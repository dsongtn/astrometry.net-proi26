/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/*
 * Developer navigation: pass-wide control
 * ---------------------------------------
 * This module owns worker TLS, global-limit sampling, terminal-state queries,
 * and cooperative stop propagation. TLS lets callbacks reached through the
 * legacy solver stack find their current worker/pass without adding mutable
 * pool state to upstream function signatures.
 *
 * Control code may request unwind, but it does not arbitrate candidates or
 * transfer results. Those are reducer responsibilities. All waits and state
 * changes follow the locks and terminal precedence in index_shard_private.h.
 */
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "index_shard_private.h"
#include "astrometry/bl.h"
#include "astrometry/errors.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/fitsbin.h"
#include "astrometry/fitsioutils.h"

static pthread_key_t index_shard_tls_key;
static pthread_once_t index_shard_tls_once = PTHREAD_ONCE_INIT;
static int index_shard_tls_key_status = EAGAIN;

/*
 * SECTION INDEX-SHARD: tls - thread logical singleton
 *
 * TLS links a running worker callback back to its pool state.
 *
 * solve_fields()/solver_run() call into onefield callbacks.  Those callbacks
 * receive only local onefield_t, so TLS is used to check the global shard stop
 * state and set local_bp->solver.quit_now.
 */
static void index_shard_make_tls_key(void) {
  index_shard_tls_key_status =
      pthread_key_create(&index_shard_tls_key, NULL);
}

int index_shard_tls_ensure(void) {
  int once_status =
      pthread_once(&index_shard_tls_once,
                   index_shard_make_tls_key);

  if (once_status) {
    return once_status;
  }
  return index_shard_tls_key_status;
}

int index_shard_set_tls(index_shard_worker_context_t *ctx) {
  int status = index_shard_tls_ensure();

  if (status) {
    return status;
  }
  return pthread_setspecific(index_shard_tls_key, ctx);
}

index_shard_worker_context_t *index_shard_get_tls(void) {
  if (index_shard_tls_ensure()) {
    return NULL;
  }
  return pthread_getspecific(index_shard_tls_key);
}

anbool index_shard_worker_context_active(void) {
  return index_shard_get_tls() != NULL;
}

anbool index_shard_worker_stop_requested(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();

  if (!ctx || !ctx->pool) {
    return FALSE;
  }

  return __atomic_load_n(
      &ctx->pool->shared.worker_stop_requested,
      __ATOMIC_ACQUIRE) != 0;
}

/*
 * SECTION INDEX-SHARD: configuration
 *
 * The engine resolves config, environment, and per-job overrides before
 * engine_run_job(). Hot execution paths consume the immutable job value only.
 */
anbool index_shard_pthread_enabled(const onefield_t *bp) {
  return bp && bp->index_shard_workers > 1;
}

anbool index_shard_trace_enabled(void) {
  /*
   * Detailed scheduler traces follow the ordinary command-line verbosity
   * model. Two -v flags select LOG_ALL without another environment control.
   */
  return log_get_level() >= LOG_ALL;
}

int index_shard_get_worker_count(const onefield_t *bp,
                                        size_t nindexes) {
  (void)nindexes;

  if (!bp) {
    return 1;
  }

  return index_shard_config_effective_workers(
      bp->index_shard_workers,
      0);
}

/* state_mutex must be held. */
int index_shard_publish_terminal_locked(
    index_shard_thread_state_t *shared,
    index_shard_terminal_cause_t cause) {
  int changed = FALSE;

  assert(shared);
  if (cause == INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  if (cause == INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY) {
    changed = !shared->fatal_error ||
        shared->terminal_cause != cause;
    if (shared->terminal_cause == INDEX_SHARD_TERMINAL_NONE) {
      shared->first_stop_wall_since_pass =
          monotonic_seconds() - shared->pass_wall_start;
    }
    shared->terminal_cause = cause;
    shared->fatal_error = TRUE;
    shared->stop_requested = TRUE;
    return changed;
  }

  if (shared->terminal_cause != INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  shared->terminal_cause = cause;
  shared->first_stop_wall_since_pass =
      monotonic_seconds() - shared->pass_wall_start;
  shared->stop_requested = TRUE;
  return TRUE;
}

void index_shard_publish_worker_stop(
    index_shard_thread_state_t *shared) {
  __atomic_store_n(
      &shared->worker_stop_requested,
      TRUE,
      __ATOMIC_RELEASE);
  index_shard_wake_pass_waiters(shared);
  index_shard_wake_queue_waiters(shared);
}

// ANCHOR INDEX-SHARD: request-fatal-stop
/*
 * Publish a hard worker/module failure and stop the current pass.
 */
void index_shard_request_fatal_stop(index_shard_thread_state_t *shared) {
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  (void)index_shard_publish_terminal_locked(
      shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
  stopped = shared->stop_requested;
  pthread_mutex_unlock(&shared->state_mutex);

  if (stopped) {
    index_shard_publish_worker_stop(shared);
  }
}

// ANCHOR INDEX-SHARD: publish-committed-solve
/*
 * Publish that the reducer committed a valid solved result.
 *
 * First-valid selection already stopped new work. This publication records
 * that the reducer transferred the selected result into master state.
 */
void index_shard_publish_committed_solve(
    index_shard_thread_state_t *shared) {
  int valid;

  pthread_mutex_lock(&shared->state_mutex);
  valid = shared->winner_selected &&
      shared->terminal_cause == INDEX_SHARD_TERMINAL_WINNER &&
      !shared->fatal_error;
  if (valid) {
    shared->solved_published = TRUE;
  } else {
    (void)index_shard_publish_terminal_locked(
        shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
  }
  pthread_mutex_unlock(&shared->state_mutex);

  index_shard_publish_worker_stop(shared);
  if (!valid) {
    logerr("[index-shard] invalid committed-solution terminal state\n");
  }
}
/*
 * Mark the point after which serial fallback is no longer safe.
 *
 * The merge hook is not transactional.  Once it begins transferring a
 * non-empty worker result into master-visible state, a later failure must not
 * cause the caller to rerun the original serial pass.
 */
void index_shard_mark_master_committed(
    index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->state_mutex);
  shared->master_committed = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);
}

// ANCHOR INDEX-SHARD: master-limit-state
int index_shard_master_limit_or_cancel_requested(
    index_shard_thread_state_t *shared,
    anbool *hit_total_cpulimit,
    anbool *hit_total_timelimit,
    anbool *cancelled) {
  onefield_t *bp = shared->bp;
  int stop;

  pthread_mutex_lock(&shared->limit_mutex);

  if (hit_total_cpulimit) {
    *hit_total_cpulimit = bp->hit_total_cpulimit;
  }

  if (hit_total_timelimit) {
    *hit_total_timelimit = bp->hit_total_timelimit;
  }

  if (cancelled) {
    *cancelled = bp->cancelled;
  }

  stop = bp->hit_total_cpulimit || bp->hit_total_timelimit || bp->cancelled;

  pthread_mutex_unlock(&shared->limit_mutex);

  return stop;
}

// ANCHOR INDEX-SHARD: master-stop-check
/*
 * Read-only stop predicate used by workers before expensive work.
 *
 * First-valid selection is mirrored through stop_requested. Workers therefore
 * do not read master bp->single_field_solved concurrently with the reducer.
 */
int index_shard_master_stop_requested(index_shard_thread_state_t *shared) {
  index_shard_pass_state_snapshot_t state;

  index_shard_pass_state_snapshot(shared, &state);

  if (state.stop_requested || state.fatal_error || state.solved_published) {
    return TRUE;
  }

  return index_shard_master_limit_or_cancel_requested(
      shared, NULL, NULL, NULL);
}

/* state_mutex and limit_mutex must be held. */
index_shard_terminal_cause_t
index_shard_sample_terminal_locked(
    index_shard_thread_state_t *shared,
    const index_shard_result_t *result,
    double *elapsed,
    int *report) {
  onefield_t *bp;
  index_shard_terminal_cause_t cause =
      INDEX_SHARD_TERMINAL_NONE;

  assert(shared);
  bp = shared->bp;
  assert(bp);

  if (elapsed) {
    *elapsed = 0.0;
  }
  if (report) {
    *report = FALSE;
  }

  if (result && result->cancelled) {
    bp->cancelled = TRUE;
  }
  if (result && result->hit_total_timelimit) {
    bp->hit_total_timelimit = TRUE;
  }
  if (result && result->hit_total_cpulimit) {
    bp->hit_total_cpulimit = TRUE;
  }

  if (bp->cancelled) {
    cause = INDEX_SHARD_TERMINAL_CANCELLED;
  } else if (bp->hit_total_timelimit) {
    cause = INDEX_SHARD_TERMINAL_WALL_LIMIT;
  } else if (bp->hit_total_cpulimit) {
    cause = INDEX_SHARD_TERMINAL_CPU_LIMIT;
  } else if (bp->total_timelimit > 0.0) {
    double now = monotonic_seconds();
    double sampled = now - bp->time_total_start;

    if (elapsed) {
      *elapsed = sampled;
    }
    if (now >= 0.0 && sampled >= bp->total_timelimit) {
      bp->hit_total_timelimit = TRUE;
      cause = INDEX_SHARD_TERMINAL_WALL_LIMIT;
    }
  }

  if (cause == INDEX_SHARD_TERMINAL_NONE &&
      bp->total_cpulimit > 0.0) {
    double sampled =
        (double)(get_cpu_usage() - bp->cpu_total_start);

    if (elapsed) {
      *elapsed = sampled;
    }
    if (sampled >= bp->total_cpulimit) {
      bp->hit_total_cpulimit = TRUE;
      cause = INDEX_SHARD_TERMINAL_CPU_LIMIT;
    }
  }

  if ((cause == INDEX_SHARD_TERMINAL_WALL_LIMIT ||
       cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) &&
      !shared->limit_reported) {
    shared->limit_reported = TRUE;
    if (report) {
      *report = TRUE;
    }
  }

  return cause;
}
// ANCHOR INDEX-SHARD: global-limits
/*
 * Process-wide elapsed-time and CPU-budget checks.
 *
 * total_timelimit is one shared monotonic wall-clock deadline and is not
 * divided by worker count. total_cpulimit is aggregate process CPU time; with
 * N active threads it can be consumed roughly N times faster than wall time.
 */
int index_shard_check_global_limits(index_shard_thread_state_t *shared) {
  onefield_t *bp = shared->bp;
  index_shard_terminal_cause_t cause = INDEX_SHARD_TERMINAL_NONE;
  double elapsed = 0.0;
  int report = FALSE;

  /*
   * state_mutex is the terminal-event arbiter. Holding it while the master
   * limit flags are inspected and updated closes the former interval between
   * deadline publication and cooperative-stop publication.
   */
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->terminal_cause != INDEX_SHARD_TERMINAL_NONE) {
    pthread_mutex_unlock(&shared->state_mutex);
    return TRUE;
  }

  pthread_mutex_lock(&shared->limit_mutex);
  cause = index_shard_sample_terminal_locked(
      shared, NULL, &elapsed, &report);
  if (cause != INDEX_SHARD_TERMINAL_NONE) {
    (void)index_shard_publish_terminal_locked(shared, cause);
  }
  pthread_mutex_unlock(&shared->limit_mutex);
  pthread_mutex_unlock(&shared->state_mutex);

  if (cause == INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  index_shard_publish_worker_stop(shared);
  if (report && cause == INDEX_SHARD_TERMINAL_WALL_LIMIT) {
    logmsg("Total wall-clock time limit reached!\n");
    logverb("[index-shard] wall-limit reached total_timelimit=%g "
            "elapsed=%.6f\n",
            bp->total_timelimit,
            elapsed);
  } else if (report && cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) {
    logmsg("Total CPU time limit reached!\n");
    logverb("[index-shard] cpu-budget reached total_cpulimit=%g "
            "elapsed=%.6f\n",
            bp->total_cpulimit,
            elapsed);
  }
  return TRUE;
}

// ANCHOR INDEX-SHARD: callback-poll
/*
 * Called from onefield callbacks/timer paths while solver_run() is active.
 *
 * This is the fast-stop path for workers already inside solver code.
 * It maps global stop -> local solver.quit_now.
 */
void index_shard_poll_from_callback(onefield_t *bp) {
  // non-worker callback, nothing to do
  // local solver should unwind normally through existing quit path
  index_shard_worker_context_t *ctx = index_shard_get_tls();

  if (!ctx || !ctx->pool)
    return;

  if (index_shard_check_global_limits(&ctx->pool->shared)) {
    bp->solver.quit_now = TRUE;
    return;
  }

  if (index_shard_master_stop_requested(&ctx->pool->shared))
    bp->solver.quit_now = TRUE;
}
