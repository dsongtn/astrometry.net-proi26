/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/*
 * Private implementation module for the index-shard subsystem.
 * See index_shard_private.h for ownership and lock-order invariants.
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

/*
 * Developer navigation: result and terminal authority
 * --------------------------------------------------
 * Workers freeze results into per-index slots. This module classifies hook
 * failures, arbitrates completed scientific candidates, establishes terminal
 * precedence, disposes losers, and transfers at most one winner to the master
 * onefield/solver state.
 *
 * Outer quick commit and inner canonical retirement are distinct: a valid
 * completed index may win before an earlier configured index, but work inside
 * each owned index must already have retired in its required legacy order.
 * Once master-visible transfer begins, serial fallback is forbidden.
 */
/*
 * SECTION INDEX-SHARD: result
 */

// ANCHOR INDEX-SHARD: result-init
/*
 * Initialize one result slot before a worker starts solving an index.
 *
 * best_logodds starts at -HUGE_VAL so diagnostics can distinguish "no match"
 * from a real low-confidence match.
 */
static const char *index_shard_failure_class_name(
    index_shard_failure_class_t failure_class) {
  switch (failure_class) {
  case INDEX_SHARD_FAILURE_TASK_LOCAL:
    return "task-local";
  case INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY:
    return "global-integrity";
  case INDEX_SHARD_FAILURE_NONE:
  default:
    return "none";
  }
}

void index_shard_result_fail(
    index_shard_result_t *result,
    index_shard_failure_class_t failure_class,
    int rc) {
  assert(result);
  assert(failure_class != INDEX_SHARD_FAILURE_NONE);

  result->failed = TRUE;
  result->rc = rc ? rc : -1;
  if (failure_class > result->failure_class) {
    result->failure_class = failure_class;
  }
}

/*
 * Map one typed bridge result into task state.
 *
 * Return 0 for normal completion, 1 for cooperative terminal observation and
 * -1 for a classified failure. A solved outcome is accepted only from the
 * solution-analysis hook.
 */
int index_shard_apply_hook_result(
    index_shard_result_t *result,
    index_shard_hook_result_t hook_result,
    int solved_allowed) {
  assert(result);

  switch (hook_result.outcome) {
  case INDEX_SHARD_HOOK_COMPLETED_UNSOLVED:
    if (!hook_result.error_code) {
      return 0;
    }
    break;
  case INDEX_SHARD_HOOK_COMPLETED_SOLVED:
    if (solved_allowed && !hook_result.error_code) {
      result->solved = TRUE;
      return 0;
    }
    break;
  case INDEX_SHARD_HOOK_CANCELLED:
    if (!hook_result.error_code) {
      result->cancelled = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_WALL_LIMIT:
    if (!hook_result.error_code) {
      result->hit_total_timelimit = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_CPU_LIMIT:
    if (!hook_result.error_code) {
      result->hit_total_cpulimit = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_TASK_LOCAL_FAILURE:
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_TASK_LOCAL,
        hook_result.error_code);
    return -1;
  case INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE:
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        hook_result.error_code);
    return -1;
  }

  index_shard_result_fail(
      result,
      INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
      hook_result.error_code);
  return -1;
}

void index_shard_result_init(index_shard_result_t *result, size_t index_order) {
  memset(result, 0, sizeof(index_shard_result_t));

  result->index_order = index_order;
  result->best_logodds = -HUGE_VAL;
  result->best_fieldnum = -1;

  result->solutions = bl_new(4, sizeof(MatchObj));
}

void index_shard_result_dispose(index_shard_result_t *result,
                                       const index_shard_hooks_t *hooks) {
  // reducer already transferred ownership to master bp
  if (!result || !result->solutions)
    return;

  if (result->merged) {
    bl_free(result->solutions);
    result->solutions = NULL;
    return;
  }

  if (hooks && hooks->free_solutions) {
    hooks->free_solutions(result->solutions);
    result->solutions = NULL;
    return;
  }

  bl_free(result->solutions);
  result->solutions = NULL;
}

/*
 * Close the full outer-task timing interval.
 *
 * One monotonic clock call is used for both task duration and position within
 * the pass, keeping the hot-path instrumentation minimal.
 */
void index_shard_result_finish_task(
    index_shard_result_t *result,
    const index_shard_thread_state_t *shared,
    double task_wall_start) {
  double task_wall_finish;
  struct rusage resource_finish;

  assert(result);
  assert(shared);

  task_wall_finish = monotonic_seconds();

  result->task_wall_seconds =
      task_wall_finish - task_wall_start;

  result->task_finish_since_pass =
      task_wall_finish - shared->pass_wall_start;

#if defined(RUSAGE_THREAD)
  if (result->task_resource_valid &&
      getrusage(RUSAGE_THREAD, &resource_finish) == 0) {
    result->task_user_seconds =
        index_shard_timeval_delta_seconds(
            &resource_finish.ru_utime,
            &result->task_resource_start.ru_utime);
    result->task_system_seconds =
        index_shard_timeval_delta_seconds(
            &resource_finish.ru_stime,
            &result->task_resource_start.ru_stime);
    result->task_major_faults =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_majflt,
            result->task_resource_start.ru_majflt);
    result->task_input_blocks =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_inblock,
            result->task_resource_start.ru_inblock);
    result->task_voluntary_switches =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_nvcsw,
            result->task_resource_start.ru_nvcsw);
  } else {
    result->task_resource_valid = FALSE;
  }
#else
  (void)resource_finish;
  result->task_resource_valid = FALSE;
#endif
}

// ANCHOR INDEX-SHARD: analyze-result
/*
 * Inspect worker-local solutions without merging them.
 *
 * This marks solved/best_logodds before immutable result publication elects
 * a winner. Only the reducer may transfer the selected solution into master
 * state.
 */
int index_shard_capture_solution_analysis(
    index_shard_thread_state_t *shared,
    index_shard_result_t *result) {
  index_shard_hook_result_t hook_result;

  if (!shared->hooks || !shared->hooks->analyze_solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  hook_result = shared->hooks->analyze_solutions(
      shared->bp, result->solutions,
      &result->best_logodds, &result->best_fieldnum);
  return index_shard_apply_hook_result(
      result, hook_result, TRUE);
}
// ANCHOR INDEX-SHARD: reduce-one-result
static void index_shard_reducer_work_account(
    index_shard_thread_state_t *shared,
    double wall_start) {
  double wall_seconds;

  assert(shared);

  wall_seconds = monotonic_seconds() - wall_start;

  if (wall_seconds < 0.0) {
    wall_seconds = 0.0;
  }

  shared->reducer_work_calls++;
  shared->reducer_work_wall_seconds += wall_seconds;
}

/*
 * Transfer one completed worker result into master onefield state.
 *
 * Only the reducer calls this.  Workers never append directly into
 * master_bp->solutions.
 *
 * The merge hook is not transactional.  Before transferring a non-empty
 * worker result, mark the pass as master-committed so a later failure cannot
 * trigger an unsafe serial rerun.
 */
static int index_shard_reduce_one_result(index_shard_thread_state_t *shared,
                                         index_shard_result_t *result) {
  anbool solved = FALSE;
  double reducer_wall_start;
  int may_mutate_master = FALSE;
  int losing_result = FALSE;

  if (!result || result->merged) {
    return 0;
  }

  pthread_mutex_lock(&shared->state_mutex);
  losing_result =
      shared->winner_selected &&
      result->index_order != shared->selected_index_order;
  pthread_mutex_unlock(&shared->state_mutex);
  if (losing_result) {
    logerr("[index-shard] refusing to merge losing result "
           "index_order=%zu\n",
           result->index_order);
    return -1;
  }

  if (shared->have_committed_result) {
    logerr("[index-shard] refusing a reducer merge after solution commit\n");
    return -1;
  }

  reducer_wall_start = monotonic_seconds();

  if (result->failed || result->rc) {
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  if (!shared->hooks || !shared->hooks->merge_solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  if (result->solutions && bl_size(result->solutions) > 0) {
    may_mutate_master = TRUE;
  }

  if (result->solved) {
    may_mutate_master = TRUE;
  }

  if (may_mutate_master) {
    index_shard_mark_master_committed(shared);
  }

  if (shared->hooks->merge_solutions(shared->bp,
                                     result->solutions,
                                     &solved)) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  result->merged = TRUE;
  shared->results_reduced++;

  if (solved || result->solved) {
    result->solved = TRUE;
    shared->have_committed_result = TRUE;
    shared->committed_index_order = result->index_order;
    index_shard_publish_committed_solve(shared);
  }

  index_shard_reducer_work_account(shared, reducer_wall_start);
  return 0;
}
// ANCHOR INDEX-SHARD: worker-get-index
/*
 * Freeze one index-independent scientific result and arbitrate the first
 * terminal event before index cleanup or producer-slot release.
 *
 * result_mutex -> state_mutex -> limit_mutex is the only lock order used
 * here. All three are released before queue_mutex can be acquired later by
 * index_shard_finish_outer_claim(). Result lifecycle metrics may still be
 * completed by the owner, but solutions, solved state and candidate identity
 * are immutable after candidate_ready is published.
 */
int index_shard_arbitrate_candidate(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  index_shard_result_t *result;
  index_shard_terminal_cause_t cause =
      INDEX_SHARD_TERMINAL_NONE;
  double elapsed = 0.0;
  int rc = 0;
  int report = FALSE;
  int selected_now = FALSE;
  int stop_now = FALSE;

  if (!shared || !shared->results ||
      index_order >= shared->nindexes) {
    return -1;
  }

  result = &shared->results[index_order];

  pthread_mutex_lock(&shared->result_mutex);
  if (result->candidate_ready ||
      (shared->completed && shared->completed[index_order])) {
    logerr("[index-shard] duplicate or late candidate publication "
           "index_order=%zu\n",
           index_order);
    rc = -1;
  } else {
    if (result->solved &&
        (result->failed || result->rc ||
         result->failure_class != INDEX_SHARD_FAILURE_NONE)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    result->candidate_sequence =
        ++shared->next_candidate_sequence;
    result->candidate_ready = TRUE;

    pthread_mutex_lock(&shared->state_mutex);
    if (result->failure_class ==
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY) {
      (void)index_shard_publish_terminal_locked(
          shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
    } else if (shared->terminal_cause ==
               INDEX_SHARD_TERMINAL_NONE) {
      pthread_mutex_lock(&shared->limit_mutex);
      cause = index_shard_sample_terminal_locked(
          shared, result, &elapsed, &report);
      if (cause != INDEX_SHARD_TERMINAL_NONE) {
        (void)index_shard_publish_terminal_locked(
            shared, cause);
      } else if (result->solved &&
                 !result->failed &&
                 result->rc == 0 &&
                 result->failure_class ==
                     INDEX_SHARD_FAILURE_NONE &&
                 !shared->winner_selected &&
                 !shared->solved_published) {
        shared->winner_selected = TRUE;
        shared->selected_index_order = index_order;
        shared->selected_candidate_sequence =
            result->candidate_sequence;
        selected_now = index_shard_publish_terminal_locked(
            shared, INDEX_SHARD_TERMINAL_WINNER);
      }
      pthread_mutex_unlock(&shared->limit_mutex);
    }
    stop_now = shared->stop_requested;
    pthread_mutex_unlock(&shared->state_mutex);
  }

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);

  if (stop_now) {
    index_shard_publish_worker_stop(shared);
  }

  if (report && cause == INDEX_SHARD_TERMINAL_WALL_LIMIT) {
    logmsg("Total wall-clock time limit reached!\n");
    logverb("[index-shard] wall-limit reached total_timelimit=%g "
            "elapsed=%.6f\n",
            shared->bp->total_timelimit,
            elapsed);
  } else if (report && cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) {
    logmsg("Total CPU time limit reached!\n");
    logverb("[index-shard] cpu-budget reached total_cpulimit=%g "
            "elapsed=%.6f\n",
            shared->bp->total_cpulimit,
            elapsed);
  }

  if (selected_now) {
    logverb("[index-shard] winner-selected index_order=%zu worker=%i "
            "candidate_sequence=%zu field=%i best_logodds=%.17g "
            "pass_wall=%.6f\n",
            index_order,
            result->worker_id,
            result->candidate_sequence,
            result->best_fieldnum,
            result->best_logodds,
            shared->first_stop_wall_since_pass);
  }

  return rc;
}

/*
 * Publish full task completion after candidate arbitration and index cleanup.
 *
 * candidate_ready protects the immutable scientific payload. completed[]
 * protects the remaining owner lifecycle and metrics. A cleanup failure can
 * still promote an elected-but-uncommitted winner to global-integrity failure,
 * but this function never performs winner election or releases outer producer
 * capacity.
 */
static int index_shard_mark_result_completed(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  index_shard_result_t *result;
  int fatal_now = FALSE;
  int late_loser_failure = FALSE;
  int rc = 0;
  int stop_now = FALSE;
  int task_local_now = FALSE;

  if (!shared || !shared->results || !shared->completed ||
      index_order >= shared->nindexes) {
    return -1;
  }

  result = &shared->results[index_order];

  /*
   * NOTE INDEX-SHARD: claimed-task-invariant
   *
   * Every claimed task must mark completion exactly once. Otherwise the
   * reducer can wait forever.
   */
  pthread_mutex_lock(&shared->result_mutex);
  if (shared->completed[index_order]) {
    logerr("[index-shard] duplicate result completion "
           "index_order=%zu\n",
           index_order);
    rc = -1;
  } else {
    /*
     * Freeze a complete cause-based terminal classification before publishing
     * the immutable result slot. Unclassified or internally inconsistent
     * failure state is itself a pass-integrity failure.
     */
    if ((result->failed || result->rc) &&
        result->failure_class == INDEX_SHARD_FAILURE_NONE) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    } else if (result->failure_class != INDEX_SHARD_FAILURE_NONE &&
               (!result->failed || !result->rc)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    result->completion_sequence =
        ++shared->next_completion_sequence;
    shared->completed[index_order] = TRUE;

    pthread_mutex_lock(&shared->state_mutex);

    if (result->solved && !result->candidate_ready) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }
    if (shared->winner_selected &&
        shared->selected_index_order == index_order &&
        (!result->candidate_ready ||
         !result->candidate_sequence ||
         result->candidate_sequence !=
             shared->selected_candidate_sequence ||
         result->failed || result->rc ||
         result->failure_class != INDEX_SHARD_FAILURE_NONE)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    if (result->failure_class ==
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY) {
      shared->global_integrity_failures++;
      (void)index_shard_publish_terminal_locked(
          shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
      fatal_now = TRUE;
      stop_now = shared->stop_requested;
    } else if (result->failure_class ==
               INDEX_SHARD_FAILURE_TASK_LOCAL) {
      shared->task_local_failures++;
      task_local_now = TRUE;
      if (shared->winner_selected &&
          shared->selected_index_order != index_order) {
        shared->late_loser_failures++;
        late_loser_failure = TRUE;
      }
    }

    pthread_mutex_unlock(&shared->state_mutex);
  }

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);

  if (stop_now) {
    index_shard_publish_worker_stop(shared);
  }

  if (fatal_now) {
    logerr("[index-shard] global integrity failure terminated pass "
           "index_order=%zu completion_sequence=%zu class=%s\n",
           index_order,
           result->completion_sequence,
           index_shard_failure_class_name(result->failure_class));
  } else if (late_loser_failure) {
    logverb("[index-shard] preserving selected winner after task-local "
            "loser index_order=%zu completion_sequence=%zu\n",
            index_order,
            result->completion_sequence);
  } else if (task_local_now) {
    logverb("[index-shard] isolated task-local failure "
            "index_order=%zu completion_sequence=%zu\n",
            index_order,
            result->completion_sequence);
  }

  return rc;
}

void index_shard_finish_outer_claim(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t index_order) {
  int underflow = FALSE;
  int completion_failed;

  /*
   * Candidate arbitration has already frozen any scientific result and
   * published its terminal event. This final limit sample covers tasks that
   * produced no candidate or crossed a limit during cleanup. Full completion
   * remains visible before the producer slot becomes claimable.
   */
  (void)index_shard_check_global_limits(shared);
  completion_failed = index_shard_mark_result_completed(
      shared, index_order);
  if (completion_failed) {
    index_shard_request_fatal_stop(shared);
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (!worker || !shared->outer_states ||
      index_order >= shared->nindexes ||
      shared->outer_states[index_order] !=
          INDEX_SHARD_OUTER_RUNNING) {
    logerr("[index-shard] invalid outer lifecycle "
           "index_order=%zu\n",
           index_order);
    underflow = TRUE;
  } else {
    shared->outer_states[index_order] =
        INDEX_SHARD_OUTER_FINISHED;
    worker->ready_before_outer_eligible = TRUE;
  }
  if (!shared->outer_running) {
    logerr("[index-shard] outer-running underflow\n");
    underflow = TRUE;
  } else {
    shared->outer_running--;
  }
  if (!shared->outer_running) {
    index_shard_queue_broadcast_locked(shared);
  } else if (shared->outer_unclaimed &&
             shared->outer_running < shared->producer_width) {
    index_shard_queue_signal_locked(shared);
  }
  pthread_mutex_unlock(&shared->queue_mutex);

  if (underflow) {
    index_shard_request_fatal_stop(shared);
  }
}

/*
 * SECTION INDEX-SHARD: reducer
 *
 * Main-thread result publication.
 *
 * Workers fill result slots.  The reducer is the only path that transfers
 * MatchObj data into master bp->solutions and updates final master solved state.
 */

// ANCHOR INDEX-SHARD: reduce-first-valid
/*
 * First-valid reducer for one submitted pass.
 *
 * Result publication selects the first immutable verified completion and
 * starts cooperative cancellation. After all workers quiesce, reduce exactly
 * that selected result. Configured index order is used only for the clean
 * unsolved drain that preserves original diagnostic accumulation.
 */
#define INDEX_SHARD_LIMIT_POLL_NANOSECONDS 100000000L

static void index_shard_limit_poll_deadline(
    struct timespec *deadline) {
  clock_gettime(CLOCK_MONOTONIC, deadline);
  deadline->tv_nsec +=
      INDEX_SHARD_LIMIT_POLL_NANOSECONDS;
  if (deadline->tv_nsec >= 1000000000L) {
    deadline->tv_sec++;
    deadline->tv_nsec -= 1000000000L;
  }
}

int index_shard_pool_reduce_first_valid(index_shard_pool_t *pool) {
  index_shard_thread_state_t *shared = &pool->shared;
  index_shard_pass_state_snapshot_t state;
  size_t i;
  int rc = 0;

  /*
   * Winner selection already published stop. Keep all result storage and
   * shared field state alive until every owner and borrowed helper has left
   * this generation.
   */
  pthread_mutex_lock(&shared->result_mutex);
  while (shared->active_workers > 0) {
    struct timespec deadline;
    int wait_result;

    index_shard_limit_poll_deadline(&deadline);
    wait_result = pthread_cond_timedwait(
        &shared->result_cv,
        &shared->result_mutex,
        &deadline);
    if (wait_result == ETIMEDOUT) {
      pthread_mutex_unlock(&shared->result_mutex);
      (void)index_shard_check_global_limits(shared);
      pthread_mutex_lock(&shared->result_mutex);
    } else if (wait_result) {
      rc = -1;
      pthread_mutex_unlock(&shared->result_mutex);
      index_shard_request_fatal_stop(shared);
      pthread_mutex_lock(&shared->result_mutex);
    }
  }

  index_shard_pass_state_snapshot(shared, &state);
  pthread_mutex_unlock(&shared->result_mutex);

  if (state.fatal_error) {
    return -1;
  }
  if (rc) {
    return rc;
  }

  if (state.winner_selected) {
    index_shard_result_t *winner;

    if (state.selected_index_order >= shared->nindexes ||
        !shared->completed[state.selected_index_order]) {
      logerr("[index-shard] selected winner result is unavailable "
             "index_order=%zu\n",
             state.selected_index_order);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    winner = &shared->results[state.selected_index_order];
    if (!winner->candidate_ready ||
        !winner->solved ||
        winner->failed ||
        winner->rc ||
        winner->failure_class != INDEX_SHARD_FAILURE_NONE ||
        winner->cancelled ||
        winner->hit_total_timelimit ||
        winner->hit_total_cpulimit ||
        winner->merged ||
        !winner->candidate_sequence ||
        winner->candidate_sequence !=
            state.selected_candidate_sequence ||
        !winner->completion_sequence ||
        !shared->completed[state.selected_index_order]) {
      logerr("[index-shard] selected winner result is inconsistent "
             "index_order=%zu candidate_sequence=%zu\n",
             state.selected_index_order,
             state.selected_candidate_sequence);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (index_shard_trace_enabled()) {
      logmsg("[index-shard] reduce-winner index_order=%zu worker=%i "
             "candidate_sequence=%zu completion_sequence=%zu "
             "solved=%i failed=%i\n",
             winner->index_order,
             winner->worker_id,
             winner->candidate_sequence,
             winner->completion_sequence,
             winner->solved,
             winner->failed);
    }

    if (index_shard_reduce_one_result(shared, winner) ||
        !shared->have_committed_result ||
        shared->committed_index_order !=
            state.selected_index_order) {
      logerr("[index-shard] failed to commit selected winner "
             "index_order=%zu\n",
             state.selected_index_order);
      index_shard_request_fatal_stop(shared);
      return -1;
    }
    return 0;
  }

  /* Cancellation or a limit won before any immutable verified result. */
  if (state.stop_requested) {
    return 0;
  }

  if (state.task_local_failures) {
    logerr("[index-shard] pass exhausted with task-local failures "
           "count=%llu; requesting exact serial retry\n",
           state.task_local_failures);
    return -1;
  }

  /*
   * A clean unsolved pass preserves configured-order accumulation of
   * below-threshold diagnostics. No result in this path may be solved.
   */
  for (i = 0; i < shared->nindexes; i++) {
    index_shard_result_t *result = &shared->results[i];

    if (!shared->completed[i] ||
        result->failed ||
        result->rc ||
        result->failure_class != INDEX_SHARD_FAILURE_NONE ||
        result->solved) {
      logerr("[index-shard] invalid clean-unsolved result "
             "index_order=%zu completed=%i solved=%i failed=%i rc=%i\n",
             i,
             shared->completed[i] ? 1 : 0,
             result->solved,
             result->failed,
             result->rc);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (index_shard_reduce_one_result(shared, result)) {
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (shared->solved_published) {
      logerr("[index-shard] clean-unsolved reduction committed a solution "
             "index_order=%zu\n", i);
      index_shard_request_fatal_stop(shared);
      return -1;
    }
  }

  return 0;
}

/*
 * Report the reducer-owned winner after every worker has quiesced and the
 * selected result has committed. Worker callbacks deliberately suppress their
 * ordinary solved line, so this is the sole parallel success report.
 */
int index_shard_report_committed_solution(
    onefield_t *bp,
    size_t nindexes,
    const index_shard_thread_state_t *shared,
    const index_shard_hooks_t *hooks,
    const index_shard_result_t *results) {
  const index_shard_result_t *committed;

  if (!hooks || !hooks->report_committed_solution) {
    logerr("[index-shard] committed-solution reporter is unavailable\n");
    return -1;
  }

  if (!shared ||
      !shared->winner_selected ||
      !shared->have_committed_result ||
      shared->committed_index_order >= nindexes ||
      shared->committed_index_order != shared->selected_index_order) {
    logerr("[index-shard] committed-solution identity is unavailable\n");
    return -1;
  }

  committed = &results[shared->committed_index_order];
  if (!committed->merged ||
      !committed->solved ||
      committed->failed ||
      committed->rc != 0) {
    logerr("[index-shard] committed-solution result is inconsistent\n");
    return -1;
  }

  return hooks->report_committed_solution(
      bp,
      committed->index_order,
      committed->best_fieldnum,
      committed->best_logodds);
}
