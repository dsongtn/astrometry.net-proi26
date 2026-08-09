/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/*
 * Developer navigation: claim selection and wake predicates
 * ---------------------------------------------------------
 * This module decides which already-admissible unit a worker claims next:
 * an outer index while producer capacity is available, or bounded inner work
 * published by an owner. Claims occur under queue_mutex and change execution
 * state exactly once; they never transfer outer ownership or reducer
 * authority.
 *
 * Queue predicates include terminal and generation state so spurious wakeups
 * cannot claim stale work. Width policy is supplied by pass setup; this module
 * must not introduce image-, index-family-, or worker-count special cases.
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
/* queue_mutex must be held. */
static anbool index_shard_inner_work_available_locked(
    const index_shard_thread_state_t *shared) {
  if (!shared) {
    return FALSE;
  }
  return shared->helper_ready_tasks ||
      shared->staged_io_ready ||
      shared->staged_compute_ready ||
      (shared->staged_submit_backpressure
           ? shared->staged_submit_credit_ready
           : shared->staged_submit_ready) ||
      shared->staged_prepare_ready;
}

/* queue_mutex must be held. */
static anbool index_shard_queue_work_available_locked(
    const index_shard_thread_state_t *shared) {
  if (!shared || !shared->pool) {
    return FALSE;
  }
  if (shared->outer_states && shared->outer_unclaimed &&
      shared->producer_width &&
      shared->outer_running < shared->producer_width) {
    return TRUE;
  }
  return index_shard_inner_work_available_locked(shared);
}

/*
 * Outer owners wait on private condition variables while their staged group
 * remains live. They are not queue_cv waiters, but they may execute foreign
 * globally-ready work. Wake one such owner when no shared queue waiter can
 * consume the handoff; otherwise ready work can remain stranded while every
 * compute worker sleeps on a different owner_cv.
 *
 * queue_mutex must be held.
 */
static void index_shard_signal_staged_owner_locked(
    index_shard_thread_state_t *shared) {
  size_t offset;
  size_t start;

  if (!shared || !shared->pool ||
      shared->worker_count <= 0 ||
      !index_shard_inner_work_available_locked(shared)) {
    return;
  }
  start = shared->staged_owner_wake_cursor;
  if (start >= (size_t)shared->worker_count) {
    start = 0U;
  }
  for (offset = 0U;
       offset < (size_t)shared->worker_count;
       offset++) {
    size_t worker_id = start + offset;
    index_shard_worker_context_t *owner;

    if (worker_id >= (size_t)shared->worker_count) {
      worker_id -= (size_t)shared->worker_count;
    }
    owner = &shared->pool->contexts[worker_id];
    if (!owner->published_staged_group ||
        owner->published_helper_group ||
        !owner->owner_cv_ready || !owner->owner_waiting ||
        owner->owner_wake_pending) {
      continue;
    }
    shared->staged_owner_wake_cursor = worker_id + 1U;
    if (shared->staged_owner_wake_cursor >=
        (size_t)shared->worker_count) {
      shared->staged_owner_wake_cursor = 0U;
    }
    index_shard_owner_signal_locked(shared, (int)worker_id);
    return;
  }
}

/* queue_mutex must be held. */
void index_shard_queue_signal_locked(
    index_shard_thread_state_t *shared) {
  if (!shared) {
    return;
  }
  if (!shared->queue_waiters) {
    index_shard_signal_staged_owner_locked(shared);
    return;
  }
  if (!index_shard_queue_work_available_locked(shared)) {
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->queue_signals_no_work);
    }
    return;
  }
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->queue_signals);
  }
  pthread_cond_signal(&shared->queue_cv);
}

/* queue_mutex must be held. */
void index_shard_owner_signal_locked(
    index_shard_thread_state_t *shared,
    int owner_worker) {
  index_shard_worker_context_t *owner;

  if (!shared || !shared->pool || owner_worker < 0 ||
      owner_worker >= shared->worker_count) {
    return;
  }
  owner = &shared->pool->contexts[owner_worker];
  if (!owner->owner_cv_ready || !owner->owner_waiting) {
    return;
  }
  if (owner->owner_wake_pending) {
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->owner_signals_coalesced);
    }
    return;
  }
  owner->owner_wake_pending = TRUE;
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->owner_signals);
  }
  pthread_cond_signal(&owner->owner_cv);
}

/* queue_mutex must be held. */
void index_shard_notify_progress_locked(
    index_shard_thread_state_t *shared,
    int owner_worker) {
  index_shard_queue_signal_locked(shared);
  index_shard_owner_signal_locked(shared, owner_worker);
}

/*
 * Global terminal/lifecycle events must wake every worker and every owner.
 * Ordinary task transitions use targeted signals above.
 * queue_mutex must be held.
 */
void index_shard_queue_broadcast_locked(
    index_shard_thread_state_t *shared) {
  int owner_worker;

  if (!shared) {
    return;
  }
  if (shared->queue_waiters) {
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->queue_broadcasts);
    }
    pthread_cond_broadcast(&shared->queue_cv);
  }
  if (!shared->pool) {
    return;
  }
  for (owner_worker = 0;
       owner_worker < shared->worker_count;
       owner_worker++) {
    index_shard_worker_context_t *owner =
        &shared->pool->contexts[owner_worker];

    if (!owner->owner_cv_ready || !owner->owner_waiting) {
      continue;
    }
    if (owner->owner_wake_pending) {
      continue;
    }
    owner->owner_wake_pending = TRUE;
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->owner_broadcasts);
    }
    pthread_cond_broadcast(&owner->owner_cv);
  }
}

// ANCHOR INDEX-SHARD: wake-pass-waiters
void index_shard_wake_pass_waiters(index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->result_mutex);
  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);
}

void index_shard_wake_queue_waiters(
    index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->queue_mutex);
  index_shard_queue_broadcast_locked(shared);
  pthread_mutex_unlock(&shared->queue_mutex);
}

const char *index_shard_terminal_cause_name(
    index_shard_terminal_cause_t cause) {
  switch (cause) {
  case INDEX_SHARD_TERMINAL_NONE:
    return "none";
  case INDEX_SHARD_TERMINAL_WINNER:
    return "winner";
  case INDEX_SHARD_TERMINAL_WALL_LIMIT:
    return "wall-limit";
  case INDEX_SHARD_TERMINAL_CPU_LIMIT:
    return "cpu-limit";
  case INDEX_SHARD_TERMINAL_CANCELLED:
    return "cancelled";
  case INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY:
    return "global-integrity";
  }
  return "invalid";
}

/*
 * SECTION INDEX-SHARD: queue
 */
// ANCHOR INDEX-SHARD: claim-one
/*
 * Claim one shard task.
 *
 * Invariant:
 *   - each index_order is claimed at most once
 *   - stop/fatal prevents new claims
 *
 * Each participating worker owns at most one synchronous task, so the worker
 * count itself is the concurrency bound; no separate task-credit condition is
 * needed.
 */
void index_shard_worker_cleanup_pass(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared);

static void index_shard_advance_canonical_cursor_locked(
    index_shard_thread_state_t *shared) {
  while (shared->canonical_scan_cursor < shared->nindexes &&
         shared->outer_states[shared->canonical_scan_cursor] !=
             INDEX_SHARD_OUTER_UNCLAIMED) {
    shared->canonical_scan_cursor++;
  }
}

int index_shard_claim_outer_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t candidate,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice) {
  if (candidate >= shared->nindexes ||
      shared->outer_states[candidate] !=
          INDEX_SHARD_OUTER_UNCLAIMED ||
      !shared->outer_unclaimed ||
      shared->outer_running >= shared->producer_width) {
    return -1;
  }
  shared->outer_states[candidate] =
      INDEX_SHARD_OUTER_RUNNING;
  worker->ready_before_outer_eligible = FALSE;
  shared->outer_unclaimed--;
  shared->outer_running++;
  shared->outer_claims++;
  *index_order = candidate;
  *mmap_advice = shared->mmap_advice;
  if (candidate == shared->canonical_scan_cursor) {
    index_shard_advance_canonical_cursor_locked(shared);
  }

  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] claim index_order=%zu lane=producer "
           "worker=%i owners=%zu producer_width=%zu "
           "outer_unclaimed=%zu payload=%s "
           "wall_since_pass=%.6f\n",
           candidate,
           worker->worker_id,
           shared->outer_running,
           shared->producer_width,
           shared->outer_unclaimed,
           fitsbin_mmap_advice_name(*mmap_advice),
           monotonic_seconds() - shared->pass_wall_start);
  }
  /*
   * A targeted publication may wake a worker that consumes the canonical
   * outer slot before looking at already-published inner work. Hand off one
   * additional wake so that outer-first priority cannot strand that work.
   */
  index_shard_queue_signal_locked(shared);
  return 0;
}

/*
 * Select work from the current band without transferring index ownership.
 *
 * A worker that just completed one outer task may execute one already-ready
 * foreign staged computation before opening another cold index. The handoff
 * is consumed before any callback runs, so the next selection restores
 * canonical outer priority. All other inner work remains eligible only when
 * no outer task is immediately claimable. The reducer remains the only
 * authority for master result mutation.
 */
index_shard_work_selection_t
index_shard_select_work(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice,
    index_shard_inner_claim_t *inner_claim) {
  if (!worker || !shared || !index_order || !mmap_advice || !inner_claim ||
      !shared->outer_states || !shared->producer_width) {
    return INDEX_SHARD_WORK_ERROR;
  }
  if (worker->worker_id < 0 ||
      worker->worker_id >= shared->worker_count) {
    return INDEX_SHARD_WORK_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  while (1) {
    index_shard_pass_state_snapshot_t state;
    size_t candidate;
    int inner_selection;

    index_shard_pass_state_snapshot(shared, &state);
    if (state.stop_requested || state.fatal_error ||
        state.solved_published) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }

    if (worker->ready_before_outer_eligible &&
        index_shard_helper_outer_claimable_locked(shared)) {
      worker->ready_before_outer_eligible = FALSE;
      inner_selection = index_shard_staged_select_locked(
          worker,
          shared,
          INDEX_SHARD_STAGED_SELECT_COMPUTE,
          INDEX_SHARD_STAGED_SCOPE_FOREIGN,
          &inner_claim->staged);
      if (inner_selection < 0) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      if (!inner_selection) {
        inner_claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
        shared->staged_ready_before_outer_claims++;
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_HELPER;
      }
    }

    index_shard_advance_canonical_cursor_locked(shared);
    candidate = shared->canonical_scan_cursor;
    if (candidate < shared->nindexes &&
        shared->outer_running < shared->producer_width) {
      if (index_shard_claim_outer_locked(
              worker,
              shared,
              candidate,
              index_order,
              mmap_advice)) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_OUTER;
    }

    inner_selection = index_shard_inner_select_locked(
        worker, shared, FALSE, inner_claim);
    if (inner_selection < 0) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    if (!inner_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_HELPER;
    }

    if (!shared->outer_running) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }
    if (worker->local_context_ready) {
      /*
       * No new outer task can become claimable for this worker in the
       * current band. Release its index-local context before it waits for
       * helper packages from the remaining owners.
       */
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_worker_cleanup_pass(worker, shared);
      pthread_mutex_lock(&shared->queue_mutex);
      continue;
    }
    shared->queue_waiters++;
    if (shared->observability_enabled) {
      double wait_start = monotonic_seconds();

      index_shard_observability_increment(
          &shared->queue_wait_calls);
      inner_selection = pthread_cond_wait(
          &shared->queue_cv, &shared->queue_mutex);
      shared->queue_wait_seconds +=
          monotonic_seconds() - wait_start;
    } else {
      inner_selection = pthread_cond_wait(
          &shared->queue_cv, &shared->queue_mutex);
    }
    if (!shared->queue_waiters) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    shared->queue_waiters--;
    if (inner_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
  }
}
