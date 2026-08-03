/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/*
 * Developer navigation: bounded helper execution
 * ----------------------------------------------
 * Helpers execute already-published immutable work groups for an outer owner.
 * A successful claim transfers execution of that group once, not ownership of
 * the owner's index, solver, reducer cursor, or result. The owner may execute
 * unclaimed groups itself and remains the only authority that retires outputs.
 *
 * Helper callbacks must be independently runnable, bounded, and nonblocking on
 * child work in this same pool. Cancellation retires the running lifetime
 * lease and returns through the group's predicate; it never abandons storage.
 * See index_shard_private.h for lock order and index_shard_staged.c for the
 * asynchronous completion registry.
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
/* Helper and staged claim execution, cancellation, and retirement. */

anbool index_shard_helper_outer_claimable_locked(
    const index_shard_thread_state_t *shared) {
  size_t candidate;

  if (!shared || !shared->outer_states ||
      !shared->producer_width ||
      shared->outer_running >= shared->producer_width) {
    return FALSE;
  }
  candidate = shared->canonical_scan_cursor;
  while (candidate < shared->nindexes &&
         shared->outer_states[candidate] !=
             INDEX_SHARD_OUTER_UNCLAIMED) {
    candidate++;
  }
  return candidate < shared->nindexes;
}

static size_t index_shard_helper_idle_workers_locked(
    const index_shard_thread_state_t *shared) {
  size_t available = 0U;
  size_t limit;
  size_t spare = 0U;

  if (!shared || shared->worker_count < 2) {
    return 0U;
  }
  limit = (size_t)shared->worker_count - 1U;
  if (index_shard_helper_outer_claimable_locked(shared)) {
    return 0U;
  }
  if ((size_t)shared->worker_count > shared->outer_running) {
    spare = (size_t)shared->worker_count - shared->outer_running;
  }
  available = shared->queue_waiters;
  if (!shared->helper_groups_active && spare > available) {
    available = spare;
  }
  available = MIN(available, limit);
  if (shared->helper_foreign_reservations >= available) {
    return 0U;
  }
  return available - shared->helper_foreign_reservations;
}

/* queue_mutex must be held. */
static int index_shard_helper_release_foreign_reservations_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t count) {
  if (!shared || !group ||
      count > group->foreign_reservations_outstanding ||
      count > shared->helper_foreign_reservations) {
    if (group) {
      group->internal_error = TRUE;
    }
    return -1;
  }
  group->foreign_reservations_outstanding -= count;
  shared->helper_foreign_reservations -= count;
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_helper_account_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  if (!shared || !group || group->ready_accounted ||
      SIZE_MAX - shared->helper_ready_tasks < group->ready_count) {
    return -1;
  }
  shared->helper_ready_tasks += group->ready_count;
  group->ready_accounted = TRUE;
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_helper_remove_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t count) {
  if (!shared || !group || count > group->ready_count) {
    return -1;
  }
  if (group->ready_accounted) {
    if (count > shared->helper_ready_tasks) {
      return -1;
    }
    shared->helper_ready_tasks -= count;
  }
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_helper_unaccount_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  if (!shared || !group || !group->ready_accounted) {
    return -1;
  }
  if (group->ready_count > shared->helper_ready_tasks) {
    return -1;
  }
  shared->helper_ready_tasks -= group->ready_count;
  group->ready_accounted = FALSE;
  return 0;
}

/* queue_mutex must be held. */
static void index_shard_helper_cancel_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    index_shard_helper_task_status_t status) {
  size_t i;

  if (!shared || !group || !group->tasks) {
    return;
  }
  (void)index_shard_helper_release_foreign_reservations_locked(
      shared,
      group,
      group->foreign_reservations_outstanding);
  if (index_shard_helper_remove_ready_locked(
          shared, group, group->ready_count)) {
    group->internal_error = TRUE;
  }
  for (i = 0U; i < group->task_count; i++) {
    index_shard_helper_task_t *task = &group->tasks[i];

    if (task->scheduler_state !=
        INDEX_SHARD_HELPER_TASK_READY) {
      continue;
    }
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_DONE;
    task->execute_status = status;
    if (group->ready_count) {
      group->ready_count--;
    } else {
      group->internal_error = TRUE;
    }
    if (group->completed_count < group->task_count) {
      group->completed_count++;
    } else {
      group->internal_error = TRUE;
    }
  }
  if (group->ready_count) {
    group->internal_error = TRUE;
    group->ready_count = 0U;
  }
  group->ready_work = 0U;
  if (status == INDEX_SHARD_HELPER_TASK_ERROR) {
    group->task_failed = TRUE;
  } else {
    group->stop_seen = TRUE;
  }
}

/* queue_mutex must be held. */
static int index_shard_helper_claim_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t *task_index) {
  index_shard_helper_task_t *task;
  unsigned long long work;
  size_t candidate;

  if (!shared || !group || !task_index || !group->tasks) {
    return -1;
  }
  if (!group->ready_count) {
    return 1;
  }

  candidate = group->next_claim;
  while (candidate < group->task_count &&
         group->tasks[candidate].scheduler_state !=
             INDEX_SHARD_HELPER_TASK_READY) {
    candidate++;
  }
  if (candidate >= group->task_count) {
    group->internal_error = TRUE;
    return -1;
  }

  task = &group->tasks[candidate];
  work = index_shard_helper_task_work(task);
  if (group->ready_work < work || !group->ready_count) {
    group->internal_error = TRUE;
    return -1;
  }

  if (index_shard_helper_remove_ready_locked(
          shared, group, 1U)) {
    group->internal_error = TRUE;
    return -1;
  }
  task->scheduler_state = INDEX_SHARD_HELPER_TASK_RUNNING;
  group->next_claim = candidate + 1U;
  group->ready_count--;
  group->running_count++;
  group->ready_work -= work;
  if (group->running_count > group->max_running) {
    group->max_running = group->running_count;
  }
  *task_index = candidate;
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_helper_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    index_shard_helper_claim_t *claim) {
  index_shard_helper_group_t *best = NULL;
  int owner;

  if (!worker || !worker->pool || !shared || !claim) {
    return -1;
  }
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->stop_requested || shared->fatal_error ||
      shared->solved_published) {
    pthread_mutex_unlock(&shared->state_mutex);
    return 1;
  }

  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_helper_group_t *candidate =
        worker->pool->contexts[owner].published_helper_group;

    if (!candidate || !candidate->ready_count) {
      continue;
    }
    if (candidate->generation != worker->generation_seen ||
        candidate->owner_worker != owner ||
        candidate->owner_epoch !=
            worker->pool->contexts[owner].helper_group_epoch ||
        candidate->owner_worker == worker->worker_id) {
      candidate->internal_error = TRUE;
      pthread_mutex_unlock(&shared->state_mutex);
      return -1;
    }

    if (!best ||
        candidate->ready_work > best->ready_work ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count > best->ready_count) ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count == best->ready_count &&
         candidate->owner_index_order < best->owner_index_order) ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count == best->ready_count &&
         candidate->owner_index_order == best->owner_index_order &&
         candidate->owner_worker < best->owner_worker)) {
      best = candidate;
    }
  }

  if (!best) {
    pthread_mutex_unlock(&shared->state_mutex);
    return 1;
  }

  claim->group = best;
  if (best->foreign_reservations_outstanding &&
      shared->helper_foreign_reservations <
          best->foreign_reservations_outstanding) {
    best->internal_error = TRUE;
    pthread_mutex_unlock(&shared->state_mutex);
    return -1;
  }
  if (index_shard_helper_claim_locked(
          shared, best, &claim->task_index)) {
    pthread_mutex_unlock(&shared->state_mutex);
    return -1;
  }
  if (best->foreign_reservations_outstanding) {
    best->foreign_reservations_outstanding--;
    shared->helper_foreign_reservations--;
  }
  best->foreign_work += index_shard_helper_task_work(
      &best->tasks[claim->task_index]);
  best->foreign_claims++;
  shared->helper_tasks_foreign++;
  index_shard_queue_signal_locked(shared);
  pthread_mutex_unlock(&shared->state_mutex);
  return 0;
}

/*
 * Select one inner action while queue_mutex is held. Already prepared
 * computation has priority over completion collection and legacy synchronous
 * packets. New I/O submission and page-plan preparation remain last. No
 * callback is invoked under the queue lock.
 */
int index_shard_inner_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    anbool allow_owner,
    index_shard_inner_claim_t *claim) {
  index_shard_staged_select_class_t select_class;
  int rc;

  if (!worker || !shared || !claim) {
    return -1;
  }
  memset(claim, 0, sizeof(*claim));

  select_class = INDEX_SHARD_STAGED_SELECT_COMPUTE;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, FALSE,
      &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_IO;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, FALSE,
      &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  rc = shared->helper_ready_tasks
      ? index_shard_helper_select_locked(
          worker, shared, &claim->helper)
      : 1;
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_HELPER;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_SUBMIT;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, FALSE,
      &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_PREPARE;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, FALSE,
      &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }
  return 1;
}

static int index_shard_helper_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim,
    index_shard_helper_task_status_t execute_status) {
  index_shard_helper_group_t *group;
  index_shard_helper_task_t *task;
  anbool canonical_ready = FALSE;
  int rc = 0;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;

  if (execute_status != INDEX_SHARD_HELPER_TASK_OK &&
      execute_status != INDEX_SHARD_HELPER_TASK_STOPPED &&
      execute_status != INDEX_SHARD_HELPER_TASK_ERROR) {
    execute_status = INDEX_SHARD_HELPER_TASK_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  task = NULL;
  if (claim->task_index < group->task_count &&
      group->tasks) {
    task = &group->tasks[claim->task_index];
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }

  /*
   * Every returning callback retires one running lifetime lease even when
   * task bookkeeping is already inconsistent. The owner may then cancel
   * READY work and safely wait for the remaining callbacks.
   */
  if (group->running_count) {
    group->running_count--;
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }
  if (group->completed_count < group->task_count) {
    group->completed_count++;
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }

  if (!task || task->scheduler_state !=
      INDEX_SHARD_HELPER_TASK_RUNNING) {
    group->internal_error = TRUE;
    rc = -1;
  } else {
    task->execute_status = execute_status;
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_DONE;
    canonical_ready =
        group->retire &&
        claim->task_index == group->next_retire &&
        execute_status == INDEX_SHARD_HELPER_TASK_OK;
  }

  if (execute_status == INDEX_SHARD_HELPER_TASK_ERROR) {
    group->task_failed = TRUE;
    shared->helper_task_failures++;
  } else if (execute_status ==
             INDEX_SHARD_HELPER_TASK_STOPPED) {
    group->stop_seen = TRUE;
  }
  if (group->internal_error || group->task_failed) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
  } else if (group->stop_seen) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  if (group->completed_count == group->task_count ||
      !group->running_count || group->internal_error ||
      group->task_failed || group->stop_seen ||
      canonical_ready) {
    index_shard_notify_progress_locked(
        shared, group->owner_worker);
  } else if (group->ready_count) {
    index_shard_queue_signal_locked(shared);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

int index_shard_helper_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim) {
  index_shard_helper_group_t *group;
  index_shard_helper_task_t *task;
  index_shard_helper_task_status_t execute_status =
      INDEX_SHARD_HELPER_TASK_ERROR;
  int internal_error = FALSE;
  int rc;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;
  if (!group->ops || !group->ops->execute ||
      !group->tasks || claim->task_index >= group->task_count) {
    internal_error = TRUE;
  } else {
    task = &group->tasks[claim->task_index];
    execute_status = group->ops->execute(
        task->input,
        task->input_bytes,
        task->output,
        task->output_bytes);
  }
  rc = index_shard_helper_complete_claim(
      shared, claim, execute_status);
  if (internal_error) {
    pthread_mutex_lock(&shared->queue_mutex);
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
    index_shard_queue_broadcast_locked(shared);
    pthread_mutex_unlock(&shared->queue_mutex);
    return -1;
  }
  return rc;
}

/* queue_mutex must be held. */
static int index_shard_helper_cancel_for_pool_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  int fatal;
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested ||
      shared->solved_published;
  pthread_mutex_unlock(&shared->state_mutex);

  if (fatal) {
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
    return TRUE;
  }
  if (stopped) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
    return TRUE;
  }
  return FALSE;
}

/* queue_mutex must be held. */
int index_shard_helper_owner_claim_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t *task_index) {
  int fatal;
  int stopped;
  int rc;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested ||
      shared->solved_published;
  if (!fatal && !stopped) {
    rc = index_shard_helper_claim_locked(
        shared, group, task_index);
    pthread_mutex_unlock(&shared->state_mutex);
    return rc;
  }
  pthread_mutex_unlock(&shared->state_mutex);

  if (fatal) {
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
  } else {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  return 1;
}

/*
 * Retire at most one completed canonical task. The queue lock publishes the
 * helper output before the owner callback observes it. The callback runs
 * without pool locks and may mutate only the owning solver.
 */
static int index_shard_helper_retire_one(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  index_shard_helper_task_t *task;
  index_shard_helper_retire_status_t status;
  size_t task_index;

  if (!shared || !group || !group->retire) {
    return 1;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (group->next_retire >= group->task_count) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task_index = group->next_retire;
  task = &group->tasks[task_index];
  if (task->scheduler_state != INDEX_SHARD_HELPER_TASK_DONE) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  if (task->execute_status != INDEX_SHARD_HELPER_TASK_OK) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task->scheduler_state = INDEX_SHARD_HELPER_TASK_RETIRING;
  pthread_mutex_unlock(&shared->queue_mutex);

  status = group->retire(
      task, task_index, group->owner_context);
  if (status != INDEX_SHARD_HELPER_RETIRE_OK &&
      status != INDEX_SHARD_HELPER_RETIRE_STOPPED &&
      status != INDEX_SHARD_HELPER_RETIRE_ERROR) {
    status = INDEX_SHARD_HELPER_RETIRE_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (task->scheduler_state !=
      INDEX_SHARD_HELPER_TASK_RETIRING ||
      group->next_retire != task_index) {
    group->internal_error = TRUE;
    status = INDEX_SHARD_HELPER_RETIRE_ERROR;
  } else {
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_RETIRED;
    group->next_retire++;
  }
  if (status == INDEX_SHARD_HELPER_RETIRE_ERROR) {
    group->task_failed = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group,
        INDEX_SHARD_HELPER_TASK_ERROR);
  } else if (status ==
             INDEX_SHARD_HELPER_RETIRE_STOPPED) {
    group->stop_seen = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group,
        INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  index_shard_notify_progress_locked(
      shared, group->owner_worker);
  pthread_mutex_unlock(&shared->queue_mutex);
  return status == INDEX_SHARD_HELPER_RETIRE_OK ? 0 : -1;
}

/*
 * A staged owner may temporarily expose one synchronous helper group while its
 * owner callback is executing. The staged task remains the lifetime owner; the
 * child group receives only immutable inputs and disjoint output ranges.
 */
anbool index_shard_helper_staged_child_allowed(
    const index_shard_worker_context_t *ctx) {
  return ctx &&
      ctx->staged_owner_callback_active &&
      ctx->staged_owner_callback_group &&
      ctx->published_staged_group ==
          ctx->staged_owner_callback_group &&
      !ctx->published_helper_group;
}

size_t index_shard_helper_available_workers(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t available = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      (!ctx->published_staged_group ||
       index_shard_helper_staged_child_allowed(ctx)) &&
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active) {
    available = index_shard_helper_idle_workers_locked(
        shared);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return available;
}

size_t index_shard_helper_prepare_reserve(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t available = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active) {
    available = index_shard_helper_idle_workers_locked(shared);
    if (available) {
      shared->helper_preparations_active++;
      ctx->helper_preparation_active = TRUE;
      ctx->helper_preparation_generation =
          ctx->generation_seen;
      ctx->helper_preparation_index_order =
          ctx->current_index_order;
      ctx->helper_preparation_workers = available;
    }
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return available;
}
/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/*
 * Published helper and staged-run APIs for owner workers.
 *
 * Owner-facing helper publication remains in this module.
 */

/* queue_mutex must be held. */
static int index_shard_helper_prepare_clear_locked(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared) {
  int invalid = FALSE;

  if (!ctx || !shared || !ctx->helper_preparation_active) {
    return 0;
  }
  if (shared->helper_preparations_active != 1U) {
    logerr("[index-shard] invalid helper preparation count=%zu\n",
           shared->helper_preparations_active);
    invalid = TRUE;
  } else {
    shared->helper_preparations_active--;
  }
  ctx->helper_preparation_active = FALSE;
  ctx->helper_preparation_generation = 0U;
  ctx->helper_preparation_index_order = SIZE_MAX;
  ctx->helper_preparation_workers = 0U;
  index_shard_queue_signal_locked(shared);
  return invalid ? -1 : 0;
}

void index_shard_helper_prepare_cancel(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;

  if (!ctx || !ctx->pool) {
    return;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (index_shard_helper_prepare_clear_locked(
          ctx, shared)) {
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
}

static index_shard_helper_run_status_t
index_shard_helper_run_internal(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  index_shard_helper_group_t group;
  index_shard_helper_claim_t claim;
  index_shard_helper_run_status_t result;
  size_t i;
  int have_claim = TRUE;
  int wait_broken = FALSE;
  int fatal_requested = FALSE;
  int prepublish_fatal = FALSE;
  int preparation_permit = FALSE;
  int staged_child_permit = FALSE;

  if (stats) {
    memset(stats, 0, sizeof(*stats));
  }
  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 || task_count < 2U) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (!ops || !ops->execute || !tasks ||
      task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }

  memset(&group, 0, sizeof(group));
  group.ops = ops;
  group.tasks = tasks;
  group.task_count = task_count;
  group.retire = retire;
  group.owner_context = owner_context;
  group.generation = ctx->generation_seen;
  group.owner_epoch = ++ctx->helper_group_epoch;
  group.owner_worker = ctx->worker_id;
  group.owner_index_order = ctx->current_index_order;
  group.verification_group = FALSE;

  for (i = 0U; i < task_count; i++) {
    unsigned long long work;

    if ((!tasks[i].input && tasks[i].input_bytes) ||
        (!tasks[i].output && tasks[i].output_bytes)) {
      index_shard_helper_prepare_cancel();
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
    tasks[i].scheduler_state = INDEX_SHARD_HELPER_TASK_READY;
    tasks[i].execute_status = INDEX_SHARD_HELPER_TASK_ERROR;
    work = index_shard_helper_task_work(&tasks[i]);
    if (index_shard_helper_add_work(
            &group.ready_work, work)) {
      index_shard_helper_prepare_cancel();
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
  }
  group.ready_count = task_count;

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  staged_child_permit =
      index_shard_helper_staged_child_allowed(ctx);
  preparation_permit =
      ctx->helper_preparation_active &&
      ctx->helper_preparation_generation ==
          ctx->generation_seen &&
      ctx->helper_preparation_index_order ==
          ctx->current_index_order;
  if (ctx->generation_seen != ctx->pool->generation ||
      ctx->published_helper_group ||
      (ctx->published_staged_group &&
       !staged_child_permit) ||
      (!preparation_permit &&
       shared->helper_preparations_active)) {
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      return INDEX_SHARD_HELPER_FATAL;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  group.foreign_reserve =
      index_shard_helper_idle_workers_locked(shared);
  if (preparation_permit &&
      group.foreign_reserve >
          ctx->helper_preparation_workers) {
    group.foreign_reserve =
        ctx->helper_preparation_workers;
  }
  if (group.foreign_reserve > task_count - 1U) {
    group.foreign_reserve = task_count - 1U;
  }
  if (!group.foreign_reserve) {
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      return INDEX_SHARD_HELPER_FATAL;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (SIZE_MAX - shared->helper_foreign_reservations <
      group.foreign_reserve) {
    (void)index_shard_helper_prepare_clear_locked(
        ctx, shared);
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return INDEX_SHARD_HELPER_FATAL;
  }

  memset(&claim, 0, sizeof(claim));
  claim.group = &group;
  result = INDEX_SHARD_HELPER_OK;
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->fatal_error) {
    prepublish_fatal = TRUE;
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (shared->stop_requested ||
             shared->solved_published) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else if (index_shard_helper_claim_locked(
                 shared, &group, &claim.task_index)) {
    prepublish_fatal = TRUE;
    result = INDEX_SHARD_HELPER_FATAL;
  }
  if (result != INDEX_SHARD_HELPER_OK) {
    pthread_mutex_unlock(&shared->state_mutex);
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      result = INDEX_SHARD_HELPER_FATAL;
      prepublish_fatal = TRUE;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    if (prepublish_fatal) {
      index_shard_request_fatal_stop(shared);
    }
    return result;
  }
  group.owner_claims++;
  group.owner_work += index_shard_helper_task_work(
      &group.tasks[claim.task_index]);
  shared->helper_tasks_owner++;
  if (index_shard_helper_prepare_clear_locked(
          ctx, shared)) {
    pthread_mutex_unlock(&shared->state_mutex);
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return INDEX_SHARD_HELPER_FATAL;
  }
  if (index_shard_helper_account_ready_locked(
          shared, &group)) {
    pthread_mutex_unlock(&shared->state_mutex);
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return INDEX_SHARD_HELPER_FATAL;
  }
  group.foreign_reservations_outstanding =
      group.foreign_reserve;
  shared->helper_foreign_reservations +=
      group.foreign_reservations_outstanding;
  ctx->published_helper_group = &group;
  shared->helper_groups_active++;
  shared->helper_groups_published++;
  if (shared->observability_enabled) {
    group.helper_window_start = monotonic_seconds();
    group.verification_group =
        ops->name && !strcmp(ops->name, "verify-context");
    index_shard_observability_increment(
        &shared->helper_windows_started);
    if (group.verification_group) {
      index_shard_observability_increment(
          &shared->verification_helper_groups);
      index_shard_observability_add(
          &shared->verification_helper_contexts,
          (unsigned long long)task_count);
    }
  }
  pthread_mutex_unlock(&shared->state_mutex);
  index_shard_queue_signal_locked(shared);
  pthread_mutex_unlock(&shared->queue_mutex);
  fitsbin_payload_io_notify_wait_helpers();

  while (1) {
    if (have_claim) {
      (void)index_shard_helper_execute_claim(
          shared, &claim);

      have_claim = FALSE;
    }

    while (!index_shard_helper_retire_one(
               shared, &group)) {
      /* Retire every currently completed canonical packet. */
    }

    pthread_mutex_lock(&shared->queue_mutex);
    (void)index_shard_helper_cancel_for_pool_locked(
        shared, &group);
    if (group.internal_error || group.task_failed) {
      index_shard_helper_cancel_ready_locked(
          shared,
          &group, INDEX_SHARD_HELPER_TASK_ERROR);
    } else if (group.stop_seen) {
      index_shard_helper_cancel_ready_locked(
          shared,
          &group, INDEX_SHARD_HELPER_TASK_STOPPED);
    }

    if (group.internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }

    /*
     * As with staged groups, a foreign task can finish after the owner
     * retirement scan but before this lock is acquired. Consume the
     * canonical completed task before deciding that the group is quiescent
     * or waiting for another notification.
     */
    if (group.retire &&
        !group.task_failed && !group.stop_seen &&
        group.next_retire < group.task_count &&
        group.tasks[group.next_retire].scheduler_state ==
            INDEX_SHARD_HELPER_TASK_DONE &&
        group.tasks[group.next_retire].execute_status ==
            INDEX_SHARD_HELPER_TASK_OK) {
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (!group.ready_count && !group.running_count) {
      break;
    }

    if (!group.ready_count) {
      if (wait_broken) {
        struct timespec pause = { 0, 1000000L };

        pthread_mutex_unlock(&shared->queue_mutex);
        nanosleep(&pause, NULL);
        continue;
      }

      {
        int wait_status;
        double wait_start;

        shared->helper_owner_wait_calls++;
        wait_start = monotonic_seconds();
        if (!ctx->owner_cv_ready || ctx->owner_waiting ||
            ctx->owner_wake_pending) {
          wait_status = EINVAL;
        } else {
          ctx->owner_waiting = TRUE;
          wait_status = pthread_cond_wait(
              &ctx->owner_cv, &shared->queue_mutex);
          ctx->owner_waiting = FALSE;
          ctx->owner_wake_pending = FALSE;
        }
        shared->helper_owner_wait_seconds +=
            monotonic_seconds() - wait_start;
        if (wait_status) {
          group.internal_error = TRUE;
          index_shard_helper_cancel_ready_locked(
              shared,
              &group, INDEX_SHARD_HELPER_TASK_ERROR);
          wait_broken = TRUE;
        }
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (group.foreign_reserve > group.foreign_claims) {
      size_t outstanding =
          group.foreign_reserve - group.foreign_claims;

      if (group.ready_count <= outstanding &&
          !group.owner_reserve_yielded) {
        group.owner_reserve_yielded = TRUE;
        index_shard_queue_signal_locked(shared);
        pthread_mutex_unlock(&shared->queue_mutex);
        sched_yield();
        continue;
      }
      if (group.owner_reserve_yielded) {
        if (index_shard_helper_release_foreign_reservations_locked(
                shared,
                &group,
                group.foreign_reservations_outstanding)) {
          index_shard_helper_cancel_ready_locked(
              shared,
              &group,
              INDEX_SHARD_HELPER_TASK_ERROR);
        }
        group.foreign_reserve = group.foreign_claims;
      }
    }

    claim.group = &group;
    {
      int claim_status =
          index_shard_helper_owner_claim_locked(
              shared, &group, &claim.task_index);

      if (claim_status < 0) {
        group.internal_error = TRUE;
        index_shard_helper_cancel_ready_locked(
            shared,
            &group, INDEX_SHARD_HELPER_TASK_ERROR);
      } else if (!claim_status) {
        group.owner_claims++;
        group.owner_work += index_shard_helper_task_work(
            &group.tasks[claim.task_index]);
        shared->helper_tasks_owner++;
        have_claim = TRUE;
        pthread_mutex_unlock(&shared->queue_mutex);
        continue;
      }
    }
    if (group.internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
  }

  if (index_shard_helper_release_foreign_reservations_locked(
          shared,
          &group,
          group.foreign_reservations_outstanding)) {
    group.internal_error = TRUE;
  }
  if (group.ready_count || group.running_count ||
      group.completed_count != group.task_count) {
    group.internal_error = TRUE;
  }
  if (group.retire &&
      !group.task_failed && !group.stop_seen &&
      group.next_retire != group.task_count) {
    group.internal_error = TRUE;
  }
  {
    size_t lanes_cleared = 0U;

    for (i = 0U; i < (size_t)shared->worker_count; i++) {
      if (ctx->pool->contexts[i].published_helper_group ==
          &group) {
        ctx->pool->contexts[i].published_helper_group = NULL;
        lanes_cleared++;
      }
    }
    if (lanes_cleared != 1U) {
      group.internal_error = TRUE;
    } else {
      shared->helper_groups_completed++;
    }
  }
  if (ctx->published_helper_group == &group) {
    ctx->published_helper_group = NULL;
    group.internal_error = TRUE;
  }
  if (index_shard_helper_unaccount_ready_locked(
          shared, &group)) {
    group.internal_error = TRUE;
  }
  if (!shared->helper_groups_active) {
    group.internal_error = TRUE;
  } else {
    shared->helper_groups_active--;
  }
  index_shard_notify_progress_locked(
      shared, group.owner_worker);
  pthread_mutex_unlock(&shared->queue_mutex);

  if (shared->observability_enabled) {
    double helper_window_seconds = 0.0;

    helper_window_seconds =
        monotonic_seconds() - group.helper_window_start;
    if (helper_window_seconds >= 0.0) {
      pthread_mutex_lock(&shared->queue_mutex);
      shared->helper_window_seconds +=
          helper_window_seconds;
      pthread_mutex_unlock(&shared->queue_mutex);
    }
  }
  if (stats) {
    stats->owner_tasks = group.owner_claims;
    stats->foreign_tasks = group.foreign_claims;
    stats->max_concurrent_tasks = group.max_running;
    stats->owner_work_units = group.owner_work;
    stats->foreign_work_units = group.foreign_work;
  }
  if (group.internal_error) {
    if (!fatal_requested) {
      index_shard_request_fatal_stop(shared);
    }
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (group.task_failed) {
    result = INDEX_SHARD_HELPER_TASK_FAILED;
  } else if (group.stop_seen) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  return result;
}

index_shard_helper_run_status_t
index_shard_helper_run(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_run_stats_t *stats) {
  return index_shard_helper_run_internal(
      ops, tasks, task_count, NULL, NULL, stats);
}

index_shard_helper_run_status_t
index_shard_helper_run_ordered(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats) {
  if (!retire) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  return index_shard_helper_run_internal(
      ops,
      tasks,
      task_count,
      retire,
      owner_context,
      stats);
}

#ifdef TESTING_INDEX_SHARD_HELPER_LEDGER
static void index_shard_helper_ledger_test_init_group(
    index_shard_helper_group_t *group,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    const index_shard_pool_t *pool,
    const index_shard_worker_context_t *owner) {
  size_t i;

  memset(group, 0, sizeof(*group));
  memset(tasks, 0, task_count * sizeof(*tasks));
  group->tasks = tasks;
  group->task_count = task_count;
  group->generation = pool->generation;
  group->owner_epoch = owner->helper_group_epoch;
  group->owner_worker = owner->worker_id;
  group->owner_index_order = owner->current_index_order;
  group->ready_count = task_count;
  group->ready_work = (unsigned long long)task_count;
  for (i = 0U; i < task_count; i++) {
    tasks[i].work_units = 1U;
    tasks[i].scheduler_state = INDEX_SHARD_HELPER_TASK_READY;
    tasks[i].execute_status = INDEX_SHARD_HELPER_TASK_ERROR;
  }
}

int index_shard_helper_readiness_ledger_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_helper_group_t group;
  index_shard_helper_group_t mismatch_group;
  index_shard_helper_task_t tasks[4];
  index_shard_helper_task_t mismatch_tasks[2];
  index_shard_helper_claim_t first_owner_claim;
  index_shard_helper_claim_t foreign_claim;
  index_shard_helper_claim_t second_owner_claim;
  index_shard_thread_state_t *shared;
  int failures = 0;
  size_t i;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(&first_owner_claim, 0, sizeof(first_owner_claim));
  memset(&foreign_claim, 0, sizeof(foreign_claim));
  memset(&second_owner_claim, 0, sizeof(second_owner_claim));
  shared = &pool.shared;
  if (pthread_mutex_init(&shared->queue_mutex, NULL)) {
    return 1;
  }
  if (pthread_mutex_init(&shared->state_mutex, NULL)) {
    pthread_mutex_destroy(&shared->queue_mutex);
    return 1;
  }

  pool.worker_count = 2;
  pool.generation = 31U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  for (i = 0U; i < 2U; i++) {
    contexts[i].worker_id = (int)i;
    contexts[i].pool = &pool;
    contexts[i].generation_seen = pool.generation;
    contexts[i].current_outer_active = TRUE;
    contexts[i].current_index_order = i;
  }
  contexts[0].helper_group_epoch = 41U;
  index_shard_helper_ledger_test_init_group(
      &group, tasks, 4U, &pool, &contexts[0]);
  first_owner_claim.group = &group;
  second_owner_claim.group = &group;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_helper_owner_claim_locked(
      shared, &group, &first_owner_claim.task_index) != 0;
  failures += first_owner_claim.task_index != 0U;
  failures += group.ready_count != 3U;
  failures += shared->helper_ready_tasks != 0U;
  failures += index_shard_helper_account_ready_locked(
      shared, &group) != 0;
  failures += !group.ready_accounted;
  failures += shared->helper_ready_tasks != 3U;
  contexts[0].published_helper_group = &group;

  failures += index_shard_helper_select_locked(
      &contexts[1], shared, &foreign_claim) != 0;
  failures += foreign_claim.group != &group;
  failures += foreign_claim.task_index != 1U;
  failures += group.ready_count != 2U;
  failures += shared->helper_ready_tasks != 2U;

  failures += index_shard_helper_owner_claim_locked(
      shared, &group, &second_owner_claim.task_index) != 0;
  failures += second_owner_claim.task_index != 2U;
  failures += group.ready_count != 1U;
  failures += group.running_count != 3U;
  failures += shared->helper_ready_tasks != 1U;
  pthread_mutex_unlock(&shared->queue_mutex);

  failures += index_shard_helper_complete_claim(
      shared, &first_owner_claim,
      INDEX_SHARD_HELPER_TASK_OK) != 0;
  failures += index_shard_helper_complete_claim(
      shared, &foreign_claim,
      INDEX_SHARD_HELPER_TASK_OK) != 0;
  failures += index_shard_helper_complete_claim(
      shared, &second_owner_claim,
      INDEX_SHARD_HELPER_TASK_OK) != 0;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += group.ready_count != 1U;
  failures += group.running_count != 0U;
  failures += group.completed_count != 3U;
  failures += shared->helper_ready_tasks != 1U;
  index_shard_helper_cancel_ready_locked(
      shared, &group, INDEX_SHARD_HELPER_TASK_STOPPED);
  failures += group.ready_count != 0U;
  failures += group.running_count != 0U;
  failures += group.completed_count != 4U;
  failures += group.ready_work != 0U;
  failures += !group.stop_seen;
  failures += shared->helper_ready_tasks != 0U;
  for (i = 0U; i < 4U; i++) {
    failures += tasks[i].scheduler_state !=
        INDEX_SHARD_HELPER_TASK_DONE;
  }
  contexts[0].published_helper_group = NULL;
  failures += index_shard_helper_unaccount_ready_locked(
      shared, &group) != 0;
  failures += group.ready_accounted;
  failures += shared->helper_ready_tasks != 0U;
  pthread_mutex_unlock(&shared->queue_mutex);

  contexts[0].helper_group_epoch++;
  index_shard_helper_ledger_test_init_group(
      &mismatch_group, mismatch_tasks, 2U,
      &pool, &contexts[0]);
  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_helper_account_ready_locked(
      shared, &mismatch_group) != 0;
  contexts[0].published_helper_group = &mismatch_group;
  failures += shared->helper_ready_tasks != 2U;
  shared->helper_ready_tasks = 1U;
  pthread_mutex_lock(&shared->state_mutex);
  shared->fatal_error = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);
  failures += !index_shard_helper_cancel_for_pool_locked(
      shared, &mismatch_group);
  failures += !mismatch_group.internal_error;
  failures += !mismatch_group.task_failed;
  failures += mismatch_group.stop_seen;
  failures += mismatch_group.ready_count != 0U;
  failures += mismatch_group.running_count != 0U;
  failures += mismatch_group.completed_count != 2U;
  failures += mismatch_group.ready_work != 0U;
  for (i = 0U; i < 2U; i++) {
    failures += mismatch_tasks[i].scheduler_state !=
        INDEX_SHARD_HELPER_TASK_DONE;
    failures += mismatch_tasks[i].execute_status !=
        INDEX_SHARD_HELPER_TASK_ERROR;
  }
  contexts[0].published_helper_group = NULL;
  failures += index_shard_helper_unaccount_ready_locked(
      shared, &mismatch_group) != 0;
  failures += mismatch_group.ready_accounted;
  failures += shared->helper_ready_tasks != 1U;
  shared->helper_ready_tasks = 0U;
  pthread_mutex_unlock(&shared->queue_mutex);

  pthread_mutex_destroy(&shared->state_mutex);
  pthread_mutex_destroy(&shared->queue_mutex);
  return failures;
}
#endif
