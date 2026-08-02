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
#ifdef TESTING_INDEX_SHARD_STAGED
/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/*
 * Deterministic tests for the private staged scheduler state machine.
 *
 * This translation unit is compiled only by the TESTING_INDEX_SHARD_STAGED
 * build, so production visibility and state layout remain unchanged.
 */

typedef struct index_shard_staged_retire_test_context {
  index_shard_thread_state_t *shared;
  size_t order[8];
  size_t calls;
  size_t task_zero_calls;
  anbool stop_during_retire;
} index_shard_staged_retire_test_context_t;

index_shard_staged_retire_status_t
index_shard_staged_retire_test_callback(
    const index_shard_staged_task_t *task,
    size_t task_index,
    void *opaque) {
  index_shard_staged_retire_test_context_t *context = opaque;

  if (!task || !context ||
      context->calls >= sizeof(context->order) /
          sizeof(context->order[0])) {
    return INDEX_SHARD_STAGED_RETIRE_ERROR;
  }
  context->order[context->calls++] = task_index;
  if (context->stop_during_retire) {
    pthread_mutex_lock(&context->shared->state_mutex);
    context->shared->stop_requested = TRUE;
    pthread_mutex_unlock(&context->shared->state_mutex);
    return INDEX_SHARD_STAGED_RETIRE_MORE;
  }
  if (task_index == 0U && context->task_zero_calls++ < 2U) {
    return INDEX_SHARD_STAGED_RETIRE_MORE;
  }
  return INDEX_SHARD_STAGED_RETIRE_OK;
}

int index_shard_staged_retire_test_init(
    index_shard_thread_state_t *shared) {
  int rc;

  memset(shared, 0, sizeof(*shared));
  rc = pthread_mutex_init(&shared->queue_mutex, NULL);
  if (rc) {
    return -1;
  }
  rc = pthread_mutex_init(&shared->state_mutex, NULL);
  if (rc) {
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  rc = pthread_cond_init(&shared->queue_cv, NULL);
  if (rc) {
    pthread_mutex_destroy(&shared->state_mutex);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  shared->worker_count = 2;
  if (index_shard_completion_registry_init(shared, 2)) {
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->state_mutex);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  return 0;
}

void index_shard_staged_retire_test_destroy(
    index_shard_thread_state_t *shared) {
  index_shard_completion_registry_destroy(shared);
  pthread_cond_destroy(&shared->queue_cv);
  pthread_mutex_destroy(&shared->state_mutex);
  pthread_mutex_destroy(&shared->queue_mutex);
}

int index_shard_staged_retire_test_ready(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task) {
  int rc;

  pthread_mutex_lock(&shared->queue_mutex);
  rc = index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_RESULTS_READY);
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

int index_shard_staged_retire_test_order(void) {
  index_shard_thread_state_t shared;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_staged_retire_test_context_t context;
  int failures = 0;

  if (index_shard_staged_retire_test_init(&shared)) {
    return 1;
  }
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&context, 0, sizeof(context));
  group.tasks = tasks;
  group.task_count = 2U;
  group.retire = index_shard_staged_retire_test_callback;
  group.owner_context = &context;
  context.shared = &shared;
  pthread_mutex_lock(&shared.queue_mutex);
  failures += index_shard_staged_set_state_locked(
      &shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_PREPARE_READY) != 0;
  failures += index_shard_staged_set_state_locked(
      &shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_RESULTS_READY) != 0;
  pthread_mutex_unlock(&shared.queue_mutex);

  failures += index_shard_staged_retire_one(&shared, &group) != 1;
  failures += context.calls != 0U;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_PREPARE_READY;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_PREPARE_READY;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RETIRED;
  failures += group.next_retire != 1U;

  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RETIRED;
  failures += group.next_retire != 2U;
  failures += context.calls != 4U;
  failures += context.order[0] != 0U;
  failures += context.order[1] != 0U;
  failures += context.order[2] != 0U;
  failures += context.order[3] != 1U;
  failures += group.reorder_ready != 0U;
  failures += shared.staged_reorder_ready != 0U;
  failures += group.internal_error;
  failures += group.task_failed;

  index_shard_staged_retire_test_destroy(&shared);
  return failures;
}

int index_shard_staged_retire_test_terminal(void) {
  index_shard_thread_state_t shared;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_staged_retire_test_context_t context;
  int failures = 0;

  if (index_shard_staged_retire_test_init(&shared)) {
    return 1;
  }
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&context, 0, sizeof(context));
  group.tasks = tasks;
  group.task_count = 2U;
  group.retire = index_shard_staged_retire_test_callback;
  group.owner_context = &context;
  context.shared = &shared;
  context.stop_during_retire = TRUE;
  pthread_mutex_lock(&shared.queue_mutex);
  failures += index_shard_staged_set_state_locked(
      &shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_RESULTS_READY) != 0;
  failures += index_shard_staged_set_state_locked(
      &shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  pthread_mutex_unlock(&shared.queue_mutex);

  failures += index_shard_staged_retire_one(&shared, &group) != -1;
  failures += context.calls != 1U;
  failures += context.order[0] != 0U;
  failures += group.next_retire != 0U;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_STOPPED;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_STOPPED;
  failures += !group.cancelling;
  failures += !group.stop_seen;
  failures += group.reorder_ready != 0U;
  failures += group.compute_ready != 0U;
  failures += shared.staged_reorder_ready != 0U;
  failures += shared.staged_compute_ready != 0U;
  failures += group.internal_error;
  failures += group.task_failed;

  index_shard_staged_retire_test_destroy(&shared);
  return failures;
}

int index_shard_staged_child_helper_test(void) {
  index_shard_worker_context_t context;
  index_shard_staged_group_t group;
  index_shard_helper_group_t helper;
  int failures = 0;

  memset(&context, 0, sizeof(context));
  memset(&group, 0, sizeof(group));
  memset(&helper, 0, sizeof(helper));

  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.published_staged_group = &group;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.staged_owner_callback_active = TRUE;
  context.staged_owner_callback_group = &group;
  failures += index_shard_helper_staged_child_allowed(
      &context) != TRUE;
  context.published_helper_group = &helper;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.published_helper_group = NULL;
  context.staged_owner_callback_group = NULL;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  return failures;
}

int index_shard_observability_counter_test(void) {
  index_shard_thread_state_t shared;
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  unsigned char outer_states[2] = {
      INDEX_SHARD_OUTER_UNCLAIMED,
      INDEX_SHARD_OUTER_UNCLAIMED
  };
  fitsbin_mmap_advice_t mmap_advice = FITSBIN_MMAP_ADVICE_NORMAL;
  size_t index_order = SIZE_MAX;
  int failures = 0;

  if (index_shard_staged_retire_test_init(&shared)) {
    return 1;
  }
  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  pool.contexts = contexts;
  pool.worker_count = 2;
  shared.pool = &pool;
  shared.worker_count = 2;
  contexts[0].owner_cv_ready =
      pthread_cond_init(&contexts[0].owner_cv, NULL) == 0;
  if (!contexts[0].owner_cv_ready) {
    index_shard_staged_retire_test_destroy(&shared);
    return 1;
  }

  shared.observability_enabled = TRUE;
  shared.outer_states = outer_states;
  shared.nindexes = 2U;
  shared.outer_unclaimed = 2U;
  shared.producer_width = 2U;
  pthread_mutex_lock(&shared.queue_mutex);
  shared.queue_waiters = 1U;
  index_shard_queue_signal_locked(&shared);
  failures += shared.queue_signals != 1U;
  failures += index_shard_claim_outer_locked(
      &contexts[0], &shared, 0U,
      &index_order, &mmap_advice) != 0;
  failures += index_order != 0U;
  failures += shared.queue_signals != 2U;
  failures += shared.outer_running != 1U;
  failures += shared.outer_unclaimed != 1U;
  contexts[0].owner_waiting = TRUE;
  index_shard_owner_signal_locked(&shared, 0);
  failures += shared.owner_signals != 1U;
  index_shard_owner_signal_locked(&shared, 0);
  failures += shared.owner_signals != 1U;
  failures += shared.owner_signals_coalesced != 1U;
  contexts[0].owner_wake_pending = FALSE;
  index_shard_queue_broadcast_locked(&shared);
  failures += shared.queue_broadcasts != 1U;
  failures += shared.owner_broadcasts != 1U;
  contexts[0].owner_waiting = FALSE;
  contexts[0].owner_wake_pending = FALSE;
  shared.queue_waiters = 0U;
  index_shard_queue_broadcast_locked(&shared);
  failures += shared.queue_broadcasts != 1U;
  failures += shared.owner_broadcasts != 1U;
  pthread_mutex_unlock(&shared.queue_mutex);

  pthread_cond_destroy(&contexts[0].owner_cv);
  index_shard_staged_retire_test_destroy(&shared);
  return failures;
}

int index_shard_staged_mask_selection_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_staged_group_t groups[2];
  index_shard_staged_task_t tasks[2][2];
  index_shard_staged_claim_t claim;
  index_shard_thread_state_t *shared;
  int failures = 0;
  int i;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(groups, 0, sizeof(groups));
  memset(tasks, 0, sizeof(tasks));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 2;
  pool.generation = 7U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  shared->observability_enabled = TRUE;

  for (i = 0; i < 2; i++) {
    contexts[i].worker_id = i;
    contexts[i].pool = &pool;
    contexts[i].generation_seen = pool.generation;
    contexts[i].current_outer_active = TRUE;
    contexts[i].current_index_order = i ? 3U : 5U;
    contexts[i].staged_group_epoch = 11U + (unsigned long long)i;
    contexts[i].published_staged_group = &groups[i];
    groups[i].pool = &pool;
    groups[i].tasks = tasks[i];
    groups[i].task_count = 2U;
    groups[i].generation = pool.generation;
    groups[i].owner_epoch = contexts[i].staged_group_epoch;
    groups[i].owner_worker = i;
    groups[i].owner_index_order = contexts[i].current_index_order;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &groups[0], &tasks[0][0],
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &groups[1], &tasks[1][1],
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  failures += index_shard_staged_select_locked(
      &contexts[0], shared,
      INDEX_SHARD_STAGED_SELECT_COMPUTE,
      FALSE, TRUE, &claim) != 0;
  failures += claim.group != &groups[1];
  failures += claim.task_index != 1U;
  failures += claim.kind != INDEX_SHARD_STAGED_CLAIM_EXECUTE;
  failures += shared->selection_tasks_scanned != 0U;
  failures += groups[1].compute_ready_mask != UINT64_C(0);
  failures += tasks[1][1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_EXECUTING;
  pthread_mutex_unlock(&shared->queue_mutex);

  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_ready_before_outer_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_thread_state_t *shared;
  index_shard_inner_claim_t claim;
  unsigned char outer_states[3] = {
      INDEX_SHARD_OUTER_RUNNING,
      INDEX_SHARD_OUTER_UNCLAIMED,
      INDEX_SHARD_OUTER_UNCLAIMED
  };
  fitsbin_mmap_advice_t mmap_advice = FITSBIN_MMAP_ADVICE_NORMAL;
  size_t index_order = SIZE_MAX;
  index_shard_work_selection_t selection;
  int failures = 0;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&claim, 0, sizeof(claim));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }

  pool.worker_count = 2;
  pool.generation = 8U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  shared->observability_enabled = TRUE;
  shared->outer_states = outer_states;
  shared->nindexes = 3U;
  shared->outer_unclaimed = 2U;
  shared->outer_running = 1U;
  shared->producer_width = 2U;
  shared->canonical_scan_cursor = 1U;

  contexts[0].worker_id = 0;
  contexts[0].pool = &pool;
  contexts[0].generation_seen = pool.generation;
  contexts[0].current_outer_active = TRUE;
  contexts[0].current_index_order = 0U;
  contexts[0].staged_group_epoch = 31U;
  contexts[0].published_staged_group = &group;
  contexts[1].worker_id = 1;
  contexts[1].pool = &pool;
  contexts[1].generation_seen = pool.generation;

  group.pool = &pool;
  group.tasks = tasks;
  group.task_count = 2U;
  group.generation = pool.generation;
  group.owner_epoch = contexts[0].staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = 0U;

  /* No ready computation: the canonical outer task wins immediately. */
  contexts[1].ready_before_outer_eligible = TRUE;
  selection = index_shard_select_work(
      &contexts[1], shared, &index_order, &mmap_advice, &claim);
  failures += selection != INDEX_SHARD_WORK_OUTER;
  failures += index_order != 1U;
  failures += contexts[1].ready_before_outer_eligible;
  failures += outer_states[1] != INDEX_SHARD_OUTER_RUNNING;
  failures += shared->staged_ready_before_outer_claims != 0U;

  /* Model completion of that outer task and publish two foreign packets. */
  pthread_mutex_lock(&shared->queue_mutex);
  outer_states[1] = INDEX_SHARD_OUTER_FINISHED;
  shared->outer_running--;
  contexts[1].ready_before_outer_eligible = TRUE;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  pthread_mutex_unlock(&shared->queue_mutex);

  memset(&claim, 0, sizeof(claim));
  index_order = SIZE_MAX;
  selection = index_shard_select_work(
      &contexts[1], shared, &index_order, &mmap_advice, &claim);
  failures += selection != INDEX_SHARD_WORK_HELPER;
  failures += claim.kind != INDEX_SHARD_INNER_CLAIM_STAGED;
  failures += claim.staged.group != &group;
  failures += claim.staged.task_index != 0U;
  failures += claim.staged.owner_claim;
  failures += contexts[1].ready_before_outer_eligible;
  failures += outer_states[2] != INDEX_SHARD_OUTER_UNCLAIMED;
  failures += shared->staged_ready_before_outer_claims != 1U;

  failures += index_shard_staged_complete_claim(
      shared, &claim.staged,
      INDEX_SHARD_STAGED_EXECUTE_OK, 0.0, 0ULL) != 0;

  /* The credit was consumed: the next outer task beats packet one. */
  memset(&claim, 0, sizeof(claim));
  index_order = SIZE_MAX;
  selection = index_shard_select_work(
      &contexts[1], shared, &index_order, &mmap_advice, &claim);
  failures += selection != INDEX_SHARD_WORK_OUTER;
  failures += index_order != 2U;
  failures += outer_states[2] != INDEX_SHARD_OUTER_RUNNING;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY;
  failures += shared->staged_ready_before_outer_claims != 1U;

  /* Terminal state always wins without consuming or claiming work. */
  contexts[1].ready_before_outer_eligible = TRUE;
  pthread_mutex_lock(&shared->state_mutex);
  shared->stop_requested = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);
  selection = index_shard_select_work(
      &contexts[1], shared, &index_order, &mmap_advice, &claim);
  failures += selection != INDEX_SHARD_WORK_DONE;
  failures += !contexts[1].ready_before_outer_eligible;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_STOPPED) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_STOPPED) != 0;
  pthread_mutex_unlock(&shared->queue_mutex);
  contexts[0].published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_outer_cap_liveness_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[4];
  index_shard_staged_group_t group;
  index_shard_staged_task_t task;
  index_shard_thread_state_t *shared;
  index_shard_inner_claim_t claim;
  unsigned char outer_states[4] = {
      INDEX_SHARD_OUTER_RUNNING,
      INDEX_SHARD_OUTER_RUNNING,
      INDEX_SHARD_OUTER_UNCLAIMED,
      INDEX_SHARD_OUTER_UNCLAIMED
  };
  fitsbin_mmap_advice_t mmap_advice = FITSBIN_MMAP_ADVICE_NORMAL;
  size_t index_order = SIZE_MAX;
  index_shard_work_selection_t selection;
  int failures = 0;
  int i;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(&group, 0, sizeof(group));
  memset(&task, 0, sizeof(task));
  memset(&claim, 0, sizeof(claim));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }

  pool.worker_count = 4;
  pool.generation = 9U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 4;
  shared->observability_enabled = TRUE;
  shared->outer_states = outer_states;
  shared->nindexes = 4U;
  shared->outer_unclaimed = 2U;
  shared->outer_running = 2U;
  shared->producer_width = 2U;
  shared->canonical_scan_cursor = 2U;
  for (i = 0; i < 4; i++) {
    contexts[i].worker_id = i;
    contexts[i].pool = &pool;
    contexts[i].generation_seen = pool.generation;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  shared->queue_waiters = 1U;
  index_shard_queue_signal_locked(shared);
  failures += shared->queue_signals != 0U;
  failures += shared->queue_signals_no_work != 1U;

  outer_states[1] = INDEX_SHARD_OUTER_FINISHED;
  shared->outer_running--;
  index_shard_queue_signal_locked(shared);
  failures += shared->queue_signals != 1U;
  failures += index_shard_claim_outer_locked(
      &contexts[2], shared, 2U,
      &index_order, &mmap_advice) != 0;
  failures += index_order != 2U;
  failures += shared->outer_running != 2U;
  failures += shared->outer_unclaimed != 1U;

  contexts[0].current_outer_active = TRUE;
  contexts[0].current_index_order = 0U;
  contexts[0].staged_group_epoch = 41U;
  contexts[0].published_staged_group = &group;
  group.pool = &pool;
  group.tasks = &task;
  group.task_count = 1U;
  group.generation = pool.generation;
  group.owner_epoch = contexts[0].staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = 0U;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &task,
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY) != 0;
  pthread_mutex_unlock(&shared->queue_mutex);

  selection = index_shard_select_work(
      &contexts[3], shared, &index_order, &mmap_advice, &claim);
  failures += selection != INDEX_SHARD_WORK_HELPER;
  failures += claim.kind != INDEX_SHARD_INNER_CLAIM_STAGED;
  failures += claim.staged.group != &group;
  failures += claim.staged.task_index != 0U;
  failures += claim.staged.owner_claim;
  failures += outer_states[3] != INDEX_SHARD_OUTER_UNCLAIMED;
  failures += shared->outer_running != 2U;

  failures += index_shard_staged_complete_claim(
      shared, &claim.staged,
      INDEX_SHARD_STAGED_EXECUTE_OK, 0.0, 0ULL) != 0;
  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &task,
      INDEX_SHARD_STAGED_TASK_STOPPED) != 0;
  shared->queue_waiters = 0U;
  pthread_mutex_unlock(&shared->queue_mutex);
  contexts[0].published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_staged_submit_rearm_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_staged_group_t groups[2];
  index_shard_staged_task_t tasks[2][2];
  index_shard_thread_state_t *shared;
  int owner = -1;
  int failures = 0;
  int i;
  int j;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(groups, 0, sizeof(groups));
  memset(tasks, 0, sizeof(tasks));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 2;
  pool.generation = 8U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  shared->observability_enabled = TRUE;

  for (i = 0; i < 2; i++) {
    contexts[i].worker_id = i;
    contexts[i].pool = &pool;
    contexts[i].generation_seen = pool.generation;
    contexts[i].current_outer_active = TRUE;
    contexts[i].current_index_order = i ? 3U : 5U;
    contexts[i].staged_group_epoch = 21U + (unsigned long long)i;
    contexts[i].published_staged_group = &groups[i];
    groups[i].pool = &pool;
    groups[i].tasks = tasks[i];
    groups[i].task_count = 2U;
    groups[i].generation = pool.generation;
    groups[i].owner_epoch = contexts[i].staged_group_epoch;
    groups[i].owner_worker = i;
    groups[i].owner_index_order = contexts[i].current_index_order;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      failures += index_shard_staged_set_state_locked(
          shared,
          &groups[i],
          &tasks[i][j],
          INDEX_SHARD_STAGED_TASK_SUBMIT_READY) != 0;
      failures += index_shard_staged_set_submit_wait_locked(
          &groups[i], &tasks[i][j], TRUE) != 0;
    }
  }

  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &owner) != 1;
  failures += owner != 1;
  failures += groups[1].submit_ready_mask != UINT64_C(1);
  failures += groups[1].submit_wait_mask != UINT64_C(2);
  failures += groups[1].submit_credit_mask != UINT64_C(1);
  failures += groups[0].submit_ready_mask != UINT64_C(0);
  failures += groups[0].submit_wait_mask != UINT64_C(3);

  owner = -1;
  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &owner) != 1;
  failures += owner != 1;
  failures += groups[1].submit_ready_mask != UINT64_C(3);
  failures += groups[1].submit_wait_mask != UINT64_C(0);
  failures += groups[1].submit_credit_mask != UINT64_C(3);

  owner = -1;
  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &owner) != 1;
  failures += owner != 0;
  failures += groups[0].submit_ready_mask != UINT64_C(1);
  failures += groups[0].submit_wait_mask != UINT64_C(2);
  failures += groups[0].submit_credit_mask != UINT64_C(1);

  owner = -1;
  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &owner) != 1;
  failures += owner != 0;
  failures += groups[0].submit_ready_mask != UINT64_C(3);
  failures += groups[0].submit_wait_mask != UINT64_C(0);
  failures += groups[0].submit_credit_mask != UINT64_C(3);

  owner = 0;
  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &owner) != 0;
  failures += owner != -1;
  failures += shared->staged_submit_rearms != 4U;
  failures += !shared->staged_submit_backpressure;
  pthread_mutex_unlock(&shared->queue_mutex);

  contexts[0].published_staged_group = NULL;
  contexts[1].published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_staged_submit_handoff_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t context;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_staged_claim_t claim;
  index_shard_thread_state_t *shared;
  int failures = 0;

  memset(&pool, 0, sizeof(pool));
  memset(&context, 0, sizeof(context));
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&claim, 0, sizeof(claim));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 1;
  pool.generation = 9U;
  pool.contexts = &context;
  shared->pool = &pool;
  shared->worker_count = 1;
  shared->observability_enabled = TRUE;
  context.worker_id = 0;
  context.pool = &pool;
  context.generation_seen = pool.generation;
  context.current_outer_active = TRUE;
  context.current_index_order = 4U;
  context.staged_group_epoch = 31U;
  context.published_staged_group = &group;
  group.pool = &pool;
  group.tasks = tasks;
  group.task_count = 2U;
  group.generation = pool.generation;
  group.owner_epoch = context.staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = context.current_index_order;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_SUBMITTING) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_SUBMIT_READY) != 0;
  failures += index_shard_staged_set_submit_wait_locked(
      &group, &tasks[1], TRUE) != 0;
  group.running_count = 1U;
  shared->staged_submit_callbacks_active = 1U;
  pthread_mutex_unlock(&shared->queue_mutex);

  claim.group = &group;
  claim.task_index = 0U;
  claim.kind = INDEX_SHARD_STAGED_CLAIM_SUBMIT;
  claim.owner_claim = TRUE;
  claim.submit_credit = TRUE;
  failures += index_shard_staged_complete_claim(
      shared,
      &claim,
      INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY,
      0.0,
      0ULL) != 0;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY;
  failures += group.compute_ready_mask != UINT64_C(1);
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_SUBMIT_READY;
  failures += group.submit_ready_mask != UINT64_C(2);
  failures += group.submit_wait_mask != UINT64_C(0);
  failures += group.submit_credit_mask != UINT64_C(2);
  failures += group.running_count != 0U;
  failures += shared->staged_submit_callbacks_active != 0U;
  failures += shared->staged_submit_rearms != 1U;
  failures += shared->staged_submit_handoffs != 1U;
  pthread_mutex_unlock(&shared->queue_mutex);

  context.published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_staged_submit_backpressure_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t context;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[3];
  index_shard_staged_claim_t claim;
  index_shard_staged_claim_t next_claim;
  index_shard_thread_state_t *shared;
  int rearmed_owner = -1;
  int failures = 0;
  size_t i;

  memset(&pool, 0, sizeof(pool));
  memset(&context, 0, sizeof(context));
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&claim, 0, sizeof(claim));
  memset(&next_claim, 0, sizeof(next_claim));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 1;
  pool.generation = 10U;
  pool.contexts = &context;
  shared->pool = &pool;
  shared->worker_count = 1;
  shared->observability_enabled = TRUE;
  context.worker_id = 0;
  context.pool = &pool;
  context.generation_seen = pool.generation;
  context.current_outer_active = TRUE;
  context.current_index_order = 6U;
  context.staged_group_epoch = 41U;
  context.published_staged_group = &group;
  group.pool = &pool;
  group.tasks = tasks;
  group.task_count = 3U;
  group.generation = pool.generation;
  group.owner_epoch = context.staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = context.current_index_order;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_SUBMITTING) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_SUBMIT_READY) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[2],
      INDEX_SHARD_STAGED_TASK_SUBMIT_READY) != 0;
  group.running_count = 1U;
  shared->staged_submit_callbacks_active = 1U;
  pthread_mutex_unlock(&shared->queue_mutex);

  claim.group = &group;
  claim.task_index = 0U;
  claim.kind = INDEX_SHARD_STAGED_CLAIM_SUBMIT;
  claim.owner_claim = TRUE;
  failures += index_shard_staged_complete_claim(
      shared,
      &claim,
      INDEX_SHARD_STAGED_SUBMIT_RETRY,
      0.0,
      0ULL) != 0;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += !shared->staged_submit_backpressure;
  failures += group.submit_ready_mask != UINT64_C(0);
  failures += group.submit_wait_mask != UINT64_C(7);
  failures += group.submit_credit_mask != UINT64_C(0);
  failures += shared->staged_submit_retries != 1U;
  failures += shared->staged_submit_deferrals != 2U;
  failures += index_shard_staged_rearm_one_submit_waiter_locked(
      &pool, &rearmed_owner) != 1;
  failures += rearmed_owner != 0;
  failures += group.submit_ready_mask != UINT64_C(1);
  failures += group.submit_wait_mask != UINT64_C(6);
  failures += group.submit_credit_mask != UINT64_C(1);
  failures += index_shard_staged_select_locked(
      &context, shared,
      INDEX_SHARD_STAGED_SELECT_SUBMIT,
      TRUE, FALSE, &next_claim) != 0;
  failures += next_claim.task_index != 0U;
  failures += !next_claim.submit_credit;
  failures += group.submit_ready_mask != UINT64_C(0);
  failures += group.submit_credit_mask != UINT64_C(0);
  failures += index_shard_staged_select_locked(
      &context, shared,
      INDEX_SHARD_STAGED_SELECT_SUBMIT,
      TRUE, FALSE, &claim) != 1;

  group.running_count = 0U;
  shared->staged_submit_callbacks_active = 0U;
  for (i = 0U; i < 3U; i++) {
    failures += index_shard_staged_set_state_locked(
        shared, &group, &tasks[i],
        INDEX_SHARD_STAGED_TASK_STOPPED) != 0;
  }
  index_shard_staged_refresh_submit_backpressure_locked(&pool);
  failures += shared->staged_submit_backpressure;
  pthread_mutex_unlock(&shared->queue_mutex);

  context.published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_completion_registry_test(void) {
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_thread_state_t *shared;
  anbool early = FALSE;
  int failures = 0;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 2;
  pool.generation = 9U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  shared->observability_enabled = TRUE;
  contexts[0].worker_id = 0;
  contexts[0].pool = &pool;
  contexts[0].generation_seen = pool.generation;
  contexts[0].current_outer_active = TRUE;
  contexts[0].current_index_order = 4U;
  contexts[0].staged_group_epoch = 13U;
  contexts[0].published_staged_group = &group;
  group.pool = &pool;
  group.tasks = tasks;
  group.task_count = 2U;
  group.generation = pool.generation;
  group.owner_epoch = contexts[0].staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = contexts[0].current_index_order;

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[0],
      INDEX_SHARD_STAGED_TASK_IO_SUBMITTED) != 0;
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_SUBMIT_READY) != 0;
  failures += index_shard_staged_set_submit_wait_locked(
      &group, &tasks[1], TRUE) != 0;
  tasks[0].completion_id = 101U;
  failures += index_shard_completion_registry_register_locked(
      shared, &group, 0U, 101U, &early) != 0;
  failures += early;
  pthread_mutex_unlock(&shared->queue_mutex);

  index_shard_staged_completion_notify(&pool, 101U);

  pthread_mutex_lock(&shared->queue_mutex);
  failures += !tasks[0].completion_pending;
  failures += !(group.completion_pending_mask & UINT64_C(1));
  failures += group.submit_wait_mask != UINT64_C(0);
  failures += group.submit_ready_mask != UINT64_C(2);
  failures += group.submit_credit_mask != UINT64_C(2);
  failures += shared->staged_submit_rearms != 1U;
  failures += shared->completion_groups_scanned != 0U;
  failures += shared->completion_tasks_scanned != 0U;
  failures += shared->completion_matches != 1U;
  failures += index_shard_completion_registry_remove_locked(
      shared, &group, 0U, 101U) != 0;
  tasks[0].completion_id = 0U;
  failures += index_shard_staged_set_completion_pending_locked(
      &group, &tasks[0], FALSE) != 0;

  shared->staged_submit_callbacks_active = 1U;
  pthread_mutex_unlock(&shared->queue_mutex);
  index_shard_staged_completion_notify(&pool, 202U);
  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared, &group, &tasks[1],
      INDEX_SHARD_STAGED_TASK_IO_SUBMITTED) != 0;
  tasks[1].completion_id = 202U;
  early = FALSE;
  failures += index_shard_completion_registry_register_locked(
      shared, &group, 1U, 202U, &early) != 0;
  failures += !early;
  shared->staged_submit_callbacks_active = 0U;
  failures += index_shard_staged_set_completion_pending_locked(
      &group, &tasks[1], early) != 0;
  failures += index_shard_completion_registry_remove_locked(
      shared, &group, 1U, 202U) != 0;
  tasks[1].completion_id = 0U;
  failures += index_shard_staged_set_completion_pending_locked(
      &group, &tasks[1], FALSE) != 0;
  failures += shared->completion_active != 0U;
  failures += group.completion_registry_entries != 0U;
  failures += shared->completion_registry_early != 1U;
  failures += shared->completion_registry_invalid != 0U;
  pthread_mutex_unlock(&shared->queue_mutex);

  contexts[0].published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

index_shard_staged_io_status_t
index_shard_inline_poll_test_callback(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes) {
  (void)input;
  (void)input_bytes;
  if (!output || output_bytes !=
          sizeof(index_shard_staged_io_status_t)) {
    return INDEX_SHARD_STAGED_IO_ERROR;
  }
  return *(index_shard_staged_io_status_t *)output;
}

int index_shard_completion_inline_poll_case(
    index_shard_staged_io_status_t poll_status,
    index_shard_staged_task_state_t expected_state,
    anbool expect_stop,
    anbool expect_error) {
  static const index_shard_staged_ops_t ops = {
      "completion-inline-poll-test",
      NULL,
      NULL,
      index_shard_inline_poll_test_callback,
      NULL,
      NULL,
      NULL,
      TRUE};
  index_shard_pool_t pool;
  index_shard_worker_context_t contexts[2];
  index_shard_staged_group_t group;
  index_shard_staged_task_t task;
  index_shard_thread_state_t *shared;
  anbool early = FALSE;
  int failures = 0;

  memset(&pool, 0, sizeof(pool));
  memset(contexts, 0, sizeof(contexts));
  memset(&group, 0, sizeof(group));
  memset(&task, 0, sizeof(task));
  shared = &pool.shared;
  if (index_shard_staged_retire_test_init(shared)) {
    return 1;
  }
  pool.worker_count = 2;
  pool.generation = 11U;
  pool.contexts = contexts;
  shared->pool = &pool;
  shared->worker_count = 2;
  shared->observability_enabled = TRUE;
  shared->staged_tickets_active = 1U;
  shared->staged_source_leases = 1U;
  shared->staged_io_submitted = 1U;
  contexts[0].worker_id = 0;
  contexts[0].pool = &pool;
  contexts[0].generation_seen = pool.generation;
  contexts[0].current_outer_active = TRUE;
  contexts[0].current_index_order = 6U;
  contexts[0].staged_group_epoch = 17U;
  contexts[0].published_staged_group = &group;
  group.pool = &pool;
  group.ops = &ops;
  group.tasks = &task;
  group.task_count = 1U;
  group.generation = pool.generation;
  group.owner_epoch = contexts[0].staged_group_epoch;
  group.owner_worker = 0;
  group.owner_index_order = contexts[0].current_index_order;
  group.io_submitted = 1U;
  task.output = &poll_status;
  task.output_bytes = sizeof(poll_status);

  pthread_mutex_lock(&shared->queue_mutex);
  failures += index_shard_staged_set_state_locked(
      shared,
      &group,
      &task,
      INDEX_SHARD_STAGED_TASK_IO_SUBMITTED) != 0;
  task.completion_id = 303U;
  failures += index_shard_completion_registry_register_locked(
      shared, &group, 0U, task.completion_id, &early) != 0;
  failures += early;
  pthread_mutex_unlock(&shared->queue_mutex);

  index_shard_staged_completion_notify(&pool, 303U);

  pthread_mutex_lock(&shared->queue_mutex);
  failures += task.scheduler_state != expected_state;
  failures += task.completion_id != 0U;
  failures += task.completion_pending;
  failures += group.completion_pending_mask != UINT64_C(0);
  failures += group.io_submitted != 0U;
  failures += group.io_completed != 1U;
  failures += group.running_count != 0U;
  failures += group.poll_claims != 1U;
  failures += group.inline_poll_claims != 1U;
  failures += group.completion_registry_entries != 0U;
  failures += shared->staged_tickets_active != 0U;
  failures += shared->staged_source_leases != 0U;
  failures += shared->staged_io_completed != 1U;
  failures += group.stop_seen != expect_stop;
  failures += group.task_failed != expect_error;
  failures += group.internal_error != expect_error;
  failures += shared->completion_registry_error != expect_error;
  failures += shared->staged_inline_poll_failures !=
      (expect_error ? 1U : 0U);
  failures += shared->completion_matches != 1U;
  failures += shared->completion_active != 0U;
  failures += shared->completion_registry_invalid != 0U;
  pthread_mutex_unlock(&shared->queue_mutex);

  contexts[0].published_staged_group = NULL;
  index_shard_staged_retire_test_destroy(shared);
  return failures;
}

int index_shard_completion_inline_poll_test(void) {
  return index_shard_completion_inline_poll_case(
      INDEX_SHARD_STAGED_IO_READY,
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY,
      FALSE,
      FALSE) +
      index_shard_completion_inline_poll_case(
          INDEX_SHARD_STAGED_IO_FAILED,
          INDEX_SHARD_STAGED_TASK_OWNER_READY,
          FALSE,
          FALSE) +
      index_shard_completion_inline_poll_case(
          INDEX_SHARD_STAGED_IO_CANCELLED,
          INDEX_SHARD_STAGED_TASK_STOPPED,
          TRUE,
          FALSE) +
      index_shard_completion_inline_poll_case(
          INDEX_SHARD_STAGED_IO_ERROR,
          INDEX_SHARD_STAGED_TASK_FAILED,
          FALSE,
          TRUE);
}

int index_shard_test_staged_retire_more(void) {
  return index_shard_staged_retire_test_order() +
      index_shard_staged_retire_test_terminal() +
      index_shard_staged_child_helper_test() +
      index_shard_observability_counter_test() +
      index_shard_staged_mask_selection_test() +
      index_shard_ready_before_outer_test() +
      index_shard_outer_cap_liveness_test() +
      index_shard_staged_submit_rearm_test() +
      index_shard_staged_submit_handoff_test() +
      index_shard_staged_submit_backpressure_test() +
      index_shard_completion_registry_test() +
      index_shard_completion_inline_poll_test();
}
#endif
