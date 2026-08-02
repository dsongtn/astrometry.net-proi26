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
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/*
 * Inner claim types, completion registry, and staged state selection.
 *
 * This module owns staged completion routing and READY task transitions.
 */


unsigned long long index_shard_helper_task_work(
    const index_shard_helper_task_t *task) {
  return task->work_units ? task->work_units : 1U;
}

int index_shard_helper_add_work(
    unsigned long long *total,
    unsigned long long work) {
  if (ULLONG_MAX - *total < work) {
    return -1;
  }
  *total += work;
  return 0;
}

static size_t index_shard_completion_hash(
    unsigned long long completion_id,
    size_t bucket_count) {
  uint64_t value = (uint64_t)completion_id;

  value ^= value >> 33;
  value *= UINT64_C(0xff51afd7ed558ccd);
  value ^= value >> 33;
  value *= UINT64_C(0xc4ceb9fe1a85ec53);
  value ^= value >> 33;
  return (size_t)value & (bucket_count - 1U);
}

int index_shard_completion_registry_init(
    index_shard_thread_state_t *shared,
    int worker_count) {
  size_t entry_capacity;
  size_t bucket_target;
  size_t bucket_count = 1U;
  size_t i;

  if (!shared || worker_count < 1 ||
      (size_t)worker_count >
          SIZE_MAX / INDEX_SHARD_HELPER_MAX_TASKS) {
    return -1;
  }
  entry_capacity =
      (size_t)worker_count * INDEX_SHARD_HELPER_MAX_TASKS;
  if (!entry_capacity || entry_capacity > SIZE_MAX / 2U) {
    return -1;
  }
  bucket_target = entry_capacity * 2U;
  while (bucket_count < bucket_target) {
    if (bucket_count > SIZE_MAX / 2U) {
      return -1;
    }
    bucket_count *= 2U;
  }

  shared->completion_buckets = calloc(
      bucket_count, sizeof(shared->completion_buckets[0]));
  shared->completion_entries = calloc(
      entry_capacity, sizeof(shared->completion_entries[0]));
  if (!shared->completion_buckets || !shared->completion_entries) {
    free(shared->completion_buckets);
    free(shared->completion_entries);
    shared->completion_buckets = NULL;
    shared->completion_entries = NULL;
    return -1;
  }
  for (i = 0U; i < bucket_count; i++) {
    shared->completion_buckets[i] =
        INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  for (i = 0U; i < entry_capacity; i++) {
    shared->completion_entries[i].next =
        i + 1U < entry_capacity
            ? i + 1U
            : INDEX_SHARD_COMPLETION_SLOT_NONE;
    shared->completion_entries[i].state =
        INDEX_SHARD_COMPLETION_ENTRY_FREE;
  }
  shared->completion_bucket_count = bucket_count;
  shared->completion_entry_capacity = entry_capacity;
  shared->completion_free_head = 0U;
  shared->completion_active = 0U;
  shared->staged_submit_callbacks_active = 0U;
  shared->completion_registry_error = FALSE;
  return 0;
}

void index_shard_completion_registry_destroy(
    index_shard_thread_state_t *shared) {
  if (!shared) {
    return;
  }
  free(shared->completion_buckets);
  free(shared->completion_entries);
  shared->completion_buckets = NULL;
  shared->completion_entries = NULL;
  shared->completion_bucket_count = 0U;
  shared->completion_entry_capacity = 0U;
  shared->completion_free_head = INDEX_SHARD_COMPLETION_SLOT_NONE;
  shared->completion_active = 0U;
  shared->staged_submit_callbacks_active = 0U;
  shared->completion_registry_error = FALSE;
}

static size_t index_shard_completion_registry_find_locked(
    const index_shard_thread_state_t *shared,
    unsigned long long completion_id,
    size_t *previous_out,
    size_t *bucket_out) {
  size_t bucket;
  size_t previous = INDEX_SHARD_COMPLETION_SLOT_NONE;
  size_t current;
  size_t visited = 0U;

  if (previous_out) {
    *previous_out = INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  if (bucket_out) {
    *bucket_out = INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  if (!shared || !completion_id ||
      !shared->completion_buckets ||
      !shared->completion_entries ||
      !shared->completion_bucket_count ||
      (shared->completion_bucket_count &
       (shared->completion_bucket_count - 1U))) {
    return INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  bucket = index_shard_completion_hash(
      completion_id, shared->completion_bucket_count);
  current = shared->completion_buckets[bucket];
  while (current != INDEX_SHARD_COMPLETION_SLOT_NONE) {
    const index_shard_completion_entry_t *entry;

    if (current >= shared->completion_entry_capacity ||
        visited++ >= shared->completion_entry_capacity) {
      return INDEX_SHARD_COMPLETION_SLOT_NONE;
    }
    entry = &shared->completion_entries[current];
    if (entry->state != INDEX_SHARD_COMPLETION_ENTRY_FREE &&
        entry->completion_id == completion_id) {
      if (previous_out) {
        *previous_out = previous;
      }
      if (bucket_out) {
        *bucket_out = bucket;
      }
      return current;
    }
    previous = current;
    current = entry->next;
  }
  if (bucket_out) {
    *bucket_out = bucket;
  }
  return INDEX_SHARD_COMPLETION_SLOT_NONE;
}

static size_t index_shard_completion_registry_allocate_locked(
    index_shard_thread_state_t *shared,
    unsigned long long completion_id) {
  index_shard_completion_entry_t *entry;
  size_t bucket;
  size_t slot;

  if (!shared || !completion_id ||
      !shared->completion_buckets ||
      !shared->completion_entries) {
    return INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  slot = shared->completion_free_head;
  if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE ||
      slot >= shared->completion_entry_capacity ||
      shared->completion_active >=
          shared->completion_entry_capacity) {
    return INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  entry = &shared->completion_entries[slot];
  if (entry->state != INDEX_SHARD_COMPLETION_ENTRY_FREE) {
    return INDEX_SHARD_COMPLETION_SLOT_NONE;
  }
  shared->completion_free_head = entry->next;
  bucket = index_shard_completion_hash(
      completion_id, shared->completion_bucket_count);
  memset(entry, 0, sizeof(*entry));
  entry->completion_id = completion_id;
  entry->next = shared->completion_buckets[bucket];
  shared->completion_buckets[bucket] = slot;
  shared->completion_active++;
  return slot;
}

/* queue_mutex must be held. */
static int index_shard_completion_registry_unlink_locked(
    index_shard_thread_state_t *shared,
    size_t slot,
    size_t previous,
    size_t bucket) {
  index_shard_completion_entry_t *entry;

  if (!shared || slot >= shared->completion_entry_capacity ||
      bucket >= shared->completion_bucket_count ||
      !shared->completion_active) {
    return -1;
  }
  entry = &shared->completion_entries[slot];
  if (entry->state == INDEX_SHARD_COMPLETION_ENTRY_FREE) {
    return -1;
  }
  if (previous == INDEX_SHARD_COMPLETION_SLOT_NONE) {
    if (shared->completion_buckets[bucket] != slot) {
      return -1;
    }
    shared->completion_buckets[bucket] = entry->next;
  } else {
    if (previous >= shared->completion_entry_capacity ||
        shared->completion_entries[previous].next != slot) {
      return -1;
    }
    shared->completion_entries[previous].next = entry->next;
  }
  memset(entry, 0, sizeof(*entry));
  entry->state = INDEX_SHARD_COMPLETION_ENTRY_FREE;
  entry->next = shared->completion_free_head;
  shared->completion_free_head = slot;
  shared->completion_active--;
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_completion_registry_record_early_locked(
    index_shard_thread_state_t *shared,
    unsigned long long completion_id) {
  index_shard_completion_entry_t *entry;
  size_t slot;

  if (!shared || !completion_id ||
      index_shard_completion_registry_find_locked(
          shared, completion_id, NULL, NULL) !=
          INDEX_SHARD_COMPLETION_SLOT_NONE) {
    return -1;
  }
  slot = index_shard_completion_registry_allocate_locked(
      shared, completion_id);
  if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE) {
    return -1;
  }
  entry = &shared->completion_entries[slot];
  entry->state = INDEX_SHARD_COMPLETION_ENTRY_EARLY;
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->completion_registry_early);
  }
  return 0;
}

/* queue_mutex must be held. */
int index_shard_completion_registry_register_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    size_t task_index,
    unsigned long long completion_id,
    anbool *already_notified) {
  index_shard_completion_entry_t *entry;
  size_t slot;

  if (already_notified) {
    *already_notified = FALSE;
  }
  if (!shared || !group || !group->pool ||
      !completion_id || task_index >= group->task_count ||
      group->owner_worker < 0 ||
      group->owner_worker >= shared->worker_count ||
      group->generation != group->pool->generation) {
    return -1;
  }
  slot = index_shard_completion_registry_find_locked(
      shared, completion_id, NULL, NULL);
  if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE) {
    slot = index_shard_completion_registry_allocate_locked(
        shared, completion_id);
    if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE) {
      return -1;
    }
    entry = &shared->completion_entries[slot];
    entry->state = INDEX_SHARD_COMPLETION_ENTRY_REGISTERED;
  } else {
    entry = &shared->completion_entries[slot];
    if (entry->state != INDEX_SHARD_COMPLETION_ENTRY_EARLY) {
      if (shared->observability_enabled) {
        index_shard_observability_increment(
            &shared->completion_registry_duplicates);
      }
      return -1;
    }
    entry->state = INDEX_SHARD_COMPLETION_ENTRY_NOTIFIED;
    if (already_notified) {
      *already_notified = TRUE;
    }
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->completion_matches);
    }
  }
  entry->pool_generation = group->generation;
  entry->owner_epoch = group->owner_epoch;
  entry->owner_worker = (size_t)group->owner_worker;
  entry->owner_index_order = group->owner_index_order;
  entry->task_index = task_index;
  group->completion_registry_entries++;
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->completion_registry_registered);
  }
  return 0;
}

/* queue_mutex must be held. */
int index_shard_completion_registry_remove_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    size_t task_index,
    unsigned long long completion_id) {
  index_shard_completion_entry_t *entry;
  size_t previous;
  size_t bucket;
  size_t slot;

  if (!shared || !group || !completion_id) {
    return -1;
  }
  slot = index_shard_completion_registry_find_locked(
      shared, completion_id, &previous, &bucket);
  if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE ||
      slot >= shared->completion_entry_capacity ||
      bucket >= shared->completion_bucket_count) {
    return -1;
  }
  entry = &shared->completion_entries[slot];
  if ((entry->state != INDEX_SHARD_COMPLETION_ENTRY_REGISTERED &&
       entry->state != INDEX_SHARD_COMPLETION_ENTRY_NOTIFIED) ||
      entry->pool_generation != group->generation ||
      entry->owner_epoch != group->owner_epoch ||
      entry->owner_worker != (size_t)group->owner_worker ||
      entry->owner_index_order != group->owner_index_order ||
      entry->task_index != task_index ||
      !group->completion_registry_entries) {
    return -1;
  }
  if (index_shard_completion_registry_unlink_locked(
          shared, slot, previous, bucket)) {
    return -1;
  }
  group->completion_registry_entries--;
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->completion_registry_removed);
  }
  return 0;
}

static uint64_t index_shard_staged_task_mask(size_t task_count) {
  if (!task_count || task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    return UINT64_C(0);
  }
  if (task_count == INDEX_SHARD_HELPER_MAX_TASKS) {
    return UINT64_MAX;
  }
  return (UINT64_C(1) << task_count) - UINT64_C(1);
}

static size_t index_shard_staged_lowest_task(uint64_t mask) {
  size_t task_index = 0U;

  if (!mask) {
    return SIZE_MAX;
  }
  while (!(mask & UINT64_C(1))) {
    mask >>= 1;
    task_index++;
  }
  return task_index;
}

/* queue_mutex must be held. */
static int index_shard_staged_task_bit_locked(
    const index_shard_staged_group_t *group,
    const index_shard_staged_task_t *task,
    uint64_t *bit_out,
    size_t *task_index_out) {
  uintptr_t base;
  uintptr_t address;
  size_t byte_offset;
  size_t task_index;

  if (bit_out) {
    *bit_out = UINT64_C(0);
  }
  if (task_index_out) {
    *task_index_out = SIZE_MAX;
  }
  if (!group || !group->tasks || !task ||
      !group->task_count ||
      group->task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    return -1;
  }
  base = (uintptr_t)group->tasks;
  address = (uintptr_t)task;
  if (address < base) {
    return -1;
  }
  byte_offset = (size_t)(address - base);
  if (byte_offset % sizeof(group->tasks[0])) {
    return -1;
  }
  task_index = byte_offset / sizeof(group->tasks[0]);
  if (task_index >= group->task_count) {
    return -1;
  }
  if (bit_out) {
    *bit_out = UINT64_C(1) << task_index;
  }
  if (task_index_out) {
    *task_index_out = task_index;
  }
  return 0;
}

/* queue_mutex must be held. */
int index_shard_staged_set_submit_wait_locked(
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    anbool waiting) {
  uint64_t bit;

  if (!group || !task ||
      task->scheduler_state !=
          INDEX_SHARD_STAGED_TASK_SUBMIT_READY ||
      index_shard_staged_task_bit_locked(
          group, task, &bit, NULL)) {
    return -1;
  }
  if (waiting) {
    group->submit_ready_mask &= ~bit;
    group->submit_wait_mask |= bit;
    group->submit_credit_mask &= ~bit;
  } else {
    group->submit_wait_mask &= ~bit;
    group->submit_ready_mask |= bit;
  }
  return 0;
}

/* queue_mutex must be held. */
int index_shard_staged_set_completion_pending_locked(
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    anbool pending) {
  uint64_t bit;

  if (!group || !task ||
      index_shard_staged_task_bit_locked(
          group, task, &bit, NULL)) {
    return -1;
  }
  task->completion_pending = pending;
  if (pending) {
    group->completion_pending_mask |= bit;
  } else {
    group->completion_pending_mask &= ~bit;
  }
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_staged_set_cancel_sent_locked(
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    anbool sent) {
  uint64_t bit;

  if (!group || !task ||
      index_shard_staged_task_bit_locked(
          group, task, &bit, NULL)) {
    return -1;
  }
  task->cancel_sent = sent;
  if (sent) {
    group->cancel_sent_mask |= bit;
  } else {
    group->cancel_sent_mask &= ~bit;
  }
  return 0;
}

static int index_shard_staged_group_valid_locked(
    const index_shard_pool_t *pool,
    const index_shard_worker_context_t *owner,
    const index_shard_staged_group_t *group);

/* queue_mutex must be held. */
void index_shard_staged_refresh_submit_backpressure_locked(
    index_shard_pool_t *pool) {
  index_shard_thread_state_t *shared;
  int owner;

  if (!pool || !pool->contexts) {
    return;
  }
  shared = &pool->shared;
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context =
        &pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;

    if (group &&
        (group->submit_wait_mask ||
         group->submit_credit_mask)) {
      shared->staged_submit_backpressure = TRUE;
      return;
    }
  }
  shared->staged_submit_backpressure = FALSE;
}

/*
 * Once bounded provider admission refuses one packet, park every other
 * uncredited submission before it enters the allocation and range-building
 * path. Completion credits remain runnable and cannot be stolen by newly
 * prepared packets.
 *
 * queue_mutex must be held.
 */
static int index_shard_staged_park_uncredited_submitters_locked(
    index_shard_pool_t *pool) {
  index_shard_thread_state_t *shared;
  int owner;

  if (!pool || !pool->contexts) {
    return -1;
  }
  shared = &pool->shared;
  shared->staged_submit_backpressure = TRUE;
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context =
        &pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;
    uint64_t mask;

    if (!group) {
      continue;
    }
    if (index_shard_staged_group_valid_locked(
            pool, context, group)) {
      group->internal_error = TRUE;
      return -1;
    }
    mask = group->submit_ready_mask &
        ~group->submit_credit_mask;
    while (mask) {
      size_t task_index =
          index_shard_staged_lowest_task(mask);
      uint64_t bit;

      if (task_index >= group->task_count) {
        group->internal_error = TRUE;
        return -1;
      }
      bit = UINT64_C(1) << task_index;
      if (index_shard_staged_set_submit_wait_locked(
              group, &group->tasks[task_index], TRUE)) {
        group->internal_error = TRUE;
        return -1;
      }
      mask &= ~bit;
      if (shared->observability_enabled) {
        index_shard_observability_increment(
            &shared->staged_submit_deferrals);
      }
    }
  }
  return 0;
}

/*
 * Rearm one canonical waiter for one released provider admission. The
 * provider frees one job before delivering each completion callback. Waking
 * every waiter here turns that single credit into a retry herd at wider
 * worker counts.
 *
 * queue_mutex must be held.
 */
int index_shard_staged_rearm_one_submit_waiter_locked(
    index_shard_pool_t *pool,
    int *owner_out) {
  index_shard_thread_state_t *shared;
  index_shard_staged_group_t *best_group = NULL;
  size_t best_task = SIZE_MAX;
  int best_owner = -1;
  int owner;

  if (owner_out) {
    *owner_out = -1;
  }
  if (!pool || !pool->contexts || !owner_out) {
    return -1;
  }
  shared = &pool->shared;
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context = &pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;
    size_t task_index;

    if (!group) {
      continue;
    }
    if (index_shard_staged_group_valid_locked(
            pool, context, group)) {
      return -1;
    }
    if (group->cancelling || group->task_failed ||
        group->stop_seen || group->internal_error) {
      continue;
    }
    task_index = index_shard_staged_lowest_task(
        group->submit_wait_mask);
    if (task_index == SIZE_MAX) {
      continue;
    }
    if (task_index >= group->task_count) {
      group->internal_error = TRUE;
      return -1;
    }
    if (!best_group ||
        group->owner_index_order < best_group->owner_index_order ||
        (group->owner_index_order == best_group->owner_index_order &&
         task_index < best_task)) {
      best_group = group;
      best_task = task_index;
      best_owner = owner;
    }
  }
  if (!best_group) {
    index_shard_staged_refresh_submit_backpressure_locked(pool);
    return 0;
  }
  if (index_shard_staged_set_submit_wait_locked(
          best_group, &best_group->tasks[best_task], FALSE)) {
    best_group->internal_error = TRUE;
    return -1;
  }
  best_group->submit_credit_mask |= UINT64_C(1) << best_task;
  shared->staged_submit_backpressure = TRUE;
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->staged_submit_rearms);
  }
  *owner_out = best_owner;
  return 1;
}

static unsigned long long index_shard_staged_task_work(
    const index_shard_staged_task_t *task) {
  return task->work_units ? task->work_units : 1U;
}

static anbool index_shard_staged_task_terminal(
    const index_shard_staged_task_t *task) {
  if (!task) {
    return FALSE;
  }
  return task->scheduler_state == INDEX_SHARD_STAGED_TASK_RETIRED ||
      task->scheduler_state == INDEX_SHARD_STAGED_TASK_STOPPED ||
      task->scheduler_state == INDEX_SHARD_STAGED_TASK_FAILED;
}

/* queue_mutex must be held. */
int index_shard_staged_set_state_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    index_shard_staged_task_state_t state) {
  index_shard_staged_task_state_t previous;
  uint64_t bit;
  int invalid = FALSE;

  if (!shared || !group || !task ||
      index_shard_staged_task_bit_locked(
          group, task, &bit, NULL)) {
    return -1;
  }
  previous = (index_shard_staged_task_state_t)
      task->scheduler_state;
  if (previous == state) {
    return 0;
  }

  switch (previous) {
  case INDEX_SHARD_STAGED_TASK_PREPARE_READY:
    group->prepare_ready_mask &= ~bit;
    break;
  case INDEX_SHARD_STAGED_TASK_SUBMIT_READY:
    group->submit_ready_mask &= ~bit;
    group->submit_wait_mask &= ~bit;
    group->submit_credit_mask &= ~bit;
    break;
  case INDEX_SHARD_STAGED_TASK_IO_SUBMITTED:
    group->io_submitted_mask &= ~bit;
    break;
  case INDEX_SHARD_STAGED_TASK_COMPUTE_READY:
    group->compute_ready_mask &= ~bit;
    if (!group->compute_ready || !shared->staged_compute_ready) {
      invalid = TRUE;
    } else {
      group->compute_ready--;
      shared->staged_compute_ready--;
    }
    break;
  case INDEX_SHARD_STAGED_TASK_OWNER_READY:
    group->owner_ready_mask &= ~bit;
    break;
  case INDEX_SHARD_STAGED_TASK_RESULTS_READY:
    group->results_ready_mask &= ~bit;
    if (!group->reorder_ready || !shared->staged_reorder_ready) {
      invalid = TRUE;
    } else {
      group->reorder_ready--;
      shared->staged_reorder_ready--;
    }
    break;
  default:
    break;
  }

  task->scheduler_state = (unsigned char)state;
  switch (state) {
  case INDEX_SHARD_STAGED_TASK_PREPARE_READY:
    group->prepare_ready_mask |= bit;
    break;
  case INDEX_SHARD_STAGED_TASK_SUBMIT_READY:
    group->submit_ready_mask |= bit;
    break;
  case INDEX_SHARD_STAGED_TASK_IO_SUBMITTED:
    group->io_submitted_mask |= bit;
    break;
  case INDEX_SHARD_STAGED_TASK_COMPUTE_READY:
    group->compute_ready_mask |= bit;
    if (group->compute_ready == SIZE_MAX ||
        shared->staged_compute_ready == SIZE_MAX) {
      invalid = TRUE;
    } else {
      group->compute_ready++;
      shared->staged_compute_ready++;
      group->max_compute_ready = MAX(
          group->max_compute_ready, group->compute_ready);
      shared->staged_max_compute_ready = MAX(
          shared->staged_max_compute_ready,
          shared->staged_compute_ready);
    }
    break;
  case INDEX_SHARD_STAGED_TASK_OWNER_READY:
    group->owner_ready_mask |= bit;
    break;
  case INDEX_SHARD_STAGED_TASK_RESULTS_READY:
    group->results_ready_mask |= bit;
    if (group->reorder_ready == SIZE_MAX ||
        shared->staged_reorder_ready == SIZE_MAX) {
      invalid = TRUE;
    } else {
      group->reorder_ready++;
      shared->staged_reorder_ready++;
      group->max_reorder_ready = MAX(
          group->max_reorder_ready, group->reorder_ready);
      shared->staged_max_reorder_ready = MAX(
          shared->staged_max_reorder_ready,
          shared->staged_reorder_ready);
    }
    break;
  default:
    break;
  }
  if (invalid) {
    group->internal_error = TRUE;
    return -1;
  }
  if ((group->submit_credit_mask &
       ~group->submit_ready_mask) ||
      (group->submit_ready_mask &
       group->submit_wait_mask)) {
    group->internal_error = TRUE;
    return -1;
  }
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_staged_group_valid_locked(
    const index_shard_pool_t *pool,
    const index_shard_worker_context_t *owner,
    const index_shard_staged_group_t *group) {
  if (!pool || !owner || !group ||
      group->pool != pool ||
      owner->worker_id < 0 ||
      owner->worker_id >= pool->shared.worker_count ||
      owner->published_staged_group != group ||
      group->owner_worker != owner->worker_id ||
      group->owner_epoch != owner->staged_group_epoch ||
      group->generation != owner->generation_seen ||
      !owner->current_outer_active ||
      group->owner_index_order != owner->current_index_order) {
    return -1;
  }
  return 0;
}

static int index_shard_staged_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim);

/*
 * Payload completion never enters with the fitsbin mutex held. The immutable
 * completion ID selects exactly one published ticket without retaining a
 * task or group pointer in the I/O service. The shared epoch is only a retry
 * event for submissions that previously met bounded queue pressure.
 */
void index_shard_staged_completion_notify(
    void *opaque,
    unsigned long long completion_id) {
  index_shard_pool_t *pool = opaque;
  index_shard_thread_state_t *shared;
  index_shard_staged_group_t *matched_group = NULL;
  index_shard_staged_task_t *matched_task = NULL;
  index_shard_completion_entry_t *entry = NULL;
  index_shard_staged_claim_t inline_claim;
  size_t slot;
  size_t matches = 0U;
  anbool inline_poll = FALSE;
  anbool invalid = FALSE;
  int rearmed_owner = -1;
  int rearm_status = 0;
  int owner;

  if (!pool || !completion_id) {
    return;
  }
  shared = &pool->shared;
  memset(&inline_claim, 0, sizeof(inline_claim));
  pthread_mutex_lock(&shared->queue_mutex);
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->completion_notifications);
  }
  shared->staged_completion_epoch++;
  if (!shared->staged_completion_epoch) {
    shared->staged_completion_epoch++;
  }

  slot = index_shard_completion_registry_find_locked(
      shared, completion_id, NULL, NULL);
  if (slot == INDEX_SHARD_COMPLETION_SLOT_NONE) {
    /*
     * submit() executes without queue_mutex. A fast provider may complete
     * after ticket publication but before the submitter can register its
     * task identity. Retain the numeric ID until that submit callback returns.
     */
    if (shared->staged_submit_callbacks_active) {
      if (index_shard_completion_registry_record_early_locked(
              shared, completion_id)) {
        invalid = TRUE;
      }
    } else {
      invalid = TRUE;
      if (shared->observability_enabled) {
        index_shard_observability_increment(
            &shared->completion_registry_misses);
      }
    }
  } else if (slot >= shared->completion_entry_capacity) {
    invalid = TRUE;
  } else {
    entry = &shared->completion_entries[slot];
    if (entry->state != INDEX_SHARD_COMPLETION_ENTRY_REGISTERED ||
        entry->owner_worker >= (size_t)shared->worker_count) {
      invalid = TRUE;
      if (shared->observability_enabled) {
        index_shard_observability_increment(
            &shared->completion_registry_duplicates);
      }
    } else {
      index_shard_worker_context_t *context =
          &pool->contexts[entry->owner_worker];

      matched_group = context->published_staged_group;
      if (!matched_group ||
          index_shard_staged_group_valid_locked(
              pool, context, matched_group) ||
          matched_group->generation != entry->pool_generation ||
          matched_group->owner_epoch != entry->owner_epoch ||
          matched_group->owner_index_order !=
              entry->owner_index_order ||
          entry->task_index >= matched_group->task_count) {
        invalid = TRUE;
      } else {
        matched_task = &matched_group->tasks[entry->task_index];
        if (matched_task->completion_id != completion_id ||
            (matched_task->scheduler_state !=
                 INDEX_SHARD_STAGED_TASK_IO_SUBMITTED &&
             matched_task->scheduler_state !=
                 INDEX_SHARD_STAGED_TASK_IO_POLLING &&
             matched_task->scheduler_state !=
                 INDEX_SHARD_STAGED_TASK_IO_CANCELLING) ||
            matched_task->completion_pending) {
          invalid = TRUE;
        } else {
          entry->state = INDEX_SHARD_COMPLETION_ENTRY_NOTIFIED;
          matches = 1U;
          if (matched_task->scheduler_state ==
                  INDEX_SHARD_STAGED_TASK_IO_SUBMITTED &&
              !matched_group->cancelling &&
              !matched_group->task_failed &&
              !matched_group->stop_seen &&
              !matched_group->internal_error &&
              matched_group->ops && matched_group->ops->poll &&
              matched_group->ops->terminal_poll_inline_safe) {
            if (matched_group->running_count == SIZE_MAX ||
                matched_group->poll_claims == SIZE_MAX ||
                matched_group->inline_poll_claims == SIZE_MAX ||
                index_shard_staged_set_state_locked(
                    shared,
                    matched_group,
                    matched_task,
                    INDEX_SHARD_STAGED_TASK_IO_POLLING)) {
              invalid = TRUE;
            } else {
              matched_group->running_count++;
              matched_group->max_running = MAX(
                  matched_group->max_running,
                  matched_group->running_count);
              matched_group->poll_claims++;
              matched_group->inline_poll_claims++;
              inline_claim.group = matched_group;
              inline_claim.task_index = entry->task_index;
              inline_claim.kind =
                  INDEX_SHARD_STAGED_CLAIM_IO_POLL;
              inline_claim.owner_claim = FALSE;
              inline_claim.completion_inline = TRUE;
              inline_claim.observed_completion_epoch =
                  shared->staged_completion_epoch;
              inline_poll = TRUE;
            }
          } else if (
              index_shard_staged_set_completion_pending_locked(
                  matched_group, matched_task, TRUE)) {
            invalid = TRUE;
          }
        }
      }
    }
  }

  if (!invalid) {
    rearm_status =
        index_shard_staged_rearm_one_submit_waiter_locked(
            pool, &rearmed_owner);
    if (rearm_status < 0) {
      invalid = TRUE;
    }
  }

  if (invalid) {
    shared->completion_registry_error = TRUE;
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->completion_registry_invalid);
    }
    for (owner = 0; owner < shared->worker_count; owner++) {
      index_shard_staged_group_t *group =
          pool->contexts[owner].published_staged_group;

      if (group) {
        group->internal_error = TRUE;
      }
    }
  }
  if (shared->observability_enabled) {
    index_shard_observability_add(
        &shared->completion_matches,
        (unsigned long long)matches);
  }
  if (invalid) {
    index_shard_queue_broadcast_locked(shared);
  } else {
    index_shard_queue_signal_locked(shared);
    if (matched_group && !inline_poll) {
      index_shard_owner_signal_locked(
          shared, matched_group->owner_worker);
    }
    if (rearm_status > 0) {
      index_shard_owner_signal_locked(shared, rearmed_owner);
    }
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  if (inline_poll &&
      index_shard_staged_execute_claim(shared, &inline_claim)) {
    pthread_mutex_lock(&shared->queue_mutex);
    index_shard_observability_increment(
        &shared->staged_inline_poll_failures);
    pthread_mutex_unlock(&shared->queue_mutex);
  }
}

/* queue_mutex must be held. */
static uint64_t index_shard_staged_select_mask_locked(
    const index_shard_thread_state_t *shared,
    const index_shard_staged_group_t *group,
    index_shard_staged_select_class_t select_class,
    anbool owner_allowed) {
  uint64_t mask;

  if (!shared || !group) {
    return UINT64_C(0);
  }
  if (select_class == INDEX_SHARD_STAGED_SELECT_IO) {
    mask = group->io_submitted_mask &
        group->completion_pending_mask;
    if (group->cancelling) {
      mask |= group->io_submitted_mask &
          ~group->cancel_sent_mask;
    }
    return mask;
  }
  if (group->cancelling || group->task_failed ||
      group->internal_error || group->stop_seen) {
    return UINT64_C(0);
  }
  switch (select_class) {
  case INDEX_SHARD_STAGED_SELECT_COMPUTE:
    mask = group->compute_ready_mask;
    if (owner_allowed) {
      mask |= group->owner_ready_mask;
    }
    return mask;
  case INDEX_SHARD_STAGED_SELECT_SUBMIT:
    return shared->staged_submit_backpressure
        ? group->submit_credit_mask
        : group->submit_ready_mask;
  case INDEX_SHARD_STAGED_SELECT_PREPARE:
    return group->prepare_ready_mask;
  case INDEX_SHARD_STAGED_SELECT_IO:
  default:
    return UINT64_C(0);
  }
}

/* queue_mutex must be held. */
int index_shard_staged_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    index_shard_staged_select_class_t select_class,
    anbool allow_owner,
    anbool foreign_only,
    index_shard_staged_claim_t *claim) {
  index_shard_staged_group_t *best_group = NULL;
  size_t best_task = SIZE_MAX;
  index_shard_staged_claim_kind_t best_kind =
      INDEX_SHARD_STAGED_CLAIM_NONE;
  int owner;

  if (!worker || !worker->pool || !shared || !claim) {
    return -1;
  }
  memset(claim, 0, sizeof(*claim));
  if (shared->observability_enabled) {
    index_shard_observability_increment(
        &shared->selection_scans);
  }
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context =
        &worker->pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;
    uint64_t mask;
    size_t task_index;
    index_shard_staged_claim_kind_t kind;
    anbool owner_allowed = allow_owner &&
        owner == worker->worker_id;

    if (!group) {
      continue;
    }
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->selection_groups_scanned);
    }
    if (index_shard_staged_group_valid_locked(
            worker->pool, context, group)) {
      group->internal_error = TRUE;
      return -1;
    }
    if (foreign_only && owner == worker->worker_id) {
      continue;
    }
    mask = index_shard_staged_select_mask_locked(
        shared, group, select_class, owner_allowed);
    task_index = index_shard_staged_lowest_task(mask);
    if (task_index == SIZE_MAX || task_index >= group->task_count) {
      continue;
    }
    if (select_class == INDEX_SHARD_STAGED_SELECT_IO) {
      uint64_t bit = UINT64_C(1) << task_index;

      kind = (group->completion_pending_mask & bit)
          ? INDEX_SHARD_STAGED_CLAIM_IO_POLL
          : INDEX_SHARD_STAGED_CLAIM_IO_CANCEL;
    } else if (select_class == INDEX_SHARD_STAGED_SELECT_COMPUTE) {
      uint64_t bit = UINT64_C(1) << task_index;

      kind = (group->compute_ready_mask & bit)
          ? INDEX_SHARD_STAGED_CLAIM_EXECUTE
          : INDEX_SHARD_STAGED_CLAIM_OWNER;
    } else if (select_class == INDEX_SHARD_STAGED_SELECT_SUBMIT) {
      kind = INDEX_SHARD_STAGED_CLAIM_SUBMIT;
    } else {
      kind = INDEX_SHARD_STAGED_CLAIM_PREPARE;
    }
    if (!best_group ||
        group->owner_index_order < best_group->owner_index_order ||
        (group->owner_index_order == best_group->owner_index_order &&
         task_index < best_task)) {
      best_group = group;
      best_task = task_index;
      best_kind = kind;
    }
  }
  if (!best_group) {
    if (shared->observability_enabled) {
      index_shard_observability_increment(
          &shared->selection_misses);
    }
    return 1;
  }

  claim->owner_claim =
      best_group->owner_worker == worker->worker_id;
  if (best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    unsigned long long work = index_shard_staged_task_work(
        &best_group->tasks[best_task]);
    unsigned long long current = claim->owner_claim
        ? best_group->owner_work
        : best_group->foreign_work;

    if (shared->staged_compute_running_global == SIZE_MAX ||
        ULLONG_MAX - current < work) {
      best_group->internal_error = TRUE;
      return -1;
    }
  }

  switch (best_kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    if (index_shard_staged_set_state_locked(
            shared, best_group, &best_group->tasks[best_task],
            INDEX_SHARD_STAGED_TASK_PREPARING)) {
      return -1;
    }
    best_group->prepare_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    claim->submit_credit =
        (best_group->submit_credit_mask &
         (UINT64_C(1) << best_task)) != UINT64_C(0);
    if (shared->staged_submit_callbacks_active == SIZE_MAX ||
        index_shard_staged_set_state_locked(
            shared, best_group, &best_group->tasks[best_task],
            INDEX_SHARD_STAGED_TASK_SUBMITTING)) {
      return -1;
    }
    shared->staged_submit_callbacks_active++;
    best_group->submit_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    if (index_shard_staged_set_state_locked(
            shared, best_group, &best_group->tasks[best_task],
            INDEX_SHARD_STAGED_TASK_IO_POLLING) ||
        index_shard_staged_set_completion_pending_locked(
            best_group, &best_group->tasks[best_task], FALSE)) {
      return -1;
    }
    best_group->poll_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    if (index_shard_staged_set_state_locked(
            shared, best_group, &best_group->tasks[best_task],
            INDEX_SHARD_STAGED_TASK_IO_CANCELLING)) {
      return -1;
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    if (best_group->tasks[best_task].scheduler_ready_seconds > 0.0) {
      double now = monotonic_seconds();

      if (now >= best_group->tasks[best_task].scheduler_ready_seconds) {
        best_group->ready_dwell_seconds +=
            now - best_group->tasks[best_task].scheduler_ready_seconds;
      }
      best_group->tasks[best_task].scheduler_ready_seconds = 0.0;
    }
    if (index_shard_staged_set_state_locked(
            shared, best_group, &best_group->tasks[best_task],
            best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE
                ? INDEX_SHARD_STAGED_TASK_EXECUTING
                : INDEX_SHARD_STAGED_TASK_OWNER_EXECUTING)) {
      return -1;
    }
    if (best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
      best_group->execute_claims++;
    } else {
      best_group->owner_execute_claims++;
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    best_group->internal_error = TRUE;
    return -1;
  }

  best_group->running_count++;
  best_group->max_running = MAX(
      best_group->max_running, best_group->running_count);
  claim->group = best_group;
  claim->task_index = best_task;
  claim->kind = best_kind;
  claim->observed_completion_epoch =
      shared->staged_completion_epoch;
  if (claim->owner_claim) {
    best_group->owner_claims++;
    shared->staged_tasks_owner++;
  } else {
    best_group->foreign_claims++;
    shared->staged_tasks_foreign++;
  }
  if (best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    unsigned long long work = index_shard_staged_task_work(
        &best_group->tasks[best_task]);

    best_group->compute_running++;
    shared->staged_compute_running_global++;
    best_group->max_compute_running = MAX(
        best_group->max_compute_running,
        best_group->compute_running);
    shared->staged_max_compute_running_global = MAX(
        shared->staged_max_compute_running_global,
        shared->staged_compute_running_global);
    if (claim->owner_claim) {
      best_group->owner_compute_executes++;
      shared->staged_compute_owner++;
      best_group->owner_work += work;
    } else {
      best_group->foreign_compute_executes++;
      shared->staged_compute_foreign++;
      best_group->foreign_work += work;
    }
  }
  /* Work-conserving handoff: another waiter may claim remaining work. */
  index_shard_queue_signal_locked(shared);
  return 0;
}
static void index_shard_staged_set_failed_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task) {
  if (!shared || !group || !task) {
    return;
  }
  (void)index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_FAILED);
  group->task_failed = TRUE;
  shared->staged_task_failures++;
}

/* queue_mutex must be held. */
static int index_shard_staged_release_ticket_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task) {
  int invalid = FALSE;

  if (!shared || !group || !task) {
    return -1;
  }
  if (task->completion_id &&
      index_shard_completion_registry_remove_locked(
          shared,
          group,
          (size_t)(task - group->tasks),
          task->completion_id)) {
    invalid = TRUE;
  }
  if (index_shard_staged_set_completion_pending_locked(
          group, task, FALSE) ||
      index_shard_staged_set_cancel_sent_locked(
          group, task, FALSE)) {
    invalid = TRUE;
  }
  if (group->io_submitted) {
    group->io_submitted--;
  } else {
    invalid = TRUE;
  }
  group->io_completed++;
  if (shared->staged_tickets_active) {
    shared->staged_tickets_active--;
  } else {
    invalid = TRUE;
  }
  if (shared->staged_source_leases) {
    shared->staged_source_leases--;
  } else {
    invalid = TRUE;
  }
  shared->staged_io_completed++;
  if (invalid) {
    shared->completion_registry_error = TRUE;
    group->internal_error = TRUE;
  }
  return invalid ? -1 : 0;
}

int index_shard_staged_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim,
    int callback_status,
    double callback_seconds,
    unsigned long long completion_id) {
  index_shard_staged_group_t *group;
  index_shard_staged_task_t *task;
  index_shard_staged_task_state_t expected;
  double now;
  anbool task_failed_before;
  anbool handoff_submit_credit = FALSE;
  int rearmed_owner = -1;
  int rearm_status = 0;
  int rc = 0;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;
  pthread_mutex_lock(&shared->queue_mutex);
  if (claim->task_index >= group->task_count ||
      !group->tasks || !group->running_count ||
      (claim->completion_inline &&
       claim->kind != INDEX_SHARD_STAGED_CLAIM_IO_POLL)) {
    group->internal_error = TRUE;
    if (group->running_count) {
      group->running_count--;
    }
    index_shard_queue_broadcast_locked(shared);
    pthread_mutex_unlock(&shared->queue_mutex);
    return -1;
  }
  task = &group->tasks[claim->task_index];
  task_failed_before = group->task_failed;
  now = monotonic_seconds();
  switch (claim->kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    expected = INDEX_SHARD_STAGED_TASK_PREPARING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    expected = INDEX_SHARD_STAGED_TASK_SUBMITTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    expected = INDEX_SHARD_STAGED_TASK_IO_POLLING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    expected = INDEX_SHARD_STAGED_TASK_IO_CANCELLING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
    expected = INDEX_SHARD_STAGED_TASK_EXECUTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    expected = INDEX_SHARD_STAGED_TASK_OWNER_EXECUTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    expected = INDEX_SHARD_STAGED_TASK_UNUSED;
    group->internal_error = TRUE;
    rc = -1;
    break;
  }
  group->running_count--;
  if (claim->kind == INDEX_SHARD_STAGED_CLAIM_SUBMIT) {
    if (!shared->staged_submit_callbacks_active) {
      shared->completion_registry_error = TRUE;
      group->internal_error = TRUE;
      rc = -1;
    } else {
      shared->staged_submit_callbacks_active--;
    }
  }
  if (claim->kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    if (!group->compute_running ||
        !shared->staged_compute_running_global) {
      group->internal_error = TRUE;
      rc = -1;
    } else {
      group->compute_running--;
      shared->staged_compute_running_global--;
    }
  }
  task->callback_status = callback_status;
  if (task->scheduler_state != expected) {
    group->internal_error = TRUE;
    rc = -1;
  }
  if (claim->kind != INDEX_SHARD_STAGED_CLAIM_SUBMIT && completion_id) {
    group->internal_error = TRUE;
    rc = -1;
  }

  if (claim->kind == INDEX_SHARD_STAGED_CLAIM_PREPARE) {
    switch ((index_shard_staged_prepare_status_t)callback_status) {
    case INDEX_SHARD_STAGED_PREPARE_MORE:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_PREPARE_READY);
      break;
    case INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_SUBMIT_READY);
      task->scheduler_epoch = ~shared->staged_completion_epoch;
      if (shared->staged_submit_backpressure) {
        if (index_shard_staged_set_submit_wait_locked(
                group, task, TRUE)) {
          group->internal_error = TRUE;
          rc = -1;
        } else if (shared->observability_enabled) {
          index_shard_observability_increment(
              &shared->staged_submit_deferrals);
        }
      }
      break;
    case INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
      task->scheduler_ready_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_OWNER_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_OWNER_READY);
      task->scheduler_ready_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_RESULTS_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_RESULTS_READY);
      task->scheduler_result_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_STOPPED:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_PREPARE_ERROR:
    default:
      index_shard_staged_set_failed_locked(shared, group, task);
      break;
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_SUBMIT) {
    if (((index_shard_staged_submit_status_t)callback_status ==
             INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED) !=
            (completion_id != 0ULL)) {
      group->internal_error = TRUE;
      index_shard_staged_set_failed_locked(shared, group, task);
      rc = -1;
    } else {
      switch ((index_shard_staged_submit_status_t)callback_status) {
      case INDEX_SHARD_STAGED_SUBMIT_RETRY: {
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_SUBMIT_READY);
        task->scheduler_epoch = claim->observed_completion_epoch;
        if (index_shard_staged_set_submit_wait_locked(
                group, task, TRUE) ||
            !group->pool ||
            index_shard_staged_park_uncredited_submitters_locked(
                group->pool)) {
          group->internal_error = TRUE;
          rc = -1;
        }
        shared->staged_submit_retries++;
        break;
      }
      case INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED: {
        anbool already_notified = FALSE;

        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
        task->scheduler_epoch = ~shared->staged_completion_epoch;
        (void)index_shard_staged_set_cancel_sent_locked(
            group, task, FALSE);
        task->completion_id = completion_id;
        if (index_shard_completion_registry_register_locked(
                shared, group, claim->task_index,
                completion_id, &already_notified)) {
          shared->completion_registry_error = TRUE;
          group->internal_error = TRUE;
          (void)index_shard_staged_set_completion_pending_locked(
              group, task, TRUE);
          rc = -1;
        } else if (index_shard_staged_set_completion_pending_locked(
                       group, task, already_notified)) {
          group->internal_error = TRUE;
          rc = -1;
        }
        task->scheduler_submit_seconds = now;
        group->io_submitted++;
        if (group->io_submitted > group->max_io_submitted) {
          group->max_io_submitted = group->io_submitted;
        }
        shared->staged_tickets_active++;
        shared->staged_max_io_submitted = MAX(
            shared->staged_max_io_submitted,
            shared->staged_tickets_active);
        shared->staged_source_leases++;
        shared->staged_io_submitted++;
        break;
      }
      case INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
        task->scheduler_ready_seconds = now;
        handoff_submit_credit = claim->submit_credit;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_OWNER_READY:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_OWNER_READY);
        task->scheduler_ready_seconds = now;
        handoff_submit_credit = claim->submit_credit;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_STOPPED:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
        group->stop_seen = TRUE;
        handoff_submit_credit = claim->submit_credit;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_ERROR:
      default:
        index_shard_staged_set_failed_locked(shared, group, task);
        handoff_submit_credit = claim->submit_credit;
        break;
      }
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_IO_POLL) {
    switch ((index_shard_staged_io_status_t)callback_status) {
    case INDEX_SHARD_STAGED_IO_PENDING:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
      if (claim->completion_inline) {
        /*
         * The provider dispatched this poll after terminal publication, so
         * PENDING cannot be followed by another completion notification.
         * Leave the ticket claimable for deterministic cancellation while
         * treating the broken provider contract as an integrity failure.
         */
        (void)index_shard_staged_set_completion_pending_locked(
            group, task, TRUE);
        shared->completion_registry_error = TRUE;
        group->internal_error = TRUE;
        rc = -1;
      }
      /* A worker poll retains a notifier racing the unlocked callback. */
      break;
    case INDEX_SHARD_STAGED_IO_READY:
      (void)index_shard_staged_release_ticket_locked(shared, group, task);
      task->completion_id = 0ULL;
      (void)index_shard_staged_set_completion_pending_locked(
          group, task, FALSE);
      if (task->scheduler_submit_seconds > 0.0 &&
          now >= task->scheduler_submit_seconds) {
        group->submit_to_ready_seconds +=
            now - task->scheduler_submit_seconds;
      }
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      } else {
        task->scheduler_ready_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_IO_FAILED:
      (void)index_shard_staged_release_ticket_locked(shared, group, task);
      task->completion_id = 0ULL;
      (void)index_shard_staged_set_completion_pending_locked(
          group, task, FALSE);
      if (task->scheduler_submit_seconds > 0.0 &&
          now >= task->scheduler_submit_seconds) {
        group->submit_to_ready_seconds +=
            now - task->scheduler_submit_seconds;
      }
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_OWNER_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      } else {
        task->scheduler_ready_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_IO_CANCELLED:
      (void)index_shard_staged_release_ticket_locked(shared, group, task);
      task->completion_id = 0ULL;
      (void)index_shard_staged_set_completion_pending_locked(
          group, task, FALSE);
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_IO_ERROR:
    default:
      (void)index_shard_staged_release_ticket_locked(shared, group, task);
      task->completion_id = 0ULL;
      (void)index_shard_staged_set_completion_pending_locked(
          group, task, FALSE);
      task->scheduler_submit_seconds = 0.0;
      index_shard_staged_set_failed_locked(shared, group, task);
      group->internal_error = TRUE;
      shared->completion_registry_error = TRUE;
      rc = -1;
      break;
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_IO_CANCEL) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
    task->scheduler_epoch = ~shared->staged_completion_epoch;
    if (callback_status < 0) {
      /* Poll once, then leave a still-pending cancellation retriable. */
      (void)index_shard_staged_set_cancel_sent_locked(
          group, task, FALSE);
      (void)index_shard_staged_set_completion_pending_locked(
          group, task, TRUE);
      group->internal_error = TRUE;
      rc = -1;
    } else {
      (void)index_shard_staged_set_cancel_sent_locked(
          group, task, TRUE);
    }
  } else {
    group->execute_seconds += callback_seconds;
    switch ((index_shard_staged_execute_status_t)callback_status) {
    case INDEX_SHARD_STAGED_EXECUTE_MORE:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_PREPARE_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      }
      break;
    case INDEX_SHARD_STAGED_EXECUTE_OK:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_RESULTS_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
        task->scheduler_result_seconds = 0.0;
      } else {
        task->scheduler_result_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_EXECUTE_STOPPED:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_EXECUTE_ERROR:
    default:
      index_shard_staged_set_failed_locked(shared, group, task);
      break;
    }
  }
  /*
   * A completion releases one provider admission before rearming one
   * waiter. If that waiter discovers already-completed pages or falls back
   * without submitting a ticket, the released admission remains unused.
   * Hand it to one more canonical waiter. A racing consumer can make that
   * waiter retry, but this remains a one-at-a-time chain rather than a herd.
   */
  if (handoff_submit_credit && group->pool) {
    rearm_status =
        index_shard_staged_rearm_one_submit_waiter_locked(
            group->pool, &rearmed_owner);
    if (rearm_status < 0) {
      shared->completion_registry_error = TRUE;
      group->internal_error = TRUE;
      rc = -1;
    } else if (rearm_status > 0) {
      shared->staged_submit_handoffs++;
    }
  }
  if (group->pool) {
    index_shard_staged_refresh_submit_backpressure_locked(
        group->pool);
  }
  if (!task_failed_before && group->task_failed) {
    logerr("[index-shard] staged callback failed ops=%s task=%zu "
           "claim=%i status=%i state=%i\n",
           group->ops && group->ops->name
               ? group->ops->name
               : "<unnamed>",
           claim->task_index,
           (int)claim->kind,
           callback_status,
           task->scheduler_state);
  }
  index_shard_notify_progress_locked(
      shared, group->owner_worker);
  if (rearm_status > 0 && rearmed_owner != group->owner_worker) {
    index_shard_owner_signal_locked(shared, rearmed_owner);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

static int index_shard_staged_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim) {
  index_shard_staged_group_t *group;
  index_shard_staged_task_t *task;
  int callback_status = INDEX_SHARD_STAGED_EXECUTE_ERROR;
  unsigned long long completion_id = 0ULL;
  double callback_start;
  double callback_seconds;

  if (!shared || !claim || !claim->group ||
      claim->task_index >= claim->group->task_count) {
    return -1;
  }
  group = claim->group;
  task = &group->tasks[claim->task_index];
  if (!group->ops) {
    return index_shard_staged_complete_claim(
        shared, claim, callback_status, 0.0, 0ULL);
  }
  callback_start = monotonic_seconds();
  switch (claim->kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    if (group->ops->prepare) {
      callback_status = group->ops->prepare(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    if (group->ops->submit) {
      callback_status = group->ops->submit(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes,
          &completion_id);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    if (group->ops->poll) {
      callback_status = group->ops->poll(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    if (group->ops->cancel) {
      callback_status = group->ops->cancel(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
    if (group->ops->execute) {
      callback_status = group->ops->execute(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    if (group->ops->owner) {
      index_shard_worker_context_t *ctx =
          index_shard_get_tls();

      if (!ctx || ctx->pool != group->pool ||
          ctx->worker_id != group->owner_worker ||
          ctx->published_staged_group != group ||
          ctx->staged_owner_callback_active ||
          ctx->staged_owner_callback_group) {
        callback_status =
            INDEX_SHARD_STAGED_EXECUTE_ERROR;
        break;
      }
      ctx->staged_owner_callback_active = TRUE;
      ctx->staged_owner_callback_group = group;
      callback_status = group->ops->owner(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
      ctx->staged_owner_callback_group = NULL;
      ctx->staged_owner_callback_active = FALSE;
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    break;
  }
  callback_seconds = monotonic_seconds() - callback_start;
  return index_shard_staged_complete_claim(
      shared, claim, callback_status, callback_seconds,
      completion_id);
}

int index_shard_inner_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_inner_claim_t *claim) {
  if (!shared || !claim) {
    return -1;
  }
  if (claim->kind == INDEX_SHARD_INNER_CLAIM_HELPER) {
    return index_shard_helper_execute_claim(
        shared, &claim->helper);
  }
  if (claim->kind == INDEX_SHARD_INNER_CLAIM_STAGED) {
    return index_shard_staged_execute_claim(
        shared, &claim->staged);
  }
  return -1;
}

/* queue_mutex must be held. */
static void index_shard_staged_cancel_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  uint64_t released_submit_credits;
  size_t i;

  if (!shared || !group || !group->tasks) {
    return;
  }
  released_submit_credits = group->submit_credit_mask;
  group->cancelling = TRUE;
  for (i = 0U; i < group->task_count; i++) {
    index_shard_staged_task_t *task = &group->tasks[i];

    switch ((index_shard_staged_task_state_t)
                task->scheduler_state) {
    case INDEX_SHARD_STAGED_TASK_PREPARE_READY:
    case INDEX_SHARD_STAGED_TASK_SUBMIT_READY:
    case INDEX_SHARD_STAGED_TASK_COMPUTE_READY:
    case INDEX_SHARD_STAGED_TASK_OWNER_READY:
    case INDEX_SHARD_STAGED_TASK_RESULTS_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      task->scheduler_result_seconds = 0.0;
      break;
    default:
      break;
    }
  }
  while (released_submit_credits && group->pool) {
    int rearmed_owner = -1;
    int rearm_status =
        index_shard_staged_rearm_one_submit_waiter_locked(
            group->pool, &rearmed_owner);

    if (rearm_status < 0) {
      group->internal_error = TRUE;
      break;
    }
    if (!rearm_status) {
      break;
    }
    index_shard_owner_signal_locked(shared, rearmed_owner);
    released_submit_credits &= released_submit_credits - UINT64_C(1);
  }
  if (group->pool) {
    index_shard_staged_refresh_submit_backpressure_locked(
        group->pool);
  }
}

/* queue_mutex must be held. */
static int index_shard_staged_cancel_for_pool_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  int fatal;
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested || shared->solved_published;
  pthread_mutex_unlock(&shared->state_mutex);
  if (fatal) {
    group->internal_error = TRUE;
    index_shard_staged_cancel_ready_locked(shared, group);
    return TRUE;
  }
  if (stopped || group->task_failed || group->stop_seen ||
      group->internal_error) {
    index_shard_staged_cancel_ready_locked(shared, group);
    return TRUE;
  }
  return FALSE;
}

/* queue_mutex must be held. */
static anbool index_shard_staged_all_terminal_locked(
    const index_shard_staged_group_t *group) {
  size_t i;

  if (!group || !group->tasks) {
    return FALSE;
  }
  for (i = 0U; i < group->task_count; i++) {
    if (!index_shard_staged_task_terminal(
            &group->tasks[i])) {
      return FALSE;
    }
  }
  return TRUE;
}

/*
 * Retire at most one complete logical task. Only the outer owner calls this
 * function. Slice preparation and execution may finish out of order, but the
 * callback sees exactly one final task in canonical task-index order.
 */
int index_shard_staged_retire_one(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  index_shard_staged_task_t *task;
  index_shard_staged_retire_status_t status;
  size_t task_index;
  double retire_start;
  double retire_seconds;

  if (!shared || !group || !group->retire) {
    return 1;
  }
  pthread_mutex_lock(&shared->queue_mutex);
  if (group->next_retire >= group->task_count) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  if (index_shard_staged_cancel_for_pool_locked(
          shared, group)) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task_index = group->next_retire;
  task = &group->tasks[task_index];
  if (task->scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  retire_start = monotonic_seconds();
  if (task->scheduler_result_seconds > 0.0 &&
      retire_start >= task->scheduler_result_seconds) {
    group->result_to_retire_seconds +=
        retire_start - task->scheduler_result_seconds;
  }
  task->scheduler_result_seconds = 0.0;
  (void)index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_RETIRING);
  pthread_mutex_unlock(&shared->queue_mutex);

  status = group->retire(task, task_index, group->owner_context);
  retire_seconds = monotonic_seconds() - retire_start;
  if (status != INDEX_SHARD_STAGED_RETIRE_OK &&
      status != INDEX_SHARD_STAGED_RETIRE_STOPPED &&
      status != INDEX_SHARD_STAGED_RETIRE_MORE &&
      status != INDEX_SHARD_STAGED_RETIRE_ERROR) {
    status = INDEX_SHARD_STAGED_RETIRE_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  group->retire_seconds += retire_seconds;
  if (task->scheduler_state !=
          INDEX_SHARD_STAGED_TASK_RETIRING ||
      group->next_retire != task_index) {
    group->internal_error = TRUE;
    status = INDEX_SHARD_STAGED_RETIRE_ERROR;
  } else if (status == INDEX_SHARD_STAGED_RETIRE_OK) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_RETIRED);
    group->next_retire++;
  } else if (status ==
             INDEX_SHARD_STAGED_RETIRE_MORE) {
    if (index_shard_staged_cancel_for_pool_locked(
            shared, group)) {
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      status = INDEX_SHARD_STAGED_RETIRE_STOPPED;
    } else {
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_PREPARE_READY);
    }
  } else if (status ==
             INDEX_SHARD_STAGED_RETIRE_STOPPED) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
    group->stop_seen = TRUE;
    index_shard_staged_cancel_ready_locked(shared, group);
  } else {
    index_shard_staged_set_failed_locked(shared, group, task);
    index_shard_staged_cancel_ready_locked(shared, group);
  }
  index_shard_notify_progress_locked(
      shared, group->owner_worker);
  pthread_mutex_unlock(&shared->queue_mutex);
  return status == INDEX_SHARD_STAGED_RETIRE_OK ||
      status == INDEX_SHARD_STAGED_RETIRE_MORE
      ? 0
      : -1;
}

/* queue_mutex must be held. */
size_t index_shard_staged_capacity(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t capacity = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active) {
    /* This is a storage bound, not an instantaneous idle-worker snapshot. */
    capacity = INDEX_SHARD_HELPER_MAX_TASKS;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return capacity;
}

size_t index_shard_staged_compute_width(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t width = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active) {
    width = (size_t)ctx->pool->worker_count;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return width;
}

index_shard_helper_run_status_t
index_shard_staged_run_ordered(
    const index_shard_staged_ops_t *ops,
    index_shard_staged_task_t *tasks,
    size_t task_count,
    index_shard_staged_retire_fn retire,
    void *owner_context,
    index_shard_staged_run_stats_t *stats) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  index_shard_staged_group_t *group;
  index_shard_helper_run_status_t result;
  size_t i;
  int fatal_requested = FALSE;
  int wait_broken = FALSE;

  if (stats) {
    memset(stats, 0, sizeof(*stats));
  }
  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered) {
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (!ops || !ops->prepare || !ops->submit || !ops->poll ||
      !ops->cancel || !ops->execute || !ops->owner || !retire ||
      !tasks || !task_count ||
      task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  for (i = 0U; i < task_count; i++) {
    if ((!tasks[i].input && tasks[i].input_bytes) ||
        (!tasks[i].output && tasks[i].output_bytes)) {
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
  }

  group = calloc(1, sizeof(*group));
  if (!group) {
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  group->pool = ctx->pool;
  group->ops = ops;
  group->tasks = tasks;
  group->task_count = task_count;
  group->retire = retire;
  group->owner_context = owner_context;
  group->generation = ctx->generation_seen;
  group->owner_epoch = ++ctx->staged_group_epoch;
  if (!group->owner_epoch) {
    group->owner_epoch = ++ctx->staged_group_epoch;
  }
  group->owner_worker = ctx->worker_id;
  group->owner_index_order = ctx->current_index_order;

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen != ctx->pool->generation ||
      !ctx->current_outer_active ||
      ctx->published_helper_group ||
      ctx->published_staged_group ||
      ctx->helper_preparation_active ||
      !ctx->pool->payload_completion_registered) {
    pthread_mutex_unlock(&shared->queue_mutex);
    free(group);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->fatal_error) {
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (shared->stop_requested || shared->solved_published) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  pthread_mutex_unlock(&shared->state_mutex);
  if (result != INDEX_SHARD_HELPER_OK) {
    pthread_mutex_unlock(&shared->queue_mutex);
    free(group);
    return result;
  }

  for (i = 0U; i < task_count; i++) {
    tasks[i].scheduler_state = INDEX_SHARD_STAGED_TASK_UNUSED;
    tasks[i].cancel_sent = FALSE;
    tasks[i].completion_pending = FALSE;
    tasks[i].callback_status = INDEX_SHARD_STAGED_EXECUTE_ERROR;
    tasks[i].scheduler_epoch = ~shared->staged_completion_epoch;
    tasks[i].completion_id = 0ULL;
    tasks[i].scheduler_submit_seconds = 0.0;
    tasks[i].scheduler_ready_seconds = 0.0;
    tasks[i].scheduler_result_seconds = 0.0;
    if (index_shard_staged_set_state_locked(
            shared, group, &tasks[i],
            INDEX_SHARD_STAGED_TASK_PREPARE_READY)) {
      group->internal_error = TRUE;
      break;
    }
  }
  if (group->internal_error ||
      group->prepare_ready_mask !=
          index_shard_staged_task_mask(task_count)) {
    pthread_mutex_unlock(&shared->queue_mutex);
    free(group);
    return INDEX_SHARD_HELPER_FATAL;
  }
  ctx->published_staged_group = group;
  shared->staged_groups_active++;
  shared->staged_groups_published++;
  index_shard_queue_signal_locked(shared);
  pthread_mutex_unlock(&shared->queue_mutex);
  fitsbin_payload_io_notify_wait_helpers();

  while (1) {
    index_shard_inner_claim_t claim;
    int selection;

    while (!index_shard_staged_retire_one(shared, group)) {
      /* Retire every complete logical packet in canonical order. */
    }

    pthread_mutex_lock(&shared->queue_mutex);
    (void)index_shard_staged_cancel_for_pool_locked(
        shared, group);
    if (group->internal_error || group->task_failed ||
        group->stop_seen) {
      index_shard_staged_cancel_ready_locked(shared, group);
    }

    if (group->internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }

    /*
     * Retirement is owner-only and runs without the queue lock. A foreign
     * completion can publish the next canonical result after the retirement
     * scan but before this lock is acquired. Recheck that predicate here so
     * the owner never sleeps after the corresponding broadcast has passed.
     */
    if (!group->cancelling &&
        group->next_retire < group->task_count &&
        group->tasks[group->next_retire].scheduler_state ==
            INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (!group->cancelling &&
        group->next_retire == group->task_count &&
        !group->running_count && !group->io_submitted) {
      break;
    }
    if (group->cancelling && !group->running_count &&
        !group->io_submitted &&
        index_shard_staged_all_terminal_locked(group)) {
      break;
    }

    selection = index_shard_inner_select_locked(
        ctx, shared, TRUE, &claim);
    if (selection < 0) {
      group->internal_error = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }
    if (!selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      if (index_shard_inner_execute_claim(shared, &claim)) {
        pthread_mutex_lock(&shared->queue_mutex);
        group->internal_error = TRUE;
        index_shard_queue_broadcast_locked(shared);
        pthread_mutex_unlock(&shared->queue_mutex);
      }
      continue;
    }

    if (wait_broken) {
      struct timespec pause = { 0, 1000000L };

      pthread_mutex_unlock(&shared->queue_mutex);
      nanosleep(&pause, NULL);
      continue;
    }
    {
      int wait_status;
      double wait_start;

      shared->staged_owner_wait_calls++;
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
      shared->staged_owner_wait_seconds +=
          monotonic_seconds() - wait_start;
      if (wait_status) {
        group->internal_error = TRUE;
        wait_broken = TRUE;
      }
    }
    pthread_mutex_unlock(&shared->queue_mutex);
  }

  if (group->running_count || group->compute_running ||
      group->io_submitted || group->compute_ready ||
      group->reorder_ready || group->completion_registry_entries ||
      group->prepare_ready_mask || group->submit_ready_mask ||
      group->submit_wait_mask || group->submit_credit_mask ||
      group->io_submitted_mask ||
      group->completion_pending_mask || group->cancel_sent_mask ||
      group->compute_ready_mask || group->owner_ready_mask ||
      group->results_ready_mask) {
    group->internal_error = TRUE;
  }
  for (i = 0U; i < task_count; i++) {
    if (tasks[i].completion_id || tasks[i].completion_pending) {
      group->internal_error = TRUE;
    }
  }
  if (!group->cancelling &&
      group->next_retire != group->task_count) {
    group->internal_error = TRUE;
  }
  {
    size_t groups_cleared = 0U;

    for (i = 0U; i < (size_t)shared->worker_count; i++) {
      if (ctx->pool->contexts[i].published_staged_group == group) {
        ctx->pool->contexts[i].published_staged_group = NULL;
        groups_cleared++;
      }
    }
    if (groups_cleared != 1U) {
      group->internal_error = TRUE;
    } else {
      shared->staged_groups_completed++;
      index_shard_staged_refresh_submit_backpressure_locked(
          ctx->pool);
    }
  }
  if (!shared->staged_groups_active) {
    group->internal_error = TRUE;
  } else {
    shared->staged_groups_active--;
  }
  shared->staged_max_io_submitted = MAX(
      shared->staged_max_io_submitted,
      group->max_io_submitted);
  shared->staged_max_compute_ready = MAX(
      shared->staged_max_compute_ready,
      group->max_compute_ready);
  shared->staged_max_reorder_ready = MAX(
      shared->staged_max_reorder_ready,
      group->max_reorder_ready);
  shared->staged_max_compute_running = MAX(
      shared->staged_max_compute_running,
      group->max_compute_running);
  shared->staged_prepare_claims += group->prepare_claims;
  shared->staged_submit_claims += group->submit_claims;
  shared->staged_poll_claims += group->poll_claims;
  shared->staged_inline_poll_claims += group->inline_poll_claims;
  shared->staged_execute_claims += group->execute_claims;
  shared->staged_owner_execute_claims +=
      group->owner_execute_claims;
  shared->staged_submit_to_ready_seconds +=
      group->submit_to_ready_seconds;
  shared->staged_ready_dwell_seconds +=
      group->ready_dwell_seconds;
  shared->staged_execute_seconds += group->execute_seconds;
  shared->staged_result_to_retire_seconds +=
      group->result_to_retire_seconds;
  shared->staged_retire_seconds += group->retire_seconds;
  index_shard_notify_progress_locked(
      shared, group->owner_worker);
  pthread_mutex_unlock(&shared->queue_mutex);

  if (stats) {
    stats->owner_claims = group->owner_claims;
    stats->foreign_claims = group->foreign_claims;
    stats->owner_compute_executes =
        group->owner_compute_executes;
    stats->foreign_compute_executes =
        group->foreign_compute_executes;
    stats->max_concurrent_claims = group->max_running;
    stats->max_compute_running = group->max_compute_running;
    stats->io_submitted = group->io_submitted + group->io_completed;
    stats->io_completed = group->io_completed;
    stats->max_io_submitted = group->max_io_submitted;
    stats->max_compute_ready = group->max_compute_ready;
    stats->max_reorder_ready = group->max_reorder_ready;
    stats->prepare_claims = group->prepare_claims;
    stats->submit_claims = group->submit_claims;
    stats->poll_claims = group->poll_claims;
    stats->execute_claims = group->execute_claims;
    stats->owner_claims_executed = group->owner_execute_claims;
    stats->submit_to_ready_seconds = group->submit_to_ready_seconds;
    stats->ready_dwell_seconds = group->ready_dwell_seconds;
    stats->execute_seconds = group->execute_seconds;
    stats->result_to_retire_seconds =
        group->result_to_retire_seconds;
    stats->retire_seconds = group->retire_seconds;
    stats->owner_work_units = group->owner_work;
    stats->foreign_work_units = group->foreign_work;
  }
  if (group->internal_error) {
    if (!fatal_requested) {
      index_shard_request_fatal_stop(shared);
    }
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (group->task_failed) {
    result = INDEX_SHARD_HELPER_TASK_FAILED;
  } else if (group->stop_seen || group->cancelling) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  free(group);
  return result;
}
