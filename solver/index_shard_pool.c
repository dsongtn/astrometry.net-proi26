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

static index_shard_pool_t *index_shard_global_pool = NULL;
static pthread_mutex_t index_shard_global_pool_mutex =
    PTHREAD_MUTEX_INITIALIZER;

/*
 * SECTION INDEX-SHARD: pool
 *
 * Pool lifecycle and pass submission.
 *
 * The pool is created once per engine job and reused across onefield_run()
 * calls.  Each submitted pass increments generation to wake workers.
 */

// ANCHOR INDEX-SHARD: shared-init
/*
 * Initialize synchronization primitives for the reusable shared pass state.
 */
static int index_shard_shared_init(index_shard_thread_state_t *shared) {
  pthread_condattr_t condattr;

  memset(shared, 0, sizeof(index_shard_thread_state_t));

  if (pthread_mutex_init(&shared->queue_mutex, NULL)) {
    return -1;
  }
  if (pthread_cond_init(&shared->queue_cv, NULL)) {
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_mutex_init(&shared->result_mutex, NULL)) {
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_condattr_init(&condattr)) {
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  if (pthread_condattr_setclock(
          &condattr,
          CLOCK_MONOTONIC)) {
    pthread_condattr_destroy(&condattr);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  if (pthread_cond_init(
          &shared->result_cv,
          &condattr)) {
    pthread_condattr_destroy(&condattr);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  pthread_condattr_destroy(&condattr);

  if (pthread_mutex_init(&shared->state_mutex, NULL)) {
    pthread_cond_destroy(&shared->result_cv);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_mutex_init(&shared->limit_mutex, NULL)) {
    pthread_mutex_destroy(&shared->state_mutex);
    pthread_cond_destroy(&shared->result_cv);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  return 0;
}

// ANCHOR INDEX-SHARD: shared-destroy
/*
 * Destroy synchronization primitives after all workers have joined.
 */
static void index_shard_shared_destroy(index_shard_thread_state_t *shared) {
  index_shard_completion_registry_destroy(shared);
  pthread_cond_destroy(&shared->queue_cv);
  pthread_mutex_destroy(&shared->queue_mutex);

  pthread_mutex_destroy(&shared->result_mutex);
  pthread_cond_destroy(&shared->result_cv);

  pthread_mutex_destroy(&shared->state_mutex);

  pthread_mutex_destroy(&shared->limit_mutex);
}


/*
 * Reserve the persistent pool for one submitted pass.
 *
 * Lock order:
 *   index_shard_global_pool_mutex -> pool->control_mutex
 *
 * This prevents a second submission from overwriting shared pass pointers and
 * also keeps pool_stop() from destroying the pool while the caller still uses
 * it.
 */
index_shard_pool_acquire_status_t
index_shard_pool_acquire_pass(onefield_t *bp,
                              solver_t *sp,
                              index_shard_pool_t **pool_out) {
  index_shard_pool_t *pool;

  if (!pool_out) {
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  *pool_out = NULL;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE;
  }

  if (pool->owner_bp != bp || pool->owner_sp != sp) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (pool->shutdown || pool->stopping || pool->pass_active) {
    pthread_mutex_unlock(&pool->control_mutex);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  pool->pass_active = TRUE;
  *pool_out = pool;

  pthread_mutex_unlock(&pool->control_mutex);
  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  return INDEX_SHARD_POOL_ACQUIRE_OK;
}

/*
 * Release one pass reservation after all caller-side work that dereferences
 * the pool has completed.
 */
void index_shard_pool_release_pass(index_shard_pool_t *pool) {
  if (!pool) {
    return;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (!pool->pass_active) {
    logerr("[index-shard] pass release requested with no active pass\n");
  } else {
    pool->pass_active = FALSE;
  }

  /*
   * pool_stop() may be waiting for pass_active to become false. Worker waiters
   * also use work_cv, but they re-check their generation predicate in a loop.
   */
  pthread_cond_broadcast(&pool->work_cv);

  pthread_mutex_unlock(&pool->control_mutex);
}

static void index_shard_context_owner_cvs_destroy(
    index_shard_pool_t *pool) {
  int i;

  if (!pool || !pool->contexts) {
    return;
  }
  for (i = 0; i < pool->worker_count; i++) {
    if (!pool->contexts[i].owner_cv_ready) {
      continue;
    }
    pthread_cond_destroy(&pool->contexts[i].owner_cv);
    pool->contexts[i].owner_cv_ready = FALSE;
    pool->contexts[i].owner_waiting = FALSE;
    pool->contexts[i].owner_wake_pending = FALSE;
  }
}

// ANCHOR INDEX-SHARD: pool-start
/*
 * Create persistent worker pool.
 *
 * Workers are created once and sleep until the first pass is submitted.
 */
int index_shard_pool_start(onefield_t *bp, solver_t *sp) {
  index_shard_pool_t *pool;
  int i;
  int tls_status;
  int worker_count;
  int payload_io_lanes;
  int payload_io_width;
  index_shard_width_plan_t width_plan;

   // pool already active for this engine job
  if (!index_shard_pthread_enabled(bp)) {
    return 0;
  }

  if (!bp || !sp) {
    ERROR("Cannot start index-shard pool without owner state");
    return -1;
  }

  tls_status = index_shard_tls_ensure();
  if (tls_status) {
    logerr("[index-shard] cannot initialize worker TLS status=%i\n",
           tls_status);
    return -1;
  }

  /* Initialize process-wide lazy read-only state before workers start. */
  (void)fits_get_endian_string();

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  if (index_shard_global_pool) {
    int reusable;

    pool = index_shard_global_pool;

    if (pool->owner_bp != bp || pool->owner_sp != sp) {
      logerr("[index-shard] global pool already belongs to another engine job\n");
      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }

    pthread_mutex_lock(&pool->control_mutex);
    reusable = !pool->shutdown && !pool->stopping;
    pthread_mutex_unlock(&pool->control_mutex);

    pthread_mutex_unlock(&index_shard_global_pool_mutex);

    if (!reusable) {
      logerr("[index-shard] owner pool is stopping or shut down\n");
      return -1;
    }

    return 0;
  }

  worker_count = index_shard_get_worker_count(bp, 0);

  pool = calloc(1, sizeof(index_shard_pool_t));
  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    SYSERROR("Failed to allocate index-shard pool");
    return -1;
  }

  pool->owner_bp = bp;
  pool->owner_sp = sp;
  pool->worker_count = worker_count;
  pool->producer_width = 1U;
  pool->helper_width = 0U;
  pool->inverse_cache_budget =
      index_shard_inverse_cache_budget();

  if (pthread_mutex_init(&pool->inverse_cache_mutex, NULL)) {
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  // initialize shared state before workers can observe pool
  if (pthread_mutex_init(&pool->control_mutex, NULL)) {
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  if (pthread_cond_init(&pool->work_cv, NULL)) {
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  if (index_shard_shared_init(&pool->shared)) {
    pthread_cond_destroy(&pool->work_cv);
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }
  pool->shared.pool = pool;
  if (index_shard_completion_registry_init(
          &pool->shared, worker_count)) {
    index_shard_shared_destroy(&pool->shared);
    pthread_cond_destroy(&pool->work_cv);
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  pool->threads = calloc((size_t)worker_count, sizeof(pthread_t));
  pool->contexts = calloc((size_t)worker_count,
                          sizeof(index_shard_worker_context_t));

    if (!pool->threads || !pool->contexts) {
      free(pool->threads);
      free(pool->contexts);

      index_shard_shared_destroy(&pool->shared);

      pthread_cond_destroy(&pool->work_cv);
      pthread_mutex_destroy(&pool->control_mutex);
      pthread_mutex_destroy(&pool->inverse_cache_mutex);

      free(pool);

      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }

  // worker contexts and owner conditions are stable for pool lifetime.
  for (i = 0; i < worker_count; i++) {
    pool->contexts[i].worker_id = i;
    pool->contexts[i].generation_seen = 0;
    pool->contexts[i].pool = pool;
    if (pthread_cond_init(&pool->contexts[i].owner_cv, NULL)) {
      index_shard_context_owner_cvs_destroy(pool);
      free(pool->threads);
      free(pool->contexts);
      index_shard_shared_destroy(&pool->shared);
      pthread_cond_destroy(&pool->work_cv);
      pthread_mutex_destroy(&pool->control_mutex);
      pthread_mutex_destroy(&pool->inverse_cache_mutex);
      free(pool);
      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }
    pool->contexts[i].owner_cv_ready = TRUE;
  }

  for (i = 0; i < worker_count; i++) {
    if (pthread_create(
            &pool->threads[i], NULL,
            index_shard_worker_main, &pool->contexts[i])) {
      int j;

      pthread_mutex_lock(&pool->control_mutex);
      pool->shutdown = TRUE;
      pthread_cond_broadcast(&pool->work_cv);
      pthread_mutex_unlock(&pool->control_mutex);
      for (j = 0; j < i; j++) {
        pthread_join(pool->threads[j], NULL);
      }
      index_shard_context_owner_cvs_destroy(pool);
      free(pool->threads);
      free(pool->contexts);
      index_shard_shared_destroy(&pool->shared);
      pthread_cond_destroy(&pool->work_cv);
      pthread_mutex_destroy(&pool->control_mutex);
      pthread_mutex_destroy(&pool->inverse_cache_mutex);
      free(pool);
      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }
  }

  pthread_mutex_lock(&pool->control_mutex);
  while (pool->ready_workers < worker_count) {
    pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
  }
  tls_status = pool->tls_startup_error;
  if (tls_status) {
    pool->shutdown = TRUE;
    pthread_cond_broadcast(&pool->work_cv);
  }
  pthread_mutex_unlock(&pool->control_mutex);

  if (tls_status) {
    logerr("[index-shard] worker TLS startup failed status=%i\n",
           tls_status);
    for (i = 0; i < worker_count; i++) {
      pthread_join(pool->threads[i], NULL);
    }
    free(pool->threads);
    index_shard_context_owner_cvs_destroy(pool);
    free(pool->contexts);
    index_shard_shared_destroy(&pool->shared);
    pthread_cond_destroy(&pool->work_cv);
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  /*
   * Completion lanes populate only packet-planned mapped pages. They are not
   * compute workers and do not own solver, index, or result state.
   */
  payload_io_lanes = worker_count;
  if (!fitsbin_payload_io_service_width()) {
    if (fitsbin_payload_io_service_start(payload_io_lanes)) {
      logverb("[index-shard] mapped-page completion unavailable; "
              "using native mmap demand\n");
    } else {
      pool->payload_io_owned = TRUE;
    }
  }
  if (fitsbin_payload_io_service_width() &&
      fitsbin_payload_io_mapped_population_supported() &&
      !fitsbin_payload_io_set_completion_notifier(
          index_shard_staged_completion_notify, pool)) {
    pool->payload_completion_registered = TRUE;
  } else if (fitsbin_payload_io_service_width() &&
             fitsbin_payload_io_mapped_population_supported()) {
    logverb("[index-shard] detached payload completion unavailable; "
            "retaining fixed helper reservation\n");
  } else if (fitsbin_payload_io_service_width()) {
    logverb("[index-shard] mapped-page population unsupported; "
            "retaining fixed helper reservation\n");
  }

  index_shard_global_pool = pool;
  payload_io_width = fitsbin_payload_io_service_width();
  if (index_shard_config_plan_widths(
          worker_count,
          payload_io_width,
          pool->payload_completion_registered,
          FALSE,
          &width_plan)) {
    logerr("[index-shard] invalid compute/delivery width plan "
           "workers=%i payload_io_width=%i detached=%i\n",
           worker_count,
           payload_io_width,
           pool->payload_completion_registered ? 1 : 0);
    pool->helper_width = worker_count > 1 ? 1U : 0U;
    pool->producer_width =
        (size_t)worker_count - pool->helper_width;
  } else {
    /*
     * Detached completion does not reserve a compute lane. Every worker may
     * own outer work while it is immediately available, then join published
     * staged work after the outer queue is exhausted.
     */
    pool->producer_width = width_plan.producer_width;
    pool->helper_width = width_plan.helper_width;
  }
  fitsbin_payload_io_configure_workers(
      pool->payload_completion_registered
          ? (int)pool->producer_width
          : worker_count);

  logverb("[index-shard] workers=%i mode=pthread "
          "compute_width=%i producer_width=%zu helper_width=%zu "
          "payload_io_width=%i inverse_cache_budget=%zu\n",
          worker_count,
          worker_count,
          pool->producer_width,
          pool->helper_width,
          payload_io_width,
          pool->inverse_cache_budget);

  pthread_mutex_unlock(&index_shard_global_pool_mutex);
  return 0;
}
// ANCHOR INDEX-SHARD: pool-stop
/*
 * Stop the pool belonging to bp and join all workers.
 *
 * New passes are rejected as soon as stopping is published. If one pass is
 * already active, shutdown waits for that pass reservation to be released.
 */
void index_shard_pool_stop(onefield_t *bp) {
  index_shard_pool_t *pool;
  int i;
  int notifier_clear_failed = FALSE;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  if (pool->owner_bp != bp) {
    logerr("[index-shard] refusing to stop pool owned by another engine job\n");
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (pool->stopping) {
    pthread_mutex_unlock(&pool->control_mutex);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  pool->stopping = TRUE;
  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  while (pool->pass_active) {
    pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
  }

  pool->shutdown = TRUE;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  if (index_shard_global_pool == pool) {
    index_shard_global_pool = NULL;
  }

  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  for (i = 0; i < pool->worker_count; i++) {
    pthread_join(pool->threads[i], NULL);
  }

  if (pool->payload_completion_registered) {
    if (fitsbin_payload_io_clear_completion_notifier(
            index_shard_staged_completion_notify, pool)) {
      logerr("[index-shard] failed to clear payload completion notifier\n");
      notifier_clear_failed = TRUE;
    } else {
      pool->payload_completion_registered = FALSE;
    }
  }
  if (pool->payload_io_owned) {
    fitsbin_payload_io_service_stop();
  }
  if (notifier_clear_failed) {
    /*
     * The process-wide notifier may still retain pool as opaque state, or an
     * external clear may still be draining an active callback. Preserve the
     * stopped allocation rather than risk a callback use-after-free.
     */
    logerr("[index-shard] retaining stopped pool after notifier "
           "teardown failure\n");
    return;
  }
  free(pool->threads);
  index_shard_context_owner_cvs_destroy(pool);
  free(pool->contexts);
  index_shard_shared_destroy(&pool->shared);

  logverb("[index-shard] inverse-cache hits=%llu misses=%llu "
          "admitted=%llu refused=%llu evicted=%llu "
          "overcommit=%llu retained=%zu active=%zu "
          "retained_peak=%zu combined_peak=%zu budget=%zu\n",
          pool->inverse_cache_hits,
          pool->inverse_cache_misses,
          pool->inverse_cache_admitted,
          pool->inverse_cache_refused,
          pool->inverse_cache_evicted,
          pool->inverse_cache_overcommit,
          pool->inverse_cache_bytes,
          pool->inverse_active_bytes,
          pool->inverse_cache_peak_bytes,
          pool->inverse_combined_peak_bytes,
          pool->inverse_cache_budget);
  index_shard_inverse_cache_destroy(pool);
  pthread_mutex_destroy(&pool->inverse_cache_mutex);
  pthread_cond_destroy(&pool->work_cv);
  pthread_mutex_destroy(&pool->control_mutex);

  logverb("[index-shard] pthread-pool stop\n");

  free(pool);
}

// ANCHOR INDEX-SHARD: pool-active
/*
 * Return true only when the compatibility pool belongs to bp and remains
 * available for a new pass reservation.
 */
int index_shard_pool_active(onefield_t *bp) {
  index_shard_pool_t *pool;
  int active = FALSE;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (pool && pool->owner_bp == bp) {
    pthread_mutex_lock(&pool->control_mutex);
    active = !pool->shutdown && !pool->stopping;
    pthread_mutex_unlock(&pool->control_mutex);
  }

  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  return active;
}
