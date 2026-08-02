/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/*
 * Developer navigation: worker and outer-owner execution
 * -----------------------------------------------------
 * Worker threads enter here after the scheduler assigns a claim. For an outer
 * claim, one context acquires one index, builds worker-private onefield and
 * solver state, runs the native solver flow, freezes its result, and releases
 * the index through the matching hook. Inner claims execute only their bounded
 * immutable package and return it to the original owner.
 *
 * There is intentionally no persistent index_t cache in this layer. Every
 * exit path must clear TLS/advice state and produce exactly one terminal task
 * transition before the context can claim work from another generation.
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
 * Load one index for one shard task through original onefield hooks.
 *
 * No persistent index_t cache here.  Full-index caching caused unacceptable
 * RSS growth because candidate sets can contain hundreds of heavy indexes.
 */
static index_shard_hook_result_t index_shard_worker_get_index(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared,
    size_t index_order,
    index_t **index_out) {
  index_shard_hook_result_t hook_result = {
      INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE, -1};

  if (index_out) {
    *index_out = NULL;
  }
  if (!shared || !shared->hooks ||
      !shared->hooks->get_index || !index_out) {
    return hook_result;
  }

  hook_result = shared->hooks->get_index(
      shared->bp, index_order, index_out);
  if (hook_result.outcome ==
          INDEX_SHARD_HOOK_COMPLETED_UNSOLVED &&
      !*index_out) {
    hook_result.outcome =
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE;
    hook_result.error_code = -1;
  }

  if (index_shard_trace_enabled() && *index_out) {
    logmsg("[index-shard] worker=%i load index_order=%zu index=%s\n",
           ctx->worker_id, index_order,
           (*index_out)->indexname
               ? (*index_out)->indexname
               : "(null)");
  }

  return hook_result;
}

static int index_shard_apply_index_mmap_advice(
    index_t* index,
    fitsbin_mmap_advice_t advice) {
    fitsbin_t* fb;
    int failures = 0;

    if (!index) {
        return 0;
    }

    /*
     * Apply the component policy to existing mappings. fitsbin keeps compact
     * topology NORMAL and applies the pass policy only to sparse payload.
     */
    if (index->codekd &&
        index->codekd->tree &&
        index->codekd->tree->io) {
        fb = (fitsbin_t*)index->codekd->tree->io;

        if (fitsbin_set_mmap_advice(
                fb,
                advice,
                TRUE)) {
            failures++;
        }
    }

    /*
     * Star KD tree.
     */
    if (index->starkd &&
        index->starkd->tree &&
        index->starkd->tree->io) {
        fb = (fitsbin_t*)index->starkd->tree->io;

        if (fitsbin_set_mmap_advice(
                fb,
                advice,
                TRUE)) {
            failures++;
        }
    }

    /*
     * Quad table.
     */
    if (index->quads &&
        index->quads->fb) {
        if (fitsbin_set_mmap_advice(
                index->quads->fb,
                advice,
                TRUE)) {
            failures++;
        }
    }
    return failures;
}

/*
 * SECTION INDEX-SHARD: worker-context
 *
 * Worker-local context reuse.
 *
 * local_bp/local solver are prepared once per submitted pass and reused across
 * many one-index tasks.  This removes repeated xylist open/close and solver
 * index-list allocation from the per-index hot path.
 */
static int index_shard_worker_prepare_pass(index_shard_worker_context_t *ctx,
                                           index_shard_thread_state_t *shared) {
  // old local context belongs to a previous generation -> cleanup first
  if (ctx->local_context_ready && ctx->local_context_generation == ctx->generation_seen)
    return 0;

  if (ctx->local_context_ready) {
    if (shared->hooks && shared->hooks->cleanup_local_context)
      shared->hooks->cleanup_local_context(&ctx->local_bp);

    ctx->local_context_ready = FALSE;
  }

  if (!shared->hooks || !shared->hooks->prepare_local_context)
    return -1;

  if (!shared->worker_view ||
      shared->hooks->prepare_local_context(
          &ctx->local_bp, shared->worker_view))
    return -1;
  // hook copies stable master config + opens worker-local xylist
  ctx->local_context_ready = TRUE;
  ctx->local_context_generation = ctx->generation_seen;

  return 0;
}
// ANCHOR INDEX-SHARD: worker-cleanup-pass
/*
 * Release worker-local pass context after this generation is finished.
 *
 * Does not touch master bp and does not free result slots.
 */
void index_shard_worker_cleanup_pass(index_shard_worker_context_t *ctx,
                                            index_shard_thread_state_t *shared) {
  if (!ctx->local_context_ready)
    return;

  if (shared->hooks && shared->hooks->cleanup_local_context)
    shared->hooks->cleanup_local_context(&ctx->local_bp);

  memset(&ctx->local_bp, 0, sizeof(onefield_t));
  ctx->local_context_ready = FALSE;
  ctx->local_context_generation = 0;
}

static index_shard_hook_result_t index_shard_done_with_index(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared,
    size_t index_order,
    index_t *index,
    index_shard_inverse_lease_t *inverse_lease) {
  index_shard_hook_result_t hook_result = {
      INDEX_SHARD_HOOK_COMPLETED_UNSOLVED, 0};

  if (!index) {
    return hook_result;
  }

  index_shard_inverse_cache_release(
      ctx, index, inverse_lease);
  if (!shared || !shared->hooks ||
      !shared->hooks->done_with_index) {
    hook_result.outcome =
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE;
    hook_result.error_code = -1;
    return hook_result;
  }

  return shared->hooks->done_with_index(
      shared->bp, index_order, index);
}

// ANCHOR INDEX-SHARD: run-one-index
/*
 * Execute one index shard in one worker.
 *
 * This function owns only local computation:
 *   - reset local context for this result slot
 *   - load one index
 *   - run solver against local_bp
 *   - analyze local solutions
 *   - release index
 *
 * It does not merge into master bp.
 */
static int index_shard_run_one_with_worker_context(index_shard_worker_context_t *ctx,
                                                   index_shard_thread_state_t *shared,
                                                   size_t index_order,
                                                   index_shard_result_t *result,
                                                   fitsbin_mmap_advice_t mmap_advice) {
  // one result slot belongs to this task
  index_t *index = NULL;
  double task_wall_start;
  double phase_wall_start;
  double wall_start;
  float cpu_start;
  index_shard_inverse_lease_t inverse_lease;
  index_shard_hook_result_t hook_result;
  int hook_status;

  index_shard_result_init(result, index_order);
  memset(&inverse_lease, 0, sizeof(inverse_lease));
  result->mmap_advice = mmap_advice;

  if (!result->solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_TASK_LOCAL,
        -1);
    return -1;
  }

  if (!ctx->local_context_ready) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  if (!shared->hooks || !shared->hooks->reset_local_context_for_task ||
      !shared->hooks->solve_one_index) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  /*
   * Full outer-task timing starts before local reset and index acquisition.
   * Queue wait and reducer wait are deliberately excluded.
   */
  task_wall_start = monotonic_seconds();

  result->task_started = TRUE;
  result->worker_id = ctx->worker_id;
  result->task_start_since_pass =
      task_wall_start - shared->pass_wall_start;
#if defined(RUSAGE_THREAD)
  result->task_resource_valid =
      (getrusage(
           RUSAGE_THREAD,
           &result->task_resource_start) == 0);
#else
  result->task_resource_valid = FALSE;
#endif

   // local_bp reused across tasks, but solutions change per task
  phase_wall_start = monotonic_seconds();

  shared->hooks->reset_local_context_for_task(&ctx->local_bp, result->solutions);

  result->reset_seconds = monotonic_seconds() - phase_wall_start;

  phase_wall_start = monotonic_seconds();

  /*
   * get_index() opens and mmaps the index components. Install the immutable
   * pass advice before acquisition so every new mapping receives the correct
   * policy immediately.
   */
  fitsbin_mmap_set_thread_advice(
      result->mmap_advice);

  hook_result = index_shard_worker_get_index(
      ctx,
      shared,
      index_order,
      &index);
  hook_status = index_shard_apply_hook_result(
      result, hook_result, FALSE);

  fitsbin_mmap_clear_thread_advice();

  result->acquire_seconds =
      monotonic_seconds() - phase_wall_start;

  if (hook_status || !index) {
    if (!hook_status && !index) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }
    if (result->failed) {
      ERROR("Failed to load index order %zu", index_order);
    }
    if (index_shard_arbitrate_candidate(
            shared, index_order)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }

    index_shard_result_finish_task(
        result, shared, task_wall_start);
    return result->failed ? -1 : 0;
  }

  /*
   * Reapply to any component that the onefield hook reused rather than opened
   * during this acquisition.
   */
  {
    int advice_failures =
        index_shard_apply_index_mmap_advice(
            index,
            result->mmap_advice);

    if (advice_failures > 0) {
      __atomic_add_fetch(
          &shared->mmap_advice_failures,
          (unsigned long long)advice_failures,
          __ATOMIC_RELAXED);
      logerr("[index-shard] failed to apply mmap advice "
             "index_order=%zu components=%i\n",
             index_order,
             advice_failures);
    }
  }
  index_shard_inverse_cache_attach(
      ctx, index, &inverse_lease);

  result->acquire_seconds =
      monotonic_seconds() - phase_wall_start;

  /*
   * Another group can solve or exhaust the aggregate budget while this owner
   * is opening an index. Recheck before faulting solver payload.
   */
  if (index_shard_check_global_limits(shared) ||
      index_shard_master_stop_requested(shared)) {
    index_shard_master_limit_or_cancel_requested(
        shared,
        &result->hit_total_cpulimit,
        &result->hit_total_timelimit,
        &result->cancelled);

    phase_wall_start = monotonic_seconds();
    hook_result = index_shard_done_with_index(
        ctx,
        shared,
        index_order,
        index,
        &inverse_lease);
    (void)index_shard_apply_hook_result(
        result, hook_result, FALSE);
    index = NULL;
    result->release_seconds =
        monotonic_seconds() - phase_wall_start;

    if (index_shard_arbitrate_candidate(
            shared, index_order)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }
    index_shard_result_finish_task(
        result,
        shared,
        task_wall_start);
    return result->failed ? -1 : 0;
  }

  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] worker=%i start index_order=%zu index=%s\n",
           ctx->worker_id,
           index_order,
           index->indexname ? index->indexname : "(null)");
  }
  // time only the actual one-index solve section
  wall_start = monotonic_seconds();
  cpu_start = get_cpu_usage();

  // Worker-lifetime TLS lets onefield callbacks poll this pass for stop.
  ctx->current_index_order = index_order;
  ctx->current_outer_active = TRUE;
  hook_result = shared->hooks->solve_one_index(
      &ctx->local_bp, index);
  ctx->current_outer_active = FALSE;

  result->wall_seconds = monotonic_seconds() - wall_start;
  result->cpu_seconds = get_cpu_usage() - cpu_start;

  result->hit_total_cpulimit =
      ctx->local_bp.hit_total_cpulimit;
  result->hit_total_timelimit =
      ctx->local_bp.hit_total_timelimit;
  result->cancelled = ctx->local_bp.cancelled;
  result->solver_profile = ctx->local_bp.solver.profile;
  hook_status = index_shard_apply_hook_result(
      result, hook_result, FALSE);

  phase_wall_start = monotonic_seconds();
  if (!hook_status && !result->failed) {
    hook_status = index_shard_capture_solution_analysis(
        shared, result);
  }
  result->analyze_seconds =
      monotonic_seconds() - phase_wall_start;

  if (index_shard_arbitrate_candidate(
          shared, index_order)) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
  }

  /*
   * Log the index name before done_with_index() releases index ownership.
   */
  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] worker=%i finish index_order=%zu index=%s\n",
           ctx->worker_id,
           index_order,
           index->indexname ? index->indexname : "(null)");
  }

  // Release through the original onefield ownership hook.
  phase_wall_start = monotonic_seconds();
  hook_result = index_shard_done_with_index(
      ctx,
      shared,
      index_order,
      index,
      &inverse_lease);
  (void)index_shard_apply_hook_result(
      result, hook_result, FALSE);
  index = NULL;

  result->release_seconds =
      monotonic_seconds() - phase_wall_start;

  index_shard_result_finish_task(
      result, shared, task_wall_start);

  return result->failed ? -1 : 0;
}
// ANCHOR INDEX-SHARD: worker-done
/*
 * Publish that this worker is done with the submitted pass.
 *
 * Reducer waits on active_workers reaching zero before final cleanup/drain.
 */
static void index_shard_worker_done(index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->result_mutex);

  if (shared->active_workers > 0)
    shared->active_workers--;

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);
}

/*
 * A compute worker waiting for one mapped-page completion may execute one
 * already-published coarse helper task. This never claims an outer index and
 * never changes index ownership. The ticket remains responsible for waking
 * the worker when its own pages become ready.
 */
static int index_shard_payload_wait_stop(void *opaque) {
  index_shard_worker_context_t *ctx = opaque;

  if (!ctx || !ctx->pool) {
    return TRUE;
  }
  return index_shard_worker_stop_requested();
}

static int index_shard_payload_wait_help(void *opaque) {
  index_shard_worker_context_t *ctx = opaque;
  index_shard_thread_state_t *shared;
  index_shard_inner_claim_t claim;
  int selection;

  if (!ctx || !ctx->pool ||
      index_shard_worker_stop_requested()) {
    return FALSE;
  }
  shared = &ctx->pool->shared;
  memset(&claim, 0, sizeof(claim));

  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->published_helper_group) {
    index_shard_helper_group_t *group =
        ctx->published_helper_group;
    size_t reserved = 0U;

    if (group->foreign_reserve > group->foreign_claims) {
      reserved = group->foreign_reserve -
          group->foreign_claims;
    }
    if (!group->ready_count ||
        group->ready_count <= reserved) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return FALSE;
    }
    claim.kind = INDEX_SHARD_INNER_CLAIM_HELPER;
    claim.helper.group = group;
    selection = index_shard_helper_owner_claim_locked(
        shared, group, &claim.helper.task_index);
    if (!selection) {
      group->owner_claims++;
      group->owner_work += index_shard_helper_task_work(
          &group->tasks[claim.helper.task_index]);
      shared->helper_tasks_owner++;
    }
  } else {
    selection = index_shard_inner_select_locked(
        ctx, shared, TRUE, &claim);
  }
  pthread_mutex_unlock(&shared->queue_mutex);

  if (selection < 0) {
    index_shard_request_fatal_stop(shared);
    return FALSE;
  }
  if (selection > 0) {
    return FALSE;
  }
  if (index_shard_inner_execute_claim(shared, &claim)) {
    index_shard_request_fatal_stop(shared);
    return FALSE;
  }
  return TRUE;
}

// ANCHOR INDEX-SHARD: worker-main
/*
 * Persistent worker loop.
 *
 * Worker sleeps until pool generation changes, prepares local pass context,
 * claims one-index tasks, then cleans up local context when the pass ends.
 */
void *index_shard_worker_main(void *userdata) {
  index_shard_worker_context_t *ctx = userdata;
  index_shard_pool_t *pool;
  int tls_status;

  if (!ctx) {
    return NULL;
  }

  pool = ctx->pool;
  if (!pool) {
    return NULL;
  }

  /*
   * Install the immutable worker context before declaring this thread ready.
   * A failed pthread TLS setup is a pool-start failure, never a worker that
   * silently runs without global stop/cancellation visibility.
   */
  tls_status = index_shard_set_tls(ctx);
  if (!tls_status) {
    tls_status = fitsbin_payload_io_set_thread_wait_helper(
        index_shard_payload_wait_help,
        index_shard_payload_wait_stop,
        ctx);
  }
  pthread_mutex_lock(&pool->control_mutex);
  if (tls_status && !pool->tls_startup_error) {
    pool->tls_startup_error = tls_status;
  }
  pool->ready_workers++;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);
  if (tls_status) {
    return NULL;
  }

  while (1) {
    index_shard_thread_state_t *shared = NULL;
    size_t index_order = 0U;
    int run_pass = FALSE;

    pthread_mutex_lock(&pool->control_mutex);

    while (!pool->shutdown &&
           ctx->generation_seen == pool->generation) {
      pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
    }

    if (pool->shutdown) {
      pthread_mutex_unlock(&pool->control_mutex);
      break;
    }

    ctx->generation_seen = pool->generation;
    shared = &pool->shared;

    if (ctx->worker_id < shared->worker_count) {
      run_pass = TRUE;
    }

    pthread_mutex_unlock(&pool->control_mutex);

    if (!run_pass) {
      continue;
    }

    ctx->pass_prepare_seconds = 0.0;
    ctx->pass_cleanup_seconds = 0.0;
    ctx->current_outer_active = FALSE;
    ctx->current_index_order = 0U;
    ctx->ready_before_outer_eligible = FALSE;

    while (1) {
      index_shard_result_t *result = NULL;
      index_shard_inner_claim_t inner_claim;
      fitsbin_mmap_advice_t mmap_advice =
          FITSBIN_MMAP_ADVICE_NORMAL;
      index_shard_work_selection_t selection;

      memset(&inner_claim, 0, sizeof(inner_claim));
      selection = index_shard_select_work(
          ctx,
          shared,
          &index_order,
          &mmap_advice,
          &inner_claim);

      if (selection == INDEX_SHARD_WORK_ERROR) {
        index_shard_request_fatal_stop(shared);
        break;
      }
      if (selection == INDEX_SHARD_WORK_DONE) {
        break;
      }
      if (selection == INDEX_SHARD_WORK_HELPER) {
        if (index_shard_inner_execute_claim(
                shared, &inner_claim)) {
          index_shard_request_fatal_stop(shared);
          break;
        }
        continue;
      }
      if (selection != INDEX_SHARD_WORK_OUTER) {
        index_shard_request_fatal_stop(shared);
        break;
      }

      result = &shared->results[index_order];

      if (index_shard_check_global_limits(shared) ||
          index_shard_master_stop_requested(shared)) {
        index_shard_master_limit_or_cancel_requested(
            shared,
            &result->hit_total_cpulimit,
            &result->hit_total_timelimit,
            &result->cancelled);
        index_shard_finish_outer_claim(ctx, shared, index_order);
        break;
      }

      /*
       * Prepare worker-local solver and field state only after this worker has
       * acquired real outer work. This preserves pass-local reuse without
       * charging field read/preprocessing costs to idle workers when the
       * configured index set is smaller than the pool.
       */
      if (!ctx->local_context_ready) {
        double context_wall_start = monotonic_seconds();

        if (index_shard_worker_prepare_pass(ctx, shared)) {
          ctx->pass_prepare_seconds =
              monotonic_seconds() - context_wall_start;

          index_shard_result_init(result, index_order);
          index_shard_result_fail(
              result,
              INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
              -1);

          index_shard_finish_outer_claim(ctx, shared, index_order);
          break;
        }

        ctx->pass_prepare_seconds =
            monotonic_seconds() - context_wall_start;
      }

      if (index_shard_run_one_with_worker_context(ctx,
                                                  shared,
                                                  index_order,
                                                  result,
                                                  mmap_advice)) {
        index_shard_finish_outer_claim(ctx, shared, index_order);
        if (result->failure_class == INDEX_SHARD_FAILURE_TASK_LOCAL) {
          continue;
        }
        break;
      }

      if (result->solved) {
        logverb("[index-shard] verified-result-ready worker=%i index_order=%zu "
                "best_logodds=%.3f field=%i wall=%.6f cpu=%.6f\n",
                ctx->worker_id,
                index_order,
                result->best_logodds,
                result->best_fieldnum,
                result->wall_seconds,
                (double)result->cpu_seconds);
      }

      if (index_shard_trace_enabled()) {
        logmsg("[index-shard] complete worker=%i index_order=%zu solved=%i "
               "failed=%i wall=%.6f cpu=%.6f pass_wall=%.6f\n",
               ctx->worker_id,
               index_order,
               result->solved,
               result->failed,
               result->wall_seconds,
               (double)result->cpu_seconds,
               monotonic_seconds() - shared->pass_wall_start);
      }

      index_shard_finish_outer_claim(ctx, shared, index_order);
      index_shard_check_global_limits(shared);
    }

    {
      double context_wall_start = monotonic_seconds();

      index_shard_worker_cleanup_pass(ctx, shared);
      ctx->pass_cleanup_seconds =
          monotonic_seconds() - context_wall_start;
    }

    index_shard_worker_done(shared);
  }

  fitsbin_payload_io_clear_thread_wait_helper();
  tls_status = index_shard_set_tls(NULL);
  if (tls_status) {
    logerr("[index-shard] failed to clear worker TLS "
           "worker=%i status=%i\n",
           ctx->worker_id,
           tls_status);
  }

  return NULL;
}
