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
// ANCHOR INDEX-SHARD: pool-submit
/*
 * Submit one onefield_run() pass to the persistent pool.
 *
 * This resets shared pass state, publishes result arrays, then increments
 * generation to wake workers.
 */
static int index_shard_pool_submit(
    index_shard_pool_t *pool,
    onefield_t *bp,
    solver_t *base_sp,
    size_t nindexes,
    const index_shard_hooks_t *hooks,
    const void *worker_view,
    index_shard_result_t *results,
    unsigned char *completed,
    unsigned char *outer_states) {
  index_shard_thread_state_t *shared = &pool->shared;
  int worker_count = pool->worker_count;
  int i;
  int payload_io_width;
  fitsbin_mmap_advice_t pass_mmap_advice;
  index_shard_width_plan_t pass_width_plan;
  anbool pass_exact_demand;

  if (!outer_states || !worker_view || !pool->producer_width) {
    return -1;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (!pool->pass_active || pool->shutdown ||
      pool->owner_bp != bp || pool->owner_sp != base_sp) {
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  payload_io_width = fitsbin_payload_io_service_width();
  pass_mmap_advice = fitsbin_mmap_advice_state_begin_pass(
      &base_sp->index_mmap_policy);
  pass_exact_demand = index_shard_config_exact_demand_pass(
      pool->payload_completion_registered,
      payload_io_width,
      fitsbin_payload_io_mapped_population_supported(),
      pass_mmap_advice == FITSBIN_MMAP_ADVICE_RANDOM,
      bp->indexnames ? (size_t)sl_size(bp->indexnames) : 0U,
      bp->indexes ? (size_t)pl_size(bp->indexes) : 0U,
      index_residency_service_active());
  if (index_shard_config_plan_widths(
          worker_count,
          payload_io_width,
          pool->payload_completion_registered,
          pass_exact_demand,
          &pass_width_plan)) {
    logerr("[index-shard] invalid pass width plan "
           "workers=%i payload_io_width=%i detached=%i "
           "exact_demand=%i\n",
           worker_count,
           payload_io_width,
           pool->payload_completion_registered ? 1 : 0,
           pass_exact_demand ? 1 : 0);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  pthread_mutex_lock(&shared->queue_mutex);
  for (i = 0; i < worker_count; i++) {
    if (pool->contexts[i].published_helper_group ||
        pool->contexts[i].published_staged_group ||
        pool->contexts[i].owner_waiting ||
        pool->contexts[i].owner_wake_pending ||
        pool->contexts[i].staged_owner_callback_active ||
        pool->contexts[i].staged_owner_callback_group ||
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] inner state remained active "
             "before pass worker=%i\n", i);
      pthread_mutex_unlock(&shared->queue_mutex);
      pthread_mutex_unlock(&pool->control_mutex);
      return -1;
    }
  }
  if (shared->completion_active ||
      shared->staged_submit_callbacks_active) {
    logerr("[index-shard] completion registry remained before pass "
           "active=%zu submit_callbacks=%zu\n",
           shared->completion_active,
           shared->staged_submit_callbacks_active);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->queue_waiters) {
    logerr("[index-shard] queue waiters remained before pass "
           "count=%zu\n", shared->queue_waiters);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->helper_groups_published !=
          shared->helper_groups_completed ||
      shared->staged_groups_published !=
          shared->staged_groups_completed) {
    logerr("[index-shard] inner group lifecycle remained "
           "before pass helper=%llu/%llu staged=%llu/%llu\n",
           shared->helper_groups_published,
           shared->helper_groups_completed,
           shared->staged_groups_published,
           shared->staged_groups_completed);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->helper_groups_active ||
      shared->helper_preparations_active ||
      shared->helper_foreign_reservations ||
      shared->staged_groups_active ||
      shared->staged_tickets_active ||
      shared->staged_source_leases ||
      shared->staged_compute_ready ||
      shared->staged_reorder_ready ||
      shared->staged_compute_running_global ||
      shared->staged_submit_backpressure) {
    logerr("[index-shard] inner activity remained before pass "
           "helper_groups=%zu preparations=%zu reservations=%zu "
           "staged_groups=%zu tickets=%zu leases=%zu "
           "compute_ready=%zu reorder_ready=%zu "
           "compute_running=%zu submit_backpressure=%i\n",
           shared->helper_groups_active,
           shared->helper_preparations_active,
           shared->helper_foreign_reservations,
           shared->staged_groups_active,
           shared->staged_tickets_active,
           shared->staged_source_leases,
           shared->staged_compute_ready,
           shared->staged_reorder_ready,
           shared->staged_compute_running_global,
           shared->staged_submit_backpressure);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  pthread_mutex_lock(&shared->result_mutex);
  pthread_mutex_lock(&shared->state_mutex);
  pthread_mutex_lock(&shared->limit_mutex);

  shared->bp = bp;
  shared->base_sp = base_sp;
  shared->hooks = hooks;
  shared->worker_view = worker_view;
  shared->nindexes = nindexes;
  shared->canonical_scan_cursor = 0U;
  shared->outer_unclaimed = nindexes;
  shared->outer_running = 0U;
  shared->producer_width = MIN(
      pass_width_plan.producer_width, nindexes);
  /* Idle outer capacity joins the inner width for this pass. */
  shared->helper_width =
      (size_t)worker_count - shared->producer_width;
  shared->queue_waiters = 0U;
  shared->helper_groups_active = 0U;
  shared->helper_preparations_active = 0U;
  shared->helper_foreign_reservations = 0U;
  shared->staged_groups_active = 0U;
  shared->staged_tickets_active = 0U;
  shared->staged_source_leases = 0U;
  shared->staged_compute_ready = 0U;
  shared->staged_reorder_ready = 0U;
  shared->staged_compute_running_global = 0U;
  shared->staged_max_compute_running = 0U;
  shared->staged_max_compute_running_global = 0U;
  shared->staged_completion_epoch = 1U;
  shared->staged_submit_backpressure = FALSE;
  shared->staged_submit_callbacks_active = 0U;
  shared->completion_registry_error = FALSE;
  shared->outer_states = outer_states;
  shared->outer_claims = 0U;
  shared->helper_groups_published = 0U;
  shared->helper_groups_completed = 0U;
  shared->helper_tasks_owner = 0U;
  shared->helper_tasks_foreign = 0U;
  shared->helper_task_failures = 0U;
  shared->helper_owner_wait_calls = 0U;
  shared->helper_owner_wait_seconds = 0.0;
  shared->staged_groups_published = 0U;
  shared->staged_groups_completed = 0U;
  shared->staged_tasks_owner = 0U;
  shared->staged_tasks_foreign = 0U;
  shared->staged_compute_owner = 0U;
  shared->staged_compute_foreign = 0U;
  shared->staged_ready_before_outer_claims = 0U;
  shared->staged_task_failures = 0U;
  shared->staged_io_submitted = 0U;
  shared->staged_io_completed = 0U;
  shared->staged_submit_retries = 0U;
  shared->staged_submit_rearms = 0U;
  shared->staged_submit_handoffs = 0U;
  shared->staged_submit_deferrals = 0U;
  shared->staged_owner_wait_calls = 0U;
  shared->staged_owner_wait_seconds = 0.0;
  shared->staged_max_io_submitted = 0U;
  shared->staged_max_compute_ready = 0U;
  shared->staged_max_reorder_ready = 0U;
  shared->staged_prepare_claims = 0U;
  shared->staged_submit_claims = 0U;
  shared->staged_poll_claims = 0U;
  shared->staged_inline_poll_claims = 0U;
  shared->staged_inline_poll_failures = 0U;
  shared->staged_execute_claims = 0U;
  shared->staged_owner_execute_claims = 0U;
  shared->staged_submit_to_ready_seconds = 0.0;
  shared->staged_ready_dwell_seconds = 0.0;
  shared->staged_execute_seconds = 0.0;
  shared->staged_result_to_retire_seconds = 0.0;
  shared->staged_retire_seconds = 0.0;
  shared->task_local_failures = 0U;
  shared->global_integrity_failures = 0U;
  shared->late_loser_failures = 0U;
  shared->observability_enabled =
      log_get_level() >= LOG_VERB;
  shared->queue_broadcasts = 0U;
  shared->queue_signals = 0U;
  shared->queue_signals_no_work = 0U;
  shared->owner_signals = 0U;
  shared->owner_signals_coalesced = 0U;
  shared->owner_broadcasts = 0U;
  shared->queue_wait_calls = 0U;
  shared->queue_wait_seconds = 0.0;
  shared->completion_notifications = 0U;
  shared->completion_groups_scanned = 0U;
  shared->completion_tasks_scanned = 0U;
  shared->completion_matches = 0U;
  shared->completion_registry_registered = 0U;
  shared->completion_registry_removed = 0U;
  shared->completion_registry_early = 0U;
  shared->completion_registry_misses = 0U;
  shared->completion_registry_duplicates = 0U;
  shared->completion_registry_invalid = 0U;
  shared->selection_scans = 0U;
  shared->selection_groups_scanned = 0U;
  shared->selection_tasks_scanned = 0U;
  shared->selection_misses = 0U;
  shared->helper_windows_started = 0U;
  shared->helper_window_seconds = 0.0;
  shared->verification_helper_groups = 0U;
  shared->verification_helper_contexts = 0U;

  shared->results = results;
  shared->completed = completed;
  shared->results_reduced = 0U;
  shared->next_completion_sequence = 0U;
  shared->next_candidate_sequence = 0U;

  shared->worker_count = worker_count;
  shared->active_workers = worker_count;

  shared->stop_requested = FALSE;
  shared->fatal_error = FALSE;
  shared->winner_selected = FALSE;
  shared->solved_published = FALSE;
  shared->master_committed = FALSE;
  shared->terminal_cause = INDEX_SHARD_TERMINAL_NONE;
  shared->first_stop_wall_since_pass = -1.0;
  __atomic_store_n(
      &shared->worker_stop_requested,
      FALSE,
      __ATOMIC_RELEASE);
  shared->selected_index_order = nindexes;
  shared->selected_candidate_sequence = 0U;
  shared->have_committed_result = FALSE;
  shared->committed_index_order = 0U;
  shared->limit_reported = FALSE;

  shared->reducer_work_calls = 0U;
  shared->reducer_work_wall_seconds = 0.0;

  /*
   * Preserve the established mmap policy. Coarse CodeKD packets may populate
   * a bounded exact DATA/PERM plan, but there is no independent delivery
   * generation and native mmap demand remains authoritative.
   */
  shared->mmap_advice = pass_mmap_advice;
  shared->mmap_pass_number =
      base_sp->index_mmap_policy.pass_number;
  shared->mmap_advice_failures = 0U;

  // Pass timing excludes final solve-field output generation.
  shared->pass_wall_start = monotonic_seconds();
  shared->pass_cpu_start = get_cpu_usage();
  shared->pass_rusage_valid =
      (getrusage(RUSAGE_SELF,
                 &shared->pass_rusage_start) == 0);

  pthread_mutex_unlock(&shared->limit_mutex);
  pthread_mutex_unlock(&shared->state_mutex);
  pthread_mutex_unlock(&shared->result_mutex);
  pthread_mutex_unlock(&shared->queue_mutex);

  /*
   * The prior generation is quiescent and the next one is not visible yet.
   * Match explicit demand admission to its cold-owner width without holding
   * a shard mutex across the payload provider's own lock.
   */
  fitsbin_payload_io_configure_workers(
      pass_exact_demand
          ? (int)pass_width_plan.producer_width
          : worker_count);

  // Generation publication is the hard band barrier release.
  pool->generation++;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);

  logverb("[index-shard] pthread-pool submit compute_width=%i "
          "producer_width=%zu helper_width=%zu "
          "candidates=%zu engine_pass=%zu "
          "depth_index=%zu scale_index=%zu startobj=%i endobj=%i "
          "scheduler=first-valid-completion inner_scheduler=ordered-codekd-packets "
          "mmap_pass=%u mmap_advice=%s "
          "mmap_scope=payload-random-topology-normal "
          "mmap_policy=parallel-random-serial-normal "
          "page_delivery=%s payload_io=%s payload_io_width=%i "
          "outer_admission=%s\n",
          worker_count,
          shared->producer_width,
          shared->helper_width,
          nindexes,
          bp->engine_pass_ordinal,
          bp->engine_depth_index,
          bp->engine_scale_index,
          base_sp->startobj,
          base_sp->endobj,
          shared->mmap_pass_number,
          fitsbin_mmap_advice_name(shared->mmap_advice),
          pool->payload_completion_registered
              ? "detached-bounded-mapped-completion"
              : "native-mmap-demand",
          payload_io_width > 0 ? "kernel-page-cache" : "native-mmap",
          payload_io_width,
          !pool->payload_completion_registered
              ? "reserved-helper"
              : (pass_exact_demand
                     ? "provider-bounded-exact-demand"
                     : "full-producer"));

  return 0;
}

/*
 * SECTION INDEX-SHARD: entry
 */

// ANCHOR INDEX-SHARD: entry
/*
 * Execute one complete index-shard pass.
 *
 * Terminal status is classified by failure scope and the master mutation
 * boundary. Only an unavailable path or an isolated task-local failure proven
 * to occur before master commit may return control to the serial path. A
 * global-integrity failure is terminal regardless of winner timing.
 */
index_shard_solve_status_t
index_shard_solve_impl(onefield_t *bp,
                       solver_t *base_sp,
                       size_t nindexes,
                       const index_shard_hooks_t *hooks) {
  index_shard_pool_t *pool;
  unsigned char *outer_states = NULL;
  index_shard_result_t *results = NULL;
  unsigned char *completed = NULL;
  void *worker_view = NULL;
  size_t i;
  int acquire_rc;
  int rc = 0;
  index_shard_pass_state_snapshot_t state;
  index_shard_pass_metrics_snapshot_t pass_metrics;
  index_shard_task_profile_snapshot_t task_profile;
  index_shard_phase_profile_snapshot_t phase_profile;
  solver_profile_t solver_profile;
  index_shard_solve_status_t status = INDEX_SHARD_SOLVE_HANDLED;
  double context_cleanup_max_seconds = 0.0;
  double context_cleanup_wall_seconds = 0.0;
  double context_prepare_max_seconds = 0.0;
  double context_prepare_wall_seconds = 0.0;
  double reduction_ex_verify_seconds;
  double stop_to_quiescence_seconds = 0.0;
  anbool pass_completed;
  anbool pass_exhaustive;
  anbool pass_solved;
  anbool pass_cancelled;
  anbool clean_exhaustion_required;
  anbool helper_quiescence_valid = TRUE;
  anbool mmap_transitioned;
  anbool worker_cancelled = FALSE;

  if (!index_shard_pthread_enabled(bp)) {
    return INDEX_SHARD_SOLVE_UNAVAILABLE;
  }

  // no candidate indexes, nothing to do
  if (!nindexes) {
    return INDEX_SHARD_SOLVE_HANDLED;
  }

  if (!hooks) {
    ERROR("index-shard hooks are NULL");
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  acquire_rc = index_shard_pool_acquire_pass(bp, base_sp, &pool);

  if (acquire_rc == INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE) {
    logverb("[index-shard] pthread mode requested but pool inactive\n");
    return INDEX_SHARD_SOLVE_UNAVAILABLE;
  }

  if (acquire_rc != INDEX_SHARD_POOL_ACQUIRE_OK) {
    logerr("[index-shard] pool ownership or pass lifecycle conflict; "
           "serial fallback suppressed\n");
    return INDEX_SHARD_SOLVE_LIFECYCLE_CONFLICT;
  }

  if (!hooks->create_worker_view ||
      !hooks->destroy_worker_view ||
      !hooks->prepare_local_context ||
      hooks->create_worker_view(
          bp, base_sp, &worker_view) ||
      !worker_view) {
    logerr("[index-shard] failed to create immutable worker view\n");
    if (worker_view && hooks->destroy_worker_view) {
      hooks->destroy_worker_view(worker_view);
    }
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  results = calloc(nindexes, sizeof(index_shard_result_t));
  completed = calloc(nindexes, sizeof(unsigned char));
  outer_states = calloc(nindexes, sizeof(*outer_states));

  if (!results || !completed || !outer_states) {
    SYSERROR("Failed to allocate index-shard pass state");

    free(outer_states);
    free(results);
    free(completed);

    hooks->destroy_worker_view(worker_view);
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  // Submit releases the hard current-band barrier to persistent workers.
  rc = index_shard_pool_submit(
      pool,
      bp,
      base_sp,
      nindexes,
      hooks,
      worker_view,
      results,
      completed,
      outer_states);

  if (rc) {
    free(outer_states);
    free(results);
    free(completed);

    hooks->destroy_worker_view(worker_view);
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  rc = index_shard_pool_reduce_first_valid(pool);

  index_shard_pass_state_snapshot(&pool->shared, &state);
  clean_exhaustion_required =
      !state.solved_published &&
      !state.stop_requested &&
      !state.fatal_error;
  pthread_mutex_lock(&pool->shared.queue_mutex);
  for (i = 0U; i < (size_t)pool->shared.worker_count; i++) {
    if (pool->contexts[i].published_helper_group ||
        pool->contexts[i].published_staged_group ||
        pool->contexts[i].staged_owner_callback_active ||
        pool->contexts[i].staged_owner_callback_group ||
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] inner state remained active "
             "after worker quiescence worker=%zu\n", i);
      rc = -1;
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
      helper_quiescence_valid = FALSE;
    }
  }
  if (pool->shared.queue_waiters) {
    logerr("[index-shard] queue waiters remained after "
           "worker quiescence count=%zu\n",
           pool->shared.queue_waiters);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_groups_active) {
    logerr("[index-shard] helper groups remained active after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_groups_active);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_preparations_active) {
    logerr("[index-shard] helper preparations remained active after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_preparations_active);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_foreign_reservations) {
    logerr("[index-shard] helper reservations remained after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_foreign_reservations);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_groups_published !=
      pool->shared.helper_groups_completed) {
    logerr("[index-shard] helper group lifecycle mismatch "
           "after worker quiescence published=%llu completed=%llu\n",
           pool->shared.helper_groups_published,
           pool->shared.helper_groups_completed);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.staged_groups_active ||
      pool->shared.staged_tickets_active ||
      pool->shared.staged_source_leases ||
      pool->shared.staged_compute_ready ||
      pool->shared.staged_reorder_ready ||
      pool->shared.staged_compute_running_global) {
    logerr("[index-shard] staged activity remained after "
           "worker quiescence groups=%zu tickets=%zu leases=%zu "
           "compute_ready=%zu reorder_ready=%zu "
           "compute_running=%zu\n",
           pool->shared.staged_groups_active,
           pool->shared.staged_tickets_active,
           pool->shared.staged_source_leases,
           pool->shared.staged_compute_ready,
           pool->shared.staged_reorder_ready,
           pool->shared.staged_compute_running_global);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.completion_active ||
      pool->shared.staged_submit_callbacks_active ||
      pool->shared.completion_registry_error) {
    logerr("[index-shard] completion registry remained after "
           "worker quiescence active=%zu submit_callbacks=%zu "
           "error=%i\n",
           pool->shared.completion_active,
           pool->shared.staged_submit_callbacks_active,
           pool->shared.completion_registry_error ? 1 : 0);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.staged_groups_published !=
      pool->shared.staged_groups_completed) {
    logerr("[index-shard] staged group lifecycle mismatch "
           "after worker quiescence published=%llu completed=%llu\n",
           pool->shared.staged_groups_published,
           pool->shared.staged_groups_completed);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.outer_running) {
    logerr("[index-shard] pass ended with %zu outer owners\n",
           pool->shared.outer_running);
    rc = -1;
  }
  if (clean_exhaustion_required &&
      (pool->shared.outer_unclaimed ||
       pool->shared.canonical_scan_cursor != nindexes)) {
    logerr("[index-shard] clean pass ended before claim exhaustion "
           "unclaimed=%zu cursor=%zu candidates=%zu\n",
           pool->shared.outer_unclaimed,
           pool->shared.canonical_scan_cursor,
           nindexes);
    rc = -1;
  }
  for (i = 0U; i < nindexes; i++) {
    if (outer_states[i] == INDEX_SHARD_OUTER_RUNNING) {
      logerr("[index-shard] pass ended with running index_order=%zu\n",
             i);
      rc = -1;
      break;
    }
    if (clean_exhaustion_required &&
        outer_states[i] != INDEX_SHARD_OUTER_FINISHED) {
      logerr("[index-shard] clean pass ended with unfinished "
             "index_order=%zu state=%u\n",
             i,
             (unsigned int)outer_states[i]);
      rc = -1;
      break;
    }
  }
  pthread_mutex_unlock(&pool->shared.queue_mutex);

  logverb("[index-shard] ownership-pass generation=%lu "
          "scale_index=%zu canonical_claims=%llu "
          "producer_width=%zu helper_width=%zu "
          "unclaimed=%zu quiescent=%i\n",
          pool->generation,
          bp->engine_scale_index,
          pool->shared.outer_claims,
          pool->shared.producer_width,
          pool->shared.helper_width,
          pool->shared.outer_unclaimed,
          helper_quiescence_valid ? 1 : 0);

  memset(&solver_profile, 0, sizeof(solver_profile));

  /*
   * Worker-local cancellation is published only after every result slot is
   * immutable. Master onefield state remains reducer/caller-owned, while the
   * limit mutex preserves the existing synchronization discipline for the
   * cancellation flag.
   */
  for (i = 0; i < nindexes; i++) {
    if (results[i].task_started) {
      solver_profile_accumulate(
          &solver_profile,
          &results[i].solver_profile);
    }

    if (results[i].cancelled) {
      worker_cancelled = TRUE;
    }
  }

  if (worker_cancelled && !state.winner_selected) {
    pthread_mutex_lock(&pool->shared.limit_mutex);
    bp->cancelled = TRUE;
    pthread_mutex_unlock(&pool->shared.limit_mutex);
  }

  /*
   * active_workers reached zero before the first-valid reducer returned, so every
   * participating context timing is immutable here. These are work sums; the
   * maxima expose the critical per-worker prepare/cleanup contribution.
   */
  for (i = 0; i < (size_t)pool->shared.worker_count; i++) {
    double cleanup_seconds = pool->contexts[i].pass_cleanup_seconds;
    double prepare_seconds = pool->contexts[i].pass_prepare_seconds;

    context_cleanup_wall_seconds += cleanup_seconds;
    context_prepare_wall_seconds += prepare_seconds;

    if (cleanup_seconds > context_cleanup_max_seconds) {
      context_cleanup_max_seconds = cleanup_seconds;
    }

    if (prepare_seconds > context_prepare_max_seconds) {
      context_prepare_max_seconds = prepare_seconds;
    }
  }

  index_shard_pass_state_snapshot(&pool->shared, &state);

  reduction_ex_verify_seconds =
      solver_profile.resolve_wall_seconds -
      solver_profile.verify_wall_seconds;

  if (reduction_ex_verify_seconds < 0.0) {
    reduction_ex_verify_seconds = 0.0;
  }

  if (state.fatal_error || state.global_integrity_failures) {
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
  } else if (rc &&
             status != INDEX_SHARD_SOLVE_TERMINAL_FAILURE) {
    if (state.master_committed) {
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    } else {
      status = INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
    }
  }

   /*
   * Defensive ownership invariant: once the reducer has crossed the master
   * mutation boundary, a precommit failure classification is impossible.
   */
  if (state.master_committed &&
      status == INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE) {
    logerr("[index-shard] invalid precommit status after master commit; "
           "promoting to terminal failure\n");
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
  }

   /*
   * Workers have left the pass and result slots are now immutable. Snapshot
   * pass timing and task-duration distribution before destroying results.
   */
  index_shard_pass_metrics_snapshot(&pool->shared,
                                    &pass_metrics);

  if (state.first_stop_wall_since_pass >= 0.0 &&
      pass_metrics.wall_seconds >= state.first_stop_wall_since_pass) {
    stop_to_quiescence_seconds =
        pass_metrics.wall_seconds - state.first_stop_wall_since_pass;
  }

  index_shard_task_profile_snapshot(results,
                                    nindexes,
                                    pass_metrics.wall_seconds,
                                    &task_profile);

  index_shard_phase_profile_snapshot(results,
                                     nindexes,
                                     &phase_profile);

  /*
   * index_shard_pool_reduce_first_valid() returned and every participating
   * worker has left this generation. The pass outcome is now immutable.
   */
  pass_completed =
      task_profile.executed == nindexes;

  pass_exhaustive =
      pass_completed &&
      pass_metrics.reduced == nindexes;

  pass_solved =
      state.solved_published ||
      bp->single_field_solved;

  pass_cancelled =
      !pass_solved &&
      (bp->cancelled ||
       bp->hit_total_cpulimit ||
       bp->hit_total_timelimit ||
       state.stop_requested);

  if (status == INDEX_SHARD_SOLVE_HANDLED &&
      rc == 0 &&
      state.solved_published) {
    if (index_shard_report_committed_solution(
            bp,
            nindexes,
            &pool->shared,
            hooks,
            results)) {
      rc = -1;
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    }
  }

  mmap_transitioned =
      fitsbin_mmap_policy_complete_pass(
          &base_sp->index_mmap_policy,
          pass_completed,
          pass_exhaustive,
          pass_solved,
          pass_cancelled,
          rc,
          (int)status);
  // dispose unmerged worker results after all workers have left pass
  for (i = 0; i < nindexes; i++) {
    index_shard_result_dispose(&results[i], hooks);
  }

  free(results);
  free(completed);

  logverb("[index-shard] done workers=%i solved=%i "
          "wall=%.6f cpu=%.6f utilization=%.1f%% "
          "effective_concurrency=%.2f stop_to_quiescence=%.6f\n",
          pool->worker_count,
          bp->single_field_solved,
          pass_metrics.wall_seconds,
          (double)pass_metrics.cpu_seconds,
          pass_metrics.cpu_percent,
          pass_metrics.cpu_percent / 100.0,
          stop_to_quiescence_seconds);

  logverb("[index-shard] reducer-pass generation=%lu candidates=%zu "
          "calls=%llu work_wall_sum=%.6f\n",
          pool->generation,
          nindexes,
          pool->shared.reducer_work_calls,
          pool->shared.reducer_work_wall_seconds);

  logverb("[index-shard] helper-pass generation=%lu "
          "groups=%llu completed=%llu owner_tasks=%llu "
          "foreign_tasks=%llu task_failures=%llu "
          "owner_waits=%llu owner_wait_seconds=%.6f\n",
          pool->generation,
          pool->shared.helper_groups_published,
          pool->shared.helper_groups_completed,
          pool->shared.helper_tasks_owner,
          pool->shared.helper_tasks_foreign,
          pool->shared.helper_task_failures,
          pool->shared.helper_owner_wait_calls,
          pool->shared.helper_owner_wait_seconds);

  logverb("[index-shard] staged-pass generation=%lu "
          "groups=%llu completed=%llu owner_claims=%llu "
          "foreign_claims=%llu owner_compute=%llu "
          "foreign_compute=%llu ready_before_outer=%llu "
          "task_failures=%llu "
          "io_submitted=%llu io_completed=%llu submit_retries=%llu "
          "submit_rearms=%llu submit_handoffs=%llu "
          "submit_deferrals=%llu "
          "max_io_submitted=%zu max_compute_ready=%zu "
          "max_reorder_ready=%zu max_compute_running=%zu "
          "max_compute_running_global=%zu prepare_claims=%llu "
          "submit_claims=%llu poll_claims=%llu inline_polls=%llu "
          "inline_poll_failures=%llu execute_claims=%llu "
          "owner_execute_claims=%llu submit_to_ready_sum=%.6f "
          "ready_dwell_sum=%.6f execute_sum=%.6f "
          "result_to_retire_sum=%.6f retire_sum=%.6f "
          "owner_waits=%llu owner_wait_seconds=%.6f\n",
          pool->generation,
          pool->shared.staged_groups_published,
          pool->shared.staged_groups_completed,
          pool->shared.staged_tasks_owner,
          pool->shared.staged_tasks_foreign,
          pool->shared.staged_compute_owner,
          pool->shared.staged_compute_foreign,
          pool->shared.staged_ready_before_outer_claims,
          pool->shared.staged_task_failures,
          pool->shared.staged_io_submitted,
          pool->shared.staged_io_completed,
          pool->shared.staged_submit_retries,
          pool->shared.staged_submit_rearms,
          pool->shared.staged_submit_handoffs,
          pool->shared.staged_submit_deferrals,
          pool->shared.staged_max_io_submitted,
          pool->shared.staged_max_compute_ready,
          pool->shared.staged_max_reorder_ready,
          pool->shared.staged_max_compute_running,
          pool->shared.staged_max_compute_running_global,
          pool->shared.staged_prepare_claims,
          pool->shared.staged_submit_claims,
          pool->shared.staged_poll_claims,
          pool->shared.staged_inline_poll_claims,
          pool->shared.staged_inline_poll_failures,
          pool->shared.staged_execute_claims,
          pool->shared.staged_owner_execute_claims,
          pool->shared.staged_submit_to_ready_seconds,
          pool->shared.staged_ready_dwell_seconds,
          pool->shared.staged_execute_seconds,
          pool->shared.staged_result_to_retire_seconds,
          pool->shared.staged_retire_seconds,
          pool->shared.staged_owner_wait_calls,
          pool->shared.staged_owner_wait_seconds);

  logverb("[index-shard] scheduler-observability generation=%lu "
          "queue_broadcasts=%llu queue_signals=%llu "
          "queue_no_work=%llu owner_signals=%llu "
          "owner_coalesced=%llu owner_broadcasts=%llu "
          "queue_waits=%llu queue_wait_seconds=%.6f completions=%llu "
          "completion_groups=%llu completion_tasks=%llu "
          "completion_matches=%llu completion_registered=%llu "
          "completion_removed=%llu completion_early=%llu "
          "completion_misses=%llu completion_duplicates=%llu "
          "completion_invalid=%llu completion_active=%zu "
          "selection_scans=%llu "
          "selection_groups=%llu selection_tasks=%llu "
          "selection_misses=%llu helper_windows=%llu "
          "helper_window_seconds=%.6f verification_groups=%llu "
          "verification_contexts=%llu "
          "max_compute_running_global=%zu\n",
          pool->generation,
          pool->shared.queue_broadcasts,
          pool->shared.queue_signals,
          pool->shared.queue_signals_no_work,
          pool->shared.owner_signals,
          pool->shared.owner_signals_coalesced,
          pool->shared.owner_broadcasts,
          pool->shared.queue_wait_calls,
          pool->shared.queue_wait_seconds,
          pool->shared.completion_notifications,
          pool->shared.completion_groups_scanned,
          pool->shared.completion_tasks_scanned,
          pool->shared.completion_matches,
          pool->shared.completion_registry_registered,
          pool->shared.completion_registry_removed,
          pool->shared.completion_registry_early,
          pool->shared.completion_registry_misses,
          pool->shared.completion_registry_duplicates,
          pool->shared.completion_registry_invalid,
          pool->shared.completion_active,
          pool->shared.selection_scans,
          pool->shared.selection_groups_scanned,
          pool->shared.selection_tasks_scanned,
          pool->shared.selection_misses,
          pool->shared.helper_windows_started,
          pool->shared.helper_window_seconds,
          pool->shared.verification_helper_groups,
          pool->shared.verification_helper_contexts,
          pool->shared.staged_max_compute_running_global);

  logverb("[index-shard] context-pass generation=%lu candidates=%zu "
          "compute_width=%i producer_width=%zu helper_width=%zu "
          "prepare_work_wall_sum=%.6f prepare_max=%.6f "
          "cleanup_work_wall_sum=%.6f cleanup_max=%.6f\n",
          pool->generation,
          nindexes,
          pool->shared.worker_count,
          pool->shared.producer_width,
          pool->shared.helper_width,
          context_prepare_wall_seconds,
          context_prepare_max_seconds,
          context_cleanup_wall_seconds,
          context_cleanup_max_seconds);

  logverb("[index-shard] solver-pass generation=%lu candidates=%zu "
          "detailed=%i failed=%i solver_run_work_wall_sum=%.6f "
          "codekd_work_wall_sum=%.6f codekd_calls=%llu "
          "codekd_hits=%llu "
          "resolve_work_wall_sum=%.6f "
          "reduction_ex_verify_hit_work_wall_sum=%.6f "
          "resolve_calls=%llu verify_hit_work_wall_sum=%.6f "
          "verify_calls=%llu "
          "batches=%llu completed=%llu stopped=%llu "
          "batch_failed=%llu hypotheses=%llu executed=%llu reduced=%llu "
          "task_ranges=%llu tasks_executed=%llu submitted=%llu "
          "inline=%llu parallel_batches=%llu observed_parallel=%llu "
          "parallel_hypotheses=%llu alloc_failures=%llu "
          "search_failures=%llu max_batch=%zu max_tasks=%zu "
          "max_parallel=%zu helper_tasks=%llu "
          "helper_combinations=%llu "
          "hypothesis_order=%016llx kd_result_order=%016llx "
          "candidate_order=%016llx\n",
          pool->generation,
          nindexes,
          solver_profile.detailed ? 1 : 0,
          solver_profile.execution_failed ? 1 : 0,
          solver_profile.solver_run_wall_seconds,
          solver_profile.codekd_wall_seconds,
          solver_profile.codekd_calls,
          solver_profile.codekd_hits,
          solver_profile.resolve_wall_seconds,
          reduction_ex_verify_seconds,
          solver_profile.resolve_calls,
          solver_profile.verify_wall_seconds,
          solver_profile.verify_calls,
          solver_profile.hypothesis_batches,
          solver_profile.hypothesis_batches_completed,
          solver_profile.hypothesis_batches_stopped,
          solver_profile.hypothesis_batches_failed,
          solver_profile.hypotheses_generated,
          solver_profile.hypotheses_executed,
          solver_profile.hypotheses_reduced,
          solver_profile.task_ranges_planned,
          solver_profile.task_ranges_executed,
          solver_profile.task_ranges_submitted,
          solver_profile.task_ranges_inline,
          solver_profile.parallel_batches,
          solver_profile.parallel_batches_observed,
          solver_profile.parallel_hypotheses,
          solver_profile.allocation_failures,
          solver_profile.search_failures,
          solver_profile.max_batch_hypotheses,
          solver_profile.max_task_ranges,
          solver_profile.max_parallel_ranges,
          solver_profile.ab_helper_tasks,
          solver_profile.ab_helper_combinations,
          solver_profile.hypothesis_order_hash,
          solver_profile.kd_result_order_hash,
          solver_profile.candidate_order_hash);

  logverb("[index-shard] page-pipeline generation=%lu "
          "descriptors=%llu complete=%llu boundary_deferrals=%llu "
          "raw_hints=%llu unique_pages=%llu coalesced_ranges=%llu "
          "logical_bytes=%llu aligned_bytes=%llu overread_bytes=%llu "
          "refusals=%llu/%llu/%llu/%llu/%llu/%llu/%llu/%llu/%llu "
          "candidate_delivery=%llu quad=%llu/%llu/%llu "
          "star=%llu/%llu/%llu windows=%llu "
          "rows=%llu/%llu/%llu/%llu "
          "verify_pages=%llu/%llu prefixes=%llu "
          "tickets=%llu/%llu fallback=%llu ready_rows=%llu "
          "ranges=%llu bytes=%llu/%llu candidate_math=%llu/%llu "
          "verify_score=%llu/%llu/%llu/%llu/%llu/%llu "
          "work=%llu work_wall_sum=%.6f\n",
          pool->generation,
          solver_profile.page_plan_descriptors_total,
          solver_profile.page_plan_descriptors_complete,
          solver_profile.page_plan_descriptor_splits,
          solver_profile.page_plan_raw_ranges,
          solver_profile.page_plan_unique_pages,
          solver_profile.page_plan_ranges_after_dedup,
          solver_profile.page_plan_logical_bytes,
          solver_profile.page_plan_aligned_bytes,
          solver_profile.page_plan_overread_bytes,
          solver_profile.page_plan_not_applicable,
          solver_profile.page_plan_allocation_refused,
          solver_profile.page_plan_source_mismatch,
          solver_profile.page_plan_invalid_range,
          solver_profile.page_plan_byte_budget_refused,
          solver_profile.page_plan_range_capacity_refused,
          solver_profile.page_plan_service_refused,
          solver_profile.page_plan_service_errors,
          solver_profile.page_plan_cancelled,
          solver_profile.candidate_delivery_candidates,
          solver_profile.candidate_quad_submitted,
          solver_profile.candidate_quad_ready,
          solver_profile.candidate_quad_fallback,
          solver_profile.candidate_star_submitted,
          solver_profile.candidate_star_ready,
          solver_profile.candidate_star_fallback,
          solver_profile.candidate_delivery_windows,
          solver_profile.candidate_quad_ready_rows,
          solver_profile.candidate_star_ready_rows,
          solver_profile.candidate_retired_rows,
          solver_profile.candidate_native_rows,
          solver_profile.verification_page_queries,
          solver_profile.verification_page_queries_planned,
          solver_profile.verification_page_prefixes,
          solver_profile.verification_page_submitted,
          solver_profile.verification_page_ready,
          solver_profile.verification_page_fallback,
          solver_profile.verification_page_ready_rows,
          solver_profile.verification_page_ranges,
          solver_profile.verification_page_logical_bytes,
          solver_profile.verification_page_aligned_bytes,
          solver_profile.candidate_math_prepared,
          solver_profile.candidate_math_reused,
          solver_profile.verification_score_batches_prepared,
          solver_profile.verification_score_contexts_prepared,
          solver_profile.verification_score_batches_executed,
          solver_profile.verification_score_contexts_completed,
          solver_profile.verification_score_fallback_batches,
          solver_profile.verification_score_stopped_batches,
          solver_profile.verification_score_work_units_completed,
          solver_profile.verification_score_wall_seconds);

  logverb("[index-shard] pass-detail candidates=%zu reduced=%zu "
          "hit_total_cpu_limit=%i hit_total_wall_limit=%i "
          "cancelled=%i rc=%i status=%i "
          "master_committed=%i winner_selected=%i terminal=%s "
          "selected_order=%zu selected_sequence=%zu "
          "task_local_failures=%llu global_integrity_failures=%llu "
          "late_loser_failures=%llu\n",
          nindexes,
          pass_metrics.reduced,
          bp->hit_total_cpulimit,
          bp->hit_total_timelimit,
          bp->cancelled,
          rc,
          (int)status,
          state.master_committed,
          state.winner_selected,
          index_shard_terminal_cause_name(state.terminal_cause),
          state.selected_index_order,
          state.selected_candidate_sequence,
          state.task_local_failures,
          state.global_integrity_failures,
          pool->shared.late_loser_failures);

  logverb("[index-shard] mmap-policy "
         "policy=%s effective=%s scope=all-chunks pass=%u "
         "clean_unsolved_passes=%u transitions=%u "
         "transitioned=%i completed=%i exhaustive=%i "
         "solved=%i cancelled=%i advice_failures=%llu\n",
         fitsbin_mmap_policy_name(
             base_sp->index_mmap_policy.policy),
         fitsbin_mmap_advice_name(
             base_sp->index_mmap_policy.effective_advice),
         base_sp->index_mmap_policy.pass_number,
         base_sp->index_mmap_policy
             .completed_clean_unsolved_passes,
         base_sp->index_mmap_policy.transition_count,
         mmap_transitioned ? 1 : 0,
         pass_completed ? 1 : 0,
         pass_exhaustive ? 1 : 0,
         pass_solved ? 1 : 0,
         pass_cancelled ? 1 : 0,
         pool->shared.mmap_advice_failures);

   if (!task_profile.executed) {
    logverb("[index-shard] task-profile executed=0\n");
  } else if (task_profile.quantiles_available) {
    logverb("[index-shard] task-profile executed=%zu "
           "task_p50=%.6f task_p90=%.6f task_p99=%.6f "
           "task_max=%.6f max_order=%zu max_worker=%i "
           "max_solve=%.6f skew_max_p50=%.1f "
           "max_pool_pct=%.1f%% serial_tail=%.6f "
           "serial_tail_pct=%.1f%% tail_order=%zu tail_worker=%i\n",
           task_profile.executed,
           task_profile.task_p50_seconds,
           task_profile.task_p90_seconds,
           task_profile.task_p99_seconds,
           task_profile.task_max_seconds,
           task_profile.max_index_order,
           task_profile.max_worker_id,
           task_profile.max_solve_seconds,
           task_profile.max_to_p50,
           task_profile.max_pool_percent,
           task_profile.serial_tail_seconds,
           task_profile.serial_tail_percent,
           task_profile.tail_index_order,
           task_profile.tail_worker_id);
  } else {
    logverb("[index-shard] task-profile executed=%zu "
           "quantiles=unavailable task_max=%.6f "
           "max_order=%zu max_worker=%i max_solve=%.6f "
           "max_pool_pct=%.1f%% serial_tail=%.6f "
           "serial_tail_pct=%.1f%% tail_order=%zu tail_worker=%i\n",
           task_profile.executed,
           task_profile.task_max_seconds,
           task_profile.max_index_order,
           task_profile.max_worker_id,
           task_profile.max_solve_seconds,
           task_profile.max_pool_percent,
           task_profile.serial_tail_seconds,
           task_profile.serial_tail_percent,
           task_profile.tail_index_order,
           task_profile.tail_worker_id);
  }

  if (!phase_profile.executed) {
    logverb("[index-shard] phase-profile executed=0\n");
  } else if (phase_profile.quantiles_available) {
    logverb("[index-shard] phase-profile executed=%zu "
           "task_work_wall_sum=%.6f "
           "reset_work_wall_sum=%.6f reset_percent=%.1f "
           "acquire_work_wall_sum=%.6f acquire_percent=%.1f "
           "solve_work_wall_sum=%.6f solve_percent=%.1f "
           "analyze_work_wall_sum=%.6f analyze_percent=%.1f "
           "release_work_wall_sum=%.6f release_percent=%.1f "
           "other_work_wall_sum=%.6f other_percent=%.1f "
           "acquire_p50=%.6f acquire_p90=%.6f "
           "acquire_p99=%.6f acquire_max=%.6f "
           "solve_p50=%.6f solve_p90=%.6f "
           "solve_p99=%.6f solve_max=%.6f\n",
           phase_profile.executed,
           phase_profile.task_wall_total,
           phase_profile.reset_total,
           phase_profile.reset_percent,
           phase_profile.acquire_total,
           phase_profile.acquire_percent,
           phase_profile.solve_total,
           phase_profile.solve_percent,
           phase_profile.analyze_total,
           phase_profile.analyze_percent,
           phase_profile.release_total,
           phase_profile.release_percent,
           phase_profile.other_total,
           phase_profile.other_percent,
           phase_profile.acquire_p50,
           phase_profile.acquire_p90,
           phase_profile.acquire_p99,
           phase_profile.acquire_max,
           phase_profile.solve_p50,
           phase_profile.solve_p90,
           phase_profile.solve_p99,
           phase_profile.solve_max);
  } else {
    logverb("[index-shard] phase-profile executed=%zu "
           "task_work_wall_sum=%.6f "
           "reset_work_wall_sum=%.6f reset_percent=%.1f "
           "acquire_work_wall_sum=%.6f acquire_percent=%.1f "
           "solve_work_wall_sum=%.6f solve_percent=%.1f "
           "analyze_work_wall_sum=%.6f analyze_percent=%.1f "
           "release_work_wall_sum=%.6f release_percent=%.1f "
           "other_work_wall_sum=%.6f other_percent=%.1f "
           "quantiles=unavailable\n",
           phase_profile.executed,
           phase_profile.task_wall_total,
           phase_profile.reset_total,
           phase_profile.reset_percent,
           phase_profile.acquire_total,
           phase_profile.acquire_percent,
           phase_profile.solve_total,
           phase_profile.solve_percent,
           phase_profile.analyze_total,
           phase_profile.analyze_percent,
           phase_profile.release_total,
           phase_profile.release_percent,
           phase_profile.other_total,
           phase_profile.other_percent);
  }

  if (pass_metrics.resource_available) {
    logverb("[index-shard] pass-resource user=%.6f sys=%.6f "
           "minflt=%ld majflt=%ld nvcsw=%ld nivcsw=%ld "
           "inblock=%ld oublock=%ld\n",
           pass_metrics.user_seconds,
           pass_metrics.system_seconds,
           pass_metrics.minor_faults,
           pass_metrics.major_faults,
           pass_metrics.voluntary_context_switches,
           pass_metrics.involuntary_context_switches,
           pass_metrics.filesystem_input_blocks,
           pass_metrics.filesystem_output_blocks);
  } else {
    logverb("[index-shard] pass-resource unavailable\n");
  }

  pthread_mutex_lock(&pool->shared.queue_mutex);
  pool->shared.outer_states = NULL;
  pool->shared.producer_width = 0U;
  pool->shared.helper_width = 0U;
  pthread_mutex_unlock(&pool->shared.queue_mutex);
  pool->shared.worker_view = NULL;
  hooks->destroy_worker_view(worker_view);
  worker_view = NULL;
  free(outer_states);
  index_shard_pool_release_pass(pool);
  return status;
}
