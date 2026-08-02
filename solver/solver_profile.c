/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* Solver profile aggregation and trace reporting. */

#include <limits.h>

#include "solver.h"
#include "log.h"
#include "index_shard_internal.h"
#include "solver_inline_internal.h"
#include "solver_profile_internal.h"
void solver_profile_accumulate(solver_profile_t* total,
                               const solver_profile_t* profile) {
    if (!total || !profile) {
        return;
    }

    total->detailed = total->detailed || profile->detailed;
    total->execution_failed =
        total->execution_failed || profile->execution_failed;

    total->solver_run_wall_seconds +=
        profile->solver_run_wall_seconds;
    total->codekd_wall_seconds += profile->codekd_wall_seconds;
    total->resolve_wall_seconds += profile->resolve_wall_seconds;
    total->verify_wall_seconds += profile->verify_wall_seconds;
    total->verification_score_wall_seconds +=
        profile->verification_score_wall_seconds;

    total->codekd_calls += profile->codekd_calls;
    total->codekd_hits += profile->codekd_hits;
    total->resolve_calls += profile->resolve_calls;
    total->verify_calls += profile->verify_calls;
    total->hypothesis_batches += profile->hypothesis_batches;
    total->hypothesis_batches_completed +=
        profile->hypothesis_batches_completed;
    total->hypothesis_batches_stopped +=
        profile->hypothesis_batches_stopped;
    total->hypothesis_batches_failed +=
        profile->hypothesis_batches_failed;
    total->hypotheses_generated += profile->hypotheses_generated;
    total->hypotheses_executed += profile->hypotheses_executed;
    total->hypotheses_reduced += profile->hypotheses_reduced;
    total->task_ranges_planned += profile->task_ranges_planned;
    total->task_ranges_executed += profile->task_ranges_executed;
    total->task_ranges_submitted += profile->task_ranges_submitted;
    total->task_ranges_inline += profile->task_ranges_inline;
    total->parallel_batches += profile->parallel_batches;
    total->parallel_batches_observed +=
        profile->parallel_batches_observed;
    total->parallel_hypotheses += profile->parallel_hypotheses;
    total->allocation_failures += profile->allocation_failures;
    total->search_failures += profile->search_failures;
    total->ab_helper_tasks += profile->ab_helper_tasks;
    total->ab_helper_combinations +=
        profile->ab_helper_combinations;
    total->page_plan_descriptors_total +=
        profile->page_plan_descriptors_total;
    total->page_plan_descriptors_complete +=
        profile->page_plan_descriptors_complete;
    total->page_plan_descriptor_splits +=
        profile->page_plan_descriptor_splits;
    total->page_plan_raw_ranges += profile->page_plan_raw_ranges;
    total->page_plan_unique_pages += profile->page_plan_unique_pages;
    total->page_plan_ranges_after_dedup +=
        profile->page_plan_ranges_after_dedup;
    total->page_plan_logical_bytes +=
        profile->page_plan_logical_bytes;
    total->page_plan_aligned_bytes +=
        profile->page_plan_aligned_bytes;
    total->page_plan_overread_bytes +=
        profile->page_plan_overread_bytes;
    total->page_plan_not_applicable +=
        profile->page_plan_not_applicable;
    total->page_plan_allocation_refused +=
        profile->page_plan_allocation_refused;
    total->page_plan_source_mismatch +=
        profile->page_plan_source_mismatch;
    total->page_plan_invalid_range +=
        profile->page_plan_invalid_range;
    total->page_plan_byte_budget_refused +=
        profile->page_plan_byte_budget_refused;
    total->page_plan_range_capacity_refused +=
        profile->page_plan_range_capacity_refused;
    total->page_plan_service_refused +=
        profile->page_plan_service_refused;
    total->page_plan_service_errors +=
        profile->page_plan_service_errors;
    total->page_plan_cancelled += profile->page_plan_cancelled;
    total->candidate_delivery_candidates +=
        profile->candidate_delivery_candidates;
    total->candidate_quad_submitted +=
        profile->candidate_quad_submitted;
    total->candidate_quad_ready +=
        profile->candidate_quad_ready;
    total->candidate_quad_fallback +=
        profile->candidate_quad_fallback;
    total->candidate_star_submitted +=
        profile->candidate_star_submitted;
    total->candidate_star_ready +=
        profile->candidate_star_ready;
    total->candidate_star_fallback +=
        profile->candidate_star_fallback;
    total->candidate_delivery_windows +=
        profile->candidate_delivery_windows;
    total->candidate_quad_ready_rows +=
        profile->candidate_quad_ready_rows;
    total->candidate_star_ready_rows +=
        profile->candidate_star_ready_rows;
    total->candidate_retired_rows +=
        profile->candidate_retired_rows;
    total->candidate_native_rows +=
        profile->candidate_native_rows;
    total->verification_page_queries +=
        profile->verification_page_queries;
    total->verification_page_queries_planned +=
        profile->verification_page_queries_planned;
    total->verification_page_prefixes +=
        profile->verification_page_prefixes;
    total->verification_page_submitted +=
        profile->verification_page_submitted;
    total->verification_page_ready +=
        profile->verification_page_ready;
    total->verification_page_fallback +=
        profile->verification_page_fallback;
    total->verification_page_ready_rows +=
        profile->verification_page_ready_rows;
    total->verification_page_ranges +=
        profile->verification_page_ranges;
    total->verification_page_logical_bytes +=
        profile->verification_page_logical_bytes;
    total->verification_page_aligned_bytes +=
        profile->verification_page_aligned_bytes;
    total->candidate_math_prepared +=
        profile->candidate_math_prepared;
    total->candidate_math_reused +=
        profile->candidate_math_reused;
    total->verification_score_batches_prepared +=
        profile->verification_score_batches_prepared;
    total->verification_score_contexts_prepared +=
        profile->verification_score_contexts_prepared;
    total->verification_score_batches_executed +=
        profile->verification_score_batches_executed;
    total->verification_score_contexts_completed +=
        profile->verification_score_contexts_completed;
    if (ULLONG_MAX -
            total->verification_score_work_units_completed <
        profile->verification_score_work_units_completed) {
        total->verification_score_work_units_completed =
            ULLONG_MAX;
    } else {
        total->verification_score_work_units_completed +=
            profile->verification_score_work_units_completed;
    }
    total->verification_score_fallback_batches +=
        profile->verification_score_fallback_batches;
    total->verification_score_stopped_batches +=
        profile->verification_score_stopped_batches;
    total->staged_owner_claims += profile->staged_owner_claims;
    total->staged_foreign_claims += profile->staged_foreign_claims;
    total->staged_io_submitted += profile->staged_io_submitted;
    total->staged_io_completed += profile->staged_io_completed;
    if (!total->hypothesis_order_hash) {
        total->hypothesis_order_hash =
            profile->hypothesis_order_hash;
    } else if (profile->hypothesis_order_hash) {
        total->hypothesis_order_hash =
            solver_order_hash_mix(
                total->hypothesis_order_hash,
                profile->hypothesis_order_hash);
    }
    if (!total->kd_result_order_hash) {
        total->kd_result_order_hash =
            profile->kd_result_order_hash;
    } else if (profile->kd_result_order_hash) {
        total->kd_result_order_hash =
            solver_order_hash_mix(
                total->kd_result_order_hash,
                profile->kd_result_order_hash);
    }
    if (!total->candidate_order_hash) {
        total->candidate_order_hash =
            profile->candidate_order_hash;
    } else if (profile->candidate_order_hash) {
        total->candidate_order_hash =
            solver_order_hash_mix(
                total->candidate_order_hash,
                profile->candidate_order_hash);
    }

    if (profile->max_batch_hypotheses >
        total->max_batch_hypotheses) {
        total->max_batch_hypotheses =
            profile->max_batch_hypotheses;
    }

    if (profile->max_task_ranges > total->max_task_ranges) {
        total->max_task_ranges = profile->max_task_ranges;
    }

    if (profile->max_parallel_ranges >
        total->max_parallel_ranges) {
        total->max_parallel_ranges = profile->max_parallel_ranges;
    }
    if (profile->max_staged_io_submitted >
        total->max_staged_io_submitted) {
        total->max_staged_io_submitted =
            profile->max_staged_io_submitted;
    }
    if (profile->max_staged_compute_ready >
        total->max_staged_compute_ready) {
        total->max_staged_compute_ready =
            profile->max_staged_compute_ready;
    }
}

void solver_profile_report(const solver_t* solver) {
    const solver_profile_t* profile;
    double reduction_wall_seconds;

    if (!solver || !index_shard_trace_enabled()) {
        return;
    }

    profile = &solver->profile;
    reduction_wall_seconds =
        profile->resolve_wall_seconds -
        profile->verify_wall_seconds;

    if (reduction_wall_seconds < 0.0) {
        reduction_wall_seconds = 0.0;
    }

    logmsg("[solver] phase-profile detailed=%i failed=%i "
           "solver_run_elapsed=%.6f codekd_work_wall_sum=%.6f "
           "codekd_calls=%llu codekd_hits=%llu "
           "resolve_work_wall_sum=%.6f "
           "reduction_ex_verify_hit_work_wall_sum=%.6f "
           "resolve_calls=%llu verify_hit_work_wall_sum=%.6f "
           "verify_calls=%llu "
           "hypothesis_batches=%llu parallel_batches=%llu "
           "parallel_batches_observed=%llu "
           "parallel_hypotheses=%llu "
           "task_ranges_planned=%llu task_ranges_executed=%llu "
           "task_ranges_inline=%llu allocation_failures=%llu "
           "search_failures=%llu "
           "helper_tasks=%llu helper_combinations=%llu "
           "hypothesis_order=%016llx "
           "kd_result_order=%016llx candidate_order=%016llx\n",
           profile->detailed ? 1 : 0,
           profile->execution_failed ? 1 : 0,
           profile->solver_run_wall_seconds,
           profile->codekd_wall_seconds,
           profile->codekd_calls,
           profile->codekd_hits,
           profile->resolve_wall_seconds,
           reduction_wall_seconds,
           profile->resolve_calls,
           profile->verify_wall_seconds,
           profile->verify_calls,
           profile->hypothesis_batches,
           profile->parallel_batches,
           profile->parallel_batches_observed,
           profile->parallel_hypotheses,
           profile->task_ranges_planned,
           profile->task_ranges_executed,
           profile->task_ranges_inline,
           profile->allocation_failures,
           profile->search_failures,
           profile->ab_helper_tasks,
           profile->ab_helper_combinations,
           profile->hypothesis_order_hash,
           profile->kd_result_order_hash,
           profile->candidate_order_hash);
    logmsg("[solver] page-pipeline descriptors=%llu complete=%llu "
           "boundary_deferrals=%llu raw_hints=%llu unique_pages=%llu "
           "coalesced_ranges=%llu raw_hint_bytes=%llu aligned_bytes=%llu "
           "positive_alignment_delta_bytes=%llu "
           "refusals=%llu/%llu/%llu/%llu/%llu/"
           "%llu/%llu/%llu/%llu staged_claims=%llu/%llu "
           "io=%llu/%llu max_io=%zu max_ready=%zu "
           "candidate_delivery=%llu quad=%llu/%llu/%llu "
           "star=%llu/%llu/%llu windows=%llu "
           "rows=%llu/%llu/%llu/%llu "
           "verify_pages=%llu/%llu prefixes=%llu tickets=%llu/%llu "
           "fallback=%llu ready_rows=%llu ranges=%llu "
           "bytes=%llu/%llu candidate_math=%llu/%llu "
           "verify_score=%llu/%llu/%llu/%llu/%llu/%llu "
           "work=%llu work_wall_sum=%.6f\n",
           profile->page_plan_descriptors_total,
           profile->page_plan_descriptors_complete,
           profile->page_plan_descriptor_splits,
           profile->page_plan_raw_ranges,
           profile->page_plan_unique_pages,
           profile->page_plan_ranges_after_dedup,
           profile->page_plan_logical_bytes,
           profile->page_plan_aligned_bytes,
           profile->page_plan_overread_bytes,
           profile->page_plan_not_applicable,
           profile->page_plan_allocation_refused,
           profile->page_plan_source_mismatch,
           profile->page_plan_invalid_range,
           profile->page_plan_byte_budget_refused,
           profile->page_plan_range_capacity_refused,
           profile->page_plan_service_refused,
           profile->page_plan_service_errors,
           profile->page_plan_cancelled,
           profile->staged_owner_claims,
           profile->staged_foreign_claims,
           profile->staged_io_submitted,
           profile->staged_io_completed,
           profile->max_staged_io_submitted,
           profile->max_staged_compute_ready,
           profile->candidate_delivery_candidates,
           profile->candidate_quad_submitted,
           profile->candidate_quad_ready,
           profile->candidate_quad_fallback,
           profile->candidate_star_submitted,
           profile->candidate_star_ready,
           profile->candidate_star_fallback,
           profile->candidate_delivery_windows,
           profile->candidate_quad_ready_rows,
           profile->candidate_star_ready_rows,
           profile->candidate_retired_rows,
           profile->candidate_native_rows,
           profile->verification_page_queries,
           profile->verification_page_queries_planned,
           profile->verification_page_prefixes,
           profile->verification_page_submitted,
           profile->verification_page_ready,
           profile->verification_page_fallback,
           profile->verification_page_ready_rows,
           profile->verification_page_ranges,
           profile->verification_page_logical_bytes,
           profile->verification_page_aligned_bytes,
           profile->candidate_math_prepared,
           profile->candidate_math_reused,
           profile->verification_score_batches_prepared,
           profile->verification_score_contexts_prepared,
           profile->verification_score_batches_executed,
           profile->verification_score_contexts_completed,
           profile->verification_score_fallback_batches,
           profile->verification_score_stopped_batches,
           profile->verification_score_work_units_completed,
           profile->verification_score_wall_seconds);
}
