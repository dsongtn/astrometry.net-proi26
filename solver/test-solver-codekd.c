/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <string.h>

#include "solver_codekd_test_private.h"

int solver_test_candidate_rolling_windows(void) {
    static const size_t expected_counts[] = { 128U, 101U };
    solver_codekd_test_window_result_t result;

    memset(&result, 0, sizeof(result));
    if (solver_codekd_test_run_candidate_windows(&result) ||
        result.window_count !=
            sizeof(expected_counts) / sizeof(expected_counts[0]) ||
        memcmp(result.window_sizes,
               expected_counts,
               sizeof(expected_counts)) ||
        result.retired_candidates != result.candidate_count ||
        result.candidate_cursor != result.candidate_count ||
        result.retire_descriptor != 3U ||
        result.delivery_windows != result.window_count) {
        return -1;
    }
    return 0;
}

int solver_test_candidate_nonresident_zero_submit_falls_back(void) {
    solver_codekd_test_fallback_result_t result;

    memset(&result, 0, sizeof(result));
    if (solver_codekd_test_run_nonresident_fallback(&result) ||
        result.submit_status != -1 ||
        !result.results_ready ||
        result.quad_compute_ready ||
        result.quad_fallback != 1U ||
        !result.quad_delivery_disabled ||
        !result.star_delivery_disabled ||
        result.quad_submitted ||
        result.quad_ready ||
        result.has_delivery_ticket ||
        result.has_delivery_source) {
        return -1;
    }
    return 0;
}

int solver_test_verification_packet_bounds(void) {
    solver_codekd_test_reserve_result_t result;

    memset(&result, 0, sizeof(result));
    if (solver_codekd_test_run_verification_reserve(&result) ||
        result.initial_status != SOLVER_CODEKD_TEST_RESERVE_OK ||
        result.initial_capacity > result.maximum_capacity ||
        result.oversized_status != SOLVER_CODEKD_TEST_RESERVE_FULL ||
        result.allocation_failed ||
        result.final_capacity > result.maximum_capacity) {
        return -1;
    }
    return 0;
}
