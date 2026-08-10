/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "solver_codekd_internal.h"
#include "solver_codekd_test_private.h"
#include "verify_prepared_internal.h"

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

int solver_test_codekd_detached_initial_planning(void) {
    solver_codekd_test_detached_result_t result;

    memset(&result, 0, sizeof(result));
    if (solver_codekd_test_run_detached_initial_planning(&result)) {
        return -1;
    }
    if (!result.supported) {
        return 0;
    }
    if (!result.plan_prepared_before_submit ||
        !result.queued_callback ||
        !result.success_compute_ready ||
        !result.repeated_cycles ||
        !result.precomputed_quad_revalidated ||
        !result.fully_resident_compute_ready ||
        !result.empty_owner_replay ||
        !result.callback_error ||
        !result.eagain_retry ||
        !result.refusal_owner_fallback ||
        !result.cancellation_stopped) {
        return -1;
    }
    return 0;
}

int solver_test_verification_retirement_horizon(void) {
    solver_codekd_packet_task_input_t input;
    solver_codekd_search_packet_t packet;
    solver_candidate_delivery_record_t records[
        2U * SOLVER_VERIFY_RETIRE_HORIZON + 3U];
    solver_codekd_result_slot_t slots[1];
    solver_codekd_page_workspace_t workspace;
    verify_field_t field;
    startree_t starkd;
    kdtree_t tree;
    fitsbin_t source;
    u32 inds[2U * SOLVER_VERIFY_RETIRE_HORIZON + 3U];
    double sdists[2U * SOLVER_VERIFY_RETIRE_HORIZON + 3U];
    size_t full_end = sizeof(records) / sizeof(records[0]);
    size_t query_bytes = 0U;
    size_t candidate_index;
    index_shard_helper_task_status_t status;
    anbool more_work = FALSE;
    int result = -1;

    memset(&input, 0, sizeof(input));
    memset(&packet, 0, sizeof(packet));
    memset(records, 0, sizeof(records));
    memset(slots, 0, sizeof(slots));
    memset(&workspace, 0, sizeof(workspace));
    memset(&field, 0, sizeof(field));
    memset(&starkd, 0, sizeof(starkd));
    memset(&tree, 0, sizeof(tree));
    memset(&source, 0, sizeof(source));
    memset(inds, 0, sizeof(inds));
    memset(sdists, 0, sizeof(sdists));
    input.verification.enabled = TRUE;
    input.verification.field = &field;
    input.verification.verify_pix = 1.0;
    input.verification.distractor_ratio = 0.25;
    input.field_maxx = 1.0;
    input.field_maxy = 1.0;
    packet.candidate_records = records;
    packet.candidate_capacity = full_end;
    packet.candidate_count = full_end;
    packet.candidate_window_count = full_end;
    packet.slots = slots;
    packet.count = 1U;
    packet.hit_capacity = full_end;
    packet.hit_count = full_end;
    packet.inds = inds;
    packet.sdists = sdists;
    packet.dimquads = DQMAX;
    packet.verify_query_budget =
        SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES;
    packet.verify_plan_end = full_end;
    packet.verify_topology_end = full_end;
    packet.verify_plan_complete = TRUE;
    packet.state =
        SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
    packet.verification_page_submitted = 1U;
    slots[0].hit_count = (unsigned int)full_end;
    slots[0].state = SOLVER_CODEKD_RESULT_READY;
    tree.io = &source;
    tree.io_is_fitsbin = TRUE;
    starkd.tree = &tree;
    packet.starkd = &starkd;
    packet.page_workspace = &workspace;

    for (candidate_index = 0U;
         candidate_index < full_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &records[candidate_index];
        verify_index_query_t* query = calloc(1, sizeof(*query));

        if (!query) {
            goto cleanup;
        }
        query->source = &starkd;
        record->quadid = (unsigned int)candidate_index;
        record->plan_action = SOLVER_AB_CANDIDATE_VERIFY;
        record->candidate_prepared = TRUE;
        record->prepared_scale = 1.0;
        record->prepared_wcs.cd[0][0] = 1.0;
        record->prepared_wcs.cd[1][1] = 1.0;
        record->prepared_wcs.imagew = 1.0;
        record->prepared_wcs.imageh = 1.0;
        record->verify_query = query;
        record->verify_query_captured = TRUE;
        inds[candidate_index] = record->quadid;
        if (!query_bytes) {
            query_bytes = verify_index_query_bytes(query);
        }
    }
    if (!query_bytes || query_bytes == SIZE_MAX ||
        full_end > SIZE_MAX / query_bytes) {
        goto cleanup;
    }
    packet.verify_query_bytes = full_end * query_bytes;

    status = solver_codekd_packet_query_verification_ready(
        &packet, &more_work);
    if (status != INDEX_SHARD_HELPER_TASK_OK || !more_work ||
        packet.state !=
            SOLVER_CODEKD_PACKET_VERIFY_PREPARE_COMPUTE_READY ||
        packet.verify_plan_first != 0U ||
        packet.verify_plan_end != SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_topology_end != full_end ||
        packet.verify_query_bytes != full_end * query_bytes ||
        packet.verification_page_submitted != 1U ||
        packet.delivery_ticket || packet.delivery_source) {
        goto cleanup;
    }
    status = solver_codekd_packet_prepare_and_score_verification_ready(
        &input, &packet);
    if (status != INDEX_SHARD_HELPER_TASK_OK ||
        packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet.verify_plan_first != 0U ||
        packet.verify_plan_end != SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_topology_end != full_end ||
        packet.verify_prepared_count !=
            SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_query_bytes !=
            (full_end - SOLVER_VERIFY_RETIRE_HORIZON) *
                query_bytes ||
        packet.verification_page_submitted != 1U ||
        solver_codekd_packet_retired_verification_complete(&packet) ==
            0) {
        goto cleanup;
    }
    for (candidate_index = 0U;
         candidate_index < full_end;
         candidate_index++) {
        if (candidate_index < SOLVER_VERIFY_RETIRE_HORIZON) {
            if (records[candidate_index].verify_query ||
                records[candidate_index].verify_query_captured == TRUE) {
                goto cleanup;
            }
        } else if (!records[candidate_index].verify_query ||
                   records[candidate_index].verify_query_captured != TRUE) {
            goto cleanup;
        }
    }
    for (candidate_index = packet.verify_plan_first;
         candidate_index < packet.verify_plan_end;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &records[candidate_index]);
    }
    if (solver_codekd_packet_retired_verification_complete(&packet)) {
        goto cleanup;
    }
    packet.candidate_window_offset = packet.verify_plan_end;
    if (solver_codekd_packet_rearm_verification_owner(&packet) ||
        packet.state !=
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY ||
        packet.verify_plan_first != SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end != full_end ||
        packet.verify_topology_end != full_end ||
        packet.verify_query_bytes !=
            (full_end - SOLVER_VERIFY_RETIRE_HORIZON) *
                query_bytes) {
        goto cleanup;
    }
    status = solver_codekd_packet_query_verification_ready(
        &packet, &more_work);
    if (status != INDEX_SHARD_HELPER_TASK_OK || !more_work ||
        packet.state !=
            SOLVER_CODEKD_PACKET_VERIFY_PREPARE_COMPUTE_READY ||
        packet.verify_plan_first != SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end !=
            2U * SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_topology_end != full_end ||
        packet.verify_query_bytes !=
            (full_end - SOLVER_VERIFY_RETIRE_HORIZON) *
                query_bytes ||
        packet.verification_page_submitted != 1U ||
        packet.delivery_ticket || packet.delivery_source) {
        goto cleanup;
    }
    status = solver_codekd_packet_prepare_and_score_verification_ready(
        &input, &packet);
    if (status != INDEX_SHARD_HELPER_TASK_OK ||
        packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet.verify_plan_first != SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end !=
            2U * SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_topology_end != full_end ||
        packet.verify_prepared_count !=
            SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_query_bytes != 3U * query_bytes ||
        packet.verification_page_submitted != 1U) {
        goto cleanup;
    }
    for (candidate_index = 0U;
         candidate_index < full_end;
         candidate_index++) {
        if (candidate_index <
                2U * SOLVER_VERIFY_RETIRE_HORIZON) {
            if (records[candidate_index].verify_query ||
                records[candidate_index].verify_query_captured == TRUE) {
                goto cleanup;
            }
        } else if (!records[candidate_index].verify_query ||
                   records[candidate_index].verify_query_captured != TRUE) {
            goto cleanup;
        }
    }
    for (candidate_index = packet.verify_plan_first;
         candidate_index < packet.verify_plan_end;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &records[candidate_index]);
    }
    if (solver_codekd_packet_retired_verification_complete(&packet)) {
        goto cleanup;
    }
    packet.candidate_window_offset = packet.verify_plan_end;
    if (solver_codekd_packet_rearm_verification_owner(&packet) ||
        packet.state !=
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY ||
        packet.verify_plan_first !=
            2U * SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end != full_end ||
        packet.verify_topology_end != full_end ||
        packet.verify_query_bytes != 3U * query_bytes) {
        goto cleanup;
    }
    more_work = FALSE;
    status = solver_codekd_packet_query_verification_ready(
        &packet, &more_work);
    if (status != INDEX_SHARD_HELPER_TASK_OK || !more_work ||
        packet.state !=
            SOLVER_CODEKD_PACKET_VERIFY_PREPARE_COMPUTE_READY ||
        packet.verify_plan_first !=
            2U * SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end != full_end ||
        packet.verify_topology_end != full_end ||
        packet.verify_query_bytes != 3U * query_bytes ||
        packet.verification_page_submitted != 1U ||
        packet.delivery_ticket || packet.delivery_source) {
        goto cleanup;
    }
    status = solver_codekd_packet_prepare_and_score_verification_ready(
        &input, &packet);
    if (status != INDEX_SHARD_HELPER_TASK_OK ||
        packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet.verify_plan_first !=
            2U * SOLVER_VERIFY_RETIRE_HORIZON ||
        packet.verify_plan_end != full_end ||
        packet.verify_topology_end != full_end ||
        packet.verify_prepared_count != 3U ||
        packet.verify_query_bytes != 0U ||
        packet.verification_page_submitted != 1U) {
        goto cleanup;
    }
    for (candidate_index = 0U;
         candidate_index < full_end;
         candidate_index++) {
        if (records[candidate_index].verify_query ||
            records[candidate_index].verify_query_captured == TRUE) {
            goto cleanup;
        }
    }
    for (candidate_index = packet.verify_plan_first;
         candidate_index < packet.verify_plan_end;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &records[candidate_index]);
    }
    if (solver_codekd_packet_retired_verification_complete(&packet)) {
        goto cleanup;
    }
    result = 0;

cleanup:
    for (candidate_index = 0U;
         candidate_index < full_end;
         candidate_index++) {
        solver_codekd_record_clear_verification_speculation(
            &records[candidate_index]);
    }
    return result;
}

int solver_test_captured_verification_requires_owned_sweep(void) {
    verify_index_query_t* query;
    verify_prepared_hit_t* prepared = NULL;
    verify_prepared_score_t score;
    verify_field_t field;
    MatchObj match;
    double* retained_refxyz;
    int* retained_refstarid;
    uint8_t* retained_sweep;

    memset(&field, 0, sizeof(field));
    memset(&match, 0, sizeof(match));
    memset(&score, 0, sizeof(score));
    match.wcs_valid = TRUE;
    match.radius = 0.5;
    match.wcstan.crval[0] = 0.0;
    match.wcstan.crval[1] = 0.0;
    match.wcstan.crpix[0] = 0.5;
    match.wcstan.crpix[1] = 0.5;
    match.wcstan.cd[0][0] = 1.0;
    match.wcstan.cd[1][1] = 1.0;
    match.wcstan.imagew = 1.0;
    match.wcstan.imageh = 1.0;

    query = calloc(1, sizeof(*query));
    if (!query) {
        return -1;
    }
    query->source_nstars = 1;
    query->radius2 = match.radius * match.radius;
    query->nrall = 1;
    query->refxyz = calloc(3U, sizeof(*query->refxyz));
    query->refstarid = calloc(1U, sizeof(*query->refstarid));
    if (!query->refxyz || !query->refstarid) {
        verify_destroy_index_query(query);
        return -1;
    }
    if (verify_prepare_captured_hit_from_query(
            &query, 0, &match, NULL, &field,
            1.0, 0.25, 1.0, 1.0,
            -1.0, 1.0, 2.0,
            FALSE, FALSE, &prepared) != -1 ||
        !query || prepared) {
        verify_destroy_index_query(query);
        verify_destroy_prepared_hit(prepared);
        return -1;
    }
    verify_destroy_index_query(query);

    query = calloc(1, sizeof(*query));
    if (!query) {
        return -1;
    }
    query->source = (const startree_t*)query;
    query->source_nstars = 1;
    query->radius2 = match.radius * match.radius;
    query->nrall = 1;
    query->refxyz = calloc(3U, sizeof(*query->refxyz));
    query->refstarid = calloc(1U, sizeof(*query->refstarid));
    query->sweep = calloc(1U, sizeof(*query->sweep));
    if (!query->refxyz || !query->refstarid || !query->sweep) {
        verify_destroy_index_query(query);
        return -1;
    }
    query->refxyz[0] = -1.0;
    if (verify_prepare_captured_hit_from_query(
            &query, 0, &match, NULL, &field,
            1.0, 0.25, 1.0, 1.0,
            -1.0, 1.0, 2.0,
            FALSE, FALSE, &prepared) ||
        query || !prepared ||
        verify_score_prepared_hit(prepared, &score) ||
        verify_finish_prepared_hit(prepared, &score, &match)) {
        verify_destroy_index_query(query);
        verify_destroy_prepared_score(&score);
        verify_destroy_prepared_hit(prepared);
        return -1;
    }
    verify_destroy_prepared_hit(prepared);
    prepared = NULL;

    memset(&match, 0, sizeof(match));
    match.wcs_valid = TRUE;
    match.radius = 0.5;
    query = calloc(1, sizeof(*query));
    if (!query) {
        return -1;
    }
    query->source_nstars = 1;
    query->radius2 = match.radius * match.radius;
    query->nrall = 1;
    query->refxyz = calloc(3U, sizeof(*query->refxyz));
    query->refstarid = calloc(1U, sizeof(*query->refstarid));
    query->sweep = calloc(1U, sizeof(*query->sweep));
    if (!query->refxyz || !query->refstarid || !query->sweep) {
        verify_destroy_index_query(query);
        return -1;
    }
    query->refstarid[0] = query->source_nstars;
    retained_refxyz = query->refxyz;
    retained_refstarid = query->refstarid;
    retained_sweep = query->sweep;
    if (verify_prepare_captured_hit_from_query(
            &query, 0, &match, NULL, &field,
            1.0, 0.25, 1.0, 1.0,
            -1.0, 1.0, 2.0,
            FALSE, FALSE, &prepared) != -1 ||
        !query || prepared ||
        query->refxyz != retained_refxyz ||
        query->refstarid != retained_refstarid ||
        query->sweep != retained_sweep) {
        verify_destroy_index_query(query);
        verify_destroy_prepared_hit(prepared);
        return -1;
    }
    verify_destroy_index_query(query);
    return 0;
}
