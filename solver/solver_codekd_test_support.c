/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <string.h>

#include "solver_codekd_internal.h"
#include "solver_codekd_test_private.h"

int solver_codekd_test_run_candidate_windows(
    solver_codekd_test_window_result_t* result) {
    solver_codekd_search_packet_t packet;
    solver_codekd_result_slot_t slots[3];
    solver_candidate_delivery_record_t records[
        SOLVER_CANDIDATE_DELIVERY_LIMIT];
    size_t window_index;
    size_t retired = 0U;

    if (!result) {
        return -1;
    }
    memset(result, 0, sizeof(*result));
    memset(&packet, 0, sizeof(packet));
    memset(slots, 0, sizeof(slots));
    memset(records, 0, sizeof(records));
    slots[0].hit_first = 0U;
    slots[0].hit_count = 17U;
    slots[0].state = SOLVER_CODEKD_RESULT_READY;
    slots[1].hit_first = 17U;
    slots[1].hit_count = 200U;
    slots[1].state = SOLVER_CODEKD_RESULT_READY;
    slots[2].hit_first = 217U;
    slots[2].hit_count = 12U;
    slots[2].state = SOLVER_CODEKD_RESULT_READY;
    packet.slots = slots;
    packet.count = 3U;
    packet.hit_count = 229U;
    packet.candidate_count = packet.hit_count;
    packet.candidate_capacity =
        SOLVER_CANDIDATE_DELIVERY_LIMIT;
    packet.candidate_records = records;

    for (window_index = 0U; window_index < 2U; window_index++) {
        size_t slot_end;

        if (solver_codekd_packet_begin_candidate_window(&packet)) {
            return -1;
        }
        result->window_sizes[result->window_count++] =
            packet.candidate_window_count;
        retired += packet.candidate_window_count;
        packet.candidate_cursor += packet.candidate_window_count;
        packet.candidate_window_offset =
            packet.candidate_window_count;
        solver_codekd_packet_clear_candidate_window(&packet);
        while (packet.retire_descriptor < packet.count) {
            slot_end =
                slots[packet.retire_descriptor].hit_first +
                slots[packet.retire_descriptor].hit_count;
            if (packet.candidate_cursor < slot_end) {
                break;
            }
            packet.retire_descriptor++;
        }
    }

    result->retired_candidates = retired;
    result->candidate_count = packet.candidate_count;
    result->candidate_cursor = packet.candidate_cursor;
    result->retire_descriptor = packet.retire_descriptor;
    result->delivery_windows = packet.candidate_delivery_windows;
    return 0;
}

int solver_codekd_test_run_nonresident_fallback(
    solver_codekd_test_fallback_result_t* result) {
    fitsbin_t source;
    quadfile_t quads;
    solver_codekd_search_packet_t packet;
    solver_candidate_delivery_record_t record;
    uint32_t quadrow[DQMAX];
    u32 quadid = 0U;

    if (!result) {
        return -1;
    }
    memset(result, 0, sizeof(*result));
    memset(&source, 0, sizeof(source));
    memset(&quads, 0, sizeof(quads));
    memset(&packet, 0, sizeof(packet));
    memset(&record, 0, sizeof(record));
    memset(quadrow, 0, sizeof(quadrow));
    source.mmap_prefetch_enabled = FALSE;
    source.payload_fully_resident = FALSE;
    quads.numquads = 1U;
    quads.dimquads = DQMAX;
    quads.fb = &source;
    quads.quadarray = quadrow;
    packet.candidate_records = &record;
    packet.inds = &quadid;
    packet.candidate_count = 1U;
    packet.candidate_capacity = 1U;
    packet.candidate_window_count = 1U;
    packet.state = SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY;

    result->submit_status =
        solver_codekd_packet_submit_candidate_pages(&packet);
    result->results_ready =
        packet.state == SOLVER_CODEKD_PACKET_RESULTS_READY;
    result->quad_compute_ready =
        packet.state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
    result->quad_fallback = packet.candidate_quad_fallback;
    result->quad_delivery_disabled =
        packet.candidate_quad_delivery_disabled;
    result->star_delivery_disabled =
        packet.candidate_star_delivery_disabled;
    result->quad_submitted = packet.candidate_quad_submitted;
    result->quad_ready = packet.candidate_quad_ready;
    result->has_delivery_ticket = packet.delivery_ticket != NULL;
    result->has_delivery_source = packet.delivery_source != NULL;
    return 0;
}

int solver_codekd_test_run_verification_reserve(
    solver_codekd_test_reserve_result_t* result) {
    solver_verification_packet_t packet;
    size_t impossible_candidates =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES /
            sizeof(solver_ab_candidate_t) +
        1U;

    if (!result) {
        return -1;
    }
    memset(result, 0, sizeof(*result));
    memset(&packet, 0, sizeof(packet));
    result->maximum_capacity =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES /
        sizeof(*packet.candidates);
    result->initial_status =
        (solver_codekd_test_reserve_status_t)
        solver_verification_packet_reserve(&packet, 1U);
    result->initial_capacity = packet.candidate_capacity;
    if (result->initial_status != SOLVER_CODEKD_TEST_RESERVE_OK ||
        result->initial_capacity > result->maximum_capacity) {
        result->allocation_failed = packet.allocation_failed;
        solver_verification_packet_free(&packet);
        return 0;
    }
    result->oversized_status =
        (solver_codekd_test_reserve_status_t)
        solver_verification_packet_reserve(
            &packet,
            impossible_candidates);
    result->final_capacity = packet.candidate_capacity;
    result->allocation_failed = packet.allocation_failed;
    solver_verification_packet_free(&packet);
    return 0;
}
