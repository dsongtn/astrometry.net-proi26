/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* CodeKD descriptor generation and bounded payload delivery. */

#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "solver.h"
#include "verify.h"
#include "fitsbin.h"
#include "fit-wcs.h"
#include "log.h"
#include "mathutil.h"
#include "quad-utils.h"
#include "tic.h"
#include "index_shard_internal.h"
#include "solver_codekd_internal.h"
#include "solver_inline_internal.h"
#include "solver_hypothesis_internal.h"
#include "../libkd/kdtree_prefetch_internal.h"
/*
 * Return zero when the bounded output arena is ready, one for an exact native
 * fallback, and minus one for an invalid packet request.
 */
int solver_codekd_search_packet_prepare(
    const solver_t* solver,
    solver_codekd_search_packet_t* packet,
    solver_ab_descriptor_output_t* descriptors,
    const kdtree_t* tree,
    const quadfile_t* quads,
    startree_t* starkd,
    int dimquads,
    anbool use_radec,
    size_t candidate_budget_bytes,
    unsigned long long sequence,
    anbool detailed) {
    size_t metadata_bytes;
    size_t hit_bytes;

    if (!solver || !packet || !descriptors || !tree || !quads || !starkd ||
        !starkd->tree ||
        dimquads <= 0 || dimquads > DQMAX ||
        SOLVER_AB_DESCRIPTOR_CAPACITY >
            SIZE_MAX / sizeof(*packet->slots)) {
        return -1;
    }
#ifdef SOLVER_CODEKD_TEST_FORCE_PACKET_ALLOC_FAILURE
    if (solver->profile.ab_helper_tasks &&
        !solver->profile.allocation_failures) {
        return 1;
    }
#endif
    memset(packet, 0, sizeof(*packet));
    metadata_bytes =
        SOLVER_AB_DESCRIPTOR_CAPACITY * sizeof(*packet->slots);
    if (metadata_bytes >= SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES) {
        return 1;
    }
    hit_bytes =
        SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES - metadata_bytes;
    packet->hit_capacity =
        hit_bytes / (sizeof(*packet->inds) + sizeof(*packet->sdists));
    if (!packet->hit_capacity ||
        packet->hit_capacity > SIZE_MAX / sizeof(*packet->inds) ||
        packet->hit_capacity > SIZE_MAX / sizeof(*packet->sdists)) {
        return 1;
    }

    packet->slots = calloc(
        SOLVER_AB_DESCRIPTOR_CAPACITY, sizeof(*packet->slots));
    packet->inds = malloc(
        packet->hit_capacity * sizeof(*packet->inds));
    packet->sdists = malloc(
        packet->hit_capacity * sizeof(*packet->sdists));
    if (!packet->slots || !packet->inds || !packet->sdists ||
        solver_codekd_page_workspace_create(
            &packet->page_workspace)) {
        (void)solver_codekd_search_packet_cleanup(packet);
        return 1;
    }

    if (candidate_budget_bytes >=
        sizeof(*packet->candidate_records)) {
        packet->candidate_capacity = MIN(
            packet->hit_capacity,
            candidate_budget_bytes /
                sizeof(*packet->candidate_records));
        /*
         * candidate_records is a rolling scratch window, not storage for the
         * complete result stream. The owner retires this bounded prefix
         * before the same logical packet submits its next payload window.
         */
        packet->candidate_capacity = MIN(
            packet->candidate_capacity,
            (size_t)SOLVER_CANDIDATE_DELIVERY_LIMIT);
        if (packet->candidate_capacity >
            SIZE_MAX / sizeof(*packet->candidate_records)) {
            packet->candidate_capacity = 0U;
        }
        if (packet->candidate_capacity) {
            packet->candidate_records = calloc(
                packet->candidate_capacity,
                sizeof(*packet->candidate_records));
            if (!packet->candidate_records) {
                packet->candidate_capacity = 0U;
            } else {
                size_t record_bytes =
                    packet->candidate_capacity *
                        sizeof(*packet->candidate_records);

                packet->verify_query_budget =
                    candidate_budget_bytes > record_bytes
                    ? candidate_budget_bytes - record_bytes
                    : 0U;
            }
        }
    }

    packet->descriptors = descriptors;
    packet->tree = tree;
    packet->quads = quads;
    packet->starkd = starkd;
    packet->dimquads = dimquads;
    packet->use_radec = use_radec;
    packet->star_delivery_eligible =
        !starkd->tree->perm || starkd->inverse_perm != NULL;
    packet->candidate_star_delivery_disabled =
        !packet->star_delivery_eligible;
    packet->sequence = sequence;
    packet->detailed = detailed;
    packet->state = SOLVER_CODEKD_PACKET_ALLOCATED;
    return 0;
}

index_shard_helper_task_status_t
solver_codekd_packet_generate_descriptors(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    index_shard_helper_task_status_t descriptor_status;

    if (!input || input_size != sizeof(*input) ||
        !packet || output_size != sizeof(*packet) ||
        !input->tree || packet->tree != input->tree ||
        !input->quads || packet->quads != input->quads ||
        !input->starkd || packet->starkd != input->starkd ||
        !packet->descriptors || !packet->slots ||
        !packet->inds || !packet->sdists ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_ALLOCATED) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    descriptor_status = solver_ab_descriptor_helper_execute(
        &input->descriptor,
        sizeof(input->descriptor),
        packet->descriptors,
        sizeof(*packet->descriptors));
    if (descriptor_status != INDEX_SHARD_HELPER_TASK_OK) {
        packet->state =
            descriptor_status == INDEX_SHARD_HELPER_TASK_STOPPED
                ? SOLVER_CODEKD_PACKET_STOPPED
                : SOLVER_CODEKD_PACKET_FAILED;
        return descriptor_status;
    }
    if (packet->descriptors->descriptor_count >
        SOLVER_AB_DESCRIPTOR_CAPACITY) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    packet->first = 0U;
    packet->count = packet->descriptors->descriptor_count;
    packet->next_descriptor = 0U;
    packet->page_stats.descriptors_total = packet->count;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static int solver_codekd_packet_mark_plan_owner_replay(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason) {
    size_t descriptor_index;

    if (!packet || !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->count) {
        return -1;
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
            SOLVER_CODEKD_RESULT_UNUSED) {
            return -1;
        }
        packet->slots[descriptor_index].state =
            SOLVER_CODEKD_RESULT_OWNER_REPLAY;
    }
    solver_codekd_page_plan_record_refusal(packet, reason);
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return 0;
}

/*
 * A permanent pre-submission refusal invalidates the current sealed plan and
 * makes another plan attempt uneconomic. Preserve the exact logical work by
 * replaying the planned prefix and every still-unplanned descriptor at owner
 * retirement.
 */
static int solver_codekd_packet_mark_outstanding_owner_replay(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason) {
    size_t descriptor_index;

    if (!packet || !packet->plan_complete ||
        packet->state != SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->next_descriptor ||
        packet->next_descriptor > packet->count) {
        return -1;
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
            SOLVER_CODEKD_RESULT_UNUSED) {
            return -1;
        }
    }
    for (descriptor_index = packet->plan_end;
         descriptor_index < packet->count;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
                SOLVER_CODEKD_RESULT_UNUSED &&
            packet->slots[descriptor_index].state !=
                SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            return -1;
        }
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->count;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state ==
            SOLVER_CODEKD_RESULT_UNUSED) {
            packet->slots[descriptor_index].state =
                SOLVER_CODEKD_RESULT_OWNER_REPLAY;
        }
    }
    solver_codekd_page_plan_record_refusal(packet, reason);
    packet->next_descriptor = packet->count;
    packet->pending_descriptor_plan = FALSE;
    packet->pending_descriptor_raw_ranges = 0U;
    packet->pending_descriptor_logical_bytes = 0U;
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return 0;
}

/* Return one when submitted, zero on bounded capacity, and minus one when
 * the complete packet must use exact owner replay. */
int solver_codekd_packet_submit_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    int submit_status;

    if (!packet || !packet->tree || !packet->tree->io ||
        !packet->tree->io_is_fitsbin || !packet->page_workspace ||
        packet->delivery_source || packet->delivery_ticket ||
        packet->state != SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE ||
        !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->next_descriptor ||
        packet->next_descriptor > packet->count ||
        !packet->plan_range_count ||
        packet->plan_range_count >
            packet->page_workspace->sealed_range_capacity) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    source = (fitsbin_t*)packet->tree->io;
    packet->delivery_source = source;
    errno = 0;
    submit_status = fitsbin_prefetch_ranges_submit(
        source,
        packet->page_workspace->sealed_ranges,
        packet->plan_range_count,
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
        &packet->delivery_ticket);
    if (submit_status == FITSBIN_PAYLOAD_IO_SUBMIT_READY &&
        !packet->delivery_ticket) {
        packet->delivery_source = NULL;
        packet->state = SOLVER_CODEKD_PACKET_COMPUTE_READY;
        return 2;
    }
    if (submit_status > 0 && packet->delivery_ticket) {
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    if (solver_codekd_packet_mark_outstanding_owner_replay(
            packet,
            submit_status < 0
                ? SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR
                : SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
    }
    return -1;
}

int solver_codekd_packet_begin_candidate_window(
    solver_codekd_search_packet_t* packet) {
    size_t descriptor_index;
    anbool found = FALSE;

    if (!packet || !packet->candidate_records ||
        !packet->candidate_capacity || !packet->candidate_count ||
        packet->candidate_cursor >= packet->candidate_count ||
        packet->candidate_window_count ||
        packet->candidate_window_offset ||
        packet->verify_sweep_reads ||
        packet->verify_sweep_buffers ||
        packet->verify_sweep_storage ||
        packet->verify_sweep_storage_bytes) {
        return -1;
    }
    for (descriptor_index = packet->retire_descriptor;
         descriptor_index < packet->count;
         descriptor_index++) {
        const solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];
        size_t slot_end;

        if (slot->state != SOLVER_CODEKD_RESULT_READY ||
            !slot->hit_count ||
            slot->hit_first > packet->hit_count ||
            (size_t)slot->hit_count >
                packet->hit_count - slot->hit_first) {
            continue;
        }
        slot_end = slot->hit_first + (size_t)slot->hit_count;
        if (packet->candidate_cursor >= slot->hit_first &&
            packet->candidate_cursor < slot_end) {
            found = TRUE;
            break;
        }
    }
    if (!found) {
        return -1;
    }
    packet->candidate_window_first = packet->candidate_cursor;
    packet->candidate_window_count = MIN(
        MIN(packet->candidate_capacity,
            (size_t)SOLVER_CANDIDATE_DELIVERY_LIMIT),
        packet->candidate_count - packet->candidate_cursor);
    if (!packet->candidate_window_count) {
        return -1;
    }
    packet->candidate_window_offset = 0U;
    packet->candidate_star_count = 0U;
    packet->candidate_quad_ready_count = 0U;
    packet->candidate_star_ready_count = 0U;
    packet->candidate_verify_query_count = 0U;
    packet->verify_plan_first = 0U;
    packet->verify_plan_end = 0U;
    packet->verify_topology_end = 0U;
    packet->verify_plan_range_count = 0U;
    packet->verify_plan_logical_bytes = 0U;
    packet->verify_query_bytes = 0U;
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_storage_bytes = 0U;
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->verify_plan_complete = FALSE;
    packet->verify_sweep_plan_complete = FALSE;
    packet->verify_pending_query_plan = FALSE;
    packet->verification_delivery_disabled = FALSE;
    packet->candidate_delivery_windows++;
    packet->state = SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY;
    return 0;
}

void solver_codekd_packet_clear_candidate_window(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;

    if (!packet) {
        return;
    }
    assert(!packet->delivery_ticket);
    assert(!packet->delivery_source);
    solver_codekd_packet_clear_sweep_storage(packet);
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &packet->candidate_records[candidate_index]);
        verify_destroy_index_query(
            packet->candidate_records[
                candidate_index].verify_query);
        packet->candidate_records[
            candidate_index].verify_query = NULL;
        packet->candidate_records[
            candidate_index].verify_query_captured = FALSE;
    }
    packet->candidate_window_first = packet->candidate_cursor;
    packet->candidate_window_count = 0U;
    packet->candidate_window_offset = 0U;
    packet->candidate_star_count = 0U;
    packet->candidate_quad_ready_count = 0U;
    packet->candidate_star_ready_count = 0U;
    packet->candidate_verify_query_count = 0U;
    packet->verify_plan_first = 0U;
    packet->verify_plan_end = 0U;
    packet->verify_topology_end = 0U;
    packet->verify_plan_range_count = 0U;
    packet->verify_plan_logical_bytes = 0U;
    packet->verify_query_bytes = 0U;
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_storage_bytes = 0U;
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->verify_plan_complete = FALSE;
    packet->verify_sweep_plan_complete = FALSE;
    packet->verify_pending_query_plan = FALSE;
    packet->verification_delivery_disabled = FALSE;
}

/*
 * Build verification queries from already copied candidate data. The record
 * also retains the exact native gates and TAN result so ordered retirement
 * can reuse the calculation after validating the complete immutable input.
 */
static int solver_codekd_packet_build_verify_queries(
    solver_codekd_search_packet_t* packet,
    const solver_codekd_packet_task_input_t* input) {
    const solver_field_geometry_t* geometry;
    size_t candidate_index;

    if (!packet || !input || input->use_radec ||
        !packet->candidate_records ||
        !packet->candidate_window_count ||
        packet->candidate_star_ready_count !=
            packet->candidate_window_count ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return -1;
    }
    geometry = input->descriptor.field_geometry;
    if (!geometry || !geometry->fieldxy ||
        geometry->numxy <= 0) {
        return -1;
    }
    packet->candidate_verify_query_count = 0U;

    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        const solver_ab_descriptor_t* descriptor;
        double field_xy[DQMAX * 2];
        double corner[3];
        double scale;
        double arcsecperpix;
        double abscale;
        double cx;
        double cy;
        double radius;
        tan_t wcs;
        int star_index;

        record->plan_action = SOLVER_AB_CANDIDATE_SCALE_SKIP;
        record->candidate_prepared = FALSE;
        memset(record->verify_center, 0, sizeof(record->verify_center));
        record->verify_radius = 0.0;
        record->verify_radius2 = 0.0;
        if (record->descriptor_index >= packet->count) {
            return -1;
        }
        descriptor = &packet->descriptors->descriptors[
            record->descriptor_index];
        for (star_index = 0;
             star_index < packet->dimquads;
             star_index++) {
            int field_star = descriptor->stars[star_index];

            if (field_star < 0 ||
                field_star >= geometry->numxy) {
                return -1;
            }
            setx(
                field_xy,
                star_index,
                starxy_getx(geometry->fieldxy, field_star));
            sety(
                field_xy,
                star_index,
                starxy_gety(geometry->fieldxy, field_star));
            record->prepared_fieldstars[star_index] = field_star;
        }
        memcpy(
            record->prepared_fieldxy,
            field_xy,
            (size_t)packet->dimquads * 2U * sizeof(*field_xy));
        record->prepared_parity = descriptor->current_parity;

        abscale = square(distsq2rad(distsq(
            record->starxyz, record->starxyz + 3, 3))) /
            distsq(field_xy, field_xy + 2, 2);
        if (abscale > input->abscale_high ||
            abscale < input->abscale_low) {
            record->plan_action =
                SOLVER_AB_CANDIDATE_ABSCALE_SKIP;
            record->candidate_prepared = TRUE;
            packet->candidate_math_prepared++;
            continue;
        }
        if (fit_tan_wcs(
                record->starxyz,
                field_xy,
                packet->dimquads,
                &wcs,
                &scale)) {
            record->plan_action =
                SOLVER_AB_CANDIDATE_BAD_QUAD;
            record->candidate_prepared = TRUE;
            packet->candidate_math_prepared++;
            continue;
        }
        arcsecperpix = scale * 3600.0;
        memcpy(&record->prepared_wcs, &wcs, sizeof(wcs));
        record->prepared_scale = arcsecperpix;
        record->candidate_prepared = TRUE;
        packet->candidate_math_prepared++;
        if (arcsecperpix > input->funits_upper ||
            arcsecperpix < input->funits_lower) {
            continue;
        }

        record->plan_action = SOLVER_AB_CANDIDATE_VERIFY;
        cx = 0.5 * (input->field_minx + input->field_maxx);
        cy = 0.5 * (input->field_miny + input->field_maxy);
        tan_pixelxy2xyzarr(
            &wcs, cx, cy, record->verify_center);
        tan_pixelxy2xyzarr(
            &wcs, input->field_minx, input->field_miny, corner);
        radius = sqrt(distsq(
            record->verify_center, corner, 3));
        record->verify_radius = radius;
        record->verify_radius2 = square(radius);
        if (!isfinite(record->verify_radius2) ||
            record->verify_radius2 < 0.0) {
            record->verify_radius = 0.0;
            record->verify_radius2 = 0.0;
            record->verify_delivery_fallback = TRUE;
            continue;
        }
        packet->candidate_verify_query_count++;
        packet->verification_page_queries++;
    }
    return 0;
}

static anbool solver_codekd_packet_candidate_data_fully_resident(
    const solver_codekd_search_packet_t* packet) {
    fitsbin_t* star_source;

    if (!packet || !packet->quads || !packet->quads->fb ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin) {
        return FALSE;
    }
    star_source = (fitsbin_t*)packet->starkd->tree->io;
    return fitsbin_payload_is_fully_resident(packet->quads->fb) &&
        fitsbin_payload_is_fully_resident(star_source);
}

/*
 * Complete the CodeKD phase through one central transition. Resident
 * candidate payloads retain the existing verification-wave path. Otherwise,
 * the first bounded Quad/A/B window becomes independently deliverable.
 */
int solver_codekd_packet_finish_codekd(
    solver_codekd_search_packet_t* packet) {
    if (!packet ||
        (packet->state != SOLVER_CODEKD_PACKET_DESCRIPTORS_READY &&
         packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY) ||
        packet->next_descriptor != packet->count ||
        packet->plan_complete || packet->delivery_ticket ||
        packet->delivery_source) {
        return -1;
    }
    if (packet->candidate_count) {
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? 0
            : -1;
    }
    if (!packet->hit_count || !packet->candidate_records ||
        !packet->candidate_capacity || packet->use_radec ||
        !packet->quads || !packet->starkd ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX ||
        solver_codekd_packet_candidate_data_fully_resident(packet)) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }

    packet->candidate_count = packet->hit_count;
    packet->candidate_cursor = 0U;
    packet->retire_descriptor = 0U;
    packet->retire_hit_offset = 0U;
    packet->retire_descriptor_started = FALSE;
    return solver_codekd_packet_begin_candidate_window(packet);
}

/*
 * Submit one bounded post-CodeKD payload stage. A permanent refusal preserves
 * the copied CodeKD results and falls back only candidate resolution.
 */
int solver_codekd_packet_submit_candidate_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    int submit_status;

    if (!packet || packet->delivery_source || packet->delivery_ticket ||
        !packet->candidate_records || !packet->candidate_count ||
        packet->candidate_cursor >= packet->candidate_count ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }

    errno = 0;
    if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY) {
        if (packet->candidate_quad_delivery_disabled) {
            packet->candidate_quad_fallback++;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return -1;
        }
        if (!packet->quads || !packet->quads->fb) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return -1;
        }
        source = packet->quads->fb;
        submit_status = quadfile_prefetch_stars_submit(
            packet->quads,
            packet->inds + packet->candidate_window_first,
            (int)packet->candidate_window_count,
            &packet->delivery_ticket);
        if (submit_status == FITSBIN_PAYLOAD_IO_SUBMIT_READY &&
            !packet->delivery_ticket) {
            packet->state = SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
            return 2;
        }
        if (submit_status > 0 && packet->delivery_ticket) {
            packet->delivery_source = source;
            packet->candidate_quad_submitted++;
            packet->state = SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED;
            return 1;
        }
        packet->delivery_ticket = NULL;
        if (!submit_status && !errno &&
            fitsbin_payload_is_fully_resident(source)) {
            packet->state = SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
            return 2;
        }
        if (!submit_status && errno == EAGAIN) {
            return 0;
        }
        packet->candidate_quad_fallback++;
        packet->candidate_quad_delivery_disabled = TRUE;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return -1;
    }

    if (packet->state != SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY ||
        packet->candidate_star_delivery_disabled ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->candidate_star_count ||
        packet->candidate_star_count > SOLVER_CANDIDATE_STAR_LIMIT) {
        packet->candidate_star_fallback++;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return -1;
    }
    source = (fitsbin_t*)packet->starkd->tree->io;
    submit_status = startree_prefetch_stars_ready_submit(
        packet->starkd,
        packet->candidate_starids,
        (int)packet->candidate_star_count,
        &packet->delivery_ticket);
    if (submit_status == FITSBIN_PAYLOAD_IO_SUBMIT_READY &&
        !packet->delivery_ticket) {
        packet->state = SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY;
        return 2;
    }
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->delivery_source = source;
        packet->candidate_star_submitted++;
        packet->state = SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED;
        return 1;
    }
    packet->delivery_ticket = NULL;
    if (!submit_status && !errno &&
        fitsbin_payload_is_fully_resident(source)) {
        packet->state = SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY;
        return 2;
    }
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    packet->candidate_star_fallback++;
    packet->candidate_star_delivery_disabled = TRUE;
    packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
    return -1;
}

int solver_codekd_packet_submit_verification_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    size_t candidate_index;
    anbool has_query = FALSE;
    int submit_status;

    if (!packet || packet->delivery_source ||
        packet->delivery_ticket ||
        !packet->candidate_records ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY ||
        packet->verification_delivery_disabled ||
        packet->verify_plan_complete ||
        packet->candidate_window_offset >=
            packet->candidate_window_count ||
        packet->candidate_star_ready_count !=
            packet->candidate_window_count ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }
    for (candidate_index = packet->candidate_window_offset;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        const solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];

        if (record->verify_delivery_fallback) {
            solver_codekd_packet_disable_verification_delivery(
                packet);
            return -1;
        }
        if (record->plan_action !=
            SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        has_query = TRUE;
        break;
    }
    if (!has_query) {
        packet->verify_plan_first =
            packet->candidate_window_offset;
        packet->verify_plan_end =
            packet->candidate_window_count;
        packet->verify_topology_end =
            packet->candidate_window_count;
        packet->verify_plan_complete = TRUE;
        packet->verification_page_prefixes++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
        return 2;
    }

    source = (fitsbin_t*)packet->starkd->tree->io;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->delivery_source = source;
    packet->state = SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED;
    errno = 0;
    submit_status = fitsbin_prefetch_ranges_planned_submit(
        source,
        solver_codekd_packet_plan_verification_pages,
        packet,
        packet->verify_delivery_budget,
        &packet->delivery_ticket);
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->verification_page_submitted++;
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && !errno &&
        fitsbin_payload_is_fully_resident(source)) {
        packet->verify_plan_first =
            packet->candidate_window_offset;
        packet->verify_plan_end =
            packet->candidate_window_count;
        packet->verify_topology_end =
            packet->candidate_window_count;
        packet->verify_plan_complete = TRUE;
        packet->verification_page_prefixes++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
        return 2;
    }
    if (!submit_status && errno == EAGAIN) {
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY;
        return 0;
    }
    solver_codekd_packet_disable_verification_delivery(packet);
    return -1;
}

int solver_codekd_packet_submit_sweep_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
#if SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED
    size_t offset;
    size_t range_index;
    anbool direct = FALSE;
#endif
    int submit_status;

    if (!packet || packet->delivery_source ||
        packet->delivery_ticket || !packet->starkd ||
        !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->page_workspace ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY ||
        !packet->verify_plan_complete ||
        !packet->verify_sweep_plan_complete ||
        !packet->verify_sweep_range_count ||
        packet->verify_sweep_range_count >
            packet->page_workspace->sealed_range_capacity ||
        !packet->verify_sweep_aligned_bytes ||
        packet->verify_sweep_aligned_bytes >
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }

    source = (fitsbin_t*)packet->starkd->tree->io;
    if (fitsbin_payload_is_fully_resident(source)) {
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
        return 2;
    }

#if SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED
    if (packet->verify_sweep_reads ||
        packet->verify_sweep_buffers ||
        packet->verify_sweep_storage ||
        packet->verify_sweep_storage_bytes) {
        if (!packet->verify_sweep_reads ||
            !packet->verify_sweep_buffers ||
            !packet->verify_sweep_storage ||
            packet->verify_sweep_storage_bytes !=
                packet->verify_sweep_aligned_bytes) {
            solver_codekd_packet_clear_sweep_storage(packet);
        } else {
            direct = TRUE;
        }
    }
    if (!direct &&
        packet->verify_sweep_range_count <=
            FITSBIN_PREAD_ASYNC_RANGE_LIMIT &&
        packet->verify_query_bytes <=
            packet->verify_query_budget &&
        packet->verify_sweep_aligned_bytes <=
            packet->verify_query_budget -
                packet->verify_query_bytes) {
        if (packet->verify_sweep_range_count <=
                SIZE_MAX / sizeof(*packet->verify_sweep_reads) &&
            packet->verify_sweep_range_count <=
                SIZE_MAX / sizeof(*packet->verify_sweep_buffers)) {
            packet->verify_sweep_reads = calloc(
                packet->verify_sweep_range_count,
                sizeof(*packet->verify_sweep_reads));
            packet->verify_sweep_buffers = calloc(
                packet->verify_sweep_range_count,
                sizeof(*packet->verify_sweep_buffers));
            packet->verify_sweep_storage = malloc(
                packet->verify_sweep_aligned_bytes);
        }
        if (packet->verify_sweep_reads &&
            packet->verify_sweep_buffers &&
            packet->verify_sweep_storage) {
            offset = 0U;
            direct = TRUE;
            for (range_index = 0U;
                 range_index < packet->verify_sweep_range_count;
                 range_index++) {
                const fitsbin_prefetch_range_t* range =
                    &packet->page_workspace->
                        sealed_ranges[range_index];

                if (!range->data || !range->size ||
                    range->size >
                        packet->verify_sweep_aligned_bytes -
                            offset) {
                    direct = FALSE;
                    break;
                }
                packet->verify_sweep_reads[range_index].data =
                    range->data;
                packet->verify_sweep_reads[range_index].size =
                    range->size;
                packet->verify_sweep_reads[
                    range_index].logical_size = range->size;
                packet->verify_sweep_reads[
                    range_index].destination =
                        packet->verify_sweep_storage + offset;
                packet->verify_sweep_buffers[
                    range_index].mapping_data = range->data;
                packet->verify_sweep_buffers[
                    range_index].size = range->size;
                packet->verify_sweep_buffers[
                    range_index].bytes =
                        packet->verify_sweep_storage + offset;
                offset += range->size;
            }
            if (offset != packet->verify_sweep_aligned_bytes) {
                direct = FALSE;
            }
        }
        if (!direct) {
            solver_codekd_packet_clear_sweep_storage(packet);
        } else {
            packet->verify_sweep_storage_bytes =
                packet->verify_sweep_aligned_bytes;
        }
    }
#else
    /*
     * A mapped-only transport comparison must bypass the direct scratch path
     * before its range arrays and page-sized destination buffer are allocated.
     * Refusing the later fitsbin submit is not equivalent: it retains the
     * allocation, initialization, and free traffic for every sweep ticket.
     */
    solver_codekd_packet_clear_sweep_storage(packet);
#endif

    errno = 0;
#if SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED
    if (direct) {
        submit_status = fitsbin_pread_mapped_ranges_submit(
            source,
            packet->verify_sweep_reads,
            packet->verify_sweep_range_count,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
            &packet->delivery_ticket);
        if (submit_status > 0 && packet->delivery_ticket) {
            packet->delivery_source = source;
            packet->verification_page_submitted++;
            packet->state =
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED;
            return 1;
        }
        packet->delivery_ticket = NULL;
        if (!submit_status && errno == EAGAIN) {
            return 0;
        }
        /*
         * A direct translator or allocation refusal is not a reason to
         * discard the already complete mapped-page plan. Preserve the
         * previous completion provider before falling back to native mmap.
         */
        solver_codekd_packet_clear_sweep_storage(packet);
        errno = 0;
    }
#endif
    submit_status = fitsbin_prefetch_ranges_submit(
        source,
        packet->page_workspace->sealed_ranges,
        packet->verify_sweep_range_count,
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
        &packet->delivery_ticket);
    if (submit_status == FITSBIN_PAYLOAD_IO_SUBMIT_READY &&
        !packet->delivery_ticket) {
        packet->delivery_source = NULL;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
        return 2;
    }
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->delivery_source = source;
        packet->verification_page_submitted++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED;
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    solver_codekd_packet_clear_sweep_storage(packet);
    packet->verification_page_fallback++;
    packet->state = SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
    return 2;
}

/* Return zero while pending, one when compute-ready, two when replanned work
 * remains, three when the logical packet is complete, and minus one on error
 * or cooperative stop. */
int solver_codekd_packet_collect_pages(
    solver_codekd_search_packet_t* packet) {
    solver_codekd_search_packet_state_t io_state;
    int ticket_result = 0;
    int ticket_errno;
    int poll_status;

    if (!packet || !packet->delivery_source ||
        !packet->delivery_ticket) {
        return -1;
    }
    errno = 0;
    poll_status = fitsbin_payload_io_ticket_poll_and_destroy(
        packet->delivery_source,
        &packet->delivery_ticket,
        &ticket_result);
    ticket_errno = errno;
    if (!poll_status) {
        return 0;
    }
    if (poll_status < 0) {
        int drain_status =
            fitsbin_payload_io_ticket_drain_and_destroy(
                packet->delivery_source,
                &packet->delivery_ticket,
                &ticket_result);

        if (drain_status == 1 && !packet->delivery_ticket) {
            packet->delivery_source = NULL;
            solver_codekd_packet_clear_sweep_storage(packet);
        } else {
            logerr("[solver] payload ticket ownership drain "
                   "failed errno=%i state=%i\n",
                   errno,
                   (int)packet->state);
        }
        logerr("[solver] payload ticket ownership transfer "
               "failed errno=%i state=%i\n",
               ticket_errno,
               (int)packet->state);
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return -1;
    }
    /*
     * Verification planned-ticket output is owned by the I/O lane until
     * terminal publication. Poll's payload_io mutex edge makes those state
     * writes visible here. Initial CodeKD tickets use an immutable sealed
     * plan, so their packet remains PAGE_PLAN_COMPLETE while queued.
     */
    io_state = packet->state;
    assert(!packet->delivery_ticket);
    if ((ticket_result == 0 && ticket_errno == ECANCELED) ||
        index_shard_worker_stop_requested()) {
        packet->delivery_source = NULL;
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    if (ticket_result > 0 &&
        io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        packet->state = SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state == SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        packet->state = SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE &&
        packet->verify_plan_complete &&
        packet->verify_plan_first ==
            packet->candidate_window_offset &&
        packet->verify_plan_end >
            packet->verify_plan_first &&
            packet->verify_plan_end <=
            packet->candidate_window_count) {
        packet->delivery_source = NULL;
        packet->verification_page_ready++;
        packet->verification_page_ready_rows +=
            packet->verify_plan_end -
                packet->verify_plan_first;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED &&
        packet->verify_plan_complete &&
        packet->verify_sweep_plan_complete &&
        packet->verify_sweep_range_count) {
        packet->delivery_source = NULL;
        packet->verification_page_ready++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->verification_page_fallback++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 &&
        (io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED ||
         io_state == SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED)) {
        packet->delivery_source = NULL;
        if (io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
        } else {
            packet->candidate_star_fallback++;
            packet->candidate_star_delivery_disabled = TRUE;
        }
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 3;
    }
    if (io_state == SOLVER_CODEKD_PACKET_RESULTS_READY &&
        packet->verification_delivery_disabled) {
        packet->delivery_source = NULL;
        return 3;
    }
    if (ticket_result < 0 &&
        (io_state == SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED ||
         io_state ==
             SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE)) {
        packet->delivery_source = NULL;
        if (packet->state == SOLVER_CODEKD_PACKET_FAILED) {
            return -1;
        }
        solver_codekd_packet_disable_verification_delivery(packet);
        return 3;
    }
    if (ticket_result > 0 &&
        packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY &&
        !packet->plan_complete) {
        packet->delivery_source = NULL;
        return 3;
    }
    if (ticket_result > 0 && packet->plan_complete &&
        packet->plan_range_count &&
        packet->state == SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
        packet->state = SOLVER_CODEKD_PACKET_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 && packet->plan_complete &&
        packet->state == SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
        if (solver_codekd_packet_mark_plan_owner_replay(
                packet, SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return -1;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? 3
            : 2;
    }
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_FAILED;
    return -1;
}

int solver_codekd_packet_cancel_pages(
    solver_codekd_search_packet_t* packet) {
    if (!packet || !packet->delivery_ticket) {
        return 0;
    }
    return fitsbin_payload_io_ticket_cancel_async(
        packet->delivery_ticket);
}

index_shard_helper_task_status_t
solver_codekd_packet_decode_quad_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    size_t stars_per_candidate;
    size_t candidate_index;
    size_t descriptor_index;

    if (!packet || !more_work || !packet->quads || !packet->starkd ||
        !packet->candidate_records ||
        packet->state != SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY ||
        !packet->candidate_count ||
        !packet->candidate_capacity ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_first >
            packet->candidate_count ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        packet->candidate_quad_ready_count ||
        packet->candidate_star_ready_count ||
        packet->use_radec ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }

    /*
     * Phase admission already identified a dense result stream. Populate
     * every quad-star row in one bounded ticket, but preserve the owner's
     * native A/B scale gate and logical consumption order.
     */
    stars_per_candidate = (size_t)packet->dimquads;
    if (packet->candidate_window_count >
        SOLVER_CANDIDATE_STAR_LIMIT / stars_per_candidate) {
        packet->candidate_quad_fallback++;
        packet->candidate_quad_delivery_disabled = TRUE;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->candidate_star_count =
        packet->candidate_window_count * stars_per_candidate;
    descriptor_index = packet->retire_descriptor;
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        size_t global_candidate =
            packet->candidate_window_first + candidate_index;
        size_t star_index;

        while (descriptor_index < packet->count) {
            const solver_codekd_result_slot_t* slot =
                &packet->slots[descriptor_index];
            size_t slot_end;

            if (slot->state != SOLVER_CODEKD_RESULT_READY ||
                !slot->hit_count) {
                descriptor_index++;
                continue;
            }
            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first) {
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            slot_end = slot->hit_first + (size_t)slot->hit_count;
            if (global_candidate < slot->hit_first) {
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            if (global_candidate < slot_end) {
                break;
            }
            descriptor_index++;
        }
        if (descriptor_index >= packet->count) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }

        memset(record, 0, sizeof(*record));
        record->descriptor_index = descriptor_index;
        record->quadid = packet->inds[global_candidate];
        if (quadfile_get_stars(
                packet->quads,
                record->quadid,
                record->stars)) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
            packet->candidate_star_count = 0U;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return INDEX_SHARD_HELPER_TASK_OK;
        }
        for (star_index = 0U;
             star_index < stars_per_candidate;
             star_index++) {
            packet->candidate_starids[
                candidate_index * stars_per_candidate + star_index] =
                record->stars[star_index];
        }
    }
    packet->candidate_quad_ready_count =
        packet->candidate_window_count;
    packet->candidate_quad_ready_rows +=
        packet->candidate_window_count;
    packet->candidate_quad_ready++;
    if (!packet->star_delivery_eligible ||
        packet->candidate_star_delivery_disabled) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->state = SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY;
    *more_work = TRUE;
    return INDEX_SHARD_HELPER_TASK_OK;
}

index_shard_helper_task_status_t
solver_codekd_packet_copy_star_ready(
    const solver_codekd_packet_task_input_t* input,
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    size_t stars_per_candidate;
    size_t candidate_index;

    if (!input || !packet || !more_work || !packet->starkd ||
        packet->state != SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY ||
        !packet->candidate_records || !packet->candidate_count ||
        !packet->candidate_capacity ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_first >
            packet->candidate_count ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        packet->candidate_quad_ready_count !=
            packet->candidate_window_count ||
        packet->candidate_star_ready_count ||
        !packet->candidate_star_count ||
        packet->candidate_star_count > SOLVER_CANDIDATE_STAR_LIMIT ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }

    stars_per_candidate = (size_t)packet->dimquads;
    if (packet->candidate_star_count !=
        packet->candidate_window_count * stars_per_candidate) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        size_t star_index;

        for (star_index = 0U;
             star_index < stars_per_candidate;
             star_index++) {
            if (startree_get_ready(
                    packet->starkd,
                    (int)record->stars[star_index],
                    record->starxyz + 3U * star_index)) {
                packet->candidate_star_fallback++;
                packet->candidate_star_delivery_disabled = TRUE;
                packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
                return INDEX_SHARD_HELPER_TASK_OK;
            }
        }
    }
    packet->candidate_star_ready_count =
        packet->candidate_window_count;
    packet->candidate_star_ready_rows +=
        packet->candidate_window_count;
    packet->candidate_star_ready++;
    if (solver_codekd_packet_build_verify_queries(
            packet, input)) {
        for (candidate_index = 0U;
             candidate_index < packet->candidate_window_count;
             candidate_index++) {
            packet->candidate_records[
                candidate_index].candidate_prepared = FALSE;
        }
        packet->verification_delivery_disabled = TRUE;
        packet->verification_page_fallback++;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->state = packet->candidate_verify_query_count
        ? SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY
        : SOLVER_CODEKD_PACKET_RESULTS_READY;
    if (packet->state == SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY) {
        *more_work = TRUE;
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}
