/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* Staged packet callbacks and exact owner execution. */

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "solver.h"
#include "verify.h"
#include "fitsbin.h"
#include "log.h"
#include "tic.h"
#include "index_shard_internal.h"
#include "solver_codekd_internal.h"
#include "solver_inline_internal.h"
#include "solver_hypothesis_internal.h"
/*
 * Execute only a fully populated, contiguous descriptor slice. Native CodeKD
 * remains the authoritative query implementation. Page eviction between I/O
 * completion and execution is harmless; it affects performance, not results.
 */
static index_shard_helper_task_status_t
solver_codekd_packet_execute_ready(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size,
    anbool* more_work) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    kdtree_qres_t* query_result = NULL;
    size_t descriptor_index;
    anbool query_failed = FALSE;

    if (more_work) {
        *more_work = FALSE;
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY) {
        return solver_codekd_packet_decode_quad_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY) {
        return solver_codekd_packet_copy_star_ready(
            input, packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY) {
        return solver_codekd_packet_query_verification_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY) {
        return solver_codekd_packet_capture_sweep_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY) {
        return solver_codekd_packet_score_verification_ready(packet);
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY) {
        if (!input || input_size != sizeof(*input) ||
            output_size != sizeof(*packet) ||
            !packet->verify_plan_complete ||
            packet->verify_plan_first !=
                packet->candidate_window_offset ||
            packet->verify_plan_end <= packet->verify_plan_first ||
            packet->verify_plan_end >
                packet->candidate_window_count) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    if (!input || input_size != sizeof(*input) ||
        !packet || output_size != sizeof(*packet) ||
        !input->tree || packet->tree != input->tree ||
        !packet->descriptors || !packet->slots ||
        !packet->inds || !packet->sdists ||
        !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->count ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_COMPUTE_READY) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    packet->state = SOLVER_CODEKD_PACKET_EXECUTING;
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[descriptor_index];
        solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];
        kdtree_qres_t* previous_result = query_result;
        double search_wall_start = 0.0;

        if (index_shard_worker_stop_requested()) {
            kdtree_free_query(query_result);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        if (packet->detailed) {
            search_wall_start = monotonic_seconds();
        }
        query_result = solver_codekd_rangesearch(
            input->tree,
            previous_result,
            descriptor->code,
            descriptor->tol2,
            SOLVER_CODEKD_SEARCH_OPTIONS);
        slot->search_errno = errno;
        if (packet->detailed) {
            slot->search_wall_seconds =
                monotonic_seconds() - search_wall_start;
        }
        if (!query_result) {
            kdtree_free_query(previous_result);
            query_result = NULL;
            slot->state = SOLVER_CODEKD_RESULT_QUERY_FAILED;
            query_failed = TRUE;
            break;
        }
        if (query_result->nres &&
            (!query_result->inds || !query_result->sdists)) {
            slot->search_errno = EIO;
            slot->state = SOLVER_CODEKD_RESULT_QUERY_FAILED;
            query_failed = TRUE;
            break;
        }
        if (packet->hit_count > packet->hit_capacity ||
            (size_t)query_result->nres >
                packet->hit_capacity - packet->hit_count) {
            slot->state = SOLVER_CODEKD_RESULT_OWNER_REPLAY;
            continue;
        }

        slot->hit_first = packet->hit_count;
        slot->hit_count = query_result->nres;
        if (slot->hit_count) {
            memcpy(
                packet->inds + slot->hit_first,
                query_result->inds,
                (size_t)slot->hit_count * sizeof(*packet->inds));
            memcpy(
                packet->sdists + slot->hit_first,
                query_result->sdists,
                (size_t)slot->hit_count * sizeof(*packet->sdists));
        }
        packet->hit_count += slot->hit_count;
        slot->state = SOLVER_CODEKD_RESULT_READY;
    }

    kdtree_free_query(query_result);
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    if (query_failed) {
        packet->next_descriptor = packet->count;
        packet->pending_descriptor_plan = FALSE;
        packet->pending_descriptor_raw_ranges = 0U;
        packet->pending_descriptor_logical_bytes = 0U;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        if (packet->hit_count) {
            if (solver_codekd_packet_finish_codekd(packet)) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            if (packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY &&
                more_work) {
                *more_work = TRUE;
            }
        }
    } else if (packet->next_descriptor < packet->count) {
        packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
        if (more_work) {
            *more_work = TRUE;
        }
    } else {
        packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
        if (more_work) {
            *more_work = TRUE;
        }
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

static index_shard_staged_prepare_status_t
solver_codekd_packet_staged_prepare(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;
    index_shard_helper_task_status_t descriptor_status;
    int plan_status;

    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_PREPARE_ERROR;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_ALLOCATED) {
        descriptor_status = solver_codekd_packet_generate_descriptors(
            input_bytes,
            input_size,
            output_bytes,
            output_size);
        if (descriptor_status == INDEX_SHARD_HELPER_TASK_STOPPED) {
            return INDEX_SHARD_STAGED_PREPARE_STOPPED;
        }
        if (descriptor_status != INDEX_SHARD_HELPER_TASK_OK) {
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
    }
    if (packet->state != SOLVER_CODEKD_PACKET_DESCRIPTORS_READY) {
        if (packet->state ==
                SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
            return packet->plan_complete &&
                    packet->plan_first < packet->plan_end &&
                    packet->plan_end <= packet->next_descriptor &&
                    packet->next_descriptor <= packet->count &&
                    packet->plan_range_count &&
                    packet->plan_range_count <=
                        packet->page_workspace->sealed_range_capacity
                ? INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY
                : INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER) {
            return INDEX_SHARD_STAGED_PREPARE_OWNER_READY;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY) {
            return INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY) {
            if (!packet->verify_plan_complete ||
                packet->verify_plan_first !=
                    packet->candidate_window_offset ||
                packet->verify_plan_end <=
                    packet->verify_plan_first ||
                packet->verify_plan_end !=
                    packet->verify_topology_end ||
                packet->verify_topology_end >
                    packet->candidate_window_count) {
                logerr("[solver] invalid verification query packet "
                       "state=%i first=%zu end=%zu topology_end=%zu "
                       "window_offset=%zu window_count=%zu\n",
                       (int)packet->state,
                       packet->verify_plan_first,
                       packet->verify_plan_end,
                       packet->verify_topology_end,
                       packet->candidate_window_offset,
                       packet->candidate_window_count);
                return INDEX_SHARD_STAGED_PREPARE_ERROR;
            }
            return INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY;
        }
        if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY ||
            packet->state == SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY ||
            packet->state == SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY) {
            return INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
        }
        if (packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY) {
            logerr("[solver] unexpected packet state during prepare "
                   "state=%i sequence=%llu\n",
                   (int)packet->state,
                   packet->sequence);
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (solver_codekd_packet_finish_codekd(packet)) {
            logerr("[solver] packet finish failed during prepare "
                   "state=%i sequence=%llu retire=%zu/%zu "
                   "verify=%zu/%zu prepared=%zu\n",
                   (int)packet->state,
                   packet->sequence,
                   packet->retire_descriptor,
                   packet->count,
                   packet->verify_plan_first,
                   packet->verify_plan_end,
                   packet->verify_prepared_count);
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_PREPARE_RESULTS_READY
            : INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
    }
    if (packet->next_descriptor < packet->count) {
        plan_status = solver_codekd_search_packet_prepare_next_plan(
            packet,
            solver_codekd_packet_plan_cancelled,
            NULL);
        if (plan_status < 0) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (plan_status == 2) {
            return INDEX_SHARD_STAGED_PREPARE_STOPPED;
        }
        if (plan_status > 0) {
            return packet->state ==
                        SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE &&
                    packet->plan_complete &&
                    packet->plan_first < packet->plan_end &&
                    packet->plan_end <= packet->next_descriptor &&
                    packet->next_descriptor <= packet->count &&
                    packet->plan_range_count &&
                    packet->plan_range_count <=
                        packet->page_workspace->sealed_range_capacity
                ? INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY
                : INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
    }
    if (packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet->next_descriptor == packet->count) {
        if (solver_codekd_packet_finish_codekd(packet)) {
            logerr("[solver] packet finish failed after planning "
                   "state=%i sequence=%llu retire=%zu/%zu "
                   "verify=%zu/%zu prepared=%zu\n",
                   (int)packet->state,
                   packet->sequence,
                   packet->retire_descriptor,
                   packet->count,
                   packet->verify_plan_first,
                   packet->verify_plan_end,
                   packet->verify_prepared_count);
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_PREPARE_RESULTS_READY
            : INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
    }
    logerr("[solver] invalid terminal packet prepare state=%i "
           "sequence=%llu next=%zu count=%zu\n",
           (int)packet->state,
           packet->sequence,
           packet->next_descriptor,
           packet->count);
    packet->state = SOLVER_CODEKD_PACKET_FAILED;
    return INDEX_SHARD_STAGED_PREPARE_ERROR;
}

static index_shard_staged_submit_status_t
solver_codekd_packet_staged_submit(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size,
    unsigned long long* completion_id_out) {
    solver_codekd_search_packet_t* packet = output_bytes;
    int status;

    (void)input_bytes;
    (void)input_size;
    if (completion_id_out) {
        *completion_id_out = 0ULL;
    }
    if (!packet || output_size != sizeof(*packet) ||
        !completion_id_out) {
        return INDEX_SHARD_STAGED_SUBMIT_ERROR;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY ||
        packet->state == SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY) {
        status = solver_codekd_packet_submit_candidate_pages(packet);
    } else if (packet->state ==
               SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY) {
        status =
            solver_codekd_packet_submit_verification_pages(packet);
    } else if (packet->state ==
               SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY) {
        status = solver_codekd_packet_submit_sweep_pages(packet);
    } else {
        status = solver_codekd_packet_submit_pages(packet);
    }
    if (status == 2 &&
        (packet->state == SOLVER_CODEKD_PACKET_COMPUTE_READY ||
         packet->state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY ||
         packet->state == SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY ||
         packet->state ==
             SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY ||
         packet->state ==
             SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY ||
         packet->state == SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY)) {
        return INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY;
    }
    if (status > 0) {
        *completion_id_out =
            fitsbin_payload_io_ticket_completion_id(
                packet->delivery_ticket);
        if (!*completion_id_out) {
            (void)solver_codekd_search_packet_release_ticket(packet);
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_SUBMIT_ERROR;
        }
        return INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED;
    }
    if (!status) {
        return INDEX_SHARD_STAGED_SUBMIT_RETRY;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_STOPPED) {
        return INDEX_SHARD_STAGED_SUBMIT_STOPPED;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
        packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY) {
        return INDEX_SHARD_STAGED_SUBMIT_OWNER_READY;
    }
    return INDEX_SHARD_STAGED_SUBMIT_ERROR;
}

static index_shard_staged_io_status_t
solver_codekd_packet_staged_poll(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;
    int status;

    (void)input_bytes;
    (void)input_size;
    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_IO_ERROR;
    }
    status = solver_codekd_packet_collect_pages(packet);
    if (!status) {
        return INDEX_SHARD_STAGED_IO_PENDING;
    }
    if (status == 1) {
        return INDEX_SHARD_STAGED_IO_READY;
    }
    if (status == 2 || status == 3) {
        return INDEX_SHARD_STAGED_IO_FAILED;
    }
    if (packet->delivery_ticket || packet->delivery_source) {
        logerr("[solver] staged payload poll returned terminal "
               "without releasing ticket ownership state=%i\n",
               (int)packet->state);
    }
    return packet->state == SOLVER_CODEKD_PACKET_STOPPED
        ? INDEX_SHARD_STAGED_IO_CANCELLED
        : INDEX_SHARD_STAGED_IO_ERROR;
}

static int solver_codekd_packet_staged_cancel(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;

    (void)input_bytes;
    (void)input_size;
    if (!packet || output_size != sizeof(*packet)) {
        return -1;
    }
    return solver_codekd_packet_cancel_pages(packet) < 0
        ? -1
        : 0;
}

static index_shard_staged_execute_status_t
solver_codekd_packet_staged_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    anbool more_work = FALSE;
    index_shard_helper_task_status_t status =
        solver_codekd_packet_execute_ready(
            input_bytes,
            input_size,
            output_bytes,
            output_size,
            &more_work);

    if (status == INDEX_SHARD_HELPER_TASK_STOPPED) {
        return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
    }
    if (status != INDEX_SHARD_HELPER_TASK_OK) {
        return INDEX_SHARD_STAGED_EXECUTE_ERROR;
    }
    return more_work
        ? INDEX_SHARD_STAGED_EXECUTE_MORE
        : INDEX_SHARD_STAGED_EXECUTE_OK;
}

static index_shard_staged_execute_status_t
solver_codekd_packet_staged_owner(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    int prepare_status;

    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_EXECUTE_ERROR;
    }
    if (index_shard_worker_stop_requested() ||
        packet->state == SOLVER_CODEKD_PACKET_STOPPED) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_DESCRIPTORS_READY) {
        return INDEX_SHARD_STAGED_EXECUTE_MORE;
    }
    if (packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER) {
        if (!input || input_size != sizeof(*input)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        prepare_status =
            solver_codekd_packet_prepare_verification_owner(
                input, packet);
        if (prepare_status == 2) {
            return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
        }
        if (prepare_status < 0) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        return prepare_status
            ? INDEX_SHARD_STAGED_EXECUTE_MORE
            : INDEX_SHARD_STAGED_EXECUTE_OK;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY) {
        if (solver_codekd_packet_finish_codekd(packet)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_EXECUTE_OK
            : INDEX_SHARD_STAGED_EXECUTE_MORE;
    }
    return INDEX_SHARD_STAGED_EXECUTE_ERROR;
}

const index_shard_staged_ops_t solver_codekd_packet_staged_ops = {
    "codekd-packet",
    solver_codekd_packet_staged_prepare,
    solver_codekd_packet_staged_submit,
    solver_codekd_packet_staged_poll,
    solver_codekd_packet_staged_cancel,
    solver_codekd_packet_staged_execute,
    solver_codekd_packet_staged_owner,
    TRUE
};

/* Exact native execution used whenever packet admission is unavailable. */
int solver_codekd_descriptor_execute_owner(
    const solver_ab_descriptor_output_t* output,
    solver_t* solver,
    int dimquads,
    kdtree_qres_t** query_result,
    unsigned long long* reduced) {
    size_t descriptor_index;

    if (!output || !solver || !query_result || !reduced) {
        return -1;
    }
    solver->profile.max_batch_hypotheses = MAX(
        solver->profile.max_batch_hypotheses,
        output->descriptor_count);
    for (descriptor_index = 0U;
         descriptor_index < output->descriptor_count;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &output->descriptors[descriptor_index];

        if (solver_poll_worker_stop(solver)) {
            return 1;
        }
        solver->rel_field_noise2 = descriptor->rel_field_noise2;
        if (solver_ab_checked_counter_delta(
                solver,
                descriptor->numtries_delta,
                descriptor->cxdx_delta,
                descriptor->meanx_delta)) {
            return -1;
        }
        solver_execute_hypothesis_owner(
            descriptor->stars,
            descriptor->code,
            dimquads,
            solver,
            descriptor->current_parity,
            descriptor->tol2,
            query_result);
        (*reduced)++;
        if (solver->profile.execution_failed) {
            return -1;
        }
        if (solver->quit_now) {
            return 1;
        }
    }
    if (output->has_final_rel_field_noise2) {
        solver->rel_field_noise2 = output->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            solver,
            output->trailing_numtries,
            output->trailing_cxdx,
            output->trailing_meanx)) {
        return -1;
    }
    return 0;
}

void solver_codekd_packet_begin_result_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    const solver_codekd_result_slot_t* slot,
    const kdtree_qres_t* result,
    int dimquad) {
    if (!solver || !descriptor || !slot || !result ||
        dimquad <= 0 || dimquad > DQMAX) {
        return;
    }
    solver_begin_hypothesis_owner(
        descriptor->stars,
        descriptor->code,
        dimquad,
        solver,
        descriptor->current_parity);
    solver->profile.codekd_calls++;
    solver->profile.hypotheses_executed++;
    if (solver->profile.detailed) {
        solver->profile.codekd_wall_seconds +=
            slot->search_wall_seconds;
        solver->profile.kd_result_order_hash =
            solver_order_hash_mix(
                solver->profile.kd_result_order_hash,
                solver_kd_result_order_digest(result));
    }
    solver->profile.codekd_hits +=
        (unsigned long long)result->nres;
    if (result->nres) {
        solver->profile.resolve_calls++;
    }
}

void solver_codekd_packet_resolve_result_range_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    kdtree_qres_t* result,
    int dimquad,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    double pixvals[DQMAX * 2];
    double resolve_wall_start = 0.0;
    int j;

    if (!solver || !descriptor || !result ||
        dimquad <= 0 || dimquad > DQMAX ||
        candidate_first < 0 ||
        candidate_end <= candidate_first ||
        candidate_end > result->nres) {
        return;
    }
    for (j = 0; j < dimquad; j++) {
        setx(
            pixvals,
            j,
            field_getx(solver, descriptor->stars[j]));
        sety(
            pixvals,
            j,
            field_gety(solver, descriptor->stars[j]));
    }
    if (solver->profile.detailed) {
        resolve_wall_start = monotonic_seconds();
    }
    resolve_matches_native_range(
        result,
        pixvals,
        descriptor->stars,
        dimquad,
        solver->numtries,
        solver,
        descriptor->current_parity,
        candidate_first,
        candidate_end,
        prepared,
        prepared_first,
        prepared_quad_count,
        prepared_star_count);
    if (solver->profile.detailed) {
        solver->profile.resolve_wall_seconds +=
            monotonic_seconds() - resolve_wall_start;
    }
}
