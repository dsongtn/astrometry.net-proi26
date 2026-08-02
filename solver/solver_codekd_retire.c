/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* Ordered CodeKD result retirement and owner reduction. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "solver.h"
#include "verify.h"
#include "fitsbin.h"
#include "log.h"
#include "index_shard_internal.h"
#include "solver_codekd_internal.h"
#include "solver_inline_internal.h"
#include "solver_hypothesis_internal.h"

static index_shard_staged_retire_status_t
solver_codekd_packet_retire_nonwindow_descriptor(
    solver_codekd_search_packet_t* packet,
    solver_codekd_packet_retire_context_t* context,
    const solver_ab_descriptor_t* descriptor,
    const solver_codekd_result_slot_t* slot) {
    if (!packet || !context || !context->solver ||
        !context->query_result || !context->reduced ||
        !descriptor || !slot ||
        packet->retire_descriptor_started ||
        (slot->state == SOLVER_CODEKD_RESULT_READY &&
         slot->hit_count)) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    if (slot->state != SOLVER_CODEKD_RESULT_QUERY_FAILED &&
        solver_poll_worker_stop(context->solver)) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_RETIRE_STOPPED;
    }
    context->solver->rel_field_noise2 =
        descriptor->rel_field_noise2;
    if (solver_ab_checked_counter_delta(
            context->solver,
            descriptor->numtries_delta,
            descriptor->cxdx_delta,
            descriptor->meanx_delta)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    if (slot->state == SOLVER_CODEKD_RESULT_READY) {
        kdtree_qres_t result_view;

        memset(&result_view, 0, sizeof(result_view));
        solver_codekd_packet_begin_result_owner(
            context->solver,
            descriptor,
            slot,
            &result_view,
            context->dimquads);
        context->solver->profile.hypotheses_reduced++;
    } else if (slot->state ==
               SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
        solver_execute_hypothesis_owner(
            descriptor->stars,
            descriptor->code,
            context->dimquads,
            context->solver,
            descriptor->current_parity,
            descriptor->tol2,
            context->query_result);
    } else if (slot->state ==
               SOLVER_CODEKD_RESULT_QUERY_FAILED) {
        solver_retire_codekd_failure_owner(
            descriptor->stars,
            descriptor->code,
            context->dimquads,
            context->solver,
            descriptor->current_parity,
            slot->search_errno,
            slot->search_wall_seconds);
    } else {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    (*context->reduced)++;
    packet->retire_descriptor++;
    if (context->solver->profile.execution_failed) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    if (context->solver->quit_now) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_RETIRE_STOPPED;
    }
    return INDEX_SHARD_STAGED_RETIRE_OK;
}

static index_shard_staged_retire_status_t
solver_codekd_search_packet_retire_window(
    solver_codekd_search_packet_t* packet,
    solver_codekd_packet_retire_context_t* context) {
    size_t window_end;

    if (!packet || !context || !context->solver ||
        packet->state != SOLVER_CODEKD_PACKET_RETIRING ||
        packet->candidate_count != packet->hit_count ||
        !packet->candidate_count ||
        !packet->candidate_records ||
        !packet->candidate_capacity ||
        packet->candidate_cursor >= packet->candidate_count ||
        !packet->candidate_window_count ||
        packet->candidate_window_first >
            packet->candidate_cursor ||
        packet->candidate_window_offset !=
            packet->candidate_cursor -
                packet->candidate_window_first ||
        packet->candidate_window_offset >=
            packet->candidate_window_count ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_quad_ready_count >
            packet->candidate_window_count ||
        packet->candidate_star_ready_count >
            packet->candidate_quad_ready_count ||
        (packet->verify_plan_complete &&
         (packet->verification_delivery_disabled ||
          packet->verify_plan_first !=
              packet->candidate_window_offset ||
          packet->verify_plan_end <=
              packet->verify_plan_first ||
          packet->verify_plan_end >
              packet->verify_topology_end ||
          packet->verify_topology_end >
              packet->candidate_window_count)) ||
        (!packet->verification_delivery_disabled &&
         packet->candidate_verify_query_count &&
         !packet->verify_plan_complete) ||
        packet->retire_descriptor >= packet->count) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    window_end =
        packet->candidate_window_first +
            packet->candidate_window_count;

    while (packet->retire_descriptor < packet->count) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[
                packet->retire_descriptor];
        const solver_codekd_result_slot_t* slot =
            &packet->slots[packet->retire_descriptor];
        index_shard_staged_retire_status_t status;

        if (slot->state == SOLVER_CODEKD_RESULT_READY &&
            slot->hit_count) {
            kdtree_qres_t result_view;
            size_t global_first;
            size_t local_end;
            size_t prepared_quad_count;
            size_t prepared_star_count;
            size_t range_count;
            size_t processed_count;
            int nummatches_before;

            if (!packet->candidate_window_count &&
                !packet->candidate_quad_delivery_disabled) {
                if (solver_codekd_packet_begin_candidate_window(
                        packet)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first ||
                packet->retire_hit_offset >
                    (size_t)slot->hit_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            memset(&result_view, 0, sizeof(result_view));
            result_view.nres = slot->hit_count;
            result_view.capacity = slot->hit_count;
            result_view.inds = packet->inds + slot->hit_first;
            result_view.sdists = packet->sdists + slot->hit_first;

            if (!packet->retire_descriptor_started) {
                if (solver_poll_worker_stop(context->solver)) {
                    packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                    return INDEX_SHARD_STAGED_RETIRE_STOPPED;
                }
                context->solver->rel_field_noise2 =
                    descriptor->rel_field_noise2;
                if (solver_ab_checked_counter_delta(
                        context->solver,
                        descriptor->numtries_delta,
                        descriptor->cxdx_delta,
                        descriptor->meanx_delta)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                solver_codekd_packet_begin_result_owner(
                    context->solver,
                    descriptor,
                    slot,
                    &result_view,
                    context->dimquads);
                packet->retire_descriptor_started = TRUE;
            }

            global_first =
                slot->hit_first + packet->retire_hit_offset;
            if (global_first != packet->candidate_cursor) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->candidate_quad_delivery_disabled) {
                local_end = (size_t)slot->hit_count;
                range_count = local_end -
                    packet->retire_hit_offset;
                nummatches_before =
                    context->solver->nummatches;
                solver_codekd_packet_resolve_result_range_owner(
                    context->solver,
                    descriptor,
                    &result_view,
                    context->dimquads,
                    (int)packet->retire_hit_offset,
                    (int)local_end,
                    NULL,
                    0U,
                    0U,
                    0U);
                if (context->solver->nummatches <
                    nummatches_before) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                processed_count = (size_t)(
                    context->solver->nummatches -
                    nummatches_before);
                if (processed_count > range_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                packet->candidate_retired_rows +=
                    processed_count;
                packet->candidate_native_rows +=
                    processed_count;
                packet->retire_hit_offset +=
                    processed_count;
                packet->candidate_cursor +=
                    processed_count;
                solver_codekd_packet_clear_candidate_window(
                    packet);

                if (context->solver->profile.execution_failed) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                if (context->solver->quit_now ||
                    solver_poll_worker_stop(context->solver)) {
                    context->solver->profile.hypotheses_reduced++;
                    (*context->reduced)++;
                    packet->retire_descriptor_started = FALSE;
                    packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                    return INDEX_SHARD_STAGED_RETIRE_STOPPED;
                }
                if (processed_count != range_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                context->solver->profile.hypotheses_reduced++;
                (*context->reduced)++;
                packet->retire_descriptor_started = FALSE;
                packet->retire_hit_offset = 0U;
                packet->retire_descriptor++;
                continue;
            }
            if (packet->candidate_cursor >= window_end ||
                packet->candidate_window_offset >=
                    packet->candidate_window_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            range_count = MIN(
                (size_t)slot->hit_count -
                    packet->retire_hit_offset,
                window_end - packet->candidate_cursor);
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled) {
                size_t verify_remaining =
                    packet->verify_plan_end -
                        packet->candidate_window_offset;

                range_count = MIN(
                    range_count, verify_remaining);
            }
            if (!range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            local_end =
                packet->retire_hit_offset +
                    range_count;
            prepared_quad_count =
                packet->candidate_quad_ready_count >
                    packet->candidate_window_offset
                ? MIN(
                    range_count,
                    packet->candidate_quad_ready_count -
                        packet->candidate_window_offset)
                : 0U;
            prepared_star_count =
                packet->candidate_star_ready_count >
                    packet->candidate_window_offset
                ? MIN(
                    range_count,
                    packet->candidate_star_ready_count -
                        packet->candidate_window_offset)
                : 0U;
            nummatches_before = context->solver->nummatches;
            solver_codekd_packet_resolve_result_range_owner(
                context->solver,
                descriptor,
                &result_view,
                context->dimquads,
                (int)packet->retire_hit_offset,
                (int)local_end,
                packet->candidate_records +
                    packet->candidate_window_offset,
                packet->retire_hit_offset,
                prepared_quad_count,
                prepared_star_count);
            if (context->solver->nummatches < nummatches_before) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            processed_count = (size_t)(
                context->solver->nummatches - nummatches_before);
            if (processed_count > range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            packet->candidate_retired_rows += processed_count;
            packet->candidate_native_rows +=
                processed_count > prepared_star_count
                ? processed_count -
                    prepared_star_count
                : 0U;
            packet->retire_hit_offset += processed_count;
            packet->candidate_cursor += processed_count;
            packet->candidate_window_offset += processed_count;

            if (context->solver->profile.execution_failed) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (context->solver->quit_now ||
                solver_poll_worker_stop(context->solver)) {
                context->solver->profile.hypotheses_reduced++;
                (*context->reduced)++;
                packet->retire_descriptor_started = FALSE;
                packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                return INDEX_SHARD_STAGED_RETIRE_STOPPED;
            }
            if (processed_count != range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled &&
                packet->candidate_window_offset ==
                    packet->verify_plan_end &&
                solver_codekd_packet_retired_verification_complete(
                    packet)) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled &&
                packet->candidate_window_offset ==
                    packet->verify_plan_end &&
                packet->candidate_window_offset <
                    packet->candidate_window_count) {
                size_t topology_end =
                    packet->verify_topology_end;

                if (packet->retire_hit_offset ==
                    (size_t)slot->hit_count) {
                    context->solver->profile.hypotheses_reduced++;
                    (*context->reduced)++;
                    packet->retire_descriptor_started = FALSE;
                    packet->retire_hit_offset = 0U;
                    packet->retire_descriptor++;
                }
                solver_codekd_packet_reset_verify_plan(packet);
                if (packet->candidate_window_offset <
                    topology_end) {
                    packet->verify_plan_first =
                        packet->candidate_window_offset;
                    packet->verify_plan_end = topology_end;
                    packet->verify_topology_end = topology_end;
                    packet->verify_plan_complete = TRUE;
                    packet->state =
                        SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
                } else {
                    packet->state =
                        SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            if (packet->candidate_window_offset ==
                packet->candidate_window_count) {
                solver_codekd_packet_clear_candidate_window(
                    packet);
            }
            if (packet->retire_hit_offset <
                (size_t)slot->hit_count) {
                if (packet->candidate_window_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                if (solver_codekd_packet_begin_candidate_window(
                        packet)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            context->solver->profile.hypotheses_reduced++;
            (*context->reduced)++;
            packet->retire_descriptor_started = FALSE;
            packet->retire_hit_offset = 0U;
            packet->retire_descriptor++;
            continue;
        }

        status = solver_codekd_packet_retire_nonwindow_descriptor(
            packet, context, descriptor, slot);
        if (status != INDEX_SHARD_STAGED_RETIRE_OK) {
            return status;
        }
    }

    if (packet->candidate_cursor != packet->candidate_count) {
        if (packet->candidate_window_count) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        if (solver_codekd_packet_begin_candidate_window(packet)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        return INDEX_SHARD_STAGED_RETIRE_MORE;
    }
    if (packet->descriptors->has_final_rel_field_noise2) {
        context->solver->rel_field_noise2 =
            packet->descriptors->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            context->solver,
            packet->descriptors->trailing_numtries,
            packet->descriptors->trailing_cxdx,
            packet->descriptors->trailing_meanx)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    return INDEX_SHARD_STAGED_RETIRE_OK;
}

index_shard_staged_retire_status_t
solver_codekd_search_packet_retire(
    const index_shard_staged_task_t* task,
    size_t task_index,
    void* owner_context) {
    solver_codekd_packet_retire_context_t* context = owner_context;
    const solver_codekd_packet_task_input_t* input;
    solver_codekd_search_packet_t* packet;
    size_t end;
    size_t descriptor_index;

    if (!task || !context || !context->solver ||
        !context->query_result || !context->reduced ||
        !task->input || task->input_bytes != sizeof(*input) ||
        !task->output || task->output_bytes != sizeof(*packet) ||
        task_index != context->next_task_index) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    input = task->input;
    packet = task->output;
    if (input->descriptor.combination_first != context->next_sequence ||
        input->descriptor.combination_first >=
            input->descriptor.combination_end ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet->descriptors == NULL ||
        packet->count != packet->descriptors->descriptor_count ||
        packet->first != 0U ||
        packet->next_descriptor != packet->count ||
        packet->plan_complete || packet->delivery_ticket ||
        packet->pending_descriptor_plan ||
        packet->pending_descriptor_raw_ranges ||
        packet->pending_descriptor_logical_bytes ||
        packet->delivery_source ||
        (packet->candidate_count &&
         packet->candidate_count != packet->hit_count) ||
        (packet->candidate_count && !packet->candidate_records) ||
        (packet->candidate_count && !packet->candidate_capacity) ||
        packet->candidate_count > packet->hit_count ||
        packet->candidate_cursor > packet->candidate_count ||
        packet->candidate_window_first >
            packet->candidate_cursor ||
        packet->candidate_window_offset !=
            packet->candidate_cursor -
                packet->candidate_window_first ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_offset >
            packet->candidate_window_count ||
        packet->candidate_quad_ready_count >
            packet->candidate_window_count ||
        packet->candidate_star_ready_count >
            packet->candidate_quad_ready_count ||
        packet->retire_descriptor > packet->count) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    packet->state = SOLVER_CODEKD_PACKET_RETIRING;
    context->solver->profile.max_batch_hypotheses = MAX(
        context->solver->profile.max_batch_hypotheses,
        packet->count);
    if (packet->candidate_count) {
        index_shard_staged_retire_status_t status =
            solver_codekd_search_packet_retire_window(
                packet, context);

        if (status == INDEX_SHARD_STAGED_RETIRE_OK) {
            packet->state = SOLVER_CODEKD_PACKET_RETIRED;
            context->next_task_index++;
            context->next_sequence =
                input->descriptor.combination_end;
        }
        return status;
    }
    end = packet->first + packet->count;
    for (descriptor_index = packet->first;
         descriptor_index < end;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[descriptor_index];
        const solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];

        if (slot->state != SOLVER_CODEKD_RESULT_QUERY_FAILED &&
            solver_poll_worker_stop(context->solver)) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_STAGED_RETIRE_STOPPED;
        }
        context->solver->rel_field_noise2 =
            descriptor->rel_field_noise2;
        if (solver_ab_checked_counter_delta(
                context->solver,
                descriptor->numtries_delta,
                descriptor->cxdx_delta,
                descriptor->meanx_delta)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }

        if (slot->state == SOLVER_CODEKD_RESULT_READY) {
            kdtree_qres_t result_view;
            solver_candidate_delivery_record_t* prepared = NULL;
            size_t prepared_quad_count = 0U;
            size_t prepared_star_count = 0U;

            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            memset(&result_view, 0, sizeof(result_view));
            result_view.nres = slot->hit_count;
            result_view.capacity = slot->hit_count;
            result_view.inds = packet->inds + slot->hit_first;
            result_view.sdists = packet->sdists + slot->hit_first;
            if (slot->hit_first <
                packet->candidate_quad_ready_count) {
                prepared =
                    packet->candidate_records + slot->hit_first;
                prepared_quad_count = MIN(
                    (size_t)slot->hit_count,
                    packet->candidate_quad_ready_count -
                        slot->hit_first);
                if (slot->hit_first <
                    packet->candidate_star_ready_count) {
                    prepared_star_count = MIN(
                        (size_t)slot->hit_count,
                        packet->candidate_star_ready_count -
                            slot->hit_first);
                }
            }
            solver_execute_prepared_hypothesis_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                &result_view,
                slot->search_wall_seconds,
                prepared,
                prepared_quad_count,
                prepared_star_count);
        } else if (slot->state ==
                   SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            solver_execute_hypothesis_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                descriptor->tol2,
                context->query_result);
        } else if (slot->state ==
                   SOLVER_CODEKD_RESULT_QUERY_FAILED) {
            solver_retire_codekd_failure_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                slot->search_errno,
                slot->search_wall_seconds);
        } else {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }

        (*context->reduced)++;
        if (context->solver->profile.execution_failed) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        if (context->solver->quit_now) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_STAGED_RETIRE_STOPPED;
        }
    }

    if (packet->descriptors->has_final_rel_field_noise2) {
        context->solver->rel_field_noise2 =
            packet->descriptors->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            context->solver,
            packet->descriptors->trailing_numtries,
            packet->descriptors->trailing_cxdx,
            packet->descriptors->trailing_meanx)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    packet->state = SOLVER_CODEKD_PACKET_RETIRED;
    context->next_task_index++;
    context->next_sequence = input->descriptor.combination_end;
    return INDEX_SHARD_STAGED_RETIRE_OK;
}
