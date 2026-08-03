/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>

#include "bl.h"
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

typedef struct solver_codekd_test_gate {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    int entered;
    int released;
} solver_codekd_test_gate_t;

typedef struct solver_codekd_test_completion {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    unsigned long long maximum_id;
} solver_codekd_test_completion_t;

typedef struct solver_codekd_test_fixture {
    fitsbin_t* source;
    kdtree_t* tree;
    double* tree_data;
    void* mapping;
    size_t mapping_size;
    solver_ab_descriptor_output_t* descriptors;
    solver_codekd_packet_task_input_t input;
    solver_codekd_search_packet_t packet;
} solver_codekd_test_fixture_t;

static int solver_codekd_test_gate_init(
    solver_codekd_test_gate_t* gate) {
    if (!gate) {
        return -1;
    }
    memset(gate, 0, sizeof(*gate));
    if (pthread_mutex_init(&gate->mutex, NULL)) {
        return -1;
    }
    if (pthread_cond_init(&gate->condition, NULL)) {
        pthread_mutex_destroy(&gate->mutex);
        return -1;
    }
    return 0;
}

static void solver_codekd_test_gate_release(
    solver_codekd_test_gate_t* gate) {
    if (!gate) {
        return;
    }
    pthread_mutex_lock(&gate->mutex);
    gate->released = TRUE;
    pthread_cond_broadcast(&gate->condition);
    pthread_mutex_unlock(&gate->mutex);
}

static void solver_codekd_test_deadline(struct timespec* deadline) {
    if (!deadline || clock_gettime(CLOCK_REALTIME, deadline)) {
        return;
    }
    deadline->tv_sec += 5;
}

static int solver_codekd_test_gate_wait_entered(
    solver_codekd_test_gate_t* gate) {
    struct timespec deadline;
    int status = 0;

    if (!gate) {
        return -1;
    }
    memset(&deadline, 0, sizeof(deadline));
    solver_codekd_test_deadline(&deadline);
    pthread_mutex_lock(&gate->mutex);
    while (!gate->entered && !status) {
        status = pthread_cond_timedwait(
            &gate->condition, &gate->mutex, &deadline);
    }
    status = gate->entered ? 0 : -1;
    pthread_mutex_unlock(&gate->mutex);
    return status;
}

static void solver_codekd_test_gate_destroy(
    solver_codekd_test_gate_t* gate) {
    if (!gate) {
        return;
    }
    solver_codekd_test_gate_release(gate);
    pthread_cond_destroy(&gate->condition);
    pthread_mutex_destroy(&gate->mutex);
}

static int solver_codekd_test_completion_init(
    solver_codekd_test_completion_t* completion) {
    if (!completion) {
        return -1;
    }
    memset(completion, 0, sizeof(*completion));
    if (pthread_mutex_init(&completion->mutex, NULL)) {
        return -1;
    }
    if (pthread_cond_init(&completion->condition, NULL)) {
        pthread_mutex_destroy(&completion->mutex);
        return -1;
    }
    return 0;
}

static void solver_codekd_test_completion_notify(
    void* opaque,
    unsigned long long completion_id) {
    solver_codekd_test_completion_t* completion = opaque;

    if (!completion || !completion_id) {
        return;
    }
    pthread_mutex_lock(&completion->mutex);
    if (completion_id > completion->maximum_id) {
        completion->maximum_id = completion_id;
    }
    pthread_cond_broadcast(&completion->condition);
    pthread_mutex_unlock(&completion->mutex);
}

static int solver_codekd_test_completion_wait(
    solver_codekd_test_completion_t* completion,
    unsigned long long completion_id) {
    struct timespec deadline;
    int status = 0;

    if (!completion || !completion_id) {
        return -1;
    }
    memset(&deadline, 0, sizeof(deadline));
    solver_codekd_test_deadline(&deadline);
    pthread_mutex_lock(&completion->mutex);
    while (completion->maximum_id < completion_id && !status) {
        status = pthread_cond_timedwait(
            &completion->condition,
            &completion->mutex,
            &deadline);
    }
    status = completion->maximum_id >= completion_id ? 0 : -1;
    pthread_mutex_unlock(&completion->mutex);
    return status;
}

static void solver_codekd_test_completion_destroy(
    solver_codekd_test_completion_t* completion) {
    if (!completion) {
        return;
    }
    pthread_cond_destroy(&completion->condition);
    pthread_mutex_destroy(&completion->mutex);
}

static int solver_codekd_test_empty_plan(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count) {
    (void)opaque;
    (void)ranges;
    (void)range_capacity;
    if (!cancelled || !range_count) {
        errno = EINVAL;
        return -1;
    }
    *range_count = 0U;
    if (cancelled(cancel_opaque)) {
        errno = ECANCELED;
        return -1;
    }
    return 0;
}

static int solver_codekd_test_blocking_plan(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count) {
    solver_codekd_test_gate_t* gate = opaque;

    (void)ranges;
    (void)range_capacity;
    if (!gate || !cancelled || !range_count) {
        errno = EINVAL;
        return -1;
    }
    *range_count = 0U;
    pthread_mutex_lock(&gate->mutex);
    gate->entered = TRUE;
    pthread_cond_broadcast(&gate->condition);
    while (!gate->released && !cancelled(cancel_opaque)) {
        pthread_cond_wait(&gate->condition, &gate->mutex);
    }
    pthread_mutex_unlock(&gate->mutex);
    if (cancelled(cancel_opaque)) {
        errno = ECANCELED;
        return -1;
    }
    return 0;
}

static anbool solver_codekd_test_never_cancel(void* opaque) {
    (void)opaque;
    return FALSE;
}

static anbool solver_codekd_test_always_cancel(void* opaque) {
    (void)opaque;
    return TRUE;
}

static int solver_codekd_test_fixture_init(
    solver_codekd_test_fixture_t* fixture,
    size_t descriptor_count) {
    fitsbin_chunk_t chunk;
    size_t data_index;
    long page_size;

    if (!fixture || !descriptor_count ||
        descriptor_count > SOLVER_AB_DESCRIPTOR_CAPACITY) {
        return -1;
    }
    memset(fixture, 0, sizeof(*fixture));
    page_size = sysconf(_SC_PAGESIZE);
    if (page_size <= 0) {
        return -1;
    }
    fixture->mapping_size = (size_t)page_size;
    fixture->mapping = mmap(
        NULL,
        fixture->mapping_size,
        PROT_READ | PROT_WRITE,
        MAP_PRIVATE | MAP_ANONYMOUS,
        -1,
        0);
    if (fixture->mapping == MAP_FAILED) {
        fixture->mapping = NULL;
        return -1;
    }
    fixture->source = calloc(1, sizeof(*fixture->source));
    if (!fixture->source) {
        goto fail;
    }
    fixture->source->payload_fd = -1;
    fixture->source->payload_fd_initialized = TRUE;
    fixture->source->filename = strdup("codekd-detached-test");
    fixture->source->chunks = bl_new(1, sizeof(fitsbin_chunk_t));
    if (!fixture->source->filename || !fixture->source->chunks) {
        goto fail;
    }
    fitsbin_chunk_init(&chunk);
    chunk.tablename = "codekd-detached";
    chunk.map = fixture->mapping;
    chunk.mapsize = fixture->mapping_size;
    chunk.data = fixture->mapping;
    chunk.itemsize = 1;
    chunk.nrows = (int)fixture->mapping_size;
    chunk.data_file_offset = 0;
    chunk.data_file_size = fixture->mapping_size;
    if (!fitsbin_add_chunk(fixture->source, &chunk)) {
        goto fail;
    }
    fixture->mapping = NULL;
    if (fitsbin_configure_index_mmap(fixture->source)) {
        goto fail;
    }

    fixture->tree_data = calloc(4U * DCMAX, sizeof(*fixture->tree_data));
    if (!fixture->tree_data) {
        goto fail;
    }
    for (data_index = DCMAX; data_index < 4U * DCMAX; data_index++) {
        fixture->tree_data[data_index] =
            (double)(data_index / DCMAX);
    }
    fixture->tree = kdtree_build(
        NULL,
        fixture->tree_data,
        4,
        DCMAX,
        4,
        KDTT_DOUBLE,
        KD_BUILD_SPLIT);
    if (!fixture->tree) {
        goto fail;
    }
    fixture->tree->io = fixture->source;
    fixture->tree->io_is_fitsbin = TRUE;

    fixture->descriptors = calloc(1, sizeof(*fixture->descriptors));
    fixture->packet.slots = calloc(
        descriptor_count, sizeof(*fixture->packet.slots));
    fixture->packet.inds = calloc(8U, sizeof(*fixture->packet.inds));
    fixture->packet.sdists = calloc(8U, sizeof(*fixture->packet.sdists));
    if (!fixture->descriptors || !fixture->packet.slots ||
        !fixture->packet.inds || !fixture->packet.sdists ||
        solver_codekd_page_workspace_create(
            &fixture->packet.page_workspace)) {
        goto fail;
    }
    fixture->descriptors->descriptor_count = descriptor_count;
    for (data_index = 0U; data_index < descriptor_count; data_index++) {
        fixture->descriptors->descriptors[data_index].tol2 = 0.25;
    }
    fixture->packet.descriptors = fixture->descriptors;
    fixture->packet.tree = fixture->tree;
    fixture->packet.hit_capacity = 8U;
    fixture->packet.first = 0U;
    fixture->packet.count = descriptor_count;
    fixture->packet.sequence = 17ULL;
    fixture->packet.state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    fixture->input.descriptor.combination_first = fixture->packet.sequence;
    fixture->input.tree = fixture->tree;
    return 0;

fail:
    (void)solver_codekd_search_packet_cleanup(&fixture->packet);
    free(fixture->descriptors);
    fixture->descriptors = NULL;
    kdtree_free(fixture->tree);
    fixture->tree = NULL;
    free(fixture->tree_data);
    fixture->tree_data = NULL;
    if (fixture->source && fixture->source->chunks) {
        fitsbin_close(fixture->source);
        fixture->source = NULL;
    } else {
        if (fixture->source) {
            free(fixture->source->filename);
            free(fixture->source);
            fixture->source = NULL;
        }
    }
    if (fixture->mapping) {
        munmap(fixture->mapping, fixture->mapping_size);
    }
    fixture->mapping = NULL;
    return -1;
}

static void solver_codekd_test_fixture_cleanup(
    solver_codekd_test_fixture_t* fixture) {
    if (!fixture) {
        return;
    }
    (void)solver_codekd_search_packet_cleanup(&fixture->packet);
    free(fixture->descriptors);
    kdtree_free(fixture->tree);
    free(fixture->tree_data);
    if (fixture->source) {
        fitsbin_close(fixture->source);
    }
    memset(fixture, 0, sizeof(*fixture));
}

static int solver_codekd_test_result_matches_native(
    const solver_codekd_test_fixture_t* fixture,
    size_t descriptor_index,
    const kdtree_qres_t* native_result) {
    const solver_codekd_result_slot_t* slot;
    size_t hit_first;
    size_t hit_count;

    if (!fixture || !native_result || native_result->nres < 0 ||
        descriptor_index >= fixture->packet.count) {
        return FALSE;
    }
    slot = &fixture->packet.slots[descriptor_index];
    hit_first = slot->hit_first;
    hit_count = (size_t)slot->hit_count;
    if (slot->state != SOLVER_CODEKD_RESULT_READY ||
        hit_count != (size_t)native_result->nres ||
        hit_first > fixture->packet.hit_count ||
        hit_count > fixture->packet.hit_count - hit_first) {
        return FALSE;
    }
    if (!hit_count) {
        return TRUE;
    }
    return native_result->inds && native_result->sdists &&
        !memcmp(
            fixture->packet.inds + hit_first,
            native_result->inds,
            hit_count * sizeof(*fixture->packet.inds)) &&
        !memcmp(
            fixture->packet.sdists + hit_first,
            native_result->sdists,
            hit_count * sizeof(*fixture->packet.sdists));
}

static int solver_codekd_test_seed_plan(
    solver_codekd_test_fixture_t* fixture) {
    solver_codekd_page_workspace_t* workspace;
    solver_codekd_page_entry_t* entry;
    fitsbin_chunk_t* chunk;
    uintptr_t begin;

    if (!fixture || !fixture->source ||
        !fixture->packet.page_workspace) {
        return -1;
    }
    chunk = fitsbin_get_chunk(fixture->source, 0);
    if (!chunk || !chunk->map || !chunk->mapsize) {
        return -1;
    }
    workspace = fixture->packet.page_workspace;
    solver_codekd_page_set_reset(&workspace->descriptor);
    if (!workspace->descriptor.entries ||
        !workspace->descriptor.capacity) {
        return -1;
    }
    begin = (uintptr_t)chunk->map;
    entry = &workspace->descriptor.entries[0];
    memset(entry, 0, sizeof(*entry));
    entry->mapping_begin = begin;
    entry->mapping_end = begin + chunk->mapsize;
    entry->page_key = begin;
    entry->populate_begin = begin;
    entry->populate_end = begin + chunk->mapsize;
    workspace->descriptor.count = 1U;
    workspace->page_limit = 1U;
    fixture->packet.pending_descriptor_plan = TRUE;
    fixture->packet.pending_descriptor_raw_ranges = 1U;
    fixture->packet.pending_descriptor_logical_bytes = 1U;
    return 0;
}

static void solver_codekd_test_ticket_drain(
    fitsbin_t* source,
    fitsbin_payload_io_ticket_t** ticket) {
    if (!source || !ticket || !*ticket) {
        return;
    }
    (void)fitsbin_payload_io_ticket_cancel_and_wait(source, *ticket);
    fitsbin_payload_io_ticket_destroy(*ticket);
    *ticket = NULL;
}

static int solver_codekd_test_detached_success(
    int* queued,
    int* compute_ready,
    int* repeated) {
    solver_codekd_test_fixture_t fixture;
    solver_codekd_test_gate_t gate;
    solver_codekd_test_completion_t completion;
    fitsbin_payload_io_ticket_t* blocker = NULL;
    kdtree_qres_t* native_result = NULL;
    unsigned long long completion_id = 0ULL;
    size_t cycle;
    int notifier_set = FALSE;
    int service_started = FALSE;
    int gate_initialized = FALSE;
    int completion_initialized = FALSE;
    int status;
    int result = -1;

    if (!queued || !compute_ready || !repeated) {
        return -1;
    }
    *queued = FALSE;
    *compute_ready = FALSE;
    *repeated = FALSE;
    if (solver_codekd_test_fixture_init(&fixture, 2U)) {
        return -1;
    }
    native_result = solver_codekd_rangesearch(
        fixture.tree,
        NULL,
        fixture.descriptors->descriptors[0].code,
        fixture.descriptors->descriptors[0].tol2,
        SOLVER_CODEKD_SEARCH_OPTIONS);
    if (!native_result) {
        goto cleanup;
    }
    if (solver_codekd_test_gate_init(&gate)) {
        goto cleanup;
    }
    gate_initialized = TRUE;
    if (solver_codekd_test_completion_init(&completion)) {
        goto cleanup;
    }
    completion_initialized = TRUE;
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    if (fitsbin_payload_io_set_completion_notifier(
            solver_codekd_test_completion_notify, &completion)) {
        goto cleanup;
    }
    notifier_set = TRUE;
    if (fitsbin_payload_io_service_start(1)) {
        goto cleanup;
    }
    service_started = TRUE;
    if (fitsbin_prefetch_ranges_planned_submit(
            fixture.source,
            solver_codekd_test_blocking_plan,
            &gate,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            &blocker) != FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED ||
        !blocker || solver_codekd_test_gate_wait_entered(&gate)) {
        goto cleanup;
    }

    for (cycle = 0U; cycle < 2U; cycle++) {
        if (solver_codekd_test_seed_plan(&fixture)) {
            goto cleanup;
        }
        status = solver_codekd_packet_staged_ops.prepare(
            &fixture.input,
            sizeof(fixture.input),
            &fixture.packet,
            sizeof(fixture.packet));
        if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
            goto cleanup;
        }
        status = solver_codekd_packet_staged_ops.submit(
            &fixture.input,
            sizeof(fixture.input),
            &fixture.packet,
            sizeof(fixture.packet),
            &completion_id);
        if (status != INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED ||
            !completion_id) {
            goto cleanup;
        }
        if (!cycle) {
            *queued =
                fixture.packet.state ==
                    SOLVER_CODEKD_PACKET_CODEKD_IO_SUBMITTED &&
                fixture.packet.next_descriptor == 0U &&
                !fixture.packet.plan_complete &&
                fixture.packet.delivery_ticket != NULL &&
                fixture.packet.delivery_source == fixture.source;
            solver_codekd_test_gate_release(&gate);
            solver_codekd_test_ticket_drain(fixture.source, &blocker);
        }
        if (solver_codekd_test_completion_wait(
                &completion, completion_id)) {
            goto cleanup;
        }
        status = solver_codekd_packet_staged_ops.poll(
            &fixture.input,
            sizeof(fixture.input),
            &fixture.packet,
            sizeof(fixture.packet));
        if (status != INDEX_SHARD_STAGED_IO_READY ||
            fixture.packet.state !=
                SOLVER_CODEKD_PACKET_COMPUTE_READY ||
            !fixture.packet.plan_complete) {
            goto cleanup;
        }
        if (!cycle) {
            *compute_ready = TRUE;
        }
        status = solver_codekd_packet_staged_ops.execute(
            &fixture.input,
            sizeof(fixture.input),
            &fixture.packet,
            sizeof(fixture.packet));
        if (status != INDEX_SHARD_STAGED_EXECUTE_MORE ||
            fixture.packet.state !=
                SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
            fixture.packet.next_descriptor != cycle + 1U ||
            fixture.packet.slots[cycle].state !=
                SOLVER_CODEKD_RESULT_READY ||
            !solver_codekd_test_result_matches_native(
                &fixture, cycle, native_result) ||
            fixture.packet.plan_complete ||
            fixture.packet.delivery_ticket ||
            fixture.packet.delivery_source) {
            goto cleanup;
        }
    }
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    *repeated =
        status == INDEX_SHARD_STAGED_PREPARE_RESULTS_READY &&
        fixture.packet.state == SOLVER_CODEKD_PACKET_RESULTS_READY &&
        fixture.packet.next_descriptor == fixture.packet.count &&
        fixture.packet.hit_count == 2U &&
        fixture.packet.page_stats.descriptors_planned == 2U;
    result = *queued && *compute_ready && *repeated ? 0 : -1;

cleanup:
    if (gate_initialized) {
        solver_codekd_test_gate_release(&gate);
    }
    solver_codekd_test_ticket_drain(fixture.source, &blocker);
    (void)solver_codekd_search_packet_release_ticket(&fixture.packet);
    if (service_started) {
        fitsbin_payload_io_service_stop();
    }
    if (notifier_set) {
        (void)fitsbin_payload_io_clear_completion_notifier(
            solver_codekd_test_completion_notify, &completion);
    }
    fitsbin_payload_io_configure_workers(1);
    if (completion_initialized) {
        solver_codekd_test_completion_destroy(&completion);
    }
    if (gate_initialized) {
        solver_codekd_test_gate_destroy(&gate);
    }
    kdtree_free_query(native_result);
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_fully_resident(void) {
    solver_codekd_test_fixture_t fixture;
    kdtree_qres_t* native_result = NULL;
    unsigned long long completion_id = 0ULL;
    int status;
    int result = -1;

    if (solver_codekd_test_fixture_init(&fixture, 1U)) {
        return -1;
    }
    native_result = solver_codekd_rangesearch(
        fixture.tree,
        NULL,
        fixture.descriptors->descriptors[0].code,
        fixture.descriptors->descriptors[0].tol2,
        SOLVER_CODEKD_SEARCH_OPTIONS);
    if (!native_result || solver_codekd_test_seed_plan(&fixture)) {
        goto cleanup;
    }
    fixture.source->payload_fully_resident = TRUE;
    fitsbin_payload_io_service_stop();
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY ||
        completion_id ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_COMPUTE_READY ||
        !fixture.packet.plan_complete ||
        fixture.packet.plan_first != 0U ||
        fixture.packet.plan_end != 1U ||
        !fixture.packet.plan_range_count ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.execute(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_EXECUTE_MORE ||
        fixture.packet.state !=
            SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
        fixture.packet.next_descriptor != 1U ||
        !solver_codekd_test_result_matches_native(
            &fixture, 0U, native_result) ||
        fixture.packet.plan_complete ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    result = 0;

cleanup:
    kdtree_free_query(native_result);
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_fully_resident_empty(void) {
    solver_codekd_test_fixture_t fixture;
    kdtree_t empty_tree;
    unsigned long long completion_id = 0ULL;
    size_t i;
    int status;
    int result = -1;

    if (solver_codekd_test_fixture_init(&fixture, 3U)) {
        return -1;
    }
    memset(&empty_tree, 0, sizeof(empty_tree));
    empty_tree.treetype = KDTT_DOUBLE;
    empty_tree.ndim = DCMAX;
    empty_tree.io = fixture.source;
    empty_tree.io_is_fitsbin = TRUE;
    fixture.packet.tree = &empty_tree;
    fixture.input.tree = &empty_tree;
    fixture.source->payload_fully_resident = TRUE;
    fitsbin_payload_io_service_stop();

    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_OWNER_READY ||
        completion_id ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        fixture.packet.next_descriptor != fixture.packet.count ||
        fixture.packet.plan_complete ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    for (i = 0U; i < fixture.packet.count; i++) {
        if (fixture.packet.slots[i].state !=
            SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            goto cleanup;
        }
    }
    status = solver_codekd_packet_staged_ops.owner(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_EXECUTE_OK ||
        fixture.packet.page_stats.refusal_counts[
            SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE] !=
                fixture.packet.count) {
        goto cleanup;
    }
    result = 0;

cleanup:
    fixture.packet.tree = fixture.tree;
    fixture.input.tree = fixture.tree;
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_detached_empty(void) {
    solver_codekd_test_fixture_t fixture;
    solver_codekd_test_completion_t completion;
    kdtree_t empty_tree;
    unsigned long long completion_id = 0ULL;
    size_t i;
    int notifier_set = FALSE;
    int service_started = FALSE;
    int completion_initialized = FALSE;
    int status;
    int result = -1;

    if (solver_codekd_test_fixture_init(&fixture, 3U)) {
        return -1;
    }
    memset(&empty_tree, 0, sizeof(empty_tree));
    empty_tree.treetype = KDTT_DOUBLE;
    empty_tree.ndim = DCMAX;
    empty_tree.io = fixture.source;
    empty_tree.io_is_fitsbin = TRUE;
    fixture.packet.tree = &empty_tree;
    fixture.input.tree = &empty_tree;
    if (solver_codekd_test_completion_init(&completion)) {
        goto cleanup;
    }
    completion_initialized = TRUE;
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    if (fitsbin_payload_io_set_completion_notifier(
            solver_codekd_test_completion_notify, &completion)) {
        goto cleanup;
    }
    notifier_set = TRUE;
    if (fitsbin_payload_io_service_start(1)) {
        goto cleanup;
    }
    service_started = TRUE;
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED ||
        !completion_id ||
        solver_codekd_test_completion_wait(
            &completion, completion_id)) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.poll(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_IO_FAILED ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        fixture.packet.next_descriptor != fixture.packet.count ||
        fixture.packet.plan_complete ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    for (i = 0U; i < fixture.packet.count; i++) {
        if (fixture.packet.slots[i].state !=
            SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            goto cleanup;
        }
    }
    status = solver_codekd_packet_staged_ops.owner(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_EXECUTE_OK ||
        fixture.packet.page_stats.refusal_counts[
            SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE] !=
                fixture.packet.count) {
        goto cleanup;
    }
    result = 0;

cleanup:
    (void)solver_codekd_search_packet_release_ticket(&fixture.packet);
    if (service_started) {
        fitsbin_payload_io_service_stop();
    }
    if (notifier_set) {
        (void)fitsbin_payload_io_clear_completion_notifier(
            solver_codekd_test_completion_notify, &completion);
    }
    fitsbin_payload_io_configure_workers(1);
    if (completion_initialized) {
        solver_codekd_test_completion_destroy(&completion);
    }
    fixture.packet.tree = fixture.tree;
    fixture.input.tree = fixture.tree;
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_detached_callback_error(void) {
    solver_codekd_test_fixture_t fixture;
    fitsbin_prefetch_range_t range;
    size_t range_count = 0U;
    int status;
    int saved_errno;
    int result;

    if (solver_codekd_test_fixture_init(&fixture, 1U) ||
        solver_codekd_test_seed_plan(&fixture)) {
        return -1;
    }
    memset(&range, 0, sizeof(range));
    fixture.packet.state = SOLVER_CODEKD_PACKET_CODEKD_IO_SUBMITTED;
    errno = 0;
    status = solver_codekd_packet_plan_codekd_pages(
        &fixture.packet,
        solver_codekd_test_never_cancel,
        NULL,
        &range,
        0U,
        &range_count);
    saved_errno = errno;
    result = status == -1 && saved_errno == E2BIG &&
        fixture.packet.state == SOLVER_CODEKD_PACKET_FAILED &&
        range_count == 0U
        ? 0
        : -1;
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_detached_refusal(void) {
    solver_codekd_test_fixture_t fixture;
    unsigned long long completion_id = 0ULL;
    size_t i;
    int status;
    int result = -1;

    if (solver_codekd_test_fixture_init(&fixture, 2U)) {
        return -1;
    }
    fitsbin_payload_io_service_stop();
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_OWNER_READY ||
        completion_id ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        fixture.packet.next_descriptor != fixture.packet.count ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    for (i = 0U; i < fixture.packet.count; i++) {
        if (fixture.packet.slots[i].state !=
            SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            goto cleanup;
        }
    }
    status = solver_codekd_packet_staged_ops.owner(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_EXECUTE_OK ||
        fixture.packet.page_stats.refusal_counts[
            SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED] != 1U) {
        goto cleanup;
    }
    result = 0;

cleanup:
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_detached_eagain(void) {
    enum { SOLVER_CODEKD_TEST_FILLER_CAPACITY = 64 };
    solver_codekd_test_fixture_t fixture;
    solver_codekd_test_gate_t gate;
    fitsbin_payload_io_ticket_t* blocker = NULL;
    fitsbin_payload_io_ticket_t* fillers[
        SOLVER_CODEKD_TEST_FILLER_CAPACITY];
    unsigned long long completion_id = 0ULL;
    size_t filler_count = 0U;
    size_t i;
    int service_started = FALSE;
    int gate_initialized = FALSE;
    int capacity_reached = FALSE;
    int status;
    int result = -1;

    memset(fillers, 0, sizeof(fillers));
    if (solver_codekd_test_fixture_init(&fixture, 1U)) {
        return -1;
    }
    if (solver_codekd_test_gate_init(&gate)) {
        goto cleanup;
    }
    gate_initialized = TRUE;
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    if (fitsbin_payload_io_service_start(1)) {
        goto cleanup;
    }
    service_started = TRUE;
    if (fitsbin_prefetch_ranges_planned_submit(
            fixture.source,
            solver_codekd_test_blocking_plan,
            &gate,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            &blocker) != FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED ||
        !blocker || solver_codekd_test_gate_wait_entered(&gate)) {
        goto cleanup;
    }
    while (filler_count < SOLVER_CODEKD_TEST_FILLER_CAPACITY) {
        errno = 0;
        status = fitsbin_prefetch_ranges_planned_submit(
            fixture.source,
            solver_codekd_test_empty_plan,
            NULL,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            &fillers[filler_count]);
        if (status == FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED &&
            fillers[filler_count]) {
            filler_count++;
            continue;
        }
        if (status == FITSBIN_PAYLOAD_IO_SUBMIT_UNAVAILABLE &&
            errno == EAGAIN && !fillers[filler_count]) {
            capacity_reached = TRUE;
        }
        break;
    }
    if (!capacity_reached) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_RETRY ||
        completion_id ||
        fixture.packet.state !=
            SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
        fixture.packet.next_descriptor != 0U ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source) {
        goto cleanup;
    }
    result = 0;

cleanup:
    if (gate_initialized) {
        solver_codekd_test_gate_release(&gate);
    }
    solver_codekd_test_ticket_drain(fixture.source, &blocker);
    for (i = 0U; i < filler_count; i++) {
        solver_codekd_test_ticket_drain(
            fixture.source, &fillers[i]);
    }
    if (service_started) {
        fitsbin_payload_io_service_stop();
    }
    fitsbin_payload_io_configure_workers(1);
    if (gate_initialized) {
        solver_codekd_test_gate_destroy(&gate);
    }
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

static int solver_codekd_test_detached_cancellation(void) {
    solver_codekd_test_fixture_t fixture;
    solver_codekd_test_gate_t gate;
    solver_codekd_test_completion_t completion;
    fitsbin_payload_io_ticket_t* blocker = NULL;
    fitsbin_prefetch_range_t range;
    unsigned long long completion_id = 0ULL;
    size_t range_count = 0U;
    int notifier_set = FALSE;
    int service_started = FALSE;
    int gate_initialized = FALSE;
    int completion_initialized = FALSE;
    int status;
    int result = -1;

    if (solver_codekd_test_fixture_init(&fixture, 1U)) {
        return -1;
    }
    memset(&range, 0, sizeof(range));
    fixture.packet.state = SOLVER_CODEKD_PACKET_CODEKD_IO_SUBMITTED;
    errno = 0;
    status = solver_codekd_packet_plan_codekd_pages(
        &fixture.packet,
        solver_codekd_test_always_cancel,
        NULL,
        &range,
        1U,
        &range_count);
    if (status != -1 || errno != ECANCELED || range_count ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_STOPPED ||
        fixture.packet.page_stats.refusal_counts[
            SOLVER_CODEKD_PAGE_PLAN_CANCELLED] != 1U) {
        goto cleanup;
    }
    fixture.packet.state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    fixture.packet.next_descriptor = 0U;
    fixture.packet.page_plan_reason = SOLVER_CODEKD_PAGE_PLAN_NONE;
    memset(fixture.packet.page_stats.refusal_counts, 0,
           sizeof(fixture.packet.page_stats.refusal_counts));
    if (solver_codekd_test_gate_init(&gate)) {
        goto cleanup;
    }
    gate_initialized = TRUE;
    if (solver_codekd_test_completion_init(&completion)) {
        goto cleanup;
    }
    completion_initialized = TRUE;
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    if (fitsbin_payload_io_set_completion_notifier(
            solver_codekd_test_completion_notify, &completion)) {
        goto cleanup;
    }
    notifier_set = TRUE;
    if (fitsbin_payload_io_service_start(1)) {
        goto cleanup;
    }
    service_started = TRUE;
    if (fitsbin_prefetch_ranges_planned_submit(
            fixture.source,
            solver_codekd_test_blocking_plan,
            &gate,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            &blocker) != FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED ||
        !blocker || solver_codekd_test_gate_wait_entered(&gate)) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.prepare(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.submit(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet),
        &completion_id);
    if (status != INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED ||
        !completion_id ||
        solver_codekd_packet_staged_ops.cancel(
            &fixture.input,
            sizeof(fixture.input),
            &fixture.packet,
            sizeof(fixture.packet))) {
        goto cleanup;
    }
    solver_codekd_test_gate_release(&gate);
    solver_codekd_test_ticket_drain(fixture.source, &blocker);
    if (solver_codekd_test_completion_wait(
            &completion, completion_id)) {
        goto cleanup;
    }
    status = solver_codekd_packet_staged_ops.poll(
        &fixture.input,
        sizeof(fixture.input),
        &fixture.packet,
        sizeof(fixture.packet));
    if (status != INDEX_SHARD_STAGED_IO_CANCELLED ||
        fixture.packet.state != SOLVER_CODEKD_PACKET_STOPPED ||
        fixture.packet.delivery_ticket ||
        fixture.packet.delivery_source ||
        fixture.packet.page_stats.refusal_counts[
            SOLVER_CODEKD_PAGE_PLAN_CANCELLED] != 1U) {
        goto cleanup;
    }
    result = 0;

cleanup:
    if (gate_initialized) {
        solver_codekd_test_gate_release(&gate);
    }
    solver_codekd_test_ticket_drain(fixture.source, &blocker);
    (void)solver_codekd_search_packet_release_ticket(&fixture.packet);
    if (service_started) {
        fitsbin_payload_io_service_stop();
    }
    if (notifier_set) {
        (void)fitsbin_payload_io_clear_completion_notifier(
            solver_codekd_test_completion_notify, &completion);
    }
    fitsbin_payload_io_configure_workers(1);
    if (completion_initialized) {
        solver_codekd_test_completion_destroy(&completion);
    }
    if (gate_initialized) {
        solver_codekd_test_gate_destroy(&gate);
    }
    solver_codekd_test_fixture_cleanup(&fixture);
    return result;
}

int solver_codekd_test_run_detached_initial_planning(
    solver_codekd_test_detached_result_t* result) {
    if (!result) {
        return -1;
    }
    memset(result, 0, sizeof(*result));
    result->supported =
        fitsbin_payload_io_mapped_population_supported();
    if (!result->supported) {
        return 0;
    }
    (void)solver_codekd_test_detached_success(
        &result->queued_callback,
        &result->success_compute_ready,
        &result->repeated_cycles);
    result->fully_resident_compute_ready =
        solver_codekd_test_fully_resident() == 0;
    result->empty_owner_replay =
        solver_codekd_test_detached_empty() == 0 &&
        solver_codekd_test_fully_resident_empty() == 0;
    result->callback_error =
        solver_codekd_test_detached_callback_error() == 0;
    result->eagain_retry =
        solver_codekd_test_detached_eagain() == 0;
    result->refusal_owner_fallback =
        solver_codekd_test_detached_refusal() == 0;
    result->cancellation_stopped =
        solver_codekd_test_detached_cancellation() == 0;
    return 0;
}
