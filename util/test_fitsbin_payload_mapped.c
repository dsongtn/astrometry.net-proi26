/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"

void test_fitsbin_payload_deferred_mapped_plan(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_planned_state_t plan;
    int submitted;
    int waited = -1;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    memset(&plan, 0, sizeof(plan));
    plan.range.data = fixture.chunk->data;
    plan.range.size = sizeof(fixture.bytes);
    plan.result = 1;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));

    submitted = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_planned_ranges,
        &plan,
        1024U * 1024U,
        &ticket);
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, 1, submitted);
    CuAssert(ct, "deferred mapped plan failed", waited > 0);
    CuAssertIntEquals(ct, 1, plan.calls);

    submitted = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_planned_ranges,
        &plan,
        1024U * 1024U,
        &ticket);
    waited = -1;
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 1, waited);
    CuAssertIntEquals(ct, 2, plan.calls);

    plan.result = 0;
    submitted = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_planned_ranges,
        &plan,
        1024U * 1024U,
        &ticket);
    waited = -1;
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 1, waited);
    CuAssertIntEquals(ct, 3, plan.calls);
#else
    CuAssertIntEquals(ct, 0, submitted);
    CuAssertIntEquals(ct, -1, waited);
    CuAssertIntEquals(ct, 0, plan.calls);
#endif

    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, 1, (int)stats.warm_calls);
    CuAssert(ct, "deferred mapped plan reported no pages",
             stats.warm_bytes > 0U);
    CuAssert(ct, "mapped completion was not reused",
             stats.cache_hits > 0U);
    CuAssert(ct, "mapped completion recorded no first miss",
             stats.cache_misses > 0U);
    CuAssertIntEquals(ct, 1, (int)stats.cache_allocations);
#else
    CuAssertIntEquals(ct, 0, (int)stats.warm_calls);
#endif
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssert(
        ct,
        "deferred mapped plan changed payload bytes",
        !memcmp((const unsigned char*)fixture.chunk->data,
                fixture.bytes,
                sizeof(fixture.bytes)));

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_queue_gap_coalescing(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_prefetch_range_t ranges[2];
    fitsbin_prefetch_range_t gap_range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    size_t page_size;
    size_t first_offset;
    size_t first_request = 0U;
    uintptr_t data_address;
    long detected_page_size;
    int submitted;
    int waited = -1;
    int calls = 0;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    detected_page_size = sysconf(_SC_PAGESIZE);
    CuAssert(ct, "invalid test page size", detected_page_size > 0);
    page_size = (size_t)detected_page_size;
    data_address = (uintptr_t)fixture.chunk->data;
    first_offset = (size_t)(data_address % page_size);
    if (first_offset) {
        first_offset = page_size - first_offset;
    }
    CuAssert(
        ct,
        "payload fixture is too small for a two-page queue gap",
        sizeof(fixture.bytes) >
            first_offset + 18U * page_size);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_set_mmap_advice(
            fixture.fitsbin,
            FITSBIN_MMAP_ADVICE_RANDOM,
            TRUE));
    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));

    ranges[0].data =
        (const unsigned char*)fixture.chunk->data + first_offset;
    ranges[0].size = 8U * page_size;
    ranges[1].data =
        (const unsigned char*)fixture.chunk->data +
        first_offset + 10U * page_size;
    ranges[1].size = 8U * page_size;
    payload_readahead_reset(1);
    submitted = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        ranges,
        2U,
        20U * page_size,
        &ticket);
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    pthread_mutex_lock(&payload_readahead.mutex);
    calls = payload_readahead.calls;
    first_request = payload_readahead.requests[0];
    pthread_mutex_unlock(&payload_readahead.mutex);

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, submitted);
    CuAssert(ct, "queue-gap mapped population failed", waited > 0);
    CuAssertIntEquals(ct, 1, calls);
    CuAssert(
        ct,
        "file-offset queue did not cover the bounded gap",
        first_request == 18U * page_size);

    /*
     * The storage request covered this page, but mapped completion did not.
     * A later exact request must still pass through the population barrier.
     */
    gap_range.data =
        (const unsigned char*)fixture.chunk->data +
        first_offset + 8U * page_size + 1U;
    gap_range.size = 1U;
    payload_readahead_reset(1);
    submitted = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &gap_range,
        1U,
        2U * page_size,
        &ticket);
    waited = -1;
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, submitted);
    CuAssert(ct, "queue-gap exact population failed", waited > 0);
#else
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_UNAVAILABLE, submitted);
    CuAssertIntEquals(ct, -1, waited);
#endif

    payload_readahead_reset(0);
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    CuAssert(
        ct,
        "queue-gap preparation changed payload bytes",
        !memcmp((const unsigned char*)fixture.chunk->data,
                fixture.bytes,
                sizeof(fixture.bytes)));
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_exact_range_order(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t ranges[2];
    unsigned char first[11];
    unsigned char second[17];
    unsigned char overlap[8];
    off_t offsets[2];
    int calls;
    int rc;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    ranges[0].data =
        (const unsigned char*)fixture.chunk->data + 31U;
    ranges[0].size = sizeof(first);
    ranges[0].logical_size = sizeof(first);
    ranges[0].destination = first;
    ranges[1].data =
        (const unsigned char*)fixture.chunk->data + 3U;
    ranges[1].size = sizeof(second);
    ranges[1].logical_size = sizeof(second);
    ranges[1].destination = second;

    payload_wrapper_reset(PAYLOAD_WRAPPER_RECORD);
    rc = fitsbin_pread_mapped_ranges(
        fixture.fitsbin, ranges, 2U);
    pthread_mutex_lock(&payload_wrapper.mutex);
    calls = payload_wrapper.calls;
    offsets[0] = payload_wrapper.offsets[0];
    offsets[1] = payload_wrapper.offsets[1];
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    CuAssertIntEquals(ct, 0, rc);
    CuAssertIntEquals(ct, 2, calls);
    CuAssert(ct, "exact ranges were not ordered by file offset",
        offsets[0] == fixture.chunk->data_file_offset + 3 &&
        offsets[1] == fixture.chunk->data_file_offset + 31);
    CuAssert(ct, "first ordered destination has wrong bytes",
        !memcmp(first, fixture.bytes + 31U, sizeof(first)));
    CuAssert(ct, "second ordered destination has wrong bytes",
        !memcmp(second, fixture.bytes + 3U, sizeof(second)));

    ranges[0].data =
        (const unsigned char*)fixture.chunk->data + 31U;
    ranges[0].size = sizeof(overlap);
    ranges[0].logical_size = sizeof(overlap);
    ranges[0].destination = overlap;
    ranges[1].data =
        (const unsigned char*)fixture.chunk->data + 3U;
    ranges[1].size = sizeof(overlap);
    ranges[1].logical_size = sizeof(overlap);
    ranges[1].destination = overlap;

    payload_wrapper_reset(PAYLOAD_WRAPPER_RECORD);
    rc = fitsbin_pread_mapped_ranges(
        fixture.fitsbin, ranges, 2U);
    pthread_mutex_lock(&payload_wrapper.mutex);
    calls = payload_wrapper.calls;
    offsets[0] = payload_wrapper.offsets[0];
    offsets[1] = payload_wrapper.offsets[1];
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    CuAssertIntEquals(ct, 0, rc);
    CuAssertIntEquals(ct, 2, calls);
    CuAssert(ct, "overlapping destinations changed caller order",
        offsets[0] == fixture.chunk->data_file_offset + 31 &&
        offsets[1] == fixture.chunk->data_file_offset + 3);
    CuAssert(ct, "overlapping destination changed final bytes",
        !memcmp(overlap, fixture.bytes + 3U, sizeof(overlap)));

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_async_direct_destination(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t ranges[2];
    fitsbin_prefetch_range_t warm;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    unsigned char first[11];
    unsigned char second[17];
    unsigned char untouched_first[sizeof(first)];
    unsigned char untouched_second[sizeof(second)];
    int refused;
    int refused_errno;
    int submitted;
    int configured;
    int async_warm;
    int sync_warm;
    int waited = -1;
    int calls;
    off_t offsets[2];

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    fitsbin_payload_set_thread_full_resident();
    configured = fitsbin_configure_index_mmap(fixture.fitsbin);
    fitsbin_payload_clear_thread_full_resident();
    CuAssertIntEquals(ct, 0, configured);
    CuAssert(
        ct,
        "full-resident source marker was not captured",
        fitsbin_payload_is_fully_resident(fixture.fitsbin));
    memset(first, 0xa5, sizeof(first));
    memset(second, 0xa5, sizeof(second));
    memcpy(untouched_first, first, sizeof(first));
    memcpy(untouched_second, second, sizeof(second));
    ranges[0].data =
        (const unsigned char*)fixture.chunk->data + 31U;
    ranges[0].size = sizeof(first);
    ranges[0].logical_size = 7U;
    ranges[0].destination = first;
    ranges[1].data =
        (const unsigned char*)fixture.chunk->data + 3U;
    ranges[1].size = sizeof(second);
    ranges[1].logical_size = 13U;
    ranges[1].destination = second;

    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    errno = 0;
    refused = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        ranges,
        2U,
        sizeof(first) + sizeof(second) - 1U,
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    refused_errno = errno;
    CuAssertIntEquals(ct, -1, refused);
    CuAssertIntEquals(ct, E2BIG, refused_errno);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssert(
        ct,
        "over-budget direct preparation wrote first destination",
        !memcmp(first, untouched_first, sizeof(first)));
    CuAssert(
        ct,
        "over-budget direct preparation wrote second destination",
        !memcmp(second, untouched_second, sizeof(second)));

    payload_wrapper_reset(PAYLOAD_WRAPPER_RECORD);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        ranges,
        2U,
        SIZE_MAX,
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    pthread_mutex_lock(&payload_wrapper.mutex);
    calls = payload_wrapper.calls;
    offsets[0] = payload_wrapper.offsets[0];
    offsets[1] = payload_wrapper.offsets[1];
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    warm.data = fixture.chunk->data;
    warm.size = sizeof(fixture.bytes);
    async_warm = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &warm,
        1U,
        SIZE_MAX,
        &ticket);
    sync_warm = fitsbin_prefetch_ranges(
        fixture.fitsbin,
        &warm,
        1U,
        SIZE_MAX);
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);

    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 2, waited);
    CuAssertIntEquals(ct, 0, async_warm);
    CuAssertIntEquals(ct, 0, sync_warm);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssertIntEquals(ct, 2, calls);
    CuAssert(ct, "async exact ranges were not ordered by file offset",
        offsets[0] == fixture.chunk->data_file_offset + 3 &&
        offsets[1] == fixture.chunk->data_file_offset + 31);
    CuAssert(
        ct,
        "first async direct destination has wrong bytes",
        !memcmp(first, fixture.bytes + 31U, sizeof(first)));
    CuAssert(
        ct,
        "second async direct destination has wrong bytes",
        !memcmp(second, fixture.bytes + 3U, sizeof(second)));
    CuAssertIntEquals(ct, 1, (int)stats.read_batches);
    CuAssertIntEquals(ct, 2, (int)stats.read_calls);
    CuAssertIntEquals(
        ct,
        (int)(sizeof(first) + sizeof(second)),
        (int)stats.read_bytes);
    CuAssertIntEquals(ct, 20, (int)stats.read_logical_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.warm_calls);
    CuAssertIntEquals(ct, 0, (int)stats.failures);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_async_mapping_boundary(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    unsigned char destination[64];
    unsigned char expected[64];
    size_t request_size;
    size_t data_map_offset;
    off_t expected_offset;
    off_t observed_offset;
    int submitted;
    int waited = 0;
    int calls;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    CuAssert(ct, "fixture has no file-backed mapping",
             fixture.chunk->map && fixture.chunk->mapsize);
    request_size = MIN(
        sizeof(destination), fixture.chunk->mapsize);
    CuAssert(ct, "fixture mapping is empty", request_size > 0U);
    data_map_offset =
        (size_t)((const unsigned char*)fixture.chunk->data -
                 (const unsigned char*)fixture.chunk->map);
    CuAssert(
        ct,
        "fixture mapping offset exceeds file payload offset",
        fixture.chunk->data_file_offset >=
            (off_t)data_map_offset);
    expected_offset =
        fixture.chunk->data_file_offset -
            (off_t)data_map_offset;
    memcpy(expected, fixture.chunk->map, request_size);
    memset(destination, 0, sizeof(destination));
    range.data = fixture.chunk->map;
    range.size = request_size;
    range.logical_size = request_size;
    range.destination = destination;

    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_RECORD);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        request_size,
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
    pthread_mutex_lock(&payload_wrapper.mutex);
    calls = payload_wrapper.calls;
    observed_offset = payload_wrapper.offsets[0];
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    fitsbin_payload_io_service_stop();

    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 1, waited);
    CuAssertIntEquals(ct, 1, calls);
    CuAssert(
        ct,
        "mapping-boundary read used the wrong file offset",
        observed_offset == expected_offset);
    CuAssert(
        ct,
        "mapping-boundary direct read returned wrong bytes",
        !memcmp(destination, expected, request_size));

    payload_fixture_close(&fixture);
}
