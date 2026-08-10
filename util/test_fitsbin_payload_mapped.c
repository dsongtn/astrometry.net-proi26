/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"
#include "fitsbin_internal.h"

void test_fitsbin_payload_mapped_plan_keeps_unrequested_gap(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_prefetch_range_t ranges[2];
    fitsbin_mapped_span_t spans[2];
    uintptr_t data_begin;
    uintptr_t first_page;
    size_t span_count = 0U;
    size_t byte_count = 0U;
    size_t logical_byte_count = 0U;
    size_t exact_span_count = 0U;
    size_t gap_count = 0U;
    size_t gap_bytes = 0U;
    unsigned long long page_count = 0ULL;
    unsigned long long reused_pages = 0ULL;
    long detected_page_size;
    size_t page_size;

    payload_fixture_open_for_test(ct, &fixture);
    detected_page_size = sysconf(_SC_PAGESIZE);
    CuAssert(ct, "failed to detect page size", detected_page_size > 0);
    page_size = (size_t)detected_page_size;
    data_begin = (uintptr_t)fixture.chunk->data;
    first_page = data_begin;
    if (first_page % (uintptr_t)page_size) {
        first_page += (uintptr_t)page_size -
            first_page % (uintptr_t)page_size;
    }
    CuAssert(
        ct,
        "payload fixture is too small for sparse mapped ranges",
        first_page + 5U * page_size <=
            data_begin + sizeof(fixture.bytes));
    memset(ranges, 0, sizeof(ranges));
    memset(spans, 0, sizeof(spans));
    ranges[0].data = (const void*)first_page;
    ranges[0].size = 2U * page_size;
    ranges[1].data = (const void*)(first_page + 3U * page_size);
    ranges[1].size = 2U * page_size;

    CuAssertIntEquals(
        ct,
        0,
        fitsbin_prepare_mapped_spans(
            fixture.fitsbin,
            ranges,
            2U,
            5U * page_size,
            0ULL,
            FALSE,
            spans,
            2U,
            &span_count,
            &byte_count,
            &logical_byte_count,
            &page_count,
            &exact_span_count,
            &reused_pages,
            &gap_count,
            &gap_bytes));
    CuAssertIntEquals(ct, 2, (int)span_count);
    CuAssertIntEquals(ct, 2, (int)exact_span_count);
    CuAssertIntEquals(
        ct, (int)(4U * page_size), (int)logical_byte_count);
    CuAssertIntEquals(ct, (int)(4U * page_size), (int)byte_count);
    CuAssertIntEquals(ct, 4, (int)page_count);
    CuAssertIntEquals(ct, 0, (int)reused_pages);
    CuAssertIntEquals(ct, 0, (int)gap_count);
    CuAssertIntEquals(ct, 0, (int)gap_bytes);
    CuAssert(
        ct,
        "sparse mapped ranges were not kept separate",
        spans[0].end < spans[1].begin);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_deferred_mapped_plan(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_planned_state_t plan;
    int submitted;
    int waited = -1;

    payload_fixture_open_for_test(ct, &fixture);
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
    CuAssert(ct, "deferred mapped plan reused no exact pages",
             stats.cache_hits > 0U);
    CuAssert(ct, "deferred mapped plan recorded no first miss",
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

void test_fitsbin_payload_fully_resident_plan_is_advisory(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_planned_state_t plan;
    int configured;
    int submitted;

    payload_fixture_open_for_test(ct, &fixture);
    fitsbin_payload_set_thread_full_resident();
    configured = fitsbin_configure_index_mmap(fixture.fitsbin);
    fitsbin_payload_clear_thread_full_resident();
    CuAssertIntEquals(ct, 0, configured);
    CuAssert(
        ct,
        "full-resident source marker was not captured",
        fitsbin_payload_is_fully_resident(fixture.fitsbin));
    memset(&plan, 0, sizeof(plan));
    plan.range.data = fixture.chunk->data;
    plan.range.size = sizeof(fixture.bytes);
    plan.result = 1;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_service_stop();

    errno = EBUSY;
    submitted = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_planned_ranges,
        &plan,
        sizeof(fixture.bytes),
        &ticket);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);

    CuAssertIntEquals(ct, 0, submitted);
    CuAssertIntEquals(ct, 0, errno);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssertIntEquals(ct, 0, plan.calls);
    CuAssertIntEquals(ct, 0, (int)stats.warm_calls);
    CuAssertIntEquals(ct, 0, (int)stats.warm_ranges);
    CuAssertIntEquals(ct, 0, (int)stats.warm_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssert(
        ct,
        "full-resident advisory path changed payload bytes",
        !memcmp((const unsigned char*)fixture.chunk->data,
                fixture.bytes,
                sizeof(fixture.bytes)));

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_mapped_source_parallel_lanes(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* first = NULL;
    fitsbin_payload_io_ticket_t* second = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_planned_gate_t gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_blocking_planned_state_t plans[2];
    int first_submit;
    int second_submit = -1;
    int first_started;
    int second_started_concurrently = -1;
    int first_wait = -1;
    int second_wait = -1;
    int condition_status;
    int mutex_status;

    payload_fixture_open_for_test(ct, &fixture);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    memset(plans, 0, sizeof(plans));
    plans[0].gate = &gate;
    plans[0].range.data = fixture.chunk->data;
    plans[0].range.size = sizeof(fixture.bytes);
    plans[1] = plans[0];
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(2));

    first_submit = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_blocking_planned_ranges,
        &plans[0],
        1024U * 1024U,
        &first);
    first_started = payload_planned_gate_wait_for_calls(
        &gate, 1, 2);
    if (!first_started) {
        second_submit = fitsbin_prefetch_ranges_planned_submit(
            fixture.fitsbin,
            payload_blocking_planned_ranges,
            &plans[1],
            1024U * 1024U,
            &second);
        second_started_concurrently = payload_planned_gate_wait_for_calls(
            &gate, 2, 1);
    }
    payload_planned_gate_release(&gate);
    if (first) {
        first_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, first);
        fitsbin_payload_io_ticket_destroy(first);
        first = NULL;
    }
    if (second) {
        second_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, second);
        fitsbin_payload_io_ticket_destroy(second);
        second = NULL;
    }
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    condition_status = pthread_cond_destroy(&gate.condition);
    mutex_status = pthread_mutex_destroy(&gate.mutex);

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, 1, first_submit);
    CuAssertIntEquals(ct, 0, first_started);
    CuAssertIntEquals(ct, 1, second_submit);
    CuAssertIntEquals(ct, 0, second_started_concurrently);
    CuAssert(ct, "first mapped source ticket failed", first_wait > 0);
    CuAssert(ct, "second mapped source ticket failed", second_wait > 0);
    CuAssertIntEquals(ct, 1, plans[0].calls);
    CuAssertIntEquals(ct, 1, plans[1].calls);
    CuAssert(ct, "parallel mapped tickets issued invalid fill count",
             stats.warm_calls >= 1U && stats.warm_calls <= 2U);
    CuAssert(ct, "parallel mapped tickets recorded no first miss",
             stats.cache_misses > 0U);
#else
    CuAssertIntEquals(ct, 0, first_submit);
#endif
    CuAssertIntEquals(ct, 0, condition_status);
    CuAssertIntEquals(ct, 0, mutex_status);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_mapped_completion_precedes_next_plan(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* blocker = NULL;
    fitsbin_payload_io_ticket_t* first = NULL;
    fitsbin_payload_io_ticket_t* second = NULL;
    payload_planned_gate_t blocker_gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_planned_gate_t second_gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_blocking_planned_state_t blocker_plan;
    payload_planned_state_t first_plan;
    payload_blocking_planned_state_t second_plan;
    int blocker_submit;
    int first_submit = -1;
    int second_submit = -1;
    int blocker_started;
    int second_started = -1;
    int first_poll = -1;
    int first_result = -1;
    int first_cancel = -1;
    int blocker_wait = -1;
    int first_wait = -1;
    int second_wait = -1;
    int blocker_condition_status;
    int blocker_mutex_status;
    int second_condition_status;
    int second_mutex_status;

    payload_fixture_open_for_test(ct, &fixture);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    memset(&blocker_plan, 0, sizeof(blocker_plan));
    memset(&first_plan, 0, sizeof(first_plan));
    memset(&second_plan, 0, sizeof(second_plan));
    blocker_plan.gate = &blocker_gate;
    blocker_plan.range.data = fixture.chunk->data;
    blocker_plan.range.size = sizeof(fixture.bytes);
    first_plan.range = blocker_plan.range;
    first_plan.result = 1;
    second_plan.gate = &second_gate;
    second_plan.range = blocker_plan.range;
    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));

    blocker_submit = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_blocking_planned_ranges,
        &blocker_plan,
        1024U * 1024U,
        &blocker);
    blocker_started = payload_planned_gate_wait_for_calls(
        &blocker_gate, 1, 2);
    if (!blocker_started) {
        first_submit = fitsbin_prefetch_ranges_planned_submit(
            fixture.fitsbin,
            payload_planned_ranges,
            &first_plan,
            1024U * 1024U,
            &first);
        second_submit = fitsbin_prefetch_ranges_planned_submit(
            fixture.fitsbin,
            payload_blocking_planned_ranges,
            &second_plan,
            1024U * 1024U,
            &second);
    }
    payload_planned_gate_release(&blocker_gate);
    if (!blocker_started && first && second) {
        second_started = payload_planned_gate_wait_for_calls(
            &second_gate, 1, 2);
        if (!second_started) {
            first_poll = fitsbin_payload_io_ticket_poll(
                fixture.fitsbin, first, &first_result);
            first_cancel =
                fitsbin_payload_io_ticket_cancel_async(first);
        }
    }
    payload_planned_gate_release(&second_gate);
    if (blocker) {
        blocker_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, blocker);
        fitsbin_payload_io_ticket_destroy(blocker);
        blocker = NULL;
    }
    if (first) {
        first_wait = first_poll == 1
            ? first_result
            : fitsbin_payload_io_ticket_wait(
                  fixture.fitsbin, first);
        fitsbin_payload_io_ticket_destroy(first);
        first = NULL;
    }
    if (second) {
        second_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, second);
        fitsbin_payload_io_ticket_destroy(second);
        second = NULL;
    }
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    blocker_condition_status =
        pthread_cond_destroy(&blocker_gate.condition);
    blocker_mutex_status = pthread_mutex_destroy(&blocker_gate.mutex);
    second_condition_status =
        pthread_cond_destroy(&second_gate.condition);
    second_mutex_status = pthread_mutex_destroy(&second_gate.mutex);

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, blocker_submit);
    CuAssertIntEquals(ct, 0, blocker_started);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, first_submit);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, second_submit);
    CuAssertIntEquals(ct, 0, second_started);
    CuAssertIntEquals(ct, 1, first_poll);
    CuAssertIntEquals(ct, 1, first_result);
    CuAssertIntEquals(ct, 0, first_cancel);
    CuAssert(ct, "blocking mapped ticket failed", blocker_wait > 0);
    CuAssert(ct, "inline mapped ticket failed", first_wait > 0);
    CuAssert(ct, "second mapped ticket failed", second_wait > 0);
    CuAssertIntEquals(ct, 1, blocker_plan.calls);
    CuAssertIntEquals(ct, 1, first_plan.calls);
    CuAssertIntEquals(ct, 1, second_plan.calls);
#else
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_UNAVAILABLE, blocker_submit);
#endif
    CuAssertIntEquals(ct, 0, blocker_condition_status);
    CuAssertIntEquals(ct, 0, blocker_mutex_status);
    CuAssertIntEquals(ct, 0, second_condition_status);
    CuAssertIntEquals(ct, 0, second_mutex_status);
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_precomputed_plan_refresh(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_ticket_t* first = NULL;
    fitsbin_payload_io_ticket_t* second = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_planned_gate_t gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_blocking_planned_state_t first_plan;
    fitsbin_prefetch_range_t second_range;
    int first_submit;
    int second_submit = -1;
    int first_started;
    int second_poll = -1;
    int second_result = -1;
    int first_wait = -1;
    int second_wait = -1;
    int condition_status;
    int mutex_status;

    payload_fixture_open_for_test(ct, &fixture);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    memset(&first_plan, 0, sizeof(first_plan));
    first_plan.gate = &gate;
    first_plan.range.data = fixture.chunk->data;
    first_plan.range.size = sizeof(fixture.bytes);
    second_range = first_plan.range;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(2));

    first_submit = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_blocking_planned_ranges,
        &first_plan,
        1024U * 1024U,
        &first);
    first_started = payload_planned_gate_wait_for_calls(
        &gate, 1, 2);
    if (!first_started) {
        second_submit = fitsbin_prefetch_ranges_submit(
            fixture.fitsbin,
            &second_range,
            1U,
            1024U * 1024U,
            &second);
        if (second) {
            second_poll = fitsbin_payload_io_ticket_poll(
                fixture.fitsbin, second, &second_result);
            second_wait = fitsbin_payload_io_ticket_wait(
                fixture.fitsbin, second);
            fitsbin_payload_io_ticket_destroy(second);
            second = NULL;
        }
    }
    payload_planned_gate_release(&gate);
    if (first) {
        first_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, first);
        fitsbin_payload_io_ticket_destroy(first);
        first = NULL;
    }
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    condition_status = pthread_cond_destroy(&gate.condition);
    mutex_status = pthread_mutex_destroy(&gate.mutex);

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, 1, first_submit);
    CuAssertIntEquals(ct, 0, first_started);
    CuAssertIntEquals(ct, 1, second_submit);
    CuAssertIntEquals(ct, 0, second_poll);
    CuAssert(ct, "first mapped source ticket failed", first_wait > 0);
    CuAssertIntEquals(ct, 1, second_wait);
    CuAssertIntEquals(ct, 1, first_plan.calls);
    CuAssertIntEquals(ct, 1, (int)stats.warm_calls);
    CuAssert(ct, "precomputed refresh reused no exact pages",
             stats.cache_hits > 0U);
#else
    CuAssertIntEquals(ct, 0, first_submit);
#endif
    CuAssertIntEquals(ct, 0, condition_status);
    CuAssertIntEquals(ct, 0, mutex_status);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_precomputed_plan_revalidates_under_queue_pressure(
    CuTest* ct) {
    payload_fixture_t fixture;

    payload_fixture_open_for_test(ct, &fixture);
#if defined(MADV_POPULATE_READ)
    fitsbin_prefetch_range_t target_ranges[2];
    fitsbin_pread_range_t demand_range;
    fitsbin_payload_io_ticket_t* warm_ticket = NULL;
    fitsbin_payload_io_ticket_t* probe_ticket = NULL;
    fitsbin_payload_io_ticket_t* blocker_ticket = NULL;
    fitsbin_payload_io_ticket_t* target_ticket = NULL;
    fitsbin_payload_io_ticket_t* demand_tickets[
        FITSBIN_PAYLOAD_IO_MAX_JOBS - 2U];
    fitsbin_payload_io_stats_t stats;
    payload_planned_gate_t gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_blocking_planned_state_t blocker_plan;
    unsigned char demand_destinations[
        FITSBIN_PAYLOAD_IO_MAX_JOBS - 2U];
    const void* first_cover = NULL;
    const void* second_cover = NULL;
    size_t first_cover_size = 0U;
    size_t second_cover_size = 0U;
    size_t first_exact_offset = 0U;
    size_t second_exact_offset = 0U;
    off_t first_file_offset = 0;
    off_t second_file_offset = 0;
    size_t demand_index;
    size_t demand_submitted = 0U;
    int demand_status = 0;
    int demand_waits_ok = 1;
    int warm_submit;
    int warm_wait = -1;
    int probe_submit;
    int probe_ticket_present;
    int blocker_submit;
    int blocker_started = -1;
    int blocker_wait = -1;
    int target_submit = -1;
    int target_ticket_present = 0;
    int target_poll = -1;
    int target_poll_result = -1;
    int target_wait = -1;
    int condition_status;
    int mutex_status;

    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    memset(&blocker_plan, 0, sizeof(blocker_plan));
    memset(demand_tickets, 0, sizeof(demand_tickets));
    memset(demand_destinations, 0, sizeof(demand_destinations));
    target_ranges[0].data = fixture.chunk->data;
    target_ranges[0].size = 1U;
    target_ranges[1].data =
        (const unsigned char*)fixture.chunk->data +
        sizeof(fixture.bytes) - 1U;
    target_ranges[1].size = 1U;
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_mapped_range_page_cover(
            fixture.fitsbin,
            target_ranges[0].data,
            target_ranges[0].size,
            &first_cover,
            &first_cover_size,
            &first_file_offset,
            &first_exact_offset));
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_mapped_range_page_cover(
            fixture.fitsbin,
            target_ranges[1].data,
            target_ranges[1].size,
            &second_cover,
            &second_cover_size,
            &second_file_offset,
            &second_exact_offset));
    CuAssert(
        ct,
        "refresh-expansion ranges share one mapped page",
        first_cover != second_cover ||
        first_file_offset != second_file_offset);
    CuAssert(ct, "first refresh cover is empty", first_cover_size > 0U);
    CuAssert(ct, "second refresh cover is empty", second_cover_size > 0U);
    demand_range.data = target_ranges[0].data;
    demand_range.size = 1U;
    demand_range.logical_size = 1U;
    demand_range.destination = demand_destinations;
    blocker_plan.gate = &gate;
    blocker_plan.range.data =
        (const unsigned char*)fixture.chunk->data +
        sizeof(fixture.bytes) / 2U;
    blocker_plan.range.size = 1U;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));

    warm_submit = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &target_ranges[0],
        1U,
        sizeof(fixture.bytes),
        &warm_ticket);
    if (warm_ticket) {
        warm_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, warm_ticket);
        fitsbin_payload_io_ticket_destroy(warm_ticket);
        warm_ticket = NULL;
    }
    probe_submit = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &target_ranges[0],
        1U,
        sizeof(fixture.bytes),
        &probe_ticket);
    probe_ticket_present = probe_ticket != NULL;
    if (probe_ticket) {
        (void)fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, probe_ticket);
        fitsbin_payload_io_ticket_destroy(probe_ticket);
        probe_ticket = NULL;
    }

    blocker_submit = fitsbin_prefetch_ranges_planned_submit(
        fixture.fitsbin,
        payload_blocking_planned_ranges,
        &blocker_plan,
        sizeof(fixture.bytes),
        &blocker_ticket);
    if (blocker_ticket) {
        blocker_started = payload_planned_gate_wait_for_calls(
            &gate, 1, 2);
    }
    if (!blocker_started) {
        target_submit = fitsbin_prefetch_ranges_submit(
            fixture.fitsbin,
            target_ranges,
            2U,
            sizeof(fixture.bytes),
            &target_ticket);
        target_ticket_present = target_ticket != NULL;
        if (target_ticket) {
            target_poll = fitsbin_payload_io_ticket_poll(
                fixture.fitsbin,
                target_ticket,
                &target_poll_result);
        }
    }
    if (!blocker_started && target_ticket) {
        for (demand_index = 0U;
             demand_index < FITSBIN_PAYLOAD_IO_MAX_JOBS - 2U;
             demand_index++) {
            demand_range.destination =
                &demand_destinations[demand_index];
            demand_status = fitsbin_pread_mapped_ranges_submit(
                fixture.fitsbin,
                &demand_range,
                1U,
                1U,
                FITSBIN_PAYLOAD_IO_PRIORITY_DEMAND,
                &demand_tickets[demand_index]);
            if (demand_status != FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED ||
                !demand_tickets[demand_index]) {
                break;
            }
            demand_submitted++;
        }
    }

    payload_planned_gate_release(&gate);
    if (blocker_ticket) {
        blocker_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, blocker_ticket);
        fitsbin_payload_io_ticket_destroy(blocker_ticket);
        blocker_ticket = NULL;
    }
    for (demand_index = 0U;
         demand_index < demand_submitted;
         demand_index++) {
        int demand_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, demand_tickets[demand_index]);

        if (demand_wait != 1 ||
            demand_destinations[demand_index] != fixture.bytes[0]) {
            demand_waits_ok = 0;
        }
        fitsbin_payload_io_ticket_destroy(
            demand_tickets[demand_index]);
        demand_tickets[demand_index] = NULL;
    }
    if (target_ticket) {
        target_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, target_ticket);
        fitsbin_payload_io_ticket_destroy(target_ticket);
        target_ticket = NULL;
    }
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    condition_status = pthread_cond_destroy(&gate.condition);
    mutex_status = pthread_mutex_destroy(&gate.mutex);

    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, warm_submit);
    CuAssertIntEquals(ct, 1, warm_wait);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, probe_submit);
    CuAssertIntEquals(ct, 1, probe_ticket_present);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, blocker_submit);
    CuAssertIntEquals(ct, 0, blocker_started);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, target_submit);
    CuAssertIntEquals(ct, 1, target_ticket_present);
    CuAssertIntEquals(ct, 0, target_poll);
    CuAssertIntEquals(ct, 0, target_poll_result);
    CuAssertIntEquals(
        ct,
        (int)(FITSBIN_PAYLOAD_IO_MAX_JOBS - 2U),
        (int)demand_submitted);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, demand_status);
    CuAssertIntEquals(ct, 1, demand_waits_ok);
    CuAssertIntEquals(ct, 1, blocker_wait);
    CuAssertIntEquals(ct, 2, target_wait);
    CuAssert(ct, "queued exact mapping reuse was not observed",
             stats.cache_hits > 0U);
    CuAssert(ct, "expired exact mapping pages were not repopulated",
             stats.cache_misses > 0U);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(ct, 0, condition_status);
    CuAssertIntEquals(ct, 0, mutex_status);
#endif
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

    payload_fixture_open_for_test(ct, &fixture);
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
    CuAssert(
        ct,
        "mapped population issued duplicate file-offset requests",
        calls == 0 || calls == 1);
    if (calls == 1) {
        CuAssert(
            ct,
            "file-offset fallback did not cover the bounded gap",
            first_request == 18U * page_size);
    }

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

    payload_fixture_open_for_test(ct, &fixture);
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
    int async_warm_ticket_present;
    int sync_warm;
    int waited = -1;
    int calls;
    off_t offsets[2];

    payload_fixture_open_for_test(ct, &fixture);
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
    async_warm_ticket_present = ticket != NULL;
    if (ticket) {
        (void)fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
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
    CuAssertIntEquals(ct, 0, async_warm_ticket_present);
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

    payload_fixture_open_for_test(ct, &fixture);
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
