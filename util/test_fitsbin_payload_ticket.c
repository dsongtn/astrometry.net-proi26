/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"

void test_fitsbin_payload_poll_transfers_ticket_ownership(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    fitsbin_payload_io_stats_t repeated_stats;
    unsigned char destination[32];
    int submitted;
    int poll_status;
    int result = 0;
    int result_errno;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    memset(destination, 0, sizeof(destination));
    range.data = (const unsigned char*)fixture.chunk->data + 17U;
    range.size = sizeof(destination);
    range.logical_size = sizeof(destination);
    range.destination = destination;

    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_BLOCK);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        sizeof(destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    CuAssertIntEquals(ct, 1, submitted);
    CuAssertPtrNotNull(ct, ticket);
    CuAssertIntEquals(
        ct, 0, payload_wrapper_wait_for_calls(1, 2));

    poll_status = fitsbin_payload_io_ticket_poll_and_destroy(
        fixture.fitsbin,
        &ticket,
        &result);
    CuAssertIntEquals(ct, 0, poll_status);
    CuAssertPtrNotNull(ct, ticket);

    payload_wrapper_release();
    fitsbin_payload_io_service_stop();
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    errno = 0;
    poll_status = fitsbin_payload_io_ticket_poll_and_destroy(
        fixture.fitsbin,
        &ticket,
        &result);
    result_errno = errno;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &repeated_stats);
    fitsbin_payload_io_configure_workers(1);

    CuAssertIntEquals(ct, 1, poll_status);
    CuAssertIntEquals(ct, 1, result);
    CuAssertIntEquals(ct, 0, result_errno);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssert(
        ct,
        "poll-and-destroy returned incorrect direct-read bytes",
        !memcmp(
            destination,
            fixture.bytes + 17U,
            sizeof(destination)));
    CuAssertIntEquals(ct, 1, (int)stats.read_batches);
    CuAssertIntEquals(ct, 1, (int)stats.read_calls);
    CuAssertIntEquals(
        ct, (int)sizeof(destination), (int)stats.read_bytes);
    CuAssertIntEquals(
        ct, (int)sizeof(destination),
        (int)stats.read_logical_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.failures);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_poll_failure_transfers_once(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    fitsbin_payload_io_stats_t repeated_stats;
    unsigned char destination[32];
    int submitted;
    int poll_status;
    int result = 0;
    int result_errno;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    memset(destination, 0, sizeof(destination));
    range.data = (const unsigned char*)fixture.chunk->data + 23U;
    range.size = sizeof(destination);
    range.logical_size = sizeof(destination);
    range.destination = destination;

    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_EOF);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        sizeof(destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    fitsbin_payload_io_service_stop();
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    errno = 0;
    poll_status = fitsbin_payload_io_ticket_poll_and_destroy(
        fixture.fitsbin,
        &ticket,
        &result);
    result_errno = errno;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &repeated_stats);
    fitsbin_payload_io_configure_workers(1);

    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 1, poll_status);
    CuAssertIntEquals(ct, -1, result);
    CuAssertIntEquals(ct, EIO, result_errno);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssertIntEquals(ct, 0, (int)stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)stats.read_calls);
    CuAssertIntEquals(ct, 1, (int)stats.failures);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.failures);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_poll_cancel_transfers_once(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    fitsbin_payload_io_stats_t repeated_stats;
    unsigned char destination[32];
    int submitted;
    int cancelled;
    int poll_status;
    int result = -1;
    int result_errno;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    memset(destination, 0, sizeof(destination));
    range.data = (const unsigned char*)fixture.chunk->data + 29U;
    range.size = sizeof(destination);
    range.logical_size = sizeof(destination);
    range.destination = destination;

    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_BLOCK);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        sizeof(destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    CuAssertIntEquals(ct, 1, submitted);
    CuAssertPtrNotNull(ct, ticket);
    CuAssertIntEquals(
        ct, 0, payload_wrapper_wait_for_calls(1, 2));
    cancelled = fitsbin_payload_io_ticket_cancel_async(ticket);
    payload_wrapper_release();
    fitsbin_payload_io_service_stop();
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    errno = 0;
    poll_status = fitsbin_payload_io_ticket_poll_and_destroy(
        fixture.fitsbin,
        &ticket,
        &result);
    result_errno = errno;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &repeated_stats);
    fitsbin_payload_io_configure_workers(1);

    CuAssertIntEquals(ct, 1, cancelled);
    CuAssertIntEquals(ct, 1, poll_status);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(ct, ECANCELED, result_errno);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssertIntEquals(ct, 0, (int)stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.failures);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_drain_registered_waiter(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    fitsbin_payload_io_stats_t repeated_stats;
    payload_wait_helper_state_t helper = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0
    };
    payload_ticket_wait_state_t waiter;
    payload_ticket_drain_state_t drain;
    pthread_t waiter_thread;
    pthread_t drain_thread;
    unsigned char destination[32];
    int submitted;
    int provider_wait_status;
    int waiter_created = 0;
    int waiter_registered = -1;
    int busy_poll_status = 0;
    int busy_poll_errno = 0;
    int busy_poll_result = 0;
    int drain_mutex_status;
    int drain_condition_status = -1;
    int drain_created = 0;
    int drain_started = -1;
    int cleanup_result = 0;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    memset(&waiter, 0, sizeof(waiter));
    memset(&drain, 0, sizeof(drain));
    memset(destination, 0, sizeof(destination));
    range.data = (const unsigned char*)fixture.chunk->data + 37U;
    range.size = sizeof(destination);
    range.logical_size = sizeof(destination);
    range.destination = destination;

    fitsbin_payload_io_configure_workers(1);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));
    payload_wrapper_reset(PAYLOAD_WRAPPER_BLOCK);
    submitted = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        sizeof(destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &ticket);
    provider_wait_status =
        payload_wrapper_wait_for_calls(1, 2);

    waiter.fitsbin = fixture.fitsbin;
    waiter.ticket = ticket;
    waiter.helper = &helper;
    waiter.result = -1;
    if (ticket &&
        !pthread_create(
            &waiter_thread,
            NULL,
            payload_ticket_wait_thread,
            &waiter)) {
        waiter_created = 1;
        waiter_registered =
            payload_wait_helper_wait_for_calls(
                &helper, 1, 2);
    }
    if (!waiter_registered) {
        errno = 0;
        busy_poll_status =
            fitsbin_payload_io_ticket_poll_and_destroy(
                fixture.fitsbin,
                &ticket,
                &busy_poll_result);
        busy_poll_errno = errno;
    }

    drain_mutex_status = pthread_mutex_init(&drain.mutex, NULL);
    if (!drain_mutex_status) {
        drain_condition_status =
            pthread_cond_init(&drain.condition, NULL);
    }
    drain.fitsbin = fixture.fitsbin;
    drain.ticket = ticket;
    drain.status = -1;
    drain.result = -1;
    if (!waiter_registered &&
        !drain_mutex_status &&
        !drain_condition_status &&
        !pthread_create(
            &drain_thread,
            NULL,
            payload_ticket_drain_thread,
            &drain)) {
        drain_created = 1;
        drain_started =
            payload_ticket_drain_wait_started(&drain, 2);
    }

    payload_wrapper_release();
    if (drain_created) {
        pthread_join(drain_thread, NULL);
        ticket = drain.ticket;
    }
    if (waiter_created) {
        pthread_join(waiter_thread, NULL);
    }
    fitsbin_payload_io_service_stop();
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    if (ticket) {
        cleanup_result =
            fitsbin_payload_io_ticket_drain_and_destroy(
                fixture.fitsbin,
                &ticket,
                &busy_poll_result);
    }
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &repeated_stats);
    fitsbin_payload_io_configure_workers(1);

    CuAssertIntEquals(ct, 1, submitted);
    CuAssertIntEquals(ct, 0, provider_wait_status);
    CuAssertIntEquals(ct, 1, waiter_created);
    CuAssertIntEquals(ct, 0, waiter_registered);
    CuAssertIntEquals(ct, -1, busy_poll_status);
    CuAssertIntEquals(ct, EBUSY, busy_poll_errno);
    CuAssertIntEquals(ct, 0, busy_poll_result);
    CuAssertIntEquals(ct, 0, drain_mutex_status);
    CuAssertIntEquals(ct, 0, drain_condition_status);
    CuAssertIntEquals(ct, 1, drain_created);
    CuAssertIntEquals(ct, 0, drain_started);
    CuAssertIntEquals(ct, 1, waiter.configured);
    CuAssertIntEquals(ct, 1, waiter.result);
    CuAssertIntEquals(ct, 0, waiter.error);
    CuAssertIntEquals(ct, 1, drain.status);
    CuAssertIntEquals(ct, 1, drain.result);
    CuAssertIntEquals(ct, 0, drain.error);
    CuAssertPtrEquals(ct, NULL, drain.ticket);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssertIntEquals(ct, 0, cleanup_result);
    CuAssert(
        ct,
        "drained waiter returned incorrect direct-read bytes",
        !memcmp(
            destination,
            fixture.bytes + 37U,
            sizeof(destination)));
    CuAssertIntEquals(ct, 1, (int)stats.read_batches);
    CuAssertIntEquals(ct, 1, (int)stats.read_calls);
    CuAssertIntEquals(
        ct, (int)sizeof(destination), (int)stats.read_bytes);
    CuAssertIntEquals(
        ct, (int)sizeof(destination),
        (int)stats.read_logical_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_batches);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.read_bytes);
    CuAssertIntEquals(
        ct, 0, (int)repeated_stats.read_logical_bytes);
    CuAssertIntEquals(
        ct, 0, (int)repeated_stats.wait_nanoseconds);
    CuAssertIntEquals(ct, 0, (int)repeated_stats.failures);

    if (!drain_condition_status) {
        pthread_cond_destroy(&drain.condition);
    }
    if (!drain_mutex_status) {
        pthread_mutex_destroy(&drain.mutex);
    }
    payload_fixture_close(&fixture);
}
