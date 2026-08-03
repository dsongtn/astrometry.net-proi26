/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"

void test_fitsbin_payload_dynamic_service_width(CuTest* ct) {
    int start_status;
    int repeated_start_status;
    int restart_status;
    int running_width;
    int repeated_width;
    int restarted_width;
    int stopped_width;

    fitsbin_payload_io_configure_workers(6);
    start_status = fitsbin_payload_io_service_start(6);
    running_width = fitsbin_payload_io_service_width();
    repeated_start_status =
        fitsbin_payload_io_service_start(3);
    repeated_width = fitsbin_payload_io_service_width();
    if (!start_status) {
        fitsbin_payload_io_service_stop();
    }
    stopped_width = fitsbin_payload_io_service_width();
    fitsbin_payload_io_service_stop();
    restart_status = fitsbin_payload_io_service_start(8);
    restarted_width = fitsbin_payload_io_service_width();
    if (!restart_status) {
        fitsbin_payload_io_service_stop();
    }
    fitsbin_payload_io_configure_workers(1);

    CuAssertIntEquals(ct, 0, start_status);
    CuAssertIntEquals(ct, 6, running_width);
    CuAssertIntEquals(ct, 0, repeated_start_status);
    CuAssertIntEquals(ct, 6, repeated_width);
    CuAssertIntEquals(ct, 0, stopped_width);
    CuAssertIntEquals(ct, 0, restart_status);
    CuAssertIntEquals(ct, 8, restarted_width);
}

void test_fitsbin_payload_shared_reader_credit(CuTest* ct) {
    payload_fixture_t fixture;
    payload_credit_result_t two_readers;
    payload_credit_result_t one_reader;

    payload_fixture_open_for_test(ct, &fixture);

    CuAssertIntEquals(
        ct,
        0,
        payload_credit_round(
            &fixture,
            4,
            4,
            4,
            &two_readers));
    CuAssertIntEquals(
        ct, 4, two_readers.calls_before_release);
    CuAssertIntEquals(
        ct, 4, two_readers.active_before_release);
    CuAssertIntEquals(ct, 4, two_readers.max_active);
    CuAssertIntEquals(ct, 4, two_readers.total_calls);
    CuAssert(
        ct,
        "parallel payload calls did not share one descriptor",
        two_readers.first_fd >= 0 &&
        !two_readers.fd_mismatch);
    CuAssert(
        ct,
        "parallel payload reads returned wrong bytes",
        two_readers.all_reads_ok);

    CuAssertIntEquals(
        ct,
        0,
        payload_credit_round(
            &fixture,
            1,
            3,
            1,
            &one_reader));
    CuAssertIntEquals(
        ct, 1, one_reader.calls_before_release);
    CuAssertIntEquals(
        ct, 1, one_reader.active_before_release);
    CuAssertIntEquals(ct, 1, one_reader.max_active);
    CuAssertIntEquals(ct, 3, one_reader.total_calls);
    CuAssert(
        ct,
        "serial payload calls did not share one descriptor",
        one_reader.first_fd >= 0 &&
        !one_reader.fd_mismatch);
    CuAssert(
        ct,
        "serial payload reads returned wrong bytes",
        one_reader.all_reads_ok);

    fitsbin_payload_io_configure_workers(1);
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_async_direct_overlap(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_pread_range_t first_range;
    fitsbin_pread_range_t second_range;
    fitsbin_payload_io_ticket_t* first = NULL;
    fitsbin_payload_io_ticket_t* second = NULL;
    fitsbin_payload_io_stats_t stats;
    payload_wait_helper_state_t helper = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0
    };
    payload_ticket_wait_state_t waiters[2];
    pthread_t waiter_threads[2];
    unsigned char first_destination[8];
    unsigned char second_destination[8];
    int first_submit;
    int second_submit;
    int wait_status;
    int helper_start_status = -1;
    int helper_wakeup_status = -1;
    int helper_condition_destroy_status;
    int helper_mutex_destroy_status;
    int active;
    int max_active;
    int waiter_count = 0;
    int first_wait = -1;
    int second_wait = -1;
    int status;
    int i;

    payload_fixture_open_for_test(ct, &fixture);
    memset(waiters, 0, sizeof(waiters));
    memset(first_destination, 0, sizeof(first_destination));
    memset(second_destination, 0, sizeof(second_destination));
    first_range.data = fixture.chunk->data;
    first_range.size = sizeof(first_destination);
    first_range.logical_size = sizeof(first_destination);
    first_range.destination = first_destination;
    second_range.data =
        (const unsigned char*)fixture.chunk->data +
        sizeof(first_destination);
    second_range.size = sizeof(second_destination);
    second_range.logical_size = sizeof(second_destination);
    second_range.destination = second_destination;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(2));

    payload_wrapper_reset(PAYLOAD_WRAPPER_BLOCK);
    first_submit = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &first_range,
        1U,
        sizeof(first_destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &first);
    second_submit = fitsbin_pread_mapped_ranges_submit(
        fixture.fitsbin,
        &second_range,
        1U,
        sizeof(second_destination),
        FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
        &second);
    wait_status = payload_wrapper_wait_for_calls(2, 3);
    pthread_mutex_lock(&payload_wrapper.mutex);
    active = payload_wrapper.active;
    max_active = payload_wrapper.max_active;
    pthread_mutex_unlock(&payload_wrapper.mutex);

    if (first && second) {
        waiters[0].fitsbin = fixture.fitsbin;
        waiters[0].ticket = first;
        waiters[0].helper = &helper;
        waiters[0].result = -1;
        waiters[1].fitsbin = fixture.fitsbin;
        waiters[1].ticket = second;
        waiters[1].helper = &helper;
        waiters[1].result = -1;
        status = pthread_create(
            &waiter_threads[0],
            NULL,
            payload_ticket_wait_thread,
            &waiters[0]);
        if (!status) {
            waiter_count = 1;
            status = pthread_create(
                &waiter_threads[1],
                NULL,
                payload_ticket_wait_thread,
                &waiters[1]);
            if (!status) {
                waiter_count = 2;
            }
        }
    }
    if (waiter_count == 2) {
        helper_start_status =
            payload_wait_helper_wait_for_calls(
                &helper, 2, 3);
        if (!helper_start_status) {
            fitsbin_payload_io_notify_wait_helpers();
            helper_wakeup_status =
                payload_wait_helper_wait_for_calls(
                    &helper, 4, 3);
        }
    }
    payload_wrapper_release();

    for (i = 0; i < waiter_count; i++) {
        pthread_join(waiter_threads[i], NULL);
    }
    if (waiter_count >= 1 && waiters[0].configured > 0) {
        first_wait = waiters[0].result;
    } else if (first) {
        first_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, first);
    }
    if (waiter_count >= 2 && waiters[1].configured > 0) {
        second_wait = waiters[1].result;
    } else if (second) {
        second_wait = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, second);
    }
    if (first) {
        fitsbin_payload_io_ticket_destroy(first);
    }
    if (second) {
        fitsbin_payload_io_ticket_destroy(second);
    }
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    helper_condition_destroy_status =
        pthread_cond_destroy(&helper.condition);
    helper_mutex_destroy_status =
        pthread_mutex_destroy(&helper.mutex);

    CuAssertIntEquals(ct, 1, first_submit);
    CuAssertIntEquals(ct, 1, second_submit);
    CuAssertIntEquals(ct, 0, wait_status);
    CuAssertIntEquals(ct, 2, active);
    CuAssertIntEquals(ct, 2, max_active);
    CuAssertIntEquals(ct, 2, waiter_count);
    CuAssertIntEquals(ct, 0, helper_start_status);
    CuAssertIntEquals(ct, 0, helper_wakeup_status);
    CuAssertIntEquals(ct, 1, waiters[0].configured);
    CuAssertIntEquals(ct, 1, waiters[1].configured);
    CuAssertIntEquals(ct, 1, first_wait);
    CuAssertIntEquals(ct, 1, second_wait);
    CuAssertIntEquals(
        ct, 0, (int)fitsbin_payload_io_wait_helper_count());
    CuAssert(
        ct,
        "first async direct range returned wrong bytes",
        !memcmp(first_destination,
                fixture.bytes,
                sizeof(first_destination)));
    CuAssert(
        ct,
        "second async direct range returned wrong bytes",
        !memcmp(second_destination,
                fixture.bytes + sizeof(first_destination),
                sizeof(second_destination)));
    CuAssertIntEquals(ct, 2, (int)stats.read_batches);
    CuAssertIntEquals(ct, 2, (int)stats.read_calls);
    CuAssertIntEquals(
        ct,
        (int)(sizeof(first_destination) +
              sizeof(second_destination)),
        (int)stats.read_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(ct, 0, helper_condition_destroy_status);
    CuAssertIntEquals(ct, 0, helper_mutex_destroy_status);
    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_async_mapped_population(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_prefetch_range_t range;
    fitsbin_pread_range_t direct_range;
    fitsbin_payload_io_ticket_t* ticket = NULL;
    fitsbin_payload_io_ticket_t* reused_ticket = NULL;
    fitsbin_payload_io_stats_t stats;
    unsigned char direct_byte = 0U;
    int expired = -1;
    int expired_ticket_present = 0;
    int expired_waited = -1;
    int reused = -1;
    int reused_ticket_present = 0;
    int submitted;
    int waited = -1;
    unsigned int i;

    payload_fixture_open_for_test(ct, &fixture);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_configure_index_mmap(fixture.fitsbin));
    range.data = fixture.chunk->data;
    range.size = 1U;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    fitsbin_payload_io_configure_workers(2);
    CuAssertIntEquals(
        ct, 0, fitsbin_payload_io_service_start(1));

    submitted = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        SIZE_MAX,
        &ticket);
    if (ticket) {
        waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
#if defined(MADV_POPULATE_READ)
    reused = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        SIZE_MAX,
        &reused_ticket);
    reused_ticket_present = reused_ticket != NULL;
    if (reused_ticket) {
        (void)fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, reused_ticket);
        fitsbin_payload_io_ticket_destroy(reused_ticket);
        reused_ticket = NULL;
    }
    direct_range.data = fixture.chunk->data;
    direct_range.size = sizeof(direct_byte);
    direct_range.logical_size = sizeof(direct_byte);
    direct_range.destination = &direct_byte;
    for (i = 0U; i < PAYLOAD_COMPLETION_EXPIRY_TICKETS; i++) {
        int direct_submitted;
        int direct_waited = -1;

        direct_submitted = fitsbin_pread_mapped_ranges_submit(
            fixture.fitsbin,
            &direct_range,
            1U,
            sizeof(direct_byte),
            FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
            &ticket);
        if (ticket) {
            direct_waited = fitsbin_payload_io_ticket_wait(
                fixture.fitsbin, ticket);
            fitsbin_payload_io_ticket_destroy(ticket);
            ticket = NULL;
        }
        CuAssertIntEquals(ct, 1, direct_submitted);
        CuAssertIntEquals(ct, 1, direct_waited);
    }
    expired = fitsbin_prefetch_ranges_submit(
        fixture.fitsbin,
        &range,
        1U,
        SIZE_MAX,
        &ticket);
    expired_ticket_present = ticket != NULL;
    if (ticket) {
        expired_waited = fitsbin_payload_io_ticket_wait(
            fixture.fitsbin, ticket);
        fitsbin_payload_io_ticket_destroy(ticket);
        ticket = NULL;
    }
#endif
    fitsbin_payload_io_service_stop();
    fitsbin_payload_io_configure_workers(1);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);

#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, 1, submitted);
    CuAssert(ct, "async mapped population failed", waited > 0);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_READY, reused);
    CuAssertIntEquals(ct, 0, reused_ticket_present);
    CuAssertIntEquals(
        ct, FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED, expired);
    CuAssertIntEquals(ct, 1, expired_ticket_present);
    CuAssert(ct, "expired mapped population failed",
             expired_waited > 0);
    CuAssertPtrEquals(ct, NULL, ticket);
    CuAssert(ct, "mapped completion aging issued too few fills",
             stats.warm_calls >= 2U);
    CuAssert(ct, "async mapped population reported no pages",
             stats.warm_bytes > 0U);
#else
    CuAssertIntEquals(ct, 0, submitted);
    CuAssertIntEquals(ct, -1, waited);
    CuAssertIntEquals(ct, 0, (int)stats.warm_calls);
#endif
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssert(
        ct,
        "async mapped population changed payload bytes",
        !memcmp((const unsigned char*)fixture.chunk->data,
                fixture.bytes,
                sizeof(fixture.bytes)));

    payload_fixture_close(&fixture);
}
