/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"

void test_fitsbin_payload_short_read_and_eintr(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_stats_t stats;
    unsigned char destination[17];
    off_t expected_offset;
    int rc;
    int calls;
    off_t offsets[4];
    size_t requests[4];

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
    expected_offset =
        fixture.chunk->data_file_offset + 5;

    payload_wrapper_reset(PAYLOAD_WRAPPER_SHORT);
    errno = 0;
    rc = fitsbin_pread_mapped_range(
        fixture.fitsbin,
        (const unsigned char*)fixture.chunk->data + 5,
        sizeof(destination),
        destination);
    pthread_mutex_lock(&payload_wrapper.mutex);
    calls = payload_wrapper.calls;
    memcpy(offsets, payload_wrapper.offsets,
           sizeof(offsets));
    memcpy(requests, payload_wrapper.requests,
           sizeof(requests));
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    CuAssertIntEquals(ct, 0, rc);
    CuAssert(
        ct,
        "short-read accumulation returned wrong bytes",
        !memcmp(destination,
                fixture.bytes + 5,
                sizeof(destination)));

    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 1, (int)stats.read_calls);
    CuAssertIntEquals(
        ct,
        (int)sizeof(destination),
        (int)stats.read_bytes);
    CuAssertIntEquals(ct, 0, (int)stats.failures);

    CuAssertIntEquals(ct, 4, calls);
    CuAssert(
        ct,
        "EINTR retry changed the file offset",
        offsets[0] == expected_offset &&
        offsets[1] == expected_offset);
    CuAssert(
        ct,
        "positive short reads did not advance exact offsets",
        offsets[2] == expected_offset + 3 &&
        offsets[3] == expected_offset + 5);
    CuAssert(
        ct,
        "positive short reads did not reduce remaining sizes",
        requests[0] == sizeof(destination) &&
        requests[1] == sizeof(destination) &&
        requests[2] == sizeof(destination) - 3U &&
        requests[3] == sizeof(destination) - 5U);

    payload_fixture_close(&fixture);
}

void test_fitsbin_payload_eof_is_eio(CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_payload_io_stats_t stats;
    unsigned char destination[12];
    int rc;
    int saved_errno;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));

    payload_wrapper_reset(PAYLOAD_WRAPPER_EOF);
    errno = 0;
    rc = fitsbin_pread_mapped_range(
        fixture.fitsbin,
        (const unsigned char*)fixture.chunk->data + 9,
        sizeof(destination),
        destination);
    saved_errno = errno;
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    CuAssertIntEquals(ct, -1, rc);
    CuAssertIntEquals(ct, EIO, saved_errno);
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 0, (int)stats.read_calls);
    CuAssertIntEquals(ct, 1, (int)stats.failures);

    /*
     * A range-local read failure does not poison the separately opened
     * payload descriptor.  The next exact request can still complete.
     */
    rc = fitsbin_pread_mapped_range(
        fixture.fitsbin,
        (const unsigned char*)fixture.chunk->data + 9,
        sizeof(destination),
        destination);
    CuAssertIntEquals(ct, 0, rc);
    CuAssert(
        ct,
        "payload descriptor did not recover after EOF",
        !memcmp(destination,
                fixture.bytes + 9,
                sizeof(destination)));
    fitsbin_take_payload_io_stats(
        fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 1, (int)stats.read_calls);
    CuAssertIntEquals(ct, 0, (int)stats.failures);

    payload_fixture_close(&fixture);
}

void test_fitsbin_mapped_population_is_all_or_nothing(
    CuTest* ct) {
    payload_fixture_t fixture;
    fitsbin_prefetch_range_t ranges[2];
    fitsbin_payload_io_stats_t stats;
    int rc;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));
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
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fitsbin_get_mmap_advice(fixture.fitsbin));
    ranges[0].data =
        (const unsigned char*)fixture.chunk->data + 3;
    ranges[0].size = 11U;
    ranges[1].data =
        (const unsigned char*)fixture.chunk->data + 9;
    ranges[1].size = 17U;
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);

    rc = fitsbin_advise_mapped_ranges(
        fixture.fitsbin,
        ranges,
        sizeof(ranges) / sizeof(ranges[0]),
        1U);
    CuAssertIntEquals(ct, 0, rc);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 0, (int)stats.warm_calls);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fitsbin_get_mmap_advice(fixture.fitsbin));

#if defined(MADV_POPULATE_READ)
    rc = fitsbin_advise_mapped_ranges(
        fixture.fitsbin,
        ranges,
        sizeof(ranges) / sizeof(ranges[0]),
        SIZE_MAX);
    CuAssert(ct, "mapped plan was not populated", rc > 0);
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 1, (int)stats.warm_calls);
    CuAssert(ct, "mapped plan reported no pages",
             stats.warm_bytes > 0U);
    CuAssertIntEquals(ct, 0, (int)stats.failures);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fitsbin_get_mmap_advice(fixture.fitsbin));
#else
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_advise_mapped_ranges(
            fixture.fitsbin,
            ranges,
            sizeof(ranges) / sizeof(ranges[0]),
            SIZE_MAX));
#endif
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fitsbin_get_mmap_advice(fixture.fitsbin));

    CuAssert(
        ct,
        "mapped population changed payload bytes",
        !memcmp((const unsigned char*)fixture.chunk->data,
                fixture.bytes,
                sizeof(fixture.bytes)));

#if defined(MADV_POPULATE_READ)
    ranges[0].data = &rc;
    ranges[0].size = sizeof(rc);
    errno = 0;
    rc = fitsbin_advise_mapped_ranges(
        fixture.fitsbin,
        ranges,
        1U,
        SIZE_MAX);
    CuAssertIntEquals(ct, -1, rc);
    CuAssertIntEquals(ct, ERANGE, errno);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fitsbin_get_mmap_advice(fixture.fitsbin));
    fitsbin_take_payload_io_stats(fixture.fitsbin, &stats);
    CuAssertIntEquals(ct, 1, (int)stats.failures);
#endif
    payload_fixture_close(&fixture);
}
