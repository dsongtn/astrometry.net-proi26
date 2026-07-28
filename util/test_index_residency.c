/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "index_residency.h"

#include "cutest.h"

#define INDEX_RESIDENCY_TEST_BYTES 8193U

typedef struct index_residency_fixture {
    char path[128];
    unsigned char bytes[INDEX_RESIDENCY_TEST_BYTES];
    struct stat source_stat;
} index_residency_fixture_t;

static int index_residency_test_write_all(
    int fd,
    const unsigned char* data,
    size_t size) {
    while (size) {
        ssize_t written = write(fd, data, size);

        if (written < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (!written) {
            return -1;
        }
        data += (size_t)written;
        size -= (size_t)written;
    }
    return 0;
}

static int index_residency_test_read_all(
    int fd,
    unsigned char* data,
    size_t size) {
    while (size) {
        ssize_t count = read(fd, data, size);

        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (!count) {
            return -1;
        }
        data += (size_t)count;
        size -= (size_t)count;
    }
    return 0;
}

static int index_residency_fixture_open(
    index_residency_fixture_t* fixture) {
    static const char path_template[] =
        "/tmp/test-index-residency-XXXXXX";
    size_t i;
    int fd;

    if (!fixture ||
        sizeof(path_template) > sizeof(fixture->path)) {
        return -1;
    }
    memset(fixture, 0, sizeof(*fixture));
    memcpy(fixture->path, path_template,
           sizeof(path_template));
    for (i = 0U; i < sizeof(fixture->bytes); i++) {
        fixture->bytes[i] = (unsigned char)(
            (i * 37U + (i >> 3U)) & 0xffU);
    }
    fd = mkstemp(fixture->path);
    if (fd < 0) {
        return -1;
    }
    if (index_residency_test_write_all(
            fd,
            fixture->bytes,
            sizeof(fixture->bytes))) {
        close(fd);
        unlink(fixture->path);
        return -1;
    }
    if (close(fd)) {
        unlink(fixture->path);
        return -1;
    }
    if (stat(
            fixture->path,
            &fixture->source_stat)) {
        unlink(fixture->path);
        return -1;
    }
    return 0;
}

static void index_residency_fixture_close(
    index_residency_fixture_t* fixture) {
    if (!fixture || !fixture->path[0]) {
        return;
    }
    unlink(fixture->path);
    fixture->path[0] = '\0';
}

void test_index_residency_prepare_acquire(CuTest* ct) {
    index_residency_fixture_t fixture;
    index_residency_handle_t* first = NULL;
    index_residency_handle_t* second = NULL;
    index_residency_stats_t stats;
    index_residency_t* service = NULL;
    unsigned char observed[INDEX_RESIDENCY_TEST_BYTES];
    const struct stat* resident_source_stat;
    const char* first_path;
    const char* second_path;
    index_residency_result_t acquire_first;
    index_residency_result_t acquire_second;
    index_residency_result_t prepare;
    int backing_fd = -1;
    int backing_is_distinct = 0;
    int byte_match = 0;
    int drain_status = -1;
    int same_backing = 0;
    int source_stat_match = 0;
    int start_status;
    int stats_status;
    int stop_status;

    memset(&stats, 0, sizeof(stats));
    CuAssertIntEquals(
        ct, 0, index_residency_fixture_open(&fixture));
    start_status = index_residency_start(
        sizeof(fixture.bytes) * 2U, 1U, &service);
    if (start_status) {
        index_residency_fixture_close(&fixture);
        CuAssertIntEquals(ct, 0, start_status);
        return;
    }
    stats_status = index_residency_get_stats(
        service, &stats);
    if (stats_status) {
        stop_status = index_residency_stop(service);
        index_residency_fixture_close(&fixture);
        CuAssertIntEquals(ct, 0, stats_status);
        CuAssertIntEquals(ct, 0, stop_status);
        return;
    }

    if (!stats.backend_supported) {
        prepare = index_residency_prepare(
            service,
            fixture.path,
            INDEX_RESIDENCY_PRIORITY_LOOKAHEAD);
        acquire_first = index_residency_acquire(
            service, fixture.path, &first);
        stop_status = index_residency_stop(service);
        index_residency_fixture_close(&fixture);

        CuAssertIntEquals(
            ct, INDEX_RESIDENCY_FALLBACK, prepare);
        CuAssertIntEquals(
            ct, INDEX_RESIDENCY_FALLBACK,
            acquire_first);
        CuAssertPtrEquals(ct, NULL, first);
        CuAssertIntEquals(ct, 0, stop_status);
        return;
    }

    prepare = index_residency_prepare(
        service,
        fixture.path,
        INDEX_RESIDENCY_PRIORITY_LOOKAHEAD);
    drain_status = index_residency_drain(service);
    acquire_first = index_residency_acquire(
        service, fixture.path, &first);
    acquire_second = index_residency_acquire(
        service, fixture.path, &second);
    if (first && second) {
        first_path = index_residency_handle_path(first);
        second_path = index_residency_handle_path(second);
        resident_source_stat =
            index_residency_handle_source_stat(first);
        if (first_path && second_path) {
            same_backing = !strcmp(first_path, second_path);
            backing_is_distinct =
                strcmp(first_path, fixture.path) != 0;
            backing_fd = open(first_path, O_RDONLY);
        }
        if (resident_source_stat) {
            source_stat_match =
                resident_source_stat->st_dev ==
                    fixture.source_stat.st_dev &&
                resident_source_stat->st_ino ==
                    fixture.source_stat.st_ino &&
                resident_source_stat->st_size ==
                    fixture.source_stat.st_size;
        }
        if (backing_fd >= 0 &&
            !index_residency_test_read_all(
                backing_fd,
                observed,
                sizeof(observed))) {
            byte_match = !memcmp(
                observed,
                fixture.bytes,
                sizeof(observed));
        }
    }
    if (backing_fd >= 0) {
        close(backing_fd);
    }
    stats_status = index_residency_get_stats(
        service, &stats);
    index_residency_release(second);
    index_residency_release(first);
    stop_status = index_residency_stop(service);
    index_residency_fixture_close(&fixture);

    CuAssertIntEquals(
        ct, INDEX_RESIDENCY_ACCEPTED, prepare);
    CuAssertIntEquals(
        ct, INDEX_RESIDENCY_ACCEPTED,
        acquire_first);
    CuAssertIntEquals(
        ct, INDEX_RESIDENCY_ACCEPTED,
        acquire_second);
    CuAssert(ct, "resident handles did not share one backing",
             same_backing);
    CuAssert(ct, "resident backing reused the source path",
             backing_is_distinct);
    CuAssert(ct, "resident backing bytes differ from source",
             byte_match);
    CuAssert(ct, "resident source identity is incorrect",
             source_stat_match);
    CuAssertIntEquals(ct, 0, drain_status);
    CuAssertIntEquals(ct, 0, stats_status);
    CuAssertIntEquals(ct, 1, (int)stats.files_queued);
    CuAssertIntEquals(ct, 1, (int)stats.files_copied);
    CuAssertIntEquals(ct, 2, (int)stats.live_handles);
    CuAssert(
        ct,
        "resident demand did not deduplicate or hit cache",
        stats.loading_deduplications +
            stats.cache_hits >= 2U);
    CuAssertIntEquals(ct, 0, stop_status);
}

void test_index_residency_tight_budget(CuTest* ct) {
    index_residency_fixture_t fixture;
    index_residency_handle_t* handle = NULL;
    index_residency_stats_t stats;
    index_residency_t* service = NULL;
    index_residency_result_t acquire;
    index_residency_result_t prepare;
    int start_status;
    int stats_status;
    int stop_status;

    memset(&stats, 0, sizeof(stats));
    CuAssertIntEquals(
        ct, 0, index_residency_fixture_open(&fixture));
    start_status = index_residency_start(
        sizeof(fixture.bytes) - 1U, 1U, &service);
    if (start_status) {
        index_residency_fixture_close(&fixture);
        CuAssertIntEquals(ct, 0, start_status);
        return;
    }
    prepare = index_residency_prepare(
        service,
        fixture.path,
        INDEX_RESIDENCY_PRIORITY_SPECULATIVE);
    acquire = index_residency_acquire(
        service, fixture.path, &handle);
    stats_status = index_residency_get_stats(
        service, &stats);
    index_residency_release(handle);
    stop_status = index_residency_stop(service);
    index_residency_fixture_close(&fixture);

    CuAssertIntEquals(
        ct, INDEX_RESIDENCY_FALLBACK, prepare);
    CuAssertIntEquals(
        ct, INDEX_RESIDENCY_FALLBACK, acquire);
    CuAssertPtrEquals(ct, NULL, handle);
    CuAssertIntEquals(ct, 0, stats_status);
    CuAssertIntEquals(ct, 0, (int)stats.files_copied);
    CuAssertIntEquals(ct, 0, (int)stats.resident_bytes);
    if (stats.backend_supported) {
        CuAssertIntEquals(
            ct, 1, (int)stats.budget_refusals);
    }
    CuAssertIntEquals(ct, 0, stop_status);
}
