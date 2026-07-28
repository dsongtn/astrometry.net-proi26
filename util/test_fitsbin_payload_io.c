/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "fitsbin.h"

#include "cutest.h"

#define PAYLOAD_FIXTURE_BYTES 128U
#define PAYLOAD_WRAPPER_RECORDS 16U
#define PAYLOAD_CREDIT_THREADS 4

typedef enum payload_wrapper_mode {
    PAYLOAD_WRAPPER_PASS,
    PAYLOAD_WRAPPER_SHORT,
    PAYLOAD_WRAPPER_EOF,
    PAYLOAD_WRAPPER_BLOCK
} payload_wrapper_mode_t;

typedef struct payload_wrapper_state {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    payload_wrapper_mode_t mode;
    int calls;
    int active;
    int max_active;
    int release;
    int first_fd;
    int fd_mismatch;
    off_t offsets[PAYLOAD_WRAPPER_RECORDS];
    size_t requests[PAYLOAD_WRAPPER_RECORDS];
} payload_wrapper_state_t;

typedef struct payload_fixture {
    fitsbin_t* fitsbin;
    fitsbin_chunk_t* chunk;
    unsigned char bytes[PAYLOAD_FIXTURE_BYTES];
    char filename[128];
} payload_fixture_t;

typedef struct payload_start_gate {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    int ready;
    int start;
} payload_start_gate_t;

typedef struct payload_thread {
    payload_start_gate_t* gate;
    fitsbin_t* fitsbin;
    const unsigned char* source;
    unsigned char destination[8];
    int rc;
    int error;
} payload_thread_t;

typedef struct payload_credit_result {
    int calls_before_release;
    int active_before_release;
    int max_active;
    int total_calls;
    int first_fd;
    int fd_mismatch;
    int all_reads_ok;
} payload_credit_result_t;

static payload_wrapper_state_t payload_wrapper = {
    PTHREAD_MUTEX_INITIALIZER,
    PTHREAD_COND_INITIALIZER,
    PAYLOAD_WRAPPER_PASS,
    0,
    0,
    0,
    0,
    -1,
    0,
    {0},
    {0}
};

/*
 * util/fitsbin.c is compiled with _FILE_OFFSET_BITS=64 and therefore calls
 * pread64.  This test target uses --wrap=pread64 so production objects and
 * APIs remain unchanged.
 */
extern ssize_t __real_pread64(
    int fd,
    void* destination,
    size_t size,
    off_t offset);

ssize_t __wrap_pread64(
    int fd,
    void* destination,
    size_t size,
    off_t offset) {
    payload_wrapper_mode_t mode;
    int call;
    size_t request = size;
    ssize_t result;

    pthread_mutex_lock(&payload_wrapper.mutex);
    mode = payload_wrapper.mode;
    if (mode == PAYLOAD_WRAPPER_PASS) {
        pthread_mutex_unlock(&payload_wrapper.mutex);
        return __real_pread64(fd, destination, size, offset);
    }

    call = payload_wrapper.calls++;
    if ((size_t)call < PAYLOAD_WRAPPER_RECORDS) {
        payload_wrapper.offsets[call] = offset;
        payload_wrapper.requests[call] = size;
    }

    if (mode == PAYLOAD_WRAPPER_BLOCK) {
        payload_wrapper.active++;
        if (payload_wrapper.active > payload_wrapper.max_active) {
            payload_wrapper.max_active = payload_wrapper.active;
        }
        if (payload_wrapper.first_fd < 0) {
            payload_wrapper.first_fd = fd;
        } else if (payload_wrapper.first_fd != fd) {
            payload_wrapper.fd_mismatch = 1;
        }
        pthread_cond_broadcast(&payload_wrapper.condition);
        while (!payload_wrapper.release) {
            pthread_cond_wait(
                &payload_wrapper.condition,
                &payload_wrapper.mutex);
        }
        pthread_mutex_unlock(&payload_wrapper.mutex);

        result = __real_pread64(fd, destination, size, offset);

        pthread_mutex_lock(&payload_wrapper.mutex);
        payload_wrapper.active--;
        pthread_cond_broadcast(&payload_wrapper.condition);
        pthread_mutex_unlock(&payload_wrapper.mutex);
        return result;
    }
    pthread_mutex_unlock(&payload_wrapper.mutex);

    if (mode == PAYLOAD_WRAPPER_SHORT) {
        if (call == 0) {
            errno = EINTR;
            return -1;
        }
        if (call == 1 && request > 3U) {
            request = 3U;
        } else if (call == 2 && request > 2U) {
            request = 2U;
        }
        return __real_pread64(
            fd, destination, request, offset);
    }

    if (call == 0) {
        if (request > 3U) {
            request = 3U;
        }
        return __real_pread64(
            fd, destination, request, offset);
    }
    return 0;
}

static void payload_wrapper_reset(
    payload_wrapper_mode_t mode) {
    pthread_mutex_lock(&payload_wrapper.mutex);
    payload_wrapper.mode = mode;
    payload_wrapper.calls = 0;
    payload_wrapper.active = 0;
    payload_wrapper.max_active = 0;
    payload_wrapper.release = 0;
    payload_wrapper.first_fd = -1;
    payload_wrapper.fd_mismatch = 0;
    memset(payload_wrapper.offsets, 0,
           sizeof(payload_wrapper.offsets));
    memset(payload_wrapper.requests, 0,
           sizeof(payload_wrapper.requests));
    pthread_mutex_unlock(&payload_wrapper.mutex);
}

static void payload_wrapper_release(void) {
    pthread_mutex_lock(&payload_wrapper.mutex);
    payload_wrapper.release = 1;
    pthread_cond_broadcast(&payload_wrapper.condition);
    pthread_mutex_unlock(&payload_wrapper.mutex);
}

static int payload_wrapper_wait_for_calls(
    int expected,
    int timeout_seconds) {
    struct timespec deadline;
    int status = 0;

    if (clock_gettime(CLOCK_REALTIME, &deadline)) {
        return -1;
    }
    deadline.tv_sec += timeout_seconds;

    pthread_mutex_lock(&payload_wrapper.mutex);
    while (payload_wrapper.calls < expected &&
           status != ETIMEDOUT) {
        status = pthread_cond_timedwait(
            &payload_wrapper.condition,
            &payload_wrapper.mutex,
            &deadline);
    }
    expected = payload_wrapper.calls >= expected;
    pthread_mutex_unlock(&payload_wrapper.mutex);
    return expected ? 0 : -1;
}

static int payload_fixture_open(payload_fixture_t* fixture) {
    fitsbin_t* output;
    fitsbin_chunk_t output_chunk;
    fitsbin_chunk_t input_chunk;
    int fd;
    size_t i;

    memset(fixture, 0, sizeof(*fixture));
    snprintf(
        fixture->filename,
        sizeof(fixture->filename),
        "/tmp/test-fitsbin-payload-io.XXXXXX");
    fd = mkstemp(fixture->filename);
    if (fd < 0) {
        return -1;
    }
    close(fd);

    for (i = 0U; i < sizeof(fixture->bytes); i++) {
        fixture->bytes[i] =
            (unsigned char)((i * 37U + 11U) & 0xffU);
    }

    output = fitsbin_open_for_writing(
        fixture->filename);
    if (!output) {
        goto fail;
    }
    fitsbin_chunk_init(&output_chunk);
    output_chunk.tablename = "payload-io";
    output_chunk.itemsize = 1;
    output_chunk.nrows = (int)sizeof(fixture->bytes);
    output_chunk.data = fixture->bytes;
    if (fitsbin_write_primary_header(output) ||
        fitsbin_write_chunk(output, &output_chunk) ||
        fitsbin_fix_primary_header(output)) {
        fitsbin_close(output);
        output = NULL;
        fitsbin_chunk_clean(&output_chunk);
        goto fail;
    }
    if (fitsbin_close(output)) {
        output = NULL;
        fitsbin_chunk_clean(&output_chunk);
        goto fail;
    }
    output = NULL;
    fitsbin_chunk_clean(&output_chunk);

    fixture->fitsbin = fitsbin_open(fixture->filename);
    if (!fixture->fitsbin) {
        goto fail;
    }
    fitsbin_chunk_init(&input_chunk);
    input_chunk.tablename = "payload-io";
    if (fitsbin_read_chunk(
            fixture->fitsbin,
            &input_chunk)) {
        goto fail;
    }
    fixture->chunk = fitsbin_get_chunk(
        fixture->fitsbin, 0);
    if (!fixture->chunk ||
        !fixture->chunk->data ||
        fixture->chunk->data_file_size !=
            sizeof(fixture->bytes)) {
        goto fail;
    }
    return 0;

fail:
    if (output) {
        fitsbin_close(output);
    }
    if (fixture->fitsbin) {
        fitsbin_close(fixture->fitsbin);
        fixture->fitsbin = NULL;
    }
    unlink(fixture->filename);
    return -1;
}

static void payload_fixture_close(
    payload_fixture_t* fixture) {
    if (fixture->fitsbin) {
        fitsbin_close(fixture->fitsbin);
    }
    unlink(fixture->filename);
    memset(fixture, 0, sizeof(*fixture));
}

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

static void* payload_credit_thread(void* opaque) {
    payload_thread_t* thread = opaque;

    pthread_mutex_lock(&thread->gate->mutex);
    thread->gate->ready++;
    pthread_cond_broadcast(&thread->gate->condition);
    while (!thread->gate->start) {
        pthread_cond_wait(
            &thread->gate->condition,
            &thread->gate->mutex);
    }
    pthread_mutex_unlock(&thread->gate->mutex);

    errno = 0;
    thread->rc = fitsbin_pread_mapped_range(
        thread->fitsbin,
        thread->source,
        sizeof(thread->destination),
        thread->destination);
    thread->error = errno;
    return NULL;
}

static int payload_credit_round(
    payload_fixture_t* fixture,
    int worker_count,
    int thread_count,
    int expected_limit,
    payload_credit_result_t* result) {
    payload_start_gate_t gate = {
        PTHREAD_MUTEX_INITIALIZER,
        PTHREAD_COND_INITIALIZER,
        0,
        0
    };
    payload_thread_t threads[PAYLOAD_CREDIT_THREADS];
    pthread_t ids[PAYLOAD_CREDIT_THREADS];
    int created = 0;
    int i;
    int wait_ok;
    int rc = -1;

    memset(result, 0, sizeof(*result));
    memset(threads, 0, sizeof(threads));
    fitsbin_payload_io_configure_workers(worker_count);
    payload_wrapper_reset(PAYLOAD_WRAPPER_BLOCK);

    for (i = 0; i < thread_count; i++) {
        threads[i].gate = &gate;
        threads[i].fitsbin = fixture->fitsbin;
        threads[i].source =
            (const unsigned char*)fixture->chunk->data +
            (size_t)i * sizeof(threads[i].destination);
        if (pthread_create(
                &ids[i],
                NULL,
                payload_credit_thread,
                &threads[i])) {
            break;
        }
        created++;
    }

    pthread_mutex_lock(&gate.mutex);
    while (gate.ready < created) {
        pthread_cond_wait(&gate.condition, &gate.mutex);
    }
    gate.start = 1;
    pthread_cond_broadcast(&gate.condition);
    pthread_mutex_unlock(&gate.mutex);

    wait_ok =
        created == thread_count &&
        payload_wrapper_wait_for_calls(
            expected_limit, 3) == 0;

    pthread_mutex_lock(&payload_wrapper.mutex);
    result->calls_before_release =
        payload_wrapper.calls;
    result->active_before_release =
        payload_wrapper.active;
    pthread_mutex_unlock(&payload_wrapper.mutex);

    payload_wrapper_release();
    for (i = 0; i < created; i++) {
        pthread_join(ids[i], NULL);
    }

    pthread_mutex_lock(&payload_wrapper.mutex);
    result->max_active =
        payload_wrapper.max_active;
    result->total_calls =
        payload_wrapper.calls;
    result->first_fd =
        payload_wrapper.first_fd;
    result->fd_mismatch =
        payload_wrapper.fd_mismatch;
    pthread_mutex_unlock(&payload_wrapper.mutex);
    payload_wrapper_reset(PAYLOAD_WRAPPER_PASS);

    result->all_reads_ok =
        created == thread_count;
    for (i = 0; i < created; i++) {
        if (threads[i].rc ||
            memcmp(
                threads[i].destination,
                fixture->bytes +
                    (size_t)i *
                        sizeof(threads[i].destination),
                sizeof(threads[i].destination))) {
            result->all_reads_ok = 0;
        }
    }

    if (wait_ok && created == thread_count) {
        rc = 0;
    }
    pthread_cond_destroy(&gate.condition);
    pthread_mutex_destroy(&gate.mutex);
    return rc;
}

void test_fitsbin_payload_shared_reader_credit(CuTest* ct) {
    payload_fixture_t fixture;
    payload_credit_result_t two_readers;
    payload_credit_result_t one_reader;

    CuAssertIntEquals(ct, 0, payload_fixture_open(&fixture));

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
