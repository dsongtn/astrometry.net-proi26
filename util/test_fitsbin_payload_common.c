/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include "test_fitsbin_payload_common.h"

payload_wrapper_state_t payload_wrapper = {
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

payload_readahead_state_t payload_readahead = {
    PTHREAD_MUTEX_INITIALIZER,
    0,
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

    if (mode == PAYLOAD_WRAPPER_RECORD) {
        return __real_pread64(fd, destination, size, offset);
    }
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

extern ssize_t __real_readahead(
    int fd,
    off_t offset,
    size_t size);

ssize_t __wrap_readahead(
    int fd,
    off_t offset,
    size_t size) {
    int call = -1;

    pthread_mutex_lock(&payload_readahead.mutex);
    if (payload_readahead.recording) {
        call = payload_readahead.calls++;
        if ((size_t)call < PAYLOAD_WRAPPER_RECORDS) {
            payload_readahead.offsets[call] = offset;
            payload_readahead.requests[call] = size;
        }
    }
    pthread_mutex_unlock(&payload_readahead.mutex);
    return __real_readahead(fd, offset, size);
}

void payload_readahead_reset(int recording) {
    pthread_mutex_lock(&payload_readahead.mutex);
    payload_readahead.recording = recording;
    payload_readahead.calls = 0;
    memset(payload_readahead.offsets, 0,
           sizeof(payload_readahead.offsets));
    memset(payload_readahead.requests, 0,
           sizeof(payload_readahead.requests));
    pthread_mutex_unlock(&payload_readahead.mutex);
}

void payload_wrapper_reset(
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

void payload_wrapper_release(void) {
    pthread_mutex_lock(&payload_wrapper.mutex);
    payload_wrapper.release = 1;
    pthread_cond_broadcast(&payload_wrapper.condition);
    pthread_mutex_unlock(&payload_wrapper.mutex);
}

int payload_wrapper_wait_for_calls(
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

static int payload_ticket_wait_helper(void* opaque) {
    payload_wait_helper_state_t* helper = opaque;

    pthread_mutex_lock(&helper->mutex);
    helper->calls++;
    pthread_cond_broadcast(&helper->condition);
    pthread_mutex_unlock(&helper->mutex);
    return 0;
}

int payload_planned_ranges(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count) {
    payload_planned_state_t* state = opaque;

    if (!state || !cancelled || !ranges || !range_count ||
        !range_capacity) {
        errno = EINVAL;
        return -1;
    }
    state->calls++;
    *range_count = 0U;
    if (cancelled(cancel_opaque)) {
        errno = ECANCELED;
        return -1;
    }
    if (state->result < 0) {
        errno = EIO;
        return -1;
    }
    if (!state->result) {
        return 0;
    }
    ranges[0] = state->range;
    *range_count = 1U;
    return 1;
}

int payload_wait_helper_wait_for_calls(
    payload_wait_helper_state_t* helper,
    int expected,
    int timeout_seconds) {
    struct timespec deadline;
    int status = 0;

    if (clock_gettime(CLOCK_REALTIME, &deadline)) {
        return -1;
    }
    deadline.tv_sec += timeout_seconds;

    pthread_mutex_lock(&helper->mutex);
    while (helper->calls < expected &&
           status != ETIMEDOUT) {
        status = pthread_cond_timedwait(
            &helper->condition,
            &helper->mutex,
            &deadline);
    }
    expected = helper->calls >= expected;
    pthread_mutex_unlock(&helper->mutex);
    return expected ? 0 : -1;
}

void* payload_ticket_wait_thread(void* opaque) {
    payload_ticket_wait_state_t* state = opaque;

    if (fitsbin_payload_io_set_thread_wait_helper(
            payload_ticket_wait_helper,
            NULL,
            state->helper)) {
        state->configured = -1;
        return NULL;
    }
    state->configured = 1;
    errno = 0;
    state->result = fitsbin_payload_io_ticket_wait(
        state->fitsbin, state->ticket);
    state->error = errno;
    fitsbin_payload_io_clear_thread_wait_helper();
    return NULL;
}

void* payload_ticket_drain_thread(void* opaque) {
    payload_ticket_drain_state_t* state = opaque;

    pthread_mutex_lock(&state->mutex);
    state->started = 1;
    pthread_cond_broadcast(&state->condition);
    pthread_mutex_unlock(&state->mutex);
    errno = 0;
    state->status = fitsbin_payload_io_ticket_drain_and_destroy(
        state->fitsbin,
        &state->ticket,
        &state->result);
    state->error = errno;
    return NULL;
}

int payload_ticket_drain_wait_started(
    payload_ticket_drain_state_t* state,
    int timeout_seconds) {
    struct timespec deadline;
    int status = 0;

    if (clock_gettime(CLOCK_REALTIME, &deadline)) {
        return -1;
    }
    deadline.tv_sec += timeout_seconds;
    pthread_mutex_lock(&state->mutex);
    while (!state->started && status != ETIMEDOUT) {
        status = pthread_cond_timedwait(
            &state->condition,
            &state->mutex,
            &deadline);
    }
    status = state->started ? 0 : -1;
    pthread_mutex_unlock(&state->mutex);
    return status;
}

int payload_fixture_open(payload_fixture_t* fixture) {
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

void payload_fixture_close(
    payload_fixture_t* fixture) {
    if (fixture->fitsbin) {
        fitsbin_close(fixture->fitsbin);
    }
    unlink(fixture->filename);
    memset(fixture, 0, sizeof(*fixture));
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

int payload_credit_round(
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
