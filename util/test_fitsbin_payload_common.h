/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#ifndef TEST_FITSBIN_PAYLOAD_COMMON_H
#define TEST_FITSBIN_PAYLOAD_COMMON_H

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

#define PAYLOAD_FIXTURE_BYTES (128U * 1024U)
#define PAYLOAD_WRAPPER_RECORDS 16U
#define PAYLOAD_CREDIT_THREADS 4
#define PAYLOAD_COMPLETION_EXPIRY_TICKETS 32U

typedef enum payload_wrapper_mode {
    PAYLOAD_WRAPPER_PASS,
    PAYLOAD_WRAPPER_RECORD,
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

typedef struct payload_wait_helper_state {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    int calls;
} payload_wait_helper_state_t;

typedef struct payload_ticket_wait_state {
    fitsbin_t* fitsbin;
    fitsbin_payload_io_ticket_t* ticket;
    payload_wait_helper_state_t* helper;
    int configured;
    int result;
    int error;
} payload_ticket_wait_state_t;

typedef struct payload_ticket_drain_state {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    fitsbin_t* fitsbin;
    fitsbin_payload_io_ticket_t* ticket;
    int started;
    int status;
    int result;
    int error;
} payload_ticket_drain_state_t;

typedef struct payload_planned_state {
    fitsbin_prefetch_range_t range;
    int result;
    int calls;
} payload_planned_state_t;

typedef struct payload_readahead_state {
    pthread_mutex_t mutex;
    int recording;
    int calls;
    off_t offsets[PAYLOAD_WRAPPER_RECORDS];
    size_t requests[PAYLOAD_WRAPPER_RECORDS];
} payload_readahead_state_t;

extern payload_wrapper_state_t payload_wrapper;
extern payload_readahead_state_t payload_readahead;

void payload_readahead_reset(int recording);
void payload_wrapper_reset(payload_wrapper_mode_t mode);
void payload_wrapper_release(void);
int payload_wrapper_wait_for_calls(
    int expected,
    int timeout_seconds);
int payload_planned_ranges(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count);
int payload_wait_helper_wait_for_calls(
    payload_wait_helper_state_t* helper,
    int expected,
    int timeout_seconds);
void* payload_ticket_wait_thread(void* opaque);
void* payload_ticket_drain_thread(void* opaque);
int payload_ticket_drain_wait_started(
    payload_ticket_drain_state_t* state,
    int timeout_seconds);
int payload_fixture_open(payload_fixture_t* fixture);
void payload_fixture_close(payload_fixture_t* fixture);
int payload_credit_round(
    payload_fixture_t* fixture,
    int worker_count,
    int thread_count,
    int expected_limit,
    payload_credit_result_t* result);

#endif
