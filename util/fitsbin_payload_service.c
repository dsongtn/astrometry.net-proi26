/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include <errno.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

#include "keywords.h"
#include "fitsbin.h"
#include "fitsbin_internal.h"
#include "fitsioutils.h"
#include "ioutils.h"
#include "fitsfile.h"
#include "errors.h"
#include "an-endian.h"
#include "tic.h"
#include "log.h"

typedef enum fitsbin_payload_io_ticket_kind {
    FITSBIN_PAYLOAD_IO_TICKET_PREFETCH = 0,
    FITSBIN_PAYLOAD_IO_TICKET_DIRECT
} fitsbin_payload_io_ticket_kind_t;

typedef enum fitsbin_payload_io_ticket_state {
    FITSBIN_PAYLOAD_IO_PLANNED = 0,
    FITSBIN_PAYLOAD_IO_SUBMITTED,
    FITSBIN_PAYLOAD_IO_READY,
    FITSBIN_PAYLOAD_IO_FAILED,
    FITSBIN_PAYLOAD_IO_CANCELLED
} fitsbin_payload_io_ticket_state_t;

struct fitsbin_payload_io_ticket {
    fitsbin_prefetch_range_t
        planned_ranges[FITSBIN_PREFETCH_RANGE_LIMIT];
    fitsbin_mapped_span_t spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    fitsbin_file_span_t
        queued_spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    fitsbin_prepared_pread_range_t* ranges;
    fitsbin_t* source;
    size_t span_count;
    size_t exact_span_count;
    unsigned long long reused_page_count;
    size_t coalesced_gap_count;
    size_t coalesced_gap_bytes;
    size_t queued_span_count;
    size_t queued_exact_span_count;
    size_t queued_gap_count;
    size_t queued_gap_bytes;
    size_t range_count;
    size_t byte_count;
    size_t admission_byte_count;
    size_t plan_byte_budget;
    size_t logical_byte_count;
    size_t queued_byte_count;
    size_t queued_ranges_submitted;
    size_t queued_bytes_submitted;
    unsigned long long page_count;
    unsigned long long sequence;
    unsigned long long read_nanoseconds;
    fitsbin_payload_io_ticket_state_t state;
    fitsbin_payload_io_ticket_kind_t kind;
    fitsbin_payload_io_priority_t priority;
    fitsbin_payload_io_plan_fn plan;
    void* plan_opaque;
    struct fitsbin_payload_io_ticket* next;
    struct fitsbin_payload_io_ticket* wait_next;
    pthread_cond_t completion_cv;
    int fd;
    int saved_errno;
    anbool cancel_requested;
    anbool counters_applied;
    anbool queued_prepare_failed;
    anbool queued_submit_failed;
    anbool wait_registered;
    anbool waiter_active;
    anbool owner_draining;
    anbool helper_waiter;
};

static pthread_mutex_t fitsbin_payload_io_mutex =
    PTHREAD_MUTEX_INITIALIZER;
static ASTROMETRY_THREAD_LOCAL anbool
    fitsbin_payload_io_planning_active = FALSE;
/* Loader queue publication and service lifecycle changes. */
static pthread_cond_t fitsbin_payload_io_cv =
    PTHREAD_COND_INITIALIZER;
/* Reader-credit availability, kept separate from ticket waiters. */
static pthread_cond_t fitsbin_payload_io_credit_cv =
    PTHREAD_COND_INITIALIZER;
/* Completion-notifier unregister waits only for callbacks already in flight. */
static pthread_cond_t fitsbin_payload_io_completion_cv =
    PTHREAD_COND_INITIALIZER;
static int fitsbin_payload_io_capacity = 1;
static int fitsbin_payload_io_limit = 1;
static int fitsbin_payload_io_active = 0;
static size_t fitsbin_payload_io_waiters = 0U;
static size_t fitsbin_payload_io_wait_helpers = 0U;
static size_t fitsbin_payload_io_wait_helpers_active = 0U;
static unsigned int fitsbin_payload_io_helper_windows = 0U;
static unsigned long long fitsbin_payload_io_work_epoch = 0ULL;
static pthread_t
    fitsbin_payload_io_threads[FITSBIN_PAYLOAD_IO_MAX_LANES];
static fitsbin_payload_io_ticket_t*
    fitsbin_payload_io_queue_head[FITSBIN_PAYLOAD_IO_PRIORITY_COUNT];
static fitsbin_payload_io_ticket_t*
    fitsbin_payload_io_queue_tail[FITSBIN_PAYLOAD_IO_PRIORITY_COUNT];
/* Registered waiters are linked only while holding payload_io_mutex. */
static fitsbin_payload_io_ticket_t*
    fitsbin_payload_io_wait_head = NULL;
static size_t fitsbin_payload_io_service_jobs = 0U;
static size_t fitsbin_payload_io_service_bytes = 0U;
static unsigned long long fitsbin_payload_io_next_sequence = 0ULL;
static unsigned long long fitsbin_payload_io_service_submitted = 0ULL;
static unsigned long long fitsbin_payload_io_service_ready = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_direct_submitted = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_direct_ready = 0ULL;

typedef struct fitsbin_payload_transport_metrics {
    unsigned long long direct_attempts;
    unsigned long long direct_policy_refused;
    unsigned long long mapped_submitted;
    unsigned long long mapped_ready;
    unsigned long long mapped_immediate_ready;
    unsigned long long direct_ranges;
    unsigned long long direct_bytes;
    unsigned long long mapped_spans;
    unsigned long long mapped_exact_spans;
    unsigned long long mapped_reused_pages;
    unsigned long long mapped_gap_merges;
    unsigned long long mapped_gap_bytes;
    unsigned long long mapped_bytes;
    unsigned long long queued_exact_spans;
    unsigned long long queued_gap_merges;
    unsigned long long queued_gap_bytes;
    unsigned long long pread_calls;
    unsigned long long pread_bytes;
    unsigned long long readahead_calls;
    unsigned long long readahead_bytes;
    unsigned long long populate_calls;
    unsigned long long populate_bytes;
} fitsbin_payload_transport_metrics_t;

static fitsbin_payload_transport_metrics_t
    fitsbin_payload_io_transport_metrics;
static unsigned long long fitsbin_payload_io_service_cancelled = 0ULL;
static unsigned long long fitsbin_payload_io_service_failed = 0ULL;
static unsigned long long fitsbin_payload_io_service_planned_bytes = 0ULL;
static unsigned long long fitsbin_payload_io_service_plan_nanoseconds = 0ULL;
static unsigned long long fitsbin_payload_io_service_read_nanoseconds = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_queued_ranges = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_queued_bytes = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_queue_prepare_failures = 0ULL;
static unsigned long long
    fitsbin_payload_io_service_queue_submit_failures = 0ULL;
static int fitsbin_payload_io_service_lanes = 0;
static anbool fitsbin_payload_io_service_running = FALSE;
static anbool fitsbin_payload_io_transport_metrics_enabled = FALSE;
static anbool fitsbin_payload_io_service_accepting = FALSE;
static anbool fitsbin_payload_io_service_stopping = FALSE;
static fitsbin_payload_io_completion_notify_fn
    fitsbin_payload_io_completion_notify = NULL;
static void* fitsbin_payload_io_completion_opaque = NULL;
static size_t fitsbin_payload_io_completion_active = 0U;
static anbool fitsbin_payload_io_completion_clearing = FALSE;
static ASTROMETRY_THREAD_LOCAL unsigned int
    fitsbin_payload_io_completion_dispatch_depth = 0U;
static ASTROMETRY_THREAD_LOCAL
    fitsbin_payload_io_wait_helper_fn
    fitsbin_payload_io_thread_wait_helper = NULL;
static ASTROMETRY_THREAD_LOCAL
    fitsbin_payload_io_stop_check_fn
    fitsbin_payload_io_thread_stop_check = NULL;
static ASTROMETRY_THREAD_LOCAL void*
    fitsbin_payload_io_thread_wait_opaque = NULL;
static ASTROMETRY_THREAD_LOCAL anbool
    fitsbin_payload_io_thread_wait_active = FALSE;
static ASTROMETRY_THREAD_LOCAL unsigned long long
    fitsbin_payload_io_thread_work_epoch = 0ULL;

unsigned long long fitsbin_payload_io_sequence_hint(void) {
    unsigned long long sequence = __atomic_load_n(
        &fitsbin_payload_io_next_sequence, __ATOMIC_ACQUIRE);

    return sequence == ULLONG_MAX ? 0ULL : sequence + 1ULL;
}

unsigned long long fitsbin_timespec_delta_nanoseconds(
    const struct timespec* finish,
    const struct timespec* start) {
    time_t seconds;
    long nanoseconds;

    if (!finish || !start) {
        return 0U;
    }
    seconds = finish->tv_sec - start->tv_sec;
    nanoseconds = finish->tv_nsec - start->tv_nsec;
    if (nanoseconds < 0) {
        seconds--;
        nanoseconds += 1000000000L;
    }
    if (seconds < 0 ||
        (unsigned long long)seconds >
            (ULLONG_MAX - (unsigned long long)nanoseconds) /
                1000000000ULL) {
        return 0U;
    }
    return (unsigned long long)seconds * 1000000000ULL +
        (unsigned long long)nanoseconds;
}

/* fitsbin_payload_io_mutex must be held. */
static fitsbin_payload_io_completion_notify_fn
fitsbin_payload_io_completion_acquire_locked(void** opaque) {
    fitsbin_payload_io_completion_notify_fn notify =
        fitsbin_payload_io_completion_notify;

    if (opaque) {
        *opaque = NULL;
    }
    if (!notify) {
        return NULL;
    }
    fitsbin_payload_io_completion_active++;
    if (opaque) {
        *opaque = fitsbin_payload_io_completion_opaque;
    }
    return notify;
}

static anbool fitsbin_payload_io_ticket_cancelled(void* opaque) {
    const fitsbin_payload_io_ticket_t* ticket = opaque;

    return ticket && __atomic_load_n(
        &ticket->cancel_requested,
        __ATOMIC_ACQUIRE);
}

static void fitsbin_payload_io_completion_release(void) {
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    assert(fitsbin_payload_io_completion_active > 0U);
    fitsbin_payload_io_completion_active--;
    pthread_cond_broadcast(&fitsbin_payload_io_completion_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

static void fitsbin_payload_io_completion_dispatch(
    fitsbin_payload_io_completion_notify_fn notify,
    void* opaque,
    unsigned long long completion_id) {
    if (!notify) {
        return;
    }
    fitsbin_payload_io_completion_dispatch_depth++;
    notify(opaque, completion_id);
    fitsbin_payload_io_completion_dispatch_depth--;
    fitsbin_payload_io_completion_release();
}

int fitsbin_payload_io_set_completion_notifier(
    fitsbin_payload_io_completion_notify_fn notify,
    void* opaque) {
    if (!notify) {
        errno = EINVAL;
        return -1;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_completion_notify ||
        fitsbin_payload_io_completion_clearing) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    fitsbin_payload_io_completion_notify = notify;
    fitsbin_payload_io_completion_opaque = opaque;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    return 0;
}

int fitsbin_payload_io_clear_completion_notifier(
    fitsbin_payload_io_completion_notify_fn notify,
    void* opaque) {
    if (!notify) {
        errno = EINVAL;
        return -1;
    }
    if (fitsbin_payload_io_completion_dispatch_depth) {
        errno = EDEADLK;
        return -1;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_completion_notify != notify ||
        fitsbin_payload_io_completion_opaque != opaque) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = ENOENT;
        return -1;
    }
    if (fitsbin_payload_io_completion_clearing) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    fitsbin_payload_io_completion_clearing = TRUE;
    fitsbin_payload_io_completion_notify = NULL;
    fitsbin_payload_io_completion_opaque = NULL;
    while (fitsbin_payload_io_completion_active) {
        pthread_cond_wait(
            &fitsbin_payload_io_completion_cv,
            &fitsbin_payload_io_mutex);
    }
    fitsbin_payload_io_completion_clearing = FALSE;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    return 0;
}

/* fitsbin_payload_io_mutex must be held. */
static void fitsbin_payload_io_update_limit_locked(void) {
    int capacity = __atomic_load_n(
        &fitsbin_payload_io_capacity,
        __ATOMIC_ACQUIRE);
    int limit = capacity;

    if (fitsbin_payload_io_helper_windows && limit > 1) {
        limit--;
    }
    __atomic_store_n(
        &fitsbin_payload_io_limit,
        limit,
        __ATOMIC_RELEASE);
}

void fitsbin_payload_io_configure_workers(int worker_count) {
    int capacity = MAX(worker_count, 1);

#if defined(FITSBIN_TEST_PAYLOAD_IO_LIMIT)
    if (FITSBIN_TEST_PAYLOAD_IO_LIMIT > 0) {
        capacity = MIN(capacity, FITSBIN_TEST_PAYLOAD_IO_LIMIT);
        capacity = MAX(capacity, 1);
    }
#endif
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    __atomic_store_n(
        &fitsbin_payload_io_capacity,
        capacity,
        __ATOMIC_RELEASE);
    fitsbin_payload_io_update_limit_locked();
    pthread_cond_broadcast(&fitsbin_payload_io_credit_cv);
    pthread_cond_broadcast(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

int fitsbin_payload_io_demand_busy(void) {
    int active = __atomic_load_n(
        &fitsbin_payload_io_active,
        __ATOMIC_ACQUIRE);
    int limit = __atomic_load_n(
        &fitsbin_payload_io_limit,
        __ATOMIC_ACQUIRE);

    return __atomic_load_n(
        &fitsbin_payload_io_waiters,
        __ATOMIC_ACQUIRE) > 0U ||
        active >= limit;
}

int fitsbin_payload_io_set_thread_wait_helper(
    fitsbin_payload_io_wait_helper_fn helper,
    fitsbin_payload_io_stop_check_fn stop_check,
    void* opaque) {
    if (!helper) {
        errno = EINVAL;
        return -1;
    }
    if (fitsbin_payload_io_thread_wait_helper ||
        fitsbin_payload_io_thread_stop_check ||
        fitsbin_payload_io_thread_wait_active) {
        errno = EBUSY;
        return -1;
    }
    fitsbin_payload_io_thread_wait_helper = helper;
    fitsbin_payload_io_thread_stop_check = stop_check;
    fitsbin_payload_io_thread_wait_opaque = opaque;
    return 0;
}

void fitsbin_payload_io_clear_thread_wait_helper(void) {
    if (fitsbin_payload_io_thread_wait_active) {
        return;
    }
    if (!fitsbin_payload_io_thread_wait_helper &&
        !fitsbin_payload_io_thread_stop_check) {
        return;
    }
    fitsbin_payload_io_thread_wait_helper = NULL;
    fitsbin_payload_io_thread_stop_check = NULL;
    fitsbin_payload_io_thread_wait_opaque = NULL;
}

size_t fitsbin_payload_io_wait_helper_count(void) {
    size_t waiters = __atomic_load_n(
        &fitsbin_payload_io_wait_helpers,
        __ATOMIC_ACQUIRE);
    size_t active = __atomic_load_n(
        &fitsbin_payload_io_wait_helpers_active,
        __ATOMIC_ACQUIRE);

    return waiters > active ? waiters - active : 0U;
}

void fitsbin_payload_io_notify_wait_helpers(void) {
    fitsbin_payload_io_ticket_t* ticket;

    __atomic_add_fetch(
        &fitsbin_payload_io_work_epoch,
        1ULL,
        __ATOMIC_RELEASE);
    if (!fitsbin_payload_io_wait_helper_count()) {
        return;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_wait_helper_count()) {
        pthread_cond_broadcast(&fitsbin_payload_io_credit_cv);
        for (ticket = fitsbin_payload_io_wait_head;
             ticket;
             ticket = ticket->wait_next) {
            if (ticket->helper_waiter) {
                pthread_cond_signal(&ticket->completion_cv);
            }
        }
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

void fitsbin_payload_io_begin_helper_window(void) {
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_helper_windows < UINT_MAX) {
        fitsbin_payload_io_helper_windows++;
    }
    fitsbin_payload_io_update_limit_locked();
    pthread_cond_broadcast(&fitsbin_payload_io_credit_cv);
    pthread_cond_broadcast(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

void fitsbin_payload_io_end_helper_window(void) {
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_helper_windows) {
        fitsbin_payload_io_helper_windows--;
    }
    fitsbin_payload_io_update_limit_locked();
    pthread_cond_broadcast(&fitsbin_payload_io_credit_cv);
    pthread_cond_broadcast(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

static int fitsbin_payload_io_try_acquire(void) {
    int active = __atomic_load_n(
        &fitsbin_payload_io_active,
        __ATOMIC_RELAXED);

    while (1) {
        int limit = __atomic_load_n(
            &fitsbin_payload_io_limit,
            __ATOMIC_ACQUIRE);

        if (active >= limit) {
            return FALSE;
        }
        if (__atomic_compare_exchange_n(
                &fitsbin_payload_io_active,
                &active,
                active + 1,
                TRUE,
                __ATOMIC_ACQ_REL,
                __ATOMIC_RELAXED)) {
            return TRUE;
        }
    }
}

unsigned long long fitsbin_payload_io_acquire(void) {
    struct timespec wait_start;
    struct timespec wait_finish;
    anbool helper_waiter_registered = FALSE;
    anbool retry_helper = TRUE;
    anbool measured =
        clock_gettime(CLOCK_MONOTONIC, &wait_start) == 0;

    if (fitsbin_payload_io_try_acquire()) {
        return 0U;
    }

    __atomic_add_fetch(
        &fitsbin_payload_io_waiters,
        1U,
        __ATOMIC_RELEASE);
    if (fitsbin_payload_io_thread_wait_helper &&
        !fitsbin_payload_io_thread_wait_active) {
        __atomic_add_fetch(
            &fitsbin_payload_io_wait_helpers,
            1U,
            __ATOMIC_RELEASE);
        helper_waiter_registered = TRUE;
    }
    while (!fitsbin_payload_io_try_acquire()) {
        if (fitsbin_payload_io_thread_wait_helper &&
            !fitsbin_payload_io_thread_wait_active &&
            (retry_helper ||
             fitsbin_payload_io_thread_work_epoch !=
                 __atomic_load_n(
                     &fitsbin_payload_io_work_epoch,
                     __ATOMIC_ACQUIRE))) {
            fitsbin_payload_io_wait_helper_fn helper =
                fitsbin_payload_io_thread_wait_helper;
            void* opaque =
                fitsbin_payload_io_thread_wait_opaque;
            int helped;

            fitsbin_payload_io_thread_work_epoch =
                __atomic_load_n(
                    &fitsbin_payload_io_work_epoch,
                    __ATOMIC_ACQUIRE);
            fitsbin_payload_io_thread_wait_active = TRUE;
            __atomic_add_fetch(
                &fitsbin_payload_io_wait_helpers_active,
                1U,
                __ATOMIC_ACQ_REL);
            helped = helper(opaque);
            __atomic_sub_fetch(
                &fitsbin_payload_io_wait_helpers_active,
                1U,
                __ATOMIC_ACQ_REL);
            fitsbin_payload_io_thread_wait_active = FALSE;
            retry_helper = helped != 0;
            if (retry_helper) {
                continue;
            }
        }

        pthread_mutex_lock(&fitsbin_payload_io_mutex);
        if (fitsbin_payload_io_try_acquire()) {
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            break;
        }
        if (fitsbin_payload_io_thread_wait_helper &&
            fitsbin_payload_io_thread_work_epoch !=
                __atomic_load_n(
                    &fitsbin_payload_io_work_epoch,
                    __ATOMIC_ACQUIRE)) {
            retry_helper = TRUE;
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            continue;
        }
        pthread_cond_wait(
            &fitsbin_payload_io_credit_cv,
            &fitsbin_payload_io_mutex);
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    }
    if (helper_waiter_registered) {
        __atomic_sub_fetch(
            &fitsbin_payload_io_wait_helpers,
            1U,
            __ATOMIC_RELEASE);
    }
    __atomic_sub_fetch(
        &fitsbin_payload_io_waiters,
        1U,
        __ATOMIC_RELEASE);
    if (!measured ||
        clock_gettime(CLOCK_MONOTONIC, &wait_finish)) {
        return 0U;
    }
    return fitsbin_timespec_delta_nanoseconds(
        &wait_finish,
        &wait_start);
}

void fitsbin_payload_io_release(void) {
    int previous = __atomic_fetch_sub(
        &fitsbin_payload_io_active,
        1,
        __ATOMIC_ACQ_REL);

    assert(previous > 0);
    (void)previous;
    if (__atomic_load_n(
            &fitsbin_payload_io_waiters,
            __ATOMIC_ACQUIRE)) {
        pthread_mutex_lock(&fitsbin_payload_io_mutex);
        pthread_cond_signal(&fitsbin_payload_io_credit_cv);
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    }
}


anbool fitsbin_payload_io_planning_is_active(void) {
    return fitsbin_payload_io_planning_active;
}

int fitsbin_payload_io_capacity_current(void) {
    int capacity = __atomic_load_n(
        &fitsbin_payload_io_capacity,
        __ATOMIC_ACQUIRE);

    return capacity < 1 ? 1 : capacity;
}

static int fitsbin_payload_io_duplicate_fd(int fd) {
    int duplicate;

#if defined(F_DUPFD_CLOEXEC)
    duplicate = fcntl(fd, F_DUPFD_CLOEXEC, 0);
#else
    duplicate = dup(fd);
    if (duplicate >= 0) {
        (void)fcntl(duplicate, F_SETFD, FD_CLOEXEC);
    }
#endif
    return duplicate;
}

static anbool fitsbin_payload_io_priority_valid(
    fitsbin_payload_io_priority_t priority) {
    switch (priority) {
    case FITSBIN_PAYLOAD_IO_PRIORITY_DEMAND:
    case FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT:
    case FITSBIN_PAYLOAD_IO_PRIORITY_SPECULATIVE:
        return TRUE;
    }
    return FALSE;
}

/* fitsbin_payload_io_mutex must be held. */
static anbool fitsbin_payload_io_queue_empty_locked(void) {
    size_t priority;

    for (priority = 0U;
         priority < FITSBIN_PAYLOAD_IO_PRIORITY_COUNT;
         priority++) {
        if (fitsbin_payload_io_queue_head[priority]) {
            return FALSE;
        }
    }
    return TRUE;
}

/* fitsbin_payload_io_mutex must be held. */
static fitsbin_payload_io_ticket_t*
fitsbin_payload_io_dequeue_locked(void) {
    size_t priority;

    for (priority = 0U;
         priority < FITSBIN_PAYLOAD_IO_PRIORITY_COUNT;
         priority++) {
        fitsbin_payload_io_ticket_t* ticket =
            fitsbin_payload_io_queue_head[priority];

        if (!ticket) {
            continue;
        }
        fitsbin_payload_io_queue_head[priority] = ticket->next;
        if (!fitsbin_payload_io_queue_head[priority]) {
            fitsbin_payload_io_queue_tail[priority] = NULL;
        }
        ticket->next = NULL;
        return ticket;
    }
    return NULL;
}

/* fitsbin_payload_io_mutex must be held. */
static anbool fitsbin_payload_io_remove_queued_locked(
    fitsbin_payload_io_ticket_t* target) {
    size_t priority;

    if (!target) {
        return FALSE;
    }
    for (priority = 0U;
         priority < FITSBIN_PAYLOAD_IO_PRIORITY_COUNT;
         priority++) {
        fitsbin_payload_io_ticket_t* previous = NULL;
        fitsbin_payload_io_ticket_t* ticket =
            fitsbin_payload_io_queue_head[priority];

        while (ticket) {
            if (ticket == target) {
                if (previous) {
                    previous->next = ticket->next;
                } else {
                    fitsbin_payload_io_queue_head[priority] =
                        ticket->next;
                }
                if (fitsbin_payload_io_queue_tail[priority] ==
                    ticket) {
                    fitsbin_payload_io_queue_tail[priority] =
                        previous;
                }
                ticket->next = NULL;
                return TRUE;
            }
            previous = ticket;
            ticket = ticket->next;
        }
    }
    return FALSE;
}

/* fitsbin_payload_io_mutex must be held. */
static fitsbin_payload_io_completion_notify_fn
fitsbin_payload_io_cancel_queued_locked(
    fitsbin_payload_io_ticket_t* ticket,
    void** completion_opaque) {
    if (completion_opaque) {
        *completion_opaque = NULL;
    }
    if (!ticket ||
        ticket->state != FITSBIN_PAYLOAD_IO_SUBMITTED ||
        !__atomic_load_n(
            &ticket->cancel_requested,
            __ATOMIC_ACQUIRE) ||
        !fitsbin_payload_io_remove_queued_locked(ticket)) {
        return NULL;
    }
    assert(fitsbin_payload_io_service_jobs > 0U);
    assert(fitsbin_payload_io_service_bytes >=
           ticket->admission_byte_count);
    fitsbin_payload_io_service_jobs--;
    fitsbin_payload_io_service_bytes -=
        ticket->admission_byte_count;
    if (ticket->fd >= 0) {
        close(ticket->fd);
        ticket->fd = -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_CANCELLED;
    fitsbin_payload_io_service_cancelled++;
    pthread_cond_broadcast(&ticket->completion_cv);
    return fitsbin_payload_io_completion_acquire_locked(
        completion_opaque);
}

/* fitsbin_payload_io_mutex must be held. */
static int fitsbin_payload_io_admit_locked(
    const fitsbin_payload_io_ticket_t* ticket) {
    size_t max_jobs = FITSBIN_PAYLOAD_IO_MAX_JOBS;
    size_t max_bytes = FITSBIN_PAYLOAD_IO_MAX_BYTES;

    if (!ticket ||
        !fitsbin_payload_io_priority_valid(ticket->priority)) {
        return -1;
    }
    if (ticket->priority != FITSBIN_PAYLOAD_IO_PRIORITY_DEMAND) {
        max_jobs -= FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_JOBS;
        max_bytes -= FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_BYTES;
    }
    if (ticket->admission_byte_count > max_bytes) {
        return -1;
    }
    if (fitsbin_payload_io_service_jobs >= max_jobs ||
        fitsbin_payload_io_service_bytes >
            max_bytes - ticket->admission_byte_count) {
        return 0;
    }
    return 1;
}

/* fitsbin_payload_io_mutex must be held. */
static void fitsbin_payload_io_enqueue_locked(
    fitsbin_payload_io_ticket_t* ticket) {
    size_t priority = (size_t)ticket->priority;

    if (fitsbin_payload_io_queue_tail[priority]) {
        fitsbin_payload_io_queue_tail[priority]->next = ticket;
    } else {
        fitsbin_payload_io_queue_head[priority] = ticket;
    }
    fitsbin_payload_io_queue_tail[priority] = ticket;
}

static fitsbin_payload_io_ticket_t*
fitsbin_payload_io_ticket_alloc(void) {
    fitsbin_payload_io_ticket_t* ticket =
        calloc(1, sizeof(*ticket));
    int status;

    if (!ticket) {
        return NULL;
    }
    ticket->fd = -1;
    status = pthread_cond_init(&ticket->completion_cv, NULL);
    if (status) {
        free(ticket);
        errno = status;
        return NULL;
    }
    return ticket;
}

static void fitsbin_payload_io_ticket_free_storage(
    fitsbin_payload_io_ticket_t* ticket) {
    if (!ticket) {
        return;
    }
    assert(!ticket->wait_registered);
    assert(!ticket->waiter_active);
    if (ticket->fd >= 0) {
        close(ticket->fd);
    }
    free(ticket->ranges);
    pthread_cond_destroy(&ticket->completion_cv);
    free(ticket);
}

/* fitsbin_payload_io_mutex must be held. */
static int fitsbin_payload_io_register_waiter_locked(
    fitsbin_payload_io_ticket_t* ticket) {
    if (!ticket || ticket->wait_registered) {
        return -1;
    }
    ticket->wait_next = fitsbin_payload_io_wait_head;
    fitsbin_payload_io_wait_head = ticket;
    ticket->wait_registered = TRUE;
    return 0;
}

/* fitsbin_payload_io_mutex must be held. */
static int fitsbin_payload_io_unregister_waiter_locked(
    fitsbin_payload_io_ticket_t* ticket) {
    fitsbin_payload_io_ticket_t** link;

    if (!ticket || !ticket->wait_registered) {
        return -1;
    }
    link = &fitsbin_payload_io_wait_head;
    while (*link && *link != ticket) {
        link = &(*link)->wait_next;
    }
    if (!*link) {
        return -1;
    }
    *link = ticket->wait_next;
    ticket->wait_next = NULL;
    ticket->wait_registered = FALSE;
    ticket->helper_waiter = FALSE;
    return 0;
}

/*
 * Materialize a callback-owned plan into ticket-owned mapped and file spans.
 * The callback runs without payload_io_mutex and has exclusive access to its
 * opaque packet until this ticket reaches a terminal state.
 */
static int fitsbin_payload_io_prepare_planned_ticket(
    fitsbin_payload_io_ticket_t* ticket) {
    size_t range_count = 0U;
    int plan_status;

    if (!ticket || !ticket->plan || !ticket->source ||
        !ticket->plan_byte_budget) {
        errno = EINVAL;
        return -1;
    }
    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    fitsbin_payload_io_planning_active = TRUE;
    plan_status = ticket->plan(
        ticket->plan_opaque,
        fitsbin_payload_io_ticket_cancelled,
        ticket,
        ticket->planned_ranges,
        FITSBIN_PREFETCH_RANGE_LIMIT,
        &range_count);
    fitsbin_payload_io_planning_active = FALSE;
    if (fitsbin_payload_io_ticket_cancelled(ticket)) {
        errno = ECANCELED;
        return 1;
    }
    if (plan_status < 0) {
        if (!errno) {
            errno = EIO;
        }
        return -1;
    }
    if (!plan_status) {
        if (range_count) {
            errno = EINVAL;
            return -1;
        }
        return 0;
    }
    if (plan_status != 1 || !range_count ||
        range_count > FITSBIN_PREFETCH_RANGE_LIMIT) {
        errno = range_count > FITSBIN_PREFETCH_RANGE_LIMIT
            ? E2BIG
            : EINVAL;
        return -1;
    }
    if (fitsbin_prepare_mapped_spans(
            ticket->source,
            ticket->planned_ranges,
            range_count,
            ticket->plan_byte_budget,
            ticket->sequence,
            ticket->spans,
            FITSBIN_PREFETCH_RANGE_LIMIT,
            &ticket->span_count,
            &ticket->byte_count,
            &ticket->logical_byte_count,
            &ticket->page_count,
            &ticket->exact_span_count,
            &ticket->reused_page_count,
            &ticket->coalesced_gap_count,
            &ticket->coalesced_gap_bytes)) {
        return -1;
    }
    assert(ticket->byte_count <= ticket->admission_byte_count);
    ticket->range_count = range_count;
#if defined(__linux__)
    if (fitsbin_prepare_mapped_file_spans(
            ticket->source,
            ticket->spans,
            ticket->span_count,
            ticket->plan_byte_budget,
            ticket->queued_spans,
            FITSBIN_PREFETCH_RANGE_LIMIT,
            &ticket->queued_span_count,
            &ticket->queued_byte_count,
            &ticket->queued_exact_span_count,
            &ticket->queued_gap_count,
            &ticket->queued_gap_bytes) ||
        ticket->fd < 0) {
        ticket->queued_span_count = 0U;
        ticket->queued_byte_count = 0U;
        ticket->queued_exact_span_count = 0U;
        ticket->queued_gap_count = 0U;
        ticket->queued_gap_bytes = 0U;
        ticket->queued_prepare_failed = TRUE;
    }
#else
    ticket->queued_prepare_failed = TRUE;
#endif
    return 0;
}

static void* fitsbin_payload_io_service_worker(void* opaque) {
    (void)opaque;
    while (1) {
        fitsbin_payload_io_ticket_t* ticket;
        fitsbin_payload_io_completion_notify_fn completion_notify;
        void* completion_opaque = NULL;
        unsigned long long completion_id;
        struct timespec read_start;
        struct timespec read_finish;
        struct timespec plan_start;
        struct timespec plan_finish;
        unsigned long long plan_nanoseconds = 0ULL;
        unsigned long long read_nanoseconds = 0ULL;
        unsigned long long pread_calls = 0ULL;
        unsigned long long pread_bytes = 0ULL;
        unsigned long long readahead_calls = 0ULL;
        unsigned long long readahead_bytes = 0ULL;
        unsigned long long populate_calls = 0ULL;
        unsigned long long populate_bytes = 0ULL;
        anbool transport_metrics = FALSE;
        anbool plan_measured = FALSE;
        anbool measured;
        anbool cancelled;
        anbool acquired = FALSE;
        size_t work_index;
        int saved_errno = 0;
        int status = 0;

        pthread_mutex_lock(&fitsbin_payload_io_mutex);
        /*
         * Keep unstarted work visible in the priority queues until this lane
         * owns a credit. A newly queued demand ticket can then move ahead of
         * preparation, and queued cancellation can retire without waiting.
         */
        while (1) {
            while (fitsbin_payload_io_queue_empty_locked() &&
                   !fitsbin_payload_io_service_stopping) {
                pthread_cond_wait(
                    &fitsbin_payload_io_cv,
                    &fitsbin_payload_io_mutex);
            }
            if (fitsbin_payload_io_queue_empty_locked() &&
                fitsbin_payload_io_service_stopping) {
                pthread_mutex_unlock(&fitsbin_payload_io_mutex);
                return NULL;
            }
            if (fitsbin_payload_io_try_acquire()) {
                acquired = TRUE;
                break;
            }
            __atomic_add_fetch(
                &fitsbin_payload_io_waiters,
                1U,
                __ATOMIC_RELEASE);
            if (fitsbin_payload_io_try_acquire()) {
                size_t previous = __atomic_fetch_sub(
                    &fitsbin_payload_io_waiters,
                    1U,
                    __ATOMIC_ACQ_REL);

                (void)previous;
                assert(previous > 0U);
                acquired = TRUE;
                break;
            }
            pthread_cond_wait(
                &fitsbin_payload_io_credit_cv,
                &fitsbin_payload_io_mutex);
            {
                size_t previous = __atomic_fetch_sub(
                    &fitsbin_payload_io_waiters,
                    1U,
                    __ATOMIC_ACQ_REL);

                (void)previous;
                assert(previous > 0U);
            }
        }
        ticket = fitsbin_payload_io_dequeue_locked();
        assert(ticket);
        transport_metrics =
            fitsbin_payload_io_transport_metrics_enabled;
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);

        cancelled = __atomic_load_n(
            &ticket->cancel_requested,
            __ATOMIC_ACQUIRE);
        measured = acquired && !cancelled &&
            clock_gettime(CLOCK_MONOTONIC, &read_start) == 0;
        if (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
            for (work_index = 0U;
                 work_index < ticket->range_count && !cancelled;
                 work_index++) {
                fitsbin_prepared_pread_range_t* range =
                    &ticket->ranges[work_index];

                if (fitsbin_pread_all_counted(
                        ticket->fd,
                        range->destination,
                        range->size,
                        range->offset,
                        transport_metrics ? &pread_calls : NULL,
                        transport_metrics ? &pread_bytes : NULL)) {
                    saved_errno = errno ? errno : EIO;
                    status = -1;
                    break;
                }
                cancelled = __atomic_load_n(
                    &ticket->cancel_requested,
                    __ATOMIC_ACQUIRE);
            }
        } else {
#if defined(MADV_POPULATE_READ)
            if (!ticket->source ||
                !ticket->source->mmap_prefetch_enabled ||
                __atomic_load_n(
                    &ticket->source->mmap_prefetch_failed,
                    __ATOMIC_ACQUIRE)) {
                saved_errno = ENOTSUP;
                status = -1;
            }
            if (!status && ticket->plan && !cancelled) {
                int plan_status =
                    0;

                plan_measured = clock_gettime(
                    CLOCK_MONOTONIC, &plan_start) == 0;
                plan_status =
                    fitsbin_payload_io_prepare_planned_ticket(ticket);
                if (plan_measured &&
                    clock_gettime(
                        CLOCK_MONOTONIC, &plan_finish) == 0) {
                    plan_nanoseconds =
                        fitsbin_timespec_delta_nanoseconds(
                            &plan_finish, &plan_start);
                    read_start = plan_finish;
                    measured = TRUE;
                } else {
                    measured = clock_gettime(
                        CLOCK_MONOTONIC, &read_start) == 0;
                }

                if (plan_status < 0) {
                    saved_errno = errno ? errno : EIO;
                    status = -1;
                } else if (plan_status > 0) {
                    cancelled = TRUE;
                }
            }
#if defined(__linux__)
            /*
             * Queue the complete, file-offset-ordered page-cache plan before
             * establishing the exact mapped PTEs. RANDOM remains the VMA
             * policy; this explicit bounded queue supplies storage depth
             * without broad speculative readahead. MADV_POPULATE_READ below
             * remains the completion barrier and authoritative readiness
             * test. Queue failure is advisory and retains the old path.
             */
            for (work_index = 0U;
                 !status &&
                     work_index < ticket->queued_span_count &&
                     !cancelled;
                 work_index++) {
                const fitsbin_file_span_t* span =
                    &ticket->queued_spans[work_index];
                size_t span_bytes =
                    (size_t)(span->end - span->begin);
                int queue_status;

                do {
                    if (transport_metrics) {
                        fitsbin_payload_counter_add(
                            &readahead_calls, 1ULL);
                    }
                    queue_status = readahead(
                        ticket->fd, span->begin, span_bytes);
                } while (queue_status && errno == EINTR);
                if (queue_status) {
                    ticket->queued_submit_failed = TRUE;
                    break;
                }
                ticket->queued_ranges_submitted++;
                ticket->queued_bytes_submitted += span_bytes;
                if (transport_metrics) {
                    fitsbin_payload_counter_add(
                        &readahead_bytes,
                        (unsigned long long)span_bytes);
                }
                cancelled = __atomic_load_n(
                    &ticket->cancel_requested,
                    __ATOMIC_ACQUIRE);
            }
#endif
            for (work_index = 0U;
                 !status && work_index < ticket->span_count &&
                     !cancelled;
                 work_index++) {
                fitsbin_mapped_span_t* span =
                    &ticket->spans[work_index];
                size_t span_bytes =
                    (size_t)(span->end - span->begin);

                int populate_status;

                if (!span_bytes) {
                    errno = EINVAL;
                    populate_status = -1;
                } else {
                    do {
                        if (transport_metrics) {
                            fitsbin_payload_counter_add(
                                &populate_calls, 1ULL);
                        }
                        populate_status = madvise(
                            (void*)span->begin,
                            span_bytes,
                            MADV_POPULATE_READ);
                    } while (populate_status && errno == EINTR);
                }
                if (populate_status) {
                    saved_errno = errno ? errno : EIO;
                    status = -1;
                    if (ticket->source &&
                        (saved_errno == EINVAL ||
                         saved_errno == ENOSYS ||
                         saved_errno == EFAULT ||
                         saved_errno == EACCES)) {
                        __atomic_store_n(
                            &ticket->source->mmap_prefetch_failed,
                            TRUE,
                            __ATOMIC_RELEASE);
                    }
                    break;
                }
                fitsbin_payload_mark_completed_span(
                    ticket->source, span, ticket->sequence);
                if (transport_metrics) {
                    fitsbin_payload_counter_add(
                        &populate_bytes,
                        (unsigned long long)span_bytes);
                }
                cancelled = __atomic_load_n(
                    &ticket->cancel_requested,
                    __ATOMIC_ACQUIRE);
            }
#else
            saved_errno = ENOTSUP;
            status = -1;
#endif
        }
        if (measured &&
            clock_gettime(CLOCK_MONOTONIC, &read_finish) == 0) {
                read_nanoseconds =
                    fitsbin_timespec_delta_nanoseconds(
                        &read_finish,
                        &read_start);
        }
        if (acquired) {
            fitsbin_payload_io_release();
        }
        if (ticket->fd >= 0) {
            close(ticket->fd);
            ticket->fd = -1;
        }

        pthread_mutex_lock(&fitsbin_payload_io_mutex);
        assert(fitsbin_payload_io_service_jobs > 0U);
        assert(fitsbin_payload_io_service_bytes >=
               ticket->admission_byte_count);
        fitsbin_payload_io_service_jobs--;
        fitsbin_payload_io_service_bytes -=
            ticket->admission_byte_count;
        if (ticket->plan) {
            fitsbin_payload_io_service_planned_bytes +=
                (unsigned long long)ticket->byte_count;
            fitsbin_payload_io_service_plan_nanoseconds +=
                plan_nanoseconds;
        }
        ticket->read_nanoseconds = read_nanoseconds;
        fitsbin_payload_io_service_read_nanoseconds += read_nanoseconds;
        if (transport_metrics) {
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.pread_calls,
                pread_calls);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.pread_bytes,
                pread_bytes);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.readahead_calls,
                readahead_calls);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.readahead_bytes,
                readahead_bytes);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.populate_calls,
                populate_calls);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.populate_bytes,
                populate_bytes);
            if (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.direct_ranges,
                    (unsigned long long)ticket->range_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.direct_bytes,
                    (unsigned long long)ticket->byte_count);
            } else {
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.mapped_spans,
                    (unsigned long long)ticket->span_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        mapped_exact_spans,
                    (unsigned long long)ticket->exact_span_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        mapped_reused_pages,
                    ticket->reused_page_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        mapped_gap_merges,
                    (unsigned long long)ticket->coalesced_gap_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        mapped_gap_bytes,
                    (unsigned long long)ticket->coalesced_gap_bytes);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.mapped_bytes,
                    (unsigned long long)ticket->byte_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        queued_exact_spans,
                    (unsigned long long)
                        ticket->queued_exact_span_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        queued_gap_merges,
                    (unsigned long long)ticket->queued_gap_count);
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.
                        queued_gap_bytes,
                    (unsigned long long)ticket->queued_gap_bytes);
            }
        }
        fitsbin_payload_io_service_queued_ranges +=
            (unsigned long long)ticket->queued_ranges_submitted;
        fitsbin_payload_io_service_queued_bytes +=
            (unsigned long long)ticket->queued_bytes_submitted;
        if (ticket->queued_prepare_failed) {
            fitsbin_payload_io_service_queue_prepare_failures++;
        }
        if (ticket->queued_submit_failed) {
            fitsbin_payload_io_service_queue_submit_failures++;
        }
        if (__atomic_load_n(
                &ticket->cancel_requested,
                __ATOMIC_ACQUIRE) ||
            cancelled) {
            ticket->state = FITSBIN_PAYLOAD_IO_CANCELLED;
            fitsbin_payload_io_service_cancelled++;
        } else if (status) {
            ticket->saved_errno = saved_errno;
            ticket->state = FITSBIN_PAYLOAD_IO_FAILED;
            fitsbin_payload_io_service_failed++;
        } else {
            ticket->state = FITSBIN_PAYLOAD_IO_READY;
            fitsbin_payload_io_service_ready++;
            if (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
                fitsbin_payload_io_service_direct_ready++;
            } else if (fitsbin_payload_io_transport_metrics_enabled) {
                fitsbin_payload_counter_add(
                    &fitsbin_payload_io_transport_metrics.mapped_ready,
                    1ULL);
            }
        }
        pthread_cond_broadcast(&ticket->completion_cv);
        completion_notify =
            fitsbin_payload_io_completion_acquire_locked(
                &completion_opaque);
        completion_id = ticket->sequence;
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        fitsbin_payload_io_completion_dispatch(
            completion_notify,
            completion_opaque,
            completion_id);
    }
    return NULL;
}

int fitsbin_payload_io_service_start(int lane_count) {
    int created = 0;
    int status = 0;

    if (lane_count <= 0) {
        errno = EINVAL;
        return -1;
    }
    lane_count = MIN(lane_count, FITSBIN_PAYLOAD_IO_MAX_LANES);

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_service_running) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return 0;
    }
    fitsbin_payload_io_transport_metrics_enabled =
        log_get_level() >= LOG_VERB;
    fitsbin_payload_io_service_running = TRUE;
    fitsbin_payload_io_service_accepting = TRUE;
    fitsbin_payload_io_service_stopping = FALSE;
    fitsbin_payload_io_service_lanes = lane_count;
    fitsbin_payload_io_service_submitted = 0ULL;
    fitsbin_payload_io_service_ready = 0ULL;
    fitsbin_payload_io_service_direct_submitted = 0ULL;
    fitsbin_payload_io_service_direct_ready = 0ULL;
    memset(&fitsbin_payload_io_transport_metrics, 0,
           sizeof(fitsbin_payload_io_transport_metrics));
    fitsbin_payload_io_service_cancelled = 0ULL;
    fitsbin_payload_io_service_failed = 0ULL;
    fitsbin_payload_io_service_planned_bytes = 0ULL;
    fitsbin_payload_io_service_plan_nanoseconds = 0ULL;
    fitsbin_payload_io_service_read_nanoseconds = 0ULL;
    fitsbin_payload_io_service_queued_ranges = 0ULL;
    fitsbin_payload_io_service_queued_bytes = 0ULL;
    fitsbin_payload_io_service_queue_prepare_failures = 0ULL;
    fitsbin_payload_io_service_queue_submit_failures = 0ULL;
    while (created < lane_count) {
        status = pthread_create(
            &fitsbin_payload_io_threads[created],
            NULL,
            fitsbin_payload_io_service_worker,
            NULL);
        if (status) {
            break;
        }
        created++;
    }
    if (!status) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return 0;
    }
    fitsbin_payload_io_service_accepting = FALSE;
    fitsbin_payload_io_service_stopping = TRUE;
    fitsbin_payload_io_service_lanes = created;
    pthread_cond_broadcast(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    while (created > 0) {
        created--;
        pthread_join(fitsbin_payload_io_threads[created], NULL);
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_service_running = FALSE;
    fitsbin_payload_io_transport_metrics_enabled = FALSE;
    fitsbin_payload_io_service_stopping = FALSE;
    fitsbin_payload_io_service_lanes = 0;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    errno = status;
    return -1;
}

void fitsbin_payload_io_service_stop(void) {
#if FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED
    const char* transport_policy = "sparse-direct-enabled";
#else
    const char* transport_policy = "sparse-mapped-only";
#endif
    int lane_count;
    int lane;

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (!fitsbin_payload_io_service_running) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return;
    }
    fitsbin_payload_io_service_accepting = FALSE;
    fitsbin_payload_io_service_stopping = TRUE;
    lane_count = fitsbin_payload_io_service_lanes;
    pthread_cond_broadcast(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    for (lane = 0; lane < lane_count; lane++) {
        pthread_join(fitsbin_payload_io_threads[lane], NULL);
    }

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    assert(fitsbin_payload_io_queue_empty_locked());
    assert(!fitsbin_payload_io_service_jobs);
    assert(!fitsbin_payload_io_service_bytes);
    logverb("[fitsbin] payload-broker submitted=%llu ready=%llu "
            "direct_submitted=%llu direct_ready=%llu "
            "mapped_submitted=%llu mapped_ready=%llu "
            "cancelled=%llu failed=%llu planned_bytes=%llu "
            "queued_ranges=%llu queued_bytes=%llu "
            "queue_failures=%llu/%llu "
            "plan_wall=%.6f read_wall=%.6f\n",
            fitsbin_payload_io_service_submitted,
            fitsbin_payload_io_service_ready,
            fitsbin_payload_io_service_direct_submitted,
            fitsbin_payload_io_service_direct_ready,
            fitsbin_payload_io_transport_metrics.mapped_submitted,
            fitsbin_payload_io_transport_metrics.mapped_ready,
            fitsbin_payload_io_service_cancelled,
            fitsbin_payload_io_service_failed,
            fitsbin_payload_io_service_planned_bytes,
            fitsbin_payload_io_service_queued_ranges,
            fitsbin_payload_io_service_queued_bytes,
            fitsbin_payload_io_service_queue_prepare_failures,
            fitsbin_payload_io_service_queue_submit_failures,
            (double)fitsbin_payload_io_service_plan_nanoseconds /
                1000000000.0,
            (double)fitsbin_payload_io_service_read_nanoseconds /
                1000000000.0);
    logverb("[fitsbin] payload-transport policy=%s "
            "direct_attempts=%llu direct_policy_refused=%llu "
            "direct_ranges=%llu direct_bytes=%llu "
            "mapped_immediate_ready=%llu "
            "mapped_spans=%llu mapped_exact_spans=%llu "
            "mapped_reused_pages=%llu "
            "mapped_gap_merges=%llu mapped_gap_bytes=%llu "
            "mapped_bytes=%llu "
            "queue_exact_spans=%llu queue_gap_merges=%llu "
            "queue_gap_bytes=%llu "
            "pread_calls=%llu pread_bytes=%llu "
            "readahead_calls=%llu readahead_bytes=%llu "
            "populate_calls=%llu populate_bytes=%llu\n",
            transport_policy,
            fitsbin_payload_io_transport_metrics.direct_attempts,
            fitsbin_payload_io_transport_metrics.direct_policy_refused,
            fitsbin_payload_io_transport_metrics.direct_ranges,
            fitsbin_payload_io_transport_metrics.direct_bytes,
            fitsbin_payload_io_transport_metrics.mapped_immediate_ready,
            fitsbin_payload_io_transport_metrics.mapped_spans,
            fitsbin_payload_io_transport_metrics.mapped_exact_spans,
            fitsbin_payload_io_transport_metrics.mapped_reused_pages,
            fitsbin_payload_io_transport_metrics.mapped_gap_merges,
            fitsbin_payload_io_transport_metrics.mapped_gap_bytes,
            fitsbin_payload_io_transport_metrics.mapped_bytes,
            fitsbin_payload_io_transport_metrics.queued_exact_spans,
            fitsbin_payload_io_transport_metrics.queued_gap_merges,
            fitsbin_payload_io_transport_metrics.queued_gap_bytes,
            fitsbin_payload_io_transport_metrics.pread_calls,
            fitsbin_payload_io_transport_metrics.pread_bytes,
            fitsbin_payload_io_transport_metrics.readahead_calls,
            fitsbin_payload_io_transport_metrics.readahead_bytes,
            fitsbin_payload_io_transport_metrics.populate_calls,
            fitsbin_payload_io_transport_metrics.populate_bytes);
    fitsbin_payload_io_service_running = FALSE;
    fitsbin_payload_io_transport_metrics_enabled = FALSE;
    fitsbin_payload_io_service_stopping = FALSE;
    fitsbin_payload_io_service_lanes = 0;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
}

int fitsbin_payload_io_service_width(void) {
    int lane_count;

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    lane_count = fitsbin_payload_io_service_running &&
        fitsbin_payload_io_service_accepting &&
        !fitsbin_payload_io_service_stopping
        ? fitsbin_payload_io_service_lanes
        : 0;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    return lane_count;
}

int fitsbin_payload_io_mapped_population_supported(void) {
#if defined(MADV_POPULATE_READ)
    return TRUE;
#else
    return FALSE;
#endif
}

static int fitsbin_payload_io_submit_ticket(
    fitsbin_payload_io_ticket_t* ticket,
    fitsbin_payload_io_ticket_t** ticket_out) {
    anbool service_available;
    int admitted = 0;

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    service_available =
        fitsbin_payload_io_service_running &&
        fitsbin_payload_io_service_accepting;
    if (service_available) {
        admitted = fitsbin_payload_io_admit_locked(ticket);
    }
    if (!admitted) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        fitsbin_payload_io_ticket_free_storage(ticket);
        errno = service_available ? EAGAIN : ENODEV;
        return 0;
    }
    if (admitted < 0) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        fitsbin_payload_io_ticket_free_storage(ticket);
        errno = E2BIG;
        return 0;
    }
    /*
     * Completion IDs are value tokens that may outlive ticket storage while
     * a notifier is draining. Never recycle an ID after counter exhaustion.
     */
    if (fitsbin_payload_io_next_sequence == ULLONG_MAX) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        fitsbin_payload_io_ticket_free_storage(ticket);
        errno = EOVERFLOW;
        return -1;
    }
    ticket->sequence = ++fitsbin_payload_io_next_sequence;
    ticket->state = FITSBIN_PAYLOAD_IO_SUBMITTED;
    fitsbin_payload_io_enqueue_locked(ticket);
    fitsbin_payload_io_service_jobs++;
    fitsbin_payload_io_service_bytes +=
        ticket->admission_byte_count;
    fitsbin_payload_io_service_submitted++;
    if (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
        fitsbin_payload_io_service_direct_submitted++;
    } else if (fitsbin_payload_io_transport_metrics_enabled) {
        fitsbin_payload_counter_add(
            &fitsbin_payload_io_transport_metrics.mapped_submitted,
            1ULL);
    }
    if (!ticket->plan) {
        fitsbin_payload_io_service_planned_bytes +=
            (unsigned long long)ticket->byte_count;
    }
    *ticket_out = ticket;
    pthread_cond_signal(&fitsbin_payload_io_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    return 1;
}

int fitsbin_pread_mapped_ranges_submit(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_payload_io_priority_t priority,
    fitsbin_payload_io_ticket_t** ticket_out) {
    fitsbin_payload_io_ticket_t* ticket;
    int fd;
    int duplicate;

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!ticket_out) {
        errno = EINVAL;
        return -1;
    }
    *ticket_out = NULL;
    if (!fb || !ranges || !range_count || !byte_budget ||
        !fitsbin_payload_io_priority_valid(priority)) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREAD_ASYNC_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (fitsbin_payload_io_transport_metrics_enabled) {
        fitsbin_payload_counter_add(
            &fitsbin_payload_io_transport_metrics.direct_attempts,
            1ULL);
    }
    if (!fitsbin_payload_io_service_running ||
        !fitsbin_payload_io_service_accepting) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = ENODEV;
        return 0;
    }
#if !FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED
    if (fitsbin_get_mmap_advice(fb) == FITSBIN_MMAP_ADVICE_RANDOM &&
        !fitsbin_payload_is_fully_resident(fb)) {
        if (fitsbin_payload_io_transport_metrics_enabled) {
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.
                    direct_policy_refused,
                1ULL);
        }
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = ENOTSUP;
        return 0;
    }
#endif
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    ticket = fitsbin_payload_io_ticket_alloc();
    if (!ticket) {
        return -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_PLANNED;
    ticket->kind = FITSBIN_PAYLOAD_IO_TICKET_DIRECT;
    ticket->priority = priority;
    ticket->source = fb;
    if (range_count > SIZE_MAX / sizeof(ticket->ranges[0])) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        errno = EOVERFLOW;
        return -1;
    }
    ticket->ranges = calloc(
        range_count, sizeof(ticket->ranges[0]));
    if (!ticket->ranges) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return -1;
    }
    if (fitsbin_prepare_direct_ranges(
            fb,
            ranges,
            range_count,
            byte_budget,
            ticket->ranges,
            &ticket->byte_count,
            &ticket->logical_byte_count,
            &ticket->page_count)) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return -1;
    }
    ticket->range_count = range_count;
    ticket->admission_byte_count = ticket->byte_count;
    fd = fitsbin_payload_fd_get(fb);
    if (fd < 0) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return -1;
    }
    duplicate = fitsbin_payload_io_duplicate_fd(fd);
    if (duplicate < 0) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return -1;
    }
    ticket->fd = duplicate;
    return fitsbin_payload_io_submit_ticket(ticket, ticket_out);
}

int fitsbin_prefetch_ranges_submit(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_payload_io_ticket_t** ticket_out) {
    fitsbin_payload_io_ticket_t* ticket;
    fitsbin_mapped_span_t
        prepared_spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    size_t prepared_span_count;
    size_t prepared_byte_count;
    size_t prepared_logical_byte_count;
    unsigned long long prepared_page_count;
    size_t prepared_exact_span_count;
    unsigned long long prepared_reused_page_count;
    size_t prepared_coalesced_gap_count;
    size_t prepared_coalesced_gap_bytes;
#if defined(__linux__)
    int fd;
    int duplicate;
#endif

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!ticket_out) {
        errno = EINVAL;
        return -1;
    }
    *ticket_out = NULL;
    if (!fb || !ranges || !range_count || !byte_budget) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREFETCH_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    if (fitsbin_payload_is_fully_resident(fb)) {
        errno = 0;
        return 0;
    }
    if (!fb->mmap_prefetch_enabled ||
        __atomic_load_n(
            &fb->mmap_prefetch_failed,
            __ATOMIC_ACQUIRE)) {
        errno = 0;
        return 0;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (!fitsbin_payload_io_service_running ||
        !fitsbin_payload_io_service_accepting) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = ENODEV;
        return 0;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    if (fitsbin_prepare_mapped_spans(
            fb,
            ranges,
            range_count,
            byte_budget,
            fitsbin_payload_io_sequence_hint(),
            prepared_spans,
            FITSBIN_PREFETCH_RANGE_LIMIT,
            &prepared_span_count,
            &prepared_byte_count,
            &prepared_logical_byte_count,
            &prepared_page_count,
            &prepared_exact_span_count,
            &prepared_reused_page_count,
            &prepared_coalesced_gap_count,
            &prepared_coalesced_gap_bytes)) {
        return -1;
    }
    if (!prepared_span_count) {
        anbool service_available;

        pthread_mutex_lock(&fitsbin_payload_io_mutex);
        service_available =
            fitsbin_payload_io_service_running &&
            fitsbin_payload_io_service_accepting;
        if (service_available &&
            fitsbin_payload_io_transport_metrics_enabled) {
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.
                    mapped_immediate_ready,
                1ULL);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.
                    mapped_exact_spans,
                (unsigned long long)prepared_exact_span_count);
            fitsbin_payload_counter_add(
                &fitsbin_payload_io_transport_metrics.
                    mapped_reused_pages,
                prepared_reused_page_count);
        }
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        if (!service_available) {
            errno = ENODEV;
            return FITSBIN_PAYLOAD_IO_SUBMIT_UNAVAILABLE;
        }
        errno = 0;
        return FITSBIN_PAYLOAD_IO_SUBMIT_READY;
    }
    ticket = fitsbin_payload_io_ticket_alloc();
    if (!ticket) {
        return -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_PLANNED;
    ticket->kind = FITSBIN_PAYLOAD_IO_TICKET_PREFETCH;
    ticket->priority = FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT;
    ticket->source = fb;
    memcpy(
        ticket->spans,
        prepared_spans,
        prepared_span_count * sizeof(*prepared_spans));
    ticket->span_count = prepared_span_count;
    ticket->byte_count = prepared_byte_count;
    ticket->logical_byte_count = prepared_logical_byte_count;
    ticket->page_count = prepared_page_count;
    ticket->exact_span_count = prepared_exact_span_count;
    ticket->reused_page_count = prepared_reused_page_count;
    ticket->coalesced_gap_count = prepared_coalesced_gap_count;
    ticket->coalesced_gap_bytes = prepared_coalesced_gap_bytes;
    ticket->range_count = range_count;
    ticket->admission_byte_count = ticket->byte_count;
#if defined(__linux__)
    if (fitsbin_prepare_mapped_file_spans(
            fb,
            ticket->spans,
            ticket->span_count,
            byte_budget,
            ticket->queued_spans,
            FITSBIN_PREFETCH_RANGE_LIMIT,
            &ticket->queued_span_count,
            &ticket->queued_byte_count,
            &ticket->queued_exact_span_count,
            &ticket->queued_gap_count,
            &ticket->queued_gap_bytes)) {
        ticket->queued_span_count = 0U;
        ticket->queued_byte_count = 0U;
        ticket->queued_exact_span_count = 0U;
        ticket->queued_gap_count = 0U;
        ticket->queued_gap_bytes = 0U;
        ticket->queued_prepare_failed = TRUE;
    } else {
        fd = fitsbin_payload_fd_get(fb);
        duplicate = fd >= 0
            ? fitsbin_payload_io_duplicate_fd(fd)
            : -1;
        if (duplicate < 0) {
            ticket->queued_span_count = 0U;
            ticket->queued_byte_count = 0U;
            ticket->queued_exact_span_count = 0U;
            ticket->queued_gap_count = 0U;
            ticket->queued_gap_bytes = 0U;
            ticket->queued_prepare_failed = TRUE;
        } else {
            ticket->fd = duplicate;
            ticket->admission_byte_count = MAX(
                ticket->admission_byte_count,
                ticket->queued_byte_count);
        }
    }
#else
    ticket->queued_prepare_failed = TRUE;
#endif
    errno = 0;
    return fitsbin_payload_io_submit_ticket(ticket, ticket_out);
}

int fitsbin_prefetch_ranges_planned_submit(
    fitsbin_t* fb,
    fitsbin_payload_io_plan_fn plan,
    void* plan_opaque,
    size_t byte_budget,
    fitsbin_payload_io_ticket_t** ticket_out) {
    fitsbin_payload_io_ticket_t* ticket;
#if defined(__linux__)
    int fd;
    int duplicate;
#endif

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!ticket_out) {
        errno = EINVAL;
        return -1;
    }
    *ticket_out = NULL;
    if (!fb || !plan || !byte_budget) {
        errno = EINVAL;
        return -1;
    }
    if (fitsbin_payload_is_fully_resident(fb)) {
        errno = 0;
        return 0;
    }
    if (!fb->mmap_prefetch_enabled ||
        __atomic_load_n(
            &fb->mmap_prefetch_failed,
            __ATOMIC_ACQUIRE)) {
        errno = 0;
        return 0;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (!fitsbin_payload_io_service_running ||
        !fitsbin_payload_io_service_accepting) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = ENODEV;
        return 0;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    ticket = fitsbin_payload_io_ticket_alloc();
    if (!ticket) {
        return -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_PLANNED;
    ticket->kind = FITSBIN_PAYLOAD_IO_TICKET_PREFETCH;
    ticket->priority = FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT;
    ticket->source = fb;
    ticket->plan = plan;
    ticket->plan_opaque = plan_opaque;
    ticket->plan_byte_budget = byte_budget;
    ticket->admission_byte_count = byte_budget;
#if defined(__linux__)
    fd = fitsbin_payload_fd_get(fb);
    duplicate = fd >= 0
        ? fitsbin_payload_io_duplicate_fd(fd)
        : -1;
    if (duplicate < 0) {
        ticket->queued_prepare_failed = TRUE;
    } else {
        ticket->fd = duplicate;
    }
#else
    ticket->queued_prepare_failed = TRUE;
#endif
    errno = 0;
    return fitsbin_payload_io_submit_ticket(ticket, ticket_out);
}

typedef struct fitsbin_payload_io_ticket_result {
    unsigned long long read_nanoseconds;
    size_t span_count;
    size_t range_count;
    size_t byte_count;
    size_t logical_byte_count;
    unsigned long long page_count;
    fitsbin_payload_io_ticket_state_t state;
    fitsbin_payload_io_ticket_kind_t kind;
    int saved_errno;
} fitsbin_payload_io_ticket_result_t;

static anbool fitsbin_payload_io_ticket_terminal(
    fitsbin_payload_io_ticket_state_t state) {
    return state == FITSBIN_PAYLOAD_IO_READY ||
        state == FITSBIN_PAYLOAD_IO_FAILED ||
        state == FITSBIN_PAYLOAD_IO_CANCELLED;
}

/* fitsbin_payload_io_mutex must be held. */
static void fitsbin_payload_io_ticket_snapshot_locked(
    const fitsbin_payload_io_ticket_t* ticket,
    fitsbin_payload_io_ticket_result_t* result) {
    assert(ticket);
    assert(result);
    result->read_nanoseconds = ticket->read_nanoseconds;
    result->span_count = ticket->span_count;
    result->range_count = ticket->range_count;
    result->byte_count = ticket->byte_count;
    result->logical_byte_count = ticket->logical_byte_count;
    result->page_count = ticket->page_count;
    result->state = ticket->state;
    result->kind = ticket->kind;
    result->saved_errno = ticket->saved_errno;
}

/* fitsbin_payload_io_mutex must be held. */
static int fitsbin_payload_io_ticket_collect_locked(
    fitsbin_payload_io_ticket_t* ticket,
    fitsbin_payload_io_ticket_result_t* result) {
    if (!ticket || !result) {
        errno = EINVAL;
        return -1;
    }
    if (ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        return 0;
    }
    if (!fitsbin_payload_io_ticket_terminal(ticket->state)) {
        errno = EINVAL;
        return -1;
    }
    if (ticket->counters_applied) {
        errno = EALREADY;
        return -1;
    }
    fitsbin_payload_io_ticket_snapshot_locked(ticket, result);
    ticket->counters_applied = TRUE;
    return 1;
}

static int fitsbin_payload_io_ticket_result_status(
    const fitsbin_payload_io_ticket_result_t* result) {
    assert(result);
    if (result->state == FITSBIN_PAYLOAD_IO_READY &&
        result->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
        errno = 0;
        return (int)result->range_count;
    }
    if (result->state == FITSBIN_PAYLOAD_IO_READY &&
        result->span_count) {
        errno = 0;
        return (int)result->span_count;
    }
    if (result->state == FITSBIN_PAYLOAD_IO_READY) {
        errno = 0;
        return 1;
    }
    if (result->state == FITSBIN_PAYLOAD_IO_FAILED) {
        errno = result->saved_errno ? result->saved_errno : EIO;
        return -1;
    }
    errno = ECANCELED;
    return 0;
}

static int fitsbin_payload_io_ticket_apply_result(
    fitsbin_t* fb,
    const fitsbin_payload_io_ticket_result_t* result,
    unsigned long long waited) {
    assert(fb);
    assert(result);

    __atomic_add_fetch(
        &fb->payload_wait_nanoseconds,
        waited,
        __ATOMIC_RELAXED);
    if (result->state == FITSBIN_PAYLOAD_IO_READY &&
        result->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
        __atomic_add_fetch(
            &fb->payload_read_batches,
            1ULL,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_calls,
            (unsigned long long)result->range_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_bytes,
            (unsigned long long)result->byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_logical_bytes,
            (unsigned long long)result->logical_byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_pages,
            result->page_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_nanoseconds,
            result->read_nanoseconds,
            __ATOMIC_RELAXED);
        return fitsbin_payload_io_ticket_result_status(result);
    }
    if (result->state == FITSBIN_PAYLOAD_IO_READY &&
        result->span_count) {
        __atomic_add_fetch(
            &fb->payload_warm_calls,
            1ULL,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_ranges,
            (unsigned long long)result->span_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_bytes,
            (unsigned long long)result->byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_nanoseconds,
            result->read_nanoseconds,
            __ATOMIC_RELAXED);
        return fitsbin_payload_io_ticket_result_status(result);
    }
    if (result->state == FITSBIN_PAYLOAD_IO_FAILED) {
        __atomic_add_fetch(
            &fb->payload_failures,
            1ULL,
            __ATOMIC_RELAXED);
    }
    return fitsbin_payload_io_ticket_result_status(result);
}

static int fitsbin_payload_io_ticket_wait_internal(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket,
    anbool cancel) {
    struct timespec wait_start;
    struct timespec wait_finish;
    fitsbin_payload_io_ticket_result_t result;
    fitsbin_payload_io_completion_notify_fn completion_notify = NULL;
    void* completion_opaque = NULL;
    unsigned long long completion_id;
    unsigned long long waited = 0ULL;
    anbool helper_registered = FALSE;
    anbool waiter_registered = FALSE;
    anbool retry_helper = TRUE;
    anbool measured;
    int result_errno;
    int result_value = -1;
    int collected;

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!fb || !ticket || ticket->source != fb) {
        errno = EINVAL;
        return -1;
    }
    measured = clock_gettime(CLOCK_MONOTONIC, &wait_start) == 0;
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    completion_id = ticket->sequence;
    if (ticket->owner_draining || ticket->waiter_active) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    ticket->waiter_active = TRUE;
    if (cancel) {
        __atomic_store_n(
            &ticket->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        completion_notify =
            fitsbin_payload_io_cancel_queued_locked(
                ticket,
                &completion_opaque);
    }
    if (ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        if (fitsbin_payload_io_register_waiter_locked(ticket)) {
            ticket->waiter_active = FALSE;
            pthread_cond_broadcast(&ticket->completion_cv);
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            errno = EINVAL;
            return -1;
        }
        waiter_registered = TRUE;
    }
    if (fitsbin_payload_io_thread_wait_helper &&
        !fitsbin_payload_io_thread_wait_active &&
        ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        __atomic_add_fetch(
            &fitsbin_payload_io_wait_helpers,
            1U,
            __ATOMIC_RELEASE);
        helper_registered = TRUE;
        ticket->helper_waiter = TRUE;
    }
    while (ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        if (fitsbin_payload_io_thread_stop_check &&
            fitsbin_payload_io_thread_stop_check(
                fitsbin_payload_io_thread_wait_opaque)) {
            __atomic_store_n(
                &ticket->cancel_requested,
                TRUE,
                __ATOMIC_RELEASE);
            completion_notify =
                fitsbin_payload_io_cancel_queued_locked(
                    ticket,
                    &completion_opaque);
        }
        if (ticket->state != FITSBIN_PAYLOAD_IO_SUBMITTED) {
            continue;
        }
        if (fitsbin_payload_io_thread_wait_helper &&
            !__atomic_load_n(
                &ticket->cancel_requested,
                __ATOMIC_ACQUIRE) &&
            !fitsbin_payload_io_thread_wait_active &&
            (retry_helper ||
             fitsbin_payload_io_thread_work_epoch !=
                 __atomic_load_n(
                     &fitsbin_payload_io_work_epoch,
                     __ATOMIC_ACQUIRE))) {
            fitsbin_payload_io_wait_helper_fn helper =
                fitsbin_payload_io_thread_wait_helper;
            void* helper_opaque =
                fitsbin_payload_io_thread_wait_opaque;
            int helped;

            fitsbin_payload_io_thread_work_epoch =
                __atomic_load_n(
                    &fitsbin_payload_io_work_epoch,
                    __ATOMIC_ACQUIRE);
            fitsbin_payload_io_thread_wait_active = TRUE;
            __atomic_add_fetch(
                &fitsbin_payload_io_wait_helpers_active,
                1U,
                __ATOMIC_ACQ_REL);
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            helped = helper(helper_opaque);
            pthread_mutex_lock(&fitsbin_payload_io_mutex);
            __atomic_sub_fetch(
                &fitsbin_payload_io_wait_helpers_active,
                1U,
                __ATOMIC_ACQ_REL);
            fitsbin_payload_io_thread_wait_active = FALSE;
            retry_helper = helped != 0;
            if (retry_helper) {
                continue;
            }
        }
        if (ticket->state != FITSBIN_PAYLOAD_IO_SUBMITTED) {
            continue;
        }
        pthread_cond_wait(
            &ticket->completion_cv,
            &fitsbin_payload_io_mutex);
    }
    if (waiter_registered) {
        int unregister_status =
            fitsbin_payload_io_unregister_waiter_locked(ticket);

        (void)unregister_status;
        assert(!unregister_status);
    }
    if (helper_registered) {
        size_t previous = __atomic_fetch_sub(
            &fitsbin_payload_io_wait_helpers,
            1U,
            __ATOMIC_ACQ_REL);

        (void)previous;
        assert(previous > 0U);
    }
    collected = fitsbin_payload_io_ticket_collect_locked(
        ticket,
        &result);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    if (collected == 1) {
        if (measured &&
            clock_gettime(CLOCK_MONOTONIC, &wait_finish) == 0) {
            waited = fitsbin_timespec_delta_nanoseconds(
                &wait_finish,
                &wait_start);
        }
        result_value = fitsbin_payload_io_ticket_apply_result(
            fb,
            &result,
            waited);
    }
    result_errno = errno;
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    ticket->waiter_active = FALSE;
    pthread_cond_broadcast(&ticket->completion_cv);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_completion_dispatch(
        completion_notify,
        completion_opaque,
        completion_id);
    errno = result_errno;
    return collected == 1 ? result_value : -1;
}

int fitsbin_payload_io_ticket_wait(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket) {
    return fitsbin_payload_io_ticket_wait_internal(
        fb, ticket, FALSE);
}

int fitsbin_payload_io_ticket_cancel_and_wait(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket) {
    return fitsbin_payload_io_ticket_wait_internal(
        fb, ticket, TRUE);
}

int fitsbin_payload_io_ticket_poll(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket,
    int* result_out) {
    fitsbin_payload_io_ticket_result_t result;
    int collected;

    if (!fb || !ticket || !result_out || ticket->source != fb) {
        errno = EINVAL;
        return -1;
    }
    *result_out = 0;
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (ticket->owner_draining || ticket->waiter_active) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    collected = fitsbin_payload_io_ticket_collect_locked(
        ticket,
        &result);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    if (collected <= 0) {
        return collected;
    }
    *result_out = fitsbin_payload_io_ticket_apply_result(
        fb,
        &result,
        0ULL);
    return 1;
}

int fitsbin_payload_io_ticket_poll_and_destroy(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t** ticket_io,
    int* result_out) {
    fitsbin_payload_io_ticket_result_t result;
    fitsbin_payload_io_ticket_t* ticket;
    int result_errno;
    int collected;

    if (!fb || !ticket_io || !*ticket_io || !result_out) {
        errno = EINVAL;
        return -1;
    }
    ticket = *ticket_io;
    if (ticket->source != fb) {
        errno = EINVAL;
        return -1;
    }
    *result_out = 0;
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (ticket->owner_draining || ticket->waiter_active) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    collected = fitsbin_payload_io_ticket_collect_locked(
        ticket,
        &result);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    if (collected <= 0) {
        return collected;
    }

    *result_out = fitsbin_payload_io_ticket_apply_result(
        fb,
        &result,
        0ULL);
    result_errno = errno;
    *ticket_io = NULL;
    fitsbin_payload_io_ticket_free_storage(ticket);
    errno = result_errno;
    return 1;
}

int fitsbin_payload_io_ticket_drain_and_destroy(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t** ticket_io,
    int* result_out) {
    struct timespec wait_start;
    struct timespec wait_finish;
    fitsbin_payload_io_ticket_result_t result;
    fitsbin_payload_io_ticket_t* ticket;
    fitsbin_payload_io_completion_notify_fn completion_notify = NULL;
    void* completion_opaque = NULL;
    unsigned long long completion_id;
    unsigned long long waited = 0ULL;
    anbool apply_counters = FALSE;
    anbool measured;
    int result_errno;
    int result_value;
    int status;

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!fb || !ticket_io || !*ticket_io || !result_out) {
        errno = EINVAL;
        return -1;
    }
    ticket = *ticket_io;
    measured = clock_gettime(CLOCK_MONOTONIC, &wait_start) == 0;
    *result_out = 0;

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (ticket->source != fb) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EINVAL;
        return -1;
    }
    if (ticket->owner_draining) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    ticket->owner_draining = TRUE;
    completion_id = ticket->sequence;
    if (ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED &&
        !ticket->waiter_active) {
        __atomic_store_n(
            &ticket->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        completion_notify =
            fitsbin_payload_io_cancel_queued_locked(
                ticket,
                &completion_opaque);
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_completion_dispatch(
        completion_notify,
        completion_opaque,
        completion_id);

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    while (ticket->waiter_active ||
           ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        status = pthread_cond_wait(
            &ticket->completion_cv,
            &fitsbin_payload_io_mutex);
        if (status) {
            ticket->owner_draining = FALSE;
            pthread_cond_broadcast(&ticket->completion_cv);
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            errno = status;
            return -1;
        }
    }
    if (!fitsbin_payload_io_ticket_terminal(ticket->state)) {
        ticket->owner_draining = FALSE;
        pthread_cond_broadcast(&ticket->completion_cv);
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EINVAL;
        return -1;
    }
    if (ticket->counters_applied) {
        fitsbin_payload_io_ticket_snapshot_locked(ticket, &result);
    } else {
        status = fitsbin_payload_io_ticket_collect_locked(
            ticket,
            &result);
        if (status != 1) {
            ticket->owner_draining = FALSE;
            pthread_cond_broadcast(&ticket->completion_cv);
            pthread_mutex_unlock(&fitsbin_payload_io_mutex);
            return -1;
        }
        apply_counters = TRUE;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    if (apply_counters) {
        if (measured &&
            clock_gettime(CLOCK_MONOTONIC, &wait_finish) == 0) {
            waited = fitsbin_timespec_delta_nanoseconds(
                &wait_finish,
                &wait_start);
        }
        result_value = fitsbin_payload_io_ticket_apply_result(
            fb,
            &result,
            waited);
    } else {
        result_value =
            fitsbin_payload_io_ticket_result_status(&result);
    }
    result_errno = errno;
    *result_out = result_value;
    *ticket_io = NULL;
    fitsbin_payload_io_ticket_free_storage(ticket);
    errno = result_errno;
    return 1;
}

int fitsbin_payload_io_ticket_cancel_async(
    fitsbin_payload_io_ticket_t* ticket) {
    fitsbin_payload_io_completion_notify_fn completion_notify = NULL;
    void* completion_opaque = NULL;
    unsigned long long completion_id;
    anbool already_requested;

    if (!ticket) {
        errno = EINVAL;
        return -1;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (ticket->owner_draining) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    if (fitsbin_payload_io_ticket_terminal(ticket->state)) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return 0;
    }
    if (ticket->state != FITSBIN_PAYLOAD_IO_SUBMITTED) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EINVAL;
        return -1;
    }
    completion_id = ticket->sequence;
    already_requested = __atomic_exchange_n(
        &ticket->cancel_requested,
        TRUE,
        __ATOMIC_ACQ_REL);
    completion_notify =
        fitsbin_payload_io_cancel_queued_locked(
            ticket,
            &completion_opaque);
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_completion_dispatch(
        completion_notify,
        completion_opaque,
        completion_id);
    return already_requested ? 0 : 1;
}

unsigned long long fitsbin_payload_io_ticket_completion_id(
    const fitsbin_payload_io_ticket_t* ticket) {
    if (!ticket) {
        return 0ULL;
    }
    return ticket->sequence;
}

int fitsbin_payload_io_ticket_destroy_checked(
    fitsbin_payload_io_ticket_t* ticket) {
    if (!ticket) {
        errno = EINVAL;
        return -1;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (!fitsbin_payload_io_ticket_terminal(ticket->state) ||
        ticket->waiter_active || ticket->owner_draining) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    if (!ticket->counters_applied) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EAGAIN;
        return -1;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_ticket_free_storage(ticket);
    return 0;
}

void fitsbin_payload_io_ticket_destroy(
    fitsbin_payload_io_ticket_t* ticket) {
    int saved_errno = errno;

    if (ticket) {
        (void)fitsbin_payload_io_ticket_destroy_checked(ticket);
    }
    errno = saved_errno;
}

int fitsbin_prefetch_ranges(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget) {
    fitsbin_file_span_t spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    unsigned char scratch[FITSBIN_PREFETCH_COPY_CHUNK];
    struct timespec read_start;
    struct timespec read_finish;
    size_t actual_bytes = 0U;
    size_t i;
    size_t merged;
    unsigned long long waited;
    unsigned long long read_nanoseconds = 0ULL;
    anbool measured;
    int fd;
    int rc = 0;

    if (fitsbin_payload_io_planning_active) {
        errno = EDEADLK;
        return -1;
    }
    if (!fb || !ranges || !range_count || !byte_budget) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREFETCH_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    if (fitsbin_payload_is_fully_resident(fb)) {
        return 0;
    }
    if (fitsbin_prepare_prefetch_spans(
            fb,
            ranges,
            range_count,
            byte_budget,
            spans,
            &merged,
            &actual_bytes)) {
        return -1;
    }
    if (!merged) {
        return 0;
    }

    fd = fitsbin_payload_fd_get(fb);
    if (fd < 0) {
        return -1;
    }
    waited = fitsbin_payload_io_acquire();
    measured =
        clock_gettime(CLOCK_MONOTONIC, &read_start) == 0;
    for (i = 0U; i < merged && !rc; i++) {
        off_t cursor = spans[i].begin;

        while (cursor < spans[i].end) {
            size_t remaining =
                (size_t)(spans[i].end - cursor);
            size_t request =
                MIN(remaining, sizeof(scratch));

            if (fitsbin_pread_all(
                    fd,
                    scratch,
                    request,
                    cursor)) {
                rc = -1;
                break;
            }
            cursor += (off_t)request;
        }
    }
    if (measured &&
        clock_gettime(CLOCK_MONOTONIC, &read_finish) == 0) {
        read_nanoseconds =
            fitsbin_timespec_delta_nanoseconds(
                &read_finish,
                &read_start);
    }
    fitsbin_payload_io_release();

    __atomic_add_fetch(
        &fb->payload_wait_nanoseconds,
        waited,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_warm_nanoseconds,
        read_nanoseconds,
        __ATOMIC_RELAXED);
    if (rc) {
        __atomic_add_fetch(
            &fb->payload_failures,
            1ULL,
            __ATOMIC_RELAXED);
        return -1;
    }
    __atomic_add_fetch(
        &fb->payload_warm_calls,
        1ULL,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_warm_ranges,
        (unsigned long long)merged,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_warm_bytes,
        (unsigned long long)actual_bytes,
        __ATOMIC_RELAXED);
    return (int)merged;
}
