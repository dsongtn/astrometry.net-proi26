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
#include "fitsioutils.h"
#include "ioutils.h"
#include "fitsfile.h"
#include "errors.h"
#include "an-endian.h"
#include "tic.h"
#include "log.h"

// For in-memory: storage of previously-written extensions.
struct fitsext {
    qfits_header* header;
    char* tablename;
    bl* items;
};
typedef struct fitsext fitsext_t;

qfits_header* fitsbin_get_header(const fitsbin_t* fb, int ext) {
    assert(fb->fits);
    return anqfits_get_header(fb->fits, ext);
}

int fitsbin_get_datinfo(fitsbin_t* fb, int ext, off_t* pstart, off_t* psize) {
    assert(fb->fits);
    if (pstart)
        *pstart = anqfits_data_start(fb->fits, ext);
    if (psize)
        *psize = anqfits_data_size(fb->fits, ext);
    return 0;
}

const qfits_table* fitsbin_get_table_const(fitsbin_t* fb, int ext) {
    assert(fb->fits);
    return anqfits_get_table_const(fb->fits, ext);
}

int fitsbin_n_ext(const fitsbin_t* fb) {
    assert(fb->fits);
    return anqfits_n_ext(fb->fits);
}

FILE* fitsbin_get_fid(fitsbin_t* fb) {
    return fb->fid;
}

int fitsbin_get_open_file_stat(
    const fitsbin_t* fb,
    struct stat* file_stat) {
    if (!fb || !file_stat ||
        !fb->open_file_stat_valid) {
        errno = ENOENT;
        return -1;
    }
    *file_stat = fb->open_file_stat;
    return 0;
}

void fitsbin_stat_times(
    const struct stat* file_stat,
    time_t* mtime_seconds,
    long* mtime_nanoseconds,
    time_t* ctime_seconds,
    long* ctime_nanoseconds) {
    if (!file_stat || !mtime_seconds ||
        !mtime_nanoseconds || !ctime_seconds ||
        !ctime_nanoseconds) {
        return;
    }
    *mtime_seconds = file_stat->st_mtime;
    *ctime_seconds = file_stat->st_ctime;
#if defined(__APPLE__)
    *mtime_nanoseconds =
        file_stat->st_mtimespec.tv_nsec;
    *ctime_nanoseconds =
        file_stat->st_ctimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
    *mtime_nanoseconds = file_stat->st_mtim.tv_nsec;
    *ctime_nanoseconds = file_stat->st_ctim.tv_nsec;
#else
    *mtime_nanoseconds = 0L;
    *ctime_nanoseconds = 0L;
#endif
}

static int nchunks(fitsbin_t* fb) {
    return bl_size(fb->chunks);
}

static fitsbin_chunk_t* get_chunk(fitsbin_t* fb, int i) {
    if (i >= bl_size(fb->chunks)) {
        ERROR("Attempt to get chunk %i from a fitsbin with only %zu chunks",
              i, bl_size(fb->chunks));
        return NULL;
    }
    if (i < 0) {
        ERROR("Attempt to get fitsbin chunk %i", i);
        return NULL;
    }
    return bl_access_const(fb->chunks, i);
}

static fitsbin_t* new_fitsbin(const char* fn) {
    fitsbin_t* fb;
    fb = calloc(1, sizeof(fitsbin_t));
    if (!fb)
        return NULL;
    fb->payload_fd = -1;
    fb->payload_fd_initialized = TRUE;
    fb->chunks = bl_new(4, sizeof(fitsbin_chunk_t));
    if (!fn)
        // Can't make it NULL or qfits freaks out.
        fb->filename = strdup("");
    else
        fb->filename = strdup(fn);
    return fb;
}
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define ASTROMETRY_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#define ASTROMETRY_THREAD_LOCAL __thread
#else
#error "A thread-local storage implementation is required"
#endif

/*
 * A negative value means that this thread has no pass-level override.
 * Mapping code then preserves the original serial NORMAL behavior. The
 * configured production policy is consumed only by an explicit shard pass.
 */
static ASTROMETRY_THREAD_LOCAL int fitsbin_thread_mmap_advice = -1;
static ASTROMETRY_THREAD_LOCAL anbool
    fitsbin_thread_payload_fully_resident = FALSE;

static fitsbin_mmap_advice_t fitsbin_mmap_policy_initial_advice(
    fitsbin_mmap_policy_t policy) {
    switch (policy) {
    case FITSBIN_MMAP_POLICY_FIXED_NORMAL:
        return FITSBIN_MMAP_ADVICE_NORMAL;

    case FITSBIN_MMAP_POLICY_FIXED_RANDOM:
    case FITSBIN_MMAP_POLICY_ADAPTIVE:
#ifdef MADV_RANDOM
        return FITSBIN_MMAP_ADVICE_RANDOM;
#else
        return FITSBIN_MMAP_ADVICE_NORMAL;
#endif
    }

    return FITSBIN_MMAP_ADVICE_NORMAL;
}

fitsbin_mmap_policy_t fitsbin_mmap_policy_parse(
    const char* value) {
    if (!value || !value[0]) {
        return FITSBIN_MMAP_POLICY_FIXED_NORMAL;
    }

    if (!strcasecmp(value, "normal")) {
        return FITSBIN_MMAP_POLICY_FIXED_NORMAL;
    }

    if (!strcasecmp(value, "random")) {
        return FITSBIN_MMAP_POLICY_FIXED_RANDOM;
    }

    if (!strcasecmp(value, "adaptive")) {
        return FITSBIN_MMAP_POLICY_ADAPTIVE;
    }

    return FITSBIN_MMAP_POLICY_FIXED_NORMAL;
}

fitsbin_mmap_policy_t fitsbin_get_configured_mmap_policy(void) {
    /*
     * W2+ shard and preparation threads install RANDOM for every mapped index
     * chunk. Serial callers have no thread-local advice and retain NORMAL.
     */
    return FITSBIN_MMAP_POLICY_FIXED_RANDOM;
}

const char* fitsbin_mmap_policy_name(
    fitsbin_mmap_policy_t policy) {
    switch (policy) {
    case FITSBIN_MMAP_POLICY_FIXED_NORMAL:
        return "normal";

    case FITSBIN_MMAP_POLICY_FIXED_RANDOM:
        return "random";

    case FITSBIN_MMAP_POLICY_ADAPTIVE:
        return "adaptive";
    }

    return "unknown";
}

const char* fitsbin_mmap_advice_name(
    fitsbin_mmap_advice_t advice) {
    switch (advice) {
    case FITSBIN_MMAP_ADVICE_NORMAL:
        return "normal";

    case FITSBIN_MMAP_ADVICE_RANDOM:
        return "random";
    }

    return "unknown";
}

const char* fitsbin_mmap_region_name(
    fitsbin_mmap_region_t region) {
    switch (region) {
    case FITSBIN_MMAP_REGION_PAYLOAD:
        return "payload";

    case FITSBIN_MMAP_REGION_TOPOLOGY:
        return "topology";
    }

    return "unknown";
}

fitsbin_mmap_advice_t fitsbin_get_mmap_advice(
    const fitsbin_t* fb) {
    if (!fb) {
        return FITSBIN_MMAP_ADVICE_NORMAL;
    }

    return fb->mmap_advice;
}

fitsbin_mmap_advice_t fitsbin_get_chunk_mmap_advice(
    const fitsbin_t* fb,
    const fitsbin_chunk_t* chunk) {
    if (!chunk) {
        return FITSBIN_MMAP_ADVICE_NORMAL;
    }

    /*
     * Compact tree topology has useful traversal locality. Sparse payload
     * follows the selected index policy; bounded mapped-page population does
     * not change that stable demand fallback.
     */
    if (chunk->mmap_region == FITSBIN_MMAP_REGION_TOPOLOGY) {
        return FITSBIN_MMAP_ADVICE_NORMAL;
    }
    return fitsbin_get_mmap_advice(fb);
}

void fitsbin_mmap_advice_state_reset(
    fitsbin_mmap_advice_state_t* state) {
    if (!state) {
        return;
    }

    state->effective_advice =
        fitsbin_mmap_policy_initial_advice(state->policy);

    state->pass_number = 0;
    state->completed_clean_unsolved_passes = 0;
    state->transition_count = 0;
}

void fitsbin_mmap_advice_state_init(
    fitsbin_mmap_advice_state_t* state,
    fitsbin_mmap_policy_t policy) {
    if (!state) {
        return;
    }

    state->policy = policy;
    fitsbin_mmap_advice_state_reset(state);
}

fitsbin_mmap_advice_t fitsbin_mmap_advice_state_begin_pass(
    const fitsbin_mmap_advice_state_t* state) {
    if (!state) {
        return fitsbin_mmap_policy_initial_advice(
            fitsbin_get_configured_mmap_policy());
    }

    return state->effective_advice;
}

void fitsbin_mmap_set_thread_advice(
    fitsbin_mmap_advice_t advice) {
    fitsbin_thread_mmap_advice = (int)advice;
}

void fitsbin_mmap_clear_thread_advice(void) {
    fitsbin_thread_mmap_advice = -1;
}

fitsbin_mmap_advice_t fitsbin_mmap_current_advice(void) {
    if (fitsbin_thread_mmap_advice >= 0) {
        return (fitsbin_mmap_advice_t)
            fitsbin_thread_mmap_advice;
    }

    return FITSBIN_MMAP_ADVICE_NORMAL;
}

anbool fitsbin_mmap_thread_advice_active(void) {
    return fitsbin_thread_mmap_advice >= 0;
}

void fitsbin_payload_set_thread_full_resident(void) {
    fitsbin_thread_payload_fully_resident = TRUE;
}

void fitsbin_payload_clear_thread_full_resident(void) {
    fitsbin_thread_payload_fully_resident = FALSE;
}

anbool fitsbin_payload_is_fully_resident(
    const fitsbin_t* fb) {
    if (!fb) {
        return FALSE;
    }
    return fb->payload_fully_resident;
}

anbool fitsbin_mmap_policy_complete_pass(
    fitsbin_mmap_advice_state_t* state,
    anbool pass_completed,
    anbool pass_exhaustive,
    anbool pass_solved,
    anbool pass_cancelled,
    int pass_rc,
    int pass_status) {
    if (!state) {
        return FALSE;
    }

    /*
     * pass_number counts fully returned pass attempts, including solved or
     * unsuccessful attempts. Partial execution is not counted.
     */
    if (pass_completed) {
        state->pass_number++;
    }

    /*
     * Only a complete, exhaustive, clean, unsolved pass is evidence for
     * changing the following pass.
     */
    if (!pass_completed ||
        !pass_exhaustive ||
        pass_solved ||
        pass_cancelled ||
        pass_rc != 0 ||
        pass_status != 0) {
        return FALSE;
    }

    if (state->policy != FITSBIN_MMAP_POLICY_ADAPTIVE) {
        return FALSE;
    }

    state->completed_clean_unsolved_passes++;

    if (state->effective_advice ==
        FITSBIN_MMAP_ADVICE_NORMAL) {
        return FALSE;
    }

    state->effective_advice =
        FITSBIN_MMAP_ADVICE_NORMAL;

    state->transition_count++;

    return TRUE;
}

static int fitsbin_mmap_os_advice(fitsbin_mmap_advice_t advice) {

    switch (advice) {
    case FITSBIN_MMAP_ADVICE_NORMAL:
        return MADV_NORMAL;

    case FITSBIN_MMAP_ADVICE_RANDOM:
#ifdef MADV_RANDOM
        return MADV_RANDOM;
#else
        return MADV_NORMAL;
#endif
    }



    return MADV_NORMAL;
}

int fitsbin_configure_index_mmap(fitsbin_t* fb) {
    long page_size;

    if (!fb) {
        return -1;
    }

    page_size = sysconf(_SC_PAGESIZE);
    if (page_size > 0) {
        fb->mmap_page_size = (size_t)page_size;
    } else {
        fb->mmap_page_size = 0;
    }

    /*
     * A shard worker installs its pass policy before opening an index. Serial
     * callers have no thread-local override and retain NORMAL.
     */
    fb->mmap_advice = fitsbin_mmap_current_advice();
    fb->payload_fully_resident =
        fitsbin_thread_payload_fully_resident;
    fb->mmap_advice_failed = FALSE;

    /*
     * Legacy whole-file prefetch and buffered solver warming remain
     * disconnected. A successful bounded preparation must mean that all
     * planned pages have been populated, not merely queued for readahead.
     */
#if defined(MADV_POPULATE_READ)
    fb->mmap_prefetch_enabled = TRUE;
#else
    fb->mmap_prefetch_enabled = FALSE;
#endif
    fb->mmap_prefetch_failed = FALSE;

    return 0;
}

#define FITSBIN_PREFETCH_RANGE_LIMIT \
    FITSBIN_PREAD_ASYNC_RANGE_LIMIT
#define FITSBIN_PREFETCH_COPY_CHUNK (16U * 1024U)
#define FITSBIN_PAYLOAD_POPULATE_TOTAL_BUDGET \
    (16U * 1024U * 1024U)
#define FITSBIN_PAYLOAD_IO_MAX_LANES 2
#define FITSBIN_PAYLOAD_IO_MAX_JOBS 24U
#define FITSBIN_PAYLOAD_IO_MAX_BYTES \
    (64U * 1024U * 1024U)
#define FITSBIN_PAYLOAD_IO_PRIORITY_COUNT 3U
#define FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_JOBS 1U
#define FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_BYTES \
    (FITSBIN_PAYLOAD_IO_MAX_BYTES / 4U)

typedef struct fitsbin_file_span {
    off_t begin;
    off_t end;
} fitsbin_file_span_t;

typedef struct fitsbin_mapped_span {
    uintptr_t map_begin;
    uintptr_t map_end;
    uintptr_t begin;
    uintptr_t end;
} fitsbin_mapped_span_t;

typedef struct fitsbin_prepared_pread_range {
    off_t offset;
    size_t size;
    size_t logical_size;
    void* destination;
} fitsbin_prepared_pread_range_t;

static int fitsbin_compare_prepared_pread_range(
    const void* left,
    const void* right) {
    const fitsbin_prepared_pread_range_t* lhs = left;
    const fitsbin_prepared_pread_range_t* rhs = right;

    if (lhs->offset < rhs->offset) {
        return -1;
    }
    if (lhs->offset > rhs->offset) {
        return 1;
    }
    if (lhs->size < rhs->size) {
        return -1;
    }
    if (lhs->size > rhs->size) {
        return 1;
    }
    return 0;
}

static anbool fitsbin_prepared_pread_destinations_disjoint(
    const fitsbin_prepared_pread_range_t* ranges,
    size_t range_count) {
    size_t i;

    if (!ranges) {
        return FALSE;
    }
    for (i = 0U; i < range_count; i++) {
        uintptr_t begin = (uintptr_t)ranges[i].destination;
        uintptr_t end;
        size_t j;

        if (ranges[i].size > UINTPTR_MAX - begin) {
            return FALSE;
        }
        end = begin + ranges[i].size;
        for (j = 0U; j < i; j++) {
            uintptr_t other_begin =
                (uintptr_t)ranges[j].destination;
            uintptr_t other_end;

            if (ranges[j].size >
                UINTPTR_MAX - other_begin) {
                return FALSE;
            }
            other_end = other_begin + ranges[j].size;
            if (begin < other_end && other_begin < end) {
                return FALSE;
            }
        }
    }
    return TRUE;
}

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
    fitsbin_mapped_span_t spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    fitsbin_prepared_pread_range_t* ranges;
    fitsbin_t* source;
    size_t span_count;
    size_t range_count;
    size_t byte_count;
    size_t logical_byte_count;
    unsigned long long page_count;
    unsigned long long sequence;
    unsigned long long read_nanoseconds;
    fitsbin_payload_io_ticket_state_t state;
    fitsbin_payload_io_ticket_kind_t kind;
    fitsbin_payload_io_priority_t priority;
    struct fitsbin_payload_io_ticket* next;
    struct fitsbin_payload_io_ticket* wait_next;
    pthread_cond_t completion_cv;
    int fd;
    int saved_errno;
    anbool cancel_requested;
    anbool counters_applied;
    anbool wait_registered;
    anbool helper_waiter;
};

static pthread_mutex_t fitsbin_payload_io_mutex =
    PTHREAD_MUTEX_INITIALIZER;
/* Loader queue publication and service lifecycle changes. */
static pthread_cond_t fitsbin_payload_io_cv =
    PTHREAD_COND_INITIALIZER;
/* Reader-credit availability, kept separate from ticket waiters. */
static pthread_cond_t fitsbin_payload_io_credit_cv =
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
static unsigned long long fitsbin_payload_io_service_cancelled = 0ULL;
static unsigned long long fitsbin_payload_io_service_failed = 0ULL;
static unsigned long long fitsbin_payload_io_service_planned_bytes = 0ULL;
static unsigned long long fitsbin_payload_io_service_read_nanoseconds = 0ULL;
static int fitsbin_payload_io_service_lanes = 0;
static anbool fitsbin_payload_io_service_running = FALSE;
static anbool fitsbin_payload_io_service_accepting = FALSE;
static anbool fitsbin_payload_io_service_stopping = FALSE;
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

static unsigned long long fitsbin_timespec_delta_nanoseconds(
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

static unsigned long long fitsbin_payload_io_acquire(void) {
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

static void fitsbin_payload_io_release(void) {
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

static anbool fitsbin_file_identity_matches(
    const struct stat* expected,
    const struct stat* actual) {
    time_t expected_mtime_seconds;
    long expected_mtime_nanoseconds;
    time_t expected_ctime_seconds;
    long expected_ctime_nanoseconds;
    time_t actual_mtime_seconds;
    long actual_mtime_nanoseconds;
    time_t actual_ctime_seconds;
    long actual_ctime_nanoseconds;

    if (!expected || !actual) {
        return FALSE;
    }
    fitsbin_stat_times(
        expected,
        &expected_mtime_seconds,
        &expected_mtime_nanoseconds,
        &expected_ctime_seconds,
        &expected_ctime_nanoseconds);
    fitsbin_stat_times(
        actual,
        &actual_mtime_seconds,
        &actual_mtime_nanoseconds,
        &actual_ctime_seconds,
        &actual_ctime_nanoseconds);
    return expected->st_dev == actual->st_dev &&
        expected->st_ino == actual->st_ino &&
        expected->st_size == actual->st_size &&
        expected_mtime_seconds == actual_mtime_seconds &&
        expected_mtime_nanoseconds ==
            actual_mtime_nanoseconds &&
        expected_ctime_seconds == actual_ctime_seconds &&
        expected_ctime_nanoseconds ==
            actual_ctime_nanoseconds;
}

static int fitsbin_payload_fd_get(fitsbin_t* fb) {
    int fd;
    int opened;
    int expected;
    int advice_rc;
    int saved_errno;
    struct stat actual;

    if (!fb || !fb->payload_fd_initialized ||
        !fb->filename ||
        !fb->open_file_stat_valid) {
        errno = ENOTSUP;
        return -1;
    }
    fd = __atomic_load_n(
        &fb->payload_fd,
        __ATOMIC_ACQUIRE);
    if (fd >= 0) {
        return fd;
    }
    if (__atomic_load_n(
            &fb->payload_fd_failed,
            __ATOMIC_ACQUIRE)) {
        errno = ENOTSUP;
        return -1;
    }

#ifdef O_CLOEXEC
    opened = open(fb->filename, O_RDONLY | O_CLOEXEC);
#else
    opened = open(fb->filename, O_RDONLY);
#endif
    if (opened < 0) {
        goto fail;
    }
    if (fstat(opened, &actual) ||
        !fitsbin_file_identity_matches(
            &fb->open_file_stat,
            &actual)) {
        close(opened);
        errno = ESTALE;
        goto fail;
    }
#if defined(POSIX_FADV_RANDOM)
    advice_rc = posix_fadvise(
        opened,
        0,
        0,
        POSIX_FADV_RANDOM);
    if (advice_rc) {
        close(opened);
        errno = advice_rc;
        goto fail;
    }
#else
    close(opened);
    errno = ENOTSUP;
    goto fail;
#endif

    expected = -1;
    if (!__atomic_compare_exchange_n(
            &fb->payload_fd,
            &expected,
            opened,
            FALSE,
            __ATOMIC_RELEASE,
            __ATOMIC_ACQUIRE)) {
        close(opened);
        if (expected >= 0) {
            return expected;
        }
        errno = EIO;
        goto fail;
    }
    return opened;

fail:
    saved_errno = errno;
    fd = __atomic_load_n(
        &fb->payload_fd,
        __ATOMIC_ACQUIRE);
    if (fd >= 0) {
        return fd;
    }
    __atomic_store_n(
        &fb->payload_fd_failed,
        TRUE,
        __ATOMIC_RELEASE);
    __atomic_add_fetch(
        &fb->payload_failures,
        1ULL,
        __ATOMIC_RELAXED);
    errno = saved_errno;
    return -1;
}

static fitsbin_chunk_t* fitsbin_find_data_chunk(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    size_t* chunk_offset) {
    uintptr_t request;
    int i;

    if (!fb || !data || !size || !chunk_offset ||
        !fb->chunks) {
        errno = EINVAL;
        return NULL;
    }
    request = (uintptr_t)data;
    for (i = 0; i < nchunks(fb); i++) {
        fitsbin_chunk_t* chunk = get_chunk(fb, i);
        uintptr_t begin;
        uintptr_t end;

        if (!chunk || !chunk->data ||
            !chunk->data_file_size) {
            continue;
        }
        begin = (uintptr_t)chunk->data;
        if (chunk->data_file_size > UINTPTR_MAX - begin) {
            continue;
        }
        end = begin + chunk->data_file_size;
        if (request < begin || request >= end ||
            size > end - request) {
            continue;
        }
        *chunk_offset = (size_t)(request - begin);
        return chunk;
    }
    errno = ERANGE;
    return NULL;
}

static int fitsbin_pread_all(
    int fd,
    void* destination,
    size_t size,
    off_t offset) {
    unsigned char* output = destination;
    size_t complete = 0U;

    if (fd < 0 || !destination || !size ||
        offset < 0 ||
        size > (size_t)LLONG_MAX ||
        offset > (off_t)LLONG_MAX - (off_t)size) {
        errno = EINVAL;
        return -1;
    }
    while (complete < size) {
        ssize_t count = pread(
            fd,
            output + complete,
            size - complete,
            offset + (off_t)complete);

        if (count > 0) {
            complete += (size_t)count;
            continue;
        }
        if (count < 0 && errno == EINTR) {
            continue;
        }
        if (!count) {
            errno = EIO;
        }
        return -1;
    }
    return 0;
}

int fitsbin_mapped_range_page_cover(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    const void** cover_data,
    size_t* cover_size,
    off_t* cover_file_offset,
    size_t* exact_offset) {
    fitsbin_chunk_t* chunk;
    size_t chunk_offset;
    size_t page_size;
    unsigned long long chunk_begin;
    unsigned long long chunk_end;
    unsigned long long request_begin;
    unsigned long long request_end;
    unsigned long long cover_begin;
    unsigned long long cover_end;
    unsigned long long remainder;
    size_t cover_chunk_offset;
    uintptr_t cover_address;

    if (!fb || !data || !size || !cover_data ||
        !cover_size || !cover_file_offset || !exact_offset) {
        errno = EINVAL;
        return -1;
    }
    *cover_data = NULL;
    *cover_size = 0U;
    *cover_file_offset = 0;
    *exact_offset = 0U;
    chunk = fitsbin_find_data_chunk(
        fb,
        data,
        size,
        &chunk_offset);
    if (!chunk || chunk->data_file_offset < 0) {
        return -1;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }
    chunk_begin =
        (unsigned long long)chunk->data_file_offset;
    if ((unsigned long long)chunk->data_file_size >
            ULLONG_MAX - chunk_begin ||
        (unsigned long long)chunk_offset >
            ULLONG_MAX - chunk_begin) {
        errno = EOVERFLOW;
        return -1;
    }
    chunk_end = chunk_begin +
        (unsigned long long)chunk->data_file_size;
    request_begin = chunk_begin +
        (unsigned long long)chunk_offset;
    if ((unsigned long long)size >
        ULLONG_MAX - request_begin) {
        errno = EOVERFLOW;
        return -1;
    }
    request_end = request_begin + (unsigned long long)size;
    cover_begin = request_begin -
        request_begin % (unsigned long long)page_size;
    if (cover_begin < chunk_begin) {
        cover_begin = chunk_begin;
    }
    remainder = request_end % (unsigned long long)page_size;
    cover_end = request_end;
    if (remainder) {
        unsigned long long padding =
            (unsigned long long)page_size - remainder;

        if (padding > ULLONG_MAX - cover_end) {
            cover_end = chunk_end;
        } else {
            cover_end += padding;
        }
    }
    if (cover_end > chunk_end) {
        cover_end = chunk_end;
    }
    if (cover_end <= cover_begin ||
        cover_begin > (unsigned long long)LLONG_MAX ||
        cover_end - cover_begin > SIZE_MAX ||
        request_begin - cover_begin > SIZE_MAX ||
        cover_begin - chunk_begin > SIZE_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    cover_chunk_offset =
        (size_t)(cover_begin - chunk_begin);
    cover_address = (uintptr_t)chunk->data;
    if (cover_chunk_offset > UINTPTR_MAX - cover_address) {
        errno = EOVERFLOW;
        return -1;
    }
    if (request_begin - cover_begin > SIZE_MAX ||
        (size_t)(request_begin - cover_begin) >
            (size_t)(cover_end - cover_begin) ||
        size > (size_t)(cover_end - cover_begin) -
            (size_t)(request_begin - cover_begin)) {
        errno = EOVERFLOW;
        return -1;
    }
    *cover_data = (const void*)(cover_address + cover_chunk_offset);
    *cover_size = (size_t)(cover_end - cover_begin);
    *cover_file_offset = (off_t)cover_begin;
    *exact_offset = (size_t)(request_begin - cover_begin);
    return 0;
}

int fitsbin_pread_mapped_ranges(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count) {
    fitsbin_prepared_pread_range_t
        prepared[FITSBIN_PREAD_RANGE_LIMIT];
    struct timespec read_start;
    struct timespec read_finish;
    unsigned long long waited;
    unsigned long long read_nanoseconds = 0ULL;
    unsigned long long physical_bytes = 0ULL;
    unsigned long long logical_bytes = 0ULL;
    unsigned long long page_count = 0ULL;
    size_t page_size;
    anbool measured;
    int saved_errno = 0;
    int fd;
    int rc = 0;
    size_t i;

    if (!fb || !ranges || !range_count) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREAD_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }
    for (i = 0U; i < range_count; i++) {
        fitsbin_chunk_t* chunk;
        size_t chunk_offset;
        unsigned long long file_offset;
        unsigned long long last_offset;
        unsigned long long first_page;
        unsigned long long last_page;
        unsigned long long page_span;

        if (!ranges[i].data || !ranges[i].size ||
            !ranges[i].logical_size ||
            ranges[i].logical_size > ranges[i].size ||
            !ranges[i].destination) {
            errno = EINVAL;
            return -1;
        }
        chunk = fitsbin_find_data_chunk(
            fb,
            ranges[i].data,
            ranges[i].size,
            &chunk_offset);
        if (!chunk || chunk->data_file_offset < 0 ||
            chunk_offset > (size_t)LLONG_MAX ||
            chunk->data_file_offset >
                (off_t)LLONG_MAX - (off_t)chunk_offset) {
            if (chunk) {
                errno = EOVERFLOW;
            }
            return -1;
        }
        prepared[i].offset =
            chunk->data_file_offset + (off_t)chunk_offset;
        prepared[i].size = ranges[i].size;
        prepared[i].logical_size = ranges[i].logical_size;
        prepared[i].destination = ranges[i].destination;
        if (ranges[i].size > (size_t)LLONG_MAX ||
            prepared[i].offset >
                (off_t)LLONG_MAX -
                    (off_t)ranges[i].size) {
            errno = EOVERFLOW;
            return -1;
        }
        if ((unsigned long long)ranges[i].size >
                ULLONG_MAX - physical_bytes ||
            (unsigned long long)ranges[i].logical_size >
                ULLONG_MAX - logical_bytes) {
            errno = EOVERFLOW;
            return -1;
        }
        physical_bytes += (unsigned long long)ranges[i].size;
        logical_bytes +=
            (unsigned long long)ranges[i].logical_size;
        file_offset =
            (unsigned long long)prepared[i].offset;
        if ((unsigned long long)ranges[i].size - 1ULL >
            ULLONG_MAX - file_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        last_offset = file_offset +
            (unsigned long long)ranges[i].size - 1ULL;
        first_page = file_offset /
            (unsigned long long)page_size;
        last_page = last_offset /
            (unsigned long long)page_size;
        page_span = last_page - first_page;
        if (page_span == ULLONG_MAX ||
            page_span + 1ULL >
                ULLONG_MAX - page_count) {
            errno = EOVERFLOW;
            return -1;
        }
        page_count += page_span + 1ULL;
    }
    if (range_count > 1U &&
        fitsbin_prepared_pread_destinations_disjoint(
            prepared, range_count)) {
        qsort(
            prepared,
            range_count,
            sizeof(prepared[0]),
            fitsbin_compare_prepared_pread_range);
    }
    fd = fitsbin_payload_fd_get(fb);
    if (fd < 0) {
        return -1;
    }
    waited = fitsbin_payload_io_acquire();
    measured =
        clock_gettime(CLOCK_MONOTONIC, &read_start) == 0;
    for (i = 0U; i < range_count; i++) {
        if (fitsbin_pread_all(
                fd,
                prepared[i].destination,
                prepared[i].size,
                prepared[i].offset)) {
            saved_errno = errno;
            rc = -1;
            break;
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
        &fb->payload_read_nanoseconds,
        read_nanoseconds,
        __ATOMIC_RELAXED);
    if (rc) {
        __atomic_add_fetch(
            &fb->payload_failures,
            1ULL,
            __ATOMIC_RELAXED);
        errno = saved_errno;
        return -1;
    }
    __atomic_add_fetch(
        &fb->payload_read_batches,
        1ULL,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_read_calls,
        (unsigned long long)range_count,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_read_bytes,
        physical_bytes,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_read_logical_bytes,
        logical_bytes,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_read_pages,
        page_count,
        __ATOMIC_RELAXED);
    return 0;
}

int fitsbin_pread_mapped_range(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    void* destination) {
    fitsbin_pread_range_t range;

    range.data = data;
    range.size = size;
    range.logical_size = size;
    range.destination = destination;
    return fitsbin_pread_mapped_ranges(fb, &range, 1U);
}

static int fitsbin_prepare_direct_ranges(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_prepared_pread_range_t* prepared,
    size_t* physical_bytes_out,
    size_t* logical_bytes_out,
    unsigned long long* page_count_out) {
    size_t physical_bytes = 0U;
    size_t logical_bytes = 0U;
    unsigned long long page_count = 0ULL;
    size_t page_size;
    size_t i;

    if (!fb || !ranges || !range_count || !byte_budget ||
        !prepared || !physical_bytes_out ||
        !logical_bytes_out || !page_count_out) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREAD_ASYNC_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }
    for (i = 0U; i < range_count; i++) {
        fitsbin_chunk_t* chunk;
        size_t chunk_offset;
        unsigned long long file_offset;
        unsigned long long last_offset;
        unsigned long long first_page;
        unsigned long long last_page;
        unsigned long long page_span;

        if (!ranges[i].data || !ranges[i].size ||
            !ranges[i].logical_size ||
            ranges[i].logical_size > ranges[i].size ||
            !ranges[i].destination) {
            errno = EINVAL;
            return -1;
        }
        chunk = fitsbin_find_data_chunk(
            fb,
            ranges[i].data,
            ranges[i].size,
            &chunk_offset);
        if (!chunk || chunk->data_file_offset < 0 ||
            chunk_offset > (size_t)LLONG_MAX ||
            chunk->data_file_offset >
                (off_t)LLONG_MAX - (off_t)chunk_offset) {
            if (chunk) {
                errno = EOVERFLOW;
            }
            return -1;
        }
        prepared[i].offset =
            chunk->data_file_offset + (off_t)chunk_offset;
        prepared[i].size = ranges[i].size;
        prepared[i].logical_size = ranges[i].logical_size;
        prepared[i].destination = ranges[i].destination;
        if (ranges[i].size > (size_t)LLONG_MAX ||
            prepared[i].offset >
                (off_t)LLONG_MAX -
                    (off_t)ranges[i].size) {
            errno = EOVERFLOW;
            return -1;
        }
        if (ranges[i].size > SIZE_MAX - physical_bytes ||
            ranges[i].logical_size >
                SIZE_MAX - logical_bytes) {
            errno = EOVERFLOW;
            return -1;
        }
        physical_bytes += ranges[i].size;
        logical_bytes += ranges[i].logical_size;
        if (physical_bytes > byte_budget) {
            errno = E2BIG;
            return -1;
        }
        file_offset =
            (unsigned long long)prepared[i].offset;
        if ((unsigned long long)ranges[i].size - 1ULL >
            ULLONG_MAX - file_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        last_offset = file_offset +
            (unsigned long long)ranges[i].size - 1ULL;
        first_page = file_offset /
            (unsigned long long)page_size;
        last_page = last_offset /
            (unsigned long long)page_size;
        page_span = last_page - first_page;
        if (page_span == ULLONG_MAX ||
            page_span + 1ULL >
                ULLONG_MAX - page_count) {
            errno = EOVERFLOW;
            return -1;
        }
        page_count += page_span + 1ULL;
    }
    *physical_bytes_out = physical_bytes;
    *logical_bytes_out = logical_bytes;
    *page_count_out = page_count;
    if (range_count > 1U &&
        fitsbin_prepared_pread_destinations_disjoint(
            prepared, range_count)) {
        qsort(
            prepared,
            range_count,
            sizeof(prepared[0]),
            fitsbin_compare_prepared_pread_range);
    }
    return 0;
}

static int fitsbin_compare_file_span(
    const void* left,
    const void* right) {
    const fitsbin_file_span_t* lhs = left;
    const fitsbin_file_span_t* rhs = right;

    if (lhs->begin < rhs->begin) {
        return -1;
    }
    if (lhs->begin > rhs->begin) {
        return 1;
    }
    if (lhs->end < rhs->end) {
        return -1;
    }
    if (lhs->end > rhs->end) {
        return 1;
    }
    return 0;
}

static int fitsbin_compare_mapped_span(
    const void* left,
    const void* right) {
    const fitsbin_mapped_span_t* lhs = left;
    const fitsbin_mapped_span_t* rhs = right;

    if (lhs->map_begin < rhs->map_begin) {
        return -1;
    }
    if (lhs->map_begin > rhs->map_begin) {
        return 1;
    }
    if (lhs->begin < rhs->begin) {
        return -1;
    }
    if (lhs->begin > rhs->begin) {
        return 1;
    }
    if (lhs->end < rhs->end) {
        return -1;
    }
    if (lhs->end > rhs->end) {
        return 1;
    }
    return 0;
}

static int fitsbin_prepare_prefetch_spans(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_file_span_t* spans,
    size_t* span_count,
    size_t* byte_count) {
    struct stat source;
    size_t accepted = 0U;
    size_t actual_bytes = 0U;
    size_t page_size;
    size_t i;
    size_t merged;

    if (!span_count || !byte_count) {
        errno = EINVAL;
        return -1;
    }
    *span_count = 0U;
    *byte_count = 0U;
    if (!range_count) {
        return 0;
    }
    if (!fb || !ranges || !byte_budget || !spans) {
        errno = EINVAL;
        return -1;
    }
    if (range_count > FITSBIN_PREFETCH_RANGE_LIMIT) {
        errno = E2BIG;
        return -1;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }
    if (fitsbin_get_open_file_stat(fb, &source)) {
        return -1;
    }

    for (i = 0U; i < range_count; i++) {
        fitsbin_chunk_t* chunk;
        size_t chunk_offset;
        off_t begin;
        off_t end;
        off_t aligned_begin;
        off_t aligned_end;

        if (!ranges[i].data || !ranges[i].size) {
            errno = EINVAL;
            return -1;
        }
        chunk = fitsbin_find_data_chunk(
            fb,
            ranges[i].data,
            ranges[i].size,
            &chunk_offset);
        if (!chunk) {
            errno = EINVAL;
            return -1;
        }
        if (chunk_offset > (size_t)LLONG_MAX ||
            chunk->data_file_offset >
                (off_t)LLONG_MAX - (off_t)chunk_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        begin = chunk->data_file_offset +
            (off_t)chunk_offset;
        if (ranges[i].size > (size_t)LLONG_MAX ||
            begin >
                (off_t)LLONG_MAX -
                    (off_t)ranges[i].size) {
            errno = EOVERFLOW;
            return -1;
        }
        end = begin + (off_t)ranges[i].size;
        aligned_begin =
            begin - begin % (off_t)page_size;
        aligned_end = end;
        if (aligned_end % (off_t)page_size) {
            off_t padding =
                (off_t)page_size -
                aligned_end % (off_t)page_size;

            if (aligned_end > source.st_size - padding) {
                aligned_end = source.st_size;
            } else {
                aligned_end += padding;
            }
        }
        if (aligned_end > source.st_size) {
            aligned_end = source.st_size;
        }
        if (aligned_end <= aligned_begin) {
            errno = EINVAL;
            return -1;
        }
        spans[accepted].begin = aligned_begin;
        spans[accepted].end = aligned_end;
        accepted++;
    }

    qsort(
        spans,
        accepted,
        sizeof(spans[0]),
        fitsbin_compare_file_span);
    merged = 0U;
    for (i = 0U; i < accepted; i++) {
        if (merged &&
            spans[i].begin <= spans[merged - 1U].end) {
            if (spans[i].end > spans[merged - 1U].end) {
                spans[merged - 1U].end = spans[i].end;
            }
            continue;
        }
        spans[merged++] = spans[i];
    }
    for (i = 0U; i < merged; i++) {
        size_t span_bytes =
            (size_t)(spans[i].end - spans[i].begin);

        if (actual_bytes > byte_budget ||
            span_bytes > byte_budget - actual_bytes) {
            errno = E2BIG;
            return -1;
        }
        actual_bytes += span_bytes;
    }
    *span_count = merged;
    *byte_count = actual_bytes;
    return 0;
}

static int fitsbin_prepare_mapped_spans(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_mapped_span_t* spans,
    size_t span_capacity,
    size_t* span_count,
    size_t* byte_count,
    size_t* logical_byte_count,
    unsigned long long* page_count) {
    size_t accepted = 0U;
    size_t merged = 0U;
    size_t aligned_bytes = 0U;
    size_t logical_bytes = 0U;
    unsigned long long pages = 0ULL;
    size_t page_size;
    size_t i;

    if (!span_count || !byte_count || !logical_byte_count ||
        !page_count) {
        errno = EINVAL;
        return -1;
    }
    *span_count = 0U;
    *byte_count = 0U;
    *logical_byte_count = 0U;
    *page_count = 0ULL;
    if (!fb || !ranges || !range_count || !byte_budget ||
        !spans || range_count > span_capacity) {
        errno = range_count > span_capacity ? E2BIG : EINVAL;
        return -1;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }

    for (i = 0U; i < range_count; i++) {
        const void* map_base;
        const void* range_start;
        size_t map_size;
        size_t range_size;
        uintptr_t map_begin;
        uintptr_t map_end;
        uintptr_t begin;
        uintptr_t end;
        uintptr_t remainder;
        int resolved;

        if (!ranges[i].data || !ranges[i].size ||
            ranges[i].size > SIZE_MAX - logical_bytes) {
            errno = ranges[i].data && ranges[i].size
                ? EOVERFLOW
                : EINVAL;
            return -1;
        }
        resolved = fitsbin_resolve_mapped_range(
            fb,
            ranges[i].data,
            ranges[i].size,
            &map_base,
            &map_size,
            &range_start,
            &range_size);
        if (resolved != 1 || range_size != ranges[i].size) {
            if (!resolved) {
                errno = ERANGE;
            }
            return -1;
        }
        map_begin = (uintptr_t)map_base;
        if (map_size > UINTPTR_MAX - map_begin) {
            errno = EOVERFLOW;
            return -1;
        }
        map_end = map_begin + map_size;
        begin = (uintptr_t)range_start;
        if (range_size > UINTPTR_MAX - begin) {
            errno = EOVERFLOW;
            return -1;
        }
        end = begin + range_size;
        begin -= begin % (uintptr_t)page_size;
        if (begin < map_begin) {
            begin = map_begin;
        }
        remainder = end % (uintptr_t)page_size;
        if (remainder) {
            uintptr_t padding =
                (uintptr_t)page_size - remainder;

            end = padding > map_end - end
                ? map_end
                : end + padding;
        }
        if (end > map_end) {
            end = map_end;
        }
        if (end <= begin) {
            errno = ERANGE;
            return -1;
        }
        spans[accepted].map_begin = map_begin;
        spans[accepted].map_end = map_end;
        spans[accepted].begin = begin;
        spans[accepted].end = end;
        accepted++;
        logical_bytes += ranges[i].size;
    }

    qsort(
        spans,
        accepted,
        sizeof(spans[0]),
        fitsbin_compare_mapped_span);
    for (i = 0U; i < accepted; i++) {
        if (merged &&
            spans[i].map_begin == spans[merged - 1U].map_begin &&
            spans[i].map_end == spans[merged - 1U].map_end &&
            spans[i].begin <= spans[merged - 1U].end) {
            if (spans[i].end > spans[merged - 1U].end) {
                spans[merged - 1U].end = spans[i].end;
            }
            continue;
        }
        spans[merged++] = spans[i];
    }
    for (i = 0U; i < merged; i++) {
        size_t span_bytes =
            (size_t)(spans[i].end - spans[i].begin);
        size_t span_pages = span_bytes / page_size;

        if (span_bytes % page_size) {
            span_pages++;
        }
        if (aligned_bytes > byte_budget ||
            span_bytes > byte_budget - aligned_bytes ||
            (unsigned long long)span_pages > ULLONG_MAX - pages) {
            errno = E2BIG;
            return -1;
        }
        aligned_bytes += span_bytes;
        pages += (unsigned long long)span_pages;
    }
    *span_count = merged;
    *byte_count = aligned_bytes;
    *logical_byte_count = logical_bytes;
    *page_count = pages;
    return 0;
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
static void fitsbin_payload_io_cancel_queued_locked(
    fitsbin_payload_io_ticket_t* ticket) {
    if (!ticket ||
        ticket->state != FITSBIN_PAYLOAD_IO_SUBMITTED ||
        !__atomic_load_n(
            &ticket->cancel_requested,
            __ATOMIC_ACQUIRE) ||
        !fitsbin_payload_io_remove_queued_locked(ticket)) {
        return;
    }
    assert(fitsbin_payload_io_service_jobs > 0U);
    assert(fitsbin_payload_io_service_bytes >=
           ticket->byte_count);
    fitsbin_payload_io_service_jobs--;
    fitsbin_payload_io_service_bytes -= ticket->byte_count;
    if (ticket->fd >= 0) {
        close(ticket->fd);
        ticket->fd = -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_CANCELLED;
    fitsbin_payload_io_service_cancelled++;
    pthread_cond_signal(&ticket->completion_cv);
}

/* fitsbin_payload_io_mutex must be held. */
static anbool fitsbin_payload_io_admit_locked(
    const fitsbin_payload_io_ticket_t* ticket) {
    size_t max_jobs = FITSBIN_PAYLOAD_IO_MAX_JOBS;
    size_t max_bytes = FITSBIN_PAYLOAD_IO_MAX_BYTES;

    if (!ticket ||
        !fitsbin_payload_io_priority_valid(ticket->priority)) {
        return FALSE;
    }
    if (ticket->priority != FITSBIN_PAYLOAD_IO_PRIORITY_DEMAND) {
        max_jobs -= FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_JOBS;
        max_bytes -= FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_BYTES;
    }
    if (fitsbin_payload_io_service_jobs >= max_jobs ||
        ticket->byte_count > max_bytes ||
        fitsbin_payload_io_service_bytes >
            max_bytes - ticket->byte_count) {
        return FALSE;
    }
    return TRUE;
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

static void* fitsbin_payload_io_service_worker(void* opaque) {
    (void)opaque;
    while (1) {
        fitsbin_payload_io_ticket_t* ticket;
        struct timespec read_start;
        struct timespec read_finish;
        unsigned long long read_nanoseconds = 0ULL;
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

                if (fitsbin_pread_all(
                        ticket->fd,
                        range->destination,
                        range->size,
                        range->offset)) {
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
            for (work_index = 0U;
                 !status && work_index < ticket->span_count &&
                     !cancelled;
                 work_index++) {
                fitsbin_mapped_span_t* span =
                    &ticket->spans[work_index];
                size_t span_bytes =
                    (size_t)(span->end - span->begin);

                if (!span_bytes ||
                    madvise(
                        (void*)span->begin,
                        span_bytes,
                        MADV_POPULATE_READ)) {
                    saved_errno = errno ? errno : EIO;
                    status = -1;
                    if (ticket->source) {
                        __atomic_store_n(
                            &ticket->source->mmap_prefetch_failed,
                            TRUE,
                            __ATOMIC_RELEASE);
                    }
                    break;
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
               ticket->byte_count);
        fitsbin_payload_io_service_jobs--;
        fitsbin_payload_io_service_bytes -= ticket->byte_count;
        ticket->read_nanoseconds = read_nanoseconds;
        fitsbin_payload_io_service_read_nanoseconds += read_nanoseconds;
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
            }
        }
        pthread_cond_signal(&ticket->completion_cv);
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
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
    fitsbin_payload_io_service_running = TRUE;
    fitsbin_payload_io_service_accepting = TRUE;
    fitsbin_payload_io_service_stopping = FALSE;
    fitsbin_payload_io_service_lanes = lane_count;
    fitsbin_payload_io_service_submitted = 0ULL;
    fitsbin_payload_io_service_ready = 0ULL;
    fitsbin_payload_io_service_direct_submitted = 0ULL;
    fitsbin_payload_io_service_direct_ready = 0ULL;
    fitsbin_payload_io_service_cancelled = 0ULL;
    fitsbin_payload_io_service_failed = 0ULL;
    fitsbin_payload_io_service_planned_bytes = 0ULL;
    fitsbin_payload_io_service_read_nanoseconds = 0ULL;
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
    fitsbin_payload_io_service_stopping = FALSE;
    fitsbin_payload_io_service_lanes = 0;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    errno = status;
    return -1;
}

void fitsbin_payload_io_service_stop(void) {
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
            "cancelled=%llu failed=%llu planned_bytes=%llu "
            "read_wall=%.6f\n",
            fitsbin_payload_io_service_submitted,
            fitsbin_payload_io_service_ready,
            fitsbin_payload_io_service_direct_submitted,
            fitsbin_payload_io_service_direct_ready,
            fitsbin_payload_io_service_cancelled,
            fitsbin_payload_io_service_failed,
            fitsbin_payload_io_service_planned_bytes,
            (double)fitsbin_payload_io_service_read_nanoseconds /
                1000000000.0);
    fitsbin_payload_io_service_running = FALSE;
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

static int fitsbin_payload_io_submit_ticket(
    fitsbin_payload_io_ticket_t* ticket,
    fitsbin_payload_io_ticket_t** ticket_out) {
    anbool admitted;

    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    admitted =
        fitsbin_payload_io_service_running &&
        fitsbin_payload_io_service_accepting &&
        fitsbin_payload_io_admit_locked(ticket);
    if (!admitted) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        fitsbin_payload_io_ticket_free_storage(ticket);
        return 0;
    }
    ticket->sequence = ++fitsbin_payload_io_next_sequence;
    ticket->state = FITSBIN_PAYLOAD_IO_SUBMITTED;
    fitsbin_payload_io_enqueue_locked(ticket);
    fitsbin_payload_io_service_jobs++;
    fitsbin_payload_io_service_bytes += ticket->byte_count;
    fitsbin_payload_io_service_submitted++;
    if (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
        fitsbin_payload_io_service_direct_submitted++;
    }
    fitsbin_payload_io_service_planned_bytes +=
        (unsigned long long)ticket->byte_count;
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
    if (!fitsbin_payload_io_service_running ||
        !fitsbin_payload_io_service_accepting) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return 0;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    ticket = fitsbin_payload_io_ticket_alloc();
    if (!ticket) {
        return -1;
    }
    ticket->state = FITSBIN_PAYLOAD_IO_PLANNED;
    ticket->kind = FITSBIN_PAYLOAD_IO_TICKET_DIRECT;
    ticket->priority = priority;
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
        return 0;
    }
    if (!fb->mmap_prefetch_enabled ||
        __atomic_load_n(
            &fb->mmap_prefetch_failed,
            __ATOMIC_ACQUIRE)) {
        return 0;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (!fitsbin_payload_io_service_running ||
        !fitsbin_payload_io_service_accepting) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
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
    if (fitsbin_prepare_mapped_spans(
            fb,
            ranges,
            range_count,
            byte_budget,
            ticket->spans,
            FITSBIN_PREFETCH_RANGE_LIMIT,
            &ticket->span_count,
            &ticket->byte_count,
            &ticket->logical_byte_count,
            &ticket->page_count)) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return -1;
    }
    ticket->range_count = range_count;
    if (!ticket->span_count) {
        fitsbin_payload_io_ticket_free_storage(ticket);
        return 0;
    }
    return fitsbin_payload_io_submit_ticket(ticket, ticket_out);
}

static int fitsbin_payload_io_ticket_wait_internal(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket,
    anbool cancel) {
    struct timespec wait_start;
    struct timespec wait_finish;
    unsigned long long waited = 0ULL;
    unsigned long long read_nanoseconds;
    size_t span_count;
    size_t range_count;
    size_t byte_count;
    size_t logical_byte_count;
    unsigned long long page_count;
    fitsbin_payload_io_ticket_kind_t kind;
    anbool helper_registered = FALSE;
    anbool waiter_registered = FALSE;
    anbool retry_helper = TRUE;
    anbool measured;
    fitsbin_payload_io_ticket_state_t state;
    int saved_errno;

    if (!fb || !ticket ||
        (ticket->kind == FITSBIN_PAYLOAD_IO_TICKET_PREFETCH &&
         ticket->source != fb)) {
        errno = EINVAL;
        return -1;
    }
    measured = clock_gettime(CLOCK_MONOTONIC, &wait_start) == 0;
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    if (ticket->wait_registered) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EBUSY;
        return -1;
    }
    if (cancel) {
        __atomic_store_n(
            &ticket->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        fitsbin_payload_io_cancel_queued_locked(ticket);
    }
    if (ticket->state == FITSBIN_PAYLOAD_IO_SUBMITTED) {
        if (fitsbin_payload_io_register_waiter_locked(ticket)) {
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
            fitsbin_payload_io_cancel_queued_locked(ticket);
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
    state = ticket->state;
    saved_errno = ticket->saved_errno;
    read_nanoseconds = ticket->read_nanoseconds;
    span_count = ticket->span_count;
    range_count = ticket->range_count;
    byte_count = ticket->byte_count;
    logical_byte_count = ticket->logical_byte_count;
    page_count = ticket->page_count;
    kind = ticket->kind;
    if (ticket->counters_applied) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        errno = EALREADY;
        return -1;
    }
    ticket->counters_applied = TRUE;
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);

    if (measured &&
        clock_gettime(CLOCK_MONOTONIC, &wait_finish) == 0) {
        waited = fitsbin_timespec_delta_nanoseconds(
            &wait_finish,
            &wait_start);
    }
    __atomic_add_fetch(
        &fb->payload_wait_nanoseconds,
        waited,
        __ATOMIC_RELAXED);
    if (state == FITSBIN_PAYLOAD_IO_READY &&
        kind == FITSBIN_PAYLOAD_IO_TICKET_DIRECT) {
        __atomic_add_fetch(
            &fb->payload_read_batches,
            1ULL,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_calls,
            (unsigned long long)range_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_bytes,
            (unsigned long long)byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_logical_bytes,
            (unsigned long long)logical_byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_pages,
            page_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_read_nanoseconds,
            read_nanoseconds,
            __ATOMIC_RELAXED);
        return (int)range_count;
    }
    if (state == FITSBIN_PAYLOAD_IO_READY) {
        __atomic_add_fetch(
            &fb->payload_warm_calls,
            1ULL,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_ranges,
            (unsigned long long)span_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_bytes,
            (unsigned long long)byte_count,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_nanoseconds,
            read_nanoseconds,
            __ATOMIC_RELAXED);
        return (int)span_count;
    }
    if (state == FITSBIN_PAYLOAD_IO_FAILED) {
        __atomic_add_fetch(
            &fb->payload_failures,
            1ULL,
            __ATOMIC_RELAXED);
        errno = saved_errno ? saved_errno : EIO;
        return -1;
    }
    errno = ECANCELED;
    return 0;
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

void fitsbin_payload_io_ticket_destroy(
    fitsbin_payload_io_ticket_t* ticket) {
    anbool terminal;

    if (!ticket) {
        return;
    }
    pthread_mutex_lock(&fitsbin_payload_io_mutex);
    terminal = ticket->state == FITSBIN_PAYLOAD_IO_READY ||
        ticket->state == FITSBIN_PAYLOAD_IO_FAILED ||
        ticket->state == FITSBIN_PAYLOAD_IO_CANCELLED;
    if (!terminal || !ticket->counters_applied) {
        pthread_mutex_unlock(&fitsbin_payload_io_mutex);
        return;
    }
    pthread_mutex_unlock(&fitsbin_payload_io_mutex);
    fitsbin_payload_io_ticket_free_storage(ticket);
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

int fitsbin_close_payload_fd(fitsbin_t* fb) {
    int fd;

    if (!fb || !fb->payload_fd_initialized) {
        return 0;
    }
    fd = __atomic_exchange_n(
        &fb->payload_fd,
        -1,
        __ATOMIC_ACQ_REL);
    if (fd >= 0 && close(fd)) {
        return -1;
    }
    __atomic_store_n(
        &fb->payload_fd_failed,
        FALSE,
        __ATOMIC_RELEASE);
    return 0;
}

void fitsbin_take_payload_io_stats(
    fitsbin_t* fb,
    fitsbin_payload_io_stats_t* stats) {
    if (!stats) {
        return;
    }
    memset(stats, 0, sizeof(*stats));
    if (!fb) {
        return;
    }
    stats->read_calls = __atomic_exchange_n(
        &fb->payload_read_calls,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->read_batches = __atomic_exchange_n(
        &fb->payload_read_batches,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->read_logical_bytes = __atomic_exchange_n(
        &fb->payload_read_logical_bytes,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->read_pages = __atomic_exchange_n(
        &fb->payload_read_pages,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->read_bytes = __atomic_exchange_n(
        &fb->payload_read_bytes,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->read_nanoseconds = __atomic_exchange_n(
        &fb->payload_read_nanoseconds,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->warm_calls = __atomic_exchange_n(
        &fb->payload_warm_calls,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->warm_ranges = __atomic_exchange_n(
        &fb->payload_warm_ranges,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->warm_bytes = __atomic_exchange_n(
        &fb->payload_warm_bytes,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->warm_nanoseconds = __atomic_exchange_n(
        &fb->payload_warm_nanoseconds,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->cache_hits = __atomic_exchange_n(
        &fb->payload_cache_hits,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->cache_misses = __atomic_exchange_n(
        &fb->payload_cache_misses,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->cache_evictions = __atomic_exchange_n(
        &fb->payload_cache_evictions,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->cache_allocations = __atomic_exchange_n(
        &fb->payload_cache_allocations,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->wait_nanoseconds = __atomic_exchange_n(
        &fb->payload_wait_nanoseconds,
        0ULL,
        __ATOMIC_ACQ_REL);
    stats->failures = __atomic_exchange_n(
        &fb->payload_failures,
        0ULL,
        __ATOMIC_ACQ_REL);
}
int fitsbin_resolve_mapped_range(fitsbin_t* fb,
                                 const void* data,
                                 size_t size,
                                 const void** map_base,
                                 size_t* map_size,
                                 const void** range_start,
                                 size_t* range_size) {
    uintptr_t request_start;
    int i;

    if (!fb ||
        !data ||
        !size ||
        !map_base ||
        !map_size ||
        !range_start ||
        !range_size) {
        errno = EINVAL;
        return -1;
    }

    *map_base = NULL;
    *map_size = 0;
    *range_start = NULL;
    *range_size = 0;

    if (!fb->chunks) {
        return 0;
    }

    request_start = (uintptr_t)data;

    for (i = 0; i < nchunks(fb); i++) {
        fitsbin_chunk_t* chunk = get_chunk(fb, i);
        uintptr_t mapping_start;
        uintptr_t mapping_end;
        size_t clipped_size;

        if (!chunk ||
            !chunk->map ||
            !chunk->mapsize) {
            continue;
        }

        mapping_start = (uintptr_t)chunk->map;

        if (chunk->mapsize > UINTPTR_MAX - mapping_start) {
            continue;
        }

        mapping_end = mapping_start + chunk->mapsize;

        if (request_start < mapping_start ||
            request_start >= mapping_end) {
            continue;
        }

        clipped_size = size;

        if (clipped_size > mapping_end - request_start) {
            clipped_size = (size_t)(mapping_end - request_start);
        }

        if (!clipped_size) {
            return 0;
        }

        *map_base = chunk->map;
        *map_size = chunk->mapsize;
        *range_start = data;
        *range_size = clipped_size;

        return 1;
    }

    return 0;
}

int fitsbin_set_mmap_range_advice(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    fitsbin_mmap_advice_t advice) {
    const void* map_base;
    const void* range_start;
    size_t map_size;
    size_t range_size;
    size_t page_size;
    uintptr_t map_begin;
    uintptr_t map_end;
    uintptr_t begin;
    uintptr_t end;
    uintptr_t remainder;
    int resolved;
    int native_advice;

    if (!fb || !data || !size ||
        (advice != FITSBIN_MMAP_ADVICE_NORMAL &&
         advice != FITSBIN_MMAP_ADVICE_RANDOM)) {
        errno = EINVAL;
        return -1;
    }
    resolved = fitsbin_resolve_mapped_range(
        fb,
        data,
        size,
        &map_base,
        &map_size,
        &range_start,
        &range_size);
    if (resolved <= 0) {
        return resolved;
    }
    if (range_size != size) {
        errno = ERANGE;
        return -1;
    }

    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected_page_size = sysconf(_SC_PAGESIZE);

        if (detected_page_size <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected_page_size;
    }

    map_begin = (uintptr_t)map_base;
    if (map_size > UINTPTR_MAX - map_begin) {
        errno = EOVERFLOW;
        return -1;
    }
    map_end = map_begin + map_size;
    begin = (uintptr_t)range_start;
    if (range_size > UINTPTR_MAX - begin) {
        end = map_end;
    } else {
        end = begin + range_size;
    }

    begin -= begin % (uintptr_t)page_size;
    if (begin < map_begin) {
        begin = map_begin;
    }
    remainder = end % (uintptr_t)page_size;
    if (remainder) {
        uintptr_t padding = (uintptr_t)page_size - remainder;

        if (padding > map_end - end) {
            end = map_end;
        } else {
            end += padding;
        }
    }
    if (end > map_end) {
        end = map_end;
    }
    if (end <= begin) {
        return 0;
    }

    native_advice = fitsbin_mmap_os_advice(advice);
    if (madvise((void*)begin, (size_t)(end - begin), native_advice)) {
        return -1;
    }
    return 1;
}

static void fitsbin_restore_normal_mmap_advice(
    fitsbin_t* fb) {
    int i;

    if (!fb) {
        return;
    }
    fb->mmap_advice = FITSBIN_MMAP_ADVICE_NORMAL;
    fb->mmap_advice_failed = TRUE;
    if (!fb->chunks) {
        return;
    }
    for (i = 0; i < bl_size(fb->chunks); i++) {
        fitsbin_chunk_t* chunk = bl_access_const(fb->chunks, i);

        if (!chunk || !chunk->map || !chunk->mapsize) {
            continue;
        }
        if (madvise(chunk->map, chunk->mapsize, MADV_NORMAL)) {
            logmsg("Warning: madvise(MADV_NORMAL) fallback failed for %s "
                   "table %s region %s: %s\n",
                   fb->filename ? fb->filename : "(unknown)",
                   chunk->tablename ? chunk->tablename : "(unknown)",
                   fitsbin_mmap_region_name(chunk->mmap_region),
                   strerror(errno));
        }
    }
}

static int fitsbin_mapped_population_failure(
    fitsbin_t* fb,
    int saved_errno) {
    if (!saved_errno) {
        saved_errno = EIO;
    }
    fb->mmap_prefetch_failed = TRUE;
    __atomic_add_fetch(
        &fb->payload_failures,
        1ULL,
        __ATOMIC_RELAXED);
    errno = saved_errno;
    return -1;
}

int fitsbin_advise_mapped_ranges(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget) {
#if defined(MADV_POPULATE_READ)
    fitsbin_mapped_span_t stack_spans[FITSBIN_PREFETCH_RANGE_LIMIT];
    fitsbin_mapped_span_t* spans = stack_spans;
    size_t accepted = 0U;
    size_t advised = 0U;
    size_t advised_bytes = 0U;
    size_t plan_byte_limit;
    size_t merged;
    size_t page_size;
    size_t i;
    unsigned long long waited = 0U;
    unsigned long long populate_nanoseconds = 0U;
    struct timespec populate_start;
    struct timespec populate_finish;
    anbool measured = FALSE;

    if (!fb) {
        errno = EINVAL;
        return -1;
    }
    if ((!ranges && range_count) ||
        (range_count && !byte_budget)) {
        return fitsbin_mapped_population_failure(fb, EINVAL);
    }
    if (!range_count || !byte_budget ||
        fb->mmap_prefetch_failed) {
        return 0;
    }
    if (range_count > (size_t)INT_MAX) {
        return fitsbin_mapped_population_failure(fb, EOVERFLOW);
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, EINVAL);
        }
        page_size = (size_t)detected;
    }
    if (range_count > FITSBIN_PREFETCH_RANGE_LIMIT) {
        if (range_count > SIZE_MAX / sizeof(*spans)) {
            return fitsbin_mapped_population_failure(fb, EOVERFLOW);
        }
        spans = malloc(range_count * sizeof(*spans));
        if (!spans) {
            return fitsbin_mapped_population_failure(fb, ENOMEM);
        }
    }

    for (i = 0U; i < range_count; i++) {
        const void* map_base;
        const void* range_start;
        size_t map_size;
        size_t range_size;
        uintptr_t map_begin;
        uintptr_t map_end;
        uintptr_t begin;
        uintptr_t end;
        uintptr_t remainder;
        int resolved;

        if (!ranges[i].data || !ranges[i].size) {
            continue;
        }
        resolved = fitsbin_resolve_mapped_range(
            fb,
            ranges[i].data,
            ranges[i].size,
            &map_base,
            &map_size,
            &range_start,
            &range_size);
        if (resolved < 0) {
            int saved_errno = errno;

            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, saved_errno);
        }
        if (!resolved) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, ERANGE);
        }
        if (range_size != ranges[i].size) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, ERANGE);
        }

        map_begin = (uintptr_t)map_base;
        if (map_size > UINTPTR_MAX - map_begin) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, EOVERFLOW);
        }
        map_end = map_begin + map_size;
        begin = (uintptr_t)range_start;
        if (range_size > UINTPTR_MAX - begin) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(fb, EOVERFLOW);
        }
        end = begin + range_size;

        begin -= begin % (uintptr_t)page_size;
        if (begin < map_begin) {
            begin = map_begin;
        }
        remainder = end % (uintptr_t)page_size;
        if (remainder) {
            uintptr_t padding = (uintptr_t)page_size - remainder;

            if (padding > map_end - end) {
                end = map_end;
            } else {
                end += padding;
            }
        }
        if (end > map_end) {
            end = map_end;
        }
        if (end <= begin) {
            continue;
        }

        spans[accepted].map_begin = map_begin;
        spans[accepted].map_end = map_end;
        spans[accepted].begin = begin;
        spans[accepted].end = end;
        accepted++;
    }
    if (!accepted) {
        if (spans != stack_spans) {
            free(spans);
        }
        return 0;
    }

    qsort(spans,
          accepted,
          sizeof(spans[0]),
          fitsbin_compare_mapped_span);
    merged = 0U;
    for (i = 0U; i < accepted; i++) {
        if (merged &&
            spans[i].map_begin == spans[merged - 1U].map_begin &&
            spans[i].map_end == spans[merged - 1U].map_end &&
            spans[i].begin <= spans[merged - 1U].end) {
            if (spans[i].end > spans[merged - 1U].end) {
                spans[merged - 1U].end = spans[i].end;
            }
            continue;
        }
        spans[merged++] = spans[i];
    }

    for (i = 0U; i < merged; i++) {
        size_t span_bytes =
            (size_t)(spans[i].end - spans[i].begin);

        if (span_bytes > SIZE_MAX - advised_bytes) {
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(
                fb, EOVERFLOW);
        }
        advised_bytes += span_bytes;
    }
    {
        int capacity = __atomic_load_n(
            &fitsbin_payload_io_capacity,
            __ATOMIC_ACQUIRE);

        if (capacity < 1) {
            capacity = 1;
        }
        plan_byte_limit =
            FITSBIN_PAYLOAD_POPULATE_TOTAL_BUDGET /
            (size_t)capacity;
    }
    if (plan_byte_limit < page_size) {
        if (spans != stack_spans) {
            free(spans);
        }
        return 0;
    }
    byte_budget = MIN(byte_budget, plan_byte_limit);
    if (advised_bytes > byte_budget ||
        fitsbin_payload_io_demand_busy()) {
        if (spans != stack_spans) {
            free(spans);
        }
        return 0;
    }

    waited = fitsbin_payload_io_acquire();
    measured =
        clock_gettime(CLOCK_MONOTONIC,
                      &populate_start) == 0;
    for (i = 0U; i < merged; i++) {
        size_t span_bytes =
            (size_t)(spans[i].end - spans[i].begin);

        if (madvise((void*)spans[i].begin,
                    span_bytes,
                    MADV_POPULATE_READ)) {
            int saved_errno = errno;

            if (measured &&
                clock_gettime(
                    CLOCK_MONOTONIC,
                    &populate_finish) == 0) {
                populate_nanoseconds =
                    fitsbin_timespec_delta_nanoseconds(
                        &populate_finish,
                        &populate_start);
            }
            fitsbin_payload_io_release();
            __atomic_add_fetch(
                &fb->payload_wait_nanoseconds,
                waited,
                __ATOMIC_RELAXED);
            __atomic_add_fetch(
                &fb->payload_warm_nanoseconds,
                populate_nanoseconds,
                __ATOMIC_RELAXED);
            if (spans != stack_spans) {
                free(spans);
            }
            return fitsbin_mapped_population_failure(
                fb, saved_errno);
        }
        advised++;
    }
    if (measured &&
        clock_gettime(
            CLOCK_MONOTONIC,
            &populate_finish) == 0) {
        populate_nanoseconds =
            fitsbin_timespec_delta_nanoseconds(
                &populate_finish,
                &populate_start);
    }
    fitsbin_payload_io_release();
    __atomic_add_fetch(
        &fb->payload_wait_nanoseconds,
        waited,
        __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_warm_nanoseconds,
        populate_nanoseconds,
        __ATOMIC_RELAXED);

    if (advised) {
        __atomic_add_fetch(
            &fb->payload_warm_calls,
            1ULL,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_ranges,
            (unsigned long long)advised,
            __ATOMIC_RELAXED);
        __atomic_add_fetch(
            &fb->payload_warm_bytes,
            (unsigned long long)advised_bytes,
            __ATOMIC_RELAXED);
    }
    if (spans != stack_spans) {
        free(spans);
    }
    return (int)advised;
#else
    if (!fb) {
        errno = EINVAL;
        return -1;
    }
    (void)ranges;
    (void)range_count;
    (void)byte_budget;
    return 0;
#endif
}

int fitsbin_set_mmap_advice(
    fitsbin_t* fb,
    fitsbin_mmap_advice_t advice,
    anbool reapply_existing) {
    int first_error = 0;
    int i;

    if (!fb) {
        errno = EINVAL;
        return -1;
    }

    if (advice != FITSBIN_MMAP_ADVICE_NORMAL &&
        advice != FITSBIN_MMAP_ADVICE_RANDOM) {
        errno = EINVAL;
        return -1;
    }

    /*
     * Avoid repeatedly issuing madvise() when several code paths encounter
     * the same fitsbin object during one pass.
     */
    if (fb->mmap_advice == advice &&
        (!reapply_existing ||
         !fb->mmap_advice_failed)) {
        return 0;
    }

    fb->mmap_advice = advice;

    if (!reapply_existing) {
        if (!fb->chunks ||
            bl_size(fb->chunks) == 0) {
            fb->mmap_advice_failed = FALSE;
        }
        return 0;
    }

    if (!fb->chunks ||
        bl_size(fb->chunks) == 0) {
        fb->mmap_advice_failed = FALSE;
        return 0;
    }

    /*
     * Reapply the selected index policy to every existing mapped chunk.
     * Use the same chunk traversal expression used by fitsbin_close().
     */
    for (i = 0; i < bl_size(fb->chunks); i++) {
        fitsbin_chunk_t* chunk = bl_access_const(fb->chunks, i);
        fitsbin_mmap_advice_t chunk_advice;
        int native_advice;

        if (!chunk || !chunk->map || !chunk->mapsize) {
            continue;
        }

        chunk_advice = fitsbin_get_chunk_mmap_advice(
            fb, chunk);
        native_advice =
            fitsbin_mmap_os_advice(chunk_advice);

        if (madvise(chunk->map, chunk->mapsize, native_advice) != 0) {
            int saved_errno = errno;

            if (!first_error) {
                first_error = saved_errno;
            }

            logmsg("Warning: madvise(%s) failed for %s table %s "
                   "region %s: %s\n",
                   chunk_advice == FITSBIN_MMAP_ADVICE_RANDOM
                       ? "MADV_RANDOM"
                       : "MADV_NORMAL",
                   fb->filename ? fb->filename : "(unknown)",
                   chunk->tablename ? chunk->tablename : "(unknown)",
                   fitsbin_mmap_region_name(chunk->mmap_region),
                   strerror(saved_errno));
        }
    }

    if (first_error) {
        fitsbin_restore_normal_mmap_advice(fb);
        errno = first_error;
        return -1;
    }

    fb->mmap_advice_failed = FALSE;
    return 0;
}
int fitsbin_prefetch_data(fitsbin_t* fb, const void* data, size_t size) {
#ifdef MADV_WILLNEED
    uintptr_t request_start;
    uintptr_t request_end;
    size_t page_size;
    int i;

    if (!fb || !data || !size ||
        fb->mmap_advice != FITSBIN_MMAP_ADVICE_RANDOM ||
        !fb->mmap_prefetch_enabled || fb->mmap_prefetch_failed) {
        return 0;
    }

    request_start = (uintptr_t)data;
    page_size = fb->mmap_page_size;

    if (!page_size) {
        long detected_page_size = sysconf(_SC_PAGESIZE);

        if (detected_page_size <= 0) {
            fb->mmap_prefetch_failed = TRUE;
            return -1;
        }

        page_size = (size_t)detected_page_size;
    }

    for (i = 0; i < nchunks(fb); i++) {
        fitsbin_chunk_t* chunk = get_chunk(fb, i);
        uintptr_t map_start;
        uintptr_t map_end;
        uintptr_t advise_start;
        uintptr_t advise_end;
        uintptr_t remainder;

        if (!chunk || !chunk->map || !chunk->mapsize) {
            continue;
        }

        map_start = (uintptr_t)chunk->map;
        if (chunk->mapsize > UINTPTR_MAX - map_start) {
            continue;
        }
        map_end = map_start + chunk->mapsize;

        if (request_start < map_start || request_start >= map_end) {
            continue;
        }

        if (size > map_end - request_start) {
            request_end = map_end;
        } else {
            request_end = request_start + size;
        }

        advise_start =
            request_start - request_start % (uintptr_t)page_size;
        advise_end = request_end;
        remainder = advise_end % (uintptr_t)page_size;

        if (remainder) {
            uintptr_t padding = (uintptr_t)page_size - remainder;

            if (padding > map_end - advise_end) {
                advise_end = map_end;
            } else {
                advise_end += padding;
            }
        }

        if (advise_end <= advise_start) {
            return 0;
        }

        if (madvise((void*)advise_start,
                    (size_t)(advise_end - advise_start),
                    MADV_WILLNEED)) {
            const char* filename =
                fb->filename ? fb->filename : "(unknown file)";
            const char* tablename =
                chunk->tablename ? chunk->tablename : "(unknown table)";

            logmsg("Warning: madvise(MADV_WILLNEED) failed for %s "
                   "table %s: %s; disabling mmap prefetch for this file.\n",
                   filename, tablename, strerror(errno));
            fb->mmap_prefetch_failed = TRUE;
            return -1;
        }

        return 1;
    }
#else
    (void)fb;
    (void)data;
    (void)size;
#endif

    return 0;
}

// Apply the file mapping advice before the first table access.
static void apply_mmap_advice(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    fitsbin_mmap_advice_t advice;
    const char* filename;
    const char* tablename;

    advice = fitsbin_get_chunk_mmap_advice(fb, chunk);

    if (advice == FITSBIN_MMAP_ADVICE_NORMAL ||
        fb->mmap_advice_failed) {
        return;
    }

    filename = fb->filename != NULL
        ? fb->filename
        : "(unknown file)";

    tablename = chunk->tablename != NULL
        ? chunk->tablename
        : "(unknown table)";

#ifdef MADV_RANDOM
    if (advice == FITSBIN_MMAP_ADVICE_RANDOM) {
        if (madvise(chunk->map, chunk->mapsize, MADV_RANDOM) != 0) {
            logmsg("Warning: madvise(MADV_RANDOM) failed for %s table %s "
                   "region %s: %s; using normal mmap advice for this "
                   "file.\n",
                   filename,
                   tablename,
                   fitsbin_mmap_region_name(chunk->mmap_region),
                   strerror(errno));

            fitsbin_restore_normal_mmap_advice(fb);
            return;
        }

        debug("Applied MADV_RANDOM to %zu bytes for %s table %s "
              "region %s.\n",
              chunk->mapsize,
              filename,
              tablename,
              fitsbin_mmap_region_name(chunk->mmap_region));
    }
#else
    (void)filename;
    (void)tablename;
#endif
}
static anbool in_memory(fitsbin_t* fb) {
    return fb->inmemory;
}


static void free_chunk(fitsbin_chunk_t* chunk) {
    if (!chunk) return;
    free(chunk->tablename_copy);
    if (chunk->header)
        qfits_header_destroy(chunk->header);
    if (chunk->map) {
        if (munmap(chunk->map, chunk->mapsize)) {
            SYSERROR("Failed to munmap fitsbin chunk");
        }
    }
}

void fitsbin_chunk_init(fitsbin_chunk_t* chunk) {
    memset(chunk, 0, sizeof(fitsbin_chunk_t));
}

void fitsbin_chunk_clean(fitsbin_chunk_t* chunk) {
    free_chunk(chunk);
}

void fitsbin_chunk_reset(fitsbin_chunk_t* chunk) {
    fitsbin_chunk_clean(chunk);
    fitsbin_chunk_init(chunk);
}

fitsbin_chunk_t* fitsbin_get_chunk(fitsbin_t* fb, int chunk) {
    return get_chunk(fb, chunk);
}

int fitsbin_n_chunks(fitsbin_t* fb) {
    return nchunks(fb);
}

fitsbin_chunk_t* fitsbin_add_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    chunk = bl_append(fb->chunks, chunk);
    chunk->tablename_copy = strdup(chunk->tablename);
    chunk->tablename = chunk->tablename_copy;
    return chunk;
}

off_t fitsbin_get_data_start(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    return chunk->header_end;
}

int fitsbin_close_fd(fitsbin_t* fb) {
    if (!fb) return 0;
    if (fb->fid) {
        FILE* fid = fb->fid;

        /*
         * Clear first: fclose() invalidates the stream even when it reports
         * a delayed write/close error. Never leave a dangling FILE* behind.
         */
        fb->fid = NULL;
        if (fclose(fid)) {
            SYSERROR("Error closing fitsbin file");
            return -1;
        }
    }
    return 0;
}

int fitsbin_close(fitsbin_t* fb) {
    int i;
    int rtn = 0;
    if (!fb) return rtn;
    if (fitsbin_close_payload_fd(fb)) {
        rtn = -1;
    }
    if (fitsbin_close_fd(fb)) {
        rtn = -1;
    }
    if (fb->primheader)
        qfits_header_destroy(fb->primheader);
    for (i=0; i<nchunks(fb); i++) {
        if (in_memory(fb)) {
            free(get_chunk(fb, i)->data);
        }
        free_chunk(get_chunk(fb, i));
    }
    free(fb->filename);
    if (fb->chunks)
        bl_free(fb->chunks);

    if (in_memory(fb)) {
        for (i=0; i<bl_size(fb->extensions); i++) {
            fitsext_t* ext = bl_access(fb->extensions, i);
            bl_free(ext->items);
            qfits_header_destroy(ext->header);
            free(ext->tablename);
        }
        bl_free(fb->extensions);
        bl_free(fb->items);
    }

    if (fb->tables) {
        for (i=0; i<fb->Next; i++) {
            if (!fb->tables[i])
                continue;
            qfits_table_close(fb->tables[i]);
        }
        free(fb->tables);
    }

    if (fb->owns_fits && fb->fits) {
        anqfits_close(fb->fits);
        fb->fits = NULL;
    }

    free(fb);
    return rtn;
}

int fitsbin_write_primary_header(fitsbin_t* fb) {
    if (in_memory(fb)) return 0;
    return fitsfile_write_primary_header(fb->fid, fb->primheader,
                                         &fb->primheader_end, fb->filename);
}

int fitsbin_write_primary_header_to(fitsbin_t* fb, FILE* fid) {
    off_t end;
    return fitsfile_write_primary_header(fid, fb->primheader, &end, "");
}

qfits_header* fitsbin_get_primary_header(const fitsbin_t* fb) {
    return fb->primheader;
}

void fitsbin_set_primary_header(fitsbin_t* fb, const qfits_header* hdr) {
    qfits_header_destroy(fb->primheader);
    fb->primheader = qfits_header_copy(hdr);
}

int fitsbin_fix_primary_header(fitsbin_t* fb) {
    if (in_memory(fb)) return 0;
    return fitsfile_fix_primary_header(fb->fid, fb->primheader,
                                       &fb->primheader_end, fb->filename);
}

qfits_header* fitsbin_get_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    qfits_table* table;
    int tablesize;
    qfits_header* hdr;
    int ncols = 1;
    char* fn = NULL;

    if (chunk->header)
        return chunk->header;

    // Create the new header.

    if (fb)
        fn = fb->filename;
    if (!fn)
        fn = "";
    // the table header
    tablesize = chunk->itemsize * chunk->nrows * ncols;
    table = qfits_table_new(fn, QFITS_BINTABLE, tablesize, ncols, chunk->nrows);
    assert(table);
    qfits_col_fill(table->col, chunk->itemsize, 0, 1,
                   chunk->forced_type ? chunk->forced_type : TFITS_BIN_TYPE_A,
                   chunk->tablename, "", "", "", 0, 0, 0, 0, 0);
    hdr = qfits_table_ext_header_default(table);
    qfits_table_close(table);
    chunk->header = hdr;
    return hdr;
}

static int write_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk, int flipped) {
    int N;
    if (fitsbin_write_chunk_header(fb, chunk)) {
        return -1;
    }
    N = chunk->nrows;
    if (!flipped) {
        if (fitsbin_write_items(fb, chunk, chunk->data, chunk->nrows))
            return -1;
    } else {
        // endian-flip words of the data of length "flipped", write them,
        // then flip them back to the way they were.

        // NO, copy to temp array, flip it, write it.

        // this is slow, but it won't be run very often...

        int i, j;
        int nper = chunk->itemsize / flipped;
        char tempdata[chunk->itemsize];
        assert(chunk->itemsize >= flipped);
        assert(nper * flipped == chunk->itemsize);
        for (i=0; i<N; i++) {
            // copy it...
            memcpy(tempdata, chunk->data + i*chunk->itemsize, chunk->itemsize);
            // swap it...
            for (j=0; j<nper; j++)
                endian_swap(tempdata + j*flipped, flipped);
            // write it...
            fitsbin_write_item(fb, chunk, tempdata);
        }
    }
    chunk->nrows -= N;
    if (fitsbin_fix_chunk_header(fb, chunk)) {
        return -1;
    }
    return 0;
}

int fitsbin_write_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    return write_chunk(fb, chunk, 0);
}

int fitsbin_write_chunk_to(fitsbin_t* fb, fitsbin_chunk_t* chunk, FILE* fid) {
    //logmsg("fitsbin_write_chunk_to: table %s, itemsize %i, nrows %i, header_start %lu, header_end %lu\n", chunk->tablename_copy, chunk->itemsize, chunk->nrows, chunk->header_start, chunk->header_end);
    //off_t off = ftello(fid);
    //logmsg("offset: %lu\n", off);
    if (fitsbin_write_chunk_header_to(fb, chunk, fid) ||
        fitsbin_write_items_to(chunk, chunk->data, chunk->nrows, fid))
        return -1;
    return 0;
}

int fitsbin_write_chunk_flipped(fitsbin_t* fb, fitsbin_chunk_t* chunk,
                                int wordsize) {
    return write_chunk(fb, chunk, wordsize);
}

int fitsbin_write_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    qfits_header* hdr;
    hdr = fitsbin_get_chunk_header(fb, chunk);
    if (in_memory(fb)) return 0;
    if (fitsfile_write_header(fb->fid, hdr,
                              &chunk->header_start, &chunk->header_end,
                              -1, fb->filename)) {
        return -1;
    }
    return 0;
}

int fitsbin_write_chunk_header_to(fitsbin_t* fb, fitsbin_chunk_t* chunk, FILE* fid) {
    off_t start, end;
    qfits_header* hdr;
    hdr = fitsbin_get_chunk_header(fb, chunk);
    if (fitsfile_write_header(fid, hdr, &start, &end, -1, ""))
        return -1;
    return 0;
}

int fitsbin_fix_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    // update NAXIS2 to reflect the number of rows written.
    fits_header_mod_int(chunk->header, "NAXIS2", chunk->nrows, NULL);

    // HACK -- leverage the fact that this is the last function called for each chunk...
    if (in_memory(fb)) {
        // Save this chunk.
        fitsext_t ext;
        // table, header, items
        if (!fb->extensions)
            fb->extensions = bl_new(4, sizeof(fitsext_t));
        ext.header = qfits_header_copy(chunk->header);
        ext.items = fb->items;
        ext.tablename = strdup(chunk->tablename);
        bl_append(fb->extensions, &ext);
        fb->items = NULL;
        return 0;
    }

    if (fitsfile_fix_header(fb->fid, chunk->header,
                            &chunk->header_start, &chunk->header_end,
                            -1, fb->filename)) {
        return -1;
    }
    return 0;
}

int fitsbin_write_items_to(fitsbin_chunk_t* chunk, void* data, int N, FILE* fid) {
    off_t offset;
    if (fwrite(data, chunk->itemsize, N, fid) != N) {
        SYSERROR("Failed to write %i items", N);
        return -1;
    }
    offset = ftello(fid);
    fits_pad_file(fid);
    if (fseeko(fid, offset, SEEK_SET)) {
        SYSERROR("Failed to fseeko in fitsbin_write_items_to.");
        return -1;
    }
    return 0;
}

int fitsbin_write_items(fitsbin_t* fb, fitsbin_chunk_t* chunk, void* data, int N) {
    if (in_memory(fb)) {
        int i;
        char* src = data;
        if (!fb->items)
            fb->items = bl_new(1024, chunk->itemsize);
        for (i=0; i<N; i++) {
            bl_append(fb->items, src);
            src += chunk->itemsize;
        }
    } else {
        if (fitsbin_write_items_to(chunk, data, N, fb->fid))
            return -1;
    }
    chunk->nrows += N;
    return 0;
}

int fitsbin_write_item(fitsbin_t* fb, fitsbin_chunk_t* chunk, void* data) {
    return fitsbin_write_items(fb, chunk, data, 1);
}

// Like fitsioutils.c : fits_find_table_column(), but using our cache...
static int find_table_column(fitsbin_t* fb, const char* colname, off_t* pstart, off_t* psize, int* pext) {
    int i;
    for (i=1; i<fb->Next; i++) {
        int c;
        const qfits_table* table = fitsbin_get_table_const(fb, i);
        if (!table)
            continue;
        c = fits_find_column(table, colname);
        if (c == -1)
            continue;
        if (fitsbin_get_datinfo(fb, i, pstart, psize)) {
            ERROR("error getting start/size for ext %i in file %s.\n", i, fb->filename);
            return -1;
        }
        if (pext) *pext = i;
        return 0;
    }
    debug("searched %i extensions in file %s but didn't find a table with a column \"%s\".\n",
          fb->Next, fb->filename, colname);
    return -1;
}

static int read_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk,
                      anbool read_payload) {
    off_t tabstart=0, tabsize=0;
    int ext;
    size_t expected = 0;
    int mode, flags;
    off_t mapstart;
    int mapoffset;
    int table_nrows;
    int table_rowsize;
    fitsext_t* inmemext = NULL;

    if (in_memory(fb)) {
        int i;
        anbool gotit = FALSE;
        for (i=0; i<bl_size(fb->extensions); i++) {
            inmemext = bl_access(fb->extensions, i);
            if (strcasecmp(inmemext->tablename, chunk->tablename))
                continue;
            // found it!
            gotit = TRUE;
            break;
        }
        if (!gotit && chunk->required) {
            ERROR("Couldn't find table \"%s\"", chunk->tablename);
            return -1;
        }
        table_nrows = bl_size(inmemext->items);
        table_rowsize = bl_datasize(inmemext->items);
        chunk->header = qfits_header_copy(inmemext->header);

    } else {
        //double t0;
        //t0 = timenow();
        if (find_table_column(fb, chunk->tablename, &tabstart, &tabsize, &ext)) {
            if (chunk->required)
                ERROR("Couldn't find table \"%s\" in file \"%s\"",
                      chunk->tablename, fb->filename);
            return -1;
        }
        //debug("fits_find_table_column(%s) took %g ms\n", chunk->tablename, 1000 * (timenow() - t0));

        //t0 = timenow();
        chunk->header = fitsbin_get_header(fb, ext);
        if (!chunk->header) {
            ERROR("Couldn't read FITS header from file \"%s\" extension %i", fb->filename, ext);
            return -1;
        }
        //debug("reading chunk header (%s) took %g ms\n", chunk->tablename, 1000 * (timenow() - t0));
        table_nrows = fitsbin_get_table_const(fb, ext)->nr;
        table_rowsize = fitsbin_get_table_const(fb, ext)->tab_w;
    }

    if (!chunk->itemsize)
        chunk->itemsize = table_rowsize;
    if (!chunk->nrows)
        chunk->nrows = table_nrows;

    if (chunk->callback_read_header &&
        chunk->callback_read_header(fb, chunk)) {
        ERROR("fitsbin callback_read_header failed");
        return -1;
    }

    if (chunk->nrows != table_nrows) {
        ERROR("Table %s in file %s: expected %i data items (ie, rows), found %i",
              chunk->tablename, fb->filename, chunk->nrows, table_nrows);
        return -1;
    }

    if (chunk->itemsize != table_rowsize) {
        ERROR("Table %s in file %s: expected data size %i (ie, row width in bytes), found %i",
              chunk->tablename, fb->filename, chunk->itemsize, table_rowsize);
        return -1;
    }

    if (chunk->itemsize < 0 ||
        chunk->nrows < 0 ||
        (chunk->nrows &&
         (size_t)chunk->itemsize >
             SIZE_MAX / (size_t)chunk->nrows)) {
        ERROR("Table %s in file %s has an invalid or overflowing "
              "payload size: itemsize=%i nrows=%i",
              chunk->tablename
                  ? chunk->tablename
                  : "(unnamed)",
              fb->filename
                  ? fb->filename
                  : "(memory)",
              chunk->itemsize,
              chunk->nrows);
        return -1;
    }
    expected =
        (size_t)chunk->itemsize *
        (size_t)chunk->nrows;
    if (!in_memory(fb) && fits_bytes_needed(expected) != tabsize) {
        ERROR("Expected table size (%zu => %i FITS blocks) is not equal to "
              "size of table \"%s\" (%zu => %i FITS blocks).",
              expected, fits_blocks_needed(expected),
              chunk->tablename, (size_t)tabsize,
              (int)(tabsize / (off_t)FITS_BLOCK_SIZE));
        return -1;
    }
    if (!read_payload) {
        return 0;
    }

    if (!in_memory(fb)) {
        chunk->data_file_offset = tabstart;
        chunk->data_file_size = expected;
    }

    if (in_memory(fb)) {
        int i;
        chunk->data = malloc(expected);
        for (i=0; i<chunk->nrows; i++) {
            memcpy(((char*)chunk->data) + (size_t)i * (size_t)chunk->itemsize,
                   bl_access(inmemext->items, i), chunk->itemsize);
        }
        // delete inmemext->items ?

    } else {

        get_mmap_size(tabstart, tabsize, &mapstart, &(chunk->mapsize), &mapoffset);
        mode = PROT_READ;
        flags = MAP_SHARED;
        chunk->map = mmap(0, chunk->mapsize, mode, flags, fileno(fb->fid), mapstart);
        if (chunk->map == MAP_FAILED) {
            SYSERROR("Couldn't mmap file \"%s\"", fb->filename);
            chunk->map = NULL;
            return -1;
        }
        apply_mmap_advice(fb, chunk);
        chunk->data = chunk->map + mapoffset;
    }
    return 0;
}

int fitsbin_read_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    if (read_chunk(fb, chunk, TRUE)) {
        return -1;
    }
    fitsbin_add_chunk(fb, chunk);
    return 0;
}

int fitsbin_read_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
    return read_chunk(fb, chunk, FALSE);
}

int fitsbin_read(fitsbin_t* fb) {
    int i;

    for (i=0; i<nchunks(fb); i++) {
        fitsbin_chunk_t* chunk = get_chunk(fb, i);
        if (read_chunk(fb, chunk, TRUE)) {
            if (chunk->required)
                goto bailout;
        }
    }
    return 0;

 bailout:
    return -1;
}

char* fitsbin_get_filename(const fitsbin_t* fb) {
    return fb->filename;
}

fitsbin_t* fitsbin_open_fits(anqfits_t* fits) {
    fitsbin_t* fb;
    fb = new_fitsbin(fits->filename);
    if (!fb)
        return fb;
    fb->fid = fopen(fits->filename, "rb");
    if (!fb->fid) {
        SYSERROR("Failed to open file \"%s\"", fits->filename);
        goto bailout;
    }
    if (fstat(
            fileno(fb->fid),
            &fb->open_file_stat)) {
        SYSERROR(
            "Failed to identify file \"%s\"",
            fits->filename);
        goto bailout;
    }
    fb->open_file_stat_valid = TRUE;
    fb->Next = anqfits_n_ext(fits);
    debug("N ext: %i\n", fb->Next);
    fb->fits = fits;
    fb->primheader = fitsbin_get_header(fb, 0);
    if (!fb->primheader) {
        ERROR("Couldn't read primary FITS header from file \"%s\"", fits->filename);
        goto bailout;
    }
    return fb;
 bailout:
    fitsbin_close(fb);
    return NULL;
}

fitsbin_t* fitsbin_open(const char* fn) {
    anqfits_t* fits;
    fitsbin_t* fb;

    fits = anqfits_open(fn);
    if (!fits) {
        ERROR("Failed to open file \"%s\"", fn);
        return NULL;
    }
    fb = fitsbin_open_fits(fits);
    if (!fb) {
        anqfits_close(fits);
        return NULL;
    }
    fb->owns_fits = TRUE;
    return fb;
}

fitsbin_t* fitsbin_open_in_memory() {
    fitsbin_t* fb;

    fb = new_fitsbin(NULL);
    if (!fb)
        return NULL;
    fb->primheader = qfits_table_prim_header_default();
    fb->inmemory = TRUE;
    return fb;
}

int fitsbin_switch_to_reading(fitsbin_t* fb) {
    int i;

    // clear the current chunk data??
    for (i=0; i<nchunks(fb); i++) {
        fitsbin_chunk_t* chunk = get_chunk(fb, i);
        if (chunk->header)
            qfits_header_destroy(chunk->header);
    }

    return 0;
}

fitsbin_t* fitsbin_open_for_writing(const char* fn) {
    fitsbin_t* fb;

    fb = new_fitsbin(fn);
    if (!fb)
        return NULL;
    fb->primheader = qfits_table_prim_header_default();
    fb->fid = fopen(fb->filename, "wb");
    if (!fb->fid) {
        SYSERROR("Couldn't open file \"%s\" for output", fb->filename);
        fitsbin_close(fb);
        return NULL;
    }
    return fb;
}
