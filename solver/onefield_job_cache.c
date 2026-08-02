/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "errors.h"
#include "fitsbin.h"
#include "index.h"
#include "log.h"
#include "mathutil.h"
#include "onefield_internal.h"
#include "os-features.h"
#include "tic.h"

#define ONEFIELD_INDEX_PREPARE_DEFER_NANOSECONDS 2000000L
#define ONEFIELD_INDEX_PREPARE_MAX_DEFER_ATTEMPTS 8U

typedef struct onefield_job_index_cache_entry {
    char* configured_path;
    index_t* index;
    struct stat identity;
    uint64_t virtual_bytes;
} onefield_job_index_cache_entry_t;

struct onefield_job_field_cache {
    anbool valid;
    int fieldnum;
    dev_t device;
    ino_t inode;
    off_t file_size;
    time_t mtime_seconds;
    long mtime_nanoseconds;
    time_t ctime_seconds;
    long ctime_nanoseconds;
    char* xcolname;
    char* ycolname;
    double pixel_xscale;
    const sip_t* predistort;
    anbool verify_uniformize;
    anbool verify_dedup;
    anbool set_crpix;
    anbool set_crpix_center;
    double crpix[2];
    double field_minx;
    double field_maxx;
    double field_miny;
    double field_maxy;
    unsigned long long reads;
    unsigned long long preprocesses;
    unsigned long long hits;
    unsigned long long invalidations;

    /*
     * The job owns one optional prepared-index handoff. Demand ownership stays
     * with its worker and ends as soon as that index task completes.
     * Descriptors are closed immediately after a coherent prepared load.
     */
    pthread_mutex_t index_mutex;
    anbool index_mutex_ready;
    pthread_cond_t index_prepare_cond;
    anbool index_prepare_cond_ready;
    pthread_t index_prepare_thread;
    anbool index_prepare_thread_ready;
    anbool index_prepare_stop;
    char* index_prepare_active_path;
    anbool index_prepare_active_stale;
    char* index_prepare_pending_path;
    onefield_job_index_cache_entry_t* index_entry;
    size_t index_entry_budget;
    uint64_t index_virtual_budget;
    uint64_t index_virtual_bytes;
    uint64_t index_virtual_peak;
    unsigned long long index_hits;
    unsigned long long index_misses;
    unsigned long long index_admitted;
    unsigned long long index_refused;
    unsigned long long index_invalidated;
    unsigned long long index_identity_retries;
    unsigned long long index_fd_close_failures;
    unsigned long long index_prepare_requests;
    unsigned long long index_prepare_started;
    unsigned long long index_prepare_completed;
    unsigned long long index_prepare_dropped;
    unsigned long long index_prepare_failures;
    unsigned long long index_prepare_waits;
    unsigned long long index_prepare_wait_timeouts;
    unsigned long long index_prepare_deferrals;
    unsigned long long index_prepare_capacity_refusals;
};


size_t onefield_internal_index_count(const onefield_t* bp) {
    return sl_size(bp->indexnames) + pl_size(bp->indexes);
}

index_t* onefield_internal_get_index(onefield_t* bp, size_t index_order) {
    if (index_order < sl_size(bp->indexnames)) {
        char* fn = sl_get(bp->indexnames, index_order);
        index_t* ind =
            onefield_internal_job_index_cache_get(bp, fn);
        if (!ind) {
            ERROR("Failed to load index %s", fn);
            exit(-1);
        }
        return ind;
    }
    index_order -= sl_size(bp->indexnames);
    return pl_get(bp->indexes, index_order);
}

char* onefield_internal_get_index_name(onefield_t* bp,
                                       size_t index_order) {
    index_t* index;
    if (index_order < sl_size(bp->indexnames)) {
        char* fn = sl_get(bp->indexnames, index_order);
        return fn;
    }
    index_order -= sl_size(bp->indexnames);
    index = pl_get(bp->indexes, index_order);
    return index->indexname;
}

int onefield_internal_done_with_index(onefield_t* bp,
                                      size_t index_order,
                                      index_t* index) {
    if (index_order < sl_size(bp->indexnames)) {
        index_free(index);
    }
    return 0;
}

static uint64_t onefield_index_cache_budget(void) {
    uint64_t budget = UINT64_MAX;
    struct rlimit address_limit;
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    long page_count;
    long page_size;
#endif

    if (sizeof(void*) < 8U) {
        return 0U;
    }
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    page_count = sysconf(_SC_PHYS_PAGES);
    page_size = sysconf(_SC_PAGESIZE);
    if (page_count <= 0 || page_size <= 0 ||
        (uint64_t)page_count >
            UINT64_MAX / (uint64_t)page_size) {
        return 0U;
    }
    budget = (uint64_t)page_count *
        (uint64_t)page_size;
    budget /= 2U;
#endif
#if defined(RLIMIT_AS)
    if (getrlimit(RLIMIT_AS, &address_limit) == 0 &&
        address_limit.rlim_cur != RLIM_INFINITY) {
        uintmax_t finite_limit =
            (uintmax_t)address_limit.rlim_cur;

        /*
         * Leave at least half of a finite address-space allowance for the
         * executable, heap, stacks, outputs, and transient solver work.
         */
        finite_limit = MIN(finite_limit, (uintmax_t)UINT64_MAX);
        finite_limit /= 2U;
        budget = MIN(
            budget,
            (uint64_t)finite_limit);
    }
#else
    (void)address_limit;
#endif
    return budget;
}

static void onefield_job_index_cache_entry_free(
    onefield_job_index_cache_entry_t* entry) {
    if (!entry) {
        return;
    }
    if (entry->index) {
        index_free(entry->index);
    }
    free(entry->configured_path);
    free(entry);
}

static int onefield_job_index_open_identity(
    index_t* index,
    struct stat* identity) {
    return index_get_source_file_stat(index, identity);
}

static anbool onefield_job_index_path_matches(
    const onefield_job_index_cache_entry_t* entry) {
    struct stat current;

    if (!entry || !entry->index ||
        !entry->index->indexfn ||
        stat(entry->index->indexfn, &current)) {
        return FALSE;
    }
    return onefield_internal_same_source_identity(
        &entry->identity,
        &current);
}

/*
 * Load one coherent index epoch and close its descriptors after mmap setup.
 *
 * The open-file fstat must match the pathname after loading. If the path was
 * replaced during acquisition, discard the entire epoch and retry once.
 */
static index_t* onefield_job_index_load_coherent(
    onefield_job_field_cache_t* cache,
    const char* configured_path,
    int index_options,
    struct stat* identity,
    anbool* descriptors_closed) {
    int acquisition_attempt;

    if (descriptors_closed) {
        *descriptors_closed = FALSE;
    }
    for (acquisition_attempt = 1;
         acquisition_attempt <= 2;
         acquisition_attempt++) {
        struct stat current;
        index_t* index =
            index_load(
                configured_path,
                index_options,
                NULL);

        if (!index) {
            return NULL;
        }
        if (onefield_job_index_open_identity(
                index,
                identity) ||
            !index->indexfn ||
            stat(index->indexfn, &current) ||
            !onefield_internal_same_source_identity(
                identity,
                &current)) {
            index_free(index);
            if (cache) {
                __atomic_add_fetch(
                    &cache->index_identity_retries,
                    1ULL,
                    __ATOMIC_RELAXED);
            }
            continue;
        }
        if (index_close_fds(index)) {
            if (cache) {
                __atomic_add_fetch(
                    &cache->index_fd_close_failures,
                    1ULL,
                    __ATOMIC_RELAXED);
            }
            /*
             * The current task may still use the completed mappings, but a
             * partially closed epoch is not admitted for cross-pass reuse.
             */
            return index;
        }
        if (descriptors_closed) {
            *descriptors_closed = TRUE;
        }
        return index;
    }
    return NULL;
}

/* index_mutex must be held. */
static anbool onefield_job_index_cache_contains_path(
    const onefield_job_field_cache_t* cache,
    const char* configured_path) {
    const onefield_job_index_cache_entry_t* entry =
        cache ? cache->index_entry : NULL;

    return configured_path && entry &&
        entry->configured_path &&
        !strcmp(entry->configured_path, configured_path);
}

/* index_mutex must be held. */
static anbool onefield_job_index_prepare_capacity_full(
    const onefield_job_field_cache_t* cache,
    uint64_t additional_bytes) {
    uint64_t available;

    if (!cache || !cache->index_entry_budget ||
        cache->index_entry ||
        cache->index_virtual_bytes >=
            cache->index_virtual_budget) {
        return TRUE;
    }
    if (!additional_bytes) {
        return FALSE;
    }
    if (additional_bytes >
        UINT64_MAX - cache->index_virtual_bytes) {
        return TRUE;
    }
    available = cache->index_virtual_budget -
        cache->index_virtual_bytes;
    return additional_bytes > available;
}

/* index_mutex must be held. */
static void onefield_job_index_prepare_record_capacity_refusal(
    onefield_job_field_cache_t* cache) {
    if (!cache) {
        return;
    }
    cache->index_prepare_capacity_refusals++;
}

/* index_mutex must be held. */
static int onefield_job_index_prepare_timedwait(
    onefield_job_field_cache_t* cache) {
    struct timespec deadline;

    if (clock_gettime(CLOCK_REALTIME, &deadline)) {
        return errno ? errno : EINVAL;
    }
    deadline.tv_nsec +=
        ONEFIELD_INDEX_PREPARE_DEFER_NANOSECONDS;
    if (deadline.tv_nsec >= 1000000000L) {
        deadline.tv_sec++;
        deadline.tv_nsec -= 1000000000L;
    }
    return pthread_cond_timedwait(
        &cache->index_prepare_cond,
        &cache->index_mutex,
        &deadline);
}

/* index_mutex must be held. */
static int onefield_job_index_prepare_defer(
    onefield_job_field_cache_t* cache) {
    cache->index_prepare_deferrals++;
    return onefield_job_index_prepare_timedwait(cache);
}

/*
 * One optional preparation lane maps exactly one near-future index at a time.
 * The mapping remains lazy: this path does not touch sparse payload pages,
 * claim solver ownership, or change index order. Descriptors are open only
 * while the mapping is built and are closed before the handoff is published.
 */
static void* onefield_job_index_prepare_main(void* opaque) {
    onefield_job_field_cache_t* cache = opaque;
    unsigned int demand_deferrals = 0U;

    if (!cache) {
        return NULL;
    }
    while (1) {
        onefield_job_index_cache_entry_t* entry = NULL;
        index_t* index = NULL;
        char* configured_path;
        char* retained_path = NULL;
        struct stat identity;
        anbool descriptors_closed = FALSE;
        anbool admitted = FALSE;
        uint64_t retained_snapshot = 0U;

        pthread_mutex_lock(&cache->index_mutex);
        while (!cache->index_prepare_stop &&
               (!cache->index_prepare_pending_path ||
                fitsbin_payload_io_demand_busy() ||
                cache->index_entry)) {
            int wait_status;

            if (cache->index_prepare_pending_path) {
                if (demand_deferrals >=
                    ONEFIELD_INDEX_PREPARE_MAX_DEFER_ATTEMPTS) {
                    free(cache->index_prepare_pending_path);
                    cache->index_prepare_pending_path = NULL;
                    cache->index_prepare_dropped++;
                    demand_deferrals = 0U;
                    continue;
                }
                wait_status =
                    onefield_job_index_prepare_defer(
                        cache);
                demand_deferrals++;
            } else {
                demand_deferrals = 0U;
                wait_status = pthread_cond_wait(
                    &cache->index_prepare_cond,
                    &cache->index_mutex);
            }

            if (wait_status &&
                wait_status != ETIMEDOUT) {
                cache->index_prepare_failures++;
                cache->index_prepare_stop = TRUE;
            }
        }
        if (cache->index_prepare_stop) {
            pthread_mutex_unlock(&cache->index_mutex);
            break;
        }
        configured_path =
            cache->index_prepare_pending_path;
        cache->index_prepare_pending_path = NULL;
        cache->index_prepare_active_path =
            configured_path;
        cache->index_prepare_active_stale = FALSE;
        cache->index_prepare_started++;
        demand_deferrals = 0U;
        pthread_mutex_unlock(&cache->index_mutex);

        /*
         * This thread exists only for parallel index preparation. Map every
         * index chunk with the production shard policy so ownership does not
         * have to repair a freshly prepared mapping.
         */
        fitsbin_mmap_set_thread_advice(
            fitsbin_mmap_advice_state_begin_pass(NULL));
        index = onefield_job_index_load_coherent(
            cache,
            configured_path,
            0,
            &identity,
            &descriptors_closed);
        fitsbin_mmap_clear_thread_advice();
        if (index && identity.st_size > 0 &&
            descriptors_closed) {
            entry = calloc(1, sizeof(*entry));
            retained_path = strdup(configured_path);
        }

        pthread_mutex_lock(&cache->index_mutex);
        if (!cache->index_prepare_stop &&
            !cache->index_prepare_active_stale &&
            !cache->index_entry &&
            entry && retained_path &&
            index && identity.st_size > 0 &&
            descriptors_closed) {
            uint64_t virtual_bytes =
                (uint64_t)identity.st_size;

            entry->configured_path = retained_path;
            entry->index = index;
            entry->identity = identity;
            entry->virtual_bytes = virtual_bytes;
            if (!onefield_job_index_prepare_capacity_full(
                    cache, virtual_bytes)) {
                cache->index_entry = entry;
                entry = NULL;
                retained_path = NULL;
                index = NULL;
                cache->index_virtual_bytes +=
                    virtual_bytes;
                cache->index_virtual_peak = MAX(
                    cache->index_virtual_peak,
                    cache->index_virtual_bytes);
                cache->index_admitted++;
                cache->index_prepare_completed++;
                retained_snapshot =
                    cache->index_virtual_bytes;
                admitted = TRUE;
            } else {
                cache->index_refused++;
                cache->index_prepare_dropped++;
                onefield_job_index_prepare_record_capacity_refusal(
                    cache);
            }
        } else if (cache->index_prepare_stop ||
                   cache->index_prepare_active_stale ||
                   cache->index_entry) {
            cache->index_prepare_dropped++;
        } else {
            cache->index_prepare_failures++;
        }
        cache->index_prepare_active_path = NULL;
        cache->index_prepare_active_stale = FALSE;
        pthread_cond_broadcast(&cache->index_prepare_cond);
        pthread_mutex_unlock(&cache->index_mutex);

        if (admitted) {
            logverb(
                "[index-shard] job-index-prepare state=ready "
                "path=%s retained=%llu budget=%llu\n",
                configured_path,
                (unsigned long long)retained_snapshot,
                (unsigned long long)
                    cache->index_virtual_budget);
        } else {
            if (entry) {
                entry->configured_path = retained_path;
                entry->index = index;
                onefield_job_index_cache_entry_free(entry);
            } else {
                free(retained_path);
                if (index) {
                    index_free(index);
                }
            }
        }
        free(configured_path);
    }
    return NULL;
}

void onefield_internal_job_index_cache_prepare(
    onefield_t* bp,
    const char* configured_path) {
    onefield_job_field_cache_t* cache;
    char* pending_path;

    if (!bp || !configured_path) {
        return;
    }
    /*
     * Full-cohort residency already owns preparation. A second mapping lane
     * could retain a source-backed index just before the resident copy becomes
     * ready, bypassing the prepared backing for the lifetime of that handoff.
     */
    if (index_residency_service_active()) {
        return;
    }
    cache = bp->job_field_cache;
    if (!cache || !cache->index_mutex_ready ||
        !cache->index_prepare_cond_ready ||
        !cache->index_prepare_thread_ready ||
        bp->index_options != 0) {
        return;
    }
    pending_path = strdup(configured_path);
    pthread_mutex_lock(&cache->index_mutex);
    cache->index_prepare_requests++;
    if (!pending_path) {
        cache->index_prepare_failures++;
    } else if (cache->index_prepare_stop) {
        cache->index_prepare_dropped++;
    } else if (onefield_job_index_prepare_capacity_full(
                   cache, 0U)) {
        cache->index_prepare_dropped++;
        onefield_job_index_prepare_record_capacity_refusal(
            cache);
    } else if (onefield_job_index_cache_contains_path(
                   cache, configured_path) ||
               (cache->index_prepare_active_path &&
                !cache->index_prepare_active_stale &&
                !strcmp(cache->index_prepare_active_path,
                        configured_path)) ||
               (cache->index_prepare_pending_path &&
                !strcmp(cache->index_prepare_pending_path,
                        configured_path))) {
        cache->index_prepare_dropped++;
    } else if (cache->index_prepare_pending_path) {
        cache->index_prepare_dropped++;
    } else {
        cache->index_prepare_pending_path = pending_path;
        pending_path = NULL;
        pthread_cond_signal(&cache->index_prepare_cond);
    }
    pthread_mutex_unlock(&cache->index_mutex);
    free(pending_path);
}

void onefield_internal_job_index_cache_flush(onefield_t* bp) {
    onefield_job_field_cache_t* cache;
    onefield_job_index_cache_entry_t* entry = NULL;

    if (!bp) {
        return;
    }
    cache = bp->job_field_cache;
    if (!cache || !cache->index_mutex_ready) {
        return;
    }
    pthread_mutex_lock(&cache->index_mutex);
    if (cache->index_prepare_pending_path) {
        free(cache->index_prepare_pending_path);
        cache->index_prepare_pending_path = NULL;
        cache->index_prepare_dropped++;
    }
    if (cache->index_prepare_active_path) {
        cache->index_prepare_active_stale = TRUE;
    }
    if (cache->index_entry) {
        entry = cache->index_entry;
        cache->index_entry = NULL;
        cache->index_virtual_bytes -= MIN(
            cache->index_virtual_bytes,
            entry->virtual_bytes);
        cache->index_prepare_dropped++;
    }
    if (cache->index_prepare_cond_ready) {
        pthread_cond_broadcast(
            &cache->index_prepare_cond);
    }
    pthread_mutex_unlock(&cache->index_mutex);
    onefield_job_index_cache_entry_free(entry);
}

static index_t* onefield_job_index_cache_take_entry(
    onefield_job_index_cache_entry_t* entry) {
    index_t* index;

    if (!entry) {
        return NULL;
    }
    index = entry->index;
    entry->index = NULL;
    return index;
}
index_t* onefield_internal_job_index_cache_get(
    onefield_t* bp,
    const char* configured_path) {
    onefield_job_field_cache_t* cache;
    onefield_job_index_cache_entry_t* entry;
    index_t* index;
    anbool waited = FALSE;

    if (!bp || !configured_path) {
        return NULL;
    }
    cache = bp->job_field_cache;
    if (!cache || !cache->index_mutex_ready ||
        cache->index_virtual_budget == 0U ||
        bp->index_options != 0) {
        return index_load(
            configured_path,
            bp->index_options,
            NULL);
    }

retry_lookup:
    pthread_mutex_lock(&cache->index_mutex);
    entry = cache->index_entry;
    if (entry && entry->configured_path &&
        !strcmp(entry->configured_path,
                configured_path)) {
        if (!onefield_job_index_path_matches(entry)) {
            cache->index_invalidated++;
            cache->index_virtual_bytes -= MIN(
                cache->index_virtual_bytes,
                entry->virtual_bytes);
            cache->index_entry = NULL;
            pthread_cond_broadcast(&cache->index_prepare_cond);
            pthread_mutex_unlock(
                &cache->index_mutex);
            logverb(
                "[index-shard] job-index-cache "
                "state=invalidate path=%s\n",
                configured_path);
            onefield_job_index_cache_entry_free(entry);
            goto retry_lookup;
        }
        cache->index_entry = NULL;
        cache->index_virtual_bytes -= MIN(
            cache->index_virtual_bytes,
            entry->virtual_bytes);
        cache->index_hits++;
        index = onefield_job_index_cache_take_entry(
            entry);
        entry->virtual_bytes = 0U;
        {
            unsigned long long hit_count =
                cache->index_hits;

            pthread_cond_broadcast(
                &cache->index_prepare_cond);
            pthread_mutex_unlock(
                &cache->index_mutex);
            logverb(
                "[index-shard] job-index-cache "
                "state=hit path=%s hits=%llu\n",
                configured_path,
                hit_count);
        }
        free(entry->configured_path);
        free(entry);
        return index;
    }
    if (cache->index_prepare_active_path &&
        cache->index_prepare_cond_ready &&
        !cache->index_prepare_active_stale &&
        !strcmp(cache->index_prepare_active_path,
                configured_path) &&
        !waited) {
        int wait_status;

        cache->index_prepare_waits++;
        wait_status =
            onefield_job_index_prepare_timedwait(cache);
        waited = TRUE;
        if (wait_status == ETIMEDOUT) {
            cache->index_prepare_wait_timeouts++;
        } else if (wait_status) {
            cache->index_prepare_failures++;
        }
        pthread_mutex_unlock(&cache->index_mutex);
        goto retry_lookup;
    }
    if (cache->index_prepare_active_path &&
        !strcmp(cache->index_prepare_active_path,
                configured_path)) {
        cache->index_prepare_active_stale = TRUE;
    }
    if (cache->index_prepare_pending_path &&
        !strcmp(cache->index_prepare_pending_path,
                configured_path)) {
        free(cache->index_prepare_pending_path);
        cache->index_prepare_pending_path = NULL;
        cache->index_prepare_dropped++;
        pthread_cond_signal(&cache->index_prepare_cond);
    }
    cache->index_misses++;
    pthread_mutex_unlock(&cache->index_mutex);

    return index_load(
        configured_path,
        bp->index_options,
        NULL);
}

int onefield_job_index_cache_test_handoff_state(void) {
    onefield_job_index_cache_entry_t entry;
    index_t retained;

    memset(&entry, 0, sizeof(entry));
    memset(&retained, 0, sizeof(retained));
    entry.index = &retained;

    if (onefield_job_index_cache_take_entry(
            &entry) != &retained ||
        entry.index) {
        return -1;
    }
    if (onefield_job_index_cache_take_entry(&entry)) {
        return -1;
    }
    return 0;
}

static void onefield_field_cache_clear_key(
    onefield_job_field_cache_t* cache) {
    if (!cache) {
        return;
    }
    free(cache->xcolname);
    free(cache->ycolname);
    cache->xcolname = NULL;
    cache->ycolname = NULL;
    cache->valid = FALSE;
}

int onefield_job_field_cache_begin(onefield_t* bp) {
    onefield_job_field_cache_t* cache;

    if (!bp) {
        return -1;
    }
    if (bp->job_field_cache) {
        return 0;
    }
    bp->job_field_cache =
        calloc(1, sizeof(*bp->job_field_cache));
    if (!bp->job_field_cache) {
        /*
         * Retention is optional. Callers may continue through the exact
         * per-run field lifecycle if metadata allocation is unavailable.
         */
        logverb("[index-shard] job-field-cache state=disabled "
                "reason=allocation\n");
        return 0;
    }
    cache = bp->job_field_cache;
    /*
     * Dynamic index claims cannot be predicted by the old one-entry handoff.
     * Preparing a guessed path can duplicate the owner's demand load and
     * compete with current-index payload delivery. Keep this lane dormant;
     * bounded preparation must be driven by an actual reserved claim.
     */
    cache->index_entry_budget = 0U;
    cache->index_virtual_budget =
        cache->index_entry_budget
            ? onefield_index_cache_budget()
            : 0U;
    if (pthread_mutex_init(
            &cache->index_mutex,
            NULL) == 0) {
        cache->index_mutex_ready = TRUE;
    } else {
        cache->index_entry_budget = 0U;
        cache->index_virtual_budget = 0U;
    }
    if (cache->index_mutex_ready &&
        cache->index_entry_budget > 0U &&
        cache->index_virtual_budget > 0U &&
        bp->index_options == 0 &&
        bp->index_shard_workers > 1 &&
        pthread_cond_init(
            &cache->index_prepare_cond,
            NULL) == 0) {
        cache->index_prepare_cond_ready = TRUE;
        if (pthread_create(
                &cache->index_prepare_thread,
                NULL,
                onefield_job_index_prepare_main,
                cache) == 0) {
            cache->index_prepare_thread_ready = TRUE;
        } else {
            pthread_cond_destroy(
                &cache->index_prepare_cond);
            cache->index_prepare_cond_ready = FALSE;
        }
    }
    if (!cache->index_prepare_thread_ready) {
        cache->index_entry_budget = 0U;
        cache->index_virtual_budget = 0U;
    }
    logverb("[index-shard] job-field-cache state=begin "
            "index_virtual_budget=%llu index_entry_budget=%zu "
            "index_prepare=%s\n",
            (unsigned long long)
                cache->index_virtual_budget,
            cache->index_entry_budget,
            cache->index_prepare_thread_ready
                ? "enabled" : "disabled");
    return 0;
}

void onefield_job_field_cache_invalidate(onefield_t* bp) {
    onefield_job_field_cache_t* cache;

    if (!bp) {
        return;
    }
    cache = bp->job_field_cache;
    if (cache && cache->valid) {
        cache->invalidations++;
    }
    solver_cleanup_field(&bp->solver);
    if (cache) {
        onefield_field_cache_clear_key(cache);
    }
    if (bp->xyls) {
        xylist_close(bp->xyls);
        bp->xyls = NULL;
    }
}

void onefield_job_field_cache_end(onefield_t* bp) {
    onefield_job_field_cache_t* cache;
    onefield_job_index_cache_entry_t* entry = NULL;

    if (!bp || !bp->job_field_cache) {
        return;
    }
    cache = bp->job_field_cache;
    logverb("[index-shard] job-field-cache state=end "
            "reads=%llu preprocesses=%llu hits=%llu "
            "invalidations=%llu\n",
            cache->reads,
            cache->preprocesses,
            cache->hits,
            cache->invalidations);
    onefield_job_field_cache_invalidate(bp);
    if (cache->index_prepare_thread_ready) {
        int join_status;

        pthread_mutex_lock(&cache->index_mutex);
        cache->index_prepare_stop = TRUE;
        if (cache->index_prepare_pending_path) {
            free(cache->index_prepare_pending_path);
            cache->index_prepare_pending_path = NULL;
            cache->index_prepare_dropped++;
        }
        pthread_cond_broadcast(&cache->index_prepare_cond);
        pthread_mutex_unlock(&cache->index_mutex);
        join_status = pthread_join(
            cache->index_prepare_thread,
            NULL);
        if (join_status) {
            logerr("[index-shard] failed to join job-index-prepare "
                   "status=%i\n",
                   join_status);
            return;
        }
        cache->index_prepare_thread_ready = FALSE;
    }
    logverb(
        "[index-shard] job-index-cache state=end "
        "hits=%llu misses=%llu admitted=%llu "
        "refused=%llu invalidated=%llu retries=%llu "
        "fd_close_failures=%llu retained=%llu peak=%llu "
        "budget=%llu entries=%zu entry_budget=%zu "
        "prepare_requests=%llu "
        "prepare_started=%llu prepare_completed=%llu "
        "prepare_dropped=%llu prepare_failures=%llu "
        "prepare_waits=%llu prepare_wait_timeouts=%llu "
        "prepare_deferrals=%llu "
        "prepare_capacity_refusals=%llu\n",
        cache->index_hits,
        cache->index_misses,
        cache->index_admitted,
        cache->index_refused,
        cache->index_invalidated,
        cache->index_identity_retries,
        cache->index_fd_close_failures,
        (unsigned long long)
            cache->index_virtual_bytes,
        (unsigned long long)
            cache->index_virtual_peak,
        (unsigned long long)
            cache->index_virtual_budget,
        cache->index_entry ? (size_t)1U : (size_t)0U,
        cache->index_entry_budget,
        cache->index_prepare_requests,
        cache->index_prepare_started,
        cache->index_prepare_completed,
        cache->index_prepare_dropped,
        cache->index_prepare_failures,
        cache->index_prepare_waits,
        cache->index_prepare_wait_timeouts,
        cache->index_prepare_deferrals,
        cache->index_prepare_capacity_refusals);
    if (cache->index_prepare_cond_ready) {
        pthread_cond_destroy(
            &cache->index_prepare_cond);
        cache->index_prepare_cond_ready = FALSE;
    }
    if (cache->index_mutex_ready) {
        pthread_mutex_lock(&cache->index_mutex);
        entry = cache->index_entry;
        cache->index_entry = NULL;
        cache->index_virtual_bytes = 0U;
        pthread_mutex_unlock(&cache->index_mutex);
        onefield_job_index_cache_entry_free(entry);
        pthread_mutex_destroy(&cache->index_mutex);
        cache->index_mutex_ready = FALSE;
    } else {
        onefield_job_index_cache_entry_free(
            cache->index_entry);
        cache->index_entry = NULL;
    }
    free(cache);
    bp->job_field_cache = NULL;
}

int onefield_internal_open_master_xyls(onefield_t* bp) {
    if (!bp || !bp->fieldfname) {
        return -1;
    }
    if (bp->xyls) {
        return 0;
    }

    logverb("Reading fields file %s...", bp->fieldfname);
    bp->xyls = xylist_open(bp->fieldfname);
    if (!bp->xyls) {
        ERROR("Failed to read xylist.\n");
        return -1;
    }
    xylist_set_xname(bp->xyls, bp->xcolname);
    xylist_set_yname(bp->xyls, bp->ycolname);
    xylist_set_include_flux(bp->xyls, FALSE);
    xylist_set_include_background(bp->xyls, FALSE);
    logverb("found %u fields.\n", xylist_n_fields(bp->xyls));
    return 0;
}

static void onefield_discard_field_acquisition(onefield_t* bp) {
    if (!bp) {
        return;
    }
    if (bp->job_field_cache) {
        onefield_job_field_cache_invalidate(bp);
        return;
    }
    solver_cleanup_field(&bp->solver);
    if (bp->xyls) {
        xylist_close(bp->xyls);
        bp->xyls = NULL;
    }
}

int onefield_internal_validate_single_field_list(onefield_t* bp) {
    int acquisition_attempt;

    if (!bp || il_size(bp->fieldlist) != 1) {
        return 0;
    }
    for (acquisition_attempt = 1;
         acquisition_attempt <= 2;
         acquisition_attempt++) {
        struct stat source_stat;
        struct stat source_stat_after;
        xylist_t* probe;
        int field_count;

        if (stat(bp->fieldfname, &source_stat)) {
            if (acquisition_attempt < 2) {
                continue;
            }
            logerr("Failed to identify XYLS source %s.\n",
                   bp->fieldfname);
            return -1;
        }
        probe = xylist_open(bp->fieldfname);
        if (!probe) {
            if (acquisition_attempt < 2) {
                continue;
            }
            logerr("Failed to inspect XYLS source %s.\n",
                   bp->fieldfname);
            return -1;
        }
        field_count = xylist_n_fields(probe);
        xylist_close(probe);
        if (stat(bp->fieldfname, &source_stat_after) ||
            !onefield_internal_same_source_identity(
                &source_stat,
                &source_stat_after)) {
            if (acquisition_attempt < 2) {
                continue;
            }
            logerr("XYLS source changed during both field-list "
                   "validation attempts.\n");
            return -1;
        }
        onefield_internal_remove_invalid_fields(bp->fieldlist, field_count);
        return 0;
    }
    return -1;
}

static void onefield_stat_times(
    const struct stat* source_stat,
    time_t* mtime_seconds,
    long* mtime_nanoseconds,
    time_t* ctime_seconds,
    long* ctime_nanoseconds) {
    *mtime_seconds = source_stat->st_mtime;
    *ctime_seconds = source_stat->st_ctime;
#if defined(__APPLE__)
    *mtime_nanoseconds = source_stat->st_mtimespec.tv_nsec;
    *ctime_nanoseconds = source_stat->st_ctimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
    *mtime_nanoseconds = source_stat->st_mtim.tv_nsec;
    *ctime_nanoseconds = source_stat->st_ctim.tv_nsec;
#else
    *mtime_nanoseconds = 0L;
    *ctime_nanoseconds = 0L;
#endif
}

anbool onefield_internal_same_source_identity(
    const struct stat* first,
    const struct stat* second) {
    time_t first_mtime_seconds;
    time_t second_mtime_seconds;
    time_t first_ctime_seconds;
    time_t second_ctime_seconds;
    long first_mtime_nanoseconds;
    long second_mtime_nanoseconds;
    long first_ctime_nanoseconds;
    long second_ctime_nanoseconds;

    if (!first || !second) {
        return FALSE;
    }
    onefield_stat_times(
        first,
        &first_mtime_seconds,
        &first_mtime_nanoseconds,
        &first_ctime_seconds,
        &first_ctime_nanoseconds);
    onefield_stat_times(
        second,
        &second_mtime_seconds,
        &second_mtime_nanoseconds,
        &second_ctime_seconds,
        &second_ctime_nanoseconds);
    return first->st_dev == second->st_dev &&
        first->st_ino == second->st_ino &&
        first->st_size == second->st_size &&
        first_mtime_seconds == second_mtime_seconds &&
        first_mtime_nanoseconds == second_mtime_nanoseconds &&
        first_ctime_seconds == second_ctime_seconds &&
        first_ctime_nanoseconds == second_ctime_nanoseconds;
}

anbool onefield_field_cache_key_matches(
    const onefield_t* bp,
    int fieldnum,
    const struct stat* source_stat) {
    const onefield_job_field_cache_t* cache;
    const solver_t* sp;
    time_t mtime_seconds;
    time_t ctime_seconds;
    long mtime_nanoseconds;
    long ctime_nanoseconds;

    if (!bp || !source_stat || !bp->job_field_cache) {
        return FALSE;
    }
    cache = bp->job_field_cache;
    sp = &bp->solver;
    if (!cache->valid || !cache->xcolname || !cache->ycolname) {
        return FALSE;
    }
    onefield_stat_times(
        source_stat,
        &mtime_seconds,
        &mtime_nanoseconds,
        &ctime_seconds,
        &ctime_nanoseconds);
    return cache->fieldnum == fieldnum &&
        cache->device == source_stat->st_dev &&
        cache->inode == source_stat->st_ino &&
        cache->file_size == source_stat->st_size &&
        cache->mtime_seconds == mtime_seconds &&
        cache->mtime_nanoseconds == mtime_nanoseconds &&
        cache->ctime_seconds == ctime_seconds &&
        cache->ctime_nanoseconds == ctime_nanoseconds &&
        !strcmp(cache->xcolname,
                bp->xcolname ? bp->xcolname : "") &&
        !strcmp(cache->ycolname,
                bp->ycolname ? bp->ycolname : "") &&
        cache->pixel_xscale == sp->pixel_xscale &&
        cache->predistort == sp->predistort &&
        cache->verify_uniformize == sp->verify_uniformize &&
        cache->verify_dedup == sp->verify_dedup &&
        cache->set_crpix == sp->set_crpix &&
        cache->set_crpix_center == sp->set_crpix_center &&
        cache->crpix[0] == sp->crpix[0] &&
        cache->crpix[1] == sp->crpix[1] &&
        cache->field_minx == sp->field_minx &&
        cache->field_maxx == sp->field_maxx &&
        cache->field_miny == sp->field_miny &&
        cache->field_maxy == sp->field_maxy;
}

static anbool onefield_field_cache_record_key(
    onefield_t* bp,
    int fieldnum,
    const struct stat* source_stat) {
    onefield_job_field_cache_t* cache = bp->job_field_cache;
    solver_t* sp = &bp->solver;
    char* xcolname;
    char* ycolname;

    if (!cache || !source_stat) {
        return FALSE;
    }
    xcolname = strdup(bp->xcolname ? bp->xcolname : "");
    ycolname = strdup(bp->ycolname ? bp->ycolname : "");
    if (!xcolname || !ycolname) {
        free(xcolname);
        free(ycolname);
        onefield_field_cache_clear_key(cache);
        logverb("[index-shard] job-field-cache state=disabled "
                "reason=key-allocation\n");
        return FALSE;
    }
    onefield_field_cache_clear_key(cache);
    cache->xcolname = xcolname;
    cache->ycolname = ycolname;
    cache->fieldnum = fieldnum;
    cache->device = source_stat->st_dev;
    cache->inode = source_stat->st_ino;
    cache->file_size = source_stat->st_size;
    onefield_stat_times(
        source_stat,
        &cache->mtime_seconds,
        &cache->mtime_nanoseconds,
        &cache->ctime_seconds,
        &cache->ctime_nanoseconds);
    cache->pixel_xscale = sp->pixel_xscale;
    cache->predistort = sp->predistort;
    cache->verify_uniformize = sp->verify_uniformize;
    cache->verify_dedup = sp->verify_dedup;
    cache->set_crpix = sp->set_crpix;
    cache->set_crpix_center = sp->set_crpix_center;
    cache->crpix[0] = sp->crpix[0];
    cache->crpix[1] = sp->crpix[1];
    cache->field_minx = sp->field_minx;
    cache->field_maxx = sp->field_maxx;
    cache->field_miny = sp->field_miny;
    cache->field_maxy = sp->field_maxy;
    cache->valid = TRUE;
    return TRUE;
}

void onefield_internal_reset_field_pass_state(onefield_t* bp) {
    solver_t* sp;

    if (!bp) {
        return;
    }
    sp = &bp->solver;
    solver_reset_best_match(sp);
    solver_reset_counters(sp);
    sp->index = NULL;
    sp->mo_template = NULL;
    sp->record_match_callback = NULL;
    sp->timer_callback = NULL;
    sp->userdata = NULL;
    memset(&sp->profile, 0, sizeof(sp->profile));
}

int onefield_internal_prepare_field_view(
    onefield_t* bp,
    int fieldnum,
    double* field_read_seconds,
    double* preprocess_seconds) {
    onefield_job_field_cache_t* cache;
    solver_t* sp;
    int acquisition_attempt;

    if (!bp || !field_read_seconds || !preprocess_seconds) {
        return -1;
    }
    *field_read_seconds = 0.0;
    *preprocess_seconds = 0.0;
    sp = &bp->solver;
    cache = bp->job_field_cache;

    /*
     * qfits table reads reopen the source pathname. Bracket the complete
     * metadata/column/preprocess acquisition and retry once from a freshly
     * opened XYLS object if the pathname identity changes. No cache key is
     * published until the closing stat matches the opening stat.
     */
    for (acquisition_attempt = 1;
         acquisition_attempt <= 2;
         acquisition_attempt++) {
        struct stat source_stat;
        anbool retainable =
            cache && il_size(bp->fieldlist) == 1;
        anbool source_stat_valid = FALSE;
        double phase_wall_start;

        if (!stat(bp->fieldfname, &source_stat)) {
            source_stat_valid = TRUE;
        } else if (retainable) {
            logverb("[onefield] job-field-cache state=retry "
                    "reason=source-identity field=%i attempt=%i\n",
                    fieldnum,
                    acquisition_attempt);
            onefield_discard_field_acquisition(bp);
            if (acquisition_attempt < 2) {
                continue;
            }
            logerr("Failed to identify XYLS source for field %i "
                   "during both acquisition attempts.\n",
                   fieldnum);
            return -1;
        }

        if (retainable &&
            source_stat_valid &&
            onefield_field_cache_key_matches(
                bp, fieldnum, &source_stat) &&
            sp->fieldxy_orig && sp->fieldxy && sp->vf) {
            struct stat source_stat_after;

            if (onefield_internal_open_master_xyls(bp) ||
                xylist_open_field(bp->xyls, fieldnum)) {
                logerr("Failed to reopen extension %i in xylist.\n",
                       fieldnum);
                onefield_discard_field_acquisition(bp);
                if (acquisition_attempt < 2) {
                    continue;
                }
                return -1;
            }
            if (stat(bp->fieldfname, &source_stat_after) ||
                !onefield_internal_same_source_identity(
                    &source_stat,
                    &source_stat_after)) {
                logverb("[index-shard] job-field-cache state=retry "
                        "reason=source-changed-during-hit field=%i "
                        "attempt=%i\n",
                        fieldnum,
                        acquisition_attempt);
                onefield_discard_field_acquisition(bp);
                if (acquisition_attempt < 2) {
                    continue;
                }
                logerr("XYLS source changed during both acquisition "
                       "attempts for field %i.\n",
                       fieldnum);
                return -1;
            }
            cache->hits++;
            solver_release_incompatible_field_geometry(sp);
            onefield_internal_reset_field_pass_state(bp);
            logverb("[index-shard] job-field-cache state=hit field=%i "
                    "hits=%llu\n",
                    fieldnum,
                    cache->hits);
            return 0;
        }

        if (cache && cache->valid) {
            logverb("[index-shard] job-field-cache state=invalidate "
                    "reason=identity-or-preprocess-key\n");
        }
        onefield_discard_field_acquisition(bp);
        if (onefield_internal_open_master_xyls(bp)) {
            if (acquisition_attempt < 2) {
                continue;
            }
            return -1;
        }
        phase_wall_start = monotonic_seconds();
        if (xylist_open_field(bp->xyls, fieldnum)) {
            logerr("Failed to open extension %i in xylist.\n",
                   fieldnum);
            onefield_discard_field_acquisition(bp);
            if (acquisition_attempt < 2) {
                continue;
            }
            return -1;
        }
        solver_set_field(sp, xylist_read_field(bp->xyls, NULL));
        *field_read_seconds +=
            monotonic_seconds() - phase_wall_start;
        if (!sp->fieldxy_orig) {
            logerr("Failed to read xylist field.\n");
            onefield_discard_field_acquisition(bp);
            if (acquisition_attempt < 2) {
                continue;
            }
            return -1;
        }
        if (cache) {
            cache->reads++;
        }

        phase_wall_start = monotonic_seconds();
        solver_preprocess_field(sp);
        *preprocess_seconds +=
            monotonic_seconds() - phase_wall_start;
        if (!sp->fieldxy || !sp->vf) {
            logerr("Failed to preprocess xylist field.\n");
            onefield_discard_field_acquisition(bp);
            return -1;
        }
        solver_release_incompatible_field_geometry(sp);
        if (cache) {
            cache->preprocesses++;
        }
        onefield_internal_reset_field_pass_state(bp);
        if (source_stat_valid) {
            struct stat source_stat_after;

            if (stat(bp->fieldfname, &source_stat_after) ||
                !onefield_internal_same_source_identity(
                    &source_stat,
                    &source_stat_after)) {
                logverb("[index-shard] job-field-cache state=retry "
                    "reason=source-changed-during-fill field=%i\n",
                    fieldnum);
                onefield_discard_field_acquisition(bp);
                if (acquisition_attempt < 2) {
                    continue;
                }
                logerr("XYLS source changed during both acquisition "
                       "attempts for field %i.\n",
                       fieldnum);
                return -1;
            }
            if (retainable) {
                if (!onefield_field_cache_record_key(
                        bp, fieldnum, &source_stat_after)) {
                    retainable = FALSE;
                }
            }
        }
        logverb("[index-shard] job-field-cache state=fill field=%i "
                "retained=%i read=%.6f preprocess=%.6f\n",
                fieldnum,
                retainable ? 1 : 0,
                *field_read_seconds,
                *preprocess_seconds);
        return 0;
    }
    return -1;
}



anbool onefield_internal_field_cache_valid(const onefield_t* bp) {
    return bp && bp->job_field_cache &&
        bp->job_field_cache->valid;
}

anbool onefield_internal_field_cache_has_field(const onefield_t* bp,
                                               int fieldnum) {
    return onefield_internal_field_cache_valid(bp) &&
        bp->job_field_cache->fieldnum == fieldnum;
}
