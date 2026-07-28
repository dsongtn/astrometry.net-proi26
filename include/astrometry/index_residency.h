/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef ASTROMETRY_INDEX_RESIDENCY_H
#define ASTROMETRY_INDEX_RESIDENCY_H

#include <stddef.h>
#include <stdint.h>
#include <sys/stat.h>

typedef struct index_residency index_residency_t;
typedef struct index_residency_handle index_residency_handle_t;

typedef enum index_residency_priority {
    INDEX_RESIDENCY_PRIORITY_SPECULATIVE = 0,
    INDEX_RESIDENCY_PRIORITY_LOOKAHEAD = 1,
    INDEX_RESIDENCY_PRIORITY_DEMAND = 2
} index_residency_priority_t;

typedef enum index_residency_result {
    INDEX_RESIDENCY_ERROR = -1,
    INDEX_RESIDENCY_FALLBACK = 0,
    INDEX_RESIDENCY_ACCEPTED = 1,
    INDEX_RESIDENCY_SOURCE_LEASE = 2
} index_residency_result_t;

typedef struct index_residency_stats {
    size_t byte_budget;
    size_t resident_bytes;
    size_t peak_resident_bytes;
    size_t ready_bytes;
    size_t loading_bytes;
    size_t leased_bytes;
    size_t ready_entries;
    size_t loading_entries;
    size_t failed_entries;
    size_t live_handles;
    unsigned int loader_lanes;
    int backend_supported;

    uint64_t prepare_requests;
    uint64_t acquire_requests;
    uint64_t cache_hits;
    uint64_t loading_deduplications;
    uint64_t files_queued;
    uint64_t files_copied;
    uint64_t bytes_copied;
    uint64_t budget_refusals;
    uint64_t evictions;
    uint64_t evicted_bytes;
    uint64_t source_changes;
    uint64_t copy_failures;
    uint64_t cancelled_entries;
    uint64_t source_leases;
    uint64_t source_requeues;
    uint64_t wait_count;
    uint64_t wait_nanoseconds;
} index_residency_stats_t;

/*
 * Start one job-scoped full-file residency service.
 *
 * byte_budget is an exact upper bound on the combined sizes of LOADING and
 * READY files. loader_lanes must be one or two. On platforms without the
 * sealed-memory backend, start still succeeds and later requests return
 * INDEX_RESIDENCY_FALLBACK.
 */
int index_residency_start(size_t byte_budget,
                          unsigned int loader_lanes,
                          index_residency_t** service);

/*
 * Optionally queue a complete source file. This never evicts a ready file.
 * INDEX_RESIDENCY_ACCEPTED means the source is already ready, already
 * loading, or was queued. INDEX_RESIDENCY_FALLBACK leaves the original path
 * authoritative and usable.
 */
index_residency_result_t index_residency_prepare(
    index_residency_t* service,
    const char* path,
    index_residency_priority_t priority);

/*
 * Acquire a demand lease. ACCEPTED supplies a complete resident source.
 * SOURCE_LEASE means a queued or active copy was cancelled and acknowledged;
 * the caller must use the original source while retaining the handle until
 * that source-backed index epoch closes. Release then requeues preparation.
 * FALLBACK supplies no handle and leaves the exact source path authoritative.
 * The only wait is for an active copy to acknowledge cancellation at its next
 * bounded copy-chunk boundary.
 */
index_residency_result_t index_residency_acquire(
    index_residency_t* service,
    const char* path,
    index_residency_handle_t** handle);

/*
 * The returned path names an immutable, read-only backing file and remains
 * valid until release. The source stat is the identity validated after the
 * copy and is likewise immutable for the handle lifetime.
 */
const char* index_residency_handle_path(
    const index_residency_handle_t* handle);
const struct stat* index_residency_handle_source_stat(
    const index_residency_handle_t* handle);

void index_residency_release(index_residency_handle_t* handle);

/* Wait until no queued or active file copy remains. */
int index_residency_drain(index_residency_t* service);

/* Take a lock-consistent, non-destructive statistics snapshot. */
int index_residency_get_stats(index_residency_t* service,
                              index_residency_stats_t* stats);

/*
 * Stop accepting or copying new work without joining loader lanes or
 * invalidating existing handles. This is safe before worker-pool shutdown.
 */
int index_residency_quiesce(index_residency_t* service);

/*
 * Stop accepting work, cancel queued preparation, and join loader lanes.
 * The service pointer is consumed and must not be used after this call.
 * Existing handles remain valid; final service cleanup is deferred until the
 * last handle is released.
 *
 * stop must not race a newly started prepare, acquire, drain, or stats call.
 * It may race acquire calls that were already waiting and may be called while
 * handles are live.
 */
int index_residency_stop(index_residency_t* service);

#endif
