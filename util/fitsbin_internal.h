/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#ifndef FITSBIN_INTERNAL_H
#define FITSBIN_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>
#include <time.h>

#include "fitsbin.h"

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define ASTROMETRY_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#define ASTROMETRY_THREAD_LOCAL __thread
#else
#error "A thread-local storage implementation is required"
#endif

#define FITSBIN_PREFETCH_RANGE_LIMIT \
    FITSBIN_MMAP_PREFETCH_RANGE_LIMIT
#define FITSBIN_PREFETCH_COPY_CHUNK (16U * 1024U)
#define FITSBIN_PAYLOAD_POPULATE_TOTAL_BUDGET \
    (16U * 1024U * 1024U)
#define FITSBIN_PAYLOAD_IO_MAX_LANES 4
#define FITSBIN_PAYLOAD_IO_MAX_JOBS 24U
#define FITSBIN_PAYLOAD_IO_MAX_BYTES \
    (64U * 1024U * 1024U)
#define FITSBIN_PAYLOAD_IO_PRIORITY_COUNT 3U
#define FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_JOBS 1U
#define FITSBIN_PAYLOAD_IO_DEMAND_RESERVED_BYTES \
    (FITSBIN_PAYLOAD_IO_MAX_BYTES / 4U)

/*
 * Join only a one-page hole inside one exact mapping, and spend no more than
 * one quarter of the ticket's exact aligned bytes on those holes. This cuts
 * small syscall fragments without recreating broad mmap readahead.
 */
#define FITSBIN_PAYLOAD_MAPPED_COALESCE_GAP_PAGES 1U
#define FITSBIN_PAYLOAD_MAPPED_COALESCE_BUDGET_DIVISOR 4U

/*
 * File-offset readahead may bridge a slightly wider hole than mapped
 * completion. This changes only storage queue shape: the exact mapped spans
 * remain the sole MADV_POPULATE_READ and READY authority.
 */
#define FITSBIN_PAYLOAD_QUEUE_COALESCE_GAP_PAGES 2U
#define FITSBIN_PAYLOAD_QUEUE_COALESCE_BUDGET_DIVISOR 8U

/*
 * Internal transport A/B. The parent source is the direct-first comparator.
 * This candidate refuses asynchronous direct tickets only for nonresident
 * sparse payload mappings, which is the production verification-sweep case;
 * the solver then exercises its existing exact mapped-population fallback.
 * Synchronous exact reads and non-sparse direct-ticket callers are unchanged.
 * This is deliberately not a public runtime option.
 */
#ifndef FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED
#define FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED 0
#endif
#if FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED != 0 && \
    FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED != 1
#error "FITSBIN_PAYLOAD_SPARSE_ASYNC_DIRECT_ENABLED must be 0 or 1"
#endif

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

int fitsbin_compare_prepared_pread_range(
    const void* left,
    const void* right);
anbool fitsbin_prepared_pread_destinations_disjoint(
    const fitsbin_prepared_pread_range_t* ranges,
    size_t range_count);

int fitsbin_payload_fd_get(fitsbin_t* fb);
fitsbin_chunk_t* fitsbin_find_data_chunk(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    size_t* chunk_offset);
int fitsbin_mapped_file_offset(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    off_t* file_offset);
void fitsbin_payload_counter_add(
    unsigned long long* counter,
    unsigned long long value);
int fitsbin_pread_all_counted(
    int fd,
    void* destination,
    size_t size,
    off_t offset,
    unsigned long long* calls,
    unsigned long long* bytes);
int fitsbin_pread_all(
    int fd,
    void* destination,
    size_t size,
    off_t offset);

int fitsbin_prepare_direct_ranges(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_prepared_pread_range_t* prepared,
    size_t* physical_bytes_out,
    size_t* logical_bytes_out,
    unsigned long long* page_count_out);
int fitsbin_compare_mapped_span(
    const void* left,
    const void* right);
int fitsbin_prepare_prefetch_spans(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_file_span_t* spans,
    size_t* span_count,
    size_t* byte_count);
int fitsbin_prepare_mapped_spans(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    unsigned long long reuse_sequence,
    fitsbin_mapped_span_t* spans,
    size_t span_capacity,
    size_t* span_count,
    size_t* byte_count,
    size_t* logical_byte_count,
    unsigned long long* page_count,
    size_t* exact_span_count,
    unsigned long long* reused_page_count,
    size_t* coalesced_gap_count,
    size_t* coalesced_gap_bytes);
int fitsbin_prepare_mapped_file_spans(
    fitsbin_t* fb,
    const fitsbin_mapped_span_t* mapped,
    size_t mapped_count,
    size_t byte_budget,
    fitsbin_file_span_t* file_spans,
    size_t file_span_capacity,
    size_t* file_span_count,
    size_t* byte_count,
    size_t* exact_span_count,
    size_t* gap_count,
    size_t* gap_bytes);
void fitsbin_payload_mark_completed_span(
    fitsbin_t* fb,
    const fitsbin_mapped_span_t* span,
    unsigned long long sequence);

unsigned long long fitsbin_payload_io_sequence_hint(void);
unsigned long long fitsbin_timespec_delta_nanoseconds(
    const struct timespec* finish,
    const struct timespec* start);
unsigned long long fitsbin_payload_io_acquire(void);
void fitsbin_payload_io_release(void);
anbool fitsbin_payload_io_planning_is_active(void);
int fitsbin_payload_io_capacity_current(void);

void fitsbin_apply_mmap_advice(
    fitsbin_t* fb,
    fitsbin_chunk_t* chunk);

#endif
