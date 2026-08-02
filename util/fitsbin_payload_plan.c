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

int fitsbin_prepare_direct_ranges(
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
        if (fitsbin_mapped_file_offset(
            fb,
            ranges[i].data,
            ranges[i].size,
            &prepared[i].offset)) {
            return -1;
        }
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

int fitsbin_compare_mapped_span(
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

static fitsbin_chunk_t* fitsbin_payload_mapping_chunk(
    fitsbin_t* fb,
    uintptr_t map_begin,
    uintptr_t map_end) {
    int i;

    if (!fb || map_end <= map_begin) {
        return NULL;
    }
    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);

        if (!chunk || chunk->mmap_region !=
                FITSBIN_MMAP_REGION_PAYLOAD ||
            !chunk->map || !chunk->mapsize) {
            continue;
        }
        if ((uintptr_t)chunk->map == map_begin &&
            chunk->mapsize == (size_t)(map_end - map_begin)) {
            return chunk;
        }
    }
    return NULL;
}

static unsigned int* fitsbin_payload_page_sequences_get(
    fitsbin_t* fb,
    fitsbin_chunk_t* chunk,
    size_t page_size,
    anbool create) {
    unsigned int* sequences;
    unsigned int* candidate;
    unsigned int* expected;
    size_t page_count;

    if (!fb || !chunk || !chunk->map || !chunk->mapsize ||
        !page_size) {
        return NULL;
    }
    sequences = __atomic_load_n(
        &chunk->payload_page_sequences, __ATOMIC_ACQUIRE);
    if (sequences || !create) {
        return sequences;
    }
    page_count = chunk->mapsize / page_size;
    if (chunk->mapsize % page_size) {
        page_count++;
    }
    if (!page_count ||
        page_count > SIZE_MAX / sizeof(*candidate)) {
        return NULL;
    }
    candidate = calloc(page_count, sizeof(*candidate));
    if (!candidate) {
        return NULL;
    }
    expected = NULL;
    if (__atomic_compare_exchange_n(
            &chunk->payload_page_sequences,
            &expected,
            candidate,
            FALSE,
            __ATOMIC_RELEASE,
            __ATOMIC_ACQUIRE)) {
        __atomic_add_fetch(
            &fb->payload_cache_allocations,
            1ULL,
            __ATOMIC_RELAXED);
        return candidate;
    }
    free(candidate);
    return expected;
}

/*
 * Drop only pages completed recently for the exact VMA. A mapping-lifetime
 * bit is not sufficient because the kernel may reclaim a populated page.
 * Limit reuse to the service admission window so READY never relies on an
 * arbitrarily old completion. Native mapped access remains authoritative.
 * Return SIZE_MAX when split output would overflow.
 */
static size_t fitsbin_payload_filter_completed_pages(
    fitsbin_t* fb,
    fitsbin_mapped_span_t* spans,
    size_t span_count,
    size_t span_capacity,
    size_t page_size,
    unsigned long long sequence,
    unsigned long long* reused_pages) {
    fitsbin_mapped_span_t filtered[FITSBIN_PREFETCH_RANGE_LIMIT];
    unsigned long long hits = 0ULL;
    unsigned long long misses = 0ULL;
    size_t filtered_count = 0U;
    size_t i;

    if (!fb || !spans || span_count > span_capacity ||
        span_capacity > FITSBIN_PREFETCH_RANGE_LIMIT || !page_size ||
        !reused_pages) {
        return SIZE_MAX;
    }
    *reused_pages = 0ULL;
    for (i = 0U; i < span_count; i++) {
        fitsbin_chunk_t* chunk = fitsbin_payload_mapping_chunk(
            fb, spans[i].map_begin, spans[i].map_end);
        unsigned int* sequences =
            fitsbin_payload_page_sequences_get(
            fb, chunk, page_size, TRUE);
        uintptr_t cursor = spans[i].begin;
        unsigned int current_sequence = (unsigned int)sequence;

        while (cursor < spans[i].end) {
            uintptr_t next = spans[i].end - cursor < page_size
                ? spans[i].end
                : cursor + page_size;
            size_t page_index = (size_t)(
                (cursor - spans[i].map_begin) / page_size);
            unsigned int completed_sequence = sequences
                ? __atomic_load_n(
                    &sequences[page_index], __ATOMIC_ACQUIRE)
                : 0U;
            anbool completed = current_sequence &&
                completed_sequence &&
                (unsigned int)(
                    current_sequence - completed_sequence) <
                    FITSBIN_PAYLOAD_IO_MAX_JOBS;

            if (completed) {
                hits++;
            } else {
                fitsbin_mapped_span_t* output;

                misses++;
                if (filtered_count &&
                    filtered[filtered_count - 1U].map_begin ==
                        spans[i].map_begin &&
                    filtered[filtered_count - 1U].map_end ==
                        spans[i].map_end &&
                    filtered[filtered_count - 1U].end == cursor) {
                    filtered[filtered_count - 1U].end = next;
                    cursor = next;
                    continue;
                }
                if (filtered_count >= span_capacity) {
                    return SIZE_MAX;
                }
                output = &filtered[filtered_count++];
                output->map_begin = spans[i].map_begin;
                output->map_end = spans[i].map_end;
                output->begin = cursor;
                output->end = next;
            }
            cursor = next;
        }
    }
    memcpy(spans, filtered, filtered_count * sizeof(*spans));
    *reused_pages = hits;
    __atomic_add_fetch(
        &fb->payload_cache_hits, hits, __ATOMIC_RELAXED);
    __atomic_add_fetch(
        &fb->payload_cache_misses, misses, __ATOMIC_RELAXED);
    return filtered_count;
}

void fitsbin_payload_mark_completed_span(
    fitsbin_t* fb,
    const fitsbin_mapped_span_t* span,
    unsigned long long sequence) {
    fitsbin_chunk_t* chunk;
    unsigned int* sequences;
    unsigned int completed_sequence = (unsigned int)sequence;
    size_t page_size;
    uintptr_t cursor;

    if (!fb || !span || span->end <= span->begin ||
        !completed_sequence) {
        return;
    }
    page_size = fb->mmap_page_size;
    if (!page_size) {
        return;
    }
    chunk = fitsbin_payload_mapping_chunk(
        fb, span->map_begin, span->map_end);
    sequences = fitsbin_payload_page_sequences_get(
        fb, chunk, page_size, TRUE);
    if (!sequences) {
        return;
    }
    cursor = span->begin;
    while (cursor < span->end) {
        size_t page_index = (size_t)(
            (cursor - span->map_begin) / page_size);
        __atomic_store_n(
            &sequences[page_index],
            completed_sequence,
            __ATOMIC_RELEASE);
        if (span->end - cursor <= page_size) {
            break;
        }
        cursor += page_size;
    }
}

int fitsbin_prepare_prefetch_spans(
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
    size_t* coalesced_gap_bytes) {
    size_t accepted = 0U;
    size_t merged = 0U;
    size_t coalesced = 0U;
    size_t aligned_bytes = 0U;
    size_t exact_bytes = 0U;
    size_t gap_budget = 0U;
    size_t gap_bytes = 0U;
    size_t gap_merges = 0U;
    size_t logical_bytes = 0U;
    unsigned long long pages = 0ULL;
    size_t page_size;
    size_t i;
    size_t original_merged;
    size_t filtered;
    unsigned long long reused_pages = 0ULL;

    if (!span_count || !byte_count || !logical_byte_count ||
        !page_count || !exact_span_count ||
        !reused_page_count ||
        !coalesced_gap_count || !coalesced_gap_bytes) {
        errno = EINVAL;
        return -1;
    }
    *span_count = 0U;
    *byte_count = 0U;
    *logical_byte_count = 0U;
    *page_count = 0ULL;
    *exact_span_count = 0U;
    *reused_page_count = 0ULL;
    *coalesced_gap_count = 0U;
    *coalesced_gap_bytes = 0U;
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

        if (exact_bytes > byte_budget ||
            span_bytes > byte_budget - exact_bytes) {
            errno = E2BIG;
            return -1;
        }
        exact_bytes += span_bytes;
    }
    original_merged = merged;
    filtered = fitsbin_payload_filter_completed_pages(
        fb,
        spans,
        merged,
        span_capacity,
        page_size,
        reuse_sequence,
        &reused_pages);
    if (filtered != SIZE_MAX) {
        merged = filtered;
        exact_bytes = 0U;
        for (i = 0U; i < merged; i++) {
            size_t span_bytes =
                (size_t)(spans[i].end - spans[i].begin);

            if (exact_bytes > byte_budget ||
                span_bytes > byte_budget - exact_bytes) {
                errno = E2BIG;
                return -1;
            }
            exact_bytes += span_bytes;
        }
    } else {
        reused_pages = 0ULL;
    }
    gap_budget = MIN(
        byte_budget - exact_bytes,
        exact_bytes /
            FITSBIN_PAYLOAD_MAPPED_COALESCE_BUDGET_DIVISOR);
    for (i = 0U; i < merged; i++) {
        if (coalesced &&
            spans[i].map_begin ==
                spans[coalesced - 1U].map_begin &&
            spans[i].map_end ==
                spans[coalesced - 1U].map_end &&
            spans[i].begin > spans[coalesced - 1U].end) {
            size_t gap = (size_t)(
                spans[i].begin - spans[coalesced - 1U].end);

            if (gap <=
                    page_size *
                        FITSBIN_PAYLOAD_MAPPED_COALESCE_GAP_PAGES &&
                gap_bytes <= gap_budget &&
                gap <= gap_budget - gap_bytes) {
                spans[coalesced - 1U].end = spans[i].end;
                gap_bytes += gap;
                gap_merges++;
                continue;
            }
        }
        spans[coalesced++] = spans[i];
    }
    for (i = 0U; i < coalesced; i++) {
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
    if (coalesced > merged || gap_merges > merged ||
        coalesced + gap_merges != merged ||
        exact_bytes > SIZE_MAX - gap_bytes ||
        aligned_bytes != exact_bytes + gap_bytes ||
        aligned_bytes > byte_budget) {
        errno = EFAULT;
        return -1;
    }
    *span_count = coalesced;
    *byte_count = aligned_bytes;
    *logical_byte_count = logical_bytes;
    *page_count = pages;
    *exact_span_count = original_merged;
    *reused_page_count = reused_pages;
    *coalesced_gap_count = gap_merges;
    *coalesced_gap_bytes = gap_bytes;
    return 0;
}

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
    size_t* gap_bytes) {
    struct stat source;
    size_t accepted = 0U;
    size_t merged = 0U;
    size_t coalesced = 0U;
    size_t exact_bytes = 0U;
    size_t gap_budget = 0U;
    size_t added_gap_bytes = 0U;
    size_t added_gap_count = 0U;
    size_t actual_bytes = 0U;
    size_t page_size;
    size_t maximum_gap;
    size_t i;

    if (!file_span_count || !byte_count || !exact_span_count ||
        !gap_count || !gap_bytes) {
        errno = EINVAL;
        return -1;
    }
    *file_span_count = 0U;
    *byte_count = 0U;
    *exact_span_count = 0U;
    *gap_count = 0U;
    *gap_bytes = 0U;
    if (!fb || !mapped || !mapped_count || !byte_budget ||
        !file_spans || mapped_count > file_span_capacity) {
        errno = mapped_count > file_span_capacity
            ? E2BIG
            : EINVAL;
        return -1;
    }
    if (fitsbin_get_open_file_stat(fb, &source)) {
        return -1;
    }
    if (source.st_size <= 0) {
        errno = EINVAL;
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
    if (page_size > SIZE_MAX /
            FITSBIN_PAYLOAD_QUEUE_COALESCE_GAP_PAGES) {
        errno = EOVERFLOW;
        return -1;
    }
    maximum_gap = page_size *
        FITSBIN_PAYLOAD_QUEUE_COALESCE_GAP_PAGES;

    for (i = 0U; i < mapped_count; i++) {
        const fitsbin_mapped_span_t* span = &mapped[i];
        fitsbin_chunk_t* chunk = NULL;
        uintptr_t data_begin;
        size_t data_map_offset;
        size_t span_begin_offset;
        size_t span_end_offset;
        off_t map_file_begin;
        off_t file_begin;
        off_t file_end;
        int chunk_index;

        if (span->map_end <= span->map_begin ||
            span->begin < span->map_begin ||
            span->end <= span->begin ||
            span->end > span->map_end) {
            errno = ERANGE;
            return -1;
        }
        for (chunk_index = 0;
             chunk_index < fitsbin_n_chunks(fb);
             chunk_index++) {
            fitsbin_chunk_t* candidate =
                fitsbin_get_chunk(fb, chunk_index);
            uintptr_t candidate_begin;
            uintptr_t candidate_size;

            if (!candidate || !candidate->map ||
                !candidate->data || !candidate->mapsize ||
                !candidate->data_file_size) {
                continue;
            }
            candidate_begin = (uintptr_t)candidate->map;
            candidate_size =
                span->map_end - span->map_begin;
            if (candidate_size > SIZE_MAX) {
                errno = EOVERFLOW;
                return -1;
            }
            if (candidate_begin == span->map_begin &&
                candidate->mapsize ==
                    (size_t)candidate_size) {
                chunk = candidate;
                break;
            }
        }
        if (!chunk) {
            errno = ERANGE;
            return -1;
        }
        data_begin = (uintptr_t)chunk->data;
        if (data_begin < span->map_begin ||
            data_begin > span->map_end) {
            errno = ERANGE;
            return -1;
        }
        if (chunk->data_file_offset < 0 ||
            chunk->data_file_offset > source.st_size ||
            chunk->data_file_size > (size_t)LLONG_MAX ||
            chunk->data_file_size >
                (size_t)(span->map_end - data_begin) ||
            (off_t)chunk->data_file_size >
                source.st_size - chunk->data_file_offset) {
            errno = ERANGE;
            return -1;
        }
        data_map_offset =
            (size_t)(data_begin - span->map_begin);
        if (chunk->data_file_offset < 0 ||
            data_map_offset > (size_t)LLONG_MAX ||
            chunk->data_file_offset < (off_t)data_map_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        map_file_begin = chunk->data_file_offset -
            (off_t)data_map_offset;
        span_begin_offset =
            (size_t)(span->begin - span->map_begin);
        span_end_offset =
            (size_t)(span->end - span->map_begin);
        if (span_begin_offset > (size_t)LLONG_MAX ||
            span_end_offset > (size_t)LLONG_MAX ||
            map_file_begin >
                (off_t)LLONG_MAX - (off_t)span_end_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        file_begin = map_file_begin +
            (off_t)span_begin_offset;
        file_end = map_file_begin +
            (off_t)span_end_offset;
        if (file_begin < 0 || file_begin >= source.st_size) {
            errno = ERANGE;
            return -1;
        }
        file_end = MIN(file_end, source.st_size);
        if (file_end <= file_begin) {
            errno = ERANGE;
            return -1;
        }
        file_spans[accepted].begin = file_begin;
        file_spans[accepted].end = file_end;
        accepted++;
    }

    qsort(file_spans,
          accepted,
          sizeof(file_spans[0]),
          fitsbin_compare_file_span);
    for (i = 0U; i < accepted; i++) {
        if (merged &&
            file_spans[i].begin <=
                file_spans[merged - 1U].end) {
            if (file_spans[i].end >
                file_spans[merged - 1U].end) {
                file_spans[merged - 1U].end =
                    file_spans[i].end;
            }
            continue;
        }
        file_spans[merged++] = file_spans[i];
    }
    for (i = 0U; i < merged; i++) {
        off_t span_length =
            file_spans[i].end - file_spans[i].begin;
        size_t span_bytes;

        if ((uintmax_t)span_length > SIZE_MAX) {
            errno = EOVERFLOW;
            return -1;
        }
        span_bytes = (size_t)span_length;

        if (exact_bytes > byte_budget ||
            span_bytes > byte_budget - exact_bytes) {
            errno = E2BIG;
            return -1;
        }
        exact_bytes += span_bytes;
    }
    gap_budget = MIN(
        byte_budget - exact_bytes,
        exact_bytes /
            FITSBIN_PAYLOAD_QUEUE_COALESCE_BUDGET_DIVISOR);
    for (i = 0U; i < merged; i++) {
        if (coalesced &&
            file_spans[i].begin >
                file_spans[coalesced - 1U].end) {
            off_t gap_offset = file_spans[i].begin -
                file_spans[coalesced - 1U].end;

            if ((uintmax_t)gap_offset <= SIZE_MAX) {
                size_t gap = (size_t)gap_offset;

                if (gap <= maximum_gap &&
                    added_gap_bytes <= gap_budget &&
                    gap <= gap_budget - added_gap_bytes) {
                    file_spans[coalesced - 1U].end =
                        file_spans[i].end;
                    added_gap_bytes += gap;
                    added_gap_count++;
                    continue;
                }
            }
        }
        file_spans[coalesced++] = file_spans[i];
    }
    for (i = 0U; i < coalesced; i++) {
        off_t span_length =
            file_spans[i].end - file_spans[i].begin;
        size_t span_bytes;

        if ((uintmax_t)span_length > SIZE_MAX) {
            errno = EOVERFLOW;
            return -1;
        }
        span_bytes = (size_t)span_length;

        if (actual_bytes > byte_budget ||
            span_bytes > byte_budget - actual_bytes) {
            errno = E2BIG;
            return -1;
        }
        actual_bytes += span_bytes;
    }
    if (coalesced > merged || added_gap_count > merged ||
        coalesced + added_gap_count != merged ||
        exact_bytes > SIZE_MAX - added_gap_bytes ||
        actual_bytes != exact_bytes + added_gap_bytes ||
        actual_bytes > byte_budget) {
        errno = EFAULT;
        return -1;
    }
    *file_span_count = coalesced;
    *byte_count = actual_bytes;
    *exact_span_count = merged;
    *gap_count = added_gap_count;
    *gap_bytes = added_gap_bytes;
    return 0;
}
