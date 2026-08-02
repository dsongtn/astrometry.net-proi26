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
     * W2+ shard and preparation threads select RANDOM for sparse payload
     * chunks. Compact topology remains NORMAL. Serial callers have no
     * thread-local advice and retain NORMAL throughout.
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

    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);
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

    if (fitsbin_payload_io_planning_is_active()) {
        errno = EDEADLK;
        return -1;
    }
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
        int capacity = fitsbin_payload_io_capacity_current();

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

    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);
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
void fitsbin_apply_mmap_advice(fitsbin_t* fb, fitsbin_chunk_t* chunk) {
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
