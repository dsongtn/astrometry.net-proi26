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

int fitsbin_compare_prepared_pread_range(
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

anbool fitsbin_prepared_pread_destinations_disjoint(
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

int fitsbin_payload_fd_get(fitsbin_t* fb) {
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

fitsbin_chunk_t* fitsbin_find_data_chunk(
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
    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);
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

/*
 * Translate any byte range inside a live file-backed mapping, including the
 * alignment and FITS-padding bytes outside the logical table payload. Exact
 * page plans describe mapped pages, so restricting this translation to
 * chunk->data would reject valid first and last pages.
 */
int fitsbin_mapped_file_offset(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    off_t* file_offset) {
    uintptr_t request;
    int i;

    if (!fb || !data || !size || !file_offset ||
        !fb->chunks) {
        errno = EINVAL;
        return -1;
    }
    request = (uintptr_t)data;
    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);
        uintptr_t map_begin;
        uintptr_t map_end;
        uintptr_t data_begin;
        size_t data_map_offset;
        size_t request_map_offset;
        off_t map_file_offset;

        if (!chunk || !chunk->map || !chunk->data ||
            !chunk->mapsize || chunk->data_file_offset < 0) {
            continue;
        }
        map_begin = (uintptr_t)chunk->map;
        if (chunk->mapsize > UINTPTR_MAX - map_begin) {
            continue;
        }
        map_end = map_begin + chunk->mapsize;
        if (request < map_begin || request >= map_end ||
            size > map_end - request) {
            continue;
        }
        data_begin = (uintptr_t)chunk->data;
        if (data_begin < map_begin || data_begin > map_end) {
            errno = ERANGE;
            return -1;
        }
        data_map_offset = (size_t)(data_begin - map_begin);
        request_map_offset = (size_t)(request - map_begin);
        if (data_map_offset > (size_t)LLONG_MAX ||
            request_map_offset > (size_t)LLONG_MAX ||
            chunk->data_file_offset <
                (off_t)data_map_offset) {
            errno = EOVERFLOW;
            return -1;
        }
        map_file_offset =
            chunk->data_file_offset - (off_t)data_map_offset;
        if (map_file_offset >
                (off_t)LLONG_MAX -
                    (off_t)request_map_offset ||
            size > (size_t)LLONG_MAX ||
            map_file_offset +
                    (off_t)request_map_offset >
                (off_t)LLONG_MAX - (off_t)size) {
            errno = EOVERFLOW;
            return -1;
        }
        *file_offset =
            map_file_offset + (off_t)request_map_offset;
        return 0;
    }
    errno = ERANGE;
    return -1;
}

void fitsbin_payload_counter_add(
    unsigned long long* counter,
    unsigned long long value) {
    if (!counter || !value || *counter == ULLONG_MAX) {
        return;
    }
    if (value > ULLONG_MAX - *counter) {
        *counter = ULLONG_MAX;
    } else {
        *counter += value;
    }
}

int fitsbin_pread_all_counted(
    int fd,
    void* destination,
    size_t size,
    off_t offset,
    unsigned long long* calls,
    unsigned long long* bytes) {
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
        ssize_t count;

        fitsbin_payload_counter_add(calls, 1ULL);
        count = pread(
            fd,
            output + complete,
            size - complete,
            offset + (off_t)complete);

        if (count > 0) {
            complete += (size_t)count;
            fitsbin_payload_counter_add(
                bytes, (unsigned long long)count);
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

int fitsbin_pread_all(
    int fd,
    void* destination,
    size_t size,
    off_t offset) {
    return fitsbin_pread_all_counted(
        fd, destination, size, offset, NULL, NULL);
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

    if (fitsbin_payload_io_planning_is_active()) {
        errno = EDEADLK;
        return -1;
    }
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
