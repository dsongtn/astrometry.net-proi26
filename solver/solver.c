/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <pthread.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdarg.h>

#include "os-features.h"
#include "ioutils.h"
#include "mathutil.h"
#include "matchobj.h"
#include "solver.h"
#include "verify.h"
#include "tic.h"
#include "solvedfile.h"
#include "fit-wcs.h"
#include "sip-utils.h"
#include "keywords.h"
#include "log.h"
#include "pquad.h"
#include "kdtree.h"
#include "quad-utils.h"
#include "errors.h"
#include "tweak2.h"
#include "astrometry/fitsbin.h"
#include "../libkd/kdtree_direct_internal.h"
#include "../libkd/kdtree_prefetch_internal.h"
#include "index_shard_internal.h"

#ifndef SOLVER_FIELD_GEOMETRY_BUDGET_BYTES
#define SOLVER_FIELD_GEOMETRY_BUDGET_BYTES \
    (64ULL * 1024ULL * 1024ULL)
#endif

#define SOLVER_CODEKD_DIRECT_MAX_POINTS \
    KDTREE_DIRECT_DSS_TASK_POINTS
#define SOLVER_CODEKD_CACHE_ENTRY_LIMIT 2048U
#define SOLVER_CODEKD_CACHE_HASH_BUCKETS 4096U
#define SOLVER_CODEKD_CACHE_NONE UINT32_MAX
#define SOLVER_CODEKD_CACHE_PAGE_PROBE_LIMIT 64U
#define SOLVER_CODEKD_READ_BATCH_MAX \
    KDTREE_DIRECT_DSS_WAVE_TASKS
#define SOLVER_CODEKD_RANGE_CACHE_BYTES \
    (8U * 1024U * 1024U)

typedef struct solver_codekd_file_segment {
    dev_t file_device;
    ino_t file_inode;
    off_t file_size;
    time_t file_mtime_seconds;
    long file_mtime_nanoseconds;
    time_t file_ctime_seconds;
    long file_ctime_nanoseconds;
    size_t page_size;
    uint32_t protected_entry;
    size_t entry_count;
    anbool fallback_reported;
    struct solver_codekd_file_segment* next;
} solver_codekd_file_segment_t;

typedef enum solver_codekd_cache_state {
    SOLVER_CODEKD_CACHE_FREE = 0,
    SOLVER_CODEKD_CACHE_LOADING = 1,
    SOLVER_CODEKD_CACHE_VALID = 2
} solver_codekd_cache_state_t;

typedef enum solver_codekd_cache_class {
    SOLVER_CODEKD_CACHE_CLASS_NONE = 0,
    SOLVER_CODEKD_CACHE_CLASS_PROTECTED = 1,
    SOLVER_CODEKD_CACHE_CLASS_EVICTABLE = 2
} solver_codekd_cache_class_t;

typedef struct solver_codekd_cache_entry {
    solver_codekd_file_segment_t* segment;
    off_t data_file_offset;
    off_t perm_file_offset;
    size_t data_size;
    size_t perm_size;
    u16* data;
    u32* perm;
    size_t data_capacity;
    size_t perm_capacity;
    uint32_t hash_bucket;
    uint32_t hash_next;
    uint32_t lru_previous;
    uint32_t lru_next;
    uint32_t free_next;
    unsigned int pin_count;
    solver_codekd_cache_state_t state;
    solver_codekd_cache_class_t cache_class;
    anbool listed;
    anbool valid;
} solver_codekd_cache_entry_t;

typedef struct solver_codekd_reader {
    const kdtree_t* tree;
    fitsbin_t* fitsbin;
    solver_codekd_file_segment_t* segments;
    solver_codekd_file_segment_t* active_segment;
    solver_codekd_cache_entry_t
        cache[SOLVER_CODEKD_CACHE_ENTRY_LIMIT];
    uint32_t hash_buckets[SOLVER_CODEKD_CACHE_HASH_BUCKETS];
    uint32_t free_head;
    uint32_t protected_head;
    uint32_t protected_tail;
    uint32_t evictable_head;
    uint32_t evictable_tail;
    kdtree_direct_dss_task_output_t
        helper_outputs[KDTREE_DIRECT_DSS_WAVE_TASKS];
    size_t cache_bytes;
    size_t total_pins;
    anbool helper_hard_failure;
} solver_codekd_reader_t;

static pthread_key_t solver_codekd_reader_key;
static pthread_once_t solver_codekd_reader_once =
    PTHREAD_ONCE_INIT;
static int solver_codekd_reader_key_status = EAGAIN;

static void solver_codekd_cache_initialize(
    solver_codekd_reader_t* reader) {
    uint32_t index;

    memset(reader->hash_buckets,
           0xff,
           sizeof(reader->hash_buckets));
    reader->free_head = 0U;
    reader->protected_head = SOLVER_CODEKD_CACHE_NONE;
    reader->protected_tail = SOLVER_CODEKD_CACHE_NONE;
    reader->evictable_head = SOLVER_CODEKD_CACHE_NONE;
    reader->evictable_tail = SOLVER_CODEKD_CACHE_NONE;
    for (index = 0U;
         index < SOLVER_CODEKD_CACHE_ENTRY_LIMIT;
         index++) {
        solver_codekd_cache_entry_t* entry =
            &reader->cache[index];

        entry->hash_bucket = SOLVER_CODEKD_CACHE_NONE;
        entry->hash_next = SOLVER_CODEKD_CACHE_NONE;
        entry->lru_previous = SOLVER_CODEKD_CACHE_NONE;
        entry->lru_next = SOLVER_CODEKD_CACHE_NONE;
        entry->free_next =
            index + 1U < SOLVER_CODEKD_CACHE_ENTRY_LIMIT
            ? index + 1U
            : SOLVER_CODEKD_CACHE_NONE;
        entry->state = SOLVER_CODEKD_CACHE_FREE;
        entry->cache_class = SOLVER_CODEKD_CACHE_CLASS_NONE;
    }
}

static void solver_codekd_reader_destroy(void* opaque) {
    solver_codekd_reader_t* reader = opaque;
    solver_codekd_file_segment_t* segment;
    size_t i;

    if (!reader) {
        return;
    }
    for (i = 0U;
         i < SOLVER_CODEKD_CACHE_ENTRY_LIMIT;
         i++) {
        if (reader->cache[i].pin_count) {
            logerr("[solver-codekd-helper] destroying pinned "
                   "cache entry slot=%zu pins=%u\n",
                   i,
                   reader->cache[i].pin_count);
        }
        free(reader->cache[i].data);
        free(reader->cache[i].perm);
    }
    segment = reader->segments;
    while (segment) {
        solver_codekd_file_segment_t* next = segment->next;

        free(segment);
        segment = next;
    }
    free(reader);
}

static void solver_codekd_reader_make_key(void) {
    solver_codekd_reader_key_status = pthread_key_create(
        &solver_codekd_reader_key,
        solver_codekd_reader_destroy);
}

static int solver_codekd_reader_prepare(void) {
    int status = pthread_once(
        &solver_codekd_reader_once,
        solver_codekd_reader_make_key);

    if (status) {
        return status;
    }
    return solver_codekd_reader_key_status;
}

static solver_codekd_reader_t* solver_codekd_reader_get(void) {
    solver_codekd_reader_t* reader;

    if (solver_codekd_reader_prepare()) {
        return NULL;
    }
    reader = pthread_getspecific(
        solver_codekd_reader_key);
    if (reader) {
        return reader;
    }
    reader = calloc(1, sizeof(*reader));
    if (!reader) {
        return NULL;
    }
    solver_codekd_cache_initialize(reader);
    if (pthread_setspecific(
            solver_codekd_reader_key,
            reader)) {
        solver_codekd_reader_destroy(reader);
        return NULL;
    }
    return reader;
}

static anbool solver_codekd_segment_matches(
    const solver_codekd_file_segment_t* segment,
    const struct stat* identity,
    time_t mtime_seconds,
    long mtime_nanoseconds,
    time_t ctime_seconds,
    long ctime_nanoseconds,
    size_t page_size) {
    return segment && identity &&
        segment->file_device == identity->st_dev &&
        segment->file_inode == identity->st_ino &&
        segment->file_size == identity->st_size &&
        segment->file_mtime_seconds == mtime_seconds &&
        segment->file_mtime_nanoseconds == mtime_nanoseconds &&
        segment->file_ctime_seconds == ctime_seconds &&
        segment->file_ctime_nanoseconds == ctime_nanoseconds &&
        segment->page_size == page_size;
}

static int solver_codekd_reader_bind(
    solver_codekd_reader_t* reader,
    const kdtree_t* tree,
    fitsbin_t* fb) {
    const struct stat* identity;
    solver_codekd_file_segment_t* segment;
    time_t mtime_seconds;
    long mtime_nanoseconds;
    time_t ctime_seconds;
    long ctime_nanoseconds;
    size_t page_size;

    if (!reader || !tree || !fb || !fb->open_file_stat_valid) {
        errno = EINVAL;
        return -1;
    }
    identity = &fb->open_file_stat;
    page_size = fb->mmap_page_size;
    if (!page_size) {
        long detected = sysconf(_SC_PAGESIZE);

        if (detected <= 0) {
            errno = EINVAL;
            return -1;
        }
        page_size = (size_t)detected;
    }

    fitsbin_stat_times(
        identity,
        &mtime_seconds,
        &mtime_nanoseconds,
        &ctime_seconds,
        &ctime_nanoseconds);
    segment = reader->active_segment;
    if (!solver_codekd_segment_matches(
            segment,
            identity,
            mtime_seconds,
            mtime_nanoseconds,
            ctime_seconds,
            ctime_nanoseconds,
            page_size)) {
        for (segment = reader->segments;
             segment;
             segment = segment->next) {
            if (solver_codekd_segment_matches(
                    segment,
                    identity,
                    mtime_seconds,
                    mtime_nanoseconds,
                    ctime_seconds,
                    ctime_nanoseconds,
                    page_size)) {
                break;
            }
        }
    }
    if (!segment) {
        segment = calloc(1, sizeof(*segment));
        if (!segment) {
            return -1;
        }
        segment->file_device = identity->st_dev;
        segment->file_inode = identity->st_ino;
        segment->file_size = identity->st_size;
        segment->file_mtime_seconds = mtime_seconds;
        segment->file_mtime_nanoseconds = mtime_nanoseconds;
        segment->file_ctime_seconds = ctime_seconds;
        segment->file_ctime_nanoseconds = ctime_nanoseconds;
        segment->page_size = page_size;
        segment->protected_entry = SOLVER_CODEKD_CACHE_NONE;
        segment->next = reader->segments;
        reader->segments = segment;
    }

    /*
     * These pointers are borrowed only for the active direct search. Never
     * inspect a pointer retained from an index whose lifetime has ended.
     */
    reader->tree = tree;
    reader->fitsbin = fb;
    reader->active_segment = segment;
    /*
     * Only scalar file identity survives an index close. The cached bytes are
     * private copies, so a deeper affinity pass can reuse them without
     * retaining an index mapping, fitsbin pointer, or other borrowed state.
     */
    return 0;
}

static void solver_codekd_reader_unbind(
    solver_codekd_reader_t* reader) {
    if (!reader) {
        return;
    }
    if (reader->total_pins) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
    }
    reader->tree = NULL;
    reader->fitsbin = NULL;
}

static uint32_t solver_codekd_cache_index(
    const solver_codekd_reader_t* reader,
    const solver_codekd_cache_entry_t* entry) {
    return (uint32_t)(entry - reader->cache);
}

static uint64_t solver_codekd_cache_mix(
    uint64_t hash,
    uint64_t value) {
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return hash ^ (value + UINT64_C(0x9e3779b97f4a7c15) +
                   (hash << 6) + (hash >> 2));
}

static uint32_t solver_codekd_cache_bucket(
    const solver_codekd_file_segment_t* segment,
    uint64_t data_page) {
    uint64_t hash = UINT64_C(0xcbf29ce484222325);

    hash = solver_codekd_cache_mix(
        hash, (uint64_t)(uintptr_t)segment);
    hash = solver_codekd_cache_mix(hash, data_page);
    return (uint32_t)hash &
        (SOLVER_CODEKD_CACHE_HASH_BUCKETS - 1U);
}

static anbool solver_codekd_cache_range_contains(
    off_t outer_file_offset,
    size_t outer_size,
    off_t inner_file_offset,
    size_t inner_size) {
    uint64_t outer;
    uint64_t inner;
    uint64_t delta;

    if (outer_file_offset < 0 || inner_file_offset < 0 ||
        inner_file_offset < outer_file_offset) {
        return FALSE;
    }
    outer = (uint64_t)outer_file_offset;
    inner = (uint64_t)inner_file_offset;
    delta = inner - outer;
    return delta <= SIZE_MAX &&
        (size_t)delta <= outer_size &&
        inner_size <= outer_size - (size_t)delta;
}

static anbool solver_codekd_cache_entry_contains(
    const solver_codekd_cache_entry_t* entry,
    const solver_codekd_file_segment_t* segment,
    off_t data_file_offset,
    size_t data_size,
    off_t perm_file_offset,
    size_t perm_size) {
    return entry && entry->valid &&
        entry->state == SOLVER_CODEKD_CACHE_VALID &&
        entry->segment == segment &&
        solver_codekd_cache_range_contains(
            entry->data_file_offset,
            entry->data_size,
            data_file_offset,
            data_size) &&
        ((!perm_size && !entry->perm_size) ||
         (perm_size && entry->perm_size &&
          solver_codekd_cache_range_contains(
              entry->perm_file_offset,
              entry->perm_size,
              perm_file_offset,
              perm_size)));
}

static size_t solver_codekd_cache_page_probes(
    const solver_codekd_reader_t* reader,
    const solver_codekd_file_segment_t* segment) {
    size_t maximum_bytes;
    size_t pages;

    if (!reader || !reader->tree || !segment ||
        !segment->page_size || reader->tree->ndim <= 0 ||
        KDTREE_DIRECT_DSS_TASK_POINTS >
            SIZE_MAX / (size_t)reader->tree->ndim ||
        KDTREE_DIRECT_DSS_TASK_POINTS *
                (size_t)reader->tree->ndim >
            SIZE_MAX / sizeof(u16)) {
        return SOLVER_CODEKD_CACHE_PAGE_PROBE_LIMIT;
    }
    maximum_bytes = KDTREE_DIRECT_DSS_TASK_POINTS *
        (size_t)reader->tree->ndim * sizeof(u16);
    pages = maximum_bytes / segment->page_size;
    if (maximum_bytes % segment->page_size) {
        pages++;
    }
    if (pages > SOLVER_CODEKD_CACHE_PAGE_PROBE_LIMIT - 2U) {
        return SOLVER_CODEKD_CACHE_PAGE_PROBE_LIMIT;
    }
    return pages + 2U;
}

static solver_codekd_cache_entry_t* solver_codekd_cache_find(
    solver_codekd_reader_t* reader,
    solver_codekd_file_segment_t* segment,
    off_t data_file_offset,
    size_t data_size,
    off_t perm_file_offset,
    size_t perm_size) {
    uint64_t data_page;
    size_t probes;
    size_t probe;

    if (!reader || !segment || !segment->page_size ||
        data_file_offset < 0 || !data_size) {
        if (reader) {
            reader->helper_hard_failure = TRUE;
        }
        errno = EPROTO;
        return NULL;
    }
    data_page = (uint64_t)data_file_offset /
        (uint64_t)segment->page_size;
    probes = solver_codekd_cache_page_probes(reader, segment);
    for (probe = 0U;
         probe < probes && (uint64_t)probe <= data_page;
         probe++) {
        uint32_t bucket = solver_codekd_cache_bucket(
            segment, data_page - (uint64_t)probe);
        uint32_t index = reader->hash_buckets[bucket];
        size_t visited = 0U;

        while (index != SOLVER_CODEKD_CACHE_NONE) {
            solver_codekd_cache_entry_t* entry;

            if (index >= SOLVER_CODEKD_CACHE_ENTRY_LIMIT ||
                visited++ >= SOLVER_CODEKD_CACHE_ENTRY_LIMIT) {
                reader->helper_hard_failure = TRUE;
                errno = EPROTO;
                return NULL;
            }
            entry = &reader->cache[index];
            if (entry->hash_bucket != bucket) {
                reader->helper_hard_failure = TRUE;
                errno = EPROTO;
                return NULL;
            }
            if (solver_codekd_cache_entry_contains(
                    entry,
                    segment,
                    data_file_offset,
                    data_size,
                    perm_file_offset,
                    perm_size)) {
                return entry;
            }
            index = entry->hash_next;
        }
    }
    return NULL;
}

static void solver_codekd_cache_list_remove(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    uint32_t* head;
    uint32_t* tail;
    uint32_t index;

    if (!entry->listed) {
        return;
    }
    if (entry->cache_class ==
            SOLVER_CODEKD_CACHE_CLASS_PROTECTED) {
        head = &reader->protected_head;
        tail = &reader->protected_tail;
    } else if (entry->cache_class ==
                   SOLVER_CODEKD_CACHE_CLASS_EVICTABLE) {
        head = &reader->evictable_head;
        tail = &reader->evictable_tail;
    } else {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    index = solver_codekd_cache_index(reader, entry);
    if (entry->lru_previous == SOLVER_CODEKD_CACHE_NONE) {
        *head = entry->lru_next;
    } else {
        reader->cache[entry->lru_previous].lru_next =
            entry->lru_next;
    }
    if (entry->lru_next == SOLVER_CODEKD_CACHE_NONE) {
        *tail = entry->lru_previous;
    } else {
        reader->cache[entry->lru_next].lru_previous =
            entry->lru_previous;
    }
    entry->lru_previous = SOLVER_CODEKD_CACHE_NONE;
    entry->lru_next = SOLVER_CODEKD_CACHE_NONE;
    entry->listed = FALSE;
    (void)index;
}

static void solver_codekd_cache_list_push_mru(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    uint32_t* head;
    uint32_t* tail;
    uint32_t index = solver_codekd_cache_index(reader, entry);

    if (entry->listed || entry->pin_count || !entry->valid) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    if (entry->cache_class ==
            SOLVER_CODEKD_CACHE_CLASS_PROTECTED) {
        head = &reader->protected_head;
        tail = &reader->protected_tail;
    } else if (entry->cache_class ==
                   SOLVER_CODEKD_CACHE_CLASS_EVICTABLE) {
        head = &reader->evictable_head;
        tail = &reader->evictable_tail;
    } else {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    entry->lru_previous = SOLVER_CODEKD_CACHE_NONE;
    entry->lru_next = *head;
    if (*head != SOLVER_CODEKD_CACHE_NONE) {
        reader->cache[*head].lru_previous = index;
    } else {
        *tail = index;
    }
    *head = index;
    entry->listed = TRUE;
}

static void solver_codekd_cache_promote(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    solver_codekd_file_segment_t* segment = entry->segment;
    uint32_t index = solver_codekd_cache_index(reader, entry);
    uint32_t previous = segment->protected_entry;

    if (previous != SOLVER_CODEKD_CACHE_NONE &&
        previous != index) {
        solver_codekd_cache_entry_t* old =
            &reader->cache[previous];

        if (old->listed) {
            solver_codekd_cache_list_remove(reader, old);
        }
        old->cache_class = SOLVER_CODEKD_CACHE_CLASS_EVICTABLE;
        if (!old->pin_count) {
            solver_codekd_cache_list_push_mru(reader, old);
        }
    }
    if (entry->listed) {
        solver_codekd_cache_list_remove(reader, entry);
    }
    entry->cache_class = SOLVER_CODEKD_CACHE_CLASS_PROTECTED;
    segment->protected_entry = index;
    if (!entry->pin_count) {
        solver_codekd_cache_list_push_mru(reader, entry);
    }
}

static void solver_codekd_cache_hash_remove(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    uint32_t bucket = entry->hash_bucket;
    uint32_t* link;
    uint32_t index = solver_codekd_cache_index(reader, entry);
    size_t visited = 0U;

    if (bucket >= SOLVER_CODEKD_CACHE_HASH_BUCKETS) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    link = &reader->hash_buckets[bucket];
    while (*link != SOLVER_CODEKD_CACHE_NONE &&
           visited++ < SOLVER_CODEKD_CACHE_ENTRY_LIMIT) {
        solver_codekd_cache_entry_t* current;

        if (*link >= SOLVER_CODEKD_CACHE_ENTRY_LIMIT) {
            break;
        }
        current = &reader->cache[*link];
        if (*link == index) {
            *link = current->hash_next;
            entry->hash_bucket = SOLVER_CODEKD_CACHE_NONE;
            entry->hash_next = SOLVER_CODEKD_CACHE_NONE;
            return;
        }
        link = &current->hash_next;
    }
    reader->helper_hard_failure = TRUE;
    errno = EPROTO;
}

static void solver_codekd_cache_detach(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    uint32_t index;

    if (!entry->valid ||
        entry->state != SOLVER_CODEKD_CACHE_VALID) {
        return;
    }
    if (entry->pin_count) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    index = solver_codekd_cache_index(reader, entry);
    solver_codekd_cache_hash_remove(reader, entry);
    if (entry->listed) {
        solver_codekd_cache_list_remove(reader, entry);
    }
    if (entry->segment->protected_entry == index) {
        entry->segment->protected_entry = SOLVER_CODEKD_CACHE_NONE;
    }
    if (entry->segment->entry_count) {
        entry->segment->entry_count--;
    } else {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
    }
    entry->segment = NULL;
    entry->data_file_offset = 0;
    entry->perm_file_offset = 0;
    entry->data_size = 0U;
    entry->perm_size = 0U;
    entry->state = SOLVER_CODEKD_CACHE_LOADING;
    entry->cache_class = SOLVER_CODEKD_CACHE_CLASS_NONE;
    entry->valid = FALSE;
    if (reader->fitsbin) {
        __atomic_add_fetch(
            &reader->fitsbin->payload_cache_evictions,
            1ULL,
            __ATOMIC_RELAXED);
    }
}

static void solver_codekd_cache_release_storage(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    size_t retained =
        entry->data_capacity + entry->perm_capacity;

    reader->cache_bytes -= MIN(reader->cache_bytes, retained);
    free(entry->data);
    free(entry->perm);
    entry->data = NULL;
    entry->perm = NULL;
    entry->data_capacity = 0U;
    entry->perm_capacity = 0U;
}

static void solver_codekd_cache_free_push(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry) {
    uint32_t index = solver_codekd_cache_index(reader, entry);

    entry->segment = NULL;
    entry->hash_bucket = SOLVER_CODEKD_CACHE_NONE;
    entry->hash_next = SOLVER_CODEKD_CACHE_NONE;
    entry->lru_previous = SOLVER_CODEKD_CACHE_NONE;
    entry->lru_next = SOLVER_CODEKD_CACHE_NONE;
    entry->pin_count = 0U;
    entry->state = SOLVER_CODEKD_CACHE_FREE;
    entry->cache_class = SOLVER_CODEKD_CACHE_CLASS_NONE;
    entry->listed = FALSE;
    entry->valid = FALSE;
    entry->free_next = reader->free_head;
    reader->free_head = index;
}

static solver_codekd_cache_entry_t* solver_codekd_cache_take(
    solver_codekd_reader_t* reader) {
    solver_codekd_cache_entry_t* entry;
    uint32_t index;

    if (reader->free_head != SOLVER_CODEKD_CACHE_NONE) {
        index = reader->free_head;
        entry = &reader->cache[index];
        reader->free_head = entry->free_next;
        entry->free_next = SOLVER_CODEKD_CACHE_NONE;
        entry->state = SOLVER_CODEKD_CACHE_LOADING;
        return entry;
    }
    index = reader->evictable_tail != SOLVER_CODEKD_CACHE_NONE
        ? reader->evictable_tail
        : reader->protected_tail;
    if (index == SOLVER_CODEKD_CACHE_NONE) {
        errno = ENOMEM;
        return NULL;
    }
    entry = &reader->cache[index];
    solver_codekd_cache_detach(reader, entry);
    return entry;
}

static int solver_codekd_cache_reclaim_storage(
    solver_codekd_reader_t* reader) {
    uint32_t index =
        reader->evictable_tail != SOLVER_CODEKD_CACHE_NONE
        ? reader->evictable_tail
        : reader->protected_tail;
    solver_codekd_cache_entry_t* entry;

    if (index == SOLVER_CODEKD_CACHE_NONE) {
        errno = ENOMEM;
        return -1;
    }
    entry = &reader->cache[index];
    solver_codekd_cache_detach(reader, entry);
    solver_codekd_cache_release_storage(reader, entry);
    solver_codekd_cache_free_push(reader, entry);
    return 0;
}

static int solver_codekd_cache_reserve(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry,
    size_t data_bytes,
    size_t perm_bytes) {
    size_t data_growth;
    size_t perm_growth;
    size_t growth;

    if (!reader || !entry || entry->pin_count) {
        if (reader) {
            reader->helper_hard_failure = TRUE;
        }
        errno = EPROTO;
        return -1;
    }
    data_growth =
        data_bytes > entry->data_capacity
        ? data_bytes - entry->data_capacity
        : 0U;
    perm_growth =
        perm_bytes > entry->perm_capacity
        ? perm_bytes - entry->perm_capacity
        : 0U;

    if (data_growth > SIZE_MAX - perm_growth) {
        errno = EOVERFLOW;
        return -1;
    }
    growth = data_growth + perm_growth;
    if (growth > SOLVER_CODEKD_RANGE_CACHE_BYTES) {
        errno = ENOMEM;
        return -1;
    }
    while (reader->cache_bytes >
           SOLVER_CODEKD_RANGE_CACHE_BYTES - growth) {
        if (solver_codekd_cache_reclaim_storage(reader)) {
            return -1;
        }
    }
    if (data_bytes > entry->data_capacity) {
        u16* grown = realloc(entry->data, data_bytes);

        if (!grown) {
            return -1;
        }
        reader->cache_bytes +=
            data_bytes - entry->data_capacity;
        entry->data = grown;
        entry->data_capacity = data_bytes;
        __atomic_add_fetch(
            &reader->fitsbin->payload_cache_allocations,
            1ULL,
            __ATOMIC_RELAXED);
    }
    if (perm_bytes > entry->perm_capacity) {
        u32* grown = realloc(entry->perm, perm_bytes);

        if (!grown) {
            return -1;
        }
        reader->cache_bytes +=
            perm_bytes - entry->perm_capacity;
        entry->perm = grown;
        entry->perm_capacity = perm_bytes;
        __atomic_add_fetch(
            &reader->fitsbin->payload_cache_allocations,
            1ULL,
            __ATOMIC_RELAXED);
    }
    return 0;
}

static void solver_codekd_cache_publish(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry,
    solver_codekd_file_segment_t* segment,
    off_t data_file_offset,
    size_t data_size,
    off_t perm_file_offset,
    size_t perm_size) {
    uint64_t data_page;
    uint32_t bucket;

    if (!segment || !segment->page_size || data_file_offset < 0) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    data_page = (uint64_t)data_file_offset /
        (uint64_t)segment->page_size;
    bucket = solver_codekd_cache_bucket(segment, data_page);
    entry->data_file_offset = data_file_offset;
    entry->perm_file_offset = perm_file_offset;
    entry->data_size = data_size;
    entry->perm_size = perm_size;
    entry->segment = segment;
    entry->state = SOLVER_CODEKD_CACHE_VALID;
    entry->valid = TRUE;
    entry->hash_bucket = bucket;
    entry->hash_next = reader->hash_buckets[bucket];
    reader->hash_buckets[bucket] =
        solver_codekd_cache_index(reader, entry);
    segment->entry_count++;
    solver_codekd_cache_promote(reader, entry);
}

#if defined(TESTING)
int solver_test_codekd_cache_churn(void) {
    solver_codekd_reader_t* reader;
    solver_codekd_file_segment_t* segment;
    fitsbin_t fitsbin;
    size_t i;
    int rc = -1;

    reader = calloc(1, sizeof(*reader));
    if (!reader) {
        return -1;
    }
    segment = calloc(1, sizeof(*segment));
    if (!segment) {
        free(reader);
        return -1;
    }
    memset(&fitsbin, 0, sizeof(fitsbin));
    solver_codekd_cache_initialize(reader);
    segment->page_size = 4096U;
    segment->protected_entry = SOLVER_CODEKD_CACHE_NONE;
    reader->segments = segment;
    reader->active_segment = segment;
    reader->fitsbin = &fitsbin;
    for (i = 0U;
         i < SOLVER_CODEKD_CACHE_ENTRY_LIMIT + 64U;
         i++) {
        solver_codekd_cache_entry_t* entry =
            solver_codekd_cache_take(reader);

        if (!entry) {
            goto cleanup;
        }
        if (solver_codekd_cache_reserve(
                reader,
                entry,
                4U * 1024U,
                4U * 1024U)) {
            goto cleanup;
        }
        solver_codekd_cache_publish(
            reader,
            entry,
            segment,
            (off_t)(i * 8192U),
            4U * 1024U,
            (off_t)(i * 8192U + 4096U),
            4U * 1024U);
        if (reader->cache_bytes >
            SOLVER_CODEKD_RANGE_CACHE_BYTES) {
            goto cleanup;
        }
    }
    i--;
    if (!solver_codekd_cache_find(
            reader,
            segment,
            (off_t)(i * 8192U),
            4U * 1024U,
            (off_t)(i * 8192U + 4096U),
            4U * 1024U) ||
        segment->protected_entry == SOLVER_CODEKD_CACHE_NONE ||
        !fitsbin.payload_cache_evictions) {
        goto cleanup;
    }
    rc = 0;

cleanup:
    solver_codekd_reader_destroy(reader);
    return rc;
}
#endif

static int solver_codekd_pin_range(
    solver_codekd_reader_t* reader,
    solver_codekd_cache_entry_t* entry,
    kdtree_direct_range_t* range) {
    if (!reader || !entry || !range || !entry->valid ||
        entry->state != SOLVER_CODEKD_CACHE_VALID ||
        entry->pin_count == UINT_MAX ||
        reader->total_pins == SIZE_MAX) {
        if (reader) {
            reader->helper_hard_failure = TRUE;
        }
        errno = EPROTO;
        return -1;
    }
    if (!entry->pin_count && entry->listed) {
        solver_codekd_cache_list_remove(reader, entry);
    }
    entry->pin_count++;
    reader->total_pins++;
    range->lease = entry;
    return 0;
}

typedef struct solver_codekd_range_plan {
    const void* data_cover_source;
    const void* perm_cover_source;
    off_t data_cover_file_offset;
    off_t perm_cover_file_offset;
    size_t data_cover_size;
    size_t perm_cover_size;
    size_t data_exact_offset;
    size_t perm_exact_offset;
    size_t data_logical_size;
    size_t perm_logical_size;
    solver_codekd_cache_entry_t* entry;
    anbool owns_loading;
} solver_codekd_range_plan_t;

static int solver_codekd_prepare_range_plan(
    solver_codekd_reader_t* reader,
    const kdtree_direct_range_request_t* request,
    solver_codekd_range_plan_t* plan) {
    const void* exact_data_source;
    const void* exact_perm_source = NULL;
    size_t first_data_element;
    size_t data_elements;

    memset(plan, 0, sizeof(*plan));
    if (!reader || !reader->tree || !reader->fitsbin ||
        !reader->active_segment || !request ||
        request->first < 0 || request->count <= 0 ||
        reader->tree->ndata <= 0 ||
        request->count > reader->tree->ndata ||
        request->first > reader->tree->ndata - request->count ||
        reader->tree->ndim <= 0 ||
        (size_t)request->count >
            SIZE_MAX / (size_t)reader->tree->ndim ||
        (size_t)request->first >
            SIZE_MAX / (size_t)reader->tree->ndim) {
        errno = EINVAL;
        return -1;
    }
    data_elements =
        (size_t)request->count * (size_t)reader->tree->ndim;
    first_data_element =
        (size_t)request->first * (size_t)reader->tree->ndim;
    if (data_elements > SIZE_MAX / sizeof(u16)) {
        errno = EOVERFLOW;
        return -1;
    }
    plan->data_logical_size = data_elements * sizeof(u16);
    exact_data_source =
        reader->tree->data.s + first_data_element;
    if (reader->tree->perm) {
        if ((size_t)request->count > SIZE_MAX / sizeof(u32)) {
            errno = EOVERFLOW;
            return -1;
        }
        plan->perm_logical_size =
            (size_t)request->count * sizeof(u32);
        exact_perm_source =
            reader->tree->perm + request->first;
    }
    if (fitsbin_mapped_range_page_cover(
            reader->fitsbin,
            exact_data_source,
            plan->data_logical_size,
            &plan->data_cover_source,
            &plan->data_cover_size,
            &plan->data_cover_file_offset,
            &plan->data_exact_offset)) {
        return -1;
    }
    if (plan->data_exact_offset % sizeof(u16) ||
        plan->data_exact_offset > plan->data_cover_size ||
        plan->data_logical_size >
            plan->data_cover_size - plan->data_exact_offset) {
        errno = EINVAL;
        return -1;
    }
    if (exact_perm_source) {
        if (fitsbin_mapped_range_page_cover(
                reader->fitsbin,
                exact_perm_source,
                plan->perm_logical_size,
                &plan->perm_cover_source,
                &plan->perm_cover_size,
                &plan->perm_cover_file_offset,
                &plan->perm_exact_offset)) {
            return -1;
        }
        if (plan->perm_exact_offset % sizeof(u32) ||
            plan->perm_exact_offset > plan->perm_cover_size ||
            plan->perm_logical_size >
                plan->perm_cover_size - plan->perm_exact_offset) {
            errno = EINVAL;
            return -1;
        }
    }
    if (plan->data_cover_size >
            SIZE_MAX - plan->perm_cover_size ||
        plan->data_cover_size + plan->perm_cover_size >
            SOLVER_CODEKD_RANGE_CACHE_BYTES) {
        errno = E2BIG;
        return -1;
    }
    return 0;
}

static anbool solver_codekd_plan_contains(
    const solver_codekd_range_plan_t* outer,
    const solver_codekd_range_plan_t* inner) {
    return solver_codekd_cache_range_contains(
               outer->data_cover_file_offset,
               outer->data_cover_size,
               inner->data_cover_file_offset,
               inner->data_cover_size) &&
        ((!inner->perm_cover_size && !outer->perm_cover_size) ||
         (inner->perm_cover_size && outer->perm_cover_size &&
          solver_codekd_cache_range_contains(
              outer->perm_cover_file_offset,
              outer->perm_cover_size,
              inner->perm_cover_file_offset,
              inner->perm_cover_size)));
}

static int solver_codekd_plan_entry_offsets(
    const solver_codekd_cache_entry_t* entry,
    const solver_codekd_range_plan_t* plan,
    size_t* data_offset,
    size_t* perm_offset) {
    uint64_t data_cover_delta;

    if (!entry || !plan || !data_offset || !perm_offset ||
        !solver_codekd_cache_entry_contains(
            entry,
            entry->segment,
            plan->data_cover_file_offset,
            plan->data_cover_size,
            plan->perm_cover_file_offset,
            plan->perm_cover_size)) {
        errno = EPROTO;
        return -1;
    }
    data_cover_delta =
        (uint64_t)plan->data_cover_file_offset -
        (uint64_t)entry->data_file_offset;
    if (data_cover_delta > SIZE_MAX ||
        plan->data_exact_offset >
            SIZE_MAX - (size_t)data_cover_delta) {
        errno = EOVERFLOW;
        return -1;
    }
    *data_offset =
        (size_t)data_cover_delta + plan->data_exact_offset;
    if (*data_offset % sizeof(u16) ||
        *data_offset > entry->data_size ||
        plan->data_logical_size >
            entry->data_size - *data_offset) {
        errno = EPROTO;
        return -1;
    }
    *perm_offset = 0U;
    if (plan->perm_cover_size) {
        uint64_t perm_cover_delta =
            (uint64_t)plan->perm_cover_file_offset -
            (uint64_t)entry->perm_file_offset;

        if (perm_cover_delta > SIZE_MAX ||
            plan->perm_exact_offset >
                SIZE_MAX - (size_t)perm_cover_delta) {
            errno = EOVERFLOW;
            return -1;
        }
        *perm_offset =
            (size_t)perm_cover_delta + plan->perm_exact_offset;
        if (*perm_offset % sizeof(u32) ||
            *perm_offset > entry->perm_size ||
            plan->perm_logical_size >
                entry->perm_size - *perm_offset) {
            errno = EPROTO;
            return -1;
        }
    }
    return 0;
}

static int solver_codekd_assign_range(
    solver_codekd_reader_t* reader,
    const solver_codekd_range_plan_t* plan,
    kdtree_direct_range_t* range) {
    size_t data_offset;
    size_t perm_offset;

    if (solver_codekd_plan_entry_offsets(
            plan->entry,
            plan,
            &data_offset,
            &perm_offset)) {
        reader->helper_hard_failure = TRUE;
        return -1;
    }
    range->data = (const u16*)(
        (const unsigned char*)plan->entry->data + data_offset);
    range->perm = plan->perm_cover_size
        ? (const u32*)(
              (const unsigned char*)plan->entry->perm + perm_offset)
        : NULL;
    return solver_codekd_pin_range(
        reader, plan->entry, range);
}

static void solver_codekd_rollback_loading_plans(
    solver_codekd_reader_t* reader,
    solver_codekd_range_plan_t* plans,
    size_t nplans) {
    size_t plan_index;

    for (plan_index = 0U; plan_index < nplans; plan_index++) {
        solver_codekd_cache_entry_t* entry =
            plans[plan_index].entry;

        if (!plans[plan_index].owns_loading || !entry ||
            entry->state != SOLVER_CODEKD_CACHE_LOADING) {
            continue;
        }
        solver_codekd_cache_release_storage(reader, entry);
        solver_codekd_cache_free_push(reader, entry);
    }
}

static int solver_codekd_read_leased_ranges(
    void* opaque,
    const kdtree_direct_range_request_t* requests,
    size_t nrequests,
    kdtree_direct_range_t* ranges) {
    solver_codekd_reader_t* reader = opaque;
    solver_codekd_range_plan_t plans[SOLVER_CODEKD_READ_BATCH_MAX];
    fitsbin_pread_range_t
        io_ranges[2U * SOLVER_CODEKD_READ_BATCH_MAX];
    size_t io_count = 0U;
    size_t request_index;
    int saved_errno;

    if (!reader || !requests || !ranges || !nrequests ||
        nrequests > SOLVER_CODEKD_READ_BATCH_MAX ||
        2U * nrequests > FITSBIN_PREAD_RANGE_LIMIT) {
        errno = EINVAL;
        return -1;
    }
    memset(ranges, 0, nrequests * sizeof(*ranges));
    memset(plans, 0, sizeof(plans));
    memset(io_ranges, 0, sizeof(io_ranges));

    for (request_index = 0U;
         request_index < nrequests;
         request_index++) {
        solver_codekd_range_plan_t* plan = &plans[request_index];
        solver_codekd_cache_entry_t* entry;
        size_t previous;

        if (solver_codekd_prepare_range_plan(
                reader,
                &requests[request_index],
                plan)) {
            goto fail;
        }
        entry = solver_codekd_cache_find(
            reader,
            reader->active_segment,
            plan->data_cover_file_offset,
            plan->data_cover_size,
            plan->perm_cover_file_offset,
            plan->perm_cover_size);
        if (reader->helper_hard_failure) {
            errno = EPROTO;
            goto fail;
        }
        if (entry) {
            plan->entry = entry;
            solver_codekd_cache_promote(reader, entry);
            __atomic_add_fetch(
                &reader->fitsbin->payload_cache_hits,
                1ULL,
                __ATOMIC_RELAXED);
            if (solver_codekd_assign_range(
                    reader, plan, &ranges[request_index])) {
                goto fail;
            }
            continue;
        }
        __atomic_add_fetch(
            &reader->fitsbin->payload_cache_misses,
            1ULL,
            __ATOMIC_RELAXED);

        for (previous = 0U; previous < request_index; previous++) {
            if (plans[previous].owns_loading &&
                plans[previous].entry &&
                plans[previous].entry->state ==
                    SOLVER_CODEKD_CACHE_LOADING &&
                solver_codekd_plan_contains(
                    &plans[previous], plan)) {
                entry = plans[previous].entry;
                break;
            }
        }
        if (previous < request_index) {
            plan->entry = entry;
            continue;
        }

        entry = solver_codekd_cache_take(reader);
        if (!entry) {
            goto fail;
        }
        plan->entry = entry;
        plan->owns_loading = TRUE;
        if (solver_codekd_cache_reserve(
                reader,
                entry,
                plan->data_cover_size,
                plan->perm_cover_size)) {
            goto fail;
        }
        io_ranges[io_count].data = plan->data_cover_source;
        io_ranges[io_count].size = plan->data_cover_size;
        io_ranges[io_count].logical_size =
            plan->data_logical_size;
        io_ranges[io_count].destination = entry->data;
        io_count++;
        if (plan->perm_cover_size) {
            io_ranges[io_count].data = plan->perm_cover_source;
            io_ranges[io_count].size = plan->perm_cover_size;
            io_ranges[io_count].logical_size =
                plan->perm_logical_size;
            io_ranges[io_count].destination = entry->perm;
            io_count++;
        }
    }

    if (io_count && fitsbin_pread_mapped_ranges(
            reader->fitsbin,
            io_ranges,
            io_count)) {
        goto fail;
    }
    for (request_index = 0U;
         request_index < nrequests;
         request_index++) {
        solver_codekd_range_plan_t* plan = &plans[request_index];

        if (!plan->owns_loading) {
            continue;
        }
        solver_codekd_cache_publish(
            reader,
            plan->entry,
            reader->active_segment,
            plan->data_cover_file_offset,
            plan->data_cover_size,
            plan->perm_cover_file_offset,
            plan->perm_cover_size);
        if (reader->helper_hard_failure) {
            errno = EPROTO;
            goto fail;
        }
    }
    for (request_index = 0U;
         request_index < nrequests;
         request_index++) {
        solver_codekd_range_plan_t* plan = &plans[request_index];

        if (ranges[request_index].lease) {
            continue;
        }
        if (solver_codekd_assign_range(
                reader, plan, &ranges[request_index])) {
            goto fail;
        }
    }
    return 0;

fail:
    saved_errno = errno ? errno : EIO;
    solver_codekd_rollback_loading_plans(
        reader, plans, nrequests);
    errno = saved_errno;
    return -1;
}

static void solver_codekd_release_range(
    void* opaque,
    void* lease) {
    solver_codekd_reader_t* reader = opaque;
    solver_codekd_cache_entry_t* entry;
    uintptr_t cache_begin;
    uintptr_t lease_address;
    size_t offset;

    if (!reader || !lease) {
        if (reader) {
            reader->helper_hard_failure = TRUE;
        }
        errno = EPROTO;
        return;
    }
    cache_begin = (uintptr_t)&reader->cache[0];
    lease_address = (uintptr_t)lease;
    if (lease_address < cache_begin ||
        lease_address - cache_begin >= sizeof(reader->cache)) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    offset = (size_t)(lease_address - cache_begin);
    if (offset % sizeof(reader->cache[0])) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    entry = &reader->cache[
        offset / sizeof(reader->cache[0])];
    if (entry != lease || !entry->pin_count ||
        !reader->total_pins || !entry->valid ||
        entry->state != SOLVER_CODEKD_CACHE_VALID) {
        reader->helper_hard_failure = TRUE;
        errno = EPROTO;
        return;
    }
    entry->pin_count--;
    reader->total_pins--;
    if (entry->pin_count) {
        return;
    }
    if (entry->segment->protected_entry ==
        SOLVER_CODEKD_CACHE_NONE) {
        solver_codekd_cache_promote(reader, entry);
    } else {
        solver_codekd_cache_list_push_mru(reader, entry);
    }
}

static index_shard_helper_task_status_t
solver_codekd_helper_execute(
    const void* input,
    size_t input_bytes,
    void* output,
    size_t output_bytes) {
    if (!input ||
        input_bytes != sizeof(kdtree_direct_dss_task_input_t) ||
        !output ||
        output_bytes != sizeof(kdtree_direct_dss_task_output_t)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    return kdtree_direct_dss_task_execute(input, output)
        ? INDEX_SHARD_HELPER_TASK_ERROR
        : INDEX_SHARD_HELPER_TASK_OK;
}

static const index_shard_helper_ops_t solver_codekd_helper_ops = {
    "codekd-leaf",
    solver_codekd_helper_execute
};

static size_t solver_codekd_helper_available(void* opaque) {
    if (!opaque || !index_shard_worker_context_active()) {
        return 0U;
    }
    return index_shard_helper_available_workers();
}

static int solver_codekd_helper_run(
    void* opaque,
    const kdtree_direct_dss_task_input_t* inputs,
    size_t ntasks,
    const kdtree_direct_dss_task_output_t** outputs) {
    solver_codekd_reader_t* reader = opaque;
    index_shard_helper_task_t tasks[KDTREE_DIRECT_DSS_WAVE_TASKS];
    index_shard_helper_run_status_t run_status;
    size_t task_index;

    if (outputs) {
        *outputs = NULL;
    }
    if (!reader || !inputs || !outputs ||
        ntasks < KDTREE_DIRECT_DSS_WAVE_MIN_TASKS ||
        ntasks > KDTREE_DIRECT_DSS_WAVE_TASKS ||
        !index_shard_worker_context_active()) {
        goto hard_failure;
    }
    memset(tasks, 0, sizeof(tasks));
    for (task_index = 0U; task_index < ntasks; task_index++) {
        const kdtree_direct_dss_task_input_t* input =
            &inputs[task_index];
        unsigned long long points = 0ULL;
        size_t span_index;

        if (!input->data || !input->nspans ||
            input->nspans > KDTREE_DIRECT_DSS_TASK_SPANS ||
            !input->cover_points ||
            input->cover_points > KDTREE_DIRECT_DSS_TASK_POINTS) {
            goto hard_failure;
        }
        for (span_index = 0U;
             span_index < input->nspans;
             span_index++) {
            u16 first = input->span_first[span_index];
            u16 last = input->span_last[span_index];

            if (first > last || (u32)last >= input->cover_points) {
                goto hard_failure;
            }
            points += (unsigned long long)(last - first) + 1ULL;
        }
        if (!points || points > KDTREE_DIRECT_DSS_TASK_POINTS) {
            goto hard_failure;
        }
        tasks[task_index].input = input;
        tasks[task_index].input_bytes = sizeof(*input);
        tasks[task_index].output = &reader->helper_outputs[task_index];
        tasks[task_index].output_bytes =
            sizeof(reader->helper_outputs[task_index]);
        tasks[task_index].work_units = points * 4ULL;
    }

    run_status = index_shard_helper_run(
        &solver_codekd_helper_ops, tasks, ntasks, NULL);
    if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE ||
        run_status == INDEX_SHARD_HELPER_STOPPED) {
        return 1;
    }
    if (run_status != INDEX_SHARD_HELPER_OK) {
        goto hard_failure;
    }
    for (task_index = 0U; task_index < ntasks; task_index++) {
        if (tasks[task_index].execute_status !=
                INDEX_SHARD_HELPER_TASK_OK ||
            reader->helper_outputs[task_index].status) {
            goto hard_failure;
        }
    }
    *outputs = reader->helper_outputs;
    return 0;

hard_failure:
    if (reader) {
        reader->helper_hard_failure = TRUE;
    }
    errno = EPROTO;
    return -1;
}

static anbool solver_codekd_direct_capable(
    const kdtree_t* tree) {
    fitsbin_t* fb;
    const void* map_base;
    const void* range_start;
    size_t map_size;
    size_t range_size;
    size_t data_elements;
    size_t data_bytes;
    size_t perm_bytes;

    if (!tree || !tree->io || !tree->io_is_fitsbin ||
        tree->treetype != KDTT_DSS ||
        tree->ndata <= 0 || tree->ndim <= 0 ||
        (size_t)tree->ndata >
            SIZE_MAX / (size_t)tree->ndim) {
        return FALSE;
    }
    data_elements =
        (size_t)tree->ndata * (size_t)tree->ndim;
    if (data_elements > SIZE_MAX / sizeof(u16)) {
        return FALSE;
    }
    data_bytes = data_elements * sizeof(u16);
    fb = tree->io;
    if (!fb->payload_fd_initialized ||
        !fb->open_file_stat_valid) {
        return FALSE;
    }
    if (fitsbin_resolve_mapped_range(
            fb,
            tree->data.s,
            data_bytes,
            &map_base,
            &map_size,
            &range_start,
            &range_size) != 1 ||
        range_start != (const void*)tree->data.s ||
        range_size != data_bytes) {
        return FALSE;
    }
    if (!tree->perm) {
        return TRUE;
    }
    if ((size_t)tree->ndata >
        SIZE_MAX / sizeof(u32)) {
        return FALSE;
    }
    perm_bytes =
        (size_t)tree->ndata * sizeof(u32);
    return fitsbin_resolve_mapped_range(
               fb,
               tree->perm,
               perm_bytes,
               &map_base,
               &map_size,
               &range_start,
               &range_size) == 1 &&
        range_start == (const void*)tree->perm &&
        range_size == perm_bytes;
}

static kdtree_qres_t* solver_codekd_rangesearch(
    const kdtree_t* tree,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options) {
    solver_codekd_reader_t* reader;
    kdtree_qres_t* direct;
    kdtree_direct_dss_executor_t executor;
    anbool direct_owned = result == NULL;
    int direct_errno;

    if (tree && tree->io && tree->io_is_fitsbin) {
        /*
         * Keep the original zero-copy shared mapping as the normal path.
         * Advice failure does not invalidate that authoritative mapping.
         */
        return kdtree_rangesearch_options_reuse(
            tree,
            result,
            query,
            maxd2,
            options);
    }
    if (!solver_codekd_direct_capable(tree)) {
        return kdtree_rangesearch_options_reuse(
            tree,
            result,
            query,
            maxd2,
            options);
    }
    reader = solver_codekd_reader_get();
    if (!reader) {
        return kdtree_rangesearch_options_reuse(
            tree,
            result,
            query,
            maxd2,
            options);
    }
    reader->helper_hard_failure = FALSE;
    if (solver_codekd_reader_bind(
            reader,
            tree,
            (fitsbin_t*)tree->io)) {
        return kdtree_rangesearch_options_reuse(
            tree,
            result,
            query,
            maxd2,
            options);
    }
    executor.opaque = reader;
    executor.available = solver_codekd_helper_available;
    executor.run = solver_codekd_helper_run;
    errno = 0;
    direct = kdtree_rangesearch_direct_dss_leased(
        tree,
        result,
        query,
        maxd2,
        options,
        SOLVER_CODEKD_DIRECT_MAX_POINTS,
        0U,
        solver_codekd_read_leased_ranges,
        solver_codekd_release_range,
        reader,
        &executor);
    direct_errno = errno;
    solver_codekd_reader_unbind(reader);
    if (reader->helper_hard_failure || direct_errno == EPROTO) {
        if (direct_owned) {
            if (direct) {
                kdtree_free_query(direct);
            }
        } else {
            kdtree_free_query(result);
        }
        logerr("[solver-codekd-helper] hard leaf-wave failure "
               "index=%s\n",
               tree->name ? tree->name : "(unnamed)");
        errno = EPROTO;
        return NULL;
    }
    if (direct) {
        return direct;
    }
    errno = direct_errno;
    if (!reader->active_segment->fallback_reported) {
        logverb(
            "[solver-codekd-io] mode=mmap-fallback "
            "index=%s errno=%i reason=%s\n",
            tree->name ? tree->name : "(unnamed)",
            errno,
            strerror(errno));
        reader->active_segment->fallback_reported = TRUE;
    }
    return kdtree_rangesearch_options_reuse(
        tree,
        result,
        query,
        maxd2,
        options);
}

typedef struct solver_pair_geometry {
    anbool scale_ok;
    double scale;
    double costheta;
    double sintheta;
    double rel_field_noise2;
} solver_pair_geometry_t;

struct solver_field_geometry {
    const starxy_t* fieldxy;
    solver_pair_geometry_t* pair_geometry;
    int numxy;
    double codetol;
    double verify_pix;
    double quadsize_min;
    double quadsize_max;
    size_t pair_capacity;
    size_t scale_ok_pairs;
    size_t bytes;
    unsigned long long reused_solver_runs;
};

static void solver_free_field_geometry_object(
    solver_field_geometry_t* geometry) {
    if (!geometry) {
        return;
    }
    free(geometry->pair_geometry);
    free(geometry);
}

static void solver_release_field_geometry(solver_t* solver) {
    solver_field_geometry_t* geometry;

    if (!solver || !solver->field_geometry) {
        return;
    }
    if (!solver->field_geometry_owned) {
        solver->field_geometry = NULL;
        return;
    }

    geometry = solver->field_geometry;
    logverb("[solver-geometry] stats mode=compact-triangular "
            "objects=%i pairs=%zu bytes=%zu reused_solver_runs=%llu\n",
            geometry->numxy,
            geometry->scale_ok_pairs,
            geometry->bytes,
            __atomic_load_n(
                &geometry->reused_solver_runs,
                __ATOMIC_RELAXED));
    solver_free_field_geometry_object(geometry);
    solver->field_geometry = NULL;
    solver->field_geometry_owned = FALSE;
}


#if TESTING_TRYALLCODES
#define DEBUGSOLVER 1
#define TRY_ALL_CODES(pq, fieldstars, dimquad, solver, tol2, presult) \
    test_try_all_codes((pquad*)(pq), \
                       (int*)(fieldstars), \
                       (dimquad), \
                       (solver), \
                       (tol2))
void test_try_all_codes(pquad* pq,
                        int* fieldstars, int dimquad,
                        solver_t* solver, double tol2);

#else
#define TRY_ALL_CODES(pq, fieldstars, dimquad, solver, tol2, presult) \
    try_all_codes((pq), \
                  (fieldstars), \
                  (dimquad), \
                  (solver), \
                  (tol2), \
                  (presult))
#endif

#if TESTING_TRYPERMUTATIONS
#define DEBUGSOLVER 1
#define TEST_TRY_PERMUTATIONS test_try_permutations
void test_try_permutations(int* stars, double* code, int dimquad, solver_t* s);

#else
#define TEST_TRY_PERMUTATIONS(u,v,x,y)  // no-op.
#endif

static inline anbool solver_poll_worker_stop(solver_t* solver) {
    if (!index_shard_worker_stop_requested()) {
        return FALSE;
    }

    solver->quit_now = TRUE;
    return TRUE;
}

static void solver_discard_verified_match(MatchObj* mo) {
    verify_free_matchobj(mo);
    if (mo->sip) {
        sip_free(mo->sip);
        mo->sip = NULL;
    }
}




void solver_set_keep_logodds(solver_t* solver, double logodds) {
    solver->logratio_tokeep = logodds;
}

int solver_set_parity(solver_t* solver, int parity) {
    if (!((parity == PARITY_NORMAL) || (parity == PARITY_FLIP) || (parity == PARITY_BOTH))) {
        ERROR("Invalid parity value: %i", parity);
        return -1;
    }
    solver->parity = parity;
    return 0;
}

anbool solver_did_solve(const solver_t* solver) {
    return solver->best_match_solves;
}

void solver_get_quad_size_range_arcsec(const solver_t* solver, double* qmin, double* qmax) {
    if (qmin) {
        *qmin = solver->quadsize_min * solver_get_pixscale_low(solver);
    }
    if (qmax) {
        double q = solver->quadsize_max;
        if (q == 0)
            q = solver->field_diag;
        *qmax = q * solver_get_pixscale_high(solver);
    }
}

double solver_get_field_jitter(const solver_t* solver) {
    return solver->verify_pix;
}

void solver_get_field_center(const solver_t* solver, double* px, double* py) {
    if (px)
        *px = (solver->field_maxx + solver->field_minx)/2.0;
    if (py)
        *py = (solver->field_maxy + solver->field_miny)/2.0;
}

double solver_get_max_radius_arcsec(const solver_t* solver) {
    return solver->funits_upper * solver->field_diag / 2.0;
}

MatchObj* solver_get_best_match(solver_t* solver) {
    return &(solver->best_match);
}

const char* solver_get_best_match_index_name(const solver_t* solver) {
    return solver->best_index->indexname;
}

double solver_get_pixscale_low(const solver_t* solver) {
    return solver->funits_lower;
}
double solver_get_pixscale_high(const solver_t* solver) {
    return solver->funits_upper;
}

void solver_set_quad_size_range(solver_t* solver, double qmin, double qmax) {
    solver->quadsize_min = qmin;
    solver->quadsize_max = qmax;
}

void solver_set_quad_size_fraction(solver_t* solver, double qmin, double qmax) {
    solver_set_quad_size_range(solver, qmin * MIN(solver_field_width(solver), solver_field_height(solver)),
                               qmax * solver->field_diag);
}

void solver_tweak2(solver_t* sp, MatchObj* mo, int order, sip_t* verifysip) {
    double* xy = NULL;
    int Nxy;
    double indexjitter;
    // quad center
    double qc[2];
    // quad radius-squared
    double Q2;
    // initial WCS
    sip_t startsip;
    int* theta;
    double* odds;
    double* refradec;
    int i;
    double newodds;
    int nm, nc, nd;
    int besti;
    int startorder;

    indexjitter = mo->index_jitter; // ref cat positional error, in arcsec.
    xy = starxy_to_xy_array(sp->fieldxy, NULL);
    Nxy = starxy_n(sp->fieldxy);
    qc[0] = (mo->quadpix[0] + mo->quadpix[2]) / 2.0;
    qc[1] = (mo->quadpix[1] + mo->quadpix[3]) / 2.0;
    Q2 = 0.25 * distsq(mo->quadpix, mo->quadpix + 2, 2);
    if (Q2 == 0.0) {
        // can happen if we're verifying an existing WCS
        // note, this is radius-squared, so 1e6 is not crazy.
        Q2 = 1e6;
        // set qc to the image center here?  or crpix?
        logverb("solver_tweak2(): setting Q2=%g; qc=(%g,%g)\n", Q2, qc[0], qc[1]);
    }

    // mo->refradec may be NULL at this point, so get it from refxyz instead...
    refradec = malloc(3 * mo->nindex * sizeof(double));
    for (i=0; i<mo->nindex; i++)
        xyzarr2radecdegarr(mo->refxyz + i*3, refradec + i*2);

    // Verifying an existing WCS?
    if (verifysip) {
        memcpy(&startsip, verifysip, sizeof(sip_t));
        startorder = MIN(verifysip->a_order, sp->tweak_aborder);
    } else {
        startorder = 1;
        sip_wrap_tan(&(mo->wcstan), &startsip);
    }

    startsip.ap_order = startsip.bp_order = sp->tweak_abporder;
    startsip.a_order = startsip.b_order = sp->tweak_aborder;
    logverb("solver_tweak2: setting orders %i, %i\n", sp->tweak_aborder, sp->tweak_abporder);

    // for TWEAK_DEBUG_PLOTs
    theta = mo->theta;
    besti = mo->nbest-1;//mo->nmatch + mo->nconflict + mo->ndistractor;

    logverb("solver_tweak2: set_crpix %i, crpix (%.1f,%.1f)\n",
            sp->set_crpix, sp->crpix[0], sp->crpix[1]);
    mo->sip = tweak2(xy, Nxy,
                     sp->verify_pix, // pixel positional noise sigma
                     solver_field_width(sp),
                     solver_field_height(sp),
                     refradec, mo->nindex,
                     indexjitter, qc, Q2,
                     sp->distractor_ratio,
                     sp->logratio_bail_threshold,
                     order, sp->tweak_abporder,
                     &startsip, NULL, &theta, &odds,
                     sp->set_crpix ? sp->crpix : NULL,
                     &newodds, &besti, mo->testperm, startorder);
    free(refradec);

    // FIXME -- update refxy?  Nobody uses it, right?
    free(mo->refxy);
    mo->refxy = NULL;
    // FIXME -- and testperm?
    free(mo->testperm);
    mo->testperm = NULL;

    if (mo->sip) {
        // Yoink the TAN solution (?)
        memcpy(&(mo->wcstan), &(mo->sip->wcstan), sizeof(tan_t));

        // Plug in the new "theta" and "odds".
        free(mo->theta);
        free(mo->matchodds);
        mo->theta = theta;
        mo->matchodds = odds;

        mo->logodds = newodds;

        verify_count_hits(theta, besti, &nm, &nc, &nd);
        mo->nmatch = nm;
        mo->nconflict = nc;
        mo->ndistractor = nd;
        matchobj_compute_derived(mo);
    }
    free(xy);
}

void solver_log_params(const solver_t* sp) {
    int i;
    logverb("Solver:\n");
    logverb("  Arcsec per pix range: %g, %g\n", sp->funits_lower, sp->funits_upper);
    logverb("  Image size: %g x %g\n", solver_field_width(sp), solver_field_height(sp));
    logverb("  Quad size range: %g, %g\n", sp->quadsize_min, sp->quadsize_max);
    logverb("  Objs: %i, %i\n", sp->startobj, sp->endobj);
    logverb("  Parity: %i, %s\n", sp->parity, sp->parity == PARITY_NORMAL ? "normal" : (sp->parity == PARITY_FLIP ? "flip" : "both"));
    if (sp->use_radec) {
        double ra,dec,rad;
        xyzarr2radecdeg(sp->centerxyz, &ra, &dec);
        rad = distsq2deg(sp->r2);
        logverb("  Use_radec? yes, (%g, %g), radius %g deg\n", ra, dec, rad);
    } else {
        logverb("  Use_radec? no\n");
    }
    logverb("  Pixel xscale: %g\n", sp->pixel_xscale);
    logverb("  Verify_pix: %g\n", sp->verify_pix);
    logverb("  Code tol: %g\n", sp->codetol);
    logverb("  Dist from quad bonus: %s\n", sp->distance_from_quad_bonus ? "yes" : "no");
    logverb("  Distractor ratio: %g\n", sp->distractor_ratio);
    logverb("  Log tune-up threshold: %g\n", sp->logratio_totune);
    logverb("  Log bail threshold: %g\n", sp->logratio_bail_threshold);
    logverb("  Log stoplooking threshold: %g\n", sp->logratio_stoplooking);
    logverb("  Maxquads %i\n", sp->maxquads);
    logverb("  Maxmatches %i\n", sp->maxmatches);
    logverb("  Set CRPIX? %s", sp->set_crpix ? "yes" : "no\n");
    if (sp->set_crpix) {
        if (sp->set_crpix_center)
            logverb(", center\n");
        else
            logverb(", %g, %g\n", sp->crpix[0], sp->crpix[1]);
    }
    logverb("  Tweak? %s\n", sp->do_tweak ? "yes" : "no");
    if (sp->do_tweak) {
        logverb("    Forward order %i\n", sp->tweak_aborder);
        logverb("    Reverse order %i\n", sp->tweak_abporder);
    }
    logverb("  Indexes: %zu\n", pl_size(sp->indexes));
    for (i=0; i<pl_size(sp->indexes); i++) {
        index_t* ind = pl_get(sp->indexes, i);
        logverb("    %s\n", ind->indexname);
    }
    if (sp->fieldxy) {
      logverb("  Field (processed): %i stars\n", starxy_n(sp->fieldxy));
      for (i=0; i<starxy_n(sp->fieldxy); i++) {
        debug("    xy (%.1f, %.1f), flux %.1f\n",
              starxy_getx(sp->fieldxy, i), starxy_gety(sp->fieldxy, i),
              sp->fieldxy->flux ? starxy_get_flux(sp->fieldxy, i) : 0.0);
      }
    }
    if (sp->fieldxy_orig) {
      logverb("  Field (orig): %i stars\n", starxy_n(sp->fieldxy_orig));
      for (i=0; i<starxy_n(sp->fieldxy_orig); i++) {
        debug("    xy (%.1f, %.1f), flux %.1f\n",
              starxy_getx(sp->fieldxy_orig, i), starxy_gety(sp->fieldxy_orig, i),
              sp->fieldxy_orig->flux ? starxy_get_flux(sp->fieldxy_orig, i) : 0.0);
      }
    }
}


void solver_print_to(const solver_t* sp, FILE* stream) {
    //int oldlevel = log_get_level();
    FILE* oldfid = log_get_fid();
    //log_set_level(LOG_ALL);
    log_to(stream);
    solver_log_params(sp);
    //log_set_level(oldlevel);
    log_to(oldfid);
}

/*
 static MatchObj* matchobj_copy_deep(const MatchObj* mo, MatchObj* dest) {
 if (!dest)
 dest = malloc(sizeof(MatchObj));
 memcpy(dest, mo, sizeof(MatchObj));
 // various modules add things to a mo...
 onefield_matchobj_deep_copy(mo, dest);
 verify_matchobj_deep_copy(mo, dest);
 return dest;
 }

 static void matchobj_free_data(MatchObj* mo) {
 verify_free_matchobj(mo);
 onefield_free_matchobj(mo);
 }
 */

static const int A = 0, B = 1, C = 2, D = 3;

// Number of stars in the "backbone" of the quad: stars A and B.
static const int NBACK = 2;

static void find_field_boundaries(solver_t* solver);

static inline double getx(const double* d, int ind) {
    return d[ind*2];
}
static inline double gety(const double* d, int ind) {
    return d[ind*2 + 1];
}
static inline void setx(double* d, int ind, double val) {
    d[ind*2] = val;
}
static inline void sety(double* d, int ind, double val) {
    d[ind*2 + 1] = val;
}

static void field_getxy(
    const solver_t* sp,
    int index,
    double* x,
    double* y) {
    *x = starxy_getx(sp->fieldxy, index);
    *y = starxy_gety(sp->fieldxy, index);
}

static double field_getx(const solver_t* sp, int index) {
    return starxy_getx(sp->fieldxy, index);
}
static double field_gety(const solver_t* sp, int index) {
    return starxy_gety(sp->fieldxy, index);
}

static void update_timeused(solver_t* sp) {
    double usertime, systime;
    get_resource_stats(&usertime, &systime, NULL);
    sp->timeused = (usertime + systime) - sp->starttime;
    if (sp->timeused < 0.0)
        sp->timeused = 0.0;
}

static void set_matchobj_template(solver_t* solver, MatchObj* mo) {
    if (solver->mo_template)
        memcpy(mo, solver->mo_template, sizeof(MatchObj));
    else
        memset(mo, 0, sizeof(MatchObj));
}

static void get_field_center(solver_t* s, double* cx, double* cy) {
    *cx = 0.5 * (s->field_minx + s->field_maxx);
    *cy = 0.5 * (s->field_miny + s->field_maxy);
}

static void get_field_ll_corner(solver_t* s, double* lx, double* ly) {
    *lx = s->field_minx;
    *ly = s->field_miny;
}

void solver_reset_counters(solver_t* s) {
    s->quit_now = FALSE;
    s->have_best_match = FALSE;
    s->best_match_solves = FALSE;
    s->numtries = 0;
    s->nummatches = 0;
    s->numscaleok = 0;
    s->last_examined_object = 0;
    s->num_cxdx_skipped = 0;
    s->num_radec_skipped = 0;
    s->num_abscale_skipped = 0;
    s->num_verified = 0;
}

double solver_field_width(const solver_t* s) {
    return s->field_maxx - s->field_minx;
}
double solver_field_height(const solver_t* s) {
    return s->field_maxy - s->field_miny;
}

void solver_set_radec(solver_t* s, double ra, double dec, double radius_deg) {
    s->use_radec = TRUE;
    radecdeg2xyzarr(ra, dec, s->centerxyz);
    s->r2 = deg2distsq(radius_deg);
}

void solver_clear_radec(solver_t* s) {
    s->use_radec = FALSE;
}

static void set_center_and_radius(solver_t* solver, MatchObj* mo,
                                  tan_t* tan, sip_t* sip) {
    double cx, cy, lx, ly;
    double xyz[3];
    get_field_center(solver, &cx, &cy);
    get_field_ll_corner(solver, &lx, &ly);
    if (sip) {
        sip_pixelxy2xyzarr(sip, cx, cy, mo->center);
        sip_pixelxy2xyzarr(sip, lx, ly, xyz);
    } else {
        tan_pixelxy2xyzarr(tan, cx, cy, mo->center);
        tan_pixelxy2xyzarr(tan, lx, ly, xyz);
    }
    mo->radius = sqrt(distsq(mo->center, xyz, 3));
    mo->radius_deg = dist2deg(mo->radius);
}
static void set_index(solver_t* s, index_t* index) {
    s->index = index;
    s->rel_index_noise2 = square(index->index_jitter / index->index_scale_lower);
}

static void set_diag(solver_t* s) {
    s->field_diag = hypot(solver_field_width(s), solver_field_height(s));
}

void solver_set_field(solver_t* s, starxy_t* field) {
    solver_free_field(s);

    /*
     * Reset the compatibility policy counters at the field boundary.
     * Shards keep topology and CodeKD at NORMAL while sparse Star and Quad
     * payloads use RANDOM. Serial solving retains original NORMAL mappings.
     */
    fitsbin_mmap_advice_state_reset(
        &s->index_mmap_policy);

    s->fieldxy_orig = field;

    // Preprocessing happens in "solver_preprocess_field()".
}

void solver_set_field_bounds(solver_t* s, double xlo, double xhi, double ylo, double yhi) {
    s->field_minx = xlo;
    s->field_maxx = xhi;
    s->field_miny = ylo;
    s->field_maxy = yhi;
    set_diag(s);
}

void solver_cleanup_field(solver_t* solver) {
    solver_reset_best_match(solver);
    solver_free_field(solver);
    solver->fieldxy = NULL;
    solver_reset_counters(solver);
}

void solver_verify_sip_wcs(solver_t* solver, sip_t* sip) { //, MatchObj* pmo) {
    int i, nindexes;
    MatchObj mo;
    MatchObj* pmo;
    anbool olddqb;

    pmo = &mo;

    if (!solver->vf)
        solver_preprocess_field(solver);

    // fabricate a match and inject it into the solver.
    set_matchobj_template(solver, pmo);
    memcpy(&(mo.wcstan), &(sip->wcstan), sizeof(tan_t));
    mo.wcs_valid = TRUE;
    mo.scale = sip_pixel_scale(sip);
    set_center_and_radius(solver, pmo, NULL, sip);
    olddqb = solver->distance_from_quad_bonus;
    solver->distance_from_quad_bonus = FALSE;

    nindexes = pl_size(solver->indexes);
    for (i=0; i<nindexes; i++) {
        index_t* index = pl_get(solver->indexes, i);
        set_index(solver, index);
        solver_inject_match(solver, pmo, sip);
    }

    // revert
    solver->distance_from_quad_bonus = olddqb;
}

void solver_add_index(solver_t* solver, index_t* index) {
    pl_append(solver->indexes, index);
}

int solver_n_indices(const solver_t* solver) {
    return pl_size(solver->indexes);
}

index_t* solver_get_index(const solver_t* solver, int i) {
    return pl_get(solver->indexes, i);
}

void solver_reset_best_match(solver_t* sp) {
    // we don't really care about very bad best matches...
    sp->best_logodds = 0;
    memset(&(sp->best_match), 0, sizeof(MatchObj));
    sp->best_index = NULL;
    sp->best_match_solves = FALSE;
    sp->have_best_match = FALSE;
}

void solver_compute_quad_range(const solver_t* sp, const index_t* index,
                               double* minAB, double* maxAB) {
    double scalefudge; // in pixels

    // compute fudge factor for quad scale: what are the extreme
    // ranges of quad scales that should be accepted, given the
    // code tolerance?
    // -what is the maximum number of pixels a C or D star can move
    //  to singlehandedly exceed the code tolerance?
    // -largest quad
    // -smallest arcsec-per-pixel scale

    // -index_scale_upper * 1/sqrt(2) is the side length of
    //  the unit-square of code space, in arcseconds.
    // -that times the code tolerance is how far a C/D star
    //  can move before exceeding the code tolerance, in arcsec.
    // -that divided by the smallest arcsec-per-pixel scale
    //  gives the largest motion in pixels.

    //logverb("Index scale %f, %f\n",
    //index->index_scale_upper, index->index_scale_lower);

    scalefudge = index->index_scale_upper * M_SQRT1_2 *
        sp->codetol / sp->funits_upper;

    if (sp->funits_upper != 0.0) {
        *minAB = index->index_scale_lower / sp->funits_upper;
        *minAB -= scalefudge;
    }
    if (sp->funits_lower != 0.0) {
        *maxAB = index->index_scale_upper / sp->funits_lower;
        *maxAB += scalefudge;
    }
}
static void try_all_codes(
    const pquad* pq,
    const int* fieldstars,
    int dimquad,
    solver_t* solver,
    double tol2,
    kdtree_qres_t** presult);

static void try_all_codes_2(
    const int* fieldstars,
    int dimquad,
    const double* code,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    kdtree_qres_t** presult);

static void try_permutations(
    const int* origstars,
    int dimquad,
    const double* origcode,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    int* stars,
    double* code,
    int slot,
    anbool* placed,
    kdtree_qres_t** presult);

static void solver_execute_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    kdtree_qres_t** presult);

static void resolve_matches(kdtree_qres_t* krez,
                            const double* field,
                            const int* fstars,
                            int dimquads,
                            int quads_tried,
                            solver_t* solver,
                            anbool current_parity);

static int solver_handle_hit(solver_t* sp,
                             MatchObj* mo,
                             sip_t* sip,
                             anbool fake_match);

static double solver_prepare_hit_for_verify(
    solver_t* sp,
    MatchObj* mo,
    double* logaccept);

static int solver_handle_hit_after_verify(
    solver_t* sp,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match,
    double match_distance_in_pixels2);

static uint64_t solver_order_hash_mix(
    uint64_t state,
    uint64_t value);

void solver_profile_accumulate(solver_profile_t* total,
                               const solver_profile_t* profile) {
    if (!total || !profile) {
        return;
    }

    total->detailed = total->detailed || profile->detailed;
    total->execution_failed =
        total->execution_failed || profile->execution_failed;

    total->solver_run_wall_seconds +=
        profile->solver_run_wall_seconds;
    total->codekd_wall_seconds += profile->codekd_wall_seconds;
    total->resolve_wall_seconds += profile->resolve_wall_seconds;
    total->verify_wall_seconds += profile->verify_wall_seconds;
    total->hypothesis_wave_wall_seconds +=
        profile->hypothesis_wave_wall_seconds;
    total->ab_planning_wall_seconds +=
        profile->ab_planning_wall_seconds;

    total->codekd_calls += profile->codekd_calls;
    total->codekd_hits += profile->codekd_hits;
    total->resolve_calls += profile->resolve_calls;
    total->verify_calls += profile->verify_calls;
    total->hypothesis_batches += profile->hypothesis_batches;
    total->hypothesis_batches_completed +=
        profile->hypothesis_batches_completed;
    total->hypothesis_batches_stopped +=
        profile->hypothesis_batches_stopped;
    total->hypothesis_batches_failed +=
        profile->hypothesis_batches_failed;
    total->hypotheses_generated += profile->hypotheses_generated;
    total->hypotheses_executed += profile->hypotheses_executed;
    total->hypotheses_reduced += profile->hypotheses_reduced;
    total->task_ranges_planned += profile->task_ranges_planned;
    total->task_ranges_executed += profile->task_ranges_executed;
    total->task_ranges_submitted += profile->task_ranges_submitted;
    total->task_ranges_inline += profile->task_ranges_inline;
    total->parallel_batches += profile->parallel_batches;
    total->parallel_batches_observed +=
        profile->parallel_batches_observed;
    total->parallel_hypotheses += profile->parallel_hypotheses;
    total->allocation_failures += profile->allocation_failures;
    total->search_failures += profile->search_failures;
    total->ab_blocks_planned += profile->ab_blocks_planned;
    total->ab_blocks_retired += profile->ab_blocks_retired;
    total->ab_blocks_owner += profile->ab_blocks_owner;
    total->ab_segments_retired += profile->ab_segments_retired;
    total->ab_segment_payload_bytes +=
        profile->ab_segment_payload_bytes;
    total->ab_pairs_planned += profile->ab_pairs_planned;
    total->ab_combinations_planned +=
        profile->ab_combinations_planned;
    total->ab_intra_pair_splits +=
        profile->ab_intra_pair_splits;
    total->ab_max_pair_combinations =
        MAX(total->ab_max_pair_combinations,
            profile->ab_max_pair_combinations);
    total->ab_helper_tasks += profile->ab_helper_tasks;
    total->ab_helper_combinations +=
        profile->ab_helper_combinations;
    if (!total->hypothesis_order_hash) {
        total->hypothesis_order_hash =
            profile->hypothesis_order_hash;
    } else if (profile->hypothesis_order_hash) {
        total->hypothesis_order_hash =
            solver_order_hash_mix(
                total->hypothesis_order_hash,
                profile->hypothesis_order_hash);
    }
    if (!total->kd_result_order_hash) {
        total->kd_result_order_hash =
            profile->kd_result_order_hash;
    } else if (profile->kd_result_order_hash) {
        total->kd_result_order_hash =
            solver_order_hash_mix(
                total->kd_result_order_hash,
                profile->kd_result_order_hash);
    }
    if (!total->candidate_order_hash) {
        total->candidate_order_hash =
            profile->candidate_order_hash;
    } else if (profile->candidate_order_hash) {
        total->candidate_order_hash =
            solver_order_hash_mix(
                total->candidate_order_hash,
                profile->candidate_order_hash);
    }

    if (profile->max_batch_hypotheses >
        total->max_batch_hypotheses) {
        total->max_batch_hypotheses =
            profile->max_batch_hypotheses;
    }

    if (profile->max_task_ranges > total->max_task_ranges) {
        total->max_task_ranges = profile->max_task_ranges;
    }

    if (profile->max_parallel_ranges >
        total->max_parallel_ranges) {
        total->max_parallel_ranges = profile->max_parallel_ranges;
    }
}

static void solver_profile_report(const solver_t* solver) {
    const solver_profile_t* profile;
    double reduction_wall_seconds;

    if (!solver || !index_shard_trace_enabled()) {
        return;
    }

    profile = &solver->profile;
    reduction_wall_seconds =
        profile->resolve_wall_seconds -
        profile->verify_wall_seconds;

    if (reduction_wall_seconds < 0.0) {
        reduction_wall_seconds = 0.0;
    }

    logmsg("[solver] phase-profile detailed=%i failed=%i "
           "solver_run_elapsed=%.6f codekd_work_wall_sum=%.6f "
           "codekd_calls=%llu codekd_hits=%llu "
           "resolve_work_wall_sum=%.6f "
           "reduction_ex_verify_hit_work_wall_sum=%.6f "
           "resolve_calls=%llu verify_hit_work_wall_sum=%.6f "
           "verify_calls=%llu hypothesis_wave_elapsed_sum=%.6f "
           "hypothesis_batches=%llu parallel_batches=%llu "
           "parallel_batches_observed=%llu "
           "parallel_hypotheses=%llu "
           "task_ranges_planned=%llu task_ranges_executed=%llu "
           "task_ranges_inline=%llu allocation_failures=%llu "
           "search_failures=%llu "
           "ab_blocks_planned=%llu ab_blocks_retired=%llu "
           "ab_blocks_owner=%llu ab_segments_retired=%llu "
           "ab_segment_payload_bytes=%llu ab_pairs=%llu "
           "ab_combinations=%llu intra_pair_splits=%llu "
           "max_pair_combinations=%llu ab_plan_wall=%.6f "
           "helper_tasks=%llu helper_combinations=%llu "
           "hypothesis_order=%016llx "
           "kd_result_order=%016llx candidate_order=%016llx\n",
           profile->detailed ? 1 : 0,
           profile->execution_failed ? 1 : 0,
           profile->solver_run_wall_seconds,
           profile->codekd_wall_seconds,
           profile->codekd_calls,
           profile->codekd_hits,
           profile->resolve_wall_seconds,
           reduction_wall_seconds,
           profile->resolve_calls,
           profile->verify_wall_seconds,
           profile->verify_calls,
           profile->hypothesis_wave_wall_seconds,
           profile->hypothesis_batches,
           profile->parallel_batches,
           profile->parallel_batches_observed,
           profile->parallel_hypotheses,
           profile->task_ranges_planned,
           profile->task_ranges_executed,
           profile->task_ranges_inline,
           profile->allocation_failures,
           profile->search_failures,
           profile->ab_blocks_planned,
           profile->ab_blocks_retired,
           profile->ab_blocks_owner,
           profile->ab_segments_retired,
           profile->ab_segment_payload_bytes,
           profile->ab_pairs_planned,
           profile->ab_combinations_planned,
           profile->ab_intra_pair_splits,
           profile->ab_max_pair_combinations,
           profile->ab_planning_wall_seconds,
           profile->ab_helper_tasks,
           profile->ab_helper_combinations,
           profile->hypothesis_order_hash,
           profile->kd_result_order_hash,
           profile->candidate_order_hash);
}

static void check_scale(pquad* pq, solver_t* s) {
    double dx, dy;
    dx = field_getx(s, pq->fieldB) - field_getx(s, pq->fieldA);
    dy = field_gety(s, pq->fieldB) - field_gety(s, pq->fieldA);
    pq->scale = dx*dx + dy*dy;
    if ((pq->scale < s->minminAB2) ||
        (pq->scale > s->maxmaxAB2)) {
        pq->scale_ok = FALSE;
        return;
    }
    pq->costheta = (dy + dx) / pq->scale;
    pq->sintheta = (dy - dx) / pq->scale;
    pq->rel_field_noise2 = (s->verify_pix * s->verify_pix) / pq->scale;
    pq->scale_ok = TRUE;
}

/*
 * Keep the native pquad path and the compact triangular-geometry path on one
 * explicitly contracted arithmetic sequence.  With -march=native -O3 GCC
 * otherwise vectorizes check_inbox() but scalarizes the compact caller, and
 * is free to choose the opposite multiplication as the fused addend.  Both
 * expressions are mathematically equivalent, but their last bits can differ
 * and then perturb the scientific traversal hash between W1 and W>1.
 */
static inline void solver_transform_code_coordinates(
    double relative_x,
    double relative_y,
    double costheta,
    double sintheta,
    double* code_x,
    double* code_y) {
    double x_addend = relative_y * sintheta;
    double y_addend = relative_y * costheta;

    *code_x = fma(relative_x, costheta, x_addend);
    *code_y = fma(-relative_x, sintheta, y_addend);
}

static void check_inbox(pquad* pq, int start, solver_t* solver) {
    int i;
    double Ax, Ay;

    if (start == 0) {
        pq->eligible_count = 0;
    }
    field_getxy(solver, pq->fieldA, &Ax, &Ay);
    // check which C, D points are inside the circle.
    for (i = start; i < pq->ninbox; i++) {
        double r;
        double Cx, Cy;
        double tol = solver->codetol;
        if (!pq->inbox[i]) {
            if (pq->inbox_prefix) {
                pq->inbox_prefix[i] =
                    (uint16_t)pq->eligible_count;
            }
            continue;
        }
        field_getxy(solver, i, &Cx, &Cy);
        Cx -= Ax;
        Cy -= Ay;
        solver_transform_code_coordinates(
            Cx,
            Cy,
            pq->costheta,
            pq->sintheta,
            &Cx,
            &Cy);

        // make sure it's in the circle centered at (0.5, 0.5)
        // with radius 1/sqrt(2) (plus codetol for fudge):
        // (x-1/2)^2 + (y-1/2)^2   <=   (r + codetol)^2
        // x^2-x+1/4 + y^2-y+1/4   <=   (1/sqrt(2) + codetol)^2
        // x^2-x + y^2-y + 1/2     <=   1/2 + sqrt(2)*codetol + codetol^2
        // x^2-x + y^2-y           <=   sqrt(2)*codetol + codetol^2
        r = (Cx * Cx - Cx) + (Cy * Cy - Cy);
        if (r > (tol * (M_SQRT2 + tol))) {
            pq->inbox[i] = FALSE;
            if (pq->inbox_prefix) {
                pq->inbox_prefix[i] =
                    (uint16_t)pq->eligible_count;
            }
            continue;
        }
        setx(pq->xy, i, Cx);
        sety(pq->xy, i, Cy);
        pq->eligible_count++;
        if (pq->inbox_prefix) {
            pq->inbox_prefix[i] =
                (uint16_t)pq->eligible_count;
        }
    }
}

static int solver_allocate_pquad_storage(
    solver_t* solver,
    pquad* pq,
    int numxy) {
    size_t inbox_bytes;
    size_t xy_bytes;

    if (!solver || !pq || numxy <= 0 ||
        (size_t)numxy > SIZE_MAX / sizeof(anbool) ||
        (size_t)numxy > SIZE_MAX /
            (2U * sizeof(double))) {
        if (solver) {
            solver->profile.allocation_failures++;
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
        }
        return -1;
    }
    inbox_bytes = (size_t)numxy * sizeof(anbool);
    xy_bytes = (size_t)numxy * 2U * sizeof(double);
    pq->inbox = malloc(inbox_bytes);
    pq->xy = malloc(xy_bytes);
    if (!pq->inbox || !pq->xy) {
        free(pq->inbox);
        free(pq->xy);
        pq->inbox = NULL;
        pq->xy = NULL;
        solver->profile.allocation_failures++;
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        SYSERROR("Failed to allocate solver pquad storage");
        return -1;
    }
    return 0;
}

static int solver_field_geometry_numxy(const solver_t* solver) {
    int numxy;

    if (!solver || !solver->fieldxy) {
        return 0;
    }
    numxy = starxy_n(solver->fieldxy);
    if (solver->endobj && numxy > solver->endobj) {
        numxy = solver->endobj;
    }
    if (numxy >= 1000) {
        numxy = 1000;
    }
    return numxy;
}

static anbool solver_field_geometry_estimate(
    int numxy,
    size_t* pairs,
    size_t* bytes) {
    size_t n;

    if (numxy <= 0 || !pairs || !bytes) {
        return FALSE;
    }
    n = (size_t)numxy;
    if (n > 1U && n > SIZE_MAX / (n - 1U)) {
        return FALSE;
    }
    *pairs = n * (n - 1U) / 2U;
    if (*pairs >
        (SIZE_MAX - sizeof(solver_field_geometry_t)) /
            sizeof(solver_pair_geometry_t)) {
        return FALSE;
    }
    *bytes = sizeof(solver_field_geometry_t) +
        *pairs * sizeof(solver_pair_geometry_t);
    return TRUE;
}

static void solver_field_geometry_check_scale(
    solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    const solver_t* solver) {
    double dx;
    double dy;
    double minimum_scale;
    double maximum_scale;

    dx = field_getx(solver, field_b) -
        field_getx(solver, field_a);
    dy = field_gety(solver, field_b) -
        field_gety(solver, field_a);
    pair->scale = dx * dx + dy * dy;
    minimum_scale = square(solver->quadsize_min);
    maximum_scale = solver->quadsize_max > 0.0
        ? square(solver->quadsize_max)
        : LARGE_VAL;
    if (pair->scale < minimum_scale ||
        pair->scale > maximum_scale) {
        pair->scale_ok = FALSE;
        return;
    }

    pair->costheta = (dy + dx) / pair->scale;
    pair->sintheta = (dy - dx) / pair->scale;
    pair->rel_field_noise2 =
        (solver->verify_pix * solver->verify_pix) / pair->scale;
    pair->scale_ok = TRUE;
}

static size_t solver_field_geometry_pair_index(
    int field_a,
    int field_b) {
    return (size_t)field_b *
        (size_t)(field_b - 1) / 2U +
        (size_t)field_a;
}

static const solver_pair_geometry_t*
solver_field_geometry_pair(
    const solver_field_geometry_t* geometry,
    int field_a,
    int field_b) {
    size_t pair_index;

    if (!geometry || !geometry->pair_geometry ||
        field_a < 0 || field_b <= field_a ||
        field_b >= geometry->numxy) {
        return NULL;
    }
    pair_index = solver_field_geometry_pair_index(
        field_a,
        field_b);
    if (pair_index >= geometry->pair_capacity) {
        return NULL;
    }
    return geometry->pair_geometry + pair_index;
}

static anbool solver_pair_geometry_transform(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int star,
    double* x,
    double* y) {
    double ax;
    double ay;
    double cx;
    double cy;
    double radius_term;

    if (!geometry || !pair || !pair->scale_ok ||
        !geometry->fieldxy || !x || !y ||
        star < 0 || star >= geometry->numxy ||
        star == field_a || star == field_b) {
        return FALSE;
    }
    ax = starxy_getx(geometry->fieldxy, field_a);
    ay = starxy_gety(geometry->fieldxy, field_a);
    cx = starxy_getx(geometry->fieldxy, star) - ax;
    cy = starxy_gety(geometry->fieldxy, star) - ay;
    solver_transform_code_coordinates(
        cx,
        cy,
        pair->costheta,
        pair->sintheta,
        &cx,
        &cy);
    radius_term =
        (cx * cx - cx) + (cy * cy - cy);
    if (radius_term >
        geometry->codetol *
            (M_SQRT2 + geometry->codetol)) {
        return FALSE;
    }
    *x = cx;
    *y = cy;
    return TRUE;
}

static int solver_pair_geometry_eligible_before(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int fieldtop) {
    int eligible = 0;
    int star;

    for (star = 0; star < fieldtop; star++) {
        double x;
        double y;

        if (solver_pair_geometry_transform(
                geometry,
                pair,
                field_a,
                field_b,
                star,
                &x,
                &y)) {
            eligible++;
        }
    }
    return eligible;
}

static anbool solver_field_geometry_compatible(
    const solver_t* solver,
    int numxy) {
    const solver_field_geometry_t* geometry;

    if (!solver || !solver->field_geometry) {
        return FALSE;
    }
    geometry = solver->field_geometry;
    return geometry->fieldxy == solver->fieldxy &&
        geometry->numxy >= numxy &&
        geometry->codetol == solver->codetol &&
        geometry->verify_pix == solver->verify_pix &&
        geometry->quadsize_min == solver->quadsize_min &&
        geometry->quadsize_max == solver->quadsize_max;
}

void solver_release_incompatible_field_geometry(
    solver_t* solver) {
    int numxy;

    if (!solver || !solver->field_geometry ||
        !solver->fieldxy) {
        return;
    }
    numxy = solver_field_geometry_numxy(solver);
    if (!solver_field_geometry_compatible(
            solver, numxy)) {
        solver_release_field_geometry(solver);
    }
}

anbool solver_prepare_field_geometry(solver_t* solver) {
    solver_field_geometry_t* geometry;
    size_t pairs;
    size_t estimated_bytes;
    size_t built_pairs = 0;
    double wall_start;
    int numxy;
    int fieldA;
    int fieldB;

    if (!solver || !solver->fieldxy) {
        return FALSE;
    }
    wall_start = monotonic_seconds();

    numxy = starxy_n(solver->fieldxy);
    if (numxy >= 1000) {
        numxy = 1000;
    }
    if (solver_field_geometry_compatible(
            solver,
            solver_field_geometry_numxy(solver)) &&
        solver->field_geometry->numxy == numxy) {
        return TRUE;
    }
    solver_release_field_geometry(solver);

    if (solver->startobj >=
        solver_field_geometry_numxy(solver)) {
        logverb("[solver-geometry] mode=legacy reason=empty-range "
                "startobj=%i objects=%i\n",
                solver->startobj,
                numxy);
        return FALSE;
    }
    if (!solver_field_geometry_estimate(
            numxy,
            &pairs,
            &estimated_bytes)) {
        logverb("[solver-geometry] mode=legacy reason=estimate "
                "objects=%i\n",
                numxy);
        return FALSE;
    }
    if (estimated_bytes > SOLVER_FIELD_GEOMETRY_BUDGET_BYTES) {
        logverb("[solver-geometry] mode=legacy reason=budget "
                "objects=%i estimate=%zu budget=%llu\n",
                numxy,
                estimated_bytes,
                (unsigned long long)
                    SOLVER_FIELD_GEOMETRY_BUDGET_BYTES);
        return FALSE;
    }

    geometry = calloc(1, sizeof(*geometry));
    if (!geometry) {
        logverb("[solver-geometry] mode=legacy reason=allocation "
                "objects=%i estimate=%zu\n",
                numxy,
                estimated_bytes);
        return FALSE;
    }
    geometry->pair_capacity = pairs;
    if (pairs) {
        geometry->pair_geometry = calloc(
            pairs,
            sizeof(*geometry->pair_geometry));
    }
    if (pairs && !geometry->pair_geometry) {
        solver_free_field_geometry_object(geometry);
        logverb("[solver-geometry] mode=legacy reason=allocation "
                "objects=%i estimate=%zu\n",
                numxy,
                estimated_bytes);
        return FALSE;
    }

    for (fieldB = 1; fieldB < numxy; fieldB++) {
        for (fieldA = 0; fieldA < fieldB; fieldA++) {
            solver_pair_geometry_t* pair =
                geometry->pair_geometry +
                solver_field_geometry_pair_index(
                    fieldA,
                    fieldB);

            solver_field_geometry_check_scale(
                pair,
                fieldA,
                fieldB,
                solver);
            if (pair->scale_ok) {
                built_pairs++;
            }
        }
    }

    geometry->fieldxy = solver->fieldxy;
    geometry->numxy = numxy;
    geometry->codetol = solver->codetol;
    geometry->verify_pix = solver->verify_pix;
    geometry->quadsize_min = solver->quadsize_min;
    geometry->quadsize_max = solver->quadsize_max;
    geometry->scale_ok_pairs = built_pairs;
    geometry->bytes = estimated_bytes;

    solver->field_geometry = geometry;
    solver->field_geometry_owned = TRUE;

    logverb("[solver-geometry] mode=compact-triangular objects=%i "
            "pairs=%zu possible_pairs=%zu bytes=%zu "
            "per_star_payload=none "
            "prepare_wall=%.6f\n",
            numxy,
            built_pairs,
            pairs,
            estimated_bytes,
            monotonic_seconds() - wall_start);
    return TRUE;
}

#if defined DEBUGSOLVER
static void print_inbox(pquad* pq) {
    int i;
    debug("[ ");
    for (i = 0; i < pq->ninbox; i++) {
        if (pq->inbox[i])
            debug("%i ", i);
    }
    debug("] (n %i)\n", pq->ninbox);
}
#else
static void print_inbox(pquad* pq) {}
#endif


void solver_reset_field_size(solver_t* s) {
    s->field_minx = s->field_maxx = s->field_miny = s->field_maxy = 0;
    s->field_diag = 0.0;
}

static void find_field_boundaries(solver_t* solver) {
    // If the bounds haven't been set, use the bounding box.
    if ((solver->field_minx == solver->field_maxx) ||
        (solver->field_miny == solver->field_maxy)) {
        int i;
        solver->field_minx = solver->field_miny =  LARGE_VAL;
        solver->field_maxx = solver->field_maxy = -LARGE_VAL;
        for (i = 0; i < starxy_n(solver->fieldxy); i++) {
            solver->field_minx = MIN(solver->field_minx, field_getx(solver, i));
            solver->field_maxx = MAX(solver->field_maxx, field_getx(solver, i));
            solver->field_miny = MIN(solver->field_miny, field_gety(solver, i));
            solver->field_maxy = MAX(solver->field_maxy, field_gety(solver, i));
        }
    }
    set_diag(solver);
}

void solver_preprocess_field(solver_t* solver) {
    int i;

    // Make a copy of the original x,y list.
    solver->fieldxy = starxy_copy(solver->fieldxy_orig);

    if ((solver->pixel_xscale > 0) && solver->predistort) {
        logerr("Error, can't do both pixel_xscale and predistortion at the same time!");
    }
    if (solver->pixel_xscale > 0) {
        logverb("Applying x-factor of %f to %i stars\n",
                solver->pixel_xscale, starxy_n(solver->fieldxy_orig));
        for (i=0; i<starxy_n(solver->fieldxy); i++)
            solver->fieldxy->x[i] *= solver->pixel_xscale;
    } else if (solver->predistort) {
        logverb("Applying undistortion to %i stars\n", starxy_n(solver->fieldxy_orig));
        // Apply the *un*distortion
        for (i=0; i<starxy_n(solver->fieldxy); i++) {
            double dx, dy;
            sip_pixel_undistortion(solver->predistort,
                                   solver->fieldxy->x[i], solver->fieldxy->y[i],
                                   &dx, &dy);
            solver->fieldxy->x[i] = dx;
            solver->fieldxy->y[i] = dy;
        }
    }

    find_field_boundaries(solver);
    // precompute a kdtree over the field
    solver->vf = verify_field_preprocess(solver->fieldxy);

    solver->vf->do_uniformize = solver->verify_uniformize;
    solver->vf->do_dedup = solver->verify_dedup;

    if (solver->set_crpix && solver->set_crpix_center) {
        solver->crpix[0] = wcs_pixel_center_for_size(solver_field_width(solver));
        solver->crpix[1] = wcs_pixel_center_for_size(solver_field_height(solver));
        logverb("Setting CRPIX to center (%.1f, %.1f) based on image size %i x %i\n",
                solver->crpix[0], solver->crpix[1],
                (int)solver_field_width(solver), (int)solver_field_height(solver));
    }
}

void solver_free_field(solver_t* solver) {
    solver_release_field_geometry(solver);
    if (solver->fieldxy) {
        starxy_free(solver->fieldxy);
    }
    solver->fieldxy = NULL;
    if (solver->fieldxy_orig) {
        starxy_free(solver->fieldxy_orig);
    }
    solver->fieldxy_orig = NULL;
    if (solver->vf) {
        verify_field_free(solver->vf);
    }
    solver->vf = NULL;
}

starxy_t* solver_get_field(solver_t* solver) {
    return solver->fieldxy;
}

static double get_tolerance_for_noise(
    double codetol,
    double rel_field_noise2,
    double rel_index_noise2) {
    (void)rel_field_noise2;
    (void)rel_index_noise2;
    return square(codetol);
    /*
     double maxtol2 = square(codetol);
     double tol2;
     tol2 = 49.0 * (rel_field_noise2 + rel_index_noise2);
     //printf("code tolerance %g.\n", sqrt(tol2));
     if (tol2 > maxtol2)
     tol2 = maxtol2;
     return tol2;
     */
}

static double get_tolerance(solver_t* solver) {
    return get_tolerance_for_noise(
        solver->codetol,
        solver->rel_field_noise2,
        solver->rel_index_noise2);
}

/*
 A somewhat tricky recursive function: stars A and B have already been
 chosen, so the code coordinate system has been fixed, and we've
 already determined which other stars will create valid codes (ie, are
 in the "box").  Now we want to build features using all sets of valid
 stars (without permutations).

 pq - data associated with the AB pair.
 field - the array of field star numbers
 fieldoffset - offset into the field array where we should add the first star
 n_to_add - number of stars to add
 adding - the star we're currently adding; in [0, n_to_add).
 fieldtop - the maximum field star number to build quads out of.
 dimquad, solver, tol2 - passed to try_all_codes.
 */
static void add_stars(const pquad* pq, int* field, int fieldoffset,
                      int n_to_add, int adding, int fieldtop,
                      int dimquad,
                      solver_t* solver, double tol2,
                      kdtree_qres_t** presult) {
    int bottom;
    int* f = field + fieldoffset;
    // When we're adding the first star, we start from index zero.
    // When we're adding subsequent stars, we start from the previous value
    // plus one, to avoid adding permutations.
    bottom = (adding ? f[adding-1] + 1 : 0);

    // It looks funny that we're using f[adding] as a loop variable, but
    // it's required because try_all_codes needs to know which field stars
    // were used to create the quad (which are stored in the "f" array)
    for (f[adding]=bottom; f[adding]<fieldtop; f[adding]++) {
        if (!pq->inbox[f[adding]])
            continue;
        if (unlikely(solver->quit_now))
            return;

        // If we've hit the end of the recursion (we're adding the last star),
        // call try_all_codes to try the quad we've built.
        if (adding == n_to_add-1) {
            // (when not testing, TRY_ALL_CODES is just try_all_codes.)
            TRY_ALL_CODES(pq,
                          field,
                          dimquad,
                          solver,
                          tol2,
                          presult);
        } else {
            // Else recurse.
            add_stars(pq, field, fieldoffset, n_to_add, adding+1,
                      fieldtop, dimquad, solver, tol2, presult);
        }
    }
}

#define SOLVER_AB_BLOCKS_PER_WORKER 4U
#ifndef SOLVER_AB_PACKET_LIMIT_BYTES
#define SOLVER_AB_PACKET_LIMIT_BYTES (1024U * 1024U)
#endif
#ifndef SOLVER_AB_HYPOTHESIS_LIMIT_BYTES
#define SOLVER_AB_HYPOTHESIS_LIMIT_BYTES (256U * 1024U)
#endif
#ifndef SOLVER_AB_CANDIDATE_LIMIT_BYTES
#define SOLVER_AB_CANDIDATE_LIMIT_BYTES \
    (SOLVER_AB_PACKET_LIMIT_BYTES - SOLVER_AB_HYPOTHESIS_LIMIT_BYTES)
#endif
#define SOLVER_AB_REORDER_WINDOWS_PER_WORKER 2U
#ifndef SOLVER_AB_PLAN_CACHE_BUDGET_BYTES
#define SOLVER_AB_PLAN_CACHE_BUDGET_BYTES \
    (64U * 1024U * 1024U)
#endif
#ifndef SOLVER_AB_PACKET_CACHE_BUDGET_BYTES
#define SOLVER_AB_PACKET_CACHE_BUDGET_BYTES \
    (64U * 1024U * 1024U)
#endif
#ifndef SOLVER_AB_ELIGIBLE_CACHE_BUDGET_BYTES
#define SOLVER_AB_ELIGIBLE_CACHE_BUDGET_BYTES \
    (64U * 1024U * 1024U)
#endif
#ifndef SOLVER_AB_INLINE_PUBLISH_MAX
#define SOLVER_AB_INLINE_PUBLISH_MAX 64U
#endif
#define SOLVER_AB_DESCRIPTOR_MAX_TASKS 6U
#define SOLVER_AB_DESCRIPTOR_CAPACITY 4096U
#define SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS 64U
#define SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS 1000
#define SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES \
    (8U * 1024U * 1024U)
#define SOLVER_AB_WAIT_NANOSECONDS 100000000L

typedef enum solver_ab_phase_kind {
    SOLVER_AB_PHASE_DIAGONAL = 0,
    SOLVER_AB_PHASE_OFF_DIAGONAL = 1
} solver_ab_phase_kind_t;

typedef enum solver_ab_phase_mode {
    SOLVER_AB_MODE_NATIVE = 0,
    SOLVER_AB_MODE_EMPTY = 1,
    SOLVER_AB_MODE_FLATTENED_OWNER = 2,
    SOLVER_AB_MODE_ASSISTED = 3
} solver_ab_phase_mode_t;

typedef enum solver_ab_candidate_action {
    SOLVER_AB_CANDIDATE_RADEC_SKIP = 0,
    SOLVER_AB_CANDIDATE_ABSCALE_SKIP = 1,
    SOLVER_AB_CANDIDATE_BAD_QUAD = 2,
    SOLVER_AB_CANDIDATE_SCALE_SKIP = 3,
    SOLVER_AB_CANDIDATE_VERIFY = 4
} solver_ab_candidate_action_t;

typedef enum solver_ab_task_state {
    SOLVER_AB_TASK_PENDING = 0,
    SOLVER_AB_TASK_RUNNING = 1,
    SOLVER_AB_TASK_DONE = 2
} solver_ab_task_state_t;

typedef enum solver_ab_reserve_result {
    SOLVER_AB_RESERVE_ERROR = -1,
    SOLVER_AB_RESERVE_OK = 0,
    SOLVER_AB_RESERVE_FULL = 1
} solver_ab_reserve_result_t;

typedef enum solver_ab_helper_state {
    SOLVER_AB_HELPER_NEW = 0,
    SOLVER_AB_HELPER_ACTIVE = 1,
    SOLVER_AB_HELPER_EXITED = 2,
    SOLVER_AB_HELPER_FAILED = 3
} solver_ab_helper_state_t;

typedef struct solver_ab_pair {
    int field_a;
    int field_b;
    unsigned long long combination_first;
    unsigned long long combination_count;
    double tol2;
} solver_ab_pair_t;

typedef struct solver_ab_descriptor {
    int stars[DQMAX];
    double code[DCMAX];
    double tol2;
    double rel_field_noise2;
    unsigned long long numtries_delta;
    unsigned long long cxdx_delta;
    unsigned long long meanx_delta;
    anbool current_parity;
} solver_ab_descriptor_t;

typedef struct solver_ab_descriptor_output {
    size_t descriptor_count;
    unsigned long long trailing_numtries;
    unsigned long long trailing_cxdx;
    unsigned long long trailing_meanx;
    double final_rel_field_noise2;
    anbool has_final_rel_field_noise2;
    solver_ab_descriptor_t
        descriptors[SOLVER_AB_DESCRIPTOR_CAPACITY];
} solver_ab_descriptor_output_t;

typedef struct solver_ab_candidate {
    solver_ab_candidate_action_t action;
    int quadno;
    double code_err;
    tan_t wcs;
    double scale;
    anbool parity;
    int quad_npeers;
    unsigned int star[DQMAX];
    int field[DQMAX];
    double quadpix[2 * DQMAX];
    double quadxyz[3 * DQMAX];
} solver_ab_candidate_t;

typedef struct solver_ab_hypothesis {
    unsigned long long numtries_delta;
    unsigned long long cxdx_delta;
    unsigned long long meanx_delta;
    size_t candidate_first;
    size_t candidate_count;
    int nresults;
    anbool search_failed;
    anbool begins_hypothesis;
    anbool ends_hypothesis;
    uint64_t hypothesis_order_digest;
    uint64_t kd_result_order_digest;
    double codekd_wall_seconds;
    double prepare_wall_seconds;
} solver_ab_hypothesis_t;

typedef struct solver_ab_packet {
    solver_ab_hypothesis_t* hypotheses;
    size_t hypothesis_count;
    size_t hypothesis_capacity;

    solver_ab_candidate_t* candidates;
    size_t candidate_count;
    size_t candidate_capacity;

    unsigned long long trailing_numtries;
    unsigned long long trailing_cxdx;
    unsigned long long trailing_meanx;

    size_t allocated_bytes;
    anbool allocation_failed;
    anbool evaluation_failed;
    anbool cancelled;
} solver_ab_packet_t;

typedef struct solver_ab_task {
    unsigned long long combination_first;
    unsigned long long combination_end;
    solver_ab_task_state_t state;
    int worker_id;
    anbool producer_done;
} solver_ab_task_t;

typedef struct solver_ab_channel {
    solver_ab_packet_t packet;
    size_t task_index;
    anbool ready;
    anbool final;
} solver_ab_channel_t;

typedef struct solver_ab_builder {
    solver_ab_executor_t* executor;
    solver_ab_packet_t* packet;
    kdtree_qres_t** query_result;
    size_t task_index;
    int worker_id;
    double tol2;
    double rel_field_noise2;
    unsigned long long pending_numtries;
    unsigned long long pending_cxdx;
    unsigned long long pending_meanx;
    solver_ab_descriptor_output_t* descriptor_output;
    anbool rel_field_noise_valid;
    anbool fatal_error;
    anbool producer_completed;
} solver_ab_builder_t;

typedef struct solver_ab_snapshot {
    index_t* index;
    const starxy_t* fieldxy;
    anbool use_radec;
    anbool cx_less_than_dx;
    anbool meanx_less_than_half;
    int parity;
    double centerxyz[3];
    double r2;
    double abscale_low;
    double abscale_high;
    double funits_lower;
    double funits_upper;
    double cxdx_margin;
    double codetol;
    double rel_index_noise2;
} solver_ab_snapshot_t;

typedef struct solver_ab_phase_telemetry {
    anbool enabled;
    anbool resource_valid;
    double wall_start;
    struct rusage resource_start;
    int combinations_start;
    int candidates_start;
    unsigned long long codekd_calls_start;
    unsigned long long codekd_hits_start;
    unsigned long long verify_calls_start;
    unsigned long long blocks_planned_start;
    unsigned long long blocks_retired_start;
    unsigned long long blocks_owner_start;
    unsigned long long segments_retired_start;
    unsigned long long segment_payload_bytes_start;
} solver_ab_phase_telemetry_t;

typedef struct solver_ab_reduce_state {
    anbool in_hypothesis;
    int expected_candidates;
    int seen_candidates;
} solver_ab_reduce_state_t;

struct solver_ab_executor {
    pthread_mutex_t mutex;
    pthread_cond_t work_cv;
    pthread_cond_t done_cv;

    int worker_count;
    int bound;
    int join_closed;
    int stopping;
    int fatal_error;
    int cancel_requested;
    int helpers_arrived;
    int helpers_exited;
    int helpers_failed;
    int work_notify_inflight;
    unsigned long long pending_allocation_failures;
    unsigned long long pending_search_failures;
    anbool pending_evaluation_failure;
    unsigned char* helper_state;
    void (*work_notify)(void*);
    int (*available_lenders)(void*);
    void* work_notify_opaque;
    anbool detailed;
    anbool starkd_initializing;
    anbool starkd_ready;

    unsigned long generation;
    anbool phase_active;

    solver_t* owner;
    index_t* index;
    solver_ab_snapshot_t snapshot;
    anbool snapshot_ready;
    const solver_field_geometry_t* field_geometry;
    int newpoint;
    int dimquads;
    solver_ab_phase_kind_t phase;
    double min_ab2;
    double max_ab2;

    solver_ab_pair_t* pairs;
    size_t pair_count;
    size_t pair_capacity;
    solver_ab_pair_t* pair_cache;
    size_t pair_cache_capacity;
    size_t pair_cache_limit_bytes;
    anbool pairs_transient;
    solver_ab_task_t* tasks;
    size_t task_count;
    size_t task_capacity;
    solver_ab_task_t* task_cache;
    size_t task_cache_capacity;
    size_t task_cache_limit_bytes;
    size_t packet_cache_limit_bytes;
    size_t eligible_cache_limit_bytes;
    anbool tasks_transient;
    unsigned long long pair_cache_grows;
    unsigned long long pair_cache_reuses;
    unsigned long long task_cache_grows;
    unsigned long long task_cache_reuses;
    unsigned long long plan_transient_allocations;
    size_t plan_transient_peak_bytes;
    unsigned long long index_epochs_started;
    unsigned long long index_epochs_quiesced;
    unsigned long long owner_query_reuses;
    unsigned long long packet_cache_trims;
    unsigned long long eligible_cache_trims;
    size_t next_task;
    size_t next_reduce;
    size_t reorder_window;
    size_t remaining_tasks;

    kdtree_qres_t** query_results;
    int** combination_eligible;
    size_t* combination_eligible_capacity;
    solver_ab_channel_t* channels;
    solver_ab_packet_t drain_packet;
};

static void solver_ab_executor_notify_work(
    solver_ab_executor_t* executor);

/*
 * startree_get() lazily builds the index-owned inverse permutation on its
 * first real StarKD lookup. Parallel AB evaluation must retain that laziness
 * while ensuring that exactly one worker performs the construction.
 *
 * Initialization is scientific-state neutral and publishes an immutable
 * vector. The first real hit may therefore initialize it from any canonical
 * work packet; retirement order continues to be enforced independently.
 */
static int solver_ab_ensure_starkd_ready(
    solver_ab_executor_t* executor,
    int worker_id,
    size_t task_index) {
    startree_t* starkd;
    const char* index_name = NULL;
    int rc = 0;

    if (!executor ||
        worker_id < 0 ||
        worker_id >= executor->worker_count) {
        return -1;
    }

    pthread_mutex_lock(&executor->mutex);
    while (TRUE) {
        if (executor->stopping ||
            __atomic_load_n(
                &executor->cancel_requested,
                __ATOMIC_ACQUIRE)) {
            pthread_mutex_unlock(&executor->mutex);
            return 1;
        }
        if (executor->fatal_error ||
            !executor->bound ||
            !executor->index ||
            !executor->index->starkd ||
            !executor->index->starkd->tree) {
            pthread_mutex_unlock(&executor->mutex);
            return -1;
        }
        starkd = executor->index->starkd;
        if (executor->starkd_ready) {
            pthread_mutex_unlock(&executor->mutex);
            return 0;
        }
        if (!starkd->tree->perm ||
            starkd->inverse_perm) {
            executor->starkd_ready = TRUE;
            pthread_mutex_unlock(&executor->mutex);
            return 0;
        }
        if (executor->starkd_initializing) {
            pthread_cond_wait(
                &executor->done_cv,
                &executor->mutex);
            continue;
        }
        executor->starkd_initializing = TRUE;
        index_name = executor->index->indexname;
        break;
    }
    pthread_mutex_unlock(&executor->mutex);

    logverb(
        "[solver-ab] starkd-inverse-init state=start "
        "index=%s worker=%i task=%zu\n",
        index_name ? index_name : "(unnamed)",
        worker_id,
        task_index);

    /* The universal StarKD lifecycle owns exact-PERM-range advice. */
    startree_compute_inverse_perm(starkd);
    if (!starkd->inverse_perm) {
        rc = -1;
    }

    logverb(
        "[solver-ab] starkd-inverse-init state=%s "
        "index=%s worker=%i task=%zu\n",
        rc ? "failed" : "ready",
        index_name ? index_name : "(unnamed)",
        worker_id,
        task_index);

    pthread_mutex_lock(&executor->mutex);
    executor->starkd_initializing = FALSE;
    if (rc) {
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
    } else {
        executor->starkd_ready = TRUE;
    }
    pthread_cond_broadcast(&executor->done_cv);
    pthread_cond_broadcast(&executor->work_cv);
    pthread_mutex_unlock(&executor->mutex);
    return rc;
}

static int solver_ab_capture_snapshot(
    solver_ab_executor_t* executor,
    const solver_t* solver) {
    solver_ab_snapshot_t* snapshot;

    if (!executor || !solver || !executor->index ||
        !solver->fieldxy) {
        return -1;
    }
    if (executor->snapshot_ready) {
        return 0;
    }

    snapshot = &executor->snapshot;
    memset(snapshot, 0, sizeof(*snapshot));
    snapshot->index = executor->index;
    snapshot->fieldxy = solver->fieldxy;
    snapshot->use_radec = solver->use_radec;
    snapshot->cx_less_than_dx =
        executor->index->cx_less_than_dx;
    snapshot->meanx_less_than_half =
        executor->index->meanx_less_than_half;
    snapshot->parity = solver->parity;
    memcpy(
        snapshot->centerxyz,
        solver->centerxyz,
        sizeof(snapshot->centerxyz));
    snapshot->r2 = solver->r2;
    snapshot->abscale_low = solver->abscale_low;
    snapshot->abscale_high = solver->abscale_high;
    snapshot->funits_lower = solver->funits_lower;
    snapshot->funits_upper = solver->funits_upper;
    snapshot->cxdx_margin = solver->cxdx_margin;
    snapshot->codetol = solver->codetol;
    snapshot->rel_index_noise2 = solver->rel_index_noise2;
    executor->snapshot_ready = TRUE;
    return 0;
}

static uint64_t solver_order_hash_mix(
    uint64_t state,
    uint64_t value) {
    int byte_index;

    if (!state) {
        state = UINT64_C(1469598103934665603);
    }
    for (byte_index = 0; byte_index < 8; byte_index++) {
        state ^=
            (value >> (8 * byte_index)) & UINT64_C(0xff);
        state *= UINT64_C(1099511628211);
    }
    return state;
}

static uint64_t solver_order_double_bits(double value) {
    uint64_t bits;

    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static uint64_t solver_order_float_bits(float value) {
    uint32_t bits;

    memcpy(&bits, &value, sizeof(bits));
    return (uint64_t)bits;
}

static uint64_t solver_hypothesis_order_digest(
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    uint64_t digest = 0;
    int dimcode = (dimquad - NBACK) * 2;
    int i;

    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x4859504f54484553));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)dimquad);
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)current_parity);
    for (i = 0; i < dimquad; i++) {
        digest = solver_order_hash_mix(
            digest,
            (uint64_t)(unsigned int)stars[i]);
    }
    for (i = 0; i < dimcode; i++) {
        digest = solver_order_hash_mix(
            digest,
            solver_order_double_bits(code[i]));
    }
    return digest;
}

static uint64_t solver_kd_result_order_digest(
    const kdtree_qres_t* result) {
    uint64_t digest = 0;
    int i;

    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x4b44524553554c54));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)result->nres);
    for (i = 0; i < result->nres; i++) {
        digest = solver_order_hash_mix(
            digest,
            (uint64_t)(unsigned int)result->inds[i]);
        digest = solver_order_hash_mix(
            digest,
            solver_order_double_bits(result->sdists[i]));
    }
    return digest;
}

static void solver_record_candidate_order(
    solver_t* solver,
    solver_ab_candidate_action_t action,
    int quadno,
    float code_err) {
    uint64_t digest = 0;

    if (!solver->profile.detailed) {
        return;
    }
    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x43414e4449444154));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)action);
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)quadno);
    digest = solver_order_hash_mix(
        digest,
        solver_order_float_bits(code_err));
    solver->profile.candidate_order_hash =
        solver_order_hash_mix(
            solver->profile.candidate_order_hash,
            digest);
}

static double solver_ab_timeval_seconds(
    const struct timeval* value) {
    return (double)value->tv_sec +
        (double)value->tv_usec * 1.0e-6;
}

static void solver_ab_phase_telemetry_begin(
    const solver_t* solver,
    solver_ab_phase_telemetry_t* telemetry) {
    memset(telemetry, 0, sizeof(*telemetry));
    if (!solver ||
        log_get_level() < LOG_VERB ||
        pl_size(solver->indexes) != 1U) {
        return;
    }
    telemetry->enabled = TRUE;
    telemetry->wall_start = monotonic_seconds();
    telemetry->resource_valid =
        getrusage(
            RUSAGE_SELF,
            &telemetry->resource_start) == 0;
    telemetry->combinations_start = solver->numtries;
    telemetry->candidates_start = solver->nummatches;
    telemetry->codekd_calls_start =
        solver->profile.codekd_calls;
    telemetry->codekd_hits_start =
        solver->profile.codekd_hits;
    telemetry->verify_calls_start =
        solver->profile.verify_calls;
    telemetry->blocks_planned_start =
        solver->profile.ab_blocks_planned;
    telemetry->blocks_retired_start =
        solver->profile.ab_blocks_retired;
    telemetry->blocks_owner_start =
        solver->profile.ab_blocks_owner;
    telemetry->segments_retired_start =
        solver->profile.ab_segments_retired;
    telemetry->segment_payload_bytes_start =
        solver->profile.ab_segment_payload_bytes;
}

static void solver_ab_phase_telemetry_report(
    const solver_t* solver,
    const solver_ab_phase_telemetry_t* telemetry,
    int newpoint,
    solver_ab_phase_kind_t phase,
    solver_ab_phase_mode_t mode) {
    struct rusage resource_end;
    anbool resource_valid;
    double user_seconds = 0.0;
    double system_seconds = 0.0;
    long major_faults = 0;
    const char* index_name;

    if (!solver || !telemetry || !telemetry->enabled) {
        return;
    }
    resource_valid =
        telemetry->resource_valid &&
        getrusage(RUSAGE_SELF, &resource_end) == 0;
    if (resource_valid) {
        user_seconds =
            solver_ab_timeval_seconds(&resource_end.ru_utime) -
            solver_ab_timeval_seconds(
                &telemetry->resource_start.ru_utime);
        system_seconds =
            solver_ab_timeval_seconds(&resource_end.ru_stime) -
            solver_ab_timeval_seconds(
                &telemetry->resource_start.ru_stime);
        major_faults =
            resource_end.ru_majflt -
            telemetry->resource_start.ru_majflt;
    }
    index_name =
        solver->index && solver->index->indexname
            ? solver->index->indexname
            : "(unknown)";
    logverb("[solver-ab-phase] index=%s object=%i phase=%s mode=%s "
            "combinations=%llu codekd_queries=%llu codekd_hits=%llu "
            "candidates=%llu verifications=%llu "
            "blocks_planned=%llu blocks_retired=%llu "
            "blocks_owner=%llu segments_retired=%llu "
            "segment_payload_bytes=%llu "
            "wall=%.6f user=%.6f system=%.6f "
            "major_faults=%ld resource=%s "
            "hypothesis_order=%016llx kd_result_order=%016llx "
            "candidate_order=%016llx\n",
            index_name,
            newpoint + 1,
            phase == SOLVER_AB_PHASE_DIAGONAL
                ? "diagonal"
                : "off-diagonal",
            mode == SOLVER_AB_MODE_ASSISTED
                ? "assisted"
                : (mode == SOLVER_AB_MODE_FLATTENED_OWNER
                    ? "flattened-owner"
                    : (mode == SOLVER_AB_MODE_EMPTY
                        ? "empty"
                        : "native")),
            (unsigned long long)(
                solver->numtries -
                telemetry->combinations_start),
            solver->profile.codekd_calls -
                telemetry->codekd_calls_start,
            solver->profile.codekd_hits -
                telemetry->codekd_hits_start,
            (unsigned long long)(
                solver->nummatches -
                telemetry->candidates_start),
            solver->profile.verify_calls -
                telemetry->verify_calls_start,
            solver->profile.ab_blocks_planned -
                telemetry->blocks_planned_start,
            solver->profile.ab_blocks_retired -
                telemetry->blocks_retired_start,
            solver->profile.ab_blocks_owner -
                telemetry->blocks_owner_start,
            solver->profile.ab_segments_retired -
                telemetry->segments_retired_start,
            solver->profile.ab_segment_payload_bytes -
                telemetry->segment_payload_bytes_start,
            monotonic_seconds() - telemetry->wall_start,
            user_seconds,
            system_seconds,
            major_faults,
            resource_valid
                ? "process-overlap"
                : "unavailable",
            solver->profile.hypothesis_order_hash,
            solver->profile.kd_result_order_hash,
            solver->profile.candidate_order_hash);
}

static unsigned long long solver_ab_saturating_add(
    unsigned long long a,
    unsigned long long b) {
    if (ULLONG_MAX - a < b) {
        return ULLONG_MAX;
    }
    return a + b;
}

static unsigned long long solver_ab_saturating_choose(int n, int k) {
    unsigned long long result = 1;
    int i;

    if (n < 0 || k < 0 || k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    for (i = 1; i <= k; i++) {
        unsigned long long numerator =
            (unsigned long long)(n - k + i);

        if (result > ULLONG_MAX / numerator) {
            return ULLONG_MAX;
        }
        result *= numerator;
        result /= (unsigned long long)i;
    }
    return result;
}

static solver_ab_reserve_result_t
solver_ab_packet_reserve_hypotheses(
    solver_ab_packet_t* packet,
    size_t required) {
    solver_ab_hypothesis_t* resized;
    size_t capacity;
    size_t max_capacity =
        SOLVER_AB_HYPOTHESIS_LIMIT_BYTES / sizeof(*resized);
    size_t bytes;

    if (required <= packet->hypothesis_capacity) {
        return SOLVER_AB_RESERVE_OK;
    }
    if (required > max_capacity || !max_capacity) {
        return SOLVER_AB_RESERVE_FULL;
    }
    capacity = packet->hypothesis_capacity ?
        packet->hypothesis_capacity : MIN(16U, max_capacity);
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            packet->allocation_failed = TRUE;
            return SOLVER_AB_RESERVE_ERROR;
        }
        if (capacity > max_capacity / 2U) {
            capacity = max_capacity;
        } else {
            capacity *= 2U;
        }
    }
    if (capacity > SIZE_MAX / sizeof(*resized)) {
        packet->allocation_failed = TRUE;
        return SOLVER_AB_RESERVE_ERROR;
    }
    bytes = capacity * sizeof(*resized);
    if (bytes > SOLVER_AB_HYPOTHESIS_LIMIT_BYTES) {
        return SOLVER_AB_RESERVE_FULL;
    }

    resized = realloc(packet->hypotheses, bytes);
    if (!resized) {
        packet->allocation_failed = TRUE;
        return SOLVER_AB_RESERVE_ERROR;
    }
    packet->hypotheses = resized;
    packet->hypothesis_capacity = capacity;
    packet->allocated_bytes =
        bytes +
        packet->candidate_capacity * sizeof(*packet->candidates);
    return SOLVER_AB_RESERVE_OK;
}

static solver_ab_reserve_result_t
solver_ab_packet_reserve_candidates(
    solver_ab_packet_t* packet,
    size_t required) {
    solver_ab_candidate_t* resized;
    size_t capacity;
    size_t max_capacity =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES / sizeof(*resized);
    size_t bytes;

    if (required <= packet->candidate_capacity) {
        return SOLVER_AB_RESERVE_OK;
    }
    if (required > max_capacity || !max_capacity) {
        return SOLVER_AB_RESERVE_FULL;
    }
    capacity = packet->candidate_capacity ?
        packet->candidate_capacity : MIN(16U, max_capacity);
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            packet->allocation_failed = TRUE;
            return SOLVER_AB_RESERVE_ERROR;
        }
        if (capacity > max_capacity / 2U) {
            capacity = max_capacity;
        } else {
            capacity *= 2U;
        }
    }
    if (capacity > SIZE_MAX / sizeof(*resized)) {
        packet->allocation_failed = TRUE;
        return SOLVER_AB_RESERVE_ERROR;
    }
    bytes = capacity * sizeof(*resized);
    if (bytes > SOLVER_AB_CANDIDATE_LIMIT_BYTES) {
        return SOLVER_AB_RESERVE_FULL;
    }

    resized = realloc(packet->candidates, bytes);
    if (!resized) {
        packet->allocation_failed = TRUE;
        return SOLVER_AB_RESERVE_ERROR;
    }
    packet->candidates = resized;
    packet->candidate_capacity = capacity;
    packet->allocated_bytes =
        packet->hypothesis_capacity * sizeof(*packet->hypotheses) +
        bytes;
    return SOLVER_AB_RESERVE_OK;
}

static void solver_ab_packet_reset(
    solver_ab_packet_t* packet) {
    packet->hypothesis_count = 0U;
    packet->candidate_count = 0U;
    packet->trailing_numtries = 0U;
    packet->trailing_cxdx = 0U;
    packet->trailing_meanx = 0U;
    packet->allocation_failed = FALSE;
    packet->evaluation_failed = FALSE;
    packet->cancelled = FALSE;
}

static anbool solver_ab_packet_has_events(
    const solver_ab_packet_t* packet) {
    return packet->hypothesis_count > 0U ||
        packet->trailing_numtries > 0U ||
        packet->trailing_cxdx > 0U ||
        packet->trailing_meanx > 0U ||
        packet->allocation_failed ||
        packet->evaluation_failed ||
        packet->cancelled;
}

static void solver_ab_packet_free(solver_ab_packet_t* packet) {
    if (!packet) {
        return;
    }
    free(packet->hypotheses);
    free(packet->candidates);
    memset(packet, 0, sizeof(*packet));
}

static anbool solver_ab_cancelled(
    const solver_ab_executor_t* executor) {
    if (__atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE)) {
        return TRUE;
    }
    return index_shard_worker_stop_requested();
}

static void solver_ab_complete_task_locked(
    solver_ab_executor_t* executor,
    size_t task_index) {
    solver_ab_task_t* task = &executor->tasks[task_index];

    if (task->producer_done) {
        return;
    }
    task->producer_done = TRUE;
    task->state = SOLVER_AB_TASK_DONE;
    if (executor->remaining_tasks > 0U) {
        executor->remaining_tasks--;
    } else {
        logerr(
            "[solver-ab] task completion underflow for block %zu\n",
            task_index);
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
    }
}

static void solver_ab_stage_pending_deltas(
    solver_ab_builder_t* builder) {
    builder->packet->trailing_numtries +=
        builder->pending_numtries;
    builder->packet->trailing_cxdx +=
        builder->pending_cxdx;
    builder->packet->trailing_meanx +=
        builder->pending_meanx;
    builder->pending_numtries = 0U;
    builder->pending_cxdx = 0U;
    builder->pending_meanx = 0U;
}

static int solver_ab_publish_segment(
    solver_ab_builder_t* builder,
    anbool final) {
    solver_ab_executor_t* executor = builder->executor;
    solver_ab_channel_t* channel =
        &executor->channels[builder->worker_id];

    solver_ab_stage_pending_deltas(builder);
    if (!final &&
        !solver_ab_packet_has_events(builder->packet)) {
        return 0;
    }

    pthread_mutex_lock(&executor->mutex);
    if (channel->ready ||
        channel->task_index != builder->task_index ||
        executor->tasks[builder->task_index].state !=
            SOLVER_AB_TASK_RUNNING) {
        logerr(
            "[solver-ab] invalid publish state for block %zu "
            "worker %i\n",
            builder->task_index,
            builder->worker_id);
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        solver_ab_complete_task_locked(
            executor,
            builder->task_index);
        builder->producer_completed = TRUE;
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }

    channel->final = final;
    channel->ready = TRUE;
    pthread_cond_broadcast(&executor->done_cv);
    while (channel->ready &&
           !executor->stopping &&
           !__atomic_load_n(
               &executor->cancel_requested,
               __ATOMIC_ACQUIRE)) {
        pthread_cond_wait(
            &executor->work_cv,
            &executor->mutex);
    }

    if (channel->ready) {
        channel->ready = FALSE;
        channel->final = FALSE;
        solver_ab_packet_reset(&channel->packet);
    }
    if (executor->stopping ||
        __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE)) {
        solver_ab_complete_task_locked(
            executor,
            builder->task_index);
        builder->producer_completed = TRUE;
        pthread_cond_broadcast(&executor->done_cv);
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }
    if (final) {
        solver_ab_complete_task_locked(
            executor,
            builder->task_index);
        builder->producer_completed = TRUE;
        pthread_cond_broadcast(&executor->done_cv);
    }
    pthread_mutex_unlock(&executor->mutex);
    return 0;
}

static int solver_ab_candidate_prepare(
    solver_ab_candidate_t* candidate,
    const kdtree_qres_t* result,
    int result_index,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    const solver_ab_snapshot_t* snapshot,
    anbool current_parity) {
    double starxyz[DQMAX * 3];
    double scale;
    double arcsecperpix;
    double abscale;
    tan_t wcs;
    unsigned int star[DQMAX];
    int thisquadno;
    int i;
    anbool outofbounds = FALSE;

    memset(candidate, 0, sizeof(*candidate));
    candidate->action = SOLVER_AB_CANDIDATE_SCALE_SKIP;

    thisquadno = result->inds[result_index];
    candidate->quadno = thisquadno;
    candidate->code_err = result->sdists[result_index];
    if (quadfile_get_stars(
            snapshot->index->quads,
            thisquadno,
            star)) {
        return -1;
    }

    if (snapshot->use_radec) {
        for (i = 0; i < dimquads; i++) {
            if (startree_get(
                    snapshot->index->starkd,
                    star[i],
                    starxyz + 3 * i)) {
                return -1;
            }
            if (distsq(starxyz + 3 * i,
                       snapshot->centerxyz,
                       3) > snapshot->r2) {
                outofbounds = TRUE;
                break;
            }
        }
        if (outofbounds) {
            candidate->action =
                SOLVER_AB_CANDIDATE_RADEC_SKIP;
            return 0;
        }
    } else {
        if (startree_get(
                snapshot->index->starkd,
                star[0],
                starxyz) ||
            startree_get(
                snapshot->index->starkd,
                star[1],
                starxyz + 3)) {
            return -1;
        }
    }

    abscale =
        square(distsq2rad(distsq(starxyz, starxyz + 3, 3))) /
        distsq(field_xy, field_xy + 2, 2);
    if (abscale > snapshot->abscale_high ||
        abscale < snapshot->abscale_low) {
        candidate->action =
            SOLVER_AB_CANDIDATE_ABSCALE_SKIP;
        return 0;
    }

    if (!snapshot->use_radec) {
        for (i = 2; i < dimquads; i++) {
            if (startree_get(
                    snapshot->index->starkd,
                    star[i],
                    starxyz + 3 * i)) {
                return -1;
            }
        }
    }

    if (fit_tan_wcs(
            starxyz,
            field_xy,
            dimquads,
            &wcs,
            &scale)) {
        candidate->action = SOLVER_AB_CANDIDATE_BAD_QUAD;
        return 0;
    }
    arcsecperpix = scale * 3600.0;
    if (arcsecperpix > snapshot->funits_upper ||
        arcsecperpix < snapshot->funits_lower) {
        candidate->action = SOLVER_AB_CANDIDATE_SCALE_SKIP;
        return 0;
    }

    memcpy(&candidate->wcs, &wcs, sizeof(tan_t));
    candidate->scale = arcsecperpix;
    candidate->parity = current_parity;
    candidate->quad_npeers = result->nres;
    for (i = 0; i < dimquads; i++) {
        candidate->star[i] = star[i];
        candidate->field[i] = fieldstars[i];
    }
    memcpy(candidate->quadpix,
           field_xy,
           (size_t)2 * (size_t)dimquads * sizeof(double));
    memcpy(candidate->quadxyz,
           starxyz,
           (size_t)3 * (size_t)dimquads * sizeof(double));
    candidate->action = SOLVER_AB_CANDIDATE_VERIFY;
    return 0;
}

static int solver_ab_record_descriptor(
    solver_ab_builder_t* builder,
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    solver_ab_descriptor_output_t* output =
        builder->descriptor_output;
    solver_ab_descriptor_t* descriptor;
    int dimcode = (dimquad - NBACK) * 2;

    if (!output || !stars || !code ||
        dimquad < NBACK || dimquad > DQMAX ||
        dimcode < 0 || dimcode > DCMAX ||
        output->descriptor_count >=
            SOLVER_AB_DESCRIPTOR_CAPACITY) {
        builder->fatal_error = TRUE;
        if (builder->packet) {
            builder->packet->evaluation_failed = TRUE;
        }
        return -1;
    }
    descriptor = &output->descriptors[
        output->descriptor_count++];
    memcpy(
        descriptor->stars,
        stars,
        (size_t)dimquad * sizeof(*stars));
    memcpy(
        descriptor->code,
        code,
        (size_t)dimcode * sizeof(*code));
    descriptor->tol2 = builder->tol2;
    descriptor->rel_field_noise2 =
        builder->rel_field_noise2;
    descriptor->numtries_delta =
        builder->pending_numtries;
    descriptor->cxdx_delta = builder->pending_cxdx;
    descriptor->meanx_delta = builder->pending_meanx;
    descriptor->current_parity = current_parity;
    builder->pending_numtries = 0U;
    builder->pending_cxdx = 0U;
    builder->pending_meanx = 0U;
    return 0;
}

static int solver_ab_record_hypothesis(
    solver_ab_builder_t* builder,
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    solver_ab_executor_t* executor = builder->executor;
    solver_ab_packet_t* packet = builder->packet;
    solver_ab_hypothesis_t* hypothesis;
    double field_xy[DQMAX * 2];
    double search_start = 0.0;
    double prepare_start = 0.0;
    solver_ab_reserve_result_t reserve_result;
    int nresults;
    int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;
    int i;

    if (solver_ab_cancelled(executor)) {
        packet->cancelled = TRUE;
        return -1;
    }
    if (builder->descriptor_output) {
        return solver_ab_record_descriptor(
            builder,
            stars,
            code,
            dimquad,
            current_parity);
    }
    reserve_result = solver_ab_packet_reserve_hypotheses(
        packet,
        packet->hypothesis_count + 1U);
    if (reserve_result == SOLVER_AB_RESERVE_FULL) {
        if (solver_ab_publish_segment(builder, FALSE)) {
            return -1;
        }
        reserve_result = solver_ab_packet_reserve_hypotheses(
            packet,
            packet->hypothesis_count + 1U);
    }
    if (reserve_result != SOLVER_AB_RESERVE_OK) {
        packet->allocation_failed = TRUE;
        builder->fatal_error = TRUE;
        return -1;
    }

    hypothesis = &packet->hypotheses[
        packet->hypothesis_count++];
    memset(hypothesis, 0, sizeof(*hypothesis));
    hypothesis->begins_hypothesis = TRUE;
    hypothesis->numtries_delta = builder->pending_numtries;
    hypothesis->cxdx_delta = builder->pending_cxdx;
    hypothesis->meanx_delta = builder->pending_meanx;
    builder->pending_numtries = 0;
    builder->pending_cxdx = 0;
    builder->pending_meanx = 0;
    if (executor->detailed) {
        hypothesis->hypothesis_order_digest =
            solver_hypothesis_order_digest(
                stars,
                code,
                dimquad,
                current_parity);
    }

#if defined(TESTING_TRYPERMUTATIONS)
    /*
     * Route the exact-order unit test through the real packet-producing
     * evaluator without requiring a synthetic CodeKD. Production builds
     * continue directly into the native range search below.
     */
    TEST_TRY_PERMUTATIONS(
        (int*)stars,
        (double*)code,
        dimquad,
        executor->owner);
    hypothesis->ends_hypothesis = TRUE;
    return 0;
#endif

    if (executor->detailed) {
        search_start = monotonic_seconds();
    }
    *builder->query_result =
        solver_codekd_rangesearch(
            executor->snapshot.index->codekd->tree,
            *builder->query_result,
            code,
            builder->tol2,
            options);
    if (executor->detailed) {
        hypothesis->codekd_wall_seconds =
            monotonic_seconds() - search_start;
    }

    if (!*builder->query_result) {
        hypothesis->search_failed = TRUE;
        hypothesis->ends_hypothesis = TRUE;
        builder->fatal_error = TRUE;
        return -1;
    }
    if (executor->detailed) {
        hypothesis->kd_result_order_digest =
            solver_kd_result_order_digest(
                *builder->query_result);
    }
    nresults = (*builder->query_result)->nres;
    hypothesis->nresults = nresults;
    hypothesis->candidate_first = packet->candidate_count;

    if (!nresults) {
        hypothesis->ends_hypothesis = TRUE;
        return 0;
    }
    {
        int starkd_status =
            solver_ab_ensure_starkd_ready(
                executor,
                builder->worker_id,
                builder->task_index);

        if (starkd_status > 0) {
            packet->cancelled = TRUE;
            hypothesis->ends_hypothesis = TRUE;
            return -1;
        }
        if (starkd_status < 0) {
            logerr("[solver-ab] failed to initialize StarKD lookup state\n");
            packet->allocation_failed = TRUE;
            hypothesis->ends_hypothesis = TRUE;
            builder->fatal_error = TRUE;
            return -1;
        }
    }

    for (i = 0; i < dimquad; i++) {
        setx(field_xy, i,
             starxy_getx(executor->snapshot.fieldxy, stars[i]));
        sety(field_xy, i,
             starxy_gety(executor->snapshot.fieldxy, stars[i]));
    }
    if (executor->detailed) {
        prepare_start = monotonic_seconds();
    }
    for (i = 0; i < nresults; i++) {
        solver_ab_candidate_t* candidate;

        if (solver_ab_cancelled(executor)) {
            packet->cancelled = TRUE;
            return -1;
        }
        reserve_result = solver_ab_packet_reserve_candidates(
            packet,
            packet->candidate_count + 1U);
        if (reserve_result == SOLVER_AB_RESERVE_FULL) {
            if (executor->detailed) {
                hypothesis->prepare_wall_seconds +=
                    monotonic_seconds() - prepare_start;
            }
            hypothesis->candidate_count =
                packet->candidate_count -
                hypothesis->candidate_first;
            if (solver_ab_publish_segment(builder, FALSE)) {
                return -1;
            }
            reserve_result =
                solver_ab_packet_reserve_hypotheses(
                    packet,
                    1U);
            if (reserve_result != SOLVER_AB_RESERVE_OK) {
                packet->allocation_failed = TRUE;
                builder->fatal_error = TRUE;
                return -1;
            }
            hypothesis =
                &packet->hypotheses[packet->hypothesis_count++];
            memset(hypothesis, 0, sizeof(*hypothesis));
            hypothesis->candidate_first =
                packet->candidate_count;
            if (executor->detailed) {
                prepare_start = monotonic_seconds();
            }
            reserve_result =
                solver_ab_packet_reserve_candidates(
                    packet,
                    packet->candidate_count + 1U);
        }
        if (reserve_result != SOLVER_AB_RESERVE_OK) {
            packet->allocation_failed = TRUE;
            builder->fatal_error = TRUE;
            return -1;
        }
        candidate =
            &packet->candidates[packet->candidate_count++];
        if (solver_ab_candidate_prepare(
                candidate,
                *builder->query_result,
                i,
                field_xy,
                stars,
                dimquad,
                &executor->snapshot,
                current_parity)) {
            packet->evaluation_failed = TRUE;
            builder->fatal_error = TRUE;
            return -1;
        }
    }
    if (executor->detailed) {
        hypothesis->prepare_wall_seconds +=
            monotonic_seconds() - prepare_start;
    }
    hypothesis->candidate_count =
        packet->candidate_count - hypothesis->candidate_first;
    hypothesis->ends_hypothesis = TRUE;
    return 0;
}

static void solver_ab_try_permutations(
    const int* origstars,
    int dimquad,
    const double* origcode,
    solver_ab_builder_t* builder,
    anbool current_parity,
    int* stars,
    double* code,
    int slot,
    anbool* placed) {
    const solver_ab_snapshot_t* snapshot =
        &builder->executor->snapshot;
    double mycode[DCMAX];
    int nstars = dimquad - NBACK;
    int lastslot = dimquad - NBACK - 1;
    int i;

    if (!code) {
        code = mycode;
    }
    if (slot >= DCMAX / 2 ||
        solver_ab_cancelled(builder->executor) ||
        builder->fatal_error) {
        if (solver_ab_cancelled(builder->executor)) {
            builder->packet->cancelled = TRUE;
        }
        return;
    }

    for (i = 0; i < nstars; i++) {
        if (placed[i]) {
            continue;
        }
        if (slot > 0 &&
            snapshot->cx_less_than_dx &&
            code[2 * (slot - 1)] >
                origcode[2 * i] + snapshot->cxdx_margin) {
            builder->pending_cxdx++;
            continue;
        }

        stars[slot + NBACK] = origstars[i + NBACK];
        code[2 * slot] = origcode[2 * i];
        code[2 * slot + 1] = origcode[2 * i + 1];

        if (snapshot->cx_less_than_dx &&
            snapshot->meanx_less_than_half) {
            double meanx = 0.0;
            int j;

            for (j = 0; j <= slot; j++) {
                meanx += code[2 * j];
            }
            meanx /= (double)(slot + 1);
            if (meanx > 0.5 + snapshot->cxdx_margin) {
                builder->pending_meanx++;
                continue;
            }
        }

        if (slot < lastslot) {
            placed[i] = TRUE;
            solver_ab_try_permutations(
                origstars,
                dimquad,
                origcode,
                builder,
                current_parity,
                stars,
                code,
                slot + 1,
                placed);
            placed[i] = FALSE;
            if (builder->fatal_error ||
                builder->packet->cancelled) {
                return;
            }
        } else if (solver_ab_record_hypothesis(
                       builder,
                       stars,
                       code,
                       dimquad,
                       current_parity)) {
            return;
        }
    }
}

static void solver_ab_try_all_codes_2(
    const int* fieldstars,
    int dimquad,
    const double* code,
    solver_ab_builder_t* builder,
    anbool current_parity) {
    double flipcode[DCMAX];
    int stars[DQMAX];
    anbool placed[DQMAX];
    int dimcode = (dimquad - NBACK) * 2;
    int i;

    stars[0] = fieldstars[0];
    stars[1] = fieldstars[1];
    memset(placed, 0, sizeof(placed));
    solver_ab_try_permutations(
        fieldstars,
        dimquad,
        code,
        builder,
        current_parity,
        stars,
        NULL,
        0,
        placed);
    if (builder->fatal_error ||
        builder->packet->cancelled) {
        return;
    }

    stars[0] = fieldstars[1];
    stars[1] = fieldstars[0];
    for (i = 0; i < dimcode; i++) {
        flipcode[i] = 1.0 - code[i];
    }
    memset(placed, 0, sizeof(placed));
    solver_ab_try_permutations(
        fieldstars,
        dimquad,
        flipcode,
        builder,
        current_parity,
        stars,
        NULL,
        0,
        placed);
}

static void solver_ab_try_all_codes(
    const solver_pair_geometry_t* pair_geometry,
    int field_a,
    int field_b,
    const int* fieldstars,
    int dimquad,
    solver_ab_builder_t* builder) {
    const solver_field_geometry_t* geometry =
        builder->executor->field_geometry;
    const solver_ab_snapshot_t* snapshot =
        &builder->executor->snapshot;
    double code[DCMAX];
    double flipcode[DCMAX];
    int dimcode = (dimquad - 2) * 2;
    int i;

    builder->pending_numtries++;
    for (i = 0; i < dimquad - NBACK; i++) {
        if (!solver_pair_geometry_transform(
                geometry,
                pair_geometry,
                field_a,
                field_b,
                fieldstars[NBACK + i],
                &code[2 * i],
                &code[2 * i + 1])) {
            builder->fatal_error = TRUE;
            builder->packet->evaluation_failed = TRUE;
            return;
        }
    }

    if (snapshot->parity == PARITY_NORMAL ||
        snapshot->parity == PARITY_BOTH) {
        solver_ab_try_all_codes_2(
            fieldstars,
            dimquad,
            code,
            builder,
            FALSE);
    }
    if (builder->fatal_error ||
        builder->packet->cancelled) {
        return;
    }
    if (snapshot->parity == PARITY_FLIP ||
        snapshot->parity == PARITY_BOTH) {
        quad_flip_parity(code, flipcode, dimcode);
        solver_ab_try_all_codes_2(
            fieldstars,
            dimquad,
            flipcode,
            builder,
            TRUE);
    }
}

typedef anbool (*solver_ab_combination_visitor_t)(
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    int* field,
    int dimquad,
    void* opaque);

static int solver_ab_select_combination(
    const int* eligible,
    int eligible_count,
    int* field,
    int fieldoffset,
    int n_to_add,
    int* selection,
    unsigned long long ordinal) {
    int position;

    for (position = 0; position < n_to_add; position++) {
        int remaining = n_to_add - position - 1;
        int bottom = position ?
            selection[position - 1] + 1 : 0;
        int candidate_index;
        anbool selected = FALSE;

        for (candidate_index = bottom;
             candidate_index < eligible_count;
             candidate_index++) {
            unsigned long long suffix_count;

            suffix_count = solver_ab_saturating_choose(
                eligible_count - candidate_index - 1,
                remaining);
            if (ordinal >= suffix_count) {
                ordinal -= suffix_count;
                continue;
            }
            selection[position] = candidate_index;
            field[fieldoffset + position] =
                eligible[candidate_index];
            selected = TRUE;
            break;
        }
        if (!selected) {
            return -1;
        }
    }
    return ordinal == 0U ? 0 : -1;
}

static anbool solver_ab_next_combination(
    const int* eligible,
    int eligible_count,
    int* field,
    int fieldoffset,
    int n_to_add,
    int* selection) {
    int position;

    for (position = n_to_add - 1;
         position >= 0;
         position--) {
        int maximum =
            eligible_count - (n_to_add - position);
        int fill;

        if (selection[position] >= maximum) {
            continue;
        }
        selection[position]++;
        field[fieldoffset + position] =
            eligible[selection[position]];
        for (fill = position + 1;
             fill < n_to_add;
             fill++) {
            selection[fill] = selection[fill - 1] + 1;
            field[fieldoffset + fill] =
                eligible[selection[fill]];
        }
        return TRUE;
    }
    return FALSE;
}

static int* solver_ab_reserve_eligible_workspace(
    solver_ab_executor_t* executor,
    int worker_id,
    size_t required) {
    int* resized;
    size_t capacity;

    if (!executor ||
        worker_id < 0 ||
        worker_id >= executor->worker_count ||
        !executor->combination_eligible ||
        !executor->combination_eligible_capacity ||
        required > SIZE_MAX / sizeof(*resized)) {
        return NULL;
    }
    capacity =
        executor->combination_eligible_capacity[worker_id];
    if (capacity >= required) {
        return executor->combination_eligible[worker_id];
    }
    resized = realloc(
        executor->combination_eligible[worker_id],
        required * sizeof(*resized));
    if (!resized) {
        return NULL;
    }
    executor->combination_eligible[worker_id] = resized;
    executor->combination_eligible_capacity[worker_id] =
        required;
    return resized;
}

static int* solver_ab_get_eligible_workspace(
    solver_ab_executor_t* executor,
    int worker_id,
    size_t required) {
    return solver_ab_reserve_eligible_workspace(
        executor,
        worker_id,
        required);
}

static int solver_ab_visit_pair_range(
    solver_ab_executor_t* executor,
    int worker_id,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    solver_ab_combination_visitor_t visitor,
    void* opaque) {
    int field[DQMAX] = {0};
    int selection[DQMAX] = {0};
    int* eligible;
    int eligible_count = 0;
    int fieldoffset;
    int n_to_add;
    int i;
    unsigned long long ordinal;

    if (local_first >= local_end ||
        local_end > pair->combination_count) {
        return -1;
    }
    field[A] = pair->field_a;
    field[B] = pair->field_b;
    if (executor->phase == SOLVER_AB_PHASE_DIAGONAL) {
        fieldoffset = C;
        n_to_add = executor->dimquads - 2;
    } else {
        field[C] = executor->newpoint;
        fieldoffset = D;
        n_to_add = executor->dimquads - 3;
    }
    if (n_to_add < 0 ||
        fieldoffset + n_to_add > DQMAX) {
        return -1;
    }
    if (n_to_add == 0) {
        if (local_first != 0U ||
            local_end != 1U) {
            return -1;
        }
        return visitor(
            pair_geometry,
            pair,
            field,
            executor->dimquads,
            opaque) ? 1 : 0;
    }
    eligible = solver_ab_get_eligible_workspace(
        executor,
        worker_id,
        (size_t)executor->newpoint);
    if (!eligible) {
        return -1;
    }
    for (i = 0; i < executor->newpoint; i++) {
        double x;
        double y;

        if (solver_pair_geometry_transform(
                executor->field_geometry,
                pair_geometry,
                pair->field_a,
                pair->field_b,
                i,
                &x,
                &y)) {
            eligible[eligible_count++] = i;
        }
    }
    if (solver_ab_saturating_choose(
            eligible_count,
            n_to_add) != pair->combination_count) {
        return -1;
    }
    if (solver_ab_select_combination(
            eligible,
            eligible_count,
            field,
            fieldoffset,
            n_to_add,
            selection,
            local_first)) {
        return -1;
    }
    for (ordinal = local_first;
         ordinal < local_end;
         ordinal++) {
        if (visitor(
                pair_geometry,
                pair,
                field,
                executor->dimquads,
                opaque)) {
            return 1;
        }
        if (ordinal + 1U < local_end &&
            !solver_ab_next_combination(
                eligible,
                eligible_count,
                field,
                fieldoffset,
                n_to_add,
                selection)) {
            return -1;
        }
    }
    return 0;
}

static size_t solver_ab_find_pair_for_work(
    const solver_ab_pair_t* pairs,
    size_t pair_count,
    unsigned long long work) {
    size_t low = 0U;
    size_t high = pair_count;

    while (low < high) {
        size_t middle = low + (high - low) / 2U;

        if (pairs[middle].combination_first <= work) {
            low = middle + 1U;
        } else {
            high = middle;
        }
    }
    return low ? low - 1U : 0U;
}

static anbool solver_ab_builder_visit(
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    int* field,
    int dimquad,
    void* opaque) {
    solver_ab_builder_t* builder = opaque;

    if (solver_ab_cancelled(builder->executor)) {
        builder->packet->cancelled = TRUE;
        return TRUE;
    }
    solver_ab_try_all_codes(
        pair_geometry,
        pair->field_a,
        pair->field_b,
        field,
        dimquad,
        builder);
    return builder->fatal_error ||
        builder->packet->cancelled;
}

typedef int (*solver_ab_pair_range_visitor_t)(
    solver_ab_executor_t* executor,
    int worker_id,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque);

static int solver_ab_walk_task_ranges(
    solver_ab_executor_t* executor,
    const solver_ab_task_t* task,
    int worker_id,
    solver_ab_pair_range_visitor_t visitor,
    void* opaque) {
    unsigned long long cursor;
    size_t pair_index;

    if (!executor || !task || !visitor ||
        !executor->pairs || !executor->pair_count ||
        task->combination_first >=
            task->combination_end) {
        return -1;
    }
    cursor = task->combination_first;
    pair_index = solver_ab_find_pair_for_work(
        executor->pairs,
        executor->pair_count,
        cursor);
    for (;
         pair_index < executor->pair_count &&
             cursor < task->combination_end;
         pair_index++) {
        const solver_ab_pair_t* pair =
            &executor->pairs[pair_index];
        unsigned long long pair_end;
        unsigned long long slice_end;
        const solver_pair_geometry_t* pair_geometry;
        int visit_result;

        if (!pair->combination_count ||
            ULLONG_MAX - pair->combination_first <
                pair->combination_count) {
            return -1;
        }
        pair_end = pair->combination_first +
            pair->combination_count;
        if (pair_end <= cursor) {
            continue;
        }
        if (pair->combination_first > cursor) {
            return -1;
        }
        slice_end = MIN(
            task->combination_end,
            pair_end);
        pair_geometry = solver_field_geometry_pair(
            executor->field_geometry,
            pair->field_a,
            pair->field_b);
        if (!pair_geometry) {
            return -1;
        }
        visit_result = visitor(
            executor,
            worker_id,
            pair,
            pair_geometry,
            cursor - pair->combination_first,
            slice_end - pair->combination_first,
            opaque);
        if (visit_result) {
            return visit_result;
        }
        cursor = slice_end;
    }
    return cursor == task->combination_end ? 0 : -1;
}

static int solver_ab_builder_visit_pair_range(
    solver_ab_executor_t* executor,
    int worker_id,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque) {
    solver_ab_builder_t* builder = opaque;

    if (solver_ab_cancelled(executor)) {
        builder->packet->cancelled = TRUE;
        return 1;
    }
    builder->tol2 = pair->tol2;
    builder->rel_field_noise2 =
        pair_geometry->rel_field_noise2;
    builder->rel_field_noise_valid = TRUE;
    return solver_ab_visit_pair_range(
        executor,
        worker_id,
        pair,
        pair_geometry,
        local_first,
        local_end,
        solver_ab_builder_visit,
        builder);
}

static int solver_ab_evaluate_task(
    solver_ab_executor_t* executor,
    size_t task_index,
    int worker_id) {
    solver_ab_task_t* task = &executor->tasks[task_index];
    solver_ab_channel_t* channel =
        &executor->channels[worker_id];
    solver_ab_builder_t builder;
    int visit_result;

    memset(&builder, 0, sizeof(builder));
    builder.executor = executor;
    builder.packet = &channel->packet;
    builder.query_result =
        &executor->query_results[worker_id];
    builder.task_index = task_index;
    builder.worker_id = worker_id;

    visit_result = solver_ab_walk_task_ranges(
        executor,
        task,
        worker_id,
        solver_ab_builder_visit_pair_range,
        &builder);
    if (visit_result < 0) {
        builder.fatal_error = TRUE;
        builder.packet->evaluation_failed = TRUE;
    }
    if (builder.fatal_error ||
        builder.packet->evaluation_failed ||
        builder.packet->allocation_failed) {
        size_t hypothesis_index;
        unsigned long long search_failures = 0U;

        for (hypothesis_index = 0U;
             hypothesis_index <
                 builder.packet->hypothesis_count;
             hypothesis_index++) {
            if (builder.packet
                    ->hypotheses[hypothesis_index]
                    .search_failed) {
                search_failures++;
            }
        }
        pthread_mutex_lock(&executor->mutex);
        if (builder.packet->allocation_failed) {
            executor->pending_allocation_failures++;
        }
        executor->pending_search_failures +=
            search_failures;
        executor->pending_evaluation_failure =
            executor->pending_evaluation_failure ||
            builder.packet->evaluation_failed;
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        pthread_mutex_unlock(&executor->mutex);
    }
    if (builder.producer_completed) {
        return -1;
    }
    if (solver_ab_publish_segment(&builder, TRUE)) {
        return -1;
    }
    return 0;
}

static anbool solver_ab_task_claimable_locked(
    const solver_ab_executor_t* executor,
    int worker_id) {
    const solver_ab_channel_t* channel;
    size_t claim_limit = executor->next_reduce;

    if (worker_id <= 0 ||
        worker_id >= executor->worker_count) {
        return FALSE;
    }
    channel = &executor->channels[worker_id];
    if (SIZE_MAX - claim_limit <
        executor->reorder_window) {
        claim_limit = SIZE_MAX;
    } else {
        claim_limit += executor->reorder_window;
    }
    return executor->phase_active &&
        !executor->stopping &&
        !executor->fatal_error &&
        !channel->ready &&
        !__atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE) &&
        executor->next_task < executor->task_count &&
        executor->next_task < claim_limit;
}

static int solver_ab_claim_task_locked(
    solver_ab_executor_t* executor,
    int worker_id,
    size_t* task_index) {
    solver_ab_channel_t* channel;
    int claimed = FALSE;

    if (worker_id <= 0 ||
        worker_id >= executor->worker_count) {
        return FALSE;
    }
    channel = &executor->channels[worker_id];
    if (solver_ab_task_claimable_locked(
            executor,
            worker_id)) {
        *task_index = executor->next_task++;
        executor->tasks[*task_index].state =
            SOLVER_AB_TASK_RUNNING;
        executor->tasks[*task_index].worker_id =
            worker_id;
        channel->task_index = *task_index;
        channel->final = FALSE;
        solver_ab_packet_reset(&channel->packet);
        claimed = TRUE;
    }
    return claimed;
}

static anbool solver_ab_claim_owner_task_locked(
    solver_ab_executor_t* executor,
    size_t* task_index) {
    size_t next;

    if (!executor || !task_index ||
        !executor->phase_active ||
        executor->stopping ||
        executor->fatal_error ||
        __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE) ||
        executor->next_task >= executor->task_count ||
        executor->next_task != executor->next_reduce) {
        return FALSE;
    }
    next = executor->next_task++;
    if (executor->tasks[next].state !=
        SOLVER_AB_TASK_PENDING) {
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        return FALSE;
    }
    executor->tasks[next].state =
        SOLVER_AB_TASK_RUNNING;
    executor->tasks[next].worker_id = 0;
    *task_index = next;
    return TRUE;
}

static int solver_ab_counter_failure(
    solver_t* solver,
    const char* counter_name) {
    logerr(
        "[solver-ab] signed counter boundary reached: %s\n",
        counter_name ? counter_name : "(unknown)");
    if (solver) {
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
    }
    return -1;
}

static anbool solver_ab_counter_can_add(
    int current,
    unsigned long long delta) {
    if (current < 0 ||
        delta > (unsigned long long)INT_MAX) {
        return FALSE;
    }
    return delta <=
        (unsigned long long)(INT_MAX - current);
}

int solver_ab_checked_counter_delta(
    solver_t* solver,
    unsigned long long numtries,
    unsigned long long cxdx,
    unsigned long long meanx) {
    if (!solver) {
        return -1;
    }
    /*
     * Preflight every destination before mutating any of them. This preserves
     * an exact reducer prefix even at the representable counter boundary.
     */
    if (!solver_ab_counter_can_add(
            solver->numtries,
            numtries)) {
        return solver_ab_counter_failure(
            solver,
            "numtries");
    }
    if (!solver_ab_counter_can_add(
            solver->num_cxdx_skipped,
            cxdx)) {
        return solver_ab_counter_failure(
            solver,
            "num_cxdx_skipped");
    }
    if (!solver_ab_counter_can_add(
            solver->num_meanx_skipped,
            meanx)) {
        return solver_ab_counter_failure(
            solver,
            "num_meanx_skipped");
    }
    solver->numtries += (int)numtries;
    solver->num_cxdx_skipped += (int)cxdx;
    solver->num_meanx_skipped += (int)meanx;
    return 0;
}

static int solver_ab_reduce_packet(
    solver_ab_executor_t* executor,
    solver_ab_packet_t* packet,
    solver_ab_reduce_state_t* state) {
    solver_t* solver = executor->owner;
    size_t hypothesis_index;

    for (hypothesis_index = 0;
         hypothesis_index < packet->hypothesis_count;
         hypothesis_index++) {
        solver_ab_hypothesis_t* hypothesis =
            &packet->hypotheses[hypothesis_index];
        size_t candidate_index;
        double reduction_start = 0.0;

        if (solver_poll_worker_stop(solver)) {
            return 0;
        }
        if (hypothesis->begins_hypothesis) {
            if (state->in_hypothesis) {
                logerr("[solver-ab] hypothesis begin before prior end\n");
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                return -1;
            }
            if (hypothesis->nresults < 0) {
                logerr(
                    "[solver-ab] hypothesis has negative "
                    "candidate count\n");
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                return -1;
            }
            if (solver_ab_checked_counter_delta(
                solver,
                hypothesis->numtries_delta,
                hypothesis->cxdx_delta,
                hypothesis->meanx_delta)) {
                return -1;
            }
            state->in_hypothesis = TRUE;
            state->expected_candidates =
                hypothesis->nresults;
            state->seen_candidates = 0;
            if (solver->profile.detailed) {
                solver->profile.hypothesis_order_hash =
                    solver_order_hash_mix(
                        solver->profile.hypothesis_order_hash,
                        hypothesis->hypothesis_order_digest);
                solver->profile.kd_result_order_hash =
                    solver_order_hash_mix(
                        solver->profile.kd_result_order_hash,
                        hypothesis->kd_result_order_digest);
            }
            solver->profile.hypotheses_generated++;
            solver->profile.hypotheses_executed++;
            solver->profile.codekd_calls++;
            solver->profile.codekd_hits +=
                (unsigned long long)hypothesis->nresults;
            solver->profile.codekd_wall_seconds +=
                hypothesis->codekd_wall_seconds;
            if (solver->profile.max_batch_hypotheses < 1U) {
                solver->profile.max_batch_hypotheses = 1U;
            }

            if (hypothesis->search_failed) {
                solver->profile.search_failures++;
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                return -1;
            }
            if (hypothesis->nresults > 0) {
                solver->profile.resolve_calls++;
            }
        } else if (!state->in_hypothesis) {
            logerr("[solver-ab] hypothesis continuation without begin\n");
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            return -1;
        }
        if (hypothesis->candidate_first >
                packet->candidate_count ||
            hypothesis->candidate_count >
                packet->candidate_count -
                    hypothesis->candidate_first) {
            logerr("[solver-ab] candidate segment bounds are invalid\n");
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            return -1;
        }
        solver->profile.resolve_wall_seconds +=
            hypothesis->prepare_wall_seconds;
        if (solver->profile.detailed &&
            hypothesis->candidate_count > 0U) {
            reduction_start = monotonic_seconds();
        }

        for (candidate_index = hypothesis->candidate_first;
             candidate_index <
                hypothesis->candidate_first +
                    hypothesis->candidate_count;
             candidate_index++) {
            solver_ab_candidate_t* candidate =
                &packet->candidates[candidate_index];

            if (solver->nummatches < 0 ||
                solver->nummatches == INT_MAX ||
                state->seen_candidates < 0 ||
                state->seen_candidates == INT_MAX) {
                return solver_ab_counter_failure(
                    solver,
                    "nummatches/seen_candidates");
            }
            if (candidate->action ==
                    SOLVER_AB_CANDIDATE_RADEC_SKIP &&
                (solver->num_radec_skipped < 0 ||
                 solver->num_radec_skipped ==
                     INT_MAX)) {
                return solver_ab_counter_failure(
                    solver,
                    "num_radec_skipped");
            }
            if (candidate->action ==
                    SOLVER_AB_CANDIDATE_ABSCALE_SKIP &&
                (solver->num_abscale_skipped < 0 ||
                 solver->num_abscale_skipped ==
                     INT_MAX)) {
                return solver_ab_counter_failure(
                    solver,
                    "num_abscale_skipped");
            }
            if (candidate->action ==
                    SOLVER_AB_CANDIDATE_VERIFY &&
                (solver->numscaleok < 0 ||
                 solver->numscaleok == INT_MAX ||
                 solver->num_verified < 0 ||
                 solver->num_verified == INT_MAX)) {
                return solver_ab_counter_failure(
                    solver,
                    "numscaleok/num_verified");
            }
            solver->nummatches++;
            state->seen_candidates++;
            solver_record_candidate_order(
                solver,
                candidate->action,
                candidate->quadno,
                candidate->code_err);
            switch (candidate->action) {
            case SOLVER_AB_CANDIDATE_RADEC_SKIP:
                solver->num_radec_skipped++;
                break;

            case SOLVER_AB_CANDIDATE_ABSCALE_SKIP:
                solver->num_abscale_skipped++;
                break;

            case SOLVER_AB_CANDIDATE_BAD_QUAD:
                logverb("bad quad at %s:%i\n",
                        __FILE__,
                        __LINE__);
                break;

            case SOLVER_AB_CANDIDATE_SCALE_SKIP:
                break;

            case SOLVER_AB_CANDIDATE_VERIFY:
            {
                MatchObj mo;
                int i;

                solver->numscaleok++;
                set_matchobj_template(solver, &mo);
                memcpy(&mo.wcstan, &candidate->wcs, sizeof(tan_t));
                mo.wcs_valid = TRUE;
                mo.code_err = candidate->code_err;
                mo.scale = candidate->scale;
                mo.parity = candidate->parity;
                mo.quad_npeers = candidate->quad_npeers;
                mo.timeused = solver->timeused;
                mo.quadno = candidate->quadno;
                mo.dimquads = executor->dimquads;
                for (i = 0; i < executor->dimquads; i++) {
                    mo.star[i] = candidate->star[i];
                    mo.field[i] = candidate->field[i];
                    mo.ids[i] = 0;
                }
                memcpy(
                    mo.quadpix,
                    candidate->quadpix,
                    (size_t)2 * (size_t)executor->dimquads *
                        sizeof(double));
                memcpy(
                    mo.quadxyz,
                    candidate->quadxyz,
                    (size_t)3 * (size_t)executor->dimquads *
                        sizeof(double));
                set_center_and_radius(
                    solver,
                    &mo,
                    &mo.wcstan,
                    NULL);
                mo.quads_tried = solver->numtries;
                mo.quads_matched = solver->nummatches;
                mo.quads_scaleok = solver->numscaleok;
                if (solver_handle_hit(
                        solver,
                        &mo,
                        NULL,
                        FALSE)) {
                    solver->quit_now = TRUE;
                }
                break;
            }

            default:
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                return -1;
            }
            if (unlikely(solver->quit_now)) {
                break;
            }
        }
        if (solver->profile.detailed &&
            hypothesis->candidate_count > 0U) {
            solver->profile.resolve_wall_seconds +=
                monotonic_seconds() - reduction_start;
        }
        if (unlikely(solver->quit_now)) {
            if (state->in_hypothesis) {
                state->in_hypothesis = FALSE;
                state->expected_candidates = 0;
                state->seen_candidates = 0;
                solver->profile.hypotheses_reduced++;
            }
            return 0;
        }
        if (hypothesis->ends_hypothesis) {
            if (!state->in_hypothesis ||
                state->seen_candidates !=
                    state->expected_candidates) {
                logerr(
                    "[solver-ab] hypothesis ended after %i of %i "
                    "candidates\n",
                    state->seen_candidates,
                    state->expected_candidates);
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                return -1;
            }
            state->in_hypothesis = FALSE;
            state->expected_candidates = 0;
            state->seen_candidates = 0;
            solver->profile.hypotheses_reduced++;
        }
    }
    if (solver_ab_checked_counter_delta(
        solver,
        packet->trailing_numtries,
        packet->trailing_cxdx,
        packet->trailing_meanx)) {
        return -1;
    }
    if (packet->evaluation_failed) {
        logerr("[solver-ab] combination-range evaluation failed\n");
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        return -1;
    }
    if (packet->allocation_failed) {
        logerr("[solver-ab] bounded segment allocation failed\n");
        solver->profile.allocation_failures++;
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        return -1;
    }
    return 0;
}

typedef struct solver_ab_native_visitor {
    solver_t* solver;
    const solver_field_geometry_t* field_geometry;
    kdtree_qres_t* query_result;
    double tol2;
} solver_ab_native_visitor_t;

static void solver_ab_try_all_codes_native(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    const int* fieldstars,
    int dimquad,
    solver_ab_native_visitor_t* visitor) {
    solver_t* solver = visitor->solver;
    double code[DCMAX];
    double flipcode[DCMAX];
    int dimcode = (dimquad - NBACK) * 2;
    int i;

    solver->numtries++;
    for (i = 0; i < dimquad - NBACK; i++) {
        if (!solver_pair_geometry_transform(
                geometry,
                pair_geometry,
                pair->field_a,
                pair->field_b,
                fieldstars[NBACK + i],
                &code[2 * i],
                &code[2 * i + 1])) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            return;
        }
    }
    if (solver->parity == PARITY_NORMAL ||
        solver->parity == PARITY_BOTH) {
        try_all_codes_2(
            fieldstars,
            dimquad,
            code,
            solver,
            FALSE,
            visitor->tol2,
            &visitor->query_result);
    }
    if (solver->quit_now) {
        return;
    }
    if (solver->parity == PARITY_FLIP ||
        solver->parity == PARITY_BOTH) {
        quad_flip_parity(code, flipcode, dimcode);
        try_all_codes_2(
            fieldstars,
            dimquad,
            flipcode,
            solver,
            TRUE,
            visitor->tol2,
            &visitor->query_result);
    }
}

static anbool solver_ab_native_visit(
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    int* field,
    int dimquad,
    void* opaque) {
    solver_ab_native_visitor_t* visitor = opaque;

    if (solver_poll_worker_stop(visitor->solver)) {
        return TRUE;
    }
    solver_ab_try_all_codes_native(
        visitor->field_geometry,
        pair_geometry,
        pair,
        field,
        dimquad,
        visitor);
    return visitor->solver->quit_now;
}

static int solver_ab_native_visit_pair_range(
    solver_ab_executor_t* executor,
    int worker_id,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque) {
    solver_ab_native_visitor_t* visitor = opaque;

    if (solver_poll_worker_stop(visitor->solver)) {
        return 1;
    }
    visitor->solver->rel_field_noise2 =
        pair_geometry->rel_field_noise2;
    visitor->tol2 = pair->tol2;
    return solver_ab_visit_pair_range(
        executor,
        worker_id,
        pair,
        pair_geometry,
        local_first,
        local_end,
        solver_ab_native_visit,
        visitor);
}

static int solver_ab_run_task_inline(
    solver_ab_executor_t* executor,
    const solver_ab_task_t* task) {
    solver_t* solver = executor->owner;
    solver_ab_native_visitor_t visitor;
    int visit_result;

    memset(&visitor, 0, sizeof(visitor));
    visitor.solver = solver;
    visitor.field_geometry = executor->field_geometry;
    visitor.query_result = executor->query_results[0];
    if (visitor.query_result) {
        executor->owner_query_reuses++;
    }
    solver->profile.ab_blocks_owner++;
    visit_result = solver_ab_walk_task_ranges(
        executor,
        task,
        0,
        solver_ab_native_visit_pair_range,
        &visitor);
    if (visit_result < 0) {
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
    }
    executor->query_results[0] = visitor.query_result;
    return solver->profile.execution_failed ? -1 : 0;
}

static int solver_ab_retire_segment(
    solver_ab_executor_t* executor,
    size_t task_index,
    solver_ab_reduce_state_t* reduce_state,
    anbool* final_out) {
    solver_t* solver = executor->owner;
    solver_ab_task_t* task = &executor->tasks[task_index];
    solver_ab_channel_t* channel;
    solver_ab_packet_t swap;
    anbool final;
    int rc;

    pthread_mutex_lock(&executor->mutex);
    if (task->worker_id <= 0 ||
        task->worker_id >= executor->worker_count) {
        pthread_mutex_unlock(&executor->mutex);
        return 1;
    }
    channel = &executor->channels[task->worker_id];
    if (!channel->ready ||
        channel->task_index != task_index) {
        pthread_mutex_unlock(&executor->mutex);
        return 1;
    }
    swap = executor->drain_packet;
    executor->drain_packet = channel->packet;
    channel->packet = swap;
    solver_ab_packet_reset(&channel->packet);
    final = channel->final;
    channel->ready = FALSE;
    channel->final = FALSE;
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    pthread_mutex_unlock(&executor->mutex);

    solver->profile.ab_segments_retired++;
    solver->profile.ab_segment_payload_bytes +=
        executor->drain_packet.hypothesis_count *
            sizeof(*executor->drain_packet.hypotheses) +
        executor->drain_packet.candidate_count *
            sizeof(*executor->drain_packet.candidates);
    if (executor->drain_packet.cancelled) {
        solver_poll_worker_stop(solver);
        if (!solver->quit_now) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            rc = -1;
        } else {
            rc = 0;
        }
    } else {
        rc = solver_ab_reduce_packet(
            executor,
            &executor->drain_packet,
            reduce_state);
    }
    solver_ab_packet_reset(&executor->drain_packet);

    if (final && !rc) {
        if (reduce_state->in_hypothesis &&
            !solver->quit_now) {
            logerr("[solver-ab] task ended inside a hypothesis\n");
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            rc = -1;
        } else {
            solver->profile.ab_blocks_retired++;
            solver->profile.task_ranges_executed++;
        }
    }
    *final_out = final;
    return rc;
}

static anbool solver_ab_canonical_segment_ready_locked(
    const solver_ab_executor_t* executor,
    size_t task_index) {
    const solver_ab_task_t* task;
    const solver_ab_channel_t* channel;

    if (task_index >= executor->task_count) {
        return FALSE;
    }
    task = &executor->tasks[task_index];
    if (task->worker_id <= 0 ||
        task->worker_id >= executor->worker_count) {
        return FALSE;
    }
    channel = &executor->channels[task->worker_id];
    return channel->ready &&
        channel->task_index == task_index;
}

static void solver_ab_configure_plan_cache_limits(
    solver_ab_executor_t* executor) {
    size_t lane_budget;
    size_t task_limit;

    if (!executor || executor->worker_count <= 0) {
        return;
    }
    lane_budget =
        (size_t)SOLVER_AB_PLAN_CACHE_BUDGET_BYTES /
        (size_t)executor->worker_count;
    if (!executor->pair_cache_limit_bytes &&
        !executor->task_cache_limit_bytes) {
        if ((size_t)executor->worker_count >
            SIZE_MAX / SOLVER_AB_BLOCKS_PER_WORKER ||
            (size_t)executor->worker_count *
                SOLVER_AB_BLOCKS_PER_WORKER >
            SIZE_MAX / sizeof(solver_ab_task_t)) {
            task_limit = lane_budget;
        } else {
            task_limit =
                (size_t)executor->worker_count *
                SOLVER_AB_BLOCKS_PER_WORKER *
                sizeof(solver_ab_task_t);
            task_limit = MIN(task_limit, lane_budget);
        }
        executor->task_cache_limit_bytes =
            task_limit;
        executor->pair_cache_limit_bytes =
            lane_budget - task_limit;
    }
    if (!executor->packet_cache_limit_bytes) {
        executor->packet_cache_limit_bytes =
            (size_t)SOLVER_AB_PACKET_CACHE_BUDGET_BYTES /
            (size_t)executor->worker_count;
    }
    if (!executor->eligible_cache_limit_bytes) {
        executor->eligible_cache_limit_bytes =
            (size_t)SOLVER_AB_ELIGIBLE_CACHE_BUDGET_BYTES /
            (size_t)executor->worker_count;
    }
}

static int solver_ab_reserve_pair_buffer(
    solver_ab_executor_t* executor,
    size_t required) {
    solver_ab_pair_t* resized;
    size_t capacity;
    size_t retained_maximum;

    if (!executor ||
        required > SIZE_MAX / sizeof(*resized)) {
        return -1;
    }
    solver_ab_configure_plan_cache_limits(executor);
    if (required <= executor->pair_capacity) {
        return 0;
    }
    capacity = executor->pair_capacity ?
        executor->pair_capacity : 64U;
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            capacity = required;
        } else {
            capacity *= 2U;
        }
    }

    retained_maximum =
        executor->pair_cache_limit_bytes /
        sizeof(*resized);
    if (!executor->pairs_transient &&
        required <= retained_maximum) {
        if (capacity > retained_maximum) {
            capacity = retained_maximum;
        }
        resized = realloc(
            executor->pair_cache,
            capacity * sizeof(*resized));
        if (!resized) {
            return -1;
        }
        executor->pair_cache = resized;
        executor->pair_cache_capacity = capacity;
        executor->pairs = resized;
        executor->pair_capacity = capacity;
        executor->pair_cache_grows++;
        return 0;
    }

    if (capacity > SIZE_MAX / sizeof(*resized)) {
        capacity = required;
    }
    if (executor->pairs_transient) {
        resized = realloc(
            executor->pairs,
            capacity * sizeof(*resized));
    } else {
        resized = malloc(capacity * sizeof(*resized));
        if (resized && required > 1U && executor->pairs) {
            memcpy(
                resized,
                executor->pairs,
                (required - 1U) * sizeof(*resized));
        }
    }
    if (!resized) {
        return -1;
    }
    executor->pairs = resized;
    executor->pair_capacity = capacity;
    executor->pairs_transient = TRUE;
    executor->plan_transient_allocations++;
    executor->plan_transient_peak_bytes = MAX(
        executor->plan_transient_peak_bytes,
        capacity * sizeof(*resized));
    return 0;
}

static int solver_ab_reserve_task_buffer(
    solver_ab_executor_t* executor,
    size_t required) {
    solver_ab_task_t* resized;
    size_t bytes;

    if (!executor ||
        required > SIZE_MAX / sizeof(*resized)) {
        return -1;
    }
    solver_ab_configure_plan_cache_limits(executor);
    if (!executor->tasks &&
        required <= executor->task_cache_capacity) {
        executor->tasks = executor->task_cache;
        executor->task_capacity =
            executor->task_cache_capacity;
        executor->tasks_transient = FALSE;
        executor->task_cache_reuses++;
        return 0;
    }
    bytes = required * sizeof(*resized);
    if (required <= executor->task_capacity) {
        return 0;
    }
    if (bytes <= executor->task_cache_limit_bytes) {
        resized = realloc(
            executor->task_cache,
            bytes);
        if (!resized) {
            return -1;
        }
        executor->task_cache = resized;
        executor->task_cache_capacity = required;
        executor->tasks = resized;
        executor->task_capacity = required;
        executor->tasks_transient = FALSE;
        executor->task_cache_grows++;
        return 0;
    }
    resized = malloc(bytes);
    if (!resized) {
        return -1;
    }
    executor->tasks = resized;
    executor->task_capacity = required;
    executor->tasks_transient = TRUE;
    executor->plan_transient_allocations++;
    executor->plan_transient_peak_bytes = MAX(
        executor->plan_transient_peak_bytes,
        bytes);
    return 0;
}

static void solver_ab_begin_pair_buffer(
    solver_ab_executor_t* executor) {
    executor->pairs = executor->pair_cache;
    executor->pair_capacity =
        executor->pair_cache_capacity;
    executor->pairs_transient = FALSE;
    if (executor->pair_cache_capacity) {
        executor->pair_cache_reuses++;
    }
}

static int solver_ab_collect_pairs(
    solver_ab_executor_t* executor,
    solver_ab_phase_kind_t phase,
    int newpoint,
    int dimquads,
    const solver_field_geometry_t* geometry,
    double min_ab2,
    double max_ab2,
    solver_ab_pair_t** pairs_out,
    size_t* pair_count_out,
    unsigned long long* total_combinations_out) {
    size_t count = 0U;
    unsigned long long total_combinations = 0;
    int field_a;
    int field_b;

    if (!executor || executor->pairs ||
        executor->pair_count != 0U) {
        return -1;
    }
    solver_ab_begin_pair_buffer(executor);
    if (newpoint <= 0) {
        *pairs_out = NULL;
        *pair_count_out = 0U;
        *total_combinations_out = 0U;
        return 0;
    }

    if (phase == SOLVER_AB_PHASE_DIAGONAL) {
        field_b = newpoint;
        for (field_a = 0;
             field_a < newpoint;
             field_a++) {
            const solver_pair_geometry_t* pair_geometry =
                solver_field_geometry_pair(
                    geometry,
                    field_a,
                    field_b);
            unsigned long long combination_count;

            if (!pair_geometry ||
                !pair_geometry->scale_ok ||
                pair_geometry->scale < min_ab2 ||
                pair_geometry->scale > max_ab2) {
                continue;
            }
            combination_count =
                solver_ab_saturating_choose(
                solver_pair_geometry_eligible_before(
                    geometry,
                    pair_geometry,
                    field_a,
                    field_b,
                    newpoint),
                dimquads - 2);
            if (!combination_count) {
                continue;
            }
            if (solver_ab_reserve_pair_buffer(
                    executor, count + 1U)) {
                return -1;
            }
            executor->pairs[count].field_a = field_a;
            executor->pairs[count].field_b = field_b;
            executor->pairs[count].combination_first =
                total_combinations;
            executor->pairs[count].combination_count =
                combination_count;
            executor->pairs[count].tol2 = get_tolerance_for_noise(
                executor->snapshot.codetol,
                pair_geometry->rel_field_noise2,
                executor->snapshot.rel_index_noise2);
            count++;
            total_combinations =
                solver_ab_saturating_add(
                    total_combinations,
                    combination_count);
        }
    } else {
        for (field_a = 0;
             field_a < newpoint;
             field_a++) {
            for (field_b = field_a + 1;
                 field_b < newpoint;
                 field_b++) {
                const solver_pair_geometry_t* pair_geometry =
                    solver_field_geometry_pair(
                        geometry,
                        field_a,
                        field_b);
                unsigned long long combination_count;
                double newpoint_x;
                double newpoint_y;

                if (!pair_geometry ||
                    !pair_geometry->scale_ok ||
                    pair_geometry->scale < min_ab2 ||
                    pair_geometry->scale > max_ab2 ||
                    !solver_pair_geometry_transform(
                        geometry,
                        pair_geometry,
                        field_a,
                        field_b,
                        newpoint,
                        &newpoint_x,
                        &newpoint_y)) {
                    continue;
                }
                combination_count = dimquads > 3 ?
                    solver_ab_saturating_choose(
                        solver_pair_geometry_eligible_before(
                            geometry,
                            pair_geometry,
                            field_a,
                            field_b,
                            newpoint),
                        dimquads - 3) :
                    1U;
                if (!combination_count) {
                    continue;
                }
                if (solver_ab_reserve_pair_buffer(
                        executor, count + 1U)) {
                    return -1;
                }
                executor->pairs[count].field_a = field_a;
                executor->pairs[count].field_b = field_b;
                executor->pairs[count].combination_first =
                    total_combinations;
                executor->pairs[count].combination_count =
                    combination_count;
                executor->pairs[count].tol2 = get_tolerance_for_noise(
                    executor->snapshot.codetol,
                    pair_geometry->rel_field_noise2,
                    executor->snapshot.rel_index_noise2);
                count++;
                total_combinations =
                    solver_ab_saturating_add(
                        total_combinations,
                        combination_count);
            }
        }
    }

    if (!count) {
        *pairs_out = NULL;
        *pair_count_out = 0U;
        *total_combinations_out = 0U;
        return 0;
    }
    *pairs_out = executor->pairs;
    *pair_count_out = count;
    *total_combinations_out = total_combinations;
    return 0;
}

#define SOLVER_PAYLOAD_PAGE_PLAN_INLINE_RANGES 128U
#define SOLVER_PAYLOAD_PAGE_PLAN_MAX_RANGES 256U
#define SOLVER_PAYLOAD_PAGE_PLAN_MAX_BYTES \
    (4U * 1024U * 1024U)

typedef struct solver_payload_page_plan {
    fitsbin_t* fitsbin;
    fitsbin_prefetch_range_t* ranges;
    size_t range_count;
    size_t range_capacity;
    size_t logical_bytes;
    anbool failed;
    fitsbin_prefetch_range_t
        inline_ranges[SOLVER_PAYLOAD_PAGE_PLAN_INLINE_RANGES];
} solver_payload_page_plan_t;

static void solver_payload_page_plan_init(
    solver_payload_page_plan_t* plan) {
    if (!plan) {
        return;
    }
    memset(plan, 0, sizeof(*plan));
    plan->ranges = plan->inline_ranges;
    plan->range_capacity =
        SOLVER_PAYLOAD_PAGE_PLAN_INLINE_RANGES;
}

static void solver_payload_page_plan_destroy(
    solver_payload_page_plan_t* plan) {
    if (!plan) {
        return;
    }
    if (plan->ranges && plan->ranges != plan->inline_ranges) {
        free(plan->ranges);
    }
    memset(plan, 0, sizeof(*plan));
}

static void solver_payload_page_plan_fail(
    solver_payload_page_plan_t* plan) {
    if (!plan || plan->failed) {
        return;
    }
    plan->failed = TRUE;
}

static int solver_payload_page_plan_enabled(
    void* userdata,
    void* mapping) {
    solver_payload_page_plan_t* plan = userdata;

    if (!plan || !mapping || plan->failed ||
        fitsbin_payload_io_demand_busy()) {
        return FALSE;
    }
    if (plan->fitsbin && plan->fitsbin != mapping) {
        return FALSE;
    }
    return TRUE;
}

static int solver_payload_page_plan_reserve(
    solver_payload_page_plan_t* plan,
    size_t needed) {
    fitsbin_prefetch_range_t* ranges;
    size_t capacity;

    if (!plan || needed <= plan->range_capacity) {
        return 0;
    }
    if (needed > SOLVER_PAYLOAD_PAGE_PLAN_MAX_RANGES) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    capacity = plan->range_capacity;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2U) {
            solver_payload_page_plan_fail(plan);
            return -1;
        }
        capacity *= 2U;
    }
    if (capacity > SIZE_MAX / sizeof(*ranges)) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    ranges = malloc(capacity * sizeof(*ranges));
    if (!ranges) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    memcpy(ranges,
           plan->ranges,
           plan->range_count * sizeof(*ranges));
    if (plan->ranges != plan->inline_ranges) {
        free(plan->ranges);
    }
    plan->ranges = ranges;
    plan->range_capacity = capacity;
    return 0;
}

static int solver_payload_page_plan_emit(
    void* userdata,
    const kdtree_prefetch_hint_t* hint) {
    solver_payload_page_plan_t* plan = userdata;

    if (!plan || !hint || !hint->mapping ||
        !hint->address || !hint->length) {
        return 0;
    }
    if (hint->kind != KDTREE_PREFETCH_ARRAY_DATA &&
        hint->kind != KDTREE_PREFETCH_ARRAY_PERM) {
        return 0;
    }
    if (plan->failed ||
        (plan->fitsbin && plan->fitsbin != hint->mapping)) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    if (!plan->fitsbin) {
        plan->fitsbin = hint->mapping;
    }
    if (plan->range_count == SIZE_MAX ||
        solver_payload_page_plan_reserve(
            plan, plan->range_count + 1U)) {
        return -1;
    }
    plan->ranges[plan->range_count].data = hint->address;
    plan->ranges[plan->range_count].size = hint->length;
    plan->range_count++;
    if (hint->length > SIZE_MAX - plan->logical_bytes) {
        plan->logical_bytes = SIZE_MAX;
    } else {
        plan->logical_bytes += hint->length;
    }
    return 0;
}

static size_t solver_payload_page_plan_budget(
    const solver_payload_page_plan_t* plan) {
    size_t budget;
    size_t page_size;
    long detected_page_size;

    if (!plan || !plan->range_count || !plan->logical_bytes) {
        return 0U;
    }
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0) {
        return 0U;
    }
    page_size = (size_t)detected_page_size;
    budget = plan->logical_bytes;
    if (plan->range_count >
        (SIZE_MAX - budget) / (2U * page_size)) {
        budget = SIZE_MAX;
    } else {
        budget += plan->range_count * 2U * page_size;
    }

    /*
     * This is a latency bound, not a host-memory heuristic. fitsbin divides
     * its process-wide population budget among configured worker credits.
     */
    return MIN(
        budget,
        (size_t)SOLVER_PAYLOAD_PAGE_PLAN_MAX_BYTES);
}

static int solver_payload_page_plan_add_query(
    solver_payload_page_plan_t* plan,
    const kdtree_t* tree,
    const void* query,
    double maxd2,
    int options) {
    kdtree_prefetch_sink_t sink;
    int status;

    if (!plan || !tree || !query ||
        !tree->io || !tree->io_is_fitsbin) {
        return 0;
    }
    if (fitsbin_payload_io_demand_busy()) {
        if (!plan->fitsbin) {
            plan->fitsbin = tree->io;
        }
        solver_payload_page_plan_fail(plan);
        return 0;
    }
    memset(&sink, 0, sizeof(sink));
    sink.userdata = plan;
    sink.enabled = solver_payload_page_plan_enabled;
    sink.emit = solver_payload_page_plan_emit;
    status = kdtree_rangesearch_prefetch_prepare(
        tree, query, maxd2, options, &sink);
    if (status || plan->failed) {
        if (!plan->fitsbin) {
            plan->fitsbin = tree->io;
        }
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    return 0;
}

static int solver_payload_page_plan_flush(
    solver_payload_page_plan_t* plan) {
    size_t budget;
    int status;

    if (!plan) {
        return -1;
    }
    if (plan->failed || fitsbin_payload_io_demand_busy()) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    if (!plan->fitsbin || !plan->range_count) {
        return 0;
    }
    budget = solver_payload_page_plan_budget(plan);
    if (!budget) {
        solver_payload_page_plan_fail(plan);
        return -1;
    }
    status = fitsbin_advise_mapped_ranges(
        plan->fitsbin,
        plan->ranges,
        plan->range_count,
        budget);
    if (status <= 0) {
        solver_payload_page_plan_fail(plan);
        return status < 0 ? -1 : 0;
    }
    return status;
}

static void solver_payload_advise_verification_search(
    solver_t* solver,
    const MatchObj* mo) {
    solver_payload_page_plan_t page_plan;
    int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_RETURN_POINTS;

    if (!solver || !mo || !solver->index ||
        !solver->index->starkd ||
        !solver->index->starkd->tree) {
        return;
    }
    solver_payload_page_plan_init(&page_plan);
    (void)solver_payload_page_plan_add_query(
        &page_plan,
        solver->index->starkd->tree,
        mo->center,
        square(mo->radius),
        options);
    (void)solver_payload_page_plan_flush(&page_plan);
    solver_payload_page_plan_destroy(&page_plan);
}

typedef struct solver_verification_score_slot {
    verify_prepared_hit_t* prepared;
} solver_verification_score_slot_t;

typedef struct solver_verification_candidate_runtime {
    MatchObj match;
    double match_distance_in_pixels2;
} solver_verification_candidate_runtime_t;

typedef struct solver_verification_task_input {
    const solver_verification_score_slot_t* slots;
    size_t slot_count;
} solver_verification_task_input_t;

static index_shard_helper_task_status_t
solver_verification_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_verification_task_input_t* input = input_bytes;
    verify_prepared_score_t* scores = output_bytes;
    size_t i;

    if (!input || input_size != sizeof(*input) ||
        !input->slots || !input->slot_count || !scores ||
        input->slot_count >
            SIZE_MAX / sizeof(*scores) ||
        output_size !=
            input->slot_count * sizeof(*scores)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    for (i = 0U; i < input->slot_count; i++) {
        const solver_verification_score_slot_t* slot =
            &input->slots[i];

        if (index_shard_worker_stop_requested()) {
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        memset(&scores[i], 0, sizeof(scores[i]));
        if (!slot->prepared ||
            verify_score_prepared_hit(
                slot->prepared,
                &scores[i])) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

static const index_shard_helper_ops_t
solver_verification_helper_ops = {
    "verify-context",
    solver_verification_helper_execute
};

#define SOLVER_VERIFICATION_WAVE_MAX_BYTES \
    (32U * 1024U * 1024U)

static size_t solver_verification_wave_memory_budget(void) {
    size_t budget = SOLVER_VERIFICATION_WAVE_MAX_BYTES;

#if defined(_SC_AVPHYS_PAGES)
    long available_pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);

    if (available_pages > 0 && page_size > 0 &&
        (unsigned long)available_pages <=
            SIZE_MAX / (unsigned long)page_size) {
        size_t pressure_budget =
            ((size_t)available_pages * (size_t)page_size) / 64U;

        budget = MIN(budget, pressure_budget);
    }
#endif
    return budget;
}

static void solver_verification_wave_cleanup(
    solver_verification_score_slot_t* slots,
    verify_prepared_score_t* scores,
    size_t slot_count,
    solver_verification_candidate_runtime_t* runtime,
    solver_ab_packet_t* packet) {
    size_t i;

    if (slots) {
        for (i = 0U; i < slot_count; i++) {
            if (scores) {
                verify_destroy_prepared_score(&scores[i]);
            }
            verify_destroy_prepared_hit(slots[i].prepared);
            slots[i].prepared = NULL;
        }
    }
    free(slots);
    free(scores);
    free(runtime);
    solver_ab_packet_free(packet);
}

/*
 * Prepare index-backed candidate inputs on the owner, score immutable
 * verification contexts in coarse helper tasks, then retire every action in
 * the original candidate order. No helper receives an index, mapping, solver,
 * callback, or reducer pointer.
 *
 * Return zero before publication/state mutation to request the native path.
 * Return one after successful handling, stop propagation, or a hard failure.
 */
static int solver_ab_try_verification_wave(
    kdtree_qres_t* result,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    const solver_candidate_payload_workspace_t* payload_workspace) {
    solver_ab_packet_t packet;
    solver_ab_snapshot_t snapshot;
    solver_verification_candidate_runtime_t* runtime = NULL;
    solver_verification_score_slot_t* slots = NULL;
    verify_prepared_score_t* scores = NULL;
    solver_payload_page_plan_t page_plan;
    index_shard_helper_task_t tasks[INDEX_SHARD_HELPER_MAX_TASKS];
    solver_verification_task_input_t
        inputs[INDEX_SHARD_HELPER_MAX_TASKS];
    index_shard_helper_run_stats_t run_stats;
    index_shard_helper_run_status_t run_status;
    size_t available;
    size_t verify_count = 0U;
    size_t peak_bytes = 0U;
    size_t peak_budget;
    size_t task_count;
    size_t task_index;
    size_t slot_cursor = 0U;
    double verify_wall_start;
    double verify_wall_seconds;
    int candidate_index;
    int handled = 0;

    if (!result || !field_xy || !fieldstars || !solver ||
        result->nres < 2 || verify_datalog_enabled() ||
        !index_shard_worker_context_active()) {
        return 0;
    }
    available = index_shard_helper_prepare_reserve();
    if (!available) {
        return 0;
    }

    memset(&packet, 0, sizeof(packet));
    memset(&snapshot, 0, sizeof(snapshot));
    memset(tasks, 0, sizeof(tasks));
    memset(inputs, 0, sizeof(inputs));
    memset(&run_stats, 0, sizeof(run_stats));
    solver_payload_page_plan_init(&page_plan);

    if (solver_ab_packet_reserve_candidates(
            &packet, (size_t)result->nres) !=
        SOLVER_AB_RESERVE_OK) {
        goto cleanup;
    }
    runtime = calloc(
        (size_t)result->nres,
        sizeof(*runtime));
    slots = calloc(
        (size_t)result->nres,
        sizeof(*slots));
    scores = calloc(
        (size_t)result->nres,
        sizeof(*scores));
    if (!runtime || !slots || !scores) {
        goto cleanup;
    }

    snapshot.index = solver->index;
    snapshot.fieldxy = solver->fieldxy;
    snapshot.use_radec = solver->use_radec;
    snapshot.cx_less_than_dx = solver->index->cx_less_than_dx;
    snapshot.meanx_less_than_half =
        solver->index->meanx_less_than_half;
    snapshot.parity = solver->parity;
    memcpy(snapshot.centerxyz,
           solver->centerxyz,
           sizeof(snapshot.centerxyz));
    snapshot.r2 = solver->r2;
    snapshot.abscale_low = solver->abscale_low;
    snapshot.abscale_high = solver->abscale_high;
    snapshot.funits_lower = solver->funits_lower;
    snapshot.funits_upper = solver->funits_upper;
    snapshot.cxdx_margin = solver->cxdx_margin;
    snapshot.codetol = solver->codetol;
    snapshot.rel_index_noise2 = solver->rel_index_noise2;

    for (candidate_index = 0;
         candidate_index < result->nres;
         candidate_index++) {
        solver_ab_candidate_t* candidate =
            &packet.candidates[candidate_index];
        const unsigned int* prepared_stars = NULL;

        if (payload_workspace) {
            if (!payload_workspace->valid[candidate_index]) {
                goto cleanup;
            }
            prepared_stars = payload_workspace->stars +
                (size_t)candidate_index * (size_t)DQMAX;
        }
        if (solver_ab_candidate_prepare(
                candidate,
                result,
                candidate_index,
                field_xy,
                fieldstars,
                dimquads,
                &snapshot,
                current_parity,
                prepared_stars)) {
            goto cleanup;
        }
        packet.candidate_count++;
        if (candidate->action ==
            SOLVER_AB_CANDIDATE_VERIFY) {
            verify_count++;
        }
    }
    if (verify_count < 2U) {
        goto cleanup;
    }

    verify_wall_start = monotonic_seconds();
    for (candidate_index = 0;
         candidate_index < result->nres;
         candidate_index++) {
        solver_ab_candidate_t* candidate =
            &packet.candidates[candidate_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[candidate_index];
        double logaccept;
        int i;

        if (candidate->action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        set_matchobj_template(solver, &candidate_runtime->match);
        memcpy(&candidate_runtime->match.wcstan,
               &candidate->wcs,
               sizeof(tan_t));
        candidate_runtime->match.wcs_valid = TRUE;
        candidate_runtime->match.code_err = candidate->code_err;
        candidate_runtime->match.scale = candidate->scale;
        candidate_runtime->match.parity = candidate->parity;
        candidate_runtime->match.quad_npeers =
            candidate->quad_npeers;
        candidate_runtime->match.timeused = solver->timeused;
        candidate_runtime->match.quadno = candidate->quadno;
        candidate_runtime->match.dimquads = dimquads;
        for (i = 0; i < dimquads; i++) {
            candidate_runtime->match.star[i] = candidate->star[i];
            candidate_runtime->match.field[i] = candidate->field[i];
            candidate_runtime->match.ids[i] = 0;
        }
        memcpy(candidate_runtime->match.quadpix,
               candidate->quadpix,
               (size_t)2 * (size_t)dimquads * sizeof(double));
        memcpy(candidate_runtime->match.quadxyz,
               candidate->quadxyz,
               (size_t)3 * (size_t)dimquads * sizeof(double));
        set_center_and_radius(
            solver,
            &candidate_runtime->match,
            &candidate_runtime->match.wcstan,
            NULL);
        candidate_runtime->match_distance_in_pixels2 =
            solver_prepare_hit_for_verify(
                solver,
                &candidate_runtime->match,
                &logaccept);
        slot_cursor++;
        if (solver_payload_page_plan_add_query(
                &page_plan,
                solver->index->starkd->tree,
                candidate_runtime->match.center,
                square(candidate_runtime->match.radius),
                KD_OPTIONS_SMALL_RADIUS |
                    KD_OPTIONS_RETURN_POINTS)) {
            break;
        }
    }
    (void)solver_payload_page_plan_flush(&page_plan);
    if (slot_cursor != verify_count) {
        goto cleanup;
    }

    peak_budget = solver_verification_wave_memory_budget();
    slot_cursor = 0U;
    for (candidate_index = 0;
         candidate_index < result->nres;
         candidate_index++) {
        solver_ab_candidate_t* candidate =
            &packet.candidates[candidate_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[candidate_index];
        solver_verification_score_slot_t* slot;
        size_t context_peak;
        double logaccept;

        if (candidate->action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        slot = &slots[slot_cursor++];
        (void)solver_prepare_hit_for_verify(
            solver,
            &candidate_runtime->match,
            &logaccept);
        if (verify_prepare_hit(
                solver->index->starkd,
                solver->index->cutnside,
                &candidate_runtime->match,
                NULL,
                solver->vf,
                candidate_runtime->match_distance_in_pixels2,
                solver->distractor_ratio,
                solver->field_maxx,
                solver->field_maxy,
                solver->logratio_bail_threshold,
                logaccept,
                solver->logratio_stoplooking,
                solver->distance_from_quad_bonus,
                FALSE,
                &slot->prepared)) {
            goto cleanup;
        }
        context_peak =
            verify_prepared_hit_peak_bytes(slot->prepared);
        if (context_peak == SIZE_MAX ||
            peak_bytes > peak_budget ||
            context_peak > peak_budget - peak_bytes) {
            goto cleanup;
        }
        peak_bytes += context_peak;
    }
    if (slot_cursor != verify_count) {
        goto cleanup;
    }

    task_count = MIN(
        verify_count,
        MIN(available + 1U,
            (size_t)INDEX_SHARD_HELPER_MAX_TASKS));
    if (task_count < 2U) {
        goto cleanup;
    }
    slot_cursor = 0U;
    for (task_index = 0U;
         task_index < task_count;
         task_index++) {
        size_t remaining_slots = verify_count - slot_cursor;
        size_t remaining_tasks = task_count - task_index;
        size_t count =
            (remaining_slots + remaining_tasks - 1U) /
            remaining_tasks;
        size_t i;
        unsigned long long work_units = 0U;

        inputs[task_index].slots = &slots[slot_cursor];
        inputs[task_index].slot_count = count;
        for (i = 0U; i < count; i++) {
            unsigned long long work =
                verify_prepared_hit_work_units(
                    slots[slot_cursor + i].prepared);

            if (!work) {
                work = 1U;
            }
            if (ULLONG_MAX - work_units < work) {
                work_units = ULLONG_MAX;
            } else {
                work_units += work;
            }
        }
        tasks[task_index].input = &inputs[task_index];
        tasks[task_index].input_bytes = sizeof(inputs[task_index]);
        tasks[task_index].output = &scores[slot_cursor];
        tasks[task_index].output_bytes =
            count * sizeof(*scores);
        tasks[task_index].work_units = work_units;
        slot_cursor += count;
    }
    if (slot_cursor != verify_count) {
        goto cleanup;
    }

    run_status = index_shard_helper_run(
        &solver_verification_helper_ops,
        tasks,
        task_count,
        &run_stats);
    if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE ||
        run_status == INDEX_SHARD_HELPER_TASK_FAILED) {
        goto cleanup;
    }
    if (run_status == INDEX_SHARD_HELPER_STOPPED) {
        (void)solver_poll_worker_stop(solver);
        solver->quit_now = TRUE;
        handled = 1;
        goto cleanup;
    }
    if (run_status != INDEX_SHARD_HELPER_OK) {
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        handled = 1;
        goto cleanup;
    }
    verify_wall_seconds =
        monotonic_seconds() - verify_wall_start;
    if (solver->profile.detailed) {
        solver->profile.verify_wall_seconds +=
            verify_wall_seconds;
    }
    if (run_stats.foreign_tasks) {
        solver->profile.parallel_batches++;
        solver->profile.parallel_batches_observed++;
        solver->profile.ab_helper_tasks =
            solver_ab_saturating_add(
                solver->profile.ab_helper_tasks,
                run_stats.foreign_tasks);
        solver->profile.max_parallel_ranges = MAX(
            solver->profile.max_parallel_ranges,
            run_stats.max_concurrent_tasks);
    }

    slot_cursor = 0U;
    for (candidate_index = 0;
         candidate_index < result->nres;
         candidate_index++) {
        solver_ab_candidate_t* candidate =
            &packet.candidates[candidate_index];

        if (solver_poll_worker_stop(solver)) {
            handled = 1;
            goto cleanup;
        }
        if (solver->nummatches < 0 ||
            solver->nummatches == INT_MAX) {
            (void)solver_ab_counter_failure(
                solver, "nummatches");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_RADEC_SKIP &&
            (solver->num_radec_skipped < 0 ||
             solver->num_radec_skipped == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "num_radec_skipped");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_ABSCALE_SKIP &&
            (solver->num_abscale_skipped < 0 ||
             solver->num_abscale_skipped == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "num_abscale_skipped");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_VERIFY &&
            (solver->numscaleok < 0 ||
             solver->numscaleok == INT_MAX ||
             solver->num_verified < 0 ||
             solver->num_verified == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "numscaleok/num_verified");
            handled = 1;
            goto cleanup;
        }
        solver->nummatches++;
        solver_record_candidate_order(
            solver,
            candidate->action,
            candidate->quadno,
            (float)candidate->code_err);
        switch (candidate->action) {
        case SOLVER_AB_CANDIDATE_RADEC_SKIP:
            solver->num_radec_skipped++;
            break;

        case SOLVER_AB_CANDIDATE_ABSCALE_SKIP:
            solver->num_abscale_skipped++;
            break;

        case SOLVER_AB_CANDIDATE_BAD_QUAD:
            logverb("bad quad at %s:%i\n", __FILE__, __LINE__);
            break;

        case SOLVER_AB_CANDIDATE_SCALE_SKIP:
            break;

        case SOLVER_AB_CANDIDATE_VERIFY:
        {
            solver_verification_candidate_runtime_t* candidate_runtime =
                &runtime[candidate_index];
            solver_verification_score_slot_t* slot =
                &slots[slot_cursor++];

            solver->numscaleok++;
            candidate_runtime->match.quads_tried = quads_tried;
            candidate_runtime->match.quads_matched =
                solver->nummatches;
            candidate_runtime->match.quads_scaleok =
                solver->numscaleok;
            if (verify_finish_prepared_hit(
                    slot->prepared,
                    &scores[slot_cursor - 1U],
                    &candidate_runtime->match)) {
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                handled = 1;
                goto cleanup;
            }
            solver->profile.verify_calls++;
            if (solver_handle_hit_after_verify(
                    solver,
                    &candidate_runtime->match,
                    NULL,
                    FALSE,
                    candidate_runtime->match_distance_in_pixels2)) {
                solver->quit_now = TRUE;
            }
            verify_destroy_prepared_hit(slot->prepared);
            slot->prepared = NULL;
            if (unlikely(solver->quit_now)) {
                handled = 1;
                goto cleanup;
            }
            break;
        }

        default:
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
    }
    handled = 1;

cleanup:
    index_shard_helper_prepare_cancel();
    solver_payload_page_plan_destroy(&page_plan);
    solver_verification_wave_cleanup(
        slots,
        scores,
        verify_count,
        runtime,
        &packet);
    return handled;
}

typedef struct solver_ab_descriptor_task_input {
    const solver_field_geometry_t* field_geometry;
    const solver_ab_pair_t* pairs;
    size_t pair_count;
    unsigned long long combination_first;
    unsigned long long combination_end;
    solver_ab_phase_kind_t phase;
    int newpoint;
    int dimquads;
    int parity;
    anbool cx_less_than_dx;
    anbool meanx_less_than_half;
    double cxdx_margin;
} solver_ab_descriptor_task_input_t;

typedef struct solver_ab_descriptor_workspace {
    solver_ab_executor_t planner;
    solver_ab_descriptor_output_t* outputs;
} solver_ab_descriptor_workspace_t;

static pthread_key_t solver_ab_descriptor_workspace_key;
static pthread_once_t solver_ab_descriptor_workspace_once =
    PTHREAD_ONCE_INIT;
static int solver_ab_descriptor_workspace_status = EAGAIN;

static void solver_ab_descriptor_release_pairs(
    solver_ab_descriptor_workspace_t* workspace) {
    solver_ab_executor_t* planner;

    if (!workspace) {
        return;
    }
    planner = &workspace->planner;
    if (planner->pairs_transient) {
        free(planner->pairs);
    }
    planner->pairs = NULL;
    planner->pair_count = 0U;
    planner->pair_capacity = 0U;
    planner->pairs_transient = FALSE;
}

static void solver_ab_descriptor_workspace_destroy(void* opaque) {
    solver_ab_descriptor_workspace_t* workspace = opaque;

    if (!workspace) {
        return;
    }
    solver_ab_descriptor_release_pairs(workspace);
    free(workspace->planner.pair_cache);
    free(workspace->outputs);
    free(workspace);
}

static void solver_ab_descriptor_workspace_make_key(void) {
    solver_ab_descriptor_workspace_status =
        pthread_key_create(
            &solver_ab_descriptor_workspace_key,
            solver_ab_descriptor_workspace_destroy);
}

static solver_ab_descriptor_workspace_t*
solver_ab_descriptor_workspace_get(void) {
    solver_ab_descriptor_workspace_t* workspace;

    if (pthread_once(
            &solver_ab_descriptor_workspace_once,
            solver_ab_descriptor_workspace_make_key) ||
        solver_ab_descriptor_workspace_status) {
        return NULL;
    }
    workspace = pthread_getspecific(
        solver_ab_descriptor_workspace_key);
    if (workspace) {
        return workspace;
    }
    workspace = calloc(1, sizeof(*workspace));
    if (!workspace) {
        return NULL;
    }
    workspace->outputs = calloc(
        SOLVER_AB_DESCRIPTOR_MAX_TASKS,
        sizeof(*workspace->outputs));
    if (!workspace->outputs) {
        free(workspace);
        return NULL;
    }
    workspace->planner.worker_count =
        (int)SOLVER_AB_DESCRIPTOR_MAX_TASKS;
    workspace->planner.pair_cache_limit_bytes =
        SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES;
    if (pthread_setspecific(
            solver_ab_descriptor_workspace_key,
            workspace)) {
        solver_ab_descriptor_workspace_destroy(workspace);
        return NULL;
    }
    return workspace;
}

static index_shard_helper_task_status_t
solver_ab_descriptor_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_ab_descriptor_task_input_t* input = input_bytes;
    solver_ab_descriptor_output_t* output = output_bytes;
    solver_ab_executor_t executor;
    solver_ab_builder_t builder;
    solver_ab_packet_t packet;
    solver_ab_task_t task;
    int eligible[SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS];
    int* eligible_pointer = eligible;
    size_t eligible_capacity =
        SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS;
    int visit_result;

    if (!input ||
        input_size != sizeof(*input) ||
        !output ||
        output_size != sizeof(*output) ||
        !input->field_geometry ||
        !input->pairs || !input->pair_count ||
        input->combination_first >= input->combination_end ||
        input->newpoint < 0 ||
        input->newpoint >=
            SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS ||
        input->dimquads < NBACK ||
        input->dimquads > DQMAX ||
        (input->phase != SOLVER_AB_PHASE_DIAGONAL &&
         input->phase != SOLVER_AB_PHASE_OFF_DIAGONAL) ||
        (input->parity != PARITY_NORMAL &&
         input->parity != PARITY_FLIP &&
         input->parity != PARITY_BOTH)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    memset(&executor, 0, sizeof(executor));
    memset(&builder, 0, sizeof(builder));
    memset(&packet, 0, sizeof(packet));
    memset(&task, 0, sizeof(task));
    output->descriptor_count = 0U;
    output->trailing_numtries = 0U;
    output->trailing_cxdx = 0U;
    output->trailing_meanx = 0U;
    output->has_final_rel_field_noise2 = FALSE;

    executor.worker_count = 1;
    executor.field_geometry = input->field_geometry;
    executor.newpoint = input->newpoint;
    executor.dimquads = input->dimquads;
    executor.phase = input->phase;
    executor.pairs = (solver_ab_pair_t*)input->pairs;
    executor.pair_count = input->pair_count;
    executor.combination_eligible = &eligible_pointer;
    executor.combination_eligible_capacity =
        &eligible_capacity;
    executor.snapshot.parity = input->parity;
    executor.snapshot.cx_less_than_dx =
        input->cx_less_than_dx;
    executor.snapshot.meanx_less_than_half =
        input->meanx_less_than_half;
    executor.snapshot.cxdx_margin = input->cxdx_margin;

    builder.executor = &executor;
    builder.packet = &packet;
    builder.descriptor_output = output;
    builder.worker_id = 0;
    task.combination_first = input->combination_first;
    task.combination_end = input->combination_end;
    visit_result = solver_ab_walk_task_ranges(
        &executor,
        &task,
        0,
        solver_ab_builder_visit_pair_range,
        &builder);
    if (packet.cancelled ||
        index_shard_worker_stop_requested()) {
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }
    if (visit_result < 0 ||
        builder.fatal_error ||
        packet.evaluation_failed ||
        packet.allocation_failed) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    output->trailing_numtries = builder.pending_numtries;
    output->trailing_cxdx = builder.pending_cxdx;
    output->trailing_meanx = builder.pending_meanx;
    output->final_rel_field_noise2 =
        builder.rel_field_noise2;
    output->has_final_rel_field_noise2 =
        builder.rel_field_noise_valid;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static const index_shard_helper_ops_t
solver_ab_descriptor_helper_ops = {
    "ab-descriptor",
    solver_ab_descriptor_helper_execute
};

static size_t solver_ab_descriptor_expansion(
    int dimquads,
    int parity) {
    size_t expansion = 2U;
    int internal_stars = dimquads - NBACK;
    int i;

    if (dimquads < NBACK || dimquads > DQMAX ||
        (parity != PARITY_NORMAL &&
         parity != PARITY_FLIP &&
         parity != PARITY_BOTH)) {
        return 0U;
    }
    if (parity == PARITY_BOTH) {
        expansion *= 2U;
    }
    for (i = 2; i <= internal_stars; i++) {
        expansion *= (size_t)i;
    }
    return expansion;
}

static int solver_ab_descriptor_reduce_output(
    solver_t* solver,
    const solver_ab_descriptor_output_t* output,
    int dimquads,
    kdtree_qres_t** query_result,
    unsigned long long* reduced) {
    size_t descriptor_index;

    for (descriptor_index = 0U;
         descriptor_index < output->descriptor_count;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &output->descriptors[descriptor_index];

        if (solver_poll_worker_stop(solver)) {
            return 1;
        }
        solver->rel_field_noise2 =
            descriptor->rel_field_noise2;
        if (solver_ab_checked_counter_delta(
                solver,
                descriptor->numtries_delta,
                descriptor->cxdx_delta,
                descriptor->meanx_delta)) {
            return -1;
        }
        solver_execute_hypothesis_owner(
            descriptor->stars,
            descriptor->code,
            dimquads,
            solver,
            descriptor->current_parity,
            descriptor->tol2,
            query_result);
        (*reduced)++;
        if (solver->profile.execution_failed) {
            return -1;
        }
        if (solver->quit_now) {
            return 1;
        }
    }
    if (output->has_final_rel_field_noise2) {
        solver->rel_field_noise2 =
            output->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            solver,
            output->trailing_numtries,
            output->trailing_cxdx,
            output->trailing_meanx)) {
        return -1;
    }
    return 0;
}

/*
 * Generate only index-free hypothesis descriptors on helper threads. The
 * owning worker retains every CodeKD, QuadFile, StarKD, verification and
 * reducer operation and consumes completed descriptors in canonical order.
 */
static int solver_ab_descriptor_execute_phase(
    solver_t* solver,
    solver_ab_phase_kind_t phase,
    int newpoint,
    const solver_field_geometry_t* field_geometry,
    int dimquads,
    double min_ab2,
    double max_ab2,
    kdtree_qres_t** query_result,
    solver_ab_phase_mode_t* mode_out) {
    solver_ab_descriptor_workspace_t* workspace;
    solver_ab_executor_t* planner;
    solver_ab_pair_t* pairs = NULL;
    size_t pair_count = 0U;
    unsigned long long total_combinations = 0U;
    unsigned long long combination_cursor = 0U;
    unsigned long long reduced = 0U;
    unsigned long long parallel_reduced = 0U;
    size_t available;
    size_t participants;
    size_t expansion;
    size_t max_task_combinations;
    anbool assisted = FALSE;
    anbool handled = FALSE;
    int result = 0;

    if (mode_out) {
        *mode_out = SOLVER_AB_MODE_NATIVE;
    }
    if (!solver || !field_geometry || !query_result ||
        solver->ab_executor || solver->maxquads != 0 ||
        solver->maxmatches != 0 ||
        pl_size(solver->indexes) != 1U ||
        !index_shard_worker_context_active()) {
        return 0;
    }
    available = index_shard_helper_available_workers();
    if (!available) {
        return 0;
    }
    participants = MIN(
        SOLVER_AB_DESCRIPTOR_MAX_TASKS,
        available + 1U);
    if (participants < 2U) {
        return 0;
    }
    expansion = solver_ab_descriptor_expansion(
        dimquads,
        solver->parity);
    if (!expansion ||
        expansion > SOLVER_AB_DESCRIPTOR_CAPACITY) {
        return 0;
    }
    max_task_combinations =
        SOLVER_AB_DESCRIPTOR_CAPACITY / expansion;
    if (max_task_combinations <
        SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS) {
        return 0;
    }

    workspace = solver_ab_descriptor_workspace_get();
    if (!workspace) {
        return 0;
    }
    planner = &workspace->planner;
    planner->snapshot.codetol = solver->codetol;
    planner->snapshot.rel_index_noise2 =
        solver->rel_index_noise2;
    if (solver_ab_collect_pairs(
            planner,
            phase,
            newpoint,
            dimquads,
            field_geometry,
            min_ab2,
            max_ab2,
            &pairs,
            &pair_count,
            &total_combinations)) {
        solver_ab_descriptor_release_pairs(workspace);
        return 0;
    }
    planner->pair_count = pair_count;
    if (!pair_count || !total_combinations ||
        total_combinations == ULLONG_MAX ||
        total_combinations <
            2U * SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS) {
        solver_ab_descriptor_release_pairs(workspace);
        return 0;
    }

    while (combination_cursor < total_combinations) {
        solver_ab_descriptor_task_input_t
            inputs[SOLVER_AB_DESCRIPTOR_MAX_TASKS];
        index_shard_helper_task_t
            tasks[SOLVER_AB_DESCRIPTOR_MAX_TASKS];
        unsigned long long remaining =
            total_combinations - combination_cursor;
        unsigned long long wave_capacity =
            (unsigned long long)max_task_combinations *
            (unsigned long long)participants;
        unsigned long long wave_combinations =
            MIN(remaining, wave_capacity);
        unsigned long long base;
        unsigned long long remainder;
        unsigned long long task_cursor = combination_cursor;
        size_t task_count =
            (size_t)((wave_combinations +
                      max_task_combinations - 1U) /
                     max_task_combinations);
        size_t task_index;
        index_shard_helper_run_status_t run_status =
            INDEX_SHARD_HELPER_UNAVAILABLE;
        index_shard_helper_run_stats_t run_stats;
        anbool run_inline = FALSE;
        anbool wave_assisted = FALSE;
        unsigned long long reduced_before = reduced;

        if (task_count < 2U &&
            wave_combinations >=
                2U * SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS) {
            task_count = 2U;
        }
        if (task_count < 2U && !handled) {
            result = 0;
            goto cleanup;
        }
        if (!task_count) {
            result = -1;
            goto fail;
        }
        base = wave_combinations /
            (unsigned long long)task_count;
        remainder = wave_combinations %
            (unsigned long long)task_count;
        memset(inputs, 0, sizeof(inputs));
        memset(tasks, 0, sizeof(tasks));
        memset(&run_stats, 0, sizeof(run_stats));
        for (task_index = 0U;
             task_index < task_count;
             task_index++) {
            unsigned long long task_combinations =
                base + (task_index < remainder ? 1U : 0U);
            unsigned long long work_units;

            if (!task_combinations ||
                task_combinations > max_task_combinations) {
                result = -1;
                goto fail;
            }
            inputs[task_index].field_geometry = field_geometry;
            inputs[task_index].pairs = pairs;
            inputs[task_index].pair_count = pair_count;
            inputs[task_index].combination_first = task_cursor;
            task_cursor += task_combinations;
            inputs[task_index].combination_end = task_cursor;
            inputs[task_index].phase = phase;
            inputs[task_index].newpoint = newpoint;
            inputs[task_index].dimquads = dimquads;
            inputs[task_index].parity = solver->parity;
            inputs[task_index].cx_less_than_dx =
                solver->index->cx_less_than_dx;
            inputs[task_index].meanx_less_than_half =
                solver->index->meanx_less_than_half;
            inputs[task_index].cxdx_margin =
                solver->cxdx_margin;
            work_units = task_combinations *
                (unsigned long long)expansion;
            tasks[task_index].input = &inputs[task_index];
            tasks[task_index].input_bytes =
                sizeof(inputs[task_index]);
            tasks[task_index].output =
                &workspace->outputs[task_index];
            tasks[task_index].output_bytes =
                sizeof(workspace->outputs[task_index]);
            tasks[task_index].work_units = work_units;
        }
        if (task_cursor !=
            combination_cursor + wave_combinations) {
            result = -1;
            goto fail;
        }

        if (task_count > 1U) {
            run_status = index_shard_helper_run(
                &solver_ab_descriptor_helper_ops,
                tasks,
                task_count,
                &run_stats);
            if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE) {
                if (!handled) {
                    result = 0;
                    goto cleanup;
                }
                run_inline = TRUE;
            } else if (run_status == INDEX_SHARD_HELPER_STOPPED) {
                (void)solver_poll_worker_stop(solver);
                solver->quit_now = TRUE;
                result = 1;
                goto cleanup;
            } else if (run_status != INDEX_SHARD_HELPER_OK) {
                result = -1;
                goto fail;
            } else {
                if (run_stats.owner_tasks +
                        run_stats.foreign_tasks != task_count ||
                    !run_stats.owner_tasks ||
                    run_stats.foreign_work_units %
                        (unsigned long long)expansion) {
                    result = -1;
                    goto fail;
                }
                solver->profile.task_ranges_inline +=
                    run_stats.owner_tasks;
                if (run_stats.foreign_tasks) {
                    wave_assisted = TRUE;
                    assisted = TRUE;
                    solver->profile.parallel_batches++;
                    solver->profile.parallel_batches_observed++;
                    solver->profile.ab_helper_tasks =
                        solver_ab_saturating_add(
                            solver->profile.ab_helper_tasks,
                            run_stats.foreign_tasks);
                    solver->profile.ab_helper_combinations =
                        solver_ab_saturating_add(
                            solver->profile.ab_helper_combinations,
                            run_stats.foreign_work_units /
                                (unsigned long long)expansion);
                    solver->profile.max_parallel_ranges = MAX(
                        solver->profile.max_parallel_ranges,
                        run_stats.max_concurrent_tasks);
                }
            }
        } else {
            run_inline = TRUE;
        }
        if (run_inline) {
            solver->profile.task_ranges_inline += task_count;
            for (task_index = 0U;
                 task_index < task_count;
                 task_index++) {
                index_shard_helper_task_status_t task_status =
                    solver_ab_descriptor_helper_execute(
                        &inputs[task_index],
                        sizeof(inputs[task_index]),
                        &workspace->outputs[task_index],
                        sizeof(workspace->outputs[task_index]));

                if (task_status ==
                    INDEX_SHARD_HELPER_TASK_STOPPED) {
                    (void)solver_poll_worker_stop(solver);
                    result = 1;
                    goto cleanup;
                }
                if (task_status !=
                    INDEX_SHARD_HELPER_TASK_OK) {
                    result = -1;
                    goto fail;
                }
            }
        }

        solver->profile.hypothesis_batches++;
        solver->profile.task_ranges_planned += task_count;
        solver->profile.task_ranges_submitted += task_count;
        solver->profile.task_ranges_executed += task_count;
        solver->profile.max_task_ranges = MAX(
            solver->profile.max_task_ranges,
            task_count);
        for (task_index = 0U;
             task_index < task_count;
             task_index++) {
            int reduce_status =
                solver_ab_descriptor_reduce_output(
                    solver,
                    &workspace->outputs[task_index],
                    dimquads,
                    query_result,
                    &reduced);

            if (reduce_status < 0) {
                result = -1;
                goto fail;
            }
            if (reduce_status > 0) {
                solver->profile.hypothesis_batches_stopped++;
                result = 1;
                goto cleanup;
            }
        }
        if (wave_assisted) {
            parallel_reduced =
                solver_ab_saturating_add(
                    parallel_reduced,
                    reduced - reduced_before);
        }
        solver->profile.hypothesis_batches_completed++;
        handled = TRUE;
        combination_cursor += wave_combinations;
    }

    result = 1;
    goto cleanup;

fail:
    solver->profile.hypothesis_batches_failed++;
    solver->profile.execution_failed = TRUE;
    solver->quit_now = TRUE;

cleanup:
    solver->profile.parallel_hypotheses =
        solver_ab_saturating_add(
            solver->profile.parallel_hypotheses,
            parallel_reduced);
    if (mode_out && assisted) {
        *mode_out = SOLVER_AB_MODE_ASSISTED;
    }
    solver_ab_descriptor_release_pairs(workspace);
    return result;
}

static int solver_ab_plan_tasks(
    solver_ab_executor_t* executor,
    size_t pair_count,
    unsigned long long total_combinations,
    int participants,
    solver_ab_task_t** tasks_out,
    size_t* task_count_out) {
    size_t desired;
    solver_ab_task_t* tasks;
    size_t task_count;
    unsigned long long base;
    unsigned long long remainder;
    unsigned long long cursor = 0U;
    size_t task_index;

    if (!executor || executor->tasks ||
        executor->task_count != 0U ||
        participants < 2 ||
        !pair_count ||
        total_combinations == ULLONG_MAX) {
        return 0;
    }
    if ((size_t)participants >
        SIZE_MAX / SOLVER_AB_BLOCKS_PER_WORKER) {
        return 0;
    }
    desired =
        (size_t)participants *
            SOLVER_AB_BLOCKS_PER_WORKER;
#if defined(SOLVER_AB_TEST_MAX_BLOCKS)
    if (SOLVER_AB_TEST_MAX_BLOCKS > 0) {
        desired = MIN(
            desired,
            (size_t)SOLVER_AB_TEST_MAX_BLOCKS);
    }
#endif
    if (total_combinations <
        (unsigned long long)desired) {
        desired = (size_t)total_combinations;
    }
    task_count = desired;
    if (task_count < 2U ||
        total_combinations < 2U) {
        return 0;
    }
    if (solver_ab_reserve_task_buffer(
            executor, task_count)) {
        return -1;
    }
    tasks = executor->tasks;
    memset(tasks, 0, task_count * sizeof(*tasks));

    base = total_combinations /
        (unsigned long long)task_count;
    remainder = total_combinations %
        (unsigned long long)task_count;
    for (task_index = 0;
         task_index < task_count;
         task_index++) {
        unsigned long long task_weight =
            base +
            (task_index < (size_t)remainder ? 1U : 0U);

        tasks[task_index].combination_first = cursor;
        cursor += task_weight;
        tasks[task_index].combination_end = cursor;
        tasks[task_index].state = SOLVER_AB_TASK_PENDING;
        tasks[task_index].worker_id = -1;
        tasks[task_index].producer_done = FALSE;
    }
    if (cursor != total_combinations) {
        return -1;
    }

    *tasks_out = tasks;
    *task_count_out = task_count;
    return 1;
}

static anbool solver_ab_plan_is_canonical(
    const solver_ab_pair_t* pairs,
    size_t pair_count,
    const solver_ab_task_t* tasks,
    size_t task_count,
    unsigned long long total_combinations) {
    unsigned long long cursor = 0U;
    size_t i;

    if (!pairs || !pair_count ||
        !tasks || !task_count ||
        !total_combinations ||
        total_combinations == ULLONG_MAX) {
        return FALSE;
    }
    for (i = 0; i < pair_count; i++) {
        if (!pairs[i].combination_count ||
            pairs[i].combination_first != cursor ||
            ULLONG_MAX - cursor <
                pairs[i].combination_count) {
            return FALSE;
        }
        cursor += pairs[i].combination_count;
    }
    if (cursor != total_combinations) {
        return FALSE;
    }
    cursor = 0U;
    for (i = 0; i < task_count; i++) {
        if (tasks[i].combination_first != cursor ||
            tasks[i].combination_end <= cursor) {
            return FALSE;
        }
        cursor = tasks[i].combination_end;
    }
    return cursor == total_combinations;
}

static size_t solver_ab_count_intra_pair_splits(
    const solver_ab_pair_t* pairs,
    size_t pair_count,
    const solver_ab_task_t* tasks,
    size_t task_count) {
    size_t pair_index = 0U;
    size_t task_index;
    size_t splits = 0U;

    if (!pairs || !pair_count || !tasks || task_count < 2U) {
        return 0U;
    }
    for (task_index = 0U;
         task_index + 1U < task_count;
         task_index++) {
        unsigned long long boundary =
            tasks[task_index].combination_end;
        unsigned long long pair_end;

        while (pair_index < pair_count) {
            pair_end =
                pairs[pair_index].combination_first +
                pairs[pair_index].combination_count;
            if (boundary <= pair_end) {
                break;
            }
            pair_index++;
        }
        if (pair_index >= pair_count) {
            break;
        }
        pair_end =
            pairs[pair_index].combination_first +
            pairs[pair_index].combination_count;
        if (boundary < pair_end) {
            splits++;
        }
    }
    return splits;
}

/*
 * Keep persistent executor memory proportional to the fixed pool, not P^2.
 *
 * At every phase boundary all helpers have left and every channel is empty.
 * Working packets may grow to the per-packet production bound while active,
 * but only a pool-scaled share is retained. Eligible-combination workspaces
 * are allocated lazily per participating lane and trimmed by the same rule.
 */
static void solver_ab_trim_epoch_buffers(
    solver_ab_executor_t* executor) {
    size_t packet_bytes;
    size_t eligible_bytes = 0U;
    int worker_id;

    if (!executor) {
        return;
    }
    solver_ab_configure_plan_cache_limits(executor);
    packet_bytes =
        executor->drain_packet.allocated_bytes;
    for (worker_id = 0;
         worker_id < executor->worker_count;
         worker_id++) {
        if (executor->channels) {
            size_t channel_bytes =
                executor->channels[worker_id]
                    .packet.allocated_bytes;

            if (packet_bytes >
                SIZE_MAX - channel_bytes) {
                packet_bytes = SIZE_MAX;
            } else {
                packet_bytes += channel_bytes;
            }
        }
        if (executor->combination_eligible &&
            executor->combination_eligible_capacity &&
            executor->combination_eligible_capacity[
                worker_id] <=
            SIZE_MAX / sizeof(int)) {
            size_t workspace_bytes =
                executor->combination_eligible_capacity[
                    worker_id] *
                sizeof(int);

            if (eligible_bytes >
                SIZE_MAX - workspace_bytes) {
                eligible_bytes = SIZE_MAX;
            } else {
                eligible_bytes += workspace_bytes;
            }
        } else {
            eligible_bytes = SIZE_MAX;
        }
    }
    for (worker_id = executor->worker_count - 1;
         executor->channels &&
             worker_id >= 0 &&
             packet_bytes >
                 executor->packet_cache_limit_bytes;
         worker_id--) {
        solver_ab_packet_t* packet =
            &executor->channels[worker_id].packet;
        size_t released;

        if (executor->channels[worker_id].ready) {
            continue;
        }
        released = packet->allocated_bytes;
        solver_ab_packet_free(packet);
        packet_bytes -= MIN(
            packet_bytes,
            released);
        if (released) {
            executor->packet_cache_trims++;
        }
    }
    if (packet_bytes >
            executor->packet_cache_limit_bytes) {
        size_t released =
            executor->drain_packet.allocated_bytes;

        solver_ab_packet_free(
            &executor->drain_packet);
        packet_bytes -= MIN(
            packet_bytes,
            released);
        if (released) {
            executor->packet_cache_trims++;
        }
    }
    for (worker_id = executor->worker_count - 1;
         executor->combination_eligible &&
             executor->combination_eligible_capacity &&
             worker_id >= 0 &&
             eligible_bytes >
                 executor->eligible_cache_limit_bytes;
         worker_id--) {
        size_t released =
            executor->combination_eligible_capacity[
                worker_id] *
            sizeof(int);

        free(executor->combination_eligible[
            worker_id]);
        executor->combination_eligible[worker_id] =
            NULL;
        executor->combination_eligible_capacity[
            worker_id] = 0U;
        eligible_bytes -= MIN(
            eligible_bytes,
            released);
        if (released) {
            executor->eligible_cache_trims++;
        }
    }
}

static void solver_ab_free_phase(
    solver_ab_executor_t* executor) {
    if (!executor) {
        return;
    }
    if (executor->pairs_transient) {
        free(executor->pairs);
    }
    if (executor->tasks_transient) {
        free(executor->tasks);
    }
    executor->pairs = NULL;
    executor->pair_capacity = 0U;
    executor->pairs_transient = FALSE;
    executor->tasks = NULL;
    executor->task_capacity = 0U;
    executor->tasks_transient = FALSE;
    executor->task_count = 0U;
    executor->pair_count = 0U;
    solver_ab_trim_epoch_buffers(executor);
}

static void solver_ab_cancel_phase(
    solver_ab_executor_t* executor) {
    size_t i;
    int worker_id;

    pthread_mutex_lock(&executor->mutex);
    __atomic_store_n(
        &executor->cancel_requested,
        TRUE,
        __ATOMIC_RELEASE);
    for (i = executor->next_task;
         i < executor->task_count;
         i++) {
        solver_ab_task_t* task = &executor->tasks[i];

        if (task->state != SOLVER_AB_TASK_PENDING) {
            continue;
        }
        solver_ab_complete_task_locked(executor, i);
    }
    executor->next_task = executor->task_count;
    for (worker_id = 1;
         worker_id < executor->worker_count;
         worker_id++) {
        solver_ab_channel_t* channel =
            &executor->channels[worker_id];

        if (!channel->ready) {
            continue;
        }
        channel->ready = FALSE;
        channel->final = FALSE;
        solver_ab_packet_reset(&channel->packet);
    }
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    pthread_mutex_unlock(&executor->mutex);
}

static void solver_ab_timespec_after_wait(
    struct timespec* deadline) {
    clock_gettime(CLOCK_REALTIME, deadline);
    deadline->tv_nsec += SOLVER_AB_WAIT_NANOSECONDS;
    if (deadline->tv_nsec >= 1000000000L) {
        deadline->tv_sec++;
        deadline->tv_nsec -= 1000000000L;
    }
}

static anbool solver_ab_poll_phase_stop(
    solver_t* solver,
    time_t* next_timer_callback_time) {
    time_t delay;
    time_t now;

    if (solver->quit_now ||
        solver_poll_worker_stop(solver)) {
        return TRUE;
    }
    if (!solver->timer_callback ||
        !next_timer_callback_time) {
        return FALSE;
    }
    now = time(NULL);
    if (now <= *next_timer_callback_time) {
        return FALSE;
    }
    update_timeused(solver);
    delay = solver->timer_callback(solver->userdata);
    if (!delay || solver->quit_now) {
        solver->quit_now = TRUE;
        return TRUE;
    }
    *next_timer_callback_time = now + delay;
    return FALSE;
}

static void solver_ab_apply_latched_failures(
    solver_ab_executor_t* executor,
    solver_t* solver) {
    unsigned long long allocation_failures;
    unsigned long long search_failures;
    anbool evaluation_failure;

    if (!executor || !solver) {
        return;
    }
    pthread_mutex_lock(&executor->mutex);
    allocation_failures =
        executor->pending_allocation_failures;
    search_failures =
        executor->pending_search_failures;
    evaluation_failure =
        executor->pending_evaluation_failure;
    executor->pending_allocation_failures = 0U;
    executor->pending_search_failures = 0U;
    executor->pending_evaluation_failure = FALSE;
    pthread_mutex_unlock(&executor->mutex);

    solver->profile.allocation_failures +=
        allocation_failures;
    solver->profile.search_failures += search_failures;
    if (allocation_failures ||
        search_failures ||
        evaluation_failure) {
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
    }
}

static int solver_ab_executor_available_lenders(
    const solver_ab_executor_t* executor) {
    int lenders;

    if (!executor || !executor->available_lenders) {
        return 0;
    }
    /*
     * This is a performance hint only. The callback is deliberately invoked
     * without executor->mutex because the fixed pool scans lanes under its
     * queue mutex before trying the executor lock.
     */
    lenders = executor->available_lenders(
        executor->work_notify_opaque);
    return lenders > 0 ? lenders : 0;
}

static unsigned long long solver_ab_publish_threshold(
    int participants) {
    unsigned long long threshold;

    if (participants <= 0) {
        return 0U;
    }
    if ((unsigned long long)participants >
        ULLONG_MAX /
            (unsigned long long)SOLVER_AB_BLOCKS_PER_WORKER) {
        threshold = ULLONG_MAX;
    } else {
        threshold =
            (unsigned long long)participants *
            (unsigned long long)SOLVER_AB_BLOCKS_PER_WORKER;
    }
    return MIN(
        threshold,
        (unsigned long long)SOLVER_AB_INLINE_PUBLISH_MAX);
}

static anbool solver_ab_phase_is_publishable(
    unsigned long long total_combinations,
    int participants,
    int available_lenders) {
    if (!total_combinations || participants < 2) {
        return FALSE;
    }
    return available_lenders > 0 ||
        total_combinations >
            solver_ab_publish_threshold(participants);
}

static int solver_ab_execute_phase(
    solver_t* solver,
    solver_ab_phase_kind_t phase,
    int newpoint,
    const solver_field_geometry_t* field_geometry,
    int dimquads,
    double min_ab2,
    double max_ab2,
    time_t* next_timer_callback_time,
    solver_ab_phase_mode_t* mode_out) {
    solver_ab_executor_t* executor = solver->ab_executor;
    solver_ab_pair_t* pairs = NULL;
    solver_ab_task_t* tasks = NULL;
    size_t pair_count = 0U;
    size_t task_count = 0U;
    unsigned long long total_combinations = 0U;
    size_t task_index;
    solver_ab_reduce_state_t reduce_state;
    int plan_result;
    int participants;
    int rc = 0;
    anbool helper_observed = FALSE;
    size_t participants_observed = 1U;
    double phase_wall_start;
    double planning_wall_start;
    unsigned long long hypotheses_before;

    if (mode_out) {
        *mode_out = SOLVER_AB_MODE_NATIVE;
    }
    if (!executor ||
        !field_geometry ||
        solver->maxquads != 0 ||
        solver->maxmatches != 0 ||
        pl_size(solver->indexes) != 1U) {
        return 0;
    }
    pthread_mutex_lock(&executor->mutex);
    if (!executor->bound ||
        executor->stopping ||
        executor->join_closed ||
        executor->fatal_error ||
        __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE)) {
        anbool failed = executor->fatal_error;

        pthread_mutex_unlock(&executor->mutex);
        if (failed) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            return -1;
        }
        return 0;
    }
    pthread_mutex_unlock(&executor->mutex);
    /*
     * Plan against the fixed pool width, not the number of idle workers at
     * this instant. The owner can retire the same canonical task ranges
     * alone, while workers that finish outer index work later can join the
     * still-live phase. This prevents a phase from becoming permanently
     * serial merely because all workers were busy at its first instruction.
     */
    participants = executor->worker_count;
    if (participants < 2) {
        return 0;
    }
    planning_wall_start = monotonic_seconds();
    if (solver_ab_capture_snapshot(executor, solver)) {
        solver->profile.allocation_failures++;
        solver->profile.ab_planning_wall_seconds +=
            monotonic_seconds() - planning_wall_start;
        solver_ab_free_phase(executor);
        return 0;
    }
    if (solver_ab_collect_pairs(
            executor,
            phase,
            newpoint,
            dimquads,
            field_geometry,
            min_ab2,
            max_ab2,
            &pairs,
            &pair_count,
            &total_combinations)) {
        solver->profile.allocation_failures++;
        solver->profile.ab_planning_wall_seconds +=
            monotonic_seconds() - planning_wall_start;
        /*
         * A failed growth may leave either the retained cache or a transient
         * overflow buffer selected as this phase's active view.  Retire that
         * view before falling back so the next frontier remains eligible for
         * assistance.
         */
        solver_ab_free_phase(executor);
        return 0;
    }
    if (!pair_count) {
        solver->profile.ab_planning_wall_seconds +=
            monotonic_seconds() - planning_wall_start;
        solver_ab_free_phase(executor);
        if (mode_out) {
            *mode_out = SOLVER_AB_MODE_EMPTY;
        }
        return 1;
    }
    /*
     * With no worker currently scanning/waiting for a loan, a plan containing
     * at most one combination per prospective task cannot amortize task
     * publication and retirement. Run that exact work through the native
     * canonical loop. Larger measured work stays publishable so a worker that
     * finishes an outer index later can still migrate into the live phase.
     */
    if (solver_ab_phase_is_publishable(
            total_combinations,
            participants,
            solver_ab_executor_available_lenders(executor))) {
        plan_result = solver_ab_plan_tasks(
            executor,
            pair_count,
            total_combinations,
            participants,
            &tasks,
            &task_count);
    } else {
        plan_result = 0;
    }
    if (plan_result <= 0) {
        /*
         * Compact geometry has no mutable pquad fallback. Keep tiny phases
         * and task-plan allocation pressure on the same exact flattened
         * sequence, using one owner task when helper publication cannot
         * amortize itself.
         */
        if (executor->tasks ||
            solver_ab_reserve_task_buffer(executor, 1U)) {
            solver->profile.allocation_failures++;
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            solver->profile.ab_planning_wall_seconds +=
                monotonic_seconds() - planning_wall_start;
            solver_ab_free_phase(executor);
            return -1;
        }
        tasks = executor->tasks;
        task_count = 1U;
        memset(tasks, 0, sizeof(*tasks));
        tasks[0].combination_first = 0U;
        tasks[0].combination_end = total_combinations;
        tasks[0].state = SOLVER_AB_TASK_PENDING;
        tasks[0].worker_id = -1;
        tasks[0].producer_done = FALSE;
    }
    if (!solver_ab_plan_is_canonical(
            pairs,
            pair_count,
            tasks,
            task_count,
            total_combinations)) {
        logerr("[solver-ab] non-canonical combination plan\n");
        solver_ab_free_phase(executor);
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        solver->profile.ab_planning_wall_seconds +=
            monotonic_seconds() - planning_wall_start;
        return -1;
    }

    pthread_mutex_lock(&executor->mutex);
    if (!executor->bound ||
        executor->stopping ||
        executor->phase_active ||
        executor->fatal_error ||
        __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE)) {
        pthread_mutex_unlock(&executor->mutex);
        solver_ab_free_phase(executor);
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        solver->profile.ab_planning_wall_seconds +=
            monotonic_seconds() - planning_wall_start;
        return -1;
    }
    executor->field_geometry = field_geometry;
    executor->newpoint = newpoint;
    executor->dimquads = dimquads;
    executor->phase = phase;
    executor->min_ab2 = min_ab2;
    executor->max_ab2 = max_ab2;
    executor->detailed = solver->profile.detailed;
    executor->pair_count = pair_count;
    executor->task_count = task_count;
    executor->next_task = 1U;
    executor->next_reduce = 0U;
    executor->reorder_window = MIN(
        task_count,
        (size_t)participants *
            SOLVER_AB_REORDER_WINDOWS_PER_WORKER);
    executor->remaining_tasks = task_count;
    tasks[0].state = SOLVER_AB_TASK_RUNNING;
    tasks[0].worker_id = 0;
    __atomic_store_n(
        &executor->cancel_requested,
        FALSE,
        __ATOMIC_RELEASE);
    executor->phase_active = TRUE;
    executor->generation++;
    if (mode_out) {
        *mode_out = SOLVER_AB_MODE_FLATTENED_OWNER;
    }
    pthread_cond_broadcast(&executor->work_cv);
    pthread_mutex_unlock(&executor->mutex);
    solver->profile.ab_planning_wall_seconds +=
        monotonic_seconds() - planning_wall_start;
    solver_ab_executor_notify_work(executor);

    solver->profile.ab_blocks_planned += task_count;
    solver->profile.ab_pairs_planned += pair_count;
    solver->profile.ab_combinations_planned +=
        total_combinations;
    solver->profile.ab_intra_pair_splits +=
        solver_ab_count_intra_pair_splits(
            pairs,
            pair_count,
            tasks,
            task_count);
    for (task_index = 0U;
         task_index < pair_count;
         task_index++) {
        solver->profile.ab_max_pair_combinations =
            MAX(solver->profile.ab_max_pair_combinations,
                pairs[task_index].combination_count);
    }
    solver->profile.hypothesis_batches++;
    solver->profile.task_ranges_planned += task_count;
    solver->profile.task_ranges_submitted += task_count;
    solver->profile.parallel_batches++;
    solver->profile.max_task_ranges =
        MAX(solver->profile.max_task_ranges, task_count);
    hypotheses_before =
        solver->profile.hypotheses_executed;
    phase_wall_start = monotonic_seconds();
    memset(&reduce_state, 0, sizeof(reduce_state));

    rc = solver_ab_run_task_inline(executor, &tasks[0]);
    pthread_mutex_lock(&executor->mutex);
    solver_ab_complete_task_locked(executor, 0U);
    executor->next_reduce = 1U;
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    pthread_mutex_unlock(&executor->mutex);
    solver->profile.ab_blocks_retired++;
    solver->profile.task_ranges_executed++;
    solver->profile.task_ranges_inline++;
    solver_ab_executor_notify_work(executor);

    if (!rc &&
        !solver->quit_now &&
        solver_ab_poll_phase_stop(
            solver,
            next_timer_callback_time)) {
        solver_ab_cancel_phase(executor);
    } else if (rc || solver->quit_now) {
        solver_ab_cancel_phase(executor);
    }

    while (!rc &&
           !solver->quit_now &&
           executor->next_reduce < task_count) {
        size_t owner_task_index = SIZE_MAX;
        anbool owner_claimed;
        anbool final = FALSE;
        int retire_result;

        pthread_mutex_lock(&executor->mutex);
        owner_claimed =
            solver_ab_claim_owner_task_locked(
                executor,
                &owner_task_index);
        pthread_mutex_unlock(&executor->mutex);
        if (owner_claimed) {
            rc = solver_ab_run_task_inline(
                executor,
                &tasks[owner_task_index]);
            pthread_mutex_lock(&executor->mutex);
            if (owner_task_index !=
                executor->next_reduce) {
                executor->fatal_error = TRUE;
                __atomic_store_n(
                    &executor->cancel_requested,
                    TRUE,
                    __ATOMIC_RELEASE);
                rc = -1;
            }
            solver_ab_complete_task_locked(
                executor,
                owner_task_index);
            executor->next_reduce++;
            pthread_cond_broadcast(
                &executor->work_cv);
            pthread_cond_broadcast(
                &executor->done_cv);
            pthread_mutex_unlock(&executor->mutex);
            solver->profile.ab_blocks_retired++;
            solver->profile.task_ranges_executed++;
            solver->profile.task_ranges_inline++;
            solver_ab_executor_notify_work(executor);
            if (rc || solver->quit_now ||
                solver_ab_poll_phase_stop(
                    solver,
                    next_timer_callback_time)) {
                solver_ab_cancel_phase(executor);
                break;
            }
            continue;
        }

        retire_result = solver_ab_retire_segment(
            executor,
            executor->next_reduce,
            &reduce_state,
            &final);

        if (retire_result < 0) {
            rc = -1;
            solver_ab_cancel_phase(executor);
            break;
        }
        if (retire_result == 0) {
            if (solver_ab_poll_phase_stop(
                    solver,
                    next_timer_callback_time)) {
                solver_ab_cancel_phase(executor);
                break;
            }
            if (final) {
                pthread_mutex_lock(&executor->mutex);
                executor->next_reduce++;
                pthread_cond_broadcast(&executor->work_cv);
                pthread_mutex_unlock(&executor->mutex);
                solver_ab_executor_notify_work(executor);
            }
            continue;
        }

        {
            struct timespec deadline;
            anbool executor_cancelled;
            anbool executor_failed;
            anbool segment_ready;
            int wait_result;

            solver_ab_timespec_after_wait(&deadline);
            pthread_mutex_lock(&executor->mutex);
            executor_failed = executor->fatal_error;
            executor_cancelled = __atomic_load_n(
                &executor->cancel_requested,
                __ATOMIC_ACQUIRE);
            segment_ready =
                solver_ab_canonical_segment_ready_locked(
                    executor,
                    executor->next_reduce);
            if (!executor_failed &&
                !executor_cancelled &&
                !segment_ready) {
                wait_result = pthread_cond_timedwait(
                    &executor->done_cv,
                    &executor->mutex,
                    &deadline);
            } else {
                wait_result = 0;
            }
            pthread_mutex_unlock(&executor->mutex);

            if (executor_failed ||
                (executor_cancelled &&
                 !index_shard_worker_stop_requested())) {
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                rc = -1;
                solver_ab_cancel_phase(executor);
                break;
            }
            if (executor_cancelled) {
                (void)solver_poll_worker_stop(solver);
                solver_ab_cancel_phase(executor);
                break;
            }
            if (wait_result == ETIMEDOUT) {
                if (solver_ab_poll_phase_stop(
                        solver,
                        next_timer_callback_time)) {
                    solver_ab_cancel_phase(executor);
                    break;
                }
            } else if (wait_result) {
                solver->profile.execution_failed = TRUE;
                solver->quit_now = TRUE;
                rc = -1;
                solver_ab_cancel_phase(executor);
                break;
            }
        }
    }

    if (!rc && !solver->quit_now) {
        anbool executor_cancelled;
        anbool executor_failed;

        pthread_mutex_lock(&executor->mutex);
        executor_failed = executor->fatal_error;
        executor_cancelled = __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE);
        pthread_mutex_unlock(&executor->mutex);
        if (executor_failed ||
            (executor_cancelled &&
             !index_shard_worker_stop_requested())) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            rc = -1;
            solver_ab_cancel_phase(executor);
        } else if (executor_cancelled) {
            (void)solver_poll_worker_stop(solver);
            solver_ab_cancel_phase(executor);
        }
    }

    pthread_mutex_lock(&executor->mutex);
    while (executor->remaining_tasks > 0U) {
        int wait_result = pthread_cond_wait(
            &executor->done_cv,
            &executor->mutex);

        if (wait_result) {
            size_t pending;

            executor->fatal_error = TRUE;
            __atomic_store_n(
                &executor->cancel_requested,
                TRUE,
                __ATOMIC_RELEASE);
            rc = -1;
            for (pending = executor->next_task;
                 pending < executor->task_count;
                 pending++) {
                if (executor->tasks[pending].state ==
                    SOLVER_AB_TASK_PENDING) {
                    solver_ab_complete_task_locked(
                        executor,
                        pending);
                }
            }
            executor->next_task =
                executor->task_count;
            pthread_cond_broadcast(
                &executor->work_cv);
            pthread_cond_broadcast(
                &executor->done_cv);
        }
    }
    executor->phase_active = FALSE;
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    while (executor->helpers_exited <
           executor->helpers_arrived) {
        int wait_result = pthread_cond_wait(
            &executor->done_cv,
            &executor->mutex);

        if (wait_result) {
            executor->fatal_error = TRUE;
            rc = -1;
            break;
        }
    }
    pthread_mutex_unlock(&executor->mutex);
    solver_ab_apply_latched_failures(
        executor, solver);
    if (solver->profile.execution_failed) {
        rc = -1;
    }
    solver_ab_executor_notify_work(executor);

    {
        int worker_id;

        for (worker_id = 1;
             worker_id < executor->worker_count;
             worker_id++) {
            anbool worker_observed = FALSE;

            for (task_index = 0;
                 task_index < task_count;
                 task_index++) {
                if (tasks[task_index].worker_id ==
                    worker_id) {
                    unsigned long long helper_combinations =
                        tasks[task_index].combination_end -
                        tasks[task_index].combination_first;

                    helper_observed = TRUE;
                    worker_observed = TRUE;
                    solver->profile.ab_helper_tasks++;
                    solver->profile.ab_helper_combinations =
                        solver_ab_saturating_add(
                            solver->profile.ab_helper_combinations,
                            helper_combinations);
                }
            }
            if (worker_observed) {
                participants_observed++;
            }
        }
        solver->profile.max_parallel_ranges =
            MAX(
                solver->profile.max_parallel_ranges,
                participants_observed);
    }
    if (helper_observed) {
        solver->profile.parallel_batches_observed++;
        solver->profile.parallel_hypotheses +=
            solver->profile.hypotheses_executed -
            hypotheses_before;
        if (mode_out) {
            *mode_out = SOLVER_AB_MODE_ASSISTED;
        }
    } else if (mode_out) {
        *mode_out = SOLVER_AB_MODE_FLATTENED_OWNER;
    }

    if (!solver->quit_now && !rc) {
        solver->profile.hypothesis_batches_completed++;
    } else if (rc || solver->profile.execution_failed) {
        solver->profile.hypothesis_batches_failed++;
    } else {
        solver->profile.hypothesis_batches_stopped++;
    }
    solver->profile.hypothesis_wave_wall_seconds +=
        monotonic_seconds() - phase_wall_start;
    solver_ab_free_phase(executor);
    return rc ? -1 : 1;
}

solver_ab_executor_t* solver_ab_executor_create(
    int worker_count) {
    solver_ab_executor_t* executor;

    if (worker_count < 2) {
        return NULL;
    }
    executor = calloc(1, sizeof(*executor));
    if (!executor) {
        return NULL;
    }
    executor->worker_count = worker_count;
    solver_ab_configure_plan_cache_limits(executor);
    executor->query_results = calloc(
        (size_t)worker_count,
        sizeof(*executor->query_results));
    executor->combination_eligible = calloc(
        (size_t)worker_count,
        sizeof(*executor->combination_eligible));
    executor->combination_eligible_capacity = calloc(
        (size_t)worker_count,
        sizeof(*executor->combination_eligible_capacity));
    executor->channels = calloc(
        (size_t)worker_count,
        sizeof(*executor->channels));
    executor->helper_state = calloc(
        (size_t)worker_count,
        sizeof(*executor->helper_state));
    if (!executor->query_results ||
        !executor->combination_eligible ||
        !executor->combination_eligible_capacity ||
        !executor->channels ||
        !executor->helper_state) {
        free(executor->helper_state);
        free(executor->channels);
        free(executor->combination_eligible_capacity);
        free(executor->combination_eligible);
        free(executor->query_results);
        free(executor);
        return NULL;
    }
    if (pthread_mutex_init(&executor->mutex, NULL)) {
        free(executor->helper_state);
        free(executor->channels);
        free(executor->combination_eligible_capacity);
        free(executor->combination_eligible);
        free(executor->query_results);
        free(executor);
        return NULL;
    }
    if (pthread_cond_init(&executor->work_cv, NULL)) {
        pthread_mutex_destroy(&executor->mutex);
        free(executor->helper_state);
        free(executor->channels);
        free(executor->combination_eligible_capacity);
        free(executor->combination_eligible);
        free(executor->query_results);
        free(executor);
        return NULL;
    }
    if (pthread_cond_init(&executor->done_cv, NULL)) {
        pthread_cond_destroy(&executor->work_cv);
        pthread_mutex_destroy(&executor->mutex);
        free(executor->helper_state);
        free(executor->channels);
        free(executor->combination_eligible_capacity);
        free(executor->combination_eligible);
        free(executor->query_results);
        free(executor);
        return NULL;
    }
    return executor;
}

int solver_ab_executor_set_lending_callbacks(
    solver_ab_executor_t* executor,
    void (*work_notify)(void*),
    int (*available_lenders)(void*),
    void* opaque) {
    if (!executor) {
        return -1;
    }
    pthread_mutex_lock(&executor->mutex);
    if (executor->bound ||
        executor->phase_active ||
        executor->stopping ||
        executor->work_notify_inflight) {
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }
    executor->work_notify = work_notify;
    executor->available_lenders = available_lenders;
    executor->work_notify_opaque = opaque;
    pthread_mutex_unlock(&executor->mutex);
    return 0;
}

static anbool solver_ab_executor_begin_notify_locked(
    solver_ab_executor_t* executor,
    void (**notify)(void*),
    void** opaque) {
    if (!executor || !notify || !opaque ||
        !executor->work_notify) {
        return FALSE;
    }
    if (executor->work_notify_inflight == INT_MAX) {
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        return FALSE;
    }
    executor->work_notify_inflight++;
    *notify = executor->work_notify;
    *opaque = executor->work_notify_opaque;
    return TRUE;
}

static void solver_ab_executor_complete_notify(
    solver_ab_executor_t* executor,
    void (*notify)(void*),
    void* opaque) {
    if (!executor || !notify) {
        return;
    }
    notify(opaque);
    pthread_mutex_lock(&executor->mutex);
    if (executor->work_notify_inflight <= 0) {
        logerr("[solver-ab] work notification underflow\n");
        executor->fatal_error = TRUE;
    } else {
        executor->work_notify_inflight--;
    }
    pthread_cond_broadcast(&executor->done_cv);
    pthread_mutex_unlock(&executor->mutex);
}

static void solver_ab_executor_notify_work(
    solver_ab_executor_t* executor) {
    void (*notify)(void*) = NULL;
    void* opaque = NULL;

    if (!executor) {
        return;
    }
    pthread_mutex_lock(&executor->mutex);
    (void)solver_ab_executor_begin_notify_locked(
        executor, &notify, &opaque);
    pthread_mutex_unlock(&executor->mutex);
    if (notify) {
        solver_ab_executor_complete_notify(
            executor, notify, opaque);
    }
}

int solver_ab_executor_bind(
    solver_ab_executor_t* executor,
    solver_t* owner,
    index_t* index) {
    int numxy;
    int worker_id;

    if (!executor || !owner || !owner->fieldxy || !index) {
        return -1;
    }
    numxy = solver_field_geometry_numxy(owner);
    /*
     * A published lane must be able to accept compact tasks.  Geometry
     * refusal/allocation failure and owner-only maxquads/maxmatches semantics
     * are valid native fallbacks, but advertising them as lendable would park
     * queued indexes behind an executor that can never publish work.
     */
    if (numxy <= 0 ||
        !solver_field_geometry_compatible(owner, numxy) ||
        owner->maxquads != 0 ||
        owner->maxmatches != 0) {
        return -1;
    }
    pthread_mutex_lock(&executor->mutex);
    if (executor->bound ||
        executor->stopping ||
        executor->fatal_error ||
        executor->phase_active ||
        executor->pairs ||
        executor->tasks ||
        executor->remaining_tasks != 0U ||
        executor->task_count != 0U ||
        executor->pair_count != 0U ||
        executor->helpers_arrived != 0 ||
        executor->helpers_exited != 0 ||
        executor->helpers_failed != 0 ||
        executor->work_notify_inflight != 0 ||
        executor->pending_allocation_failures != 0U ||
        executor->pending_search_failures != 0U ||
        executor->pending_evaluation_failure) {
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }
    if (!solver_ab_reserve_eligible_workspace(
            executor,
            0,
            (size_t)numxy)) {
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }
    for (worker_id = 1;
         worker_id < executor->worker_count;
         worker_id++) {
        if (executor->helper_state[worker_id] !=
            SOLVER_AB_HELPER_NEW) {
            pthread_mutex_unlock(&executor->mutex);
            return -1;
        }
    }
    executor->owner = owner;
    executor->index = index;
    executor->detailed = owner->profile.detailed;
    executor->join_closed = FALSE;
    executor->starkd_initializing = FALSE;
    executor->starkd_ready = FALSE;
    __atomic_store_n(
        &executor->cancel_requested,
        FALSE,
        __ATOMIC_RELEASE);
    executor->bound = TRUE;
    executor->index_epochs_started++;
    logverb(
        "[solver-ab] starkd-permutation index=%s "
        "present=%i inverse_ready=%i\n",
        index->indexname ? index->indexname : "(unnamed)",
        index->starkd && index->starkd->tree &&
            index->starkd->tree->perm,
        index->starkd && index->starkd->inverse_perm);
    pthread_cond_broadcast(&executor->work_cv);
    pthread_mutex_unlock(&executor->mutex);
    return 0;
}

int solver_ab_executor_try_join(
    solver_ab_executor_t* executor,
    int* worker_id) {
    int slot;
    size_t task_index;

    if (!executor || !worker_id) {
        return -1;
    }
    pthread_mutex_lock(&executor->mutex);
    if (!executor->bound ||
        executor->join_closed ||
        executor->stopping ||
        executor->fatal_error ||
        !executor->phase_active ||
        __atomic_load_n(
            &executor->cancel_requested,
            __ATOMIC_ACQUIRE) ||
        executor->next_task >=
            executor->task_count) {
        pthread_mutex_unlock(&executor->mutex);
        return 0;
    }
    for (slot = 1; slot < executor->worker_count; slot++) {
        if (executor->helper_state[slot] ==
                SOLVER_AB_HELPER_NEW &&
            solver_ab_task_claimable_locked(
                executor,
                slot)) {
            if (executor->helpers_arrived == INT_MAX) {
                executor->fatal_error = TRUE;
                __atomic_store_n(
                    &executor->cancel_requested,
                    TRUE,
                    __ATOMIC_RELEASE);
                pthread_cond_broadcast(
                    &executor->done_cv);
                pthread_cond_broadcast(
                    &executor->work_cv);
                pthread_mutex_unlock(
                    &executor->mutex);
                return -1;
            }
            if (!solver_ab_claim_task_locked(
                    executor,
                    slot,
                    &task_index)) {
                executor->fatal_error = TRUE;
                __atomic_store_n(
                    &executor->cancel_requested,
                    TRUE,
                    __ATOMIC_RELEASE);
                pthread_cond_broadcast(
                    &executor->done_cv);
                pthread_cond_broadcast(
                    &executor->work_cv);
                pthread_mutex_unlock(
                    &executor->mutex);
                return -1;
            }
            executor->helper_state[slot] =
                SOLVER_AB_HELPER_ACTIVE;
            executor->helpers_arrived++;
            *worker_id = slot;
            pthread_cond_broadcast(&executor->done_cv);
            pthread_cond_broadcast(&executor->work_cv);
            pthread_mutex_unlock(&executor->mutex);
            return 1;
        }
    }
    pthread_mutex_unlock(&executor->mutex);
    return 0;
}

int solver_ab_executor_run_joined(
    solver_ab_executor_t* executor,
    int worker_id) {
    size_t task_index = SIZE_MAX;
    void (*notify)(void*) = NULL;
    void* notify_opaque = NULL;
    int rc = 0;

    if (!executor) {
        return -1;
    }
    if (worker_id <= 0 ||
        worker_id >= executor->worker_count) {
        pthread_mutex_lock(&executor->mutex);
        executor->helpers_failed++;
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }

    pthread_mutex_lock(&executor->mutex);
    if (executor->helper_state[worker_id] !=
        SOLVER_AB_HELPER_ACTIVE) {
        /*
         * A non-active slot carries no reservation owned by this call.
         * Do not guess from a stale/default channel index and do not account
         * another helper's exit; either could release phase storage while a
         * real producer is still using it.
         */
        executor->helpers_failed++;
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        pthread_mutex_unlock(&executor->mutex);
        return -1;
    }

    task_index =
        executor->channels[worker_id].task_index;
    if (task_index >= executor->task_count ||
        executor->tasks[task_index].state !=
            SOLVER_AB_TASK_RUNNING ||
        executor->tasks[task_index].worker_id !=
            worker_id ||
        executor->channels[worker_id].ready) {
        executor->helpers_failed++;
        executor->fatal_error = TRUE;
        /*
         * The ACTIVE slot proves that this call owns one reservation.  Close
         * its task only if the task and channel both prove the same ownership;
         * otherwise cancellation drains the corrupted phase without
         * completing an unrelated owner/helper task.
         */
        if (task_index < executor->task_count &&
            executor->channels[worker_id].task_index ==
                task_index &&
            executor->tasks[task_index].state ==
                SOLVER_AB_TASK_RUNNING &&
            executor->tasks[task_index].worker_id ==
                worker_id) {
            solver_ab_complete_task_locked(
                executor, task_index);
        }
        executor->helper_state[worker_id] =
            SOLVER_AB_HELPER_FAILED;
        executor->helpers_exited++;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        rc = -1;
    }
    pthread_mutex_unlock(&executor->mutex);

    if (!rc && task_index != SIZE_MAX &&
        solver_ab_evaluate_task(
            executor,
            task_index,
            worker_id) &&
        !solver_ab_cancelled(executor)) {
        pthread_mutex_lock(&executor->mutex);
        executor->helpers_failed++;
        executor->fatal_error = TRUE;
        __atomic_store_n(
            &executor->cancel_requested,
            TRUE,
            __ATOMIC_RELEASE);
        pthread_cond_broadcast(&executor->done_cv);
        pthread_cond_broadcast(&executor->work_cv);
        pthread_mutex_unlock(&executor->mutex);
        rc = -1;
    }

    pthread_mutex_lock(&executor->mutex);
    if (executor->helper_state[worker_id] ==
        SOLVER_AB_HELPER_FAILED) {
        rc = -1;
    } else if (executor->helper_state[worker_id] !=
        SOLVER_AB_HELPER_ACTIVE) {
        logerr(
            "[solver-ab] invalid leave state for helper %i\n",
            worker_id);
        executor->helpers_failed++;
        executor->fatal_error = TRUE;
        rc = -1;
    } else {
        executor->helper_state[worker_id] =
            SOLVER_AB_HELPER_NEW;
        executor->helpers_exited++;
        (void)solver_ab_executor_begin_notify_locked(
            executor,
            &notify,
            &notify_opaque);
    }
    pthread_cond_broadcast(&executor->done_cv);
    pthread_cond_broadcast(&executor->work_cv);
    if (executor->fatal_error) {
        rc = -1;
    }
    pthread_mutex_unlock(&executor->mutex);
    if (notify) {
        solver_ab_executor_complete_notify(
            executor,
            notify,
            notify_opaque);
    }
    return rc;
}

/*
 * Retire one clean index epoch without destroying tree-neutral capacity.
 *
 * The lane must already be unpublished by the outer scheduler. This closes
 * new joins, waits for every reserved helper and callback to leave, proves
 * the phase is empty, then rearms the shell for another index. KD query
 * results are tree-dependent and are therefore released at every epoch
 * boundary; planner buffers, packet capacity, and eligibility workspaces are
 * retained under their existing limits.
 */
int solver_ab_executor_quiesce(
    solver_ab_executor_t* executor) {
    int worker_id;
    int rc;

    if (!executor) {
        return -1;
    }
    pthread_mutex_lock(&executor->mutex);
    executor->join_closed = TRUE;
    __atomic_store_n(
        &executor->cancel_requested,
        TRUE,
        __ATOMIC_RELEASE);
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    while (executor->helpers_exited <
               executor->helpers_arrived ||
           executor->work_notify_inflight > 0) {
        int wait_result = pthread_cond_wait(
            &executor->done_cv,
            &executor->mutex);

        if (wait_result) {
            executor->fatal_error = TRUE;
            pthread_cond_broadcast(
                &executor->work_cv);
        }
    }
    if (!executor->fatal_error) {
        if (!executor->bound ||
            executor->stopping ||
            executor->phase_active ||
            executor->pairs ||
            executor->tasks ||
            executor->remaining_tasks != 0U ||
            executor->task_count != 0U ||
            executor->pair_count != 0U ||
            executor->starkd_initializing ||
            executor->work_notify_inflight != 0 ||
            executor->pending_allocation_failures != 0U ||
            executor->pending_search_failures != 0U ||
            executor->pending_evaluation_failure) {
            executor->fatal_error = TRUE;
        }
        for (worker_id = 0;
             !executor->fatal_error &&
                 worker_id < executor->worker_count;
             worker_id++) {
            if (executor->channels[worker_id].ready ||
                (worker_id > 0 &&
                 executor->helper_state[worker_id] !=
                     SOLVER_AB_HELPER_NEW)) {
                executor->fatal_error = TRUE;
            }
        }
    }
    rc = executor->fatal_error ? -1 : 0;
    if (rc) {
        logerr(
            "[solver-ab] executor quiesce failed: arrived=%i "
            "exited=%i failed=%i remaining=%zu\n",
            executor->helpers_arrived,
            executor->helpers_exited,
            executor->helpers_failed,
            executor->remaining_tasks);
        /*
         * No owner or helper can still execute after the wait above. Clear
         * stale phase views before the index mapping is released; the caller
         * will destroy this non-reusable shell.
         */
        solver_ab_free_phase(executor);
        executor->phase_active = FALSE;
        executor->remaining_tasks = 0U;
        executor->next_task = 0U;
        executor->next_reduce = 0U;
    }
    for (worker_id = 0;
         worker_id < executor->worker_count;
         worker_id++) {
        kdtree_free_query(
            executor->query_results[worker_id]);
        executor->query_results[worker_id] = NULL;
    }
    executor->bound = FALSE;
    executor->owner = NULL;
    executor->index = NULL;
    memset(
        &executor->snapshot,
        0,
        sizeof(executor->snapshot));
    executor->snapshot_ready = FALSE;
    executor->starkd_initializing = FALSE;
    executor->starkd_ready = FALSE;
    if (!rc) {
        executor->helpers_arrived = 0;
        executor->helpers_exited = 0;
        executor->helpers_failed = 0;
        for (worker_id = 1;
             worker_id < executor->worker_count;
             worker_id++) {
            executor->helper_state[worker_id] =
                SOLVER_AB_HELPER_NEW;
        }
        executor->join_closed = FALSE;
        __atomic_store_n(
            &executor->cancel_requested,
            FALSE,
            __ATOMIC_RELEASE);
        executor->index_epochs_quiesced++;
    }
    logverb(
        "[solver-ab] executor-epoch started=%llu quiesced=%llu "
        "owner_query_reuses=%llu pair_bytes=%zu task_bytes=%zu "
        "pair_grows=%llu pair_reuses=%llu "
        "task_grows=%llu task_reuses=%llu "
        "packet_cache_limit=%zu packet_trims=%llu "
        "eligible_cache_limit=%zu eligible_trims=%llu\n",
        executor->index_epochs_started,
        executor->index_epochs_quiesced,
        executor->owner_query_reuses,
        executor->pair_cache_capacity *
            sizeof(*executor->pair_cache),
        executor->task_cache_capacity *
            sizeof(*executor->task_cache),
        executor->pair_cache_grows,
        executor->pair_cache_reuses,
        executor->task_cache_grows,
        executor->task_cache_reuses,
        executor->packet_cache_limit_bytes,
        executor->packet_cache_trims,
        executor->eligible_cache_limit_bytes,
        executor->eligible_cache_trims);
    pthread_mutex_unlock(&executor->mutex);
    return rc;
}

int solver_ab_executor_finish(
    solver_ab_executor_t* executor) {
    int worker_id;
    int rc;

    if (!executor) {
        return -1;
    }
    pthread_mutex_lock(&executor->mutex);
    executor->join_closed = TRUE;
    executor->stopping = TRUE;
    __atomic_store_n(
        &executor->cancel_requested,
        TRUE,
        __ATOMIC_RELEASE);
    pthread_cond_broadcast(&executor->work_cv);
    pthread_cond_broadcast(&executor->done_cv);
    while (executor->helpers_exited <
               executor->helpers_arrived ||
           executor->work_notify_inflight > 0) {
        int wait_result = pthread_cond_wait(
            &executor->done_cv,
            &executor->mutex);

        if (wait_result) {
            executor->fatal_error = TRUE;
            pthread_cond_broadcast(
                &executor->work_cv);
        }
    }
    if (!executor->fatal_error) {
        if (executor->phase_active ||
            executor->pairs ||
            executor->tasks ||
            executor->remaining_tasks != 0U ||
            executor->task_count != 0U ||
            executor->pair_count != 0U ||
            executor->starkd_initializing ||
            executor->work_notify_inflight != 0 ||
            executor->pending_allocation_failures != 0U ||
            executor->pending_search_failures != 0U ||
            executor->pending_evaluation_failure) {
            executor->fatal_error = TRUE;
        }
        for (worker_id = 0;
             !executor->fatal_error &&
                 worker_id < executor->worker_count;
             worker_id++) {
            if (executor->channels[worker_id].ready ||
                (worker_id > 0 &&
                 executor->helper_state[worker_id] !=
                    SOLVER_AB_HELPER_NEW)) {
                executor->fatal_error = TRUE;
            }
        }
    }
    rc = executor->fatal_error ? -1 : 0;
    if (executor->fatal_error) {
        logerr(
            "[solver-ab] executor finish failed: arrived=%i "
            "exited=%i failed=%i remaining=%zu\n",
            executor->helpers_arrived,
            executor->helpers_exited,
            executor->helpers_failed,
            executor->remaining_tasks);
    }
    logverb(
        "[solver-ab] plan-cache pair_bytes=%zu task_bytes=%zu "
        "pair_limit=%zu task_limit=%zu pair_grows=%llu "
        "pair_reuses=%llu task_grows=%llu task_reuses=%llu "
        "transient_allocations=%llu transient_peak=%zu\n",
        executor->pair_cache_capacity *
            sizeof(*executor->pair_cache),
        executor->task_cache_capacity *
            sizeof(*executor->task_cache),
        executor->pair_cache_limit_bytes,
        executor->task_cache_limit_bytes,
        executor->pair_cache_grows,
        executor->pair_cache_reuses,
        executor->task_cache_grows,
        executor->task_cache_reuses,
        executor->plan_transient_allocations,
        executor->plan_transient_peak_bytes);
    executor->bound = FALSE;
    executor->owner = NULL;
    executor->index = NULL;
    memset(
        &executor->snapshot,
        0,
        sizeof(executor->snapshot));
    executor->snapshot_ready = FALSE;
    executor->starkd_initializing = FALSE;
    executor->starkd_ready = FALSE;
    executor->work_notify = NULL;
    executor->available_lenders = NULL;
    executor->work_notify_opaque = NULL;
    pthread_mutex_unlock(&executor->mutex);
    return rc;
}

void solver_ab_executor_abort(
    solver_ab_executor_t* executor) {
    if (!executor) {
        return;
    }
    pthread_mutex_lock(&executor->mutex);
    executor->fatal_error = TRUE;
    pthread_mutex_unlock(&executor->mutex);
    (void)solver_ab_executor_finish(executor);
}

void solver_ab_executor_destroy(
    solver_ab_executor_t* executor) {
    int i;

    if (!executor) {
        return;
    }
    solver_ab_free_phase(executor);
    free(executor->task_cache);
    free(executor->pair_cache);
    for (i = 0; i < executor->worker_count; i++) {
        kdtree_free_query(executor->query_results[i]);
        free(executor->combination_eligible[i]);
        solver_ab_packet_free(
            &executor->channels[i].packet);
    }
    solver_ab_packet_free(&executor->drain_packet);
    free(executor->helper_state);
    free(executor->channels);
    free(executor->combination_eligible_capacity);
    free(executor->combination_eligible);
    free(executor->query_results);
    pthread_cond_destroy(&executor->done_cv);
    pthread_cond_destroy(&executor->work_cv);
    pthread_mutex_destroy(&executor->mutex);
    free(executor);
}

#if defined(TESTING_TRYPERMUTATIONS)
int solver_test_ab_packet_bounds(void) {
    solver_ab_packet_t packet;
    size_t impossible_candidates =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES /
            sizeof(solver_ab_candidate_t) +
        1U;
    solver_ab_reserve_result_t reserve_result;

    memset(&packet, 0, sizeof(packet));
    if (solver_ab_packet_reserve_hypotheses(
            &packet,
            1U) != SOLVER_AB_RESERVE_OK ||
        solver_ab_packet_reserve_candidates(
            &packet,
            1U) != SOLVER_AB_RESERVE_OK ||
        packet.allocated_bytes >
            SOLVER_AB_PACKET_LIMIT_BYTES) {
        solver_ab_packet_free(&packet);
        return -1;
    }
    reserve_result = solver_ab_packet_reserve_candidates(
        &packet,
        impossible_candidates);
    if (reserve_result != SOLVER_AB_RESERVE_FULL ||
        packet.allocation_failed ||
        packet.allocated_bytes >
            SOLVER_AB_PACKET_LIMIT_BYTES) {
        solver_ab_packet_free(&packet);
        return -1;
    }
    solver_ab_packet_free(&packet);
    return 0;
}

int solver_test_ab_planner_invariants(void) {
    static const unsigned long long totals[] = {
        1U, 2U, 3U, 15U, 16U, 17U
    };
    static const size_t expected_tasks[] = {
        0U, 2U, 3U, 15U, 16U, 16U
    };
    solver_ab_executor_t executor;
    solver_ab_pair_t pair;
    solver_ab_pair_t split_pairs[3];
    solver_ab_task_t* tasks;
    size_t task_count;
    size_t i;
    int plan_result;
    int rc = -1;

    memset(&executor, 0, sizeof(executor));
    executor.worker_count = 4;
    memset(&pair, 0, sizeof(pair));
    for (i = 0U; i < sizeof(totals) / sizeof(totals[0]); i++) {
        pair.combination_count = totals[i];
        tasks = NULL;
        task_count = 0U;
        plan_result = solver_ab_plan_tasks(
            &executor,
            1U,
            totals[i],
            executor.worker_count,
            &tasks,
            &task_count);
        if ((totals[i] == 1U && plan_result != 0) ||
            (totals[i] > 1U && plan_result != 1) ||
            task_count != expected_tasks[i] ||
            (plan_result == 1 &&
             !solver_ab_plan_is_canonical(
                 &pair,
                 1U,
                 tasks,
                 task_count,
                 totals[i]))) {
            goto cleanup;
        }
        solver_ab_free_phase(&executor);
    }

    memset(split_pairs, 0, sizeof(split_pairs));
    split_pairs[0].combination_first = 0U;
    split_pairs[0].combination_count = 1U;
    split_pairs[1].combination_first = 1U;
    split_pairs[1].combination_count = 100U;
    split_pairs[2].combination_first = 101U;
    split_pairs[2].combination_count = 2U;
    tasks = NULL;
    task_count = 0U;
    if (solver_ab_plan_tasks(
            &executor,
            3U,
            103U,
            executor.worker_count,
            &tasks,
            &task_count) != 1 ||
        task_count != 16U ||
        !solver_ab_plan_is_canonical(
            split_pairs,
            3U,
            tasks,
            task_count,
            103U) ||
        solver_ab_count_intra_pair_splits(
            split_pairs,
            3U,
            tasks,
            task_count) != 15U) {
        goto cleanup;
    }
    solver_ab_free_phase(&executor);

    if (solver_ab_publish_threshold(4) != 16U ||
        solver_ab_phase_is_publishable(16U, 4, 0) ||
        !solver_ab_phase_is_publishable(17U, 4, 0) ||
        !solver_ab_phase_is_publishable(2U, 4, 1) ||
        solver_ab_publish_threshold(1024) !=
            SOLVER_AB_INLINE_PUBLISH_MAX) {
        goto cleanup;
    }

    executor.pair_cache_limit_bytes =
        2U * sizeof(solver_ab_pair_t);
    executor.task_cache_limit_bytes =
        16U * sizeof(solver_ab_task_t);
    solver_ab_begin_pair_buffer(&executor);
    if (solver_ab_reserve_pair_buffer(&executor, 2U) ||
        executor.pairs_transient ||
        executor.pair_cache_capacity != 2U ||
        solver_ab_reserve_pair_buffer(&executor, 3U) ||
        !executor.pairs_transient ||
        executor.pair_cache_capacity != 2U) {
        goto cleanup;
    }
    solver_ab_free_phase(&executor);
    solver_ab_begin_pair_buffer(&executor);
    if (executor.pairs != executor.pair_cache ||
        executor.pair_capacity != 2U ||
        executor.pairs_transient) {
        goto cleanup;
    }
    rc = 0;

cleanup:
    solver_ab_free_phase(&executor);
    free(executor.task_cache);
    free(executor.pair_cache);
    return rc;
}

int solver_test_ab_blocks_serial(
    solver_t* solver,
    int worker_count) {
    solver_ab_executor_t executor;
    const solver_field_geometry_t* geometry;
    index_t* index;
    double min_ab;
    double max_ab;
    int numxy;
    int newpoint;
    int worker_id;

    if (!solver ||
        worker_count < 1 ||
        pl_size(solver->indexes) != 1U) {
        return -1;
    }
    numxy = solver_field_geometry_numxy(solver);
    if (!solver_field_geometry_compatible(solver, numxy)) {
        return -1;
    }
    geometry = solver->field_geometry;
    index = pl_get(solver->indexes, 0);
    solver_compute_quad_range(
        solver,
        index,
        &min_ab,
        &max_ab);
    solver->minminAB2 =
        MAX(square(min_ab), square(solver->quadsize_min));
    solver->maxmaxAB2 = square(max_ab);
    if (solver->quadsize_max != 0.0) {
        solver->maxmaxAB2 =
            MIN(
                solver->maxmaxAB2,
                square(solver->quadsize_max));
    }
    if (index->cx_less_than_dx) {
        solver->cxdx_margin =
            1.5 * solver->codetol;
    }
    set_index(solver, index);

    memset(&executor, 0, sizeof(executor));
    executor.worker_count = worker_count;
    executor.query_results = calloc(
        (size_t)worker_count,
        sizeof(*executor.query_results));
    executor.combination_eligible = calloc(
        (size_t)worker_count,
        sizeof(*executor.combination_eligible));
    executor.combination_eligible_capacity = calloc(
        (size_t)worker_count,
        sizeof(*executor.combination_eligible_capacity));
    if (!executor.query_results ||
        !executor.combination_eligible ||
        !executor.combination_eligible_capacity) {
        free(executor.combination_eligible_capacity);
        free(executor.combination_eligible);
        free(executor.query_results);
        return -1;
    }
    executor.owner = solver;
    executor.index = index;
    executor.field_geometry = geometry;
    executor.dimquads = index_dimquads(index);
    executor.min_ab2 = square(min_ab);
    executor.max_ab2 = square(max_ab);
    if (solver_ab_capture_snapshot(&executor, solver)) {
        free(executor.combination_eligible_capacity);
        free(executor.combination_eligible);
        free(executor.query_results);
        return -1;
    }
    for (worker_id = 0;
         worker_id < worker_count;
         worker_id++) {
        if (!solver_ab_reserve_eligible_workspace(
                &executor,
                worker_id,
                (size_t)numxy)) {
            goto fail;
        }
    }

    for (newpoint = solver->startobj;
         newpoint < numxy;
         newpoint++) {
        solver_ab_phase_kind_t phase;

        solver->last_examined_object = newpoint;
        for (phase = SOLVER_AB_PHASE_DIAGONAL;
             phase <= SOLVER_AB_PHASE_OFF_DIAGONAL;
             phase++) {
            solver_ab_pair_t* pairs = NULL;
            solver_ab_task_t* tasks = NULL;
            size_t pair_count = 0U;
            size_t task_count = 0U;
            unsigned long long total_combinations = 0U;
            int plan_result;
            size_t task_index;
            solver_ab_reduce_state_t reduce_state;

            if (solver_ab_collect_pairs(
                    &executor,
                    phase,
                    newpoint,
                    executor.dimquads,
                    executor.field_geometry,
                    executor.min_ab2,
                    executor.max_ab2,
                    &pairs,
                    &pair_count,
                    &total_combinations)) {
                goto fail;
            }
            if (!pair_count) {
                solver_ab_free_phase(&executor);
                continue;
            }
            plan_result = solver_ab_plan_tasks(
                &executor,
                pair_count,
                total_combinations,
                worker_count,
                &tasks,
                &task_count);
            if (plan_result < 0) {
                goto fail;
            }
            if (!plan_result) {
                if (solver_ab_reserve_task_buffer(
                        &executor, 1U)) {
                    goto fail;
                }
                tasks = executor.tasks;
                memset(tasks, 0, sizeof(*tasks));
                tasks[0].combination_first = 0U;
                tasks[0].combination_end =
                    total_combinations;
                task_count = 1U;
            }
            if (!solver_ab_plan_is_canonical(
                    pairs,
                    pair_count,
                    tasks,
                    task_count,
                    total_combinations)) {
                goto fail;
            }

            executor.newpoint = newpoint;
            executor.phase = phase;
            executor.pair_count = pair_count;
            executor.task_count = task_count;
            memset(&reduce_state, 0, sizeof(reduce_state));
            for (task_index = 0;
                 task_index < task_count;
                 task_index++) {
                solver_ab_packet_t packet;
                solver_ab_builder_t builder;
                int visit_result;

                memset(&packet, 0, sizeof(packet));
                memset(&builder, 0, sizeof(builder));
                builder.executor = &executor;
                builder.packet = &packet;
                builder.query_result =
                    &executor.query_results[0];
                visit_result = solver_ab_walk_task_ranges(
                    &executor,
                    &tasks[task_index],
                    0,
                    solver_ab_builder_visit_pair_range,
                    &builder);
                if (visit_result < 0) {
                    builder.fatal_error = TRUE;
                    packet.evaluation_failed = TRUE;
                }
                solver_ab_stage_pending_deltas(&builder);
                if (builder.fatal_error ||
                    packet.cancelled ||
                    solver_ab_reduce_packet(
                        &executor,
                        &packet,
                        &reduce_state) ||
                    reduce_state.in_hypothesis) {
                    solver_ab_packet_free(&packet);
                    solver_ab_free_phase(&executor);
                    goto fail;
                }
                solver_ab_packet_free(&packet);
            }
            solver_ab_free_phase(&executor);
        }
    }
    for (worker_id = 0;
         worker_id < worker_count;
         worker_id++) {
        free(executor.combination_eligible[worker_id]);
    }
    free(executor.combination_eligible_capacity);
    free(executor.combination_eligible);
    free(executor.query_results);
    free(executor.task_cache);
    free(executor.pair_cache);
    return 0;

fail:
    solver_ab_free_phase(&executor);
    for (worker_id = 0;
         worker_id < worker_count;
         worker_id++) {
        free(executor.combination_eligible[worker_id]);
    }
    free(executor.combination_eligible_capacity);
    free(executor.combination_eligible);
    free(executor.query_results);
    free(executor.task_cache);
    free(executor.pair_cache);
    return -1;
}
#endif


// The real deal
int solver_run(solver_t* solver) {
    int numxy, newpoint;
    double usertime, systime;
    double run_wall_start;
    // first timer callback is called after 1 second
    time_t next_timer_callback_time = time(NULL) + 1;
    pquad* pquads = NULL;
    const solver_field_geometry_t* field_geometry = NULL;
    anbool shared_field_geometry = FALSE;
    size_t i, num_indexes;
    double tol2;
    int field[DQMAX];
    // Reuse CodeKD result capacity across all hypotheses in this run.
    kdtree_qres_t* codekd_result = NULL;

    run_wall_start = monotonic_seconds();
    memset(&solver->profile, 0, sizeof(solver->profile));
    solver->profile.detailed = index_shard_trace_enabled();

    get_resource_stats(&usertime, &systime, NULL);

    if (!solver->vf)
        solver_preprocess_field(solver);

    memset(field, 0, sizeof(field));

    solver->starttime = usertime + systime;

    numxy = starxy_n(solver->fieldxy);
    if (solver->endobj && (numxy > solver->endobj))
        numxy = solver->endobj;
    if (solver->startobj >= numxy) {
        solver->profile.solver_run_wall_seconds =
            monotonic_seconds() - run_wall_start;
        solver_profile_report(solver);
        return 0;
    }


    if (numxy >= 1000) {
        logverb("Limiting search to first 1000 objects\n");
        numxy = 1000;
    }

    if (solver->field_geometry &&
        !solver_field_geometry_compatible(
            solver, numxy)) {
        solver_release_field_geometry(solver);
    }
    num_indexes = pl_size(solver->indexes);
    if (solver_field_geometry_compatible(solver, numxy)) {
        field_geometry = solver->field_geometry;
        /*
         * Compact geometry is advisory for index-free descriptor work. Only
         * the dormant flattened executor replaces native pquad state with
         * it; descriptor and fallback paths retain owner-private pquads.
         */
        shared_field_geometry =
            solver->ab_executor &&
            num_indexes == 1U &&
            solver->maxquads == 0 &&
            solver->maxmatches == 0;
        if (shared_field_geometry) {
            __atomic_fetch_add(
                &solver->field_geometry->reused_solver_runs,
                1,
                __ATOMIC_RELAXED);
        }
    }
    {
        double minAB2s[num_indexes];
        double maxAB2s[num_indexes];
        solver->minminAB2 = LARGE_VAL;
        solver->maxmaxAB2 = -LARGE_VAL;
        for (i = 0; i < num_indexes; i++) {
            double minAB=0, maxAB=0;
            index_t* index = pl_get(solver->indexes, i);
            // The limits on the size of quads that we try to match, in pixels.
            // Derived from index_scale_* and funits_*.
            solver_compute_quad_range(solver, index, &minAB, &maxAB);
            //logverb("Index \"%s\" quad range %f to %f\n", index->indexname,
            //minAB, maxAB);
            minAB2s[i] = square(minAB);
            maxAB2s[i] = square(maxAB);
            solver->minminAB2 = MIN(solver->minminAB2, minAB2s[i]);
            solver->maxmaxAB2 = MAX(solver->maxmaxAB2, maxAB2s[i]);

            if (index->cx_less_than_dx) {
                solver->cxdx_margin = 1.5 * solver->codetol;
                // FIXME die horribly if the indexes have differing cx_less_than_dx
            }
        }
        solver->minminAB2 = MAX(solver->minminAB2, square(solver->quadsize_min));
        if (solver->quadsize_max != 0.0)
            solver->maxmaxAB2 = MIN(solver->maxmaxAB2, square(solver->quadsize_max));
        logverb("Quad scale range: [%g, %g] pixels\n", sqrt(solver->minminAB2), sqrt(solver->maxmaxAB2));

        // quick-n-dirty scale estimate using stars A,B.
        solver->abscale_high = square(arcsec2rad(solver->funits_upper) * (1.0 + solver->codetol));
        solver->abscale_low  = square(arcsec2rad(solver->funits_lower) * (1.0 - solver->codetol));

        /** Ugh, I want to avoid doing distsq2rad when checking scale,
         but that means correcting for the difference between the
         distance along the curve of the sphere vs the chord distance.
         This affects the lower bound for the largest quads in a messy way...
         This below isn't right.
         solver->abscale_high = square(arcsec2rad(solver->funits_upper) * (1.0 + solver->codetol));
         solver->abscale_low = arcsec2rad(solver->funits_lower) * (1.0 - solver->codetol) *
         MIN(M_PI, arcsec2rad(field_diag * solver->funits_upper)) ...
         */

        if (!shared_field_geometry) {
            pquads = calloc(
                (size_t)numxy * (size_t)numxy,
                sizeof(pquad));
            if (!pquads) {
                SYSERROR("Failed to allocate solver pquad array");
                solver->profile.execution_failed = TRUE;
                goto finish;
            }
        }

        /* We maintain an array of "potential quads" (pquad) structs, where
         * each struct corresponds to one choice of stars A and B; the struct
         * at index (B * numxy + A) holds information about quads that could be
         * created using stars A,B.
         *
         * (We only use the above-diagonal elements of this 2D array because
         * A<B.)
         *
         * For each AB pair, we cache the scale and the rotation parameters,
         * and we keep an array "inbox" of length "numxy" of booleans, one for
         * each star, which say whether that star is eligible to be star C or D
         * of a quad with AB at the corners.  (Obviously A and B aren't
         * eligible).
         *
         * The "ninbox" parameter is somewhat misnamed - it says that "inbox"
         * elements in the range [0, ninbox) have been initialized.
         */

        /* (See explanatory paragraph below) If "solver->startobj" isn't zero,
         * then we need to initialize the triangle of "pquads" up to
         * A=startobj-2, B=startobj-1. */
        if (!shared_field_geometry && solver->startobj) {
            debug("startobj > 0; priming pquad arrays.\n");
            for (field[B] = 0; field[B] < solver->startobj; field[B]++) {
                for (field[A] = 0; field[A] < field[B]; field[A]++) {
                    pquad* pq = pquads + field[B] * numxy + field[A];
                    pq->fieldA = field[A];
                    pq->fieldB = field[B];
                    debug("trying A=%i, B=%i\n", field[A], field[B]);
                    check_scale(pq, solver);
                    if (!pq->scale_ok) {
                        debug("  bad scale for A=%i, B=%i\n", field[A], field[B]);
                        continue;
                    }
                    if (solver_allocate_pquad_storage(
                            solver, pq, numxy)) {
                        goto quitnow;
                    }
                    memset(pq->inbox, TRUE, solver->startobj);
                    pq->ninbox = solver->startobj;
                    pq->inbox[field[A]] = FALSE;
                    pq->inbox[field[B]] = FALSE;
                    check_inbox(pq, 0, solver);
                    debug("  inbox(A=%i, B=%i): ", field[A], field[B]);
                    print_inbox(pq);
                }
            }
        }

        /* Each time through the "for" loop below, we consider a new star
         * ("newpoint").  First, we try building all quads that have the new
         * star on the diagonal (star B).  Then, we try building all quads that
         * have the star not on the diagonal (star D).
         *
         * For each AB pair, we have a "potential_quad" or "pquad" struct.
         * This caches the computation we need to do: deciding whether the
         * scale is acceptable, computing the transformation to code
         * coordinates, and deciding which C,D stars are in the circle.
         */
        for (newpoint = solver->startobj; newpoint < numxy; newpoint++) {

            debug("Trying newpoint=%i (%.1f,%.1f)\n", newpoint,
                  field_getx(solver,newpoint), field_gety(solver,newpoint));

            if (solver_poll_worker_stop(solver)) {
                break;
            }

            // Give our caller a chance to cancel us midway. The callback
            // returns how long to wait before calling again.

            if (solver->timer_callback) {
                time_t delay;
                time_t now = time(NULL);
                if (now > next_timer_callback_time) {
                    update_timeused(solver);
                    delay = solver->timer_callback(solver->userdata);
                    if (delay == 0) // Canceled
                        break;
                    next_timer_callback_time = now + delay;
                }
            }

            solver->last_examined_object = newpoint;
            // quads with the new star on the diagonal:
            field[B] = newpoint;
            debug("Trying quads with B=%i\n", newpoint);

            // first do an index-independent scale check...
            if (!shared_field_geometry) {
                for (field[A] = 0; field[A] < newpoint; field[A]++) {
                    // initialize the "pquad" struct for this AB combo.
                    pquad* pq = pquads + field[B] * numxy + field[A];
                    pq->fieldA = field[A];
                    pq->fieldB = field[B];
                    debug("  trying A=%i, B=%i\n", field[A], field[B]);
                    check_scale(pq, solver);
                    if (!pq->scale_ok) {
                        debug("    bad scale for A=%i, B=%i\n", field[A], field[B]);
                        continue;
                    }
                    // initialize the "inbox" array:
                    if (solver_allocate_pquad_storage(
                            solver, pq, numxy)) {
                        goto quitnow;
                    }
                    // -try all stars up to "newpoint"...
                    assert(sizeof(anbool) == 1);
                    memset(pq->inbox, TRUE, newpoint + 1);
                    pq->ninbox = newpoint + 1;
                    // -except A and B.
                    pq->inbox[field[A]] = FALSE;
                    pq->inbox[field[B]] = FALSE;
                    check_inbox(pq, 0, solver);
                    debug("    inbox(A=%i, B=%i): ", field[A], field[B]);
                    print_inbox(pq);
                }
            }

            // Now iterate through the different indices
            for (i = 0; i < num_indexes; i++) {
                index_t* index = pl_get(solver->indexes, i);
                int dimquads;
                int ab_phase_result = 0;
                solver_ab_phase_mode_t phase_mode =
                    SOLVER_AB_MODE_NATIVE;
                solver_ab_phase_telemetry_t phase_telemetry;

                set_index(solver, index);
                dimquads = index_dimquads(index);
                solver_ab_phase_telemetry_begin(
                    solver,
                    &phase_telemetry);
                if (num_indexes == 1U &&
                    !solver->ab_executor &&
                    field_geometry) {
                    ab_phase_result =
                        solver_ab_descriptor_execute_phase(
                            solver,
                            SOLVER_AB_PHASE_DIAGONAL,
                            newpoint,
                            field_geometry,
                            dimquads,
                            minAB2s[i],
                            maxAB2s[i],
                            &codekd_result,
                            &phase_mode);
                    if (ab_phase_result < 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_DIAGONAL,
                            phase_mode);
                        goto quitnow;
                    }
                    if (ab_phase_result > 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_DIAGONAL,
                            phase_mode);
                        continue;
                    }
                }
                if (num_indexes == 1U &&
                    solver->ab_executor) {
                    ab_phase_result = solver_ab_execute_phase(
                        solver,
                        SOLVER_AB_PHASE_DIAGONAL,
                        newpoint,
                        field_geometry,
                        dimquads,
                        minAB2s[i],
                        maxAB2s[i],
                        &next_timer_callback_time,
                        &phase_mode);
                    if (ab_phase_result < 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_DIAGONAL,
                            phase_mode);
                        goto quitnow;
                    }
                    if (shared_field_geometry &&
                        ab_phase_result == 0) {
                        logerr("[solver-ab] compact diagonal phase "
                               "declined after geometry commit\n");
                        solver->profile.execution_failed = TRUE;
                        solver->quit_now = TRUE;
                        goto quitnow;
                    }
                    if (ab_phase_result > 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_DIAGONAL,
                            phase_mode);
                        continue;
                    }
                }
                for (field[A] = 0; field[A] < newpoint; field[A]++) {
                    // initialize the "pquad" struct for this AB combo.
                    pquad* pq = pquads + field[B] * numxy + field[A];
                    if (!pq->scale_ok)
                        continue;
                    if ((pq->scale < minAB2s[i]) ||
                        (pq->scale > maxAB2s[i]))
                        continue;
                    // set code tolerance for this index and AB pair...
                    solver->rel_field_noise2 = pq->rel_field_noise2;
                    tol2 = get_tolerance(solver);
                    // Now look at all sets of (C, D, ...) stars (subject to field[C] < field[D] < ...)
                    // ("dimquads - 2" because we've set stars A and B at this point)
                    add_stars(pq,
                              field,
                              C,
                              dimquads - 2,
                              0,
                              newpoint,
                              dimquads,
                              solver,
                              tol2,
                              &codekd_result);
                    if (solver->quit_now) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_DIAGONAL,
                            SOLVER_AB_MODE_NATIVE);
                        goto quitnow;
                    }
                }
                solver_ab_phase_telemetry_report(
                    solver,
                    &phase_telemetry,
                    newpoint,
                    SOLVER_AB_PHASE_DIAGONAL,
                    SOLVER_AB_MODE_NATIVE);
            }

            if (solver->quit_now)
                goto quitnow;

            // Now try building quads with the new star not on the diagonal:
            field[C] = newpoint;
            // (in this loop field[C] > field[D])
            debug("Trying quads with C=%i\n", newpoint);
            {
                int ab_phase_result = 0;
                solver_ab_phase_mode_t phase_mode =
                    SOLVER_AB_MODE_NATIVE;
                anbool phase_pquads_prepared =
                    shared_field_geometry;
                solver_ab_phase_telemetry_t phase_telemetry;

                if (num_indexes == 1U) {
                    set_index(
                        solver,
                        pl_get(solver->indexes, 0));
                }
                solver_ab_phase_telemetry_begin(
                    solver,
                    &phase_telemetry);

                if (num_indexes == 1U &&
                    !solver->ab_executor &&
                    field_geometry) {
                    index_t* index = pl_get(solver->indexes, 0);
                    int dimquads;

                    set_index(solver, index);
                    dimquads = index_dimquads(index);
                    ab_phase_result =
                        solver_ab_descriptor_execute_phase(
                            solver,
                            SOLVER_AB_PHASE_OFF_DIAGONAL,
                            newpoint,
                            field_geometry,
                            dimquads,
                            minAB2s[0],
                            maxAB2s[0],
                            &codekd_result,
                            &phase_mode);
                    if (ab_phase_result < 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_OFF_DIAGONAL,
                            phase_mode);
                        goto quitnow;
                    }
                    /*
                     * A later phase may need the native fallback path.
                     * Preserve its incremental pquad state only after the
                     * index-free descriptor path accepted this phase.
                     */
                    if (ab_phase_result > 0 &&
                        !solver->quit_now &&
                        !phase_pquads_prepared) {
                        for (field[A] = 0;
                             field[A] < newpoint;
                             field[A]++) {
                            if (solver_ab_poll_phase_stop(
                                    solver,
                                    &next_timer_callback_time)) {
                                goto quitnow;
                            }
                            for (field[B] = field[A] + 1;
                                 field[B] < newpoint;
                                 field[B]++) {
                                pquad* pq =
                                    pquads +
                                    field[B] * numxy +
                                    field[A];

                                if (!pq->scale_ok) {
                                    continue;
                                }
                                pq->inbox[field[C]] = TRUE;
                                pq->ninbox = field[C] + 1;
                                check_inbox(
                                    pq,
                                    field[C],
                                    solver);
                            }
                        }
                        phase_pquads_prepared = TRUE;
                    }
                }
                if (num_indexes == 1U &&
                    solver->ab_executor) {
                    index_t* index = pl_get(solver->indexes, 0);
                    int dimquads;

                    set_index(solver, index);
                    dimquads = index_dimquads(index);
                    /*
                     * Off-diagonal AB assistance consumes the current C
                     * membership as an immutable phase snapshot. Shared
                     * whole-field geometry already contains it. The native
                     * owner-private pquad path prepares the same rows once
                     * before helpers are exposed, so geometry allocation is
                     * a cache decision rather than a scheduler gate.
                     */
                    if (!phase_pquads_prepared) {
                        for (field[A] = 0;
                             field[A] < newpoint;
                             field[A]++) {
                            if (solver_ab_poll_phase_stop(
                                    solver,
                                    &next_timer_callback_time)) {
                                goto quitnow;
                            }
                            for (field[B] = field[A] + 1;
                                 field[B] < newpoint;
                                 field[B]++) {
                                pquad* pq =
                                    pquads +
                                    field[B] * numxy +
                                    field[A];

                                if (!pq->scale_ok) {
                                    continue;
                                }
                                pq->inbox[field[C]] = TRUE;
                                pq->ninbox = field[C] + 1;
                                check_inbox(
                                    pq,
                                    field[C],
                                    solver);
                            }
                        }
                        phase_pquads_prepared = TRUE;
                    }
                    ab_phase_result = solver_ab_execute_phase(
                        solver,
                        SOLVER_AB_PHASE_OFF_DIAGONAL,
                        newpoint,
                        field_geometry,
                        dimquads,
                        minAB2s[0],
                        maxAB2s[0],
                        &next_timer_callback_time,
                        &phase_mode);
                    if (ab_phase_result < 0) {
                        solver_ab_phase_telemetry_report(
                            solver,
                            &phase_telemetry,
                            newpoint,
                            SOLVER_AB_PHASE_OFF_DIAGONAL,
                            phase_mode);
                        goto quitnow;
                    }
                    if (shared_field_geometry &&
                        ab_phase_result == 0) {
                        logerr("[solver-ab] compact off-diagonal phase "
                               "declined after geometry commit\n");
                        solver->profile.execution_failed = TRUE;
                        solver->quit_now = TRUE;
                        goto quitnow;
                    }
                }
                if (ab_phase_result <= 0) {
                    for (field[A] = 0;
                         field[A] < newpoint;
                         field[A]++) {
                        for (field[B] = field[A] + 1;
                             field[B] < newpoint;
                             field[B]++) {
                            // grab the "pquad" for this AB combo
                            pquad* pq =
                                pquads +
                                field[B] * numxy +
                                field[A];
                            if (!pq->scale_ok) {
                                debug("  bad scale for A=%i, B=%i\n",
                                      field[A],
                                      field[B]);
                                continue;
                            }
                            // test if this C is in the box:
                            if (!phase_pquads_prepared) {
                                pq->inbox[field[C]] = TRUE;
                                pq->ninbox = field[C] + 1;
                                check_inbox(pq, field[C], solver);
                            }
                            if (!pq->inbox[field[C]]) {
                                debug("  C is not in the box for A=%i, "
                                      "B=%i\n",
                                      field[A],
                                      field[B]);
                                continue;
                            }
                            debug("  C is in the box for A=%i, B=%i\n",
                                  field[A],
                                  field[B]);
                            debug("    box now:");
                            print_inbox(pq);
                            debug("\n");

                            solver->rel_field_noise2 =
                                pq->rel_field_noise2;

                            for (i = 0;
                                 i < pl_size(solver->indexes);
                                 i++) {
                                int dimquads;
                                index_t* index =
                                    pl_get(solver->indexes, i);
                                if ((pq->scale < minAB2s[i]) ||
                                    (pq->scale > maxAB2s[i])) {
                                    continue;
                                }

                                set_index(solver, index);
                                dimquads =
                                    index_dimquads(index);

                                tol2 = get_tolerance(solver);

                                if (dimquads > 3) {
                                    /*
                                     * dimquads - 3 because A, B, and C
                                     * are already fixed.
                                     */
                                    add_stars(
                                        pq,
                                        field,
                                        D,
                                        dimquads - 3,
                                        0,
                                        newpoint,
                                        dimquads,
                                        solver,
                                        tol2,
                                        &codekd_result);
                                } else {
                                    TRY_ALL_CODES(
                                        pq,
                                        field,
                                        dimquads,
                                        solver,
                                        tol2,
                                        &codekd_result);
                                }
                                if (solver->quit_now) {
                                    solver_ab_phase_telemetry_report(
                                        solver,
                                        &phase_telemetry,
                                        newpoint,
                                        SOLVER_AB_PHASE_OFF_DIAGONAL,
                                        SOLVER_AB_MODE_NATIVE);
                                    goto quitnow;
                                }
                            }
                        }
                    }
                }
                solver_ab_phase_telemetry_report(
                    solver,
                    &phase_telemetry,
                    newpoint,
                    SOLVER_AB_PHASE_OFF_DIAGONAL,
                    phase_mode);
            }

            logverb("object %u of %u: %i quads tried, %i matched.\n",
                    newpoint + 1, numxy, solver->numtries, solver->nummatches);

            if ((solver->maxquads && (solver->numtries >= solver->maxquads))
                || (solver->maxmatches && (solver->nummatches >= solver->maxmatches))
                || solver->quit_now)
                break;
        }

    quitnow:
        kdtree_free_query(codekd_result);

        if (!shared_field_geometry) {
            for (i = 0; i < (numxy*numxy); i++) {
                pquad* pq = pquads + i;
                free(pq->inbox);
                free(pq->xy);
            }
            free(pquads);
        }
    }

finish:
    solver->profile.solver_run_wall_seconds =
        monotonic_seconds() - run_wall_start;
    solver_profile_report(solver);
    return solver->profile.execution_failed ? -1 : 0;
}

/**
 All the stars in this quad have been chosen.  Figure out which
 permutations of stars CDE are valid and search for matches.
 */
static void try_all_codes(const pquad* pq,
                          const int* fieldstars, int dimquad,
                          solver_t* solver, double tol2,
                          kdtree_qres_t** presult) {
    int dimcode = (dimquad - 2) * 2;
    double code[DCMAX];
    double flipcode[DCMAX];
    int i;

    solver->numtries++;

    debug("  trying quad [");
    for (i=0; i<dimquad; i++) {
        debug("%s%i", (i?" ":""), fieldstars[i]);
    }
    debug("]\n");

    for (i=0; i<dimquad-NBACK; i++) {
        code[2*i  ] = getx(pq->xy, fieldstars[NBACK+i]);
        code[2*i+1] = gety(pq->xy, fieldstars[NBACK+i]);
    }

    if (solver->parity == PARITY_NORMAL ||
        solver->parity == PARITY_BOTH) {

        debug("    trying normal parity: code=[");
        for (i=0; i<dimcode; i++)
            debug("%s%g", (i?", ":""), code[i]);
        debug("].\n");

        try_all_codes_2(fieldstars, dimquad, code, solver, FALSE,
                        tol2, presult);
    }

    if (unlikely(solver->quit_now))
        return;

    if (solver->parity == PARITY_FLIP ||
        solver->parity == PARITY_BOTH) {

        quad_flip_parity(code, flipcode, dimcode);

        debug("    trying reverse parity: code=[");
        for (i=0; i<dimcode; i++)
            debug("%s%g", (i?", ":""), flipcode[i]);
        debug("].\n");

        try_all_codes_2(fieldstars, dimquad, flipcode, solver, TRUE,
                        tol2, presult);
    }
}

/**
 This function tries the quad with the "backbone" stars A and B in
 normal and flipped configurations.
 */
static void try_all_codes_2(const int* fieldstars, int dimquad,
                            const double* code, solver_t* solver,
                            anbool current_parity, double tol2,
                            kdtree_qres_t** presult) {
    int i;
    int dimcode = (dimquad - NBACK) * 2;
    int stars[DQMAX];
    double flipcode[DCMAX];
    anbool placed[DQMAX];

    if (unlikely(solver->quit_now))
        return;

    // We actually only use elements up to dimquads-2.

    // Un-flipped:
    stars[0] = fieldstars[0];
    stars[1] = fieldstars[1];

    for (i=0; i<DQMAX; i++)
        placed[i] = FALSE;

    try_permutations(fieldstars, dimquad, code, solver,
                     current_parity, tol2, stars, NULL, 0,
                     placed, presult);

    if (unlikely(solver->quit_now))
        return;

    // Flipped:
    stars[0] = fieldstars[1];
    stars[1] = fieldstars[0];

    for (i=0; i<dimcode; i++)
        flipcode[i] = 1.0 - code[i];

    for (i=0; i<DQMAX; i++)
        placed[i] = FALSE;

    try_permutations(fieldstars, dimquad, flipcode, solver,
                     current_parity, tol2, stars, NULL, 0,
                     placed, presult);
}

/**
 This functions tries different permutations of the non-backbone
 stars C [, D [,E ] ]

 origstars: [0] and [1] are the "backbone" stars

 stars: only elements [0] and [1] are set; they will be equal to origstars [0],[1] or [1],[0].
 code: may be NULL, in which case use a local variable
 slot: 0 on initial call; incremented on recursive calls
 */
static void try_permutations(const int* origstars, int dimquad,
                             const double* origcode,
                             solver_t* solver, anbool current_parity,
                             double tol2,
                             int* stars, double* code,
                             int slot, anbool* placed,
                             kdtree_qres_t** presult) {
    int i;
    double mycode[DCMAX];
    int Nstars = dimquad - NBACK;
    int lastslot = dimquad - NBACK - 1;
    /*
     This is a recursive function that tries all combinations of the
     "internal" stars (ie, not stars A,B that form the "backbone" of
     the quad).

     We fill the "stars" array with the star IDs (from "origstars") of
     the stars that form the quad, while simultaneously filling the
     "code" array with the corresponding code coordinates (from
     "origcode").

     For example, if "dimquad" is 5, and "origstars" contains
     A,B,C,D,E, we want to call "resolve_matches" with the following
     combinations in "stars":

     AB CDE
     AB CED
     AB DCE
     AB DEC
     AB ECD
     AB EDC

     This call will try to put each star in "slot" in turn, then for
     each one recurse to "slot" in the rest of the stars.

     Note that we are filling stars[2], stars[3], etc; the first two
     elements are already filled by stars A and B.
     */

    if (code == NULL)
        code = mycode;

    // try to convince the compiler that this is okay
    if (slot >= DCMAX/2)
        return;

    // We try putting each star that hasn't already been placed in
    // this "slot".
    for (i=0; i<Nstars; i++) {
        if (placed[i])
            continue;

        // Check cx <= dx, if we're a "dx".
        if (slot > 0 && solver->index->cx_less_than_dx) {
            if (code[2*(slot - 1) +0] >
                origcode[2*i +0] + solver->cxdx_margin) {
                debug("cx <= dx check failed: %g > %g + %g\n",
                      code[2*(slot - 1) +0], origcode[2*i +0],
                      solver->cxdx_margin);
                solver->num_cxdx_skipped++;
                continue;
            }
        }

        // Slot in this star...
        stars[slot + NBACK] = origstars[i + NBACK];
        code[2*slot +0] = origcode[2*i +0];
        code[2*slot +1] = origcode[2*i +1];

        // Check meanx <= 1/2.
        if (solver->index->cx_less_than_dx &&
            solver->index->meanx_less_than_half) {
            // Check the "cx + dx <= 1" condition (for quads); in general,
            // combined with the "cx <= dx" condition, this means that the
            // mean(x) <= 1/2.
            int j;
            double meanx = 0;
            for (j=0; j<=slot; j++)
                meanx += code[2*j];
            meanx /= (slot+1);
            if (meanx > 0.5 + solver->cxdx_margin) {
                debug("meanx <= 0.5 check failed: %g > 0.5 + %g\n",
                      meanx, solver->cxdx_margin);
                solver->num_meanx_skipped++;
                continue;
            }
        }

        // If we have more slots to fill...
        if (slot < lastslot) {
            placed[i] = TRUE;
            try_permutations(origstars, dimquad, origcode, solver,
                             current_parity, tol2, stars, code,
                             slot+1, placed, presult);
            placed[i] = FALSE;

            if (unlikely(solver->quit_now))
                return;
        } else {
#if defined(TESTING_TRYPERMUTATIONS)
            TEST_TRY_PERMUTATIONS(stars, code, dimquad, solver);
            continue;
#endif
            solver_execute_hypothesis_owner(
                stars,
                code,
                dimquad,
                solver,
                current_parity,
                tol2,
                presult);
            if (unlikely(solver->quit_now))
                return;
        }
    }
}

static void solver_execute_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    kdtree_qres_t** presult) {
    int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;
    double search_wall_start = 0.0;

    if (solver_poll_worker_stop(solver)) {
        return;
    }

    // Search with the code we've built.
    solver->profile.hypotheses_generated++;
    if (solver->profile.max_batch_hypotheses < 1)
        solver->profile.max_batch_hypotheses = 1;

    if (solver->profile.detailed) {
        uint64_t hypothesis_digest =
            solver_hypothesis_order_digest(
                stars,
                code,
                dimquad,
                current_parity);

        solver->profile.hypothesis_order_hash =
            solver_order_hash_mix(
                solver->profile.hypothesis_order_hash,
                hypothesis_digest);
        search_wall_start = monotonic_seconds();
    }

    *presult = solver_codekd_rangesearch(
        solver->index->codekd->tree,
        *presult,
        code,
        tol2,
        options);

    solver->profile.codekd_calls++;
    solver->profile.hypotheses_executed++;

    if (solver->profile.detailed) {
        solver->profile.codekd_wall_seconds +=
            monotonic_seconds() - search_wall_start;
    }

    if (!*presult) {
        solver->profile.search_failures++;
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        return;
    }
    if (solver->profile.detailed) {
        uint64_t kd_result_digest =
            solver_kd_result_order_digest(*presult);

        solver->profile.kd_result_order_hash =
            solver_order_hash_mix(
                solver->profile.kd_result_order_hash,
                kd_result_digest);
    }
    solver->profile.codekd_hits +=
        (unsigned long long)(*presult)->nres;

    if (solver_poll_worker_stop(solver)) {
        return;
    }

    if ((*presult)->nres) {
        double pixvals[DQMAX*2];
        double resolve_wall_start = 0.0;
        int j;
        for (j=0; j<dimquad; j++) {
            setx(pixvals, j, field_getx(solver, stars[j]));
            sety(pixvals, j, field_gety(solver, stars[j]));
        }

        if (solver->profile.detailed)
            resolve_wall_start = monotonic_seconds();

        resolve_matches(*presult, pixvals, stars, dimquad,
                        solver->numtries, solver,
                        current_parity);
        solver->profile.resolve_calls++;

        if (solver->profile.detailed) {
            solver->profile.resolve_wall_seconds +=
                monotonic_seconds() - resolve_wall_start;
        }
    }

    solver->profile.hypotheses_reduced++;
}

static void solver_index_payload_failure(
    solver_t* solver,
    const char* component) {
    logerr("[solver-io] failed to decode %s payload\n",
           component);
    solver->profile.execution_failed = TRUE;
    solver->quit_now = TRUE;
}

static void resolve_matches(kdtree_qres_t* krez, const double *field_xy,
                            const int* fieldstars, int dimquads,
                            int quads_tried,
                            solver_t* solver, anbool current_parity) {
    // "field_xy" contains the xy pixel coordinates of stars A,B,C,D forming the quad
    //    [x_A,y_A, x_B,y_B, x_C,y_C, ...]
    int jj, thisquadno;
    MatchObj mo;
    solver_candidate_payload_workspace_t* payload_workspace = NULL;

    assert(krez);
    assert(dimquads > 0);
    assert(dimquads <= DQMAX);

    if (krez->nres > 0) {
        (void)quadfile_advise_rows(
            solver->index->quads,
            krez->inds,
            krez->nres);
        payload_workspace =
            solver_candidate_payload_workspace_get(
                (size_t)krez->nres);
        if (payload_workspace) {
            size_t advice_count = 0U;

            for (jj = 0; jj < krez->nres; jj++) {
                unsigned int* prepared_stars =
                    payload_workspace->stars +
                    (size_t)jj * (size_t)DQMAX;

                payload_workspace->valid[jj] =
                    quadfile_get_stars(
                        solver->index->quads,
                        krez->inds[jj],
                        prepared_stars) == 0;
                if (!payload_workspace->valid[jj]) {
                    continue;
                }
                memcpy(
                    payload_workspace->advice_starids + advice_count,
                    prepared_stars,
                    (size_t)dimquads * sizeof(*prepared_stars));
                advice_count += (size_t)dimquads;
            }
            if (advice_count && advice_count <= (size_t)INT_MAX) {
                (void)startree_advise_rows(
                    solver->index->starkd,
                    payload_workspace->advice_starids,
                    (int)advice_count);
            }
        }
    }

    if (krez->nres && solver->ab_executor) {
        int starkd_status =
            solver_ab_ensure_starkd_ready(
                solver->ab_executor,
                0,
                0U);

        if (starkd_status > 0) {
            solver->quit_now = TRUE;
            return;
        }
        if (starkd_status < 0) {
            logerr("[solver-ab] failed to initialize StarKD lookup state\n");
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            return;
        }
    }

    if (solver_ab_try_verification_wave(
            krez,
            field_xy,
            fieldstars,
            dimquads,
            quads_tried,
            solver,
            current_parity,
            payload_workspace)) {
        return;
    }

    for (jj = 0; jj < krez->nres; jj++) {
        unsigned int star_storage[DQMAX];
        const unsigned int* star;
        double starxyz[DQMAX*3];
        double scale;
        double arcsecperpix;
        tan_t wcs;
        int i;
        anbool outofbounds = FALSE;
        double abscale;

        if (solver_poll_worker_stop(solver)) {
            return;
        }

        solver->nummatches++;
        thisquadno = krez->inds[jj];

        if (payload_workspace) {
            star = payload_workspace->stars +
                (size_t)jj * (size_t)DQMAX;
        } else {
            star = star_storage;
        }
        if ((payload_workspace &&
             !payload_workspace->valid[jj]) ||
            (!payload_workspace &&
             quadfile_get_stars(
                 solver->index->quads,
                 thisquadno,
                 star_storage))) {
            solver_index_payload_failure(
                solver, "QuadFile");
            return;
        }


        if (solver->use_radec) {
            /*
             * Preserve the original all-star loading and rejection order
             * when the sky-position constraint depends on every quad star.
             */
            for (i = 0; i < dimquads; i++) {
                if (startree_get(
                        solver->index->starkd,
                        star[i],
                        starxyz + 3 * i)) {
                    solver_index_payload_failure(
                        solver, "StarKD");
                    return;
                }

                if (distsq(starxyz + 3 * i,
                           solver->centerxyz,
                           3) > solver->r2) {
                    outofbounds = TRUE;
                    break;
                }
            }
            if (outofbounds) {
                debug("Quad match is out of bounds.\n");
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_RADEC_SKIP,
                    thisquadno,
                    (float)krez->sdists[jj]);
                solver->num_radec_skipped++;
                continue;
            }
        } else {
            /*
             * The quick scale gate uses only A and B. Delay C/D/E decoding
             * until the candidate has passed that gate.
             */
            if (startree_get(
                    solver->index->starkd,
                    star[0],
                    starxyz) ||
                startree_get(
                    solver->index->starkd,
                    star[1],
                    starxyz + 3)) {
                solver_index_payload_failure(
                    solver, "StarKD");
                return;
            }
        }

        debug("        stars [");
        for (i=0; i<dimquads; i++)
            debug("%s%i", (i?" ":""), star[i]);
        debug("]\n");

        // Quick-n-dirty scale estimate based on two stars.
        // in (rad per pix)**2
        abscale = square(distsq2rad(distsq(starxyz, starxyz+3, 3))) /
            distsq(field_xy, field_xy+2, 2);
        if (abscale > solver->abscale_high ||
            abscale < solver->abscale_low) {
            solver_record_candidate_order(
                solver,
                SOLVER_AB_CANDIDATE_ABSCALE_SKIP,
                thisquadno,
                (float)krez->sdists[jj]);
            solver->num_abscale_skipped++;
            continue;
        }

        if (!solver->use_radec) {
            for (i = 2; i < dimquads; i++) {
                if (startree_get(
                        solver->index->starkd,
                        star[i],
                        starxyz + 3 * i)) {
                    solver_index_payload_failure(
                        solver, "StarKD");
                    return;
                }
            }
        }

        // compute TAN projection from the matching quad alone.
        if (fit_tan_wcs(starxyz, field_xy, dimquads, &wcs, &scale)) {
            // bad quad.
            solver_record_candidate_order(
                solver,
                SOLVER_AB_CANDIDATE_BAD_QUAD,
                thisquadno,
                (float)krez->sdists[jj]);
            logverb("bad quad at %s:%i\n", __FILE__, __LINE__);
            continue;
        }
        arcsecperpix = scale * 3600.0;

        // FIXME - should there be scale fudge here?
        if (arcsecperpix > solver->funits_upper ||
            arcsecperpix < solver->funits_lower) {
            solver_record_candidate_order(
                solver,
                SOLVER_AB_CANDIDATE_SCALE_SKIP,
                thisquadno,
                (float)krez->sdists[jj]);
            debug("          bad scale (%g arcsec/pix, range %g %g)\n",
                  arcsecperpix, solver->funits_lower, solver->funits_upper);
            continue;
        }
        solver_record_candidate_order(
            solver,
            SOLVER_AB_CANDIDATE_VERIFY,
            thisquadno,
            (float)krez->sdists[jj]);
        solver->numscaleok++;

        set_matchobj_template(solver, &mo);
        memcpy(&(mo.wcstan), &wcs, sizeof(tan_t));
        mo.wcs_valid = TRUE;
        mo.code_err = krez->sdists[jj];
        mo.scale = arcsecperpix;
        mo.parity = current_parity;
        mo.quads_tried = quads_tried;
        mo.quads_matched = solver->nummatches;
        mo.quads_scaleok = solver->numscaleok;
        mo.quad_npeers = krez->nres;
        mo.timeused = solver->timeused;
        mo.quadno = thisquadno;
        mo.dimquads = dimquads;
        for (i=0; i<dimquads; i++) {
            mo.star[i] = star[i];
            mo.field[i] = fieldstars[i];
            mo.ids[i] = 0;
        }

        memcpy(mo.quadpix, field_xy, 2 * dimquads * sizeof(double));
        memcpy(mo.quadxyz, starxyz, 3 * dimquads * sizeof(double));

        set_center_and_radius(solver, &mo, &(mo.wcstan), NULL);

        if (solver_handle_hit(solver, &mo, NULL, FALSE))
            solver->quit_now = TRUE;

        if (unlikely(solver->quit_now))
            return;
    }
}

void solver_inject_match(solver_t* solver, MatchObj* mo, sip_t* sip) {
    solver_handle_hit(solver, mo, sip, TRUE);
}

static double solver_prepare_hit_for_verify(
    solver_t* sp,
    MatchObj* mo,
    double* logaccept) {
    assert(sp);
    assert(mo);
    assert(logaccept);

    mo->indexid = sp->index->indexid;
    mo->healpix = sp->index->healpix;
    mo->hpnside = sp->index->hpnside;
    mo->wcstan.imagew = sp->field_maxx;
    mo->wcstan.imageh = sp->field_maxy;
    mo->dimquads = quadfile_dimquads(sp->index->quads);
    *logaccept = MIN(
        sp->logratio_tokeep,
        sp->logratio_totune);
    return square(sp->verify_pix) +
        square(sp->index->index_jitter / mo->scale);
}

static int solver_handle_hit(solver_t* sp, MatchObj* mo, sip_t* verifysip,
                             anbool fake_match) {
    double match_distance_in_pixels2;
    double verify_wall_start = 0.0;
    double logaccept;

    assert(sp);
    assert(mo);

    if (solver_poll_worker_stop(sp)) {
        return FALSE;
    }

    match_distance_in_pixels2 =
        solver_prepare_hit_for_verify(
            sp, mo, &logaccept);

    if (sp->profile.detailed) {
        verify_wall_start = monotonic_seconds();
    }

    solver_payload_advise_verification_search(sp, mo);
    verify_hit(sp->index->starkd, sp->index->cutnside,
               mo, verifysip, sp->vf, match_distance_in_pixels2,
               sp->distractor_ratio, sp->field_maxx, sp->field_maxy,
               sp->logratio_bail_threshold, logaccept,
               sp->logratio_stoplooking,
               sp->distance_from_quad_bonus, fake_match);
    sp->profile.verify_calls++;

    if (sp->profile.detailed) {
        sp->profile.verify_wall_seconds +=
            monotonic_seconds() - verify_wall_start;
    }

    return solver_handle_hit_after_verify(
        sp,
        mo,
        verifysip,
        fake_match,
        match_distance_in_pixels2);
}

static int solver_handle_hit_after_verify(
    solver_t* sp,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match,
    double match_distance_in_pixels2) {
    double verify_wall_start = 0.0;
    anbool solved;

    mo->nverified = sp->num_verified++;

    if (mo->logodds >= sp->best_logodds) {
        sp->best_logodds = mo->logodds;
        logverb("Got a new best match: logodds %g.\n", mo->logodds);
    }

    if (mo->logodds >= sp->logratio_totune &&
        mo->logodds < sp->logratio_tokeep) {
        if (solver_poll_worker_stop(sp)) {
            solver_discard_verified_match(mo);
            return FALSE;
        }

        logverb("Trying to tune up this solution (logodds = %g; %g)...\n",
                mo->logodds, exp(mo->logodds));
        solver_tweak2(sp, mo, 1, NULL);
        logverb("After tuning, logodds = %g (%g)\n",
                mo->logodds, exp(mo->logodds));

        // Since we tuned up this solution, we can't just accept the
        // resulting log-odds at face value.
        if (!fake_match) {
            if (solver_poll_worker_stop(sp)) {
                solver_discard_verified_match(mo);
                return FALSE;
            }

            if (sp->profile.detailed) {
                verify_wall_start = monotonic_seconds();
            }

            solver_payload_advise_verification_search(sp, mo);
            verify_hit(sp->index->starkd, sp->index->cutnside,
                       mo, mo->sip, sp->vf, match_distance_in_pixels2,
                       sp->distractor_ratio,
                       sp->field_maxx, sp->field_maxy,
                       sp->logratio_bail_threshold,
                       sp->logratio_tokeep,
                       sp->logratio_stoplooking,
                       sp->distance_from_quad_bonus,
                       fake_match);
            sp->profile.verify_calls++;

            if (sp->profile.detailed) {
                sp->profile.verify_wall_seconds +=
                    monotonic_seconds() - verify_wall_start;
            }

            logverb("Checking tuned result: logodds = %g (%g)\n",
                    mo->logodds, exp(mo->logodds));
        }
    }

    if (mo->logodds < sp->logratio_toprint)
        return FALSE;

    // Also copy original field star coordinates
    //mo.quadpix_orig
    logverb("mo field stars:\n");
    int i;
    for (i=0; i<mo->dimquads; i++) {
        logverb("  star %i; field_xy %.1f,%.1f, field_orig %.1f,%.1f\n",
               mo->field[i], mo->quadpix[2*i+0], mo->quadpix[2*i+1],
               starxy_getx(sp->fieldxy_orig, mo->field[i]),
               starxy_gety(sp->fieldxy_orig, mo->field[i]));
        mo->quadpix_orig[2*i+0] = starxy_getx(sp->fieldxy_orig, mo->field[i]);
        mo->quadpix_orig[2*i+1] = starxy_gety(sp->fieldxy_orig, mo->field[i]);
    }

    update_timeused(sp);
    mo->timeused = sp->timeused;

    matchobj_print(mo, log_get_level());

    if (mo->logodds < sp->logratio_tokeep)
        return FALSE;

    logverb("Pixel scale: %g arcsec/pix.\n", mo->scale);
    logverb("Parity: %s.\n", (mo->parity ? "neg" : "pos"));

    mo->index = sp->index;
    mo->index_jitter = sp->index->index_jitter;

    if (sp->predistort || (sp->pixel_xscale > 0)) {
        int i;
        double* matchxy;
        double* matchxyz;
        double* weights;
        int N;
        int Ngood;
        double dx,dy;

        // Apply the distortion.
        if (sp->predistort)
            logverb("Applying the distortion pattern and recomputing WCS...\n");
        else
            logverb("Applying pixel scaling and recomputing WCS...\n");

        if (log_get_level() >= LOG_VERB) {
            printf("Initial WCS:\n");
            tan_print(&(mo->wcstan));
        }

        // this includes conflicts and distractors; we won't fill these arrays.
        N = mo->nbest;
        matchxy = malloc(N * 2 * sizeof(double));
        matchxyz = malloc(N * 3 * sizeof(double));
        weights = malloc(N * sizeof(double));

        Ngood = 0;
        for (i=0; i<N; i++) {
            if (mo->theta[i] < 0)
                continue;
            // Plug in the original (distorted) coordinates
            dx = starxy_get_x(sp->fieldxy_orig, i);
            dy = starxy_get_y(sp->fieldxy_orig, i);
            matchxy[2*Ngood + 0] = dx;
            matchxy[2*Ngood + 1] = dy;
            memcpy(matchxyz + 3*Ngood, mo->refxyz + 3*mo->theta[i],
                   3*sizeof(double));
            weights[Ngood] = verify_logodds_to_weight(mo->matchodds[i]);

            double xx,yy;
            Unused anbool ok;
            ok = tan_xyzarr2pixelxy(&mo->wcstan, matchxyz+3*Ngood, &xx, &yy);
            assert(ok);
            logverb("match: ref(%.1f, %.1f) -- undist(%.1f, %.1f) --> dist(%.1f, %.1f)\n",
                    xx, yy, starxy_get_x(sp->fieldxy, i), starxy_get_y(sp->fieldxy, i), dx, dy);
            Ngood++;
        }

        if (sp->do_tweak) {
            // Compute the SIP solution using the correspondences
            // found during verify(), but with the original (distorted) positions.
            sip_t* sip = sip_create();
            memset(sip, 0, sizeof(sip_t));
            memcpy(&(sip->wcstan), &(mo->wcstan), sizeof(tan_t));
            sip->a_order = sip->b_order = sp->tweak_aborder;
            sip->ap_order = sip->bp_order = sp->tweak_abporder;
            sip->wcstan.imagew = solver_field_width(sp);
            sip->wcstan.imageh = solver_field_height(sp);
            if (sp->set_crpix) {
                sip->wcstan.crpix[0] = sp->crpix[0];
                sip->wcstan.crpix[1] = sp->crpix[1];
                if (sp->predistort) {
                    // find matching crval...
                    sip_pixel_undistortion(sp->predistort,
                                           sp->crpix[0], sp->crpix[1], &dx, &dy);
                } else {
                    dx = sp->crpix[0] / sp->pixel_xscale;
                    dy = sp->crpix[1];
                }
                tan_pixelxy2radecarr(&mo->wcstan, dx, dy, sip->wcstan.crval);

            } else {
                // keep TAN WCS's crval but distort the crpix.
                if (sp->predistort) {
                    sip_pixel_distortion(sp->predistort,
                                         mo->wcstan.crpix[0], mo->wcstan.crpix[1],
                                         sip->wcstan.crpix+0, sip->wcstan.crpix+1);
                } else {
                    sip->wcstan.crpix[0] = mo->wcstan.crpix[0] / sp->pixel_xscale;
                    sip->wcstan.crpix[1] = mo->wcstan.crpix[1];
                }
            }

            if (log_get_level() >= LOG_VERB) {
                printf("Initial SIP on distorted positions:\n");
                sip_print(sip);
            }

            int doshift = 1;
            fit_sip_wcs(matchxyz, matchxy, weights, Ngood, &(sip->wcstan),
                        sp->tweak_aborder, sp->tweak_abporder, doshift,
                        sip);

            if (log_get_level() >= LOG_VERB) {
                printf("Final SIP on distorted positions:\n");
                sip_print(sip);
            }

            for (i=0; i<Ngood; i++) {
                double xx,yy;
                Unused anbool ok;
                ok = sip_xyzarr2pixelxy(sip, matchxyz+3*i, &xx, &yy);
                assert(ok);
                logverb("match: ref(%.1f, %.1f) -- dist(%.1f, %.1f)\n",
                        xx, yy, matchxy[2*i+0], matchxy[2*i+1]);
            }
            mo->sip = sip;

        } else {
            // Take the coordinates after applying the --predistort,
            // fit a TAN to those, and include the --predistort in the output
            // SIP WCS.
            if (sp->predistort) {
                Ngood = 0;
                for (i=0; i<N; i++) {
                    if (mo->theta[i] < 0)
                        continue;
                    // Plug in the (undistorted) coordinates
                    dx = starxy_get_x(sp->fieldxy, i);
                    dy = starxy_get_y(sp->fieldxy, i);
                    matchxy[2*Ngood + 0] = dx;
                    matchxy[2*Ngood + 1] = dy;
                }
            }

            // Compute new TAN WCS...?
            fit_tan_wcs_weighted(matchxyz, matchxy, weights, Ngood,
                                 &mo->wcstan, NULL);
            if (sp->set_crpix) {
                tan_t wcs2;
                fit_tan_wcs_move_tangent_point(matchxyz, matchxy, Ngood,
                                               sp->crpix, &mo->wcstan, &wcs2);
                fit_tan_wcs_move_tangent_point(matchxyz, matchxy, Ngood,
                                               sp->crpix, &wcs2, &mo->wcstan);
            }
            if (sp->predistort) {
                // Copy the distortion
                sip_t* sip = sip_create();
                memcpy(sip, sp->predistort, sizeof(sip_t));
                memcpy(&sip->wcstan, &mo->wcstan, sizeof(tan_t));
                mo->sip = sip;
            }
        }

        free(matchxy);
        free(matchxyz);
        free(weights);

    } else if (sp->do_tweak) {
        solver_tweak2(sp, mo, sp->tweak_aborder, verifysip);

    } else if (!verifysip && sp->set_crpix) {
        tan_t wcs2;
        tan_t wcs3;
        fit_tan_wcs_move_tangent_point(mo->quadxyz, mo->quadpix, mo->dimquads,
                                       sp->crpix, &(mo->wcstan), &wcs2);
        fit_tan_wcs_move_tangent_point(mo->quadxyz, mo->quadpix, mo->dimquads,
                                       sp->crpix, &wcs2, &wcs3);
        memcpy(&(mo->wcstan), &wcs3, sizeof(tan_t));
        /*
         Good test case:
         http://antwrp.gsfc.nasa.gov/apod/image/0912/Geminid2007_pacholka850wp.jpg
         solve-field --config backend.cfg Geminid2007_pacholka850wp.xy \
         --scale-low 10 --scale-units degwidth -v --no-tweak --continue --new-fits none \
         -o 4 --crpix-center --depth 40-45

         printf("Original WCS:\n");
         tan_print_to(&(mo->wcstan), stdout);
         printf("\n");
         printf("Moved WCS:\n");
         tan_print_to(&wcs2, stdout);
         printf("\n");
         printf("Moved again WCS:\n");
         tan_print_to(&wcs3, stdout);
         printf("\n");
         */
    }

    // If the user didn't supply a callback, or if the callback
    // returns TRUE, consider it solved.
    solved = (!sp->record_match_callback ||
              sp->record_match_callback(mo, sp->userdata));

    // New best match?
    if (!sp->have_best_match || (mo->logodds > sp->best_match.logodds)) {
        if (sp->have_best_match)
            verify_free_matchobj(&sp->best_match);
        memcpy(&sp->best_match, mo, sizeof(MatchObj));
        sp->have_best_match = TRUE;
        sp->best_index = sp->index;
    } else {
        verify_free_matchobj(mo);
    }

    if (solved) {
        sp->best_match_solves = TRUE;
        return TRUE;
    }
    return FALSE;
}

solver_t* solver_new() {
    solver_t* solver = calloc(1, sizeof(solver_t));
    solver_set_default_values(solver);
    return solver;
}

void solver_set_default_values(solver_t* solver) {
    memset(solver, 0, sizeof(solver_t));

    fitsbin_mmap_advice_state_init(
        &solver->index_mmap_policy,
        fitsbin_get_configured_mmap_policy());

    solver->indexes = pl_new(16);
    solver->funits_upper = LARGE_VAL;
    solver->logratio_bail_threshold = log(1e-100);
    solver->logratio_stoplooking = LARGE_VAL;
    solver->logratio_totune = LARGE_VAL;
    solver->parity = DEFAULT_PARITY;
    solver->codetol = DEFAULT_CODE_TOL;
    solver->distractor_ratio = DEFAULT_DISTRACTOR_RATIO;
    solver->verify_pix = DEFAULT_VERIFY_PIX;
    solver->verify_uniformize = TRUE;
    solver->verify_dedup = TRUE;
    solver->distance_from_quad_bonus = TRUE;
    solver->tweak_aborder = DEFAULT_TWEAK_ABORDER;
    solver->tweak_abporder = DEFAULT_TWEAK_ABPORDER;
}

void solver_clear_indexes(solver_t* solver) {
    pl_remove_all(solver->indexes);
    solver->index = NULL;
}

void solver_cleanup(solver_t* solver) {
    solver_free_field(solver);
    pl_free(solver->indexes);
    solver->indexes = NULL;
    if (solver->have_best_match) {
        verify_free_matchobj(&solver->best_match);
        solver->have_best_match = FALSE;
    }
    if (solver->predistort)
        sip_free(solver->predistort);
    solver->predistort = NULL;
}

void solver_free(solver_t* solver) {
    if (!solver) return;
    solver_cleanup(solver);
    free(solver);
}
