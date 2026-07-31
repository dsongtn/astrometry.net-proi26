/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "index_residency.h"

#if defined(__linux__) && defined(MFD_CLOEXEC) && \
    defined(MFD_ALLOW_SEALING) && defined(F_ADD_SEALS) && \
    defined(F_SEAL_SEAL) && defined(F_SEAL_SHRINK) && \
    defined(F_SEAL_GROW) && defined(F_SEAL_WRITE)
#define INDEX_RESIDENCY_HAVE_MEMFD 1
#else
#define INDEX_RESIDENCY_HAVE_MEMFD 0
#endif

#define INDEX_RESIDENCY_COPY_BUFFER (1024U * 1024U)
#define INDEX_RESIDENCY_BACKING_PATH 64U

static uint64_t index_residency_monotonic_nanoseconds(void) {
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now)) {
        return 0U;
    }
    if ((uint64_t)now.tv_sec >
        UINT64_MAX / 1000000000ULL) {
        return UINT64_MAX;
    }
    return (uint64_t)now.tv_sec * 1000000000ULL +
        (uint64_t)now.tv_nsec;
}

typedef enum index_residency_entry_state {
    INDEX_RESIDENCY_ENTRY_LOADING = 0,
    INDEX_RESIDENCY_ENTRY_READY = 1,
    INDEX_RESIDENCY_ENTRY_FAILED = 2
} index_residency_entry_state_t;

typedef struct index_residency_entry {
    char* source_path;
    struct stat source_stat;
    size_t bytes;
    int source_fd;
    int backing_fd;
    char backing_path[INDEX_RESIDENCY_BACKING_PATH];

    index_residency_entry_state_t state;
    index_residency_priority_t priority;
    int queued;
    int copy_was_speculative;
    int demand_cancelled;
    int permanent_source_fallback;
    uint64_t queue_sequence;
    uint64_t last_used;
    size_t leases;
    size_t waiters;

    struct index_residency_entry* next;
} index_residency_entry_t;

struct index_residency {
    pthread_mutex_t mutex;
    pthread_cond_t work_cond;
    pthread_cond_t state_cond;
    int mutex_ready;
    int work_cond_ready;
    int state_cond_ready;

    pthread_t threads[2];
    unsigned int loader_lanes;
    unsigned int threads_started;
    int backend_supported;
    int stopping;
    int cancel_requested;
    int workers_joined;
    int destroy_pending;

    index_residency_entry_t* entries;
    size_t byte_budget;
    size_t resident_bytes;
    size_t peak_resident_bytes;
    size_t queue_count;
    size_t inflight;
    size_t speculative_inflight;
    size_t live_handles;
    size_t active_acquires;
    uint64_t queue_sequence;
    uint64_t access_tick;

    index_residency_stats_t counters;
};

struct index_residency_handle {
    index_residency_t* service;
    index_residency_entry_t* entry;
    int source_fallback;
};

static long index_residency_stat_mtime_nsec(
    const struct stat* source_stat) {
#if defined(__APPLE__)
    return source_stat->st_mtimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
    return source_stat->st_mtim.tv_nsec;
#else
    (void)source_stat;
    return 0L;
#endif
}

static long index_residency_stat_ctime_nsec(
    const struct stat* source_stat) {
#if defined(__APPLE__)
    return source_stat->st_ctimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
    return source_stat->st_ctim.tv_nsec;
#else
    (void)source_stat;
    return 0L;
#endif
}

static int index_residency_same_source(
    const struct stat* left,
    const struct stat* right) {
    if (!left || !right) {
        return 0;
    }
    return left->st_dev == right->st_dev &&
        left->st_ino == right->st_ino &&
        left->st_size == right->st_size &&
        left->st_mtime == right->st_mtime &&
        left->st_ctime == right->st_ctime &&
        index_residency_stat_mtime_nsec(left) ==
            index_residency_stat_mtime_nsec(right) &&
        index_residency_stat_ctime_nsec(left) ==
            index_residency_stat_ctime_nsec(right);
}

/*
 * Open through the configured path, then confirm that the pathname and open
 * description identify one coherent source epoch. A transient replacement is
 * retried once.
 */
static int index_residency_open_source(
    const char* path,
    int* source_fd,
    struct stat* source_stat,
    int* source_changed) {
    int attempt;

    if (!path || !path[0] || !source_fd || !source_stat) {
        return -1;
    }
    *source_fd = -1;
    if (source_changed) {
        *source_changed = 0;
    }
    for (attempt = 0; attempt < 2; attempt++) {
        struct stat path_stat;
        int flags = O_RDONLY;
        int fd;

#ifdef O_CLOEXEC
        flags |= O_CLOEXEC;
#endif
        fd = open(path, flags);
        if (fd < 0) {
            return -1;
        }
        if (fstat(fd, source_stat) ||
            !S_ISREG(source_stat->st_mode) ||
            source_stat->st_size < 0 ||
            stat(path, &path_stat)) {
            close(fd);
            return -1;
        }
        if (index_residency_same_source(
                source_stat, &path_stat)) {
            *source_fd = fd;
            return 0;
        }
        if (source_changed) {
            *source_changed = 1;
        }
        close(fd);
    }
    return -1;
}

static index_residency_entry_t* index_residency_entry_new(
    const char* path,
    int source_fd,
    const struct stat* source_stat,
    index_residency_priority_t priority) {
    index_residency_entry_t* entry;

    if (!path || source_fd < 0 || !source_stat ||
        source_stat->st_size < 0 ||
        (uintmax_t)source_stat->st_size >
            (uintmax_t)SIZE_MAX) {
        return NULL;
    }
    entry = calloc(1, sizeof(*entry));
    if (!entry) {
        return NULL;
    }
    entry->source_path = strdup(path);
    if (!entry->source_path) {
        free(entry);
        return NULL;
    }
    entry->source_stat = *source_stat;
    entry->bytes = (size_t)source_stat->st_size;
    entry->source_fd = source_fd;
    entry->backing_fd = -1;
    entry->state = INDEX_RESIDENCY_ENTRY_LOADING;
    entry->priority = priority;
    return entry;
}

static void index_residency_entry_free(
    index_residency_entry_t* entry) {
    if (!entry) {
        return;
    }
    if (entry->source_fd >= 0) {
        close(entry->source_fd);
    }
    if (entry->backing_fd >= 0) {
        close(entry->backing_fd);
    }
    free(entry->source_path);
    free(entry);
}

static index_residency_entry_t* index_residency_find_locked(
    index_residency_t* service,
    const struct stat* source_stat) {
    index_residency_entry_t* entry;

    for (entry = service->entries;
         entry;
         entry = entry->next) {
        if (index_residency_same_source(
                &entry->source_stat, source_stat)) {
            return entry;
        }
    }
    return NULL;
}

static void index_residency_remove_locked(
    index_residency_t* service,
    index_residency_entry_t* entry) {
    index_residency_entry_t** cursor;

    if (!service || !entry) {
        return;
    }
    cursor = &service->entries;
    while (*cursor && *cursor != entry) {
        cursor = &(*cursor)->next;
    }
    if (*cursor == entry) {
        *cursor = entry->next;
        entry->next = NULL;
    }
}

static int index_residency_can_reserve_locked(
    const index_residency_t* service,
    size_t bytes) {
    return bytes <= service->byte_budget &&
        service->resident_bytes <=
            service->byte_budget - bytes;
}

static void index_residency_queue_locked(
    index_residency_t* service,
    index_residency_entry_t* entry) {
    entry->state = INDEX_RESIDENCY_ENTRY_LOADING;
    if (entry->source_fd >= 0) {
        close(entry->source_fd);
        entry->source_fd = -1;
    }
    entry->queued = 1;
    entry->queue_sequence = ++service->queue_sequence;
    entry->next = service->entries;
    service->entries = entry;
    service->queue_count++;
    service->resident_bytes += entry->bytes;
    if (service->resident_bytes >
        service->peak_resident_bytes) {
        service->peak_resident_bytes =
            service->resident_bytes;
    }
    service->counters.files_queued++;
    pthread_cond_signal(&service->work_cond);
}

static index_residency_entry_t*
index_residency_select_work_locked(
    index_residency_t* service) {
    index_residency_entry_t* selected = NULL;
    index_residency_entry_t* entry;

    for (entry = service->entries;
         entry;
         entry = entry->next) {
        if (entry->state != INDEX_RESIDENCY_ENTRY_LOADING ||
            !entry->queued) {
            continue;
        }
        if (service->loader_lanes > 1U &&
            entry->priority ==
                INDEX_RESIDENCY_PRIORITY_SPECULATIVE &&
            service->speculative_inflight >=
                (size_t)service->loader_lanes - 1U) {
            continue;
        }
        if (!selected ||
            entry->priority > selected->priority ||
            (entry->priority == selected->priority &&
             entry->queue_sequence <
                 selected->queue_sequence)) {
            selected = entry;
        }
    }
    return selected;
}

static int index_residency_validate_after_copy(
    const index_residency_entry_t* entry,
    int source_fd) {
    struct stat open_stat;
    struct stat path_stat;

    if (!entry || source_fd < 0 ||
        fstat(source_fd, &open_stat) ||
        stat(entry->source_path, &path_stat)) {
        return -1;
    }
    if (!index_residency_same_source(
            &entry->source_stat, &open_stat) ||
        !index_residency_same_source(
            &entry->source_stat, &path_stat)) {
        errno = ESTALE;
        return 1;
    }
    return 0;
}

static int index_residency_write_all(
    int fd,
    const void* data,
    size_t size,
    off_t offset) {
    const unsigned char* cursor = data;

    while (size) {
        ssize_t written = pwrite(fd, cursor, size, offset);

        if (written < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (!written) {
            errno = EIO;
            return -1;
        }
        cursor += (size_t)written;
        size -= (size_t)written;
        offset += written;
    }
    return 0;
}

static int index_residency_read_all(
    int fd,
    void* data,
    size_t size,
    off_t offset) {
    unsigned char* cursor = data;

    while (size) {
        ssize_t count = pread(fd, cursor, size, offset);

        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (!count) {
            errno = EIO;
            return -1;
        }
        cursor += (size_t)count;
        size -= (size_t)count;
        offset += count;
    }
    return 0;
}

static int index_residency_copy_source(
    index_residency_t* service,
    index_residency_entry_t* entry,
    int source_fd,
    int* backing_fd,
    char* backing_path,
    size_t backing_path_size,
    int* source_changed) {
#if INDEX_RESIDENCY_HAVE_MEMFD
    unsigned char* buffer = NULL;
    off_t offset = 0;
    int writable_fd = -1;
    int readonly_fd = -1;
    int result = -1;

    if (!service || !entry || source_fd < 0 ||
        !backing_fd || !backing_path ||
        !backing_path_size || !source_changed) {
        errno = EINVAL;
        return -1;
    }
    *backing_fd = -1;
    *source_changed = 0;
    backing_path[0] = '\0';

    writable_fd = memfd_create(
        "astrometry-index-resident",
        MFD_CLOEXEC | MFD_ALLOW_SEALING);
    if (writable_fd < 0 ||
        ftruncate(writable_fd, entry->source_stat.st_size)) {
        goto cleanup;
    }
    buffer = malloc(INDEX_RESIDENCY_COPY_BUFFER);
    if (!buffer) {
        goto cleanup;
    }
    while (offset < entry->source_stat.st_size) {
        off_t remaining = entry->source_stat.st_size - offset;
        size_t request = remaining >
            (off_t)INDEX_RESIDENCY_COPY_BUFFER
                ? INDEX_RESIDENCY_COPY_BUFFER
                : (size_t)remaining;
        if (__atomic_load_n(
                &service->cancel_requested,
                __ATOMIC_ACQUIRE) ||
            __atomic_load_n(
                &entry->demand_cancelled,
                __ATOMIC_ACQUIRE)) {
            errno = ECANCELED;
            goto cleanup;
        }
        if (index_residency_read_all(
                source_fd,
                buffer,
                request,
                offset) ||
            index_residency_write_all(
                writable_fd,
                buffer,
                request,
                offset)) {
            goto cleanup;
        }
        offset += (off_t)request;
    }
    result = index_residency_validate_after_copy(
        entry, source_fd);
    if (result) {
        if (result > 0) {
            *source_changed = 1;
        }
        result = -1;
        goto cleanup;
    }
    if (fchmod(writable_fd, S_IRUSR) ||
        fcntl(
            writable_fd,
            F_ADD_SEALS,
            F_SEAL_SEAL |
                F_SEAL_SHRINK |
                F_SEAL_GROW |
                F_SEAL_WRITE)) {
        goto cleanup;
    }
    {
        int path_length = snprintf(
            backing_path,
            backing_path_size,
            "/proc/self/fd/%i",
            writable_fd);

        if (path_length < 0 ||
            (size_t)path_length >= backing_path_size) {
            errno = ENAMETOOLONG;
            goto cleanup;
        }
    }
    readonly_fd = open(
        backing_path,
        O_RDONLY | O_CLOEXEC);
    if (readonly_fd < 0) {
        goto cleanup;
    }
    {
        int path_length = snprintf(
            backing_path,
            backing_path_size,
            "/proc/self/fd/%i",
            readonly_fd);

        if (path_length < 0 ||
            (size_t)path_length >= backing_path_size) {
            errno = ENAMETOOLONG;
            goto cleanup;
        }
    }
    close(writable_fd);
    writable_fd = -1;
    *backing_fd = readonly_fd;
    readonly_fd = -1;
    result = 0;

cleanup:
    free(buffer);
    if (readonly_fd >= 0) {
        close(readonly_fd);
    }
    if (writable_fd >= 0) {
        close(writable_fd);
    }
    return result;
#else
    (void)entry;
    (void)backing_fd;
    (void)service;
    (void)source_fd;
    (void)backing_path;
    (void)backing_path_size;
    (void)source_changed;
    errno = ENOTSUP;
    return -1;
#endif
}

static void* index_residency_loader_main(void* opaque) {
    index_residency_t* service = opaque;

    while (1) {
        index_residency_entry_t* entry;
        index_residency_entry_t* discard = NULL;
        char backing_path[INDEX_RESIDENCY_BACKING_PATH] = {0};
        struct stat source_stat;
        int backing_fd = -1;
        int source_changed = 0;
        int source_fd = -1;
        int copy_status;

        pthread_mutex_lock(&service->mutex);
        while (!service->stopping &&
               !service->queue_count) {
            pthread_cond_wait(
                &service->work_cond,
                &service->mutex);
        }
        if (service->stopping) {
            pthread_mutex_unlock(&service->mutex);
            break;
        }
        entry = index_residency_select_work_locked(
            service);
        if (!entry) {
            pthread_cond_wait(
                &service->work_cond,
                &service->mutex);
            pthread_mutex_unlock(&service->mutex);
            continue;
        }
        entry->queued = 0;
        service->queue_count--;
        service->inflight++;
        entry->copy_was_speculative =
            entry->priority ==
                INDEX_RESIDENCY_PRIORITY_SPECULATIVE;
        if (entry->copy_was_speculative) {
            service->speculative_inflight++;
        }
        pthread_mutex_unlock(&service->mutex);

        copy_status = index_residency_open_source(
            entry->source_path,
            &source_fd,
            &source_stat,
            &source_changed);
        if (!copy_status &&
            !index_residency_same_source(
                &entry->source_stat,
                &source_stat)) {
            source_changed = 1;
            errno = ESTALE;
            copy_status = -1;
        }
        if (!copy_status) {
            copy_status = index_residency_copy_source(
                service,
                entry,
                source_fd,
                &backing_fd,
                backing_path,
                sizeof(backing_path),
                &source_changed);
        }

        pthread_mutex_lock(&service->mutex);
        service->inflight--;
        if (entry->copy_was_speculative) {
            if (service->speculative_inflight) {
                service->speculative_inflight--;
            }
            entry->copy_was_speculative = 0;
        }
        if (!copy_status && !service->stopping &&
            !__atomic_load_n(
                &entry->demand_cancelled,
                __ATOMIC_ACQUIRE)) {
            entry->backing_fd = backing_fd;
            backing_fd = -1;
            memcpy(
                entry->backing_path,
                backing_path,
                sizeof(entry->backing_path));
            entry->state = INDEX_RESIDENCY_ENTRY_READY;
            entry->last_used = ++service->access_tick;
            service->counters.files_copied++;
            service->counters.bytes_copied += entry->bytes;
        } else {
            entry->state = INDEX_RESIDENCY_ENTRY_FAILED;
            service->resident_bytes -= entry->bytes;
            if (source_changed) {
                service->counters.source_changes++;
            } else if (service->stopping ||
                       __atomic_load_n(
                           &entry->demand_cancelled,
                           __ATOMIC_ACQUIRE)) {
                service->counters.cancelled_entries++;
            } else {
                service->counters.copy_failures++;
            }
            if (!entry->waiters) {
                index_residency_remove_locked(
                    service, entry);
                discard = entry;
            }
        }
        pthread_cond_broadcast(&service->state_cond);
        pthread_cond_broadcast(&service->work_cond);
        pthread_mutex_unlock(&service->mutex);

        if (source_fd >= 0) {
            close(source_fd);
        }
        if (backing_fd >= 0) {
            close(backing_fd);
        }
        index_residency_entry_free(discard);
    }
    return NULL;
}

static int index_residency_service_can_destroy_locked(
    const index_residency_t* service) {
    return service->destroy_pending &&
        service->workers_joined &&
        !service->live_handles &&
        !service->active_acquires;
}

static void index_residency_service_destroy(
    index_residency_t* service) {
    index_residency_entry_t* entry;

    if (!service) {
        return;
    }
    entry = service->entries;
    while (entry) {
        index_residency_entry_t* next = entry->next;

        index_residency_entry_free(entry);
        entry = next;
    }
    if (service->state_cond_ready) {
        pthread_cond_destroy(&service->state_cond);
    }
    if (service->work_cond_ready) {
        pthread_cond_destroy(&service->work_cond);
    }
    if (service->mutex_ready) {
        pthread_mutex_destroy(&service->mutex);
    }
    free(service);
}

static int index_residency_priority_valid(
    index_residency_priority_t priority) {
    return priority >=
            INDEX_RESIDENCY_PRIORITY_SPECULATIVE &&
        priority <= INDEX_RESIDENCY_PRIORITY_DEMAND;
}

int index_residency_start(size_t byte_budget,
                          unsigned int loader_lanes,
                          index_residency_t** service_out) {
    index_residency_t* service;
    unsigned int i;

    if (!service_out || loader_lanes < 1U ||
        loader_lanes > 2U) {
        return -1;
    }
    *service_out = NULL;
    service = calloc(1, sizeof(*service));
    if (!service) {
        return -1;
    }
    service->byte_budget = byte_budget;
    service->loader_lanes = loader_lanes;
    service->backend_supported =
        INDEX_RESIDENCY_HAVE_MEMFD;
    service->counters.byte_budget = byte_budget;
    service->counters.loader_lanes = loader_lanes;
    service->counters.backend_supported =
        service->backend_supported;

    if (pthread_mutex_init(&service->mutex, NULL)) {
        index_residency_service_destroy(service);
        return -1;
    }
    service->mutex_ready = 1;
    if (pthread_cond_init(&service->work_cond, NULL)) {
        index_residency_service_destroy(service);
        return -1;
    }
    service->work_cond_ready = 1;
    if (pthread_cond_init(&service->state_cond, NULL)) {
        index_residency_service_destroy(service);
        return -1;
    }
    service->state_cond_ready = 1;

    if (!service->backend_supported) {
        service->workers_joined = 1;
        *service_out = service;
        return 0;
    }
    for (i = 0U; i < loader_lanes; i++) {
        if (pthread_create(
                &service->threads[i],
                NULL,
                index_residency_loader_main,
                service)) {
            unsigned int j;

            pthread_mutex_lock(&service->mutex);
            service->stopping = 1;
            pthread_cond_broadcast(&service->work_cond);
            pthread_mutex_unlock(&service->mutex);
            for (j = 0U;
                 j < service->threads_started;
                 j++) {
                pthread_join(service->threads[j], NULL);
            }
            service->workers_joined = 1;
            index_residency_service_destroy(service);
            return -1;
        }
        service->threads_started++;
    }
    *service_out = service;
    return 0;
}

index_residency_result_t index_residency_prepare(
    index_residency_t* service,
    const char* path,
    index_residency_priority_t priority) {
    index_residency_entry_t* candidate;
    index_residency_entry_t* entry;
    struct stat source_stat;
    int source_changed = 0;
    int source_fd = -1;

    if (!service || !path || !path[0] ||
        !index_residency_priority_valid(priority)) {
        return INDEX_RESIDENCY_ERROR;
    }
    pthread_mutex_lock(&service->mutex);
    service->counters.prepare_requests++;
    if (service->stopping ||
        !service->backend_supported) {
        pthread_mutex_unlock(&service->mutex);
        return INDEX_RESIDENCY_FALLBACK;
    }
    pthread_mutex_unlock(&service->mutex);

    if (index_residency_open_source(
            path,
            &source_fd,
            &source_stat,
            &source_changed)) {
        if (source_changed) {
            pthread_mutex_lock(&service->mutex);
            service->counters.source_changes++;
            pthread_mutex_unlock(&service->mutex);
        }
        return INDEX_RESIDENCY_FALLBACK;
    }
    candidate = index_residency_entry_new(
        path, source_fd, &source_stat, priority);
    if (!candidate) {
        close(source_fd);
        return INDEX_RESIDENCY_FALLBACK;
    }

    pthread_mutex_lock(&service->mutex);
    if (service->stopping) {
        pthread_mutex_unlock(&service->mutex);
        index_residency_entry_free(candidate);
        return INDEX_RESIDENCY_FALLBACK;
    }
    entry = index_residency_find_locked(
        service, &source_stat);
    if (entry &&
        entry->state == INDEX_RESIDENCY_ENTRY_FAILED &&
        !entry->demand_cancelled &&
        !entry->waiters && !entry->leases) {
        index_residency_remove_locked(service, entry);
        index_residency_entry_free(entry);
        entry = NULL;
    }
    if (entry &&
        entry->state == INDEX_RESIDENCY_ENTRY_READY) {
        entry->last_used = ++service->access_tick;
        pthread_mutex_unlock(&service->mutex);
        index_residency_entry_free(candidate);
        return INDEX_RESIDENCY_ACCEPTED;
    }
    if (entry &&
        entry->state == INDEX_RESIDENCY_ENTRY_LOADING) {
        if (priority > entry->priority) {
            entry->priority = priority;
            pthread_cond_signal(&service->work_cond);
        }
        service->counters.loading_deduplications++;
        pthread_mutex_unlock(&service->mutex);
        index_residency_entry_free(candidate);
        return INDEX_RESIDENCY_ACCEPTED;
    }
    if (entry ||
        !index_residency_can_reserve_locked(
            service, candidate->bytes)) {
        service->counters.budget_refusals++;
        pthread_mutex_unlock(&service->mutex);
        index_residency_entry_free(candidate);
        return INDEX_RESIDENCY_FALLBACK;
    }
    index_residency_queue_locked(service, candidate);
    pthread_mutex_unlock(&service->mutex);
    return INDEX_RESIDENCY_ACCEPTED;
}

static void index_residency_finish_acquire(
    index_residency_t* service) {
    int destroy;

    pthread_mutex_lock(&service->mutex);
    if (service->active_acquires) {
        service->active_acquires--;
    }
    destroy =
        index_residency_service_can_destroy_locked(
            service);
    pthread_mutex_unlock(&service->mutex);
    if (destroy) {
        index_residency_service_destroy(service);
    }
}

index_residency_result_t index_residency_acquire(
    index_residency_t* service,
    const char* path,
    index_residency_handle_t** handle_out) {
    index_residency_handle_t* handle = NULL;
    index_residency_entry_t* discard = NULL;
    index_residency_entry_t* entry;
    struct stat source_stat;
    int destroy;

    if (!service || !path || !path[0] || !handle_out) {
        return INDEX_RESIDENCY_ERROR;
    }
    *handle_out = NULL;
    pthread_mutex_lock(&service->mutex);
    service->counters.acquire_requests++;
    if (service->stopping ||
        !service->backend_supported) {
        pthread_mutex_unlock(&service->mutex);
        return INDEX_RESIDENCY_FALLBACK;
    }
    service->active_acquires++;
    pthread_mutex_unlock(&service->mutex);

    if (stat(path, &source_stat) ||
        !S_ISREG(source_stat.st_mode)) {
        index_residency_finish_acquire(service);
        return INDEX_RESIDENCY_FALLBACK;
    }

    pthread_mutex_lock(&service->mutex);
    entry = index_residency_find_locked(
        service, &source_stat);
    if (entry &&
        entry->state == INDEX_RESIDENCY_ENTRY_FAILED &&
        !entry->demand_cancelled &&
        !entry->waiters && !entry->leases) {
        index_residency_remove_locked(service, entry);
        index_residency_entry_free(entry);
        entry = NULL;
    }
    if (!entry) {
        if (service->active_acquires) {
            service->active_acquires--;
        }
        destroy =
            index_residency_service_can_destroy_locked(service);
        pthread_mutex_unlock(&service->mutex);
        free(handle);
        if (destroy) {
            index_residency_service_destroy(service);
        }
        return INDEX_RESIDENCY_FALLBACK;
    }
    if (entry->state == INDEX_RESIDENCY_ENTRY_LOADING) {
        uint64_t wait_start = 0U;

        service->counters.loading_deduplications++;
        __atomic_store_n(
            &entry->demand_cancelled,
            1, __ATOMIC_RELEASE);
        if (entry->queued) {
            entry->queued = 0;
            entry->state =
                INDEX_RESIDENCY_ENTRY_FAILED;
            service->queue_count--;
            service->resident_bytes -= entry->bytes;
            service->counters.cancelled_entries++;
        } else {
            wait_start =
                index_residency_monotonic_nanoseconds();
            entry->waiters++;
            service->counters.wait_count++;
            pthread_cond_broadcast(
                &service->work_cond);
            while (entry->state ==
                   INDEX_RESIDENCY_ENTRY_LOADING) {
                (void)pthread_cond_wait(
                    &service->state_cond,
                    &service->mutex);
            }
            if (entry->waiters) {
                entry->waiters--;
            }
            if (wait_start) {
                uint64_t wait_end =
                    index_residency_monotonic_nanoseconds();

                if (wait_end >= wait_start) {
                    service->counters.wait_nanoseconds +=
                        wait_end - wait_start;
                }
            }
        }
        pthread_cond_broadcast(&service->work_cond);
        pthread_cond_broadcast(&service->state_cond);
    }
    handle = calloc(1, sizeof(*handle));
    if (!handle) {
        if (!service->stopping &&
            entry->state ==
                INDEX_RESIDENCY_ENTRY_FAILED &&
            __atomic_load_n(
                &entry->demand_cancelled,
                __ATOMIC_ACQUIRE)) {
            entry->permanent_source_fallback = 1;
        }
        if (service->active_acquires) {
            service->active_acquires--;
        }
        destroy =
            index_residency_service_can_destroy_locked(
                service);
        pthread_mutex_unlock(&service->mutex);
        if (destroy) {
            index_residency_service_destroy(service);
        }
        return INDEX_RESIDENCY_FALLBACK;
    }
    if (entry->state == INDEX_RESIDENCY_ENTRY_READY &&
        !service->stopping) {
        entry->leases++;
        entry->last_used = ++service->access_tick;
        service->live_handles++;
        service->counters.cache_hits++;
        handle->service = service;
        handle->entry = entry;
        *handle_out = handle;
        if (service->active_acquires) {
            service->active_acquires--;
        }
        pthread_mutex_unlock(&service->mutex);
        return INDEX_RESIDENCY_ACCEPTED;
    }
    if (entry->state == INDEX_RESIDENCY_ENTRY_FAILED &&
        __atomic_load_n(
            &entry->demand_cancelled,
            __ATOMIC_ACQUIRE) &&
        !service->stopping) {
        entry->leases++;
        service->live_handles++;
        service->counters.source_leases++;
        handle->service = service;
        handle->entry = entry;
        handle->source_fallback = 1;
        *handle_out = handle;
        if (service->active_acquires) {
            service->active_acquires--;
        }
        pthread_mutex_unlock(&service->mutex);
        return INDEX_RESIDENCY_SOURCE_LEASE;
    }
    if (entry->state == INDEX_RESIDENCY_ENTRY_FAILED &&
        !entry->waiters && !entry->leases) {
        index_residency_remove_locked(service, entry);
        discard = entry;
    }
    if (service->active_acquires) {
        service->active_acquires--;
    }
    destroy =
        index_residency_service_can_destroy_locked(service);
    pthread_mutex_unlock(&service->mutex);
    free(handle);
    index_residency_entry_free(discard);
    if (destroy) {
        index_residency_service_destroy(service);
    }
    return INDEX_RESIDENCY_FALLBACK;
}

const char* index_residency_handle_path(
    const index_residency_handle_t* handle) {
    if (!handle || !handle->entry ||
        handle->entry->state !=
            INDEX_RESIDENCY_ENTRY_READY) {
        return NULL;
    }
    return handle->entry->backing_path;
}

const struct stat* index_residency_handle_source_stat(
    const index_residency_handle_t* handle) {
    if (!handle || !handle->entry ||
        handle->entry->state !=
            INDEX_RESIDENCY_ENTRY_READY) {
        return NULL;
    }
    return &handle->entry->source_stat;
}

static int index_residency_requeue_source_locked(
    index_residency_t* service,
    index_residency_entry_t* entry) {
    if (!service || !entry || service->stopping ||
        entry->state != INDEX_RESIDENCY_ENTRY_FAILED ||
        !__atomic_load_n(
            &entry->demand_cancelled,
            __ATOMIC_ACQUIRE)) {
        return -1;
    }
    if (!index_residency_can_reserve_locked(
            service, entry->bytes)) {
        return -1;
    }
    __atomic_store_n(
        &entry->demand_cancelled,
        0, __ATOMIC_RELEASE);
    entry->state = INDEX_RESIDENCY_ENTRY_LOADING;
    entry->priority = INDEX_RESIDENCY_PRIORITY_LOOKAHEAD;
    entry->queued = 1;
    entry->queue_sequence = ++service->queue_sequence;
    service->queue_count++;
    service->resident_bytes += entry->bytes;
    if (service->resident_bytes >
        service->peak_resident_bytes) {
        service->peak_resident_bytes =
            service->resident_bytes;
    }
    service->counters.files_queued++;
    service->counters.source_requeues++;
    pthread_cond_signal(&service->work_cond);
    return 0;
}

void index_residency_release(
    index_residency_handle_t* handle) {
    index_residency_entry_t* discard = NULL;
    index_residency_t* service;
    int destroy;

    if (!handle || !handle->service ||
        !handle->entry) {
        free(handle);
        return;
    }
    service = handle->service;
    pthread_mutex_lock(&service->mutex);
    if (handle->entry->leases) {
        handle->entry->leases--;
    }
    if (service->live_handles) {
        service->live_handles--;
    }
    if (handle->source_fallback &&
        !handle->entry->leases &&
        !handle->entry->waiters &&
        !service->stopping &&
        !handle->entry->permanent_source_fallback &&
        handle->entry->state ==
            INDEX_RESIDENCY_ENTRY_FAILED) {
        if (index_residency_requeue_source_locked(
                service, handle->entry)) {
            index_residency_remove_locked(
                service, handle->entry);
            discard = handle->entry;
        }
    }
    if (service->stopping &&
        !handle->entry->leases &&
        !handle->entry->waiters &&
        handle->entry->state !=
            INDEX_RESIDENCY_ENTRY_LOADING) {
        index_residency_remove_locked(
            service, handle->entry);
        if (handle->entry->state ==
            INDEX_RESIDENCY_ENTRY_READY) {
            service->resident_bytes -=
                handle->entry->bytes;
        }
        discard = handle->entry;
    }
    destroy =
        index_residency_service_can_destroy_locked(
            service);
    pthread_mutex_unlock(&service->mutex);

    free(handle);
    index_residency_entry_free(discard);
    if (destroy) {
        index_residency_service_destroy(service);
    }
}

int index_residency_drain(
    index_residency_t* service) {
    int status = 0;

    if (!service) {
        return -1;
    }
    pthread_mutex_lock(&service->mutex);
    while (!service->stopping &&
           (service->queue_count ||
            service->inflight) &&
           !status) {
        status = pthread_cond_wait(
            &service->state_cond,
            &service->mutex);
    }
    if (service->stopping) {
        status = -1;
    }
    pthread_mutex_unlock(&service->mutex);
    return status ? -1 : 0;
}

int index_residency_get_stats(
    index_residency_t* service,
    index_residency_stats_t* stats) {
    index_residency_entry_t* entry;

    if (!service || !stats) {
        return -1;
    }
    pthread_mutex_lock(&service->mutex);
    *stats = service->counters;
    stats->byte_budget = service->byte_budget;
    stats->resident_bytes =
        service->resident_bytes;
    stats->peak_resident_bytes =
        service->peak_resident_bytes;
    stats->live_handles =
        service->live_handles;
    stats->loader_lanes =
        service->loader_lanes;
    stats->backend_supported =
        service->backend_supported;
    for (entry = service->entries;
         entry;
         entry = entry->next) {
        switch (entry->state) {
        case INDEX_RESIDENCY_ENTRY_LOADING:
            stats->loading_entries++;
            stats->loading_bytes += entry->bytes;
            break;
        case INDEX_RESIDENCY_ENTRY_READY:
            stats->ready_entries++;
            stats->ready_bytes += entry->bytes;
            if (entry->leases) {
                stats->leased_bytes += entry->bytes;
            }
            break;
        case INDEX_RESIDENCY_ENTRY_FAILED:
            stats->failed_entries++;
            break;
        }
    }
    pthread_mutex_unlock(&service->mutex);
    return 0;
}

int index_residency_quiesce(
    index_residency_t* service) {
    index_residency_entry_t* entry;

    if (!service) {
        return 0;
    }
    pthread_mutex_lock(&service->mutex);
    if (service->stopping) {
        pthread_mutex_unlock(&service->mutex);
        return 0;
    }
    service->stopping = 1;
    __atomic_store_n(
        &service->cancel_requested,
        1, __ATOMIC_RELEASE);
    for (entry = service->entries;
         entry;
         entry = entry->next) {
        if (entry->state !=
                INDEX_RESIDENCY_ENTRY_LOADING ||
            !entry->queued) {
            continue;
        }
        entry->queued = 0;
        entry->state =
            INDEX_RESIDENCY_ENTRY_FAILED;
        service->queue_count--;
        service->resident_bytes -= entry->bytes;
        service->counters.cancelled_entries++;
        if (entry->source_fd >= 0) {
            close(entry->source_fd);
            entry->source_fd = -1;
        }
    }
    pthread_cond_broadcast(&service->work_cond);
    pthread_cond_broadcast(&service->state_cond);
    pthread_mutex_unlock(&service->mutex);
    return 0;
}

int index_residency_stop(
    index_residency_t* service) {
    unsigned int i;
    int join_error = 0;
    int destroy;

    if (!service) {
        return 0;
    }
    (void)index_residency_quiesce(service);

    for (i = 0U; i < service->threads_started; i++) {
        if (pthread_join(service->threads[i], NULL)) {
            join_error = -1;
        }
    }

    pthread_mutex_lock(&service->mutex);
    service->workers_joined = !join_error;
    service->destroy_pending = 1;
    destroy =
        index_residency_service_can_destroy_locked(
            service);
    pthread_mutex_unlock(&service->mutex);
    if (destroy) {
        index_residency_service_destroy(service);
    }
    return join_error;
}
