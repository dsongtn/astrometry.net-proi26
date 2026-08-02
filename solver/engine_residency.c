/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "engine_private.h"
#include "log.h"
#include "mathutil.h"

static anbool engine_index_residency_eligible(
    const index_t* index) {
    return index && index->indexfn &&
        !index->codekd && !index->quads && !index->starkd;
}

static int engine_index_cohort_measure(
    const engine_t* engine,
    size_t* cohort_bytes,
    size_t* cohort_files) {
    struct stat* sources;
    size_t source_count = 0U;
    size_t bytes = 0U;
    int index_count;
    int i;

    if (!engine || !engine->indexes ||
        !cohort_bytes || !cohort_files) {
        return -1;
    }
    index_count = pl_size(engine->indexes);
    sources = calloc(
        index_count ? (size_t)index_count : 1U,
        sizeof(*sources));
    if (!sources) {
        return -1;
    }
    for (i = 0; i < index_count; i++) {
        const index_t* index = pl_get(engine->indexes, i);
        struct stat source;
        size_t j;
        anbool duplicate = FALSE;

        if (!engine_index_residency_eligible(index) ||
            stat(index->indexfn, &source) ||
            !S_ISREG(source.st_mode) ||
            source.st_size < 0 ||
            (uintmax_t)source.st_size > (uintmax_t)SIZE_MAX) {
            free(sources);
            return -1;
        }
        for (j = 0U; j < source_count; j++) {
            if (sources[j].st_dev == source.st_dev &&
                sources[j].st_ino == source.st_ino) {
                duplicate = TRUE;
                break;
            }
        }
        if (duplicate) {
            continue;
        }
        if ((size_t)source.st_size > SIZE_MAX - bytes) {
            free(sources);
            return -1;
        }
        sources[source_count++] = source;
        bytes += (size_t)source.st_size;
    }
    free(sources);
    *cohort_bytes = bytes;
    *cohort_files = source_count;
    return 0;
}

static int engine_available_memory(size_t* available_bytes) {
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE)
    long pages;
    long page_size;

    if (!available_bytes) {
        return -1;
    }
    pages = sysconf(_SC_AVPHYS_PAGES);
    page_size = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || page_size <= 0 ||
        (uintmax_t)pages >
            (uintmax_t)SIZE_MAX / (uintmax_t)page_size) {
        return -1;
    }
    *available_bytes = (size_t)pages * (size_t)page_size;
    return 0;
#else
    (void)available_bytes;
    return -1;
#endif
}

static int engine_read_memory_limit(
    const char* path,
    size_t* value) {
    char buffer[64];
    char* end;
    char* token;
    FILE* file;
    uintmax_t parsed;

    if (!path || !value) {
        errno = EINVAL;
        return -1;
    }
    file = fopen(path, "r");
    if (!file) {
        return -1;
    }
    if (!fgets(buffer, sizeof(buffer), file)) {
        int read_error = errno;

        fclose(file);
        errno = read_error ? read_error : EIO;
        return -1;
    }
    fclose(file);
    token = buffer;
    while (*token == ' ' || *token == '\t' ||
           *token == '\r' || *token == '\n') {
        token++;
    }
    end = token + strlen(token);
    while (end > token &&
           (end[-1] == ' ' || end[-1] == '\t' ||
            end[-1] == '\r' || end[-1] == '\n')) {
        end--;
    }
    *end = '\0';
    if (!strcmp(token, "max")) {
        return 1;
    }
    if (*token < '0' || *token > '9') {
        errno = EINVAL;
        return -1;
    }
    errno = 0;
    parsed = strtoumax(token, &end, 10);
    if (errno || end == token || *end ||
        parsed > SIZE_MAX) {
        errno = EINVAL;
        return -1;
    }
    *value = (size_t)parsed;
    return 0;
}

#ifdef __linux__
#define ENGINE_CGROUP_PATH_SIZE 4096U
#define ENGINE_CGROUP_LINE_SIZE (ENGINE_CGROUP_PATH_SIZE * 4U)

static anbool engine_cgroup_list_contains(
    const char* list,
    const char* item) {
    size_t item_length;

    if (!list || !item) {
        return FALSE;
    }
    item_length = strlen(item);
    while (*list) {
        const char* end = strchr(list, ',');
        size_t length = end
            ? (size_t)(end - list) : strlen(list);

        if (length == item_length &&
            !strncmp(list, item, length)) {
            return TRUE;
        }
        if (!end) {
            break;
        }
        list = end + 1;
    }
    return FALSE;
}

static int engine_cgroup_decode_path(
    const char* source,
    char* destination,
    size_t destination_size) {
    size_t source_length;
    size_t input = 0U;
    size_t output = 0U;

    if (!source || !destination || !destination_size) {
        return -1;
    }
    source_length = strlen(source);
    while (input < source_length) {
        unsigned int value;

        if (source[input] == '\\' &&
            input + 3U < source_length &&
            source[input + 1U] >= '0' &&
            source[input + 1U] <= '7' &&
            source[input + 2U] >= '0' &&
            source[input + 2U] <= '7' &&
            source[input + 3U] >= '0' &&
            source[input + 3U] <= '7') {
            value = (unsigned int)(source[input + 1U] - '0') * 64U +
                (unsigned int)(source[input + 2U] - '0') * 8U +
                (unsigned int)(source[input + 3U] - '0');
            input += 4U;
        } else {
            value = (unsigned char)source[input++];
        }
        if (!value || output + 1U >= destination_size) {
            return -1;
        }
        destination[output++] = (char)value;
    }
    if (!output || destination[0] != '/') {
        return -1;
    }
    destination[output] = '\0';
    return 0;
}

static int engine_cgroup_membership(
    char* hierarchy_path,
    size_t hierarchy_size,
    anbool* unified) {
    char line[ENGINE_CGROUP_PATH_SIZE + 256U];
    char unified_path[ENGINE_CGROUP_PATH_SIZE] = {0};
    FILE* file;

    if (!hierarchy_path || !hierarchy_size || !unified) {
        return -1;
    }
    hierarchy_path[0] = '\0';
    file = fopen("/proc/self/cgroup", "r");
    if (!file) {
        return -1;
    }
    while (fgets(line, sizeof(line), file)) {
        char* first = strchr(line, ':');
        char* second = first ? strchr(first + 1, ':') : NULL;
        char* newline = strchr(line, '\n');
        char* selected = NULL;

        if (!newline && !feof(file)) {
            fclose(file);
            return -1;
        }
        if (newline) {
            *newline = '\0';
        }
        if (!first || !second || second[1] != '/') {
            continue;
        }
        *first = '\0';
        *second = '\0';
        if (!first[1] && !strcmp(line, "0")) {
            selected = unified_path;
        } else if (engine_cgroup_list_contains(
                       first + 1, "memory")) {
            selected = hierarchy_path;
        }
        if (selected) {
            size_t length = strlen(second + 1);

            if (!length || length >= ENGINE_CGROUP_PATH_SIZE) {
                fclose(file);
                return -1;
            }
            memcpy(selected, second + 1, length + 1U);
        }
    }
    fclose(file);
    if (hierarchy_path[0]) {
        *unified = FALSE;
        return 1;
    }
    if (unified_path[0]) {
        size_t length = strlen(unified_path);

        if (length >= hierarchy_size) {
            return -1;
        }
        memcpy(hierarchy_path, unified_path, length + 1U);
        *unified = TRUE;
        return 1;
    }
    return 0;
}

static anbool engine_cgroup_path_contains(
    const char* root,
    const char* path) {
    size_t length;

    if (!root || !path || root[0] != '/' || path[0] != '/') {
        return FALSE;
    }
    if (!strcmp(root, "/")) {
        return TRUE;
    }
    length = strlen(root);
    return !strncmp(root, path, length) &&
        (path[length] == '\0' || path[length] == '/');
}

static int engine_cgroup_mount(
    const char* hierarchy_path,
    anbool unified,
    char* mount_point,
    char* leaf_path) {
    char line[ENGINE_CGROUP_LINE_SIZE];
    size_t best_root_length = 0U;
    FILE* file;

    file = fopen("/proc/self/mountinfo", "r");
    if (!file) {
        return -1;
    }
    mount_point[0] = '\0';
    leaf_path[0] = '\0';
    while (fgets(line, sizeof(line), file)) {
        char encoded_root[ENGINE_CGROUP_PATH_SIZE];
        char encoded_mount[ENGINE_CGROUP_PATH_SIZE];
        char root[ENGINE_CGROUP_PATH_SIZE];
        char mount[ENGINE_CGROUP_PATH_SIZE];
        char filesystem[32];
        char super_options[ENGINE_CGROUP_PATH_SIZE];
        char candidate[ENGINE_CGROUP_PATH_SIZE];
        char* separator;
        const char* relative;
        size_t root_length;
        int prefix_length = 0;
        int candidate_length;

        if (!strchr(line, '\n') && !feof(file)) {
            continue;
        }
        if (sscanf(line, "%*s %*s %*s %4095s %4095s %n",
                   encoded_root, encoded_mount, &prefix_length) != 2) {
            continue;
        }
        separator = strstr(line + prefix_length, " - ");
        if (!separator ||
            sscanf(separator + 3, "%31s %*s %4095s",
                   filesystem, super_options) != 2) {
            continue;
        }
        if ((unified && strcmp(filesystem, "cgroup2")) ||
            (!unified &&
             (strcmp(filesystem, "cgroup") ||
              !engine_cgroup_list_contains(
                  super_options, "memory"))) ||
            engine_cgroup_decode_path(
                encoded_root, root, sizeof(root)) ||
            engine_cgroup_decode_path(
                encoded_mount, mount, sizeof(mount)) ||
            strcmp(root, "/") ||
            !engine_cgroup_path_contains(root, hierarchy_path)) {
            continue;
        }
        root_length = strlen(root);
        relative = !strcmp(root, "/")
            ? hierarchy_path : hierarchy_path + root_length;
        if (!relative[0] || !strcmp(relative, "/")) {
            candidate_length = snprintf(
                candidate, sizeof(candidate), "%s", mount);
        } else if (!strcmp(mount, "/")) {
            candidate_length = snprintf(
                candidate, sizeof(candidate), "%s", relative);
        } else {
            candidate_length = snprintf(
                candidate, sizeof(candidate), "%s%s", mount, relative);
        }
        if (candidate_length < 0 ||
            (size_t)candidate_length >= sizeof(candidate) ||
            root_length < best_root_length) {
            continue;
        }
        memcpy(mount_point, mount, strlen(mount) + 1U);
        memcpy(leaf_path, candidate, strlen(candidate) + 1U);
        best_root_length = root_length;
    }
    fclose(file);
    return mount_point[0] && leaf_path[0] ? 0 : -1;
}

static int engine_cgroup_apply_limits(
    const char* mount_point,
    const char* leaf_path,
    anbool unified,
    size_t* capacity_bytes,
    size_t* available_bytes) {
    char current[ENGINE_CGROUP_PATH_SIZE];
    const char* limit_name = unified
        ? "memory.max" : "memory.limit_in_bytes";
    const char* usage_name = unified
        ? "memory.current" : "memory.usage_in_bytes";
    size_t mount_length = strlen(mount_point);
    int found = 0;

    if (strlen(leaf_path) >= sizeof(current) ||
        !engine_cgroup_path_contains(mount_point, leaf_path)) {
        return -1;
    }
    memcpy(current, leaf_path, strlen(leaf_path) + 1U);
    while (1) {
        char limit_path[ENGINE_CGROUP_PATH_SIZE];
        char usage_path[ENGINE_CGROUP_PATH_SIZE];
        size_t limit;
        size_t usage;
        int limit_status;
        int usage_status;
        int limit_length = snprintf(
            limit_path, sizeof(limit_path),
            "%s/%s", current, limit_name);
        int usage_length = snprintf(
            usage_path, sizeof(usage_path),
            "%s/%s", current, usage_name);

        if (limit_length <= 0 || usage_length <= 0 ||
            (size_t)limit_length >= sizeof(limit_path) ||
            (size_t)usage_length >= sizeof(usage_path)) {
            return -1;
        }
        errno = 0;
        limit_status =
            engine_read_memory_limit(limit_path, &limit);
        if (limit_status < 0) {
            if (errno != ENOENT) {
                return -1;
            }
        } else {
            found = 1;
            if (!limit_status) {
                errno = 0;
                usage_status = engine_read_memory_limit(
                    usage_path, &usage);
                if (usage_status) {
                    return -1;
                }
                *capacity_bytes =
                    MIN(*capacity_bytes, limit);
                *available_bytes = MIN(*available_bytes,
                    usage < limit ? limit - usage : 0U);
            }
        }
        if (!strcmp(current, mount_point)) {
            break;
        }
        {
            char* slash = strrchr(current, '/');

            if (!slash) {
                return -1;
            }
            if (!strcmp(mount_point, "/") && slash == current) {
                current[1] = '\0';
                continue;
            }
            if ((size_t)(slash - current) < mount_length) {
                return -1;
            }
            *slash = '\0';
        }
    }
    return found ? 1 : -1;
}
#endif

static int engine_limit_memory_by_cgroup(
    size_t* capacity_bytes,
    size_t* available_bytes) {
    if (!capacity_bytes || !available_bytes) {
        return -1;
    }
#ifdef __linux__
    {
        char hierarchy_path[ENGINE_CGROUP_PATH_SIZE];
        char mount_point[ENGINE_CGROUP_PATH_SIZE];
        char leaf_path[ENGINE_CGROUP_PATH_SIZE];
        anbool unified;
        int status;

        status = engine_cgroup_membership(
            hierarchy_path,
            sizeof(hierarchy_path),
            &unified);
        if (status <= 0) {
            return status;
        }
        if (engine_cgroup_mount(
                hierarchy_path,
                unified,
                mount_point,
                leaf_path)) {
            return -1;
        }
        status = engine_cgroup_apply_limits(
            mount_point,
            leaf_path,
            unified,
            capacity_bytes,
            available_bytes);
        return status < 0 ? -1 : 0;
    }
#else
    return 0;
#endif
}

static void engine_limit_memory_by_address_space(
    size_t page_size,
    size_t* capacity_bytes,
    size_t* available_bytes) {
    struct rlimit address_limit;
    uintmax_t pages = 0U;
    size_t current_bytes = 0U;
    size_t limit_bytes;
    FILE* file;

    if (!page_size || !capacity_bytes || !available_bytes ||
        getrlimit(RLIMIT_AS, &address_limit) ||
        address_limit.rlim_cur == RLIM_INFINITY ||
        (uintmax_t)address_limit.rlim_cur > SIZE_MAX) {
        return;
    }
    limit_bytes = (size_t)address_limit.rlim_cur;
    file = fopen("/proc/self/statm", "r");
    if (file) {
        if (fscanf(file, "%ju", &pages) == 1 &&
            pages <= SIZE_MAX / page_size) {
            current_bytes = (size_t)pages * page_size;
        }
        fclose(file);
    }
    *capacity_bytes = MIN(*capacity_bytes, limit_bytes);
    *available_bytes = MIN(
        *available_bytes,
        current_bytes < limit_bytes
            ? limit_bytes - current_bytes : 0U);
}

static anbool engine_job_local_residency_enabled(void) {
    return FALSE;
}

index_residency_t* engine_index_residency_begin(
    engine_t* engine,
    const onefield_t* bp) {
    index_residency_t* service = NULL;
    index_residency_stats_t residency_stats;
    size_t cohort_bytes;
    size_t cohort_files;
    size_t available_bytes;
    size_t physical_bytes;
    size_t physical_headroom;
    size_t physical_full_limit;
    size_t available_headroom;
    size_t available_full_limit;
    size_t worker_headroom;
    unsigned int lanes;
    long physical_pages;
    long page_size;
    int i;

    if (!engine || !bp || bp->index_shard_workers <= 1) {
        return NULL;
    }

    /*
     * The current service copies the complete cohort inside the solve wall
     * clock, blocks without solver-limit polling, and marks reclaimable memfd
     * pages as permanently resident. Keep exact demand delivery authoritative
     * until residency has a persistent pre-job lifecycle and recoverable page
     * state.
     */
    if (!engine_job_local_residency_enabled()) {
        logverb("[index-residency] mode=exact-demand "
                "reason=job-local-residency-quarantined\n");
        return NULL;
    }

    if (engine_index_cohort_measure(
            engine, &cohort_bytes, &cohort_files) ||
        !cohort_files || !cohort_bytes ||
        engine_available_memory(&available_bytes)) {
        return NULL;
    }
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    physical_pages = sysconf(_SC_PHYS_PAGES);
    page_size = sysconf(_SC_PAGESIZE);
    if (physical_pages <= 0 || page_size <= 0 ||
        (uintmax_t)physical_pages >
            (uintmax_t)SIZE_MAX / (uintmax_t)page_size) {
        return NULL;
    }
    physical_bytes =
        (size_t)physical_pages * (size_t)page_size;
#else
    (void)physical_pages;
    (void)page_size;
    return NULL;
#endif

    /*
     * Whole-file residency is useful only when every eligible source can
     * remain resident for the job. Partial whole-file LRU would copy and
     * discard broad data for sparse queries, recreating I/O amplification.
     * Capacity, current availability, cgroups and address-space limits are
     * all advisory admission guards; refusal keeps exact delivery unchanged.
     */
    if (engine_limit_memory_by_cgroup(
            &physical_bytes, &available_bytes)) {
        logverb("[index-residency] mode=exact-demand "
                "reason=cgroup-admission-unavailable\n");
        return NULL;
    }
    engine_limit_memory_by_address_space(
        (size_t)page_size,
        &physical_bytes,
        &available_bytes);
    physical_headroom = physical_bytes / 4U;
    physical_full_limit = physical_bytes - physical_headroom;
    worker_headroom = 512U * 1024U * 1024U;
    if ((size_t)bp->index_shard_workers <=
        (SIZE_MAX - worker_headroom) /
            (128U * 1024U * 1024U)) {
        worker_headroom +=
            (size_t)bp->index_shard_workers *
            (128U * 1024U * 1024U);
    }
    available_headroom = MAX(
        available_bytes / 5U,
        worker_headroom);
    available_full_limit =
        available_headroom < available_bytes
            ? available_bytes - available_headroom : 0U;
    if (cohort_bytes > physical_full_limit ||
        cohort_bytes > available_full_limit) {
        logverb(
            "[index-residency] mode=exact-demand "
            "reason=cohort-does-not-fit files=%zu "
            "cohort_bytes=%zu capacity_bytes=%zu "
            "available_bytes=%zu\n",
            cohort_files,
            cohort_bytes,
            physical_bytes,
            available_bytes);
        return NULL;
    }

    lanes = bp->index_shard_workers >= 4 ? 2U : 1U;
    if (index_residency_start(
            cohort_bytes, lanes, &service)) {
        return NULL;
    }
    for (i = 0; i < pl_size(engine->indexes); i++) {
        const index_t* index = pl_get(engine->indexes, i);
        index_residency_result_t prepare_status;

        if (!engine_index_residency_eligible(index)) {
            continue;
        }
        prepare_status = index_residency_prepare(
            service,
            index->indexfn,
            INDEX_RESIDENCY_PRIORITY_LOOKAHEAD);
        if (prepare_status != INDEX_RESIDENCY_ACCEPTED) {
            (void)index_residency_stop(service);
            logverb(
                "[index-residency] mode=exact-demand "
                "reason=prepare-fallback\n");
            return NULL;
        }
    }
    if (index_residency_drain(service) ||
        index_residency_get_stats(
            service, &residency_stats) ||
        residency_stats.ready_entries != cohort_files ||
        residency_stats.ready_bytes != cohort_bytes ||
        residency_stats.resident_bytes != cohort_bytes ||
        residency_stats.loading_entries ||
        residency_stats.loading_bytes ||
        residency_stats.failed_entries) {
        (void)index_residency_stop(service);
        logverb(
            "[index-residency] mode=exact-demand "
            "reason=full-cohort-not-ready\n");
        return NULL;
    }
    if (index_bind_residency_service(service)) {
        (void)index_residency_stop(service);
        logverb(
            "[index-residency] mode=exact-demand "
            "reason=concurrent-binding\n");
        return NULL;
    }
    logverb(
        "[index-residency] mode=full-cohort files=%zu "
        "cohort_bytes=%zu budget_bytes=%zu capacity_bytes=%zu "
        "available_bytes=%zu lanes=%u\n",
        cohort_files,
        cohort_bytes,
        cohort_bytes,
        physical_bytes,
        available_bytes,
        lanes);
    return service;
}
