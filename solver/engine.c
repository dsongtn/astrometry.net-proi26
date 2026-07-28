/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/**
 * Accepts an augmented xylist that describes a field or set of fields to solve.
 * Reads a config file to find local indices, and merges information about the
 * indices with the job description to create an input file for 'onefield'.
 * Runs and merges the results.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <libgen.h>
#include <getopt.h>
#include <dirent.h>
#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <sys/resource.h>
#include <unistd.h>

#include "math.h"

#include "an-bool.h"
#include "anqfits.h"
#include "astrometry/index_shard.h"
#include "astrometry/index_residency.h"
#include "bl.h"
#include "engine.h"
#include "errors.h"
#include "fileutils.h"
#include "fitsioutils.h"
#include "healpix.h"
#include "indexset.h"
#include "ioutils.h"
#include "log.h"
#include "mathutil.h"
#include "multiindex.h"
#include "onefield.h"
#include "os-features.h"
#include "sip-utils.h"
#include "solver.h"
#include "solverutils.h"
#include "tic.h"
#include "index_shard_config.h"
#include "engine_internal.h"

void engine_pass_cursor_init(engine_pass_cursor_t* cursor) {
    if (!cursor) {
        return;
    }
    memset(cursor, 0, sizeof(*cursor));
}

anbool engine_pass_cursor_next(const job_t* job,
                               double default_lower,
                               double default_upper,
                               engine_pass_cursor_t* cursor,
                               engine_pass_t* pass) {
    size_t depth_count;
    size_t scale_count;
    int raw_start;
    int raw_end;
    double raw_lower;
    double raw_upper;

    if (!job || !job->depths || !job->scales || !cursor || !pass) {
        return FALSE;
    }
    depth_count = (size_t)il_size(job->depths) / 2U;
    scale_count = (size_t)dl_size(job->scales) / 2U;
    if (!depth_count || !scale_count ||
        cursor->next_depth_index >= depth_count) {
        return FALSE;
    }

    memset(pass, 0, sizeof(*pass));
    pass->ordinal = cursor->next_ordinal;
    pass->depth_index = cursor->next_depth_index;
    pass->scale_index = cursor->next_scale_index;

    raw_start = il_get(job->depths, pass->depth_index * 2U);
    raw_end = il_get(job->depths, pass->depth_index * 2U + 1U);
    if (raw_start < 0 || raw_end < 0) {
        return FALSE;
    }
    pass->startobj = raw_start ? raw_start - 1 : 0;
    /*
     * The user-facing upper bound is inclusive and one-based. Its numeric
     * value is therefore already the zero-based exclusive bound. Zero is the
     * native open-upper sentinel and must be written on every pass.
     */
    pass->endobj = raw_end;

    raw_lower = dl_get(job->scales, pass->scale_index * 2U);
    raw_upper = dl_get(job->scales, pass->scale_index * 2U + 1U);
    pass->funits_lower =
        raw_lower == 0.0 ? default_lower : raw_lower;
    pass->funits_upper =
        raw_upper == 0.0 ? default_upper : raw_upper;

    cursor->next_scale_index++;
    cursor->next_ordinal++;
    if (cursor->next_scale_index >= scale_count) {
        cursor->next_scale_index = 0U;
        cursor->next_depth_index++;
    }
    return TRUE;
}

void engine_pass_apply(solver_t* solver, const engine_pass_t* pass) {
    if (!solver || !pass) {
        return;
    }
    solver->startobj = pass->startobj;
    solver->endobj = pass->endobj;
    solver->funits_lower = pass->funits_lower;
    solver->funits_upper = pass->funits_upper;
}

void engine_add_search_path(engine_t* engine, const char* path) {
    sl_append(engine->index_paths, path);
}

char* engine_find_index(engine_t* engine, const char* name) {
    int j;

    for (j=-1; j<(int)sl_size(engine->index_paths); j++) {
        char* path;
        if (j == -1)
            if (strlen(name) && name[0] == '/') {
                // try as an absolute filename.
                path = strdup(name);
            } else {
                continue;
            }
        else
            asprintf_safe(&path, "%s/%s", sl_get(engine->index_paths, j), name);

        logverb("Trying path %s...\n", path);
        if (index_is_file_index(path))
            return path;
        free(path);
    }
    return NULL;
}

int engine_autoindex_search_paths(engine_t* engine) {
    int i;
    // Search the paths specified and add any indexes that are found.
    for (i=0; i<sl_size(engine->index_paths); i++) {
        char* path = sl_get(engine->index_paths, i);
        DIR* dir = opendir(path);
        sl* tryinds;
        int j;
        if (!dir) {
            SYSERROR("Warning: failed to open index directory: \"%s\"\n", path);
            continue;
        }
        logverb("Auto-indexing directory \"%s\" ...\n", path);
        tryinds = sl_new(16);
        while (1) {
            struct dirent* de;
            char* name;
            char* fullpath;
            char* err;
            anbool ok;
            errno = 0;
            de = readdir(dir);
            if (!de) {
                if (errno)
                    SYSERROR("Failed to read entry from directory \"%s\"", path);
                break;
            }
            name = de->d_name;
            asprintf_safe(&fullpath, "%s/%s", path, name);
            if (path_is_dir(fullpath)) {
                logverb("Skipping directory %s\n", fullpath);
                free(fullpath);
                continue;
            }

            logverb("Checking file \"%s\"\n", fullpath);
            errors_start_logging_to_string();
            ok = index_is_file_index(fullpath);
            err = errors_stop_logging_to_string(": ");
            if (!ok) {
                logverb("File is not an index: %s\n", err);
                free(err);
                free(fullpath);
                continue;
            }
            free(err);

            sl_insert_sorted_nocopy(tryinds, fullpath);
        }
        closedir(dir);

        // add them in reverse order... (why?)
        for (j=sl_size(tryinds)-1; j>=0; j--) {
            char* path = sl_get(tryinds, j);
            logverb("Trying to add index \"%s\".\n", path);
            if (engine_add_index(engine, path))
                logmsg("Failed to add index \"%s\".\n", path);
        }
        sl_free2(tryinds);
    }
    return 0;
}

static int add_index(engine_t* engine, index_t* ind) {
    int k;
    // check that an index with the same id and healpix isn't already listed.
    for (k=0; k<pl_size(engine->indexes); k++) {
        index_t* m = pl_get(engine->indexes, k);
        if (m->indexid == ind->indexid &&
            m->healpix == ind->healpix) {
            logmsg("Warning: encountered two index files with the same INDEXID = %i and HEALPIX = %i: \"%s\" and \"%s\".  Keeping both.\n",
                   m->indexid, m->healpix, m->indexname, ind->indexname);
            //index_free(ind);
            //return 0;
        }
    }

    pl_append(engine->indexes, ind);

    // <= smallest we've seen?
    if (ind->index_scale_lower < engine->sizesmallest) {
        engine->sizesmallest = ind->index_scale_lower;
        bl_remove_all(engine->ismallest);
        il_append(engine->ismallest, pl_size(engine->indexes) - 1);
    } else if (ind->index_scale_lower == engine->sizesmallest) {
        il_append(engine->ismallest, pl_size(engine->indexes) - 1);
    }

    // >= largest we've seen?
    if (ind->index_scale_upper > engine->sizebiggest) {
        engine->sizebiggest = ind->index_scale_upper;
        bl_remove_all(engine->ibiggest);
        il_append(engine->ibiggest, pl_size(engine->indexes) - 1);
    } else if (ind->index_scale_upper == engine->sizebiggest) {
        il_append(engine->ibiggest, pl_size(engine->indexes) - 1);
    }
    return 0;
}

int engine_add_index(engine_t* engine, char* path) {
    int k;
    index_t* ind = NULL;
    char* quadpath = index_get_quad_filename(path);
    char* base = basename_safe(quadpath);
    double t0;
    free(quadpath);

    // check that an index with the same filename hasn't already been added.
    for (k=0; k<pl_size(engine->indexes); k++) {
        ind = pl_get(engine->indexes, k);
        // ind->indexname is a path to the quad filename; strip off directory component.
        char* mbase = basename_safe(ind->indexname);
        anbool eq = streq(base, mbase);
        free(mbase);
        if (eq) {
            logmsg("Warning: we've already seen an index with the same name: \"%s\".  Adding it anyway...\n", ind->indexname);
            //free(base);
            //return 0;
        }
    }
    free(base);

    t0 = timenow();
    /*
     * Ordinary registration is always metadata-only. Legacy grouped mode
     * still loads all selected filename-owned indexes together inside
     * onefield; it no longer needs every configured payload resident before
     * scale and sky selection.
     */
    ind = index_load(path, INDEX_ONLY_LOAD_METADATA, NULL);
    debug("index_load(\"%s\") took %g ms\n", path, 1000 * (timenow() - t0));
    if (!ind) {
        ERROR("Failed to load index from path %s", path);
        return -1;
    }
    if (add_index(engine, ind)) {
        ERROR("Failed to add index \"%s\"", path);
        return -1;
    }
    pl_append(engine->free_indexes, ind);
    return 0;
}
int engine_parse_config_file(engine_t* engine, const char* fn) {
    FILE* fconf;
    int rtn;
    fconf = fopen(fn, "r");
    if (!fconf) {
        SYSERROR("Failed to open config file \"%s\"", fn);
        return -1;
    }
    rtn = engine_parse_config_file_stream(engine, fconf);
    fclose(fconf);
    return rtn;
}

int engine_parse_config_file_stream(engine_t* engine, FILE* fconf) {
    sl* indices = sl_new(16);
    sl* indexsets = sl_new(16);
    sl* mindices = sl_new(16);
    anbool auto_index = FALSE;
    int i;
    int rtn = 0;

    while (1) {
        char buffer[10240];
        char* nextword;
        char* line;
        if (!fgets(buffer, sizeof(buffer), fconf)) {
            if (feof(fconf))
                break;
            SYSERROR("Failed to read a line from the config file");
            rtn = -1;
            goto done;
        }
        line = buffer;
        // strip off newline
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';
        // skip leading whitespace:
        while (*line && isspace((unsigned)(*line)))
            line++;
        // skip comments
        if (line[0] == '#')
            continue;
        // skip blank lines.
        if (line[0] == '\0')
            continue;

        if (is_word(line, "index ", &nextword)) {
            // don't try to find the index yet - because search paths may be
            // added later.
            sl_append(indices, nextword);
        } else if (is_word(line, "indexset ", &nextword)) {
            // don't try to find the index yet - because search paths may be
            // added later.
            sl_append(indexsets, nextword);
        } else if (is_word(line, "multiindex ", &nextword)) {
            // don't try to find the index yet - because search paths may be
            // added later.
            sl_append(mindices, nextword);
        } else if (is_word(line, "autoindex", &nextword)) {
            auto_index = TRUE;
        } else if (is_word(line, "inparallel", &nextword)) {
            engine->inparallel = TRUE;
        } else if (is_word(line, "minwidth ", &nextword)) {
            engine->minwidth = atof(nextword);
        } else if (is_word(line, "maxwidth ", &nextword)) {
            engine->maxwidth = atof(nextword);
        } else if (is_word(line, "cpulimit ", &nextword)) {
            engine->cpulimit = atof(nextword);
        } else if (is_word(line, "p_workers ", &nextword) ||
                   is_word(line, "index_shard_workers ", &nextword)) {
            int available_cpus = index_shard_config_available_cpus();
            int requested_workers;

            if (index_shard_config_parse_workers(nextword,
                                                 available_cpus,
                                                 &requested_workers)) {
                ERROR("Invalid p_workers value \"%s\": "
                      "expected \"auto\" or an integer from 1 through %i",
                      nextword,
                      available_cpus);
                rtn = -1;
                goto done;
            }

            engine->index_shard_workers_config = requested_workers;
            engine->index_shard_workers_config_set = TRUE;
        } else if (is_word(line, "depths ", &nextword)) {
            if (parse_depth_string(engine->default_depths, nextword)) {
                rtn = -1;
                goto done;
            }
        } else if (is_word(line, "add_path ", &nextword)) {
            engine_add_search_path(engine, nextword);
        } else {
            ERROR("Didn't understand this config file line: \"%s\"", line);
            // unknown config line is a firing offense
            rtn = -1;
            goto done;
        }
    }

    for (i=0; i<sl_size(indices); i++) {
        char* ind = sl_get(indices, i);
        char* path;
        logverb("Trying index %s...\n", ind);

        path = engine_find_index(engine, ind);
        if (!path) {
            logmsg("Couldn't find index \"%s\".\n", ind);
            rtn = -1;
            goto done;
        }
        if (engine_add_index(engine, path))
            logmsg("Failed to add index \"%s\".\n", path);
        free(path);
    }

    for (i=0; i<sl_size(indexsets); i++) {
        char* ind = sl_get(indexsets, i);
        pl* indexes = pl_new(16);
        int i, j;

        logverb("Trying index-set %s...\n", ind);
        indexset_get(ind, indexes);
        if (bl_size(indexes) == 0) {
            ERROR("Unknown index-set \"%s\"", ind);
            rtn = -1;
            goto done;
        }
        // See which index files in the set exist
        // NOTE, no i++ here -- we only advance conditionally
        for (i=0; i<bl_size(indexes);) {
            index_t* indx = pl_get(indexes, i);
            for (j=0; j<(int)sl_size(engine->index_paths); j++) {
                char* path;
                asprintf_safe(&path, "%s/%s", sl_get(engine->index_paths, j), indx->indexname);
                if (file_readable(path)) {
                    indx->indexfn = path;
                    break;
                }
                free(path);
            }
            if (!indx->indexfn) {
                logverb("Did not find file for index name \"%s\"\n", indx->indexname);
                index_free(indx);
                bl_remove_index(indexes, i);
                continue;
            }
            // Found an index where the file exists!
            // bypass engine_add_index(), which opens the file...
            if (add_index(engine, indx)) {
                ERROR("Failed to add index \"%s\"", indx->indexfn);
                rtn = -1;
                goto done;
            }
            pl_append(engine->free_indexes, indx);
            logverb("Added index %s from indexset %s\n", indx->indexfn, ind);

            i++;
        }
        pl_free(indexes);
    }

    for (i=0; i<sl_size(mindices); i++) {
        char* ind = sl_get(mindices, i);
        char* path;
        char* skdt;
        char* skdtpath;
        int j;
        sl* words = sl_split(NULL, ind, " ");
        multiindex_t* mi;

        if (sl_size(words) < 2) {
            logmsg("Config line 'multiindex' must be followed by skdt and inds\n");
            rtn = -1;
            goto done;
        }
        skdt = sl_get(words, 0);
        sl_remove(words, 0);
        {
            char* s = sl_join(words, " / ");
            logverb("Trying multi-index %s + %s...\n", skdt, s);
            free(s);
        }
        skdtpath = engine_find_index(engine, skdt);
        if (!skdtpath) {
            logmsg("Couldn't find skdt \"%s\".\n", skdt);
            rtn = -1;
            goto done;
        }
        for (j=0; j<sl_size(words); j++) {
            ind = sl_get(words, j);
            path = engine_find_index(engine, ind);
            if (!path) {
                logmsg("Couldn't find index \"%s\".\n", ind);
                rtn = -1;
                goto done;
            }
            sl_set(words, j, path);
            // sl_set makes a copy.
            free(path);
        }

        mi = multiindex_open(skdtpath, words, 0);
        if (!mi) {
            char* s = sl_join(words, " / ");
            logerr("Failed to open multiindex: %s + %s\n", skdt, s);
            free(s);
            rtn = -1;
            goto done;
        }
        for (j=0; j<multiindex_n(mi); j++) {
            index_t* ind = multiindex_get(mi, j);
            if (add_index(engine, ind)) {
                ERROR("Failed to add index \"%s\"", sl_get(words, j));
                return -1;
            }
        }
        pl_append(engine->free_mindexes, mi);
        sl_free2(words);
        free(skdt);
        free(skdtpath);
    }

    if (auto_index) {
        engine_autoindex_search_paths(engine);
    }

 done:
    sl_free2(indices);
    sl_free2(mindices);
    sl_free2(indexsets);
    return rtn;
}

static job_t* job_new() {
    job_t* job = calloc(1, sizeof(job_t));
    if (!job) {
        SYSERROR("Failed to allocate a new job_t.");
        return NULL;
    }
    job->scales = dl_new(8);
    job->depths = il_new(8);
    job->index_shard_workers_override = INDEX_SHARD_WORKERS_UNSET;
    return job;
}

static anbool engine_index_residency_eligible(
    const index_t* index) {
    return index && index->indexfn &&
        !index->codekd && !index->quads && !index->starkd;
}

void job_free(job_t* job) {
    if (!job)
        return;
    dl_free(job->scales);
    il_free(job->depths);
    free(job);
}

static double job_imagew(job_t* job) {
    return job->bp.solver.field_maxx;
}
static double job_imageh(job_t* job) {
    return job->bp.solver.field_maxy;
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

static index_residency_t* engine_index_residency_begin(
    engine_t* engine,
    const onefield_t* bp) {
    index_residency_t* service = NULL;
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

    if (!engine || !bp || bp->index_shard_workers <= 1 ||
        engine_index_cohort_measure(
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
            INDEX_RESIDENCY_PRIORITY_SPECULATIVE);
        if (prepare_status != INDEX_RESIDENCY_ACCEPTED) {
            (void)index_residency_stop(service);
            logverb(
                "[index-residency] mode=exact-demand "
                "reason=prepare-fallback\n");
            return NULL;
        }
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

int engine_run_job(engine_t* engine, job_t* job) {
    onefield_t* bp = &(job->bp);
    solver_t* sp = &(bp->solver);

    int rtn = 0;
    double app_min_default;
    double app_max_default;
    double engine_wall_start = monotonic_seconds();
    double pool_start_seconds = 0.0;
    double pool_stop_seconds = 0.0;
    anbool index_shard_pool_started = FALSE;
    anbool legacy_grouped =
        engine->inparallel && !job->index_shard_workers_controlled;
    index_residency_t* residency = NULL;
    engine_pass_cursor_t pass_cursor;
    engine_pass_t pass;

    if (onefield_is_run_obsolete(bp, sp)) {
        goto finish;
    }
    // SECTION INDEX-SHARD: engine-lifecycle
    bp->time_total_start = monotonic_seconds();
    bp->cpu_total_start = get_cpu_usage();
    bp->indexes_inparallel = legacy_grouped;

    app_min_default = deg2arcsec(engine->minwidth) / job_imagew(job);
    app_max_default = deg2arcsec(engine->maxwidth) / job_imagew(job);

    if (job->use_radec_center) {
        logmsg("Only searching for solutions within %g degrees of RA,Dec (%g,%g)\n",
               job->search_radius, job->ra_center, job->dec_center);
        solver_set_radec(sp, job->ra_center, job->dec_center, job->search_radius);
    }

    if (onefield_job_field_cache_begin(bp)) {
        ERROR("Failed to initialize job field cache");
        rtn = -1;
        goto finish;
    }

    residency = engine_index_residency_begin(engine, bp);

    if (index_shard_pthread_enabled(bp) && !legacy_grouped) {
        double pool_wall_start = monotonic_seconds();

        if (index_shard_pool_start(bp, sp)) {
            ERROR("Failed to start parallel solver pool");
            rtn = -1;
            goto finish;
        }

        pool_start_seconds =
            monotonic_seconds() - pool_wall_start;
        index_shard_pool_started = TRUE;
    }

    engine_pass_cursor_init(&pass_cursor);
    while (engine_pass_cursor_next(
               job,
               app_min_default,
               app_max_default,
               &pass_cursor,
               &pass)) {
            double fmin, fmax;
            double app_max, app_min;
            int k;
            il* indexlist;
            il* selectedlist;
            anbool selected_loaded_index = FALSE;
            anbool pass_limit_reached = FALSE;

            /*
             * Index selection and materialization can be expensive and fault
             * mapped metadata.  A job budget is terminal across the whole
             * pass sequence; never start another pass after it expires.
             */
            if (onefield_check_total_limits(bp)) {
                break;
            }

            // arcsec per pixel range
            app_min = pass.funits_lower;
            app_max = pass.funits_upper;
            engine_pass_apply(sp, &pass);
            bp->engine_pass_ordinal = pass.ordinal;
            bp->engine_depth_index = pass.depth_index;
            bp->engine_scale_index = pass.scale_index;
            logverb("[engine-pass] state=begin ordinal=%zu "
                    "depth_index=%zu scale_index=%zu "
                    "startobj=%i endobj=%i "
                    "funits_lower=%.17g funits_upper=%.17g\n",
                    pass.ordinal,
                    pass.depth_index,
                    pass.scale_index,
                    pass.startobj,
                    pass.endobj,
                    pass.funits_lower,
                    pass.funits_upper);

            // minimum quad size to try (in pixels)
            sp->quadsize_min = bp->quad_size_fraction_lo *
                MIN(job_imagew(job), job_imageh(job));

            // range of quad sizes that could be found in the field,
            // in arcsec.
            // the hypotenuse...
            fmax = bp->quad_size_fraction_hi *
                hypot(job_imagew(job), job_imageh(job)) * app_max;
            fmin = sp->quadsize_min * app_min;

            // Select the indices that should be checked.
            indexlist = il_new(16);
            for (k = 0; k < pl_size(engine->indexes); k++) {
                index_t* index = pl_get(engine->indexes, k);
                if (!index_overlaps_scale_range(index, fmin, fmax))
                    continue;
                il_append(indexlist, k);
            }

            // Use the (list of) smallest or largest indices if no other one fits.
            if (!il_size(indexlist)) {
                il* list = NULL;
                if (fmin > engine->sizebiggest) {
                    list = engine->ibiggest;
                } else if (fmax < engine->sizesmallest) {
                    list = engine->ismallest;
                } else {
                    assert(0);
                }
                il_append_list(indexlist, list);
            }

            selectedlist = il_new(il_size(indexlist));
            for (k=0; k<il_size(indexlist); k++) {
                int ii = il_get(indexlist, k);
                index_t* index = pl_get(engine->indexes, ii);
                anbool inrange = TRUE;
                if (job->use_radec_center) {
                    inrange = index_is_within_range(index, job->ra_center, job->dec_center, job->search_radius);
                }
                if (!inrange) {
                    logverb("Not using index %s because it's not within %g degrees of (RA,Dec) = (%g,%g)\n",
                            index->indexname, job->search_radius, job->ra_center, job->dec_center);
                    continue;
                }
                il_append(selectedlist, ii);
                if (index->starkd && index->quads && index->codekd) {
                    selected_loaded_index = TRUE;
                }
            }

            il_free(indexlist);
            if (onefield_check_total_limits(bp)) {
                il_free(selectedlist);
                logverb("[engine-pass] state=end ordinal=%zu "
                        "reason=limit-before-materialization\n",
                        pass.ordinal);
                break;
            }
            /*
             * onefield keeps filename and loaded handles in separate lists.
             * If a pass contains a borrowed multiindex component, materialize
             * every ordinary member into the loaded list so their original
             * interleaved order is preserved exactly.
             */
            for (k = 0; k < il_size(selectedlist); k++) {
                int ii = il_get(selectedlist, k);
                index_t* index = pl_get(engine->indexes, ii);

                if (!selected_loaded_index) {
                    onefield_add_index(bp, index->indexfn);
                } else if (index->starkd &&
                           index->quads &&
                           index->codekd) {
                    onefield_add_loaded_index(bp, index);
                } else {
                    index_t* owned_index =
                        index_load(index->indexfn, 0, NULL);

                    if (!owned_index) {
                        ERROR("Failed to load selected index %s",
                              index->indexfn);
                        il_free(selectedlist);
                        rtn = -1;
                        goto finish;
                    }
                    onefield_add_owned_index(bp, owned_index);
                }

                if (onefield_check_total_limits(bp)) {
                    pass_limit_reached = TRUE;
                    break;
                }
            }
            il_free(selectedlist);
            if (pass_limit_reached) {
                onefield_clear_indexes(bp);
                solver_clear_indexes(sp);
                logverb("[engine-pass] state=end ordinal=%zu "
                        "reason=limit-during-materialization\n",
                        pass.ordinal);
                break;
            }

            logverb("Running solver:\n");
            onefield_log_run_parameters(bp);

            onefield_run(bp);

            if (bp->solver_failed) {
                rtn = -1;
                goto finish;
            }

            // we only want to try using the verify_wcses the first time.
            onefield_clear_verify_wcses(bp);
            onefield_clear_indexes(bp);
            onefield_clear_solutions(bp);
            solver_clear_indexes(sp);

            logverb("[engine-pass] state=end ordinal=%zu "
                    "solved=%i cancelled=%i "
                    "hit_total_cpu_limit=%i "
                    "hit_total_wall_limit=%i failed=%i\n",
                    pass.ordinal,
                    bp->single_field_solved ? 1 : 0,
                    bp->cancelled ? 1 : 0,
                    bp->hit_total_cpulimit ? 1 : 0,
                    bp->hit_total_timelimit ? 1 : 0,
                    bp->solver_failed ? 1 : 0);

            if (onefield_check_total_limits(bp)) {
                break;
            }
            if (onefield_is_run_obsolete(bp, sp)) {
                break;
            }
    }

    logverb("cx<=dx constraints: %i\n", sp->num_cxdx_skipped);
    logverb("meanx constraints: %i\n", sp->num_meanx_skipped);
    logverb("RA,Dec constraints: %i\n", sp->num_radec_skipped);
    logverb("AB scale constraints: %i\n", sp->num_abscale_skipped);

 finish:
   // SECTION INDEX-SHARD: engine-lifecycle
   if (index_shard_pool_started) {
     double pool_wall_start = monotonic_seconds();

     index_shard_pool_stop(bp);
     pool_stop_seconds =
         monotonic_seconds() - pool_wall_start;
   }
   if (residency) {
     (void)index_residency_quiesce(residency);
   }
   onefield_job_field_cache_end(bp);
   if (residency) {
     index_unbind_residency_service(residency);
   }

   logverb("[engine-profile] pool_start=%.6f pool_stop=%.6f "
           "engine_total=%.6f solver_failed=%i\n",
           pool_start_seconds,
           pool_stop_seconds,
           monotonic_seconds() - engine_wall_start,
           bp->solver_failed ? 1 : 0);

   solver_cleanup(sp);
   onefield_cleanup(bp);
   if (residency) {
     index_residency_stats_t stats;

     if (!index_residency_get_stats(residency, &stats)) {
       logverb(
           "[index-residency] copied_files=%llu copied_bytes=%llu "
           "hits=%llu deduplicated=%llu waits=%llu wait_ms=%.3f "
           "source_leases=%llu source_requeues=%llu "
           "cancellations=%llu "
           "peak_bytes=%zu ready_bytes=%zu live_handles=%zu "
           "failures=%llu source_changes=%llu\n",
           (unsigned long long)stats.files_copied,
           (unsigned long long)stats.bytes_copied,
           (unsigned long long)stats.cache_hits,
           (unsigned long long)stats.loading_deduplications,
           (unsigned long long)stats.wait_count,
           (double)stats.wait_nanoseconds / 1000000.0,
           (unsigned long long)stats.source_leases,
           (unsigned long long)stats.source_requeues,
           (unsigned long long)stats.cancelled_entries,
           stats.peak_resident_bytes,
           stats.ready_bytes,
           stats.live_handles,
           (unsigned long long)stats.copy_failures,
           (unsigned long long)stats.source_changes);
     }
     (void)index_residency_stop(residency);
   }
   return rtn;
}

static void parse_sip_coeffs(const qfits_header* hdr, const char* prefix, sip_t* wcs) {
    char key[64];
    int order, i, j;
    sprintf(key, "%sSAO", prefix);
    order = qfits_header_getint(hdr, key, -1);
    if (order >= 2) {
        if (order > 9)
            order = 9;
        wcs->a_order = order;
        wcs->b_order = order;
        for (i=0; i<=order; i++) {
            for (j=0; (i+j)<=order; j++) {
                if (i+j < 1)
                    continue;
                sprintf(key, "%sA%i%i", prefix, i, j);
                wcs->a[i][j] = qfits_header_getdouble(hdr, key, 0.0);
                sprintf(key, "%sB%i%i", prefix, i, j);
                wcs->b[i][j] = qfits_header_getdouble(hdr, key, 0.0);
            }
        }
    }
    sprintf(key, "%sSAPO", prefix);
    order = qfits_header_getint(hdr, key, -1);
    if (order >= 2) {
        if (order > 9)
            order = 9;
        wcs->ap_order = order;
        wcs->bp_order = order;
        for (i=0; i<=order; i++) {
            for (j=0; (i+j)<=order; j++) {
                if (i+j < 1)
                    continue;
                sprintf(key, "%sAP%i%i", prefix, i, j);
                wcs->ap[i][j] = qfits_header_getdouble(hdr, key, 0.0);
                sprintf(key, "%sBP%i%i", prefix, i, j);
                wcs->bp[i][j] = qfits_header_getdouble(hdr, key, 0.0);
            }
        }
    }
}

static anbool parse_job_from_qfits_header(const qfits_header* hdr, job_t* job) {
    onefield_t* bp = &(job->bp);
    solver_t* sp = &(bp->solver);

    double dnil = -LARGE_VAL;
    char *pstr;
    int n;
    anbool run;

    anbool default_tweak = TRUE;
    int default_tweakorder = 2;
    double default_odds_toprint = 1e6;
    double default_odds_tokeep = 1e9;
    double default_odds_tosolve = 1e9;
    double default_odds_totune = 1e6;
    //double default_image_fraction = 1.0;
    char* fn;
    double val;
    char pretty[FITS_LINESZ+1];

    onefield_init(bp);
    // must be in this order because init_parameters handily zeros out sp
    solver_set_default_values(sp);

    // Here we assume that the field's pixel coordinataes go from zero to IMAGEW,H.
    sp->field_maxx = qfits_header_getdouble(hdr, "IMAGEW", dnil);
    sp->field_maxy = qfits_header_getdouble(hdr, "IMAGEH", dnil);
    if ((sp->field_maxx == dnil) || (sp->field_maxy == dnil) ||
        (sp->field_maxx <= 0.0) || (sp->field_maxy <= 0.0)) {
        logerr("Must specify positive \"IMAGEW\" and \"IMAGEH\".\n");
        goto bailout;
    }

    sp->verify_uniformize = qfits_header_getboolean(hdr, "ANVERUNI", sp->verify_uniformize);
    sp->verify_dedup = qfits_header_getboolean(hdr, "ANVERDUP", sp->verify_dedup);

    val = qfits_header_getdouble(hdr, "ANPOSERR", 0.0);
    if (val > 0.0)
        sp->verify_pix = val;
    val = qfits_header_getdouble(hdr, "ANCTOL", 0.0);
    if (val > 0.0)
        sp->codetol = val;
    val = qfits_header_getdouble(hdr, "ANDISTR", 0.0);
    if (val > 0.0)
        sp->distractor_ratio = val;

    onefield_set_solvedout_file  (bp, fn=fits_get_long_string(hdr, "ANSOLVED"));
    free(fn);
    onefield_set_solvedin_file  (bp, fn=fits_get_long_string(hdr, "ANSOLVIN"));
    free(fn);
    onefield_set_match_file   (bp, fn=fits_get_long_string(hdr, "ANMATCH" ));
    free(fn);
    onefield_set_rdls_file    (bp, fn=fits_get_long_string(hdr, "ANRDLS"  ));
    free(fn);
    onefield_set_scamp_file   (bp, fn=fits_get_long_string(hdr, "ANSCAMP" ));
    free(fn);
    onefield_set_wcs_file     (bp, fn=fits_get_long_string(hdr, "ANWCS"   ));
    free(fn);
    onefield_set_corr_file    (bp, fn=fits_get_long_string(hdr, "ANCORR"  ));
    free(fn);
    onefield_set_cancel_file  (bp, fn=fits_get_long_string(hdr, "ANCANCEL"));
    free(fn);

    onefield_set_xcol(bp, fn=fits_get_dupstring(hdr, "ANXCOL"));
    free(fn);
    onefield_set_ycol(bp, fn=fits_get_dupstring(hdr, "ANYCOL"));
    free(fn);

    bp->timelimit = qfits_header_getdouble(hdr, "ANTLIM", 0.0);
    bp->cpulimit = qfits_header_getdouble(hdr, "ANCLIM", 0.0);
    if (qfits_header_getstr(hdr, "ANSHWRK")) {
        int requested_workers =
            qfits_header_getint(hdr, "ANSHWRK", INT_MIN);

        if (requested_workers == INT_MIN) {
            logerr("Invalid ANSHWRK worker value in augmented job header.\n");
            goto bailout;
        }

        job->index_shard_workers_override = requested_workers;
        job->index_shard_workers_override_set = TRUE;
    }
    bp->logratio_tosolve = log(qfits_header_getdouble(hdr, "ANODDSSL", default_odds_tosolve));
    logverb("Set odds ratio to solve to %g (log = %g)\n", exp(bp->logratio_tosolve), bp->logratio_tosolve);


    sp->logratio_toprint = log(qfits_header_getdouble(hdr, "ANODDSPR", default_odds_toprint));
    sp->logratio_tokeep = log(qfits_header_getdouble(hdr, "ANODDSKP", default_odds_tokeep));
    sp->logratio_totune = log(qfits_header_getdouble(hdr, "ANODDSTU", default_odds_totune));
    sp->logratio_bail_threshold = log(qfits_header_getdouble(hdr, "ANODDSBL", DEFAULT_BAIL_THRESHOLD));
    val = qfits_header_getdouble(hdr, "ANODDSST", 0.0);
    if (val > 0.0)
        sp->logratio_stoplooking = log(val);
    bp->best_hit_only = TRUE;

    // gotta keep it to solve it!
    sp->logratio_tokeep = MIN(sp->logratio_tokeep, bp->logratio_tosolve);
    // gotta print it to keep it (so what if that doesn't make sense)!
    sp->logratio_toprint = MIN(sp->logratio_toprint, sp->logratio_tokeep);

    // job->image_fraction = qfits_header_getdouble(hdr, "ANIMFRAC", job->image_fraction);
    job->include_default_scales = qfits_header_getboolean(hdr, "ANAPPDEF", 0);

    sp->parity = PARITY_BOTH;
    pstr = qfits_pretty_string_r(qfits_header_getstr(hdr, "ANPARITY"), pretty);
    if (pstr && streq(pstr, "NEG"))
        sp->parity = PARITY_FLIP;
    else if (pstr && streq(pstr, "POS"))
        sp->parity = PARITY_NORMAL;

    sp->set_crpix_center = qfits_header_getboolean(hdr, "ANCRPIXC", FALSE);
    sp->crpix[0] = qfits_header_getdouble(hdr, "ANCRPIX1", sp->crpix[0]);
    sp->crpix[1] = qfits_header_getdouble(hdr, "ANCRPIX2", sp->crpix[1]);
    sp->set_crpix = (sp->set_crpix_center ||
                     // were the values set?
                     qfits_header_getstr(hdr, "ANCRPIX1") ||
                     qfits_header_getstr(hdr, "ANCRPIX2"));

    if (qfits_header_getboolean(hdr, "ANTWEAK", default_tweak)) {
        int order = qfits_header_getint(hdr, "ANTWEAKO", default_tweakorder);
        //bp->do_tweak = TRUE;
        sp->do_tweak = TRUE;
        sp->tweak_aborder = order;
        sp->tweak_abporder = order;
    }

    if (!sp->do_tweak) {
        // No tweak: set tweak order to linear, because the tweak alg
        // can still be invoked via tune-up.
        sp->tweak_aborder = sp->tweak_abporder = 1;
    }

    val = qfits_header_getdouble(hdr, "ANQSFMIN", 0.0);
    if (val > 0.0)
        bp->quad_size_fraction_lo = val;
    val = qfits_header_getdouble(hdr, "ANQSFMAX", 0.0);
    if (val > 0.0)
        bp->quad_size_fraction_hi = val;

    job->ra_center = qfits_header_getdouble(hdr, "ANERA", dnil);
    job->dec_center = qfits_header_getdouble(hdr, "ANEDEC", dnil);
    job->search_radius = qfits_header_getdouble(hdr, "ANERAD", dnil);
    job->use_radec_center = ((job->ra_center     != dnil) &&
                             (job->dec_center    != dnil) &&
                             (job->search_radius != dnil));

    // tag-along columns
    bp->rdls_tagalong_all = qfits_header_getboolean(hdr, "ANTAGALL", FALSE);
    if (!bp->rdls_tagalong_all) {
        n = 1;
        while (1) {
            char key[64];
            char* val;
            sprintf(key, "ANTAG%i", n);
            val = fits_get_dupstring(hdr, key);
            if (!val)
                break;
            if (!bp->rdls_tagalong)
                bp->rdls_tagalong = sl_new(16);
            sl_append_nocopy(bp->rdls_tagalong, val);
            n++;
        }
    }

    // sort RDLS column
    bp->sort_rdls = fits_get_dupstring(hdr, "ANRDSORT");

    n = 1;
    while (1) {
        char key[64];
        double lo, hi;
        sprintf(key, "ANAPPL%i", n);
        lo = qfits_header_getdouble(hdr, key, 0.);
        sprintf(key, "ANAPPU%i", n);
        hi = qfits_header_getdouble(hdr, key, 0.);
        if ((hi == 0.) && (lo == 0.))
            break;
        if ((lo != 0.) && (hi != 0.)) {
            if ((lo < 0) || (lo > hi)) {
                logerr("Scale range %g to %g is invalid: min must be >= 0, max must be >= min.\n", lo, hi);
                goto bailout;
            }
        }
        dl_append(job->scales, lo);
        dl_append(job->scales, hi);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        int dlo, dhi;
        sprintf(key, "ANDPL%i", n);
        dlo = qfits_header_getint(hdr, key, 0);
        sprintf(key, "ANDPU%i", n);
        dhi = qfits_header_getint(hdr, key, 0);
        if (dlo == 0 && dhi == 0)
            break;
        if ((dlo < 1) || (dlo > dhi)) {
            logerr("Depth range %i to %i is invalid: min must be >= 1, max must be >= min.\n", dlo, dhi);
            goto bailout;
        }
        il_append(job->depths, dlo);
        il_append(job->depths, dhi);
        n++;
    }

    n = 1;
    while (1) {
        char lokey[64];
        char hikey[64];
        int lo, hi;
        sprintf(lokey, "ANFDL%i", n);
        lo = qfits_header_getint(hdr, lokey, -1);
        if (lo == -1)
            break;
        sprintf(hikey, "ANFDU%i", n);
        hi = qfits_header_getint(hdr, hikey, -1);
        if (hi == -1)
            break;
        if ((lo <= 0) || (lo > hi)) {
            char pretty1[FITS_LINESZ+1];
            char pretty2[FITS_LINESZ+1];
            logerr("Field range %i to %i is invalid: min must be >= 1, max must be >= min.\n", lo, hi);
            qfits_pretty_string_r(qfits_header_getstr(hdr, lokey), pretty1);
            qfits_pretty_string_r(qfits_header_getstr(hdr, hikey), pretty2);
            logmsg("  (FITS headers: \"%s = %s\", \"%s = %s\")\n",
                   lokey, pretty1, hikey, pretty2);
            goto bailout;
        }

        onefield_add_field_range(bp, lo, hi);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        int fld;
        sprintf(key, "ANFD%i", n);
        fld = qfits_header_getint(hdr, key, -1);
        if (fld == -1)
            break;
        if (fld <= 0) {
            qfits_pretty_string_r(qfits_header_getstr(hdr, key), pretty);
            logerr("Field %i is invalid: must be >= 1.  (FITS header: \"%s = %s\")\n", fld, key, pretty);
            goto bailout;
        }

        onefield_add_field(bp, fld);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        sip_t wcs;
        char* keys[] = { "ANW%iPIX1", "ANW%iPIX2", "ANW%iVAL1", "ANW%iVAL2",
                         "ANW%iCD11", "ANW%iCD12", "ANW%iCD21", "ANW%iCD22" };
        double* vals[] = { &(wcs.wcstan. crval[0]), &(wcs.wcstan.crval[1]),
                           &(wcs.wcstan.crpix[0]), &(wcs.wcstan.crpix[1]),
                           &(wcs.wcstan.cd[0][0]), &(wcs.wcstan.cd[0][1]),
                           &(wcs.wcstan.cd[1][0]), &(wcs.wcstan.cd[1][1]) };
        int j;
        int bail = 0;
        memset(&wcs, 0, sizeof(wcs));
        for (j = 0; j < 8; j++) {
            sprintf(key, keys[j], n);
            *(vals[j]) = qfits_header_getdouble(hdr, key, dnil);
            if (*(vals[j]) == dnil) {
                bail = 1;
                break;
            }
        }
        if (bail)
            break;

        // SIP terms
        sprintf(key, "ANW%i", n);
        parse_sip_coeffs(hdr, key, &wcs);

        sip_ensure_inverse_polynomials(&wcs);

        onefield_add_verify_wcs(bp, &wcs);
        n++;
    }

    // Distortion to apply before matching...
    do {
        sip_t dsip;
        double p0, p1;
        memset(&dsip, 0, sizeof(sip_t));
        p0 = qfits_header_getdouble(hdr, "ANDPIX0", dnil);
        if (p0 == dnil)
            break;
        p1 = qfits_header_getdouble(hdr, "ANDPIX1", dnil);
        if (p1 == dnil)
            break;
        dsip.wcstan.crpix[0] = p0;
        dsip.wcstan.crpix[1] = p1;
        parse_sip_coeffs(hdr, "AND", &dsip);
        if ((dsip.a_order > 1 && dsip.b_order > 1) ||
            (dsip.ap_order > 1 && dsip.bp_order > 1)) {
            sp->predistort = malloc(sizeof(sip_t));
            memcpy(sp->predistort, &dsip, sizeof(sip_t));
        }
    } while (0);

    sp->pixel_xscale = qfits_header_getdouble(hdr, "ANPXSCAL", 0.);

    run = qfits_header_getboolean(hdr, "ANRUN", FALSE);

    // Default: solve first field.
    if (run && !il_size(bp->fieldlist)) {
        onefield_add_field(bp, 1);
    }

    return TRUE;

 bailout:
    return FALSE;
}



engine_t* engine_new() {
    engine_t* engine = calloc(1, sizeof(engine_t));
    engine->index_paths = sl_new(10);
    engine->indexes = pl_new(16);
    engine->free_indexes = pl_new(16);
    engine->free_mindexes = pl_new(16);
    engine->ismallest = il_new(4);
    engine->ibiggest = il_new(4);
    engine->default_depths = il_new(4);
    engine->sizesmallest = LARGE_VAL;
    engine->sizebiggest = -LARGE_VAL;

    // Default scale estimate: field width, in degrees:
    engine->minwidth = 0.1;
    engine->maxwidth = 180.0;
    engine->cpulimit = 600.0;
    engine->index_shard_workers_config = INDEX_SHARD_WORKERS_AUTO;
    return engine;
}

void engine_free(engine_t* engine) {
    int i;
    if (!engine)
        return;
    if (engine->free_indexes) {
        for (i=0; i<pl_size(engine->free_indexes); i++) {
            index_t* ind = pl_get(engine->free_indexes, i);
            index_free(ind);
        }
        pl_free(engine->free_indexes);
    }
    if (engine->free_mindexes) {
        for (i=0; i<pl_size(engine->free_mindexes); i++) {
            multiindex_t* mi = pl_get(engine->free_mindexes, i);
            multiindex_free(mi);
        }
        pl_free(engine->free_mindexes);
    }
    pl_free(engine->indexes);
    if (engine->ismallest)
        il_free(engine->ismallest);
    if (engine->ibiggest)
        il_free(engine->ibiggest);
    if (engine->default_depths)
        il_free(engine->default_depths);
    if (engine->index_paths)
        sl_free2(engine->index_paths);
    free(engine);
}

static int engine_resolve_index_shard_workers(engine_t *engine,
                                              job_t *job) {
    const char *environment_value;
    const char *source;
    int available_cpus;
    int requested_workers;
    int resolved_workers;
    char requested_text[32];

    if (!engine || !job) {
        ERROR("Cannot resolve parallel workers without engine and job state");
        return -1;
    }

    available_cpus = index_shard_config_available_cpus();
    requested_workers = engine->index_shard_workers_config;
    source = engine->index_shard_workers_config_set
        ? "config"
        : "built-in";

    if (job->index_shard_workers_override_set) {
        requested_workers = job->index_shard_workers_override;
        if (index_shard_config_validate_workers(requested_workers,
                                                available_cpus)) {
            ERROR("Invalid ANSHWRK worker override %i: "
                  "expected automatic selection or an integer from "
                  "1 through %i",
                  requested_workers,
                  available_cpus);
            return -1;
        }
        source = "solve-field";
    } else {
        /*
         * Retain the legacy environment override for existing measurement
         * harnesses. New production commands should use the per-job
         * solve-field option, whose AXY header has higher precedence.
         */
        environment_value = getenv("ASTROMETRY_P_WORKERS");
        if (!environment_value || !environment_value[0]) {
            environment_value = getenv("ASTROMETRY_INDEX_SHARD_WORKERS");
        }
        if (environment_value && environment_value[0]) {
            if (index_shard_config_parse_workers(environment_value,
                                                 available_cpus,
                                                 &requested_workers)) {
                ERROR("Invalid parallel worker environment value \"%s\": "
                      "expected \"auto\" or an integer from 1 through %i",
                      environment_value,
                      available_cpus);
                return -1;
            }
            source = "environment";
        }
    }

    resolved_workers =
        index_shard_config_resolve_workers(requested_workers,
                                           available_cpus);
    if (resolved_workers < 1) {
        ERROR("Failed to resolve parallel worker count");
        return -1;
    }

    job->bp.index_shard_workers = resolved_workers;
    job->index_shard_workers_controlled =
        strcmp(source, "built-in") != 0;

    if (requested_workers == INDEX_SHARD_WORKERS_AUTO) {
        snprintf(requested_text, sizeof(requested_text), "auto");
    } else {
        snprintf(requested_text,
                 sizeof(requested_text),
                 "%i",
                 requested_workers);
    }

    logverb("[parallel] worker-config source=%s requested=%s "
            "available=%i effective=%i mode=%s\n",
            source,
            requested_text,
            available_cpus,
            resolved_workers,
            resolved_workers > 1 ? "pthread" : "serial");

    return 0;
}

job_t* engine_read_job_file(engine_t* engine, const char* jobfn) {
    qfits_header* hdr;
    job_t* job;
    onefield_t* bp;

    // Read primary header.
    hdr = anqfits_get_header2(jobfn, 0);
    if (!hdr) {
        ERROR("Failed to parse FITS header from file \"%s\"", jobfn);
        return NULL;
    }
    job = job_new();
    if (!parse_job_from_qfits_header(hdr, job)) {
        job_free(job);
        qfits_header_destroy(hdr);
        return NULL;
    }
    qfits_header_destroy(hdr);

    bp = &(job->bp);

    if (engine_resolve_index_shard_workers(engine, job)) {
        solver_cleanup(&bp->solver);
        onefield_cleanup(bp);
        job_free(job);
        return NULL;
    }

    onefield_set_field_file(bp, jobfn);

    // If the job has no scale estimate, search everything provided
    // by the engine
    if (!dl_size(job->scales) || job->include_default_scales) {
        double arcsecperpix;
        arcsecperpix = deg2arcsec(engine->minwidth) / job_imagew(job);
        dl_append(job->scales, arcsecperpix);
        arcsecperpix = deg2arcsec(engine->maxwidth) / job_imagew(job);
        dl_append(job->scales, arcsecperpix);
    }

    // The job can only decrease the CPU limit.
    // SECTION INDEX-SHARD: cpu-limit-precedence
    /*
     * Upstream-compatible CPU-limit handling.
     *
     * Sources:
     *   - bp->cpulimit      job/ANCLIM limit, normally produced by solve-field
     *                       --cpulimit
     *   - engine->cpulimit  backend/config limit from astrometry.cfg
     *
     * Semantics:
     *   - if both exist, use the smaller one
     *   - if only one exists, use that one
     *   - if neither exists, run without a CPU limit
     *
     * This preserves the original astrometry.net behavior: solve-field may
     * reduce a backend/config CPU limit but must not increase it.
     *
     * pthread index-sharding detail:
     *   In pthread mode the effective limit is stored in bp->total_cpulimit and
     *   bp->cpulimit is cleared.  This prevents worker-local onefield copies
     * from treating the same budget as an independent per-index/per-worker
     * limit.
     */
    {
      double job_cpulimit = bp->cpulimit;
      double cfg_cpulimit = engine->cpulimit;
      double effective_cpulimit = 0.0;

      if (job_cpulimit > 0.0 && cfg_cpulimit > 0.0) {
        effective_cpulimit =
            (job_cpulimit < cfg_cpulimit) ? job_cpulimit : cfg_cpulimit;
      } else if (job_cpulimit > 0.0) {
        effective_cpulimit = job_cpulimit;
      } else if (cfg_cpulimit > 0.0) {
        effective_cpulimit = cfg_cpulimit;
      }

      if (effective_cpulimit > 0.0) {
        logverb("Using effective CPU time limit of %g seconds "
                "(job=%g, config=%g)\n",
                effective_cpulimit, job_cpulimit, cfg_cpulimit);
      } else {
        logverb("No CPU time limit set for this job "
                "(job=%g, config=%g)\n",
                job_cpulimit, cfg_cpulimit);
      }

      if (index_shard_pthread_enabled(bp)) {
        /*
         * pthread path:
         * total_cpulimit is the process-wide budget checked by
         * index_shard_check_global_cpu_limit().
         */
        bp->total_cpulimit = effective_cpulimit;
        bp->cpulimit = 0.0;
      } else {
        /*
         * Original/non-pthread path:
         * keep bp->cpulimit as the ordinary effective run limit.
         */
        bp->cpulimit = effective_cpulimit;
        bp->total_cpulimit = effective_cpulimit;
      }

      bp->total_timelimit = bp->timelimit;
    }

    logverb("[index-shard] engine limits after setup: "
            "cpulimit=%f total_cpulimit=%f timelimit=%g total_timelimit=%g\n",
            bp->cpulimit, bp->total_cpulimit, bp->timelimit,
            bp->total_timelimit);

    // If the job didn't specify depths, set defaults.
    if (il_size(job->depths) == 0) {
        if (il_size(engine->default_depths) != 0) {
            il_append_list(job->depths, engine->default_depths);
        } else {
            /*
             * An empty site default means the original unbounded depth
             * interval. Keep this scientific search space independent of
             * worker count and of the legacy "inparallel" token.
             */
            il_append(job->depths, 0);
            il_append(job->depths, 0);
        }
    }

    if (engine->cancelfn)
        onefield_set_cancel_file(bp, engine->cancelfn);
    if (engine->solvedfn)
        onefield_set_solved_file(bp, engine->solvedfn);

    return job;
}

void job_set_cancel_file(job_t* job, const char* fn) {
    onefield_set_cancel_file(&(job->bp), fn);
}

void job_set_solved_file(job_t* job, const char* fn) {
    onefield_set_solved_file(&(job->bp), fn);
}

// Modify all filenames to be relative to "dir".
int job_set_base_dir(job_t* job, const char* dir) {
    return job_set_output_base_dir(job, dir) ||
        job_set_input_base_dir(job, dir);
}

int job_set_input_base_dir(job_t* job, const char* dir) {
    char* path;
    onefield_t* bp = &(job->bp);
    logverb("Changing input file base dir to %s\n", dir);
    if (bp->fieldfname) {
        path = resolve_path(bp->fieldfname, dir);
        logverb("Changing %s to %s\n", bp->fieldfname, path);
        onefield_set_field_file(bp, path);
    }
    return 0;
}

int job_set_output_base_dir(job_t* job, const char* dir) {
    char* path;
    onefield_t* bp = &(job->bp);
    logverb("Changing output file base dir to %s\n", dir);
    if (bp->cancelfname) {
        path = resolve_path(bp->cancelfname, dir);
        logverb("Cancel file was %s, changing to %s.\n", bp->cancelfname, path);
        onefield_set_cancel_file(bp, path);
    }
    if (bp->solved_in) {
        path = resolve_path(bp->solved_in, dir);
        logverb("Changing %s to %s\n", bp->solved_in, path);
        onefield_set_solvedin_file(bp, path);
    }
    if (bp->solved_out) {
        path = resolve_path(bp->solved_out, dir);
        logverb("Changing %s to %s\n", bp->solved_out, path);
        onefield_set_solvedout_file(bp, path);
    }
    if (bp->matchfname) {
        path = resolve_path(bp->matchfname, dir);
        logverb("Changing %s to %s\n", bp->matchfname, path);
        onefield_set_match_file(bp, path);
    }
    if (bp->indexrdlsfname) {
        path = resolve_path(bp->indexrdlsfname, dir);
        logverb("Changing %s to %s\n", bp->indexrdlsfname, path);
        onefield_set_rdls_file(bp, path);
    }
    if (bp->scamp_fname) {
        path = resolve_path(bp->scamp_fname, dir);
        logverb("Changing %s to %s\n", bp->scamp_fname, path);
        onefield_set_scamp_file(bp, path);
    }
    if (bp->corr_fname) {
        path = resolve_path(bp->corr_fname, dir);
        logverb("Changing %s to %s\n", bp->corr_fname, path);
        onefield_set_corr_file(bp, path);
    }
    if (bp->wcs_template) {
        path = resolve_path(bp->wcs_template, dir);
        logverb("Changing %s to %s\n", bp->wcs_template, path);
        onefield_set_wcs_file(bp, path);
    }
    return 0;
}
