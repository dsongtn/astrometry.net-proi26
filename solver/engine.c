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
#include "engine_private.h"

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
        } else if (is_word(line, "walllimit ", &nextword) ||
                   is_word(line, "wall_limit ", &nextword)) {
            engine->walllimit = atof(nextword);
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

    app_min_default = deg2arcsec(engine->minwidth) /
        engine_job_imagew(job);
    app_max_default = deg2arcsec(engine->maxwidth) /
        engine_job_imagew(job);

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
                MIN(engine_job_imagew(job), engine_job_imageh(job));

            // range of quad sizes that could be found in the field,
            // in arcsec.
            // the hypotenuse...
            fmax = bp->quad_size_fraction_hi *
                hypot(engine_job_imagew(job), engine_job_imageh(job)) *
                app_max;
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
    engine->walllimit = 300.0;
    engine->cpulimit = 0.0;
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
