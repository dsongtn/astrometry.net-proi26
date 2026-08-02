/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

/**
 * Solve a single field
 *
 * Inputs: .ckdt .quad .skdt
 * Output: .match .rdls .wcs, ...
 */

#include <sys/types.h>
#include <sys/time.h>
#include <libgen.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "anqfits.h"
#include "index_shard_config.h"
#include "index_shard_internal.h"
#include "bl-sort.h"
#include "boilerplate.h"
#include "codekd.h"
#include "errors.h"
#include "fitsbin.h"
#include "fitsioutils.h"
#include "index.h"
#include "log.h"
#include "mathutil.h"
#include "onefield.h"
#include "onefield_internal.h"
#include "os-features.h"
#include "permutedsort.h"
#include "quadfile.h"
#include "scamp-catalog.h"
#include "sip_qfits.h"
#include "solvedfile.h"
#include "starkd.h"
#include "starutil.h"
#include "tic.h"
#include "tweak.h"
#include "tweak2.h"
#include "verify.h"

static void add_onefield_params(onefield_t* bp, qfits_header* hdr);
static void load_and_parse_wcsfiles(onefield_t* bp);
static void solve_fields(onefield_t* bp, sip_t* verify_wcs);
static anbool is_field_solved(onefield_t* bp, int fieldnum);
static int write_solutions(onefield_t* bp);
static int publish_solved_fields(onefield_t* bp);
static void remove_duplicate_solutions(onefield_t* bp);
// A tag-along column for index rdls / correspondence file.
struct tagalong {
    tfits_type type;
    int arraysize;
    char* name;
    char* units;
    void* data;
    // size in bytes of one item.
    int itemsize;
    int Ndata;
    // assigned by rdlist_add_tagalong_column
    int colnum;
};
typedef struct tagalong tagalong_t;

static anbool grab_tagalong_data(startree_t* starkd, MatchObj* mo, onefield_t* bp,
                                 const int* starinds, int N) {
    fitstable_t* tagalong;
    int i;
    tagalong = startree_get_tagalong(starkd);
    if (!tagalong) {
        ERROR("Failed to find tag-along table in index");
        return FALSE;
    }
    if (!mo->tagalong)
        mo->tagalong = bl_new(16, sizeof(tagalong_t));

    if (bp->rdls_tagalong_all) { // && ! bp->done_rdls_tagalong_all
        char* cols;
        // retrieve all column names.
        bp->rdls_tagalong = fitstable_get_fits_column_names(tagalong, bp->rdls_tagalong);
        cols = sl_join(bp->rdls_tagalong, ", ");
        logverb("Found tag-along columns: %s\n", cols);
        free(cols);
        //
        sl_remove_duplicates(bp->rdls_tagalong);
        cols = sl_join(bp->rdls_tagalong, ", ");
        logverb("After removing duplicates: %s\n", cols);
        free(cols);
    }
    for (i=0; i<sl_size(bp->rdls_tagalong); i++) {
        const char* col = sl_get(bp->rdls_tagalong, i);
        tagalong_t tag;
        if (fitstable_find_fits_column(tagalong, col, &(tag.units), &(tag.type), &(tag.arraysize))) {
            ERROR("Failed to find column \"%s\" in index", col);
            continue;
        }
        tag.data = fitstable_read_column_array_inds(tagalong, col, tag.type, starinds, N, NULL);
        if (!tag.data) {
            ERROR("Failed to read data for column \"%s\" in index", col);
            continue;
        }
        if ((strcaseeq(col, "ra") || strcaseeq(col, "dec")))
            asprintf_safe(&(tag.name), "%s_ref", col);
        else
            tag.name = strdup(col);
        tag.units = strdup(tag.units);
        tag.itemsize = fits_get_atom_size(tag.type) * tag.arraysize;
        tag.Ndata = N;
        bl_append(mo->tagalong, &tag);
    }
    return TRUE;
}

static anbool grab_field_tagalong_data(MatchObj* mo, xylist_t* xy, int N) {
    fitstable_t* tagalong;
    int i;
    sl* lst;
    if (!mo->field_tagalong)
        mo->field_tagalong = bl_new(16, sizeof(tagalong_t));
    tagalong = xy->table;
    lst = xylist_get_tagalong_column_names(xy, NULL);
    {
        char* txt = sl_join(lst, " ");
        logverb("Found tag-along columns from field: %s\n", txt);
        free(txt);
    }
    for (i=0; i<sl_size(lst); i++) {
        const char* col = sl_get(lst, i);
        tagalong_t tag;
        if (fitstable_find_fits_column(tagalong, col, &(tag.units), &(tag.type), &(tag.arraysize))) {
            ERROR("Failed to find column \"%s\" in index", col);
            continue;
        }
        tag.data = fitstable_read_column_array(tagalong, col, tag.type);
        if (!tag.data) {
            ERROR("Failed to read data for column \"%s\" in index", col);
            continue;
        }
        tag.name = strdup(col);
        tag.units = strdup(tag.units);
        tag.itemsize = fits_get_atom_size(tag.type) * tag.arraysize;
        tag.Ndata = N;
        bl_append(mo->field_tagalong, &tag);
    }
    sl_free2(lst);
    return TRUE;
}


void onefield_clear_verify_wcses(onefield_t* bp) {
    bl_remove_all(bp->verify_wcs_list);
}

void onefield_clear_solutions(onefield_t* bp) {
    bl_remove_all(bp->solutions);
    il_remove_all(bp->solved_fields_pending);
}

void onefield_clear_indexes(onefield_t* bp) {
    int i;

    for (i = 0; i < pl_size(bp->owned_indexes); i++) {
        index_free(pl_get(bp->owned_indexes, i));
    }
    pl_remove_all(bp->owned_indexes);
    sl_remove_all(bp->indexnames);
    pl_remove_all(bp->indexes);
}

void onefield_set_field_file(onefield_t* bp, const char* fn) {
    onefield_job_field_cache_invalidate(bp);
    free(bp->fieldfname);
    bp->fieldfname = strdup_safe(fn);
}

void onefield_set_solved_file(onefield_t* bp, const char* fn) {
    onefield_set_solvedin_file (bp, fn);
    onefield_set_solvedout_file(bp, fn);
}

void onefield_set_solvedin_file(onefield_t* bp, const char* fn) {
    free(bp->solved_in);
    bp->solved_in = strdup_safe(fn);
}

void onefield_set_solvedout_file(onefield_t* bp, const char* fn) {
    free(bp->solved_out);
    bp->solved_out = strdup_safe(fn);
}

void onefield_set_cancel_file(onefield_t* bp, const char* fn) {
    free(bp->cancelfname);
    bp->cancelfname = strdup_safe(fn);
}

void onefield_set_match_file(onefield_t* bp, const char* fn) {
    free(bp->matchfname);
    bp->matchfname = strdup_safe(fn);
}

void onefield_set_rdls_file(onefield_t* bp, const char* fn) {
    free(bp->indexrdlsfname);
    bp->indexrdlsfname = strdup_safe(fn);
}

void onefield_set_scamp_file(onefield_t* bp, const char* fn) {
    free(bp->scamp_fname);
    bp->scamp_fname = strdup_safe(fn);
}

void onefield_set_corr_file(onefield_t* bp, const char* fn) {
    free(bp->corr_fname);
    bp->corr_fname = strdup_safe(fn);
}

void onefield_set_wcs_file(onefield_t* bp, const char* fn) {
    free(bp->wcs_template);
    bp->wcs_template = strdup_safe(fn);
}

void onefield_set_xcol(onefield_t* bp, const char* x) {
    onefield_job_field_cache_invalidate(bp);
    free(bp->xcolname);
    if (!x)
        x = "X";
    bp->xcolname = strdup(x);
}

void onefield_set_ycol(onefield_t* bp, const char* y) {
    onefield_job_field_cache_invalidate(bp);
    free(bp->ycolname);
    if (!y)
        y = "Y";
    bp->ycolname = strdup_safe(y);
}

void onefield_add_index(onefield_t* bp, const char* index) {
    sl_append(bp->indexnames, index);
}

void onefield_add_loaded_index(onefield_t* bp, index_t* ind) {
    pl_append(bp->indexes, ind);
}

void onefield_add_owned_index(onefield_t* bp, index_t* ind) {
    pl_append(bp->indexes, ind);
    pl_append(bp->owned_indexes, ind);
}

void onefield_add_verify_wcs(onefield_t* bp, sip_t* wcs) {
    bl_append(bp->verify_wcs_list, wcs);
}

void onefield_add_field(onefield_t* bp, int field) {
    onefield_job_field_cache_invalidate(bp);
    il_insert_unique_ascending(bp->fieldlist, field);
}

void onefield_add_field_range(onefield_t* bp, int lo, int hi) {
    int i;
    onefield_job_field_cache_invalidate(bp);
    for (i=lo; i<=hi; i++) {
        il_insert_unique_ascending(bp->fieldlist, i);
    }
}

anbool onefield_check_total_limits(onefield_t* bp) {
    if (!bp) {
        return TRUE;
    }

    if (bp->total_timelimit > 0.0 && !bp->hit_total_timelimit) {
        double now = monotonic_seconds();
        if (now >= 0.0 &&
            now - bp->time_total_start >= bp->total_timelimit) {
            logmsg("Total wall-clock time limit reached!\n");
            bp->hit_total_timelimit = TRUE;
        }
    }

    if (bp->total_cpulimit > 0.0 && !bp->hit_total_cpulimit) {
        float now = get_cpu_usage();
        if (now - bp->cpu_total_start >= bp->total_cpulimit) {
            logmsg("Total CPU time limit reached!\n");
            bp->hit_total_cpulimit = TRUE;
        }
    }

    if (bp->hit_total_timelimit || bp->hit_total_cpulimit) {
        bp->solver.quit_now = TRUE;
        return TRUE;
    }

    return FALSE;
}

static void check_time_limits(onefield_t* bp) {
    onefield_check_total_limits(bp);

    if (bp->timelimit > 0.0 && !bp->hit_timelimit) {
        double now = monotonic_seconds();
        if (now >= 0.0 && now - bp->time_start >= bp->timelimit) {
            logmsg("Wall-clock time limit reached!\n");
            bp->hit_timelimit = TRUE;
        }
    }

    if (bp->cpulimit > 0.0 && !bp->hit_cpulimit) {
        float now = get_cpu_usage();
        if (now - bp->cpu_start > bp->cpulimit) {
            logmsg("CPU time limit reached!\n");
            bp->hit_cpulimit = TRUE;
        }
    }

    if (bp->hit_total_timelimit ||
        bp->hit_total_cpulimit ||
        bp->hit_timelimit ||
        bp->hit_cpulimit) {
        bp->solver.quit_now = TRUE;
    }
}
void onefield_run(onefield_t* bp) {
    solver_t* sp = &(bp->solver);
    size_t i, I;
    size_t Nindexes = 0U;
    size_t profile_indexes_executed = 0;
    double profile_wall_start = monotonic_seconds();
    double profile_acquire_seconds = 0.0;
    double profile_solver_seconds = 0.0;
    double profile_release_seconds = 0.0;
    double profile_output_seconds = 0.0;
    const char* profile_mode = "serial";
    anbool verification_datalog;
    anbool job_field_prepared = FALSE;
    anbool shard_candidate;

    /*
     * engine_run_job() initializes total-job limits once. Direct onefield
     * callers arrive with zeroed start values and initialize them here.
     */
    if (bp->time_total_start <= 0.0) {
        bp->time_total_start = monotonic_seconds();
    }
    if (bp->cpu_total_start <= 0.0) {
        bp->cpu_total_start = get_cpu_usage();
    }
    if (onefield_check_total_limits(bp)) {
        goto cleanup;
    }

    // Parse WCS files submitted for verification.
    load_and_parse_wcsfiles(bp);

    /*
     * A single-field cached job validates field existence inside the stable
     * identity acquisition. Opening here would create an unguarded metadata
     * epoch before that acquisition and could remove a field based on a
     * pathname that is replaced immediately afterward. Multi-field and
     * uncached jobs retain the original up-front enumeration.
     */
    if (bp->job_field_cache &&
        il_size(bp->fieldlist) == 1) {
        if (onefield_internal_validate_single_field_list(bp)) {
            bp->solver_failed = TRUE;
            goto cleanup;
        }
    } else {
        if (onefield_internal_open_master_xyls(bp)) {
            exit(-1);
        }
        onefield_internal_remove_invalid_fields(
            bp->fieldlist,
            xylist_n_fields(bp->xyls));
    }

    Nindexes = onefield_internal_index_count(bp);
    verification_datalog = verify_datalog_enabled();
    if (onefield_check_total_limits(bp)) {
        goto cleanup;
    }

    // Verify any WCS estimates we have.
    if (bl_size(bp->verify_wcs_list)) {
        int i;
        int w;

        // We want to get the best logodds out of all the indices, so we set the
        // logodds-to-solve impossibly high so that a "good enough" solution doesn't
        // stop us from continuing to search...
        double oldodds = bp->logratio_tosolve;
        bp->logratio_tosolve = LARGE_VAL;

        for (w = 0; w < bl_size(bp->verify_wcs_list); w++) {
            double pixscale;
            double quadlo, quadhi;
            sip_t* wcs = bl_access(bp->verify_wcs_list, w);

            // We don't want to try to verify a wide-field image using a narrow-
            // field index, because it will contain a TON of index stars in the
            // field.  We therefore only try to verify using indices that contain
            // quads that could have been found in the image.
            if (wcs->wcstan.imagew == 0.0 && sp->field_maxx > 0.0)
                wcs->wcstan.imagew = sp->field_maxx;
            if (wcs->wcstan.imageh == 0.0 && sp->field_maxy > 0.0)
                wcs->wcstan.imageh = sp->field_maxy;

            if ((wcs->wcstan.imagew == 0) ||
                (wcs->wcstan.imageh == 0)) {
                logmsg("Verifying WCS: image width or height is zero / unknown.\n");
                continue;
            }
            pixscale = sip_pixel_scale(wcs);
            quadlo = bp->quad_size_fraction_lo
                * MIN(wcs->wcstan.imagew, wcs->wcstan.imageh)
                * pixscale;
            quadhi = bp->quad_size_fraction_hi
                * MAX(wcs->wcstan.imagew, wcs->wcstan.imageh)
                * pixscale;
            logmsg("Verifying WCS using indices with quads of size [%g, %g] arcmin\n",
                   arcsec2arcmin(quadlo), arcsec2arcmin(quadhi));

            for (I=0; I<Nindexes; I++) {
                index_t* index = onefield_internal_get_index(bp, I);
                if (!index_overlaps_scale_range(index, quadlo, quadhi)) {
                    if (onefield_internal_done_with_index(bp, I, index)) {
                        bp->solver_failed = TRUE;
                        break;
                    }
                    continue;
                }
                solver_add_index(sp, index);
                sp->index = index;
                logmsg("Verifying WCS with index %zu of %zu (%s)\n",  I + 1, Nindexes, index->indexname);
                // Do it!
                solve_fields(bp, wcs);
                // Clean up this index...
                if (onefield_internal_done_with_index(bp, I, index)) {
                    bp->solver_failed = TRUE;
                }
                solver_clear_indexes(sp);

                if (bp->solver_failed) {
                    break;
                }
            }

            if (bp->solver_failed) {
                break;
            }
        }

        bp->logratio_tosolve = oldodds;

        if (bp->solver_failed) {
            goto cleanup;
        }

        logmsg("Got %zu solutions.\n", bl_size(bp->solutions));

        if (bp->best_hit_only)
            remove_duplicate_solutions(bp);

        for (i=0; i<bl_size(bp->solutions); i++) {
            MatchObj* mo = bl_access(bp->solutions, i);
            if (mo->logodds >= bp->logratio_tosolve)
                onefield_internal_solved_field(bp, mo->fieldnum);
        }
    }

    if (bp->single_field_solved)
        goto cleanup;

    /*
     * A nonzero lower bound reconstructs the same field prefix for every
     * selected index. Prepare that immutable geometry once whenever there is
     * cross-index reuse, independent of whether this run uses one worker or
     * the pthread pool. Scheduler topology only determines who borrows it.
     */
    if (sp->startobj > 0 &&
        Nindexes > 1U &&
        il_size(bp->fieldlist) == 1) {
      int prepare_status =
          onefield_index_shard_prepare_job_field_for_run(bp);

      if (prepare_status < 0) {
        bp->solver_failed = TRUE;
        goto cleanup;
      }
      job_field_prepared = (prepare_status == 0);
      if (onefield_check_total_limits(bp)) {
        goto cleanup;
      }
    }

    // SECTION INDEX-SHARD: onefield-entry
    /*
     * The outer reducer currently commits the first valid field and stops the
     * pass. That is exact for the production single-field, first-solution
     * workload, but not for multi-extension XYLS input where other fields
     * must continue or for nsolves > 1 when hits can span indexes. Keep those
     * runs on the legacy serial path until the reducer tracks the missing
     * completion state authoritatively.
     */
    shard_candidate =
        index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        Nindexes > 0 &&
        pl_size(bp->indexes) == 0 &&
        il_size(bp->fieldlist) == 1 &&
        bp->nsolves <= 1 &&
        sp->maxquads == 0 &&
        sp->maxmatches == 0 &&
        !bp->rdls_tagalong_all &&
        !verification_datalog;
    if (shard_candidate && !job_field_prepared) {
      int prepare_status =
          onefield_index_shard_prepare_job_field_for_run(bp);

      if (prepare_status < 0) {
        bp->solver_failed = TRUE;
        goto cleanup;
      }
      if (prepare_status > 0) {
        logmsg("[index-shard] pthread path unavailable without "
               "a stable shared-field identity; using original "
               "serial path\n");
        profile_mode = "serial-unavailable-field-identity";
        shard_candidate = FALSE;
      } else {
        job_field_prepared = TRUE;
      }
    }
    if (shard_candidate) {
      index_shard_solve_status_t shard_status;
      double shard_wall_start = monotonic_seconds();
      if (onefield_check_total_limits(bp)) {
        goto cleanup;
      }

      profile_mode = "pthread-data-ready";
      shard_status =
          onefield_index_shard_solve(
              bp,
              sp,
              Nindexes);
      profile_solver_seconds +=
          monotonic_seconds() - shard_wall_start;
      onefield_internal_job_index_cache_flush(bp);

      switch (shard_status) {
      case INDEX_SHARD_SOLVE_HANDLED:
        goto cleanup;

      case INDEX_SHARD_SOLVE_UNAVAILABLE:
        logmsg("[index-shard] pthread path unavailable; using original "
               "serial path\n");
        profile_mode = "serial-unavailable";
        break;

      case INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE:
        logmsg("[index-shard] pthread solve failed before master commit; "
               "using original serial path\n");
        profile_mode = "serial-precommit-retry";
        break;

      case INDEX_SHARD_SOLVE_TERMINAL_FAILURE:
        logerr("[index-shard] pthread solve hit a global-integrity or "
               "post-commit failure; serial fallback suppressed\n");
        bp->solver_failed = TRUE;
        goto cleanup;

      case INDEX_SHARD_SOLVE_LIFECYCLE_CONFLICT:
        logerr("[index-shard] pthread lifecycle conflict; "
               "serial fallback suppressed\n");
        bp->solver_failed = TRUE;
        goto cleanup;

      default:
        logerr("[index-shard] unexpected solve status %i; "
               "serial fallback suppressed\n",
               (int)shard_status);
        bp->solver_failed = TRUE;
        goto cleanup;
      }
    }

    if (index_shard_pthread_enabled(bp) &&
        pl_size(bp->indexes) != 0) {
        logverb("[index-shard] loaded multiindex component uses exact "
                "serial path\n");
        profile_mode = "serial-loaded-index";
    }

    if (index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        il_size(bp->fieldlist) != 1) {
        logverb("[index-shard] multi-field input uses exact serial path "
                "until field-aware reduction is available\n");
        profile_mode = "serial-multi-field";
    }

    if (index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        il_size(bp->fieldlist) == 1 &&
        bp->nsolves > 1) {
        logverb("[index-shard] nsolves=%i uses exact serial path until "
                "cross-index solve counting is available\n",
                bp->nsolves);
        profile_mode = "serial-nsolves";
    }

    if (index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        il_size(bp->fieldlist) == 1 &&
        bp->nsolves <= 1 &&
        (sp->maxquads != 0 || sp->maxmatches != 0)) {
        logverb("[index-shard] maxquads=%i maxmatches=%i use exact serial "
                "path until process-wide hypothesis limits are available\n",
                sp->maxquads,
                sp->maxmatches);
        profile_mode = "serial-hypothesis-limits";
    }

    if (index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        il_size(bp->fieldlist) == 1 &&
        bp->nsolves <= 1 &&
        sp->maxquads == 0 &&
        sp->maxmatches == 0 &&
        bp->rdls_tagalong_all) {
        logverb("[index-shard] automatic RDLS tag-along discovery uses exact "
                "serial path until its column list is worker-private\n");
        profile_mode = "serial-rdls-tagalong-all";
    }

    if (index_shard_pthread_enabled(bp) &&
        index_shard_pool_active(bp) &&
        verification_datalog) {
        logverb("[index-shard] process-global verification datalog uses "
                "the exact serial path\n");
        profile_mode = "serial-verification-datalog";
    }

    // Start solving...
    if (bp->indexes_inparallel) {

        // Add all the indexes...
        for (I=0; I<Nindexes; I++) {
            double phase_wall_start = monotonic_seconds();
            index_t* index = onefield_internal_get_index(bp, I);

            profile_acquire_seconds +=
                monotonic_seconds() - phase_wall_start;
            solver_add_index(sp, index);
        }

        // Record current CPU usage.
        bp->cpu_start = get_cpu_usage();
        // Record current wall-clock time.
        bp->time_start = monotonic_seconds();

        // Do it!
        {
            double phase_wall_start = monotonic_seconds();

            solve_fields(bp, NULL);
            profile_solver_seconds +=
                monotonic_seconds() - phase_wall_start;
        }

        profile_indexes_executed = Nindexes;

        if (bp->solver_failed || sp->profile.execution_failed) {
            bp->solver_failed = TRUE;
        }

        // Clean up the indices...
        for (I=0; I<Nindexes; I++) {
            double phase_wall_start = monotonic_seconds();
            index_t* index = solver_get_index(sp, I);

            /*
             * Release the exact handle acquired above.  Re-entering
             * onefield_internal_get_index() here would open a second copy when indexnames owns
             * the admission list, leaking the mappings actually used by the
             * grouped solver.
             */
            if (onefield_internal_done_with_index(bp, I, index)) {
                bp->solver_failed = TRUE;
            }
            profile_release_seconds +=
                monotonic_seconds() - phase_wall_start;
        }
        solver_clear_indexes(sp);

    } else {

        for (I=0; I<Nindexes; I++) {
            index_t* index;

            /*
             * The serial path creates a fresh solver_run() for every index.
             * Poll the process-wide limits here so short per-index runs cannot
             * indefinitely postpone the native deadline callback.
             */
            if (onefield_check_total_limits(bp)) {
                break;
            }
            if (bp->single_field_solved) {
                break;
            }
            if (bp->cancelled) {
                break;
            }

            // Load the index...
            {
                double phase_wall_start = monotonic_seconds();

                index = onefield_internal_get_index(bp, I);
                profile_acquire_seconds +=
                    monotonic_seconds() - phase_wall_start;
            }

            /*
             * Index loading can fault substantial mapped data. If the native
             * deadline expired during acquisition, release the index without
             * entering the scalar solver.
             */
            if (onefield_check_total_limits(bp)) {
                if (onefield_internal_done_with_index(bp, I, index)) {
                    bp->solver_failed = TRUE;
                }
                break;
            }

            solver_add_index(sp, index);
            logverb("Trying index %s...\n", index->indexname);

            // Record current CPU usage.
            bp->cpu_start = get_cpu_usage();
            // Record current wall-clock time.
            bp->time_start = monotonic_seconds();

            // Do it!
            {
                double phase_wall_start = monotonic_seconds();

                solve_fields(bp, NULL);
                profile_solver_seconds +=
                    monotonic_seconds() - phase_wall_start;
            }

            profile_indexes_executed++;

            if (bp->solver_failed || sp->profile.execution_failed) {
                bp->solver_failed = TRUE;
                logerr("Solver execution failed while trying index %s\n",
                       index->indexname ? index->indexname : "(null)");
            }

            /*
             * Persist any total-limit decision across the next solver reset.
             * This check runs once per completed index, not in CodeKD.
             */
            onefield_check_total_limits(bp);

            // Clean up this index...
            {
                double phase_wall_start = monotonic_seconds();

                if (onefield_internal_done_with_index(bp, I, index)) {
                    bp->solver_failed = TRUE;
                }
                profile_release_seconds +=
                    monotonic_seconds() - phase_wall_start;
            }
            solver_clear_indexes(sp);

            if (bp->solver_failed) {
                break;
            }

            if (bp->hit_total_timelimit || bp->hit_total_cpulimit) {
                break;
            }
        }
    }

 cleanup:
    // Clean up.
    if (!onefield_internal_field_cache_valid(bp) ||
        il_size(bp->fieldlist) != 1) {
        if (bp->xyls) {
            xylist_close(bp->xyls);
            bp->xyls = NULL;
        }
    }

    if (bp->solver_failed) {
        logerr("Suppressing solution output after solver execution failure\n");
    } else {
        double phase_wall_start = monotonic_seconds();

        if (write_solutions(bp)) {
            exit(-1);
        }

        if (publish_solved_fields(bp)) {
            bp->solver_failed = TRUE;
            logerr("Solution output completed, but solved-marker publication "
                   "failed\n");
        }

        profile_output_seconds =
            monotonic_seconds() - phase_wall_start;
    }

    logverb("[onefield-profile] mode=%s candidates=%zu serial_executed=%zu "
            "acquire=%.6f solver=%.6f release=%.6f output=%.6f "
            "total=%.6f failed=%i cancelled=%i\n",
            profile_mode,
            Nindexes,
            profile_indexes_executed,
            profile_acquire_seconds,
            profile_solver_seconds,
            profile_release_seconds,
            profile_output_seconds,
            monotonic_seconds() - profile_wall_start,
            bp->solver_failed ? 1 : 0,
            bp->cancelled ? 1 : 0);

    for (i=0; i<bl_size(bp->solutions); i++) {
        MatchObj* mo = bl_access(bp->solutions, i);
        verify_free_matchobj(mo);
        onefield_free_matchobj(mo);
    }
    bl_remove_all(bp->solutions);
}

void onefield_init(onefield_t* bp) {
    // Reset params.
    memset(bp, 0, sizeof(onefield_t));

    bp->fieldlist = il_new(256);
    bp->solutions = bl_new(16, sizeof(MatchObj));
    bp->solved_fields_pending = il_new(4);
    bp->indexnames = sl_new(16);
    bp->indexes = pl_new(16);
    bp->owned_indexes = pl_new(4);
    bp->verify_wcs_list = bl_new(1, sizeof(sip_t));
    bp->verify_wcsfiles = sl_new(1);
    bp->fieldid_key = strdup("FIELDID");
    onefield_set_xcol(bp, NULL);
    onefield_set_ycol(bp, NULL);
    bp->quad_size_fraction_lo = DEFAULT_QSF_LO;
    bp->quad_size_fraction_hi = DEFAULT_QSF_HI;
    bp->nsolves = 1;
    bp->index_shard_workers = 1;

    bp->xyls_tagalong_all = TRUE;
    // don't set sp-> here because solver_set_default_values()
    // will get called next and wipe it out...
}

int onefield_parameters_are_okay(onefield_t* bp, solver_t* sp) {
    if (sp->distractor_ratio == 0) {
        logerr("You must set a \"distractors\" proportion.\n");
        return 0;
    }
    if (!(sl_size(bp->indexnames) || pl_size(bp->indexes))) {
        logerr("You must specify one or more indexes.\n");
        return 0;
    }
    if (!bp->fieldfname) {
        logerr("You must specify a field filename (xylist).\n");
        return 0;
    }
    if (sp->codetol < 0.0) {
        logerr("You must specify codetol > 0\n");
        return 0;
    }
    if (sp->verify_pix <= 0.0) {
        logerr("You must specify a positive verify_pix.\n");
        return 0;
    }
    if ((sp->funits_lower != 0.0) && (sp->funits_upper != 0.0) &&
        (sp->funits_lower > sp->funits_upper)) {
        logerr("fieldunits_lower MUST be less than fieldunits_upper.\n");
        logerr("\n(in other words, the lower-bound of scale estimate must "
               "be less than the upper-bound!)\n\n");
        return 0;
    }
    return 1;
}

int onefield_is_run_obsolete(onefield_t* bp, solver_t* sp) {
  if (bp->single_field_solved)
    return 1;

  if (bp->cancelled)
    return 1;

  if (bp->hit_total_cpulimit || bp->hit_total_timelimit)
    return 1;
  // If we're just solving one field, check to see if it's already
  // solved before doing a bunch of work and spewing tons of output.
  if ((il_size(bp->fieldlist) == 1) && bp->solved_in) {
    if (is_field_solved(bp, il_get(bp->fieldlist, 0)))
      return 1;
    }
    // Early check to see if this job was cancelled.
    if (bp->cancelfname) {
        if (file_exists(bp->cancelfname)) {
            logerr("Run cancelled.\n");
            return 1;
        }
    }

    return 0;
}

static void load_and_parse_wcsfiles(onefield_t* bp) {
    int i;
    for (i = 0; i < sl_size(bp->verify_wcsfiles); i++) {
        sip_t wcs;
        char* fn = sl_get(bp->verify_wcsfiles, i);
        logmsg("Reading WCS header to verify from file %s\n", fn);
        memset(&wcs, 0, sizeof(sip_t));
        if (!sip_read_header_file(fn, &wcs)) {
            logerr("Failed to parse WCS header from file %s\n", fn);
            continue;
        }
        bl_append(bp->verify_wcs_list, &wcs);
    }
}

void onefield_log_run_parameters(onefield_t* bp) {
    solver_t* sp = &(bp->solver);
    int i, N;

    logverb("solver run parameters:\n");
    logverb("indexes:\n");
    N = onefield_internal_index_count(bp);
    for (i=0; i<N; i++)
        logverb("  %s\n", onefield_internal_get_index_name(bp, i));
    if (bp->fieldfname)
        logverb("fieldfname %s\n", bp->fieldfname);
    logverb("fields ");
    for (i = 0; i < il_size(bp->fieldlist); i++)
        logverb("%i ", il_get(bp->fieldlist, i));
    logverb("\n");
    for (i = 0; i < sl_size(bp->verify_wcsfiles); i++)
        logverb("verify %s\n", sl_get(bp->verify_wcsfiles, i));
    logverb("fieldid %i\n", bp->fieldid);
    if (bp->matchfname)
        logverb("matchfname %s\n", bp->matchfname);
    if (bp->solved_in)
        logverb("solved_in %s\n", bp->solved_in);
    if (bp->solved_out)
        logverb("solved_out %s\n", bp->solved_out);
    if (bp->cancelfname)
        logverb("cancel %s\n", bp->cancelfname);
    if (bp->wcs_template)
        logverb("wcs %s\n", bp->wcs_template);
    if (bp->fieldid_key)
        logverb("fieldid_key %s\n", bp->fieldid_key);
    if (bp->indexrdlsfname)
        logverb("indexrdlsfname %s\n", bp->indexrdlsfname);
    logverb("parity %i\n", sp->parity);
    logverb("codetol %g\n", sp->codetol);
    logverb("startdepth %i\n", sp->startobj);
    logverb("enddepth %i\n", sp->endobj);
    logverb("fieldunits_lower %g\n", sp->funits_lower);
    logverb("fieldunits_upper %g\n", sp->funits_upper);
    logverb("verify_pix %g\n", sp->verify_pix);
    if (bp->xcolname)
        logverb("xcolname %s\n", bp->xcolname);
    if (bp->ycolname)
        logverb("ycolname %s\n", bp->ycolname);
    logverb("maxquads %i\n", sp->maxquads);
    logverb("maxmatches %i\n", sp->maxmatches);
    logverb("p_workers %i\n", bp->index_shard_workers);
    logverb("cpulimit %f\n", bp->cpulimit);
    logverb("timelimit %g\n", bp->timelimit);
    logverb("total_timelimit %g\n", bp->total_timelimit);
    logverb("total_cpulimit %f\n", bp->total_cpulimit);
}

void onefield_cleanup(onefield_t* bp) {
    onefield_job_field_cache_end(bp);
    if (bp->xyls) {
        xylist_close(bp->xyls);
        bp->xyls = NULL;
    }
    onefield_clear_indexes(bp);
    il_free(bp->fieldlist);
    il_free(bp->solved_fields_pending);
    bl_free(bp->solutions);
    sl_free2(bp->indexnames);
    pl_free(bp->indexes);
    pl_free(bp->owned_indexes);
    sl_free2(bp->verify_wcsfiles);
    bl_free(bp->verify_wcs_list);
    sl_free2(bp->rdls_tagalong);

    free(bp->cancelfname);
    free(bp->fieldfname);
    free(bp->fieldid_key);
    free(bp->indexrdlsfname);
    free(bp->scamp_fname);
    free(bp->corr_fname);
    free(bp->matchfname);
    free(bp->solved_in);
    free(bp->solved_out);
    free(bp->wcs_template);
    free(bp->xcolname);
    free(bp->ycolname);
    free(bp->sort_rdls);
}

static int sort_rdls(MatchObj* mymo, onefield_t* bp) {
    const solver_t* sp = &(bp->solver);
    anbool asc = TRUE;
    char* colname = bp->sort_rdls;
    double* sortdata;
    fitstable_t* tagalong;
    int* perm;
    int i;
    logverb("Sorting RDLS by column \"%s\"\n", bp->sort_rdls);
    if (colname[0] == '-') {
        colname++;
        asc = FALSE;
    }
    tagalong = startree_get_tagalong(sp->index->starkd);
    if (!tagalong) {
        ERROR("Failed to find tag-along table in index");
        return -1;
    }
    sortdata = fitstable_read_column_inds(tagalong, colname, fitscolumn_double_type(),
                                          mymo->refstarid, mymo->nindex);
    if (!sortdata) {
        ERROR("Failed to read data for column \"%s\" in index", colname);
        return -1;
    }
    perm = permutation_init(NULL, mymo->nindex);
    permuted_sort(sortdata, sizeof(double), asc ? compare_doubles_asc : compare_doubles_desc,
                  perm, mymo->nindex);
    free(sortdata);

    if (mymo->refxyz)
        permutation_apply(perm, mymo->nindex, mymo->refxyz, mymo->refxyz, 3*sizeof(double));
    // probably not set yet, but what the heck...
    if (mymo->refradec)
        permutation_apply(perm, mymo->nindex, mymo->refradec,  mymo->refradec, 2*sizeof(double));
    if (mymo->refxy)
        permutation_apply(perm, mymo->nindex, mymo->refxy,     mymo->refxy,    2*sizeof(double));
    if (mymo->refstarid)
        permutation_apply(perm, mymo->nindex, mymo->refstarid, mymo->refstarid,  sizeof(int));
    if (mymo->theta)
        for (i=0; i<mymo->nfield; i++) {
            if (mymo->theta[i] < 0)
                continue;
            mymo->theta[i] = perm[mymo->theta[i]];
        }
    free(perm);
    return 0;
}

anbool onefield_internal_record_match_callback(MatchObj* mo, void* userdata) {
    onefield_t* bp = userdata;
    // shard poll
    index_shard_poll_from_callback(bp);
    solver_t* sp = &(bp->solver);
    MatchObj* mymo;
    int ind;

    check_time_limits(bp);

    // Copy "mo" to "mymo".
    ind = bl_insert_sorted(bp->solutions, mo, onefield_internal_compare_matchobjs);
    mymo = bl_access(bp->solutions, ind);

    // steal these arrays from "mo" (prevent them from being free()'d
    // by the caller)
    mo->theta = NULL;
    mo->matchodds = NULL;
    mo->refxyz = NULL;
    mo->refxy = NULL;
    mo->refstarid = NULL;
    mo->testperm = NULL;

    // We have no guarantee that the index will still be open when it
    // comes time to write our output files, so we've got to grab everything
    // we need now while it's at hand.

    if (bp->indexrdlsfname || bp->scamp_fname || bp->corr_fname) {
        int i;

        // This must happen first, because it reorders the "ref" arrays,
        // and we want that to be done before more data are integrated.
        if (bp->sort_rdls) {
            if (sort_rdls(mymo, bp)) {
                ERROR("Failed to sort RDLS file by column \"%s\"", bp->sort_rdls);
            }
        }

        logdebug("Converting %i reference stars from xyz to radec\n", mymo->nindex);
        mymo->refradec = malloc(mymo->nindex * 2 * sizeof(double));
        for (i=0; i<mymo->nindex; i++) {
            xyzarr2radecdegarr(mymo->refxyz+i*3, mymo->refradec+i*2);
            logdebug("  %i: radec %.2f,%.2f\n", i, mymo->refradec[i*2], mymo->refradec[i*2+1]);
        }

        mymo->fieldxy = malloc(mymo->nfield * 2 * sizeof(double));
        // whew!
        memcpy(mymo->fieldxy, bp->solver.vf->xy, mymo->nfield * 2 * sizeof(double));

        // Tweak was here...

        // FIXME -- add MAG, MAGERR, and positional errors for SCAMP catalog.

        if (bp->rdls_tagalong || bp->rdls_tagalong_all)
            grab_tagalong_data(sp->index->starkd, mymo, bp, mymo->refstarid, mymo->nindex);

        // FIXME -- we don't support specifying individual fields (yet)
        assert(bp->xyls_tagalong_all);
        assert(!bp->xyls_tagalong);
        if (bp->xyls_tagalong_all)
            grab_field_tagalong_data(mymo, bp->xyls, mymo->nfield);
    }

    if (mymo->logodds < bp->logratio_tosolve)
        return FALSE;

    // this match is considered a solution.

    bp->nsolves_sofar++;
    if (bp->nsolves_sofar < bp->nsolves) {
        logmsg("Found a quad that solves the image; that makes %i of %i required.\n",
               bp->nsolves_sofar, bp->nsolves);
    } else {
        if (!index_shard_worker_context_active()) {
            if (bp->solver.index) {
                char* base = basename_safe(bp->solver.index->indexname);
                logmsg("Field %i: solved with index %s.\n", mymo->fieldnum, base);
                free(base);
            } else {
                logmsg("Field %i: solved with index %i", mymo->fieldnum, mymo->indexid);
                if (mymo->healpix >= 0) {
                    logmsg(", healpix %i\n", mymo->healpix);
                } else {
                    logmsg("\n");
                }
            }
        }
        return TRUE;
    }
    return FALSE;
}

time_t onefield_internal_timer_callback(void* user_data) {
    onefield_t* bp = user_data;
    // shard poll
    index_shard_poll_from_callback(bp);
    check_time_limits(bp);

    // check if the field has already been solved...
    if (is_field_solved(bp, bp->fieldnum))
        return 0;
    if (bp->cancelfname && file_exists(bp->cancelfname)) {
        bp->cancelled = TRUE;
        logmsg("File \"%s\" exists: cancelling.\n", bp->cancelfname);
        return 0;
    }
    return 1; // wait 1 second... FIXME config?
}

static void add_onefield_params(onefield_t* bp, qfits_header* hdr) {
    solver_t* sp = &(bp->solver);
    int i;
    int Nindexes;
    fits_add_long_comment(hdr, "-- onefield solver parameters: --");
    if (sp->index) {
        fits_add_long_comment(hdr, "Index name: %s", sp->index->indexname?sp->index->indexname:"(null)");
        fits_add_long_comment(hdr, "Index id: %i", sp->index->indexid);
        fits_add_long_comment(hdr, "Index healpix: %i", sp->index->healpix);
        fits_add_long_comment(hdr, "Index healpix nside: %i", sp->index->hpnside);
        fits_add_long_comment(hdr, "Index scale lower: %g arcsec", sp->index->index_scale_lower);
        fits_add_long_comment(hdr, "Index scale upper: %g arcsec", sp->index->index_scale_upper);
        fits_add_long_comment(hdr, "Index jitter: %g", sp->index->index_jitter);
        fits_add_long_comment(hdr, "Circle: %s", sp->index->circle ? "yes" : "no");
        fits_add_long_comment(hdr, "Cxdx margin: %g", sp->cxdx_margin);
    }
    Nindexes = onefield_internal_index_count(bp);
    for (i = 0; i < Nindexes; i++)
        fits_add_long_comment(hdr, "Index(%i): %s", i, onefield_internal_get_index_name(bp, i)?onefield_internal_get_index_name(bp, i):"(null)");

    fits_add_long_comment(hdr, "Field name: %s", bp->fieldfname?bp->fieldfname:"(null)");
    fits_add_long_comment(hdr, "Field scale lower: %g arcsec/pixel", sp->funits_lower);
    fits_add_long_comment(hdr, "Field scale upper: %g arcsec/pixel", sp->funits_upper);
    fits_add_long_comment(hdr, "X col name: %s", bp->xcolname?bp->xcolname:"(null)");
    fits_add_long_comment(hdr, "Y col name: %s", bp->ycolname?bp->ycolname:"(null)");
    fits_add_long_comment(hdr, "Start obj: %i", sp->startobj);
    fits_add_long_comment(hdr, "End obj: %i", sp->endobj);

    // 'Solved_in' is often a NULL pointer.
    // If %s is a NULL pointer, vasprintf() causes a segmentation fault (due to
    // strlen()) on Solaris -> added treatment of this case for portability.
    // GNU/Linux implementation of vasprintf() catches NULL pointer and prints
    // "(null)" in header. Seems to be an issue on Solaris only.
    fits_add_long_comment(hdr, "Solved_in: %s", bp->solved_in?bp->solved_in:"(null)");
    fits_add_long_comment(hdr, "Solved_out: %s", bp->solved_out?bp->solved_out:"(null)");

    fits_add_long_comment(hdr, "Parity: %i", sp->parity);
    fits_add_long_comment(hdr, "Codetol: %g", sp->codetol);
    fits_add_long_comment(hdr, "Verify pixels: %g pix", sp->verify_pix);

    fits_add_long_comment(hdr, "Maxquads: %i", sp->maxquads);
    fits_add_long_comment(hdr, "Maxmatches: %i", sp->maxmatches);
    fits_add_long_comment(hdr, "Cpu limit: %f s", bp->cpulimit);
    fits_add_long_comment(hdr, "Time limit: %g s", bp->timelimit);
    fits_add_long_comment(hdr, "Total time limit: %g s", bp->total_timelimit);
    fits_add_long_comment(hdr, "Total CPU limit: %f s", bp->total_cpulimit);

    fits_add_long_comment(hdr, "Tweak: %s", (sp->do_tweak ? "yes" : "no"));
    if (sp->do_tweak) {
        fits_add_long_comment(hdr, "Tweak AB order: %i", sp->tweak_aborder);
        fits_add_long_comment(hdr, "Tweak ABP order: %i", sp->tweak_abporder);
    }

    fits_add_long_comment(hdr, "--");
}

void onefield_internal_remove_invalid_fields(il* fieldlist, int maxfield) {
    int i;
    for (i=0; i<il_size(fieldlist); i++) {
        int fieldnum = il_get(fieldlist, i);
        if (fieldnum >= 1 && fieldnum <= maxfield)
            continue;
        if (fieldnum > maxfield) {
            logerr("Field %i does not exist (max=%i).\n", fieldnum, maxfield);
        }
        if (fieldnum < 1) {
            logerr("Field %i is invalid (must be >= 1).\n", fieldnum);
        }
        il_remove(fieldlist, i);
        i--;
    }
}

static void solve_fields(onefield_t* bp, sip_t* verify_wcs) {
    solver_t* sp = &(bp->solver);
    solver_profile_t profile_total;
    double last_utime, last_stime;
    double utime, stime;
    struct timeval wtime, last_wtime;
    int fi;

    memset(&profile_total, 0, sizeof(profile_total));
    memset(&sp->profile, 0, sizeof(sp->profile));

    get_resource_stats(&last_utime, &last_stime, NULL);
    gettimeofday(&last_wtime, NULL);

    for (fi = 0; fi < il_size(bp->fieldlist); fi++) {
        int fieldnum;
        MatchObj template ;
        qfits_header* fieldhdr = NULL;
        double field_wall_start = monotonic_seconds();
        double field_read_seconds = 0.0;
        double preprocess_seconds = 0.0;

        fieldnum = il_get(bp->fieldlist, fi);

        if (onefield_check_total_limits(bp)) {
            break;
        }

        memset(&template, 0, sizeof(MatchObj));
        template.fieldnum = fieldnum;
        template.fieldfile = bp->fieldid;

        // Has the field already been solved?
        if (is_field_solved(bp, fieldnum)) {
            goto cleanup;
        }

        if (onefield_internal_prepare_field_view(
                bp,
                fieldnum,
                &field_read_seconds,
                &preprocess_seconds)) {
            bp->solver_failed = TRUE;
            profile_total.execution_failed = TRUE;
            goto cleanup;
        }
        if (onefield_check_total_limits(bp)) {
            goto cleanup;
        }

        /*
         * Read metadata from the same currently prepared XYLS epoch. Cache
         * invalidation may close and reopen the handle, so header access must
         * follow field preparation rather than precede it.
         */
        fieldhdr = xylist_get_header(bp->xyls);
        if (fieldhdr) {
            char* idstr =
                fits_get_dupstring(fieldhdr, bp->fieldid_key);
            if (idstr) {
                strncpy(
                    template.fieldname,
                    idstr,
                    sizeof(template.fieldname) - 1);
            }
            free(idstr);
        }

        sp->mo_template = &template;
        sp->record_match_callback = onefield_internal_record_match_callback;
        sp->timer_callback = onefield_internal_timer_callback;
        sp->userdata = bp;

        bp->fieldnum = fieldnum;
        bp->nsolves_sofar = 0;

        if (verify_wcs) {
            //MatchObj mo;
            logmsg("Verifying WCS of field %i.\n", fieldnum);
            solver_verify_sip_wcs(sp, verify_wcs); //, &mo);
            logmsg(" --> log-odds %g\n", sp->best_logodds);

        } else {
            logverb("Solving field %i.\n", fieldnum);
            sp->distance_from_quad_bonus = TRUE;
            solver_log_params(sp);

            // The real thing
            if (solver_run(sp)) {
                bp->solver_failed = TRUE;
            }

            solver_profile_accumulate(&profile_total, &sp->profile);

            logverb("[onefield-field-profile] field=%i read=%.6f "
                    "preprocess=%.6f solver_run=%.6f total=%.6f "
                    "failed=%i\n",
                    fieldnum,
                    field_read_seconds,
                    preprocess_seconds,
                    sp->profile.solver_run_wall_seconds,
                    monotonic_seconds() - field_wall_start,
                    sp->profile.execution_failed ? 1 : 0);

            if (bp->solver_failed || sp->profile.execution_failed) {
                bp->solver_failed = TRUE;
                logerr("Solver execution failed for field %i\n",
                       fieldnum);
                goto cleanup;
            }

            logverb("Field %i: tried %i quads, matched %i codes.\n",
                    fieldnum, sp->numtries, sp->nummatches);

            if (sp->maxquads && sp->numtries >= sp->maxquads)
                logmsg("  exceeded the number of quads to try: %i >= %i.\n",
                       sp->numtries, sp->maxquads);
            if (sp->maxmatches && sp->nummatches >= sp->maxmatches)
                logmsg("  exceeded the number of quads to match: %i >= %i.\n",
                       sp->nummatches, sp->maxmatches);
            if (bp->cancelled)
                logmsg("  cancelled at user request.\n");
        }


        if (sp->best_match_solves) {
            onefield_internal_solved_field(bp, fieldnum);
        } else if (!verify_wcs) {
            // Field unsolved.
            if (bp->solver.index && bp->solver.index->indexname) {
                char* copy;
                char* base;
                copy = strdup_safe(bp->solver.index->indexname);
                base = basename(copy);
                if (bp->solver.endobj) {
                    logerr("Field %i did not solve (index %s, field objects %i-%i).\n",
                           fieldnum, base, bp->solver.startobj+1, bp->solver.endobj);
                } else {
                    logerr("Field %i did not solve (index %s).\n", fieldnum, base);
                }
                free(copy);
            } else {
                logerr("Field %i did not solve.\n", fieldnum);
            }
            if (sp->have_best_match) {
                logverb("Best match encountered: ");
                matchobj_print(&(sp->best_match), log_get_level());
            } else {
                logverb("Best odds encountered: %g\n", exp(sp->best_logodds));
            }
        }

        get_resource_stats(&utime, &stime, NULL);
        gettimeofday(&wtime, NULL);
        logverb("Spent %g s user, %g s system, %g s total, %g s wall time.\n",
                (utime - last_utime), (stime - last_stime),
                (stime - last_stime + utime - last_utime),
                millis_between(&last_wtime, &wtime) * 0.001);

        last_utime = utime;
        last_stime = stime;
        last_wtime = wtime;

    cleanup:
        sp->mo_template = NULL;
        sp->record_match_callback = NULL;
        sp->timer_callback = NULL;
        sp->userdata = NULL;
        if (onefield_internal_field_cache_has_field(bp, fieldnum) &&
            il_size(bp->fieldlist) == 1 &&
            sp->fieldxy_orig && sp->fieldxy && sp->vf) {
            onefield_internal_reset_field_pass_state(bp);
        } else {
            solver_cleanup_field(sp);
        }

        if (bp->solver_failed) {
            break;
        }
    }

    sp->profile = profile_total;
}

static anbool is_field_solved(onefield_t* bp, int fieldnum) {
    anbool solved = FALSE;

    if (bp->solved_fields_pending &&
        il_sorted_contains(bp->solved_fields_pending, fieldnum)) {
        logverb("Field %i has already been solved in this run.\n", fieldnum);
        return TRUE;
    }

    if (bp->solved_in) {
        solved = solvedfile_get(bp->solved_in, fieldnum);
        logverb("Checking %s file %i to see if the field is solved: %s.\n",
                bp->solved_in, fieldnum, (solved ? "yes" : "no"));
    }
    if (solved) {
        // file exists; field has already been solved.
        logmsg("Field %i: solvedfile %s: field has been solved.\n", fieldnum, bp->solved_in);
        return TRUE;
    }
    return FALSE;
}

void onefield_internal_solved_field(onefield_t* bp, int fieldnum) {
    if (bp->solved_fields_pending) {
        il_insert_unique_ascending(bp->solved_fields_pending, fieldnum);
    }

    // If we're just solving a single field, and we solved it...
    if (il_size(bp->fieldlist) == 1)
        bp->single_field_solved = TRUE;
}

/*
 * Publish external solved markers only after the complete solve pass has
 * quiesced and final solution output succeeded. First-valid selection stops
 * workers before the reducer commits in-memory solved state, and neither event
 * may leave a marker that makes a failed output retry skip the input.
 */
static int publish_solved_fields(onefield_t* bp) {
    int i;
    int failed = FALSE;

    if (!bp || !bp->solved_out || !bp->solved_fields_pending) {
        return 0;
    }

    for (i = 0; i < il_size(bp->solved_fields_pending); i++) {
        int fieldnum = il_get(bp->solved_fields_pending, i);

        logmsg("Field %i solved: writing to file %s to indicate this.\n",
               fieldnum,
               bp->solved_out);

        if (solvedfile_set(bp->solved_out, fieldnum)) {
            logerr("Failed to write solvedfile %s.\n", bp->solved_out);
            failed = TRUE;
        }
    }

    return failed ? -1 : 0;
}

void onefield_matchobj_deep_copy(const MatchObj* mo, MatchObj* dest) {
    if (!mo || !dest)
        return;
    if (mo->sip) {
        dest->sip = sip_create();
        memcpy(dest->sip, mo->sip, sizeof(sip_t));
    }
    if (mo->refradec) {
        dest->refradec = malloc(mo->nindex * 2 * sizeof(double));
        memcpy(dest->refradec, mo->refradec, mo->nindex * 2 * sizeof(double));
    }
    if (mo->fieldxy) {
        dest->fieldxy = malloc(mo->nfield * 2 * sizeof(double));
        memcpy(dest->fieldxy, mo->fieldxy, mo->nfield * 2 * sizeof(double));
    }
    if (mo->tagalong) {
        int i;
        dest->tagalong = bl_new(16, sizeof(tagalong_t));
        for (i=0; i<bl_size(mo->tagalong); i++) {
            tagalong_t* tag = bl_access(mo->tagalong, i);
            tagalong_t tagcopy;
            memcpy(&tagcopy, tag, sizeof(tagalong_t));
            tagcopy.name = strdup_safe(tag->name);
            tagcopy.units = strdup_safe(tag->units);
            if (tag->data) {
                tagcopy.data = malloc((size_t)tag->Ndata * (size_t)tag->itemsize);
                memcpy(tagcopy.data, tag->data, (size_t)tag->Ndata * (size_t)tag->itemsize);
            }
            bl_append(dest->tagalong, &tagcopy);
        }
    }
    // NOT SUPPORTED (yet)
    assert(!mo->field_tagalong);
}

// Free the things I added to the mo.
void onefield_free_matchobj(MatchObj* mo) {
    if (!mo) return;
    if (mo->sip) {
        sip_free(mo->sip);
        mo->sip = NULL;
    }
    free(mo->refradec);
    free(mo->fieldxy);
    free(mo->theta);
    free(mo->matchodds);
    free(mo->refxyz);
    free(mo->refxy);
    free(mo->refstarid);
    free(mo->testperm);
    mo->refradec = NULL;
    mo->fieldxy = NULL;
    mo->theta = NULL;
    mo->matchodds = NULL;
    mo->refxyz = NULL;
    mo->refxy = NULL;
    mo->refstarid = NULL;
    mo->testperm = NULL;

    if (mo->tagalong) {
        int i;
        for (i=0; i<bl_size(mo->tagalong); i++) {
            tagalong_t* tag = bl_access(mo->tagalong, i);
            free(tag->name);
            free(tag->units);
            free(tag->data);
        }
        bl_free(mo->tagalong);
        mo->tagalong = NULL;
    }
    if (mo->field_tagalong) {
        int i;
        for (i=0; i<bl_size(mo->field_tagalong); i++) {
            tagalong_t* tag = bl_access(mo->field_tagalong, i);
            free(tag->name);
            free(tag->units);
            free(tag->data);
        }
        bl_free(mo->field_tagalong);
        mo->field_tagalong = NULL;
    }
}

static void remove_duplicate_solutions(onefield_t* bp) {
    int i, j;
    // The solutions can fall out of order because tweak2() updates their logodds.
    bl_sort(bp->solutions, onefield_internal_compare_matchobjs);

    for (i=0; i<bl_size(bp->solutions); i++) {
        MatchObj* mo = bl_access(bp->solutions, i);
        j = i+1;
        while (j < bl_size(bp->solutions)) {
            MatchObj* mo2 = bl_access(bp->solutions, j);
            if (mo->fieldfile != mo2->fieldfile)
                break;
            if (mo->fieldnum != mo2->fieldnum)
                break;
            assert(mo2->logodds <= mo->logodds);
            onefield_free_matchobj(mo2);
            verify_free_matchobj(mo2);
            bl_remove_index(bp->solutions, j);
        }
    }
}

static int write_match_file(onefield_t* bp) {
    int i;
    bp->mf = matchfile_open_for_writing(bp->matchfname);
    if (!bp->mf) {
        logerr("Failed to open file %s to write match file.\n", bp->matchfname);
        return -1;
    }
    BOILERPLATE_ADD_FITS_HEADERS(bp->mf->header);
    qfits_header_add(bp->mf->header, "DATE", qfits_get_datetime_iso8601(), "Date this file was created.", NULL);
    add_onefield_params(bp, bp->mf->header);
    if (matchfile_write_headers(bp->mf)) {
        logerr("Failed to write matchfile header.\n");
        return -1;
    }
    for (i=0; i<bl_size(bp->solutions); i++) {
        MatchObj* mo = bl_access(bp->solutions, i);
        if (matchfile_write_match(bp->mf, mo)) {
            logerr("Field %i: error writing a match.\n", mo->fieldnum);
            return -1;
        }
    }
    if (matchfile_fix_headers(bp->mf) ||
        matchfile_close(bp->mf)) {
        logerr("Error closing matchfile.\n");
        return -1;
    }
    bp->mf = NULL;
    return 0;
}

static int write_rdls_file(onefield_t* bp) {
    int i;
    qfits_header* h;
    bp->indexrdls = rdlist_open_for_writing(bp->indexrdlsfname);
    if (!bp->indexrdls) {
        logerr("Failed to open index RDLS file %s for writing.\n",
               bp->indexrdlsfname);
        return -1;
    }
    h = rdlist_get_primary_header(bp->indexrdls);

    BOILERPLATE_ADD_FITS_HEADERS(h);
    fits_add_long_history(h, "This \"indexrdls\" file contains the RA/DEC of index objects that were found inside a solved field.");
    qfits_header_add(h, "DATE", qfits_get_datetime_iso8601(), "Date this file was created.", NULL);
    add_onefield_params(bp, h);
    if (rdlist_write_primary_header(bp->indexrdls)) {
        logerr("Failed to write index RDLS header.\n");
        return -1;
    }

    for (i=0; i<bl_size(bp->solutions); i++) {
        MatchObj* mo = bl_access(bp->solutions, i);
        rd_t rd;
        if (strlen(mo->fieldname)) {
            qfits_header* hdr = rdlist_get_header(bp->indexrdls);
            qfits_header_add(hdr, "FIELDID", mo->fieldname, "Name of this field", NULL);
        }
        if (mo->tagalong) {
            int j;
            for (j=0; j<bl_size(mo->tagalong); j++) {
                tagalong_t* tag = bl_access(mo->tagalong, j);
                tag->colnum = rdlist_add_tagalong_column(bp->indexrdls, tag->type, tag->arraysize,
                                                         tag->type, tag->name, tag->units);
            }
        }
        if (rdlist_write_header(bp->indexrdls)) {
            logerr("Failed to write index RDLS field header.\n");
            return -1;
        }
        assert(mo->refradec);

        rd_from_array(&rd, mo->refradec, mo->nindex);
        if (rdlist_write_field(bp->indexrdls, &rd)) {
            logerr("Failed to write index RDLS entry.\n");
            return -1;
        }
        rd_free_data(&rd);

        if (mo->tagalong) {
            int j;
            for (j=0; j<bl_size(mo->tagalong); j++) {
                tagalong_t* tag = bl_access(mo->tagalong, j);
                if (rdlist_write_tagalong_column(bp->indexrdls, tag->colnum,
                                                 0, mo->nindex, tag->data, tag->itemsize)) {
                    ERROR("Failed to write tag-along data column %s", tag->name);
                    return -1;
                }
            }
        }

        if (rdlist_fix_header(bp->indexrdls)) {
            logerr("Failed to fix index RDLS field header.\n");
            return -1;
        }
        rdlist_next_field(bp->indexrdls);
    }

    if (rdlist_fix_primary_header(bp->indexrdls) ||
        rdlist_close(bp->indexrdls)) {
        logerr("Failed to close index RDLS file.\n");
        return -1;
    }
    bp->indexrdls = NULL;
    return 0;
}

static int write_wcs_file(onefield_t* bp) {
    int i;
    for (i=0; i<bl_size(bp->solutions); i++) {
        char wcs_fn[1024];
        FILE* fout;
        qfits_header* hdr;
        char* tm;

        MatchObj* mo = bl_access(bp->solutions, i);
        snprintf(wcs_fn, sizeof(wcs_fn), bp->wcs_template, mo->fieldnum);
        fout = fopen(wcs_fn, "wb");
        if (!fout) {
            logerr("Failed to open WCS output file %s: %s\n", wcs_fn, strerror(errno));
            return -1;
        }
        assert(mo->wcs_valid);

        if (mo->sip)
            hdr = sip_create_header(mo->sip);
        else
            hdr = tan_create_header(&(mo->wcstan));

        BOILERPLATE_ADD_FITS_HEADERS(hdr);
        qfits_header_add(hdr, "HISTORY", "This is a WCS header was created by Astrometry.net.", NULL, NULL);
        tm = qfits_get_datetime_iso8601();
        qfits_header_add(hdr, "DATE", tm, "Date this file was created.", NULL);
        add_onefield_params(bp, hdr);
        fits_add_long_comment(hdr, "-- properties of the matching quad: --");
        fits_add_long_comment(hdr, "index id: %i", mo->indexid);
        fits_add_long_comment(hdr, "index healpix: %i", mo->healpix);
        fits_add_long_comment(hdr, "index hpnside: %i", mo->hpnside);
        fits_add_long_comment(hdr, "log odds: %g", mo->logodds);
        fits_add_long_comment(hdr, "odds: %g", exp(mo->logodds));
        fits_add_long_comment(hdr, "quadno: %i", mo->quadno);
        fits_add_long_comment(hdr, "stars: %i,%i,%i,%i", mo->star[0], mo->star[1], mo->star[2], mo->star[3]);
        fits_add_long_comment(hdr, "field: %i,%i,%i,%i", mo->field[0], mo->field[1], mo->field[2], mo->field[3]);
        fits_add_long_comment(hdr, "code error: %g", sqrt(mo->code_err));
        fits_add_long_comment(hdr, "nmatch: %i", mo->nmatch);
        fits_add_long_comment(hdr, "nconflict: %i", mo->nconflict);
        fits_add_long_comment(hdr, "nfield: %i", mo->nfield);
        fits_add_long_comment(hdr, "nindex: %i", mo->nindex);
        fits_add_long_comment(hdr, "scale: %g arcsec/pix", mo->scale);
        fits_add_long_comment(hdr, "parity: %i", (int)mo->parity);
        fits_add_long_comment(hdr, "quads tried: %i", mo->quads_tried);
        fits_add_long_comment(hdr, "quads matched: %i", mo->quads_matched);
        fits_add_long_comment(hdr, "quads verified: %i", mo->nverified);
        fits_add_long_comment(hdr, "objs tried: %i", mo->objs_tried);
        fits_add_long_comment(hdr, "cpu time: %g", mo->timeused);
        fits_add_long_comment(hdr, "--");

        if (strlen(mo->fieldname))
            qfits_header_add(hdr, bp->fieldid_key, mo->fieldname, "Field name (copied from input field)", NULL);

        if (qfits_header_dump(hdr, fout)) {
            logerr("Failed to write FITS WCS header.\n");
            return -1;
        }
        fits_pad_file(fout);
        qfits_header_destroy(hdr);
        fclose(fout);
    }
    return 0;
}

static int write_scamp_file(onefield_t* bp) {
    int i;
    scamp_cat_t* scamp;
    qfits_header* hdr = NULL;
    MatchObj* mo;
    tan_t fakewcs;

    // HACK -- just hdr = NULL?
    hdr = qfits_header_default();
    fits_header_add_int(hdr, "BITPIX", 0, NULL);
    fits_header_add_int(hdr, "NAXIS", 2, NULL);
    fits_header_add_int(hdr, "NAXIS1", 0, NULL);
    fits_header_add_int(hdr, "NAXIS2", 0, NULL);
    qfits_header_add(hdr, "EXTEND", "T", "", NULL);
    memset(&fakewcs, 0, sizeof(tan_t));
    tan_add_to_header(hdr, &fakewcs);

    scamp = scamp_catalog_open_for_writing(bp->scamp_fname, TRUE);
    if (!scamp) {
        logerr("Failed to open SCAMP reference catalog for writing.\n");
        return -1;
    }
    if (scamp_catalog_write_field_header(scamp, hdr)) {
        logerr("Failed to write SCAMP headers.\n");
        return -1;
    }
    mo = bl_access(bp->solutions, 0);
    for (i=0; i<mo->nindex; i++) {
        scamp_ref_t ref;
        ref.ra  = mo->refradec[2*i + 0];
        ref.dec = mo->refradec[2*i + 1];
        ref.err_a = ref.err_b = arcsec2deg(mo->index_jitter);
        // HACK
        ref.mag = 10.0;
        ref.err_mag = 0.1;

        if (scamp_catalog_write_reference(scamp, &ref)) {
            logerr("Failed to write SCAMP object.\n");
            return -1;
        }
    }
    if (scamp_catalog_close(scamp)) {
        logerr("Failed to close SCAMP reference catalog.\n");
        return -1;
    }
    return 0;
}

static int write_corr_file(onefield_t* bp) {
    int i;
    fitstable_t* tab;
    tab = fitstable_open_for_writing(bp->corr_fname);
    if (!tab) {
        ERROR("Failed to open correspondences file \"%s\" for writing", bp->corr_fname);
        return -1;
    }
    // FIXME -- add header boilerplate.

    if (fitstable_write_primary_header(tab)) {
        ERROR("Failed to write primary header for corr file \"%s\"", bp->corr_fname);
        return -1;
    }

    for (i=0; i<bl_size(bp->solutions); i++) {
        MatchObj* mo;
        sip_t thesip;
        sip_t* wcs;
        int j;
        tfits_type dubl = fitscolumn_double_type();
        tfits_type itype = fitscolumn_int_type();

        mo = bl_access(bp->solutions, i);

        if (mo->sip)
            wcs = mo->sip;
        else {
            sip_wrap_tan(&mo->wcstan, &thesip);
            wcs = &thesip;
        }

        fitstable_add_write_column(tab, dubl, "field_x",   "pixels");
        fitstable_add_write_column(tab, dubl, "field_y",   "pixels");
        fitstable_add_write_column(tab, dubl, "field_ra",  "degrees");
        fitstable_add_write_column(tab, dubl, "field_dec", "degrees");
        fitstable_add_write_column(tab, dubl, "index_x",   "pixels");
        fitstable_add_write_column(tab, dubl, "index_y",   "pixels");
        fitstable_add_write_column(tab, dubl, "index_ra",  "degrees");
        fitstable_add_write_column(tab, dubl, "index_dec", "degrees");
        fitstable_add_write_column(tab, itype, "index_id", "none");
        fitstable_add_write_column(tab, itype, "field_id", "none");
        fitstable_add_write_column(tab, dubl, "match_weight", "none");

        if (mo->tagalong) {
            for (j=0; j<bl_size(mo->tagalong); j++) {
                tagalong_t* tag = bl_access(mo->tagalong, j);
                fitstable_add_write_column_struct(tab, tag->type, tag->arraysize, 0, tag->type, tag->name, tag->units);
                tag->colnum = fitstable_ncols(tab)-1;
            }
        }

        // FIXME -- check for duplicate column names
        if (mo->field_tagalong) {
            int j;
            for (j=0; j<bl_size(mo->field_tagalong); j++) {
                tagalong_t* tag = bl_access(mo->field_tagalong, j);
                fitstable_add_write_column_struct(tab, tag->type, tag->arraysize, 0, tag->type, tag->name, tag->units);
                tag->colnum = fitstable_ncols(tab)-1;
            }
        }

        if (fitstable_write_header(tab)) {
            ERROR("Failed to write correspondence file header.");
            return -1;
        }

        {
            int rows = 0;
            for (j=0; j<mo->nfield; j++) {
                if (mo->theta[j] < 0)
                    continue;
                rows++;
            }
            logverb("Writing %i rows (of %i field and %i index objects) to correspondence file.\n", rows, mo->nfield, mo->nindex);
        }
        for (j=0; j<mo->nfield; j++) {
            double fx,fy,fra,fdec;
            double rx,ry,rra,rdec;
            double weight;
            int ti, ri;
            ri = mo->theta[j];
            if (ri < 0)
                continue;
            ti = j;
            rra  = mo->refradec[2*ri+0];
            rdec = mo->refradec[2*ri+1];
            if (!sip_radec2pixelxy(wcs, rra, rdec, &rx, &ry))
                continue;
            fx = mo->fieldxy[2*ti+0];
            fy = mo->fieldxy[2*ti+1];
            sip_pixelxy2radec(wcs, fx, fy, &fra, &fdec);
            logdebug("Writing field xy %.1f,%.1f, radec %.2f,%.2f; index xy %.1f,%.1f, radec %.2f,%.2f\n", fx, fy, fra, fdec, rx, ry, rra, rdec);
            weight = verify_logodds_to_weight(mo->matchodds[j]);
            if (fitstable_write_row(tab, &fx, &fy, &fra, &fdec, &rx, &ry, &rra, &rdec, &ri, &ti, &weight)) {
                ERROR("Failed to write coordinates to correspondences file \"%s\"", bp->corr_fname);
                return -1;
            }
        }

        if (mo->tagalong) {
            for (j=0; j<bl_size(mo->tagalong); j++) {
                tagalong_t* tag = bl_access(mo->tagalong, j);
                int row = 0;
                int k;
                // Ugh, we write each datum individually...
                for (k=0; k<mo->nfield; k++) {
                    int ri = mo->theta[k];
                    if (ri < 0)
                        continue;
                    fitstable_write_one_column(tab, tag->colnum, row, 1,
                                               (char*)tag->data + ri*tag->itemsize, 0);
                    row++;
                }
            }
        }
        if (mo->field_tagalong) {
            for (j=0; j<bl_size(mo->field_tagalong); j++) {
                tagalong_t* tag = bl_access(mo->field_tagalong, j);
                int row = 0;
                int k;
                // Ugh, we write each datum individually...
                for (k=0; k<mo->nfield; k++) {
                    if (mo->theta[k] < 0)
                        continue;
                    fitstable_write_one_column(tab, tag->colnum, row, 1,
                                               (char*)tag->data + k*tag->itemsize, 0);
                    row++;
                }
            }
        }

        if (fitstable_fix_header(tab)) {
            ERROR("Failed to fix correspondence file header.");
            return -1;
        }

        fitstable_next_extension(tab);
        fitstable_clear_table(tab);
    }

    if (fitstable_close(tab)) {
        ERROR("Failed to close correspondence file");
        return -1;
    }

    return 0;
}

static int write_solutions(onefield_t* bp) {
    anbool got_solutions = (bl_size(bp->solutions) > 0);

    // If we found no solution, don't write empty output files!
    if (!got_solutions)
        return 0;

    // The solutions can fall out of order because tweak2() updates their logodds.
    bl_sort(bp->solutions, onefield_internal_compare_matchobjs);

    if (bp->matchfname) {
        if (write_match_file(bp))
            return -1;
    }
    if (bp->indexrdlsfname) {
        if (write_rdls_file(bp))
            return -1;
    }

    // We only want the best solution for each field in the following outputs:
    remove_duplicate_solutions(bp);

    if (bp->wcs_template) {
        if (write_wcs_file(bp))
            return -1;
    }
    if (bp->scamp_fname) {
        if (write_scamp_file(bp))
            return -1;
    }
    if (bp->corr_fname) {
        if (write_corr_file(bp))
            return -1;
    }
    return 0;
}

int onefield_internal_compare_matchobjs(const void* v1, const void* v2) {
    int diff;
    float fdiff;
    const MatchObj* mo1 = v1;
    const MatchObj* mo2 = v2;
    diff = mo1->fieldfile - mo2->fieldfile;
    if (diff) return diff;
    diff = mo1->fieldnum - mo2->fieldnum;
    if (diff) return diff;
    fdiff = mo1->logodds - mo2->logodds;
    if (fdiff == 0.0)
        return 0;
    if (fdiff > 0.0)
        return -1;
    return 1;
}
