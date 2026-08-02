/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <libgen.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "bl-sort.h"
#include "boilerplate.h"
#include "errors.h"
#include "fitsioutils.h"
#include "index.h"
#include "ioutils.h"
#include "log.h"
#include "mathutil.h"
#include "onefield_internal.h"
#include "os-features.h"
#include "permutedsort.h"
#include "tic.h"
#include "verify.h"

// SECTION INDEX-SHARD: bridge

static index_shard_hook_result_t onefield_index_shard_hook_result(
    index_shard_hook_outcome_t outcome,
    int error_code) {
  index_shard_hook_result_t result = {outcome, error_code};

  return result;
}

// ANCHOR INDEX-SHARD: bridge-get-index
static index_shard_hook_result_t onefield_index_shard_get_index(
    onefield_t *bp,
    size_t index_order,
    index_t **index_out) {
  index_t *index;

  if (index_out) {
    *index_out = NULL;
  }
  if (!bp || !index_out) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  if (index_order < (size_t)sl_size(bp->indexnames)) {
    const char *index_name = sl_get(bp->indexnames, index_order);
    size_t worker_stride =
        bp->index_shard_workers > 1
            ? (size_t)bp->index_shard_workers
            : 1U;

    index = onefield_internal_job_index_cache_get(
        bp,
        index_name);
    if (!index) {
      ERROR("Failed to load index %s", index_name);
      return onefield_index_shard_hook_result(
          INDEX_SHARD_HOOK_TASK_LOCAL_FAILURE,
          -1);
    }

    /*
     * Preparing one worker-width ahead overlaps mapping setup without claiming
     * that future index. The single prepared handoff transfers to the first
     * exact claimant and is freed when that worker finishes the index task.
     */
    if (index_order <= SIZE_MAX - worker_stride) {
      size_t prepare_order =
          index_order + worker_stride;

      if (prepare_order <
          (size_t)sl_size(bp->indexnames)) {
        onefield_internal_job_index_cache_prepare(
            bp,
            sl_get(bp->indexnames,
                   prepare_order));
      }
    }

    *index_out = index;
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_COMPLETED_UNSOLVED,
        0);
  }

  index_order -= (size_t)sl_size(bp->indexnames);
  if (index_order >= (size_t)pl_size(bp->indexes)) {
    ERROR("Index order %zu is outside the loaded index list", index_order);
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  index = pl_get(bp->indexes, index_order);
  if (!index) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  *index_out = index;
  return onefield_index_shard_hook_result(
      INDEX_SHARD_HOOK_COMPLETED_UNSOLVED,
      0);
}

// ANCHOR INDEX-SHARD: bridge-done-with-index
static index_shard_hook_result_t
onefield_index_shard_done_with_index(
    onefield_t *bp,
    size_t index_order,
    index_t *index) {
  if (!bp || !index) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }
  if (onefield_internal_done_with_index(bp, index_order, index)) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }
  return onefield_index_shard_hook_result(
      INDEX_SHARD_HOOK_COMPLETED_UNSOLVED,
      0);
}

// ANCHOR INDEX-SHARD: bridge-report-committed-solution
static int onefield_index_shard_report_committed_solution(
    onefield_t *bp,
    size_t index_order,
    int fieldnum,
    double best_logodds) {
  const char *index_name;
  char *index_base;

  if (!bp || fieldnum < 0) {
    logerr("[index-shard] invalid committed-solution metadata\n");
    return -1;
  }

  index_name = onefield_internal_get_index_name(bp, index_order);
  if (!index_name) {
    logerr("[index-shard] committed index order %zu has no filename\n",
           index_order);
    return -1;
  }

  index_base = basename_safe(index_name);
  if (!index_base) {
    SYSERROR("Failed to allocate committed index basename");
    return -1;
  }

  logmsg("Field %i: solved with index %s.\n", fieldnum, index_base);
  logverb("[index-shard] committed-solution index_order=%zu "
          "field=%i best_logodds=%.17g index_path=%s\n",
          index_order,
          fieldnum,
          best_logodds,
          index_name);

  free(index_base);
  return 0;
}

// ANCHOR INDEX-SHARD: bridge-prepare-shared-field
/*
 * Prepare the immutable field representation once on the pass owner. Worker
 * solvers borrow these pointers and release only their private task state.
 */
int onefield_index_shard_prepare_job_field_for_run(
    onefield_t *bp) {
  double field_read_seconds = 0.0;
  double preprocess_seconds = 0.0;
  int fieldnum;

  if (!bp || il_size(bp->fieldlist) != 1) {
    logerr("[index-shard] shared field preparation requires one field\n");
    return -1;
  }

  fieldnum = il_get(bp->fieldlist, 0);
  bp->fieldnum = fieldnum;

  if (onefield_internal_prepare_field_view(
          bp,
          fieldnum,
          &field_read_seconds,
          &preprocess_seconds)) {
    return -1;
  }
  if (!onefield_internal_field_cache_valid(bp)) {
    logverb("[index-shard] shared-field-cache state=unavailable "
            "reason=identity-key\n");
    return 1;
  }
  if (onefield_check_total_limits(bp)) {
    return 0;
  }

  if (bp->index_shard_workers > 1 &&
      !solver_prepare_field_geometry(&bp->solver)) {
    logverb("[index-shard] shared-field-geometry state=unavailable "
            "fallback=native\n");
  }

  logverb("[onefield] shared-field-cache=job-owned field=%i "
          "read=%.6f preprocess=%.6f\n",
          fieldnum,
          field_read_seconds,
          preprocess_seconds);

  return 0;
}

typedef struct onefield_index_shard_worker_view {
  solver_t solver;
  struct stat source_identity;
  char *fieldfname;
  char *indexrdlsfname;
  char *corr_fname;
  char *scamp_fname;
  char *solved_in;
  char *xcolname;
  char *ycolname;
  char *fieldid_key;
  char *sort_rdls;
  char *cancelfname;
  const sl *rdls_tagalong;
  double logratio_tosolve;
  int nsolves;
  int fieldnum;
  int fieldid;
  anbool xyls_tagalong_all;
} onefield_index_shard_worker_view_t;

/*
 * Initialize one worker from an allowlist of immutable configuration and
 * pass-bounded field views. Mutable solver state, output ownership, open file
 * handles, and task results are deliberately not cloned from the master.
 */
static void onefield_index_shard_initialize_local_solver(
    solver_t *local,
    const solver_t *base) {
  memset(local, 0, sizeof(*local));

  local->fieldxy = base->fieldxy;
  local->pixel_xscale = base->pixel_xscale;
  local->predistort = base->predistort;
  local->fieldxy_orig = base->fieldxy_orig;
  local->funits_lower = base->funits_lower;
  local->funits_upper = base->funits_upper;
  local->logratio_toprint = base->logratio_toprint;
  local->logratio_tokeep = base->logratio_tokeep;
  local->logratio_totune = base->logratio_totune;
  local->distance_from_quad_bonus = base->distance_from_quad_bonus;
  local->verify_uniformize = base->verify_uniformize;
  local->verify_dedup = base->verify_dedup;
  local->do_tweak = base->do_tweak;
  local->tweak_aborder = base->tweak_aborder;
  local->tweak_abporder = base->tweak_abporder;
  local->verify_pix = base->verify_pix;
  local->distractor_ratio = base->distractor_ratio;
  local->codetol = base->codetol;
  local->quadsize_min = base->quadsize_min;
  local->quadsize_max = base->quadsize_max;
  local->startobj = base->startobj;
  local->endobj = base->endobj;
  local->parity = base->parity;
  local->use_radec = base->use_radec;
  memcpy(local->centerxyz, base->centerxyz, sizeof(local->centerxyz));
  local->r2 = base->r2;
  local->logratio_bail_threshold = base->logratio_bail_threshold;
  local->logratio_stoplooking = base->logratio_stoplooking;
  local->maxquads = base->maxquads;
  local->maxmatches = base->maxmatches;
  local->set_crpix = base->set_crpix;
  local->set_crpix_center = base->set_crpix_center;
  memcpy(local->crpix, base->crpix, sizeof(local->crpix));

  local->minminAB2 = base->minminAB2;
  local->maxmaxAB2 = base->maxmaxAB2;
  local->rel_index_noise2 = base->rel_index_noise2;
  local->rel_field_noise2 = base->rel_field_noise2;
  local->abscale_low = base->abscale_low;
  local->abscale_high = base->abscale_high;
  local->field_minx = base->field_minx;
  local->field_maxx = base->field_maxx;
  local->field_miny = base->field_miny;
  local->field_maxy = base->field_maxy;
  local->field_diag = base->field_diag;
  local->cxdx_margin = base->cxdx_margin;
  local->vf = base->vf;
  local->field_geometry = base->field_geometry;
  local->field_geometry_owned = FALSE;
  local->index_mmap_policy = base->index_mmap_policy;

  solver_reset_counters(local);
  local->num_meanx_skipped = 0;
  solver_reset_best_match(local);
}

static int onefield_index_shard_duplicate_optional(
    char **destination,
    const char *source) {
  if (!destination) {
    return -1;
  }
  *destination = NULL;
  if (!source) {
    return 0;
  }
  *destination = strdup(source);
  return *destination ? 0 : -1;
}

static void onefield_index_shard_destroy_worker_view(
    void *opaque) {
  onefield_index_shard_worker_view_t *view = opaque;

  if (!view) {
    return;
  }
  free(view->fieldfname);
  free(view->indexrdlsfname);
  free(view->corr_fname);
  free(view->scamp_fname);
  free(view->solved_in);
  free(view->xcolname);
  free(view->ycolname);
  free(view->fieldid_key);
  free(view->sort_rdls);
  free(view->cancelfname);
  free(view);
}

static int onefield_index_shard_create_worker_view(
    onefield_t *master,
    const solver_t *base_solver,
    void **worker_view_out) {
  onefield_index_shard_worker_view_t *view;
  struct stat source_identity;
  int fieldnum;

  if (!worker_view_out) {
    return -1;
  }
  *worker_view_out = NULL;
  if (!master || !base_solver ||
      il_size(master->fieldlist) != 1 ||
      !base_solver->fieldxy_orig ||
      !base_solver->fieldxy ||
      !base_solver->vf ||
      master->rdls_tagalong_all ||
      master->xyls_tagalong ||
      !master->xyls_tagalong_all ||
      !master->fieldfname ||
      stat(master->fieldfname, &source_identity)) {
    return -1;
  }
  fieldnum = il_get_const(master->fieldlist, 0);
  if (!onefield_internal_field_cache_valid(master) ||
      !onefield_field_cache_key_matches(
          master, fieldnum, &source_identity)) {
    return -1;
  }

  view = calloc(1, sizeof(*view));
  if (!view) {
    return -1;
  }
  onefield_index_shard_initialize_local_solver(
      &view->solver, base_solver);
  view->source_identity = source_identity;
  view->rdls_tagalong = master->rdls_tagalong;
  view->logratio_tosolve = master->logratio_tosolve;
  view->nsolves = master->nsolves;
  view->fieldnum = fieldnum;
  view->fieldid = master->fieldid;
  view->xyls_tagalong_all = master->xyls_tagalong_all;

  if (onefield_index_shard_duplicate_optional(
          &view->fieldfname, master->fieldfname) ||
      onefield_index_shard_duplicate_optional(
          &view->indexrdlsfname, master->indexrdlsfname) ||
      onefield_index_shard_duplicate_optional(
          &view->corr_fname, master->corr_fname) ||
      onefield_index_shard_duplicate_optional(
          &view->scamp_fname, master->scamp_fname) ||
      onefield_index_shard_duplicate_optional(
          &view->solved_in, master->solved_in) ||
      onefield_index_shard_duplicate_optional(
          &view->xcolname, master->xcolname) ||
      onefield_index_shard_duplicate_optional(
          &view->ycolname, master->ycolname) ||
      onefield_index_shard_duplicate_optional(
          &view->fieldid_key, master->fieldid_key) ||
      onefield_index_shard_duplicate_optional(
          &view->sort_rdls, master->sort_rdls) ||
      onefield_index_shard_duplicate_optional(
          &view->cancelfname, master->cancelfname)) {
    onefield_index_shard_destroy_worker_view(view);
    return -1;
  }

  *worker_view_out = view;
  return 0;
}

static void onefield_index_shard_initialize_local_params(
    onefield_t *local,
    const onefield_index_shard_worker_view_t *view) {
  memset(local, 0, sizeof(*local));
  onefield_index_shard_initialize_local_solver(
      &local->solver, &view->solver);

  local->logratio_tosolve = view->logratio_tosolve;
  local->nsolves = view->nsolves;
  local->fieldfname = view->fieldfname;
  local->indexrdlsfname = view->indexrdlsfname;
  local->corr_fname = view->corr_fname;
  local->scamp_fname = view->scamp_fname;
  local->solved_in = view->solved_in;
  local->fieldnum = view->fieldnum;
  local->fieldid = view->fieldid;
  local->xcolname = view->xcolname;
  local->ycolname = view->ycolname;
  local->fieldid_key = view->fieldid_key;
  local->rdls_tagalong = (sl*)view->rdls_tagalong;
  local->rdls_tagalong_all = FALSE;
  local->sort_rdls = view->sort_rdls;
  local->xyls_tagalong = NULL;
  local->xyls_tagalong_all = view->xyls_tagalong_all;
  local->cancelfname = view->cancelfname;
}

// ANCHOR INDEX-SHARD: bridge-prepare-local-context
static int onefield_index_shard_prepare_local_context(onefield_t *local_bp,
                                                      const void *opaque) {
  const onefield_index_shard_worker_view_t *view = opaque;
  struct stat source_stat_after;
  int fieldnum;

  if (!local_bp || !view) {
    return -1;
  }
  onefield_index_shard_initialize_local_params(
      local_bp, view);

  local_bp->solver.indexes = pl_new(1);
  if (!local_bp->solver.indexes) {
    SYSERROR("Failed to allocate worker-local solver index list");
    return -1;
  }

  local_bp->solver.index = NULL;
  local_bp->solver.mo_template = NULL;
  local_bp->solver.record_match_callback = NULL;
  local_bp->solver.timer_callback = NULL;
  local_bp->solver.userdata = NULL;
  local_bp->solver.quit_now = FALSE;
  memset(&local_bp->solver.profile,
         0,
         sizeof(local_bp->solver.profile));

  solver_reset_counters(&local_bp->solver);
  solver_reset_best_match(&local_bp->solver);

  local_bp->solutions = NULL;
  local_bp->solved_out = NULL;
  local_bp->solved_fields_pending = NULL;

  local_bp->single_field_solved = FALSE;
  local_bp->solver_failed = FALSE;
  local_bp->nsolves_sofar = 0;

  local_bp->hit_cpulimit = FALSE;
  local_bp->hit_total_cpulimit = FALSE;
  local_bp->hit_timelimit = FALSE;
  local_bp->hit_total_timelimit = FALSE;
  local_bp->cancelled = FALSE;

  local_bp->cpulimit = 0.0;
  local_bp->total_cpulimit = 0.0;
  local_bp->timelimit = 0.0;
  local_bp->total_timelimit = 0.0;

  if (!view->solver.fieldxy_orig ||
      !view->solver.fieldxy ||
      !view->solver.vf) {
    logerr("[index-shard] shared field view is not prepared\n");
    goto fail;
  }

  fieldnum = view->fieldnum;
  local_bp->fieldnum = fieldnum;
  if (stat(view->fieldfname, &source_stat_after) ||
      !onefield_internal_same_source_identity(
          &view->source_identity, &source_stat_after)) {
    logerr("[index-shard] worker field source changed before "
           "local view preparation\n");
    goto fail;
  }

  /*
   * fieldxy_orig, fieldxy and vf were copied from the immutable pass view.
   * The master solver retains ownership until every worker is quiescent.
   */
  local_bp->xyls = xylist_open(view->fieldfname);
  if (!local_bp->xyls) {
    ERROR("Failed to open worker-local xylist %s", view->fieldfname);
    goto fail;
  }

  xylist_set_xname(local_bp->xyls, local_bp->xcolname);
  xylist_set_yname(local_bp->xyls, local_bp->ycolname);
  xylist_set_include_flux(local_bp->xyls, FALSE);
  xylist_set_include_background(local_bp->xyls, FALSE);

  if (xylist_open_field(local_bp->xyls, fieldnum)) {
    logerr("Failed to open extension %i in worker-local xylist.\n",
           fieldnum);
    goto fail;
  }
  if (stat(view->fieldfname, &source_stat_after) ||
      !onefield_internal_same_source_identity(
          &view->source_identity, &source_stat_after)) {
    logerr("[index-shard] worker field source changed during "
           "local view preparation\n");
    goto fail;
  }

  logverb("[index-shard] worker-field-view=borrowed field=%i\n",
          fieldnum);

  return 0;

fail:
  /* Never release master-owned field data from a worker error path. */
  local_bp->solver.fieldxy_orig = NULL;
  local_bp->solver.fieldxy = NULL;
  local_bp->solver.vf = NULL;
  local_bp->solver.field_geometry = NULL;
  local_bp->solver.field_geometry_owned = FALSE;
  solver_cleanup_field(&local_bp->solver);

  if (local_bp->xyls) {
    xylist_close(local_bp->xyls);
    local_bp->xyls = NULL;
  }

  if (local_bp->solver.indexes) {
    pl_free(local_bp->solver.indexes);
    local_bp->solver.indexes = NULL;
  }

  return -1;
}

// ANCHOR INDEX-SHARD: bridge-reset-local-context
static void
onefield_index_shard_reset_local_context_for_task(onefield_t *local_bp,
                                                  bl *local_solutions) {
  local_bp->solutions = local_solutions;

  local_bp->single_field_solved = FALSE;
  local_bp->solver_failed = FALSE;
  local_bp->nsolves_sofar = 0;

  local_bp->hit_cpulimit = FALSE;
  local_bp->hit_total_cpulimit = FALSE;
  local_bp->hit_timelimit = FALSE;
  local_bp->hit_total_timelimit = FALSE;
  local_bp->cancelled = FALSE;

  local_bp->solver.quit_now = FALSE;
  local_bp->solver.index = NULL;
  memset(&local_bp->solver.profile,
         0,
         sizeof(local_bp->solver.profile));

  solver_reset_counters(&local_bp->solver);
  local_bp->solver.num_meanx_skipped = 0;
  solver_reset_best_match(&local_bp->solver);

  solver_clear_indexes(&local_bp->solver);
}

// ANCHOR INDEX-SHARD: bridge-cleanup-local-context
static void onefield_index_shard_cleanup_local_context(onefield_t *local_bp) {
  if (!local_bp) {
    return;
  }

  local_bp->solver.mo_template = NULL;
  local_bp->solver.record_match_callback = NULL;
  local_bp->solver.timer_callback = NULL;
  local_bp->solver.userdata = NULL;

  solver_clear_indexes(&local_bp->solver);

  /* Field storage is owned by the master solver for the pass. */
  local_bp->solver.fieldxy_orig = NULL;
  local_bp->solver.fieldxy = NULL;
  local_bp->solver.vf = NULL;
  local_bp->solver.field_geometry = NULL;
  local_bp->solver.field_geometry_owned = FALSE;
  solver_cleanup_field(&local_bp->solver);

  if (local_bp->xyls) {
    xylist_close(local_bp->xyls);
    local_bp->xyls = NULL;
  }

  if (local_bp->solver.indexes) {
    pl_free(local_bp->solver.indexes);
    local_bp->solver.indexes = NULL;
  }

  local_bp->solutions = NULL;
}

/*
 * Run one index against the immutable field representation prepared once by
 * the pass owner. The owner-visible MatchObj template and every solver counter
 * remain task-local; only the XYLS field, star-list copy and verification
 * KD-tree are retained until the outer pass has quiesced.
 */
static int onefield_index_shard_solve_preprocessed_field(onefield_t *local_bp) {
  solver_t *sp = &local_bp->solver;
  MatchObj template;
  qfits_header *fieldhdr;
  double field_wall_start;
  anbool interrupted_by_parallel_stop = FALSE;
  int fieldnum = local_bp->fieldnum;

  if (!sp->fieldxy_orig || !sp->fieldxy || !sp->vf) {
    logerr("[index-shard] shared field view is not prepared\n");
    sp->profile.execution_failed = TRUE;
    local_bp->solver_failed = TRUE;
    return -1;
  }

  memset(&template, 0, sizeof(MatchObj));
  template.fieldnum = fieldnum;
  template.fieldfile = local_bp->fieldid;

  fieldhdr = xylist_get_header(local_bp->xyls);
  if (fieldhdr) {
    char *idstr = fits_get_dupstring(fieldhdr, local_bp->fieldid_key);

    if (idstr) {
      strncpy(template.fieldname,
              idstr,
              sizeof(template.fieldname) - 1);
    }
    free(idstr);
  }

  sp->mo_template = &template;
  sp->record_match_callback = onefield_internal_record_match_callback;
  sp->timer_callback = onefield_internal_timer_callback;
  sp->userdata = local_bp;
  sp->distance_from_quad_bonus = TRUE;

  local_bp->nsolves_sofar = 0;
  field_wall_start = monotonic_seconds();

  logverb("Solving field %i.\n", fieldnum);
  solver_log_params(sp);

  if (solver_run(sp)) {
    local_bp->solver_failed = TRUE;
  }

  /*
   * A losing owner can observe the pass stop while unwinding from solver_run.
   * Its partial traversal is not a completed scientific "did not solve"
   * result. Preserve the ordinary per-index report for completed work and
   * local user/limit cancellation, but omit this one misleading observation.
   */
  interrupted_by_parallel_stop =
      index_shard_worker_stop_requested() &&
      !local_bp->cancelled &&
      !local_bp->hit_total_timelimit &&
      !local_bp->hit_total_cpulimit;

  sp->mo_template = NULL;
  sp->record_match_callback = NULL;
  sp->timer_callback = NULL;
  sp->userdata = NULL;

  logverb("[onefield-field-profile] field=%i read=0.000000 "
          "preprocess=0.000000 solver_run=%.6f total=%.6f "
          "field_view=job-borrowed failed=%i\n",
          fieldnum,
          sp->profile.solver_run_wall_seconds,
          monotonic_seconds() - field_wall_start,
          sp->profile.execution_failed ? 1 : 0);

  if (local_bp->solver_failed || sp->profile.execution_failed) {
    local_bp->solver_failed = TRUE;
    logerr("Solver execution failed for field %i\n", fieldnum);
    return -1;
  }

  logverb("Field %i: tried %i quads, matched %i codes.\n",
          fieldnum,
          sp->numtries,
          sp->nummatches);

  if (sp->maxquads && sp->numtries >= sp->maxquads) {
    logmsg("  exceeded the number of quads to try: %i >= %i.\n",
           sp->numtries,
           sp->maxquads);
  }
  if (sp->maxmatches && sp->nummatches >= sp->maxmatches) {
    logmsg("  exceeded the number of quads to match: %i >= %i.\n",
           sp->nummatches,
           sp->maxmatches);
  }
  if (local_bp->cancelled) {
    logmsg("  cancelled at user request.\n");
  }

  if (sp->best_match_solves) {
    local_bp->single_field_solved = TRUE;
  } else if (!interrupted_by_parallel_stop &&
             sp->index && sp->index->indexname) {
    char *copy = strdup_safe(sp->index->indexname);
    char *base = basename(copy);

    if (sp->endobj) {
      logerr("Field %i did not solve (index %s, field objects %i-%i).\n",
             fieldnum,
             base,
             sp->startobj + 1,
             sp->endobj);
    } else {
      logerr("Field %i did not solve (index %s).\n",
             fieldnum,
             base);
    }
    free(copy);

    if (sp->have_best_match) {
      logverb("Best match encountered: ");
      matchobj_print(&sp->best_match, log_get_level());
    } else {
      logverb("Best odds encountered: %g\n", exp(sp->best_logodds));
    }
  } else if (!interrupted_by_parallel_stop) {
    logerr("Field %i did not solve.\n", fieldnum);
  }

  return 0;
}

// ANCHOR INDEX-SHARD: bridge-solve-one-index
static index_shard_hook_result_t
onefield_index_shard_solve_one_index(
    onefield_t *local_bp,
    index_t *index) {
  index_shard_hook_result_t hook_result;
  int rc;

  if (!local_bp || !index) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  solver_add_index(&local_bp->solver, index);

  local_bp->cpu_start = get_cpu_usage();
  local_bp->time_start = monotonic_seconds();

  rc = onefield_index_shard_solve_preprocessed_field(local_bp);
  if (rc || local_bp->solver_failed ||
      local_bp->solver.profile.execution_failed) {
    hook_result = onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        rc ? rc : -1);
  } else if (local_bp->cancelled) {
    hook_result = onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_CANCELLED,
        0);
  } else if (local_bp->hit_total_timelimit) {
    hook_result = onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_WALL_LIMIT,
        0);
  } else if (local_bp->hit_total_cpulimit) {
    hook_result = onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_CPU_LIMIT,
        0);
  } else {
    hook_result = onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_COMPLETED_UNSOLVED,
        0);
  }

  solver_clear_indexes(&local_bp->solver);
  return hook_result;
}

// ANCHOR INDEX-SHARD: bridge-analyze-solutions
static index_shard_hook_result_t
onefield_index_shard_analyze_solutions(
    onefield_t *master_bp,
    bl *solutions,
    double *best_logodds,
    int *best_fieldnum) {
  int required_solutions;
  int solution_count = 0;
  int i;

  if (best_logodds) {
    *best_logodds = -HUGE_VAL;
  }

  if (best_fieldnum) {
    *best_fieldnum = -1;
  }

  if (!master_bp || !solutions) {
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  if (il_size(master_bp->fieldlist) != 1) {
    logerr("[index-shard] refusing non-single-field solution analysis\n");
    return onefield_index_shard_hook_result(
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE,
        -1);
  }

  required_solutions = MAX(1, master_bp->nsolves);

  for (i = 0; i < bl_size(solutions); i++) {
    MatchObj *mo = bl_access(solutions, i);

    if (mo->logodds >= master_bp->logratio_tosolve) {
      solution_count++;
    }

    if (best_logodds && mo->logodds > *best_logodds) {
      *best_logodds = mo->logodds;

      if (best_fieldnum) {
        *best_fieldnum = mo->fieldnum;
      }
    }
  }

  return onefield_index_shard_hook_result(
      solution_count >= required_solutions
          ? INDEX_SHARD_HOOK_COMPLETED_SOLVED
          : INDEX_SHARD_HOOK_COMPLETED_UNSOLVED,
      0);
}

// ANCHOR INDEX-SHARD: bridge-disown-matchobj
static void onefield_index_shard_disown_matchobj(MatchObj *mo) {
  if (!mo)
    return;

  mo->sip = NULL;
  mo->refradec = NULL;
  mo->fieldxy = NULL;
  mo->theta = NULL;
  mo->matchodds = NULL;
  mo->refxyz = NULL;
  mo->refxy = NULL;
  mo->refstarid = NULL;
  mo->testperm = NULL;
  mo->tagalong = NULL;
  mo->field_tagalong = NULL;
}

// ANCHOR INDEX-SHARD: bridge-merge-solutions
static int onefield_index_shard_merge_solutions(onefield_t *master_bp,
                                                bl *solutions,
                                                anbool *solved_out) {
  int required_solutions;
  int solution_count = 0;
  int i;
  anbool solved = FALSE;

  if (solved_out) {
    *solved_out = FALSE;
  }

  if (!master_bp || !solutions) {
    return 0;
  }

  if (il_size(master_bp->fieldlist) != 1) {
    logerr("[index-shard] refusing non-single-field solution merge\n");
    return -1;
  }

  required_solutions = MAX(1, master_bp->nsolves);

  for (i = 0; i < bl_size(solutions); i++) {
    MatchObj *src = bl_access(solutions, i);

    bl_insert_sorted(master_bp->solutions, src, onefield_internal_compare_matchobjs);

    if (src->logodds >= master_bp->logratio_tosolve) {
      solution_count++;
    }

    if (solution_count == required_solutions) {
      /*
       * The serial callback declares a field solved on exactly the Nth
       * above-threshold match. Preserve that nsolves contract per field at the
       * authoritative master commit point. The worker already emitted the
       * chronological MatchObj diagnostics; do not print a sorted-list entry
       * here and misrepresent it as the Nth chronological hit.
       */
      onefield_internal_solved_field(master_bp, src->fieldnum);
      solved = TRUE;
    }

    onefield_index_shard_disown_matchobj(src);
  }

  bl_remove_all(solutions);

  if (solved_out) {
    *solved_out = solved;
  }

  return 0;
}

// ANCHOR INDEX-SHARD: bridge-free-solutions
static void onefield_index_shard_free_solutions(bl *solutions) {
  int i;

  if (!solutions)
    return;

  for (i = 0; i < bl_size(solutions); i++) {
    MatchObj *mo = bl_access(solutions, i);
    verify_free_matchobj(mo);
    onefield_free_matchobj(mo);
  }

  bl_free(solutions);
}

// ANCHOR INDEX-SHARD: bridge-hooks
static const index_shard_hooks_t onefield_index_shard_hooks = {
    .get_index = onefield_index_shard_get_index,
    .done_with_index = onefield_index_shard_done_with_index,
    .report_committed_solution =
        onefield_index_shard_report_committed_solution,

    .create_worker_view =
        onefield_index_shard_create_worker_view,
    .destroy_worker_view =
        onefield_index_shard_destroy_worker_view,
    .prepare_local_context =
        onefield_index_shard_prepare_local_context,
    .reset_local_context_for_task =
        onefield_index_shard_reset_local_context_for_task,
    .cleanup_local_context =
        onefield_index_shard_cleanup_local_context,

    .solve_one_index = onefield_index_shard_solve_one_index,
    .analyze_solutions =
        onefield_index_shard_analyze_solutions,
    .merge_solutions = onefield_index_shard_merge_solutions,
    .free_solutions = onefield_index_shard_free_solutions};



index_shard_solve_status_t onefield_index_shard_solve(
    onefield_t* bp,
    solver_t* solver,
    size_t index_count) {
  return index_shard_solve(
      bp,
      solver,
      index_count,
      &onefield_index_shard_hooks);
}
