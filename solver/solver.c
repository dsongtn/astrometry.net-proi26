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
#include "../libkd/kdtree_prefetch_internal.h"
#include "tweak2.h"
#include "astrometry/fitsbin.h"
#include "index_shard_internal.h"

#ifndef SOLVER_FIELD_GEOMETRY_BUDGET_BYTES
#define SOLVER_FIELD_GEOMETRY_BUDGET_BYTES \
    (64ULL * 1024ULL * 1024ULL)
#endif

#define SOLVER_CODEKD_SEARCH_OPTIONS \
    (KD_OPTIONS_SMALL_RADIUS | \
     KD_OPTIONS_COMPUTE_DISTS | \
     KD_OPTIONS_NO_RESIZE_RESULTS | \
     KD_OPTIONS_USE_SPLIT)

#define SOLVER_STARKD_VERIFY_SEARCH_OPTIONS \
    (KD_OPTIONS_SMALL_RADIUS | KD_OPTIONS_RETURN_POINTS)

static kdtree_qres_t* solver_codekd_rangesearch(
    const kdtree_t* tree,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options) {
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
    int* theta = NULL;
    double* odds = NULL;
    sip_t* tuned_sip = NULL;
    double* refradec;
    int i;
    double newodds;
    int nm, nc, nd;
    int besti;
    int startorder;

    if (!sp || !mo || !sp->fieldxy ||
        !mo->refxyz || mo->nindex <= 0 ||
        (size_t)mo->nindex > SIZE_MAX / (2U * sizeof(double))) {
        return;
    }

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
    refradec = malloc(2 * mo->nindex * sizeof(double));
    if (!xy || !refradec) {
        free(refradec);
        free(xy);
        logverb("solver_tweak2: allocation failed; "
                "preserving verified TAN solution\n");
        return;
    }
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
    besti = mo->nbest-1;//mo->nmatch + mo->nconflict + mo->ndistractor;

    logverb("solver_tweak2: set_crpix %i, crpix (%.1f,%.1f)\n",
            sp->set_crpix, sp->crpix[0], sp->crpix[1]);
    tuned_sip = tweak2(xy, Nxy,
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
    free(xy);
    xy = NULL;

    if (!tuned_sip || !theta || !odds ||
        besti < 0 || besti >= Nxy || !isfinite(newodds)) {
        sip_free(tuned_sip);
        free(theta);
        free(odds);
        logverb("solver_tweak2: tune failed; "
                "preserving verified TAN solution\n");
        return;
    }

    // Commit the tuned result only after every output is complete.
    if (mo->sip) {
        sip_free(mo->sip);
    }
    mo->sip = tuned_sip;
    free(mo->refxy);
    mo->refxy = NULL;
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
     * Shard workers apply RANDOM only to sparse payload chunks; compact
     * topology and serial solving retain original NORMAL mappings.
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

static void solver_begin_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity);

typedef struct solver_candidate_delivery_record
    solver_candidate_delivery_record_t;

static void solver_execute_prepared_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    kdtree_qres_t* result,
    double search_wall_seconds,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_quad_count,
    size_t prepared_star_count);

static void solver_retire_codekd_failure_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    int search_errno,
    double search_wall_seconds);

static void resolve_matches_with_delivery(
    kdtree_qres_t* krez,
    const double* field,
    const int* fstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count);

static void resolve_matches_native_range(
    kdtree_qres_t* krez,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count);

static int solver_handle_hit(solver_t* sp,
                             MatchObj* mo,
                             sip_t* sip,
                             anbool fake_match);

static int solver_handle_query_hit(
    solver_t* sp,
    MatchObj* mo,
    sip_t* sip,
    anbool fake_match,
    verify_index_query_t** query);

static int solver_handle_scored_query_hit(
    solver_t* sp,
    MatchObj* mo,
    sip_t* sip,
    anbool fake_match,
    solver_candidate_delivery_record_t* record);

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
    total->verification_score_wall_seconds +=
        profile->verification_score_wall_seconds;

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
    total->ab_helper_tasks += profile->ab_helper_tasks;
    total->ab_helper_combinations +=
        profile->ab_helper_combinations;
    total->page_plan_descriptors_total +=
        profile->page_plan_descriptors_total;
    total->page_plan_descriptors_complete +=
        profile->page_plan_descriptors_complete;
    total->page_plan_descriptor_splits +=
        profile->page_plan_descriptor_splits;
    total->page_plan_raw_ranges += profile->page_plan_raw_ranges;
    total->page_plan_unique_pages += profile->page_plan_unique_pages;
    total->page_plan_ranges_after_dedup +=
        profile->page_plan_ranges_after_dedup;
    total->page_plan_logical_bytes +=
        profile->page_plan_logical_bytes;
    total->page_plan_aligned_bytes +=
        profile->page_plan_aligned_bytes;
    total->page_plan_overread_bytes +=
        profile->page_plan_overread_bytes;
    total->page_plan_not_applicable +=
        profile->page_plan_not_applicable;
    total->page_plan_allocation_refused +=
        profile->page_plan_allocation_refused;
    total->page_plan_source_mismatch +=
        profile->page_plan_source_mismatch;
    total->page_plan_invalid_range +=
        profile->page_plan_invalid_range;
    total->page_plan_byte_budget_refused +=
        profile->page_plan_byte_budget_refused;
    total->page_plan_range_capacity_refused +=
        profile->page_plan_range_capacity_refused;
    total->page_plan_service_refused +=
        profile->page_plan_service_refused;
    total->page_plan_service_errors +=
        profile->page_plan_service_errors;
    total->page_plan_cancelled += profile->page_plan_cancelled;
    total->candidate_delivery_candidates +=
        profile->candidate_delivery_candidates;
    total->candidate_quad_submitted +=
        profile->candidate_quad_submitted;
    total->candidate_quad_ready +=
        profile->candidate_quad_ready;
    total->candidate_quad_fallback +=
        profile->candidate_quad_fallback;
    total->candidate_star_submitted +=
        profile->candidate_star_submitted;
    total->candidate_star_ready +=
        profile->candidate_star_ready;
    total->candidate_star_fallback +=
        profile->candidate_star_fallback;
    total->candidate_delivery_windows +=
        profile->candidate_delivery_windows;
    total->candidate_quad_ready_rows +=
        profile->candidate_quad_ready_rows;
    total->candidate_star_ready_rows +=
        profile->candidate_star_ready_rows;
    total->candidate_retired_rows +=
        profile->candidate_retired_rows;
    total->candidate_native_rows +=
        profile->candidate_native_rows;
    total->verification_page_queries +=
        profile->verification_page_queries;
    total->verification_page_queries_planned +=
        profile->verification_page_queries_planned;
    total->verification_page_prefixes +=
        profile->verification_page_prefixes;
    total->verification_page_submitted +=
        profile->verification_page_submitted;
    total->verification_page_ready +=
        profile->verification_page_ready;
    total->verification_page_fallback +=
        profile->verification_page_fallback;
    total->verification_page_ready_rows +=
        profile->verification_page_ready_rows;
    total->verification_page_ranges +=
        profile->verification_page_ranges;
    total->verification_page_logical_bytes +=
        profile->verification_page_logical_bytes;
    total->verification_page_aligned_bytes +=
        profile->verification_page_aligned_bytes;
    total->candidate_math_prepared +=
        profile->candidate_math_prepared;
    total->candidate_math_reused +=
        profile->candidate_math_reused;
    total->verification_score_batches_prepared +=
        profile->verification_score_batches_prepared;
    total->verification_score_contexts_prepared +=
        profile->verification_score_contexts_prepared;
    total->verification_score_batches_executed +=
        profile->verification_score_batches_executed;
    total->verification_score_contexts_completed +=
        profile->verification_score_contexts_completed;
    if (ULLONG_MAX -
            total->verification_score_work_units_completed <
        profile->verification_score_work_units_completed) {
        total->verification_score_work_units_completed =
            ULLONG_MAX;
    } else {
        total->verification_score_work_units_completed +=
            profile->verification_score_work_units_completed;
    }
    total->verification_score_fallback_batches +=
        profile->verification_score_fallback_batches;
    total->verification_score_stopped_batches +=
        profile->verification_score_stopped_batches;
    total->staged_owner_claims += profile->staged_owner_claims;
    total->staged_foreign_claims += profile->staged_foreign_claims;
    total->staged_io_submitted += profile->staged_io_submitted;
    total->staged_io_completed += profile->staged_io_completed;
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
    if (profile->max_staged_io_submitted >
        total->max_staged_io_submitted) {
        total->max_staged_io_submitted =
            profile->max_staged_io_submitted;
    }
    if (profile->max_staged_compute_ready >
        total->max_staged_compute_ready) {
        total->max_staged_compute_ready =
            profile->max_staged_compute_ready;
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
           "verify_calls=%llu "
           "hypothesis_batches=%llu parallel_batches=%llu "
           "parallel_batches_observed=%llu "
           "parallel_hypotheses=%llu "
           "task_ranges_planned=%llu task_ranges_executed=%llu "
           "task_ranges_inline=%llu allocation_failures=%llu "
           "search_failures=%llu "
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
           profile->hypothesis_batches,
           profile->parallel_batches,
           profile->parallel_batches_observed,
           profile->parallel_hypotheses,
           profile->task_ranges_planned,
           profile->task_ranges_executed,
           profile->task_ranges_inline,
           profile->allocation_failures,
           profile->search_failures,
           profile->ab_helper_tasks,
           profile->ab_helper_combinations,
           profile->hypothesis_order_hash,
           profile->kd_result_order_hash,
           profile->candidate_order_hash);
    logmsg("[solver] page-pipeline descriptors=%llu complete=%llu "
           "boundary_deferrals=%llu raw_hints=%llu unique_pages=%llu "
           "coalesced_ranges=%llu raw_hint_bytes=%llu aligned_bytes=%llu "
           "positive_alignment_delta_bytes=%llu "
           "refusals=%llu/%llu/%llu/%llu/%llu/"
           "%llu/%llu/%llu/%llu staged_claims=%llu/%llu "
           "io=%llu/%llu max_io=%zu max_ready=%zu "
           "candidate_delivery=%llu quad=%llu/%llu/%llu "
           "star=%llu/%llu/%llu windows=%llu "
           "rows=%llu/%llu/%llu/%llu "
           "verify_pages=%llu/%llu prefixes=%llu tickets=%llu/%llu "
           "fallback=%llu ready_rows=%llu ranges=%llu "
           "bytes=%llu/%llu candidate_math=%llu/%llu "
           "verify_score=%llu/%llu/%llu/%llu/%llu/%llu "
           "work=%llu work_wall_sum=%.6f\n",
           profile->page_plan_descriptors_total,
           profile->page_plan_descriptors_complete,
           profile->page_plan_descriptor_splits,
           profile->page_plan_raw_ranges,
           profile->page_plan_unique_pages,
           profile->page_plan_ranges_after_dedup,
           profile->page_plan_logical_bytes,
           profile->page_plan_aligned_bytes,
           profile->page_plan_overread_bytes,
           profile->page_plan_not_applicable,
           profile->page_plan_allocation_refused,
           profile->page_plan_source_mismatch,
           profile->page_plan_invalid_range,
           profile->page_plan_byte_budget_refused,
           profile->page_plan_range_capacity_refused,
           profile->page_plan_service_refused,
           profile->page_plan_service_errors,
           profile->page_plan_cancelled,
           profile->staged_owner_claims,
           profile->staged_foreign_claims,
           profile->staged_io_submitted,
           profile->staged_io_completed,
           profile->max_staged_io_submitted,
           profile->max_staged_compute_ready,
           profile->candidate_delivery_candidates,
           profile->candidate_quad_submitted,
           profile->candidate_quad_ready,
           profile->candidate_quad_fallback,
           profile->candidate_star_submitted,
           profile->candidate_star_ready,
           profile->candidate_star_fallback,
           profile->candidate_delivery_windows,
           profile->candidate_quad_ready_rows,
           profile->candidate_star_ready_rows,
           profile->candidate_retired_rows,
           profile->candidate_native_rows,
           profile->verification_page_queries,
           profile->verification_page_queries_planned,
           profile->verification_page_prefixes,
           profile->verification_page_submitted,
           profile->verification_page_ready,
           profile->verification_page_fallback,
           profile->verification_page_ready_rows,
           profile->verification_page_ranges,
           profile->verification_page_logical_bytes,
           profile->verification_page_aligned_bytes,
           profile->candidate_math_prepared,
           profile->candidate_math_reused,
           profile->verification_score_batches_prepared,
           profile->verification_score_contexts_prepared,
           profile->verification_score_batches_executed,
           profile->verification_score_contexts_completed,
           profile->verification_score_fallback_batches,
           profile->verification_score_stopped_batches,
           profile->verification_score_work_units_completed,
           profile->verification_score_wall_seconds);
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

#ifndef SOLVER_AB_CANDIDATE_LIMIT_BYTES
#define SOLVER_AB_CANDIDATE_LIMIT_BYTES (768U * 1024U)
#endif
#define SOLVER_AB_DESCRIPTOR_CAPACITY 4096U
/* Keep task splitting coarse, but admit one-task delivery by actual work. */
#define SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS 64U
#define SOLVER_AB_DESCRIPTOR_DELIVERY_MIN_HYPOTHESES 64U
#define SOLVER_AB_DESCRIPTOR_TASK_MIN_HYPOTHESES 32U
#define SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS 1000
#define SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES \
    (8U * 1024U * 1024U)

typedef enum solver_ab_phase_kind {
    SOLVER_AB_PHASE_DIAGONAL = 0,
    SOLVER_AB_PHASE_OFF_DIAGONAL = 1
} solver_ab_phase_kind_t;

typedef enum solver_ab_phase_mode {
    SOLVER_AB_MODE_NATIVE = 0,
    SOLVER_AB_MODE_FLATTENED_OWNER = 1,
    SOLVER_AB_MODE_ASSISTED = 2
} solver_ab_phase_mode_t;

typedef enum solver_ab_candidate_action {
    SOLVER_AB_CANDIDATE_RADEC_SKIP = 0,
    SOLVER_AB_CANDIDATE_ABSCALE_SKIP = 1,
    SOLVER_AB_CANDIDATE_BAD_QUAD = 2,
    SOLVER_AB_CANDIDATE_SCALE_SKIP = 3,
    SOLVER_AB_CANDIDATE_VERIFY = 4
} solver_ab_candidate_action_t;

typedef enum solver_packet_reserve_result {
    SOLVER_PACKET_RESERVE_ERROR = -1,
    SOLVER_PACKET_RESERVE_OK = 0,
    SOLVER_PACKET_RESERVE_FULL = 1
} solver_packet_reserve_result_t;

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

typedef struct solver_verification_packet {
    solver_ab_candidate_t* candidates;
    size_t candidate_capacity;
    anbool allocation_failed;
} solver_verification_packet_t;

typedef struct solver_descriptor_status {
    anbool evaluation_failed;
    anbool cancelled;
} solver_descriptor_status_t;

typedef struct solver_ab_task {
    unsigned long long combination_first;
    unsigned long long combination_end;
} solver_ab_task_t;

typedef struct solver_ab_descriptor_planner
    solver_ab_descriptor_planner_t;

typedef struct solver_ab_builder {
    solver_ab_descriptor_planner_t* planner;
    solver_descriptor_status_t* status;
    double tol2;
    double rel_field_noise2;
    unsigned long long pending_numtries;
    unsigned long long pending_cxdx;
    unsigned long long pending_meanx;
    solver_ab_descriptor_output_t* descriptor_output;
    anbool rel_field_noise_valid;
    anbool fatal_error;
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
} solver_ab_phase_telemetry_t;

struct solver_ab_descriptor_planner {
    solver_ab_snapshot_t snapshot;
    const solver_field_geometry_t* field_geometry;
    int newpoint;
    int dimquads;
    solver_ab_phase_kind_t phase;

    solver_ab_pair_t* pairs;
    size_t pair_count;
    size_t pair_capacity;
    solver_ab_pair_t* pair_cache;
    size_t pair_cache_capacity;
    size_t pair_cache_limit_bytes;
    anbool pairs_transient;

    int* combination_eligible;
    size_t combination_eligible_capacity;
};

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
    double code_err) {
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
        solver_order_double_bits(code_err));
    solver->profile.candidate_order_hash =
        solver_order_hash_mix(
            solver->profile.candidate_order_hash,
            digest);
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

/*
 * Descriptor planning is index-free. A synchronous packet may borrow a const
 * CodeKD tree view; QuadFile, StarKD, lazy index initialization, verification,
 * counters, callbacks, and reducer state remain owner-local.
 */
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
                    : "native"),
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

static solver_packet_reserve_result_t
solver_verification_packet_reserve(
    solver_verification_packet_t* packet,
    size_t required) {
    solver_ab_candidate_t* resized;
    size_t capacity;
    size_t max_capacity =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES / sizeof(*resized);
    size_t bytes;

    if (required <= packet->candidate_capacity) {
        return SOLVER_PACKET_RESERVE_OK;
    }
    if (required > max_capacity || !max_capacity) {
        return SOLVER_PACKET_RESERVE_FULL;
    }
    capacity = packet->candidate_capacity ?
        packet->candidate_capacity : MIN(16U, max_capacity);
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            packet->allocation_failed = TRUE;
            return SOLVER_PACKET_RESERVE_ERROR;
        }
        if (capacity > max_capacity / 2U) {
            capacity = max_capacity;
        } else {
            capacity *= 2U;
        }
    }
    if (capacity > SIZE_MAX / sizeof(*resized)) {
        packet->allocation_failed = TRUE;
        return SOLVER_PACKET_RESERVE_ERROR;
    }
    bytes = capacity * sizeof(*resized);
    if (bytes > SOLVER_AB_CANDIDATE_LIMIT_BYTES) {
        return SOLVER_PACKET_RESERVE_FULL;
    }

    resized = realloc(packet->candidates, bytes);
    if (!resized) {
        packet->allocation_failed = TRUE;
        return SOLVER_PACKET_RESERVE_ERROR;
    }
    packet->candidates = resized;
    packet->candidate_capacity = capacity;
    return SOLVER_PACKET_RESERVE_OK;
}

static void solver_verification_packet_free(solver_verification_packet_t* packet) {
    if (!packet) {
        return;
    }
    free(packet->candidates);
    memset(packet, 0, sizeof(*packet));
}

static anbool solver_ab_cancelled(
    const solver_ab_descriptor_planner_t* executor) {
    (void)executor;
    return index_shard_worker_stop_requested();
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
        if (builder->status) {
            builder->status->evaluation_failed = TRUE;
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
    if (!builder || !builder->status ||
        index_shard_worker_stop_requested()) {
        if (builder && builder->status) {
            builder->status->cancelled = TRUE;
        }
        return -1;
    }
    return solver_ab_record_descriptor(
        builder,
        stars,
        code,
        dimquad,
        current_parity);
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
        &builder->planner->snapshot;
    double mycode[DCMAX];
    int nstars = dimquad - NBACK;
    int lastslot = dimquad - NBACK - 1;
    int i;

    if (!code) {
        code = mycode;
    }
    if (slot >= DCMAX / 2 ||
        solver_ab_cancelled(builder->planner) ||
        builder->fatal_error) {
        if (solver_ab_cancelled(builder->planner)) {
            builder->status->cancelled = TRUE;
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
                builder->status->cancelled) {
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
        builder->status->cancelled) {
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
        builder->planner->field_geometry;
    const solver_ab_snapshot_t* snapshot =
        &builder->planner->snapshot;
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
            builder->status->evaluation_failed = TRUE;
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
        builder->status->cancelled) {
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

static int* solver_ab_get_eligible_workspace(
    solver_ab_descriptor_planner_t* planner,
    size_t required) {
    if (!planner || !planner->combination_eligible ||
        required > planner->combination_eligible_capacity) {
        return NULL;
    }
    return planner->combination_eligible;
}

static int solver_ab_visit_pair_range(
    solver_ab_descriptor_planner_t* executor,
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

    if (solver_ab_cancelled(builder->planner)) {
        builder->status->cancelled = TRUE;
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
        builder->status->cancelled;
}

typedef int (*solver_ab_pair_range_visitor_t)(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque);

static int solver_ab_walk_task_ranges(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_task_t* task,
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
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque) {
    solver_ab_builder_t* builder = opaque;

    if (solver_ab_cancelled(executor)) {
        builder->status->cancelled = TRUE;
        return 1;
    }
    builder->tol2 = pair->tol2;
    builder->rel_field_noise2 =
        pair_geometry->rel_field_noise2;
    builder->rel_field_noise_valid = TRUE;
    return solver_ab_visit_pair_range(
        executor,
        pair,
        pair_geometry,
        local_first,
        local_end,
        solver_ab_builder_visit,
        builder);
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

static int solver_ab_reserve_pair_buffer(
    solver_ab_descriptor_planner_t* executor,
    size_t required) {
    solver_ab_pair_t* resized;
    size_t capacity;
    size_t retained_maximum;

    if (!executor ||
        required > SIZE_MAX / sizeof(*resized)) {
        return -1;
    }
    if (!executor->pair_cache_limit_bytes) {
        executor->pair_cache_limit_bytes =
            SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES;
    }
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
    return 0;
}

static void solver_ab_begin_pair_buffer(
    solver_ab_descriptor_planner_t* executor) {
    executor->pairs = executor->pair_cache;
    executor->pair_capacity =
        executor->pair_cache_capacity;
    executor->pairs_transient = FALSE;
}

static int solver_ab_collect_pairs(
    solver_ab_descriptor_planner_t* executor,
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

#define SOLVER_PAYLOAD_CANDIDATE_BATCH 32U
#define SOLVER_VERIFICATION_WINDOW_CANDIDATES \
    (2U * SOLVER_PAYLOAD_CANDIDATE_BATCH)

/*
 * Inner compute assistance never starts or waits for a cold page provider.
 * Only an already resident Quad/Star payload may enter prepared verification;
 * every other case remains on the native owner-local mmap path.
 */
static anbool solver_payload_candidate_data_fully_resident(
    const solver_t* solver) {
    if (!solver || !solver->index || !solver->index->quads ||
        !solver->index->quads->fb || !solver->index->starkd ||
        !solver->index->starkd->tree ||
        !solver->index->starkd->tree->io ||
        !solver->index->starkd->tree->io_is_fitsbin) {
        return FALSE;
    }
    return fitsbin_payload_is_fully_resident(
               solver->index->quads->fb) &&
        fitsbin_payload_is_fully_resident(
               (const fitsbin_t*)solver->index->starkd->tree->io);
}

typedef struct solver_verification_score_slot {
    verify_prepared_hit_t* prepared;
} solver_verification_score_slot_t;

typedef struct solver_verification_candidate_runtime {
    MatchObj match;
    double match_distance_in_pixels2;
    double logaccept;
} solver_verification_candidate_runtime_t;

typedef struct solver_verification_task_input {
    const solver_verification_score_slot_t* slots;
    size_t slot_first;
    size_t slot_count;
} solver_verification_task_input_t;

typedef struct solver_verification_retire_context {
    size_t next_slot;
} solver_verification_retire_context_t;

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

static index_shard_helper_retire_status_t
solver_verification_packet_retire(
    const index_shard_helper_task_t* task,
    size_t task_index,
    void* owner_context) {
    const solver_verification_task_input_t* input;
    solver_verification_retire_context_t* context = owner_context;

    (void)task_index;
    if (!task || !context || !task->input ||
        task->input_bytes != sizeof(*input)) {
        return INDEX_SHARD_HELPER_RETIRE_ERROR;
    }
    input = task->input;
    if (!input->slot_count ||
        input->slot_first != context->next_slot ||
        input->slot_count > SIZE_MAX - context->next_slot) {
        return INDEX_SHARD_HELPER_RETIRE_ERROR;
    }
    context->next_slot += input->slot_count;
    return INDEX_SHARD_HELPER_RETIRE_OK;
}

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
    solver_verification_packet_t* packet) {
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
    solver_verification_packet_free(packet);
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
    int candidate_first,
    int candidate_end) {
    solver_verification_packet_t packet;
    solver_ab_snapshot_t snapshot;
    solver_verification_candidate_runtime_t* runtime = NULL;
    solver_verification_score_slot_t* slots = NULL;
    verify_prepared_score_t* scores = NULL;
    index_shard_helper_task_t tasks[INDEX_SHARD_HELPER_MAX_TASKS];
    solver_verification_task_input_t
        inputs[INDEX_SHARD_HELPER_MAX_TASKS];
    index_shard_helper_run_stats_t run_stats;
    index_shard_helper_run_status_t run_status;
    solver_verification_retire_context_t retire_context;
    size_t available = 0U;
    size_t candidate_count;
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
        candidate_first < 0 ||
        candidate_end <= candidate_first ||
        candidate_end > result->nres ||
        candidate_end - candidate_first >
            (int)SOLVER_VERIFICATION_WINDOW_CANDIDATES ||
        verify_datalog_enabled() ||
        !solver_payload_candidate_data_fully_resident(solver) ||
        !index_shard_worker_context_active()) {
        return 0;
    }
    candidate_count = (size_t)(candidate_end - candidate_first);
    available = index_shard_helper_available_workers();
    if (!available) {
        return 0;
    }

    memset(&packet, 0, sizeof(packet));
    memset(&snapshot, 0, sizeof(snapshot));
    memset(tasks, 0, sizeof(tasks));
    memset(inputs, 0, sizeof(inputs));
    memset(&run_stats, 0, sizeof(run_stats));
    memset(&retire_context, 0, sizeof(retire_context));

    if (solver_verification_packet_reserve(
            &packet, candidate_count) !=
        SOLVER_PACKET_RESERVE_OK) {
        goto cleanup;
    }
    runtime = calloc(candidate_count, sizeof(*runtime));
    slots = calloc(candidate_count, sizeof(*slots));
    scores = calloc(candidate_count, sizeof(*scores));
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

    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];

        if (solver_poll_worker_stop(solver)) {
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
        if (solver_ab_candidate_prepare(
                candidate,
                result,
                candidate_index,
                field_xy,
                fieldstars,
                dimquads,
                &snapshot,
                current_parity)) {
            goto cleanup;
        }
        if (candidate->action ==
            SOLVER_AB_CANDIDATE_VERIFY) {
            verify_count++;
        }
    }
    if (verify_count < 2U) {
        goto cleanup;
    }
    available = index_shard_helper_prepare_reserve();
    if (!available) {
        goto cleanup;
    }
    verify_wall_start = monotonic_seconds();
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[packet_index];
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
                &candidate_runtime->logaccept);
    }
    peak_budget = solver_verification_wave_memory_budget();
    slot_cursor = 0U;
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[packet_index];
        solver_verification_score_slot_t* slot;
        size_t context_peak;

        if (candidate->action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        slot = &slots[slot_cursor++];
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
                candidate_runtime->logaccept,
                solver->logratio_stoplooking,
                solver->distance_from_quad_bonus,
                FALSE,
                &slot->prepared)) {
            goto cleanup;
        }
        if (solver_poll_worker_stop(solver)) {
            solver->quit_now = TRUE;
            handled = 1;
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

    if (verify_count) {
        task_count = MIN(
            verify_count,
            MIN(available + 1U,
                (size_t)INDEX_SHARD_HELPER_MAX_TASKS));
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
            inputs[task_index].slot_first = slot_cursor;
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

        run_status = index_shard_helper_run_ordered(
            &solver_verification_helper_ops,
            tasks,
            task_count,
            solver_verification_packet_retire,
            &retire_context,
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
        if (retire_context.next_slot != verify_count) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
    } else {
        task_count = 0U;
        run_status = INDEX_SHARD_HELPER_OK;
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
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];

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
            candidate->code_err);
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
                &runtime[packet_index];
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
                verify_destroy_prepared_score(
                    &scores[slot_cursor - 1U]);
                verify_destroy_prepared_hit(slot->prepared);
                slot->prepared = NULL;
                if (solver_handle_hit(
                        solver,
                        &candidate_runtime->match,
                        NULL,
                        FALSE)) {
                    solver->quit_now = TRUE;
                }
                if (unlikely(solver->quit_now)) {
                    handled = 1;
                    goto cleanup;
                }
                break;
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
    solver_ab_descriptor_planner_t planner;
    solver_ab_descriptor_output_t* outputs;
    size_t output_capacity;
} solver_ab_descriptor_workspace_t;

static pthread_key_t solver_ab_descriptor_workspace_key;
static pthread_once_t solver_ab_descriptor_workspace_once =
    PTHREAD_ONCE_INIT;
static int solver_ab_descriptor_workspace_status = EAGAIN;

static void solver_ab_descriptor_release_pairs(
    solver_ab_descriptor_workspace_t* workspace) {
    solver_ab_descriptor_planner_t* planner;

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

static int solver_ab_descriptor_workspace_reserve_outputs(
    solver_ab_descriptor_workspace_t* workspace,
    size_t output_count) {
    solver_ab_descriptor_output_t* outputs;

    if (!workspace || !output_count ||
        output_count > SIZE_MAX / sizeof(*outputs)) {
        return -1;
    }
    if (output_count <= workspace->output_capacity) {
        return 0;
    }
    outputs = calloc(output_count, sizeof(*outputs));
    if (!outputs) {
        return -1;
    }
    free(workspace->outputs);
    workspace->outputs = outputs;
    workspace->output_capacity = output_count;
    return 0;
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
    solver_ab_descriptor_planner_t executor;
    solver_ab_builder_t builder;
    solver_descriptor_status_t status;
    solver_ab_task_t task;
    int eligible[SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS];
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
    memset(&status, 0, sizeof(status));
    memset(&task, 0, sizeof(task));
    output->descriptor_count = 0U;
    output->trailing_numtries = 0U;
    output->trailing_cxdx = 0U;
    output->trailing_meanx = 0U;
    output->has_final_rel_field_noise2 = FALSE;

    executor.field_geometry = input->field_geometry;
    executor.newpoint = input->newpoint;
    executor.dimquads = input->dimquads;
    executor.phase = input->phase;
    executor.pairs = (solver_ab_pair_t*)input->pairs;
    executor.pair_count = input->pair_count;
    executor.combination_eligible = eligible;
    executor.combination_eligible_capacity =
        SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS;
    executor.snapshot.parity = input->parity;
    executor.snapshot.cx_less_than_dx =
        input->cx_less_than_dx;
    executor.snapshot.meanx_less_than_half =
        input->meanx_less_than_half;
    executor.snapshot.cxdx_margin = input->cxdx_margin;

    builder.planner = &executor;
    builder.status = &status;
    builder.descriptor_output = output;
    task.combination_first = input->combination_first;
    task.combination_end = input->combination_end;
    visit_result = solver_ab_walk_task_ranges(
        &executor,
        &task,
        solver_ab_builder_visit_pair_range,
        &builder);
    if (status.cancelled ||
        index_shard_worker_stop_requested()) {
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }
    if (visit_result < 0 ||
        builder.fatal_error ||
        status.evaluation_failed) {
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

#define SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES \
    (2U * 1024U * 1024U)

#define SOLVER_CODEKD_DELIVERY_RANGE_CAPACITY \
    FITSBIN_MMAP_PREFETCH_RANGE_LIMIT
#define SOLVER_CODEKD_DELIVERY_BUDGET_BYTES \
    (2U * 1024U * 1024U)
#define SOLVER_CANDIDATE_DELIVERY_LIMIT \
    (FITSBIN_PREAD_ASYNC_RANGE_LIMIT / 2U)
#define SOLVER_VERIFY_QUERY_LOOKAHEAD 16U
#define SOLVER_CANDIDATE_STAR_LIMIT \
    (SOLVER_CANDIDATE_DELIVERY_LIMIT * DQMAX)
#if SOLVER_CANDIDATE_STAR_LIMIT > FITSBIN_MMAP_PREFETCH_RANGE_LIMIT
#error "candidate Star delivery exceeds mapped-prefetch range capacity"
#endif

typedef enum solver_codekd_result_state {
    SOLVER_CODEKD_RESULT_UNUSED = 0,
    SOLVER_CODEKD_RESULT_READY,
    SOLVER_CODEKD_RESULT_OWNER_REPLAY,
    SOLVER_CODEKD_RESULT_QUERY_FAILED
} solver_codekd_result_state_t;

typedef struct solver_codekd_result_slot {
    size_t hit_first;
    unsigned int hit_count;
    double search_wall_seconds;
    int search_errno;
    solver_codekd_result_state_t state;
} solver_codekd_result_slot_t;

typedef enum solver_codekd_search_packet_state {
    SOLVER_CODEKD_PACKET_ALLOCATED = 0,
    SOLVER_CODEKD_PACKET_DESCRIPTORS_READY,
    SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE,
    SOLVER_CODEKD_PACKET_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE,
    SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER,
    SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_EXECUTING,
    SOLVER_CODEKD_PACKET_RESULTS_READY,
    SOLVER_CODEKD_PACKET_RETIRING,
    SOLVER_CODEKD_PACKET_RETIRED,
    SOLVER_CODEKD_PACKET_STOPPED,
    SOLVER_CODEKD_PACKET_FAILED
} solver_codekd_search_packet_state_t;

typedef enum solver_codekd_page_plan_reason {
    SOLVER_CODEKD_PAGE_PLAN_NONE = 0,
    SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE,
    SOLVER_CODEKD_PAGE_PLAN_ALLOCATION,
    SOLVER_CODEKD_PAGE_PLAN_SOURCE_MISMATCH,
    SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE,
    SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET,
    SOLVER_CODEKD_PAGE_PLAN_RANGE_CAPACITY,
    SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED,
    SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR,
    SOLVER_CODEKD_PAGE_PLAN_CANCELLED
} solver_codekd_page_plan_reason_t;

typedef struct solver_codekd_page_entry {
    uintptr_t mapping_begin;
    uintptr_t mapping_end;
    uintptr_t page_key;
    uintptr_t populate_begin;
    uintptr_t populate_end;
} solver_codekd_page_entry_t;

typedef struct solver_codekd_page_set {
    solver_codekd_page_entry_t* entries;
    size_t* hash_indices;
    unsigned int* hash_generations;
    size_t count;
    size_t capacity;
    size_t hash_capacity;
    unsigned int generation;
} solver_codekd_page_set_t;

typedef struct solver_codekd_page_workspace {
    solver_codekd_page_set_t descriptor;
    solver_codekd_page_set_t group;
    solver_codekd_page_entry_t* sort_entries;
    fitsbin_prefetch_range_t* sealed_ranges;
    size_t page_size;
    size_t page_limit;
    size_t sealed_range_capacity;
} solver_codekd_page_workspace_t;

typedef struct solver_codekd_page_plan_stats {
    size_t descriptors_total;
    size_t descriptors_planned;
    size_t descriptor_splits;
    size_t raw_ranges;
    size_t unique_pages;
    size_t ranges_before_dedup;
    size_t ranges_after_dedup;
    size_t logical_bytes;
    size_t aligned_bytes;
    size_t overread_bytes;
    size_t refusal_counts[SOLVER_CODEKD_PAGE_PLAN_CANCELLED + 1U];
} solver_codekd_page_plan_stats_t;

typedef struct solver_codekd_verification_snapshot {
    MatchObj match_template;
    const verify_field_t* field;
    int index_cutnside;
    int indexid;
    int healpix;
    int hpnside;
    double index_jitter;
    double verify_pix;
    double distractor_ratio;
    double logratio_bail_threshold;
    double logaccept;
    double logratio_stoplooking;
    anbool distance_from_quad_bonus;
    anbool enabled;
} solver_codekd_verification_snapshot_t;

struct solver_candidate_delivery_record {
    unsigned int quadid;
    unsigned int stars[DQMAX];
    double starxyz[DQMAX * 3];
    int prepared_fieldstars[DQMAX];
    double prepared_fieldxy[DQMAX * 2];
    tan_t prepared_wcs;
    double prepared_scale;
    verify_index_query_t* verify_query;
    size_t descriptor_index;
    solver_ab_candidate_action_t plan_action;
    double verify_center[3];
    double verify_radius;
    double verify_radius2;
    double prepared_verify_pix2;
    double prepared_logaccept;
    double prepared_distractor_ratio;
    double prepared_logratio_bail_threshold;
    double prepared_logratio_stoplooking;
    double prepared_field_maxx;
    double prepared_field_maxy;
    verify_prepared_hit_t* prepared_verification;
    verify_prepared_score_t prepared_score;
    anbool prepared_parity;
    anbool candidate_prepared;
    anbool verify_delivery_fallback;
    anbool verify_query_captured;
    anbool prepared_distance_from_quad_bonus;
    anbool verification_score_ready;
};

typedef struct solver_codekd_search_packet {
    solver_ab_descriptor_output_t* descriptors;
    const kdtree_t* tree;
    const quadfile_t* quads;
    startree_t* starkd;
    solver_codekd_result_slot_t* slots;
    u32* inds;
    double* sdists;
    size_t hit_capacity;
    size_t hit_count;
    size_t first;
    size_t count;
    size_t next_descriptor;
    size_t plan_first;
    size_t plan_end;
    size_t plan_range_count;
    size_t plan_logical_bytes;
    size_t pending_descriptor_raw_ranges;
    size_t pending_descriptor_logical_bytes;
    size_t candidate_count;
    size_t candidate_capacity;
    size_t candidate_cursor;
    size_t candidate_window_first;
    size_t candidate_window_count;
    size_t candidate_window_offset;
    size_t candidate_star_count;
    size_t candidate_quad_ready_count;
    size_t candidate_star_ready_count;
    size_t candidate_verify_query_count;
    size_t verify_plan_first;
    size_t verify_plan_end;
    size_t verify_topology_end;
    size_t verify_plan_range_count;
    size_t verify_plan_logical_bytes;
    size_t verify_query_budget;
    size_t verify_query_bytes;
    size_t verify_prepared_count;
    size_t verify_prepared_retained_bytes;
    size_t verify_prepared_transient_bytes;
    size_t verify_sweep_range_count;
    size_t verify_sweep_aligned_bytes;
    size_t verify_sweep_storage_bytes;
    size_t verify_pending_query_index;
    size_t verify_pending_query_raw_ranges;
    size_t verify_pending_query_logical_bytes;
    size_t verify_delivery_budget;
    size_t retire_descriptor;
    size_t retire_hit_offset;
    unsigned long long sequence;
    unsigned int candidate_starids[SOLVER_CANDIDATE_STAR_LIMIT];
    solver_candidate_delivery_record_t* candidate_records;
    fitsbin_pread_range_t* verify_sweep_reads;
    verify_mapped_page_buffer_t* verify_sweep_buffers;
    unsigned char* verify_sweep_storage;
    fitsbin_t* delivery_source;
    fitsbin_payload_io_ticket_t* delivery_ticket;
    solver_codekd_page_workspace_t* page_workspace;
    solver_codekd_page_plan_stats_t page_stats;
    solver_codekd_page_plan_reason_t page_plan_reason;
    unsigned int candidate_quad_submitted;
    unsigned int candidate_quad_ready;
    unsigned int candidate_quad_fallback;
    unsigned int candidate_star_submitted;
    unsigned int candidate_star_ready;
    unsigned int candidate_star_fallback;
    unsigned long long candidate_delivery_windows;
    unsigned long long candidate_quad_ready_rows;
    unsigned long long candidate_star_ready_rows;
    unsigned long long candidate_retired_rows;
    unsigned long long candidate_native_rows;
    unsigned long long verification_page_queries;
    unsigned long long verification_page_queries_planned;
    unsigned long long verification_page_prefixes;
    unsigned long long verification_page_submitted;
    unsigned long long verification_page_ready;
    unsigned long long verification_page_fallback;
    unsigned long long verification_page_ready_rows;
    unsigned long long verification_page_ranges;
    unsigned long long verification_page_logical_bytes;
    unsigned long long verification_page_aligned_bytes;
    unsigned long long candidate_math_prepared;
    unsigned long long verification_score_batches_prepared;
    unsigned long long verification_score_contexts_prepared;
    unsigned long long verification_score_batches_executed;
    unsigned long long verification_score_contexts_completed;
    unsigned long long verification_score_work_units_completed;
    unsigned long long verification_score_fallback_batches;
    unsigned long long verification_score_stopped_batches;
    double verification_score_wall_seconds;
    int dimquads;
    anbool plan_complete;
    anbool pending_descriptor_plan;
    anbool detailed;
    anbool use_radec;
    anbool star_delivery_eligible;
    anbool candidate_quad_delivery_disabled;
    anbool candidate_star_delivery_disabled;
    anbool verification_delivery_disabled;
    anbool verify_plan_complete;
    anbool verify_sweep_plan_complete;
    anbool verify_pending_query_plan;
    anbool retire_descriptor_started;
    solver_codekd_search_packet_state_t state;
} solver_codekd_search_packet_t;

typedef struct solver_codekd_packet_task_input {
    solver_ab_descriptor_task_input_t descriptor;
    const kdtree_t* tree;
    const quadfile_t* quads;
    startree_t* starkd;
    anbool detailed;
    anbool use_radec;
    double field_minx;
    double field_maxx;
    double field_miny;
    double field_maxy;
    double abscale_low;
    double abscale_high;
    double funits_lower;
    double funits_upper;
    solver_codekd_verification_snapshot_t verification;
} solver_codekd_packet_task_input_t;

typedef struct solver_codekd_packet_wave {
    solver_codekd_packet_task_input_t* inputs;
    solver_codekd_search_packet_t* packets;
    index_shard_staged_task_t* tasks;
    size_t capacity;
} solver_codekd_packet_wave_t;

static solver_codekd_packet_wave_t*
solver_codekd_packet_wave_create(size_t capacity) {
    solver_codekd_packet_wave_t* wave;

    if (!capacity ||
        capacity > SIZE_MAX / sizeof(*wave->inputs) ||
        capacity > SIZE_MAX / sizeof(*wave->packets) ||
        capacity > SIZE_MAX / sizeof(*wave->tasks)) {
        return NULL;
    }
    wave = calloc(1, sizeof(*wave));
    if (!wave) {
        return NULL;
    }
    wave->inputs = calloc(capacity, sizeof(*wave->inputs));
    wave->packets = calloc(capacity, sizeof(*wave->packets));
    wave->tasks = calloc(capacity, sizeof(*wave->tasks));
    if (!wave->inputs || !wave->packets || !wave->tasks) {
        free(wave->tasks);
        free(wave->packets);
        free(wave->inputs);
        free(wave);
        return NULL;
    }
    wave->capacity = capacity;
    return wave;
}

static void solver_codekd_packet_wave_destroy(
    solver_codekd_packet_wave_t* wave) {
    if (!wave) {
        return;
    }
    free(wave->tasks);
    free(wave->packets);
    free(wave->inputs);
    free(wave);
}

static int solver_codekd_packet_finish_codekd(
    solver_codekd_search_packet_t* packet);

typedef struct solver_codekd_page_plan {
    fitsbin_t* source;
    solver_codekd_page_workspace_t* workspace;
    fitsbin_payload_io_cancel_check_fn cancelled;
    void* cancel_opaque;
    size_t raw_ranges;
    size_t logical_bytes;
    solver_codekd_page_plan_reason_t reason;
    anbool cancellation_observed;
    anbool enabled;
} solver_codekd_page_plan_t;

static uint64_t solver_codekd_page_hash(
    uintptr_t mapping_begin,
    uintptr_t page_key) {
    uint64_t value = (uint64_t)mapping_begin;

    value ^= (uint64_t)page_key + UINT64_C(0x9e3779b97f4a7c15) +
        (value << 6) + (value >> 2);
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return value;
}

static void solver_codekd_page_set_cleanup(
    solver_codekd_page_set_t* set) {
    if (!set) {
        return;
    }
    free(set->entries);
    free(set->hash_indices);
    free(set->hash_generations);
    memset(set, 0, sizeof(*set));
}

static int solver_codekd_page_set_init(
    solver_codekd_page_set_t* set,
    size_t capacity,
    size_t hash_capacity) {
    if (!set || !capacity || hash_capacity < capacity ||
        (hash_capacity & (hash_capacity - 1U)) != 0U ||
        capacity > SIZE_MAX / sizeof(*set->entries) ||
        hash_capacity > SIZE_MAX / sizeof(*set->hash_indices) ||
        hash_capacity > SIZE_MAX / sizeof(*set->hash_generations)) {
        return -1;
    }
    memset(set, 0, sizeof(*set));
    set->entries = malloc(capacity * sizeof(*set->entries));
    set->hash_indices = malloc(
        hash_capacity * sizeof(*set->hash_indices));
    set->hash_generations = calloc(
        hash_capacity, sizeof(*set->hash_generations));
    if (!set->entries || !set->hash_indices ||
        !set->hash_generations) {
        solver_codekd_page_set_cleanup(set);
        return 1;
    }
    set->capacity = capacity;
    set->hash_capacity = hash_capacity;
    return 0;
}

static void solver_codekd_page_set_reset(
    solver_codekd_page_set_t* set) {
    if (!set) {
        return;
    }
    set->count = 0U;
    set->generation++;
    if (!set->generation) {
        memset(set->hash_generations, 0,
               set->hash_capacity *
                   sizeof(*set->hash_generations));
        set->generation = 1U;
    }
}

static int solver_codekd_page_set_find(
    const solver_codekd_page_set_t* set,
    uintptr_t mapping_begin,
    uintptr_t page_key,
    size_t* slot_out) {
    size_t slot;
    size_t probes;

    if (!set || !set->hash_capacity || !set->generation) {
        return 0;
    }
    slot = (size_t)solver_codekd_page_hash(
        mapping_begin, page_key) & (set->hash_capacity - 1U);
    for (probes = 0U; probes < set->hash_capacity; probes++) {
        size_t entry_index;

        if (set->hash_generations[slot] != set->generation) {
            if (slot_out) {
                *slot_out = slot;
            }
            return 0;
        }
        entry_index = set->hash_indices[slot];
        if (entry_index < set->count &&
            set->entries[entry_index].mapping_begin == mapping_begin &&
            set->entries[entry_index].page_key == page_key) {
            if (slot_out) {
                *slot_out = slot;
            }
            return 1;
        }
        slot = (slot + 1U) & (set->hash_capacity - 1U);
    }
    return -1;
}

static int solver_codekd_page_set_add(
    solver_codekd_page_set_t* set,
    const solver_codekd_page_entry_t* entry) {
    size_t slot = 0U;
    int found;

    if (!set || !entry) {
        return -1;
    }
    found = solver_codekd_page_set_find(
        set, entry->mapping_begin, entry->page_key, &slot);
    if (found) {
        return found > 0 ? 0 : -1;
    }
    if (set->count >= set->capacity) {
        return 1;
    }
    set->entries[set->count] = *entry;
    set->hash_indices[slot] = set->count;
    set->hash_generations[slot] = set->generation;
    set->count++;
    return 0;
}

static void solver_codekd_page_workspace_cleanup(
    solver_codekd_page_workspace_t* workspace) {
    if (!workspace) {
        return;
    }
    solver_codekd_page_set_cleanup(&workspace->descriptor);
    solver_codekd_page_set_cleanup(&workspace->group);
    free(workspace->sort_entries);
    free(workspace->sealed_ranges);
    free(workspace);
}

static int solver_codekd_page_workspace_create(
    solver_codekd_page_workspace_t** workspace_out) {
    solver_codekd_page_workspace_t* workspace;
    size_t hash_capacity = 1U;
    size_t page_limit;
    long detected_page_size;
    int status;

    if (!workspace_out) {
        return -1;
    }
    *workspace_out = NULL;
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0) {
        return 1;
    }
    page_limit = SOLVER_CODEKD_DELIVERY_BUDGET_BYTES /
        (size_t)detected_page_size;
    if (!page_limit || page_limit > SIZE_MAX / 2U) {
        return 1;
    }
    while (hash_capacity < page_limit * 2U) {
        if (hash_capacity > SIZE_MAX / 2U) {
            return 1;
        }
        hash_capacity *= 2U;
    }
    workspace = calloc(1, sizeof(*workspace));
    if (!workspace) {
        return 1;
    }
    status = solver_codekd_page_set_init(
        &workspace->descriptor, page_limit, hash_capacity);
    if (!status) {
        status = solver_codekd_page_set_init(
            &workspace->group, page_limit, hash_capacity);
    }
    if (!status) {
        workspace->sort_entries = malloc(
            page_limit * sizeof(*workspace->sort_entries));
        workspace->sealed_ranges = malloc(
            SOLVER_CODEKD_DELIVERY_RANGE_CAPACITY *
                sizeof(*workspace->sealed_ranges));
        if (!workspace->sort_entries || !workspace->sealed_ranges) {
            status = 1;
        }
    }
    if (status) {
        solver_codekd_page_workspace_cleanup(workspace);
        return status;
    }
    workspace->page_size = (size_t)detected_page_size;
    workspace->page_limit = page_limit;
    workspace->sealed_range_capacity =
        SOLVER_CODEKD_DELIVERY_RANGE_CAPACITY;
    solver_codekd_page_set_reset(&workspace->descriptor);
    solver_codekd_page_set_reset(&workspace->group);
    *workspace_out = workspace;
    return 0;
}

static int solver_codekd_page_plan_enabled(
    void* opaque,
    void* mapping) {
    solver_codekd_page_plan_t* plan = opaque;
    fitsbin_t* source = mapping;

    if (!plan || !plan->enabled || !source ||
        plan->source != source) {
        return FALSE;
    }
    if (plan->cancelled &&
        plan->cancelled(plan->cancel_opaque)) {
        plan->cancellation_observed = TRUE;
        plan->reason = SOLVER_CODEKD_PAGE_PLAN_CANCELLED;
        return FALSE;
    }
    return TRUE;
}

static void solver_codekd_page_plan_add_size(
    size_t* total,
    size_t increment) {
    if (!total) {
        return;
    }
    if (increment > SIZE_MAX - *total) {
        *total = SIZE_MAX;
    } else {
        *total += increment;
    }
}

static int solver_codekd_page_plan_emit(
    void* opaque,
    const kdtree_prefetch_hint_t* hint) {
    solver_codekd_page_plan_t* plan = opaque;
    solver_codekd_page_workspace_t* workspace;
    const void* map_base;
    const void* range_start;
    size_t map_size;
    size_t range_size;
    uintptr_t map_begin;
    uintptr_t map_end;
    uintptr_t request_begin;
    uintptr_t request_end;
    uintptr_t page_key;
    int resolved;

    if (!plan || !hint || !hint->mapping ||
        !hint->address || !hint->length) {
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    if (plan->cancelled &&
        plan->cancelled(plan->cancel_opaque)) {
        plan->cancellation_observed = TRUE;
        plan->reason = SOLVER_CODEKD_PAGE_PLAN_CANCELLED;
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    if (hint->priority != KDTREE_PREFETCH_PRIORITY_LEAF ||
        (hint->kind != KDTREE_PREFETCH_ARRAY_DATA &&
         hint->kind != KDTREE_PREFETCH_ARRAY_PERM)) {
        return KDTREE_PREFETCH_EMIT_CONTINUE;
    }
    if (plan->source != hint->mapping || !plan->workspace) {
        plan->reason = SOLVER_CODEKD_PAGE_PLAN_SOURCE_MISMATCH;
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    workspace = plan->workspace;
    if (plan->raw_ranges < SIZE_MAX) {
        plan->raw_ranges++;
    }
    solver_codekd_page_plan_add_size(
        &plan->logical_bytes, hint->length);
    resolved = fitsbin_resolve_mapped_range(
        plan->source,
        hint->address,
        hint->length,
        &map_base,
        &map_size,
        &range_start,
        &range_size);
    if (resolved <= 0 || range_size != hint->length) {
        plan->reason = resolved < 0
            ? SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE
            : SOLVER_CODEKD_PAGE_PLAN_SOURCE_MISMATCH;
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    map_begin = (uintptr_t)map_base;
    request_begin = (uintptr_t)range_start;
    if (map_size > UINTPTR_MAX - map_begin ||
        range_size > UINTPTR_MAX - request_begin) {
        plan->reason = SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    map_end = map_begin + map_size;
    request_end = request_begin + range_size;
    if (request_begin < map_begin || request_end > map_end) {
        plan->reason = SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    page_key = request_begin -
        request_begin % (uintptr_t)workspace->page_size;
    while (page_key < request_end) {
        solver_codekd_page_entry_t entry;
        uintptr_t next_page;
        int add_status;

        if ((uintptr_t)workspace->page_size > UINTPTR_MAX - page_key) {
            plan->reason = SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
            return KDTREE_PREFETCH_EMIT_ERROR;
        }
        next_page = page_key + (uintptr_t)workspace->page_size;
        entry.mapping_begin = map_begin;
        entry.mapping_end = map_end;
        entry.page_key = page_key;
        entry.populate_begin = MAX(page_key, map_begin);
        entry.populate_end = MIN(next_page, map_end);
        if (entry.populate_end <= entry.populate_begin) {
            plan->reason = SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
            return KDTREE_PREFETCH_EMIT_ERROR;
        }
        add_status = solver_codekd_page_set_add(
            &workspace->descriptor, &entry);
        if (add_status > 0) {
            plan->reason = SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET;
            return KDTREE_PREFETCH_EMIT_REFUSED;
        }
        if (add_status < 0) {
            plan->reason = SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
            return KDTREE_PREFETCH_EMIT_ERROR;
        }
        page_key = next_page;
    }
    return KDTREE_PREFETCH_EMIT_CONTINUE;
}

static int solver_codekd_page_entry_compare(
    const void* left_opaque,
    const void* right_opaque) {
    const solver_codekd_page_entry_t* left = left_opaque;
    const solver_codekd_page_entry_t* right = right_opaque;

    if (left->mapping_begin < right->mapping_begin) {
        return -1;
    }
    if (left->mapping_begin > right->mapping_begin) {
        return 1;
    }
    if (left->mapping_end < right->mapping_end) {
        return -1;
    }
    if (left->mapping_end > right->mapping_end) {
        return 1;
    }
    if (left->populate_begin < right->populate_begin) {
        return -1;
    }
    if (left->populate_begin > right->populate_begin) {
        return 1;
    }
    return 0;
}

/*
 * Materialize the physical union only after page-key deduplication. Sorting
 * this private page list cannot change descriptor or KD-hit order.
 */
static int solver_codekd_page_plan_seal_union(
    solver_codekd_page_workspace_t* workspace,
    anbool include_descriptor,
    size_t* unique_pages_out,
    size_t* range_count_out,
    size_t* aligned_bytes_out) {
    size_t unique_pages = 0U;
    size_t range_count = 0U;
    size_t aligned_bytes = 0U;
    size_t i;

    if (!workspace || !unique_pages_out || !range_count_out ||
        !aligned_bytes_out) {
        return -1;
    }
    if (workspace->group.count > workspace->page_limit) {
        return -1;
    }
    for (i = 0U; i < workspace->group.count; i++) {
        workspace->sort_entries[unique_pages++] =
            workspace->group.entries[i];
    }
    if (include_descriptor) {
        for (i = 0U; i < workspace->descriptor.count; i++) {
            const solver_codekd_page_entry_t* entry =
                &workspace->descriptor.entries[i];
            int found = solver_codekd_page_set_find(
                &workspace->group,
                entry->mapping_begin,
                entry->page_key,
                NULL);

            if (found < 0) {
                return -1;
            }
            if (found) {
                continue;
            }
            if (unique_pages >= workspace->page_limit) {
                return 1;
            }
            workspace->sort_entries[unique_pages++] = *entry;
        }
    }
    if (!unique_pages) {
        *unique_pages_out = 0U;
        *range_count_out = 0U;
        *aligned_bytes_out = 0U;
        return 0;
    }
    qsort(workspace->sort_entries,
          unique_pages,
          sizeof(*workspace->sort_entries),
          solver_codekd_page_entry_compare);
    for (i = 0U; i < unique_pages; i++) {
        const solver_codekd_page_entry_t* entry =
            &workspace->sort_entries[i];
        fitsbin_prefetch_range_t* previous = range_count
            ? &workspace->sealed_ranges[range_count - 1U]
            : NULL;
        uintptr_t previous_begin = previous
            ? (uintptr_t)previous->data
            : 0U;
        uintptr_t previous_end = previous_begin;

        if (previous && previous->size <=
                UINTPTR_MAX - previous_begin) {
            previous_end += previous->size;
        }
        if (previous &&
            entry->mapping_begin ==
                workspace->sort_entries[i - 1U].mapping_begin &&
            entry->mapping_end ==
                workspace->sort_entries[i - 1U].mapping_end &&
            entry->populate_begin <= previous_end) {
            if (entry->populate_end > previous_end) {
                previous->size = (size_t)(
                    entry->populate_end - previous_begin);
            }
            continue;
        }
        if (range_count >= workspace->sealed_range_capacity) {
            return 2;
        }
        workspace->sealed_ranges[range_count].data =
            (const void*)entry->populate_begin;
        workspace->sealed_ranges[range_count].size =
            (size_t)(entry->populate_end - entry->populate_begin);
        range_count++;
    }
    for (i = 0U; i < range_count; i++) {
        if (workspace->sealed_ranges[i].size >
            SIZE_MAX - aligned_bytes) {
            return -1;
        }
        aligned_bytes += workspace->sealed_ranges[i].size;
    }
    if (aligned_bytes > SOLVER_CODEKD_DELIVERY_BUDGET_BYTES) {
        return 1;
    }
    *unique_pages_out = unique_pages;
    *range_count_out = range_count;
    *aligned_bytes_out = aligned_bytes;
    return 0;
}

static int solver_codekd_page_set_merge_descriptor(
    solver_codekd_page_workspace_t* workspace) {
    size_t i;

    if (!workspace) {
        return -1;
    }
    for (i = 0U; i < workspace->descriptor.count; i++) {
        int add_status = solver_codekd_page_set_add(
            &workspace->group,
            &workspace->descriptor.entries[i]);

        if (add_status) {
            return -1;
        }
    }
    return 0;
}

/*
 * Check whether the current descriptor page set fits in the group without
 * sorting or materializing physical ranges. The group hash already defines
 * the exact mapped-page identity used by the final seal.
 */
static int solver_codekd_page_set_union_fits(
    const solver_codekd_page_workspace_t* workspace) {
    size_t additional = 0U;
    size_t i;

    if (!workspace ||
        workspace->group.count > workspace->page_limit) {
        return -1;
    }
    for (i = 0U; i < workspace->descriptor.count; i++) {
        const solver_codekd_page_entry_t* entry =
            &workspace->descriptor.entries[i];
        int found = solver_codekd_page_set_find(
            &workspace->group,
            entry->mapping_begin,
            entry->page_key,
            NULL);

        if (found < 0) {
            return -1;
        }
        if (found) {
            continue;
        }
        if (additional >=
            workspace->page_limit - workspace->group.count) {
            return 0;
        }
        additional++;
    }
    return 1;
}

static void solver_codekd_page_plan_record_refusal(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason) {
    if (!packet || reason <= SOLVER_CODEKD_PAGE_PLAN_NONE ||
        reason > SOLVER_CODEKD_PAGE_PLAN_CANCELLED) {
        return;
    }
    packet->page_plan_reason = reason;
    packet->page_stats.refusal_counts[reason]++;
}

static anbool solver_codekd_packet_plan_cancelled(void* opaque) {
    (void)opaque;
    return index_shard_worker_stop_requested();
}

/*
 * Plan the next contiguous descriptor slice transactionally. A positive
 * result exposes one complete sealed plan. Zero means every descriptor now
 * has either a populated-query result pending or exact owner replay. No
 * incomplete prefix is ever submitted.
 */
static int solver_codekd_search_packet_prepare_next_plan(
    solver_codekd_search_packet_t* packet,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque) {
    solver_codekd_page_workspace_t* workspace;
    fitsbin_t* source;
    size_t group_first = 0U;
    size_t group_end = 0U;
    size_t group_logical_bytes = 0U;

    if (!packet || !packet->tree ||
        !packet->tree->io || !packet->tree->io_is_fitsbin ||
        !packet->descriptors || !packet->slots ||
        !packet->page_workspace || !cancelled ||
        packet->state != SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
        packet->next_descriptor > packet->count) {
        return -1;
    }
    workspace = packet->page_workspace;
    source = (fitsbin_t*)packet->tree->io;
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    solver_codekd_page_set_reset(&workspace->group);

    while (packet->next_descriptor < packet->count) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[
                packet->next_descriptor];
        solver_codekd_page_plan_t plan;
        kdtree_prefetch_sink_t sink;
        solver_codekd_page_plan_reason_t reason =
            SOLVER_CODEKD_PAGE_PLAN_NONE;
        size_t descriptor_raw_ranges;
        size_t descriptor_logical_bytes;
        size_t unique_pages;
        size_t range_count;
        size_t aligned_bytes;
        int prepare_status;
        int seal_status;

        if (cancelled(cancel_opaque)) {
            solver_codekd_page_plan_record_refusal(
                packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return 2;
        }
        if (packet->pending_descriptor_plan) {
            descriptor_raw_ranges =
                packet->pending_descriptor_raw_ranges;
            descriptor_logical_bytes =
                packet->pending_descriptor_logical_bytes;
        } else {
            solver_codekd_page_set_reset(&workspace->descriptor);
            memset(&plan, 0, sizeof(plan));
            memset(&sink, 0, sizeof(sink));
            plan.source = source;
            plan.workspace = workspace;
            plan.cancelled = cancelled;
            plan.cancel_opaque = cancel_opaque;
            plan.enabled = TRUE;
            sink.userdata = &plan;
            sink.enabled = solver_codekd_page_plan_enabled;
            sink.emit = solver_codekd_page_plan_emit;
            prepare_status = kdtree_rangesearch_prefetch_prepare(
                packet->tree,
                descriptor->code,
                descriptor->tol2,
                SOLVER_CODEKD_SEARCH_OPTIONS,
                &sink);
            if (plan.cancellation_observed) {
                solver_codekd_page_plan_record_refusal(
                    packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
                packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                return 2;
            }
            descriptor_raw_ranges = plan.raw_ranges;
            descriptor_logical_bytes = plan.logical_bytes;
            if (prepare_status !=
                    KDTREE_PREFETCH_PREPARE_COMPLETE ||
                !workspace->descriptor.count) {
                if (prepare_status ==
                    KDTREE_PREFETCH_PREPARE_NOT_APPLICABLE) {
                    reason = SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE;
                } else if (prepare_status ==
                           KDTREE_PREFETCH_PREPARE_REFUSED) {
                    reason = plan.reason !=
                            SOLVER_CODEKD_PAGE_PLAN_NONE
                        ? plan.reason
                        : SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET;
                } else if (prepare_status !=
                           KDTREE_PREFETCH_PREPARE_COMPLETE) {
                    reason = plan.reason !=
                            SOLVER_CODEKD_PAGE_PLAN_NONE
                        ? plan.reason
                        : SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE;
                } else {
                    reason = SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE;
                }
                solver_codekd_page_plan_record_refusal(packet, reason);
                packet->slots[packet->next_descriptor].state =
                    SOLVER_CODEKD_RESULT_OWNER_REPLAY;
                packet->next_descriptor++;
                if (workspace->group.count) {
                    break;
                }
                continue;
            }
        }

        seal_status = solver_codekd_page_plan_seal_union(
            workspace,
            TRUE,
            &unique_pages,
            &range_count,
            &aligned_bytes);
        if (seal_status) {
            reason = seal_status == 2
                ? SOLVER_CODEKD_PAGE_PLAN_RANGE_CAPACITY
                : SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET;
            if (workspace->group.count) {
                packet->pending_descriptor_plan = TRUE;
                packet->pending_descriptor_raw_ranges =
                    descriptor_raw_ranges;
                packet->pending_descriptor_logical_bytes =
                    descriptor_logical_bytes;
                packet->page_stats.descriptor_splits++;
                break;
            }
            solver_codekd_page_plan_record_refusal(packet, reason);
            packet->slots[packet->next_descriptor].state =
                SOLVER_CODEKD_RESULT_OWNER_REPLAY;
            packet->pending_descriptor_plan = FALSE;
            packet->pending_descriptor_raw_ranges = 0U;
            packet->pending_descriptor_logical_bytes = 0U;
            packet->next_descriptor++;
            continue;
        }
        if (!workspace->group.count) {
            group_first = packet->next_descriptor;
        }
        if (solver_codekd_page_set_merge_descriptor(workspace)) {
            return -1;
        }
        solver_codekd_page_plan_add_size(
            &packet->page_stats.raw_ranges,
            descriptor_raw_ranges);
        solver_codekd_page_plan_add_size(
            &packet->page_stats.logical_bytes,
            descriptor_logical_bytes);
        solver_codekd_page_plan_add_size(
            &group_logical_bytes,
            descriptor_logical_bytes);
        packet->pending_descriptor_plan = FALSE;
        packet->pending_descriptor_raw_ranges = 0U;
        packet->pending_descriptor_logical_bytes = 0U;
        packet->next_descriptor++;
        group_end = packet->next_descriptor;
        if (workspace->group.count >= workspace->page_limit) {
            break;
        }
    }

    if (!workspace->group.count) {
        if (packet->next_descriptor != packet->count) {
            return -1;
        }
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    {
        size_t unique_pages;
        size_t range_count;
        size_t aligned_bytes;
        int seal_status = solver_codekd_page_plan_seal_union(
            workspace,
            FALSE,
            &unique_pages,
            &range_count,
            &aligned_bytes);

        if (seal_status || !unique_pages || !range_count) {
            return -1;
        }
        packet->plan_first = group_first;
        packet->plan_end = group_end;
        packet->plan_range_count = range_count;
        packet->plan_logical_bytes = group_logical_bytes;
        packet->plan_complete = TRUE;
        packet->page_stats.descriptors_planned +=
            packet->plan_end - packet->plan_first;
        packet->page_stats.unique_pages += unique_pages;
        packet->page_stats.ranges_before_dedup += unique_pages;
        packet->page_stats.ranges_after_dedup += range_count;
        solver_codekd_page_plan_add_size(
            &packet->page_stats.aligned_bytes, aligned_bytes);
        /*
         * Raw hints can overlap. This positive delta is diagnostic only and
         * is not an exact physical-overread or coverage measurement.
         */
        if (aligned_bytes > group_logical_bytes) {
            solver_codekd_page_plan_add_size(
                &packet->page_stats.overread_bytes,
                aligned_bytes - group_logical_bytes);
        }
        packet->state = SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE;
    }
    return 1;
}

static void solver_codekd_packet_reset_verify_plan(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;
    size_t retained_bytes = 0U;

    if (!packet) {
        return;
    }
    packet->verify_plan_first = 0U;
    packet->verify_plan_end = 0U;
    packet->verify_topology_end = 0U;
    packet->verify_plan_range_count = 0U;
    packet->verify_plan_logical_bytes = 0U;
    packet->verify_plan_complete = FALSE;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_plan_complete = FALSE;
    for (candidate_index = 0U;
         packet->candidate_records &&
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        size_t query_bytes = verify_index_query_bytes(
            packet->candidate_records[
                candidate_index].verify_query);

        if (query_bytes == SIZE_MAX ||
            query_bytes > SIZE_MAX - retained_bytes) {
            retained_bytes = SIZE_MAX;
            break;
        }
        retained_bytes += query_bytes;
    }
    packet->verify_query_bytes = retained_bytes;
}

static void solver_codekd_packet_disable_verification_delivery(
    solver_codekd_search_packet_t* packet) {
    if (!packet) {
        return;
    }
    solver_codekd_packet_reset_verify_plan(packet);
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    packet->verify_pending_query_plan = FALSE;
    packet->verification_delivery_disabled = TRUE;
    packet->verification_page_fallback++;
    packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
}

/*
 * Build the maximal complete canonical verification-query prefix. Each query
 * is transactional: no page from an incomplete query enters a READY prefix.
 */
static int solver_codekd_search_packet_prepare_verify_plan(
    solver_codekd_search_packet_t* packet,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque) {
    solver_codekd_page_workspace_t* workspace;
    fitsbin_t* source;
    size_t segment_first;
    size_t segment_end;
    size_t group_logical_bytes = 0U;
    size_t candidate_index;

    if (!packet || !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->candidate_records || !packet->page_workspace ||
        !cancelled ||
        packet->state != SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED ||
        packet->candidate_window_offset >=
            packet->candidate_window_count ||
        packet->candidate_star_ready_count !=
            packet->candidate_window_count ||
        packet->verification_delivery_disabled) {
        return -1;
    }
    workspace = packet->page_workspace;
    source = (fitsbin_t*)packet->starkd->tree->io;
    segment_first = packet->candidate_window_offset;
    segment_end = segment_first;
    solver_codekd_packet_reset_verify_plan(packet);
    solver_codekd_page_set_reset(&workspace->group);

    for (candidate_index = segment_first;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        solver_codekd_page_plan_t plan;
        kdtree_prefetch_sink_t sink;
        size_t query_raw_ranges;
        size_t query_logical_bytes;
        size_t unique_pages;
        size_t range_count;
        size_t aligned_bytes;
        int prepare_status;
        int seal_status;

        if (cancelled(cancel_opaque)) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return 2;
        }
        if (record->verify_delivery_fallback) {
            if (workspace->group.count) {
                break;
            }
            solver_codekd_packet_disable_verification_delivery(packet);
            return 0;
        }
        if (record->plan_action !=
            SOLVER_AB_CANDIDATE_VERIFY) {
            segment_end = candidate_index + 1U;
            continue;
        }

        if (packet->verify_pending_query_plan) {
            if (packet->verify_pending_query_index != candidate_index ||
                !workspace->descriptor.count) {
                return -1;
            }
            query_raw_ranges =
                packet->verify_pending_query_raw_ranges;
            query_logical_bytes =
                packet->verify_pending_query_logical_bytes;
        } else {
            solver_codekd_page_set_reset(&workspace->descriptor);
            memset(&plan, 0, sizeof(plan));
            memset(&sink, 0, sizeof(sink));
            plan.source = source;
            plan.workspace = workspace;
            plan.cancelled = cancelled;
            plan.cancel_opaque = cancel_opaque;
            plan.enabled = TRUE;
            sink.userdata = &plan;
            sink.enabled = solver_codekd_page_plan_enabled;
            sink.emit = solver_codekd_page_plan_emit;
            prepare_status = kdtree_rangesearch_prefetch_prepare(
                packet->starkd->tree,
                record->verify_center,
                record->verify_radius2,
                SOLVER_STARKD_VERIFY_SEARCH_OPTIONS,
                &sink);
            if (plan.cancellation_observed) {
                packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                return 2;
            }
            query_raw_ranges = plan.raw_ranges;
            query_logical_bytes = plan.logical_bytes;
            if (prepare_status !=
                KDTREE_PREFETCH_PREPARE_COMPLETE) {
                record->verify_delivery_fallback = TRUE;
                if (workspace->group.count) {
                    break;
                }
                solver_codekd_packet_disable_verification_delivery(
                    packet);
                return 0;
            }
            if (!workspace->descriptor.count) {
                packet->verification_page_queries_planned++;
                segment_end = candidate_index + 1U;
                continue;
            }
        }

        seal_status = solver_codekd_page_plan_seal_union(
            workspace,
            TRUE,
            &unique_pages,
            &range_count,
            &aligned_bytes);
        if (seal_status) {
            if (workspace->group.count) {
                packet->verify_pending_query_plan = TRUE;
                packet->verify_pending_query_index = candidate_index;
                packet->verify_pending_query_raw_ranges =
                    query_raw_ranges;
                packet->verify_pending_query_logical_bytes =
                    query_logical_bytes;
                break;
            }
            record->verify_delivery_fallback = TRUE;
            solver_codekd_packet_disable_verification_delivery(packet);
            return 0;
        }
        if (solver_codekd_page_set_merge_descriptor(workspace)) {
            return -1;
        }
        packet->verification_page_queries_planned++;
        solver_codekd_page_plan_add_size(
            &group_logical_bytes, query_logical_bytes);
        packet->verify_pending_query_plan = FALSE;
        packet->verify_pending_query_index = 0U;
        packet->verify_pending_query_raw_ranges = 0U;
        packet->verify_pending_query_logical_bytes = 0U;
        segment_end = candidate_index + 1U;
        if (workspace->group.count >= workspace->page_limit) {
            break;
        }
    }

    if (segment_end <= segment_first) {
        return -1;
    }
    packet->verify_plan_first = segment_first;
    packet->verify_plan_end = segment_end;
    packet->verify_topology_end = segment_end;
    packet->verify_plan_logical_bytes = group_logical_bytes;
    packet->verify_plan_complete = TRUE;
    packet->verification_page_prefixes++;
    if (!workspace->group.count) {
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE;
        return 0;
    }

    {
        size_t unique_pages;
        size_t range_count;
        size_t aligned_bytes;
        int seal_status = solver_codekd_page_plan_seal_union(
            workspace,
            FALSE,
            &unique_pages,
            &range_count,
            &aligned_bytes);

        if (seal_status || !unique_pages || !range_count) {
            return -1;
        }
        packet->verify_plan_range_count = range_count;
        packet->verification_page_ranges += range_count;
        packet->verification_page_logical_bytes =
            solver_ab_saturating_add(
                packet->verification_page_logical_bytes,
                (unsigned long long)group_logical_bytes);
        packet->verification_page_aligned_bytes =
            solver_ab_saturating_add(
                packet->verification_page_aligned_bytes,
                (unsigned long long)aligned_bytes);
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE;
    }
    return 1;
}

static int solver_codekd_packet_plan_verification_pages(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count) {
    solver_codekd_search_packet_t* packet = opaque;
    int plan_status;

    if (!packet || !cancelled || !ranges || !range_count) {
        errno = EINVAL;
        return -1;
    }
    *range_count = 0U;
    plan_status = solver_codekd_search_packet_prepare_verify_plan(
        packet, cancelled, cancel_opaque);
    if (plan_status < 0) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        errno = EIO;
        return -1;
    }
    if (plan_status == 2) {
        errno = ECANCELED;
        return -1;
    }
    if (!packet->verify_plan_complete ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE ||
        !packet->page_workspace ||
        packet->verify_plan_range_count > range_capacity) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        errno = packet->verify_plan_range_count > range_capacity
            ? E2BIG
            : EINVAL;
        return -1;
    }
    if (!plan_status) {
        return 0;
    }
    if (!packet->verify_plan_range_count) {
        return 1;
    }
    memcpy(
        ranges,
        packet->page_workspace->sealed_ranges,
        packet->verify_plan_range_count * sizeof(*ranges));
    *range_count = packet->verify_plan_range_count;
    return 1;
}

static int solver_codekd_search_packet_release_ticket(
    solver_codekd_search_packet_t* packet) {
    if (!packet) {
        return -1;
    }
    if (!packet->delivery_ticket) {
        return 0;
    }
    if (!packet->delivery_source) {
        return -1;
    }
    /*
     * A terminal READY ticket returns a positive work count. Checked destroy,
     * not the sign of the collection result, is the ownership authority.
     */
    (void)fitsbin_payload_io_ticket_cancel_and_wait(
        packet->delivery_source,
        packet->delivery_ticket);
    if (fitsbin_payload_io_ticket_destroy_checked(
            packet->delivery_ticket)) {
        return -1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    return 0;
}

static void solver_codekd_packet_clear_sweep_storage(
    solver_codekd_search_packet_t* packet) {
    if (!packet) {
        return;
    }
    assert(!packet->delivery_ticket);
    free(packet->verify_sweep_reads);
    free(packet->verify_sweep_buffers);
    free(packet->verify_sweep_storage);
    packet->verify_sweep_reads = NULL;
    packet->verify_sweep_buffers = NULL;
    packet->verify_sweep_storage = NULL;
    packet->verify_sweep_storage_bytes = 0U;
}

static void solver_codekd_record_clear_prepared_verification(
    solver_candidate_delivery_record_t* record) {
    if (!record) {
        return;
    }
    verify_destroy_prepared_score(&record->prepared_score);
    verify_destroy_prepared_hit(record->prepared_verification);
    record->prepared_verification = NULL;
    record->prepared_verify_pix2 = 0.0;
    record->prepared_logaccept = 0.0;
    record->prepared_distractor_ratio = 0.0;
    record->prepared_logratio_bail_threshold = 0.0;
    record->prepared_logratio_stoplooking = 0.0;
    record->prepared_field_maxx = 0.0;
    record->prepared_field_maxy = 0.0;
    record->prepared_distance_from_quad_bonus = FALSE;
    record->verification_score_ready = FALSE;
}

static void solver_codekd_record_clear_verification_speculation(
    solver_candidate_delivery_record_t* record) {
    if (!record) {
        return;
    }
    verify_destroy_index_query(record->verify_query);
    record->verify_query = NULL;
    record->verify_query_captured = FALSE;
    solver_codekd_record_clear_prepared_verification(record);
}

static int solver_codekd_search_packet_cleanup(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;

    if (!packet ||
        solver_codekd_search_packet_release_ticket(packet)) {
        return -1;
    }
    solver_codekd_packet_clear_sweep_storage(packet);
    for (candidate_index = 0U;
         packet->candidate_records &&
         candidate_index < packet->candidate_capacity;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &packet->candidate_records[candidate_index]);
        verify_destroy_index_query(
            packet->candidate_records[
                candidate_index].verify_query);
        packet->candidate_records[
            candidate_index].verify_query = NULL;
        packet->candidate_records[
            candidate_index].verify_query_captured = FALSE;
    }
    solver_codekd_page_workspace_cleanup(packet->page_workspace);
    free(packet->slots);
    free(packet->inds);
    free(packet->sdists);
    free(packet->candidate_records);
    memset(packet, 0, sizeof(*packet));
    return 0;
}

static void solver_codekd_packet_profile_accumulate(
    solver_t* solver,
    const solver_codekd_search_packet_t* packet) {
    const solver_codekd_page_plan_stats_t* stats;

    if (!solver || !packet) {
        return;
    }
    stats = &packet->page_stats;
    solver->profile.page_plan_descriptors_total +=
        stats->descriptors_total;
    solver->profile.page_plan_descriptors_complete +=
        stats->descriptors_planned;
    solver->profile.page_plan_descriptor_splits +=
        stats->descriptor_splits;
    solver->profile.page_plan_raw_ranges += stats->raw_ranges;
    solver->profile.page_plan_unique_pages += stats->unique_pages;
    solver->profile.page_plan_ranges_after_dedup +=
        stats->ranges_after_dedup;
    solver->profile.page_plan_logical_bytes += stats->logical_bytes;
    solver->profile.page_plan_aligned_bytes += stats->aligned_bytes;
    solver->profile.page_plan_overread_bytes += stats->overread_bytes;
    solver->profile.page_plan_not_applicable +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE];
    solver->profile.page_plan_allocation_refused +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_ALLOCATION];
    solver->profile.page_plan_source_mismatch +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_SOURCE_MISMATCH];
    solver->profile.page_plan_invalid_range +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE];
    solver->profile.page_plan_byte_budget_refused +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET];
    solver->profile.page_plan_range_capacity_refused +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_RANGE_CAPACITY];
    solver->profile.page_plan_service_refused +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED];
    solver->profile.page_plan_service_errors +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR];
    solver->profile.page_plan_cancelled +=
        stats->refusal_counts[SOLVER_CODEKD_PAGE_PLAN_CANCELLED];
    solver->profile.candidate_delivery_candidates +=
        packet->candidate_count;
    solver->profile.candidate_quad_submitted +=
        packet->candidate_quad_submitted;
    solver->profile.candidate_quad_ready +=
        packet->candidate_quad_ready;
    solver->profile.candidate_quad_fallback +=
        packet->candidate_quad_fallback;
    solver->profile.candidate_star_submitted +=
        packet->candidate_star_submitted;
    solver->profile.candidate_star_ready +=
        packet->candidate_star_ready;
    solver->profile.candidate_star_fallback +=
        packet->candidate_star_fallback;
    solver->profile.candidate_delivery_windows +=
        packet->candidate_delivery_windows;
    solver->profile.candidate_quad_ready_rows +=
        packet->candidate_quad_ready_rows;
    solver->profile.candidate_star_ready_rows +=
        packet->candidate_star_ready_rows;
    solver->profile.candidate_retired_rows +=
        packet->candidate_retired_rows;
    solver->profile.candidate_native_rows +=
        packet->candidate_native_rows;
    solver->profile.verification_page_queries +=
        packet->verification_page_queries;
    solver->profile.verification_page_queries_planned +=
        packet->verification_page_queries_planned;
    solver->profile.verification_page_prefixes +=
        packet->verification_page_prefixes;
    solver->profile.verification_page_submitted +=
        packet->verification_page_submitted;
    solver->profile.verification_page_ready +=
        packet->verification_page_ready;
    solver->profile.verification_page_fallback +=
        packet->verification_page_fallback;
    solver->profile.verification_page_ready_rows +=
        packet->verification_page_ready_rows;
    solver->profile.verification_page_ranges +=
        packet->verification_page_ranges;
    solver->profile.verification_page_logical_bytes +=
        packet->verification_page_logical_bytes;
    solver->profile.verification_page_aligned_bytes +=
        packet->verification_page_aligned_bytes;
    solver->profile.candidate_math_prepared +=
        packet->candidate_math_prepared;
    solver->profile.verification_score_batches_prepared +=
        packet->verification_score_batches_prepared;
    solver->profile.verification_score_contexts_prepared +=
        packet->verification_score_contexts_prepared;
    solver->profile.verification_score_batches_executed +=
        packet->verification_score_batches_executed;
    solver->profile.verification_score_contexts_completed +=
        packet->verification_score_contexts_completed;
    if (ULLONG_MAX -
            solver->profile.verification_score_work_units_completed <
        packet->verification_score_work_units_completed) {
        solver->profile.verification_score_work_units_completed =
            ULLONG_MAX;
    } else {
        solver->profile.verification_score_work_units_completed +=
            packet->verification_score_work_units_completed;
    }
    solver->profile.verification_score_fallback_batches +=
        packet->verification_score_fallback_batches;
    solver->profile.verification_score_stopped_batches +=
        packet->verification_score_stopped_batches;
    solver->profile.verification_score_wall_seconds +=
        packet->verification_score_wall_seconds;
    if (packet->detailed) {
        solver->profile.verify_wall_seconds +=
            packet->verification_score_wall_seconds;
        solver->profile.resolve_wall_seconds +=
            packet->verification_score_wall_seconds;
    }
}

/*
 * Return zero when the bounded output arena is ready, one for an exact native
 * fallback, and minus one for an invalid packet request.
 */
static int solver_codekd_search_packet_prepare(
    solver_codekd_search_packet_t* packet,
    solver_ab_descriptor_output_t* descriptors,
    const kdtree_t* tree,
    const quadfile_t* quads,
    startree_t* starkd,
    int dimquads,
    anbool use_radec,
    size_t candidate_budget_bytes,
    unsigned long long sequence,
    anbool detailed) {
    size_t metadata_bytes;
    size_t hit_bytes;

    if (!packet || !descriptors || !tree || !quads || !starkd ||
        !starkd->tree ||
        dimquads <= 0 || dimquads > DQMAX ||
        SOLVER_AB_DESCRIPTOR_CAPACITY >
            SIZE_MAX / sizeof(*packet->slots)) {
        return -1;
    }
    memset(packet, 0, sizeof(*packet));
    metadata_bytes =
        SOLVER_AB_DESCRIPTOR_CAPACITY * sizeof(*packet->slots);
    if (metadata_bytes >= SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES) {
        return 1;
    }
    hit_bytes =
        SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES - metadata_bytes;
    packet->hit_capacity =
        hit_bytes / (sizeof(*packet->inds) + sizeof(*packet->sdists));
    if (!packet->hit_capacity ||
        packet->hit_capacity > SIZE_MAX / sizeof(*packet->inds) ||
        packet->hit_capacity > SIZE_MAX / sizeof(*packet->sdists)) {
        return 1;
    }

    packet->slots = calloc(
        SOLVER_AB_DESCRIPTOR_CAPACITY, sizeof(*packet->slots));
    packet->inds = malloc(
        packet->hit_capacity * sizeof(*packet->inds));
    packet->sdists = malloc(
        packet->hit_capacity * sizeof(*packet->sdists));
    if (!packet->slots || !packet->inds || !packet->sdists ||
        solver_codekd_page_workspace_create(
            &packet->page_workspace)) {
        (void)solver_codekd_search_packet_cleanup(packet);
        return 1;
    }

    if (candidate_budget_bytes >=
        sizeof(*packet->candidate_records)) {
        packet->candidate_capacity = MIN(
            packet->hit_capacity,
            candidate_budget_bytes /
                sizeof(*packet->candidate_records));
        /*
         * candidate_records is a rolling scratch window, not storage for the
         * complete result stream. The owner retires this bounded prefix
         * before the same logical packet submits its next payload window.
         */
        packet->candidate_capacity = MIN(
            packet->candidate_capacity,
            (size_t)SOLVER_CANDIDATE_DELIVERY_LIMIT);
        if (packet->candidate_capacity >
            SIZE_MAX / sizeof(*packet->candidate_records)) {
            packet->candidate_capacity = 0U;
        }
        if (packet->candidate_capacity) {
            packet->candidate_records = calloc(
                packet->candidate_capacity,
                sizeof(*packet->candidate_records));
            if (!packet->candidate_records) {
                packet->candidate_capacity = 0U;
            } else {
                size_t record_bytes =
                    packet->candidate_capacity *
                        sizeof(*packet->candidate_records);

                packet->verify_query_budget =
                    candidate_budget_bytes > record_bytes
                    ? candidate_budget_bytes - record_bytes
                    : 0U;
            }
        }
    }

    packet->descriptors = descriptors;
    packet->tree = tree;
    packet->quads = quads;
    packet->starkd = starkd;
    packet->dimquads = dimquads;
    packet->use_radec = use_radec;
    packet->star_delivery_eligible =
        !starkd->tree->perm || starkd->inverse_perm != NULL;
    packet->candidate_star_delivery_disabled =
        !packet->star_delivery_eligible;
    packet->sequence = sequence;
    packet->detailed = detailed;
    packet->state = SOLVER_CODEKD_PACKET_ALLOCATED;
    return 0;
}

static index_shard_helper_task_status_t
solver_codekd_packet_generate_descriptors(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    index_shard_helper_task_status_t descriptor_status;

    if (!input || input_size != sizeof(*input) ||
        !packet || output_size != sizeof(*packet) ||
        !input->tree || packet->tree != input->tree ||
        !input->quads || packet->quads != input->quads ||
        !input->starkd || packet->starkd != input->starkd ||
        !packet->descriptors || !packet->slots ||
        !packet->inds || !packet->sdists ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_ALLOCATED) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    descriptor_status = solver_ab_descriptor_helper_execute(
        &input->descriptor,
        sizeof(input->descriptor),
        packet->descriptors,
        sizeof(*packet->descriptors));
    if (descriptor_status != INDEX_SHARD_HELPER_TASK_OK) {
        packet->state =
            descriptor_status == INDEX_SHARD_HELPER_TASK_STOPPED
                ? SOLVER_CODEKD_PACKET_STOPPED
                : SOLVER_CODEKD_PACKET_FAILED;
        return descriptor_status;
    }
    if (packet->descriptors->descriptor_count >
        SOLVER_AB_DESCRIPTOR_CAPACITY) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    packet->first = 0U;
    packet->count = packet->descriptors->descriptor_count;
    packet->next_descriptor = 0U;
    packet->page_stats.descriptors_total = packet->count;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static int solver_codekd_packet_mark_plan_owner_replay(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason) {
    size_t descriptor_index;

    if (!packet || !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->count) {
        return -1;
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
            SOLVER_CODEKD_RESULT_UNUSED) {
            return -1;
        }
        packet->slots[descriptor_index].state =
            SOLVER_CODEKD_RESULT_OWNER_REPLAY;
    }
    solver_codekd_page_plan_record_refusal(packet, reason);
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return 0;
}

/*
 * A permanent pre-submission refusal invalidates the current sealed plan and
 * makes another plan attempt uneconomic. Preserve the exact logical work by
 * replaying the planned prefix and every still-unplanned descriptor at owner
 * retirement.
 */
static int solver_codekd_packet_mark_outstanding_owner_replay(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason) {
    size_t descriptor_index;

    if (!packet || !packet->plan_complete ||
        packet->state != SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->next_descriptor ||
        packet->next_descriptor > packet->count) {
        return -1;
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
            SOLVER_CODEKD_RESULT_UNUSED) {
            return -1;
        }
    }
    for (descriptor_index = packet->plan_end;
         descriptor_index < packet->count;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state !=
                SOLVER_CODEKD_RESULT_UNUSED &&
            packet->slots[descriptor_index].state !=
                SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            return -1;
        }
    }
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->count;
         descriptor_index++) {
        if (packet->slots[descriptor_index].state ==
            SOLVER_CODEKD_RESULT_UNUSED) {
            packet->slots[descriptor_index].state =
                SOLVER_CODEKD_RESULT_OWNER_REPLAY;
        }
    }
    solver_codekd_page_plan_record_refusal(packet, reason);
    packet->next_descriptor = packet->count;
    packet->pending_descriptor_plan = FALSE;
    packet->pending_descriptor_raw_ranges = 0U;
    packet->pending_descriptor_logical_bytes = 0U;
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
    return 0;
}

/* Return one when submitted, zero on bounded capacity, and minus one when
 * the complete packet must use exact owner replay. */
static int solver_codekd_packet_submit_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    int submit_status;

    if (!packet || !packet->tree || !packet->tree->io ||
        !packet->tree->io_is_fitsbin || !packet->page_workspace ||
        packet->delivery_source || packet->delivery_ticket ||
        packet->state != SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE ||
        !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->next_descriptor ||
        packet->next_descriptor > packet->count ||
        !packet->plan_range_count ||
        packet->plan_range_count >
            packet->page_workspace->sealed_range_capacity) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    source = (fitsbin_t*)packet->tree->io;
    packet->delivery_source = source;
    errno = 0;
    submit_status = fitsbin_prefetch_ranges_submit(
        source,
        packet->page_workspace->sealed_ranges,
        packet->plan_range_count,
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
        &packet->delivery_ticket);
    if (submit_status > 0 && packet->delivery_ticket) {
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    if (solver_codekd_packet_mark_outstanding_owner_replay(
            packet,
            submit_status < 0
                ? SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR
                : SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
    }
    return -1;
}

static int solver_codekd_packet_begin_candidate_window(
    solver_codekd_search_packet_t* packet) {
    size_t descriptor_index;
    anbool found = FALSE;

    if (!packet || !packet->candidate_records ||
        !packet->candidate_capacity || !packet->candidate_count ||
        packet->candidate_cursor >= packet->candidate_count ||
        packet->candidate_window_count ||
        packet->candidate_window_offset ||
        packet->verify_sweep_reads ||
        packet->verify_sweep_buffers ||
        packet->verify_sweep_storage ||
        packet->verify_sweep_storage_bytes) {
        return -1;
    }
    for (descriptor_index = packet->retire_descriptor;
         descriptor_index < packet->count;
         descriptor_index++) {
        const solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];
        size_t slot_end;

        if (slot->state != SOLVER_CODEKD_RESULT_READY ||
            !slot->hit_count ||
            slot->hit_first > packet->hit_count ||
            (size_t)slot->hit_count >
                packet->hit_count - slot->hit_first) {
            continue;
        }
        slot_end = slot->hit_first + (size_t)slot->hit_count;
        if (packet->candidate_cursor >= slot->hit_first &&
            packet->candidate_cursor < slot_end) {
            found = TRUE;
            break;
        }
    }
    if (!found) {
        return -1;
    }
    packet->candidate_window_first = packet->candidate_cursor;
    packet->candidate_window_count = MIN(
        MIN(packet->candidate_capacity,
            (size_t)SOLVER_CANDIDATE_DELIVERY_LIMIT),
        packet->candidate_count - packet->candidate_cursor);
    if (!packet->candidate_window_count) {
        return -1;
    }
    packet->candidate_window_offset = 0U;
    packet->candidate_star_count = 0U;
    packet->candidate_quad_ready_count = 0U;
    packet->candidate_star_ready_count = 0U;
    packet->candidate_verify_query_count = 0U;
    packet->verify_plan_first = 0U;
    packet->verify_plan_end = 0U;
    packet->verify_topology_end = 0U;
    packet->verify_plan_range_count = 0U;
    packet->verify_plan_logical_bytes = 0U;
    packet->verify_query_bytes = 0U;
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_storage_bytes = 0U;
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->verify_plan_complete = FALSE;
    packet->verify_sweep_plan_complete = FALSE;
    packet->verify_pending_query_plan = FALSE;
    packet->verification_delivery_disabled = FALSE;
    packet->candidate_delivery_windows++;
    packet->state = SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY;
    return 0;
}

static void solver_codekd_packet_clear_candidate_window(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;

    if (!packet) {
        return;
    }
    assert(!packet->delivery_ticket);
    assert(!packet->delivery_source);
    solver_codekd_packet_clear_sweep_storage(packet);
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &packet->candidate_records[candidate_index]);
        verify_destroy_index_query(
            packet->candidate_records[
                candidate_index].verify_query);
        packet->candidate_records[
            candidate_index].verify_query = NULL;
        packet->candidate_records[
            candidate_index].verify_query_captured = FALSE;
    }
    packet->candidate_window_first = packet->candidate_cursor;
    packet->candidate_window_count = 0U;
    packet->candidate_window_offset = 0U;
    packet->candidate_star_count = 0U;
    packet->candidate_quad_ready_count = 0U;
    packet->candidate_star_ready_count = 0U;
    packet->candidate_verify_query_count = 0U;
    packet->verify_plan_first = 0U;
    packet->verify_plan_end = 0U;
    packet->verify_topology_end = 0U;
    packet->verify_plan_range_count = 0U;
    packet->verify_plan_logical_bytes = 0U;
    packet->verify_query_bytes = 0U;
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_storage_bytes = 0U;
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->verify_plan_complete = FALSE;
    packet->verify_sweep_plan_complete = FALSE;
    packet->verify_pending_query_plan = FALSE;
    packet->verification_delivery_disabled = FALSE;
}

/*
 * Build verification queries from already copied candidate data. The record
 * also retains the exact native gates and TAN result so ordered retirement
 * can reuse the calculation after validating the complete immutable input.
 */
static int solver_codekd_packet_build_verify_queries(
    solver_codekd_search_packet_t* packet,
    const solver_codekd_packet_task_input_t* input) {
    const solver_field_geometry_t* geometry;
    size_t candidate_index;

    if (!packet || !input || input->use_radec ||
        !packet->candidate_records ||
        !packet->candidate_window_count ||
        packet->candidate_star_ready_count !=
            packet->candidate_window_count ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return -1;
    }
    geometry = input->descriptor.field_geometry;
    if (!geometry || !geometry->fieldxy ||
        geometry->numxy <= 0) {
        return -1;
    }
    packet->candidate_verify_query_count = 0U;

    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        const solver_ab_descriptor_t* descriptor;
        double field_xy[DQMAX * 2];
        double corner[3];
        double scale;
        double arcsecperpix;
        double abscale;
        double cx;
        double cy;
        double radius;
        tan_t wcs;
        int star_index;

        record->plan_action = SOLVER_AB_CANDIDATE_SCALE_SKIP;
        record->candidate_prepared = FALSE;
        memset(record->verify_center, 0, sizeof(record->verify_center));
        record->verify_radius = 0.0;
        record->verify_radius2 = 0.0;
        if (record->descriptor_index >= packet->count) {
            return -1;
        }
        descriptor = &packet->descriptors->descriptors[
            record->descriptor_index];
        for (star_index = 0;
             star_index < packet->dimquads;
             star_index++) {
            int field_star = descriptor->stars[star_index];

            if (field_star < 0 ||
                field_star >= geometry->numxy) {
                return -1;
            }
            setx(
                field_xy,
                star_index,
                starxy_getx(geometry->fieldxy, field_star));
            sety(
                field_xy,
                star_index,
                starxy_gety(geometry->fieldxy, field_star));
            record->prepared_fieldstars[star_index] = field_star;
        }
        memcpy(
            record->prepared_fieldxy,
            field_xy,
            (size_t)packet->dimquads * 2U * sizeof(*field_xy));
        record->prepared_parity = descriptor->current_parity;

        abscale = square(distsq2rad(distsq(
            record->starxyz, record->starxyz + 3, 3))) /
            distsq(field_xy, field_xy + 2, 2);
        if (abscale > input->abscale_high ||
            abscale < input->abscale_low) {
            record->plan_action =
                SOLVER_AB_CANDIDATE_ABSCALE_SKIP;
            record->candidate_prepared = TRUE;
            packet->candidate_math_prepared++;
            continue;
        }
        if (fit_tan_wcs(
                record->starxyz,
                field_xy,
                packet->dimquads,
                &wcs,
                &scale)) {
            record->plan_action =
                SOLVER_AB_CANDIDATE_BAD_QUAD;
            record->candidate_prepared = TRUE;
            packet->candidate_math_prepared++;
            continue;
        }
        arcsecperpix = scale * 3600.0;
        memcpy(&record->prepared_wcs, &wcs, sizeof(wcs));
        record->prepared_scale = arcsecperpix;
        record->candidate_prepared = TRUE;
        packet->candidate_math_prepared++;
        if (arcsecperpix > input->funits_upper ||
            arcsecperpix < input->funits_lower) {
            continue;
        }

        record->plan_action = SOLVER_AB_CANDIDATE_VERIFY;
        cx = 0.5 * (input->field_minx + input->field_maxx);
        cy = 0.5 * (input->field_miny + input->field_maxy);
        tan_pixelxy2xyzarr(
            &wcs, cx, cy, record->verify_center);
        tan_pixelxy2xyzarr(
            &wcs, input->field_minx, input->field_miny, corner);
        radius = sqrt(distsq(
            record->verify_center, corner, 3));
        record->verify_radius = radius;
        record->verify_radius2 = square(radius);
        if (!isfinite(record->verify_radius2) ||
            record->verify_radius2 < 0.0) {
            record->verify_radius = 0.0;
            record->verify_radius2 = 0.0;
            record->verify_delivery_fallback = TRUE;
            continue;
        }
        packet->candidate_verify_query_count++;
        packet->verification_page_queries++;
    }
    return 0;
}

static anbool solver_codekd_packet_candidate_data_fully_resident(
    const solver_codekd_search_packet_t* packet) {
    fitsbin_t* star_source;

    if (!packet || !packet->quads || !packet->quads->fb ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin) {
        return FALSE;
    }
    star_source = (fitsbin_t*)packet->starkd->tree->io;
    return fitsbin_payload_is_fully_resident(packet->quads->fb) &&
        fitsbin_payload_is_fully_resident(star_source);
}

/*
 * Complete the CodeKD phase through one central transition. Resident
 * candidate payloads retain the existing verification-wave path. Otherwise,
 * the first bounded Quad/A/B window becomes independently deliverable.
 */
static int solver_codekd_packet_finish_codekd(
    solver_codekd_search_packet_t* packet) {
    if (!packet ||
        (packet->state != SOLVER_CODEKD_PACKET_DESCRIPTORS_READY &&
         packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY) ||
        packet->next_descriptor != packet->count ||
        packet->plan_complete || packet->delivery_ticket ||
        packet->delivery_source) {
        return -1;
    }
    if (packet->candidate_count) {
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? 0
            : -1;
    }
    if (!packet->hit_count || !packet->candidate_records ||
        !packet->candidate_capacity || packet->use_radec ||
        !packet->quads || !packet->starkd ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX ||
        solver_codekd_packet_candidate_data_fully_resident(packet)) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }

    packet->candidate_count = packet->hit_count;
    packet->candidate_cursor = 0U;
    packet->retire_descriptor = 0U;
    packet->retire_hit_offset = 0U;
    packet->retire_descriptor_started = FALSE;
    return solver_codekd_packet_begin_candidate_window(packet);
}

/*
 * Submit one bounded post-CodeKD payload stage. A permanent refusal preserves
 * the copied CodeKD results and falls back only candidate resolution.
 */
static int solver_codekd_packet_submit_candidate_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    int submit_status;

    if (!packet || packet->delivery_source || packet->delivery_ticket ||
        !packet->candidate_records || !packet->candidate_count ||
        packet->candidate_cursor >= packet->candidate_count ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }

    errno = 0;
    if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY) {
        if (packet->candidate_quad_delivery_disabled) {
            packet->candidate_quad_fallback++;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return -1;
        }
        if (!packet->quads || !packet->quads->fb) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return -1;
        }
        source = packet->quads->fb;
        submit_status = quadfile_prefetch_stars_submit(
            packet->quads,
            packet->inds + packet->candidate_window_first,
            (int)packet->candidate_window_count,
            &packet->delivery_ticket);
        if (submit_status > 0 && packet->delivery_ticket) {
            packet->delivery_source = source;
            packet->candidate_quad_submitted++;
            packet->state = SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED;
            return 1;
        }
        packet->delivery_ticket = NULL;
        if (!submit_status && !errno &&
            fitsbin_payload_is_fully_resident(source)) {
            packet->state = SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
            return 2;
        }
        if (!submit_status && errno == EAGAIN) {
            return 0;
        }
        packet->candidate_quad_fallback++;
        packet->candidate_quad_delivery_disabled = TRUE;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return -1;
    }

    if (packet->state != SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY ||
        packet->candidate_star_delivery_disabled ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->candidate_star_count ||
        packet->candidate_star_count > SOLVER_CANDIDATE_STAR_LIMIT) {
        packet->candidate_star_fallback++;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return -1;
    }
    source = (fitsbin_t*)packet->starkd->tree->io;
    submit_status = startree_prefetch_stars_ready_submit(
        packet->starkd,
        packet->candidate_starids,
        (int)packet->candidate_star_count,
        &packet->delivery_ticket);
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->delivery_source = source;
        packet->candidate_star_submitted++;
        packet->state = SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED;
        return 1;
    }
    packet->delivery_ticket = NULL;
    if (!submit_status && !errno &&
        fitsbin_payload_is_fully_resident(source)) {
        packet->state = SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY;
        return 2;
    }
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    packet->candidate_star_fallback++;
    packet->candidate_star_delivery_disabled = TRUE;
    packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
    return -1;
}

static int solver_codekd_packet_submit_verification_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    size_t candidate_index;
    anbool has_query = FALSE;
    int submit_status;

    if (!packet || packet->delivery_source ||
        packet->delivery_ticket ||
        !packet->candidate_records ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY ||
        packet->verification_delivery_disabled ||
        packet->verify_plan_complete ||
        packet->candidate_window_offset >=
            packet->candidate_window_count ||
        packet->candidate_star_ready_count !=
            packet->candidate_window_count ||
        !packet->starkd || !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }
    for (candidate_index = packet->candidate_window_offset;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        const solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];

        if (record->verify_delivery_fallback) {
            solver_codekd_packet_disable_verification_delivery(
                packet);
            return -1;
        }
        if (record->plan_action !=
            SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        has_query = TRUE;
        break;
    }
    if (!has_query) {
        packet->verify_plan_first =
            packet->candidate_window_offset;
        packet->verify_plan_end =
            packet->candidate_window_count;
        packet->verify_topology_end =
            packet->candidate_window_count;
        packet->verify_plan_complete = TRUE;
        packet->verification_page_prefixes++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
        return 2;
    }

    source = (fitsbin_t*)packet->starkd->tree->io;
    packet->verify_delivery_budget =
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES;
    packet->delivery_source = source;
    packet->state = SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED;
    errno = 0;
    submit_status = fitsbin_prefetch_ranges_planned_submit(
        source,
        solver_codekd_packet_plan_verification_pages,
        packet,
        packet->verify_delivery_budget,
        &packet->delivery_ticket);
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->verification_page_submitted++;
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && !errno &&
        fitsbin_payload_is_fully_resident(source)) {
        packet->verify_plan_first =
            packet->candidate_window_offset;
        packet->verify_plan_end =
            packet->candidate_window_count;
        packet->verify_topology_end =
            packet->candidate_window_count;
        packet->verify_plan_complete = TRUE;
        packet->verification_page_prefixes++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
        return 2;
    }
    if (!submit_status && errno == EAGAIN) {
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY;
        return 0;
    }
    solver_codekd_packet_disable_verification_delivery(packet);
    return -1;
}

static int solver_codekd_packet_submit_sweep_pages(
    solver_codekd_search_packet_t* packet) {
    fitsbin_t* source;
    size_t offset;
    size_t range_index;
    anbool direct = FALSE;
    int submit_status;

    if (!packet || packet->delivery_source ||
        packet->delivery_ticket || !packet->starkd ||
        !packet->starkd->tree ||
        !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->page_workspace ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY ||
        !packet->verify_plan_complete ||
        !packet->verify_sweep_plan_complete ||
        !packet->verify_sweep_range_count ||
        packet->verify_sweep_range_count >
            packet->page_workspace->sealed_range_capacity ||
        !packet->verify_sweep_aligned_bytes ||
        packet->verify_sweep_aligned_bytes >
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES) {
        return -1;
    }
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return -1;
    }

    source = (fitsbin_t*)packet->starkd->tree->io;
    if (fitsbin_payload_is_fully_resident(source)) {
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
        return 2;
    }

    if (packet->verify_sweep_reads ||
        packet->verify_sweep_buffers ||
        packet->verify_sweep_storage ||
        packet->verify_sweep_storage_bytes) {
        if (!packet->verify_sweep_reads ||
            !packet->verify_sweep_buffers ||
            !packet->verify_sweep_storage ||
            packet->verify_sweep_storage_bytes !=
                packet->verify_sweep_aligned_bytes) {
            solver_codekd_packet_clear_sweep_storage(packet);
        } else {
            direct = TRUE;
        }
    }
    if (!direct &&
        packet->verify_sweep_range_count <=
            FITSBIN_PREAD_ASYNC_RANGE_LIMIT &&
        packet->verify_query_bytes <=
            packet->verify_query_budget &&
        packet->verify_sweep_aligned_bytes <=
            packet->verify_query_budget -
                packet->verify_query_bytes) {
        if (packet->verify_sweep_range_count <=
                SIZE_MAX / sizeof(*packet->verify_sweep_reads) &&
            packet->verify_sweep_range_count <=
                SIZE_MAX / sizeof(*packet->verify_sweep_buffers)) {
            packet->verify_sweep_reads = calloc(
                packet->verify_sweep_range_count,
                sizeof(*packet->verify_sweep_reads));
            packet->verify_sweep_buffers = calloc(
                packet->verify_sweep_range_count,
                sizeof(*packet->verify_sweep_buffers));
            packet->verify_sweep_storage = malloc(
                packet->verify_sweep_aligned_bytes);
        }
        if (packet->verify_sweep_reads &&
            packet->verify_sweep_buffers &&
            packet->verify_sweep_storage) {
            offset = 0U;
            direct = TRUE;
            for (range_index = 0U;
                 range_index < packet->verify_sweep_range_count;
                 range_index++) {
                const fitsbin_prefetch_range_t* range =
                    &packet->page_workspace->
                        sealed_ranges[range_index];

                if (!range->data || !range->size ||
                    range->size >
                        packet->verify_sweep_aligned_bytes -
                            offset) {
                    direct = FALSE;
                    break;
                }
                packet->verify_sweep_reads[range_index].data =
                    range->data;
                packet->verify_sweep_reads[range_index].size =
                    range->size;
                packet->verify_sweep_reads[
                    range_index].logical_size = range->size;
                packet->verify_sweep_reads[
                    range_index].destination =
                        packet->verify_sweep_storage + offset;
                packet->verify_sweep_buffers[
                    range_index].mapping_data = range->data;
                packet->verify_sweep_buffers[
                    range_index].size = range->size;
                packet->verify_sweep_buffers[
                    range_index].bytes =
                        packet->verify_sweep_storage + offset;
                offset += range->size;
            }
            if (offset != packet->verify_sweep_aligned_bytes) {
                direct = FALSE;
            }
        }
        if (!direct) {
            solver_codekd_packet_clear_sweep_storage(packet);
        } else {
            packet->verify_sweep_storage_bytes =
                packet->verify_sweep_aligned_bytes;
        }
    }

    errno = 0;
    if (direct) {
        submit_status = fitsbin_pread_mapped_ranges_submit(
            source,
            packet->verify_sweep_reads,
            packet->verify_sweep_range_count,
            SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
            FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
            &packet->delivery_ticket);
        if (submit_status > 0 && packet->delivery_ticket) {
            packet->delivery_source = source;
            packet->verification_page_submitted++;
            packet->state =
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED;
            return 1;
        }
        packet->delivery_ticket = NULL;
        if (!submit_status && errno == EAGAIN) {
            return 0;
        }
        /*
         * A direct translator or allocation refusal is not a reason to
         * discard the already complete mapped-page plan. Preserve the
         * previous completion provider before falling back to native mmap.
         */
        solver_codekd_packet_clear_sweep_storage(packet);
        errno = 0;
    }
    submit_status = fitsbin_prefetch_ranges_submit(
        source,
        packet->page_workspace->sealed_ranges,
        packet->verify_sweep_range_count,
        SOLVER_CODEKD_DELIVERY_BUDGET_BYTES,
        &packet->delivery_ticket);
    if (submit_status > 0 && packet->delivery_ticket) {
        packet->delivery_source = source;
        packet->verification_page_submitted++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED;
        return 1;
    }
    packet->delivery_ticket = NULL;
    packet->delivery_source = NULL;
    if (!submit_status && errno == EAGAIN) {
        return 0;
    }
    solver_codekd_packet_clear_sweep_storage(packet);
    packet->verification_page_fallback++;
    packet->state = SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
    return 2;
}

/* Return zero while pending, one when compute-ready, two when replanned work
 * remains, three when the logical packet is complete, and minus one on error
 * or cooperative stop. */
static int solver_codekd_packet_collect_pages(
    solver_codekd_search_packet_t* packet) {
    solver_codekd_search_packet_state_t io_state;
    int ticket_result = 0;
    int ticket_errno;
    int poll_status;

    if (!packet || !packet->delivery_source ||
        !packet->delivery_ticket) {
        return -1;
    }
    errno = 0;
    poll_status = fitsbin_payload_io_ticket_poll(
        packet->delivery_source,
        packet->delivery_ticket,
        &ticket_result);
    ticket_errno = errno;
    if (!poll_status) {
        return 0;
    }
    if (poll_status < 0) {
        if (!solver_codekd_search_packet_release_ticket(packet)) {
            solver_codekd_packet_clear_sweep_storage(packet);
        }
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return -1;
    }
    /*
     * Verification planned-ticket output is owned by the I/O lane until
     * terminal publication. Poll's payload_io mutex edge makes those state
     * writes visible here. Initial CodeKD tickets use an immutable sealed
     * plan, so their packet remains PAGE_PLAN_COMPLETE while queued.
     */
    io_state = packet->state;
    if (fitsbin_payload_io_ticket_destroy_checked(
            packet->delivery_ticket)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return -1;
    }
    packet->delivery_ticket = NULL;
    if ((ticket_result == 0 && ticket_errno == ECANCELED) ||
        index_shard_worker_stop_requested()) {
        packet->delivery_source = NULL;
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        solver_codekd_page_plan_record_refusal(
            packet, SOLVER_CODEKD_PAGE_PLAN_CANCELLED);
        return -1;
    }
    if (ticket_result > 0 &&
        io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        packet->state = SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state == SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        packet->state = SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE &&
        packet->verify_plan_complete &&
        packet->verify_plan_first ==
            packet->candidate_window_offset &&
        packet->verify_plan_end >
            packet->verify_plan_first &&
            packet->verify_plan_end <=
            packet->candidate_window_count) {
        packet->delivery_source = NULL;
        packet->verification_page_ready++;
        packet->verification_page_ready_rows +=
            packet->verify_plan_end -
                packet->verify_plan_first;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
        return 1;
    }
    if (ticket_result > 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED &&
        packet->verify_plan_complete &&
        packet->verify_sweep_plan_complete &&
        packet->verify_sweep_range_count) {
        packet->delivery_source = NULL;
        packet->verification_page_ready++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 &&
        io_state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED) {
        packet->delivery_source = NULL;
        solver_codekd_packet_clear_sweep_storage(packet);
        packet->verification_page_fallback++;
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 &&
        (io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED ||
         io_state == SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED)) {
        packet->delivery_source = NULL;
        if (io_state == SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
        } else {
            packet->candidate_star_fallback++;
            packet->candidate_star_delivery_disabled = TRUE;
        }
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 3;
    }
    if (io_state == SOLVER_CODEKD_PACKET_RESULTS_READY &&
        packet->verification_delivery_disabled) {
        packet->delivery_source = NULL;
        return 3;
    }
    if (ticket_result < 0 &&
        (io_state == SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED ||
         io_state ==
             SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE)) {
        packet->delivery_source = NULL;
        if (packet->state == SOLVER_CODEKD_PACKET_FAILED) {
            return -1;
        }
        solver_codekd_packet_disable_verification_delivery(packet);
        return 3;
    }
    if (ticket_result > 0 &&
        packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY &&
        !packet->plan_complete) {
        packet->delivery_source = NULL;
        return 3;
    }
    if (ticket_result > 0 && packet->plan_complete &&
        packet->plan_range_count &&
        packet->state == SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
        packet->state = SOLVER_CODEKD_PACKET_COMPUTE_READY;
        return 1;
    }
    if (ticket_result < 0 && packet->plan_complete &&
        packet->state == SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
        if (solver_codekd_packet_mark_plan_owner_replay(
                packet, SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return -1;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? 3
            : 2;
    }
    packet->delivery_source = NULL;
    packet->state = SOLVER_CODEKD_PACKET_FAILED;
    return -1;
}

static int solver_codekd_packet_cancel_pages(
    solver_codekd_search_packet_t* packet) {
    if (!packet || !packet->delivery_ticket) {
        return 0;
    }
    return fitsbin_payload_io_ticket_cancel_async(
        packet->delivery_ticket);
}

static index_shard_helper_task_status_t
solver_codekd_packet_decode_quad_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    size_t stars_per_candidate;
    size_t candidate_index;
    size_t descriptor_index;

    if (!packet || !more_work || !packet->quads || !packet->starkd ||
        !packet->candidate_records ||
        packet->state != SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY ||
        !packet->candidate_count ||
        !packet->candidate_capacity ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_first >
            packet->candidate_count ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        packet->candidate_quad_ready_count ||
        packet->candidate_star_ready_count ||
        packet->use_radec ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }

    /*
     * Phase admission already identified a dense result stream. Populate
     * every quad-star row in one bounded ticket, but preserve the owner's
     * native A/B scale gate and logical consumption order.
     */
    stars_per_candidate = (size_t)packet->dimquads;
    if (packet->candidate_window_count >
        SOLVER_CANDIDATE_STAR_LIMIT / stars_per_candidate) {
        packet->candidate_quad_fallback++;
        packet->candidate_quad_delivery_disabled = TRUE;
        packet->candidate_star_delivery_disabled = TRUE;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->candidate_star_count =
        packet->candidate_window_count * stars_per_candidate;
    descriptor_index = packet->retire_descriptor;
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        size_t global_candidate =
            packet->candidate_window_first + candidate_index;
        size_t star_index;

        while (descriptor_index < packet->count) {
            const solver_codekd_result_slot_t* slot =
                &packet->slots[descriptor_index];
            size_t slot_end;

            if (slot->state != SOLVER_CODEKD_RESULT_READY ||
                !slot->hit_count) {
                descriptor_index++;
                continue;
            }
            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first) {
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            slot_end = slot->hit_first + (size_t)slot->hit_count;
            if (global_candidate < slot->hit_first) {
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            if (global_candidate < slot_end) {
                break;
            }
            descriptor_index++;
        }
        if (descriptor_index >= packet->count) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }

        memset(record, 0, sizeof(*record));
        record->descriptor_index = descriptor_index;
        record->quadid = packet->inds[global_candidate];
        if (quadfile_get_stars(
                packet->quads,
                record->quadid,
                record->stars)) {
            packet->candidate_quad_fallback++;
            packet->candidate_quad_delivery_disabled = TRUE;
            packet->candidate_star_delivery_disabled = TRUE;
            packet->candidate_star_count = 0U;
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return INDEX_SHARD_HELPER_TASK_OK;
        }
        for (star_index = 0U;
             star_index < stars_per_candidate;
             star_index++) {
            packet->candidate_starids[
                candidate_index * stars_per_candidate + star_index] =
                record->stars[star_index];
        }
    }
    packet->candidate_quad_ready_count =
        packet->candidate_window_count;
    packet->candidate_quad_ready_rows +=
        packet->candidate_window_count;
    packet->candidate_quad_ready++;
    if (!packet->star_delivery_eligible ||
        packet->candidate_star_delivery_disabled) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->state = SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY;
    *more_work = TRUE;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static index_shard_helper_task_status_t
solver_codekd_packet_copy_star_ready(
    const solver_codekd_packet_task_input_t* input,
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    size_t stars_per_candidate;
    size_t candidate_index;

    if (!input || !packet || !more_work || !packet->starkd ||
        packet->state != SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY ||
        !packet->candidate_records || !packet->candidate_count ||
        !packet->candidate_capacity ||
        !packet->candidate_window_count ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_first >
            packet->candidate_count ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_first != packet->candidate_cursor ||
        packet->candidate_window_offset ||
        packet->candidate_quad_ready_count !=
            packet->candidate_window_count ||
        packet->candidate_star_ready_count ||
        !packet->candidate_star_count ||
        packet->candidate_star_count > SOLVER_CANDIDATE_STAR_LIMIT ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }

    stars_per_candidate = (size_t)packet->dimquads;
    if (packet->candidate_star_count !=
        packet->candidate_window_count * stars_per_candidate) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        size_t star_index;

        for (star_index = 0U;
             star_index < stars_per_candidate;
             star_index++) {
            if (startree_get_ready(
                    packet->starkd,
                    (int)record->stars[star_index],
                    record->starxyz + 3U * star_index)) {
                packet->candidate_star_fallback++;
                packet->candidate_star_delivery_disabled = TRUE;
                packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
                return INDEX_SHARD_HELPER_TASK_OK;
            }
        }
    }
    packet->candidate_star_ready_count =
        packet->candidate_window_count;
    packet->candidate_star_ready_rows +=
        packet->candidate_window_count;
    packet->candidate_star_ready++;
    if (solver_codekd_packet_build_verify_queries(
            packet, input)) {
        for (candidate_index = 0U;
             candidate_index < packet->candidate_window_count;
             candidate_index++) {
            packet->candidate_records[
                candidate_index].candidate_prepared = FALSE;
        }
        packet->verification_delivery_disabled = TRUE;
        packet->verification_page_fallback++;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    packet->state = packet->candidate_verify_query_count
        ? SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY
        : SOLVER_CODEKD_PACKET_RESULTS_READY;
    if (packet->state == SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY) {
        *more_work = TRUE;
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

/*
 * Execute each native StarKD query once after its exact DATA/PERM pages are
 * ready. Retain the logical result order, then derive only the mapped sweep
 * pages that those results will dereference. Physical page keys may be
 * deduplicated and sorted; logical query results must never be reordered.
 */
static index_shard_helper_task_status_t
solver_codekd_packet_query_verification_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    solver_codekd_page_workspace_t* workspace;
    fitsbin_t* source;
    size_t plan_end;
    size_t segment_end;
    size_t candidate_index;
    size_t logical_bytes = 0U;
    size_t unique_pages = 0U;
    size_t range_count = 0U;
    size_t aligned_bytes = 0U;
    size_t captured_queries = 0U;
    int seal_status;

    if (!packet || !more_work || !packet->starkd ||
        !packet->starkd->tree || !packet->starkd->tree->io ||
        !packet->starkd->tree->io_is_fitsbin ||
        !packet->candidate_records || !packet->page_workspace ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY ||
        !packet->verify_plan_complete ||
        packet->verify_plan_first !=
            packet->candidate_window_offset ||
        packet->verify_plan_end <= packet->verify_plan_first ||
        packet->verify_plan_end >
            packet->candidate_window_count ||
        packet->verify_topology_end !=
            packet->verify_plan_end) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }
    solver_codekd_packet_clear_sweep_storage(packet);

    workspace = packet->page_workspace;
    source = (fitsbin_t*)packet->starkd->tree->io;
    /*
     * The descriptor scratch is now reused for exact sweep pages. Any
     * topology plan that was deferred beyond this prefix must be rebuilt
     * when owner retirement advances to it.
     */
    packet->verify_pending_query_plan = FALSE;
    packet->verify_pending_query_index = 0U;
    packet->verify_pending_query_raw_ranges = 0U;
    packet->verify_pending_query_logical_bytes = 0U;
    solver_codekd_page_set_reset(&workspace->descriptor);
    solver_codekd_page_set_reset(&workspace->group);
    plan_end = MIN(
        packet->verify_plan_end,
        packet->verify_plan_first +
            (size_t)SOLVER_VERIFY_QUERY_LOOKAHEAD);
    segment_end = packet->verify_plan_first;
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_plan_complete = FALSE;

    for (candidate_index = packet->verify_plan_first;
         candidate_index < plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        solver_codekd_page_plan_t plan;
        verify_index_query_t* query = NULL;
        size_t query_logical_bytes = 0U;
        size_t query_bytes;
        size_t result_index;
        anbool query_plan_failed = FALSE;

        if (record->plan_action !=
            SOLVER_AB_CANDIDATE_VERIFY) {
            segment_end = candidate_index + 1U;
            continue;
        }
        if (record->verify_query_captured) {
            if (!record->verify_query) {
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            captured_queries++;
            segment_end = candidate_index + 1U;
            continue;
        }
        if (index_shard_worker_stop_requested()) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        query = record->verify_query;
        if (!query) {
            if (verify_query_hit(
                    packet->starkd,
                    record->verify_center,
                    record->verify_radius2,
                    &query) ||
                !query) {
                segment_end = candidate_index + 1U;
                continue;
            }
            query_bytes = verify_index_query_bytes(query);
            if (query_bytes == SIZE_MAX ||
                query_bytes > packet->verify_query_budget ||
                packet->verify_query_bytes >
                    packet->verify_query_budget -
                        query_bytes) {
                verify_destroy_index_query(query);
                segment_end = candidate_index + 1U;
                continue;
            }
            record->verify_query = query;
            packet->verify_query_bytes += query_bytes;
        }

        solver_codekd_page_set_reset(&workspace->descriptor);
        memset(&plan, 0, sizeof(plan));
        plan.source = source;
        plan.workspace = workspace;
        plan.enabled = TRUE;
        for (result_index = 0U;
             result_index < verify_index_query_count(query);
             result_index++) {
            kdtree_prefetch_hint_t hint;
            const void* data = NULL;
            size_t size = 0U;
            int emit_status;

            if (!(result_index & 255U) &&
                index_shard_worker_stop_requested()) {
                packet->state =
                    SOLVER_CODEKD_PACKET_STOPPED;
                return INDEX_SHARD_HELPER_TASK_STOPPED;
            }
            if (verify_index_query_sweep_range(
                    packet->starkd,
                    query,
                    result_index,
                    &data,
                    &size)) {
                query_plan_failed = TRUE;
                break;
            }
            memset(&hint, 0, sizeof(hint));
            hint.mapping = source;
            hint.address = data;
            hint.length = size;
            hint.kind = KDTREE_PREFETCH_ARRAY_DATA;
            hint.priority = KDTREE_PREFETCH_PRIORITY_LEAF;
            emit_status =
                solver_codekd_page_plan_emit(&plan, &hint);
            if (emit_status !=
                KDTREE_PREFETCH_EMIT_CONTINUE) {
                query_plan_failed = TRUE;
                break;
            }
            solver_codekd_page_plan_add_size(
                &query_logical_bytes, size);
        }
        if (query_plan_failed) {
            if (segment_end > packet->verify_plan_first) {
                break;
            }
            packet->verify_plan_end = candidate_index + 1U;
            packet->verification_page_fallback++;
            packet->state =
                SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
            *more_work = TRUE;
            return INDEX_SHARD_HELPER_TASK_OK;
        }

        seal_status = solver_codekd_page_set_union_fits(workspace);
        if (seal_status <= 0) {
            if (segment_end > packet->verify_plan_first) {
                break;
            }
            packet->verify_plan_end = candidate_index + 1U;
            packet->verification_page_fallback++;
            packet->state =
                SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
            *more_work = TRUE;
            return INDEX_SHARD_HELPER_TASK_OK;
        }
        if (solver_codekd_page_set_merge_descriptor(workspace)) {
            goto sweep_fallback;
        }
        solver_codekd_page_plan_add_size(
            &logical_bytes, query_logical_bytes);
        segment_end = candidate_index + 1U;
        if (workspace->group.count >= workspace->page_limit) {
            break;
        }
    }

    if (segment_end <= packet->verify_plan_first) {
        goto sweep_fallback;
    }
    if (!workspace->group.count) {
        packet->verify_plan_end = segment_end;
        packet->verify_sweep_range_count = 0U;
        packet->verify_sweep_aligned_bytes = 0U;
        packet->verify_sweep_plan_complete = FALSE;
        packet->state = captured_queries
            ? SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER
            : SOLVER_CODEKD_PACKET_RESULTS_READY;
        *more_work = captured_queries != 0U;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    seal_status = solver_codekd_page_plan_seal_union(
        workspace,
        FALSE,
        &unique_pages,
        &range_count,
        &aligned_bytes);
    if (seal_status) {
        goto sweep_fallback;
    }
    packet->verify_plan_end = segment_end;
    packet->verify_sweep_range_count = range_count;
    packet->verify_sweep_aligned_bytes = aligned_bytes;
    packet->verify_sweep_plan_complete = TRUE;
    packet->verification_page_ranges += range_count;
    packet->verification_page_logical_bytes =
        solver_ab_saturating_add(
            packet->verification_page_logical_bytes,
            (unsigned long long)logical_bytes);
    packet->verification_page_aligned_bytes =
        solver_ab_saturating_add(
            packet->verification_page_aligned_bytes,
            (unsigned long long)aligned_bytes);
    packet->state = range_count
        ? SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY
        : SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY;
    *more_work = TRUE;
    return INDEX_SHARD_HELPER_TASK_OK;

sweep_fallback:
    packet->verify_sweep_range_count = 0U;
    packet->verify_sweep_aligned_bytes = 0U;
    packet->verify_sweep_plan_complete = FALSE;
    packet->verification_page_fallback++;
    packet->state = SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY;
    *more_work = TRUE;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static index_shard_helper_task_status_t
solver_codekd_packet_capture_sweep_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work) {
    size_t candidate_index;
    size_t captured_queries = 0U;
    anbool direct;
    anbool storage_valid = TRUE;
    anbool capture_failed = FALSE;

    if (!packet || !more_work || !packet->starkd ||
        !packet->candidate_records ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY ||
        !packet->verify_plan_complete ||
        !packet->verify_sweep_plan_complete ||
        packet->verify_plan_first !=
            packet->candidate_window_offset ||
        packet->verify_plan_end <= packet->verify_plan_first ||
        packet->verify_plan_end >
            packet->candidate_window_count ||
        packet->verify_topology_end <
            packet->verify_plan_end ||
        packet->verify_topology_end >
            packet->candidate_window_count) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    *more_work = FALSE;
    direct = packet->verify_sweep_storage ||
        packet->verify_sweep_buffers ||
        packet->verify_sweep_reads ||
        packet->verify_sweep_storage_bytes;
    if (direct &&
        (!packet->verify_sweep_storage ||
         !packet->verify_sweep_buffers ||
         !packet->verify_sweep_reads ||
         packet->verify_sweep_storage_bytes !=
             packet->verify_sweep_aligned_bytes)) {
        storage_valid = FALSE;
        capture_failed = TRUE;
    }
    for (candidate_index = packet->verify_plan_first;
         candidate_index < packet->verify_plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];

        if (index_shard_worker_stop_requested()) {
            solver_codekd_packet_clear_sweep_storage(packet);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        if (record->verify_query_captured) {
            if (!record->verify_query) {
                solver_codekd_packet_clear_sweep_storage(packet);
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            captured_queries++;
            continue;
        }
        if (record->verify_query &&
            (!direct || storage_valid)) {
            int capture_status = direct
                ? verify_index_query_capture_sweep_buffers(
                    packet->starkd,
                    record->verify_query,
                    packet->verify_sweep_buffers,
                    packet->verify_sweep_range_count)
                : verify_index_query_capture_sweep(
                    packet->starkd,
                    record->verify_query);

            if (capture_status) {
                capture_failed = TRUE;
            } else {
                record->verify_query_captured = TRUE;
                captured_queries++;
            }
        }
    }
    solver_codekd_packet_clear_sweep_storage(packet);
    if (index_shard_worker_stop_requested()) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }
    if (capture_failed) {
        packet->verification_page_fallback++;
    }
    packet->state = captured_queries
        ? SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER
        : SOLVER_CODEKD_PACKET_RESULTS_READY;
    *more_work =
        packet->state == SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER;
    return INDEX_SHARD_HELPER_TASK_OK;
}

static int solver_codekd_packet_build_verification_match(
    const solver_codekd_packet_task_input_t* input,
    const solver_codekd_search_packet_t* packet,
    const solver_candidate_delivery_record_t* record,
    size_t candidate_index,
    MatchObj* match,
    double* verify_pix2,
    double* logaccept) {
    const solver_codekd_verification_snapshot_t* snapshot;
    const solver_codekd_result_slot_t* slot;
    size_t global_candidate;
    size_t slot_end;
    double radius2;
    int star_index;

    if (!input || !packet || !record || !match ||
        !verify_pix2 || !logaccept ||
        candidate_index >= packet->candidate_window_count ||
        packet->candidate_window_first >
            packet->candidate_count ||
        candidate_index >=
            packet->candidate_count -
                packet->candidate_window_first ||
        record->descriptor_index >= packet->count ||
        record->plan_action != SOLVER_AB_CANDIDATE_VERIFY ||
        !record->candidate_prepared ||
        packet->dimquads <= 0 || packet->dimquads > DQMAX) {
        return -1;
    }
    snapshot = &input->verification;
    if (!snapshot->enabled || !snapshot->field ||
        !isfinite(record->prepared_scale) ||
        record->prepared_scale == 0.0 ||
        !isfinite(record->verify_radius) ||
        record->verify_radius < 0.0 ||
        !isfinite(snapshot->verify_pix) ||
        !isfinite(snapshot->index_jitter) ||
        !isfinite(snapshot->logaccept) ||
        !isfinite(snapshot->logratio_bail_threshold)) {
        return -1;
    }
    global_candidate =
        packet->candidate_window_first + candidate_index;
    slot = &packet->slots[record->descriptor_index];
    if (slot->state != SOLVER_CODEKD_RESULT_READY ||
        !slot->hit_count ||
        slot->hit_first > packet->hit_count ||
        (size_t)slot->hit_count >
            packet->hit_count - slot->hit_first) {
        return -1;
    }
    slot_end = slot->hit_first + (size_t)slot->hit_count;
    if (global_candidate < slot->hit_first ||
        global_candidate >= slot_end ||
        packet->inds[global_candidate] != record->quadid) {
        return -1;
    }
    radius2 = square(record->verify_radius);
    if (memcmp(&radius2, &record->verify_radius2, sizeof(radius2))) {
        return -1;
    }

    memcpy(match, &snapshot->match_template, sizeof(*match));
    memcpy(&match->wcstan, &record->prepared_wcs,
           sizeof(match->wcstan));
    match->wcs_valid = TRUE;
    match->code_err = packet->sdists[global_candidate];
    match->scale = record->prepared_scale;
    match->parity = record->prepared_parity;
    match->quad_npeers = slot->hit_count;
    match->quadno = record->quadid;
    match->dimquads = packet->dimquads;
    for (star_index = 0;
         star_index < packet->dimquads;
         star_index++) {
        match->star[star_index] = record->stars[star_index];
        match->field[star_index] =
            record->prepared_fieldstars[star_index];
        match->ids[star_index] = 0;
    }
    memcpy(match->quadpix, record->prepared_fieldxy,
           (size_t)packet->dimquads * 2U *
               sizeof(*match->quadpix));
    memcpy(match->quadxyz, record->starxyz,
           (size_t)packet->dimquads * 3U *
               sizeof(*match->quadxyz));
    memcpy(match->center, record->verify_center,
           sizeof(match->center));
    match->radius = record->verify_radius;
    match->radius_deg = dist2deg(match->radius);
    match->indexid = snapshot->indexid;
    match->healpix = snapshot->healpix;
    match->hpnside = snapshot->hpnside;
    match->wcstan.imagew = input->field_maxx;
    match->wcstan.imageh = input->field_maxy;
    *verify_pix2 = square(snapshot->verify_pix) +
        square(snapshot->index_jitter / match->scale);
    *logaccept = snapshot->logaccept;
    return isfinite(*verify_pix2) ? 0 : -1;
}

static void solver_codekd_packet_discard_prepared_verification(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;

    if (!packet || !packet->candidate_records) {
        return;
    }
    for (candidate_index = 0U;
         candidate_index < packet->candidate_window_count;
         candidate_index++) {
        solver_codekd_record_clear_prepared_verification(
            &packet->candidate_records[candidate_index]);
    }
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
}

static int solver_codekd_packet_retired_verification_complete(
    solver_codekd_search_packet_t* packet) {
    size_t candidate_index;

    if (!packet || !packet->candidate_records ||
        !packet->verify_plan_complete ||
        packet->verify_plan_end > packet->candidate_window_count) {
        return -1;
    }
    for (candidate_index = packet->verify_plan_first;
         candidate_index < packet->verify_plan_end;
         candidate_index++) {
        const solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];

        if (record->prepared_verification ||
            record->verification_score_ready ||
            record->prepared_score.complete ||
            record->prepared_score.theta ||
            record->prepared_score.allodds) {
            return -1;
        }
    }
    packet->verify_prepared_count = 0U;
    packet->verify_prepared_retained_bytes = 0U;
    packet->verify_prepared_transient_bytes = 0U;
    return 0;
}

/*
 * Saturating score-work accounting is diagnostic only; it must never alter
 * packet control flow.
 */
static void solver_codekd_packet_add_verification_score_work(
    solver_codekd_search_packet_t* packet,
    unsigned long long work_units) {
    if (!packet) {
        return;
    }
    if (packet->verification_score_work_units_completed ==
            ULLONG_MAX ||
        work_units == ULLONG_MAX ||
        ULLONG_MAX -
            packet->verification_score_work_units_completed <
                work_units) {
        packet->verification_score_work_units_completed =
            ULLONG_MAX;
    } else {
        packet->verification_score_work_units_completed +=
            work_units;
    }
}

static void solver_codekd_packet_destroy_score_array(
    verify_prepared_score_t* scores,
    size_t count) {
    size_t score_index;

    if (!scores) {
        return;
    }
    for (score_index = 0U; score_index < count; score_index++) {
        verify_destroy_prepared_score(&scores[score_index]);
    }
    free(scores);
}

typedef enum solver_codekd_verification_window_status {
    SOLVER_CODEKD_VERIFY_WINDOW_ERROR = -1,
    SOLVER_CODEKD_VERIFY_WINDOW_UNAVAILABLE = 0,
    SOLVER_CODEKD_VERIFY_WINDOW_SCORED = 1,
    SOLVER_CODEKD_VERIFY_WINDOW_STOPPED = 2,
    SOLVER_CODEKD_VERIFY_WINDOW_FAILED = 3
} solver_codekd_verification_window_status_t;

/*
 * Score one bounded immutable window on the existing shard pool. Each helper
 * receives only prepared verification contexts and a disjoint score range.
 * Candidate records remain owner-only until the synchronous group is fully
 * quiescent, after which scores are moved into their canonical record slots.
 */
static solver_codekd_verification_window_status_t
solver_codekd_packet_score_verification_window_parallel(
    solver_codekd_search_packet_t* packet) {
    solver_verification_score_slot_t* slots = NULL;
    verify_prepared_score_t* scores = NULL;
    index_shard_helper_task_t tasks[INDEX_SHARD_HELPER_MAX_TASKS];
    solver_verification_task_input_t
        inputs[INDEX_SHARD_HELPER_MAX_TASKS];
    index_shard_helper_run_stats_t run_stats;
    index_shard_helper_run_status_t run_status;
    size_t task_count;
    size_t task_index;
    size_t slot_index = 0U;
    size_t candidate_index;
    double wall_start = 0.0;

    if (!packet || !packet->candidate_records ||
        packet->state != SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER ||
        !packet->verify_plan_complete ||
        packet->verify_plan_first !=
            packet->candidate_window_offset ||
        packet->verify_plan_end <= packet->verify_plan_first ||
        packet->verify_plan_end > packet->candidate_window_count ||
        packet->verify_prepared_count < 2U ||
        packet->verify_prepared_count >
            INDEX_SHARD_HELPER_MAX_TASKS) {
        return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
    }
    slots = calloc(packet->verify_prepared_count, sizeof(*slots));
    scores = calloc(packet->verify_prepared_count, sizeof(*scores));
    if (!slots || !scores) {
        free(slots);
        free(scores);
        return SOLVER_CODEKD_VERIFY_WINDOW_UNAVAILABLE;
    }
    memset(tasks, 0, sizeof(tasks));
    memset(inputs, 0, sizeof(inputs));
    memset(&run_stats, 0, sizeof(run_stats));

    for (candidate_index = packet->verify_plan_first;
         candidate_index < packet->verify_plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];

        if (!record->prepared_verification) {
            continue;
        }
        if (record->plan_action != SOLVER_AB_CANDIDATE_VERIFY ||
            record->verification_score_ready ||
            slot_index >= packet->verify_prepared_count) {
            solver_codekd_packet_destroy_score_array(
                scores, packet->verify_prepared_count);
            free(slots);
            return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
        }
        slots[slot_index].prepared =
            record->prepared_verification;
        slot_index++;
    }
    if (slot_index != packet->verify_prepared_count) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
    }

    /*
     * One prepared context is the smallest safe score unit. The window is
     * capped by SOLVER_VERIFY_QUERY_LOOKAHEAD and the packet memory budget, so
     * the bounded task count cannot flood the shared queue.
     */
    task_count = packet->verify_prepared_count;
    for (task_index = 0U; task_index < task_count; task_index++) {
        unsigned long long work_units =
            verify_prepared_hit_work_units(
                slots[task_index].prepared);

        if (!work_units) {
            work_units = 1U;
        }
        inputs[task_index].slots = &slots[task_index];
        inputs[task_index].slot_first = task_index;
        inputs[task_index].slot_count = 1U;
        tasks[task_index].input = &inputs[task_index];
        tasks[task_index].input_bytes = sizeof(inputs[task_index]);
        tasks[task_index].output = &scores[task_index];
        tasks[task_index].output_bytes = sizeof(scores[task_index]);
        tasks[task_index].work_units = work_units;
    }

    if (packet->detailed) {
        wall_start = monotonic_seconds();
    }
    run_status = index_shard_helper_run(
        &solver_verification_helper_ops,
        tasks,
        task_count,
        &run_stats);
    if (packet->detailed &&
        run_status != INDEX_SHARD_HELPER_UNAVAILABLE) {
        packet->verification_score_wall_seconds +=
            monotonic_seconds() - wall_start;
    }
    if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_UNAVAILABLE;
    }
    if (run_status == INDEX_SHARD_HELPER_STOPPED) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_STOPPED;
    }
    if (run_status == INDEX_SHARD_HELPER_FATAL) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
    }
    if (run_status != INDEX_SHARD_HELPER_OK) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_FAILED;
    }

    slot_index = 0U;
    for (candidate_index = packet->verify_plan_first;
         candidate_index < packet->verify_plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        unsigned long long work_units;

        if (!record->prepared_verification) {
            continue;
        }
        if (slot_index >= packet->verify_prepared_count ||
            !scores[slot_index].complete) {
            solver_codekd_packet_destroy_score_array(
                scores, packet->verify_prepared_count);
            free(slots);
            return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
        }
        record->prepared_score = scores[slot_index];
        memset(&scores[slot_index], 0, sizeof(scores[slot_index]));
        record->verification_score_ready = TRUE;
        work_units = verify_prepared_hit_work_units(
            record->prepared_verification);
        solver_codekd_packet_add_verification_score_work(
            packet, work_units);
        slot_index++;
    }
    if (slot_index != packet->verify_prepared_count) {
        solver_codekd_packet_destroy_score_array(
            scores, packet->verify_prepared_count);
        free(slots);
        return SOLVER_CODEKD_VERIFY_WINDOW_ERROR;
    }
    packet->verification_score_batches_executed++;
    packet->verification_score_contexts_completed +=
        packet->verify_prepared_count;
    solver_codekd_packet_destroy_score_array(
        scores, packet->verify_prepared_count);
    free(slots);
    return SOLVER_CODEKD_VERIFY_WINDOW_SCORED;
}

/*
 * The owner consumes only fully captured native query results and creates a
 * bounded contiguous set of immutable, index-free score contexts. No solver
 * or reducer state is touched. Foreign scoring is attempted only after every
 * context in the selected prefix has been completely prepared.
 */
static int solver_codekd_packet_prepare_verification_owner(
    const solver_codekd_packet_task_input_t* input,
    solver_codekd_search_packet_t* packet) {
    const solver_codekd_verification_snapshot_t* snapshot;
    size_t original_plan_end;
    size_t candidate_index;
    size_t verify_index;
    solver_codekd_verification_window_status_t window_status;

    if (!input || !packet || !packet->candidate_records ||
        packet->state != SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER ||
        !packet->verify_plan_complete ||
        packet->verify_plan_first !=
            packet->candidate_window_offset ||
        packet->verify_plan_end <= packet->verify_plan_first ||
        packet->verify_plan_end > packet->candidate_window_count ||
        packet->verify_prepared_count ||
        packet->verify_prepared_retained_bytes ||
        packet->verify_prepared_transient_bytes) {
        return -1;
    }
    snapshot = &input->verification;
    if (!snapshot->enabled || !snapshot->field ||
        verify_datalog_enabled()) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    original_plan_end = packet->verify_plan_end;
    verify_index = original_plan_end;
    for (candidate_index = packet->verify_plan_first;
         candidate_index < original_plan_end;
         candidate_index++) {
        if (packet->candidate_records[candidate_index].plan_action ==
                SOLVER_AB_CANDIDATE_VERIFY) {
            verify_index = candidate_index;
            break;
        }
    }
    if (verify_index == original_plan_end) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    if (verify_index > packet->verify_plan_first) {
        packet->verify_plan_end = verify_index;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    if (!index_shard_helper_available_workers()) {
        original_plan_end = verify_index + 1U;
        packet->verify_plan_end = original_plan_end;
    }

    for (candidate_index = verify_index;
         candidate_index < original_plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        verify_prepared_hit_t* prepared = NULL;
        MatchObj match;
        size_t query_bytes;
        size_t remaining_query_bytes;
        size_t prepared_bytes;
        size_t score_bytes;
        size_t peak_bytes;
        size_t retained_bytes;
        size_t transient_bytes;
        size_t retained_total;
        size_t transient_max;
        double verify_pix2;
        double logaccept;

        if (record->plan_action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        if (index_shard_worker_stop_requested()) {
            solver_codekd_packet_discard_prepared_verification(packet);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return 2;
        }
        if (!record->candidate_prepared ||
            !record->verify_query ||
            !record->verify_query_captured) {
            record->verify_delivery_fallback = TRUE;
            goto close_prefix;
        }
        query_bytes = verify_index_query_bytes(record->verify_query);
        if (query_bytes == SIZE_MAX ||
            query_bytes > packet->verify_query_bytes ||
            solver_codekd_packet_build_verification_match(
                input, packet, record, candidate_index,
                &match, &verify_pix2, &logaccept) ||
            verify_prepare_hit_from_query(
                packet->starkd,
                &record->verify_query,
                snapshot->index_cutnside,
                &match,
                NULL,
                snapshot->field,
                verify_pix2,
                snapshot->distractor_ratio,
                input->field_maxx,
                input->field_maxy,
                snapshot->logratio_bail_threshold,
                logaccept,
                snapshot->logratio_stoplooking,
                snapshot->distance_from_quad_bonus,
                FALSE,
                &prepared)) {
            verify_destroy_prepared_hit(prepared);
            record->verify_delivery_fallback = TRUE;
            goto close_prefix;
        }
        record->verify_query_captured = FALSE;
        remaining_query_bytes =
            packet->verify_query_bytes - query_bytes;
        prepared_bytes = verify_prepared_hit_bytes(prepared);
        score_bytes = verify_prepared_score_bytes(prepared);
        peak_bytes = verify_prepared_hit_peak_bytes(prepared);
        if (prepared_bytes == SIZE_MAX ||
            score_bytes == SIZE_MAX ||
            peak_bytes == SIZE_MAX ||
            prepared_bytes > SIZE_MAX - score_bytes ||
            peak_bytes < prepared_bytes + score_bytes ||
            remaining_query_bytes > packet->verify_query_budget ||
            packet->verify_prepared_retained_bytes >
                SIZE_MAX - prepared_bytes - score_bytes) {
            verify_destroy_prepared_hit(prepared);
            record->verify_delivery_fallback = TRUE;
            packet->verify_query_bytes = remaining_query_bytes;
            goto close_prefix;
        }
        retained_bytes = prepared_bytes + score_bytes;
        transient_bytes = peak_bytes - retained_bytes;
        retained_total =
            packet->verify_prepared_retained_bytes +
            retained_bytes;
        transient_max = MAX(
            packet->verify_prepared_transient_bytes,
            transient_bytes);
        if (retained_total >
                packet->verify_query_budget -
                    remaining_query_bytes ||
            transient_max >
                packet->verify_query_budget -
                    remaining_query_bytes -
                    retained_total) {
            verify_destroy_prepared_hit(prepared);
            record->verify_delivery_fallback = TRUE;
            packet->verify_query_bytes = remaining_query_bytes;
            goto close_prefix;
        }
        packet->verify_query_bytes = remaining_query_bytes;
        record->prepared_verification = prepared;
        record->prepared_verify_pix2 = verify_pix2;
        record->prepared_logaccept = logaccept;
        record->prepared_distractor_ratio =
            snapshot->distractor_ratio;
        record->prepared_logratio_bail_threshold =
            snapshot->logratio_bail_threshold;
        record->prepared_logratio_stoplooking =
            snapshot->logratio_stoplooking;
        record->prepared_field_maxx = input->field_maxx;
        record->prepared_field_maxy = input->field_maxy;
        record->prepared_distance_from_quad_bonus =
            snapshot->distance_from_quad_bonus;
        packet->verify_prepared_retained_bytes =
            retained_total;
        packet->verify_prepared_transient_bytes =
            transient_max;
        packet->verify_prepared_count++;
        continue;

close_prefix:
        if (!packet->verify_prepared_count) {
            packet->verify_plan_end = candidate_index + 1U;
            packet->verification_score_fallback_batches++;
            solver_codekd_packet_discard_prepared_verification(packet);
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            return 0;
        }
        packet->verify_plan_end = candidate_index;
        packet->verify_topology_end = candidate_index;
        break;
    }

    if (!packet->verify_prepared_count ||
        packet->verify_plan_end <= packet->verify_plan_first) {
        packet->verification_score_fallback_batches++;
        solver_codekd_packet_discard_prepared_verification(packet);
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    packet->verification_score_batches_prepared++;
    packet->verification_score_contexts_prepared +=
        packet->verify_prepared_count;
    if (packet->verify_prepared_count == 1U) {
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY;
        return 1;
    }

    window_status =
        solver_codekd_packet_score_verification_window_parallel(packet);
    if (window_status == SOLVER_CODEKD_VERIFY_WINDOW_SCORED) {
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    if (window_status == SOLVER_CODEKD_VERIFY_WINDOW_UNAVAILABLE) {
        packet->state =
            SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY;
        return 1;
    }
    if (window_status == SOLVER_CODEKD_VERIFY_WINDOW_STOPPED) {
        packet->verification_score_stopped_batches++;
        solver_codekd_packet_discard_prepared_verification(packet);
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return 2;
    }
    if (window_status == SOLVER_CODEKD_VERIFY_WINDOW_FAILED) {
        packet->verification_score_fallback_batches++;
        solver_codekd_packet_discard_prepared_verification(packet);
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return 0;
    }
    solver_codekd_packet_discard_prepared_verification(packet);
    packet->state = SOLVER_CODEKD_PACKET_FAILED;
    return -1;
}

static index_shard_helper_task_status_t
solver_codekd_packet_score_verification_ready(
    solver_codekd_search_packet_t* packet) {
    index_shard_helper_task_status_t status =
        INDEX_SHARD_HELPER_TASK_OK;
    size_t candidate_index;
    size_t completed = 0U;
    double wall_start = 0.0;

    if (!packet || !packet->candidate_records ||
        packet->state !=
            SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY ||
        !packet->verify_plan_complete ||
        packet->verify_plan_first !=
            packet->candidate_window_offset ||
        packet->verify_plan_first == SIZE_MAX ||
        packet->verify_plan_end <= packet->verify_plan_first ||
        packet->verify_plan_end > packet->candidate_window_count ||
        !packet->verify_prepared_count ||
        !packet->verify_prepared_retained_bytes) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    if (packet->detailed) {
        wall_start = monotonic_seconds();
    }
    packet->verification_score_batches_executed++;
    for (candidate_index = packet->verify_plan_first;
         candidate_index < packet->verify_plan_end;
         candidate_index++) {
        solver_candidate_delivery_record_t* record =
            &packet->candidate_records[candidate_index];
        unsigned long long work_units;

        if (!record->prepared_verification) {
            continue;
        }
        if (record->plan_action != SOLVER_AB_CANDIDATE_VERIFY ||
            record->verification_score_ready) {
            status = INDEX_SHARD_HELPER_TASK_ERROR;
            goto done;
        }
        if (index_shard_worker_stop_requested()) {
            packet->verification_score_stopped_batches++;
            solver_codekd_packet_discard_prepared_verification(packet);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            status = INDEX_SHARD_HELPER_TASK_STOPPED;
            goto done;
        }
        if (verify_score_prepared_hit(
                record->prepared_verification,
                &record->prepared_score)) {
            packet->verification_score_fallback_batches++;
            solver_codekd_packet_discard_prepared_verification(packet);
            packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
            goto done;
        }
        record->verification_score_ready = TRUE;
        completed++;
        work_units = verify_prepared_hit_work_units(
            record->prepared_verification);
        solver_codekd_packet_add_verification_score_work(
            packet, work_units);
    }
    if (completed != packet->verify_prepared_count) {
        status = INDEX_SHARD_HELPER_TASK_ERROR;
        goto done;
    }
    packet->verification_score_contexts_completed += completed;
    packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;

done:
    if (packet->detailed) {
        packet->verification_score_wall_seconds +=
            monotonic_seconds() - wall_start;
    }
    return status;
}

/*
 * Execute only a fully populated, contiguous descriptor slice. Native CodeKD
 * remains the authoritative query implementation. Page eviction between I/O
 * completion and execution is harmless; it affects performance, not results.
 */
static index_shard_helper_task_status_t
solver_codekd_packet_execute_ready(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size,
    anbool* more_work) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    kdtree_qres_t* query_result = NULL;
    size_t descriptor_index;
    anbool query_failed = FALSE;

    if (more_work) {
        *more_work = FALSE;
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY) {
        return solver_codekd_packet_decode_quad_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY) {
        return solver_codekd_packet_copy_star_ready(
            input, packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY) {
        return solver_codekd_packet_query_verification_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY) {
        return solver_codekd_packet_capture_sweep_ready(
            packet, more_work);
    }
    if (packet && more_work &&
        packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY) {
        return solver_codekd_packet_score_verification_ready(packet);
    }
    if (packet && more_work &&
        packet->state == SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY) {
        if (!input || input_size != sizeof(*input) ||
            output_size != sizeof(*packet) ||
            !packet->verify_plan_complete ||
            packet->verify_plan_first !=
                packet->candidate_window_offset ||
            packet->verify_plan_end <= packet->verify_plan_first ||
            packet->verify_plan_end >
                packet->candidate_window_count) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        return INDEX_SHARD_HELPER_TASK_OK;
    }
    if (!input || input_size != sizeof(*input) ||
        !packet || output_size != sizeof(*packet) ||
        !input->tree || packet->tree != input->tree ||
        !packet->descriptors || !packet->slots ||
        !packet->inds || !packet->sdists ||
        !packet->plan_complete ||
        packet->plan_first >= packet->plan_end ||
        packet->plan_end > packet->count ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_COMPUTE_READY) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    packet->state = SOLVER_CODEKD_PACKET_EXECUTING;
    for (descriptor_index = packet->plan_first;
         descriptor_index < packet->plan_end;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[descriptor_index];
        solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];
        kdtree_qres_t* previous_result = query_result;
        double search_wall_start = 0.0;

        if (index_shard_worker_stop_requested()) {
            kdtree_free_query(query_result);
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        if (packet->detailed) {
            search_wall_start = monotonic_seconds();
        }
        query_result = solver_codekd_rangesearch(
            input->tree,
            previous_result,
            descriptor->code,
            descriptor->tol2,
            SOLVER_CODEKD_SEARCH_OPTIONS);
        slot->search_errno = errno;
        if (packet->detailed) {
            slot->search_wall_seconds =
                monotonic_seconds() - search_wall_start;
        }
        if (!query_result) {
            kdtree_free_query(previous_result);
            query_result = NULL;
            slot->state = SOLVER_CODEKD_RESULT_QUERY_FAILED;
            query_failed = TRUE;
            break;
        }
        if (query_result->nres &&
            (!query_result->inds || !query_result->sdists)) {
            slot->search_errno = EIO;
            slot->state = SOLVER_CODEKD_RESULT_QUERY_FAILED;
            query_failed = TRUE;
            break;
        }
        if (packet->hit_count > packet->hit_capacity ||
            (size_t)query_result->nres >
                packet->hit_capacity - packet->hit_count) {
            slot->state = SOLVER_CODEKD_RESULT_OWNER_REPLAY;
            continue;
        }

        slot->hit_first = packet->hit_count;
        slot->hit_count = query_result->nres;
        if (slot->hit_count) {
            memcpy(
                packet->inds + slot->hit_first,
                query_result->inds,
                (size_t)slot->hit_count * sizeof(*packet->inds));
            memcpy(
                packet->sdists + slot->hit_first,
                query_result->sdists,
                (size_t)slot->hit_count * sizeof(*packet->sdists));
        }
        packet->hit_count += slot->hit_count;
        slot->state = SOLVER_CODEKD_RESULT_READY;
    }

    kdtree_free_query(query_result);
    packet->plan_complete = FALSE;
    packet->plan_first = 0U;
    packet->plan_end = 0U;
    packet->plan_range_count = 0U;
    packet->plan_logical_bytes = 0U;
    packet->delivery_source = NULL;
    if (query_failed) {
        packet->next_descriptor = packet->count;
        packet->pending_descriptor_plan = FALSE;
        packet->pending_descriptor_raw_ranges = 0U;
        packet->pending_descriptor_logical_bytes = 0U;
        packet->state = SOLVER_CODEKD_PACKET_RESULTS_READY;
        if (packet->hit_count) {
            if (solver_codekd_packet_finish_codekd(packet)) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_HELPER_TASK_ERROR;
            }
            if (packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY &&
                more_work) {
                *more_work = TRUE;
            }
        }
    } else if (packet->next_descriptor < packet->count) {
        packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
        if (more_work) {
            *more_work = TRUE;
        }
    } else {
        packet->state = SOLVER_CODEKD_PACKET_DESCRIPTORS_READY;
        if (more_work) {
            *more_work = TRUE;
        }
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

static index_shard_staged_prepare_status_t
solver_codekd_packet_staged_prepare(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;
    index_shard_helper_task_status_t descriptor_status;
    int plan_status;

    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_PREPARE_ERROR;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_ALLOCATED) {
        descriptor_status = solver_codekd_packet_generate_descriptors(
            input_bytes,
            input_size,
            output_bytes,
            output_size);
        if (descriptor_status == INDEX_SHARD_HELPER_TASK_STOPPED) {
            return INDEX_SHARD_STAGED_PREPARE_STOPPED;
        }
        if (descriptor_status != INDEX_SHARD_HELPER_TASK_OK) {
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
    }
    if (packet->state != SOLVER_CODEKD_PACKET_DESCRIPTORS_READY) {
        if (packet->state ==
                SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE) {
            return packet->plan_complete &&
                    packet->plan_first < packet->plan_end &&
                    packet->plan_end <= packet->next_descriptor &&
                    packet->next_descriptor <= packet->count &&
                    packet->plan_range_count &&
                    packet->plan_range_count <=
                        packet->page_workspace->sealed_range_capacity
                ? INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY
                : INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER) {
            return INDEX_SHARD_STAGED_PREPARE_OWNER_READY;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY) {
            return INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY;
        }
        if (packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY) {
            if (!packet->verify_plan_complete ||
                packet->verify_plan_first !=
                    packet->candidate_window_offset ||
                packet->verify_plan_end <=
                    packet->verify_plan_first ||
                packet->verify_plan_end !=
                    packet->verify_topology_end ||
                packet->verify_topology_end >
                    packet->candidate_window_count) {
                logerr("[solver] invalid verification query packet "
                       "state=%i first=%zu end=%zu topology_end=%zu "
                       "window_offset=%zu window_count=%zu\n",
                       (int)packet->state,
                       packet->verify_plan_first,
                       packet->verify_plan_end,
                       packet->verify_topology_end,
                       packet->candidate_window_offset,
                       packet->candidate_window_count);
                return INDEX_SHARD_STAGED_PREPARE_ERROR;
            }
            return INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY;
        }
        if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY ||
            packet->state == SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY ||
            packet->state == SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY ||
            packet->state ==
                SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY) {
            return INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
        }
        if (packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY) {
            logerr("[solver] unexpected packet state during prepare "
                   "state=%i sequence=%llu\n",
                   (int)packet->state,
                   packet->sequence);
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (solver_codekd_packet_finish_codekd(packet)) {
            logerr("[solver] packet finish failed during prepare "
                   "state=%i sequence=%llu retire=%zu/%zu "
                   "verify=%zu/%zu prepared=%zu\n",
                   (int)packet->state,
                   packet->sequence,
                   packet->retire_descriptor,
                   packet->count,
                   packet->verify_plan_first,
                   packet->verify_plan_end,
                   packet->verify_prepared_count);
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_PREPARE_RESULTS_READY
            : INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
    }
    if (packet->next_descriptor < packet->count) {
        plan_status = solver_codekd_search_packet_prepare_next_plan(
            packet,
            solver_codekd_packet_plan_cancelled,
            NULL);
        if (plan_status < 0) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        if (plan_status == 2) {
            return INDEX_SHARD_STAGED_PREPARE_STOPPED;
        }
        if (plan_status > 0) {
            return packet->state ==
                        SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE &&
                    packet->plan_complete &&
                    packet->plan_first < packet->plan_end &&
                    packet->plan_end <= packet->next_descriptor &&
                    packet->next_descriptor <= packet->count &&
                    packet->plan_range_count &&
                    packet->plan_range_count <=
                        packet->page_workspace->sealed_range_capacity
                ? INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY
                : INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
    }
    if (packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet->next_descriptor == packet->count) {
        if (solver_codekd_packet_finish_codekd(packet)) {
            logerr("[solver] packet finish failed after planning "
                   "state=%i sequence=%llu retire=%zu/%zu "
                   "verify=%zu/%zu prepared=%zu\n",
                   (int)packet->state,
                   packet->sequence,
                   packet->retire_descriptor,
                   packet->count,
                   packet->verify_plan_first,
                   packet->verify_plan_end,
                   packet->verify_prepared_count);
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_PREPARE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_PREPARE_RESULTS_READY
            : INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY;
    }
    logerr("[solver] invalid terminal packet prepare state=%i "
           "sequence=%llu next=%zu count=%zu\n",
           (int)packet->state,
           packet->sequence,
           packet->next_descriptor,
           packet->count);
    packet->state = SOLVER_CODEKD_PACKET_FAILED;
    return INDEX_SHARD_STAGED_PREPARE_ERROR;
}

static index_shard_staged_submit_status_t
solver_codekd_packet_staged_submit(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size,
    unsigned long long* completion_id_out) {
    solver_codekd_search_packet_t* packet = output_bytes;
    int status;

    (void)input_bytes;
    (void)input_size;
    if (completion_id_out) {
        *completion_id_out = 0ULL;
    }
    if (!packet || output_size != sizeof(*packet) ||
        !completion_id_out) {
        return INDEX_SHARD_STAGED_SUBMIT_ERROR;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY ||
        packet->state == SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY) {
        status = solver_codekd_packet_submit_candidate_pages(packet);
    } else if (packet->state ==
               SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY) {
        status =
            solver_codekd_packet_submit_verification_pages(packet);
    } else if (packet->state ==
               SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY) {
        status = solver_codekd_packet_submit_sweep_pages(packet);
    } else {
        status = solver_codekd_packet_submit_pages(packet);
    }
    if (status == 2 &&
        (packet->state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY ||
         packet->state == SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY ||
         packet->state ==
             SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY ||
         packet->state ==
             SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY ||
         packet->state == SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY)) {
        return INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY;
    }
    if (status > 0) {
        *completion_id_out =
            fitsbin_payload_io_ticket_completion_id(
                packet->delivery_ticket);
        if (!*completion_id_out) {
            (void)solver_codekd_search_packet_release_ticket(packet);
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_SUBMIT_ERROR;
        }
        return INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED;
    }
    if (!status) {
        return INDEX_SHARD_STAGED_SUBMIT_RETRY;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_STOPPED) {
        return INDEX_SHARD_STAGED_SUBMIT_STOPPED;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_DESCRIPTORS_READY ||
        packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY) {
        return INDEX_SHARD_STAGED_SUBMIT_OWNER_READY;
    }
    return INDEX_SHARD_STAGED_SUBMIT_ERROR;
}

static index_shard_staged_io_status_t
solver_codekd_packet_staged_poll(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;
    int status;

    (void)input_bytes;
    (void)input_size;
    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_IO_ERROR;
    }
    status = solver_codekd_packet_collect_pages(packet);
    if (!status) {
        return INDEX_SHARD_STAGED_IO_PENDING;
    }
    if (status == 1) {
        return INDEX_SHARD_STAGED_IO_READY;
    }
    if (status == 2 || status == 3) {
        return INDEX_SHARD_STAGED_IO_FAILED;
    }
    return packet->state == SOLVER_CODEKD_PACKET_STOPPED
        ? INDEX_SHARD_STAGED_IO_CANCELLED
        : INDEX_SHARD_STAGED_IO_ERROR;
}

static int solver_codekd_packet_staged_cancel(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    solver_codekd_search_packet_t* packet = output_bytes;

    (void)input_bytes;
    (void)input_size;
    if (!packet || output_size != sizeof(*packet)) {
        return -1;
    }
    return solver_codekd_packet_cancel_pages(packet) < 0
        ? -1
        : 0;
}

static index_shard_staged_execute_status_t
solver_codekd_packet_staged_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    anbool more_work = FALSE;
    index_shard_helper_task_status_t status =
        solver_codekd_packet_execute_ready(
            input_bytes,
            input_size,
            output_bytes,
            output_size,
            &more_work);

    if (status == INDEX_SHARD_HELPER_TASK_STOPPED) {
        return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
    }
    if (status != INDEX_SHARD_HELPER_TASK_OK) {
        return INDEX_SHARD_STAGED_EXECUTE_ERROR;
    }
    return more_work
        ? INDEX_SHARD_STAGED_EXECUTE_MORE
        : INDEX_SHARD_STAGED_EXECUTE_OK;
}

static index_shard_staged_execute_status_t
solver_codekd_packet_staged_owner(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_codekd_packet_task_input_t* input = input_bytes;
    solver_codekd_search_packet_t* packet = output_bytes;
    int prepare_status;

    if (!packet || output_size != sizeof(*packet)) {
        return INDEX_SHARD_STAGED_EXECUTE_ERROR;
    }
    if (index_shard_worker_stop_requested() ||
        packet->state == SOLVER_CODEKD_PACKET_STOPPED) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_DESCRIPTORS_READY) {
        return INDEX_SHARD_STAGED_EXECUTE_MORE;
    }
    if (packet->state ==
            SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER) {
        if (!input || input_size != sizeof(*input)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        prepare_status =
            solver_codekd_packet_prepare_verification_owner(
                input, packet);
        if (prepare_status == 2) {
            return INDEX_SHARD_STAGED_EXECUTE_STOPPED;
        }
        if (prepare_status < 0) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        return prepare_status
            ? INDEX_SHARD_STAGED_EXECUTE_MORE
            : INDEX_SHARD_STAGED_EXECUTE_OK;
    }
    if (packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY) {
        if (solver_codekd_packet_finish_codekd(packet)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_EXECUTE_ERROR;
        }
        return packet->state == SOLVER_CODEKD_PACKET_RESULTS_READY
            ? INDEX_SHARD_STAGED_EXECUTE_OK
            : INDEX_SHARD_STAGED_EXECUTE_MORE;
    }
    return INDEX_SHARD_STAGED_EXECUTE_ERROR;
}

static const index_shard_staged_ops_t solver_codekd_packet_staged_ops = {
    "codekd-packet",
    solver_codekd_packet_staged_prepare,
    solver_codekd_packet_staged_submit,
    solver_codekd_packet_staged_poll,
    solver_codekd_packet_staged_cancel,
    solver_codekd_packet_staged_execute,
    solver_codekd_packet_staged_owner
};

/* Exact native execution used whenever packet admission is unavailable. */
static int solver_codekd_descriptor_execute_owner(
    const solver_ab_descriptor_output_t* output,
    solver_t* solver,
    int dimquads,
    kdtree_qres_t** query_result,
    unsigned long long* reduced) {
    size_t descriptor_index;

    if (!output || !solver || !query_result || !reduced) {
        return -1;
    }
    solver->profile.max_batch_hypotheses = MAX(
        solver->profile.max_batch_hypotheses,
        output->descriptor_count);
    for (descriptor_index = 0U;
         descriptor_index < output->descriptor_count;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &output->descriptors[descriptor_index];

        if (solver_poll_worker_stop(solver)) {
            return 1;
        }
        solver->rel_field_noise2 = descriptor->rel_field_noise2;
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
        solver->rel_field_noise2 = output->final_rel_field_noise2;
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

static void solver_codekd_packet_begin_result_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    const solver_codekd_result_slot_t* slot,
    const kdtree_qres_t* result,
    int dimquad) {
    if (!solver || !descriptor || !slot || !result ||
        dimquad <= 0 || dimquad > DQMAX) {
        return;
    }
    solver_begin_hypothesis_owner(
        descriptor->stars,
        descriptor->code,
        dimquad,
        solver,
        descriptor->current_parity);
    solver->profile.codekd_calls++;
    solver->profile.hypotheses_executed++;
    if (solver->profile.detailed) {
        solver->profile.codekd_wall_seconds +=
            slot->search_wall_seconds;
        solver->profile.kd_result_order_hash =
            solver_order_hash_mix(
                solver->profile.kd_result_order_hash,
                solver_kd_result_order_digest(result));
    }
    solver->profile.codekd_hits +=
        (unsigned long long)result->nres;
    if (result->nres) {
        solver->profile.resolve_calls++;
    }
}

static void solver_codekd_packet_resolve_result_range_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    kdtree_qres_t* result,
    int dimquad,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    double pixvals[DQMAX * 2];
    double resolve_wall_start = 0.0;
    int j;

    if (!solver || !descriptor || !result ||
        dimquad <= 0 || dimquad > DQMAX ||
        candidate_first < 0 ||
        candidate_end <= candidate_first ||
        candidate_end > result->nres) {
        return;
    }
    for (j = 0; j < dimquad; j++) {
        setx(
            pixvals,
            j,
            field_getx(solver, descriptor->stars[j]));
        sety(
            pixvals,
            j,
            field_gety(solver, descriptor->stars[j]));
    }
    if (solver->profile.detailed) {
        resolve_wall_start = monotonic_seconds();
    }
    resolve_matches_native_range(
        result,
        pixvals,
        descriptor->stars,
        dimquad,
        solver->numtries,
        solver,
        descriptor->current_parity,
        candidate_first,
        candidate_end,
        prepared,
        prepared_first,
        prepared_quad_count,
        prepared_star_count);
    if (solver->profile.detailed) {
        solver->profile.resolve_wall_seconds +=
            monotonic_seconds() - resolve_wall_start;
    }
}

typedef struct solver_codekd_packet_retire_context {
    solver_t* solver;
    int dimquads;
    kdtree_qres_t** query_result;
    unsigned long long* reduced;
    size_t next_task_index;
    unsigned long long next_sequence;
} solver_codekd_packet_retire_context_t;

static index_shard_staged_retire_status_t
solver_codekd_packet_retire_nonwindow_descriptor(
    solver_codekd_search_packet_t* packet,
    solver_codekd_packet_retire_context_t* context,
    const solver_ab_descriptor_t* descriptor,
    const solver_codekd_result_slot_t* slot) {
    if (!packet || !context || !context->solver ||
        !context->query_result || !context->reduced ||
        !descriptor || !slot ||
        packet->retire_descriptor_started ||
        (slot->state == SOLVER_CODEKD_RESULT_READY &&
         slot->hit_count)) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    if (slot->state != SOLVER_CODEKD_RESULT_QUERY_FAILED &&
        solver_poll_worker_stop(context->solver)) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_RETIRE_STOPPED;
    }
    context->solver->rel_field_noise2 =
        descriptor->rel_field_noise2;
    if (solver_ab_checked_counter_delta(
            context->solver,
            descriptor->numtries_delta,
            descriptor->cxdx_delta,
            descriptor->meanx_delta)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    if (slot->state == SOLVER_CODEKD_RESULT_READY) {
        kdtree_qres_t result_view;

        memset(&result_view, 0, sizeof(result_view));
        solver_codekd_packet_begin_result_owner(
            context->solver,
            descriptor,
            slot,
            &result_view,
            context->dimquads);
        context->solver->profile.hypotheses_reduced++;
    } else if (slot->state ==
               SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
        solver_execute_hypothesis_owner(
            descriptor->stars,
            descriptor->code,
            context->dimquads,
            context->solver,
            descriptor->current_parity,
            descriptor->tol2,
            context->query_result);
    } else if (slot->state ==
               SOLVER_CODEKD_RESULT_QUERY_FAILED) {
        solver_retire_codekd_failure_owner(
            descriptor->stars,
            descriptor->code,
            context->dimquads,
            context->solver,
            descriptor->current_parity,
            slot->search_errno,
            slot->search_wall_seconds);
    } else {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    (*context->reduced)++;
    packet->retire_descriptor++;
    if (context->solver->profile.execution_failed) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    if (context->solver->quit_now) {
        packet->state = SOLVER_CODEKD_PACKET_STOPPED;
        return INDEX_SHARD_STAGED_RETIRE_STOPPED;
    }
    return INDEX_SHARD_STAGED_RETIRE_OK;
}

static index_shard_staged_retire_status_t
solver_codekd_search_packet_retire_window(
    solver_codekd_search_packet_t* packet,
    solver_codekd_packet_retire_context_t* context) {
    size_t window_end;

    if (!packet || !context || !context->solver ||
        packet->state != SOLVER_CODEKD_PACKET_RETIRING ||
        packet->candidate_count != packet->hit_count ||
        !packet->candidate_count ||
        !packet->candidate_records ||
        !packet->candidate_capacity ||
        packet->candidate_cursor >= packet->candidate_count ||
        !packet->candidate_window_count ||
        packet->candidate_window_first >
            packet->candidate_cursor ||
        packet->candidate_window_offset !=
            packet->candidate_cursor -
                packet->candidate_window_first ||
        packet->candidate_window_offset >=
            packet->candidate_window_count ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            SOLVER_CANDIDATE_DELIVERY_LIMIT ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_quad_ready_count >
            packet->candidate_window_count ||
        packet->candidate_star_ready_count >
            packet->candidate_quad_ready_count ||
        (packet->verify_plan_complete &&
         (packet->verification_delivery_disabled ||
          packet->verify_plan_first !=
              packet->candidate_window_offset ||
          packet->verify_plan_end <=
              packet->verify_plan_first ||
          packet->verify_plan_end >
              packet->verify_topology_end ||
          packet->verify_topology_end >
              packet->candidate_window_count)) ||
        (!packet->verification_delivery_disabled &&
         packet->candidate_verify_query_count &&
         !packet->verify_plan_complete) ||
        packet->retire_descriptor >= packet->count) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    window_end =
        packet->candidate_window_first +
            packet->candidate_window_count;

    while (packet->retire_descriptor < packet->count) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[
                packet->retire_descriptor];
        const solver_codekd_result_slot_t* slot =
            &packet->slots[packet->retire_descriptor];
        index_shard_staged_retire_status_t status;

        if (slot->state == SOLVER_CODEKD_RESULT_READY &&
            slot->hit_count) {
            kdtree_qres_t result_view;
            size_t global_first;
            size_t local_end;
            size_t prepared_quad_count;
            size_t prepared_star_count;
            size_t range_count;
            size_t processed_count;
            int nummatches_before;

            if (!packet->candidate_window_count &&
                !packet->candidate_quad_delivery_disabled) {
                if (solver_codekd_packet_begin_candidate_window(
                        packet)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first ||
                packet->retire_hit_offset >
                    (size_t)slot->hit_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            memset(&result_view, 0, sizeof(result_view));
            result_view.nres = slot->hit_count;
            result_view.capacity = slot->hit_count;
            result_view.inds = packet->inds + slot->hit_first;
            result_view.sdists = packet->sdists + slot->hit_first;

            if (!packet->retire_descriptor_started) {
                if (solver_poll_worker_stop(context->solver)) {
                    packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                    return INDEX_SHARD_STAGED_RETIRE_STOPPED;
                }
                context->solver->rel_field_noise2 =
                    descriptor->rel_field_noise2;
                if (solver_ab_checked_counter_delta(
                        context->solver,
                        descriptor->numtries_delta,
                        descriptor->cxdx_delta,
                        descriptor->meanx_delta)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                solver_codekd_packet_begin_result_owner(
                    context->solver,
                    descriptor,
                    slot,
                    &result_view,
                    context->dimquads);
                packet->retire_descriptor_started = TRUE;
            }

            global_first =
                slot->hit_first + packet->retire_hit_offset;
            if (global_first != packet->candidate_cursor) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->candidate_quad_delivery_disabled) {
                local_end = (size_t)slot->hit_count;
                range_count = local_end -
                    packet->retire_hit_offset;
                nummatches_before =
                    context->solver->nummatches;
                solver_codekd_packet_resolve_result_range_owner(
                    context->solver,
                    descriptor,
                    &result_view,
                    context->dimquads,
                    (int)packet->retire_hit_offset,
                    (int)local_end,
                    NULL,
                    0U,
                    0U,
                    0U);
                if (context->solver->nummatches <
                    nummatches_before) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                processed_count = (size_t)(
                    context->solver->nummatches -
                    nummatches_before);
                if (processed_count > range_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                packet->candidate_retired_rows +=
                    processed_count;
                packet->candidate_native_rows +=
                    processed_count;
                packet->retire_hit_offset +=
                    processed_count;
                packet->candidate_cursor +=
                    processed_count;
                solver_codekd_packet_clear_candidate_window(
                    packet);

                if (context->solver->profile.execution_failed) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                if (context->solver->quit_now ||
                    solver_poll_worker_stop(context->solver)) {
                    context->solver->profile.hypotheses_reduced++;
                    (*context->reduced)++;
                    packet->retire_descriptor_started = FALSE;
                    packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                    return INDEX_SHARD_STAGED_RETIRE_STOPPED;
                }
                if (processed_count != range_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                context->solver->profile.hypotheses_reduced++;
                (*context->reduced)++;
                packet->retire_descriptor_started = FALSE;
                packet->retire_hit_offset = 0U;
                packet->retire_descriptor++;
                continue;
            }
            if (packet->candidate_cursor >= window_end ||
                packet->candidate_window_offset >=
                    packet->candidate_window_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            range_count = MIN(
                (size_t)slot->hit_count -
                    packet->retire_hit_offset,
                window_end - packet->candidate_cursor);
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled) {
                size_t verify_remaining =
                    packet->verify_plan_end -
                        packet->candidate_window_offset;

                range_count = MIN(
                    range_count, verify_remaining);
            }
            if (!range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            local_end =
                packet->retire_hit_offset +
                    range_count;
            prepared_quad_count =
                packet->candidate_quad_ready_count >
                    packet->candidate_window_offset
                ? MIN(
                    range_count,
                    packet->candidate_quad_ready_count -
                        packet->candidate_window_offset)
                : 0U;
            prepared_star_count =
                packet->candidate_star_ready_count >
                    packet->candidate_window_offset
                ? MIN(
                    range_count,
                    packet->candidate_star_ready_count -
                        packet->candidate_window_offset)
                : 0U;
            nummatches_before = context->solver->nummatches;
            solver_codekd_packet_resolve_result_range_owner(
                context->solver,
                descriptor,
                &result_view,
                context->dimquads,
                (int)packet->retire_hit_offset,
                (int)local_end,
                packet->candidate_records +
                    packet->candidate_window_offset,
                packet->retire_hit_offset,
                prepared_quad_count,
                prepared_star_count);
            if (context->solver->nummatches < nummatches_before) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            processed_count = (size_t)(
                context->solver->nummatches - nummatches_before);
            if (processed_count > range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            packet->candidate_retired_rows += processed_count;
            packet->candidate_native_rows +=
                processed_count > prepared_star_count
                ? processed_count -
                    prepared_star_count
                : 0U;
            packet->retire_hit_offset += processed_count;
            packet->candidate_cursor += processed_count;
            packet->candidate_window_offset += processed_count;

            if (context->solver->profile.execution_failed) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (context->solver->quit_now ||
                solver_poll_worker_stop(context->solver)) {
                context->solver->profile.hypotheses_reduced++;
                (*context->reduced)++;
                packet->retire_descriptor_started = FALSE;
                packet->state = SOLVER_CODEKD_PACKET_STOPPED;
                return INDEX_SHARD_STAGED_RETIRE_STOPPED;
            }
            if (processed_count != range_count) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled &&
                packet->candidate_window_offset ==
                    packet->verify_plan_end &&
                solver_codekd_packet_retired_verification_complete(
                    packet)) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            if (packet->verify_plan_complete &&
                !packet->verification_delivery_disabled &&
                packet->candidate_window_offset ==
                    packet->verify_plan_end &&
                packet->candidate_window_offset <
                    packet->candidate_window_count) {
                size_t topology_end =
                    packet->verify_topology_end;

                if (packet->retire_hit_offset ==
                    (size_t)slot->hit_count) {
                    context->solver->profile.hypotheses_reduced++;
                    (*context->reduced)++;
                    packet->retire_descriptor_started = FALSE;
                    packet->retire_hit_offset = 0U;
                    packet->retire_descriptor++;
                }
                solver_codekd_packet_reset_verify_plan(packet);
                if (packet->candidate_window_offset <
                    topology_end) {
                    packet->verify_plan_first =
                        packet->candidate_window_offset;
                    packet->verify_plan_end = topology_end;
                    packet->verify_topology_end = topology_end;
                    packet->verify_plan_complete = TRUE;
                    packet->state =
                        SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY;
                } else {
                    packet->state =
                        SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            if (packet->candidate_window_offset ==
                packet->candidate_window_count) {
                solver_codekd_packet_clear_candidate_window(
                    packet);
            }
            if (packet->retire_hit_offset <
                (size_t)slot->hit_count) {
                if (packet->candidate_window_count) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                if (solver_codekd_packet_begin_candidate_window(
                        packet)) {
                    packet->state = SOLVER_CODEKD_PACKET_FAILED;
                    return INDEX_SHARD_STAGED_RETIRE_ERROR;
                }
                return INDEX_SHARD_STAGED_RETIRE_MORE;
            }
            context->solver->profile.hypotheses_reduced++;
            (*context->reduced)++;
            packet->retire_descriptor_started = FALSE;
            packet->retire_hit_offset = 0U;
            packet->retire_descriptor++;
            continue;
        }

        status = solver_codekd_packet_retire_nonwindow_descriptor(
            packet, context, descriptor, slot);
        if (status != INDEX_SHARD_STAGED_RETIRE_OK) {
            return status;
        }
    }

    if (packet->candidate_cursor != packet->candidate_count) {
        if (packet->candidate_window_count) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        if (solver_codekd_packet_begin_candidate_window(packet)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        return INDEX_SHARD_STAGED_RETIRE_MORE;
    }
    if (packet->descriptors->has_final_rel_field_noise2) {
        context->solver->rel_field_noise2 =
            packet->descriptors->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            context->solver,
            packet->descriptors->trailing_numtries,
            packet->descriptors->trailing_cxdx,
            packet->descriptors->trailing_meanx)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    return INDEX_SHARD_STAGED_RETIRE_OK;
}

static index_shard_staged_retire_status_t
solver_codekd_search_packet_retire(
    const index_shard_staged_task_t* task,
    size_t task_index,
    void* owner_context) {
    solver_codekd_packet_retire_context_t* context = owner_context;
    const solver_codekd_packet_task_input_t* input;
    solver_codekd_search_packet_t* packet;
    size_t end;
    size_t descriptor_index;

    if (!task || !context || !context->solver ||
        !context->query_result || !context->reduced ||
        !task->input || task->input_bytes != sizeof(*input) ||
        !task->output || task->output_bytes != sizeof(*packet) ||
        task_index != context->next_task_index) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }
    input = task->input;
    packet = task->output;
    if (input->descriptor.combination_first != context->next_sequence ||
        input->descriptor.combination_first >=
            input->descriptor.combination_end ||
        packet->sequence != input->descriptor.combination_first ||
        packet->state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet->descriptors == NULL ||
        packet->count != packet->descriptors->descriptor_count ||
        packet->first != 0U ||
        packet->next_descriptor != packet->count ||
        packet->plan_complete || packet->delivery_ticket ||
        packet->pending_descriptor_plan ||
        packet->pending_descriptor_raw_ranges ||
        packet->pending_descriptor_logical_bytes ||
        packet->delivery_source ||
        (packet->candidate_count &&
         packet->candidate_count != packet->hit_count) ||
        (packet->candidate_count && !packet->candidate_records) ||
        (packet->candidate_count && !packet->candidate_capacity) ||
        packet->candidate_count > packet->hit_count ||
        packet->candidate_cursor > packet->candidate_count ||
        packet->candidate_window_first >
            packet->candidate_cursor ||
        packet->candidate_window_offset !=
            packet->candidate_cursor -
                packet->candidate_window_first ||
        packet->candidate_window_count >
            packet->candidate_capacity ||
        packet->candidate_window_count >
            packet->candidate_count -
                packet->candidate_window_first ||
        packet->candidate_window_offset >
            packet->candidate_window_count ||
        packet->candidate_quad_ready_count >
            packet->candidate_window_count ||
        packet->candidate_star_ready_count >
            packet->candidate_quad_ready_count ||
        packet->retire_descriptor > packet->count) {
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    packet->state = SOLVER_CODEKD_PACKET_RETIRING;
    context->solver->profile.max_batch_hypotheses = MAX(
        context->solver->profile.max_batch_hypotheses,
        packet->count);
    if (packet->candidate_count) {
        index_shard_staged_retire_status_t status =
            solver_codekd_search_packet_retire_window(
                packet, context);

        if (status == INDEX_SHARD_STAGED_RETIRE_OK) {
            packet->state = SOLVER_CODEKD_PACKET_RETIRED;
            context->next_task_index++;
            context->next_sequence =
                input->descriptor.combination_end;
        }
        return status;
    }
    end = packet->first + packet->count;
    for (descriptor_index = packet->first;
         descriptor_index < end;
         descriptor_index++) {
        const solver_ab_descriptor_t* descriptor =
            &packet->descriptors->descriptors[descriptor_index];
        const solver_codekd_result_slot_t* slot =
            &packet->slots[descriptor_index];

        if (slot->state != SOLVER_CODEKD_RESULT_QUERY_FAILED &&
            solver_poll_worker_stop(context->solver)) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_STAGED_RETIRE_STOPPED;
        }
        context->solver->rel_field_noise2 =
            descriptor->rel_field_noise2;
        if (solver_ab_checked_counter_delta(
                context->solver,
                descriptor->numtries_delta,
                descriptor->cxdx_delta,
                descriptor->meanx_delta)) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }

        if (slot->state == SOLVER_CODEKD_RESULT_READY) {
            kdtree_qres_t result_view;
            solver_candidate_delivery_record_t* prepared = NULL;
            size_t prepared_quad_count = 0U;
            size_t prepared_star_count = 0U;

            if (slot->hit_first > packet->hit_count ||
                (size_t)slot->hit_count >
                    packet->hit_count - slot->hit_first) {
                packet->state = SOLVER_CODEKD_PACKET_FAILED;
                return INDEX_SHARD_STAGED_RETIRE_ERROR;
            }
            memset(&result_view, 0, sizeof(result_view));
            result_view.nres = slot->hit_count;
            result_view.capacity = slot->hit_count;
            result_view.inds = packet->inds + slot->hit_first;
            result_view.sdists = packet->sdists + slot->hit_first;
            if (slot->hit_first <
                packet->candidate_quad_ready_count) {
                prepared =
                    packet->candidate_records + slot->hit_first;
                prepared_quad_count = MIN(
                    (size_t)slot->hit_count,
                    packet->candidate_quad_ready_count -
                        slot->hit_first);
                if (slot->hit_first <
                    packet->candidate_star_ready_count) {
                    prepared_star_count = MIN(
                        (size_t)slot->hit_count,
                        packet->candidate_star_ready_count -
                            slot->hit_first);
                }
            }
            solver_execute_prepared_hypothesis_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                &result_view,
                slot->search_wall_seconds,
                prepared,
                prepared_quad_count,
                prepared_star_count);
        } else if (slot->state ==
                   SOLVER_CODEKD_RESULT_OWNER_REPLAY) {
            solver_execute_hypothesis_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                descriptor->tol2,
                context->query_result);
        } else if (slot->state ==
                   SOLVER_CODEKD_RESULT_QUERY_FAILED) {
            solver_retire_codekd_failure_owner(
                descriptor->stars,
                descriptor->code,
                context->dimquads,
                context->solver,
                descriptor->current_parity,
                slot->search_errno,
                slot->search_wall_seconds);
        } else {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }

        (*context->reduced)++;
        if (context->solver->profile.execution_failed) {
            packet->state = SOLVER_CODEKD_PACKET_FAILED;
            return INDEX_SHARD_STAGED_RETIRE_ERROR;
        }
        if (context->solver->quit_now) {
            packet->state = SOLVER_CODEKD_PACKET_STOPPED;
            return INDEX_SHARD_STAGED_RETIRE_STOPPED;
        }
    }

    if (packet->descriptors->has_final_rel_field_noise2) {
        context->solver->rel_field_noise2 =
            packet->descriptors->final_rel_field_noise2;
    }
    if (solver_ab_checked_counter_delta(
            context->solver,
            packet->descriptors->trailing_numtries,
            packet->descriptors->trailing_cxdx,
            packet->descriptors->trailing_meanx)) {
        packet->state = SOLVER_CODEKD_PACKET_FAILED;
        return INDEX_SHARD_STAGED_RETIRE_ERROR;
    }

    packet->state = SOLVER_CODEKD_PACKET_RETIRED;
    context->next_task_index++;
    context->next_sequence = input->descriptor.combination_end;
    return INDEX_SHARD_STAGED_RETIRE_OK;
}

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
    solver_ab_descriptor_planner_t* planner;
    solver_codekd_packet_wave_t* wave = NULL;
    kdtree_t* code_tree;
    fitsbin_t* delivery_source;
    solver_ab_pair_t* pairs = NULL;
    size_t pair_count = 0U;
    unsigned long long total_combinations = 0U;
    unsigned long long combination_cursor = 0U;
    unsigned long long reduced = 0U;
    unsigned long long parallel_reduced = 0U;
    size_t staged_capacity;
    size_t compute_width;
    size_t participants;
    size_t expansion;
    size_t max_task_combinations;
    anbool assisted = FALSE;
    anbool packet_cleanup_failed = FALSE;
    int result = 0;

    if (mode_out) {
        *mode_out = SOLVER_AB_MODE_NATIVE;
    }
    if (!solver || !field_geometry || !query_result ||
        !solver->index || !solver->index->codekd ||
        !solver->index->codekd->tree ||
        solver->maxquads != 0 ||
        solver->maxmatches != 0 ||
        pl_size(solver->indexes) != 1U ||
        !index_shard_worker_context_active()) {
        return 0;
    }
    code_tree = solver->index->codekd->tree;
    if (!code_tree->io || !code_tree->io_is_fitsbin) {
        return 0;
    }
    delivery_source = (fitsbin_t*)code_tree->io;
    if (fitsbin_payload_io_service_width() <= 0 ||
        fitsbin_payload_is_fully_resident(delivery_source) ||
        fitsbin_get_mmap_advice(delivery_source) !=
            FITSBIN_MMAP_ADVICE_RANDOM) {
        return 0;
    }
    staged_capacity = index_shard_staged_capacity();
    compute_width = index_shard_staged_compute_width();
    participants = MIN(staged_capacity, compute_width);
    if (!participants) {
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
    if (!workspace ||
        solver_ab_descriptor_workspace_reserve_outputs(
            workspace, participants)) {
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
        total_combinations >
            ULLONG_MAX / (unsigned long long)expansion ||
        total_combinations * (unsigned long long)expansion <
            SOLVER_AB_DESCRIPTOR_DELIVERY_MIN_HYPOTHESES) {
        solver_ab_descriptor_release_pairs(workspace);
        return 0;
    }
    wave = solver_codekd_packet_wave_create(participants);
    if (!wave) {
        solver_ab_descriptor_release_pairs(workspace);
        return 0;
    }

    while (combination_cursor < total_combinations) {
        solver_codekd_packet_task_input_t* inputs = wave->inputs;
        solver_codekd_search_packet_t* packets = wave->packets;
        index_shard_staged_task_t* tasks = wave->tasks;
        unsigned long long remaining =
            total_combinations - combination_cursor;
        unsigned long long wave_capacity =
            (unsigned long long)max_task_combinations *
            (unsigned long long)participants;
        unsigned long long wave_combinations =
            MIN(remaining, wave_capacity);
        unsigned long long expanded_work =
            wave_combinations * (unsigned long long)expansion;
        unsigned long long work_task_count =
            expanded_work /
                SOLVER_AB_DESCRIPTOR_TASK_MIN_HYPOTHESES +
            (expanded_work %
                 SOLVER_AB_DESCRIPTOR_TASK_MIN_HYPOTHESES
             ? 1ULL
             : 0ULL);
        unsigned long long base;
        unsigned long long remainder;
        unsigned long long task_cursor = combination_cursor;
        size_t task_count =
            (size_t)((wave_combinations +
                      max_task_combinations - 1U) /
                     max_task_combinations);
        size_t task_index;
        size_t candidate_budget_per_task;
        size_t prepared_packets = 0U;
        index_shard_helper_run_status_t run_status =
            INDEX_SHARD_HELPER_UNAVAILABLE;
        index_shard_staged_run_stats_t run_stats;
        solver_codekd_packet_retire_context_t retire_context;
        anbool run_inline = FALSE;
        anbool wave_assisted = FALSE;
        unsigned long long reduced_before = reduced;

        work_task_count = MIN(
            work_task_count,
            (unsigned long long)participants);
        work_task_count = MIN(
            work_task_count,
            wave_combinations);
        if (task_count < (size_t)work_task_count) {
            task_count = (size_t)work_task_count;
        }
        if (!task_count) {
            result = -1;
            goto fail;
        }
        candidate_budget_per_task =
            solver_verification_wave_memory_budget() / task_count;
        base = wave_combinations /
            (unsigned long long)task_count;
        remainder = wave_combinations %
            (unsigned long long)task_count;
        memset(
            inputs, 0,
            participants * sizeof(*inputs));
        memset(
            packets, 0,
            participants * sizeof(*packets));
        memset(
            tasks, 0,
            participants * sizeof(*tasks));
        memset(&run_stats, 0, sizeof(run_stats));
        memset(&retire_context, 0, sizeof(retire_context));
        retire_context.solver = solver;
        retire_context.dimquads = dimquads;
        retire_context.query_result = query_result;
        retire_context.reduced = &reduced;
        retire_context.next_sequence = combination_cursor;

        for (task_index = 0U;
             task_index < task_count;
             task_index++) {
            solver_ab_descriptor_task_input_t* descriptor_input =
                &inputs[task_index].descriptor;
            unsigned long long task_combinations =
                base + (task_index < remainder ? 1U : 0U);
            unsigned long long work_units;

            if (!task_combinations ||
                task_combinations > max_task_combinations) {
                result = -1;
                goto fail;
            }
            descriptor_input->field_geometry = field_geometry;
            descriptor_input->pairs = pairs;
            descriptor_input->pair_count = pair_count;
            descriptor_input->combination_first = task_cursor;
            task_cursor += task_combinations;
            descriptor_input->combination_end = task_cursor;
            descriptor_input->phase = phase;
            descriptor_input->newpoint = newpoint;
            descriptor_input->dimquads = dimquads;
            descriptor_input->parity = solver->parity;
            descriptor_input->cx_less_than_dx =
                solver->index->cx_less_than_dx;
            descriptor_input->meanx_less_than_half =
                solver->index->meanx_less_than_half;
            descriptor_input->cxdx_margin = solver->cxdx_margin;
            inputs[task_index].tree = code_tree;
            inputs[task_index].quads = solver->index->quads;
            inputs[task_index].starkd = solver->index->starkd;
            inputs[task_index].detailed = solver->profile.detailed;
            inputs[task_index].use_radec = solver->use_radec;
            inputs[task_index].field_minx = solver->field_minx;
            inputs[task_index].field_maxx = solver->field_maxx;
            inputs[task_index].field_miny = solver->field_miny;
            inputs[task_index].field_maxy = solver->field_maxy;
            inputs[task_index].abscale_low = solver->abscale_low;
            inputs[task_index].abscale_high = solver->abscale_high;
            inputs[task_index].funits_lower = solver->funits_lower;
            inputs[task_index].funits_upper = solver->funits_upper;
            if (solver->mo_template) {
                memcpy(
                    &inputs[task_index].verification.match_template,
                    solver->mo_template,
                    sizeof(*solver->mo_template));
            }
            inputs[task_index].verification.field = solver->vf;
            inputs[task_index].verification.index_cutnside =
                solver->index->cutnside;
            inputs[task_index].verification.indexid =
                solver->index->indexid;
            inputs[task_index].verification.healpix =
                solver->index->healpix;
            inputs[task_index].verification.hpnside =
                solver->index->hpnside;
            inputs[task_index].verification.index_jitter =
                solver->index->index_jitter;
            inputs[task_index].verification.verify_pix =
                solver->verify_pix;
            inputs[task_index].verification.distractor_ratio =
                solver->distractor_ratio;
            inputs[task_index].verification.logratio_bail_threshold =
                solver->logratio_bail_threshold;
            inputs[task_index].verification.logaccept = MIN(
                solver->logratio_tokeep,
                solver->logratio_totune);
            inputs[task_index].verification.logratio_stoplooking =
                solver->logratio_stoplooking;
            inputs[task_index].verification.distance_from_quad_bonus =
                solver->distance_from_quad_bonus;
            inputs[task_index].verification.enabled =
                solver->vf && !verify_datalog_enabled();

            work_units = task_combinations *
                (unsigned long long)expansion;
            tasks[task_index].input = &inputs[task_index];
            tasks[task_index].input_bytes = sizeof(inputs[task_index]);
            tasks[task_index].output = &packets[task_index];
            tasks[task_index].output_bytes = sizeof(packets[task_index]);
            tasks[task_index].work_units = work_units;
        }
        if (task_cursor !=
            combination_cursor + wave_combinations) {
            result = -1;
            goto fail;
        }

        solver->profile.hypothesis_batches++;
        solver->profile.task_ranges_planned += task_count;
        solver->profile.task_ranges_submitted += task_count;
        solver->profile.task_ranges_executed += task_count;
        solver->profile.max_task_ranges = MAX(
            solver->profile.max_task_ranges,
            task_count);

        if (task_count) {
            int packet_prepare_status = 0;

            for (task_index = 0U;
                 task_index < task_count;
                 task_index++) {
                packet_prepare_status =
                    solver_codekd_search_packet_prepare(
                        &packets[task_index],
                        &workspace->outputs[task_index],
                        inputs[task_index].tree,
                        inputs[task_index].quads,
                        inputs[task_index].starkd,
                        inputs[task_index].descriptor.dimquads,
                        inputs[task_index].use_radec,
                        candidate_budget_per_task,
                        inputs[task_index].descriptor.combination_first,
                        inputs[task_index].detailed);
                if (packet_prepare_status) {
                    break;
                }
                prepared_packets++;
            }
            if (packet_prepare_status < 0) {
                for (task_index = 0U;
                     task_index < prepared_packets;
                     task_index++) {
                    if (solver_codekd_search_packet_cleanup(
                            &packets[task_index])) {
                        packet_cleanup_failed = TRUE;
                    }
                }
                result = -1;
                goto fail;
            }
            if (packet_prepare_status > 0) {
                solver->profile.allocation_failures++;
                solver->profile.page_plan_allocation_refused++;
                for (task_index = 0U;
                     task_index < prepared_packets;
                     task_index++) {
                    if (solver_codekd_search_packet_cleanup(
                            &packets[task_index])) {
                        packet_cleanup_failed = TRUE;
                    }
                }
                if (packet_cleanup_failed) {
                    result = -1;
                    goto fail;
                }
                run_inline = TRUE;
            } else {
                run_status = index_shard_staged_run_ordered(
                    &solver_codekd_packet_staged_ops,
                    tasks,
                    task_count,
                    solver_codekd_search_packet_retire,
                    &retire_context,
                    &run_stats);
                for (task_index = 0U;
                     task_index < prepared_packets;
                     task_index++) {
                    solver_codekd_packet_profile_accumulate(
                        solver, &packets[task_index]);
                    if (solver_codekd_search_packet_cleanup(
                            &packets[task_index])) {
                        packet_cleanup_failed = TRUE;
                    }
                }
                solver->profile.staged_owner_claims +=
                    run_stats.owner_claims;
                solver->profile.staged_foreign_claims +=
                    run_stats.foreign_claims;
                solver->profile.staged_io_submitted +=
                    run_stats.io_submitted;
                solver->profile.staged_io_completed +=
                    run_stats.io_completed;
                solver->profile.max_staged_io_submitted = MAX(
                    solver->profile.max_staged_io_submitted,
                    run_stats.max_io_submitted);
                solver->profile.max_staged_compute_ready = MAX(
                    solver->profile.max_staged_compute_ready,
                    run_stats.max_compute_ready);

                if (packet_cleanup_failed) {
                    result = -1;
                    goto fail;
                }

                if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE) {
                    run_inline = TRUE;
                } else if (run_status == INDEX_SHARD_HELPER_STOPPED) {
                    (void)solver_poll_worker_stop(solver);
                    solver->quit_now = TRUE;
                    solver->profile.hypothesis_batches_stopped++;
                    result = 1;
                    goto cleanup;
                } else if (run_status != INDEX_SHARD_HELPER_OK) {
                    result = -1;
                    goto fail;
                } else {
                    if (run_stats.io_completed >
                            run_stats.io_submitted ||
                        run_stats.owner_claims +
                            run_stats.foreign_claims == 0U) {
                        result = -1;
                        goto fail;
                    }
                    if (run_stats.foreign_compute_executes) {
                        wave_assisted = TRUE;
                        assisted = TRUE;
                        solver->profile.parallel_batches++;
                        solver->profile.parallel_batches_observed++;
                        solver->profile.ab_helper_tasks =
                            solver_ab_saturating_add(
                                solver->profile.ab_helper_tasks,
                                run_stats.foreign_compute_executes);
                        solver->profile.max_parallel_ranges = MAX(
                            solver->profile.max_parallel_ranges,
                            run_stats.max_compute_running);
                    }
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
                const solver_ab_descriptor_task_input_t* descriptor_input =
                    &inputs[task_index].descriptor;
                index_shard_helper_task_status_t task_status;
                int owner_status;

                task_status = solver_ab_descriptor_helper_execute(
                    descriptor_input,
                    sizeof(*descriptor_input),
                    &workspace->outputs[task_index],
                    sizeof(workspace->outputs[task_index]));
                if (task_status == INDEX_SHARD_HELPER_TASK_STOPPED) {
                    (void)solver_poll_worker_stop(solver);
                    solver->profile.hypothesis_batches_stopped++;
                    result = 1;
                    goto cleanup;
                }
                if (task_status != INDEX_SHARD_HELPER_TASK_OK) {
                    result = -1;
                    goto fail;
                }
                owner_status = solver_codekd_descriptor_execute_owner(
                    &workspace->outputs[task_index],
                    solver,
                    dimquads,
                    query_result,
                    &reduced);
                if (owner_status < 0) {
                    result = -1;
                    goto fail;
                }
                if (owner_status > 0) {
                    solver->profile.hypothesis_batches_stopped++;
                    result = 1;
                    goto cleanup;
                }
                retire_context.next_task_index++;
                retire_context.next_sequence =
                    descriptor_input->combination_end;
            }
        }

        if (retire_context.next_task_index != task_count ||
            retire_context.next_sequence != task_cursor) {
            result = -1;
            goto fail;
        }

        if (wave_assisted) {
            parallel_reduced = solver_ab_saturating_add(
                parallel_reduced,
                reduced - reduced_before);
        }
        solver->profile.hypothesis_batches_completed++;
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
    if (mode_out && result > 0) {
        *mode_out = assisted
            ? SOLVER_AB_MODE_ASSISTED
            : SOLVER_AB_MODE_FLATTENED_OWNER;
    }
    solver_ab_descriptor_release_pairs(workspace);
    if (!packet_cleanup_failed) {
        solver_codekd_packet_wave_destroy(wave);
    }
    return result;
}

#if defined(TESTING_TRYPERMUTATIONS)
int solver_test_candidate_rolling_windows(void) {
    static const size_t expected_counts[] = { 128U, 101U };
    solver_codekd_search_packet_t packet;
    solver_codekd_result_slot_t slots[3];
    solver_candidate_delivery_record_t records[
        SOLVER_CANDIDATE_DELIVERY_LIMIT];
    size_t window_index;
    size_t retired = 0U;

    memset(&packet, 0, sizeof(packet));
    memset(slots, 0, sizeof(slots));
    memset(records, 0, sizeof(records));
    slots[0].hit_first = 0U;
    slots[0].hit_count = 17U;
    slots[0].state = SOLVER_CODEKD_RESULT_READY;
    slots[1].hit_first = 17U;
    slots[1].hit_count = 200U;
    slots[1].state = SOLVER_CODEKD_RESULT_READY;
    slots[2].hit_first = 217U;
    slots[2].hit_count = 12U;
    slots[2].state = SOLVER_CODEKD_RESULT_READY;
    packet.slots = slots;
    packet.count = 3U;
    packet.hit_count = 229U;
    packet.candidate_count = packet.hit_count;
    packet.candidate_capacity =
        SOLVER_CANDIDATE_DELIVERY_LIMIT;
    packet.candidate_records = records;

    for (window_index = 0U;
         window_index <
             sizeof(expected_counts) / sizeof(expected_counts[0]);
         window_index++) {
        size_t slot_end;

        if (solver_codekd_packet_begin_candidate_window(&packet) ||
            packet.candidate_window_count !=
                expected_counts[window_index] ||
            packet.candidate_window_count >
                packet.candidate_capacity) {
            return -1;
        }
        retired += packet.candidate_window_count;
        packet.candidate_cursor += packet.candidate_window_count;
        packet.candidate_window_offset =
            packet.candidate_window_count;
        solver_codekd_packet_clear_candidate_window(&packet);
        while (packet.retire_descriptor < packet.count) {
            slot_end =
                slots[packet.retire_descriptor].hit_first +
                slots[packet.retire_descriptor].hit_count;
            if (packet.candidate_cursor < slot_end) {
                break;
            }
            packet.retire_descriptor++;
        }
    }
    if (retired != packet.hit_count ||
        packet.candidate_cursor != packet.candidate_count ||
        packet.retire_descriptor != packet.count ||
        packet.candidate_delivery_windows !=
            sizeof(expected_counts) / sizeof(expected_counts[0])) {
        return -1;
    }
    return 0;
}

int solver_test_candidate_nonresident_zero_submit_falls_back(void) {
    fitsbin_t source;
    quadfile_t quads;
    solver_codekd_search_packet_t packet;
    solver_candidate_delivery_record_t record;
    uint32_t quadrow[DQMAX];
    u32 quadid = 0U;
    int status;

    memset(&source, 0, sizeof(source));
    memset(&quads, 0, sizeof(quads));
    memset(&packet, 0, sizeof(packet));
    memset(&record, 0, sizeof(record));
    memset(quadrow, 0, sizeof(quadrow));
    source.mmap_prefetch_enabled = FALSE;
    source.payload_fully_resident = FALSE;
    quads.numquads = 1U;
    quads.dimquads = DQMAX;
    quads.fb = &source;
    quads.quadarray = quadrow;
    packet.candidate_records = &record;
    packet.inds = &quadid;
    packet.candidate_count = 1U;
    packet.candidate_capacity = 1U;
    packet.candidate_window_count = 1U;
    packet.state = SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY;

    status = solver_codekd_packet_submit_candidate_pages(&packet);
    if (status != -1 ||
        packet.state != SOLVER_CODEKD_PACKET_RESULTS_READY ||
        packet.state == SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY ||
        packet.candidate_quad_fallback != 1U ||
        !packet.candidate_quad_delivery_disabled ||
        !packet.candidate_star_delivery_disabled ||
        packet.candidate_quad_submitted ||
        packet.candidate_quad_ready ||
        packet.delivery_ticket ||
        packet.delivery_source) {
        return -1;
    }
    return 0;
}

int solver_test_verification_packet_bounds(void) {
    solver_verification_packet_t packet;
    size_t impossible_candidates =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES /
            sizeof(solver_ab_candidate_t) +
        1U;
    solver_packet_reserve_result_t reserve_result;

    memset(&packet, 0, sizeof(packet));
    if (solver_verification_packet_reserve(
            &packet,
            1U) != SOLVER_PACKET_RESERVE_OK ||
        packet.candidate_capacity >
            SOLVER_AB_CANDIDATE_LIMIT_BYTES /
                sizeof(*packet.candidates)) {
        solver_verification_packet_free(&packet);
        return -1;
    }
    reserve_result = solver_verification_packet_reserve(
        &packet,
        impossible_candidates);
    if (reserve_result != SOLVER_PACKET_RESERVE_FULL ||
        packet.allocation_failed ||
        packet.candidate_capacity >
            SOLVER_AB_CANDIDATE_LIMIT_BYTES /
                sizeof(*packet.candidates)) {
        solver_verification_packet_free(&packet);
        return -1;
    }
    solver_verification_packet_free(&packet);
    return 0;
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

        {
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
        if (solver->startobj) {
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
            {
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
                anbool phase_pquads_prepared = FALSE;
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

        {
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

static void solver_begin_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity) {
    solver->profile.hypotheses_generated++;
    if (solver->profile.max_batch_hypotheses < 1U) {
        solver->profile.max_batch_hypotheses = 1U;
    }
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
    }
}

static void solver_reduce_codekd_result_owner(
    const int* stars,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    kdtree_qres_t* result,
    double search_wall_seconds,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    solver->profile.codekd_calls++;
    solver->profile.hypotheses_executed++;
    if (solver->profile.detailed) {
        uint64_t kd_result_digest =
            solver_kd_result_order_digest(result);

        solver->profile.codekd_wall_seconds += search_wall_seconds;
        solver->profile.kd_result_order_hash =
            solver_order_hash_mix(
                solver->profile.kd_result_order_hash,
                kd_result_digest);
    }
    solver->profile.codekd_hits +=
        (unsigned long long)result->nres;

    if (solver_poll_worker_stop(solver)) {
        return;
    }
    if (result->nres) {
        double pixvals[DQMAX * 2];
        double resolve_wall_start = 0.0;
        int j;

        for (j = 0; j < dimquad; j++) {
            setx(pixvals, j, field_getx(solver, stars[j]));
            sety(pixvals, j, field_gety(solver, stars[j]));
        }
        if (solver->profile.detailed) {
            resolve_wall_start = monotonic_seconds();
        }
        resolve_matches_with_delivery(
            result,
            pixvals,
            stars,
            dimquad,
            solver->numtries,
            solver,
            current_parity,
            prepared,
            0U,
            prepared_quad_count,
            prepared_star_count);
        solver->profile.resolve_calls++;
        if (solver->profile.detailed) {
            solver->profile.resolve_wall_seconds +=
                monotonic_seconds() - resolve_wall_start;
        }
    }
    solver->profile.hypotheses_reduced++;
}

static void solver_execute_prepared_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    kdtree_qres_t* result,
    double search_wall_seconds,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    if (!result || solver_poll_worker_stop(solver)) {
        return;
    }
    solver_begin_hypothesis_owner(
        stars, code, dimquad, solver, current_parity);
    solver_reduce_codekd_result_owner(
        stars,
        dimquad,
        solver,
        current_parity,
        result,
        search_wall_seconds,
        prepared,
        prepared_quad_count,
        prepared_star_count);
}

static void solver_retire_codekd_failure_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    int search_errno,
    double search_wall_seconds) {
    if (!solver) {
        return;
    }
    solver_begin_hypothesis_owner(
        stars, code, dimquad, solver, current_parity);
    solver->profile.codekd_calls++;
    solver->profile.hypotheses_executed++;
    if (solver->profile.detailed) {
        solver->profile.codekd_wall_seconds += search_wall_seconds;
    }
    errno = search_errno;
    if (search_errno == ECANCELED &&
        solver_poll_worker_stop(solver)) {
        solver->quit_now = TRUE;
        return;
    }
    solver->profile.search_failures++;
    solver->profile.execution_failed = TRUE;
    solver->quit_now = TRUE;
}

static void solver_execute_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    kdtree_qres_t** presult) {
    double search_wall_seconds = 0.0;
    double search_wall_start = 0.0;

    if (solver_poll_worker_stop(solver)) {
        return;
    }
    solver_begin_hypothesis_owner(
        stars, code, dimquad, solver, current_parity);
    if (solver->profile.detailed) {
        search_wall_start = monotonic_seconds();
    }
    *presult = solver_codekd_rangesearch(
        solver->index->codekd->tree,
        *presult,
        code,
        tol2,
        SOLVER_CODEKD_SEARCH_OPTIONS);
    if (solver->profile.detailed) {
        search_wall_seconds =
            monotonic_seconds() - search_wall_start;
    }

    if (!*presult) {
        solver->profile.codekd_calls++;
        solver->profile.hypotheses_executed++;
        if (solver->profile.detailed) {
            solver->profile.codekd_wall_seconds +=
                search_wall_seconds;
        }
        if (errno == ECANCELED &&
            solver_poll_worker_stop(solver)) {
            solver->quit_now = TRUE;
            return;
        }
        solver->profile.search_failures++;
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
        return;
    }

    solver_reduce_codekd_result_owner(
        stars,
        dimquad,
        solver,
        current_parity,
        *presult,
        search_wall_seconds,
        NULL,
        0U,
        0U);
}

static void solver_index_payload_failure(
    solver_t* solver,
    const char* component) {
    logerr("[solver-io] failed to decode %s payload\n",
           component);
    solver->profile.execution_failed = TRUE;
    solver->quit_now = TRUE;
}

static void resolve_matches_native_range(
    kdtree_qres_t* krez,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    int jj, thisquadno;
    MatchObj mo;

    for (jj = candidate_first; jj < candidate_end; jj++) {
        unsigned int star[DQMAX];
        double starxyz[DQMAX*3];
        double scale;
        double arcsecperpix;
        tan_t wcs;
        int i;
        anbool outofbounds = FALSE;
        anbool prepared_stars = FALSE;
        anbool prepared_candidate = FALSE;
        double abscale;
        solver_candidate_delivery_record_t* record = NULL;

        if (solver_poll_worker_stop(solver)) {
            return;
        }

        solver->nummatches++;
        thisquadno = krez->inds[jj];

        if (prepared && (size_t)jj >= prepared_first &&
            (size_t)jj - prepared_first < prepared_quad_count &&
            prepared[(size_t)jj - prepared_first].quadid ==
                (unsigned int)thisquadno) {
            record = &prepared[(size_t)jj - prepared_first];
            memcpy(star, record->stars,
                   (size_t)dimquads * sizeof(*star));
        } else if (quadfile_get_stars(
                       solver->index->quads,
                       thisquadno,
                       star)) {
            solver_index_payload_failure(
                solver, "QuadFile");
            return;
        }
        if (record && (size_t)jj >= prepared_first &&
            (size_t)jj - prepared_first <
                prepared_star_count) {
            prepared_stars = TRUE;
            if (record->candidate_prepared &&
                record->prepared_parity == current_parity &&
                !memcmp(
                    record->prepared_fieldstars,
                    fieldstars,
                    (size_t)dimquads *
                        sizeof(*fieldstars)) &&
                !memcmp(
                    record->prepared_fieldxy,
                    field_xy,
                    (size_t)dimquads * 2U *
                        sizeof(*field_xy))) {
                prepared_candidate = TRUE;
                memcpy(
                    starxyz,
                    record->starxyz,
                    (size_t)dimquads * 3U *
                        sizeof(*starxyz));
            }
        }
        if (record && !prepared_candidate) {
            solver_codekd_record_clear_verification_speculation(
                record);
        }


        if (solver->use_radec) {
            /*
             * Preserve the original all-star loading and rejection order
             * when the sky-position constraint depends on every quad star.
             */
            for (i = 0; i < dimquads; i++) {
                if (prepared_stars) {
                    memcpy(starxyz + 3 * i,
                           record->starxyz + 3 * i,
                           3U * sizeof(*starxyz));
                } else if (startree_get(
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
                    krez->sdists[jj]);
                solver->num_radec_skipped++;
                continue;
            }
        } else {
            /*
             * The quick scale gate uses only A and B. Consume C/D/E only
             * after the candidate has passed that gate.
             */
            if (prepared_stars) {
                memcpy(starxyz, record->starxyz,
                       6U * sizeof(*starxyz));
            } else if (startree_get(
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

        if (prepared_candidate) {
            switch (record->plan_action) {
            case SOLVER_AB_CANDIDATE_RADEC_SKIP:
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_RADEC_SKIP,
                    thisquadno,
                    krez->sdists[jj]);
                solver->num_radec_skipped++;
                solver->profile.candidate_math_reused++;
                continue;

            case SOLVER_AB_CANDIDATE_ABSCALE_SKIP:
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_ABSCALE_SKIP,
                    thisquadno,
                    krez->sdists[jj]);
                solver->num_abscale_skipped++;
                solver->profile.candidate_math_reused++;
                continue;

            case SOLVER_AB_CANDIDATE_BAD_QUAD:
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_BAD_QUAD,
                    thisquadno,
                    krez->sdists[jj]);
                solver->profile.candidate_math_reused++;
                logverb("bad quad at %s:%i\n", __FILE__, __LINE__);
                continue;

            case SOLVER_AB_CANDIDATE_SCALE_SKIP:
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_SCALE_SKIP,
                    thisquadno,
                    krez->sdists[jj]);
                solver->profile.candidate_math_reused++;
                debug("          bad scale (%g arcsec/pix, range %g %g)\n",
                      record->prepared_scale,
                      solver->funits_lower,
                      solver->funits_upper);
                continue;

            case SOLVER_AB_CANDIDATE_VERIFY:
                memcpy(&wcs, &record->prepared_wcs, sizeof(wcs));
                arcsecperpix = record->prepared_scale;
                solver->profile.candidate_math_reused++;
                break;

            default:
                prepared_candidate = FALSE;
                break;
            }
        }
        if (record && !prepared_candidate) {
            solver_codekd_record_clear_verification_speculation(
                record);
        }
        if (!prepared_candidate) {
            // Quick-n-dirty scale estimate based on two stars.
            // in (rad per pix)**2
            abscale = square(distsq2rad(distsq(
                starxyz, starxyz + 3, 3))) /
                distsq(field_xy, field_xy + 2, 2);
            if (abscale > solver->abscale_high ||
                abscale < solver->abscale_low) {
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_ABSCALE_SKIP,
                    thisquadno,
                    krez->sdists[jj]);
                solver->num_abscale_skipped++;
                continue;
            }

            if (!solver->use_radec) {
                for (i = 2; i < dimquads; i++) {
                    if (prepared_stars) {
                        memcpy(starxyz + 3 * i,
                               record->starxyz + 3 * i,
                               3U * sizeof(*starxyz));
                    } else if (startree_get(
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
            if (fit_tan_wcs(
                    starxyz,
                    field_xy,
                    dimquads,
                    &wcs,
                    &scale)) {
                // bad quad.
                solver_record_candidate_order(
                    solver,
                    SOLVER_AB_CANDIDATE_BAD_QUAD,
                    thisquadno,
                    krez->sdists[jj]);
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
                    krez->sdists[jj]);
                debug("          bad scale (%g arcsec/pix, range %g %g)\n",
                      arcsecperpix,
                      solver->funits_lower,
                      solver->funits_upper);
                continue;
            }
        }
        solver_record_candidate_order(
            solver,
            SOLVER_AB_CANDIDATE_VERIFY,
            thisquadno,
            krez->sdists[jj]);
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

        if (record && record->verify_query &&
            (!prepared_candidate ||
             record->plan_action !=
                 SOLVER_AB_CANDIDATE_VERIFY)) {
            verify_destroy_index_query(record->verify_query);
            record->verify_query = NULL;
        }
        if (record && prepared_candidate &&
            record->plan_action ==
                SOLVER_AB_CANDIDATE_VERIFY &&
            record->prepared_verification &&
            record->verification_score_ready) {
            if (solver_handle_scored_query_hit(
                    solver,
                    &mo,
                    NULL,
                    FALSE,
                    record)) {
                solver->quit_now = TRUE;
            }
        } else if (record && prepared_candidate &&
            record->plan_action ==
                SOLVER_AB_CANDIDATE_VERIFY &&
            record->verify_query) {
            if (solver_handle_query_hit(
                    solver,
                    &mo,
                    NULL,
                    FALSE,
                    &record->verify_query)) {
                solver->quit_now = TRUE;
            }
        } else if (solver_handle_hit(
                       solver, &mo, NULL, FALSE)) {
            solver->quit_now = TRUE;
        }

        if (unlikely(solver->quit_now))
            return;
    }
}

static void resolve_matches_with_delivery(
    kdtree_qres_t* krez,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count) {
    size_t window_first;

    /*
     * field_xy contains [x_A, y_A, x_B, y_B, ...]. Cold payloads remain on
     * the original owner-local path. Prepared verification is permitted only
     * when its complete Quad/Star input is already resident.
     */
    assert(krez);
    assert(dimquads > 0);
    assert(dimquads <= DQMAX);

    if (!solver_payload_candidate_data_fully_resident(solver)) {
        resolve_matches_native_range(
            krez,
            field_xy,
            fieldstars,
            dimquads,
            quads_tried,
            solver,
            current_parity,
            0,
            krez->nres,
            prepared,
            prepared_first,
            prepared_quad_count,
            prepared_star_count);
        return;
    }

    window_first = 0U;
    while (window_first < (size_t)krez->nres) {
        size_t window_end = MIN(
            (size_t)krez->nres,
            window_first +
                (size_t)SOLVER_VERIFICATION_WINDOW_CANDIDATES);
        int handled = solver_ab_try_verification_wave(
            krez,
            field_xy,
            fieldstars,
            dimquads,
            quads_tried,
            solver,
            current_parity,
            (int)window_first,
            (int)window_end);

        if (!handled) {
            resolve_matches_native_range(
                krez,
                field_xy,
                fieldstars,
                dimquads,
                quads_tried,
                solver,
                current_parity,
                (int)window_first,
                (int)window_end,
                NULL,
                0U,
                0U,
                0U);
        }
        if (solver->quit_now || solver_poll_worker_stop(solver)) {
            return;
        }
        window_first = window_end;
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

/*
 * Continue one canonical owner verification from the exact native StarKD
 * query retained by the delivery packet. Allocation or preparation failure
 * falls back before MatchObj publication to the original verify_hit() path.
 */
static int solver_handle_query_hit(
    solver_t* sp,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match,
    verify_index_query_t** query) {
    verify_prepared_hit_t* prepared = NULL;
    verify_prepared_score_t score;
    double match_distance_in_pixels2;
    double verify_wall_start = 0.0;
    double logaccept;
    int status;

    assert(sp);
    assert(mo);

    if (!query || !*query) {
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }
    if (solver_poll_worker_stop(sp)) {
        verify_destroy_index_query(*query);
        *query = NULL;
        return FALSE;
    }

    match_distance_in_pixels2 =
        solver_prepare_hit_for_verify(
            sp, mo, &logaccept);
    memset(&score, 0, sizeof(score));
    if (sp->profile.detailed) {
        verify_wall_start = monotonic_seconds();
    }
    status = verify_prepare_hit_from_query(
        sp->index->starkd,
        query,
        sp->index->cutnside,
        mo,
        verifysip,
        sp->vf,
        match_distance_in_pixels2,
        sp->distractor_ratio,
        sp->field_maxx,
        sp->field_maxy,
        sp->logratio_bail_threshold,
        logaccept,
        sp->logratio_stoplooking,
        sp->distance_from_quad_bonus,
        fake_match,
        &prepared);
    if (status ||
        verify_score_prepared_hit(prepared, &score)) {
        verify_destroy_prepared_score(&score);
        verify_destroy_prepared_hit(prepared);
        verify_destroy_index_query(*query);
        *query = NULL;
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }
    if (verify_finish_prepared_hit(
            prepared, &score, mo)) {
        verify_destroy_prepared_score(&score);
        verify_destroy_prepared_hit(prepared);
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }
    verify_destroy_prepared_hit(prepared);
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

/*
 * Finish one already scored immutable verification context on the solver
 * owner. Any changed live verification parameter invalidates speculation and
 * replays the exact native path before scientific state is published.
 */
static int solver_handle_scored_query_hit(
    solver_t* sp,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match,
    solver_candidate_delivery_record_t* record) {
    double match_distance_in_pixels2;
    double verify_wall_start = 0.0;
    double logaccept;

    assert(sp);
    assert(mo);

    if (!record || fake_match ||
        record->verify_query ||
        !record->prepared_verification ||
        !record->verification_score_ready) {
        if (record) {
            if (record->verification_score_ready) {
                sp->profile.verification_score_fallback_batches++;
            }
            solver_codekd_record_clear_prepared_verification(record);
        }
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }
    if (solver_poll_worker_stop(sp)) {
        solver_codekd_record_clear_prepared_verification(record);
        return FALSE;
    }

    match_distance_in_pixels2 =
        solver_prepare_hit_for_verify(
            sp, mo, &logaccept);
    if (memcmp(
            mo->center,
            record->verify_center,
            sizeof(mo->center)) ||
        memcmp(
            &mo->radius,
            &record->verify_radius,
            sizeof(mo->radius)) ||
        memcmp(
            &match_distance_in_pixels2,
            &record->prepared_verify_pix2,
            sizeof(match_distance_in_pixels2)) ||
        memcmp(
            &logaccept,
            &record->prepared_logaccept,
            sizeof(logaccept)) ||
        memcmp(
            &sp->distractor_ratio,
            &record->prepared_distractor_ratio,
            sizeof(sp->distractor_ratio)) ||
        memcmp(
            &sp->logratio_bail_threshold,
            &record->prepared_logratio_bail_threshold,
            sizeof(sp->logratio_bail_threshold)) ||
        memcmp(
            &sp->logratio_stoplooking,
            &record->prepared_logratio_stoplooking,
            sizeof(sp->logratio_stoplooking)) ||
        memcmp(
            &sp->field_maxx,
            &record->prepared_field_maxx,
            sizeof(sp->field_maxx)) ||
        memcmp(
            &sp->field_maxy,
            &record->prepared_field_maxy,
            sizeof(sp->field_maxy)) ||
        sp->distance_from_quad_bonus !=
            record->prepared_distance_from_quad_bonus) {
        sp->profile.verification_score_fallback_batches++;
        solver_codekd_record_clear_prepared_verification(record);
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }

    if (sp->profile.detailed) {
        verify_wall_start = monotonic_seconds();
    }
    if (verify_finish_prepared_hit(
            record->prepared_verification,
            &record->prepared_score,
            mo)) {
        sp->profile.verification_score_fallback_batches++;
        if (sp->profile.detailed) {
            sp->profile.verify_wall_seconds +=
                monotonic_seconds() - verify_wall_start;
        }
        solver_codekd_record_clear_prepared_verification(record);
        return solver_handle_hit(
            sp, mo, verifysip, fake_match);
    }
    solver_codekd_record_clear_prepared_verification(record);
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
