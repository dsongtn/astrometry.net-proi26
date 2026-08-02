/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "os-features.h"
#include "verify.h"
#include "permutedsort.h"
#include "mathutil.h"
#include "keywords.h"
#include "log.h"
#include "sip-utils.h"
#include "healpix.h"
#include "datalog.h"
#include "verify_internal.h"
#include "verify_theta_tail.h"

#define DEBUGVERIFY 0

#if DEBUGVERIFY
#define debug2(args...) logdebug(args)
#else
#define debug2(args...)
#endif

#define DATALOG_MASK_VERIFY 0x1

// level
#define DLOG_ODDS 10

anbool verify_datalog_enabled(void) {
    return data_log_passes(DATALOG_MASK_VERIFY, DLOG_ODDS);
}

static anbool* verify_deduplicate_field_stars(
    verify_t* v, const verify_field_t* vf, double nsigmas);

static int verify_uniformize_field_checked(
    const double* xy,
    int* perm,
    int N,
    double fieldW,
    double fieldH,
    int nw,
    int nh,
    int** p_bincounts,
    int** p_binids);

verify_field_t* verify_field_preprocess(const starxy_t* fieldxy) {
    verify_field_t* vf;
    int Nleaf = 5;

    vf = calloc(1, sizeof(verify_field_t));
    if (!vf) {
        fprintf(stderr, "Failed to allocate space for a verify_field_t().\n");
        return NULL;
    }
    vf->field = fieldxy;
    // Note on kdtree type: I tried U32 (duu) but it was marginally slower.
    // I didn't try U16 (dss) because we need a fair bit of accuracy here.
    // Make a copy of the field objects, because we're going to build a
    // kdtree out of them and that shuffles their order.
    vf->fieldcopy = starxy_copy_xy(fieldxy);
    vf->xy = starxy_copy_xy(fieldxy);
    if (!vf->fieldcopy || !vf->xy) {
        fprintf(stderr, "Failed to copy the field.\n");
        verify_field_free(vf);
        return NULL;
    }
    // Build a tree out of the field objects (in pixel space)
    vf->ftree = kdtree_build(NULL, vf->fieldcopy, starxy_n(vf->field),
                             2, Nleaf, KDTT_DOUBLE, KD_BUILD_SPLIT);
    if (!vf->ftree) {
        fprintf(stderr, "Failed to build the verification field tree.\n");
        verify_field_free(vf);
        return NULL;
    }

    vf->do_uniformize = TRUE;
    vf->do_dedup = TRUE;
    vf->do_ror = TRUE;

    return vf;
}

void verify_field_free(verify_field_t* vf) {
    if (!vf) {
        return;
    }
    kdtree_free(vf->ftree);
    free(vf->xy);
    free(vf->fieldcopy);
    free(vf);
}

static double get_sigma2_at_radius(double verify_pix2, double r2, double quadr2) {
    return verify_pix2 * (1.0 + r2/quadr2);
}

static double* compute_sigma2s(const verify_field_t* vf,
                               const double* xy, int NF,
                               const double* qc, double Q2,
                               double verify_pix2, anbool do_gamma) {
    double* sigma2s;
    int i;
    double R2;

    if (NF < 0 ||
        (size_t)NF > SIZE_MAX / sizeof(double)) {
        return NULL;
    }
    sigma2s = malloc((size_t)NF * sizeof(double));
    if (NF && !sigma2s) {
        return NULL;
    }
    if (!do_gamma) {
        for (i=0; i<NF; i++)
            sigma2s[i] = verify_pix2;
    } else {
        // Compute individual positional variances for every field
        // star.
        for (i=0; i<NF; i++) {
            if (vf) {
                double sxy[2];
                starxy_get(vf->field, i, sxy);
                // Distance from the quad center of this field star:
                R2 = distsq(sxy, qc, 2);
            } else
                R2 = distsq(xy + 2*i, qc, 2);

            // Variance of a field star at that distance from the quad center:
            sigma2s[i] = get_sigma2_at_radius(verify_pix2, R2, Q2);
        }
    }
    return sigma2s;
}

double* verify_compute_sigma2s(const verify_field_t* vf, const MatchObj* mo,
                               double verify_pix2, anbool do_gamma) {
    int NF;
    double qc[2];
    double Q2=0;
    NF = starxy_n(vf->field);
    if (do_gamma) {
        verify_get_quad_center(vf, mo, qc, &Q2);
        debug2("Quad radius = %g pixels\n", sqrt(Q2));
    }
    return compute_sigma2s(vf, NULL, NF, qc, Q2, verify_pix2, do_gamma);
}

double* verify_compute_sigma2s_arr(const double* xy, int NF,
                                   const double* qc, double Q2,
                                   double verify_pix2, anbool do_gamma) {
    return compute_sigma2s(NULL, xy, NF, qc, Q2, verify_pix2, do_gamma);
}

static int get_xy_bin(const double* xy,
                      double fieldW, double fieldH,
                      int nw, int nh) {
    int ix, iy;
    ix = (int)floor(nw * xy[0] / fieldW);
    ix = MAX(0, MIN(nw-1, ix));
    iy = (int)floor(nh * xy[1] / fieldH);
    iy = MAX(0, MIN(nh-1, iy));
    return iy * nw + ix;
}

static void print_test_perm(verify_t* v) {
    int i;
    for (i=0; i<v->NTall; i++) {
        if (i == v->NT)
            debug2("(NT)");
        debug2("%i ", v->testperm[i]);
    }
}

int verify_internal_get_test_stars(
    verify_t* v,
    const verify_field_t* vf,
    MatchObj* mo,
    double pix2,
    anbool do_gamma,
    anbool fake_match) {
    anbool* keepers = NULL;
    int i;
    int ibad=0, igood=0;

    v->NTall = starxy_n(vf->field);
    if (v->NTall < 0 ||
        (size_t)v->NTall > SIZE_MAX / sizeof(int)) {
        return -1;
    }
    v->testxy = vf->xy;
    v->NT = v->NTall;
    v->testsigma = verify_compute_sigma2s(vf, mo, pix2, do_gamma);
    v->testperm = permutation_init(NULL, v->NTall);
    v->tbadguys = malloc((size_t)v->NTall * sizeof(int));
    if ((v->NTall && !v->testsigma) ||
        (v->NTall && !v->testperm) ||
        (v->NTall && !v->tbadguys)) {
        return -1;
    }

    if (DEBUGVERIFY) {
        debug2("start:\n");
        print_test_perm(v);
        debug2("\n");
    }

    if (vf->do_dedup) {
        // Deduplicate test stars.  This could be done (approximately) in preprocessing.
        // FIXME -- this should be at the reference deduplication radius, not relative to sigma!
        // -- this requires the match scale
        // -- can perhaps discretize dedup to nearest power-of-sqrt(2) pixel radius and cache it.
        // -- we can compute sigma much later
        keepers = verify_deduplicate_field_stars(v, vf, 1.0);
        if (v->NTall && !keepers) {
            return -1;
        }

        // Remove test quad stars.  Do this after deduplication so we
        // don't end up with (duplicate) test stars near the quad stars.
        if (!fake_match) {
            for (i=0; i<mo->dimquads; i++) {
                assert(mo->field[i] >= 0);
                assert(mo->field[i] < v->NTall);
                keepers[mo->field[i]] = FALSE;
            }
        }

        ibad = igood = 0;
        for (i=0; i<v->NT; i++) {
            int ti = v->testperm[i];
            if (keepers[ti]) {
                v->testperm[igood] = ti;
                igood++;
            } else {
                v->tbadguys[ibad] = ti;
                ibad++;
            }
        }
    } else {
        // Remove the quad.
        if (!fake_match) {
            int j;
            for (i=0; i<mo->dimquads; i++) {
                assert(mo->field[i] >= 0);
                assert(mo->field[i] < v->NTall);
            }
            ibad = igood = 0;
            for (i=0; i<v->NT; i++) {
                int ti = v->testperm[i];
                anbool isquad = FALSE;
                for (j=0; j<mo->dimquads; j++) {
                    if (ti == mo->field[j]) {
                        isquad = TRUE;
                        break;
                    }
                }
                if (!isquad) {
                    v->testperm[igood] = ti;
                    igood++;
                } else {
                    v->tbadguys[ibad] = ti;
                    ibad++;
                }
            }
        } else {
            igood = v->NT;
        }
    }

    v->NT = igood;
    // remember the bad guys
    if (ibad) {
        memcpy(v->testperm + igood, v->tbadguys,
               (size_t)ibad * sizeof(int));
    }
    free(keepers);

    if (DEBUGVERIFY) {
        debug2("after dedup and removing quad:\n");
        print_test_perm(v);
        debug2("\n");
    }

    return 0;
}

double verify_get_ror2(double Q2, double area,
                       double distractors, int NR, double pix2) {
    return Q2 * MAX(1, (area*(1 - distractors) / (4. * M_PI * NR * pix2) - 1));
}

int verify_internal_apply_ror(verify_t* v,
                              int index_cutnside,
                              MatchObj* mo,
                              const verify_field_t* vf,
                              double pix2,
                              double distractors,
                              double fieldW,
                              double fieldH,
                              anbool do_gamma, anbool fake_match,
                              double* p_effA,
                              int* p_uninw, int* p_uninh) {
    int i;
    int uni_nw = 0, uni_nh = 0;
    double effA = fieldW * fieldH;
    double qc[2], Q2=0;
    int igood, ibad;
    int* binids = NULL;
    double* bincenters = NULL;

    // If we're verifying an existing WCS solution, then don't increase the variance
    // away from the center of the matched quad.
    if (fake_match)
        do_gamma = FALSE;

    if (verify_internal_get_test_stars(
            v, vf, mo, pix2, do_gamma, fake_match)) {
        goto fail;
    }
    debug2("Number of test stars: %i\n", v->NT);
    debug2("Number of reference stars: %i\n", v->NR);

    if (!fake_match)
        verify_get_quad_center(vf, mo, qc, &Q2);

    // Uniformize test stars
    // FIXME - can do this (possibly at several scales) in preprocessing.
    if (vf->do_uniformize) {
        // -get uniformization scale.
        verify_get_uniformize_scale(index_cutnside, mo->scale, fieldW, fieldH, &uni_nw, &uni_nh);
        debug2("uniformizing into %i x %i blocks.\n", uni_nw, uni_nh);

        // uniformize!
        if (uni_nw > 1 || uni_nh > 1) {
            if (verify_uniformize_field_checked(
                    vf->xy, v->testperm, v->NT,
                    fieldW, fieldH, uni_nw, uni_nh,
                    NULL, &binids)) {
                goto fail;
            }
            bincenters = verify_uniformize_bin_centers(fieldW, fieldH, uni_nw, uni_nh);
            if (!bincenters) {
                goto fail;
            }

            if (DEBUGVERIFY) {
                debug2("after uniformizing:\n");
                print_test_perm(v);
                debug2("\n");
            }
        }
    }
    if (vf->do_ror && !fake_match) {
        anbool* goodbins = NULL;
        int Ngoodbins;
        double ror2;

        debug2("Quad radius = %g\n", sqrt(Q2));
        ror2 = verify_get_ror2(Q2, fieldW*fieldH, distractors, v->NR, pix2);
        debug2("(strong) Radius of relevance is %.1f\n", sqrt(ror2));

        if (binids) {
            assert(uni_nw);
            goodbins = malloc((size_t)uni_nw * (size_t)uni_nh * sizeof(anbool));
            if (!goodbins) {
                goto fail;
            }
            Ngoodbins = 0;
            for (i=0; i<(uni_nw * uni_nh); i++) {
                double binr2 = distsq(bincenters + 2*i, qc, 2);
                goodbins[i] = (binr2 < ror2);
                if (goodbins[i])
                    Ngoodbins++;
            }
            // Remove test stars in irrelevant bins...
            igood = ibad = 0;
            for (i=0; i<v->NT; i++) {
                int ti = v->testperm[i];
                if (goodbins[binids[i]]) {
                    v->testperm[igood] = ti;
                    igood++;
                } else {
                    v->tbadguys[ibad] = ti;
                    ibad++;
                }
            }
        } else {
            // Remove test stars outside the RoR.
            igood = ibad = 0;
            for (i=0; i<v->NT; i++) {
                int ti = v->testperm[i];
                double r2 = distsq(qc, vf->xy + 2*ti, 2);
                if (r2 < ror2) {
                    v->testperm[igood] = ti;
                    igood++;
                } else {
                    v->tbadguys[ibad] = ti;
                    ibad++;
                }
            }
            // Count good bins to find effective area... (ugh)
            assert(!bincenters);
            if (!uni_nw)
                verify_get_uniformize_scale(index_cutnside, mo->scale, fieldW, fieldH, &uni_nw, &uni_nh);
            bincenters = verify_uniformize_bin_centers(fieldW, fieldH, uni_nw, uni_nh);
            if (!bincenters) {
                goto fail;
            }
            Ngoodbins = 0;
            for (i=0; i<(uni_nw * uni_nh); i++) {
                double binr2 = distsq(bincenters + 2*i, qc, 2);
                if (binr2 < ror2)
                    Ngoodbins++;
            }
        }

        v->NT = igood;
        if (ibad) {
            memcpy(v->testperm + igood, v->tbadguys,
                   (size_t)ibad * sizeof(int));
        }
        debug2("After removing %i/%i irrelevant bins: %i test stars.\n", (uni_nw*uni_nh)-Ngoodbins, uni_nw*uni_nh, v->NT);

        if (DEBUGVERIFY) {
            debug2("after applying RoR:\n");
            print_test_perm(v);
            debug2("\n");
        }

        // Effective area: A * proportion of good bins.
        effA *= Ngoodbins / (double)(uni_nw * uni_nh);

        // Remove reference stars in bad bins.
        igood = ibad = 0;
        if (goodbins) {
            assert(uni_nw);
            for (i=0; i<v->NR; i++) {
                int ri = v->refperm[i];
                int binid = get_xy_bin(v->refxy + 2*ri, fieldW, fieldH, uni_nw, uni_nh);
                if (goodbins[binid]) {
                    v->refperm[igood] = ri;
                    igood++;
                } else {
                    v->badguys[ibad] = ri;
                    ibad++;
                }
            }
        } else {
            for (i=0; i<v->NR; i++) {
                int ri = v->refperm[i];
                if (distsq(qc, v->refxy + 2*ri, 2) < ror2) {
                    v->refperm[igood] = ri;
                    igood++;
                } else {
                    v->badguys[ibad] = ri;
                    ibad++;
                }
            }
        }
        // remember the bad guys
        if (ibad) {
            memcpy(v->refperm + igood, v->badguys,
                   (size_t)ibad * sizeof(int));
        }
        v->NR = igood;
        debug2("After removing irrelevant ref stars: %i ref stars.\n", v->NR);

        // New ROR is...
        debug2("ROR changed from %g to %g\n", sqrt(ror2),
               sqrt(verify_get_ror2(Q2, effA, distractors, v->NR, pix2)));

        free(goodbins);
    }
    free(bincenters);
    free(binids);

    *p_effA = effA;
    if (p_uninw)
        *p_uninw = uni_nw;
    if (p_uninh)
        *p_uninh = uni_nh;
    return 0;

fail:
    free(bincenters);
    free(binids);
    return -1;
}

void verify_get_index_stars(const double* fieldcenter, double fieldr2,
                            const startree_t* skdt, const sip_t* sip, const tan_t* tan,
                            double fieldW, double fieldH,
                            double** p_indexradec,
                            double** indexpix, int** p_starids, int* p_nindex) {
    double* indxyz;
    int i, N, NI;
    int* sweep;
    int* starid;
    int* inbounds;
    int* perm;
    double* radec = NULL;

    assert(skdt->sweep);
    assert(p_nindex);
    assert(sip || tan);

    // Find all index stars within the bounding circle of the field.
    startree_search_for(skdt, fieldcenter, fieldr2, &indxyz, NULL, &starid, &N);

    if (!indxyz) {
        // no stars in range.
        *p_nindex = 0;
        return;
    }

    // Find index stars within the rectangular field.
    if (!verify_internal_filter_stars_in_field_parallel(
            sip, tan, indxyz, N, indexpix, &inbounds, &NI)) {
        inbounds = sip_filter_stars_in_field(
            sip, tan, indxyz, NULL, N, indexpix, NULL, &NI);
    }
    // Apply the permutation now, so that "indexpix" and "starid" stay in sync:
    // indexpix is already in the "inbounds" ordering.
    permutation_apply(inbounds, NI, starid, starid, sizeof(int));

    // Compute index RA,Decs if requested.
    if (p_indexradec) {
        radec = malloc(2 * NI * sizeof(double));
        for (i=0; i<NI; i++)
            // note that the "inbounds" permutation is applied to "indxyz" here.
            // we will apply the sweep permutation below.
            xyzarr2radecdegarr(indxyz + 3*inbounds[i], radec + 2*i);
        *p_indexradec = radec;
    }
    free(indxyz);
    free(inbounds);

    // Each index star has a "sweep number" assigned during index building;
    // it roughly represents a local brightness ordering.  Use this to sort the
    // index stars.
    sweep = malloc(NI * sizeof(int));
    for (i=0; i<NI; i++)
        sweep[i] = skdt->sweep[starid[i]];
    perm = permuted_sort(sweep, sizeof(int), compare_ints_asc, NULL, NI);
    free(sweep);

    if (indexpix) {
        permutation_apply(perm, NI, *indexpix, *indexpix, 2 * sizeof(double));
        *indexpix = realloc(*indexpix, NI * 2 * sizeof(double));
    }

    if (p_starids) {
        permutation_apply(perm, NI, starid, starid, sizeof(int));
        starid = realloc(starid, NI * sizeof(int));
        *p_starids = starid;
    } else
        free(starid);

    if (p_indexradec)
        permutation_apply(perm, NI, radec, radec, 2 * sizeof(double));

    free(perm);

    *p_nindex = NI;
}

/**
 If field objects are within "sigma" of each other (where sigma depends on the
 distance from the matched quad), then they are not very useful for verification.
 We filter out field stars within sigma of each other, taking only the brightest.

 Returns an array indicating which field stars should be kept.
 */
static anbool* verify_deduplicate_field_stars(verify_t* v, const verify_field_t* vf, double nsigmas) {
    anbool* keepers = NULL;
    int i, j, ti;
    kdtree_qres_t* res = NULL;
    double nsig2 = nsigmas*nsigmas;
    int options = KD_OPTIONS_NO_RESIZE_RESULTS | KD_OPTIONS_SMALL_RADIUS;

    if (!v || !vf || v->NTall < 0 || v->NT < 0 ||
        v->NT > v->NTall ||
        (v->NT && (!v->testperm || !v->testsigma))) {
        return NULL;
    }
    // default to FALSE
    keepers = calloc((size_t)v->NTall, sizeof(anbool));
    if (v->NTall && !keepers) {
        return NULL;
    }
    for (i=0; i<v->NT; i++) {
        ti = v->testperm[i];
        keepers[ti] = TRUE;
    }
    for (i=0; i<v->NT; i++) {
        double sxy[2];
        ti = v->testperm[i];
        if (!keepers[ti])
            continue;
        starxy_get(vf->field, ti, sxy);
        res = kdtree_rangesearch_options_reuse(vf->ftree, res, sxy, nsig2 * v->testsigma[ti], options);
        if (!res || res->nres < 0 ||
            (res->nres && !res->inds)) {
            kdtree_free_query(res);
            free(keepers);
            return NULL;
        }
        for (j=0; j<res->nres; j++) {
            int ind = res->inds[j];
            if (ind > i) {
                keepers[ind] = FALSE;
                if (DEBUGVERIFY) {
                    double otherxy[2];
                    starxy_get(vf->field, ind, otherxy);
                    logdebug("Field star %i at %g,%g: is close to field star %i at %g,%g.  dist is %g, sigma is %g\n",
                             i, sxy[0], sxy[1], ind, otherxy[0], otherxy[1],
                             sqrt(distsq(sxy, otherxy, 2)), sqrt(nsig2 * v->testsigma[ti]));
                }
            }
        }
    }
    kdtree_free_query(res);
    return keepers;
}

void verify_get_quad_center(const verify_field_t* vf, const MatchObj* mo, double* centerpix,
                            double* quadr2) {
    double Axy[2], Bxy[2];
    // Find the midpoint of AB of the quad in pixel space.
    starxy_get(vf->field, mo->field[0], Axy);
    starxy_get(vf->field, mo->field[1], Bxy);
    centerpix[0] = 0.5 * (Axy[0] + Bxy[0]);
    centerpix[1] = 0.5 * (Axy[1] + Bxy[1]);
    // Find the radius-squared of the quad = distsq(qc, A)
    *quadr2 = distsq(Axy, centerpix, 2);
}

void verify_get_uniformize_scale(int cutnside, double scale, int W, int H, int* cutnw, int* cutnh) {
    double cutarcsec, cutpix;
    cutarcsec = healpix_side_length_arcmin(cutnside) * 60.0;
    cutpix = cutarcsec / scale;
    debug2("cut nside: %i\n", cutnside);
    debug2("cut scale: %g arcsec\n", cutarcsec);
    debug2("match scale: %g arcsec/pix\n", scale);
    debug2("cut scale: %g pixels\n", cutpix);
    if (cutnw)
        *cutnw = MAX(1, (int)round(W / cutpix));
    if (cutnh)
        *cutnh = MAX(1, (int)round(H / cutpix));
}

static int verify_uniformize_field_checked(
    const double* xy,
    int* perm,
    int N,
    double fieldW,
    double fieldH,
    int nw,
    int nh,
    int** p_bincounts,
    int** p_binids) {
    int* workspace = NULL;
    int* bincount;
    int* binoffset;
    int* binwrite;
    int* binmembers;
    int* inputbins;
    int i,k,p;
    int activecount;
    int nbins;
    int* bincounts = NULL;
    int* binids = NULL;
    size_t workspace_count;
    size_t nbins_size;
    size_t n_size;

    if (p_bincounts) {
        *p_bincounts = NULL;
    }
    if (p_binids) {
        *p_binids = NULL;
    }
    if (N < 0 || nw <= 0 || nh <= 0 ||
        (N && (!xy || !perm)) ||
        !isfinite(fieldW) || !isfinite(fieldH) ||
        fieldW <= 0.0 || fieldH <= 0.0 ||
        nw > INT_MAX / nh) {
        return -1;
    }
    nbins = nw * nh;
    nbins_size = (size_t)nbins;
    n_size = (size_t)N;
    if (n_size > SIZE_MAX / 2U ||
        nbins_size > (SIZE_MAX - 2U * n_size) / 3U ||
        3U * nbins_size + 2U * n_size >
            SIZE_MAX / sizeof(int)) {
        return -1;
    }
    workspace_count = 3U * nbins_size + 2U * n_size;

    if (p_binids) {
        if (n_size > SIZE_MAX / sizeof(int)) {
            return -1;
        }
        binids = malloc(n_size * sizeof(int));
        if (n_size && !binids) {
            return -1;
        }
    }

    workspace = malloc(workspace_count * sizeof(int));
    if (workspace_count && !workspace) {
        goto fail;
    }
    bincount = workspace;
    binoffset = bincount + nbins;
    binwrite = binoffset + nbins;
    binmembers = binwrite + nbins;
    inputbins = binmembers + N;
    memset(bincount, 0, (size_t)nbins * sizeof(int));

    // Count stars in each bin.
    debug2("Test star bins:\n");
    for (i=0; i<N; i++) {
        int ind;
        int bin;
        ind = perm[i];
        bin = get_xy_bin(xy + 2*ind, fieldW, fieldH, nw, nh);
        debug2("%i ", bin);
        inputbins[i] = bin;
        bincount[bin]++;
    }
    debug2("\n");

    if (p_bincounts) {
        // note the bin occupancies.
        bincounts = malloc(nbins_size * sizeof(int));
        if (nbins_size && !bincounts) {
            goto fail;
        }
        memcpy(bincounts, bincount, nbins_size * sizeof(int));
    }

    // Lay out each bin contiguously while preserving input permutation order.
    p = 0;
    for (i=0; i<nbins; i++) {
        binoffset[i] = p;
        binwrite[i] = p;
        p += bincount[i];
    }
    assert(p == N);
    for (i=0; i<N; i++) {
        int ind;
        int bin;
        ind = perm[i];
        bin = inputbins[i];
        binmembers[binwrite[bin]] = ind;
        binwrite[bin]++;
    }

    /*
     * Make the same round-robin sweeps as the legacy nested loops, but keep
     * only non-empty bins in the active list. The old maxcount * nbins scan
     * repeatedly visited empty/exhausted bins and dominated deep verification
     * when a field occupied only a small part of a fine grid.
     *
     * binwrite is dead after the contiguous layout above, so reuse it as the
     * ordered active-bin list without another allocation. Stable in-place
     * compaction preserves ascending bin order in every sweep and therefore
     * preserves the exact output permutation.
     */
    activecount = 0;
    for (i=0; i<nbins; i++) {
        if (bincount[i] > 0) {
            binwrite[activecount++] = i;
        }
    }

    p=0;
    for (k=0; activecount > 0; k++) {
        int nextactive = 0;

        for (i=0; i<activecount; i++) {
            int binid = binwrite[i];

            assert(k < bincount[binid]);
            perm[p] = binmembers[binoffset[binid] + k];
            if (binids) {
                binids[p] = binid;
            }
            p++;

            if (k + 1 < bincount[binid]) {
                binwrite[nextactive++] = binid;
            }
        }
        activecount = nextactive;
    }
    assert(p == N);

    free(workspace);
    if (p_bincounts) {
        *p_bincounts = bincounts;
    }
    if (p_binids) {
        *p_binids = binids;
    }
    return 0;

fail:
    free(workspace);
    free(bincounts);
    free(binids);
    return -1;
}

void verify_uniformize_field(const double* xy,
                             int* perm,
                             int N,
                             double fieldW, double fieldH,
                             int nw, int nh,
                             int** p_bincounts,
                             int** p_binids) {
    (void)verify_uniformize_field_checked(
        xy, perm, N, fieldW, fieldH, nw, nh,
        p_bincounts, p_binids);
}

double* verify_uniformize_bin_centers(double fieldW, double fieldH,
                                      int nw, int nh) {
    int i,j;
    size_t count;
    double* bxy;

    if (nw <= 0 || nh <= 0 ||
        !isfinite(fieldW) || !isfinite(fieldH) ||
        nw > INT_MAX / nh ||
        (size_t)nw * (size_t)nh >
            SIZE_MAX / (2U * sizeof(double))) {
        return NULL;
    }
    count = (size_t)nw * (size_t)nh * 2U;
    bxy = malloc(count * sizeof(double));
    if (!bxy) {
        return NULL;
    }
    for (j=0; j<nh; j++)
        for (i=0; i<nw; i++) {
            bxy[(j * nw + i)*2 +0] = (i + 0.5) * fieldW / (double)nw;
            bxy[(j * nw + i)*2 +1] = (j + 0.5) * fieldH / (double)nh;
        }
    return bxy;
}

void verify_wcs(const startree_t* skdt,
                int index_cutnside,
                const sip_t* sip,
                const verify_field_t* vf,
                double verify_pix2,
                double distractors,
                double fieldW,
                double fieldH,
                double logbail,
                double logaccept,
                double logstoplooking,

                double* logodds,
                int* nfield, int* nindex,
                int* nmatch, int* nconflict, int* ndistractor
                // int** theta ?
                ) {
    MatchObj mo;

    memset(&mo, 0, sizeof(MatchObj));

    radecdeg2xyzarr(sip->wcstan.crval[0], sip->wcstan.crval[1], mo.center);
    mo.radius = arcsec2dist(hypot(fieldW, fieldH)/2.0 * sip_pixel_scale(sip));
    memcpy(&(mo.wcstan), &(sip->wcstan), sizeof(tan_t));
    mo.wcs_valid = TRUE;

    verify_hit(skdt, index_cutnside, &mo, sip, vf, verify_pix2,
               distractors, fieldW, fieldH, logbail, logaccept,
               logstoplooking, FALSE, TRUE);

    if (logodds)
        *logodds = mo.logodds;
    if (nfield)
        *nfield = mo.nfield;
    if (nindex)
        *nindex = mo.nindex;
    if (nmatch)
        *nmatch = mo.nmatch;
    if (nconflict)
        *nconflict = mo.nconflict;
    if (ndistractor)
        *ndistractor = mo.ndistractor;
}


void verify_internal_set_null_mo(MatchObj* mo) {
    mo->nfield = 0;
    mo->nmatch = 0;
    matchobj_compute_derived(mo);
    mo->logodds = -LARGE_VAL;
}

static void check_permutation(const int* perm, int N) {
    int i;
    int* counts = calloc(N, sizeof(int));
    for (i=0; i<N; i++) {
        assert(perm[i] >= 0);
        assert(perm[i] < N);
        counts[perm[i]]++;
    }
    for (i=0; i<N; i++) {
        assert(counts[i] == 1);
    }
    free(counts);
}

static void verify_permutation_apply_workspace(
    const int* perm,
    int count,
    void* array,
    size_t element_size,
    void* workspace) {
    const unsigned char* input = array;
    unsigned char* output = workspace;
    int i;

    if (!count) {
        return;
    }
    for (i = 0; i < count; i++) {
        memcpy(
            output + (size_t)i * element_size,
            input + (size_t)perm[i] * element_size,
            element_size);
    }
    memcpy(array, workspace, (size_t)count * element_size);
}

int verify_internal_fixup_theta(int* theta, double* allodds,
                                int ibailed, int istopped, verify_t* v,
                                int besti, int NRimage, double* refxyz,
                                int** p_etheta, double** p_eodds) {
    int* etheta = NULL;
    double* eodds = NULL;
    int* invrperm = NULL;
    unsigned char* permutation_workspace = NULL;
    size_t permutation_stride;
    int i, ti;

    if (!p_etheta || !p_eodds) {
        return -1;
    }
    *p_etheta = NULL;
    *p_eodds = NULL;
    if (!theta || !allodds || !v ||
        v->NT < 0 || v->NTall < 0 ||
        v->NRall < 0 || NRimage < 0 ||
        v->NT > v->NTall || NRimage > v->NRall ||
        (v->NTall && !v->testperm) ||
        (NRimage && (!v->refperm || !v->refxy)) ||
        ibailed < -1 || ibailed >= v->NT ||
        istopped < -1 || istopped >= v->NT ||
        (size_t)v->NTall > SIZE_MAX / sizeof(*etheta) ||
        (size_t)v->NTall > SIZE_MAX / sizeof(*eodds) ||
        (size_t)v->NRall > SIZE_MAX / sizeof(*invrperm)) {
        return -1;
    }
    if (verify_theta_mark_unprocessed(
            theta, v->NT, ibailed, istopped)) {
        return -1;
    }
    permutation_stride = refxyz
        ? 3U * sizeof(double)
        : 2U * sizeof(double);
    if ((size_t)NRimage >
        SIZE_MAX / permutation_stride) {
        return -1;
    }
    for (i = 0; i < NRimage; i++) {
        if (v->refperm[i] < 0 ||
            v->refperm[i] >= v->NRall) {
            return -1;
        }
    }
    for (i = 0; i < v->NTall; i++) {
        if (v->testperm[i] < 0 ||
            v->testperm[i] >= v->NTall) {
            return -1;
        }
    }
    for (i = 0; i < v->NT; i++) {
        if (theta[i] >= v->NRall) {
            return -1;
        }
    }

    if (v->NTall) {
        etheta = malloc((size_t)v->NTall * sizeof(*etheta));
        eodds = malloc((size_t)v->NTall * sizeof(*eodds));
    }
    if (v->NRall) {
        invrperm = malloc((size_t)v->NRall * sizeof(*invrperm));
    }
    if (NRimage) {
        permutation_workspace = malloc(
            (size_t)NRimage * permutation_stride);
    }
    if ((v->NTall && (!etheta || !eodds)) ||
        (v->NRall && !invrperm) ||
        (NRimage && !permutation_workspace)) {
        free(permutation_workspace);
        free(invrperm);
        free(eodds);
        free(etheta);
        return -1;
    }

#define BAD_PERM -1000000
    for (i = 0; i < v->NRall; i++) {
        invrperm[i] = BAD_PERM;
    }
    for (i = 0; i < NRimage; i++) {
        invrperm[v->refperm[i]] = i;
    }
    for (i = 0; i < v->NT; i++) {
        if (theta[i] >= 0 &&
            invrperm[theta[i]] == BAD_PERM) {
            free(permutation_workspace);
            free(invrperm);
            free(eodds);
            free(etheta);
            return -1;
        }
    }

    if (DEBUGVERIFY) {
        // The "testperm" permutation should be "complete".
        check_permutation(v->testperm, v->NTall);
        // "refperm" has vals < NRall in elements < NRimage.
        //check_permutation(v->refperm, NRimage);
        for (i=0; i<NRimage; i++) {
            assert(v->refperm[i] >= 0);
            assert(v->refperm[i] < v->NRall);
        }
    }

    // At this point, "theta[0]" is the *reference* star index
    // that was matched by the test star "v->testperm[0]".
    // Meanwhile, "v->refperm" lists all the valid reference stars.

    // We want to produce "etheta", which has elements parallel to
    // the test stars in their original (brightness) ordering; that is,
    // we want to eliminate the need for "v->testperm".

    if (DEBUGVERIFY) {
        for (i=0; i<v->NT; i++) {
            Unused int ri;
            if (i == besti)
                debug2("* ");
            debug2("Theta[%i] = %i", i, theta[i]);
            if (theta[i] < 0) {
                debug2("\n");
                continue;
            }
            ri = theta[i];
            ti = v->testperm[i];
            debug2(" (starid %i), testxy=(%.1f, %.1f), refxy=(%.1f, %.1f)\n",
                   (v->refstarid ? v->refstarid[ri] : -1000), v->testxy[ti*2+0], v->testxy[ti*2+1], v->refxy[ri*2+0], v->refxy[ri*2+1]);
        }
    }
    // Apply the "refperm" permutation, mostly to cut out the stars that
    // aren't in the image (we want to have "nindex" = "NRimage" = "NRall").
    // This requires computing the inverse perm so we can fix theta to match.

    // The reference stars include stars that are actually outside
    // the field; we want to collapse the reference star list,
    // which will renumber them.

    if (v->refstarid) {
        verify_permutation_apply_workspace(
            v->refperm, NRimage, v->refstarid,
            sizeof(int), permutation_workspace);
    }
    verify_permutation_apply_workspace(
        v->refperm, NRimage, v->refxy,
        2U * sizeof(double), permutation_workspace);
    if (refxyz) {
        verify_permutation_apply_workspace(
            v->refperm, NRimage, refxyz,
            3U * sizeof(double), permutation_workspace);
    }

    // New v->refstarid[i] is old v->refstarid[ v->refperm[i] ]

    if (DEBUGVERIFY) {
        for (i=0; i<v->NTall; i++)
            etheta[i] = BAD_PERM;
    }

    for (i=0; i<v->NT; i++) {
        ti = v->testperm[i];
        if (DEBUGVERIFY)
            // assert that we haven't touched this element yet.
            assert(etheta[ti] == BAD_PERM);
        if (theta[i] < 0) {
            etheta[ti] = theta[i];
            // No match -> no weight.
            eodds[ti] = -LARGE_VAL;
        } else {
            if (DEBUGVERIFY)
                assert(invrperm[theta[i]] != BAD_PERM);
            etheta[ti] = invrperm[theta[i]];
            eodds[ti] = allodds[i];
        }
    }

    free(permutation_workspace);
    free(invrperm);

    for (i=v->NT; i<v->NTall; i++) {
        ti = v->testperm[i];
        etheta[ti] = THETA_FILTERED;
        eodds[ti] = -LARGE_VAL;
    }

    if (DEBUGVERIFY) {
        // We should touch every element.
        for (i=0; i<v->NTall; i++)
            assert(etheta[i] != BAD_PERM);
        for (i=0; i<v->NTall; i++)
            if (etheta[i] >= 0)
                assert(etheta[i] < NRimage);
            else
                assert(etheta[i] == THETA_FILTERED ||
                       etheta[i] == THETA_DISTRACTOR ||
                       etheta[i] == THETA_CONFLICT ||
                       etheta[i] == THETA_BAILEDOUT ||
                       etheta[i] == THETA_STOPPEDLOOKING);

    }

    *p_etheta = etheta;
    *p_eodds = eodds;
    return 0;
}

void verify_count_hits(int* theta, int besti, int* p_nmatch, int* p_nconflict, int* p_ndistractor) {
    int i;
    int d, c, m;
    d = 0;
    c = 0;
    m = 0;
    for (i=0; i<=besti; i++) {
        if (theta[i] == THETA_DISTRACTOR)
            d++;
        else if (theta[i] == THETA_CONFLICT)
            c++;
        else
            m++;
    }
    if (p_nconflict) *p_nconflict = c;
    if (p_ndistractor) *p_ndistractor = d;
    if (p_nmatch) *p_nmatch = m;
}


static void verify_hit_original(const startree_t* skdt, int index_cutnside,
                                MatchObj* mo, const sip_t* sip,
                                const verify_field_t* vf,
                                double pix2, double distractors,
                                double fieldW, double fieldH,
                                double logbail, double logaccept,
                                double logstoplooking,
                                anbool do_gamma, anbool fake_match) {
    int i,j;
    double* fieldcenter;
    double fieldr2;
    double effA, K, worst;
    int besti;
    int* theta = NULL;
    double* allodds = NULL;
    sip_t thewcs;
    int ibad, igood;
    double* refxyz = NULL;
    int* sweep = NULL;
    verify_t the_v;
    verify_t* v = &the_v;
    int NRimage;
    int ibailed, istopped;
    anbool score_completed;

    assert(mo->wcs_valid || sip);
    assert(isfinite(logaccept));
    assert(isfinite(logbail));

    memset(v, 0, sizeof(verify_t));

    if (sip)
        v->wcs = sip;
    else {
        sip_wrap_tan(&mo->wcstan, &thewcs);
        v->wcs = &thewcs;
    }

    // center and radius of the field in xyz space:
    fieldcenter = mo->center;
    fieldr2 = square(mo->radius);

    // find index stars and project them into pixel coordinates.
    /*
     verify_get_index_stars(fieldcenter, fieldr2, skdt, sip, &(mo->wcstan),
     fieldW, fieldH, NULL, &refxy, &starids, &NR);
     */
    /*
     Gotta be a bit careful with reference stars:

     We want to be able to return a list of all the reference
     stars in the image, but during the verification process we
     want to apply some filtering of reference stars.  We
     therefore keep an int array ("refperm") of indices into the
     arrays of reference star quantities.  There are "NR" good
     stars, but "NRall" in total.  Thus operations on all the
     stars must go to "NRall" in the original arrays, but
     operations on good stars must go to "NR", using "refperm" to
     redirect.

     This means that "refperm" should remain a permutation array (ie,
     no duplicates), and each value should be less than "NRall"; when
     filtering out an index, it should get moved to the part of the
     array between "NR" and "NRall".  We use the "badguys" array to
     hold these indices temporarily.
     */
    assert(skdt->sweep);
    // Find all index stars within the bounding circle of the field.
    startree_search_for(skdt, fieldcenter, fieldr2, &refxyz, NULL, &v->refstarid, &v->NRall);
    debug2("%i reference stars in the bounding circle\n", v->NRall);
    if (!refxyz) {
        // no stars in range.
        logverb("No reference stars in the bounding circle\n");
        goto bailout;
    }
    //logverb("Found %i reference stars in the bounding circle\n", v->NRall);
    // Find index stars within the rectangular field.
    v->refxy = malloc(v->NRall * 2 * sizeof(double));
    v->refperm = malloc(v->NRall * sizeof(int));
    igood = 0;
    for (i=0; i<v->NRall; i++) {
        if (!sip_xyzarr2pixelxy(v->wcs, refxyz+i*3, v->refxy+i*2, v->refxy+i*2 +1) ||
            !sip_pixel_is_inside_image(v->wcs, v->refxy[i*2], v->refxy[i*2+1])) {
            continue;
        }
        v->refperm[igood] = i;
        igood++;
    }
    v->NR = igood;
    // We sort of want to forget about stars not within the image...
    // but we don't want to change NRall...
    NRimage = v->NR;
    // NOTE that at this point, v->refperm elements past NRimage are invalid
    // (ie, may contain repeats)

    // Sort by sweep #.
    // Each index star has a "sweep number" assigned during index building;
    // it roughly represents a local brightness ordering.  Use this to sort the
    // index stars.
    // (NOTE that here we do want "sweep" to be size "NRall"; only the
    // bottom "NRimage" of the "refperm" array will be accessed in the
    // permuted_sort below, so none of
    // the elements between NRimage and NRall will be touched.)
    sweep = malloc(v->NRall * sizeof(int));
    for (i=0; i<v->NRall; i++)
        sweep[i] = skdt->sweep[v->refstarid[i]];
    // Note here that we're passing in an existing permutation array; it
    // gets re-permuted during this call.
    permuted_sort(sweep, sizeof(int), compare_ints_asc, v->refperm, v->NR);
    free(sweep);
    sweep = NULL;
    debug2("Found %i reference stars.\n", v->NR);

    // "refstarids" are indices into the star kdtree and could be used to
    // retrieve "tag-along" data with, eg, startree_get_data_column().

    v->badguys = malloc(v->NR * sizeof(int));

    // remove reference stars that are part of the quad.
    if (!fake_match) {
        ibad = 0;
        igood = 0;
        for (i=0; i<v->NR; i++) {
            anbool inquad = FALSE;
            int ri = v->refperm[i];
            for (j=0; j<mo->dimquads; j++) {
                if (v->refstarid[ri] == mo->star[j]) {
                    inquad = TRUE;
                    //debug2("Skipping ref star index %i, starid %i: quad star %i\n", ri, v->refstarid[ri], j);
                    v->badguys[ibad] = ri;
                    ibad++;
                    break;
                }
            }
            if (inquad)
                continue;
            v->refperm[igood] = ri;
            igood++;
        }
        // remember the bad guys
        if (ibad) {
            memcpy(v->refperm + igood, v->badguys,
                   (size_t)ibad * sizeof(int));
        }
        v->NR = igood;
        debug2("After removing stars in the quad: %i reference stars.\n", v->NR);
    }

    if (!v->NR) {
        logverb("After removing quad stars: no reference stars\n");
        goto bailout;
    }

    ///// FIXME -- we could compute the RoR and search for ref stars
    // based on the quad center and RoR rather than the image center
    // and image radius.

    if (!fake_match) {
        if (verify_internal_apply_ror(v, index_cutnside, mo,
                             vf, pix2, distractors,
                             fieldW, fieldH,
                             do_gamma, fake_match,
                             &effA, NULL, NULL)) {
            goto bailout;
        }
        if (!v->NR) {
            logerr("After applying ROR, NR = 0!\n");
            goto bailout;
        }
    } else {
        if (verify_internal_get_test_stars(
                v, vf, mo, pix2, do_gamma, fake_match)) {
            goto bailout;
        }
        effA = fieldW * fieldH;
        debug2("Number of test stars: %i\n", v->NT);
    }
    if (!v->NR || !v->NT) {
        logverb("After applying RoR, NR=%i, NT=%i\n", v->NR, v->NT);
        goto bailout;
    }

    worst = -LARGE_VAL;
    K = verify_internal_star_lists(v, effA, distractors,
                               logbail, logstoplooking, &besti, &allodds, &theta, &worst,
                               &ibailed, &istopped, &score_completed);
    if (!score_completed) {
        goto bailout;
    }
    mo->logodds = K;
    mo->worstlogodds = worst;
    // NTall so that caller knows how big 'etheta' is.
    mo->nfield = v->NTall;
    // NRimage: only the stars inside the image bounds.
    mo->nindex = NRimage;

    if (log_get_level() >= LOG_ALL) {
        int nm, nc, nd;
        verify_count_hits(theta, besti, &nm, &nc, &nd);
        debug("verify: logodds %g, %i matches, %i conflicts, %i distractors after %i field objects.\n",
              K, nm, nc, nd, besti);
    }

    if (K >= logaccept) {
        int ri, ti;
        int* etheta;
        double* eodds;
        int nm, nc, nd;
        verify_count_hits(theta, besti, &nm, &nc, &nd);
        mo->nmatch = nm;
        mo->nconflict = nc;
        mo->ndistractor = nd;

        if (verify_internal_fixup_theta(
                theta, allodds, ibailed, istopped,
                v, besti, NRimage, refxyz,
                &etheta, &eodds)) {
            goto bailout;
        }

        // Reinsert the matched quad...
        if (!fake_match) {
            for (j=0; j<mo->dimquads; j++) {
                // the ref star should have been eliminated, so it
                // should be in the "bad" part of the array, but
                // search the whole thing anyway.
                for (i=0; i<NRimage; i++) {
                    ri = i;
                    if (v->refstarid[ri] == mo->star[j]) {
                        ti = mo->field[j];
                        assert(etheta[ti] == THETA_FILTERED);
                        etheta[ti] = ri;
                        eodds[ti] = LARGE_VAL;
                        debug2("Matched ref index %i (star %i) to test index %i; ref pos=(%.1f, %.1f), test pos=(%.1f, %.1f)\n",
                               ri, v->refstarid[ri], ti, v->refxy[ri*2+0], v->refxy[ri*2+1], v->testxy[ti*2+0], v->testxy[ti*2+1]);
                        break;
                    }
                }
            }
        }

        if (DEBUGVERIFY) {
            debug2("\n");
            for (i=0; i<v->NTall; i++) {
                debug2("ETheta[%i] = %i", i, etheta[i]);
                if (etheta[i] < 0) {
                    debug2(" (w=%g)\n", verify_logodds_to_weight(eodds[i]));
                    continue;
                }
                ri = etheta[i];
                ti = i;
                debug2(" (starid %i), testxy=(%.1f, %.1f), refxy=(%.1f, %.1f), logodds=%g, w=%g\n",
                       v->refstarid[ri], v->testxy[ti*2+0], v->testxy[ti*2+1], v->refxy[ri*2+0], v->refxy[ri*2+1],
                       eodds[i], verify_logodds_to_weight(eodds[i]));
            }
        }

        mo->theta = etheta;
        mo->matchodds = eodds;
        mo->refxyz = refxyz;
        refxyz = NULL;
        mo->refxy = v->refxy;
        v->refxy = NULL;
        mo->refstarid = v->refstarid;
        v->refstarid = NULL;
        mo->testperm = v->testperm;
        v->testperm = NULL;

        matchobj_compute_derived(mo);
    }

 cleanup:
    free(refxyz);
    free(theta);
    free(allodds);
    free(v->testperm);
    free(v->testsigma);
    free(v->tbadguys);
    free(v->refperm);
    free(v->refxy);
    free(v->refstarid);
    free(v->badguys);
    return;

 bailout:
    verify_internal_set_null_mo(mo);
    // uh oh, spaghetti-code-oh!
    goto cleanup;
}

void verify_hit(const startree_t* skdt, int index_cutnside, MatchObj* mo,
                const sip_t* sip, const verify_field_t* vf,
                double pix2, double distractors,
                double fieldW, double fieldH,
                double logbail, double logaccept, double logstoplooking,
                anbool do_gamma, anbool fake_match) {
    double fieldr2;

    fieldr2 = square(mo->radius);
    debug("Field center %g,%g,%g, radius2 %g\n",
          mo->center[0], mo->center[1], mo->center[2], fieldr2);
    if (log_get_level() >= LOG_VERB) {
        double ra;
        double dec;
        double r;

        xyzarr2radecdeg(mo->center, &ra, &dec);
        r = distsq2deg(fieldr2);
        debug("Field center RA,Dec %g,%g, radius %g deg\n",
              ra, dec, r);
    }

    verify_hit_original(skdt, index_cutnside, mo, sip, vf,
                        pix2, distractors, fieldW, fieldH,
                        logbail, logaccept, logstoplooking,
                        do_gamma, fake_match);
}

// Free the things we added to this mo.
void verify_free_matchobj(MatchObj* mo) {
    free(mo->refxyz);
    free(mo->refstarid);
    free(mo->refxy);
    free(mo->theta);
    free(mo->matchodds);
    free(mo->testperm);
    mo->testperm = NULL;
    mo->refxyz = NULL;
    mo->refstarid = NULL;
    mo->refxy = NULL;
    mo->theta = NULL;
    mo->matchodds = NULL;
}

void verify_matchobj_deep_copy(const MatchObj* mo, MatchObj* dest) {
    if (mo->refxyz) {
        dest->refxyz = malloc(mo->nindex * 3 * sizeof(double));
        memcpy(dest->refxyz, mo->refxyz, mo->nindex * 3 * sizeof(double));
    }
    if (mo->refxy) {
        dest->refxy = malloc(mo->nindex * 2 * sizeof(double));
        memcpy(dest->refxy, mo->refxy, mo->nindex * 2 * sizeof(double));
    }
    if (mo->refstarid) {
        dest->refstarid = malloc(mo->nindex * sizeof(int));
        memcpy(dest->refstarid, mo->refstarid, mo->nindex * sizeof(int));
    }
    if (mo->matchodds) {
        dest->matchodds = malloc(mo->nfield * sizeof(double));
        memcpy(dest->matchodds, mo->matchodds, mo->nfield * sizeof(double));
    }
    if (mo->theta) {
        dest->theta = malloc(mo->nfield * sizeof(int));
        memcpy(dest->theta, mo->theta, mo->nfield * sizeof(int));
    }
}

double verify_logodds_to_weight(double lodds) {
    if (lodds > 40.)
        return 1.0;
    if (lodds < -700)
        return 0.0;
    return exp(lodds) / (1.0 + exp(lodds));
}


double verify_star_lists(double* refxys, int NR,
                         const double* testxys, const double* testsigma2s, int NT,
                         double effective_area,
                         double distractors,
                         double logodds_bail,
                         double logodds_stoplooking,
                         int* p_besti,
                         double** p_all_logodds, int** p_theta,
                         double* p_worstlogodds,
                         int** p_testperm) {
    double X;
    verify_t v;
    double* eodds;
    int* etheta;
    int ibailed, istopped;
    int besti;
    int* theta;
    double* allodds;
    anbool score_completed;

    memset(&v, 0, sizeof(verify_t));
    v.NRall = v.NR = NR;
    v.NTall = v.NT = NT;
    // discard const here...
    v.refxy = (double*)refxys;
    v.testxy = (double*)testxys;
    v.testsigma = (double*)testsigma2s;

    v.refperm = permutation_init(NULL, NR);
    v.testperm = permutation_init(NULL, NT);

    X = verify_internal_star_lists(&v, effective_area, distractors,
                               logodds_bail, logodds_stoplooking, &besti,
                               &allodds, &theta,
                               p_worstlogodds, &ibailed, &istopped,
                               &score_completed);
    if (!score_completed) {
        if (p_all_logodds) {
            *p_all_logodds = NULL;
        }
        if (p_theta) {
            *p_theta = NULL;
        }
        if (p_testperm) {
            *p_testperm = NULL;
        }
        free(v.testperm);
        free(v.refperm);
        free(v.badguys);
        return -LARGE_VAL;
    }
    if (verify_internal_fixup_theta(
            theta, allodds, ibailed, istopped,
            &v, besti, NR, NULL,
            &etheta, &eodds)) {
        free(theta);
        free(allodds);
        free(v.testperm);
        free(v.refperm);
        free(v.badguys);
        if (p_all_logodds) {
            *p_all_logodds = NULL;
        }
        if (p_theta) {
            *p_theta = NULL;
        }
        if (p_testperm) {
            *p_testperm = NULL;
        }
        return -LARGE_VAL;
    }
    free(theta);
    free(allodds);

    if (p_all_logodds)
        *p_all_logodds = eodds;
    else
        free(eodds);
    if (p_theta)
        *p_theta = etheta;
    else
        free(etheta);

    if (p_besti)
        *p_besti = besti;

    if (p_testperm)
        *p_testperm = v.testperm;
    else
        free(v.testperm);

    free(v.refperm);
    free(v.badguys);
    return X;
}









double verify_star_lists_ror(double* refxys, int NR,
                             const double* testxys, const double* testsigma2s, int NT,
                             double pix2, double gamma,
                             const double* qc, double Q2,
                             double W, double H,
                             double distractors,
                             double logodds_bail,
                             double logodds_stoplooking,
                             int* p_besti,
                             double** p_all_logodds, int** p_theta,
                             double* p_worstlogodds,
                             int** p_testperm, int** p_refperm) {
    double X = -LARGE_VAL;
    verify_t v;
    double* eodds = NULL;
    int* etheta = NULL;
    int ibailed, istopped;
    int besti = -1;
    int* theta = NULL;
    double* allodds = NULL;
    // RoR
    double ror2;
    int igood, ibad;
    int NB = 100;
    int NBx, NBy;
    double bx0, by0;
    double stepx, stepy;
    int i, j;
    int Ngood;
    double effective_area;
    anbool score_completed;

    if (p_besti) {
        *p_besti = -1;
    }
    if (p_all_logodds) {
        *p_all_logodds = NULL;
    }
    if (p_theta) {
        *p_theta = NULL;
    }
    if (p_worstlogodds) {
        *p_worstlogodds = -LARGE_VAL;
    }
    if (p_testperm) {
        *p_testperm = NULL;
    }
    if (p_refperm) {
        *p_refperm = NULL;
    }

    memset(&v, 0, sizeof(verify_t));
    v.NRall = v.NR = NR;
    v.NTall = v.NT = NT;
    v.refxy = refxys;
    // instead of verify_internal_get_test_stars()...
    // (so we don't do:
    // --dedup
    // --remove quad stars
    // --uniformize
    // )

    // discard const here...
    v.testxy = (double*)testxys;
    v.testsigma = (double*)testsigma2s;
    v.refperm = permutation_init(NULL, NR);
    v.testperm = permutation_init(NULL, NT);
    v.tbadguys = malloc(v.NTall * sizeof(int));
    v.badguys = malloc(v.NRall * sizeof(int));

    ror2 = verify_get_ror2(Q2, W*H, distractors, NR, pix2);
    logverb("RoR: %g\n", sqrt(ror2));

    // Remove test stars outside the RoR.
    igood = ibad = 0;
    for (i=0; i<v.NT; i++) {
        int ti = v.testperm[i];
        double r2 = distsq(qc, v.testxy + 2*ti, 2);
        if (r2 < ror2) {
            v.testperm[igood] = ti;
            igood++;
        } else {
            v.tbadguys[ibad] = ti;
            ibad++;
        }
    }
    v.NT = igood;
    // remember the bad guys
    memcpy(v.testperm + igood, v.tbadguys, ibad * sizeof(int));
    logverb("Test stars in RoR: %i of %i\n", v.NT, v.NTall);

    // Count good bins to find effective area...
    NBx = ceil((double)W / sqrt(W*H) * sqrt(NB));
    NBy = ceil((double)H / sqrt(W*H) * sqrt(NB));
    NB = NBx * NBy;
    stepx = (double)W / (double)NBx;
    stepy = (double)H / (double)NBy;
    bx0 = stepx/2.0;
    by0 = stepy/2.0;
    Ngood = 0;
    for (i=0; i<NBy; i++) {
        double bxy[2];
        bxy[1] = by0 + i*stepy;
        for (j=0; j<NBx; j++) {
            double r2;
            bxy[0] = bx0 + j*stepx;
            r2 = distsq(bxy, qc, 2);
            if (r2 < ror2)
                Ngood++;
        }
    }
    effective_area = W*H * (double)Ngood / (double)NB;
    logverb("Good bins: %i / %i; effA %g of %g\n", Ngood, NB, W*H, effective_area);

    // Remove ref stars outside RoR.
    igood = ibad = 0;
    for (i=0; i<v.NR; i++) {
        int ri = v.refperm[i];
        if (distsq(qc, v.refxy + 2*ri, 2) < ror2) {
            v.refperm[igood] = ri;
            igood++;
        } else {
            v.badguys[ibad] = ri;
            ibad++;
        }
    }
    // remember the bad guys
    memcpy(v.refperm + igood, v.badguys, ibad * sizeof(int));
    v.NR = igood;
    logverb("Ref stars in RoR: %i of %i\n", v.NR, v.NRall);

    if (v.NR) {
        X = verify_internal_star_lists(&v, effective_area, distractors,
                                   logodds_bail, logodds_stoplooking, &besti,
                                   &allodds, &theta,
                                   p_worstlogodds, &ibailed, &istopped,
                                   &score_completed);
        if (!score_completed) {
            X = -LARGE_VAL;
            if (p_all_logodds) {
                *p_all_logodds = NULL;
            }
            if (p_theta) {
                *p_theta = NULL;
            }
            goto cleanup;
        }
        if (verify_internal_fixup_theta(
                theta, allodds, ibailed, istopped,
                &v, besti, NR, NULL,
                &etheta, &eodds)) {
            X = -LARGE_VAL;
            if (p_all_logodds) {
                *p_all_logodds = NULL;
            }
            if (p_theta) {
                *p_theta = NULL;
            }
            goto cleanup;
        }
        if (p_all_logodds)
            *p_all_logodds = eodds;
        else
            free(eodds);
        if (p_theta)
            *p_theta = etheta;
        else
            free(etheta);

        if (p_besti)
            *p_besti = besti;

    } else {
        X = -LARGE_VAL;
    }


cleanup:
    free(theta);
    free(allodds);
    if (p_testperm)
        *p_testperm = v.testperm;
    else
        free(v.testperm);


    if (p_refperm)
        *p_refperm = v.refperm;
    else
        free(v.refperm);

    free(v.badguys);
    free(v.tbadguys);

    return X;
}
