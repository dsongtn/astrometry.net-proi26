/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "datalog.h"
#include "log.h"
#include "mathutil.h"
#include "verify_internal.h"
#include "verify_theta_tail.h"

#define DEBUGVERIFY 0

#if DEBUGVERIFY
#define debug2(args...) logdebug(args)
#else
#define debug2(args...)
#endif

#define DATALOG_MASK_VERIFY 0x1
#define DLOG_ODDS 10
#define DLOG_ODDS_MIN log(1e6)
#define dlog(lev, fmt, ...) \
    data_log(DATALOG_MASK_VERIFY, lev, fmt, ##__VA_ARGS__)

/*
 * Verification projects a different reference-star set for every candidate.
 * Building a general-purpose KD tree for each set is expensive, particularly
 * for candidates that bail out after only a few test stars.  Use a bounded
 * flat spatial hash for the common case and retain the original KD path for
 * exact ties and pathological geometry.
 */
#define VERIFY_NN_LINEAR_WORK_LIMIT 256
#define VERIFY_NN_SMALL_KD_LIMIT 24
#define VERIFY_NN_MAX_CELL_OCCUPANCY 64
#define VERIFY_NN_MIN_SLOTS 64

typedef enum verify_nn_mode {
    VERIFY_NN_LINEAR,
    VERIFY_NN_GRID,
    VERIFY_NN_KDTREE
} verify_nn_mode_t;

typedef struct verify_nn {
    verify_nn_mode_t mode;
    double* points;
    int npoints;
    double cellsize;
    size_t nslots;
    unsigned char* workspace;
    int* heads;
    int* next;
    int64_t* cellx;
    int64_t* celly;
    kdtree_t* tree;
} verify_nn_t;

static uint64_t verify_nn_hash_word(uint64_t value) {
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static size_t verify_nn_hash_cell(int64_t x, int64_t y, size_t mask) {
    uint64_t hx = verify_nn_hash_word((uint64_t)x);
    uint64_t hy = verify_nn_hash_word((uint64_t)y);
    return (size_t)(hx ^ ((hy << 32) | (hy >> 32))) & mask;
}

static anbool verify_nn_cell_coord(double value, double cellsize,
                                  int64_t* cell) {
    double scaled;

    scaled = value / cellsize;
    /*
     * Beyond 2^52, adjacent double-precision cell coordinates are no longer
     * reliable.  The legacy KD implementation remains exact for that case.
     */
    if (!isfinite(scaled) || fabs(scaled) > 0x1p52) {
        return FALSE;
    }
    *cell = (int64_t)floor(scaled);
    return TRUE;
}

static void verify_nn_free_grid(verify_nn_t* nn) {
    free(nn->workspace);
    nn->workspace = NULL;
    nn->heads = NULL;
    nn->next = NULL;
    nn->cellx = NULL;
    nn->celly = NULL;
    nn->nslots = 0;
}

static anbool verify_nn_grid_workspace_size(int npoints,
                                            size_t* p_nslots,
                                            size_t* p_coordoffset,
                                            size_t* p_nbytes) {
    size_t coordoffset;
    size_t intcount;
    size_t nbytes;
    size_t nslots;
    size_t points;
    size_t target;

    if (npoints < 0 || !p_nslots || !p_coordoffset || !p_nbytes) {
        return FALSE;
    }
    points = (size_t)npoints;
    if (points > SIZE_MAX / 2U) {
        return FALSE;
    }
    target = points * 2U;
    nslots = VERIFY_NN_MIN_SLOTS;
    while (nslots < target) {
        if (nslots > SIZE_MAX / 2U) {
            return FALSE;
        }
        nslots *= 2U;
    }

    if (nslots > SIZE_MAX - points) {
        return FALSE;
    }
    intcount = nslots + points;
    if (intcount > SIZE_MAX / sizeof(int)) {
        return FALSE;
    }
    coordoffset = intcount * sizeof(int);
    if (coordoffset > SIZE_MAX - (sizeof(int64_t) - 1U)) {
        return FALSE;
    }
    coordoffset = (coordoffset + sizeof(int64_t) - 1U) &
                  ~(sizeof(int64_t) - 1U);
    if (nslots > (SIZE_MAX - coordoffset) /
                     (2U * sizeof(int64_t))) {
        return FALSE;
    }
    nbytes = coordoffset + 2U * nslots * sizeof(int64_t);

    *p_nslots = nslots;
    *p_coordoffset = coordoffset;
    *p_nbytes = nbytes;
    return TRUE;
}

static anbool verify_nn_build_grid(verify_nn_t* nn,
                                  const double* sigma2,
                                  const int* testperm,
                                  int ntest) {
    double maxd2 = 0.0;
    size_t coordoffset;
    size_t nbytes;
    size_t mask;
    size_t nslots;
    int maxoccupancy = 0;
    int i;
    size_t slot;

    for (i=0; i<ntest; i++) {
        double d2 = sigma2[testperm[i]] * 25.0;
        if (!isfinite(d2) || d2 <= 0.0) {
            return FALSE;
        }
        maxd2 = MAX(maxd2, d2);
    }
    if (!(maxd2 > 0.0)) {
        return FALSE;
    }
    nn->cellsize = nextafter(sqrt(maxd2), HUGE_VAL);
    if (!isfinite(nn->cellsize) || nn->cellsize <= 0.0) {
        return FALSE;
    }

    if (!verify_nn_grid_workspace_size(
            nn->npoints, &nslots, &coordoffset, &nbytes)) {
        return FALSE;
    }

    nn->workspace = calloc(1, nbytes);
    if (!nn->workspace) {
        return FALSE;
    }
    nn->nslots = nslots;
    nn->heads = (int*)nn->workspace;
    nn->next = nn->heads + nslots;
    nn->cellx = (int64_t*)(nn->workspace + coordoffset);
    nn->celly = nn->cellx + nslots;
    mask = nslots - 1U;

    for (i=0; i<nn->npoints; i++) {
        int64_t x;
        int64_t y;
        size_t slot;

        if (!verify_nn_cell_coord(nn->points[2*i], nn->cellsize, &x) ||
            !verify_nn_cell_coord(nn->points[2*i+1], nn->cellsize, &y)) {
            verify_nn_free_grid(nn);
            return FALSE;
        }
        slot = verify_nn_hash_cell(x, y, mask);
        while (nn->heads[slot] &&
               (nn->cellx[slot] != x || nn->celly[slot] != y)) {
            slot = (slot + 1U) & mask;
        }
        if (!nn->heads[slot]) {
            nn->cellx[slot] = x;
            nn->celly[slot] = y;
        }
        nn->next[i] = nn->heads[slot] - 1;
        nn->heads[slot] = i + 1;
    }

    for (slot=0; slot<nslots; slot++) {
        int occupancy = 0;
        int point;
        for (point=nn->heads[slot] - 1;
             point >= 0;
             point=nn->next[point]) {
            occupancy++;
        }
        maxoccupancy = MAX(maxoccupancy, occupancy);
    }
    if (maxoccupancy > VERIFY_NN_MAX_CELL_OCCUPANCY) {
        verify_nn_free_grid(nn);
        return FALSE;
    }
    return TRUE;
}

static anbool verify_nn_init(verify_nn_t* nn,
                             double* points,
                             int npoints,
                             const double* sigma2,
                             const int* testperm,
                             int ntest) {
    memset(nn, 0, sizeof(*nn));
    nn->mode = VERIFY_NN_LINEAR;
    nn->points = points;
    nn->npoints = npoints;
    if (npoints > 0 && ntest > 0 &&
        (size_t)npoints <=
        (size_t)VERIFY_NN_LINEAR_WORK_LIMIT / (size_t)ntest) {
        return TRUE;
    }
    if (npoints <= VERIFY_NN_SMALL_KD_LIMIT) {
        nn->tree = kdtree_build(NULL, nn->points, nn->npoints, 2, 10,
                                KDTT_DOUBLE, KD_BUILD_SPLIT);
        if (!nn->tree) {
            return FALSE;
        }
        nn->mode = VERIFY_NN_KDTREE;
        return TRUE;
    }
    if (verify_nn_build_grid(nn, sigma2, testperm, ntest)) {
        nn->mode = VERIFY_NN_GRID;
        return TRUE;
    }
    nn->tree = kdtree_build(NULL, nn->points, nn->npoints, 2, 10,
                            KDTT_DOUBLE, KD_BUILD_SPLIT);
    if (!nn->tree) {
        return FALSE;
    }
    nn->mode = VERIFY_NN_KDTREE;
    return TRUE;
}

static size_t verify_nn_find_slot(const verify_nn_t* nn,
                                  int64_t x, int64_t y) {
    size_t mask = nn->nslots - 1U;
    size_t slot = verify_nn_hash_cell(x, y, mask);

    while (nn->heads[slot]) {
        if (nn->cellx[slot] == x && nn->celly[slot] == y) {
            return slot;
        }
        slot = (slot + 1U) & mask;
    }
    return SIZE_MAX;
}

static void verify_nn_consider_point(const verify_nn_t* nn,
                                     int point,
                                     const double* query,
                                     double* bestd2,
                                     int* best,
                                     anbool* tied) {
    const double* ref = nn->points + 2*point;
    double delta;
    double d2;

    /*
     * Force a scalar rounding point between dimensions.  This matches the KD
     * kernel even when the build permits floating-point contraction.
     */
    delta = query[0] - ref[0];
    d2 = delta * delta;
    {
        volatile double rounded = d2;
        d2 = rounded;
    }
    delta = query[1] - ref[1];
    {
        volatile double squared = delta * delta;
        d2 += squared;
    }
    if (d2 > *bestd2) {
        return;
    }
    if (*best == -1 || d2 < *bestd2) {
        *best = point;
        *bestd2 = d2;
        *tied = FALSE;
    } else if (d2 == *bestd2) {
        *tied = TRUE;
    }
}

static anbool verify_nn_flat_query(const verify_nn_t* nn,
                                   const double* query,
                                   double maxd2,
                                   int* best,
                                   double* bestd2,
                                   anbool* tied) {
    int i;

    if (!isfinite(query[0]) || !isfinite(query[1]) ||
        !isfinite(maxd2) || maxd2 < 0.0) {
        return FALSE;
    }
    *best = -1;
    *bestd2 = maxd2;
    *tied = FALSE;

    if (nn->mode == VERIFY_NN_LINEAR) {
        for (i=0; i<nn->npoints; i++) {
            verify_nn_consider_point(nn, i, query, bestd2, best, tied);
        }
        return TRUE;
    }

    {
        int64_t qx;
        int64_t qy;
        int dx;
        int dy;

        if (!verify_nn_cell_coord(query[0], nn->cellsize, &qx) ||
            !verify_nn_cell_coord(query[1], nn->cellsize, &qy)) {
            return FALSE;
        }
        for (dy=-1; dy<=1; dy++) {
            for (dx=-1; dx<=1; dx++) {
                size_t slot = verify_nn_find_slot(nn, qx + dx, qy + dy);
                int point;
                if (slot == SIZE_MAX) {
                    continue;
                }
                for (point=nn->heads[slot] - 1;
                     point >= 0;
                     point=nn->next[point]) {
                    verify_nn_consider_point(nn, point, query,
                                             bestd2, best, tied);
                }
            }
        }
    }
    return TRUE;
}

static anbool verify_nn_promote(verify_nn_t* nn) {
    if (nn->mode == VERIFY_NN_KDTREE) {
        return nn->tree != NULL;
    }
    verify_nn_free_grid(nn);
    nn->tree = kdtree_build(NULL, nn->points, nn->npoints, 2, 10,
                            KDTT_DOUBLE, KD_BUILD_SPLIT);
    if (!nn->tree) {
        return FALSE;
    }
    nn->mode = VERIFY_NN_KDTREE;
    return TRUE;
}

static int verify_nn_query(verify_nn_t* nn,
                           const double* query,
                           double maxd2,
                           double* bestd2,
                           anbool* failed) {
    int best;

    if (failed) {
        *failed = FALSE;
    }
    if (nn->mode != VERIFY_NN_KDTREE) {
        anbool tied;
        if (verify_nn_flat_query(nn, query, maxd2,
                                 &best, bestd2, &tied) && !tied) {
            return best;
        }
        if (!verify_nn_promote(nn)) {
            if (failed) {
                *failed = TRUE;
            }
            return -1;
        }
    }

    if (!nn->tree) {
        if (failed) {
            *failed = TRUE;
        }
        return -1;
    }
    best = kdtree_nearest_neighbour_within(nn->tree, query, maxd2, bestd2);
    if (best == -1) {
        return -1;
    }
    return kdtree_permute(nn->tree, best);
}

static void verify_nn_cleanup(verify_nn_t* nn) {
    verify_nn_free_grid(nn);
    kdtree_free(nn->tree);
    nn->tree = NULL;
}



static double verify_logd_at(double distractor, int mu, int NR,
                             double logbg) {
    return log(distractor + (1.0-distractor)*mu / (double)NR) + logbg;
}

static double verify_logd_cached(double* values, unsigned char* ready,
                                 int count, double distractor, int mu,
                                 int NR, double logbg) {
    assert(mu >= 0);
    assert(mu < count);
    (void)count;
    if (!ready[mu]) {
        values[mu] = verify_logd_at(distractor, mu, NR, logbg);
        ready[mu] = TRUE;
    }
    return values[mu];
}

double verify_internal_star_lists(verify_t* v,
                                     double effective_area,
                                     double distractors,
                                     double logodds_bail,
                                     double logodds_stoplooking,
                                     int* p_besti,
                                     double** p_logodds, int** p_theta,
                                     double* p_worstlogodds,
                                     int* p_ibailed, int* p_istopped,
                                     anbool* p_completed) {
    int i, j;
    double worstlogodds;
    double bestworstlogodds;
    double bestlogodds;
    int besti;
    double logodds;
    double logbg;
    double logd;
    //double matchnsigma = 5.0;
    unsigned char* arrays = NULL;
    double* refcopy;
    verify_nn_t nn;
    int* rmatches;
    double* rprobs;
    double* logdcache;
    unsigned char* logdready;
    double* all_logodds = NULL;
    int* theta = NULL;
    int mu;
    int* rperm;
    size_t arrays_bytes;
    size_t double_count;
    size_t nr;
    size_t nt;
    anbool allocated_badguys = FALSE;

    memset(&nn, 0, sizeof(nn));
    if (p_completed) {
        *p_completed = FALSE;
    }
    if (p_besti) {
        *p_besti = -1;
    }
    if (p_logodds) {
        *p_logodds = NULL;
    }
    if (p_theta) {
        *p_theta = NULL;
    }
    if (p_worstlogodds) {
        *p_worstlogodds = -LARGE_VAL;
    }
    if (p_ibailed) {
        *p_ibailed = -1;
    }
    if (p_istopped) {
        *p_istopped = -1;
    }
    if (!v->NR || !v->NT) {
        logerr("real_verify_star_lists: NR=%i, NT=%i\n", v->NR, v->NT);
        return -LARGE_VAL;
    }
    if (v->NR < 0 || v->NT < 0) {
        return -LARGE_VAL;
    }
    nr = (size_t)v->NR;
    nt = (size_t)v->NT;

    /*
     * Keep the fixed-size per-candidate arrays in one allocation.  The
     * fallback KD builder may scramble refcopy, so it remains a packed copy.
     */
    if (nr > (SIZE_MAX - 1U) / 3U) {
        goto fail;
    }
    double_count = 3U * nr + 1U;
    if (nt > SIZE_MAX - double_count) {
        goto fail;
    }
    double_count += nt;
    if (double_count > SIZE_MAX / sizeof(double)) {
        goto fail;
    }
    arrays_bytes = double_count * sizeof(double);
    if (nr > (SIZE_MAX - arrays_bytes) / sizeof(int)) {
        goto fail;
    }
    arrays_bytes += nr * sizeof(int);
    if (nt == SIZE_MAX ||
        nt + 1U > SIZE_MAX - arrays_bytes) {
        goto fail;
    }
    arrays_bytes += (nt + 1U) * sizeof(unsigned char);
    arrays = malloc(arrays_bytes);
    if (!arrays) {
        goto fail;
    }
    refcopy = (double*)arrays;
    rprobs = refcopy + 2 * v->NR;
    logdcache = rprobs + v->NR;
    rmatches = (int*)(logdcache + v->NT + 1);
    logdready = (unsigned char*)(rmatches + v->NR);
    memset(logdready, 0, ((size_t)v->NT + 1U) * sizeof(unsigned char));

    // we must pack/unpermute the refxys; remember this packing order in "rperm".
    // we borrow storage for "rperm"...
    if (!v->badguys) {
        if (nr > SIZE_MAX / sizeof(int)) {
            goto fail;
        }
        v->badguys = malloc(nr * sizeof(int));
        if (!v->badguys) {
            goto fail;
        }
        allocated_badguys = TRUE;
    }
    rperm = v->badguys;
    for (i=0; i<v->NR; i++) {
        int ri = v->refperm[i];
        rperm[i] = ri;
        refcopy[2*i+0] = v->refxy[2*ri+0];
        refcopy[2*i+1] = v->refxy[2*ri+1];
    }
    if (!verify_nn_init(&nn, refcopy, v->NR,
                        v->testsigma, v->testperm, v->NT)) {
        goto fail;
    }

    for (i=0; i<v->NR; i++) {
        rmatches[i] = -1;
        rprobs[i] = -LARGE_VAL;
    }

    if (p_logodds || data_log_passes(DATALOG_MASK_VERIFY, DLOG_ODDS)) {
        if (nt > SIZE_MAX / sizeof(double)) {
            goto fail;
        }
        all_logodds = calloc(v->NT, sizeof(double));
        if (!all_logodds) {
            goto fail;
        }
    }
    if (nt > SIZE_MAX / sizeof(int)) {
        goto fail;
    }
    theta = malloc(nt * sizeof(int));
    if (!theta) {
        goto fail;
    }

    logbg = log(1.0 / effective_area);

    worstlogodds = 0;
    bestlogodds = -LARGE_VAL;
    bestworstlogodds = -LARGE_VAL;
    besti = -1;
    logodds = 0.0;
    mu = 0;
    for (i=0; i<v->NT; i++) {
        const double* testxy;
        double sig2;
        int refi;
        double d2;
        anbool query_failed;
        //double reallogfg;
        double logfg;
        int ti;

        ti = v->testperm[i];
        testxy = v->testxy + 2*ti;
        sig2 = v->testsigma[ti];

        logd = verify_logd_cached(logdcache, logdready, v->NT + 1,
                           distractors, mu, v->NR, logbg);

        debug2("\n");
        debug2("test star %i: (%.1f,%.1f), sigma: %.1f\n", i, testxy[0], testxy[1], sqrt(sig2));

        // find nearest ref star (within 5 sigma)
        refi = verify_nn_query(
            &nn, testxy, sig2 * 25.0, &d2, &query_failed);
        if (query_failed) {
            goto fail;
        }
        if (refi == -1) {
            // no nearest neighbour within range.
            debug2("  No nearest neighbour.\n");
            logfg = -LARGE_VAL;
        } else {
            double loggmax;
            // peak value of the Gaussian
            loggmax = log((1.0 - distractors) / (2.0 * M_PI * sig2 * v->NR));
            // FIXME - do something with uninformative hits?
            // these should be eliminated by RoR filtering...
            if (loggmax < logbg)
                debug2("  This star is uninformative: peak %.1f, bg %.1f.\n", loggmax, logbg);

            // value of the foreground Gaussian
            logfg = loggmax - d2 / (2.0 * sig2);

            debug2("  NN: ref star %i, dist %.2f, sigmas: %.3f, logfg: %.1f (%.1f above distractor, %.1f above bg)\n",
                   refi, sqrt(d2), sqrt(d2 / sig2), logfg, logfg - logd, logfg - logbg);
        }

        if (logfg < logd) {
            //reallogfg =
            logfg = logd;
            debug2("  Distractor.\n");
            theta[i] = THETA_DISTRACTOR;
        } else {
            // duplicate match?
            if (rmatches[refi] != -1) {
                double oldfg = rprobs[refi];
                //debug2("Conflict: odds was %g, now %g.\n", oldfg, logfg);
                // Conflict.  Compute probabilities of old vs new theta.
                // if we keep the old one: the new star is a distractor
                double keepfg = logd;

                // if we switch to the new one: the new star is a match...
                double switchfg = logfg;
                // ... and the old one becomes a distractor...
                int oldj = rmatches[refi];
                int muj = 0;
                //reallogfg = logfg;
                for (j=0; j<oldj; j++)
                    if (theta[j] >= 0)
                        muj++;
                switchfg +=
                    verify_logd_cached(logdcache, logdready, v->NT + 1,
                                distractors, muj, v->NR, logbg) - oldfg;
                // FIXME - could estimate/bound the distractor change and avoid computing it...

                // ... and the intervening distractors become worse.
                debug2("  oldj is %i, muj is %i.\n", oldj, muj);
                debug2("  changing old point to distractor: %.1f change in logodds\n",
                       (verify_logd_cached(logdcache, logdready, v->NT + 1,
                                    distractors, muj, v->NR, logbg) - oldfg));
                for (; j<i; j++)
                    if (theta[j] < 0) {
                        double current_logd =
                            verify_logd_cached(logdcache, logdready, v->NT + 1,
                                        distractors, muj, v->NR, logbg);
                        double next_logd =
                            verify_logd_cached(logdcache, logdready, v->NT + 1,
                                        distractors, muj+1, v->NR, logbg);
                        switchfg += current_logd - next_logd;
                        debug2("  adjusting distractor %i: %g change in logodds\n",
                               j, current_logd - next_logd);
                    } else
                        muj++;
                debug2("  Conflict: keeping   old match, logfg would be %.1f\n", keepfg);
                debug2("  Conflict: accepting new match, logfg would be %.1f\n", switchfg);

                if (switchfg > keepfg) {
                    // upgrade: old match becomes a distractor.
                    debug2("  Conflict: upgrading.\n");
                    theta[oldj] = THETA_CONFLICT;
                    // Note that here we want the entries in "theta" to be
                    // indices into "v->refxy" et al, so apply the "rperm" permutation.
                    theta[i] = rperm[refi];
                    // record this new match.
                    rmatches[refi] = i;
                    rprobs[refi] = logfg;

                    // "switchfg" incorporates the cost of adjusting the previous probabilities.
                    logfg = switchfg;

                    // FIXME -- Do we need to repeat the distractor-adjustment
                    // loop above, updating all_logodds entries??
                    // No, not really -- we update "logfg" in this loop, and record it below
                    // and that's sort of right -- it's THIS star that resulting in all the changes.
                    /*
                     if (all_logodds) {
                     muj = 0;
                     for (j=0; j<oldj; j++)
                     if (theta[j] >= 0)
                     muj++;
                     all_logodds[oldj] = verify_logd_at(distractors, muj, v->NR, logbg) - logbg;
                     for (j=oldj; j<i; j++)
                     if (theta[j] < 0) {
                     all_logodds[j] = verify_logd_at(distractors, muj, v->NR, logbg) - logbg;
                     } else {
                     muj++;
                     }
                     double logp = 0.;
                     for (j=0; j<i; j++)
                     logp += all_logodds[j];
                     logverb("updated all_logodds = %g, vs logodds %g\n",
                     logp, logodds);
                     }
                     */


                } else {
                    // old match was better: this match becomes a distractor.
                    debug2("  Conflict: not upgrading.\n"); //  logprob was %.1f, now %.1f.\n", oldfg, logfg);
                    logfg = keepfg;
                    theta[i] = THETA_CONFLICT;
                }
                // no change in mu.

            } else {
                // new match.
                rmatches[refi] = i;
                rprobs[refi] = logfg;
                theta[i] = rperm[refi];
                mu++;
            }
        }

        logodds += (logfg - logbg);
        debug2("  Logodds: change %.1f, now %.1f\n", (logfg - logbg), logodds);

        if (all_logodds)
            all_logodds[i] = logfg - logbg;

        if (logodds < logodds_bail) {
            debug2("  logodds %g less than bailout %g\n", logodds, logodds_bail);
            if (p_ibailed)
                *p_ibailed = i;
            break;
        }

        worstlogodds = MIN(worstlogodds, logodds);

        if (logodds > bestlogodds) {
            bestlogodds = logodds;
            besti = i;
            // Record the worst log-odds we've seen up to this point.
            bestworstlogodds = worstlogodds;
        }

        if (logodds > logodds_stoplooking) {
            if (p_istopped)
                *p_istopped = i;
            break;
        }
    }

    if (bestlogodds > DLOG_ODDS_MIN) {
        // when the loop stopped...
        int iend = i;
        data_log_start_item(DATALOG_MASK_VERIFY, DLOG_ODDS, "logodds");
        dlog(DLOG_ODDS, "[");
        for (i=0; i<iend; i++)
            dlog(DLOG_ODDS, "%s%g", (i ? ", ":""), all_logodds[i]);
        dlog(DLOG_ODDS, "]");
        data_log_end_item(DATALOG_MASK_VERIFY, DLOG_ODDS);

        data_log_start_item(DATALOG_MASK_VERIFY, DLOG_ODDS, "bestlogodds");
        dlog(DLOG_ODDS, "%g", bestlogodds);
        data_log_end_item(DATALOG_MASK_VERIFY, DLOG_ODDS);

        /*
         double lnp = 0.0;
         for (i=0; i<5; i++)
         lnp += all_logodds[i];
         if (lnp > 4.) {
         printf("lnp at step 5: %g\n", lnp);
         printf("test perm:");
         for (i=0; i<10; i++)
         printf(" %i", v->testperm[i]);
         printf("\n");
         printf("theta:");
         for (i=0; i<10; i++)
         printf(" %i", theta[i]);
         printf("\n");

         data_log_start_item(DATALOG_MASK_VERIFY, DLOG_ODDS, "match");
         dlog(DLOG_ODDS, "{ 'refxy': [");
         for (i=0; i<v->NRall; i++)
         dlog(DLOG_ODDS, "(%.3f,%.3f),", v->refxy[2*i+0], v->refxy[2*i+1]);
         dlog(DLOG_ODDS, "], 'testxy': [");
         for (i=0; i<v->NTall; i++)
         dlog(DLOG_ODDS, "(%.3f,%.3f),", v->testxy[2*i+0], v->testxy[2*i+1]);
         dlog(DLOG_ODDS, "], 'testperm': [");
         for (i=0; i<v->NT; i++)
         dlog(DLOG_ODDS, "%i,", v->testperm[i]);
         dlog(DLOG_ODDS, "], 'refperm': [");
         for (i=0; i<v->NR; i++)
         dlog(DLOG_ODDS, "%i,", v->refperm[i]);
         dlog(DLOG_ODDS, "], 'theta': [");
         for (i=0; i<v->NT; i++)
         dlog(DLOG_ODDS, "%i,", theta[i]);
         dlog(DLOG_ODDS, "], 'logodds5': %g, 'all_logodds': [", lnp);
         for (i=0; i<iend; i++)
         dlog(DLOG_ODDS, "%g,", all_logodds[i]);
         dlog(DLOG_ODDS, "] }");
         data_log_end_item(DATALOG_MASK_VERIFY, DLOG_ODDS);
         }
         */
    }

    if (p_theta) {
        *p_theta = theta;
        theta = NULL;
    } else {
        free(theta);
        theta = NULL;
    }

    if (p_besti) {
        *p_besti = besti;
    }

    if (p_worstlogodds) {
        *p_worstlogodds = bestworstlogodds;
    }

    if (p_logodds) {
        *p_logodds = all_logodds;
        all_logodds = NULL;
    } else {
        free(all_logodds);
        all_logodds = NULL;
    }

    verify_nn_cleanup(&nn);
    free(arrays);
    if (p_completed) {
        *p_completed = TRUE;
    }

    return bestlogodds;

fail:
    verify_nn_cleanup(&nn);
    free(theta);
    free(all_logodds);
    free(arrays);
    if (allocated_badguys) {
        free(v->badguys);
        v->badguys = NULL;
    }
    if (p_logodds) {
        *p_logodds = NULL;
    }
    if (p_theta) {
        *p_theta = NULL;
    }
    return -LARGE_VAL;
}


static size_t verify_score_add_bytes(size_t total,
                                     size_t count,
                                     size_t element_size) {
    size_t bytes;

    if (count > SIZE_MAX / element_size) {
        return SIZE_MAX;
    }
    bytes = count * element_size;
    if (total > SIZE_MAX - bytes) {
        return SIZE_MAX;
    }
    return total + bytes;
}

static size_t verify_nn_kdtree_workspace_bytes(int npoints) {
    size_t bottom = 1U;
    size_t interior;
    size_t quotient;
    size_t total;

    if (npoints <= 0) {
        return 0U;
    }
    quotient = (size_t)npoints / 10U;
    while (quotient) {
        if (bottom > SIZE_MAX / 2U) {
            return SIZE_MAX;
        }
        bottom *= 2U;
        quotient >>= 1U;
    }
    interior = bottom - 1U;

    total = sizeof(kdtree_t);
    total = verify_score_add_bytes(
        total, (size_t)npoints, sizeof(u32));
    total = verify_score_add_bytes(
        total, bottom, sizeof(int32_t));
    total = verify_score_add_bytes(
        total, interior, sizeof(double) + sizeof(u8));
    return total;
}

size_t verify_internal_score_workspace_bytes(int npoints) {
    size_t coordoffset;
    size_t grid_bytes;
    size_t kd_bytes;
    size_t slots;

    if (!verify_nn_grid_workspace_size(
            npoints, &slots, &coordoffset, &grid_bytes)) {
        return SIZE_MAX;
    }
    kd_bytes = verify_nn_kdtree_workspace_bytes(npoints);
    if (kd_bytes == SIZE_MAX) {
        return SIZE_MAX;
    }
    (void)slots;
    (void)coordoffset;
    return MAX(grid_bytes, kd_bytes);
}
