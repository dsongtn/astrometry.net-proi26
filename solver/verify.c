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
#include "index_shard_internal.h"

#define DEBUGVERIFY 0

#define VERIFY_PROJECTION_MAX_TASKS 6U
#define VERIFY_PROJECTION_MIN_STARS_PER_TASK 256U

/*
 * verify.o is also linked into legacy tools that do not include
 * index_shard.o. Weak references keep those tools on the original inline
 * path while the solver engine supplies the strong helper scheduler.
 */
#if defined(__GNUC__) && !defined(_WIN32)
extern size_t index_shard_helper_available_workers(void)
    __attribute__((weak));
extern index_shard_helper_run_status_t index_shard_helper_run(
    const index_shard_helper_ops_t* ops,
    index_shard_helper_task_t* tasks,
    size_t task_count,
    index_shard_helper_run_stats_t* stats) __attribute__((weak));
extern anbool index_shard_worker_stop_requested(void)
    __attribute__((weak));
#define VERIFY_PROJECTION_HELPERS_LINKED 1
#define VERIFY_STOP_CHECK_LINKED 1
#else
#define VERIFY_PROJECTION_HELPERS_LINKED 0
#define VERIFY_STOP_CHECK_LINKED 0
#endif

#if DEBUGVERIFY
#define debug2(args...) logdebug(args)
#else
#define debug2(args...)
#endif

#define DATALOG_MASK_VERIFY 0x1

// level
#define DLOG_ODDS 10

#define DLOG_ODDS_MIN log(1e6)

#define dlog(lev, fmt, ...) data_log(DATALOG_MASK_VERIFY, lev, fmt, ##__VA_ARGS__)

anbool verify_datalog_enabled(void) {
    return data_log_passes(DATALOG_MASK_VERIFY, DLOG_ODDS);
}

struct verify_s {
    const sip_t* wcs;

    // Reference stars:
    int NR;
    int NRall;
    int* refperm;
    int* refstarid;
    double* refxy;
    // temp storage when filtering
    int* badguys;

    // Image stars:
    int NT;
    int NTall;
    int* testperm;
    double* testxy;
    double* testsigma; // actually sigma**2.
    // temp storage
    int* tbadguys;

};
typedef struct verify_s verify_t;

typedef enum verify_prepared_state {
    VERIFY_PREPARED_READY = 0,
    VERIFY_PREPARED_NO_REFERENCE = 1,
    VERIFY_PREPARED_NO_QUAD_REFERENCE = 2,
    VERIFY_PREPARED_NO_ROR_REFERENCE = 3,
    VERIFY_PREPARED_EMPTY_LISTS = 4
} verify_prepared_state_t;

struct verify_index_query {
    const startree_t* source;
    double center[3];
    double radius2;
    double* refxyz;
    int* refstarid;
    uint8_t* sweep;
    int nrall;
};

struct verify_prepared_hit {
    verify_t verify;
    sip_t wcs;
    double* refxyz;
    double effective_area;
    double distractors;
    double logbail;
    double logaccept;
    double logstoplooking;
    int nrimage;
    verify_prepared_state_t state;
    anbool fake_match;
};

int verify_query_hit(const startree_t* skdt,
                     const double center[3],
                     double radius2,
                     verify_index_query_t** query) {
    verify_index_query_t* result;

    if (!query) {
        return -1;
    }
    *query = NULL;
    if (!skdt || !center || !skdt->tree || !skdt->sweep) {
        return -1;
    }
    result = calloc(1, sizeof(*result));
    if (!result) {
        return -1;
    }
    result->source = skdt;
    memcpy(result->center, center, sizeof(result->center));
    result->radius2 = radius2;
    startree_search_for(skdt, center, radius2,
                        &result->refxyz, NULL,
                        &result->refstarid, &result->nrall);
    if (result->nrall < 0 ||
        (result->nrall &&
         (!result->refxyz || !result->refstarid)) ||
        (!result->nrall &&
         (result->refxyz || result->refstarid))) {
        verify_destroy_index_query(result);
        return -1;
    }
    *query = result;
    return 0;
}

size_t verify_index_query_count(const verify_index_query_t* query) {
    if (!query || query->nrall <= 0) {
        return 0U;
    }
    return (size_t)query->nrall;
}

size_t verify_index_query_bytes(const verify_index_query_t* query) {
    size_t count;
    size_t total;

    if (!query || query->nrall < 0) {
        return 0U;
    }
    count = (size_t)query->nrall;
    total = sizeof(*query);
    if (count > (SIZE_MAX - total) / (3U * sizeof(double))) {
        return SIZE_MAX;
    }
    total += count * 3U * sizeof(double);
    if (count > (SIZE_MAX - total) / sizeof(int)) {
        return SIZE_MAX;
    }
    total += count * sizeof(int);
    if (count > (SIZE_MAX - total) / sizeof(*query->sweep)) {
        return SIZE_MAX;
    }
    total += count * sizeof(*query->sweep);
    return total;
}

int verify_index_query_sweep_range(const startree_t* skdt,
                                   const verify_index_query_t* query,
                                   size_t index,
                                   const void** data,
                                   size_t* size) {
    int starid;

    if (!data || !size) {
        return -1;
    }
    *data = NULL;
    *size = 0U;
    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        index >= (size_t)query->nrall ||
        !query->refstarid) {
        return -1;
    }
    starid = query->refstarid[index];
    if (starid < 0 || starid >= startree_N(skdt)) {
        return -1;
    }
    *data = skdt->sweep + starid;
    *size = sizeof(*skdt->sweep);
    return 0;
}

int verify_index_query_capture_sweep(const startree_t* skdt,
                                     verify_index_query_t* query) {
    uint8_t* sweep;
    size_t count;
    int nstars;
    int i;

    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        (query->nrall && !query->refstarid)) {
        return -1;
    }
    if (query->sweep || !query->nrall) {
        return 0;
    }
    count = (size_t)query->nrall;
    nstars = startree_N(skdt);
    for (i = 0; i < query->nrall; i++) {
        if (query->refstarid[i] < 0 ||
            query->refstarid[i] >= nstars) {
            return -1;
        }
    }
    sweep = malloc(count * sizeof(*sweep));
    if (!sweep) {
        return -1;
    }
    for (i = 0; i < query->nrall; i++) {
#if VERIFY_STOP_CHECK_LINKED
        if (!(i & 255) &&
            index_shard_worker_stop_requested &&
            index_shard_worker_stop_requested()) {
            free(sweep);
            return -1;
        }
#endif
        sweep[i] = skdt->sweep[query->refstarid[i]];
    }
    query->sweep = sweep;
    return 0;
}

static int verify_mapped_page_buffers_validate(
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count) {
    uintptr_t previous_end = 0U;
    size_t i;

    if (!buffers || !buffer_count) {
        return -1;
    }
    for (i = 0U; i < buffer_count; i++) {
        uintptr_t begin;
        uintptr_t end;

        if (!buffers[i].mapping_data || !buffers[i].bytes ||
            !buffers[i].size) {
            return -1;
        }
        begin = (uintptr_t)buffers[i].mapping_data;
        if (buffers[i].size > UINTPTR_MAX - begin) {
            return -1;
        }
        end = begin + buffers[i].size;
        if (i && begin < previous_end) {
            return -1;
        }
        previous_end = end;
    }
    return 0;
}

static int verify_mapped_page_buffer_find(
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count,
    uintptr_t target,
    size_t target_size,
    const unsigned char** bytes) {
    size_t low = 0U;
    size_t high = buffer_count;

    if (!bytes || !target_size) {
        return -1;
    }
    *bytes = NULL;
    while (low < high) {
        size_t middle = low + (high - low) / 2U;
        uintptr_t begin = (uintptr_t)buffers[middle].mapping_data;
        uintptr_t end = begin + buffers[middle].size;

        if (target < begin) {
            high = middle;
            continue;
        }
        if (target >= end) {
            low = middle + 1U;
            continue;
        }
        if (target_size > end - target) {
            return -1;
        }
        *bytes = buffers[middle].bytes + (size_t)(target - begin);
        return 0;
    }
    return -1;
}

int verify_index_query_capture_sweep_buffers(
    const startree_t* skdt,
    verify_index_query_t* query,
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count) {
    uint8_t* sweep;
    size_t count;
    int nstars;
    int i;

    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        (query->nrall && !query->refstarid)) {
        return -1;
    }
    if (query->sweep || !query->nrall) {
        return 0;
    }
    if (verify_mapped_page_buffers_validate(
            buffers, buffer_count)) {
        return -1;
    }
    count = (size_t)query->nrall;
    if (count > SIZE_MAX / sizeof(*sweep)) {
        return -1;
    }
    sweep = malloc(count * sizeof(*sweep));
    if (!sweep) {
        return -1;
    }
    nstars = startree_N(skdt);
    for (i = 0; i < query->nrall; i++) {
        const unsigned char* bytes;
        uintptr_t target;
        int starid = query->refstarid[i];

#if VERIFY_STOP_CHECK_LINKED
        if (!(i & 255) &&
            index_shard_worker_stop_requested &&
            index_shard_worker_stop_requested()) {
            free(sweep);
            return -1;
        }
#endif
        if (starid < 0 || starid >= nstars) {
            free(sweep);
            return -1;
        }
        target = (uintptr_t)(skdt->sweep + starid);
        if (verify_mapped_page_buffer_find(
                buffers, buffer_count,
                target, sizeof(*sweep), &bytes)) {
            free(sweep);
            return -1;
        }
        memcpy(&sweep[i], bytes, sizeof(*sweep));
    }
    query->sweep = sweep;
    return 0;
}

void verify_destroy_index_query(verify_index_query_t* query) {
    if (!query) {
        return;
    }
    free(query->refxyz);
    free(query->refstarid);
    free(query->sweep);
    memset(query, 0, sizeof(*query));
    free(query);
}

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

static double logd_at(double distractor, int mu, int NR, double logbg) {
    return log(distractor + (1.0-distractor)*mu / (double)NR) + logbg;
}

static double logd_cached(double* values, unsigned char* ready, int count,
                          double distractor, int mu, int NR, double logbg) {
    assert(mu >= 0);
    assert(mu < count);
    (void)count;
    if (!ready[mu]) {
        values[mu] = logd_at(distractor, mu, NR, logbg);
        ready[mu] = TRUE;
    }
    return values[mu];
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

static int verify_get_test_stars(
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

static int verify_apply_ror(verify_t* v,
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

    if (verify_get_test_stars(
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

static double real_verify_star_lists(verify_t* v,
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

        logd = logd_cached(logdcache, logdready, v->NT + 1,
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
                    logd_cached(logdcache, logdready, v->NT + 1,
                                distractors, muj, v->NR, logbg) - oldfg;
                // FIXME - could estimate/bound the distractor change and avoid computing it...

                // ... and the intervening distractors become worse.
                debug2("  oldj is %i, muj is %i.\n", oldj, muj);
                debug2("  changing old point to distractor: %.1f change in logodds\n",
                       (logd_cached(logdcache, logdready, v->NT + 1,
                                    distractors, muj, v->NR, logbg) - oldfg));
                for (; j<i; j++)
                    if (theta[j] < 0) {
                        double current_logd =
                            logd_cached(logdcache, logdready, v->NT + 1,
                                        distractors, muj, v->NR, logbg);
                        double next_logd =
                            logd_cached(logdcache, logdready, v->NT + 1,
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
                     all_logodds[oldj] = logd_at(distractors, muj, v->NR, logbg) - logbg;
                     for (j=oldj; j<i; j++)
                     if (theta[j] < 0) {
                     all_logodds[j] = logd_at(distractors, muj, v->NR, logbg) - logbg;
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

typedef struct verify_projection_context {
    const double* xyz;
    int npoints;
    int width;
    int height;
    anbool use_sip;
    union {
        sip_t sip;
        tan_t tan;
    } wcs;
} verify_projection_context_t;

typedef struct verify_projection_task_input {
    const verify_projection_context_t* context;
    size_t first;
    size_t count;
} verify_projection_task_input_t;

typedef struct verify_projection_task_output {
    double x;
    double y;
    anbool inside;
} verify_projection_task_output_t;

static index_shard_helper_task_status_t
verify_projection_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const verify_projection_task_input_t* input = input_bytes;
    const verify_projection_context_t* context;
    verify_projection_task_output_t* output = output_bytes;
    size_t point;

    if (!input || input_size != sizeof(*input) ||
        !output ||
        input->count > SIZE_MAX / sizeof(*output) ||
        output_size != input->count * sizeof(*output)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    context = input->context;
    if (!context || !context->xyz || context->npoints < 0 ||
        input->first > (size_t)context->npoints ||
        input->count >
            (size_t)context->npoints - input->first) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    for (point = 0U; point < input->count; point++) {
        const double* xyz = context->xyz +
            3U * (input->first + point);
        double x = 0.0;
        double y = 0.0;
        anbool projected;

        if (context->use_sip) {
            projected = sip_xyzarr2pixelxy(
                &context->wcs.sip, xyz, &x, &y);
        } else {
            projected = tan_xyzarr2pixelxy(
                &context->wcs.tan, xyz, &x, &y);
        }
        output[point].x = x;
        output[point].y = y;
        output[point].inside =
            projected &&
            x >= 0.0 && y >= 0.0 &&
            x < context->width && y < context->height;
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

static const index_shard_helper_ops_t
verify_projection_helper_ops = {
    "verify-projection",
    verify_projection_helper_execute
};

/*
 * Project only the already-copied startree_search_for() result. The helper
 * package is index-free and immutable; stable owner compaction preserves the
 * exact input ordering of sip_filter_stars_in_field().
 */
static anbool verify_filter_stars_in_field_parallel(
    const sip_t* sip,
    const tan_t* tan,
    const double* xyz,
    int npoints,
    double** p_xy,
    int** p_inbounds,
    int* p_ngood) {
#if VERIFY_PROJECTION_HELPERS_LINKED
    verify_projection_context_t context;
    verify_projection_task_input_t
        inputs[VERIFY_PROJECTION_MAX_TASKS];
    index_shard_helper_task_t
        tasks[VERIFY_PROJECTION_MAX_TASKS];
    verify_projection_task_output_t* projected = NULL;
    int* inbounds = NULL;
    double* xy = NULL;
    size_t available;
    size_t task_count;
    size_t base;
    size_t remainder;
    size_t cursor = 0U;
    size_t task_index;
    int ngood = 0;
    int point;
    index_shard_helper_run_status_t run_status;

    if (!index_shard_helper_available_workers ||
        !index_shard_helper_run ||
        (!sip && !tan) || !xyz || npoints <= 0 ||
        !p_inbounds || !p_ngood ||
        (size_t)npoints <
            2U * VERIFY_PROJECTION_MIN_STARS_PER_TASK) {
        return FALSE;
    }
    available = index_shard_helper_available_workers();
    if (!available) {
        return FALSE;
    }
    task_count = (size_t)npoints /
        VERIFY_PROJECTION_MIN_STARS_PER_TASK;
    if (task_count > VERIFY_PROJECTION_MAX_TASKS) {
        task_count = VERIFY_PROJECTION_MAX_TASKS;
    }
    if (available < task_count - 1U) {
        task_count = available + 1U;
    }
    if (task_count < 2U ||
        (size_t)npoints > SIZE_MAX / sizeof(*projected) ||
        (size_t)npoints > SIZE_MAX / sizeof(*inbounds) ||
        (p_xy &&
         (size_t)npoints > SIZE_MAX / (2U * sizeof(*xy)))) {
        return FALSE;
    }

    projected = malloc((size_t)npoints * sizeof(*projected));
    inbounds = malloc((size_t)npoints * sizeof(*inbounds));
    if (p_xy) {
        xy = malloc((size_t)npoints * 2U * sizeof(*xy));
    }
    if (!projected || !inbounds || (p_xy && !xy)) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    memset(&context, 0, sizeof(context));
    context.xyz = xyz;
    context.npoints = npoints;
    context.use_sip = sip != NULL;
    if (sip) {
        context.wcs.sip = *sip;
        context.width = sip->wcstan.imagew;
        context.height = sip->wcstan.imageh;
    } else {
        context.wcs.tan = *tan;
        context.width = tan->imagew;
        context.height = tan->imageh;
    }

    memset(inputs, 0, sizeof(inputs));
    memset(tasks, 0, sizeof(tasks));
    base = (size_t)npoints / task_count;
    remainder = (size_t)npoints % task_count;
    for (task_index = 0U;
         task_index < task_count;
         task_index++) {
        size_t count = base +
            (task_index < remainder ? 1U : 0U);

        inputs[task_index].context = &context;
        inputs[task_index].first = cursor;
        inputs[task_index].count = count;
        tasks[task_index].input = &inputs[task_index];
        tasks[task_index].input_bytes = sizeof(inputs[task_index]);
        tasks[task_index].output = projected + cursor;
        tasks[task_index].output_bytes =
            count * sizeof(*projected);
        tasks[task_index].work_units =
            (unsigned long long)count;
        cursor += count;
    }
    if (cursor != (size_t)npoints) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    run_status = index_shard_helper_run(
        &verify_projection_helper_ops,
        tasks,
        task_count,
        NULL);
    if (run_status != INDEX_SHARD_HELPER_OK) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    for (point = 0; point < npoints; point++) {
        if (!projected[point].inside) {
            continue;
        }
        inbounds[ngood] = point;
        if (xy) {
            xy[2 * ngood] = projected[point].x;
            xy[2 * ngood + 1] = projected[point].y;
        }
        ngood++;
    }
    free(projected);

    if (!ngood) {
        free(xy);
        free(inbounds);
        xy = NULL;
        inbounds = NULL;
    }
    if (p_xy) {
        *p_xy = xy;
    }
    *p_inbounds = inbounds;
    *p_ngood = ngood;
    return TRUE;
#else
    (void)sip;
    (void)tan;
    (void)xyz;
    (void)npoints;
    (void)p_xy;
    (void)p_inbounds;
    (void)p_ngood;
    return FALSE;
#endif
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
    if (!verify_filter_stars_in_field_parallel(
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


static void set_null_mo(MatchObj* mo) {
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

static int fixup_theta(int* theta, double* allodds,
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

    // "theta" has length v->NT.

    if (ibailed != -1)
        for (i=ibailed+1; i<v->NT; i++)
            theta[i] = THETA_BAILEDOUT;

    if (istopped != -1)
        for (i=istopped+1; i<v->NT; i++)
            theta[i] = THETA_STOPPEDLOOKING;

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
        if (verify_apply_ror(v, index_cutnside, mo,
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
        if (verify_get_test_stars(
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
    K = real_verify_star_lists(v, effA, distractors,
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

        if (fixup_theta(
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
    set_null_mo(mo);
    // uh oh, spaghetti-code-oh!
    goto cleanup;
}

int verify_prepare_hit_from_query(const startree_t* skdt,
                                  verify_index_query_t** query,
                                  int index_cutnside,
                                  const MatchObj* mo, const sip_t* sip,
                                  const verify_field_t* vf,
                                  double pix2, double distractors,
                                  double fieldW, double fieldH,
                                  double logbail, double logaccept,
                                  double logstoplooking,
                                  anbool do_gamma, anbool fake_match,
                                  verify_prepared_hit_t** prepared) {
    verify_index_query_t* query_context;
    verify_prepared_hit_t* context;
    verify_t* v;
    double fieldr2;
    int* sweep = NULL;
    int nstars;
    int i;
    int j;
    int ibad;
    int igood;

    if (!prepared) {
        return -1;
    }
    *prepared = NULL;
    if (!query || !*query || !skdt || !mo || !vf ||
        (!mo->wcs_valid && !sip) ||
        !isfinite(logaccept) || !isfinite(logbail)) {
        return -1;
    }
    query_context = *query;
    fieldr2 = square(mo->radius);
    if (!skdt->tree ||
        query_context->source != skdt ||
        memcmp(query_context->center, mo->center,
               sizeof(query_context->center)) ||
        memcmp(&query_context->radius2, &fieldr2, sizeof(fieldr2)) ||
        query_context->nrall < 0 ||
        (query_context->nrall &&
         (!query_context->refxyz || !query_context->refstarid)) ||
        (!query_context->nrall &&
         (query_context->refxyz || query_context->refstarid)) ||
        (query_context->nrall && !query_context->sweep &&
         !skdt->sweep)) {
        return -1;
    }
    context = calloc(1, sizeof(*context));
    if (!context) {
        return -1;
    }
    context->distractors = distractors;
    context->logbail = logbail;
    context->logaccept = logaccept;
    context->logstoplooking = logstoplooking;
    context->fake_match = fake_match;
    context->state = VERIFY_PREPARED_READY;
    v = &context->verify;

    if (sip) {
        memcpy(&context->wcs, sip, sizeof(context->wcs));
    } else {
        sip_wrap_tan(&mo->wcstan, &context->wcs);
    }
    v->wcs = &context->wcs;
    context->refxyz = query_context->refxyz;
    v->refstarid = query_context->refstarid;
    v->NRall = query_context->nrall;
    if (!context->refxyz) {
        context->state = VERIFY_PREPARED_NO_REFERENCE;
        goto done;
    }

    if ((size_t)v->NRall > SIZE_MAX / (2U * sizeof(double)) ||
        (size_t)v->NRall > SIZE_MAX / sizeof(int)) {
        goto fail;
    }
    v->refxy = malloc((size_t)v->NRall * 2U * sizeof(double));
    v->refperm = malloc((size_t)v->NRall * sizeof(int));
    if (!v->refxy || !v->refperm) {
        goto fail;
    }
    igood = 0;
    for (i = 0; i < v->NRall; i++) {
        if (!sip_xyzarr2pixelxy(v->wcs,
                                context->refxyz + 3 * i,
                                v->refxy + 2 * i,
                                v->refxy + 2 * i + 1) ||
            !sip_pixel_is_inside_image(v->wcs,
                                       v->refxy[2 * i],
                                       v->refxy[2 * i + 1])) {
            continue;
        }
        v->refperm[igood++] = i;
    }
    v->NR = igood;
    context->nrimage = v->NR;

    sweep = malloc((size_t)v->NRall * sizeof(int));
    if (!sweep) {
        goto fail;
    }
    nstars = startree_N(skdt);
    for (i = 0; i < v->NRall; i++) {
        int starid = v->refstarid[i];

        if (starid < 0 || starid >= nstars) {
            goto fail;
        }
        if (query_context->sweep) {
            sweep[i] = query_context->sweep[i];
        } else {
            sweep[i] = skdt->sweep[starid];
        }
    }
    permuted_sort(sweep, sizeof(int), compare_ints_asc,
                  v->refperm, v->NR);
    free(sweep);
    sweep = NULL;

    if (v->NR) {
        v->badguys = malloc((size_t)v->NR * sizeof(int));
        if (!v->badguys) {
            goto fail;
        }
    }
    if (!fake_match) {
        ibad = 0;
        igood = 0;
        for (i = 0; i < v->NR; i++) {
            anbool inquad = FALSE;
            int ri = v->refperm[i];

            for (j = 0; j < mo->dimquads; j++) {
                if (v->refstarid[ri] == (int)mo->star[j]) {
                    inquad = TRUE;
                    v->badguys[ibad++] = ri;
                    break;
                }
            }
            if (!inquad) {
                v->refperm[igood++] = ri;
            }
        }
        if (ibad) {
            memcpy(v->refperm + igood, v->badguys,
                   (size_t)ibad * sizeof(int));
        }
        v->NR = igood;
    }
    if (!v->NR) {
        context->state = VERIFY_PREPARED_NO_QUAD_REFERENCE;
        goto done;
    }

    if (!fake_match) {
        if (verify_apply_ror(
                v, index_cutnside, (MatchObj*)mo,
                vf, pix2, distractors, fieldW, fieldH,
                do_gamma, fake_match,
                &context->effective_area, NULL, NULL)) {
            goto fail;
        }
        if (!v->NR) {
            context->state = VERIFY_PREPARED_NO_ROR_REFERENCE;
            goto done;
        }
    } else {
        if (verify_get_test_stars(
                v, vf, (MatchObj*)mo,
                pix2, do_gamma, fake_match)) {
            goto fail;
        }
        context->effective_area = fieldW * fieldH;
    }
    if (!v->NR || !v->NT) {
        context->state = VERIFY_PREPARED_EMPTY_LISTS;
        goto done;
    }

done:
    free(v->badguys);
    v->badguys = NULL;
    free(v->tbadguys);
    v->tbadguys = NULL;
    query_context->refxyz = NULL;
    query_context->refstarid = NULL;
    verify_destroy_index_query(query_context);
    *query = NULL;
    *prepared = context;
    return 0;

fail:
    free(sweep);
    context->refxyz = NULL;
    v->refstarid = NULL;
    verify_destroy_prepared_hit(context);
    return -1;
}

int verify_prepare_hit(const startree_t* skdt, int index_cutnside,
                       const MatchObj* mo, const sip_t* sip,
                       const verify_field_t* vf,
                       double pix2, double distractors,
                       double fieldW, double fieldH,
                       double logbail, double logaccept,
                       double logstoplooking,
                       anbool do_gamma, anbool fake_match,
                       verify_prepared_hit_t** prepared) {
    verify_index_query_t* query = NULL;
    double fieldr2;
    int status;

    if (!prepared || !skdt || !mo || !vf ||
        (!mo->wcs_valid && !sip) ||
        !isfinite(logaccept) || !isfinite(logbail)) {
        return -1;
    }
    *prepared = NULL;
    fieldr2 = square(mo->radius);
    if (verify_query_hit(skdt, mo->center, fieldr2, &query)) {
        return -1;
    }
    status = verify_prepare_hit_from_query(
        skdt, &query, index_cutnside, mo, sip, vf,
        pix2, distractors, fieldW, fieldH,
        logbail, logaccept, logstoplooking,
        do_gamma, fake_match, prepared);
    verify_destroy_index_query(query);
    return status;
}

int verify_score_prepared_hit(const verify_prepared_hit_t* prepared,
                              verify_prepared_score_t* score) {
    verify_t local;
    anbool score_completed;

    if (!prepared || !score || score->theta || score->allodds ||
        score->complete) {
        return -1;
    }
    memset(score, 0, sizeof(*score));
    score->logodds = -LARGE_VAL;
    score->worstlogodds = -LARGE_VAL;
    score->besti = -1;
    score->ibailed = -1;
    score->istopped = -1;
    if (prepared->state != VERIFY_PREPARED_READY) {
        score->complete = TRUE;
        return 0;
    }

    local = prepared->verify;
    if (local.NR <= 0 ||
        (size_t)local.NR > SIZE_MAX / sizeof(int)) {
        return -1;
    }
    local.badguys = malloc((size_t)local.NR * sizeof(int));
    if (!local.badguys) {
        return -1;
    }
    local.tbadguys = NULL;
    score->logodds = real_verify_star_lists(
        &local,
        prepared->effective_area,
        prepared->distractors,
        prepared->logbail,
        prepared->logstoplooking,
        &score->besti,
        &score->allodds,
        &score->theta,
        &score->worstlogodds,
        &score->ibailed,
        &score->istopped,
        &score_completed);
    free(local.badguys);
    if (!score_completed) {
        verify_destroy_prepared_score(score);
        return -1;
    }
    score->complete = TRUE;
    return 0;
}

int verify_finish_prepared_hit(verify_prepared_hit_t* prepared,
                               verify_prepared_score_t* score,
                               MatchObj* mo) {
    verify_t* v;
    double* refxyz;
    double K;
    int besti;
    int i;
    int j;

    if (!prepared || !score || !score->complete || !mo) {
        return -1;
    }
    switch (prepared->state) {
    case VERIFY_PREPARED_NO_REFERENCE:
        logverb("No reference stars in the bounding circle\n");
        set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_NO_QUAD_REFERENCE:
        logverb("After removing quad stars: no reference stars\n");
        set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_NO_ROR_REFERENCE:
        logerr("After applying ROR, NR = 0!\n");
        set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_EMPTY_LISTS:
        logverb("After applying RoR, NR=%i, NT=%i\n",
                prepared->verify.NR, prepared->verify.NT);
        set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_READY:
        break;
    default:
        return -1;
    }

    v = &prepared->verify;
    refxyz = prepared->refxyz;
    K = score->logodds;
    besti = score->besti;

    if (log_get_level() >= LOG_ALL) {
        int nm;
        int nc;
        int nd;

        verify_count_hits(score->theta, besti, &nm, &nc, &nd);
        debug("verify: logodds %g, %i matches, %i conflicts, "
              "%i distractors after %i field objects.\n",
              K, nm, nc, nd, besti);
    }

    if (K >= prepared->logaccept) {
        int* etheta;
        double* eodds;
        int nm;
        int nc;
        int nd;

        verify_count_hits(score->theta, besti, &nm, &nc, &nd);
        if (fixup_theta(
                score->theta, score->allodds,
                score->ibailed, score->istopped,
                v, besti, prepared->nrimage, refxyz,
                &etheta, &eodds)) {
            return -1;
        }
        mo->logodds = K;
        mo->worstlogodds = score->worstlogodds;
        mo->nfield = v->NTall;
        mo->nindex = prepared->nrimage;
        mo->nmatch = nm;
        mo->nconflict = nc;
        mo->ndistractor = nd;

        if (!prepared->fake_match) {
            for (j = 0; j < mo->dimquads; j++) {
                for (i = 0; i < prepared->nrimage; i++) {
                    if (v->refstarid[i] == (int)mo->star[j]) {
                        int ti = mo->field[j];

                        assert(etheta[ti] == THETA_FILTERED);
                        etheta[ti] = i;
                        eodds[ti] = LARGE_VAL;
                        break;
                    }
                }
            }
        }

        mo->theta = etheta;
        mo->matchodds = eodds;
        mo->refxyz = prepared->refxyz;
        prepared->refxyz = NULL;
        mo->refxy = v->refxy;
        v->refxy = NULL;
        mo->refstarid = v->refstarid;
        v->refstarid = NULL;
        mo->testperm = v->testperm;
        v->testperm = NULL;
        matchobj_compute_derived(mo);
    } else {
        mo->logodds = K;
        mo->worstlogodds = score->worstlogodds;
        mo->nfield = v->NTall;
        mo->nindex = prepared->nrimage;
    }
    verify_destroy_prepared_score(score);
    return 0;
}

static size_t verify_prepared_add_bytes(size_t total,
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
    total = verify_prepared_add_bytes(
        total, (size_t)npoints, sizeof(u32));
    total = verify_prepared_add_bytes(
        total, bottom, sizeof(int32_t));
    total = verify_prepared_add_bytes(
        total, interior, sizeof(double) + sizeof(u8));
    return total;
}

size_t verify_prepared_hit_bytes(const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t total;

    if (!prepared) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NRall < 0 || v->NTall < 0) {
        return SIZE_MAX;
    }
    total = sizeof(*prepared);
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, 3U * sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, 2U * sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, sizeof(int));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, sizeof(int));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NTall, sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NTall, sizeof(int));
    return total;
}

size_t verify_prepared_score_bytes(
    const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t total = 0U;

    if (!prepared) {
        return 0U;
    }
    if (prepared->state != VERIFY_PREPARED_READY) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NT < 0) {
        return SIZE_MAX;
    }
    total = verify_prepared_add_bytes(
        total, (size_t)v->NT, sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NT, sizeof(int));
    return total;
}

size_t verify_prepared_hit_peak_bytes(
    const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t finish_peak = 0U;
    size_t grid_bytes;
    size_t grid_coordoffset;
    size_t grid_slots;
    size_t kd_bytes;
    size_t nn_bytes;
    size_t retained;
    size_t score_peak = 0U;
    size_t transient_peak;

    if (!prepared) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NR < 0 || v->NRall < 0 ||
        v->NT < 0 || v->NTall < 0 ||
        prepared->nrimage < 0) {
        return SIZE_MAX;
    }
    retained = verify_prepared_hit_bytes(prepared);
    if (retained == SIZE_MAX) {
        return SIZE_MAX;
    }
    switch (prepared->state) {
    case VERIFY_PREPARED_NO_REFERENCE:
    case VERIFY_PREPARED_NO_QUAD_REFERENCE:
    case VERIFY_PREPARED_NO_ROR_REFERENCE:
    case VERIFY_PREPARED_EMPTY_LISTS:
        return retained;
    case VERIFY_PREPARED_READY:
        break;
    default:
        return SIZE_MAX;
    }
    if (!verify_nn_grid_workspace_size(
            v->NR, &grid_slots, &grid_coordoffset, &grid_bytes)) {
        return SIZE_MAX;
    }
    kd_bytes = verify_nn_kdtree_workspace_bytes(v->NR);
    if (kd_bytes == SIZE_MAX) {
        return SIZE_MAX;
    }
    nn_bytes = MAX(grid_bytes, kd_bytes);
    (void)grid_slots;
    (void)grid_coordoffset;

    /*
     * Scoring retains one conflict array, the packed nearest-neighbor
     * inputs, the result vectors, and at most one nearest-neighbor
     * workspace. Grid and legacy KD storage never coexist, so use the
     * larger exact payload bound.
     */
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NR, sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NR,
        3U * sizeof(double) + sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NT,
        sizeof(double) + sizeof(unsigned char));
    score_peak = verify_prepared_add_bytes(
        score_peak, 1U, sizeof(double) + sizeof(unsigned char));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NT,
        sizeof(double) + sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, nn_bytes, 1U);

    /*
     * Owner retirement keeps the score vectors while fixup_theta() builds
     * its expanded result and permutation workspaces. These allocations do
     * not overlap the scoring workspace, so admission uses the larger peak.
     */
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NT,
        sizeof(double) + sizeof(int));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NTall,
        sizeof(int) + sizeof(double));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NRall, sizeof(int));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)prepared->nrimage,
        3U * sizeof(double));
    transient_peak = MAX(score_peak, finish_peak);
    return verify_prepared_add_bytes(retained, transient_peak, 1U);
}

unsigned long long
verify_prepared_hit_work_units(const verify_prepared_hit_t* prepared) {
    unsigned long long nr;
    unsigned long long nt;

    if (!prepared || prepared->state != VERIFY_PREPARED_READY) {
        return 0U;
    }
    nr = (unsigned long long)prepared->verify.NR;
    nt = (unsigned long long)prepared->verify.NT;
    if (nr && nt > ULLONG_MAX / nr) {
        return ULLONG_MAX;
    }
    return nr * nt;
}

void verify_destroy_prepared_score(verify_prepared_score_t* score) {
    if (!score) {
        return;
    }
    free(score->allodds);
    free(score->theta);
    memset(score, 0, sizeof(*score));
}

void verify_destroy_prepared_hit(verify_prepared_hit_t* prepared) {
    verify_t* v;

    if (!prepared) {
        return;
    }
    v = &prepared->verify;
    free(prepared->refxyz);
    free(v->testperm);
    free(v->testsigma);
    free(v->tbadguys);
    free(v->refperm);
    free(v->refxy);
    free(v->refstarid);
    free(v->badguys);
    memset(prepared, 0, sizeof(*prepared));
    free(prepared);
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

    X = real_verify_star_lists(&v, effective_area, distractors,
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
    if (fixup_theta(
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
    double X;
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

    memset(&v, 0, sizeof(verify_t));
    v.NRall = v.NR = NR;
    v.NTall = v.NT = NT;
    v.refxy = refxys;
    // instead of verify_get_test_stars()...
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
        X = real_verify_star_lists(&v, effective_area, distractors,
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
        if (fixup_theta(
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
