/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>


#include "os-features.h"
#include "kdtree.h"
#include "kdtree_internal.h"
#include "kdtree_prefetch_internal.h"
#include "kdtree_continuation_internal.h"
#include "kdtree_mem.h"
#include "keywords.h"
#include "errors.h"
#include "mathutil.h"

#define KDTREE_MAX_RESULTS 1000
#define KDTREE_MAX_DIM 100

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define KDTREE_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#define KDTREE_THREAD_LOCAL __thread
#else
#error "A thread-local storage implementation is required"
#endif

#ifndef KDTREE_PREFETCH_STACK_MAX
#define KDTREE_PREFETCH_STACK_MAX 4096
#endif


#define WARNING(x, ...) fprintf(stderr, x, ## __VA_ARGS__)

#define MANGLE(x) KDMANGLE(x, ETYPE, DTYPE, TTYPE)

/*
 The "external" type is the data type that the outside world works in.

 The "data" type is the type in which we store the points.

 The "tree" type is the type in which we store the bounding boxes or splitting planes.

 Recall that    etype >= dtype >= ttype:
 .   (etype >= dtype: because there's no point storing more precision than needed by the
 .                    outside world;
 .    dtype >= ttype: because there's no point keeping more precision in the splitting
 .                    plane than exists in the data.)


 The following will be defined:

 - etype:           typedef, the "external" type.
 - ETYPE_INTEGER:   1 if the "external" space is integral (not floating-point)
 - ETYPE_MIN
 - ETYPE_MAX:       the limits of the "external" space.
 - ETYPE:           the "external" type, in a form that the preprocessor can make C
 .                  identifiers out of it: eg 'd' for double; use this to get to the
 .                  particular data type in the kdtree unions:
 .         kdtree_qres_t* query_results = ...;
 .         etype* mydata = query_results.results.ETYPE;

 - dtype:           typedef, the "data" type.
 - DTYPE_INTEGER:   1 if the "data" space is integral (not floating-point)
 - DTYPE_DOUBLE:    1 if the "data" space has type 'double'.
 - DTYPE_MIN
 - DTYPE_MAX:       the limits of the "data" space.
 - DTYPE_KDT_DATA:  eg KDT_DATA_DOUBLE
 - DTYPE:           the "data" type, in a form that the preprocessor can make C
 .                  identifiers out of it: eg 'd' for double:
 .          dtype* mydata = kd->data.DTYPE;

 - ttype:           typedef, the "tree" type.
 - TTYPE_INTEGER:   1 if the "tree" space is integral (not floating-point)
 - TTYPE_MIN
 - TTYPE_MAX:       the limits of the "tree" space.
 - TTYPE_SQRT_MAX:  the square root of the maximum value of the "tree" space.
 - TTYPE:           the "tree" type, in a form that the preprocessor can make C
 .                  identifiers out of it: eg 'd' for double.
 .          ttype* mybb = kd->bb.TTYPE;

 - bigttype:        typedef, a type that can hold a "ttype" squared.
 - BIGTTYPE:        #define the bigttype; used for "STRINGIFY" macro.
 - BIGTTYPE_MAX:    maximum value of a 'bigttype'.


 - EQUAL_ED:   1 if "external" and "data" spaces are the same.
 - EQUAL_DT:   1 if "data" and "tree" spaces are the same.
 - EQUAL_ET:   1 if "external" and "tree" spaces are the same; implies EQUAL_ED && EQUAL_DT.


 - POINT_ED(kd, d, c, func):  d: dimension;
 .                            c: coordinate in that dimension,
 .                            func: function to apply (may be empty)
 .          Converts a coordinate from "external" to "data" space.
 .          eg      POINT_ED(kd, 0, 0.1, floor);
 .                  POINT_ED(kd, 1, 0.0, );

 - POINT_DT(kd, d, c, func):
 .          Converts a coordinate from "data" to "tree" space.

 - POINT_ET(kd, d, c, func):
 .          Converts a coordinate from "external" to "tree" space.

 - POINT_TD(kd, d, c):
 .          Converts a coordinate from "tree" to "data" space.

 - POINT_DE(kd, d, c):
 .          Converts a coordinate from "data" to "external" space.

 - POINT_TE(kd, d, c):
 .          Converts a coordinate from "tree" to "external" space.


 - DIST_ED(kd, dist, func):
 - DIST2_ED(kd, dist2, func):
 .          Convert distance or distance-squared from "external" to "data" space,
 .          optionally applying a function to the result.

 - DIST_DT(kd, dist, func):
 - DIST2_DT(kd, dist2, func):
 .          Convert distance or distance-squared from "data" to "tree" space.

 - DIST_ET(kd, dist, func):
 - DIST2_ET(kd, dist2, func):
 .          Convert distance or distance-squared from "external" to "tree" space.

 - DIST_TD(kd, dist):
 - DIST2_TD(kd, dist2):
 .          Convert distance or distance-squared from "tree" to "data" space.

 - DIST_DE(kd, dist):
 - DIST2_DE(kd, dist2):
 .          Convert distance or distance-squared from "data" to "external" space.

 - DIST_TE(kd, dist):
 - DIST2_TE(kd, dist2):
 .          Convert distance or distance-squared from "tree" to "external" space.
 */

// Which function do we use for rounding?
#define KD_ROUND rint

// Get the low corner of the bounding box
#define LOW_HR( kd, D, i) ((kd)->bb.TTYPE + (2*(size_t)(i)*(size_t)(D)))

// Get the high corner of the bounding box
#define HIGH_HR(kd, D, i) ((kd)->bb.TTYPE + ((2*(size_t)(i)+1)*(size_t)(D)))

// Get the splitting-plane position
#define KD_SPLIT(kd, i) ((kd)->split.TTYPE + (size_t)(i))

// Get a pointer to the 'i'-th data point.
#define KD_DATA(kd, D, i) ((kd)->data.DTYPE + ((size_t)(D)*(size_t)(i)))

// Dereference an array of dimension D at element i, dimension d.
#define KD_ARRAY_VAL(arr, D, i, d) ((arr)[(((size_t)(D)*(size_t)(i)) + (size_t)(d))])

#define KD_ARRAY_REF(arr, D, i, d) ((arr) + (((size_t)(D)*(size_t)(i)) + (size_t)(d)))

// Get the 'i'-th element of the permutation vector, or 'i' if there is no permutation vector.
#define KD_PERM(kd, i) ((kd)->perm ? (kd)->perm[i] : i)

// Get the dimension of this tree.
#if defined(KD_DIM)
#define DIMENSION(kd)   (KD_DIM)
#else
#define DIMENSION(kd)   (kd->ndim)
#endif

// Get the size of a single point in the tree.
#define SIZEOF_PT(kd)  (sizeof(dtype)*DIMENSION(kd))

// compatibility macros (for DEPRECATED old-fashioned trees)
#define COMPAT_NODE_SIZE(kd)    (sizeof(kdtree_node_t) + (SIZEOF_PT(kd) * 2))
#define COMPAT_HIGH_HR(kd, i)  ((ttype*)(((char*)(kd)->nodes)           \
                                         + COMPAT_NODE_SIZE(kd)*(i)     \
                                         + sizeof(kdtree_node_t)        \
                                         + SIZEOF_PT(kd)))
#define COMPAT_LOW_HR(kd, i)   ((ttype*) (((char*)(kd)->nodes)          \
                                          + COMPAT_NODE_SIZE(kd)*(i)    \
                                          + sizeof(kdtree_node_t)))

// "converted" <-> "dtype" conversion.
#define POINT_CTOR(kd, d, val)  (((val) - ((kd)->minval[(d)])) * (kd)->scale)
#define POINT_RTOC(kd, d, val)  (((val) * ((kd)->invscale)) + (kd)->minval[d])

#define DIST_CTOR(kd, dist)       ((dist) * (kd)->scale)
#define DIST2_CTOR(kd, dist2)     ((dist2) * (kd)->scale * (kd)->scale)
#define DIST2_RTOC(kd, dist2)     ((dist2) * (kd)->invscale * (kd)->invscale)


// this must be at least as big as the biggest integer TTYPE.
typedef u32 bigint;

#define ACTUAL_STRINGIFY(x) (#x)
#define STRINGIFY(x) ACTUAL_STRINGIFY(x)

#define MYGLUE2(a, b) a ## b
#define DIST_FUNC_MANGLE(x, suffix) MYGLUE2(x, suffix)

/*
 static void split_dim_and_value(kdtree_t* kd, int node,
 uint8_t* splitdim, ttype* splitval) {
 if (kd->splitdim) {
 *splitdim = kd->splitdim[node];
 *splitval = *KD_SPLIT(kd, node);
 } else {
 bigint tmpsplit = *KD_SPLIT(kd, node);
 *splitdim = (uint8_t)(tmpsplit & kd->dimmask);
 *splitval = (ttype)(tmpsplit & kd->splitmask);
 }
 }
 */

/* min/maxdist functions. */
#define CAN_OVERFLOW 0
#undef  DELTAMAX

#define PTYPE etype
#define DISTTYPE double
#define FUNC_SUFFIX
#include "kdtree_internal_dists.c"
#undef PTYPE
#undef DISTTYPE
#undef FUNC_SUFFIX

#undef CAN_OVERFLOW
#define CAN_OVERFLOW 1

#define PTYPE ttype
#define DISTTYPE ttype
#define FUNC_SUFFIX _ttype
#define DELTAMAX TTYPE_SQRT_MAX
#include "kdtree_internal_dists.c"
#undef PTYPE
#undef DISTTYPE
#undef FUNC_SUFFIX
#undef DELTAMAX

#define PTYPE ttype
#define DISTTYPE bigttype
#define FUNC_SUFFIX _bigttype
#undef  DELTAMAX
#include "kdtree_internal_dists.c"
#undef PTYPE
#undef DISTTYPE
#undef FUNC_SUFFIX

#undef CAN_OVERFLOW


void MANGLE(kdtree_update_funcs)(kdtree_t* kd);

                                static anbool bboxes(const kdtree_t* kd, int node,
                                                     ttype** p_tlo, ttype** p_thi, int D) {
                                    if (kd->bb.any) {
                                        // bb trees
                                        *p_tlo =  LOW_HR(kd, D, node);
                                        *p_thi = HIGH_HR(kd, D, node);
                                        return TRUE;
                                    } else {
                                        return FALSE;
                                    }
                                }



static inline double dist2(const kdtree_t* kd, const etype* q, const dtype* p, int D) {
    int d;
    double d2 = 0.0;
#if defined(KD_DIM)
    D = KD_DIM;
#endif
    for (d=0; d<D; d++) {
        etype pp = POINT_DE(kd, d, p[d]);
        double delta;
        if (TTYPE_INTEGER) {
            if (q[d] > pp)
                delta = q[d] - pp;
            else
                delta = pp - q[d];
        } else {
            delta = q[d] - pp;
        }
        d2 += delta * delta;
    }
    return d2;
}

static inline void dist2_bailout(const kdtree_t *kd,
                                 const etype *q,
                                 const dtype *p,
                                 int D,
                                 double maxd2,
                                 anbool *bailedout,
                                 double *d2res) {
  int d;
  double d2 = 0.0;

#if defined(KD_DIM)
  D = KD_DIM;
#endif

  /*
   * Phase C: specialized 4D CodeKD distance kernel.
   *
   * The dominant solver CodeKD workload is four-dimensional.  Keep the
   * original POINT_DE() decoding, accumulation order, bailout predicate and
   * result semantics exactly unchanged while removing the generic loop
   * control from this hot path.
   */
  if (D == 4) {
    double delta;
    etype pp;

    pp = POINT_DE(kd, 0, p[0]);
    delta = q[0] - pp;
    d2 += delta * delta;

    if (d2 > maxd2) {
      *bailedout = TRUE;
      return;
    }

    pp = POINT_DE(kd, 1, p[1]);
    delta = q[1] - pp;
    d2 += delta * delta;

    if (d2 > maxd2) {
      *bailedout = TRUE;
      return;
    }

    pp = POINT_DE(kd, 2, p[2]);
    delta = q[2] - pp;
    d2 += delta * delta;

    if (d2 > maxd2) {
      *bailedout = TRUE;
      return;
    }

    pp = POINT_DE(kd, 3, p[3]);
    delta = q[3] - pp;
    d2 += delta * delta;

    if (d2 > maxd2) {
      *bailedout = TRUE;
      return;
    }

    *d2res = d2;
    return;
  }

  /*
   * Original generic implementation for every non-4D specialization.
   */
  for (d = 0; d < D; d++) {
    double delta;
    etype pp = POINT_DE(kd, d, p[d]);

    delta = q[d] - pp;
    d2 += delta * delta;

    if (d2 > maxd2) {
      *bailedout = TRUE;
      return;
    }
  }

  *d2res = d2;
}
#if defined(KDTREE_CODEKD_DSS_U16_FAST_PATH)

/*
 * Specialized four-dimensional CodeKD leaf-distance kernel for the
 * double/U16/U16 kdtree representation.
 *
 * Numerical invariants:
 *   - POINT_DE() remains the authoritative stored-to-external conversion.
 *   - Dimensions are evaluated in original order: 0, 1, 2, 3.
 *   - Squared terms are accumulated sequentially.
 *   - Bailout remains strictly d2 > maxd2.
 *   - d2res is written only for an accepted point.
 */
static inline anbool MANGLE(codekd_dist2_4d_dss_bailout)
     (const kdtree_t *kd,
      const etype *query,
      const dtype *point,
      double maxd2,
      double *d2res) {
  double delta;
  double d2;

  delta = (double)query[0] -
      (double)POINT_DE(kd, 0, point[0]);

  d2 = delta * delta;

  if (d2 > maxd2) {
    return TRUE;
  }

  delta = (double)query[1] -
      (double)POINT_DE(kd, 1, point[1]);

  d2 += delta * delta;

  if (d2 > maxd2) {
    return TRUE;
  }

  delta = (double)query[2] -
      (double)POINT_DE(kd, 2, point[2]);

  d2 += delta * delta;

  if (d2 > maxd2) {
    return TRUE;
  }

  delta = (double)query[3] -
      (double)POINT_DE(kd, 3, point[3]);

  d2 += delta * delta;

  if (d2 > maxd2) {
    return TRUE;
  }

  *d2res = d2;
  return FALSE;
}

#endif
static inline void ddist2_bailout(const kdtree_t* kd,
                                  const dtype* q, const dtype* p,
                                  int D, bigttype maxd2, anbool* bailedout,
                                  bigttype* d2res) {
    int d;
    bigttype d2 = 0;
#if defined(KD_DIM)
    D = KD_DIM;
#endif
    for (d=0; d<D; d++) {
        dtype delta;
        if (q[d] > p[d])
            delta = q[d] - p[d];
        else
            delta = p[d] - q[d];

        d2 += (bigttype)delta * (bigttype)delta;
        if (d2 > maxd2) {
            *bailedout = TRUE;
            return;
        }
    }
    *d2res = d2;
}

//ANCHOR - dist2 bailout
static inline anbool MANGLE(codekd_dist2_2d_bailout)
     (const kdtree_t *kd,
      const etype *q,
      const dtype *p,
      double maxd2,
      double *d2res) {
  double delta0;
  double delta1;
  double d2;

  /*
   * Specialized CodeKD leaf test for 2D code trees.
   *
   * This preserves the same external-coordinate distance semantics as
   * dist2_bailout(): stored dtype values are converted through POINT_DE().
   */
  delta0 = (double)q[0] - (double)POINT_DE(kd, 0, p[0]);
  d2 = delta0 * delta0;
  if (d2 > maxd2)
    return TRUE;

  delta1 = (double)q[1] - (double)POINT_DE(kd, 1, p[1]);
  d2 += delta1 * delta1;
  if (d2 > maxd2)
    return TRUE;

  *d2res = d2;
  return FALSE;
}

static inline anbool MANGLE(codekd_can_use_2d_leaf_kernel)
     (const kdtree_t *kd,
      int D,
      int options,
      anbool do_dists) {
  if (!kd)
    return FALSE;

  if (D != 2)
    return FALSE;

  if (!do_dists)
    return FALSE;

  if (!(options & KD_OPTIONS_COMPUTE_DISTS))
    return FALSE;

  if (!(options & KD_OPTIONS_SMALL_RADIUS))
    return FALSE;

  if (!(options & KD_OPTIONS_USE_SPLIT))
    return FALSE;

  if (options & KD_OPTIONS_SORT_DISTS)
    return FALSE;

  return TRUE;
}


static inline anbool dist2_exceeds(const kdtree_t* kd, const etype* q, const dtype* p, int D, double maxd2) {
    int d;
    double d2 = 0.0;
#if defined(KD_DIM)
    D = KD_DIM;
#endif
    for (d=0; d<D; d++) {
        double delta;
        double pp = POINT_DE(kd, d, p[d]);
        if (TTYPE_INTEGER) {
            if (q[d] > pp)
                delta = q[d] - pp;
            else
                delta = pp - q[d];
        } else {
            delta = q[d] - pp;
        }
        d2 += delta * delta;
        if (d2 > maxd2)
            return 1;
    }
    return 0;
}

static anbool bb_point_l1mindist_exceeds_ttype(ttype* lo, ttype* hi,
                                               ttype* query, int D,
                                               ttype maxl1, ttype maxlinf) {
    ttype dist = 0;
    ttype newdist;
    ttype delta;
    int d;
#if defined(KD_DIM)
    D = KD_DIM;
#endif
    for (d=0; d<D; d++) {
        if (query[d] < lo[d])
            delta = lo[d] - query[d];
        else if (query[d] > hi[d])
            delta = query[d] - hi[d];
        else
            continue;
        if (delta > maxlinf) {
            //printf("linf: %u > %u.\n", (unsigned int)delta, (unsigned int)maxlinf);
            return TRUE;
        }
        newdist = dist + delta;
        if (newdist < dist) {
            // overflow!
            return TRUE;
        }
        if (newdist > maxl1) {
            //printf("l1: %u > %u\n", (unsigned int)newdist, (unsigned int)maxl1);
            return TRUE;
        }
        dist = newdist;
    }
    return FALSE;
}

static void compute_splitbits(kdtree_t* kd) {
    int D;
    int bits;
    u32 val;
    D = kd->ndim;
    bits = 0;
    val = 1;
    while (val < D) {
        bits++;
        val *= 2;
    }
    kd->dimmask = val - 1;
    kd->dimbits = bits;
    kd->splitmask = ~kd->dimmask;
}

/* Sorts results by kq->sdists */
static int kdtree_qsort_results(kdtree_qres_t *kq, int D) {
    int beg[KDTREE_MAX_RESULTS], end[KDTREE_MAX_RESULTS], i = 0, j, L, R;
    etype piv_vec[KDTREE_MAX_DIM];
    anbool do_points = (kq->results.any != NULL);
    unsigned int piv_perm;
    double piv;

    beg[0] = 0;
    end[0] = kq->nres - 1;
    while (i >= 0) {
        L = beg[i];
        R = end[i];
        if (L < R) {
            piv = kq->sdists[L];
            if (do_points) {
                for (j=0; j<D; j++)
                    piv_vec[j] = kq->results.ETYPE[L*D + j];
            }
            piv_perm = kq->inds[L];
            if (i == KDTREE_MAX_RESULTS - 1) /* Sanity */
                assert(0);
            while (L < R) {
                while (kq->sdists[R] >= piv && L < R)
                    R--;
                if (L < R) {
                    if (do_points) {
                        for (j=0; j<D; j++)
                            kq->results.ETYPE[L*D + j] = kq->results.ETYPE[R*D + j];
                    }
                    kq->inds  [L] = kq->inds  [R];
                    kq->sdists[L] = kq->sdists[R];
                    L++;
                }
                while (kq->sdists[L] <= piv && L < R)
                    L++;
                if (L < R) {
                    if (do_points) {
                        for (j=0; j<D; j++)
                            kq->results.ETYPE[R*D + j] = kq->results.ETYPE[L*D + j];
                    }
                    kq->inds  [R] = kq->inds  [L];
                    kq->sdists[R] = kq->sdists[L];
                    R--;
                }
            }
            if (do_points) {
                for (j=0; j<D; j++)
                    kq->results.ETYPE[D*L + j] = piv_vec[j];
            }
            kq->inds  [L] = piv_perm;
            kq->sdists[L] = piv;
            beg[i + 1] = L + 1;
            end[i + 1] = end[i];
            end[i++] = L;
        } else
            i--;
    }
    return 1;
}

static void print_results(kdtree_qres_t* res, int D) {
    if (TRUE) {
        int i, d;
        printf("%i results.\n", res->nres);
        for (i=0; i<res->nres; i++) {
            printf("  ind %i", res->inds[i]);
            if (res->sdists)
                printf(", dist %g", res->sdists[i]);
            if (res->results.any) {
                printf(", pt [ ");
                for (d=0; d<D; d++)
                    printf("%g ", (double)res->results.ETYPE[i*D + d]);
                printf("]");
            }
            printf("\n");
        }
        printf("\n");
    }
}

static
anbool resize_results(kdtree_qres_t* res, int newsize, int D,
                      anbool do_dists, anbool do_points) {

    if (FALSE) {
        printf("resize results: before:\n");
        print_results(res, D);
    }

    if (do_dists)
        res->sdists  = REALLOC(res->sdists , newsize * sizeof(double));
    if (do_points) {
        res->results.any = REALLOC(res->results.any, (size_t)newsize * (size_t)D * sizeof(etype));
    } else {
        FREE(res->results.any);
        res->results.any = NULL;
    }
    res->inds = REALLOC(res->inds, newsize * sizeof(u32));
    if (newsize && ((do_points && !res->results.any) ||
                    (do_dists && !res->sdists) || !res->inds))
        SYSERROR("Failed to resize kdtree results arrays");
    res->capacity = newsize;

    if (FALSE) {
        printf("resize results: after:\n");
        print_results(res, D);
    }

    return TRUE;
}

static
anbool add_result(const kdtree_t* kd, kdtree_qres_t* res, double sdist,
                  unsigned int ind, const dtype* pt,
                  int D, anbool do_dists, anbool do_points) {

    if (FALSE) {
        printf("Before adding new result:\n");
        print_results(res, D);
    }

    if (do_dists)
        res->sdists[res->nres] = sdist;
    res->inds  [res->nres] = ind;
    if (do_points) {
        int d;
        for (d=0; d<D; d++)
            res->results.ETYPE[res->nres * D + d] = POINT_DE(kd, d, pt[d]);
    }
    res->nres++;
    if (res->nres == res->capacity) {
        // enlarge arrays.
        return resize_results(res, res->capacity * 2, D, do_dists, do_points);
    }

    if (FALSE) {
        printf("After adding new result:\n");
        print_results(res, D);
    }

    return TRUE;
}

/*
 Can the query be represented as a ttype?

 If so, place the converted value in "tquery".
 */
static anbool ttype_query(const kdtree_t* kd, const etype* query, ttype* tquery) {
    etype val;
    int d, D=kd->ndim;
    for (d=0; d<D; d++) {
        val = POINT_ET(kd, d, query[d], );
        if (val < TTYPE_MIN || val > TTYPE_MAX)
            return FALSE;
        tquery[d] = (ttype)val;
    }
    return TRUE;
}



double MANGLE(kdtree_get_splitval)(const kdtree_t* kd, int nodeid) {
    Unused int dim;
    ttype split = *KD_SPLIT(kd, nodeid);
    if (EQUAL_ET) {
        return split;
    }
    if (!kd->splitdim && TTYPE_INTEGER) {
        bigint tmpsplit = split;
        dim = tmpsplit & kd->dimmask;
        return POINT_TE(kd, dim, tmpsplit & kd->splitmask);
    } else {
        dim = kd->splitdim[nodeid];
    }
    return POINT_TE(kd, dim, split);
}


static void kdtree_nn_bb(const kdtree_t* kd, const etype* query,
                         double* p_bestd2, int* p_ibest) {
    int nodestack[100];
    double dist2stack[100];
    int stackpos = 0;
    int D = (kd ? kd->ndim : 0);
    anbool use_tquery = FALSE;
    anbool use_tmath = FALSE;
    anbool use_bigtmath = FALSE;
    ttype tquery[D];
    double bestd2 = *p_bestd2;
    int ibest = *p_ibest;
    ttype tl2 = 0;
    bigttype bigtl2 = 0;

#if defined(KD_DIM)
    assert(kd->ndim == KD_DIM);
    D = KD_DIM;
#else
    D = kd->ndim;
#endif

    if (TTYPE_INTEGER) {
        use_tquery = ttype_query(kd, query, tquery);
    }
    if (TTYPE_INTEGER && use_tquery) {
        double dtl2 = DIST2_ET(kd, bestd2, );
        if (dtl2 < TTYPE_MAX) {
            use_tmath = TRUE;
        } else if (dtl2 < BIGTTYPE_MAX) {
            use_bigtmath = TRUE;
        }
        bigtl2 = ceil(dtl2);
        tl2    = bigtl2;
    }

    // queue root.
    nodestack[0] = 0;
    dist2stack[0] = 0.0;
    if (kd->fun.nn_enqueue)
        kd->fun.nn_enqueue(kd, 0, 1);

    while (stackpos >= 0) {
        int nodeid;
        int i;
        int L, R;
        ttype *tlo=NULL, *thi=NULL;
        int child;
        double childd2[2];
        double firstd2, secondd2;
        int firstid, secondid;

        if (dist2stack[stackpos] > bestd2) {
            // pruned!
            if (kd->fun.nn_prune)
                kd->fun.nn_prune(kd, nodestack[stackpos], dist2stack[stackpos], bestd2, 1);
            stackpos--;
            continue;
        }
        nodeid = nodestack[stackpos];
        stackpos--;
        if (kd->fun.nn_explore)
            kd->fun.nn_explore(kd, nodeid, dist2stack[stackpos+1], bestd2);

        if (KD_IS_LEAF(kd, nodeid)) {
            // Back when leaf nodes didn't have BBoxes:
            //|| KD_IS_LEAF(kd, KD_CHILD_LEFT(nodeid)))
            dtype* data;
            L = kdtree_left(kd, nodeid);
            R = kdtree_right(kd, nodeid);
            for (i=L; i<=R; i++) {
                anbool bailedout = FALSE;
                double dsqd;
                if (kd->fun.nn_point)
                    kd->fun.nn_point(kd, nodeid, i);
                data = KD_DATA(kd, D, i);
                dist2_bailout(kd, query, data, D, bestd2, &bailedout, &dsqd);
                if (bailedout)
                    continue;
                // new best
                ibest = i;
                bestd2 = dsqd;
                if (kd->fun.nn_new_best)
                    kd->fun.nn_new_best(kd, nodeid, i, bestd2);
            }
            continue;
        }

        childd2[0] = childd2[1] = LARGE_VAL;
        for (child=0; child<2; child++) {
            anbool bailed;
            double dist2;
            int childid = (child ? KD_CHILD_RIGHT(nodeid) : KD_CHILD_LEFT(nodeid));

            bboxes(kd, childid, &tlo, &thi, D);

            bailed = FALSE;
            if (TTYPE_INTEGER && use_tmath) {
                ttype newd2 = 0;
                bb_point_mindist2_bailout_ttype(tlo, thi, tquery, D, tl2, &bailed, &newd2);
                if (bailed) {
                    if (kd->fun.nn_prune)
                        kd->fun.nn_prune(kd, nodeid, newd2, bestd2, 2);
                    continue;
                }
                dist2 = DIST2_TE(kd, newd2);
            } else if (TTYPE_INTEGER && use_bigtmath) {
                bigttype newd2 = 0;
                bb_point_mindist2_bailout_bigttype(tlo, thi, tquery, D, bigtl2, &bailed, &newd2);
                if (bailed) {
                    if (kd->fun.nn_prune)
                        kd->fun.nn_prune(kd, nodeid, newd2, bestd2, 3);
                    continue;
                }
                dist2 = DIST2_TE(kd, newd2);
            } else {
                etype bblo, bbhi;
                int d;
                // this is just bb_point_mindist2_bailout...
                dist2 = 0.0;
                for (d=0; d<D; d++) {
                    bblo = POINT_TE(kd, d, tlo[d]);
                    if (query[d] < bblo) {
                        dist2 += (bblo - query[d])*(bblo - query[d]);
                    } else {
                        bbhi = POINT_TE(kd, d, thi[d]);
                        if (query[d] > bbhi) {
                            dist2 += (query[d] - bbhi)*(query[d] - bbhi);
                        } else
                            continue;
                    }
                    if (dist2 > bestd2) {
                        bailed = TRUE;
                        break;
                    }
                }
                if (bailed) {
                    if (kd->fun.nn_prune)
                        kd->fun.nn_prune(kd, childid, dist2, bestd2, 4);
                    continue;
                }
            }
            childd2[child] = dist2;
        }

        if (childd2[0] <= childd2[1]) {
            firstd2 = childd2[0];
            secondd2 = childd2[1];
            firstid = KD_CHILD_LEFT(nodeid);
            secondid = KD_CHILD_RIGHT(nodeid);
        } else {
            firstd2 = childd2[1];
            secondd2 = childd2[0];
            firstid = KD_CHILD_RIGHT(nodeid);
            secondid = KD_CHILD_LEFT(nodeid);
        }

        if (firstd2 == LARGE_VAL)
            continue;

        // it's a stack, so put the "second" one on first.
        if (secondd2 != LARGE_VAL) {
            stackpos++;
            nodestack[stackpos] = secondid;
            dist2stack[stackpos] = secondd2;
            if (kd->fun.nn_enqueue)
                kd->fun.nn_enqueue(kd, secondid, 2);
        }

        stackpos++;
        nodestack[stackpos] = firstid;
        dist2stack[stackpos] = firstd2;
        if (kd->fun.nn_enqueue)
            kd->fun.nn_enqueue(kd, firstid, 2);

    }
    *p_bestd2 = bestd2;
    *p_ibest = ibest;
}

static void kdtree_nn_int_split(const kdtree_t* kd, const etype* query,
                                const ttype* tquery,
                                double* p_bestd2, int* p_ibest) {
    int nodestack[100];
    ttype mindists[100];

    int stackpos = 0;
    int D = kd->ndim;

    ttype closest_so_far;
    bigttype closest2;

    int ibest = -1;

    dtype* data;
    dtype* dquery = (dtype*)tquery;

    /** FIXME **/
    assert(sizeof(dtype) == sizeof(ttype));

    {
        double closest;
        closest = DIST_ET(kd, sqrt(*p_bestd2), );
        if (closest > TTYPE_MAX) {
            closest_so_far = TTYPE_MAX;
            closest2 = BIGTTYPE_MAX;
        } else {
            closest_so_far = ceil(closest);
            closest2 = (bigttype)closest_so_far * (bigttype)closest_so_far;
        }
    }

    // queue root.
    nodestack[0] = 0;
    mindists[0] = 0;

    while (stackpos >= 0) {
        int nodeid;
        int i;
        int dim = -1;
        int L, R;
        ttype split = 0;

        if (mindists[stackpos] > closest_so_far) {
            // pruned!
            stackpos--;
            continue;
        }
        nodeid = nodestack[stackpos];
        stackpos--;

        if (KD_IS_LEAF(kd, nodeid)) {
            int oldbest = ibest;

            L = kdtree_left(kd, nodeid);
            R = kdtree_right(kd, nodeid);
            for (i=L; i<=R; i++) {
                anbool bailedout = FALSE;
                bigttype dsqd;
                data = KD_DATA(kd, D, i);
                ddist2_bailout(kd, dquery, data, D, closest2, &bailedout, &dsqd);
                if (bailedout)
                    continue;
                // new best
                ibest = i;
                closest2 = dsqd;
            }

            if (oldbest != ibest) {
                // FIXME - replace with int sqrt
                closest_so_far = ceil(sqrt((double)closest2));
            }
            continue;
        }

        // split/dim trees
        split = *KD_SPLIT(kd, nodeid);

        if (kd->splitdim)
            dim = kd->splitdim[nodeid];
        else {
            bigint tmpsplit;
            tmpsplit = split;
            dim = tmpsplit & kd->dimmask;
            split = tmpsplit & kd->splitmask;
        }


        if (tquery[dim] < split) {
            // query is on the "left" side of the split.
            assert(query[dim] < POINT_TE(kd, dim, split));
            // is the right child within range?
            // look mum, no int overflow!
            if (split - tquery[dim] <= closest_so_far) {
                // visit right child - it is within range.
                assert(POINT_TE(kd, dim, split) - query[dim] > 0.0);
                //assert(POINT_TE(kd, dim, split) - query[dim] <= bestdist);
                stackpos++;
                nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
                mindists[stackpos] = split - tquery[dim];
            }
            stackpos++;
            nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
            mindists[stackpos] = 0;

        } else {
            // query is on "right" side.
            assert(POINT_TE(kd, dim, split) <= query[dim]);
            // is the left child within range?
            if (tquery[dim] - split < closest_so_far) {
                assert(query[dim] - POINT_TE(kd, dim, split) >= 0.0);
                //assert(query[dim] - POINT_TE(kd, dim, split) < bestdist);
                stackpos++;
                nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
                mindists[stackpos] = tquery[dim] - split;
            }
            stackpos++;

            nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
            mindists[stackpos] = 0;
        }
    }
    if (ibest != -1) {
        //*p_bestd2 = DIST2_TE(kd, closest2);
        // Recompute the d2 more precisely in "etype":
        data = KD_DATA(kd, D, ibest);
        *p_bestd2 = dist2(kd, query, data, D);
        *p_ibest = ibest;
    }
}

void MANGLE(kdtree_nn)(const kdtree_t* kd, const void* vquery,
                       double* p_bestd2, int* p_ibest) {
    int nodestack[100];
    double dist2stack[100];
    int stackpos = 0;
    int D = (kd ? kd->ndim : 0);

    double bestd2 = *p_bestd2;
    int ibest = *p_ibest;
    const etype* query = vquery;

    if (!kd) {
        WARNING("kdtree_nn: null tree!\n");
        return;
    }

    // Bounding boxes
    if (!kd->split.any) {
        kdtree_nn_bb(kd, query, p_bestd2, p_ibest);
        return;
    }

#if defined(KD_DIM)
    assert(kd->ndim == KD_DIM);
    D = KD_DIM;
#else
    D = kd->ndim;
#endif

    // Integers.
    if (TTYPE_INTEGER) {
        ttype tquery[D];
        if (ttype_query(kd, query, tquery)) {
            kdtree_nn_int_split(kd, query, tquery, p_bestd2, p_ibest);
            return;
        }
    }

    // We got splitting planes, and the splits are either doubles, or ints
    // but the query doesn't find into the integer range.

    // queue root.
    nodestack[0] = 0;
    dist2stack[0] = 0.0;
    if (kd->fun.nn_enqueue)
        kd->fun.nn_enqueue(kd, 0, 1);

    while (stackpos >= 0) {
        int nodeid;
        int i;
        int dim = -1;
        int L, R;
        ttype split = 0;
        double del;
        etype rsplit;

        int nearchild;
        int farchild;
        double fard2;

        if (dist2stack[stackpos] > bestd2) {
            // pruned!
            if (kd->fun.nn_prune)
                kd->fun.nn_prune(kd, nodestack[stackpos], dist2stack[stackpos], bestd2, 1);
            stackpos--;
            continue;
        }
        nodeid = nodestack[stackpos];
        stackpos--;

        if (kd->fun.nn_explore)
            kd->fun.nn_explore(kd, nodeid, dist2stack[stackpos+1], bestd2);

        if (KD_IS_LEAF(kd, nodeid)) {
            dtype* data;
            L = kdtree_left(kd, nodeid);
            R = kdtree_right(kd, nodeid);
            for (i=L; i<=R; i++) {
                anbool bailedout = FALSE;
                double dsqd;

                if (kd->fun.nn_point)
                    kd->fun.nn_point(kd, nodeid, i);

                data = KD_DATA(kd, D, i);
                dist2_bailout(kd, query, data, D, bestd2, &bailedout, &dsqd);
                if (bailedout)
                    continue;
                // new best
                ibest = i;
                bestd2 = dsqd;

                if (kd->fun.nn_new_best)
                    kd->fun.nn_new_best(kd, nodeid, i, bestd2);
            }
            continue;
        }

        // split/dim trees
        split = *KD_SPLIT(kd, nodeid);
        if (kd->splitdim) {
            dim = kd->splitdim[nodeid];
        } else {
            // packed int
            bigint tmpsplit = split;
            dim = tmpsplit & kd->dimmask;
            split = tmpsplit & kd->splitmask;
        }
        rsplit = POINT_TE(kd, dim, split);
        del = query[dim] - rsplit;
        fard2 = del*del;
        if (query[dim] < rsplit) {
            nearchild = KD_CHILD_LEFT (nodeid);
            farchild  = KD_CHILD_RIGHT(nodeid);
        } else {
            nearchild = KD_CHILD_RIGHT(nodeid);
            farchild  = KD_CHILD_LEFT (nodeid);
        }

        if (fard2 <= bestd2) {
            // is the far child within range?
            stackpos++;
            nodestack[stackpos] = farchild;
            dist2stack[stackpos] = fard2;
            if (kd->fun.nn_enqueue)
                kd->fun.nn_enqueue(kd, farchild, 8);
        } else {
            if (kd->fun.nn_prune)
                kd->fun.nn_prune(kd, farchild, fard2, bestd2, 7);
        }

        // stack near child.
        stackpos++;
        nodestack[stackpos] = nearchild;
        dist2stack[stackpos] = 0.0;
        if (kd->fun.nn_enqueue)
            kd->fun.nn_enqueue(kd, nearchild, 9);
    }
    *p_bestd2 = bestd2;
    *p_ibest = ibest;
}


static inline int MANGLE(kdtree_rangesearch_continuation_push)(
    kdtree_rangesearch_continuation_t* continuation,
    int nodeid) {
    if (!continuation) {
        return -1;
    }

    if ((continuation->stackpos + 1) >=
        KDTREE_RANGESEARCH_CONTINUATION_STACK_MAX) {
        ERROR("KD range-search continuation stack overflow at node %i",
              nodeid);
        return -1;
    }

    continuation->stackpos++;
    continuation->nodestack[continuation->stackpos] = nodeid;

    return 0;
}

/*
 * Process exactly one real split-tree traversal node.
 *
 * This implementation mirrors the established scalar split path while the
 * scalar function remains an independent correctness fallback during the
 * continuation gate.  Child push order matches the historical stack exactly,
 * preserving result order.
 */
static inline int MANGLE(kdtree_rangesearch_continuation_visit_node)(
    kdtree_rangesearch_continuation_t* continuation,
    int nodeid) {
    const kdtree_t* kd;
    const etype* query;
    ttype* tquery;
    int D;
    int dim = -1;
    ttype split = 0;

    if (!continuation || !continuation->tree || !continuation->query ||
        !continuation->result) {
        return -1;
    }

    kd = continuation->tree;
    query = continuation->query;
    tquery = continuation->tquery.TTYPE;
    D = continuation->dimension;

    if (nodeid < 0 ||
        (nodeid >= kd->ninterior &&
         nodeid >= kd->ninterior + kd->nbottom)) {
        ERROR("KD range-search continuation encountered invalid node %i",
              nodeid);
        return -1;
    }

    if (KD_IS_LEAF(kd, nodeid)) {
        int L;
        int R;
        int i;

        L = kdtree_left(kd, nodeid);
        R = kdtree_right(kd, nodeid);

        if (L < 0 || R < L || R >= kd->ndata) {
            ERROR("KD range-search continuation encountered invalid leaf "
                  "range node=%i L=%i R=%i ndata=%i",
                  nodeid,
                  L,
                  R,
                  kd->ndata);
            return -1;
        }

        if (continuation->do_dists) {
            for (i = L; i <= R; i++) {
                anbool bailedout = FALSE;
                double dsqd;
                dtype* data;

                data = KD_DATA(kd, D, i);

                dist2_bailout(kd,
                              query,
                              data,
                              D,
                              continuation->maxd2,
                              &bailedout,
                              &dsqd);

                if (bailedout) {
                    continue;
                }

                if (!add_result(kd,
                                continuation->result,
                                dsqd,
                                KD_PERM(kd, i),
                                data,
                                D,
                                continuation->do_dists,
                                continuation->do_points)) {
                    return -1;
                }
            }
        } else {
            for (i = L; i <= R; i++) {
                dtype* data;

                data = KD_DATA(kd, D, i);

                if (dist2_exceeds(kd,
                                  query,
                                  data,
                                  D,
                                  continuation->maxd2)) {
                    continue;
                }

                if (!add_result(kd,
                                continuation->result,
                                LARGE_VAL,
                                KD_PERM(kd, i),
                                data,
                                D,
                                continuation->do_dists,
                                continuation->do_points)) {
                    return -1;
                }
            }
        }

        return 0;
    }

    if (kd->splitdim) {
        dim = kd->splitdim[nodeid];
    }

    split = *KD_SPLIT(kd, nodeid);

    if (!kd->splitdim && TTYPE_INTEGER) {
        bigint tmpsplit;

        tmpsplit = split;
        dim = tmpsplit & kd->dimmask;
        split = tmpsplit & kd->splitmask;
    }

    if (dim < 0 || dim >= D) {
        ERROR("KD range-search continuation encountered invalid split "
              "dimension node=%i dim=%i ndim=%i",
              nodeid,
              dim,
              D);
        return -1;
    }

    if (TTYPE_INTEGER && continuation->use_tsplit) {
        if (tquery[dim] < split) {
            if (MANGLE(kdtree_rangesearch_continuation_push)(
                    continuation,
                    KD_CHILD_LEFT(nodeid))) {
                return -1;
            }

            if (split - tquery[dim] <=
                continuation->tlinf.TTYPE) {
                if (MANGLE(kdtree_rangesearch_continuation_push)(
                        continuation,
                        KD_CHILD_RIGHT(nodeid))) {
                    return -1;
                }
            }
        } else {
            if (MANGLE(kdtree_rangesearch_continuation_push)(
                    continuation,
                    KD_CHILD_RIGHT(nodeid))) {
                return -1;
            }

            if (tquery[dim] - split <=
                continuation->tlinf.TTYPE) {
                if (MANGLE(kdtree_rangesearch_continuation_push)(
                        continuation,
                        KD_CHILD_LEFT(nodeid))) {
                    return -1;
                }
            }
        }
    } else {
        dtype rsplit;

        rsplit = POINT_TE(kd, dim, split);

        if (query[dim] < rsplit) {
            if (MANGLE(kdtree_rangesearch_continuation_push)(
                    continuation,
                    KD_CHILD_LEFT(nodeid))) {
                return -1;
            }

            if (rsplit - query[dim] <= continuation->maxdist) {
                if (MANGLE(kdtree_rangesearch_continuation_push)(
                        continuation,
                        KD_CHILD_RIGHT(nodeid))) {
                    return -1;
                }
            }
        } else {
            if (MANGLE(kdtree_rangesearch_continuation_push)(
                    continuation,
                    KD_CHILD_RIGHT(nodeid))) {
                return -1;
            }

            if (query[dim] - rsplit <= continuation->maxdist) {
                if (MANGLE(kdtree_rangesearch_continuation_push)(
                        continuation,
                        KD_CHILD_LEFT(nodeid))) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

int MANGLE(kdtree_rangesearch_continuation_init)(
    kdtree_rangesearch_continuation_t* continuation,
    const kdtree_t* kd,
    kdtree_qres_t* res,
    const void* vquery,
    double maxd2,
    int options) {
    const etype* query = vquery;
    int D;
    double dtlinf = 0.0;

    if (!continuation) {
        return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
    }

    /*
     * Initialize only lifecycle fields before validating caller input.  The
     * previous implementation cleared the complete continuation object in
     * execute(), init(), and cleanup().  That wrote more than one kilobyte
     * per CodeKD query even when only four query dimensions were active.
     */
    continuation->tree = NULL;
    continuation->query = NULL;
    continuation->result = NULL;
    continuation->stackpos = -1;
    continuation->owns_result = FALSE;
    continuation->result_released = FALSE;
    continuation->nodes_visited = 0;
    continuation->status =
        KDTREE_RANGESEARCH_CONTINUATION_UNINITIALIZED;

    if (!kd || !query) {
        return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
    }

#if defined(KD_DIM)
    assert(kd->ndim == KD_DIM);
    D = KD_DIM;
#else
    D = kd->ndim;
#endif

    if (D <= 0 || D > KDTREE_RANGESEARCH_CONTINUATION_DIM_MAX) {
        return KDTREE_RANGESEARCH_CONTINUATION_INIT_UNSUPPORTED;
    }

    /*
     * The first continuation gate deliberately supports the real split-tree
     * scalar path used by CodeKD.  Bounding-box traversal remains on the
     * established scalar implementation until it receives its own parity
     * gate.
     */
    if (!kd->split.any ||
        (kd->bb.any && !(options & KD_OPTIONS_USE_SPLIT))) {
        return KDTREE_RANGESEARCH_CONTINUATION_INIT_UNSUPPORTED;
    }

    if (options & KD_OPTIONS_SORT_DISTS) {
        options |= KD_OPTIONS_COMPUTE_DISTS;
    }

    continuation->tree = kd;
    continuation->query = query;
    continuation->maxd2 = maxd2;
    continuation->maxdist = sqrt(maxd2);
    continuation->options = options;
    continuation->dimension = D;
    continuation->do_dists =
        (options & KD_OPTIONS_COMPUTE_DISTS) != 0;
    continuation->do_points =
        (options & KD_OPTIONS_RETURN_POINTS) != 0;
    continuation->use_tquery = FALSE;
    continuation->use_tsplit = FALSE;
    continuation->tlinf.TTYPE = 0;

    if (TTYPE_INTEGER) {
        continuation->use_tquery =
            ttype_query(kd,
                        query,
                        continuation->tquery.TTYPE);
    }

    if (TTYPE_INTEGER && continuation->use_tquery) {
        dtlinf = DIST_ET(kd, continuation->maxdist, );
        continuation->tlinf.TTYPE = ceil(dtlinf);
    }

    continuation->use_tsplit =
        continuation->use_tquery &&
        (dtlinf < TTYPE_MAX);

    continuation->owns_result = (res == NULL);

    if (res) {
        if (!res->capacity) {
            if (!resize_results(res,
                                KDTREE_MAX_RESULTS,
                                D,
                                continuation->do_dists,
                                continuation->do_points)) {
                return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
            }
        } else {
            if (!resize_results(res,
                                res->capacity,
                                D,
                                continuation->do_dists,
                                continuation->do_points)) {
                return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
            }
        }

        res->nres = 0;
    } else {
        res = CALLOC(1, sizeof(kdtree_qres_t));
        if (!res) {
            SYSERROR("Failed to allocate kdtree_qres_t struct");
            return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
        }

        if (!resize_results(res,
                            KDTREE_MAX_RESULTS,
                            D,
                            continuation->do_dists,
                            continuation->do_points)) {
            kdtree_free_query(res);
            return KDTREE_RANGESEARCH_CONTINUATION_INIT_ERROR;
        }
    }

    continuation->result = res;
    continuation->stackpos = 0;
    continuation->nodestack[0] = 0;
    continuation->status = KDTREE_RANGESEARCH_CONTINUATION_MORE;

    return KDTREE_RANGESEARCH_CONTINUATION_INIT_OK;
}

int MANGLE(kdtree_rangesearch_continuation_step)(
    kdtree_rangesearch_continuation_t* continuation,
    size_t node_budget) {
    size_t processed = 0;

    if (!continuation || !continuation->tree ||
        !continuation->result || !node_budget) {
        if (continuation) {
            continuation->status = KDTREE_RANGESEARCH_CONTINUATION_ERROR;
        }
        return KDTREE_RANGESEARCH_CONTINUATION_ERROR;
    }

    if (continuation->status !=
        KDTREE_RANGESEARCH_CONTINUATION_MORE) {
        return continuation->status;
    }

    while (continuation->stackpos >= 0 && processed < node_budget) {
        int nodeid;

        nodeid = continuation->nodestack[continuation->stackpos];
        continuation->stackpos--;

        if (MANGLE(kdtree_rangesearch_continuation_visit_node)(
                continuation,
                nodeid)) {
            continuation->status = KDTREE_RANGESEARCH_CONTINUATION_ERROR;
            return continuation->status;
        }

        continuation->nodes_visited++;
        processed++;
    }

    if (continuation->stackpos < 0) {
        continuation->status = KDTREE_RANGESEARCH_CONTINUATION_DONE;
    }

    return continuation->status;
}

kdtree_qres_t* MANGLE(kdtree_rangesearch_continuation_finish)(
    kdtree_rangesearch_continuation_t* continuation) {
    kdtree_qres_t* result;

    if (!continuation || !continuation->tree ||
        !continuation->result ||
        continuation->status != KDTREE_RANGESEARCH_CONTINUATION_DONE) {
        return NULL;
    }

    if (!(continuation->options & KD_OPTIONS_NO_RESIZE_RESULTS)) {
        if (!resize_results(continuation->result,
                            continuation->result->nres,
                            continuation->dimension,
                            continuation->do_dists,
                            continuation->do_points)) {
            continuation->status = KDTREE_RANGESEARCH_CONTINUATION_ERROR;
            return NULL;
        }
    }

    if (continuation->options & KD_OPTIONS_SORT_DISTS) {
        kdtree_qsort_results(continuation->result,
                             continuation->tree->ndim);
    }

    result = continuation->result;
    continuation->result_released = TRUE;
    continuation->status = KDTREE_RANGESEARCH_CONTINUATION_FINISHED;

    return result;
}

void MANGLE(kdtree_rangesearch_continuation_cleanup)(
    kdtree_rangesearch_continuation_t* continuation) {
    if (!continuation) {
        return;
    }

    if (continuation->result && !continuation->result_released) {
        if (continuation->owns_result) {
            kdtree_free_query(continuation->result);
        } else {
            continuation->result->nres = 0;
        }
    }

    continuation->tree = NULL;
    continuation->query = NULL;
    continuation->result = NULL;
    continuation->stackpos = -1;
    continuation->owns_result = FALSE;
    continuation->result_released = FALSE;
    continuation->nodes_visited = 0;
    continuation->status =
        KDTREE_RANGESEARCH_CONTINUATION_UNINITIALIZED;
}

kdtree_qres_t* MANGLE(kdtree_rangesearch_options)
     (const kdtree_t* kd, kdtree_qres_t* res, const void* vquery,
      double maxd2, int options)
{
    int nodestack[100];
    int stackpos = 0;
    int D = (kd ? kd->ndim : 0);

    anbool do_dists;
    anbool do_points;
    anbool do_wholenode_check;
    double maxdist = 0.0;
    ttype tlinf = 0;
    ttype tl1 = 0;
    ttype tl2 = 0;
    bigttype bigtl2 = 0;

    anbool use_tquery = FALSE;
    anbool use_tsplit = FALSE;
    anbool use_tmath = FALSE;
    anbool use_bigtmath = FALSE;

    anbool do_precheck = FALSE;
    anbool do_l1precheck = FALSE;

    anbool use_bboxes = FALSE;
    Unused anbool use_splits = FALSE;

    double dtl1=0.0, dtl2=0.0, dtlinf=0.0;

    const etype* query = vquery;

    //dtype dquery[D];
    ttype tquery[D];

    if (!kd || !query)
        return NULL;
#if defined(KD_DIM)
    assert(kd->ndim == KD_DIM);
    D = KD_DIM;
#else
    D = kd->ndim;
#endif

    if (options & KD_OPTIONS_SORT_DISTS)
        // gotta compute 'em if ya wanna sort 'em!
        options |= KD_OPTIONS_COMPUTE_DISTS;
    do_dists = options & KD_OPTIONS_COMPUTE_DISTS;
    do_points = options & KD_OPTIONS_RETURN_POINTS;
    do_wholenode_check = !(options & KD_OPTIONS_SMALL_RADIUS);

    if ((options & KD_OPTIONS_SPLIT_PRECHECK) &&
        kd->bb.any && kd->splitdim) {
        do_precheck = TRUE;
    }

    if ((options & KD_OPTIONS_L1_PRECHECK) &&
        kd->bb.any) {
        do_l1precheck = TRUE;
    }

    if (!kd->split.any) {
        assert(kd->bb.any);
        use_bboxes = TRUE;
    } else {
        if (kd->bb.any) {
            // Got both BBoxes and Splits.
            if (options & KD_OPTIONS_USE_SPLIT) {
                use_splits = TRUE;
            } else {
                // Use bboxes by default.
                use_bboxes = TRUE;
            }
        } else {
            use_splits = TRUE;
            assert(kd->splitdim || TTYPE_INTEGER);
        }
    }

     assert(use_splits || use_bboxes);

    maxdist = sqrt(maxd2);

    if (TTYPE_INTEGER &&
        (kd->split.any || do_precheck || do_l1precheck)) {
        use_tquery = ttype_query(kd, query, tquery);
    }

    if (TTYPE_INTEGER && use_tquery) {
        dtl1   = DIST_ET(kd, maxdist * sqrt(D),);
        dtl2   = DIST2_ET(kd, maxd2, );
        dtlinf = DIST_ET(kd, maxdist, );
        tl1    = ceil(dtl1);
        tlinf  = ceil(dtlinf);
        bigtl2 = ceil(dtl2);
        tl2    = bigtl2;
    }

    use_tsplit = use_tquery && (dtlinf < TTYPE_MAX);

    if (do_l1precheck)
        if (dtl1 > TTYPE_MAX) {
            //printf("L1 maxdist %g overflows ttype representation.  L1 precheck disabled.\n", dtl1);
            do_l1precheck = FALSE;
        }

    if (TTYPE_INTEGER && use_tquery && kd->bb.any) {
        if (dtl2 < TTYPE_MAX) {
            use_tmath = TRUE;
            /*
             printf("Using %s integer math.\n", STRINGIFY(TTYPE));
             printf("(tl2 = %u).\n", (unsigned int)tl2);
             */
        } else if (dtl2 < BIGTTYPE_MAX) {
            use_bigtmath = TRUE;
        } else {
            /*
             printf("L2 maxdist overflows u16 and u32 representation; not using int math.  %g -> %g > %u\n",
             maxd2, dtl2, UINT32_MAX);
             */
        }
        if (use_bigtmath) {
            if (options & KD_OPTIONS_NO_BIG_INT_MATH)
                use_bigtmath = FALSE;
            else {
                /*
                 printf("Using %s/%s integer math.\n", STRINGIFY(TTYPE), STRINGIFY(BIGTTYPE));
                 printf("(bigtl2 = %llu).\n", (long long unsigned int)bigtl2);
                 */
            }
        }
    }


    if (res) {
        if (!res->capacity) {
            resize_results(res, KDTREE_MAX_RESULTS, D, do_dists, do_points);
        } else {
            // call the resize routine just in case the old result struct was
            // from a tree of different type or dimensionality.
            resize_results(res, res->capacity, D, do_dists, do_points);
        }
        res->nres = 0;
    } else {
        res = CALLOC(1, sizeof(kdtree_qres_t));
        if (!res) {
            SYSERROR("Failed to allocate kdtree_qres_t struct");
            return NULL;
        }
        resize_results(res, KDTREE_MAX_RESULTS, D, do_dists, do_points);
    }

    // queue root.
    nodestack[0] = 0;

    while (stackpos >= 0) {
        int nodeid;
        int i;
        int dim = -1;
        int L, R;
        ttype split = 0;
        ttype *tlo=NULL, *thi=NULL;

        nodeid = nodestack[stackpos];
        stackpos--;


        if (KD_IS_LEAF(kd, nodeid)) {
            dtype* data;

            L = kdtree_left(kd, nodeid);
            R = kdtree_right(kd, nodeid);

            if (do_dists) {
                for (i = L; i <= R; i++) {
                    anbool bailedout = FALSE;
                    double dsqd;

                    data = KD_DATA(kd, D, i);

                    dist2_bailout(kd,
                                  query,
                                  data,
                                  D,
                                  maxd2,
                                  &bailedout,
                                  &dsqd);

                    if (bailedout) {
                        continue;
                    }

                    if (!add_result(kd,
                                    res,
                                    dsqd,
                                    KD_PERM(kd, i),
                                    data,
                                    D,
                                    do_dists,
                                    do_points)) {
                        return NULL;
                    }
                }
            } else {
                for (i = L; i <= R; i++) {
                    data = KD_DATA(kd, D, i);

                    if (dist2_exceeds(kd,
                                      query,
                                      data,
                                      D,
                                      maxd2)) {
                        continue;
                    }

                    if (!add_result(kd,
                                    res,
                                    LARGE_VAL,
                                    KD_PERM(kd, i),
                                    data,
                                    D,
                                    do_dists,
                                    do_points)) {
                        return NULL;
                    }
                }
            }
            continue;
        }

        if (kd->splitdim)
            dim = kd->splitdim[nodeid];

        if (use_bboxes) {
            anbool wholenode = FALSE;

            bboxes(kd, nodeid, &tlo, &thi, D);
            assert(tlo && thi);

            if (do_precheck && nodeid) {
                anbool isleftchild = KD_IS_LEFT_CHILD(nodeid);
                // we need to use the dimension our _parent_ split on, not ours!
                int pdim;
                anbool cut;
                if (kd->splitdim)
                    pdim = kd->splitdim[KD_PARENT(nodeid)];
                else {
                    pdim = kd->split.TTYPE[KD_PARENT(nodeid)];
                    pdim &= kd->dimmask;
                }
                if (TTYPE_INTEGER && use_tquery) {
                    if (isleftchild)
                        cut = ((tquery[pdim] > thi[pdim]) &&
                               (tquery[pdim] - thi[pdim] > tlinf));
                    else
                        cut = ((tlo[pdim] > tquery[pdim]) &&
                               (tlo[pdim] - tquery[pdim] > tlinf));
                } else {
                    etype bb;
                    if (isleftchild) {
                        bb = POINT_TE(kd, pdim, thi[pdim]);
                        cut = (query[pdim] - bb > maxdist);
                    } else {
                        bb = POINT_TE(kd, pdim, tlo[pdim]);
                        cut = (bb - query[pdim] > maxdist);
                    }
                }
                if (cut) {
                    // precheck failed!
                    //printf("precheck failed!\n");
                    continue;
                }
            }

            if (TTYPE_INTEGER && do_l1precheck && use_tquery)
                if (bb_point_l1mindist_exceeds_ttype(tlo, thi, tquery, D, tl1, tlinf)) {
                    //printf("l1 precheck failed!\n");
                    continue;
                }

            if (TTYPE_INTEGER && use_tmath) {
                if (bb_point_mindist2_exceeds_ttype(tlo, thi, tquery, D, tl2))
                    continue;
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds_ttype(tlo, thi, tquery, D, tl2);
            } else if (TTYPE_INTEGER && use_bigtmath) {
                if (bb_point_mindist2_exceeds_bigttype(tlo, thi, tquery, D, bigtl2))
                    continue;
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds_bigttype(tlo, thi, tquery, D, bigtl2);
            } else {
                etype bblo[D], bbhi[D];
                int d;
                for (d=0; d<D; d++) {
                    bblo[d] = POINT_TE(kd, d, tlo[d]);
                    bbhi[d] = POINT_TE(kd, d, thi[d]);
                }
                if (bb_point_mindist2_exceeds(bblo, bbhi, query, D, maxd2))
                    continue;
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds(bblo, bbhi, query, D, maxd2);
            }

            if (wholenode) {
                L = kdtree_left(kd, nodeid);
                R = kdtree_right(kd, nodeid);
                if (do_dists) {
                    for (i=L; i<=R; i++) {
                        double dsqd = dist2(kd, query, KD_DATA(kd, D, i), D);
                        if (!add_result(kd, res, dsqd, KD_PERM(kd, i),
                                        KD_DATA(kd, D, i), D,
                                        do_dists, do_points))
                            return NULL;
                    }
                } else {
                    for (i=L; i<=R; i++)
                        if (!add_result(kd, res, LARGE_VAL, KD_PERM(kd, i),
                                        KD_DATA(kd, D, i), D,
                                        do_dists, do_points))
                            return NULL;
                }
                continue;
            }

            stackpos++;
            nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
            stackpos++;
            nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);

        } else {
            // use_splits.

            split = *KD_SPLIT(kd, nodeid);
            if (!kd->splitdim && TTYPE_INTEGER) {
                bigint tmpsplit;
                tmpsplit = split;
                dim = tmpsplit & kd->dimmask;
                split = tmpsplit & kd->splitmask;
            }

            if (TTYPE_INTEGER && use_tsplit) {

                if (tquery[dim] < split) {
                    // query is on the "left" side of the split.
                    assert(query[dim] < POINT_TE(kd, dim, split));
                    stackpos++;
                    nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
                    // look mum, no int overflow!
                    if (split - tquery[dim] <= tlinf) {
                        // right child is okay.
                        assert(POINT_TE(kd, dim, split) - query[dim] >= 0.0);
                        // This may fail due to rounding?
                        //assert(POINT_TE(kd, dim, split) - query[dim] <= maxdist);
                        stackpos++;
                        nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
                    }

                } else {
                    // query is on "right" side.
                    //assert(POINT_TE(kd, dim, split) <= query[dim]);
                    assert(POINT_TE(kd, dim, split) <= query[dim] + kd->invscale/2.0);
                    stackpos++;
                    nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
                    if (tquery[dim] - split <= tlinf) {
                        //assert(query[dim] - POINT_TE(kd, dim, split) >= 0.0);
                        assert((query[dim] - POINT_TE(kd, dim, split)) >= -kd->invscale/2.0);
                        // This may fail due to rounding?
                        //assert(query[dim] - POINT_TE(kd, dim, split) <= maxdist);
                        stackpos++;
                        nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
                    }
                }
            } else {
                dtype rsplit = POINT_TE(kd, dim, split);
                if (query[dim] < rsplit) {
                    // query is on the "left" side of the split.
                    stackpos++;
                    nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
                    if (rsplit - query[dim] <= maxdist) {
                        stackpos++;
                        nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
                    }
                } else {
                    // query is on the "right" side
                    stackpos++;
                    nodestack[stackpos] = KD_CHILD_RIGHT(nodeid);
                    if (query[dim] - rsplit <= maxdist) {
                        stackpos++;
                        nodestack[stackpos] = KD_CHILD_LEFT(nodeid);
                    }
                }
            }
        }
    }

    /* Resize result arrays. */
    if (!(options & KD_OPTIONS_NO_RESIZE_RESULTS))
        resize_results(res, res->nres, D, do_dists, do_points);

    /* Sort by ascending distance away from target point before returning */
    if (options & KD_OPTIONS_SORT_DISTS) {
        if (FALSE) {
            printf("before sorting results:\n");
            print_results(res, D);
        }
        kdtree_qsort_results(res, kd->ndim);
        if (FALSE) {
            printf("after sorting results:\n");
            print_results(res, D);
        }
    }

    return res;
}


static void* get_data(const kdtree_t* kd, int i) {
    return KD_DATA(kd, kd->ndim, i);
}

static void copy_data_double(const kdtree_t* kd, int start, int N,
                             double* dest) {
    Unused int i, j, d;
    int D;

    D = kd->ndim;
#if DTYPE_DOUBLE
    //#warning "Data type is double; just copying data."
    memcpy(dest, kd->data.DTYPE + start*D, (size_t)N*(size_t)D*sizeof(etype));
#elif (!DTYPE_INTEGER && !ETYPE_INTEGER)
    //#warning "Etype and Dtype are both reals; just casting values."
    for (i=0; i<(N * D); i++)
        dest[i] = kd->data.DTYPE[start*D + i];
#else
    //#warning "Using POINT_DE to convert data."
    j=0;
    for (i=0; i<N; i++)
        for (d=0; d<D; d++) {
            dest[j] = POINT_DE(kd, D, kd->data.DTYPE[(start + i)*D + d]);
            j++;
        }
#endif
}

static KDTREE_THREAD_LOCAL dtype* kdqsort_arr;
static KDTREE_THREAD_LOCAL int kdqsort_D;

static int kdqsort_compare(const void* v1, const void* v2)
{
    int i1, i2;
    dtype val1, val2;

    i1 = *((int*)v1);
    i2 = *((int*)v2);
    val1 = kdqsort_arr[(size_t)i1 * (size_t)kdqsort_D];
    val2 = kdqsort_arr[(size_t)i2 * (size_t)kdqsort_D];
    if (val1 < val2) {
        return -1;
    } else if (val1 > val2) {
        return 1;
    }
    return 0;
}

static int kdtree_qsort(dtype *arr, unsigned int *parr, int l, int r, int D, int d)
{
    int* permute;
    int i, j, N;
    dtype* tmparr;
    int* tmpparr;

    N = r - l + 1;
    permute = MALLOC((size_t)N * sizeof(int));
    if (!permute) {
        SYSERROR("Failed to allocate extra permutation array");
        return -1;
    }
    for (i = 0; i < N; i++)
        permute[i] = i;

    kdqsort_arr = arr + (size_t)l * (size_t)D + (size_t)d;
    kdqsort_D = D;
    qsort(permute, N, sizeof(int), kdqsort_compare);

    // permute the data one dimension at a time...
    tmparr = MALLOC((size_t)N * sizeof(dtype));
    if (!tmparr) {
        SYSERROR("Failed to allocate temp permutation array");
        FREE(permute);
        return -1;
    }
    for (j = 0; j < D; j++) {
        for (i = 0; i < N; i++) {
            int pi = permute[i];
            tmparr[i] = arr[(size_t)(l + pi) * (size_t)D + (size_t)j];
        }
        for (i = 0; i < N; i++)
            arr[(size_t)(l + i) * (size_t)D + (size_t)j] = tmparr[i];
    }
    FREE(tmparr);
    tmpparr = MALLOC((size_t)N * sizeof(int));
    if (!tmpparr) {
        SYSERROR("Failed to allocate temp permutation array");
        FREE(permute);
        return -1;
    }
    for (i = 0; i < N; i++) {
        int pi = permute[i];
        tmpparr[i] = parr[l + pi];
    }
    memcpy(parr + l, tmpparr, (size_t)N*(size_t)sizeof(int));
    FREE(tmpparr);
    FREE(permute);
    return 0;
}


#define GET(x) (arr[(size_t)(x)*(size_t)D+(size_t)d])
#if defined(KD_DIM)
#define ELEM_SWAP(il, ir) {                                             \
        if ((il) != (ir)) {                                             \
            tmpperm  = parr[il];                                        \
            assert(tmpperm != -1);                                      \
            parr[il] = parr[ir];                                        \
            parr[ir] = tmpperm;                                         \
            { int d; for (d=0; d<D; d++) {                              \
                    tmpdata[0] = KD_ARRAY_VAL(arr, D, il, d);           \
                    *KD_ARRAY_REF(arr, D, il, d) = KD_ARRAY_VAL(arr, D, ir, d); \
                    *KD_ARRAY_REF(arr, D, ir, d) = tmpdata[0]; }}}}
#else
#define ELEM_SWAP(il, ir) {                                             \
        if ((il) != (ir)) {                                             \
            tmpperm  = parr[il];                                        \
            assert(tmpperm != -1);                                      \
            parr[il] = parr[ir];                                        \
            parr[ir] = tmpperm;                                         \
            memcpy(tmpdata,                     KD_ARRAY_REF(arr, D, il, 0), D*sizeof(dtype)); \
            memcpy(KD_ARRAY_REF(arr, D, il, 0), KD_ARRAY_REF(arr, D, ir, 0), D*sizeof(dtype)); \
            memcpy(KD_ARRAY_REF(arr, D, ir, 0), tmpdata,                     D*sizeof(dtype)); }}
#endif
#define ELEM_ROT(iA, iB, iC) {                                          \
        tmpperm  = parr[iC];                                            \
        parr[iC] = parr[iB];                                            \
        parr[iB] = parr[iA];                                            \
        parr[iA] = tmpperm;                                             \
        assert(tmpperm != -1);                                          \
        memcpy(tmpdata,                     KD_ARRAY_REF(arr, D, iC, 0), D*sizeof(dtype)); \
        memcpy(KD_ARRAY_REF(arr, D, iC, 0), KD_ARRAY_REF(arr, D, iB, 0), D*sizeof(dtype)); \
        memcpy(KD_ARRAY_REF(arr, D, iB, 0), KD_ARRAY_REF(arr, D, iA, 0), D*sizeof(dtype)); \
        memcpy(KD_ARRAY_REF(arr, D, iA, 0), tmpdata,                     D*sizeof(dtype)); }

static void kdtree_quickselect_partition(dtype *arr, unsigned int *parr,
                                         int L, int R, int D, int d,
                                         int rank) {
    int i;
    int low, high;

#if defined(KD_DIM)
    // tell the compiler this is a constant...
    D = KD_DIM;
#endif

    /* sanity is good */
    assert(R >= L);

    /* Find the "rank"th point and partition the data. */
    /* For us, "rank" is usually the median of L and R. */
    low = L;
    high = R;
    while(1) {
        dtype vals[3];
        dtype tmp;
        dtype pivot;
        int i,j;
        int iless, iequal, igreater;
        int endless, endequal, endgreater;
        int middle;
        int nless, nequal;
        // temp storage for ELEM_SWAP and ELEM_ROT macros.
        dtype tmpdata[D];
        int tmpperm;

        if (high == low)
            break;

        /* Choose the pivot: find the median of the values in low,
         middle, and high positions. */
        middle = (low + high) / 2;
        vals[0] = GET(low);
        vals[1] = GET(middle);
        vals[2] = GET(high);
        /* (Bubblesort the three elements.) */
        for (i=0; i<2; i++)
            for (j=0; j<(2-i); j++)
                if (vals[j] > vals[j+1]) {
                    tmp = vals[j];
                    vals[j] = vals[j+1];
                    vals[j+1] = tmp;
                }
        assert(vals[0] <= vals[1]);
        assert(vals[1] <= vals[2]);
        pivot = vals[1];

        /* Count the number of items that are less than, and equal to, the pivot. */
        nless = nequal = 0;
        for (i=low; i<=high; i++) {
            if (GET(i) < pivot)
                nless++;
            else if (GET(i) == pivot)
                nequal++;
        }

        /* These are the indices where the <, =, and > entries will start. */
        iless = low;
        iequal = low + nless;
        igreater = low + nless + nequal;

        /* These are the indices where they will end; ie the elements less than the
         pivot will live in [iless, endless).  (But note that we'll be incrementing
         "iequal" et al in the loop below.) */
        endless = iequal;
        endequal = igreater;
        endgreater = high+1;

        while (1) {
            /* Find an element in the "less" section that is out of place. */
            while ( (iless < endless) && (GET(iless) < pivot) )
                iless++;

            /* Find an element in the "equal" section that is out of place. */
            while ( (iequal < endequal) && (GET(iequal) == pivot) )
                iequal++;

            /* Find an element in the "greater" section that is out of place. */
            while ( (igreater < endgreater) && (GET(igreater) > pivot)  )
                igreater++;


            /* We're looking at three positions, and each one has three cases:
             we're finished that segment, or the element we're looking at belongs in
             one of the other two segments.  This yields 27 cases, but many of them
             are ruled out because, eg, if the element at "iequal" belongs in the "less"
             segment, then we can't be done the "less" segment.

             It turns out there are only 6 cases to handle:

             ---------------------------------------------
             case   iless    iequal   igreater    action
             ---------------------------------------------
             1      D        D        D           done
             2      G        ?        L           swap l,g
             3      E        L        ?           swap l,e
             4      ?        G        E           swap e,g
             5      E        G        L           rotate A
             6      G        L        E           rotate B
             ---------------------------------------------

             legend:
             D: done
             ?: don't care
             L: (element < pivot)
             E: (element == pivot)
             G: (element > pivot)
             */

            /* case 1: done? */
            if ((iless == endless) && (iequal == endequal) && (igreater == endgreater))
                break;

            /* case 2: swap l,g */
            if ((iless < endless) && (igreater < endgreater) &&
                (GET(iless) > pivot) && (GET(igreater) < pivot)) {
                ELEM_SWAP(iless, igreater);
                assert(GET(iless) < pivot);
                assert(GET(igreater) > pivot);
                continue;
            }

            /* cases 3,4,5,6 */
            assert(iequal < endequal);
            if (GET(iequal) > pivot) {
                /* cases 4,5: */
                assert(igreater < endgreater);
                if (GET(igreater) == pivot) {
                    /* case 4: swap e,g */
                    ELEM_SWAP(iequal, igreater);
                    assert(GET(iequal) == pivot);
                    assert(GET(igreater) > pivot);
                } else {
                    /* case 5: rotate. */
                    assert(GET(iless) == pivot);
                    assert(GET(iequal) > pivot);
                    assert(GET(igreater) < pivot);
                    ELEM_ROT(iless, iequal, igreater);
                    assert(GET(iless) < pivot);
                    assert(GET(iequal) == pivot);
                    assert(GET(igreater) > pivot);
                }
            } else {
                /* cases 3,6 */
                assert(GET(iequal) < pivot);
                assert(iless < endless);
                if (GET(iless) == pivot) {
                    /* case 3: swap l,e */
                    ELEM_SWAP(iless, iequal);
                    assert(GET(iless) < pivot);
                    assert(GET(iequal) == pivot);
                } else {
                    /* case 6: rotate. */
                    assert(GET(iless) > pivot);
                    assert(GET(iequal) < pivot);
                    assert(GET(igreater) == pivot);
                    ELEM_ROT(igreater, iequal, iless);
                    assert(GET(iless) < pivot);
                    assert(GET(iequal) == pivot);
                    assert(GET(igreater) > pivot);
                }
            }
        }

        /* Reset the indices of where the segments start. */
        iless = low;
        iequal = low + nless;
        igreater = low + nless + nequal;

        /* Assert that "<" values are in the "less" partition, "=" values are in the
         "equal" partition, and ">" values are in the "greater" partition. */
        for (i=iless; i<iequal; i++)
            assert(GET(i) < pivot);
        for (i=iequal; i<igreater; i++)
            assert(GET(i) == pivot);
        for (i=igreater; i<=high; i++)
            assert(GET(i) > pivot);

        /* Is the median in the "<", "=", or ">" partition? */
        if (rank < iequal)
            /* median is in the "<" partition.  low is unchanged. */
            high = iequal - 1;
        else if (rank < igreater)
            /* the median is inside the "=" partition; we're done! */
            break;
        else
            /* median is in the ">" partition.  high is unchanged. */
            low = igreater;
    }

    /* check that it worked. */
    for (i=L; i<rank; i++)
        assert(GET(i) <= GET(rank));
    for (i=rank; i<=R; i++)
        assert(GET(i) >= GET(rank));
}
#undef ELEM_SWAP
#undef ELEM_ROT
#undef GET



static int kdtree_check_node(const kdtree_t* kd, int nodeid) {
    int D = kd->ndim;
    int L, R;
    int i, d;

    L = kdtree_left (kd, nodeid);
    R = kdtree_right(kd, nodeid);

    if (kdtree_is_node_empty(kd, nodeid)) {
        assert(L == (R+1));
        assert(L >= 0);
        assert(L <= kd->ndata);
        assert(R >= -1);
        assert(R < kd->ndata);
        if (!((L == (R+1)) && (L >= -1) && (L <= kd->ndata) && (R >= -1) && (R < kd->ndata))) {
            ERROR("kdtree_check: L,R out of range for empty node");
            return -1;
        }
    } else {
        assert(L < kd->ndata);
        assert(R < kd->ndata);
        assert(L >= 0);
        assert(R >= 0);
        assert(L <= R);
        if (!((L < kd->ndata) && (R < kd->ndata) && (L >= 0) && (R >= 0) && (L <= R))) {
            ERROR("kdtree_check: L,R out of range for non-empty node");
            return -1;
        }
    }

    // if it's the root node, make sure that each value in the permutation
    // array is present exactly once.
    if (!nodeid && kd->perm) {
        unsigned char* counts = CALLOC(kd->ndata, 1);
        for (i=0; i<kd->ndata; i++)
            counts[kd->perm[i]]++;
        for (i=0; i<kd->ndata; i++)
            assert(counts[i] == 1);
        for (i=0; i<kd->ndata; i++)
            if (counts[i] != 1) {
                ERROR("kdtree_check: permutation vector failure");
                return -1;
            }
        FREE(counts);
    }

    if (kd->perm) {
        for (i=L; i<=R; i++) {
            assert(kd->perm[i] >= 0);
            assert(kd->perm[i] < kd->ndata);
            if (kd->perm[i] >= kd->ndata) {
                ERROR("kdtree_check: permutation vector range failure");
                return -1;
            }
        }
    }

    if (KD_IS_LEAF(kd, nodeid)) {
        if ((kd->minval && !kd->maxval) ||
            (!kd->minval && kd->maxval)) {
            ERROR("kdtree_check: minval but no maxval (or vice versa)");
            return -1;
        }
        if (kd->minval && kd->maxval) {
            for (i=L; i<=R; i++) {
                dtype* dat = KD_DATA(kd, D, i);
                for (d=0; d<D; d++) {
                    Unused etype e = POINT_DE(kd, d, dat[d]);
                    assert(e >= kd->minval[d]);
                    assert(e <= kd->maxval[d]);
                }
            }
        }
        return 0;
    }

    if (kd->bb.any) {
        ttype* bb;
        ttype *plo, *phi;
        anbool ok = FALSE;
        plo = LOW_HR( kd, D, nodeid);
        phi = HIGH_HR(kd, D, nodeid);

        // bounding-box sanity.
        for (d=0; d<D; d++) {
            assert(plo[d] <= phi[d]);
            if (plo[d] > phi[d]) {
                ERROR("kdtree_check: bounding-box sanity failure");
                return -1;
            }
        }
        // check that the points owned by this node are within its bounding box.
        for (i=L; i<=R; i++) {
            dtype* dat = KD_DATA(kd, D, i);
            for (d=0; d<D; d++) {
                Unused ttype t = POINT_DT(kd, d, dat[d], KD_ROUND);
                assert(plo[d] <= t);
                assert(t <= phi[d]);
                if (plo[d] > t || t > phi[d]) {
                    ERROR("kdtree_check: bounding-box failure");
                    return -1;
                }
            }
        }

        if (!KD_IS_LEAF(kd, nodeid)) {
            // check that the children's bounding box corners are within
            // the parent's bounding box.
            bb = LOW_HR(kd, D, KD_CHILD_LEFT(nodeid));
            for (d=0; d<D; d++) {
                assert(plo[d] <= bb[d]);
                assert(bb[d] <= phi[d]);
                if (plo[d] > bb[d] || bb[d] > phi[d]) {
                    ERROR("kdtree_check: bounding-box nesting failure");
                    return -1;
                }
            }
            bb = HIGH_HR(kd, D, KD_CHILD_LEFT(nodeid));
            for (d=0; d<D; d++) {
                assert(plo[d] <= bb[d]);
                assert(bb[d] <= phi[d]);
                if (plo[d] > bb[d] || bb[d] > phi[d]) {
                    ERROR("kdtree_check: bounding-box nesting failure");
                    return -1;
                }
            }
            bb = LOW_HR(kd, D, KD_CHILD_RIGHT(nodeid));
            for (d=0; d<D; d++) {
                assert(plo[d] <= bb[d]);
                assert(bb[d] <= phi[d]);
                if (plo[d] > bb[d] || bb[d] > phi[d]) {
                    ERROR("kdtree_check: bounding-box nesting failure");
                    return -1;
                }
            }
            bb = HIGH_HR(kd, D, KD_CHILD_RIGHT(nodeid));
            for (d=0; d<D; d++) {
                assert(plo[d] <= bb[d]);
                assert(bb[d] <= phi[d]);
                if (plo[d] > bb[d] || bb[d] > phi[d]) {
                    ERROR("kdtree_check: bounding-box nesting failure");
                    return -1;
                }
            }
            // check that the children don't overlap with positive volume.
            for (d=0; d<D; d++) {
                ttype* bb1 = HIGH_HR(kd, D, KD_CHILD_LEFT(nodeid));
                ttype* bb2 = LOW_HR(kd, D, KD_CHILD_RIGHT(nodeid));
                if (bb2[d] >= bb1[d]) {
                    ok = TRUE;
                    break;
                }
            }
            assert(ok);
            if (!ok) {
                ERROR("kdtree_check: peer overlap failure");
                return -1;
            }
        }
    }
    if (kd->split.any) {

        if (!KD_IS_LEAF(kd, nodeid)) {
            // check split/dim.
            ttype split;
            int dim = 0;
            int cL, cR;
            dtype dsplit;

            split = *KD_SPLIT(kd, nodeid);
            if (kd->splitdim)
                dim = kd->splitdim[nodeid];
            else {
                if (TTYPE_INTEGER) {
                    bigint tmpsplit;
                    tmpsplit = split;
                    dim = tmpsplit & kd->dimmask;
                    split = tmpsplit & kd->splitmask;
                }
            }

            dsplit = POINT_TD(kd, dim, split);

            cL = kdtree_left (kd, KD_CHILD_LEFT(nodeid));
            cR = kdtree_right(kd, KD_CHILD_LEFT(nodeid));
            for (i=cL; i<=cR; i++) {
                Unused dtype* dat = KD_DATA(kd, D, i);
                assert(dat[dim] <= dsplit);
                if (dat[dim] > dsplit) {
                    ERROR("kdtree_check: split-plane failure (1)");
                    printf("Data item %i, dim %i: %g vs %g\n", i, dim, (double)dat[dim], (double)dsplit);
                    return -1;
                }
            }

            cL = kdtree_left (kd, KD_CHILD_RIGHT(nodeid));
            cR = kdtree_right(kd, KD_CHILD_RIGHT(nodeid));
            for (i=cL; i<=cR; i++) {
                Unused dtype* dat = KD_DATA(kd, D, i);
                assert(dat[dim] >= dsplit);
                if (dat[dim] < dsplit) {
                    ERROR("kdtree_check: split-plane failure (2)");
                    return -1;
                }
            }
        }
    }

    return 0;
}

int MANGLE(kdtree_check)(const kdtree_t* kd) {
    int i;
    if (kd->split.any) {
        assert(kd->splitmask);
        if (!kd->splitdim) {
            assert(kd->dimmask);
        }
    }
    for (i=0; i<kd->nnodes; i++) {
        if (kdtree_check_node(kd, i))
            return -1;
    }
    return 0;
}

static double maxrange(double* lo, double* hi, int D) {
    double range;
    int d;
    range = 0.0;
    for (d=0; d<D; d++)
        if (hi[d] - lo[d] > range)
            range = hi[d] - lo[d];
    if (range == 0.0)
        return 1.0;
    return range;
}

static double compute_scale(dtype* ddata, int N, int D,
                            double* lo, double* hi) {
    double range;
    int i, d;
    for (d=0; d<D; d++) {
        lo[d] = DTYPE_MAX;
        hi[d] = DTYPE_MIN;
    }
    for (i=0; i<N; i++) {
        for (d=0; d<D; d++) {
            if (*ddata > hi[d]) hi[d] = *ddata;
            if (*ddata < lo[d]) lo[d] = *ddata;
            ddata++;
        }
    }
    range = maxrange(lo, hi, D);
    return (double)(DTYPE_MAX) / range;
}

// same as "compute_scale" but takes data of "etype".
static double compute_scale_ext(etype* edata, int N, int D,
                                double* lo, double* hi) {
    double range;
    int i, d;
    for (d=0; d<D; d++) {
        lo[d] = ETYPE_MAX;
        hi[d] = ETYPE_MIN;
    }
    for (i=0; i<N; i++) {
        for (d=0; d<D; d++) {
            if (*edata > hi[d]) hi[d] = *edata;
            if (*edata < lo[d]) lo[d] = *edata;
            edata++;
        }
    }
    range = maxrange(lo, hi, D);
    return (double)DTYPE_MAX / range;
}

static void convert_data(kdtree_t* kd, etype* edata, int N, int D, int Nleaf) {
    dtype* ddata;
    int i, d;

    if (!kd->minval || !kd->maxval) {
        kd->minval = MALLOC((size_t)D * sizeof(double));
        kd->maxval = MALLOC((size_t)D * sizeof(double));
        kd->scale = compute_scale_ext(edata, N, D, kd->minval, kd->maxval);
    } else {
        // limits were pre-set by the user.  just compute scale.
        double range;
        range = maxrange(kd->minval, kd->maxval, D);
        kd->scale = (double)DTYPE_MAX / range;
    }
    kd->invscale = 1.0 / kd->scale;

    ddata = kd->data.any = MALLOC((size_t)N * (size_t)D * (size_t)sizeof(dtype));
    if (!ddata) {
        ERROR("Failed to malloc %i x %i x %i\n", N, D, (int)sizeof(dtype));
        return;
    }
    kd->free_data = TRUE;
    for (i=0; i<N; i++) {
        for (d=0; d<D; d++) {
            etype dd = POINT_ED(kd, d, *edata, KD_ROUND);
            if (dd > DTYPE_MAX) {
                WARNING("Clamping value %.12g -> %.12g to %u", (double)*edata, (double)dd, (unsigned int)DTYPE_MAX);
                dd = DTYPE_MAX;
            }
            if (dd < DTYPE_MIN) {
                WARNING("Clamping value %.12g -> %.12g to %u.\n", (double)*edata, (double)dd, (unsigned int)DTYPE_MIN);
                dd = DTYPE_MIN;
            }
            // Right place for this?  Not really....
            if (!ETYPE_INTEGER) {
                // to avoid compiler warnings about int types, even though this will never happen at runtime.
                double ddd = (double)dd;
                // NaN and Inf detection...
                if (!isfinite(ddd) || isnan(ddd)) {
                    WARNING("Replacing inf/nan value (element %i,%i) = %g with %g\n", i, d, ddd, (double)DTYPE_MAX);
                    dd = DTYPE_MAX;
                }
            }
            *ddata = (dtype)dd;

            ddata++;
            edata++;
        }
    }

    // tune up kd->maxval so that it is larger than the external-type
    // value of all points in the tree.
    for (d=0; d<D; d++) {
        dtype dmax = POINT_ED(kd, d, kd->maxval[d], KD_ROUND);
        etype emax = POINT_DE(kd, d, dmax);
        kd->maxval[d] = MAX(kd->maxval[d], emax);
    }
#ifndef NDEBUG
    for (i=0; i<N; i++) {
        for (d=0; d<D; d++) {
            etype e = POINT_DE(kd, d, KD_DATA(kd, D, i)[d]);
            assert(e <= kd->maxval[d]);
            assert(e >= kd->minval[d]);
        }
    }
#endif
}

static void compute_bb(const dtype* data, int D, int N, dtype* lo, dtype* hi) {
    int d, i;

    for (d=0; d<D; d++) {
        hi[d] = DTYPE_MIN;
        lo[d] = DTYPE_MAX;
    }
    /* (since data is stored lexicographically we can just iterate through it) */
    /* (avoid doing kd->data[NODE(i)*D + d] many times; just ++ the pointer) */
    for (i=0; i<N; i++) {
        for (d=0; d<D; d++) {
            if (*data > hi[d]) hi[d] = *data;
            if (*data < lo[d]) lo[d] = *data;
            data++;
        }
    }
}

static void save_bb(kdtree_t* kd, int i, const dtype* lo, const dtype* hi) {
    int D = kd->ndim;
    int d;
    for (d=0; d<D; d++) {
        (LOW_HR (kd, D, i))[d] = POINT_DT(kd, d, lo[d], floor);
        (HIGH_HR(kd, D, i))[d] = POINT_DT(kd, d, hi[d], ceil);
    }
}

static int needs_data_conversion() {
    return DTYPE_INTEGER && !ETYPE_INTEGER;
}

kdtree_t* MANGLE(kdtree_build_2)
     (kdtree_t* kd, etype* indata, int N, int D, int Nleaf, int treetype, unsigned int options, double* minval, double* maxval) {
    int i;
    int xx;
    int lnext, level;
    int maxlevel;
    dtype hi[D], lo[D];
    dtype* data = NULL;

    dtype nullbb[D];
    for (i=0; i<D; i++)
        nullbb[i] = 0;

    maxlevel = kdtree_compute_levels(N, Nleaf);

    assert(maxlevel > 0);
    assert(D <= KDTREE_MAX_DIM);
#if defined(KD_DIM)
    assert(D == KD_DIM);
    // let the compiler know that D is a constant...
    D = KD_DIM;
#endif

    /* Parameters checking */
    if (!indata || !N || !D) {
        ERROR("Data, N, or D is zero");
        return NULL;
    }
    /* Make sure we have enough data */
    if ((1 << maxlevel) - 1 > N) {
        ERROR("Too few data points for number of tree levels (%i < %i)!", N, (1 << maxlevel) - 1);
        return NULL;
    }

    if ((options & KD_BUILD_SPLIT) &&
        (options & KD_BUILD_NO_LR) &&
        !(options & KD_BUILD_SPLITDIM) &&
        TTYPE_INTEGER) {
        ERROR("Currently you can't set KD_BUILD_NO_LR for int trees "
              "unless you also set KD_BUILD_SPLITDIM");
        return NULL;
    }

    if (!kd)
        kd = kdtree_new(N, D, Nleaf);

    kd->treetype = treetype;
    if (minval) {
        kd->minval = CALLOC(D, sizeof(double));
        memcpy(kd->minval, minval, D*sizeof(double));
    }
    if (maxval) {
        kd->maxval = CALLOC(D, sizeof(double));
        memcpy(kd->maxval, maxval, D*sizeof(double));
    }

    if (!kd->data.any) {
        if (needs_data_conversion()) {
            // need to convert the data from "etype" to "dtype".
            convert_data(kd, indata, N, D, Nleaf);
        } else {
            kd->data.any = indata;

            // ???
            if (!ETYPE_INTEGER) {
                int i,d;
                etype* edata = indata;
                for (i=0; i<N; i++) {
                    for (d=0; d<D; d++) {
                        etype dd = edata[i*D + d];
                        // to avoid compiler warnings about int types, even though this will never happen at runtime.
                        double ddd = (double)dd;
                        // NaN and Inf detection...
                        if (!isfinite(ddd) || isnan(ddd)) {
                            WARNING("Replacing inf/nan value (element %i,%i) = %g with %g\n", i, d, (double)dd, (double)DTYPE_MAX);
                            edata[i*D + d] = DTYPE_MAX;
                        }
                    }
                }
            }
        }
    }
    if (needs_data_conversion()) {
        // compute scaling params
        if (!kd->minval || !kd->maxval) {
            free(kd->minval);
            free(kd->maxval);
            kd->minval = MALLOC((size_t)D * sizeof(double));
            kd->maxval = MALLOC((size_t)D * sizeof(double));
            assert(kd->minval);
            assert(kd->maxval);
            kd->scale = compute_scale(kd->data.any, N, D, kd->minval, kd->maxval);
        } else {
            // limits were pre-set by the user.  just compute scale.
            double range;
            range = maxrange(kd->minval, kd->maxval, D);
            kd->scale = (double)TTYPE_MAX / range;
        }
        kd->invscale = 1.0 / kd->scale;
    }

    /* perm stores the permutation indexes. This gets shuffled around during
     * sorts to keep track of the original index. */
    kd->perm = MALLOC((size_t)N * sizeof(u32));
    assert(kd->perm);
    for (i = 0; i < N; i++)
        kd->perm[i] = i;

    kd->lr = MALLOC((size_t)kd->nbottom * sizeof(int32_t));
    assert(kd->lr);

    if (options & KD_BUILD_BBOX) {
        kd->bb.any = MALLOC((size_t)kd->nnodes * (size_t)2 * (size_t)D * sizeof(ttype));
        kd->n_bb = kd->nnodes;
        assert(kd->bb.any);
    }
    if (options & KD_BUILD_SPLIT) {
        kd->split.any = MALLOC((size_t)kd->ninterior * sizeof(ttype));
        assert(kd->split.any);
    }
    if (((options & KD_BUILD_SPLIT) && !TTYPE_INTEGER) ||
        (options & KD_BUILD_SPLITDIM)) {
        kd->splitdim = MALLOC((size_t)kd->ninterior * sizeof(u8));
        kd->splitmask = UINT32_MAX;
        kd->dimmask = 0;
    } else if (options & KD_BUILD_SPLIT)
        compute_splitbits(kd);

    if (options & KD_BUILD_LINEAR_LR)
        kd->has_linear_lr = TRUE;

    /* Use the lr array as a stack while building. In place in your face! */
    kd->lr[0] = N - 1;
    lnext = 1;
    level = 0;

    // shorthand
    data = kd->data.DTYPE;

    /* And in one shot, make the kdtree. Because the lr pointers
     * are only stored for the bottom layer, we use the lr array as a
     * stack. At finish, it contains the r pointers for the bottom nodes.
     * The l pointer is simply +1 of the previous right pointer, or 0 if we
     * are at the first element of the lr array. */
    for (i = 0; i < kd->ninterior; i++) {
        unsigned int d;
        int left, right;
        dtype maxrange;
        ttype s;
        unsigned int c;
        int dim = 0;
        int m;
        dtype qsplit = 0;

        /* Have we reached the next level in the tree? */
        if (i == lnext) {
            level++;
            lnext = lnext * 2 + 1;
        }

        /* Since we're not storing the L pointers, we have to infer L */
        if (i == (1<<level)-1) {
            left = 0;
        } else {
            left = kd->lr[i-1] + 1;
        }
        right = kd->lr[i];

        assert(right != (unsigned int)-1);

        if (left >= right) {
            //debug("Empty node %i: left=right=%i\n", i, left);
            if (options & KD_BUILD_BBOX)
                save_bb(kd, i, nullbb, nullbb);
            if (kd->splitdim)
                kd->splitdim[i] = 0;
            c = 2*i;
            if (level == maxlevel - 2)
                c -= kd->ninterior;
            kd->lr[c+1] = right;
            kd->lr[c+2] = right;
            continue;
        }

        /* More sanity */
        assert(0 <= left);
        assert(left <= right);
        assert(right < N);

        /* Find the bounding-box for this node. */
        compute_bb(KD_DATA(kd, D, left), D, right - left + 1, lo, hi);

        if (options & KD_BUILD_BBOX)
            save_bb(kd, i, lo, hi);

        /* Split along dimension with largest range */
        maxrange = DTYPE_MIN;
        for (d=0; d<D; d++)
            if ((hi[d] - lo[d]) >= maxrange) {
                maxrange = hi[d] - lo[d];
                dim = d;
            }
        d = dim;
        assert (d < D);

        if ((options & KD_BUILD_FORCE_SORT) ||
            (TTYPE_INTEGER && !(options & KD_BUILD_SPLITDIM))) {

            /* We're packing dimension and split location into an int. */

            /* Sort the data. */

            /* Because the nature of the inttree is to bin the split
             * planes, we have to be careful. Here, we MUST sort instead
             * of merely partitioning, because we may not be able to
             * properly represent the median as a split plane. Imagine the
             * following on the dtype line:
             *
             *    |P P   | P M  | P    |P     |  PP |  ------> X
             *           1      2
             * The |'s are possible split positions. If M is selected to
             * split on, we actually cannot select the split 1 or 2
             * immediately, because if we selected 2, then M would be on
             * the wrong side (the medians always go to the right) and we
             * can't select 1 because then P would be on the wrong side.
             * So, the solution is to try split 2, and if point M-1 is on
             * the correct side, great. Otherwise, we have to move shift
             * point M-1 into the right side and only then chose plane 1. */


            /* FIXME but qsort allocates a 2nd perm array GAH */
            if (kdtree_qsort(data, kd->perm, left, right, D, dim)) {
                ERROR("kdtree_qsort failed");
                // FIXME: memleak mania!
                return NULL;
            }
            m = (1 + (size_t)left + (size_t)right)/2;
            assert(m >= 0);
            assert(m >= left);
            assert(m <= right);

            /* Make sure sort works */
            for(xx=left; xx<=right-1; xx++) {
                assert(KD_ARRAY_VAL(data, D, xx,   d) <=
                       KD_ARRAY_VAL(data, D, xx+1, d));
            }

            /* Encode split dimension and value. */
            /* "s" is the location of the splitting plane in the "tree"
             data type. */
            s = POINT_DT(kd, d, KD_ARRAY_VAL(data, D, m, d), KD_ROUND);

            if (kd->split.any) {
                /* If we are using the "split" array to store both the
                 splitting plane and the splitting dimension, then we
                 truncate a few bits from "s" here. */
                bigint tmps = s;
                tmps &= kd->splitmask;
                assert((tmps & kd->dimmask) == 0);
                s = tmps;
            }
            /* "qsplit" is the location of the splitting plane in the "data"
             type. */
            qsplit = POINT_TD(kd, d, s);

            /* Play games to make sure we properly partition the data */
            while (m < right && KD_ARRAY_VAL(data, D, m, d) < qsplit) m++;
            while (left < m  && qsplit < KD_ARRAY_VAL(data, D, m-1, d)) m--;

            /* Even more sanity */
            assert(m >= -1);
            assert(left <= m);
            assert(m <= right);
            for (xx=left; m && xx<=m-1; xx++)
                assert(KD_ARRAY_VAL(data, D, xx, d) <= qsplit);
            for (xx=m; xx<=right; xx++)
                assert(qsplit <= KD_ARRAY_VAL(data, D, xx, d));

        } else {
            /* "m-1" becomes R of the left child;
             "m" becomes L of the right child. */
            if (kd->has_linear_lr) {
                m = kdtree_left(kd, KD_CHILD_RIGHT(i));
            } else {
                /* Pivot the data at the median */
                m = (1 + (size_t)left + (size_t)right) / 2;
            }
            assert(m >= 0);
            assert(m >= left);
            assert(m <= right);
            kdtree_quickselect_partition(data, kd->perm, left, right, D, dim, m);

            s = POINT_DT(kd, d, KD_ARRAY_VAL(data, D, m, d), KD_ROUND);

            assert(m != 0);
            assert(left <= (m-1));
            assert(m <= right);
            for (xx=left; xx<=m-1; xx++)
                assert(KD_ARRAY_VAL(data, D, xx, d) <=
                       KD_ARRAY_VAL(data, D, m, d));
            for (xx=left; xx<=m-1; xx++)
                assert(KD_ARRAY_VAL(data, D, xx, d) <= s);
            for (xx=m; xx<=right; xx++)
                assert(KD_ARRAY_VAL(data, D, m, d) <=
                       KD_ARRAY_VAL(data, D, xx, d));
            for (xx=m; xx<=right; xx++)
                assert(s <= KD_ARRAY_VAL(data, D, xx, d));
        }

        if (kd->split.any) {
            if (kd->splitdim)
                *KD_SPLIT(kd, i) = s;
            else {
                bigint tmps = s;
                *KD_SPLIT(kd, i) = tmps | dim;
            }
        }
        if (kd->splitdim)
            kd->splitdim[i] = dim;

        /* Store the R pointers for each child */
        c = 2*i;
        if (level == maxlevel - 2)
            c -= kd->ninterior;

        kd->lr[c+1] = m-1;
        kd->lr[c+2] = right;

        assert(c+2 < kd->nbottom);
    }

    for (i=0; i<kd->nbottom-1; i++)
        assert(kd->lr[i] <= kd->lr[i+1]);

    if (options & KD_BUILD_BBOX) {
        // Compute bounding boxes for leaf nodes.
        int L, R = -1;
        for (i=0; i<kd->nbottom; i++) {
            L = R + 1;
            R = kd->lr[i];
            assert(L == kdtree_leaf_left(kd, i + kd->ninterior));
            assert(R == kdtree_leaf_right(kd, i + kd->ninterior));
            compute_bb(KD_DATA(kd, D, L), D, R - L + 1, lo, hi);
            save_bb(kd, i + kd->ninterior, lo, hi);
        }

        // check that it worked...
#ifndef NDEBUG
        for (i=0; i<kd->nbottom; i++) {
            ttype* lo;
            ttype* hi;
            int j, d;
            int nodeid = i + kd->ninterior;
            lo = LOW_HR(kd, kd->ndim, nodeid);
            hi = HIGH_HR(kd, kd->ndim, nodeid);
            for (j=kdtree_leaf_left(kd, nodeid); j<=kdtree_leaf_right(kd, nodeid); j++) {
                dtype* data = KD_DATA(kd, kd->ndim, j);
                for (d=0; d<kd->ndim; d++) {
                    assert(POINT_TD(kd, d, lo[d]) <= data[d]);
                    assert(POINT_TD(kd, d, hi[d]) >= data[d]);
                    assert(POINT_DT(kd, d, data[d], KD_ROUND) <= hi[d]);
                    assert(POINT_DT(kd, d, data[d], KD_ROUND) >= lo[d]);
                }
            }
        }
#endif

    }

    if (options & KD_BUILD_NO_LR) {
        FREE(kd->lr);
        kd->lr = NULL;
    }

    // set function table pointers.
    MANGLE(kdtree_update_funcs)(kd);

    return kd;
}

void MANGLE(kdtree_fix_bounding_boxes)(kdtree_t* kd) {
    // FIXME - do this log(N) times more efficiently by propagating
    // bounding boxes up the levels of the tree...
    int i;
    int D = kd->ndim;
    kd->bb.any = MALLOC((size_t)kd->nnodes * sizeof(ttype) * (size_t)(D * 2));
    assert(kd->bb.any);
    for (i=0; i<kd->nnodes; i++) {
        unsigned int left, right;
        dtype hi[D], lo[D];
        left = kdtree_left(kd, i);
        right = kdtree_right(kd, i);
        compute_bb(KD_DATA(kd, D, left), D, right - left + 1, lo, hi);
        save_bb(kd, i, lo, hi);
    }
}

double MANGLE(kdtree_node_point_mindist2)
     (const kdtree_t* kd, int node, const etype* query) {
    int D = kd->ndim;
    int d;
    ttype* tlo, *thi;
    double d2 = 0.0;
    if (!bboxes(kd, node, &tlo, &thi, D)) {
        ERROR("Error: kdtree does not have bounding boxes!");
        return LARGE_VAL;
    }
    for (d=0; d<D; d++) {
        etype delta;
        etype lo = POINT_TE(kd, d, tlo[d]);
        if (query[d] < lo)
            delta = lo - query[d];
        else {
            etype hi = POINT_TE(kd, d, thi[d]);
            if (query[d] > hi)
                delta = query[d] - hi;
            else
                continue;
        }
        d2 += delta * delta;
    }
    return d2;
}

double MANGLE(kdtree_node_point_maxdist2)
     (const kdtree_t* kd, int node, const etype* query) {
    int D = kd->ndim;
    int d;
    ttype* tlo=NULL, *thi=NULL;
    double d2 = 0.0;
    if (!bboxes(kd, node, &tlo, &thi, D)) {
        ERROR("Error: kdtree_node_point_maxdist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }
    for (d=0; d<D; d++) {
        etype delta1, delta2, delta;
        etype lo = POINT_TE(kd, d, tlo[d]);
        etype hi = POINT_TE(kd, d, thi[d]);
        if (query[d] < lo)
            delta = hi - query[d];
        else if (query[d] > hi)
            delta = query[d] - lo;
        else {
            delta1 = hi - query[d];
            delta2 = query[d] - lo;
            delta = MAX(delta1, delta2);
        }
        d2 += delta*delta;
    }
    return d2;
}

anbool MANGLE(kdtree_node_point_mindist2_exceeds)
     (const kdtree_t* kd, int node, const etype* query, double maxd2) {
    int D = kd->ndim;
    int d;
    ttype* tlo, *thi;
    double d2 = 0.0;

    if (!bboxes(kd, node, &tlo, &thi, D)) {
        //ERROR("Error: kdtree does not have bounding boxes!");
        return FALSE;
    }
    for (d=0; d<D; d++) {
        etype delta;
        etype lo = POINT_TE(kd, d, tlo[d]);
        if (query[d] < lo)
            delta = lo - query[d];
        else {
            etype hi = POINT_TE(kd, d, thi[d]);
            if (query[d] > hi)
                delta = query[d] - hi;
            else
                continue;
        }
        d2 += delta * delta;
        if (d2 > maxd2)
            return TRUE;
    }
    return FALSE;
}

anbool MANGLE(kdtree_node_point_maxdist2_exceeds)
     (const kdtree_t* kd, int node, const etype* query, double maxd2) {
    int D = kd->ndim;
    int d;
    ttype* tlo=NULL, *thi=NULL;
    double d2 = 0.0;

    if (!bboxes(kd, node, &tlo, &thi, D)) {
        ERROR("Error: kdtree_node_point_maxdist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }
    for (d=0; d<D; d++) {
        etype delta1, delta2, delta;
        etype lo = POINT_TE(kd, d, tlo[d]);
        etype hi = POINT_TE(kd, d, thi[d]);
        if (query[d] < lo)
            delta = hi - query[d];
        else if (query[d] > hi)
            delta = query[d] - lo;
        else {
            delta1 = hi - query[d];
            delta2 = query[d] - lo;
            delta = (delta1 > delta2 ? delta1 : delta2);
        }
        d2 += delta*delta;
        if (d2 > maxd2)
            return TRUE;
    }
    return FALSE;
}

double MANGLE(kdtree_node_node_maxdist2)
     (const kdtree_t* kd1, int node1,
      const kdtree_t* kd2, int node2) {
    ttype *tlo1=NULL, *tlo2=NULL, *thi1=NULL, *thi2=NULL;
    double d2 = 0.0;
    int d, D = kd1->ndim;

    assert(kd1->ndim == kd2->ndim);
    if (!bboxes(kd1, node1, &tlo1, &thi1, D)) {
        ERROR("Error: kdtree_node_node_maxdist2: kdtree does not have bounding boxes!");
        return FALSE;
    }
    if (!bboxes(kd2, node2, &tlo2, &thi2, D)) {
        ERROR("Error: kdtree_node_node_maxdist2: kdtree does not have bounding boxes!");
        return FALSE;
    }
    // Since the two trees can have different conversion factors,
    // we have to convert both to the external type.
    // FIXME - we do assume that POINT_TE works for both of them --
    // ie, ~we assume they are the same treetype.
    for (d=0; d<D; d++) {
        etype alo, ahi, blo, bhi;
        etype delta1, delta2, delta;
        alo = POINT_TE(kd1, d, tlo1[d]);
        ahi = POINT_TE(kd1, d, thi1[d]);
        blo = POINT_TE(kd2, d, tlo2[d]);
        bhi = POINT_TE(kd2, d, thi2[d]);
        if (ETYPE_INTEGER)
            WARNING("HACK - int overflow is possible here.");
        delta1 = bhi - alo;
        delta2 = ahi - blo;
        delta = MAX(delta1, delta2);
        d2 += delta*delta;
    }
    return d2;
}

double MANGLE(kdtree_node_node_mindist2)
     (const kdtree_t* kd1, int node1,
      const kdtree_t* kd2, int node2) {
    ttype *tlo1=NULL, *tlo2=NULL, *thi1=NULL, *thi2=NULL;
    double d2 = 0.0;
    int d, D = kd1->ndim;
    assert(kd1->ndim == kd2->ndim);
    if (!bboxes(kd1, node1, &tlo1, &thi1, D)) {
        ERROR("Error: kdtree_node_node_mindist2: kdtree does not have bounding boxes!");
        return FALSE;
    }
    if (!bboxes(kd2, node2, &tlo2, &thi2, D)) {
        ERROR("Error: kdtree_node_node_mindist2: kdtree does not have bounding boxes!");
        return FALSE;
    }
    for (d=0; d<D; d++) {
        etype alo, ahi, blo, bhi;
        etype delta;
        ahi = POINT_TE(kd1, d, thi1[d]);
        blo = POINT_TE(kd2, d, tlo2[d]);
        if (ahi < blo)
            delta = blo - ahi;
        else {
            alo = POINT_TE(kd1, d, tlo1[d]);
            bhi = POINT_TE(kd2, d, thi2[d]);
            if (bhi < alo)
                delta = alo - bhi;
            else
                continue;
        }
        d2 += delta*delta;
    }
    return d2;
}

anbool MANGLE(kdtree_node_node_maxdist2_exceeds)
     (const kdtree_t* kd1, int node1,
      const kdtree_t* kd2, int node2,
      double maxd2) {
    ttype *tlo1=NULL, *tlo2=NULL, *thi1=NULL, *thi2=NULL;
    double d2 = 0.0;
    int d, D = kd1->ndim;

    //assert(kd1->treetype == kd2->treetype);
    assert(kd1->ndim == kd2->ndim);

    if (!bboxes(kd1, node1, &tlo1, &thi1, D)) {
        ERROR("Error: kdtree_node_node_maxdist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }

    if (!bboxes(kd2, node2, &tlo2, &thi2, D)) {
        ERROR("Error: kdtree_node_node_maxdist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }

    for (d=0; d<D; d++) {
        etype alo, ahi, blo, bhi;
        etype delta1, delta2, delta;
        alo = POINT_TE(kd1, d, tlo1[d]);
        ahi = POINT_TE(kd1, d, thi1[d]);
        blo = POINT_TE(kd2, d, tlo2[d]);
        bhi = POINT_TE(kd2, d, thi2[d]);
        // HACK - if etype is integer...
        if (ETYPE_INTEGER)
            WARNING("HACK - int overflow is possible here.");
        delta1 = bhi - alo;
        delta2 = ahi - blo;
        delta = (delta1 > delta2 ? delta1 : delta2);
        d2 += delta*delta;
        if (d2 > maxd2)
            return TRUE;
    }
    return FALSE;
}

anbool MANGLE(kdtree_node_node_mindist2_exceeds)
     (const kdtree_t* kd1, int node1,
      const kdtree_t* kd2, int node2,
      double maxd2) {
    ttype *tlo1=NULL, *tlo2=NULL, *thi1=NULL, *thi2=NULL;
    double d2 = 0.0;
    int d, D = kd1->ndim;

    //assert(kd1->treetype == kd2->treetype);
    assert(kd1->ndim == kd2->ndim);

    if (!bboxes(kd1, node1, &tlo1, &thi1, D)) {
        //ERROR("Error: kdtree_node_node_mindist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }

    if (!bboxes(kd2, node2, &tlo2, &thi2, D)) {
        //ERROR("Error: kdtree_node_node_mindist2_exceeds: kdtree does not have bounding boxes!");
        return FALSE;
    }

    for (d=0; d<D; d++) {
        etype alo, ahi, blo, bhi;
        etype delta;
        ahi = POINT_TE(kd1, d, thi1[d]);
        blo = POINT_TE(kd2, d, tlo2[d]);
        if (ahi < blo)
            delta = blo - ahi;
        else {
            alo = POINT_TE(kd1, d, tlo1[d]);
            bhi = POINT_TE(kd2, d, thi2[d]);
            if (bhi < alo)
                delta = alo - bhi;
            else
                continue;
        }
        d2 += delta*delta;
        if (d2 > maxd2)
            return TRUE;
    }
    return FALSE;
}

static anbool do_boxes_overlap(const ttype* lo1, const ttype* hi1,
                               const ttype* lo2, const ttype* hi2, int D) {
    int d;
    for (d=0; d<D; d++) {
        if (lo1[d] > hi2[d])
            return FALSE;
        if (lo2[d] > hi1[d])
            return FALSE;
    }
    return TRUE;
}

/* Is the first box contained within the second? */
static anbool is_box_contained(const ttype* lo1, const ttype* hi1,
                               const ttype* lo2, const ttype* hi2, int D) {
    int d;
    for (d=0; d<D; d++) {
        if (lo1[d] < lo2[d])
            return FALSE;
        if (hi1[d] > hi2[d])
            return FALSE;
    }
    return TRUE;
}

static void nodes_contained_rec(const kdtree_t* kd,
                                int nodeid,
                                const ttype* qlo, const ttype* qhi,
                                void (*cb_contained)(const kdtree_t* kd, int node, void* extra),
                                void (*cb_overlap)(const kdtree_t* kd, int node, void* extra),
                                void* cb_extra) {
    ttype *tlo=NULL, *thi=NULL;
    int D = kd->ndim;

    // leaf nodes don't have bounding boxes, so we have to do this check first!
    if (KD_IS_LEAF(kd, nodeid)) {
        cb_overlap(kd, nodeid, cb_extra);
        return;
    }

    if (!bboxes(kd, nodeid, &tlo, &thi, D)) {
        ERROR("Error: kdtree_nodes_contained: node %i doesn't have a bounding box", nodeid);
        return;
    }

    if (!do_boxes_overlap(tlo, thi, qlo, qhi, D))
        return;

    if (is_box_contained(tlo, thi, qlo, qhi, D)) {
        cb_contained(kd, nodeid, cb_extra);
        return;
    }

    nodes_contained_rec(kd,  KD_CHILD_LEFT(nodeid), qlo, qhi,
                        cb_contained, cb_overlap, cb_extra);
    nodes_contained_rec(kd, KD_CHILD_RIGHT(nodeid), qlo, qhi,
                        cb_contained, cb_overlap, cb_extra);
}

void MANGLE(kdtree_nodes_contained)
     (const kdtree_t* kd,
      const void* vquerylow, const void* vqueryhi,
      void (*cb_contained)(const kdtree_t* kd, int node, void* extra),
      void (*cb_overlap)(const kdtree_t* kd, int node, void* extra),
      void* cb_extra) {
    int D = kd->ndim;
    int d;
    ttype qlo[D], qhi[D];
    const etype* querylow = vquerylow;
    const etype* queryhi = vqueryhi;

    for (d=0; d<D; d++) {
        double q;
        qlo[d] = q = POINT_ET(kd, d, querylow[d], floor);
        if (q < TTYPE_MIN) {
            //WARNING("Error: query value %g is below the minimum range of the tree %g.\n", q, (double)TTYPE_MIN);
            qlo[d] = TTYPE_MIN;
        } else if (q > TTYPE_MAX) {
            // query's low position is more than the tree's max: no overlap is possible.
            return;
        }
        qhi[d] = q = POINT_ET(kd, d, queryhi [d], ceil );
        if (q > TTYPE_MAX) {
            //WARNING("Error: query value %g is above the maximum range of the tree %g.\n", q, (double)TTYPE_MAX);
            qhi[d] = TTYPE_MAX;
        } else if (q < TTYPE_MIN) {
            // query's high position is less than the tree's min: no overlap is possible.
            return;
        }
    }

    nodes_contained_rec(kd, 0, qlo, qhi, cb_contained, cb_overlap, cb_extra);
}

int MANGLE(kdtree_get_bboxes)(const kdtree_t* kd, int node,
                              void* vbblo, void* vbbhi) {
    etype* bblo = vbblo;
    etype* bbhi = vbbhi;
    ttype *tlo=NULL, *thi=NULL;
    int D = kd->ndim;
    int d;

    if (!bboxes(kd, node, &tlo, &thi, D))
        return FALSE;

    for (d=0; d<D; d++) {
        bblo[d] = POINT_TE(kd, d, tlo[d]);
        bbhi[d] = POINT_TE(kd, d, thi[d]);
    }
    return TRUE;
}

void MANGLE(kdtree_update_funcs)(kdtree_t* kd) {
    kd->fun.get_data = get_data;
    kd->fun.copy_data_double = copy_data_double;
    kd->fun.get_splitval = MANGLE(kdtree_get_splitval);
    kd->fun.get_bboxes = MANGLE(kdtree_get_bboxes);
    kd->fun.check = MANGLE(kdtree_check);
    kd->fun.fix_bounding_boxes = MANGLE(kdtree_fix_bounding_boxes);
    kd->fun.nearest_neighbour_internal = MANGLE(kdtree_nn);
    kd->fun.rangesearch = MANGLE(kdtree_rangesearch_options);
    kd->fun.rangesearch_continuation_init =
        MANGLE(kdtree_rangesearch_continuation_init);
    kd->fun.rangesearch_continuation_step =
        MANGLE(kdtree_rangesearch_continuation_step);
    kd->fun.rangesearch_continuation_finish =
        MANGLE(kdtree_rangesearch_continuation_finish);
    kd->fun.rangesearch_continuation_cleanup =
        MANGLE(kdtree_rangesearch_continuation_cleanup);
    kd->fun.nodes_contained = MANGLE(kdtree_nodes_contained);
}

static inline anbool MANGLE(kdtree_prefetch_valid_node)
     (const kdtree_t *kd,
      int nodeid) {
  if (!kd)
    return FALSE;

  if (nodeid < 0)
    return FALSE;

  /*
   * Astrometry.net KD nodes are used as a packed binary tree:
   *
   *   interior nodes: 0 ... kd->ninterior - 1
   *   bottom leaves:  kd->ninterior ... kd->ninterior + kd->nbottom - 1
   *
   * Do not use kd->nnodes alone as the structural validity test for product
   * traversal, because kdtree_leaf_left/right() index kd->lr by:
   *
   *   leafid = nodeid - kd->ninterior
   */
  if (nodeid < kd->ninterior)
    return TRUE;

  if (nodeid < kd->ninterior + kd->nbottom)
    return TRUE;

  return FALSE;
}

static inline anbool MANGLE(kdtree_prefetch_valid_leaf_node)
     (const kdtree_t *kd,
      int nodeid) {
  if (!kd)
    return FALSE;

  if (nodeid < kd->ninterior)
    return FALSE;

  if (nodeid >= kd->ninterior + kd->nbottom)
    return FALSE;

  return TRUE;
}
static inline anbool MANGLE(kdtree_prefetch_valid_data_range)
     (const kdtree_t *kd,
      int L,
      int R) {
  if (!kd)
    return FALSE;

  if (L < 0)
    return FALSE;

  if (R < L)
    return FALSE;

  if (R >= kd->ndata)
    return FALSE;

  return TRUE;
}

static inline int MANGLE(kdtree_prefetch_push_node)
     (const kdtree_t *kd,
      int *nodestack,
      int *stackpos,
      int nodeid) {
    if (!MANGLE(kdtree_prefetch_valid_node)(kd, nodeid)) {
            fprintf(stderr,
                    "[kd-prefetch] refusing invalid child node=%i "
                    "nnodes=%i stackpos=%i\n",
                    nodeid, kd ? kd->nnodes : -1, stackpos ? *stackpos : -999);
            return -1;
        }

    if ((*stackpos + 1) >= KDTREE_PREFETCH_STACK_MAX) {
        fprintf(stderr,
                "[kd-prefetch] subtree stack overflow next=%i max=%i "
                "node=%i nnodes=%i\n",
                *stackpos + 1, KDTREE_PREFETCH_STACK_MAX,
                nodeid, kd ? kd->nnodes : -1);
        return -1;
    }

    (*stackpos)++;
    nodestack[*stackpos] = nodeid;

    return 0;
}

static anbool MANGLE(kdtree_prefetch_size_mul)
     (size_t left,
      size_t right,
      size_t *result) {
    if (!result) {
        return FALSE;
    }

    if (left && right > SIZE_MAX / left) {
        return FALSE;
    }

    *result = left * right;
    return TRUE;
}

static int MANGLE(kdtree_prefetch_emit_hint)
     (const kdtree_t *kd,
      const void *address,
      size_t length,
      kdtree_prefetch_array_kind_t kind,
      unsigned int priority,
      const kdtree_prefetch_sink_t *sink) {
    kdtree_prefetch_hint_t hint;
    int status;

    if (!kd ||
        !kd->io ||
        !address ||
        !length ||
        !sink ||
        !sink->emit) {
        return KDTREE_PREFETCH_EMIT_ERROR;
    }

    memset(&hint, 0, sizeof(hint));

    hint.mapping = kd->io;
    hint.address = address;
    hint.length = length;
    hint.kind = kind;
    hint.priority = priority;

    status = sink->emit(sink->userdata, &hint);
    if (status < 0) {
        return KDTREE_PREFETCH_EMIT_ERROR;
    }
    if (status > 0) {
        return KDTREE_PREFETCH_EMIT_REFUSED;
    }
    return KDTREE_PREFETCH_EMIT_CONTINUE;
}

static int MANGLE(kdtree_prefetch_emit_leaf_metadata)
     (const kdtree_t *kd,
      int nodeid,
      const kdtree_prefetch_sink_t *sink) {
    int first_lr;
    int leafid;
    int lr_count;
    size_t lr_bytes;

    if (!kd || !kd->lr || kd->has_linear_lr ||
        !MANGLE(kdtree_prefetch_valid_leaf_node)(kd, nodeid)) {
        return KDTREE_PREFETCH_EMIT_CONTINUE;
    }

    leafid = nodeid - kd->ninterior;
    first_lr = leafid > 0 ? leafid - 1 : leafid;
    lr_count = leafid > 0 ? 2 : 1;
    if (first_lr < 0 || first_lr >= kd->nbottom ||
        lr_count <= 0 || lr_count > kd->nbottom - first_lr ||
        !MANGLE(kdtree_prefetch_size_mul)
            ((size_t)lr_count, sizeof(*kd->lr), &lr_bytes)) {
        return KDTREE_PREFETCH_EMIT_CONTINUE;
    }

    return MANGLE(kdtree_prefetch_emit_hint)
        (kd,
         kd->lr + first_lr,
         lr_bytes,
         KDTREE_PREFETCH_ARRAY_LR,
         KDTREE_PREFETCH_PRIORITY_METADATA,
         sink);
}

static int MANGLE(kdtree_prefetch_emit_payload)
     (const kdtree_t *kd,
      int D,
      int L,
      int R,
      const kdtree_prefetch_sink_t *sink) {
    size_t count;
    size_t coordinates;
    size_t data_bytes;
    size_t perm_bytes;

    if (!MANGLE(kdtree_prefetch_valid_data_range)(kd, L, R) ||
        D <= 0) {
        return KDTREE_PREFETCH_EMIT_ERROR;
    }

    count = (size_t)(R - L + 1);
    if (kd->data.any) {
        if (!MANGLE(kdtree_prefetch_size_mul)
                (count, (size_t)D, &coordinates) ||
            !MANGLE(kdtree_prefetch_size_mul)
                (coordinates, sizeof(dtype), &data_bytes)) {
            return KDTREE_PREFETCH_EMIT_ERROR;
        }
        int emit_status = MANGLE(kdtree_prefetch_emit_hint)
            (kd,
             KD_DATA(kd, D, L),
             data_bytes,
             KDTREE_PREFETCH_ARRAY_DATA,
             KDTREE_PREFETCH_PRIORITY_LEAF,
             sink);
        if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
            return emit_status;
        }
    }

    if (kd->perm) {
        if (!MANGLE(kdtree_prefetch_size_mul)
                (count, sizeof(*kd->perm), &perm_bytes)) {
            return KDTREE_PREFETCH_EMIT_ERROR;
        }
        int emit_status = MANGLE(kdtree_prefetch_emit_hint)
            (kd,
             kd->perm + L,
             perm_bytes,
             KDTREE_PREFETCH_ARRAY_PERM,
             KDTREE_PREFETCH_PRIORITY_LEAF,
             sink);
        if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
            return emit_status;
        }
    }
    return KDTREE_PREFETCH_EMIT_CONTINUE;
}

static int MANGLE(kdtree_prefetch_emit_split_metadata)
     (const kdtree_t *kd,
      int nodeid,
      const kdtree_prefetch_sink_t *sink) {
    int emit_status;

    if (!kd || !kd->split.any || nodeid < 0 ||
        nodeid >= kd->ninterior) {
        return KDTREE_PREFETCH_EMIT_CONTINUE;
    }

    emit_status = MANGLE(kdtree_prefetch_emit_hint)
        (kd,
         KD_SPLIT(kd, nodeid),
         sizeof(ttype),
         KDTREE_PREFETCH_ARRAY_SPLIT,
         KDTREE_PREFETCH_PRIORITY_METADATA,
         sink);
    if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
        return emit_status;
    }
    if (kd->splitdim) {
        emit_status = MANGLE(kdtree_prefetch_emit_hint)
            (kd,
             kd->splitdim + nodeid,
             sizeof(*kd->splitdim),
             KDTREE_PREFETCH_ARRAY_SPLITDIM,
             KDTREE_PREFETCH_PRIORITY_METADATA,
             sink);
    }
    return emit_status;
}

static int MANGLE(kdtree_prefetch_emit_bbox_metadata)
     (const kdtree_t *kd,
      int D,
      int nodeid,
      const kdtree_prefetch_sink_t *sink) {
    size_t coordinates;
    size_t bbox_bytes;

    if (!kd || !kd->bb.any || D <= 0 || nodeid < 0 ||
        nodeid >= kd->ninterior ||
        !MANGLE(kdtree_prefetch_size_mul)
            ((size_t)D, 2U, &coordinates) ||
        !MANGLE(kdtree_prefetch_size_mul)
            (coordinates, sizeof(ttype), &bbox_bytes)) {
        return KDTREE_PREFETCH_EMIT_CONTINUE;
    }

    return MANGLE(kdtree_prefetch_emit_hint)
        (kd,
         LOW_HR(kd, D, nodeid),
         bbox_bytes,
         KDTREE_PREFETCH_ARRAY_BBOX,
         KDTREE_PREFETCH_PRIORITY_METADATA,
         sink);
}

/*
 * Traverse only KD topology and emit the DATA and PERM intervals that the
 * scalar range search will inspect. Query pruning and child order mirror the
 * established solver path. No result is allocated and payload is never
 * dereferenced here. Page alignment, merging, and delivery remain the
 * caller's responsibility.
 */
int MANGLE(kdtree_rangesearch_prefetch_prepare)
     (const kdtree_t *kd,
      const void *vquery,
      double maxd2,
      int options,
      const kdtree_prefetch_sink_t *sink) {
    int nodestack[KDTREE_PREFETCH_STACK_MAX];
    ttype tquery[KDTREE_MAX_DIM];
    int stackpos = 0;
    int D;
    anbool do_wholenode_check;
    anbool do_precheck = FALSE;
    anbool do_l1precheck = FALSE;
    anbool use_bboxes = FALSE;
    anbool use_splits = FALSE;
    anbool use_tquery = FALSE;
    anbool use_tsplit = FALSE;
    anbool use_tmath = FALSE;
    anbool use_bigtmath = FALSE;
    double maxdist;
    double dtl1 = 0.0;
    double dtl2 = 0.0;
    double dtlinf = 0.0;
    ttype tl1 = 0;
    ttype tl2 = 0;
    ttype tlinf = 0;
    bigttype bigtl2 = 0;
    const etype *query = vquery;

    if (!kd || !query || !sink ||
        !sink->emit || !sink->enabled) {
        return KDTREE_PREFETCH_PREPARE_ERROR;
    }
    if (!kd->io ||
        !sink->enabled(sink->userdata, kd->io)) {
        return KDTREE_PREFETCH_PREPARE_NOT_APPLICABLE;
    }

#if defined(KD_DIM)
    assert(kd->ndim == KD_DIM);
    D = KD_DIM;
#else
    D = kd->ndim;
#endif
    if (D <= 0 || D > KDTREE_MAX_DIM) {
        return KDTREE_PREFETCH_PREPARE_ERROR;
    }

    do_wholenode_check = !(options & KD_OPTIONS_SMALL_RADIUS);
    if ((options & KD_OPTIONS_SPLIT_PRECHECK) &&
        kd->bb.any && kd->splitdim) {
        do_precheck = TRUE;
    }
    if ((options & KD_OPTIONS_L1_PRECHECK) && kd->bb.any) {
        do_l1precheck = TRUE;
    }

    if (!kd->split.any) {
        if (!kd->bb.any) {
            return KDTREE_PREFETCH_PREPARE_NOT_APPLICABLE;
        }
        use_bboxes = TRUE;
    } else if (kd->bb.any &&
               !(options & KD_OPTIONS_USE_SPLIT)) {
        use_bboxes = TRUE;
    } else {
        use_splits = TRUE;
    }

    maxdist = sqrt(maxd2);
    if (TTYPE_INTEGER &&
        (kd->split.any || do_precheck || do_l1precheck)) {
        use_tquery = ttype_query(kd, query, tquery);
    }
    if (TTYPE_INTEGER && use_tquery) {
        dtl1 = DIST_ET(kd, maxdist * sqrt(D),);
        dtl2 = DIST2_ET(kd, maxd2, );
        dtlinf = DIST_ET(kd, maxdist, );
        tl1 = ceil(dtl1);
        tlinf = ceil(dtlinf);
        bigtl2 = ceil(dtl2);
        tl2 = bigtl2;
    }
    use_tsplit = use_tquery && (dtlinf < TTYPE_MAX);
    if (do_l1precheck && dtl1 > TTYPE_MAX) {
        do_l1precheck = FALSE;
    }
    if (TTYPE_INTEGER && use_tquery && kd->bb.any) {
        if (dtl2 < TTYPE_MAX) {
            use_tmath = TRUE;
        } else if (dtl2 < BIGTTYPE_MAX) {
            use_bigtmath = TRUE;
        }
        if (use_bigtmath &&
            (options & KD_OPTIONS_NO_BIG_INT_MATH)) {
            use_bigtmath = FALSE;
        }
    }

    nodestack[0] = 0;
    while (stackpos >= 0) {
        int nodeid = nodestack[stackpos--];
        int dim = -1;
        int L;
        int R;
        ttype split = 0;
        int emit_status;

        if (!MANGLE(kdtree_prefetch_valid_node)(kd, nodeid)) {
            return KDTREE_PREFETCH_PREPARE_ERROR;
        }

        if (KD_IS_LEAF(kd, nodeid)) {
            emit_status = MANGLE(kdtree_prefetch_emit_leaf_metadata)
                (kd, nodeid, sink);
            if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
                return emit_status < 0
                    ? KDTREE_PREFETCH_PREPARE_ERROR
                    : KDTREE_PREFETCH_PREPARE_REFUSED;
            }
            L = kdtree_left(kd, nodeid);
            R = kdtree_right(kd, nodeid);
            emit_status = MANGLE(kdtree_prefetch_emit_payload)
                (kd, D, L, R, sink);
            if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
                return emit_status < 0
                    ? KDTREE_PREFETCH_PREPARE_ERROR
                    : KDTREE_PREFETCH_PREPARE_REFUSED;
            }
            continue;
        }

        if (kd->splitdim) {
            dim = kd->splitdim[nodeid];
        }

        if (use_bboxes) {
            ttype *tlo = NULL;
            ttype *thi = NULL;
            anbool wholenode = FALSE;

            emit_status = MANGLE(kdtree_prefetch_emit_bbox_metadata)
                (kd, D, nodeid, sink);
            if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
                return emit_status < 0
                    ? KDTREE_PREFETCH_PREPARE_ERROR
                    : KDTREE_PREFETCH_PREPARE_REFUSED;
            }
            if (!bboxes(kd, nodeid, &tlo, &thi, D) ||
                !tlo || !thi) {
                return KDTREE_PREFETCH_PREPARE_ERROR;
            }

            if (do_precheck && nodeid) {
                anbool isleftchild = KD_IS_LEFT_CHILD(nodeid);
                int pdim;
                anbool cut;

                if (kd->splitdim) {
                    pdim = kd->splitdim[KD_PARENT(nodeid)];
                } else {
                    pdim = kd->split.TTYPE[KD_PARENT(nodeid)];
                    pdim &= kd->dimmask;
                }
                if (TTYPE_INTEGER && use_tquery) {
                    if (isleftchild) {
                        cut = ((tquery[pdim] > thi[pdim]) &&
                               (tquery[pdim] - thi[pdim] > tlinf));
                    } else {
                        cut = ((tlo[pdim] > tquery[pdim]) &&
                               (tlo[pdim] - tquery[pdim] > tlinf));
                    }
                } else {
                    etype bb;

                    if (isleftchild) {
                        bb = POINT_TE(kd, pdim, thi[pdim]);
                        cut = (query[pdim] - bb > maxdist);
                    } else {
                        bb = POINT_TE(kd, pdim, tlo[pdim]);
                        cut = (bb - query[pdim] > maxdist);
                    }
                }
                if (cut) {
                    continue;
                }
            }

            if (TTYPE_INTEGER && do_l1precheck && use_tquery &&
                bb_point_l1mindist_exceeds_ttype(
                    tlo, thi, tquery, D, tl1, tlinf)) {
                continue;
            }

            if (TTYPE_INTEGER && use_tmath) {
                if (bb_point_mindist2_exceeds_ttype(
                        tlo, thi, tquery, D, tl2)) {
                    continue;
                }
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds_ttype(
                        tlo, thi, tquery, D, tl2);
            } else if (TTYPE_INTEGER && use_bigtmath) {
                if (bb_point_mindist2_exceeds_bigttype(
                        tlo, thi, tquery, D, bigtl2)) {
                    continue;
                }
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds_bigttype(
                        tlo, thi, tquery, D, bigtl2);
            } else {
                etype bblo[KDTREE_MAX_DIM];
                etype bbhi[KDTREE_MAX_DIM];
                int d;

                for (d = 0; d < D; d++) {
                    bblo[d] = POINT_TE(kd, d, tlo[d]);
                    bbhi[d] = POINT_TE(kd, d, thi[d]);
                }
                if (bb_point_mindist2_exceeds(
                        bblo, bbhi, query, D, maxd2)) {
                    continue;
                }
                wholenode = do_wholenode_check &&
                    !bb_point_maxdist2_exceeds(
                        bblo, bbhi, query, D, maxd2);
            }

            if (wholenode) {
                L = kdtree_left(kd, nodeid);
                R = kdtree_right(kd, nodeid);
                emit_status = MANGLE(kdtree_prefetch_emit_payload)
                    (kd, D, L, R, sink);
                if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
                    return emit_status < 0
                        ? KDTREE_PREFETCH_PREPARE_ERROR
                        : KDTREE_PREFETCH_PREPARE_REFUSED;
                }
                continue;
            }

            if (MANGLE(kdtree_prefetch_push_node)
                    (kd, nodestack, &stackpos,
                     KD_CHILD_LEFT(nodeid)) ||
                MANGLE(kdtree_prefetch_push_node)
                    (kd, nodestack, &stackpos,
                     KD_CHILD_RIGHT(nodeid))) {
                return KDTREE_PREFETCH_PREPARE_ERROR;
            }
            continue;
        }

        if (!use_splits) {
            return KDTREE_PREFETCH_PREPARE_ERROR;
        }
        emit_status = MANGLE(kdtree_prefetch_emit_split_metadata)
            (kd, nodeid, sink);
        if (emit_status != KDTREE_PREFETCH_EMIT_CONTINUE) {
            return emit_status < 0
                ? KDTREE_PREFETCH_PREPARE_ERROR
                : KDTREE_PREFETCH_PREPARE_REFUSED;
        }
        split = *KD_SPLIT(kd, nodeid);
        if (!kd->splitdim && TTYPE_INTEGER) {
            bigint tmpsplit = split;

            dim = tmpsplit & kd->dimmask;
            split = tmpsplit & kd->splitmask;
        }
        if (dim < 0 || dim >= D) {
            return KDTREE_PREFETCH_PREPARE_ERROR;
        }

        if (TTYPE_INTEGER && use_tsplit) {
            if (tquery[dim] < split) {
                if (MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_LEFT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
                if (split - tquery[dim] <= tlinf &&
                    MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_RIGHT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
            } else {
                if (MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_RIGHT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
                if (tquery[dim] - split <= tlinf &&
                    MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_LEFT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
            }
        } else {
            dtype rsplit = POINT_TE(kd, dim, split);

            if (query[dim] < rsplit) {
                if (MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_LEFT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
                if (rsplit - query[dim] <= maxdist &&
                    MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_RIGHT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
            } else {
                if (MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_RIGHT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
                if (query[dim] - rsplit <= maxdist &&
                    MANGLE(kdtree_prefetch_push_node)
                        (kd, nodestack, &stackpos,
                         KD_CHILD_LEFT(nodeid))) {
                    return KDTREE_PREFETCH_PREPARE_ERROR;
                }
            }
        }
    }

    return KDTREE_PREFETCH_PREPARE_COMPLETE;
}
