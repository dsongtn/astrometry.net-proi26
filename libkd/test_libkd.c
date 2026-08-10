/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <strings.h>

#include "os-features.h"
#include "errors.h"
#include "cutest.h"
#include "kdtree.h"
#include "kdtree_continuation_internal.h"
#include "kdtree_prefetch_internal.h"
#include "mathutil.h"
#include "an-fls.h"

#include "test_libkd_common.c"

typedef struct span_test_sink {
    void* mapping;
    kdtree_rangesearch_span_t spans[4096];
    size_t count;
    size_t capacity;
} span_test_sink_t;

static int span_test_enabled(void* opaque, void* mapping) {
    span_test_sink_t* sink = opaque;

    return sink && sink->mapping == mapping;
}

static int span_test_emit(
    void* opaque,
    const kdtree_prefetch_hint_t* hint) {
    (void)opaque;
    return hint ? KDTREE_PREFETCH_EMIT_CONTINUE
                : KDTREE_PREFETCH_EMIT_ERROR;
}

static int span_test_emit_span(
    void* opaque,
    const kdtree_rangesearch_span_t* span) {
    span_test_sink_t* sink = opaque;

    if (!sink || !span || sink->count >= sink->capacity ||
        sink->capacity >
            sizeof(sink->spans) / sizeof(sink->spans[0])) {
        return KDTREE_PREFETCH_EMIT_REFUSED;
    }
    sink->spans[sink->count++] = *span;
    return KDTREE_PREFETCH_EMIT_CONTINUE;
}

static void span_test_compare_results(
    CuTest* ct,
    const char* case_name,
    const kdtree_qres_t* scalar,
    const kdtree_qres_t* compiled,
    int dimension,
    int options) {
    size_t result_count;

    CuAssertPtrNotNull(ct, scalar);
    CuAssertPtrNotNull(ct, compiled);
    CuAssertIntEquals(ct, scalar->nres, compiled->nres);
    result_count = (size_t)scalar->nres;
    if (result_count) {
        CuAssert(ct, case_name,
                 scalar->inds && compiled->inds &&
                 !memcmp(scalar->inds, compiled->inds,
                         result_count * sizeof(*scalar->inds)));
    }
    if ((options & KD_OPTIONS_COMPUTE_DISTS) ||
        (options & KD_OPTIONS_SORT_DISTS)) {
        CuAssert(ct, "compiled span distances must be identical",
                 (!result_count ||
                  (scalar->sdists && compiled->sdists &&
                   !memcmp(scalar->sdists, compiled->sdists,
                           result_count * sizeof(*scalar->sdists)))));
    }
    if (options & KD_OPTIONS_RETURN_POINTS) {
        CuAssert(ct, "compiled span points must be identical",
                 (!result_count ||
                  (scalar->results.d && compiled->results.d &&
                   !memcmp(scalar->results.d, compiled->results.d,
                           result_count * (size_t)dimension *
                               sizeof(*scalar->results.d)))));
    }
}

static void span_test_run_case(
    CuTest* ct,
    const char* case_name,
    int treetype,
    int build_options,
    int search_options,
    const double* query,
    double maxd2,
    anbool remove_perm,
    anbool expect_accept_all,
    anbool exercise_reuse) {
    const int N = 257;
    const int D = 4;
    double second_query[4] = {0.19, 0.73, 0.41, 0.27};
    double* data;
    kdtree_t* kd;
    kdtree_qres_t* scalar = NULL;
    kdtree_qres_t* compiled = NULL;
    kdtree_prefetch_sink_t prefetch;
    span_test_sink_t spans;
    size_t span_index;
    anbool saw_accept_all = FALSE;
    int i;
    int d;

    data = malloc((size_t)N * (size_t)D * sizeof(double));
    CuAssertPtrNotNull(ct, data);
    for (i = 0; i < N; i++) {
        for (d = 0; d < D; d++) {
            unsigned int value =
                (unsigned int)(i * 37 + d * 53 + i * d * 11);

            data[(size_t)i * (size_t)D + (size_t)d] =
                (double)(value % 997U) / 996.0;
        }
    }
    kd = kdtree_build(
        NULL, data, N, D, 8, treetype, build_options);
    CuAssertPtrNotNull(ct, kd);
    if (remove_perm) {
        free(kd->perm);
        kd->perm = NULL;
    }

    memset(&spans, 0, sizeof(spans));
    memset(&prefetch, 0, sizeof(prefetch));
    spans.mapping = &spans;
    spans.capacity = sizeof(spans.spans) / sizeof(spans.spans[0]);
    kd->io = spans.mapping;
    prefetch.userdata = &spans;
    prefetch.enabled = span_test_enabled;
    prefetch.emit = span_test_emit;
    prefetch.emit_span = span_test_emit_span;
    CuAssertIntEquals(
        ct,
        KDTREE_PREFETCH_PREPARE_COMPLETE,
        kdtree_rangesearch_prefetch_prepare(
            kd, query, maxd2, search_options, &prefetch));
    for (span_index = 0U; span_index < spans.count; span_index++) {
        saw_accept_all |= spans.spans[span_index].mode ==
            KDTREE_RANGESEARCH_SPAN_ACCEPT_ALL;
    }
    CuAssertIntEquals(ct, expect_accept_all, saw_accept_all);

    scalar = kdtree_rangesearch_options_reuse(
        kd, scalar, query, maxd2, search_options);
    compiled = kdtree_rangesearch_execute_spans(
        kd, compiled, query, maxd2, search_options,
        spans.spans, spans.count);
    span_test_compare_results(
        ct, case_name, scalar, compiled, D, search_options);

    if (exercise_reuse) {
        spans.count = 0U;
        CuAssertIntEquals(
            ct,
            KDTREE_PREFETCH_PREPARE_COMPLETE,
            kdtree_rangesearch_prefetch_prepare(
                kd, second_query, 0.12,
                search_options, &prefetch));
        scalar = kdtree_rangesearch_options_reuse(
            kd, scalar, second_query, 0.12, search_options);
        compiled = kdtree_rangesearch_execute_spans(
            kd, compiled, second_query, 0.12, search_options,
            spans.spans, spans.count);
        span_test_compare_results(
            ct, case_name, scalar, compiled, D, search_options);
    }

    kd->io = NULL;
    kdtree_free_query(compiled);
    kdtree_free_query(scalar);
    kdtree_free(kd);
    free(data);
}

static int calculate_R(int leafid, int nlevels, int N) {
    int l;
    unsigned int mask, L;

    int nbottom = 1 << (nlevels - 1);
    

    mask = (1 << (nlevels-1));
    L = 0;
    // Compute the L index of the node one to the right of this node.
    int nextguy = leafid + 1;
    if (nextguy == nbottom)
        return N-1;
    for (l=0; l<(nlevels-1); l++) {
        mask /= 2;
        if (nextguy & mask) {
            L += N/2;
            N = (N+1)/2;
        } else {
            N = N/2;
        }
    }
    L--;
    return L;
}

int linearR(int leafid, int nbottom, int N) {
    int64_t res = leafid + 1;
    res *= N;
    res /= nbottom;
    return res - 1;
}

double linearRF(int leafid, int nbottom, int N) {
    double res = leafid + 1.0;
    res *= N;
    res /= (double)nbottom;
    return res - 1.0;
}

void tst_1(CuTest* ct) {
    kdtree_t* kd;
    double * data;
    int N = 1000;
    int Nleaf = 5;
    int D = 3;
    int i;

    data = random_points_d(N, D);
    kd = build_tree(ct, data, N, D, Nleaf, KDTT_DOUBLE, KD_BUILD_SPLIT);

    printf("kd->nbottom = %i, kd->nlevels = %i.\n", kd->nbottom, kd->nlevels);

    for (i=0; i<kd->nbottom; i++) {
        int R1 = kdtree_right(kd, kd->ninterior + i);
        int R2 = calculate_R(i, kd->nlevels, N);
        int R3 = linearR(i, kd->nbottom, N);
        double d3 = linearRF(i, kd->nbottom, N);

        printf("%i %i %i %g\n", R1, R2, R3, d3);
        printf("                               %s   %g\n",
               (R1 != R3) ? "***" : "   ",
               (double)R3 - d3);
        /*
         CuAssertIntEquals(ct, R1, R2);
         CuAssertIntEquals(ct, R1, R3);
         */
    }
    kdtree_free(kd);
}

static void compute_splitbits(int ndim, uint32_t* dimmask, uint32_t* dimbits, uint32_t* splitmask) {
    int D;
    int bits;
    uint32_t val;
    D = ndim;
    bits = 0;
    val = 1;
    while (val < D) {
        bits++;
        val *= 2;
    }
    *dimmask = val - 1;
    *dimbits = bits;
    *splitmask = ~(*dimmask);
}

void test_splitbits(CuTest* ct) {
    uint32_t dmask, dbits, smask;
    int dim;

    compute_splitbits(1, &dmask, &dbits, &smask);
    CuAssertIntEquals(ct, 0x0, dbits);
    CuAssertIntEquals(ct, 0x0, dmask);
    CuAssertIntEquals(ct,~0x0, smask);

    compute_splitbits(2, &dmask, &dbits, &smask);
    CuAssertIntEquals(ct, 0x1, dbits);
    CuAssertIntEquals(ct, 0x1, dmask);
    CuAssertIntEquals(ct,~0x1, smask);

    for (dim=3; dim<=4; dim++) {
        compute_splitbits(dim, &dmask, &dbits, &smask);
        CuAssertIntEquals(ct, 0x2, dbits);
        CuAssertIntEquals(ct, 0x3, dmask);
        CuAssertIntEquals(ct,~0x3, smask);
    }

    for (dim=5; dim<=8; dim++) {
        compute_splitbits(dim, &dmask, &dbits, &smask);
        CuAssertIntEquals(ct, 0x3, dbits);
        CuAssertIntEquals(ct, 0x7, dmask);
        CuAssertIntEquals(ct,~0x7, smask);
    }

    for (dim=9; dim<=16; dim++) {
        compute_splitbits(dim, &dmask, &dbits, &smask);
        CuAssertIntEquals(ct, 0x4, dbits);
        CuAssertIntEquals(ct, 0xF, dmask);
        CuAssertIntEquals(ct,~0xF, smask);
    }
}

void test_short_partition(CuTest* ct) {
    kdtree_t* kd;
    double* data;
    int N = 21;
    int Nleaf = 16;
    int D = 2;
    int i;
    double minval[D], maxval[D];
    uint16_t cdata[] = { 12669, 12669, 12669, 12669, 12669, 12669, 
                         12669, 12669, 12669, 12669, 12669, 12669,
                         13860, 13913, 14164, 14557, 15283, 17130,
                         17130, 17130, 17130 };


    minval[0] = -0.20710678118654757;
    minval[1] = -0.20710678118654757;
    maxval[0] =  1.2071067811865475;
    maxval[1] =  1.2071067811865475;

    data = calloc(N*D, sizeof(double));
    // convert from "cdata" to "data" space (I got the test data above
    // from the internal data representation of a problem tree).
    double scale = (maxval[0] - minval[0]) / (double)UINT16_MAX;
    for (i=0; i<N; i++) {
        data[2*i+0] = minval[0];
        data[2*i+1] = minval[0] + (double)cdata[i] * scale;
    }
    kd = kdtree_build_2(NULL, data, N, D, Nleaf, KDTT_DSS, KD_BUILD_SPLIT,
                        minval, maxval);
    CuAssertPtrNotNull(ct, kd);
    free(data);

    /*
     printf("kd:\n");
     kdtree_print(kd);
     */

    uint16_t* kdata = kd->data.s;
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, kdata[2*i+0], 0);
        CuAssertIntEquals(ct, kdata[2*i+1], cdata[i]);
    }
    CuAssertIntEquals(ct, 0, kdtree_check(kd));
    kdtree_free(kd);
}

void test_empty_node(CuTest* ct) {
    kdtree_t* kd;
    double* data;
    int N = 21;
    int Nleaf = 16;
    int D = 2;
    int i, ok;
    double minval[D], maxval[D];

    minval[0] = 0;
    minval[1] = 0;
    maxval[0] = 1;
    maxval[1] = 1;

    data = calloc(N*D, sizeof(double));
    // convert from "cdata" to "data" space
    double scale = (maxval[0] - minval[0]) / (double)UINT16_MAX;
    for (i=0; i<N; i++) {
        data[2*i+0] = minval[0];
        data[2*i+1] = minval[0] + 3 * scale;
    }

    kd = kdtree_build_2(NULL, data, N, D, Nleaf, KDTT_DSS, KD_BUILD_SPLIT,
                        minval, maxval);
    CuAssertPtrNotNull(ct, kd);
    free(data);

    uint16_t* kdata = kd->data.s;
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, 0, kdata[2*i+0]);
        CuAssertIntEquals(ct, 3, kdata[2*i+1]);
    }

    ok = kdtree_check(kd);
    CuAssertIntEquals(ct, 0, ok);
    kdtree_free(kd);
}

static inline u8 node_level(int nodeid) {
    int val = (nodeid + 1) >> 1;
    u8 level = 0;
    while (val) {
        val = val >> 1;
        level++;
    }
    return level;
}

void test_2(CuTest* ct) {
    int N = 1024;
    int i;

    for (i=0; i<N; i++) {
        int L1 = node_level(i);
        int L3 = an_flsB(i+1);
        CuAssertIntEquals(ct, L1, L3);
    }
}

void test_nlevels(CuTest* ct) {
    // No nodes: no levels!
    CuAssertIntEquals(ct, 0, kdtree_nnodes_to_nlevels(0));
    // Single node.
    CuAssertIntEquals(ct, 1, kdtree_nnodes_to_nlevels(1));
    // sort of invalid input - incomplete level...
    CuAssertIntEquals(ct, 1, kdtree_nnodes_to_nlevels(2));
    CuAssertIntEquals(ct, 2, kdtree_nnodes_to_nlevels(3));
    CuAssertIntEquals(ct, 10, kdtree_nnodes_to_nlevels(1023));
}

static void run_continuation_to_completion(
    CuTest* ct,
    const kdtree_t* kd,
    const double* query,
    double maxd2,
    int options,
    size_t node_budget,
    kdtree_qres_t* reuse,
    kdtree_qres_t** result_out) {
    kdtree_rangesearch_continuation_t continuation;
    kdtree_rangesearch_continuation_init_status_t init_status;
    kdtree_rangesearch_continuation_status_t status;

    CuAssertPtrNotNull(ct, result_out);
    *result_out = NULL;

    kdtree_rangesearch_continuation_zero(&continuation);

    init_status = kdtree_rangesearch_continuation_init(
        &continuation,
        kd,
        reuse,
        query,
        maxd2,
        options);

    CuAssertIntEquals(
        ct,
        KDTREE_RANGESEARCH_CONTINUATION_INIT_OK,
        init_status);

    do {
        status = kdtree_rangesearch_continuation_step(
            &continuation,
            node_budget);
    } while (status == KDTREE_RANGESEARCH_CONTINUATION_MORE);

    CuAssertIntEquals(
        ct,
        KDTREE_RANGESEARCH_CONTINUATION_DONE,
        status);

    *result_out = kdtree_rangesearch_continuation_finish(
        &continuation);

    CuAssertPtrNotNull(ct, *result_out);
    CuAssert(
        ct,
        "continuation must visit at least the root",
        continuation.nodes_visited > 0);

    kdtree_rangesearch_continuation_cleanup(&continuation);
}

static void assert_query_indices_and_distances_equal(
    CuTest* ct,
    const kdtree_qres_t* expected,
    const kdtree_qres_t* actual) {
    unsigned int i;

    CuAssertPtrNotNull(ct, expected);
    CuAssertPtrNotNull(ct, actual);
    CuAssertIntEquals(ct, (int)expected->nres, (int)actual->nres);

    for (i = 0; i < expected->nres; i++) {
        CuAssertIntEquals(
            ct,
            (int)expected->inds[i],
            (int)actual->inds[i]);

        CuAssert(
            ct,
            "range-search squared distance differs",
            expected->sdists[i] == actual->sdists[i]);
    }
}

static void assert_query_results_equal(
    CuTest* ct,
    const kdtree_qres_t* expected,
    const kdtree_qres_t* actual,
    int dimension) {
    unsigned int i;
    int d;

    CuAssertPtrNotNull(ct, expected);
    CuAssertPtrNotNull(ct, actual);
    CuAssertIntEquals(ct, (int)expected->nres, (int)actual->nres);

    for (i = 0; i < expected->nres; i++) {
        CuAssertIntEquals(
            ct,
            (int)expected->inds[i],
            (int)actual->inds[i]);

        CuAssert(
            ct,
            "continuation squared distance differs",
            expected->sdists[i] == actual->sdists[i]);

        for (d = 0; d < dimension; d++) {
            CuAssert(
                ct,
                "continuation result point differs",
                expected->results.d[(size_t)i * (size_t)dimension +
                                    (size_t)d] ==
                actual->results.d[(size_t)i * (size_t)dimension +
                                  (size_t)d]);
        }
    }
}

void test_rangesearch_index_distance_only(CuTest* ct) {
    const int N = 257;
    const int D = 4;
    const int Nleaf = 8;
    const int base_options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_SORT_DISTS |
        KD_OPTIONS_USE_SPLIT;
    double query[4] = {0.47, 0.51, 0.39, 0.62};
    double* data;
    kdtree_t* kd;
    kdtree_qres_t* with_points;
    kdtree_qres_t* without_points;
    kdtree_qres_t* reuse;
    kdtree_qres_t* continuation;
    kdtree_qres_t* defaults;
    int i;
    int d;

    data = malloc((size_t)N * (size_t)D * sizeof(double));
    CuAssertPtrNotNull(ct, data);

    for (i = 0; i < N; i++) {
        for (d = 0; d < D; d++) {
            unsigned int value =
                (unsigned int)(i * 37 + d * 53 + i * d * 11);

            data[(size_t)i * (size_t)D + (size_t)d] =
                (double)(value % 997U) / 996.0;
        }
    }

    kd = kdtree_build(
        NULL,
        data,
        N,
        D,
        Nleaf,
        KDTT_DSS,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM);

    CuAssertPtrNotNull(ct, kd);
    free(data);

    with_points = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        query,
        0.08,
        base_options | KD_OPTIONS_RETURN_POINTS);

    CuAssertPtrNotNull(ct, with_points);
    CuAssertPtrNotNull(ct, with_points->results.any);
    CuAssert(ct, "expected non-empty query result", with_points->nres > 0);

    without_points = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        query,
        0.08,
        base_options);

    CuAssertPtrNotNull(ct, without_points);
    CuAssert(
        ct,
        "point storage must remain absent when not requested",
        without_points->results.any == NULL);
    CuAssertPtrNotNull(ct, without_points->inds);
    CuAssertPtrNotNull(ct, without_points->sdists);
    assert_query_indices_and_distances_equal(
        ct,
        with_points,
        without_points);

    reuse = kdtree_rangesearch_options_reuse(
        kd,
        with_points,
        query,
        0.08,
        base_options);

    CuAssert(
        ct,
        "reuse query must preserve the result container",
        reuse == with_points);
    CuAssert(
        ct,
        "reuse query must release stale point storage",
        reuse->results.any == NULL);
    assert_query_indices_and_distances_equal(
        ct,
        without_points,
        reuse);

    continuation = kdtree_rangesearch_continuation_execute(
        kd,
        NULL,
        query,
        0.08,
        base_options,
        KDTREE_RANGESEARCH_CONTINUATION_RUN_TO_COMPLETION);

    CuAssertPtrNotNull(ct, continuation);
    CuAssert(
        ct,
        "continuation must honor index-and-distance-only output",
        continuation->results.any == NULL);
    assert_query_indices_and_distances_equal(
        ct,
        without_points,
        continuation);

    defaults = kdtree_rangesearch(kd, query, 0.08);
    CuAssertPtrNotNull(ct, defaults);
    CuAssertPtrNotNull(ct, defaults->results.any);
    assert_query_indices_and_distances_equal(
        ct,
        without_points,
        defaults);

    kdtree_free_query(defaults);
    kdtree_free_query(continuation);
    kdtree_free_query(reuse);
    kdtree_free_query(without_points);
    kdtree_free(kd);
}

void test_rangesearch_compiled_span_parity(CuTest* ct) {
    double query[4] = {0.47, 0.51, 0.39, 0.62};
    double boundary_query[4] = {
        0.0, 53.0 / 996.0, 106.0 / 996.0, 159.0 / 996.0};
    double empty_query[4] = {2.0, 2.0, 2.0, 2.0};
    double* data;
    kdtree_t* kd;
    kdtree_qres_t* preserved;
    kdtree_prefetch_sink_t prefetch;
    span_test_sink_t spans;

    span_test_run_case(
        ct, "compiled split double indices must be identical",
        KDTT_DOUBLE,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_NO_RESIZE_RESULTS |
            KD_OPTIONS_USE_SPLIT,
        query, 0.08, FALSE, FALSE, TRUE);
    span_test_run_case(
        ct, "compiled no-perm indices must be identical",
        KDTT_DOUBLE,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS | KD_OPTIONS_USE_SPLIT,
        query, 0.08, TRUE, FALSE, FALSE);
    span_test_run_case(
        ct, "compiled bbox accept-all indices must be identical",
        KDTT_DOUBLE, KD_BUILD_BBOX,
        KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_RETURN_POINTS |
            KD_OPTIONS_SORT_DISTS,
        query, 4.0, FALSE, TRUE, TRUE);
    span_test_run_case(
        ct, "compiled integer bbox indices must be identical",
        KDTT_DSS, KD_BUILD_BBOX,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_SORT_DISTS |
            KD_OPTIONS_SPLIT_PRECHECK |
            KD_OPTIONS_L1_PRECHECK,
        query, 0.08, FALSE, FALSE, FALSE);
    span_test_run_case(
        ct, "compiled unsigned split indices must be identical",
        KDTT_DUU,
        KD_BUILD_BBOX | KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_USE_SPLIT,
        query, 0.08, FALSE, FALSE, FALSE);
    span_test_run_case(
        ct, "compiled mixed bbox split indices must be identical",
        KDTT_DOUBLE,
        KD_BUILD_BBOX | KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_SPLIT_PRECHECK,
        query, 0.08, FALSE, FALSE, FALSE);
    span_test_run_case(
        ct, "compiled boundary indices must be identical",
        KDTT_DOUBLE,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_USE_SPLIT,
        boundary_query, 0.0, FALSE, FALSE, FALSE);
    span_test_run_case(
        ct, "compiled empty indices must be identical",
        KDTT_DOUBLE,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_USE_SPLIT,
        empty_query, 0.0, FALSE, FALSE, FALSE);

    data = calloc(32U * 4U, sizeof(*data));
    CuAssertPtrNotNull(ct, data);
    kd = kdtree_build(
        NULL, data, 32, 4, 4, KDTT_DOUBLE,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM);
    CuAssertPtrNotNull(ct, kd);
    memset(&spans, 0, sizeof(spans));
    memset(&prefetch, 0, sizeof(prefetch));
    spans.mapping = &spans;
    spans.capacity = 0U;
    kd->io = spans.mapping;
    prefetch.userdata = &spans;
    prefetch.enabled = span_test_enabled;
    prefetch.emit = span_test_emit;
    prefetch.emit_span = span_test_emit_span;
    CuAssertIntEquals(
        ct,
        KDTREE_PREFETCH_PREPARE_REFUSED,
        kdtree_rangesearch_prefetch_prepare(
            kd, query, 4.0,
            KD_OPTIONS_USE_SPLIT,
            &prefetch));
    CuAssertIntEquals(ct, 0, (int)spans.count);

    preserved = kdtree_rangesearch_options_reuse(
        kd, NULL, query, 0.08,
        KD_OPTIONS_SMALL_RADIUS |
            KD_OPTIONS_COMPUTE_DISTS |
            KD_OPTIONS_USE_SPLIT);
    CuAssertPtrNotNull(ct, preserved);
    spans.spans[0].left = -1;
    spans.spans[0].right = 0;
    spans.spans[0].mode = KDTREE_RANGESEARCH_SPAN_FILTER;
    {
        unsigned int preserved_count = preserved->nres;

        CuAssertPtrEquals(
            ct, NULL,
            kdtree_rangesearch_execute_spans(
                kd, preserved, query, 0.08,
                KD_OPTIONS_SMALL_RADIUS |
                    KD_OPTIONS_COMPUTE_DISTS |
                    KD_OPTIONS_USE_SPLIT,
                spans.spans, 1U));
        CuAssertIntEquals(ct, preserved_count, preserved->nres);
    }
    kd->io = NULL;
    kdtree_free_query(preserved);
    kdtree_free(kd);
    free(data);
}

void test_rangesearch_continuation_split_parity(CuTest* ct) {
    const int N = 257;
    const int D = 4;
    const int Nleaf = 8;
    const int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_RETURN_POINTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;

    double query[4] = {0.47, 0.51, 0.39, 0.62};
    double* data;
    kdtree_t* kd;
    kdtree_qres_t* one_shot;
    kdtree_qres_t* fast_path;
    kdtree_qres_t* stepped_one;
    kdtree_qres_t* stepped_seven;
    kdtree_qres_t* stepped_reuse;
    kdtree_qres_t* reuse;
    double reuse_seed_query[4] = {0.12, 0.18, 0.24, 0.30};
    int i;
    int d;

    data = malloc((size_t)N * (size_t)D * sizeof(double));
    CuAssertPtrNotNull(ct, data);

    for (i = 0; i < N; i++) {
        for (d = 0; d < D; d++) {
            unsigned int value =
                (unsigned int)(i * 37 + d * 53 + i * d * 11);

            data[(size_t)i * (size_t)D + (size_t)d] =
                (double)(value % 997U) / 996.0;
        }
    }

    kd = kdtree_build(
        NULL,
        data,
        N,
        D,
        Nleaf,
        KDTT_DSS,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM);

    CuAssertPtrNotNull(ct, kd);
    free(data);

    one_shot = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        query,
        0.08,
        options);

    CuAssertPtrNotNull(ct, one_shot);

    fast_path = kdtree_rangesearch_continuation_execute(
        kd,
        NULL,
        query,
        0.08,
        options,
        KDTREE_RANGESEARCH_CONTINUATION_RUN_TO_COMPLETION);

    CuAssertPtrNotNull(ct, fast_path);

    run_continuation_to_completion(
        ct,
        kd,
        query,
        0.08,
        options,
        1,
        NULL,
        &stepped_one);

    run_continuation_to_completion(
        ct,
        kd,
        query,
        0.08,
        options,
        7,
        NULL,
        &stepped_seven);

    reuse = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        reuse_seed_query,
        0.01,
        options);

    CuAssertPtrNotNull(ct, reuse);

    run_continuation_to_completion(
        ct,
        kd,
        query,
        0.08,
        options,
        3,
        reuse,
        &stepped_reuse);

    CuAssert(
        ct,
        "continuation must preserve the caller's reuse object",
        stepped_reuse == reuse);

    assert_query_results_equal(ct, one_shot, fast_path, D);
    assert_query_results_equal(ct, one_shot, stepped_one, D);
    assert_query_results_equal(ct, one_shot, stepped_seven, D);
    assert_query_results_equal(ct, one_shot, stepped_reuse, D);

    kdtree_free_query(stepped_reuse);
    kdtree_free_query(fast_path);
    kdtree_free_query(stepped_seven);
    kdtree_free_query(stepped_one);
    kdtree_free_query(one_shot);
    kdtree_free(kd);
}

static void run_test_nn(CuTest* tc, int treetype, int treeopts,
                        double eps) {
    int N = 1000;
    int Nleaf = 10;
    int D = 3;
    int Q = 10;
    kdtree_t* kd;
    double* origdata;
    double* treedata;
    double query[D];
    int i, q, d;

    srand(0);

    origdata = random_points_d(N, D);
    treedata = malloc(N * D * sizeof(double));
    memcpy(treedata, origdata, N*D*sizeof(double));

    kd = build_tree(tc, treedata, N, D, Nleaf, treetype, treeopts);

    CuAssert(tc, "kd", kd != NULL);

    if (treeopts & KD_BUILD_NO_LR)
        CuAssert(tc, "no lr", kd->lr == NULL);

    for (q=0; q<Q; q++) {
        int ind;
        double d2;
        double trued2;
        int trueind;
        for (d=0; d<D; d++)
            query[d] = rand() / (double)RAND_MAX;

        ind = kdtree_nearest_neighbour(kd, query, &d2);

        trued2 = LARGE_VAL;
        trueind = -1;
        for (i=0; i<N; i++) {
            double d2 = distsq(query, origdata + i*D, D);
            if (d2 < trued2) {
                trueind = i;
                trued2 = d2;
            }
        }

        /*
         printf("Naive : ind %i, dist %g.\n", trueind, sqrt(trued2));
         printf("Kdtree: ind %i, dist %g.\n", kd->perm[ind], sqrt(d2));
         */
        CuAssertIntEquals(tc, kd->perm[ind], trueind);

        if (fabs(sqrt(d2) - sqrt(trued2)) >= eps) {
            printf("Naive : %.12g\n", sqrt(trued2));
            printf("Kdtree: %.12g\n", sqrt(d2));
        }
        
        CuAssertDblEquals(tc, sqrt(d2), sqrt(trued2), eps);
    }

    kdtree_free(kd);
    free(treedata);
    free(origdata);
}

static void run_test_rs_ND(CuTest* tc, int treetype, int treeopts,
                           double eps, int N, int D) {
    int Nleaf = 10;
    int Q = 10;
    double rad2 = 0.01;
    double* origdata;
    double* treedata;
    kdtree_t* kd;
    double query[D];
    int i, q, d;

    srand(0);

    origdata = random_points_d(N, D);
    treedata = malloc(N * D * sizeof(double));
    memcpy(treedata, origdata, N*D*sizeof(double));

    kd = build_tree(tc, treedata, N, D, Nleaf, treetype, treeopts);
    CuAssert(tc, "kd", kd != NULL);

    for (q=0; q<Q; q++) {
        int ind;
        double d2;
        double trued2;
        int ntrue;
        kdtree_qres_t* res;

        for (d=0; d<D; d++)
            query[d] = rand() / (double)RAND_MAX;

        res = kdtree_rangesearch(kd, query, rad2);

        ntrue = 0;
        for (i=0; i<N; i++) {
            double d2 = distsq(query, origdata + i*D, D);
            if (d2 <= rad2) {
                ntrue++;
            }
        }

        CuAssertIntEquals(tc, res->nres, ntrue);

        for (i=0; i<res->nres; i++) {
            ind = res->inds[i];
            d2 = res->sdists[i];
            trued2 = distsq(query, origdata + ind*D, D);
        }

        CuAssert(tc, "res", res != NULL);
        for (i=0; i<res->nres; i++) {
            ind = res->inds[i];
            d2 = res->sdists[i];
            trued2 = distsq(query, origdata + ind*D, D);

            CuAssert(tc, "ind pos", ind >= 0);
            CuAssert(tc, "ind pos", ind < N);
            CuAssert(tc, "inrange", d2 <= rad2);
            CuAssert(tc, "inrange", trued2 <= rad2);
            CuAssert(tc, "d2pos", d2 >= 0.0);
            CuAssert(tc, "trued2pos", trued2 >= 0.0);
            CuAssertDblEquals(tc, sqrt(d2), sqrt(trued2), sqrt(eps));
        }
        /*
         printf("Naive : ind %i, dist %g.\n", trueind, sqrt(trued2));
         printf("Kdtree: ind %i, dist %g.\n", kd->perm[ind], sqrt(d2));
         */

        kdtree_free_query(res);
    }

    kdtree_free(kd);
    free(treedata);
    free(origdata);
}

static void run_test_rs(CuTest* tc, int treetype, int treeopts,
                        double eps) {
    int N = 1000;
    int D = 3;
    run_test_rs_ND(tc, treetype, treeopts, eps, N, D);
}

void test_rs_bb_duu(CuTest* tc) {
    run_test_rs(tc, KDTT_DUU, KD_BUILD_BBOX, 1e-9);
}

void test_rs_bb_ddd_small(CuTest* tc) {
    run_test_rs_ND(tc, KDTT_DOUBLE, KD_BUILD_BBOX, 1e-9, 10, 1);
}

void test_rs_bb_ddd(CuTest* tc) {
    run_test_rs(tc, KDTT_DOUBLE, KD_BUILD_BBOX, 1e-9);
}
void test_rs_split_ddd(CuTest* tc) {
    run_test_rs(tc, KDTT_DOUBLE, KD_BUILD_SPLIT, 1e-9);
}
void test_rs_both_ddd(CuTest* tc) {
    run_test_rs(tc, KDTT_DOUBLE, KD_BUILD_BBOX | KD_BUILD_SPLIT, 1e-9);
}

void test_rs_split_duu(CuTest* tc) {
    run_test_rs(tc, KDTT_DUU, KD_BUILD_SPLIT, 1e-9);
}

/**
 Sadly, does not work.

 void test_rs_split_ddu(CuTest* tc) {
 run_test_rs(tc, KDTT_DDU, KD_BUILD_SPLIT, 1e-9);
 }
 */

void test_rs_bb_dss(CuTest* tc) {
    run_test_rs(tc, KDTT_DSS, KD_BUILD_BBOX, 1e-5);
}
void test_rs_split_dss(CuTest* tc) {
    run_test_rs(tc, KDTT_DSS, KD_BUILD_SPLIT, 1e-5);
}




void test_nn_bb_ddd(CuTest* tc) {
    run_test_nn(tc, KDTT_DOUBLE, KD_BUILD_BBOX, 1e-9);
}

void test_nn_split_ddd(CuTest* tc) {
    run_test_nn(tc, KDTT_DOUBLE, KD_BUILD_SPLIT, 1e-9);
}

void test_nn_both_ddd(CuTest* tc) {
    run_test_nn(tc, KDTT_DOUBLE, KD_BUILD_SPLIT | KD_BUILD_BBOX, 1e-9);
}

void test_nn_split_duu(CuTest* tc) {
    run_test_nn(tc, KDTT_DUU, KD_BUILD_SPLIT, 1e-9);
}

void test_nn_bb_duu(CuTest* tc) {
    run_test_nn(tc, KDTT_DUU, KD_BUILD_BBOX, 1e-9);
}

void test_nn_split_dss(CuTest* tc) {
    run_test_nn(tc, KDTT_DSS, KD_BUILD_SPLIT, 1e-5);
}

void test_nn_split_dssB(CuTest* tc) {
    run_test_nn(tc, KDTT_DSS, KD_BUILD_SPLIT | KD_BUILD_NO_LR | KD_BUILD_SPLITDIM, 1e-5);
}

void test_nn_bb_dss(CuTest* tc) {
    run_test_nn(tc, KDTT_DSS, KD_BUILD_BBOX, 1e-5);
}

void test_nn_bb_dssB(CuTest* tc) {
    run_test_nn(tc, KDTT_DSS, KD_BUILD_BBOX | KD_BUILD_NO_LR, 1e-5);
}



void test_nn_split_ddd_linearlr(CuTest* tc) {
    run_test_nn(tc, KDTT_DOUBLE, KD_BUILD_SPLIT | KD_BUILD_NO_LR | KD_BUILD_LINEAR_LR, 1e-9);
}
void test_nn_split_duu_linearlr(CuTest* tc) {
    run_test_nn(tc, KDTT_DUU, KD_BUILD_SPLIT | KD_BUILD_SPLITDIM | KD_BUILD_NO_LR | KD_BUILD_LINEAR_LR, 1e-9);
}
void test_nn_split_dss_linearlr(CuTest* tc) {
    run_test_nn(tc, KDTT_DSS, KD_BUILD_SPLIT | KD_BUILD_SPLITDIM | KD_BUILD_NO_LR | KD_BUILD_LINEAR_LR, 1e-5);
}

void run_test_lr(CuTest* tc, int D, int Nleaf, int treetype, int treeopts) {
    int i;
    kdtree_t* kd;
    double* treedata;
    int32_t* lr;
    int N;
    for (N=100; N<=1000; N+=9) {
        treedata = random_points_d(N, D);
        kd = build_tree(tc, treedata, N, D, Nleaf, treetype, treeopts);
        CuAssert(tc, "kd", kd != NULL);

        lr = kd->lr;
        kd->lr = NULL;

        for (i=0; i<kd->nbottom; i++) {
            if (i)
                CuAssertIntEquals(tc, lr[i-1]+1, kdtree_left(kd, i + kd->ninterior));
            CuAssertIntEquals(tc, lr[i], kdtree_right(kd, i + kd->ninterior));
        }

        kd->lr = lr;
        kdtree_free(kd);
        free(treedata);
    }
}

void test_lr_ddd(CuTest* tc) {
    run_test_lr(tc, 3, 10, KDTT_DOUBLE, KD_BUILD_SPLIT);
}

void test_no_lr_with_ints(CuTest* tc) {
    double* data;
    kdtree_t* kd;
    int N = 1000;
    int D = 3;
    int Nleaf = 10;
    data = random_points_d(N, D);
    kd = build_tree(tc, data, N, D, Nleaf, KDTT_DSS, KD_BUILD_SPLIT | KD_BUILD_NO_LR);
    CuAssert(tc, "no kd", kd == NULL);
    free(data);

    errors_free();
}
