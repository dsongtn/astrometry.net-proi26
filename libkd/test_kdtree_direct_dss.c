/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cutest.h"
#include "kdtree.h"
#include "kdtree_direct_internal.h"

typedef struct direct_reader {
    const kdtree_t* kd;
    u16* data;
    u32* perm;
    size_t capacity_points;
    int calls;
    int fail_on_call;
    int largest_count;
} direct_reader_t;

static uint32_t direct_test_prng(uint32_t* state) {
    uint32_t value = *state;

    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    *state = value;
    return value;
}

static double* direct_test_points(int npoints, int ndim) {
    double* data;
    uint32_t state = UINT32_C(0x6d2b79f5);
    int i;

    data = malloc((size_t)npoints * (size_t)ndim * sizeof(double));
    if (!data) {
        return NULL;
    }
    for (i = 0; i < npoints * ndim; i++) {
        uint32_t value = direct_test_prng(&state);

        /*
         * Keep points away from an exactly uniform grid while retaining
         * deterministic values at both ends of the quantized U16 range.
         */
        data[i] = (double)(value & UINT32_C(0xffff)) / 65535.0;
    }
    data[0] = 0.0;
    data[1] = 1.0;
    data[2] = 0.0;
    data[3] = 1.0;
    return data;
}

static int direct_reader_init(
    direct_reader_t* reader,
    const kdtree_t* kd,
    size_t capacity_points) {
    memset(reader, 0, sizeof(*reader));
    reader->kd = kd;
    reader->capacity_points = capacity_points;
    reader->data = malloc(capacity_points * (size_t)kd->ndim * sizeof(u16));
    reader->perm = malloc(capacity_points * sizeof(u32));
    if (!reader->data || !reader->perm) {
        free(reader->perm);
        free(reader->data);
        memset(reader, 0, sizeof(*reader));
        errno = ENOMEM;
        return -1;
    }
    return 0;
}

static void direct_reader_cleanup(direct_reader_t* reader) {
    free(reader->perm);
    free(reader->data);
    memset(reader, 0, sizeof(*reader));
}

static int direct_reader_read(
    void* opaque,
    int first,
    int count,
    const u16** data,
    const u32** perm) {
    direct_reader_t* reader = opaque;
    size_t data_count;

    reader->calls++;
    if (reader->fail_on_call == reader->calls) {
        errno = EIO;
        return -1;
    }
    if (first < 0 || count <= 0 ||
        first > reader->kd->ndata - count ||
        (size_t)count > reader->capacity_points) {
        errno = ERANGE;
        return -1;
    }

    if (count > reader->largest_count) {
        reader->largest_count = count;
    }
    data_count = (size_t)count * (size_t)reader->kd->ndim;
    memcpy(reader->data,
           reader->kd->data.s +
               (size_t)first * (size_t)reader->kd->ndim,
           data_count * sizeof(u16));
    if (reader->kd->perm) {
        memcpy(reader->perm,
               reader->kd->perm + first,
               (size_t)count * sizeof(u32));
    }

    *data = reader->data;
    *perm = reader->kd->perm ? reader->perm : NULL;
    return 0;
}

static void direct_assert_result_equal(
    CuTest* ct,
    const kdtree_qres_t* expected,
    const kdtree_qres_t* actual) {
    size_t result_count;

    CuAssertPtrNotNull(ct, expected);
    CuAssertPtrNotNull(ct, actual);
    CuAssertIntEquals(ct, (int)expected->nres, (int)actual->nres);

    result_count = (size_t)expected->nres;
    CuAssert(
        ct,
        "direct search must preserve exact result index order",
        !memcmp(expected->inds,
                actual->inds,
                result_count * sizeof(expected->inds[0])));
    CuAssert(
        ct,
        "direct search must preserve squared-distance bits",
        !memcmp(expected->sdists,
                actual->sdists,
                result_count * sizeof(expected->sdists[0])));
    CuAssert(
        ct,
        "index-and-distance search must not materialize points",
        actual->results.any == NULL);
}

static void direct_compare_query(
    CuTest* ct,
    const kdtree_t* kd,
    direct_reader_t* reader,
    kdtree_qres_t** direct_reuse,
    const double* query,
    double maxd2,
    size_t max_points,
    size_t merge_gap_points) {
    const int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;
    kdtree_qres_t* legacy;
    kdtree_qres_t* direct;
    kdtree_qres_t* original_reuse = *direct_reuse;

    legacy = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        query,
        maxd2,
        options);
    CuAssertPtrNotNull(ct, legacy);

    reader->calls = 0;
    reader->largest_count = 0;
    reader->fail_on_call = 0;
    errno = 0;
    direct = kdtree_rangesearch_direct_dss(
        kd,
        *direct_reuse,
        query,
        maxd2,
        options,
        max_points,
        merge_gap_points,
        direct_reader_read,
        reader);
    CuAssertPtrNotNull(ct, direct);
    if (original_reuse) {
        CuAssert(
            ct,
            "direct search must preserve the caller's result container",
            direct == original_reuse);
    }
    CuAssert(
        ct,
        "direct reader request must honor max_points",
        reader->largest_count <= (int)max_points);
    direct_assert_result_equal(ct, legacy, direct);

    *direct_reuse = direct;
    kdtree_free_query(legacy);
}

static void direct_run_matrix(
    CuTest* ct,
    kdtree_t* kd,
    direct_reader_t* reader) {
    static const double fixed_queries[][4] = {
        {0.0, 1.0, 0.0, 1.0},
        {0.5, 0.5, 0.5, 0.5},
        {-0.001, 0.25, 0.75, 1.001},
        {1.0, 1.0, 1.0, 1.0},
        {0.125, 0.875, 0.375, 0.625}
    };
    static const double fixed_maxd2[] = {
        0.0,
        1.0e-12,
        0.0025,
        0.125,
        4.0
    };
    kdtree_qres_t* reuse = NULL;
    double query[4];
    uint32_t state = UINT32_C(0xc001d00d);
    int i;
    int d;

    for (i = 0; i < (int)(sizeof(fixed_queries) /
                           sizeof(fixed_queries[0])); i++) {
        direct_compare_query(
            ct,
            kd,
            reader,
            &reuse,
            fixed_queries[i],
            fixed_maxd2[i],
            (i & 1) ? 7U : 19U,
            (i & 1) ? 0U : 3U);
    }

    /*
     * Deterministic randomized queries exercise different leaf batching,
     * including max_points smaller than a leaf and a reusable result object.
     */
    for (i = 0; i < 96; i++) {
        double radius;

        for (d = 0; d < 4; d++) {
            uint32_t value = direct_test_prng(&state);
            query[d] = ((double)(value & UINT32_C(0xffff)) / 65535.0) *
                1.04 - 0.02;
        }
        radius = (double)(direct_test_prng(&state) & UINT32_C(0x1fff)) /
            65535.0;
        direct_compare_query(
            ct,
            kd,
            reader,
            &reuse,
            query,
            radius * radius,
            (i % 3 == 0) ? 3U : ((i % 3 == 1) ? 11U : 37U),
            (i % 4 == 0) ? 0U : 5U);
    }
    kdtree_free_query(reuse);
}

static void direct_test_tree(
    CuTest* ct,
    int build_options,
    int remove_permutation) {
    const int npoints = 513;
    const int ndim = 4;
    const int nleaf = 9;
    double* data;
    kdtree_t* kd;
    direct_reader_t reader;
    u32* saved_permutation = NULL;

    data = direct_test_points(npoints, ndim);
    CuAssertPtrNotNull(ct, data);
    kd = kdtree_build(
        NULL,
        data,
        npoints,
        ndim,
        nleaf,
        KDTT_DSS,
        build_options);
    free(data);
    CuAssertPtrNotNull(ct, kd);
    CuAssertIntEquals(ct, 0, kdtree_check(kd));

    if (remove_permutation) {
        saved_permutation = kd->perm;
        kd->perm = NULL;
    }

    CuAssertIntEquals(
        ct,
        0,
        direct_reader_init(&reader, kd, (size_t)npoints));
    direct_run_matrix(ct, kd, &reader);
    direct_reader_cleanup(&reader);

    kd->perm = saved_permutation ? saved_permutation : kd->perm;
    kdtree_free(kd);
}

void test_kdtree_direct_dss_splitdim_permutation(CuTest* ct) {
    direct_test_tree(
        ct,
        KD_BUILD_SPLIT |
            KD_BUILD_SPLITDIM |
            KD_BUILD_LINEAR_LR,
        0);
}

void test_kdtree_direct_dss_packed_split_no_permutation(CuTest* ct) {
    direct_test_tree(ct, KD_BUILD_SPLIT, 1);
}

void test_kdtree_direct_dss_callback_failure_resets_reuse(CuTest* ct) {
    const int options =
        KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;
    const int npoints = 513;
    const int ndim = 4;
    const double seed_query[4] = {0.10, 0.20, 0.30, 0.40};
    const double broad_query[4] = {0.50, 0.50, 0.50, 0.50};
    double* data;
    kdtree_t* kd;
    direct_reader_t reader;
    kdtree_qres_t* reuse;
    kdtree_qres_t* result;
    kdtree_qres_t* expected;

    data = direct_test_points(npoints, ndim);
    CuAssertPtrNotNull(ct, data);
    kd = kdtree_build(
        NULL,
        data,
        npoints,
        ndim,
        9,
        KDTT_DSS,
        KD_BUILD_SPLIT | KD_BUILD_SPLITDIM);
    free(data);
    CuAssertPtrNotNull(ct, kd);
    CuAssertIntEquals(
        ct,
        0,
        direct_reader_init(&reader, kd, (size_t)npoints));

    reuse = kdtree_rangesearch_direct_dss(
        kd,
        NULL,
        seed_query,
        0.05,
        options,
        3,
        0,
        direct_reader_read,
        &reader);
    CuAssertPtrNotNull(ct, reuse);

    reader.calls = 0;
    reader.fail_on_call = 2;
    errno = 0;
    result = kdtree_rangesearch_direct_dss(
        kd,
        reuse,
        broad_query,
        4.0,
        options,
        3,
        0,
        direct_reader_read,
        &reader);
    CuAssertPtrEquals(ct, NULL, result);
    CuAssertIntEquals(ct, EIO, errno);
    CuAssertIntEquals(ct, 0, (int)reuse->nres);

    expected = kdtree_rangesearch_options_reuse(
        kd,
        NULL,
        broad_query,
        4.0,
        options);
    CuAssertPtrNotNull(ct, expected);

    reader.calls = 0;
    reader.fail_on_call = 0;
    result = kdtree_rangesearch_direct_dss(
        kd,
        reuse,
        broad_query,
        4.0,
        options,
        3,
        0,
        direct_reader_read,
        &reader);
    CuAssert(
        ct,
        "reuse remains valid after callback failure",
        result == reuse);
    direct_assert_result_equal(ct, expected, result);

    kdtree_free_query(expected);
    kdtree_free_query(reuse);
    direct_reader_cleanup(&reader);
    kdtree_free(kd);
}
