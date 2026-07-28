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
#include <errno.h>

#include "cutest.h"
#include "kdtree.h"
#include "kdtree_fits_io.h"

#include "test_libkd_common.c"

static void assert_kdtrees_equal(CuTest* ct, const kdtree_t* kd, const kdtree_t* kd2) {
    double del = 1e-10;
    size_t sz, sz2;

    if (!kd) {
        CuAssertPtrEquals(ct, NULL, (kdtree_t*)kd2);
        return;
    }
    CuAssertPtrNotNull(ct, kd2);

    CuAssertIntEquals(ct, kd->treetype, kd2->treetype);
    CuAssertIntEquals(ct, kd->dimbits, kd2->dimbits);
    CuAssertIntEquals(ct, kd->dimmask, kd2->dimmask);
    CuAssertIntEquals(ct, kd->splitmask, kd2->splitmask);
    CuAssertIntEquals(ct, kd->ndata, kd2->ndata);
    CuAssertIntEquals(ct, kd->ndim, kd2->ndim);
    CuAssertIntEquals(ct, kd->nnodes, kd2->nnodes);
    CuAssertIntEquals(ct, kd->nbottom, kd2->nbottom);
    CuAssertIntEquals(ct, kd->ninterior, kd2->ninterior);
    CuAssertIntEquals(ct, kd->nlevels, kd2->nlevels);
    CuAssertIntEquals(ct, kd->has_linear_lr, kd2->has_linear_lr);
    CuAssertDblEquals(ct, kd->scale,    kd2->scale,    del);
    CuAssertDblEquals(ct, kd->invscale, kd2->invscale, del);

    if (kd->lr) {
        CuAssertPtrNotNull(ct, kd2->lr);
        sz  = kdtree_sizeof_lr(kd );
        sz2 = kdtree_sizeof_lr(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "lr equal", memcmp(kd->lr, kd2->lr, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->lr);
    }

    if (kd->perm) {
        CuAssertPtrNotNull(ct, kd2->perm);
        sz  = kdtree_sizeof_perm(kd );
        sz2 = kdtree_sizeof_perm(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "perm equal", memcmp(kd->perm, kd2->perm, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->perm);
    }

    if (kd->data.any) {
        CuAssertPtrNotNull(ct, kd2->data.any);
        sz  = kdtree_sizeof_data(kd );
        sz2 = kdtree_sizeof_data(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "data equal", memcmp(kd->data.any, kd2->data.any, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->data.any);
    }

    if (kd->splitdim) {
        CuAssertPtrNotNull(ct, kd2->splitdim);
        sz  = kdtree_sizeof_splitdim(kd );
        sz2 = kdtree_sizeof_splitdim(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "splitdim equal", memcmp(kd->splitdim, kd2->splitdim, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->splitdim);
    }

    if (kd->split.any) {
        CuAssertPtrNotNull(ct, kd2->split.any);
        sz  = kdtree_sizeof_split(kd );
        sz2 = kdtree_sizeof_split(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "split equal", memcmp(kd->split.any, kd2->split.any, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->split.any);
    }

    if (kd->bb.any) {
        CuAssertPtrNotNull(ct, kd2->bb.any);
        sz  = kdtree_sizeof_bb(kd );
        sz2 = kdtree_sizeof_bb(kd2);
        CuAssertIntEquals(ct, sz, sz2);
        CuAssert(ct, "bb equal", memcmp(kd->bb.any, kd2->bb.any, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->bb.any);
    }

    if (kd->minval) {
        sz  = kd->ndim * sizeof(double);
        CuAssertPtrNotNull(ct, kd2->minval);
        CuAssert(ct, "minval equal", memcmp(kd->minval, kd2->minval, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->minval);
    }

    if (kd->maxval) {
        CuAssertPtrNotNull(ct, kd2->maxval);
        sz  = kd->ndim * sizeof(double);
        CuAssert(ct, "maxval equal", memcmp(kd->maxval, kd2->maxval, sz) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->maxval);
    }

    if (kd->name) {
        CuAssertPtrNotNull(ct, kd2->name);
        CuAssert(ct, "name equal", strcmp(kd->name, kd2->name) == 0);
    } else {
        CuAssertPtrEquals(ct, NULL, kd2->name);
    }
}

void test_read_write_single_tree_unnamed(CuTest* ct) {
    kdtree_t* kd;
    double * data;
    int N = 1000;
    int Nleaf = 5;
    int D = 3;
    char fn[1024];
    int rtn;
    kdtree_t* kd2;
    int fd;

    data = random_points_d(N, D);
    kd = build_tree(ct, data, N, D, Nleaf, KDTT_DOUBLE, KD_BUILD_SPLIT);
    kd->name = NULL;

    sprintf(fn, "/tmp/test_libkd_io_single_tree_unnamed.XXXXXX");
    fd = mkstemp(fn);
    if (fd == -1) {
        fprintf(stderr, "Failed to generate a temp filename: %s\n", strerror(errno));
        CuFail(ct, "mkstemp");
    }
    close(fd);
    printf("Single tree unnamed: writing to file %s.\n", fn);

    rtn = kdtree_fits_write(kd, fn, NULL);
    CuAssertIntEquals(ct, 0, rtn);

    kd2 = kdtree_fits_read(fn, NULL, NULL);

    assert_kdtrees_equal(ct, kd, kd2);

    free(data);
    kdtree_free(kd);
    kdtree_fits_close(kd2);
}

void test_read_write_single_tree_named(CuTest* ct) {
    kdtree_t* kd;
    double * data;
    int N = 1000;
    int Nleaf = 5;
    int D = 3;
    char fn[1024];
    int rtn;
    kdtree_t* kd2;
    int fd;

    data = random_points_d(N, D);
    kd = build_tree(ct, data, N, D, Nleaf, KDTT_DOUBLE,
                    KD_BUILD_SPLIT | KD_BUILD_BBOX | KD_BUILD_LINEAR_LR);
    kd->name = strdup("christmas");

    sprintf(fn, "/tmp/test_libkd_io_single_tree_named.XXXXXX");
    fd = mkstemp(fn);
    if (fd == -1) {
        fprintf(stderr, "Failed to generate a temp filename: %s\n", strerror(errno));
        CuFail(ct, "mkstemp");
    }
    close(fd);
    printf("Single tree named: writing to file %s.\n", fn);

    rtn = kdtree_fits_write(kd, fn, NULL);
    CuAssertIntEquals(ct, 0, rtn);

    // Loading any tree should succeed.
    kd2 = kdtree_fits_read(fn, NULL, NULL);
    assert_kdtrees_equal(ct, kd, kd2);
    kdtree_fits_close(kd2);

    // Attempting to load a nonexist named tree should fail.
    kd2 = kdtree_fits_read(fn, "none", NULL);
    CuAssertPtrEquals(ct, NULL, kd2);

    // Loading by its correct name should work.
    kd2 = kdtree_fits_read(fn, "christmas", NULL);
    assert_kdtrees_equal(ct, kd, kd2);
    kdtree_fits_close(kd2);

    free(data);
    kdtree_free(kd);
}

static anbool test_chunk_name_is(
    const char* tablename,
    const char* component) {
    size_t component_length;

    if (!tablename || !component) {
        return FALSE;
    }

    component_length = strlen(component);
    return !strncmp(tablename, component, component_length) &&
        (tablename[component_length] == '\0' ||
         tablename[component_length] == '_');
}

void test_production_mmap_component_classification(CuTest* ct) {
    kdtree_t* kd;
    kdtree_t* mapped_kd;
    kdtree_fits_t* io;
    fitsbin_t* fb;
    double* data;
    char fn[1024];
    int fd;
    int i;
    int payload_chunks = 0;
    int topology_chunks = 0;
    int saw_data = 0;
    int saw_split = 0;

    data = random_points_d(1000, 3);
    kd = build_tree(
        ct,
        data,
        1000,
        3,
        5,
        KDTT_DSS,
        KD_BUILD_SPLIT |
            KD_BUILD_SPLITDIM |
            KD_BUILD_BBOX |
            KD_BUILD_LINEAR_LR);
    kd->name = strdup("production");
    CuAssertPtrNotNull(ct, kd->name);

    sprintf(fn, "/tmp/test_libkd_io_mmap_components.XXXXXX");
    fd = mkstemp(fn);
    if (fd == -1) {
        free(data);
        kdtree_free(kd);
        CuFail(ct, "mkstemp");
        return;
    }
    close(fd);

    CuAssertIntEquals(ct, 0, kdtree_fits_write(kd, fn, NULL));

    /*
     * Use a controlled RANDOM mapping policy. The reader still assigns
     * distinct region classes: sparse DATA/PERM payload stays RANDOM while
     * repeatedly traversed topology stays NORMAL.
     */
    fitsbin_mmap_set_thread_advice(FITSBIN_MMAP_ADVICE_RANDOM);
    io = kdtree_fits_open(fn);
    CuAssertPtrNotNull(ct, io);
    mapped_kd = kdtree_fits_read_tree(io, "production", NULL);
    fitsbin_mmap_clear_thread_advice();
    CuAssertPtrNotNull(ct, mapped_kd);

    fb = kdtree_fits_get_fitsbin(io);
    for (i = 0; i < fitsbin_n_chunks(fb); i++) {
        fitsbin_chunk_t* chunk = fitsbin_get_chunk(fb, i);

        CuAssertPtrNotNull(ct, chunk);
        CuAssertPtrNotNull(ct, chunk->tablename);

        if (test_chunk_name_is(chunk->tablename, KD_STR_DATA) ||
            test_chunk_name_is(chunk->tablename, KD_STR_PERM)) {
            payload_chunks++;
            if (test_chunk_name_is(chunk->tablename, KD_STR_DATA)) {
                saw_data = 1;
            }
            CuAssertIntEquals(
                ct,
                FITSBIN_MMAP_REGION_PAYLOAD,
                chunk->mmap_region);
            CuAssertIntEquals(
                ct,
                FITSBIN_MMAP_ADVICE_RANDOM,
                fitsbin_get_chunk_mmap_advice(fb, chunk));
        } else {
            topology_chunks++;
            if (test_chunk_name_is(chunk->tablename, KD_STR_SPLIT)) {
                saw_split = 1;
            }
            CuAssertIntEquals(
                ct,
                FITSBIN_MMAP_REGION_TOPOLOGY,
                chunk->mmap_region);
            CuAssertIntEquals(
                ct,
                FITSBIN_MMAP_ADVICE_NORMAL,
                fitsbin_get_chunk_mmap_advice(fb, chunk));
        }
    }

    CuAssert(ct, "production data payload chunk", saw_data);
    CuAssert(ct, "production split topology chunk", saw_split);
    CuAssert(ct, "production payload chunks", payload_chunks >= 1);
    CuAssert(ct, "production topology chunks", topology_chunks >= 1);

    kdtree_fits_close(mapped_kd);
    unlink(fn);
    free(data);
    kdtree_free(kd);
}

void test_read_write_two_trees(CuTest* ct) {
    kdtree_t* kd;
    kdtree_t* kdB;
    double * data;
    double * dataB;
    int N = 1000;
    int Nleaf = 5;
    int D = 3;
    char fn[1024];
    int rtn;
    kdtree_t* kd2;
    kdtree_t* kd2B;
    int fd;
    kdtree_fits_t* io;

    data = random_points_d(N, D);
    kd = build_tree(ct, data, N, D, Nleaf, KDTT_DOUBLE,
                    KD_BUILD_SPLIT | KD_BUILD_BBOX | KD_BUILD_LINEAR_LR);
    kd->name = strdup("christmas");

    dataB = random_points_d(N, D);
    kdB = build_tree(ct, dataB, N, D, Nleaf, KDTT_DUU,
                     KD_BUILD_SPLIT | KD_BUILD_SPLITDIM | KD_BUILD_LINEAR_LR);
    kdB->name = strdup("watermelon");

    sprintf(fn, "/tmp/test_libkd_io_two_trees.XXXXXX");
    fd = mkstemp(fn);
    if (fd == -1) {
        fprintf(stderr, "Failed to generate a temp filename: %s\n", strerror(errno));
        CuFail(ct, "mkstemp");
    }
    printf("Two trees: writing to file %s.\n", fn);

    close(fd);

    io = kdtree_fits_open_for_writing(fn);
    if (!io) {
        fprintf(stderr, "Failed to open temp file: %s\n", strerror(errno));
        CuFail(ct, "fdopen");
    }
    rtn = kdtree_fits_write_primary_header(io, NULL);
    CuAssertIntEquals(ct, 0, rtn);
    rtn = kdtree_fits_append_tree(io, kd, NULL);
    CuAssertIntEquals(ct, 0, rtn);
    rtn = kdtree_fits_append_tree(io, kdB, NULL);
    CuAssertIntEquals(ct, 0, rtn);

    if (kdtree_fits_io_close(io)) {
        fprintf(stderr, "Failed to close temp file: %s\n", strerror(errno));
        CuFail(ct, "fclose");
    }

    // Loading any tree should return the first one.
    kd2 = kdtree_fits_read(fn, NULL, NULL);
    assert_kdtrees_equal(ct, kd, kd2);
    kdtree_fits_close(kd2);

    // Attempting to load a nonexist named tree should fail.
    kd2 = kdtree_fits_read(fn, "none", NULL);
    CuAssertPtrEquals(ct, NULL, kd2);

    // Loading by the correct names should work.
    kd2 = kdtree_fits_read(fn, "christmas", NULL);
    assert_kdtrees_equal(ct, kd, kd2);

    kd2B = kdtree_fits_read(fn, "watermelon", NULL);
    assert_kdtrees_equal(ct, kdB, kd2B);

    kdtree_fits_close(kd2);
    kdtree_fits_close(kd2B);

    free(data);
    kdtree_free(kd);

    free(dataB);
    kdtree_free(kdB);
}
