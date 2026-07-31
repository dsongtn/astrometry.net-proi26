/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <string.h>
#include <unistd.h>

#include "cutest.h"
#include "index.h"
#include "starutil.h"

static const char* find_index_fixture(const char* from_root,
                                      const char* from_util) {
    if (access(from_root, R_OK) == 0) {
        return from_root;
    }
    if (access(from_util, R_OK) == 0) {
        return from_util;
    }
    return NULL;
}

static void assert_metadata_equal(CuTest* ct,
                                  const index_t* expected,
                                  const index_t* actual) {
    CuAssertIntEquals(ct, expected->indexid, actual->indexid);
    CuAssertIntEquals(ct, expected->healpix, actual->healpix);
    CuAssertIntEquals(ct, expected->hpnside, actual->hpnside);
    CuAssertDblEquals(ct, expected->index_jitter,
                      actual->index_jitter, 0.0);
    CuAssertIntEquals(ct, expected->cutnside, actual->cutnside);
    CuAssertIntEquals(ct, expected->cutnsweep, actual->cutnsweep);
    CuAssertDblEquals(ct, expected->cutdedup, actual->cutdedup, 0.0);
    CuAssertIntEquals(ct, expected->cutmargin, actual->cutmargin);
    CuAssertIntEquals(ct, expected->circle, actual->circle);
    CuAssertIntEquals(ct, expected->cx_less_than_dx,
                      actual->cx_less_than_dx);
    CuAssertIntEquals(ct, expected->meanx_less_than_half,
                      actual->meanx_less_than_half);
    CuAssertDblEquals(ct, expected->index_scale_upper,
                      actual->index_scale_upper, 0.0);
    CuAssertDblEquals(ct, expected->index_scale_lower,
                      actual->index_scale_lower, 0.0);
    CuAssertIntEquals(ct, expected->dimquads, actual->dimquads);
    CuAssertIntEquals(ct, expected->nstars, actual->nstars);
    CuAssertIntEquals(ct, expected->nquads, actual->nquads);

    CuAssertIntEquals(ct, expected->cutband == NULL,
                      actual->cutband == NULL);
    if (expected->cutband) {
        CuAssertStrEquals(ct, expected->cutband, actual->cutband);
    }
    CuAssertStrEquals(ct, expected->indexname, actual->indexname);
    CuAssertStrEquals(ct, expected->indexfn, actual->indexfn);
}

static void check_fixture(CuTest* ct, const char* filename) {
    double code[DCMAX];
    double xyz[3];
    unsigned int stars[DQMAX];
    index_t* full;
    index_t* metadata;

    full = index_load(filename, 0, NULL);
    CuAssertPtrNotNull(ct, full);

    metadata = index_load(filename, INDEX_ONLY_LOAD_METADATA, NULL);
    CuAssertPtrNotNull(ct, metadata);
    assert_metadata_equal(ct, full, metadata);
    CuAssertPtrNotNull(ct, metadata->fits);
    CuAssertPtrEquals(ct, NULL, metadata->starkd);
    CuAssertPtrEquals(ct, NULL, metadata->quads);
    CuAssertPtrEquals(ct, NULL, metadata->codekd);

    CuAssertIntEquals(ct, 0, index_reload(metadata));
    CuAssertPtrNotNull(ct, metadata->starkd);
    CuAssertPtrNotNull(ct, metadata->quads);
    CuAssertPtrNotNull(ct, metadata->codekd);
    assert_metadata_equal(ct, full, metadata);
    CuAssertIntEquals(ct, startree_N(full->starkd),
                      startree_N(metadata->starkd));
    CuAssertIntEquals(ct, codetree_N(full->codekd),
                      codetree_N(metadata->codekd));

    CuAssertIntEquals(ct, 0,
                      quadfile_get_stars(metadata->quads, 0, stars));
    CuAssertIntEquals(ct, 0,
                      startree_get(metadata->starkd, stars[0], xyz));
    CuAssertIntEquals(ct, 0, codetree_get(metadata->codekd, 0, code));

    index_free(metadata);
    index_free(full);
}

void test_index_metadata_header_only(CuTest* ct) {
    const char* fixture_4119 =
        find_index_fixture("demo/index-4119.fits",
                           "../demo/index-4119.fits");
    const char* fixture_9918 =
        find_index_fixture("solver/index-9918.fits",
                           "../solver/index-9918.fits");

    CuAssertPtrNotNull(ct, fixture_4119);
    CuAssertPtrNotNull(ct, fixture_9918);
    check_fixture(ct, fixture_4119);
    check_fixture(ct, fixture_9918);
}
