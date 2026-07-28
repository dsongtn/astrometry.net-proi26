/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <unistd.h>

#include "solver.h"
#include "index.h"
#include "pquad.h"
#include "permutedsort.h"
#include "bl-sort.h"
#include "index_shard_internal.h"
#include "log.h"

extern int solver_test_codekd_cache_churn(void);

static int compare_n(const void* v1, const void* v2, int N) {
    const int* u1 = v1;
    const int* u2 = v2;
    int i;
    for (i=0; i<N; i++) {
        if (u1[i] < u2[i]) return -1;
        if (u1[i] > u2[i]) return 1;
    }
    return 0;
}

static int compare_tri(const void* v1, const void* v2) {
    return compare_n(v1, v2, 3);
}
static int compare_quad(const void* v1, const void* v2) {
    return compare_n(v1, v2, 4);
}
static int compare_quint(const void* v1, const void* v2) {
    return compare_n(v1, v2, 5);
}




bl* quadlist;

void test_try_all_codes(pquad* pq,
                        unsigned int* fieldstars, int dimquad,
                        solver_t* solver, double tol2) {
    int sorted[dimquad];
    int i;
    fflush(NULL);
    printf("test_try_all_codes: [");
    for (i=0; i<dimquad; i++) {
        printf("%s%i", (i?" ":""), fieldstars[i]);
    }
    printf("]");

    // sort AB and C[DE]...
    memcpy(sorted, fieldstars, dimquad * sizeof(int));
    qsort(sorted, 2, sizeof(int), compare_ints_asc);
    qsort(sorted+2, dimquad-2, sizeof(int), compare_ints_asc);

    printf(" -> [");
    for (i=0; i<dimquad; i++) {
        printf("%s%i", (i?" ":""), sorted[i]);
    }
    printf("]\n");
    fflush(NULL);

    bl_append(quadlist, sorted);
}

static starxy_t* field1() {
    starxy_t* starxy;
    double field[14];
    int i=0, N;
    // star0 A: (0,0)
    field[i++] = 0.0;
    field[i++] = 0.0;
    // star1 B: (2,2)
    field[i++] = 2.0;
    field[i++] = 2.0;
    // star2
    field[i++] = -1.0;
    field[i++] = 3.0;
    // star3
    field[i++] = 0.5;
    field[i++] = 1.5;
    // star4
    field[i++] = 1.0;
    field[i++] = 1.0;
    // star5
    field[i++] = 1.5;
    field[i++] = 0.5;
    // star6
    field[i++] = 3.0;
    field[i++] = -1.0;

    N = i/2;
    starxy = starxy_new(N, FALSE, FALSE);
    for (i=0; i<N; i++) {
        starxy_setx(starxy, i, field[i*2+0]);
        starxy_sety(starxy, i, field[i*2+1]);
    }
    return starxy;
}

static bl* collect_frontier_quads(
    int dimquads,
    int startobj,
    anbool use_geometry_cache) {
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    bl* collected;

    starxy = field1();
    collected = bl_new(16, (size_t)dimquads * sizeof(uint));
    quadlist = collected;

    solver = solver_new();
    memset(&index, 0, sizeof(index_t));
    index.index_scale_lower = 1;
    index.index_scale_upper = 10;
    index.dimquads = dimquads;

    solver->funits_lower = 0.1;
    solver->funits_upper = 10;
    solver->startobj = startobj;
    solver->endobj = starxy_n(starxy);

    solver_add_index(solver, &index);
    solver_set_field(solver, starxy);
    solver_preprocess_field(solver);
    if (use_geometry_cache) {
        assert(solver_prepare_field_geometry(solver));
    }

    assert(!solver_run(solver));
    solver_free_field(solver);
    solver_free(solver);
    return collected;
}

static void test_geometry_cache_exact(void) {
    static const int starts[] = {0, 1, 3, 5};
    int dimquads;
    int start_index;

    for (dimquads = 3; dimquads <= 5; dimquads++) {
        for (start_index = 0;
             start_index < (int)(sizeof(starts) / sizeof(starts[0]));
             start_index++) {
            int startobj = starts[start_index];
            bl* legacy =
                collect_frontier_quads(dimquads, startobj, FALSE);
            bl* cached =
                collect_frontier_quads(dimquads, startobj, TRUE);
            int i;

            assert(bl_size(legacy) == bl_size(cached));
            for (i = 0; i < bl_size(legacy); i++) {
                assert(!memcmp(
                    bl_access(legacy, i),
                    bl_access(cached, i),
                    (size_t)dimquads * sizeof(uint)));
            }
            bl_free(legacy);
            bl_free(cached);
        }
    }
}

static void test_geometry_cache_deep_admission(void) {
    enum { STAR_COUNT = 200 };
    solver_t* solver;
    starxy_t* starxy;
    int i;

    starxy = starxy_new(STAR_COUNT, FALSE, FALSE);
    assert(starxy);
    for (i = 0; i < STAR_COUNT; i++) {
        starxy_setx(starxy, i, (double)(i % 20));
        starxy_sety(starxy, i, (double)(i / 20));
    }

    solver = solver_new();
    solver->startobj = 7;
    solver->endobj = STAR_COUNT;
    solver_set_field(solver, starxy);
    solver_preprocess_field(solver);

    /*
     * V19's per-pair per-star payload crossed the 64 MiB cache budget
     * around this field size. V20's triangular transform table is O(n^2)
     * and must remain admitted for the same later-band field.
     */
    assert(solver_prepare_field_geometry(solver));
    assert(solver->field_geometry);

    solver_free_field(solver);
    solver_free(solver);
}

void test1() {
    int i;
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    int wanted[][4] = { { 0,1,3,4 },
                        { 0,2,3,4 },
                        { 1,2,3,4 },
                        { 2,5,0,1 },
                        { 2,5,0,3 },
                        { 2,5,0,4 },
                        { 2,5,1,3 },
                        { 2,5,1,4 },
                        { 2,5,3,4 },
                        { 0,1,3,5 },
                        { 0,1,4,5 },
                        { 0,6,4,5 },
                        { 1,6,4,5 },
                        { 2,6,0,1 },
                        { 2,6,0,3 },
                        { 2,6,0,4 },
                        { 2,6,0,5 },
                        { 2,6,1,3 },
                        { 2,6,1,4 },
                        { 2,6,1,5 },
                        { 2,6,3,4 },
                        { 2,6,3,5 },
                        { 2,6,4,5 },
                        { 3,6,0,1 },
                        { 3,6,0,4 },
                        { 3,6,0,5 },
                        { 3,6,1,4 },
                        { 3,6,1,5 },
                        { 3,6,4,5 },
    };

    starxy = field1();

    quadlist = bl_new(16, 4*sizeof(uint));

    solver = solver_new();

    memset(&index, 0, sizeof(index_t));
    index.index_scale_lower = 1;
    index.index_scale_upper = 10;
    index.dimquads = 4;

    solver->funits_lower = 0.1;
    solver->funits_upper = 10;

    solver_add_index(solver, &index);
    solver_set_field(solver, starxy);
    solver_preprocess_field(solver);

    solver_run(solver);

    solver_free_field(solver);
    solver_free(solver);

    //
    assert(bl_size(quadlist) == (sizeof(wanted) / (4*sizeof(uint))));
    for (i=0; i<bl_size(quadlist); i++) {
        assert(compare_quad(bl_access(quadlist, i), wanted[i]) == 0);
    }

    bl_free(quadlist);
}

void test2() {
    int i;
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    int wanted[][3] = { { 0, 1, 3 },
                        { 0, 1, 4 },
                        { 0, 1, 5 },
                        { 0, 2, 3 },
                        { 0, 2, 4 },
                        { 0, 3, 4 },
                        { 0, 5, 4 },
                        { 0, 6, 4 },
                        { 0, 6, 5 },
                        { 1, 2, 3 },
                        { 1, 2, 4 },
                        { 1, 3, 4 },
                        { 1, 5, 4 },
                        { 1, 6, 4 },
                        { 1, 6, 5 },
                        { 2, 4, 3 },
                        { 2, 5, 0 },
                        { 2, 5, 1 },
                        { 2, 5, 3 },
                        { 2, 5, 4 },
                        { 2, 6, 0 },
                        { 2, 6, 1 },
                        { 2, 6, 3 },
                        { 2, 6, 4 },
                        { 2, 6, 5 },
                        { 3, 5, 4 },
                        { 3, 6, 0 },
                        { 3, 6, 1 },
                        { 3, 6, 4 },
                        { 3, 6, 5 },
                        { 4, 6, 5 },
    };

    starxy = field1();
    quadlist = bl_new(16, 3*sizeof(uint));
    solver = solver_new();
    memset(&index, 0, sizeof(index_t));
    index.index_scale_lower = 1;
    index.index_scale_upper = 10;
    index.dimquads = 3;

    solver->funits_lower = 0.1;
    solver->funits_upper = 10;

    solver_add_index(solver, &index);
    solver_set_field(solver, starxy);
    solver_preprocess_field(solver);

    solver_run(solver);

    solver_free_field(solver);
    solver_free(solver);

    //
    assert(bl_size(quadlist) == (sizeof(wanted) / (3*sizeof(uint))));
    bl_sort(quadlist, compare_tri);
    for (i=0; i<bl_size(quadlist); i++) {
        assert(compare_tri(bl_access(quadlist, i), wanted[i]) == 0);
    }
    bl_free(quadlist);
}


char* OPTIONS = "v";

static void test_ab_counter_boundaries(void) {
    solver_t solver;
    int cxdx_before;
    int meanx_before;

    memset(&solver, 0, sizeof(solver));
    assert(solver_ab_checked_counter_delta(
        &solver,
        3ULL,
        5ULL,
        7ULL) == 0);
    assert(solver.numtries == 3);
    assert(solver.num_cxdx_skipped == 5);
    assert(solver.num_meanx_skipped == 7);

    solver.numtries = INT_MAX;
    cxdx_before = solver.num_cxdx_skipped;
    meanx_before = solver.num_meanx_skipped;
    assert(solver_ab_checked_counter_delta(
        &solver,
        1ULL,
        1ULL,
        1ULL) != 0);
    assert(solver.numtries == INT_MAX);
    assert(solver.num_cxdx_skipped == cxdx_before);
    assert(solver.num_meanx_skipped == meanx_before);
    assert(solver.profile.execution_failed);
    assert(solver.quit_now);

    memset(&solver, 0, sizeof(solver));
    solver.numtries = INT_MAX - 1;
    assert(solver_ab_checked_counter_delta(
        &solver,
        1ULL,
        0ULL,
        0ULL) == 0);
    assert(solver.numtries == INT_MAX);

    memset(&solver, 0, sizeof(solver));
    solver.num_cxdx_skipped = INT_MAX - 1;
    assert(solver_ab_checked_counter_delta(
        &solver,
        2ULL,
        2ULL,
        0ULL) != 0);
    assert(solver.numtries == 0);
    assert(solver.num_cxdx_skipped == INT_MAX - 1);
    assert(solver.num_meanx_skipped == 0);

    memset(&solver, 0, sizeof(solver));
    solver.num_meanx_skipped = INT_MAX;
    assert(solver_ab_checked_counter_delta(
        &solver,
        1ULL,
        1ULL,
        1ULL) != 0);
    assert(solver.numtries == 0);
    assert(solver.num_cxdx_skipped == 0);
    assert(solver.num_meanx_skipped == INT_MAX);
}

static void test_index_close_fds_failure_state(void) {
    index_t index;
    quadfile_t quads;
    codetree_t codes;
    startree_t stars;
    kdtree_t code_tree;
    kdtree_t star_tree;
    fitsbin_t quad_fits;
    fitsbin_t code_fits;
    fitsbin_t star_fits;

    memset(&index, 0, sizeof(index));
    memset(&quads, 0, sizeof(quads));
    memset(&codes, 0, sizeof(codes));
    memset(&stars, 0, sizeof(stars));
    memset(&code_tree, 0, sizeof(code_tree));
    memset(&star_tree, 0, sizeof(star_tree));
    memset(&quad_fits, 0, sizeof(quad_fits));
    memset(&code_fits, 0, sizeof(code_fits));
    memset(&star_fits, 0, sizeof(star_fits));

    quad_fits.fid = tmpfile();
    code_fits.fid = tmpfile();
    star_fits.fid = tmpfile();
    assert(quad_fits.fid);
    assert(code_fits.fid);
    assert(star_fits.fid);
    assert(close(fileno(quad_fits.fid)) == 0);

    quads.fb = &quad_fits;
    code_tree.io = &code_fits;
    star_tree.io = &star_fits;
    codes.tree = &code_tree;
    stars.tree = &star_tree;
    index.quads = &quads;
    index.codekd = &codes;
    index.starkd = &stars;

    /*
     * fclose() reports EBADF for the deliberately invalid first stream.
     * Every fitsbin must nevertheless relinquish its invalid FILE* and later
     * components must still be closed.
     */
    printf("test_index_close_fds_failure_state: "
           "begin expected EBADF injection\n");
    assert(index_close_fds(&index) != 0);
    assert(!quad_fits.fid);
    assert(!code_fits.fid);
    assert(!star_fits.fid);
    printf("test_index_close_fds_failure_state: PASS\n");
}

static void test_zero_initialized_payload_fd_is_unowned(void) {
    fitsbin_t fits;
    int opened_stdin = 0;

    memset(&fits, 0, sizeof(fits));
    errno = 0;
    if (fcntl(STDIN_FILENO, F_GETFD) < 0 &&
        errno == EBADF) {
        assert(open("/dev/null", O_RDONLY) ==
               STDIN_FILENO);
        opened_stdin = 1;
    }
    assert(fitsbin_close_payload_fd(&fits) == 0);
    assert(fcntl(STDIN_FILENO, F_GETFD) >= 0);
    if (opened_stdin) {
        assert(close(STDIN_FILENO) == 0);
    }
}

int main(int argc, char** args) {
    int argchar;

    while ((argchar = getopt(argc, args, OPTIONS)) != -1)
        switch (argchar) {
        case 'v':
            log_init(LOG_ALL+1);
            break;
        }

    test1();
    test2();
    test_geometry_cache_exact();
    test_geometry_cache_deep_admission();
    test_ab_counter_boundaries();
    assert(!solver_test_codekd_cache_churn());
    assert(!onefield_job_index_cache_test_release_state());
    test_zero_initialized_payload_fd_is_unowned();
    test_index_close_fds_failure_state();
    return 0;
}
