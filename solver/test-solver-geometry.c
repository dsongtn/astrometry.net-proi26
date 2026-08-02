/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "solver.h"
#include "index.h"
#include "pquad.h"
#include "permutedsort.h"
#include "test_solver_private.h"

static bl* collected_quads;

starxy_t* test_solver_geometry_field(void) {
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

void test_solver_geometry_use_quadlist(bl* quadlist) {
    collected_quads = quadlist;
}

void test_try_all_codes(pquad* pq,
                        int* fieldstars, int dimquad,
                        solver_t* solver, double tol2) {
    int sorted[dimquad];
    int i;

    (void)pq;
    (void)solver;
    (void)tol2;
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

    bl_append(collected_quads, sorted);
}

static bl* collect_frontier_quads(
    int dimquads,
    int startobj,
    anbool use_geometry_cache) {
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    bl* collected;

    starxy = test_solver_geometry_field();
    collected = bl_new(16, (size_t)dimquads * sizeof(uint));
    test_solver_geometry_use_quadlist(collected);

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
    test_solver_geometry_use_quadlist(NULL);
    return collected;
}

void test_solver_geometry_cache_exact(void) {
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

void test_solver_geometry_cache_deep_admission(void) {
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
     * The former per-pair per-star payload crossed the 64 MiB cache budget
     * around this field size. The triangular transform table is O(n^2) and
     * must remain admitted for the same later-band field.
     */
    assert(solver_prepare_field_geometry(solver));
    assert(solver->field_geometry);

    solver_free_field(solver);
    solver_free(solver);
}
