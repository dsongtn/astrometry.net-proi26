/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "solver.h"
#include "index.h"
#include "bl-sort.h"
#include "log.h"
#include "onefield_test_private.h"
#include "test_solver_private.h"

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



void test1() {
    int i;
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    bl* quadlist;
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

    starxy = test_solver_geometry_field();

    quadlist = bl_new(16, 4*sizeof(uint));
    test_solver_geometry_use_quadlist(quadlist);

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
    test_solver_geometry_use_quadlist(NULL);
}

void test2() {
    int i;
    solver_t* solver;
    index_t index;
    starxy_t* starxy;
    bl* quadlist;
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

    starxy = test_solver_geometry_field();
    quadlist = bl_new(16, 3*sizeof(uint));
    test_solver_geometry_use_quadlist(quadlist);
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
    test_solver_geometry_use_quadlist(NULL);
}


char* OPTIONS = "v";

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
    test_solver_geometry_cache_exact();
    test_solver_geometry_cache_deep_admission();
    test_solver_ab_counter_boundaries();
    test_solver_ab_descriptor_partition_count();
    assert(!onefield_job_index_cache_test_handoff_state());
    test_solver_zero_initialized_payload_fd_is_unowned();
    test_solver_index_close_fds_failure_state();
    return 0;
}
