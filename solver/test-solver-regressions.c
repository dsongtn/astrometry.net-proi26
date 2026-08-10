/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "solver.h"
#include "index.h"
#include "solver_hypothesis_internal.h"
#include "solver_inline_internal.h"
#include "test_solver_private.h"

void test_solver_ab_counter_boundaries(void) {
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

void test_solver_ab_descriptor_partition_count(void) {
    assert(solver_ab_descriptor_partition_count(
        0ULL, 2U, 100U, 4U) == 0U);
    assert(solver_ab_descriptor_partition_count(
        1ULL, 0U, 100U, 4U) == 0U);
    assert(solver_ab_descriptor_partition_count(
        1ULL, 2U, 0U, 4U) == 0U);
    assert(solver_ab_descriptor_partition_count(
        1ULL, 2U, 100U, 0U) == 0U);
    assert(solver_ab_descriptor_partition_count(
        ULLONG_MAX, 2U, 100U, 4U) == 0U);

    assert(solver_ab_descriptor_partition_count(
        64ULL, 2U, 100U, 8U) == 1U);
    assert(solver_ab_descriptor_partition_count(
        65ULL, 2U, 100U, 8U) == 2U);
    assert(solver_ab_descriptor_partition_count(
        256ULL, 2U, 256U, 8U) == 4U);
    assert(solver_ab_descriptor_partition_count(
        256ULL, 2U, 256U, 4U) == 4U);

    assert(solver_ab_descriptor_partition_count(
        300ULL, 1U, 100U, 8U) == 3U);
    assert(solver_ab_descriptor_partition_count(
        400ULL, 1U, 100U, 4U) == 4U);
    assert(solver_ab_descriptor_partition_count(
        401ULL, 1U, 100U, 4U) == 0U);

    assert(solver_ab_descriptor_lead_combinations(
        0ULL, 2U) == 0ULL);
    assert(solver_ab_descriptor_lead_combinations(
        64ULL, 0U) == 0ULL);
    assert(solver_ab_descriptor_lead_combinations(
        63ULL, 2U) == 0ULL);
    assert(solver_ab_descriptor_lead_combinations(
        64ULL, 2U) == 32ULL);
    assert(solver_ab_descriptor_lead_combinations(
        128ULL, 1U) == 64ULL);
    assert(solver_ab_descriptor_lead_combinations(
        129ULL, 1U) == 64ULL);
    assert(solver_ab_descriptor_lead_combinations(
        1ULL, 128U) == 0ULL);
    assert(solver_ab_descriptor_lead_combinations(
        2ULL, 128U) == 1ULL);
}

static void assert_descriptor_scientific_equal(
    const solver_ab_descriptor_t* expected,
    const solver_ab_descriptor_t* actual,
    int dimquads) {
    size_t code_count = (size_t)(dimquads - NBACK) * 2U;

    assert(expected);
    assert(actual);
    assert(!memcmp(
        expected->stars,
        actual->stars,
        (size_t)dimquads * sizeof(*expected->stars)));
    assert(!memcmp(
        expected->code,
        actual->code,
        code_count * sizeof(*expected->code)));
    assert(!memcmp(
        &expected->tol2,
        &actual->tol2,
        sizeof(expected->tol2)));
    assert(!memcmp(
        &expected->rel_field_noise2,
        &actual->rel_field_noise2,
        sizeof(expected->rel_field_noise2)));
    assert(expected->current_parity ==
        actual->current_parity);
}

static void assert_descriptor_split_equivalent(
    const solver_ab_descriptor_output_t* complete,
    const solver_ab_descriptor_output_t* lead,
    const solver_ab_descriptor_output_t* tail,
    int dimquads) {
    const solver_ab_descriptor_output_t* split[2] = {
        lead, tail
    };
    unsigned long long complete_numtries = 0ULL;
    unsigned long long complete_cxdx = 0ULL;
    unsigned long long complete_meanx = 0ULL;
    unsigned long long split_numtries = 0ULL;
    unsigned long long split_cxdx = 0ULL;
    unsigned long long split_meanx = 0ULL;
    double complete_noise = 0.0;
    double split_noise = 0.0;
    size_t split_packet = 0U;
    size_t split_descriptor = 0U;
    size_t descriptor_index;

    assert(complete);
    assert(lead);
    assert(tail);
    assert(lead->descriptor_count + tail->descriptor_count ==
        complete->descriptor_count);

    for (descriptor_index = 0U;
         descriptor_index < complete->descriptor_count;
         descriptor_index++) {
        const solver_ab_descriptor_t* expected =
            &complete->descriptors[descriptor_index];
        const solver_ab_descriptor_t* actual;

        while (split_packet < 2U &&
               split_descriptor ==
                   split[split_packet]->descriptor_count) {
            split_numtries +=
                split[split_packet]->trailing_numtries;
            split_cxdx += split[split_packet]->trailing_cxdx;
            split_meanx += split[split_packet]->trailing_meanx;
            if (split[split_packet]
                    ->has_final_rel_field_noise2) {
                split_noise = split[split_packet]
                    ->final_rel_field_noise2;
            }
            split_packet++;
            split_descriptor = 0U;
        }
        assert(split_packet < 2U);
        actual = &split[split_packet]
            ->descriptors[split_descriptor++];

        complete_numtries += expected->numtries_delta;
        complete_cxdx += expected->cxdx_delta;
        complete_meanx += expected->meanx_delta;
        split_numtries += actual->numtries_delta;
        split_cxdx += actual->cxdx_delta;
        split_meanx += actual->meanx_delta;
        complete_noise = expected->rel_field_noise2;
        split_noise = actual->rel_field_noise2;

        assert_descriptor_scientific_equal(
            expected, actual, dimquads);
        assert(complete_numtries == split_numtries);
        assert(complete_cxdx == split_cxdx);
        assert(complete_meanx == split_meanx);
        assert(!memcmp(
            &complete_noise,
            &split_noise,
            sizeof(complete_noise)));
    }

    complete_numtries += complete->trailing_numtries;
    complete_cxdx += complete->trailing_cxdx;
    complete_meanx += complete->trailing_meanx;
    if (complete->has_final_rel_field_noise2) {
        complete_noise = complete->final_rel_field_noise2;
    }
    while (split_packet < 2U) {
        assert(split_descriptor ==
            split[split_packet]->descriptor_count);
        split_numtries += split[split_packet]->trailing_numtries;
        split_cxdx += split[split_packet]->trailing_cxdx;
        split_meanx += split[split_packet]->trailing_meanx;
        if (split[split_packet]->has_final_rel_field_noise2) {
            split_noise =
                split[split_packet]->final_rel_field_noise2;
        }
        split_packet++;
        split_descriptor = 0U;
    }
    assert(complete_numtries == split_numtries);
    assert(complete_cxdx == split_cxdx);
    assert(complete_meanx == split_meanx);
    assert(complete->has_final_rel_field_noise2 ==
        tail->has_final_rel_field_noise2);
    assert(!memcmp(
        &complete_noise,
        &split_noise,
        sizeof(complete_noise)));
}

void test_solver_ab_descriptor_split_equivalence(void) {
    solver_ab_descriptor_workspace_t workspace;
    solver_ab_descriptor_planner_t* planner;
    solver_ab_descriptor_task_input_t input;
    solver_ab_descriptor_output_t* complete;
    solver_ab_descriptor_output_t* lead;
    solver_ab_descriptor_output_t* tail;
    solver_ab_pair_t* pairs = NULL;
    solver_t* solver;
    starxy_t* starxy;
    index_t index;
    size_t pair_count = 0U;
    size_t expansion;
    unsigned long long total_combinations = 0ULL;
    unsigned long long split_combination;

    memset(&workspace, 0, sizeof(workspace));
    starxy = test_solver_geometry_field();
    assert(starxy);
    solver = solver_new();
    assert(solver);
    memset(&index, 0, sizeof(index));
    index.index_scale_lower = 1;
    index.index_scale_upper = 10;
    index.dimquads = 4;
    index.cx_less_than_dx = TRUE;
    index.meanx_less_than_half = TRUE;
    solver->funits_lower = 0.1;
    solver->funits_upper = 10;
    solver->parity = PARITY_BOTH;
    solver->endobj = starxy_n(starxy);
    solver_add_index(solver, &index);
    solver_set_field(solver, starxy);
    solver_preprocess_field(solver);
    assert(solver_prepare_field_geometry(solver));

    planner = &workspace.planner;
    planner->snapshot.codetol = solver->codetol;
    planner->snapshot.rel_index_noise2 =
        solver->rel_index_noise2;
    assert(!solver_ab_collect_pairs(
        planner,
        SOLVER_AB_PHASE_DIAGONAL,
        6,
        4,
        solver->field_geometry,
        0.0,
        DBL_MAX,
        &pairs,
        &pair_count,
        &total_combinations));
    planner->pair_count = pair_count;
    expansion = solver_ab_descriptor_expansion(
        4, solver->parity);
    assert(pair_count);
    assert(total_combinations > 1ULL);
    split_combination =
        solver_ab_descriptor_lead_combinations(
            total_combinations - 1ULL, expansion);
    assert(split_combination > 0ULL);
    split_combination++;
    assert(split_combination > 1ULL);
    assert(split_combination < total_combinations);

    complete = calloc(1U, sizeof(*complete));
    lead = calloc(1U, sizeof(*lead));
    tail = calloc(1U, sizeof(*tail));
    assert(complete);
    assert(lead);
    assert(tail);
    memset(&input, 0, sizeof(input));
    input.field_geometry = solver->field_geometry;
    input.pairs = pairs;
    input.pair_count = pair_count;
    input.combination_first = 1ULL;
    input.combination_end = total_combinations;
    input.phase = SOLVER_AB_PHASE_DIAGONAL;
    input.newpoint = 6;
    input.dimquads = 4;
    input.parity = solver->parity;
    input.cx_less_than_dx = index.cx_less_than_dx;
    input.meanx_less_than_half = index.meanx_less_than_half;
    input.cxdx_margin = solver->cxdx_margin;
    assert(solver_ab_descriptor_helper_execute(
        &input, sizeof(input), complete, sizeof(*complete)) ==
        INDEX_SHARD_HELPER_TASK_OK);
    input.combination_end = split_combination;
    assert(solver_ab_descriptor_helper_execute(
        &input, sizeof(input), lead, sizeof(*lead)) ==
        INDEX_SHARD_HELPER_TASK_OK);
    input.combination_first = split_combination;
    input.combination_end = total_combinations;
    assert(solver_ab_descriptor_helper_execute(
        &input, sizeof(input), tail, sizeof(*tail)) ==
        INDEX_SHARD_HELPER_TASK_OK);
    assert_descriptor_split_equivalent(
        complete, lead, tail, 4);

    free(tail);
    free(lead);
    free(complete);
    solver_ab_descriptor_release_pairs(&workspace);
    free(workspace.planner.pair_cache);
    solver_free_field(solver);
    solver_free(solver);
}

void test_solver_index_close_fds_failure_state(void) {
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

void test_solver_zero_initialized_payload_fd_is_unowned(void) {
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
