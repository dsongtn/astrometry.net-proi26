/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "solver.h"
#include "index.h"
#include "solver_hypothesis_internal.h"
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
