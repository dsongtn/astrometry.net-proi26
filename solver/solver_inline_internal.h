/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef SOLVER_INLINE_INTERNAL_H
#define SOLVER_INLINE_INTERNAL_H

#include <math.h>
#include <stdint.h>

#include "solver.h"
#include "tic.h"
#include "index_shard_internal.h"

enum {
    A = 0,
    B = 1,
    C = 2,
    D = 3,
    NBACK = 2
};

static inline double getx(const double* values, int index) {
    return values[index * 2];
}

static inline double gety(const double* values, int index) {
    return values[index * 2 + 1];
}

static inline void setx(double* values, int index, double value) {
    values[index * 2] = value;
}

static inline void sety(double* values, int index, double value) {
    values[index * 2 + 1] = value;
}

static inline anbool solver_poll_worker_stop(solver_t* solver) {
    if (!index_shard_worker_stop_requested()) {
        return FALSE;
    }

    solver->quit_now = TRUE;
    return TRUE;
}

static inline double field_getx(const solver_t* solver, int index) {
    return starxy_getx(solver->fieldxy, index);
}

static inline double field_gety(const solver_t* solver, int index) {
    return starxy_gety(solver->fieldxy, index);
}

static inline void update_timeused(solver_t* solver) {
    double usertime;
    double systime;

    get_resource_stats(&usertime, &systime, NULL);
    solver->timeused = (usertime + systime) - solver->starttime;
    if (solver->timeused < 0.0) {
        solver->timeused = 0.0;
    }
}

/*
 * Keep the native pquad path and the compact triangular-geometry path on one
 * explicitly contracted arithmetic sequence. With -march=native -O3 GCC
 * otherwise vectorizes check_inbox() but scalarizes the compact caller, and
 * is free to choose the opposite multiplication as the fused addend. Both
 * expressions are mathematically equivalent, but their last bits can differ
 * and then perturb the scientific traversal hash between W1 and W>1.
 */
static inline void solver_transform_code_coordinates(
    double relative_x,
    double relative_y,
    double costheta,
    double sintheta,
    double* code_x,
    double* code_y) {
    double x_addend = relative_y * sintheta;
    double y_addend = relative_y * costheta;

    *code_x = fma(relative_x, costheta, x_addend);
    *code_y = fma(-relative_x, sintheta, y_addend);
}

static inline uint64_t solver_order_hash_mix(
    uint64_t state,
    uint64_t value) {
    int byte_index;

    if (!state) {
        state = UINT64_C(1469598103934665603);
    }
    for (byte_index = 0; byte_index < 8; byte_index++) {
        state ^= (value >> (8 * byte_index)) & UINT64_C(0xff);
        state *= UINT64_C(1099511628211);
    }
    return state;
}

#endif
