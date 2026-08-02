/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef SOLVER_FIELD_GEOMETRY_INTERNAL_H
#define SOLVER_FIELD_GEOMETRY_INTERNAL_H

#include <stddef.h>

#include "solver.h"

/*
 * Bounded native-geometry cache contract.
 *
 * solver_t owns solver_field_geometry when field_geometry_owned is true.
 * Pair pointers borrow the table and are valid only while that owner remains
 * compatible and alive. Allocation or compatibility failure selects the
 * original per-pair calculation; it is not a scientific failure.
 */

#ifndef SOLVER_FIELD_GEOMETRY_BUDGET_BYTES
#define SOLVER_FIELD_GEOMETRY_BUDGET_BYTES \
    (64ULL * 1024ULL * 1024ULL)
#endif

typedef struct solver_pair_geometry {
    anbool scale_ok;
    double scale;
    double costheta;
    double sintheta;
    double rel_field_noise2;
} solver_pair_geometry_t;

struct solver_field_geometry {
    const starxy_t* fieldxy;
    solver_pair_geometry_t* pair_geometry;
    int numxy;
    double codetol;
    double verify_pix;
    double quadsize_min;
    double quadsize_max;
    size_t pair_capacity;
    size_t scale_ok_pairs;
    size_t bytes;
    unsigned long long reused_solver_runs;
};

void solver_release_field_geometry(solver_t* solver);

anbool solver_field_geometry_compatible(
    const solver_t* solver,
    int numxy);

const solver_pair_geometry_t* solver_field_geometry_pair(
    const solver_field_geometry_t* geometry,
    int field_a,
    int field_b);

anbool solver_pair_geometry_transform(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int star,
    double* x,
    double* y);

int solver_pair_geometry_eligible_before(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int fieldtop);

#endif
