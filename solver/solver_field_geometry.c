/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* Bounded field-geometry cache construction and pair lookup. */

#include <stdint.h>
#include <stdlib.h>

#include "solver.h"
#include "log.h"
#include "mathutil.h"
#include "tic.h"
#include "solver_field_geometry_internal.h"
#include "solver_inline_internal.h"

static void solver_free_field_geometry_object(
    solver_field_geometry_t* geometry) {
    if (!geometry) {
        return;
    }
    free(geometry->pair_geometry);
    free(geometry);
}

void solver_release_field_geometry(solver_t* solver) {
    solver_field_geometry_t* geometry;

    if (!solver || !solver->field_geometry) {
        return;
    }
    if (!solver->field_geometry_owned) {
        solver->field_geometry = NULL;
        return;
    }

    geometry = solver->field_geometry;
    logverb("[solver-geometry] stats mode=compact-triangular "
            "objects=%i pairs=%zu bytes=%zu reused_solver_runs=%llu\n",
            geometry->numxy,
            geometry->scale_ok_pairs,
            geometry->bytes,
            __atomic_load_n(
                &geometry->reused_solver_runs,
                __ATOMIC_RELAXED));
    solver_free_field_geometry_object(geometry);
    solver->field_geometry = NULL;
    solver->field_geometry_owned = FALSE;
}

static int solver_field_geometry_numxy(const solver_t* solver) {
    int numxy;

    if (!solver || !solver->fieldxy) {
        return 0;
    }
    numxy = starxy_n(solver->fieldxy);
    if (solver->endobj && numxy > solver->endobj) {
        numxy = solver->endobj;
    }
    if (numxy >= 1000) {
        numxy = 1000;
    }
    return numxy;
}

static anbool solver_field_geometry_estimate(
    int numxy,
    size_t* pairs,
    size_t* bytes) {
    size_t n;

    if (numxy <= 0 || !pairs || !bytes) {
        return FALSE;
    }
    n = (size_t)numxy;
    if (n > 1U && n > SIZE_MAX / (n - 1U)) {
        return FALSE;
    }
    *pairs = n * (n - 1U) / 2U;
    if (*pairs >
        (SIZE_MAX - sizeof(solver_field_geometry_t)) /
            sizeof(solver_pair_geometry_t)) {
        return FALSE;
    }
    *bytes = sizeof(solver_field_geometry_t) +
        *pairs * sizeof(solver_pair_geometry_t);
    return TRUE;
}

static void solver_field_geometry_check_scale(
    solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    const solver_t* solver) {
    double dx;
    double dy;
    double minimum_scale;
    double maximum_scale;

    dx = field_getx(solver, field_b) -
        field_getx(solver, field_a);
    dy = field_gety(solver, field_b) -
        field_gety(solver, field_a);
    pair->scale = dx * dx + dy * dy;
    minimum_scale = square(solver->quadsize_min);
    maximum_scale = solver->quadsize_max > 0.0
        ? square(solver->quadsize_max)
        : LARGE_VAL;
    if (pair->scale < minimum_scale ||
        pair->scale > maximum_scale) {
        pair->scale_ok = FALSE;
        return;
    }

    pair->costheta = (dy + dx) / pair->scale;
    pair->sintheta = (dy - dx) / pair->scale;
    pair->rel_field_noise2 =
        (solver->verify_pix * solver->verify_pix) / pair->scale;
    pair->scale_ok = TRUE;
}

static size_t solver_field_geometry_pair_index(
    int field_a,
    int field_b) {
    return (size_t)field_b *
        (size_t)(field_b - 1) / 2U +
        (size_t)field_a;
}

const solver_pair_geometry_t*
solver_field_geometry_pair(
    const solver_field_geometry_t* geometry,
    int field_a,
    int field_b) {
    size_t pair_index;

    if (!geometry || !geometry->pair_geometry ||
        field_a < 0 || field_b <= field_a ||
        field_b >= geometry->numxy) {
        return NULL;
    }
    pair_index = solver_field_geometry_pair_index(
        field_a,
        field_b);
    if (pair_index >= geometry->pair_capacity) {
        return NULL;
    }
    return geometry->pair_geometry + pair_index;
}

anbool solver_pair_geometry_transform(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int star,
    double* x,
    double* y) {
    double ax;
    double ay;
    double cx;
    double cy;
    double radius_term;

    if (!geometry || !pair || !pair->scale_ok ||
        !geometry->fieldxy || !x || !y ||
        star < 0 || star >= geometry->numxy ||
        star == field_a || star == field_b) {
        return FALSE;
    }
    ax = starxy_getx(geometry->fieldxy, field_a);
    ay = starxy_gety(geometry->fieldxy, field_a);
    cx = starxy_getx(geometry->fieldxy, star) - ax;
    cy = starxy_gety(geometry->fieldxy, star) - ay;
    solver_transform_code_coordinates(
        cx,
        cy,
        pair->costheta,
        pair->sintheta,
        &cx,
        &cy);
    radius_term =
        (cx * cx - cx) + (cy * cy - cy);
    if (radius_term >
        geometry->codetol *
            (M_SQRT2 + geometry->codetol)) {
        return FALSE;
    }
    *x = cx;
    *y = cy;
    return TRUE;
}

int solver_pair_geometry_eligible_before(
    const solver_field_geometry_t* geometry,
    const solver_pair_geometry_t* pair,
    int field_a,
    int field_b,
    int fieldtop) {
    int eligible = 0;
    int star;

    for (star = 0; star < fieldtop; star++) {
        double x;
        double y;

        if (solver_pair_geometry_transform(
                geometry,
                pair,
                field_a,
                field_b,
                star,
                &x,
                &y)) {
            eligible++;
        }
    }
    return eligible;
}

anbool solver_field_geometry_compatible(
    const solver_t* solver,
    int numxy) {
    const solver_field_geometry_t* geometry;

    if (!solver || !solver->field_geometry) {
        return FALSE;
    }
    geometry = solver->field_geometry;
    return geometry->fieldxy == solver->fieldxy &&
        geometry->numxy >= numxy &&
        geometry->codetol == solver->codetol &&
        geometry->verify_pix == solver->verify_pix &&
        geometry->quadsize_min == solver->quadsize_min &&
        geometry->quadsize_max == solver->quadsize_max;
}

void solver_release_incompatible_field_geometry(
    solver_t* solver) {
    int numxy;

    if (!solver || !solver->field_geometry ||
        !solver->fieldxy) {
        return;
    }
    numxy = solver_field_geometry_numxy(solver);
    if (!solver_field_geometry_compatible(
            solver, numxy)) {
        solver_release_field_geometry(solver);
    }
}

anbool solver_prepare_field_geometry(solver_t* solver) {
    solver_field_geometry_t* geometry;
    size_t pairs;
    size_t estimated_bytes;
    size_t built_pairs = 0;
    double wall_start;
    int numxy;
    int fieldA;
    int fieldB;

    if (!solver || !solver->fieldxy) {
        return FALSE;
    }
    wall_start = monotonic_seconds();

    numxy = starxy_n(solver->fieldxy);
    if (numxy >= 1000) {
        numxy = 1000;
    }
    if (solver_field_geometry_compatible(
            solver,
            solver_field_geometry_numxy(solver)) &&
        solver->field_geometry->numxy == numxy) {
        return TRUE;
    }
    solver_release_field_geometry(solver);

    if (solver->startobj >=
        solver_field_geometry_numxy(solver)) {
        logverb("[solver-geometry] mode=native reason=empty-range "
                "startobj=%i objects=%i\n",
                solver->startobj,
                numxy);
        return FALSE;
    }
    if (!solver_field_geometry_estimate(
            numxy,
            &pairs,
            &estimated_bytes)) {
        logverb("[solver-geometry] mode=native reason=estimate "
                "objects=%i\n",
                numxy);
        return FALSE;
    }
    if (estimated_bytes > SOLVER_FIELD_GEOMETRY_BUDGET_BYTES) {
        logverb("[solver-geometry] mode=native reason=budget "
                "objects=%i estimate=%zu budget=%llu\n",
                numxy,
                estimated_bytes,
                (unsigned long long)
                    SOLVER_FIELD_GEOMETRY_BUDGET_BYTES);
        return FALSE;
    }

    geometry = calloc(1, sizeof(*geometry));
    if (!geometry) {
        logverb("[solver-geometry] mode=native reason=allocation "
                "objects=%i estimate=%zu\n",
                numxy,
                estimated_bytes);
        return FALSE;
    }
    geometry->pair_capacity = pairs;
    if (pairs) {
        geometry->pair_geometry = calloc(
            pairs,
            sizeof(*geometry->pair_geometry));
    }
    if (pairs && !geometry->pair_geometry) {
        solver_free_field_geometry_object(geometry);
        logverb("[solver-geometry] mode=native reason=allocation "
                "objects=%i estimate=%zu\n",
                numxy,
                estimated_bytes);
        return FALSE;
    }

    for (fieldB = 1; fieldB < numxy; fieldB++) {
        for (fieldA = 0; fieldA < fieldB; fieldA++) {
            solver_pair_geometry_t* pair =
                geometry->pair_geometry +
                solver_field_geometry_pair_index(
                    fieldA,
                    fieldB);

            solver_field_geometry_check_scale(
                pair,
                fieldA,
                fieldB,
                solver);
            if (pair->scale_ok) {
                built_pairs++;
            }
        }
    }

    geometry->fieldxy = solver->fieldxy;
    geometry->numxy = numxy;
    geometry->codetol = solver->codetol;
    geometry->verify_pix = solver->verify_pix;
    geometry->quadsize_min = solver->quadsize_min;
    geometry->quadsize_max = solver->quadsize_max;
    geometry->scale_ok_pairs = built_pairs;
    geometry->bytes = estimated_bytes;

    solver->field_geometry = geometry;
    solver->field_geometry_owned = TRUE;

    logverb("[solver-geometry] mode=compact-triangular objects=%i "
            "pairs=%zu possible_pairs=%zu bytes=%zu "
            "per_star_payload=none "
            "prepare_wall=%.6f\n",
            numxy,
            built_pairs,
            pairs,
            estimated_bytes,
            monotonic_seconds() - wall_start);
    return TRUE;
}
