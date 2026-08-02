/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef ASTROMETRY_ENGINE_INTERNAL_H
#define ASTROMETRY_ENGINE_INTERNAL_H

#include <stddef.h>

#include "astrometry/an-bool.h"
#include "astrometry/engine.h"
#include "astrometry/solver.h"

/*
 * Value-only engine control structures.
 *
 * A cursor produces the canonical depth/scale pass sequence. A limit policy
 * records both the effective values and their origin so terminal reporting can
 * distinguish job and configuration limits. These structures do not own
 * solver, index, pool, or payload resources.
 */

typedef struct engine_pass {
    size_t ordinal;
    size_t depth_index;
    size_t scale_index;
    int startobj;
    int endobj;
    double funits_lower;
    double funits_upper;
} engine_pass_t;

typedef struct engine_pass_cursor {
    size_t next_ordinal;
    size_t next_depth_index;
    size_t next_scale_index;
} engine_pass_cursor_t;

void engine_pass_cursor_init(engine_pass_cursor_t* cursor);

anbool engine_pass_cursor_next(const job_t* job,
                               double default_lower,
                               double default_upper,
                               engine_pass_cursor_t* cursor,
                               engine_pass_t* pass);

void engine_pass_apply(solver_t* solver, const engine_pass_t* pass);

typedef struct engine_limit_policy {
    double wall_seconds;
    double cpu_seconds;
    anbool wall_job_clamped;
    anbool wall_from_job;
    anbool cpu_from_job;
} engine_limit_policy_t;

void engine_limit_policy_resolve(double job_wall_seconds,
                                 double config_wall_seconds,
                                 double job_cpu_seconds,
                                 double config_cpu_seconds,
                                 engine_limit_policy_t* policy);

#endif
