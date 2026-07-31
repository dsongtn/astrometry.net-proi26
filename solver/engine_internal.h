#ifndef ASTROMETRY_ENGINE_INTERNAL_H
#define ASTROMETRY_ENGINE_INTERNAL_H

#include <stddef.h>

#include "astrometry/an-bool.h"
#include "astrometry/engine.h"
#include "astrometry/solver.h"

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

#endif
