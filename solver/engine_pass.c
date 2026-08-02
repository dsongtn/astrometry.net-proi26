/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <string.h>

#include "engine_internal.h"

void engine_pass_cursor_init(engine_pass_cursor_t* cursor) {
    if (!cursor) {
        return;
    }
    memset(cursor, 0, sizeof(*cursor));
}

anbool engine_pass_cursor_next(const job_t* job,
                               double default_lower,
                               double default_upper,
                               engine_pass_cursor_t* cursor,
                               engine_pass_t* pass) {
    size_t depth_count;
    size_t scale_count;
    int raw_start;
    int raw_end;
    double raw_lower;
    double raw_upper;

    if (!job || !job->depths || !job->scales || !cursor || !pass) {
        return FALSE;
    }
    depth_count = (size_t)il_size(job->depths) / 2U;
    scale_count = (size_t)dl_size(job->scales) / 2U;
    if (!depth_count || !scale_count ||
        cursor->next_depth_index >= depth_count) {
        return FALSE;
    }

    memset(pass, 0, sizeof(*pass));
    pass->ordinal = cursor->next_ordinal;
    pass->depth_index = cursor->next_depth_index;
    pass->scale_index = cursor->next_scale_index;

    raw_start = il_get(job->depths, pass->depth_index * 2U);
    raw_end = il_get(job->depths, pass->depth_index * 2U + 1U);
    if (raw_start < 0 || raw_end < 0) {
        return FALSE;
    }
    pass->startobj = raw_start ? raw_start - 1 : 0;
    /*
     * The user-facing upper bound is inclusive and one-based. Its numeric
     * value is therefore already the zero-based exclusive bound. Zero is the
     * native open-upper sentinel and must be written on every pass.
     */
    pass->endobj = raw_end;

    raw_lower = dl_get(job->scales, pass->scale_index * 2U);
    raw_upper = dl_get(job->scales, pass->scale_index * 2U + 1U);
    pass->funits_lower =
        raw_lower == 0.0 ? default_lower : raw_lower;
    pass->funits_upper =
        raw_upper == 0.0 ? default_upper : raw_upper;

    cursor->next_scale_index++;
    cursor->next_ordinal++;
    if (cursor->next_scale_index >= scale_count) {
        cursor->next_scale_index = 0U;
        cursor->next_depth_index++;
    }
    return TRUE;
}

void engine_pass_apply(solver_t* solver, const engine_pass_t* pass) {
    if (!solver || !pass) {
        return;
    }
    solver->startobj = pass->startobj;
    solver->endobj = pass->endobj;
    solver->funits_lower = pass->funits_lower;
    solver->funits_upper = pass->funits_upper;
}
