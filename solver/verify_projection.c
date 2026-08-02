/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sip-utils.h"
#include "verify_internal.h"

#define VERIFY_PROJECTION_MAX_TASKS 6U
#define VERIFY_PROJECTION_MIN_STARS_PER_TASK 256U

typedef struct verify_projection_context {
    const double* xyz;
    int npoints;
    int width;
    int height;
    anbool use_sip;
    union {
        sip_t sip;
        tan_t tan;
    } wcs;
} verify_projection_context_t;

typedef struct verify_projection_task_input {
    const verify_projection_context_t* context;
    size_t first;
    size_t count;
} verify_projection_task_input_t;

typedef struct verify_projection_task_output {
    double x;
    double y;
    anbool inside;
} verify_projection_task_output_t;

static index_shard_helper_task_status_t
verify_projection_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const verify_projection_task_input_t* input = input_bytes;
    const verify_projection_context_t* context;
    verify_projection_task_output_t* output = output_bytes;
    size_t point;

    if (!input || input_size != sizeof(*input) ||
        !output ||
        input->count > SIZE_MAX / sizeof(*output) ||
        output_size != input->count * sizeof(*output)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    context = input->context;
    if (!context || !context->xyz || context->npoints < 0 ||
        input->first > (size_t)context->npoints ||
        input->count >
            (size_t)context->npoints - input->first) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    for (point = 0U; point < input->count; point++) {
        const double* xyz = context->xyz +
            3U * (input->first + point);
        double x = 0.0;
        double y = 0.0;
        anbool projected;

        if (context->use_sip) {
            projected = sip_xyzarr2pixelxy(
                &context->wcs.sip, xyz, &x, &y);
        } else {
            projected = tan_xyzarr2pixelxy(
                &context->wcs.tan, xyz, &x, &y);
        }
        output[point].x = x;
        output[point].y = y;
        output[point].inside =
            projected &&
            x >= 0.0 && y >= 0.0 &&
            x < context->width && y < context->height;
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

static const index_shard_helper_ops_t
verify_projection_helper_ops = {
    "verify-projection",
    verify_projection_helper_execute
};

/*
 * Project only the already-copied startree_search_for() result. The helper
 * package is index-free and immutable; stable owner compaction preserves the
 * exact input ordering of sip_filter_stars_in_field().
 */
anbool verify_internal_filter_stars_in_field_parallel(
    const sip_t* sip,
    const tan_t* tan,
    const double* xyz,
    int npoints,
    double** p_xy,
    int** p_inbounds,
    int* p_ngood) {
#if VERIFY_INTERNAL_HELPERS_LINKED
    verify_projection_context_t context;
    verify_projection_task_input_t
        inputs[VERIFY_PROJECTION_MAX_TASKS];
    index_shard_helper_task_t
        tasks[VERIFY_PROJECTION_MAX_TASKS];
    verify_projection_task_output_t* projected = NULL;
    int* inbounds = NULL;
    double* xy = NULL;
    size_t available;
    size_t task_count;
    size_t base;
    size_t remainder;
    size_t cursor = 0U;
    size_t task_index;
    int ngood = 0;
    int point;
    index_shard_helper_run_status_t run_status;

    if (!index_shard_helper_available_workers ||
        !index_shard_helper_run ||
        (!sip && !tan) || !xyz || npoints <= 0 ||
        !p_inbounds || !p_ngood ||
        (size_t)npoints <
            2U * VERIFY_PROJECTION_MIN_STARS_PER_TASK) {
        return FALSE;
    }
    available = index_shard_helper_available_workers();
    if (!available) {
        return FALSE;
    }
    task_count = (size_t)npoints /
        VERIFY_PROJECTION_MIN_STARS_PER_TASK;
    if (task_count > VERIFY_PROJECTION_MAX_TASKS) {
        task_count = VERIFY_PROJECTION_MAX_TASKS;
    }
    if (available < task_count - 1U) {
        task_count = available + 1U;
    }
    if (task_count < 2U ||
        (size_t)npoints > SIZE_MAX / sizeof(*projected) ||
        (size_t)npoints > SIZE_MAX / sizeof(*inbounds) ||
        (p_xy &&
         (size_t)npoints > SIZE_MAX / (2U * sizeof(*xy)))) {
        return FALSE;
    }

    projected = malloc((size_t)npoints * sizeof(*projected));
    inbounds = malloc((size_t)npoints * sizeof(*inbounds));
    if (p_xy) {
        xy = malloc((size_t)npoints * 2U * sizeof(*xy));
    }
    if (!projected || !inbounds || (p_xy && !xy)) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    memset(&context, 0, sizeof(context));
    context.xyz = xyz;
    context.npoints = npoints;
    context.use_sip = sip != NULL;
    if (sip) {
        context.wcs.sip = *sip;
        context.width = sip->wcstan.imagew;
        context.height = sip->wcstan.imageh;
    } else {
        context.wcs.tan = *tan;
        context.width = tan->imagew;
        context.height = tan->imageh;
    }

    memset(inputs, 0, sizeof(inputs));
    memset(tasks, 0, sizeof(tasks));
    base = (size_t)npoints / task_count;
    remainder = (size_t)npoints % task_count;
    for (task_index = 0U;
         task_index < task_count;
         task_index++) {
        size_t count = base +
            (task_index < remainder ? 1U : 0U);

        inputs[task_index].context = &context;
        inputs[task_index].first = cursor;
        inputs[task_index].count = count;
        tasks[task_index].input = &inputs[task_index];
        tasks[task_index].input_bytes = sizeof(inputs[task_index]);
        tasks[task_index].output = projected + cursor;
        tasks[task_index].output_bytes =
            count * sizeof(*projected);
        tasks[task_index].work_units =
            (unsigned long long)count;
        cursor += count;
    }
    if (cursor != (size_t)npoints) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    run_status = index_shard_helper_run(
        &verify_projection_helper_ops,
        tasks,
        task_count,
        NULL);
    if (run_status != INDEX_SHARD_HELPER_OK) {
        free(xy);
        free(inbounds);
        free(projected);
        return FALSE;
    }

    for (point = 0; point < npoints; point++) {
        if (!projected[point].inside) {
            continue;
        }
        inbounds[ngood] = point;
        if (xy) {
            xy[2 * ngood] = projected[point].x;
            xy[2 * ngood + 1] = projected[point].y;
        }
        ngood++;
    }
    free(projected);

    if (!ngood) {
        free(xy);
        free(inbounds);
        xy = NULL;
        inbounds = NULL;
    }
    if (p_xy) {
        *p_xy = xy;
    }
    *p_inbounds = inbounds;
    *p_ngood = ngood;
    return TRUE;
#else
    (void)sip;
    (void)tan;
    (void)xyz;
    (void)npoints;
    (void)p_xy;
    (void)p_inbounds;
    (void)p_ngood;
    return FALSE;
#endif
}
