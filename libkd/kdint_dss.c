/*
# This file is part of libkd.
# Licensed under a 3-clause BSD style license - see LICENSE
*/

#include <errno.h>
#include <limits.h>

#include "kdtree.h"
#include "kdtree_internal_common.h"

#include "kdint_etype_d.h"
#include "kdint_dtype_s.h"
#include "kdint_ttype_s.h"

#define POINT_ED(kd, d, r, func)    (func(POINT_SCALE(kd, d, r)))
#define POINT_DT(kd, d, r, func)    (r)
#define POINT_ET(kd, d, r, func)    (func(POINT_SCALE(kd, d, r)))

#define POINT_TD(kd, d, r)          (r)
#define POINT_DE(kd, d, r)          POINT_INVSCALE(kd, d, r)
#define POINT_TE(kd, d, r)          POINT_INVSCALE(kd, d, r)

#define DIST_ED(kd, dist, func)     (func(DIST_SCALE(kd, dist)))
#define DIST_DT(kd, dist, func)     (dist)
#define DIST_ET(kd, dist, func)     (func(DIST_SCALE(kd, dist)))

#define DIST_TD(kd, dist)           (dist)
#define DIST_DE(kd, dist)           DIST_INVSCALE(kd, dist)
#define DIST_TE(kd, dist)           DIST_INVSCALE(kd, dist)

#define DIST2_ED(kd, dist2, func)   (func(DIST2_SCALE(kd, dist2)))
#define DIST2_DT(kd, dist2, func)   (dist2)
#define DIST2_ET(kd, dist2, func)   (func(DIST2_SCALE(kd, dist2)))

#define DIST2_TD(kd, dist2)         (dist2)
#define DIST2_DE(kd, dist2)         DIST2_INVSCALE(kd, dist2)
#define DIST2_TE(kd, dist2)         DIST2_INVSCALE(kd, dist2)

#define EQUAL_ED 0
#define EQUAL_DT 1
#define EQUAL_ET 0

/*
 * Enable the specialized 4D CodeKD leaf kernel only for the dominant
 * double/U16/U16 specialization.
 */
#define KDTREE_CODEKD_DSS_U16_FAST_PATH 1

#include "kdtree_internal.c"

#include "kdtree_direct_internal.h"

#define KDTREE_DIRECT_SPAN_CAPACITY KDTREE_DIRECT_DSS_TASK_SPANS

typedef struct kdtree_direct_span {
    int left;
    int right;
} kdtree_direct_span_t;

typedef int (*kdtree_direct_flush_fn)(
    void* opaque,
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last,
    anbool final);

int kdtree_direct_dss_task_execute(
    const kdtree_direct_dss_task_input_t* input,
    kdtree_direct_dss_task_output_t* output) {
    size_t span_index;

    if (!output) {
        return -1;
    }
    output->nres = 0U;
    output->status = -1;
    if (!input || !input->data || !input->nspans ||
        input->nspans > KDTREE_DIRECT_DSS_TASK_SPANS ||
        !input->cover_points ||
        input->cover_points > KDTREE_DIRECT_DSS_TASK_POINTS) {
        return -1;
    }
    for (span_index = 0U;
         span_index < input->nspans;
         span_index++) {
        u16 first = input->span_first[span_index];
        u16 last = input->span_last[span_index];
        unsigned int offset;

        if (first > last || (u32)last >= input->cover_points) {
            return -1;
        }
        for (offset = first; offset <= (unsigned int)last; offset++) {
            const u16* point = input->data + (size_t)offset * 4U;
            double pp;
            double delta;
            double d2 = 0.0;

            pp = ((double)point[0] * input->invscale) +
                input->minval[0];
            delta = input->query[0] - pp;
            d2 += delta * delta;
            if (d2 > input->maxd2) {
                continue;
            }
            pp = ((double)point[1] * input->invscale) +
                input->minval[1];
            delta = input->query[1] - pp;
            d2 += delta * delta;
            if (d2 > input->maxd2) {
                continue;
            }
            pp = ((double)point[2] * input->invscale) +
                input->minval[2];
            delta = input->query[2] - pp;
            d2 += delta * delta;
            if (d2 > input->maxd2) {
                continue;
            }
            pp = ((double)point[3] * input->invscale) +
                input->minval[3];
            delta = input->query[3] - pp;
            d2 += delta * delta;
            if (d2 > input->maxd2) {
                continue;
            }
            if (output->nres >= KDTREE_DIRECT_DSS_TASK_POINTS) {
                return -1;
            }
            output->offsets[output->nres] = (u16)offset;
            output->sdists[output->nres] = d2;
            output->nres++;
        }
    }
    output->status = 0;
    return 0;
}

static int kdtree_direct_flush_spans(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last,
    kdtree_direct_read_range_fn read_range,
    void* read_opaque) {
    const u16* loaded_data = NULL;
    const u32* loaded_perm = NULL;
    size_t span_index;
    int D = kd->ndim;

    if (!nspans) {
        return 0;
    }
    if (first < 0 || last < first || last >= kd->ndata) {
        errno = EINVAL;
        return -1;
    }
    if (read_range(read_opaque,
                   first,
                   last - first + 1,
                   &loaded_data,
                   &loaded_perm)) {
        if (!errno) {
            errno = EIO;
        }
        return -1;
    }
    if (!loaded_data || (kd->perm && !loaded_perm)) {
        errno = EIO;
        return -1;
    }

    for (span_index = 0; span_index < nspans; span_index++) {
        int i;
        const kdtree_direct_span_t* span = spans + span_index;

        for (i = span->left; i <= span->right; i++) {
            anbool bailedout = FALSE;
            const dtype* point = loaded_data +
                (size_t)(i - first) * (size_t)D;
            double dsqd = 0.0;
            unsigned int index = loaded_perm
                ? loaded_perm[i - first]
                : (unsigned int)i;

            dist2_bailout(kd,
                          query,
                          point,
                          D,
                          maxd2,
                          &bailedout,
                          &dsqd);
            if (bailedout) {
                continue;
            }
            if (!add_result(kd,
                            result,
                            dsqd,
                            index,
                            point,
                            D,
                            TRUE,
                            FALSE)) {
                errno = ENOMEM;
                return -1;
            }
        }
    }
    return 0;
}

typedef struct kdtree_direct_inline_reader {
    kdtree_direct_read_range_fn read_range;
    void* read_opaque;
} kdtree_direct_inline_reader_t;

static int kdtree_direct_flush_inline(
    void* opaque,
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last,
    anbool final) {
    kdtree_direct_inline_reader_t* reader = opaque;

    (void)final;
    return kdtree_direct_flush_spans(
        kd,
        result,
        query,
        maxd2,
        spans,
        nspans,
        first,
        last,
        reader->read_range,
        reader->read_opaque);
}

static int kdtree_direct_emit_span(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    int left,
    int right,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_flush_fn flush,
    void* flush_opaque,
    kdtree_direct_span_t* spans,
    size_t* nspans,
    int* batch_first,
    int* batch_last) {
    size_t span_points;

    if (left > right) {
        return 0;
    }
    span_points = (size_t)(right - left) + 1U;

    if (span_points > max_points) {
        int first = left;

        if (flush(flush_opaque,
                  kd,
                  result,
                  query,
                  maxd2,
                  spans,
                  *nspans,
                  *batch_first,
                  *batch_last,
                  FALSE)) {
            return -1;
        }
        *nspans = 0;
        *batch_first = -1;
        *batch_last = -1;

        while (first <= right) {
            int count = right - first + 1;
            kdtree_direct_span_t large_span;

            if ((size_t)count > max_points) {
                count = (int)max_points;
            }
            large_span.left = first;
            large_span.right = first + count - 1;
            if (flush(flush_opaque,
                      kd,
                      result,
                      query,
                      maxd2,
                      &large_span,
                      1U,
                      large_span.left,
                      large_span.right,
                      FALSE)) {
                return -1;
            }
            first += count;
        }
        return 0;
    }

    if (*nspans) {
        int new_first = left < *batch_first ? left : *batch_first;
        int new_last = right > *batch_last ? right : *batch_last;
        size_t union_points = (size_t)(new_last - new_first) + 1U;
        size_t gap_points = 0;

        if (right < *batch_first) {
            gap_points = (size_t)(*batch_first - right - 1);
        } else if (left > *batch_last) {
            gap_points = (size_t)(left - *batch_last - 1);
        }

        if (*nspans == KDTREE_DIRECT_SPAN_CAPACITY ||
            union_points > max_points ||
            gap_points > merge_gap_points) {
            if (flush(flush_opaque,
                      kd,
                      result,
                      query,
                      maxd2,
                      spans,
                      *nspans,
                      *batch_first,
                      *batch_last,
                      FALSE)) {
                return -1;
            }
            *nspans = 0;
            *batch_first = -1;
            *batch_last = -1;
        }
    }

    spans[*nspans].left = left;
    spans[*nspans].right = right;
    (*nspans)++;
    if (*batch_first < 0 || left < *batch_first) {
        *batch_first = left;
    }
    if (*batch_last < 0 || right > *batch_last) {
        *batch_last = right;
    }
    return 0;
}

typedef struct kdtree_direct_leased_group {
    kdtree_direct_span_t spans[KDTREE_DIRECT_SPAN_CAPACITY];
    size_t nspans;
    size_t selected_points;
    int first;
    int last;
} kdtree_direct_leased_group_t;

typedef struct kdtree_direct_leased_state {
    const kdtree_t* kd;
    kdtree_qres_t* result;
    const etype* query;
    double maxd2;
    kdtree_direct_read_leased_ranges_fn read_ranges;
    kdtree_direct_release_range_fn release_range;
    void* read_opaque;
    const kdtree_direct_dss_executor_t* executor;
    const kdtree_direct_dss_lookahead_t* lookahead;
    kdtree_direct_leased_group_t groups[KDTREE_DIRECT_DSS_WAVE_TASKS];
    size_t ngroups;
    kdtree_direct_leased_group_t
        pending_groups[KDTREE_DIRECT_DSS_WAVE_TASKS];
    size_t npending_groups;
    void* pending_handle;
} kdtree_direct_leased_state_t;

static void kdtree_direct_release_leased_range(
    kdtree_direct_leased_state_t* state,
    kdtree_direct_range_t* range) {
    if (range->lease) {
        state->release_range(state->read_opaque, range->lease);
    }
    memset(range, 0, sizeof(*range));
}

static int kdtree_direct_build_group_requests(
    const kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* groups,
    size_t ngroups,
    kdtree_direct_range_request_t* requests) {
    size_t group_index;

    if (!state || !state->kd || !groups || !requests ||
        !ngroups || ngroups > KDTREE_DIRECT_DSS_WAVE_TASKS) {
        errno = EINVAL;
        return -1;
    }
    for (group_index = 0U; group_index < ngroups; group_index++) {
        const kdtree_direct_leased_group_t* group =
            &groups[group_index];
        int count = group->last - group->first + 1;

        if (group->first < 0 || group->last < group->first ||
            group->last >= state->kd->ndata || count <= 0) {
            errno = EINVAL;
            return -1;
        }
        requests[group_index].first = group->first;
        requests[group_index].count = count;
    }
    return 0;
}

static int kdtree_direct_read_leased_groups(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* groups,
    size_t ngroups,
    kdtree_direct_range_t* ranges) {
    kdtree_direct_range_request_t
        requests[KDTREE_DIRECT_DSS_WAVE_TASKS];
    size_t group_index;
    int saved_errno;

    memset(ranges,
           0,
           KDTREE_DIRECT_DSS_WAVE_TASKS * sizeof(*ranges));
    if (!groups || !ngroups ||
        ngroups > KDTREE_DIRECT_DSS_WAVE_TASKS) {
        errno = EINVAL;
        return -1;
    }
    if (kdtree_direct_build_group_requests(
            state, groups, ngroups, requests)) {
        return -1;
    }
    if (state->read_ranges(state->read_opaque,
                           requests,
                           ngroups,
                           ranges)) {
        saved_errno = errno ? errno : EIO;
        goto fail;
    }
    for (group_index = 0U;
         group_index < ngroups;
         group_index++) {
        kdtree_direct_range_t* range = &ranges[group_index];

        if (!range->lease) {
            saved_errno = EPROTO;
            goto fail;
        }
        if (!range->data || (state->kd->perm && !range->perm)) {
            saved_errno = EIO;
            goto fail;
        }
    }
    return 0;

fail:
    for (group_index = 0U;
         group_index < ngroups;
         group_index++) {
        kdtree_direct_release_leased_range(
            state, &ranges[group_index]);
    }
    errno = saved_errno;
    return -1;
}

static int kdtree_direct_reduce_inline_group(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* group,
    const kdtree_direct_range_t* range) {
    size_t span_index;
    int D = state->kd->ndim;

    for (span_index = 0U;
         span_index < group->nspans;
         span_index++) {
        const kdtree_direct_span_t* span = &group->spans[span_index];
        int i;

        for (i = span->left; i <= span->right; i++) {
            anbool bailedout = FALSE;
            size_t offset = (size_t)(i - group->first);
            const dtype* point = range->data + offset * (size_t)D;
            double dsqd = 0.0;
            unsigned int index = range->perm
                ? range->perm[offset]
                : (unsigned int)i;

            dist2_bailout(state->kd,
                          state->query,
                          point,
                          D,
                          state->maxd2,
                          &bailedout,
                          &dsqd);
            if (bailedout) {
                continue;
            }
            if (!add_result(state->kd,
                            state->result,
                            dsqd,
                            index,
                            point,
                            D,
                            TRUE,
                            FALSE)) {
                errno = ENOMEM;
                return -1;
            }
        }
    }
    return 0;
}

static int kdtree_direct_prepare_task_input(
    const kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* group,
    const kdtree_direct_range_t* range,
    kdtree_direct_dss_task_input_t* input) {
    size_t cover_points = (size_t)(group->last - group->first) + 1U;
    size_t span_index;
    int d;

    if (state->kd->ndim != 4 ||
        cover_points > KDTREE_DIRECT_DSS_TASK_POINTS ||
        group->nspans > KDTREE_DIRECT_DSS_TASK_SPANS) {
        return -1;
    }
    memset(input, 0, sizeof(*input));
    input->data = range->data;
    input->invscale = state->kd->invscale;
    input->maxd2 = state->maxd2;
    input->cover_points = (u32)cover_points;
    input->nspans = (u16)group->nspans;
    for (d = 0; d < 4; d++) {
        input->query[d] = state->query[d];
        input->minval[d] = state->kd->minval[d];
    }
    for (span_index = 0U;
         span_index < group->nspans;
         span_index++) {
        int first = group->spans[span_index].left - group->first;
        int last = group->spans[span_index].right - group->first;

        if (first < 0 || last < first ||
            (size_t)last >= cover_points) {
            return -1;
        }
        input->span_first[span_index] = (u16)first;
        input->span_last[span_index] = (u16)last;
    }
    return 0;
}

static int kdtree_direct_reduce_task_output(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* group,
    const kdtree_direct_range_t* range,
    const kdtree_direct_dss_task_input_t* input,
    const kdtree_direct_dss_task_output_t* output) {
    size_t span_cursor = 0U;
    u16 previous_offset = 0U;
    anbool have_previous = FALSE;
    size_t result_index;

    if (!output || output->status ||
        output->nres > KDTREE_DIRECT_DSS_TASK_POINTS ||
        output->nres > group->selected_points) {
        errno = EPROTO;
        return -1;
    }
    for (result_index = 0U;
         result_index < output->nres;
         result_index++) {
        u16 offset = output->offsets[result_index];
        const dtype* point;
        unsigned int index;

        if ((u32)offset >= input->cover_points) {
            errno = EPROTO;
            return -1;
        }
        while (span_cursor < input->nspans) {
            u16 first = input->span_first[span_cursor];
            u16 last = input->span_last[span_cursor];

            if (offset >= first && offset <= last) {
                break;
            }
            span_cursor++;
            have_previous = FALSE;
        }
        if (span_cursor == input->nspans ||
            (have_previous && offset <= previous_offset)) {
            errno = EPROTO;
            return -1;
        }
        previous_offset = offset;
        have_previous = TRUE;
        point = range->data + (size_t)offset * 4U;
        index = range->perm
            ? range->perm[offset]
            : (unsigned int)(group->first + (int)offset);
        if (!add_result(state->kd,
                        state->result,
                        output->sdists[result_index],
                        index,
                        point,
                        4,
                        TRUE,
                        FALSE)) {
            errno = ENOMEM;
            return -1;
        }
    }
    return 0;
}

static int kdtree_direct_process_leased_wave(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* groups,
    size_t ngroups,
    kdtree_direct_range_t* prepared_ranges) {
    kdtree_direct_range_t ranges[KDTREE_DIRECT_DSS_WAVE_TASKS];
    kdtree_direct_dss_task_input_t inputs[KDTREE_DIRECT_DSS_WAVE_TASKS];
    const kdtree_direct_dss_task_output_t* outputs = NULL;
    size_t selected_points = 0U;
    size_t loaded = 0U;
    size_t group_index;
    anbool task_eligible = TRUE;
    anbool use_executor;
    int status = 0;

    if (!ngroups) {
        return 0;
    }
    if (!groups || ngroups > KDTREE_DIRECT_DSS_WAVE_TASKS) {
        errno = EINVAL;
        return -1;
    }
    memset(ranges, 0, sizeof(ranges));
    if (prepared_ranges) {
        memcpy(ranges,
               prepared_ranges,
               ngroups * sizeof(ranges[0]));
        memset(prepared_ranges,
               0,
               ngroups * sizeof(prepared_ranges[0]));
    }
    for (group_index = 0U;
         group_index < ngroups;
         group_index++) {
        const kdtree_direct_leased_group_t* group =
            &groups[group_index];
        size_t cover_points =
            (size_t)(group->last - group->first) + 1U;

        selected_points += group->selected_points;
        if (cover_points > KDTREE_DIRECT_DSS_TASK_POINTS ||
            group->nspans > KDTREE_DIRECT_DSS_TASK_SPANS) {
            task_eligible = FALSE;
        }
    }
    use_executor = task_eligible && state->kd->ndim == 4 &&
        ngroups >= KDTREE_DIRECT_DSS_WAVE_MIN_TASKS &&
        selected_points >= KDTREE_DIRECT_DSS_WAVE_MIN_POINTS &&
        state->executor && state->executor->available &&
        state->executor->run &&
        state->executor->available(state->executor->opaque) > 0U;

    if (!prepared_ranges &&
        kdtree_direct_read_leased_groups(
            state, groups, ngroups, ranges)) {
        status = -1;
        goto cleanup;
    }
    loaded = ngroups;
    for (group_index = 0U;
         group_index < ngroups;
         group_index++) {
        if (use_executor && kdtree_direct_prepare_task_input(
                state,
                &groups[group_index],
                &ranges[group_index],
                &inputs[group_index])) {
            use_executor = FALSE;
        }
    }
    if (!use_executor) {
        for (group_index = 0U;
             group_index < ngroups;
             group_index++) {
            if (kdtree_direct_reduce_inline_group(
                    state,
                    &groups[group_index],
                    &ranges[group_index])) {
                status = -1;
                goto cleanup;
            }
            kdtree_direct_release_leased_range(
                state, &ranges[group_index]);
        }
        goto cleanup;
    }

    status = state->executor->run(
        state->executor->opaque,
        inputs,
        ngroups,
        &outputs);
    if (status > 0) {
        status = 0;
        for (group_index = 0U;
             group_index < ngroups;
             group_index++) {
            if (kdtree_direct_reduce_inline_group(
                    state,
                    &groups[group_index],
                    &ranges[group_index])) {
                status = -1;
                goto cleanup;
            }
            kdtree_direct_release_leased_range(
                state, &ranges[group_index]);
        }
        goto cleanup;
    }
    if (status < 0 || !outputs) {
        errno = EPROTO;
        status = -1;
        goto cleanup;
    }
    for (group_index = 0U;
         group_index < ngroups;
         group_index++) {
        if (kdtree_direct_reduce_task_output(
                state,
                &groups[group_index],
                &ranges[group_index],
                &inputs[group_index],
                &outputs[group_index])) {
            status = -1;
            goto cleanup;
        }
    }

cleanup:
    for (group_index = 0U; group_index < loaded; group_index++) {
        kdtree_direct_release_leased_range(
            state, &ranges[group_index]);
    }
    return status;
}

static int kdtree_direct_submit_leased_lookahead(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_leased_group_t* groups,
    size_t ngroups,
    void** handle) {
    kdtree_direct_range_request_t
        requests[KDTREE_DIRECT_DSS_WAVE_TASKS];
    int saved_errno = errno;
    int status;

    if (!handle) {
        errno = EINVAL;
        return -1;
    }
    *handle = NULL;
    if (!state->lookahead) {
        return 0;
    }
    if (!state->lookahead->submit ||
        !state->lookahead->finish ||
        !state->lookahead->cancel) {
        errno = EPROTO;
        return -1;
    }
    if (kdtree_direct_build_group_requests(
            state, groups, ngroups, requests)) {
        return -1;
    }
    status = state->lookahead->submit(
        state->lookahead->opaque,
        requests,
        ngroups,
        handle);
    if (status > 0 && *handle) {
        return 1;
    }
    if (*handle) {
        state->lookahead->cancel(
            state->lookahead->opaque, *handle);
        *handle = NULL;
    }
    if (status > 0) {
        errno = EPROTO;
        return -1;
    }
    if (!status) {
        errno = saved_errno;
    }
    return status;
}

static int kdtree_direct_finish_leased_lookahead(
    kdtree_direct_leased_state_t* state,
    void** handle,
    const kdtree_direct_leased_group_t* groups,
    size_t ngroups,
    kdtree_direct_range_t* ranges) {
    kdtree_direct_range_request_t
        requests[KDTREE_DIRECT_DSS_WAVE_TASKS];
    void* pending;
    size_t range_index;
    int status;

    memset(ranges,
           0,
           KDTREE_DIRECT_DSS_WAVE_TASKS * sizeof(*ranges));
    if (!handle || !*handle) {
        return 0;
    }
    pending = *handle;
    *handle = NULL;
    if (kdtree_direct_build_group_requests(
            state, groups, ngroups, requests)) {
        state->lookahead->cancel(
            state->lookahead->opaque, pending);
        return -1;
    }
    status = state->lookahead->finish(
        state->lookahead->opaque,
        pending,
        requests,
        ngroups,
        ranges);
    if (status > 0) {
        for (range_index = 0U;
             range_index < ngroups;
             range_index++) {
            if (!ranges[range_index].lease ||
                !ranges[range_index].data ||
                (state->kd->perm &&
                 !ranges[range_index].perm)) {
                errno = EPROTO;
                status = -1;
                break;
            }
        }
    }
    if (status <= 0) {
        for (range_index = 0U;
             range_index < ngroups;
             range_index++) {
            kdtree_direct_release_leased_range(
                state, &ranges[range_index]);
        }
    }
    return status;
}

/*
 * Publish the previous canonical wave, submit the next one, then reduce the
 * previous wave while its successor is loading. Publishing first lets the
 * successor share overlapping cache entries instead of issuing duplicate I/O.
 */
static int kdtree_direct_advance_leased_pipeline(
    kdtree_direct_leased_state_t* state) {
    kdtree_direct_range_t
        pending_ranges[KDTREE_DIRECT_DSS_WAVE_TASKS];
    void* next_handle = NULL;
    int next_status;
    int pending_status = 0;

    if (!state->ngroups) {
        return 0;
    }
    if (state->npending_groups) {
        pending_status = kdtree_direct_finish_leased_lookahead(
            state,
            &state->pending_handle,
            state->pending_groups,
            state->npending_groups,
            pending_ranges);
        if (pending_status < 0) {
            state->npending_groups = 0U;
            state->ngroups = 0U;
            return -1;
        }
    }
    next_status = kdtree_direct_submit_leased_lookahead(
        state,
        state->groups,
        state->ngroups,
        &next_handle);
    if (next_status < 0) {
        if (pending_status > 0) {
            size_t range_index;

            for (range_index = 0U;
                 range_index < state->npending_groups;
                 range_index++) {
                kdtree_direct_release_leased_range(
                    state, &pending_ranges[range_index]);
            }
        }
        state->npending_groups = 0U;
        state->ngroups = 0U;
        return -1;
    }
    if (state->npending_groups) {
        if (kdtree_direct_process_leased_wave(
                state,
                state->pending_groups,
                state->npending_groups,
                pending_status > 0 ? pending_ranges : NULL)) {
            if (next_handle) {
                state->lookahead->cancel(
                    state->lookahead->opaque, next_handle);
            }
            state->npending_groups = 0U;
            state->ngroups = 0U;
            return -1;
        }
        state->npending_groups = 0U;
    }
    if (next_status > 0) {
        memcpy(state->pending_groups,
               state->groups,
               state->ngroups * sizeof(state->groups[0]));
        state->npending_groups = state->ngroups;
        state->pending_handle = next_handle;
    } else if (kdtree_direct_process_leased_wave(
                   state,
                   state->groups,
                   state->ngroups,
                   NULL)) {
        state->ngroups = 0U;
        return -1;
    }
    state->ngroups = 0U;
    return 0;
}

static int kdtree_direct_finish_leased_pipeline(
    kdtree_direct_leased_state_t* state) {
    kdtree_direct_range_t
        pending_ranges[KDTREE_DIRECT_DSS_WAVE_TASKS];
    int pending_status;

    if (state->npending_groups) {
        if (state->ngroups &&
            kdtree_direct_advance_leased_pipeline(state)) {
            return -1;
        }
        if (!state->npending_groups) {
            return 0;
        }
        pending_status = kdtree_direct_finish_leased_lookahead(
            state,
            &state->pending_handle,
            state->pending_groups,
            state->npending_groups,
            pending_ranges);
        if (pending_status < 0) {
            state->npending_groups = 0U;
            return -1;
        }
        if (kdtree_direct_process_leased_wave(
                state,
                state->pending_groups,
                state->npending_groups,
                pending_status > 0 ? pending_ranges : NULL)) {
            state->npending_groups = 0U;
            return -1;
        }
        state->npending_groups = 0U;
        return 0;
    }
    if (state->ngroups) {
        int status = kdtree_direct_process_leased_wave(
            state,
            state->groups,
            state->ngroups,
            NULL);

        state->ngroups = 0U;
        return status;
    }
    return 0;
}

static void kdtree_direct_abort_leased_pipeline(
    kdtree_direct_leased_state_t* state) {
    if (state->pending_handle && state->lookahead) {
        state->lookahead->cancel(
            state->lookahead->opaque,
            state->pending_handle);
    }
    state->pending_handle = NULL;
    state->npending_groups = 0U;
    state->ngroups = 0U;
}

static int kdtree_direct_enqueue_leased_group(
    kdtree_direct_leased_state_t* state,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last) {
    kdtree_direct_leased_group_t* group;
    size_t span_index;

    if (!nspans) {
        return 0;
    }
    if (nspans > KDTREE_DIRECT_SPAN_CAPACITY ||
        state->ngroups >= KDTREE_DIRECT_DSS_WAVE_TASKS) {
        errno = EOVERFLOW;
        return -1;
    }
    if (first < 0 || last < first) {
        errno = EINVAL;
        return -1;
    }
    group = &state->groups[state->ngroups];
    memset(group, 0, sizeof(*group));
    memcpy(group->spans,
           spans,
           nspans * sizeof(kdtree_direct_span_t));
    group->nspans = nspans;
    group->first = first;
    group->last = last;
    for (span_index = 0U; span_index < nspans; span_index++) {
        const kdtree_direct_span_t* span = &spans[span_index];

        if (span->left < first || span->right < span->left ||
            span->right > last) {
            errno = EINVAL;
            return -1;
        }
        group->selected_points +=
            (size_t)(span->right - span->left) + 1U;
    }
    state->ngroups++;
    if (state->ngroups == KDTREE_DIRECT_DSS_WAVE_TASKS) {
        return kdtree_direct_advance_leased_pipeline(state);
    }
    return 0;
}

static int kdtree_direct_flush_leased_wave(
    void* opaque,
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last,
    anbool final) {
    kdtree_direct_leased_state_t* state = opaque;

    state->kd = kd;
    state->result = result;
    state->query = query;
    state->maxd2 = maxd2;
    if (nspans &&
        kdtree_direct_enqueue_leased_group(
            state, spans, nspans, first, last)) {
        return -1;
    }
    if (final) {
        return kdtree_direct_finish_leased_pipeline(state);
    }
    return 0;
}

typedef struct kdtree_direct_plan_state {
    kdtree_direct_plan_ranges_fn emit_ranges;
    void* emit_opaque;
    kdtree_direct_range_request_t
        requests[KDTREE_DIRECT_DSS_WAVE_TASKS];
    size_t nrequests;
} kdtree_direct_plan_state_t;

static int kdtree_direct_emit_plan_batch(
    kdtree_direct_plan_state_t* state) {
    int saved_errno;

    if (!state->nrequests) {
        return 0;
    }
    errno = 0;
    if (state->emit_ranges(
            state->emit_opaque,
            state->requests,
            state->nrequests)) {
        saved_errno = errno ? errno : EIO;
        state->nrequests = 0U;
        errno = saved_errno;
        return -1;
    }
    state->nrequests = 0U;
    return 0;
}

static int kdtree_direct_flush_plan(
    void* opaque,
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const etype* query,
    double maxd2,
    const kdtree_direct_span_t* spans,
    size_t nspans,
    int first,
    int last,
    anbool final) {
    kdtree_direct_plan_state_t* state = opaque;
    size_t count;
    size_t span_index;

    (void)result;
    (void)query;
    (void)maxd2;
    if (nspans) {
        if (!state || !kd || !spans ||
            nspans > KDTREE_DIRECT_SPAN_CAPACITY ||
            state->nrequests >= KDTREE_DIRECT_DSS_WAVE_TASKS ||
            first < 0 || last < first || last >= kd->ndata) {
            errno = EINVAL;
            return -1;
        }
        count = (size_t)last - (size_t)first + 1U;
        if (count > (size_t)INT_MAX) {
            errno = EOVERFLOW;
            return -1;
        }
        for (span_index = 0U; span_index < nspans; span_index++) {
            if (spans[span_index].left < first ||
                spans[span_index].right < spans[span_index].left ||
                spans[span_index].right > last) {
                errno = EINVAL;
                return -1;
            }
        }
        state->requests[state->nrequests].first = first;
        state->requests[state->nrequests].count = (int)count;
        state->nrequests++;
        if (state->nrequests == KDTREE_DIRECT_DSS_WAVE_TASKS &&
            kdtree_direct_emit_plan_batch(state)) {
            return -1;
        }
    }
    if (final && kdtree_direct_emit_plan_batch(state)) {
        return -1;
    }
    return 0;
}

static kdtree_qres_t* kdtree_rangesearch_direct_dss_core(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_flush_fn flush,
    void* flush_opaque,
    anbool produce_result,
    int* walk_status) {
    const int supported_options = KD_OPTIONS_SMALL_RADIUS |
        KD_OPTIONS_COMPUTE_DISTS |
        KD_OPTIONS_NO_RESIZE_RESULTS |
        KD_OPTIONS_USE_SPLIT;
    int nodestack[100];
    int stackpos = 0;
    int D;
    double maxdist;
    double dtlinf;
    ttype tlinf;
    ttype tquery[KDTREE_MAX_DIM];
    anbool use_tquery;
    anbool use_tsplit;
    kdtree_direct_span_t spans[KDTREE_DIRECT_SPAN_CAPACITY];
    size_t nspans = 0;
    int batch_first = -1;
    int batch_last = -1;
    anbool owns_result = FALSE;

    if (walk_status) {
        *walk_status = -1;
    }
    if (!produce_result && result) {
        errno = EINVAL;
        return NULL;
    }
    if (!kd || !query || !flush || !max_points || maxd2 < 0.0 ||
        options != supported_options || kd->treetype != KDTT_DSS ||
        !kd->split.any || kd->ndim <= 0 ||
        kd->ndim > KDTREE_MAX_DIM) {
        errno = ENOTSUP;
        return NULL;
    }

    D = kd->ndim;
    maxdist = sqrt(maxd2);
    use_tquery = ttype_query(kd, query, tquery);
    dtlinf = use_tquery ? DIST_ET(kd, maxdist, ) : 0.0;
    tlinf = use_tquery ? (ttype)ceil(dtlinf) : 0;
    use_tsplit = use_tquery && (dtlinf < TTYPE_MAX);

    if (produce_result) {
        if (result) {
            if (!result->capacity) {
                if (!resize_results(result,
                                    KDTREE_MAX_RESULTS,
                                    D,
                                    TRUE,
                                    FALSE)) {
                    errno = ENOMEM;
                    return NULL;
                }
            } else if (!resize_results(result,
                                       result->capacity,
                                       D,
                                       TRUE,
                                       FALSE)) {
                errno = ENOMEM;
                return NULL;
            }
            result->nres = 0;
        } else {
            result = CALLOC(1, sizeof(kdtree_qres_t));
            if (!result) {
                errno = ENOMEM;
                return NULL;
            }
            owns_result = TRUE;
            if (!resize_results(result,
                                KDTREE_MAX_RESULTS,
                                D,
                                TRUE,
                                FALSE)) {
                kdtree_free_query(result);
                errno = ENOMEM;
                return NULL;
            }
        }
    }

    nodestack[0] = 0;
    while (stackpos >= 0) {
        int nodeid = nodestack[stackpos--];

        if (KD_IS_LEAF(kd, nodeid)) {
            int left = kdtree_left(kd, nodeid);
            int right = kdtree_right(kd, nodeid);

            if (kdtree_direct_emit_span(kd,
                                        result,
                                        query,
                                        maxd2,
                                        left,
                                        right,
                                        max_points,
                                        merge_gap_points,
                                        flush,
                                        flush_opaque,
                                        spans,
                                        &nspans,
                                        &batch_first,
                                        &batch_last)) {
                goto bailout;
            }
            continue;
        }

        {
            int dim = -1;
            ttype split = *KD_SPLIT(kd, nodeid);

            if (kd->splitdim) {
                dim = kd->splitdim[nodeid];
            } else {
                bigint tmpsplit = split;
                dim = (int)(tmpsplit & kd->dimmask);
                split = (ttype)(tmpsplit & kd->splitmask);
            }
            if (dim < 0 || dim >= D) {
                errno = EINVAL;
                goto bailout;
            }

            if (use_tsplit) {
                if (tquery[dim] < split) {
                    if (stackpos >= 98) {
                        errno = EOVERFLOW;
                        goto bailout;
                    }
                    nodestack[++stackpos] = KD_CHILD_LEFT(nodeid);
                    if (split - tquery[dim] <= tlinf) {
                        nodestack[++stackpos] = KD_CHILD_RIGHT(nodeid);
                    }
                } else {
                    if (stackpos >= 98) {
                        errno = EOVERFLOW;
                        goto bailout;
                    }
                    nodestack[++stackpos] = KD_CHILD_RIGHT(nodeid);
                    if (tquery[dim] - split <= tlinf) {
                        nodestack[++stackpos] = KD_CHILD_LEFT(nodeid);
                    }
                }
            } else {
                dtype real_split = POINT_TE(kd, dim, split);

                if (query[dim] < real_split) {
                    if (stackpos >= 98) {
                        errno = EOVERFLOW;
                        goto bailout;
                    }
                    nodestack[++stackpos] = KD_CHILD_LEFT(nodeid);
                    if (real_split - query[dim] <= maxdist) {
                        nodestack[++stackpos] = KD_CHILD_RIGHT(nodeid);
                    }
                } else {
                    if (stackpos >= 98) {
                        errno = EOVERFLOW;
                        goto bailout;
                    }
                    nodestack[++stackpos] = KD_CHILD_RIGHT(nodeid);
                    if (query[dim] - real_split <= maxdist) {
                        nodestack[++stackpos] = KD_CHILD_LEFT(nodeid);
                    }
                }
            }
        }
    }

    if (flush(flush_opaque,
              kd,
              result,
              query,
              maxd2,
              spans,
              nspans,
              batch_first,
              batch_last,
              TRUE)) {
        goto bailout;
    }
    if (walk_status) {
        *walk_status = 0;
    }
    return result;

bailout:
    if (owns_result) {
        kdtree_free_query(result);
    } else if (result) {
        result->nres = 0;
    }
    return NULL;
}

int kdtree_rangesearch_direct_dss_plan(
    const kdtree_t* kd,
    const double* query,
    double maxd2,
    int options,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_plan_ranges_fn emit_ranges,
    void* emit_opaque) {
    kdtree_direct_plan_state_t state;
    int walk_status = -1;

    if (!emit_ranges) {
        errno = ENOTSUP;
        return -1;
    }
    memset(&state, 0, sizeof(state));
    state.emit_ranges = emit_ranges;
    state.emit_opaque = emit_opaque;
    (void)kdtree_rangesearch_direct_dss_core(
        kd,
        NULL,
        query,
        maxd2,
        options,
        max_points,
        merge_gap_points,
        kdtree_direct_flush_plan,
        &state,
        FALSE,
        &walk_status);
    return walk_status;
}

kdtree_qres_t* kdtree_rangesearch_direct_dss(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_read_range_fn read_range,
    void* read_opaque) {
    kdtree_direct_inline_reader_t reader;

    if (!read_range) {
        errno = ENOTSUP;
        return NULL;
    }
    reader.read_range = read_range;
    reader.read_opaque = read_opaque;
    return kdtree_rangesearch_direct_dss_core(
        kd,
        result,
        query,
        maxd2,
        options,
        max_points,
        merge_gap_points,
        kdtree_direct_flush_inline,
        &reader,
        TRUE,
        NULL);
}

kdtree_qres_t* kdtree_rangesearch_direct_dss_leased(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_read_leased_ranges_fn read_ranges,
    kdtree_direct_release_range_fn release_range,
    void* read_opaque,
    const kdtree_direct_dss_executor_t* executor,
    const kdtree_direct_dss_lookahead_t* lookahead) {
    kdtree_direct_leased_state_t state;
    kdtree_qres_t* searched;

    if (!read_ranges || !release_range) {
        errno = ENOTSUP;
        return NULL;
    }
    memset(&state, 0, sizeof(state));
    state.read_ranges = read_ranges;
    state.release_range = release_range;
    state.read_opaque = read_opaque;
    state.executor = executor;
    state.lookahead = lookahead;
    searched = kdtree_rangesearch_direct_dss_core(
        kd,
        result,
        query,
        maxd2,
        options,
        max_points,
        merge_gap_points,
        kdtree_direct_flush_leased_wave,
        &state,
        TRUE,
        NULL);
    kdtree_direct_abort_leased_pipeline(&state);
    return searched;
}

#undef KDTREE_DIRECT_SPAN_CAPACITY

#undef KDTREE_CODEKD_DSS_U16_FAST_PATH

#include "kdtree_internal_fits.c"
