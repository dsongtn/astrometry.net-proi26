/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef KDTREE_DIRECT_INTERNAL_H
#define KDTREE_DIRECT_INTERNAL_H

#include <stddef.h>

#include "astrometry/kdtree.h"

/*
 * Reads a contiguous range of points from the original CodeKD payload.
 * The returned pointers remain valid until the next callback on this thread.
 */
typedef int (*kdtree_direct_read_range_fn)(
    void* opaque,
    int first,
    int count,
    const u16** data,
    const u32** perm);

#define KDTREE_DIRECT_DSS_TASK_POINTS 4096U
#define KDTREE_DIRECT_DSS_TASK_SPANS 128U
#define KDTREE_DIRECT_DSS_WAVE_TASKS 6U
#define KDTREE_DIRECT_DSS_WAVE_MIN_TASKS 2U
#define KDTREE_DIRECT_DSS_WAVE_MIN_POINTS 8192U

typedef struct kdtree_direct_range {
    const u16* data;
    const u32* perm;
    void* lease;
} kdtree_direct_range_t;

typedef struct kdtree_direct_range_request {
    int first;
    int count;
} kdtree_direct_range_request_t;

/*
 * Success fills every range with non-NULL data and lease, plus perm when the
 * tree has one. The pointers remain valid until release_range() receives each
 * lease. Failure may return partial non-NULL leases; libkd releases each one
 * exactly once. The callback completes all demand I/O before it returns.
 */
typedef int (*kdtree_direct_read_leased_ranges_fn)(
    void* opaque,
    const kdtree_direct_range_request_t* requests,
    size_t nrequests,
    kdtree_direct_range_t* ranges);

typedef void (*kdtree_direct_release_range_fn)(
    void* opaque,
    void* lease);

typedef struct kdtree_direct_dss_task_input {
    const u16* data;
    double query[4];
    double minval[4];
    double invscale;
    double maxd2;
    u16 span_first[KDTREE_DIRECT_DSS_TASK_SPANS];
    u16 span_last[KDTREE_DIRECT_DSS_TASK_SPANS];
    u16 nspans;
    u32 cover_points;
} kdtree_direct_dss_task_input_t;

typedef struct kdtree_direct_dss_task_output {
    u32 nres;
    u16 offsets[KDTREE_DIRECT_DSS_TASK_POINTS];
    double sdists[KDTREE_DIRECT_DSS_TASK_POINTS];
    int status;
} kdtree_direct_dss_task_output_t;

typedef struct kdtree_direct_dss_executor {
    void* opaque;
    size_t (*available)(void* opaque);
    /*
     * Zero supplies ordered outputs valid through synchronous reduction.
     * Positive requests exact owner-inline recomputation. Negative reports a
     * hard execution or invariant failure.
     */
    int (*run)(
        void* opaque,
        const kdtree_direct_dss_task_input_t* inputs,
        size_t ntasks,
        const kdtree_direct_dss_task_output_t** outputs);
} kdtree_direct_dss_executor_t;

int kdtree_direct_dss_task_execute(
    const kdtree_direct_dss_task_input_t* input,
    kdtree_direct_dss_task_output_t* output);

/*
 * Exact split-tree search for the double/U16/U16 CodeKD representation.
 * Unsupported calls return NULL with errno=ENOTSUP. I/O failures return NULL
 * with the callback's errno and leave the caller free to execute the legacy
 * mapped search from the beginning.
 */
kdtree_qres_t* kdtree_rangesearch_direct_dss(
    const kdtree_t* kd,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options,
    size_t max_points,
    size_t merge_gap_points,
    kdtree_direct_read_range_fn read_range,
    void* read_opaque);

/*
 * The lease-aware path uses the same exact traversal and releases every
 * acquired lease before return. Positive executor status recomputes the wave
 * owner-inline; negative status is returned as a hard direct-path failure.
 */
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
    const kdtree_direct_dss_executor_t* executor);

#endif
