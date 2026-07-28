#ifndef KDTREE_PREFETCH_INTERNAL_H
#define KDTREE_PREFETCH_INTERNAL_H

#include <stddef.h>

#include "astrometry/kdtree.h"

/*
 * Logical source array for one prefetch hint.
 *
 * The coordinator uses this only for priority/accounting. The actual mmap
 * identity is resolved independently from the hinted address.
 */
typedef enum kdtree_prefetch_array_kind {
    KDTREE_PREFETCH_ARRAY_SPLIT = 0,
    KDTREE_PREFETCH_ARRAY_SPLITDIM,
    KDTREE_PREFETCH_ARRAY_BBOX,
    KDTREE_PREFETCH_ARRAY_LR,
    KDTREE_PREFETCH_ARRAY_DATA,
    KDTREE_PREFETCH_ARRAY_PERM
} kdtree_prefetch_array_kind_t;

/*
 * Lower numeric values are more urgent.
 */
typedef enum kdtree_prefetch_priority {
    KDTREE_PREFETCH_PRIORITY_METADATA = 0,
    KDTREE_PREFETCH_PRIORITY_LEAF = 1
} kdtree_prefetch_priority_t;

/*
 * One immutable array-specific hint.
 *
 * mapping is an opaque mapping-owner identity. For normal index-backed
 * kdtree objects this is kd->io, currently backed by fitsbin_t.
 *
 * The emitter copies the hint before returning; the caller retains ownership.
 */
typedef struct kdtree_prefetch_hint {
    void *mapping;
    const void *address;
    size_t length;

    kdtree_prefetch_array_kind_t kind;
    unsigned int priority;
} kdtree_prefetch_hint_t;

#define KDTREE_PREFETCH_ARRAY_MASK(kind) \
    (1U << (unsigned int)(kind))

typedef struct kdtree_prefetch_sink {
    void *userdata;

    int (*enabled)(void *userdata, void *mapping);

    int (*emit)(void *userdata,
                const kdtree_prefetch_hint_t *hint);

    int (*flush)(void *userdata);
} kdtree_prefetch_sink_t;

/*
 * Mirror the scalar topology traversal and emit array-specific hints for the
 * selected leaf DATA/PERM ranges and traversal metadata. No query result is
 * produced, payload is not dereferenced, and libkd issues no I/O advice.
 *
 * Return:
 *   0  normal completion or prefetch not applicable
 *  -1  invalid arguments
 */
int kdtree_rangesearch_prefetch_prepare(
    const kdtree_t *kd,
    const void *query,
    double maxd2,
    int options,
    const kdtree_prefetch_sink_t *sink);

#endif
