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

typedef enum kdtree_prefetch_emit_status {
    KDTREE_PREFETCH_EMIT_ERROR = -1,
    KDTREE_PREFETCH_EMIT_CONTINUE = 0,
    KDTREE_PREFETCH_EMIT_REFUSED = 1
} kdtree_prefetch_emit_status_t;

typedef enum kdtree_prefetch_prepare_status {
    KDTREE_PREFETCH_PREPARE_ERROR = -1,
    KDTREE_PREFETCH_PREPARE_NOT_APPLICABLE = 0,
    KDTREE_PREFETCH_PREPARE_COMPLETE = 1,
    KDTREE_PREFETCH_PREPARE_REFUSED = 2
} kdtree_prefetch_prepare_status_t;

typedef enum kdtree_rangesearch_span_mode {
    KDTREE_RANGESEARCH_SPAN_FILTER = 0,
    KDTREE_RANGESEARCH_SPAN_ACCEPT_ALL = 1
} kdtree_rangesearch_span_mode_t;

typedef struct kdtree_rangesearch_span {
    int left;
    int right;
    unsigned char mode;
} kdtree_rangesearch_span_t;

typedef struct kdtree_prefetch_sink {
    void *userdata;

    int (*enabled)(void *userdata, void *mapping);

    int (*emit)(void *userdata,
                const kdtree_prefetch_hint_t *hint);

    int (*flush)(void *userdata);

    /*
     * Optional logical output. It is emitted only after every DATA/PERM hint
     * covering the span has been accepted. The callback must copy the span.
     */
    int (*emit_span)(void *userdata,
                     const kdtree_rangesearch_span_t *span);
} kdtree_prefetch_sink_t;

/*
 * Mirror the scalar topology traversal and emit complete covering DATA/PERM
 * hints plus traversal-local metadata hints. No query result is produced,
 * payload is not dereferenced, and libkd issues no I/O advice.
 *
 * Return:
 *  KDTREE_PREFETCH_PREPARE_COMPLETE when every emitted hint was accepted.
 *  KDTREE_PREFETCH_PREPARE_NOT_APPLICABLE when preparation is disabled.
 *  KDTREE_PREFETCH_PREPARE_REFUSED when the sink refuses a hint before full
 *  traversal completion; previously accepted hints are an incomplete prefix.
 *  KDTREE_PREFETCH_PREPARE_ERROR for invalid input or traversal failure.
 */
int kdtree_rangesearch_prefetch_prepare(
    const kdtree_t *kd,
    const void *query,
    double maxd2,
    int options,
    const kdtree_prefetch_sink_t *sink);

/*
 * Execute one canonical span program without traversing KD topology again.
 * All spans are validated before result storage is changed. The caller must
 * supply spans emitted by a complete preparation using the same tree, query,
 * radius, and options. On unsupported input this returns NULL, allowing the
 * authoritative scalar search to run unchanged.
 */
kdtree_qres_t *kdtree_rangesearch_execute_spans(
    const kdtree_t *kd,
    kdtree_qres_t *result,
    const void *query,
    double maxd2,
    int options,
    const kdtree_rangesearch_span_t *spans,
    size_t span_count);

#endif
