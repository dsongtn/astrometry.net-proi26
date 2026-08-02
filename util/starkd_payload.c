/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "starkd.h"
#include "starkd_internal.h"
#include "kdtree.h"
#include "kdtree_fits_io.h"
#include "starutil.h"
#include "fitsbin.h"
#include "fitstable.h"
#include "errors.h"
#include "tic.h"
#include "log.h"
#include "ioutils.h"
#include "fitsioutils.h"

int startree_prepare_stars(startree_t* s,
                            const unsigned int* starids,
                            int nstars) {
    fitsbin_prefetch_range_t ranges[160];
    fitsbin_t* fb;
    size_t row_size;
    size_t byte_budget;
    size_t per_range_budget;
    long detected_page_size;
    int i;

    if (!s || !s->tree || !starids || nstars <= 0 ||
        !s->tree->io || !s->tree->io_is_fitsbin ||
        !s->tree->data.any || startree_data_count_internal(s) <= 0) {
        return 0;
    }
    if ((size_t)nstars > sizeof(ranges) / sizeof(ranges[0])) {
        errno = E2BIG;
        return -1;
    }
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0) {
        errno = EINVAL;
        return -1;
    }
    fb = s->tree->io;
    row_size = kdtree_sizeof_data(s->tree) / (size_t)startree_data_count_internal(s);
    if (!row_size || (size_t)detected_page_size >
        (SIZE_MAX - row_size) / 2U) {
        errno = EOVERFLOW;
        return -1;
    }
    per_range_budget = row_size +
        2U * (size_t)detected_page_size;
    if ((size_t)nstars > SIZE_MAX / per_range_budget) {
        errno = EOVERFLOW;
        return -1;
    }
    byte_budget = (size_t)nstars * per_range_budget;

    for (i = 0; i < nstars; i++) {
        int data_index;

        if (starids[i] >= (unsigned int)startree_data_count_internal(s)) {
            errno = EINVAL;
            return -1;
        }
        data_index = startree_data_index_internal(s, (int)starids[i]);
        if (data_index < 0) {
            return -1;
        }
        ranges[i].data = kdtree_get_data(s->tree, data_index);
        ranges[i].size = row_size;
    }
    return fitsbin_prefetch_ranges(
        fb,
        ranges,
        (size_t)nstars,
        byte_budget);
}

int startree_prefetch_stars(startree_t* s,
                            const unsigned int* starids,
                            int nstars) {
    int status = startree_prepare_stars(
        s, starids, nstars);

    return status < 0 ? -1 : 0;
}

int startree_prefetch_stars_submit(
    startree_t* s,
    const unsigned int* starids,
    int nstars,
    fitsbin_payload_io_ticket_t** ticket) {
    fitsbin_prefetch_range_t ranges[160];
    fitsbin_t* fb;
    size_t row_size;
    size_t byte_budget;
    size_t per_range_budget;
    long detected_page_size;
    int i;

    if (!ticket) {
        errno = EINVAL;
        return -1;
    }
    *ticket = NULL;
    if (!s || !s->tree || !starids || nstars <= 0 ||
        !s->tree->io || !s->tree->io_is_fitsbin ||
        !s->tree->data.any || startree_data_count_internal(s) <= 0) {
        return 0;
    }
    if ((size_t)nstars > sizeof(ranges) / sizeof(ranges[0])) {
        errno = E2BIG;
        return -1;
    }
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0) {
        errno = EINVAL;
        return -1;
    }
    row_size = kdtree_sizeof_data(s->tree) / (size_t)startree_data_count_internal(s);
    if (!row_size || (size_t)detected_page_size >
        (SIZE_MAX - row_size) / 2U) {
        errno = EOVERFLOW;
        return -1;
    }
    per_range_budget = row_size +
        2U * (size_t)detected_page_size;
    if ((size_t)nstars > SIZE_MAX / per_range_budget) {
        errno = EOVERFLOW;
        return -1;
    }
    byte_budget = (size_t)nstars * per_range_budget;
    fb = s->tree->io;

    for (i = 0; i < nstars; i++) {
        int data_index;

        if (starids[i] >= (unsigned int)startree_data_count_internal(s)) {
            errno = EINVAL;
            return -1;
        }
        data_index = startree_data_index_internal(s, (int)starids[i]);
        if (data_index < 0) {
            return -1;
        }
        ranges[i].data = kdtree_get_data(s->tree, data_index);
        ranges[i].size = row_size;
    }
    return fitsbin_prefetch_ranges_submit(
        fb,
        ranges,
        (size_t)nstars,
        byte_budget,
        ticket);
}

int startree_prefetch_stars_ready_submit(
    const startree_t* s,
    const unsigned int* starids,
    int nstars,
    fitsbin_payload_io_ticket_t** ticket) {
    fitsbin_prefetch_range_t
        ranges[FITSBIN_MMAP_PREFETCH_RANGE_LIMIT];
    fitsbin_t* fb;
    size_t row_size;
    size_t byte_budget;
    size_t per_range_budget;
    long detected_page_size;
    int ndata;
    int i;

    if (!ticket) {
        errno = EINVAL;
        return -1;
    }
    *ticket = NULL;
    if (!s || !s->tree || !starids || nstars <= 0 ||
        !s->tree->io || !s->tree->io_is_fitsbin ||
        !s->tree->data.any) {
        return 0;
    }
    ndata = startree_data_count_internal(s);
    if (ndata <= 0) {
        return 0;
    }
    if ((size_t)nstars > sizeof(ranges) / sizeof(ranges[0])) {
        errno = E2BIG;
        return -1;
    }
    if (s->tree->perm && !s->inverse_perm) {
        errno = EAGAIN;
        return 0;
    }
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0) {
        errno = EINVAL;
        return -1;
    }
    row_size = kdtree_sizeof_data(s->tree) / (size_t)ndata;
    if (!row_size || (size_t)detected_page_size >
        (SIZE_MAX - row_size) / 2U) {
        errno = EOVERFLOW;
        return -1;
    }
    per_range_budget = row_size +
        2U * (size_t)detected_page_size;
    if ((size_t)nstars > SIZE_MAX / per_range_budget) {
        errno = EOVERFLOW;
        return -1;
    }
    byte_budget = (size_t)nstars * per_range_budget;
    fb = s->tree->io;

    for (i = 0; i < nstars; i++) {
        unsigned int starid = starids[i];
        int data_index;

        if (starid >= (unsigned int)ndata) {
            errno = EINVAL;
            return -1;
        }
        data_index = s->inverse_perm
            ? s->inverse_perm[starid]
            : (int)starid;
        if (data_index < 0 || data_index >= ndata) {
            errno = EINVAL;
            return -1;
        }
        ranges[i].data = kdtree_get_data(s->tree, data_index);
        ranges[i].size = row_size;
    }
    return fitsbin_prefetch_ranges_submit(
        fb,
        ranges,
        (size_t)nstars,
        byte_budget,
        ticket);
}

int startree_get_ready(
    const startree_t* s,
    int starid,
    double* posn) {
    int data_index;
    int ndata;

    if (!s || !s->tree || !posn) {
        errno = EINVAL;
        return -1;
    }
    ndata = startree_data_count_internal(s);
    if (starid < 0 || starid >= ndata) {
        errno = EINVAL;
        return -1;
    }
    if (s->tree->perm && !s->inverse_perm) {
        errno = EAGAIN;
        return -1;
    }
    data_index = s->inverse_perm
        ? s->inverse_perm[starid]
        : starid;
    if (data_index < 0 || data_index >= ndata) {
        errno = EINVAL;
        return -1;
    }
    kdtree_copy_data_double(s->tree, data_index, 1, posn);
    return 0;
}

#define STARTREE_MAPPED_ADVICE_RANGE_LIMIT 256U

int startree_advise_rows(startree_t* s,
                         const unsigned int* starids,
                         int nstars) {
    fitsbin_prefetch_range_t
        ranges[STARTREE_MAPPED_ADVICE_RANGE_LIMIT];
    fitsbin_t* fb;
    size_t accepted;
    size_t byte_budget;
    size_t per_range_budget;
    size_t row_size;
    long detected_page_size;
    int advised;
    int ndata;
    size_t i;

    if (!s || !s->tree || !starids || nstars <= 0 ||
        !s->tree->io || !s->tree->io_is_fitsbin ||
        !s->tree->data.any) {
        return 0;
    }
    ndata = startree_data_count_internal(s);
    if (ndata <= 0) {
        return 0;
    }

    /*
     * Advice must not turn into a compulsory full PERM sweep. The normal
     * data lookup remains responsible for constructing the inverse mapping
     * when it is genuinely needed.
     */
    if (s->tree->perm && !s->inverse_perm) {
        return 0;
    }
    row_size =
        kdtree_sizeof_data(s->tree) / (size_t)ndata;
    if (!row_size) {
        return 0;
    }
    accepted = MIN(
        (size_t)nstars,
        (size_t)STARTREE_MAPPED_ADVICE_RANGE_LIMIT);
    detected_page_size = sysconf(_SC_PAGESIZE);
    if (detected_page_size <= 0 ||
        (size_t)detected_page_size >
            (SIZE_MAX - row_size) / 2U) {
        return 0;
    }
    per_range_budget = row_size +
        2U * (size_t)detected_page_size;
    if (accepted > SIZE_MAX / per_range_budget) {
        return 0;
    }
    byte_budget = accepted * per_range_budget;
    fb = s->tree->io;

    /*
     * Populate only a fixed canonical lookahead. Later rows retain the base
     * mapping policy and are consumed by the original loop on demand.
     */
    for (i = 0U; i < accepted; i++) {
        unsigned int starid = starids[i];
        int data_index;

        if (starid >= (unsigned int)ndata) {
            return -1;
        }
        data_index = s->inverse_perm
            ? s->inverse_perm[starid]
            : (int)starid;
        if (data_index < 0 || data_index >= ndata) {
            return -1;
        }
        ranges[i].data =
            kdtree_get_data(s->tree, data_index);
        ranges[i].size = row_size;
    }
    advised = fitsbin_advise_mapped_ranges(
        fb,
        ranges,
        accepted,
        byte_budget);
    return advised;
}
