/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef VERIFY_PREPARED_INTERNAL_H
#define VERIFY_PREPARED_INTERNAL_H

#include <stdint.h>

#include "verify_internal.h"

/*
 * Prepared-object ownership summary.
 *
 * verify_index_query owns refxyz, refstarid, and captured sweep storage while
 * borrowing source. verify_prepared_hit owns its index-derived arrays and WCS
 * value; verify.testxy may borrow the pass field. Destruction must therefore
 * complete before the field view or source lifetime ends.
 */

typedef enum verify_prepared_state {
    VERIFY_PREPARED_READY = 0,
    VERIFY_PREPARED_NO_REFERENCE = 1,
    VERIFY_PREPARED_NO_QUAD_REFERENCE = 2,
    VERIFY_PREPARED_NO_ROR_REFERENCE = 3,
    VERIFY_PREPARED_EMPTY_LISTS = 4
} verify_prepared_state_t;

struct verify_index_query {
    /* source is borrowed until this query is consumed or destroyed. */
    const startree_t* source;
    /* Captured preparation validates ids without dereferencing source. */
    int source_nstars;
    double center[3];
    double radius2;
    double* refxyz;
    int* refstarid;
    uint8_t* sweep;
    int nrall;
};

/*
 * Continue only from a fully captured query. This private path never reads
 * source, its tree, or mapped sweep storage, so a foreign compute worker may
 * prepare the packet-private context while the owner retains index lifetime.
 */
int verify_prepare_captured_hit_from_query(
    verify_index_query_t** query,
    int index_cutnside,
    const MatchObj* mo,
    const sip_t* sip,
    const verify_field_t* vf,
    double pix2,
    double distractors,
    double fieldW,
    double fieldH,
    double logbail,
    double logaccept,
    double logstoplooking,
    anbool do_gamma,
    anbool fake_match,
    verify_prepared_hit_t** prepared);

struct verify_prepared_hit {
    /* Index-backed arrays are owned; verify.testxy borrows verify_field data. */
    verify_t verify;
    sip_t wcs;
    double* refxyz;
    double effective_area;
    double distractors;
    double logbail;
    double logaccept;
    double logstoplooking;
    int nrimage;
    verify_prepared_state_t state;
    anbool fake_match;
};

#endif
