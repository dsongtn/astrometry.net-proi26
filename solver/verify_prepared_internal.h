/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef VERIFY_PREPARED_INTERNAL_H
#define VERIFY_PREPARED_INTERNAL_H

#include <stdint.h>

#include "verify_internal.h"

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
    double center[3];
    double radius2;
    double* refxyz;
    int* refstarid;
    uint8_t* sweep;
    int nrall;
};

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
