/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef VERIFY_THETA_TAIL_H
#define VERIFY_THETA_TAIL_H

#include "verify.h"

/*
 * verify_internal_star_lists() stops immediately after the star that
 * crosses a bailout or stop-looking threshold. The remaining theta
 * entries have never been written and must receive their semantic
 * sentinel before any full-array validation or permutation pass
 * reads them.
 *
 * When both markers are present, preserve the historical ordering:
 * bailout markers are written first and stop-looking markers
 * overwrite the later suffix.
 */
static inline int verify_theta_mark_unprocessed(
    int* theta,
    int count,
    int ibailed,
    int istopped) {
    int i;

    if (count < 0 ||
        (count && !theta) ||
        ibailed < -1 || ibailed >= count ||
        istopped < -1 || istopped >= count) {
        return -1;
    }

    if (ibailed != -1) {
        for (i = ibailed + 1; i < count; i++) {
            theta[i] = THETA_BAILEDOUT;
        }
    }

    if (istopped != -1) {
        for (i = istopped + 1; i < count; i++) {
            theta[i] = THETA_STOPPEDLOOKING;
        }
    }

    return 0;
}

#endif
