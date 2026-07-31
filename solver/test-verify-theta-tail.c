/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "verify_theta_tail.h"

#define ARRAY_SIZE(array) (sizeof(array) / sizeof((array)[0]))

static int expect_value(const char* case_name,
                        const int* values,
                        size_t index,
                        int expected) {
    if (values[index] == expected) {
        return 0;
    }
    fprintf(stderr,
            "%s: theta[%zu]=%d, expected %d\n",
            case_name, index, values[index], expected);
    return -1;
}

static int test_bailout_tail(void) {
    int theta[] = { 0, 1, 2, INT_MAX, INT_MAX, INT_MAX };
    size_t i;

    if (verify_theta_mark_unprocessed(
            theta, (int)ARRAY_SIZE(theta), 2, -1)) {
        return -1;
    }
    for (i = 0; i <= 2; i++) {
        if (expect_value(
                "bailout-prefix", theta, i, (int)i)) {
            return -1;
        }
    }
    for (i = 3; i < ARRAY_SIZE(theta); i++) {
        if (expect_value(
                "bailout-tail", theta, i, THETA_BAILEDOUT)) {
            return -1;
        }
    }
    return 0;
}

static int test_stop_tail(void) {
    int theta[] = { 3, 2, 1, 0, INT_MAX, INT_MAX };
    const int prefix[] = { 3, 2, 1, 0 };
    size_t i;

    if (verify_theta_mark_unprocessed(
            theta, (int)ARRAY_SIZE(theta), -1, 3)) {
        return -1;
    }
    for (i = 0; i < ARRAY_SIZE(prefix); i++) {
        if (expect_value("stop-prefix", theta, i, prefix[i])) {
            return -1;
        }
    }
    for (i = ARRAY_SIZE(prefix); i < ARRAY_SIZE(theta); i++) {
        if (expect_value(
                "stop-tail", theta, i,
                THETA_STOPPEDLOOKING)) {
            return -1;
        }
    }
    return 0;
}

static int test_historical_overlap_order(void) {
    int theta[] = { 7, 6, INT_MAX, INT_MAX, INT_MAX, INT_MAX };

    if (verify_theta_mark_unprocessed(
            theta, (int)ARRAY_SIZE(theta), 1, 3)) {
        return -1;
    }
    if (expect_value("overlap-prefix-0", theta, 0, 7) ||
        expect_value("overlap-prefix-1", theta, 1, 6) ||
        expect_value(
            "overlap-bail-2", theta, 2, THETA_BAILEDOUT) ||
        expect_value(
            "overlap-bail-3", theta, 3, THETA_BAILEDOUT) ||
        expect_value(
            "overlap-stop-4", theta, 4,
            THETA_STOPPEDLOOKING) ||
        expect_value(
            "overlap-stop-5", theta, 5,
            THETA_STOPPEDLOOKING)) {
        return -1;
    }
    return 0;
}

static int test_no_termination_marker(void) {
    int theta[] = { 4, 3, 2, 1 };
    int expected[ARRAY_SIZE(theta)];

    memcpy(expected, theta, sizeof(theta));
    if (verify_theta_mark_unprocessed(
            theta, (int)ARRAY_SIZE(theta), -1, -1)) {
        return -1;
    }
    if (memcmp(theta, expected, sizeof(theta))) {
        fprintf(stderr,
                "no-marker: theta changed unexpectedly\n");
        return -1;
    }
    return 0;
}

static int test_invalid_arguments(void) {
    int theta[] = { 0, 1 };

    if (!verify_theta_mark_unprocessed(NULL, 1, -1, -1) ||
        !verify_theta_mark_unprocessed(theta, -1, -1, -1) ||
        !verify_theta_mark_unprocessed(theta, 2, 2, -1) ||
        !verify_theta_mark_unprocessed(theta, 2, -1, 2) ||
        verify_theta_mark_unprocessed(NULL, 0, -1, -1)) {
        fprintf(stderr,
                "invalid-arguments: validation contract failed\n");
        return -1;
    }
    return 0;
}

int main(void) {
    if (test_bailout_tail() ||
        test_stop_tail() ||
        test_historical_overlap_order() ||
        test_no_termination_marker() ||
        test_invalid_arguments()) {
        return 1;
    }
    printf("VERIFY_THETA_TAIL_OK\n");
    return 0;
}
