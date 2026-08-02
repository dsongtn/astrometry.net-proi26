/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef TEST_SOLVER_PRIVATE_H
#define TEST_SOLVER_PRIVATE_H

#include "bl.h"
#include "starxy.h"

/* Shared geometry fixture used by the try-all-codes test binary. */
starxy_t* test_solver_geometry_field(void);
/* The caller owns the list and keeps it alive until solver_run() returns. */
void test_solver_geometry_use_quadlist(bl* quadlist);

void test_solver_geometry_cache_exact(void);
void test_solver_geometry_cache_deep_admission(void);

void test_solver_ab_counter_boundaries(void);
void test_solver_ab_descriptor_partition_count(void);
void test_solver_index_close_fds_failure_state(void);
void test_solver_zero_initialized_payload_fd_is_unowned(void);

#endif
