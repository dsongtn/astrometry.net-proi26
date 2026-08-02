/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef SOLVER_PROFILE_INTERNAL_H
#define SOLVER_PROFILE_INTERNAL_H

#include "solver.h"

/* Reporting-only boundary; solver_profile_t remains owned by its solver. */

void solver_profile_report(const solver_t* solver);

#endif
