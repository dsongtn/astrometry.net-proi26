/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef ASTROMETRY_ENGINE_PRIVATE_H
#define ASTROMETRY_ENGINE_PRIVATE_H

#include "astrometry/engine.h"
#include "astrometry/index_residency.h"
#include "astrometry/onefield.h"

/*
 * Private engine cross-module boundary.
 *
 * engine_job.c owns job parsing and image dimensions; engine_policy.c resolves
 * limits and requested worker width; engine_residency.c owns optional cohort
 * admission. None of these declarations transfers index or solver ownership.
 * Keep pass-only types and functions in engine_internal.h.
 */

double engine_job_imagew(job_t* job);
double engine_job_imageh(job_t* job);

int engine_resolve_index_shard_workers(engine_t* engine, job_t* job);

index_residency_t* engine_index_residency_begin(
    engine_t* engine,
    const onefield_t* bp);

#endif
