#ifndef ASTROMETRY_ENGINE_PRIVATE_H
#define ASTROMETRY_ENGINE_PRIVATE_H

#include "astrometry/engine.h"
#include "astrometry/index_residency.h"
#include "astrometry/onefield.h"

double engine_job_imagew(job_t* job);
double engine_job_imageh(job_t* job);

int engine_resolve_index_shard_workers(engine_t* engine, job_t* job);

index_residency_t* engine_index_residency_begin(
    engine_t* engine,
    const onefield_t* bp);

#endif
