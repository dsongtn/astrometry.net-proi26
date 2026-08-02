/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine_internal.h"
#include "engine_private.h"
#include "errors.h"
#include "index_shard_config.h"
#include "log.h"
#include "mathutil.h"

void engine_limit_policy_resolve(double job_wall_seconds,
                                 double config_wall_seconds,
                                 double job_cpu_seconds,
                                 double config_cpu_seconds,
                                 engine_limit_policy_t* policy) {
    if (!policy) {
        return;
    }
    memset(policy, 0, sizeof(*policy));

    if (job_wall_seconds > 0.0 && config_wall_seconds > 0.0) {
        policy->wall_seconds = MIN(job_wall_seconds, config_wall_seconds);
        policy->wall_from_job = job_wall_seconds <= config_wall_seconds;
        policy->wall_job_clamped = job_wall_seconds > config_wall_seconds;
    } else if (job_wall_seconds > 0.0) {
        policy->wall_seconds = job_wall_seconds;
        policy->wall_from_job = TRUE;
    } else if (config_wall_seconds > 0.0) {
        policy->wall_seconds = config_wall_seconds;
    }

    if (job_cpu_seconds > 0.0) {
        policy->cpu_seconds = job_cpu_seconds;
        policy->cpu_from_job = TRUE;
    } else if (config_cpu_seconds > 0.0) {
        policy->cpu_seconds = config_cpu_seconds;
    }
}

int engine_resolve_index_shard_workers(engine_t *engine,
                                       job_t *job) {
    const char *environment_value;
    const char *source;
    int available_cpus;
    int requested_workers;
    int resolved_workers;
    char requested_text[32];

    if (!engine || !job) {
        ERROR("Cannot resolve parallel workers without engine and job state");
        return -1;
    }

    available_cpus = index_shard_config_available_cpus();
    requested_workers = engine->index_shard_workers_config;
    source = engine->index_shard_workers_config_set
        ? "config"
        : "built-in";

    if (job->index_shard_workers_override_set) {
        requested_workers = job->index_shard_workers_override;
        if (index_shard_config_validate_workers(requested_workers,
                                                available_cpus)) {
            ERROR("Invalid ANSHWRK worker override %i: "
                  "expected automatic selection or an integer from "
                  "1 through %i",
                  requested_workers,
                  available_cpus);
            return -1;
        }
        source = "solve-field";
    } else {
        /*
         * Retain the legacy environment override for existing measurement
         * harnesses. New production commands should use the per-job
         * solve-field option, whose AXY header has higher precedence.
         */
        environment_value = getenv("ASTROMETRY_P_WORKERS");
        if (!environment_value || !environment_value[0]) {
            environment_value = getenv("ASTROMETRY_INDEX_SHARD_WORKERS");
        }
        if (environment_value && environment_value[0]) {
            if (index_shard_config_parse_workers(environment_value,
                                                 available_cpus,
                                                 &requested_workers)) {
                ERROR("Invalid parallel worker environment value \"%s\": "
                      "expected \"auto\" or an integer from 1 through %i",
                      environment_value,
                      available_cpus);
                return -1;
            }
            source = "environment";
        }
    }

    resolved_workers =
        index_shard_config_resolve_workers(requested_workers,
                                           available_cpus);
    if (resolved_workers < 1) {
        ERROR("Failed to resolve parallel worker count");
        return -1;
    }

    job->bp.index_shard_workers = resolved_workers;
    job->index_shard_workers_controlled =
        strcmp(source, "built-in") != 0;

    if (requested_workers == INDEX_SHARD_WORKERS_AUTO) {
        snprintf(requested_text, sizeof(requested_text), "auto");
    } else {
        snprintf(requested_text,
                 sizeof(requested_text),
                 "%i",
                 requested_workers);
    }

    logverb("[parallel] worker-config source=%s requested=%s "
            "available=%i effective=%i mode=%s\n",
            source,
            requested_text,
            available_cpus,
            resolved_workers,
            resolved_workers > 1 ? "pthread" : "serial");

    return 0;
}
