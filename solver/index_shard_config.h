#ifndef INDEX_SHARD_CONFIG_H
#define INDEX_SHARD_CONFIG_H

#include <stddef.h>

/*
 * Public text uses "auto".  Zero is only the internal representation carried
 * through the AXY job header; zero and negative numeric user input are invalid.
 */
#define INDEX_SHARD_WORKERS_AUTO 0
#define INDEX_SHARD_WORKERS_UNSET (-1)

typedef struct index_shard_width_plan {
  size_t producer_width;
  size_t helper_width;
} index_shard_width_plan_t;

/*
 * Return the number of logical CPUs currently available to this process.
 * Linux process affinity is preferred; portable online-CPU detection is the
 * fallback.  The function always returns at least one.
 */
int index_shard_config_available_cpus(void);

/*
 * Parse exactly "auto" or a decimal integer in [1, available_cpus].
 * On success, requested_workers is AUTO or a positive explicit count.
 */
int index_shard_config_parse_workers(const char *value,
                                     int available_cpus,
                                     int *requested_workers);

/*
 * Validate an internally represented request.  AUTO is valid; UNSET and all
 * other non-positive values are invalid.
 */
int index_shard_config_validate_workers(int requested_workers,
                                        int available_cpus);

/*
 * Resolve AUTO to the affinity-visible CPU count.
 */
int index_shard_config_resolve_workers(int requested_workers,
                                       int available_cpus);

/*
 * Return the configured fixed-pool width. nindexes is retained for source
 * compatibility with the preceding implementation.
 */
int index_shard_config_effective_workers(int configured_workers,
                                         size_t nindexes);

/*
 * Bound simultaneous cold outer owners by detached delivery capacity.
 * Surplus compute threads remain available to execute staged work.
 * Without detached completion, retain the conservative fixed-helper policy.
 */
int index_shard_config_plan_widths(
    int worker_count,
    int payload_io_width,
    int detached_completion,
    index_shard_width_plan_t *plan);

#endif
