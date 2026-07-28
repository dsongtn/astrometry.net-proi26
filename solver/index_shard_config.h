#ifndef INDEX_SHARD_CONFIG_H
#define INDEX_SHARD_CONFIG_H

#include <stddef.h>

/*
 * Public text uses "auto".  Zero is only the internal representation carried
 * through the AXY job header; zero and negative numeric user input are invalid.
 */
#define INDEX_SHARD_WORKERS_AUTO 0
#define INDEX_SHARD_WORKERS_UNSET (-1)

/*
 * Compatibility topology retained for focused configuration tests. Production
 * assigns every physical worker a full outer-index owner role and does not
 * create phase-assistance groups.
 */
typedef struct index_shard_phase_group_config {
  size_t index_order;
  int first_worker;
  int worker_count;
} index_shard_phase_group_config_t;

/*
 * Compatibility state for the measured owner-credit controller.
 *
 * Owner-only resource samples do not steer production topology because they
 * omit work performed elsewhere. Production now uses full outer ownership;
 * these pure transitions remain available only to focused tests and possible
 * future whole-pool measurement work.
 */
typedef struct index_shard_resource_sample {
  double user_seconds;
  double system_seconds;
  double wall_seconds;
  unsigned long long major_faults;
  unsigned long long input_blocks;
  unsigned long long voluntary_switches;
  int completed;
  int cancelled;
  int failed;
} index_shard_resource_sample_t;

typedef struct index_shard_admission_controller {
  int worker_count;
  int owner_target;
  int owner_pressure_score;
  unsigned long long samples;
  unsigned long long io_samples;
  unsigned long long cpu_samples;
  unsigned long long owner_target_decreases;
  unsigned long long owner_target_increases;
} index_shard_admission_controller_t;

typedef struct index_shard_admission_update {
  int owner_target_changed;
  int io_bound;
  int cpu_bound;
  int read_amplified;
  int fault_stalled;
} index_shard_admission_update_t;

void index_shard_admission_init(
    index_shard_admission_controller_t *controller,
    int worker_count);

/*
 * Record one non-overlapping outer-owner thread sample.
 *
 * The controller uses hysteresis. A malformed, cancelled, failed, or
 * incomplete sample is ignored; the fixed NORMAL payload policy is therefore
 * never stranded behind an incomplete feedback transition.
 */
index_shard_admission_update_t index_shard_admission_record(
    index_shard_admission_controller_t *controller,
    const index_shard_resource_sample_t *sample);

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
 * Phase assistance is a pool-capacity decision only. Field-object ordinals,
 * configured range order, index identity, and cache state are deliberately
 * absent from this predicate.
 */
int index_shard_config_phase_assist_eligible(int worker_count,
                                             size_t nindexes);

/*
 * Divide the complete fixed pool into one disjoint group per index.  The
 * quotient/remainder rule is deterministic and satisfies:
 *
 *   nindexes + sum(group.worker_count - 1) == worker_count.
 *
 * Returns zero on success and -1 for an ineligible or invalid request.
 */
int index_shard_config_build_phase_groups(
    int worker_count,
    size_t nindexes,
    index_shard_phase_group_config_t *groups,
    size_t group_capacity,
    size_t *worker_group,
    int *worker_local_rank);

#endif
