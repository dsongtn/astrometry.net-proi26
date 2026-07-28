#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef __linux__
#include <sched.h>
#endif

#include "index_shard_config.h"

#define INDEX_SHARD_MIN_SAMPLE_BLOCKS 8192ULL
#define INDEX_SHARD_AMPLIFIED_BLOCKS_PER_FAULT 512ULL
#define INDEX_SHARD_EXACT_BLOCKS_PER_FAULT 128ULL
#define INDEX_SHARD_MIN_FAULT_STALLS 32ULL
#define INDEX_SHARD_PRESSURE_LIMIT 4

static int index_shard_saturating_add(
    int value,
    int delta,
    int limit) {
  if (delta > 0 && value > limit - delta) {
    return limit;
  }
  if (delta < 0 && value < -limit - delta) {
    return -limit;
  }
  value += delta;
  if (value > limit) {
    return limit;
  }
  if (value < -limit) {
    return -limit;
  }
  return value;
}

static int index_shard_sample_read_amplified(
    const index_shard_resource_sample_t *sample) {
  if (!sample ||
      sample->major_faults == 0ULL ||
      sample->input_blocks < INDEX_SHARD_MIN_SAMPLE_BLOCKS ||
      sample->system_seconds <= sample->user_seconds) {
    return 0;
  }
  return sample->input_blocks / sample->major_faults >=
      INDEX_SHARD_AMPLIFIED_BLOCKS_PER_FAULT;
}

static int index_shard_sample_fault_stalled(
    const index_shard_resource_sample_t *sample) {
  if (!sample ||
      sample->major_faults < INDEX_SHARD_MIN_FAULT_STALLS ||
      sample->input_blocks / sample->major_faults >
          INDEX_SHARD_EXACT_BLOCKS_PER_FAULT) {
    return 0;
  }
  return sample->voluntary_switches >=
      sample->major_faults / 2ULL;
}

void index_shard_admission_init(
    index_shard_admission_controller_t *controller,
    int worker_count) {
  if (!controller) {
    return;
  }
  memset(controller, 0, sizeof(*controller));
  if (worker_count < 1) {
    worker_count = 1;
  }
  controller->worker_count = worker_count;
  controller->owner_target = worker_count;
}

index_shard_admission_update_t index_shard_admission_record(
    index_shard_admission_controller_t *controller,
    const index_shard_resource_sample_t *sample) {
  index_shard_admission_update_t update;
  int owner_delta = 0;

  memset(&update, 0, sizeof(update));
  if (!controller || !sample ||
      controller->worker_count < 1 ||
      !sample->completed ||
      sample->cancelled ||
      sample->failed ||
      sample->user_seconds < 0.0 ||
      sample->system_seconds < 0.0 ||
      sample->wall_seconds < 0.0) {
    return update;
  }

  controller->samples++;
  update.read_amplified =
      index_shard_sample_read_amplified(sample);
  update.fault_stalled =
      index_shard_sample_fault_stalled(sample);
  update.io_bound =
      update.read_amplified ||
      update.fault_stalled ||
      (sample->input_blocks >=
           INDEX_SHARD_MIN_SAMPLE_BLOCKS &&
       sample->system_seconds >
           sample->user_seconds);
  update.cpu_bound =
      !update.io_bound &&
      sample->user_seconds >
          sample->system_seconds * 2.0;

  if (update.fault_stalled) {
    /*
     * When NORMAL's own feedback has converged to exact page faults, more
     * independent owners can overlap storage latency. Reducing ownership in
     * this state would compound the stall.
     */
    controller->io_samples++;
    owner_delta = -2;
  } else if (update.io_bound) {
    /*
     * NORMAL read amplification benefits from fewer simultaneous readaround
     * streams and more compute lending into already-open index epochs.
     */
    controller->io_samples++;
    owner_delta = 2;
  } else if (update.cpu_bound) {
    controller->cpu_samples++;
    owner_delta = -1;
  } else if (controller->owner_pressure_score > 0) {
    owner_delta = -1;
  } else if (controller->owner_pressure_score < 0) {
    owner_delta = 1;
  }

  controller->owner_pressure_score =
      index_shard_saturating_add(
          controller->owner_pressure_score,
          owner_delta,
          INDEX_SHARD_PRESSURE_LIMIT);
  if (controller->owner_pressure_score >=
          INDEX_SHARD_PRESSURE_LIMIT) {
    if (controller->owner_target > 1) {
      controller->owner_target--;
      controller->owner_target_decreases++;
      update.owner_target_changed = 1;
    }
    controller->owner_pressure_score = 0;
  } else if (controller->owner_pressure_score <=
                 -INDEX_SHARD_PRESSURE_LIMIT) {
    if (controller->owner_target <
        controller->worker_count) {
      controller->owner_target++;
      controller->owner_target_increases++;
      update.owner_target_changed = 1;
    }
    controller->owner_pressure_score = 0;
  }

  return update;
}

int index_shard_config_available_cpus(void) {
  long online_cpus = -1;

#ifdef __linux__
  {
    cpu_set_t affinity;

    CPU_ZERO(&affinity);
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
      int affinity_cpus = CPU_COUNT(&affinity);

      if (affinity_cpus > 0) {
        return affinity_cpus;
      }
    }
  }
#endif

#ifdef _SC_NPROCESSORS_ONLN
  online_cpus = sysconf(_SC_NPROCESSORS_ONLN);
#endif

  if (online_cpus < 1) {
    return 1;
  }

  if (online_cpus > INT_MAX) {
    return INT_MAX;
  }

  return (int)online_cpus;
}

int index_shard_config_validate_workers(int requested_workers,
                                        int available_cpus) {
  if (available_cpus < 1) {
    return -1;
  }

  if (requested_workers == INDEX_SHARD_WORKERS_AUTO) {
    return 0;
  }

  if (requested_workers < 1 || requested_workers > available_cpus) {
    return -1;
  }

  return 0;
}

int index_shard_config_parse_workers(const char *value,
                                     int available_cpus,
                                     int *requested_workers) {
  char *end = NULL;
  long parsed;

  if (!value || !requested_workers || available_cpus < 1) {
    return -1;
  }

  if (!strcmp(value, "auto")) {
    *requested_workers = INDEX_SHARD_WORKERS_AUTO;
    return 0;
  }

  /*
   * Reject signs, leading zeroes, whitespace, and partially parsed values.
   * This keeps 0 and every negative spelling unambiguously invalid.
   */
  if (value[0] < '1' || value[0] > '9') {
    return -1;
  }

  errno = 0;
  parsed = strtol(value, &end, 10);

  if (errno == ERANGE ||
      end == value ||
      *end != '\0' ||
      parsed > INT_MAX) {
    return -1;
  }

  if (index_shard_config_validate_workers((int)parsed, available_cpus)) {
    return -1;
  }

  *requested_workers = (int)parsed;
  return 0;
}

int index_shard_config_resolve_workers(int requested_workers,
                                       int available_cpus) {
  if (index_shard_config_validate_workers(requested_workers,
                                          available_cpus)) {
    return -1;
  }

  if (requested_workers == INDEX_SHARD_WORKERS_AUTO) {
    return available_cpus;
  }

  return requested_workers;
}

int index_shard_config_effective_workers(int configured_workers,
                                         size_t nindexes) {
  int workers = configured_workers;

  if (workers < 1) {
    workers = 1;
  }

  (void)nindexes;

  return workers;
}

int index_shard_config_phase_assist_eligible(int worker_count,
                                             size_t nindexes) {
  return worker_count > 1 &&
      nindexes > 0U &&
      nindexes < (size_t)worker_count;
}

int index_shard_config_build_phase_groups(
    int worker_count,
    size_t nindexes,
    index_shard_phase_group_config_t *groups,
    size_t group_capacity,
    size_t *worker_group,
    int *worker_local_rank) {
  int group_count;
  int workers_per_group;
  int groups_with_extra;
  int next_worker = 0;
  int group_index;

  if (!index_shard_config_phase_assist_eligible(
          worker_count,
          nindexes) ||
      !groups ||
      group_capacity < nindexes ||
      !worker_group ||
      !worker_local_rank) {
    return -1;
  }

  group_count = (int)nindexes;
  workers_per_group = worker_count / group_count;
  groups_with_extra = worker_count % group_count;

  for (group_index = 0;
       group_index < group_count;
       group_index++) {
    index_shard_phase_group_config_t *group =
        &groups[group_index];
    int local_rank;

    group->index_order = (size_t)group_index;
    group->first_worker = next_worker;
    group->worker_count =
        workers_per_group +
        (group_index < groups_with_extra ? 1 : 0);
    if (group->worker_count < 1 ||
        group->first_worker < 0 ||
        group->first_worker >
            worker_count - group->worker_count) {
      return -1;
    }

    for (local_rank = 0;
         local_rank < group->worker_count;
         local_rank++) {
      int worker_id = group->first_worker + local_rank;

      worker_group[worker_id] = (size_t)group_index;
      worker_local_rank[worker_id] = local_rank;
    }
    next_worker += group->worker_count;
  }

  return next_worker == worker_count ? 0 : -1;
}
