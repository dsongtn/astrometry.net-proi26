#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef __linux__
#include <sched.h>
#endif

#include "index_shard_config.h"

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
