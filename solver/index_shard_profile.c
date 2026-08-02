/*
 * Private implementation module for the index-shard subsystem.
 * See index_shard_private.h for ownership and lock-order invariants.
 */
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "index_shard_private.h"
#include "astrometry/bl.h"
#include "astrometry/errors.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/fitsbin.h"
#include "astrometry/fitsioutils.h"
/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/*
 * Pass, task, and phase profiling snapshots and aggregate counters.
 *
 * This module owns immutable snapshots and aggregate observability helpers.
 */


// ANCHOR INDEX-SHARD: pass-state-snapshot
/*
 * Take one synchronized snapshot of pass termination state.
 *
 * state_mutex is the only lock that protects these fields.  Callers may hold
 * queue_mutex or result_mutex while taking this snapshot; no function may hold
 * state_mutex and then acquire either of those locks.
 */
void index_shard_pass_state_snapshot(index_shard_thread_state_t *shared,
                                            index_shard_pass_state_snapshot_t *snapshot) {
  assert(shared);
  assert(snapshot);

  pthread_mutex_lock(&shared->state_mutex);

  snapshot->stop_requested = shared->stop_requested;
  snapshot->fatal_error = shared->fatal_error;
  snapshot->winner_selected = shared->winner_selected;
  snapshot->solved_published = shared->solved_published;
  snapshot->master_committed = shared->master_committed;
  snapshot->terminal_cause = shared->terminal_cause;
  snapshot->selected_index_order = shared->selected_index_order;
  snapshot->selected_candidate_sequence =
      shared->selected_candidate_sequence;
  snapshot->task_local_failures =
      shared->task_local_failures;
  snapshot->global_integrity_failures =
      shared->global_integrity_failures;
  snapshot->first_stop_wall_since_pass =
      shared->first_stop_wall_since_pass;

  pthread_mutex_unlock(&shared->state_mutex);
}

/*
 * Snapshot completed-pass metrics.
 *
 * This function is called only after index_shard_pool_reduce_first_valid() has
 * returned and all participating workers have left the pass.
 */
double index_shard_timeval_delta_seconds(
    const struct timeval *finish,
    const struct timeval *start) {
  double seconds;

  assert(finish);
  assert(start);

  seconds =
      (double)(finish->tv_sec - start->tv_sec) +
      ((double)(finish->tv_usec - start->tv_usec) / 1000000.0);

  if (seconds < 0.0) {
    return 0.0;
  }

  return seconds;
}

long index_shard_nonnegative_long_delta(long finish,
                                                long start) {
  if (finish < start) {
    return 0;
  }

  return finish - start;
}

/*
 * Snapshot completed-pass timing and process resource usage.
 *
 * The rusage counters cover the complete process while the pthread pass is
 * active. They are suitable for pass-level attribution but deliberately not
 * treated as exclusive per-task measurements.
 */
void index_shard_pass_metrics_snapshot(
    index_shard_thread_state_t *shared,
    index_shard_pass_metrics_snapshot_t *snapshot) {
  struct rusage finish;

  assert(shared);
  assert(snapshot);

  memset(snapshot, 0, sizeof(*snapshot));

  snapshot->reduced = shared->results_reduced;
  snapshot->wall_seconds =
      monotonic_seconds() - shared->pass_wall_start;
  snapshot->cpu_seconds =
      get_cpu_usage() - shared->pass_cpu_start;

  if (snapshot->wall_seconds > 0.0) {
    snapshot->cpu_percent =
        (100.0 * (double)snapshot->cpu_seconds) /
        snapshot->wall_seconds;
  }

  if (!shared->pass_rusage_valid ||
      getrusage(RUSAGE_SELF, &finish)) {
    return;
  }

  snapshot->resource_available = TRUE;

  snapshot->user_seconds =
      index_shard_timeval_delta_seconds(
          &finish.ru_utime,
          &shared->pass_rusage_start.ru_utime);

  snapshot->system_seconds =
      index_shard_timeval_delta_seconds(
          &finish.ru_stime,
          &shared->pass_rusage_start.ru_stime);

  snapshot->minor_faults =
      index_shard_nonnegative_long_delta(
          finish.ru_minflt,
          shared->pass_rusage_start.ru_minflt);

  snapshot->major_faults =
      index_shard_nonnegative_long_delta(
          finish.ru_majflt,
          shared->pass_rusage_start.ru_majflt);

  snapshot->voluntary_context_switches =
      index_shard_nonnegative_long_delta(
          finish.ru_nvcsw,
          shared->pass_rusage_start.ru_nvcsw);

  snapshot->involuntary_context_switches =
      index_shard_nonnegative_long_delta(
          finish.ru_nivcsw,
          shared->pass_rusage_start.ru_nivcsw);

  snapshot->filesystem_input_blocks =
      index_shard_nonnegative_long_delta(
          finish.ru_inblock,
          shared->pass_rusage_start.ru_inblock);

  snapshot->filesystem_output_blocks =
      index_shard_nonnegative_long_delta(
          finish.ru_oublock,
          shared->pass_rusage_start.ru_oublock);
}
/*
 * Compare task durations for percentile calculation.
 */
static int index_shard_compare_double(const void *left,
                                      const void *right) {
  const double lhs = *(const double *)left;
  const double rhs = *(const double *)right;

  if (lhs < rhs) {
    return -1;
  }

  if (lhs > rhs) {
    return 1;
  }

  return 0;
}

/*
 * Return the zero-based nearest-rank percentile index.
 *
 * The calculation avoids multiplying the full sample count by the percentile,
 * so it remains safe for large size_t values.
 */
static size_t index_shard_percentile_index(size_t count,
                                           unsigned int percentile) {
  size_t rank;
  size_t whole;
  size_t remainder;

  assert(count > 0);
  assert(percentile >= 1);
  assert(percentile <= 100);

  whole = (count / 100) * percentile;
  remainder = (count % 100) * percentile;
  rank = whole + ((remainder + 99) / 100);

  if (rank == 0) {
    rank = 1;
  }

  if (rank > count) {
    rank = count;
  }

  return rank - 1;
}

/*
 * Build a completed-pass profile from immutable worker result slots.
 *
 * This runs only after index_shard_pool_reduce_first_valid() has waited for every
 * participating worker to leave the pass. No task modifies result storage at
 * this point.
 *
 * The terminal serial tail is the interval after every other measured task
 * has finished while the latest-finishing task is still running. If that task
 * started after the second-latest completion, its own start time is used as
 * the lower bound.
 */
void index_shard_task_profile_snapshot(
    const index_shard_result_t *results,
    size_t nresults,
    double pool_wall_seconds,
    index_shard_task_profile_snapshot_t *snapshot) {
  double *durations = NULL;
  size_t sample_count = 0;
  size_t i;

  anbool have_max = FALSE;
  anbool have_latest = FALSE;
  anbool have_second_latest = FALSE;

  double latest_finish = 0.0;
  double latest_start = 0.0;
  double second_latest_finish = 0.0;

  memset(snapshot, 0, sizeof(*snapshot));

  snapshot->max_worker_id = -1;
  snapshot->tail_worker_id = -1;

  if (!results || !nresults) {
    return;
  }

  if (nresults <= ((size_t)-1) / sizeof(*durations)) {
    durations = malloc(nresults * sizeof(*durations));
  }

  for (i = 0; i < nresults; i++) {
    const index_shard_result_t *result = &results[i];

    if (!result->task_started) {
      continue;
    }

    if (!isfinite(result->task_wall_seconds) ||
        !isfinite(result->task_start_since_pass) ||
        !isfinite(result->task_finish_since_pass) ||
        result->task_wall_seconds < 0.0 ||
        result->task_finish_since_pass < result->task_start_since_pass) {
      continue;
    }

    snapshot->executed++;

    if (durations) {
      durations[sample_count++] = result->task_wall_seconds;
    }

    if (!have_max ||
        result->task_wall_seconds > snapshot->task_max_seconds) {
      have_max = TRUE;

      snapshot->task_max_seconds = result->task_wall_seconds;
      snapshot->max_solve_seconds = result->wall_seconds;
      snapshot->max_index_order = result->index_order;
      snapshot->max_worker_id = result->worker_id;
    }

    if (!have_latest ||
        result->task_finish_since_pass > latest_finish) {
      if (have_latest) {
        second_latest_finish = latest_finish;
        have_second_latest = TRUE;
      }

      have_latest = TRUE;
      latest_finish = result->task_finish_since_pass;
      latest_start = result->task_start_since_pass;

      snapshot->tail_index_order = result->index_order;
      snapshot->tail_worker_id = result->worker_id;
    } else if (!have_second_latest ||
               result->task_finish_since_pass > second_latest_finish) {
      second_latest_finish = result->task_finish_since_pass;
      have_second_latest = TRUE;
    }
  }

  if (durations &&
      sample_count > 0 &&
      sample_count == snapshot->executed) {
    qsort(durations,
          sample_count,
          sizeof(*durations),
          index_shard_compare_double);

    snapshot->quantiles_available = TRUE;

    snapshot->task_p50_seconds =
        durations[index_shard_percentile_index(sample_count, 50)];

    snapshot->task_p90_seconds =
        durations[index_shard_percentile_index(sample_count, 90)];

    snapshot->task_p99_seconds =
        durations[index_shard_percentile_index(sample_count, 99)];

    if (snapshot->task_p50_seconds > 0.0) {
      snapshot->max_to_p50 =
          snapshot->task_max_seconds /
          snapshot->task_p50_seconds;
    }
  }

  if (pool_wall_seconds > 0.0) {
    snapshot->max_pool_percent =
        (100.0 * snapshot->task_max_seconds) /
        pool_wall_seconds;
  }

  if (have_latest) {
    double tail_start = latest_start;

    if (have_second_latest &&
        second_latest_finish > tail_start) {
      tail_start = second_latest_finish;
    }

    if (latest_finish > tail_start) {
      snapshot->serial_tail_seconds =
          latest_finish - tail_start;
    }

    if (pool_wall_seconds > 0.0) {
      snapshot->serial_tail_percent =
          (100.0 * snapshot->serial_tail_seconds) /
          pool_wall_seconds;
    }
  }

  free(durations);
}

/*
 * Attribute completed outer-task wall time to reset, index acquisition,
 * solving, result analysis, and index release.
 *
 * Phase totals are sums of per-task wall durations and may exceed pool wall
 * time because pthread workers execute concurrently. Percentages are therefore
 * relative to summed task wall, not elapsed pass wall.
 */
void index_shard_phase_profile_snapshot(
    const index_shard_result_t *results,
    size_t nresults,
    index_shard_phase_profile_snapshot_t *snapshot) {
  double *acquire_durations = NULL;
  double *solve_durations = NULL;
  size_t sample_count = 0;
  size_t i;
  double measured_total;

  memset(snapshot, 0, sizeof(*snapshot));

  if (!results || !nresults) {
    return;
  }

  if (nresults <= ((size_t)-1) / sizeof(*acquire_durations)) {
    acquire_durations =
        malloc(nresults * sizeof(*acquire_durations));

    solve_durations =
        malloc(nresults * sizeof(*solve_durations));
  }

  if (!acquire_durations || !solve_durations) {
    free(acquire_durations);
    free(solve_durations);

    acquire_durations = NULL;
    solve_durations = NULL;
  }

  for (i = 0; i < nresults; i++) {
    const index_shard_result_t *result = &results[i];

    if (!result->task_started ||
        !isfinite(result->task_wall_seconds) ||
        result->task_wall_seconds < 0.0) {
      continue;
    }

    snapshot->executed++;
    snapshot->task_wall_total += result->task_wall_seconds;

    if (isfinite(result->reset_seconds) &&
        result->reset_seconds >= 0.0) {
      snapshot->reset_total += result->reset_seconds;
    }

    if (isfinite(result->acquire_seconds) &&
        result->acquire_seconds >= 0.0) {
      snapshot->acquire_total += result->acquire_seconds;
    }

    if (isfinite(result->wall_seconds) &&
        result->wall_seconds >= 0.0) {
      snapshot->solve_total += result->wall_seconds;
    }

    if (isfinite(result->analyze_seconds) &&
        result->analyze_seconds >= 0.0) {
      snapshot->analyze_total += result->analyze_seconds;
    }

    if (isfinite(result->release_seconds) &&
        result->release_seconds >= 0.0) {
      snapshot->release_total += result->release_seconds;
    }

    if (acquire_durations &&
        solve_durations &&
        isfinite(result->acquire_seconds) &&
        result->acquire_seconds >= 0.0 &&
        isfinite(result->wall_seconds) &&
        result->wall_seconds >= 0.0) {
      acquire_durations[sample_count] =
          result->acquire_seconds;

      solve_durations[sample_count] =
          result->wall_seconds;

      sample_count++;
    }
  }

  measured_total =
      snapshot->reset_total +
      snapshot->acquire_total +
      snapshot->solve_total +
      snapshot->analyze_total +
      snapshot->release_total;

  if (snapshot->task_wall_total > measured_total) {
    snapshot->other_total =
        snapshot->task_wall_total - measured_total;
  }

  if (snapshot->task_wall_total > 0.0) {
    snapshot->reset_percent =
        (100.0 * snapshot->reset_total) /
        snapshot->task_wall_total;

    snapshot->acquire_percent =
        (100.0 * snapshot->acquire_total) /
        snapshot->task_wall_total;

    snapshot->solve_percent =
        (100.0 * snapshot->solve_total) /
        snapshot->task_wall_total;

    snapshot->analyze_percent =
        (100.0 * snapshot->analyze_total) /
        snapshot->task_wall_total;

    snapshot->release_percent =
        (100.0 * snapshot->release_total) /
        snapshot->task_wall_total;

    snapshot->other_percent =
        (100.0 * snapshot->other_total) /
        snapshot->task_wall_total;
  }

  if (acquire_durations &&
      solve_durations &&
      sample_count > 0 &&
      sample_count == snapshot->executed) {
    qsort(acquire_durations,
          sample_count,
          sizeof(*acquire_durations),
          index_shard_compare_double);

    qsort(solve_durations,
          sample_count,
          sizeof(*solve_durations),
          index_shard_compare_double);

    snapshot->quantiles_available = TRUE;

    snapshot->acquire_p50 =
        acquire_durations[
            index_shard_percentile_index(sample_count, 50)];

    snapshot->acquire_p90 =
        acquire_durations[
            index_shard_percentile_index(sample_count, 90)];

    snapshot->acquire_p99 =
        acquire_durations[
            index_shard_percentile_index(sample_count, 99)];

    snapshot->acquire_max =
        acquire_durations[sample_count - 1];

    snapshot->solve_p50 =
        solve_durations[
            index_shard_percentile_index(sample_count, 50)];

    snapshot->solve_p90 =
        solve_durations[
            index_shard_percentile_index(sample_count, 90)];

    snapshot->solve_p99 =
        solve_durations[
            index_shard_percentile_index(sample_count, 99)];

    snapshot->solve_max =
        solve_durations[sample_count - 1];
  }

  free(acquire_durations);
  free(solve_durations);
}
void index_shard_observability_increment(
    unsigned long long *counter) {
  if (counter && *counter != ULLONG_MAX) {
    (*counter)++;
  }
}

void index_shard_observability_add(
    unsigned long long *counter,
    unsigned long long value) {
  if (!counter || !value || *counter == ULLONG_MAX) {
    return;
  }
  if (ULLONG_MAX - *counter < value) {
    *counter = ULLONG_MAX;
  } else {
    *counter += value;
  }
}
