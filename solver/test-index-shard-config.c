#include <stdio.h>
#include <stdlib.h>

#include "index_shard_config.h"

static int failures = 0;

#define CHECK(condition)                                                     \
  do {                                                                       \
    if (!(condition)) {                                                      \
      fprintf(stderr,                                                        \
              "FAIL %s:%i: %s\n",                                           \
              __FILE__,                                                      \
              __LINE__,                                                      \
              #condition);                                                   \
      failures++;                                                            \
    }                                                                        \
  } while (0)

static void check_parse_success(const char *text, int expected) {
  int parsed = INDEX_SHARD_WORKERS_UNSET;

  CHECK(index_shard_config_parse_workers(text, 8, &parsed) == 0);
  CHECK(parsed == expected);
}

static void check_parse_failure(const char *text) {
  int parsed = 77;

  CHECK(index_shard_config_parse_workers(text, 8, &parsed) != 0);
  CHECK(parsed == 77);
}

static void check_phase_topology(
    int workers,
    size_t nindexes,
    const int *expected_group_widths) {
  index_shard_phase_group_config_t groups[8];
  size_t worker_group[8];
  int worker_local_rank[8];
  int expected_first_worker = 0;
  int helper_count = 0;
  int worker_id;
  size_t index_order;

  CHECK(workers <= 8);
  CHECK(nindexes <= 8U);
  if (!expected_group_widths) {
    CHECK(!index_shard_config_phase_assist_eligible(
        workers,
        nindexes));
    CHECK(index_shard_config_build_phase_groups(
        workers,
        nindexes,
        groups,
        8U,
        worker_group,
        worker_local_rank) != 0);
    return;
  }

  CHECK(index_shard_config_phase_assist_eligible(
      workers,
      nindexes));
  CHECK(index_shard_config_build_phase_groups(
      workers,
      nindexes,
      groups,
      8U,
      worker_group,
      worker_local_rank) == 0);

  for (index_order = 0U;
       index_order < nindexes;
       index_order++) {
    CHECK(groups[index_order].index_order == index_order);
    CHECK(groups[index_order].first_worker ==
          expected_first_worker);
    CHECK(groups[index_order].worker_count ==
          expected_group_widths[index_order]);
    CHECK(groups[index_order].worker_count >= 1);
    CHECK(worker_group[expected_first_worker] == index_order);
    CHECK(worker_local_rank[expected_first_worker] == 0);
    expected_first_worker += groups[index_order].worker_count;
    helper_count += groups[index_order].worker_count - 1;
  }

  CHECK(expected_first_worker == workers);
  CHECK((int)nindexes + helper_count == workers);
  for (worker_id = 0;
       worker_id < workers;
       worker_id++) {
    size_t group_index = worker_group[worker_id];

    CHECK(group_index < nindexes);
    CHECK(worker_local_rank[worker_id] >= 0);
    CHECK(worker_local_rank[worker_id] <
          groups[group_index].worker_count);
    CHECK(groups[group_index].first_worker +
          worker_local_rank[worker_id] == worker_id);
  }
}

static index_shard_resource_sample_t resource_sample(
    double user_seconds,
    double system_seconds,
    unsigned long long major_faults,
    unsigned long long input_blocks,
    unsigned long long voluntary_switches) {
  index_shard_resource_sample_t sample;

  sample.user_seconds = user_seconds;
  sample.system_seconds = system_seconds;
  sample.wall_seconds = user_seconds + system_seconds;
  sample.major_faults = major_faults;
  sample.input_blocks = input_blocks;
  sample.voluntary_switches = voluntary_switches;
  sample.completed = 1;
  sample.cancelled = 0;
  sample.failed = 0;
  return sample;
}

static void check_admission_controller(void) {
  index_shard_admission_controller_t controller;
  index_shard_admission_update_t update;
  index_shard_resource_sample_t amplified;
  index_shard_resource_sample_t stalled;
  index_shard_resource_sample_t cpu;
  int i;

  amplified = resource_sample(
      0.25,
      1.00,
      64ULL,
      64ULL * 1024ULL,
      70ULL);
  stalled = resource_sample(
      1.00,
      0.50,
      128ULL,
      128ULL * 8ULL,
      120ULL);
  cpu = resource_sample(
      1.00,
      0.10,
      0ULL,
      0ULL,
      1ULL);

  index_shard_admission_init(&controller, 4);
  CHECK(controller.owner_target == 4);

  update = index_shard_admission_record(
      &controller, &amplified);
  CHECK(update.read_amplified);
  CHECK(update.io_bound);
  update = index_shard_admission_record(
      &controller, &amplified);
  CHECK(controller.owner_target == 3);

  for (i = 0; i < 2; i++) {
    update = index_shard_admission_record(
        &controller, &stalled);
  }
  update = index_shard_admission_record(
      &controller, &stalled);
  update = index_shard_admission_record(
      &controller, &stalled);
  CHECK(controller.owner_target == 4);

  for (i = 0; i < 4; i++) {
    (void)index_shard_admission_record(
        &controller, &cpu);
  }
  CHECK(controller.owner_target == 4);

  amplified.completed = 0;
  update = index_shard_admission_record(
      &controller, &amplified);
  CHECK(!update.owner_target_changed);
  CHECK(controller.samples == 10ULL);
}

int main(void) {
  static const int p2_n1[] = {2};
  static const int p4_n1[] = {4};
  static const int p4_n2[] = {2, 2};
  static const int p4_n3[] = {2, 1, 1};
  static const int p5_n2[] = {3, 2};
  static const int p5_n3[] = {2, 2, 1};
  const char *expected_available =
      getenv("TEST_EXPECTED_AVAILABLE_CPUS");
  int detected_available =
      index_shard_config_available_cpus();

  check_parse_success("auto", INDEX_SHARD_WORKERS_AUTO);
  check_parse_success("1", 1);
  check_parse_success("8", 8);

  check_parse_failure("");
  check_parse_failure("0");
  check_parse_failure("-1");
  check_parse_failure("+1");
  check_parse_failure("01");
  check_parse_failure("9");
  check_parse_failure("1x");
  check_parse_failure(" 1");
  check_parse_failure("1 ");
  check_parse_failure("999999999999999999999999999999999999");

  CHECK(index_shard_config_parse_workers(NULL, 8, NULL) != 0);
  CHECK(index_shard_config_parse_workers("1", 0, NULL) != 0);

  CHECK(index_shard_config_validate_workers(INDEX_SHARD_WORKERS_AUTO, 8) == 0);
  CHECK(index_shard_config_validate_workers(1, 8) == 0);
  CHECK(index_shard_config_validate_workers(8, 8) == 0);
  CHECK(index_shard_config_validate_workers(INDEX_SHARD_WORKERS_UNSET, 8) != 0);
  CHECK(index_shard_config_validate_workers(-2, 8) != 0);
  CHECK(index_shard_config_validate_workers(9, 8) != 0);

  CHECK(index_shard_config_resolve_workers(INDEX_SHARD_WORKERS_AUTO, 2) == 2);
  CHECK(index_shard_config_resolve_workers(INDEX_SHARD_WORKERS_AUTO, 8) == 8);
  CHECK(index_shard_config_resolve_workers(4, 8) == 4);
  CHECK(index_shard_config_resolve_workers(9, 8) == -1);

  CHECK(index_shard_config_effective_workers(8, 0) == 8);
  CHECK(index_shard_config_effective_workers(8, 3) == 8);
  CHECK(index_shard_config_effective_workers(1, 9) == 1);
  CHECK(index_shard_config_effective_workers(0, 9) == 1);
  CHECK(detected_available >= 1);
  if (expected_available) {
    CHECK(detected_available == atoi(expected_available));
  }

  check_phase_topology(2, 1U, p2_n1);
  check_phase_topology(4, 1U, p4_n1);
  check_phase_topology(4, 2U, p4_n2);
  check_phase_topology(4, 3U, p4_n3);
  check_phase_topology(5, 2U, p5_n2);
  check_phase_topology(5, 3U, p5_n3);
  check_phase_topology(4, 4U, NULL);
  check_phase_topology(4, 5U, NULL);
  check_phase_topology(1, 1U, NULL);
  check_phase_topology(4, 0U, NULL);
  check_admission_controller();

  if (failures) {
    fprintf(stderr, "%i index-shard configuration test(s) failed\n", failures);
    return 1;
  }

  printf("PASS: index-shard worker configuration\n");
  return 0;
}
