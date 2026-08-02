#include <stdint.h>
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

static void check_width_plan(int workers,
                             int io_width,
                             int detached,
                             int exact_demand,
                             size_t expected_producers,
                             size_t expected_helpers) {
  index_shard_width_plan_t plan = { SIZE_MAX, SIZE_MAX };

  CHECK(index_shard_config_plan_widths(
      workers, io_width, detached, exact_demand, &plan) == 0);
  CHECK(plan.producer_width == expected_producers);
  CHECK(plan.helper_width == expected_helpers);
  CHECK(plan.producer_width + plan.helper_width ==
        (size_t)workers);
}

int main(void) {
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

  CHECK(index_shard_config_exact_demand_pass(
      1, 4, 1, 1, 349U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      0, 4, 1, 1, 349U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 0, 1, 1, 349U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 4, 0, 1, 349U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 4, 1, 0, 349U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 4, 1, 1, 0U, 0U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 4, 1, 1, 349U, 1U, 0));
  CHECK(!index_shard_config_exact_demand_pass(
      1, 4, 1, 1, 349U, 0U, 1));
  CHECK(!index_shard_config_exact_demand_pass(
      2, 4, 1, 1, 349U, 0U, 0));

  check_width_plan(1, 0, 0, 0, 1U, 0U);
  check_width_plan(2, 1, 1, 0, 2U, 0U);
  check_width_plan(3, 1, 1, 0, 3U, 0U);
  check_width_plan(4, 2, 1, 0, 4U, 0U);
  check_width_plan(8, 4, 1, 0, 8U, 0U);
  check_width_plan(2, 2, 1, 1, 2U, 0U);
  check_width_plan(4, 4, 1, 1, 4U, 0U);
  check_width_plan(5, 4, 1, 1, 4U, 1U);
  check_width_plan(6, 4, 1, 1, 4U, 2U);
  check_width_plan(8, 4, 1, 1, 4U, 4U);
  check_width_plan(8, 6, 1, 1, 6U, 2U);
  check_width_plan(8, 8, 1, 1, 8U, 0U);
  check_width_plan(4, 0, 0, 0, 3U, 1U);
  {
    index_shard_width_plan_t invalid_plan = { 9U, 9U };

    CHECK(index_shard_config_plan_widths(
        0, 1, 1, 0, &invalid_plan) != 0);
    CHECK(index_shard_config_plan_widths(
        4, -1, 1, 0, &invalid_plan) != 0);
    CHECK(index_shard_config_plan_widths(
        4, 1, 1, 0, NULL) != 0);
    CHECK(index_shard_config_plan_widths(
        4, 0, 1, 1, &invalid_plan) != 0);
    CHECK(index_shard_config_plan_widths(
        4, 2, 0, 1, &invalid_plan) != 0);
    CHECK(index_shard_config_plan_widths(
        4, 2, 1, 2, &invalid_plan) != 0);
    CHECK(index_shard_config_plan_widths(
        4, 2, 2, 0, &invalid_plan) != 0);
  }

  CHECK(detected_available >= 1);
  if (expected_available) {
    CHECK(detected_available == atoi(expected_available));
  }

  if (failures) {
    fprintf(stderr, "%i index-shard configuration test(s) failed\n", failures);
    return 1;
  }

  printf("PASS: index-shard worker configuration\n");
  return 0;
}
