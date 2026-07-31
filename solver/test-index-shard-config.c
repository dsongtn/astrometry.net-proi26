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
