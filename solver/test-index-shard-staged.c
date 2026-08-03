/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>

int index_shard_test_staged_retire_more(void);

int main(void) {
  int failures = index_shard_test_staged_retire_more();

  if (failures) {
    fprintf(stderr,
            "INDEX_SHARD_STAGED_TEST_FAILED failures=%i\n",
            failures);
    return 1;
  }
  printf("INDEX_SHARD_STAGED_TEST_OK cases=19\n");
  return 0;
}
