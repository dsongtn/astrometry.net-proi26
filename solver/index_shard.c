/*
 * SECTION INDEX-SHARD: module-overview
 *
 * pthread index-sharding for onefield_run()
 *
 * This module executes one candidate index as one outer shard task.  An outer
 * owner may also publish bounded immutable helper tasks.  It does not split
 * the image, xylist, field stars, quads, or verification math.
 *
 * Ownership model:
 *   - worker threads compute local shard results
 *   - reducer thread merges results into the master onefield_t
 *   - master onefield_t remains the final source of truth
 *
 * Threading model:
 *   - one persistent worker pool per engine job
 *   - one submitted pass per onefield_run() call
 *   - one outer task = one candidate index
 *   - helper tasks carry only bounded, disposable byte ranges
 *   - no pthread_cancel
 *   - stop is cooperative through shared flags + solver.quit_now
 *
 * Safety constraints:
 *   - no shared solver_t between workers
 *   - no worker writes directly into master bp->solutions
 *   - no persistent full index_t cache in production path
 *   - index load/release follows the original onefield ownership hooks
 */
#include "index_shard_private.h"

index_shard_solve_status_t
index_shard_solve(onefield_t *bp,
                  solver_t *base_sp,
                  size_t nindexes,
                  const index_shard_hooks_t *hooks) {
  return index_shard_solve_impl(
      bp,
      base_sp,
      nindexes,
      hooks);
}
