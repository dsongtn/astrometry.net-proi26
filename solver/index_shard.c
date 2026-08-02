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
 *
 * Navigation:
 *   - index_shard_private.h defines shared state, locks, and the module map
 *   - index_shard_pass.c owns one generation from submit through quiescence
 *   - index_shard_scheduler.c and worker.c own claims and execution
 *   - index_shard_staged.c and helper.c own inner-package completion
 *   - index_shard_reducer.c is the sole master-publication boundary
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
