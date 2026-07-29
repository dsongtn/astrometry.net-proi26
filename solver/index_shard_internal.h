#ifndef ASTROMETRY_INDEX_SHARD_INTERNAL_H
#define ASTROMETRY_INDEX_SHARD_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "astrometry/bl.h"
#include "astrometry/index.h"
#include "astrometry/index_shard.h"
/*
 * Terminal status and ownership contract.
 *
 * INDEX_SHARD_SOLVE_LIFECYCLE_CONFLICT:
 *   The requested pass could not acquire the pool lifecycle. No serial
 *   fallback is permitted because a competing or incompatible execution
 *   context may still own the pool.
 *
 * INDEX_SHARD_SOLVE_TERMINAL_FAILURE:
 *   Execution failed after master-visible state may have been mutated.
 *   Serial fallback is forbidden.
 *
 * INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE:
 *   Execution failed before any worker result was transferred into
 *   master-visible state. Serial fallback is permitted.
 *
 * INDEX_SHARD_SOLVE_HANDLED:
 *   The parallel pass completed normally, either solved or unsolved.
 *
 * INDEX_SHARD_SOLVE_UNAVAILABLE:
 *   The parallel path was not active or available. The original serial
 *   path may be used.
 *
 * Ownership:
 *
 * - Workers own private result slots during computation.
 * - Workers never commit externally visible solver state.
 * - The reducer is the only writer of master-visible solution state.
 * - Losing worker products are destroyed during result disposal.
 * - Winning products transfer ownership exactly once during reduction.
 * - master_committed is published before any potentially mutating merge.
 * - Once master_committed is true, serial fallback is permanently forbidden
 *   for that pass.
 */
typedef enum index_shard_solve_status {
  INDEX_SHARD_SOLVE_LIFECYCLE_CONFLICT = -3,
  INDEX_SHARD_SOLVE_TERMINAL_FAILURE = -2,
  INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE = -1,
  INDEX_SHARD_SOLVE_HANDLED = 0,
  INDEX_SHARD_SOLVE_UNAVAILABLE = 1
} index_shard_solve_status_t;

typedef struct index_shard_hooks {
  index_t *(*get_index)(onefield_t *bp, size_t index_order);
  int (*done_with_index)(onefield_t *bp,
                         size_t index_order,
                         index_t *index);

  /*
   * Report the one reducer-owned solution only after the pass is quiescent
   * and its terminal status is known to be successful.
   */
  int (*report_committed_solution)(onefield_t *bp,
                                   size_t index_order,
                                   int fieldnum,
                                   double best_logodds);

  /*
   * Worker-local context lifecycle.
   *
   * prepare_local_context() runs once per worker per submitted pass.
   * reset_local_context_for_task() runs before every one-index solve.
   * cleanup_local_context() runs once when the worker finishes the pass.
   */
  int (*prepare_local_context)(onefield_t *local_bp,
                               onefield_t *master_bp,
                               const solver_t *base_sp);

  void (*reset_local_context_for_task)(onefield_t *local_bp,
                                       bl *local_solutions);

  void (*cleanup_local_context)(onefield_t *local_bp);

  int (*solve_one_index)(onefield_t *local_bp, index_t *index);

  anbool (*analyze_solutions)(onefield_t *master_bp,
                              bl *solutions,
                              double *best_logodds,
                              int *best_fieldnum);

  int (*merge_solutions)(onefield_t *master_bp,
                         bl *solutions,
                         anbool *solved_out);

  void (*free_solutions)(bl *solutions);
} index_shard_hooks_t;

anbool index_shard_trace_enabled(void);

/* True only while the calling thread is executing a shard worker task. */
anbool index_shard_worker_context_active(void);

/* Focused unit seam for exact prepared-index ownership transfer. */
int onefield_job_index_cache_test_handoff_state(void);

/*
 * Lock-free cooperative-stop check for hot solver boundaries.
 * This is meaningful only while the calling thread owns a shard task.
 */
anbool index_shard_worker_stop_requested(void);

/*
 * Generic bounded helper work.
 *
 * Input bytes are immutable until index_shard_helper_run() returns. Output
 * ranges are pairwise disjoint, do not alias any input, and are invalid when
 * the group returns TASK_FAILED or STOPPED. Inputs may point into explicitly
 * immutable, index-free package arenas whose lifetime covers the group. No
 * range may point to mutable owner state, an index, a FITS mapping, a solver,
 * a callback, or a reducer. execute() must be reentrant and thread-safe and
 * must not mutate global state or retain either byte-range pointer after it
 * returns.
 */
typedef enum index_shard_helper_run_status {
  INDEX_SHARD_HELPER_FATAL = -3,
  INDEX_SHARD_HELPER_TASK_FAILED = -2,
  INDEX_SHARD_HELPER_STOPPED = -1,
  INDEX_SHARD_HELPER_OK = 0,
  INDEX_SHARD_HELPER_UNAVAILABLE = 1
} index_shard_helper_run_status_t;

typedef enum index_shard_helper_task_status {
  INDEX_SHARD_HELPER_TASK_ERROR = -1,
  INDEX_SHARD_HELPER_TASK_OK = 0,
  INDEX_SHARD_HELPER_TASK_STOPPED = 1
} index_shard_helper_task_status_t;

#define INDEX_SHARD_HELPER_MAX_TASKS 64U

typedef index_shard_helper_task_status_t
(*index_shard_helper_execute_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes);

typedef struct index_shard_helper_ops {
  const char *name;
  index_shard_helper_execute_fn execute;
} index_shard_helper_ops_t;

typedef struct index_shard_helper_task {
  const void *input;
  size_t input_bytes;
  void *output;
  size_t output_bytes;
  unsigned long long work_units;

  /* Scheduler-owned while index_shard_helper_run() is active. */
  unsigned char scheduler_state;
  int execute_status;
} index_shard_helper_task_t;

typedef struct index_shard_helper_run_stats {
  size_t owner_tasks;
  size_t foreign_tasks;
  size_t max_concurrent_tasks;
  unsigned long long owner_work_units;
  unsigned long long foreign_work_units;
} index_shard_helper_run_stats_t;

/*
 * Return an advisory count of workers that are outer-idle. The value may
 * become stale immediately after return and must never be used as a
 * correctness predicate.
 */
size_t index_shard_helper_available_workers(void);

/*
 * Acquire an owner-local preparation permit before constructing an expensive
 * helper package. The returned worker count is an admission snapshot, not a
 * worker reservation. A successful permit is consumed by helper_run() or
 * must be released explicitly on every native-fallback path.
 */
size_t index_shard_helper_prepare_reserve(void);
void index_shard_helper_prepare_cancel(void);

/*
 * Publish one synchronous, fixed helper group from the current outer owner.
 * Task zero is reserved for that owner. The call returns only after every
 * task has completed and the stack-backed group is no longer published.
 * Foreign admission reservations remain pool-accounted until claimed,
 * cancelled, yielded to the owner, or released at group quiescence.
 * UNAVAILABLE publishes nothing and permits the native inline path. STOPPED
 * and TASK_FAILED invalidate every output. FATAL has already requested pool
 * fatal state and must be propagated by the caller.
 */
index_shard_helper_run_status_t
index_shard_helper_run(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_run_stats_t *stats);

/*
 * Publish a speculative worker-local solution. It narrows the claim ceiling
 * without cancelling earlier canonical indexes; only the reducer may commit
 * the winner and stop the pool.
 */
void index_shard_worker_publish_solution_candidate(void);

index_shard_solve_status_t
index_shard_solve(onefield_t *bp,
                  solver_t *base_sp,
                  size_t nindexes,
                  const index_shard_hooks_t *hooks);

/*
 * Dormant compatibility executor retained for focused solver seams. The
 * production bounded-producer scheduler never creates, binds, or joins it.
 */
solver_ab_executor_t* solver_ab_executor_create(int worker_count);
void solver_ab_executor_destroy(solver_ab_executor_t* executor);
int solver_ab_executor_set_lending_callbacks(
    solver_ab_executor_t* executor,
    void (*work_notify)(void*),
    int (*available_lenders)(void*),
    void* opaque);
int solver_ab_executor_bind(solver_ab_executor_t* executor,
                            solver_t* owner,
                            index_t* index);
int solver_ab_executor_try_join(solver_ab_executor_t* executor,
                                int* worker_id);
int solver_ab_executor_run_joined(solver_ab_executor_t* executor,
                                  int worker_id);
int solver_ab_executor_quiesce(solver_ab_executor_t* executor);
int solver_ab_executor_finish(solver_ab_executor_t* executor);
void solver_ab_executor_abort(solver_ab_executor_t* executor);

/*
 * Apply reducer-owned traversal deltas atomically. A signed counter boundary
 * is a deterministic execution failure; no counter is partially updated.
 */
int solver_ab_checked_counter_delta(
    solver_t* solver,
    unsigned long long numtries,
    unsigned long long cxdx,
    unsigned long long meanx);

#endif
