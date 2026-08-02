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
 *   A global-integrity failure occurred, or execution failed after
 *   master-visible state may have been mutated. Serial fallback is forbidden.
 *
 * INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE:
 *   One or more task-local executions failed before any worker result was
 *   transferred into master-visible state. Serial fallback is permitted.
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

/*
 * Typed worker-hook outcome.
 *
 * The hook that observes an error owns its scope classification. The shard
 * scheduler must not infer task-local versus global-integrity failure from
 * winner timing or from a generic nonzero return value.
 */
typedef enum index_shard_hook_outcome {
  INDEX_SHARD_HOOK_COMPLETED_UNSOLVED = 0,
  INDEX_SHARD_HOOK_COMPLETED_SOLVED,
  INDEX_SHARD_HOOK_CANCELLED,
  INDEX_SHARD_HOOK_WALL_LIMIT,
  INDEX_SHARD_HOOK_CPU_LIMIT,
  INDEX_SHARD_HOOK_TASK_LOCAL_FAILURE,
  INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE
} index_shard_hook_outcome_t;

typedef struct index_shard_hook_result {
  index_shard_hook_outcome_t outcome;
  int error_code;
} index_shard_hook_result_t;

typedef struct index_shard_hooks {
  index_shard_hook_result_t (*get_index)(
      onefield_t *bp,
      size_t index_order,
      index_t **index_out);

  index_shard_hook_result_t (*done_with_index)(
      onefield_t *bp,
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
   * create_worker_view() freezes one pass-owned immutable view before the
   * generation is published. prepare_local_context() receives only that view,
   * never the reducer-owned onefield_t or its mutable solver. Narrow index
   * acquisition and terminal-limit services remain separate synchronized
   * pass boundaries.
   *
   * prepare_local_context() runs once per worker per submitted pass.
   * reset_local_context_for_task() runs before every one-index solve.
   * cleanup_local_context() runs once when the worker finishes the pass.
   * destroy_worker_view() runs only after every worker has quiesced.
   */
  int (*create_worker_view)(onefield_t *master_bp,
                            const solver_t *base_sp,
                            void **worker_view_out);

  void (*destroy_worker_view)(void *worker_view);

  int (*prepare_local_context)(onefield_t *local_bp,
                               const void *worker_view);

  void (*reset_local_context_for_task)(onefield_t *local_bp,
                                       bl *local_solutions);

  void (*cleanup_local_context)(onefield_t *local_bp);

  index_shard_hook_result_t (*solve_one_index)(
      onefield_t *local_bp,
      index_t *index);

  index_shard_hook_result_t (*analyze_solutions)(
      onefield_t *master_bp,
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
 * range may point to an immutable package arena. No range may point to
 * mutable owner state, an index_t, a solver, a callback, or a reducer. A task
 * may borrow a const read-only payload view when the
 * outer owner retains its source lease until the synchronous group is fully
 * quiescent. execute() must be reentrant and thread-safe and must not mutate
 * global state or retain any task or source pointer after it returns.
 *
 * A staged task owner may call index_shard_helper_run() only from its typed
 * owner callback. This publishes one same-pool synchronous child group while
 * the staged task retains every input lifetime. Recursive helper publication
 * and calls from staged prepare/submit/poll/foreign-execute callbacks remain
 * unavailable.
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

typedef enum index_shard_helper_retire_status {
  INDEX_SHARD_HELPER_RETIRE_ERROR = -1,
  INDEX_SHARD_HELPER_RETIRE_OK = 0,
  INDEX_SHARD_HELPER_RETIRE_STOPPED = 1
} index_shard_helper_retire_status_t;

/*
 * Called only by the outer owner, never by a foreign helper. Completed tasks
 * retire in increasing task_index order while later tasks may still execute.
 */
typedef index_shard_helper_retire_status_t
(*index_shard_helper_retire_fn)(
    const index_shard_helper_task_t *task,
    size_t task_index,
    void *owner_context);

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
 * Ordered continuation form of index_shard_helper_run(). The owner retires a
 * completed canonical prefix before claiming more local work. READY work is
 * never abandoned while the owner waits for a later completion. Once retire()
 * mutates owner state, a later task failure is terminal for that packet group;
 * callers must not replay the group through a native fallback.
 */
index_shard_helper_run_status_t
index_shard_helper_run_ordered(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats);

/*
 * Persistent staged helper work.
 *
 * Unlike the synchronous helper task above, one logical staged task may
 * release its compute claim while a payload ticket is pending. The outer
 * owner still waits for the complete group, retains every borrowed input and
 * mapping until quiescence, and alone retires results in task order.
 *
 * All operation callbacks run without the shard queue mutex. prepare(),
 * submit(), execute(), and owner() may mutate only their task-private output.
 * poll() must not wait: a terminal return must also collect and destroy the
 * task's payload ticket. cancel() requests cancellation without waiting and
 * returns any nonnegative success value or a negative error.
 * No callback may mutate solver_t, reducer state, or final publication state.
 */
typedef enum index_shard_staged_prepare_status {
  INDEX_SHARD_STAGED_PREPARE_ERROR = -2,
  INDEX_SHARD_STAGED_PREPARE_STOPPED = -1,
  INDEX_SHARD_STAGED_PREPARE_MORE = 0,
  INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY = 1,
  INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY = 2,
  INDEX_SHARD_STAGED_PREPARE_OWNER_READY = 3,
  INDEX_SHARD_STAGED_PREPARE_RESULTS_READY = 4
} index_shard_staged_prepare_status_t;

typedef enum index_shard_staged_submit_status {
  INDEX_SHARD_STAGED_SUBMIT_ERROR = -2,
  INDEX_SHARD_STAGED_SUBMIT_STOPPED = -1,
  INDEX_SHARD_STAGED_SUBMIT_RETRY = 0,
  INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED = 1,
  INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY = 2,
  INDEX_SHARD_STAGED_SUBMIT_OWNER_READY = 3
} index_shard_staged_submit_status_t;

typedef enum index_shard_staged_io_status {
  INDEX_SHARD_STAGED_IO_ERROR = -2,
  INDEX_SHARD_STAGED_IO_CANCELLED = -1,
  INDEX_SHARD_STAGED_IO_PENDING = 0,
  INDEX_SHARD_STAGED_IO_READY = 1,
  INDEX_SHARD_STAGED_IO_FAILED = 2
} index_shard_staged_io_status_t;

typedef enum index_shard_staged_execute_status {
  INDEX_SHARD_STAGED_EXECUTE_ERROR = -1,
  INDEX_SHARD_STAGED_EXECUTE_OK = 0,
  INDEX_SHARD_STAGED_EXECUTE_STOPPED = 1,
  INDEX_SHARD_STAGED_EXECUTE_MORE = 2
} index_shard_staged_execute_status_t;

typedef index_shard_staged_prepare_status_t
(*index_shard_staged_prepare_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes);

typedef index_shard_staged_submit_status_t
(*index_shard_staged_submit_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes,
    unsigned long long *completion_id_out);

/*
 * A successful IO_SUBMITTED return must publish one nonzero immutable
 * completion ID through completion_id_out. Every other return leaves it zero.
 */
typedef index_shard_staged_io_status_t
(*index_shard_staged_poll_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes);

typedef int
(*index_shard_staged_cancel_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes);

typedef index_shard_staged_execute_status_t
(*index_shard_staged_execute_fn)(
    const void *input,
    size_t input_bytes,
    void *output,
    size_t output_bytes);

typedef struct index_shard_staged_ops {
  const char *name;
  index_shard_staged_prepare_fn prepare;
  index_shard_staged_submit_fn submit;
  index_shard_staged_poll_fn poll;
  index_shard_staged_cancel_fn cancel;
  index_shard_staged_execute_fn execute;
  index_shard_staged_execute_fn owner;
  /* Safe to collect a notified terminal ticket on its delivery lane. */
  anbool terminal_poll_inline_safe;
} index_shard_staged_ops_t;

typedef struct index_shard_staged_task {
  const void *input;
  size_t input_bytes;
  void *output;
  size_t output_bytes;
  unsigned long long work_units;

  /* Scheduler-owned while index_shard_staged_run_ordered() is active. */
  unsigned char scheduler_state;
  unsigned char cancel_sent;
  unsigned char completion_pending;
  int callback_status;
  unsigned long long scheduler_epoch;
  unsigned long long completion_id;
  double scheduler_submit_seconds;
  double scheduler_ready_seconds;
  double scheduler_result_seconds;
} index_shard_staged_task_t;

typedef enum index_shard_staged_retire_status {
  INDEX_SHARD_STAGED_RETIRE_ERROR = -1,
  INDEX_SHARD_STAGED_RETIRE_OK = 0,
  INDEX_SHARD_STAGED_RETIRE_STOPPED = 1,
  /*
   * The canonical owner retired one bounded prefix. Requeue this same task
   * without advancing the ordered retirement cursor.
   */
  INDEX_SHARD_STAGED_RETIRE_MORE = 2
} index_shard_staged_retire_status_t;

typedef index_shard_staged_retire_status_t
(*index_shard_staged_retire_fn)(
    const index_shard_staged_task_t *task,
    size_t task_index,
    void *owner_context);

typedef struct index_shard_staged_run_stats {
  size_t owner_claims;
  size_t foreign_claims;
  size_t owner_compute_executes;
  size_t foreign_compute_executes;
  size_t max_concurrent_claims;
  size_t max_compute_running;
  size_t io_submitted;
  size_t io_completed;
  size_t max_io_submitted;
  size_t max_compute_ready;
  size_t max_reorder_ready;
  size_t prepare_claims;
  size_t submit_claims;
  size_t poll_claims;
  size_t execute_claims;
  size_t owner_claims_executed;
  double submit_to_ready_seconds;
  double ready_dwell_seconds;
  double execute_seconds;
  double result_to_retire_seconds;
  double retire_seconds;
  unsigned long long owner_work_units;
  unsigned long long foreign_work_units;
} index_shard_staged_run_stats_t;

/*
 * Return the maximum useful logical packet width for the current pool.
 * W1 and a pool without a completion notifier return zero.
 */
size_t index_shard_staged_capacity(void);

/*
 * Return the stable compute-lane width for a staged group. Unlike
 * index_shard_staged_capacity(), this is not a task-storage bound.
 */
size_t index_shard_staged_compute_width(void);

/*
 * Publish one bounded owner-scoped group. The group object is heap-backed and
 * remains published while its outwardly synchronous owner call schedules
 * preparation, submission, completion, execution, and ordered retirement.
 * The caller retains tasks, inputs, outputs, and owner_context until return.
 */
index_shard_helper_run_status_t
index_shard_staged_run_ordered(
    const index_shard_staged_ops_t *ops,
    index_shard_staged_task_t *tasks,
    size_t task_count,
    index_shard_staged_retire_fn retire,
    void *owner_context,
    index_shard_staged_run_stats_t *stats);

index_shard_solve_status_t
index_shard_solve(onefield_t *bp,
                  solver_t *base_sp,
                  size_t nindexes,
                  const index_shard_hooks_t *hooks);

#endif
