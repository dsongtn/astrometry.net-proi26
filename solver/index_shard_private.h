/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef ASTROMETRY_INDEX_SHARD_PRIVATE_H
#define ASTROMETRY_INDEX_SHARD_PRIVATE_H

#include <pthread.h>
#include <stdint.h>
#include <stddef.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "index_shard_internal.h"
#include "index_shard_config.h"
#include "astrometry/fitsbin.h"

/*
 * Private index-shard ownership contract.
 *
 * Lock order:
 *   global pool -> control -> queue -> result -> state -> limit
 *
 * The inverse-cache mutex is independent. A path holding state_mutex must
 * not acquire queue_mutex or result_mutex. queue_mutex is the predicate
 * mutex for queue_cv and every owner_cv; result_mutex owns result_cv.
 */
/*
 * File map for state transitions:
 *
 *   control.c   - TLS, limits, terminal observation, cooperative stop
 *   scheduler.c - queue predicates and exactly-once claims
 *   worker.c    - outer-owner and claimed-package execution
 *   helper.c    - bounded helper groups and their lifetime leases
 *   staged.c    - asynchronous completion registry and READY promotion
 *   reducer.c   - failure classification, arbitration, and master transfer
 *   pass.c      - generation setup, submission, quiescence, and reduction
 *   pool.c      - persistent thread and synchronization-object lifetime
 *   inverse.c   - separately locked inverse-permutation cache
 *   profile.c   - aggregate observation only
 *
 * Add shared mutable state here only when its owner, protecting mutex,
 * generation behavior, and teardown order are explicit.
 */
/*
 * SECTION INDEX-SHARD: types
 */
/*
 * Failure scope is determined by the originating operation, never by whether
 * a winner happened to publish first. Task-local failures are isolated to one
 * index execution and require an exact retry only when the pass finds no
 * winner. Global-integrity failures invalidate the pass at any time.
 */
typedef enum index_shard_failure_class {
  INDEX_SHARD_FAILURE_NONE = 0,
  INDEX_SHARD_FAILURE_TASK_LOCAL,
  INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY
} index_shard_failure_class_t;

/*
 * The first non-fatal terminal event is the pass linearization point.
 * A later global-integrity failure remains terminal and invalidates an
 * elected-but-uncommitted winner.
 */
typedef enum index_shard_terminal_cause {
  INDEX_SHARD_TERMINAL_NONE = 0,
  INDEX_SHARD_TERMINAL_WINNER,
  INDEX_SHARD_TERMINAL_WALL_LIMIT,
  INDEX_SHARD_TERMINAL_CPU_LIMIT,
  INDEX_SHARD_TERMINAL_CANCELLED,
  INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY
} index_shard_terminal_cause_t;

// ANCHOR INDEX-SHARD: result-state
/*
 * Result produced by exactly one shard task.
 *
 * The worker owns all mutable state until completed[index_order] is published.
 * candidate_ready freezes solutions, solved state, best-match metadata and
 * candidate identity for election before index cleanup. The reducer still
 * waits for full task completion and pass quiescence before transferring the
 * selected MatchObj payload into master bp->solutions.
 *
 * Important:
 *   - solutions is worker-local and immutable after candidate_ready
 *   - merged prevents double-free / double-merge
 *   - solved means this shard contains an accepted verified MatchObj
 */
typedef struct index_shard_result {
  bl *solutions; // worker-local MatchObj list for this index

  /* Fixed before candidate publication; aggregated only after quiescence. */
  solver_profile_t solver_profile;

  int failed; // classified failure, not normal "did not solve"
  int rc;
  index_shard_failure_class_t failure_class;

  anbool solved; // accepted solution detected for this shard

  double best_logodds; // diagnostic + future usefulness hint
  int best_fieldnum;

  /*
   * solve wall covers solve_one_index() only.
   *
   * cpu_seconds is a process-wide CPU delta and is diagnostic only when
   * pthread workers overlap; it is not exclusive CPU time for this task.
   */
  double wall_seconds;
  float cpu_seconds;

   /*
   * Full outer-task timing covers local reset, index acquisition, solving,
   * result analysis, and index release. It excludes queue wait and reducer
   * wait.
   */
  anbool task_started;
  int worker_id;
  double task_wall_seconds;
  double task_start_since_pass;
  double task_finish_since_pass;

  /*
   * Non-overlapping wall-time attribution inside the outer task.
   *
   * solve_seconds is represented by wall_seconds above to preserve the
   * existing result layout and diagnostics.
   */
  double reset_seconds;
  double acquire_seconds;
  double analyze_seconds;
  double release_seconds;

  fitsbin_mmap_advice_t mmap_advice;
  anbool task_resource_valid;
  struct rusage task_resource_start;
  double task_user_seconds;
  double task_system_seconds;
  unsigned long long task_major_faults;
  unsigned long long task_input_blocks;
  unsigned long long task_voluntary_switches;

  anbool hit_total_cpulimit;
  anbool hit_total_timelimit;

  anbool cancelled;

  size_t index_order; // original candidate index order in onefield pass
  anbool candidate_ready;
  size_t candidate_sequence;
  size_t completion_sequence;
  int merged;         // reducer already consumed/transferred this result
} index_shard_result_t;

// ANCHOR INDEX-SHARD: shared-pass-state
/*
 * Shared state for one submitted onefield_run() pass.
 *
 * Lifetime:
 *   - initialized by index_shard_pool_submit()
 *   - read/updated by workers + reducer during one pass
 *   - result arrays are owned by index_shard_solve()
 *
 * Locking:
 *   - queue_mutex protects claim state and owner count
 *   - result_mutex protects completed slots, completion sequence, and active
 *     worker count
 *   - state_mutex protects stop/fatal/selected/committed-solve pass state
 *   - limit_mutex protects process-wide CPU-limit publication
 *
 * Do not store per-worker heavy data here.  Per-worker context belongs in
 * index_shard_worker_context_t.
 */
typedef struct index_shard_pool index_shard_pool_t;
typedef struct index_shard_helper_group index_shard_helper_group_t;
typedef struct index_shard_staged_group index_shard_staged_group_t;

#define INDEX_SHARD_COMPLETION_SLOT_NONE SIZE_MAX

#if INDEX_SHARD_HELPER_MAX_TASKS > 64U
#error "staged runnable masks require at most 64 tasks"
#endif

typedef enum index_shard_completion_entry_state {
  INDEX_SHARD_COMPLETION_ENTRY_FREE = 0,
  INDEX_SHARD_COMPLETION_ENTRY_EARLY,
  INDEX_SHARD_COMPLETION_ENTRY_REGISTERED,
  INDEX_SHARD_COMPLETION_ENTRY_NOTIFIED
} index_shard_completion_entry_state_t;

/*
 * Scheduler-owned completion identity. The payload provider retains only its
 * immutable numeric completion ID; no task, group, solver, index, or mapping
 * pointer crosses into the provider. queue_mutex protects this registry.
 */
typedef struct index_shard_completion_entry {
  unsigned long long completion_id;
  unsigned long pool_generation;
  unsigned long long owner_epoch;
  size_t owner_worker;
  size_t owner_index_order;
  size_t task_index;
  size_t next;
  unsigned char state;
} index_shard_completion_entry_t;

typedef struct index_shard_thread_state {
  index_shard_pool_t *pool;          // persistent pool, never task-owned
  onefield_t *bp;                   // master bp, reducer-owned for writes
  const solver_t *base_sp;          // read-only template for local solvers
  const index_shard_hooks_t *hooks; // bridge back into onefield.c
  const void *worker_view;          // immutable, pass-owned worker snapshot

  size_t nindexes;
  size_t canonical_scan_cursor;
  size_t outer_unclaimed;
  size_t outer_running;
  size_t producer_width;
  size_t helper_width;
  size_t queue_waiters;
  size_t helper_groups_active;
  size_t helper_preparations_active;
  size_t helper_foreign_reservations;
  size_t staged_groups_active;
  size_t staged_tickets_active;
  /* Logical borrows kept live by the synchronous outer index owner. */
  size_t staged_source_leases;
  size_t staged_compute_ready;
  size_t staged_reorder_ready;
  size_t staged_compute_running_global;
  size_t staged_max_compute_running;
  size_t staged_max_compute_running_global;
  unsigned long long staged_completion_epoch;
  anbool staged_submit_backpressure;

  /* O(1) completion routing, bounded by workers times staged-task limit. */
  index_shard_completion_entry_t *completion_entries;
  size_t *completion_buckets;
  size_t completion_entry_capacity;
  size_t completion_bucket_count;
  size_t completion_free_head;
  size_t completion_active;
  size_t staged_submit_callbacks_active;
  anbool completion_registry_error;

  unsigned char *outer_states;

  index_shard_result_t *results;
  unsigned char *completed; // result slot is visible to reducer
  size_t results_reduced;
  size_t next_completion_sequence;

  size_t next_candidate_sequence;
  pthread_mutex_t queue_mutex;
  pthread_cond_t queue_cv;

  pthread_mutex_t result_mutex;
  pthread_cond_t result_cv;

  pthread_mutex_t state_mutex;

  pthread_mutex_t limit_mutex;

  int worker_count;
  int active_workers; // workers still participating in pass

  int stop_requested;   // cooperative stop, no new claims
  int fatal_error;      // hard worker/module failure
  int winner_selected;  // first immutable verified result won the pass
  int solved_published; // reducer committed a valid solved result
  int master_committed;
  index_shard_terminal_cause_t terminal_cause;
  double first_stop_wall_since_pass;

  /*
   * Hot solvers read this with an atomic load through worker TLS. The locked
   * state above remains authoritative; this flag only shortens their unwind.
   */
  int worker_stop_requested;

  /* First-valid selection identity, immutable after publication. */
  size_t selected_index_order;
  size_t selected_candidate_sequence;

  /* Reducer-owned identity of the first and only master solution commit. */
  int have_committed_result;
  size_t committed_index_order;

  int limit_reported; // avoid repeated CPU-limit log spam

  double pass_wall_start;
  float pass_cpu_start;

  fitsbin_mmap_advice_t mmap_advice;
  unsigned int mmap_pass_number;
  unsigned long long mmap_advice_failures;

  struct rusage pass_rusage_start;
  int pass_rusage_valid;

  unsigned long long reducer_work_calls;
  double reducer_work_wall_seconds;
  unsigned long long outer_claims;
  unsigned long long helper_groups_published;
  unsigned long long helper_groups_completed;
  unsigned long long helper_tasks_owner;
  unsigned long long helper_tasks_foreign;
  unsigned long long helper_task_failures;
  unsigned long long helper_owner_wait_calls;
  double helper_owner_wait_seconds;
  unsigned long long staged_groups_published;
  unsigned long long staged_groups_completed;
  unsigned long long staged_tasks_owner;
  unsigned long long staged_tasks_foreign;
  unsigned long long staged_compute_owner;
  unsigned long long staged_compute_foreign;
  unsigned long long staged_ready_before_outer_claims;
  unsigned long long staged_task_failures;
  unsigned long long staged_io_submitted;
  unsigned long long staged_io_completed;
  unsigned long long staged_submit_retries;
  unsigned long long staged_submit_rearms;
  unsigned long long staged_submit_handoffs;
  unsigned long long staged_submit_deferrals;
  unsigned long long staged_owner_wait_calls;
  double staged_owner_wait_seconds;
  size_t staged_max_io_submitted;
  size_t staged_max_compute_ready;
  size_t staged_max_reorder_ready;
  unsigned long long staged_prepare_claims;
  unsigned long long staged_submit_claims;
  unsigned long long staged_poll_claims;
  unsigned long long staged_inline_poll_claims;
  unsigned long long staged_inline_poll_failures;
  unsigned long long staged_execute_claims;
  unsigned long long staged_owner_execute_claims;
  double staged_submit_to_ready_seconds;
  double staged_ready_dwell_seconds;
  double staged_execute_seconds;
  double staged_result_to_retire_seconds;
  double staged_retire_seconds;
  unsigned long long task_local_failures;
  unsigned long long global_integrity_failures;
  unsigned long long late_loser_failures;

  /* Aggregate scheduler observability; updated only while profiling is on. */
  anbool observability_enabled;
  unsigned long long queue_broadcasts;
  unsigned long long queue_signals;
  unsigned long long queue_signals_no_work;
  unsigned long long owner_signals;
  unsigned long long owner_signals_coalesced;
  unsigned long long owner_broadcasts;
  unsigned long long queue_wait_calls;
  double queue_wait_seconds;
  unsigned long long completion_notifications;
  unsigned long long completion_groups_scanned;
  unsigned long long completion_tasks_scanned;
  unsigned long long completion_matches;
  unsigned long long completion_registry_registered;
  unsigned long long completion_registry_removed;
  unsigned long long completion_registry_early;
  unsigned long long completion_registry_misses;
  unsigned long long completion_registry_duplicates;
  unsigned long long completion_registry_invalid;
  unsigned long long selection_scans;
  unsigned long long selection_groups_scanned;
  unsigned long long selection_tasks_scanned;
  unsigned long long selection_misses;
  unsigned long long helper_windows_started;
  double helper_window_seconds;
  unsigned long long verification_helper_groups;
  unsigned long long verification_helper_contexts;
} index_shard_thread_state_t;

typedef enum index_shard_outer_state {
  INDEX_SHARD_OUTER_UNCLAIMED = 0,
  INDEX_SHARD_OUTER_RUNNING = 1,
  INDEX_SHARD_OUTER_FINISHED = 2
} index_shard_outer_state_t;

typedef enum index_shard_helper_task_state {
  INDEX_SHARD_HELPER_TASK_UNUSED = 0,
  INDEX_SHARD_HELPER_TASK_READY = 1,
  INDEX_SHARD_HELPER_TASK_RUNNING = 2,
  INDEX_SHARD_HELPER_TASK_DONE = 3,
  INDEX_SHARD_HELPER_TASK_RETIRING = 4,
  INDEX_SHARD_HELPER_TASK_RETIRED = 5
} index_shard_helper_task_state_t;

/*
 * One bounded synchronous group published from an outer-index owner stack.
 * queue_mutex protects the lane pointer, task states, and every counter.
 */
struct index_shard_helper_group {
  const index_shard_helper_ops_t *ops;
  index_shard_helper_task_t *tasks;
  size_t task_count;
  index_shard_helper_retire_fn retire;
  void *owner_context;

  unsigned long generation;
  unsigned long long owner_epoch;
  int owner_worker;
  size_t owner_index_order;

  size_t next_claim;
  size_t next_retire;
  size_t ready_count;
  size_t running_count;
  size_t completed_count;
  size_t foreign_claims;
  size_t owner_claims;
  size_t foreign_reserve;
  size_t foreign_reservations_outstanding;
  anbool owner_reserve_yielded;
  size_t max_running;
  unsigned long long ready_work;
  unsigned long long foreign_work;
  unsigned long long owner_work;

  anbool task_failed;
  anbool stop_seen;
  anbool internal_error;
  anbool verification_group;
  double helper_window_start;
};

typedef enum index_shard_staged_task_state {
  INDEX_SHARD_STAGED_TASK_UNUSED = 0,
  INDEX_SHARD_STAGED_TASK_PREPARE_READY,
  INDEX_SHARD_STAGED_TASK_PREPARING,
  INDEX_SHARD_STAGED_TASK_SUBMIT_READY,
  INDEX_SHARD_STAGED_TASK_SUBMITTING,
  INDEX_SHARD_STAGED_TASK_IO_SUBMITTED,
  INDEX_SHARD_STAGED_TASK_IO_POLLING,
  INDEX_SHARD_STAGED_TASK_IO_CANCELLING,
  INDEX_SHARD_STAGED_TASK_COMPUTE_READY,
  INDEX_SHARD_STAGED_TASK_EXECUTING,
  INDEX_SHARD_STAGED_TASK_OWNER_READY,
  INDEX_SHARD_STAGED_TASK_OWNER_EXECUTING,
  INDEX_SHARD_STAGED_TASK_RESULTS_READY,
  INDEX_SHARD_STAGED_TASK_RETIRING,
  INDEX_SHARD_STAGED_TASK_RETIRED,
  INDEX_SHARD_STAGED_TASK_STOPPED,
  INDEX_SHARD_STAGED_TASK_FAILED
} index_shard_staged_task_state_t;

/*
 * Heap-backed group whose task, input, output, and owner-context storage is
 * retained by the outwardly synchronous caller. queue_mutex protects every
 * field below. No operation callback runs while that mutex is held.
 */
struct index_shard_staged_group {
  index_shard_pool_t *pool;
  const index_shard_staged_ops_t *ops;
  index_shard_staged_task_t *tasks;
  size_t task_count;
  index_shard_staged_retire_fn retire;
  void *owner_context;

  unsigned long generation;
  unsigned long long owner_epoch;
  int owner_worker;
  size_t owner_index_order;

  size_t next_retire;
  size_t running_count;
  size_t compute_running;
  size_t owner_claims;
  size_t foreign_claims;
  size_t owner_compute_executes;
  size_t foreign_compute_executes;
  size_t max_running;
  size_t max_compute_running;
  size_t io_submitted;
  size_t io_completed;
  size_t compute_ready;
  size_t reorder_ready;
  size_t max_io_submitted;
  size_t max_compute_ready;
  size_t max_reorder_ready;
  size_t prepare_claims;
  size_t submit_claims;
  size_t poll_claims;
  size_t inline_poll_claims;
  size_t execute_claims;
  size_t owner_execute_claims;
  double submit_to_ready_seconds;
  double ready_dwell_seconds;
  double execute_seconds;
  double result_to_retire_seconds;
  double retire_seconds;
  unsigned long long owner_work;
  unsigned long long foreign_work;

  size_t completion_registry_entries;

  /* Exact runnable indexes; task_count is bounded to 64. */
  uint64_t prepare_ready_mask;
  uint64_t submit_ready_mask;
  uint64_t submit_wait_mask;
  uint64_t submit_credit_mask;
  uint64_t io_submitted_mask;
  uint64_t completion_pending_mask;
  uint64_t cancel_sent_mask;
  uint64_t compute_ready_mask;
  uint64_t owner_ready_mask;
  uint64_t results_ready_mask;

  anbool cancelling;
  anbool task_failed;
  anbool stop_seen;
  anbool internal_error;
};

typedef struct index_shard_inverse_cache_entry {
  char *filename;
  dev_t device;
  ino_t inode;
  off_t file_size;
  time_t mtime_seconds;
  long mtime_nanoseconds;
  time_t ctime_seconds;
  long ctime_nanoseconds;
  int ndata;
  int ndim;
  u32 treetype;
  int *inverse_perm;
  size_t bytes;
  unsigned int users;
  unsigned long long last_used_tick;
  struct index_shard_inverse_cache_entry *next;
} index_shard_inverse_cache_entry_t;

typedef struct index_shard_inverse_source {
  const char *filename;
  dev_t device;
  ino_t inode;
  off_t file_size;
  time_t mtime_seconds;
  long mtime_nanoseconds;
  time_t ctime_seconds;
  long ctime_nanoseconds;
  int ndata;
  int ndim;
  u32 treetype;
  size_t bytes;
} index_shard_inverse_source_t;

typedef struct index_shard_inverse_lease {
  index_shard_pool_t *pool;
  index_shard_inverse_cache_entry_t *entry;
  startree_t *starkd;
  index_shard_inverse_source_t source;
  anbool source_valid;
  anbool borrowed;
  anbool initially_empty;
  anbool callbacks_registered;
  anbool active_reserved;
  anbool admission_reserved;
  anbool allocation_completed;
  size_t reserved_bytes;
  unsigned long generation_seen;
} index_shard_inverse_lease_t;

// ANCHOR INDEX-SHARD: worker-context-state
/*
 * Private state for one pthread worker.
 *
 * local_bp is reused across all tasks within one submitted pass.  This avoids
 * repeated xylist open/close + local solver allocation per index.
 *
 * Important:
 *   - local_bp must never publish directly into master bp->solutions
 *   - local_context_generation ties local_bp to the active pool generation
 *   - no persistent full index_t cache here in the production path
 */
typedef struct index_shard_worker_context {
  int worker_id;
  unsigned long generation_seen;
  struct index_shard_pool *pool;

  onefield_t local_bp; // worker-local onefield copy
  int local_context_ready;
  unsigned long local_context_generation;
  double pass_prepare_seconds;
  double pass_cleanup_seconds;
  anbool current_outer_active;
  size_t current_index_order;
  /*
   * Set after this worker finishes an outer index. The next selection may
   * consume one foreign COMPUTE_READY staged packet before another outer
   * claim; it never admits preparation, submission, polling, or owner work.
   */
  anbool ready_before_outer_eligible;
  index_shard_helper_group_t *published_helper_group;
  unsigned long long helper_group_epoch;
  index_shard_staged_group_t *published_staged_group;
  unsigned long long staged_group_epoch;
  /*
   * True only while this worker is running the owner callback for its
   * published staged group. That callback may publish one bounded synchronous
   * helper group for immutable child computation; recursive publication is
   * forbidden.
   */
  anbool staged_owner_callback_active;
  index_shard_staged_group_t *staged_owner_callback_group;
  anbool helper_preparation_active;
  unsigned long helper_preparation_generation;
  size_t helper_preparation_index_order;
  size_t helper_preparation_workers;

  /* queue_mutex is the predicate mutex for this owner-only condition. */
  pthread_cond_t owner_cv;
  anbool owner_cv_ready;
  anbool owner_waiting;
  anbool owner_wake_pending;

} index_shard_worker_context_t;

// ANCHOR INDEX-SHARD: pool-state
/*
 * Persistent worker pool for one engine job.
 *
 * The pool survives across multiple onefield_run() submissions.  Workers sleep
 * between generations and wake when index_shard_pool_submit() increments
 * generation.
 */
typedef enum index_shard_work_selection {
  INDEX_SHARD_WORK_ERROR = -1,
  INDEX_SHARD_WORK_DONE = 0,
  INDEX_SHARD_WORK_OUTER = 1,
  INDEX_SHARD_WORK_HELPER = 2
} index_shard_work_selection_t;

typedef struct index_shard_helper_claim {
  index_shard_helper_group_t *group;
  size_t task_index;
} index_shard_helper_claim_t;

typedef enum index_shard_staged_claim_kind {
  INDEX_SHARD_STAGED_CLAIM_NONE = 0,
  INDEX_SHARD_STAGED_CLAIM_PREPARE,
  INDEX_SHARD_STAGED_CLAIM_SUBMIT,
  INDEX_SHARD_STAGED_CLAIM_IO_POLL,
  INDEX_SHARD_STAGED_CLAIM_IO_CANCEL,
  INDEX_SHARD_STAGED_CLAIM_EXECUTE,
  INDEX_SHARD_STAGED_CLAIM_OWNER
} index_shard_staged_claim_kind_t;

typedef struct index_shard_staged_claim {
  index_shard_staged_group_t *group;
  size_t task_index;
  index_shard_staged_claim_kind_t kind;
  anbool owner_claim;
  anbool completion_inline;
  anbool submit_credit;
  unsigned long long observed_completion_epoch;
} index_shard_staged_claim_t;

typedef enum index_shard_inner_claim_kind {
  INDEX_SHARD_INNER_CLAIM_NONE = 0,
  INDEX_SHARD_INNER_CLAIM_HELPER,
  INDEX_SHARD_INNER_CLAIM_STAGED
} index_shard_inner_claim_kind_t;

typedef struct index_shard_inner_claim {
  index_shard_inner_claim_kind_t kind;
  index_shard_helper_claim_t helper;
  index_shard_staged_claim_t staged;
} index_shard_inner_claim_t;

typedef struct index_shard_pass_state_snapshot {
  int stop_requested;
  int fatal_error;
  int winner_selected;
  int solved_published;
  int master_committed;
  index_shard_terminal_cause_t terminal_cause;
  size_t selected_index_order;
  size_t selected_candidate_sequence;
  unsigned long long task_local_failures;
  unsigned long long global_integrity_failures;
  double first_stop_wall_since_pass;
} index_shard_pass_state_snapshot_t;

typedef struct index_shard_pass_metrics_snapshot {
  size_t reduced;

  double wall_seconds;
  float cpu_seconds;
  double cpu_percent;

  int resource_available;
  double user_seconds;
  double system_seconds;

  long minor_faults;
  long major_faults;
  long voluntary_context_switches;
  long involuntary_context_switches;
  long filesystem_input_blocks;
  long filesystem_output_blocks;
} index_shard_pass_metrics_snapshot_t;

typedef struct index_shard_task_profile_snapshot {
  size_t executed;
  int quantiles_available;

  double task_p50_seconds;
  double task_p90_seconds;
  double task_p99_seconds;
  double task_max_seconds;

  double max_solve_seconds;
  size_t max_index_order;
  int max_worker_id;

  double max_to_p50;
  double max_pool_percent;

  double serial_tail_seconds;
  double serial_tail_percent;
  size_t tail_index_order;
  int tail_worker_id;
} index_shard_task_profile_snapshot_t;

typedef struct index_shard_phase_profile_snapshot {
  size_t executed;
  int quantiles_available;

  double task_wall_total;

  double reset_total;
  double acquire_total;
  double solve_total;
  double analyze_total;
  double release_total;
  double other_total;

  double reset_percent;
  double acquire_percent;
  double solve_percent;
  double analyze_percent;
  double release_percent;
  double other_percent;

  double acquire_p50;
  double acquire_p90;
  double acquire_p99;
  double acquire_max;

  double solve_p50;
  double solve_p90;
  double solve_p99;
  double solve_max;
} index_shard_phase_profile_snapshot_t;

typedef enum index_shard_pool_acquire_status {
  INDEX_SHARD_POOL_ACQUIRE_CONFLICT = -1,
  INDEX_SHARD_POOL_ACQUIRE_OK = 0,
  INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE = 1
} index_shard_pool_acquire_status_t;

typedef enum index_shard_staged_select_class {
  INDEX_SHARD_STAGED_SELECT_COMPUTE = 0,
  INDEX_SHARD_STAGED_SELECT_IO,
  INDEX_SHARD_STAGED_SELECT_SUBMIT,
  INDEX_SHARD_STAGED_SELECT_PREPARE
} index_shard_staged_select_class_t;

struct index_shard_pool {

  onefield_t *owner_bp;
  solver_t *owner_sp;

  int worker_count;
  size_t producer_width;
  size_t helper_width;
  pthread_t *threads;
  index_shard_worker_context_t *contexts;

  pthread_mutex_t control_mutex;
  pthread_cond_t work_cv;

  pthread_mutex_t inverse_cache_mutex;
  index_shard_inverse_cache_entry_t *inverse_cache;
  size_t inverse_cache_budget;
  size_t inverse_cache_bytes;
  size_t inverse_active_bytes;
  unsigned long long inverse_cache_hits;
  unsigned long long inverse_cache_misses;
  unsigned long long inverse_cache_admitted;
  unsigned long long inverse_cache_refused;
  unsigned long long inverse_cache_evicted;
  unsigned long long inverse_cache_overcommit;
  unsigned long long inverse_cache_access_tick;
  size_t inverse_cache_peak_bytes;
  size_t inverse_combined_peak_bytes;

  int shutdown;
  int stopping;
  int pass_active;
  int payload_io_owned;
  int payload_completion_registered;
  int ready_workers;
  int tls_startup_error;
  unsigned long generation; // pass submission counter

  index_shard_thread_state_t shared;

} ;


/* index_shard_pass.c */
index_shard_solve_status_t
index_shard_solve_impl(onefield_t *bp,
                       solver_t *base_sp,
                       size_t nindexes,
                       const index_shard_hooks_t *hooks);

/* index_shard.c */
index_shard_solve_status_t
index_shard_solve(onefield_t *bp,
                  solver_t *base_sp,
                  size_t nindexes,
                  const index_shard_hooks_t *hooks);

/* index_shard_control.c */
int index_shard_tls_ensure(void);
int index_shard_set_tls(index_shard_worker_context_t *ctx);
index_shard_worker_context_t *index_shard_get_tls(void);
anbool index_shard_worker_context_active(void);
anbool index_shard_worker_stop_requested(void);
anbool index_shard_pthread_enabled(const onefield_t *bp);
anbool index_shard_trace_enabled(void);
int index_shard_get_worker_count(const onefield_t *bp,
                                        size_t nindexes);
int index_shard_publish_terminal_locked(
    index_shard_thread_state_t *shared,
    index_shard_terminal_cause_t cause);
void index_shard_publish_worker_stop(
    index_shard_thread_state_t *shared);
void index_shard_request_fatal_stop(index_shard_thread_state_t *shared);
void index_shard_publish_committed_solve(
    index_shard_thread_state_t *shared);
void index_shard_mark_master_committed(
    index_shard_thread_state_t *shared);
int index_shard_master_limit_or_cancel_requested(
    index_shard_thread_state_t *shared,
    anbool *hit_total_cpulimit,
    anbool *hit_total_timelimit,
    anbool *cancelled);
int index_shard_master_stop_requested(index_shard_thread_state_t *shared);
index_shard_terminal_cause_t
index_shard_sample_terminal_locked(
    index_shard_thread_state_t *shared,
    const index_shard_result_t *result,
    double *elapsed,
    int *report);
int index_shard_check_global_limits(index_shard_thread_state_t *shared);
void index_shard_poll_from_callback(onefield_t *bp);

/* index_shard_helper.c */
anbool index_shard_helper_outer_claimable_locked(
    const index_shard_thread_state_t *shared);
int index_shard_inner_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    anbool allow_owner,
    index_shard_inner_claim_t *claim);
int index_shard_helper_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim);
int index_shard_helper_owner_claim_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t *task_index);
anbool index_shard_helper_staged_child_allowed(
    const index_shard_worker_context_t *ctx);
size_t index_shard_helper_available_workers(void);
size_t index_shard_helper_prepare_reserve(void);
void index_shard_helper_prepare_cancel(void);
index_shard_helper_run_status_t
index_shard_helper_run(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_run_stats_t *stats);
index_shard_helper_run_status_t
index_shard_helper_run_ordered(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats);

/* index_shard_inverse.c */
size_t index_shard_inverse_cache_budget(void);
void index_shard_inverse_cache_attach(
    index_shard_worker_context_t *ctx,
    index_t *index,
    index_shard_inverse_lease_t *lease);
void index_shard_inverse_cache_release(
    index_shard_worker_context_t *ctx,
    index_t *index,
    index_shard_inverse_lease_t *lease);
void index_shard_inverse_cache_destroy(
    index_shard_pool_t *pool);

/* index_shard_pool.c */
index_shard_pool_acquire_status_t
index_shard_pool_acquire_pass(onefield_t *bp,
                              solver_t *sp,
                              index_shard_pool_t **pool_out);
void index_shard_pool_release_pass(index_shard_pool_t *pool);
int index_shard_pool_start(onefield_t *bp, solver_t *sp);
void index_shard_pool_stop(onefield_t *bp);
int index_shard_pool_active(onefield_t *bp);

/* index_shard_profile.c */
void index_shard_pass_state_snapshot(index_shard_thread_state_t *shared,
                                            index_shard_pass_state_snapshot_t *snapshot);
double index_shard_timeval_delta_seconds(
    const struct timeval *finish,
    const struct timeval *start);
long index_shard_nonnegative_long_delta(long finish,
                                                long start);
void index_shard_pass_metrics_snapshot(
    index_shard_thread_state_t *shared,
    index_shard_pass_metrics_snapshot_t *snapshot);
void index_shard_task_profile_snapshot(
    const index_shard_result_t *results,
    size_t nresults,
    double pool_wall_seconds,
    index_shard_task_profile_snapshot_t *snapshot);
void index_shard_phase_profile_snapshot(
    const index_shard_result_t *results,
    size_t nresults,
    index_shard_phase_profile_snapshot_t *snapshot);
void index_shard_observability_increment(
    unsigned long long *counter);
void index_shard_observability_add(
    unsigned long long *counter,
    unsigned long long value);

/* index_shard_reducer.c */
void index_shard_result_fail(
    index_shard_result_t *result,
    index_shard_failure_class_t failure_class,
    int rc);
int index_shard_apply_hook_result(
    index_shard_result_t *result,
    index_shard_hook_result_t hook_result,
    int solved_allowed);
void index_shard_result_init(index_shard_result_t *result, size_t index_order);
void index_shard_result_dispose(index_shard_result_t *result,
                                       const index_shard_hooks_t *hooks);
void index_shard_result_finish_task(
    index_shard_result_t *result,
    const index_shard_thread_state_t *shared,
    double task_wall_start);
int index_shard_capture_solution_analysis(
    index_shard_thread_state_t *shared,
    index_shard_result_t *result);
int index_shard_arbitrate_candidate(
    index_shard_thread_state_t *shared,
    size_t index_order);
void index_shard_finish_outer_claim(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t index_order);
int index_shard_pool_reduce_first_valid(index_shard_pool_t *pool);
int index_shard_report_committed_solution(
    onefield_t *bp,
    size_t nindexes,
    const index_shard_thread_state_t *shared,
    const index_shard_hooks_t *hooks,
    const index_shard_result_t *results);

/* index_shard_scheduler.c */
void index_shard_queue_signal_locked(
    index_shard_thread_state_t *shared);
void index_shard_owner_signal_locked(
    index_shard_thread_state_t *shared,
    int owner_worker);
void index_shard_notify_progress_locked(
    index_shard_thread_state_t *shared,
    int owner_worker);
void index_shard_queue_broadcast_locked(
    index_shard_thread_state_t *shared);
void index_shard_wake_pass_waiters(index_shard_thread_state_t *shared);
void index_shard_wake_queue_waiters(
    index_shard_thread_state_t *shared);
const char *index_shard_terminal_cause_name(
    index_shard_terminal_cause_t cause);
int index_shard_claim_outer_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t candidate,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice);
index_shard_work_selection_t
index_shard_select_work(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice,
    index_shard_inner_claim_t *inner_claim);

/* index_shard_staged.c */
unsigned long long index_shard_helper_task_work(
    const index_shard_helper_task_t *task);
int index_shard_helper_add_work(
    unsigned long long *total,
    unsigned long long work);
int index_shard_completion_registry_init(
    index_shard_thread_state_t *shared,
    int worker_count);
void index_shard_completion_registry_destroy(
    index_shard_thread_state_t *shared);
int index_shard_completion_registry_register_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    size_t task_index,
    unsigned long long completion_id,
    anbool *already_notified);
int index_shard_completion_registry_remove_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    size_t task_index,
    unsigned long long completion_id);
int index_shard_staged_set_submit_wait_locked(
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    anbool waiting);
int index_shard_staged_set_completion_pending_locked(
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    anbool pending);
void index_shard_staged_refresh_submit_backpressure_locked(
    index_shard_pool_t *pool);
int index_shard_staged_rearm_one_submit_waiter_locked(
    index_shard_pool_t *pool,
    int *owner_out);
int index_shard_staged_set_state_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    index_shard_staged_task_state_t state);
void index_shard_staged_completion_notify(
    void *opaque,
    unsigned long long completion_id);
int index_shard_staged_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    index_shard_staged_select_class_t select_class,
    anbool allow_owner,
    anbool foreign_only,
    index_shard_staged_claim_t *claim);
int index_shard_staged_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim,
    int callback_status,
    double callback_seconds,
    unsigned long long completion_id);
int index_shard_inner_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_inner_claim_t *claim);
int index_shard_staged_retire_one(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group);
size_t index_shard_staged_capacity(void);
size_t index_shard_staged_compute_width(void);
index_shard_helper_run_status_t
index_shard_staged_run_ordered(
    const index_shard_staged_ops_t *ops,
    index_shard_staged_task_t *tasks,
    size_t task_count,
    index_shard_staged_retire_fn retire,
    void *owner_context,
    index_shard_staged_run_stats_t *stats);

/* index_shard_worker.c */
void index_shard_worker_cleanup_pass(index_shard_worker_context_t *ctx,
                                            index_shard_thread_state_t *shared);
void *index_shard_worker_main(void *userdata);

#endif
