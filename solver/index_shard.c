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

#include "index_shard_internal.h"
#include "index_shard_config.h"
#include "astrometry/bl.h"
#include "astrometry/errors.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/fitsbin.h"
#include "astrometry/fitsioutils.h"

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

typedef struct index_shard_thread_state {
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
  size_t staged_max_compute_running;
  unsigned long long staged_completion_epoch;
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
  unsigned long long staged_task_failures;
  unsigned long long staged_io_submitted;
  unsigned long long staged_io_completed;
  unsigned long long staged_submit_retries;
  unsigned long long staged_owner_wait_calls;
  double staged_owner_wait_seconds;
  size_t staged_max_io_submitted;
  size_t staged_max_compute_ready;
  size_t staged_max_reorder_ready;
  unsigned long long staged_prepare_claims;
  unsigned long long staged_submit_claims;
  unsigned long long staged_poll_claims;
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
  size_t execute_claims;
  size_t owner_execute_claims;
  double submit_to_ready_seconds;
  double ready_dwell_seconds;
  double execute_seconds;
  double result_to_retire_seconds;
  double retire_seconds;
  unsigned long long owner_work;
  unsigned long long foreign_work;

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

} index_shard_worker_context_t;

// ANCHOR INDEX-SHARD: pool-state
/*
 * Persistent worker pool for one engine job.
 *
 * The pool survives across multiple onefield_run() submissions.  Workers sleep
 * between generations and wake when index_shard_pool_submit() increments
 * generation.
 */
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

static index_shard_pool_t *index_shard_global_pool = NULL;

static pthread_mutex_t index_shard_global_pool_mutex = PTHREAD_MUTEX_INITIALIZER;

static pthread_key_t index_shard_tls_key;
static pthread_once_t index_shard_tls_once = PTHREAD_ONCE_INIT;
static int index_shard_tls_key_status = EAGAIN;

static size_t index_shard_inverse_cache_budget(void) {
  struct rlimit address_limit;
  long available_pages;
  long page_size;
  size_t available_bytes;
  size_t budget;

  if (sizeof(void*) < 8U) {
    return 0U;
  }
#if defined(_SC_AVPHYS_PAGES)
  available_pages = sysconf(_SC_AVPHYS_PAGES);
#else
  available_pages = -1;
#endif
  page_size = sysconf(_SC_PAGESIZE);
  if (available_pages <= 0 || page_size <= 0 ||
      (unsigned long)available_pages >
          SIZE_MAX / (unsigned long)page_size) {
    return 0U;
  }
  available_bytes =
      (size_t)available_pages * (size_t)page_size;

  /*
   * Inverse permutations replace compulsory full PERM sweeps in deeper bands,
   * but remain recomputable heap state. Use a bounded share of current memory
   * rather than the old fixed 128 MiB ceiling, which rejected useful entries
   * independently of host size and configured cohort.
   */
  budget = available_bytes / 8U;
#if defined(RLIMIT_AS)
  if (getrlimit(RLIMIT_AS, &address_limit) == 0 &&
      address_limit.rlim_cur != RLIM_INFINITY) {
    uintmax_t finite_limit =
        (uintmax_t)address_limit.rlim_cur;

    finite_limit = MIN(
        finite_limit,
        (uintmax_t)SIZE_MAX);
    budget = MIN(
        budget,
        (size_t)finite_limit / 8U);
  }
#else
  (void)address_limit;
#endif
  return budget;
}

static int index_shard_inverse_source(
    const startree_t *starkd,
    index_shard_inverse_source_t *source) {
  fitsbin_t *fb;
  struct stat source_stat;

  if (!starkd || !starkd->tree || !starkd->tree->io ||
      !starkd->tree->perm || !source) {
    return -1;
  }
  memset(source, 0, sizeof(*source));
  fb = (fitsbin_t*)starkd->tree->io;
  if (fitsbin_get_open_file_stat(
          fb, &source_stat)) {
    return -1;
  }
  source->filename = fitsbin_get_filename(fb);
  source->device = source_stat.st_dev;
  source->inode = source_stat.st_ino;
  source->file_size = source_stat.st_size;
  source->mtime_seconds = source_stat.st_mtime;
  source->ctime_seconds = source_stat.st_ctime;
#if defined(__APPLE__)
  source->mtime_nanoseconds =
      source_stat.st_mtimespec.tv_nsec;
  source->ctime_nanoseconds =
      source_stat.st_ctimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
  source->mtime_nanoseconds = source_stat.st_mtim.tv_nsec;
  source->ctime_nanoseconds = source_stat.st_ctim.tv_nsec;
#endif
  source->ndata = startree_N(starkd);
  source->ndim = starkd->tree->ndim;
  source->treetype = starkd->tree->treetype;
  if (source->ndata <= 0 ||
      (size_t)source->ndata > SIZE_MAX / sizeof(int)) {
    return -1;
  }
  source->bytes = (size_t)source->ndata * sizeof(int);
  return 0;
}

static anbool index_shard_inverse_entry_matches(
    const index_shard_inverse_cache_entry_t *entry,
    const index_shard_inverse_source_t *source) {
  return entry && source &&
      entry->ndata == source->ndata &&
      entry->ndim == source->ndim &&
      entry->treetype == source->treetype &&
      entry->bytes == source->bytes &&
      entry->device == source->device &&
      entry->inode == source->inode &&
      entry->file_size == source->file_size &&
      entry->mtime_seconds == source->mtime_seconds &&
      entry->mtime_nanoseconds == source->mtime_nanoseconds &&
      entry->ctime_seconds == source->ctime_seconds &&
      entry->ctime_nanoseconds == source->ctime_nanoseconds;
}

static anbool index_shard_inverse_source_path_unchanged(
    const index_shard_inverse_source_t *source) {
  struct stat current;
  time_t mtime_seconds;
  time_t ctime_seconds;
  long mtime_nanoseconds = 0L;
  long ctime_nanoseconds = 0L;

  if (!source || !source->filename ||
      stat(source->filename, &current)) {
    return FALSE;
  }
  mtime_seconds = current.st_mtime;
  ctime_seconds = current.st_ctime;
#if defined(__APPLE__)
  mtime_nanoseconds = current.st_mtimespec.tv_nsec;
  ctime_nanoseconds = current.st_ctimespec.tv_nsec;
#elif defined(__linux__) || defined(__FreeBSD__)
  mtime_nanoseconds = current.st_mtim.tv_nsec;
  ctime_nanoseconds = current.st_ctim.tv_nsec;
#endif
  return source->device == current.st_dev &&
      source->inode == current.st_ino &&
      source->file_size == current.st_size &&
      source->mtime_seconds == mtime_seconds &&
      source->mtime_nanoseconds == mtime_nanoseconds &&
      source->ctime_seconds == ctime_seconds &&
      source->ctime_nanoseconds == ctime_nanoseconds;
}

static index_shard_inverse_cache_entry_t*
index_shard_inverse_cache_find(
    index_shard_pool_t *pool,
    const index_shard_inverse_source_t *source) {
  index_shard_inverse_cache_entry_t *entry;

  for (entry = pool->inverse_cache;
       entry;
       entry = entry->next) {
    if (index_shard_inverse_entry_matches(
            entry, source)) {
      return entry;
    }
  }
  return NULL;
}

static void index_shard_inverse_cache_free_entry(
    index_shard_inverse_cache_entry_t *entry) {
  if (!entry) {
    return;
  }
  free(entry->inverse_perm);
  free(entry->filename);
  free(entry);
}

static anbool index_shard_inverse_cache_make_room(
    index_shard_pool_t *pool,
    size_t bytes) {
  size_t evictable = 0U;
  size_t retained_limit;
  index_shard_inverse_cache_entry_t *entry;

  /*
   * If currently active, non-evictable inverses already make this request
   * impossible, leave the retained LRU intact. Evicting hot state cannot
   * create enough room in that case.
   */
  if (bytes > pool->inverse_cache_budget ||
      pool->inverse_active_bytes >
          pool->inverse_cache_budget - bytes) {
    return FALSE;
  }
  retained_limit =
      pool->inverse_cache_budget -
      pool->inverse_active_bytes -
      bytes;
  if (pool->inverse_cache_bytes <= retained_limit) {
    return TRUE;
  }

  /*
   * Preflight reclaimable bytes before unlinking anything. A pinned cache can
   * reject a new admission, but that refusal must not partially destroy the
   * useful unpinned LRU.
   */
  for (entry = pool->inverse_cache;
       entry;
       entry = entry->next) {
    if (!entry->users) {
      if (SIZE_MAX - evictable < entry->bytes) {
        evictable = SIZE_MAX;
        break;
      }
      evictable += entry->bytes;
    }
  }
  if (evictable <
      pool->inverse_cache_bytes - retained_limit) {
    return FALSE;
  }

  while (pool->inverse_cache_bytes > retained_limit) {
    index_shard_inverse_cache_entry_t *previous = NULL;
    index_shard_inverse_cache_entry_t *oldest = NULL;
    index_shard_inverse_cache_entry_t *oldest_previous = NULL;

    for (entry = pool->inverse_cache;
         entry;
         entry = entry->next) {
      if (!entry->users &&
          (!oldest ||
           entry->last_used_tick <
               oldest->last_used_tick)) {
        oldest = entry;
        oldest_previous = previous;
      }
      previous = entry;
    }
    if (!oldest) {
      return FALSE;
    }
    if (oldest_previous) {
      oldest_previous->next = oldest->next;
    } else {
      pool->inverse_cache = oldest->next;
    }
    assert(pool->inverse_cache_bytes >= oldest->bytes);
    pool->inverse_cache_bytes -= oldest->bytes;
    pool->inverse_cache_evicted++;
    logverb("[index-shard] inverse-cache state=evict index=%s "
            "bytes=%zu used=%zu budget=%zu\n",
            oldest->filename ? oldest->filename : "(unnamed)",
            oldest->bytes,
            pool->inverse_cache_bytes,
            pool->inverse_cache_budget);
    index_shard_inverse_cache_free_entry(oldest);
  }
  return TRUE;
}

static void index_shard_inverse_active_release_locked(
    index_shard_pool_t *pool,
    index_shard_inverse_lease_t *lease);

static void index_shard_inverse_prepare_callback(
    void *opaque,
    size_t bytes) {
  index_shard_inverse_lease_t *lease = opaque;
  index_shard_pool_t *pool;
  index_shard_inverse_cache_entry_t *entry;

  if (!lease || !lease->pool || !bytes) {
    return;
  }
  pool = lease->pool;
  pthread_mutex_lock(&pool->inverse_cache_mutex);
  /*
   * The descriptor snapshot identifies the mapping that is open in this
   * worker, but the retained inverse is reusable only while the configured
   * pathname still names the same, unchanged file.  Revalidate immediately
   * before a lazy borrow as well as before admission: release-time checking
   * alone cannot protect an old entry from an in-place replacement that
   * happens between index open and the first StarKD hit.
   */
  entry = NULL;
  if (index_shard_inverse_source_path_unchanged(
          &lease->source)) {
    entry = index_shard_inverse_cache_find(
        pool, &lease->source);
  }
  if (entry &&
      !startree_borrow_inverse_perm(
          lease->starkd,
          entry->inverse_perm,
          entry->ndata)) {
    entry->users++;
    entry->last_used_tick =
        ++pool->inverse_cache_access_tick;
    lease->entry = entry;
    lease->borrowed = TRUE;
    pool->inverse_cache_hits++;
    logverb("[index-shard] inverse-cache state=hit index=%s "
            "bytes=%zu generation=%lu\n",
            lease->source.filename ?
                lease->source.filename : "(unnamed)",
            entry->bytes,
            lease->generation_seen);
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }
  pool->inverse_cache_misses++;
  if (lease->active_reserved ||
      lease->allocation_completed ||
      lease->reserved_bytes) {
    logerr("[index-shard] duplicate inverse allocation prepare\n");
    pool->inverse_cache_overcommit++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }

  lease->reserved_bytes = bytes;
  if (bytes == lease->source.bytes &&
      index_shard_inverse_cache_make_room(pool, bytes)) {
    lease->admission_reserved = TRUE;
  } else {
    pool->inverse_cache_overcommit++;
  }

  if (SIZE_MAX - pool->inverse_active_bytes >= bytes) {
    pool->inverse_active_bytes += bytes;
    lease->active_reserved = TRUE;
    if (SIZE_MAX - pool->inverse_cache_bytes >=
        pool->inverse_active_bytes) {
      pool->inverse_combined_peak_bytes =
          MAX(pool->inverse_combined_peak_bytes,
              pool->inverse_cache_bytes +
                  pool->inverse_active_bytes);
    } else {
      pool->inverse_combined_peak_bytes = SIZE_MAX;
    }
  } else {
    lease->admission_reserved = FALSE;
    pool->inverse_cache_overcommit++;
  }
  pthread_mutex_unlock(&pool->inverse_cache_mutex);
}

static void index_shard_inverse_complete_callback(
    void *opaque,
    size_t bytes,
    anbool allocated) {
  index_shard_inverse_lease_t *lease = opaque;
  index_shard_pool_t *pool;

  if (!lease || !lease->pool) {
    return;
  }
  pool = lease->pool;
  pthread_mutex_lock(&pool->inverse_cache_mutex);
  if (lease->borrowed) {
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }
  if (bytes != lease->reserved_bytes) {
    logerr("[index-shard] inverse allocation byte mismatch "
           "prepared=%zu completed=%zu\n",
           lease->reserved_bytes,
           bytes);
    lease->admission_reserved = FALSE;
    pool->inverse_cache_overcommit++;
  }
  lease->allocation_completed = allocated;
  if (!allocated) {
    index_shard_inverse_active_release_locked(
        pool, lease);
    lease->admission_reserved = FALSE;
    /*
     * StarKD is allowed to retry a failed lazy allocation on a later real
     * hit.  Reset the prepare/completion handshake completely; otherwise the
     * retry is misclassified as a duplicate and the cache policy silently
     * remains disabled for this live index.
     */
    lease->allocation_completed = FALSE;
    lease->reserved_bytes = 0U;
  }
  pthread_mutex_unlock(&pool->inverse_cache_mutex);
}

static void index_shard_inverse_cache_attach(
    index_shard_worker_context_t *ctx,
    index_t *index,
    index_shard_inverse_lease_t *lease) {
  index_shard_pool_t *pool;

  if (!lease) {
    return;
  }
  memset(lease, 0, sizeof(*lease));
  if (!ctx || !ctx->pool || !index || !index->starkd ||
      !index->starkd->tree || !index->starkd->tree->perm ||
      index->starkd->inverse_perm) {
    return;
  }
  pool = ctx->pool;
  /*
   * Loaded multiindex components can share one persistent startree_t. Their
   * inverse is naturally retained by that owner and must never enter this
   * per-ephemeral-handle transfer cache.
   */
  if (!pool->shared.bp ||
      pl_size(pool->shared.bp->indexes) != 0) {
    return;
  }
  if (!pool->inverse_cache_budget ||
      index_shard_inverse_source(
          index->starkd,
          &lease->source)) {
    return;
  }
  lease->source_valid = TRUE;
  lease->initially_empty = TRUE;
  lease->pool = pool;
  lease->starkd = index->starkd;
  lease->generation_seen = ctx->generation_seen;
  if (!startree_set_inverse_perm_callbacks(
          index->starkd,
          index_shard_inverse_prepare_callback,
          index_shard_inverse_complete_callback,
          lease)) {
    lease->callbacks_registered = TRUE;
  } else {
    pthread_mutex_lock(&pool->inverse_cache_mutex);
    pool->inverse_cache_refused++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    lease->source_valid = FALSE;
    lease->pool = NULL;
    lease->starkd = NULL;
  }
}

static void index_shard_inverse_active_release_locked(
    index_shard_pool_t *pool,
    index_shard_inverse_lease_t *lease) {
  if (!pool || !lease || !lease->active_reserved) {
    return;
  }
  if (pool->inverse_active_bytes <
      lease->reserved_bytes) {
    logerr("[index-shard] inverse active-byte underflow\n");
    pool->inverse_active_bytes = 0U;
  } else {
    pool->inverse_active_bytes -=
        lease->reserved_bytes;
  }
  lease->active_reserved = FALSE;
}

static void index_shard_inverse_cache_release(
    index_shard_worker_context_t *ctx,
    index_t *index,
    index_shard_inverse_lease_t *lease) {
  index_shard_pool_t *pool;
  startree_t *starkd;
  index_shard_inverse_cache_entry_t *entry;
  index_shard_inverse_source_t source_now;
  const char *log_name;
  size_t admitted_bytes;
  size_t cache_bytes;

  if (!lease) {
    return;
  }
  pool = lease->pool;
  starkd = lease->starkd;

  if (lease->callbacks_registered && starkd) {
    if (startree_clear_inverse_perm_callbacks(
            starkd, lease)) {
      logerr("[index-shard] failed to clear inverse callbacks\n");
      lease->admission_reserved = FALSE;
    }
    lease->callbacks_registered = FALSE;
  }

  if (!ctx || !ctx->pool || ctx->pool != pool ||
      !pool || !index || !index->starkd ||
      !lease->source_valid ||
      starkd != index->starkd) {
    if (pool) {
      pthread_mutex_lock(&pool->inverse_cache_mutex);
      index_shard_inverse_active_release_locked(
          pool, lease);
      pthread_mutex_unlock(&pool->inverse_cache_mutex);
    }
    return;
  }

  if (lease->borrowed) {
    int *borrowed =
        startree_release_borrowed_inverse_perm(starkd);

    pthread_mutex_lock(&pool->inverse_cache_mutex);
    entry = lease->entry;
    if (!entry || entry->inverse_perm != borrowed ||
        !entry->users) {
      logerr("[index-shard] inverse-cache lease mismatch\n");
    } else {
      entry->users--;
      entry->last_used_tick =
          ++pool->inverse_cache_access_tick;
    }
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }

  if (!lease->initially_empty ||
      !starkd->inverse_perm ||
      !starkd->inverse_perm_owned ||
      !lease->allocation_completed ||
      !lease->active_reserved ||
      !lease->admission_reserved ||
      index_shard_inverse_source(starkd, &source_now) ||
      !index_shard_inverse_source_path_unchanged(
          &lease->source) ||
      !index_shard_inverse_entry_matches(
          &(index_shard_inverse_cache_entry_t){
              .device = lease->source.device,
              .inode = lease->source.inode,
              .file_size = lease->source.file_size,
              .mtime_seconds = lease->source.mtime_seconds,
              .mtime_nanoseconds =
                  lease->source.mtime_nanoseconds,
              .ctime_seconds = lease->source.ctime_seconds,
              .ctime_nanoseconds =
                  lease->source.ctime_nanoseconds,
              .ndata = lease->source.ndata,
              .ndim = lease->source.ndim,
              .treetype = lease->source.treetype,
              .bytes = lease->source.bytes},
          &source_now)) {
    pthread_mutex_lock(&pool->inverse_cache_mutex);
    index_shard_inverse_active_release_locked(
        pool, lease);
    if (lease->allocation_completed) {
      pool->inverse_cache_refused++;
    }
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }

  entry = calloc(1, sizeof(*entry));
  if (!entry) {
    pthread_mutex_lock(&pool->inverse_cache_mutex);
    index_shard_inverse_active_release_locked(
        pool, lease);
    pool->inverse_cache_refused++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }
  entry->filename = strdup(
      source_now.filename ? source_now.filename : "(unnamed)");
  if (!entry->filename) {
    free(entry);
    pthread_mutex_lock(&pool->inverse_cache_mutex);
    index_shard_inverse_active_release_locked(
        pool, lease);
    pool->inverse_cache_refused++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    return;
  }
  entry->device = source_now.device;
  entry->inode = source_now.inode;
  entry->file_size = source_now.file_size;
  entry->mtime_seconds = source_now.mtime_seconds;
  entry->mtime_nanoseconds = source_now.mtime_nanoseconds;
  entry->ctime_seconds = source_now.ctime_seconds;
  entry->ctime_nanoseconds = source_now.ctime_nanoseconds;
  entry->ndata = source_now.ndata;
  entry->ndim = source_now.ndim;
  entry->treetype = source_now.treetype;
  entry->bytes = source_now.bytes;

  pthread_mutex_lock(&pool->inverse_cache_mutex);
  entry->last_used_tick =
      ++pool->inverse_cache_access_tick;
  if (index_shard_inverse_cache_find(
          pool, &source_now) ||
      entry->bytes != lease->reserved_bytes ||
      entry->bytes > pool->inverse_cache_budget) {
    index_shard_inverse_active_release_locked(
        pool, lease);
    pool->inverse_cache_refused++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    index_shard_inverse_cache_free_entry(entry);
    return;
  }
  entry->inverse_perm =
      startree_take_inverse_perm(starkd);
  if (!entry->inverse_perm) {
    index_shard_inverse_active_release_locked(
        pool, lease);
    pool->inverse_cache_refused++;
    pthread_mutex_unlock(&pool->inverse_cache_mutex);
    index_shard_inverse_cache_free_entry(entry);
    return;
  }
  index_shard_inverse_active_release_locked(
      pool, lease);
  entry->next = pool->inverse_cache;
  pool->inverse_cache = entry;
  pool->inverse_cache_bytes += entry->bytes;
  pool->inverse_cache_peak_bytes =
      MAX(pool->inverse_cache_peak_bytes,
          pool->inverse_cache_bytes);
  pool->inverse_cache_admitted++;
  admitted_bytes = entry->bytes;
  cache_bytes = pool->inverse_cache_bytes;
  log_name = index->indexname ? index->indexname :
      source_now.filename;
  pthread_mutex_unlock(&pool->inverse_cache_mutex);

  logverb("[index-shard] inverse-cache state=admit index=%s "
          "bytes=%zu used=%zu budget=%zu generation=%lu\n",
          log_name ? log_name : "(unnamed)",
          admitted_bytes,
          cache_bytes,
          pool->inverse_cache_budget,
          ctx->generation_seen);
}

static void index_shard_inverse_cache_destroy(
    index_shard_pool_t *pool) {
  index_shard_inverse_cache_entry_t *entry;

  if (!pool) {
    return;
  }
  entry = pool->inverse_cache;
  while (entry) {
    index_shard_inverse_cache_entry_t *next = entry->next;

    if (entry->users) {
      logerr("[index-shard] inverse-cache destroyed with %u users\n",
             entry->users);
    }
    index_shard_inverse_cache_free_entry(entry);
    entry = next;
  }
  pool->inverse_cache = NULL;
  pool->inverse_cache_bytes = 0U;
  if (pool->inverse_active_bytes) {
    logerr("[index-shard] inverse-cache destroyed with "
           "%zu active bytes\n",
           pool->inverse_active_bytes);
    pool->inverse_active_bytes = 0U;
  }
}

/*
 * SECTION INDEX-SHARD: tls - thread logical singleton
 *
 * TLS links a running worker callback back to its pool state.
 *
 * solve_fields()/solver_run() call into onefield callbacks.  Those callbacks
 * receive only local onefield_t, so TLS is used to check the global shard stop
 * state and set local_bp->solver.quit_now.
 */
static void index_shard_make_tls_key(void) {
  index_shard_tls_key_status =
      pthread_key_create(&index_shard_tls_key, NULL);
}

static int index_shard_tls_ensure(void) {
  int once_status =
      pthread_once(&index_shard_tls_once,
                   index_shard_make_tls_key);

  if (once_status) {
    return once_status;
  }
  return index_shard_tls_key_status;
}

static int index_shard_set_tls(index_shard_worker_context_t *ctx) {
  int status = index_shard_tls_ensure();

  if (status) {
    return status;
  }
  return pthread_setspecific(index_shard_tls_key, ctx);
}

static index_shard_worker_context_t *index_shard_get_tls(void) {
  if (index_shard_tls_ensure()) {
    return NULL;
  }
  return pthread_getspecific(index_shard_tls_key);
}

anbool index_shard_worker_context_active(void) {
  return index_shard_get_tls() != NULL;
}

anbool index_shard_worker_stop_requested(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();

  if (!ctx || !ctx->pool) {
    return FALSE;
  }

  return __atomic_load_n(
      &ctx->pool->shared.worker_stop_requested,
      __ATOMIC_ACQUIRE) != 0;
}

/*
 * SECTION INDEX-SHARD: configuration
 *
 * The engine resolves config, environment, and per-job overrides before
 * engine_run_job(). Hot execution paths consume the immutable job value only.
 */
anbool index_shard_pthread_enabled(const onefield_t *bp) {
  return bp && bp->index_shard_workers > 1;
}

anbool index_shard_trace_enabled(void) {
  /*
   * Detailed scheduler traces follow the ordinary command-line verbosity
   * model. Two -v flags select LOG_ALL without another environment control.
   */
  return log_get_level() >= LOG_ALL;
}

static int index_shard_get_worker_count(const onefield_t *bp,
                                        size_t nindexes) {
  (void)nindexes;

  if (!bp) {
    return 1;
  }

  return index_shard_config_effective_workers(
      bp->index_shard_workers,
      0);
}

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

// ANCHOR INDEX-SHARD: pass-state-snapshot
/*
 * Take one synchronized snapshot of pass termination state.
 *
 * state_mutex is the only lock that protects these fields.  Callers may hold
 * queue_mutex or result_mutex while taking this snapshot; no function may hold
 * state_mutex and then acquire either of those locks.
 */
static void index_shard_pass_state_snapshot(index_shard_thread_state_t *shared,
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
static double index_shard_timeval_delta_seconds(
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

static long index_shard_nonnegative_long_delta(long finish,
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
static void index_shard_pass_metrics_snapshot(
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
static void index_shard_task_profile_snapshot(
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
static void index_shard_phase_profile_snapshot(
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
// ANCHOR INDEX-SHARD: wake-pass-waiters
static void index_shard_wake_pass_waiters(index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->result_mutex);
  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);
}

static void index_shard_wake_queue_waiters(
    index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->queue_mutex);
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
}

static const char *index_shard_terminal_cause_name(
    index_shard_terminal_cause_t cause) {
  switch (cause) {
  case INDEX_SHARD_TERMINAL_NONE:
    return "none";
  case INDEX_SHARD_TERMINAL_WINNER:
    return "winner";
  case INDEX_SHARD_TERMINAL_WALL_LIMIT:
    return "wall-limit";
  case INDEX_SHARD_TERMINAL_CPU_LIMIT:
    return "cpu-limit";
  case INDEX_SHARD_TERMINAL_CANCELLED:
    return "cancelled";
  case INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY:
    return "global-integrity";
  }
  return "invalid";
}

/* state_mutex must be held. */
static int index_shard_publish_terminal_locked(
    index_shard_thread_state_t *shared,
    index_shard_terminal_cause_t cause) {
  int changed = FALSE;

  assert(shared);
  if (cause == INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  if (cause == INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY) {
    changed = !shared->fatal_error ||
        shared->terminal_cause != cause;
    if (shared->terminal_cause == INDEX_SHARD_TERMINAL_NONE) {
      shared->first_stop_wall_since_pass =
          monotonic_seconds() - shared->pass_wall_start;
    }
    shared->terminal_cause = cause;
    shared->fatal_error = TRUE;
    shared->stop_requested = TRUE;
    return changed;
  }

  if (shared->terminal_cause != INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  shared->terminal_cause = cause;
  shared->first_stop_wall_since_pass =
      monotonic_seconds() - shared->pass_wall_start;
  shared->stop_requested = TRUE;
  return TRUE;
}

static void index_shard_publish_worker_stop(
    index_shard_thread_state_t *shared) {
  __atomic_store_n(
      &shared->worker_stop_requested,
      TRUE,
      __ATOMIC_RELEASE);
  index_shard_wake_pass_waiters(shared);
  index_shard_wake_queue_waiters(shared);
}


// ANCHOR INDEX-SHARD: request-fatal-stop
/*
 * Publish a hard worker/module failure and stop the current pass.
 */
static void index_shard_request_fatal_stop(index_shard_thread_state_t *shared) {
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  (void)index_shard_publish_terminal_locked(
      shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
  stopped = shared->stop_requested;
  pthread_mutex_unlock(&shared->state_mutex);

  if (stopped) {
    index_shard_publish_worker_stop(shared);
  }
}

// ANCHOR INDEX-SHARD: publish-committed-solve
/*
 * Publish that the reducer committed a valid solved result.
 *
 * First-valid selection already stopped new work. This publication records
 * that the reducer transferred the selected result into master state.
 */
static void index_shard_publish_committed_solve(
    index_shard_thread_state_t *shared) {
  int valid;

  pthread_mutex_lock(&shared->state_mutex);
  valid = shared->winner_selected &&
      shared->terminal_cause == INDEX_SHARD_TERMINAL_WINNER &&
      !shared->fatal_error;
  if (valid) {
    shared->solved_published = TRUE;
  } else {
    (void)index_shard_publish_terminal_locked(
        shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
  }
  pthread_mutex_unlock(&shared->state_mutex);

  index_shard_publish_worker_stop(shared);
  if (!valid) {
    logerr("[index-shard] invalid committed-solution terminal state\n");
  }
}
/*
 * Mark the point after which serial fallback is no longer safe.
 *
 * The merge hook is not transactional.  Once it begins transferring a
 * non-empty worker result into master-visible state, a later failure must not
 * cause the caller to rerun the original serial pass.
 */
static void index_shard_mark_master_committed(
    index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->state_mutex);
  shared->master_committed = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);
}

// ANCHOR INDEX-SHARD: master-limit-state
static int index_shard_master_limit_or_cancel_requested(
    index_shard_thread_state_t *shared,
    anbool *hit_total_cpulimit,
    anbool *hit_total_timelimit,
    anbool *cancelled) {
  onefield_t *bp = shared->bp;
  int stop;

  pthread_mutex_lock(&shared->limit_mutex);

  if (hit_total_cpulimit) {
    *hit_total_cpulimit = bp->hit_total_cpulimit;
  }

  if (hit_total_timelimit) {
    *hit_total_timelimit = bp->hit_total_timelimit;
  }

  if (cancelled) {
    *cancelled = bp->cancelled;
  }

  stop = bp->hit_total_cpulimit || bp->hit_total_timelimit || bp->cancelled;

  pthread_mutex_unlock(&shared->limit_mutex);

  return stop;
}

// ANCHOR INDEX-SHARD: master-stop-check
/*
 * Read-only stop predicate used by workers before expensive work.
 *
 * First-valid selection is mirrored through stop_requested. Workers therefore
 * do not read master bp->single_field_solved concurrently with the reducer.
 */
static int index_shard_master_stop_requested(index_shard_thread_state_t *shared) {
  index_shard_pass_state_snapshot_t state;

  index_shard_pass_state_snapshot(shared, &state);

  if (state.stop_requested || state.fatal_error || state.solved_published) {
    return TRUE;
  }

  return index_shard_master_limit_or_cancel_requested(
      shared, NULL, NULL, NULL);
}

/* state_mutex and limit_mutex must be held. */
static index_shard_terminal_cause_t
index_shard_sample_terminal_locked(
    index_shard_thread_state_t *shared,
    const index_shard_result_t *result,
    double *elapsed,
    int *report) {
  onefield_t *bp;
  index_shard_terminal_cause_t cause =
      INDEX_SHARD_TERMINAL_NONE;

  assert(shared);
  bp = shared->bp;
  assert(bp);

  if (elapsed) {
    *elapsed = 0.0;
  }
  if (report) {
    *report = FALSE;
  }

  if (result && result->cancelled) {
    bp->cancelled = TRUE;
  }
  if (result && result->hit_total_timelimit) {
    bp->hit_total_timelimit = TRUE;
  }
  if (result && result->hit_total_cpulimit) {
    bp->hit_total_cpulimit = TRUE;
  }

  if (bp->cancelled) {
    cause = INDEX_SHARD_TERMINAL_CANCELLED;
  } else if (bp->hit_total_timelimit) {
    cause = INDEX_SHARD_TERMINAL_WALL_LIMIT;
  } else if (bp->hit_total_cpulimit) {
    cause = INDEX_SHARD_TERMINAL_CPU_LIMIT;
  } else if (bp->total_timelimit > 0.0) {
    double now = monotonic_seconds();
    double sampled = now - bp->time_total_start;

    if (elapsed) {
      *elapsed = sampled;
    }
    if (now >= 0.0 && sampled >= bp->total_timelimit) {
      bp->hit_total_timelimit = TRUE;
      cause = INDEX_SHARD_TERMINAL_WALL_LIMIT;
    }
  }

  if (cause == INDEX_SHARD_TERMINAL_NONE &&
      bp->total_cpulimit > 0.0) {
    double sampled =
        (double)(get_cpu_usage() - bp->cpu_total_start);

    if (elapsed) {
      *elapsed = sampled;
    }
    if (sampled >= bp->total_cpulimit) {
      bp->hit_total_cpulimit = TRUE;
      cause = INDEX_SHARD_TERMINAL_CPU_LIMIT;
    }
  }

  if ((cause == INDEX_SHARD_TERMINAL_WALL_LIMIT ||
       cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) &&
      !shared->limit_reported) {
    shared->limit_reported = TRUE;
    if (report) {
      *report = TRUE;
    }
  }

  return cause;
}
// ANCHOR INDEX-SHARD: global-limits
/*
 * Process-wide elapsed-time and CPU-budget checks.
 *
 * total_timelimit is one shared monotonic wall-clock deadline and is not
 * divided by worker count. total_cpulimit is aggregate process CPU time; with
 * N active threads it can be consumed roughly N times faster than wall time.
 */
static int index_shard_check_global_limits(index_shard_thread_state_t *shared) {
  onefield_t *bp = shared->bp;
  index_shard_terminal_cause_t cause = INDEX_SHARD_TERMINAL_NONE;
  double elapsed = 0.0;
  int report = FALSE;

  /*
   * state_mutex is the terminal-event arbiter. Holding it while the master
   * limit flags are inspected and updated closes the former interval between
   * deadline publication and cooperative-stop publication.
   */
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->terminal_cause != INDEX_SHARD_TERMINAL_NONE) {
    pthread_mutex_unlock(&shared->state_mutex);
    return TRUE;
  }

  pthread_mutex_lock(&shared->limit_mutex);
  cause = index_shard_sample_terminal_locked(
      shared, NULL, &elapsed, &report);
  if (cause != INDEX_SHARD_TERMINAL_NONE) {
    (void)index_shard_publish_terminal_locked(shared, cause);
  }
  pthread_mutex_unlock(&shared->limit_mutex);
  pthread_mutex_unlock(&shared->state_mutex);

  if (cause == INDEX_SHARD_TERMINAL_NONE) {
    return FALSE;
  }

  index_shard_publish_worker_stop(shared);
  if (report && cause == INDEX_SHARD_TERMINAL_WALL_LIMIT) {
    logmsg("Total wall-clock time limit reached!\n");
    logverb("[index-shard] wall-limit reached total_timelimit=%g "
            "elapsed=%.6f\n",
            bp->total_timelimit,
            elapsed);
  } else if (report && cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) {
    logmsg("Total CPU time limit reached!\n");
    logverb("[index-shard] cpu-budget reached total_cpulimit=%g "
            "elapsed=%.6f\n",
            bp->total_cpulimit,
            elapsed);
  }
  return TRUE;
}

// ANCHOR INDEX-SHARD: callback-poll
/*
 * Called from onefield callbacks/timer paths while solver_run() is active.
 *
 * This is the fast-stop path for workers already inside solver code.
 * It maps global stop -> local solver.quit_now.
 */
void index_shard_poll_from_callback(onefield_t *bp) {
  // non-worker callback, nothing to do
  // local solver should unwind normally through existing quit path
  index_shard_worker_context_t *ctx = index_shard_get_tls();

  if (!ctx || !ctx->pool)
    return;

  if (index_shard_check_global_limits(&ctx->pool->shared)) {
    bp->solver.quit_now = TRUE;
    return;
  }

  if (index_shard_master_stop_requested(&ctx->pool->shared))
    bp->solver.quit_now = TRUE;
}

/*
 * SECTION INDEX-SHARD: result
 */

// ANCHOR INDEX-SHARD: result-init
/*
 * Initialize one result slot before a worker starts solving an index.
 *
 * best_logodds starts at -HUGE_VAL so diagnostics can distinguish "no match"
 * from a real low-confidence match.
 */
static const char *index_shard_failure_class_name(
    index_shard_failure_class_t failure_class) {
  switch (failure_class) {
  case INDEX_SHARD_FAILURE_TASK_LOCAL:
    return "task-local";
  case INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY:
    return "global-integrity";
  case INDEX_SHARD_FAILURE_NONE:
  default:
    return "none";
  }
}

static void index_shard_result_fail(
    index_shard_result_t *result,
    index_shard_failure_class_t failure_class,
    int rc) {
  assert(result);
  assert(failure_class != INDEX_SHARD_FAILURE_NONE);

  result->failed = TRUE;
  result->rc = rc ? rc : -1;
  if (failure_class > result->failure_class) {
    result->failure_class = failure_class;
  }
}

/*
 * Map one typed bridge result into task state.
 *
 * Return 0 for normal completion, 1 for cooperative terminal observation and
 * -1 for a classified failure. A solved outcome is accepted only from the
 * solution-analysis hook.
 */
static int index_shard_apply_hook_result(
    index_shard_result_t *result,
    index_shard_hook_result_t hook_result,
    int solved_allowed) {
  assert(result);

  switch (hook_result.outcome) {
  case INDEX_SHARD_HOOK_COMPLETED_UNSOLVED:
    if (!hook_result.error_code) {
      return 0;
    }
    break;
  case INDEX_SHARD_HOOK_COMPLETED_SOLVED:
    if (solved_allowed && !hook_result.error_code) {
      result->solved = TRUE;
      return 0;
    }
    break;
  case INDEX_SHARD_HOOK_CANCELLED:
    if (!hook_result.error_code) {
      result->cancelled = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_WALL_LIMIT:
    if (!hook_result.error_code) {
      result->hit_total_timelimit = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_CPU_LIMIT:
    if (!hook_result.error_code) {
      result->hit_total_cpulimit = TRUE;
      return 1;
    }
    break;
  case INDEX_SHARD_HOOK_TASK_LOCAL_FAILURE:
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_TASK_LOCAL,
        hook_result.error_code);
    return -1;
  case INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE:
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        hook_result.error_code);
    return -1;
  }

  index_shard_result_fail(
      result,
      INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
      hook_result.error_code);
  return -1;
}

static void index_shard_result_init(index_shard_result_t *result, size_t index_order) {
  memset(result, 0, sizeof(index_shard_result_t));

  result->index_order = index_order;
  result->best_logodds = -HUGE_VAL;
  result->best_fieldnum = -1;

  result->solutions = bl_new(4, sizeof(MatchObj));
}

static void index_shard_result_dispose(index_shard_result_t *result,
                                       const index_shard_hooks_t *hooks) {
  // reducer already transferred ownership to master bp
  if (!result || !result->solutions)
    return;

  if (result->merged) {
    bl_free(result->solutions);
    result->solutions = NULL;
    return;
  }

  if (hooks && hooks->free_solutions) {
    hooks->free_solutions(result->solutions);
    result->solutions = NULL;
    return;
  }

  bl_free(result->solutions);
  result->solutions = NULL;
}

/*
 * Close the full outer-task timing interval.
 *
 * One monotonic clock call is used for both task duration and position within
 * the pass, keeping the hot-path instrumentation minimal.
 */
static void index_shard_result_finish_task(
    index_shard_result_t *result,
    const index_shard_thread_state_t *shared,
    double task_wall_start) {
  double task_wall_finish;
  struct rusage resource_finish;

  assert(result);
  assert(shared);

  task_wall_finish = monotonic_seconds();

  result->task_wall_seconds =
      task_wall_finish - task_wall_start;

  result->task_finish_since_pass =
      task_wall_finish - shared->pass_wall_start;

#if defined(RUSAGE_THREAD)
  if (result->task_resource_valid &&
      getrusage(RUSAGE_THREAD, &resource_finish) == 0) {
    result->task_user_seconds =
        index_shard_timeval_delta_seconds(
            &resource_finish.ru_utime,
            &result->task_resource_start.ru_utime);
    result->task_system_seconds =
        index_shard_timeval_delta_seconds(
            &resource_finish.ru_stime,
            &result->task_resource_start.ru_stime);
    result->task_major_faults =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_majflt,
            result->task_resource_start.ru_majflt);
    result->task_input_blocks =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_inblock,
            result->task_resource_start.ru_inblock);
    result->task_voluntary_switches =
        (unsigned long long)
        index_shard_nonnegative_long_delta(
            resource_finish.ru_nvcsw,
            result->task_resource_start.ru_nvcsw);
  } else {
    result->task_resource_valid = FALSE;
  }
#else
  (void)resource_finish;
  result->task_resource_valid = FALSE;
#endif
}

// ANCHOR INDEX-SHARD: analyze-result
/*
 * Inspect worker-local solutions without merging them.
 *
 * This marks solved/best_logodds before immutable result publication elects
 * a winner. Only the reducer may transfer the selected solution into master
 * state.
 */
static int index_shard_capture_solution_analysis(
    index_shard_thread_state_t *shared,
    index_shard_result_t *result) {
  index_shard_hook_result_t hook_result;

  if (!shared->hooks || !shared->hooks->analyze_solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  hook_result = shared->hooks->analyze_solutions(
      shared->bp, result->solutions,
      &result->best_logodds, &result->best_fieldnum);
  return index_shard_apply_hook_result(
      result, hook_result, TRUE);
}
// ANCHOR INDEX-SHARD: reduce-one-result
static void index_shard_reducer_work_account(
    index_shard_thread_state_t *shared,
    double wall_start) {
  double wall_seconds;

  assert(shared);

  wall_seconds = monotonic_seconds() - wall_start;

  if (wall_seconds < 0.0) {
    wall_seconds = 0.0;
  }

  shared->reducer_work_calls++;
  shared->reducer_work_wall_seconds += wall_seconds;
}

/*
 * Transfer one completed worker result into master onefield state.
 *
 * Only the reducer calls this.  Workers never append directly into
 * master_bp->solutions.
 *
 * The merge hook is not transactional.  Before transferring a non-empty
 * worker result, mark the pass as master-committed so a later failure cannot
 * trigger an unsafe serial rerun.
 */
static int index_shard_reduce_one_result(index_shard_thread_state_t *shared,
                                         index_shard_result_t *result) {
  anbool solved = FALSE;
  double reducer_wall_start;
  int may_mutate_master = FALSE;
  int losing_result = FALSE;

  if (!result || result->merged) {
    return 0;
  }


  pthread_mutex_lock(&shared->state_mutex);
  losing_result =
      shared->winner_selected &&
      result->index_order != shared->selected_index_order;
  pthread_mutex_unlock(&shared->state_mutex);
  if (losing_result) {
    logerr("[index-shard] refusing to merge losing result "
           "index_order=%zu\n",
           result->index_order);
    return -1;
  }

  if (shared->have_committed_result) {
    logerr("[index-shard] refusing a reducer merge after solution commit\n");
    return -1;
  }

  reducer_wall_start = monotonic_seconds();

  if (result->failed || result->rc) {
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  if (!shared->hooks || !shared->hooks->merge_solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  if (result->solutions && bl_size(result->solutions) > 0) {
    may_mutate_master = TRUE;
  }

  if (result->solved) {
    may_mutate_master = TRUE;
  }

  if (may_mutate_master) {
    index_shard_mark_master_committed(shared);
  }

  if (shared->hooks->merge_solutions(shared->bp,
                                     result->solutions,
                                     &solved)) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  result->merged = TRUE;
  shared->results_reduced++;

  if (solved || result->solved) {
    result->solved = TRUE;
    shared->have_committed_result = TRUE;
    shared->committed_index_order = result->index_order;
    index_shard_publish_committed_solve(shared);
  }

  index_shard_reducer_work_account(shared, reducer_wall_start);
  return 0;
}
// ANCHOR INDEX-SHARD: worker-get-index
/*
 * Load one index for one shard task through original onefield hooks.
 *
 * No persistent index_t cache here.  Full-index caching caused unacceptable
 * RSS growth because candidate sets can contain hundreds of heavy indexes.
 */
static index_shard_hook_result_t index_shard_worker_get_index(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared,
    size_t index_order,
    index_t **index_out) {
  index_shard_hook_result_t hook_result = {
      INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE, -1};

  if (index_out) {
    *index_out = NULL;
  }
  if (!shared || !shared->hooks ||
      !shared->hooks->get_index || !index_out) {
    return hook_result;
  }

  hook_result = shared->hooks->get_index(
      shared->bp, index_order, index_out);
  if (hook_result.outcome ==
          INDEX_SHARD_HOOK_COMPLETED_UNSOLVED &&
      !*index_out) {
    hook_result.outcome =
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE;
    hook_result.error_code = -1;
  }

  if (index_shard_trace_enabled() && *index_out) {
    logmsg("[index-shard] worker=%i load index_order=%zu index=%s\n",
           ctx->worker_id, index_order,
           (*index_out)->indexname
               ? (*index_out)->indexname
               : "(null)");
  }

  return hook_result;
}

static int index_shard_apply_index_mmap_advice(
    index_t* index,
    fitsbin_mmap_advice_t advice) {
    fitsbin_t* fb;
    int failures = 0;

    if (!index) {
        return 0;
    }

    /*
     * Apply the component policy to existing mappings. fitsbin keeps compact
     * topology NORMAL and applies the pass policy only to sparse payload.
     */
    if (index->codekd &&
        index->codekd->tree &&
        index->codekd->tree->io) {
        fb = (fitsbin_t*)index->codekd->tree->io;

        if (fitsbin_set_mmap_advice(
                fb,
                advice,
                TRUE)) {
            failures++;
        }
    }

    /*
     * Star KD tree.
     */
    if (index->starkd &&
        index->starkd->tree &&
        index->starkd->tree->io) {
        fb = (fitsbin_t*)index->starkd->tree->io;

        if (fitsbin_set_mmap_advice(
                fb,
                advice,
                TRUE)) {
            failures++;
        }
    }

    /*
     * Quad table.
     */
    if (index->quads &&
        index->quads->fb) {
        if (fitsbin_set_mmap_advice(
                index->quads->fb,
                advice,
                TRUE)) {
            failures++;
        }
    }
    return failures;
}

/*
 * SECTION INDEX-SHARD: queue
 */
// ANCHOR INDEX-SHARD: claim-one
/*
 * Claim one shard task.
 *
 * Invariant:
 *   - each index_order is claimed at most once
 *   - stop/fatal prevents new claims
 *
 * Each participating worker owns at most one synchronous task, so the worker
 * count itself is the concurrency bound; no separate task-credit condition is
 * needed.
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

static unsigned long long index_shard_helper_task_work(
    const index_shard_helper_task_t *task) {
  return task->work_units ? task->work_units : 1U;
}

static int index_shard_helper_add_work(
    unsigned long long *total,
    unsigned long long work) {
  if (ULLONG_MAX - *total < work) {
    return -1;
  }
  *total += work;
  return 0;
}

static unsigned long long index_shard_staged_task_work(
    const index_shard_staged_task_t *task) {
  return task->work_units ? task->work_units : 1U;
}

static anbool index_shard_staged_task_terminal(
    const index_shard_staged_task_t *task) {
  if (!task) {
    return FALSE;
  }
  return task->scheduler_state == INDEX_SHARD_STAGED_TASK_RETIRED ||
      task->scheduler_state == INDEX_SHARD_STAGED_TASK_STOPPED ||
      task->scheduler_state == INDEX_SHARD_STAGED_TASK_FAILED;
}

/* queue_mutex must be held. */
static int index_shard_staged_set_state_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task,
    index_shard_staged_task_state_t state) {
  index_shard_staged_task_state_t previous;
  int invalid = FALSE;

  if (!shared || !group || !task) {
    return -1;
  }
  previous = (index_shard_staged_task_state_t)
      task->scheduler_state;
  if (previous == state) {
    return 0;
  }
  if (previous == INDEX_SHARD_STAGED_TASK_COMPUTE_READY) {
    if (!group->compute_ready || !shared->staged_compute_ready) {
      invalid = TRUE;
    } else {
      group->compute_ready--;
      shared->staged_compute_ready--;
    }
  } else if (previous == INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
    if (!group->reorder_ready || !shared->staged_reorder_ready) {
      invalid = TRUE;
    } else {
      group->reorder_ready--;
      shared->staged_reorder_ready--;
    }
  }

  task->scheduler_state = (unsigned char)state;
  if (state == INDEX_SHARD_STAGED_TASK_COMPUTE_READY) {
    if (group->compute_ready == SIZE_MAX ||
        shared->staged_compute_ready == SIZE_MAX) {
      invalid = TRUE;
    } else {
      group->compute_ready++;
      shared->staged_compute_ready++;
      group->max_compute_ready = MAX(
          group->max_compute_ready, group->compute_ready);
      shared->staged_max_compute_ready = MAX(
          shared->staged_max_compute_ready,
          shared->staged_compute_ready);
    }
  } else if (state == INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
    if (group->reorder_ready == SIZE_MAX ||
        shared->staged_reorder_ready == SIZE_MAX) {
      invalid = TRUE;
    } else {
      group->reorder_ready++;
      shared->staged_reorder_ready++;
      group->max_reorder_ready = MAX(
          group->max_reorder_ready, group->reorder_ready);
      shared->staged_max_reorder_ready = MAX(
          shared->staged_max_reorder_ready,
          shared->staged_reorder_ready);
    }
  }
  if (invalid) {
    group->internal_error = TRUE;
    return -1;
  }
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_staged_group_valid_locked(
    const index_shard_pool_t *pool,
    const index_shard_worker_context_t *owner,
    const index_shard_staged_group_t *group) {
  if (!pool || !owner || !group ||
      group->pool != pool ||
      owner->worker_id < 0 ||
      owner->worker_id >= pool->shared.worker_count ||
      owner->published_staged_group != group ||
      group->owner_worker != owner->worker_id ||
      group->owner_epoch != owner->staged_group_epoch ||
      group->generation != owner->generation_seen ||
      !owner->current_outer_active ||
      group->owner_index_order != owner->current_index_order) {
    return -1;
  }
  return 0;
}

/*
 * Payload completion never enters with the fitsbin mutex held. The immutable
 * completion ID selects exactly one published ticket without retaining a
 * task or group pointer in the I/O service. The shared epoch is only a retry
 * event for submissions that previously met bounded queue pressure.
 */
static void index_shard_staged_completion_notify(
    void *opaque,
    unsigned long long completion_id) {
  index_shard_pool_t *pool = opaque;
  index_shard_thread_state_t *shared;
  index_shard_staged_group_t *matched_group = NULL;
  size_t matches = 0U;
  anbool wake = FALSE;
  int owner;

  if (!pool || !completion_id) {
    return;
  }
  shared = &pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  shared->staged_completion_epoch++;
  if (!shared->staged_completion_epoch) {
    shared->staged_completion_epoch++;
  }
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context =
        &pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;
    size_t task_index;

    if (!group) {
      continue;
    }
    if (index_shard_staged_group_valid_locked(
            pool, context, group)) {
      group->internal_error = TRUE;
      wake = TRUE;
      continue;
    }
    for (task_index = 0U;
         task_index < group->task_count;
         task_index++) {
      index_shard_staged_task_t *task =
          &group->tasks[task_index];

      if (task->scheduler_state ==
          INDEX_SHARD_STAGED_TASK_SUBMIT_READY) {
        wake = TRUE;
      }
      if (task->completion_id != completion_id) {
        continue;
      }
      matches++;
      if (!matched_group) {
        matched_group = group;
      } else if (matched_group != group || matches > 1U) {
        matched_group->internal_error = TRUE;
        group->internal_error = TRUE;
      }
      if (task->scheduler_state !=
              INDEX_SHARD_STAGED_TASK_IO_SUBMITTED &&
          task->scheduler_state !=
              INDEX_SHARD_STAGED_TASK_IO_POLLING &&
          task->scheduler_state !=
              INDEX_SHARD_STAGED_TASK_IO_CANCELLING) {
        group->internal_error = TRUE;
      } else {
        task->completion_pending = TRUE;
      }
      wake = TRUE;
    }
  }
  if (wake) {
    pthread_cond_broadcast(&shared->queue_cv);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
}

typedef enum index_shard_staged_select_class {
  INDEX_SHARD_STAGED_SELECT_COMPUTE = 0,
  INDEX_SHARD_STAGED_SELECT_IO,
  INDEX_SHARD_STAGED_SELECT_SUBMIT,
  INDEX_SHARD_STAGED_SELECT_PREPARE
} index_shard_staged_select_class_t;

/* queue_mutex must be held. */
static int index_shard_staged_task_claimable_locked(
    const index_shard_thread_state_t *shared,
    const index_shard_staged_group_t *group,
    const index_shard_staged_task_t *task,
    index_shard_staged_select_class_t select_class,
    anbool owner_allowed,
    index_shard_staged_claim_kind_t *kind) {
  if (!shared || !group || !task || !kind) {
    return -1;
  }
  *kind = INDEX_SHARD_STAGED_CLAIM_NONE;
  if (select_class == INDEX_SHARD_STAGED_SELECT_IO) {
    if (task->scheduler_state ==
            INDEX_SHARD_STAGED_TASK_IO_SUBMITTED &&
        task->completion_pending) {
      *kind = INDEX_SHARD_STAGED_CLAIM_IO_POLL;
      return 0;
    }
    if (group->cancelling &&
        task->scheduler_state ==
            INDEX_SHARD_STAGED_TASK_IO_SUBMITTED &&
        !task->cancel_sent) {
      *kind = INDEX_SHARD_STAGED_CLAIM_IO_CANCEL;
      return 0;
    }
    return 1;
  }
  if (group->cancelling || group->task_failed ||
      group->internal_error || group->stop_seen) {
    return 1;
  }
  if (select_class == INDEX_SHARD_STAGED_SELECT_COMPUTE) {
    if (task->scheduler_state ==
        INDEX_SHARD_STAGED_TASK_COMPUTE_READY) {
      *kind = INDEX_SHARD_STAGED_CLAIM_EXECUTE;
      return 0;
    }
    if (owner_allowed &&
        task->scheduler_state ==
            INDEX_SHARD_STAGED_TASK_OWNER_READY) {
      *kind = INDEX_SHARD_STAGED_CLAIM_OWNER;
      return 0;
    }
    return 1;
  }
  if (select_class == INDEX_SHARD_STAGED_SELECT_SUBMIT) {
    if (task->scheduler_state ==
            INDEX_SHARD_STAGED_TASK_SUBMIT_READY &&
        task->scheduler_epoch !=
            shared->staged_completion_epoch) {
      *kind = INDEX_SHARD_STAGED_CLAIM_SUBMIT;
      return 0;
    }
    return 1;
  }
  if (select_class == INDEX_SHARD_STAGED_SELECT_PREPARE &&
      task->scheduler_state ==
          INDEX_SHARD_STAGED_TASK_PREPARE_READY) {
    *kind = INDEX_SHARD_STAGED_CLAIM_PREPARE;
    return 0;
  }
  return 1;
}

/* queue_mutex must be held. */
static int index_shard_staged_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    index_shard_staged_select_class_t select_class,
    anbool allow_owner,
    index_shard_staged_claim_t *claim) {
  index_shard_staged_group_t *best_group = NULL;
  size_t best_task = 0U;
  index_shard_staged_claim_kind_t best_kind =
      INDEX_SHARD_STAGED_CLAIM_NONE;
  int owner;

  if (!worker || !worker->pool || !shared || !claim) {
    return -1;
  }
  memset(claim, 0, sizeof(*claim));
  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_worker_context_t *context =
        &worker->pool->contexts[owner];
    index_shard_staged_group_t *group =
        context->published_staged_group;
    size_t task_index;

    if (!group) {
      continue;
    }
    if (index_shard_staged_group_valid_locked(
            worker->pool, context, group)) {
      group->internal_error = TRUE;
      return -1;
    }
    for (task_index = 0U;
         task_index < group->task_count;
         task_index++) {
      index_shard_staged_claim_kind_t kind;
      anbool owner_allowed = allow_owner &&
          owner == worker->worker_id;
      int eligible = index_shard_staged_task_claimable_locked(
          shared,
          group,
          &group->tasks[task_index],
          select_class,
          owner_allowed,
          &kind);

      if (eligible < 0) {
        group->internal_error = TRUE;
        return -1;
      }
      if (eligible) {
        continue;
      }
      if (!best_group ||
          group->owner_index_order <
              best_group->owner_index_order ||
          (group->owner_index_order ==
               best_group->owner_index_order &&
           task_index < best_task)) {
        best_group = group;
        best_task = task_index;
        best_kind = kind;
      }
      break;
    }
  }
  if (!best_group) {
    return 1;
  }

  claim->owner_claim =
      best_group->owner_worker == worker->worker_id;
  if (best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    unsigned long long work = index_shard_staged_task_work(
        &best_group->tasks[best_task]);
    unsigned long long current = claim->owner_claim
        ? best_group->owner_work
        : best_group->foreign_work;

    if (ULLONG_MAX - current < work) {
      best_group->internal_error = TRUE;
      return -1;
    }
  }

  switch (best_kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_PREPARING)) {
      return -1;
    }
    best_group->prepare_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_SUBMITTING)) {
      return -1;
    }
    best_group->submit_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_IO_POLLING)) {
      return -1;
    }
    best_group->tasks[best_task].completion_pending = FALSE;
    best_group->poll_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_IO_CANCELLING)) {
      return -1;
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
    if (best_group->tasks[best_task].scheduler_ready_seconds > 0.0) {
      double now = monotonic_seconds();

      if (now >= best_group->tasks[best_task].scheduler_ready_seconds) {
        best_group->ready_dwell_seconds +=
            now - best_group->tasks[best_task].scheduler_ready_seconds;
      }
      best_group->tasks[best_task].scheduler_ready_seconds = 0.0;
    }
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_EXECUTING)) {
      return -1;
    }
    best_group->execute_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    if (best_group->tasks[best_task].scheduler_ready_seconds > 0.0) {
      double now = monotonic_seconds();

      if (now >= best_group->tasks[best_task].scheduler_ready_seconds) {
        best_group->ready_dwell_seconds +=
            now - best_group->tasks[best_task].scheduler_ready_seconds;
      }
      best_group->tasks[best_task].scheduler_ready_seconds = 0.0;
    }
    if (index_shard_staged_set_state_locked(
        shared, best_group, &best_group->tasks[best_task],
        INDEX_SHARD_STAGED_TASK_OWNER_EXECUTING)) {
      return -1;
    }
    best_group->owner_execute_claims++;
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    best_group->internal_error = TRUE;
    return -1;
  }
  best_group->running_count++;
  if (best_group->running_count > best_group->max_running) {
    best_group->max_running = best_group->running_count;
  }
  claim->group = best_group;
  claim->task_index = best_task;
  claim->kind = best_kind;
  claim->observed_completion_epoch =
      shared->staged_completion_epoch;
  if (claim->owner_claim) {
    best_group->owner_claims++;
    shared->staged_tasks_owner++;
  } else {
    best_group->foreign_claims++;
    shared->staged_tasks_foreign++;
  }
  if (best_kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    unsigned long long work = index_shard_staged_task_work(
        &best_group->tasks[best_task]);

    best_group->compute_running++;
    best_group->max_compute_running = MAX(
        best_group->max_compute_running,
        best_group->compute_running);
    if (claim->owner_claim) {
      best_group->owner_compute_executes++;
      shared->staged_compute_owner++;
      best_group->owner_work += work;
    } else {
      best_group->foreign_compute_executes++;
      shared->staged_compute_foreign++;
      best_group->foreign_work += work;
    }
  }
  return 0;
}

/* queue_mutex must be held. */
static anbool index_shard_helper_outer_claimable_locked(
    const index_shard_thread_state_t *shared) {
  size_t candidate;

  if (!shared || !shared->outer_states ||
      !shared->producer_width ||
      shared->outer_running >= shared->producer_width) {
    return FALSE;
  }
  candidate = shared->canonical_scan_cursor;
  while (candidate < shared->nindexes &&
         shared->outer_states[candidate] !=
             INDEX_SHARD_OUTER_UNCLAIMED) {
    candidate++;
  }
  return candidate < shared->nindexes;
}

/* queue_mutex must be held. */
static size_t index_shard_helper_idle_workers_locked(
    const index_shard_thread_state_t *shared) {
  size_t available = 0U;
  size_t limit;
  size_t spare = 0U;

  if (!shared || shared->worker_count < 2) {
    return 0U;
  }
  limit = (size_t)shared->worker_count - 1U;
  if (index_shard_helper_outer_claimable_locked(shared)) {
    return 0U;
  }
  if ((size_t)shared->worker_count > shared->outer_running) {
    spare = (size_t)shared->worker_count - shared->outer_running;
  }
  available = shared->queue_waiters;
  if (!shared->helper_groups_active && spare > available) {
    available = spare;
  }
  available = MIN(available, limit);
  if (shared->helper_foreign_reservations >= available) {
    return 0U;
  }
  return available - shared->helper_foreign_reservations;
}

/* queue_mutex must be held. */
static int index_shard_helper_release_foreign_reservations_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t count) {
  if (!shared || !group ||
      count > group->foreign_reservations_outstanding ||
      count > shared->helper_foreign_reservations) {
    if (group) {
      group->internal_error = TRUE;
    }
    return -1;
  }
  group->foreign_reservations_outstanding -= count;
  shared->helper_foreign_reservations -= count;
  return 0;
}

/* queue_mutex must be held. */
static void index_shard_helper_cancel_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    index_shard_helper_task_status_t status) {
  size_t i;

  if (!shared || !group || !group->tasks) {
    return;
  }
  (void)index_shard_helper_release_foreign_reservations_locked(
      shared,
      group,
      group->foreign_reservations_outstanding);
  for (i = 0U; i < group->task_count; i++) {
    index_shard_helper_task_t *task = &group->tasks[i];

    if (task->scheduler_state !=
        INDEX_SHARD_HELPER_TASK_READY) {
      continue;
    }
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_DONE;
    task->execute_status = status;
    if (group->ready_count) {
      group->ready_count--;
    } else {
      group->internal_error = TRUE;
    }
    if (group->completed_count < group->task_count) {
      group->completed_count++;
    } else {
      group->internal_error = TRUE;
    }
  }
  if (group->ready_count) {
    group->internal_error = TRUE;
    group->ready_count = 0U;
  }
  group->ready_work = 0U;
  if (status == INDEX_SHARD_HELPER_TASK_ERROR) {
    group->task_failed = TRUE;
  } else {
    group->stop_seen = TRUE;
  }
}

/* queue_mutex must be held. */
static int index_shard_helper_claim_locked(
    index_shard_helper_group_t *group,
    size_t *task_index) {
  index_shard_helper_task_t *task;
  unsigned long long work;
  size_t candidate;

  if (!group || !task_index || !group->tasks) {
    return -1;
  }
  if (!group->ready_count) {
    return 1;
  }

  candidate = group->next_claim;
  while (candidate < group->task_count &&
         group->tasks[candidate].scheduler_state !=
             INDEX_SHARD_HELPER_TASK_READY) {
    candidate++;
  }
  if (candidate >= group->task_count) {
    group->internal_error = TRUE;
    return -1;
  }

  task = &group->tasks[candidate];
  work = index_shard_helper_task_work(task);
  if (group->ready_work < work || !group->ready_count) {
    group->internal_error = TRUE;
    return -1;
  }

  task->scheduler_state = INDEX_SHARD_HELPER_TASK_RUNNING;
  group->next_claim = candidate + 1U;
  group->ready_count--;
  group->running_count++;
  group->ready_work -= work;
  if (group->running_count > group->max_running) {
    group->max_running = group->running_count;
  }
  *task_index = candidate;
  return 0;
}

/* queue_mutex must be held. */
static int index_shard_helper_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    index_shard_helper_claim_t *claim) {
  index_shard_helper_group_t *best = NULL;
  int owner;

  if (!worker || !worker->pool || !shared || !claim) {
    return -1;
  }
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->stop_requested || shared->fatal_error ||
      shared->solved_published) {
    pthread_mutex_unlock(&shared->state_mutex);
    return 1;
  }

  for (owner = 0; owner < shared->worker_count; owner++) {
    index_shard_helper_group_t *candidate =
        worker->pool->contexts[owner].published_helper_group;

    if (!candidate || !candidate->ready_count) {
      continue;
    }
    if (candidate->generation != worker->generation_seen ||
        candidate->owner_worker != owner ||
        candidate->owner_epoch !=
            worker->pool->contexts[owner].helper_group_epoch ||
        candidate->owner_worker == worker->worker_id) {
      candidate->internal_error = TRUE;
      pthread_mutex_unlock(&shared->state_mutex);
      return -1;
    }

    if (!best ||
        candidate->ready_work > best->ready_work ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count > best->ready_count) ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count == best->ready_count &&
         candidate->owner_index_order < best->owner_index_order) ||
        (candidate->ready_work == best->ready_work &&
         candidate->ready_count == best->ready_count &&
         candidate->owner_index_order == best->owner_index_order &&
         candidate->owner_worker < best->owner_worker)) {
      best = candidate;
    }
  }

  if (!best) {
    pthread_mutex_unlock(&shared->state_mutex);
    return 1;
  }

  claim->group = best;
  if (best->foreign_reservations_outstanding &&
      shared->helper_foreign_reservations <
          best->foreign_reservations_outstanding) {
    best->internal_error = TRUE;
    pthread_mutex_unlock(&shared->state_mutex);
    return -1;
  }
  if (index_shard_helper_claim_locked(
          best, &claim->task_index)) {
    pthread_mutex_unlock(&shared->state_mutex);
    return -1;
  }
  if (best->foreign_reservations_outstanding) {
    best->foreign_reservations_outstanding--;
    shared->helper_foreign_reservations--;
  }
  best->foreign_work += index_shard_helper_task_work(
      &best->tasks[claim->task_index]);
  best->foreign_claims++;
  shared->helper_tasks_foreign++;
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->state_mutex);
  return 0;
}

/*
 * Select one inner action while queue_mutex is held. Already prepared
 * computation has priority over completion collection and legacy synchronous
 * packets. New I/O submission and page-plan preparation remain last. No
 * callback is invoked under the queue lock.
 */
static int index_shard_inner_select_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    anbool allow_owner,
    index_shard_inner_claim_t *claim) {
  index_shard_staged_select_class_t select_class;
  int rc;

  if (!worker || !shared || !claim) {
    return -1;
  }
  memset(claim, 0, sizeof(*claim));

  select_class = INDEX_SHARD_STAGED_SELECT_COMPUTE;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_IO;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  rc = index_shard_helper_select_locked(
      worker, shared, &claim->helper);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_HELPER;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_SUBMIT;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }

  select_class = INDEX_SHARD_STAGED_SELECT_PREPARE;
  rc = index_shard_staged_select_locked(
      worker, shared, select_class, allow_owner, &claim->staged);
  if (rc <= 0) {
    if (!rc) {
      claim->kind = INDEX_SHARD_INNER_CLAIM_STAGED;
    }
    return rc;
  }
  return 1;
}

static int index_shard_helper_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim,
    index_shard_helper_task_status_t execute_status) {
  index_shard_helper_group_t *group;
  index_shard_helper_task_t *task;
  anbool canonical_ready = FALSE;
  int rc = 0;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;

  if (execute_status != INDEX_SHARD_HELPER_TASK_OK &&
      execute_status != INDEX_SHARD_HELPER_TASK_STOPPED &&
      execute_status != INDEX_SHARD_HELPER_TASK_ERROR) {
    execute_status = INDEX_SHARD_HELPER_TASK_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  task = NULL;
  if (claim->task_index < group->task_count &&
      group->tasks) {
    task = &group->tasks[claim->task_index];
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }

  /*
   * Every returning callback retires one running lifetime lease even when
   * task bookkeeping is already inconsistent. The owner may then cancel
   * READY work and safely wait for the remaining callbacks.
   */
  if (group->running_count) {
    group->running_count--;
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }
  if (group->completed_count < group->task_count) {
    group->completed_count++;
  } else {
    group->internal_error = TRUE;
    rc = -1;
  }

  if (!task || task->scheduler_state !=
      INDEX_SHARD_HELPER_TASK_RUNNING) {
    group->internal_error = TRUE;
    rc = -1;
  } else {
    task->execute_status = execute_status;
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_DONE;
    canonical_ready =
        group->retire &&
        claim->task_index == group->next_retire &&
        execute_status == INDEX_SHARD_HELPER_TASK_OK;
  }

  if (execute_status == INDEX_SHARD_HELPER_TASK_ERROR) {
    group->task_failed = TRUE;
    shared->helper_task_failures++;
  } else if (execute_status ==
             INDEX_SHARD_HELPER_TASK_STOPPED) {
    group->stop_seen = TRUE;
  }
  if (group->internal_error || group->task_failed) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
  } else if (group->stop_seen) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  if (group->completed_count == group->task_count ||
      !group->running_count || group->internal_error ||
      group->task_failed || group->stop_seen ||
      canonical_ready) {
    pthread_cond_broadcast(&shared->queue_cv);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

static int index_shard_helper_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim) {
  index_shard_helper_group_t *group;
  index_shard_helper_task_t *task;
  index_shard_helper_task_status_t execute_status =
      INDEX_SHARD_HELPER_TASK_ERROR;
  int internal_error = FALSE;
  int rc;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;
  if (!group->ops || !group->ops->execute ||
      !group->tasks || claim->task_index >= group->task_count) {
    internal_error = TRUE;
  } else {
    task = &group->tasks[claim->task_index];
    execute_status = group->ops->execute(
        task->input,
        task->input_bytes,
        task->output,
        task->output_bytes);
  }
  rc = index_shard_helper_complete_claim(
      shared, claim, execute_status);
  if (internal_error) {
    pthread_mutex_lock(&shared->queue_mutex);
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
    pthread_cond_broadcast(&shared->queue_cv);
    pthread_mutex_unlock(&shared->queue_mutex);
    return -1;
  }
  return rc;
}

/* queue_mutex must be held. */
static void index_shard_staged_set_failed_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task) {
  if (!shared || !group || !task) {
    return;
  }
  (void)index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_FAILED);
  group->task_failed = TRUE;
  shared->staged_task_failures++;
}

/* queue_mutex must be held. */
static int index_shard_staged_release_ticket_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  int invalid = FALSE;

  if (!shared || !group) {
    return -1;
  }
  if (group->io_submitted) {
    group->io_submitted--;
  } else {
    invalid = TRUE;
  }
  group->io_completed++;
  if (shared->staged_tickets_active) {
    shared->staged_tickets_active--;
  } else {
    invalid = TRUE;
  }
  if (shared->staged_source_leases) {
    shared->staged_source_leases--;
  } else {
    invalid = TRUE;
  }
  shared->staged_io_completed++;
  if (invalid) {
    group->internal_error = TRUE;
  }
  return invalid ? -1 : 0;
}

/* queue_mutex must be held. */
static int index_shard_staged_completion_id_active_locked(
    const index_shard_staged_group_t *current_group,
    const index_shard_staged_task_t *current_task,
    unsigned long long completion_id) {
  index_shard_pool_t *pool;
  int owner;

  if (!current_group || !current_group->pool || !completion_id) {
    return -1;
  }
  pool = current_group->pool;
  for (owner = 0; owner < pool->shared.worker_count; owner++) {
    const index_shard_staged_group_t *group =
        pool->contexts[owner].published_staged_group;
    size_t task_index;

    if (!group) {
      continue;
    }
    for (task_index = 0U; task_index < group->task_count; task_index++) {
      const index_shard_staged_task_t *task = &group->tasks[task_index];

      if (group == current_group && task == current_task) {
        continue;
      }
      if (task->completion_id == completion_id) {
        return 1;
      }
    }
  }
  return 0;
}

static int index_shard_staged_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim,
    int callback_status,
    double callback_seconds,
    unsigned long long completion_id) {
  index_shard_staged_group_t *group;
  index_shard_staged_task_t *task;
  index_shard_staged_task_state_t expected;
  double now;
  anbool task_failed_before;
  int rc = 0;

  if (!shared || !claim || !claim->group) {
    return -1;
  }
  group = claim->group;
  pthread_mutex_lock(&shared->queue_mutex);
  if (claim->task_index >= group->task_count ||
      !group->tasks || !group->running_count) {
    group->internal_error = TRUE;
    if (group->running_count) {
      group->running_count--;
    }
    pthread_cond_broadcast(&shared->queue_cv);
    pthread_mutex_unlock(&shared->queue_mutex);
    return -1;
  }
  task = &group->tasks[claim->task_index];
  task_failed_before = group->task_failed;
  now = monotonic_seconds();
  switch (claim->kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    expected = INDEX_SHARD_STAGED_TASK_PREPARING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    expected = INDEX_SHARD_STAGED_TASK_SUBMITTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    expected = INDEX_SHARD_STAGED_TASK_IO_POLLING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    expected = INDEX_SHARD_STAGED_TASK_IO_CANCELLING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
    expected = INDEX_SHARD_STAGED_TASK_EXECUTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    expected = INDEX_SHARD_STAGED_TASK_OWNER_EXECUTING;
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    expected = INDEX_SHARD_STAGED_TASK_UNUSED;
    group->internal_error = TRUE;
    rc = -1;
    break;
  }
  group->running_count--;
  if (claim->kind == INDEX_SHARD_STAGED_CLAIM_EXECUTE) {
    if (!group->compute_running) {
      group->internal_error = TRUE;
      rc = -1;
    } else {
      group->compute_running--;
    }
  }
  task->callback_status = callback_status;
  if (task->scheduler_state != expected) {
    group->internal_error = TRUE;
    rc = -1;
  }
  if (claim->kind != INDEX_SHARD_STAGED_CLAIM_SUBMIT && completion_id) {
    group->internal_error = TRUE;
    rc = -1;
  }

  if (claim->kind == INDEX_SHARD_STAGED_CLAIM_PREPARE) {
    switch ((index_shard_staged_prepare_status_t)callback_status) {
    case INDEX_SHARD_STAGED_PREPARE_MORE:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_PREPARE_READY);
      break;
    case INDEX_SHARD_STAGED_PREPARE_SUBMIT_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_SUBMIT_READY);
      task->scheduler_epoch = ~shared->staged_completion_epoch;
      break;
    case INDEX_SHARD_STAGED_PREPARE_COMPUTE_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
      task->scheduler_ready_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_OWNER_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_OWNER_READY);
      task->scheduler_ready_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_RESULTS_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_RESULTS_READY);
      task->scheduler_result_seconds = now;
      break;
    case INDEX_SHARD_STAGED_PREPARE_STOPPED:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_PREPARE_ERROR:
    default:
      index_shard_staged_set_failed_locked(shared, group, task);
      break;
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_SUBMIT) {
    if (((index_shard_staged_submit_status_t)callback_status ==
             INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED) !=
            (completion_id != 0ULL) ||
        (completion_id &&
         index_shard_staged_completion_id_active_locked(
             group, task, completion_id) != 0)) {
      group->internal_error = TRUE;
      index_shard_staged_set_failed_locked(shared, group, task);
      rc = -1;
    } else {
      switch ((index_shard_staged_submit_status_t)callback_status) {
      case INDEX_SHARD_STAGED_SUBMIT_RETRY:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_SUBMIT_READY);
        /*
         * Preserve the epoch observed when this callback was claimed. If a
         * completion arrived while the queue mutex was released, the mismatch
         * remains visible and this retry is immediately eligible.
         */
        task->scheduler_epoch = claim->observed_completion_epoch;
        shared->staged_submit_retries++;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_IO_SUBMITTED:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
        task->scheduler_epoch = ~shared->staged_completion_epoch;
        task->cancel_sent = FALSE;
        task->completion_id = completion_id;
        task->completion_pending =
            shared->staged_completion_epoch !=
                claim->observed_completion_epoch;
        task->scheduler_submit_seconds = now;
        group->io_submitted++;
        if (group->io_submitted > group->max_io_submitted) {
          group->max_io_submitted = group->io_submitted;
        }
        shared->staged_tickets_active++;
        shared->staged_max_io_submitted = MAX(
            shared->staged_max_io_submitted,
            shared->staged_tickets_active);
        shared->staged_source_leases++;
        shared->staged_io_submitted++;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_COMPUTE_READY:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
        task->scheduler_ready_seconds = now;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_OWNER_READY:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_OWNER_READY);
        task->scheduler_ready_seconds = now;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_STOPPED:
        (void)index_shard_staged_set_state_locked(
            shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
        group->stop_seen = TRUE;
        break;
      case INDEX_SHARD_STAGED_SUBMIT_ERROR:
      default:
        index_shard_staged_set_failed_locked(shared, group, task);
        break;
      }
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_IO_POLL) {
    switch ((index_shard_staged_io_status_t)callback_status) {
    case INDEX_SHARD_STAGED_IO_PENDING:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
      /* completion_pending preserves a notifier racing this unlocked poll. */
      break;
    case INDEX_SHARD_STAGED_IO_READY:
      (void)index_shard_staged_release_ticket_locked(shared, group);
      task->completion_id = 0ULL;
      task->completion_pending = FALSE;
      if (task->scheduler_submit_seconds > 0.0 &&
          now >= task->scheduler_submit_seconds) {
        group->submit_to_ready_seconds +=
            now - task->scheduler_submit_seconds;
      }
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_COMPUTE_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      } else {
        task->scheduler_ready_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_IO_FAILED:
      (void)index_shard_staged_release_ticket_locked(shared, group);
      task->completion_id = 0ULL;
      task->completion_pending = FALSE;
      if (task->scheduler_submit_seconds > 0.0 &&
          now >= task->scheduler_submit_seconds) {
        group->submit_to_ready_seconds +=
            now - task->scheduler_submit_seconds;
      }
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_OWNER_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      } else {
        task->scheduler_ready_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_IO_CANCELLED:
      (void)index_shard_staged_release_ticket_locked(shared, group);
      task->completion_id = 0ULL;
      task->completion_pending = FALSE;
      task->scheduler_submit_seconds = 0.0;
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_IO_ERROR:
    default:
      (void)index_shard_staged_release_ticket_locked(shared, group);
      task->completion_id = 0ULL;
      task->completion_pending = FALSE;
      task->scheduler_submit_seconds = 0.0;
      index_shard_staged_set_failed_locked(shared, group, task);
      break;
    }
  } else if (claim->kind == INDEX_SHARD_STAGED_CLAIM_IO_CANCEL) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_IO_SUBMITTED);
    task->scheduler_epoch = ~shared->staged_completion_epoch;
    if (callback_status < 0) {
      /* Poll once, then leave a still-pending cancellation retriable. */
      task->cancel_sent = FALSE;
      task->completion_pending = TRUE;
      group->internal_error = TRUE;
      rc = -1;
    } else {
      task->cancel_sent = TRUE;
    }
  } else {
    group->execute_seconds += callback_seconds;
    switch ((index_shard_staged_execute_status_t)callback_status) {
    case INDEX_SHARD_STAGED_EXECUTE_MORE:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_PREPARE_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
      }
      break;
    case INDEX_SHARD_STAGED_EXECUTE_OK:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, group->cancelling
              ? INDEX_SHARD_STAGED_TASK_STOPPED
              : INDEX_SHARD_STAGED_TASK_RESULTS_READY);
      if (group->cancelling) {
        group->stop_seen = TRUE;
        task->scheduler_result_seconds = 0.0;
      } else {
        task->scheduler_result_seconds = now;
      }
      break;
    case INDEX_SHARD_STAGED_EXECUTE_STOPPED:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      break;
    case INDEX_SHARD_STAGED_EXECUTE_ERROR:
    default:
      index_shard_staged_set_failed_locked(shared, group, task);
      break;
    }
  }
  if (!task_failed_before && group->task_failed) {
    logerr("[index-shard] staged callback failed ops=%s task=%zu "
           "claim=%i status=%i state=%i\n",
           group->ops && group->ops->name
               ? group->ops->name
               : "<unnamed>",
           claim->task_index,
           (int)claim->kind,
           callback_status,
           task->scheduler_state);
  }
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

static int index_shard_staged_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_staged_claim_t *claim) {
  index_shard_staged_group_t *group;
  index_shard_staged_task_t *task;
  int callback_status = INDEX_SHARD_STAGED_EXECUTE_ERROR;
  unsigned long long completion_id = 0ULL;
  double callback_start;
  double callback_seconds;

  if (!shared || !claim || !claim->group ||
      claim->task_index >= claim->group->task_count) {
    return -1;
  }
  group = claim->group;
  task = &group->tasks[claim->task_index];
  if (!group->ops) {
    return index_shard_staged_complete_claim(
        shared, claim, callback_status, 0.0, 0ULL);
  }
  callback_start = monotonic_seconds();
  switch (claim->kind) {
  case INDEX_SHARD_STAGED_CLAIM_PREPARE:
    if (group->ops->prepare) {
      callback_status = group->ops->prepare(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_SUBMIT:
    if (group->ops->submit) {
      callback_status = group->ops->submit(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes,
          &completion_id);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_POLL:
    if (group->ops->poll) {
      callback_status = group->ops->poll(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_IO_CANCEL:
    if (group->ops->cancel) {
      callback_status = group->ops->cancel(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_EXECUTE:
    if (group->ops->execute) {
      callback_status = group->ops->execute(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_OWNER:
    if (group->ops->owner) {
      index_shard_worker_context_t *ctx =
          index_shard_get_tls();

      if (!ctx || ctx->pool != group->pool ||
          ctx->worker_id != group->owner_worker ||
          ctx->published_staged_group != group ||
          ctx->staged_owner_callback_active ||
          ctx->staged_owner_callback_group) {
        callback_status =
            INDEX_SHARD_STAGED_EXECUTE_ERROR;
        break;
      }
      ctx->staged_owner_callback_active = TRUE;
      ctx->staged_owner_callback_group = group;
      callback_status = group->ops->owner(
          task->input,
          task->input_bytes,
          task->output,
          task->output_bytes);
      ctx->staged_owner_callback_group = NULL;
      ctx->staged_owner_callback_active = FALSE;
    }
    break;
  case INDEX_SHARD_STAGED_CLAIM_NONE:
  default:
    break;
  }
  callback_seconds = monotonic_seconds() - callback_start;
  return index_shard_staged_complete_claim(
      shared, claim, callback_status, callback_seconds,
      completion_id);
}

static int index_shard_inner_execute_claim(
    index_shard_thread_state_t *shared,
    const index_shard_inner_claim_t *claim) {
  if (!shared || !claim) {
    return -1;
  }
  if (claim->kind == INDEX_SHARD_INNER_CLAIM_HELPER) {
    return index_shard_helper_execute_claim(
        shared, &claim->helper);
  }
  if (claim->kind == INDEX_SHARD_INNER_CLAIM_STAGED) {
    return index_shard_staged_execute_claim(
        shared, &claim->staged);
  }
  return -1;
}

/* queue_mutex must be held. */
static void index_shard_staged_cancel_ready_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  size_t i;

  if (!shared || !group || !group->tasks) {
    return;
  }
  group->cancelling = TRUE;
  for (i = 0U; i < group->task_count; i++) {
    index_shard_staged_task_t *task = &group->tasks[i];

    switch ((index_shard_staged_task_state_t)
                task->scheduler_state) {
    case INDEX_SHARD_STAGED_TASK_PREPARE_READY:
    case INDEX_SHARD_STAGED_TASK_SUBMIT_READY:
    case INDEX_SHARD_STAGED_TASK_COMPUTE_READY:
    case INDEX_SHARD_STAGED_TASK_OWNER_READY:
    case INDEX_SHARD_STAGED_TASK_RESULTS_READY:
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      task->scheduler_result_seconds = 0.0;
      break;
    default:
      break;
    }
  }
}

/* queue_mutex must be held. */
static int index_shard_staged_cancel_for_pool_locked(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  int fatal;
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested || shared->solved_published;
  pthread_mutex_unlock(&shared->state_mutex);
  if (fatal) {
    group->internal_error = TRUE;
    index_shard_staged_cancel_ready_locked(shared, group);
    return TRUE;
  }
  if (stopped || group->task_failed || group->stop_seen ||
      group->internal_error) {
    index_shard_staged_cancel_ready_locked(shared, group);
    return TRUE;
  }
  return FALSE;
}

/* queue_mutex must be held. */
static anbool index_shard_staged_all_terminal_locked(
    const index_shard_staged_group_t *group) {
  size_t i;

  if (!group || !group->tasks) {
    return FALSE;
  }
  for (i = 0U; i < group->task_count; i++) {
    if (!index_shard_staged_task_terminal(
            &group->tasks[i])) {
      return FALSE;
    }
  }
  return TRUE;
}

/*
 * Retire at most one complete logical task. Only the outer owner calls this
 * function. Slice preparation and execution may finish out of order, but the
 * callback sees exactly one final task in canonical task-index order.
 */
static int index_shard_staged_retire_one(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group) {
  index_shard_staged_task_t *task;
  index_shard_staged_retire_status_t status;
  size_t task_index;
  double retire_start;
  double retire_seconds;

  if (!shared || !group || !group->retire) {
    return 1;
  }
  pthread_mutex_lock(&shared->queue_mutex);
  if (group->next_retire >= group->task_count) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  if (index_shard_staged_cancel_for_pool_locked(
          shared, group)) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task_index = group->next_retire;
  task = &group->tasks[task_index];
  if (task->scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  retire_start = monotonic_seconds();
  if (task->scheduler_result_seconds > 0.0 &&
      retire_start >= task->scheduler_result_seconds) {
    group->result_to_retire_seconds +=
        retire_start - task->scheduler_result_seconds;
  }
  task->scheduler_result_seconds = 0.0;
  (void)index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_RETIRING);
  pthread_mutex_unlock(&shared->queue_mutex);

  status = group->retire(task, task_index, group->owner_context);
  retire_seconds = monotonic_seconds() - retire_start;
  if (status != INDEX_SHARD_STAGED_RETIRE_OK &&
      status != INDEX_SHARD_STAGED_RETIRE_STOPPED &&
      status != INDEX_SHARD_STAGED_RETIRE_MORE &&
      status != INDEX_SHARD_STAGED_RETIRE_ERROR) {
    status = INDEX_SHARD_STAGED_RETIRE_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  group->retire_seconds += retire_seconds;
  if (task->scheduler_state !=
          INDEX_SHARD_STAGED_TASK_RETIRING ||
      group->next_retire != task_index) {
    group->internal_error = TRUE;
    status = INDEX_SHARD_STAGED_RETIRE_ERROR;
  } else if (status == INDEX_SHARD_STAGED_RETIRE_OK) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_RETIRED);
    group->next_retire++;
  } else if (status ==
             INDEX_SHARD_STAGED_RETIRE_MORE) {
    if (index_shard_staged_cancel_for_pool_locked(
            shared, group)) {
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
      group->stop_seen = TRUE;
      status = INDEX_SHARD_STAGED_RETIRE_STOPPED;
    } else {
      (void)index_shard_staged_set_state_locked(
          shared, group, task, INDEX_SHARD_STAGED_TASK_PREPARE_READY);
    }
  } else if (status ==
             INDEX_SHARD_STAGED_RETIRE_STOPPED) {
    (void)index_shard_staged_set_state_locked(
        shared, group, task, INDEX_SHARD_STAGED_TASK_STOPPED);
    group->stop_seen = TRUE;
    index_shard_staged_cancel_ready_locked(shared, group);
  } else {
    index_shard_staged_set_failed_locked(shared, group, task);
    index_shard_staged_cancel_ready_locked(shared, group);
  }
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
  return status == INDEX_SHARD_STAGED_RETIRE_OK ||
      status == INDEX_SHARD_STAGED_RETIRE_MORE
      ? 0
      : -1;
}

/* queue_mutex must be held. */
static int index_shard_helper_cancel_for_pool_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  int fatal;
  int stopped;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested ||
      shared->solved_published;
  pthread_mutex_unlock(&shared->state_mutex);

  if (fatal) {
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
    return TRUE;
  }
  if (stopped) {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
    return TRUE;
  }
  return FALSE;
}

/* queue_mutex must be held. */
static int index_shard_helper_owner_claim_locked(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group,
    size_t *task_index) {
  int fatal;
  int stopped;
  int rc;

  pthread_mutex_lock(&shared->state_mutex);
  fatal = shared->fatal_error;
  stopped = shared->stop_requested ||
      shared->solved_published;
  if (!fatal && !stopped) {
    rc = index_shard_helper_claim_locked(
        group, task_index);
    pthread_mutex_unlock(&shared->state_mutex);
    return rc;
  }
  pthread_mutex_unlock(&shared->state_mutex);

  if (fatal) {
    group->internal_error = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_ERROR);
  } else {
    index_shard_helper_cancel_ready_locked(
        shared,
        group, INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  return 1;
}

/*
 * Retire at most one completed canonical task. The queue lock publishes the
 * helper output before the owner callback observes it. The callback runs
 * without pool locks and may mutate only the owning solver.
 */
static int index_shard_helper_retire_one(
    index_shard_thread_state_t *shared,
    index_shard_helper_group_t *group) {
  index_shard_helper_task_t *task;
  index_shard_helper_retire_status_t status;
  size_t task_index;

  if (!shared || !group || !group->retire) {
    return 1;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (group->next_retire >= group->task_count) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task_index = group->next_retire;
  task = &group->tasks[task_index];
  if (task->scheduler_state != INDEX_SHARD_HELPER_TASK_DONE) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  if (task->execute_status != INDEX_SHARD_HELPER_TASK_OK) {
    pthread_mutex_unlock(&shared->queue_mutex);
    return 1;
  }
  task->scheduler_state = INDEX_SHARD_HELPER_TASK_RETIRING;
  pthread_mutex_unlock(&shared->queue_mutex);

  status = group->retire(
      task, task_index, group->owner_context);
  if (status != INDEX_SHARD_HELPER_RETIRE_OK &&
      status != INDEX_SHARD_HELPER_RETIRE_STOPPED &&
      status != INDEX_SHARD_HELPER_RETIRE_ERROR) {
    status = INDEX_SHARD_HELPER_RETIRE_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (task->scheduler_state !=
      INDEX_SHARD_HELPER_TASK_RETIRING ||
      group->next_retire != task_index) {
    group->internal_error = TRUE;
    status = INDEX_SHARD_HELPER_RETIRE_ERROR;
  } else {
    task->scheduler_state = INDEX_SHARD_HELPER_TASK_RETIRED;
    group->next_retire++;
  }
  if (status == INDEX_SHARD_HELPER_RETIRE_ERROR) {
    group->task_failed = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group,
        INDEX_SHARD_HELPER_TASK_ERROR);
  } else if (status ==
             INDEX_SHARD_HELPER_RETIRE_STOPPED) {
    group->stop_seen = TRUE;
    index_shard_helper_cancel_ready_locked(
        shared,
        group,
        INDEX_SHARD_HELPER_TASK_STOPPED);
  }
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
  return status == INDEX_SHARD_HELPER_RETIRE_OK ? 0 : -1;
}

/*
 * A staged owner may temporarily expose one synchronous helper group while its
 * owner callback is executing. The staged task remains the lifetime owner; the
 * child group receives only immutable inputs and disjoint output ranges.
 */
static anbool index_shard_helper_staged_child_allowed(
    const index_shard_worker_context_t *ctx) {
  return ctx &&
      ctx->staged_owner_callback_active &&
      ctx->staged_owner_callback_group &&
      ctx->published_staged_group ==
          ctx->staged_owner_callback_group &&
      !ctx->published_helper_group;
}

size_t index_shard_helper_available_workers(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t available = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      (!ctx->published_staged_group ||
       index_shard_helper_staged_child_allowed(ctx)) &&
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active) {
    available = index_shard_helper_idle_workers_locked(
        shared);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return available;
}

size_t index_shard_helper_prepare_reserve(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t available = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active) {
    available = index_shard_helper_idle_workers_locked(shared);
    if (available) {
      shared->helper_preparations_active++;
      ctx->helper_preparation_active = TRUE;
      ctx->helper_preparation_generation =
          ctx->generation_seen;
      ctx->helper_preparation_index_order =
          ctx->current_index_order;
      ctx->helper_preparation_workers = available;
    }
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return available;
}

/* queue_mutex must be held. */
static int index_shard_helper_prepare_clear_locked(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared) {
  int invalid = FALSE;

  if (!ctx || !shared || !ctx->helper_preparation_active) {
    return 0;
  }
  if (shared->helper_preparations_active != 1U) {
    logerr("[index-shard] invalid helper preparation count=%zu\n",
           shared->helper_preparations_active);
    invalid = TRUE;
  } else {
    shared->helper_preparations_active--;
  }
  ctx->helper_preparation_active = FALSE;
  ctx->helper_preparation_generation = 0U;
  ctx->helper_preparation_index_order = SIZE_MAX;
  ctx->helper_preparation_workers = 0U;
  pthread_cond_broadcast(&shared->queue_cv);
  return invalid ? -1 : 0;
}

void index_shard_helper_prepare_cancel(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;

  if (!ctx || !ctx->pool) {
    return;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (index_shard_helper_prepare_clear_locked(
          ctx, shared)) {
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
}

static index_shard_helper_run_status_t
index_shard_helper_run_internal(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  index_shard_helper_group_t group;
  index_shard_helper_claim_t claim;
  index_shard_helper_run_status_t result;
  size_t i;
  int have_claim = TRUE;
  int wait_broken = FALSE;
  int fatal_requested = FALSE;
  int prepublish_fatal = FALSE;
  int preparation_permit = FALSE;
  int staged_child_permit = FALSE;
  int helper_window_active = FALSE;

  if (stats) {
    memset(stats, 0, sizeof(*stats));
  }
  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 || task_count < 2U) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (!ops || !ops->execute || !tasks ||
      task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }

  memset(&group, 0, sizeof(group));
  group.ops = ops;
  group.tasks = tasks;
  group.task_count = task_count;
  group.retire = retire;
  group.owner_context = owner_context;
  group.generation = ctx->generation_seen;
  group.owner_epoch = ++ctx->helper_group_epoch;
  group.owner_worker = ctx->worker_id;
  group.owner_index_order = ctx->current_index_order;

  for (i = 0U; i < task_count; i++) {
    unsigned long long work;

    if ((!tasks[i].input && tasks[i].input_bytes) ||
        (!tasks[i].output && tasks[i].output_bytes)) {
      index_shard_helper_prepare_cancel();
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
    tasks[i].scheduler_state = INDEX_SHARD_HELPER_TASK_READY;
    tasks[i].execute_status = INDEX_SHARD_HELPER_TASK_ERROR;
    work = index_shard_helper_task_work(&tasks[i]);
    if (index_shard_helper_add_work(
            &group.ready_work, work)) {
      index_shard_helper_prepare_cancel();
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
  }
  group.ready_count = task_count;

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  staged_child_permit =
      index_shard_helper_staged_child_allowed(ctx);
  preparation_permit =
      ctx->helper_preparation_active &&
      ctx->helper_preparation_generation ==
          ctx->generation_seen &&
      ctx->helper_preparation_index_order ==
          ctx->current_index_order;
  if (ctx->generation_seen != ctx->pool->generation ||
      ctx->published_helper_group ||
      (ctx->published_staged_group &&
       !staged_child_permit) ||
      (!preparation_permit &&
       shared->helper_preparations_active)) {
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      return INDEX_SHARD_HELPER_FATAL;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  group.foreign_reserve =
      index_shard_helper_idle_workers_locked(shared);
  if (preparation_permit &&
      group.foreign_reserve >
          ctx->helper_preparation_workers) {
    group.foreign_reserve =
        ctx->helper_preparation_workers;
  }
  if (group.foreign_reserve > task_count - 1U) {
    group.foreign_reserve = task_count - 1U;
  }
  if (!group.foreign_reserve) {
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      return INDEX_SHARD_HELPER_FATAL;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (SIZE_MAX - shared->helper_foreign_reservations <
      group.foreign_reserve) {
    (void)index_shard_helper_prepare_clear_locked(
        ctx, shared);
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return INDEX_SHARD_HELPER_FATAL;
  }

  memset(&claim, 0, sizeof(claim));
  claim.group = &group;
  result = INDEX_SHARD_HELPER_OK;
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->fatal_error) {
    prepublish_fatal = TRUE;
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (shared->stop_requested ||
             shared->solved_published) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else if (index_shard_helper_claim_locked(
                 &group, &claim.task_index)) {
    prepublish_fatal = TRUE;
    result = INDEX_SHARD_HELPER_FATAL;
  }
  if (result != INDEX_SHARD_HELPER_OK) {
    pthread_mutex_unlock(&shared->state_mutex);
    if (index_shard_helper_prepare_clear_locked(
            ctx, shared)) {
      result = INDEX_SHARD_HELPER_FATAL;
      prepublish_fatal = TRUE;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
    if (prepublish_fatal) {
      index_shard_request_fatal_stop(shared);
    }
    return result;
  }
  group.owner_claims++;
  group.owner_work += index_shard_helper_task_work(
      &group.tasks[claim.task_index]);
  shared->helper_tasks_owner++;
  if (index_shard_helper_prepare_clear_locked(
          ctx, shared)) {
    pthread_mutex_unlock(&shared->state_mutex);
    pthread_mutex_unlock(&shared->queue_mutex);
    index_shard_request_fatal_stop(shared);
    return INDEX_SHARD_HELPER_FATAL;
  }
  group.foreign_reservations_outstanding =
      group.foreign_reserve;
  shared->helper_foreign_reservations +=
      group.foreign_reservations_outstanding;
  ctx->published_helper_group = &group;
  shared->helper_groups_active++;
  shared->helper_groups_published++;
  pthread_mutex_unlock(&shared->state_mutex);
  fitsbin_payload_io_begin_helper_window();
  helper_window_active = TRUE;
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
  fitsbin_payload_io_notify_wait_helpers();

  while (1) {
    if (have_claim) {
      (void)index_shard_helper_execute_claim(
          shared, &claim);

      have_claim = FALSE;
    }

    while (!index_shard_helper_retire_one(
               shared, &group)) {
      /* Retire every currently completed canonical packet. */
    }

    pthread_mutex_lock(&shared->queue_mutex);
    (void)index_shard_helper_cancel_for_pool_locked(
        shared, &group);
    if (group.internal_error || group.task_failed) {
      index_shard_helper_cancel_ready_locked(
          shared,
          &group, INDEX_SHARD_HELPER_TASK_ERROR);
    } else if (group.stop_seen) {
      index_shard_helper_cancel_ready_locked(
          shared,
          &group, INDEX_SHARD_HELPER_TASK_STOPPED);
    }

    if (group.internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }

    /*
     * As with staged groups, a foreign task can finish after the owner
     * retirement scan but before this lock is acquired. Consume the
     * canonical completed task before deciding that the group is quiescent
     * or waiting for another notification.
     */
    if (group.retire &&
        !group.task_failed && !group.stop_seen &&
        group.next_retire < group.task_count &&
        group.tasks[group.next_retire].scheduler_state ==
            INDEX_SHARD_HELPER_TASK_DONE &&
        group.tasks[group.next_retire].execute_status ==
            INDEX_SHARD_HELPER_TASK_OK) {
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (!group.ready_count && !group.running_count) {
      break;
    }

    if (!group.ready_count) {
      if (wait_broken) {
        struct timespec pause = { 0, 1000000L };

        pthread_mutex_unlock(&shared->queue_mutex);
        nanosleep(&pause, NULL);
        continue;
      }

      {
        int wait_status;
        double wait_start;

        shared->helper_owner_wait_calls++;
        wait_start = monotonic_seconds();
        wait_status = pthread_cond_wait(
            &shared->queue_cv, &shared->queue_mutex);
        shared->helper_owner_wait_seconds +=
            monotonic_seconds() - wait_start;
        if (wait_status) {
          group.internal_error = TRUE;
          index_shard_helper_cancel_ready_locked(
              shared,
              &group, INDEX_SHARD_HELPER_TASK_ERROR);
          wait_broken = TRUE;
        }
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (group.foreign_reserve > group.foreign_claims) {
      size_t outstanding =
          group.foreign_reserve - group.foreign_claims;

      if (group.ready_count <= outstanding &&
          !group.owner_reserve_yielded) {
        group.owner_reserve_yielded = TRUE;
        pthread_cond_broadcast(&shared->queue_cv);
        pthread_mutex_unlock(&shared->queue_mutex);
        sched_yield();
        continue;
      }
      if (group.owner_reserve_yielded) {
        if (index_shard_helper_release_foreign_reservations_locked(
                shared,
                &group,
                group.foreign_reservations_outstanding)) {
          index_shard_helper_cancel_ready_locked(
              shared,
              &group,
              INDEX_SHARD_HELPER_TASK_ERROR);
        }
        group.foreign_reserve = group.foreign_claims;
      }
    }

    claim.group = &group;
    {
      int claim_status =
          index_shard_helper_owner_claim_locked(
              shared, &group, &claim.task_index);

      if (claim_status < 0) {
        group.internal_error = TRUE;
        index_shard_helper_cancel_ready_locked(
            shared,
            &group, INDEX_SHARD_HELPER_TASK_ERROR);
      } else if (!claim_status) {
        group.owner_claims++;
              group.owner_work += index_shard_helper_task_work(
            &group.tasks[claim.task_index]);
        shared->helper_tasks_owner++;
        have_claim = TRUE;
        pthread_mutex_unlock(&shared->queue_mutex);
        continue;
      }
    }
    if (group.internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }
    pthread_mutex_unlock(&shared->queue_mutex);
  }

  if (index_shard_helper_release_foreign_reservations_locked(
          shared,
          &group,
          group.foreign_reservations_outstanding)) {
    group.internal_error = TRUE;
  }
  if (group.ready_count || group.running_count ||
      group.completed_count != group.task_count) {
    group.internal_error = TRUE;
  }
  if (group.retire &&
      !group.task_failed && !group.stop_seen &&
      group.next_retire != group.task_count) {
    group.internal_error = TRUE;
  }
  {
    size_t lanes_cleared = 0U;

    for (i = 0U; i < (size_t)shared->worker_count; i++) {
      if (ctx->pool->contexts[i].published_helper_group ==
          &group) {
        ctx->pool->contexts[i].published_helper_group = NULL;
        lanes_cleared++;
      }
    }
    if (lanes_cleared != 1U) {
      group.internal_error = TRUE;
    } else {
      shared->helper_groups_completed++;
    }
  }
  if (ctx->published_helper_group == &group) {
    ctx->published_helper_group = NULL;
    group.internal_error = TRUE;
  }
  if (!shared->helper_groups_active) {
    group.internal_error = TRUE;
  } else {
    shared->helper_groups_active--;
  }
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);

  if (helper_window_active) {
    fitsbin_payload_io_end_helper_window();
    helper_window_active = FALSE;
  }
  if (stats) {
    stats->owner_tasks = group.owner_claims;
    stats->foreign_tasks = group.foreign_claims;
    stats->max_concurrent_tasks = group.max_running;
    stats->owner_work_units = group.owner_work;
    stats->foreign_work_units = group.foreign_work;
  }
  if (group.internal_error) {
    if (!fatal_requested) {
      index_shard_request_fatal_stop(shared);
    }
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (group.task_failed) {
    result = INDEX_SHARD_HELPER_TASK_FAILED;
  } else if (group.stop_seen) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  return result;
}

index_shard_helper_run_status_t
index_shard_helper_run(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_run_stats_t *stats) {
  return index_shard_helper_run_internal(
      ops, tasks, task_count, NULL, NULL, stats);
}

index_shard_helper_run_status_t
index_shard_helper_run_ordered(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
    index_shard_helper_retire_fn retire,
    void *owner_context,
    index_shard_helper_run_stats_t *stats) {
  if (!retire) {
    index_shard_helper_prepare_cancel();
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  return index_shard_helper_run_internal(
      ops,
      tasks,
      task_count,
      retire,
      owner_context,
      stats);
}

size_t index_shard_staged_capacity(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t capacity = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active) {
    /* This is a storage bound, not an instantaneous idle-worker snapshot. */
    capacity = INDEX_SHARD_HELPER_MAX_TASKS;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return capacity;
}

size_t index_shard_staged_compute_width(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  size_t width = 0U;

  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered ||
      index_shard_worker_stop_requested()) {
    return 0U;
  }
  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen == ctx->pool->generation &&
      !ctx->published_helper_group &&
      !ctx->published_staged_group &&
      !ctx->helper_preparation_active) {
    width = (size_t)ctx->pool->worker_count;
  }
  pthread_mutex_unlock(&shared->queue_mutex);
  return width;
}

index_shard_helper_run_status_t
index_shard_staged_run_ordered(
    const index_shard_staged_ops_t *ops,
    index_shard_staged_task_t *tasks,
    size_t task_count,
    index_shard_staged_retire_fn retire,
    void *owner_context,
    index_shard_staged_run_stats_t *stats) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();
  index_shard_thread_state_t *shared;
  index_shard_staged_group_t *group;
  index_shard_helper_run_status_t result;
  size_t i;
  int fatal_requested = FALSE;
  int wait_broken = FALSE;

  if (stats) {
    memset(stats, 0, sizeof(*stats));
  }
  if (!ctx || !ctx->pool || !ctx->current_outer_active ||
      ctx->pool->worker_count < 2 ||
      !ctx->pool->payload_completion_registered) {
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  if (!ops || !ops->prepare || !ops->submit || !ops->poll ||
      !ops->cancel || !ops->execute || !ops->owner || !retire ||
      !tasks || !task_count ||
      task_count > INDEX_SHARD_HELPER_MAX_TASKS) {
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  for (i = 0U; i < task_count; i++) {
    if ((!tasks[i].input && tasks[i].input_bytes) ||
        (!tasks[i].output && tasks[i].output_bytes)) {
      return INDEX_SHARD_HELPER_TASK_FAILED;
    }
  }

  group = calloc(1, sizeof(*group));
  if (!group) {
    return INDEX_SHARD_HELPER_TASK_FAILED;
  }
  group->pool = ctx->pool;
  group->ops = ops;
  group->tasks = tasks;
  group->task_count = task_count;
  group->retire = retire;
  group->owner_context = owner_context;
  group->generation = ctx->generation_seen;
  group->owner_epoch = ++ctx->staged_group_epoch;
  if (!group->owner_epoch) {
    group->owner_epoch = ++ctx->staged_group_epoch;
  }
  group->owner_worker = ctx->worker_id;
  group->owner_index_order = ctx->current_index_order;

  shared = &ctx->pool->shared;
  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->generation_seen != ctx->pool->generation ||
      !ctx->current_outer_active ||
      ctx->published_helper_group ||
      ctx->published_staged_group ||
      ctx->helper_preparation_active ||
      !ctx->pool->payload_completion_registered) {
    pthread_mutex_unlock(&shared->queue_mutex);
    free(group);
    return INDEX_SHARD_HELPER_UNAVAILABLE;
  }
  pthread_mutex_lock(&shared->state_mutex);
  if (shared->fatal_error) {
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (shared->stop_requested || shared->solved_published) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  pthread_mutex_unlock(&shared->state_mutex);
  if (result != INDEX_SHARD_HELPER_OK) {
    pthread_mutex_unlock(&shared->queue_mutex);
    free(group);
    return result;
  }

  for (i = 0U; i < task_count; i++) {
    tasks[i].scheduler_state =
        INDEX_SHARD_STAGED_TASK_PREPARE_READY;
    tasks[i].cancel_sent = FALSE;
    tasks[i].completion_pending = FALSE;
    tasks[i].callback_status = INDEX_SHARD_STAGED_EXECUTE_ERROR;
    tasks[i].scheduler_epoch = ~shared->staged_completion_epoch;
    tasks[i].completion_id = 0ULL;
    tasks[i].scheduler_submit_seconds = 0.0;
    tasks[i].scheduler_ready_seconds = 0.0;
    tasks[i].scheduler_result_seconds = 0.0;
  }
  ctx->published_staged_group = group;
  shared->staged_groups_active++;
  shared->staged_groups_published++;
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);
  fitsbin_payload_io_notify_wait_helpers();

  while (1) {
    index_shard_inner_claim_t claim;
    int selection;

    while (!index_shard_staged_retire_one(shared, group)) {
      /* Retire every complete logical packet in canonical order. */
    }

    pthread_mutex_lock(&shared->queue_mutex);
    (void)index_shard_staged_cancel_for_pool_locked(
        shared, group);
    if (group->internal_error || group->task_failed ||
        group->stop_seen) {
      index_shard_staged_cancel_ready_locked(shared, group);
    }

    if (group->internal_error && !fatal_requested) {
      fatal_requested = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_request_fatal_stop(shared);
      continue;
    }

    /*
     * Retirement is owner-only and runs without the queue lock. A foreign
     * completion can publish the next canonical result after the retirement
     * scan but before this lock is acquired. Recheck that predicate here so
     * the owner never sleeps after the corresponding broadcast has passed.
     */
    if (!group->cancelling &&
        group->next_retire < group->task_count &&
        group->tasks[group->next_retire].scheduler_state ==
            INDEX_SHARD_STAGED_TASK_RESULTS_READY) {
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }

    if (!group->cancelling &&
        group->next_retire == group->task_count &&
        !group->running_count && !group->io_submitted) {
      break;
    }
    if (group->cancelling && !group->running_count &&
        !group->io_submitted &&
        index_shard_staged_all_terminal_locked(group)) {
      break;
    }

    selection = index_shard_inner_select_locked(
        ctx, shared, TRUE, &claim);
    if (selection < 0) {
      group->internal_error = TRUE;
      pthread_mutex_unlock(&shared->queue_mutex);
      continue;
    }
    if (!selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      if (index_shard_inner_execute_claim(shared, &claim)) {
        pthread_mutex_lock(&shared->queue_mutex);
        group->internal_error = TRUE;
        pthread_cond_broadcast(&shared->queue_cv);
        pthread_mutex_unlock(&shared->queue_mutex);
      }
      continue;
    }

    if (wait_broken) {
      struct timespec pause = { 0, 1000000L };

      pthread_mutex_unlock(&shared->queue_mutex);
      nanosleep(&pause, NULL);
      continue;
    }
    {
      int wait_status;
      double wait_start;

      shared->staged_owner_wait_calls++;
      wait_start = monotonic_seconds();
      wait_status = pthread_cond_wait(
          &shared->queue_cv, &shared->queue_mutex);
      shared->staged_owner_wait_seconds +=
          monotonic_seconds() - wait_start;
      if (wait_status) {
        group->internal_error = TRUE;
        wait_broken = TRUE;
      }
    }
    pthread_mutex_unlock(&shared->queue_mutex);
  }

  if (group->running_count || group->compute_running ||
      group->io_submitted || group->compute_ready ||
      group->reorder_ready) {
    group->internal_error = TRUE;
  }
  for (i = 0U; i < task_count; i++) {
    if (tasks[i].completion_id || tasks[i].completion_pending) {
      group->internal_error = TRUE;
    }
  }
  if (!group->cancelling &&
      group->next_retire != group->task_count) {
    group->internal_error = TRUE;
  }
  {
    size_t groups_cleared = 0U;

    for (i = 0U; i < (size_t)shared->worker_count; i++) {
      if (ctx->pool->contexts[i].published_staged_group == group) {
        ctx->pool->contexts[i].published_staged_group = NULL;
        groups_cleared++;
      }
    }
    if (groups_cleared != 1U) {
      group->internal_error = TRUE;
    } else {
      shared->staged_groups_completed++;
    }
  }
  if (!shared->staged_groups_active) {
    group->internal_error = TRUE;
  } else {
    shared->staged_groups_active--;
  }
  shared->staged_max_io_submitted = MAX(
      shared->staged_max_io_submitted,
      group->max_io_submitted);
  shared->staged_max_compute_ready = MAX(
      shared->staged_max_compute_ready,
      group->max_compute_ready);
  shared->staged_max_reorder_ready = MAX(
      shared->staged_max_reorder_ready,
      group->max_reorder_ready);
  shared->staged_max_compute_running = MAX(
      shared->staged_max_compute_running,
      group->max_compute_running);
  shared->staged_prepare_claims += group->prepare_claims;
  shared->staged_submit_claims += group->submit_claims;
  shared->staged_poll_claims += group->poll_claims;
  shared->staged_execute_claims += group->execute_claims;
  shared->staged_owner_execute_claims +=
      group->owner_execute_claims;
  shared->staged_submit_to_ready_seconds +=
      group->submit_to_ready_seconds;
  shared->staged_ready_dwell_seconds +=
      group->ready_dwell_seconds;
  shared->staged_execute_seconds += group->execute_seconds;
  shared->staged_result_to_retire_seconds +=
      group->result_to_retire_seconds;
  shared->staged_retire_seconds += group->retire_seconds;
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);

  if (stats) {
    stats->owner_claims = group->owner_claims;
    stats->foreign_claims = group->foreign_claims;
    stats->owner_compute_executes =
        group->owner_compute_executes;
    stats->foreign_compute_executes =
        group->foreign_compute_executes;
    stats->max_concurrent_claims = group->max_running;
    stats->max_compute_running = group->max_compute_running;
    stats->io_submitted = group->io_submitted + group->io_completed;
    stats->io_completed = group->io_completed;
    stats->max_io_submitted = group->max_io_submitted;
    stats->max_compute_ready = group->max_compute_ready;
    stats->max_reorder_ready = group->max_reorder_ready;
    stats->prepare_claims = group->prepare_claims;
    stats->submit_claims = group->submit_claims;
    stats->poll_claims = group->poll_claims;
    stats->execute_claims = group->execute_claims;
    stats->owner_claims_executed = group->owner_execute_claims;
    stats->submit_to_ready_seconds = group->submit_to_ready_seconds;
    stats->ready_dwell_seconds = group->ready_dwell_seconds;
    stats->execute_seconds = group->execute_seconds;
    stats->result_to_retire_seconds =
        group->result_to_retire_seconds;
    stats->retire_seconds = group->retire_seconds;
    stats->owner_work_units = group->owner_work;
    stats->foreign_work_units = group->foreign_work;
  }
  if (group->internal_error) {
    if (!fatal_requested) {
      index_shard_request_fatal_stop(shared);
    }
    result = INDEX_SHARD_HELPER_FATAL;
  } else if (group->task_failed) {
    result = INDEX_SHARD_HELPER_TASK_FAILED;
  } else if (group->stop_seen || group->cancelling) {
    result = INDEX_SHARD_HELPER_STOPPED;
  } else {
    result = INDEX_SHARD_HELPER_OK;
  }
  free(group);
  return result;
}

static void index_shard_worker_cleanup_pass(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared);

static void index_shard_advance_canonical_cursor_locked(
    index_shard_thread_state_t *shared) {
  while (shared->canonical_scan_cursor < shared->nindexes &&
         shared->outer_states[shared->canonical_scan_cursor] !=
             INDEX_SHARD_OUTER_UNCLAIMED) {
    shared->canonical_scan_cursor++;
  }
}

static int index_shard_claim_outer_locked(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t candidate,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice) {
  if (candidate >= shared->nindexes ||
      shared->outer_states[candidate] !=
          INDEX_SHARD_OUTER_UNCLAIMED ||
      !shared->outer_unclaimed ||
      shared->outer_running >= shared->producer_width) {
    return -1;
  }
  shared->outer_states[candidate] =
      INDEX_SHARD_OUTER_RUNNING;
  shared->outer_unclaimed--;
  shared->outer_running++;
  shared->outer_claims++;
  *index_order = candidate;
  *mmap_advice = shared->mmap_advice;
  if (candidate == shared->canonical_scan_cursor) {
    index_shard_advance_canonical_cursor_locked(shared);
  }

  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] claim index_order=%zu lane=producer "
           "worker=%i owners=%zu producer_width=%zu "
           "outer_unclaimed=%zu payload=%s "
           "wall_since_pass=%.6f\n",
           candidate,
           worker->worker_id,
           shared->outer_running,
           shared->producer_width,
           shared->outer_unclaimed,
           fitsbin_mmap_advice_name(*mmap_advice),
           monotonic_seconds() - shared->pass_wall_start);
  }
  return 0;
}

/*
 * Select work from the current band without transferring index ownership.
 *
 * New canonical outer work has priority whenever a producer slot is free.
 * Helper packages become eligible only when no outer task is immediately
 * claimable. The reducer remains the only authority for master result
 * mutation.
 */
static index_shard_work_selection_t
index_shard_select_work(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice,
    index_shard_inner_claim_t *inner_claim) {
  if (!worker || !shared || !index_order || !mmap_advice || !inner_claim ||
      !shared->outer_states || !shared->producer_width) {
    return INDEX_SHARD_WORK_ERROR;
  }
  if (worker->worker_id < 0 ||
      worker->worker_id >= shared->worker_count) {
    return INDEX_SHARD_WORK_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  while (1) {
    index_shard_pass_state_snapshot_t state;
    size_t candidate;
    int inner_selection;

    index_shard_pass_state_snapshot(shared, &state);
    if (state.stop_requested || state.fatal_error ||
        state.solved_published) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }

    index_shard_advance_canonical_cursor_locked(shared);
    candidate = shared->canonical_scan_cursor;
    if (candidate < shared->nindexes &&
        shared->outer_running < shared->producer_width) {
      if (index_shard_claim_outer_locked(
              worker,
              shared,
              candidate,
              index_order,
              mmap_advice)) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_OUTER;
    }

    inner_selection = index_shard_inner_select_locked(
        worker, shared, FALSE, inner_claim);
    if (inner_selection < 0) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    if (!inner_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_HELPER;
    }

    if (!shared->outer_running) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }
    if (worker->local_context_ready) {
      /*
       * No new outer task can become claimable for this worker in the
       * current band. Release its index-local context before it waits for
       * helper packages from the remaining owners.
       */
      pthread_mutex_unlock(&shared->queue_mutex);
      index_shard_worker_cleanup_pass(worker, shared);
      pthread_mutex_lock(&shared->queue_mutex);
      continue;
    }
    shared->queue_waiters++;
    inner_selection = pthread_cond_wait(
        &shared->queue_cv, &shared->queue_mutex);
    if (!shared->queue_waiters) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    shared->queue_waiters--;
    if (!shared->queue_waiters) {
      pthread_cond_broadcast(&shared->queue_cv);
    }
    if (inner_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
  }
}

/*
 * Freeze one index-independent scientific result and arbitrate the first
 * terminal event before index cleanup or producer-slot release.
 *
 * result_mutex -> state_mutex -> limit_mutex is the only lock order used
 * here. All three are released before queue_mutex can be acquired later by
 * index_shard_finish_outer_claim(). Result lifecycle metrics may still be
 * completed by the owner, but solutions, solved state and candidate identity
 * are immutable after candidate_ready is published.
 */
static int index_shard_arbitrate_candidate(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  index_shard_result_t *result;
  index_shard_terminal_cause_t cause =
      INDEX_SHARD_TERMINAL_NONE;
  double elapsed = 0.0;
  int rc = 0;
  int report = FALSE;
  int selected_now = FALSE;
  int stop_now = FALSE;

  if (!shared || !shared->results ||
      index_order >= shared->nindexes) {
    return -1;
  }

  result = &shared->results[index_order];

  pthread_mutex_lock(&shared->result_mutex);
  if (result->candidate_ready ||
      (shared->completed && shared->completed[index_order])) {
    logerr("[index-shard] duplicate or late candidate publication "
           "index_order=%zu\n",
           index_order);
    rc = -1;
  } else {
    if (result->solved &&
        (result->failed || result->rc ||
         result->failure_class != INDEX_SHARD_FAILURE_NONE)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    result->candidate_sequence =
        ++shared->next_candidate_sequence;
    result->candidate_ready = TRUE;

    pthread_mutex_lock(&shared->state_mutex);
    if (result->failure_class ==
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY) {
      (void)index_shard_publish_terminal_locked(
          shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
    } else if (shared->terminal_cause ==
               INDEX_SHARD_TERMINAL_NONE) {
      pthread_mutex_lock(&shared->limit_mutex);
      cause = index_shard_sample_terminal_locked(
          shared, result, &elapsed, &report);
      if (cause != INDEX_SHARD_TERMINAL_NONE) {
        (void)index_shard_publish_terminal_locked(
            shared, cause);
      } else if (result->solved &&
                 !result->failed &&
                 result->rc == 0 &&
                 result->failure_class ==
                     INDEX_SHARD_FAILURE_NONE &&
                 !shared->winner_selected &&
                 !shared->solved_published) {
        shared->winner_selected = TRUE;
        shared->selected_index_order = index_order;
        shared->selected_candidate_sequence =
            result->candidate_sequence;
        selected_now = index_shard_publish_terminal_locked(
            shared, INDEX_SHARD_TERMINAL_WINNER);
      }
      pthread_mutex_unlock(&shared->limit_mutex);
    }
    stop_now = shared->stop_requested;
    pthread_mutex_unlock(&shared->state_mutex);
  }

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);

  if (stop_now) {
    index_shard_publish_worker_stop(shared);
  }

  if (report && cause == INDEX_SHARD_TERMINAL_WALL_LIMIT) {
    logmsg("Total wall-clock time limit reached!\n");
    logverb("[index-shard] wall-limit reached total_timelimit=%g "
            "elapsed=%.6f\n",
            shared->bp->total_timelimit,
            elapsed);
  } else if (report && cause == INDEX_SHARD_TERMINAL_CPU_LIMIT) {
    logmsg("Total CPU time limit reached!\n");
    logverb("[index-shard] cpu-budget reached total_cpulimit=%g "
            "elapsed=%.6f\n",
            shared->bp->total_cpulimit,
            elapsed);
  }

  if (selected_now) {
    logverb("[index-shard] winner-selected index_order=%zu worker=%i "
            "candidate_sequence=%zu field=%i best_logodds=%.17g "
            "pass_wall=%.6f\n",
            index_order,
            result->worker_id,
            result->candidate_sequence,
            result->best_fieldnum,
            result->best_logodds,
            shared->first_stop_wall_since_pass);
  }

  return rc;
}

/*
 * Publish full task completion after candidate arbitration and index cleanup.
 *
 * candidate_ready protects the immutable scientific payload. completed[]
 * protects the remaining owner lifecycle and metrics. A cleanup failure can
 * still promote an elected-but-uncommitted winner to global-integrity failure,
 * but this function never performs winner election or releases outer producer
 * capacity.
 */
static int index_shard_mark_result_completed(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  index_shard_result_t *result;
  int fatal_now = FALSE;
  int late_loser_failure = FALSE;
  int rc = 0;
  int stop_now = FALSE;
  int task_local_now = FALSE;

  if (!shared || !shared->results || !shared->completed ||
      index_order >= shared->nindexes) {
    return -1;
  }

  result = &shared->results[index_order];

  /*
   * NOTE INDEX-SHARD: claimed-task-invariant
   *
   * Every claimed task must mark completion exactly once. Otherwise the
   * reducer can wait forever.
   */
  pthread_mutex_lock(&shared->result_mutex);
  if (shared->completed[index_order]) {
    logerr("[index-shard] duplicate result completion "
           "index_order=%zu\n",
           index_order);
    rc = -1;
  } else {
    /*
     * Freeze a complete cause-based terminal classification before publishing
     * the immutable result slot. Unclassified or internally inconsistent
     * failure state is itself a pass-integrity failure.
     */
    if ((result->failed || result->rc) &&
        result->failure_class == INDEX_SHARD_FAILURE_NONE) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    } else if (result->failure_class != INDEX_SHARD_FAILURE_NONE &&
               (!result->failed || !result->rc)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    result->completion_sequence =
        ++shared->next_completion_sequence;
    shared->completed[index_order] = TRUE;

    pthread_mutex_lock(&shared->state_mutex);

    if (result->solved && !result->candidate_ready) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }
    if (shared->winner_selected &&
        shared->selected_index_order == index_order &&
        (!result->candidate_ready ||
         !result->candidate_sequence ||
         result->candidate_sequence !=
             shared->selected_candidate_sequence ||
         result->failed || result->rc ||
         result->failure_class != INDEX_SHARD_FAILURE_NONE)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          result->rc);
    }

    if (result->failure_class ==
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY) {
      shared->global_integrity_failures++;
      (void)index_shard_publish_terminal_locked(
          shared, INDEX_SHARD_TERMINAL_GLOBAL_INTEGRITY);
      fatal_now = TRUE;
      stop_now = shared->stop_requested;
    } else if (result->failure_class ==
               INDEX_SHARD_FAILURE_TASK_LOCAL) {
      shared->task_local_failures++;
      task_local_now = TRUE;
      if (shared->winner_selected &&
          shared->selected_index_order != index_order) {
        shared->late_loser_failures++;
        late_loser_failure = TRUE;
      }
    }

    pthread_mutex_unlock(&shared->state_mutex);
  }

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);

  if (stop_now) {
    index_shard_publish_worker_stop(shared);
  }

  if (fatal_now) {
    logerr("[index-shard] global integrity failure terminated pass "
           "index_order=%zu completion_sequence=%zu class=%s\n",
           index_order,
           result->completion_sequence,
           index_shard_failure_class_name(result->failure_class));
  } else if (late_loser_failure) {
    logverb("[index-shard] preserving selected winner after task-local "
            "loser index_order=%zu completion_sequence=%zu\n",
            index_order,
            result->completion_sequence);
  } else if (task_local_now) {
    logverb("[index-shard] isolated task-local failure "
            "index_order=%zu completion_sequence=%zu\n",
            index_order,
            result->completion_sequence);
  }

  return rc;
}

static void index_shard_finish_outer_claim(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  int underflow = FALSE;
  int completion_failed;

  /*
   * Candidate arbitration has already frozen any scientific result and
   * published its terminal event. This final limit sample covers tasks that
   * produced no candidate or crossed a limit during cleanup. Full completion
   * remains visible before the producer slot becomes claimable.
   */
  (void)index_shard_check_global_limits(shared);
  completion_failed = index_shard_mark_result_completed(
      shared, index_order);
  if (completion_failed) {
    index_shard_request_fatal_stop(shared);
  }

  pthread_mutex_lock(&shared->queue_mutex);
  if (!shared->outer_states ||
      index_order >= shared->nindexes ||
      shared->outer_states[index_order] !=
          INDEX_SHARD_OUTER_RUNNING) {
    logerr("[index-shard] invalid outer lifecycle "
           "index_order=%zu\n",
           index_order);
    underflow = TRUE;
  } else {
    shared->outer_states[index_order] =
        INDEX_SHARD_OUTER_FINISHED;
  }
  if (!shared->outer_running) {
    logerr("[index-shard] outer-running underflow\n");
    underflow = TRUE;
  } else {
    shared->outer_running--;
  }
  if (!shared->outer_running) {
    pthread_cond_broadcast(&shared->queue_cv);
  }
  pthread_mutex_unlock(&shared->queue_mutex);

  if (underflow) {
    index_shard_request_fatal_stop(shared);
  }
}

/*
 * SECTION INDEX-SHARD: worker-context
 *
 * Worker-local context reuse.
 *
 * local_bp/local solver are prepared once per submitted pass and reused across
 * many one-index tasks.  This removes repeated xylist open/close and solver
 * index-list allocation from the per-index hot path.
 */
static int index_shard_worker_prepare_pass(index_shard_worker_context_t *ctx,
                                           index_shard_thread_state_t *shared) {
  // old local context belongs to a previous generation -> cleanup first
  if (ctx->local_context_ready && ctx->local_context_generation == ctx->generation_seen)
    return 0;

  if (ctx->local_context_ready) {
    if (shared->hooks && shared->hooks->cleanup_local_context)
      shared->hooks->cleanup_local_context(&ctx->local_bp);

    ctx->local_context_ready = FALSE;
  }

  if (!shared->hooks || !shared->hooks->prepare_local_context)
    return -1;

  if (!shared->worker_view ||
      shared->hooks->prepare_local_context(
          &ctx->local_bp, shared->worker_view))
    return -1;
  // hook copies stable master config + opens worker-local xylist
  ctx->local_context_ready = TRUE;
  ctx->local_context_generation = ctx->generation_seen;

  return 0;
}
// ANCHOR INDEX-SHARD: worker-cleanup-pass
/*
 * Release worker-local pass context after this generation is finished.
 *
 * Does not touch master bp and does not free result slots.
 */
static void index_shard_worker_cleanup_pass(index_shard_worker_context_t *ctx,
                                            index_shard_thread_state_t *shared) {
  if (!ctx->local_context_ready)
    return;

  if (shared->hooks && shared->hooks->cleanup_local_context)
    shared->hooks->cleanup_local_context(&ctx->local_bp);

  memset(&ctx->local_bp, 0, sizeof(onefield_t));
  ctx->local_context_ready = FALSE;
  ctx->local_context_generation = 0;
}

static index_shard_hook_result_t index_shard_done_with_index(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared,
    size_t index_order,
    index_t *index,
    index_shard_inverse_lease_t *inverse_lease) {
  index_shard_hook_result_t hook_result = {
      INDEX_SHARD_HOOK_COMPLETED_UNSOLVED, 0};

  if (!index) {
    return hook_result;
  }

  index_shard_inverse_cache_release(
      ctx, index, inverse_lease);
  if (!shared || !shared->hooks ||
      !shared->hooks->done_with_index) {
    hook_result.outcome =
        INDEX_SHARD_HOOK_GLOBAL_INTEGRITY_FAILURE;
    hook_result.error_code = -1;
    return hook_result;
  }

  return shared->hooks->done_with_index(
      shared->bp, index_order, index);
}

// ANCHOR INDEX-SHARD: run-one-index
/*
 * Execute one index shard in one worker.
 *
 * This function owns only local computation:
 *   - reset local context for this result slot
 *   - load one index
 *   - run solver against local_bp
 *   - analyze local solutions
 *   - release index
 *
 * It does not merge into master bp.
 */
static int index_shard_run_one_with_worker_context(index_shard_worker_context_t *ctx,
                                                   index_shard_thread_state_t *shared,
                                                   size_t index_order,
                                                   index_shard_result_t *result,
                                                   fitsbin_mmap_advice_t mmap_advice) {
  // one result slot belongs to this task
  index_t *index = NULL;
  double task_wall_start;
  double phase_wall_start;
  double wall_start;
  float cpu_start;
  index_shard_inverse_lease_t inverse_lease;
  index_shard_hook_result_t hook_result;
  int hook_status;

  index_shard_result_init(result, index_order);
  memset(&inverse_lease, 0, sizeof(inverse_lease));
  result->mmap_advice = mmap_advice;

  if (!result->solutions) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_TASK_LOCAL,
        -1);
    return -1;
  }

  if (!ctx->local_context_ready) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  if (!shared->hooks || !shared->hooks->reset_local_context_for_task ||
      !shared->hooks->solve_one_index) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
    return -1;
  }

  /*
   * Full outer-task timing starts before local reset and index acquisition.
   * Queue wait and reducer wait are deliberately excluded.
   */
  task_wall_start = monotonic_seconds();

  result->task_started = TRUE;
  result->worker_id = ctx->worker_id;
  result->task_start_since_pass =
      task_wall_start - shared->pass_wall_start;
#if defined(RUSAGE_THREAD)
  result->task_resource_valid =
      (getrusage(
           RUSAGE_THREAD,
           &result->task_resource_start) == 0);
#else
  result->task_resource_valid = FALSE;
#endif

   // local_bp reused across tasks, but solutions change per task
  phase_wall_start = monotonic_seconds();

  shared->hooks->reset_local_context_for_task(&ctx->local_bp, result->solutions);

  result->reset_seconds = monotonic_seconds() - phase_wall_start;

  phase_wall_start = monotonic_seconds();

  /*
   * get_index() opens and mmaps the index components. Install the immutable
   * pass advice before acquisition so every new mapping receives the correct
   * policy immediately.
   */
  fitsbin_mmap_set_thread_advice(
      result->mmap_advice);

  hook_result = index_shard_worker_get_index(
      ctx,
      shared,
      index_order,
      &index);
  hook_status = index_shard_apply_hook_result(
      result, hook_result, FALSE);

  fitsbin_mmap_clear_thread_advice();

  result->acquire_seconds =
      monotonic_seconds() - phase_wall_start;

  if (hook_status || !index) {
    if (!hook_status && !index) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }
    if (result->failed) {
      ERROR("Failed to load index order %zu", index_order);
    }
    if (index_shard_arbitrate_candidate(
            shared, index_order)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }

    index_shard_result_finish_task(
        result, shared, task_wall_start);
    return result->failed ? -1 : 0;
  }

  /*
   * Reapply to any component that the onefield hook reused rather than opened
   * during this acquisition.
   */
  {
    int advice_failures =
        index_shard_apply_index_mmap_advice(
            index,
            result->mmap_advice);

    if (advice_failures > 0) {
      __atomic_add_fetch(
          &shared->mmap_advice_failures,
          (unsigned long long)advice_failures,
          __ATOMIC_RELAXED);
      logerr("[index-shard] failed to apply mmap advice "
             "index_order=%zu components=%i\n",
             index_order,
             advice_failures);
    }
  }
  index_shard_inverse_cache_attach(
      ctx, index, &inverse_lease);

  result->acquire_seconds =
      monotonic_seconds() - phase_wall_start;

  /*
   * Another group can solve or exhaust the aggregate budget while this owner
   * is opening an index. Recheck before faulting solver payload.
   */
  if (index_shard_check_global_limits(shared) ||
      index_shard_master_stop_requested(shared)) {
    index_shard_master_limit_or_cancel_requested(
        shared,
        &result->hit_total_cpulimit,
        &result->hit_total_timelimit,
        &result->cancelled);

    phase_wall_start = monotonic_seconds();
    hook_result = index_shard_done_with_index(
        ctx,
        shared,
        index_order,
        index,
        &inverse_lease);
    (void)index_shard_apply_hook_result(
        result, hook_result, FALSE);
    index = NULL;
    result->release_seconds =
        monotonic_seconds() - phase_wall_start;

    if (index_shard_arbitrate_candidate(
            shared, index_order)) {
      index_shard_result_fail(
          result,
          INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
          -1);
    }
    index_shard_result_finish_task(
        result,
        shared,
        task_wall_start);
    return result->failed ? -1 : 0;
  }

  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] worker=%i start index_order=%zu index=%s\n",
           ctx->worker_id,
           index_order,
           index->indexname ? index->indexname : "(null)");
  }
  // time only the actual one-index solve section
  wall_start = monotonic_seconds();
  cpu_start = get_cpu_usage();

  // Worker-lifetime TLS lets onefield callbacks poll this pass for stop.
  ctx->current_index_order = index_order;
  ctx->current_outer_active = TRUE;
  hook_result = shared->hooks->solve_one_index(
      &ctx->local_bp, index);
  ctx->current_outer_active = FALSE;

  result->wall_seconds = monotonic_seconds() - wall_start;
  result->cpu_seconds = get_cpu_usage() - cpu_start;

  result->hit_total_cpulimit =
      ctx->local_bp.hit_total_cpulimit;
  result->hit_total_timelimit =
      ctx->local_bp.hit_total_timelimit;
  result->cancelled = ctx->local_bp.cancelled;
  result->solver_profile = ctx->local_bp.solver.profile;
  hook_status = index_shard_apply_hook_result(
      result, hook_result, FALSE);

  phase_wall_start = monotonic_seconds();
  if (!hook_status && !result->failed) {
    hook_status = index_shard_capture_solution_analysis(
        shared, result);
  }
  result->analyze_seconds =
      monotonic_seconds() - phase_wall_start;

  if (index_shard_arbitrate_candidate(
          shared, index_order)) {
    index_shard_result_fail(
        result,
        INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
        -1);
  }

  /*
   * Log the index name before done_with_index() releases index ownership.
   */
  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] worker=%i finish index_order=%zu index=%s\n",
           ctx->worker_id,
           index_order,
           index->indexname ? index->indexname : "(null)");
  }

  // Release through the original onefield ownership hook.
  phase_wall_start = monotonic_seconds();
  hook_result = index_shard_done_with_index(
      ctx,
      shared,
      index_order,
      index,
      &inverse_lease);
  (void)index_shard_apply_hook_result(
      result, hook_result, FALSE);
  index = NULL;

  result->release_seconds =
      monotonic_seconds() - phase_wall_start;

  index_shard_result_finish_task(
      result, shared, task_wall_start);

  return result->failed ? -1 : 0;
}
// ANCHOR INDEX-SHARD: worker-done
/*
 * Publish that this worker is done with the submitted pass.
 *
 * Reducer waits on active_workers reaching zero before final cleanup/drain.
 */
static void index_shard_worker_done(index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->result_mutex);

  if (shared->active_workers > 0)
    shared->active_workers--;

  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);
}

/*
 * A compute worker waiting for one mapped-page completion may execute one
 * already-published coarse helper task. This never claims an outer index and
 * never changes index ownership. The ticket remains responsible for waking
 * the worker when its own pages become ready.
 */
static int index_shard_payload_wait_stop(void *opaque) {
  index_shard_worker_context_t *ctx = opaque;

  if (!ctx || !ctx->pool) {
    return TRUE;
  }
  return index_shard_worker_stop_requested();
}

static int index_shard_payload_wait_help(void *opaque) {
  index_shard_worker_context_t *ctx = opaque;
  index_shard_thread_state_t *shared;
  index_shard_inner_claim_t claim;
  int selection;

  if (!ctx || !ctx->pool ||
      index_shard_worker_stop_requested()) {
    return FALSE;
  }
  shared = &ctx->pool->shared;
  memset(&claim, 0, sizeof(claim));

  pthread_mutex_lock(&shared->queue_mutex);
  if (ctx->published_helper_group) {
    index_shard_helper_group_t *group =
        ctx->published_helper_group;
    size_t reserved = 0U;

    if (group->foreign_reserve > group->foreign_claims) {
      reserved = group->foreign_reserve -
          group->foreign_claims;
    }
    if (!group->ready_count ||
        group->ready_count <= reserved) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return FALSE;
    }
    claim.kind = INDEX_SHARD_INNER_CLAIM_HELPER;
    claim.helper.group = group;
    selection = index_shard_helper_owner_claim_locked(
        shared, group, &claim.helper.task_index);
    if (!selection) {
      group->owner_claims++;
      group->owner_work += index_shard_helper_task_work(
          &group->tasks[claim.helper.task_index]);
      shared->helper_tasks_owner++;
    }
  } else {
    selection = index_shard_inner_select_locked(
        ctx, shared, TRUE, &claim);
  }
  pthread_mutex_unlock(&shared->queue_mutex);

  if (selection < 0) {
    index_shard_request_fatal_stop(shared);
    return FALSE;
  }
  if (selection > 0) {
    return FALSE;
  }
  if (index_shard_inner_execute_claim(shared, &claim)) {
    index_shard_request_fatal_stop(shared);
    return FALSE;
  }
  return TRUE;
}

// ANCHOR INDEX-SHARD: worker-main
/*
 * Persistent worker loop.
 *
 * Worker sleeps until pool generation changes, prepares local pass context,
 * claims one-index tasks, then cleans up local context when the pass ends.
 */
static void *index_shard_worker_main(void *userdata) {
  index_shard_worker_context_t *ctx = userdata;
  index_shard_pool_t *pool;
  int tls_status;

  if (!ctx) {
    return NULL;
  }

  pool = ctx->pool;
  if (!pool) {
    return NULL;
  }

  /*
   * Install the immutable worker context before declaring this thread ready.
   * A failed pthread TLS setup is a pool-start failure, never a worker that
   * silently runs without global stop/cancellation visibility.
   */
  tls_status = index_shard_set_tls(ctx);
  if (!tls_status) {
    tls_status = fitsbin_payload_io_set_thread_wait_helper(
        index_shard_payload_wait_help,
        index_shard_payload_wait_stop,
        ctx);
  }
  pthread_mutex_lock(&pool->control_mutex);
  if (tls_status && !pool->tls_startup_error) {
    pool->tls_startup_error = tls_status;
  }
  pool->ready_workers++;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);
  if (tls_status) {
    return NULL;
  }

  while (1) {
    index_shard_thread_state_t *shared = NULL;
    size_t index_order = 0U;
    int run_pass = FALSE;

    pthread_mutex_lock(&pool->control_mutex);

    while (!pool->shutdown &&
           ctx->generation_seen == pool->generation) {
      pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
    }

    if (pool->shutdown) {
      pthread_mutex_unlock(&pool->control_mutex);
      break;
    }

    ctx->generation_seen = pool->generation;
    shared = &pool->shared;

    if (ctx->worker_id < shared->worker_count) {
      run_pass = TRUE;
    }

    pthread_mutex_unlock(&pool->control_mutex);

    if (!run_pass) {
      continue;
    }

    ctx->pass_prepare_seconds = 0.0;
    ctx->pass_cleanup_seconds = 0.0;
    ctx->current_outer_active = FALSE;
    ctx->current_index_order = 0U;

    while (1) {
      index_shard_result_t *result = NULL;
      index_shard_inner_claim_t inner_claim;
      fitsbin_mmap_advice_t mmap_advice =
          FITSBIN_MMAP_ADVICE_NORMAL;
      index_shard_work_selection_t selection;

      memset(&inner_claim, 0, sizeof(inner_claim));
      selection = index_shard_select_work(
          ctx,
          shared,
          &index_order,
          &mmap_advice,
          &inner_claim);

      if (selection == INDEX_SHARD_WORK_ERROR) {
        index_shard_request_fatal_stop(shared);
        break;
      }
      if (selection == INDEX_SHARD_WORK_DONE) {
        break;
      }
      if (selection == INDEX_SHARD_WORK_HELPER) {
        if (index_shard_inner_execute_claim(
                shared, &inner_claim)) {
          index_shard_request_fatal_stop(shared);
          break;
        }
        continue;
      }
      if (selection != INDEX_SHARD_WORK_OUTER) {
        index_shard_request_fatal_stop(shared);
        break;
      }

      result = &shared->results[index_order];

      if (index_shard_check_global_limits(shared) ||
          index_shard_master_stop_requested(shared)) {
        index_shard_master_limit_or_cancel_requested(
            shared,
            &result->hit_total_cpulimit,
            &result->hit_total_timelimit,
            &result->cancelled);
        index_shard_finish_outer_claim(shared, index_order);
        break;
      }

      /*
       * Prepare worker-local solver and field state only after this worker has
       * acquired real outer work. This preserves pass-local reuse without
       * charging field read/preprocessing costs to idle workers when the
       * configured index set is smaller than the pool.
       */
      if (!ctx->local_context_ready) {
        double context_wall_start = monotonic_seconds();

        if (index_shard_worker_prepare_pass(ctx, shared)) {
          ctx->pass_prepare_seconds =
              monotonic_seconds() - context_wall_start;

          index_shard_result_init(result, index_order);
          index_shard_result_fail(
              result,
              INDEX_SHARD_FAILURE_GLOBAL_INTEGRITY,
              -1);

          index_shard_finish_outer_claim(shared, index_order);
          break;
        }

        ctx->pass_prepare_seconds =
            monotonic_seconds() - context_wall_start;
      }

      if (index_shard_run_one_with_worker_context(ctx,
                                                  shared,
                                                  index_order,
                                                  result,
                                                  mmap_advice)) {
        index_shard_finish_outer_claim(shared, index_order);
        if (result->failure_class == INDEX_SHARD_FAILURE_TASK_LOCAL) {
          continue;
        }
        break;
      }

      if (result->solved) {
        logverb("[index-shard] verified-result-ready worker=%i index_order=%zu "
                "best_logodds=%.3f field=%i wall=%.6f cpu=%.6f\n",
                ctx->worker_id,
                index_order,
                result->best_logodds,
                result->best_fieldnum,
                result->wall_seconds,
                (double)result->cpu_seconds);
      }

      if (index_shard_trace_enabled()) {
        logmsg("[index-shard] complete worker=%i index_order=%zu solved=%i "
               "failed=%i wall=%.6f cpu=%.6f pass_wall=%.6f\n",
               ctx->worker_id,
               index_order,
               result->solved,
               result->failed,
               result->wall_seconds,
               (double)result->cpu_seconds,
               monotonic_seconds() - shared->pass_wall_start);
      }

      index_shard_finish_outer_claim(shared, index_order);
      index_shard_check_global_limits(shared);
    }

    {
      double context_wall_start = monotonic_seconds();

      index_shard_worker_cleanup_pass(ctx, shared);
      ctx->pass_cleanup_seconds =
          monotonic_seconds() - context_wall_start;
    }

    index_shard_worker_done(shared);
  }

  fitsbin_payload_io_clear_thread_wait_helper();
  tls_status = index_shard_set_tls(NULL);
  if (tls_status) {
    logerr("[index-shard] failed to clear worker TLS "
           "worker=%i status=%i\n",
           ctx->worker_id,
           tls_status);
  }

  return NULL;
}

/*
 * SECTION INDEX-SHARD: reducer
 *
 * Main-thread result publication.
 *
 * Workers fill result slots.  The reducer is the only path that transfers
 * MatchObj data into master bp->solutions and updates final master solved state.
 */

// ANCHOR INDEX-SHARD: reduce-first-valid
/*
 * First-valid reducer for one submitted pass.
 *
 * Result publication selects the first immutable verified completion and
 * starts cooperative cancellation. After all workers quiesce, reduce exactly
 * that selected result. Configured index order is used only for the clean
 * unsolved drain that preserves original diagnostic accumulation.
 */
#define INDEX_SHARD_LIMIT_POLL_NANOSECONDS 100000000L

static void index_shard_limit_poll_deadline(
    struct timespec *deadline) {
  clock_gettime(CLOCK_MONOTONIC, deadline);
  deadline->tv_nsec +=
      INDEX_SHARD_LIMIT_POLL_NANOSECONDS;
  if (deadline->tv_nsec >= 1000000000L) {
    deadline->tv_sec++;
    deadline->tv_nsec -= 1000000000L;
  }
}

static int index_shard_pool_reduce_first_valid(index_shard_pool_t *pool) {
  index_shard_thread_state_t *shared = &pool->shared;
  index_shard_pass_state_snapshot_t state;
  size_t i;
  int rc = 0;

  /*
   * Winner selection already published stop. Keep all result storage and
   * shared field state alive until every owner and borrowed helper has left
   * this generation.
   */
  pthread_mutex_lock(&shared->result_mutex);
  while (shared->active_workers > 0) {
    struct timespec deadline;
    int wait_result;

    index_shard_limit_poll_deadline(&deadline);
    wait_result = pthread_cond_timedwait(
        &shared->result_cv,
        &shared->result_mutex,
        &deadline);
    if (wait_result == ETIMEDOUT) {
      pthread_mutex_unlock(&shared->result_mutex);
      (void)index_shard_check_global_limits(shared);
      pthread_mutex_lock(&shared->result_mutex);
    } else if (wait_result) {
      rc = -1;
      pthread_mutex_unlock(&shared->result_mutex);
      index_shard_request_fatal_stop(shared);
      pthread_mutex_lock(&shared->result_mutex);
    }
  }

  index_shard_pass_state_snapshot(shared, &state);
  pthread_mutex_unlock(&shared->result_mutex);

  if (state.fatal_error) {
    return -1;
  }
  if (rc) {
    return rc;
  }

  if (state.winner_selected) {
    index_shard_result_t *winner;

    if (state.selected_index_order >= shared->nindexes ||
        !shared->completed[state.selected_index_order]) {
      logerr("[index-shard] selected winner result is unavailable "
             "index_order=%zu\n",
             state.selected_index_order);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    winner = &shared->results[state.selected_index_order];
    if (!winner->candidate_ready ||
        !winner->solved ||
        winner->failed ||
        winner->rc ||
        winner->failure_class != INDEX_SHARD_FAILURE_NONE ||
        winner->cancelled ||
        winner->hit_total_timelimit ||
        winner->hit_total_cpulimit ||
        winner->merged ||
        !winner->candidate_sequence ||
        winner->candidate_sequence !=
            state.selected_candidate_sequence ||
        !winner->completion_sequence ||
        !shared->completed[state.selected_index_order]) {
      logerr("[index-shard] selected winner result is inconsistent "
             "index_order=%zu candidate_sequence=%zu\n",
             state.selected_index_order,
             state.selected_candidate_sequence);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (index_shard_trace_enabled()) {
      logmsg("[index-shard] reduce-winner index_order=%zu worker=%i "
             "candidate_sequence=%zu completion_sequence=%zu "
             "solved=%i failed=%i\n",
             winner->index_order,
             winner->worker_id,
             winner->candidate_sequence,
             winner->completion_sequence,
             winner->solved,
             winner->failed);
    }

    if (index_shard_reduce_one_result(shared, winner) ||
        !shared->have_committed_result ||
        shared->committed_index_order !=
            state.selected_index_order) {
      logerr("[index-shard] failed to commit selected winner "
             "index_order=%zu\n",
             state.selected_index_order);
      index_shard_request_fatal_stop(shared);
      return -1;
    }
    return 0;
  }

  /* Cancellation or a limit won before any immutable verified result. */
  if (state.stop_requested) {
    return 0;
  }

  if (state.task_local_failures) {
    logerr("[index-shard] pass exhausted with task-local failures "
           "count=%llu; requesting exact serial retry\n",
           state.task_local_failures);
    return -1;
  }

  /*
   * A clean unsolved pass preserves configured-order accumulation of
   * below-threshold diagnostics. No result in this path may be solved.
   */
  for (i = 0; i < shared->nindexes; i++) {
    index_shard_result_t *result = &shared->results[i];

    if (!shared->completed[i] ||
        result->failed ||
        result->rc ||
        result->failure_class != INDEX_SHARD_FAILURE_NONE ||
        result->solved) {
      logerr("[index-shard] invalid clean-unsolved result "
             "index_order=%zu completed=%i solved=%i failed=%i rc=%i\n",
             i,
             shared->completed[i] ? 1 : 0,
             result->solved,
             result->failed,
             result->rc);
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (index_shard_reduce_one_result(shared, result)) {
      index_shard_request_fatal_stop(shared);
      return -1;
    }

    if (shared->solved_published) {
      logerr("[index-shard] clean-unsolved reduction committed a solution "
             "index_order=%zu\n", i);
      index_shard_request_fatal_stop(shared);
      return -1;
    }
  }

  return 0;
}

/*
 * Report the reducer-owned winner after every worker has quiesced and the
 * selected result has committed. Worker callbacks deliberately suppress their
 * ordinary solved line, so this is the sole parallel success report.
 */
static int index_shard_report_committed_solution(
    onefield_t *bp,
    size_t nindexes,
    const index_shard_thread_state_t *shared,
    const index_shard_hooks_t *hooks,
    const index_shard_result_t *results) {
  const index_shard_result_t *committed;

  if (!hooks || !hooks->report_committed_solution) {
    logerr("[index-shard] committed-solution reporter is unavailable\n");
    return -1;
  }

  if (!shared ||
      !shared->winner_selected ||
      !shared->have_committed_result ||
      shared->committed_index_order >= nindexes ||
      shared->committed_index_order != shared->selected_index_order) {
    logerr("[index-shard] committed-solution identity is unavailable\n");
    return -1;
  }

  committed = &results[shared->committed_index_order];
  if (!committed->merged ||
      !committed->solved ||
      committed->failed ||
      committed->rc != 0) {
    logerr("[index-shard] committed-solution result is inconsistent\n");
    return -1;
  }

  return hooks->report_committed_solution(
      bp,
      committed->index_order,
      committed->best_fieldnum,
      committed->best_logodds);
}
/*
 * SECTION INDEX-SHARD: pool
 *
 * Pool lifecycle and pass submission.
 *
 * The pool is created once per engine job and reused across onefield_run()
 * calls.  Each submitted pass increments generation to wake workers.
 */

// ANCHOR INDEX-SHARD: shared-init
/*
 * Initialize synchronization primitives for the reusable shared pass state.
 */
static int index_shard_shared_init(index_shard_thread_state_t *shared) {
  pthread_condattr_t condattr;

  memset(shared, 0, sizeof(index_shard_thread_state_t));

  if (pthread_mutex_init(&shared->queue_mutex, NULL)) {
    return -1;
  }
  if (pthread_cond_init(&shared->queue_cv, NULL)) {
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_mutex_init(&shared->result_mutex, NULL)) {
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_condattr_init(&condattr)) {
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  if (pthread_condattr_setclock(
          &condattr,
          CLOCK_MONOTONIC)) {
    pthread_condattr_destroy(&condattr);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  if (pthread_cond_init(
          &shared->result_cv,
          &condattr)) {
    pthread_condattr_destroy(&condattr);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  pthread_condattr_destroy(&condattr);

  if (pthread_mutex_init(&shared->state_mutex, NULL)) {
    pthread_cond_destroy(&shared->result_cv);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  if (pthread_mutex_init(&shared->limit_mutex, NULL)) {
    pthread_mutex_destroy(&shared->state_mutex);
    pthread_cond_destroy(&shared->result_cv);
    pthread_mutex_destroy(&shared->result_mutex);
    pthread_cond_destroy(&shared->queue_cv);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }

  return 0;
}

// ANCHOR INDEX-SHARD: shared-destroy
/*
 * Destroy synchronization primitives after all workers have joined.
 */
static void index_shard_shared_destroy(index_shard_thread_state_t *shared) {
  pthread_cond_destroy(&shared->queue_cv);
  pthread_mutex_destroy(&shared->queue_mutex);

  pthread_mutex_destroy(&shared->result_mutex);
  pthread_cond_destroy(&shared->result_cv);

  pthread_mutex_destroy(&shared->state_mutex);

  pthread_mutex_destroy(&shared->limit_mutex);
}

typedef enum index_shard_pool_acquire_status {
  INDEX_SHARD_POOL_ACQUIRE_CONFLICT = -1,
  INDEX_SHARD_POOL_ACQUIRE_OK = 0,
  INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE = 1
} index_shard_pool_acquire_status_t;

/*
 * Reserve the persistent pool for one submitted pass.
 *
 * Lock order:
 *   index_shard_global_pool_mutex -> pool->control_mutex
 *
 * This prevents a second submission from overwriting shared pass pointers and
 * also keeps pool_stop() from destroying the pool while the caller still uses
 * it.
 */
static index_shard_pool_acquire_status_t
index_shard_pool_acquire_pass(onefield_t *bp,
                              solver_t *sp,
                              index_shard_pool_t **pool_out) {
  index_shard_pool_t *pool;

  if (!pool_out) {
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  *pool_out = NULL;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE;
  }

  if (pool->owner_bp != bp || pool->owner_sp != sp) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (pool->shutdown || pool->stopping || pool->pass_active) {
    pthread_mutex_unlock(&pool->control_mutex);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return INDEX_SHARD_POOL_ACQUIRE_CONFLICT;
  }

  pool->pass_active = TRUE;
  *pool_out = pool;

  pthread_mutex_unlock(&pool->control_mutex);
  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  return INDEX_SHARD_POOL_ACQUIRE_OK;
}

/*
 * Release one pass reservation after all caller-side work that dereferences
 * the pool has completed.
 */
static void index_shard_pool_release_pass(index_shard_pool_t *pool) {
  if (!pool) {
    return;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (!pool->pass_active) {
    logerr("[index-shard] pass release requested with no active pass\n");
  } else {
    pool->pass_active = FALSE;
  }

  /*
   * pool_stop() may be waiting for pass_active to become false. Worker waiters
   * also use work_cv, but they re-check their generation predicate in a loop.
   */
  pthread_cond_broadcast(&pool->work_cv);

  pthread_mutex_unlock(&pool->control_mutex);
}

// ANCHOR INDEX-SHARD: pool-start
/*
 * Create persistent worker pool.
 *
 * Workers are created once and sleep until the first pass is submitted.
 */
int index_shard_pool_start(onefield_t *bp, solver_t *sp) {
  index_shard_pool_t *pool;
  int i;
  int tls_status;
  int worker_count;
  int payload_io_lanes;
  int payload_io_width;

   // pool already active for this engine job
  if (!index_shard_pthread_enabled(bp)) {
    return 0;
  }

  if (!bp || !sp) {
    ERROR("Cannot start index-shard pool without owner state");
    return -1;
  }

  tls_status = index_shard_tls_ensure();
  if (tls_status) {
    logerr("[index-shard] cannot initialize worker TLS status=%i\n",
           tls_status);
    return -1;
  }

  /* Initialize process-wide lazy read-only state before workers start. */
  (void)fits_get_endian_string();

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  if (index_shard_global_pool) {
    int reusable;

    pool = index_shard_global_pool;

    if (pool->owner_bp != bp || pool->owner_sp != sp) {
      logerr("[index-shard] global pool already belongs to another engine job\n");
      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }

    pthread_mutex_lock(&pool->control_mutex);
    reusable = !pool->shutdown && !pool->stopping;
    pthread_mutex_unlock(&pool->control_mutex);

    pthread_mutex_unlock(&index_shard_global_pool_mutex);

    if (!reusable) {
      logerr("[index-shard] owner pool is stopping or shut down\n");
      return -1;
    }

    return 0;
  }

  worker_count = index_shard_get_worker_count(bp, 0);

  pool = calloc(1, sizeof(index_shard_pool_t));
  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    SYSERROR("Failed to allocate index-shard pool");
    return -1;
  }

  pool->owner_bp = bp;
  pool->owner_sp = sp;
  pool->worker_count = worker_count;
  pool->producer_width = 1U;
  pool->helper_width = 0U;
  pool->inverse_cache_budget =
      index_shard_inverse_cache_budget();

  if (pthread_mutex_init(&pool->inverse_cache_mutex, NULL)) {
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  // initialize shared state before workers can observe pool
  if (pthread_mutex_init(&pool->control_mutex, NULL)) {
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  if (pthread_cond_init(&pool->work_cv, NULL)) {
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  if (index_shard_shared_init(&pool->shared)) {
    pthread_cond_destroy(&pool->work_cv);
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  pool->threads = calloc((size_t)worker_count, sizeof(pthread_t));
  pool->contexts = calloc((size_t)worker_count,
                          sizeof(index_shard_worker_context_t));

    if (!pool->threads || !pool->contexts) {
      free(pool->threads);
      free(pool->contexts);

      index_shard_shared_destroy(&pool->shared);

      pthread_cond_destroy(&pool->work_cv);
      pthread_mutex_destroy(&pool->control_mutex);
      pthread_mutex_destroy(&pool->inverse_cache_mutex);

      free(pool);

      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }

  // worker contexts are stable for lifetime of the pool
  for (i = 0; i < worker_count; i++) {
    pool->contexts[i].worker_id = i;
    pool->contexts[i].generation_seen = 0;
    pool->contexts[i].pool = pool;

    if (pthread_create(&pool->threads[i], NULL, index_shard_worker_main, &pool->contexts[i])) {
      int j;

      pthread_mutex_lock(&pool->control_mutex);
      pool->shutdown = TRUE;
      pthread_cond_broadcast(&pool->work_cv);
      pthread_mutex_unlock(&pool->control_mutex);

            for (j = 0; j < i; j++) {
        pthread_join(pool->threads[j], NULL);
      }

       free(pool->threads);
      free(pool->contexts);

      index_shard_shared_destroy(&pool->shared);

      pthread_cond_destroy(&pool->work_cv);
      pthread_mutex_destroy(&pool->control_mutex);
      pthread_mutex_destroy(&pool->inverse_cache_mutex);

      free(pool);

      pthread_mutex_unlock(&index_shard_global_pool_mutex);
      return -1;
    }
  }

  /* Include worker scheduling/entry latency in pool startup, not first pass. */
  pthread_mutex_lock(&pool->control_mutex);
  while (pool->ready_workers < worker_count) {
    pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
  }
  tls_status = pool->tls_startup_error;
  if (tls_status) {
    pool->shutdown = TRUE;
    pthread_cond_broadcast(&pool->work_cv);
  }
  pthread_mutex_unlock(&pool->control_mutex);

  if (tls_status) {
    logerr("[index-shard] worker TLS startup failed status=%i\n",
           tls_status);
    for (i = 0; i < worker_count; i++) {
      pthread_join(pool->threads[i], NULL);
    }
    free(pool->threads);
    free(pool->contexts);
    index_shard_shared_destroy(&pool->shared);
    pthread_cond_destroy(&pool->work_cv);
    pthread_mutex_destroy(&pool->control_mutex);
    pthread_mutex_destroy(&pool->inverse_cache_mutex);
    free(pool);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return -1;
  }

  /*
   * Completion lanes populate only packet-planned mapped pages. They are not
   * compute workers and do not own solver, index, or result state.
   */
  payload_io_lanes = MIN(2, MAX(1, worker_count / 2));
  if (!fitsbin_payload_io_service_width()) {
    if (fitsbin_payload_io_service_start(payload_io_lanes)) {
      logverb("[index-shard] mapped-page completion unavailable; "
              "using native mmap demand\n");
    } else {
      pool->payload_io_owned = TRUE;
    }
  }
  if (fitsbin_payload_io_service_width() &&
      fitsbin_payload_io_mapped_population_supported() &&
      !fitsbin_payload_io_set_completion_notifier(
          index_shard_staged_completion_notify, pool)) {
    pool->payload_completion_registered = TRUE;
  } else if (fitsbin_payload_io_service_width() &&
             fitsbin_payload_io_mapped_population_supported()) {
    logverb("[index-shard] detached payload completion unavailable; "
            "retaining fixed helper reservation\n");
  } else if (fitsbin_payload_io_service_width()) {
    logverb("[index-shard] mapped-page population unsupported; "
            "retaining fixed helper reservation\n");
  }

  index_shard_global_pool = pool;
  payload_io_width = fitsbin_payload_io_service_width();
  if (pool->payload_completion_registered) {
    /*
     * Mapped-page delivery width and outer compute width are independent.
     * Every compute worker may own an outer index. A worker with no
     * immediately claimable outer work may still execute READY staged work
     * from any published owner.
     */
    pool->producer_width = (size_t)worker_count;
    pool->helper_width = 0U;
  } else {
    pool->helper_width = worker_count > 1 ? 1U : 0U;
    pool->producer_width =
        (size_t)worker_count - pool->helper_width;
  }
  fitsbin_payload_io_configure_workers(
      pool->payload_completion_registered
          ? (int)pool->producer_width
          : worker_count);

  logverb("[index-shard] workers=%i mode=pthread "
          "compute_width=%i producer_width=%zu helper_width=%zu "
          "payload_io_width=%i inverse_cache_budget=%zu\n",
          worker_count,
          worker_count,
          pool->producer_width,
          pool->helper_width,
          payload_io_width,
          pool->inverse_cache_budget);

  pthread_mutex_unlock(&index_shard_global_pool_mutex);
  return 0;
}
// ANCHOR INDEX-SHARD: pool-stop
/*
 * Stop the pool belonging to bp and join all workers.
 *
 * New passes are rejected as soon as stopping is published. If one pass is
 * already active, shutdown waits for that pass reservation to be released.
 */
void index_shard_pool_stop(onefield_t *bp) {
  index_shard_pool_t *pool;
  int i;
  int notifier_clear_failed = FALSE;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (!pool) {
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  if (pool->owner_bp != bp) {
    logerr("[index-shard] refusing to stop pool owned by another engine job\n");
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (pool->stopping) {
    pthread_mutex_unlock(&pool->control_mutex);
    pthread_mutex_unlock(&index_shard_global_pool_mutex);
    return;
  }

  pool->stopping = TRUE;
  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  while (pool->pass_active) {
    pthread_cond_wait(&pool->work_cv, &pool->control_mutex);
  }

  pool->shutdown = TRUE;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  if (index_shard_global_pool == pool) {
    index_shard_global_pool = NULL;
  }

  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  for (i = 0; i < pool->worker_count; i++) {
    pthread_join(pool->threads[i], NULL);
  }

  if (pool->payload_completion_registered) {
    if (fitsbin_payload_io_clear_completion_notifier(
            index_shard_staged_completion_notify, pool)) {
      logerr("[index-shard] failed to clear payload completion notifier\n");
      notifier_clear_failed = TRUE;
    } else {
      pool->payload_completion_registered = FALSE;
    }
  }
  if (pool->payload_io_owned) {
    fitsbin_payload_io_service_stop();
  }
  if (notifier_clear_failed) {
    /*
     * The process-wide notifier may still retain pool as opaque state, or an
     * external clear may still be draining an active callback. Preserve the
     * stopped allocation rather than risk a callback use-after-free.
     */
    logerr("[index-shard] retaining stopped pool after notifier "
           "teardown failure\n");
    return;
  }
  free(pool->threads);
  free(pool->contexts);
  index_shard_shared_destroy(&pool->shared);

  logverb("[index-shard] inverse-cache hits=%llu misses=%llu "
          "admitted=%llu refused=%llu evicted=%llu "
          "overcommit=%llu retained=%zu active=%zu "
          "retained_peak=%zu combined_peak=%zu budget=%zu\n",
          pool->inverse_cache_hits,
          pool->inverse_cache_misses,
          pool->inverse_cache_admitted,
          pool->inverse_cache_refused,
          pool->inverse_cache_evicted,
          pool->inverse_cache_overcommit,
          pool->inverse_cache_bytes,
          pool->inverse_active_bytes,
          pool->inverse_cache_peak_bytes,
          pool->inverse_combined_peak_bytes,
          pool->inverse_cache_budget);
  index_shard_inverse_cache_destroy(pool);
  pthread_mutex_destroy(&pool->inverse_cache_mutex);
  pthread_cond_destroy(&pool->work_cv);
  pthread_mutex_destroy(&pool->control_mutex);

  logverb("[index-shard] pthread-pool stop\n");

  free(pool);
}

// ANCHOR INDEX-SHARD: pool-active
/*
 * Return true only when the compatibility pool belongs to bp and remains
 * available for a new pass reservation.
 */
int index_shard_pool_active(onefield_t *bp) {
  index_shard_pool_t *pool;
  int active = FALSE;

  pthread_mutex_lock(&index_shard_global_pool_mutex);

  pool = index_shard_global_pool;

  if (pool && pool->owner_bp == bp) {
    pthread_mutex_lock(&pool->control_mutex);
    active = !pool->shutdown && !pool->stopping;
    pthread_mutex_unlock(&pool->control_mutex);
  }

  pthread_mutex_unlock(&index_shard_global_pool_mutex);

  return active;
}

// ANCHOR INDEX-SHARD: pool-submit
/*
 * Submit one onefield_run() pass to the persistent pool.
 *
 * This resets shared pass state, publishes result arrays, then increments
 * generation to wake workers.
 */
static int index_shard_pool_submit(
    index_shard_pool_t *pool,
    onefield_t *bp,
    solver_t *base_sp,
    size_t nindexes,
    const index_shard_hooks_t *hooks,
    const void *worker_view,
    index_shard_result_t *results,
    unsigned char *completed,
    unsigned char *outer_states) {
  index_shard_thread_state_t *shared = &pool->shared;
  int worker_count = pool->worker_count;
  int i;
  int payload_io_width;

  if (!outer_states || !worker_view || !pool->producer_width) {
    return -1;
  }

  pthread_mutex_lock(&pool->control_mutex);

  if (!pool->pass_active || pool->shutdown ||
      pool->owner_bp != bp || pool->owner_sp != base_sp) {
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  pthread_mutex_lock(&shared->queue_mutex);
  for (i = 0; i < worker_count; i++) {
    if (pool->contexts[i].published_helper_group ||
        pool->contexts[i].published_staged_group ||
        pool->contexts[i].staged_owner_callback_active ||
        pool->contexts[i].staged_owner_callback_group ||
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] inner state remained active "
             "before pass worker=%i\n", i);
      pthread_mutex_unlock(&shared->queue_mutex);
      pthread_mutex_unlock(&pool->control_mutex);
      return -1;
    }
  }
  if (shared->queue_waiters) {
    logerr("[index-shard] queue waiters remained before pass "
           "count=%zu\n", shared->queue_waiters);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->helper_groups_published !=
          shared->helper_groups_completed ||
      shared->staged_groups_published !=
          shared->staged_groups_completed) {
    logerr("[index-shard] inner group lifecycle remained "
           "before pass helper=%llu/%llu staged=%llu/%llu\n",
           shared->helper_groups_published,
           shared->helper_groups_completed,
           shared->staged_groups_published,
           shared->staged_groups_completed);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->helper_groups_active ||
      shared->helper_preparations_active ||
      shared->helper_foreign_reservations ||
      shared->staged_groups_active ||
      shared->staged_tickets_active ||
      shared->staged_source_leases ||
      shared->staged_compute_ready ||
      shared->staged_reorder_ready) {
    logerr("[index-shard] inner activity remained before pass "
           "helper_groups=%zu preparations=%zu reservations=%zu "
           "staged_groups=%zu tickets=%zu leases=%zu "
           "compute_ready=%zu reorder_ready=%zu\n",
           shared->helper_groups_active,
           shared->helper_preparations_active,
           shared->helper_foreign_reservations,
           shared->staged_groups_active,
           shared->staged_tickets_active,
           shared->staged_source_leases,
           shared->staged_compute_ready,
           shared->staged_reorder_ready);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  pthread_mutex_lock(&shared->result_mutex);
  pthread_mutex_lock(&shared->state_mutex);
  pthread_mutex_lock(&shared->limit_mutex);

  shared->bp = bp;
  shared->base_sp = base_sp;
  shared->hooks = hooks;
  shared->worker_view = worker_view;
  shared->nindexes = nindexes;
  shared->canonical_scan_cursor = 0U;
  shared->outer_unclaimed = nindexes;
  shared->outer_running = 0U;
  shared->producer_width = MIN(pool->producer_width, nindexes);
  /* Idle outer capacity joins the inner width for this pass. */
  shared->helper_width =
      (size_t)worker_count - shared->producer_width;
  shared->queue_waiters = 0U;
  shared->helper_groups_active = 0U;
  shared->helper_preparations_active = 0U;
  shared->helper_foreign_reservations = 0U;
  shared->staged_groups_active = 0U;
  shared->staged_tickets_active = 0U;
  shared->staged_source_leases = 0U;
  shared->staged_compute_ready = 0U;
  shared->staged_reorder_ready = 0U;
  shared->staged_max_compute_running = 0U;
  shared->staged_completion_epoch = 1U;
  shared->outer_states = outer_states;
  shared->outer_claims = 0U;
  shared->helper_groups_published = 0U;
  shared->helper_groups_completed = 0U;
  shared->helper_tasks_owner = 0U;
  shared->helper_tasks_foreign = 0U;
  shared->helper_task_failures = 0U;
  shared->helper_owner_wait_calls = 0U;
  shared->helper_owner_wait_seconds = 0.0;
  shared->staged_groups_published = 0U;
  shared->staged_groups_completed = 0U;
  shared->staged_tasks_owner = 0U;
  shared->staged_tasks_foreign = 0U;
  shared->staged_compute_owner = 0U;
  shared->staged_compute_foreign = 0U;
  shared->staged_task_failures = 0U;
  shared->staged_io_submitted = 0U;
  shared->staged_io_completed = 0U;
  shared->staged_submit_retries = 0U;
  shared->staged_owner_wait_calls = 0U;
  shared->staged_owner_wait_seconds = 0.0;
  shared->staged_max_io_submitted = 0U;
  shared->staged_max_compute_ready = 0U;
  shared->staged_max_reorder_ready = 0U;
  shared->staged_prepare_claims = 0U;
  shared->staged_submit_claims = 0U;
  shared->staged_poll_claims = 0U;
  shared->staged_execute_claims = 0U;
  shared->staged_owner_execute_claims = 0U;
  shared->staged_submit_to_ready_seconds = 0.0;
  shared->staged_ready_dwell_seconds = 0.0;
  shared->staged_execute_seconds = 0.0;
  shared->staged_result_to_retire_seconds = 0.0;
  shared->staged_retire_seconds = 0.0;
  shared->task_local_failures = 0U;
  shared->global_integrity_failures = 0U;
  shared->late_loser_failures = 0U;

  shared->results = results;
  shared->completed = completed;
  shared->results_reduced = 0U;
  shared->next_completion_sequence = 0U;
  shared->next_candidate_sequence = 0U;

  shared->worker_count = worker_count;
  shared->active_workers = worker_count;

  shared->stop_requested = FALSE;
  shared->fatal_error = FALSE;
  shared->winner_selected = FALSE;
  shared->solved_published = FALSE;
  shared->master_committed = FALSE;
  shared->terminal_cause = INDEX_SHARD_TERMINAL_NONE;
  shared->first_stop_wall_since_pass = -1.0;
  __atomic_store_n(
      &shared->worker_stop_requested,
      FALSE,
      __ATOMIC_RELEASE);
  shared->selected_index_order = nindexes;
  shared->selected_candidate_sequence = 0U;
  shared->have_committed_result = FALSE;
  shared->committed_index_order = 0U;
  shared->limit_reported = FALSE;

  shared->reducer_work_calls = 0U;
  shared->reducer_work_wall_seconds = 0.0;

  /*
   * Preserve the established mmap policy. Coarse CodeKD packets may populate
   * a bounded exact DATA/PERM plan, but there is no independent delivery
   * generation and native mmap demand remains authoritative.
   */
  shared->mmap_advice =
      fitsbin_mmap_advice_state_begin_pass(
          &base_sp->index_mmap_policy);
  shared->mmap_pass_number =
      base_sp->index_mmap_policy.pass_number;
  shared->mmap_advice_failures = 0U;

  // Pass timing excludes final solve-field output generation.
  shared->pass_wall_start = monotonic_seconds();
  shared->pass_cpu_start = get_cpu_usage();
  shared->pass_rusage_valid =
      (getrusage(RUSAGE_SELF,
                 &shared->pass_rusage_start) == 0);

  pthread_mutex_unlock(&shared->limit_mutex);
  pthread_mutex_unlock(&shared->state_mutex);
  pthread_mutex_unlock(&shared->result_mutex);
  pthread_mutex_unlock(&shared->queue_mutex);

  // Generation publication is the hard band barrier release.
  pool->generation++;
  pthread_cond_broadcast(&pool->work_cv);
  pthread_mutex_unlock(&pool->control_mutex);
  payload_io_width = fitsbin_payload_io_service_width();

  logverb("[index-shard] pthread-pool submit compute_width=%i "
          "producer_width=%zu helper_width=%zu "
          "candidates=%zu engine_pass=%zu "
          "depth_index=%zu scale_index=%zu startobj=%i endobj=%i "
          "scheduler=first-valid-completion inner_scheduler=ordered-codekd-packets "
          "mmap_pass=%u mmap_advice=%s "
          "mmap_scope=payload-random-topology-normal "
          "mmap_policy=parallel-random-serial-normal "
          "page_delivery=%s payload_io=%s payload_io_width=%i "
          "outer_admission=%s\n",
          worker_count,
          shared->producer_width,
          shared->helper_width,
          nindexes,
          bp->engine_pass_ordinal,
          bp->engine_depth_index,
          bp->engine_scale_index,
          base_sp->startobj,
          base_sp->endobj,
          shared->mmap_pass_number,
          fitsbin_mmap_advice_name(shared->mmap_advice),
          pool->payload_completion_registered
              ? "detached-bounded-mapped-completion"
              : "native-mmap-demand",
          payload_io_width > 0 ? "kernel-page-cache" : "native-mmap",
          payload_io_width,
          pool->payload_completion_registered
              ? (shared->producer_width ==
                         (size_t)worker_count
                     ? "full-producer"
                     : "bounded-delivery-window")
              : "reserved-helper");

  return 0;
}

/*
 * SECTION INDEX-SHARD: entry
 */

// ANCHOR INDEX-SHARD: entry
/*
 * Execute one complete index-shard pass.
 *
 * Terminal status is classified by failure scope and the master mutation
 * boundary. Only an unavailable path or an isolated task-local failure proven
 * to occur before master commit may return control to the serial path. A
 * global-integrity failure is terminal regardless of winner timing.
 */
static index_shard_solve_status_t
index_shard_solve_impl(onefield_t *bp,
                       solver_t *base_sp,
                       size_t nindexes,
                       const index_shard_hooks_t *hooks) {
  index_shard_pool_t *pool;
  unsigned char *outer_states = NULL;
  index_shard_result_t *results = NULL;
  unsigned char *completed = NULL;
  void *worker_view = NULL;
  size_t i;
  int acquire_rc;
  int rc = 0;
  index_shard_pass_state_snapshot_t state;
  index_shard_pass_metrics_snapshot_t pass_metrics;
  index_shard_task_profile_snapshot_t task_profile;
  index_shard_phase_profile_snapshot_t phase_profile;
  solver_profile_t solver_profile;
  index_shard_solve_status_t status = INDEX_SHARD_SOLVE_HANDLED;
  double context_cleanup_max_seconds = 0.0;
  double context_cleanup_wall_seconds = 0.0;
  double context_prepare_max_seconds = 0.0;
  double context_prepare_wall_seconds = 0.0;
  double reduction_ex_verify_seconds;
  double stop_to_quiescence_seconds = 0.0;
  anbool pass_completed;
  anbool pass_exhaustive;
  anbool pass_solved;
  anbool pass_cancelled;
  anbool clean_exhaustion_required;
  anbool helper_quiescence_valid = TRUE;
  anbool mmap_transitioned;
  anbool worker_cancelled = FALSE;

  if (!index_shard_pthread_enabled(bp)) {
    return INDEX_SHARD_SOLVE_UNAVAILABLE;
  }

  // no candidate indexes, nothing to do
  if (!nindexes) {
    return INDEX_SHARD_SOLVE_HANDLED;
  }

  if (!hooks) {
    ERROR("index-shard hooks are NULL");
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  acquire_rc = index_shard_pool_acquire_pass(bp, base_sp, &pool);

  if (acquire_rc == INDEX_SHARD_POOL_ACQUIRE_UNAVAILABLE) {
    logverb("[index-shard] pthread mode requested but pool inactive\n");
    return INDEX_SHARD_SOLVE_UNAVAILABLE;
  }

  if (acquire_rc != INDEX_SHARD_POOL_ACQUIRE_OK) {
    logerr("[index-shard] pool ownership or pass lifecycle conflict; "
           "serial fallback suppressed\n");
    return INDEX_SHARD_SOLVE_LIFECYCLE_CONFLICT;
  }

  if (!hooks->create_worker_view ||
      !hooks->destroy_worker_view ||
      !hooks->prepare_local_context ||
      hooks->create_worker_view(
          bp, base_sp, &worker_view) ||
      !worker_view) {
    logerr("[index-shard] failed to create immutable worker view\n");
    if (worker_view && hooks->destroy_worker_view) {
      hooks->destroy_worker_view(worker_view);
    }
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  results = calloc(nindexes, sizeof(index_shard_result_t));
  completed = calloc(nindexes, sizeof(unsigned char));
  outer_states = calloc(nindexes, sizeof(*outer_states));

  if (!results || !completed || !outer_states) {
    SYSERROR("Failed to allocate index-shard pass state");

    free(outer_states);
    free(results);
    free(completed);

    hooks->destroy_worker_view(worker_view);
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  // Submit releases the hard current-band barrier to persistent workers.
  rc = index_shard_pool_submit(
      pool,
      bp,
      base_sp,
      nindexes,
      hooks,
      worker_view,
      results,
      completed,
      outer_states);

  if (rc) {
    free(outer_states);
    free(results);
    free(completed);

    hooks->destroy_worker_view(worker_view);
    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  rc = index_shard_pool_reduce_first_valid(pool);

  index_shard_pass_state_snapshot(&pool->shared, &state);
  clean_exhaustion_required =
      !state.solved_published &&
      !state.stop_requested &&
      !state.fatal_error;
  pthread_mutex_lock(&pool->shared.queue_mutex);
  for (i = 0U; i < (size_t)pool->shared.worker_count; i++) {
    if (pool->contexts[i].published_helper_group ||
        pool->contexts[i].published_staged_group ||
        pool->contexts[i].staged_owner_callback_active ||
        pool->contexts[i].staged_owner_callback_group ||
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] inner state remained active "
             "after worker quiescence worker=%zu\n", i);
      rc = -1;
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
      helper_quiescence_valid = FALSE;
    }
  }
  if (pool->shared.queue_waiters) {
    logerr("[index-shard] queue waiters remained after "
           "worker quiescence count=%zu\n",
           pool->shared.queue_waiters);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_groups_active) {
    logerr("[index-shard] helper groups remained active after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_groups_active);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_preparations_active) {
    logerr("[index-shard] helper preparations remained active after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_preparations_active);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_foreign_reservations) {
    logerr("[index-shard] helper reservations remained after "
           "worker quiescence count=%zu\n",
           pool->shared.helper_foreign_reservations);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.helper_groups_published !=
      pool->shared.helper_groups_completed) {
    logerr("[index-shard] helper group lifecycle mismatch "
           "after worker quiescence published=%llu completed=%llu\n",
           pool->shared.helper_groups_published,
           pool->shared.helper_groups_completed);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.staged_groups_active ||
      pool->shared.staged_tickets_active ||
      pool->shared.staged_source_leases ||
      pool->shared.staged_compute_ready ||
      pool->shared.staged_reorder_ready) {
    logerr("[index-shard] staged activity remained after "
           "worker quiescence groups=%zu tickets=%zu leases=%zu "
           "compute_ready=%zu reorder_ready=%zu\n",
           pool->shared.staged_groups_active,
           pool->shared.staged_tickets_active,
           pool->shared.staged_source_leases,
           pool->shared.staged_compute_ready,
           pool->shared.staged_reorder_ready);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.staged_groups_published !=
      pool->shared.staged_groups_completed) {
    logerr("[index-shard] staged group lifecycle mismatch "
           "after worker quiescence published=%llu completed=%llu\n",
           pool->shared.staged_groups_published,
           pool->shared.staged_groups_completed);
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    helper_quiescence_valid = FALSE;
  }
  if (pool->shared.outer_running) {
    logerr("[index-shard] pass ended with %zu outer owners\n",
           pool->shared.outer_running);
    rc = -1;
  }
  if (clean_exhaustion_required &&
      (pool->shared.outer_unclaimed ||
       pool->shared.canonical_scan_cursor != nindexes)) {
    logerr("[index-shard] clean pass ended before claim exhaustion "
           "unclaimed=%zu cursor=%zu candidates=%zu\n",
           pool->shared.outer_unclaimed,
           pool->shared.canonical_scan_cursor,
           nindexes);
    rc = -1;
  }
  for (i = 0U; i < nindexes; i++) {
    if (outer_states[i] == INDEX_SHARD_OUTER_RUNNING) {
      logerr("[index-shard] pass ended with running index_order=%zu\n",
             i);
      rc = -1;
      break;
    }
    if (clean_exhaustion_required &&
        outer_states[i] != INDEX_SHARD_OUTER_FINISHED) {
      logerr("[index-shard] clean pass ended with unfinished "
             "index_order=%zu state=%u\n",
             i,
             (unsigned int)outer_states[i]);
      rc = -1;
      break;
    }
  }
  pthread_mutex_unlock(&pool->shared.queue_mutex);

  logverb("[index-shard] ownership-pass generation=%lu "
          "scale_index=%zu canonical_claims=%llu "
          "producer_width=%zu helper_width=%zu "
          "unclaimed=%zu quiescent=%i\n",
          pool->generation,
          bp->engine_scale_index,
          pool->shared.outer_claims,
          pool->shared.producer_width,
          pool->shared.helper_width,
          pool->shared.outer_unclaimed,
          helper_quiescence_valid ? 1 : 0);

  memset(&solver_profile, 0, sizeof(solver_profile));

  /*
   * Worker-local cancellation is published only after every result slot is
   * immutable. Master onefield state remains reducer/caller-owned, while the
   * limit mutex preserves the existing synchronization discipline for the
   * cancellation flag.
   */
  for (i = 0; i < nindexes; i++) {
    if (results[i].task_started) {
      solver_profile_accumulate(
          &solver_profile,
          &results[i].solver_profile);
    }

    if (results[i].cancelled) {
      worker_cancelled = TRUE;
    }
  }

  if (worker_cancelled && !state.winner_selected) {
    pthread_mutex_lock(&pool->shared.limit_mutex);
    bp->cancelled = TRUE;
    pthread_mutex_unlock(&pool->shared.limit_mutex);
  }

  /*
   * active_workers reached zero before the first-valid reducer returned, so every
   * participating context timing is immutable here. These are work sums; the
   * maxima expose the critical per-worker prepare/cleanup contribution.
   */
  for (i = 0; i < (size_t)pool->shared.worker_count; i++) {
    double cleanup_seconds = pool->contexts[i].pass_cleanup_seconds;
    double prepare_seconds = pool->contexts[i].pass_prepare_seconds;

    context_cleanup_wall_seconds += cleanup_seconds;
    context_prepare_wall_seconds += prepare_seconds;

    if (cleanup_seconds > context_cleanup_max_seconds) {
      context_cleanup_max_seconds = cleanup_seconds;
    }

    if (prepare_seconds > context_prepare_max_seconds) {
      context_prepare_max_seconds = prepare_seconds;
    }
  }

  index_shard_pass_state_snapshot(&pool->shared, &state);

  reduction_ex_verify_seconds =
      solver_profile.resolve_wall_seconds -
      solver_profile.verify_wall_seconds;

  if (reduction_ex_verify_seconds < 0.0) {
    reduction_ex_verify_seconds = 0.0;
  }

  if (state.fatal_error || state.global_integrity_failures) {
    rc = -1;
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
  } else if (rc &&
             status != INDEX_SHARD_SOLVE_TERMINAL_FAILURE) {
    if (state.master_committed) {
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    } else {
      status = INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
    }
  }

   /*
   * Defensive ownership invariant: once the reducer has crossed the master
   * mutation boundary, a precommit failure classification is impossible.
   */
  if (state.master_committed &&
      status == INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE) {
    logerr("[index-shard] invalid precommit status after master commit; "
           "promoting to terminal failure\n");
    status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
  }

   /*
   * Workers have left the pass and result slots are now immutable. Snapshot
   * pass timing and task-duration distribution before destroying results.
   */
  index_shard_pass_metrics_snapshot(&pool->shared,
                                    &pass_metrics);

  if (state.first_stop_wall_since_pass >= 0.0 &&
      pass_metrics.wall_seconds >= state.first_stop_wall_since_pass) {
    stop_to_quiescence_seconds =
        pass_metrics.wall_seconds - state.first_stop_wall_since_pass;
  }

  index_shard_task_profile_snapshot(results,
                                    nindexes,
                                    pass_metrics.wall_seconds,
                                    &task_profile);

  index_shard_phase_profile_snapshot(results,
                                     nindexes,
                                     &phase_profile);


  /*
   * index_shard_pool_reduce_first_valid() returned and every participating
   * worker has left this generation. The pass outcome is now immutable.
   */
  pass_completed =
      task_profile.executed == nindexes;

  pass_exhaustive =
      pass_completed &&
      pass_metrics.reduced == nindexes;

  pass_solved =
      state.solved_published ||
      bp->single_field_solved;

  pass_cancelled =
      !pass_solved &&
      (bp->cancelled ||
       bp->hit_total_cpulimit ||
       bp->hit_total_timelimit ||
       state.stop_requested);

  if (status == INDEX_SHARD_SOLVE_HANDLED &&
      rc == 0 &&
      state.solved_published) {
    if (index_shard_report_committed_solution(
            bp,
            nindexes,
            &pool->shared,
            hooks,
            results)) {
      rc = -1;
      status = INDEX_SHARD_SOLVE_TERMINAL_FAILURE;
    }
  }

  mmap_transitioned =
      fitsbin_mmap_policy_complete_pass(
          &base_sp->index_mmap_policy,
          pass_completed,
          pass_exhaustive,
          pass_solved,
          pass_cancelled,
          rc,
          (int)status);
  // dispose unmerged worker results after all workers have left pass
  for (i = 0; i < nindexes; i++) {
    index_shard_result_dispose(&results[i], hooks);
  }

  free(results);
  free(completed);

  logverb("[index-shard] done workers=%i solved=%i "
          "wall=%.6f cpu=%.6f utilization=%.1f%% "
          "effective_concurrency=%.2f stop_to_quiescence=%.6f\n",
          pool->worker_count,
          bp->single_field_solved,
          pass_metrics.wall_seconds,
          (double)pass_metrics.cpu_seconds,
          pass_metrics.cpu_percent,
          pass_metrics.cpu_percent / 100.0,
          stop_to_quiescence_seconds);

  logverb("[index-shard] reducer-pass generation=%lu candidates=%zu "
          "calls=%llu work_wall_sum=%.6f\n",
          pool->generation,
          nindexes,
          pool->shared.reducer_work_calls,
          pool->shared.reducer_work_wall_seconds);

  logverb("[index-shard] helper-pass generation=%lu "
          "groups=%llu completed=%llu owner_tasks=%llu "
          "foreign_tasks=%llu task_failures=%llu "
          "owner_waits=%llu owner_wait_seconds=%.6f\n",
          pool->generation,
          pool->shared.helper_groups_published,
          pool->shared.helper_groups_completed,
          pool->shared.helper_tasks_owner,
          pool->shared.helper_tasks_foreign,
          pool->shared.helper_task_failures,
          pool->shared.helper_owner_wait_calls,
          pool->shared.helper_owner_wait_seconds);

  logverb("[index-shard] staged-pass generation=%lu "
          "groups=%llu completed=%llu owner_claims=%llu "
          "foreign_claims=%llu owner_compute=%llu "
          "foreign_compute=%llu task_failures=%llu "
          "io_submitted=%llu io_completed=%llu submit_retries=%llu "
          "max_io_submitted=%zu max_compute_ready=%zu "
          "max_reorder_ready=%zu max_compute_running=%zu "
          "prepare_claims=%llu "
          "submit_claims=%llu poll_claims=%llu execute_claims=%llu "
          "owner_execute_claims=%llu submit_to_ready_sum=%.6f "
          "ready_dwell_sum=%.6f execute_sum=%.6f "
          "result_to_retire_sum=%.6f retire_sum=%.6f "
          "owner_waits=%llu owner_wait_seconds=%.6f\n",
          pool->generation,
          pool->shared.staged_groups_published,
          pool->shared.staged_groups_completed,
          pool->shared.staged_tasks_owner,
          pool->shared.staged_tasks_foreign,
          pool->shared.staged_compute_owner,
          pool->shared.staged_compute_foreign,
          pool->shared.staged_task_failures,
          pool->shared.staged_io_submitted,
          pool->shared.staged_io_completed,
          pool->shared.staged_submit_retries,
          pool->shared.staged_max_io_submitted,
          pool->shared.staged_max_compute_ready,
          pool->shared.staged_max_reorder_ready,
          pool->shared.staged_max_compute_running,
          pool->shared.staged_prepare_claims,
          pool->shared.staged_submit_claims,
          pool->shared.staged_poll_claims,
          pool->shared.staged_execute_claims,
          pool->shared.staged_owner_execute_claims,
          pool->shared.staged_submit_to_ready_seconds,
          pool->shared.staged_ready_dwell_seconds,
          pool->shared.staged_execute_seconds,
          pool->shared.staged_result_to_retire_seconds,
          pool->shared.staged_retire_seconds,
          pool->shared.staged_owner_wait_calls,
          pool->shared.staged_owner_wait_seconds);

  logverb("[index-shard] context-pass generation=%lu candidates=%zu "
          "compute_width=%i producer_width=%zu helper_width=%zu "
          "prepare_work_wall_sum=%.6f prepare_max=%.6f "
          "cleanup_work_wall_sum=%.6f cleanup_max=%.6f\n",
          pool->generation,
          nindexes,
          pool->shared.worker_count,
          pool->shared.producer_width,
          pool->shared.helper_width,
          context_prepare_wall_seconds,
          context_prepare_max_seconds,
          context_cleanup_wall_seconds,
          context_cleanup_max_seconds);

  logverb("[index-shard] solver-pass generation=%lu candidates=%zu "
          "detailed=%i failed=%i solver_run_work_wall_sum=%.6f "
          "codekd_work_wall_sum=%.6f codekd_calls=%llu "
          "codekd_hits=%llu "
          "resolve_work_wall_sum=%.6f "
          "reduction_ex_verify_hit_work_wall_sum=%.6f "
          "resolve_calls=%llu verify_hit_work_wall_sum=%.6f "
          "verify_calls=%llu "
          "batches=%llu completed=%llu stopped=%llu "
          "batch_failed=%llu hypotheses=%llu executed=%llu reduced=%llu "
          "task_ranges=%llu tasks_executed=%llu submitted=%llu "
          "inline=%llu parallel_batches=%llu observed_parallel=%llu "
          "parallel_hypotheses=%llu alloc_failures=%llu "
          "search_failures=%llu max_batch=%zu max_tasks=%zu "
          "max_parallel=%zu helper_tasks=%llu "
          "helper_combinations=%llu "
          "hypothesis_order=%016llx kd_result_order=%016llx "
          "candidate_order=%016llx\n",
          pool->generation,
          nindexes,
          solver_profile.detailed ? 1 : 0,
          solver_profile.execution_failed ? 1 : 0,
          solver_profile.solver_run_wall_seconds,
          solver_profile.codekd_wall_seconds,
          solver_profile.codekd_calls,
          solver_profile.codekd_hits,
          solver_profile.resolve_wall_seconds,
          reduction_ex_verify_seconds,
          solver_profile.resolve_calls,
          solver_profile.verify_wall_seconds,
          solver_profile.verify_calls,
          solver_profile.hypothesis_batches,
          solver_profile.hypothesis_batches_completed,
          solver_profile.hypothesis_batches_stopped,
          solver_profile.hypothesis_batches_failed,
          solver_profile.hypotheses_generated,
          solver_profile.hypotheses_executed,
          solver_profile.hypotheses_reduced,
          solver_profile.task_ranges_planned,
          solver_profile.task_ranges_executed,
          solver_profile.task_ranges_submitted,
          solver_profile.task_ranges_inline,
          solver_profile.parallel_batches,
          solver_profile.parallel_batches_observed,
          solver_profile.parallel_hypotheses,
          solver_profile.allocation_failures,
          solver_profile.search_failures,
          solver_profile.max_batch_hypotheses,
          solver_profile.max_task_ranges,
          solver_profile.max_parallel_ranges,
          solver_profile.ab_helper_tasks,
          solver_profile.ab_helper_combinations,
          solver_profile.hypothesis_order_hash,
          solver_profile.kd_result_order_hash,
          solver_profile.candidate_order_hash);

  logverb("[index-shard] page-pipeline generation=%lu "
          "descriptors=%llu complete=%llu boundary_deferrals=%llu "
          "raw_hints=%llu unique_pages=%llu coalesced_ranges=%llu "
          "logical_bytes=%llu aligned_bytes=%llu overread_bytes=%llu "
          "refusals=%llu/%llu/%llu/%llu/%llu/%llu/%llu/%llu/%llu "
          "candidate_delivery=%llu quad=%llu/%llu/%llu "
          "star=%llu/%llu/%llu windows=%llu "
          "rows=%llu/%llu/%llu/%llu "
          "verify_pages=%llu/%llu prefixes=%llu "
          "tickets=%llu/%llu fallback=%llu ready_rows=%llu "
          "ranges=%llu bytes=%llu/%llu candidate_math=%llu/%llu "
          "verify_score=%llu/%llu/%llu/%llu/%llu/%llu "
          "work=%llu work_wall_sum=%.6f\n",
          pool->generation,
          solver_profile.page_plan_descriptors_total,
          solver_profile.page_plan_descriptors_complete,
          solver_profile.page_plan_descriptor_splits,
          solver_profile.page_plan_raw_ranges,
          solver_profile.page_plan_unique_pages,
          solver_profile.page_plan_ranges_after_dedup,
          solver_profile.page_plan_logical_bytes,
          solver_profile.page_plan_aligned_bytes,
          solver_profile.page_plan_overread_bytes,
          solver_profile.page_plan_not_applicable,
          solver_profile.page_plan_allocation_refused,
          solver_profile.page_plan_source_mismatch,
          solver_profile.page_plan_invalid_range,
          solver_profile.page_plan_byte_budget_refused,
          solver_profile.page_plan_range_capacity_refused,
          solver_profile.page_plan_service_refused,
          solver_profile.page_plan_service_errors,
          solver_profile.page_plan_cancelled,
          solver_profile.candidate_delivery_candidates,
          solver_profile.candidate_quad_submitted,
          solver_profile.candidate_quad_ready,
          solver_profile.candidate_quad_fallback,
          solver_profile.candidate_star_submitted,
          solver_profile.candidate_star_ready,
          solver_profile.candidate_star_fallback,
          solver_profile.candidate_delivery_windows,
          solver_profile.candidate_quad_ready_rows,
          solver_profile.candidate_star_ready_rows,
          solver_profile.candidate_retired_rows,
          solver_profile.candidate_native_rows,
          solver_profile.verification_page_queries,
          solver_profile.verification_page_queries_planned,
          solver_profile.verification_page_prefixes,
          solver_profile.verification_page_submitted,
          solver_profile.verification_page_ready,
          solver_profile.verification_page_fallback,
          solver_profile.verification_page_ready_rows,
          solver_profile.verification_page_ranges,
          solver_profile.verification_page_logical_bytes,
          solver_profile.verification_page_aligned_bytes,
          solver_profile.candidate_math_prepared,
          solver_profile.candidate_math_reused,
          solver_profile.verification_score_batches_prepared,
          solver_profile.verification_score_contexts_prepared,
          solver_profile.verification_score_batches_executed,
          solver_profile.verification_score_contexts_completed,
          solver_profile.verification_score_fallback_batches,
          solver_profile.verification_score_stopped_batches,
          solver_profile.verification_score_work_units_completed,
          solver_profile.verification_score_wall_seconds);

  logverb("[index-shard] pass-detail candidates=%zu reduced=%zu "
          "hit_total_cpu_limit=%i hit_total_wall_limit=%i "
          "cancelled=%i rc=%i status=%i "
          "master_committed=%i winner_selected=%i terminal=%s "
          "selected_order=%zu selected_sequence=%zu "
          "task_local_failures=%llu global_integrity_failures=%llu "
          "late_loser_failures=%llu\n",
          nindexes,
          pass_metrics.reduced,
          bp->hit_total_cpulimit,
          bp->hit_total_timelimit,
          bp->cancelled,
          rc,
          (int)status,
          state.master_committed,
          state.winner_selected,
          index_shard_terminal_cause_name(state.terminal_cause),
          state.selected_index_order,
          state.selected_candidate_sequence,
          state.task_local_failures,
          state.global_integrity_failures,
          pool->shared.late_loser_failures);

  logverb("[index-shard] mmap-policy "
         "policy=%s effective=%s scope=all-chunks pass=%u "
         "clean_unsolved_passes=%u transitions=%u "
         "transitioned=%i completed=%i exhaustive=%i "
         "solved=%i cancelled=%i advice_failures=%llu\n",
         fitsbin_mmap_policy_name(
             base_sp->index_mmap_policy.policy),
         fitsbin_mmap_advice_name(
             base_sp->index_mmap_policy.effective_advice),
         base_sp->index_mmap_policy.pass_number,
         base_sp->index_mmap_policy
             .completed_clean_unsolved_passes,
         base_sp->index_mmap_policy.transition_count,
         mmap_transitioned ? 1 : 0,
         pass_completed ? 1 : 0,
         pass_exhaustive ? 1 : 0,
         pass_solved ? 1 : 0,
         pass_cancelled ? 1 : 0,
         pool->shared.mmap_advice_failures);

   if (!task_profile.executed) {
    logverb("[index-shard] task-profile executed=0\n");
  } else if (task_profile.quantiles_available) {
    logverb("[index-shard] task-profile executed=%zu "
           "task_p50=%.6f task_p90=%.6f task_p99=%.6f "
           "task_max=%.6f max_order=%zu max_worker=%i "
           "max_solve=%.6f skew_max_p50=%.1f "
           "max_pool_pct=%.1f%% serial_tail=%.6f "
           "serial_tail_pct=%.1f%% tail_order=%zu tail_worker=%i\n",
           task_profile.executed,
           task_profile.task_p50_seconds,
           task_profile.task_p90_seconds,
           task_profile.task_p99_seconds,
           task_profile.task_max_seconds,
           task_profile.max_index_order,
           task_profile.max_worker_id,
           task_profile.max_solve_seconds,
           task_profile.max_to_p50,
           task_profile.max_pool_percent,
           task_profile.serial_tail_seconds,
           task_profile.serial_tail_percent,
           task_profile.tail_index_order,
           task_profile.tail_worker_id);
  } else {
    logverb("[index-shard] task-profile executed=%zu "
           "quantiles=unavailable task_max=%.6f "
           "max_order=%zu max_worker=%i max_solve=%.6f "
           "max_pool_pct=%.1f%% serial_tail=%.6f "
           "serial_tail_pct=%.1f%% tail_order=%zu tail_worker=%i\n",
           task_profile.executed,
           task_profile.task_max_seconds,
           task_profile.max_index_order,
           task_profile.max_worker_id,
           task_profile.max_solve_seconds,
           task_profile.max_pool_percent,
           task_profile.serial_tail_seconds,
           task_profile.serial_tail_percent,
           task_profile.tail_index_order,
           task_profile.tail_worker_id);
  }

  if (!phase_profile.executed) {
    logverb("[index-shard] phase-profile executed=0\n");
  } else if (phase_profile.quantiles_available) {
    logverb("[index-shard] phase-profile executed=%zu "
           "task_work_wall_sum=%.6f "
           "reset_work_wall_sum=%.6f reset_percent=%.1f "
           "acquire_work_wall_sum=%.6f acquire_percent=%.1f "
           "solve_work_wall_sum=%.6f solve_percent=%.1f "
           "analyze_work_wall_sum=%.6f analyze_percent=%.1f "
           "release_work_wall_sum=%.6f release_percent=%.1f "
           "other_work_wall_sum=%.6f other_percent=%.1f "
           "acquire_p50=%.6f acquire_p90=%.6f "
           "acquire_p99=%.6f acquire_max=%.6f "
           "solve_p50=%.6f solve_p90=%.6f "
           "solve_p99=%.6f solve_max=%.6f\n",
           phase_profile.executed,
           phase_profile.task_wall_total,
           phase_profile.reset_total,
           phase_profile.reset_percent,
           phase_profile.acquire_total,
           phase_profile.acquire_percent,
           phase_profile.solve_total,
           phase_profile.solve_percent,
           phase_profile.analyze_total,
           phase_profile.analyze_percent,
           phase_profile.release_total,
           phase_profile.release_percent,
           phase_profile.other_total,
           phase_profile.other_percent,
           phase_profile.acquire_p50,
           phase_profile.acquire_p90,
           phase_profile.acquire_p99,
           phase_profile.acquire_max,
           phase_profile.solve_p50,
           phase_profile.solve_p90,
           phase_profile.solve_p99,
           phase_profile.solve_max);
  } else {
    logverb("[index-shard] phase-profile executed=%zu "
           "task_work_wall_sum=%.6f "
           "reset_work_wall_sum=%.6f reset_percent=%.1f "
           "acquire_work_wall_sum=%.6f acquire_percent=%.1f "
           "solve_work_wall_sum=%.6f solve_percent=%.1f "
           "analyze_work_wall_sum=%.6f analyze_percent=%.1f "
           "release_work_wall_sum=%.6f release_percent=%.1f "
           "other_work_wall_sum=%.6f other_percent=%.1f "
           "quantiles=unavailable\n",
           phase_profile.executed,
           phase_profile.task_wall_total,
           phase_profile.reset_total,
           phase_profile.reset_percent,
           phase_profile.acquire_total,
           phase_profile.acquire_percent,
           phase_profile.solve_total,
           phase_profile.solve_percent,
           phase_profile.analyze_total,
           phase_profile.analyze_percent,
           phase_profile.release_total,
           phase_profile.release_percent,
           phase_profile.other_total,
           phase_profile.other_percent);
  }

  if (pass_metrics.resource_available) {
    logverb("[index-shard] pass-resource user=%.6f sys=%.6f "
           "minflt=%ld majflt=%ld nvcsw=%ld nivcsw=%ld "
           "inblock=%ld oublock=%ld\n",
           pass_metrics.user_seconds,
           pass_metrics.system_seconds,
           pass_metrics.minor_faults,
           pass_metrics.major_faults,
           pass_metrics.voluntary_context_switches,
           pass_metrics.involuntary_context_switches,
           pass_metrics.filesystem_input_blocks,
           pass_metrics.filesystem_output_blocks);
  } else {
    logverb("[index-shard] pass-resource unavailable\n");
  }

  pthread_mutex_lock(&pool->shared.queue_mutex);
  pool->shared.outer_states = NULL;
  pool->shared.producer_width = 0U;
  pool->shared.helper_width = 0U;
  pthread_mutex_unlock(&pool->shared.queue_mutex);
  pool->shared.worker_view = NULL;
  hooks->destroy_worker_view(worker_view);
  worker_view = NULL;
  free(outer_states);
  index_shard_pool_release_pass(pool);
  return status;
}

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

#ifdef TESTING_INDEX_SHARD_STAGED

typedef struct index_shard_staged_retire_test_context {
  index_shard_thread_state_t *shared;
  size_t order[8];
  size_t calls;
  size_t task_zero_calls;
  anbool stop_during_retire;
} index_shard_staged_retire_test_context_t;

static index_shard_staged_retire_status_t
index_shard_staged_retire_test_callback(
    const index_shard_staged_task_t *task,
    size_t task_index,
    void *opaque) {
  index_shard_staged_retire_test_context_t *context = opaque;

  if (!task || !context ||
      context->calls >= sizeof(context->order) /
          sizeof(context->order[0])) {
    return INDEX_SHARD_STAGED_RETIRE_ERROR;
  }
  context->order[context->calls++] = task_index;
  if (context->stop_during_retire) {
    pthread_mutex_lock(&context->shared->state_mutex);
    context->shared->stop_requested = TRUE;
    pthread_mutex_unlock(&context->shared->state_mutex);
    return INDEX_SHARD_STAGED_RETIRE_MORE;
  }
  if (task_index == 0U && context->task_zero_calls++ < 2U) {
    return INDEX_SHARD_STAGED_RETIRE_MORE;
  }
  return INDEX_SHARD_STAGED_RETIRE_OK;
}

static int index_shard_staged_retire_test_init(
    index_shard_thread_state_t *shared) {
  int rc;

  memset(shared, 0, sizeof(*shared));
  rc = pthread_mutex_init(&shared->queue_mutex, NULL);
  if (rc) {
    return -1;
  }
  rc = pthread_mutex_init(&shared->state_mutex, NULL);
  if (rc) {
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  rc = pthread_cond_init(&shared->queue_cv, NULL);
  if (rc) {
    pthread_mutex_destroy(&shared->state_mutex);
    pthread_mutex_destroy(&shared->queue_mutex);
    return -1;
  }
  return 0;
}

static void index_shard_staged_retire_test_destroy(
    index_shard_thread_state_t *shared) {
  pthread_cond_destroy(&shared->queue_cv);
  pthread_mutex_destroy(&shared->state_mutex);
  pthread_mutex_destroy(&shared->queue_mutex);
}

static int index_shard_staged_retire_test_ready(
    index_shard_thread_state_t *shared,
    index_shard_staged_group_t *group,
    index_shard_staged_task_t *task) {
  int rc;

  pthread_mutex_lock(&shared->queue_mutex);
  rc = index_shard_staged_set_state_locked(
      shared, group, task, INDEX_SHARD_STAGED_TASK_RESULTS_READY);
  pthread_mutex_unlock(&shared->queue_mutex);
  return rc;
}

static int index_shard_staged_retire_test_order(void) {
  index_shard_thread_state_t shared;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_staged_retire_test_context_t context;
  int failures = 0;

  if (index_shard_staged_retire_test_init(&shared)) {
    return 1;
  }
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&context, 0, sizeof(context));
  group.tasks = tasks;
  group.task_count = 2U;
  group.retire = index_shard_staged_retire_test_callback;
  group.owner_context = &context;
  context.shared = &shared;
  tasks[0].scheduler_state =
      INDEX_SHARD_STAGED_TASK_PREPARE_READY;
  tasks[1].scheduler_state =
      INDEX_SHARD_STAGED_TASK_RESULTS_READY;
  group.reorder_ready = 1U;
  shared.staged_reorder_ready = 1U;

  failures += index_shard_staged_retire_one(&shared, &group) != 1;
  failures += context.calls != 0U;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_PREPARE_READY;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_PREPARE_READY;
  failures += group.next_retire != 0U;

  failures += index_shard_staged_retire_test_ready(
      &shared, &group, &tasks[0]) != 0;
  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RETIRED;
  failures += group.next_retire != 1U;

  failures += index_shard_staged_retire_one(&shared, &group) != 0;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_RETIRED;
  failures += group.next_retire != 2U;
  failures += context.calls != 4U;
  failures += context.order[0] != 0U;
  failures += context.order[1] != 0U;
  failures += context.order[2] != 0U;
  failures += context.order[3] != 1U;
  failures += group.reorder_ready != 0U;
  failures += shared.staged_reorder_ready != 0U;
  failures += group.internal_error;
  failures += group.task_failed;

  index_shard_staged_retire_test_destroy(&shared);
  return failures;
}

static int index_shard_staged_retire_test_terminal(void) {
  index_shard_thread_state_t shared;
  index_shard_staged_group_t group;
  index_shard_staged_task_t tasks[2];
  index_shard_staged_retire_test_context_t context;
  int failures = 0;

  if (index_shard_staged_retire_test_init(&shared)) {
    return 1;
  }
  memset(&group, 0, sizeof(group));
  memset(tasks, 0, sizeof(tasks));
  memset(&context, 0, sizeof(context));
  group.tasks = tasks;
  group.task_count = 2U;
  group.retire = index_shard_staged_retire_test_callback;
  group.owner_context = &context;
  context.shared = &shared;
  context.stop_during_retire = TRUE;
  tasks[0].scheduler_state =
      INDEX_SHARD_STAGED_TASK_RESULTS_READY;
  tasks[1].scheduler_state =
      INDEX_SHARD_STAGED_TASK_COMPUTE_READY;
  group.reorder_ready = 1U;
  group.compute_ready = 1U;
  shared.staged_reorder_ready = 1U;
  shared.staged_compute_ready = 1U;

  failures += index_shard_staged_retire_one(&shared, &group) != -1;
  failures += context.calls != 1U;
  failures += context.order[0] != 0U;
  failures += group.next_retire != 0U;
  failures += tasks[0].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_STOPPED;
  failures += tasks[1].scheduler_state !=
      INDEX_SHARD_STAGED_TASK_STOPPED;
  failures += !group.cancelling;
  failures += !group.stop_seen;
  failures += group.reorder_ready != 0U;
  failures += group.compute_ready != 0U;
  failures += shared.staged_reorder_ready != 0U;
  failures += shared.staged_compute_ready != 0U;
  failures += group.internal_error;
  failures += group.task_failed;

  index_shard_staged_retire_test_destroy(&shared);
  return failures;
}

static int index_shard_staged_child_helper_test(void) {
  index_shard_worker_context_t context;
  index_shard_staged_group_t group;
  index_shard_helper_group_t helper;
  int failures = 0;

  memset(&context, 0, sizeof(context));
  memset(&group, 0, sizeof(group));
  memset(&helper, 0, sizeof(helper));

  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.published_staged_group = &group;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.staged_owner_callback_active = TRUE;
  context.staged_owner_callback_group = &group;
  failures += index_shard_helper_staged_child_allowed(
      &context) != TRUE;
  context.published_helper_group = &helper;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  context.published_helper_group = NULL;
  context.staged_owner_callback_group = NULL;
  failures += index_shard_helper_staged_child_allowed(
      &context) != FALSE;
  return failures;
}

int index_shard_test_staged_retire_more(void) {
  return index_shard_staged_retire_test_order() +
      index_shard_staged_retire_test_terminal() +
      index_shard_staged_child_helper_test();
}

#endif
