/*
 * SECTION INDEX-SHARD: module-overview
 *
 * pthread index-sharding for onefield_run()
 *
 * This module executes one candidate index as one outer shard task.  An outer
 * owner may also publish bounded, index-free helper tasks.  It does not split
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
// ANCHOR INDEX-SHARD: result-state
/*
 * Result produced by exactly one shard task.
 *
 * The worker owns this object until completed[index_order] is published.
 * The reducer may then transfer MatchObj payloads into master bp->solutions.
 *
 * Important:
 *   - solutions is worker-local until merge
 *   - merged prevents double-free / double-merge
 *   - solved means this shard produced an accepted candidate or local solve flag
 */
typedef struct index_shard_result {
  bl *solutions; // worker-local MatchObj list for this index

  /* Immutable after task publication; aggregated only after pass quiescence. */
  solver_profile_t solver_profile;

  int failed; // hard task failure, not normal "did not solve"
  int rc;

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

  anbool cancelled;

  size_t index_order; // original candidate index order in onefield pass
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
 *   - queue_mutex protects claim state, owner count, and solved frontier
 *   - result_mutex protects completed slots + active worker count
 *   - state_mutex protects stop/fatal/committed-solve pass state
 *   - limit_mutex protects process-wide CPU-limit publication
 *
 * Do not store per-worker heavy data here.  Per-worker context belongs in
 * index_shard_worker_context_t.
 */
typedef struct index_shard_affinity_entry {
  char *identity;
  size_t canonical_order;
  int owner_worker;
} index_shard_affinity_entry_t;

typedef struct index_shard_affinity_domain {
  size_t scale_index;
  index_shard_affinity_entry_t *entries;
  size_t count;
  size_t capacity;
  struct index_shard_affinity_domain *next;
} index_shard_affinity_domain_t;

typedef struct index_shard_pool index_shard_pool_t;
typedef struct index_shard_helper_group index_shard_helper_group_t;

typedef struct index_shard_thread_state {
  onefield_t *bp;                   // master bp, reducer-owned for writes
  const solver_t *base_sp;          // read-only template for local solvers
  const index_shard_hooks_t *hooks; // bridge back into onefield.c

  size_t nindexes;
  size_t canonical_scan_cursor;
  size_t outer_unclaimed;
  size_t outer_running;
  size_t queue_waiters;
  size_t helper_groups_active;
  size_t helper_preparations_active;
  size_t helper_foreign_reservations;
  unsigned char *outer_states;
  int *preferred_owner;
  size_t *preferred_next;
  size_t *preferred_head;

  index_shard_result_t *results;
  unsigned char *completed; // result slot is visible to reducer
  size_t next_reduce;       // ordered prefix reducer cursor

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
  int solved_published; // reducer committed a valid solved result
  int master_committed;
  double first_stop_wall_since_pass;

  /*
   * Hot solvers read this with an atomic load through worker TLS. The locked
   * state above remains authoritative; this flag only shortens their unwind.
   */
  int worker_stop_requested;

  /* Reducer-owned identity of the first and only master solution commit. */
  int have_committed_result;
  size_t committed_index_order;

  int have_solved_order;
  size_t solved_index_order;

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
  unsigned long long affinity_claims;
  unsigned long long fallback_claims;
  unsigned long long affinity_reassignments;
  unsigned long long helper_groups_published;
  unsigned long long helper_groups_completed;
  unsigned long long helper_tasks_owner;
  unsigned long long helper_tasks_foreign;
  unsigned long long helper_task_failures;
  unsigned long long helper_owner_wait_calls;
  double helper_owner_wait_seconds;
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
  INDEX_SHARD_HELPER_TASK_DONE = 3
} index_shard_helper_task_state_t;

/*
 * One bounded synchronous group published from an outer-index owner stack.
 * queue_mutex protects the lane pointer, task states, and every counter.
 */
struct index_shard_helper_group {
  const index_shard_helper_ops_t *ops;
  index_shard_helper_task_t *tasks;
  size_t task_count;

  unsigned long generation;
  unsigned long long owner_epoch;
  int owner_worker;
  size_t owner_index_order;

  size_t next_claim;
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
  index_shard_affinity_domain_t *affinity_domains;

  int shutdown;
  int stopping;
  int pass_active;
  int ready_workers;
  int tls_startup_error;
  unsigned long generation; // pass submission counter

  index_shard_thread_state_t shared;

} ;

static index_shard_pool_t *index_shard_global_pool = NULL;

static index_shard_affinity_domain_t*
index_shard_affinity_get_domain(
    index_shard_pool_t* pool,
    size_t scale_index,
    anbool create) {
  index_shard_affinity_domain_t* domain;

  if (!pool) {
    return NULL;
  }
  for (domain = pool->affinity_domains;
       domain;
       domain = domain->next) {
    if (domain->scale_index == scale_index) {
      return domain;
    }
  }
  if (!create) {
    return NULL;
  }
  domain = calloc(1, sizeof(*domain));
  if (!domain) {
    return NULL;
  }
  domain->scale_index = scale_index;
  domain->next = pool->affinity_domains;
  pool->affinity_domains = domain;
  return domain;
}

static size_t index_shard_affinity_get_entry(
    index_shard_affinity_domain_t* domain,
    const char* identity,
    size_t canonical_order) {
  index_shard_affinity_entry_t* grown;
  char* identity_copy;
  size_t capacity;
  size_t i;

  if (!domain || !identity || !identity[0]) {
    return SIZE_MAX;
  }
  for (i = 0U; i < domain->count; i++) {
    if (domain->entries[i].canonical_order ==
            canonical_order &&
        !strcmp(domain->entries[i].identity, identity)) {
      return i;
    }
  }
  identity_copy = strdup(identity);
  if (!identity_copy) {
    return SIZE_MAX;
  }
  if (domain->count == domain->capacity) {
    capacity = domain->capacity ? 2U * domain->capacity : 32U;
    if (capacity < domain->capacity ||
        capacity > SIZE_MAX / sizeof(*grown)) {
      free(identity_copy);
      return SIZE_MAX;
    }
    grown = realloc(
        domain->entries,
        capacity * sizeof(*grown));
    if (!grown) {
      free(identity_copy);
      return SIZE_MAX;
    }
    domain->entries = grown;
    domain->capacity = capacity;
  }
  i = domain->count++;
  domain->entries[i].identity = identity_copy;
  domain->entries[i].canonical_order = canonical_order;
  domain->entries[i].owner_worker = -1;
  return i;
}

static index_shard_affinity_domain_t*
index_shard_affinity_prepare_pass(
    index_shard_pool_t* pool,
    onefield_t* bp,
    const index_shard_hooks_t* hooks,
    size_t nindexes,
    size_t* affinity_slots,
    int* preferred_owner,
    size_t* preferred_next,
    size_t* preferred_head,
    size_t* preferred_tail) {
  index_shard_affinity_domain_t* domain = NULL;
  size_t i;

  if (!pool || !affinity_slots || !preferred_owner ||
      !preferred_next || !preferred_head || !preferred_tail) {
    return NULL;
  }
  for (i = 0U; i < nindexes; i++) {
    affinity_slots[i] = SIZE_MAX;
    preferred_owner[i] = -1;
    preferred_next[i] = SIZE_MAX;
  }
  for (i = 0U; i < (size_t)pool->worker_count; i++) {
    preferred_head[i] = SIZE_MAX;
    preferred_tail[i] = SIZE_MAX;
  }
  if (!bp || !hooks ||
      !hooks->get_index_identity) {
    return NULL;
  }
  domain = index_shard_affinity_get_domain(
      pool,
      bp->engine_scale_index,
      TRUE);
  if (!domain) {
    return NULL;
  }
  for (i = 0U; i < nindexes; i++) {
    const char* identity =
        hooks->get_index_identity(bp, i);
    size_t slot = index_shard_affinity_get_entry(
        domain,
        identity,
        i);
    int owner;

    if (slot == SIZE_MAX) {
      continue;
    }
    affinity_slots[i] = slot;
    owner = domain->entries[slot].owner_worker;
    if (owner < 0 || owner >= pool->worker_count) {
      continue;
    }
    preferred_owner[i] = owner;
    if (preferred_tail[owner] == SIZE_MAX) {
      preferred_head[owner] = i;
    } else {
      preferred_next[preferred_tail[owner]] = i;
    }
    preferred_tail[owner] = i;
  }
  return domain;
}

static void index_shard_affinity_commit_pass(
    index_shard_thread_state_t* shared,
    index_shard_affinity_domain_t* domain,
    const size_t* affinity_slots) {
  size_t i;

  if (!shared || !domain || !affinity_slots) {
    return;
  }
  for (i = 0U; i < shared->nindexes; i++) {
    const index_shard_result_t* result =
        &shared->results[i];
    size_t slot = affinity_slots[i];
    int previous;

    if (!result->task_started ||
        !shared->completed[i] ||
        !shared->outer_states ||
        shared->outer_states[i] !=
            INDEX_SHARD_OUTER_FINISHED ||
        result->failed || result->rc ||
        result->cancelled ||
        result->worker_id < 0 ||
        result->worker_id >= shared->worker_count ||
        slot == SIZE_MAX || slot >= domain->count) {
      continue;
    }
    previous = domain->entries[slot].owner_worker;
    if (previous >= 0 && previous != result->worker_id) {
      shared->affinity_reassignments++;
    }
    domain->entries[slot].owner_worker = result->worker_id;
  }
}

static void index_shard_affinity_destroy(
    index_shard_pool_t* pool) {
  index_shard_affinity_domain_t* domain;

  if (!pool) {
    return;
  }
  domain = pool->affinity_domains;
  while (domain) {
    index_shard_affinity_domain_t* next = domain->next;
    size_t i;

    for (i = 0U; i < domain->count; i++) {
      free(domain->entries[i].identity);
    }
    free(domain->entries);
    free(domain);
    domain = next;
  }
  pool->affinity_domains = NULL;
}

static pthread_mutex_t index_shard_global_pool_mutex = PTHREAD_MUTEX_INITIALIZER;

static pthread_key_t index_shard_tls_key;
static pthread_once_t index_shard_tls_once = PTHREAD_ONCE_INIT;
static int index_shard_tls_key_status = EAGAIN;

static size_t index_shard_inverse_cache_budget(void) {
  const size_t ceiling = 128U * 1024U * 1024U;
  long available_pages;
  long page_size;
  size_t available_bytes;
  size_t budget;

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
  budget = available_bytes / 64U;
  return MIN(budget, ceiling);
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
  int solved_published;
  int master_committed;
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
  snapshot->solved_published = shared->solved_published;
  snapshot->master_committed = shared->master_committed;
  snapshot->first_stop_wall_since_pass =
      shared->first_stop_wall_since_pass;

  pthread_mutex_unlock(&shared->state_mutex);
}

/*
 * Snapshot completed-pass metrics.
 *
 * This function is called only after index_shard_pool_reduce_online() has
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

  snapshot->reduced = shared->next_reduce;
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
 * This runs only after index_shard_pool_reduce_online() has waited for every
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

// ANCHOR INDEX-SHARD: request-stop
/*
 * Cooperative global stop.
 *
 * Stop means:
 *   - workers should not claim new tasks
 *   - reducer should wake and merge any solved/completed result
 *   - already-running solver calls must exit through callback polling
 *
 * This is intentionally not pthread_cancel.
 */
static void index_shard_request_stop(index_shard_thread_state_t *shared) {
  int was_stopped;

  pthread_mutex_lock(&shared->state_mutex);
  was_stopped = shared->stop_requested;
  if (!was_stopped) {
    shared->first_stop_wall_since_pass =
        monotonic_seconds() - shared->pass_wall_start;
  }
  shared->stop_requested = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);

  __atomic_store_n(
      &shared->worker_stop_requested,
      TRUE,
      __ATOMIC_RELEASE);

  if (!was_stopped && index_shard_trace_enabled()) {
    logmsg("[index-shard] stop-request pass_wall=%.6f\n",
           monotonic_seconds() - shared->pass_wall_start);
  }

  index_shard_wake_pass_waiters(shared);
  index_shard_wake_queue_waiters(shared);
}

static void index_shard_publish_solved(
    index_shard_thread_state_t *shared,
    size_t index_order);

/*
 * Publish a speculative solved frontier without cancelling earlier indexes.
 *
 * Affinity claims may be out of order. A solved index narrows the claim
 * ceiling, but every earlier unclaimed index remains eligible. The ordered
 * reducer is still the only path that can publish the canonical winner and
 * stop the complete pool.
 */
void index_shard_worker_publish_solution_candidate(void) {
  index_shard_worker_context_t *ctx = index_shard_get_tls();

  if (!ctx || !ctx->pool ||
      !ctx->current_outer_active) {
    return;
  }

  index_shard_publish_solved(
      &ctx->pool->shared,
      ctx->current_index_order);
}

// ANCHOR INDEX-SHARD: request-fatal-stop
/*
 * Publish a hard worker/module failure and stop the current pass.
 */
static void index_shard_request_fatal_stop(index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->state_mutex);
  if (!shared->stop_requested) {
    shared->first_stop_wall_since_pass =
        monotonic_seconds() - shared->pass_wall_start;
  }
  shared->fatal_error = TRUE;
  shared->stop_requested = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);

  __atomic_store_n(
      &shared->worker_stop_requested,
      TRUE,
      __ATOMIC_RELEASE);

  index_shard_wake_pass_waiters(shared);
  index_shard_wake_queue_waiters(shared);
}

// ANCHOR INDEX-SHARD: publish-committed-solve
/*
 * Publish that the reducer committed a valid solved result.
 *
 * Quick commit is intentional project policy. Once a solved result becomes
 * master-visible, no new shard tasks should be claimed.
 */
static void index_shard_publish_committed_solve(
    index_shard_thread_state_t *shared) {
  pthread_mutex_lock(&shared->state_mutex);
  if (!shared->stop_requested) {
    shared->first_stop_wall_since_pass =
        monotonic_seconds() - shared->pass_wall_start;
  }
  shared->solved_published = TRUE;
  shared->stop_requested = TRUE;
  pthread_mutex_unlock(&shared->state_mutex);

  __atomic_store_n(
      &shared->worker_stop_requested,
      TRUE,
      __ATOMIC_RELEASE);

  index_shard_wake_pass_waiters(shared);
  index_shard_wake_queue_waiters(shared);
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
static int index_shard_master_limit_or_cancel_requested(index_shard_thread_state_t *shared,
                                                        anbool *hit_total_cpulimit,
                                                        anbool *cancelled) {
  onefield_t *bp = shared->bp;
  int stop;

  pthread_mutex_lock(&shared->limit_mutex);

  if (hit_total_cpulimit) {
    *hit_total_cpulimit = bp->hit_total_cpulimit;
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
 * Solver completion is mirrored through solved_published.  Workers therefore
 * do not read master bp->single_field_solved concurrently with the reducer.
 */
static int index_shard_master_stop_requested(index_shard_thread_state_t *shared) {
  index_shard_pass_state_snapshot_t state;

  index_shard_pass_state_snapshot(shared, &state);

  if (state.stop_requested || state.fatal_error || state.solved_published) {
    return TRUE;
  }

  return index_shard_master_limit_or_cancel_requested(shared, NULL, NULL);
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
  index_shard_pass_state_snapshot_t state;
  int hit = FALSE;

  index_shard_pass_state_snapshot(shared, &state);

  if (state.stop_requested || state.fatal_error || state.solved_published) {
    return TRUE;
  }

  pthread_mutex_lock(&shared->limit_mutex);

  if (bp->cancelled || bp->hit_total_cpulimit || bp->hit_total_timelimit) {
    hit = TRUE;
  } else {
    if (bp->total_timelimit > 0.0) {
      double now = monotonic_seconds();

      if (now >= 0.0 &&
          now - bp->time_total_start >= bp->total_timelimit) {
        bp->hit_total_timelimit = TRUE;
        hit = TRUE;

        if (!shared->limit_reported) {
          shared->limit_reported = TRUE;
          logmsg("Total wall-clock time limit reached!\n");
          logverb("[index-shard] wall-limit reached total_timelimit=%g "
                  "elapsed=%.6f\n",
                  bp->total_timelimit,
                  now - bp->time_total_start);
        }
      }
    }

    if (!hit && bp->total_cpulimit > 0.0) {
      float now = get_cpu_usage();
      double elapsed = (double)(now - bp->cpu_total_start);

      if (elapsed >= bp->total_cpulimit) {
        bp->hit_total_cpulimit = TRUE;
        hit = TRUE;

        if (!shared->limit_reported) {
          shared->limit_reported = TRUE;
          logmsg("Total CPU time limit reached!\n");
          logverb("[index-shard] cpu-budget reached total_cpulimit=%g "
                  "elapsed=%.6f\n",
                  bp->total_cpulimit,
                  elapsed);
        }
      }
    }
  }

  pthread_mutex_unlock(&shared->limit_mutex);

  if (hit) {
    index_shard_request_stop(shared);
  }

  return hit;
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
 * This marks solved/best_logodds before the ordered reducer reaches this
 * slot.  The worker may stop later claims, but only the reducer can commit a
 * solution and cancel earlier-index work.
 */
static void index_shard_capture_solution_analysis(index_shard_thread_state_t *shared,
                                                  index_shard_result_t *result) {
  if (!shared->hooks || !shared->hooks->analyze_solutions)
    return;

  result->solved = shared->hooks->analyze_solutions(shared->bp, result->solutions,
                                                    &result->best_logodds, &result->best_fieldnum);
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

  if (!result || result->merged) {
    return 0;
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
    result->failed = TRUE;
    result->rc = -1;
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
    result->failed = TRUE;
    result->rc = -1;
    index_shard_reducer_work_account(shared, reducer_wall_start);
    return -1;
  }

  result->merged = TRUE;

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
static index_t *index_shard_worker_get_index(index_shard_worker_context_t *ctx,
                                             index_shard_thread_state_t *shared,
                                             size_t index_order) {
  // original onefield ownership path, released after task
  index_t *index;

  if (!shared->hooks || !shared->hooks->get_index)
    return NULL;

  index = shared->hooks->get_index(shared->bp, index_order);

  if (index_shard_trace_enabled() && index) {
    logmsg("[index-shard] worker=%i load index_order=%zu index=%s\n", ctx->worker_id, index_order,
           index->indexname ? index->indexname : "(null)");
  }

  return index;
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
     * Apply the pass policy to payload chunks of retained or newly prepared
     * index components. fitsbin keeps topology chunks under NORMAL.
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

/* queue_mutex must be held. */
static anbool index_shard_helper_outer_claimable_locked(
    const index_shard_thread_state_t *shared) {
  size_t limit;
  size_t candidate;

  if (!shared || !shared->outer_states) {
    return FALSE;
  }
  limit = shared->have_solved_order ?
      shared->solved_index_order : shared->nindexes;
  candidate = shared->canonical_scan_cursor;
  while (candidate < limit &&
         shared->outer_states[candidate] !=
             INDEX_SHARD_OUTER_UNCLAIMED) {
    candidate++;
  }
  return candidate < limit;
}

/* queue_mutex must be held. */
static size_t index_shard_helper_idle_workers_locked(
    const index_shard_thread_state_t *shared) {
  size_t available;
  size_t limit;
  size_t spare = 0U;

  if (!shared || shared->worker_count < 2 ||
      index_shard_helper_outer_claimable_locked(shared)) {
    return 0U;
  }
  limit = (size_t)shared->worker_count - 1U;
  if ((size_t)shared->worker_count >
      shared->outer_running) {
    spare = (size_t)shared->worker_count -
        shared->outer_running;
  }
  available = shared->queue_waiters;
  if (!shared->helper_groups_active &&
      spare > available) {
    available = spare;
  }
  if (shared->helper_foreign_reservations >= available) {
    return 0U;
  }
  available -= shared->helper_foreign_reservations;
  return available < limit ? available : limit;
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
    if (shared->have_solved_order &&
        candidate->owner_index_order >=
            shared->solved_index_order) {
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

static int index_shard_helper_complete_claim(
    index_shard_thread_state_t *shared,
    const index_shard_helper_claim_t *claim,
    index_shard_helper_task_status_t execute_status) {
  index_shard_helper_group_t *group;
  index_shard_helper_task_t *task;
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
      group->task_failed || group->stop_seen) {
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
  if (stopped ||
      (shared->have_solved_order &&
       group->owner_index_order >=
           shared->solved_index_order)) {
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
  if (!fatal && !stopped &&
      (!shared->have_solved_order ||
       group->owner_index_order <
           shared->solved_index_order)) {
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
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active &&
      (!shared->have_solved_order ||
       ctx->current_index_order < shared->solved_index_order)) {
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
      !ctx->helper_preparation_active &&
      !shared->helper_preparations_active &&
      (!shared->have_solved_order ||
       ctx->current_index_order < shared->solved_index_order)) {
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

index_shard_helper_run_status_t
index_shard_helper_run(
    const index_shard_helper_ops_t *ops,
    index_shard_helper_task_t *tasks,
    size_t task_count,
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
  int outer_claimable;
  int preparation_permit = FALSE;

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
  preparation_permit =
      ctx->helper_preparation_active &&
      ctx->helper_preparation_generation ==
          ctx->generation_seen &&
      ctx->helper_preparation_index_order ==
          ctx->current_index_order;
  outer_claimable =
      index_shard_helper_outer_claimable_locked(shared);
  if (ctx->generation_seen != ctx->pool->generation ||
      ctx->published_helper_group ||
      (!preparation_permit &&
       shared->helper_preparations_active) ||
      outer_claimable ||
      (shared->have_solved_order &&
       group.owner_index_order >=
           shared->solved_index_order)) {
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
  pthread_cond_broadcast(&shared->queue_cv);
  pthread_mutex_unlock(&shared->queue_mutex);

  while (1) {
    if (have_claim) {
      (void)index_shard_helper_execute_claim(
          shared, &claim);

      have_claim = FALSE;
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

static size_t index_shard_claim_limit_locked(
    const index_shard_thread_state_t *shared) {
  if (shared->have_solved_order) {
    return shared->solved_index_order;
  }
  return shared->nindexes;
}

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
    anbool affinity_claim,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice) {
  int preferred;

  if (candidate >= shared->nindexes ||
      shared->outer_states[candidate] !=
          INDEX_SHARD_OUTER_UNCLAIMED ||
      !shared->outer_unclaimed ||
      shared->outer_running >=
          (size_t)shared->worker_count) {
    return -1;
  }
  preferred = shared->preferred_owner[candidate];
  shared->outer_states[candidate] =
      INDEX_SHARD_OUTER_RUNNING;
  shared->outer_unclaimed--;
  shared->outer_running++;
  if (affinity_claim) {
    shared->affinity_claims++;
  } else {
    shared->fallback_claims++;
  }
  *index_order = candidate;
  *mmap_advice = shared->mmap_advice;
  if (candidate == shared->canonical_scan_cursor) {
    index_shard_advance_canonical_cursor_locked(shared);
  }

  if (index_shard_trace_enabled()) {
    logmsg("[index-shard] claim index_order=%zu lane=%s "
           "worker=%i preferred=%i owners=%zu "
           "outer_unclaimed=%zu payload=%s "
           "wall_since_pass=%.6f\n",
           candidate,
           affinity_claim ? "affinity" : "fallback",
           worker->worker_id,
           preferred,
           shared->outer_running,
           shared->outer_unclaimed,
           fitsbin_mmap_advice_name(*mmap_advice),
           monotonic_seconds() - shared->pass_wall_start);
  }
  return 0;
}

/*
 * Select one whole-index owner task from the current band.
 *
 * A worker first consumes its canonical-order list of exact indexes retained
 * from the previous depth at this scale. If none remains claimable, it steals
 * the canonical-lowest unstarted index. The preference is never a reservation:
 * a free worker starts useful outer work immediately. The explicit state map
 * permits affinity claims out of order while the reducer remains canonical.
 *
 * The old cross-owner AB executor is deliberately absent here. A worker with
 * no claimable outer index may claim a fixed, index-free package published by
 * a current outer owner. Neither path transfers index ownership or advances
 * the engine band independently.
 */
static index_shard_work_selection_t
index_shard_select_work(
    index_shard_worker_context_t *worker,
    index_shard_thread_state_t *shared,
    size_t *index_order,
    fitsbin_mmap_advice_t *mmap_advice,
    index_shard_helper_claim_t *helper_claim) {
  int worker_id;

  if (!worker || !shared || !index_order || !mmap_advice || !helper_claim ||
      !shared->outer_states || !shared->preferred_owner ||
      !shared->preferred_next || !shared->preferred_head) {
    return INDEX_SHARD_WORK_ERROR;
  }
  worker_id = worker->worker_id;
  if (worker_id < 0 || worker_id >= shared->worker_count) {
    return INDEX_SHARD_WORK_ERROR;
  }

  pthread_mutex_lock(&shared->queue_mutex);
  while (1) {
    index_shard_pass_state_snapshot_t state;
    size_t claim_limit;
    size_t candidate;
    int helper_selection;

    index_shard_pass_state_snapshot(shared, &state);
    if (state.stop_requested || state.fatal_error ||
        state.solved_published) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }

    claim_limit = index_shard_claim_limit_locked(shared);
    while (shared->preferred_head[worker_id] != SIZE_MAX) {
      candidate = shared->preferred_head[worker_id];
      if (candidate >= shared->nindexes) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      shared->preferred_head[worker_id] =
          shared->preferred_next[candidate];
      if (candidate >= claim_limit ||
          shared->outer_states[candidate] !=
              INDEX_SHARD_OUTER_UNCLAIMED) {
        continue;
      }
      if (index_shard_claim_outer_locked(
              worker,
              shared,
              candidate,
              TRUE,
              index_order,
              mmap_advice)) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_OUTER;
    }

    index_shard_advance_canonical_cursor_locked(shared);
    candidate = shared->canonical_scan_cursor;
    if (candidate < claim_limit) {
      if (index_shard_claim_outer_locked(
              worker,
              shared,
              candidate,
              FALSE,
              index_order,
              mmap_advice)) {
        pthread_mutex_unlock(&shared->queue_mutex);
        return INDEX_SHARD_WORK_ERROR;
      }
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_OUTER;
    }

    helper_selection = index_shard_helper_select_locked(
        worker, shared, helper_claim);
    if (helper_selection < 0) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    if (!helper_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_HELPER;
    }

    if (!shared->outer_running) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_DONE;
    }
    shared->queue_waiters++;
    helper_selection = pthread_cond_wait(
        &shared->queue_cv, &shared->queue_mutex);
    if (!shared->queue_waiters) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
    shared->queue_waiters--;
    if (!shared->queue_waiters) {
      pthread_cond_broadcast(&shared->queue_cv);
    }
    if (helper_selection) {
      pthread_mutex_unlock(&shared->queue_mutex);
      return INDEX_SHARD_WORK_ERROR;
    }
  }
}

static int index_shard_mark_result_completed(index_shard_thread_state_t *shared,
                                             size_t index_order) {
  int rc = 0;

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
    shared->completed[index_order] = TRUE;
  }
  pthread_cond_broadcast(&shared->result_cv);
  pthread_mutex_unlock(&shared->result_mutex);
  return rc;
}

static void index_shard_finish_outer_claim(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  int underflow = FALSE;
  int completion_failed;

  /*
   * Publish the terminal owner state before completed[index_order] makes the
   * immutable result visible to the canonical reducer.
   */
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

  completion_failed = index_shard_mark_result_completed(
      shared, index_order);
  if (underflow || completion_failed) {
    index_shard_request_fatal_stop(shared);
  }
}

// ANCHOR INDEX-SHARD: publish-solved
/*
 * Publish solved frontier.
 *
 * This prevents workers from claiming at or above the earliest known solved
 * order. Affinity can leave earlier holes, so those indexes remain claimable.
 * The actual solution merge and global stop still belong to the canonical
 * ordered reducer.
 */
static void index_shard_publish_solved(
    index_shard_thread_state_t *shared,
    size_t index_order) {
  pthread_mutex_lock(&shared->queue_mutex);
  if (index_order < shared->nindexes &&
      (!shared->have_solved_order ||
       index_order < shared->solved_index_order)) {
    shared->have_solved_order = TRUE;
    shared->solved_index_order = index_order;
    pthread_cond_broadcast(&shared->queue_cv);
  }
  pthread_mutex_unlock(&shared->queue_mutex);
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

  if (shared->hooks->prepare_local_context(&ctx->local_bp, shared->bp, shared->base_sp))
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

static int index_shard_done_with_index(
    index_shard_worker_context_t *ctx,
    index_shard_thread_state_t *shared,
    size_t index_order,
    index_t *index,
    index_shard_inverse_lease_t *inverse_lease) {
  if (!index) {
    return 0;
  }
  index_shard_inverse_cache_release(
      ctx, index, inverse_lease);
  if (shared && shared->hooks &&
      shared->hooks->done_with_index) {
    return shared->hooks->done_with_index(
        shared->bp,
        index_order,
        index);
  }
  return 0;
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
  int rc = 0;
  index_shard_inverse_lease_t inverse_lease;

  index_shard_result_init(result, index_order);
  memset(&inverse_lease, 0, sizeof(inverse_lease));
  result->mmap_advice = mmap_advice;

  if (!result->solutions) {
    result->failed = TRUE;
    result->rc = -1;
    return -1;
  }

  if (!ctx->local_context_ready) {
    result->failed = TRUE;
    result->rc = -1;
    return -1;
  }

  if (!shared->hooks || !shared->hooks->reset_local_context_for_task ||
      !shared->hooks->solve_one_index) {
    result->failed = TRUE;
    result->rc = -1;
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

  index = index_shard_worker_get_index(
      ctx,
      shared,
      index_order);

  fitsbin_mmap_clear_thread_advice();

  if (!index) {
    result->acquire_seconds =
        monotonic_seconds() - phase_wall_start;

    ERROR("Failed to load index order %zu", index_order);
    result->failed = TRUE;
    result->rc = -1;

    index_shard_result_finish_task(result,
                                   shared,
                                   task_wall_start);
    return -1;
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
        &result->cancelled);

    phase_wall_start = monotonic_seconds();
    if (index_shard_done_with_index(
            ctx,
            shared,
            index_order,
            index,
            &inverse_lease)) {
      result->failed = TRUE;
      result->rc = -1;
    }
    index = NULL;
    result->release_seconds =
        monotonic_seconds() - phase_wall_start;
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

  // Worker-lifetime TLS lets onefield callbacks publish this exact order.
  ctx->current_index_order = index_order;
  ctx->current_outer_active = TRUE;
  rc = shared->hooks->solve_one_index(&ctx->local_bp, index);
  ctx->current_outer_active = FALSE;

  result->wall_seconds = monotonic_seconds() - wall_start;
  result->cpu_seconds = get_cpu_usage() - cpu_start;

  result->hit_total_cpulimit = ctx->local_bp.hit_total_cpulimit;
  result->cancelled = ctx->local_bp.cancelled;
  result->solver_profile = ctx->local_bp.solver.profile;
  result->rc = rc;

  if (rc) {
    result->failed = TRUE;
  }

  if (result->cancelled) {
    index_shard_request_stop(shared);
  }

  // analyze before reducer so worker can trigger fast stop
  phase_wall_start = monotonic_seconds();

  if (!result->failed) {
    index_shard_capture_solution_analysis(shared, result);

    if (ctx->local_bp.single_field_solved) {
      result->solved = TRUE;
    }
  }

  result->analyze_seconds =
      monotonic_seconds() - phase_wall_start;

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
  if (index_shard_done_with_index(
          ctx,
          shared,
          index_order,
          index,
          &inverse_lease)) {
    result->failed = TRUE;
    result->rc = -1;
    rc = -1;
  }
  index = NULL;

  result->release_seconds =
      monotonic_seconds() - phase_wall_start;

  index_shard_result_finish_task(result,
                                 shared,
                                 task_wall_start);

  if (rc) {
    return rc;
  }

  return 0;
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
      index_shard_helper_claim_t helper_claim;
      fitsbin_mmap_advice_t mmap_advice =
          FITSBIN_MMAP_ADVICE_NORMAL;
      index_shard_work_selection_t selection;

      memset(&helper_claim, 0, sizeof(helper_claim));
      selection = index_shard_select_work(
          ctx,
          shared,
          &index_order,
          &mmap_advice,
          &helper_claim);

      if (selection == INDEX_SHARD_WORK_ERROR) {
        index_shard_request_fatal_stop(shared);
        break;
      }
      if (selection == INDEX_SHARD_WORK_DONE) {
        break;
      }
      if (selection == INDEX_SHARD_WORK_HELPER) {
        if (index_shard_helper_execute_claim(
                shared, &helper_claim)) {
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
          result->failed = TRUE;
          result->rc = -1;

          index_shard_finish_outer_claim(shared, index_order);
          index_shard_request_fatal_stop(shared);
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
        index_shard_request_fatal_stop(shared);
        break;
      }

      if (result->solved) {
        logverb("[index-shard] solved-candidate worker=%i index_order=%zu "
                "best_logodds=%.3f field=%i wall=%.6f cpu=%.6f\n",
                ctx->worker_id,
                index_order,
                result->best_logodds,
                result->best_fieldnum,
                result->wall_seconds,
                (double)result->cpu_seconds);

        index_shard_publish_solved(shared, index_order);
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

      index_shard_check_global_limits(shared);
      index_shard_finish_outer_claim(shared, index_order);
    }

    {
      double context_wall_start = monotonic_seconds();

      index_shard_worker_cleanup_pass(ctx, shared);
      ctx->pass_cleanup_seconds =
          monotonic_seconds() - context_wall_start;
    }

    index_shard_worker_done(shared);
  }

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

// ANCHOR INDEX-SHARD: reduce-online
/*
 * Online reducer for one submitted pass.
 *
 * Reduce only the completed canonical prefix.  A solved result becomes
 * master-visible only after every earlier configured index has completed and
 * been reduced, preserving the original serial index priority independently
 * of worker completion order.
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

static int index_shard_pool_reduce_online(index_shard_pool_t *pool) {
  index_shard_thread_state_t *shared = &pool->shared;
  size_t i;
  int rc = 0;

  while (shared->next_reduce < shared->nindexes) {
    index_shard_pass_state_snapshot_t state;
    int can_reduce = FALSE;
    int workers_done = FALSE;
    int fatal = FALSE;
    int wait_failed = FALSE;

    pthread_mutex_lock(&shared->result_mutex);

    while (!shared->completed[shared->next_reduce] &&
           shared->active_workers > 0) {
      index_shard_pass_state_snapshot(shared, &state);

      if (state.fatal_error) {
        break;
      }

      {
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
          wait_failed = TRUE;
          break;
        }
      }
    }

    can_reduce = shared->completed[shared->next_reduce];
    workers_done = (shared->active_workers == 0);

    index_shard_pass_state_snapshot(shared, &state);
    fatal = state.fatal_error;

    pthread_mutex_unlock(&shared->result_mutex);

    if (wait_failed) {
      rc = -1;
      index_shard_request_fatal_stop(shared);
      break;
    }

    if (fatal) {
      rc = -1;
      index_shard_request_stop(shared);
      break;
    }

    if (can_reduce) {
      index_shard_result_t *result =
          &shared->results[shared->next_reduce];

      if (index_shard_trace_enabled()) {
        logmsg("[index-shard] reduce index_order=%zu solved=%i failed=%i\n",
               shared->next_reduce,
               result->solved,
               result->failed);
      }

      if (result->failed) {
        rc = -1;
        index_shard_request_fatal_stop(shared);
        break;
      }

      if (index_shard_reduce_one_result(shared, result)) {
        rc = -1;
        index_shard_request_fatal_stop(shared);
        break;
      }

      shared->next_reduce++;

      if (index_shard_master_stop_requested(shared)) {
        index_shard_request_stop(shared);
        break;
      }

      continue;
    }

    if (workers_done) {
      break;
    }
  }

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

  pthread_mutex_unlock(&shared->result_mutex);

  /*
   * A worker failure can race with an earlier valid-solution commit.  Never
   * convert that hard failure into a successful pass: after quiescence every
   * result slot is immutable, so perform one authoritative failure scan.
   */
  for (i = 0; i < shared->nindexes; i++) {
    if (shared->completed[i] &&
        (shared->results[i].failed || shared->results[i].rc)) {
      rc = -1;
      index_shard_request_fatal_stop(shared);
      break;
    }
  }

  {
    index_shard_pass_state_snapshot_t state;

    index_shard_pass_state_snapshot(shared, &state);
    if (state.fatal_error) {
      rc = -1;
    }
  }

  while (rc == 0 &&
         shared->next_reduce < shared->nindexes &&
         shared->completed[shared->next_reduce]) {
    index_shard_pass_state_snapshot_t state;
    index_shard_result_t *result =
        &shared->results[shared->next_reduce];

    index_shard_pass_state_snapshot(shared, &state);

    if (state.solved_published) {
      break;
    }

    if (!result->failed && !result->merged) {
      if (index_shard_reduce_one_result(shared, result)) {
        rc = -1;
        break;
      }
    }

    shared->next_reduce++;

    index_shard_pass_state_snapshot(shared, &state);

    if (state.solved_published) {
      break;
    }
  }

  return rc;
}

/*
 * Report the reducer-owned winner after every worker has quiesced and the
 * late hard-failure scan has passed. Worker callbacks deliberately suppress
 * their ordinary solved line, so this is the sole parallel success report.
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
      !shared->have_committed_result ||
      shared->committed_index_order >= nindexes) {
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

  index_shard_global_pool = pool;
  fitsbin_payload_io_configure_workers(worker_count);

  logverb("[index-shard] workers=%i mode=pthread "
          "inverse_cache_budget=%zu\n",
          worker_count,
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
  fitsbin_payload_io_configure_workers(1);

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
  index_shard_affinity_destroy(pool);
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
    index_shard_result_t *results,
    unsigned char *completed,
    unsigned char *outer_states,
    int *preferred_owner,
    size_t *preferred_next,
    size_t *preferred_head) {
  index_shard_thread_state_t *shared = &pool->shared;
  int worker_count = pool->worker_count;
  int i;

  if (!outer_states || !preferred_owner ||
      !preferred_next || !preferred_head) {
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
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] helper state remained active "
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
      shared->helper_groups_completed) {
    logerr("[index-shard] helper group lifecycle remained "
           "before pass published=%llu completed=%llu\n",
           shared->helper_groups_published,
           shared->helper_groups_completed);
    pthread_mutex_unlock(&shared->queue_mutex);
    pthread_mutex_unlock(&pool->control_mutex);
    return -1;
  }
  if (shared->helper_groups_active ||
      shared->helper_preparations_active ||
      shared->helper_foreign_reservations) {
    logerr("[index-shard] helper activity remained before pass "
           "groups=%zu preparations=%zu reservations=%zu\n",
           shared->helper_groups_active,
           shared->helper_preparations_active,
           shared->helper_foreign_reservations);
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
  shared->nindexes = nindexes;
  shared->canonical_scan_cursor = 0U;
  shared->outer_unclaimed = nindexes;
  shared->outer_running = 0U;
  shared->queue_waiters = 0U;
  shared->helper_groups_active = 0U;
  shared->helper_preparations_active = 0U;
  shared->helper_foreign_reservations = 0U;
  shared->outer_states = outer_states;
  shared->preferred_owner = preferred_owner;
  shared->preferred_next = preferred_next;
  shared->preferred_head = preferred_head;
  shared->affinity_claims = 0U;
  shared->fallback_claims = 0U;
  shared->affinity_reassignments = 0U;
  shared->helper_groups_published = 0U;
  shared->helper_groups_completed = 0U;
  shared->helper_tasks_owner = 0U;
  shared->helper_tasks_foreign = 0U;
  shared->helper_task_failures = 0U;
  shared->helper_owner_wait_calls = 0U;
  shared->helper_owner_wait_seconds = 0.0;

  shared->results = results;
  shared->completed = completed;
  shared->next_reduce = 0U;

  shared->worker_count = worker_count;
  shared->active_workers = worker_count;

  shared->stop_requested = FALSE;
  shared->fatal_error = FALSE;
  shared->solved_published = FALSE;
  shared->master_committed = FALSE;
  shared->first_stop_wall_since_pass = -1.0;
  __atomic_store_n(
      &shared->worker_stop_requested,
      FALSE,
      __ATOMIC_RELEASE);
  shared->have_committed_result = FALSE;
  shared->committed_index_order = 0U;
  shared->have_solved_order = FALSE;
  shared->solved_index_order = nindexes;
  shared->limit_reported = FALSE;

  shared->reducer_work_calls = 0U;
  shared->reducer_work_wall_seconds = 0.0;

  /*
   * Payload mappings follow the parallel pass policy while topology remains
   * NORMAL. Bounded exact page population is additive and never changes that
   * base policy.
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

  logverb("[index-shard] pthread-pool submit workers=%i pool_workers=%i "
          "candidates=%zu engine_pass=%zu depth_index=%zu scale_index=%zu "
          "startobj=%i endobj=%i scheduler=affinity-first chunk=1 "
          "inner_scheduler=owner-helper-groups "
          "mmap_pass=%u mmap_advice=%s "
          "mmap_policy=payload-random-topology-normal "
          "page_delivery=bounded-populate "
          "payload_io=mmap-zero-copy credits=%i "
          "outer_admission=full-owner-affinity\n",
          worker_count,
          pool->worker_count,
          nindexes,
          bp->engine_pass_ordinal,
          bp->engine_depth_index,
          bp->engine_scale_index,
          base_sp->startobj,
          base_sp->endobj,
          shared->mmap_pass_number,
          fitsbin_mmap_advice_name(shared->mmap_advice),
          worker_count);

  return 0;
}

/*
 * SECTION INDEX-SHARD: entry
 */

// ANCHOR INDEX-SHARD: entry
/*
 * Execute one complete index-shard pass.
 *
 * Terminal status is classified according to whether master-visible solver
 * state has already been mutated. Only an unavailable path or a failure
 * proven to occur before master commit may return control to the serial path.
 */
static index_shard_solve_status_t
index_shard_solve_impl(onefield_t *bp,
                       solver_t *base_sp,
                       size_t nindexes,
                       const index_shard_hooks_t *hooks) {
  index_shard_pool_t *pool;
  index_shard_affinity_domain_t *affinity_domain = NULL;
  size_t *affinity_slots = NULL;
  int *preferred_owner = NULL;
  size_t *preferred_next = NULL;
  size_t *preferred_head = NULL;
  size_t *preferred_tail = NULL;
  unsigned char *outer_states = NULL;
  index_shard_result_t *results = NULL;
  unsigned char *completed = NULL;
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

  results = calloc(nindexes, sizeof(index_shard_result_t));
  completed = calloc(nindexes, sizeof(unsigned char));
  outer_states = calloc(nindexes, sizeof(*outer_states));
  affinity_slots = calloc(nindexes, sizeof(*affinity_slots));
  preferred_owner = calloc(nindexes, sizeof(*preferred_owner));
  preferred_next = calloc(nindexes, sizeof(*preferred_next));
  preferred_head = calloc(
      (size_t)pool->worker_count,
      sizeof(*preferred_head));
  preferred_tail = calloc(
      (size_t)pool->worker_count,
      sizeof(*preferred_tail));

  if (!results || !completed || !outer_states ||
      !affinity_slots || !preferred_owner || !preferred_next ||
      !preferred_head || !preferred_tail) {
    SYSERROR("Failed to allocate index-shard pass state");

    free(preferred_tail);
    free(preferred_head);
    free(preferred_next);
    free(preferred_owner);
    free(affinity_slots);
    free(outer_states);
    free(results);
    free(completed);

    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  affinity_domain = index_shard_affinity_prepare_pass(
      pool,
      bp,
      hooks,
      nindexes,
      affinity_slots,
      preferred_owner,
      preferred_next,
      preferred_head,
      preferred_tail);

  // Submit releases the hard current-band barrier to persistent workers.
  rc = index_shard_pool_submit(
      pool,
      bp,
      base_sp,
      nindexes,
      hooks,
      results,
      completed,
      outer_states,
      preferred_owner,
      preferred_next,
      preferred_head);

  if (rc) {
    free(preferred_tail);
    free(preferred_head);
    free(preferred_next);
    free(preferred_owner);
    free(affinity_slots);
    free(outer_states);
    free(results);
    free(completed);

    index_shard_pool_release_pass(pool);
    return INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
  }

  rc = index_shard_pool_reduce_online(pool);

  index_shard_pass_state_snapshot(&pool->shared, &state);
  clean_exhaustion_required =
      !state.solved_published &&
      !state.stop_requested &&
      !state.fatal_error;
  pthread_mutex_lock(&pool->shared.queue_mutex);
  for (i = 0U; i < (size_t)pool->shared.worker_count; i++) {
    if (pool->contexts[i].published_helper_group ||
        pool->contexts[i].helper_preparation_active) {
      logerr("[index-shard] helper state remained active "
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

  if (helper_quiescence_valid) {
    index_shard_affinity_commit_pass(
        &pool->shared,
        affinity_domain,
        affinity_slots);
  } else {
    logerr("[index-shard] affinity ownership not committed "
           "after helper quiescence failure\n");
  }
  logverb("[index-shard] affinity-pass generation=%lu "
          "scale_index=%zu familiar=%llu fallback=%llu "
          "reassigned=%llu unclaimed=%zu committed=%i\n",
          pool->generation,
          bp->engine_scale_index,
          pool->shared.affinity_claims,
          pool->shared.fallback_claims,
          pool->shared.affinity_reassignments,
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

  if (worker_cancelled) {
    pthread_mutex_lock(&pool->shared.limit_mutex);
    bp->cancelled = TRUE;
    pthread_mutex_unlock(&pool->shared.limit_mutex);
  }

  /*
   * active_workers reached zero before reduce_online() returned, so every
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

  if (rc && status != INDEX_SHARD_SOLVE_TERMINAL_FAILURE) {
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
   * index_shard_pool_reduce_online() has returned and every participating
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

  logverb("[index-shard] context-pass generation=%lu candidates=%zu "
          "outer_workers=%i prepare_work_wall_sum=%.6f "
          "prepare_max=%.6f cleanup_work_wall_sum=%.6f "
          "cleanup_max=%.6f\n",
          pool->generation,
          nindexes,
          pool->shared.worker_count,
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
          "verify_calls=%llu hypothesis_wave_elapsed_sum=%.6f "
          "batches=%llu completed=%llu stopped=%llu "
          "batch_failed=%llu hypotheses=%llu executed=%llu reduced=%llu "
          "task_ranges=%llu tasks_executed=%llu submitted=%llu "
          "inline=%llu parallel_batches=%llu observed_parallel=%llu "
          "parallel_hypotheses=%llu alloc_failures=%llu "
          "search_failures=%llu max_batch=%zu max_tasks=%zu "
          "max_parallel=%zu ab_blocks_planned=%llu "
          "ab_blocks_retired=%llu ab_blocks_owner=%llu "
          "ab_segments_retired=%llu "
          "ab_segment_payload_bytes=%llu ab_pairs=%llu "
          "ab_combinations=%llu intra_pair_splits=%llu "
          "max_pair_combinations=%llu ab_plan_wall=%.6f "
          "helper_tasks=%llu helper_combinations=%llu "
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
          solver_profile.hypothesis_wave_wall_seconds,
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
          solver_profile.ab_blocks_planned,
          solver_profile.ab_blocks_retired,
          solver_profile.ab_blocks_owner,
          solver_profile.ab_segments_retired,
          solver_profile.ab_segment_payload_bytes,
          solver_profile.ab_pairs_planned,
          solver_profile.ab_combinations_planned,
          solver_profile.ab_intra_pair_splits,
          solver_profile.ab_max_pair_combinations,
          solver_profile.ab_planning_wall_seconds,
          solver_profile.ab_helper_tasks,
          solver_profile.ab_helper_combinations,
          solver_profile.hypothesis_order_hash,
          solver_profile.kd_result_order_hash,
          solver_profile.candidate_order_hash);

  logverb("[index-shard] pass-detail candidates=%zu reduced=%zu "
          "hit_total_cpu_limit=%i hit_total_wall_limit=%i "
          "cancelled=%i rc=%i status=%i "
          "master_committed=%i\n",
          nindexes,
          pass_metrics.reduced,
          bp->hit_total_cpulimit,
          bp->hit_total_timelimit,
          bp->cancelled,
          rc,
          (int)status,
          state.master_committed);

  logverb("[index-shard] mmap-policy "
         "policy=%s effective=%s pass=%u "
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
  pool->shared.preferred_owner = NULL;
  pool->shared.preferred_next = NULL;
  pool->shared.preferred_head = NULL;
  pthread_mutex_unlock(&pool->shared.queue_mutex);
  free(preferred_tail);
  free(preferred_head);
  free(preferred_next);
  free(preferred_owner);
  free(affinity_slots);
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
