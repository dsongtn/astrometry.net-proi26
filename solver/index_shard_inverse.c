/*
 * Private implementation module for the index-shard subsystem.
 * See index_shard_private.h for ownership and lock-order invariants.
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

#include "index_shard_private.h"
#include "astrometry/bl.h"
#include "astrometry/errors.h"
#include "astrometry/log.h"
#include "astrometry/tic.h"
#include "astrometry/fitsbin.h"
#include "astrometry/fitsioutils.h"
/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/*
 * Bounded inverse-permutation cache and asynchronous preparation.
 *
 * This module owns the bounded inverse-permutation cache and its leases.
 */

size_t index_shard_inverse_cache_budget(void) {
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

void index_shard_inverse_cache_attach(
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

void index_shard_inverse_cache_release(
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

void index_shard_inverse_cache_destroy(
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
