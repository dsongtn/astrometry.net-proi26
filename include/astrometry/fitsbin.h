/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef FITSBIN_H
#define FITSBIN_H

#include <stdio.h>
#include <sys/stat.h>

#include "astrometry/anqfits.h"
#include "astrometry/bl.h"
#include "astrometry/an-bool.h"
/**
 "fitsbin" is our abuse of FITS binary tables to hold raw binary data,
 *without endian flips*, by storing the data as characters/bytes.
 This has the advantage that they can be directly mmap()'d, but of
 course means that they aren't endian-independent.  We accept that
 tradeoff in the interest of speed and the recognition that x86 is
 pretty much all we care about.


 Typical usage patterns:

 -Reading:
 fitsbin_open();
 fitsbin_add_chunk();
 fitsbin_add_chunk();
 ...
 fitsbin_read();
 ...
 fitsbin_close();

 OR:
 fb = fitsbin_open();
 fitsbin_chunk_init(&chunk);
 chunk.tablename = "hello";
 fitsbin_read_chunk(fb, &chunk);
 // chunk.data;
 //NO fitsbin_add_chunk(fb, &chunk);
 fitsbin_close();

 -Writing:
 fitsbin_open_for_writing();
 fitsbin_add_chunk();
 fitsbin_add_chunk();
 ...
 fitsbin_write_primary_header();
 ...
 fitsbin_write_chunk_header();
 fitsbin_write_items();
 ...
 fitsbin_fix_chunk_header();

 fitsbin_write_chunk_header();
 fitsbin_write_items();
 ...
 fitsbin_fix_chunk_header();
 ...
 fitsbin_fix_primary_header();
 fitsbin_close();

 OR:

 fb = fitsbin_open_for_writing();
 fitsbin_write_primary_header();

 fitsbin_chunk_init(&chunk);
 chunk.tablename = "whatever";
 chunk.data = ...;
 chunk.itemsize = 4;
 chunk.nrows = 1000;
 fitsbin_write_chunk(fb, &chunk);
 fitsbin_chunk_clean(&chunk);

 fitsbin_fix_primary_header();
 fitsbin_close(fb);

 */

struct fitsbin_t;
// mmap access policy for solver index tables.
typedef enum fitsbin_mmap_advice {
    FITSBIN_MMAP_ADVICE_NORMAL = 0,
    FITSBIN_MMAP_ADVICE_RANDOM = 1
} fitsbin_mmap_advice_t;

/*
 * Region metadata distinguishes topology and payload for page planning and
 * telemetry. Payload remains the default so zero-initialized third-party
 * chunks preserve compatibility. Topology uses NORMAL; payload follows the
 * file/pass policy.
 */
typedef enum fitsbin_mmap_region {
    FITSBIN_MMAP_REGION_PAYLOAD = 0,
    FITSBIN_MMAP_REGION_TOPOLOGY = 1
} fitsbin_mmap_region_t;


struct fitsbin_chunk_t {
    char* tablename;

    // internal use: pointer to strdup'd name.
    char* tablename_copy;

    // The data (NULL if the table was not found)
    void* data;

    // The size of a single row in bytes.
    int itemsize;

    // The number of items (rows)
    int nrows;

    // abort if this table isn't found?
    int required;

    // Reading:
    int (*callback_read_header)(struct fitsbin_t* fb, struct fitsbin_chunk_t* chunk);
    void* userdata;

    qfits_header* header;

    // Writing:
    off_t header_start;
    off_t header_end;

    // on output, force a type other than A?
    tfits_type forced_type;

    // Internal use:
    // The mmap'ed address
    char* map;
    // The mmap'ed size.
    size_t mapsize;

    /*
     * Recent bounded-delivery sequence for each exact VMA page. These
     * entries contain no payload data and expire after one service window.
     */
    unsigned int* payload_page_sequences;

    // Access-pattern class used when applying mmap advice.
    fitsbin_mmap_region_t mmap_region;

    /*
     * Exact file interval occupied by the table payload.  Unlike map/mapsize,
     * this excludes mmap alignment and FITS padding. The payload provider uses
     * it to translate proven mapped addresses into bounded page intervals and
     * test reads.
     */
    off_t data_file_offset;
    size_t data_file_size;
};
typedef struct fitsbin_chunk_t fitsbin_chunk_t;


struct fitsbin_t {
    char* filename;

    anqfits_t* fits;
    anbool owns_fits;

    bl* chunks;

    // Writing:
    FILE* fid;

    // only used for in_memory():
    anbool inmemory;
    bl* items;
    bl* extensions;

    // The primary FITS header
    qfits_header* primheader;
    off_t primheader_end;

    // for use when reading (not in_memory()): cache the tables in this FITS file.
    // ideally this would be pushed down to the qfits layer...
    qfits_table** tables;
    // number of extensions, include the primary; extensions < Next are valid.
    int Next;

    // for use by callback_read_header().
    void* userdata;
     // mmap policy applied to chunks read after opening this FITS file.
    fitsbin_mmap_advice_t mmap_advice;

    // Prevent repeated attempts after the kernel rejects an advice request.
    anbool mmap_advice_failed;

    // Cached system page size for bounded mmap prefetch requests.
    size_t mmap_page_size;

    // Enables bounded, caller-directed population of mapped payload pages.
    anbool mmap_prefetch_enabled;
    anbool mmap_prefetch_failed;

    /*
     * Captured source policy. This marks a source whose payload is already
     * fully resident; it does not change mapping or solver behavior itself.
     */
    anbool payload_fully_resident;

    /*
     * A separately opened buffered-I/O description for exact payload reads.
     * It receives POSIX_FADV_RANDOM but never backs a compute mapping, so its
     * file-description advice cannot change the mmap VMA policy.
     */
    int payload_fd;
    anbool payload_fd_initialized;
    anbool payload_fd_failed;

    /* Lock-free, per-file payload-I/O telemetry. */
    unsigned long long payload_read_calls;
    unsigned long long payload_read_batches;
    unsigned long long payload_read_logical_bytes;
    unsigned long long payload_read_pages;
    unsigned long long payload_read_bytes;
    unsigned long long payload_read_nanoseconds;
    unsigned long long payload_warm_calls;
    unsigned long long payload_warm_ranges;
    unsigned long long payload_warm_bytes;
    unsigned long long payload_warm_nanoseconds;
    unsigned long long payload_cache_hits;
    unsigned long long payload_cache_misses;
    unsigned long long payload_cache_evictions;
    unsigned long long payload_cache_allocations;
    unsigned long long payload_wait_nanoseconds;
    unsigned long long payload_failures;

    /*
     * Immutable identity captured from the descriptor that backs read-time
     * mappings. It remains valid after fitsbin_close_fd() so mapped data can
     * be identified without reopening a pathname or touching a closed FILE*.
     */
    anbool open_file_stat_valid;
    struct stat open_file_stat;
};
typedef struct fitsbin_t fitsbin_t;

// Initializes a chunk to default values
void fitsbin_chunk_init(fitsbin_chunk_t* chunk);

// Frees contents of this chunk.
void fitsbin_chunk_clean(fitsbin_chunk_t* chunk);

// clean + init
void fitsbin_chunk_reset(fitsbin_chunk_t* chunk);

char* fitsbin_get_filename(const fitsbin_t* fb);

// Reading: returns a new copy of the given FITS extension header.
// (-> *qfits_get_header)
qfits_header* fitsbin_get_header(const fitsbin_t* fb, int ext);

// Reading: how many extensions in this file?  (-> *qfits_query_n_ext)
int fitsbin_n_ext(const fitsbin_t* fb);

fitsbin_t* fitsbin_open(const char* fn);

fitsbin_t* fitsbin_open_fits(anqfits_t* fits);

fitsbin_t* fitsbin_open_for_writing(const char* fn);

fitsbin_t* fitsbin_open_in_memory(void);

int fitsbin_close_fd(fitsbin_t* fb);

/**
 Configures the mmap policy for solver index data.

 A W2+ shard worker supplies RANDOM before mapping every index chunk. A caller
 without shard-local policy retains the original serial NORMAL behavior.
 Bounded exact page population is advisory and does not change either mapping
 policy.
 */
int fitsbin_configure_index_mmap(fitsbin_t* fb);

/**
 Requests population of the mapped pages covering this data range.

 The request is a no-op unless bounded mmap population is enabled for the
 fitsbin.
 */
int fitsbin_prefetch_data(fitsbin_t* fb, const void* data, size_t size);

typedef struct fitsbin_prefetch_range {
    const void* data;
    size_t size;
} fitsbin_prefetch_range_t;

typedef struct fitsbin_payload_io_ticket
    fitsbin_payload_io_ticket_t;

typedef anbool (*fitsbin_payload_io_cancel_check_fn)(void* opaque);

/*
 * Construct one complete mapped-range plan on an I/O lane in provider-owned
 * storage. The callback is bounded and nonblocking: it must not enter or wait
 * on the payload provider. Return one for a complete nonempty plan, zero for
 * a complete logical operation requiring no mapped population, and minus one
 * for failure. Cancellation is authoritative only through cancelled().
 */
typedef int (*fitsbin_payload_io_plan_fn)(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count);

/*
 * Loader dequeue order is demand, current-index preparation, then
 * speculation. Non-demand admission leaves one job and one quarter of the
 * in-flight byte ceiling available for demand.
 */
typedef enum fitsbin_payload_io_priority {
    FITSBIN_PAYLOAD_IO_PRIORITY_DEMAND = 0,
    FITSBIN_PAYLOAD_IO_PRIORITY_CURRENT,
    FITSBIN_PAYLOAD_IO_PRIORITY_SPECULATIVE
} fitsbin_payload_io_priority_t;

typedef enum fitsbin_payload_io_submit_status {
    FITSBIN_PAYLOAD_IO_SUBMIT_ERROR = -1,
    FITSBIN_PAYLOAD_IO_SUBMIT_UNAVAILABLE = 0,
    FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED = 1,
    FITSBIN_PAYLOAD_IO_SUBMIT_READY = 2
} fitsbin_payload_io_submit_status_t;

typedef struct fitsbin_payload_io_stats {
    unsigned long long read_calls;
    unsigned long long read_batches;
    unsigned long long read_logical_bytes;
    unsigned long long read_pages;
    unsigned long long read_bytes;
    unsigned long long read_nanoseconds;
    unsigned long long warm_calls;
    unsigned long long warm_ranges;
    unsigned long long warm_bytes;
    unsigned long long warm_nanoseconds;
    unsigned long long cache_hits;
    unsigned long long cache_misses;
    unsigned long long cache_evictions;
    unsigned long long cache_allocations;
    unsigned long long wait_nanoseconds;
    unsigned long long failures;
} fitsbin_payload_io_stats_t;

/*
 * Configure the process-wide I/O admission width. Production mapped-page
 * population and focused buffered-read tests share this bounded mechanism.
 */
void fitsbin_payload_io_configure_workers(int worker_count);

/*
 * Start or stop the bounded payload loader. The service is advisory: callers
 * retain their original synchronous path when startup or submission fails.
 */
int fitsbin_payload_io_service_start(int lane_count);
void fitsbin_payload_io_service_stop(void);

/* Return the live loader width, or zero while the service is unavailable. */
int fitsbin_payload_io_service_width(void);

/* Return nonzero when mapped-page population is compiled into the loader. */
int fitsbin_payload_io_mapped_population_supported(void);

/*
 * One process-wide scheduler wakeup for terminal payload tickets. The
 * immutable completion ID is unique for the lifetime of the process and lets
 * a scheduler wake only the task that owns the completed ticket. The notifier
 * runs after the payload mutex is released, must remain bounded, and must not
 * perform payload I/O or wait for the payload service. It may briefly
 * synchronize with scheduler state. Clearing waits for active notifier calls
 * and must not be called by the notifier itself.
 */
typedef void (*fitsbin_payload_io_completion_notify_fn)(
    void* opaque,
    unsigned long long completion_id);

int fitsbin_payload_io_set_completion_notifier(
    fitsbin_payload_io_completion_notify_fn notify,
    void* opaque);
int fitsbin_payload_io_clear_completion_notifier(
    fitsbin_payload_io_completion_notify_fn notify,
    void* opaque);

/*
 * Cheap advisory predicate for optional work that may compete with payload
 * demand. It is nonzero while all reader credits are occupied or any demand
 * reader is waiting. A zero result grants no credit and may become stale
 * immediately; speculative callers must remain optional and bounded.
 */
int fitsbin_payload_io_demand_busy(void);

/*
 * Optional worker callback used only while a buffered demand reader has no
 * payload credit. It may execute one bounded, index-free helper package; a
 * nonzero return requests one immediate credit/helper retry. The optional
 * stop callback is sampled while a submitted ticket is pending; a nonzero
 * result cooperatively cancels that ticket.
 */
typedef int (*fitsbin_payload_io_wait_helper_fn)(void* opaque);
typedef int (*fitsbin_payload_io_stop_check_fn)(void* opaque);

int fitsbin_payload_io_set_thread_wait_helper(
    fitsbin_payload_io_wait_helper_fn helper,
    fitsbin_payload_io_stop_check_fn stop_check,
    void* opaque);
void fitsbin_payload_io_clear_thread_wait_helper(void);

/*
 * Count payload-credit waiters with a worker-local callback and wake them
 * after bounded helper work is published.
 */
size_t fitsbin_payload_io_wait_helper_count(void);
void fitsbin_payload_io_notify_wait_helpers(void);

/*
 * Delimit a production helper window. Nested windows share one reserved
 * reader slot so demand stays bounded while helper work is runnable.
 */
void fitsbin_payload_io_begin_helper_window(void);
void fitsbin_payload_io_end_helper_window(void);

typedef struct fitsbin_pread_range {
    const void* data;
    size_t size;
    size_t logical_size;
    void* destination;
} fitsbin_pread_range_t;

#define FITSBIN_PREAD_RANGE_LIMIT 16U
/* Direct async tickets own up to this many prepared range records. */
#define FITSBIN_PREAD_ASYNC_RANGE_LIMIT 256U
/* Mapped population accepts this many logical ranges before page deduplication. */
#define FITSBIN_MMAP_PREFETCH_RANGE_LIMIT 640U

/*
 * Resolve the file-page covering interval for one exact mapped request.
 * Alignment is calculated from immutable file offsets and clipped to the
 * containing FITS payload. The returned mapped address and file offset name
 * the same cover; the mapped address is only a token and is never
 * dereferenced or faulted by this function.
 */
int fitsbin_mapped_range_page_cover(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    const void** cover_data,
    size_t* cover_size,
    off_t* cover_file_offset,
    size_t* exact_offset);

/*
 * Read one through FITSBIN_PREAD_RANGE_LIMIT fully validated mapped ranges
 * under one demand-I/O credit. All ranges are resolved before I/O begins. A
 * failure invalidates the whole batch from the caller's perspective.
 * Disjoint destinations may be read in file-offset order. Any overlapping
 * destinations retain caller order.
 */
int fitsbin_pread_mapped_ranges(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count);

/*
 * Submit a fully validated direct read into caller-owned unpublished storage.
 * All mapped ranges are resolved before queueing and no mapped pointer is
 * retained. The complete physical byte count must fit byte_budget.
 *
 * Return 1 when queued, zero when bounded service capacity is unavailable,
 * and -1 for invalid input or failed preparation. After a return of 1, every
 * destination must remain allocated and inaccessible to the caller until
 * ticket_wait(), ticket_cancel_and_wait(), or a terminal ticket_poll(). Only
 * a positive collected result permits publication; failure or cancellation
 * invalidates the whole batch even when some destination bytes were written.
 * Disjoint destinations may be read in file-offset order; overlapping
 * destinations retain caller order. The source fitsbin must remain live
 * through collection so operation counters can be applied.
 *
 * On a zero return, errno is ENOTSUP when the exact source candidate disables
 * asynchronous direct transport for a nonresident sparse payload mapping,
 * ENODEV when the service is unavailable,
 * EAGAIN when current queue occupancy prevents admission, and E2BIG when the
 * unchanged ticket can never fit the service byte ceiling.  Callers must keep
 * their authoritative mapped fallback available for every zero return.
 */
int fitsbin_pread_mapped_ranges_submit(
    fitsbin_t* fb,
    const fitsbin_pread_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_payload_io_priority_t priority,
    fitsbin_payload_io_ticket_t** ticket);

/*
 * Read the exact mapped payload bytes through the dedicated buffered-I/O
 * description.  The mapped address is used only to resolve the immutable file
 * offset; bytes are returned in caller storage.
 *
 * Return 0 on success and -1 on identity, range, allocation, or I/O failure.
 */
int fitsbin_pread_mapped_range(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    void* destination);

/*
 * Prepare a complete bounded group of mapped ranges through exact buffered
 * reads. Ranges are page aligned, sorted, deduplicated, and merged only when
 * overlapping or adjacent. Oversized, unresolvable, or over-budget groups
 * fail as a unit; callers must split larger plans explicitly. Correctness
 * never depends on this operation. A fully resident source returns zero
 * without issuing I/O or changing warm counters.
 */
int fitsbin_prefetch_ranges(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget);

/*
 * Submit one complete current-index mapped-page population to the bounded
 * loader. The source must first be initialized by
 * fitsbin_configure_index_mmap(). The ticket copies validated page-aligned
 * spans but borrows the fitsbin owner and its mappings.
 *
 * Return FITSBIN_PAYLOAD_IO_SUBMIT_QUEUED with a ticket when queued,
 * FITSBIN_PAYLOAD_IO_SUBMIT_READY without a ticket when the exact live-mapping
 * completion record already covers every requested page, zero when the
 * optional service or bounded capacity is unavailable, and -1 for an invalid
 * or failed preparation. READY does not pin pages; the native mapped read
 * remains authoritative if the kernel has reclaimed one. The caller must keep
 * every source mapping live through blocking or polled collection before
 * closing the source fitsbin. A fully resident source returns zero without
 * creating a ticket. On a zero return caused by the service, errno
 * distinguishes ENODEV, transient admission refusal EAGAIN, and permanent
 * size refusal E2BIG.
 */
int fitsbin_prefetch_ranges_submit(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget,
    fitsbin_payload_io_ticket_t** ticket);

/*
 * Submit a bounded mapped-population ticket whose exact ranges are planned on
 * the I/O lane. byte_budget is reserved at admission; actual mapped bytes are
 * recorded after the callback completes. plan_opaque and every object it can
 * reach must remain alive and exclusively owned by the ticket until terminal
 * collection. An empty successful plan collects as a positive logical result
 * and does not increment mapped-warm counters. Refusal and failure leave the
 * native mapped path authoritative.
 */
int fitsbin_prefetch_ranges_planned_submit(
    fitsbin_t* fb,
    fitsbin_payload_io_plan_fn plan,
    void* plan_opaque,
    size_t byte_budget,
    fitsbin_payload_io_ticket_t** ticket);

/*
 * Wait for one ticket and apply operation-specific I/O counters to the
 * still-live fitsbin. A failed or cancelled operation leaves the original
 * synchronous solver path valid. A direct ticket returns its range count only
 * after every destination is complete. destroy() is valid only after blocking
 * or polled collection. The exact source passed at submission must remain
 * live and must be passed again when collecting the ticket.
 */
int fitsbin_payload_io_ticket_wait(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket);
int fitsbin_payload_io_ticket_cancel_and_wait(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket);

/*
 * Poll one ticket without waiting. Return one when terminal and place the
 * ordinary wait result in result_out, zero while pending, and minus one for
 * invalid or already-collected input. A terminal poll applies counters and
 * must be followed by exactly one checked destroy. A terminal failure is
 * reported as a successful poll with result_out set to -1 and errno set to
 * the I/O failure.
 */
int fitsbin_payload_io_ticket_poll(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t* ticket,
    int* result_out);

/*
 * Poll a ticket and, when terminal, atomically transfer its result and destroy
 * its storage before returning success. A zero return leaves *ticket_io
 * unchanged. A positive return sets *ticket_io to NULL and stores the
 * provider result in *result_out. A negative return is an ownership or API
 * integrity error and leaves *ticket_io unchanged; it is not a retryable
 * capacity refusal. This is the preferred operation for a scheduler that must
 * not release its source lease before ticket ownership is gone.
 */
int fitsbin_payload_io_ticket_poll_and_destroy(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t** ticket_io,
    int* result_out);

/*
 * Take exclusive owner responsibility for one ticket. A registered waiter is
 * drained through its final result accounting; otherwise pending work is
 * cancelled and drained. Then destroy the ticket and set *ticket_io to NULL.
 * Return one after terminal transfer and minus one for invalid ownership. The
 * source must remain live until this function returns.
 */
int fitsbin_payload_io_ticket_drain_and_destroy(
    fitsbin_t* fb,
    fitsbin_payload_io_ticket_t** ticket_io,
    int* result_out);

/*
 * Request cancellation without waiting or applying counters. Return one for
 * a new request, zero when already requested or terminal, and minus one for
 * invalid input. Queued tickets become terminal immediately; running tickets
 * complete cancellation at their next bounded polling point.
 */
int fitsbin_payload_io_ticket_cancel_async(
    fitsbin_payload_io_ticket_t* ticket);

/*
 * Return the immutable nonzero completion ID assigned at successful
 * submission, or zero for invalid or unpublished input. The caller must
 * still own the live ticket.
 */
unsigned long long fitsbin_payload_io_ticket_completion_id(
    const fitsbin_payload_io_ticket_t* ticket);

/*
 * Destroy a terminal, collected ticket. A failure returns minus one with
 * errno set to EBUSY for pending or actively waited tickets, or EAGAIN when
 * terminal counters have not been collected. Ticket operations, collection,
 * and destruction require exclusive caller serialization. The legacy void
 * wrapper preserves its original silent no-op behavior.
 */
int fitsbin_payload_io_ticket_destroy_checked(
    fitsbin_payload_io_ticket_t* ticket);
void fitsbin_payload_io_ticket_destroy(
    fitsbin_payload_io_ticket_t* ticket);

/*
 * Populate every page-aligned mapped range without copying payload bytes.
 * Ranges are sorted, deduplicated, and merged within each mapping. The whole
 * aligned plan is preflighted against byte_budget and the process-wide byte
 * ceiling before any page is touched. This operation never changes the base
 * mapping advice. Refusal or failure therefore leaves the caller's NORMAL or
 * RANDOM mapping policy intact.
 *
 * Return the number of fully populated spans, zero when the complete plan is
 * not applicable, or -1 on invalid input or kernel failure.
 */
int fitsbin_advise_mapped_ranges(
    fitsbin_t* fb,
    const fitsbin_prefetch_range_t* ranges,
    size_t range_count,
    size_t byte_budget);

/* Close only the dedicated payload reader; compute mappings remain valid. */
int fitsbin_close_payload_fd(fitsbin_t* fb);

/* Atomically take and reset this fitsbin's payload-I/O counters. */
void fitsbin_take_payload_io_stats(
    fitsbin_t* fb,
    fitsbin_payload_io_stats_t* stats);

/**
 Resolves a requested data range to the actual mmap region containing it.

 On success:
   - map_base/map_size identify the exact mmap mapping.
   - range_start/range_size identify the requested interval clipped to that
     mapping.

 Return:
    1  containing mmap region found
    0  address is not backed by a mapped fitsbin chunk
   -1  invalid arguments
 */
int fitsbin_resolve_mapped_range(fitsbin_t* fb,
                                 const void* data,
                                 size_t size,
                                 const void** map_base,
                                 size_t* map_size,
                                 const void** range_start,
                                 size_t* range_size);

/*
 * Apply advice only to the page-aligned mapped interval containing the
 * requested data. This is used for access phases whose intent is narrower
 * than the containing FITS payload chunk, such as a sequential PERM sweep.
 *
 * Return:
 *   1  advice applied
 *   0  range is not mmap-backed
 *  -1  invalid input or kernel rejection
 */
int fitsbin_set_mmap_range_advice(
    fitsbin_t* fb,
    const void* data,
    size_t size,
    fitsbin_mmap_advice_t advice);

int fitsbin_switch_to_reading(fitsbin_t* fb);

int fitsbin_read(fitsbin_t* fb);

fitsbin_chunk_t* fitsbin_get_chunk(fitsbin_t* fb, int chunk);

off_t fitsbin_get_data_start(fitsbin_t* fb, fitsbin_chunk_t* chunk);

int fitsbin_n_chunks(fitsbin_t* fb);

/**
 Appends the given chunk -- makes a copy of the contents of "chunk" and
 returns a pointer to the stored location.
 */
fitsbin_chunk_t* fitsbin_add_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk);

/**
 Immediately tries to read a chunk.  If the chunk is not found, -1 is returned
 and the chunk is not added to this fitsbin's list.  If it's found, 0 is
 returned, a copy of the chunk is stored, and the results are placed in
 "chunk".
 */
int fitsbin_read_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk);

/*
 * Reads and validates a chunk's FITS/table headers without mapping or
 * copying its payload.  The chunk is not added to fb.
 */
int fitsbin_read_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk);

FILE* fitsbin_get_fid(fitsbin_t* fb);

int fitsbin_get_open_file_stat(
    const fitsbin_t* fb,
    struct stat* file_stat);

void fitsbin_stat_times(
    const struct stat* file_stat,
    time_t* mtime_seconds,
    long* mtime_nanoseconds,
    time_t* ctime_seconds,
    long* ctime_nanoseconds);

int fitsbin_close(fitsbin_t* fb);

qfits_header* fitsbin_get_primary_header(const fitsbin_t* fb);

void fitsbin_set_primary_header(fitsbin_t* fb, const qfits_header* hdr);

// (pads to FITS block size)
int fitsbin_write_primary_header(fitsbin_t* fb);

// (pads to FITS block size)
int fitsbin_fix_primary_header(fitsbin_t* fb);

qfits_header* fitsbin_get_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk);

int fitsbin_write_chunk(fitsbin_t* fb, fitsbin_chunk_t* chunk);

int fitsbin_write_chunk_flipped(fitsbin_t* fb, fitsbin_chunk_t* chunk,
                                int wordsize);

// (pads to FITS block size)
int fitsbin_write_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk);

// (pads to FITS block size)
int fitsbin_fix_chunk_header(fitsbin_t* fb, fitsbin_chunk_t* chunk);

int fitsbin_write_item(fitsbin_t* fb, fitsbin_chunk_t* chunk, void* data);

int fitsbin_write_items(fitsbin_t* fb, fitsbin_chunk_t* chunk, void* data, int N);


// direct FILE* output:

int fitsbin_write_primary_header_to(fitsbin_t* fb, FILE* fid);

int fitsbin_write_chunk_header_to(fitsbin_t* fb, fitsbin_chunk_t* chunk, FILE* fid);

int fitsbin_write_items_to(fitsbin_chunk_t* chunk, void* data, int N, FILE* fid);

int fitsbin_write_chunk_to(fitsbin_t* fb, fitsbin_chunk_t* chunk, FILE* fid);

/*
 * Index mmap policy.
 *
 * Fixed policies preserve one base mapping advice for the complete field.
 * Parallel production uses FIXED_RANDOM for every mapped index chunk, while
 * serial callers remain NORMAL. Bounded exact page population is additive and
 * policy-neutral.
 */
typedef enum fitsbin_mmap_policy {
    FITSBIN_MMAP_POLICY_FIXED_NORMAL = 0,
    FITSBIN_MMAP_POLICY_FIXED_RANDOM,
    FITSBIN_MMAP_POLICY_ADAPTIVE
} fitsbin_mmap_policy_t;

typedef struct fitsbin_mmap_advice_state {
    fitsbin_mmap_policy_t policy;
    fitsbin_mmap_advice_t effective_advice;

    /*
     * Number of fully quiesced passes recorded for the current field.
     */
    unsigned int pass_number;

    /*
     * Adaptive-policy evidence and transitions for the current field.
     */
    unsigned int completed_clean_unsolved_passes;
    unsigned int transition_count;
} fitsbin_mmap_advice_state_t;

fitsbin_mmap_policy_t fitsbin_mmap_policy_parse(const char* value);

/*
 * Return the automatic production policy. Runtime environment overrides are
 * intentionally not part of the public UX.
 */
fitsbin_mmap_policy_t fitsbin_get_configured_mmap_policy(void);

const char* fitsbin_mmap_policy_name(
    fitsbin_mmap_policy_t policy);

const char* fitsbin_mmap_advice_name(
    fitsbin_mmap_advice_t advice);

const char* fitsbin_mmap_region_name(
    fitsbin_mmap_region_t region);

fitsbin_mmap_advice_t fitsbin_get_mmap_advice(
    const fitsbin_t* fb);

/*
 * Topology chunks use NORMAL. Payload chunks follow file/pass advice.
 */
fitsbin_mmap_advice_t fitsbin_get_chunk_mmap_advice(
    const fitsbin_t* fb,
    const fitsbin_chunk_t* chunk);

void fitsbin_mmap_advice_state_init(
    fitsbin_mmap_advice_state_t* state,
    fitsbin_mmap_policy_t policy);

void fitsbin_mmap_advice_state_reset(
    fitsbin_mmap_advice_state_t* state);

fitsbin_mmap_advice_t fitsbin_mmap_advice_state_begin_pass(
    const fitsbin_mmap_advice_state_t* state);

/*
 * Records one fully quiesced pass.
 *
 * Returns TRUE only when adaptive policy changes from RANDOM to NORMAL.
 */
anbool fitsbin_mmap_policy_complete_pass(
    fitsbin_mmap_advice_state_t* state,
    anbool pass_completed,
    anbool pass_exhaustive,
    anbool pass_solved,
    anbool pass_cancelled,
    int pass_rc,
    int pass_status);

/*
 * Worker-local pass advice used while index components are opened and mmaped.
 */
void fitsbin_mmap_set_thread_advice(
    fitsbin_mmap_advice_t advice);

void fitsbin_mmap_clear_thread_advice(void);

fitsbin_mmap_advice_t fitsbin_mmap_current_advice(void);

/* TRUE only while an outer shard worker owns the current index open. */
anbool fitsbin_mmap_thread_advice_active(void);

/*
 * Mark index components opened by this thread as fully resident. The marker
 * is captured by fitsbin_configure_index_mmap() and has no effect by itself.
 */
void fitsbin_payload_set_thread_full_resident(void);
void fitsbin_payload_clear_thread_full_resident(void);
anbool fitsbin_payload_is_fully_resident(
    const fitsbin_t* fb);

/*
 * Changes the logical advice and optionally reapplies it to existing chunks.
 */
int fitsbin_set_mmap_advice(
    fitsbin_t* fb,
    fitsbin_mmap_advice_t advice,
    anbool reapply_existing);

#endif
