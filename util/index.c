/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include "index.h"
#include "log.h"
#include "errors.h"
#include "ioutils.h"
#include "healpix.h"
#include "tic.h"

#include "anqfits.h"
#include "qfits_rw.h"
#include "starutil.h"

#include <pthread.h>

static index_residency_t* active_index_residency = NULL;
static size_t active_index_residency_readers = 0U;
static pthread_mutex_t active_index_residency_mutex =
    PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t active_index_residency_cond =
    PTHREAD_COND_INITIALIZER;

int index_bind_residency_service(index_residency_t* service) {
    int status = 0;

    if (!service) {
        return -1;
    }
    pthread_mutex_lock(&active_index_residency_mutex);
    if (active_index_residency) {
        status = -1;
    } else {
        active_index_residency = service;
    }
    pthread_mutex_unlock(&active_index_residency_mutex);
    return status;
}

void index_unbind_residency_service(index_residency_t* service) {
    pthread_mutex_lock(&active_index_residency_mutex);
    if (active_index_residency == service) {
        active_index_residency = NULL;
    }
    while (active_index_residency_readers) {
        pthread_cond_wait(
            &active_index_residency_cond,
            &active_index_residency_mutex);
    }
    pthread_mutex_unlock(&active_index_residency_mutex);
}

anbool index_residency_service_active(void) {
    anbool active;

    pthread_mutex_lock(&active_index_residency_mutex);
    active = active_index_residency != NULL;
    pthread_mutex_unlock(&active_index_residency_mutex);
    return active;
}

static index_residency_t* index_borrow_residency_service(void) {
    index_residency_t* service;

    pthread_mutex_lock(&active_index_residency_mutex);
    service = active_index_residency;
    if (service) {
        active_index_residency_readers++;
    }
    pthread_mutex_unlock(&active_index_residency_mutex);
    return service;
}

static void index_return_residency_service(
    index_residency_t* service) {
    if (!service) {
        return;
    }
    pthread_mutex_lock(&active_index_residency_mutex);
    if (active_index_residency_readers) {
        active_index_residency_readers--;
    }
    if (!active_index_residency_readers) {
        pthread_cond_broadcast(
            &active_index_residency_cond);
    }
    pthread_mutex_unlock(&active_index_residency_mutex);
}

static anbool index_try_resident_open(index_t* index) {
    index_residency_t* service;
    index_residency_handle_t* handle = NULL;
    index_residency_result_t result;
    const char* backing_path;
    anqfits_t* resident_fits;

    if (!index || !index->indexfn ||
        index->residency || index->residency_source) {
        return FALSE;
    }
    service = index_borrow_residency_service();
    if (!service) {
        return FALSE;
    }
    result = index_residency_acquire(
        service,
        index->indexfn,
        &handle);
    index_return_residency_service(service);
    if (result == INDEX_RESIDENCY_SOURCE_LEASE &&
        handle) {
        index->residency_source = handle;
        return FALSE;
    }
    if (result != INDEX_RESIDENCY_ACCEPTED || !handle) {
        index_residency_release(handle);
        return FALSE;
    }
    backing_path = index_residency_handle_path(handle);
    resident_fits = backing_path ? anqfits_open(backing_path) : NULL;
    if (!resident_fits) {
        index_residency_release(handle);
        return FALSE;
    }
    if (index->fits) {
        anqfits_close(index->fits);
    }
    index->fits = resident_fits;
    index->residency = handle;
    return TRUE;
}

int index_get_source_file_stat(
    const index_t* index,
    struct stat* identity) {
    const struct stat* resident_identity;

    if (!index || !identity) {
        return -1;
    }
    if (index->residency) {
        resident_identity =
            index_residency_handle_source_stat(index->residency);
        if (!resident_identity) {
            return -1;
        }
        *identity = *resident_identity;
        return 0;
    }
    if (index->quads && index->quads->fb &&
        !fitsbin_get_open_file_stat(
            index->quads->fb,
            identity)) {
        return 0;
    }
    return index->indexfn ? stat(index->indexfn, identity) : -1;
}

anbool index_overlaps_scale_range(index_t* meta,
                                  double quadlo, double quadhi) {
    anbool rtn =
        !((quadlo > meta->index_scale_upper) ||
          (quadhi < meta->index_scale_lower));
    debug("index_overlaps_scale_range: index %s has quads [%g, %g] arcsec; image has quads [%g, %g] arcsec.  In range? %s\n",
          meta->indexname, meta->index_scale_lower, meta->index_scale_upper, quadlo, quadhi, rtn ? "yes" : "no");
    return rtn;
}

anbool index_is_within_range(index_t* meta, double ra, double dec, double radius_deg) {
    if (meta->healpix == -1) {
        // allsky; tautology
        return TRUE;
    }
    return (healpix_distance_to_radec(meta->healpix, meta->hpnside, ra, dec, NULL) <= radius_deg);
}

int index_get_quad_dim(const index_t* index) {
    return quadfile_dimquads(index->quads);
}

int index_get_code_dim(const index_t* index) {
    return dimquad2dimcode(index_get_quad_dim(index));
}

int index_nquads(const index_t* index) {
    return quadfile_nquads(index->quads);
}

int index_nstars(const index_t* index) {
    return startree_N(index->starkd);
}

// Returns a newly-allocated string containing the filename, or NULL if not found.
// Tries the index name, or name + ".fits" as a filename.
static char* get_filename(const char* indexname) {
    char* fits;
    if (file_readable(indexname)) {
        logverb("Index name \"%s\" is readable, using as index filename\n", indexname);
        return strdup(indexname);
    }
    asprintf_safe(&fits, "%s.fits", indexname);
    if (file_readable(fits)) {
        // assume single-file index.
        logverb("Index name \"%s\" with .fits suffix, \"%s\", is readable, using as index filename.\n", indexname, fits);
        return fits;
    }
    free(fits);
    return NULL;
}

char* index_get_quad_filename(const char* indexname) {
    if (!index_is_file_index(indexname))
        return NULL;
    return get_filename(indexname);
}

char* index_get_qidx_filename(const char* indexname) {
    char* indexfn;
    char* qidxfn = NULL;
    if (!index_is_file_index(indexname))
        return NULL;
    indexfn = get_filename(indexname);
    if (ends_with(indexfn, ".fits")) {
        asprintf_safe(&qidxfn, "%.*s.qidx.fits", (int)(strlen(indexfn)-5), indexfn);
    } else {
        asprintf_safe(&qidxfn, "%s.qidx.fits", indexfn);
    }
    free(indexfn);
    return qidxfn;
}

anbool index_is_file_index(const char* filename) {
    char* indexfn;
    anbool rtn = TRUE;
    indexfn = get_filename(filename);
    if (!file_readable(indexfn)) {
        logdebug("Index file %s -> %s is not readable.", filename, indexfn);
        rtn = FALSE;
        goto finish;
    }
    if (!qfits_is_fits(indexfn)) {
        ERROR("Index file %s is not FITS.\n", indexfn);
        rtn = FALSE;
        goto finish;
    }
    /* This is a bit expensive...
     if (index_get_meta(filename, &meta)) {
     if (singlefile)
     ERROR("File %s does not contain an index.\n", quadfn);
     else
     ERROR("Files %s , %s , and %s do not contain an index.\n",
     quadfn, skdtfn, ckdtfn);
     rtn = FALSE;
     }
     */

 finish:
    free(indexfn);
    return rtn;
}

int index_get_meta(const char* filename, index_t* meta) {
    index_t* ind = index_load(filename, INDEX_ONLY_LOAD_METADATA, meta);
    if (!ind)
        return -1;
    return 0;
}

int index_get_missing_cut_params(int indexid, int* hpnside, int* nsweep,
                                 double* dedup, int* margin, char** pband) {
    // The 200-series indices use cut 100 (usnob)
    // The 500-series indices use cut 100 (usnob)
    // The 600-series indices use cut 300 (2mass)
    // The 700-series indices use cut 400 (usnob)
    int i = -1;
    int ns, n, marg;
    double dd;
    char* band;

    if ((indexid >= 200 && indexid < 220) ||
        (indexid >= 500 && indexid < 520)) {
        // Cut 100 params:
        int cut100hp[] = { 1760, 1245, 880, 622, 440, 312, 220, 156, 110, 78, 55, 39, 28, 20, 14, 10, 7, 5, 4, 3 };
        int cut100n[] = { 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 };
        double cut100dd[] = { 8, 8, 8, 8, 8, 9.6, 13.2, 18.0, 25.2, 36, 51, 72, 102, 144, 204, 288, 408, 600, 840, 1200 };
        int cut100margin = 5;

        i = indexid % 100;
        ns = cut100hp[i];
        n = cut100n[i];
        dd = cut100dd[i];
        marg = cut100margin;
        band = "R";

    } else if (indexid >= 602 && indexid < 620) {
        // Cut 300 params:
        int cut300hp[] = { 0, 0, 880, 624, 440, 312, 220, 156, 110, 78, 56, 40, 28, 20, 14, 10, 8, 6, 4, 4 };
        int cut300n = 10;
        double cut300dd = 8.0;
        //double cut300dd[] = { 8, 8, 8, 8, 8, 9.6, 13.2, 18.0, 25.2, 36, 51, 72, 102, 144, 204, 288, 408, 600, 840, 1200 };
        int cut300margin = 10;

        i = indexid % 100;
        ns = cut300hp[i];
        n = cut300n;
        dd = cut300dd;
        marg = cut300margin;
        band = "J";

    } else if (indexid >= 700 && indexid < 720) {
        // Cut 400 params:
        // (cut 400 used cut 200 as input: it had dedup=8, and n=6,10,10,...)
        int cut400hp[] = { 1760, 1246, 880, 624, 440, 312, 220, 156, 110, 78, 55, 39, 28, 20, 14, 10, 7, 5, 4, 3 };
        int cut400n[] = { 6, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
        double cut400dd = 8.0;
        //double cut400dd[] = { 8, 8, 8, 8, 8, 9.6, 13.2, 18.0, 25.2, 36, 51, 72, 102, 144, 204, 288, 408, 600, 840, 1200 };
        int cut400margin = 10;

        i = indexid % 100;
        ns = cut400hp[i];
        n = cut400n[i];
        dd = cut400dd;
        marg = cut400margin;
        band = "R";
    } else {
        return -1;
    }

    if (hpnside)
        *hpnside = ns;
    if (nsweep)
        *nsweep = n;
    if (dedup)
        *dedup = dd;
    if (margin)
        *margin = marg;
    if (pband)
        *pband = strdup(band);
    return 0;
}

static void get_cut_params(index_t* index) {
    index->index_jitter = startree_get_jitter(index->starkd);
    if (index->index_jitter == 0.0)
        index->index_jitter = DEFAULT_INDEX_JITTER;

    index->cutnside = startree_get_cut_nside(index->starkd);
    index->cutnsweep = startree_get_cut_nsweeps(index->starkd);
    index->cutdedup = startree_get_cut_dedup(index->starkd);
    index->cutband = strdup_safe(startree_get_cut_band(index->starkd));
    index->cutmargin = startree_get_cut_margin(index->starkd);

    // HACK - fill in values that are missing in old index files.
    {
        int *nside = NULL, *nsweep = NULL, *margin = NULL;
        char** band = NULL;
        double* dedup = NULL;

        if (index->cutnside == -1)
            nside = &(index->cutnside);
        if (index->cutnsweep == 0)
            nsweep = &(index->cutnsweep);
        if (index->cutmargin == -1)
            margin = &(index->cutmargin);
        if (index->cutdedup == 0)
            dedup = &(index->cutdedup);
        if (!index->cutband)
            band = &(index->cutband);

        index_get_missing_cut_params(index->indexid, nside, nsweep, dedup, margin, band);
    }

}

static void set_meta(index_t* index) {
    index->index_scale_upper = quadfile_get_index_scale_upper_arcsec(index->quads);
    index->index_scale_lower = quadfile_get_index_scale_lower_arcsec(index->quads);
    index->indexid = index->quads->indexid;
    index->healpix = index->quads->healpix;
    index->hpnside = index->quads->hpnside;
    index->dimquads = index->quads->dimquads;
    index->nquads = index->quads->numquads;
    index->nstars = index->quads->numstars;

    // This must get called after meta.indexid is set: otherwise we won't be
    // able to fill in values that are missing in old index files.
    get_cut_params(index);

    // check for CIRCLE field in ckdt header...
    index->circle = qfits_header_getboolean(index->codekd->header, "CIRCLE", 0);

    // New indexes are cooked such that cx < dx for all codes, but not
    // all of the old ones are like this.
    index->cx_less_than_dx = qfits_header_getboolean(index->codekd->header, "CXDX", FALSE);

    index->meanx_less_than_half = qfits_header_getboolean(index->codekd->header, "CXDXLT1", FALSE);
}

static int index_reload_internal(index_t* index, anbool metadata_only) {
    anbool full_resident = index && index->residency;

    if (full_resident)
        fitsbin_payload_set_thread_full_resident();
    // Read .skdt file...
    if (!index->starkd) {
        if (metadata_only) {
            index->starkd = startree_open_fits_metadata(index->fits);
        } else {
            index->starkd = startree_open_fits(index->fits);
        }
        if (!index->starkd) {
            ERROR("Failed to read star kdtree from file %s", index->indexfn);
            goto bailout;
        }
    }

    // Read .quad file...
    if (!index->quads) {
        if (metadata_only) {
            index->quads = quadfile_open_fits_metadata(index->fits);
        } else {
            index->quads = quadfile_open_fits(index->fits);
        }
        if (!index->quads) {
            ERROR("Failed to read quads from %s", index->indexfn);
            goto bailout;
        }
    }

    // Read .ckdt file...
    if (!index->codekd) {
        if (metadata_only) {
            index->codekd = codetree_open_fits_metadata(index->fits);
        } else {
            index->codekd = codetree_open_fits(index->fits);
        }
        if (!index->codekd) {
            ERROR("Failed to read code kdtree from file %s", index->indexfn);
            goto bailout;
        }
    }
    if (full_resident)
        fitsbin_payload_clear_thread_full_resident();
    return 0;

 bailout:
    if (full_resident)
        fitsbin_payload_clear_thread_full_resident();
    return -1;
}

static int index_reopen_original_backing(index_t* index) {
    if (!index || !index->residency || !index->indexfn) {
        return -1;
    }
    index_unload(index);
    if (index->fits) {
        anqfits_close(index->fits);
        index->fits = NULL;
    }
    index_residency_release(index->residency);
    index->residency = NULL;
    index->fits = anqfits_open(index->indexfn);
    return index->fits ? 0 : -1;
}

int index_dimquads(index_t* indx) {
    return indx->dimquads;
}

index_t* index_build_from(codetree_t* codekd, quadfile_t* quads, startree_t* starkd) {
    index_t* index = calloc(1, sizeof(index_t));
    index->codekd = codekd;
    index->quads = quads;
    index->starkd = starkd;
    set_meta(index);
    return index;
}

index_t* index_load(const char* indexname, int flags, index_t* dest) {
    index_t* allocd = NULL;

    if (flags & INDEX_ONLY_LOAD_METADATA)
        logverb("Loading metadata for %s...\n", indexname);

    if (!dest)
        allocd = dest = calloc(1, sizeof(index_t));
    else
        memset(dest, 0, sizeof(index_t));

    dest->indexname = strdup(indexname);

    dest->indexfn = get_filename(indexname);
    if (!dest->indexfn) {
        ERROR("Did not find file for index named %s", dest->indexname);
        goto bailout;
    }
    if (!(flags & INDEX_ONLY_LOAD_METADATA))
        (void)index_try_resident_open(dest);
    if (!dest->fits)
        dest->fits = anqfits_open(dest->indexfn);
    if (!dest->fits) {
        ERROR("Failed to open FITS file %s", dest->indexfn);
        goto bailout;
    }
    if (index_reload_internal(
            dest, flags & INDEX_ONLY_LOAD_METADATA)) {
        if (!dest->residency ||
            index_reopen_original_backing(dest) ||
            index_reload_internal(dest, FALSE)) {
            goto bailout;
        }
    }

    free(dest->indexname);
    dest->indexname = strdup(
        dest->residency ? dest->indexfn :
        quadfile_get_filename(dest->quads));
    set_meta(dest);

    logverb("Index scale: [%g, %g] arcmin, [%g, %g] arcsec\n",
            dest->index_scale_lower / 60.0, dest->index_scale_upper / 60.0,
            dest->index_scale_lower, dest->index_scale_upper);
    logverb("Index has %i quads and %i stars\n", dest->nquads, dest->nstars);

    if (!dest->circle) {
        ERROR("Code kdtree does not contain the CIRCLE header.");
        goto bailout;
    }

    if (flags & INDEX_ONLY_LOAD_METADATA) {
        index_unload(dest);
        // If we're using anqfits_t (dest->fits), keep that open for
        // fast reopening.  anqfits_t doesn't keep a FILE* or anything
        // open, so that's fine.
    }
    return dest;

 bailout:
    index_close(dest);
    free(allocd);
    return NULL;
}

int index_reload(index_t* index) {
    int status;

    if (!index) {
        return -1;
    }
    if (!index->codekd && !index->quads && !index->starkd) {
        (void)index_try_resident_open(index);
    }
    status = index_reload_internal(index, FALSE);
    if (!status || !index->residency) {
        return status;
    }
    if (index_reopen_original_backing(index)) {
        return -1;
    }
    return index_reload_internal(index, FALSE);
}

void index_unload(index_t* index) {
    if (index->starkd) {
        startree_close(index->starkd);
        index->starkd = NULL;
    }
    if (index->codekd) {
        codetree_close(index->codekd);
        index->codekd = NULL;
    }
    if (index->quads) {
        quadfile_close(index->quads);
        index->quads = NULL;
    }
}

int index_close_fds(index_t* ind) {
    kdtree_fits_t* io;
    int rc = 0;

    if (!ind || !ind->quads || !ind->quads->fb ||
        !ind->codekd || !ind->codekd->tree ||
        !ind->codekd->tree->io ||
        !ind->starkd || !ind->starkd->tree ||
        !ind->starkd->tree->io) {
        ERROR("Cannot close descriptors for an incomplete index");
        return -1;
    }
    if (fitsbin_close_fd(ind->quads->fb)) {
        ERROR("Failed to close an astrometry_net_data quadfile");
        rc = -1;
    }
    io = ind->codekd->tree->io;
    if (fitsbin_close_fd(io)) {
        ERROR("Failed to close an astrometry_net_data code kdtree");
        rc = -1;
    }
    io = (kdtree_fits_t*)ind->starkd->tree->io;
    if (fitsbin_close_fd(io)) {
        ERROR("Failed to close an astrometry_net_data star kdtree");
        rc = -1;
    }
    return rc;
}

static int index_close_one_payload_fd(
    fitsbin_t* fb,
    const char* component,
    const char* index_name) {
    fitsbin_payload_io_stats_t stats;
    int rc;

    if (!fb) {
        return 0;
    }
    fitsbin_take_payload_io_stats(fb, &stats);
    if (stats.read_calls ||
        stats.warm_calls ||
        stats.cache_hits ||
        stats.cache_misses ||
        stats.failures) {
        logverb(
            "[index-payload-io] index=%s component=%s "
            "direct_calls=%llu direct_batches=%llu "
            "logical_bytes=%llu aligned_bytes=%llu page_coverage=%llu "
            "direct_ms=%.3f "
            "warm_calls=%llu warm_ranges=%llu warm_bytes=%llu "
            "warm_ms=%.3f cache_hits=%llu cache_misses=%llu "
            "cache_evictions=%llu cache_allocations=%llu "
            "credit_wait_ms=%.3f failures=%llu\n",
            index_name ? index_name : "(unnamed)",
            component ? component : "(unknown)",
            stats.read_calls,
            stats.read_batches,
            stats.read_logical_bytes,
            stats.read_bytes,
            stats.read_pages,
            (double)stats.read_nanoseconds / 1000000.0,
            stats.warm_calls,
            stats.warm_ranges,
            stats.warm_bytes,
            (double)stats.warm_nanoseconds / 1000000.0,
            stats.cache_hits,
            stats.cache_misses,
            stats.cache_evictions,
            stats.cache_allocations,
            (double)stats.wait_nanoseconds / 1000000.0,
            stats.failures);
    }
    rc = fitsbin_close_payload_fd(fb);
    if (rc) {
        ERROR("Failed to close %s exact-payload descriptor",
              component ? component : "index");
    }
    return rc;
}

int index_close_payload_fds(index_t* ind) {
    int rc = 0;

    if (!ind) {
        return 0;
    }
    if (ind->quads &&
        index_close_one_payload_fd(
            ind->quads->fb,
            "quad",
            ind->indexname)) {
        rc = -1;
    }
    if (ind->codekd && ind->codekd->tree &&
        index_close_one_payload_fd(
            (fitsbin_t*)ind->codekd->tree->io,
            "codekd",
            ind->indexname)) {
        rc = -1;
    }
    if (ind->starkd && ind->starkd->tree &&
        index_close_one_payload_fd(
            (fitsbin_t*)ind->starkd->tree->io,
            "starkd",
            ind->indexname)) {
        rc = -1;
    }
    return rc;
}

void index_close(index_t* index) {
    if (!index) return;
    index_unload(index);
    if (index->fits)
        anqfits_close(index->fits);
    index->fits = NULL;
    if (index->residency)
        index_residency_release(index->residency);
    index->residency = NULL;
    if (index->residency_source)
        index_residency_release(index->residency_source);
    index->residency_source = NULL;
    free(index->indexname);
    free(index->indexfn);
    free(index->cutband);
    index->indexname = index->indexfn = NULL;
}

void index_free(index_t* index) {
    index_close(index);
    free(index);
}
