/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
#include <errno.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fitsbin.h"
#include "fitsioutils.h"

#include "cutest.h"

static char* get_tmpfile(int i) {
    static char fn[256];
    sprintf(fn, "/tmp/test-fitsbin-%i", i);
    return fn;
}

static fitsbin_t* make_mmap_test_fitsbin(void* mapping,
                                         size_t mapping_size) {
    fitsbin_t* fb;
    fitsbin_chunk_t chunk;

    fb = calloc(1, sizeof(*fb));
    if (!fb) {
        return NULL;
    }
    fb->chunks = bl_new(1, sizeof(fitsbin_chunk_t));
    if (!fb->chunks) {
        free(fb);
        return NULL;
    }
    fitsbin_chunk_init(&chunk);
    chunk.tablename = "mmap-test";
    chunk.map = mapping;
    chunk.mapsize = mapping_size;
    fitsbin_add_chunk(fb, &chunk);
    fitsbin_configure_index_mmap(fb);
    return fb;
}

// Verify that legacy mmap environment variables no longer alter production.
void test_fitsbin_mmap_advice(CuTest* ct) {
    static const char* envname = "ASTROMETRY_INDEX_MMAP_ADVICE";
    const char* oldvalue;
    char* saved = NULL;
    fitsbin_t fb;
    fitsbin_chunk_t payload;
    fitsbin_chunk_t topology;
    int result;
    int restore_result;

    oldvalue = getenv(envname);
    if (oldvalue) {
        saved = strdup(oldvalue);
        CuAssertPtrNotNull(ct, saved);
    }

    result = setenv(envname, "normal", 1);
    CuAssertIntEquals(ct, 0, result);

    /*
     * The parser is retained for controlled callers, but its fallback must
     * agree with the serial no-override NORMAL behavior. An absent or
     * misspelled value must not select RANDOM through the parser.
     */
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_NORMAL,
        fitsbin_mmap_policy_parse(NULL));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_NORMAL,
        fitsbin_mmap_policy_parse(""));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_NORMAL,
        fitsbin_mmap_policy_parse("unknown"));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_NORMAL,
        fitsbin_mmap_policy_parse("normal"));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_RANDOM,
        fitsbin_mmap_policy_parse("random"));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_ADAPTIVE,
        fitsbin_mmap_policy_parse("adaptive"));

    fitsbin_mmap_clear_thread_advice();
    memset(&fb, 0, sizeof(fb));

    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_POLICY_FIXED_RANDOM,
        fitsbin_get_configured_mmap_policy());

    result = fitsbin_configure_index_mmap(&fb);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        fb.mmap_advice);

#ifdef MADV_RANDOM
    fitsbin_mmap_set_thread_advice(FITSBIN_MMAP_ADVICE_RANDOM);
#else
    fitsbin_mmap_set_thread_advice(FITSBIN_MMAP_ADVICE_NORMAL);
#endif

    memset(&fb, 0, sizeof(fb));
    result = fitsbin_configure_index_mmap(&fb);
    CuAssertIntEquals(ct, 0, result);

#ifdef MADV_RANDOM
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        fb.mmap_advice);
#else
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        fb.mmap_advice);
#endif

    fitsbin_chunk_init(&payload);
    fitsbin_chunk_init(&topology);
    topology.mmap_region = FITSBIN_MMAP_REGION_TOPOLOGY;

    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_REGION_PAYLOAD,
        payload.mmap_region);
    CuAssertIntEquals(
        ct,
        fb.mmap_advice,
        fitsbin_get_mmap_advice(&fb));
    CuAssertIntEquals(
        ct,
        fb.mmap_advice,
        fitsbin_get_chunk_mmap_advice(&fb, &payload));
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        fitsbin_get_chunk_mmap_advice(&fb, &topology));

    fitsbin_mmap_clear_thread_advice();
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        fitsbin_mmap_current_advice());

    if (saved) {
        restore_result = setenv(envname, saved, 1);
    } else {
        restore_result = unsetenv(envname);
    }

    free(saved);
    CuAssertIntEquals(ct, 0, restore_result);
}

void test_fitsbin_mmap_prefetch(CuTest* ct) {
    static const char* envname = "ASTROMETRY_INDEX_MMAP_PREFETCH";
    const char* oldvalue;
    char* saved = NULL;
    fitsbin_t fb;
    int result;
    int restore_result;

    oldvalue = getenv(envname);
    if (oldvalue) {
        saved = strdup(oldvalue);
        CuAssertPtrNotNull(ct, saved);
    }

    result = setenv(envname, "1", 1);
    CuAssertIntEquals(ct, 0, result);

    memset(&fb, 0, sizeof(fb));
    result = fitsbin_configure_index_mmap(&fb);

    CuAssertIntEquals(ct, 0, result);
#if defined(MADV_POPULATE_READ)
    CuAssertIntEquals(ct, TRUE, fb.mmap_prefetch_enabled);
#else
    CuAssertIntEquals(ct, FALSE, fb.mmap_prefetch_enabled);
#endif

    if (saved) {
        restore_result = setenv(envname, saved, 1);
    } else {
        restore_result = unsetenv(envname);
    }

    free(saved);
    CuAssertIntEquals(ct, 0, restore_result);
}

void test_fitsbin_mmap_range_advice(CuTest* ct) {
    fitsbin_t empty_fb;
    fitsbin_t* fb;
    fitsbin_chunk_t* chunk;
    const void* map_base;
    const void* range_start;
    size_t map_size;
    size_t range_size;
    long detected_page_size;
    size_t page_size;
    char* mapping;
    int result;

    memset(&empty_fb, 0, sizeof(empty_fb));
    map_base = (const void*)1;
    map_size = 1U;
    range_start = (const void*)1;
    range_size = 1U;
    result = fitsbin_resolve_mapped_range(
        &empty_fb,
        &empty_fb,
        1U,
        &map_base,
        &map_size,
        &range_start,
        &range_size);
    CuAssertIntEquals(ct, 0, result);
    CuAssertPtrEquals(ct, NULL, map_base);
    CuAssertIntEquals(ct, 0, (int)map_size);
    CuAssertPtrEquals(ct, NULL, range_start);
    CuAssertIntEquals(ct, 0, (int)range_size);

    detected_page_size = sysconf(_SC_PAGESIZE);
    CuAssert(ct, "system page size unavailable", detected_page_size > 0);
    page_size = (size_t)detected_page_size;

    mapping = mmap(
        NULL,
        page_size * 2U,
        PROT_READ | PROT_WRITE,
        MAP_PRIVATE | MAP_ANONYMOUS,
        -1,
        0);
    CuAssert(ct, "anonymous mmap failed", mapping != MAP_FAILED);

    fb = make_mmap_test_fitsbin(mapping, page_size * 2U);
    CuAssertPtrNotNull(ct, fb);
    chunk = fitsbin_get_chunk(fb, 0);
    CuAssertPtrNotNull(ct, chunk);

    result = fitsbin_resolve_mapped_range(
        fb,
        mapping + 17,
        113U,
        &map_base,
        &map_size,
        &range_start,
        &range_size);
    CuAssertIntEquals(ct, 1, result);
    CuAssertPtrEquals(ct, mapping, map_base);
    CuAssertIntEquals(ct, (int)(page_size * 2U), (int)map_size);
    CuAssertPtrEquals(ct, mapping + 17, range_start);
    CuAssertIntEquals(ct, 113, (int)range_size);

    result = fitsbin_set_mmap_range_advice(
        fb,
        mapping + 17,
        113U,
        FITSBIN_MMAP_ADVICE_NORMAL);
    CuAssertIntEquals(ct, 1, result);

    result = fitsbin_resolve_mapped_range(
        fb,
        mapping + page_size * 2U - 32U,
        64U,
        &map_base,
        &map_size,
        &range_start,
        &range_size);
    CuAssertIntEquals(ct, 1, result);
    CuAssertIntEquals(ct, 32, (int)range_size);

    errno = 0;
    result = fitsbin_set_mmap_range_advice(
        fb,
        mapping + page_size * 2U - 32U,
        64U,
        FITSBIN_MMAP_ADVICE_NORMAL);
    CuAssertIntEquals(ct, -1, result);
    CuAssertIntEquals(ct, ERANGE, errno);

    /*
     * Force the actual madvise() call to reject an unaligned mapping. A
     * range-local failure must not poison the whole-file policy state.
     */
    chunk->map = mapping + 1;
    chunk->mapsize = page_size * 2U - 1U;
    errno = 0;
    result = fitsbin_set_mmap_range_advice(
        fb,
        mapping + 2,
        1U,
        FITSBIN_MMAP_ADVICE_NORMAL);
    CuAssertIntEquals(ct, -1, result);
    CuAssertTrue(ct, errno != 0);
    CuAssertIntEquals(ct, FALSE, fb->mmap_advice_failed);

    chunk->map = mapping;
    chunk->mapsize = page_size * 2U;
    CuAssertIntEquals(ct, 0, fitsbin_close(fb));
}

void test_fitsbin_mmap_advice_retry(CuTest* ct) {
    fitsbin_t empty_fb;
    fitsbin_t* fb;
    fitsbin_chunk_t* chunk;
    long detected_page_size;
    size_t page_size;
    char* mapping;
    int result;

    memset(&empty_fb, 0, sizeof(empty_fb));
    empty_fb.mmap_advice = FITSBIN_MMAP_ADVICE_NORMAL;
    empty_fb.mmap_advice_failed = TRUE;
    result = fitsbin_set_mmap_advice(
        &empty_fb,
        FITSBIN_MMAP_ADVICE_NORMAL,
        TRUE);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(ct, FALSE, empty_fb.mmap_advice_failed);

    detected_page_size = sysconf(_SC_PAGESIZE);
    CuAssert(ct, "system page size unavailable", detected_page_size > 0);
    page_size = (size_t)detected_page_size;

    mapping = mmap(
        NULL,
        page_size,
        PROT_READ | PROT_WRITE,
        MAP_PRIVATE | MAP_ANONYMOUS,
        -1,
        0);
    CuAssert(ct, "anonymous mmap failed", mapping != MAP_FAILED);

    /*
     * Keep the real allocation alive for cleanup, but advertise an unaligned
     * VMA address so the first full-chunk madvise() deterministically fails.
     */
    fb = make_mmap_test_fitsbin(mapping + 1, page_size - 1U);
    CuAssertPtrNotNull(ct, fb);
    chunk = fitsbin_get_chunk(fb, 0);
    CuAssertPtrNotNull(ct, chunk);

    errno = 0;
    result = fitsbin_set_mmap_advice(
        fb,
        FITSBIN_MMAP_ADVICE_RANDOM,
        TRUE);
    CuAssertIntEquals(ct, -1, result);
    CuAssertTrue(ct, errno != 0);
    CuAssertIntEquals(ct, TRUE, fb->mmap_advice_failed);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        fb->mmap_advice);

    /*
     * The fallback records NORMAL but cannot certify a mapping when madvise()
     * itself failed. Preserve the failure until a full reapplication works.
     */
    result = fitsbin_set_mmap_advice(
        fb,
        FITSBIN_MMAP_ADVICE_NORMAL,
        FALSE);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(ct, TRUE, fb->mmap_advice_failed);

    chunk->map = mapping;
    chunk->mapsize = page_size;
    result = fitsbin_set_mmap_advice(
        fb,
        FITSBIN_MMAP_ADVICE_NORMAL,
        TRUE);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(ct, FALSE, fb->mmap_advice_failed);

    /*
     * A successful same-advice call is idempotent and remains healthy.
     */
    result = fitsbin_set_mmap_advice(
        fb,
        FITSBIN_MMAP_ADVICE_NORMAL,
        TRUE);
    CuAssertIntEquals(ct, 0, result);
    CuAssertIntEquals(ct, FALSE, fb->mmap_advice_failed);
    CuAssertIntEquals(ct, 0, fitsbin_close(fb));
}

void test_fitsbin_open_file_identity(CuTest* ct) {
    fitsbin_t* in;
    fitsbin_t* out;
    struct stat before_close;
    struct stat after_close;
    char* fn;

    fn = get_tmpfile(99);
    out = fitsbin_open_for_writing(fn);
    CuAssertPtrNotNull(ct, out);
    CuAssertIntEquals(ct, 0, fitsbin_write_primary_header(out));
    CuAssertIntEquals(ct, 0, fitsbin_fix_primary_header(out));
    CuAssertIntEquals(ct, 0, fitsbin_close(out));

    in = fitsbin_open(fn);
    CuAssertPtrNotNull(ct, in);
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_get_open_file_stat(in, &before_close));
    CuAssertIntEquals(ct, 0, fitsbin_close_fd(in));
    CuAssertPtrEquals(ct, NULL, fitsbin_get_fid(in));
    CuAssertIntEquals(
        ct,
        0,
        fitsbin_get_open_file_stat(in, &after_close));
    CuAssertTrue(ct, before_close.st_dev == after_close.st_dev);
    CuAssertTrue(ct, before_close.st_ino == after_close.st_ino);
    CuAssertTrue(ct, before_close.st_size == after_close.st_size);
    CuAssertTrue(ct, before_close.st_mtime == after_close.st_mtime);
    CuAssertTrue(ct, before_close.st_ctime == after_close.st_ctime);
    CuAssertIntEquals(ct, 0, fitsbin_close(in));
    CuAssertIntEquals(ct, 0, unlink(fn));
}

void test_fitsbin_mmap_adaptive_policy(CuTest* ct) {
    fitsbin_mmap_advice_state_t state;
    anbool transitioned;

    fitsbin_mmap_advice_state_init(
        &state,
        FITSBIN_MMAP_POLICY_ADAPTIVE);

    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    CuAssertIntEquals(ct, 0, state.pass_number);
    CuAssertIntEquals(
        ct,
        0,
        state.completed_clean_unsolved_passes);

    CuAssertIntEquals(ct, 0, state.transition_count);

    /*
     * Partial pass.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    /*
     * Cancelled pass.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        TRUE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    /*
     * Failed pass.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        -1,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    /*
     * Solved pass.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        FALSE,
        TRUE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    /*
     * First clean, exhaustive, unsolved pass.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, TRUE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        state.effective_advice);

    CuAssertIntEquals(
        ct,
        1,
        state.completed_clean_unsolved_passes);

    CuAssertIntEquals(ct, 1, state.transition_count);

    /*
     * Later clean pass remains NORMAL and does not transition again.
     */
    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        state.effective_advice);

    CuAssertIntEquals(ct, 1, state.transition_count);

    /*
     * New field.
     */
    fitsbin_mmap_advice_state_reset(&state);

    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);

    CuAssertIntEquals(ct, 0, state.pass_number);
    CuAssertIntEquals(
        ct,
        0,
        state.completed_clean_unsolved_passes);

    CuAssertIntEquals(ct, 0, state.transition_count);

    /*
     * Fixed policies never transition.
     */
    fitsbin_mmap_advice_state_init(
        &state,
        FITSBIN_MMAP_POLICY_FIXED_NORMAL);

    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_NORMAL,
        state.effective_advice);

    fitsbin_mmap_advice_state_init(
        &state,
        FITSBIN_MMAP_POLICY_FIXED_RANDOM);

    transitioned = fitsbin_mmap_policy_complete_pass(
        &state,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        0,
        0);

    CuAssertIntEquals(ct, FALSE, transitioned);
    CuAssertIntEquals(
        ct,
        FITSBIN_MMAP_ADVICE_RANDOM,
        state.effective_advice);
}

void test_fitsbin_1(CuTest* ct) {
    fitsbin_t* in, *out;
    int i;
    int N = 6;
    double outdata[6];
    double* indata;
    char* fn;
    fitsbin_chunk_t chunk;

    fn = get_tmpfile(0);
    out = fitsbin_open_for_writing(fn);
    CuAssertPtrNotNull(ct, out);

    CuAssertIntEquals(ct, 0, fitsbin_write_primary_header(out));

    for (i=0; i<N; i++) {
        outdata[i] = i*i;
    }

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test1";
    chunk.itemsize = sizeof(double);
    chunk.nrows = N;
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk(out, &chunk));
    CuAssertIntEquals(ct, fitsbin_fix_primary_header(out), 0);
    CuAssertIntEquals(ct, fitsbin_close(out), 0);

    fitsbin_chunk_clean(&chunk);

    // writing shouldn't affect the data values
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, outdata[i], i*i);
    }

    in = fitsbin_open(fn);
    CuAssertPtrNotNull(ct, in);

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test1";

    CuAssertIntEquals(ct, 0, fitsbin_read_chunk(in, &chunk));
    CuAssertIntEquals(ct, sizeof(double), chunk.itemsize);
    CuAssertIntEquals(ct, N, chunk.nrows);
    indata = chunk.data;
    CuAssertPtrNotNull(ct, indata);
    CuAssertIntEquals(ct, 0, memcmp(outdata, indata, sizeof(outdata)));
    CuAssertIntEquals(ct, 0, fitsbin_close(in));
}

void test_fitsbin_2(CuTest* ct) {
    fitsbin_t* in, *out;
    int i;
    int N = 6;
    double outdata[6];
    double* indata;
    char* fn;
    fitsbin_chunk_t chunk;
    fitsbin_chunk_t* ch;

    fn = get_tmpfile(0);
    printf("Writing to %s\n", fn);
    out = fitsbin_open_for_writing(fn);
    CuAssertPtrNotNull(ct, out);

    CuAssertIntEquals(ct, 0, fitsbin_write_primary_header(out));

    for (i=0; i<N; i++) {
        outdata[i] = i*i;
    }

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test2";
    chunk.itemsize = 1;
    //chunk.nrows = N * sizeof(double);
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk_header(out, &chunk));
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, 0, fitsbin_write_items(out, &chunk, &outdata[i], sizeof(double)));
    }
    CuAssertIntEquals(ct, 0, fitsbin_fix_chunk_header(out, &chunk));

    fitsbin_chunk_reset(&chunk);
    chunk.tablename = "test2B";
    chunk.itemsize = sizeof(double);
    //chunk.nrows = N;
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk_header(out, &chunk));
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, 0, fitsbin_write_items(out, &chunk, outdata + (N-1) - i, 1));
    }
    CuAssertIntEquals(ct, 0, fitsbin_fix_chunk_header(out, &chunk));

    CuAssertIntEquals(ct, fitsbin_close(out), 0);

    fitsbin_chunk_clean(&chunk);

    // writing shouldn't affect the data values
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, outdata[i], i*i);
    }

    in = fitsbin_open(fn);
    CuAssertPtrNotNull(ct, in);

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test2";
    fitsbin_add_chunk(in, &chunk);
    chunk.tablename = "test2B";
    fitsbin_add_chunk(in, &chunk);

    CuAssertIntEquals(ct, 0, fitsbin_read(in));

    ch = fitsbin_get_chunk(in, 0);
    CuAssertIntEquals(ct, 1, ch->itemsize);
    CuAssertIntEquals(ct, N * sizeof(double), ch->nrows);
    indata = ch->data;
    CuAssertPtrNotNull(ct, indata);
    CuAssertIntEquals(ct, 0, memcmp(outdata, indata, sizeof(outdata)));

    ch = fitsbin_get_chunk(in, 1);
    CuAssertIntEquals(ct, sizeof(double), ch->itemsize);
    CuAssertIntEquals(ct, N, ch->nrows);
    indata = ch->data;
    CuAssertPtrNotNull(ct, indata);
    for (i=0; i<N; i++) {
        CuAssertDblEquals(ct, outdata[N-1 - i], indata[i], 1e-10);
    }

    CuAssertIntEquals(ct, 0, fitsbin_close(in));
}



void test_inmemory_fitsbin_1(CuTest* ct) {
    fitsbin_t* fb;
    int i;
    int N = 6;
    double outdata[6];
    double* indata;
    fitsbin_chunk_t chunk;

    fb = fitsbin_open_in_memory();
    CuAssertPtrNotNull(ct, fb);

    CuAssertIntEquals(ct, 0, fitsbin_write_primary_header(fb));

    for (i=0; i<N; i++) {
        outdata[i] = i*i;
    }

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test1";
    chunk.itemsize = sizeof(double);
    chunk.nrows = N;
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk(fb, &chunk));
    CuAssertIntEquals(ct, 0, fitsbin_fix_primary_header(fb));

    fitsbin_chunk_clean(&chunk);

    // writing shouldn't affect the data values
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, outdata[i], i*i);
    }

    CuAssertIntEquals(ct, 0, fitsbin_switch_to_reading(fb));

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test1";

    CuAssertIntEquals(ct, 0, fitsbin_read_chunk(fb, &chunk));
    CuAssertIntEquals(ct, sizeof(double), chunk.itemsize);
    CuAssertIntEquals(ct, N, chunk.nrows);
    indata = chunk.data;
    CuAssertPtrNotNull(ct, indata);
    CuAssertIntEquals(ct, 0, memcmp(outdata, indata, sizeof(outdata)));

    CuAssertIntEquals(ct, 0, fitsbin_close(fb));
}

void test_inmemory_fitsbin_2(CuTest* ct) {
    fitsbin_t* fb;
    int i;
    int N = 6;
    double outdata[6];
    double* indata;
    fitsbin_chunk_t chunk;
    fitsbin_chunk_t* ch;

    fb = fitsbin_open_in_memory();
    CuAssertPtrNotNull(ct, fb);

    CuAssertIntEquals(ct, 0, fitsbin_write_primary_header(fb));

    for (i=0; i<N; i++) {
        outdata[i] = i*i;
    }

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test2";
    chunk.itemsize = 1;
    //chunk.nrows = N * sizeof(double);
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk_header(fb, &chunk));
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, 0, fitsbin_write_items(fb, &chunk, &outdata[i], sizeof(double)));
    }
    CuAssertIntEquals(ct, 0, fitsbin_fix_chunk_header(fb, &chunk));

    fitsbin_chunk_reset(&chunk);
    chunk.tablename = "test2B";
    chunk.itemsize = sizeof(double);
    //chunk.nrows = N;
    chunk.data = outdata;

    CuAssertIntEquals(ct, 0, fitsbin_write_chunk_header(fb, &chunk));
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, 0, fitsbin_write_items(fb, &chunk, outdata + (N-1) - i, 1));
    }
    CuAssertIntEquals(ct, 0, fitsbin_fix_chunk_header(fb, &chunk));

    fitsbin_chunk_clean(&chunk);

    // writing shouldn't affect the data values
    for (i=0; i<N; i++) {
        CuAssertIntEquals(ct, outdata[i], i*i);
    }

    CuAssertIntEquals(ct, 0, fitsbin_switch_to_reading(fb));

    fitsbin_chunk_init(&chunk);
    chunk.tablename = "test2";
    fitsbin_add_chunk(fb, &chunk);
    chunk.tablename = "test2B";
    fitsbin_add_chunk(fb, &chunk);

    CuAssertIntEquals(ct, 0, fitsbin_read(fb));

    ch = fitsbin_get_chunk(fb, 0);
    CuAssertIntEquals(ct, 1, ch->itemsize);
    CuAssertIntEquals(ct, N * sizeof(double), ch->nrows);
    indata = ch->data;
    CuAssertPtrNotNull(ct, indata);
    CuAssertIntEquals(ct, 0, memcmp(outdata, indata, sizeof(outdata)));

    ch = fitsbin_get_chunk(fb, 1);
    CuAssertIntEquals(ct, sizeof(double), ch->itemsize);
    CuAssertIntEquals(ct, N, ch->nrows);
    indata = ch->data;
    CuAssertPtrNotNull(ct, indata);
    for (i=0; i<N; i++) {
        CuAssertDblEquals(ct, outdata[N-1 - i], indata[i], 1e-10);
    }

    CuAssertIntEquals(ct, 0, fitsbin_close(fb));
}
