/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "mathutil.h"
#include "permutedsort.h"
#include "sip-utils.h"
#include "verify_prepared_internal.h"

int verify_query_hit(const startree_t* skdt,
                     const double center[3],
                     double radius2,
                     verify_index_query_t** query) {
    verify_index_query_t* result;

    if (!query) {
        return -1;
    }
    *query = NULL;
    if (!skdt || !center || !skdt->tree || !skdt->sweep) {
        return -1;
    }
    result = calloc(1, sizeof(*result));
    if (!result) {
        return -1;
    }
    result->source = skdt;
    memcpy(result->center, center, sizeof(result->center));
    result->radius2 = radius2;
    startree_search_for(skdt, center, radius2,
                        &result->refxyz, NULL,
                        &result->refstarid, &result->nrall);
    if (result->nrall < 0 ||
        (result->nrall &&
         (!result->refxyz || !result->refstarid)) ||
        (!result->nrall &&
         (result->refxyz || result->refstarid))) {
        verify_destroy_index_query(result);
        return -1;
    }
    *query = result;
    return 0;
}

size_t verify_index_query_count(const verify_index_query_t* query) {
    if (!query || query->nrall <= 0) {
        return 0U;
    }
    return (size_t)query->nrall;
}

size_t verify_index_query_bytes(const verify_index_query_t* query) {
    size_t count;
    size_t total;

    if (!query || query->nrall < 0) {
        return 0U;
    }
    count = (size_t)query->nrall;
    total = sizeof(*query);
    if (count > (SIZE_MAX - total) / (3U * sizeof(double))) {
        return SIZE_MAX;
    }
    total += count * 3U * sizeof(double);
    if (count > (SIZE_MAX - total) / sizeof(int)) {
        return SIZE_MAX;
    }
    total += count * sizeof(int);
    if (count > (SIZE_MAX - total) / sizeof(*query->sweep)) {
        return SIZE_MAX;
    }
    total += count * sizeof(*query->sweep);
    return total;
}

int verify_index_query_sweep_range(const startree_t* skdt,
                                   const verify_index_query_t* query,
                                   size_t index,
                                   const void** data,
                                   size_t* size) {
    int starid;

    if (!data || !size) {
        return -1;
    }
    *data = NULL;
    *size = 0U;
    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        index >= (size_t)query->nrall ||
        !query->refstarid) {
        return -1;
    }
    starid = query->refstarid[index];
    if (starid < 0 || starid >= startree_N(skdt)) {
        return -1;
    }
    *data = skdt->sweep + starid;
    *size = sizeof(*skdt->sweep);
    return 0;
}

int verify_index_query_capture_sweep(const startree_t* skdt,
                                     verify_index_query_t* query) {
    uint8_t* sweep;
    size_t count;
    int nstars;
    int i;

    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        (query->nrall && !query->refstarid)) {
        return -1;
    }
    if (query->sweep || !query->nrall) {
        return 0;
    }
    count = (size_t)query->nrall;
    nstars = startree_N(skdt);
    for (i = 0; i < query->nrall; i++) {
        if (query->refstarid[i] < 0 ||
            query->refstarid[i] >= nstars) {
            return -1;
        }
    }
    sweep = malloc(count * sizeof(*sweep));
    if (!sweep) {
        return -1;
    }
    for (i = 0; i < query->nrall; i++) {
        if (!(i & 255) &&
            verify_internal_worker_stop_requested()) {
            free(sweep);
            return -1;
        }
        sweep[i] = skdt->sweep[query->refstarid[i]];
    }
    query->sweep = sweep;
    return 0;
}

static int verify_mapped_page_buffers_validate(
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count) {
    uintptr_t previous_end = 0U;
    size_t i;

    if (!buffers || !buffer_count) {
        return -1;
    }
    for (i = 0U; i < buffer_count; i++) {
        uintptr_t begin;
        uintptr_t end;

        if (!buffers[i].mapping_data || !buffers[i].bytes ||
            !buffers[i].size) {
            return -1;
        }
        begin = (uintptr_t)buffers[i].mapping_data;
        if (buffers[i].size > UINTPTR_MAX - begin) {
            return -1;
        }
        end = begin + buffers[i].size;
        if (i && begin < previous_end) {
            return -1;
        }
        previous_end = end;
    }
    return 0;
}

static int verify_mapped_page_buffer_find(
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count,
    uintptr_t target,
    size_t target_size,
    const unsigned char** bytes) {
    size_t low = 0U;
    size_t high = buffer_count;

    if (!bytes || !target_size) {
        return -1;
    }
    *bytes = NULL;
    while (low < high) {
        size_t middle = low + (high - low) / 2U;
        uintptr_t begin = (uintptr_t)buffers[middle].mapping_data;
        uintptr_t end = begin + buffers[middle].size;

        if (target < begin) {
            high = middle;
            continue;
        }
        if (target >= end) {
            low = middle + 1U;
            continue;
        }
        if (target_size > end - target) {
            return -1;
        }
        *bytes = buffers[middle].bytes + (size_t)(target - begin);
        return 0;
    }
    return -1;
}

int verify_index_query_capture_sweep_buffers(
    const startree_t* skdt,
    verify_index_query_t* query,
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count) {
    uint8_t* sweep;
    size_t count;
    int nstars;
    int i;

    if (!skdt || !query || query->source != skdt ||
        !skdt->tree || !skdt->sweep ||
        query->nrall < 0 ||
        (query->nrall && !query->refstarid)) {
        return -1;
    }
    if (query->sweep || !query->nrall) {
        return 0;
    }
    if (verify_mapped_page_buffers_validate(
            buffers, buffer_count)) {
        return -1;
    }
    count = (size_t)query->nrall;
    if (count > SIZE_MAX / sizeof(*sweep)) {
        return -1;
    }
    sweep = malloc(count * sizeof(*sweep));
    if (!sweep) {
        return -1;
    }
    nstars = startree_N(skdt);
    for (i = 0; i < query->nrall; i++) {
        const unsigned char* bytes;
        uintptr_t target;
        int starid = query->refstarid[i];

        if (!(i & 255) &&
            verify_internal_worker_stop_requested()) {
            free(sweep);
            return -1;
        }
        if (starid < 0 || starid >= nstars) {
            free(sweep);
            return -1;
        }
        target = (uintptr_t)(skdt->sweep + starid);
        if (verify_mapped_page_buffer_find(
                buffers, buffer_count,
                target, sizeof(*sweep), &bytes)) {
            free(sweep);
            return -1;
        }
        memcpy(&sweep[i], bytes, sizeof(*sweep));
    }
    query->sweep = sweep;
    return 0;
}

void verify_destroy_index_query(verify_index_query_t* query) {
    if (!query) {
        return;
    }
    free(query->refxyz);
    free(query->refstarid);
    free(query->sweep);
    memset(query, 0, sizeof(*query));
    free(query);
}


int verify_prepare_hit_from_query(const startree_t* skdt,
                                  verify_index_query_t** query,
                                  int index_cutnside,
                                  const MatchObj* mo, const sip_t* sip,
                                  const verify_field_t* vf,
                                  double pix2, double distractors,
                                  double fieldW, double fieldH,
                                  double logbail, double logaccept,
                                  double logstoplooking,
                                  anbool do_gamma, anbool fake_match,
                                  verify_prepared_hit_t** prepared) {
    verify_index_query_t* query_context;
    verify_prepared_hit_t* context;
    verify_t* v;
    double fieldr2;
    int* sweep = NULL;
    int nstars;
    int i;
    int j;
    int ibad;
    int igood;

    if (!prepared) {
        return -1;
    }
    *prepared = NULL;
    if (!query || !*query || !skdt || !mo || !vf ||
        (!mo->wcs_valid && !sip) ||
        !isfinite(logaccept) || !isfinite(logbail)) {
        return -1;
    }
    query_context = *query;
    fieldr2 = square(mo->radius);
    if (!skdt->tree ||
        query_context->source != skdt ||
        memcmp(query_context->center, mo->center,
               sizeof(query_context->center)) ||
        memcmp(&query_context->radius2, &fieldr2, sizeof(fieldr2)) ||
        query_context->nrall < 0 ||
        (query_context->nrall &&
         (!query_context->refxyz || !query_context->refstarid)) ||
        (!query_context->nrall &&
         (query_context->refxyz || query_context->refstarid)) ||
        (query_context->nrall && !query_context->sweep &&
         !skdt->sweep)) {
        return -1;
    }
    context = calloc(1, sizeof(*context));
    if (!context) {
        return -1;
    }
    context->distractors = distractors;
    context->logbail = logbail;
    context->logaccept = logaccept;
    context->logstoplooking = logstoplooking;
    context->fake_match = fake_match;
    context->state = VERIFY_PREPARED_READY;
    v = &context->verify;

    if (sip) {
        memcpy(&context->wcs, sip, sizeof(context->wcs));
    } else {
        sip_wrap_tan(&mo->wcstan, &context->wcs);
    }
    v->wcs = &context->wcs;
    context->refxyz = query_context->refxyz;
    v->refstarid = query_context->refstarid;
    v->NRall = query_context->nrall;
    if (!context->refxyz) {
        context->state = VERIFY_PREPARED_NO_REFERENCE;
        goto done;
    }

    if ((size_t)v->NRall > SIZE_MAX / (2U * sizeof(double)) ||
        (size_t)v->NRall > SIZE_MAX / sizeof(int)) {
        goto fail;
    }
    v->refxy = malloc((size_t)v->NRall * 2U * sizeof(double));
    v->refperm = malloc((size_t)v->NRall * sizeof(int));
    if (!v->refxy || !v->refperm) {
        goto fail;
    }
    igood = 0;
    for (i = 0; i < v->NRall; i++) {
        if (!sip_xyzarr2pixelxy(v->wcs,
                                context->refxyz + 3 * i,
                                v->refxy + 2 * i,
                                v->refxy + 2 * i + 1) ||
            !sip_pixel_is_inside_image(v->wcs,
                                       v->refxy[2 * i],
                                       v->refxy[2 * i + 1])) {
            continue;
        }
        v->refperm[igood++] = i;
    }
    v->NR = igood;
    context->nrimage = v->NR;

    sweep = malloc((size_t)v->NRall * sizeof(int));
    if (!sweep) {
        goto fail;
    }
    nstars = startree_N(skdt);
    for (i = 0; i < v->NRall; i++) {
        int starid = v->refstarid[i];

        if (starid < 0 || starid >= nstars) {
            goto fail;
        }
        if (query_context->sweep) {
            sweep[i] = query_context->sweep[i];
        } else {
            sweep[i] = skdt->sweep[starid];
        }
    }
    permuted_sort(sweep, sizeof(int), compare_ints_asc,
                  v->refperm, v->NR);
    free(sweep);
    sweep = NULL;

    if (v->NR) {
        v->badguys = malloc((size_t)v->NR * sizeof(int));
        if (!v->badguys) {
            goto fail;
        }
    }
    if (!fake_match) {
        ibad = 0;
        igood = 0;
        for (i = 0; i < v->NR; i++) {
            anbool inquad = FALSE;
            int ri = v->refperm[i];

            for (j = 0; j < mo->dimquads; j++) {
                if (v->refstarid[ri] == (int)mo->star[j]) {
                    inquad = TRUE;
                    v->badguys[ibad++] = ri;
                    break;
                }
            }
            if (!inquad) {
                v->refperm[igood++] = ri;
            }
        }
        if (ibad) {
            memcpy(v->refperm + igood, v->badguys,
                   (size_t)ibad * sizeof(int));
        }
        v->NR = igood;
    }
    if (!v->NR) {
        context->state = VERIFY_PREPARED_NO_QUAD_REFERENCE;
        goto done;
    }

    if (!fake_match) {
        if (verify_internal_apply_ror(
                v, index_cutnside, (MatchObj*)mo,
                vf, pix2, distractors, fieldW, fieldH,
                do_gamma, fake_match,
                &context->effective_area, NULL, NULL)) {
            goto fail;
        }
        if (!v->NR) {
            context->state = VERIFY_PREPARED_NO_ROR_REFERENCE;
            goto done;
        }
    } else {
        if (verify_internal_get_test_stars(
                v, vf, (MatchObj*)mo,
                pix2, do_gamma, fake_match)) {
            goto fail;
        }
        context->effective_area = fieldW * fieldH;
    }
    if (!v->NR || !v->NT) {
        context->state = VERIFY_PREPARED_EMPTY_LISTS;
        goto done;
    }

done:
    free(v->badguys);
    v->badguys = NULL;
    free(v->tbadguys);
    v->tbadguys = NULL;
    query_context->refxyz = NULL;
    query_context->refstarid = NULL;
    verify_destroy_index_query(query_context);
    *query = NULL;
    *prepared = context;
    return 0;

fail:
    free(sweep);
    context->refxyz = NULL;
    v->refstarid = NULL;
    verify_destroy_prepared_hit(context);
    return -1;
}

int verify_prepare_hit(const startree_t* skdt, int index_cutnside,
                       const MatchObj* mo, const sip_t* sip,
                       const verify_field_t* vf,
                       double pix2, double distractors,
                       double fieldW, double fieldH,
                       double logbail, double logaccept,
                       double logstoplooking,
                       anbool do_gamma, anbool fake_match,
                       verify_prepared_hit_t** prepared) {
    verify_index_query_t* query = NULL;
    double fieldr2;
    int status;

    if (!prepared || !skdt || !mo || !vf ||
        (!mo->wcs_valid && !sip) ||
        !isfinite(logaccept) || !isfinite(logbail)) {
        return -1;
    }
    *prepared = NULL;
    fieldr2 = square(mo->radius);
    if (verify_query_hit(skdt, mo->center, fieldr2, &query)) {
        return -1;
    }
    status = verify_prepare_hit_from_query(
        skdt, &query, index_cutnside, mo, sip, vf,
        pix2, distractors, fieldW, fieldH,
        logbail, logaccept, logstoplooking,
        do_gamma, fake_match, prepared);
    verify_destroy_index_query(query);
    return status;
}

int verify_score_prepared_hit(const verify_prepared_hit_t* prepared,
                              verify_prepared_score_t* score) {
    verify_t local;
    anbool score_completed;

    if (!prepared || !score || score->theta || score->allodds ||
        score->complete) {
        return -1;
    }
    memset(score, 0, sizeof(*score));
    score->logodds = -LARGE_VAL;
    score->worstlogodds = -LARGE_VAL;
    score->besti = -1;
    score->ibailed = -1;
    score->istopped = -1;
    if (prepared->state != VERIFY_PREPARED_READY) {
        score->complete = TRUE;
        return 0;
    }

    local = prepared->verify;
    if (local.NR <= 0 ||
        (size_t)local.NR > SIZE_MAX / sizeof(int)) {
        return -1;
    }
    local.badguys = malloc((size_t)local.NR * sizeof(int));
    if (!local.badguys) {
        return -1;
    }
    local.tbadguys = NULL;
    score->logodds = verify_internal_star_lists(
        &local,
        prepared->effective_area,
        prepared->distractors,
        prepared->logbail,
        prepared->logstoplooking,
        &score->besti,
        &score->allodds,
        &score->theta,
        &score->worstlogodds,
        &score->ibailed,
        &score->istopped,
        &score_completed);
    free(local.badguys);
    if (!score_completed) {
        verify_destroy_prepared_score(score);
        return -1;
    }
    score->complete = TRUE;
    return 0;
}

int verify_finish_prepared_hit(verify_prepared_hit_t* prepared,
                               verify_prepared_score_t* score,
                               MatchObj* mo) {
    verify_t* v;
    double* refxyz;
    double K;
    int besti;
    int i;
    int j;

    if (!prepared || !score || !score->complete || !mo) {
        return -1;
    }
    switch (prepared->state) {
    case VERIFY_PREPARED_NO_REFERENCE:
        logverb("No reference stars in the bounding circle\n");
        verify_internal_set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_NO_QUAD_REFERENCE:
        logverb("After removing quad stars: no reference stars\n");
        verify_internal_set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_NO_ROR_REFERENCE:
        logerr("After applying ROR, NR = 0!\n");
        verify_internal_set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_EMPTY_LISTS:
        logverb("After applying RoR, NR=%i, NT=%i\n",
                prepared->verify.NR, prepared->verify.NT);
        verify_internal_set_null_mo(mo);
        verify_destroy_prepared_score(score);
        return 0;
    case VERIFY_PREPARED_READY:
        break;
    default:
        return -1;
    }

    v = &prepared->verify;
    refxyz = prepared->refxyz;
    K = score->logodds;
    besti = score->besti;

    if (log_get_level() >= LOG_ALL) {
        int nm;
        int nc;
        int nd;

        verify_count_hits(score->theta, besti, &nm, &nc, &nd);
        debug("verify: logodds %g, %i matches, %i conflicts, "
              "%i distractors after %i field objects.\n",
              K, nm, nc, nd, besti);
    }

    if (K >= prepared->logaccept) {
        int* etheta;
        double* eodds;
        int nm;
        int nc;
        int nd;

        verify_count_hits(score->theta, besti, &nm, &nc, &nd);
        if (verify_internal_fixup_theta(
                score->theta, score->allodds,
                score->ibailed, score->istopped,
                v, besti, prepared->nrimage, refxyz,
                &etheta, &eodds)) {
            return -1;
        }
        mo->logodds = K;
        mo->worstlogodds = score->worstlogodds;
        mo->nfield = v->NTall;
        mo->nindex = prepared->nrimage;
        mo->nmatch = nm;
        mo->nconflict = nc;
        mo->ndistractor = nd;

        if (!prepared->fake_match) {
            for (j = 0; j < mo->dimquads; j++) {
                for (i = 0; i < prepared->nrimage; i++) {
                    if (v->refstarid[i] == (int)mo->star[j]) {
                        int ti = mo->field[j];

                        assert(etheta[ti] == THETA_FILTERED);
                        etheta[ti] = i;
                        eodds[ti] = LARGE_VAL;
                        break;
                    }
                }
            }
        }

        mo->theta = etheta;
        mo->matchodds = eodds;
        mo->refxyz = prepared->refxyz;
        prepared->refxyz = NULL;
        mo->refxy = v->refxy;
        v->refxy = NULL;
        mo->refstarid = v->refstarid;
        v->refstarid = NULL;
        mo->testperm = v->testperm;
        v->testperm = NULL;
        matchobj_compute_derived(mo);
    } else {
        mo->logodds = K;
        mo->worstlogodds = score->worstlogodds;
        mo->nfield = v->NTall;
        mo->nindex = prepared->nrimage;
    }
    verify_destroy_prepared_score(score);
    return 0;
}

static size_t verify_prepared_add_bytes(size_t total,
                                        size_t count,
                                        size_t element_size) {
    size_t bytes;

    if (count > SIZE_MAX / element_size) {
        return SIZE_MAX;
    }
    bytes = count * element_size;
    if (total > SIZE_MAX - bytes) {
        return SIZE_MAX;
    }
    return total + bytes;
}

size_t verify_prepared_hit_bytes(const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t total;

    if (!prepared) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NRall < 0 || v->NTall < 0) {
        return SIZE_MAX;
    }
    total = sizeof(*prepared);
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, 3U * sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, 2U * sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, sizeof(int));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NRall, sizeof(int));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NTall, sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NTall, sizeof(int));
    return total;
}

size_t verify_prepared_score_bytes(
    const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t total = 0U;

    if (!prepared) {
        return 0U;
    }
    if (prepared->state != VERIFY_PREPARED_READY) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NT < 0) {
        return SIZE_MAX;
    }
    total = verify_prepared_add_bytes(
        total, (size_t)v->NT, sizeof(double));
    total = verify_prepared_add_bytes(
        total, (size_t)v->NT, sizeof(int));
    return total;
}

size_t verify_prepared_hit_peak_bytes(
    const verify_prepared_hit_t* prepared) {
    const verify_t* v;
    size_t finish_peak = 0U;
    size_t nn_bytes;
    size_t retained;
    size_t score_peak = 0U;
    size_t transient_peak;

    if (!prepared) {
        return 0U;
    }
    v = &prepared->verify;
    if (v->NR < 0 || v->NRall < 0 ||
        v->NT < 0 || v->NTall < 0 ||
        prepared->nrimage < 0) {
        return SIZE_MAX;
    }
    retained = verify_prepared_hit_bytes(prepared);
    if (retained == SIZE_MAX) {
        return SIZE_MAX;
    }
    switch (prepared->state) {
    case VERIFY_PREPARED_NO_REFERENCE:
    case VERIFY_PREPARED_NO_QUAD_REFERENCE:
    case VERIFY_PREPARED_NO_ROR_REFERENCE:
    case VERIFY_PREPARED_EMPTY_LISTS:
        return retained;
    case VERIFY_PREPARED_READY:
        break;
    default:
        return SIZE_MAX;
    }
    nn_bytes = verify_internal_score_workspace_bytes(v->NR);
    if (nn_bytes == SIZE_MAX) {
        return SIZE_MAX;
    }

    /*
     * Scoring retains one conflict array, the packed nearest-neighbor
     * inputs, the result vectors, and at most one nearest-neighbor
     * workspace. Grid and legacy KD storage never coexist, so use the
     * larger exact payload bound.
     */
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NR, sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NR,
        3U * sizeof(double) + sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NT,
        sizeof(double) + sizeof(unsigned char));
    score_peak = verify_prepared_add_bytes(
        score_peak, 1U, sizeof(double) + sizeof(unsigned char));
    score_peak = verify_prepared_add_bytes(
        score_peak, (size_t)v->NT,
        sizeof(double) + sizeof(int));
    score_peak = verify_prepared_add_bytes(
        score_peak, nn_bytes, 1U);

    /*
     * Owner retirement keeps the score vectors while
     * verify_internal_fixup_theta() builds its expanded result and permutation
     * workspaces. These allocations do not overlap the scoring workspace, so
     * admission uses the larger peak.
     */
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NT,
        sizeof(double) + sizeof(int));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NTall,
        sizeof(int) + sizeof(double));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)v->NRall, sizeof(int));
    finish_peak = verify_prepared_add_bytes(
        finish_peak, (size_t)prepared->nrimage,
        3U * sizeof(double));
    transient_peak = MAX(score_peak, finish_peak);
    return verify_prepared_add_bytes(retained, transient_peak, 1U);
}

unsigned long long
verify_prepared_hit_work_units(const verify_prepared_hit_t* prepared) {
    unsigned long long nr;
    unsigned long long nt;

    if (!prepared || prepared->state != VERIFY_PREPARED_READY) {
        return 0U;
    }
    nr = (unsigned long long)prepared->verify.NR;
    nt = (unsigned long long)prepared->verify.NT;
    if (nr && nt > ULLONG_MAX / nr) {
        return ULLONG_MAX;
    }
    return nr * nt;
}

void verify_destroy_prepared_score(verify_prepared_score_t* score) {
    if (!score) {
        return;
    }
    free(score->allodds);
    free(score->theta);
    memset(score, 0, sizeof(*score));
}

void verify_destroy_prepared_hit(verify_prepared_hit_t* prepared) {
    verify_t* v;

    if (!prepared) {
        return;
    }
    v = &prepared->verify;
    free(prepared->refxyz);
    free(v->testperm);
    free(v->testsigma);
    free(v->tbadguys);
    free(v->refperm);
    free(v->refxy);
    free(v->refstarid);
    free(v->badguys);
    memset(prepared, 0, sizeof(*prepared));
    free(prepared);
}
