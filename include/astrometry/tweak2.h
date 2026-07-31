/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef TWEAK_2_H
#define TWEAK_2_H

#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "astrometry/log.h"
#include "astrometry/sip.h"

typedef enum tweak2_status {
    TWEAK2_STATUS_SUCCESS = 0,
    TWEAK2_STATUS_FAILED_UNCLASSIFIED
} tweak2_status_t;

static inline const char* tweak2_status_name(tweak2_status_t status) {
    switch (status) {
    case TWEAK2_STATUS_SUCCESS:
        return "success";
    case TWEAK2_STATUS_FAILED_UNCLASSIFIED:
        return "failed-unclassified";
    }
    return "unknown";
}

/**
 Given an initial WCS solution, compute SIP polynomial distortions
 using an annealing-like strategy.  That is, it finds matches between
 image and reference catalog by searching within a radius, and that
 radius is small near a region of confidence, and grows as you move
 away.  That makes it possible to pick up more distant matches, but
 they are downweighted in the fit.  The annealing process reduces the
 slope of the growth of the matching radius with respect to the
 distance from the region of confidence.
 
 In Astrometry.net, the confidence region is the center of the quad
 that matched.

 fieldxy: source positions, x0,y0,x1,y1,...
 Nfield: number of source positions
 fieldjitter: standard deviation of field sources, in pixels;
 (FIXME: this should be per-source!)
 W, H: size of the image.
 indexradec: reference star positions, in degrees, ra0,dec0,ra1,dec1,...
 Nindex: number of reference stars
 indexjitter: standard deviation of reference sources, in arcsec.
 (FIXME: this should be per-source)
 quadcenter: the center of the quad that matched, in pixels.
 quadR2: the radius-squared, in pixels, of the quad that matched.
 distractors: how often you find unexpected test stars; we've always kept this fixed at 0.25
 logodds_bail: totally uninteresting; shouldn't really be here.  Set to, say, -100.
 sip_order: polynomial distortion order you want.  1=linear works.
 startwcs: initial WCS solution.  See util/sip.h : sip_wrap_tan() if you have a tan_t rather than a sip_t.
 destwcs: where to put the solution; NULL to allocate a new sip_t.
 newtheta: "theta" maps field stars to reference stars in the final matching that we produce.  Set this non-NULL to pull it out.
 newodds: this tells the confidence in the matches.  Use verify_logodds_to_weight() to turn these into a weight in [0,1].
 crpix: if you want to keep the reference point fixed, set this to a (2-element) array of the image reference position.
 */
sip_t* tweak2(const double* fieldxy, int Nfield,
              double fieldjitter,
              int W, int H,
              const double* indexradec, int Nindex,
              double indexjitter,
              const double* quadcenter, double quadR2,
              double distractors,
              double logodds_bail,
              int sip_order,
              int sip_invorder,
              const sip_t* startwcs,
              sip_t* destwcs,
              int** newtheta, double** newodds,
              double* crpix,
              double* p_logodds,
              int* p_besti,
              int* testperm, int startorder);

#define TWEAK2_TRACE_FNV1A_OFFSET UINT64_C(1469598103934665603)
#define TWEAK2_TRACE_FNV1A_PRIME UINT64_C(1099511628211)

static inline uint64_t tweak2_trace_hash_bytes(uint64_t hash,
                                               const void* data,
                                               size_t size) {
    const unsigned char* bytes = data;
    size_t i;

    if (!data || !size) {
        return hash;
    }
    for (i = 0; i < size; i++) {
        hash ^= bytes[i];
        hash *= TWEAK2_TRACE_FNV1A_PRIME;
    }
    return hash;
}

static inline uint64_t tweak2_trace_hash_int(uint64_t hash, int value) {
    return tweak2_trace_hash_bytes(hash, &value, sizeof(value));
}

static inline uint64_t tweak2_trace_hash_double(uint64_t hash, double value) {
    uint64_t bits;

    memcpy(&bits, &value, sizeof(bits));
    return tweak2_trace_hash_bytes(hash, &bits, sizeof(bits));
}

static inline uint64_t tweak2_trace_hash_int_array(const int* values,
                                                   int count) {
    uint64_t hash = TWEAK2_TRACE_FNV1A_OFFSET;

    hash = tweak2_trace_hash_int(hash, values != NULL);
    hash = tweak2_trace_hash_int(hash, count);
    if (values && count > 0 &&
        (size_t)count <= SIZE_MAX / sizeof(*values)) {
        hash = tweak2_trace_hash_bytes(
            hash, values, (size_t)count * sizeof(*values));
    }
    return hash;
}

static inline uint64_t tweak2_trace_hash_double_array(const double* values,
                                                      int count) {
    uint64_t hash = TWEAK2_TRACE_FNV1A_OFFSET;

    hash = tweak2_trace_hash_int(hash, values != NULL);
    hash = tweak2_trace_hash_int(hash, count);
    if (values && count > 0 &&
        (size_t)count <= SIZE_MAX / sizeof(*values)) {
        hash = tweak2_trace_hash_bytes(
            hash, values, (size_t)count * sizeof(*values));
    }
    return hash;
}

static inline uint64_t tweak2_trace_hash_polynomial(
    uint64_t hash,
    const double coefficients[SIP_MAXORDER][SIP_MAXORDER],
    int order) {
    int p;
    int q;

    hash = tweak2_trace_hash_int(hash, order);
    if (order < 0) {
        return hash;
    }
    for (p = 0; p <= order && p < SIP_MAXORDER; p++) {
        for (q = 0; q <= order - p && q < SIP_MAXORDER; q++) {
            hash = tweak2_trace_hash_double(hash, coefficients[p][q]);
        }
    }
    return hash;
}

static inline uint64_t tweak2_trace_hash_sip(const sip_t* sip) {
    uint64_t hash = TWEAK2_TRACE_FNV1A_OFFSET;
    int row;
    int column;

    hash = tweak2_trace_hash_int(hash, sip != NULL);
    if (!sip) {
        return hash;
    }
    for (row = 0; row < 2; row++) {
        hash = tweak2_trace_hash_double(hash, sip->wcstan.crval[row]);
        hash = tweak2_trace_hash_double(hash, sip->wcstan.crpix[row]);
        for (column = 0; column < 2; column++) {
            hash = tweak2_trace_hash_double(
                hash, sip->wcstan.cd[row][column]);
        }
    }
    hash = tweak2_trace_hash_double(hash, sip->wcstan.imagew);
    hash = tweak2_trace_hash_double(hash, sip->wcstan.imageh);
    hash = tweak2_trace_hash_int(hash, sip->wcstan.sin);
    hash = tweak2_trace_hash_polynomial(
        hash, sip->a, sip->a_order);
    hash = tweak2_trace_hash_polynomial(
        hash, sip->b, sip->b_order);
    hash = tweak2_trace_hash_polynomial(
        hash, sip->ap, sip->ap_order);
    hash = tweak2_trace_hash_polynomial(
        hash, sip->bp, sip->bp_order);
    return hash;
}

static inline uint64_t tweak2_trace_candidate_fingerprint(
    const double* fieldxy,
    int Nfield,
    double fieldjitter,
    int W,
    int H,
    const double* indexradec,
    int Nindex,
    double indexjitter,
    const double* quadcenter,
    double quadR2,
    double distractors,
    double logodds_bail,
    int sip_order,
    int sip_invorder,
    const sip_t* startwcs,
    const double* crpix,
    const int* testperm,
    int startorder,
    uint64_t* fieldxy_hash,
    uint64_t* indexradec_hash,
    uint64_t* startwcs_hash,
    uint64_t* testperm_hash) {
    uint64_t hash = TWEAK2_TRACE_FNV1A_OFFSET;
    int pair_count;

    pair_count = Nfield > 0 && Nfield <= INT32_MAX / 2
        ? Nfield * 2 : -1;
    *fieldxy_hash = tweak2_trace_hash_double_array(fieldxy, pair_count);
    pair_count = Nindex > 0 && Nindex <= INT32_MAX / 2
        ? Nindex * 2 : -1;
    *indexradec_hash = tweak2_trace_hash_double_array(
        indexradec, pair_count);
    *startwcs_hash = tweak2_trace_hash_sip(startwcs);
    *testperm_hash = tweak2_trace_hash_int_array(testperm, Nfield);

    hash = tweak2_trace_hash_int(hash, Nfield);
    hash = tweak2_trace_hash_int(hash, Nindex);
    hash = tweak2_trace_hash_int(hash, W);
    hash = tweak2_trace_hash_int(hash, H);
    hash = tweak2_trace_hash_int(hash, sip_order);
    hash = tweak2_trace_hash_int(hash, sip_invorder);
    hash = tweak2_trace_hash_int(hash, startorder);
    hash = tweak2_trace_hash_double(hash, fieldjitter);
    hash = tweak2_trace_hash_double(hash, indexjitter);
    hash = tweak2_trace_hash_double(hash, quadR2);
    hash = tweak2_trace_hash_double(hash, distractors);
    hash = tweak2_trace_hash_double(hash, logodds_bail);
    hash = tweak2_trace_hash_bytes(
        hash, fieldxy_hash, sizeof(*fieldxy_hash));
    hash = tweak2_trace_hash_bytes(
        hash, indexradec_hash, sizeof(*indexradec_hash));
    hash = tweak2_trace_hash_bytes(
        hash, startwcs_hash, sizeof(*startwcs_hash));
    hash = tweak2_trace_hash_bytes(
        hash, testperm_hash, sizeof(*testperm_hash));
    hash = tweak2_trace_hash_int(hash, quadcenter != NULL);
    if (quadcenter) {
        hash = tweak2_trace_hash_double(hash, quadcenter[0]);
        hash = tweak2_trace_hash_double(hash, quadcenter[1]);
    }
    hash = tweak2_trace_hash_int(hash, crpix != NULL);
    if (crpix) {
        hash = tweak2_trace_hash_double(hash, crpix[0]);
        hash = tweak2_trace_hash_double(hash, crpix[1]);
    }
    return hash;
}

static inline sip_t* tweak2_with_status(
    const double* fieldxy,
    int Nfield,
    double fieldjitter,
    int W,
    int H,
    const double* indexradec,
    int Nindex,
    double indexjitter,
    const double* quadcenter,
    double quadR2,
    double distractors,
    double logodds_bail,
    int sip_order,
    int sip_invorder,
    const sip_t* startwcs,
    sip_t* destwcs,
    int** newtheta,
    double** newodds,
    double* crpix,
    double* p_logodds,
    int* p_besti,
    int* testperm,
    int startorder,
    tweak2_status_t* status_out) {
    uint64_t candidate = 0;
    uint64_t fieldxy_hash = 0;
    uint64_t indexradec_hash = 0;
    uint64_t startwcs_hash = 0;
    uint64_t testperm_hash = 0;
    uint64_t outputwcs_hash = 0;
    uint64_t theta_hash = 0;
    uint64_t odds_hash = 0;
    tweak2_status_t status;
    sip_t* result;
    int trace_enabled = log_get_level() >= LOG_VERB;

    if (trace_enabled) {
        candidate = tweak2_trace_candidate_fingerprint(
            fieldxy, Nfield, fieldjitter, W, H,
            indexradec, Nindex, indexjitter,
            quadcenter, quadR2, distractors, logodds_bail,
            sip_order, sip_invorder, startwcs, crpix,
            testperm, startorder,
            &fieldxy_hash, &indexradec_hash,
            &startwcs_hash, &testperm_hash);
        logverb("[tune-fingerprint] phase=enter candidate=%016" PRIx64
                " fieldxy=%016" PRIx64
                " indexradec=%016" PRIx64
                " startwcs=%016" PRIx64
                " testperm=%016" PRIx64
                " nfield=%i nindex=%i order=%i invorder=%i"
                " startorder=%i\n",
                candidate, fieldxy_hash, indexradec_hash,
                startwcs_hash, testperm_hash,
                Nfield, Nindex, sip_order, sip_invorder, startorder);
    }

    result = tweak2(fieldxy, Nfield,
                    fieldjitter,
                    W, H,
                    indexradec, Nindex,
                    indexjitter,
                    quadcenter, quadR2,
                    distractors,
                    logodds_bail,
                    sip_order,
                    sip_invorder,
                    startwcs,
                    destwcs,
                    newtheta, newodds,
                    crpix,
                    p_logodds,
                    p_besti,
                    testperm, startorder);
    status = result
        ? TWEAK2_STATUS_SUCCESS
        : TWEAK2_STATUS_FAILED_UNCLASSIFIED;
    if (status_out) {
        *status_out = status;
    }

    if (trace_enabled) {
        outputwcs_hash = tweak2_trace_hash_sip(result);
        if (result && newtheta && *newtheta) {
            theta_hash = tweak2_trace_hash_int_array(*newtheta, Nfield);
        }
        if (result && newodds && *newodds) {
            odds_hash = tweak2_trace_hash_double_array(*newodds, Nfield);
        }
        logverb("[tune-fingerprint] phase=exit candidate=%016" PRIx64
                " status=%s outputwcs=%016" PRIx64
                " theta=%016" PRIx64
                " odds=%016" PRIx64
                " besti=%i logodds=%.17g\n",
                candidate, tweak2_status_name(status),
                outputwcs_hash, theta_hash, odds_hash,
                p_besti ? *p_besti : -1,
                p_logodds ? *p_logodds : 0.0);
    }
    return result;
}

static inline sip_t* tweak2_traced(
    const double* fieldxy,
    int Nfield,
    double fieldjitter,
    int W,
    int H,
    const double* indexradec,
    int Nindex,
    double indexjitter,
    const double* quadcenter,
    double quadR2,
    double distractors,
    double logodds_bail,
    int sip_order,
    int sip_invorder,
    const sip_t* startwcs,
    sip_t* destwcs,
    int** newtheta,
    double** newodds,
    double* crpix,
    double* p_logodds,
    int* p_besti,
    int* testperm,
    int startorder) {
    return tweak2_with_status(
        fieldxy, Nfield, fieldjitter, W, H,
        indexradec, Nindex, indexjitter,
        quadcenter, quadR2, distractors, logodds_bail,
        sip_order, sip_invorder, startwcs, destwcs,
        newtheta, newodds, crpix, p_logodds, p_besti,
        testperm, startorder, NULL);
}

/*
 * Existing callers are redirected only after the raw external declaration and
 * wrapper bodies above have been parsed.  solver/tweak2.c does not include
 * this header, so the implementation symbol remains unchanged.
 */
#define tweak2 tweak2_traced

#endif
