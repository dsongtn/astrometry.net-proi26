/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef VERIFY_H
#define VERIFY_H

#include <stddef.h>

#include "astrometry/kdtree.h"
#include "astrometry/matchobj.h"
#include "astrometry/bl.h"
#include "astrometry/starkd.h"
#include "astrometry/sip.h"
#include "astrometry/bl.h"
#include "astrometry/starxy.h"

struct verify_field_t {
    const starxy_t* field;
    // this copy is normal.
    double* xy;
    // this copy is permuted by the kdtree
    double* fieldcopy;
    kdtree_t* ftree;

    // should this field be spatially uniformized at the index's scale?
    anbool do_uniformize;
    // should this field be de-duplicated (have nearby sources removed)?
    anbool do_dedup;
    // apply radius-of-relevance filtering
    anbool do_ror;
};
typedef struct verify_field_t verify_field_t;


/*
 This function must be called once for each field before verification
 begins.  We build a kdtree out of the field stars (in pixel space)
 which will be used during deduplication.
 */
verify_field_t* verify_field_preprocess(const starxy_t* fieldxy);

/*
 This function must be called after all verification calls for a field
 are finished; we clean up the data structures we created in the
 verify_field_preprocess() function.
 */
void verify_field_free(verify_field_t* vf);




void verify_count_hits(int* theta, int besti, int* p_nmatch, int* p_nconflict, int* p_ndistractor);

void verify_wcs(const startree_t* skdt,
                int index_cutnside,
                const sip_t* sip,
                const verify_field_t* vf,
                double verify_pix2,
                double distractors,
                double fieldW,
                double fieldH,
                double logratio_tobail,
                double logratio_toaccept,
                double logratio_tostoplooking,

                double* logodds,
                int* nfield, int* nindex,
                int* nmatch, int* nconflict, int* ndistractor
                // int** theta ?
                );

/*
 Uses the following entries in the "mo" struct:
 -wcs_valid
 -wcstan
 -center
 -radius
 -field[]
 -star[]
 -dimquads

 Sets the following:
 -nfield
 -noverlap
 -nconflict
 -nindex
 -(matchobj_compute_derived() values)
 -logodds
 -corr_field
 -corr_index
 */
void verify_hit(const startree_t* skdt,
                int index_cutnside,
                // input/output param.
                MatchObj* mo,
                const sip_t* sip, // if non-NULL, verify this SIP WCS.
                const verify_field_t* vf,
                double verify_pix2,
                double distractors,
                double fieldW,
                double fieldH,
                double logratio_tobail,
                double logratio_toaccept,
                double logratio_tostoplooking,
                anbool distance_from_quad_bonus,
                anbool fake_match);

/*
 * Heap-owned result of the native StarKD traversal used by verification. The
 * logical result preserves the exact native reference-star order. Optional
 * sweep capture only attaches an owned copy; it does not change that result.
 * The source startree and its mappings must outlive the query. The query owns
 * no mapped payload pointer. A successful zero-result query still returns a
 * non-NULL context with count zero.
 */
typedef struct verify_index_query verify_index_query_t;

int verify_query_hit(const startree_t* skdt,
                     const double center[3],
                     double radius2,
                     verify_index_query_t** query);

size_t verify_index_query_count(const verify_index_query_t* query);
size_t verify_index_query_bytes(const verify_index_query_t* query);

/*
 * Return the exact mapped sweep element used by the indexed native result.
 * The returned mapping remains owned by skdt. Callers may reorder or
 * deduplicate physical ranges, but must not reorder the logical query result.
 */
int verify_index_query_sweep_range(const startree_t* skdt,
                                   const verify_index_query_t* query,
                                   size_t index,
                                   const void** data,
                                   size_t* size);

/*
 * Copy sweep values in native query-result order. This is intended to run
 * after the ranges above have completed delivery. The copy is optional:
 * continuation uses the authoritative mapped sweep when it is absent.
 */
int verify_index_query_capture_sweep(const startree_t* skdt,
                                     verify_index_query_t* query);

/*
 * One caller-owned direct-read buffer for an exact range of the compute
 * mapping. mapping_data is an address token in that mapping and is never
 * dereferenced by the capture routine. bytes contains the corresponding
 * contents. A buffer array must be ordered by mapping_data and must not
 * overlap.
 */
typedef struct verify_mapped_page_buffer {
    const void* mapping_data;
    size_t size;
    const unsigned char* bytes;
} verify_mapped_page_buffer_t;

/*
 * Copy sweep values from fully read mapping-relative buffers in native query
 * result order. The buffers remain caller-owned and are not retained. Every
 * requested sweep byte must be covered; failure leaves query unchanged so
 * the authoritative mapped fallback remains available.
 */
int verify_index_query_capture_sweep_buffers(
    const startree_t* skdt,
    verify_index_query_t* query,
    const verify_mapped_page_buffer_t* buffers,
    size_t buffer_count);

void verify_destroy_index_query(verify_index_query_t* query);

/*
 * Prepared verification separates index-backed input collection from the
 * index-free scoring kernel and the ordered MatchObj update.
 *
 * verify_prepare_hit() must run on the index owner. It performs every StarKD
 * lookup and field preparation step and retains no index, solver, callback,
 * or FITS mapping pointer. The resulting context is immutable until
 * verify_finish_prepared_hit() and may be scored by a foreign helper.
 *
 * verify_score_prepared_hit() calls the original serial verification
 * mathematics on prepared immutable data. The verify_field_t supplied during
 * preparation must outlive the prepared context. Independent contexts may be
 * scored concurrently only when verify_datalog_enabled() is false. Callers
 * must zero-initialize each verify_prepared_score_t before its first score.
 *
 * verify_finish_prepared_hit() must run on the owner in canonical candidate
 * order. It transfers accepted correspondence arrays to MatchObj exactly
 * once and consumes the score. Destroy functions accept partially consumed
 * objects and are always safe after a successful finish.
 */
typedef struct verify_prepared_hit verify_prepared_hit_t;

typedef struct verify_prepared_score {
    double logodds;
    double worstlogodds;
    int besti;
    int ibailed;
    int istopped;
    double* allodds;
    int* theta;
    anbool complete;
} verify_prepared_score_t;

int verify_prepare_hit(const startree_t* skdt,
                       int index_cutnside,
                       const MatchObj* mo,
                       const sip_t* sip,
                       const verify_field_t* vf,
                       double verify_pix2,
                       double distractors,
                       double fieldW,
                       double fieldH,
                       double logratio_tobail,
                       double logratio_toaccept,
                       double logratio_tostoplooking,
                       anbool distance_from_quad_bonus,
                       anbool fake_match,
                       verify_prepared_hit_t** prepared);

/*
 * Continue verification from one exact native StarKD query result.
 *
 * On success, this function consumes *query, sets it to NULL, and transfers
 * its arrays into *prepared. On validation or allocation failure it leaves
 * *query unchanged and caller-owned, sets *prepared to NULL, and returns -1.
 */
int verify_prepare_hit_from_query(const startree_t* skdt,
                                  verify_index_query_t** query,
                                  int index_cutnside,
                                  const MatchObj* mo,
                                  const sip_t* sip,
                                  const verify_field_t* vf,
                                  double verify_pix2,
                                  double distractors,
                                  double fieldW,
                                  double fieldH,
                                  double logratio_tobail,
                                  double logratio_toaccept,
                                  double logratio_tostoplooking,
                                  anbool distance_from_quad_bonus,
                                  anbool fake_match,
                                  verify_prepared_hit_t** prepared);

int verify_score_prepared_hit(const verify_prepared_hit_t* prepared,
                              verify_prepared_score_t* score);

int verify_finish_prepared_hit(verify_prepared_hit_t* prepared,
                               verify_prepared_score_t* score,
                               MatchObj* mo);

size_t verify_prepared_hit_bytes(const verify_prepared_hit_t* prepared);

size_t verify_prepared_score_bytes(
    const verify_prepared_hit_t* prepared);

size_t verify_prepared_hit_peak_bytes(
    const verify_prepared_hit_t* prepared);

unsigned long long
verify_prepared_hit_work_units(const verify_prepared_hit_t* prepared);

void verify_destroy_prepared_score(verify_prepared_score_t* score);
void verify_destroy_prepared_hit(verify_prepared_hit_t* prepared);

/*
 * Verification datalog records are process-global and ordered.  Callers that
 * would otherwise verify independent MatchObj instances concurrently must
 * retain serial execution while this stream is enabled.
 */
anbool verify_datalog_enabled(void);

// Distractor
#define THETA_DISTRACTOR -1
// Conflict
#define THETA_CONFLICT -2
// Filtered out
#define THETA_FILTERED -3
// Not examined because the bail-out threshold was reached.
#define THETA_BAILEDOUT -4
// Not examined because the stop-looking threshold was reached.
#define THETA_STOPPEDLOOKING -5

/*
 void verify_apply_ror(double* refxy, int* starids, int* p_NR,
 int index_cutnside,
 MatchObj* mo,
 const verify_field_t* vf,
 double pix2,
 double distractors,
 double fieldW,
 double fieldH,
 anbool do_gamma, anbool fake_match,
 double** p_testxy, double** p_sigma2s,
 int* p_NT, int** p_perm, double* p_effA,
 int* p_uninw, int* p_uninh);
 */

/**
 Returns the best log-odds encountered.
 */
double verify_star_lists(double* refxys, int NR,
                         const double* testxys, const double* testsigma2s, int NT,
                         double effective_area,
                         double distractors,
                         double logodds_bail,
                         double logodds_accept,
                         int* p_besti,
                         double** p_all_logodds, int** p_theta,
                         double* p_worstlogodds,
                         int** p_testperm);

void verify_get_uniformize_scale(int cutnside, double scale, int W, int H, int* uni_nw, int* uni_nh);

void verify_uniformize_field(const double* xy, int* perm, int N,
                             double fieldW, double fieldH,
                             int nw, int nh,
                             int** p_bincounts,
                             int** p_binids);

double* verify_uniformize_bin_centers(double fieldW, double fieldH,
                                      int nw, int nh);

void verify_get_quad_center(const verify_field_t* vf, const MatchObj* mo, double* centerpix,
                            double* quadr2);

/*
 int verify_get_test_stars(const verify_field_t* vf, MatchObj* mo,
 double pix2, anbool do_gamma,
 anbool fake_match,
 double** p_sigma2s, int** p_perm);
 */

void verify_get_index_stars(const double* fieldcenter, double fieldr2,
                            const startree_t* skdt, const sip_t* sip, const tan_t* tan,
                            double fieldW, double fieldH,
                            double** p_indexradec,
                            double** p_indexpix, int** p_starids, int* p_nindex);

/*
 anbool* verify_deduplicate_field_stars(const verify_field_t* vf, double* sigma2s, double nsigmas);
 */
/*
 double* verify_compute_sigma2s_arr(const double* xy, int NF,
 const double* qc, double Q2,
 double verify_pix2, anbool do_gamma);
 */

// For use with matchobj.h : matchodds
double verify_logodds_to_weight(double lodds);

void verify_free_matchobj(MatchObj* mo);

void verify_matchobj_deep_copy(const MatchObj* mo, MatchObj* dest);

double verify_get_ror2(double Q2, double area,
                       double distractors, int NR, double pix2);



double verify_star_lists_ror(double* refxys, int NR,
                             const double* testxys, const double* testsigma2s, int NT,
                             double pix2, double gamma,
                             const double* qc, double Q2,
                             double W, double H,
                             double distractors,
                             double logodds_bail,
                             double logodds_stoplooking,
                             int* p_besti,
                             double** p_all_logodds, int** p_theta,
                             double* p_worstlogodds,
                             int** p_testperm, int** p_refperm);

#endif
