/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef VERIFY_INTERNAL_H
#define VERIFY_INTERNAL_H

#include <stddef.h>

#include "index_shard_internal.h"
#include "verify.h"

typedef struct verify_s {
    const sip_t* wcs;

    int NR;
    int NRall;
    int* refperm;
    int* refstarid;
    double* refxy;
    int* badguys;

    int NT;
    int NTall;
    int* testperm;
    double* testxy;
    double* testsigma;
    int* tbadguys;
} verify_t;

/*
 * verify.o is also linked into legacy tools that do not include
 * index_shard.o. Weak references keep those tools on the original inline
 * path while the solver engine supplies the strong helper scheduler.
 */
#if defined(__GNUC__) && !defined(_WIN32)
extern size_t index_shard_helper_available_workers(void)
    __attribute__((weak));
extern index_shard_helper_run_status_t index_shard_helper_run(
    const index_shard_helper_ops_t* ops,
    index_shard_helper_task_t* tasks,
    size_t task_count,
    index_shard_helper_run_stats_t* stats) __attribute__((weak));
extern anbool index_shard_worker_stop_requested(void)
    __attribute__((weak));

static inline anbool verify_internal_worker_stop_requested(void) {
    return index_shard_worker_stop_requested &&
        index_shard_worker_stop_requested();
}
#define VERIFY_INTERNAL_HELPERS_LINKED 1
#else
static inline anbool verify_internal_worker_stop_requested(void) {
    return FALSE;
}
#define VERIFY_INTERNAL_HELPERS_LINKED 0
#endif

int verify_internal_get_test_stars(
    verify_t* v,
    const verify_field_t* vf,
    MatchObj* mo,
    double pix2,
    anbool do_gamma,
    anbool fake_match);

int verify_internal_apply_ror(verify_t* v,
                              int index_cutnside,
                              MatchObj* mo,
                              const verify_field_t* vf,
                              double pix2,
                              double distractors,
                              double fieldW,
                              double fieldH,
                              anbool do_gamma,
                              anbool fake_match,
                              double* p_effA,
                              int* p_uninw,
                              int* p_uninh);

double verify_internal_star_lists(verify_t* v,
                                  double effective_area,
                                  double distractors,
                                  double logodds_bail,
                                  double logodds_stoplooking,
                                  int* p_besti,
                                  double** p_logodds,
                                  int** p_theta,
                                  double* p_worstlogodds,
                                  int* p_ibailed,
                                  int* p_istopped,
                                  anbool* p_completed);

/* Exact transient bound for the NN implementation above. */
size_t verify_internal_score_workspace_bytes(int npoints);

anbool verify_internal_filter_stars_in_field_parallel(
    const sip_t* sip,
    const tan_t* tan,
    const double* xyz,
    int npoints,
    double** p_xy,
    int** p_inbounds,
    int* p_ngood);

void verify_internal_set_null_mo(MatchObj* mo);

int verify_internal_fixup_theta(int* theta,
                                double* allodds,
                                int ibailed,
                                int istopped,
                                verify_t* v,
                                int besti,
                                int nrimage,
                                double* refxyz,
                                int** p_etheta,
                                double** p_eodds);

#endif
