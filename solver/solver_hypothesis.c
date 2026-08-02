/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */
/* AB-pair enumeration, descriptors, and bounded verification waves. */

#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#include "solver.h"
#include "verify.h"
#include "fit-wcs.h"
#include "fitsbin.h"
#include "log.h"
#include "mathutil.h"
#include "quad-utils.h"
#include "tic.h"
#include "index_shard_internal.h"
#include "solver_field_geometry_internal.h"
#include "solver_inline_internal.h"
#include "solver_hypothesis_internal.h"
static uint64_t solver_order_double_bits(double value) {
    uint64_t bits;

    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

uint64_t solver_hypothesis_order_digest(
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    uint64_t digest = 0;
    int dimcode = (dimquad - NBACK) * 2;
    int i;

    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x4859504f54484553));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)dimquad);
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)current_parity);
    for (i = 0; i < dimquad; i++) {
        digest = solver_order_hash_mix(
            digest,
            (uint64_t)(unsigned int)stars[i]);
    }
    for (i = 0; i < dimcode; i++) {
        digest = solver_order_hash_mix(
            digest,
            solver_order_double_bits(code[i]));
    }
    return digest;
}

uint64_t solver_kd_result_order_digest(
    const kdtree_qres_t* result) {
    uint64_t digest = 0;
    int i;

    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x4b44524553554c54));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)result->nres);
    for (i = 0; i < result->nres; i++) {
        digest = solver_order_hash_mix(
            digest,
            (uint64_t)(unsigned int)result->inds[i]);
        digest = solver_order_hash_mix(
            digest,
            solver_order_double_bits(result->sdists[i]));
    }
    return digest;
}

void solver_record_candidate_order(
    solver_t* solver,
    solver_ab_candidate_action_t action,
    int quadno,
    double code_err) {
    uint64_t digest = 0;

    if (!solver->profile.detailed) {
        return;
    }
    digest = solver_order_hash_mix(
        digest,
        UINT64_C(0x43414e4449444154));
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)action);
    digest = solver_order_hash_mix(
        digest,
        (uint64_t)(unsigned int)quadno);
    digest = solver_order_hash_mix(
        digest,
        solver_order_double_bits(code_err));
    solver->profile.candidate_order_hash =
        solver_order_hash_mix(
            solver->profile.candidate_order_hash,
            digest);
}

anbool solver_ab_poll_phase_stop(
    solver_t* solver,
    time_t* next_timer_callback_time) {
    time_t delay;
    time_t now;

    if (solver->quit_now ||
        solver_poll_worker_stop(solver)) {
        return TRUE;
    }
    if (!solver->timer_callback ||
        !next_timer_callback_time) {
        return FALSE;
    }
    now = time(NULL);
    if (now <= *next_timer_callback_time) {
        return FALSE;
    }
    update_timeused(solver);
    delay = solver->timer_callback(solver->userdata);
    if (!delay || solver->quit_now) {
        solver->quit_now = TRUE;
        return TRUE;
    }
    *next_timer_callback_time = now + delay;
    return FALSE;
}

/*
 * Descriptor planning is index-free. A synchronous packet may borrow a const
 * CodeKD tree view; QuadFile, StarKD, lazy index initialization, verification,
 * counters, callbacks, and reducer state remain owner-local.
 */
static double solver_ab_timeval_seconds(
    const struct timeval* value) {
    return (double)value->tv_sec +
        (double)value->tv_usec * 1.0e-6;
}

void solver_ab_phase_telemetry_begin(
    const solver_t* solver,
    solver_ab_phase_telemetry_t* telemetry) {
    memset(telemetry, 0, sizeof(*telemetry));
    if (!solver ||
        log_get_level() < LOG_VERB ||
        pl_size(solver->indexes) != 1U) {
        return;
    }
    telemetry->enabled = TRUE;
    telemetry->wall_start = monotonic_seconds();
    telemetry->resource_valid =
        getrusage(
            RUSAGE_SELF,
            &telemetry->resource_start) == 0;
    telemetry->combinations_start = solver->numtries;
    telemetry->candidates_start = solver->nummatches;
    telemetry->codekd_calls_start =
        solver->profile.codekd_calls;
    telemetry->codekd_hits_start =
        solver->profile.codekd_hits;
    telemetry->verify_calls_start =
        solver->profile.verify_calls;
}

void solver_ab_phase_telemetry_report(
    const solver_t* solver,
    const solver_ab_phase_telemetry_t* telemetry,
    int newpoint,
    solver_ab_phase_kind_t phase,
    solver_ab_phase_mode_t mode) {
    struct rusage resource_end;
    anbool resource_valid;
    double user_seconds = 0.0;
    double system_seconds = 0.0;
    long major_faults = 0;
    const char* index_name;

    if (!solver || !telemetry || !telemetry->enabled) {
        return;
    }
    resource_valid =
        telemetry->resource_valid &&
        getrusage(RUSAGE_SELF, &resource_end) == 0;
    if (resource_valid) {
        user_seconds =
            solver_ab_timeval_seconds(&resource_end.ru_utime) -
            solver_ab_timeval_seconds(
                &telemetry->resource_start.ru_utime);
        system_seconds =
            solver_ab_timeval_seconds(&resource_end.ru_stime) -
            solver_ab_timeval_seconds(
                &telemetry->resource_start.ru_stime);
        major_faults =
            resource_end.ru_majflt -
            telemetry->resource_start.ru_majflt;
    }
    index_name =
        solver->index && solver->index->indexname
            ? solver->index->indexname
            : "(unknown)";
    logverb("[solver-ab-phase] index=%s object=%i phase=%s mode=%s "
            "combinations=%llu codekd_queries=%llu codekd_hits=%llu "
            "candidates=%llu verifications=%llu "
            "wall=%.6f user=%.6f system=%.6f "
            "major_faults=%ld resource=%s "
            "hypothesis_order=%016llx kd_result_order=%016llx "
            "candidate_order=%016llx\n",
            index_name,
            newpoint + 1,
            phase == SOLVER_AB_PHASE_DIAGONAL
                ? "diagonal"
                : "off-diagonal",
            mode == SOLVER_AB_MODE_ASSISTED
                ? "assisted"
                : (mode == SOLVER_AB_MODE_FLATTENED_OWNER
                    ? "flattened-owner"
                    : "native"),
            (unsigned long long)(
                solver->numtries -
                telemetry->combinations_start),
            solver->profile.codekd_calls -
                telemetry->codekd_calls_start,
            solver->profile.codekd_hits -
                telemetry->codekd_hits_start,
            (unsigned long long)(
                solver->nummatches -
                telemetry->candidates_start),
            solver->profile.verify_calls -
                telemetry->verify_calls_start,
            monotonic_seconds() - telemetry->wall_start,
            user_seconds,
            system_seconds,
            major_faults,
            resource_valid
                ? "process-overlap"
                : "unavailable",
            solver->profile.hypothesis_order_hash,
            solver->profile.kd_result_order_hash,
            solver->profile.candidate_order_hash);
}

unsigned long long solver_ab_saturating_add(
    unsigned long long a,
    unsigned long long b) {
    if (ULLONG_MAX - a < b) {
        return ULLONG_MAX;
    }
    return a + b;
}

static unsigned long long solver_ab_saturating_choose(int n, int k) {
    unsigned long long result = 1;
    int i;

    if (n < 0 || k < 0 || k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    for (i = 1; i <= k; i++) {
        unsigned long long numerator =
            (unsigned long long)(n - k + i);

        if (result > ULLONG_MAX / numerator) {
            return ULLONG_MAX;
        }
        result *= numerator;
        result /= (unsigned long long)i;
    }
    return result;
}

solver_packet_reserve_result_t
solver_verification_packet_reserve(
    solver_verification_packet_t* packet,
    size_t required) {
    solver_ab_candidate_t* resized;
    size_t capacity;
    size_t max_capacity =
        SOLVER_AB_CANDIDATE_LIMIT_BYTES / sizeof(*resized);
    size_t bytes;

    if (required <= packet->candidate_capacity) {
        return SOLVER_PACKET_RESERVE_OK;
    }
    if (required > max_capacity || !max_capacity) {
        return SOLVER_PACKET_RESERVE_FULL;
    }
    capacity = packet->candidate_capacity ?
        packet->candidate_capacity : MIN(16U, max_capacity);
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            packet->allocation_failed = TRUE;
            return SOLVER_PACKET_RESERVE_ERROR;
        }
        if (capacity > max_capacity / 2U) {
            capacity = max_capacity;
        } else {
            capacity *= 2U;
        }
    }
    if (capacity > SIZE_MAX / sizeof(*resized)) {
        packet->allocation_failed = TRUE;
        return SOLVER_PACKET_RESERVE_ERROR;
    }
    bytes = capacity * sizeof(*resized);
    if (bytes > SOLVER_AB_CANDIDATE_LIMIT_BYTES) {
        return SOLVER_PACKET_RESERVE_FULL;
    }

    resized = realloc(packet->candidates, bytes);
    if (!resized) {
        packet->allocation_failed = TRUE;
        return SOLVER_PACKET_RESERVE_ERROR;
    }
    packet->candidates = resized;
    packet->candidate_capacity = capacity;
    return SOLVER_PACKET_RESERVE_OK;
}

void solver_verification_packet_free(solver_verification_packet_t* packet) {
    if (!packet) {
        return;
    }
    free(packet->candidates);
    memset(packet, 0, sizeof(*packet));
}

static anbool solver_ab_cancelled(
    const solver_ab_descriptor_planner_t* executor) {
    (void)executor;
    return index_shard_worker_stop_requested();
}

static int solver_ab_candidate_prepare(
    solver_ab_candidate_t* candidate,
    const kdtree_qres_t* result,
    int result_index,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    const solver_ab_snapshot_t* snapshot,
    anbool current_parity) {
    double starxyz[DQMAX * 3];
    double scale;
    double arcsecperpix;
    double abscale;
    tan_t wcs;
    unsigned int star[DQMAX];
    int thisquadno;
    int i;
    anbool outofbounds = FALSE;

    memset(candidate, 0, sizeof(*candidate));
    candidate->action = SOLVER_AB_CANDIDATE_SCALE_SKIP;

    thisquadno = result->inds[result_index];
    candidate->quadno = thisquadno;
    candidate->code_err = result->sdists[result_index];
    if (quadfile_get_stars(
            snapshot->index->quads,
            thisquadno,
            star)) {
        return -1;
    }

    if (snapshot->use_radec) {
        for (i = 0; i < dimquads; i++) {
            if (startree_get(
                    snapshot->index->starkd,
                    star[i],
                    starxyz + 3 * i)) {
                return -1;
            }
            if (distsq(starxyz + 3 * i,
                       snapshot->centerxyz,
                       3) > snapshot->r2) {
                outofbounds = TRUE;
                break;
            }
        }
        if (outofbounds) {
            candidate->action =
                SOLVER_AB_CANDIDATE_RADEC_SKIP;
            return 0;
        }
    } else {
        if (startree_get(
                snapshot->index->starkd,
                star[0],
                starxyz) ||
            startree_get(
                snapshot->index->starkd,
                star[1],
                starxyz + 3)) {
            return -1;
        }
    }

    abscale =
        square(distsq2rad(distsq(starxyz, starxyz + 3, 3))) /
        distsq(field_xy, field_xy + 2, 2);
    if (abscale > snapshot->abscale_high ||
        abscale < snapshot->abscale_low) {
        candidate->action =
            SOLVER_AB_CANDIDATE_ABSCALE_SKIP;
        return 0;
    }

    if (!snapshot->use_radec) {
        for (i = 2; i < dimquads; i++) {
            if (startree_get(
                    snapshot->index->starkd,
                    star[i],
                    starxyz + 3 * i)) {
                return -1;
            }
        }
    }

    if (fit_tan_wcs(
            starxyz,
            field_xy,
            dimquads,
            &wcs,
            &scale)) {
        candidate->action = SOLVER_AB_CANDIDATE_BAD_QUAD;
        return 0;
    }
    arcsecperpix = scale * 3600.0;
    if (arcsecperpix > snapshot->funits_upper ||
        arcsecperpix < snapshot->funits_lower) {
        candidate->action = SOLVER_AB_CANDIDATE_SCALE_SKIP;
        return 0;
    }

    memcpy(&candidate->wcs, &wcs, sizeof(tan_t));
    candidate->scale = arcsecperpix;
    candidate->parity = current_parity;
    candidate->quad_npeers = result->nres;
    for (i = 0; i < dimquads; i++) {
        candidate->star[i] = star[i];
        candidate->field[i] = fieldstars[i];
    }
    memcpy(candidate->quadpix,
           field_xy,
           (size_t)2 * (size_t)dimquads * sizeof(double));
    memcpy(candidate->quadxyz,
           starxyz,
           (size_t)3 * (size_t)dimquads * sizeof(double));
    candidate->action = SOLVER_AB_CANDIDATE_VERIFY;
    return 0;
}

static int solver_ab_record_descriptor(
    solver_ab_builder_t* builder,
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    solver_ab_descriptor_output_t* output =
        builder->descriptor_output;
    solver_ab_descriptor_t* descriptor;
    int dimcode = (dimquad - NBACK) * 2;

    if (!output || !stars || !code ||
        dimquad < NBACK || dimquad > DQMAX ||
        dimcode < 0 || dimcode > DCMAX ||
        output->descriptor_count >=
            SOLVER_AB_DESCRIPTOR_CAPACITY) {
        builder->fatal_error = TRUE;
        if (builder->status) {
            builder->status->evaluation_failed = TRUE;
        }
        return -1;
    }
    descriptor = &output->descriptors[
        output->descriptor_count++];
    memcpy(
        descriptor->stars,
        stars,
        (size_t)dimquad * sizeof(*stars));
    memcpy(
        descriptor->code,
        code,
        (size_t)dimcode * sizeof(*code));
    descriptor->tol2 = builder->tol2;
    descriptor->rel_field_noise2 =
        builder->rel_field_noise2;
    descriptor->numtries_delta =
        builder->pending_numtries;
    descriptor->cxdx_delta = builder->pending_cxdx;
    descriptor->meanx_delta = builder->pending_meanx;
    descriptor->current_parity = current_parity;
    builder->pending_numtries = 0U;
    builder->pending_cxdx = 0U;
    builder->pending_meanx = 0U;
    return 0;
}

static int solver_ab_record_hypothesis(
    solver_ab_builder_t* builder,
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity) {
    if (!builder || !builder->status ||
        index_shard_worker_stop_requested()) {
        if (builder && builder->status) {
            builder->status->cancelled = TRUE;
        }
        return -1;
    }
    return solver_ab_record_descriptor(
        builder,
        stars,
        code,
        dimquad,
        current_parity);
}

static void solver_ab_try_permutations(
    const int* origstars,
    int dimquad,
    const double* origcode,
    solver_ab_builder_t* builder,
    anbool current_parity,
    int* stars,
    double* code,
    int slot,
    anbool* placed) {
    const solver_ab_snapshot_t* snapshot =
        &builder->planner->snapshot;
    double mycode[DCMAX];
    int nstars = dimquad - NBACK;
    int lastslot = dimquad - NBACK - 1;
    int i;

    if (!code) {
        code = mycode;
    }
    if (slot >= DCMAX / 2 ||
        solver_ab_cancelled(builder->planner) ||
        builder->fatal_error) {
        if (solver_ab_cancelled(builder->planner)) {
            builder->status->cancelled = TRUE;
        }
        return;
    }

    for (i = 0; i < nstars; i++) {
        if (placed[i]) {
            continue;
        }
        if (slot > 0 &&
            snapshot->cx_less_than_dx &&
            code[2 * (slot - 1)] >
                origcode[2 * i] + snapshot->cxdx_margin) {
            builder->pending_cxdx++;
            continue;
        }

        stars[slot + NBACK] = origstars[i + NBACK];
        code[2 * slot] = origcode[2 * i];
        code[2 * slot + 1] = origcode[2 * i + 1];

        if (snapshot->cx_less_than_dx &&
            snapshot->meanx_less_than_half) {
            double meanx = 0.0;
            int j;

            for (j = 0; j <= slot; j++) {
                meanx += code[2 * j];
            }
            meanx /= (double)(slot + 1);
            if (meanx > 0.5 + snapshot->cxdx_margin) {
                builder->pending_meanx++;
                continue;
            }
        }

        if (slot < lastslot) {
            placed[i] = TRUE;
            solver_ab_try_permutations(
                origstars,
                dimquad,
                origcode,
                builder,
                current_parity,
                stars,
                code,
                slot + 1,
                placed);
            placed[i] = FALSE;
            if (builder->fatal_error ||
                builder->status->cancelled) {
                return;
            }
        } else if (solver_ab_record_hypothesis(
                       builder,
                       stars,
                       code,
                       dimquad,
                       current_parity)) {
            return;
        }
    }
}

static void solver_ab_try_all_codes_2(
    const int* fieldstars,
    int dimquad,
    const double* code,
    solver_ab_builder_t* builder,
    anbool current_parity) {
    double flipcode[DCMAX];
    int stars[DQMAX];
    anbool placed[DQMAX];
    int dimcode = (dimquad - NBACK) * 2;
    int i;

    stars[0] = fieldstars[0];
    stars[1] = fieldstars[1];
    memset(placed, 0, sizeof(placed));
    solver_ab_try_permutations(
        fieldstars,
        dimquad,
        code,
        builder,
        current_parity,
        stars,
        NULL,
        0,
        placed);
    if (builder->fatal_error ||
        builder->status->cancelled) {
        return;
    }

    stars[0] = fieldstars[1];
    stars[1] = fieldstars[0];
    for (i = 0; i < dimcode; i++) {
        flipcode[i] = 1.0 - code[i];
    }
    memset(placed, 0, sizeof(placed));
    solver_ab_try_permutations(
        fieldstars,
        dimquad,
        flipcode,
        builder,
        current_parity,
        stars,
        NULL,
        0,
        placed);
}

static void solver_ab_try_all_codes(
    const solver_pair_geometry_t* pair_geometry,
    int field_a,
    int field_b,
    const int* fieldstars,
    int dimquad,
    solver_ab_builder_t* builder) {
    const solver_field_geometry_t* geometry =
        builder->planner->field_geometry;
    const solver_ab_snapshot_t* snapshot =
        &builder->planner->snapshot;
    double code[DCMAX];
    double flipcode[DCMAX];
    int dimcode = (dimquad - 2) * 2;
    int i;

    builder->pending_numtries++;
    for (i = 0; i < dimquad - NBACK; i++) {
        if (!solver_pair_geometry_transform(
                geometry,
                pair_geometry,
                field_a,
                field_b,
                fieldstars[NBACK + i],
                &code[2 * i],
                &code[2 * i + 1])) {
            builder->fatal_error = TRUE;
            builder->status->evaluation_failed = TRUE;
            return;
        }
    }

    if (snapshot->parity == PARITY_NORMAL ||
        snapshot->parity == PARITY_BOTH) {
        solver_ab_try_all_codes_2(
            fieldstars,
            dimquad,
            code,
            builder,
            FALSE);
    }
    if (builder->fatal_error ||
        builder->status->cancelled) {
        return;
    }
    if (snapshot->parity == PARITY_FLIP ||
        snapshot->parity == PARITY_BOTH) {
        quad_flip_parity(code, flipcode, dimcode);
        solver_ab_try_all_codes_2(
            fieldstars,
            dimquad,
            flipcode,
            builder,
            TRUE);
    }
}

typedef anbool (*solver_ab_combination_visitor_t)(
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    int* field,
    int dimquad,
    void* opaque);

static int solver_ab_select_combination(
    const int* eligible,
    int eligible_count,
    int* field,
    int fieldoffset,
    int n_to_add,
    int* selection,
    unsigned long long ordinal) {
    int position;

    for (position = 0; position < n_to_add; position++) {
        int remaining = n_to_add - position - 1;
        int bottom = position ?
            selection[position - 1] + 1 : 0;
        int candidate_index;
        anbool selected = FALSE;

        for (candidate_index = bottom;
             candidate_index < eligible_count;
             candidate_index++) {
            unsigned long long suffix_count;

            suffix_count = solver_ab_saturating_choose(
                eligible_count - candidate_index - 1,
                remaining);
            if (ordinal >= suffix_count) {
                ordinal -= suffix_count;
                continue;
            }
            selection[position] = candidate_index;
            field[fieldoffset + position] =
                eligible[candidate_index];
            selected = TRUE;
            break;
        }
        if (!selected) {
            return -1;
        }
    }
    return ordinal == 0U ? 0 : -1;
}

static anbool solver_ab_next_combination(
    const int* eligible,
    int eligible_count,
    int* field,
    int fieldoffset,
    int n_to_add,
    int* selection) {
    int position;

    for (position = n_to_add - 1;
         position >= 0;
         position--) {
        int maximum =
            eligible_count - (n_to_add - position);
        int fill;

        if (selection[position] >= maximum) {
            continue;
        }
        selection[position]++;
        field[fieldoffset + position] =
            eligible[selection[position]];
        for (fill = position + 1;
             fill < n_to_add;
             fill++) {
            selection[fill] = selection[fill - 1] + 1;
            field[fieldoffset + fill] =
                eligible[selection[fill]];
        }
        return TRUE;
    }
    return FALSE;
}

static int* solver_ab_get_eligible_workspace(
    solver_ab_descriptor_planner_t* planner,
    size_t required) {
    if (!planner || !planner->combination_eligible ||
        required > planner->combination_eligible_capacity) {
        return NULL;
    }
    return planner->combination_eligible;
}

static int solver_ab_visit_pair_range(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    solver_ab_combination_visitor_t visitor,
    void* opaque) {
    int field[DQMAX] = {0};
    int selection[DQMAX] = {0};
    int* eligible;
    int eligible_count = 0;
    int fieldoffset;
    int n_to_add;
    int i;
    unsigned long long ordinal;

    if (local_first >= local_end ||
        local_end > pair->combination_count) {
        return -1;
    }
    field[A] = pair->field_a;
    field[B] = pair->field_b;
    if (executor->phase == SOLVER_AB_PHASE_DIAGONAL) {
        fieldoffset = C;
        n_to_add = executor->dimquads - 2;
    } else {
        field[C] = executor->newpoint;
        fieldoffset = D;
        n_to_add = executor->dimquads - 3;
    }
    if (n_to_add < 0 ||
        fieldoffset + n_to_add > DQMAX) {
        return -1;
    }
    if (n_to_add == 0) {
        if (local_first != 0U ||
            local_end != 1U) {
            return -1;
        }
        return visitor(
            pair_geometry,
            pair,
            field,
            executor->dimquads,
            opaque) ? 1 : 0;
    }
    eligible = solver_ab_get_eligible_workspace(
        executor,
        (size_t)executor->newpoint);
    if (!eligible) {
        return -1;
    }
    for (i = 0; i < executor->newpoint; i++) {
        double x;
        double y;

        if (solver_pair_geometry_transform(
                executor->field_geometry,
                pair_geometry,
                pair->field_a,
                pair->field_b,
                i,
                &x,
                &y)) {
            eligible[eligible_count++] = i;
        }
    }
    if (solver_ab_saturating_choose(
            eligible_count,
            n_to_add) != pair->combination_count) {
        return -1;
    }
    if (solver_ab_select_combination(
            eligible,
            eligible_count,
            field,
            fieldoffset,
            n_to_add,
            selection,
            local_first)) {
        return -1;
    }
    for (ordinal = local_first;
         ordinal < local_end;
         ordinal++) {
        if (visitor(
                pair_geometry,
                pair,
                field,
                executor->dimquads,
                opaque)) {
            return 1;
        }
        if (ordinal + 1U < local_end &&
            !solver_ab_next_combination(
                eligible,
                eligible_count,
                field,
                fieldoffset,
                n_to_add,
                selection)) {
            return -1;
        }
    }
    return 0;
}

static size_t solver_ab_find_pair_for_work(
    const solver_ab_pair_t* pairs,
    size_t pair_count,
    unsigned long long work) {
    size_t low = 0U;
    size_t high = pair_count;

    while (low < high) {
        size_t middle = low + (high - low) / 2U;

        if (pairs[middle].combination_first <= work) {
            low = middle + 1U;
        } else {
            high = middle;
        }
    }
    return low ? low - 1U : 0U;
}

static anbool solver_ab_builder_visit(
    const solver_pair_geometry_t* pair_geometry,
    const solver_ab_pair_t* pair,
    int* field,
    int dimquad,
    void* opaque) {
    solver_ab_builder_t* builder = opaque;

    if (solver_ab_cancelled(builder->planner)) {
        builder->status->cancelled = TRUE;
        return TRUE;
    }
    solver_ab_try_all_codes(
        pair_geometry,
        pair->field_a,
        pair->field_b,
        field,
        dimquad,
        builder);
    return builder->fatal_error ||
        builder->status->cancelled;
}

typedef int (*solver_ab_pair_range_visitor_t)(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque);

static int solver_ab_walk_task_ranges(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_task_t* task,
    solver_ab_pair_range_visitor_t visitor,
    void* opaque) {
    unsigned long long cursor;
    size_t pair_index;

    if (!executor || !task || !visitor ||
        !executor->pairs || !executor->pair_count ||
        task->combination_first >=
            task->combination_end) {
        return -1;
    }
    cursor = task->combination_first;
    pair_index = solver_ab_find_pair_for_work(
        executor->pairs,
        executor->pair_count,
        cursor);
    for (;
         pair_index < executor->pair_count &&
             cursor < task->combination_end;
         pair_index++) {
        const solver_ab_pair_t* pair =
            &executor->pairs[pair_index];
        unsigned long long pair_end;
        unsigned long long slice_end;
        const solver_pair_geometry_t* pair_geometry;
        int visit_result;

        if (!pair->combination_count ||
            ULLONG_MAX - pair->combination_first <
                pair->combination_count) {
            return -1;
        }
        pair_end = pair->combination_first +
            pair->combination_count;
        if (pair_end <= cursor) {
            continue;
        }
        if (pair->combination_first > cursor) {
            return -1;
        }
        slice_end = MIN(
            task->combination_end,
            pair_end);
        pair_geometry = solver_field_geometry_pair(
            executor->field_geometry,
            pair->field_a,
            pair->field_b);
        if (!pair_geometry) {
            return -1;
        }
        visit_result = visitor(
            executor,
            pair,
            pair_geometry,
            cursor - pair->combination_first,
            slice_end - pair->combination_first,
            opaque);
        if (visit_result) {
            return visit_result;
        }
        cursor = slice_end;
    }
    return cursor == task->combination_end ? 0 : -1;
}

static int solver_ab_builder_visit_pair_range(
    solver_ab_descriptor_planner_t* executor,
    const solver_ab_pair_t* pair,
    const solver_pair_geometry_t* pair_geometry,
    unsigned long long local_first,
    unsigned long long local_end,
    void* opaque) {
    solver_ab_builder_t* builder = opaque;

    if (solver_ab_cancelled(executor)) {
        builder->status->cancelled = TRUE;
        return 1;
    }
    builder->tol2 = pair->tol2;
    builder->rel_field_noise2 =
        pair_geometry->rel_field_noise2;
    builder->rel_field_noise_valid = TRUE;
    return solver_ab_visit_pair_range(
        executor,
        pair,
        pair_geometry,
        local_first,
        local_end,
        solver_ab_builder_visit,
        builder);
}

static int solver_ab_counter_failure(
    solver_t* solver,
    const char* counter_name) {
    logerr(
        "[solver-ab] signed counter boundary reached: %s\n",
        counter_name ? counter_name : "(unknown)");
    if (solver) {
        solver->profile.execution_failed = TRUE;
        solver->quit_now = TRUE;
    }
    return -1;
}

static anbool solver_ab_counter_can_add(
    int current,
    unsigned long long delta) {
    if (current < 0 ||
        delta > (unsigned long long)INT_MAX) {
        return FALSE;
    }
    return delta <=
        (unsigned long long)(INT_MAX - current);
}

int solver_ab_checked_counter_delta(
    solver_t* solver,
    unsigned long long numtries,
    unsigned long long cxdx,
    unsigned long long meanx) {
    if (!solver) {
        return -1;
    }
    /*
     * Preflight every destination before mutating any of them. This preserves
     * an exact reducer prefix even at the representable counter boundary.
     */
    if (!solver_ab_counter_can_add(
            solver->numtries,
            numtries)) {
        return solver_ab_counter_failure(
            solver,
            "numtries");
    }
    if (!solver_ab_counter_can_add(
            solver->num_cxdx_skipped,
            cxdx)) {
        return solver_ab_counter_failure(
            solver,
            "num_cxdx_skipped");
    }
    if (!solver_ab_counter_can_add(
            solver->num_meanx_skipped,
            meanx)) {
        return solver_ab_counter_failure(
            solver,
            "num_meanx_skipped");
    }
    solver->numtries += (int)numtries;
    solver->num_cxdx_skipped += (int)cxdx;
    solver->num_meanx_skipped += (int)meanx;
    return 0;
}

static int solver_ab_reserve_pair_buffer(
    solver_ab_descriptor_planner_t* executor,
    size_t required) {
    solver_ab_pair_t* resized;
    size_t capacity;
    size_t retained_maximum;

    if (!executor ||
        required > SIZE_MAX / sizeof(*resized)) {
        return -1;
    }
    if (!executor->pair_cache_limit_bytes) {
        executor->pair_cache_limit_bytes =
            SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES;
    }
    if (required <= executor->pair_capacity) {
        return 0;
    }
    capacity = executor->pair_capacity ?
        executor->pair_capacity : 64U;
    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            capacity = required;
        } else {
            capacity *= 2U;
        }
    }

    retained_maximum =
        executor->pair_cache_limit_bytes /
        sizeof(*resized);
    if (!executor->pairs_transient &&
        required <= retained_maximum) {
        if (capacity > retained_maximum) {
            capacity = retained_maximum;
        }
        resized = realloc(
            executor->pair_cache,
            capacity * sizeof(*resized));
        if (!resized) {
            return -1;
        }
        executor->pair_cache = resized;
        executor->pair_cache_capacity = capacity;
        executor->pairs = resized;
        executor->pair_capacity = capacity;
        return 0;
    }

    if (capacity > SIZE_MAX / sizeof(*resized)) {
        capacity = required;
    }
    if (executor->pairs_transient) {
        resized = realloc(
            executor->pairs,
            capacity * sizeof(*resized));
    } else {
        resized = malloc(capacity * sizeof(*resized));
        if (resized && required > 1U && executor->pairs) {
            memcpy(
                resized,
                executor->pairs,
                (required - 1U) * sizeof(*resized));
        }
    }
    if (!resized) {
        return -1;
    }
    executor->pairs = resized;
    executor->pair_capacity = capacity;
    executor->pairs_transient = TRUE;
    return 0;
}

static void solver_ab_begin_pair_buffer(
    solver_ab_descriptor_planner_t* executor) {
    executor->pairs = executor->pair_cache;
    executor->pair_capacity =
        executor->pair_cache_capacity;
    executor->pairs_transient = FALSE;
}

int solver_ab_collect_pairs(
    solver_ab_descriptor_planner_t* executor,
    solver_ab_phase_kind_t phase,
    int newpoint,
    int dimquads,
    const solver_field_geometry_t* geometry,
    double min_ab2,
    double max_ab2,
    solver_ab_pair_t** pairs_out,
    size_t* pair_count_out,
    unsigned long long* total_combinations_out) {
    size_t count = 0U;
    unsigned long long total_combinations = 0;
    int field_a;
    int field_b;

    if (!executor || executor->pairs ||
        executor->pair_count != 0U) {
        return -1;
    }
    solver_ab_begin_pair_buffer(executor);
    if (newpoint <= 0) {
        *pairs_out = NULL;
        *pair_count_out = 0U;
        *total_combinations_out = 0U;
        return 0;
    }

    if (phase == SOLVER_AB_PHASE_DIAGONAL) {
        field_b = newpoint;
        for (field_a = 0;
             field_a < newpoint;
             field_a++) {
            const solver_pair_geometry_t* pair_geometry =
                solver_field_geometry_pair(
                    geometry,
                    field_a,
                    field_b);
            unsigned long long combination_count;

            if (!pair_geometry ||
                !pair_geometry->scale_ok ||
                pair_geometry->scale < min_ab2 ||
                pair_geometry->scale > max_ab2) {
                continue;
            }
            combination_count =
                solver_ab_saturating_choose(
                solver_pair_geometry_eligible_before(
                    geometry,
                    pair_geometry,
                    field_a,
                    field_b,
                    newpoint),
                dimquads - 2);
            if (!combination_count) {
                continue;
            }
            if (solver_ab_reserve_pair_buffer(
                    executor, count + 1U)) {
                return -1;
            }
            executor->pairs[count].field_a = field_a;
            executor->pairs[count].field_b = field_b;
            executor->pairs[count].combination_first =
                total_combinations;
            executor->pairs[count].combination_count =
                combination_count;
            executor->pairs[count].tol2 = get_tolerance_for_noise(
                executor->snapshot.codetol,
                pair_geometry->rel_field_noise2,
                executor->snapshot.rel_index_noise2);
            count++;
            total_combinations =
                solver_ab_saturating_add(
                    total_combinations,
                    combination_count);
        }
    } else {
        for (field_a = 0;
             field_a < newpoint;
             field_a++) {
            for (field_b = field_a + 1;
                 field_b < newpoint;
                 field_b++) {
                const solver_pair_geometry_t* pair_geometry =
                    solver_field_geometry_pair(
                        geometry,
                        field_a,
                        field_b);
                unsigned long long combination_count;
                double newpoint_x;
                double newpoint_y;

                if (!pair_geometry ||
                    !pair_geometry->scale_ok ||
                    pair_geometry->scale < min_ab2 ||
                    pair_geometry->scale > max_ab2 ||
                    !solver_pair_geometry_transform(
                        geometry,
                        pair_geometry,
                        field_a,
                        field_b,
                        newpoint,
                        &newpoint_x,
                        &newpoint_y)) {
                    continue;
                }
                combination_count = dimquads > 3 ?
                    solver_ab_saturating_choose(
                        solver_pair_geometry_eligible_before(
                            geometry,
                            pair_geometry,
                            field_a,
                            field_b,
                            newpoint),
                        dimquads - 3) :
                    1U;
                if (!combination_count) {
                    continue;
                }
                if (solver_ab_reserve_pair_buffer(
                        executor, count + 1U)) {
                    return -1;
                }
                executor->pairs[count].field_a = field_a;
                executor->pairs[count].field_b = field_b;
                executor->pairs[count].combination_first =
                    total_combinations;
                executor->pairs[count].combination_count =
                    combination_count;
                executor->pairs[count].tol2 = get_tolerance_for_noise(
                    executor->snapshot.codetol,
                    pair_geometry->rel_field_noise2,
                    executor->snapshot.rel_index_noise2);
                count++;
                total_combinations =
                    solver_ab_saturating_add(
                        total_combinations,
                        combination_count);
            }
        }
    }

    if (!count) {
        *pairs_out = NULL;
        *pair_count_out = 0U;
        *total_combinations_out = 0U;
        return 0;
    }
    *pairs_out = executor->pairs;
    *pair_count_out = count;
    *total_combinations_out = total_combinations;
    return 0;
}

/*
 * Inner compute assistance never starts or waits for a cold page provider.
 * Only an already resident Quad/Star payload may enter prepared verification;
 * every other case remains on the native owner-local mmap path.
 */
anbool solver_payload_candidate_data_fully_resident(
    const solver_t* solver) {
    if (!solver || !solver->index || !solver->index->quads ||
        !solver->index->quads->fb || !solver->index->starkd ||
        !solver->index->starkd->tree ||
        !solver->index->starkd->tree->io ||
        !solver->index->starkd->tree->io_is_fitsbin) {
        return FALSE;
    }
    return fitsbin_payload_is_fully_resident(
               solver->index->quads->fb) &&
        fitsbin_payload_is_fully_resident(
               (const fitsbin_t*)solver->index->starkd->tree->io);
}

static index_shard_helper_task_status_t
solver_verification_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_verification_task_input_t* input = input_bytes;
    verify_prepared_score_t* scores = output_bytes;
    size_t i;

    if (!input || input_size != sizeof(*input) ||
        !input->slots || !input->slot_count || !scores ||
        input->slot_count >
            SIZE_MAX / sizeof(*scores) ||
        output_size !=
            input->slot_count * sizeof(*scores)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    for (i = 0U; i < input->slot_count; i++) {
        const solver_verification_score_slot_t* slot =
            &input->slots[i];

        if (index_shard_worker_stop_requested()) {
            return INDEX_SHARD_HELPER_TASK_STOPPED;
        }
        memset(&scores[i], 0, sizeof(scores[i]));
        if (!slot->prepared ||
            verify_score_prepared_hit(
                slot->prepared,
                &scores[i])) {
            return INDEX_SHARD_HELPER_TASK_ERROR;
        }
    }
    return INDEX_SHARD_HELPER_TASK_OK;
}

const index_shard_helper_ops_t
solver_verification_helper_ops = {
    "verify-context",
    solver_verification_helper_execute
};

static index_shard_helper_retire_status_t
solver_verification_packet_retire(
    const index_shard_helper_task_t* task,
    size_t task_index,
    void* owner_context) {
    const solver_verification_task_input_t* input;
    solver_verification_retire_context_t* context = owner_context;

    (void)task_index;
    if (!task || !context || !task->input ||
        task->input_bytes != sizeof(*input)) {
        return INDEX_SHARD_HELPER_RETIRE_ERROR;
    }
    input = task->input;
    if (!input->slot_count ||
        input->slot_first != context->next_slot ||
        input->slot_count > SIZE_MAX - context->next_slot) {
        return INDEX_SHARD_HELPER_RETIRE_ERROR;
    }
    context->next_slot += input->slot_count;
    return INDEX_SHARD_HELPER_RETIRE_OK;
}

#define SOLVER_VERIFICATION_WAVE_MAX_BYTES \
    (32U * 1024U * 1024U)
size_t solver_verification_wave_memory_budget(void) {
    size_t budget = SOLVER_VERIFICATION_WAVE_MAX_BYTES;

#if defined(_SC_AVPHYS_PAGES)
    long available_pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);

    if (available_pages > 0 && page_size > 0 &&
        (unsigned long)available_pages <=
            SIZE_MAX / (unsigned long)page_size) {
        size_t pressure_budget =
            ((size_t)available_pages * (size_t)page_size) / 64U;

        budget = MIN(budget, pressure_budget);
    }
#endif
    return budget;
}

static void solver_verification_wave_cleanup(
    solver_verification_score_slot_t* slots,
    verify_prepared_score_t* scores,
    size_t slot_count,
    solver_verification_candidate_runtime_t* runtime,
    solver_verification_packet_t* packet) {
    size_t i;

    if (slots) {
        for (i = 0U; i < slot_count; i++) {
            if (scores) {
                verify_destroy_prepared_score(&scores[i]);
            }
            verify_destroy_prepared_hit(slots[i].prepared);
            slots[i].prepared = NULL;
        }
    }
    free(slots);
    free(scores);
    free(runtime);
    solver_verification_packet_free(packet);
}

/*
 * Prepare index-backed candidate inputs on the owner, score immutable
 * verification contexts in coarse helper tasks, then retire every action in
 * the original candidate order. No helper receives an index, mapping, solver,
 * callback, or reducer pointer.
 *
 * Return zero before publication/state mutation to request the native path.
 * Return one after successful handling, stop propagation, or a hard failure.
 */
int solver_ab_try_verification_wave(
    kdtree_qres_t* result,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    int candidate_first,
    int candidate_end) {
    solver_verification_packet_t packet;
    solver_ab_snapshot_t snapshot;
    solver_verification_candidate_runtime_t* runtime = NULL;
    solver_verification_score_slot_t* slots = NULL;
    verify_prepared_score_t* scores = NULL;
    index_shard_helper_task_t tasks[INDEX_SHARD_HELPER_MAX_TASKS];
    solver_verification_task_input_t
        inputs[INDEX_SHARD_HELPER_MAX_TASKS];
    index_shard_helper_run_stats_t run_stats;
    index_shard_helper_run_status_t run_status;
    solver_verification_retire_context_t retire_context;
    size_t available = 0U;
    size_t candidate_count;
    size_t verify_count = 0U;
    size_t peak_bytes = 0U;
    size_t peak_budget;
    size_t task_count;
    size_t task_index;
    size_t slot_cursor = 0U;
    double verify_wall_start;
    double verify_wall_seconds;
    int candidate_index;
    int handled = 0;

    if (!result || !field_xy || !fieldstars || !solver ||
        candidate_first < 0 ||
        candidate_end <= candidate_first ||
        candidate_end > result->nres ||
        candidate_end - candidate_first >
            (int)SOLVER_VERIFICATION_WINDOW_CANDIDATES ||
        verify_datalog_enabled() ||
        !solver_payload_candidate_data_fully_resident(solver) ||
        !index_shard_worker_context_active()) {
        return 0;
    }
    candidate_count = (size_t)(candidate_end - candidate_first);
    available = index_shard_helper_available_workers();
    if (!available) {
        return 0;
    }

    memset(&packet, 0, sizeof(packet));
    memset(&snapshot, 0, sizeof(snapshot));
    memset(tasks, 0, sizeof(tasks));
    memset(inputs, 0, sizeof(inputs));
    memset(&run_stats, 0, sizeof(run_stats));
    memset(&retire_context, 0, sizeof(retire_context));

    {
        solver_packet_reserve_result_t reserve_status;

        reserve_status = solver_verification_packet_reserve(
            &packet, candidate_count);
        if (reserve_status != SOLVER_PACKET_RESERVE_OK) {
            if (packet.allocation_failed) {
                solver->profile.allocation_failures++;
            }
            goto cleanup;
        }
    }
    runtime = calloc(candidate_count, sizeof(*runtime));
    slots = calloc(candidate_count, sizeof(*slots));
    scores = calloc(candidate_count, sizeof(*scores));
    if (!runtime || !slots || !scores) {
        goto cleanup;
    }

    snapshot.index = solver->index;
    snapshot.fieldxy = solver->fieldxy;
    snapshot.use_radec = solver->use_radec;
    snapshot.cx_less_than_dx = solver->index->cx_less_than_dx;
    snapshot.meanx_less_than_half =
        solver->index->meanx_less_than_half;
    snapshot.parity = solver->parity;
    memcpy(snapshot.centerxyz,
           solver->centerxyz,
           sizeof(snapshot.centerxyz));
    snapshot.r2 = solver->r2;
    snapshot.abscale_low = solver->abscale_low;
    snapshot.abscale_high = solver->abscale_high;
    snapshot.funits_lower = solver->funits_lower;
    snapshot.funits_upper = solver->funits_upper;
    snapshot.cxdx_margin = solver->cxdx_margin;
    snapshot.codetol = solver->codetol;
    snapshot.rel_index_noise2 = solver->rel_index_noise2;

    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];

        if (solver_poll_worker_stop(solver)) {
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
        if (solver_ab_candidate_prepare(
                candidate,
                result,
                candidate_index,
                field_xy,
                fieldstars,
                dimquads,
                &snapshot,
                current_parity)) {
            goto cleanup;
        }
        if (candidate->action ==
            SOLVER_AB_CANDIDATE_VERIFY) {
            verify_count++;
        }
    }
    if (verify_count < 2U) {
        goto cleanup;
    }
    available = index_shard_helper_prepare_reserve();
    if (!available) {
        goto cleanup;
    }
    verify_wall_start = monotonic_seconds();
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[packet_index];
        int i;

        if (candidate->action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        set_matchobj_template(solver, &candidate_runtime->match);
        memcpy(&candidate_runtime->match.wcstan,
               &candidate->wcs,
               sizeof(tan_t));
        candidate_runtime->match.wcs_valid = TRUE;
        candidate_runtime->match.code_err = candidate->code_err;
        candidate_runtime->match.scale = candidate->scale;
        candidate_runtime->match.parity = candidate->parity;
        candidate_runtime->match.quad_npeers =
            candidate->quad_npeers;
        candidate_runtime->match.timeused = solver->timeused;
        candidate_runtime->match.quadno = candidate->quadno;
        candidate_runtime->match.dimquads = dimquads;
        for (i = 0; i < dimquads; i++) {
            candidate_runtime->match.star[i] = candidate->star[i];
            candidate_runtime->match.field[i] = candidate->field[i];
            candidate_runtime->match.ids[i] = 0;
        }
        memcpy(candidate_runtime->match.quadpix,
               candidate->quadpix,
               (size_t)2 * (size_t)dimquads * sizeof(double));
        memcpy(candidate_runtime->match.quadxyz,
               candidate->quadxyz,
               (size_t)3 * (size_t)dimquads * sizeof(double));
        set_center_and_radius(
            solver,
            &candidate_runtime->match,
            &candidate_runtime->match.wcstan,
            NULL);
        candidate_runtime->match_distance_in_pixels2 =
            solver_prepare_hit_for_verify(
                solver,
                &candidate_runtime->match,
                &candidate_runtime->logaccept);
    }
    peak_budget = solver_verification_wave_memory_budget();
    slot_cursor = 0U;
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];
        solver_verification_candidate_runtime_t* candidate_runtime =
            &runtime[packet_index];
        solver_verification_score_slot_t* slot;
        size_t context_peak;

        if (candidate->action != SOLVER_AB_CANDIDATE_VERIFY) {
            continue;
        }
        slot = &slots[slot_cursor++];
        if (verify_prepare_hit(
                solver->index->starkd,
                solver->index->cutnside,
                &candidate_runtime->match,
                NULL,
                solver->vf,
                candidate_runtime->match_distance_in_pixels2,
                solver->distractor_ratio,
                solver->field_maxx,
                solver->field_maxy,
                solver->logratio_bail_threshold,
                candidate_runtime->logaccept,
                solver->logratio_stoplooking,
                solver->distance_from_quad_bonus,
                FALSE,
                &slot->prepared)) {
            goto cleanup;
        }
        if (solver_poll_worker_stop(solver)) {
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
        context_peak =
            verify_prepared_hit_peak_bytes(slot->prepared);
        if (context_peak == SIZE_MAX ||
            peak_bytes > peak_budget ||
            context_peak > peak_budget - peak_bytes) {
            goto cleanup;
        }
        peak_bytes += context_peak;
    }
    if (slot_cursor != verify_count) {
        goto cleanup;
    }

    if (verify_count) {
        task_count = MIN(
            verify_count,
            MIN(available + 1U,
                (size_t)INDEX_SHARD_HELPER_MAX_TASKS));
        slot_cursor = 0U;
        for (task_index = 0U;
             task_index < task_count;
             task_index++) {
            size_t remaining_slots = verify_count - slot_cursor;
            size_t remaining_tasks = task_count - task_index;
            size_t count =
                (remaining_slots + remaining_tasks - 1U) /
                remaining_tasks;
            size_t i;
            unsigned long long work_units = 0U;

            inputs[task_index].slots = &slots[slot_cursor];
            inputs[task_index].slot_first = slot_cursor;
            inputs[task_index].slot_count = count;
            for (i = 0U; i < count; i++) {
                unsigned long long work =
                    verify_prepared_hit_work_units(
                        slots[slot_cursor + i].prepared);

                if (!work) {
                    work = 1U;
                }
                if (ULLONG_MAX - work_units < work) {
                    work_units = ULLONG_MAX;
                } else {
                    work_units += work;
                }
            }
            tasks[task_index].input = &inputs[task_index];
            tasks[task_index].input_bytes = sizeof(inputs[task_index]);
            tasks[task_index].output = &scores[slot_cursor];
            tasks[task_index].output_bytes =
                count * sizeof(*scores);
            tasks[task_index].work_units = work_units;
            slot_cursor += count;
        }
        if (slot_cursor != verify_count) {
            goto cleanup;
        }

        run_status = index_shard_helper_run_ordered(
            &solver_verification_helper_ops,
            tasks,
            task_count,
            solver_verification_packet_retire,
            &retire_context,
            &run_stats);
        if (run_status == INDEX_SHARD_HELPER_UNAVAILABLE ||
            run_status == INDEX_SHARD_HELPER_TASK_FAILED) {
            goto cleanup;
        }
        if (run_status == INDEX_SHARD_HELPER_STOPPED) {
            (void)solver_poll_worker_stop(solver);
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
        if (run_status != INDEX_SHARD_HELPER_OK) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
        if (retire_context.next_slot != verify_count) {
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
    } else {
        task_count = 0U;
        run_status = INDEX_SHARD_HELPER_OK;
    }
    verify_wall_seconds =
        monotonic_seconds() - verify_wall_start;
    if (solver->profile.detailed) {
        solver->profile.verify_wall_seconds +=
            verify_wall_seconds;
    }
    if (run_stats.foreign_tasks) {
        solver->profile.parallel_batches++;
        solver->profile.parallel_batches_observed++;
        solver->profile.ab_helper_tasks =
            solver_ab_saturating_add(
                solver->profile.ab_helper_tasks,
                run_stats.foreign_tasks);
        solver->profile.max_parallel_ranges = MAX(
            solver->profile.max_parallel_ranges,
            run_stats.max_concurrent_tasks);
    }

    slot_cursor = 0U;
    for (candidate_index = candidate_first;
         candidate_index < candidate_end;
         candidate_index++) {
        int packet_index = candidate_index - candidate_first;
        solver_ab_candidate_t* candidate =
            &packet.candidates[packet_index];

        if (solver_poll_worker_stop(solver)) {
            handled = 1;
            goto cleanup;
        }
        if (solver->nummatches < 0 ||
            solver->nummatches == INT_MAX) {
            (void)solver_ab_counter_failure(
                solver, "nummatches");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_RADEC_SKIP &&
            (solver->num_radec_skipped < 0 ||
             solver->num_radec_skipped == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "num_radec_skipped");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_ABSCALE_SKIP &&
            (solver->num_abscale_skipped < 0 ||
             solver->num_abscale_skipped == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "num_abscale_skipped");
            handled = 1;
            goto cleanup;
        }
        if (candidate->action ==
                SOLVER_AB_CANDIDATE_VERIFY &&
            (solver->numscaleok < 0 ||
             solver->numscaleok == INT_MAX ||
             solver->num_verified < 0 ||
             solver->num_verified == INT_MAX)) {
            (void)solver_ab_counter_failure(
                solver, "numscaleok/num_verified");
            handled = 1;
            goto cleanup;
        }
        solver->nummatches++;
        solver_record_candidate_order(
            solver,
            candidate->action,
            candidate->quadno,
            candidate->code_err);
        switch (candidate->action) {
        case SOLVER_AB_CANDIDATE_RADEC_SKIP:
            solver->num_radec_skipped++;
            break;

        case SOLVER_AB_CANDIDATE_ABSCALE_SKIP:
            solver->num_abscale_skipped++;
            break;

        case SOLVER_AB_CANDIDATE_BAD_QUAD:
            logverb("bad quad at %s:%i\n", __FILE__, __LINE__);
            break;

        case SOLVER_AB_CANDIDATE_SCALE_SKIP:
            break;

        case SOLVER_AB_CANDIDATE_VERIFY:
        {
            solver_verification_candidate_runtime_t* candidate_runtime =
                &runtime[packet_index];
            solver_verification_score_slot_t* slot =
                &slots[slot_cursor++];

            solver->numscaleok++;
            candidate_runtime->match.quads_tried = quads_tried;
            candidate_runtime->match.quads_matched =
                solver->nummatches;
            candidate_runtime->match.quads_scaleok =
                solver->numscaleok;
            if (verify_finish_prepared_hit(
                    slot->prepared,
                    &scores[slot_cursor - 1U],
                    &candidate_runtime->match)) {
                verify_destroy_prepared_score(
                    &scores[slot_cursor - 1U]);
                verify_destroy_prepared_hit(slot->prepared);
                slot->prepared = NULL;
                if (solver_handle_hit(
                        solver,
                        &candidate_runtime->match,
                        NULL,
                        FALSE)) {
                    solver->quit_now = TRUE;
                }
                if (unlikely(solver->quit_now)) {
                    handled = 1;
                    goto cleanup;
                }
                break;
            }
            solver->profile.verify_calls++;
            if (solver_handle_hit_after_verify(
                    solver,
                    &candidate_runtime->match,
                    NULL,
                    FALSE,
                    candidate_runtime->match_distance_in_pixels2)) {
                solver->quit_now = TRUE;
            }
            verify_destroy_prepared_hit(slot->prepared);
            slot->prepared = NULL;
            if (unlikely(solver->quit_now)) {
                handled = 1;
                goto cleanup;
            }
            break;
        }

        default:
            solver->profile.execution_failed = TRUE;
            solver->quit_now = TRUE;
            handled = 1;
            goto cleanup;
        }
    }
    handled = 1;

cleanup:
    index_shard_helper_prepare_cancel();
    solver_verification_wave_cleanup(
        slots,
        scores,
        verify_count,
        runtime,
        &packet);
    return handled;
}

static pthread_key_t solver_ab_descriptor_workspace_key;
static pthread_once_t solver_ab_descriptor_workspace_once =
    PTHREAD_ONCE_INIT;
static int solver_ab_descriptor_workspace_status = EAGAIN;

void solver_ab_descriptor_release_pairs(
    solver_ab_descriptor_workspace_t* workspace) {
    solver_ab_descriptor_planner_t* planner;

    if (!workspace) {
        return;
    }
    planner = &workspace->planner;
    if (planner->pairs_transient) {
        free(planner->pairs);
    }
    planner->pairs = NULL;
    planner->pair_count = 0U;
    planner->pair_capacity = 0U;
    planner->pairs_transient = FALSE;
}

static void solver_ab_descriptor_workspace_destroy(void* opaque) {
    solver_ab_descriptor_workspace_t* workspace = opaque;

    if (!workspace) {
        return;
    }
    solver_ab_descriptor_release_pairs(workspace);
    free(workspace->planner.pair_cache);
    free(workspace->outputs);
    free(workspace);
}

static void solver_ab_descriptor_workspace_make_key(void) {
    solver_ab_descriptor_workspace_status =
        pthread_key_create(
            &solver_ab_descriptor_workspace_key,
            solver_ab_descriptor_workspace_destroy);
}

int solver_ab_descriptor_workspace_reserve_outputs(
    solver_ab_descriptor_workspace_t* workspace,
    size_t output_count) {
    solver_ab_descriptor_output_t* outputs;

    if (!workspace || !output_count ||
        output_count > SIZE_MAX / sizeof(*outputs)) {
        return -1;
    }
    if (output_count <= workspace->output_capacity) {
        return 0;
    }
    outputs = calloc(output_count, sizeof(*outputs));
    if (!outputs) {
        return -1;
    }
    free(workspace->outputs);
    workspace->outputs = outputs;
    workspace->output_capacity = output_count;
    return 0;
}

solver_ab_descriptor_workspace_t*
solver_ab_descriptor_workspace_get(void) {
    solver_ab_descriptor_workspace_t* workspace;

    if (pthread_once(
            &solver_ab_descriptor_workspace_once,
            solver_ab_descriptor_workspace_make_key) ||
        solver_ab_descriptor_workspace_status) {
        return NULL;
    }
    workspace = pthread_getspecific(
        solver_ab_descriptor_workspace_key);
    if (workspace) {
        return workspace;
    }
    workspace = calloc(1, sizeof(*workspace));
    if (!workspace) {
        return NULL;
    }
    workspace->planner.pair_cache_limit_bytes =
        SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES;
    if (pthread_setspecific(
            solver_ab_descriptor_workspace_key,
            workspace)) {
        solver_ab_descriptor_workspace_destroy(workspace);
        return NULL;
    }
    return workspace;
}

index_shard_helper_task_status_t
solver_ab_descriptor_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size) {
    const solver_ab_descriptor_task_input_t* input = input_bytes;
    solver_ab_descriptor_output_t* output = output_bytes;
    solver_ab_descriptor_planner_t executor;
    solver_ab_builder_t builder;
    solver_descriptor_status_t status;
    solver_ab_task_t task;
    int eligible[SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS];
    int visit_result;

    if (!input ||
        input_size != sizeof(*input) ||
        !output ||
        output_size != sizeof(*output) ||
        !input->field_geometry ||
        !input->pairs || !input->pair_count ||
        input->combination_first >= input->combination_end ||
        input->newpoint < 0 ||
        input->newpoint >=
            SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS ||
        input->dimquads < NBACK ||
        input->dimquads > DQMAX ||
        (input->phase != SOLVER_AB_PHASE_DIAGONAL &&
         input->phase != SOLVER_AB_PHASE_OFF_DIAGONAL) ||
        (input->parity != PARITY_NORMAL &&
         input->parity != PARITY_FLIP &&
         input->parity != PARITY_BOTH)) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }

    memset(&executor, 0, sizeof(executor));
    memset(&builder, 0, sizeof(builder));
    memset(&status, 0, sizeof(status));
    memset(&task, 0, sizeof(task));
    output->descriptor_count = 0U;
    output->trailing_numtries = 0U;
    output->trailing_cxdx = 0U;
    output->trailing_meanx = 0U;
    output->has_final_rel_field_noise2 = FALSE;

    executor.field_geometry = input->field_geometry;
    executor.newpoint = input->newpoint;
    executor.dimquads = input->dimquads;
    executor.phase = input->phase;
    executor.pairs = (solver_ab_pair_t*)input->pairs;
    executor.pair_count = input->pair_count;
    executor.combination_eligible = eligible;
    executor.combination_eligible_capacity =
        SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS;
    executor.snapshot.parity = input->parity;
    executor.snapshot.cx_less_than_dx =
        input->cx_less_than_dx;
    executor.snapshot.meanx_less_than_half =
        input->meanx_less_than_half;
    executor.snapshot.cxdx_margin = input->cxdx_margin;

    builder.planner = &executor;
    builder.status = &status;
    builder.descriptor_output = output;
    task.combination_first = input->combination_first;
    task.combination_end = input->combination_end;
    visit_result = solver_ab_walk_task_ranges(
        &executor,
        &task,
        solver_ab_builder_visit_pair_range,
        &builder);
    if (status.cancelled ||
        index_shard_worker_stop_requested()) {
        return INDEX_SHARD_HELPER_TASK_STOPPED;
    }
    if (visit_result < 0 ||
        builder.fatal_error ||
        status.evaluation_failed) {
        return INDEX_SHARD_HELPER_TASK_ERROR;
    }
    output->trailing_numtries = builder.pending_numtries;
    output->trailing_cxdx = builder.pending_cxdx;
    output->trailing_meanx = builder.pending_meanx;
    output->final_rel_field_noise2 =
        builder.rel_field_noise2;
    output->has_final_rel_field_noise2 =
        builder.rel_field_noise_valid;
    return INDEX_SHARD_HELPER_TASK_OK;
}

size_t solver_ab_descriptor_expansion(
    int dimquads,
    int parity) {
    size_t expansion = 2U;
    int internal_stars = dimquads - NBACK;
    int i;

    if (dimquads < NBACK || dimquads > DQMAX ||
        (parity != PARITY_NORMAL &&
         parity != PARITY_FLIP &&
         parity != PARITY_BOTH)) {
        return 0U;
    }
    if (parity == PARITY_BOTH) {
        expansion *= 2U;
    }
    for (i = 2; i <= internal_stars; i++) {
        expansion *= (size_t)i;
    }
    return expansion;
}

size_t solver_ab_descriptor_partition_count(
    unsigned long long wave_combinations,
    size_t expansion,
    size_t max_task_combinations,
    size_t participants) {
    unsigned long long expanded_work;
    unsigned long long capacity_tasks;
    unsigned long long work_tasks;
    unsigned long long task_count;

    if (!wave_combinations || !expansion ||
        !max_task_combinations || !participants ||
        wave_combinations > ULLONG_MAX /
            (unsigned long long)expansion) {
        return 0U;
    }
    expanded_work = wave_combinations *
        (unsigned long long)expansion;
    capacity_tasks = wave_combinations /
        (unsigned long long)max_task_combinations;
    if (wave_combinations %
            (unsigned long long)max_task_combinations) {
        capacity_tasks++;
    }
    if (!capacity_tasks ||
        capacity_tasks > (unsigned long long)participants) {
        return 0U;
    }
    work_tasks = expanded_work /
        SOLVER_AB_DESCRIPTOR_TASK_TARGET_HYPOTHESES;
    if (expanded_work %
            SOLVER_AB_DESCRIPTOR_TASK_TARGET_HYPOTHESES) {
        work_tasks++;
    }
    work_tasks = MIN(
        work_tasks,
        (unsigned long long)participants);
    work_tasks = MIN(work_tasks, wave_combinations);
    task_count = MAX(capacity_tasks, work_tasks);
    if (!task_count ||
        task_count > (unsigned long long)participants ||
        task_count > (unsigned long long)SIZE_MAX) {
        return 0U;
    }
    return (size_t)task_count;
}
