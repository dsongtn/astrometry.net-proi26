/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "index_shard_internal.h"

#include "astrometry/index_shard.h"
#include "astrometry/ioutils.h"
#include "astrometry/log.h"
#include "astrometry/matchfile.h"
#include "astrometry/matchobj.h"
#include "astrometry/onefield.h"
#include "astrometry/sip_qfits.h"
#include "astrometry/solver.h"

static int parse_positive_int(
    const char* text,
    const char* label,
    int* value_out) {
    char* end = NULL;
    long value;

    if (!text || !text[0] || !value_out) {
        return -1;
    }
    errno = 0;
    value = strtol(text, &end, 10);
    if (errno == ERANGE ||
        end == text ||
        *end != '\0' ||
        value < 1 ||
        value > INT_MAX) {
        fprintf(stderr, "%s must be a positive integer\n", label);
        return -1;
    }
    *value_out = (int)value;
    return 0;
}

static int parse_positive_double_env(
    const char* name,
    double* value_out) {
    const char* text = getenv(name);
    char* end = NULL;
    double value;

    if (!text || !text[0]) {
        return 0;
    }
    errno = 0;
    value = strtod(text, &end);
    if (errno == ERANGE ||
        end == text ||
        *end != '\0' ||
        !isfinite(value) ||
        value <= 0.0) {
        fprintf(stderr, "%s must be a positive finite number\n", name);
        return -1;
    }
    *value_out = value;
    return 1;
}

static int touch_before_second_pass(const char* env_name) {
    const char* path = getenv(env_name);

    if (!path || !path[0]) {
        return 0;
    }
    if (utimensat(AT_FDCWD, path, NULL, 0)) {
        fprintf(
            stderr,
            "failed to update %s source identity for %s: %s\n",
            env_name,
            path,
            strerror(errno));
        return -1;
    }
    return 0;
}

static void signature_hash_bytes(
    uint64_t* hash,
    const void* data,
    size_t length) {
    const unsigned char* bytes = data;
    size_t i;

    for (i = 0U; i < length; i++) {
        *hash ^= (uint64_t)bytes[i];
        *hash *= UINT64_C(1099511628211);
    }
}

typedef struct result_signature {
    uint64_t hash;
    int solution_count;
    anbool best_valid;
    MatchObj best;
} result_signature_t;

/*
 * Hash every persisted scientific value in one match row. Traversal counters
 * are added separately to the exact-result signature because they reset at a
 * configured pass boundary and therefore are not part of continuous-versus-
 * split match-set identity. TIMEUSED is deliberately excluded from both
 * hashes because it is measurement data, not a scientific result.
 */
static void signature_hash_science(
    uint64_t* hash,
    const MatchObj* match) {
    unsigned int logodds_bits;
    unsigned int worstlogodds_bits;

    memcpy(
        &logodds_bits,
        &match->logodds,
        sizeof(logodds_bits));
    memcpy(
        &worstlogodds_bits,
        &match->worstlogodds,
        sizeof(worstlogodds_bits));
    signature_hash_bytes(hash, &match->quadno, sizeof(match->quadno));
    signature_hash_bytes(hash, &match->dimquads, sizeof(match->dimquads));
    signature_hash_bytes(hash, match->star, sizeof(match->star));
    signature_hash_bytes(hash, match->field, sizeof(match->field));
    signature_hash_bytes(hash, match->ids, sizeof(match->ids));
    signature_hash_bytes(hash, &match->code_err, sizeof(match->code_err));
    signature_hash_bytes(hash, match->quadpix, sizeof(match->quadpix));
    signature_hash_bytes(
        hash,
        match->quadpix_orig,
        sizeof(match->quadpix_orig));
    signature_hash_bytes(hash, match->quadxyz, sizeof(match->quadxyz));
    signature_hash_bytes(hash, match->center, sizeof(match->center));
    signature_hash_bytes(
        hash,
        &match->radius_deg,
        sizeof(match->radius_deg));
    signature_hash_bytes(hash, &match->nmatch, sizeof(match->nmatch));
    signature_hash_bytes(
        hash,
        &match->ndistractor,
        sizeof(match->ndistractor));
    signature_hash_bytes(
        hash,
        &match->nconflict,
        sizeof(match->nconflict));
    signature_hash_bytes(hash, &match->nfield, sizeof(match->nfield));
    signature_hash_bytes(hash, &match->nindex, sizeof(match->nindex));
    signature_hash_bytes(
        hash,
        match->wcstan.crval,
        sizeof(match->wcstan.crval));
    signature_hash_bytes(
        hash,
        match->wcstan.crpix,
        sizeof(match->wcstan.crpix));
    signature_hash_bytes(hash, match->wcstan.cd, sizeof(match->wcstan.cd));
    signature_hash_bytes(
        hash,
        &match->wcs_valid,
        sizeof(match->wcs_valid));
    signature_hash_bytes(hash, &match->fieldnum, sizeof(match->fieldnum));
    signature_hash_bytes(hash, &match->fieldfile, sizeof(match->fieldfile));
    signature_hash_bytes(hash, &match->indexid, sizeof(match->indexid));
    signature_hash_bytes(hash, &match->healpix, sizeof(match->healpix));
    signature_hash_bytes(hash, &match->hpnside, sizeof(match->hpnside));
    signature_hash_bytes(
        hash,
        match->fieldname,
        sizeof(match->fieldname) - 1U);
    signature_hash_bytes(hash, &match->parity, sizeof(match->parity));
    signature_hash_bytes(
        hash,
        &match->quad_npeers,
        sizeof(match->quad_npeers));
    signature_hash_bytes(hash, &logodds_bits, sizeof(logodds_bits));
    signature_hash_bytes(
        hash,
        &worstlogodds_bits,
        sizeof(worstlogodds_bits));
}

static void signature_hash_match(
    uint64_t* hash,
    const MatchObj* match) {
    signature_hash_science(hash, match);
    signature_hash_bytes(hash, &match->nagree, sizeof(match->nagree));
    signature_hash_bytes(
        hash,
        &match->quads_tried,
        sizeof(match->quads_tried));
    signature_hash_bytes(
        hash,
        &match->quads_matched,
        sizeof(match->quads_matched));
    signature_hash_bytes(
        hash,
        &match->quads_scaleok,
        sizeof(match->quads_scaleok));
    signature_hash_bytes(
        hash,
        &match->nverified,
        sizeof(match->nverified));
}

static uint64_t scientific_match_digest(
    const MatchObj* match) {
    uint64_t hash = UINT64_C(1469598103934665603);

    signature_hash_science(&hash, match);
    return hash;
}

/*
 * onefield_run() intentionally releases its in-memory MatchObj list before
 * returning. Read the match artifact written before that cleanup so this
 * signature covers reducer-owned scientific results rather than an empty
 * post-run container.
 */
static int read_result_signature(
    const char* match_path,
    result_signature_t* result) {
    matchfile* matches;
    int row;

    if (!match_path || !result) {
        return -1;
    }
    memset(result, 0, sizeof(*result));
    result->hash = UINT64_C(1469598103934665603);

    if (!file_exists(match_path)) {
        return 0;
    }
    matches = matchfile_open(match_path);
    if (!matches) {
        return -1;
    }
    result->solution_count = matchfile_count(matches);
    if (result->solution_count < 0) {
        matchfile_close(matches);
        return -1;
    }

    for (row = 0; row < result->solution_count; row++) {
        const MatchObj* match = matchfile_read_match(matches);

        if (!match) {
            matchfile_close(matches);
            return -1;
        }
        signature_hash_match(&result->hash, match);
        if (!result->best_valid ||
            match->logodds > result->best.logodds) {
            result->best = *match;
            result->best_valid = TRUE;
        }
    }
    if (matchfile_close(matches)) {
        return -1;
    }
    return 0;
}

static int print_match_records(
    const char* match_path,
    const char* mode,
    int pass_number,
    int workers,
    int first_object,
    int last_object) {
    matchfile* matches;
    int solution_count;
    int row;

    if (!file_exists(match_path)) {
        return 0;
    }
    matches = matchfile_open(match_path);
    if (!matches) {
        return -1;
    }
    solution_count = matchfile_count(matches);
    if (solution_count < 0) {
        matchfile_close(matches);
        return -1;
    }
    for (row = 0; row < solution_count; row++) {
        const MatchObj* match = matchfile_read_match(matches);

        if (!match) {
            matchfile_close(matches);
            return -1;
        }
        printf(
            "AB_MATCH mode=%s pass=%i workers=%i "
            "first_object=%i last_object=%i "
            "digest=%016llx indexid=%i healpix=%i hpnside=%i "
            "parity=%i quad=%u max_field_object=%i\n",
            mode,
            pass_number,
            workers,
            first_object,
            last_object,
            (unsigned long long)scientific_match_digest(match),
            match->indexid,
            match->healpix,
            match->hpnside,
            match->parity,
            match->quadno,
            match->objs_tried);
    }
    return matchfile_close(matches);
}

static void print_result_record(
    const onefield_t* bp,
    const result_signature_t* result,
    const char* mode,
    int pass_number,
    int workers,
    int first_object,
    int last_object,
    int index_count) {
    const MatchObj* best =
        result->best_valid ? &result->best : NULL;

    printf(
        "AB_RESULT mode=%s pass=%i workers=%i "
        "first_object=%i last_object=%i "
        "indexes=%i solutions=%i cancelled=%i wall_limit=%i "
        "cpu_limit=%i failed=%i signature=%016llx",
        mode,
        pass_number,
        workers,
        first_object,
        last_object,
        index_count,
        result->solution_count,
        bp->cancelled,
        bp->hit_total_timelimit,
        bp->hit_total_cpulimit,
        bp->solver_failed,
        (unsigned long long)result->hash);
    if (best) {
        unsigned int logodds_bits;

        memcpy(
            &logodds_bits,
            &best->logodds,
            sizeof(logodds_bits));
        printf(
            " best_indexid=%i best_healpix=%i best_hpnside=%i "
            "best_parity=%i best_quad=%u best_max_field_object=%i "
            "best_logodds_bits=%08x",
            best->indexid,
            best->healpix,
            best->hpnside,
            best->parity,
            best->quadno,
            best->objs_tried,
            logodds_bits);
    }
    printf("\n");
}

static int print_wcs_record(
    const char* wcs_path,
    const char* mode,
    int pass_number,
    int workers,
    int first_object,
    int last_object,
    anbool allow_missing) {
    tan_t wcs;

    if (!tan_read_header_file(wcs_path, &wcs)) {
        if (allow_missing) {
            printf(
                "AB_WCS mode=%s pass=%i workers=%i "
                "first_object=%i last_object=%i none=1\n",
                mode,
                pass_number,
                workers,
                first_object,
                last_object);
            return 0;
        }
        fprintf(stderr, "failed to read integration WCS output\n");
        return -1;
    }
    printf(
        "AB_WCS mode=%s pass=%i workers=%i "
        "first_object=%i last_object=%i "
        "crval=%.17g,%.17g crpix=%.17g,%.17g "
        "cd=%.17g,%.17g,%.17g,%.17g\n",
        mode,
        pass_number,
        workers,
        first_object,
        last_object,
        wcs.crval[0],
        wcs.crval[1],
        wcs.crpix[0],
        wcs.crpix[1],
        wcs.cd[0][0],
        wcs.cd[0][1],
        wcs.cd[1][0],
        wcs.cd[1][1]);
    return 0;
}

typedef struct canonical_index_order_test {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    index_t indexes[2];
    int later_completed;
    int completion_count;
    int completion_order[2];
    int analyzed[2];
    int merge_count;
    int merged_order;
    int report_count;
    size_t reported_order;
    int loser_freed;
    int failed;
} canonical_index_order_test_t;

static canonical_index_order_test_t canonical_index_order_test;

static index_t* canonical_index_order_get_index(
    onefield_t* bp,
    size_t index_order) {
    (void)bp;

    if (index_order >= 2U) {
        return NULL;
    }
    return &canonical_index_order_test.indexes[index_order];
}

static int canonical_index_order_done_with_index(
    onefield_t* bp,
    size_t index_order,
    index_t* index) {
    (void)bp;
    (void)index_order;
    (void)index;
    return 0;
}

static int canonical_index_order_report_solution(
    onefield_t* bp,
    size_t index_order,
    int fieldnum,
    double best_logodds) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;

    (void)bp;
    (void)best_logodds;

    pthread_mutex_lock(&test->mutex);
    test->report_count++;
    test->reported_order = index_order;
    if (fieldnum != 1) {
        test->failed = TRUE;
    }
    pthread_mutex_unlock(&test->mutex);
    return 0;
}

static int canonical_index_order_prepare_local(
    onefield_t* local_bp,
    onefield_t* master_bp,
    const solver_t* base_sp) {
    (void)master_bp;
    (void)base_sp;

    memset(local_bp, 0, sizeof(*local_bp));
    local_bp->fieldnum = -1;
    return 0;
}

static void canonical_index_order_reset_local(
    onefield_t* local_bp,
    bl* local_solutions) {
    local_bp->solutions = local_solutions;
    local_bp->single_field_solved = FALSE;
    local_bp->solver_failed = FALSE;
    local_bp->fieldnum = -1;
    memset(
        &local_bp->solver.profile,
        0,
        sizeof(local_bp->solver.profile));
}

/*
 * Worker cleanup follows index_shard_finish_outer_claim(), so observing the
 * later lane here proves its result slot was published before the earlier
 * lane is allowed to return from solve_one_index().
 */
static void canonical_index_order_cleanup_local(
    onefield_t* local_bp) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;
    int index_order = local_bp->fieldnum;

    pthread_mutex_lock(&test->mutex);
    if (index_order < 0 ||
        index_order >= 2 ||
        test->completion_count >= 2) {
        test->failed = TRUE;
    } else {
        test->completion_order[
            test->completion_count++] = index_order;
    }
    if (index_order == 1) {
        test->later_completed = TRUE;
        pthread_cond_broadcast(&test->condition);
    }
    pthread_mutex_unlock(&test->mutex);

    local_bp->solutions = NULL;
}

static int canonical_index_order_solve_one(
    onefield_t* local_bp,
    index_t* index) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;
    MatchObj match;
    int index_order;
    int wait_status = 0;

    if (!local_bp || !local_bp->solutions || !index) {
        return -1;
    }

    index_order = index->indexid;
    if (index_order < 0 || index_order >= 2) {
        return -1;
    }
    local_bp->fieldnum = index_order;

    if (index_order == 0) {
        pthread_mutex_lock(&test->mutex);
        while (!test->later_completed && !wait_status) {
            wait_status =
                pthread_cond_wait(
                    &test->condition,
                    &test->mutex);
        }
        if (wait_status) {
            test->failed = TRUE;
        }
        pthread_mutex_unlock(&test->mutex);
        if (wait_status) {
            return -1;
        }
    }

    memset(&match, 0, sizeof(match));
    match.fieldnum = 1;
    match.indexid = index_order;
    match.logodds = 100.0;
    if (!bl_append(local_bp->solutions, &match)) {
        return -1;
    }
    return 0;
}

static anbool canonical_index_order_analyze(
    onefield_t* master_bp,
    bl* solutions,
    double* best_logodds,
    int* best_fieldnum) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;
    MatchObj* match;
    int index_order;

    (void)master_bp;

    if (!solutions || bl_size(solutions) != 1U) {
        return FALSE;
    }
    match = bl_access(solutions, 0);
    index_order = match->indexid;
    if (index_order < 0 || index_order >= 2) {
        return FALSE;
    }

    pthread_mutex_lock(&test->mutex);
    test->analyzed[index_order]++;
    pthread_mutex_unlock(&test->mutex);

    if (best_logodds) {
        *best_logodds = match->logodds;
    }
    if (best_fieldnum) {
        *best_fieldnum = match->fieldnum;
    }
    return TRUE;
}

static int canonical_index_order_merge(
    onefield_t* master_bp,
    bl* solutions,
    anbool* solved_out) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;
    MatchObj* match;

    if (solved_out) {
        *solved_out = FALSE;
    }
    if (!master_bp ||
        !solutions ||
        bl_size(solutions) != 1U) {
        return -1;
    }
    match = bl_access(solutions, 0);

    pthread_mutex_lock(&test->mutex);
    test->merge_count++;
    test->merged_order = match->indexid;
    pthread_mutex_unlock(&test->mutex);

    master_bp->single_field_solved = TRUE;
    bl_remove_all(solutions);
    if (solved_out) {
        *solved_out = TRUE;
    }
    return 0;
}

static void canonical_index_order_free_solutions(
    bl* solutions) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;

    if (solutions && bl_size(solutions) == 1U) {
        MatchObj* match = bl_access(solutions, 0);

        pthread_mutex_lock(&test->mutex);
        if (match->indexid == 1) {
            test->loser_freed++;
        } else {
            test->failed = TRUE;
        }
        pthread_mutex_unlock(&test->mutex);
    } else {
        pthread_mutex_lock(&test->mutex);
        test->failed = TRUE;
        pthread_mutex_unlock(&test->mutex);
    }
    bl_free(solutions);
}

static const index_shard_hooks_t canonical_index_order_hooks = {
    canonical_index_order_get_index,
    canonical_index_order_done_with_index,
    canonical_index_order_report_solution,
    canonical_index_order_prepare_local,
    canonical_index_order_reset_local,
    canonical_index_order_cleanup_local,
    canonical_index_order_solve_one,
    canonical_index_order_analyze,
    canonical_index_order_merge,
    canonical_index_order_free_solutions};

static int run_canonical_index_order_test(void) {
    canonical_index_order_test_t* test =
        &canonical_index_order_test;
    index_shard_solve_status_t status =
        INDEX_SHARD_SOLVE_PRECOMMIT_FAILURE;
    onefield_t bp;
    anbool pool_started = FALSE;
    int result = 1;

    memset(test, 0, sizeof(*test));
    test->merged_order = -1;
    test->reported_order = SIZE_MAX;
    if (pthread_mutex_init(&test->mutex, NULL)) {
        fprintf(stderr, "failed to initialize canonical-order mutex\n");
        return 1;
    }
    if (pthread_cond_init(&test->condition, NULL)) {
        fprintf(
            stderr,
            "failed to initialize canonical-order condition\n");
        pthread_mutex_destroy(&test->mutex);
        return 1;
    }

    test->indexes[0].indexname = "canonical-index-order-0";
    test->indexes[0].indexid = 0;
    test->indexes[1].indexname = "canonical-index-order-1";
    test->indexes[1].indexid = 1;

    log_init(LOG_ALL);
    onefield_init(&bp);
    solver_set_default_values(&bp.solver);
    /* Two producer lanes plus one index-free helper lane. */
    bp.index_shard_workers = 3;

    if (index_shard_pool_start(&bp, &bp.solver)) {
        fprintf(stderr, "failed to start canonical-order pool\n");
        goto cleanup;
    }
    pool_started = TRUE;

    status = index_shard_solve(
        &bp,
        &bp.solver,
        2U,
        &canonical_index_order_hooks);

    index_shard_pool_stop(&bp);
    pool_started = FALSE;

    if (status != INDEX_SHARD_SOLVE_HANDLED ||
        test->failed ||
        test->completion_count != 2 ||
        test->completion_order[0] != 1 ||
        test->completion_order[1] != 0 ||
        test->analyzed[0] != 1 ||
        test->analyzed[1] != 1 ||
        test->merge_count != 1 ||
        test->merged_order != 0 ||
        test->report_count != 1 ||
        test->reported_order != 0U ||
        test->loser_freed != 1 ||
        !bp.single_field_solved) {
        fprintf(
            stderr,
            "canonical index-order test failed: "
            "status=%i failed=%i completion=%i:%i,%i "
            "analyzed=%i,%i merged=%i:%i reported=%i:%zu "
            "loser_freed=%i solved=%i\n",
            (int)status,
            test->failed,
            test->completion_count,
            test->completion_order[0],
            test->completion_order[1],
            test->analyzed[0],
            test->analyzed[1],
            test->merge_count,
            test->merged_order,
            test->report_count,
            test->reported_order,
            test->loser_freed,
            bp.single_field_solved);
        goto cleanup;
    }

    printf(
        "INDEX_SHARD_CANONICAL_ORDER_OK "
        "completion=1,0 analyzed=1,1 "
        "merged=0 reported=0 loser_freed=1\n");
    result = 0;

cleanup:
    if (pool_started) {
        index_shard_pool_stop(&bp);
    }
    solver_cleanup(&bp.solver);
    onefield_cleanup(&bp);
    pthread_cond_destroy(&test->condition);
    pthread_mutex_destroy(&test->mutex);
    return result;
}

int main(int argc, char** argv) {
    onefield_t bp;
    const char* cancel_file;
    const char* mode;
    char* match_path = NULL;
    result_signature_t signature;
    double total_cpu_limit = 0.0;
    double total_wall_limit = 0.0;
    int cpu_limit_status;
    int wall_limit_status;
    int first_object;
    int last_object;
    int second_first_object = 0;
    int second_last_object = 0;
    int pass_count = 1;
    int pass_number;
    int workers;
    int index_arg;
    int result = 0;
    anbool pool_started = FALSE;

    if (argc == 2 &&
        !strcmp(argv[1], "--canonical-index-order")) {
        return run_canonical_index_order_test();
    }

    if (argc < 8) {
        fprintf(
            stderr,
            "usage: %s FIELD_XYLS WORKERS FIRST_OBJECT LAST_OBJECT "
            "WCS_OUT exhaustive|probe|winner|cancel|limit|multipass "
            "INDEX_FITS...\n",
            argv[0]);
        return 2;
    }
    if (parse_positive_int(argv[2], "WORKERS", &workers) ||
        parse_positive_int(
            argv[3],
            "FIRST_OBJECT",
            &first_object) ||
        parse_positive_int(
            argv[4],
            "LAST_OBJECT",
            &last_object)) {
        return 2;
    }
    if (last_object < first_object) {
        fprintf(stderr, "LAST_OBJECT must not precede FIRST_OBJECT\n");
        return 2;
    }
    mode = argv[6];
    if (strcmp(mode, "exhaustive") &&
        strcmp(mode, "probe") &&
        strcmp(mode, "winner") &&
        strcmp(mode, "cancel") &&
        strcmp(mode, "limit") &&
        strcmp(mode, "multipass")) {
        fprintf(stderr, "unknown mode: %s\n", mode);
        return 2;
    }
    cancel_file = getenv("AB_CANCEL_FILE");
    if (!strcmp(mode, "cancel") &&
        (!cancel_file || !cancel_file[0])) {
        fprintf(
            stderr,
            "cancel mode requires AB_CANCEL_FILE\n");
        return 2;
    }
    wall_limit_status = parse_positive_double_env(
        "AB_TOTAL_WALL_LIMIT",
        &total_wall_limit);
    cpu_limit_status = parse_positive_double_env(
        "AB_TOTAL_CPU_LIMIT",
        &total_cpu_limit);
    if (wall_limit_status < 0 || cpu_limit_status < 0) {
        return 2;
    }
    if (!strcmp(mode, "limit") &&
        !wall_limit_status &&
        !cpu_limit_status) {
        fprintf(
            stderr,
            "limit mode requires AB_TOTAL_WALL_LIMIT or "
            "AB_TOTAL_CPU_LIMIT\n");
        return 2;
    }
    if (!strcmp(mode, "multipass")) {
        if (parse_positive_int(
                getenv("AB_SECOND_FIRST_OBJECT"),
                "AB_SECOND_FIRST_OBJECT",
                &second_first_object) ||
            parse_positive_int(
                getenv("AB_SECOND_LAST_OBJECT"),
                "AB_SECOND_LAST_OBJECT",
                &second_last_object) ||
            second_last_object < second_first_object) {
            fprintf(
                stderr,
                "multipass mode requires a valid second object range\n");
            return 2;
        }
        pass_count = 2;
    }

    log_init(LOG_ALL);
    onefield_init(&bp);
    solver_set_default_values(&bp.solver);
    if (onefield_job_field_cache_begin(&bp)) {
        fprintf(stderr, "failed to initialize integration field cache\n");
        onefield_cleanup(&bp);
        return 1;
    }

    onefield_set_field_file(&bp, argv[1]);
    onefield_add_field(&bp, 1);
    for (index_arg = 7; index_arg < argc; index_arg++) {
        onefield_add_index(&bp, argv[index_arg]);
    }
    onefield_set_wcs_file(&bp, argv[5]);
    match_path = malloc(strlen(argv[5]) + sizeof(".match"));
    if (!match_path) {
        result = 1;
        goto cleanup;
    }
    snprintf(
        match_path,
        strlen(argv[5]) + sizeof(".match"),
        "%s.match",
        argv[5]);
    if (unlink(match_path) && errno != ENOENT) {
        fprintf(stderr, "failed to remove stale match artifact\n");
        result = 1;
        goto cleanup;
    }
    onefield_set_match_file(&bp, match_path);
    if (cancel_file) {
        onefield_set_cancel_file(&bp, cancel_file);
    }

    bp.index_shard_workers = workers;
    if (!strcmp(mode, "winner")) {
        /*
         * Accept a verified astronomical match, not the first arbitrary
         * negative-log-odds CodeKD candidate.  This makes a valid nonmatching
         * index usable ahead of the known APOD4 winner.
         */
        bp.logratio_tosolve = log(1.0e6);
        bp.solver.logratio_toprint = log(1.0e6);
        bp.solver.logratio_tokeep = log(1.0e6);
    } else {
        bp.logratio_tosolve = HUGE_VAL;
        bp.solver.logratio_toprint = log(1.0e6);
        bp.solver.logratio_tokeep = log(1.0e9);
    }
    bp.solver.startobj = first_object - 1;
    bp.solver.endobj = last_object;
    bp.solver.logratio_totune = HUGE_VAL;
    bp.solver.quadsize_min = 0.1 * 507.0;
    bp.solver.funits_lower = 3600.0 * 30.0 / 719.0;
    bp.solver.funits_upper = 3600.0 * 40.0 / 719.0;
    bp.solver.do_tweak = FALSE;
    bp.total_timelimit = total_wall_limit;
    bp.total_cpulimit = (float)total_cpu_limit;
    solver_set_field_bounds(
        &bp.solver,
        0.0,
        719.0,
        0.0,
        507.0);

    if (workers > 1) {
        if (index_shard_pool_start(&bp, &bp.solver)) {
            result = 1;
            goto cleanup;
        }
        pool_started = TRUE;
    }

    for (pass_number = 1;
         pass_number <= pass_count;
         pass_number++) {
        if (pass_number == 2) {
            if (touch_before_second_pass(
                    "AB_TOUCH_FIELD_BEFORE_SECOND") ||
                touch_before_second_pass(
                    "AB_TOUCH_INDEX_BEFORE_SECOND")) {
                result = 1;
                goto cleanup;
            }
            first_object = second_first_object;
            last_object = second_last_object;
            bp.solver.startobj = first_object - 1;
            bp.solver.endobj = last_object;
        }
        if ((unlink(match_path) && errno != ENOENT) ||
            (unlink(argv[5]) && errno != ENOENT)) {
            fprintf(stderr, "failed to remove stale pass artifact\n");
            result = 1;
            goto cleanup;
        }

        onefield_run(&bp);
        if (read_result_signature(match_path, &signature)) {
            fprintf(stderr, "failed to read integration match artifact\n");
            result = 1;
            goto cleanup;
        }
        print_result_record(
            &bp,
            &signature,
            mode,
            pass_number,
            workers,
            first_object,
            last_object,
            argc - 7);
        if (print_match_records(
                match_path,
                mode,
                pass_number,
                workers,
                first_object,
                last_object)) {
            fprintf(stderr, "failed to read integration match records\n");
            result = 1;
            goto cleanup;
        }
        if (bp.solver_failed) {
            result = 1;
            goto cleanup;
        }
        if (!strcmp(mode, "cancel")) {
            if (!bp.cancelled) {
                fprintf(
                    stderr,
                    "integration cancellation was not observed\n");
                result = 1;
            }
            goto cleanup;
        }
        if (!strcmp(mode, "limit")) {
            if (!bp.hit_total_timelimit &&
                !bp.hit_total_cpulimit) {
                fprintf(
                    stderr,
                    "integration limit was not observed\n");
                result = 1;
            }
            goto cleanup;
        }
        if (print_wcs_record(
                argv[5],
                mode,
                pass_number,
                workers,
                first_object,
                last_object,
                !strcmp(mode, "probe") ||
                    !strcmp(mode, "multipass"))) {
            result = 1;
            goto cleanup;
        }
    }

cleanup:
    if (pool_started) {
        index_shard_pool_stop(&bp);
    }
    solver_cleanup(&bp.solver);
    onefield_cleanup(&bp);
    free(match_path);
    return result;
}
