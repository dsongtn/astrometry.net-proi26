#include <math.h>
#include <stdio.h>
#include <string.h>

#include "astrometry/engine.h"
#include "astrometry/solverutils.h"
#include "engine_internal.h"

static int failures = 0;

#define CHECK(expression)                                                \
    do {                                                                 \
        if (!(expression)) {                                             \
            fprintf(stderr, "CHECK failed at %s:%i: %s\n",             \
                    __FILE__, __LINE__, #expression);                    \
            failures++;                                                  \
        }                                                                \
    } while (0)

static anbool close_enough(double left, double right) {
    return fabs(left - right) < 1.0e-12;
}

static void bounded_then_open_resets_endobj(void) {
    job_t job;
    solver_t solver;
    engine_pass_cursor_t cursor;
    engine_pass_t pass;

    memset(&job, 0, sizeof(job));
    memset(&solver, 0, sizeof(solver));
    job.depths = il_new(4);
    job.scales = dl_new(2);
    CHECK(job.depths != NULL);
    CHECK(job.scales != NULL);
    CHECK(parse_depth_string(job.depths, "5-10 30-") == 0);
    dl_append(job.scales, 0.5);
    dl_append(job.scales, 1.0);

    solver.endobj = 777;
    engine_pass_cursor_init(&cursor);
    CHECK(engine_pass_cursor_next(
        &job, 7.5, 9.5, &cursor, &pass));
    CHECK(pass.ordinal == 0U);
    engine_pass_apply(&solver, &pass);
    CHECK(solver.startobj == 4);
    CHECK(solver.endobj == 10);

    CHECK(engine_pass_cursor_next(
        &job, 7.5, 9.5, &cursor, &pass));
    CHECK(pass.ordinal == 1U);
    engine_pass_apply(&solver, &pass);
    CHECK(solver.startobj == 29);
    CHECK(solver.endobj == 0);
    CHECK(!engine_pass_cursor_next(
        &job, 7.5, 9.5, &cursor, &pass));

    il_free(job.depths);
    dl_free(job.scales);
}

static void range_major_scale_minor_order(void) {
    static const int expected_start[] = {
        4, 4, 4, 1, 1, 1, 29, 29, 29
    };
    static const int expected_end[] = {
        10, 10, 10, 3, 3, 3, 0, 0, 0
    };
    static const double expected_lower[] = {
        0.5, 2.0, 7.5, 0.5, 2.0, 7.5, 0.5, 2.0, 7.5
    };
    static const double expected_upper[] = {
        1.0, 4.0, 9.5, 1.0, 4.0, 9.5, 1.0, 4.0, 9.5
    };
    job_t job;
    solver_t solver;
    engine_pass_cursor_t cursor;
    engine_pass_t pass;
    size_t ordinal;

    memset(&job, 0, sizeof(job));
    memset(&solver, 0, sizeof(solver));
    job.depths = il_new(6);
    job.scales = dl_new(6);
    CHECK(job.depths != NULL);
    CHECK(job.scales != NULL);
    CHECK(parse_depth_string(job.depths, "5-10 2-3 30-") == 0);
    dl_append(job.scales, 0.5);
    dl_append(job.scales, 1.0);
    dl_append(job.scales, 2.0);
    dl_append(job.scales, 4.0);
    dl_append(job.scales, 0.0);
    dl_append(job.scales, 0.0);

    engine_pass_cursor_init(&cursor);
    for (ordinal = 0U; ordinal < 9U; ordinal++) {
        CHECK(engine_pass_cursor_next(
            &job, 7.5, 9.5, &cursor, &pass));
        CHECK(pass.ordinal == ordinal);
        CHECK(pass.depth_index == ordinal / 3U);
        CHECK(pass.scale_index == ordinal % 3U);
        CHECK(pass.startobj == expected_start[ordinal]);
        CHECK(pass.endobj == expected_end[ordinal]);
        CHECK(close_enough(
            pass.funits_lower, expected_lower[ordinal]));
        CHECK(close_enough(
            pass.funits_upper, expected_upper[ordinal]));

        solver.startobj = -1;
        solver.endobj = -1;
        engine_pass_apply(&solver, &pass);
        CHECK(solver.startobj == expected_start[ordinal]);
        CHECK(solver.endobj == expected_end[ordinal]);
        CHECK(close_enough(
            solver.funits_lower, expected_lower[ordinal]));
        CHECK(close_enough(
            solver.funits_upper, expected_upper[ordinal]));
    }
    CHECK(!engine_pass_cursor_next(
        &job, 7.5, 9.5, &cursor, &pass));

    il_free(job.depths);
    dl_free(job.scales);
}

static void primary_wall_and_optional_cpu_policy(void) {
    engine_limit_policy_t limits;
    engine_t* engine;

    engine = engine_new();
    CHECK(engine != NULL);
    if (engine) {
        CHECK(close_enough(engine->walllimit, 300.0));
        CHECK(close_enough(engine->cpulimit, 0.0));
        engine_free(engine);
    }

    {
        FILE* config = tmpfile();

        CHECK(config != NULL);
        if (config) {
            engine = engine_new();
            CHECK(engine != NULL);
            CHECK(fputs("walllimit 240\ncpulimit 0\n", config) >= 0);
            rewind(config);
            if (engine) {
                CHECK(engine_parse_config_file_stream(engine, config) == 0);
                CHECK(close_enough(engine->walllimit, 240.0));
                CHECK(close_enough(engine->cpulimit, 0.0));
                engine_free(engine);
            }
            fclose(config);
        }
    }

    engine_limit_policy_resolve(0.0, 300.0, 0.0, 0.0, &limits);
    CHECK(close_enough(limits.wall_seconds, 300.0));
    CHECK(close_enough(limits.cpu_seconds, 0.0));
    CHECK(!limits.wall_from_job);
    CHECK(!limits.wall_job_clamped);
    CHECK(!limits.cpu_from_job);

    engine_limit_policy_resolve(60.0, 300.0, 720.0, 0.0, &limits);
    CHECK(close_enough(limits.wall_seconds, 60.0));
    CHECK(close_enough(limits.cpu_seconds, 720.0));
    CHECK(limits.wall_from_job);
    CHECK(!limits.wall_job_clamped);
    CHECK(limits.cpu_from_job);

    engine_limit_policy_resolve(600.0, 300.0, 720.0, 120.0, &limits);
    CHECK(close_enough(limits.wall_seconds, 300.0));
    CHECK(close_enough(limits.cpu_seconds, 720.0));
    CHECK(!limits.wall_from_job);
    CHECK(limits.wall_job_clamped);
    CHECK(limits.cpu_from_job);

    engine_limit_policy_resolve(600.0, 0.0, 0.0, 120.0, &limits);
    CHECK(close_enough(limits.wall_seconds, 600.0));
    CHECK(close_enough(limits.cpu_seconds, 120.0));
    CHECK(limits.wall_from_job);
    CHECK(!limits.wall_job_clamped);
    CHECK(!limits.cpu_from_job);
}

int main(void) {
    bounded_then_open_resets_endobj();
    range_major_scale_minor_order();
    primary_wall_and_optional_cpu_policy();
    if (failures) {
        fprintf(stderr, "ENGINE_PASS_TEST_FAILED failures=%i\n", failures);
        return 1;
    }
    printf("ENGINE_PASS_TEST_OK cases=3 ordered_passes=9\n");
    return 0;
}
