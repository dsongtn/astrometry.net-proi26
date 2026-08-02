/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef SOLVER_HYPOTHESIS_INTERNAL_H
#define SOLVER_HYPOTHESIS_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include <sys/resource.h>
#include <time.h>

#include "solver.h"
#include "verify.h"
#include "pquad.h"
#include "index_shard_internal.h"
#include "solver_field_geometry_internal.h"

#ifndef SOLVER_AB_CANDIDATE_LIMIT_BYTES
#define SOLVER_AB_CANDIDATE_LIMIT_BYTES (768U * 1024U)
#endif
#define SOLVER_AB_DESCRIPTOR_CAPACITY 4096U
#define SOLVER_AB_DESCRIPTOR_MIN_COMBINATIONS 64U
#define SOLVER_AB_DESCRIPTOR_DELIVERY_MIN_HYPOTHESES 64U
#define SOLVER_AB_DESCRIPTOR_TASK_TARGET_HYPOTHESES \
    (2U * SOLVER_AB_DESCRIPTOR_DELIVERY_MIN_HYPOTHESES)
#define SOLVER_AB_DESCRIPTOR_MAX_FIELD_OBJECTS 1000
#define SOLVER_AB_DESCRIPTOR_PAIR_CACHE_BYTES \
    (8U * 1024U * 1024U)
#define SOLVER_PAYLOAD_CANDIDATE_BATCH 32U
#define SOLVER_VERIFICATION_WINDOW_CANDIDATES \
    (2U * SOLVER_PAYLOAD_CANDIDATE_BATCH)

typedef enum solver_ab_phase_kind {
    SOLVER_AB_PHASE_DIAGONAL = 0,
    SOLVER_AB_PHASE_OFF_DIAGONAL = 1
} solver_ab_phase_kind_t;

typedef enum solver_ab_phase_mode {
    SOLVER_AB_MODE_NATIVE = 0,
    SOLVER_AB_MODE_FLATTENED_OWNER = 1,
    SOLVER_AB_MODE_ASSISTED = 2
} solver_ab_phase_mode_t;

typedef enum solver_ab_candidate_action {
    SOLVER_AB_CANDIDATE_RADEC_SKIP = 0,
    SOLVER_AB_CANDIDATE_ABSCALE_SKIP = 1,
    SOLVER_AB_CANDIDATE_BAD_QUAD = 2,
    SOLVER_AB_CANDIDATE_SCALE_SKIP = 3,
    SOLVER_AB_CANDIDATE_VERIFY = 4
} solver_ab_candidate_action_t;

typedef enum solver_packet_reserve_result {
    SOLVER_PACKET_RESERVE_ERROR = -1,
    SOLVER_PACKET_RESERVE_OK = 0,
    SOLVER_PACKET_RESERVE_FULL = 1
} solver_packet_reserve_result_t;

typedef struct solver_ab_pair {
    int field_a;
    int field_b;
    unsigned long long combination_first;
    unsigned long long combination_count;
    double tol2;
} solver_ab_pair_t;

typedef struct solver_ab_descriptor {
    int stars[DQMAX];
    double code[DCMAX];
    double tol2;
    double rel_field_noise2;
    unsigned long long numtries_delta;
    unsigned long long cxdx_delta;
    unsigned long long meanx_delta;
    anbool current_parity;
} solver_ab_descriptor_t;

typedef struct solver_ab_descriptor_output {
    size_t descriptor_count;
    unsigned long long trailing_numtries;
    unsigned long long trailing_cxdx;
    unsigned long long trailing_meanx;
    double final_rel_field_noise2;
    anbool has_final_rel_field_noise2;
    solver_ab_descriptor_t descriptors[SOLVER_AB_DESCRIPTOR_CAPACITY];
} solver_ab_descriptor_output_t;

typedef struct solver_ab_candidate {
    solver_ab_candidate_action_t action;
    int quadno;
    double code_err;
    tan_t wcs;
    double scale;
    anbool parity;
    int quad_npeers;
    unsigned int star[DQMAX];
    int field[DQMAX];
    double quadpix[2 * DQMAX];
    double quadxyz[3 * DQMAX];
} solver_ab_candidate_t;

typedef struct solver_verification_packet {
    solver_ab_candidate_t* candidates;
    size_t candidate_capacity;
    anbool allocation_failed;
} solver_verification_packet_t;

typedef struct solver_descriptor_status {
    anbool evaluation_failed;
    anbool cancelled;
} solver_descriptor_status_t;

typedef struct solver_ab_task {
    unsigned long long combination_first;
    unsigned long long combination_end;
} solver_ab_task_t;

typedef struct solver_ab_descriptor_planner
    solver_ab_descriptor_planner_t;

typedef struct solver_ab_builder {
    solver_ab_descriptor_planner_t* planner;
    solver_descriptor_status_t* status;
    double tol2;
    double rel_field_noise2;
    unsigned long long pending_numtries;
    unsigned long long pending_cxdx;
    unsigned long long pending_meanx;
    solver_ab_descriptor_output_t* descriptor_output;
    anbool rel_field_noise_valid;
    anbool fatal_error;
} solver_ab_builder_t;

typedef struct solver_ab_snapshot {
    index_t* index;
    const starxy_t* fieldxy;
    anbool use_radec;
    anbool cx_less_than_dx;
    anbool meanx_less_than_half;
    int parity;
    double centerxyz[3];
    double r2;
    double abscale_low;
    double abscale_high;
    double funits_lower;
    double funits_upper;
    double cxdx_margin;
    double codetol;
    double rel_index_noise2;
} solver_ab_snapshot_t;

typedef struct solver_ab_phase_telemetry {
    anbool enabled;
    anbool resource_valid;
    double wall_start;
    struct rusage resource_start;
    int combinations_start;
    int candidates_start;
    unsigned long long codekd_calls_start;
    unsigned long long codekd_hits_start;
    unsigned long long verify_calls_start;
} solver_ab_phase_telemetry_t;

struct solver_ab_descriptor_planner {
    solver_ab_snapshot_t snapshot;
    const solver_field_geometry_t* field_geometry;
    int newpoint;
    int dimquads;
    solver_ab_phase_kind_t phase;
    solver_ab_pair_t* pairs;
    size_t pair_count;
    size_t pair_capacity;
    solver_ab_pair_t* pair_cache;
    size_t pair_cache_capacity;
    size_t pair_cache_limit_bytes;
    anbool pairs_transient;
    int* combination_eligible;
    size_t combination_eligible_capacity;
};

typedef struct solver_verification_score_slot {
    verify_prepared_hit_t* prepared;
} solver_verification_score_slot_t;

typedef struct solver_verification_candidate_runtime {
    MatchObj match;
    double match_distance_in_pixels2;
    double logaccept;
} solver_verification_candidate_runtime_t;

typedef struct solver_verification_task_input {
    const solver_verification_score_slot_t* slots;
    size_t slot_first;
    size_t slot_count;
} solver_verification_task_input_t;

typedef struct solver_verification_retire_context {
    size_t next_slot;
} solver_verification_retire_context_t;

typedef struct solver_ab_descriptor_task_input {
    const solver_field_geometry_t* field_geometry;
    const solver_ab_pair_t* pairs;
    size_t pair_count;
    unsigned long long combination_first;
    unsigned long long combination_end;
    solver_ab_phase_kind_t phase;
    int newpoint;
    int dimquads;
    int parity;
    anbool cx_less_than_dx;
    anbool meanx_less_than_half;
    double cxdx_margin;
} solver_ab_descriptor_task_input_t;

typedef struct solver_ab_descriptor_workspace {
    solver_ab_descriptor_planner_t planner;
    solver_ab_descriptor_output_t* outputs;
    size_t output_capacity;
} solver_ab_descriptor_workspace_t;

extern const index_shard_helper_ops_t solver_verification_helper_ops;

uint64_t solver_hypothesis_order_digest(
    const int* stars,
    const double* code,
    int dimquad,
    anbool current_parity);
uint64_t solver_kd_result_order_digest(const kdtree_qres_t* result);
void solver_record_candidate_order(
    solver_t* solver,
    solver_ab_candidate_action_t action,
    int quadno,
    double code_err);
anbool solver_ab_poll_phase_stop(
    solver_t* solver,
    time_t* next_timer_callback_time);
void solver_ab_phase_telemetry_begin(
    const solver_t* solver,
    solver_ab_phase_telemetry_t* telemetry);
void solver_ab_phase_telemetry_report(
    const solver_t* solver,
    const solver_ab_phase_telemetry_t* telemetry,
    int newpoint,
    solver_ab_phase_kind_t phase,
    solver_ab_phase_mode_t mode);
unsigned long long solver_ab_saturating_add(
    unsigned long long a,
    unsigned long long b);
solver_packet_reserve_result_t solver_verification_packet_reserve(
    solver_verification_packet_t* packet,
    size_t required);
void solver_verification_packet_free(
    solver_verification_packet_t* packet);
int solver_ab_try_verification_wave(
    kdtree_qres_t* result,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    int candidate_first,
    int candidate_end);
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
    unsigned long long* total_combinations_out);
void solver_ab_descriptor_release_pairs(
    solver_ab_descriptor_workspace_t* workspace);
int solver_ab_descriptor_workspace_reserve_outputs(
    solver_ab_descriptor_workspace_t* workspace,
    size_t output_count);
solver_ab_descriptor_workspace_t* solver_ab_descriptor_workspace_get(void);
index_shard_helper_task_status_t solver_ab_descriptor_helper_execute(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size);
size_t solver_ab_descriptor_expansion(int dimquads, int parity);
int solver_ab_checked_counter_delta(
    solver_t* solver,
    unsigned long long numtries,
    unsigned long long cxdx,
    unsigned long long meanx);
size_t solver_ab_descriptor_partition_count(
    unsigned long long wave_combinations,
    size_t expansion,
    size_t max_task_combinations,
    size_t participants);
anbool solver_payload_candidate_data_fully_resident(
    const solver_t* solver);
size_t solver_verification_wave_memory_budget(void);

void set_matchobj_template(solver_t* solver, MatchObj* mo);
void set_center_and_radius(
    solver_t* solver,
    MatchObj* mo,
    tan_t* tan,
    sip_t* sip);
double get_tolerance_for_noise(
    double codetol,
    double rel_field_noise2,
    double rel_index_noise2);
int solver_handle_hit(
    solver_t* solver,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match);
double solver_prepare_hit_for_verify(
    solver_t* solver,
    MatchObj* mo,
    double* logaccept);
int solver_handle_hit_after_verify(
    solver_t* solver,
    MatchObj* mo,
    sip_t* verifysip,
    anbool fake_match,
    double match_distance_in_pixels2);

#endif
