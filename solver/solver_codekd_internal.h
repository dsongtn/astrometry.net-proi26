/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef SOLVER_CODEKD_INTERNAL_H
#define SOLVER_CODEKD_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "solver.h"
#include "verify.h"
#include "fitsbin.h"
#include "index_shard_internal.h"
#include "solver_hypothesis_internal.h"
#include "../libkd/kdtree_prefetch_internal.h"

#define SOLVER_CODEKD_SEARCH_OPTIONS \
    (KD_OPTIONS_SMALL_RADIUS | \
     KD_OPTIONS_COMPUTE_DISTS | \
     KD_OPTIONS_NO_RESIZE_RESULTS | \
     KD_OPTIONS_USE_SPLIT)
#define SOLVER_STARKD_VERIFY_SEARCH_OPTIONS \
    (KD_OPTIONS_SMALL_RADIUS | KD_OPTIONS_RETURN_POINTS)

#define SOLVER_CODEKD_PACKET_RESULT_LIMIT_BYTES \
    (2U * 1024U * 1024U)
#define SOLVER_CODEKD_DELIVERY_RANGE_CAPACITY \
    FITSBIN_MMAP_PREFETCH_RANGE_LIMIT
#define SOLVER_CODEKD_DELIVERY_BUDGET_BYTES \
    (2U * 1024U * 1024U)
#define SOLVER_CANDIDATE_DELIVERY_LIMIT \
    (FITSBIN_PREAD_ASYNC_RANGE_LIMIT / 2U)
#define SOLVER_VERIFY_QUERY_LOOKAHEAD 16U
#ifndef SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED
#define SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED 0
#endif
#if SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED != 0 && \
    SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED != 1
#error "SOLVER_VERIFY_SWEEP_ASYNC_DIRECT_ENABLED must be 0 or 1"
#endif
#define SOLVER_CANDIDATE_STAR_LIMIT \
    (SOLVER_CANDIDATE_DELIVERY_LIMIT * DQMAX)
#if SOLVER_CANDIDATE_STAR_LIMIT > FITSBIN_MMAP_PREFETCH_RANGE_LIMIT
#error "candidate Star delivery exceeds mapped-prefetch range capacity"
#endif

typedef enum solver_codekd_result_state {
    SOLVER_CODEKD_RESULT_UNUSED = 0,
    SOLVER_CODEKD_RESULT_READY,
    SOLVER_CODEKD_RESULT_OWNER_REPLAY,
    SOLVER_CODEKD_RESULT_QUERY_FAILED
} solver_codekd_result_state_t;

typedef struct solver_codekd_result_slot {
    size_t hit_first;
    unsigned int hit_count;
    double search_wall_seconds;
    int search_errno;
    solver_codekd_result_state_t state;
} solver_codekd_result_slot_t;

typedef enum solver_codekd_search_packet_state {
    SOLVER_CODEKD_PACKET_ALLOCATED = 0,
    SOLVER_CODEKD_PACKET_DESCRIPTORS_READY,
    SOLVER_CODEKD_PACKET_PAGE_PLAN_COMPLETE,
    SOLVER_CODEKD_PACKET_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_QUAD_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_QUAD_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_QUAD_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_STAR_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_STAR_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_STAR_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_VERIFY_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_VERIFY_PAGE_PLAN_COMPLETE,
    SOLVER_CODEKD_PACKET_VERIFY_QUERY_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_SUBMIT_READY,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_IO_SUBMITTED,
    SOLVER_CODEKD_PACKET_VERIFY_SWEEP_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_PREPARE_OWNER,
    SOLVER_CODEKD_PACKET_VERIFY_SCORE_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_VERIFY_COMPUTE_READY,
    SOLVER_CODEKD_PACKET_EXECUTING,
    SOLVER_CODEKD_PACKET_RESULTS_READY,
    SOLVER_CODEKD_PACKET_RETIRING,
    SOLVER_CODEKD_PACKET_RETIRED,
    SOLVER_CODEKD_PACKET_STOPPED,
    SOLVER_CODEKD_PACKET_FAILED
} solver_codekd_search_packet_state_t;

typedef enum solver_codekd_page_plan_reason {
    SOLVER_CODEKD_PAGE_PLAN_NONE = 0,
    SOLVER_CODEKD_PAGE_PLAN_NOT_APPLICABLE,
    SOLVER_CODEKD_PAGE_PLAN_ALLOCATION,
    SOLVER_CODEKD_PAGE_PLAN_SOURCE_MISMATCH,
    SOLVER_CODEKD_PAGE_PLAN_INVALID_RANGE,
    SOLVER_CODEKD_PAGE_PLAN_BYTE_BUDGET,
    SOLVER_CODEKD_PAGE_PLAN_RANGE_CAPACITY,
    SOLVER_CODEKD_PAGE_PLAN_SERVICE_REFUSED,
    SOLVER_CODEKD_PAGE_PLAN_SERVICE_ERROR,
    SOLVER_CODEKD_PAGE_PLAN_CANCELLED
} solver_codekd_page_plan_reason_t;

typedef struct solver_codekd_page_entry {
    uintptr_t mapping_begin;
    uintptr_t mapping_end;
    uintptr_t page_key;
    uintptr_t populate_begin;
    uintptr_t populate_end;
} solver_codekd_page_entry_t;

typedef struct solver_codekd_page_set {
    solver_codekd_page_entry_t* entries;
    size_t* hash_indices;
    unsigned int* hash_generations;
    size_t count;
    size_t capacity;
    size_t hash_capacity;
    unsigned int generation;
} solver_codekd_page_set_t;

typedef struct solver_codekd_page_workspace {
    solver_codekd_page_set_t descriptor;
    solver_codekd_page_set_t group;
    solver_codekd_page_entry_t* sort_entries;
    solver_codekd_page_entry_t* sort_scratch;
    fitsbin_prefetch_range_t* sealed_ranges;
    size_t page_size;
    size_t page_limit;
    size_t sealed_range_capacity;
} solver_codekd_page_workspace_t;

typedef struct solver_codekd_page_plan_stats {
    size_t descriptors_total;
    size_t descriptors_planned;
    size_t descriptor_splits;
    size_t raw_ranges;
    size_t unique_pages;
    size_t ranges_before_dedup;
    size_t ranges_after_dedup;
    size_t logical_bytes;
    size_t aligned_bytes;
    size_t overread_bytes;
    size_t refusal_counts[SOLVER_CODEKD_PAGE_PLAN_CANCELLED + 1U];
} solver_codekd_page_plan_stats_t;

typedef struct solver_codekd_verification_snapshot {
    MatchObj match_template;
    const verify_field_t* field;
    int index_cutnside;
    int indexid;
    int healpix;
    int hpnside;
    double index_jitter;
    double verify_pix;
    double distractor_ratio;
    double logratio_bail_threshold;
    double logaccept;
    double logratio_stoplooking;
    anbool distance_from_quad_bonus;
    anbool enabled;
} solver_codekd_verification_snapshot_t;

typedef struct solver_candidate_delivery_record
    solver_candidate_delivery_record_t;

struct solver_candidate_delivery_record {
    unsigned int quadid;
    unsigned int stars[DQMAX];
    double starxyz[DQMAX * 3];
    int prepared_fieldstars[DQMAX];
    double prepared_fieldxy[DQMAX * 2];
    tan_t prepared_wcs;
    double prepared_scale;
    verify_index_query_t* verify_query;
    size_t descriptor_index;
    solver_ab_candidate_action_t plan_action;
    double verify_center[3];
    double verify_radius;
    double verify_radius2;
    double prepared_verify_pix2;
    double prepared_logaccept;
    double prepared_distractor_ratio;
    double prepared_logratio_bail_threshold;
    double prepared_logratio_stoplooking;
    double prepared_field_maxx;
    double prepared_field_maxy;
    verify_prepared_hit_t* prepared_verification;
    verify_prepared_score_t prepared_score;
    anbool prepared_parity;
    anbool candidate_prepared;
    anbool verify_delivery_fallback;
    anbool verify_query_captured;
    anbool prepared_distance_from_quad_bonus;
    anbool verification_score_ready;
};

typedef struct solver_codekd_search_packet {
    solver_ab_descriptor_output_t* descriptors;
    const kdtree_t* tree;
    const quadfile_t* quads;
    startree_t* starkd;
    solver_codekd_result_slot_t* slots;
    u32* inds;
    double* sdists;
    size_t hit_capacity;
    size_t hit_count;
    size_t first;
    size_t count;
    size_t next_descriptor;
    size_t plan_first;
    size_t plan_end;
    size_t plan_range_count;
    size_t plan_logical_bytes;
    size_t pending_descriptor_raw_ranges;
    size_t pending_descriptor_logical_bytes;
    size_t candidate_count;
    size_t candidate_capacity;
    size_t candidate_cursor;
    size_t candidate_window_first;
    size_t candidate_window_count;
    size_t candidate_window_offset;
    size_t candidate_star_count;
    size_t candidate_quad_ready_count;
    size_t candidate_star_ready_count;
    size_t candidate_verify_query_count;
    size_t verify_plan_first;
    size_t verify_plan_end;
    size_t verify_topology_end;
    size_t verify_plan_range_count;
    size_t verify_plan_logical_bytes;
    size_t verify_query_budget;
    size_t verify_query_bytes;
    size_t verify_prepared_count;
    size_t verify_prepared_retained_bytes;
    size_t verify_prepared_transient_bytes;
    size_t verify_sweep_range_count;
    size_t verify_sweep_aligned_bytes;
    size_t verify_sweep_storage_bytes;
    size_t verify_pending_query_index;
    size_t verify_pending_query_raw_ranges;
    size_t verify_pending_query_logical_bytes;
    size_t verify_delivery_budget;
    size_t retire_descriptor;
    size_t retire_hit_offset;
    unsigned long long sequence;
    unsigned int candidate_starids[SOLVER_CANDIDATE_STAR_LIMIT];
    solver_candidate_delivery_record_t* candidate_records;
    fitsbin_pread_range_t* verify_sweep_reads;
    verify_mapped_page_buffer_t* verify_sweep_buffers;
    unsigned char* verify_sweep_storage;
    fitsbin_t* delivery_source;
    fitsbin_payload_io_ticket_t* delivery_ticket;
    solver_codekd_page_workspace_t* page_workspace;
    solver_codekd_page_plan_stats_t page_stats;
    solver_codekd_page_plan_reason_t page_plan_reason;
    unsigned int candidate_quad_submitted;
    unsigned int candidate_quad_ready;
    unsigned int candidate_quad_fallback;
    unsigned int candidate_star_submitted;
    unsigned int candidate_star_ready;
    unsigned int candidate_star_fallback;
    unsigned long long candidate_delivery_windows;
    unsigned long long candidate_quad_ready_rows;
    unsigned long long candidate_star_ready_rows;
    unsigned long long candidate_retired_rows;
    unsigned long long candidate_native_rows;
    unsigned long long verification_page_queries;
    unsigned long long verification_page_queries_planned;
    unsigned long long verification_page_prefixes;
    unsigned long long verification_page_submitted;
    unsigned long long verification_page_ready;
    unsigned long long verification_page_fallback;
    unsigned long long verification_page_ready_rows;
    unsigned long long verification_page_ranges;
    unsigned long long verification_page_logical_bytes;
    unsigned long long verification_page_aligned_bytes;
    unsigned long long candidate_math_prepared;
    unsigned long long verification_score_batches_prepared;
    unsigned long long verification_score_contexts_prepared;
    unsigned long long verification_score_batches_executed;
    unsigned long long verification_score_contexts_completed;
    unsigned long long verification_score_work_units_completed;
    unsigned long long verification_score_fallback_batches;
    unsigned long long verification_score_stopped_batches;
    double verification_score_wall_seconds;
    int dimquads;
    anbool plan_complete;
    anbool pending_descriptor_plan;
    anbool detailed;
    anbool use_radec;
    anbool star_delivery_eligible;
    anbool candidate_quad_delivery_disabled;
    anbool candidate_star_delivery_disabled;
    anbool verification_delivery_disabled;
    anbool verify_plan_complete;
    anbool verify_sweep_plan_complete;
    anbool verify_pending_query_plan;
    anbool retire_descriptor_started;
    solver_codekd_search_packet_state_t state;
} solver_codekd_search_packet_t;

typedef struct solver_codekd_packet_task_input {
    solver_ab_descriptor_task_input_t descriptor;
    const kdtree_t* tree;
    const quadfile_t* quads;
    startree_t* starkd;
    anbool detailed;
    anbool use_radec;
    double field_minx;
    double field_maxx;
    double field_miny;
    double field_maxy;
    double abscale_low;
    double abscale_high;
    double funits_lower;
    double funits_upper;
    solver_codekd_verification_snapshot_t verification;
} solver_codekd_packet_task_input_t;

typedef struct solver_codekd_packet_wave {
    solver_codekd_packet_task_input_t* inputs;
    solver_codekd_search_packet_t* packets;
    index_shard_staged_task_t* tasks;
    size_t capacity;
} solver_codekd_packet_wave_t;

typedef struct solver_codekd_page_plan {
    fitsbin_t* source;
    solver_codekd_page_workspace_t* workspace;
    fitsbin_payload_io_cancel_check_fn cancelled;
    void* cancel_opaque;
    size_t raw_ranges;
    size_t logical_bytes;
    solver_codekd_page_plan_reason_t reason;
    anbool cancellation_observed;
    anbool enabled;
} solver_codekd_page_plan_t;

typedef struct solver_codekd_packet_retire_context {
    solver_t* solver;
    int dimquads;
    kdtree_qres_t** query_result;
    unsigned long long* reduced;
    size_t next_task_index;
    unsigned long long next_sequence;
} solver_codekd_packet_retire_context_t;

kdtree_qres_t* solver_codekd_rangesearch(
    const kdtree_t* tree,
    kdtree_qres_t* result,
    const double* query,
    double maxd2,
    int options);
void solver_execute_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    double tol2,
    kdtree_qres_t** presult);
void solver_begin_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity);
void solver_execute_prepared_hypothesis_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    kdtree_qres_t* result,
    double search_wall_seconds,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_quad_count,
    size_t prepared_star_count);
void solver_retire_codekd_failure_owner(
    const int* stars,
    const double* code,
    int dimquad,
    solver_t* solver,
    anbool current_parity,
    int search_errno,
    double search_wall_seconds);
void resolve_matches_native_range(
    kdtree_qres_t* result,
    const double* field_xy,
    const int* fieldstars,
    int dimquads,
    int quads_tried,
    solver_t* solver,
    anbool current_parity,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count);

solver_codekd_packet_wave_t* solver_codekd_packet_wave_create(
    size_t capacity);
void solver_codekd_packet_wave_destroy(
    solver_codekd_packet_wave_t* wave);
int solver_codekd_packet_finish_codekd(
    solver_codekd_search_packet_t* packet);
int solver_codekd_page_workspace_create(
    solver_codekd_page_workspace_t** workspace_out);
void solver_codekd_page_set_reset(solver_codekd_page_set_t* set);
int solver_codekd_page_set_merge_descriptor(
    solver_codekd_page_workspace_t* workspace);
int solver_codekd_page_set_union_fits(
    const solver_codekd_page_workspace_t* workspace);
void solver_codekd_page_plan_add_size(
    size_t* destination,
    size_t value);
int solver_codekd_page_plan_emit(
    void* opaque,
    const kdtree_prefetch_hint_t* hint);
int solver_codekd_page_plan_seal_union(
    solver_codekd_page_workspace_t* workspace,
    anbool include_descriptor,
    size_t* unique_pages_out,
    size_t* range_count_out,
    size_t* aligned_bytes_out);
void solver_codekd_page_plan_record_refusal(
    solver_codekd_search_packet_t* packet,
    solver_codekd_page_plan_reason_t reason);
int solver_codekd_search_packet_prepare_next_plan(
    solver_codekd_search_packet_t* packet,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque);
anbool solver_codekd_packet_plan_cancelled(void* opaque);
void solver_codekd_packet_disable_verification_delivery(
    solver_codekd_search_packet_t* packet);
void solver_codekd_packet_reset_verify_plan(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_plan_verification_pages(
    void* opaque,
    fitsbin_payload_io_cancel_check_fn cancelled,
    void* cancel_opaque,
    fitsbin_prefetch_range_t* ranges,
    size_t range_capacity,
    size_t* range_count);
void solver_codekd_packet_clear_sweep_storage(
    solver_codekd_search_packet_t* packet);
void solver_codekd_record_clear_prepared_verification(
    solver_candidate_delivery_record_t* record);
void solver_codekd_record_clear_verification_speculation(
    solver_candidate_delivery_record_t* record);
int solver_codekd_search_packet_cleanup(
    solver_codekd_search_packet_t* packet);
int solver_codekd_search_packet_release_ticket(
    solver_codekd_search_packet_t* packet);
void solver_codekd_packet_profile_accumulate(
    solver_t* solver,
    const solver_codekd_search_packet_t* packet);

int solver_codekd_search_packet_prepare(
    const solver_t* solver,
    solver_codekd_search_packet_t* packet,
    solver_ab_descriptor_output_t* descriptors,
    const kdtree_t* tree,
    const quadfile_t* quads,
    startree_t* starkd,
    int dimquads,
    anbool use_radec,
    size_t candidate_budget_bytes,
    unsigned long long sequence,
    anbool detailed);
index_shard_helper_task_status_t
solver_codekd_packet_generate_descriptors(
    const void* input_bytes,
    size_t input_size,
    void* output_bytes,
    size_t output_size);
int solver_codekd_packet_begin_candidate_window(
    solver_codekd_search_packet_t* packet);
void solver_codekd_packet_clear_candidate_window(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_submit_pages(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_submit_candidate_pages(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_submit_verification_pages(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_submit_sweep_pages(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_collect_pages(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_cancel_pages(
    solver_codekd_search_packet_t* packet);
index_shard_helper_task_status_t solver_codekd_packet_decode_quad_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work);
index_shard_helper_task_status_t solver_codekd_packet_copy_star_ready(
    const solver_codekd_packet_task_input_t* input,
    solver_codekd_search_packet_t* packet,
    anbool* more_work);

index_shard_helper_task_status_t
solver_codekd_packet_query_verification_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work);
index_shard_helper_task_status_t
solver_codekd_packet_capture_sweep_ready(
    solver_codekd_search_packet_t* packet,
    anbool* more_work);
int solver_codekd_packet_prepare_verification_owner(
    const solver_codekd_packet_task_input_t* input,
    solver_codekd_search_packet_t* packet);
index_shard_helper_task_status_t
solver_codekd_packet_score_verification_ready(
    solver_codekd_search_packet_t* packet);
int solver_codekd_packet_retired_verification_complete(
    solver_codekd_search_packet_t* packet);

int solver_codekd_descriptor_execute_owner(
    const solver_ab_descriptor_output_t* output,
    solver_t* solver,
    int dimquads,
    kdtree_qres_t** query_result,
    unsigned long long* reduced);
void solver_codekd_packet_begin_result_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    const solver_codekd_result_slot_t* slot,
    const kdtree_qres_t* result,
    int dimquad);
void solver_codekd_packet_resolve_result_range_owner(
    solver_t* solver,
    const solver_ab_descriptor_t* descriptor,
    kdtree_qres_t* result,
    int dimquad,
    int candidate_first,
    int candidate_end,
    solver_candidate_delivery_record_t* prepared,
    size_t prepared_first,
    size_t prepared_quad_count,
    size_t prepared_star_count);
extern const index_shard_staged_ops_t solver_codekd_packet_staged_ops;

index_shard_staged_retire_status_t solver_codekd_search_packet_retire(
    const index_shard_staged_task_t* task,
    size_t task_index,
    void* owner_context);

#endif
