/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef ONEFIELD_INTERNAL_H
#define ONEFIELD_INTERNAL_H

#include <stddef.h>
#include <sys/stat.h>
#include <time.h>

#include "index_shard_internal.h"
#include "onefield.h"

/* Filename-loaded indexes transfer to one task and are released once. */
size_t onefield_internal_index_count(const onefield_t* bp);
index_t* onefield_internal_get_index(onefield_t* bp, size_t index_order);
char* onefield_internal_get_index_name(onefield_t* bp, size_t index_order);
int onefield_internal_done_with_index(onefield_t* bp,
                                      size_t index_order,
                                      index_t* index);

index_t* onefield_internal_job_index_cache_get(onefield_t* bp,
                                               const char* configured_path);
void onefield_internal_job_index_cache_prepare(
    onefield_t* bp,
    const char* configured_path);
void onefield_internal_job_index_cache_flush(onefield_t* bp);

/*
 * The master owns cached field storage. Worker views borrow it only while the
 * enclosing shard pass is live; the opaque cache is never exposed here.
 */
int onefield_internal_open_master_xyls(onefield_t* bp);
int onefield_internal_validate_single_field_list(onefield_t* bp);
anbool onefield_internal_same_source_identity(
    const struct stat* first,
    const struct stat* second);
int onefield_internal_prepare_field_view(onefield_t* bp,
                                         int fieldnum,
                                         double* field_read_seconds,
                                         double* preprocess_seconds);
void onefield_internal_reset_field_pass_state(onefield_t* bp);
anbool onefield_internal_field_cache_valid(const onefield_t* bp);
anbool onefield_internal_field_cache_has_field(const onefield_t* bp,
                                               int fieldnum);
anbool onefield_field_cache_key_matches(
    const onefield_t* bp,
    int fieldnum,
    const struct stat* source_stat);

int onefield_index_shard_prepare_job_field_for_run(onefield_t* bp);
index_shard_solve_status_t onefield_index_shard_solve(onefield_t* bp,
                                                      solver_t* solver,
                                                      size_t index_count);

/* Owner-only callbacks used by worker-local onefield adapters and reduction. */
anbool onefield_internal_record_match_callback(MatchObj* mo,
                                               void* userdata);
time_t onefield_internal_timer_callback(void* user_data);
void onefield_internal_solved_field(onefield_t* bp, int fieldnum);
int onefield_internal_compare_matchobjs(const void* v1, const void* v2);
void onefield_internal_remove_invalid_fields(il* fieldlist, int maxfield);

#endif
