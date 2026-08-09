/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include "test_fitsbin_payload_common.h"

void test_fitsbin_payload_short_read_and_eintr(CuTest* ct);
void test_fitsbin_payload_eof_is_eio(CuTest* ct);
void test_fitsbin_mapped_population_is_all_or_nothing(CuTest* ct);
void test_fitsbin_payload_shared_reader_credit(CuTest* ct);
void test_fitsbin_payload_dynamic_service_width(CuTest* ct);
void test_fitsbin_payload_async_direct_overlap(CuTest* ct);
void test_fitsbin_payload_async_mapped_population(CuTest* ct);
void test_fitsbin_payload_deferred_mapped_plan(CuTest* ct);
void test_fitsbin_payload_fully_resident_plan_is_advisory(CuTest* ct);
void test_fitsbin_payload_mapped_source_parallel_lanes(CuTest* ct);
void test_fitsbin_payload_mapped_completion_precedes_next_plan(CuTest* ct);
void test_fitsbin_payload_precomputed_plan_refresh(CuTest* ct);
void test_fitsbin_payload_precomputed_plan_expands_after_reuse_expiry(
    CuTest* ct);
void test_fitsbin_payload_queue_gap_coalescing(CuTest* ct);
void test_fitsbin_payload_exact_range_order(CuTest* ct);
void test_fitsbin_payload_async_direct_destination(CuTest* ct);
void test_fitsbin_payload_async_mapping_boundary(CuTest* ct);
void test_fitsbin_payload_poll_transfers_ticket_ownership(CuTest* ct);
void test_fitsbin_payload_poll_failure_transfers_once(CuTest* ct);
void test_fitsbin_payload_poll_cancel_transfers_once(CuTest* ct);
void test_fitsbin_payload_drain_registered_waiter(CuTest* ct);
