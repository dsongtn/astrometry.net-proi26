#!/usr/bin/env bash
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

solver_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" ||
    die "cannot resolve the solver directory"
source_root="$(cd "$solver_dir/.." && pwd)" ||
    die "cannot resolve the source root"
field="${1:-$source_root/demo/apod4.xyls}"
winner_index="${2:-$solver_dir/index-9918.fits}"
nonwinner_index="${3:-$source_root/demo/index-4119.fits}"
permuted_index="${4:-}"
output_dir="${SOLVER_TEST_OUTPUT_DIR:-}"
temporary_output=0
output_claim=""
claim_acquired=0

if [[ -z "$permuted_index" ]]; then
    printf '%s\n' \
        "a genuine permuted-StarKD fixture is required as argument four" >&2
    exit 2
fi

if [[ -n "$output_dir" ]]; then
    if [[ "$output_dir" != /* ]]; then
        printf '%s\n' \
            "SOLVER_TEST_OUTPUT_DIR must be absolute" >&2
        exit 2
    fi
    mkdir -p "$output_dir" ||
        die "cannot create SOLVER_TEST_OUTPUT_DIR"
    output_claim="$output_dir/.solver-test-validation-claim"
    if ! mkdir "$output_claim" 2>/dev/null; then
        printf '%s\n' \
            "SOLVER_TEST_OUTPUT_DIR is already claimed" >&2
        exit 2
    fi
    claim_acquired=1
    unexpected_path="$(
        find "$output_dir" \
            -mindepth 1 \
            -maxdepth 1 \
            ! -name '.solver-test-validation-claim' \
            -print -quit
    )" || die "cannot inspect SOLVER_TEST_OUTPUT_DIR"
    if [[ -n "$unexpected_path" ]]; then
        rmdir "$output_claim" ||
            die "cannot release the output-directory claim"
        output_claim=""
        claim_acquired=0
        printf '%s\n' \
            "SOLVER_TEST_OUTPUT_DIR must be empty" >&2
        exit 2
    fi
    work_dir="$output_dir"
else
    work_dir="$(mktemp -d)" ||
        die "cannot create a temporary validation directory"
    temporary_output=1
fi

cleanup() {
    local original_status=$?
    local cleanup_failed=0

    if [[ "$temporary_output" -eq 1 ]]; then
        if ! rm -rf -- "$work_dir"; then
            printf 'ERROR: cannot remove temporary validation directory: %s\n' \
                "$work_dir" >&2
            cleanup_failed=1
        fi
    elif [[ "$claim_acquired" -eq 1 ]]; then
        if ! rmdir "$output_claim" 2>/dev/null; then
            printf 'ERROR: cannot release validation claim: %s\n' \
                "$output_claim" >&2
            cleanup_failed=1
        fi
    fi

    trap - EXIT
    if [[ "$original_status" -ne 0 ]]; then
        exit "$original_status"
    fi
    if [[ "$cleanup_failed" -ne 0 ]]; then
        exit 1
    fi
}
trap cleanup EXIT

verify_default_fixture() {
    local path="$1"
    local default_path="$2"
    local expected="$3"
    local label="$4"
    local actual

    if [[ "$path" != "$default_path" ]]; then
        return
    fi
    actual="$(sha256sum "$path")" ||
        die "cannot hash $label fixture: $path"
    actual="${actual%% *}"
    if [[ "$actual" != "$expected" ]]; then
        printf '%s fixture hash mismatch\nexpected: %s\nactual:   %s\n' \
            "$label" "$expected" "$actual" >&2
        exit 1
    fi
}

verify_default_fixture \
    "$field" \
    "$source_root/demo/apod4.xyls" \
    "b2b380360fa7fc910a796f63a8d161aede7085190b8b45ad0c58e9ff1162570d" \
    "APOD4 field"
verify_default_fixture \
    "$winner_index" \
    "$solver_dir/index-9918.fits" \
    "69d89681a8f65618d999eb4a1feaefc29374508b72f018e7a60f6cb12d1cac77" \
    "winner index"
verify_default_fixture \
    "$nonwinner_index" \
    "$source_root/demo/index-4119.fits" \
    "79b36eea45b72448c8471f6e8af0e1c8635a3821d04f1a23d6cf5ecd5f59d31c" \
    "nonwinner index"

make -C "$source_root/util" -j2 \
    test_fitsbin_payload_io \
    >"$work_dir/payload-io-build.log" 2>&1 ||
    {
        tail -80 "$work_dir/payload-io-build.log" >&2
        die "payload I/O unit binary did not build"
    }
"$source_root/util/test_fitsbin_payload_io" \
    >"$work_dir/payload-io-unit.log" 2>&1 ||
    {
        tail -80 "$work_dir/payload-io-unit.log" >&2
        die "payload I/O unit test failed"
    }

make -C "$solver_dir" -j2 \
    test-solver \
    test-solver-parallel-integration \
    test-solver-streaming-integration \
    test-solver-allocation-failure \
    test-solver-geometry-fallback \
    >"$work_dir/build.log" 2>&1 ||
    {
        tail -80 "$work_dir/build.log" >&2
        die "integration binaries did not build"
    }
mkdir "$work_dir/bin" ||
    die "cannot create the validation binary directory"
install -m 755 \
    "$solver_dir/test-solver" \
    "$solver_dir/test-solver-parallel-integration" \
    "$solver_dir/test-solver-streaming-integration" \
    "$solver_dir/test-solver-allocation-failure" \
    "$solver_dir/test-solver-geometry-fallback" \
    "$work_dir/bin/" ||
    die "cannot install validation binaries"

private_fixture_dir="$work_dir/private-fixtures"
mkdir -m 700 "$private_fixture_dir" ||
    die "cannot create the private fixture directory"
field_touch_w1="$private_fixture_dir/field-touch-w1.xyls"
field_touch_w4="$private_fixture_dir/field-touch-w4.xyls"
index_touch_w1="$private_fixture_dir/index-touch-w1.fits"
index_touch_w4="$private_fixture_dir/index-touch-w4.fits"
index_nohit_w4="$private_fixture_dir/index-nohit-w4.fits"
install -m 600 "$field" "$field_touch_w1" ||
    die "cannot copy the W1 field-mutation fixture"
install -m 600 "$field" "$field_touch_w4" ||
    die "cannot copy the W4 field-mutation fixture"
install -m 600 "$permuted_index" "$index_touch_w1" ||
    die "cannot copy the W1 index-mutation fixture"
install -m 600 "$permuted_index" "$index_touch_w4" ||
    die "cannot copy the W4 index-mutation fixture"
install -m 600 "$permuted_index" "$index_nohit_w4" ||
    die "cannot copy the no-hit permuted-index fixture"

"$work_dir/bin/test-solver" \
    >"$work_dir/solver-unit.log" 2>&1 ||
    {
        tail -80 "$work_dir/solver-unit.log" >&2
        die "solver unit test failed"
    }

run_case() {
    local name="$1"
    local binary="$2"
    local workers="$3"
    local first_object="$4"
    local last_object="$5"
    local mode="$6"
    shift 6
    local status
    local case_field="${SOLVER_TEST_CASE_FIELD:-$field}"

    timeout 45 \
        "$work_dir/bin/$binary" \
        "$case_field" \
        "$workers" \
        "$first_object" \
        "$last_object" \
        "$work_dir/$name.wcs" \
        "$mode" \
        "$@" \
        >"$work_dir/$name.log" 2>&1
    status=$?
    if [[ "$status" -ne 0 ]]; then
        printf 'integration case %s failed with status %i\n' \
            "$name" "$status" >&2
        tail -60 "$work_dir/$name.log" >&2
        die "integration case $name failed"
    fi
}

profile_key() {
    awk '
        /\[solver\] phase-profile/ {
            line = $0
        }
        END {
            if (line == "") {
                printf "<missing-profile:%s>", FILENAME
                exit
            }
            count = split(line, fields)
            found = 0
            for (i = 1; i <= count; i++) {
                if (fields[i] ~ /^(codekd_calls|codekd_hits|resolve_calls|verify_calls|hypothesis_order|kd_result_order|candidate_order)=/) {
                    printf "%s ", fields[i]
                    found++
                }
            }
            if (found != 7) {
                printf "<incomplete-profile:%s:%i>", FILENAME, found
            }
        }
    ' "$1"
}

result_key() {
    awk '
        /^SOLVER_TEST_RESULT / {
            line = $0
        }
        END {
            if (line == "") {
                printf "<missing-result:%s>", FILENAME
                exit
            }
            sub(/workers=[0-9]+ /, "", line)
            print line
        }
    ' "$1"
}

wcs_key() {
    awk '
        /^SOLVER_TEST_WCS / {
            line = $0
        }
        END {
            if (line == "") {
                printf "<missing-wcs:%s>", FILENAME
                exit
            }
            sub(/workers=[0-9]+ /, "", line)
            print line
        }
    ' "$1"
}

profile_sequence_key() {
    awk '
        /\[solver\] phase-profile/ {
            count = split($0, fields)
            found = 0
            for (i = 1; i <= count; i++) {
                if (fields[i] ~ /^(codekd_calls|codekd_hits|resolve_calls|verify_calls|hypothesis_order|kd_result_order|candidate_order)=/) {
                    printf "%s ", fields[i]
                    found++
                }
            }
            if (found != 7) {
                printf "<incomplete-profile:%s:%i>", FILENAME, found
            }
            printf "\n"
            rows++
        }
        END {
            if (!rows) {
                printf "<missing-profile-sequence:%s>\n", FILENAME
            }
        }
    ' "$1"
}

result_sequence_key() {
    awk '
        /^SOLVER_TEST_RESULT / {
            line = $0
            sub(/workers=[0-9]+ /, "", line)
            print line
            rows++
        }
        END {
            if (!rows) {
                printf "<missing-result-sequence:%s>\n", FILENAME
            }
        }
    ' "$1"
}

wcs_sequence_key() {
    awk '
        /^SOLVER_TEST_WCS / {
            line = $0
            sub(/workers=[0-9]+ /, "", line)
            print line
            rows++
        }
        END {
            if (!rows) {
                printf "<missing-wcs-sequence:%s>\n", FILENAME
            }
        }
    ' "$1"
}

match_set_key() {
    awk '
        /^SOLVER_TEST_MATCH / {
            line = $0
            sub(/^SOLVER_TEST_MATCH mode=[^ ]+ pass=[0-9]+ workers=[0-9]+ first_object=[0-9]+ last_object=[0-9]+ /, "", line)
            print line
            rows++
        }
        END {
            if (!rows) {
                printf "<missing-match-set:%s>\n", FILENAME
            }
        }
    ' "$1" |
        LC_ALL=C sort
}

assert_equal() {
    local label="$1"
    local expected="$2"
    local actual="$3"

    if [[ "$expected" != "$actual" ]]; then
        printf '%s mismatch\nexpected: %s\nactual:   %s\n' \
            "$label" "$expected" "$actual" >&2
        die "$label differs"
    fi
}

assert_contains() {
    local label="$1"
    local pattern="$2"
    local file="$3"

    if ! grep -qE "$pattern" "$file"; then
        printf '%s: missing pattern %s in %s\n' \
            "$label" "$pattern" "$file" >&2
        tail -60 "$file" >&2
        die "$label is missing"
    fi
}

assert_not_contains() {
    local label="$1"
    local pattern="$2"
    local file="$3"

    if grep -qE "$pattern" "$file"; then
        printf '%s: unexpected pattern %s in %s\n' \
            "$label" "$pattern" "$file" >&2
        tail -60 "$file" >&2
        die "$label is present"
    fi
}

source "$solver_dir/check-solver-integration-common.sh"
compare_single_index_case() {
    local label="$1"
    local first_object="$2"
    local last_object="$3"
    local w1="${label}_w1"
    local w4="${label}_w4"

    run_case \
        "$w1" \
        test-solver-parallel-integration \
        1 \
        "$first_object" \
        "$last_object" \
        probe \
        "$winner_index"
    run_case \
        "$w4" \
        test-solver-parallel-integration \
        4 \
        "$first_object" \
        "$last_object" \
        probe \
        "$winner_index"
    assert_equal \
        "$label profile" \
        "$(profile_key "$work_dir/$w1.log")" \
        "$(profile_key "$work_dir/$w4.log")"
    assert_equal \
        "$label result" \
        "$(result_key "$work_dir/$w1.log")" \
        "$(result_key "$work_dir/$w4.log")"
    assert_equal \
        "$label WCS" \
        "$(wcs_key "$work_dir/$w1.log")" \
        "$(wcs_key "$work_dir/$w4.log")"
    assert_assist_lifecycle \
        "$label dynamic lending" \
        "$work_dir/$w4.log" \
        1 \
        1
}

timeout 45 \
    "$work_dir/bin/test-solver-parallel-integration" \
    --canonical-index-order \
    >"$work_dir/canonical-index-order.log" 2>&1 ||
    {
        tail -60 "$work_dir/canonical-index-order.log" >&2
        die "canonical index-order regression failed"
    }
assert_contains \
    "canonical later-completion winner" \
    '^INDEX_SHARD_CANONICAL_ORDER_OK completion=1,0 analyzed=1,1 merged=1 reported=1 loser_freed=1$' \
    "$work_dir/canonical-index-order.log"
assert_contains \
    "canonical ordered reduction" \
    '\[index-shard\] pass-detail candidates=2 reduced=1 .*rc=0 status=0 master_committed=1 winner_selected=1 ' \
    "$work_dir/canonical-index-order.log"
assert_not_contains \
    "canonical losing-index reduction" \
    '\[index-shard\] reduce index_order=1 ' \
    "$work_dir/canonical-index-order.log"

# Selection is invariant across representative frontier ranges.  Actual scalar
# versus assisted execution remains a live pair-count/weight planner decision.
compare_single_index_case range_1_8 1 8
compare_single_index_case range_1_10 1 10
compare_single_index_case range_1_11 1 11
compare_single_index_case range_5_15 5 15
compare_single_index_case range_11_20 11 20

permuted_hash="$(sha256sum "$permuted_index")" ||
    die "cannot hash the permuted StarKD fixture"
permuted_hash="${permuted_hash%% *}"
if [[ "$permuted_hash" != \
    "04b9e86038eb99fe4037738b497caa14d7777c6553b1f631dbc05bd8cc2be73b" ]]; then
    printf '%s\nexpected: %s\nactual:   %s\n' \
        "permuted StarKD fixture hash mismatch" \
        "04b9e86038eb99fe4037738b497caa14d7777c6553b1f631dbc05bd8cc2be73b" \
        "$permuted_hash" >&2
    exit 1
fi

run_case \
    permuted_w1 \
    test-solver-parallel-integration \
    1 1 20 probe \
    "$permuted_index"
run_case \
    permuted_w4 \
    test-solver-parallel-integration \
    4 1 20 probe \
    "$permuted_index"
assert_equal \
    "permuted StarKD profile" \
    "$(profile_key "$work_dir/permuted_w1.log")" \
    "$(profile_key "$work_dir/permuted_w4.log")"
assert_equal \
    "permuted StarKD result" \
    "$(result_key "$work_dir/permuted_w1.log")" \
    "$(result_key "$work_dir/permuted_w4.log")"
assert_equal \
    "permuted StarKD WCS" \
    "$(wcs_key "$work_dir/permuted_w1.log")" \
    "$(wcs_key "$work_dir/permuted_w4.log")"
assert_assist_lifecycle \
    "permuted StarKD dynamic lending" \
    "$work_dir/permuted_w4.log" \
    1 \
    1

# Reopen the same stable permuted index in a second fixed-pool generation.
# The native initializer still owns each StarKD acquisition bracket, while
# the job-scoped inverse cache must admit generation one, transfer on
# generation two, and end with no active lease.
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    permuted_cache_w4 \
    test-solver-parallel-integration \
    4 1 20 multipass \
    "$permuted_index"
assert_assist_lifecycle \
    "permuted StarKD cache multipass" \
    "$work_dir/permuted_cache_w4.log" \
    2 \
    2
if [[ "$(line_count \
        '\[index-shard\] inverse-cache state=admit ' \
        "$work_dir/permuted_cache_w4.log")" -ne 1 ||
    "$(line_count \
        '\[index-shard\] inverse-cache state=hit ' \
        "$work_dir/permuted_cache_w4.log")" -ne 1 ]]; then
    die "permuted StarKD cache did not complete one admit-to-hit transfer"
fi
assert_contains \
    "permuted StarKD cache accounting" \
    '\[index-shard\] inverse-cache hits=1 misses=1 admitted=1 .* active=0 ' \
    "$work_dir/permuted_cache_w4.log"
assert_contains \
    "stable index mapping release" \
    '\[index-shard\] job-index-cache state=end hits=0 misses=0 admitted=0 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* budget=0 entries=0 ' \
    "$work_dir/permuted_cache_w4.log"

# Mutate only private index copies between identical passes. A retained
# inverse from generation one must not match the new source identity.
SOLVER_TEST_TOUCH_INDEX_BEFORE_SECOND="$index_touch_w1" \
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    index_touch_w1 \
    test-solver-parallel-integration \
    1 1 20 multipass \
    "$index_touch_w1"
SOLVER_TEST_TOUCH_INDEX_BEFORE_SECOND="$index_touch_w4" \
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    index_touch_w4 \
    test-solver-parallel-integration \
    4 1 20 multipass \
    "$index_touch_w4"
assert_equal \
    "index-identity mutation profile sequence" \
    "$(profile_sequence_key "$work_dir/index_touch_w1.log")" \
    "$(profile_sequence_key "$work_dir/index_touch_w4.log")"
assert_equal \
    "index-identity mutation result sequence" \
    "$(result_sequence_key "$work_dir/index_touch_w1.log")" \
    "$(result_sequence_key "$work_dir/index_touch_w4.log")"
assert_equal \
    "index-identity mutation WCS sequence" \
    "$(wcs_sequence_key "$work_dir/index_touch_w1.log")" \
    "$(wcs_sequence_key "$work_dir/index_touch_w4.log")"
assert_equal \
    "index-identity mutation match set" \
    "$(match_set_key "$work_dir/index_touch_w1.log")" \
    "$(match_set_key "$work_dir/index_touch_w4.log")"
assert_assist_lifecycle \
    "index-identity mutation dynamic lending" \
    "$work_dir/index_touch_w4.log" \
    2 \
    2
if [[ "$(line_count \
        '\[index-shard\] inverse-cache state=hit ' \
        "$work_dir/index_touch_w4.log")" -ne 0 ||
    "$(line_count \
        '\[index-shard\] inverse-cache state=admit ' \
        "$work_dir/index_touch_w4.log")" -ne 2 ]]; then
    die "touched permuted index reused stale inverse state"
fi
assert_contains \
    "touched permuted-index cache accounting" \
    '\[index-shard\] inverse-cache hits=0 misses=2 admitted=2 .* active=0 ' \
    "$work_dir/index_touch_w4.log"
assert_contains \
    "W1 touched index mapping release" \
    '\[index-shard\] job-index-cache state=end hits=0 misses=0 admitted=0 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* budget=0 entries=0 ' \
    "$work_dir/index_touch_w1.log"
assert_contains \
    "W4 touched index mapping release" \
    '\[index-shard\] job-index-cache state=end hits=0 misses=0 admitted=0 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* budget=0 entries=0 ' \
    "$work_dir/index_touch_w4.log"

# Merely opening a permuted StarKD must not initialize or pin its inverse.
# Objects 1--2 contain no real CodeKD hit in this pinned fixture.
run_case \
    permuted_nohit_w4 \
    test-solver-parallel-integration \
    4 1 2 probe \
    "$index_nohit_w4"
assert_contains \
    "permuted no-hit pinned result" \
    '^SOLVER_TEST_RESULT .*indexes=1 solutions=0 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=14650fb0739d0383$' \
    "$work_dir/permuted_nohit_w4.log"
assert_assist_lifecycle \
    "permuted no-hit dynamic lending" \
    "$work_dir/permuted_nohit_w4.log" \
    1 \
    1
assert_contains \
    "permuted no-hit cache accounting" \
    '\[index-shard\] inverse-cache hits=0 misses=0 admitted=0 .* active=0 ' \
    "$work_dir/permuted_nohit_w4.log"

# The production engine reuses one fixed pool across configured frontier
# passes. Exercise generation reset in both ascending and deliberately
# reordered range sequences; range order is user configuration, not a
# scheduler definition.
SOLVER_TEST_SECOND_FIRST_OBJECT=11 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    multipass_forward_w1 \
    test-solver-parallel-integration \
    1 1 10 multipass \
    "$winner_index"
SOLVER_TEST_SECOND_FIRST_OBJECT=11 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    multipass_forward_w4 \
    test-solver-parallel-integration \
    4 1 10 multipass \
    "$winner_index"
assert_equal \
    "forward multipass profile sequence" \
    "$(profile_sequence_key "$work_dir/multipass_forward_w1.log")" \
    "$(profile_sequence_key "$work_dir/multipass_forward_w4.log")"
assert_equal \
    "forward multipass result sequence" \
    "$(result_sequence_key "$work_dir/multipass_forward_w1.log")" \
    "$(result_sequence_key "$work_dir/multipass_forward_w4.log")"
assert_equal \
    "forward multipass WCS sequence" \
    "$(wcs_sequence_key "$work_dir/multipass_forward_w1.log")" \
    "$(wcs_sequence_key "$work_dir/multipass_forward_w4.log")"
assert_assist_lifecycle \
    "forward multipass dynamic lending" \
    "$work_dir/multipass_forward_w4.log" \
    2 \
    2
if [[ "$(line_count \
        '^\[index-shard\] workers=4 mode=pthread ' \
        "$work_dir/multipass_forward_w4.log")" -ne 1 ||
    "$(line_count \
        '\[index-shard\] solver-pass generation=1 ' \
        "$work_dir/multipass_forward_w4.log")" -ne 1 ||
    "$(line_count \
        '\[index-shard\] solver-pass generation=2 ' \
        "$work_dir/multipass_forward_w4.log")" -ne 1 ||
    "$(line_count \
        '\[index-shard\] pthread-pool stop' \
        "$work_dir/multipass_forward_w4.log")" -ne 1 ]]; then
    die "forward multipass fixed-pool lifecycle is not exactly start-submit-submit-stop"
fi
assert_contains \
    "ordinary two-pass field cache" \
    '\[index-shard\] job-field-cache state=end reads=1 preprocesses=1 hits=1 invalidations=0$' \
    "$work_dir/multipass_forward_w4.log"
assert_contains \
    "W1 ordinary two-pass index mapping release" \
    '\[index-shard\] job-index-cache state=end hits=0 misses=0 admitted=0 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* budget=0 entries=0 ' \
    "$work_dir/multipass_forward_w1.log"
assert_contains \
    "W4 ordinary two-pass index mapping release" \
    '\[index-shard\] job-index-cache state=end hits=0 misses=0 admitted=0 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* budget=0 entries=0 ' \
    "$work_dir/multipass_forward_w4.log"

# Mutate only private XYLS copies between identical passes. The cache must
# invalidate on source identity, reread and preprocess exactly once, while
# W1 and W4 retain identical persisted science.
SOLVER_TEST_CASE_FIELD="$field_touch_w1" \
SOLVER_TEST_TOUCH_FIELD_BEFORE_SECOND="$field_touch_w1" \
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    field_touch_w1 \
    test-solver-parallel-integration \
    1 1 20 multipass \
    "$winner_index"
SOLVER_TEST_CASE_FIELD="$field_touch_w4" \
SOLVER_TEST_TOUCH_FIELD_BEFORE_SECOND="$field_touch_w4" \
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=20 \
run_case \
    field_touch_w4 \
    test-solver-parallel-integration \
    4 1 20 multipass \
    "$winner_index"
assert_equal \
    "field-identity mutation profile sequence" \
    "$(profile_sequence_key "$work_dir/field_touch_w1.log")" \
    "$(profile_sequence_key "$work_dir/field_touch_w4.log")"
assert_equal \
    "field-identity mutation result sequence" \
    "$(result_sequence_key "$work_dir/field_touch_w1.log")" \
    "$(result_sequence_key "$work_dir/field_touch_w4.log")"
assert_equal \
    "field-identity mutation WCS sequence" \
    "$(wcs_sequence_key "$work_dir/field_touch_w1.log")" \
    "$(wcs_sequence_key "$work_dir/field_touch_w4.log")"
assert_equal \
    "field-identity mutation match set" \
    "$(match_set_key "$work_dir/field_touch_w1.log")" \
    "$(match_set_key "$work_dir/field_touch_w4.log")"
assert_contains \
    "W1 field-identity invalidation accounting" \
    '\[index-shard\] job-field-cache state=end reads=2 preprocesses=2 hits=0 invalidations=1$' \
    "$work_dir/field_touch_w1.log"
assert_contains \
    "W4 field-identity invalidation accounting" \
    '\[index-shard\] job-field-cache state=end reads=2 preprocesses=2 hits=0 invalidations=1$' \
    "$work_dir/field_touch_w4.log"
assert_assist_lifecycle \
    "field-identity mutation dynamic lending" \
    "$work_dir/field_touch_w4.log" \
    2 \
    2

SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=8 \
run_case \
    multipass_reordered_w1 \
    test-solver-parallel-integration \
    1 11 20 multipass \
    "$winner_index"
SOLVER_TEST_SECOND_FIRST_OBJECT=1 \
SOLVER_TEST_SECOND_LAST_OBJECT=8 \
run_case \
    multipass_reordered_w4 \
    test-solver-parallel-integration \
    4 11 20 multipass \
    "$winner_index"
assert_equal \
    "reordered multipass profile sequence" \
    "$(profile_sequence_key "$work_dir/multipass_reordered_w1.log")" \
    "$(profile_sequence_key "$work_dir/multipass_reordered_w4.log")"
assert_equal \
    "reordered multipass result sequence" \
    "$(result_sequence_key "$work_dir/multipass_reordered_w1.log")" \
    "$(result_sequence_key "$work_dir/multipass_reordered_w4.log")"
assert_equal \
    "reordered multipass WCS sequence" \
    "$(wcs_sequence_key "$work_dir/multipass_reordered_w1.log")" \
    "$(wcs_sequence_key "$work_dir/multipass_reordered_w4.log")"
assert_assist_lifecycle \
    "reordered multipass dynamic lending" \
    "$work_dir/multipass_reordered_w4.log" \
    2 \
    2

run_case \
    regular_w1_exhaustive \
    test-solver-parallel-integration \
    1 1 20 exhaustive \
    "$winner_index"
run_case \
    regular_w4_exhaustive \
    test-solver-parallel-integration \
    4 1 20 exhaustive \
    "$winner_index"
run_case \
    stream_w4_exhaustive \
    test-solver-streaming-integration \
    4 1 20 exhaustive \
    "$winner_index"

assert_equal \
    "regular exhaustive profile" \
    "$(profile_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(profile_key "$work_dir/regular_w4_exhaustive.log")"
assert_equal \
    "streaming exhaustive profile" \
    "$(profile_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(profile_key "$work_dir/stream_w4_exhaustive.log")"
assert_equal \
    "regular exhaustive result" \
    "$(result_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(result_key "$work_dir/regular_w4_exhaustive.log")"
assert_equal \
    "streaming exhaustive result" \
    "$(result_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(result_key "$work_dir/stream_w4_exhaustive.log")"
assert_equal \
    "regular exhaustive WCS" \
    "$(wcs_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(wcs_key "$work_dir/regular_w4_exhaustive.log")"
assert_equal \
    "streaming exhaustive WCS" \
    "$(wcs_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(wcs_key "$work_dir/stream_w4_exhaustive.log")"
assert_assist_lifecycle \
    "regular exhaustive dynamic lending" \
    "$work_dir/regular_w4_exhaustive.log" \
    1 \
    1
assert_assist_lifecycle \
    "streaming exhaustive dynamic lending" \
    "$work_dir/stream_w4_exhaustive.log" \
    1 \
    1

# A deterministic bounded CodeKD packet allocation refusal occurs after one
# assisted phase completes and before the next phase mutates owner state. It
# must use the exact native phase path without failing or replaying the
# complete index task.
run_case \
    allocation_failure_recovery \
    test-solver-allocation-failure \
    4 1 20 exhaustive \
    "$winner_index"
assert_equal \
    "allocation-failure recovery profile" \
    "$(profile_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(profile_key "$work_dir/allocation_failure_recovery.log")"
assert_equal \
    "allocation-failure recovery result" \
    "$(result_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(result_key "$work_dir/allocation_failure_recovery.log")"
assert_equal \
    "allocation-failure recovery WCS" \
    "$(wcs_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(wcs_key "$work_dir/allocation_failure_recovery.log")"
assert_assist_lifecycle \
    "allocation-failure dynamic lending" \
    "$work_dir/allocation_failure_recovery.log" \
    1 \
    1
if [[ "$(
        grep -c '\[onefield-profile\] mode=serial-precommit-retry ' \
            "$work_dir/allocation_failure_recovery.log"
    )" -ne 0 ||
    "$(
        grep -c '^SOLVER_TEST_RESULT ' \
            "$work_dir/allocation_failure_recovery.log"
    )" -ne 1 ]]; then
    printf '%s\n' \
        "allocation-failure recovery lifecycle is not unique" >&2
    exit 1
fi
assert_contains \
    "allocation-failure assisted profile" \
    '\[solver\] phase-profile detailed=1 failed=0 .*allocation_failures=1 .*helper_tasks=[1-9][0-9]* ' \
    "$work_dir/allocation_failure_recovery.log"
assert_contains \
    "allocation-failure shard profile" \
    '\[index-shard\] solver-pass generation=1 .*failed=0 .*batch_failed=0 ' \
    "$work_dir/allocation_failure_recovery.log"
assert_equal \
    "continuous versus split match set W1" \
    "$(match_set_key "$work_dir/regular_w1_exhaustive.log")" \
    "$(match_set_key "$work_dir/multipass_forward_w1.log")"
assert_equal \
    "continuous versus split match set W4" \
    "$(match_set_key "$work_dir/regular_w4_exhaustive.log")" \
    "$(match_set_key "$work_dir/multipass_forward_w4.log")"

segments="$(
    awk '
        /\[solver\] page-pipeline/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^windows=/) {
                    split($i, value, "=")
                    result = value[2]
                }
            }
        }
        END {
            print result
        }
    ' "$work_dir/stream_w4_exhaustive.log"
)"
parallel_batches="$(
    awk '
        /\[solver\] phase-profile/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^parallel_batches_observed=/) {
                    split($i, value, "=")
                    result = value[2]
                }
            }
        }
        END {
            print result
        }
    ' "$work_dir/stream_w4_exhaustive.log"
)"
helper_tasks="$(
    awk '
        /\[solver\] phase-profile/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^helper_tasks=/) {
                    split($i, value, "=")
                    result = value[2]
                }
            }
        }
        END {
            print result
        }
    ' "$work_dir/stream_w4_exhaustive.log"
)"
if [[ -z "$segments" || "$segments" -le 1 ]]; then
    printf 'tiny-buffer test did not force multiple delivery windows\n' >&2
    exit 1
fi
if [[ -z "$parallel_batches" || "$parallel_batches" -le 0 ]]; then
    printf 'tiny-buffer test did not observe helper execution\n' >&2
    exit 1
fi
if [[ -z "$helper_tasks" || "$helper_tasks" -le 0 ]]; then
    die "tiny-buffer test did not retire any dynamically lent helper task"
fi

# Two, three, and four configured indexes preserve one owner per index.
# Workers exhaust the outer queue first, then lend bounded tasks to any
# still-live lane. Identical indexes make the result signature an exact,
# order-sensitive parity check without prescribing which lane receives a
# late worker.
run_case \
    multi2_w1 \
    test-solver-parallel-integration \
    1 5 15 probe \
    "$winner_index" "$winner_index"
run_case \
    multi2_w4 \
    test-solver-parallel-integration \
    4 5 15 probe \
    "$winner_index" "$winner_index"
assert_equal \
    "two-index result" \
    "$(result_key "$work_dir/multi2_w1.log")" \
    "$(result_key "$work_dir/multi2_w4.log")"
assert_contains \
    "two-index pinned known answer" \
    '^SOLVER_TEST_RESULT .*indexes=2 solutions=2 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=082e8d6419d6f45f ' \
    "$work_dir/multi2_w1.log"
assert_assist_lifecycle \
    "two-index dynamic lending" \
    "$work_dir/multi2_w4.log" \
    1 \
    2

run_case \
    multi3_w1 \
    test-solver-parallel-integration \
    1 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index"
run_case \
    multi3_w4 \
    test-solver-parallel-integration \
    4 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index"
assert_equal \
    "three-index result" \
    "$(result_key "$work_dir/multi3_w1.log")" \
    "$(result_key "$work_dir/multi3_w4.log")"
assert_contains \
    "three-index pinned known answer" \
    '^SOLVER_TEST_RESULT .*indexes=3 solutions=3 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=2bd3e32be2102808 ' \
    "$work_dir/multi3_w1.log"
assert_assist_lifecycle \
    "three-index dynamic lending" \
    "$work_dir/multi3_w4.log" \
    1 \
    3

# At pool-width index concurrency every worker first owns an index. Tail
# lending becomes legal only after an owner exhausts the outer queue.
run_case \
    multi4_w1 \
    test-solver-parallel-integration \
    1 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index"
run_case \
    multi4_w4 \
    test-solver-parallel-integration \
    4 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index"
assert_equal \
    "four-index result" \
    "$(result_key "$work_dir/multi4_w1.log")" \
    "$(result_key "$work_dir/multi4_w4.log")"
assert_contains \
    "four-index pinned known answer" \
    '^SOLVER_TEST_RESULT .*indexes=4 solutions=4 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=957afd8b0473076b ' \
    "$work_dir/multi4_w1.log"
assert_assist_lifecycle \
    "four-index dynamic lending" \
    "$work_dir/multi4_w4.log" \
    1 \
    4
assert_contains \
    "four-index compact prefix geometry" \
    '\[solver-geometry\] mode=compact-triangular ' \
    "$work_dir/multi4_w4.log"

# The tiny-packet build retains the same eight-index scientific result while
# exercising real assistance only after every current-band index has been
# claimed by an outer owner. More indexes than workers also proves executor
# and owner-query reuse without reducing the production owner width.
run_case \
    owner_credit_w1 \
    test-solver-streaming-integration \
    1 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index" \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index"
run_case \
    owner_credit_w4 \
    test-solver-streaming-integration \
    4 5 15 probe \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index" \
    "$winner_index" "$winner_index" "$winner_index" "$winner_index"
assert_equal \
    "full-owner eight-index result" \
    "$(result_key "$work_dir/owner_credit_w1.log")" \
    "$(result_key "$work_dir/owner_credit_w4.log")"
assert_assist_lifecycle \
    "full-owner tail lending" \
    "$work_dir/owner_credit_w4.log" \
    1 \
    2
assert_full_producer_assistance \
    "full-owner tail lending" \
    "$work_dir/owner_credit_w4.log"
assert_contains \
    "full-owner test admission" \
    'outer_admission=(full-producer|provider-bounded-exact-demand)$' \
    "$work_dir/owner_credit_w4.log"

# A compact-geometry budget refusal is nonfatal. It deliberately falls back
# to the original owner-private pquad loop for that pass; helper publication
# requires an immutable compact snapshot.
run_case \
    geometry_fallback_w4 \
    test-solver-geometry-fallback \
    4 5 15 probe \
    "$winner_index"
assert_equal \
    "geometry refusal result" \
    "$(result_key "$work_dir/range_5_15_w1.log")" \
    "$(result_key "$work_dir/geometry_fallback_w4.log")"
assert_contains \
    "geometry refusal" \
    '\[solver-geometry\] mode=native reason=budget' \
    "$work_dir/geometry_fallback_w4.log"
assert_contains \
    "geometry native-lane refusal" \
    '\[solver\] phase-profile detailed=1 failed=0 .*hypothesis_batches=0 .*helper_tasks=0 ' \
    "$work_dir/geometry_fallback_w4.log"
assert_contains \
    "geometry refusal helper quiescence" \
    '\[index-shard\] helper-pass generation=1 groups=0 completed=0 .*task_failures=0 ' \
    "$work_dir/geometry_fallback_w4.log"

# Geometry-refused owners cannot publish helper lanes. All four queued indexes
# must still be claimed through the ordered outer path without parking behind
# owners that can never publish work.
run_case \
    geometry_fallback_overflow_w4 \
    test-solver-geometry-fallback \
    4 5 15 probe \
    "$nonwinner_index" "$nonwinner_index" \
    "$nonwinner_index" "$nonwinner_index"
[[ "$(grep -cE \
    'claim index_order=[0-9]+ lane=producer ' \
    "$work_dir/geometry_fallback_overflow_w4.log")" -eq 4 ]] ||
    die "geometry refusal did not claim four producer indexes"
# Both accepted outer-scheduler winner directions remain legal.
run_case \
    first_winner_w1 \
    test-solver-parallel-integration \
    1 1 20 winner \
    "$winner_index" "$nonwinner_index"
run_case \
    first_winner_w4 \
    test-solver-parallel-integration \
    4 1 20 winner \
    "$winner_index" "$nonwinner_index"
assert_equal \
    "first-index winner WCS" \
    "$(wcs_key "$work_dir/first_winner_w1.log")" \
    "$(wcs_key "$work_dir/first_winner_w4.log")"
assert_equal \
    "first-index winner result" \
    "$(result_key "$work_dir/first_winner_w1.log")" \
    "$(result_key "$work_dir/first_winner_w4.log")"
assert_contains \
    "first-index winner order" \
    'committed-solution index_order=0 ' \
    "$work_dir/first_winner_w4.log"
assert_assist_lifecycle \
    "first-index winner dynamic lending" \
    "$work_dir/first_winner_w4.log" \
    1 \
    2

run_case \
    later_winner_w1 \
    test-solver-parallel-integration \
    1 1 20 winner \
    "$nonwinner_index" "$winner_index"
run_case \
    later_winner_w4 \
    test-solver-parallel-integration \
    4 1 20 winner \
    "$nonwinner_index" "$winner_index"
assert_equal \
    "later-index winner WCS" \
    "$(wcs_key "$work_dir/later_winner_w1.log")" \
    "$(wcs_key "$work_dir/later_winner_w4.log")"
assert_equal \
    "later-index winner result" \
    "$(result_key "$work_dir/later_winner_w1.log")" \
    "$(result_key "$work_dir/later_winner_w4.log")"
assert_contains \
    "later-index winner order" \
    'committed-solution index_order=1 ' \
    "$work_dir/later_winner_w4.log"
assert_assist_lifecycle \
    "later-index winner dynamic lending" \
    "$work_dir/later_winner_w4.log" \
    1 \
    2

# An already-expired aggregate wall budget exercises pre-bind helper
# quiescence without relying on a timing race.
SOLVER_TEST_TOTAL_WALL_LIMIT=0.000001 \
run_case \
    expired_limit_w4 \
    test-solver-parallel-integration \
    4 1 20 limit \
    "$winner_index" "$winner_index" "$winner_index"
assert_contains \
    "expired wall limit" \
    '^SOLVER_TEST_RESULT .*wall_limit=1 ' \
    "$work_dir/expired_limit_w4.log"

if grep --exclude='allocation_failure_recovery.log' -qE \
    '\[solver-ab\].*(failed|invalid|underflow)|failed=1|serial-precommit-retry|releasing leased job-index-cache' \
    "$work_dir"/*.log; then
    grep --exclude='allocation_failure_recovery.log' -nE \
        '\[solver-ab\].*(failed|invalid|underflow)|failed=1|serial-precommit-retry|releasing leased job-index-cache' \
        "$work_dir"/*.log >&2
    exit 1
fi

printf \
    'SOLVER_PARALLEL_INTEGRATION_OK ranges=5 multipass=2 retained_state=field-hit,field-invalidate,index-release,inverse-invalidate,index-nohit multi=2,3,4 permuted=%s inverse_cache=admit-hit windows=%s parallel_batches=%s helper_tasks=%s output=%s\n' \
    "yes" \
    "$segments" \
    "$parallel_batches" \
    "$helper_tasks" \
    "$([[ "$temporary_output" -eq 1 ]] && printf transient || printf '%s' "$work_dir")" ||
    die "cannot write ordinary integration summary"
