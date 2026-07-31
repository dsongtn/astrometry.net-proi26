#!/usr/bin/env bash

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
output_dir="${AB_VALIDATION_OUTPUT_DIR:-}"
temporary_output=0
claim_dir=
claim_acquired=0

cleanup() {
    local original_status=$?
    local cleanup_failed=0

    if [[ "$claim_acquired" -eq 1 ]]; then
        if ! rmdir -- "$claim_dir" 2>/dev/null; then
            printf 'ERROR: cannot release validation claim: %s\n' \
                "$claim_dir" >&2
            cleanup_failed=1
        fi
    fi
    if [[ "$temporary_output" -eq 1 ]]; then
        if ! rm -rf -- "$output_dir"; then
            printf 'ERROR: cannot remove temporary validation directory: %s\n' \
                "$output_dir" >&2
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

if [[ -z "$permuted_index" ]]; then
    die "a genuine permuted-StarKD fixture is required as argument four"
fi

if [[ -n "$output_dir" ]]; then
    [[ "$output_dir" == /* ]] ||
        die "AB_VALIDATION_OUTPUT_DIR must be absolute"
    mkdir -p "$output_dir" ||
        die "cannot create AB_VALIDATION_OUTPUT_DIR"
else
    output_dir="$(mktemp -d)" ||
        die "cannot create a temporary validation directory"
    temporary_output=1
fi

claim_dir="$output_dir/.ab-validation-claim"
mkdir "$claim_dir" 2>/dev/null ||
    die "AB_VALIDATION_OUTPUT_DIR is already claimed by another runner"
claim_acquired=1
unexpected_path="$(
    find "$output_dir" \
        -mindepth 1 \
        -maxdepth 1 \
        ! -path "$claim_dir" \
        -print -quit
)" || die "cannot inspect AB_VALIDATION_OUTPUT_DIR"
if [[ -n "$unexpected_path" ]]; then
    die "AB_VALIDATION_OUTPUT_DIR must be empty"
fi

verify_hash() {
    local path="$1"
    local expected="$2"
    local label="$3"
    local actual

    [[ -f "$path" ]] || die "$label is not a regular file: $path"
    actual="$(sha256sum "$path")" ||
        die "cannot hash $label: $path"
    actual="${actual%% *}"
    [[ "$actual" == "$expected" ]] ||
        die "$label hash mismatch: expected $expected, got $actual"
}

verify_hash \
    "$field" \
    "b2b380360fa7fc910a796f63a8d161aede7085190b8b45ad0c58e9ff1162570d" \
    "APOD4 field"
verify_hash \
    "$winner_index" \
    "69d89681a8f65618d999eb4a1feaefc29374508b72f018e7a60f6cb12d1cac77" \
    "winner index"
verify_hash \
    "$nonwinner_index" \
    "79b36eea45b72448c8471f6e8af0e1c8635a3821d04f1a23d6cf5ecd5f59d31c" \
    "nonwinner index"
verify_hash \
    "$permuted_index" \
    "04b9e86038eb99fe4037738b497caa14d7777c6553b1f631dbc05bd8cc2be73b" \
    "permuted StarKD index"

make -C "$solver_dir" -j2 \
    test-ab-block-integration \
    test-ab-stream-asan \
    test-ab-stream-tsan \
    test-ab-allocfail-asan \
    test-ab-allocfail-tsan \
    >"$output_dir/build.log" 2>&1 ||
    {
        tail -80 "$output_dir/build.log" >&2
        die "sanitizer integration binaries did not build"
    }

mkdir "$output_dir/bin" ||
    die "cannot create the sanitizer binary directory"
install -m 755 \
    "$solver_dir/test-ab-block-integration" \
    "$solver_dir/test-ab-stream-asan" \
    "$solver_dir/test-ab-stream-tsan" \
    "$solver_dir/test-ab-allocfail-asan" \
    "$solver_dir/test-ab-allocfail-tsan" \
    "$output_dir/bin/" ||
    die "cannot install sanitizer integration binaries"
install -m 644 \
    "$solver_dir/test-ab-stream-asan.map" \
    "$solver_dir/test-ab-stream-tsan.map" \
    "$solver_dir/test-ab-allocfail-asan.map" \
    "$solver_dir/test-ab-allocfail-tsan.map" \
    "$output_dir/" ||
    die "cannot install sanitizer linker maps"

assert_map() {
    local flavor="$1"
    local stem="$2"
    local archive="$3"
    local map="$output_dir/${stem}.map"
    local solver_member="$4"
    local selected
    local plain

    selected="$(
        grep -Ec \
            "^${archive}\\.a\\((${solver_member}|onefield_ab_stream_${flavor}|index_shard_ab_stream_${flavor}|index_shard_config_ab_stream_${flavor}|starkd_ab_stream_${flavor}|fitsbin_ab_stream_${flavor})\\.o\\)" \
            "$map"
    )"
    plain="$(
        grep -Ec \
            "^${archive}\\.a\\((solver|onefield|index_shard|index_shard_config|starkd|fitsbin)\\.o\\)" \
            "$map"
    )"
    [[ "$selected" -eq 6 ]] ||
        die "$stem map selected $selected sanitizer objects, expected 6"
    [[ "$plain" -eq 0 ]] ||
        die "$stem map selected $plain conflicting plain objects"
}

assert_map \
    asan \
    test-ab-stream-asan \
    libastrometry-ab-stream-asan \
    solver_ab_stream_asan
assert_map \
    tsan \
    test-ab-stream-tsan \
    libastrometry-ab-stream-tsan \
    solver_ab_stream_tsan
assert_map \
    asan \
    test-ab-allocfail-asan \
    libastrometry-ab-allocfail-asan \
    solver_ab_allocfail_asan
assert_map \
    tsan \
    test-ab-allocfail-tsan \
    libastrometry-ab-allocfail-tsan \
    solver_ab_allocfail_tsan

for object in solver onefield index_shard index_shard_config starkd fitsbin; do
    nm -u "$solver_dir/${object}_ab_stream_asan.o" |
        grep -Eq '__(asan|ubsan)_' ||
        die "$object ASan/UBSan object has no sanitizer references"
    nm -u "$solver_dir/${object}_ab_stream_tsan.o" |
        grep -q '__tsan_' ||
        die "$object TSan object has no sanitizer references"
done
nm -u "$solver_dir/solver_ab_allocfail_asan.o" |
    grep -Eq '__(asan|ubsan)_' ||
    die "allocation-failure ASan/UBSan solver is not instrumented"
nm -u "$solver_dir/solver_ab_allocfail_tsan.o" |
    grep -q '__tsan_' ||
    die "allocation-failure TSan solver is not instrumented"

run_logged() {
    local name="$1"
    local timeout_seconds="$2"
    shift 2
    local status

    timeout "$timeout_seconds" "$@" \
        >"$output_dir/$name.log" 2>&1
    status=$?
    if [[ "$status" -ne 0 ]]; then
        printf '%s failed with status %i\n' "$name" "$status" >&2
        tail -80 "$output_dir/$name.log" >&2
        exit "$status"
    fi
}

run_cancel_logged() {
    local name="$1"
    local timeout_seconds="$2"
    shift 2
    local log="$output_dir/$name.log"
    local run_pid
    local status
    local polls=0
    local triggered=0

    timeout "$timeout_seconds" "$@" >"$log" 2>&1 &
    run_pid=$!
    while [[ "$polls" -lt 2000 ]]; do
        # Trigger from an observed dynamic loan, not from elapsed time or an
        # empty early-frontier phase. This makes cancellation an actual
        # in-flight lane-quiescence test.
        if grep -Eq \
            '\[index-shard\] assist-lane state=join ' \
            "$log" 2>/dev/null; then
            touch "$output_dir/$name.cancel-request" ||
                die "$name cannot create its cancellation request"
            triggered=1
            break
        fi
        if ! kill -0 "$run_pid" 2>/dev/null; then
            break
        fi
        polls=$((polls + 1))
        sleep 0.01 ||
            die "$name cancellation poll sleep failed"
    done
    if [[ "$triggered" -ne 1 ]]; then
        wait "$run_pid"
        status=$?
        die "$name ended before the log-driven cancellation trigger"
    fi
    wait "$run_pid"
    status=$?
    if [[ "$status" -ne 0 ]]; then
        printf '%s failed with status %i\n' "$name" "$status" >&2
        tail -80 "$log" >&2
        exit "$status"
    fi
}

result_key() {
    awk '
        /^AB_RESULT / {
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

line_count() {
    PATTERN="$1" awk '
        $0 ~ ENVIRON["PATTERN"] {
            count++
        }
        END {
            print count + 0
        }
    ' "$2"
}

assert_assist_lifecycle() {
    local label="$1"
    local file="$2"
    local generations="$3"
    local minimum_publishes="$4"
    local allow_failed_leaves="${5:-0}"
    local submits
    local passes
    local all_passes
    local publishes
    local unpublishes
    local joins
    local leaves
    local bad_leaves
    local loans

    submits="$(line_count \
        'pthread-pool submit .*inner_scheduler=dynamic-pool-lending .*mmap_advice=normal' \
        "$file")"
    passes="$(line_count \
        'assist-pass generation=[0-9]+ .*waiters=0' \
        "$file")"
    all_passes="$(line_count 'assist-pass generation=' "$file")"
    publishes="$(line_count 'assist-lane state=publish ' "$file")"
    unpublishes="$(line_count 'assist-lane state=unpublish ' "$file")"
    joins="$(line_count 'assist-lane state=join ' "$file")"
    leaves="$(line_count 'assist-lane state=leave ' "$file")"
    bad_leaves="$(awk '
        /assist-lane state=leave / && $0 !~ / rc=0$/ {
            count++
        }
        END {
            print count + 0
        }
    ' "$file")"
    loans="$(awk '
        /assist-pass generation=/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^loans=/) {
                    split($i, value, "=")
                    total += value[2]
                }
            }
        }
        END {
            print total + 0
        }
    ' "$file")"

    [[ "$submits" -eq "$generations" ]] ||
        die "$label has $submits dynamic NORMAL submits, expected $generations"
    [[ "$passes" -eq "$generations" ]] ||
        die "$label has $passes clean assist-pass records, expected $generations"
    [[ "$all_passes" -eq "$generations" ]] ||
        die "$label has $all_passes total assist-pass records, expected $generations"
    [[ "$publishes" -ge "$minimum_publishes" ]] ||
        die "$label published $publishes lanes, expected at least $minimum_publishes"
    [[ "$publishes" -eq "$unpublishes" ]] ||
        die "$label lane publish/unpublish count differs ($publishes/$unpublishes)"
    [[ "$joins" -eq "$leaves" ]] ||
        die "$label lane join/leave count differs ($joins/$leaves)"
    [[ "$joins" -eq "$loans" ]] ||
        die "$label helper joins differ from claimed loans ($joins/$loans)"
    if [[ "$allow_failed_leaves" -eq 0 &&
        "$bad_leaves" -ne 0 ]]; then
        die "$label contains $bad_leaves failed helper loans"
    fi
    if grep -qE 'phase-group|static-phase-assist' "$file"; then
        die "$label contains obsolete static phase-group telemetry"
    fi
}

assert_tail_only_lending() {
    local label="$1"
    local file="$2"
    local totals
    local prequeue
    local loans

    totals="$(awk '
        /assist-pass generation=/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^prequeue_loans=/) {
                    split($i, value, "=")
                    prequeue += value[2]
                } else if ($i ~ /^loans=/) {
                    split($i, value, "=")
                    loans += value[2]
                }
            }
        }
        END {
            print prequeue + 0, loans + 0
        }
    ' "$file")"
    read -r prequeue loans <<<"$totals"

    [[ "$prequeue" -eq 0 ]] ||
        die "$label lent $prequeue tasks before outer queue exhaustion"
    [[ "$loans" -gt 0 ]] ||
        die "$label did not exercise tail assistance"
}

assert_executor_quiescence() {
    local label="$1"
    local file="$2"

    awk '
        /\[solver-ab\] executor-epoch / {
            started = 0
            quiesced = 0
            queries = 0
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^started=/) {
                    split($i, value, "=")
                    started = value[2] + 0
                } else if ($i ~ /^quiesced=/) {
                    split($i, value, "=")
                    quiesced = value[2] + 0
                } else if ($i ~ /^owner_query_reuses=/) {
                    split($i, value, "=")
                    queries = value[2] + 0
                }
            }
            if (started <= 0 ||
                quiesced != started ||
                queries <= 0) {
                failed = 1
            }
            found = 1
        }
        END {
            exit found && !failed ? 0 : 1
        }
    ' "$file" ||
        die "$label did not prove clean executor/query quiescence"
}

run_adaptive_limit_logged() {
    local flavor="$1"
    local kind="$2"
    local timeout_seconds="$3"
    local options_name="$4"
    local options_value="$5"
    local binary="$output_dir/bin/test-ab-stream-$flavor"
    local limit_name
    local expected_result
    local attempt=0
    local accepted=0
    local budget
    local log
    local status
    local result_count
    local all_result_count

    case "$kind" in
        wall)
            limit_name=AB_TOTAL_WALL_LIMIT
            expected_result='wall_limit=1 cpu_limit=0 failed=0'
            ;;
        cpu)
            limit_name=AB_TOTAL_CPU_LIMIT
            expected_result='wall_limit=0 cpu_limit=1 failed=0'
            ;;
        *)
            die "unknown adaptive limit kind: $kind"
            ;;
    esac

    : >"$output_dir/${flavor}_${kind}.attempts" ||
        die "cannot create $flavor $kind attempt record"
    for budget in 0.10 0.20 0.40 0.80 1.60 3.20 6.40; do
        attempt=$((attempt + 1))
        log="$output_dir/${flavor}_${kind}_attempt_${attempt}.log"
        timeout "$timeout_seconds" \
            env "$options_name=$options_value" "$limit_name=$budget" \
            "$binary" "$field" 4 11 100 \
            "$output_dir/${flavor}_${kind}_attempt_${attempt}.wcs" \
            limit "$winner_index" \
            >"$log" 2>&1
        status=$?

        if grep -Eq \
            'ERROR: AddressSanitizer|runtime error:|WARNING: ThreadSanitizer|SUMMARY: ThreadSanitizer|data race' \
            "$log"; then
            printf \
                'attempt=%i budget=%s range=11-100 status=%i accepted=0\n' \
                "$attempt" "$budget" "$status" \
                >>"$output_dir/${flavor}_${kind}.attempts" ||
                die "cannot record $flavor $kind sanitizer failure"
            die "$flavor $kind limit attempt $attempt produced a sanitizer diagnostic"
        fi
        if grep -Eq '(^|[[:space:]])failed=1([[:space:]]|$)' "$log"; then
            printf \
                'attempt=%i budget=%s range=11-100 status=%i accepted=0\n' \
                "$attempt" "$budget" "$status" \
                >>"$output_dir/${flavor}_${kind}.attempts" ||
                die "cannot record $flavor $kind solver failure"
            die "$flavor $kind limit attempt $attempt reported solver failure"
        fi
        if [[ "$status" -ne 0 ]]; then
            printf \
                'attempt=%i budget=%s range=11-100 status=%i accepted=0\n' \
                "$attempt" "$budget" "$status" \
                >>"$output_dir/${flavor}_${kind}.attempts" ||
                die "cannot record $flavor $kind process failure"
            printf \
                '%s %s limit attempt %i failed with status %i\n' \
                "$flavor" "$kind" "$attempt" "$status" >&2
            tail -80 "$log" >&2
            exit "$status"
        fi

        result_count="$(
            grep -Ec \
                "^AB_RESULT mode=limit .*${expected_result} " \
                "$log" ||
                true
        )"
        all_result_count="$(line_count '^AB_RESULT ' "$log")"
        if [[ "$result_count" -ne 1 ||
            "$all_result_count" -ne 1 ]]; then
            printf \
                'attempt=%i budget=%s range=11-100 status=%i accepted=0\n' \
                "$attempt" "$budget" "$status" \
                >>"$output_dir/${flavor}_${kind}.attempts" ||
                die "cannot record $flavor $kind invalid result"
            die "$flavor $kind limit attempt $attempt has an incomplete or invalid result footer"
        fi

        if grep -Eq \
            '\[index-shard\] solver-pass generation=1 .*helper_tasks=[1-9][0-9]* ' \
            "$log" &&
            grep -Eq \
                '\[index-shard\] assist-pass generation=1 .*waiters=0' \
                "$log"; then
            accepted=1
            install -m 644 \
                "$log" \
                "$output_dir/${flavor}_${kind}.log" ||
                die "cannot retain selected $flavor $kind log"
            printf \
                'attempt=%i budget=%s range=11-100 status=%i accepted=1\n' \
                "$attempt" "$budget" "$status" \
                >"$output_dir/${flavor}_${kind}.selected" ||
                die "cannot write selected $flavor $kind record"
        fi

        printf \
            'attempt=%i budget=%s range=11-100 status=%i accepted=%i\n' \
            "$attempt" "$budget" "$status" "$accepted" \
            >>"$output_dir/${flavor}_${kind}.attempts" ||
            die "cannot append $flavor $kind attempt record"
        if [[ "$accepted" -eq 1 ]]; then
            break
        fi
    done

    [[ "$accepted" -eq 1 ]] ||
        die "$flavor $kind limit did not find an in-flight assisted-work bracket"
}

run_flavor() {
    local flavor="$1"
    local timeout_seconds="$2"
    local options_name="$3"
    local options_value="$4"
    local binary="$output_dir/bin/test-ab-stream-$flavor"
    local cancel_path="$output_dir/${flavor}_cancel.cancel-request"

    run_logged \
        "${flavor}_canonical_index_order" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" --canonical-index-order
    grep -qx \
        'INDEX_SHARD_CANONICAL_ORDER_OK completion=1,0 analyzed=1,1 merged=0 reported=0 loser_freed=1' \
        "$output_dir/${flavor}_canonical_index_order.log" ||
        die "$flavor canonical index-order regression failed"
    grep -Eq \
        '\[index-shard\] pass-detail candidates=2 reduced=1 .*rc=0 status=0 master_committed=1$' \
        "$output_dir/${flavor}_canonical_index_order.log" ||
        die "$flavor canonical index-order reduction is invalid"
    if grep -q \
        '\[index-shard\] reduce index_order=1 ' \
        "$output_dir/${flavor}_canonical_index_order.log"; then
        die "$flavor canonical losing index was reduced"
    fi

    run_logged \
        "${flavor}_single" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 1 20 \
        "$output_dir/${flavor}_single.wcs" exhaustive \
        "$winner_index"
    run_logged \
        "${flavor}_multipass" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        AB_SECOND_FIRST_OBJECT=11 AB_SECOND_LAST_OBJECT=20 \
        "$binary" "$field" 4 1 10 \
        "$output_dir/${flavor}_multipass.wcs" multipass \
        "$winner_index"
    run_logged \
        "${flavor}_n2" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 5 15 \
        "$output_dir/${flavor}_n2.wcs" probe \
        "$winner_index" "$winner_index"
    run_logged \
        "${flavor}_n3" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 5 15 \
        "$output_dir/${flavor}_n3.wcs" probe \
        "$winner_index" "$winner_index" "$winner_index"
    run_logged \
        "${flavor}_n4" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 5 15 \
        "$output_dir/${flavor}_n4.wcs" probe \
        "$winner_index" "$winner_index" "$winner_index" "$winner_index"
    run_logged \
        "${flavor}_later_winner" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 1 20 \
        "$output_dir/${flavor}_later_winner.wcs" winner \
        "$nonwinner_index" "$winner_index"
    run_logged \
        "${flavor}_permuted" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$binary" "$field" 4 1 20 \
        "$output_dir/${flavor}_permuted.wcs" probe \
        "$permuted_index"
    run_logged \
        "${flavor}_permuted_cache" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        AB_SECOND_FIRST_OBJECT=1 AB_SECOND_LAST_OBJECT=20 \
        "$binary" "$field" 4 1 20 \
        "$output_dir/${flavor}_permuted_cache.wcs" multipass \
        "$permuted_index"

    run_cancel_logged \
        "${flavor}_cancel" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        AB_CANCEL_FILE="$cancel_path" \
        "$binary" "$field" 4 1 100 \
        "$output_dir/${flavor}_cancel.wcs" cancel \
        "$winner_index"

    run_adaptive_limit_logged \
        "$flavor" wall "$timeout_seconds" \
        "$options_name" "$options_value"
    run_adaptive_limit_logged \
        "$flavor" cpu "$timeout_seconds" \
        "$options_name" "$options_value"

    [[ "$(line_count \
        '\[index-shard\] solver-pass generation=' \
        "$output_dir/${flavor}_multipass.log")" -eq 2 ]] ||
        die "$flavor multipass did not complete two pool generations"
    [[ "$(line_count \
        '^AB_RESULT mode=multipass .*failed=0 ' \
        "$output_dir/${flavor}_multipass.log")" -eq 2 ]] ||
        die "$flavor multipass did not publish two clean results"
    assert_assist_lifecycle \
        "$flavor multipass dynamic lending" \
        "$output_dir/${flavor}_multipass.log" \
        2 \
        2
    grep -Eq \
        '\[index-shard\] job-index-cache state=end hits=1 misses=1 admitted=1 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* entries=1$' \
        "$output_dir/${flavor}_multipass.log" ||
        die "$flavor multipass did not retain one coherent index mmap epoch"
    assert_executor_quiescence \
        "$flavor multipass executor quiescence" \
        "$output_dir/${flavor}_multipass.log"
    [[ "$(line_count \
        '^\[index-shard\] workers=4 mode=pthread ' \
        "$output_dir/${flavor}_multipass.log")" -eq 1 &&
        "$(line_count \
            '\[index-shard\] solver-pass generation=1 ' \
            "$output_dir/${flavor}_multipass.log")" -eq 1 &&
        "$(line_count \
            '\[index-shard\] solver-pass generation=2 ' \
            "$output_dir/${flavor}_multipass.log")" -eq 1 &&
        "$(line_count \
            '\[index-shard\] pthread-pool stop' \
            "$output_dir/${flavor}_multipass.log")" -eq 1 ]] ||
        die "$flavor multipass pool lifecycle is not start-submit-submit-stop"
    assert_assist_lifecycle \
        "$flavor single-index dynamic lending" \
        "$output_dir/${flavor}_single.log" \
        1 \
        1
    grep -Eq \
        '^AB_RESULT mode=probe .*indexes=2 solutions=2 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=082e8d6419d6f45f ' \
        "$output_dir/${flavor}_n2.log" ||
        die "$flavor N2 result differs from its pinned fixture answer"
    assert_assist_lifecycle \
        "$flavor N2 dynamic lending" \
        "$output_dir/${flavor}_n2.log" \
        1 \
        2
    grep -Eq \
        '^AB_RESULT mode=probe .*indexes=3 solutions=3 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=2bd3e32be2102808 ' \
        "$output_dir/${flavor}_n3.log" ||
        die "$flavor N3 result differs from its pinned fixture answer"
    assert_assist_lifecycle \
        "$flavor N3 dynamic lending" \
        "$output_dir/${flavor}_n3.log" \
        1 \
        3
    grep -Eq \
        '^AB_RESULT mode=probe .*indexes=4 solutions=4 cancelled=0 wall_limit=0 cpu_limit=0 failed=0 signature=957afd8b0473076b ' \
        "$output_dir/${flavor}_n4.log" ||
        die "$flavor N4 result differs from its pinned fixture answer"
    assert_assist_lifecycle \
        "$flavor N4 dynamic lending" \
        "$output_dir/${flavor}_n4.log" \
        1 \
        4
    assert_tail_only_lending \
        "$flavor N4 tail lending" \
        "$output_dir/${flavor}_n4.log"
    grep -q \
        '\[solver-geometry\] mode=compact-triangular' \
        "$output_dir/${flavor}_n4.log" ||
        die "$flavor N4 outer scheduler did not receive shared prefix geometry"
    grep -q \
        'committed-solution index_order=1 ' \
        "$output_dir/${flavor}_later_winner.log" ||
        die "$flavor later-index winner was not committed"
    assert_assist_lifecycle \
        "$flavor later-index winner dynamic lending" \
        "$output_dir/${flavor}_later_winner.log" \
        1 \
        2
    [[ "$(
        grep -c '\[solver-ab\] starkd-inverse-init state=start ' \
            "$output_dir/${flavor}_permuted.log"
    )" -eq 1 ]] ||
        die "$flavor permuted index did not start exactly one lazy initializer"
    [[ "$(
        grep -c '\[solver-ab\] starkd-inverse-init state=ready ' \
            "$output_dir/${flavor}_permuted.log"
    )" -eq 1 ]] ||
        die "$flavor permuted index did not finish exactly one lazy initializer"
    assert_assist_lifecycle \
        "$flavor permuted dynamic lending" \
        "$output_dir/${flavor}_permuted.log" \
        1 \
        1
    assert_assist_lifecycle \
        "$flavor permuted cache multipass" \
        "$output_dir/${flavor}_permuted_cache.log" \
        2 \
        2
    [[ "$(line_count \
        '\[index-shard\] inverse-cache state=admit ' \
        "$output_dir/${flavor}_permuted_cache.log")" -eq 1 &&
        "$(line_count \
            '\[index-shard\] inverse-cache state=hit ' \
            "$output_dir/${flavor}_permuted_cache.log")" -eq 1 ]] ||
        die "$flavor permuted cache did not complete one admit-to-hit transfer"
    grep -Eq \
        '\[index-shard\] inverse-cache hits=1 misses=1 admitted=1 .* active=0 ' \
        "$output_dir/${flavor}_permuted_cache.log" ||
        die "$flavor permuted cache accounting did not end quiescent"
    grep -Eq \
        '\[index-shard\] job-index-cache state=end hits=1 misses=1 admitted=1 refused=0 invalidated=0 retries=0 fd_close_failures=0 .* entries=1$' \
        "$output_dir/${flavor}_permuted_cache.log" ||
        die "$flavor permuted mmap cache did not end coherent"
    grep -Eq \
        '^AB_RESULT mode=cancel .*cancelled=1 .*failed=0 ' \
        "$output_dir/${flavor}_cancel.log" ||
        die "$flavor active cancellation result is invalid"
    grep -Eq \
        '\[index-shard\] solver-pass generation=1 .*batches=[1-9][0-9]* ' \
        "$output_dir/${flavor}_cancel.log" ||
        die "$flavor cancellation occurred before phase work began"
    assert_assist_lifecycle \
        "$flavor active-cancellation dynamic lending" \
        "$output_dir/${flavor}_cancel.log" \
        1 \
        1
    grep -Eq \
        '^AB_RESULT mode=limit .*wall_limit=1 cpu_limit=0 failed=0 ' \
        "$output_dir/${flavor}_wall.log" ||
        die "$flavor in-flight wall limit was not observed"
    grep -Eq \
        '\[index-shard\] solver-pass generation=1 .*helper_tasks=[1-9][0-9]* ' \
        "$output_dir/${flavor}_wall.log" ||
        die "$flavor wall-limit bracket has no dynamically lent work"
    assert_assist_lifecycle \
        "$flavor wall-limit dynamic lending" \
        "$output_dir/${flavor}_wall.log" \
        1 \
        1
    grep -Eq \
        '^AB_RESULT mode=limit .*wall_limit=0 cpu_limit=1 failed=0 ' \
        "$output_dir/${flavor}_cpu.log" ||
        die "$flavor in-flight CPU limit was not observed"
    grep -Eq \
        '\[index-shard\] solver-pass generation=1 .*helper_tasks=[1-9][0-9]* ' \
        "$output_dir/${flavor}_cpu.log" ||
        die "$flavor CPU-limit bracket has no dynamically lent work"
    assert_assist_lifecycle \
        "$flavor CPU-limit dynamic lending" \
        "$output_dir/${flavor}_cpu.log" \
        1 \
        1
}

asan_options='detect_leaks=0:halt_on_error=1:abort_on_error=1'
tsan_options='halt_on_error=1:abort_on_error=1:history_size=7'
export UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'
run_flavor asan 90 ASAN_OPTIONS "$asan_options"
run_flavor tsan 120 TSAN_OPTIONS "$tsan_options"

run_logged \
    allocfail_baseline \
    90 \
    "$output_dir/bin/test-ab-block-integration" \
    "$field" 1 1 20 "$output_dir/allocfail_baseline.wcs" \
    probe "$winner_index"

for flavor in asan tsan; do
    if [[ "$flavor" == asan ]]; then
        timeout_seconds=90
        options_name=ASAN_OPTIONS
        options_value="$asan_options"
    else
        timeout_seconds=120
        options_name=TSAN_OPTIONS
        options_value="$tsan_options"
    fi
    run_logged \
        "allocfail_$flavor" \
        "$timeout_seconds" \
        env "$options_name=$options_value" \
        "$output_dir/bin/test-ab-allocfail-$flavor" \
        "$field" 4 1 20 "$output_dir/allocfail_$flavor.wcs" \
        probe "$winner_index"
    log="$output_dir/allocfail_$flavor.log"
    assert_assist_lifecycle \
        "$flavor allocation-failure dynamic lending" \
        "$log" \
        1 \
        1 \
        1
    grep -Eq \
        '\[solver\] phase-profile detailed=1 failed=1 .*allocation_failures=1 .*helper_tasks=[1-9][0-9]* ' \
        "$log" ||
        die "$flavor assisted failure profile is missing"
    grep -Eq \
        '\[onefield-profile\] mode=serial-precommit-retry .*failed=0 ' \
        "$log" ||
        die "$flavor precommit serial recovery is missing"
    [[ "$(line_count '^AB_RESULT ' "$log")" -eq 1 ]] ||
        die "$flavor allocation recovery published duplicate final results"
    [[ "$(result_key "$log")" == "$(
        result_key "$output_dir/allocfail_baseline.log"
    )" ]] ||
        die "$flavor allocation recovery differs from W1"
done

if grep -Eq \
    'ERROR: AddressSanitizer|runtime error:|WARNING: ThreadSanitizer|SUMMARY: ThreadSanitizer|data race|releasing leased job-index-cache' \
    "$output_dir"/*.log; then
    die "a sanitizer diagnostic was found"
fi

lsan_status=not-run
if [[ "${RUN_LSAN:-0}" == 1 ]]; then
    run_logged \
        lsan \
        90 \
        env \
        ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:abort_on_error=1 \
        UBSAN_OPTIONS="$UBSAN_OPTIONS" \
        "$output_dir/bin/test-ab-stream-asan" \
        "$field" 4 1 20 "$output_dir/lsan.wcs" exhaustive \
        "$winner_index"
    assert_assist_lifecycle \
        "LSan dynamic lending" \
        "$output_dir/lsan.log" \
        1 \
        1
    lsan_status=passed
fi

selected_budget() {
    awk '
        {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^budget=/) {
                    split($i, value, "=")
                    result = value[2]
                }
            }
        }
        END {
            print result
        }
    ' "$1"
}

asan_wall_budget="$(selected_budget "$output_dir/asan_wall.selected")"
asan_cpu_budget="$(selected_budget "$output_dir/asan_cpu.selected")"
tsan_wall_budget="$(selected_budget "$output_dir/tsan_wall.selected")"
tsan_cpu_budget="$(selected_budget "$output_dir/tsan_cpu.selected")"
if [[ -z "$asan_wall_budget" ||
    -z "$asan_cpu_budget" ||
    -z "$tsan_wall_budget" ||
    -z "$tsan_cpu_budget" ]]; then
    die "one or more selected adaptive-limit budgets are missing"
fi

printf \
    'AB_SANITIZER_INTEGRATION_OK asan=passed ubsan=passed tsan=passed lsan=%s retained_state=index-mmap-hit,executor-quiescence full_owner=asan-tail,tsan-tail inverse_cache=asan-admit-hit,tsan-admit-hit adaptive_limits=asan-wall:%s,asan-cpu:%s,tsan-wall:%s,tsan-cpu:%s output=%s\n' \
    "$lsan_status" \
    "$asan_wall_budget" \
    "$asan_cpu_budget" \
    "$tsan_wall_budget" \
    "$tsan_cpu_budget" \
    "$output_dir" ||
    die "cannot write sanitizer integration summary"
