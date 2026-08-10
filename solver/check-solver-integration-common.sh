#!/usr/bin/env bash
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE

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
    # Retained for call compatibility with the retired lending checker.
    local legacy_minimum_foreign="$4"
    local allow_task_failures="${5:-0}"
    local submits
    local lifecycle
    local helper_records
    local staged_records
    local scheduler_records
    local foreign_work
    local lifecycle_errors
    local task_failures

    submits="$(line_count \
        'pthread-pool submit .*inner_scheduler=ordered-codekd-packets .*page_delivery=detached-bounded-mapped-completion' \
        "$file")"
    lifecycle="$(awk '
        /^SOLVER_TEST_RESULT / &&
                ($0 ~ /cancelled=1/ ||
                 $0 ~ /wall_limit=1/ ||
                 $0 ~ /cpu_limit=1/) {
            terminal = 1
        }
        /\[index-shard\] helper-pass generation=/ {
            helper_records++
            groups = completed = foreign = failures = 0
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^groups=/) {
                    split($i, value, "=")
                    groups = value[2] + 0
                } else if ($i ~ /^completed=/) {
                    split($i, value, "=")
                    completed = value[2] + 0
                } else if ($i ~ /^foreign_tasks=/) {
                    split($i, value, "=")
                    foreign = value[2] + 0
                } else if ($i ~ /^task_failures=/) {
                    split($i, value, "=")
                    failures = value[2] + 0
                }
            }
            foreign_work += foreign
            task_failures += failures
            if (groups != completed) {
                balance_errors++
            }
        }
        /\[index-shard\] staged-pass generation=/ {
            staged_records++
            groups = completed = foreign = failures = 0
            submitted = io_completed = 0
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^groups=/) {
                    split($i, value, "=")
                    groups = value[2] + 0
                } else if ($i ~ /^completed=/) {
                    split($i, value, "=")
                    completed = value[2] + 0
                } else if ($i ~ /^foreign_claims=/) {
                    split($i, value, "=")
                    foreign = value[2] + 0
                } else if ($i ~ /^task_failures=/) {
                    split($i, value, "=")
                    failures = value[2] + 0
                } else if ($i ~ /^io_submitted=/) {
                    split($i, value, "=")
                    submitted = value[2] + 0
                } else if ($i ~ /^io_completed=/) {
                    split($i, value, "=")
                    io_completed = value[2] + 0
                }
            }
            foreign_work += foreign
            task_failures += failures
            if (groups != completed || submitted != io_completed) {
                balance_errors++
            }
        }
        /\[index-shard\] scheduler-observability generation=/ {
            scheduler_records++
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^completion_active=/) {
                    split($i, value, "=")
                    if (value[2] + 0 != 0) {
                        lifecycle_errors++
                    }
                }
            }
        }
        END {
            if (!terminal) {
                lifecycle_errors += balance_errors
            }
            print helper_records + 0, staged_records + 0,
                scheduler_records + 0, foreign_work + 0,
                lifecycle_errors + 0, task_failures + 0
        }
    ' "$file")"
    read -r helper_records staged_records scheduler_records \
        foreign_work lifecycle_errors task_failures <<<"$lifecycle"

    [[ "$submits" -eq "$generations" ]] ||
        die "$label has $submits current-architecture submits, expected $generations"
    [[ "$helper_records" -eq "$generations" ]] ||
        die "$label has $helper_records helper-pass records, expected $generations"
    [[ "$staged_records" -eq "$generations" ]] ||
        die "$label has $staged_records staged-pass records, expected $generations"
    [[ "$scheduler_records" -eq "$generations" ]] ||
        die "$label has $scheduler_records scheduler records, expected $generations"
    : "$legacy_minimum_foreign" "$foreign_work"
    [[ "$lifecycle_errors" -eq 0 ]] ||
        die "$label contains $lifecycle_errors unbalanced lifecycle records"
    if [[ "$allow_task_failures" -eq 0 &&
        "$task_failures" -ne 0 ]]; then
        die "$label contains $task_failures failed helper tasks"
    fi
    if grep -qE \
        'assist-pass|assist-lane state=|phase-group|static-phase-assist' \
        "$file"; then
        die "$label contains obsolete assist telemetry"
    fi
}

assert_full_producer_assistance() {
    local label="$1"
    local file="$2"
    local values
    local full_width
    local ownership_complete
    local foreign_work

    values="$(awk '
        /\[index-shard\] pthread-pool submit / {
            compute = producer = helper = candidates = -1
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^compute_width=/) {
                    split($i, value, "=")
                    compute = value[2] + 0
                } else if ($i ~ /^producer_width=/) {
                    split($i, value, "=")
                    producer = value[2] + 0
                } else if ($i ~ /^helper_width=/) {
                    split($i, value, "=")
                    helper = value[2] + 0
                } else if ($i ~ /^candidates=/) {
                    split($i, value, "=")
                    candidates = value[2] + 0
                }
            }
            if (compute > 0 && producer == compute && helper == 0 &&
                candidates >= compute) {
                full_width++
            }
            submitted_candidates = candidates
        }
        /\[index-shard\] ownership-pass generation=/ {
            claims = unclaimed = quiescent = -1
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^canonical_claims=/) {
                    split($i, value, "=")
                    claims = value[2] + 0
                } else if ($i ~ /^unclaimed=/) {
                    split($i, value, "=")
                    unclaimed = value[2] + 0
                } else if ($i ~ /^quiescent=/) {
                    split($i, value, "=")
                    quiescent = value[2] + 0
                }
            }
            if (claims == submitted_candidates &&
                unclaimed == 0 && quiescent == 1) {
                ownership_complete++
            }
        }
        /\[index-shard\] helper-pass generation=/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^foreign_tasks=/) {
                    split($i, value, "=")
                    foreign_work += value[2]
                }
            }
        }
        /\[index-shard\] staged-pass generation=/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^foreign_claims=/) {
                    split($i, value, "=")
                    foreign_work += value[2]
                }
            }
        }
        END {
            print full_width + 0, ownership_complete + 0,
                foreign_work + 0
        }
    ' "$file")"
    read -r full_width ownership_complete foreign_work <<<"$values"

    [[ "$full_width" -eq 1 ]] ||
        die "$label did not retain full outer producer width"
    [[ "$ownership_complete" -eq 1 ]] ||
        die "$label did not complete every outer ownership claim"
    [[ "$foreign_work" -gt 0 ]] ||
        die "$label did not exercise dynamic foreign assistance"
}
