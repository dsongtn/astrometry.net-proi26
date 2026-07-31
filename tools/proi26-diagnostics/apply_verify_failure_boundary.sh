#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    printf 'Usage: %s [--revert] [REPO]\n' "$0" >&2
}

MODE=apply
if [[ ${1:-} == --revert ]]; then
    MODE=revert
    shift
fi
if [[ $# -gt 1 ]]; then
    usage
    exit 2
fi

REPO=${1:-$PWD}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PATCH="$SCRIPT_DIR/verify_failure_boundary.patch"
EXPECTED_BASE=a07ba885ee41ec6d141b7457e053a438a4036ad1
EXPECTED_PATCHED=3f349817984fb7b19e70225eff1e37479dfc8bd2
EXPECTED_BRANCH=assistant/verification-reorder-buffer

if [[ ! -f $PATCH ]]; then
    printf 'ERROR: patch not found: %s\n' "$PATCH" >&2
    exit 1
fi
if ! git -C "$REPO" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    printf 'ERROR: not a Git worktree: %s\n' "$REPO" >&2
    exit 1
fi

BRANCH=$(git -C "$REPO" branch --show-current)
if [[ $BRANCH != "$EXPECTED_BRANCH" ]]; then
    printf 'ERROR: expected branch %s, found %s\n' "$EXPECTED_BRANCH" "$BRANCH" >&2
    exit 1
fi

if ! git -C "$REPO" diff --cached --quiet -- solver/verify.c; then
    printf 'ERROR: solver/verify.c has staged modifications\n' >&2
    exit 1
fi
if [[ $MODE == apply ]] &&
   ! git -C "$REPO" diff --quiet -- solver/verify.c; then
    printf 'ERROR: solver/verify.c has unstaged modifications\n' >&2
    exit 1
fi

CURRENT=$(git -C "$REPO" hash-object solver/verify.c)

if [[ $MODE == apply ]]; then
    if [[ $CURRENT != "$EXPECTED_BASE" ]]; then
        printf 'ERROR: unexpected solver/verify.c blob before apply\n' >&2
        printf 'expected=%s\nactual=%s\n' "$EXPECTED_BASE" "$CURRENT" >&2
        exit 1
    fi
    git -C "$REPO" apply --check "$PATCH"
    git -C "$REPO" apply "$PATCH"
    AFTER=$(git -C "$REPO" hash-object solver/verify.c)
    if [[ $AFTER != "$EXPECTED_PATCHED" ]]; then
        printf 'ERROR: patched blob mismatch\nexpected=%s\nactual=%s\n' \
            "$EXPECTED_PATCHED" "$AFTER" >&2
        exit 1
    fi
    rm -f \
        "$REPO/solver/verify.o" \
        "$REPO/solver/libastrometry.a" \
        "$REPO/solver/astrometry-engine"
    printf 'VERIFY_FAILURE_BOUNDARY_APPLIED\n'
    printf 'branch=%s\nverify_blob=%s\n' "$BRANCH" "$AFTER"
else
    if [[ $CURRENT != "$EXPECTED_PATCHED" ]]; then
        printf 'ERROR: unexpected solver/verify.c blob before revert\n' >&2
        printf 'expected=%s\nactual=%s\n' "$EXPECTED_PATCHED" "$CURRENT" >&2
        exit 1
    fi
    git -C "$REPO" apply -R --check "$PATCH"
    git -C "$REPO" apply -R "$PATCH"
    AFTER=$(git -C "$REPO" hash-object solver/verify.c)
    if [[ $AFTER != "$EXPECTED_BASE" ]]; then
        printf 'ERROR: reverted blob mismatch\nexpected=%s\nactual=%s\n' \
            "$EXPECTED_BASE" "$AFTER" >&2
        exit 1
    fi
    rm -f \
        "$REPO/solver/verify.o" \
        "$REPO/solver/libastrometry.a" \
        "$REPO/solver/astrometry-engine"
    printf 'VERIFY_FAILURE_BOUNDARY_REVERTED\n'
    printf 'branch=%s\nverify_blob=%s\n' "$BRANCH" "$AFTER"
fi
