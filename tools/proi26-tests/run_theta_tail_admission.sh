#!/usr/bin/env bash
set -euo pipefail

fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
ROOT="$(dirname -- "$REPO")"
INSTALL="$ROOT/astrometry-install-solver-tasks"

BRANCH=assistant/verification-reorder-buffer
SOURCE_IMAGE="$REPO/demo/apod4.jpg"
BASE_CONFIG="$INSTALL/etc/astrometry.cfg"
INDEX_FILE="$ROOT/astrometry-data/index-4116.fits"
ENGINE="$INSTALL/bin/astrometry-engine"
SOLVE_FIELD="$INSTALL/bin/solve-field"
LISTHEAD="$INSTALL/bin/listhead"
BUILT_ENGINE="$REPO/solver/astrometry-engine"

cd "$REPO"
[ "$(git branch --show-current)" = "$BRANCH" ] || fail "switch to $BRANCH first"

for path in "$SOURCE_IMAGE" "$BASE_CONFIG" "$INDEX_FILE" \
            "$SOLVE_FIELD" "$LISTHEAD"; do
    [ -e "$path" ] || fail "missing required path: $path"
done

# Search production source only; do not match this runner's own marker strings.
if git grep -n -E '\[tune-fingerprint\]|\[verify-failure\]|tweak2_status_t' \
        -- include solver >/dev/null; then
    fail "temporary diagnostic code remains in production source"
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
CASE_NAME="apod4_theta_tail_fix_${STAMP}"
CASE_ROOT="$REPO/run_outputs/$CASE_NAME"
ARCHIVE="$REPO/run_outputs/${CASE_NAME}.tar.gz"
CASE_CONFIG="$CASE_ROOT/astrometry-4116.cfg"
mkdir -p "$CASE_ROOT"

{
    printf 'date=%s\n' "$(date --iso-8601=seconds)"
    printf 'branch=%s\n' "$(git branch --show-current)"
    printf 'head=%s\n' "$(git rev-parse HEAD)"
    printf 'verify_blob=%s\n' "$(git hash-object solver/verify.c)"
    printf 'helper_blob=%s\n' "$(git hash-object solver/verify_theta_tail.h)"
    printf 'test_blob=%s\n' "$(git hash-object solver/test-verify-theta-tail.c)"
} >"$CASE_ROOT/manifest.txt"

git status --short >"$CASE_ROOT/git-status-before.txt"
git diff 4b85abbdd235072117178fec387577243fa88a19..HEAD \
    >"$CASE_ROOT/branch.diff"
cp "$0" "$CASE_ROOT/admission-runner.sh"

JOBS="$(command -v nproc >/dev/null 2>&1 && nproc || printf 1)"
make -C "$REPO" -j"$JOBS" solver >"$CASE_ROOT/build.log" 2>&1 ||
    fail "solver build failed; see $CASE_ROOT/build.log"
make -C "$REPO/solver" test-verify-theta-tail \
    >"$CASE_ROOT/regression-build.log" 2>&1 || fail "regression build failed"
"$REPO/solver/test-verify-theta-tail" \
    >"$CASE_ROOT/regression-run.log" 2>&1 || fail "regression test failed"
grep -Fxq VERIFY_THETA_TAIL_OK "$CASE_ROOT/regression-run.log" ||
    fail "regression success marker missing"

[ -x "$BUILT_ENGINE" ] || fail "built engine missing"
install -m 0755 "$BUILT_ENGINE" "${ENGINE}.new.$$"
mv -f "${ENGINE}.new.$$" "$ENGINE"

for marker in '[verify-failure]' '[tune-fingerprint]' 'tweak2_status_t'; do
    LC_ALL=C grep -aFq "$marker" "$ENGINE" &&
        fail "installed engine contains temporary marker: $marker"
done

awk '
    /^[[:space:]]*(index|autoindex)([[:space:]]|$)/ { next }
    { print }
' "$BASE_CONFIG" >"$CASE_CONFIG"
printf '\nindex %s\n' "$INDEX_FILE" >>"$CASE_CONFIG"
[ "$(grep -Ec '^[[:space:]]*index[[:space:]]+' "$CASE_CONFIG" || true)" -eq 1 ] ||
    fail "generated config does not contain exactly one index"

grep -E '^[[:space:]]*(add_path|index|autoindex)([[:space:]]|$)' \
    "$CASE_CONFIG" >"$CASE_ROOT/config-directives.txt" || true

extract_wcs() {
    awk -F= '
        /^[[:space:]]*(CTYPE1|CTYPE2|CRVAL1|CRVAL2|CRPIX1|CRPIX2|CD1_1|CD1_2|CD2_1|CD2_2|IMAGEW|IMAGEH|A_ORDER|B_ORDER|AP_ORDER|BP_ORDER|A_[0-9]+_[0-9]+|B_[0-9]+_[0-9]+|AP_[0-9]+_[0-9]+|BP_[0-9]+_[0-9]+)[[:space:]]*=/ {
            key=$1; gsub(/[[:space:]]/, "", key)
            value=$2; sub(/[[:space:]]*\/.*/, "", value)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            print key "=" value
        }
    ' "$1" | sort >"$2"
}

for workers in 1 2 4; do
    out="$CASE_ROOT/w${workers}"
    log="$out/run.log"
    mkdir -p "$out"
    printf 'workers=%s\n' "$workers" >"$out/case.txt"

    rc=0
    ASTROMETRY_INDEX_SHARDING=pthread \
    ASTROMETRY_INDEX_SHARD_WORKERS="$workers" \
    LC_ALL=C /usr/bin/time -v "$SOLVE_FIELD" \
        --config "$CASE_CONFIG" --overwrite --wall-limit 120 \
        --scale-low 10 --scale-units degwidth --p-workers "$workers" \
        --no-plots -v --dir "$out" "$SOURCE_IMAGE" \
        >"$log" 2>&1 || rc=$?

    printf '%s\n' "$rc" >"$out/exit-status.txt"
    [ "$rc" -eq 0 ] || fail "W${workers} returned $rc"
    grep -Fq 'Field 1: solved with index index-4116.fits.' "$log" ||
        fail "W${workers} did not solve with index-4116"
    grep -Fq 'Tweak2: final WCS:' "$log" ||
        fail "W${workers} did not complete tuning"
    if grep -Eq '\[verify-failure\]|\[tune-fingerprint\]|verification produced no usable correspondence set|solver_tweak2: tune failed|preserving verified TAN solution' "$log"; then
        fail "W${workers} entered a forbidden fallback path"
    fi

    [ -f "$out/apod4.wcs" ] || fail "W${workers} WCS missing"
    [ -f "$out/apod4.axy" ] || fail "W${workers} AXY missing"
    "$LISTHEAD" "$out/apod4.wcs" >"$out/wcs-header.txt" 2>&1 ||
        fail "W${workers} WCS header read failed"
    extract_wcs "$out/wcs-header.txt" "$out/wcs-scientific-keys.txt"
    grep -Eq '^CTYPE1=.*SIP' "$out/wcs-scientific-keys.txt" ||
        fail "W${workers} CTYPE1 is not SIP"
    grep -Eq '^A_ORDER=' "$out/wcs-scientific-keys.txt" ||
        fail "W${workers} lacks A_ORDER"
    grep -Eq '^B_ORDER=' "$out/wcs-scientific-keys.txt" ||
        fail "W${workers} lacks B_ORDER"

    grep -E 'Tweak2: final WCS:|Field 1: solved with index|Field center:|Field size:|Field rotation angle:|Field parity:|User time \(seconds\):|System time \(seconds\):|Percent of CPU this job got:|Elapsed \(wall clock\) time' \
        "$log" >"$out/summary.txt" || true
    sha256sum "$out/apod4.axy" "$out/apod4.wcs" >"$out/output-sha256.txt"
    printf 'W%s_OK\n' "$workers"
done

cmp -s "$CASE_ROOT/w1/wcs-scientific-keys.txt" \
       "$CASE_ROOT/w2/wcs-scientific-keys.txt" || fail "W1/W2 WCS mismatch"
cmp -s "$CASE_ROOT/w1/wcs-scientific-keys.txt" \
       "$CASE_ROOT/w4/wcs-scientific-keys.txt" || fail "W1/W4 WCS mismatch"

cp "$CASE_ROOT/w1/wcs-scientific-keys.txt" \
   "$CASE_ROOT/common-wcs-scientific-keys.txt"
{
    printf 'result=PASS\n'
    printf 'head=%s\n' "$(git rev-parse HEAD)"
    printf 'regression=VERIFY_THETA_TAIL_OK\n'
    printf 'workers=1,2,4\n'
    printf 'sip_output=present\n'
    printf 'wcs_parity=exact\n'
    printf 'temporary_diagnostics=absent\n'
} >"$CASE_ROOT/result-summary.txt"

git restore -- __init__.py 2>/dev/null || true
git status --short >"$CASE_ROOT/git-status-final.txt"
tar -C "$REPO/run_outputs" -czf "$ARCHIVE" "$CASE_NAME"
sha256sum "$ARCHIVE" >"${ARCHIVE}.sha256"

printf '\nPOST_FIX_ARCHIVE_READY\n'
printf 'upload_required=%s\n' "$ARCHIVE"
printf 'archive_sha256=%s\n' "$(cat "${ARCHIVE}.sha256")"
