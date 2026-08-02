# Validation and field-test checklist

## Build and smoke checks

Repeat these from the repository root:

```sh
make -C catalogs libcatalogs.a
make -C solver -j4 astrometry-engine
python3 -m unittest -v \
  solver/benchmarks/narrow-index/test_harness.py
python3 -m py_compile \
  solver/benchmarks/narrow-index/run_matrix.py \
  solver/benchmarks/narrow-index/summarize.py \
  solver/benchmarks/narrow-index/fits_wcs_signature.py \
  solver/benchmarks/narrow-index/test_harness.py
git diff --check
```

Use a clean rebuild in the target installation. The specialized integration
and sanitizer objects are deliberately listed in `NODEP_OBJS`, so an
incremental build after an internal-header change is not admissible evidence.
After the clean build, run both bundled gates and provide the mandatory
genuinely permuted StarKD fixture to each:

```sh
solver/check-solver-parallel-integration.sh \
  demo/apod4.xyls \
  solver/index-9918.fits \
  demo/index-4119.fits \
  /absolute/path/to/permuted-sweep.fits

SOLVER_TEST_OUTPUT_DIR=/absolute/path/to/empty-sanitizer-results \
solver/check-solver-parallel-sanitizers.sh \
  demo/apod4.xyls \
  solver/index-9918.fits \
  demo/index-4119.fits \
  /absolute/path/to/permuted-sweep.fits
```

## Real matrix

Create one timing config per production image and run:

```sh
python3 solver/benchmarks/narrow-index/run_matrix.py \
  --config /absolute/path/to/mobile-case-timing.json
```

Repeat in a different empty output root with
`"measurement_mode": "phase"`. The quiet timing suite decides crossover; the
one-verbose phase suite validates the observed candidate/pass shape, attributes
the coarse phases and records outer-worker activity. If hot-call timing or
queue-wait evidence is required, run a third empty root with
`"measurement_mode": "detailed"` (`-vv`). Phase and detailed latency
comparisons are diagnostic-only. Each seven-repetition suite contains 48
unique cells and 336 timed processes when two field-object ranges are
configured. In general, R configured ranges produce 24R cells per repetition
and 168R timed processes across seven repetitions.

The copied configs must differ only in `measurement_mode` and `output_root`.
After the phase root completes, run the mandatory cross-mode gate:

```sh
python3 solver/benchmarks/narrow-index/summarize.py \
  /absolute/path/to/timing-results \
  --companion /absolute/path/to/phase-results
```

If a detailed root exists, append another `--companion`. Do not accept a timing
crossover unless the resulting `companion_validation.status` is `PASS` and
`timing_acceptance_ready` is `true`. The comparator revalidates the companion's
matrix, evidence, profiles and correctness and requires the stable workload
fingerprint and result signatures to match across modes.

Every `common_args` value-bearing option must use one token with `=`, such as
`--scale-units=arcsecperpix`; positional and split option/value tokens are
rejected.

`--no-verify` in the example means "ignore a pre-existing WCS header in the
input FITS image." It prevents a header-verification shortcut; it does not
disable `verify_hit()` for newly generated candidates.

Use at least two representative images, including one solved within an earlier
configured field-object range and one requiring a later and/or larger
configured range in the same declared sequence. "Earlier" and "later" are
relative to that experiment's range order, not intrinsic kinds of astrometry.net
work. Preserve result roots unmodified when returning evidence.

The field-test pthread path requires exactly one requested XYLS field and
`nsolves <= 1`. Multi-field input is intentionally routed through the exact
legacy serial path until a field-aware all-fields-complete reducer is
implemented and validated. `nsolves > 1` also stays serial until qualifying
matches can be counted exactly across private index result slots.
Nonzero `maxquads` or `maxmatches` also selects serial execution until those
process-wide limits have an authoritative cross-worker counter.
Automatic RDLS tag-along discovery (`rdls_tagalong_all`) selects serial
execution because its discovered column list is not yet worker-private.
Every safety-gated or precommit-fallback run uses the ordinary one-index-at-a-
time serial solver. In particular, a resolved worker count of one never starts
the pthread pool or enters legacy grouped-index execution.

## Production-admission failure matrix

Use a test-only fault injector or debugger breakpoint; do not add runtime
production environment hooks merely to inject failures.

The bundled local gates cover a deterministic active-phase allocation failure
with safe precommit retry, active cancellation, in-flight wall and CPU limits,
ordinary and later-index winners, W1/W4 result and WCS equality, topology
boundaries, and lazy StarKD initialization under ASan/UBSan and TSan. They do
not claim to cover every engine/output failure below. Items without retained
evidence remain production-admission work; passing the local gates alone does
not activate this candidate as the production default.

1. Fail hypothesis allocation/search before any reducer commit. Expect one
   safe serial retry and no partial output from the failed pass.
2. Fail a different shard after a valid result has committed. Expect terminal
   nonzero status, no serial fallback, no final solution output and no solved
   marker.
3. Place a valid-looking solution in a result slot and then mark that task
   failed. The reducer must reject it.
4. Race one valid solved shard against one hard failure. Hard failure must not
   become apparent success, regardless of completion order.
5. Trigger the cancel file from a worker. Master cancellation and process
   status/output must match the documented cancellation contract.
6. Fail a later field read in a multi-field XYLS. Expect immediate stop,
   aggregated failure status, nonzero engine exit and no final output/marker.
7. Compare serial and parallel accepted WCS/correctness signatures for every
   solved fixture. Investigate any difference even if both WCSes appear valid.
   Confirm a changed SIP coefficient is rejected after canonical signature
   recomputation, while DATE/HISTORY/COMMENT changes leave only the raw WCS
   digest different.
8. Run I1/W4, I2/W4 and I>=W under ThreadSanitizer with a representative
   shifted configured field-object range. Confirm every pthread pass reports
   `dynamic-pool-lending`; outer index claims retain priority; workers that
   exhaust the outer queue borrow one bounded task at a time from any published
   lane; and all lanes, loans and waiters quiesce before index release. Exercise
   solve, cancellation, limits and an injected helper-arrival failure. Dynamic
   lending must not depend on bounded prefix geometry, an absolute field-object
   number, or `N < W`.
9. Make the solved-marker destination unwritable after output setup. Expect a
   nonzero engine status; a completed solution output may already exist because
   marker publication deliberately occurs after output succeeds.

Solved markers are intentionally staged and written only after final solution
output succeeds. This prevents a failed concurrent run from poisoning a retry.
The trade-off is that a multi-field process killed before successful final
output no longer retains markers for fields solved earlier in that process.
Failure to publish a staged marker is now a hard run failure. Publication is
not atomic with the separate solution files, so callers must treat the process
status""”not mere output-file existence""”as authoritative.

## Evidence interpretation

- `_work_wall_sum`, `_elapsed_sum`, `_cpu_sum` and `_queue_wait_sum` fields are
  summed work intervals, not critical-path elapsed time.
- Per-frontier `[solver-ab-phase]` `RUSAGE_SELF` deltas are labelled
  `resource=process-overlap`: other owners/helpers can run inside that
  interval, so those deltas are diagnostic and neither lane-local nor
  additive. Use `[index-shard] pass-resource` for the process-wide submitted
  pass and `resource-usage.json` for the complete reaped process tree.
- `native`, `empty`, `flattened-owner`, and `assisted` are exact observed phase
  modes. A range may contain any mixture; do not turn those outcomes into
  permanent policies for a field-object ordinal or configured range ID.
- Quiet timing runs cannot observe internal pass shape; pair them
  with a phase suite before accepting the experiment definition.
- Phase (`-v`) and detailed (`-vv`) runs are diagnostic and cannot set the
  crossover; detailed mode additionally instruments hot calls.
- A page-cache-cold-labelled run is valid only when Linux page-cache control
  succeeded. That label says nothing about filesystem-server, hypervisor,
  controller, or device caches. Zero major faults means that run did not test
  storage-bound behavior.
- CPU work is the high-resolution `wait4` child-tree rusage, not GNU-time text.
  Compare millisecond wall results with the recorded 31-run taskset/null-launch
  distribution. No launcher value is subtracted; a result near that distribution
  needs a longer-lived server benchmark before adopting a crossover policy.
- Never compare image-input and pre-extracted-AXY results as the same latency
  target; provenance must identify the input stage.
