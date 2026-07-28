# Narrow-index solve-field benchmark

This harness measures the production regime where scale and other constraints
leave one or a few candidate indexes and per-index hypothesis work can grow as
more or larger field objects are admitted. It produces the complete
I1/I2/I4/I8 x configured field-object range ID x W1/W2/W4 x OS-page-cache
state matrix without running two solves concurrently.

The timing suite is an evidence gate, not an adaptive-policy implementation.
It preserves raw logs and rejects timing conclusions when the exact constructed
command, immutable candidate prefix, field-object interval, run completeness,
or exact W1 result differs from the declared experiment. A companion phase
suite validates the verbose log's observed candidate/pass shape and provides
attribution.

## What one suite does

Let R be the number of configured field-object ranges (R must be at least two).
For each of at least seven repetitions, the runner executes all 12R logical
index/range/worker conditions, each with two OS-page-cache states: exactly 24R
timed cells per repetition. A seeded shuffle fixes the initial condition order;
each repetition rotates that order by a coprime stride. Every condition is run
as an immediate page-cache-cold/reuse pair:

1. synchronize dirty pages;
2. invoke the configured privileged page-cache drop;
3. start a fresh `solve-field` process labelled `cold`;
4. immediately start the identical fresh process labelled `warm`, without
   another page-cache drop;
5. parse and validate both runs.

Seven repetitions produce 84R pairs and 168R timed processes. With the two
ranges in the example, that is 168 pairs and 336 processes. The `warm` label
means immediate OS-page-cache reuse with a new process and a new solver pool;
it does not mean persistent process or persistent index residency. The `cold`
label means a Linux kernel page-cache drop was completed before the run. It
does not establish that filesystem-server, hypervisor, controller, or
device-internal caches are cold, so the recorded storage topology remains part
of the experiment definition.

The authoritative measurement launches `taskset` directly and reaps that exact
child with Linux `wait4`. Wall time comes from nanosecond `perf_counter_ns`;
user/system CPU, RSS, faults, context switches and filesystem I/O come from the
returned child-tree rusage without GNU `time`'s 0.01-second display
quantization. The raw record is preserved in `resource-usage.json`. The wall
boundary includes Python launch/reap and the `taskset` exec, while CPU includes
`taskset`, solve-field and descendants reaped by solve-field. At suite startup,
31 identical `taskset /usr/bin/true` null launches are recorded in provenance
as a fixed-boundary estimate. This estimate is never subtracted: it bounds when
millisecond results approach launcher resolution without pretending to be an
exact per-run correction. The runner pins the complete process tree to one
configured CPU set with `taskset` and forces OpenMP, OpenBLAS, MKL, BLIS,
ATLAS, GotoBLAS, vecLib and NumExpr
thread counts to one. The runner uses the backward-compatible
`ASTROMETRY_INDEX_SHARD_WORKERS` environment input to vary worker counts within
the matrix; production commands may use `--p-workers` instead.

## Prepare an experiment

Use one output root per input image and per measurement mode. Copy
`indexes.example.txt` and replace every line with an exact absolute path. The
order must be the real production admission order; do not move a known winning
index forward. The runner requires exactly eight distinct files, rejects globs,
hashes every index and declared pipeline file, and forms I1/I2/I4/I8 only from
immutable prefixes.

Copy `config.example.json` and set:

- the exact feature-branch `solve-field` and matching `wcsinfo` binaries;
- every executable or script that the selected input path can launch in
  `pipeline_files` (at minimum the `astrometry-engine` adjacent to
  `solve-field`; use a pilot log to include image conversion, PNM,
  source-extractor, line-removal and uniformization helpers actually used);
- one production input image (or the exact pre-extracted input stage being
  studied);
- the absolute index-manifest and output paths;
- the production scale, sky and other constraints in `common_args`;
- an explicitly ordered `field_object_ranges` array with at least two named
  inclusive intervals;
- a fixed CPU set containing at least four suitable physical cores;
- a cache-control method that succeeds non-interactively.

Each range entry has a filesystem-safe `id` and an `objects` value such as
`"1-10"`. IDs are labels, not solver classes: the first declared range can
contain any field objects, and no ID carries an intrinsic "early", "late",
size, difficulty, or scheduler meaning. Declaration order is preserved in the
plan and reports. The harness passes each exact `objects` interval through
solve-field's native `--depth` option; it does not change astrometry.net's
field-object ordering or interpretation. Range specifications and IDs must be
unique, but the harness does not infer a physical distinction between ranges.
Set `rotation_stride` coprime to the resulting 12R logical conditions; the
example value 5 is valid for its two ranges.

The example uses `sudo -n tee /proc/sys/vm/drop_caches`. A narrowly scoped,
audited cache-drop helper is preferable on a shared machine. Alternatively,
set `cache_control.method` to `proc` and run the harness as root. There is no
uncontrolled page-cache-cold fallback: inability to drop the page cache aborts
the suite. This operation affects the whole host, so use an otherwise idle
dedicated test machine and coordinate with its operator.

The harness owns verbosity, output names, the native `--depth` field-object
interval, index and config options. It rejects those options in `common_args`.
Every `common_args` entry must be a single long-option token; attach values with
`=`, for example `--scale-low=1.0`. Split option/value entries and positional
tokens are refused so an undeclared second input cannot enter a quiet timing
run. Keep exactly one scale interval and one input field. Run the companion
`phase` suite before accepting a timing result: if scale filtering removes a
declared prefix member, or if the command creates more than one solver pass,
its exact verbose-log assertion aborts instead of silently attributing a
different workload. Quiet timing mode cannot observe that internal log shape
and records this limitation explicitly.

Before a real collection, stabilize the machine: disable opportunistic jobs,
record the intended CPU governor/turbo policy, avoid SMT siblings in the CPU
set when possible, and keep binaries, input, indexes and results on their
production-representative filesystems. The provenance bundle captures file
hashes, stat identities, CPU/memory information, selected environment values,
binary versions, `ldd`, and git state. It deliberately does not dump the full
ambient environment because that can contain credentials.

## Run and resume

From this directory:

```sh
python3 run_matrix.py --config /absolute/path/to/config.json
```

The output root must be empty. A completed pair has an atomic `pair.json`
marker. After a clean interruption between pairs, resume with:

```sh
python3 run_matrix.py --config /absolute/path/to/config.json --resume
```

Resume recomputes input and binary hashes and accepts only an identical
experiment fingerprint. Every completed pair and run record is bound to that
fingerprint and to a cross-mode workload fingerprint, then revalidated against
its deterministic plan position before it
can be skipped. Command, solver log, raw wait4 resource output, `wcsinfo`, solved marker,
WCS, and page-cache-control evidence are hashed; resume and summarization
reject missing or altered evidence. They also recompute the canonical scientific WCS
signature, `wcsinfo` digest and raw-WCS digest, and require solved metadata, marker, WCS, and the
unsolved sentinel to agree. An interrupted or failed partial pair is retained for
diagnosis and is never overwritten; start a new output root after inspecting
it. The runner never deletes evidence.

Each state directory contains `command.json`, `solve.log`,
`resource-usage.json`, `run.json`, `wcsinfo.txt`, solver outputs and, for the
page-cache-cold state, `cache-control.json`. The plan and full provenance are
at the suite root.

## Timing, phase, and detailed diagnostics

Use `measurement_mode: "timing"` for crossover decisions. It invokes no
`--verbose` option, so the authoritative latency numbers contain neither
verbose logging nor phase-record overhead. It validates the exact persisted
and executed argv, ordered immutable index prefix, field-object interval,
single final input, exit status, solved/WCS agreement, correctness signatures
and evidence hashes.
It does not pretend that candidate/pass log shape or phase attribution is
observable: run records say `profile_validation: not_observable_quiet`, and
helper-use output says `NOT_OBSERVABLE_QUIET`.
Any structured profile record in a quiet run is treated as contamination and
invalidates the suite.

For low-overhead phase attribution and observed candidate/pass assertions,
copy the timing config to a new output root and set `measurement_mode: "phase"`.
This invokes one `--verbose` and captures the source, engine, onefield, pool,
context, reducer and solver aggregates in `phase_summary.tsv`. Phase-suite
latency comparisons are `DIAGNOSTIC_ONLY`; they do not set the crossover.

CodeKD, exact reduction and verification timing instruments hot calls and are
therefore enabled only by `-vv`. To collect them, copy the same config, change
the output root and set `measurement_mode: "detailed"`. That suite repeats the
same controlled matrix with two `--verbose` options and captures
`[solver] phase-profile`, including `codekd_work_wall_sum`,
`resolve_work_wall_sum`, `reduction_ex_verify_hit_work_wall_sum`, and
`verify_hit_work_wall_sum`. It also captures `[solver] verification-context`
so bounded verification-tile submission, observed overlap, ordered commits,
discarded speculative candidates, and serial fallbacks remain auditable. Its
comparisons are explicitly marked
`DIAGNOSTIC_ONLY`; never use phase or detailed suites to set the crossover.
Diagnostic reports leave `evidence_selected_workers` blank; their independently
reported `raw_fastest_workers` observation is not a production selection.

Each suite records a `workload_fingerprint` over the normalized configuration,
plan, binaries, input, indexes, declared pipeline files and harness. Only
`measurement_mode` and `output_root` are excluded. After both roots complete,
make the timing report validate the phase evidence explicitly:

```sh
python3 summarize.py /absolute/path/to/timing-results \
  --companion /absolute/path/to/phase-results
```

Add a detailed root with another `--companion` when collected. The comparator
fully revalidates companion completeness, evidence hashes, profiles and exact
correctness, requires matching workload fingerprints and cross-mode result
signatures, and refuses duplicate modes. A timing result is eligible for
acceptance only when `companion_validation.status` is `PASS` and
`timing_acceptance_ready` is `true`; the automatic summary produced before the
phase suite exists deliberately reports `NOT_CHECKED` and is not an acceptance
artifact.

Fields ending in `_work_wall_sum`, `_elapsed_sum`, `_cpu_sum`, or
`_queue_wait_sum` are sums of work intervals. Parallel sums can exceed pass
elapsed time and must not be added together or interpreted as critical-path
shares. End-to-end elapsed time comes from the process-tree timer and the
explicit `total`/pass records.

When a profile prefix occurs more than once in one process (for example, one
`onefield-field-profile` and one detailed solver profile per attempted index),
the generic phase table sums each numeric field within that run before taking
the cross-run median/MAD. Its `record_count` exposes that multiplicity. These
are work aggregates, not process critical-path elapsed time.

The authoritative reducer aggregate from `[index-shard] reducer-pass` and
solver aggregate from `[index-shard] solver-pass` are retained.
`[index-shard] context-pass` separates worker-local onefield/solver preparation
and cleanup from index-task work; its sums are parallel work totals and its
maxima approximate the largest per-worker contribution to pass latency.
`[index-shard] phase-profile` reports `task_work_wall_sum` plus the
`reset_work_wall_sum`, `acquire_work_wall_sum`, `solve_work_wall_sum`,
`analyze_work_wall_sum`, `release_work_wall_sum`, and `other_work_wall_sum`
components. Separate `*_percent` fields are descriptive; the phase table keeps
both them and the numeric work sums under their explicit names.

Dynamic fixed-pool lending has its own evidence. Every owner publishes and
unpublishes a temporary `assist-lane`; a pool worker joins only after the outer
index queue is exhausted, borrows one bounded phase task, leaves, and can then
migrate to another live lane. `[index-shard] assist-pass` gives authoritative
loan, notification and final-waiter counts. The harness balances every
publish/unpublish and join/leave pair and requires the aggregate loan count to
equal observed joins. There is no predefined owner/helper partition and no
`N < W` eligibility rule.

Detailed mode also records every `[solver-ab-phase]` with the exact mode:
`native`, `empty`, `flattened-owner`, or `assisted`. These modes describe
observed work at that frontier; they are not scheduler classes and need not be
the same in repeated runs. The phase record's `user`, `system`, and
`major_faults` values come from `RUSAGE_SELF` around a lane-local interval while
other owners and helpers may run. They are therefore explicitly labelled
`resource=process-overlap` and must not be summed or attributed to that lane.
`[index-shard] pass-resource` is the authoritative process-wide resource delta
for the submitted pass. End-to-end `resource-usage.json`, measured by `wait4`,
remains authoritative for the complete process tree.

The parser still recognizes historical `aux-pass` records so old evidence
archives remain readable, but current reports derive lending evidence from
`assist-pass`, exact phase modes, `helper_tasks`, and `helper_combinations`.

## Assertions and correctness

For every run, the runner verifies that `command.json` exactly matches the argv
given to the timed process, with the declared verbosity count, ordered immutable
index prefix, exact field-object interval supplied to native `--depth`, and one
final configured input argument.

For phase and detailed runs it additionally verifies from the solver log:

- exactly one `solver run parameters` block;
- exactly one input field in that block;
- the exact ordered index prefix and candidate count;
- the exact internal `startobj`/`endobj` for the declared pass;
- no pthread submission for W1;
- exactly one full-width pthread submission for W2/W4 whenever at least one
  index is eligible, with the configured pool size, candidate count and pass
  bounds;
- `inner_scheduler=dynamic-pool-lending` for every pthread submission,
  independent of index count, field-object ordinal, declared range order, and
  bounded-geometry outcome;
- no bounded-prefix geometry outcome for a zero-based range, and exactly one
  shared or refused outcome when a nonzero lower bound is reused across
  pthread owners or multiple W1 indexes;
- balanced dynamic-lane publish/unpublish and join/leave lifecycles;
- exactly one `assist-pass` record with `waiters=0` and `loans` equal to the
  observed join count;
- in detailed mode, only the four exact AB phase modes listed above, with every
  per-phase `RUSAGE_SELF` scope declared as `process-overlap` or `unavailable`.

Repeated per-index field profiles must match `serial_executed` for every
native-serial pass. For parallel runs their count must not exceed the outer-task
`executed` count, which itself is bounded by the declared candidate count;
cooperatively cancelled outer tasks may execute without entering
`solve_fields`, so equality is not required in that mode.

Timing runs deliberately record observed log shape and winning index as
`not_observable_quiet`; these are not inferred from missing quiet-log lines.

Solved status is derived from the solved marker and must agree with WCS output.
The authoritative winning index order is extracted from the reducer for
parallel phase/detailed runs and from the verbose serial attempt stream for W1,
bounds-checked against the immutable prefix, and reported with each diagnostic
condition.
The canonical correctness signature is the SHA-256 of a Python-standard-library
parse of `solve.wcs`. It includes the TAN core, image geometry, all standard
WCS matrix/projection keys, and every SIP A/B/AP/BP order and coefficient card.
Card order, formatting, DATE, HISTORY and COMMENT do not affect scientific
equality. Complete `wcsinfo` output and raw WCS bytes retain independent hashes
as diagnostic and chain-of-custody evidence, and resume/summarization recompute
all three. W1 must be deterministic across all repetitions and both cache states;
W2/W4 must match that one W1 signature exactly. A different but potentially
valid out-of-order solution is reported as a correctness failure for review,
not silently accepted into performance results.

## Statistics and crossover output

`summarize.py` is Python-standard-library only and runs automatically after a
complete collection. It can also be rerun:

```sh
python3 summarize.py /absolute/path/to/results/root
```

It reports unscaled median and MAD for timing/resource metrics and all numeric
phase fields. W2 and W4 are paired with W1 by repetition. A deterministic
paired bootstrap resamples the per-repetition speedup ratios and reports a 95%
interval (20,000 resamples by default).

The crossover flags are conservative:

- `PARALLEL`: at least one parallel width has a speedup interval wholly above
  1.0 and exact W1 correctness passed;
- `SERIAL`: both W2 and W4 intervals are wholly below 1.0;
- `UNCERTAIN`: neither conclusion is supported;
- `INVALID`: evidence completeness or exact correctness failed;
- `DIAGNOSTIC_ONLY`: phase logging or hot-call instrumentation was enabled.

Raw fastest width is reported separately and does not override uncertainty.
CPU-work ratios and major faults accompany wall-time speedups so a latency win
that spends substantially more CPU or changes the I/O regime remains visible.
No multiplicity correction is applied; treat the intervals as engineering
evidence for this declared workload, not population-wide statistical claims.

The summary directory contains:

- `report.md`  - compact human-readable decisions;
- `condition_summary.tsv`  - median/MAD time and resource metrics;
- `phase_summary.tsv`  - median/MAD structured instrumentation fields;
- `helper_usage.tsv`  - dynamic loans, notification counts, exact phase modes,
  helper task/combinations, and observed-assistance classifications;
- `parallel_vs_w1.tsv`  - paired speedup intervals and CPU-work ratios;
- `crossover.tsv`  - serial/parallel decision per index/range/OS-page-cache cell;
- `correctness.tsv`  - exact W1 signature validation;
- `summary.json`  - the same information in a machine-readable form.

Run the self-tests before deploying a modified harness:

```sh
python3 -m unittest -v test_harness.py
```

## Related document

- [Validation and field-test checklist](VALIDATION.md)
