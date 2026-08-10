# Consumption-Bounded Exact-Demand Result

Status: `REJECTED / RESTORED`

```text
date: 2026-08-10
branch: test/delivery-width-decoupling
head: 9cbba4738642c6a3ade7bbdf8675a7dfed836e8b
source_bearing_commit: 2db39ae1ede4602e893754770df4372cb538b398
phase_00_tree: 5e486e9e0c6e0fc90d4942faf94ed3b82a7f70fa
candidate_tree: bf65bf23d0d2a45f51b442e94b941854a9828286
restored_tree: 5e486e9e0c6e0fc90d4942faf94ed3b82a7f70fa
decision: REJECT
```

## Iterations

```text
production_source_iterations: 2
  1: complete D1 compiled exact-demand pipeline
  2: telemetry, fallback, lifecycle, and verification correction

causal_build_variants: 2
  1: D2 delivery width
  2: compiled CodeKD plus restored downstream delivery

rollback_iterations: 1
total_cycles_after_phase_00: 5
```

The verified-install checkpoint did not introduce another source delta. It
recorded the same candidate tree after the clean install and focused gates.

## Candidate

The candidate implemented one CodeKD exact-demand pipeline:

```text
compute_width: requested W1/W2/W4
outer_producer_width: 1
mapped_delivery_width: 1
live_aligned_byte_limit: 16777216
lookahead_waves: 1
post_codekd_delivery: disabled
native_downstream: Quad, Star, verification, sweep
```

Each admitted descriptor compiled its canonical leaf spans and exact mapped
page plan in one topology traversal. Execution reused the retained spans,
held a generation-safe byte lease through the final mapped dereference, and
fell back to exact owner replay on refusal or failure.

## Structural and scientific result

The mechanism itself passed its focused gates:

- one topology traversal per planned descriptor;
- zero topology replay during retained-span execution;
- exact scalar/span result IDs, order, and squared distances;
- zero post-CodeKD tickets in the primary build;
- balanced live-byte acquires and releases below 16 MiB;
- exact owner replay on refusal and allocation failure;
- no observed WCS or ordered-digest mismatch in the focused fixtures;
- ASan, UBSan, Valgrind, and Helgrind checks passed;
- TSan could not start on this VM and remains `UNSUPPORTED`.

No critical scientific, ownership, lifetime, or memory-safety defect caused
the rejection.

## Decisive performance result

The causal probe used one APOD1 AXY, only `index-4107.fits`, W4, a five-second
engine wall limit, a 60 CPU-second limit, and file-cache eviction advice before
each execution. Four counterbalanced repetitions were used. This was a
bounded source-exact falsifier, not canonical 349-index admission.

| Variant | Runs | Median wall | Median CPU | Median major faults | Median voluntary switches |
|---|---:|---:|---:|---:|---:|
| Phase 00 | 4 | 2.930 s | 164.5% | 3,720 | 34,762.5 |
| Compiled CodeKD, D1, native downstream | 4 | 3.905 s | 94.0% | 16,769.5 | 20,740.5 |
| Compiled CodeKD, D2, native downstream | 4 | 4.035 s | 97.5% | 17,025 | 21,908 |
| Compiled CodeKD plus old downstream delivery | 4 | 2.980 s | 149.5% | 910 | 24,994 |
| Restored engine | 1 | 2.870 s | 155.0% | 3,699 | 31,798 |

Every D1 pair regressed by 32.5% to 35.1%; the median regression was 33.3%.
D2 did not recover the loss. Restoring the old downstream delivery recovered
approximately baseline performance but did not produce a useful improvement.

## Root cause and decision

The falsified assumption was:

```text
compiled CodeKD delivery plus native downstream access supplies enough useful
parallel work to repay planning, population, scheduling, and retirement cost
```

It did not. The removed Quad, Star, and verification packages supplied much
of the useful helper work. With those packages absent, owners performed more
downstream work locally while helper availability fell. Effective CPU use
dropped, sparse mapped faults increased, and result-to-retire delay grew.
Eliminating CodeKD topology replay was correct but too small a saving to repay
the complete pipeline.

The complete candidate was therefore rejected. No subcomponent was retained
because the hybrid only returned to baseline while preserving additional
complexity. The live source and installed engine were restored to Phase 00.

## Evidence and rollback

```text
candidate_checkpoint:
  refs/local-checkpoints/astrometry.net-solver-tasks/
  consumption-bounded-exact-demand/phase-02-focused-gates

restored_checkpoint:
  refs/local-checkpoints/astrometry.net-solver-tasks/
  consumption-bounded-exact-demand/phase-03-rejected-restored

rejected_patch_sha256:
  9c35c280c9feecb556118a2872811da563a9ad013f33f3197c6e616d175eea0c

local_evidence_directory:
  run_outputs/phase12_consumption_bounded_exact_demand_20260810T020512Z

summary_sha256:
  5a66fd174ff201a5731ae60c6eaf28faabc7b184c4eb2e98e99d7959f9c67097

restored_engine_sha256:
  3bc1567d1287795a1acb4146a136089b2e74dcd5624ce769680acb618a9b8a34
```

Private checkpoint refs and raw run outputs are local evidence and are not
part of normal branch publication. Reconsider this design only with new
evidence that preserves useful downstream parallelism and predicts a material
critical-path gain.
