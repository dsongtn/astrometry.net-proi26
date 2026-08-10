# Final Reference Continuation

Status: `EXPERIMENT / FIELD ADMISSION REQUIRED`

This is branch-local communication material for
`test/delivery-width-decoupling`. It is not release documentation and must not
be merged into an official branch.

## Decision

The previous compiled-span candidate proved that helpers were active, but it
did not produce a repeatable full-field gain. This continuation therefore does
not restore that machinery. It makes two smaller changes against the restored
Phase 00 source boundary:

1. sparse mapped payload plans no longer bridge a positive page gap;
2. the first sufficiently large descriptor task is split into one 64-work-unit
   lead range and one canonical contiguous tail range.

The second change is intended to publish useful downstream CodeKD, Quad, Star,
and verification work earlier. It does not introduce a new queue, packet state,
worker class, I/O policy, reducer, or user-facing option.

The Phase 03 continuation adds one bounded scheduler experiment. A worker that
has just completed a valid outer index may execute a short burst of foreign
`COMPUTE_READY` packets before it claims another outer index. The experiment
does not revive a fixed helper reservation and does not make I/O-pending or
owner-only work eligible.

## FITS-side correction

Mapped exact-demand plans retain overlap and adjacency merging. They no longer
read an otherwise untouched page only to combine two physical ranges. File-
offset fallback behavior is unchanged. Native mmap remains the authoritative
reader, and allocation, mapping, advice, or delivery failure retains the exact
native fallback.

This change has a deterministic one-page-gap test. Its measured effect is an
upper bound on avoided mapped bytes, not evidence of an end-to-end wall-time
gain.

## Lead-packet mechanism

For the existing one-wave packet path only:

```text
original first task [A, B)
  -> lead [A, M)
  -> tail [M, B)
  -> all later tasks unchanged
```

`M - A` is the smallest combination count that expands to at least 64 native
hypothesis work units. A split occurs only when both lead and tail meet that
same minimum. At most one extra physical task is created per wave.

The existing scheduler may finish the lead through mapped completion and
downstream execution while the tail or later packets remain in CodeKD, I/O, or
compute stages. Owner-only retirement still requires the exact next task index
and `combination_first` sequence. If the extra workspace or wave slot cannot be
allocated, the original unsplit partition is used.

W1 remains on the native serial path. The mechanism contains no worker-count,
image, band, or index-family special case.

## Rejected alternative

An in-place rolling prefix inside one packet was rejected before production
editing. The present packet has one scheduler state and one result boundary.
Publishing a partial prefix would pause that packet's remaining CodeKD work;
it would not create true overlap without adding produced and retired
watermarks, dynamic result counts, partial trailing-counter ownership, and a
new cancellation protocol. That cost and risk were not justified.

## Bounded READY-first outer admission

The scheduler normally gives immediately claimable outer indexes priority.
That rule can leave useful downstream work READY while a recently freed lane
opens another cold index. Phase 03 permits only the following exception:

```text
finish one valid outer index
  -> grant a bounded READY-first budget
  -> claim foreign COMPUTE_READY work while inside the dynamic target
  -> reconsider outer work after at most two successful claims
```

The dynamic target is the smaller of the current `COMPUTE_READY` supply and
half the configured pool. A finishing lane outside that target immediately
returns to outer production. A miss, terminal state, selection error, or new
outer claim clears its budget. Preparation, submission, I/O polling, and
owner-only reduction cannot consume the budget.

The two-claim maximum is a liveness bound for this experiment, not a claimed
general optimum. It prevents one finished owner from draining an unbounded
foreign queue while canonical outer work remains. The final field decision
must consider READY dwell, owner wait, outer progress, page faults, input,
system time, and time to solution together.

## Exactness and safety evidence

The deterministic split test executes one nonzero native combination interval
both unsplit and as lead plus tail. It requires:

- identical hypothesis count and order;
- bit-identical active stars, codes, tolerance, and noise values;
- identical parity;
- identical cumulative `numtries`, `cxdx`, and `meanx` state at every retired
  hypothesis;
- identical final trailing counters and noise state.

Focused results against the current dirty experiment overlay:

| Gate | Result |
| --- | --- |
| Solver enumeration and split equivalence | pass |
| CodeKD lifecycle and cancellation | pass |
| Staged READY ordering and quiescence | pass |
| FITS payload exact-range suite | `OK (22 tests)` |
| Small W1/W4 scientific integration | pass |
| Forced packet-allocation fallback | pass |
| ASan and UBSan streaming path | pass |
| ASan and UBSan allocation-failure path | pass |
| Helgrind small W4 path | zero reported errors |

The Phase 03 scheduler extension also passed:

- deterministic W1, W2, W3, W4, W6, and W8 budget boundaries;
- two simultaneous W4 foreign claims with a half-pool cap;
- per-worker burst exhaustion and forced return to outer work;
- terminal-state budget clearing and complete counter closure;
- the focused solver, engine-pass, worker-configuration, payload, and staged
  suites;
- the small genuine permuted-StarKD W1/W4 integration;
- ASan and UBSan lifecycle and allocation-failure fallback runs;
- Helgrind on the changed W4 lifecycle with zero reported errors.

TSan is `UNSUPPORTED` on this VM because the runtime aborts before the test
path starts with an unexpected memory mapping. This is not recorded as a pass.

The clean installed Phase 03 engine identity is:

```text
astrometry-engine SHA-256:
  16505075f259242c0535a3d92fb839dcdca539c234ed8dbc4742b82ffa360d3c
build ID:
  6b8d935bb5308da7f3f62f1ca698763e79cb0665
solve-field SHA-256:
  7b65f75a1fea2c1d76df5c3f2593616be26114d4767760defe9c8a5c6215d648
configuration SHA-256:
  23e9c877b20bad9955ae8d0a799ee97a6cf096a83207072ce074f71731aa6c2a
ordered index count and hash:
  349
  7ee54821c993f3a8ac16e212256d428bcf3a7b2f3308571866d3f2d88a27337e
```

The small W4 fixture recorded one active lead split with 64 nominal work units.
The same W1 fixture recorded zero. These counters prove activation only. They
do not establish an APOD1 performance gain.

## Admission gate

Run one source-exact APOD1 W4 comparison against the restored Phase 00 control.
The candidate must preserve scientific output and ordered digests, show active
lead splits with non-inline staged work, and improve useful progress or solve
latency without materially increasing faults, filesystem input, system time,
or context switches.

For Phase 03, require `staged_ready_before_outer_claims` to be nonzero. If the
counter remains zero, the scheduler mechanism was dormant and cannot explain
the field result. Reject the experiment if READY dwell decreases only by
delaying the next outer claim, or if it moves blocking or mapped faults to
another lane without improving useful progress.

Reject and restore the private pre-change checkpoint if the field result is
neutral, noisy in opposite directions, or slower. Do not tune the 64-work-unit
lead using APOD1-specific data.

## Private rollback boundaries

```text
series: final-reference-continuation
phase-00-phase00-overlay
phase-01b-exact-mapped-ranges-clean
phase-02-lead-packet
phase-03-balanced-ready
```

These refs are local-only snapshots outside normal branch history.

## Phase 06 completion-aware CodeKD delivery

Status: `EXPERIMENT / FIELD ADMISSION REQUIRED`

The final continuation addresses two measured structural problems in the
mapped CodeKD path without changing native query execution:

1. KD topology page discovery now runs as bounded scheduler PREPARE work.
   The I/O lane receives an already sealed exact-mapping range plan and only
   refreshes the mapping and performs population. It no longer serializes KD
   topology traversal before population.
2. A queued ticket may reuse a recent completion only for the same live
   mapping generation and only when a service-lane `mincore()` check reports
   that page resident. A missing marker table, an old marker, a nonresident
   page, allocation pressure, or a failed residency check retains the page in
   the fresh population request.

Descriptor generation yields once before page planning. This keeps descriptor
construction and KD plan traversal separately claimable, allowing the
scheduler to choose already-READY compute between those coarse operations.
Cold first-use mappings do not allocate reuse metadata or call `mincore()`;
the metadata is created only after successful population.

`mincore()` is a point-in-time cache-residency observation. It does not pin a
page or prove that the mapping will remain fault-free. Native mmap access is
still authoritative, and eviction after READY produces an ordinary native
fault rather than a changed result. Service refusal, capacity pressure,
invalid mapping state, cancellation, or terminal I/O failure preserves exact
owner replay for every unfinished descriptor.

The following remain unchanged:

- native CodeKD range search and result order;
- descriptor and hypothesis sequence numbers;
- owner-only ordered retirement;
- Quad, Star, and verification mathematics;
- index order, winner publication, limits, and cancellation precedence;
- source, mapping, ticket, and index lifetime rules.

## Phase 06 bounded evidence

Focused validation completed on the dirty branch overlay:

| Gate | Result |
| --- | --- |
| FITS payload suite | `OK (22 tests)` |
| CodeKD packet suite | pass |
| Small W1/W4 scientific integration | pass |
| ASan and UBSan streaming matrix | pass |
| ASan and UBSan allocation-failure fallback | pass |
| Helgrind small W4 path | zero reported errors |
| TSan | `UNSUPPORTED`, runtime failed before the test path |

The clean installed candidate identity is:

```text
astrometry-engine SHA-256:
  1fb6d5007309a8cc7de2ef8b453d4e3925aff763d39c97a91b8d7a7c417eb801
build ID:
  52eb3ee99732feebad97707b0c944225375096f3
solve-field SHA-256:
  7b65f75a1fea2c1d76df5c3f2593616be26114d4767760defe9c8a5c6215d648
configuration SHA-256:
  23e9c877b20bad9955ae8d0a799ee97a6cf096a83207072ce074f71731aa6c2a
ordered index count and hash:
  349
  7ee54821c993f3a8ac16e212256d428bcf3a7b2f3308571866d3f2d88a27337e
```

Two quick, non-admission smoke observations were recorded after the clean
build. They were not cache-controlled repetitions:

| Case | Result | Wall | CPU | Major faults | Input |
| --- | --- | ---: | ---: | ---: | ---: |
| APOD4 W4, scale-low 10, first run | solved | 2.51 s | 77% | 505 | 342688 blocks |
| APOD4 W4, immediate repeat | solved | 1.39 s | 111% | 1 | 59256 blocks |
| APOD1 W4, 20 s screen | wall limit | 21.35 s | 237% | 40439 | 18234240 blocks |

The APOD1 screen reached the index-5200-33 area in field objects 11-20. This
is useful activation evidence only. The high input volume and the APOD4 result
show that cache and storage sensitivity remain material. The two APOD4 runs
also selected different finite-budget valid-looking solutions; neither was
admitted against an independent normalized WCS reference in this smoke step.

Raw smoke evidence:

```text
run_outputs/final-reference-apod4-smoke.Qc5YAL
run_outputs/final-reference-apod4-repeat.5gP6EU
run_outputs/final-reference-apod1-20s.ZO8JLT
```

The exact Phase 06 source is preserved at:

```text
refs/local-checkpoints/astrometry.net-solver-tasks/
  final-reference-continuation/phase-06-compute-plan-resident-reuse
checkpoint commit: 1494bf901d809f9491c6d5a3abbccf8812c19ec2
checkpoint tree:   200aaa3521020e4d55073eeb257f83c0151601dc
```

The remaining admission gate is one source-exact, guest-page-cache-cold APOD1
W4 run against the retained control, followed by normalized WCS validation for
any solution. A short screen cannot establish the under-120-second goal or a
reproducible wall-time speedup.

## Phase 09 single-pass exact-demand continuation

Status: `EXPERIMENT / NOT ADMITTED`

The local prior-art review in
`local-research/parallel-maturity-reduction/PRIOR_ART_AND_RESCUE_DESIGN.md`
identified one defensible rescue attempt: retain the existing bounded mapped
provider and complete downstream helper chain, but replace the second full KD
planning traversal with exact resumable native traversals. The implemented
path keeps a bounded ring of CodeKD descriptors. Each continuation advances in
native order until its next leaf, publishes only that leaf's mapped DATA/PERM
page cover, and resumes after mapped completion. Results are copied into
canonical descriptor slots and only the owner retires the contiguous prefix.

The implementation preserves:

- the W1 native authority;
- native split predicates, child order, result IDs, and squared distances;
- descriptor and candidate sequence order;
- the existing Quad, Star, and verification helper stages;
- generation-safe mapping leases, cancellation, and exact owner replay;
- native mmap as the authoritative reader after provider refusal or eviction.

### Causal results

The initial eight-continuation implementation amplified mapped submissions to
196,528 and took 7.91 seconds on the one-index APOD1 control. A packet-local,
page-keyed record of completed mapped pages reduced that to 83,449 submissions
and 5.12 seconds. The remaining bounded width study was:

| Active continuations | File-evicted wall | Mapped submissions | Voluntary switches | Peak RSS |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 3.38 s | 19,723 | 127,680 | 206,464 KiB |
| 128 | 3.32 s | 13,324 | 86,954 | 180,448 KiB |
| 256 | 3.32 s | 9,639 | 64,935 | 196,256 KiB |
| 512 | 3.24 s | 7,681 | 51,922 | 218,180 KiB |
| 1024 | 3.34 s | 6,725 | 47,891 | 245,492 KiB |

The table is single-run causal evidence, not benchmark admission. It proves
that wider rings reduce tickets and wakeups, but that the reduction stops
improving wall time and increases retained memory. The balanced source bound
therefore remains 128 continuations.

A global delivery-width-one build took 3.70 seconds. It throttled every
payload class and therefore does not test the reviewed CodeKD-only
consume-before-refill policy. A safe CodeKD-only credit must be scheduler-owned
from submit through leaf consumption; it cannot be attached to a ticket that
is destroyed at I/O completion. That larger lifecycle change was not rushed
into this experiment.

A typed bulk internal-node advance was also rejected. Four immediate-repeat
measurements produced a median of approximately 2.80 seconds versus
approximately 2.72 seconds for the simpler peek-and-step loop. The primitive
and its ABI surface were removed.

The completed-page identity history was enlarged from one 512-page delivery
batch to a bounded 8,192-entry metadata window. This stores mapping/page
identities only, not FITS payload. It reduced file-evicted submissions from
13,324 to 12,516 and reduced immediate-repeat voluntary switches by about 8%.
The file-evicted one-index observation improved from 3.32 to 3.24 seconds, but
the change alone does not make the full continuation candidate faster than the
retained Phase 08 control.

Extending provider reuse across the entire mapping lifetime was rejected. The
I/O lane correctly used `mincore()` to reduce mapped population from about
1.87 GiB to 0.33 GiB, but planning time rose and the file-evicted wall time
regressed to 3.30 seconds. This is direct evidence that fewer populated bytes
do not by themselves prove a faster solver.

### Decision and rollback

The exact-demand continuation is scientifically intact in focused tests, but
it did not cross its one-index causal gate and must not be described as a
performance improvement. The best measured, now-removed experimental shape
was:

```text
128 active exact native continuations
8-node compute quantum
8,192 packet-local completed-page identities
existing bounded mapped provider
existing downstream Quad/Star/verification packages
owner-only canonical retirement
```

The pre-wavefront source remains the accepted comparison boundary:

```text
series: final-reference-continuation
phase-08-pre-wavefront
checkpoint commit: eee1eda2683e7660a4ed900a92af2458dc24a9b2
checkpoint tree:   1e21efc99b610ba7fdb2e83d9b0fea467b5b104e
```

The rejected bulk-advance source is preserved only for diagnosis:

```text
phase-09-bulk-topology-advance
checkpoint commit: b0f2901516b51e02559cdbe858d645fd5b61017f
checkpoint tree:   ba6bb0f241bc15511b8180ce3ec18f610a94f94e
```

The final falsifying gate reconstructed Phase 08 from its exact checkpoint in
an isolated worktree. Its quiet one-index W4 immediate-repeat wall times were
2.50, 2.53, and 2.56 seconds, with a 2.53-second median. Its index-advisory-
evicted result was 3.01 seconds, 1,987 major faults, 29,010 voluntary context
switches, and 237,160 input blocks. The wavefront's best comparable quiet
results remained approximately 2.67 seconds resident and 3.10 seconds after
the same index-only eviction. Every control run solved with index 4107 and
produced the same normalized WCS digest.

This falsifies the continuation candidate at its mandatory first gate. The
production source and installed engine were restored exactly to Phase 08; the
wavefront source is not present in the final working tree or final link. The
installed restored engine has SHA-256
`1fb6d5007309a8cc7de2ef8b453d4e3925aff763d39c97a91b8d7a7c417eb801`
and build ID `52eb3ee99732feebad97707b0c944225375096f3`.

The exact control evidence is preserved locally under:

```text
run_outputs/phase08_exact_control_20260810T151616Z
```

The complete review-relevant provenance and result table are committed with
this branch under:

```text
development-review/delivery-width-decoupling/evidence/
  phase08-exact-control/
```

No APOD1 field claim is made for the rejected wavefront. Reconsideration would
require a materially different mechanism that first beats the Phase 08
one-index control while preserving the same normalized WCS; further ring,
lookahead, or page-history tuning is not justified by these results.

The mapping-lifetime `mincore()` reuse policy was also isolated directly on
the restored Phase 08 source. Its immediate-repeat results were 2.63 to 2.93
seconds and its advisory-evicted result was 3.04 seconds, versus the exact
Phase 08 ranges above. It was therefore reverted as neutral-to-negative; the
original 24-job recent-completion horizon remains installed.
