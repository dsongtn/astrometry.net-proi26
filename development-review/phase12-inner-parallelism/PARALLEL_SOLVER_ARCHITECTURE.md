# Parallel Solver Architecture

Branch-local review copy for `test/inner-parallelism-and-scaling`. Do not merge
or publish this file with an official branch. Review responses belong in the
external `proi26-documentation/engineering-workspace`.

Status: `KEEP`

Applies to source tree:
`acc451dad31db38582ef35dc22ca4fb83db26cf1`

This document describes the public baseline at commit `e07c44c2`. It records
ownership and lifetime boundaries, not every internal function.

Active private candidate deviation:

```text
branch commit: adc027322558f05b6063f9e62a9a33ae903d3ea3
branch tree: c29c6039c1e25db6b4a546f1956fcad2516069cd
original checkpoint: phase-12-two-stage-mapped-prime
status: PROVISIONAL
```

The public architecture remains the release authority. Subsections explicitly
labelled Phase 12 describe only the private candidate.

## 1. End-to-end structure

```text
solve-field
  -> astrometry-engine job and pass control
  -> onefield pass
  -> W1: native serial solver
     W2+: fixed index-shard pthread pool
       -> outer index owner
       -> worker-private solver_run()
       -> bounded CodeKD and verification packages
       -> bounded payload I/O service
       -> owner-only ordered retirement
       -> reducer-only master publication
```

Parallel execution is internal to one image solve. It does not depend on
running unrelated images concurrently.

## 2. Ownership model

| State | Owner | Sharing rule |
| --- | --- | --- |
| Engine job and pass policy | Engine/main thread | Creates and tears down pass-scoped services |
| Mutable `onefield_t` and `solver_t` | One worker | Never concurrently mutated by another worker |
| Index task | One outer owner | Executed at most once; reaches one terminal state |
| Immutable field view and index metadata | Pass/job | Borrowed only through bounded lifetimes |
| CodeKD packet inputs | Packet owner | Helpers receive immutable descriptors and disjoint outputs |
| Payload ticket | Provider until completion, then scheduler/owner | Exactly one terminal transfer or drain |
| Candidate result | Worker, then reducer | Frozen before publication and transferred once |
| Final master result | Reducer | At most one authoritative commit |

Cooperative stop predicates govern cancellation. The design does not use
`pthread_cancel` and does not create a nested worker pool.

## 3. Scientific order

Outer and inner order have different contracts:

- Outer quick commit may publish the first completed scientifically valid WCS,
  even when its index is not first in the configured manifest.
- Within one index-owned pipeline, descriptors, hits, candidates, prepared
  verification results, and reducer-visible outputs retire in their required
  legacy sequence.
- Helpers may complete out of order but cannot retire out of order.
- A different valid parallel winner under finite limits is reported and
  independently validated; it is not silently called identical.

W1 remains the reference for enumeration, mathematics, and scientific output.

## 4. Scheduler widths

The public control is `--p-workers auto|count`.

- `auto` uses CPUs visible through process affinity.
- `1` selects the native serial path.
- Invalid, zero, negative, or over-affinity counts are rejected.

The fixed pool has two internal roles:

1. Outer producers own indexes and run the native solver flow.
2. Surplus workers claim already-published staged or verification work.

For frozen public production exact-demand passes, outer producer width is
bounded by live payload-delivery width. The public I/O lane ceiling is four.
Workers beyond the producer width are not permanently idle; they are eligible
for inner packages. This width separation reduces simultaneous cold mappings,
but its wider-width performance is provisional.

No worker may block waiting for child work submitted to the same saturated
pool. Owner progress, helper publication, completion routing, and reduction
use explicit predicates and generation checks.

## 5. CodeKD packet pipeline

The CodeKD path converts a bounded contiguous descriptor range into a packet.
The detailed state machine includes separate Quad, Star, verification-query,
verification-sweep, scoring, execution, result, and retirement states. Its
logical shape is:

```text
descriptors ready
  -> complete bounded page plan
  -> I/O submitted, or exact owner fallback
  -> compute ready
  -> executing
  -> results ready
  -> canonical owner retirement
```

A packet is not considered compute-ready merely because some ranges were
prepared. Incomplete, refused, unsupported, cancelled, or failed preparation
routes to the exact owner/native behavior or an explicit terminal error.

Current packet bounds include:

```text
CodeKD result bytes per packet: 2 MiB
CodeKD mapped-delivery bytes per packet: 2 MiB
provider in-flight bytes: 64 MiB
provider in-flight jobs: 24
public provider lane ceiling: 4
```

These are implementation safety bounds, not workload-specific performance
promises. Any change requires allocation, overflow, cancellation, and
backpressure tests.

### 5.1 Private Phase 12 packet window

Phase 12 permits one following descriptor wave in the same ordered staged
group. The maximum logical lookahead is two waves. This changes physical
overlap opportunity, not descriptor grain, sequence assignment, or owner-only
retirement. It is not part of the frozen public baseline.

## 6. Mapped-page delivery

Native FITS mmap remains authoritative. The payload provider is an overlap and
readiness mechanism, not a second general-purpose page cache.

The current policy is deliberately mixed:

- Serial callers retain `MADV_NORMAL`.
- Parallel compact KD topology retains `MADV_NORMAL`.
- Parallel sparse payload mappings use `MADV_RANDOM`.
- Selected exact mapped ranges may be populated with
  `MADV_POPULATE_READ` when the platform supports it.
- Advice or provider failure returns to the native path.

Mapped completion is keyed to the mapping that will later be dereferenced,
not merely to an underlying file offset. File-page identity can avoid
duplicate storage work, but it cannot replace mapping identity when readiness
means that the exact compute mapping has been populated.

Prepared pages are not pinned. Memory pressure can reclaim them. Therefore the
READY queue and in-flight byte budget remain bounded, and readiness is
revalidated through ticket and generation ownership rather than assumed to be
permanent.

The provider:

- sorts and coalesces physical file spans within strict gap and byte bounds;
- keeps exact mapped spans as the completion authority;
- reserves demand capacity ahead of speculative work;
- has explicit completion notification, cancellation, drain, and shutdown;
- never carries solver, index, or mapping pointers in the scheduler completion
  registry;
- records immutable numeric completion identities and validates owner epochs.

### 6.1 Private Phase 12 two-stage mapped delivery

The private candidate separates advisory priming from authoritative mapped
population:

```text
SUBMITTED
  -> exact plan
  -> bounded WILLNEED/readahead prime
  -> requeue when another ticket is waiting
  -> authoritative MADV_POPULATE_READ
  -> READY, FAILED, or CANCELLED
```

The first pass does not publish READY. If no other ticket is queued, it falls
through to population instead of requeueing. One ticket retains the duplicated
descriptor, source lease, mapping generation, priority, sequence, and
completion identity across both passes.

The candidate also allows the payload service to start at the requested live
compute width, bounded by the existing provider job and byte ceilings. For an
exact-demand pass, producer width follows this live delivery width. W4 therefore
uses four producers and four delivery lanes; W8 may use eight. Wider behavior
is an experiment and can increase cold mappings, page pressure, and storage
contention.

Priming is advisory. It does not prove that storage completed before the
authoritative population pass, and prepared pages remain reclaimable. Phase 12
therefore improves ordering opportunity but is not a true asynchronous storage
completion interface.

## 7. Whole-cohort residency

The source contains a bounded whole-file residency service, but job-local use
is explicitly disabled. It copied complete cohorts inside the solve wall,
could block without solver-limit polling, and treated reclaimable memfd pages
as if they were permanently resident.

Production therefore logs exact-demand mode and returns without starting the
service. Do not enable this code as a performance shortcut. Reopening it would
require a persistent pre-job lifecycle, recoverable residency state, source
identity proof, strict admission against physical/cgroup/address-space limits,
and a separate field experiment.

## 8. Module map

The recognizable upstream facade files remain:

```text
solver/engine.c
solver/onefield.c
solver/index_shard.c
solver/solver.c
solver/verify.c
util/fitsbin.c
util/starkd.c
```

Project-owned responsibilities are split into 32 translation units:

| Subsystem | Modules | Responsibility |
| --- | --- | --- |
| Engine | `engine_job.c`, `engine_pass.c`, `engine_policy.c`, `engine_residency.c` | Job, pass, worker policy, and residency boundary |
| Onefield | `onefield_job_cache.c`, `onefield_index_shard.c` | Field/index view lifetime and shard bridge |
| Index shard | `index_shard_control.c`, `index_shard_helper.c`, `index_shard_inverse.c`, `index_shard_pass.c`, `index_shard_pool.c`, `index_shard_profile.c`, `index_shard_reducer.c`, `index_shard_scheduler.c`, `index_shard_staged.c`, `index_shard_worker.c` | Pool, task, completion, helper, reduction, and terminal state |
| Solver core | `solver_profile.c`, `solver_field_geometry.c`, `solver_hypothesis.c` | Profiling, geometry, and hypothesis generation |
| CodeKD | `solver_codekd_plan.c`, `solver_codekd_delivery.c`, `solver_codekd_verification.c`, `solver_codekd_staged.c`, `solver_codekd_retire.c` | Packet planning, delivery, verification preparation, completion, and ordered retirement |
| Verification | `verify_score.c`, `verify_projection.c`, `verify_prepared.c` | Scoring, projection, and prepared-hit ownership |
| FITSBIN | `fitsbin_mmap.c`, `fitsbin_payload_source.c`, `fitsbin_payload_plan.c`, `fitsbin_payload_service.c` | Mapping policy, source identity, range planning, and bounded service |
| StarKD | `starkd_payload.c` | Star payload preparation |

Private headers live beside their owners and are not installed. Public ABI
continues through `include/astrometry/`.

## 9. Required lifecycle invariants

```text
executions[task] <= 1
terminal_transitions[task] == 1
commits <= 1
post_terminal_claims == 0
queued_tasks == 0 at quiescence
running_tasks == 0 at quiescence
live_contexts == 0 at quiescence
reference_acquires == reference_releases
```

Also require:

- no mapped pointer survives index release;
- no pass-N completion enters pass N+1;
- no condition wait occurs outside a predicate loop;
- no late completion publishes after terminal stop;
- partial initialization has one complete cleanup path;
- the first non-fatal terminal event is the pass linearization point, while a
  later global-integrity failure can still invalidate an elected but
  uncommitted winner;
- fallback never replays or drops scientific work.

## 10. Rejected architecture directions

Do not reintroduce these without new contrary evidence and a bounded gate:

- shared mutable `solver_t`;
- nested pools or owner tasks waiting on children in the same pool;
- tiny per-permutation or per-candidate tasks;
- same-index affinity as the governing scheduling policy;
- a second general-purpose user-space page cache;
- whole-corpus preload on an oversized corpus;
- broad `MADV_WILLNEED` or global blanket advice changes;
- millions of synchronous tiny `pread()` operations;
- foreign traversal of an owner's index state;
- persistent index-family routing sidecars;
- hard-coded APOD, index-family, scale, or search-order shortcuts;
- candidate dropping, changed tolerances, or reduced search coverage;
- extra threads whose only benefit is a higher reported CPU percentage.

## 11. Safe extension rule

Before adding an inner module, identify:

1. creator and owner;
2. immutable input boundary;
3. task grain and measured ceiling;
4. output and canonical retirement order;
5. cancellation and generation behavior;
6. mapping, source, and index lease lifetime;
7. bounded memory and queue use;
8. native fallback;
9. final linked object and test coverage;
10. one falsifying field experiment.

Move code only when the destination owns a coherent responsibility. Do not
split state machines merely to reduce line count.
