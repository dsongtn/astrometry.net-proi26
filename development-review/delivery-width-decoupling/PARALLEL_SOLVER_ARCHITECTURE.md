# Parallel Solver Architecture

Branch-local review copy for `test/delivery-width-decoupling`. Do not merge
or publish this file with an official branch. The external
`proi26-documentation/engineering-workspace` is inbound review material;
branch implementation updates and responses to that advice are recorded only
in this directory.

Status: baseline `KEEP`; active branch candidate `EXPERIMENT / NOT ADMITTED`

Public baseline source tree:
`acc451dad31db38582ef35dc22ca4fb83db26cf1`

This document describes the public baseline at commit `e07c44c2`. It records
ownership and lifetime boundaries, not every internal function.

Active private candidate deviation:

```text
source-bearing commit: 2db39ae1ede4602e893754770df4372cb538b398
source tree: 06fe375c644aaf1ae1d1fe48e4a574e4e16d8a7c
source checkpoint: phase-05-owner-progress
mapped prime/requeue: REJECTED
CodeKD lookahead: one wave
owner progress and width policy: EXPERIMENT
```

The public architecture remains the release authority. Historical Phase 12
subsections are retained as audit records and explicitly labelled rejected or
superseded where the active branch no longer implements them.

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

The active branch derives payload completion width independently from outer
ownership. By default the requested payload width follows compute width and is
bounded by the provider's 24-job safety capacity. Internal build-only limits
can reduce producer or payload width for attribution; they are not public
command-line or environment controls.

For detached exact-demand passes, outer producer width is bounded by the live
payload-delivery width. Remaining compute workers are helpers for published
staged work. Resident, loaded-index, and other non-exact-demand passes retain
full outer width. Without detached completion, the legacy fixed-helper policy
remains the fallback. No tested width is currently admitted as the production
optimum.

No worker may block waiting for child work submitted to the same saturated
pool. Owner progress, helper publication, completion routing, and reduction
use explicit predicates and generation checks.

### 4.1 Bounded owner-local progress - EXPERIMENT

An owner with a published staged group receives one local selection
opportunity before it borrows global work. The bounded selection order is:

```text
COMPUTE -> IO completion -> SUBMIT -> PREPARE -> global fallback
```

This advances only the owner's already-published staged group. It does not let
a helper traverse mutable owner state, change index ownership, or retire the
owner's scientific results. Helpers retain global compute-first selection.

The active call path is:

```text
index_shard_staged_run_ordered()
  -> index_shard_owner_or_global_select_locked()
  -> index_shard_owner_progress_select_locked()
  -> index_shard_staged_select_locked(..., INDEX_SHARD_STAGED_SCOPE_OWNER)
```

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
default requested lanes: compute width
attribution lane limit: internal build control only
```

These are implementation safety bounds, not workload-specific performance
promises. Any change requires allocation, overflow, cancellation, and
backpressure tests.

### 5.1 Active packet window

The active branch publishes one descriptor wave in an ordered staged group.
Two waves remain the compile-time supported maximum only so the rejected
lookahead experiment can be reproduced against an exact source. The second
wave is not enabled in the candidate and must not be described as active.

This reduction does not change descriptor grain, sequence assignment, or
owner-only retirement.

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

Prepared pages are not pinned. Memory pressure can reclaim them. The provider
queue and I/O admission are bounded, and readiness is revalidated through
ticket and generation ownership rather than assumed to be permanent. The
existing provider-byte admission ends when the ticket becomes terminal; it
does not currently account for all READY pages through compute consumption.

The provider:

- sorts and coalesces physical file spans within strict gap and byte bounds;
- keeps exact mapped spans as the completion authority;
- reserves demand capacity ahead of speculative work;
- has explicit completion notification, cancellation, drain, and shutdown;
- never carries solver, index, or mapping pointers in the scheduler completion
  registry;
- records immutable numeric completion identities and validates owner epochs.

### 6.1 Active mapped delivery

The Phase 12 two-pass `prime -> requeue -> populate` lifecycle is `REJECTED`
and has been removed. The active mapped ticket remains within one service
execution:

```text
SUBMITTED
  -> exact plan and mapping validation
  -> bounded WILLNEED/readahead advice
  -> authoritative MADV_POPULATE_READ
  -> READY, FAILED, or CANCELLED
```

There is no prime requeue boundary and no detached READY interval between
advice and mapped population. Advice remains best-effort; authoritative
population may block one provider lane. The ticket retains its source lease,
duplicated descriptor, mapping generation, priority, sequence, and numeric
completion identity until one terminal transfer.

Provider width can exceed producer width in a controlled attribution build,
but configured width does not prove that lanes overlap usefully. A provider
lane may be waiting in mapped population while compute workers exhaust READY
work. Native mapped dereference remains the exact fallback on refusal, failure,
unsupported input, or cancellation.

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

## 11. Current evidence boundary

The source-exact Phase 05 APOD1 W4 `ABC` screen compared P1/D1, P1/D2, and
P2/D2 with one CodeKD wave and no mapped prime/requeue. All three candidates
reached the 120-second wall limit without a solution.

P1/D2 showed a severe single-screen regression: effective process CPU usage
fell to 0.25 cores and generation 2 contained an 88.90-second flattened-owner
phase.
P2/D2 made substantially more diagnostic progress than P1/D2, but neither a
joint timeout nor one ordering admits a width. The result does not prove that
data delivery, worker scaling, or owner progress is solved.

Current decisions are:

```text
ownership, generation checks, ordered retirement, native fallback: KEEP
mapped prime/requeue: REJECTED
second CodeKD wave: REJECTED for the active candidate
bounded owner progress: EXPERIMENT
producer/payload width policy: EXPERIMENT
aggregate branch candidate: NOT ADMITTED
```

The full source identities and results are recorded in
`ENGINEERING_PROGRESS_AND_EVALUATION.md`.

## 12. Reviewed next mechanisms - NOT IMPLEMENTED

The FITS/catalog-side review identified three separable mechanisms. They are
proposals, not descriptions of active production behavior.

### 12.1 Width factorization

Keep compute width, outer producer width, and delivery width independently
fixed during attribution. Select an absolute producer/delivery pair before
measuring C2/C4/C6. Do not add adaptive width control until fixed-width
evidence establishes a stable direction.

### 12.2 Delivery-stage ablation

The staged CodeKD pipeline includes leaf, Quad, candidate-Star, verification,
and sweep delivery. Remove one tail stage per experiment and use the native
owner path as the authoritative fallback. Verification-query construction
depends on candidate Quad/Star data, so not every stage can be enabled or
disabled independently without adding a new bridge.

The intended decision metric is useful scientific progress per wall second,
not ticket count, fault count, or CPU percentage alone.

### 12.3 Prepared-page lifetime budget

The provider already limits admitted jobs and bytes. A future memory-pressure
correction, if justified, must extend exact aligned-byte ownership from
admission through READY execution, cancellation, or retirement. Releasing the
lease at I/O completion would not bound READY-page pressure. This change must
include generation-safe cancellation and single-release accounting.

### 12.4 CodeKD inspector/executor fusion

Physical range order and logical solver order are separate. Physical page
ranges may be coalesced, but per-descriptor leaf spans must retain native child
and result order. Fusion is deferred until the reduced pipeline shows that the
second topology traversal remains material. W1 stays on the native serial path.

## 13. Safe extension rule

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
