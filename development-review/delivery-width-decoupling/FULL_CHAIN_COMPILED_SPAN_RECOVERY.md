# Full-chain Compiled-span Recovery

Status: `REJECTED`; local Rung 4 winner that failed repeatable Rung 5 field
admission. Active source returned to the Phase 00 boundary.

This is branch-local communication material for
`test/delivery-width-decoupling`. It is not release documentation.

## Decision being tested

The rejected CodeKD-only/native-downstream candidate removed useful Quad,
Star, verification, and sweep helper packages together with downstream page
delivery. It therefore could not isolate the economics of compiled CodeKD
spans. The corrected candidate tests one narrower question:

> Can a canonical compiled span program replace duplicate CodeKD topology
> traversal while the complete work-conserving downstream helper chain remains
> active?

The unchanged CodeKD-only candidate remains `REJECTED`.

## Production path

```text
DESCRIPTORS_READY
  -> owner or helper compiles canonical spans and exact mapped pages once
  -> PAGE_PLAN_COMPLETE
  -> delivery submission without a compute claim
  -> CODEKD_IO_SUBMITTED
  -> COMPUTE_READY
  -> any eligible worker executes retained spans
  -> DESCRIPTORS_READY for the next bounded slice, or RESULTS_READY
  -> candidate window
  -> Quad READY
  -> Star READY
  -> verification query / sweep / preparation / score READY
  -> owner-only canonical retirement
```

Physical page ranges may be sorted, deduplicated, and coalesced. Descriptor,
span, hit, candidate, and retirement order remain canonical.

## Ownership and lifetime

- One worker owns each mutable `solver_t` and its outer index task.
- Packet inputs and retained span programs are immutable while claimable.
- The mapped source lease and generation remain valid through the last mapped
  dereference.
- An I/O-pending packet holds no compute claim.
- A successful completion returns to the global READY set.
- A 16 MiB scheduler ledger bounds aligned mapped bytes from submit through
  execution or cancellation.
- Live-byte acquisition and release counts must balance at quiescence.
- Only the owner retires ordered results and mutates solver state.
- Late or cancelled completions release resources without publication.

## Exact fallback

The native owner path remains authoritative. It is used when:

- the KD tree or mapping is unsupported;
- a complete page or span plan cannot be produced;
- a byte, page, range, span, result, or candidate bound is exceeded;
- allocation or delivery fails;
- cancellation is observed;
- execution cannot safely consume the retained span program.

Fallback preserves the original query, tolerance, options, candidate order,
stop semantics, and owner-only publication boundary.

## Causal widths

The first probe uses fixed absolute internal widths:

```text
C = requested W4 compute width
P = 1 exact-demand outer producer
D = 1 mapped delivery lane
B = 16 MiB total READY-to-consumption live mapped bytes
N = 1 packet lookahead wave
```

P1/D1 is an attribution contract, not a general optimum. The same absolute
P, D, and B must be retained when W1, W2, and W4 are compared at Rung 5.

## Local result

Four counterbalanced pairs used the same APOD1 AXY, index 4107, W4, a five
second engine wall limit, and a 60 CPU-second ceiling.

| Pair | Control | Candidate | Gain |
| --- | ---: | ---: | ---: |
| 1 | 3.66 s | 2.98 s | 18.58% |
| 2 | 3.69 s | 2.96 s | 19.78% |
| 3 | 3.61 s | 3.01 s | 16.62% |
| 4 | 3.65 s | 2.97 s | 18.63% |

Median wall improved from 3.655 seconds to 2.975 seconds. CPU utilization
rose from a 120.0 percent median to 145.5 percent while the normalized WCS
digest remained identical. Filesystem input was essentially unchanged, so the
measured benefit is attributed to removal of repeated topology work with
preserved useful helper computation, not to less physical input.

## Verification

Completed against the exact postimage:

- libkd scalar/span parity: 30 tests;
- staged lifecycle: 20 deterministic cases;
- solver and worker-configuration suites;
- full downstream integration and helper-work accounting;
- allocation-failure and exact fallback paths;
- ASan and UBSan changed-path integration;
- Helgrind full-chain W4 path: zero reported errors;
- active clean build and installed-binary identity;
- APOD4 W1/W4 scientific smoke with identical normalized WCS.

TSan is `UNSUPPORTED` on the current VM because it aborts at startup before
the test path runs.

## Rung 5 field decision

The first quiet APOD1 candidate-then-control pair retired 483 candidate tasks
and 551 control tasks, a 12.34 percent candidate loss. The exact-AXY reversed
confirmation retired 577 candidate tasks and 553 control tasks, a 4.34 percent
candidate gain. Combined progress was 3.99 percent lower for the candidate,
and the paired directions disagreed materially.

Verbose attribution favored the candidate by 5.50 percent, but its log volume
changed execution economics and could not admit latency. It nevertheless
closed the structural questions:

- all 2,027,792 retained spans executed;
- there were no execution or hit-capacity replays;
- there was no non-cancellation plan refusal;
- foreign workers executed 61,003 compute packages versus 23,037 owner
  packages;
- Quad, Star, verification, and sweep stages remained active;
- submit-to-ready and owner-wait cost per event matched the control.

The candidate therefore did not fail because helper packages were removed,
and attribution did not identify a stable localized correction. APOD5 was
neutral with additional resource cost. APOD4 improved cold but regressed
slightly on its immediate repeat. The required repeatable full-field gain was
not established, so scaling was not admitted and the bounded-persistence rule
ended the mechanism.

The validated compiler, span-order, lifetime, fallback, and deterministic test
material remain archived for possible future reuse. They are not production
enabled. Active source was mechanically restored to the Phase 00 code boundary
without changing the two pre-existing user-owned worktree edits.

The restored boundary was then clean-built and installed. Every candidate
production path matches the Phase 00 checkpoint, with only one unrelated final
blank line omitted from `libkd/test_libkd.c`. Candidate identifiers are absent
from the source, `libkd.a`, the built engine, and the installed engine. The
restored integration path passed payload I/O, libkd, solver, worker-configuration,
staged-lifecycle, engine-pass, canonical-order, W1/W4 scientific, cancellation,
allocation-failure, and fallback checks. No post-rollback sanitizer rerun was
required because the rollback re-established the previously validated Phase 00
implementation without introducing production logic.

## Exact artifacts

```text
base HEAD:
  a75ff4757f486d24666b8d97356f882038ace072

clean review postimage tree:
  2cfb0c9aefac5949abe943ebc5b0c807cdb59d37

private checkpoint:
  refs/local-checkpoints/astrometry.net-solver-tasks/
  full-chain-compiled-span-recovery/phase-01-full-chain-local-winner

mechanism-stop checkpoint:
  refs/local-checkpoints/astrometry.net-solver-tasks/
  full-chain-compiled-span-recovery/phase-02-rung5-mechanism-stop

restored-boundary checkpoint:
  refs/local-checkpoints/astrometry.net-solver-tasks/
  full-chain-compiled-span-recovery/phase-03-phase00-restored

patch:
  run_outputs/full_chain_compiled_span_recovery_20260810T092905Z/
  full-chain-compiled-span-candidate.patch

patch SHA-256:
  c59b2a8085311855d39a91e4d151413f8db54074ba3accd57c77997174554f4c

compact evidence:
  run_outputs/full_chain_compiled_span_recovery_20260810T092905Z/
  EVIDENCE_SUMMARY.md

final quiet confirmation:
  run_outputs/full_chain_compiled_span_recovery_20260810T092905Z/
  rung5-field-confirm-20260810T104453Z-pid816072/
```
