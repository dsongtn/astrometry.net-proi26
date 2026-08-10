# Delivery Width Decoupling Branch Review Notes

Status: branch-local communication material only.

These files belong only to `test/delivery-width-decoupling`. They are not
release documentation and must not be merged or published with an official
branch.

This directory is the outbound development record for this branch. The
external `proi26-documentation/engineering-workspace` is the inbound review
and advice workspace; implementation progress must not be written there.

Read:

1. `ENGINEERING_PROGRESS_AND_EVALUATION.md` for dates, phases, implementation
   steps, results, failures, and current candidate status.
2. `PARALLEL_SOLVER_ARCHITECTURE.md` for design, ownership, data delivery,
   scheduling, fallback, and delivery-width deviations from the published
   baseline.
3. `CONSUMPTION_BOUNDED_EXACT_DEMAND_RESULT.md` for the final implementation,
   causal measurement, rejection, and exact restoration result from the
   2026-08-10 bounded pipeline attempt.
4. `FULL_CHAIN_COMPILED_SPAN_RECOVERY.md` for the corrected integration that
   retained the complete downstream helper chain, won locally, failed
   repeatable field admission, and was restored to Phase 00.
5. `FINAL_REFERENCE_CONTINUATION.md` for the exact-demand continuation attempt,
   bounded causal measurements, rejected alternatives, and final exact Phase
   08 restoration.
6. `evidence/phase08-exact-control/` for the source-pinned one-index control
   provenance, complete timing table, scientific digest, and capsule hashes.

Communication is simplex from this branch: implementation status and evidence
are updated here, while reviewer responses remain in
`proi26-documentation/engineering-workspace`.

Implementation snapshot documented here:

```text
branch: test/delivery-width-decoupling
documentation HEAD before this result: a947b992e9cdbe757ab74c6b690e903272952bc1
restored source checkpoint tree: 1e21efc99b610ba7fdb2e83d9b0fea467b5b104e
rejected wavefront checkpoint tree: ba6bb0f241bc15511b8180ce3ec18f610a94f94e
base: e07c44c2cd7995c684bfa20ab3dc7b861038e1b8
review copy prepared: 2026-08-09
last branch update: 2026-08-10
```

Current implementation decisions:

- mapped prime/requeue: `REJECTED` and removed;
- two-wave CodeKD lookahead: `REDUCED` to one wave;
- bounded owner-local progress: `EXPERIMENT`, focused-test verified;
- producer and payload width: `EXPERIMENT`, one `ABC` screen complete and no
  width admitted;
- delivery-stage ablation: `REJECTED`; removing downstream packages caused a
  repeatable 33.3% median regression in the source-exact bounded probe;
- READY-to-consumption page-byte lease: implemented and verified only in the
  rejected candidate, then removed by exact restoration;
- CodeKD inspector/executor fusion: structurally correct and locally faster,
  but not admitted because its full-field effect was not repeatable;
- unchanged CodeKD-only/native-downstream candidate: `REJECTED`;
- full-chain compiled-span replacement: `REJECTED`, despite its local Rung 4
  win, because the two quiet APOD1 field pairs disagreed and combined progress
  was lower;
- positive-gap mapped-range bridging: `REMOVED`; overlap and adjacency merging
  remain, while sparse mapped demand no longer reads an otherwise untouched
  page solely to join two ranges;
- in-place rolling descriptor prefix: `REJECTED` before implementation because
  the existing all-or-nothing packet state cannot overlap one packet's CodeKD
  tail with its own downstream retirement;
- bounded first-packet lead split: `EXPERIMENT`; exact traversal equivalence,
  lifecycle, sanitizer, allocation-failure, and small-fixture checks pass, but
  no APOD1 field-performance claim exists yet;
- single-pass CodeKD exact-demand wavefront: `REJECTED`; its best bounded
  one-index measurements remained slower than the exact Phase 08 control;
- mapping-lifetime `mincore()` reuse: `REJECTED`; it reduced population bytes
  but added enough planning and system cost to remain neutral-to-negative;
- aggregate field candidate: `NOT ADMITTED`; active source restored to the
  exact Phase 08 boundary.
