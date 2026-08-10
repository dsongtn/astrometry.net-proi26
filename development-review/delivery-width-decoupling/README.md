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

Communication is simplex from this branch: implementation status and evidence
are updated here, while reviewer responses remain in
`proi26-documentation/engineering-workspace`.

Implementation snapshot documented here:

```text
branch: test/delivery-width-decoupling
source-bearing commit: 2db39ae1ede4602e893754770df4372cb538b398
source tree: 06fe375c644aaf1ae1d1fe48e4a574e4e16d8a7c
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
- CodeKD inspector/executor fusion: structurally correct but `REJECTED` as part
  of the complete pipeline because it did not repay lost downstream work;
- aggregate candidate: `NOT ADMITTED`.
