# Inner Parallelism Branch Review Notes

Status: branch-local communication material only.

These files belong only to `test/inner-parallelism-and-scaling`. They are not
release documentation and must not be merged or published with an official
branch.

This directory is the outbound development record for this branch. The
external `proi26-documentation/engineering-workspace` is the inbound review
and advice workspace; implementation progress must not be written there.

Read:

1. `ENGINEERING_PROGRESS_AND_EVALUATION.md` for dates, phases, implementation
   steps, results, failures, and current candidate status.
2. `PARALLEL_SOLVER_ARCHITECTURE.md` for design, ownership, data delivery,
   scheduling, fallback, and Phase 12 deviations from the published baseline.

Communication is simplex from this branch: implementation status and evidence
are updated here, while reviewer responses remain in
`proi26-documentation/engineering-workspace`.

Implementation snapshot documented here:

```text
branch: test/inner-parallelism-and-scaling
commit: d0b5a0c66d6d3bdeebd51949bab56558091bf324
tree: 86044a499518ceb37267ef66d94160f718034ad3
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
- aggregate candidate: `NOT ADMITTED`.
