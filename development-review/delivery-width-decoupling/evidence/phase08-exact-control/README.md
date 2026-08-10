# Phase 08 Exact One-Index Control Evidence

Status: `CONTROL / SOURCE-EXACT / NOT CANONICAL COLD`

This compact evidence capsule belongs only to
`test/delivery-width-decoupling`. It records the exact control used to decide
the final single-pass continuation experiment. It is not release admission
evidence and must not be generalized to the complete 349-index APOD1 field.

## Purpose

The control answers one bounded causal question: does the proposed CodeKD
continuation wavefront beat the restored Phase 08 implementation when both
solve the same APOD1 AXY against the same single index with four workers?

The source was reconstructed from private checkpoint commit
`eee1eda2683e7660a4ed900a92af2458dc24a9b2`, built in an isolated worktree,
and executed on CPUs 0-3. The source and binary identities are recorded in
`provenance.txt`; every measured result is in `summary.tsv`.

Equivalent engine arguments:

```text
astrometry-engine -p -c none \
  -i astrometry-data/index-4107.fits \
  -d RUN_DIR apod1.axy
```

No verbose tracing, scale hint, downsampling, or index shortcut was used. The
input AXY and index hashes are frozen in `provenance.txt`.

## Results

| Case | Wall | User | System | CPU | Major faults | Voluntary switches | Input blocks |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Warmup | 2.58 s | 3.48 s | 0.57 s | 157% | 0 | 26,300 | 0 |
| Immediate 1 | 2.50 s | 3.39 s | 0.53 s | 156% | 0 | 26,531 | 0 |
| Immediate 2 | 2.53 s | 3.43 s | 0.56 s | 157% | 0 | 26,647 | 0 |
| Immediate 3 | 2.56 s | 3.45 s | 0.54 s | 156% | 0 | 26,297 | 0 |
| Index advisory-evicted | 3.01 s | 3.50 s | 1.27 s | 158% | 1,987 | 29,010 | 237,160 |

The immediate-repeat median is 2.53 seconds. The advisory-evicted run used
successful `POSIX_FADV_DONTNEED` only on `index-4107.fits`. It did not clear
the guest kernel page cache, hypervisor cache, AXY pages, executable pages, or
physical storage cache and therefore is not canonical cold evidence.

All five runs:

- exited successfully;
- solved with index 4107 at field objects 41-50;
- preserved parity, center, scale, rotation, match statistics, and TAN/SIP
  coefficients;
- produced normalized WCS SHA-256
  `c739fee6dcd87f51aea100a7382b696ba1257fc3d9727db37146a9140ccd55ce`.

Raw FITS WCS hashes differ because runtime-dependent FITS metadata differs.
The normalized scientific header digest is identical.

## Decision supported by this evidence

The best continuation measurements remained approximately 2.67 seconds with
resident index pages and 3.10 seconds after the same index-only advisory
eviction. The continuation therefore failed its mandatory one-index causal
gate and was removed. The branch source and installed candidate were restored
to the Phase 08 boundary before publication.

The original raw evidence remains locally at:

```text
run_outputs/phase08_exact_control_20260810T151616Z
```

Its raw checksum manifest has SHA-256
`e2a5c9610ae76982c8a8df93dd5ee4062c7f3d4d58c35d2a99b29dff5bd3f868`.
The raw directory is intentionally not committed because it consists largely
of generated FITS products. The complete review-relevant measurements and
provenance are preserved in this capsule.
