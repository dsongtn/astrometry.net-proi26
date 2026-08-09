# Engineering Progress and Evaluation

Branch-local review copy for `test/inner-parallelism-and-scaling`. Do not merge
or publish this file with an official branch. Review responses belong in the
external `proi26-documentation/engineering-workspace`.

Status: `KEEP AS AUDIT RECORD`; current candidate status: `PROVISIONAL`

Development record last reconciled: 2026-08-03

Branch review identity added: 2026-08-09

## 1. Purpose

This is the ongoing engineering-accountability record for substantial private
development cycles. It records what was changed, why it was changed, what was
actually verified, what failed, and whether the decisions were effective.

It is not a replacement for source, private checkpoints, or raw run evidence.
When this document conflicts with an exact source tree or raw log, the source
and raw evidence take precedence.

This entry covers the reported six-hour development session on 2026-08-02.
The checkpointed implementation interval was 18:47 through 22:55 local time;
the last field run completed at 23:02. Earlier investigation and later result
analysis were outside that checkpoint interval.

## 2. Exact development identity

The final local checkpoint tree described below was subsequently captured as
the branch snapshot used for this review:

```text
branch: test/inner-parallelism-and-scaling
commit: adc027322558f05b6063f9e62a9a33ae903d3ea3
tree: c29c6039c1e25db6b4a546f1956fcad2516069cd
parent: 5bddecc9a8dc78f3ac27d6a258f92eb6cb992e5a
published base: e07c44c2cd7995c684bfa20ab3dc7b861038e1b8
```

The following normal-HEAD and checkpoint details describe the original
2026-08-02 development session and are retained as historical provenance.

The normal Git branch and HEAD did not move during the cycle:

```text
repository: /home/dsongt/HAW/26/PROINF/astrometry.net-solver-tasks
branch: test/inner-parallelism-and-scaling
HEAD: 5bddecc9a8dc78f3ac27d6a258f92eb6cb992e5a
HEAD tree: c985657faa32adc91fcf9767a2c367c524b4adee
normal commits created: 0
normal index/staging changes created by checkpoints: 0
```

The active source is an uncommitted overlay captured by private local-only
checkpoints. The exact final snapshot is:

```text
checkpoint series: final-performance
checkpoint: phase-12-two-stage-mapped-prime
checkpoint commit: 93d26b9871b75c92cfff035c9ead2a9d2b4c9e23
checkpoint tree: c29c6039c1e25db6b4a546f1956fcad2516069cd
predecessor checkpoint: f451af43e491e801e276bd268478d8b855b9ddc3
predecessor tree: c6e2613595d4df3e5767c7fc073b826d3358fce2
```

At reconciliation, the live source matched the Phase 12 checkpoint and
`git diff --check` was clean.

The clean-installed Phase 12 runtime at reconciliation was:

| Artifact | SHA-256 |
| --- | --- |
| `bin/astrometry-engine` | `8787f35543f516be6637f76d7239afb9648bbb58658df26d957c1b6007875b9a` |
| `bin/solve-field` | `8fbd2b20c2b7bd879c38233482c729abc273bc0e518ed20a7de1505d766aabac` |
| `etc/astrometry.cfg` | `23e9c877b20bad9955ae8d0a799ee97a6cf096a83207072ce074f71731aa6c2a` |

These hashes are historical identities for the Phase 12 field observations.
The shared workspace installation may later be rebuilt from another branch;
reviewers must not infer its current contents from this table.

The configuration retained 349 unique indexes in the required order with
manifest SHA-256:

```text
7ee54821c993f3a8ac16e212256d428bcf3a7b2f3308571866d3f2d88a27337e
```

## 3. Starting problem and hypothesis

The accepted architecture already separated compute workers from mapped-page
completion lanes, but the measured delivery path still behaved synchronously:

1. a CodeKD packet was planned;
2. an I/O lane issued page advice;
3. the same I/O lane immediately entered authoritative mapped population;
4. the lane could block in `MADV_POPULATE_READ` before preparing later work;
5. compute workers eventually exhausted ready packets and waited.

The cycle tested this hypothesis:

> More accurate readiness accounting, detached planning, bounded lookahead,
> fair scheduling, and a separate prime-before-populate pass can create useful
> compute/I/O overlap without changing scientific order or native fallback.

The intended counter movement was higher useful completion throughput, lower
system CPU and filesystem input for the same scientific progress, and fewer
periods with no compute-ready packet. Lower major-fault count by itself was
not the objective.

## 4. Scope of the cycle

From Phase 00 to Phase 12, the source-only delta covered 31 C, header, test,
and build-support files:

```text
5,219 insertions
1,237 deletions
```

Approximately 2,785 additions were deterministic tests or test support and
approximately 2,434 additions were production or production-facing code.
This is a large private candidate, even though the final Phase 11 to Phase 12
mechanism was much smaller:

```text
4 source/test files changed
253 insertions
6 deletions
```

The full cycle therefore must not be described as a small four-file change.
The four-file statement applies only to the final Phase 12 mechanism relative
to the Phase 11 checkpoint.

## 5. Checkpoint-by-checkpoint work

| Phase | Time | Main action | Decision or observed result |
| --- | --- | --- | --- |
| 00 baseline | 18:47 | Captured the exact pre-cycle overlay | Rollback boundary only |
| 01 delivery refresh batch | 19:07 | Refreshed validated mapped spans on the I/O lane; added bounded vectored population with exact per-span fallback | Retained as infrastructure; verbose APOD1 run remained unsolved |
| 02 readiness ledger | 20:06 | Added explicit claimable-ready counts and state-specific scheduler accounting for staged and helper work | Reduced blind scanning and made readiness auditable; broad and high-risk delta |
| 03 local validation | 20:31 | Added deterministic readiness tests and bounded all-span `MADV_WILLNEED` submission before authoritative population | Retained; label understates that production delivery code also changed |
| 04 dynamic delivery width | 20:50 | Made the payload service start at the requested live width with lifecycle-safe start/stop transitions | Functional, but the first following APOD1 run exposed severe starvation |
| 05 bounded retirement | 21:05 | Let the owner retire bounded verification slices under an already captured topology horizon rather than resubmitting equivalent delivery | Retained with focused tests |
| 05 logical-retirement correction | 21:09 | Repaired owner wake selection and atomic provider capability state | Corrected a work-conservation defect found after Phase 04 |
| 06 work-conserving wake | 21:10 | Validation checkpoint over the corrected Phase 05 source | No source delta from the preceding correction checkpoint |
| 07 detached CodeKD planning | 21:54 | Moved page-plan construction behind detached completion identities; added owner fallback and round-robin owner wakes | Architecturally important; large test and state-machine expansion |
| 08 CodeKD resident contract | 22:03 | Defined fully resident behavior separately from mapped population and expanded detached success/refusal/cancellation tests | Retained as correctness contract, not field-performance proof |
| 09 bounded lookahead | 22:23 | Published a current and one following descriptor wave and rotated submit/selection across owners | Later reduced because lookahead without a separate prime stage was not justified |
| 10 fair cross-owner | 22:31 | Preferred immediately ready local work and used fair foreign rotation only when local work was absent; disabled the second lookahead wave | Retained local-first/fairness rules; avoided leaving unproven speculation enabled |
| 11 validated retirement/width balance | 22:38 | Consolidated and reviewed the Phase 10 source | No source delta from Phase 10; status was validation, not a new mechanism |
| 11 clean-built fair delivery | 22:41 | Recorded the clean-build state | No source delta |
| 12 two-stage mapped prime | 22:55 | Split mapped delivery into prime/requeue and authoritative population passes; re-enabled one following CodeKD wave | Current `PROVISIONAL` candidate |

### 5.1 Important Phase 04 regression

The first quiet APOD1 W4 run after dynamic service-width work averaged only
89 percent CPU, completed 702 reported index-pass attempts, and timed out at
122.84 seconds. It had fewer faults and less input than later successful runs,
but it was slower because workers were not making enough scientific progress.

This was a useful falsification: low I/O counters are not automatically good.
The immediate response was to repair work-conserving owner wakes and bounded
retirement rather than interpreting the low fault count as success.

### 5.2 Why plain lookahead was not retained

Phase 09 doubled the logical packet window so later descriptors could enter
the pipeline. Phase 10 reduced the active window back to one because the
provider still performed prime and blocking population in one service pass.
That made the extra wave additional in-flight work without a reliable overlap
point.

Phase 12 re-enabled the second wave only after separating the advisory prime
pass from authoritative mapped population. This is a more defensible use of
lookahead, but its field effectiveness remains provisional.

## 6. Final Phase 12 mechanism

Phase 12 changed these source files relative to Phase 11:

```text
solver/solver_codekd_internal.h
util/fitsbin_payload_service.c
util/test_fitsbin_payload_io.c
util/test_fitsbin_payload_mapped.c
```

The mapped ticket lifecycle is now:

```text
SUBMITTED
  -> exact plan and mapping validation
  -> bounded WILLNEED/readahead prime
  -> requeue at the same priority when other work is queued
  -> later authoritative MADV_POPULATE_READ
  -> READY, FAILED, or CANCELLED
```

If no other ticket is queued, the lane falls through to population instead of
performing a useless requeue. Direct payload tickets are unchanged.

The same ticket retains its source lease, duplicated file descriptor,
mapping generation, sequence, plan state, and completion identity across the
requeue. Scientific objects do not enter the I/O queue. Native mapped access
remains the fallback and authoritative reader.

New aggregate counters distinguish:

```text
mapped_prime_passes
mapped_prime_requeues
```

The CodeKD logical lookahead remains strictly bounded at two waves. Descriptor
grain, sequence numbers, and owner-only retirement do not change.

## 7. Review findings during implementation

The most important review finding was that the new deterministic mapped-prime
test existed in `util/test_fitsbin_payload_mapped.c` but was initially missing
from the declarations consumed by the generated CuTest runner. The suite
reported 20 passing tests while silently omitting the new case.

The omission was repaired before the clean build and final checkpoint. The
suite then reported 21 tests. This incident matters because a green generated
test suite is not meaningful unless new tests are proven reachable.

The lifecycle review also checked:

- the source lease and duplicated descriptor survive both service passes;
- queued cancellation can remove a primed ticket before population;
- service stop drains queued tickets without requeue livelock;
- terminal accounting is applied once after the final pass;
- direct-ticket behavior is unchanged;
- the native mapped fallback remains available;
- descriptor and result retirement order remains owner-canonical.

No normal commit, branch, tag, or staging operation was used for phase
snapshots.

## 8. Verification actually performed

### 8.1 Build and focused tests

The authoritative clean build completed successfully:

```bash
cd /home/dsongt/HAW/26/PROINF
env SYSTEM_GSL=no workflow_scripts/build_solver_tasks.sh --clean
```

The following focused results were observed after the build:

| Check | Result |
| --- | --- |
| Payload I/O suite | `OK (21 tests)` |
| Solver permutation/enumeration suite | exit 0 |
| Staged scheduler suite | `INDEX_SHARD_STAGED_TEST_OK cases=19` |
| `git diff --check` | clean |
| Clean installed binaries | hashes recorded in Section 2 |

The interactive build and focused-test transcript was not saved as a dedicated
raw log. The result is therefore an observed engineering check, not an
independently replayable evidence package.

### 8.2 Deterministic mapped-prime test

The added test uses one provider lane and three tickets:

1. one blocking planner holds the lane;
2. a first mapped ticket and a second blocking ticket are queued;
3. releasing the blocker proves the first ticket is primed and requeued;
4. the second planner starts before the first ticket is populated;
5. the first ticket remains nonterminal and can be cancelled while queued;
6. the blocker and remaining ticket complete with single ownership transfer.

This proves queue ordering and cancellation mechanics. It does not prove disk
latency improvement.

### 8.3 Deep one-index fixture

A temporary one-index APOD1 fixture for `index-5206-47` completed the intended
deep CodeKD path with:

```text
hypotheses generated/executed/reduced: 37,944 / 37,944 / 37,944
maximum packet tasks: 8
I/O submitted/completed: 387 / 387
mapped prime passes: 387
mapped prime requeues: 85
mapped immediate-ready tickets: 13
WILLNEED failures: 0
ticket failures/cancellations: 0 / 0
page-pipeline fallback: 0
global integrity failures: 0
```

This showed that the mechanism was reachable and that two-wave publication
occurred. The fixture was warm and temporary, so its 0.06-second wall time was
not performance evidence. The `/tmp` evidence is no longer present.

### 8.4 APOD4 functional smoke

Temporary APOD4 W1 and W4 runs with `--scale-low 10` both solved with
`index-4116.fits` and produced the same normalized WCS signature:

```text
27e3f651184a6ee96f32504acc2e7ff6c18ac23b51b6fcc416f959dd5050e423
```

Observed immediate-cache walls were 1.26 seconds for W1 and 1.21 seconds for
W4. These were functional checks only. Their temporary outputs are no longer
present and they are not latency admission evidence.

### 8.5 Sanitizer gap

No new ASan, UBSan, or supported race-checker run was completed specifically
for the Phase 12 requeue path. Earlier project hardening does not substitute
for sanitizer evidence on newly changed ownership code. ThreadSanitizer was
previously unsupported on this host, but that status is not a Phase 12 pass.

## 9. Field evidence from the cycle

The five preserved APOD1 logs are summarized below. Directory phase labels do
not identify binaries; attribution uses checkpoint and file times and remains
qualified where a runtime hash was not embedded in the log.

| Evidence | Mode and likely source window | Result | Wall | CPU | Major faults | Filesystem input | Reported unsolved completions |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `apod1_phase01_delivery_w4_20260802T171800Z` | Verbose attribution after Phase 01 | timeout | 123.78 s | 287% | 425,377 | 57,765,568 | 884 |
| `final_phase03_apod1_w4_20260802T185459Z` | Quiet run after Phase 04 | timeout | 122.84 s | 89% | 177,337 | 22,389,488 | 702 |
| `final_phase03_apod1_w4_20260802T192220Z` | Quiet run after Phase 05/06 correction | solved | 101.44 s | 289% | 516,624 | 50,327,080 | 925 |
| `phase11_apod1_w4_20260802T205709Z` | Quiet run made after Phase 12 install | timeout | 123.64 s | 264% | 468,138 | 67,117,456 | 913 |
| `apod1_phase12_w4_20260802T210100Z` | Quiet Phase 12 run | solved | 95.39 s | 290% | 538,607 | 42,631,000 | 926 |

The first run was verbose and cannot support a latency comparison. The two
files named `final_phase03` were not source-enrolled and are historical
orientation only.

### 9.1 Mislabeled same-install evidence

The directory named `phase11_apod1_w4_20260802T205709Z` is not a valid Phase
11 versus Phase 12 comparator. The final installed binaries have unchanged
mtime and ctime from 22:52, while that run started around 22:57 and ended at
22:59. The final Phase 12 run started at 23:01. Both therefore used the same
installed candidate unless an unrecorded runtime substitution occurred.

This supplies two adjacent observations of the Phase 12 installation:

- one timed out at 123.64 seconds;
- one solved in 95.39 seconds.

The 28.25-second difference cannot be attributed to a source change. It
demonstrates substantial runtime variance from storage state, scheduling,
read amplification, or a combination of them.

### 9.2 Phase 12 positive result

The 95.39-second run:

- passed the minimum APOD1 W4 objective of less than 120 seconds;
- missed the 80-second challenge by 15.39 seconds;
- missed the 60-second stretch objective by 35.39 seconds;
- solved with `index-5204-23.fits` in field objects 21-30;
- reported the expected center, scale, parity, rotation, and field dimensions;
- matched the stored original reference's normalized TAN/SIP WCS exactly.

The normalized scientific WCS SHA-256 for both reference and Phase 12 output
was:

```text
83dd7fdd41b390b5c625ff5b2c2f330e82976c36e33699d4297144bb07f6ede4
```

### 9.3 Resource interpretation

Compared with the adjacent timeout from the same install, the successful run:

```text
wall time:                 -22.8 percent
user CPU:                   -7.5 percent
system CPU:                -28.6 percent
filesystem input:          -36.5 percent
reported completion rate:  +31.6 percent
major faults:              +15.1 percent
maximum RSS:                +6.9 percent
voluntary switches:         -1.1 percent
```

This is consistent with more precise I/O and less read amplification, not with
eliminating major faults. More faults accompanied less total input and faster
scientific progress.

The result does not show that the delivery problem is solved. Total CPU was
290 percent, consisting of about 202 percent user CPU and 89 percent system
CPU. About 1.1 of four cores remained unused on average, while 538,607 major
faults and 1,531,127 voluntary context switches remained. The provider still
uses advisory priming followed by a potentially blocking readiness authority.

The successful run also immediately followed a full-corpus traversal. Guest
page-cache dropping cannot clear hypervisor or host-storage caches. It is a
guest-page-cache-cold screen, not independently physical-cold admission.

## 10. Effectiveness and decision-quality assessment

### 10.1 Decisions that were responsible

- Private checkpoints preserved every major rollback boundary without
  modifying publishable history.
- W1 behavior, index order, solver mathematics, owner-private state, and
  reducer-only publication were kept outside the optimization.
- The Phase 04 low-CPU regression was treated as a failure rather than being
  misrepresented as an I/O improvement.
- Plain two-wave lookahead was disabled when its overlap mechanism was not
  convincing and re-enabled only after adding a separate prime stage.
- The final mechanism retained bounded tickets, native fallback, generation
  checks, cancellation, and owner-ordered retirement.
- A missing test-runner declaration was found and fixed before the final
  checkpoint.
- The final APOD1 WCS was checked against the original scientific reference,
  not merely against another parallel run.

### 10.2 Decisions or process weaknesses

- The Phase 00 to Phase 12 overlay was large: 31 files and more than 6,400
  changed lines. Several middle phases changed scheduler, delivery, and
  verification behavior before a source-enrolled field gate. This reduces
  causal clarity and increases review cost.
- Runtime hashes were not embedded in the five field logs. Binary attribution
  had to be reconstructed from install timestamps and checkpoint order.
- One output directory was named Phase 11 even though the Phase 12 install was
  already active. Phase labels should never substitute for recorded hashes.
- The final quiet run omitted internal prime/requeue counters, so the field
  result cannot quantify how much of the new path was used.
- Focused build/test results and temporary APOD4/deep-fixture outputs were not
  preserved in a durable evidence directory.
- The new Phase 12 ownership path did not receive a specific sanitizer run.
- No APOD5 field check or wider-width scaling check was performed for the
  final candidate.
- The adjacent same-binary timeout proves that one favorable screen is not a
  repeatability result.

### 10.3 Confidence by category

| Category | Confidence | Basis |
| --- | --- | --- |
| Scientific correctness for exercised APOD1/APOD4 cases | High | Same winner and exact normalized WCS; no scientific formulas changed |
| Phase 12 lifecycle correctness | Medium-high | Deterministic cancellation/requeue test and focused suites; no Phase 12 sanitizer evidence |
| Architectural soundness | Medium | Bounded and fallback-safe, but advisory priming does not guarantee completed storage I/O |
| APOD1 performance improvement | Low-to-medium | One 95.39-second solve and one 123.64-second timeout from the same install |
| Wider-worker scaling | Low | Not measured for Phase 12; dynamic lane width can increase cold concurrency |
| Release readiness | Low | Missing three-run field admission, APOD5, final sanitizer coverage, and reproducible provenance package |

## 11. Current decision

Phase 12 is classified as:

```text
correctness/lifecycle mechanism: KEEP AS PRIVATE CANDIDATE
performance status: PROVISIONAL
release status: NOT ADMITTED
```

Do not add another scheduler, lookahead, or payload mechanism before resolving
the current candidate's repeatability. Do not use the 95.39-second result as a
median or claim that the data-delivery problem is complete.

The first rollback boundary is Phase 11:

```text
refs/local-checkpoints/astrometry.net-solver-tasks/final-performance/
phase-11-validated-retirement-width-balance
```

The final Phase 12 delta can be reviewed independently because it affects four
source/test files relative to that checkpoint.

## 12. Single next gate

Run one later isolated APOD1 W4 guest-page-cache-cold screen using the exact
Phase 12 binary, configuration, manifest, and input identities, with no prior
corpus scan or unrelated storage workload. Capture the hashes in the run
directory before launch.

Decision rule:

- If it solves below 120 seconds with the expected normalized WCS, preserve
  Phase 12 and obtain only the remaining repetition needed for admission.
- If it times out or materially regresses, do not increase lookahead or lane
  width. Use source-exact counters to determine whether the failure was prime
  coverage, population wait, ready-queue starvation, or storage variation.
- Any next code change must target only the measured failed component and use
  Phase 11 as the rollback boundary.

## 13. Durable evidence

| Evidence | Purpose |
| --- | --- |
| `run_outputs/apod1_phase01_delivery_w4_20260802T171800Z/run.log` | Verbose Phase 01 delivery attribution |
| `run_outputs/final_phase03_apod1_w4_20260802T185459Z/run.log` | Low-CPU starvation regression |
| `run_outputs/final_phase03_apod1_w4_20260802T192220Z/run.log` | Restored work-conserving solved screen |
| `run_outputs/phase11_apod1_w4_20260802T205709Z/run.log` | Same-install timeout, despite misleading directory label |
| `run_outputs/apod1_phase12_w4_20260802T210100Z/run.log` | Phase 12 solved screen and resource metrics |
| `run_outputs/apod1_phase12_w4_20260802T210100Z/apod1.wcs` | Phase 12 scientific result |
| `run_outputs/v20-performance/20260726T231601Z-baseline-wcs/apod1/output/apod1.wcs` | Stored original APOD1 scientific reference |

Raw logs remain immutable. Any future conclusion must state whether it uses
the checkpoint identity, the installed binary identity, or only a historical
directory label.
