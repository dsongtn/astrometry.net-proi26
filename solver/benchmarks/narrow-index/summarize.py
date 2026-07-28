#!/usr/bin/env python3
"""Summarize a narrow-index matrix without third-party Python packages."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from fits_wcs_signature import SIGNATURE_METHOD, WcsSignatureError, wcs_signature


TIME_METRICS = (
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "cpu_seconds",
    "cpu_percent",
    "max_rss_kbytes",
    "major_faults",
    "minor_faults",
    "voluntary_context_switches",
    "involuntary_context_switches",
    "filesystem_inputs",
    "filesystem_outputs",
)
INDEX_COUNTS = (1, 2, 4, 8)
WORKER_COUNTS = (1, 2, 4)
CACHE_STATES = ("cold", "warm")
QUIET_SCHEDULER_CLASSIFICATION = "not-observable-quiet"
SCHEDULER_CLASSIFICATIONS = frozenset(
    {
        "serial",
        "dynamic-pool-lending",
    }
)
EVIDENCE_FILE_NAMES = frozenset(
    {
        "cache-control.json",
        "command.json",
        "solve.log",
        "solve.solved",
        "solve.wcs",
        "resource-usage.json",
        "wcsinfo.txt",
    }
)
SHA256_RE = re.compile(r"[0-9a-f]{64}")
FIELD_OBJECT_RANGE_RE = re.compile(
    r"^(?P<low>[1-9][0-9]*)-(?P<high>[1-9][0-9]*)$"
)
RANGE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class SummaryError(RuntimeError):
    """The result set is incomplete, inconsistent, or malformed."""


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            block = stream.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SummaryError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise SummaryError(f"top-level JSON in {path} is not an object")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(canonical_json(value))
    temporary.replace(path)


def median(values: Sequence[float]) -> float:
    if not values:
        raise SummaryError("median requested for an empty sequence")
    return float(statistics.median(values))


def median_mad(values: Sequence[float]) -> Tuple[float, float]:
    center = median(values)
    return center, median([abs(value - center) for value in values])


def percentile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise SummaryError("percentile requested for an empty sequence")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = probability * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(sorted_values[lower])
    fraction = position - lower
    return float(sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction)


def paired_bootstrap_median_ratio(
    baseline: Sequence[float],
    candidate: Sequence[float],
    samples: int,
    seed: int,
) -> Tuple[float, float, float]:
    if len(baseline) != len(candidate) or not baseline:
        raise SummaryError("paired bootstrap requires equal, non-empty samples")
    if any(value <= 0.0 for value in baseline) or any(value <= 0.0 for value in candidate):
        raise SummaryError("paired bootstrap wall times must be positive")
    ratios = [base / contender for base, contender in zip(baseline, candidate)]
    observed = median(ratios)
    generator = random.Random(seed)
    count = len(ratios)
    bootstrap = [
        median([ratios[generator.randrange(count)] for _ in range(count)])
        for _ in range(samples)
    ]
    bootstrap.sort()
    return observed, percentile(bootstrap, 0.025), percentile(bootstrap, 0.975)


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            delimiter="\t",
            fieldnames=list(fieldnames),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def condition_key(run: Mapping[str, Any]) -> Tuple[int, str, int, str]:
    return (
        int(run["index_count"]),
        str(run["range_id"]),
        int(run["workers"]),
        str(run["cache_state"]),
    )


def correctness_key(run: Mapping[str, Any]) -> Tuple[int, str, int]:
    return int(run["index_count"]), str(run["range_id"]), int(run["workers"])


def expected_scheduler_classifications(
    index_count: int, workers: int
) -> frozenset[str]:
    if workers == 1:
        return frozenset({"serial"})
    return frozenset({"dynamic-pool-lending"})


def scheduler_classification(
    run: Mapping[str, Any], measurement_mode: str
) -> Optional[str]:
    if measurement_mode == "timing":
        return QUIET_SCHEDULER_CLASSIFICATION
    assertions = run.get("assertions")
    if not isinstance(assertions, dict):
        return None
    value = assertions.get("scheduler_classification")
    return value if isinstance(value, str) else None


def scheduler_classification_summary(
    runs: Sequence[Mapping[str, Any]], measurement_mode: str
) -> str:
    values = {
        value
        for value in (
            scheduler_classification(run, measurement_mode) for run in runs
        )
        if value is not None
    }
    if not values:
        return "missing"
    if len(values) == 1:
        return next(iter(values))
    return "inconsistent:" + ",".join(sorted(values))


def collect_scheduler_classifications(
    runs: Sequence[Mapping[str, Any]], measurement_mode: str
) -> Tuple[Dict[Tuple[int, str, int], str], List[str]]:
    values: Dict[Tuple[int, str, int], set[str]] = {}
    errors: List[str] = []
    for run in runs:
        path = run.get("_record_path")
        try:
            index_count = int(run["index_count"])
            range_id = str(run["range_id"])
            workers = int(run["workers"])
        except (KeyError, TypeError, ValueError):
            errors.append(
                f"scheduler classification identity is malformed in {path}"
            )
            continue
        classification = scheduler_classification(run, measurement_mode)
        if classification is None:
            errors.append(
                f"missing assertions.scheduler_classification in {path}"
            )
            continue
        if measurement_mode == "timing":
            explicit = run.get("assertions", {}).get(
                "scheduler_classification"
            )
            if explicit not in (None, QUIET_SCHEDULER_CLASSIFICATION):
                errors.append(
                    f"quiet timing scheduler classification must be absent or "
                    f"{QUIET_SCHEDULER_CLASSIFICATION!r} in {path}, got "
                    f"{explicit!r}"
                )
        elif classification not in SCHEDULER_CLASSIFICATIONS:
            errors.append(
                f"unknown assertions.scheduler_classification "
                f"{classification!r} in {path}"
            )
            continue
        expected = (
            frozenset({QUIET_SCHEDULER_CLASSIFICATION})
            if measurement_mode == "timing"
            else expected_scheduler_classifications(
                index_count, workers
            )
        )
        if classification not in expected:
            errors.append(
                f"scheduler classification mismatch in {path}: "
                f"N={index_count} W={workers} allows "
                f"{sorted(expected)}, got {classification!r}"
            )
        values.setdefault((index_count, range_id, workers), set()).add(
            classification
        )

    classifications: Dict[Tuple[int, str, int], str] = {}
    for key, observed in sorted(values.items()):
        if len(observed) != 1:
            errors.append(
                f"scheduler classification is inconsistent across repetitions "
                f"and cache states for i{key[0]}/{key[1]}/W{key[2]}: "
                f"{sorted(observed)}"
            )
            continue
        classifications[key] = next(iter(observed))
    return classifications, errors


def scheduler_classification_rows(
    classifications: Mapping[Tuple[int, str, int], str]
) -> List[Dict[str, Any]]:
    return [
        {
            "index_count": index_count,
            "range_id": range_id,
            "workers": workers,
            "scheduler_classification": classification,
        }
        for (
            index_count,
            range_id,
            workers,
        ), classification in sorted(classifications.items())
    ]


def load_runs(root: Path, allow_incomplete: bool) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    suite = read_json(root / "suite.json")
    if not suite.get("complete") and not allow_incomplete:
        raise SummaryError("suite is not marked complete; use --allow-incomplete only for diagnosis")
    runs: List[Dict[str, Any]] = []
    for path in sorted((root / "runs").glob("r*/*/*/run.json")):
        run = read_json(path)
        run["_record_path"] = str(path)
        runs.append(run)
    expected = int(suite.get("run_count", 0))
    if len(runs) != expected and not allow_incomplete:
        raise SummaryError(f"expected {expected} run records, found {len(runs)}")
    if not runs:
        raise SummaryError("no run records found")
    return suite, runs


def validate_field_object_ranges(value: Any) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    if not isinstance(value, list) or len(value) < 2:
        return [], [
            "suite field_object_ranges must be an ordered array with at least two entries"
        ]
    range_ids: List[str] = []
    seen_ids: set[str] = set()
    seen_specs: set[str] = set()
    for position, entry in enumerate(value):
        name = f"suite field_object_ranges[{position}]"
        if not isinstance(entry, dict) or set(entry) != {"id", "objects"}:
            errors.append(
                f"{name} must contain exactly 'id' and 'objects'"
            )
            continue
        range_id = entry["id"]
        field_objects = entry["objects"]
        if not isinstance(range_id, str) or not RANGE_ID_RE.fullmatch(range_id):
            errors.append(f"{name}.id is invalid: {range_id!r}")
        elif range_id in seen_ids:
            errors.append(f"duplicate suite field-object range id {range_id!r}")
        else:
            seen_ids.add(range_id)
            range_ids.append(range_id)
        match = (
            FIELD_OBJECT_RANGE_RE.fullmatch(field_objects)
            if isinstance(field_objects, str)
            else None
        )
        if not match:
            errors.append(
                f"{name}.objects must be one explicit inclusive LOW-HIGH range"
            )
        else:
            low = int(match.group("low"))
            high = int(match.group("high"))
            if low > high:
                errors.append(f"{name}.objects has LOW greater than HIGH")
            elif field_objects in seen_specs:
                errors.append(
                    f"duplicate suite field-object range specification "
                    f"{field_objects!r}"
                )
            else:
                seen_specs.add(field_objects)
    return range_ids, errors


def validate_suite_shape(
    suite: Mapping[str, Any], range_ids: Sequence[str]
) -> List[str]:
    errors: List[str] = []
    repetitions = int(suite.get("repetitions", 0))
    if suite.get("range_ids") != list(range_ids):
        errors.append(
            "suite range_ids must exactly preserve field_object_ranges declaration order"
        )
    expected_cells = (
        len(INDEX_COUNTS)
        * len(range_ids)
        * len(WORKER_COUNTS)
        * len(CACHE_STATES)
    )
    expected_pairs = expected_cells // len(CACHE_STATES) * repetitions
    expected_runs = expected_cells * repetitions
    expected_metadata = {
        "cell_count": expected_cells,
        "pair_count": expected_pairs,
        "run_count": expected_runs,
        "index_counts": list(INDEX_COUNTS),
        "workers": list(WORKER_COUNTS),
        "cache_states": list(CACHE_STATES),
    }
    for key, expected in expected_metadata.items():
        if suite.get(key) != expected:
            errors.append(f"suite {key} mismatch: expected {expected}, got {suite.get(key)}")
    fingerprint = suite.get("fingerprint")
    if not isinstance(fingerprint, str) or not fingerprint:
        errors.append("suite fingerprint is missing")
    workload_fingerprint = suite.get("workload_fingerprint")
    if not isinstance(workload_fingerprint, str) or not SHA256_RE.fullmatch(
        workload_fingerprint
    ):
        errors.append("suite workload fingerprint is missing or malformed")
    if suite.get("measurement_mode") not in ("timing", "phase", "detailed"):
        errors.append(
            f"suite measurement_mode is invalid: {suite.get('measurement_mode')!r}"
        )
    if suite.get("resource_measurement") != "os.wait4":
        errors.append("suite resource_measurement must be os.wait4")
    launcher_median = suite.get("launcher_overhead_wall_median_seconds")
    if not isinstance(launcher_median, (int, float)) or launcher_median <= 0.0:
        errors.append("suite launcher-overhead median is missing or non-positive")
    return errors


def validate_runs(
    runs: Sequence[Mapping[str, Any]],
    repetitions: int,
    suite_fingerprint: str,
    workload_fingerprint: str,
    range_ids: Sequence[str],
) -> List[str]:
    errors: List[str] = []
    seen: set[Tuple[int, str, int, str, int]] = set()
    for run in runs:
        identity = (*condition_key(run), int(run["repetition"]))
        if identity in seen:
            errors.append(f"duplicate run identity {identity}")
        seen.add(identity)
        if not run.get("complete"):
            errors.append(f"incomplete run {run.get('_record_path')}")
        if run.get("suite_fingerprint") != suite_fingerprint:
            errors.append(f"suite fingerprint mismatch in {run.get('_record_path')}")
        if run.get("workload_fingerprint") != workload_fingerprint:
            errors.append(f"workload fingerprint mismatch in {run.get('_record_path')}")
        evidence = run.get("evidence_sha256")
        if not isinstance(evidence, dict) or not evidence:
            errors.append(f"missing evidence hashes in {run.get('_record_path')}")
        else:
            run_dir = Path(str(run.get("_record_path"))).parent
            required_evidence = {
                "command.json",
                "resource-usage.json",
                "solve.log",
                "wcsinfo.txt",
            }
            if run.get("cache_state") == "cold":
                required_evidence.add("cache-control.json")
            if run.get("solved"):
                required_evidence.update(("solve.solved", "solve.wcs"))
            missing_evidence = required_evidence - set(evidence)
            if missing_evidence:
                errors.append(
                    f"missing evidence hashes {sorted(missing_evidence)} in {run.get('_record_path')}"
                )
            for name, expected_digest in evidence.items():
                if name not in EVIDENCE_FILE_NAMES:
                    errors.append(f"unexpected evidence name {name!r} in {run.get('_record_path')}")
                    continue
                if not isinstance(expected_digest, str) or len(expected_digest) != 64:
                    errors.append(f"invalid evidence digest for {name!r} in {run.get('_record_path')}")
                    continue
                evidence_path = run_dir / name
                if not evidence_path.is_file():
                    errors.append(f"missing evidence file {evidence_path}")
                    continue
                if sha256_file(evidence_path) != expected_digest:
                    errors.append(f"evidence hash mismatch for {evidence_path}")
            resource_path = run_dir / "resource-usage.json"
            if resource_path.is_file():
                try:
                    resource_record = read_json(resource_path)
                except SummaryError as error:
                    errors.append(str(error))
                else:
                    if resource_record.get("measurement_method") != "os.wait4":
                        errors.append(
                            f"non-wait4 resource evidence in {run.get('_record_path')}"
                        )
                    timing_record = run.get("time", {})
                    for metric in (
                        *TIME_METRICS,
                        "returncode",
                        "wait_status",
                        "timed_out",
                    ):
                        if resource_record.get(metric) != timing_record.get(metric):
                            errors.append(
                                f"time.{metric} differs from raw resource evidence in "
                                f"{run.get('_record_path')}"
                            )
                    if resource_record.get("timed_out") is not False:
                        errors.append(f"timed-out run marked complete in {run.get('_record_path')}")
                    if resource_record.get("returncode") != 0:
                        errors.append(f"nonzero wait4 result in {run.get('_record_path')}")
            solved = run.get("solved")
            signature = run.get("correctness_signature")
            signature_method = run.get("correctness_signature_method")
            wcsinfo_sha256 = run.get("wcsinfo_sha256")
            raw_wcs_sha256 = run.get("raw_wcs_sha256")
            solved_path = run_dir / "solve.solved"
            wcs_path = run_dir / "solve.wcs"
            wcsinfo_path = run_dir / "wcsinfo.txt"
            if signature_method != SIGNATURE_METHOD:
                errors.append(
                    f"unsupported correctness signature method in {run.get('_record_path')}"
                )
            if not isinstance(wcsinfo_sha256, str) or not SHA256_RE.fullmatch(
                wcsinfo_sha256
            ):
                errors.append(f"invalid wcsinfo digest in {run.get('_record_path')}")
            elif wcsinfo_path.is_file() and sha256_file(wcsinfo_path) != wcsinfo_sha256:
                errors.append(f"wcsinfo digest mismatch in {run.get('_record_path')}")
            if not isinstance(solved, bool):
                errors.append(f"non-boolean solved status in {run.get('_record_path')}")
            elif solved:
                if not solved_path.is_file() or solved_path.stat().st_size == 0:
                    errors.append(f"missing non-empty solved marker in {run.get('_record_path')}")
                if not wcs_path.is_file():
                    errors.append(f"missing WCS in solved run {run.get('_record_path')}")
                if not isinstance(signature, str) or len(signature) != 64:
                    errors.append(f"invalid correctness signature in {run.get('_record_path')}")
                elif wcs_path.is_file():
                    try:
                        actual_signature = wcs_signature(wcs_path)
                    except WcsSignatureError as error:
                        errors.append(
                            f"cannot canonicalize WCS in {run.get('_record_path')}: {error}"
                        )
                    else:
                        if actual_signature != signature:
                            errors.append(
                                f"canonical WCS signature mismatch in {run.get('_record_path')}"
                            )
                if not isinstance(raw_wcs_sha256, str) or len(raw_wcs_sha256) != 64:
                    errors.append(f"invalid raw WCS digest in {run.get('_record_path')}")
                elif wcs_path.is_file() and sha256_file(wcs_path) != raw_wcs_sha256:
                    errors.append(f"raw WCS digest mismatch in {run.get('_record_path')}")
            else:
                if solved_path.exists() or wcs_path.exists():
                    errors.append(f"unsolved run retains solved artifacts in {run.get('_record_path')}")
                if signature != "UNSOLVED" or raw_wcs_sha256 is not None:
                    errors.append(f"inconsistent unsolved metadata in {run.get('_record_path')}")
                if wcsinfo_path.is_file() and wcsinfo_path.read_bytes() != b"UNSOLVED\n":
                    errors.append(f"invalid unsolved wcsinfo sentinel in {run.get('_record_path')}")
        if not run.get("assertions", {}).get("ok"):
            errors.append(f"failed candidate/pass assertion {run.get('_record_path')}")
        timing = run.get("time", {})
        for metric in TIME_METRICS:
            if not isinstance(timing.get(metric), (int, float)):
                errors.append(f"missing numeric time.{metric} in {run.get('_record_path')}")
        if timing.get("wall_seconds", 0) <= 0:
            errors.append(f"non-positive wall time in {run.get('_record_path')}")
    expected_identities = {
        (index_count, range_id, workers, cache_state, repetition)
        for index_count in INDEX_COUNTS
        for range_id in range_ids
        for workers in WORKER_COUNTS
        for cache_state in CACHE_STATES
        for repetition in range(1, repetitions + 1)
    }
    missing = sorted(expected_identities - seen)
    unexpected = sorted(seen - expected_identities)
    if missing:
        errors.append(
            f"matrix is missing {len(missing)} expected run identities; first={missing[:5]}"
        )
    if unexpected:
        errors.append(
            f"matrix has {len(unexpected)} unexpected run identities; first={unexpected[:5]}"
        )
    return errors


def validate_profile_evidence(
    runs: Sequence[Mapping[str, Any]], measurement_mode: str
) -> List[str]:
    errors: List[str] = []
    _, scheduler_errors = collect_scheduler_classifications(
        runs, measurement_mode
    )
    errors.extend(scheduler_errors)

    if measurement_mode == "timing":
        for run in runs:
            if run.get("profile_validation") != "not_observable_quiet":
                errors.append(
                    f"quiet timing profile-validation marker mismatch in "
                    f"{run.get('_record_path')}"
                )
            assertions = run.get("assertions", {})
            if assertions.get("observed_log_shape") != "not_observable_quiet":
                errors.append(
                    f"quiet timing log-shape marker mismatch in {run.get('_record_path')}"
                )
            if assertions.get("verbose_count") != 0:
                errors.append(
                    f"quiet timing run used verbosity in {run.get('_record_path')}"
                )
            profiles = run.get("profiles")
            if not isinstance(profiles, dict) or profiles:
                errors.append(
                    f"quiet timing run contains structured profiles in "
                    f"{run.get('_record_path')}"
                )
        return errors

    def require_keys(
        record: Mapping[str, Any], profile_name: str, keys: Sequence[str], path: Any
    ) -> None:
        missing = [key for key in keys if key not in record]
        if missing:
            errors.append(f"{profile_name} profile in {path} is missing keys {missing}")

    for run in runs:
        profiles = run.get("profiles", {})
        path = run.get("_record_path")
        if not isinstance(profiles, dict):
            errors.append(f"profiles object is malformed in {path}")
            continue
        if run.get("profile_validation") != "validated_verbose":
            errors.append(f"verbose profile-validation marker mismatch in {path}")
        assertions = run.get("assertions", {})
        expected_verbose = 1 if measurement_mode == "phase" else 2
        if assertions.get("observed_log_shape") != "validated_verbose":
            errors.append(f"verbose log-shape marker mismatch in {path}")
        if assertions.get("verbose_count") != expected_verbose:
            errors.append(
                f"verbosity mismatch in {path}: expected {expected_verbose}, "
                f"got {assertions.get('verbose_count')}"
            )
        for name in ("solve_field", "engine", "onefield"):
            records = profiles.get(name, []) if isinstance(profiles, dict) else []
            if not isinstance(records, list) or len(records) != 1:
                errors.append(f"expected one {name} profile in {path}, found {len(records) if isinstance(records, list) else 0}")
        field_records = profiles.get("onefield_field", [])
        if not isinstance(field_records, list) or not field_records:
            errors.append(f"expected one or more onefield_field profiles in {path}")
            field_records = []

        solve_records = profiles.get("solve_field", [])
        if isinstance(solve_records, list) and len(solve_records) == 1:
            require_keys(
                solve_records[0],
                "solve_field",
                ("source", "engine", "output", "total", "just_augment"),
                path,
            )
            if solve_records[0].get("just_augment") != 0:
                errors.append(f"solve-field unexpectedly ran in just-augment mode in {path}")

        engine_records = profiles.get("engine", [])
        if isinstance(engine_records, list) and len(engine_records) == 1:
            require_keys(
                engine_records[0],
                "engine",
                ("pool_start", "pool_stop", "engine_total", "solver_failed"),
                path,
            )
            if engine_records[0].get("solver_failed") != 0:
                errors.append(f"engine reported solver failure in {path}")

        for record in field_records:
            require_keys(
                record,
                "onefield_field",
                ("field", "read", "preprocess", "solver_run", "total", "failed"),
                path,
            )
            if record.get("failed") != 0:
                errors.append(f"onefield field execution failed in {path}")

        workers = int(run["workers"])
        index_count = int(run["index_count"])
        use_pthread = workers > 1
        onefield_records = profiles.get("onefield", []) if isinstance(profiles, dict) else []
        if isinstance(onefield_records, list) and len(onefield_records) == 1:
            require_keys(
                onefield_records[0],
                "onefield",
                (
                    "mode",
                    "candidates",
                    "serial_executed",
                    "acquire",
                    "solver",
                    "release",
                    "output",
                    "total",
                    "failed",
                    "cancelled",
                ),
                path,
            )
            expected_mode = (
                "serial"
                if not use_pthread
                else "pthread-dynamic-lending"
            )
            if onefield_records[0].get("mode") != expected_mode:
                errors.append(
                    f"onefield mode mismatch in {path}: expected {expected_mode}, "
                    f"got {onefield_records[0].get('mode')}"
                )
            if onefield_records[0].get("candidates") != index_count:
                errors.append(
                    f"onefield candidate mismatch in {path}: expected {index_count}, "
                    f"got {onefield_records[0].get('candidates')}"
                )
            if onefield_records[0].get("failed") != 0:
                errors.append(f"onefield reported failure in {path}")
            if onefield_records[0].get("cancelled") != 0:
                errors.append(f"onefield reported external cancellation in {path}")
            expected_serial_executed = (
                0 if use_pthread else len(field_records)
            )
            if (
                onefield_records[0].get("serial_executed")
                != expected_serial_executed
            ):
                errors.append(
                    f"serial execution/profile count mismatch in {path}: "
                    f"executed={onefield_records[0].get('serial_executed')} "
                    f"expected={expected_serial_executed} "
                    f"field_profiles={len(field_records)}"
                )
        if use_pthread:
            for name in (
                "index_shard_assist",
                "index_shard_reducer",
                "index_shard_context",
                "index_shard_solver",
                "index_shard_phase",
                "index_shard_resource",
            ):
                records = profiles.get(name, []) if isinstance(profiles, dict) else []
                if not isinstance(records, list) or len(records) != 1:
                    errors.append(
                        f"expected one {name} profile in parallel run {path}, "
                        f"found {len(records) if isinstance(records, list) else 0}"
                    )
            assist_records = profiles.get("index_shard_assist", [])
            if isinstance(assist_records, list) and len(assist_records) == 1:
                require_keys(
                    assist_records[0],
                    "index_shard_assist",
                    (
                        "generation",
                        "loans",
                        "notify_broadcasts",
                        "notify_skipped",
                        "waiters",
                    ),
                    path,
                )
                if assist_records[0].get("waiters") != 0:
                    errors.append(
                        f"index-shard assist aggregate ended with waiters in {path}"
                    )
                asserted_assist = assertions.get("assist_pass")
                if isinstance(asserted_assist, dict):
                    for key in (
                        "generation",
                        "loans",
                        "notify_broadcasts",
                        "notify_skipped",
                        "waiters",
                    ):
                        if assist_records[0].get(key) != asserted_assist.get(key):
                            errors.append(
                                f"assist assertion/profile {key} mismatch in "
                                f"{path}"
                            )
            resource_records = profiles.get("index_shard_resource", [])
            if isinstance(resource_records, list) and len(resource_records) == 1:
                require_keys(
                    resource_records[0],
                    "index_shard_resource",
                    (
                        "user",
                        "sys",
                        "minflt",
                        "majflt",
                        "nvcsw",
                        "nivcsw",
                        "inblock",
                        "oublock",
                    ),
                    path,
                )
            aux_records = profiles.get("index_shard_aux", [])
            if not isinstance(aux_records, list) or len(aux_records) > 1:
                errors.append(
                    f"expected at most one historical index_shard_aux profile "
                    f"in {path}, found "
                    f"{len(aux_records) if isinstance(aux_records, list) else 0}"
                )
            if isinstance(aux_records, list) and len(aux_records) == 1:
                require_keys(
                    aux_records[0],
                    "index_shard_aux",
                    (
                        "pool_workers",
                        "outer_workers",
                        "candidates",
                        "detailed",
                        "submitted",
                        "executed",
                        "idle_tasks",
                        "idle_queue_wait_sum",
                        "idle_work_wall_sum",
                        "lender_tasks",
                        "lender_queue_wait_sum",
                        "lender_work_wall_sum",
                        "owner_tasks",
                        "owner_queue_wait_sum",
                        "owner_work_wall_sum",
                        "total_queue_wait_sum",
                        "total_work_wall_sum",
                    ),
                    path,
                )
                expected_aux = {
                    "pool_workers": workers,
                    "outer_workers": workers,
                    "candidates": index_count,
                }
                for key, expected in expected_aux.items():
                    if aux_records[0].get(key) != expected:
                        errors.append(
                            f"aux {key} mismatch in {path}: expected {expected}, "
                            f"got {aux_records[0].get(key)}"
                        )
                expected_detailed = 1 if measurement_mode == "detailed" else 0
                if aux_records[0].get("detailed") != expected_detailed:
                    errors.append(
                        f"index_shard_aux detailed mismatch in {path}: "
                        f"expected {expected_detailed}, got {aux_records[0].get('detailed')}"
                    )
            for profile_name in ("index_shard_reducer", "index_shard_solver"):
                records = profiles.get(profile_name, [])
                if isinstance(records, list) and len(records) == 1:
                    if records[0].get("candidates") != index_count:
                        errors.append(
                            f"{profile_name} candidate mismatch in {path}: expected {index_count}, "
                            f"got {records[0].get('candidates')}"
                        )
            reducer_records = profiles.get("index_shard_reducer", [])
            if isinstance(reducer_records, list) and len(reducer_records) == 1:
                require_keys(
                    reducer_records[0],
                    "index_shard_reducer",
                    ("generation", "candidates", "calls", "work_wall_sum"),
                    path,
                )
            context_records = profiles.get("index_shard_context", [])
            if isinstance(context_records, list) and len(context_records) == 1:
                require_keys(
                    context_records[0],
                    "index_shard_context",
                    (
                        "generation",
                        "candidates",
                        "outer_workers",
                        "prepare_work_wall_sum",
                        "prepare_max",
                        "cleanup_work_wall_sum",
                        "cleanup_max",
                    ),
                    path,
                )
                if context_records[0].get("candidates") != index_count:
                    errors.append(f"index_shard_context candidate mismatch in {path}")
                if context_records[0].get("outer_workers") != workers:
                    errors.append(
                        f"index_shard_context outer-worker mismatch in {path}: "
                        f"expected {workers}, got "
                        f"{context_records[0].get('outer_workers')}"
                    )
            phase_records = profiles.get("index_shard_phase", [])
            if isinstance(phase_records, list) and len(phase_records) == 1:
                require_keys(
                    phase_records[0],
                    "index_shard_phase",
                    (
                        "executed",
                        "task_work_wall_sum",
                        "reset_work_wall_sum",
                        "reset_percent",
                        "acquire_work_wall_sum",
                        "acquire_percent",
                        "solve_work_wall_sum",
                        "solve_percent",
                        "analyze_work_wall_sum",
                        "analyze_percent",
                        "release_work_wall_sum",
                        "release_percent",
                        "other_work_wall_sum",
                        "other_percent",
                    ),
                    path,
                )
                outer_executed = phase_records[0].get("executed")
                if not isinstance(outer_executed, int) or not (
                    len(field_records) <= outer_executed <= index_count
                ):
                    errors.append(
                        f"outer/field profile count mismatch in {path}: "
                        f"outer_executed={outer_executed} field_profiles={len(field_records)} "
                        f"candidates={index_count}"
                    )
            solver_records = profiles.get("index_shard_solver", [])
            if isinstance(solver_records, list) and len(solver_records) == 1:
                require_keys(
                    solver_records[0],
                    "index_shard_solver",
                    (
                        "generation",
                        "candidates",
                        "detailed",
                        "failed",
                        "solver_run_work_wall_sum",
                        "codekd_work_wall_sum",
                        "resolve_work_wall_sum",
                        "reduction_ex_verify_hit_work_wall_sum",
                        "verify_hit_work_wall_sum",
                        "hypothesis_wave_elapsed_sum",
                        "observed_parallel",
                        "parallel_hypotheses",
                        "helper_tasks",
                        "helper_combinations",
                    ),
                    path,
                )
                expected_detailed = 1 if measurement_mode == "detailed" else 0
                if solver_records[0].get("detailed") != expected_detailed:
                    errors.append(
                        f"index_shard_solver detailed mismatch in {path}: "
                        f"expected {expected_detailed}, got {solver_records[0].get('detailed')}"
                    )
                if solver_records[0].get("failed") != 0:
                    errors.append(f"index-shard solver aggregate reported failure in {path}")
        if measurement_mode == "detailed":
            ab_phase_records = (
                profiles.get("solver_ab_phase", [])
                if isinstance(profiles, dict)
                else []
            )
            if not isinstance(ab_phase_records, list) or not ab_phase_records:
                errors.append(
                    f"missing exact solver_ab_phase telemetry in {path}"
                )
                ab_phase_records = []
            observed_modes = {
                mode: 0
                for mode in (
                    "assisted",
                    "empty",
                    "flattened-owner",
                    "native",
                )
            }
            for record in ab_phase_records:
                require_keys(
                    record,
                    "solver_ab_phase",
                    (
                        "index",
                        "object",
                        "phase",
                        "mode",
                        "combinations",
                        "codekd_queries",
                        "codekd_hits",
                        "candidates",
                        "verifications",
                        "blocks_planned",
                        "blocks_retired",
                        "blocks_owner",
                        "segments_retired",
                        "segment_payload_bytes",
                        "wall",
                        "user",
                        "system",
                        "major_faults",
                        "resource",
                    ),
                    path,
                )
                mode = record.get("mode")
                if mode not in observed_modes:
                    errors.append(
                        f"unknown exact AB phase mode {mode!r} in {path}"
                    )
                else:
                    observed_modes[mode] += 1
                if record.get("resource") not in {
                    "process-overlap",
                    "unavailable",
                }:
                    errors.append(
                        f"AB phase RUSAGE_SELF scope is not declared as "
                        f"process-overlap/unavailable in {path}"
                    )
            asserted_modes = assertions.get("phase_mode_counts")
            if isinstance(asserted_modes, dict) and asserted_modes != observed_modes:
                errors.append(
                    f"exact AB phase-mode assertion/profile mismatch in {path}: "
                    f"asserted={asserted_modes}, observed={observed_modes}"
                )
            records = profiles.get("solver_phase", []) if isinstance(profiles, dict) else []
            if not isinstance(records, list) or not records:
                errors.append(f"missing detailed solver_phase profile in {path}")
            else:
                if len(records) != len(field_records):
                    errors.append(
                        f"solver/field profile count mismatch in {path}: "
                        f"solver={len(records)} field={len(field_records)}"
                    )
                for record in records:
                    require_keys(
                        record,
                        "solver_phase",
                        (
                            "detailed",
                            "failed",
                            "solver_run_elapsed",
                            "codekd_work_wall_sum",
                            "codekd_calls",
                            "resolve_work_wall_sum",
                            "reduction_ex_verify_hit_work_wall_sum",
                            "resolve_calls",
                            "verify_hit_work_wall_sum",
                            "verify_calls",
                            "hypothesis_wave_elapsed_sum",
                        ),
                        path,
                    )
                    if record.get("detailed") != 1:
                        errors.append(f"non-detailed solver_phase record in detailed run {path}")
                    if record.get("failed") != 0:
                        errors.append(f"solver phase reported failure in {path}")
        elif isinstance(profiles, dict) and any(
            record.get("detailed") == 1 for record in profiles.get("solver_phase", [])
        ):
            errors.append(f"hot-call detailed instrumentation contaminated phase run {path}")
    return errors


def profile_values(run: Mapping[str, Any]) -> Dict[Tuple[str, str], float]:
    values: Dict[Tuple[str, str], float] = {}
    profiles = run.get("profiles", {})
    if not isinstance(profiles, dict):
        return values
    for profile_name, records in profiles.items():
        if not isinstance(records, list):
            continue
        per_metric: Dict[str, List[float]] = {}
        for record in records:
            if not isinstance(record, dict):
                continue
            for metric, value in record.items():
                if metric == "raw" or isinstance(value, bool) or not isinstance(value, (int, float)):
                    continue
                per_metric.setdefault(metric, []).append(float(value))
        for metric, observations in per_metric.items():
            values[(str(profile_name), metric)] = sum(observations)
        values[(str(profile_name), "record_count")] = float(len(records))
    return values


def condition_group_sort_key(
    item: Tuple[
        Tuple[int, str, int, str],
        Sequence[Mapping[str, Any]],
    ],
    range_ids: Sequence[str],
) -> Tuple[int, int, int, int]:
    index_count, range_id, workers, cache_state = item[0]
    range_positions = {name: position for position, name in enumerate(range_ids)}
    cache_positions = {
        name: position for position, name in enumerate(CACHE_STATES)
    }
    return (
        index_count,
        range_positions.get(range_id, len(range_positions)),
        workers,
        cache_positions.get(cache_state, len(cache_positions)),
    )


def make_condition_rows(
    groups: Mapping[Tuple[int, str, int, str], Sequence[Mapping[str, Any]]],
    range_ids: Sequence[str],
    measurement_mode: str = "timing",
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for (index_count, range_id, workers, cache_state), observations in sorted(
        groups.items(),
        key=lambda item: condition_group_sort_key(item, range_ids),
    ):
        row: Dict[str, Any] = {
            "index_count": index_count,
            "range_id": range_id,
            "workers": workers,
            "cache_state": cache_state,
            "scheduler_classification": scheduler_classification_summary(
                observations, measurement_mode
            ),
            "n": len(observations),
            "solved_count": sum(bool(run["solved"]) for run in observations),
            "winner_positions": ",".join(
                str(position)
                for position in sorted(
                    {
                        run.get("winner", {}).get("position_one_based")
                        for run in observations
                        if run.get("winner", {}).get("position_one_based") is not None
                    }
                )
            ),
        }
        for metric in TIME_METRICS:
            center, spread = median_mad(
                [float(run["time"][metric]) for run in observations]
            )
            row[f"{metric}_median"] = center
            row[f"{metric}_mad"] = spread
        rows.append(row)
    return rows


def make_phase_rows(
    groups: Mapping[Tuple[int, str, int, str], Sequence[Mapping[str, Any]]],
    range_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for (index_count, range_id, workers, cache_state), observations in sorted(
        groups.items(),
        key=lambda item: condition_group_sort_key(item, range_ids),
    ):
        per_metric: Dict[Tuple[str, str], List[float]] = {}
        for run in observations:
            for key, value in profile_values(run).items():
                per_metric.setdefault(key, []).append(value)
        for (profile, metric), values in sorted(per_metric.items()):
            center, spread = median_mad(values)
            rows.append(
                {
                    "index_count": index_count,
                    "range_id": range_id,
                    "workers": workers,
                    "cache_state": cache_state,
                    "profile": profile,
                    "metric": metric,
                    "n": len(values),
                    "median": center,
                    "mad": spread,
                }
            )
    return rows


def make_helper_rows(
    groups: Mapping[Tuple[int, str, int, str], Sequence[Mapping[str, Any]]],
    measurement_mode: str,
    range_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    metric_names = (
        ("index_shard_assist", "loans"),
        ("index_shard_assist", "notify_broadcasts"),
        ("index_shard_assist", "notify_skipped"),
        ("index_shard_solver", "observed_parallel"),
        ("index_shard_solver", "parallel_hypotheses"),
        ("index_shard_solver", "helper_tasks"),
        ("index_shard_solver", "helper_combinations"),
    )
    for (index_count, range_id, workers, cache_state), observations in sorted(
        groups.items(),
        key=lambda item: condition_group_sort_key(item, range_ids),
    ):
        if workers == 1:
            continue
        observable = measurement_mode != "timing"
        classification = scheduler_classification_summary(
            observations, measurement_mode
        )
        run_classifications = [
            scheduler_classification(run, measurement_mode)
            for run in observations
        ]
        per_run = (
            [profile_values(run) for run in observations]
            if observable
            else [{} for _ in observations]
        )
        row: Dict[str, Any] = {
            "index_count": index_count,
            "range_id": range_id,
            "workers": workers,
            "cache_state": cache_state,
            "scheduler_classification": classification,
            "runs": len(observations),
            "runs_dynamic_pool_lending": sum(
                value == "dynamic-pool-lending"
                for value in run_classifications
            ),
            "runs_with_lane_publication": sum(
                bool(
                    run.get("assertions", {}).get(
                        "assist_lane_publications"
                    )
                )
                for run in observations
            ),
            "runs_with_loans": sum(
                values.get(
                    ("index_shard_assist", "loans"), 0.0
                )
                > 0.0
                for values in per_run
            ),
            "runs_with_assisted_phase": sum(
                run.get("assertions", {})
                .get("phase_mode_counts", {})
                .get("assisted", 0)
                > 0
                for run in observations
            ),
            "runs_with_flattened_owner_phase": sum(
                run.get("assertions", {})
                .get("phase_mode_counts", {})
                .get("flattened-owner", 0)
                > 0
                for run in observations
            ),
            "runs_with_exact_phase_modes": sum(
                sum(
                    int(value)
                    for value in run.get("assertions", {})
                    .get("phase_mode_counts", {})
                    .values()
                    if isinstance(value, int)
                )
                > 0
                for run in observations
            ),
            "runs_with_observed_assistance": sum(
                values.get(
                    ("index_shard_solver", "observed_parallel"), 0.0
                )
                > 0.0
                or values.get(
                    ("index_shard_assist", "loans"), 0.0
                )
                > 0.0
                for values in per_run
            ),
        }
        for profile, metric in metric_names:
            values = [
                record[(profile, metric)]
                for record in per_run
                if (profile, metric) in record
            ]
            row[f"{metric}_median"] = median(values) if values else None
            row[f"{metric}_max"] = max(values) if values else None
        if not observable:
            flag = "NOT_OBSERVABLE_QUIET"
        elif classification != "dynamic-pool-lending":
            flag = "INVALID_SCHEDULER_CLASSIFICATION"
        elif row["runs_with_observed_assistance"] > 0:
            flag = "DYNAMIC_ASSIST_OBSERVED"
        elif row["runs_with_flattened_owner_phase"] > 0:
            flag = "FLATTENED_OWNER_ONLY"
        elif measurement_mode == "detailed":
            flag = "NATIVE_OR_EMPTY_ONLY"
        else:
            flag = "NO_LOAN_OBSERVED"
        row["flag"] = flag
        rows.append(row)
    return rows


def make_correctness_rows(
    runs: Sequence[Mapping[str, Any]],
    repetitions: int,
    range_ids: Sequence[str],
) -> Tuple[List[Dict[str, Any]], Dict[Tuple[int, str, int], bool], List[str]]:
    signatures: Dict[Tuple[int, str, int], set[str]] = {}
    solved_states: Dict[Tuple[int, str, int], set[bool]] = {}
    counts: Dict[Tuple[int, str, int], int] = {}
    for run in runs:
        key = correctness_key(run)
        signatures.setdefault(key, set()).add(str(run["correctness_signature"]))
        solved_states.setdefault(key, set()).add(bool(run["solved"]))
        counts[key] = counts.get(key, 0) + 1
    rows: List[Dict[str, Any]] = []
    validity: Dict[Tuple[int, str, int], bool] = {}
    errors: List[str] = []
    for index_count, range_id in (
        (index_count, range_id)
        for index_count in INDEX_COUNTS
        for range_id in range_ids
    ):
        baseline_key = (index_count, range_id, 1)
        baseline_signatures = signatures.get(baseline_key, set())
        baseline_deterministic = len(baseline_signatures) == 1
        if not baseline_deterministic:
            errors.append(
                f"W1 correctness is non-deterministic for "
                f"i{index_count}/{range_id}: "
                f"{sorted(baseline_signatures)}"
            )
        for workers in (1, 2, 4):
            key = (index_count, range_id, workers)
            candidate_signatures = signatures.get(key, set())
            expected_count = repetitions * 2
            complete = counts.get(key, 0) == expected_count
            deterministic = len(candidate_signatures) == 1
            exact_match = (
                baseline_deterministic
                and deterministic
                and candidate_signatures == baseline_signatures
            )
            valid = complete and exact_match and len(solved_states.get(key, set())) == 1
            validity[key] = valid
            status = "PASS" if valid else "FAIL"
            if not valid:
                errors.append(
                    f"correctness {status} for "
                    f"i{index_count}/{range_id}/W{workers}: "
                    f"baseline={sorted(baseline_signatures)} candidate={sorted(candidate_signatures)} "
                    f"runs={counts.get(key, 0)}/{expected_count}"
                )
            rows.append(
                {
                    "index_count": index_count,
                    "range_id": range_id,
                    "workers": workers,
                    "runs": counts.get(key, 0),
                    "expected_runs": expected_count,
                    "solved_states": ",".join(str(int(value)) for value in sorted(solved_states.get(key, set()))),
                    "w1_signature_count": len(baseline_signatures),
                    "signature_count": len(candidate_signatures),
                    "exact_w1_match": int(exact_match),
                    "status": status,
                }
            )
    return rows, validity, errors


def correctness_signatures(
    runs: Sequence[Mapping[str, Any]],
) -> Dict[Tuple[int, str, int], Tuple[str, ...]]:
    signatures: Dict[Tuple[int, str, int], set[str]] = {}
    for run in runs:
        signatures.setdefault(correctness_key(run), set()).add(
            str(run.get("correctness_signature"))
        )
    return {key: tuple(sorted(values)) for key, values in signatures.items()}


def validate_companion_roots(
    current_root: Path,
    current_suite: Mapping[str, Any],
    current_runs: Sequence[Mapping[str, Any]],
    companion_paths: Sequence[Path],
) -> Tuple[Dict[str, Any], List[str]]:
    current_mode = str(current_suite.get("measurement_mode"))
    expected_fingerprint = str(current_suite.get("workload_fingerprint", ""))
    current_signatures = correctness_signatures(current_runs)
    records: List[Dict[str, Any]] = []
    errors: List[str] = []
    seen_roots: set[Path] = set()
    seen_modes: set[str] = set()
    current_classifications, current_scheduler_errors = (
        collect_scheduler_classifications(current_runs, current_mode)
    )
    observed_classifications: Dict[
        str, Dict[Tuple[int, str, int], str]
    ] = {}
    if (
        current_mode in ("phase", "detailed")
        and not current_scheduler_errors
    ):
        observed_classifications[current_mode] = current_classifications

    for configured_path in companion_paths:
        try:
            root = configured_path.expanduser().resolve(strict=True)
        except OSError as error:
            errors.append(f"companion root does not resolve: {configured_path}: {error}")
            continue
        if root == current_root:
            errors.append(f"companion root is the current suite root: {root}")
            continue
        if root in seen_roots:
            errors.append(f"duplicate companion root: {root}")
            continue
        seen_roots.add(root)

        record: Dict[str, Any] = {"root": str(root)}
        companion_errors: List[str] = []
        try:
            suite, runs = load_runs(root, allow_incomplete=False)
            provenance = read_json(root / "provenance" / "provenance.json")
        except SummaryError as error:
            companion_errors.append(str(error))
            record.update({"mode": None, "valid": False, "errors": companion_errors})
            records.append(record)
            errors.extend(f"companion {root}: {item}" for item in companion_errors)
            continue

        mode = str(suite.get("measurement_mode"))
        fingerprint = suite.get("workload_fingerprint")
        record.update({"mode": mode, "workload_fingerprint": fingerprint})
        duplicate_mode = mode in seen_modes
        if mode == current_mode:
            companion_errors.append(
                f"companion mode duplicates current measurement mode {current_mode!r}"
            )
        if duplicate_mode:
            companion_errors.append(f"duplicate companion measurement mode {mode!r}")
        seen_modes.add(mode)
        if fingerprint != expected_fingerprint:
            companion_errors.append(
                "workload fingerprint differs from the current suite: "
                f"{fingerprint!r} != {expected_fingerprint!r}"
            )
        if provenance.get("workload_fingerprint") != fingerprint:
            companion_errors.append(
                "companion provenance and suite workload fingerprints differ"
            )

        companion_range_ids, companion_range_errors = (
            validate_field_object_ranges(suite.get("field_object_ranges"))
        )
        companion_errors.extend(companion_range_errors)
        companion_errors.extend(
            validate_suite_shape(suite, companion_range_ids)
        )
        companion_errors.extend(
            validate_runs(
                runs,
                int(suite.get("repetitions", 0)),
                str(suite.get("fingerprint", "")),
                str(fingerprint or ""),
                companion_range_ids,
            )
        )
        companion_errors.extend(validate_profile_evidence(runs, mode))
        companion_classifications, companion_scheduler_errors = (
            collect_scheduler_classifications(runs, mode)
        )
        if (
            mode in ("phase", "detailed")
            and not duplicate_mode
            and not companion_scheduler_errors
        ):
            observed_classifications[mode] = companion_classifications
        _, _, correctness_errors = make_correctness_rows(
            runs,
            int(suite.get("repetitions", 0)),
            companion_range_ids,
        )
        companion_errors.extend(correctness_errors)
        if correctness_signatures(runs) != current_signatures:
            companion_errors.append(
                "correctness signatures differ from the current suite"
            )

        record.update(
            {
                "scheduler_classifications": scheduler_classification_rows(
                    companion_classifications
                ),
                "valid": not companion_errors,
                "errors": companion_errors,
            }
        )
        records.append(record)
        errors.extend(f"companion {root}: {error}" for error in companion_errors)

    if companion_paths and current_mode == "timing" and "phase" not in seen_modes:
        errors.append("a timing-suite comparison requires a phase companion")

    scheduler_consistency = "NOT_CHECKED"
    if (
        "phase" in observed_classifications
        and "detailed" in observed_classifications
    ):
        if (
            observed_classifications["phase"]
            != observed_classifications["detailed"]
        ):
            errors.append(
                "phase/detailed scheduler classifications differ"
            )
            scheduler_consistency = "FAIL"
        else:
            scheduler_consistency = "PASS"

    if "phase" in observed_classifications:
        classification_source = "phase"
        surfaced_classifications = observed_classifications["phase"]
    elif current_mode in observed_classifications:
        classification_source = current_mode
        surfaced_classifications = observed_classifications[current_mode]
    elif "detailed" in observed_classifications:
        classification_source = "detailed"
        surfaced_classifications = observed_classifications["detailed"]
    else:
        classification_source = QUIET_SCHEDULER_CLASSIFICATION
        surfaced_classifications = current_classifications

    if not companion_paths:
        status = "NOT_CHECKED"
    elif errors:
        status = "FAIL"
    else:
        status = "PASS"
    return (
        {
            "status": status,
            "current_mode": current_mode,
            "workload_fingerprint": expected_fingerprint,
            "scheduler_classification_source": classification_source,
            "scheduler_classification_consistency": scheduler_consistency,
            "scheduler_classifications": scheduler_classification_rows(
                surfaced_classifications
            ),
            "timing_acceptance_ready": (
                current_mode == "timing" and status == "PASS" and "phase" in seen_modes
            ),
            "companions": records,
        },
        errors,
    )


def derived_seed(base_seed: int, key: Sequence[Any]) -> int:
    material = canonical_json([base_seed, *key])
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big")


def paired_by_repetition(
    observations: Sequence[Mapping[str, Any]], workers: int
) -> Tuple[List[Mapping[str, Any]], List[Mapping[str, Any]]]:
    baseline = {
        int(run["repetition"]): run for run in observations if int(run["workers"]) == 1
    }
    candidate = {
        int(run["repetition"]): run for run in observations if int(run["workers"]) == workers
    }
    repetitions = sorted(set(baseline) & set(candidate))
    return [baseline[item] for item in repetitions], [candidate[item] for item in repetitions]


def make_parallel_rows(
    runs: Sequence[Mapping[str, Any]],
    validity: Mapping[Tuple[int, str, int], bool],
    repetitions: int,
    bootstrap_samples: int,
    bootstrap_seed: int,
    measurement_mode: str,
    range_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for index_count, range_id, cache_state in (
        (index_count, range_id, cache_state)
        for index_count in INDEX_COUNTS
        for range_id in range_ids
        for cache_state in CACHE_STATES
    ):
        observations = [
            run
            for run in runs
            if int(run["index_count"]) == index_count
            and str(run["range_id"]) == range_id
            and str(run["cache_state"]) == cache_state
        ]
        for workers in (2, 4):
            baseline, candidate = paired_by_repetition(observations, workers)
            pair_count = len(baseline)
            correctness_ok = validity.get((index_count, range_id, workers), False)
            complete = pair_count == repetitions
            row: Dict[str, Any] = {
                "index_count": index_count,
                "range_id": range_id,
                "cache_state": cache_state,
                "workers": workers,
                "pairs": pair_count,
                "correctness": "PASS" if correctness_ok else "FAIL",
                "measurement_mode": measurement_mode,
            }
            if not complete:
                row["flag"] = "INVALID_INCOMPLETE"
                rows.append(row)
                continue
            base_wall = [float(run["time"]["wall_seconds"]) for run in baseline]
            candidate_wall = [float(run["time"]["wall_seconds"]) for run in candidate]
            speedup, ci_low, ci_high = paired_bootstrap_median_ratio(
                base_wall,
                candidate_wall,
                bootstrap_samples,
                derived_seed(
                    bootstrap_seed, (index_count, range_id, cache_state, workers)
                ),
            )
            base_work = [
                float(run["time"]["user_seconds"]) + float(run["time"]["system_seconds"])
                for run in baseline
            ]
            candidate_work = [
                float(run["time"]["user_seconds"]) + float(run["time"]["system_seconds"])
                for run in candidate
            ]
            work_ratios = [
                contender / base if base > 0.0 else float("nan")
                for base, contender in zip(base_work, candidate_work)
            ]
            finite_work_ratios = [value for value in work_ratios if math.isfinite(value)]
            row.update(
                {
                    "w1_wall_median": median(base_wall),
                    "candidate_wall_median": median(candidate_wall),
                    "paired_speedup_median": speedup,
                    "speedup_ci95_low": ci_low,
                    "speedup_ci95_high": ci_high,
                    "paired_wall_delta_median": median(
                        [contender - base for base, contender in zip(base_wall, candidate_wall)]
                    ),
                    "candidate_over_w1_cpu_work_median": (
                        median(finite_work_ratios) if finite_work_ratios else None
                    ),
                    "w1_major_faults_median": median(
                        [float(run["time"]["major_faults"]) for run in baseline]
                    ),
                    "candidate_major_faults_median": median(
                        [float(run["time"]["major_faults"]) for run in candidate]
                    ),
                }
            )
            if not correctness_ok:
                flag = "INVALID_CORRECTNESS"
            elif measurement_mode != "timing":
                flag = "DIAGNOSTIC_ONLY"
            elif ci_low > 1.0:
                flag = "PARALLEL_FASTER"
            elif ci_high < 1.0:
                flag = "PARALLEL_SLOWER"
            else:
                flag = "INCONCLUSIVE"
            row["flag"] = flag
            rows.append(row)
    return rows


def make_crossover_rows(
    condition_rows: Sequence[Mapping[str, Any]],
    parallel_rows: Sequence[Mapping[str, Any]],
    range_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    condition_lookup = {
        (
            int(row["index_count"]),
            str(row["range_id"]),
            str(row["cache_state"]),
            int(row["workers"]),
        ): row
        for row in condition_rows
    }
    rows: List[Dict[str, Any]] = []
    for index_count, range_id, cache_state in (
        (index_count, range_id, cache_state)
        for index_count in INDEX_COUNTS
        for range_id in range_ids
        for cache_state in CACHE_STATES
    ):
        comparisons = [
            row
            for row in parallel_rows
            if int(row["index_count"]) == index_count
            and str(row["range_id"]) == range_id
            and str(row["cache_state"]) == cache_state
        ]
        medians: Dict[int, Optional[float]] = {}
        for workers in (1, 2, 4):
            condition = condition_lookup.get(
                (index_count, range_id, cache_state, workers)
            )
            medians[workers] = (
                float(condition["wall_seconds_median"]) if condition is not None else None
            )
        if any(value is None for value in medians.values()):
            rows.append(
                {
                    "index_count": index_count,
                    "range_id": range_id,
                    "cache_state": cache_state,
                    "w1_wall_median": medians[1],
                    "w2_wall_median": medians[2],
                    "w4_wall_median": medians[4],
                    "raw_fastest_workers": None,
                    "crossover_flag": "INVALID",
                    "evidence_selected_workers": 1,
                }
            )
            continue
        complete_medians = {workers: float(value) for workers, value in medians.items() if value is not None}
        best_workers = min(complete_medians, key=lambda workers: complete_medians[workers])
        flags = {str(row["flag"]) for row in comparisons}
        confident = [row for row in comparisons if row.get("flag") == "PARALLEL_FASTER"]
        if any(flag.startswith("INVALID") for flag in flags):
            crossover = "INVALID"
            selected = 1
        elif "DIAGNOSTIC_ONLY" in flags:
            crossover = "DIAGNOSTIC_ONLY"
            selected = None
        elif confident:
            selected_row = max(
                confident, key=lambda row: float(row["paired_speedup_median"])
            )
            selected = int(selected_row["workers"])
            crossover = "PARALLEL"
        elif comparisons and all(row.get("flag") == "PARALLEL_SLOWER" for row in comparisons):
            crossover = "SERIAL"
            selected = 1
        else:
            crossover = "UNCERTAIN"
            selected = 1
        rows.append(
            {
                "index_count": index_count,
                "range_id": range_id,
                "cache_state": cache_state,
                "w1_wall_median": complete_medians[1],
                "w2_wall_median": complete_medians[2],
                "w4_wall_median": complete_medians[4],
                "raw_fastest_workers": best_workers,
                "crossover_flag": crossover,
                "evidence_selected_workers": selected,
            }
        )
    return rows


def format_number(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_report(
    path: Path,
    suite: Mapping[str, Any],
    condition_rows: Sequence[Mapping[str, Any]],
    crossover_rows: Sequence[Mapping[str, Any]],
    correctness_rows: Sequence[Mapping[str, Any]],
    helper_rows: Sequence[Mapping[str, Any]],
    companion_validation: Mapping[str, Any],
    errors: Sequence[str],
    bootstrap_samples: int,
) -> None:
    lines = [
        "# Narrow-index benchmark report",
        "",
        f"Measurement mode: `{suite.get('measurement_mode', 'timing')}`. "
        f"Repetitions per condition/state: {suite.get('repetitions')}. "
        f"Paired bootstrap resamples: {bootstrap_samples}.",
        "",
        "Median absolute deviation (MAD) is the unscaled median absolute deviation. "
        "A timing crossover is asserted only when the paired 95% bootstrap interval "
        "lies wholly on one side of 1.0 and exact W1 correctness passes.",
        "",
        (
            "Quiet timing mode intentionally has no verbose profile or observed "
            "pass-shape evidence; use the companion phase suite for those diagnostics."
            if suite.get("measurement_mode") == "timing"
            else "This verbose suite is diagnostic-only and cannot set the production crossover."
        ),
        "",
        "Companion validation: `{status}`. Timing acceptance ready: `{ready}`.".format(
            status=companion_validation.get("status"),
            ready=str(bool(companion_validation.get("timing_acceptance_ready"))).lower(),
        ),
        "Scheduler classification evidence source: `{source}`.".format(
            source=companion_validation.get(
                "scheduler_classification_source",
                QUIET_SCHEDULER_CLASSIFICATION,
            )
        ),
        "",
        "## Crossover decisions",
        "",
        "| Indexes | Field-object range ID | OS page cache | W1 median (s) | W2 median (s) | W4 median (s) | Decision | Selected W |",
        "|---:|---|---|---:|---:|---:|---|---:|",
    ]
    for row in crossover_rows:
        lines.append(
            "| {index_count} | {range_id} | {cache_state} | {w1} | {w2} | {w4} | {flag} | {selected} |".format(
                index_count=row["index_count"],
                range_id=row["range_id"],
                cache_state=row["cache_state"],
                w1=format_number(row["w1_wall_median"]),
                w2=format_number(row["w2_wall_median"]),
                w4=format_number(row["w4_wall_median"]),
                flag=row["crossover_flag"],
                selected=format_number(row["evidence_selected_workers"]),
            )
        )
    lines.extend(
        [
            "",
            "## Correctness",
            "",
            "| Indexes | Field-object range ID | Workers | Runs | Exact W1 match | Status |",
            "|---:|---|---:|---:|---:|---|",
        ]
    )
    for row in correctness_rows:
        lines.append(
            f"| {row['index_count']} | {row['range_id']} | {row['workers']} | "
            f"{row['runs']} | {row['exact_w1_match']} | {row['status']} |"
        )
    lines.extend(
        [
            "",
            "## Dynamic pool lending",
            "",
            "| Indexes | Field-object range ID | OS page cache | Workers | Scheduler | Observed parallel runs | Flag |",
            "|---:|---|---|---:|---|---:|---|",
        ]
    )
    for row in helper_rows:
        lines.append(
            f"| {row['index_count']} | {row['range_id']} | {row['cache_state']} | "
            f"{row['workers']} | {row['scheduler_classification']} | "
            f"{row['runs_with_observed_assistance']} | "
            f"{row['flag']} |"
        )
    if errors:
        lines.extend(["", "## Validation failures", ""])
        lines.extend(f"- {error}" for error in errors)
    else:
        lines.extend(["", "All completeness, pass-shape, and exact W1 correctness checks passed."])
    lines.extend(
        [
            "",
            "Machine-readable detail is in `condition_summary.tsv`, `phase_summary.tsv`, "
            "`helper_usage.tsv`, `parallel_vs_w1.tsv`, `crossover.tsv`, and `correctness.tsv`.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument(
        "--companion",
        action="append",
        default=[],
        type=Path,
        metavar="ROOT",
        help=(
            "validate a companion timing/phase/detailed suite against the same "
            "workload; repeat for more than one companion"
        ),
    )
    args = parser.parse_args(argv)
    root = args.root.expanduser().resolve(strict=True)
    suite, runs = load_runs(root, args.allow_incomplete)
    repetitions = int(suite.get("repetitions", 0))
    if repetitions < 7 and not args.allow_incomplete:
        raise SummaryError(f"suite has only {repetitions} repetitions; at least 7 are required")

    provenance = read_json(root / "provenance" / "provenance.json")
    config = provenance.get("config", {})
    bootstrap_samples = int(config.get("bootstrap_samples", 20000))
    bootstrap_seed = int(config.get("bootstrap_seed", 20260721))
    if bootstrap_samples < 5000:
        raise SummaryError("bootstrap_samples must be >= 5000")
    measurement_mode = str(suite.get("measurement_mode", "timing"))

    range_ids, range_errors = validate_field_object_ranges(
        suite.get("field_object_ranges")
    )
    validation_errors = range_errors + validate_suite_shape(suite, range_ids)
    workload_fingerprint = str(suite.get("workload_fingerprint", ""))
    if provenance.get("workload_fingerprint") != workload_fingerprint:
        validation_errors.append(
            "provenance and suite workload fingerprints differ"
        )
    validation_errors.extend(
        validate_runs(
            runs,
            repetitions,
            str(suite.get("fingerprint", "")),
            workload_fingerprint,
            range_ids,
        )
    )
    validation_errors.extend(validate_profile_evidence(runs, measurement_mode))
    scheduler_classifications, _ = collect_scheduler_classifications(
        runs, measurement_mode
    )
    companion_validation, companion_errors = validate_companion_roots(
        root, suite, runs, args.companion
    )
    validation_errors.extend(companion_errors)
    groups: Dict[Tuple[int, str, int, str], List[Mapping[str, Any]]] = {}
    for run in runs:
        groups.setdefault(condition_key(run), []).append(run)
    condition_rows = make_condition_rows(
        groups, range_ids, measurement_mode
    )
    phase_rows = make_phase_rows(groups, range_ids)
    helper_rows = make_helper_rows(groups, measurement_mode, range_ids)
    correctness_rows, validity, correctness_errors = make_correctness_rows(
        runs, repetitions, range_ids
    )
    all_errors = validation_errors + correctness_errors
    parallel_rows = make_parallel_rows(
        runs,
        validity,
        repetitions,
        bootstrap_samples,
        bootstrap_seed,
        measurement_mode,
        range_ids,
    )
    crossover_rows = make_crossover_rows(
        condition_rows, parallel_rows, range_ids
    )

    summary_dir = root / "summary"
    condition_fields = [
        "index_count",
        "range_id",
        "workers",
        "cache_state",
        "scheduler_classification",
        "n",
        "solved_count",
        "winner_positions",
    ] + [f"{metric}_{suffix}" for metric in TIME_METRICS for suffix in ("median", "mad")]
    write_tsv(summary_dir / "condition_summary.tsv", condition_rows, condition_fields)
    write_tsv(
        summary_dir / "phase_summary.tsv",
        phase_rows,
        (
            "index_count",
            "range_id",
            "workers",
            "cache_state",
            "profile",
            "metric",
            "n",
            "median",
            "mad",
        ),
    )
    helper_fields = (
        "index_count",
        "range_id",
        "workers",
        "cache_state",
        "scheduler_classification",
        "runs",
        "runs_dynamic_pool_lending",
        "runs_with_lane_publication",
        "runs_with_loans",
        "runs_with_assisted_phase",
        "runs_with_flattened_owner_phase",
        "runs_with_exact_phase_modes",
        "runs_with_observed_assistance",
        "loans_median",
        "loans_max",
        "notify_broadcasts_median",
        "notify_broadcasts_max",
        "notify_skipped_median",
        "notify_skipped_max",
        "observed_parallel_median",
        "observed_parallel_max",
        "parallel_hypotheses_median",
        "parallel_hypotheses_max",
        "helper_tasks_median",
        "helper_tasks_max",
        "helper_combinations_median",
        "helper_combinations_max",
        "flag",
    )
    write_tsv(summary_dir / "helper_usage.tsv", helper_rows, helper_fields)
    parallel_fields = (
        "index_count",
        "range_id",
        "cache_state",
        "workers",
        "pairs",
        "correctness",
        "measurement_mode",
        "w1_wall_median",
        "candidate_wall_median",
        "paired_speedup_median",
        "speedup_ci95_low",
        "speedup_ci95_high",
        "paired_wall_delta_median",
        "candidate_over_w1_cpu_work_median",
        "w1_major_faults_median",
        "candidate_major_faults_median",
        "flag",
    )
    write_tsv(summary_dir / "parallel_vs_w1.tsv", parallel_rows, parallel_fields)
    write_tsv(
        summary_dir / "crossover.tsv",
        crossover_rows,
        (
            "index_count",
            "range_id",
            "cache_state",
            "w1_wall_median",
            "w2_wall_median",
            "w4_wall_median",
            "raw_fastest_workers",
            "crossover_flag",
            "evidence_selected_workers",
        ),
    )
    write_tsv(
        summary_dir / "correctness.tsv",
        correctness_rows,
        (
            "index_count",
            "range_id",
            "workers",
            "runs",
            "expected_runs",
            "solved_states",
            "w1_signature_count",
            "signature_count",
            "exact_w1_match",
            "status",
        ),
    )
    write_report(
        summary_dir / "report.md",
        suite,
        condition_rows,
        crossover_rows,
        correctness_rows,
        helper_rows,
        companion_validation,
        all_errors,
        bootstrap_samples,
    )
    write_json(
        summary_dir / "summary.json",
        {
            "schema_version": 1,
            "valid": not all_errors,
            "measurement_mode": measurement_mode,
            "workload_fingerprint": workload_fingerprint,
            "scheduler_classifications": scheduler_classification_rows(
                scheduler_classifications
            ),
            "companion_validation": companion_validation,
            "bootstrap_samples": bootstrap_samples,
            "validation_errors": all_errors,
            "condition_summary": condition_rows,
            "phase_summary": phase_rows,
            "helper_usage": helper_rows,
            "parallel_vs_w1": parallel_rows,
            "crossover": crossover_rows,
            "correctness": correctness_rows,
        },
    )
    print(summary_dir / "report.md")
    return 0 if not all_errors else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SummaryError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
