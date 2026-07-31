#!/usr/bin/env python3
"""Run the reproducible narrow-index solve-field benchmark matrix.

This runner deliberately executes one solve-field process at a time.  Every
condition is an OS-page-cache cold/reuse pair, and the reuse process starts
immediately after the page-cache-cold process.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import select
import shlex
import signal
import statistics
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from fits_wcs_signature import SIGNATURE_METHOD, WcsSignatureError, wcs_signature


SCHEMA_VERSION = 1
INDEX_COUNTS = (1, 2, 4, 8)
WORKER_COUNTS = (1, 2, 4)
CACHE_STATES = ("cold", "warm")
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
PERFORMANCE_ENV = {
    "ATLAS_NUM_THREADS": "1",
    "BLAS_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
    "MKL_DYNAMIC": "FALSE",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
RESERVED_ARGUMENTS = {
    "--axy",
    "--backend-config",
    "--batch",
    "--config",
    "--continue",
    "--depth",
    "--dir",
    "--files-on-stdin",
    "--index-dir",
    "--index-file",
    "--out",
    "--overwrite",
    "--skip-solved",
    "--solved",
    "--temp-axy",
    "--timestamp",
    "--wcs",
    "--verbose",
    "-D",
    "-J",
    "-K",
    "-M",
    "-N",
    "-O",
    "-P",
    "-R",
    "-U",
    "-Z",
    "-B",
    "-b",
    "-d",
    "-f",
    "-i",
    "-j",
    "-k",
    "-n",
    "-o",
    "-v",
}
DISABLE_ONLY_ARGUMENTS = {
    "--corr",
    "--index-xyls",
    "--keep-xylist",
    "--kmz",
    "--match",
    "--new-fits",
    "--pnm",
    "--rdls",
    "--scamp",
    "--scamp-config",
    "--scamp-ref",
}
SUBMIT_RE = re.compile(
    r"\[index-shard\] pthread-pool submit "
    r"workers=(?P<workers>\d+) pool_workers=(?P<pool_workers>\d+) "
    r"candidates=(?P<candidates>\d+) "
    r"engine_pass=(?P<engine_pass>\d+) "
    r"depth_index=(?P<depth_index>\d+) "
    r"scale_index=(?P<scale_index>\d+) "
    r"startobj=(?P<startobj>\d+) endobj=(?P<endobj>\d+) "
    r"scheduler=ordered chunk=1 "
    r"inner_scheduler=(?P<inner_scheduler>[A-Za-z0-9_-]+) "
    r"mmap_pass=(?P<mmap_pass>\d+) "
    r"mmap_advice=(?P<mmap_advice>[A-Za-z0-9_-]+)"
)
ASSIST_LANE_PUBLISH_RE = re.compile(
    r"\[index-shard\] assist-lane state=(?P<state>publish|unpublish) "
    r"owner_worker=(?P<owner_worker>\d+) "
    r"index_order=(?P<index_order>\d+)"
)
ASSIST_LANE_LOAN_RE = re.compile(
    r"\[index-shard\] assist-lane state=(?P<state>join|leave) "
    r"helper_worker=(?P<helper_worker>\d+) "
    r"helper_slot=(?P<helper_slot>\d+) "
    r"index_order=(?P<index_order>\d+)"
    r"(?: rc=(?P<rc>-?\d+))?"
)
ASSIST_PASS_RE = re.compile(
    r"\[index-shard\] assist-pass "
    r"generation=(?P<generation>\d+) "
    r"loans=(?P<loans>\d+) "
    r"notify_broadcasts=(?P<notify_broadcasts>\d+) "
    r"notify_skipped=(?P<notify_skipped>\d+) "
    r"waiters=(?P<waiters>\d+)"
)
AB_PHASE_RE = re.compile(
    r"^\[solver-ab-phase\] (?P<payload>.*)$",
    re.MULTILINE,
)
GEOMETRY_FALLBACK_RE = re.compile(
    r"\[solver-geometry\] mode=legacy "
    r"reason=(?P<reason>budget|allocation|empty-range|estimate)"
)
GEOMETRY_SHARED_RE = re.compile(
    r"\[solver-geometry\] mode=shared-readonly "
)
REDUCED_WINNER_RE = re.compile(r"\[index-shard\] reduce solved index_order=(?P<order>\d+)")
SERIAL_TRY_RE = re.compile(r"^Trying index (?P<path>.+)\.\.\.$")
SERIAL_VERIFY_RE = re.compile(
    r"^Verifying WCS with index \d+ of \d+ \((?P<path>.+)\)$"
)
SOLVED_INDEX_RE = re.compile(r"^Field \d+: solved with index ")
FIELD_OBJECT_RANGE_RE = re.compile(
    r"^(?P<low>[1-9][0-9]*)-(?P<high>[1-9][0-9]*)$"
)
RANGE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
KEY_VALUE_RE = re.compile(r"(?P<key>[A-Za-z][A-Za-z0-9_.-]*)=(?P<value>[^\s]+)")


class HarnessError(RuntimeError):
    """A benchmark precondition or evidence assertion failed."""


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds")


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            block = stream.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(canonical_json(value))
    os.replace(temporary, path)


def read_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise HarnessError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise HarnessError(f"top-level JSON value in {path} must be an object")
    return value


def require_regular_file(value: Any, name: str, executable: bool = False) -> Path:
    if not isinstance(value, str) or not value:
        raise HarnessError(f"{name} must be a non-empty path string")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise HarnessError(f"{name} must be absolute: {path}")
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise HarnessError(f"{name} does not resolve: {path}: {error}") from error
    if not resolved.is_file():
        raise HarnessError(f"{name} is not a regular file: {resolved}")
    if executable and not os.access(resolved, os.X_OK):
        raise HarnessError(f"{name} is not executable: {resolved}")
    return resolved


def parse_field_object_range(value: Any, range_id: str) -> Tuple[str, int, int]:
    if not isinstance(value, str):
        raise HarnessError(f"field-object range {range_id!r} must be a string")
    match = FIELD_OBJECT_RANGE_RE.fullmatch(value)
    if not match:
        raise HarnessError(
            f"field-object range {range_id!r} must be one explicit inclusive "
            f"LOW-HIGH range, got {value!r}"
        )
    low = int(match.group("low"))
    high = int(match.group("high"))
    if low > high:
        raise HarnessError(
            f"field-object range {range_id!r} has LOW greater than HIGH: {value}"
        )
    return value, low - 1, high


def parse_field_object_ranges(
    value: Any,
) -> List[Tuple[str, str, int, int]]:
    if not isinstance(value, list) or len(value) < 2:
        raise HarnessError(
            "field_object_ranges must be an ordered array with at least two entries"
        )
    parsed: List[Tuple[str, str, int, int]] = []
    seen_ids: set[str] = set()
    seen_specs: set[str] = set()
    for position, entry in enumerate(value):
        name = f"field_object_ranges[{position}]"
        if not isinstance(entry, dict) or set(entry) != {"id", "objects"}:
            raise HarnessError(
                f"{name} must be an object containing exactly 'id' and 'objects'"
            )
        range_id = entry["id"]
        if not isinstance(range_id, str) or not RANGE_ID_RE.fullmatch(range_id):
            raise HarnessError(
                f"{name}.id must match {RANGE_ID_RE.pattern!r}, got {range_id!r}"
            )
        if range_id in seen_ids:
            raise HarnessError(f"duplicate field-object range id {range_id!r}")
        field_objects, startobj, endobj = parse_field_object_range(
            entry["objects"], range_id
        )
        if field_objects in seen_specs:
            raise HarnessError(
                f"duplicate field-object range specification {field_objects!r}"
            )
        seen_ids.add(range_id)
        seen_specs.add(field_objects)
        parsed.append((range_id, field_objects, startobj, endobj))
    return parsed


def validate_common_args(value: Any) -> List[str]:
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise HarnessError("common_args must be an array of strings")
    args = list(value)
    for argument in args:
        if not argument.startswith("--") or argument == "--":
            raise HarnessError(
                "every common_args entry must be one long option token; "
                "attach option values with '=' so positional input files cannot be injected: "
                f"{argument!r}"
            )
        option = argument.split("=", 1)[0]
        if option in RESERVED_ARGUMENTS:
            raise HarnessError(
                f"common_args contains harness-owned option {argument!r}; remove it"
            )
        if option in DISABLE_ONLY_ARGUMENTS:
            if "=" not in argument:
                raise HarnessError(
                    f"common_args option {argument!r} must attach its value with '='"
                )
            output_value = argument.split("=", 1)[1]
            if output_value != "none":
                raise HarnessError(
                    f"common_args may only set {option} to 'none'; per-run output paths are harness-owned"
                )
    return args


def workload_config(config: Mapping[str, Any]) -> Dict[str, Any]:
    """Return the exact cross-mode workload definition.

    A companion suite must differ only in where evidence is written and in the
    requested instrumentation mode.  Everything else, including the human
    label, cache control, repetitions and common solver arguments, is part of
    the workload identity.
    """
    normalized = dict(config)
    normalized.pop("output_root", None)
    normalized.pop("measurement_mode", None)
    return normalized


def compute_workload_fingerprint(
    config: Mapping[str, Any],
    provenance: Mapping[str, Any],
    plan: Sequence[Mapping[str, Any]],
) -> str:
    value = {
        "config": workload_config(config),
        "solve_field": provenance["solve_field"],
        "wcsinfo": provenance["wcsinfo"],
        "input": provenance["input"],
        "manifest_file": provenance["manifest_file"],
        "indexes": provenance["indexes"],
        "pipeline_files": provenance["pipeline_files"],
        "taskset_binary": provenance["taskset_binary"],
        "launcher_probe_binary": provenance["launcher_probe_binary"],
        "harness": provenance["harness"],
        "plan": list(plan),
    }
    return sha256_bytes(canonical_json(value))


def load_indexes(path: Path) -> List[Path]:
    indexes: List[Path] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if any(character in line for character in "*?[]{}"):
            raise HarnessError(
                f"index manifest {path}:{line_number} contains a glob; list an exact path"
            )
        index = require_regular_file(line, f"index manifest line {line_number}")
        indexes.append(index)
    if len(indexes) != 8:
        raise HarnessError(f"index manifest must contain exactly 8 paths, found {len(indexes)}")
    if len(set(indexes)) != len(indexes):
        raise HarnessError("index manifest contains duplicate resolved paths")
    return indexes


def load_pipeline_files(value: Any, solve_field: Path) -> List[Path]:
    if not isinstance(value, list) or not value:
        raise HarnessError("pipeline_files must be a non-empty array of absolute paths")
    files = [
        require_regular_file(item, f"pipeline_files[{position}]")
        for position, item in enumerate(value)
    ]
    if len(files) != len(set(files)):
        raise HarnessError("pipeline_files contains duplicate resolved paths")
    engine_path = solve_field.parent / "astrometry-engine"
    try:
        expected_engine = engine_path.resolve(strict=True)
    except OSError as error:
        raise HarnessError(
            f"cannot resolve the astrometry-engine adjacent to solve_field: {engine_path}: {error}"
        ) from error
    if expected_engine not in files:
        raise HarnessError(
            f"pipeline_files must include the engine selected by solve-field: {expected_engine}"
        )
    if not os.access(expected_engine, os.X_OK):
        raise HarnessError(f"astrometry-engine is not executable: {expected_engine}")
    return files


def file_record(path: Path) -> Dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "mode": stat.st_mode,
    }


def stat_identity(path: Path) -> Dict[str, int]:
    stat = path.stat()
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
    }


def run_capture(command: Sequence[str], timeout: float = 30.0) -> Dict[str, Any]:
    started = utc_now()
    try:
        completed = subprocess.run(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
        return {
            "command": list(command),
            "started_utc": started,
            "finished_utc": utc_now(),
            "returncode": completed.returncode,
            "output": completed.stdout,
        }
    except (OSError, subprocess.TimeoutExpired) as error:
        return {
            "command": list(command),
            "started_utc": started,
            "finished_utc": utc_now(),
            "error": repr(error),
        }


def collect_provenance(
    config_path: Path,
    config: Mapping[str, Any],
    solve_field: Path,
    wcsinfo: Path,
    input_path: Path,
    manifest_path: Path,
    indexes: Sequence[Path],
    pipeline_files: Sequence[Path],
    taskset_binary: Path,
    launcher_probe_binary: Path,
) -> Dict[str, Any]:
    script_dir = Path(__file__).resolve().parent
    runner = Path(__file__).resolve()
    summarizer = script_dir / "summarize.py"
    wcs_signer = script_dir / "fits_wcs_signature.py"
    index_records = [file_record(path) for path in indexes]
    provenance: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "captured_utc": utc_now(),
        "config_path": str(config_path),
        "config_sha256": sha256_file(config_path),
        "config_file": file_record(config_path),
        "config": dict(config),
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_file": file_record(manifest_path),
        "solve_field": file_record(solve_field),
        "wcsinfo": file_record(wcsinfo),
        "input": file_record(input_path),
        "indexes": index_records,
        "pipeline_files": [file_record(path) for path in pipeline_files],
        "taskset_binary": file_record(taskset_binary),
        "launcher_probe_binary": file_record(launcher_probe_binary),
        "harness": {
            "run_matrix.py": file_record(runner),
            "summarize.py": file_record(summarizer),
            "fits_wcs_signature.py": file_record(wcs_signer),
        },
        "python": {
            "executable": sys.executable,
            "version": sys.version,
        },
        "uname": list(os.uname()),
        "performance_environment": dict(PERFORMANCE_ENV),
        "ambient_environment": {
            name: os.environ.get(name)
            for name in sorted(set(PERFORMANCE_ENV) | {"LANG", "LC_ALL", "PATH", "TZ"})
        },
        "commands": {},
        "system_files": {},
    }
    commands = {
        "solve_field_version": [str(solve_field), "--version"],
        "git_revision": ["git", "-C", str(script_dir.parents[2]), "rev-parse", "HEAD"],
        "git_status": [
            "git",
            "-C",
            str(script_dir.parents[2]),
            "status",
            "--short",
            "--branch",
        ],
        "lscpu": ["lscpu"],
        "numactl": ["numactl", "--hardware"],
        "ldd_solve_field": ["ldd", str(solve_field)],
        "taskset_version": [str(taskset_binary), "--version"],
    }
    for name, command in commands.items():
        provenance["commands"][name] = run_capture(command)
    for position, path in enumerate(pipeline_files):
        provenance["commands"][f"ldd_pipeline_{position:02d}"] = run_capture(
            ["ldd", str(path)]
        )
    for system_path in (
        Path("/proc/cpuinfo"),
        Path("/proc/meminfo"),
        Path("/proc/cmdline"),
        Path("/sys/devices/system/cpu/online"),
        Path("/sys/devices/system/cpu/smt/active"),
        Path("/sys/devices/system/cpu/intel_pstate/no_turbo"),
    ):
        try:
            provenance["system_files"][str(system_path)] = system_path.read_text(
                encoding="utf-8", errors="replace"
            )
        except OSError as error:
            provenance["system_files"][str(system_path)] = {"error": repr(error)}
    return provenance


def condition_plan(
    range_ids: Sequence[str], repetitions: int, seed: int, rotation_stride: int
) -> List[Dict[str, Any]]:
    base: List[Tuple[int, str, int]] = [
        (index_count, range_id, workers)
        for index_count in INDEX_COUNTS
        for range_id in range_ids
        for workers in WORKER_COUNTS
    ]
    random.Random(seed).shuffle(base)
    if math.gcd(rotation_stride, len(base)) != 1:
        raise HarnessError(
            f"rotation_stride={rotation_stride} must be coprime to {len(base)} conditions"
        )
    plan: List[Dict[str, Any]] = []
    sequence = 0
    for repetition in range(1, repetitions + 1):
        offset = ((repetition - 1) * rotation_stride) % len(base)
        rotated = base[offset:] + base[:offset]
        for order_in_repetition, (index_count, range_id, workers) in enumerate(
            rotated, 1
        ):
            sequence += 1
            plan.append(
                {
                    "sequence": sequence,
                    "repetition": repetition,
                    "order_in_repetition": order_in_repetition,
                    "index_count": index_count,
                    "range_id": range_id,
                    "workers": workers,
                }
            )
    return plan


def safe_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-")
    if not label:
        raise HarnessError(f"cannot make a filesystem label from {value!r}")
    return label


def pair_relative_path(item: Mapping[str, Any]) -> Path:
    condition = (
        f"o{item['order_in_repetition']:02d}_i{item['index_count']}_"
        f"{safe_label(str(item['range_id']))}_w{item['workers']}"
    )
    return Path("runs") / f"r{item['repetition']:02d}" / condition


def benchmark_environment(workers: int) -> Dict[str, str]:
    environment = os.environ.copy()
    environment.update(PERFORMANCE_ENV)
    environment.update(
        {
            "ASTROMETRY_INDEX_SHARD_WORKERS": str(workers),
            "LANG": "C",
            "LC_ALL": "C",
            "TZ": "UTC",
        }
    )
    return environment


def build_solve_command(
    solve_field: Path,
    run_dir: Path,
    input_path: Path,
    indexes: Sequence[Path],
    field_objects: str,
    common_args: Sequence[str],
    measurement_mode: str,
) -> List[str]:
    command = [str(solve_field)]
    if measurement_mode == "phase":
        command.append("--verbose")
    elif measurement_mode == "detailed":
        command.extend(("--verbose", "--verbose"))
    command.extend(
        (
            "--overwrite",
            "--no-plots",
            "--config",
            "none",
            "--dir",
            str(run_dir),
            "--out",
            "solve",
            "--depth",
            field_objects,
        )
    )
    for index in indexes:
        command.extend(("--index-file", str(index)))
    command.extend(common_args)
    command.append(str(input_path))
    return command


def render_command(command: Sequence[str]) -> str:
    return shlex.join(list(command))


def run_drop_caches(cache_config: Mapping[str, Any]) -> Dict[str, Any]:
    method = cache_config.get("method")
    started = utc_now()
    monotonic_start = time.perf_counter_ns()
    if hasattr(os, "sync"):
        os.sync()
    if method == "proc":
        if os.geteuid() != 0:
            raise HarnessError("cache_control method 'proc' requires root")
        proc_path = Path("/proc/sys/vm/drop_caches")
        try:
            proc_path.write_text("3\n", encoding="ascii")
        except OSError as error:
            raise HarnessError(f"failed to write {proc_path}: {error}") from error
        result: Dict[str, Any] = {"method": method, "path": str(proc_path), "returncode": 0}
    elif method == "command":
        command = cache_config.get("command")
        stdin_text = cache_config.get("stdin", "")
        if not isinstance(command, list) or not command or not all(
            isinstance(item, str) and item for item in command
        ):
            raise HarnessError("cache_control.command must be a non-empty string array")
        if not isinstance(stdin_text, str):
            raise HarnessError("cache_control.stdin must be a string")
        try:
            completed = subprocess.run(
                command,
                input=stdin_text,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise HarnessError(f"cache-control command failed to execute: {error}") from error
        result = {
            "method": method,
            "command": command,
            "stdin_sha256": sha256_bytes(stdin_text.encode("utf-8")),
            "returncode": completed.returncode,
            "output": completed.stdout,
        }
        if completed.returncode != 0:
            raise HarnessError(
                f"cache-control command returned {completed.returncode}: {completed.stdout.strip()}"
            )
    else:
        raise HarnessError(
            "cache_control.method must be 'proc' or 'command'; "
            "an uncontrolled page-cache-cold state is refused"
        )
    result.update(
        {
            "started_utc": started,
            "finished_utc": utc_now(),
            "wall_seconds": (time.perf_counter_ns() - monotonic_start) / 1.0e9,
        }
    )
    return result


def _wait_until_pidfd_ready(pidfd: int, timeout_seconds: float) -> bool:
    deadline = time.monotonic() + timeout_seconds
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0.0:
            return False
        try:
            readable, _, _ = select.select([pidfd], [], [], remaining)
        except InterruptedError:
            continue
        return bool(readable)


def _wait4_with_timeout(
    pid: int, timeout_seconds: float
) -> Tuple[int, int, Any, bool]:
    """Wait for exactly *pid* and return its high-resolution child rusage."""
    pidfd: Optional[int] = None
    if hasattr(os, "pidfd_open"):
        try:
            pidfd = os.pidfd_open(pid)
        except OSError:
            pidfd = None
    timed_out = False
    try:
        if pidfd is not None:
            ready = _wait_until_pidfd_ready(pidfd, timeout_seconds)
        else:
            deadline = time.monotonic() + timeout_seconds
            ready = False
            while time.monotonic() < deadline:
                waited_pid, status, usage = os.wait4(pid, os.WNOHANG)
                if waited_pid == pid:
                    return waited_pid, status, usage, False
                time.sleep(min(0.0005, max(0.0, deadline - time.monotonic())))
        if not ready:
            timed_out = True
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            if pidfd is not None:
                ready = _wait_until_pidfd_ready(pidfd, 5.0)
            else:
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline:
                    waited_pid, status, usage = os.wait4(pid, os.WNOHANG)
                    if waited_pid == pid:
                        return waited_pid, status, usage, True
                    time.sleep(0.0005)
                ready = False
            if not ready:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        waited_pid, status, usage = os.wait4(pid, 0)
        return waited_pid, status, usage, timed_out
    finally:
        if pidfd is not None:
            os.close(pidfd)


def _resource_record(
    command: Sequence[str],
    wall_seconds: float,
    wait_status: int,
    usage: Any,
    timed_out: bool,
) -> Dict[str, Any]:
    returncode = os.waitstatus_to_exitcode(wait_status)
    user_seconds = float(usage.ru_utime)
    system_seconds = float(usage.ru_stime)
    cpu_seconds = user_seconds + system_seconds
    return {
        "schema_version": SCHEMA_VERSION,
        "measurement_method": "os.wait4",
        "measurement_boundary": "python-launch -> taskset-exec -> target-tree -> wait4-reap",
        "wall_clock": "time.perf_counter_ns",
        "command": list(command),
        "wait_status": wait_status,
        "returncode": returncode,
        "timed_out": timed_out,
        "wall_seconds": wall_seconds,
        "user_seconds": user_seconds,
        "system_seconds": system_seconds,
        "cpu_seconds": cpu_seconds,
        "cpu_percent": 100.0 * cpu_seconds / wall_seconds if wall_seconds > 0.0 else 0.0,
        "max_rss_kbytes": int(usage.ru_maxrss),
        "major_faults": int(usage.ru_majflt),
        "minor_faults": int(usage.ru_minflt),
        "voluntary_context_switches": int(usage.ru_nvcsw),
        "involuntary_context_switches": int(usage.ru_nivcsw),
        "filesystem_inputs": int(usage.ru_inblock),
        "filesystem_outputs": int(usage.ru_oublock),
    }


def _launch_and_measure(
    command: Sequence[str],
    output: Any,
    environment: Mapping[str, str],
    timeout_seconds: float,
) -> Dict[str, Any]:
    monotonic_start = time.perf_counter_ns()
    try:
        process = subprocess.Popen(
            list(command),
            stdout=output,
            stderr=subprocess.STDOUT,
            env=dict(environment),
            start_new_session=True,
        )
    except OSError as error:
        raise HarnessError(f"failed to start measured process: {error}") from error
    waited_pid, wait_status, usage, timed_out = _wait4_with_timeout(
        process.pid, timeout_seconds
    )
    if waited_pid != process.pid:
        raise HarnessError(
            f"wait4 reaped unexpected pid {waited_pid}; expected {process.pid}"
        )
    process.returncode = os.waitstatus_to_exitcode(wait_status)
    wall_seconds = (time.perf_counter_ns() - monotonic_start) / 1.0e9
    return _resource_record(command, wall_seconds, wait_status, usage, timed_out)


def measure_launcher_overhead(
    taskset_binary: Path,
    cpu_list: str,
    probe_binary: Path,
    repetitions: int = 31,
) -> Dict[str, Any]:
    """Measure, but never subtract, the taskset/launch/reap null boundary."""
    command = [str(taskset_binary), "-c", cpu_list, str(probe_binary)]
    environment = benchmark_environment(1)
    samples = [
        _launch_and_measure(command, subprocess.DEVNULL, environment, 5.0)
        for _ in range(repetitions)
    ]
    failures = [
        sample
        for sample in samples
        if sample["timed_out"] or sample["returncode"] != 0
    ]
    if failures:
        raise HarnessError(f"launcher overhead calibration failed: {failures[0]}")
    walls = [float(sample["wall_seconds"]) for sample in samples]
    cpus = [float(sample["cpu_seconds"]) for sample in samples]
    return {
        "method": "taskset-null-wait4",
        "interpretation": "fixed-boundary estimate only; no value is subtracted from solve measurements",
        "command": command,
        "repetitions": repetitions,
        "wall_seconds": {
            "minimum": min(walls),
            "median": float(statistics.median(walls)),
            "maximum": max(walls),
        },
        "cpu_seconds": {
            "minimum": min(cpus),
            "median": float(statistics.median(cpus)),
            "maximum": max(cpus),
        },
        "samples": samples,
    }


def execute_timed(
    solve_command: Sequence[str],
    run_dir: Path,
    workers: int,
    cpu_list: str,
    timeout_seconds: float,
    taskset_binary: Path,
) -> Dict[str, Any]:
    log_path = run_dir / "solve.log"
    resource_path = run_dir / "resource-usage.json"
    measured_command = [str(taskset_binary), "-c", cpu_list, *solve_command]
    environment = benchmark_environment(workers)
    started_utc = utc_now()
    with log_path.open("wb") as log_stream:
        resource_record = _launch_and_measure(
            measured_command, log_stream, environment, timeout_seconds
        )
    write_json(resource_path, resource_record)
    return {
        "started_utc": started_utc,
        "finished_utc": utc_now(),
        "harness_wall_seconds": resource_record["wall_seconds"],
        "returncode": resource_record["returncode"],
        "timed_out": resource_record["timed_out"],
        "measured_command": measured_command,
        "solve_command": list(solve_command),
        "solve_command_shell": render_command(solve_command),
        "resource_measurement": {
            "method": resource_record["measurement_method"],
            "boundary": resource_record["measurement_boundary"],
            "raw_evidence": resource_path.name,
        },
        "cpu_list": cpu_list,
        "environment": {
            key: environment[key]
            for key in sorted(set(PERFORMANCE_ENV) | {"ASTROMETRY_INDEX_SHARD_WORKERS", "LANG", "LC_ALL", "TZ"})
        },
    }


def parse_resource_usage(path: Path) -> Dict[str, Any]:
    record = read_json(path)
    if record.get("measurement_method") != "os.wait4":
        raise HarnessError(f"{path}: resource measurement is not os.wait4")
    required_numeric = (
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
        "returncode",
        "wait_status",
    )
    missing = [
        name for name in required_numeric if not isinstance(record.get(name), (int, float))
    ]
    if missing:
        raise HarnessError(f"{path}: resource measurement lacks numeric fields {missing}")
    if not isinstance(record.get("timed_out"), bool):
        raise HarnessError(f"{path}: resource measurement lacks boolean timed_out")
    if record["wall_seconds"] <= 0.0:
        raise HarnessError(f"{path}: resource measurement has non-positive wall time")
    return {
        **{name: record[name] for name in required_numeric},
        "timed_out": record["timed_out"],
    }


def parameter_blocks(log_text: str) -> List[Dict[str, Any]]:
    lines = log_text.splitlines()
    blocks: List[Dict[str, Any]] = []
    position = 0
    while position < len(lines):
        if lines[position].strip() != "solver run parameters:":
            position += 1
            continue
        block: Dict[str, Any] = {"indexes": []}
        position += 1
        if position >= len(lines) or lines[position].strip() != "indexes:":
            block["malformed"] = "missing indexes header"
            blocks.append(block)
            continue
        position += 1
        while position < len(lines) and lines[position].startswith("  "):
            block["indexes"].append(lines[position].strip())
            position += 1
        while position < len(lines):
            stripped = lines[position].strip()
            if stripped == "solver run parameters:":
                break
            if stripped.startswith("startdepth "):
                block["startobj"] = int(stripped.split()[1])
            elif stripped.startswith("enddepth "):
                block["endobj"] = int(stripped.split()[1])
            elif stripped.startswith("fields "):
                block["fields"] = [int(value) for value in stripped.split()[1:]]
            elif stripped.startswith("[index-shard] pthread-pool submit"):
                break
            position += 1
        blocks.append(block)
    return blocks


def assert_log_shape(
    log_path: Path,
    expected_indexes: Sequence[Path],
    workers: int,
    expected_startobj: int,
    expected_endobj: int,
) -> Dict[str, Any]:
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    blocks = parameter_blocks(log_text)
    if len(blocks) != 1:
        raise HarnessError(f"{log_path}: expected one solver parameter block, found {len(blocks)}")
    block = blocks[0]
    try:
        logged_indexes = [str(Path(value).resolve(strict=True)) for value in block["indexes"]]
    except OSError as error:
        raise HarnessError(f"{log_path}: a logged index no longer resolves: {error}") from error
    expected = [str(path) for path in expected_indexes]
    if logged_indexes != expected:
        raise HarnessError(
            f"{log_path}: logged candidate list differs from immutable prefix; "
            f"expected {expected}, got {logged_indexes}"
        )
    if block.get("startobj") != expected_startobj or block.get("endobj") != expected_endobj:
        raise HarnessError(
            f"{log_path}: expected pass startobj={expected_startobj} endobj={expected_endobj}, "
            f"got startobj={block.get('startobj')} endobj={block.get('endobj')}"
        )
    if len(block.get("fields", [])) != 1:
        raise HarnessError(
            f"{log_path}: expected exactly one input field, got {block.get('fields')}"
        )
    submissions = [match.groupdict() for match in SUBMIT_RE.finditer(log_text)]
    lane_publications = [
        {
            "state": match.group("state"),
            "owner_worker": int(match.group("owner_worker")),
            "index_order": int(match.group("index_order")),
        }
        for match in ASSIST_LANE_PUBLISH_RE.finditer(log_text)
    ]
    lane_loans = [
        {
            "state": match.group("state"),
            "helper_worker": int(match.group("helper_worker")),
            "helper_slot": int(match.group("helper_slot")),
            "index_order": int(match.group("index_order")),
            "rc": (
                int(match.group("rc"))
                if match.group("rc") is not None
                else None
            ),
        }
        for match in ASSIST_LANE_LOAN_RE.finditer(log_text)
    ]
    assist_passes = [
        {key: int(value) for key, value in match.groupdict().items()}
        for match in ASSIST_PASS_RE.finditer(log_text)
    ]
    phase_records = [
        {
            match.group("key"): parse_profile_value(match.group("value"))
            for match in KEY_VALUE_RE.finditer(phase_match.group("payload"))
        }
        for phase_match in AB_PHASE_RE.finditer(log_text)
    ]
    geometry_fallbacks = [
        match.groupdict() for match in GEOMETRY_FALLBACK_RE.finditer(log_text)
    ]
    geometry_shared = len(GEOMETRY_SHARED_RE.findall(log_text))
    use_pthread = workers > 1 and len(expected_indexes) > 0
    geometry_requested = (
        expected_startobj > 0
        and (use_pthread or len(expected_indexes) > 1)
    )
    allowed_phase_modes = {
        "native",
        "empty",
        "flattened-owner",
        "assisted",
    }
    phase_mode_counts = {
        mode: sum(record.get("mode") == mode for record in phase_records)
        for mode in sorted(allowed_phase_modes)
    }
    for record in phase_records:
        if record.get("mode") not in allowed_phase_modes:
            raise HarnessError(
                f"{log_path}: unknown solver AB phase mode "
                f"{record.get('mode')!r}"
            )
        if record.get("phase") not in {"diagonal", "off-diagonal"}:
            raise HarnessError(
                f"{log_path}: unknown solver AB phase kind "
                f"{record.get('phase')!r}"
            )
        if record.get("resource") not in {
            "process-overlap",
            "unavailable",
        }:
            raise HarnessError(
                f"{log_path}: solver AB phase resource scope must be "
                f"'process-overlap' or 'unavailable', got "
                f"{record.get('resource')!r}"
            )
    if len(geometry_fallbacks) > 1 or geometry_shared > 1:
        raise HarnessError(
            f"{log_path}: expected at most one geometry outcome, "
            f"found shared={geometry_shared}, refused={geometry_fallbacks}"
        )
    if geometry_fallbacks and geometry_shared:
        raise HarnessError(
            f"{log_path}: geometry was both prepared and refused"
        )
    if geometry_requested and not (geometry_fallbacks or geometry_shared):
        raise HarnessError(
            f"{log_path}: pass omitted its required prefix-geometry outcome"
        )
    if not geometry_requested and (geometry_fallbacks or geometry_shared):
        raise HarnessError(
            f"{log_path}: pass unexpectedly prepared prefix geometry"
        )
    if not use_pthread:
        if (
            submissions
            or lane_publications
            or lane_loans
            or assist_passes
        ):
            raise HarnessError(
                f"{log_path}: serial pass unexpectedly submitted pthread work"
            )
        unexpected_modes = {
            mode: count
            for mode, count in phase_mode_counts.items()
            if mode != "native" and count
        }
        if unexpected_modes:
            raise HarnessError(
                f"{log_path}: serial pass reported executor-only phase modes "
                f"{unexpected_modes}"
            )
    else:
        if len(submissions) != 1:
            raise HarnessError(
                f"{log_path}: expected one pthread pass submission, found {len(submissions)}"
            )
        submission = {
            key: (
                value
                if key in {"inner_scheduler", "mmap_advice"}
                else int(value)
            )
            for key, value in submissions[0].items()
        }
        expected_submission_subset = {
            "workers": workers,
            "pool_workers": workers,
            "candidates": len(expected_indexes),
            "startobj": expected_startobj,
            "endobj": expected_endobj,
            "inner_scheduler": "dynamic-pool-lending",
            "mmap_advice": "normal",
        }
        mismatched_submission = {
            key: (expected, submission.get(key))
            for key, expected in expected_submission_subset.items()
            if submission.get(key) != expected
        }
        if mismatched_submission:
            raise HarnessError(
                f"{log_path}: pthread submission differs: "
                f"{mismatched_submission}"
            )
        if len(assist_passes) != 1:
            raise HarnessError(
                f"{log_path}: expected one dynamic assist-pass record, "
                f"found {len(assist_passes)}"
            )
        assist_pass = assist_passes[0]
        if assist_pass["waiters"] != 0:
            raise HarnessError(
                f"{log_path}: assist pass ended with "
                f"{assist_pass['waiters']} waiter(s)"
            )

        publishes = [
            (record["owner_worker"], record["index_order"])
            for record in lane_publications
            if record["state"] == "publish"
        ]
        unpublishes = [
            (record["owner_worker"], record["index_order"])
            for record in lane_publications
            if record["state"] == "unpublish"
        ]
        if sorted(publishes) != sorted(unpublishes):
            raise HarnessError(
                f"{log_path}: dynamic lane publish/unpublish lifecycle differs; "
                f"publish={publishes}, unpublish={unpublishes}"
            )
        if any(
            owner < 0
            or owner >= workers
            or index_order < 0
            or index_order >= len(expected_indexes)
            for owner, index_order in publishes
        ):
            raise HarnessError(
                f"{log_path}: dynamic lane publication is out of bounds"
            )

        joins = [
            (
                record["helper_worker"],
                record["helper_slot"],
                record["index_order"],
            )
            for record in lane_loans
            if record["state"] == "join"
        ]
        leaves = [
            (
                record["helper_worker"],
                record["helper_slot"],
                record["index_order"],
            )
            for record in lane_loans
            if record["state"] == "leave"
        ]
        if sorted(joins) != sorted(leaves):
            raise HarnessError(
                f"{log_path}: dynamic lane join/leave lifecycle differs; "
                f"join={joins}, leave={leaves}"
            )
        if any(
            record["rc"] is None
            for record in lane_loans
            if record["state"] == "leave"
        ):
            raise HarnessError(
                f"{log_path}: dynamic lane leave omitted its result"
            )
        failed_leaves = [
            record
            for record in lane_loans
            if record["state"] == "leave" and record["rc"] != 0
        ]
        if failed_leaves:
            raise HarnessError(
                f"{log_path}: dynamic lane helper failed: {failed_leaves}"
            )
        if any(
            helper < 0
            or helper >= workers
            or slot <= 0
            or slot >= workers
            or index_order < 0
            or index_order >= len(expected_indexes)
            for helper, slot, index_order in joins
        ):
            raise HarnessError(
                f"{log_path}: dynamic lane loan is out of bounds"
            )
        if assist_pass["loans"] != len(joins):
            raise HarnessError(
                f"{log_path}: assist-pass loans={assist_pass['loans']} "
                f"but observed joins={len(joins)}"
            )
        if phase_records and phase_mode_counts["assisted"] and not joins:
            raise HarnessError(
                f"{log_path}: assisted phase mode has no observed helper loan"
            )
        if phase_records and joins and not phase_mode_counts["assisted"]:
            raise HarnessError(
                f"{log_path}: observed helper loans without an assisted phase"
            )

    assistance_observation = (
        "not-applicable-serial"
        if not use_pthread
        else "assisted"
        if phase_mode_counts["assisted"]
        else "loan-observed"
        if lane_loans
        else "flattened-owner"
        if phase_mode_counts["flattened-owner"]
        else "native-or-empty"
    )
    return {
        "ok": True,
        "candidate_count": len(logged_indexes),
        "candidate_paths": logged_indexes,
        "startobj": block["startobj"],
        "endobj": block["endobj"],
        "field": block["fields"][0],
        "pthread_submissions": submissions,
        "assist_lane_publications": lane_publications,
        "assist_lane_loans": lane_loans,
        "assist_pass": assist_passes[0] if assist_passes else None,
        "phase_mode_counts": phase_mode_counts,
        "phase_resource_scope": (
            "process-overlap"
            if any(
                record.get("resource") == "process-overlap"
                for record in phase_records
            )
            else "unavailable"
            if phase_records
            else "not-observed"
        ),
        "assistance_observation": assistance_observation,
        "geometry_fallbacks": geometry_fallbacks,
        "geometry_shared": geometry_shared,
        "geometry_classification": (
            "not-requested"
            if not geometry_requested
            else "shared-readonly"
            if geometry_shared
            else f"refused-{geometry_fallbacks[0]['reason']}"
        ),
        "scheduler_classification": (
            "serial"
            if not use_pthread
            else "dynamic-pool-lending"
        ),
    }


def assert_constructed_command(
    command_path: Path,
    execution: Mapping[str, Any],
    expected_indexes: Sequence[Path],
    expected_field_objects: str,
    expected_startobj: int,
    expected_endobj: int,
    expected_input: Path,
    measurement_mode: str,
) -> Dict[str, Any]:
    command_record = read_json(command_path)
    command = command_record.get("argv")
    executed_command = execution.get("solve_command")
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        raise HarnessError(f"{command_path}: argv must be a string array")
    if command != executed_command:
        raise HarnessError(f"{command_path}: persisted and executed solve argv differ")
    expected_verbose = {"timing": 0, "phase": 1, "detailed": 2}[measurement_mode]
    if command.count("--verbose") != expected_verbose:
        raise HarnessError(
            f"{command_path}: expected {expected_verbose} --verbose arguments, "
            f"found {command.count('--verbose')}"
        )
    depth_positions = [position for position, value in enumerate(command) if value == "--depth"]
    if len(depth_positions) != 1 or depth_positions[0] + 1 >= len(command):
        raise HarnessError(f"{command_path}: expected exactly one complete --depth option")
    if command[depth_positions[0] + 1] != expected_field_objects:
        raise HarnessError(
            f"{command_path}: expected --depth field-object range "
            f"{expected_field_objects!r}, "
            f"got {command[depth_positions[0] + 1]!r}"
        )
    index_values: List[str] = []
    for position, value in enumerate(command):
        if value != "--index-file":
            continue
        if position + 1 >= len(command):
            raise HarnessError(f"{command_path}: incomplete --index-file option")
        index_values.append(command[position + 1])
    expected_paths = [str(path) for path in expected_indexes]
    if index_values != expected_paths:
        raise HarnessError(
            f"{command_path}: constructed index prefix differs; "
            f"expected {expected_paths}, got {index_values}"
        )
    last_index_position = max(
        position for position, value in enumerate(command) if value == "--index-file"
    )
    common_tail = command[last_index_position + 2 : -1]
    invalid_common_tokens = [
        value for value in common_tail if not value.startswith("--") or value == "--"
    ]
    if invalid_common_tokens:
        raise HarnessError(
            f"{command_path}: common-argument tail contains positional tokens: "
            f"{invalid_common_tokens}"
        )
    if not command or command[-1] != str(expected_input) or command.count(str(expected_input)) != 1:
        raise HarnessError(f"{command_path}: expected exactly one final input argument")
    return {
        "ok": True,
        "validation": "constructed_command",
        "candidate_count": len(index_values),
        "candidate_paths": index_values,
        "field_objects": expected_field_objects,
        "startobj": expected_startobj,
        "endobj": expected_endobj,
        "input": str(expected_input),
        "verbose_count": expected_verbose,
    }


def parse_profile_value(value: str) -> Any:
    candidate = value.rstrip(",")
    if candidate.endswith("%"):
        candidate = candidate[:-1]
    work_with_percent = re.fullmatch(
        r"(?P<work>[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
        r"(?:[eE][-+]?[0-9]+)?)"
        r"\([-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)%\)",
        candidate,
    )
    if work_with_percent:
        return float(work_with_percent.group("work"))
    try:
        if re.fullmatch(r"[-+]?[0-9]+", candidate):
            return int(candidate)
        if re.fullmatch(
            r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?",
            candidate,
        ):
            return float(candidate)
    except ValueError:
        pass
    return value


def parse_profiles(log_path: Path) -> Dict[str, List[Dict[str, Any]]]:
    """Extract stable key/value profile records without assuming key order."""
    profiles: Dict[str, List[Dict[str, Any]]] = {}
    for raw_line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        profile_name: Optional[str] = None
        payload = ""
        for prefix, name in (
            ("[solve-field-profile]", "solve_field"),
            ("[engine-profile]", "engine"),
            ("[onefield-profile]", "onefield"),
            ("[onefield-field-profile]", "onefield_field"),
            ("[solver-ab-phase]", "solver_ab_phase"),
            ("[solver] phase-profile", "solver_phase"),
            ("[solver] verification-context", "solver_verification"),
            ("[index-shard] aux-pass", "index_shard_aux"),
            ("[index-shard] assist-pass", "index_shard_assist"),
            ("[index-shard] reducer-pass", "index_shard_reducer"),
            ("[index-shard] context-pass", "index_shard_context"),
            ("[index-shard] solver-pass", "index_shard_solver"),
            ("[index-shard] phase-profile", "index_shard_phase"),
            ("[index-shard] pass-detail", "index_shard_pass"),
            ("[index-shard] task-profile", "index_shard_task"),
            ("[index-shard] pass-resource", "index_shard_resource"),
        ):
            if line.startswith(prefix):
                profile_name = name
                payload = line[len(prefix) :].strip()
                break
        if profile_name is None:
            continue
        record = {
            match.group("key"): parse_profile_value(match.group("value"))
            for match in KEY_VALUE_RE.finditer(payload)
        }
        record["raw"] = line
        profiles.setdefault(profile_name, []).append(record)
    return profiles


def extract_winner(
    log_path: Path,
    solved: bool,
    expected_indexes: Sequence[Path],
    workers: int,
) -> Dict[str, Any]:
    lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    use_pthread = workers > 1 and len(expected_indexes) > 0
    if use_pthread:
        orders = [int(match.group("order")) for match in REDUCED_WINNER_RE.finditer("\n".join(lines))]
        if solved and len(orders) != 1:
            raise HarnessError(
                f"{log_path}: expected one authoritative reduced winner, found {orders}"
            )
        if not solved and orders:
            raise HarnessError(f"{log_path}: unsolved run logged reduced winner(s) {orders}")
        order = orders[0] if orders else None
    else:
        current_index: Optional[Path] = None
        order = None
        index_lookup = {path: position for position, path in enumerate(expected_indexes)}
        for line in lines:
            stripped = line.strip()
            match = SERIAL_TRY_RE.match(stripped) or SERIAL_VERIFY_RE.match(stripped)
            if match:
                try:
                    current_index = Path(match.group("path")).resolve(strict=True)
                except OSError as error:
                    raise HarnessError(
                        f"{log_path}: winning-index candidate no longer resolves: {error}"
                    ) from error
            if SOLVED_INDEX_RE.match(stripped):
                if current_index not in index_lookup:
                    raise HarnessError(
                        f"{log_path}: serial solution could not be mapped to the immutable prefix"
                    )
                order = index_lookup[current_index]
                break
        if solved and order is None:
            raise HarnessError(f"{log_path}: solved serial run has no identifiable winning index")
        if not solved and order is not None:
            raise HarnessError(f"{log_path}: unsolved serial run has an identified winning index")
    if order is not None and not 0 <= order < len(expected_indexes):
        raise HarnessError(
            f"{log_path}: winning index order {order} is outside {len(expected_indexes)} candidates"
        )
    return {
        "order_zero_based": order,
        "position_one_based": order + 1 if order is not None else None,
        "path": str(expected_indexes[order]) if order is not None else None,
    }


def run_wcsinfo(wcsinfo: Path, wcs_path: Path, destination: Path) -> str:
    try:
        completed = subprocess.run(
            [str(wcsinfo), str(wcs_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=benchmark_environment(1),
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise HarnessError(f"failed to run wcsinfo for {wcs_path}: {error}") from error
    destination.write_bytes(completed.stdout)
    if completed.returncode != 0:
        raise HarnessError(f"wcsinfo returned {completed.returncode} for {wcs_path}")
    return sha256_bytes(completed.stdout)


def output_inventory(run_dir: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for path in sorted(run_dir.iterdir(), key=lambda item: item.name):
        if path.is_file():
            stat = path.stat()
            records.append({"name": path.name, "size": stat.st_size})
    return records


def evidence_hashes(run_dir: Path) -> Dict[str, str]:
    names = ["command.json", "resource-usage.json", "solve.log", "wcsinfo.txt"]
    for optional_name in ("cache-control.json", "solve.solved", "solve.wcs"):
        if (run_dir / optional_name).is_file():
            names.append(optional_name)
    return {name: sha256_file(run_dir / name) for name in names}


def verify_evidence_hashes(run_dir: Path, expected: Any) -> None:
    if not isinstance(expected, dict) or not expected:
        raise HarnessError(f"run evidence hashes are missing: {run_dir}")
    for name, expected_digest in expected.items():
        if name not in EVIDENCE_FILE_NAMES:
            raise HarnessError(f"unexpected run evidence name {name!r}: {run_dir}")
        if not isinstance(expected_digest, str) or not re.fullmatch(
            r"[0-9a-f]{64}", expected_digest
        ):
            raise HarnessError(f"invalid evidence digest for {name!r}: {run_dir}")
        path = run_dir / name
        if not path.is_file():
            raise HarnessError(f"run evidence file is missing: {path}")
        actual_digest = sha256_file(path)
        if actual_digest != expected_digest:
            raise HarnessError(
                f"run evidence hash mismatch for {path}: {expected_digest} -> {actual_digest}"
            )


def required_evidence_names(cache_state: str, solved: bool) -> set[str]:
    required = {"command.json", "resource-usage.json", "solve.log", "wcsinfo.txt"}
    if cache_state == "cold":
        required.add("cache-control.json")
    if solved:
        required.update(("solve.solved", "solve.wcs"))
    return required


def verify_result_evidence(run_dir: Path, record: Mapping[str, Any]) -> None:
    resource_timing = parse_resource_usage(run_dir / "resource-usage.json")
    recorded_timing = record.get("time")
    if not isinstance(recorded_timing, dict):
        raise HarnessError(f"run has no timing object: {run_dir}")
    for name, value in resource_timing.items():
        if recorded_timing.get(name) != value:
            raise HarnessError(
                f"run timing differs from raw resource evidence at {name}: {run_dir}"
            )
    solved = record.get("solved")
    if not isinstance(solved, bool):
        raise HarnessError(f"run solved status is not boolean: {run_dir}")
    signature = record.get("correctness_signature")
    signature_method = record.get("correctness_signature_method")
    wcsinfo_sha256 = record.get("wcsinfo_sha256")
    raw_wcs_sha256 = record.get("raw_wcs_sha256")
    solved_path = run_dir / "solve.solved"
    wcs_path = run_dir / "solve.wcs"
    wcsinfo_path = run_dir / "wcsinfo.txt"
    if signature_method != SIGNATURE_METHOD:
        raise HarnessError(f"run has unsupported correctness signature method: {run_dir}")
    if not isinstance(wcsinfo_sha256, str) or not re.fullmatch(
        r"[0-9a-f]{64}", wcsinfo_sha256
    ):
        raise HarnessError(f"run has invalid wcsinfo digest: {run_dir}")
    if sha256_file(wcsinfo_path) != wcsinfo_sha256:
        raise HarnessError(f"wcsinfo digest mismatch: {run_dir}")
    if solved:
        if not solved_path.is_file() or solved_path.stat().st_size == 0:
            raise HarnessError(f"solved run has no non-empty solved marker: {run_dir}")
        if not wcs_path.is_file():
            raise HarnessError(f"solved run has no WCS: {run_dir}")
        if not isinstance(signature, str) or not re.fullmatch(r"[0-9a-f]{64}", signature):
            raise HarnessError(f"solved run has invalid correctness signature: {run_dir}")
        try:
            actual_signature = wcs_signature(wcs_path)
        except WcsSignatureError as error:
            raise HarnessError(f"cannot validate canonical WCS signature: {error}") from error
        if actual_signature != signature:
            raise HarnessError(f"canonical WCS signature mismatch: {run_dir}")
        if not isinstance(raw_wcs_sha256, str) or not re.fullmatch(
            r"[0-9a-f]{64}", raw_wcs_sha256
        ):
            raise HarnessError(f"solved run has invalid raw WCS digest: {run_dir}")
        if sha256_file(wcs_path) != raw_wcs_sha256:
            raise HarnessError(f"raw WCS digest mismatch: {run_dir}")
    else:
        if solved_path.exists() or wcs_path.exists():
            raise HarnessError(f"unsolved run retains solved artifacts: {run_dir}")
        if signature != "UNSOLVED" or raw_wcs_sha256 is not None:
            raise HarnessError(f"unsolved result metadata is inconsistent: {run_dir}")
        if wcsinfo_path.read_bytes() != b"UNSOLVED\n":
            raise HarnessError(f"unsolved wcsinfo sentinel is inconsistent: {run_dir}")


def finalize_run(
    run_dir: Path,
    execution: Mapping[str, Any],
    plan_item: Mapping[str, Any],
    cache_state: str,
    expected_indexes: Sequence[Path],
    expected_field_objects: str,
    expected_startobj: int,
    expected_endobj: int,
    expected_input: Path,
    wcsinfo: Path,
    measurement_mode: str,
    suite_fingerprint: str,
    workload_fingerprint: str,
) -> Dict[str, Any]:
    if execution["timed_out"]:
        raise HarnessError(f"timed solve exceeded timeout in {run_dir}")
    if execution["returncode"] != 0:
        raise HarnessError(f"timed solve returned {execution['returncode']} in {run_dir}")
    timing = parse_resource_usage(run_dir / "resource-usage.json")
    if timing["returncode"] != 0:
        raise HarnessError(f"wait4 recorded exit status {timing['returncode']} in {run_dir}")
    if not math.isclose(
        float(timing["wall_seconds"]),
        float(execution["harness_wall_seconds"]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise HarnessError(f"resource/effective wall-time mismatch in {run_dir}")
    assertions = assert_constructed_command(
        run_dir / "command.json",
        execution,
        expected_indexes,
        expected_field_objects,
        expected_startobj,
        expected_endobj,
        expected_input,
        measurement_mode,
    )
    if measurement_mode == "timing":
        assertions["observed_log_shape"] = "not_observable_quiet"
        profile_validation = "not_observable_quiet"
    else:
        assertions.update(
            assert_log_shape(
                run_dir / "solve.log",
                expected_indexes,
                int(plan_item["workers"]),
                expected_startobj,
                expected_endobj,
            )
        )
        assertions["observed_log_shape"] = "validated_verbose"
        profile_validation = "validated_verbose"
    solved_path = run_dir / "solve.solved"
    wcs_path = run_dir / "solve.wcs"
    solved_marker_exists = solved_path.is_file()
    if solved_marker_exists and solved_path.stat().st_size == 0:
        raise HarnessError(f"{run_dir}: solved marker is empty")
    solved = solved_marker_exists
    if solved != wcs_path.is_file():
        raise HarnessError(
            f"{run_dir}: solved marker/WCS disagreement (solved={solved}, wcs={wcs_path.is_file()})"
        )
    if solved:
        wcsinfo_sha256 = run_wcsinfo(wcsinfo, wcs_path, run_dir / "wcsinfo.txt")
        try:
            signature = wcs_signature(wcs_path)
        except WcsSignatureError as error:
            raise HarnessError(f"cannot create canonical WCS signature: {error}") from error
        raw_wcs_sha256: Optional[str] = sha256_file(wcs_path)
    else:
        (run_dir / "wcsinfo.txt").write_text("UNSOLVED\n", encoding="ascii")
        signature = "UNSOLVED"
        wcsinfo_sha256 = sha256_file(run_dir / "wcsinfo.txt")
        raw_wcs_sha256 = None
    if measurement_mode == "timing":
        winner = {
            "observation": "not_observable_quiet",
            "order_zero_based": None,
            "position_one_based": None,
            "path": None,
        }
    else:
        winner = extract_winner(
            run_dir / "solve.log", solved, expected_indexes, int(plan_item["workers"])
        )
    profiles = parse_profiles(run_dir / "solve.log")
    if measurement_mode == "timing" and profiles:
        raise HarnessError(
            f"{run_dir}: quiet timing run unexpectedly emitted structured profiles"
        )
    record: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "complete": True,
        "suite_fingerprint": suite_fingerprint,
        "workload_fingerprint": workload_fingerprint,
        **dict(plan_item),
        "cache_state": cache_state,
        "run_dir": str(run_dir),
        "execution": dict(execution),
        "time": timing,
        "assertions": assertions,
        "profile_validation": profile_validation,
        "profiles": profiles,
        "solved": solved,
        "correctness_signature": signature,
        "correctness_signature_method": SIGNATURE_METHOD,
        "wcsinfo_sha256": wcsinfo_sha256,
        "raw_wcs_sha256": raw_wcs_sha256,
        "winner": winner,
        "outputs": output_inventory(run_dir),
        "evidence_sha256": evidence_hashes(run_dir),
    }
    write_json(run_dir / "run.json", record)
    return record


def load_config(config_path: Path) -> Dict[str, Any]:
    config = read_json(config_path)
    if config.get("schema_version") != SCHEMA_VERSION:
        raise HarnessError(
            f"config schema_version must be {SCHEMA_VERSION}, got {config.get('schema_version')!r}"
        )
    return config


def ensure_output_root(output_root: Path, resume: bool) -> None:
    if output_root.exists():
        if not output_root.is_dir():
            raise HarnessError(f"output root exists and is not a directory: {output_root}")
        if not resume and any(output_root.iterdir()):
            raise HarnessError(f"output root is not empty (use --resume if appropriate): {output_root}")
    else:
        output_root.mkdir(parents=True)


def validate_affinity(taskset_binary: Path, cpu_list: str) -> int:
    result = run_capture(
        [
            str(taskset_binary),
            "-c",
            cpu_list,
            sys.executable,
            "-c",
            "import os; print(len(os.sched_getaffinity(0)))",
        ]
    )
    if result.get("returncode") != 0:
        raise HarnessError(
            f"CPU affinity {cpu_list!r} is not usable: {result.get('output', result.get('error'))}"
        )
    try:
        available = int(str(result.get("output", "")).strip())
    except ValueError as error:
        raise HarnessError(f"could not determine CPU count for affinity {cpu_list!r}") from error
    if available < max(WORKER_COUNTS):
        raise HarnessError(
            f"CPU affinity {cpu_list!r} exposes {available} CPUs; W4 requires at least 4"
        )
    return available


def pair_complete(
    pair_dir: Path,
    expected_item: Mapping[str, Any],
    suite_fingerprint: str,
    workload_fingerprint: str,
) -> bool:
    marker = pair_dir / "pair.json"
    if not marker.exists():
        return False
    record = read_json(marker)
    if not record.get("complete"):
        return False
    if record.get("suite_fingerprint") != suite_fingerprint:
        raise HarnessError(f"resume marker {marker} belongs to a different suite fingerprint")
    if record.get("workload_fingerprint") != workload_fingerprint:
        raise HarnessError(f"resume marker {marker} belongs to a different workload fingerprint")
    for key in (
        "sequence",
        "repetition",
        "order_in_repetition",
        "index_count",
        "range_id",
        "workers",
    ):
        if record.get(key) != expected_item.get(key):
            raise HarnessError(f"resume marker {marker} does not match deterministic plan at {key}")
    for state in CACHE_STATES:
        run_record = read_json(pair_dir / state / "run.json")
        if not run_record.get("complete"):
            raise HarnessError(f"resume run is incomplete: {pair_dir / state}")
        if run_record.get("suite_fingerprint") != suite_fingerprint:
            raise HarnessError(f"resume run belongs to a different suite: {pair_dir / state}")
        if run_record.get("workload_fingerprint") != workload_fingerprint:
            raise HarnessError(
                f"resume run belongs to a different workload: {pair_dir / state}"
            )
        if run_record.get("cache_state") != state:
            raise HarnessError(f"resume run has wrong cache state: {pair_dir / state}")
        if not run_record.get("assertions", {}).get("ok"):
            raise HarnessError(f"resume run has failed assertions: {pair_dir / state}")
        evidence_names = set(run_record.get("evidence_sha256", {}))
        missing_evidence = required_evidence_names(
            state, bool(run_record.get("solved"))
        ) - evidence_names
        if missing_evidence:
            raise HarnessError(
                f"resume run is missing evidence hashes {sorted(missing_evidence)}: "
                f"{pair_dir / state}"
            )
        verify_evidence_hashes(pair_dir / state, run_record.get("evidence_sha256"))
        verify_result_evidence(pair_dir / state, run_record)
        for key in (
            "sequence",
            "repetition",
            "order_in_repetition",
            "index_count",
            "range_id",
            "workers",
        ):
            if run_record.get(key) != expected_item.get(key):
                raise HarnessError(
                    f"resume run {pair_dir / state} does not match deterministic plan at {key}"
                )
    return True


def assert_no_partial_pair(pair_dir: Path) -> None:
    if pair_dir.exists() and any(pair_dir.iterdir()):
        raise HarnessError(
            f"incomplete pair directory exists: {pair_dir}; preserve it for diagnosis and use a new output root"
        )


def verify_file_identities(paths: Sequence[Path], original: Sequence[Mapping[str, Any]]) -> None:
    for path, record in zip(paths, original):
        current = stat_identity(path)
        expected = {key: record[key] for key in current}
        if current != expected:
            raise HarnessError(f"input file changed during benchmark: {path}: {expected} -> {current}")


def write_plan(path: Path, plan: Sequence[Mapping[str, Any]]) -> None:
    lines = [
        "sequence\trepetition\torder_in_repetition\tindex_count\trange_id\tworkers"
    ]
    for item in plan:
        lines.append(
            "\t".join(
                str(item[key])
                for key in (
                    "sequence",
                    "repetition",
                    "order_in_repetition",
                    "index_count",
                    "range_id",
                    "workers",
                )
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="override config output_root")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)

    config_path = args.config.expanduser().resolve(strict=True)
    config = load_config(config_path)
    solve_field = require_regular_file(config.get("solve_field"), "solve_field", executable=True)
    wcsinfo = require_regular_file(config.get("wcsinfo"), "wcsinfo", executable=True)
    input_path = require_regular_file(config.get("input"), "input")
    manifest_path = require_regular_file(config.get("index_manifest"), "index_manifest")
    taskset_binary = require_regular_file(config.get("taskset_binary", "/usr/bin/taskset"), "taskset_binary", executable=True)
    launcher_probe_binary = require_regular_file(
        config.get("launcher_probe_binary", "/usr/bin/true"),
        "launcher_probe_binary",
        executable=True,
    )
    indexes = load_indexes(manifest_path)
    pipeline_files = load_pipeline_files(config.get("pipeline_files"), solve_field)

    output_value = args.output if args.output is not None else config.get("output_root")
    if output_value is None:
        raise HarnessError("set output_root in config or pass --output")
    output_root = Path(output_value).expanduser().resolve()
    ensure_output_root(output_root, args.resume)

    repetitions = config.get("repetitions", 7)
    if not isinstance(repetitions, int) or repetitions < 7:
        raise HarnessError("repetitions must be an integer >= 7")
    seed = config.get("random_seed", 20260721)
    rotation_stride = config.get("rotation_stride", 5)
    if not isinstance(seed, int) or not isinstance(rotation_stride, int):
        raise HarnessError("random_seed and rotation_stride must be integers")
    cpu_list = config.get("cpu_list")
    if not isinstance(cpu_list, str) or not cpu_list.strip():
        raise HarnessError("cpu_list must be a non-empty taskset CPU-list string")
    timeout_seconds = config.get("timeout_seconds", 300)
    if not isinstance(timeout_seconds, (int, float)) or timeout_seconds <= 0:
        raise HarnessError("timeout_seconds must be positive")
    cache_config = config.get("cache_control")
    if not isinstance(cache_config, dict):
        raise HarnessError("cache_control must be an object")
    common_args = validate_common_args(config.get("common_args"))
    measurement_mode = config.get("measurement_mode", "timing")
    if measurement_mode not in ("timing", "phase", "detailed"):
        raise HarnessError("measurement_mode must be 'timing', 'phase', or 'detailed'")

    field_object_ranges_value = config.get("field_object_ranges")
    parsed_ranges = parse_field_object_ranges(field_object_ranges_value)
    field_object_ranges: Dict[str, Tuple[str, int, int]] = {
        range_id: (objects, startobj, endobj)
        for range_id, objects, startobj, endobj in parsed_ranges
    }
    range_ids = [range_id for range_id, _, _, _ in parsed_ranges]

    affinity_cpu_count = validate_affinity(taskset_binary, cpu_list)
    launcher_calibration = measure_launcher_overhead(
        taskset_binary, cpu_list, launcher_probe_binary
    )
    plan = condition_plan(range_ids, repetitions, seed, rotation_stride)
    provenance = collect_provenance(
        config_path,
        config,
        solve_field,
        wcsinfo,
        input_path,
        manifest_path,
        indexes,
        pipeline_files,
        taskset_binary,
        launcher_probe_binary,
    )
    provenance["launcher_calibration"] = launcher_calibration
    workload_fingerprint = compute_workload_fingerprint(config, provenance, plan)
    provenance["workload_config"] = workload_config(config)
    provenance["workload_fingerprint"] = workload_fingerprint
    fingerprint_value = {
        "config": config,
        "solve_field": provenance["solve_field"],
        "wcsinfo": provenance["wcsinfo"],
        "input": provenance["input"],
        "manifest_sha256": provenance["manifest_sha256"],
        "config_file": provenance["config_file"],
        "manifest_file": provenance["manifest_file"],
        "indexes": provenance["indexes"],
        "pipeline_files": provenance["pipeline_files"],
        "taskset_binary": provenance["taskset_binary"],
        "launcher_probe_binary": provenance["launcher_probe_binary"],
        "harness": provenance["harness"],
        "plan": plan,
    }
    fingerprint = sha256_bytes(canonical_json(fingerprint_value))
    suite_path = output_root / "suite.json"
    if args.resume:
        if not suite_path.exists():
            raise HarnessError(f"--resume requested but {suite_path} does not exist")
        old_suite = read_json(suite_path)
        if old_suite.get("fingerprint") != fingerprint:
            raise HarnessError("resume fingerprint differs from current binaries, inputs, config, or plan")
        if old_suite.get("workload_fingerprint") != workload_fingerprint:
            raise HarnessError("resume workload fingerprint differs from the current experiment")
    else:
        (output_root / "provenance").mkdir(parents=True, exist_ok=True)
        write_json(output_root / "provenance" / "provenance.json", provenance)
        (output_root / "provenance" / "config.json").write_bytes(canonical_json(config))
        (output_root / "provenance" / "indexes.manifest.txt").write_text(
            "".join(f"{path}\n" for path in indexes), encoding="utf-8"
        )
        write_plan(output_root / "plan.tsv", plan)
        write_json(
            suite_path,
            {
                "schema_version": SCHEMA_VERSION,
                "complete": False,
                "created_utc": utc_now(),
                "fingerprint": fingerprint,
                "workload_fingerprint": workload_fingerprint,
                "pair_count": len(plan),
                "cell_count": len(INDEX_COUNTS)
                * len(range_ids)
                * len(WORKER_COUNTS)
                * len(CACHE_STATES),
                "run_count": len(plan) * len(CACHE_STATES),
                "repetitions": repetitions,
                "index_counts": list(INDEX_COUNTS),
                "workers": list(WORKER_COUNTS),
                "cache_states": list(CACHE_STATES),
                "field_object_ranges": field_object_ranges_value,
                "range_ids": range_ids,
                "measurement_mode": measurement_mode,
                "cpu_list": cpu_list,
                "affinity_cpu_count": affinity_cpu_count,
                "resource_measurement": "os.wait4",
                "launcher_overhead_wall_median_seconds": launcher_calibration[
                    "wall_seconds"
                ]["median"],
            },
        )

    runner_path = Path(__file__).resolve()
    summarizer_path = runner_path.with_name("summarize.py")
    original_files = [
        provenance["config_file"],
        provenance["manifest_file"],
        provenance["solve_field"],
        provenance["wcsinfo"],
        provenance["input"],
        *provenance["indexes"],
        *provenance["pipeline_files"],
        provenance["taskset_binary"],
        provenance["launcher_probe_binary"],
        provenance["harness"]["run_matrix.py"],
        provenance["harness"]["summarize.py"],
        provenance["harness"]["fits_wcs_signature.py"],
    ]
    immutable_paths = [
        config_path,
        manifest_path,
        solve_field,
        wcsinfo,
        input_path,
        *indexes,
        *pipeline_files,
        taskset_binary,
        launcher_probe_binary,
        runner_path,
        summarizer_path,
        runner_path.with_name("fits_wcs_signature.py"),
    ]
    for item in plan:
        pair_dir = output_root / pair_relative_path(item)
        if args.resume and pair_complete(
            pair_dir, item, fingerprint, workload_fingerprint
        ):
            print(f"skip complete pair {item['sequence']}/{len(plan)}: {pair_dir}", flush=True)
            continue
        assert_no_partial_pair(pair_dir)
        cold_dir = pair_dir / "cold"
        warm_dir = pair_dir / "warm"
        cold_dir.mkdir(parents=True)
        warm_dir.mkdir(parents=True)
        selected_indexes = indexes[: int(item["index_count"])]
        (
            field_objects,
            expected_startobj,
            expected_endobj,
        ) = field_object_ranges[str(item["range_id"])]
        cold_command = build_solve_command(
            solve_field,
            cold_dir,
            input_path,
            selected_indexes,
            field_objects,
            common_args,
            measurement_mode,
        )
        warm_command = build_solve_command(
            solve_field,
            warm_dir,
            input_path,
            selected_indexes,
            field_objects,
            common_args,
            measurement_mode,
        )
        write_json(cold_dir / "command.json", {"argv": cold_command})
        write_json(warm_dir / "command.json", {"argv": warm_command})
        print(
            f"run pair {item['sequence']}/{len(plan)} r{item['repetition']} "
            f"i{item['index_count']} range={item['range_id']} W{item['workers']}",
            flush=True,
        )
        cache_evidence = run_drop_caches(cache_config)
        cold_execution = execute_timed(
            cold_command,
            cold_dir,
            int(item["workers"]),
            cpu_list,
            float(timeout_seconds),
            taskset_binary,
        )
        if cold_execution["timed_out"] or cold_execution["returncode"] != 0:
            write_json(cold_dir / "failed-execution.json", cold_execution)
            raise HarnessError(f"cold execution failed in {cold_dir}")
        warm_execution = execute_timed(
            warm_command,
            warm_dir,
            int(item["workers"]),
            cpu_list,
            float(timeout_seconds),
            taskset_binary,
        )
        write_json(cold_dir / "cache-control.json", cache_evidence)
        if warm_execution["timed_out"] or warm_execution["returncode"] != 0:
            write_json(warm_dir / "failed-execution.json", warm_execution)
            raise HarnessError(f"warm execution failed in {warm_dir}")
        finalize_run(
            cold_dir,
            cold_execution,
            item,
            "cold",
            selected_indexes,
            field_objects,
            expected_startobj,
            expected_endobj,
            input_path,
            wcsinfo,
            measurement_mode,
            fingerprint,
            workload_fingerprint,
        )
        finalize_run(
            warm_dir,
            warm_execution,
            item,
            "warm",
            selected_indexes,
            field_objects,
            expected_startobj,
            expected_endobj,
            input_path,
            wcsinfo,
            measurement_mode,
            fingerprint,
            workload_fingerprint,
        )
        write_json(
            pair_dir / "pair.json",
            {
                "complete": True,
                "suite_fingerprint": fingerprint,
                "workload_fingerprint": workload_fingerprint,
                **dict(item),
            },
        )

    verify_file_identities(immutable_paths, original_files)
    completed_suite = read_json(suite_path)
    completed_suite.update({"complete": True, "completed_utc": utc_now()})
    write_json(suite_path, completed_suite)
    result = subprocess.run([sys.executable, str(summarizer_path), str(output_root)], check=False)
    if result.returncode != 0:
        print(
            f"matrix completed, but correctness/statistical validation failed (summary exit {result.returncode})",
            file=sys.stderr,
        )
        return result.returncode
    print(f"matrix and summary complete: {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except HarnessError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
