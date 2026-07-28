#!/usr/bin/env python3
"""Unit and synthetic integration tests for the narrow-index harness."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run_matrix  # noqa: E402
import summarize  # noqa: E402
from fits_wcs_signature import SIGNATURE_METHOD, wcs_signature  # noqa: E402


def fits_card(keyword: str, value: object = None, comment: str = "") -> str:
    if value is None:
        text = f"{keyword:<8}"
        if comment:
            text += comment
    else:
        if isinstance(value, str):
            rendered = "'" + value.replace("'", "''") + "'"
        elif isinstance(value, bool):
            rendered = "T" if value else "F"
        else:
            rendered = str(value)
        text = f"{keyword:<8}= {rendered:>20}"
        if comment:
            text += f" / {comment}"
    return text[:80].ljust(80)


def synthetic_wcs_bytes(
    sip_a_2_0: float = 1.25e-5,
    date: str = "2026-07-21T00:00:00",
    history: str = "first run",
) -> bytes:
    cards = [
        fits_card("SIMPLE", True),
        fits_card("BITPIX", 8),
        fits_card("NAXIS", 0),
        fits_card("WCSAXES", 2),
        fits_card("CTYPE1", "RA---TAN-SIP"),
        fits_card("CTYPE2", "DEC--TAN-SIP"),
        fits_card("CRPIX1", 512.0),
        fits_card("CRPIX2", 384.0),
        fits_card("CRVAL1", 123.456),
        fits_card("CRVAL2", -12.5),
        fits_card("CD1_1", -2.5e-4),
        fits_card("CD1_2", 1.0e-7),
        fits_card("CD2_1", 2.0e-7),
        fits_card("CD2_2", 2.5e-4),
        fits_card("A_ORDER", 2),
        fits_card("A_2_0", sip_a_2_0),
        fits_card("A_1_1", -2e-7),
        fits_card("B_ORDER", 2),
        fits_card("B_0_2", 3e-6),
        fits_card("AP_ORDER", 2),
        fits_card("AP_2_0", -1.25e-5),
        fits_card("BP_ORDER", 2),
        fits_card("BP_0_2", -3e-6),
        fits_card("IMAGEW", 1024),
        fits_card("IMAGEH", 768),
        fits_card("DATE", date),
        fits_card("HISTORY", comment=history),
        fits_card("COMMENT", comment="volatile comment"),
        fits_card("END"),
    ]
    header = "".join(cards).encode("ascii")
    return header + b" " * ((-len(header)) % 2880)


def synthetic_resource_record(wall: float = 0.01) -> dict[str, object]:
    return {
        "schema_version": 1,
        "measurement_method": "os.wait4",
        "measurement_boundary": "python-launch -> taskset-exec -> target-tree -> wait4-reap",
        "wall_clock": "time.perf_counter_ns",
        "command": ["/usr/bin/true"],
        "wait_status": 0,
        "returncode": 0,
        "timed_out": False,
        "wall_seconds": wall,
        "user_seconds": wall * 0.4,
        "system_seconds": wall * 0.1,
        "cpu_seconds": wall * 0.5,
        "cpu_percent": 50.0,
        "max_rss_kbytes": 1024,
        "major_faults": 0,
        "minor_faults": 10,
        "voluntary_context_switches": 1,
        "involuntary_context_switches": 0,
        "filesystem_inputs": 0,
        "filesystem_outputs": 0,
    }


class RunnerTests(unittest.TestCase):
    def test_common_argument_ownership(self) -> None:
        self.assertEqual(
            run_matrix.validate_common_args(["--new-fits=none", "--scale-low=1"]),
            ["--new-fits=none", "--scale-low=1"],
        )
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.validate_common_args(["--new-fits=/tmp/shared.new"])
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.validate_common_args(["-vv"])
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.validate_common_args(
                ["--scale-low=1", "/tmp/second-input.jpg"]
            )
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.validate_common_args(["--scale-low", "1"])

    def test_field_object_range_schema(self) -> None:
        self.assertEqual(
            run_matrix.parse_field_object_range("11-50", "objects-11-50"),
            ("11-50", 10, 50),
        )
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.parse_field_object_range("1,10", "bad")
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.parse_field_object_range("50-11", "bad")
        ranges = [
            {"id": "any-first-label", "objects": "21-30"},
            {"id": "another.range", "objects": "1-10"},
            {"id": "third_range", "objects": "101-120"},
        ]
        self.assertEqual(
            run_matrix.parse_field_object_ranges(ranges),
            [
                ("any-first-label", "21-30", 20, 30),
                ("another.range", "1-10", 0, 10),
                ("third_range", "101-120", 100, 120),
            ],
        )
        invalid_values = (
            [{"id": "only", "objects": "1-10"}],
            [
                {"id": "same", "objects": "1-10"},
                {"id": "same", "objects": "11-20"},
            ],
            [
                {"id": "first", "objects": "1-10"},
                {"id": "second", "objects": "1-10"},
            ],
            [
                {"id": "../unsafe", "objects": "1-10"},
                {"id": "safe", "objects": "11-20"},
            ],
            [
                {"id": "first", "objects": "1-10", "meaning": "early"},
                {"id": "second", "objects": "11-20"},
            ],
        )
        for value in invalid_values:
            with self.subTest(value=value), self.assertRaises(
                run_matrix.HarnessError
            ):
                run_matrix.parse_field_object_ranges(value)

    def test_measurement_mode_verbosity_and_constructed_command(self) -> None:
        solve_field = Path("/synthetic/solve-field")
        input_path = Path("/synthetic/input.jpg")
        indexes = [Path("/synthetic/index-1.fits"), Path("/synthetic/index-2.fits")]
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            for mode, verbose_count in (("timing", 0), ("phase", 1), ("detailed", 2)):
                command = run_matrix.build_solve_command(
                    solve_field,
                    run_dir,
                    input_path,
                    indexes,
                    "1-10",
                    ["--scale-low=1"],
                    mode,
                )
                self.assertEqual(command.count("--verbose"), verbose_count)
                command_path = run_dir / f"command-{mode}.json"
                run_matrix.write_json(command_path, {"argv": command})
                assertion = run_matrix.assert_constructed_command(
                    command_path,
                    {"solve_command": command},
                    indexes,
                    "1-10",
                    0,
                    10,
                    input_path,
                    mode,
                )
                self.assertEqual(assertion["candidate_count"], 2)
                self.assertEqual(assertion["verbose_count"], verbose_count)
            injected = run_matrix.build_solve_command(
                solve_field,
                run_dir,
                input_path,
                indexes,
                "1-10",
                ["/synthetic/second-input.jpg"],
                "timing",
            )
            injected_path = run_dir / "command-injected.json"
            run_matrix.write_json(injected_path, {"argv": injected})
            with self.assertRaises(run_matrix.HarnessError):
                run_matrix.assert_constructed_command(
                    injected_path,
                    {"solve_command": injected},
                    indexes,
                    "1-10",
                    0,
                    10,
                    input_path,
                    "timing",
                )

    def test_rotated_plan_is_deterministic_and_complete(self) -> None:
        range_ids = ["objects-1-10", "objects-11-50", "objects-51-100"]
        first = run_matrix.condition_plan(range_ids, 7, 1234, 5)
        second = run_matrix.condition_plan(range_ids, 7, 1234, 5)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 7 * 36)
        expected = {
            (indexes, range_id, workers)
            for indexes in run_matrix.INDEX_COUNTS
            for range_id in range_ids
            for workers in run_matrix.WORKER_COUNTS
        }
        for repetition in range(1, 8):
            observed = {
                (item["index_count"], item["range_id"], item["workers"])
                for item in first
                if item["repetition"] == repetition
            }
            self.assertEqual(observed, expected)
        with self.assertRaises(run_matrix.HarnessError):
            run_matrix.condition_plan(range_ids, 7, 1234, 6)

    def test_workload_fingerprint_excludes_only_mode_and_output(self) -> None:
        provenance = {
            name: {"path": f"/synthetic/{name}", "sha256": name * 4}
            for name in (
                "solve_field",
                "wcsinfo",
                "input",
                "manifest_file",
                "taskset_binary",
                "launcher_probe_binary",
            )
        }
        provenance.update(
            {
                "indexes": [{"path": "/synthetic/index", "sha256": "i" * 64}],
                "pipeline_files": [
                    {"path": "/synthetic/engine", "sha256": "e" * 64}
                ],
                "harness": {
                    "run_matrix.py": {"sha256": "r" * 64},
                    "summarize.py": {"sha256": "s" * 64},
                    "fits_wcs_signature.py": {"sha256": "w" * 64},
                },
            }
        )
        plan = [{"index_count": 1, "range_id": "objects-1-10", "workers": 1}]
        timing = {
            "label": "same-workload",
            "measurement_mode": "timing",
            "output_root": "/results/timing",
            "field_object_ranges": [
                {"id": "objects-1-10", "objects": "1-10"},
                {"id": "objects-11-50", "objects": "11-50"},
            ],
        }
        phase = {
            **timing,
            "measurement_mode": "phase",
            "output_root": "/results/phase",
        }
        self.assertEqual(
            run_matrix.compute_workload_fingerprint(timing, provenance, plan),
            run_matrix.compute_workload_fingerprint(phase, provenance, plan),
        )
        changed = {**phase, "label": "different-workload"}
        self.assertNotEqual(
            run_matrix.compute_workload_fingerprint(timing, provenance, plan),
            run_matrix.compute_workload_fingerprint(changed, provenance, plan),
        )

    def test_wait4_resource_parser_and_live_child(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "resource-usage.json"
            run_matrix.write_json(path, synthetic_resource_record(0.004321))
            parsed = run_matrix.parse_resource_usage(path)
        self.assertEqual(parsed["wall_seconds"], 0.004321)
        self.assertEqual(parsed["cpu_seconds"], 0.0021605)
        live = run_matrix._launch_and_measure(
            ["/usr/bin/true"],
            subprocess.DEVNULL,
            run_matrix.benchmark_environment(1),
            2.0,
        )
        self.assertEqual(live["measurement_method"], "os.wait4")
        self.assertFalse(live["timed_out"])
        self.assertEqual(live["returncode"], 0)
        self.assertGreater(live["wall_seconds"], 0.0)
        timed_out = run_matrix._launch_and_measure(
            ["/bin/sh", "-c", "sleep 5"],
            subprocess.DEVNULL,
            run_matrix.benchmark_environment(1),
            0.02,
        )
        self.assertTrue(timed_out["timed_out"])
        self.assertNotEqual(timed_out["returncode"], 0)

    def test_canonical_wcs_signature_tracks_sip_not_volatile_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline.wcs"
            metadata_changed = root / "metadata.wcs"
            sip_changed = root / "sip.wcs"
            baseline.write_bytes(synthetic_wcs_bytes())
            metadata_changed.write_bytes(
                synthetic_wcs_bytes(
                    date="2027-01-01T12:34:56", history="different invocation"
                )
            )
            sip_changed.write_bytes(synthetic_wcs_bytes(sip_a_2_0=1.30e-5))
            self.assertEqual(wcs_signature(baseline), wcs_signature(metadata_changed))
            self.assertNotEqual(wcs_signature(baseline), wcs_signature(sip_changed))

    def test_resume_evidence_recomputes_canonical_wcs_signature(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            resource_path = run_dir / "resource-usage.json"
            run_matrix.write_json(resource_path, synthetic_resource_record())
            (run_dir / "solve.solved").write_bytes(b"\x01")
            (run_dir / "wcsinfo.txt").write_bytes(b"stable diagnostic text\n")
            wcs_path = run_dir / "solve.wcs"
            wcs_path.write_bytes(synthetic_wcs_bytes())
            record = {
                "solved": True,
                "time": run_matrix.parse_resource_usage(resource_path),
                "correctness_signature": wcs_signature(wcs_path),
                "correctness_signature_method": SIGNATURE_METHOD,
                "wcsinfo_sha256": run_matrix.sha256_file(run_dir / "wcsinfo.txt"),
                "raw_wcs_sha256": run_matrix.sha256_file(wcs_path),
            }
            run_matrix.verify_result_evidence(run_dir, record)

            wcs_path.write_bytes(
                synthetic_wcs_bytes(
                    date="2030-12-31T23:59:59", history="new volatile metadata"
                )
            )
            record["raw_wcs_sha256"] = run_matrix.sha256_file(wcs_path)
            run_matrix.verify_result_evidence(run_dir, record)

            wcs_path.write_bytes(synthetic_wcs_bytes(sip_a_2_0=9.9e-5))
            record["raw_wcs_sha256"] = run_matrix.sha256_file(wcs_path)
            with self.assertRaises(run_matrix.HarnessError):
                run_matrix.verify_result_evidence(run_dir, record)

    def test_serial_and_parallel_log_assertions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "timing"
            root.mkdir()
            indexes = []
            for number in range(4):
                index = root / f"index-{number}.fits"
                index.write_bytes(b"index")
                indexes.append(index.resolve())
            index_lines = "".join(f"  {path}\n" for path in indexes[:2])
            common = (
                "solver run parameters:\n"
                "indexes:\n"
                f"{index_lines}"
                "fieldfname field.axy\n"
                "fields 1\n"
                "startdepth 0\n"
                "enddepth 10\n"
            )
            serial = root / "serial.log"
            serial.write_text(common, encoding="utf-8")
            assertion = run_matrix.assert_log_shape(
                serial, indexes[:2], 1, 0, 10
            )
            self.assertEqual(assertion["candidate_count"], 2)
            self.assertEqual(
                assertion["scheduler_classification"], "serial"
            )
            shifted_serial = root / "shifted-serial.log"
            shifted_serial.write_text(
                common.replace("startdepth 0\n", "startdepth 4\n")
                + "[solver-geometry] mode=shared-readonly objects=10 "
                "valid_pairs=1 inbox_entries=1 estimated=1 allocated=1\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                shifted_serial, indexes[:2], 1, 4, 10
            )
            self.assertEqual(
                assertion["geometry_classification"], "shared-readonly"
            )

            def submission(
                candidates: int, startobj: int = 0, endobj: int = 10
            ) -> str:
                return (
                    "[index-shard] pthread-pool submit "
                    "workers=4 pool_workers=4 "
                    f"candidates={candidates} "
                    "engine_pass=0 depth_index=0 scale_index=0 "
                    f"startobj={startobj} endobj={endobj} "
                    "scheduler=ordered chunk=1 "
                    "inner_scheduler=dynamic-pool-lending "
                    "mmap_pass=0 mmap_advice=normal\n"
                )

            parallel = root / "parallel.log"
            parallel.write_text(
                common
                + submission(2)
                + "[index-shard] assist-lane state=publish "
                "owner_worker=0 index_order=0\n"
                + "[index-shard] assist-lane state=publish "
                "owner_worker=1 index_order=1\n"
                + "[index-shard] assist-lane state=join "
                "helper_worker=2 helper_slot=1 index_order=1\n"
                + "[index-shard] assist-lane state=leave "
                "helper_worker=2 helper_slot=1 index_order=1 rc=0\n"
                + "[index-shard] assist-lane state=unpublish "
                "owner_worker=0 index_order=0\n"
                + "[index-shard] assist-lane state=unpublish "
                "owner_worker=1 index_order=1\n"
                + "[index-shard] assist-pass generation=1 loans=1 "
                "notify_broadcasts=2 notify_skipped=9 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes[:2], 4, 0, 10
            )
            self.assertEqual(len(assertion["pthread_submissions"]), 1)
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )
            self.assertEqual(
                assertion["geometry_classification"], "not-requested"
            )
            self.assertEqual(assertion["assist_pass"]["loans"], 1)
            self.assertEqual(
                assertion["assistance_observation"], "loan-observed"
            )

            single_common = (
                "solver run parameters:\n"
                "indexes:\n"
                f"  {indexes[0]}\n"
                "fieldfname field.axy\n"
                "fields 1\n"
                "startdepth 0\n"
                "enddepth 10\n"
            )
            parallel.write_text(
                single_common
                + submission(1)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes[:1], 4, 0, 10
            )
            self.assertEqual(len(assertion["pthread_submissions"]), 1)
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )

            parallel.write_text(
                single_common
                + submission(1)
                + "[index-shard] assist-lane state=publish "
                "owner_worker=0 index_order=0\n"
                + "[index-shard] assist-lane state=join "
                "helper_worker=2 helper_slot=1 index_order=0\n"
                + "[index-shard] assist-lane state=leave "
                "helper_worker=2 helper_slot=1 index_order=0 rc=0\n"
                + "[index-shard] assist-lane state=unpublish "
                "owner_worker=0 index_order=0\n"
                + "[index-shard] assist-pass generation=1 loans=1 "
                "notify_broadcasts=2 notify_skipped=9 waiters=0\n"
                + "[solver-ab-phase] index=index.fits object=10 "
                "phase=diagonal mode=assisted combinations=32 "
                "codekd_queries=32 codekd_hits=4 candidates=2 "
                "verifications=1 blocks_planned=4 blocks_retired=4 "
                "blocks_owner=2 segments_retired=8 "
                "segment_payload_bytes=512 wall=0.1 user=0.12 "
                "system=0.01 major_faults=2 "
                "resource=process-overlap\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes[:1], 4, 0, 10
            )
            self.assertEqual(
                assertion["phase_mode_counts"]["assisted"], 1
            )
            self.assertEqual(
                assertion["phase_resource_scope"], "process-overlap"
            )
            self.assertEqual(
                assertion["assistance_observation"], "assisted"
            )

            shifted_common = single_common.replace(
                "startdepth 0\n", "startdepth 4\n"
            )
            parallel.write_text(
                shifted_common
                + "[solver-geometry] mode=legacy reason=budget\n"
                + submission(1, 4)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes[:1], 4, 4, 10
            )
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )
            self.assertEqual(
                assertion["geometry_classification"], "refused-budget"
            )

            outer_common = (
                "solver run parameters:\n"
                "indexes:\n"
                + "".join(f"  {path}\n" for path in indexes)
                + "fieldfname field.axy\n"
                "fields 1\n"
                "startdepth 0\n"
                "enddepth 10\n"
            )
            parallel.write_text(
                outer_common
                + submission(4)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes, 4, 0, 10
            )
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )
            self.assertEqual(
                assertion["geometry_classification"], "not-requested"
            )

            outer_prefix_common = outer_common.replace(
                "startdepth 0\n", "startdepth 4\n"
            )
            parallel.write_text(
                outer_prefix_common
                + "[solver-geometry] mode=shared-readonly objects=10 "
                "valid_pairs=1 inbox_entries=1 estimated=1 allocated=1\n"
                + submission(4, 4)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes, 4, 4, 10
            )
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )
            self.assertEqual(
                assertion["geometry_classification"], "shared-readonly"
            )

            parallel.write_text(
                outer_prefix_common
                + "[solver-geometry] mode=legacy reason=budget\n"
                + submission(4, 4)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n",
                encoding="utf-8",
            )
            assertion = run_matrix.assert_log_shape(
                parallel, indexes, 4, 4, 10
            )
            self.assertEqual(
                assertion["scheduler_classification"],
                "dynamic-pool-lending",
            )
            self.assertEqual(
                assertion["geometry_classification"], "refused-budget"
            )

            serial.write_text(
                common
                + f"Trying index {indexes[0]}...\n"
                + "Field 1: solved with index index-0.fits.\n",
                encoding="utf-8",
            )
            winner = run_matrix.extract_winner(
                serial, True, indexes[:2], 1
            )
            self.assertEqual(winner["position_one_based"], 1)
            parallel.write_text(
                single_common
                + submission(1)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n"
                + "[index-shard] reduce solved index_order=0 "
                "best_logodds=42 field=1\n",
                encoding="utf-8",
            )
            winner = run_matrix.extract_winner(
                parallel, True, indexes[:1], 4
            )
            self.assertEqual(winner["position_one_based"], 1)
            parallel.write_text(
                common
                + submission(2)
                + "[index-shard] assist-pass generation=1 loans=0 "
                "notify_broadcasts=0 notify_skipped=2 waiters=0\n"
                + "[index-shard] reduce solved index_order=1 best_logodds=42 field=1\n",
                encoding="utf-8",
            )
            winner = run_matrix.extract_winner(
                parallel, True, indexes[:2], 4
            )
            self.assertEqual(winner["path"], str(indexes[1]))

    def test_profile_parser(self) -> None:
        content = """\
[solve-field-profile] source=0.1 engine=0.2 output=0.3 total=0.6 just_augment=0
[onefield-profile] mode=pthread-dynamic-lending candidates=2 serial_executed=0 acquire=0.01 solver=0.19 release=0.01 output=0 total=0.21 failed=0 cancelled=0
[solver-ab-phase] index=index.fits object=11 phase=diagonal mode=assisted combinations=32 codekd_queries=32 codekd_hits=4 candidates=2 verifications=1 blocks_planned=4 blocks_retired=4 blocks_owner=2 segments_retired=8 segment_payload_bytes=512 wall=0.1 user=0.12 system=0.01 major_faults=2 resource=process-overlap hypothesis_order=0000000000000001 kd_result_order=0000000000000002 candidate_order=0000000000000003
[solver] phase-profile detailed=1 failed=0 solver_run_elapsed=0.19 codekd_work_wall_sum=0.1 codekd_calls=40 resolve_work_wall_sum=0.05 reduction_ex_verify_hit_work_wall_sum=0.03 resolve_calls=2 verify_hit_work_wall_sum=0.02 verify_calls=1 hypothesis_wave_elapsed_sum=1
[solver] verification-context enabled=1 waves=2 tiles=8 parallel_tiles=6 observed_parallel_tiles=5 serial_tiles=2 single_search_waves=1 opt_in_disabled_waves=0 datalog_serial_waves=0 raw_candidates=900 candidates=700 executed=700 committed=700 discarded=0 prepared_off_owner_thread=500 committed_from_off_owner_preparation=500 commit_thread_violations=0 winning_candidates=1 winning_from_parallel_tile=1 winning_from_off_owner_preparation=0 tune_reverify_commits=0 tune_reverify_from_off_owner_preparation=0 discarded_after_solve=12 task_ranges=26 tasks_executed=26 submitted=18 inline=8 serial_fallbacks=0 failures=0 max_candidates=400 max_tile_candidates=256 max_tasks=4 max_parallel=4
[index-shard] aux-pass pool_workers=4 outer_workers=4 candidates=2 detailed=1 idle_tasks=3 lender_tasks=0 owner_tasks=0 total_work_wall_sum=0.001
[index-shard] assist-pass generation=1 loans=1 notify_broadcasts=2 notify_skipped=9 waiters=0
[index-shard] reducer-pass generation=1 candidates=2 calls=1 work_wall_sum=0.001
[index-shard] context-pass generation=1 candidates=2 outer_workers=4 prepare_work_wall_sum=0.002 prepare_max=0.001 cleanup_work_wall_sum=0.001 cleanup_max=0.0005
[index-shard] solver-pass generation=1 candidates=2 detailed=1 codekd_work_wall_sum=0.1 hypotheses=40 helper_tasks=1 helper_combinations=32
[index-shard] phase-profile executed=2 task_work_wall_sum=0.2 reset_work_wall_sum=0.01 reset_percent=5.0 acquire_work_wall_sum=0.02 acquire_percent=10.0 solve_work_wall_sum=0.14 solve_percent=70.0 analyze_work_wall_sum=0.01 analyze_percent=5.0 release_work_wall_sum=0.01 release_percent=5.0 other_work_wall_sum=0.01 other_percent=5.0 acquire_p50=0.01
[index-shard] pass-resource user=0.2 sys=0.03 minflt=20 majflt=2 nvcsw=4 nivcsw=1 inblock=16 oublock=0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "solve.log"
            path.write_text(content, encoding="utf-8")
            profiles = run_matrix.parse_profiles(path)
        self.assertEqual(profiles["solver_phase"][0]["codekd_calls"], 40)
        self.assertEqual(
            profiles["solver_ab_phase"][0]["mode"], "assisted"
        )
        self.assertEqual(
            profiles["solver_ab_phase"][0]["resource"],
            "process-overlap",
        )
        self.assertEqual(profiles["index_shard_assist"][0]["loans"], 1)
        self.assertEqual(profiles["index_shard_resource"][0]["majflt"], 2)
        self.assertEqual(
            profiles["solver_phase"][0]["reduction_ex_verify_hit_work_wall_sum"],
            0.03,
        )
        self.assertEqual(
            profiles["solver_verification"][0]["parallel_tiles"], 6
        )
        self.assertEqual(
            profiles["solver_verification"][0]["max_tile_candidates"], 256
        )
        self.assertEqual(
            profiles["solver_verification"][0]["single_search_waves"], 1
        )
        self.assertEqual(
            profiles["solver_verification"][0][
                "committed_from_off_owner_preparation"
            ],
            500,
        )
        self.assertEqual(
            profiles["solver_verification"][0]["commit_thread_violations"],
            0,
        )
        self.assertEqual(profiles["index_shard_aux"][0]["idle_tasks"], 3)
        self.assertEqual(profiles["index_shard_reducer"][0]["calls"], 1)
        self.assertEqual(profiles["index_shard_reducer"][0]["work_wall_sum"], 0.001)
        self.assertEqual(profiles["index_shard_context"][0]["outer_workers"], 4)
        self.assertEqual(profiles["index_shard_solver"][0]["hypotheses"], 40)
        self.assertEqual(
            profiles["index_shard_phase"][0]["acquire_work_wall_sum"], 0.02
        )
        self.assertEqual(profiles["index_shard_phase"][0]["acquire_percent"], 10.0)
        self.assertEqual(run_matrix.parse_profile_value("0.02(10.0%)"), 0.02)

    def test_resume_pair_is_bound_to_suite_fingerprint(self) -> None:
        item = {
            "sequence": 1,
            "repetition": 1,
            "order_in_repetition": 1,
            "index_count": 1,
            "range_id": "objects-1-10",
            "workers": 1,
        }
        with tempfile.TemporaryDirectory() as directory:
            pair_dir = Path(directory) / "pair"
            for state in run_matrix.CACHE_STATES:
                state_dir = pair_dir / state
                state_dir.mkdir(parents=True)
                evidence_contents = {
                    "command.json": b"{}\n",
                    "solve.log": b"evidence\n",
                    "resource-usage.json": run_matrix.canonical_json(
                        synthetic_resource_record()
                    ),
                    "wcsinfo.txt": b"UNSOLVED\n",
                }
                if state == "cold":
                    evidence_contents["cache-control.json"] = b"{}\n"
                for name, content in evidence_contents.items():
                    (state_dir / name).write_bytes(content)
                (state_dir / "run.json").write_text(
                    json.dumps(
                        {
                            "complete": True,
                            "suite_fingerprint": "suite-a",
                            "workload_fingerprint": "workload-a",
                            "cache_state": state,
                            "assertions": {"ok": True},
                            "solved": False,
                            "correctness_signature": "UNSOLVED",
                            "correctness_signature_method": SIGNATURE_METHOD,
                            "wcsinfo_sha256": run_matrix.sha256_file(
                                state_dir / "wcsinfo.txt"
                            ),
                            "raw_wcs_sha256": None,
                            "time": run_matrix.parse_resource_usage(
                                state_dir / "resource-usage.json"
                            ),
                            "evidence_sha256": {
                                name: run_matrix.sha256_file(state_dir / name)
                                for name in evidence_contents
                            },
                            **item,
                        }
                    ),
                    encoding="utf-8",
                )
            marker = pair_dir / "pair.json"
            marker.write_text(
                json.dumps(
                    {
                        "complete": True,
                        "suite_fingerprint": "suite-a",
                        "workload_fingerprint": "workload-a",
                        **item,
                    }
                ),
                encoding="utf-8",
            )
            self.assertTrue(
                run_matrix.pair_complete(
                    pair_dir, item, "suite-a", "workload-a"
                )
            )
            with self.assertRaises(run_matrix.HarnessError):
                run_matrix.pair_complete(
                    pair_dir, item, "suite-b", "workload-a"
                )
            with self.assertRaises(run_matrix.HarnessError):
                run_matrix.pair_complete(
                    pair_dir, item, "suite-a", "workload-b"
                )
            (pair_dir / "warm" / "solve.log").write_text(
                "altered evidence\n", encoding="utf-8"
            )
            with self.assertRaises(run_matrix.HarnessError):
                run_matrix.pair_complete(
                    pair_dir, item, "suite-a", "workload-a"
                )


class SummaryTests(unittest.TestCase):
    def test_median_mad_and_paired_bootstrap(self) -> None:
        self.assertEqual(summarize.median_mad([1.0, 2.0, 100.0]), (2.0, 1.0))
        observed, low, high = summarize.paired_bootstrap_median_ratio(
            [2.0] * 7, [1.0] * 7, 5000, 99
        )
        self.assertEqual((observed, low, high), (2.0, 2.0, 2.0))

    def test_summary_field_object_range_schema(self) -> None:
        value = [
            {"id": "later-label-first", "objects": "21-30"},
            {"id": "earlier-label-second", "objects": "1-10"},
            {"id": "third", "objects": "31-40"},
        ]
        range_ids, errors = summarize.validate_field_object_ranges(value)
        self.assertEqual(
            range_ids,
            ["later-label-first", "earlier-label-second", "third"],
        )
        self.assertEqual(errors, [])
        _, errors = summarize.validate_field_object_ranges(value[:1])
        self.assertTrue(errors)
        _, errors = summarize.validate_field_object_ranges(
            [
                {"id": "duplicate-spec-a", "objects": "1-10"},
                {"id": "duplicate-spec-b", "objects": "1-10"},
            ]
        )
        self.assertTrue(any("duplicate" in error for error in errors))

    def test_quiet_timing_rejects_structured_profiles(self) -> None:
        run = {
            "_record_path": "/synthetic/run.json",
            "profile_validation": "not_observable_quiet",
            "assertions": {
                "observed_log_shape": "not_observable_quiet",
                "verbose_count": 0,
            },
            "profiles": {"solve_field": [{"source": 0.1}]},
        }
        errors = summarize.validate_profile_evidence([run], "timing")
        self.assertTrue(any("structured profiles" in error for error in errors))

    def test_scheduler_classification_is_topology_bound_and_repeatable(
        self,
    ) -> None:
        def classified(
            index_count: int,
            workers: int,
            classification: str,
            repetition: int,
        ) -> dict:
            return {
                "_record_path": f"/synthetic/r{repetition}/run.json",
                "index_count": index_count,
                "range_id": "arbitrary-band",
                "workers": workers,
                "assertions": {
                    "scheduler_classification": classification,
                },
            }

        classifications, errors = (
            summarize.collect_scheduler_classifications(
                [
                    classified(1, 4, "dynamic-pool-lending", 1),
                    classified(1, 4, "dynamic-pool-lending", 2),
                    classified(4, 4, "dynamic-pool-lending", 1),
                ],
                "phase",
            )
        )
        self.assertEqual(errors, [])
        self.assertEqual(
            classifications[(1, "arbitrary-band", 4)],
            "dynamic-pool-lending",
        )

        _, errors = summarize.collect_scheduler_classifications(
            [
                classified(1, 4, "dynamic-pool-lending", 1),
                classified(1, 4, "serial", 2),
            ],
            "detailed",
        )
        self.assertTrue(
            any("inconsistent across repetitions" in error for error in errors)
        )

        _, errors = summarize.collect_scheduler_classifications(
            [classified(4, 4, "serial", 1)],
            "phase",
        )
        self.assertTrue(
            any("scheduler classification mismatch" in error for error in errors)
        )

    def test_complete_synthetic_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "timing"
            root.mkdir()
            workload_fingerprint = "a" * 64
            saved_profiles = {}
            (root / "provenance").mkdir()
            (root / "provenance" / "provenance.json").write_text(
                json.dumps(
                    {
                        "workload_fingerprint": workload_fingerprint,
                        "config": {
                            "bootstrap_samples": 5000,
                            "bootstrap_seed": 77,
                        }
                    }
                ),
                encoding="utf-8",
            )
            run_count = 0
            order = 0
            field_object_ranges = (
                ("objects-11-50", "11-50", 10, 50),
                ("objects-1-10", "1-10", 0, 10),
            )
            for repetition in range(1, 8):
                for index_count in (1, 2, 4, 8):
                    for (
                        range_id,
                        field_objects,
                        startobj,
                        endobj,
                    ) in field_object_ranges:
                        for workers in (1, 2, 4):
                            order += 1
                            for cache_state in ("cold", "warm"):
                                run_count += 1
                                run_dir = (
                                    root
                                    / "runs"
                                    / f"r{repetition:02d}"
                                    / f"o{order:03d}"
                                    / cache_state
                                )
                                run_dir.mkdir(parents=True)
                                wall = {1: 2.0, 2: 1.5, 4: 1.0}[workers]
                                wall += repetition * 0.001
                                time_record = {
                                    "wall_seconds": wall,
                                    "user_seconds": wall * workers * 0.8,
                                    "system_seconds": 0.1,
                                    "cpu_seconds": wall * workers * 0.8 + 0.1,
                                    "cpu_percent": workers * 90,
                                    "max_rss_kbytes": 10000 + workers,
                                    "major_faults": 0,
                                    "minor_faults": 100,
                                    "voluntary_context_switches": 2,
                                    "involuntary_context_switches": 1,
                                    "filesystem_inputs": 0,
                                    "filesystem_outputs": 8,
                                    "returncode": 0,
                                    "wait_status": 0,
                                    "timed_out": False,
                                }
                                use_pthread = workers > 1
                                scheduler_classification = (
                                    "serial"
                                    if workers == 1
                                    else "dynamic-pool-lending"
                                )
                                assistance_observed = use_pthread
                                field_profile_count = (
                                    1
                                    if use_pthread
                                    else min(index_count, 2)
                                )
                                record = {
                                    "complete": True,
                                    "suite_fingerprint": "synthetic-suite",
                                    "workload_fingerprint": workload_fingerprint,
                                    "repetition": repetition,
                                    "index_count": index_count,
                                    "range_id": range_id,
                                    "workers": workers,
                                    "cache_state": cache_state,
                                    "assertions": {
                                        "ok": True,
                                        "validation": "constructed_command",
                                        "candidate_count": index_count,
                                        "candidate_paths": [
                                            f"/synthetic/index-{number}.fits"
                                            for number in range(1, index_count + 1)
                                        ],
                                        "field_objects": field_objects,
                                        "startobj": startobj,
                                        "endobj": endobj,
                                        "input": "/synthetic/input.jpg",
                                        "verbose_count": 0,
                                        "observed_log_shape": "not_observable_quiet",
                                    },
                                    "profile_validation": "not_observable_quiet",
                                    "solved": True,
                                    "time": time_record,
                                    "profiles": {
                                        "solve_field": [
                                            {
                                                "source": 0.2,
                                                "engine": wall - 0.2,
                                                "output": 0.0,
                                                "total": wall,
                                                "just_augment": 0,
                                                "raw": "ignored",
                                            }
                                        ],
                                        "engine": [
                                            {
                                                "pool_start": 0.0,
                                                "pool_stop": 0.0,
                                                "engine_total": wall,
                                                "solver_failed": 0,
                                            }
                                        ],
                                        "onefield": [
                                            {
                                                "mode": (
                                                    "serial"
                                                    if not use_pthread
                                                    else "pthread-dynamic-lending"
                                                ),
                                                "candidates": index_count,
                                                "serial_executed": (
                                                    0
                                                    if use_pthread
                                                    else field_profile_count
                                                ),
                                                "acquire": 0.01,
                                                "solver": wall - 0.2,
                                                "release": 0.01,
                                                "output": 0.01,
                                                "total": wall - 0.1,
                                                "failed": 0,
                                                "cancelled": 0,
                                            }
                                        ],
                                        "onefield_field": [
                                            {
                                                "field": 1,
                                                "read": 0.01,
                                                "preprocess": 0.01,
                                                "solver_run": wall - 0.2,
                                                "total": wall - 0.1,
                                                "failed": 0,
                                            }
                                            for _ in range(field_profile_count)
                                        ],
                                    },
                                }
                                if use_pthread:
                                    record["profiles"].update(
                                        {
                                            "index_shard_assist": [
                                                {
                                                    "generation": repetition,
                                                    "loans": (
                                                        1
                                                        if assistance_observed
                                                        else 0
                                                    ),
                                                    "notify_broadcasts": 1,
                                                    "notify_skipped": 3,
                                                    "waiters": 0,
                                                }
                                            ],
                                            "index_shard_reducer": [
                                                {
                                                    "generation": repetition,
                                                    "candidates": index_count,
                                                    "calls": 1,
                                                    "work_wall_sum": 0.01,
                                                }
                                            ],
                                            "index_shard_context": [
                                                {
                                                    "generation": repetition,
                                                    "candidates": index_count,
                                                    "outer_workers": workers,
                                                    "prepare_work_wall_sum": 0.01,
                                                    "prepare_max": 0.01,
                                                    "cleanup_work_wall_sum": 0.01,
                                                    "cleanup_max": 0.01,
                                                }
                                            ],
                                            "index_shard_phase": [
                                                {
                                                    "executed": 1,
                                                    "task_work_wall_sum": wall - 0.1,
                                                    "reset_work_wall_sum": 0.01,
                                                    "reset_percent": 1.0,
                                                    "acquire_work_wall_sum": 0.01,
                                                    "acquire_percent": 1.0,
                                                    "solve_work_wall_sum": wall - 0.2,
                                                    "solve_percent": 90.0,
                                                    "analyze_work_wall_sum": 0.01,
                                                    "analyze_percent": 1.0,
                                                    "release_work_wall_sum": 0.01,
                                                    "release_percent": 1.0,
                                                    "other_work_wall_sum": 0.06,
                                                    "other_percent": 6.0,
                                                }
                                            ],
                                            "index_shard_solver": [
                                                {
                                                    "generation": repetition,
                                                    "candidates": index_count,
                                                    "detailed": 0,
                                                    "failed": 0,
                                                    "solver_run_work_wall_sum": wall - 0.2,
                                                    "codekd_work_wall_sum": 0.0,
                                                    "resolve_work_wall_sum": 0.0,
                                                    "reduction_ex_verify_hit_work_wall_sum": 0.0,
                                                    "verify_hit_work_wall_sum": 0.0,
                                                    "hypothesis_wave_elapsed_sum": 0.0,
                                                    "observed_parallel": (
                                                        1
                                                        if assistance_observed
                                                        else 0
                                                    ),
                                                    "parallel_hypotheses": (
                                                        32
                                                        if assistance_observed
                                                        else 0
                                                    ),
                                                    "helper_tasks": (
                                                        1
                                                        if assistance_observed
                                                        else 0
                                                    ),
                                                    "helper_combinations": (
                                                        32
                                                        if assistance_observed
                                                        else 0
                                                    ),
                                                }
                                            ],
                                            "index_shard_resource": [
                                                {
                                                    "user": wall * workers * 0.8,
                                                    "sys": 0.1,
                                                    "minflt": 100,
                                                    "majflt": 0,
                                                    "nvcsw": 2,
                                                    "nivcsw": 1,
                                                    "inblock": 0,
                                                    "oublock": 8,
                                                }
                                            ],
                                        }
                                    )
                                evidence_contents = {
                                    "command.json": b"{}\n",
                                    "solve.log": b"synthetic log\n",
                                    "resource-usage.json": run_matrix.canonical_json(
                                        {
                                            **synthetic_resource_record(wall),
                                            **time_record,
                                        }
                                    ),
                                    "wcsinfo.txt": (
                                        f"synthetic wcsinfo i{index_count} "
                                        f"{range_id}\n".encode()
                                    ),
                                    "solve.solved": b"\x01",
                                    "solve.wcs": synthetic_wcs_bytes(),
                                }
                                if cache_state == "cold":
                                    evidence_contents["cache-control.json"] = b"{}\n"
                                for name, content in evidence_contents.items():
                                    (run_dir / name).write_bytes(content)
                                record["evidence_sha256"] = {
                                    name: summarize.sha256_file(run_dir / name)
                                    for name in evidence_contents
                                }
                                record["correctness_signature"] = wcs_signature(
                                    run_dir / "solve.wcs"
                                )
                                record["correctness_signature_method"] = SIGNATURE_METHOD
                                record["wcsinfo_sha256"] = record["evidence_sha256"][
                                    "wcsinfo.txt"
                                ]
                                record["raw_wcs_sha256"] = record[
                                    "evidence_sha256"
                                ]["solve.wcs"]
                                record_path = run_dir / "run.json"
                                saved_profiles[str(record_path)] = record["profiles"]
                                record["profiles"] = {}
                                record_path.write_text(
                                    json.dumps(record), encoding="utf-8"
                                )
            (root / "suite.json").write_text(
                json.dumps(
                    {
                        "complete": True,
                        "fingerprint": "synthetic-suite",
                        "workload_fingerprint": workload_fingerprint,
                        "cell_count": 48,
                        "pair_count": 168,
                        "run_count": run_count,
                        "repetitions": 7,
                        "index_counts": [1, 2, 4, 8],
                        "workers": [1, 2, 4],
                        "cache_states": ["cold", "warm"],
                        "field_object_ranges": [
                            {"id": range_id, "objects": field_objects}
                            for range_id, field_objects, _, _ in field_object_ranges
                        ],
                        "range_ids": [
                            range_id
                            for range_id, _, _, _ in field_object_ranges
                        ],
                        "measurement_mode": "timing",
                        "resource_measurement": "os.wait4",
                        "launcher_overhead_wall_median_seconds": 0.0005,
                    }
                ),
                encoding="utf-8",
            )
            returncode = summarize.main([str(root)])
            self.assertEqual(returncode, 0)
            summary = json.loads(
                (root / "summary" / "summary.json").read_text(encoding="utf-8")
            )
            self.assertTrue(summary["valid"])
            self.assertEqual(
                [
                    row["range_id"]
                    for row in summary["condition_summary"]
                    if row["index_count"] == 1
                    and row["workers"] == 1
                    and row["cache_state"] == "cold"
                ],
                ["objects-11-50", "objects-1-10"],
            )
            self.assertTrue(
                all(
                    row["crossover_flag"] == "PARALLEL"
                    for row in summary["crossover"]
                )
            )
            self.assertTrue(
                all(
                    row["scheduler_classification"]
                    == "not-observable-quiet"
                    for row in summary["condition_summary"]
                )
            )
            self.assertEqual(
                summary["companion_validation"]["status"], "NOT_CHECKED"
            )
            self.assertFalse(
                summary["companion_validation"]["timing_acceptance_ready"]
            )
            dynamic_rows = summary["helper_usage"]
            self.assertTrue(dynamic_rows)
            self.assertTrue(
                all(
                    row["flag"] == "NOT_OBSERVABLE_QUIET"
                    for row in dynamic_rows
                )
            )

            phase_root = Path(directory) / "phase"
            shutil.copytree(root, phase_root)
            suite_path = phase_root / "suite.json"
            phase_suite = json.loads(suite_path.read_text(encoding="utf-8"))
            phase_suite["measurement_mode"] = "phase"
            suite_path.write_text(json.dumps(phase_suite), encoding="utf-8")
            for record_path in (phase_root / "runs").glob("r*/*/*/run.json"):
                phase_record = json.loads(record_path.read_text(encoding="utf-8"))
                phase_record["profile_validation"] = "validated_verbose"
                phase_record["assertions"]["observed_log_shape"] = "validated_verbose"
                phase_record["assertions"]["verbose_count"] = 1
                phase_record["assertions"]["scheduler_classification"] = (
                    "serial"
                    if phase_record["workers"] == 1
                    else "dynamic-pool-lending"
                )
                phase_record["assertions"]["phase_mode_counts"] = {
                    "assisted": 0,
                    "empty": 0,
                    "flattened-owner": 0,
                    "native": 0,
                }
                timing_record_path = root / record_path.relative_to(phase_root)
                phase_record["profiles"] = saved_profiles[str(timing_record_path)]
                if phase_record["workers"] > 1:
                    phase_record["assertions"]["assist_pass"] = dict(
                        phase_record["profiles"]["index_shard_assist"][0]
                    )
                    phase_record["assertions"][
                        "assist_lane_publications"
                    ] = [
                        {
                            "state": "publish",
                            "owner_worker": 0,
                            "index_order": 0,
                        },
                        {
                            "state": "unpublish",
                            "owner_worker": 0,
                            "index_order": 0,
                        },
                    ]
                record_path.write_text(json.dumps(phase_record), encoding="utf-8")
            self.assertEqual(summarize.main([str(phase_root)]), 0)
            phase_summary = json.loads(
                (phase_root / "summary" / "summary.json").read_text(encoding="utf-8")
            )
            phase_snapshot_root = Path(directory) / "phase-snapshot"
            shutil.copytree(phase_root, phase_snapshot_root)
            self.assertTrue(
                all(
                    row["crossover_flag"] == "DIAGNOSTIC_ONLY"
                    for row in phase_summary["crossover"]
                )
            )
            self.assertTrue(
                all(
                    row["evidence_selected_workers"] is None
                    for row in phase_summary["crossover"]
                )
            )
            phase_dynamic_rows = phase_summary["helper_usage"]
            self.assertTrue(
                all(
                    row["flag"] == "DYNAMIC_ASSIST_OBSERVED"
                    for row in phase_dynamic_rows
                )
            )
            single_index_rows = [
                row
                for row in phase_summary["helper_usage"]
                if row["index_count"] == 1 and row["workers"] > 1
            ]
            self.assertTrue(
                all(
                    row["scheduler_classification"]
                    == "dynamic-pool-lending"
                    and row["flag"] == "DYNAMIC_ASSIST_OBSERVED"
                    for row in single_index_rows
                )
            )
            outer_rows = [
                row
                for row in phase_summary["helper_usage"]
                if row["index_count"] >= row["workers"]
            ]
            self.assertTrue(outer_rows)
            self.assertTrue(
                all(
                    row["scheduler_classification"]
                    == "dynamic-pool-lending"
                    and row["flag"] == "DYNAMIC_ASSIST_OBSERVED"
                    for row in outer_rows
                )
            )
            phase_text = (phase_root / "summary" / "phase_summary.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("solve_field\tsource", phase_text)
            self.assertEqual(
                summarize.main(
                    [str(root), "--companion", str(phase_root)]
                ),
                0,
            )
            compared_timing = json.loads(
                (root / "summary" / "summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                compared_timing["companion_validation"]["status"], "PASS"
            )
            self.assertTrue(
                compared_timing["companion_validation"]["timing_acceptance_ready"]
            )
            self.assertEqual(
                compared_timing["companion_validation"][
                    "scheduler_classification_source"
                ],
                "phase",
            )
            self.assertTrue(
                any(
                    row["index_count"] == 1
                    and row["workers"] == 4
                    and row["scheduler_classification"]
                    == "dynamic-pool-lending"
                    for row in compared_timing["companion_validation"][
                        "scheduler_classifications"
                    ]
                )
            )
            mismatched_suite = json.loads(suite_path.read_text(encoding="utf-8"))
            mismatched_suite["workload_fingerprint"] = "b" * 64
            suite_path.write_text(json.dumps(mismatched_suite), encoding="utf-8")
            current_suite, current_runs = summarize.load_runs(
                root, allow_incomplete=False
            )
            companion_validation, companion_errors = (
                summarize.validate_companion_roots(
                    root,
                    current_suite,
                    current_runs,
                    [phase_root],
                )
            )
            self.assertEqual(companion_validation["status"], "FAIL")
            self.assertTrue(
                any("workload fingerprint differs" in error for error in companion_errors)
            )
            mismatched_suite["workload_fingerprint"] = workload_fingerprint
            suite_path.write_text(json.dumps(mismatched_suite), encoding="utf-8")

            detailed_suite = json.loads(suite_path.read_text(encoding="utf-8"))
            detailed_suite["measurement_mode"] = "detailed"
            suite_path.write_text(json.dumps(detailed_suite), encoding="utf-8")
            for record_path in (phase_root / "runs").glob("r*/*/*/run.json"):
                detailed_record = json.loads(record_path.read_text(encoding="utf-8"))
                detailed_record["assertions"]["verbose_count"] = 2
                profiles = detailed_record["profiles"]
                field_count = len(profiles["onefield_field"])
                profiles["solver_phase"] = [
                    {
                        "detailed": 1,
                        "failed": 0,
                        "solver_run_elapsed": 0.5,
                        "codekd_work_wall_sum": 0.1,
                        "codekd_calls": 40,
                        "resolve_work_wall_sum": 0.2,
                        "reduction_ex_verify_hit_work_wall_sum": 0.1,
                        "resolve_calls": 4,
                        "verify_hit_work_wall_sum": 0.1,
                        "verify_calls": 2,
                        "hypothesis_wave_elapsed_sum": 0.3,
                    }
                    for _ in range(field_count)
                ]
                ab_mode = (
                    "native"
                    if detailed_record["workers"] == 1
                    else "assisted"
                )
                profiles["solver_ab_phase"] = [
                    {
                        "index": "synthetic-index.fits",
                        "object": 11,
                        "phase": "diagonal",
                        "mode": ab_mode,
                        "combinations": 32,
                        "codekd_queries": 32,
                        "codekd_hits": 4,
                        "candidates": 2,
                        "verifications": 1,
                        "blocks_planned": (
                            4 if ab_mode == "assisted" else 0
                        ),
                        "blocks_retired": (
                            4 if ab_mode == "assisted" else 0
                        ),
                        "blocks_owner": (
                            2 if ab_mode == "assisted" else 0
                        ),
                        "segments_retired": (
                            8 if ab_mode == "assisted" else 0
                        ),
                        "segment_payload_bytes": (
                            512 if ab_mode == "assisted" else 0
                        ),
                        "wall": 0.1,
                        "user": 0.1,
                        "system": 0.01,
                        "major_faults": 0,
                        "resource": "process-overlap",
                    }
                    for _ in range(field_count)
                ]
                detailed_record["assertions"]["phase_mode_counts"] = {
                    "assisted": (
                        field_count if ab_mode == "assisted" else 0
                    ),
                    "empty": 0,
                    "flattened-owner": 0,
                    "native": field_count if ab_mode == "native" else 0,
                }
                if "index_shard_solver" in profiles:
                    profiles["index_shard_solver"][0]["detailed"] = 1
                record_path.write_text(json.dumps(detailed_record), encoding="utf-8")
            self.assertEqual(summarize.main([str(phase_root)]), 0)
            detailed_summary = json.loads(
                (phase_root / "summary" / "summary.json").read_text(encoding="utf-8")
            )
            self.assertTrue(
                all(
                    row["crossover_flag"] == "DIAGNOSTIC_ONLY"
                    and row["evidence_selected_workers"] is None
                    for row in detailed_summary["crossover"]
                )
            )
            self.assertTrue(
                all(
                    row["runs_with_exact_phase_modes"] == row["runs"]
                    for row in detailed_summary["helper_usage"]
                )
            )
            current_suite, current_runs = summarize.load_runs(
                root, allow_incomplete=False
            )
            companion_validation, companion_errors = (
                summarize.validate_companion_roots(
                    root,
                    current_suite,
                    current_runs,
                    [phase_snapshot_root, phase_root],
                )
            )
            self.assertEqual(companion_errors, [])
            self.assertEqual(companion_validation["status"], "PASS")
            self.assertEqual(
                companion_validation[
                    "scheduler_classification_consistency"
                ],
                "PASS",
            )

            reclassified_paths = []
            for record_path in (phase_root / "runs").glob(
                "r*/*/*/run.json"
            ):
                record = json.loads(
                    record_path.read_text(encoding="utf-8")
                )
                if (
                    record["index_count"] == 1
                    and record["workers"] == 4
                ):
                    record["assertions"][
                        "scheduler_classification"
                    ] = "serial"
                    record["profiles"]["onefield"][0]["mode"] = "serial"
                    record_path.write_text(
                        json.dumps(record), encoding="utf-8"
                    )
                    reclassified_paths.append(record_path)
            companion_validation, companion_errors = (
                summarize.validate_companion_roots(
                    root,
                    current_suite,
                    current_runs,
                    [phase_snapshot_root, phase_root],
                )
            )
            self.assertEqual(companion_validation["status"], "FAIL")
            self.assertEqual(
                companion_validation[
                    "scheduler_classification_consistency"
                ],
                "NOT_CHECKED",
            )
            self.assertTrue(
                any(
                    "scheduler classification mismatch"
                    in error
                    for error in companion_errors
                )
            )
            for record_path in reclassified_paths:
                record = json.loads(
                    record_path.read_text(encoding="utf-8")
                )
                record["assertions"][
                    "scheduler_classification"
                ] = "dynamic-pool-lending"
                record["profiles"]["onefield"][0][
                    "mode"
                ] = "pthread-dynamic-lending"
                record_path.write_text(
                    json.dumps(record), encoding="utf-8"
                )
            mismatch_path = next(
                path
                for path in (phase_root / "runs").glob("r*/*/*/run.json")
                if json.loads(path.read_text(encoding="utf-8"))["workers"] == 2
            )
            mismatch_record = json.loads(mismatch_path.read_text(encoding="utf-8"))
            mismatch_wcs = mismatch_path.parent / "solve.wcs"
            mismatch_wcs.write_bytes(synthetic_wcs_bytes(sip_a_2_0=9.9e-5))
            mismatch_record["raw_wcs_sha256"] = summarize.sha256_file(mismatch_wcs)
            mismatch_record["evidence_sha256"]["solve.wcs"] = summarize.sha256_file(
                mismatch_wcs
            )
            mismatch_path.write_text(json.dumps(mismatch_record), encoding="utf-8")
            self.assertEqual(summarize.main([str(phase_root)]), 2)
            invalid_summary = json.loads(
                (phase_root / "summary" / "summary.json").read_text(encoding="utf-8")
            )
            self.assertFalse(invalid_summary["valid"])
            self.assertTrue(
                any(
                    "canonical WCS signature mismatch" in error
                    for error in invalid_summary["validation_errors"]
                )
            )
            suite, runs = summarize.load_runs(phase_root, allow_incomplete=False)
            runs[0]["index_count"] = 16
            errors = summarize.validate_runs(
                runs,
                7,
                suite["fingerprint"],
                suite["workload_fingerprint"],
                suite["range_ids"],
            )
            self.assertTrue(any("missing" in error for error in errors))
            self.assertTrue(any("unexpected" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
