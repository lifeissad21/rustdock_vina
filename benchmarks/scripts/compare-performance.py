#!/usr/bin/env python3
"""Compare two release builds on five targets; require identical seeded outputs."""

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--full-docking",
        action="store_true",
        help="Use normal exhaustiveness-1 docking without an evaluation cap",
    )
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    engines = {"before": args.before.resolve(), "after": args.after.resolve()}
    cases = json.loads(
        (ROOT / "benchmarks/datasets/dude-100/manifest.json").read_text()
    )["cases"]
    selected = []
    seen = set()
    for case in cases:
        if case["target"] not in seen:
            selected.append(case)
            seen.add(case["target"])
    report = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "repeats": args.repeats,
        "binaries": {
            k: {"path": str(v), "sha256": hashlib.sha256(v.read_bytes()).hexdigest()}
            for k, v in engines.items()
        },
        "results": [],
    }
    with tempfile.TemporaryDirectory() as temp:
        output = Path(temp) / "pose.pdbqt"
        for case in selected:
            for workload in (
                "map-preparation",
                "full-docking" if args.full_docking else "limited-docking",
            ):
                common = [
                    "--receptor",
                    case["receptor"],
                    "--ligand",
                    case["pdbqt"],
                    "--config",
                    case["config"],
                    "--metal",
                    "off",
                    "--cpu",
                    "1",
                    "--seed",
                    "42",
                    "--verbosity",
                    "0",
                    "--out",
                    str(output),
                ]
                common += (
                    ["--randomize_only"]
                    if workload == "map-preparation"
                    else ["--exhaustiveness", "1", "--num_modes", "3"]
                )
                if workload == "limited-docking":
                    common += ["--max_evals", "2000"]
                times = {k: [] for k in engines}
                hashes = {}
                # Warm each binary, then alternate order to reduce temperature/order bias.
                for repeat in range(-1, args.repeats):
                    order = (
                        list(engines) if repeat % 2 == 0 else list(reversed(engines))
                    )
                    for label in order:
                        start = time.perf_counter()
                        completed = subprocess.run(
                            [str(engines[label]), *common],
                            cwd=ROOT,
                            capture_output=True,
                            check=True,
                        )
                        elapsed = time.perf_counter() - start
                        digest = hashlib.sha256(
                            output.read_bytes() + completed.stdout
                        ).hexdigest()
                        if hashes and digest != next(iter(hashes.values())):
                            raise RuntimeError(
                                f"Output drift: {case['target']} {workload} {label}"
                            )
                        hashes[label] = digest
                        if repeat >= 0:
                            times[label].append(elapsed)
                before, after = (statistics.median(times[k]) for k in engines)
                row = {
                    "target": case["target"],
                    "ligand": case["pdbqt"],
                    "workload": workload,
                    "arguments": common,
                    "seconds": times,
                    "before_median": before,
                    "after_median": after,
                    "speedup": before / after,
                    "output_sha256": hashes["before"],
                }
                report["results"].append(row)
                args.report.parent.mkdir(parents=True, exist_ok=True)
                args.report.write_text(json.dumps(report, indent=2) + "\n")
                print(
                    f"{case['target']:6} {workload:16} {before:.3f}s -> {after:.3f}s "
                    f"({before / after:.2f}x), identical output",
                    flush=True,
                )


if __name__ == "__main__":
    main()
