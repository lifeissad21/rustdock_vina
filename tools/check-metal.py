#!/usr/bin/env python3
"""Check real Metal search, reproducibility and CPU fallback; record wall times."""

import argparse
import json
import re
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustiveness", type=int, default=1)
    args = parser.parse_args()
    basic = ROOT / "reference/example/basic_docking/solution"
    common = [
        str(ROOT / "target/release/vina"),
        "--receptor",
        str(basic / "1iep_receptor.pdbqt"),
        "--ligand",
        str(basic / "1iep_ligand.pdbqt"),
        "--config",
        str(basic / "1iep_receptor.box.txt"),
        "--seed",
        "42",
        "--cpu",
        "8",
        "--num_modes",
        "3",
        "--exhaustiveness",
        str(args.exhaustiveness),
    ]
    results = {}
    with tempfile.TemporaryDirectory(prefix="rustdock-metal-check-") as tmp:

        def run(name, mode, extra=()):
            out = Path(tmp) / (name + ".pdbqt")
            start = time.monotonic()
            process = subprocess.run(
                common + ["--metal", mode, "--out", str(out), *extra],
                capture_output=True,
                text=True,
                check=True,
            )
            elapsed = time.monotonic() - start
            text = out.read_text()
            score = float(re.search(r"REMARK VINA RESULT:\s*([-+\d.]+)", text).group(1))
            assert text.count("MODEL ") >= 1 and text.count("ENDMDL") == text.count(
                "MODEL "
            )
            results[name] = {
                "wall_seconds": elapsed,
                "best_affinity": score,
                "stdout": process.stdout,
                "stderr": process.stderr,
            }
            print(f"{name}: {elapsed:.3f} s, {score:.3f} kcal/mol", flush=True)
            return text

        run("cpu", "off")
        first = run("metal", "on")
        assert "Metal:" in results["metal"]["stdout"]
        assert run("metal_repeat", "on") == first, "Same-seed Metal output changed"
        off = run("limited_cpu", "off", ["--max_evals", "2000"])
        assert run("auto_fallback", "auto", ["--max_evals", "2000"]) == off
        assert "using CPU" in results["auto_fallback"]["stderr"]
    report = {
        "exhaustiveness": args.exhaustiveness,
        "seed": 42,
        "cpu": 8,
        "results": results,
    }
    (ROOT / "benchmarks/results/validation/metal-results.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
