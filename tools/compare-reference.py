#!/usr/bin/env python3
"""Run native Rust and official Vina on preserved examples; fail on score/local parity drift."""

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "reference" / "example"


def run(binary, args):
    result = subprocess.run(
        [str(binary), *map(str, args)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    output = result.stdout + result.stderr
    match = re.search(
        r"(?:Affinity|Estimated Free Energy of Binding)\s*:\s*([-+\d.]+)", output
    )
    if match:
        return float(match.group(1))
    poses = Path(args[args.index("--out") + 1]).read_text()
    return float(re.search(r"REMARK VINA RESULT:\s*([-+\d.]+)", poses).group(1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference",
        type=Path,
        default=ROOT.parent / "vina-multicore-benchmark" / "bin" / "vina",
    )
    parser.add_argument(
        "--rust", type=Path, default=ROOT / "target" / "release" / "vina"
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "docs" / "porting" / "regression-results.json",
    )
    parser.add_argument(
        "--full-docking",
        action="store_true",
        help="Also run exhaustiveness=8 docking; scores may diverge along stochastic search paths",
    )
    options = parser.parse_args()
    basic = EXAMPLES / "basic_docking" / "solution"
    rigid = [
        "--receptor",
        basic / "1iep_receptor.pdbqt",
        "--ligand",
        basic / "1iep_ligand.pdbqt",
        "--config",
        basic / "1iep_receptor.box.txt",
    ]
    flex = EXAMPLES / "flexible_docking" / "solution"
    flexible = [
        "--receptor",
        flex / "1fpu_receptor_rigid.pdbqt",
        "--flex",
        flex / "1fpu_receptor_flex.pdbqt",
        "--ligand",
        flex / "1iep_ligand.pdbqt",
        "--config",
        flex / "1fpu_receptor.box.txt",
    ]
    macro = EXAMPLES / "docking_with_macrocycles" / "solution"
    cases = [
        ("vina_score", rigid + ["--score_only"], 0.001),
        ("vina_grid_score", rigid + ["--score_only", "--no_refine"], 0.001),
        ("vinardo_score", rigid + ["--scoring", "vinardo", "--score_only"], 0.001),
        (
            "ad4_score",
            [
                "--ligand",
                basic / "1iep_ligand.pdbqt",
                "--maps",
                basic / "1iep_receptor",
                "--scoring",
                "ad4",
                "--score_only",
            ],
            0.001,
        ),
        ("flexible_score", flexible + ["--score_only"], 0.001),
        (
            "macrocycle_score",
            [
                "--receptor",
                macro / "BACE_1_receptor.pdbqt",
                "--ligand",
                macro / "BACE_1_ligand.pdbqt",
                "--config",
                macro / "BACE_1_receptor_vina_box.txt",
                "--score_only",
            ],
            0.001,
        ),
        ("vina_local", rigid + ["--local_only"], 0.001),
        ("vinardo_local", rigid + ["--scoring", "vinardo", "--local_only"], 0.001),
        ("flexible_local", flexible + ["--local_only"], 0.001),
        (
            "seeded_docking",
            rigid
            + ["--exhaustiveness", "2", "--num_modes", "3", "--max_evals", "20000"],
            None,
        ),
    ]
    if options.full_docking:
        cases.append(
            (
                "full_docking",
                rigid + ["--exhaustiveness", "8", "--num_modes", "9"],
                None,
            )
        )
    results = []
    with tempfile.TemporaryDirectory(prefix="rustdock-regression-") as tmp:
        for name, args, tolerance in cases:
            common = args + ["--cpu", "2", "--seed", "42", "--verbosity", "0"]
            scores = {}
            for label, binary in [
                ("official", options.reference),
                ("rust", options.rust),
            ]:
                out = Path(tmp) / f"{name}-{label}.pdbqt"
                scores[label] = run(binary.resolve(), common + ["--out", out])
            delta = abs(scores["official"] - scores["rust"])
            passed = tolerance is None or delta <= tolerance + 1e-12
            results.append(
                dict(
                    case=name,
                    **scores,
                    absolute_delta=delta,
                    tolerance=tolerance,
                    passed=passed,
                )
            )
            print(
                f"{name}: official={scores['official']:.3f}, Rust={scores['rust']:.3f}, delta={delta:.3f}",
                flush=True,
            )
    report = dict(
        reference_version=subprocess.check_output(
            [str(options.reference.resolve()), "--version"], text=True
        ).strip(),
        seed=42,
        cpu=2,
        results=results,
    )
    options.report.parent.mkdir(parents=True, exist_ok=True)
    options.report.write_text(json.dumps(report, indent=2) + "\n")
    if not all(case["passed"] for case in results):
        raise SystemExit(
            "Reference parity regression failed; inspect " + str(options.report)
        )


if __name__ == "__main__":
    main()
