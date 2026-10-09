#!/usr/bin/env python3
"""Publish the recorded 2026-10-09 Apple M3 benchmark snapshot and figures."""

import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODES = ["vina", "off", "on"]
LABELS = {"vina": "Official Vina 1.2.7", "off": "Rust CPU", "on": "Rust Metal"}
COLORS = {"vina": "#64748b", "off": "#2563eb", "on": "#0d9488"}


def first_pose(path):
    score, atoms = None, {}
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("REMARK VINA RESULT:"):
            score = float(line.split(":", 1)[1].split()[0])
        if line.startswith(("ATOM  ", "HETATM")) and line.split()[-1] not in (
            "H",
            "HD",
            "HS",
        ):
            serial = int(line[6:11])
            assert serial not in atoms, f"Duplicate serial: {path}"
            atoms[serial] = [float(line[i : i + 8]) for i in (30, 38, 46)]
    assert score is not None and math.isfinite(score) and atoms, f"Invalid pose: {path}"
    return score, atoms


def validate(run):
    manifest = json.loads((ROOT / "benchmark-100/manifest.json").read_text())
    cases = manifest["cases"]
    assert len(cases) == 100
    summaries = sorted(run.glob("case-*/summary.json"))
    assert len(summaries) == 100, "Run must contain all 100 cases"
    rows, digests = [], {}
    for i, (case, path) in enumerate(zip(cases, summaries)):
        assert path.parent.name == f"case-{i + 1:03d}"
        raw = path.read_bytes()
        digests[path.parent.name] = hashlib.sha256(raw).hexdigest()
        data = json.loads(raw)
        settings = data["settings"]
        for key, expected in {
            "trials": "1",
            "cpu": "8",
            "exhaustiveness": "8",
            "num_modes": "1",
            "seed": str(20260717 + i),
            "ligand": case["pdbqt"],
            "receptor": case["receptor"],
        }.items():
            assert settings[key] == [expected], f"Unexpected {key}: {path}"
        assert not set(settings).intersection(
            {
                "lanes",
                "steps",
                "local-steps",
                "max_evals",
                "no_refine",
                "maps",
                "scoring",
            }
        ), f"Unexpected protocol override: {path}"
        assert data["reference"]["version"] == "AutoDock Vina v1.2.7"
        runs = data["runs"]
        assert len(runs) == 3 and [r["requested_backend"] for r in runs] == MODES
        reference_time = runs[0]["seconds"]
        reference_score, reference_atoms = first_pose(
            path.parent / "trial-001-vina.pdbqt"
        )
        for record in runs:
            mode = record["requested_backend"]
            assert (
                record["actual_backend"] == mode
                and record["trial"] == 1
                and record["seed"] == 20260717 + i
            )
            assert math.isfinite(record["seconds"]) and record["seconds"] > 0
            score, atoms = first_pose(path.parent / f"trial-001-{mode}.pdbqt")
            assert math.isclose(score, record["best_affinity"], abs_tol=1e-9)
            assert atoms.keys() == reference_atoms.keys()
            rmsd = math.sqrt(
                sum(
                    sum((xyz[j] - reference_atoms[k][j]) ** 2 for j in range(3))
                    for k, xyz in atoms.items()
                )
                / len(atoms)
            )
            delta, speedup = score - reference_score, reference_time / record["seconds"]
            if mode != "vina":
                for key, calculated in [
                    ("rmsd", rmsd),
                    ("score_delta", delta),
                    ("speedup_vs_baseline", speedup),
                ]:
                    assert math.isclose(
                        record[key], calculated, rel_tol=1e-9, abs_tol=1e-9
                    ), f"Metric mismatch: {path} {key}"
            args = record["command_args"]
            assert "--maps" not in args
            assert args[args.index("--receptor") + 1] == case["receptor"]
            assert args[args.index("--ligand") + 1] == case["pdbqt"]
            log = (path.parent / f"trial-001-{mode}.log").read_text()
            assert "Vina error:" not in log and "using CPU:" not in log
            if mode == "on":
                assert "Metal: Apple M3" in log, f"Unexpected GPU: {path}"
            rows.append(
                {
                    "case": i + 1,
                    "target": case["target"],
                    "ligand": case["ligand"],
                    "heavy_atoms": case["heavy_atoms"],
                    "torsions": case["rotatable_bonds"],
                    "backend": mode,
                    "seed": record["seed"],
                    "seconds": record["seconds"],
                    "affinity_kcal_mol": score,
                    "delta_vs_vina_kcal_mol": delta,
                    "rmsd_vs_vina_angstrom": rmsd,
                    "speedup_vs_vina": speedup,
                }
            )
    return manifest, rows, digests


def metrics(rows):
    result = {}
    reference_total = sum(r["seconds"] for r in rows if r["backend"] == "vina")
    for mode in MODES:
        rr = [r for r in rows if r["backend"] == mode]
        total = sum(r["seconds"] for r in rr)
        result[mode] = {
            "runs": len(rr),
            "total_seconds": total,
            "mean_seconds": total / len(rr),
            "median_seconds": statistics.median(r["seconds"] for r in rr),
            "aggregate_speedup_vs_vina": reference_total / total,
            "mean_affinity_kcal_mol": statistics.mean(
                r["affinity_kcal_mol"] for r in rr
            ),
            "mean_delta_kcal_mol": statistics.mean(
                r["delta_vs_vina_kcal_mol"] for r in rr
            ),
            "mean_absolute_delta_kcal_mol": statistics.mean(
                abs(r["delta_vs_vina_kcal_mol"]) for r in rr
            ),
            "max_absolute_delta_kcal_mol": max(
                abs(r["delta_vs_vina_kcal_mol"]) for r in rr
            ),
            "median_rmsd_angstrom": statistics.median(
                r["rmsd_vs_vina_angstrom"] for r in rr
            ),
            "max_rmsd_angstrom": max(r["rmsd_vs_vina_angstrom"] for r in rr),
            "cases_faster_than_vina": sum(r["speedup_vs_vina"] > 1 for r in rr),
        }
    return result


def figures(rows, targets, out):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.titleweight": "bold",
            "svg.fonttype": "none",
            "svg.hashsalt": "rustdock-vina-benchmark",
        }
    )

    def save(fig, name):
        fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
        svg = out / f"{name}.svg"
        fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
        svg.write_text(
            "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n"
        )
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.2), layout="constrained")
    x = np.arange(len(targets))
    for i, mode in enumerate(MODES):
        means = [
            statistics.mean(
                r["seconds"]
                for r in rows
                if r["backend"] == mode and r["target"] == target
            )
            for target in targets
        ]
        bars = ax.bar(
            x + (i - 1) * 0.25,
            means,
            width=0.25,
            label=LABELS[mode],
            color=COLORS[mode],
        )
        ax.bar_label(bars, fmt="%.2f", fontsize=9, padding=3)
    ax.set_xticks(x, [t.upper() for t in targets])
    ax.set_ylabel("Mean wall time per docking (seconds)")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
    ax.set_title("Docking runtime across five targets", loc="left", pad=18)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)
    fig.supxlabel(
        "Apple M3 · 20 cases per target · 8 CPU threads · exhaustiveness 8", fontsize=10
    )
    save(fig, "runtime")
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 5.2), layout="constrained", sharex=True, sharey=True
    )
    target_colors = dict(
        zip(targets, ["#2563eb", "#0d9488", "#9333ea", "#ea580c", "#be185d"])
    )
    bounds = [r["affinity_kcal_mol"] for r in rows]
    low, high = min(bounds) - 0.4, max(bounds) + 0.4
    reference = {
        r["case"]: r["affinity_kcal_mol"] for r in rows if r["backend"] == "vina"
    }
    for ax, mode in zip(axes, ["off", "on"]):
        for target in targets:
            rr = [r for r in rows if r["backend"] == mode and r["target"] == target]
            ax.scatter(
                [reference[r["case"]] for r in rr],
                [r["affinity_kcal_mol"] for r in rr],
                s=23,
                alpha=0.8,
                color=target_colors[target],
                label=target.upper(),
            )
        ax.plot([low, high], [low, high], "--", color="#475569", lw=1)
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.set_aspect("equal")
        mae = statistics.mean(
            abs(r["delta_vs_vina_kcal_mol"]) for r in rows if r["backend"] == mode
        )
        ax.set_title(
            f"{LABELS[mode]}\nMean absolute difference: {mae:.3f} kcal/mol", loc="left"
        )
        ax.set_xlabel("Official Vina affinity (kcal/mol)")
        ax.grid(alpha=0.12)
    axes[0].set_ylabel("Rust affinity (kcal/mol)")
    axes[1].legend(frameon=False, fontsize=9, loc="lower right")
    fig.suptitle(
        "Best-pose affinity estimates",
        fontsize=16,
        fontweight="bold",
        ha="left",
        x=0.02,
    )
    save(fig, "affinity")
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for mode in ["off", "on"]:
        values = sorted(
            r["rmsd_vs_vina_angstrom"] for r in rows if r["backend"] == mode
        )
        ax.step(
            values,
            np.arange(1, len(values) + 1) / len(values) * 100,
            where="post",
            lw=2,
            color=COLORS[mode],
            label=f"{LABELS[mode]} · median {statistics.median(values):.3f} Å",
        )
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xticks([0, 0.1, 1, 10], labels=["0", "0.1", "1", "10"])
    ax.set_ylim(0, 102)
    ax.set_xlabel(
        "Heavy-atom RMSD from official Vina best pose (Å; scale changes above 0.1)"
    )
    ax.set_ylabel("Cases at or below this RMSD (%)")
    ax.set_title("Predicted-pose differences", loc="left", pad=15)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.15)
    fig.supxlabel(
        "Direct atom-serial comparison · no alignment or symmetry correction · not crystal-pose accuracy",
        fontsize=9,
    )
    save(fig, "pose-rmsd")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: uv run --no-project --with matplotlib python tools/publish-benchmark.py RUN_DIR OUTPUT_DIR"
        )
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    manifest, rows, digests = validate(run)
    overall = metrics(rows)
    by_target = {
        target: metrics([r for r in rows if r["target"] == target])
        for target in manifest["targets"]
    }
    out.mkdir(parents=True, exist_ok=True)
    with (out / "results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "source_run": run.name,
        "benchmark_version": manifest["benchmark_version"],
        "run_date": "2026-10-09",
        "gpu_from_run_logs": "Apple M3",
        "host_at_analysis": {"chip": "Apple M3", "memory_gib": 8, "macos": "26.7.1"},
        "source_commit_at_analysis": "7d706e5",
        "protocol": {
            "cases": 100,
            "targets": 5,
            "trials_per_case_per_backend": 1,
            "cpu_threads": 8,
            "exhaustiveness": 8,
            "num_modes": 1,
            "first_seed": 20260717,
            "last_seed": 20260816,
            "reference_version": "AutoDock Vina v1.2.7",
            "metal_controls": "automatic",
            "maps": "generated per run",
            "timing": "wall time including startup, map generation, search and backend refinement",
            "backend_order": MODES,
        },
        "overall": overall,
        "by_target": by_target,
        "source_summary_sha256": digests,
        "limitations": [
            "One trial per case; no repeated-run confidence intervals.",
            "Backend order was fixed, not randomized.",
            "RMSD is against official predictions, not crystal poses; no fitting or symmetry correction.",
            "Score differences mix scoring and stochastic pose-search differences.",
            "Hardware state and binary hashes were not recorded by the runner at execution time; host metadata and repository HEAD were collected at analysis, not verified as build provenance.",
            "These measurements do not establish universal speedup or full scientific parity.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures(rows, manifest["targets"], out)
    print(json.dumps(overall, indent=2))


if __name__ == "__main__":
    main()
