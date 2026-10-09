# RustDock Vina

A native Rust AutoDock Vina engine with a CLI, Python API and optional Apple Metal GPU search. Supports Vina, Vinardo and AutoDock4 scoring; final refinement runs in Rust.

## Benchmark results

**Metal was 1.36× faster than official Vina and 1.40× faster than Rust CPU** on the completed Apple M3 benchmark. Rust CPU was 2.7% slower than official Vina overall. Metal was faster in 94/100 cases.

100 ligand cases across five targets, 300 runs. Eight CPU threads, exhaustiveness 8, one pose and one trial per backend/case. Timings include startup, map generation, search and refinement.

| Backend | Total time (s) | Mean / case (s) | Speedup vs Vina | Mean absolute affinity difference (kcal/mol) | Median pose RMSD (Å) |
|---|---:|---:|---:|---:|---:|
| Official Vina 1.2.7 | 814.28 | 8.14 | 1.00× | — | — |
| Rust CPU | 835.95 | 8.36 | 0.97× | 0.062 | 0.118 |
| Rust Metal | 598.87 | 5.99 | 1.36× | 0.104 | 0.178 |

Affinity and pose differences are relative to official Vina's best prediction, not experimental accuracy. This is one run on one Mac; speedup is the ratio of total times.

![Mean docking time by target](docs/benchmarks/2026-10-09/runtime.png)

![Rust affinity estimates compared with official Vina](docs/benchmarks/2026-10-09/affinity.png)

[Full report and pose-RMSD graph](docs/benchmarks/2026-10-09/REPORT.md) · [300 measurements](docs/benchmarks/2026-10-09/results.csv) · [Protocol and provenance](docs/benchmarks/2026-10-09/summary.json)

## Build and dock

```sh
cargo build --release --workspace
./target/release/vina \
  --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --config box.txt --out poses.pdbqt --metal on
```

Use `--metal off` for CPU (default), `on` to require Metal, or `auto` for CPU fallback. Metal needs macOS 13+, Apple Command Line Tools, a rigid receptor, default Vina scoring, and a ligand with at most 64 heavy atoms and eight torsions. Other supported docking workflows use CPU.

## Run your own benchmark

All 100 included cases:

```sh
./tools/run-benchmark-100.sh \
  --reference-vina /path/to/official/vina \
  --cpu 8 --exhaustiveness 8 --num_modes 1 --seed 20260717 --trials 1
```

One receptor–ligand pair:

```sh
./target/release/vina-benchmark \
  --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
  --reference-vina /path/to/official/vina \
  --backends vina,off,on --trials 5 --benchmark-dir benchmark-results
```

The Rust benchmark tool reports timing, affinity differences and pose RMSD, and saves CSV/JSON, poses and logs. The shell script creates a new results folder per run. Use `--backends off,on` to compare only Rust CPU and Metal.

[Python API, experimental controls and validation](docs/USAGE.md) · [Benchmark instructions](benchmark-100/README.md) · [Port status and limits](docs/porting/COMPLETE_PORT_STATUS.md)
