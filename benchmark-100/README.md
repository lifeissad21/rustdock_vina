# Rust CPU / Metal benchmarks

Both single-pair trials and the 100-case runner use the Rust `vina-benchmark` command in its own `rustdock-vina-benchmark` crate. No Bun or TypeScript is required.

```sh
cargo build --release --locked --bin vina --bin vina-benchmark
./target/release/vina-benchmark --help
```

For your own pair, pass regular CLI inputs (`--receptor`, `--ligand`, `--config` or box coordinates, optional `--maps`) and search settings. `--trials 5 --backends vina,off,on --reference-vina /path/to/official/vina` runs five paired official Vina / Rust CPU / Rust Metal trials. Defaults are `vina,off,on`, using `../vina-multicore-benchmark/bin/vina` when no reference path is supplied. Use `--backends off,on` for Rust-only comparisons. `--metal off` selects CPU alone, `--metal on` selects required Metal alone, and `--backends off,on,auto` includes automatic fallback. CPU and Metal off are the same backend. GPU controls apply only to on. Default trials: three; first seed: 20260717, incremented each trial with matching seeds across backends.

The included dataset can be run from the repository root:

```sh
./tools/run-benchmark-100.sh --check
./tools/run-benchmark-100.sh --limit 3 --benchmark-dir benchmark-100/trial-run
./tools/run-benchmark-100.sh
./tools/run-benchmark-100.sh --trials 5 --benchmark-dir benchmark-100/repeated-run
```

The script builds the release binary and selects the manifest, one trial per case/backend, and a fresh dated directory under `benchmark-100/results/` as its default output directory. The shell script preserves previous results and creates a new directory on each invocation. Explicit output directories must be new; existing directories are rejected. Run from the repository root when supplying the manifest directly: its paths are relative to the working directory. `--check` checks input presence and box options, not complete map contents or hardware compatibility. The full benchmark has not been run during setup.

Outputs per pair: `results.csv`, `summary.json`, and `trial-NNN-off/on/auto.pdbqt` and `.log` files. The 100-case runner places each pair in `case-001` through `case-100` subdirectories. The terminal prints each pair's per-trial measurements and backend averages. Reports are written after each successful run; a failure stops with its log preserved. There is no automatic resume or overwrite option.

Wall time includes process startup, map loading or generation, GPU compilation when applicable, search and each backend’s final refinement. Affinity delta is backend minus official Vina (or Rust CPU when official Vina is omitted) in kcal/mol. Pose RMSD compares heavy atoms by serial in the receptor coordinate frame, with no fitting or symmetry correction; it is not crystal-pose accuracy. Different stochastic trajectories can produce different minima. Experimental GPU settings affect runtime and convergence without guaranteeing better poses. Auto results can include CPU fallback; inspect actual_backend in saved records.

The manifest and 155 prepared receptor, ligand, box and affinity-map files were copied from the sibling multicore benchmark: DUD-E Diverse, version dude-diverse-100-v1. AMPC, CXCR4, GCR, HIVPR and HIVRT each have 20 cases, at most 38 heavy atoms and eight torsions. Inputs are included here; the sibling directory is not needed. Preparation has not been independently audited. This Rust protocol differs from the sibling stock-Vina/Swift timing protocol.

CPU-only benchmarks work without Metal. Metal requires a compatible Mac, macOS 13+, and Apple Command Line Tools when building the helper.

Official Vina's CLI cannot combine receptor and maps. For reference comparisons with a receptor, the runner ignores precomputed maps for all backends and generates maps from the same box each run. All timings therefore include map generation and final refinement. Reference version/path, per-run binary/arguments, paired speedup, signed affinity differences and direct pose RMSD are saved. Different minima can give different estimates; this is not crystal-pose accuracy or proof of identical scoring. The reference executable is checked with `--version` and must differ from the Rust binary.
