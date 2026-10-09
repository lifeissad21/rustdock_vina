# Port status

All `reference/src` sources map to compiled Rust implementations in `crates/rustdock-vina-core/src/porting.rs`; the coverage test verifies the mapping. Placeholder source files have been replaced. `Ported` means an initial implementation exists, not universal numerical equivalence.

## Implemented

- PDBQT parsing, ligand/flexible-residue trees, mobility, bond inference, atom typing, model append and coordinate output.
- Vina/Vinardo/AD4 scoring, map generation/read/write, local BFGS optimization, native threaded Monte Carlo docking, refinement and pose selection.
- MT19937 with Boost-compatible uniform and normal distributions, checked against C++ oracle values.
- CLI configuration, batch/multiple ligands, scoring/search actions, map/output options and `vina_split`.
- Compatible Python frontend backed by a persistent JSON-lines Rust server; platform wheel builds bundle the server.
- Optional macOS Metal search with `--metal on|off|auto`; binary Rust maps/topology, runtime calibration and native Rust final refinement. See root README for supported inputs.

## Evidence and limits

After integration, all 97 Rust tests and four Python workflow tests passed on this Mac. Formatting, release builds and CPU-only checks passed. A platform wheel was built, installed into a temporary directory and used for the Python tests. Clippy completes with style warnings; it is not warning-free.

`metal-results.json` records a real Apple M3 hardware check: exhaustiveness 8 CPU wall time 34.18 seconds versus Metal 12.74 seconds; best affinities -13.183 versus -13.299. Same-seed Metal poses repeat exactly on this device, and auto fallback preserves CPU output. Single-fixture timing is not a general performance or scientific-equivalence claim. Platform CI and binary artifacts are configured but have not run on Linux or Windows here.

`regression-results.json` records comparison with official Vina 1.2.7. Nine scoring/local examples (Vina, grid-only, Vinardo, AD4, flexible receptor and macrocycle) match within 0.001 kcal/mol. Seed 42 full 1IEP docking at exhaustiveness 8 produced official -13.237 and Rust -13.183. Limited docking at 20,000 evaluations produced -10.819 and -8.565 respectively. Docking comparisons are observations, not passing equality assertions: distribution matching does not establish identical stochastic trajectories or convergence.

The generic legacy `parallel.rs` FnMut helper remains sequential; production Monte Carlo uses scoped threads in `parallel_mc.rs`. Python uses CPU search. Exact Boost CLI help/error formatting is not replicated. Original documentation, examples and scientific assets remain preserved under `reference/`.

Linux/Windows runtime behavior and release archives require their platform CI runs; local checks on this Mac cannot establish those results. Scientific coverage is fixture-based and should be expanded before treating this rewrite as a validated replacement for every official Vina use case.

New-target 3PTB trypsin–benzamidine runs are recorded in `results/3ptb-benzamidine/`: five Metal seeds and one CPU run recovered best poses within 0.38 Angstrom crystal RMSD. This new translated box exposed an upper-bound reconstruction rounding failure in explicit refinement; the shared spatial-grid check now uses scale-adjusted machine-epsilon tolerance, with tests for valid boundary and invalid outside queries. The 96-test workspace suite and formatting check passed after the fix.

CLI presentation now groups settings and stage messages, labels pose-table units, reports saved outputs and elapsed time, and uses terminal-only colors respecting NO_COLOR. The 97-test workspace suite, formatting and release build passed; the new 1STP docking inputs remain unrun.

Experimental Metal controls `--lanes`, `--steps`, and `--local-steps` are exposed with validation and require `--metal on`. They control GPU search; final refinement remains native Rust. Calibration now checks map interpolation at the box center without boundary penalties, using separate synthetic coordinates that do not change ligand search geometry. Synthetic hardware checks passed against all five benchmark targets' maps.

The self-contained `benchmark-100/` dataset and `tools/run-benchmark-100.sh` compare this Rust binary's CPU search against Metal search with Rust refinement. Input checks verified 100 cases, five targets, and 155 prepared files. The benchmark itself has not been run. Its protocol differs from the original multicore benchmark, so timings are not directly interchangeable. Current verification passed: 101 Rust tests, three Bun runner tests, formatting, release build, shell syntax, and CPU-only compilation.

Benchmark orchestration is now native Rust via `vina-benchmark`, accepting standard docking inputs, paired seeded trials, off/on/auto backend selections, CSV/JSON measurements, pose files and logs. The 100-case shell wrapper invokes this Rust command; the TypeScript runner and tests were removed. Existing output directories are rejected; the new runner does not resume previous results. Full benchmark docking remains unrun.

The separate `rustdock-vina-benchmark` crate launches the `vina` executable, which calls `rustdock-vina-core`; shared CLI parsing is exported by `rustdock-vina-cli`. Verification passed: all 105 workspace tests, formatting, release builds of both executables, CPU-only compilation, shell syntax and all 100 manifest input checks. A synthetic subprocess test verifies paired seeds, actual auto fallback, score/RMSD comparisons, saved records and refusal to overwrite existing output directories. No new real molecular benchmark trials were run during this replacement.

Benchmark presentation now shares the engine's terminal output helper: grouped settings, progress, aligned trial/average tables with units, actual backend labels and saved paths. The shell wrapper creates a fresh dated run folder under `benchmark-100/results/`, preserves existing results, and accepts trial/directory overrides without duplicate arguments. All 106 workspace tests, formatting, release builds, shell syntax and wrapper input-only checks passed, including explicit overrides on macOS Bash. Full benchmark docking remains unrun.

The benchmark now defaults to official Vina / Rust CPU / Rust Metal, accepting `--reference-vina FILE` and recording the reference version/path. Paired per-trial speedup, signed affinity differences and direct heavy-atom pose RMSD are computed against official Vina, or Rust CPU when reference runs are omitted. Official CLI boolean switches and GPU-option filtering are tested. Reference comparisons with a receptor generate maps from the same box for every backend because official Vina rejects receptor plus precomputed maps; timing includes map generation and each backend's refinement. All 107 workspace tests, formatting, release builds, CPU-only compilation and a reference-enabled wrapper input check passed. The local reference reports AutoDock Vina v1.2.7. New molecular benchmark comparisons remain unrun; synthetic subprocess tests cover reference orchestration and metrics.

The user's completed 2026-10-09 100-case run is published under `docs/benchmarks/2026-10-09/`, with 300 portable measurements, summary/provenance, and runtime, affinity and pose-RMSD figures. Official Vina / Rust CPU / Metal total times were 814.284 / 835.953 / 598.871 seconds (Metal 1.36× versus official). All saved pose scores and pairwise metrics were independently checked. One trial per case, fixed backend order and stochastic differences limit interpretation; this supersedes earlier statements that the full benchmark remained unrun.
