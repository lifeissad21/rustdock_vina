# Performance review — 2026-10-09

## What the project does

`rustdock-vina-core` implements PDBQT parsing, receptor/ligand topology, Vina,
Vinardo and AD4 scoring, affinity maps, BFGS local optimization, seeded Monte
Carlo search and final pose refinement. CPU searches share immutable grids and
precalculated pair tables across scoped workers; each task owns its mutable
model. The CLI exposes these operations, the Python frontend uses a persistent
JSON-lines server, and the benchmark crate launches isolated CLI processes.
Metal is an optional candidate-search backend with native Rust refinement.

The port registry lists every upstream source/header as `Ported`. Source
coverage and the workspace suite pass, with official score/local golden tests
for the main workflows. This remains fixture-based parity, not universal
scientific validation or identical stochastic convergence.

## Changes tested

- Reuse the affinity accumulation buffer across grid points, clearing it before
  each point. This removes one heap allocation per point.
- Traverse generated maps with X as the inner loop, matching `Array3d` storage.
  Each point retains the original receptor-atom summation order.

ThinLTO with one codegen unit was also tested, then **rejected**. It reduced
preparation/limited-docking medians by roughly 2–5%, but uncapped docking showed
slower medians on some targets (about 10% on CXCR4 in the extended comparison).
The existing release profile is preserved. Its experimental measurements are
saved separately in `performance-experimental-lto-limited-results.json` and
`performance-experimental-lto-full-results.json`; the latter is partial because
the repeat was stopped after detecting the regression. Host-load variation
limits causal interpretation; these data do not establish that ThinLTO always
causes a regression. There was insufficient benefit to keep it.

## Measurement protocol

`benchmarks/scripts/compare-performance.py` compares preserved release binaries on the first
ligand from each of AMPC, CXCR4, GCR, HIVPR and HIVRT. Every run generates maps
from the included receptor and box. CPU=1, seed=42 and Metal=off are fixed.
Preparation runs use `--randomize_only` because some supplied ligand coordinates
are outside the box. This includes parsing, precalculation, map generation,
randomization and output, rather than timing map generation alone.

Limited docking uses exhaustiveness=1, max_evals=2000 and three output modes.
The uncapped comparison removes max_evals, retaining exhaustiveness=1. Timing
includes startup and final refinement. Each workload warms both binaries once,
then alternates execution order. Validation/build processes are stopped before
measuring. Output files and captured stdout must match byte-for-byte across
both binaries and all repetitions; a mismatch stops the comparison.

The final comparison uses seven measured repeats per binary/workload, plus
warmups, with uncapped docking. Raw per-run seconds, medians, arguments, platform
and executable SHA-256 hashes are recorded in `performance-results.json`.
The benchmark script was also checked with the same baseline executable on both
sides across all five targets; all equality checks passed.

To repeat, save the original release executable before applying the changes,
build the updated executable, then run:

```sh
python3 benchmarks/scripts/compare-performance.py \
  --before /tmp/rustdock-performance/vina-before \
  --after target/release/vina --repeats 7 --full-docking \
  --report /tmp/performance-results.json
```

## Final results

Seven measured repeats per binary/workload, one ligand per target, on this Mac:

| Target | Preparation before → after (s) | Uncapped docking before → after (s) |
|---|---:|---:|
| AMPC | 0.461 → 0.460 | 1.377 → 1.382 |
| CXCR4 | 0.629 → 0.625 | 3.983 → 3.975 |
| GCR | 0.714 → 0.703 | 1.855 → 1.855 |
| HIVPR | 0.634 → 0.628 | 3.328 → 3.338 |
| HIVRT | 0.781 → 0.774 | 3.705 → 3.654 |

The sums of the five case medians declined by **0.89% for preparation** and
**0.30% for uncapped docking**. These small differences, with some slower case
medians and observed run-to-run variability, do **not** establish a material
runtime improvement. This is a negative result for substantial speedup, not a
claim that the engine is now generally faster.

The retained map change removes a known allocation cost without altering
numerics: the AMPC box at spacing 0.375 has 52 × 50 × 61 = 158,600 points, so
scratch allocations drop from 158,600 to one for that map-population call.
The new grid test compares both generated atom-type maps at every point against
direct receptor-pair sums and verifies repeat population leaves maps unchanged.
All warmup and measured preparation/docking outputs matched byte-for-byte.

## Verification

- `cargo fmt --all --check` passed.
- `cargo test --workspace` passed: 108 tests across 12 suites.
- The final release workspace build passed.
- The official comparison harness passed all nine Vina score/local comparisons,
  with zero differences at the reported three-decimal precision. Its existing
  limited seeded docking difference remains: official -10.819 versus Rust
  -8.565 kcal/mol. That observation is not an equality assertion.
- The benchmark script compiles and its same-executable self-check passed.

## Further opportunities

1. **Reuse receptor maps and scoring tables across batch ligands.** The CLI
   constructs a new engine and prepares maps for each batch entry. For a fixed
   receptor, scoring weights, box and spacing, reuse could amortize preparation
   across many ligands. Handle newly required atom types and preserve per-ligand
   seed behavior. Measure a real multi-ligand batch before changing this flow.
2. **Profile BFGS allocation and gradient access.** BFGS clones nested
   configuration/gradient vectors repeatedly, and its matrix operations use
   flattened accessors that walk ligand/residue structures. Reusable buffers
   or direct iteration may help long CPU searches. Preserve update order,
   rollback behavior and evaluation counts; measure search separately first.
3. **Buffer map-file output.** Map writing currently calls `writeln!` on a
   `File` for every grid value. `BufWriter` is a small potential win for
   `--write_maps`, with flush errors propagated. It would not speed normal
   docking that does not export maps.

These are code-inspection candidates, not measured speedup promises. The retained
change leaves search budgets, precision, scoring equations, dependencies and
user-visible options unchanged.
Metal hardware, Python workflows, Linux/Windows, higher-exhaustiveness runs and
all 100 cases require separate validation before generalizing these timings.
