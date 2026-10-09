# 100-case benchmark — 2026-10-09

Complete user-run dataset: `run-20261009-125526-42917`. All 100 case summaries contain one successful official Vina, Rust CPU and Rust Metal run: 300 measurements. No partial or failed earlier run was included. The GPU logs identify Apple M3; this host reported 8 GiB memory and macOS 26.7.1 at analysis. The source checkout at analysis was commit `7d706e5`; the runner did not capture execution-time binary hashes or system load, so the build revision cannot be independently verified from these records.

## Overall results

| Backend | Total time (s) | Mean / case (s) | Speedup vs Vina | Mean absolute affinity difference (kcal/mol) | Median pose RMSD (Å) |
|---|---:|---:|---:|---:|---:|
| Official Vina 1.2.7 | 814.28 | 8.14 | 1.00× | — | — |
| Rust CPU | 835.95 | 8.36 | 0.97× | 0.062 | 0.118 |
| Rust Metal | 598.87 | 5.99 | 1.36× | 0.104 | 0.178 |

Metal / Rust CPU aggregate speedup: 1.40×. Metal was faster than official Vina in 94/100 cases; Rust CPU in 30/100.

## Results by target

| Target | Cases | Official mean (s) | Rust CPU mean (s) | Metal mean (s) | Metal speedup vs Vina | CPU mean absolute affinity difference | Metal mean absolute affinity difference |
|---|---:|---:|---:|---:|---:|---:|---:|
| AMPC | 20 | 3.403 | 3.626 | 2.610 | 1.30× | 0.025 | 0.043 |
| CXCR4 | 20 | 9.154 | 10.042 | 6.629 | 1.38× | 0.104 | 0.124 |
| GCR | 20 | 5.373 | 5.452 | 3.746 | 1.43× | 0.033 | 0.023 |
| HIVPR | 20 | 13.853 | 13.351 | 10.975 | 1.26× | 0.072 | 0.161 |
| HIVRT | 20 | 8.931 | 9.327 | 5.984 | 1.49× | 0.077 | 0.167 |

Affinity differences are in kcal/mol. Global means weight all 100 cases equally. Target speedup uses each target's sum of 20 runtimes.

## Protocol and verification

- DUD-E Diverse prepared manifest `dude-diverse-100-v1`: five targets, 20 prepared ligand cases each. Cases can repeat compound identities with different conformers; they are not 100 independent unique compounds.
- Official AutoDock Vina 1.2.7, native Rust CPU search, and Metal GPU search with native Rust final refinement.
- Eight CPU threads, exhaustiveness 8, one returned pose, one trial per case/backend, no GPU overrides. Seed is 20260717 + case index, matched across the three backends.
- Sequential fixed order: official Vina, Rust CPU, Metal. All generate maps using the same receptor and search-box inputs. Wall time includes startup, map generation, search, GPU initialization/compilation where applicable, and each backend's refinement.
- Best affinity was re-read from all 300 PDBQT outputs. All 200 score differences, direct heavy-atom RMSDs and paired speedup ratios were recomputed and checked against the saved JSON values. Backend, seed, settings, inputs and GPU identity were checked for every case.
- Signed affinity difference is Rust minus official Vina. Lower affinity is more negative. The maximum absolute differences are 0.842 kcal/mol (CPU) and 0.866 kcal/mol (Metal).
- RMSD compares heavy atoms by serial in the receptor frame, without fitting or symmetry correction. Maximum RMSDs are 7.039 Å (CPU) and 7.828 Å (Metal).

## Figures

![Mean time per target](runtime.png)

![Best-pose affinity comparison](affinity.png)

![Direct RMSD cumulative distributions](pose-rmsd.png)

PNG images are used in the README; SVG versions are included for export. The RMSD graph uses a linear scale to 0.1 Å and a logarithmic scale above that value to show both small differences and outliers.

## Limits

These are observations on one Mac and one run per case, with no repeated-run confidence intervals. Fixed backend order, warm caches, background load and thermal state can affect timing. Hardware metadata beyond GPU identity was collected after the run. Matching seeds do not imply matching stochastic trajectories or search budgets. Comparing best-pose affinity estimates mixes search and scoring differences; official predictions are a reference, not experimental truth. Neither affinity agreement nor direct RMSD proves crystal-pose accuracy, complete parity, or a universal performance improvement.

## Reproduce

Run from the repository root with a compatible Mac and an official Vina binary:

```sh
./tools/run-benchmark-100.sh \
  --reference-vina /path/to/official/vina \
  --cpu 8 --exhaustiveness 8 --num_modes 1 --seed 20260717 --trials 1
```

The publication script validates this protocol and the recorded Apple M3 GPU logs, exports portable CSV/JSON, and regenerates the figures. Its snapshot metadata describes this publication; it is not a generic hardware inventory collector.

```sh
uv run --no-project --with matplotlib python tools/publish-benchmark.py \
  benchmark-100/results/run-20261009-125526-42917 docs/benchmarks/2026-10-09
```

The ignored raw source directory remains local. [CSV](results.csv) contains all 300 measurements. [JSON](summary.json) contains aggregate metrics, protocol, limitations and SHA-256 digests of the source summaries. Personal absolute paths and raw logs are excluded from these published artifacts.
