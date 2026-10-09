#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
benchmark_arguments=(--manifest benchmarks/datasets/dude-100/manifest.json)
benchmark_has_dir=false
benchmark_has_trials=false
for argument in "$@"; do
  case "$argument" in
    --benchmark-dir|--benchmark-dir=*) benchmark_has_dir=true ;;
    --trials|--trials=*) benchmark_has_trials=true ;;
  esac
done
if ! "$benchmark_has_dir"; then
  mkdir -p benchmarks/runs/100-case
  benchmark_arguments+=(--benchmark-dir "benchmarks/runs/100-case/run-$(date +%Y%m%d-%H%M%S)-$$")
fi
if ! "$benchmark_has_trials"; then
  benchmark_arguments+=(--trials 1)
fi
cargo build --release --locked --bin vina --bin vina-benchmark
exec ./target/release/vina-benchmark "${benchmark_arguments[@]}" "$@"
