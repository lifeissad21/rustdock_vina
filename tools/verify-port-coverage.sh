#!/usr/bin/env sh
set -eu

cargo test -p rustdock-vina-core every_reference_source_file_has_a_rust_counterpart
