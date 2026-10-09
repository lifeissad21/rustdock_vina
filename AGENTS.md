# AI Porting Guidance

Before continuing any Rust porting work, read:

- `docs/porting/COMPLETE_PORT_STATUS.md`
- `docs/porting/README.md`
- `crates/rustdock-vina-core/src/porting.rs`

Rules for this rewrite:

1. Do not start optimization work until every `reference/src` entry is marked `Ported` in `SOURCE_PORTS` and has parity tests.
2. Every original source/header under `reference/src` must keep a Rust counterpart.
3. When a scaffold is replaced by real Rust code, update `SOURCE_PORTS` from `Scaffolded` to `Ported`.
4. Add or update tests for each ported module before moving on.
5. Run `cargo fmt --all --check` and `cargo test --workspace` before handing off.
6. Preserve behavior first. Improve performance only after the full port and regression harness are complete.
