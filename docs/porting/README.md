# RustDock Vina Porting Notes

This rewrite is coverage-first.

Before optimization work, every C++ source/header under `reference/src` must have a Rust counterpart in the workspace and must be listed in `rustdock_vina_core::porting::SOURCE_PORTS`.

Status meanings:

- `Ported`: an initial Rust implementation exists and is compiled/tested.
- `Scaffolded`: a deliberate placeholder exists so the original file cannot be lost during the port.

Non-source assets under `reference/docs`, `reference/example`, `reference/data`, and `reference/build/python` remain preserved fixtures until their replacement surfaces are implemented.
