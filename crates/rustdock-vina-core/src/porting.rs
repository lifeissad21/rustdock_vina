#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortStatus {
    Ported,
    Scaffolded,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourcePort {
    pub reference_path: &'static str,
    pub rust_path: &'static str,
    pub status: PortStatus,
}

pub const SOURCE_PORTS: &[SourcePort] = &[
    SourcePort {
        reference_path: "reference/src/lib/ad4cache.cpp",
        rust_path: "crates/rustdock-vina-core/src/ad4cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/ad4cache.h",
        rust_path: "crates/rustdock-vina-core/src/ad4cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/array3d.h",
        rust_path: "crates/rustdock-vina-core/src/array3d.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/atom.h",
        rust_path: "crates/rustdock-vina-core/src/atom.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/atom_base.h",
        rust_path: "crates/rustdock-vina-core/src/atom_base.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/atom_constants.h",
        rust_path: "crates/rustdock-vina-core/src/atom_constants.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/atom_type.h",
        rust_path: "crates/rustdock-vina-core/src/atom_type.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/bfgs.h",
        rust_path: "crates/rustdock-vina-core/src/bfgs.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/brick.h",
        rust_path: "crates/rustdock-vina-core/src/brick.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/cache.cpp",
        rust_path: "crates/rustdock-vina-core/src/cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/cache.h",
        rust_path: "crates/rustdock-vina-core/src/cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/common.h",
        rust_path: "crates/rustdock-vina-core/src/common.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/conf.h",
        rust_path: "crates/rustdock-vina-core/src/conf.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/conf_independent.cpp",
        rust_path: "crates/rustdock-vina-core/src/conf_independent.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/conf_independent.h",
        rust_path: "crates/rustdock-vina-core/src/conf_independent.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/convert_substring.h",
        rust_path: "crates/rustdock-vina-core/src/convert_substring.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/coords.cpp",
        rust_path: "crates/rustdock-vina-core/src/coords.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/coords.h",
        rust_path: "crates/rustdock-vina-core/src/coords.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/curl.h",
        rust_path: "crates/rustdock-vina-core/src/curl.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/file.h",
        rust_path: "crates/rustdock-vina-core/src/file.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/grid.cpp",
        rust_path: "crates/rustdock-vina-core/src/grid.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/grid.h",
        rust_path: "crates/rustdock-vina-core/src/grid.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/grid_dim.h",
        rust_path: "crates/rustdock-vina-core/src/grid_dim.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/igrid.h",
        rust_path: "crates/rustdock-vina-core/src/igrid.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/incrementable.h",
        rust_path: "crates/rustdock-vina-core/src/incrementable.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/int_pow.h",
        rust_path: "crates/rustdock-vina-core/src/int_pow.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/macros.h",
        rust_path: "crates/rustdock-vina-core/src/macros.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/matrix.h",
        rust_path: "crates/rustdock-vina-core/src/matrix.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/model.cpp",
        rust_path: "crates/rustdock-vina-core/src/model.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/model.h",
        rust_path: "crates/rustdock-vina-core/src/model.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/monte_carlo.cpp",
        rust_path: "crates/rustdock-vina-core/src/monte_carlo.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/monte_carlo.h",
        rust_path: "crates/rustdock-vina-core/src/monte_carlo.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/mutate.cpp",
        rust_path: "crates/rustdock-vina-core/src/mutate.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/mutate.h",
        rust_path: "crates/rustdock-vina-core/src/mutate.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/non_cache.cpp",
        rust_path: "crates/rustdock-vina-core/src/non_cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/non_cache.h",
        rust_path: "crates/rustdock-vina-core/src/non_cache.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parallel.h",
        rust_path: "crates/rustdock-vina-core/src/parallel.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parallel_mc.cpp",
        rust_path: "crates/rustdock-vina-core/src/parallel_mc.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parallel_mc.h",
        rust_path: "crates/rustdock-vina-core/src/parallel_mc.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parallel_progress.cpp",
        rust_path: "crates/rustdock-vina-core/src/parallel_progress.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parallel_progress.h",
        rust_path: "crates/rustdock-vina-core/src/parallel_progress.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parse_error.h",
        rust_path: "crates/rustdock-vina-core/src/parse_error.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parse_pdbqt.cpp",
        rust_path: "crates/rustdock-vina-core/src/parse_pdbqt.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/parse_pdbqt.h",
        rust_path: "crates/rustdock-vina-core/src/parse_pdbqt.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/potentials.h",
        rust_path: "crates/rustdock-vina-core/src/potentials.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/precalculate.h",
        rust_path: "crates/rustdock-vina-core/src/precalculate.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/quasi_newton.cpp",
        rust_path: "crates/rustdock-vina-core/src/quasi_newton.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/quasi_newton.h",
        rust_path: "crates/rustdock-vina-core/src/quasi_newton.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/quaternion.cpp",
        rust_path: "crates/rustdock-vina-core/src/quaternion.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/quaternion.h",
        rust_path: "crates/rustdock-vina-core/src/quaternion.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/random.cpp",
        rust_path: "crates/rustdock-vina-core/src/random.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/random.h",
        rust_path: "crates/rustdock-vina-core/src/random.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/scoring_function.h",
        rust_path: "crates/rustdock-vina-core/src/scoring_function.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/szv_grid.cpp",
        rust_path: "crates/rustdock-vina-core/src/szv_grid.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/szv_grid.h",
        rust_path: "crates/rustdock-vina-core/src/szv_grid.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/tree.h",
        rust_path: "crates/rustdock-vina-core/src/tree.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/triangular_matrix_index.h",
        rust_path: "crates/rustdock-vina-core/src/triangular_matrix_index.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/utils.cpp",
        rust_path: "crates/rustdock-vina-core/src/utils.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/utils.h",
        rust_path: "crates/rustdock-vina-core/src/utils.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/vina.cpp",
        rust_path: "crates/rustdock-vina-core/src/vina.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/lib/vina.h",
        rust_path: "crates/rustdock-vina-core/src/vina.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/main/main.cpp",
        rust_path: "crates/rustdock-vina-cli/src/main.rs",
        status: PortStatus::Ported,
    },
    SourcePort {
        reference_path: "reference/src/split/split.cpp",
        rust_path: "crates/rustdock-vina-split/src/lib.rs",
        status: PortStatus::Ported,
    },
];
