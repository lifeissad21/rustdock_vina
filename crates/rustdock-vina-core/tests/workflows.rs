//! Golden scores captured from official AutoDock Vina 1.2.7, CPU=1, seed=42.
use rustdock_vina_core::{common::Vec3, vina::Vina};
use std::path::PathBuf;

fn fixture(path: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../reference/example")
        .join(path)
        .to_str()
        .unwrap()
        .into()
}
fn box_dims(path: &str) -> (Vec3, Vec3) {
    let text = std::fs::read_to_string(fixture(path)).unwrap();
    let value = |key: &str| {
        text.lines()
            .find_map(|line| {
                line.split_once('=')
                    .filter(|(k, _)| k.trim() == key)
                    .map(|(_, v)| v.trim().parse::<f64>().unwrap())
            })
            .unwrap()
    };
    (
        Vec3::new(value("center_x"), value("center_y"), value("center_z")),
        Vec3::new(value("size_x"), value("size_y"), value("size_z")),
    )
}
fn basic(scoring: &str) -> Vina {
    basic_with_refinement(scoring, false)
}
fn basic_with_refinement(scoring: &str, no_refine: bool) -> Vina {
    let mut v = Vina::new(scoring, 1, 42, no_refine).unwrap();
    v.set_receptor(&fixture("basic_docking/solution/1iep_receptor.pdbqt"), "")
        .unwrap();
    v.set_ligand_from_file(&fixture("basic_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    let (center, size) = box_dims("basic_docking/solution/1iep_receptor.box.txt");
    v.compute_vina_maps(center, size, 0.375, true).unwrap();
    v
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 0.0006, "{actual} != {expected}");
}
#[test]
fn vina_score_and_local_optimization_match_official_binary() {
    let mut v = basic("vina");
    assert!(v.metal_grid().is_ok());
    let mut invalid = rustdock_vina_core::conf::OutputType::new(v.model.get_initial_conf(), 0.0);
    invalid.c.ligands[0].torsions.clear();
    assert!(v.finalize_poses(vec![invalid]).is_err());
    let es = v.score().unwrap();
    near(es[0], -12.513);
    near(es[1], -17.634);
    near(es[5], -0.485);
    near(es[6], 5.121);
    near(v.optimize(0).unwrap()[0], -13.170);
}
#[test]
fn vinardo_score_matches_official_binary() {
    let mut v = basic("vinardo");
    let es = v.score().unwrap();
    near(es[0], -9.180);
    near(es[1], -12.936);
    near(es[5], -0.457);
}
#[test]
fn ad4_maps_score_matches_official_binary() {
    let mut v = Vina::new("ad4", 1, 42, false).unwrap();
    v.set_ligand_from_file(&fixture("basic_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    v.load_maps(&fixture("basic_docking/solution/1iep_receptor"))
        .unwrap();
    let es = v.score().unwrap();
    near(es[0], -14.375);
    near(es[1], -16.463);
    near(es[5], -1.147);
    near(es[6], 2.088);
    assert!(v.optimize(5).unwrap().iter().all(|e| e.is_finite()));
}
#[test]
fn flexible_receptor_score_and_coordinate_assembly_match_reference() {
    let mut v = Vina::new("vina", 1, 42, false).unwrap();
    v.set_receptor(
        &fixture("flexible_docking/solution/1fpu_receptor_rigid.pdbqt"),
        &fixture("flexible_docking/solution/1fpu_receptor_flex.pdbqt"),
    )
    .unwrap();
    v.set_ligand_from_file(&fixture("flexible_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    let coords = v.model.coords.clone();
    v.model.set(&v.model.get_initial_conf());
    for (a, b) in coords.iter().zip(&v.model.coords) {
        assert!((*a - *b).norm() < 1e-10);
    }
    let (center, size) = box_dims("flexible_docking/solution/1fpu_receptor.box.txt");
    v.compute_vina_maps(center, size, 0.375, false).unwrap();
    let es = v.score().unwrap();
    near(es[0], -4.275);
    near(es[1] + es[2], -5.951);
    near(es[3] + es[4] + es[5], 0.772);
    near(es[7], 0.844);
    assert!(v.model.num_atoms() > v.model.num_movable_atoms());
}
#[test]
fn macrocycle_score_matches_official_binary() {
    let mut v = Vina::new("vina", 1, 42, false).unwrap();
    v.set_receptor(
        &fixture("docking_with_macrocycles/solution/BACE_1_receptor.pdbqt"),
        "",
    )
    .unwrap();
    v.set_ligand_from_file(&fixture(
        "docking_with_macrocycles/solution/BACE_1_ligand.pdbqt",
    ))
    .unwrap();
    let (center, size) = box_dims("docking_with_macrocycles/solution/BACE_1_receptor_vina_box.txt");
    v.compute_vina_maps(center, size, 0.375, false).unwrap();
    let es = v.score().unwrap();
    near(es[0], -7.628);
    near(es[1], -17.216);
    near(es[5], -0.878);
    assert!(!v.model.glue_pairs.is_empty());
}
#[test]
fn docking_is_reproducible_and_emits_parseable_poses() {
    let mut v = basic("vina");
    v.global_search(2, 3, 1.0, 2000).unwrap();
    let poses = v.get_poses(3, 100.0).unwrap();
    let energies = v.get_poses_energies(3, 100.0).unwrap();
    assert!(!energies.is_empty());
    assert!(energies.iter().all(|e| e.iter().all(|x| x.is_finite())));
    assert!(energies.windows(2).all(|e| e[0][0] <= e[1][0]));
    assert_eq!(poses.matches("MODEL ").count(), energies.len());
    assert_eq!(
        v.get_poses_coordinates(3, 100.0).unwrap()[0].len(),
        3 * v.model.num_atoms()
    );
    let first = poses
        .lines()
        .skip(1)
        .take_while(|l| *l != "ENDMDL")
        .collect::<Vec<_>>()
        .join("\n");
    rustdock_vina_core::parse_pdbqt::parse_ligand_pdbqt_from_string(
        &first,
        v.model.atom_typing_used(),
    )
    .unwrap();
    v.set_ligand_from_file(&fixture("basic_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    v.global_search(2, 3, 1.0, 2000).unwrap();
    assert_eq!(v.get_poses(3, 100.0).unwrap(), poses);
}
#[test]
fn multiple_ligands_keep_indexed_bonds_and_contexts() {
    let mut v = Vina::new("vina", 1, 42, false).unwrap();
    let path = fixture("basic_docking/solution/1iep_ligand.pdbqt");
    v.set_ligands_from_files(&[path.clone(), path]).unwrap();
    let n = v.model.num_atoms() / 2;
    assert_eq!(v.model.ligands[1].range.begin, n);
    for a in &v.model.atoms[n..] {
        assert!(a.bonds.iter().all(|b| b.connected_atom_index.i >= n));
    }
    assert!(!v.model.inter_pairs.is_empty());
    let coords = v.model.coords.clone();
    v.model.set(&v.model.get_initial_conf());
    for (a, b) in coords.iter().zip(&v.model.coords) {
        assert!((*a - *b).norm() < 1e-10);
    }
}
#[test]
fn map_write_read_preserves_scoring_and_rejects_corrupt_maps() {
    let v = basic_with_refinement("vina", true);
    let dir = std::env::temp_dir().join(format!(
        "rustdock-map-test-{}",
        rustdock_vina_core::random::auto_seed()
    ));
    std::fs::create_dir(&dir).unwrap();
    let prefix = dir.join("receptor").to_str().unwrap().to_string();
    v.write_maps(&prefix).unwrap();
    let mut loaded = Vina::new("vina", 1, 42, true).unwrap();
    loaded
        .set_ligand_from_file(&fixture("basic_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    loaded.load_maps(&prefix).unwrap();
    let mut grid_only = basic_with_refinement("vina", true);
    // Four-decimal map serialization changes each atom contribution by at most 0.00005.
    let a = loaded.score().unwrap()[0];
    let b = grid_only.score().unwrap()[0];
    assert!((a - b).abs() < 0.005);
    let map = dir.join("receptor.C_H.map");
    std::fs::write(map, "broken map\n").unwrap();
    assert!(loaded.load_maps(&prefix).is_err());
    assert!(loaded.score().is_ok());
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn invalid_input_does_not_replace_valid_ligand_and_errors_are_explicit() {
    let mut v = Vina::new("vina", 1, 42, false).unwrap();
    assert!(v.score().is_err());
    assert!(v.global_search(1, 1, 1.0, 0).is_err());
    v.set_ligand_from_file(&fixture("basic_docking/solution/1iep_ligand.pdbqt"))
        .unwrap();
    let model = v.model.clone();
    assert!(v.set_ligand_from_string("not pdbqt").is_err());
    assert_eq!(model, v.model);
    assert!(v
        .compute_vina_maps(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0),
            0.0,
            false
        )
        .is_err());
}
