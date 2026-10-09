use rustdock_vina_core::vina::Vina;

#[derive(Default, Debug, Clone, Copy)]
pub struct Controls {
    pub lanes: Option<u32>,
    pub steps: Option<u32>,
    pub local_steps: Option<u32>,
}
#[derive(Debug, PartialEq, Eq)]
struct SearchSettings {
    lanes: u32,
    steps: u32,
    local_steps: u32,
}
impl Controls {
    pub fn is_set(&self) -> bool {
        self.lanes.is_some() || self.steps.is_some() || self.local_steps.is_some()
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.lanes.is_some_and(|n| n == 0 || n > 65_536) {
            return Err("--lanes must be between 1 and 65536".into());
        }
        if self.steps == Some(0) || self.local_steps == Some(0) {
            return Err("--steps and --local-steps must be positive integers".into());
        }
        Ok(())
    }
    fn resolve(
        &self,
        exhaustiveness: usize,
        atoms: usize,
        torsions: usize,
    ) -> Result<SearchSettings, String> {
        self.validate()?;
        if exhaustiveness == 0 {
            return Err("Exhaustiveness must be positive".into());
        }
        let global_steps = 105u64 * (50 + atoms as u64 + 10 * (6 + torsions as u64));
        let budget = (exhaustiveness as u64)
            .checked_mul(global_steps)
            .ok_or("Metal search budget overflow")?;
        let target = (budget / 32).clamp(256, 4096) as u32;
        let automatic_lanes = 1 << (31 - target.leading_zeros());
        let lanes = self.lanes.unwrap_or(automatic_lanes);
        let steps = match self.steps {
            Some(n) => n,
            None => u32::try_from(budget.div_ceil(lanes as u64))
                .map_err(|_| "Metal step count exceeds 32-bit limit")?,
        };
        let local_steps = self.local_steps.unwrap_or(
            u32::try_from((25 + atoms) / 3).map_err(|_| "Metal local-step count overflow")?,
        );
        Ok(SearchSettings {
            lanes,
            steps,
            local_steps,
        })
    }
}

pub fn dock(
    engine: &mut Vina,
    exhaustiveness: usize,
    modes: usize,
    min_rmsd: f64,
    max_evals: u32,
    verbosity: u32,
    controls: Controls,
) -> Result<(), String> {
    if max_evals != 0 {
        return Err("Metal does not support --max_evals; use 0 or --metal off".into());
    }
    if exhaustiveness == 0 || modes == 0 || !min_rmsd.is_finite() || min_rmsd <= 0.0 {
        return Err("Exhaustiveness, pose count and minimum RMSD must be positive".into());
    }
    engine.metal_grid()?;
    let settings = controls.resolve(
        exhaustiveness,
        engine.model.num_movable_atoms(),
        engine.model.get_size().ligands[0],
    )?;
    #[cfg(has_metal)]
    {
        let model = engine.model.clone();
        let poses = engine.poses.clone();
        let result = native::dock(engine, modes, min_rmsd, verbosity, settings);
        if result.is_err() {
            engine.model = model;
            engine.poses = poses;
        }
        result
    }
    #[cfg(not(has_metal))]
    {
        let _ = (verbosity, settings);
        Err("Metal is unavailable in this build. It requires macOS, a Metal device, and Apple Command Line Tools at build time. CPU builds use --metal off".into())
    }
}

#[cfg(has_metal)]
mod native {
    use super::*;
    use rustdock_vina_core::{
        cache::canonical_xs_map_type, common::Vec3, conf::OutputType,
        coords::add_to_output_container, quaternion::Quaternion, tree::Branch,
    };
    use serde_json::json;
    use std::{
        fs,
        io::{BufWriter, Write},
        path::PathBuf,
        process::Command,
    };
    struct Temp(PathBuf);
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    pub fn dock(
        engine: &mut Vina,
        modes: usize,
        min_rmsd: f64,
        verbosity: u32,
        settings: SearchSettings,
    ) -> Result<(), String> {
        let directory = std::env::temp_dir().join(format!(
            "rustdock-metal-{}-{}",
            std::process::id(),
            rustdock_vina_core::random::auto_seed()
        ));
        fs::create_dir(&directory).map_err(|e| e.to_string())?;
        let temp = Temp(directory);
        let cache = engine.metal_grid()?;
        let gd = cache.gd();
        let model = &engine.model;
        let heavy: Vec<_> = model
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(i, a)| (!a.base.atom_type.is_hydrogen()).then_some(i))
            .collect();
        let index = |i: usize| {
            heavy
                .iter()
                .position(|&a| a == i)
                .ok_or("Torsion axis requires a scored heavy atom".to_string())
        };
        let mut types = Vec::new();
        for &i in &heavy {
            let t = canonical_xs_map_type(model.atoms[i].base.atom_type.xs)
                .ok_or("Unsupported Metal atom type")?;
            if !types.contains(&t) {
                types.push(t);
            }
        }
        let root = model.ligands[0].body.node.atom_frame.frame.origin();
        let atoms: Vec<_> = heavy
            .iter()
            .map(|&i| {
                let p = model.coords[i] - root;
                let t = model.atoms[i].base.atom_type.xs;
                [
                    p[0],
                    p[1],
                    p[2],
                    types.iter().position(|&s| s == t).unwrap() as f64,
                ]
            })
            .collect();
        fn subtree_indices(branch: &Branch, out: &mut Vec<usize>) {
            out.extend(
                branch.node.axis_frame.atom_frame.range.begin
                    ..branch.node.axis_frame.atom_frame.range.end,
            );
            for child in &branch.children {
                subtree_indices(child, out);
            }
        }
        fn torsions(
            branch: &Branch,
            model: &rustdock_vina_core::model::Model,
            heavy: &[usize],
            out: &mut Vec<[u32; 4]>,
        ) -> Result<(), String> {
            let frame = &branch.node.axis_frame;
            let child = model
                .coords
                .iter()
                .position(|p| (*p - frame.atom_frame.frame.origin()).norm() < 1e-7)
                .ok_or("Metal torsion endpoint missing")?;
            let parent = model.atoms[child]
                .bonds
                .iter()
                .find_map(|bond| {
                    if !bond.rotatable || bond.connected_atom_index.in_grid {
                        return None;
                    }
                    let i = bond.connected_atom_index.i;
                    let delta = model.coords[child] - model.coords[i];
                    (delta.norm() > 1e-8
                        && (delta * (1.0 / delta.norm()) - frame.axis).norm() < 1e-5)
                        .then_some(i)
                })
                .ok_or("Metal torsion parent bond missing")?;
            let idx = |i| {
                heavy
                    .iter()
                    .position(|&a| a == i)
                    .ok_or("Metal requires heavy-atom torsion axes")
            };
            let mut descendants = vec![];
            subtree_indices(branch, &mut descendants);
            let mut mask = 0u64;
            for i in descendants {
                if let Some(h) = heavy.iter().position(|&a| a == i) {
                    mask |= 1u64 << h;
                }
            }
            out.push([
                idx(parent)? as u32,
                idx(child)? as u32,
                mask as u32,
                (mask >> 32) as u32,
            ]);
            for child in &branch.children {
                torsions(child, model, heavy, out)?;
            }
            Ok(())
        }
        let mut branches = vec![];
        for branch in &model.ligands[0].body.children {
            torsions(branch, model, &heavy, &mut branches)?;
        }
        let pairs = model.ligands[0]
            .pairs
            .iter()
            .filter_map(|p| match (index(p.a), index(p.b)) {
                (Ok(a), Ok(b)) => Some([
                    a as u32,
                    b as u32,
                    model.atoms[p.a].base.atom_type.xs as u32,
                    model.atoms[p.b].base.atom_type.xs as u32,
                ]),
                _ => None,
            })
            .collect::<Vec<_>>();
        let path = temp.0.join("maps.bin");
        let mut maps = BufWriter::new(fs::File::create(path).map_err(|e| e.to_string())?);
        for &t in &types {
            let g = cache.grid_for_type(t).ok_or("Metal map missing")?;
            for z in 0..g.data.dim2() {
                for y in 0..g.data.dim1() {
                    for x in 0..g.data.dim0() {
                        maps.write_all(&(*g.data.get(x, y, z) as f32).to_le_bytes())
                            .map_err(|e| e.to_string())?;
                    }
                }
            }
        }
        maps.flush().map_err(|e| e.to_string())?;
        // Calibrate inside the box even when an input ligand starts far from its pocket.
        let calibration_position = Vec3::new(
            (gd[0].begin + gd[0].end) / 2.0,
            (gd[1].begin + gd[1].end) / 2.0,
            (gd[2].begin + gd[2].end) / 2.0,
        );
        let mut calibration = model.clone();
        calibration.coords.fill(calibration_position);
        let input = json!({"dims":gd.map(|g|g.n_voxels+1),"origin":gd.map(|g|g.begin),"root":calibration_position.data,"spacing":(gd[0].end-gd[0].begin)/gd[0].n_voxels as f64,"atoms":atoms,"torsions":branches,"pairs":pairs,"lanes":settings.lanes,"steps":settings.steps,"localSteps":settings.local_steps,"seed":engine.seed as u32,"gyrationRadius":model.gyration_radius(0),"expectedGridEnergy":cache.eval(&calibration,1000.0)});
        fs::write(temp.0.join("input.json"), input.to_string()).map_err(|e| e.to_string())?;
        fs::write(
            temp.0.join("grid_score.metal"),
            include_bytes!("../metal/grid_score.metal"),
        )
        .map_err(|e| e.to_string())?;
        let helper = temp.0.join("rustdock-vina-metal");
        fs::write(
            &helper,
            include_bytes!(concat!(env!("OUT_DIR"), "/rustdock-vina-metal")),
        )
        .map_err(|e| e.to_string())?;
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&helper, fs::Permissions::from_mode(0o700))
            .map_err(|e| e.to_string())?;
        let output = Command::new(&helper)
            .arg(&temp.0)
            .output()
            .map_err(|e| e.to_string())?;
        if !output.status.success() {
            return Err(format!(
                "Metal search failed: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        let report: serde_json::Value =
            serde_json::from_slice(&output.stdout).map_err(|e| e.to_string())?;
        let bytes = fs::read(temp.0.join("candidates.bin")).map_err(|e| e.to_string())?;
        if bytes.is_empty() || bytes.len() % 64 != 0 {
            return Err("Invalid Metal candidate buffer".into());
        }
        let template = engine.model.get_initial_conf();
        let mut candidates = Vec::new();
        for record in bytes.chunks_exact(64) {
            let values: Vec<_> = record
                .chunks_exact(4)
                .map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64)
                .collect();
            if values.iter().any(|v| !v.is_finite()) {
                continue;
            }
            let mut c = template.clone();
            let q = Quaternion::new(values[7], values[4], values[5], values[6]);
            if q.norm() < 1e-8 {
                continue;
            }
            c.ligands[0].rigid.position = Vec3::new(values[0], values[1], values[2]);
            c.ligands[0].rigid.orientation = q.normalized();
            for (i, angle) in c.ligands[0].torsions.iter_mut().enumerate() {
                *angle = values[8 + i];
            }
            engine.model.set(&c);
            let mut pose = OutputType::new(c, values[3]);
            pose.coords = engine.model.get_heavy_atom_movable_coords();
            candidates.push(pose);
        }
        candidates.sort_by(rustdock_vina_core::conf::output_compare);
        let mut clustered = vec![];
        for pose in candidates {
            add_to_output_container(&mut clustered, pose, min_rmsd, modes);
        }
        engine.finalize_poses(clustered)?;
        if verbosity > 0 {
            println!(
                "  Metal: {} · {} lanes × {} steps · {} local steps · GPU search {:.2} s · Rust refinement",
                report["device"].as_str().unwrap_or("Unknown device"),
                report["lanes"],
                report["steps"],
                report["local_steps"],
                report["search_ms"].as_f64().unwrap_or(0.0) / 1000.0
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn default_budget_matches_existing_scheduler() {
        assert_eq!(
            Controls::default().resolve(8, 40, 8).unwrap(),
            SearchSettings {
                lanes: 4096,
                steps: 48,
                local_steps: 21
            }
        );
    }
    #[test]
    fn independent_overrides_preserve_unspecified_defaults() {
        let c = Controls {
            lanes: Some(1000),
            ..Default::default()
        };
        assert_eq!(
            c.resolve(8, 40, 8).unwrap(),
            SearchSettings {
                lanes: 1000,
                steps: 194,
                local_steps: 21
            }
        );
        let c = Controls {
            steps: Some(10),
            local_steps: Some(4),
            ..Default::default()
        };
        assert_eq!(
            c.resolve(8, 40, 8).unwrap(),
            SearchSettings {
                lanes: 4096,
                steps: 10,
                local_steps: 4
            }
        );
        let c = Controls {
            lanes: Some(8192),
            steps: Some(20),
            local_steps: Some(5),
        };
        assert_eq!(
            c.resolve(8, 40, 8).unwrap(),
            SearchSettings {
                lanes: 8192,
                steps: 20,
                local_steps: 5
            }
        );
    }
    #[test]
    fn invalid_controls_and_budget_overflow_return_errors() {
        for lanes in [0, 65537] {
            assert!(Controls {
                lanes: Some(lanes),
                ..Default::default()
            }
            .validate()
            .is_err());
        }
        assert!(Controls {
            steps: Some(0),
            ..Default::default()
        }
        .validate()
        .is_err());
        assert!(Controls {
            local_steps: Some(0),
            ..Default::default()
        }
        .validate()
        .is_err());
        assert!(Controls::default().resolve(usize::MAX, 40, 8).is_err());
    }
}
