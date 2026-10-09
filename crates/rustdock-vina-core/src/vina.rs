//! Native docking workflow. No external Vina executable is used by this engine.
use crate::ad4cache::Ad4Cache;
use crate::cache::Cache;
use crate::common::{Fl, Vec3, MAX_FL};
use crate::conf::{Change, OutputContainer, OutputType};
use crate::conf_independent::ConfIndependentInputs;
use crate::coords::add_to_output_container;
use crate::grid_dim::GridDim;
use crate::igrid::IGrid;
use crate::model::{Model, SearchModel};
use crate::monte_carlo::MonteCarlo;
use crate::non_cache::NonCache;
use crate::parallel_mc::ParallelMc;
use crate::parse_pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt};
use crate::precalculate::{Precalculate, PrecalculateByAtom};
use crate::quasi_newton::QuasiNewton;
use crate::random::{auto_seed, Rng64};
use crate::scoring_function::{ScoringFunction, ScoringFunctionChoice};

pub struct Vina {
    pub model: Model,
    receptor: Model,
    sf: ScoringFunction,
    choice: ScoringFunctionChoice,
    precalculate: Option<PrecalculateByAtom>,
    grid: Option<DockingGrid>,
    explicit: Option<NonCache>,
    pub poses: OutputContainer,
    pub seed: u64,
    pub cpu: usize,
    no_refine: bool,
}
#[derive(Clone)]
enum DockingGrid {
    Vina(Cache),
    Ad4(Ad4Cache),
}
impl DockingGrid {
    fn dims(&self) -> crate::grid_dim::GridDims {
        match self {
            Self::Vina(g) => g.gd(),
            Self::Ad4(g) => g.gd(),
        }
    }
    fn corners(&self) -> (Vec3, Vec3) {
        match self {
            Self::Vina(g) => (g.corner1(), g.corner2()),
            Self::Ad4(g) => (g.corner1(), g.corner2()),
        }
    }
    fn within(&self, m: &Model) -> bool {
        match self {
            Self::Vina(g) => g.is_in_grid(m, 0.0001),
            Self::Ad4(g) => g.is_in_grid(m, 0.0001),
        }
    }
    fn complete(&self, m: &Model) -> bool {
        let types = m.get_movable_atom_types(m.atom_typing_used());
        match self {
            Self::Vina(g) => g.are_atom_types_grid_initialized(&types),
            Self::Ad4(g) => g.are_atom_types_grid_initialized(&types),
        }
    }
}
impl IGrid<Model> for DockingGrid {
    fn eval(&self, m: &Model, v: Fl) -> Fl {
        match self {
            Self::Vina(g) => g.eval(m, v),
            Self::Ad4(g) => g.eval(m, v),
        }
    }
    fn eval_intra(&self, m: &mut Model, v: Fl) -> Fl {
        match self {
            Self::Vina(g) => g.eval_intra(m, v),
            Self::Ad4(g) => g.eval_intra(m, v),
        }
    }
    fn eval_deriv(&self, m: &mut Model, v: Fl) -> Fl {
        match self {
            Self::Vina(g) => g.eval_deriv(m, v),
            Self::Ad4(g) => g.eval_deriv(m, v),
        }
    }
}
const AUTHENTIC: Vec3 = Vec3::new(1000.0, 1000.0, 1000.0);
impl Vina {
    pub fn new(scoring: &str, cpu: usize, seed: u64, no_refine: bool) -> Result<Self, String> {
        let (choice, weights) = match scoring {
            "vina" => (
                ScoringFunctionChoice::Vina,
                vec![
                    -0.035579,
                    -0.005156,
                    0.840245,
                    -0.035069,
                    -0.587439,
                    50.0,
                    5.0 * 0.05846 / 0.1 - 1.0,
                ],
            ),
            "vinardo" => (
                ScoringFunctionChoice::Vinardo,
                vec![-0.045, 0.8, -0.035, -0.600, 50.0, 5.0 * 0.05846 / 0.1 - 1.0],
            ),
            "ad4" => (
                ScoringFunctionChoice::Ad42,
                vec![0.1662, 0.1209, 0.1406, 0.1322, 50.0, 0.2983],
            ),
            _ => return Err("Scoring function must be vina, vinardo or ad4".into()),
        };
        let sf = ScoringFunction::new(choice, weights);
        Ok(Self {
            model: Model::new(sf.atom_typing()),
            receptor: Model::new(sf.atom_typing()),
            sf,
            choice,
            precalculate: None,
            grid: None,
            explicit: None,
            poses: vec![],
            seed: if seed == 0 { auto_seed() } else { seed },
            cpu: if cpu == 0 {
                std::thread::available_parallelism().map_or(1, usize::from)
            } else {
                cpu
            },
            no_refine,
        })
    }
    pub fn set_weights(&mut self, mut weights: Vec<Fl>) -> Result<(), String> {
        let count = match self.choice {
            ScoringFunctionChoice::Vina => 7,
            _ => 6,
        };
        if weights.len() != count || weights.iter().any(|w| !w.is_finite()) {
            return Err(format!("Expected {count} finite weights"));
        }
        if self.choice != ScoringFunctionChoice::Ad42 {
            weights[count - 1] = 5.0 * weights[count - 1] / 0.1 - 1.0;
        }
        self.sf = ScoringFunction::new(self.choice, weights);
        self.grid = None;
        self.explicit = None;
        self.poses.clear();
        self.recalculate();
        Ok(())
    }
    pub fn set_receptor(&mut self, rigid: &str, flex: &str) -> Result<(), String> {
        if rigid.is_empty() && flex.is_empty() {
            return Err("Specify a rigid receptor or flexible residues".into());
        }
        if self.choice == ScoringFunctionChoice::Ad42 && !rigid.is_empty() {
            return Err("AD4 uses precomputed maps; omit the rigid receptor".into());
        }
        let receptor = parse_receptor_pdbqt(rigid, flex, self.sf.atom_typing())?;
        self.model = receptor.clone();
        self.receptor = receptor;
        self.precalculate = None;
        self.grid = None;
        self.explicit = None;
        self.poses.clear();
        Ok(())
    }
    pub fn set_ligands_from_strings(&mut self, ligands: &[String]) -> Result<(), String> {
        if ligands.is_empty() {
            return Err("No ligands supplied".into());
        }
        let mut model = self.receptor.clone();
        for text in ligands {
            model.append(
                &parse_ligand_pdbqt_from_string(text, self.sf.atom_typing())
                    .map_err(|e| e.to_string())?,
            );
        }
        if self.grid.as_ref().is_some_and(|g| !g.complete(&model)) {
            return Err("Loaded maps do not cover all ligand atom types".into());
        }
        self.model = model;
        self.poses.clear();
        self.recalculate();
        Ok(())
    }
    pub fn set_ligand_from_string(&mut self, text: &str) -> Result<(), String> {
        self.set_ligands_from_strings(&[text.into()])
    }
    pub fn set_ligand_from_file(&mut self, path: &str) -> Result<(), String> {
        self.set_ligand_from_string(
            &std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?,
        )
    }
    pub fn set_ligands_from_files(&mut self, paths: &[String]) -> Result<(), String> {
        let texts = paths
            .iter()
            .map(|p| std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}")))
            .collect::<Result<Vec<_>, _>>()?;
        self.set_ligands_from_strings(&texts)
    }
    fn recalculate(&mut self) {
        self.precalculate = if self.model.num_ligands() > 0 {
            Some(PrecalculateByAtom::new(&self.sf, &self.model, MAX_FL, 32.0))
        } else {
            None
        };
    }
    pub fn compute_vina_maps(
        &mut self,
        center: Vec3,
        size: Vec3,
        spacing: Fl,
        even: bool,
    ) -> Result<(), String> {
        if self.choice == ScoringFunctionChoice::Ad42 {
            return Err("AD4 requires precomputed AutoGrid maps".into());
        }
        if self.model.grid_atoms.is_empty() {
            return Err("Rigid receptor has not been initialized".into());
        }
        if !spacing.is_finite()
            || spacing <= 0.0
            || center.data.iter().any(|v| !v.is_finite())
            || size.data.iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(
                "Grid center must be finite; sizes and spacing must be finite and positive".into(),
            );
        }
        let mut dims = [GridDim::disabled(); 3];
        let mut points = 1usize;
        for i in 0..3 {
            let voxels = (size[i] / spacing).ceil();
            if !voxels.is_finite() || voxels > usize::MAX as f64 - 2.0 {
                return Err("Grid size overflow".into());
            }
            let mut n = voxels as usize;
            if even && n % 2 != 0 {
                n += 1;
            }
            points = points.checked_mul(n + 1).ok_or("Grid size overflow")?;
            let half = n as Fl * spacing / 2.0;
            dims[i] = GridDim::new(center[i] - half, center[i] + half, n);
        }
        if points > 100_000_000 {
            return Err(
                "Grid exceeds 100 million points; increase spacing or reduce box size".into(),
            );
        }
        let p = Precalculate::new(&self.sf, MAX_FL, 32.0);
        let mut grid = Cache::with_dims(dims, 1e6);
        let types = if self.model.num_ligands() > 0 {
            self.model.get_movable_atom_types(self.sf.atom_typing())
        } else {
            self.sf.atom_types()
        };
        grid.populate(&self.model, &p, &types);
        self.explicit = if self.no_refine {
            None
        } else {
            Some(NonCache::new(&self.model, dims, p, 1e6))
        };
        self.grid = Some(DockingGrid::Vina(grid));
        self.poses.clear();
        Ok(())
    }
    pub fn load_maps(&mut self, prefix: &str) -> Result<(), String> {
        let grid = if self.choice == ScoringFunctionChoice::Ad42 {
            let mut g = Ad4Cache::new(1e6);
            g.read(prefix)?;
            DockingGrid::Ad4(g)
        } else {
            let mut g = Cache::new(1e6);
            g.read(prefix)?;
            DockingGrid::Vina(g)
        };
        if self.model.num_ligands() > 0 && !grid.complete(&self.model) {
            return Err("Maps are missing required ligand atom types".into());
        }
        self.explicit = if !self.no_refine && !self.model.grid_atoms.is_empty() {
            Some(NonCache::new(
                &self.model,
                grid.dims(),
                Precalculate::new(&self.sf, MAX_FL, 32.0),
                1e6,
            ))
        } else {
            None
        };
        self.grid = Some(grid);
        self.poses.clear();
        Ok(())
    }
    pub fn write_maps(&self, prefix: &str) -> Result<(), String> {
        self.write_maps_with_metadata(prefix, "NULL", "NULL", "NULL")
    }
    pub fn write_maps_with_metadata(
        &self,
        prefix: &str,
        gpf: &str,
        fld: &str,
        receptor: &str,
    ) -> Result<(), String> {
        let g = self.grid.as_ref().ok_or("Maps have not been initialized")?;
        let types = if self.model.num_ligands() > 0 {
            self.model.get_movable_atom_types(self.sf.atom_typing())
        } else {
            self.sf.atom_types()
        };
        match g {
            DockingGrid::Vina(g) => g.write(prefix, &types, gpf, fld, receptor),
            DockingGrid::Ad4(g) => g.write(prefix, &types, gpf, fld, receptor),
        }
    }
    fn ready(&self) -> Result<(), String> {
        if self.model.num_ligands() == 0 {
            return Err("Ligands have not been initialized".into());
        }
        let g = self.grid.as_ref().ok_or("Maps have not been initialized")?;
        if !g.complete(&self.model) {
            return Err("Required atom-type maps are missing".into());
        }
        Ok(())
    }
    pub fn score(&mut self) -> Result<[Fl; 8], String> {
        self.ready()?;
        let g = self.grid.as_ref().unwrap();
        if !g.within(&self.model) {
            return Err("Ligand is outside the grid box".into());
        }
        let intra = if self.choice == ScoringFunctionChoice::Ad42 {
            0.0
        } else {
            self.model
                .eval_intramolecular(self.precalculate.as_ref().unwrap(), g, AUTHENTIC)
        };
        Ok(self.score_with_intra(intra))
    }
    pub fn score_with_unbound(&mut self, unbound: Fl) -> Result<[Fl; 8], String> {
        self.ready()?;
        if !unbound.is_finite() {
            return Err("Unbound energy must be finite".into());
        }
        if !self.grid.as_ref().unwrap().within(&self.model) {
            return Err("Ligand is outside the grid box".into());
        }
        Ok(self.score_with_intra(unbound))
    }
    fn score_with_intra(&mut self, unbound: Fl) -> [Fl; 8] {
        let p = self.precalculate.as_ref().unwrap();
        let (all, flex) = if let Some(g) = &self.explicit {
            (
                g.eval(&self.model, 1000.0),
                g.eval_intra(&mut self.model, 1000.0),
            )
        } else {
            let g = self.grid.as_ref().unwrap();
            (
                g.eval(&self.model, 1000.0),
                g.eval_intra(&mut self.model, 1000.0),
            )
        };
        let inter_pairs = self.model.eval_inter(p, AUTHENTIC);
        let other = self.model.evalo(p, AUTHENTIC);
        let lig_intra = self.model.evali(p, AUTHENTIC);
        let intra = flex + other + lig_intra;
        let inter = all - flex + inter_pairs;
        let inputs = ConfIndependentInputs::from_model(&self.model);
        let (total, torsion, unbound) = if self.choice == ScoringFunctionChoice::Ad42 {
            let t = self.sf.conf_independent_from_inputs(&inputs, 0.0);
            (inter + t, t, intra)
        } else {
            let raw = inter + intra - unbound;
            let total = self.sf.conf_independent_from_inputs(&inputs, raw);
            (total, total - raw, unbound)
        };
        [
            total,
            all - flex,
            inter_pairs,
            flex,
            other,
            lig_intra,
            torsion,
            unbound,
        ]
    }
    pub fn optimize(&mut self, max_steps: u32) -> Result<[Fl; 8], String> {
        self.ready()?;
        if !self.grid.as_ref().unwrap().within(&self.model) {
            return Err("Ligand is outside the grid box".into());
        }
        let c = self
            .poses
            .first()
            .map_or_else(|| self.model.get_initial_conf(), |p| p.c.clone());
        let mut out = OutputType::new(c, 0.0);
        let steps = if max_steps == 0 {
            (25 + self.model.num_movable_atoms()) as u32 / 3
        } else {
            max_steps
        };
        let grid = self.grid.as_ref().unwrap();
        let p = self.precalculate.as_ref().unwrap();
        for _ in 0..5 {
            optimize_model(&mut self.model, p, grid, &mut out, steps);
            if grid.within(&self.model) {
                break;
            }
        }
        self.poses.clear();
        self.score()
    }
    pub fn randomize(&mut self, max_steps: usize) -> Result<(), String> {
        self.ready()?;
        if max_steps == 0 {
            return Err("Randomization steps must be positive".into());
        }
        let (a, b) = self.grid.as_ref().unwrap().corners();
        let mut rng = Rng64::new(self.seed);
        let mut best = self.model.get_initial_conf();
        let mut penalty = MAX_FL;
        for _ in 0..max_steps {
            let mut c = best.clone();
            c.randomize(a, b, &mut rng);
            self.model.set(&c);
            let e = self.model.clash_penalty();
            if e < penalty {
                best = c;
                penalty = e;
            }
            if penalty == 0.0 {
                break;
            }
        }
        self.model.set(&best);
        self.poses.clear();
        Ok(())
    }
    pub fn global_search(
        &mut self,
        exhaustiveness: usize,
        n_poses: usize,
        min_rmsd: Fl,
        max_evals: u32,
    ) -> Result<(), String> {
        self.ready()?;
        if exhaustiveness == 0 || n_poses == 0 || !min_rmsd.is_finite() || min_rmsd <= 0.0 {
            return Err("Exhaustiveness, pose count and minimum RMSD must be positive".into());
        }
        let heuristic =
            self.model.num_movable_atoms() + 10 * self.model.get_size().num_degrees_of_freedom();
        let local_steps = (25 + self.model.num_movable_atoms()) as u32 / 3;
        let pmc = ParallelMc {
            mc: MonteCarlo {
                global_steps: (105 * (50 + heuristic)) as u32,
                local_steps,
                max_evals,
                min_rmsd,
                num_saved_mins: n_poses,
                hunt_cap: Vec3::new(10.0, 10.0, 10.0),
                ..Default::default()
            },
            num_tasks: exhaustiveness,
            num_threads: self.cpu,
            display_progress: false,
        };
        let grid = self.grid.as_ref().unwrap();
        let p = self.precalculate.as_ref().unwrap();
        let model = SearchModel {
            model: self.model.clone(),
            precalculate: p,
            grid,
        };
        let (a, b) = grid.corners();
        let mut poses = vec![];
        pmc.run(
            &model,
            &mut poses,
            a,
            b,
            &mut Rng64::new(self.seed),
            Option::<fn(f64)>::None,
        );
        let mut unique = vec![];
        for pose in poses {
            add_to_output_container(&mut unique, pose, min_rmsd, n_poses);
        }
        self.finalize_poses(unique)
    }
    pub fn finalize_poses(&mut self, mut unique: OutputContainer) -> Result<(), String> {
        self.ready()?;
        let size = self.model.get_size();
        if unique.is_empty()
            || unique.iter().any(|pose| {
                !pose.e.is_finite()
                    || pose.c.ligands.len() != size.ligands.len()
                    || pose.c.flex.len() != size.flex.len()
                    || pose.c.ligands.iter().zip(&size.ligands).any(|(l, n)| {
                        let q = l.rigid.orientation;
                        l.torsions.len() != *n
                            || l.torsions
                                .iter()
                                .chain(l.rigid.position.data.iter())
                                .any(|v| !v.is_finite())
                            || [q.r, q.i, q.j, q.k].iter().any(|v| !v.is_finite())
                            || (q.norm() - 1.0).abs() > 1e-5
                    })
                    || pose.c.flex.iter().zip(&size.flex).any(|(f, n)| {
                        f.torsions.len() != *n || f.torsions.iter().any(|v| !v.is_finite())
                    })
            })
        {
            return Err("No valid search candidates".into());
        }
        let local_steps = (25 + self.model.num_movable_atoms()) as u32 / 3;
        let p = self.precalculate.as_ref().unwrap();
        let grid = self.grid.as_ref().unwrap();
        if let Some(explicit) = &mut self.explicit {
            for pose in &mut unique {
                for pass in 0..5 {
                    explicit.slope = 100.0 * 10.0_f64.powi(2 * pass);
                    optimize_model(&mut self.model, p, explicit, pose, local_steps);
                    if explicit.within(&self.model, 0.0001) {
                        break;
                    }
                }
                explicit.slope = 1e6;
                self.model.set(&pose.c);
                pose.coords = self.model.get_heavy_atom_movable_coords();
                pose.e = explicit.eval(&self.model, 1000.0)
                    + self.model.eval_inter(p, AUTHENTIC)
                    + self.model.evalo(p, AUTHENTIC)
                    + self.model.evali(p, AUTHENTIC);
            }
        }
        unique.sort_by(crate::conf::output_compare);
        let best = unique.first().ok_or("No docking poses found")?;
        self.model.set(&best.c);
        let unbound = if self.choice == ScoringFunctionChoice::Ad42 {
            0.0
        } else if let Some(g) = &self.explicit {
            self.model.eval_intramolecular(p, g, AUTHENTIC)
        } else {
            self.model.eval_intramolecular(p, grid, AUTHENTIC)
        };
        for pose in &mut unique {
            self.model.set(&pose.c);
            let es = self.score_with_intra(unbound);
            pose.e = es[0];
            pose.inter = es[1] + es[2];
            pose.intra = es[3] + es[4] + es[5];
            pose.conf_independent = es[6];
            pose.unbound = es[7];
            pose.total = pose.inter + pose.intra;
        }
        unique.sort_by(crate::conf::output_compare);
        let mut reference = self.model.clone();
        reference.set(&unique[0].c);
        for pose in &mut unique {
            self.model.set(&pose.c);
            pose.lb = self.model.rmsd_lower_bound(&reference);
            pose.ub = self.model.rmsd_upper_bound_model(&reference);
        }
        self.model.set(&unique[0].c);
        self.poses = unique;
        Ok(())
    }
    pub fn metal_grid(&self) -> Result<&Cache, String> {
        self.ready()?;
        if self.choice != ScoringFunctionChoice::Vina
            || self.sf.weights()
                != [
                    -0.035579,
                    -0.005156,
                    0.840245,
                    -0.035069,
                    -0.587439,
                    50.0,
                    5.0 * 0.05846 / 0.1 - 1.0,
                ]
        {
            return Err("Metal currently supports the default Vina force field only".into());
        }
        if self.model.num_ligands() != 1
            || !self.model.flex.is_empty()
            || !self.model.glue_pairs.is_empty()
        {
            return Err("Metal currently supports one ligand and a rigid receptor, without macrocycle glue atoms".into());
        }
        if self.model.get_size().ligands[0] > 8 {
            return Err("Metal currently supports at most 8 torsions".into());
        }
        if self
            .model
            .atoms
            .iter()
            .any(|a| !a.base.atom_type.is_hydrogen() && a.base.atom_type.xs > 15)
        {
            return Err("Metal currently supports standard organic Vina atom types only".into());
        }
        let count = self.model.get_heavy_atom_movable_coords().len();
        if count == 0 || count > 64 {
            return Err("Metal currently supports 1 to 64 heavy atoms".into());
        }
        match self.grid.as_ref().unwrap() {
            DockingGrid::Vina(g) => Ok(g),
            _ => Err("Metal requires Vina affinity maps".into()),
        }
    }
    pub fn get_poses(&self, how_many: usize, energy_range: Fl) -> Result<String, String> {
        if how_many == 0 || !energy_range.is_finite() || energy_range < 0.0 {
            return Err(
                "Pose count must be positive and energy range finite and nonnegative".into(),
            );
        }
        let best = self.poses.first().ok_or("Docking has not been run")?.e;
        let mut model = self.model.clone();
        let mut text = String::new();
        for (i, pose) in self
            .poses
            .iter()
            .take(how_many)
            .take_while(|p| p.e <= best + energy_range)
            .enumerate()
        {
            model.set(&pose.c);
            let mut remark=format!("REMARK VINA RESULT: {:9.3}  {:9.3}  {:9.3}\nREMARK INTER + INTRA:    {:12.3}\nREMARK INTER:            {:12.3}\nREMARK INTRA:            {:12.3}\n",pose.e,pose.lb,pose.ub,pose.total,pose.inter,pose.intra);
            if self.choice == ScoringFunctionChoice::Ad42 {
                remark.push_str(&format!(
                    "REMARK CONF_INDEPENDENT: {:12.3}\n",
                    pose.conf_independent
                ));
            }
            remark.push_str(&format!("REMARK UNBOUND:          {:12.3}\n", pose.unbound));
            text.push_str(&model.write_model(i + 1, &remark));
        }
        Ok(text)
    }
    pub fn write_poses(&self, path: &str, how_many: usize, energy_range: Fl) -> Result<(), String> {
        std::fs::write(path, self.get_poses(how_many, energy_range)?)
            .map_err(|e| format!("{path}: {e}"))
    }
    pub fn get_poses_coordinates(
        &self,
        how_many: usize,
        energy_range: Fl,
    ) -> Result<Vec<Vec<Fl>>, String> {
        self.get_poses(how_many, energy_range)?;
        let mut model = self.model.clone();
        let best = self.poses[0].e;
        Ok(self
            .poses
            .iter()
            .take(how_many)
            .take_while(|p| p.e <= best + energy_range)
            .map(|pose| {
                model.set(&pose.c);
                model
                    .ligands
                    .iter()
                    .flat_map(|l| {
                        model.coords[l.range.begin..l.range.end]
                            .iter()
                            .flat_map(|v| v.data)
                    })
                    .collect()
            })
            .collect())
    }
    pub fn get_poses_energies(
        &self,
        how_many: usize,
        energy_range: Fl,
    ) -> Result<Vec<[Fl; 5]>, String> {
        self.get_poses(how_many, energy_range)?;
        let best = self.poses[0].e;
        Ok(self
            .poses
            .iter()
            .take(how_many)
            .take_while(|p| p.e <= best + energy_range)
            .map(|p| [p.e, p.inter, p.intra, p.conf_independent, p.unbound])
            .collect())
    }
    pub fn write_pose(&self, path: &str, remark: &str) -> Result<(), String> {
        self.ready()?;
        std::fs::write(path, self.model.write_model(1, remark)).map_err(|e| format!("{path}: {e}"))
    }
    pub fn grid_dimensions_from_ligand(&self, buffer: Fl) -> Result<(Vec3, Vec3), String> {
        if self.model.num_ligands() == 0 || !buffer.is_finite() || buffer <= 0.0 {
            return Err("Initialize ligands and use a finite positive buffer".into());
        }
        let center = self.model.center();
        let mut size = Vec3::new(0.0, 0.0, 0.0);
        for coord in &self.model.coords[..self.model.num_movable_atoms()] {
            for i in 0..3 {
                size[i] = size[i].max((coord[i] - center[i]).abs());
            }
        }
        for i in 0..3 {
            size[i] = (2.0 * (size[i] + buffer)).ceil();
        }
        Ok((center, size))
    }
}
fn optimize_model<G: IGrid<Model>>(
    model: &mut Model,
    p: &PrecalculateByAtom,
    grid: &G,
    out: &mut OutputType,
    steps: u32,
) {
    let mut adapter = SearchModel {
        model: model.clone(),
        precalculate: p,
        grid,
    };
    let mut gradient = Change::new(&model.get_size());
    QuasiNewton {
        max_steps: steps,
        average_required_improvement: 0.0,
    }
    .optimize(&mut adapter, out, &mut gradient, AUTHENTIC, &mut 0);
    *model = adapter.model;
}
