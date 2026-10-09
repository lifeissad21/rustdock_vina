use std::fs;
use std::io::Write;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_constants::*;
use crate::atom_type::{num_atom_types, AtomTyping};
use crate::common::{Fl, Vec3};
use crate::grid::Grid;
use crate::grid_dim::{GridDim, GridDims};
use crate::precalculate::Precalculate;
use crate::szv_grid::{szv_grid_dims, SzvGrid, SzvGridModel};

pub fn convert_xs_to_string(t: usize) -> &'static str {
    match t {
        XS_TYPE_C_H => "C_H",
        XS_TYPE_C_P => "C_P",
        XS_TYPE_N_P => "N_P",
        XS_TYPE_N_D => "N_D",
        XS_TYPE_N_A => "N_A",
        XS_TYPE_N_DA => "N_DA",
        XS_TYPE_O_P => "O_P",
        XS_TYPE_O_D => "O_D",
        XS_TYPE_O_A => "O_A",
        XS_TYPE_O_DA => "O_DA",
        XS_TYPE_S_P => "S_P",
        XS_TYPE_P_P => "P_P",
        XS_TYPE_F_H => "F_H",
        XS_TYPE_CL_H => "Cl_H",
        XS_TYPE_BR_H => "Br_H",
        XS_TYPE_I_H => "I_H",
        XS_TYPE_SI => "Si",
        XS_TYPE_AT => "At",
        XS_TYPE_MET_D => "Met_D",
        XS_TYPE_W => "W",
        _ => panic!("unsupported XS map type {t}"),
    }
}

pub fn canonical_xs_map_type(t: usize) -> Option<usize> {
    match t {
        XS_TYPE_G0 | XS_TYPE_G1 | XS_TYPE_G2 | XS_TYPE_G3 => None,
        XS_TYPE_C_H_CG0 | XS_TYPE_C_H_CG1 | XS_TYPE_C_H_CG2 | XS_TYPE_C_H_CG3 => Some(XS_TYPE_C_H),
        XS_TYPE_C_P_CG0 | XS_TYPE_C_P_CG1 | XS_TYPE_C_P_CG2 | XS_TYPE_C_P_CG3 => Some(XS_TYPE_C_P),
        t if t < XS_TYPE_SIZE => Some(t),
        _ => None,
    }
}

pub trait CacheModel: SzvGridModel {
    fn num_movable_atoms(&self) -> usize;
    fn atoms(&self) -> &[Atom];
    fn coords(&self) -> &[Vec3];
    fn minus_forces_mut(&mut self) -> &mut [Vec3];
    fn is_atom_in_ligand(&self, atom_index: usize) -> bool;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cache {
    gd: GridDims,
    slope: Fl,
    grids: Vec<Grid>,
}

impl Default for Cache {
    fn default() -> Self {
        Self::new(1e6)
    }
}

impl Cache {
    pub fn new(slope: Fl) -> Self {
        Self {
            gd: [
                GridDim::disabled(),
                GridDim::disabled(),
                GridDim::disabled(),
            ],
            slope,
            grids: vec![Grid::default(); XS_TYPE_SIZE],
        }
    }

    pub fn with_dims(gd: GridDims, slope: Fl) -> Self {
        Self {
            gd,
            slope,
            grids: vec![Grid::default(); XS_TYPE_SIZE],
        }
    }

    pub fn grid_for_type(&self, atom_type: usize) -> Option<&Grid> {
        self.grids.get(atom_type).filter(|g| g.initialized())
    }
    pub fn gd(&self) -> GridDims {
        self.gd
    }

    pub fn corner1(&self) -> Vec3 {
        Vec3::new(self.gd[0].begin, self.gd[1].begin, self.gd[2].begin)
    }

    pub fn corner2(&self) -> Vec3 {
        Vec3::new(self.gd[0].end, self.gd[1].end, self.gd[2].end)
    }

    pub fn is_atom_type_grid_initialized(&self, t: usize) -> bool {
        self.grids[t].initialized()
    }

    pub fn eval<M: CacheModel>(&self, model: &M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            let atom = &model.atoms()[i];
            let Some(t) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                continue;
            };
            if t >= nat {
                continue;
            }
            e += self.grids[t].evaluate(model.coords()[i], self.slope, v);
        }
        e
    }

    pub fn eval_intra<M: CacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            if model.is_atom_in_ligand(i) {
                continue;
            }
            let atom = &model.atoms()[i];
            let Some(t) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                continue;
            };
            if t >= nat {
                continue;
            }
            e += self.grids[t].evaluate(model.coords()[i], self.slope, v);
        }
        e
    }

    pub fn eval_deriv<M: CacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            let atom = &model.atoms()[i];
            let Some(t) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                model.minus_forces_mut()[i] = Vec3::new(0.0, 0.0, 0.0);
                continue;
            };
            if t >= nat {
                model.minus_forces_mut()[i] = Vec3::new(0.0, 0.0, 0.0);
                continue;
            }
            let mut deriv = Vec3::new(0.0, 0.0, 0.0);
            e += self.grids[t].evaluate_with_deriv(model.coords()[i], self.slope, v, &mut deriv);
            model.minus_forces_mut()[i] = deriv;
        }
        e
    }

    pub fn is_in_grid<M: CacheModel>(&self, model: &M, margin: Fl) -> bool {
        for i in 0..model.num_movable_atoms() {
            if model.atoms()[i].base.atom_type.is_hydrogen() {
                continue;
            }
            let coords = model.coords()[i];
            for j in 0..3 {
                if self.gd[j].n_voxels > 0
                    && (coords[j] < self.gd[j].begin - margin
                        || coords[j] > self.gd[j].end + margin)
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn are_atom_types_grid_initialized(&self, atom_types: &[usize]) -> bool {
        let nat = num_atom_types(AtomTyping::Xs);
        for t in atom_types {
            let Some(t) = canonical_xs_map_type(*t) else {
                continue;
            };
            if t >= nat {
                continue;
            }
            if !self.is_atom_type_grid_initialized(t) {
                return false;
            }
        }
        true
    }

    pub fn read(&mut self, map_prefix: &str) -> Result<(), String> {
        let mut gds = Vec::new();
        let mut found_at_least_one_map = false;
        let mut got_c_h = false;
        let mut got_c_p = false;

        for atom_type in 0..XS_TYPE_SIZE {
            let mut t = atom_type;
            match t {
                XS_TYPE_G0 | XS_TYPE_G1 | XS_TYPE_G2 | XS_TYPE_G3 => continue,
                XS_TYPE_C_H_CG0 | XS_TYPE_C_H_CG1 | XS_TYPE_C_H_CG2 | XS_TYPE_C_H_CG3 => {
                    if got_c_h {
                        continue;
                    }
                    t = XS_TYPE_C_H;
                    got_c_h = true;
                }
                XS_TYPE_C_P_CG0 | XS_TYPE_C_P_CG1 | XS_TYPE_C_P_CG2 | XS_TYPE_C_P_CG3 => {
                    if got_c_p {
                        continue;
                    }
                    t = XS_TYPE_C_P;
                    got_c_p = true;
                }
                _ => {}
            }
            let filename = format!("{map_prefix}.{}.map", convert_xs_to_string(t));
            if Path::new(&filename).exists() {
                read_vina_map(&filename, &mut gds, &mut self.grids[t])?;
                found_at_least_one_map = true;
            }
        }
        if !found_at_least_one_map {
            return Err(format!("No *.map files with prefix \"{map_prefix}\""));
        }
        self.gd = gds[0];
        Ok(())
    }

    pub fn write(
        &self,
        out_prefix: &str,
        atom_types: &[usize],
        gpf_filename: &str,
        fld_filename: &str,
        receptor_filename: &str,
    ) -> Result<(), String> {
        let mut got_c_h = false;
        let mut got_c_p = false;
        for t in atom_types {
            let Some(t) = canonical_xs_map_type(*t) else {
                continue;
            };
            if t == XS_TYPE_C_H {
                if got_c_h {
                    continue;
                }
                got_c_h = true;
            }
            if t == XS_TYPE_C_P {
                if got_c_p {
                    continue;
                }
                got_c_p = true;
            }
            if !self.grids[t].initialized() {
                continue;
            }
            let grid = &self.grids[t];
            let nx = grid.data.dim0() - 1;
            let ny = grid.data.dim1() - 1;
            let nz = grid.data.dim2() - 1;
            if nx % 2 == 1 || ny % 2 == 1 || nz % 2 == 1 {
                return Err("Can't write maps. Number of voxels (NELEMENTS) is odd. Use --force_even_voxels.".to_string());
            }
            let filename = format!("{out_prefix}.{}.map", convert_xs_to_string(t));
            let mut out = fs::File::create(&filename).map_err(|err| err.to_string())?;
            writeln!(out, "GRID_PARAMETER_FILE {gpf_filename}").map_err(|err| err.to_string())?;
            writeln!(out, "GRID_DATA_FILE {fld_filename}").map_err(|err| err.to_string())?;
            writeln!(out, "MACROMOLECULE {receptor_filename}").map_err(|err| err.to_string())?;
            writeln!(out, "SPACING {}", grid.factor_inv[0]).map_err(|err| err.to_string())?;
            writeln!(out, "NELEMENTS {nx} {ny} {nz}").map_err(|err| err.to_string())?;
            let cx = grid.init[0] + grid.range[0] * 0.5;
            let cy = grid.init[1] + grid.range[1] * 0.5;
            let cz = grid.init[2] + grid.range[2] * 0.5;
            writeln!(out, "CENTER {cx} {cy} {cz}").map_err(|err| err.to_string())?;
            for z in 0..grid.data.dim2() {
                for y in 0..grid.data.dim1() {
                    for x in 0..grid.data.dim0() {
                        writeln!(out, "{:.4}", grid.data.get(x, y, z))
                            .map_err(|err| err.to_string())?;
                    }
                }
            }
        }
        Ok(())
    }

    pub fn populate<M: CacheModel>(
        &mut self,
        model: &M,
        precalculate: &Precalculate,
        atom_types_needed: &[usize],
    ) {
        let mut needed = Vec::new();
        let mut got_c_h = false;
        let mut got_c_p = false;
        for t in atom_types_needed {
            let Some(t) = canonical_xs_map_type(*t) else {
                continue;
            };
            if t == XS_TYPE_C_H {
                if got_c_h {
                    continue;
                }
                got_c_h = true;
            }
            if t == XS_TYPE_C_P {
                if got_c_p {
                    continue;
                }
                got_c_p = true;
            }
            if !self.grids[t].initialized() {
                needed.push(t);
                self.grids[t].init(&self.gd);
            }
        }
        if needed.is_empty() {
            return;
        }

        let nat = num_atom_types(AtomTyping::Xs);
        let reduced = szv_grid_dims(&self.gd);
        let ig = SzvGrid::new(model, &reduced, precalculate.cutoff_sqr());
        let dimensions = (
            self.grids[needed[0]].data.dim0(),
            self.grids[needed[0]].data.dim1(),
            self.grids[needed[0]].data.dim2(),
        );

        for x in 0..dimensions.0 {
            for y in 0..dimensions.1 {
                for z in 0..dimensions.2 {
                    let mut affinities = vec![0.0; needed.len()];
                    let probe_coords = self.grids[needed[0]].index_to_argument(x, y, z);
                    for i in ig.possibilities(probe_coords) {
                        let atom = &model.grid_atoms()[*i];
                        let t1 = atom.base.atom_type.xs;
                        if t1 >= nat {
                            continue;
                        }
                        let r2 = crate::common::vec_distance_sqr(atom.coords, probe_coords);
                        if r2 <= precalculate.cutoff_sqr() {
                            for (j, t2) in needed.iter().enumerate() {
                                let type_pair_index = crate::triangular_matrix_index::triangular_matrix_index_permissive(nat, t1, *t2);
                                affinities[j] += precalculate.eval_fast(type_pair_index, r2);
                            }
                        }
                    }
                    for (j, t) in needed.iter().enumerate() {
                        *self.grids[*t].data.get_mut(x, y, z) = affinities[j];
                    }
                }
            }
        }
    }
}

pub fn vina_split_fields(value: &str) -> Vec<&str> {
    value.split(' ').collect()
}

pub fn read_vina_map(
    filename: &str,
    gds: &mut Vec<GridDims>,
    grid: &mut Grid,
) -> Result<(), String> {
    let input = fs::read_to_string(filename).map_err(|err| format!("{filename}: {err}"))?;
    let lines: Vec<_> = input.lines().collect();
    if lines.len() < 6 {
        return Err(format!("{filename}: incomplete map header"));
    }
    let fields = |index: usize, label: &str, count: usize| -> Result<Vec<&str>, String> {
        let f: Vec<_> = lines[index].split_whitespace().collect();
        if f.len() != count || f[0] != label {
            return Err(format!("{filename}: invalid {label} header"));
        }
        Ok(f)
    };
    let spacing = fields(3, "SPACING", 2)?[1]
        .parse::<Fl>()
        .map_err(|e| e.to_string())?;
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err("Map spacing must be finite and positive".into());
    }
    let counts = fields(4, "NELEMENTS", 4)?;
    let centers = fields(5, "CENTER", 4)?;
    let mut gd = [GridDim::disabled(); 3];
    let mut points = 1usize;
    for i in 0..3 {
        let n = counts[i + 1].parse::<usize>().map_err(|e| e.to_string())?;
        if n == 0 || n % 2 != 0 {
            return Err("Map NELEMENTS must be positive and even".into());
        }
        let center = centers[i + 1].parse::<Fl>().map_err(|e| e.to_string())?;
        let halfspan = n as Fl * spacing / 2.0;
        if !center.is_finite() || !halfspan.is_finite() {
            return Err("Nonfinite map dimensions".into());
        }
        gd[i] = GridDim::new(center - halfspan, center + halfspan, n);
        points = points
            .checked_mul(n.checked_add(1).ok_or("Map size overflow")?)
            .ok_or("Map size overflow")?;
    }
    if lines.len() - 6 != points {
        return Err(format!(
            "{filename}: expected {points} grid values, found {}",
            lines.len() - 6
        ));
    }
    if gds.first().is_some_and(|previous| *previous != gd) {
        return Err(format!("{filename}: inconsistent map dimensions"));
    }
    let mut parsed = Grid::default();
    parsed.init(&gd);
    let nx = gd[0].n_voxels + 1;
    let ny = gd[1].n_voxels + 1;
    for (i, line) in lines[6..].iter().enumerate() {
        let value = line
            .trim()
            .parse::<Fl>()
            .map_err(|e| format!("{filename}: {e}"))?;
        if !value.is_finite() {
            return Err(format!("{filename}: nonfinite map value"));
        }
        *parsed.data.get_mut(i % nx, (i / nx) % ny, i / (nx * ny)) = value;
    }
    gds.push(gd);
    *grid = parsed;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_constants::{AD_TYPE_C, EL_TYPE_C};
    use crate::atom_type::AtomType;

    struct FakeModel {
        atoms: Vec<Atom>,
        coords: Vec<Vec3>,
        minus_forces: Vec<Vec3>,
    }

    impl SzvGridModel for FakeModel {
        fn atom_typing_used(&self) -> AtomTyping {
            AtomTyping::Xs
        }

        fn grid_atoms(&self) -> &[Atom] {
            &self.atoms
        }
    }

    impl CacheModel for FakeModel {
        fn num_movable_atoms(&self) -> usize {
            self.atoms.len()
        }

        fn atoms(&self) -> &[Atom] {
            &self.atoms
        }

        fn coords(&self) -> &[Vec3] {
            &self.coords
        }

        fn minus_forces_mut(&mut self) -> &mut [Vec3] {
            &mut self.minus_forces
        }

        fn is_atom_in_ligand(&self, _atom_index: usize) -> bool {
            false
        }
    }

    #[test]
    fn canonicalizes_macrocycle_types() {
        assert_eq!(canonical_xs_map_type(XS_TYPE_C_H_CG0), Some(XS_TYPE_C_H));
        assert_eq!(canonical_xs_map_type(XS_TYPE_G0), None);
        assert_eq!(convert_xs_to_string(XS_TYPE_C_H), "C_H");
    }

    #[test]
    fn evaluates_initialized_grid() {
        let dims = [
            GridDim::new(0.0, 2.0, 2),
            GridDim::new(0.0, 2.0, 2),
            GridDim::new(0.0, 2.0, 2),
        ];
        let mut cache = Cache::with_dims(dims, 0.0);
        cache.grids[XS_TYPE_C_H].init(&dims);
        for x in 0..cache.grids[XS_TYPE_C_H].data.dim0() {
            for y in 0..cache.grids[XS_TYPE_C_H].data.dim1() {
                for z in 0..cache.grids[XS_TYPE_C_H].data.dim2() {
                    *cache.grids[XS_TYPE_C_H].data.get_mut(x, y, z) = 1.0;
                }
            }
        }
        let mut atom = Atom::default();
        atom.base.atom_type = AtomType {
            el: EL_TYPE_C,
            ad: AD_TYPE_C,
            xs: XS_TYPE_C_H,
            sy: 0,
        };
        let model = FakeModel {
            atoms: vec![atom],
            coords: vec![Vec3::new(1.0, 1.0, 1.0)],
            minus_forces: vec![Vec3::new(0.0, 0.0, 0.0)],
        };
        assert_eq!(cache.eval(&model, crate::common::MAX_FL), 1.0);
        assert!(cache.is_in_grid(&model, 0.0001));
    }
}
