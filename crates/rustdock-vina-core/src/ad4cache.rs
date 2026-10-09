use std::fs;
use std::io::Write;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_constants::*;
use crate::cache::read_vina_map;
use crate::common::{Fl, Vec3, ZERO_VEC};
use crate::grid::Grid;
use crate::grid_dim::{GridDim, GridDims};
use crate::szv_grid::SzvGridModel;

pub const AD4_ELECTROSTATIC_MAP: usize = AD_TYPE_SIZE;
pub const AD4_DESOLVATION_MAP: usize = AD_TYPE_SIZE + 1;
pub const AD4_NUM_MAPS: usize = AD_TYPE_SIZE + 2;

pub fn convert_ad_to_string(t: usize) -> &'static str {
    match t {
        AD_TYPE_C => "C",
        AD_TYPE_A => "A",
        AD_TYPE_N => "N",
        AD_TYPE_O => "O",
        AD_TYPE_P => "P",
        AD_TYPE_S => "S",
        AD_TYPE_H => "H",
        AD_TYPE_F => "F",
        AD_TYPE_I => "I",
        AD_TYPE_NA => "NA",
        AD_TYPE_OA => "OA",
        AD_TYPE_SA => "SA",
        AD_TYPE_HD => "HD",
        AD_TYPE_MG => "Mg",
        AD_TYPE_MN => "Mn",
        AD_TYPE_ZN => "Zn",
        AD_TYPE_CA => "Ca",
        AD_TYPE_FE => "Fe",
        AD_TYPE_CL => "Cl",
        AD_TYPE_BR => "Br",
        AD_TYPE_SI => "Si",
        AD_TYPE_AT => "At",
        AD_TYPE_W => "W",
        AD4_ELECTROSTATIC_MAP => "e",
        AD4_DESOLVATION_MAP => "d",
        _ => panic!("unsupported AD4 map type {t}"),
    }
}

pub fn canonical_ad4_map_type(t: usize) -> Option<usize> {
    match t {
        AD_TYPE_G0 | AD_TYPE_G1 | AD_TYPE_G2 | AD_TYPE_G3 => None,
        AD_TYPE_CG0 | AD_TYPE_CG1 | AD_TYPE_CG2 | AD_TYPE_CG3 => Some(AD_TYPE_C),
        t if t < AD_TYPE_SIZE => Some(t),
        _ => None,
    }
}

pub trait Ad4CacheModel: SzvGridModel {
    fn num_movable_atoms(&self) -> usize;
    fn atoms(&self) -> &[Atom];
    fn coords(&self) -> &[Vec3];
    fn minus_forces_mut(&mut self) -> &mut [Vec3];
    fn is_atom_in_ligand(&self, atom_index: usize) -> bool;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Ad4Cache {
    gd: GridDims,
    slope: Fl,
    grids: Vec<Grid>,
}

impl Default for Ad4Cache {
    fn default() -> Self {
        Self::new(1e6)
    }
}

impl Ad4Cache {
    pub fn new(slope: Fl) -> Self {
        Self {
            gd: [
                GridDim::disabled(),
                GridDim::disabled(),
                GridDim::disabled(),
            ],
            slope,
            grids: vec![Grid::default(); AD4_NUM_MAPS],
        }
    }

    pub fn with_dims(gd: GridDims, slope: Fl) -> Self {
        Self {
            gd,
            slope,
            grids: vec![Grid::default(); AD4_NUM_MAPS],
        }
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

    pub fn eval<M: Ad4CacheModel>(&self, model: &M, v: Fl) -> Fl {
        let mut e = 0.0;
        for i in 0..model.num_movable_atoms() {
            if !model.is_atom_in_ligand(i) {
                continue;
            }
            e += self.eval_atom(&model.atoms()[i], model.coords()[i], v);
        }
        e
    }

    pub fn eval_intra<M: Ad4CacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        for i in 0..model.num_movable_atoms() {
            if model.is_atom_in_ligand(i) {
                continue;
            }
            e += self.eval_atom(&model.atoms()[i], model.coords()[i], v);
        }
        e
    }

    pub fn eval_deriv<M: Ad4CacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        for i in 0..model.num_movable_atoms() {
            let atom = model.atoms()[i].clone();
            let coords = model.coords()[i];
            let Some(t) = canonical_ad4_map_type(atom.base.atom_type.ad) else {
                model.minus_forces_mut()[i] = ZERO_VEC;
                continue;
            };

            let mut minus_force = ZERO_VEC;

            let mut deriv = ZERO_VEC;
            e += self.grids[t].evaluate_with_deriv(coords, self.slope, v, &mut deriv);
            minus_force += deriv;

            e += self.grids[AD4_ELECTROSTATIC_MAP]
                .evaluate_with_deriv(coords, self.slope, v, &mut deriv)
                * atom.base.charge;
            minus_force += atom.base.charge * deriv;

            let charge_abs = atom.base.charge.abs();
            e += self.grids[AD4_DESOLVATION_MAP]
                .evaluate_with_deriv(coords, self.slope, v, &mut deriv)
                * charge_abs;
            minus_force += charge_abs * deriv;

            model.minus_forces_mut()[i] = minus_force;
        }
        e
    }

    pub fn is_in_grid<M: Ad4CacheModel>(&self, model: &M, margin: Fl) -> bool {
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
        for t in atom_types {
            let Some(t) = canonical_ad4_map_type(*t) else {
                continue;
            };
            if !self.is_atom_type_grid_initialized(t) {
                return false;
            }
        }
        self.is_atom_type_grid_initialized(AD4_ELECTROSTATIC_MAP)
            && self.is_atom_type_grid_initialized(AD4_DESOLVATION_MAP)
    }

    pub fn read(&mut self, map_prefix: &str) -> Result<(), String> {
        let mut gds = Vec::new();
        let mut got_c = false;

        for atom_type in 0..AD_TYPE_SIZE {
            let Some(t) = canonical_ad4_map_type(atom_type) else {
                continue;
            };
            if t == AD_TYPE_C {
                if got_c {
                    continue;
                }
                got_c = true;
            }
            let filename = format!("{map_prefix}.{}.map", convert_ad_to_string(t));
            if Path::new(&filename).exists() {
                read_vina_map(&filename, &mut gds, &mut self.grids[t])?;
            }
        }

        read_vina_map(
            &format!("{map_prefix}.e.map"),
            &mut gds,
            &mut self.grids[AD4_ELECTROSTATIC_MAP],
        )?;
        read_vina_map(
            &format!("{map_prefix}.d.map"),
            &mut gds,
            &mut self.grids[AD4_DESOLVATION_MAP],
        )?;

        if gds.is_empty() {
            return Err(format!("No AD4 map files with prefix \"{map_prefix}\""));
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
        let mut got_c = false;
        for t in atom_types {
            let Some(t) = canonical_ad4_map_type(*t) else {
                continue;
            };
            if t == AD_TYPE_C {
                if got_c {
                    continue;
                }
                got_c = true;
            }
            self.write_one_map(t, out_prefix, gpf_filename, fld_filename, receptor_filename)?;
        }
        self.write_one_map(
            AD4_ELECTROSTATIC_MAP,
            out_prefix,
            gpf_filename,
            fld_filename,
            receptor_filename,
        )?;
        self.write_one_map(
            AD4_DESOLVATION_MAP,
            out_prefix,
            gpf_filename,
            fld_filename,
            receptor_filename,
        )?;
        Ok(())
    }

    fn eval_atom(&self, atom: &Atom, coords: Vec3, v: Fl) -> Fl {
        let Some(t) = canonical_ad4_map_type(atom.base.atom_type.ad) else {
            return 0.0;
        };
        self.grids[t].evaluate(coords, self.slope, v)
            + self.grids[AD4_ELECTROSTATIC_MAP].evaluate(coords, self.slope, v) * atom.base.charge
            + self.grids[AD4_DESOLVATION_MAP].evaluate(coords, self.slope, v)
                * atom.base.charge.abs()
    }

    fn write_one_map(
        &self,
        t: usize,
        out_prefix: &str,
        gpf_filename: &str,
        fld_filename: &str,
        receptor_filename: &str,
    ) -> Result<(), String> {
        if !self.grids[t].initialized() {
            return Ok(());
        }
        let grid = &self.grids[t];
        let size_x = if grid.data.dim0() % 2 == 0 {
            grid.data.dim0()
        } else {
            grid.data.dim0() - 1
        };
        let size_y = if grid.data.dim1() % 2 == 0 {
            grid.data.dim1()
        } else {
            grid.data.dim1() - 1
        };
        let size_z = if grid.data.dim2() % 2 == 0 {
            grid.data.dim2()
        } else {
            grid.data.dim2() - 1
        };
        let filename = format!("{out_prefix}.{}.map", convert_ad_to_string(t));
        let mut out = fs::File::create(&filename).map_err(|err| err.to_string())?;
        writeln!(out, "GRID_PARAMETER_FILE {gpf_filename}").map_err(|err| err.to_string())?;
        writeln!(out, "GRID_DATA_FILE {fld_filename}").map_err(|err| err.to_string())?;
        writeln!(out, "MACROMOLECULE {receptor_filename}").map_err(|err| err.to_string())?;
        writeln!(out, "SPACING {}", grid.factor_inv[0]).map_err(|err| err.to_string())?;
        writeln!(out, "NELEMENTS {size_x} {size_y} {size_z}").map_err(|err| err.to_string())?;
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
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_type::AtomType;
    use crate::atom_type::AtomTyping;
    use crate::common::{eq_fl, eq_vec, MAX_FL};
    use crate::grid_dim::GridDim;

    struct FakeModel {
        atoms: Vec<Atom>,
        coords: Vec<Vec3>,
        minus_forces: Vec<Vec3>,
        in_ligand: Vec<bool>,
    }

    impl SzvGridModel for FakeModel {
        fn atom_typing_used(&self) -> AtomTyping {
            AtomTyping::Ad
        }

        fn grid_atoms(&self) -> &[Atom] {
            &self.atoms
        }
    }

    impl Ad4CacheModel for FakeModel {
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

        fn is_atom_in_ligand(&self, atom_index: usize) -> bool {
            self.in_ligand[atom_index]
        }
    }

    fn carbon(charge: Fl) -> Atom {
        let mut atom = Atom::default();
        atom.base.atom_type = AtomType {
            el: EL_TYPE_C,
            ad: AD_TYPE_C,
            xs: XS_TYPE_C_H,
            sy: 0,
        };
        atom.base.charge = charge;
        atom
    }

    fn filled_cache() -> Ad4Cache {
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let mut cache = Ad4Cache::with_dims(dims, 0.0);
        for (t, value) in [
            (AD_TYPE_C, 1.0),
            (AD4_ELECTROSTATIC_MAP, 2.0),
            (AD4_DESOLVATION_MAP, 3.0),
        ] {
            cache.grids[t].init(&dims);
            for x in 0..cache.grids[t].data.dim0() {
                for y in 0..cache.grids[t].data.dim1() {
                    for z in 0..cache.grids[t].data.dim2() {
                        *cache.grids[t].data.get_mut(x, y, z) = value;
                    }
                }
            }
        }
        cache
    }

    #[test]
    fn canonicalizes_ad4_types() {
        assert_eq!(canonical_ad4_map_type(AD_TYPE_CG0), Some(AD_TYPE_C));
        assert_eq!(canonical_ad4_map_type(AD_TYPE_G0), None);
        assert_eq!(convert_ad_to_string(AD_TYPE_HD), "HD");
        assert_eq!(convert_ad_to_string(AD4_ELECTROSTATIC_MAP), "e");
    }

    #[test]
    fn evaluates_ligand_maps_with_charge_terms() {
        let cache = filled_cache();
        let atom = carbon(-0.5);
        let model = FakeModel {
            atoms: vec![atom],
            coords: vec![Vec3::new(0.5, 0.5, 0.5)],
            minus_forces: vec![ZERO_VEC],
            in_ligand: vec![true],
        };
        assert!(eq_fl(cache.eval(&model, MAX_FL), 1.0 - 1.0 + 1.5));
        assert!(cache.is_in_grid(&model, 0.0001));
    }

    #[test]
    fn eval_intra_skips_ligand_atoms() {
        let cache = filled_cache();
        let atom = carbon(1.0);
        let mut model = FakeModel {
            atoms: vec![atom],
            coords: vec![Vec3::new(0.5, 0.5, 0.5)],
            minus_forces: vec![ZERO_VEC],
            in_ligand: vec![true],
        };
        assert!(eq_fl(cache.eval_intra(&mut model, MAX_FL), 0.0));
    }

    #[test]
    fn eval_deriv_combines_map_derivatives() {
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let mut cache = Ad4Cache::with_dims(dims, 0.0);
        for t in [AD_TYPE_C, AD4_ELECTROSTATIC_MAP, AD4_DESOLVATION_MAP] {
            cache.grids[t].init(&dims);
            *cache.grids[t].data.get_mut(0, 0, 0) = 0.0;
            *cache.grids[t].data.get_mut(1, 0, 0) = 1.0;
            *cache.grids[t].data.get_mut(0, 1, 0) = 0.0;
            *cache.grids[t].data.get_mut(1, 1, 0) = 1.0;
            *cache.grids[t].data.get_mut(0, 0, 1) = 0.0;
            *cache.grids[t].data.get_mut(1, 0, 1) = 1.0;
            *cache.grids[t].data.get_mut(0, 1, 1) = 0.0;
            *cache.grids[t].data.get_mut(1, 1, 1) = 1.0;
        }
        let atom = carbon(0.5);
        let mut model = FakeModel {
            atoms: vec![atom],
            coords: vec![Vec3::new(0.5, 0.5, 0.5)],
            minus_forces: vec![ZERO_VEC],
            in_ligand: vec![true],
        };
        let e = cache.eval_deriv(&mut model, MAX_FL);
        assert!(eq_fl(e, 1.0));
        assert!(eq_vec(model.minus_forces[0], Vec3::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn required_maps_include_electrostatic_and_desolvation() {
        let mut cache = filled_cache();
        assert!(cache.are_atom_types_grid_initialized(&[AD_TYPE_CG0]));
        cache.grids[AD4_DESOLVATION_MAP] = Grid::default();
        assert!(!cache.are_atom_types_grid_initialized(&[AD_TYPE_C]));
    }
}
