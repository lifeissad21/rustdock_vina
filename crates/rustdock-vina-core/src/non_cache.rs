use crate::atom::Atom;
use crate::atom_type::{get_type_pair_index, num_atom_types, AtomTyping};
use crate::cache::canonical_xs_map_type;
use crate::common::{Fl, Vec3, ZERO_VEC};
use crate::curl::{curl, curl_with_derivative};
use crate::grid_dim::GridDims;
use crate::precalculate::Precalculate;
use crate::szv_grid::{szv_grid_dims, SzvGrid, SzvGridModel};

pub trait NonCacheModel: SzvGridModel {
    fn num_movable_atoms(&self) -> usize;
    fn atoms(&self) -> &[Atom];
    fn coords(&self) -> &[Vec3];
    fn minus_forces_mut(&mut self) -> &mut [Vec3];
    fn is_atom_in_ligand(&self, atom_index: usize) -> bool;
}

#[derive(Debug, Clone, PartialEq)]
pub struct NonCache {
    sgrid: SzvGrid,
    gd: GridDims,
    p: Precalculate,
    pub slope: Fl,
}

impl NonCache {
    pub fn new<M: NonCacheModel>(model: &M, gd: GridDims, p: Precalculate, slope: Fl) -> Self {
        Self {
            sgrid: SzvGrid::new(model, &szv_grid_dims(&gd), p.cutoff_sqr()),
            gd,
            p,
            slope,
        }
    }

    pub fn eval<M: NonCacheModel>(&self, model: &M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            let mut this_e = 0.0;
            let mut out_of_bounds_penalty = 0.0;
            let atom = &model.atoms()[i];
            let Some(t1) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                continue;
            };
            if t1 >= nat {
                continue;
            }
            let a_coords = model.coords()[i];
            let adjusted = self.adjusted_coords(a_coords, &mut out_of_bounds_penalty, None);
            out_of_bounds_penalty *= self.slope;
            this_e += self.eval_adjusted_atom(atom, adjusted, model);
            curl(&mut this_e, v);
            e += this_e + out_of_bounds_penalty;
        }
        e
    }

    pub fn eval_intra<M: NonCacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            if model.is_atom_in_ligand(i) {
                continue;
            }
            let mut this_e = 0.0;
            let mut out_of_bounds_penalty = 0.0;
            let atom = &model.atoms()[i];
            let Some(t1) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                continue;
            };
            if t1 >= nat {
                continue;
            }
            let adjusted =
                self.adjusted_coords(model.coords()[i], &mut out_of_bounds_penalty, None);
            out_of_bounds_penalty *= self.slope;
            this_e += self.eval_adjusted_atom(atom, adjusted, model);
            curl(&mut this_e, v);
            e += this_e + out_of_bounds_penalty;
        }
        e
    }

    pub fn eval_deriv<M: NonCacheModel>(&self, model: &mut M, v: Fl) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for i in 0..model.num_movable_atoms() {
            let atom = &model.atoms()[i];
            let Some(t1) = canonical_xs_map_type(atom.base.atom_type.xs) else {
                model.minus_forces_mut()[i] = ZERO_VEC;
                continue;
            };
            if t1 >= nat {
                model.minus_forces_mut()[i] = ZERO_VEC;
                continue;
            }

            let mut this_e = 0.0;
            let mut deriv = ZERO_VEC;
            let mut out_of_bounds_deriv = ZERO_VEC;
            let mut out_of_bounds_penalty = 0.0;
            let adjusted = self.adjusted_coords(
                model.coords()[i],
                &mut out_of_bounds_penalty,
                Some(&mut out_of_bounds_deriv),
            );
            out_of_bounds_penalty *= self.slope;
            out_of_bounds_deriv *= self.slope;

            for j in self.sgrid.possibilities(adjusted) {
                let b = &model.grid_atoms()[*j];
                if canonical_xs_map_type(b.base.atom_type.xs).is_none_or(|t2| t2 >= nat) {
                    continue;
                }
                let r_ba = adjusted - b.coords;
                let r2 = r_ba.norm_sqr();
                if r2 < self.p.cutoff_sqr() {
                    let type_pair_index =
                        get_type_pair_index(AtomTyping::Xs, atom.base.atom_type, b.base.atom_type);
                    let (pair_e, dor) = self.p.eval_deriv(type_pair_index, r2);
                    this_e += pair_e;
                    deriv += dor * r_ba;
                }
            }
            curl_with_derivative(&mut this_e, &mut deriv, v);
            model.minus_forces_mut()[i] = deriv + out_of_bounds_deriv;
            e += this_e + out_of_bounds_penalty;
        }
        e
    }

    pub fn within<M: NonCacheModel>(&self, model: &M, margin: Fl) -> bool {
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

    fn adjusted_coords(
        &self,
        coords: Vec3,
        out_of_bounds_penalty: &mut Fl,
        mut out_of_bounds_deriv: Option<&mut Vec3>,
    ) -> Vec3 {
        let mut adjusted = coords;
        for j in 0..3 {
            if self.gd[j].n_voxels == 0 {
                continue;
            }
            if coords[j] < self.gd[j].begin {
                adjusted[j] = self.gd[j].begin;
                if let Some(deriv) = out_of_bounds_deriv.as_deref_mut() {
                    deriv[j] = -1.0;
                }
                *out_of_bounds_penalty += (coords[j] - self.gd[j].begin).abs();
            } else if coords[j] > self.gd[j].end {
                adjusted[j] = self.gd[j].end;
                if let Some(deriv) = out_of_bounds_deriv.as_deref_mut() {
                    deriv[j] = 1.0;
                }
                *out_of_bounds_penalty += (coords[j] - self.gd[j].end).abs();
            }
        }
        adjusted
    }

    fn eval_adjusted_atom<M: NonCacheModel>(&self, atom: &Atom, adjusted: Vec3, model: &M) -> Fl {
        let mut e = 0.0;
        let nat = num_atom_types(AtomTyping::Xs);
        for j in self.sgrid.possibilities(adjusted) {
            let b = &model.grid_atoms()[*j];
            if canonical_xs_map_type(b.base.atom_type.xs).is_none_or(|t2| t2 >= nat) {
                continue;
            }
            let r_ba = adjusted - b.coords;
            let r2 = r_ba.norm_sqr();
            if r2 < self.p.cutoff_sqr() {
                let type_pair_index =
                    get_type_pair_index(AtomTyping::Xs, atom.base.atom_type, b.base.atom_type);
                e += self.p.eval_fast(type_pair_index, r2);
            }
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_constants::{AD_TYPE_C, EL_TYPE_C, XS_TYPE_C_H};
    use crate::atom_type::AtomType;
    use crate::common::{eq_fl, eq_vec};
    use crate::grid_dim::GridDim;
    use crate::precalculate::Precalculate;
    use crate::scoring_function::{ScoringFunction, ScoringFunctionChoice};

    #[derive(Clone)]
    struct FakeModel {
        atoms: Vec<Atom>,
        coords: Vec<Vec3>,
        grid_atoms: Vec<Atom>,
        minus_forces: Vec<Vec3>,
        atom_in_ligand: Vec<bool>,
    }

    impl SzvGridModel for FakeModel {
        fn atom_typing_used(&self) -> AtomTyping {
            AtomTyping::Xs
        }

        fn grid_atoms(&self) -> &[Atom] {
            &self.grid_atoms
        }
    }

    impl NonCacheModel for FakeModel {
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
            self.atom_in_ligand[atom_index]
        }
    }

    fn carbon_atom(coords: Vec3) -> Atom {
        let mut atom = Atom {
            coords,
            ..Atom::default()
        };
        atom.base.atom_type = AtomType {
            el: EL_TYPE_C,
            ad: AD_TYPE_C,
            xs: XS_TYPE_C_H,
            sy: 0,
        };
        atom
    }

    fn precalculate() -> Precalculate {
        let sf = ScoringFunction::new(ScoringFunctionChoice::Vina, vec![0.0; 6]);
        Precalculate::new(&sf, 1000.0, 16.0)
    }

    #[test]
    fn within_rejects_heavy_atom_outside_grid() {
        let atom = carbon_atom(Vec3::new(2.0, 0.5, 0.5));
        let model = FakeModel {
            atoms: vec![atom.clone()],
            coords: vec![atom.coords],
            grid_atoms: vec![carbon_atom(Vec3::new(0.5, 0.5, 0.5))],
            minus_forces: vec![ZERO_VEC],
            atom_in_ligand: vec![false],
        };
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let non_cache = NonCache::new(&model, dims, precalculate(), 10.0);
        assert!(!non_cache.within(&model, 0.0001));
    }

    #[test]
    fn eval_deriv_adds_out_of_bounds_force() {
        let atom = carbon_atom(Vec3::new(-1.0, 0.5, 0.5));
        let mut model = FakeModel {
            atoms: vec![atom.clone()],
            coords: vec![atom.coords],
            grid_atoms: vec![carbon_atom(Vec3::new(0.5, 0.5, 0.5))],
            minus_forces: vec![ZERO_VEC],
            atom_in_ligand: vec![false],
        };
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let non_cache = NonCache::new(&model, dims, precalculate(), 10.0);
        let e = non_cache.eval_deriv(&mut model, 1000.0);
        assert!(eq_fl(e, 10.0));
        assert!(eq_vec(model.minus_forces[0], Vec3::new(-10.0, 0.0, 0.0)));
    }

    #[test]
    fn eval_intra_skips_ligand_atoms() {
        let atom = carbon_atom(Vec3::new(-1.0, 0.5, 0.5));
        let mut model = FakeModel {
            atoms: vec![atom.clone()],
            coords: vec![atom.coords],
            grid_atoms: vec![carbon_atom(Vec3::new(0.5, 0.5, 0.5))],
            minus_forces: vec![ZERO_VEC],
            atom_in_ligand: vec![true],
        };
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let non_cache = NonCache::new(&model, dims, precalculate(), 10.0);
        assert!(eq_fl(non_cache.eval_intra(&mut model, 1000.0), 0.0));
    }
}
