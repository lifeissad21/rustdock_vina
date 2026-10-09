use crate::ad4cache::Ad4Cache;
use crate::ad4cache::Ad4CacheModel;
use crate::atom::{Atom, AtomIndex, Bond};
use crate::atom_constants::*;
use crate::atom_type::{num_atom_types, AtomTyping};
use crate::cache::Cache;
use crate::cache::CacheModel;
use crate::common::{sqr, vec_distance_sqr, Fl, Vec3, MAX_FL, ZERO_VEC};
use crate::conf::{Change, Conf, ConfSize};
use crate::coords::rmsd_upper_bound;
use crate::curl::{curl, curl_with_derivative};
use crate::igrid::IGrid;
use crate::matrix::StrictlyTriangularMatrix;
use crate::monte_carlo::MonteCarloModel;
use crate::mutate::GyrationRadiusModel;
use crate::non_cache::NonCache;
use crate::non_cache::NonCacheModel;
use crate::precalculate::{PrecalculateByAtom, PrecalculateByAtomModel};
use crate::quasi_newton::QuasiNewtonModel;
use crate::szv_grid::SzvGridModel;
use crate::tree::{count_torsions, FlexibleBody, MainBranch};
use crate::triangular_matrix_index::triangular_matrix_index_permissive;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InteractingPair {
    pub type_pair_index: usize,
    pub a: usize,
    pub b: usize,
}

impl InteractingPair {
    pub const fn new(type_pair_index: usize, a: usize, b: usize) -> Self {
        Self {
            type_pair_index,
            a,
            b,
        }
    }
}

pub type InteractingPairs = Vec<InteractingPair>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedLine {
    pub text: String,
    pub atom_index: Option<usize>,
}

impl ParsedLine {
    pub fn new(text: impl Into<String>, atom_index: Option<usize>) -> Self {
        Self {
            text: text.into(),
            atom_index,
        }
    }
}

pub type Context = Vec<ParsedLine>;

#[derive(Debug, Clone, PartialEq)]
pub struct Ligand {
    pub body: FlexibleBody,
    pub range: crate::tree::AtomRange,
    pub degrees_of_freedom: usize,
    pub pairs: InteractingPairs,
    pub context: Context,
}

impl Ligand {
    pub fn new(body: FlexibleBody, degrees_of_freedom: usize) -> Self {
        Self {
            body,
            range: crate::tree::AtomRange::new(0, 0),
            degrees_of_freedom,
            pairs: Vec::new(),
            context: Vec::new(),
        }
    }

    pub fn set_range(&mut self) {
        self.range = flexible_body_atom_range(&self.body);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub branch: MainBranch,
}

impl Residue {
    pub fn new(branch: MainBranch) -> Self {
        Self { branch }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceType {
    Fixed,
    Rotor,
    Variable,
}

pub type DistanceTypeMatrix = StrictlyTriangularMatrix<DistanceType>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelPairClass {
    LigandInternal,
    Other,
    Inter,
    Glue,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub coords: Vec<Vec3>,
    pub minus_forces: Vec<Vec3>,
    pub grid_atoms: Vec<Atom>,
    pub atoms: Vec<Atom>,
    pub ligands: Vec<Ligand>,
    pub flex: Vec<Residue>,
    pub flex_context: Context,
    pub other_pairs: InteractingPairs,
    pub inter_pairs: InteractingPairs,
    pub glue_pairs: InteractingPairs,
    num_movable_atoms: usize,
    atom_typing_used: AtomTyping,
}

impl Default for Model {
    fn default() -> Self {
        Self::new(AtomTyping::Xs)
    }
}

impl Model {
    pub fn new(atom_typing_used: AtomTyping) -> Self {
        Self {
            coords: Vec::new(),
            minus_forces: Vec::new(),
            grid_atoms: Vec::new(),
            atoms: Vec::new(),
            ligands: Vec::new(),
            flex: Vec::new(),
            flex_context: Vec::new(),
            other_pairs: Vec::new(),
            inter_pairs: Vec::new(),
            glue_pairs: Vec::new(),
            num_movable_atoms: 0,
            atom_typing_used,
        }
    }

    pub fn set_num_movable_atoms(&mut self, value: usize) {
        assert!(value <= self.atoms.len());
        self.num_movable_atoms = value;
    }

    pub fn push_atom(&mut self, atom: Atom, movable: bool) {
        if movable {
            assert_eq!(self.num_movable_atoms, self.atoms.len());
            self.num_movable_atoms += 1;
        }
        self.coords.push(atom.coords);
        self.minus_forces.push(ZERO_VEC);
        self.atoms.push(atom);
    }

    pub fn push_grid_atom(&mut self, atom: Atom) {
        self.grid_atoms.push(atom);
    }

    pub fn atom_typing_used(&self) -> AtomTyping {
        self.atom_typing_used
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn num_movable_atoms(&self) -> usize {
        self.num_movable_atoms
    }

    pub fn num_ligands(&self) -> usize {
        self.ligands.len()
    }

    pub fn num_flex(&self) -> usize {
        self.flex.len()
    }

    pub fn num_internal_pairs(&self) -> usize {
        self.ligands.iter().map(|ligand| ligand.pairs.len()).sum()
    }

    pub fn num_other_pairs(&self) -> usize {
        self.other_pairs.len()
    }

    pub fn ligand_degrees_of_freedom(&self, ligand_number: usize) -> usize {
        self.ligands[ligand_number].degrees_of_freedom
    }

    pub fn find_ligand(&self, atom_index: usize) -> Option<usize> {
        self.ligands
            .iter()
            .position(|ligand| atom_index >= ligand.range.begin && atom_index < ligand.range.end)
    }

    pub fn is_atom_in_ligand(&self, atom_index: usize) -> bool {
        self.find_ligand(atom_index).is_some()
    }

    pub fn is_movable_atom(&self, atom_index: usize) -> bool {
        atom_index < self.num_movable_atoms
    }

    pub fn center(&self) -> Vec3 {
        if self.num_movable_atoms == 0 {
            return ZERO_VEC;
        }
        let mut center = ZERO_VEC;
        for coord in &self.coords[..self.num_movable_atoms] {
            center += *coord;
        }
        (1.0 / self.num_movable_atoms as Fl) * center
    }

    pub fn get_movable_atom_types(&self, atom_typing_used: AtomTyping) -> Vec<usize> {
        let n = num_atom_types(atom_typing_used);
        let mut out = Vec::new();
        for atom in &self.atoms[..self.num_movable_atoms] {
            let t = atom.base.atom_type.get(atom_typing_used);
            if t < n && !out.contains(&t) {
                out.push(t);
            }
        }
        out
    }

    pub fn get_ligand_coords_vec3(&self) -> Vec<Vec3> {
        assert_eq!(self.ligands.len(), 1);
        let ligand = &self.ligands[0];
        self.coords[ligand.range.begin..ligand.range.end].to_vec()
    }

    pub fn get_ligand_coords_flat(&self) -> Vec<Fl> {
        let mut out = Vec::new();
        for coord in self.get_ligand_coords_vec3() {
            out.extend_from_slice(&coord.data);
        }
        out
    }

    pub fn get_heavy_atom_movable_coords(&self) -> Vec<Vec3> {
        self.atoms[..self.num_movable_atoms]
            .iter()
            .zip(&self.coords)
            .filter_map(|(atom, coord)| {
                if atom.base.atom_type.el != EL_TYPE_H {
                    Some(*coord)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_size(&self) -> ConfSize {
        ConfSize {
            ligands: self
                .ligands
                .iter()
                .map(|ligand| count_torsions(&ligand.body))
                .collect(),
            flex: self
                .flex
                .iter()
                .map(|residue| count_torsions(&residue.branch))
                .collect(),
        }
    }

    pub fn get_initial_conf(&self) -> Conf {
        let mut conf = Conf::new(&self.get_size());
        conf.set_to_null();
        for (out, ligand) in conf.ligands.iter_mut().zip(&self.ligands) {
            out.rigid.position = ligand.body.node.atom_frame.frame.origin();
        }
        conf
    }

    pub fn set(&mut self, c: &Conf) {
        assert_eq!(self.ligands.len(), c.ligands.len());
        assert_eq!(self.flex.len(), c.flex.len());
        for (ligand, conf) in self.ligands.iter_mut().zip(&c.ligands) {
            ligand
                .body
                .set_ligand_conf(&self.atoms, &mut self.coords, conf);
        }
        for (residue, conf) in self.flex.iter_mut().zip(&c.flex) {
            residue
                .branch
                .set_residue_conf(&self.atoms, &mut self.coords, conf);
        }
    }

    pub fn gyration_radius(&self, ligand_number: usize) -> Fl {
        let ligand = &self.ligands[ligand_number];
        let origin = ligand.body.node.atom_frame.frame.origin();
        let mut acc = 0.0;
        let mut counter = 0;
        for i in ligand.range.begin..ligand.range.end {
            if self.atoms[i].base.atom_type.el != EL_TYPE_H {
                acc += vec_distance_sqr(self.coords[i], origin);
                counter += 1;
            }
        }
        if counter > 0 {
            (acc / counter as Fl).sqrt()
        } else {
            0.0
        }
    }

    pub fn evalo(&self, p: &PrecalculateByAtom, v: Vec3) -> Fl {
        eval_interacting_pairs(p, v[2], &self.other_pairs, &self.coords, false)
    }

    pub fn eval_inter(&self, p: &PrecalculateByAtom, v: Vec3) -> Fl {
        eval_interacting_pairs(p, v[2], &self.inter_pairs, &self.coords, false)
    }

    pub fn evali(&self, p: &PrecalculateByAtom, v: Vec3) -> Fl {
        self.ligands
            .iter()
            .map(|ligand| eval_interacting_pairs(p, v[0], &ligand.pairs, &self.coords, false))
            .sum()
    }

    pub fn eval_deriv_with_grid<G>(
        &mut self,
        p: &PrecalculateByAtom,
        ig: &G,
        v: Vec3,
        g: &mut Change,
    ) -> Fl
    where
        G: IGrid<Model>,
    {
        let mut e = ig.eval_deriv(self, v[1]);
        for ligand in &self.ligands {
            e += eval_interacting_pairs_deriv(
                p,
                v[0],
                &ligand.pairs,
                &self.coords,
                &mut self.minus_forces,
                false,
            );
        }
        e += eval_interacting_pairs_deriv(
            p,
            v[2],
            &self.inter_pairs,
            &self.coords,
            &mut self.minus_forces,
            false,
        );
        e += eval_interacting_pairs_deriv(
            p,
            v[2],
            &self.other_pairs,
            &self.coords,
            &mut self.minus_forces,
            false,
        );
        e += eval_interacting_pairs_deriv(
            p,
            v[2],
            &self.glue_pairs,
            &self.coords,
            &mut self.minus_forces,
            true,
        );
        self.write_derivatives(g);
        e
    }

    pub fn eval_intramolecular<G>(&mut self, p: &PrecalculateByAtom, ig: &G, v: Vec3) -> Fl
    where
        G: IGrid<Model>,
    {
        let mut e = self.evali(p, v);
        e += ig.eval_intra(self, v[1]);
        e += eval_interacting_pairs(p, v[2], &self.other_pairs, &self.coords, false);
        e
    }

    pub fn rmsd_lower_bound(&self, other: &Self) -> Fl {
        self.rmsd_lower_bound_asymmetric(other)
            .max(other.rmsd_lower_bound_asymmetric(self))
    }

    pub fn rmsd_upper_bound_model(&self, other: &Self) -> Fl {
        assert_eq!(self.num_movable_atoms, other.num_movable_atoms);
        let a = self.get_heavy_atom_movable_coords();
        let b = other.get_heavy_atom_movable_coords();
        rmsd_upper_bound(&a, &b)
    }

    pub fn rmsd_ligands_upper_bound(&self, other: &Self) -> Fl {
        assert_eq!(self.ligands.len(), other.ligands.len());
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (ligand, other_ligand) in self.ligands.iter().zip(&other.ligands) {
            assert_eq!(ligand.range, other_ligand.range);
            for i in ligand.range.begin..ligand.range.end {
                assert_eq!(
                    self.atoms[i].base.atom_type.ad,
                    other.atoms[i].base.atom_type.ad
                );
                assert_eq!(
                    self.atoms[i].base.atom_type.xs,
                    other.atoms[i].base.atom_type.xs
                );
                if self.atoms[i].base.atom_type.el != EL_TYPE_H {
                    a.push(self.coords[i]);
                    b.push(other.coords[i]);
                }
            }
        }
        rmsd_upper_bound(&a, &b)
    }

    pub fn clash_penalty(&self) -> Fl {
        self.clash_penalty_aux(&self.other_pairs)
            + self.clash_penalty_aux(&self.inter_pairs)
            + self
                .ligands
                .iter()
                .map(|ligand| self.clash_penalty_aux(&ligand.pairs))
                .sum::<Fl>()
    }

    pub fn atom_coords(&self, index: AtomIndex) -> Vec3 {
        if index.in_grid {
            self.grid_atoms[index.i].coords
        } else {
            self.coords[index.i]
        }
    }

    pub fn distance_sqr_between(&self, a: AtomIndex, b: AtomIndex) -> Fl {
        vec_distance_sqr(self.atom_coords(a), self.atom_coords(b))
    }

    pub fn distance_type_between(
        &self,
        mobility: &DistanceTypeMatrix,
        i: AtomIndex,
        j: AtomIndex,
    ) -> DistanceType {
        if i.in_grid && j.in_grid {
            return DistanceType::Fixed;
        }
        if i.in_grid {
            return if j.i < self.num_movable_atoms {
                DistanceType::Variable
            } else {
                DistanceType::Fixed
            };
        }
        if j.in_grid {
            return if i.i < self.num_movable_atoms {
                DistanceType::Variable
            } else {
                DistanceType::Fixed
            };
        }
        if i.i == j.i {
            DistanceType::Fixed
        } else {
            *mobility.get(mobility_i(i.i, j.i).0, mobility_i(i.i, j.i).1)
        }
    }

    pub fn sz_to_atom_index(&self, i: usize) -> AtomIndex {
        if i < self.grid_atoms.len() {
            AtomIndex::new(i, true)
        } else {
            AtomIndex::new(i - self.grid_atoms.len(), false)
        }
    }

    pub fn initialize(&mut self, mobility: &DistanceTypeMatrix) {
        for ligand in &mut self.ligands {
            ligand.set_range();
        }
        self.assign_bonds(mobility);
        self.assign_types();
        self.initialize_pairs(mobility);
    }

    pub fn assign_bonds(&mut self, mobility: &DistanceTypeMatrix) {
        for atom in self.grid_atoms.iter_mut().chain(&mut self.atoms) {
            atom.bonds.clear();
        }
        let n = self.grid_atoms.len() + self.atoms.len();
        // Spatial buckets only narrow candidates; bond criteria match the reference.
        let mut buckets = std::collections::HashMap::<(i32, i32, i32), Vec<usize>>::new();
        let key = |v: Vec3| {
            (
                (v[0] / 4.0).floor() as i32,
                (v[1] / 4.0).floor() as i32,
                (v[2] / 4.0).floor() as i32,
            )
        };
        for i in 0..n {
            buckets
                .entry(key(self.atom_coords(self.sz_to_atom_index(i))))
                .or_default()
                .push(i);
        }
        for i in 0..n {
            let a = self.sz_to_atom_index(i);
            let atom = self.get_atom_ref(a);
            let radius = if atom.base.atom_type.ad < AD_TYPE_SIZE {
                ad_type_property(atom.base.atom_type.ad).covalent_radius
            } else {
                max_covalent_radius()
            };
            let cutoff = sqr(1.1 * (radius + max_covalent_radius()));
            let (x, y, z) = key(self.atom_coords(a));
            let mut relevant = Vec::new();
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(indices) = buckets.get(&(x + dx, y + dy, z + dz)) {
                            relevant.extend(indices.iter().copied().filter(|&j| {
                                j != i
                                    && self.distance_type_between(
                                        mobility,
                                        a,
                                        self.sz_to_atom_index(j),
                                    ) != DistanceType::Variable
                                    && self.distance_sqr_between(a, self.sz_to_atom_index(j))
                                        < cutoff
                            }));
                        }
                    }
                }
            }
            for &j in &relevant {
                if j <= i {
                    continue;
                }
                let b = self.sz_to_atom_index(j);
                let length = self
                    .get_atom_ref(a)
                    .base
                    .atom_type
                    .optimal_covalent_bond_length(self.get_atom_ref(b).base.atom_type);
                let r2 = self.distance_sqr_between(a, b);
                let between = relevant.iter().any(|&k| {
                    let c = self.sz_to_atom_index(k);
                    c != a
                        && c != b
                        && self.distance_type_between(mobility, a, c) != DistanceType::Variable
                        && self.distance_type_between(mobility, b, c) != DistanceType::Variable
                        && self.distance_sqr_between(a, c) < r2
                        && self.distance_sqr_between(b, c) < r2
                });
                if r2 < sqr(1.1 * length) && !between {
                    let rotor = self.distance_type_between(mobility, a, b) == DistanceType::Rotor;
                    self.get_atom_mut(a)
                        .bonds
                        .push(Bond::new(b, r2.sqrt(), rotor));
                    self.get_atom_mut(b)
                        .bonds
                        .push(Bond::new(a, r2.sqrt(), rotor));
                }
            }
        }
    }

    pub fn append(&mut self, other: &Model) {
        assert_eq!(self.atom_typing_used, other.atom_typing_used);
        let a_m = self.num_movable_atoms;
        let b_m = other.num_movable_atoms;
        let a_n = self.atoms.len();
        let a_g = self.grid_atoms.len();
        let map_a = |i: usize| if i < a_m { i } else { i + b_m };
        let map_b = |i: usize| if i < b_m { i + a_m } else { i + a_n };
        fn remap(model: &mut Model, f: impl Fn(usize) -> usize + Copy, grid_offset: usize) {
            for atom in model.atoms.iter_mut().chain(&mut model.grid_atoms) {
                for bond in &mut atom.bonds {
                    bond.connected_atom_index.i = if bond.connected_atom_index.in_grid {
                        bond.connected_atom_index.i + grid_offset
                    } else {
                        f(bond.connected_atom_index.i)
                    };
                }
            }
            let pairs = |pairs: &mut InteractingPairs| {
                for p in pairs {
                    p.a = f(p.a);
                    p.b = f(p.b);
                }
            };
            let context = |context: &mut Context| {
                for line in context {
                    line.atom_index = line.atom_index.map(f);
                }
            };
            for ligand in &mut model.ligands {
                crate::tree::transform_ranges(&mut ligand.body, f);
                ligand.range.transform(f);
                pairs(&mut ligand.pairs);
                context(&mut ligand.context);
            }
            for residue in &mut model.flex {
                crate::tree::transform_ranges(&mut residue.branch, f);
            }
            context(&mut model.flex_context);
            pairs(&mut model.other_pairs);
            pairs(&mut model.inter_pairs);
            pairs(&mut model.glue_pairs);
        }
        remap(self, map_a, 0);
        let mut b = other.clone();
        remap(&mut b, map_b, a_g);
        self.atoms.splice(a_m..a_m, b.atoms.drain(..b_m));
        self.atoms.extend(b.atoms);
        self.coords.splice(a_m..a_m, b.coords.drain(..b_m));
        self.coords.extend(b.coords);
        self.minus_forces
            .splice(a_m..a_m, b.minus_forces.drain(..b_m));
        self.minus_forces.extend(b.minus_forces);
        self.grid_atoms.extend(b.grid_atoms);
        self.ligands.extend(b.ligands);
        self.flex.extend(b.flex);
        self.flex_context.extend(b.flex_context);
        self.other_pairs.extend(b.other_pairs);
        self.inter_pairs.extend(b.inter_pairs);
        self.glue_pairs.extend(b.glue_pairs);
        for i in 0..a_m {
            for j in a_m..a_m + b_m {
                if self.is_closure_clash(i, j) || self.is_unmatched_closure_dummy(i, j) {
                    continue;
                }
                let t1 = self.atoms[i].base.atom_type.get(self.atom_typing_used);
                let t2 = self.atoms[j].base.atom_type.get(self.atom_typing_used);
                let n = num_atom_types(self.atom_typing_used);
                if t1 >= n || t2 >= n {
                    continue;
                }
                let p = InteractingPair::new(triangular_matrix_index_permissive(n, t1, t2), i, j);
                if self.is_glue_pair(i, j) {
                    self.glue_pairs.push(p);
                } else if self.is_atom_in_ligand(i) || self.is_atom_in_ligand(j) {
                    self.inter_pairs.push(p);
                } else {
                    self.other_pairs.push(p);
                }
            }
        }
        self.num_movable_atoms += b_m;
    }

    pub fn assign_types(&mut self) {
        for i in 0..(self.grid_atoms.len() + self.atoms.len()) {
            let ai = self.sz_to_atom_index(i);
            let bonded_to_hd = self.bonded_to_hd(ai);
            let bonded_to_heteroatom = self.bonded_to_heteroatom(ai);
            let atom = self.get_atom_mut(ai);
            atom.base.atom_type.assign_el();
            let acceptor =
                atom.base.atom_type.ad == AD_TYPE_OA || atom.base.atom_type.ad == AD_TYPE_NA;
            let donor_n_or_o = atom.base.atom_type.el == EL_TYPE_MET || bonded_to_hd;
            atom.base.atom_type.xs = match atom.base.atom_type.el {
                EL_TYPE_H => atom.base.atom_type.xs,
                EL_TYPE_C => match atom.base.atom_type.ad {
                    AD_TYPE_CG0 => {
                        if bonded_to_heteroatom {
                            XS_TYPE_C_P_CG0
                        } else {
                            XS_TYPE_C_H_CG0
                        }
                    }
                    AD_TYPE_CG1 => {
                        if bonded_to_heteroatom {
                            XS_TYPE_C_P_CG1
                        } else {
                            XS_TYPE_C_H_CG1
                        }
                    }
                    AD_TYPE_CG2 => {
                        if bonded_to_heteroatom {
                            XS_TYPE_C_P_CG2
                        } else {
                            XS_TYPE_C_H_CG2
                        }
                    }
                    AD_TYPE_CG3 => {
                        if bonded_to_heteroatom {
                            XS_TYPE_C_P_CG3
                        } else {
                            XS_TYPE_C_H_CG3
                        }
                    }
                    _ => {
                        if bonded_to_heteroatom {
                            XS_TYPE_C_P
                        } else {
                            XS_TYPE_C_H
                        }
                    }
                },
                EL_TYPE_N => {
                    if acceptor && donor_n_or_o {
                        XS_TYPE_N_DA
                    } else if acceptor {
                        XS_TYPE_N_A
                    } else if donor_n_or_o {
                        XS_TYPE_N_D
                    } else {
                        XS_TYPE_N_P
                    }
                }
                EL_TYPE_O => {
                    if acceptor && donor_n_or_o {
                        XS_TYPE_O_DA
                    } else if acceptor {
                        XS_TYPE_O_A
                    } else if donor_n_or_o {
                        XS_TYPE_O_D
                    } else {
                        XS_TYPE_O_P
                    }
                }
                EL_TYPE_S => XS_TYPE_S_P,
                EL_TYPE_P => XS_TYPE_P_P,
                EL_TYPE_F => XS_TYPE_F_H,
                EL_TYPE_CL => XS_TYPE_CL_H,
                EL_TYPE_BR => XS_TYPE_BR_H,
                EL_TYPE_I => XS_TYPE_I_H,
                EL_TYPE_SI => XS_TYPE_SI,
                EL_TYPE_AT => XS_TYPE_AT,
                EL_TYPE_MET => XS_TYPE_MET_D,
                EL_TYPE_DUMMY => match atom.base.atom_type.ad {
                    AD_TYPE_G0 => XS_TYPE_G0,
                    AD_TYPE_G1 => XS_TYPE_G1,
                    AD_TYPE_G2 => XS_TYPE_G2,
                    AD_TYPE_G3 => XS_TYPE_G3,
                    AD_TYPE_W => XS_TYPE_SIZE,
                    _ => panic!("unexpected dummy atom type"),
                },
                EL_TYPE_SIZE => atom.base.atom_type.xs,
                _ => panic!("unexpected element type"),
            };
        }
    }

    pub fn initialize_pairs(&mut self, mobility: &DistanceTypeMatrix) {
        self.other_pairs.clear();
        self.inter_pairs.clear();
        self.glue_pairs.clear();
        for ligand in &mut self.ligands {
            ligand.pairs.clear();
        }

        for i in 0..self.atoms.len() {
            let i_lig = self.find_ligand(i);
            let bonded_atoms = self.bonded_to(i, 3);
            for j in (i + 1)..self.atoms.len() {
                if *mobility.get(i, j) != DistanceType::Variable || bonded_atoms.contains(&j) {
                    continue;
                }
                if self.is_closure_clash(i, j) || self.is_unmatched_closure_dummy(i, j) {
                    continue;
                }
                let t1 = self.atoms[i].base.atom_type.get(self.atom_typing_used);
                let t2 = self.atoms[j].base.atom_type.get(self.atom_typing_used);
                let n = num_atom_types(self.atom_typing_used);
                if t1 >= n || t2 >= n {
                    continue;
                }
                let ip = InteractingPair::new(triangular_matrix_index_permissive(n, t1, t2), i, j);
                if self.is_glue_pair(i, j) {
                    self.glue_pairs.push(ip);
                } else if let Some(ligand_index) = i_lig {
                    if self.find_ligand(j) == Some(ligand_index) {
                        self.ligands[ligand_index].pairs.push(ip);
                    }
                } else if !self.is_atom_in_ligand(j) {
                    self.other_pairs.push(ip);
                }
            }
        }
    }

    pub fn write_model(&self, model_number: usize, remark: &str) -> String {
        let mut out = format!("MODEL {model_number}\n{remark}");
        for ligand in &self.ligands {
            out.push_str(&self.write_context(&ligand.context));
        }
        if !self.flex.is_empty() {
            out.push_str(&self.write_context(&self.flex_context));
        }
        out.push_str("ENDMDL\n");
        out
    }

    fn write_context(&self, context: &Context) -> String {
        let mut out = String::new();
        for line in context {
            if let Some(index) = line.atom_index {
                out.push_str(&coords_to_pdbqt_string(self.coords[index], &line.text));
            } else {
                out.push_str(&line.text);
            }
            out.push('\n');
        }
        out
    }

    fn write_derivatives(&self, g: &mut Change) {
        for (ligand, change) in self.ligands.iter().zip(&mut g.ligands) {
            ligand
                .body
                .derivative_ligand(&self.coords, &self.minus_forces, change);
        }
        for (residue, change) in self.flex.iter().zip(&mut g.flex) {
            residue
                .branch
                .derivative_residue(&self.coords, &self.minus_forces, change);
        }
    }

    fn get_atom_ref(&self, index: AtomIndex) -> &Atom {
        if index.in_grid {
            &self.grid_atoms[index.i]
        } else {
            &self.atoms[index.i]
        }
    }

    fn get_atom_mut(&mut self, index: AtomIndex) -> &mut Atom {
        if index.in_grid {
            &mut self.grid_atoms[index.i]
        } else {
            &mut self.atoms[index.i]
        }
    }

    fn bonded_to_hd(&self, index: AtomIndex) -> bool {
        self.get_atom_ref(index).bonds.iter().any(|bond| {
            self.get_atom_ref(bond.connected_atom_index)
                .base
                .atom_type
                .ad
                == AD_TYPE_HD
        })
    }

    fn bonded_to_heteroatom(&self, index: AtomIndex) -> bool {
        self.get_atom_ref(index).bonds.iter().any(|bond| {
            self.get_atom_ref(bond.connected_atom_index)
                .base
                .atom_type
                .is_heteroatom()
        })
    }

    fn bonded_to(&self, atom_index: usize, depth: usize) -> Vec<usize> {
        let mut out = Vec::new();
        self.bonded_to_rec(atom_index, depth, &mut out);
        out
    }

    fn bonded_to_rec(&self, atom_index: usize, depth: usize, out: &mut Vec<usize>) {
        if out.contains(&atom_index) {
            return;
        }
        out.push(atom_index);
        if depth == 0 {
            return;
        }
        for bond in &self.atoms[atom_index].bonds {
            if !bond.connected_atom_index.in_grid {
                self.bonded_to_rec(bond.connected_atom_index.i, depth - 1, out);
            }
        }
    }

    fn is_glue_pair(&self, i: usize, j: usize) -> bool {
        matches_glue_pair(
            self.atoms[i].base.atom_type.ad,
            self.atoms[j].base.atom_type.ad,
        )
    }

    fn is_unmatched_closure_dummy(&self, i: usize, j: usize) -> bool {
        unmatched_closure_dummy(
            self.atoms[i].base.atom_type.ad,
            self.atoms[j].base.atom_type.ad,
        )
    }

    fn is_closure_clash(&self, i: usize, j: usize) -> bool {
        let t1 = self.atoms[i].base.atom_type.ad;
        let t2 = self.atoms[j].base.atom_type.ad;
        if matches_glue_pair(t1, t2) {
            return false;
        }
        let neighbors_i = self.bonded_to(i, 1);
        let neighbors_j = self.bonded_to(j, 1);
        let mut i_has = [false; 4];
        for index in neighbors_i {
            match self.atoms[index].base.atom_type.ad {
                AD_TYPE_CG0 => i_has[0] = true,
                AD_TYPE_CG1 => i_has[1] = true,
                AD_TYPE_CG2 => i_has[2] = true,
                AD_TYPE_CG3 => i_has[3] = true,
                _ => {}
            }
        }
        neighbors_j.into_iter().any(|index| {
            matches!(
                (self.atoms[index].base.atom_type.ad, i_has),
                (AD_TYPE_CG0, [true, _, _, _])
                    | (AD_TYPE_CG1, [_, true, _, _])
                    | (AD_TYPE_CG2, [_, _, true, _])
                    | (AD_TYPE_CG3, [_, _, _, true])
            )
        })
    }

    fn clash_penalty_aux(&self, pairs: &[InteractingPair]) -> Fl {
        let mut e = 0.0;
        for pair in pairs {
            let a = &self.atoms[pair.a];
            let b = &self.atoms[pair.b];
            let covalent_r = a
                .base
                .atom_type
                .optimal_covalent_bond_length(b.base.atom_type);
            e += pairwise_clash_penalty(
                vec_distance_sqr(self.coords[pair.a], self.coords[pair.b]).sqrt(),
                covalent_r,
            );
        }
        e
    }

    fn rmsd_lower_bound_asymmetric(&self, other: &Self) -> Fl {
        assert_eq!(self.num_movable_atoms, other.num_movable_atoms);
        let mut sum = 0.0;
        let mut counter = 0;
        for i in 0..self.num_movable_atoms {
            let a = &self.atoms[i];
            if a.base.atom_type.el == EL_TYPE_H {
                continue;
            }
            let mut best = MAX_FL;
            for j in 0..self.num_movable_atoms {
                let b = &other.atoms[j];
                if a.base.atom_type.same_element(b.base.atom_type)
                    && !b.base.atom_type.is_hydrogen()
                {
                    best = best.min(vec_distance_sqr(self.coords[i], other.coords[j]));
                }
            }
            assert!(best < 0.1 * MAX_FL);
            sum += best;
            counter += 1;
        }
        if counter == 0 {
            0.0
        } else {
            (sum / counter as Fl).sqrt()
        }
    }
}

impl SzvGridModel for Model {
    fn atom_typing_used(&self) -> AtomTyping {
        self.atom_typing_used
    }

    fn grid_atoms(&self) -> &[Atom] {
        &self.grid_atoms
    }
}

impl PrecalculateByAtomModel for Model {
    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}

impl CacheModel for Model {
    fn num_movable_atoms(&self) -> usize {
        self.num_movable_atoms
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
        Model::is_atom_in_ligand(self, atom_index)
    }
}

impl NonCacheModel for Model {
    fn num_movable_atoms(&self) -> usize {
        self.num_movable_atoms
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
        Model::is_atom_in_ligand(self, atom_index)
    }
}

impl Ad4CacheModel for Model {
    fn num_movable_atoms(&self) -> usize {
        self.num_movable_atoms
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
        Model::is_atom_in_ligand(self, atom_index)
    }
}

impl IGrid<Model> for Cache {
    fn eval(&self, model: &Model, v: Fl) -> Fl {
        Cache::eval(self, model, v)
    }

    fn eval_intra(&self, model: &mut Model, v: Fl) -> Fl {
        Cache::eval_intra(self, model, v)
    }

    fn eval_deriv(&self, model: &mut Model, v: Fl) -> Fl {
        Cache::eval_deriv(self, model, v)
    }
}

impl IGrid<Model> for NonCache {
    fn eval(&self, model: &Model, v: Fl) -> Fl {
        NonCache::eval(self, model, v)
    }

    fn eval_intra(&self, model: &mut Model, v: Fl) -> Fl {
        NonCache::eval_intra(self, model, v)
    }

    fn eval_deriv(&self, model: &mut Model, v: Fl) -> Fl {
        NonCache::eval_deriv(self, model, v)
    }
}

impl IGrid<Model> for Ad4Cache {
    fn eval(&self, model: &Model, v: Fl) -> Fl {
        Ad4Cache::eval(self, model, v)
    }

    fn eval_intra(&self, model: &mut Model, v: Fl) -> Fl {
        Ad4Cache::eval_intra(self, model, v)
    }

    fn eval_deriv(&self, model: &mut Model, v: Fl) -> Fl {
        Ad4Cache::eval_deriv(self, model, v)
    }
}

impl GyrationRadiusModel for Model {
    fn gyration_radius(&self, ligand_index: usize) -> Fl {
        Model::gyration_radius(self, ligand_index)
    }
}

impl crate::conf_independent::ConfIndependentModel for Model {
    fn num_ligands(&self) -> usize {
        self.ligands.len()
    }
    fn ligand(&self, index: usize) -> crate::conf_independent::LigandInfo {
        let l = &self.ligands[index];
        crate::conf_independent::LigandInfo {
            begin: l.range.begin,
            end: l.range.end,
            degrees_of_freedom: l.degrees_of_freedom as Fl,
        }
    }
    fn ligand_length(&self, index: usize) -> usize {
        fn depth(branch: &crate::tree::Branch) -> usize {
            1 + branch.children.iter().map(depth).max().unwrap_or(0)
        }
        self.ligands[index]
            .body
            .children
            .iter()
            .map(depth)
            .max()
            .unwrap_or(0)
    }
    fn atom(&self, index: AtomIndex) -> &Atom {
        self.get_atom_ref(index)
    }
}

pub struct SearchModel<'a, G> {
    pub model: Model,
    pub precalculate: &'a PrecalculateByAtom,
    pub grid: &'a G,
}
impl<G> Clone for SearchModel<'_, G> {
    fn clone(&self) -> Self {
        Self {
            model: self.model.clone(),
            precalculate: self.precalculate,
            grid: self.grid,
        }
    }
}
impl<G: IGrid<Model>> QuasiNewtonModel for SearchModel<'_, G> {
    fn set_conf(&mut self, conf: &Conf) {
        self.model.set(conf);
    }
    fn eval_deriv(&mut self, v: Vec3, gradient: &mut Change) -> Fl {
        self.model
            .eval_deriv_with_grid(self.precalculate, self.grid, v, gradient)
    }
}
impl<G> GyrationRadiusModel for SearchModel<'_, G> {
    fn gyration_radius(&self, index: usize) -> Fl {
        self.model.gyration_radius(index)
    }
}
impl<G: IGrid<Model>> MonteCarloModel for SearchModel<'_, G> {
    fn size(&self) -> ConfSize {
        self.model.get_size()
    }
    fn heavy_atom_movable_coords(&self) -> Vec<Vec3> {
        self.model.get_heavy_atom_movable_coords()
    }
}

pub fn eval_interacting_pairs(
    p: &PrecalculateByAtom,
    v: Fl,
    pairs: &[InteractingPair],
    coords: &[Vec3],
    with_max_cutoff: bool,
) -> Fl {
    let cutoff_sqr = if with_max_cutoff {
        p.max_cutoff_sqr()
    } else {
        p.cutoff_sqr()
    };
    let mut e = 0.0;
    for pair in pairs {
        let r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
        if r2 < cutoff_sqr {
            let mut tmp = p.eval_fast(pair.a, pair.b, r2);
            curl(&mut tmp, v);
            e += tmp;
        }
    }
    e
}

pub fn eval_interacting_pairs_deriv(
    p: &PrecalculateByAtom,
    v: Fl,
    pairs: &[InteractingPair],
    coords: &[Vec3],
    forces: &mut [Vec3],
    with_max_cutoff: bool,
) -> Fl {
    let cutoff_sqr = if with_max_cutoff {
        p.max_cutoff_sqr()
    } else {
        p.cutoff_sqr()
    };
    let mut e = 0.0;
    for pair in pairs {
        let r = coords[pair.b] - coords[pair.a];
        let r2 = r.norm_sqr();
        if r2 < cutoff_sqr {
            let (mut pair_e, dor) = p.eval_deriv(pair.a, pair.b, r2);
            let mut force = dor * r;
            curl_with_derivative(&mut pair_e, &mut force, v);
            e += pair_e;
            forces[pair.a] -= force;
            forces[pair.b] += force;
        }
    }
    e
}

pub fn pairwise_clash_penalty(r: Fl, covalent_r: Fl) -> Fl {
    assert!(r >= 0.0);
    assert!(covalent_r > crate::common::EPSILON_FL);
    let x = r / covalent_r;
    if x > 2.0 {
        0.0
    } else {
        1.0 - sqr(x) / 4.0
    }
}

pub fn coords_to_pdbqt_string(coords: Vec3, line: &str) -> String {
    let mut out = line.to_string();
    write_pdbqt_coord(&mut out, 31, coords[0]);
    write_pdbqt_coord(&mut out, 39, coords[1]);
    write_pdbqt_coord(&mut out, 47, coords[2]);
    out
}

fn write_pdbqt_coord(line: &mut String, one_based: usize, value: Fl) {
    assert!(one_based > 0);
    let start = one_based - 1;
    assert!(line.len() >= start + 8);
    let formatted = format!("{value:8.3}");
    line.replace_range(start..start + 8, &formatted);
}

fn flexible_body_atom_range(body: &FlexibleBody) -> crate::tree::AtomRange {
    let mut range = body.node.atom_frame.range;
    for child in &body.children {
        let child_range = branch_atom_range(child);
        range.begin = range.begin.min(child_range.begin);
        range.end = range.end.max(child_range.end);
    }
    range
}

fn branch_atom_range(branch: &crate::tree::Branch) -> crate::tree::AtomRange {
    let mut range = branch.node.axis_frame.atom_frame.range;
    for child in &branch.children {
        let child_range = branch_atom_range(child);
        range.begin = range.begin.min(child_range.begin);
        range.end = range.end.max(child_range.end);
    }
    range
}

fn matches_glue_pair(t1: usize, t2: usize) -> bool {
    matches!(
        (t1, t2),
        (AD_TYPE_CG0, AD_TYPE_G0)
            | (AD_TYPE_G0, AD_TYPE_CG0)
            | (AD_TYPE_CG1, AD_TYPE_G1)
            | (AD_TYPE_G1, AD_TYPE_CG1)
            | (AD_TYPE_CG2, AD_TYPE_G2)
            | (AD_TYPE_G2, AD_TYPE_CG2)
            | (AD_TYPE_CG3, AD_TYPE_G3)
            | (AD_TYPE_G3, AD_TYPE_CG3)
    )
}

fn unmatched_closure_dummy(t1: usize, t2: usize) -> bool {
    matches!(
        (t1, t2),
        (AD_TYPE_G0, t) if t != AD_TYPE_CG0
    ) || matches!((t1, t2), (t, AD_TYPE_G0) if t != AD_TYPE_CG0)
        || matches!((t1, t2), (AD_TYPE_G1, t) if t != AD_TYPE_CG1)
        || matches!((t1, t2), (t, AD_TYPE_G1) if t != AD_TYPE_CG1)
        || matches!((t1, t2), (AD_TYPE_G2, t) if t != AD_TYPE_CG2)
        || matches!((t1, t2), (t, AD_TYPE_G2) if t != AD_TYPE_CG2)
        || matches!((t1, t2), (AD_TYPE_G3, t) if t != AD_TYPE_CG3)
        || matches!((t1, t2), (t, AD_TYPE_G3) if t != AD_TYPE_CG3)
}

fn mobility_i(i: usize, j: usize) -> (usize, usize) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_type::AtomType;
    use crate::common::{eq_fl, eq_vec};
    use crate::conf::RigidConf;
    use crate::quaternion::QT_IDENTITY;
    use crate::tree::{Heterotree, RigidBody};

    fn typed_atom(el: usize, ad: usize, xs: usize, coords: Vec3) -> Atom {
        let mut atom = Atom {
            coords,
            ..Atom::default()
        };
        atom.base.atom_type = AtomType { el, ad, xs, sy: 0 };
        atom
    }

    #[test]
    fn model_basic_queries_match_reference_shape() {
        let mut model = Model::default();
        model.push_atom(
            typed_atom(EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, Vec3::new(0.0, 0.0, 0.0)),
            true,
        );
        model.push_atom(
            typed_atom(EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_A, Vec3::new(2.0, 0.0, 0.0)),
            true,
        );
        model.ligands.push(Ligand {
            body: Heterotree::new(RigidBody::new(ZERO_VEC, 0, 2)),
            range: crate::tree::AtomRange::new(0, 2),
            degrees_of_freedom: 6,
            pairs: Vec::new(),
            context: Vec::new(),
        });
        assert_eq!(model.num_movable_atoms(), 2);
        assert!(model.is_atom_in_ligand(1));
        assert_eq!(model.find_ligand(1), Some(0));
        assert!(eq_vec(model.center(), Vec3::new(1.0, 0.0, 0.0)));
        assert_eq!(
            model.get_movable_atom_types(AtomTyping::Xs),
            vec![XS_TYPE_C_H, XS_TYPE_O_A]
        );
    }

    #[test]
    fn set_updates_ligand_coordinates_from_conf() {
        let mut model = Model::default();
        model.push_atom(
            typed_atom(EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, Vec3::new(1.0, 0.0, 0.0)),
            true,
        );
        model.ligands.push(Ligand {
            body: Heterotree::new(RigidBody::new(ZERO_VEC, 0, 1)),
            range: crate::tree::AtomRange::new(0, 1),
            degrees_of_freedom: 6,
            pairs: Vec::new(),
            context: Vec::new(),
        });
        let mut conf = model.get_initial_conf();
        conf.ligands[0].rigid = RigidConf {
            position: Vec3::new(5.0, 0.0, 0.0),
            orientation: QT_IDENTITY,
        };
        model.set(&conf);
        assert!(eq_vec(model.coords[0], Vec3::new(6.0, 0.0, 0.0)));
    }

    #[test]
    fn interaction_pair_evaluation_uses_atom_precalculate_indices() {
        let pair = InteractingPair::new(0, 0, 1);
        let coords = [ZERO_VEC, Vec3::new(2.0, 0.0, 0.0)];
        let mut forces = vec![ZERO_VEC; 2];
        let mut fake = Model::default();
        fake.push_atom(
            typed_atom(EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, coords[0]),
            true,
        );
        fake.push_atom(
            typed_atom(EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, coords[1]),
            true,
        );
        let sf = crate::scoring_function::ScoringFunction::new(
            crate::scoring_function::ScoringFunctionChoice::Vina,
            vec![0.0; 6],
        );
        let p = PrecalculateByAtom::new(&sf, &fake, 1000.0, 16.0);
        assert!(eq_fl(
            eval_interacting_pairs(&p, 1000.0, &[pair], &coords, false),
            0.0
        ));
        assert!(eq_fl(
            eval_interacting_pairs_deriv(&p, 1000.0, &[pair], &coords, &mut forces, false),
            0.0
        ));
    }

    #[test]
    fn pdbqt_coord_rewrite_matches_fixed_columns() {
        let line = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00    +0.000 C";
        let out = coords_to_pdbqt_string(Vec3::new(1.25, -2.5, 3.75), line);
        assert_eq!(&out[30..54], "   1.250  -2.500   3.750");
    }
}
