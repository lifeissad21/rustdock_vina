use crate::atom::{Atom, AtomIndex};
use crate::atom_constants::{xs_is_acceptor, xs_is_donor, xs_is_hydrophobic, EL_TYPE_H};
use crate::common::{sqr, Fl, EPSILON_FL, MAX_FL};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LigandInfo {
    pub begin: usize,
    pub end: usize,
    pub degrees_of_freedom: Fl,
}

pub trait ConfIndependentModel {
    fn num_ligands(&self) -> usize;
    fn ligand(&self, index: usize) -> LigandInfo;
    fn ligand_length(&self, index: usize) -> usize;
    fn atom(&self, index: AtomIndex) -> &Atom;
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConfIndependentInputs {
    pub torsdof: Fl,
    pub num_tors: Fl,
    pub num_rotors: Fl,
    pub num_heavy_atoms: Fl,
    pub num_hydrophobic_atoms: Fl,
    pub ligand_max_num_h_bonds: Fl,
    pub num_ligands: Fl,
    pub ligand_lengths_sum: Fl,
}

impl Default for ConfIndependentInputs {
    fn default() -> Self {
        Self {
            torsdof: 0.0,
            num_tors: 0.0,
            num_rotors: 0.0,
            num_heavy_atoms: 0.0,
            num_hydrophobic_atoms: 0.0,
            ligand_max_num_h_bonds: 0.0,
            num_ligands: 0.0,
            ligand_lengths_sum: 0.0,
        }
    }
}

impl ConfIndependentInputs {
    pub fn from_model<M: ConfIndependentModel>(model: &M) -> Self {
        let mut inputs = Self {
            num_ligands: model.num_ligands() as Fl,
            ..Self::default()
        };

        for i in 0..model.num_ligands() {
            let ligand = model.ligand(i);
            inputs.ligand_lengths_sum += model.ligand_length(i) as Fl;
            inputs.torsdof += ligand.degrees_of_freedom;
            for j in ligand.begin..ligand.end {
                let atom = model.atom(AtomIndex::new(j, false));
                if atom.base.atom_type.el != EL_TYPE_H {
                    let ar = atom_rotors(model, AtomIndex::new(j, false));
                    inputs.num_tors += 0.5 * ar as Fl;
                    if ar > 2 {
                        inputs.num_rotors += 0.5;
                    } else {
                        inputs.num_rotors += 0.5 * ar as Fl;
                    }
                    if xs_is_hydrophobic(atom.base.atom_type.xs) {
                        inputs.num_hydrophobic_atoms += 1.0;
                    }
                    if xs_is_acceptor(atom.base.atom_type.xs) || xs_is_donor(atom.base.atom_type.xs)
                    {
                        inputs.ligand_max_num_h_bonds += 1.0;
                    }
                    inputs.num_heavy_atoms += 1.0;
                }
            }
        }
        inputs
    }

    pub fn as_vec(&self) -> Vec<Fl> {
        vec![
            self.num_tors,
            self.num_rotors,
            self.num_heavy_atoms,
            self.num_hydrophobic_atoms,
            self.ligand_max_num_h_bonds,
            self.num_ligands,
            self.ligand_lengths_sum,
        ]
    }

    pub fn names() -> Vec<&'static str> {
        vec![
            "num_tors",
            "num_rotors",
            "num_heavy_atoms",
            "num_hydrophobic_atoms",
            "ligand_max_num_h_bonds",
            "num_ligands",
            "ligand_lengths_sum",
        ]
    }
}

fn num_bonded_heavy_atoms<M: ConfIndependentModel>(model: &M, i: AtomIndex) -> u32 {
    model
        .atom(i)
        .bonds
        .iter()
        .filter(|bond| {
            !model
                .atom(bond.connected_atom_index)
                .base
                .atom_type
                .is_hydrogen()
        })
        .count() as u32
}

fn atom_rotors<M: ConfIndependentModel>(model: &M, i: AtomIndex) -> u32 {
    model
        .atom(i)
        .bonds
        .iter()
        .filter(|bond| {
            bond.rotatable
                && !model
                    .atom(bond.connected_atom_index)
                    .base
                    .atom_type
                    .is_hydrogen()
                && num_bonded_heavy_atoms(model, bond.connected_atom_index) > 1
        })
        .count() as u32
}

pub fn read_weight(weights: &[Fl], index: &mut usize) -> Fl {
    let value = weights[*index];
    *index += 1;
    value
}

pub fn conf_smooth_div(x: Fl, y: Fl) -> Fl {
    if x.abs() < EPSILON_FL {
        0.0
    } else if y.abs() < EPSILON_FL {
        if x * y > 0.0 {
            MAX_FL
        } else {
            -MAX_FL
        }
    } else {
        x / y
    }
}

pub fn num_tors_sqr(input: &ConfIndependentInputs, x: Fl, weights: &[Fl], index: &mut usize) -> Fl {
    let weight = 0.1 * read_weight(weights, index);
    x + weight * sqr(input.num_tors) / 5.0
}

pub fn num_tors_sqrt(
    input: &ConfIndependentInputs,
    x: Fl,
    weights: &[Fl],
    index: &mut usize,
) -> Fl {
    let weight = 0.1 * read_weight(weights, index);
    x + weight * input.num_tors.sqrt() / 5.0_f64.sqrt()
}

pub fn num_tors_div(input: &ConfIndependentInputs, x: Fl, weights: &[Fl], index: &mut usize) -> Fl {
    let weight = 0.1 * (read_weight(weights, index) + 1.0);
    conf_smooth_div(x, 1.0 + weight * input.num_tors / 5.0)
}

pub fn ligand_length(
    input: &ConfIndependentInputs,
    x: Fl,
    weights: &[Fl],
    index: &mut usize,
) -> Fl {
    x + read_weight(weights, index) * input.ligand_lengths_sum
}

pub fn num_ligands(input: &ConfIndependentInputs, x: Fl, weights: &[Fl], index: &mut usize) -> Fl {
    x + read_weight(weights, index) * input.num_ligands
}

pub fn num_heavy_atoms_div(
    input: &ConfIndependentInputs,
    x: Fl,
    weights: &[Fl],
    index: &mut usize,
) -> Fl {
    let weight = 0.05 * read_weight(weights, index);
    conf_smooth_div(x, 1.0 + weight * input.num_heavy_atoms)
}

pub fn num_heavy_atoms(
    input: &ConfIndependentInputs,
    x: Fl,
    weights: &[Fl],
    index: &mut usize,
) -> Fl {
    x + 0.05 * read_weight(weights, index) * input.num_heavy_atoms
}

pub fn num_hydrophobic_atoms(
    input: &ConfIndependentInputs,
    x: Fl,
    weights: &[Fl],
    index: &mut usize,
) -> Fl {
    x + 0.05 * read_weight(weights, index) * input.num_hydrophobic_atoms
}

pub fn ad4_tors_add(input: &ConfIndependentInputs, x: Fl, weights: &[Fl], index: &mut usize) -> Fl {
    x + read_weight(weights, index) * input.torsdof
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_fl;

    #[test]
    fn exposes_input_vector_and_names() {
        let input = ConfIndependentInputs {
            num_tors: 1.0,
            num_rotors: 2.0,
            num_heavy_atoms: 3.0,
            num_hydrophobic_atoms: 4.0,
            ligand_max_num_h_bonds: 5.0,
            num_ligands: 6.0,
            ligand_lengths_sum: 7.0,
            torsdof: 8.0,
        };
        assert_eq!(input.as_vec(), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert_eq!(ConfIndependentInputs::names().len(), input.as_vec().len());
    }

    #[test]
    fn evaluates_vina_torsion_div_formula() {
        let input = ConfIndependentInputs {
            num_tors: 5.0,
            ..Default::default()
        };
        let mut index = 0;
        assert!(eq_fl(
            num_tors_div(&input, 10.0, &[0.0], &mut index),
            10.0 / 1.1
        ));
    }
}
