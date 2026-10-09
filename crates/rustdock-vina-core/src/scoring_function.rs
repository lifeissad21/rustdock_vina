use crate::atom::Atom;
use crate::atom_type::{num_atom_types, AtomTyping};
use crate::common::Fl;
use crate::conf_independent::{ad4_tors_add, num_tors_div, ConfIndependentInputs};
use crate::potentials::{
    Ad4Electrostatic, Ad4Hb, Ad4Solvation, Ad4Vdw, Gaussian, Hydrophobic, LinearAttraction,
    NonDirHBond, Potential, Repulsion,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringFunctionChoice {
    Vina,
    Ad42,
    Vinardo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfIndependentTerm {
    NumTorsDiv,
    Ad4TorsAdd,
}

pub struct ScoringFunction {
    potentials: Vec<Box<dyn Potential>>,
    conf_independents: Vec<ConfIndependentTerm>,
    weights: Vec<Fl>,
    cutoff: Fl,
    max_cutoff: Fl,
    atom_typing: AtomTyping,
}

impl ScoringFunction {
    pub fn new(choice: ScoringFunctionChoice, weights: Vec<Fl>) -> Self {
        let mut potentials: Vec<Box<dyn Potential>> = Vec::new();
        let mut conf_independents = Vec::new();
        let cutoff;
        let max_cutoff;
        let atom_typing;

        match choice {
            ScoringFunctionChoice::Vina => {
                potentials.push(Box::new(Gaussian {
                    offset: 0.0,
                    width: 0.5,
                    cutoff: 8.0,
                    vinardo: false,
                }));
                potentials.push(Box::new(Gaussian {
                    offset: 3.0,
                    width: 2.0,
                    cutoff: 8.0,
                    vinardo: false,
                }));
                potentials.push(Box::new(Repulsion {
                    offset: 0.0,
                    cutoff: 8.0,
                    vinardo: false,
                }));
                potentials.push(Box::new(Hydrophobic {
                    good: 0.5,
                    bad: 1.5,
                    cutoff: 8.0,
                    vinardo: false,
                }));
                potentials.push(Box::new(NonDirHBond {
                    good: -0.7,
                    bad: 0.0,
                    cutoff: 8.0,
                    vinardo: false,
                }));
                potentials.push(Box::new(LinearAttraction { cutoff: 20.0 }));
                conf_independents.push(ConfIndependentTerm::NumTorsDiv);
                atom_typing = AtomTyping::Xs;
                cutoff = 8.0;
                max_cutoff = 20.0;
            }
            ScoringFunctionChoice::Vinardo => {
                potentials.push(Box::new(Gaussian {
                    offset: 0.0,
                    width: 0.8,
                    cutoff: 8.0,
                    vinardo: true,
                }));
                potentials.push(Box::new(Repulsion {
                    offset: 0.0,
                    cutoff: 8.0,
                    vinardo: true,
                }));
                potentials.push(Box::new(Hydrophobic {
                    good: 0.0,
                    bad: 2.5,
                    cutoff: 8.0,
                    vinardo: true,
                }));
                potentials.push(Box::new(NonDirHBond {
                    good: -0.6,
                    bad: 0.0,
                    cutoff: 8.0,
                    vinardo: true,
                }));
                potentials.push(Box::new(LinearAttraction { cutoff: 20.0 }));
                conf_independents.push(ConfIndependentTerm::NumTorsDiv);
                atom_typing = AtomTyping::Xs;
                cutoff = 8.0;
                max_cutoff = 20.0;
            }
            ScoringFunctionChoice::Ad42 => {
                potentials.push(Box::new(Ad4Vdw {
                    smoothing: 0.5,
                    cap: 100000.0,
                    cutoff: 8.0,
                }));
                potentials.push(Box::new(Ad4Hb {
                    smoothing: 0.5,
                    cap: 100000.0,
                    cutoff: 8.0,
                }));
                potentials.push(Box::new(Ad4Electrostatic {
                    cap: 100.0,
                    cutoff: 20.48,
                }));
                potentials.push(Box::new(Ad4Solvation {
                    desolvation_sigma: 3.6,
                    solvation_q: 0.01097,
                    charge_dependent: true,
                    cutoff: 20.48,
                }));
                potentials.push(Box::new(LinearAttraction { cutoff: 20.0 }));
                conf_independents.push(ConfIndependentTerm::Ad4TorsAdd);
                atom_typing = AtomTyping::Ad;
                cutoff = 20.48;
                max_cutoff = 20.48;
            }
        }

        Self {
            potentials,
            conf_independents,
            weights,
            cutoff,
            max_cutoff,
            atom_typing,
        }
    }

    pub fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        self.potentials
            .iter()
            .enumerate()
            .map(|(i, potential)| {
                self.weights.get(i).copied().unwrap_or(0.0) * potential.eval_atoms(a, b, r)
            })
            .sum()
    }

    pub fn eval_atom_refs(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        self.eval_atoms(a, b, r)
    }

    pub fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        self.potentials
            .iter()
            .enumerate()
            .map(|(i, potential)| {
                self.weights.get(i).copied().unwrap_or(0.0) * potential.eval_types(t1, t2, r)
            })
            .sum()
    }

    pub fn conf_independent_from_inputs(&self, input: &ConfIndependentInputs, mut e: Fl) -> Fl {
        let mut index = self.potentials.len();
        for term in &self.conf_independents {
            e = match term {
                ConfIndependentTerm::NumTorsDiv => {
                    num_tors_div(input, e, &self.weights, &mut index)
                }
                ConfIndependentTerm::Ad4TorsAdd => {
                    ad4_tors_add(input, e, &self.weights, &mut index)
                }
            };
        }
        assert_eq!(index, self.weights.len());
        e
    }

    pub fn cutoff(&self) -> Fl {
        self.cutoff
    }

    pub fn max_cutoff(&self) -> Fl {
        self.max_cutoff
    }

    pub fn atom_typing(&self) -> AtomTyping {
        self.atom_typing
    }

    pub fn get_atom_typing(&self) -> AtomTyping {
        self.atom_typing()
    }

    pub fn atom_types(&self) -> Vec<usize> {
        (0..num_atom_types(self.atom_typing)).collect()
    }

    pub fn num_atom_types(&self) -> usize {
        num_atom_types(self.atom_typing)
    }

    pub fn weights(&self) -> &[Fl] {
        &self.weights
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_constants::XS_TYPE_C_H;

    #[test]
    fn vina_scoring_function_has_reference_shape() {
        let sf = ScoringFunction::new(
            ScoringFunctionChoice::Vina,
            vec![
                -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50.0, 0.05846,
            ],
        );
        assert_eq!(sf.atom_typing(), AtomTyping::Xs);
        assert_eq!(sf.cutoff(), 8.0);
        assert_eq!(sf.max_cutoff(), 20.0);
        assert_eq!(sf.atom_types().len(), sf.num_atom_types());
        let _ = sf.eval_types(XS_TYPE_C_H, XS_TYPE_C_H, 4.0);
    }
}
