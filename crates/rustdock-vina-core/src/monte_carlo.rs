use crate::common::{Fl, Vec3, MAX_FL};
use crate::conf::{Change, Conf, ConfSize, OutputContainer, OutputType};
use crate::coords::add_to_output_container;
use crate::incrementable::Incrementable;
use crate::mutate::{mutate_conf, GyrationRadiusModel};
use crate::quasi_newton::{QuasiNewton, QuasiNewtonModel};
use crate::random::{random_fl, Rng64};

pub trait MonteCarloModel: QuasiNewtonModel + GyrationRadiusModel {
    fn size(&self) -> ConfSize;
    fn heavy_atom_movable_coords(&self) -> Vec<Vec3>;
}

pub fn metropolis_accept(old_f: Fl, new_f: Fl, temperature: Fl, generator: &mut Rng64) -> bool {
    if new_f < old_f {
        return true;
    }
    let acceptance_probability = ((old_f - new_f) / temperature).exp();
    random_fl(0.0, 1.0, generator) < acceptance_probability
}

#[derive(Debug, Clone, PartialEq)]
pub struct MonteCarlo {
    pub max_evals: u32,
    pub global_steps: u32,
    pub temperature: Fl,
    pub hunt_cap: Vec3,
    pub min_rmsd: Fl,
    pub num_saved_mins: usize,
    pub mutation_amplitude: Fl,
    pub local_steps: u32,
}

impl Default for MonteCarlo {
    fn default() -> Self {
        Self {
            max_evals: 0,
            global_steps: 2500,
            temperature: 1.2,
            hunt_cap: Vec3::new(10.0, 1.5, 10.0),
            min_rmsd: 0.5,
            num_saved_mins: 50,
            mutation_amplitude: 2.0,
            local_steps: 0,
        }
    }
}

impl MonteCarlo {
    pub fn run_best<M: MonteCarloModel>(
        &self,
        model: &mut M,
        corner1: Vec3,
        corner2: Vec3,
        increment_me: Option<&mut dyn Incrementable>,
        generator: &mut Rng64,
    ) -> OutputType {
        let mut out = OutputContainer::new();
        self.run(model, &mut out, corner1, corner2, increment_me, generator);
        assert!(!out.is_empty());
        out.remove(0)
    }

    pub fn run<M: MonteCarloModel>(
        &self,
        model: &mut M,
        out: &mut OutputContainer,
        corner1: Vec3,
        corner2: Vec3,
        mut increment_me: Option<&mut dyn Incrementable>,
        generator: &mut Rng64,
    ) {
        let mut evalcount = 0;
        let authentic_v = Vec3::new(1000.0, 1000.0, 1000.0);
        let size = model.size();
        let mut gradient = Change::new(&size);
        let mut tmp = OutputType::new(Conf::new(&size), 0.0);
        tmp.c.randomize(corner1, corner2, generator);
        let mut best_e = MAX_FL;
        let quasi_newton = QuasiNewton {
            max_steps: self.local_steps,
            average_required_improvement: 0.0,
        };

        for step in 0..self.global_steps {
            if let Some(increment_me) = increment_me.as_deref_mut() {
                increment_me.increment();
            }
            if self.max_evals > 0 && evalcount > self.max_evals as i32 {
                break;
            }
            let mut candidate = tmp.clone();
            mutate_conf(&mut candidate.c, model, self.mutation_amplitude, generator);
            quasi_newton.optimize(
                model,
                &mut candidate,
                &mut gradient,
                self.hunt_cap,
                &mut evalcount,
            );
            if step == 0 || metropolis_accept(tmp.e, candidate.e, self.temperature, generator) {
                tmp = candidate;
                if tmp.e < best_e || out.len() < self.num_saved_mins {
                    quasi_newton.optimize(
                        model,
                        &mut tmp,
                        &mut gradient,
                        authentic_v,
                        &mut evalcount,
                    );
                    tmp.coords = model.heavy_atom_movable_coords();
                    add_to_output_container(out, tmp.clone(), self.min_rmsd, self.num_saved_mins);
                    if tmp.e < best_e {
                        best_e = tmp.e;
                    }
                }
            }
        }
        assert!(!out.is_empty());
        assert!(out.first().unwrap().e <= out.last().unwrap().e);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfgs::ChangeVector;

    #[derive(Clone)]
    struct FakeModel {
        size: ConfSize,
    }

    impl GyrationRadiusModel for FakeModel {
        fn gyration_radius(&self, _ligand_index: usize) -> Fl {
            1.0
        }
    }

    impl QuasiNewtonModel for FakeModel {
        fn set_conf(&mut self, _conf: &Conf) {}

        fn eval_deriv(&mut self, _v: Vec3, gradient: &mut Change) -> Fl {
            for i in 0..gradient.num_floats() {
                gradient.set(i, 0.0);
            }
            0.0
        }
    }

    impl MonteCarloModel for FakeModel {
        fn size(&self) -> ConfSize {
            self.size.clone()
        }

        fn heavy_atom_movable_coords(&self) -> Vec<Vec3> {
            vec![Vec3::new(0.0, 0.0, 0.0)]
        }
    }

    #[test]
    fn accepts_lower_energy_and_runs_smoke() {
        let mut rng = Rng64::new(5);
        assert!(metropolis_accept(2.0, 1.0, 1.2, &mut rng));

        let mc = MonteCarlo {
            global_steps: 2,
            num_saved_mins: 2,
            ..Default::default()
        };
        let mut model = FakeModel {
            size: ConfSize {
                ligands: vec![1],
                flex: vec![],
            },
        };
        let mut out = Vec::new();
        mc.run(
            &mut model,
            &mut out,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0),
            None,
            &mut rng,
        );
        assert!(!out.is_empty());
    }
}
