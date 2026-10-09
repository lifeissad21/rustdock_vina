use crate::common::{Fl, Vec3};
use crate::conf::{OutputContainer, OutputType};
use crate::coords::add_to_output_container;
use crate::monte_carlo::{MonteCarlo, MonteCarloModel};
use crate::parallel_progress::ParallelProgress;
use crate::random::{random_int, Rng64};

#[derive(Debug, Clone, PartialEq)]
pub struct ParallelMc {
    pub mc: MonteCarlo,
    pub num_tasks: usize,
    pub num_threads: usize,
    pub display_progress: bool,
}

impl Default for ParallelMc {
    fn default() -> Self {
        Self {
            mc: MonteCarlo::default(),
            num_tasks: 8,
            num_threads: 1,
            display_progress: true,
        }
    }
}

pub fn merge_output_containers(
    incoming: &[OutputType],
    out: &mut OutputContainer,
    min_rmsd: Fl,
    max_size: usize,
) {
    for output in incoming {
        add_to_output_container(out, output.clone(), min_rmsd, max_size);
    }
}

pub fn merge_many_output_containers(
    many: &[OutputContainer],
    out: &mut OutputContainer,
    _min_rmsd: Fl,
    max_size: usize,
) {
    let min_rmsd = 2.0;
    for incoming in many {
        merge_output_containers(incoming, out, min_rmsd, max_size);
    }
    out.sort_by(crate::conf::output_compare);
}

impl ParallelMc {
    pub fn run<M, F>(
        &self,
        model: &M,
        out: &mut OutputContainer,
        corner1: Vec3,
        corner2: Vec3,
        generator: &mut Rng64,
        progress_callback: Option<F>,
    ) where
        M: MonteCarloModel + Clone + Send + Sync,
        F: FnMut(f64),
    {
        let mut progress = ParallelProgress::new(
            progress_callback.map(|callback| Box::new(callback) as Box<dyn FnMut(f64)>),
        );
        if self.display_progress {
            progress.init((self.num_tasks as u32 * self.mc.global_steps) as u64);
        }

        assert!(self.num_tasks > 0 && self.num_threads > 0);
        let seeds: Vec<_> = (0..self.num_tasks)
            .map(|_| random_int(0, 1_000_000, generator) as u64)
            .collect();
        let mut task_outputs = vec![Vec::new(); self.num_tasks];
        let workers = self.num_threads.min(self.num_tasks);
        let chunk_size = self.num_tasks.div_ceil(workers);
        let (sender, receiver) = std::sync::mpsc::channel();
        std::thread::scope(|scope| {
            for (chunk, outputs) in task_outputs.chunks_mut(chunk_size).enumerate() {
                let sender = sender.clone();
                let seeds = &seeds;
                scope.spawn(move || {
                    for (offset, out) in outputs.iter_mut().enumerate() {
                        let seed = seeds[chunk * chunk_size + offset];
                        let mut task_model = model.clone();
                        let mut rng = Rng64::new(seed);
                        let mut steps = 0u64;
                        let mut tick = || {
                            steps += 1;
                            if steps % 100 == 0 {
                                let _ = sender.send(100);
                            }
                        };
                        self.mc.run(
                            &mut task_model,
                            out,
                            corner1,
                            corner2,
                            if self.display_progress {
                                Some(&mut tick)
                            } else {
                                None
                            },
                            &mut rng,
                        );
                        if self.display_progress && steps % 100 != 0 {
                            let _ = sender.send(steps % 100);
                        }
                    }
                });
            }
            drop(sender);
            for steps in receiver {
                progress.advance(steps);
            }
        });

        merge_many_output_containers(&task_outputs, out, self.mc.min_rmsd, self.mc.num_saved_mins);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfgs::ChangeVector;
    use crate::conf::{Change, Conf, ConfSize};
    use crate::mutate::GyrationRadiusModel;
    use crate::quasi_newton::QuasiNewtonModel;

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
    fn parallel_mc_runs_tasks_with_stable_merge_order() {
        let pmc = ParallelMc {
            mc: MonteCarlo {
                global_steps: 1,
                num_saved_mins: 2,
                ..Default::default()
            },
            num_tasks: 2,
            num_threads: 1,
            display_progress: false,
        };
        let model = FakeModel {
            size: ConfSize {
                ligands: vec![1],
                flex: vec![],
            },
        };
        let mut out = Vec::new();
        pmc.run(
            &model,
            &mut out,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0),
            &mut Rng64::new(3),
            Option::<fn(f64)>::None,
        );
        assert!(!out.is_empty());
    }
}
