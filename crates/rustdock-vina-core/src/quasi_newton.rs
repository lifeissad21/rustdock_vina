use crate::bfgs::bfgs;
use crate::common::{Fl, Vec3};
use crate::conf::{Change, Conf, OutputType};

pub trait QuasiNewtonModel {
    fn set_conf(&mut self, conf: &Conf);
    fn eval_deriv(&mut self, v: Vec3, gradient: &mut Change) -> Fl;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuasiNewton {
    pub max_steps: u32,
    pub average_required_improvement: Fl,
}

impl Default for QuasiNewton {
    fn default() -> Self {
        Self {
            max_steps: 1000,
            average_required_improvement: 0.0,
        }
    }
}

impl QuasiNewton {
    pub fn optimize<M: QuasiNewtonModel>(
        &self,
        model: &mut M,
        out: &mut OutputType,
        gradient: &mut Change,
        v: Vec3,
        evalcount: &mut i32,
    ) {
        let mut aux = |conf: &Conf, gradient: &mut Change| {
            model.set_conf(conf);
            model.eval_deriv(v, gradient)
        };
        let result = bfgs(
            &mut aux,
            &mut out.c,
            gradient,
            self.max_steps,
            self.average_required_improvement,
            10,
            evalcount,
        );
        model.set_conf(&out.c);
        out.e = result;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfgs::ChangeVector;
    use crate::conf::{ConfSize, OutputType};

    struct FakeModel;

    impl QuasiNewtonModel for FakeModel {
        fn set_conf(&mut self, _conf: &Conf) {}

        fn eval_deriv(&mut self, _v: Vec3, gradient: &mut Change) -> Fl {
            gradient.set(0, 0.0);
            1.5
        }
    }

    #[test]
    fn updates_output_energy() {
        let size = ConfSize {
            ligands: vec![],
            flex: vec![1],
        };
        let mut out = OutputType::new(Conf::new(&size), 0.0);
        let mut gradient = Change::new(&size);
        let mut evalcount = 0;
        QuasiNewton {
            max_steps: 0,
            average_required_improvement: 0.0,
        }
        .optimize(
            &mut FakeModel,
            &mut out,
            &mut gradient,
            Vec3::new(1.0, 1.0, 1.0),
            &mut evalcount,
        );
        assert_eq!(out.e, 1.5);
        assert_eq!(evalcount, 1);
    }
}
