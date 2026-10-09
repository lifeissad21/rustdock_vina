use crate::common::Fl;

pub trait IGrid<M> {
    fn eval(&self, model: &M, v: Fl) -> Fl;
    fn eval_intra(&self, model: &mut M, v: Fl) -> Fl;
    fn eval_deriv(&self, model: &mut M, v: Fl) -> Fl;
}
