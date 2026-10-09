use crate::atom::Atom;
use crate::common::{sqr, Fl, Pair, EPSILON_FL};
use crate::matrix::TriangularMatrix;
use crate::scoring_function::ScoringFunction;

#[derive(Debug, Clone, PartialEq)]
pub struct PrecalculateElement {
    pub smooth: Vec<Pair>,
    fast: Vec<Fl>,
    factor: Fl,
}

impl PrecalculateElement {
    pub fn new(n: usize, factor: Fl) -> Self {
        Self {
            smooth: vec![(0.0, 0.0); n],
            fast: vec![0.0; n],
            factor,
        }
    }

    pub fn eval_fast(&self, r2: Fl) -> Fl {
        assert!(r2 * self.factor < self.fast.len() as Fl);
        let i = (self.factor * r2) as usize;
        self.fast[i]
    }

    pub fn eval_deriv(&self, r2: Fl) -> Pair {
        let r2_factored = self.factor * r2;
        assert!(r2_factored + 1.0 < self.smooth.len() as Fl);
        let i1 = r2_factored as usize;
        let i2 = i1 + 1;
        let rem = r2_factored - i1 as Fl;
        assert!(rem >= -EPSILON_FL);
        assert!(rem < 1.0 + EPSILON_FL);
        let p1 = self.smooth[i1];
        let p2 = self.smooth[i2];
        let e = p1.0 + rem * (p2.0 - p1.0);
        let dor = p1.1 + rem * (p2.1 - p1.1);
        (e, dor)
    }

    pub fn init_from_smooth_fst(&mut self, rs: &[Fl]) {
        let n = self.smooth.len();
        assert_eq!(rs.len(), n);
        assert_eq!(self.fast.len(), n);
        for i in 0..n {
            if i == 0 || i == n - 1 {
                self.smooth[i].1 = 0.0;
            } else {
                let delta = rs[i + 1] - rs[i - 1];
                let r = rs[i];
                self.smooth[i].1 = (self.smooth[i + 1].0 - self.smooth[i - 1].0) / (delta * r);
            }
            let f1 = self.smooth[i].0;
            let f2 = if i + 1 >= n {
                0.0
            } else {
                self.smooth[i + 1].0
            };
            self.fast[i] = (f2 + f1) / 2.0;
        }
    }

    pub fn min_smooth_fst(&self) -> usize {
        let mut tmp = 0;
        for i_inv in 0..self.smooth.len() {
            let i = self.smooth.len() - i_inv - 1;
            if i_inv == 0 || self.smooth[i].0 < self.smooth[tmp].0 {
                tmp = i;
            }
        }
        tmp
    }

    pub fn widen_smooth_fst(&mut self, rs: &[Fl], left: Fl, right: Fl) {
        let mut tmp = vec![0.0; self.smooth.len()];
        let min_index = self.min_smooth_fst();
        assert!(min_index < rs.len());
        assert_eq!(rs.len(), self.smooth.len());
        let optimal_r = rs[min_index];
        for i in 0..self.smooth.len() {
            let mut r = rs[i];
            if r < optimal_r - left {
                r += left;
            } else if r > optimal_r + right {
                r -= right;
            } else {
                r = optimal_r;
            }
            if r < 0.0 {
                r = 0.0;
            }
            if r > *rs.last().expect("rs non-empty") {
                r = *rs.last().expect("rs non-empty");
            }
            tmp[i] = self.eval_deriv(sqr(r)).0;
        }
        for (pair, value) in self.smooth.iter_mut().zip(tmp) {
            pair.0 = value;
        }
    }

    pub fn widen(&mut self, rs: &[Fl], left: Fl, right: Fl) {
        self.widen_smooth_fst(rs, left, right);
        self.init_from_smooth_fst(rs);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Precalculate {
    cutoff_sqr: Fl,
    max_cutoff_sqr: Fl,
    n: usize,
    factor: Fl,
    data: TriangularMatrix<PrecalculateElement>,
}

impl Precalculate {
    pub fn new(sf: &ScoringFunction, v: Fl, factor: Fl) -> Self {
        let cutoff_sqr = sqr(sf.cutoff());
        let max_cutoff_sqr = sqr(sf.max_cutoff());
        let n = (factor * max_cutoff_sqr) as usize + 3;
        assert!(factor > EPSILON_FL);
        assert!((max_cutoff_sqr * factor) as usize + 1 < n);
        assert!(max_cutoff_sqr * factor + 1.0 < n as Fl);
        let mut data =
            TriangularMatrix::new(sf.num_atom_types(), PrecalculateElement::new(n, factor));
        let rs = calculate_rs(n, factor);

        for t1 in 0..data.dim() {
            for t2 in t1..data.dim() {
                let p = data.get_mut(t1, t2);
                for (i, smooth) in p.smooth.iter_mut().enumerate() {
                    smooth.0 = v.min(sf.eval_types(t1, t2, rs[i]));
                }
                p.init_from_smooth_fst(&rs);
            }
        }

        Self {
            cutoff_sqr,
            max_cutoff_sqr,
            n,
            factor,
            data,
        }
    }

    pub fn eval_fast(&self, type_pair_index: usize, r2: Fl) -> Fl {
        assert!(r2 <= self.max_cutoff_sqr);
        self.data.get_flat(type_pair_index).eval_fast(r2)
    }

    pub fn eval_deriv(&self, type_pair_index: usize, r2: Fl) -> Pair {
        assert!(r2 <= self.max_cutoff_sqr);
        self.data.get_flat(type_pair_index).eval_deriv(r2)
    }

    pub fn index_permissive(&self, t1: usize, t2: usize) -> usize {
        self.data.index_permissive(t1, t2)
    }

    pub fn cutoff_sqr(&self) -> Fl {
        self.cutoff_sqr
    }

    pub fn max_cutoff_sqr(&self) -> Fl {
        self.max_cutoff_sqr
    }

    pub fn widen(&mut self, left: Fl, right: Fl) {
        let rs = calculate_rs(self.n, self.factor);
        for t1 in 0..self.data.dim() {
            for t2 in t1..self.data.dim() {
                self.data.get_mut(t1, t2).widen(&rs, left, right);
            }
        }
    }
}

pub trait PrecalculateByAtomModel {
    fn num_atoms(&self) -> usize;
    fn atoms(&self) -> &[Atom];
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrecalculateByAtom {
    cutoff_sqr: Fl,
    max_cutoff_sqr: Fl,
    n: usize,
    factor: Fl,
    data: TriangularMatrix<PrecalculateElement>,
}

impl PrecalculateByAtom {
    pub fn new<M: PrecalculateByAtomModel>(
        sf: &ScoringFunction,
        model: &M,
        v: Fl,
        factor: Fl,
    ) -> Self {
        let cutoff_sqr = sqr(sf.cutoff());
        let max_cutoff_sqr = sqr(sf.max_cutoff());
        let n = (factor * max_cutoff_sqr) as usize + 3;
        let atoms = model.atoms();
        assert_eq!(atoms.len(), model.num_atoms());
        assert!(factor > EPSILON_FL);
        assert!((max_cutoff_sqr * factor) as usize + 1 < n);
        assert!(max_cutoff_sqr * factor + 1.0 < n as Fl);
        let mut data =
            TriangularMatrix::new(model.num_atoms(), PrecalculateElement::new(n, factor));
        let rs = calculate_rs(n, factor);

        for i in 0..data.dim() {
            for j in i..data.dim() {
                let p = data.get_mut(i, j);
                for (k, smooth) in p.smooth.iter_mut().enumerate() {
                    smooth.0 = v.min(sf.eval_atom_refs(&atoms[i], &atoms[j], rs[k]));
                }
                p.init_from_smooth_fst(&rs);
            }
        }

        Self {
            cutoff_sqr,
            max_cutoff_sqr,
            n,
            factor,
            data,
        }
    }

    pub fn eval_fast(&self, i: usize, j: usize, r2: Fl) -> Fl {
        assert!(r2 <= self.max_cutoff_sqr);
        self.data.get(i, j).eval_fast(r2)
    }

    pub fn eval_deriv(&self, i: usize, j: usize, r2: Fl) -> Pair {
        assert!(r2 <= self.max_cutoff_sqr);
        self.data.get(i, j).eval_deriv(r2)
    }

    pub fn index_permissive(&self, t1: usize, t2: usize) -> usize {
        self.data.index_permissive(t1, t2)
    }

    pub fn cutoff_sqr(&self) -> Fl {
        self.cutoff_sqr
    }

    pub fn max_cutoff_sqr(&self) -> Fl {
        self.max_cutoff_sqr
    }

    pub fn factor(&self) -> Fl {
        self.factor
    }

    pub fn widen(&mut self, left: Fl, right: Fl) {
        let rs = calculate_rs(self.n, self.factor);
        for t1 in 0..self.data.dim() {
            for t2 in t1..self.data.dim() {
                self.data.get_mut(t1, t2).widen(&rs, left, right);
            }
        }
    }
}

fn calculate_rs(n: usize, factor: Fl) -> Vec<Fl> {
    (0..n).map(|i| (i as Fl / factor).sqrt()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::atom_constants::{AD_TYPE_C, EL_TYPE_C, XS_TYPE_C_H};
    use crate::atom_type::AtomType;
    use crate::common::MAX_FL;
    use crate::scoring_function::ScoringFunctionChoice;

    #[test]
    fn element_interpolates_derivative_table() {
        let mut element = PrecalculateElement::new(4, 1.0);
        element.smooth[0].0 = 0.0;
        element.smooth[1].0 = 1.0;
        element.smooth[2].0 = 4.0;
        element.smooth[3].0 = 9.0;
        element.init_from_smooth_fst(&[0.0, 1.0, 2.0, 3.0]);
        assert_eq!(element.eval_fast(1.0), 2.5);
        assert_eq!(element.eval_deriv(1.5).0, 2.5);
    }

    #[test]
    fn builds_type_precalculate_table() {
        let sf = ScoringFunction::new(
            ScoringFunctionChoice::Vina,
            vec![
                -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50.0, 0.05846,
            ],
        );
        let p = Precalculate::new(&sf, MAX_FL, 4.0);
        let index = p.index_permissive(XS_TYPE_C_H, XS_TYPE_C_H);
        let _ = p.eval_fast(index, 1.0);
    }

    struct FakeModel {
        atoms: Vec<Atom>,
    }

    impl PrecalculateByAtomModel for FakeModel {
        fn num_atoms(&self) -> usize {
            self.atoms.len()
        }

        fn atoms(&self) -> &[Atom] {
            &self.atoms
        }
    }

    #[test]
    fn builds_atom_precalculate_table() {
        let sf = ScoringFunction::new(
            ScoringFunctionChoice::Vina,
            vec![
                -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50.0, 0.05846,
            ],
        );
        let mut atom = Atom::default();
        atom.base.atom_type = AtomType {
            el: EL_TYPE_C,
            ad: AD_TYPE_C,
            xs: XS_TYPE_C_H,
            sy: 0,
        };
        let p = PrecalculateByAtom::new(&sf, &FakeModel { atoms: vec![atom] }, MAX_FL, 4.0);
        let _ = p.eval_fast(0, 0, 1.0);
    }
}
