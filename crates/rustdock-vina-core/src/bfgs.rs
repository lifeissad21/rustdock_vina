use crate::common::{Fl, EPSILON_FL};
use crate::matrix::TriangularMatrix;

pub trait ChangeVector: Clone {
    fn num_floats(&self) -> usize;
    fn get(&self, index: usize) -> Fl;
    fn set(&mut self, index: usize, value: Fl);
}

pub trait IncrementBy<P> {
    fn increment(&mut self, change: &P, factor: Fl);
}

impl ChangeVector for crate::conf::Change {
    fn num_floats(&self) -> usize {
        self.num_floats()
    }

    fn get(&self, index: usize) -> Fl {
        self.get(index)
    }

    fn set(&mut self, index: usize, value: Fl) {
        *self.get_mut(index) = value;
    }
}

impl IncrementBy<crate::conf::Change> for crate::conf::Conf {
    fn increment(&mut self, change: &crate::conf::Change, factor: Fl) {
        self.increment(change, factor);
    }
}

pub type FlMat = TriangularMatrix<Fl>;

pub fn minus_mat_vec_product<C: ChangeVector>(m: &FlMat, input: &C, output: &mut C) {
    let n = m.dim();
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += *m.get_permissive(i, j) * input.get(j);
        }
        output.set(i, -sum);
    }
}

pub fn scalar_product<C: ChangeVector>(a: &C, b: &C, n: usize) -> Fl {
    let mut tmp = 0.0;
    for i in 0..n {
        tmp += a.get(i) * b.get(i);
    }
    tmp
}

pub fn bfgs_update<C: ChangeVector>(h: &mut FlMat, p: &C, y: &C, alpha: Fl) -> bool {
    let yp = scalar_product(y, p, h.dim());
    if alpha * yp < EPSILON_FL {
        return false;
    }
    let mut minus_hy = y.clone();
    minus_mat_vec_product(h, y, &mut minus_hy);
    let yhy = -scalar_product(y, &minus_hy, h.dim());
    let r = 1.0 / (alpha * yp);
    let n = p.num_floats();
    for i in 0..n {
        for j in i..n {
            let value = *h.get(i, j)
                + alpha * r * (minus_hy.get(i) * p.get(j) + minus_hy.get(j) * p.get(i))
                + alpha * alpha * (r * r * yhy + r) * p.get(i) * p.get(j);
            *h.get_mut(i, j) = value;
        }
    }
    true
}

pub fn line_search<F, X, C>(
    f: &mut F,
    n: usize,
    x: &X,
    g: &C,
    f0: Fl,
    p: &C,
    x_new: &mut X,
    g_new: &mut C,
    f1: &mut Fl,
    evalcount: &mut i32,
) -> Fl
where
    F: FnMut(&X, &mut C) -> Fl,
    X: Clone + IncrementBy<C>,
    C: ChangeVector,
{
    let c0 = 0.0001;
    let max_trials = 10;
    let multiplier = 0.5;
    let mut alpha = 1.0;
    let pg = scalar_product(p, g, n);

    for _ in 0..max_trials {
        *x_new = x.clone();
        x_new.increment(p, alpha);
        *f1 = f(x_new, g_new);
        *evalcount += 1;
        if *f1 - f0 < c0 * alpha * pg {
            break;
        }
        alpha *= multiplier;
    }
    alpha
}

pub fn set_diagonal(m: &mut FlMat, x: Fl) {
    for i in 0..m.dim() {
        *m.get_mut(i, i) = x;
    }
}

pub fn subtract_change<C: ChangeVector>(b: &mut C, a: &C, n: usize) {
    for i in 0..n {
        b.set(i, b.get(i) - a.get(i));
    }
}

pub fn bfgs<F, X, C>(
    f: &mut F,
    x: &mut X,
    g: &mut C,
    max_steps: u32,
    _average_required_improvement: Fl,
    _over: usize,
    evalcount: &mut i32,
) -> Fl
where
    F: FnMut(&X, &mut C) -> Fl,
    X: Clone + IncrementBy<C>,
    C: ChangeVector,
{
    let n = g.num_floats();
    let mut h = TriangularMatrix::new(n, 0.0);
    set_diagonal(&mut h, 1.0);

    let mut g_new = g.clone();
    let mut x_new = x.clone();
    let mut f0 = f(x, g);
    *evalcount += 1;

    let f_orig = f0;
    let g_orig = g.clone();
    let x_orig = x.clone();
    let mut p = g.clone();

    for step in 0..max_steps {
        minus_mat_vec_product(&h, g, &mut p);
        let mut f1 = 0.0;
        let alpha = line_search(
            f, n, x, g, f0, &p, &mut x_new, &mut g_new, &mut f1, evalcount,
        );
        let mut y = g_new.clone();
        subtract_change(&mut y, g, n);

        f0 = f1;
        *x = x_new.clone();
        if scalar_product(g, g, n).sqrt() < 1e-5 {
            break;
        }
        *g = g_new.clone();

        if step == 0 {
            let yy = scalar_product(&y, &y, n);
            if yy.abs() > EPSILON_FL {
                set_diagonal(&mut h, alpha * scalar_product(&y, &p, n) / yy);
            }
        }
        bfgs_update(&mut h, &p, &y, alpha);
    }

    if f0 > f_orig {
        f0 = f_orig;
        *x = x_orig;
        *g = g_orig;
    }
    f0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conf::{Change, Conf, ConfSize};

    #[test]
    fn scalar_product_reads_change_flattening() {
        let size = ConfSize {
            ligands: vec![1],
            flex: vec![],
        };
        let mut a = Change::new(&size);
        let mut b = Change::new(&size);
        a.set(0, 2.0);
        b.set(0, 3.0);
        assert_eq!(scalar_product(&a, &b, a.num_floats()), 6.0);
    }

    #[test]
    fn bfgs_runs_on_quadratic_smoke_case() {
        let size = ConfSize {
            ligands: vec![],
            flex: vec![1],
        };
        let mut x = Conf::new(&size);
        x.flex[0].torsions[0] = 2.0;
        let mut g = Change::new(&size);
        let mut evals = 0;
        let mut f = |conf: &Conf, grad: &mut Change| {
            let value = conf.flex[0].torsions[0];
            grad.set(0, 2.0 * value);
            value * value
        };
        let result = bfgs(&mut f, &mut x, &mut g, 10, 0.0, 10, &mut evals);
        assert!(result <= 4.0);
        assert!(evals > 0);
    }
}
