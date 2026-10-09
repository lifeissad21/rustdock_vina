use crate::atom::Atom;
use crate::atom_constants::*;
use crate::atom_type::AtomType;
use crate::common::{not_max, sqr, Fl, EPSILON_FL, MAX_FL, PI};
use crate::int_pow::int_pow;

pub fn slope_step(x_bad: Fl, x_good: Fl, x: Fl) -> Fl {
    if x_bad < x_good {
        if x <= x_bad {
            return 0.0;
        }
        if x >= x_good {
            return 1.0;
        }
    } else {
        if x >= x_bad {
            return 0.0;
        }
        if x <= x_good {
            return 1.0;
        }
    }
    (x - x_bad) / (x_good - x_bad)
}

pub fn is_glue_type(xs_t: usize) -> bool {
    matches!(xs_t, XS_TYPE_G0 | XS_TYPE_G1 | XS_TYPE_G2 | XS_TYPE_G3)
}

pub fn optimal_distance(xs_t1: usize, xs_t2: usize) -> Fl {
    if is_glue_type(xs_t1) || is_glue_type(xs_t2) {
        0.0
    } else {
        xs_radius(xs_t1) + xs_radius(xs_t2)
    }
}

pub fn smooth_div(x: Fl, y: Fl) -> Fl {
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

pub fn optimal_distance_vinardo(xs_t1: usize, xs_t2: usize) -> Fl {
    if is_glue_type(xs_t1) || is_glue_type(xs_t2) {
        0.0
    } else {
        xs_vinardo_radius(xs_t1) + xs_vinardo_radius(xs_t2)
    }
}

pub fn smoothen(r: Fl, rij: Fl, mut smoothing: Fl) -> Fl {
    smoothing *= 0.5;
    if r > rij + smoothing {
        r - smoothing
    } else if r < rij - smoothing {
        r + smoothing
    } else {
        rij
    }
}

pub fn ad4_hb_eps(a: usize) -> Fl {
    ad_type_property(a).hb_depth
}

pub fn ad4_hb_radius(t: usize) -> Fl {
    ad_type_property(t).hb_radius
}

pub fn ad4_vdw_eps(a: usize) -> Fl {
    ad_type_property(a).depth
}

pub fn ad4_vdw_radius(t: usize) -> Fl {
    ad_type_property(t).radius
}

pub fn is_glued(xs_t1: usize, xs_t2: usize) -> bool {
    matches!(
        (xs_t1, xs_t2),
        (XS_TYPE_G0, XS_TYPE_C_H_CG0)
            | (XS_TYPE_G0, XS_TYPE_C_P_CG0)
            | (XS_TYPE_C_H_CG0, XS_TYPE_G0)
            | (XS_TYPE_C_P_CG0, XS_TYPE_G0)
            | (XS_TYPE_G1, XS_TYPE_C_H_CG1)
            | (XS_TYPE_G1, XS_TYPE_C_P_CG1)
            | (XS_TYPE_C_H_CG1, XS_TYPE_G1)
            | (XS_TYPE_C_P_CG1, XS_TYPE_G1)
            | (XS_TYPE_G2, XS_TYPE_C_H_CG2)
            | (XS_TYPE_G2, XS_TYPE_C_P_CG2)
            | (XS_TYPE_C_H_CG2, XS_TYPE_G2)
            | (XS_TYPE_C_P_CG2, XS_TYPE_G2)
            | (XS_TYPE_G3, XS_TYPE_C_H_CG3)
            | (XS_TYPE_G3, XS_TYPE_C_P_CG3)
            | (XS_TYPE_C_H_CG3, XS_TYPE_G3)
            | (XS_TYPE_C_P_CG3, XS_TYPE_G3)
    )
}

pub trait Potential {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl;
    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl;
    fn cutoff(&self) -> Fl;
}

#[derive(Debug, Clone, Copy)]
pub struct Gaussian {
    pub offset: Fl,
    pub width: Fl,
    pub cutoff: Fl,
    pub vinardo: bool,
}

impl Gaussian {
    fn gauss(self, x: Fl) -> Fl {
        (-(sqr(x / self.width))).exp()
    }
}

impl Potential for Gaussian {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        let xs_a = a.base.atom_type.xs;
        let xs_b = b.base.atom_type.xs;
        if r >= self.cutoff || xs_a >= XS_TYPE_SIZE || xs_b >= XS_TYPE_SIZE {
            return 0.0;
        }
        self.eval_types(xs_a, xs_b, r)
    }

    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let optimal = if self.vinardo {
            optimal_distance_vinardo(t1, t2)
        } else {
            optimal_distance(t1, t2)
        };
        self.gauss(r - (optimal + self.offset))
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Repulsion {
    pub offset: Fl,
    pub cutoff: Fl,
    pub vinardo: bool,
}

impl Potential for Repulsion {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        let xs_a = a.base.atom_type.xs;
        let xs_b = b.base.atom_type.xs;
        if r >= self.cutoff || xs_a >= XS_TYPE_SIZE || xs_b >= XS_TYPE_SIZE {
            return 0.0;
        }
        self.eval_types(xs_a, xs_b, r)
    }

    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let optimal = if self.vinardo {
            optimal_distance_vinardo(t1, t2)
        } else {
            optimal_distance(t1, t2)
        };
        let d = r - (optimal + self.offset);
        if d > 0.0 {
            0.0
        } else {
            d * d
        }
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Hydrophobic {
    pub good: Fl,
    pub bad: Fl,
    pub cutoff: Fl,
    pub vinardo: bool,
}

impl Potential for Hydrophobic {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        let xs_a = a.base.atom_type.xs;
        let xs_b = b.base.atom_type.xs;
        if r >= self.cutoff || xs_a >= XS_TYPE_SIZE || xs_b >= XS_TYPE_SIZE {
            return 0.0;
        }
        self.eval_types(xs_a, xs_b, r)
    }

    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        if r >= self.cutoff || !xs_is_hydrophobic(t1) || !xs_is_hydrophobic(t2) {
            return 0.0;
        }
        let optimal = if self.vinardo {
            optimal_distance_vinardo(t1, t2)
        } else {
            optimal_distance(t1, t2)
        };
        slope_step(self.bad, self.good, r - optimal)
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NonDirHBond {
    pub good: Fl,
    pub bad: Fl,
    pub cutoff: Fl,
    pub vinardo: bool,
}

impl Potential for NonDirHBond {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        let xs_a = a.base.atom_type.xs;
        let xs_b = b.base.atom_type.xs;
        if r >= self.cutoff || xs_a >= XS_TYPE_SIZE || xs_b >= XS_TYPE_SIZE {
            return 0.0;
        }
        self.eval_types(xs_a, xs_b, r)
    }

    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        if r >= self.cutoff || !xs_h_bond_possible(t1, t2) {
            return 0.0;
        }
        let optimal = if self.vinardo {
            optimal_distance_vinardo(t1, t2)
        } else {
            optimal_distance(t1, t2)
        };
        slope_step(self.bad, self.good, r - optimal)
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ad4Electrostatic {
    pub cap: Fl,
    pub cutoff: Fl,
}

impl Potential for Ad4Electrostatic {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let q1q2 = a.base.charge * b.base.charge * 332.0;
        let big_b = 78.4 + 8.5525;
        let lb = -big_b * 0.003627;
        let diel = -8.5525 + (big_b / (1.0 + 7.7839 * (lb * r).exp()));
        if r < EPSILON_FL {
            q1q2 * self.cap / diel
        } else {
            q1q2 * self.cap.min(1.0 / (r * diel))
        }
    }

    fn eval_types(&self, _t1: usize, _t2: usize, _r: Fl) -> Fl {
        0.0
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ad4Solvation {
    pub desolvation_sigma: Fl,
    pub solvation_q: Fl,
    pub charge_dependent: bool,
    pub cutoff: Fl,
}

impl Ad4Solvation {
    fn volume(atom_type: AtomType) -> Fl {
        if atom_type.ad < AD_TYPE_SIZE {
            ad_type_property(atom_type.ad).volume
        } else if atom_type.xs < XS_TYPE_SIZE {
            4.0 * PI / 3.0 * int_pow::<3>(xs_radius(atom_type.xs))
        } else {
            panic!("atom type lacks volume")
        }
    }

    fn solvation_parameter(atom_type: AtomType) -> Fl {
        if atom_type.ad < AD_TYPE_SIZE {
            ad_type_property(atom_type.ad).solvation
        } else if atom_type.xs == XS_TYPE_MET_D {
            METAL_SOLVATION_PARAMETER
        } else {
            panic!("atom type lacks solvation parameter")
        }
    }
}

impl Potential for Ad4Solvation {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let q1 = a.base.charge;
        let q2 = b.base.charge;
        assert!(not_max(q1));
        assert!(not_max(q2));
        let solv1 = Self::solvation_parameter(a.base.atom_type);
        let solv2 = Self::solvation_parameter(b.base.atom_type);
        let volume1 = Self::volume(a.base.atom_type);
        let volume2 = Self::volume(b.base.atom_type);
        let my_solv = if self.charge_dependent {
            self.solvation_q
        } else {
            0.0
        };
        let tmp = ((solv1 + my_solv * q1.abs()) * volume2 + (solv2 + my_solv * q2.abs()) * volume1)
            * (-0.5 * sqr(r / self.desolvation_sigma)).exp();
        assert!(not_max(tmp));
        tmp
    }

    fn eval_types(&self, _t1: usize, _t2: usize, _r: Fl) -> Fl {
        0.0
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ad4Vdw {
    pub smoothing: Fl,
    pub cap: Fl,
    pub cutoff: Fl,
}

impl Potential for Ad4Vdw {
    fn eval_atoms(&self, a: &Atom, b: &Atom, mut r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let t1 = a.base.atom_type.ad;
        let t2 = b.base.atom_type.ad;
        let hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        let vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        let vdw_depth = (ad4_vdw_eps(t1) * ad4_vdw_eps(t2)).sqrt();
        if hb_depth < 0.0 {
            return 0.0;
        }
        r = smoothen(r, vdw_rij, self.smoothing);
        let c_12 = int_pow::<12>(vdw_rij) * vdw_depth;
        let c_6 = int_pow::<6>(vdw_rij) * vdw_depth * 2.0;
        let r6 = int_pow::<6>(r);
        let r12 = int_pow::<12>(r);
        if r12 > EPSILON_FL && r6 > EPSILON_FL {
            self.cap.min(c_12 / r12 - c_6 / r6)
        } else {
            self.cap
        }
    }

    fn eval_types(&self, _t1: usize, _t2: usize, _r: Fl) -> Fl {
        0.0
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ad4Hb {
    pub smoothing: Fl,
    pub cap: Fl,
    pub cutoff: Fl,
}

impl Potential for Ad4Hb {
    fn eval_atoms(&self, a: &Atom, b: &Atom, mut r: Fl) -> Fl {
        if r >= self.cutoff {
            return 0.0;
        }
        let t1 = a.base.atom_type.ad;
        let t2 = b.base.atom_type.ad;
        let hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2);
        let hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        if hb_depth >= 0.0 {
            return 0.0;
        }
        r = smoothen(r, hb_rij, self.smoothing);
        let c_12 = int_pow::<12>(hb_rij) * -hb_depth * 10.0 / 2.0;
        let c_10 = int_pow::<10>(hb_rij) * -hb_depth * 12.0 / 2.0;
        let r10 = int_pow::<10>(r);
        let r12 = int_pow::<12>(r);
        if r12 > EPSILON_FL && r10 > EPSILON_FL {
            self.cap.min(c_12 / r12 - c_10 / r10)
        } else {
            self.cap
        }
    }

    fn eval_types(&self, _t1: usize, _t2: usize, _r: Fl) -> Fl {
        0.0
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LinearAttraction {
    pub cutoff: Fl,
}

impl Potential for LinearAttraction {
    fn eval_atoms(&self, a: &Atom, b: &Atom, r: Fl) -> Fl {
        self.eval_types(a.base.atom_type.xs, b.base.atom_type.xs, r)
    }

    fn eval_types(&self, t1: usize, t2: usize, r: Fl) -> Fl {
        if r >= self.cutoff {
            0.0
        } else if is_glued(t1, t2) {
            r
        } else {
            0.0
        }
    }

    fn cutoff(&self) -> Fl {
        self.cutoff
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_fl;

    #[test]
    fn glue_types_have_zero_optimal_distance() {
        assert_eq!(optimal_distance(XS_TYPE_G0, XS_TYPE_C_H), 0.0);
        assert_eq!(xs_radius(XS_TYPE_G0), 0.0);
    }

    #[test]
    fn vina_repulsion_matches_reference_formula() {
        let repulsion = Repulsion {
            offset: 0.0,
            cutoff: 8.0,
            vinardo: false,
        };
        let optimal = optimal_distance(XS_TYPE_C_H, XS_TYPE_C_H);
        assert!(eq_fl(
            repulsion.eval_types(XS_TYPE_C_H, XS_TYPE_C_H, optimal - 1.0),
            1.0
        ));
        assert_eq!(
            repulsion.eval_types(XS_TYPE_C_H, XS_TYPE_C_H, optimal + 1.0),
            0.0
        );
    }

    #[test]
    fn macrocycle_linear_attraction_only_for_glued_pairs() {
        let attraction = LinearAttraction { cutoff: 10.0 };
        assert_eq!(attraction.eval_types(XS_TYPE_G0, XS_TYPE_C_H_CG0, 2.0), 2.0);
        assert_eq!(attraction.eval_types(XS_TYPE_G0, XS_TYPE_C_H, 2.0), 0.0);
    }
}
