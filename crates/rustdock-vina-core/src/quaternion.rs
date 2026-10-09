use crate::common::{eq_fl, normalized_angle, Mat3, Vec3, EPSILON_FL, PI, ZERO_VEC};
use crate::random::{random_normal, Rng64};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    pub r: f64,
    pub i: f64,
    pub j: f64,
    pub k: f64,
}

pub const QT_IDENTITY: Quaternion = Quaternion::new(1.0, 0.0, 0.0, 0.0);

impl Quaternion {
    pub const fn new(r: f64, i: f64, j: f64, k: f64) -> Self {
        Self { r, i, j, k }
    }

    pub fn norm_sqr(self) -> f64 {
        self.r * self.r + self.i * self.i + self.j * self.j + self.k * self.k
    }

    pub fn norm(self) -> f64 {
        self.norm_sqr().sqrt()
    }

    pub fn normalized(mut self) -> Self {
        let norm = self.norm();
        assert!(norm > EPSILON_FL);
        self.r /= norm;
        self.i /= norm;
        self.j /= norm;
        self.k /= norm;
        self
    }

    pub fn normalize_approx(&mut self, tolerance: f64) {
        let s = self.norm_sqr();
        if (s - 1.0).abs() >= tolerance {
            *self = self.normalized();
        }
    }

    pub fn conjugate(self) -> Self {
        Self::new(self.r, -self.i, -self.j, -self.k)
    }
}

impl core::ops::Mul for Quaternion {
    type Output = Quaternion;

    fn mul(self, rhs: Self) -> Self::Output {
        Quaternion::new(
            self.r * rhs.r - self.i * rhs.i - self.j * rhs.j - self.k * rhs.k,
            self.r * rhs.i + self.i * rhs.r + self.j * rhs.k - self.k * rhs.j,
            self.r * rhs.j - self.i * rhs.k + self.j * rhs.r + self.k * rhs.i,
            self.r * rhs.k + self.i * rhs.j - self.j * rhs.i + self.k * rhs.r,
        )
    }
}

pub fn quaternion_approx_eq(a: Quaternion, b: Quaternion) -> bool {
    eq_fl(a.r, b.r) && eq_fl(a.i, b.i) && eq_fl(a.j, b.j) && eq_fl(a.k, b.k)
}

pub fn quaternion_is_normalized(q: Quaternion) -> bool {
    eq_fl(q.norm_sqr(), 1.0) && eq_fl(q.norm(), 1.0)
}

pub fn angle_axis_to_quaternion(axis: Vec3, angle: f64) -> Quaternion {
    assert!(eq_fl(axis.norm(), 1.0));
    let angle = normalized_angle(angle);
    let c = (angle / 2.0).cos();
    let s = (angle / 2.0).sin();
    Quaternion::new(c, s * axis[0], s * axis[1], s * axis[2])
}

pub fn rotation_to_quaternion(rotation: Vec3) -> Quaternion {
    let angle = rotation.norm();
    if angle > EPSILON_FL {
        angle_axis_to_quaternion((1.0 / angle) * rotation, angle)
    } else {
        QT_IDENTITY
    }
}

pub fn quaternion_to_angle(q: Quaternion) -> Vec3 {
    assert!(quaternion_is_normalized(q));
    let c = q.r;
    if c > -1.0 && c < 1.0 {
        let mut angle = 2.0 * c.acos();
        if angle > PI {
            angle -= 2.0 * PI;
        }
        let mut axis = Vec3::new(q.i, q.j, q.k);
        let s = (angle / 2.0).sin();
        if s.abs() < EPSILON_FL {
            return ZERO_VEC;
        }
        axis *= angle / s;
        axis
    } else {
        ZERO_VEC
    }
}

pub fn quaternion_to_r3(q: Quaternion) -> Mat3 {
    assert!(quaternion_is_normalized(q));

    let a = q.r;
    let b = q.i;
    let c = q.j;
    let d = q.k;

    let aa = a * a;
    let ab = a * b;
    let ac = a * c;
    let ad = a * d;
    let bb = b * b;
    let bc = b * c;
    let bd = b * d;
    let cc = c * c;
    let cd = c * d;
    let dd = d * d;

    assert!(eq_fl(aa + bb + cc + dd, 1.0));

    Mat3::new(
        aa + bb - cc - dd,
        2.0 * (-ad + bc),
        2.0 * (ac + bd),
        2.0 * (ad + bc),
        aa - bb + cc - dd,
        2.0 * (-ab + cd),
        2.0 * (-ac + bd),
        2.0 * (ab + cd),
        aa - bb - cc + dd,
    )
}

pub fn random_orientation(generator: &mut Rng64) -> Quaternion {
    loop {
        let q = Quaternion::new(
            random_normal(0.0, 1.0, generator),
            random_normal(0.0, 1.0, generator),
            random_normal(0.0, 1.0, generator),
            random_normal(0.0, 1.0, generator),
        );
        if q.norm() > EPSILON_FL {
            return q.normalized();
        }
    }
}

pub fn quaternion_increment(q: &mut Quaternion, rotation: Vec3) {
    assert!(quaternion_is_normalized(*q));
    *q = rotation_to_quaternion(rotation) * *q;
    q.normalize_approx(1e-6);
}

pub fn quaternion_difference(b: Quaternion, a: Quaternion) -> Vec3 {
    assert!(quaternion_is_normalized(a));
    assert!(quaternion_is_normalized(b));
    quaternion_to_angle(b * a.conjugate())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_vec;

    #[test]
    fn identity_roundtrips_to_zero_angle() {
        assert!(eq_vec(quaternion_to_angle(QT_IDENTITY), ZERO_VEC));
        assert!(quaternion_is_normalized(QT_IDENTITY));
    }

    #[test]
    fn axis_angle_roundtrip_matches_reference_behavior() {
        let q = angle_axis_to_quaternion(Vec3::new(0.0, 0.0, 1.0), PI / 2.0);
        let angle = quaternion_to_angle(q);
        assert!(eq_fl(angle[0], 0.0));
        assert!(eq_fl(angle[1], 0.0));
        assert!(eq_fl(angle[2], PI / 2.0));
    }

    #[test]
    fn increment_applies_left_multiplication() {
        let mut q = QT_IDENTITY;
        quaternion_increment(&mut q, Vec3::new(0.0, 0.0, PI / 2.0));
        assert!(eq_fl(quaternion_to_angle(q)[2], PI / 2.0));
    }
}
