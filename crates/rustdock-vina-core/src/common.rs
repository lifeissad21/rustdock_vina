use core::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

pub type Fl = f64;
pub type Sz = usize;
pub type Pair = (Fl, Fl);

pub const PI: Fl = std::f64::consts::PI;
pub const FL_TOLERANCE: Fl = 0.001;
pub const EPSILON_FL: Fl = Fl::EPSILON;
pub const MAX_FL: Fl = Fl::MAX;
pub const MAX_SZ: Sz = Sz::MAX;
pub const PK_TO_ENERGY_FACTOR: Fl = -8.31 * 0.001 * 300.0 / 4.184 * core::f64::consts::LN_10;
pub const ZERO_VEC: Vec3 = Vec3::new(0.0, 0.0, 0.0);
pub const MAX_VEC: Vec3 = Vec3::new(MAX_FL, MAX_FL, MAX_FL);

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub data: [Fl; 3],
}

impl Vec3 {
    pub const fn new(x: Fl, y: Fl, z: Fl) -> Self {
        Self { data: [x, y, z] }
    }

    pub fn splat(value: Fl) -> Self {
        Self::new(value, value, value)
    }

    pub fn norm_sqr(self) -> Fl {
        sqr(self.data[0]) + sqr(self.data[1]) + sqr(self.data[2])
    }

    pub fn norm(self) -> Fl {
        self.norm_sqr().sqrt()
    }

    pub fn dot(self, other: Self) -> Fl {
        self.data[0] * other[0] + self.data[1] * other[1] + self.data[2] * other[2]
    }

    pub const fn len(self) -> usize {
        3
    }

    pub const fn is_empty(self) -> bool {
        false
    }
}

impl Index<usize> for Vec3 {
    type Output = Fl;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2])
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2])
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        self.data[0] += rhs[0];
        self.data[1] += rhs[1];
        self.data[2] += rhs[2];
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        self.data[0] -= rhs[0];
        self.data[1] -= rhs[1];
        self.data[2] -= rhs[2];
    }
}

impl Mul<Fl> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: Fl) -> Self::Output {
        Self::new(self[0] * rhs, self[1] * rhs, self[2] * rhs)
    }
}

impl Mul<Vec3> for Fl {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<Fl> for Vec3 {
    fn mul_assign(&mut self, rhs: Fl) {
        self.data[0] *= rhs;
        self.data[1] *= rhs;
        self.data[2] *= rhs;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    data: [Fl; 9],
}

impl Mat3 {
    pub const fn new(
        xx: Fl,
        xy: Fl,
        xz: Fl,
        yx: Fl,
        yy: Fl,
        yz: Fl,
        zx: Fl,
        zy: Fl,
        zz: Fl,
    ) -> Self {
        Self {
            data: [xx, yx, zx, xy, yy, zy, xz, yz, zz],
        }
    }

    pub const fn zero() -> Self {
        Self { data: [0.0; 9] }
    }

    pub fn get(&self, i: usize, j: usize) -> Fl {
        debug_assert!(i < 3);
        debug_assert!(j < 3);
        self.data[i + 3 * j]
    }

    pub fn set(&mut self, i: usize, j: usize, value: Fl) {
        debug_assert!(i < 3);
        debug_assert!(j < 3);
        self.data[i + 3 * j] = value;
    }

    pub fn mul_vec(self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.data[0] * v[0] + self.data[3] * v[1] + self.data[6] * v[2],
            self.data[1] * v[0] + self.data[4] * v[1] + self.data[7] * v[2],
            self.data[2] * v[0] + self.data[5] * v[1] + self.data[8] * v[2],
        )
    }
}

impl MulAssign<Fl> for Mat3 {
    fn mul_assign(&mut self, rhs: Fl) {
        for value in &mut self.data {
            *value *= rhs;
        }
    }
}

pub fn sqr<T>(x: T) -> T
where
    T: Copy + Mul<Output = T>,
{
    x * x
}

pub fn cross_product(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )
}

pub fn elementwise_product(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a[0] * b[0], a[1] * b[1], a[2] * b[2])
}

pub fn fl_to_sz(x: Fl, max_sz: Sz) -> Sz {
    if x <= 0.0 {
        0
    } else if x >= max_sz as Fl {
        max_sz
    } else {
        (x as Sz).min(max_sz)
    }
}

pub fn eq_fl(a: Fl, b: Fl) -> bool {
    (a - b).abs() < FL_TOLERANCE
}

pub fn eq_vec(a: Vec3, b: Vec3) -> bool {
    eq_fl(a[0], b[0]) && eq_fl(a[1], b[1]) && eq_fl(a[2], b[2])
}

pub fn not_max(x: Fl) -> bool {
    x < 0.1 * MAX_FL
}

pub fn vec_distance_sqr(a: Vec3, b: Vec3) -> Fl {
    sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2])
}

pub fn normalize_angle(x: &mut Fl) {
    if *x > 3.0 * PI {
        let n = (*x - PI) / (2.0 * PI);
        *x -= 2.0 * PI * n.ceil();
        normalize_angle(x);
    } else if *x < -3.0 * PI {
        let n = (-*x - PI) / (2.0 * PI);
        *x += 2.0 * PI * n.ceil();
        normalize_angle(x);
    } else if *x > PI {
        *x -= 2.0 * PI;
    } else if *x < -PI {
        *x += 2.0 * PI;
    }
    debug_assert!(*x >= -PI && *x <= PI);
}

pub fn normalized_angle(mut x: Fl) -> Fl {
    normalize_angle(&mut x);
    x
}

pub fn sum<T>(values: &[T]) -> T
where
    T: Copy + Default + AddAssign,
{
    let mut acc = T::default();
    for value in values {
        acc += *value;
    }
    acc
}

pub fn pk_to_energy(pk: Fl) -> Fl {
    PK_TO_ENERGY_FACTOR * pk
}

pub fn starts_with(value: &str, start: &str) -> bool {
    value.starts_with(start)
}

pub fn has<T: PartialEq>(values: &[T], element: &T) -> bool {
    values.contains(element)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InternalError {
    pub file: &'static str,
    pub line: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_math_matches_reference_layout() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.norm_sqr(), 14.0);
        assert_eq!(v.dot(Vec3::new(4.0, 5.0, 6.0)), 32.0);
        assert_eq!(
            cross_product(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
            Vec3::new(0.0, 0.0, 1.0)
        );
    }

    #[test]
    fn matrix_is_column_major_like_cpp() {
        let m = Mat3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(0, 1), 2.0);
        assert_eq!(m.get(2, 2), 9.0);
        assert_eq!(
            m.mul_vec(Vec3::new(1.0, 1.0, 1.0)),
            Vec3::new(6.0, 15.0, 24.0)
        );
    }

    #[test]
    fn angle_normalization_matches_cpp_bounds() {
        assert!(eq_fl(normalized_angle(4.0 * PI), 0.0));
        assert!(eq_fl(normalized_angle(-4.0 * PI), 0.0));
        assert!(normalized_angle(1.5 * PI) < 0.0);
    }
}
