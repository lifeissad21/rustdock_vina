use crate::common::Vec3;

// MT19937 and distributions follow Boost.Random, as used by reference Vina.
use crate::random_tables::{EXPONENTIAL_X, EXPONENTIAL_Y, NORMAL_X, NORMAL_Y};
#[derive(Debug, Clone)]
pub struct Rng64 {
    state: [u32; 624],
    index: usize,
}
impl Rng64 {
    pub fn new(seed: u64) -> Self {
        let mut state = [0u32; 624];
        state[0] = seed as u32;
        for i in 1..624 {
            state[i] = 1812433253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }
    pub fn next_u32(&mut self) -> u32 {
        if self.index == 624 {
            for i in 0..624 {
                let y = (self.state[i] & 0x80000000) | (self.state[(i + 1) % 624] & 0x7fffffff);
                self.state[i] = self.state[(i + 397) % 624]
                    ^ (y >> 1)
                    ^ if y & 1 != 0 { 0x9908b0df } else { 0 };
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }
    pub fn next_u64(&mut self) -> u64 {
        (u64::from(self.next_u32()) << 32) | u64::from(self.next_u32())
    }
    pub fn next_f64_open01(&mut self) -> f64 {
        self.uniform01().max(f64::MIN_POSITIVE)
    }
    fn uniform01(&mut self) -> f64 {
        f64::from(self.next_u32()) / 4294967296.0
    }
    fn int_float_pair(&mut self) -> (f64, usize) {
        let u = self.next_u32();
        let bucket = (u & 255) as usize;
        let r = f64::from(u >> 8) / 16777216.0;
        let remaining = self.next_u32() & ((1 << 29) - 1);
        ((r + f64::from(remaining)) / 536870912.0, bucket)
    }
}
pub fn random_fl(a: f64, b: f64, generator: &mut Rng64) -> f64 {
    assert!(a < b);
    loop {
        let value = generator.uniform01() * (b - a) + a;
        if value < b {
            return value;
        }
    }
}
fn exponential(generator: &mut Rng64) -> f64 {
    let mut shift = 0.0;
    loop {
        let (r, i) = generator.int_float_pair();
        let x = r * EXPONENTIAL_X[i];
        if x < EXPONENTIAL_X[i + 1] {
            return shift + x;
        }
        if i == 0 {
            shift += EXPONENTIAL_X[1];
            continue;
        }
        let y01 = generator.uniform01();
        let y = EXPONENTIAL_Y[i] + y01 * (EXPONENTIAL_Y[i + 1] - EXPONENTIAL_Y[i]);
        let upper = (EXPONENTIAL_X[i] - EXPONENTIAL_X[i + 1]) * y01 - (EXPONENTIAL_X[i] - x);
        let lower = y - (EXPONENTIAL_Y[i + 1] + (EXPONENTIAL_X[i + 1] - x) * EXPONENTIAL_Y[i + 1]);
        if upper < 0.0 && (lower < 0.0 || y < (-x).exp()) {
            return x + shift;
        }
    }
}
pub fn random_normal(mean: f64, sigma: f64, generator: &mut Rng64) -> f64 {
    assert!(sigma >= 0.0);
    loop {
        let (r, bucket) = generator.int_float_pair();
        let i = bucket >> 1;
        let sign = if bucket & 1 == 0 { -1.0 } else { 1.0 };
        let x = r * NORMAL_X[i];
        if x < NORMAL_X[i + 1] {
            return mean + sigma * x * sign;
        }
        if i == 0 {
            loop {
                let x = exponential(generator) / NORMAL_X[1];
                let y = exponential(generator);
                if 2.0 * y > x * x {
                    return mean + sigma * (x + NORMAL_X[1]) * sign;
                }
            }
        }
        let y01 = generator.uniform01();
        let y = NORMAL_Y[i] + y01 * (NORMAL_Y[i + 1] - NORMAL_Y[i]);
        let diagonal = (NORMAL_X[i] - NORMAL_X[i + 1]) * y01 - (NORMAL_X[i] - x);
        let tangent = y - (NORMAL_Y[i] + (NORMAL_X[i] - x) * NORMAL_Y[i] * NORMAL_X[i]);
        let (upper, lower) = if NORMAL_X[i] >= 1.0 {
            (diagonal, tangent)
        } else {
            (tangent, diagonal)
        };
        if upper < 0.0 && (lower < 0.0 || y < (-(x * x / 2.0)).exp()) {
            return mean + sigma * x * sign;
        }
    }
}
pub fn random_int(a: i32, b: i32, generator: &mut Rng64) -> i32 {
    assert!(a <= b);
    let width = (i64::from(b) - i64::from(a) + 1) as u64;
    if width == 1 {
        return a;
    }
    let bucket = 4294967296u64 / width;
    loop {
        let value = u64::from(generator.next_u32()) / bucket;
        if value < width {
            return (value as i64 + i64::from(a)) as i32;
        }
    }
}

pub fn random_sz(a: usize, b: usize, generator: &mut Rng64) -> usize {
    assert!(a <= b);
    random_int(a as i32, b as i32, generator) as usize
}

pub fn random_inside_sphere(generator: &mut Rng64) -> Vec3 {
    loop {
        let tmp = Vec3::new(
            random_fl(-1.0, 1.0, generator),
            random_fl(-1.0, 1.0, generator),
            random_fl(-1.0, 1.0, generator),
        );
        if tmp.norm_sqr() < 1.0 {
            return tmp;
        }
    }
}

pub fn random_in_box(corner1: Vec3, corner2: Vec3, generator: &mut Rng64) -> Vec3 {
    Vec3::new(
        random_fl(corner1[0], corner2[0], generator),
        random_fl(corner1[1], corner2[1], generator),
        random_fl(corner1[2], corner2[2], generator),
    )
}

pub fn auto_seed() -> u64 {
    use std::hash::{Hash, Hasher};
    use std::time::{SystemTime, UNIX_EPOCH};

    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or_default();
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    now.hash(&mut hasher);
    std::process::id().hash(&mut hasher);
    hasher.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mt19937_matches_reference_sequence() {
        let mut rng = Rng64::new(5489);
        for expected in [3499211612, 581869302, 3890346734, 3586334585, 545404204] {
            assert_eq!(rng.next_u32(), expected);
        }
    }

    #[test]
    fn boost_distributions_match_cpp_oracle() {
        let mut rng = Rng64::new(42);
        let expected = [
            (0.62178080016747117, 796712, -0.83680812597664633),
            (3.1239575683139265, 779856, -0.86923168979247489),
            (-0.9078695303760469, 445923, -0.75893171670396486),
            (-1.5934147224761546, 459342, 1.2299577240166031),
        ];
        for (real, integer, normal) in expected {
            assert!((random_fl(-2.0, 5.0, &mut rng) - real).abs() < 1e-14);
            assert_eq!(random_int(-10, 1000000, &mut rng), integer);
            assert!((random_normal(0.0, 1.0, &mut rng) - normal).abs() < 1e-14);
        }
        let mut rng = Rng64::new(7);
        let mut sum = 0.0;
        let mut square = 0.0;
        for _ in 0..10000 {
            let value = random_normal(0.0, 1.0, &mut rng);
            sum += value;
            square += value * value;
        }
        assert!((sum - (-39.743039124672883)).abs() < 1e-10);
        assert!((square - 9991.8533422858727).abs() < 1e-8);
    }

    #[test]
    fn generated_values_stay_in_requested_bounds() {
        let mut rng = Rng64::new(42);
        for _ in 0..100 {
            let value = random_fl(-2.0, 5.0, &mut rng);
            assert!((-2.0..=5.0).contains(&value));
        }
    }

    #[test]
    fn deterministic_seed_repeats_sequence() {
        let mut a = Rng64::new(7);
        let mut b = Rng64::new(7);
        assert_eq!(a.next_u64(), b.next_u64());
        assert_eq!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn inside_sphere_returns_valid_vector() {
        let mut rng = Rng64::new(11);
        let value = random_inside_sphere(&mut rng);
        assert!(value.norm_sqr() < 1.0);
    }
}
