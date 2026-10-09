use crate::common::{not_max, sqr, Fl, EPSILON_FL};

pub trait CurlDerivative {
    fn scale_by(&mut self, value: Fl);
}

impl CurlDerivative for Fl {
    fn scale_by(&mut self, value: Fl) {
        *self *= value;
    }
}

impl CurlDerivative for crate::common::Vec3 {
    fn scale_by(&mut self, value: Fl) {
        *self *= value;
    }
}

pub fn curl_with_derivative<T: CurlDerivative>(energy: &mut Fl, derivative: &mut T, v: Fl) {
    if *energy > 0.0 && not_max(v) {
        let tmp = if v < EPSILON_FL {
            0.0
        } else {
            v / (v + *energy)
        };
        *energy *= tmp;
        derivative.scale_by(sqr(tmp));
    }
}

pub fn curl(energy: &mut Fl, v: Fl) {
    if *energy > 0.0 && not_max(v) {
        let tmp = if v < EPSILON_FL {
            0.0
        } else {
            v / (v + *energy)
        };
        *energy *= tmp;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{eq_fl, Vec3};

    #[test]
    fn curls_positive_energy_and_derivative() {
        let mut energy = 2.0;
        let mut deriv = Vec3::new(1.0, 2.0, 3.0);
        curl_with_derivative(&mut energy, &mut deriv, 2.0);
        assert!(eq_fl(energy, 1.0));
        assert!(eq_fl(deriv[0], 0.25));
        assert!(eq_fl(deriv[1], 0.5));
        assert!(eq_fl(deriv[2], 0.75));
    }
}
