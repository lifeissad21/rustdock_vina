use crate::array3d::Array3d;
use crate::common::{elementwise_product, fl_to_sz, Fl, Vec3, EPSILON_FL};
use crate::curl::{curl, curl_with_derivative};
use crate::grid_dim::GridDims;

#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    pub init: Vec3,
    pub range: Vec3,
    pub factor_inv: Vec3,
    pub data: Array3d<Fl>,
    factor: Vec3,
    dim_fl_minus_1: Vec3,
}

impl Default for Grid {
    fn default() -> Self {
        Self {
            init: Vec3::new(0.0, 0.0, 0.0),
            range: Vec3::new(1.0, 1.0, 1.0),
            factor: Vec3::new(1.0, 1.0, 1.0),
            dim_fl_minus_1: Vec3::new(-1.0, -1.0, -1.0),
            factor_inv: Vec3::new(1.0, 1.0, 1.0),
            data: Array3d::default(),
        }
    }
}

impl Grid {
    pub fn new(gd: &GridDims) -> Self {
        let mut grid = Grid::default();
        grid.init(gd);
        grid
    }

    pub fn init(&mut self, gd: &GridDims) {
        self.data = Array3d::new(
            gd[0].n_voxels + 1,
            gd[1].n_voxels + 1,
            gd[2].n_voxels + 1,
            0.0,
        )
        .expect("grid dimensions should fit usize");
        self.init = Vec3::new(gd[0].begin, gd[1].begin, gd[2].begin);
        self.range = Vec3::new(gd[0].span(), gd[1].span(), gd[2].span());
        assert!(self.range[0] > 0.0);
        assert!(self.range[1] > 0.0);
        assert!(self.range[2] > 0.0);
        self.dim_fl_minus_1 = Vec3::new(
            (self.data.dim0() - 1) as Fl,
            (self.data.dim1() - 1) as Fl,
            (self.data.dim2() - 1) as Fl,
        );
        for i in 0..3 {
            self.factor[i] = self.dim_fl_minus_1[i] / self.range[i];
            self.factor_inv[i] = 1.0 / self.factor[i];
        }
    }

    pub fn index_to_argument(&self, x: usize, y: usize, z: usize) -> Vec3 {
        Vec3::new(
            self.init[0] + self.factor_inv[0] * x as Fl,
            self.init[1] + self.factor_inv[1] * y as Fl,
            self.init[2] + self.factor_inv[2] * z as Fl,
        )
    }

    pub fn initialized(&self) -> bool {
        self.data.dim0() > 0 && self.data.dim1() > 0 && self.data.dim2() > 0
    }

    pub fn evaluate(&self, location: Vec3, slope: Fl, c: Fl) -> Fl {
        self.evaluate_aux(location, slope, c, None)
    }

    pub fn evaluate_with_deriv(&self, location: Vec3, slope: Fl, c: Fl, deriv: &mut Vec3) -> Fl {
        self.evaluate_aux(location, slope, c, Some(deriv))
    }

    fn evaluate_aux(&self, location: Vec3, slope: Fl, v: Fl, deriv: Option<&mut Vec3>) -> Fl {
        let mut s = elementwise_product(location - self.init, self.factor);
        let mut miss = Vec3::new(0.0, 0.0, 0.0);
        let mut region = [0_i32; 3];
        let mut a = [0_usize; 3];

        for i in 0..3 {
            if s[i] < 0.0 {
                miss[i] = -s[i];
                region[i] = -1;
                a[i] = 0;
                s[i] = 0.0;
            } else if s[i] >= self.dim_fl_minus_1[i] {
                miss[i] = s[i] - self.dim_fl_minus_1[i];
                region[i] = 1;
                assert!(self.data.dim(i) >= 2);
                a[i] = self.data.dim(i) - 2;
                s[i] = 1.0;
            } else {
                region[i] = 0;
                a[i] = fl_to_sz(s[i], usize::MAX);
                s[i] -= a[i] as Fl;
            }
            assert!(s[i] >= 0.0);
            assert!(s[i] <= 1.0);
            assert!(a[i] + 1 < self.data.dim(i));
        }

        let penalty = slope * miss.dot(self.factor_inv);
        assert!(penalty > -EPSILON_FL);

        let (x0, y0, z0) = (a[0], a[1], a[2]);
        let (x1, y1, z1) = (x0 + 1, y0 + 1, z0 + 1);
        let f000 = *self.data.get(x0, y0, z0);
        let f100 = *self.data.get(x1, y0, z0);
        let f010 = *self.data.get(x0, y1, z0);
        let f110 = *self.data.get(x1, y1, z0);
        let f001 = *self.data.get(x0, y0, z1);
        let f101 = *self.data.get(x1, y0, z1);
        let f011 = *self.data.get(x0, y1, z1);
        let f111 = *self.data.get(x1, y1, z1);

        let (x, y, z) = (s[0], s[1], s[2]);
        let (mx, my, mz) = (1.0 - x, 1.0 - y, 1.0 - z);

        let mut f = f000 * mx * my * mz
            + f100 * x * my * mz
            + f010 * mx * y * mz
            + f110 * x * y * mz
            + f001 * mx * my * z
            + f101 * x * my * z
            + f011 * mx * y * z
            + f111 * x * y * z;

        if let Some(deriv) = deriv {
            let x_g = f000 * -my * mz
                + f100 * my * mz
                + f010 * -y * mz
                + f110 * y * mz
                + f001 * -my * z
                + f101 * my * z
                + f011 * -y * z
                + f111 * y * z;
            let y_g = f000 * mx * -mz
                + f100 * x * -mz
                + f010 * mx * mz
                + f110 * x * mz
                + f001 * mx * -z
                + f101 * x * -z
                + f011 * mx * z
                + f111 * x * z;
            let z_g = f000 * mx * my * -1.0
                + f100 * x * my * -1.0
                + f010 * mx * y * -1.0
                + f110 * x * y * -1.0
                + f001 * mx * my
                + f101 * x * my
                + f011 * mx * y
                + f111 * x * y;
            let mut gradient = Vec3::new(x_g, y_g, z_g);
            curl_with_derivative(&mut f, &mut gradient, v);
            for i in 0..3 {
                let gradient_everywhere = if region[i] == 0 { gradient[i] } else { 0.0 };
                deriv[i] = self.factor[i] * gradient_everywhere + slope * region[i] as Fl;
            }
            f + penalty
        } else {
            curl(&mut f, v);
            f + penalty
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{eq_fl, MAX_FL};
    use crate::grid_dim::GridDim;

    #[test]
    fn trilinear_interpolation_matches_corner_average() {
        let dims = [
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
            GridDim::new(0.0, 1.0, 1),
        ];
        let mut grid = Grid::new(&dims);
        *grid.data.get_mut(0, 0, 0) = 0.0;
        *grid.data.get_mut(1, 0, 0) = 1.0;
        *grid.data.get_mut(0, 1, 0) = 1.0;
        *grid.data.get_mut(1, 1, 0) = 2.0;
        *grid.data.get_mut(0, 0, 1) = 1.0;
        *grid.data.get_mut(1, 0, 1) = 2.0;
        *grid.data.get_mut(0, 1, 1) = 2.0;
        *grid.data.get_mut(1, 1, 1) = 3.0;
        assert!(eq_fl(
            grid.evaluate(Vec3::new(0.5, 0.5, 0.5), 0.0, MAX_FL),
            1.5
        ));
    }
}
