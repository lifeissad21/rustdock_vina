use crate::common::{eq_fl, Fl, Vec3};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridDim {
    pub begin: Fl,
    pub end: Fl,
    pub n_voxels: usize,
}

impl GridDim {
    pub const fn new(begin: Fl, end: Fl, n_voxels: usize) -> Self {
        Self {
            begin,
            end,
            n_voxels,
        }
    }

    pub const fn disabled() -> Self {
        Self::new(0.0, 0.0, 0)
    }

    pub fn span(self) -> Fl {
        self.end - self.begin
    }

    pub fn enabled(self) -> bool {
        self.n_voxels > 0
    }

    pub fn approx_eq(self, other: Self) -> bool {
        self.n_voxels == other.n_voxels
            && eq_fl(self.begin, other.begin)
            && eq_fl(self.end, other.end)
    }
}

pub type GridDims = [GridDim; 3];

pub fn grid_dims_approx_eq(a: &GridDims, b: &GridDims) -> bool {
    a[0].approx_eq(b[0]) && a[1].approx_eq(b[1]) && a[2].approx_eq(b[2])
}

pub fn grid_dims_begin(gd: &GridDims) -> Vec3 {
    Vec3::new(gd[0].begin, gd[1].begin, gd[2].begin)
}

pub fn grid_dims_end(gd: &GridDims) -> Vec3 {
    Vec3::new(gd[0].end, gd[1].end, gd[2].end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimensions_expose_reference_helpers() {
        let dims = [
            GridDim::new(-1.0, 1.0, 4),
            GridDim::new(-2.0, 2.0, 8),
            GridDim::new(-3.0, 3.0, 12),
        ];
        assert_eq!(dims[0].span(), 2.0);
        assert!(dims[0].enabled());
        assert_eq!(grid_dims_begin(&dims), Vec3::new(-1.0, -2.0, -3.0));
        assert_eq!(grid_dims_end(&dims), Vec3::new(1.0, 2.0, 3.0));
    }
}
