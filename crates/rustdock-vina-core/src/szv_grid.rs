use crate::array3d::Array3d;
use crate::atom::Atom;
use crate::atom_type::{num_atom_types, AtomTyping};
use crate::brick::brick_distance_sqr;
use crate::common::{fl_to_sz, Fl, Vec3, EPSILON_FL};
use crate::grid_dim::{GridDim, GridDims};

pub trait SzvGridModel {
    fn atom_typing_used(&self) -> AtomTyping;
    fn grid_atoms(&self) -> &[Atom];
}

#[derive(Debug, Clone, PartialEq)]
pub struct SzvGrid {
    data: Array3d<Vec<usize>>,
    init: Vec3,
    range: Vec3,
}

impl Default for SzvGrid {
    fn default() -> Self {
        Self {
            data: Array3d::default(),
            init: Vec3::new(0.0, 0.0, 0.0),
            range: Vec3::new(0.0, 0.0, 0.0),
        }
    }
}

impl SzvGrid {
    pub fn new<M: SzvGridModel>(model: &M, gd: &GridDims, cutoff_sqr: Fl) -> Self {
        let data = Array3d::new(gd[0].n_voxels, gd[1].n_voxels, gd[2].n_voxels, Vec::new())
            .expect("szv grid dimensions should fit usize");
        let init = Vec3::new(gd[0].begin, gd[1].begin, gd[2].begin);
        let end = Vec3::new(gd[0].end, gd[1].end, gd[2].end);
        let range = end - init;
        let nat = num_atom_types(model.atom_typing_used());

        let relevant_indexes = model
            .grid_atoms()
            .iter()
            .enumerate()
            .filter_map(|(i, atom)| {
                let atom_type = atom.base.atom_type.get(model.atom_typing_used());
                if atom_type < nat && brick_distance_sqr(init, end, atom.coords) < cutoff_sqr {
                    Some(i)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        let mut grid = Self { data, init, range };
        for x in 0..grid.data.dim0() {
            for y in 0..grid.data.dim1() {
                for z in 0..grid.data.dim2() {
                    let begin = grid.index_to_coord(x, y, z);
                    let end = grid.index_to_coord(x + 1, y + 1, z + 1);
                    for i in &relevant_indexes {
                        if brick_distance_sqr(begin, end, model.grid_atoms()[*i].coords)
                            < cutoff_sqr
                        {
                            grid.data.get_mut(x, y, z).push(*i);
                        }
                    }
                }
            }
        }
        grid
    }

    pub fn average_num_possibilities(&self) -> Fl {
        let mut counter = 0;
        for x in 0..self.data.dim0() {
            for y in 0..self.data.dim1() {
                for z in 0..self.data.dim2() {
                    counter += self.data.get(x, y, z).len();
                }
            }
        }
        counter as Fl / (self.data.dim0() * self.data.dim1() * self.data.dim2()) as Fl
    }

    pub fn possibilities(&self, coords: Vec3) -> &[usize] {
        let mut index = [0_usize; 3];
        for i in 0..3 {
            // Reconstructing the upper bound can lose an ulp at translated boxes.
            let tolerance = 4.0
                * EPSILON_FL
                * self.init[i]
                    .abs()
                    .max((self.init[i] + self.range[i]).abs())
                    .max(1.0);
            assert!(coords[i] + tolerance >= self.init[i]);
            assert!(coords[i] <= self.init[i] + self.range[i] + tolerance);
            let tmp = (coords[i] - self.init[i]) * self.data.dim(i) as Fl / self.range[i];
            index[i] = fl_to_sz(tmp, self.data.dim(i) - 1);
        }
        self.data.get(index[0], index[1], index[2])
    }

    pub fn index_to_coord(&self, i: usize, j: usize, k: usize) -> Vec3 {
        let index = Vec3::new(i as Fl, j as Fl, k as Fl);
        Vec3::new(
            self.init[0] + self.range[0] * index[0] / self.data.dim(0) as Fl,
            self.init[1] + self.range[1] * index[1] / self.data.dim(1) as Fl,
            self.init[2] + self.range[2] * index[2] / self.data.dim(2) as Fl,
        )
    }
}

pub fn szv_grid_dims(gd: &GridDims) -> GridDims {
    [single_dim(gd[0]), single_dim(gd[1]), single_dim(gd[2])]
}

fn single_dim(dim: GridDim) -> GridDim {
    let n_fl = (dim.end - dim.begin) / 3.0;
    let n_int = n_fl as i32;
    GridDim::new(
        dim.begin,
        dim.end,
        if n_int < 1 { 1 } else { n_int as usize },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::atom_constants::{AD_TYPE_C, EL_TYPE_C, XS_TYPE_C_H};
    use crate::atom_type::AtomType;

    struct FakeModel {
        atoms: Vec<Atom>,
    }

    impl SzvGridModel for FakeModel {
        fn atom_typing_used(&self) -> AtomTyping {
            AtomTyping::Xs
        }

        fn grid_atoms(&self) -> &[Atom] {
            &self.atoms
        }
    }

    #[test]
    fn translated_upper_boundary_survives_rounding() {
        let begin = -10.759111;
        let end = 7.240889;
        assert!(begin + (end - begin) < end);
        let dims = [GridDim::new(begin, end, 6); 3];
        let grid = SzvGrid::new(&FakeModel { atoms: vec![] }, &dims, 4.0);
        assert!(grid.possibilities(Vec3::new(end, end, end)).is_empty());
        assert!(grid
            .possibilities(Vec3::new(begin, begin, begin))
            .is_empty());
    }

    #[test]
    #[should_panic]
    fn still_rejects_coordinates_outside_the_box() {
        let dims = [GridDim::new(-10.759111, 7.240889, 6); 3];
        let grid = SzvGrid::new(&FakeModel { atoms: vec![] }, &dims, 4.0);
        grid.possibilities(Vec3::new(7.241, 0.0, 0.0));
    }

    #[test]
    fn builds_possibility_grid() {
        let mut atom = Atom::default();
        atom.coords = Vec3::new(0.5, 0.5, 0.5);
        atom.base.atom_type = AtomType {
            el: EL_TYPE_C,
            ad: AD_TYPE_C,
            xs: XS_TYPE_C_H,
            sy: 0,
        };
        let dims = [
            GridDim::new(0.0, 3.0, 1),
            GridDim::new(0.0, 3.0, 1),
            GridDim::new(0.0, 3.0, 1),
        ];
        let grid = SzvGrid::new(&FakeModel { atoms: vec![atom] }, &dims, 4.0);
        assert_eq!(grid.possibilities(Vec3::new(0.1, 0.1, 0.1)), &[0]);
        assert_eq!(szv_grid_dims(&dims)[0].n_voxels, 1);
    }
}
