use crate::triangular_matrix_index::{
    strictly_triangular_matrix_index, strictly_triangular_matrix_index_permissive,
    triangular_matrix_index, triangular_matrix_index_permissive,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<T> {
    data: Vec<T>,
    i: usize,
    j: usize,
}

impl<T> Default for Matrix<T> {
    fn default() -> Self {
        Self {
            data: Vec::new(),
            i: 0,
            j: 0,
        }
    }
}

impl<T: Clone> Matrix<T> {
    pub fn new(i: usize, j: usize, filler: T) -> Self {
        Self {
            data: vec![filler; i * j],
            i,
            j,
        }
    }

    pub fn resize(&mut self, m: usize, n: usize, filler: T) {
        if m == self.dim_1() && n == self.dim_2() {
            return;
        }
        assert!(m >= self.dim_1());
        assert!(n >= self.dim_2());
        let mut tmp = vec![filler; m * n];
        for i in 0..self.i {
            for j in 0..self.j {
                tmp[i + m * j] = self.get(i, j).clone();
            }
        }
        self.data = tmp;
        self.i = m;
        self.j = n;
    }

    pub fn append(&mut self, x: &Matrix<T>, filler: T) {
        let m = self.dim_1();
        let n = self.dim_2();
        self.resize(m + x.dim_1(), n + x.dim_2(), filler);
        for i in 0..x.dim_1() {
            for j in 0..x.dim_2() {
                *self.get_mut(i + m, j + n) = x.get(i, j).clone();
            }
        }
    }
}

impl<T> Matrix<T> {
    pub fn index(&self, i: usize, j: usize) -> usize {
        debug_assert!(j < self.j);
        debug_assert!(i < self.i);
        i + self.i * j
    }

    pub fn get_flat(&self, i: usize) -> &T {
        &self.data[i]
    }

    pub fn get_flat_mut(&mut self, i: usize) -> &mut T {
        &mut self.data[i]
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[self.index(i, j)]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        let index = self.index(i, j);
        &mut self.data[index]
    }

    pub fn dim_1(&self) -> usize {
        self.i
    }

    pub fn dim_2(&self) -> usize {
        self.j
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TriangularMatrix<T> {
    data: Vec<T>,
    dim: usize,
}

impl<T: Clone> TriangularMatrix<T> {
    pub fn new(dim: usize, filler: T) -> Self {
        Self {
            data: vec![filler; dim * (dim + 1) / 2],
            dim,
        }
    }
}

impl<T> TriangularMatrix<T> {
    pub fn index(&self, i: usize, j: usize) -> usize {
        triangular_matrix_index(self.dim, i, j)
    }

    pub fn index_permissive(&self, i: usize, j: usize) -> usize {
        triangular_matrix_index_permissive(self.dim, i, j)
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[self.index(i, j)]
    }

    pub fn get_flat(&self, index: usize) -> &T {
        &self.data[index]
    }

    pub fn get_permissive(&self, i: usize, j: usize) -> &T {
        &self.data[self.index_permissive(i, j)]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        let index = self.index(i, j);
        &mut self.data[index]
    }

    pub fn dim(&self) -> usize {
        self.dim
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StrictlyTriangularMatrix<T> {
    data: Vec<T>,
    dim: usize,
}

impl<T: Clone> StrictlyTriangularMatrix<T> {
    pub fn new(dim: usize, filler: T) -> Self {
        Self {
            data: vec![filler; dim * dim.saturating_sub(1) / 2],
            dim,
        }
    }

    pub fn resize(&mut self, dim: usize, filler: T) {
        if dim == self.dim {
            return;
        }
        assert!(dim > self.dim);
        self.dim = dim;
        self.data.resize(dim * dim.saturating_sub(1) / 2, filler);
    }

    pub fn append(&mut self, matrix: &StrictlyTriangularMatrix<T>, filler: T) {
        let n = self.dim();
        self.resize(n + matrix.dim(), filler);
        for i in 0..matrix.dim() {
            for j in (i + 1)..matrix.dim() {
                *self.get_mut(i + n, j + n) = matrix.get(i, j).clone();
            }
        }
    }

    pub fn append_rectangular(
        &mut self,
        rectangular: &Matrix<T>,
        triangular: &StrictlyTriangularMatrix<T>,
    ) {
        assert_eq!(self.dim(), rectangular.dim_1());
        assert_eq!(rectangular.dim_2(), triangular.dim());

        if rectangular.dim_2() == 0 {
            return;
        }
        if rectangular.dim_1() == 0 {
            *self = triangular.clone();
            return;
        }
        let filler = rectangular.get(0, 0).clone();
        let n = self.dim();
        self.append(triangular, filler);
        for i in 0..rectangular.dim_1() {
            for j in 0..rectangular.dim_2() {
                *self.get_mut(i, n + j) = rectangular.get(i, j).clone();
            }
        }
    }
}

impl<T> StrictlyTriangularMatrix<T> {
    pub fn index(&self, i: usize, j: usize) -> usize {
        strictly_triangular_matrix_index(self.dim, i, j)
    }

    pub fn index_permissive(&self, i: usize, j: usize) -> usize {
        strictly_triangular_matrix_index_permissive(self.dim, i, j)
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[self.index(i, j)]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        let index = self.index(i, j);
        &mut self.data[index]
    }

    pub fn dim(&self) -> usize {
        self.dim
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_is_column_major() {
        let mut matrix = Matrix::new(2, 3, 0);
        *matrix.get_mut(1, 2) = 9;
        assert_eq!(matrix.index(1, 2), 5);
        assert_eq!(*matrix.get(1, 2), 9);
    }

    #[test]
    fn triangular_matrix_uses_reference_index() {
        let mut matrix = TriangularMatrix::new(3, 0);
        *matrix.get_mut(1, 2) = 7;
        assert_eq!(matrix.index(1, 2), 4);
        assert_eq!(*matrix.get(1, 2), 7);
    }

    #[test]
    fn strict_append_rectangular_matches_original_layout() {
        let mut base = StrictlyTriangularMatrix::new(1, 0);
        let mut rectangular = Matrix::new(1, 2, 0);
        *rectangular.get_mut(0, 0) = 3;
        *rectangular.get_mut(0, 1) = 4;
        let mut triangular = StrictlyTriangularMatrix::new(2, 0);
        *triangular.get_mut(0, 1) = 5;

        base.append_rectangular(&rectangular, &triangular);

        assert_eq!(base.dim(), 3);
        assert_eq!(*base.get(0, 1), 3);
        assert_eq!(*base.get(0, 2), 4);
        assert_eq!(*base.get(1, 2), 5);
    }
}
