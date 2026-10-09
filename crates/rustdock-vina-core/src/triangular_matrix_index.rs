pub fn triangular_matrix_index(n: usize, i: usize, j: usize) -> usize {
    debug_assert!(j < n);
    debug_assert!(i <= j);
    i + j * (j + 1) / 2
}

pub fn triangular_matrix_index_permissive(n: usize, i: usize, j: usize) -> usize {
    if i <= j {
        triangular_matrix_index(n, i, j)
    } else {
        triangular_matrix_index(n, j, i)
    }
}

pub fn strictly_triangular_matrix_index(n: usize, i: usize, j: usize) -> usize {
    debug_assert!(j < n);
    debug_assert!(i < j);
    debug_assert!(j >= 1);
    i + j * (j - 1) / 2
}

pub fn strictly_triangular_matrix_index_permissive(n: usize, i: usize, j: usize) -> usize {
    if i < j {
        strictly_triangular_matrix_index(n, i, j)
    } else {
        strictly_triangular_matrix_index(n, j, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mirrors_reference_indices() {
        assert_eq!(triangular_matrix_index(4, 0, 0), 0);
        assert_eq!(triangular_matrix_index(4, 0, 1), 1);
        assert_eq!(triangular_matrix_index(4, 1, 1), 2);
        assert_eq!(triangular_matrix_index_permissive(4, 2, 1), 4);
        assert_eq!(strictly_triangular_matrix_index(4, 0, 1), 0);
        assert_eq!(strictly_triangular_matrix_index(4, 1, 2), 2);
    }
}
