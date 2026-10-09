#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CapacityOverflow;

pub fn checked_multiply2(i: usize, j: usize) -> Result<usize, CapacityOverflow> {
    if i == 0 || j == 0 {
        return Ok(0);
    }
    i.checked_mul(j).ok_or(CapacityOverflow)
}

pub fn checked_multiply3(i: usize, j: usize, k: usize) -> Result<usize, CapacityOverflow> {
    checked_multiply2(checked_multiply2(i, j)?, k)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Array3d<T> {
    i: usize,
    j: usize,
    k: usize,
    data: Vec<T>,
}

impl<T> Default for Array3d<T> {
    fn default() -> Self {
        Self {
            i: 0,
            j: 0,
            k: 0,
            data: Vec::new(),
        }
    }
}

impl<T: Clone> Array3d<T> {
    pub fn new(i: usize, j: usize, k: usize, filler: T) -> Result<Self, CapacityOverflow> {
        Ok(Self {
            i,
            j,
            k,
            data: vec![filler; checked_multiply3(i, j, k)?],
        })
    }

    pub fn resize(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        filler: T,
    ) -> Result<(), CapacityOverflow> {
        self.i = i;
        self.j = j;
        self.k = k;
        self.data.resize(checked_multiply3(i, j, k)?, filler);
        Ok(())
    }
}

impl<T> Array3d<T> {
    pub fn dim0(&self) -> usize {
        self.i
    }

    pub fn dim1(&self) -> usize {
        self.j
    }

    pub fn dim2(&self) -> usize {
        self.k
    }

    pub fn dim(&self, index: usize) -> usize {
        match index {
            0 => self.i,
            1 => self.j,
            2 => self.k,
            _ => panic!("array3d dimension index out of range"),
        }
    }

    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        debug_assert!(i < self.i);
        debug_assert!(j < self.j);
        debug_assert!(k < self.k);
        i + self.i * (j + self.j * k)
    }

    pub fn get(&self, i: usize, j: usize, k: usize) -> &T {
        &self.data[self.index(i, j, k)]
    }

    pub fn get_mut(&mut self, i: usize, j: usize, k: usize) -> &mut T {
        let index = self.index(i, j, k);
        &mut self.data[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uses_reference_index_layout() {
        let mut values = Array3d::new(2, 3, 4, 0).unwrap();
        *values.get_mut(1, 2, 3) = 17;
        assert_eq!(values.index(1, 2, 3), 23);
        assert_eq!(*values.get(1, 2, 3), 17);
    }

    #[test]
    fn checked_multiply_detects_overflow() {
        assert!(checked_multiply2(usize::MAX, 2).is_err());
        assert_eq!(checked_multiply3(2, 3, 4).unwrap(), 24);
    }
}
