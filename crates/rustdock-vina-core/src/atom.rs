use crate::atom_base::AtomBase;
use crate::common::{Fl, Vec3, MAX_SZ, MAX_VEC};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomIndex {
    pub i: usize,
    pub in_grid: bool,
}

impl Default for AtomIndex {
    fn default() -> Self {
        Self {
            i: MAX_SZ,
            in_grid: false,
        }
    }
}

impl AtomIndex {
    pub const fn new(i: usize, in_grid: bool) -> Self {
        Self { i, in_grid }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub connected_atom_index: AtomIndex,
    pub length: Fl,
    pub rotatable: bool,
}

impl Default for Bond {
    fn default() -> Self {
        Self {
            connected_atom_index: AtomIndex::default(),
            length: 0.0,
            rotatable: false,
        }
    }
}

impl Bond {
    pub const fn new(connected_atom_index: AtomIndex, length: Fl, rotatable: bool) -> Self {
        Self {
            connected_atom_index,
            length,
            rotatable,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub base: AtomBase,
    pub coords: Vec3,
    pub bonds: Vec<Bond>,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            base: AtomBase::default(),
            coords: MAX_VEC,
            bonds: Vec::new(),
        }
    }
}

pub type AtomVec = Vec<Atom>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_reference_constructors() {
        assert_eq!(AtomIndex::default().i, MAX_SZ);
        assert!(!AtomIndex::default().in_grid);
        assert_eq!(Bond::default().length, 0.0);
        assert!(!Bond::default().rotatable);
        assert_eq!(Atom::default().coords, MAX_VEC);
    }
}
