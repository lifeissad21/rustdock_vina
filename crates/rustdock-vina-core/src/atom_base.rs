use crate::atom_type::AtomType;
use crate::common::Fl;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomBase {
    pub atom_type: AtomType,
    pub charge: Fl,
}

impl Default for AtomBase {
    fn default() -> Self {
        Self {
            atom_type: AtomType::default(),
            charge: 0.0,
        }
    }
}
