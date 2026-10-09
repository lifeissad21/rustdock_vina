use crate::atom_constants::*;
use crate::common::{Fl, Sz, MAX_SZ};
use crate::triangular_matrix_index::triangular_matrix_index;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomTyping {
    El,
    Ad,
    Xs,
    Sy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomType {
    pub el: Sz,
    pub ad: Sz,
    pub xs: Sz,
    pub sy: Sz,
}

impl Default for AtomType {
    fn default() -> Self {
        Self {
            el: EL_TYPE_SIZE,
            ad: AD_TYPE_SIZE,
            xs: XS_TYPE_SIZE,
            sy: SY_TYPE_SIZE,
        }
    }
}

impl AtomType {
    pub fn get(self, atom_typing_used: AtomTyping) -> Sz {
        match atom_typing_used {
            AtomTyping::El => self.el,
            AtomTyping::Ad => self.ad,
            AtomTyping::Xs => self.xs,
            AtomTyping::Sy => self.sy,
        }
    }

    pub fn is_hydrogen(self) -> bool {
        ad_is_hydrogen(self.ad)
    }

    pub fn is_heteroatom(self) -> bool {
        ad_is_heteroatom(self.ad) || self.xs == XS_TYPE_MET_D
    }

    pub fn acceptable_type(self) -> bool {
        self.ad < AD_TYPE_SIZE || self.xs == XS_TYPE_MET_D
    }

    pub fn assign_el(&mut self) {
        self.el = ad_type_to_el_type(self.ad);
        if self.ad == AD_TYPE_SIZE && self.xs == XS_TYPE_MET_D {
            self.el = EL_TYPE_MET;
        }
    }

    pub fn same_element(self, other: Self) -> bool {
        self.el == other.el
    }

    pub fn covalent_radius(self) -> Fl {
        if self.ad < AD_TYPE_SIZE {
            ad_type_property(self.ad).covalent_radius
        } else if self.xs == XS_TYPE_MET_D {
            METAL_COVALENT_RADIUS
        } else {
            panic!("atom type has no covalent radius")
        }
    }

    pub fn optimal_covalent_bond_length(self, other: Self) -> Fl {
        self.covalent_radius() + other.covalent_radius()
    }
}

pub fn num_atom_types(atom_typing_used: AtomTyping) -> Sz {
    match atom_typing_used {
        AtomTyping::El => EL_TYPE_SIZE,
        AtomTyping::Ad => AD_TYPE_SIZE,
        AtomTyping::Xs => XS_TYPE_SIZE,
        AtomTyping::Sy => SY_TYPE_SIZE,
    }
}

pub fn get_type_pair_index(atom_typing_used: AtomTyping, a: AtomType, b: AtomType) -> Sz {
    let n = num_atom_types(atom_typing_used);
    let i = a.get(atom_typing_used);
    let j = b.get(atom_typing_used);
    assert!(i < n, "unassigned atom type {}", MAX_SZ);
    assert!(j < n, "unassigned atom type {}", MAX_SZ);

    if i <= j {
        triangular_matrix_index(n, i, j)
    } else {
        triangular_matrix_index(n, j, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_unassigned_reference_state() {
        let atom_type = AtomType::default();
        assert_eq!(atom_type.el, EL_TYPE_SIZE);
        assert_eq!(atom_type.ad, AD_TYPE_SIZE);
        assert_eq!(atom_type.xs, XS_TYPE_SIZE);
        assert_eq!(atom_type.sy, SY_TYPE_SIZE);
    }

    #[test]
    fn assigns_element_and_bond_lengths() {
        let mut carbon = AtomType {
            ad: AD_TYPE_C,
            ..Default::default()
        };
        carbon.assign_el();
        let mut oxygen = AtomType {
            ad: AD_TYPE_OA,
            ..Default::default()
        };
        oxygen.assign_el();

        assert_eq!(carbon.el, EL_TYPE_C);
        assert_eq!(oxygen.el, EL_TYPE_O);
        assert!(oxygen.is_heteroatom());
        assert_eq!(
            get_type_pair_index(AtomTyping::Ad, carbon, oxygen),
            triangular_matrix_index(AD_TYPE_SIZE, AD_TYPE_C, AD_TYPE_OA)
        );
    }
}
