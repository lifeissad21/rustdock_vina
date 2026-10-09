use crate::common::{Fl, Sz};

pub const EL_TYPE_H: Sz = 0;
pub const EL_TYPE_C: Sz = 1;
pub const EL_TYPE_N: Sz = 2;
pub const EL_TYPE_O: Sz = 3;
pub const EL_TYPE_S: Sz = 4;
pub const EL_TYPE_P: Sz = 5;
pub const EL_TYPE_F: Sz = 6;
pub const EL_TYPE_CL: Sz = 7;
pub const EL_TYPE_BR: Sz = 8;
pub const EL_TYPE_I: Sz = 9;
pub const EL_TYPE_SI: Sz = 10;
pub const EL_TYPE_AT: Sz = 11;
pub const EL_TYPE_MET: Sz = 12;
pub const EL_TYPE_DUMMY: Sz = 13;
pub const EL_TYPE_SIZE: Sz = 14;

pub const AD_TYPE_C: Sz = 0;
pub const AD_TYPE_A: Sz = 1;
pub const AD_TYPE_N: Sz = 2;
pub const AD_TYPE_O: Sz = 3;
pub const AD_TYPE_P: Sz = 4;
pub const AD_TYPE_S: Sz = 5;
pub const AD_TYPE_H: Sz = 6;
pub const AD_TYPE_F: Sz = 7;
pub const AD_TYPE_I: Sz = 8;
pub const AD_TYPE_NA: Sz = 9;
pub const AD_TYPE_OA: Sz = 10;
pub const AD_TYPE_SA: Sz = 11;
pub const AD_TYPE_HD: Sz = 12;
pub const AD_TYPE_MG: Sz = 13;
pub const AD_TYPE_MN: Sz = 14;
pub const AD_TYPE_ZN: Sz = 15;
pub const AD_TYPE_CA: Sz = 16;
pub const AD_TYPE_FE: Sz = 17;
pub const AD_TYPE_CL: Sz = 18;
pub const AD_TYPE_BR: Sz = 19;
pub const AD_TYPE_SI: Sz = 20;
pub const AD_TYPE_AT: Sz = 21;
pub const AD_TYPE_G0: Sz = 22;
pub const AD_TYPE_G1: Sz = 23;
pub const AD_TYPE_G2: Sz = 24;
pub const AD_TYPE_G3: Sz = 25;
pub const AD_TYPE_CG0: Sz = 26;
pub const AD_TYPE_CG1: Sz = 27;
pub const AD_TYPE_CG2: Sz = 28;
pub const AD_TYPE_CG3: Sz = 29;
pub const AD_TYPE_W: Sz = 30;
pub const AD_TYPE_SIZE: Sz = 31;

pub const XS_TYPE_C_H: Sz = 0;
pub const XS_TYPE_C_P: Sz = 1;
pub const XS_TYPE_N_P: Sz = 2;
pub const XS_TYPE_N_D: Sz = 3;
pub const XS_TYPE_N_A: Sz = 4;
pub const XS_TYPE_N_DA: Sz = 5;
pub const XS_TYPE_O_P: Sz = 6;
pub const XS_TYPE_O_D: Sz = 7;
pub const XS_TYPE_O_A: Sz = 8;
pub const XS_TYPE_O_DA: Sz = 9;
pub const XS_TYPE_S_P: Sz = 10;
pub const XS_TYPE_P_P: Sz = 11;
pub const XS_TYPE_F_H: Sz = 12;
pub const XS_TYPE_CL_H: Sz = 13;
pub const XS_TYPE_BR_H: Sz = 14;
pub const XS_TYPE_I_H: Sz = 15;
pub const XS_TYPE_SI: Sz = 16;
pub const XS_TYPE_AT: Sz = 17;
pub const XS_TYPE_MET_D: Sz = 18;
pub const XS_TYPE_C_H_CG0: Sz = 19;
pub const XS_TYPE_C_P_CG0: Sz = 20;
pub const XS_TYPE_G0: Sz = 21;
pub const XS_TYPE_C_H_CG1: Sz = 22;
pub const XS_TYPE_C_P_CG1: Sz = 23;
pub const XS_TYPE_G1: Sz = 24;
pub const XS_TYPE_C_H_CG2: Sz = 25;
pub const XS_TYPE_C_P_CG2: Sz = 26;
pub const XS_TYPE_G2: Sz = 27;
pub const XS_TYPE_C_H_CG3: Sz = 28;
pub const XS_TYPE_C_P_CG3: Sz = 29;
pub const XS_TYPE_G3: Sz = 30;
pub const XS_TYPE_W: Sz = 31;
pub const XS_TYPE_SIZE: Sz = 32;

pub const SY_TYPE_C_3: Sz = 0;
pub const SY_TYPE_C_2: Sz = 1;
pub const SY_TYPE_C_AR: Sz = 2;
pub const SY_TYPE_C_CAT: Sz = 3;
pub const SY_TYPE_N_3: Sz = 4;
pub const SY_TYPE_N_AR: Sz = 5;
pub const SY_TYPE_N_AM: Sz = 6;
pub const SY_TYPE_N_PL3: Sz = 7;
pub const SY_TYPE_O_3: Sz = 8;
pub const SY_TYPE_O_2: Sz = 9;
pub const SY_TYPE_O_CO2: Sz = 10;
pub const SY_TYPE_S: Sz = 11;
pub const SY_TYPE_P: Sz = 12;
pub const SY_TYPE_F: Sz = 13;
pub const SY_TYPE_CL: Sz = 14;
pub const SY_TYPE_BR: Sz = 15;
pub const SY_TYPE_I: Sz = 16;
pub const SY_TYPE_MET: Sz = 17;
pub const SY_TYPE_SIZE: Sz = 18;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomKind {
    pub name: &'static str,
    pub radius: Fl,
    pub depth: Fl,
    pub hb_depth: Fl,
    pub hb_radius: Fl,
    pub solvation: Fl,
    pub volume: Fl,
    pub covalent_radius: Fl,
}

pub const ATOM_KIND_DATA: [AtomKind; AD_TYPE_SIZE] = [
    AtomKind {
        name: "C",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "A",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00052,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "N",
        radius: 1.75000,
        depth: 0.16000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00162,
        volume: 22.44930,
        covalent_radius: 0.75,
    },
    AtomKind {
        name: "O",
        radius: 1.60000,
        depth: 0.20000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00251,
        volume: 17.15730,
        covalent_radius: 0.73,
    },
    AtomKind {
        name: "P",
        radius: 2.10000,
        depth: 0.20000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 38.79240,
        covalent_radius: 1.06,
    },
    AtomKind {
        name: "S",
        radius: 2.00000,
        depth: 0.20000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00214,
        volume: 33.51030,
        covalent_radius: 1.02,
    },
    AtomKind {
        name: "H",
        radius: 1.00000,
        depth: 0.02000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00051,
        volume: 0.00000,
        covalent_radius: 0.37,
    },
    AtomKind {
        name: "F",
        radius: 1.54500,
        depth: 0.08000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 15.44800,
        covalent_radius: 0.71,
    },
    AtomKind {
        name: "I",
        radius: 2.36000,
        depth: 0.55000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 55.05850,
        covalent_radius: 1.33,
    },
    AtomKind {
        name: "NA",
        radius: 1.75000,
        depth: 0.16000,
        hb_depth: -5.0,
        hb_radius: 1.9,
        solvation: -0.00162,
        volume: 22.44930,
        covalent_radius: 0.75,
    },
    AtomKind {
        name: "OA",
        radius: 1.60000,
        depth: 0.20000,
        hb_depth: -5.0,
        hb_radius: 1.9,
        solvation: -0.00251,
        volume: 17.15730,
        covalent_radius: 0.73,
    },
    AtomKind {
        name: "SA",
        radius: 2.00000,
        depth: 0.20000,
        hb_depth: -1.0,
        hb_radius: 2.5,
        solvation: -0.00214,
        volume: 33.51030,
        covalent_radius: 1.02,
    },
    AtomKind {
        name: "HD",
        radius: 1.00000,
        depth: 0.02000,
        hb_depth: 1.0,
        hb_radius: 0.0,
        solvation: 0.00051,
        volume: 0.00000,
        covalent_radius: 0.37,
    },
    AtomKind {
        name: "Mg",
        radius: 0.65000,
        depth: 0.87500,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 1.56000,
        covalent_radius: 1.30,
    },
    AtomKind {
        name: "Mn",
        radius: 0.65000,
        depth: 0.87500,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 2.14000,
        covalent_radius: 1.39,
    },
    AtomKind {
        name: "Zn",
        radius: 0.74000,
        depth: 0.55000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 1.70000,
        covalent_radius: 1.31,
    },
    AtomKind {
        name: "Ca",
        radius: 0.99000,
        depth: 0.55000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 2.77000,
        covalent_radius: 1.74,
    },
    AtomKind {
        name: "Fe",
        radius: 0.65000,
        depth: 0.01000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 1.84000,
        covalent_radius: 1.25,
    },
    AtomKind {
        name: "Cl",
        radius: 2.04500,
        depth: 0.27600,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 35.82350,
        covalent_radius: 0.99,
    },
    AtomKind {
        name: "Br",
        radius: 2.16500,
        depth: 0.38900,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 42.56610,
        covalent_radius: 1.14,
    },
    AtomKind {
        name: "Si",
        radius: 2.30000,
        depth: 0.20000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 50.96500,
        covalent_radius: 1.11,
    },
    AtomKind {
        name: "At",
        radius: 2.40000,
        depth: 0.55000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00110,
        volume: 57.90580,
        covalent_radius: 1.44,
    },
    AtomKind {
        name: "G0",
        radius: 0.00000,
        depth: 0.00000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00000,
        volume: 0.00000,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "G1",
        radius: 0.00000,
        depth: 0.00000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00000,
        volume: 0.00000,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "G2",
        radius: 0.00000,
        depth: 0.00000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00000,
        volume: 0.00000,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "G3",
        radius: 0.00000,
        depth: 0.00000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00000,
        volume: 0.00000,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "CG0",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "CG1",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "CG2",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "CG3",
        radius: 2.00000,
        depth: 0.15000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: -0.00143,
        volume: 33.51030,
        covalent_radius: 0.77,
    },
    AtomKind {
        name: "W",
        radius: 0.00000,
        depth: 0.00000,
        hb_depth: 0.0,
        hb_radius: 0.0,
        solvation: 0.00000,
        volume: 0.00000,
        covalent_radius: 0.00,
    },
];

pub const METAL_SOLVATION_PARAMETER: Fl = -0.00110;
pub const METAL_COVALENT_RADIUS: Fl = 1.75;
pub const ATOM_EQUIVALENCE_DATA: [(&str, &str); 1] = [("Se", "S")];
pub const ACCEPTOR_KIND_DATA: [(Sz, Fl, Fl); 3] = [
    (AD_TYPE_NA, 1.9, 5.0),
    (AD_TYPE_OA, 1.9, 5.0),
    (AD_TYPE_SA, 2.5, 1.0),
];
pub const XS_VDW_RADII: [Fl; XS_TYPE_SIZE] = [
    1.9, 1.9, 1.8, 1.8, 1.8, 1.8, 1.7, 1.7, 1.7, 1.7, 2.0, 2.1, 1.5, 1.8, 2.0, 2.2, 2.2, 2.3, 1.2,
    1.9, 1.9, 0.0, 1.9, 1.9, 0.0, 1.9, 1.9, 0.0, 1.9, 1.9, 0.0, 0.0,
];
pub const XS_VINARDO_VDW_RADII: [Fl; XS_TYPE_SIZE] = [
    2.0, 2.0, 1.7, 1.7, 1.7, 1.7, 1.6, 1.6, 1.6, 1.6, 2.0, 2.1, 1.5, 1.8, 2.0, 2.2, 2.2, 2.3, 1.2,
    2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0,
];
pub const NON_AD_METAL_NAMES: [&str; 9] = ["Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"];

pub fn ad_is_hydrogen(ad: Sz) -> bool {
    ad == AD_TYPE_H || ad == AD_TYPE_HD
}

pub fn ad_is_heteroatom(ad: Sz) -> bool {
    ad != AD_TYPE_A && ad != AD_TYPE_C && ad != AD_TYPE_H && ad != AD_TYPE_HD && ad < AD_TYPE_SIZE
}

pub fn ad_type_to_el_type(t: Sz) -> Sz {
    match t {
        AD_TYPE_C | AD_TYPE_A | AD_TYPE_CG0 | AD_TYPE_CG1 | AD_TYPE_CG2 | AD_TYPE_CG3 => EL_TYPE_C,
        AD_TYPE_N | AD_TYPE_NA => EL_TYPE_N,
        AD_TYPE_O | AD_TYPE_OA => EL_TYPE_O,
        AD_TYPE_P => EL_TYPE_P,
        AD_TYPE_S | AD_TYPE_SA => EL_TYPE_S,
        AD_TYPE_H | AD_TYPE_HD => EL_TYPE_H,
        AD_TYPE_F => EL_TYPE_F,
        AD_TYPE_I => EL_TYPE_I,
        AD_TYPE_MG | AD_TYPE_MN | AD_TYPE_ZN | AD_TYPE_CA | AD_TYPE_FE => EL_TYPE_MET,
        AD_TYPE_CL => EL_TYPE_CL,
        AD_TYPE_BR => EL_TYPE_BR,
        AD_TYPE_SI => EL_TYPE_SI,
        AD_TYPE_AT => EL_TYPE_AT,
        AD_TYPE_G0 | AD_TYPE_G1 | AD_TYPE_G2 | AD_TYPE_G3 | AD_TYPE_W => EL_TYPE_DUMMY,
        AD_TYPE_SIZE => EL_TYPE_SIZE,
        _ => panic!("unknown AutoDock atom type {t}"),
    }
}

pub fn xs_radius(t: Sz) -> Fl {
    XS_VDW_RADII[t]
}

pub fn xs_vinardo_radius(t: Sz) -> Fl {
    XS_VINARDO_VDW_RADII[t]
}

pub fn is_non_ad_metal_name(name: &str) -> bool {
    NON_AD_METAL_NAMES.contains(&name)
}

pub fn xs_is_hydrophobic(xs: Sz) -> bool {
    matches!(
        xs,
        XS_TYPE_C_H | XS_TYPE_F_H | XS_TYPE_CL_H | XS_TYPE_BR_H | XS_TYPE_I_H
    )
}

pub fn xs_is_acceptor(xs: Sz) -> bool {
    matches!(xs, XS_TYPE_N_A | XS_TYPE_N_DA | XS_TYPE_O_A | XS_TYPE_O_DA)
}

pub fn xs_is_donor(xs: Sz) -> bool {
    matches!(
        xs,
        XS_TYPE_N_D | XS_TYPE_N_DA | XS_TYPE_O_D | XS_TYPE_O_DA | XS_TYPE_MET_D
    )
}

pub fn xs_donor_acceptor(t1: Sz, t2: Sz) -> bool {
    xs_is_donor(t1) && xs_is_acceptor(t2)
}

pub fn xs_h_bond_possible(t1: Sz, t2: Sz) -> bool {
    xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1)
}

pub fn ad_type_property(i: Sz) -> &'static AtomKind {
    &ATOM_KIND_DATA[i]
}

pub fn string_to_ad_type(name: &str) -> Sz {
    if let Some((index, _)) = ATOM_KIND_DATA
        .iter()
        .enumerate()
        .find(|(_, kind)| kind.name == name)
    {
        return index;
    }
    if let Some((_, to)) = ATOM_EQUIVALENCE_DATA.iter().find(|(from, _)| *from == name) {
        return string_to_ad_type(to);
    }
    AD_TYPE_SIZE
}

pub fn max_covalent_radius() -> Fl {
    ATOM_KIND_DATA
        .iter()
        .map(|kind| kind.covalent_radius)
        .fold(0.0, Fl::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_fl;

    #[test]
    fn maps_ad_types_to_elements_like_reference() {
        assert_eq!(ad_type_to_el_type(AD_TYPE_C), EL_TYPE_C);
        assert_eq!(ad_type_to_el_type(AD_TYPE_NA), EL_TYPE_N);
        assert_eq!(ad_type_to_el_type(AD_TYPE_G0), EL_TYPE_DUMMY);
        assert_eq!(ad_type_to_el_type(AD_TYPE_SIZE), EL_TYPE_SIZE);
    }

    #[test]
    fn resolves_names_and_equivalences() {
        assert_eq!(string_to_ad_type("OA"), AD_TYPE_OA);
        assert_eq!(string_to_ad_type("Se"), AD_TYPE_S);
        assert_eq!(string_to_ad_type("Unknown"), AD_TYPE_SIZE);
    }

    #[test]
    fn exposes_scoring_class_helpers() {
        assert!(ad_is_hydrogen(AD_TYPE_HD));
        assert!(ad_is_heteroatom(AD_TYPE_ZN));
        assert!(xs_h_bond_possible(XS_TYPE_N_D, XS_TYPE_O_A));
        assert!(is_non_ad_metal_name("Cu"));
        assert!(eq_fl(max_covalent_radius(), 1.74));
    }
}
