use crate::common::{
    normalize_angle, normalized_angle, sum, vec_distance_sqr, Fl, Vec3, PI, ZERO_VEC,
};
use crate::quaternion::{
    quaternion_difference, quaternion_increment, quaternion_to_r3, random_orientation, Quaternion,
    QT_IDENTITY,
};
use crate::random::{random_fl, random_in_box, random_inside_sphere, Rng64};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scale {
    pub position: Fl,
    pub orientation: Fl,
    pub torsion: Fl,
}

impl Scale {
    pub const fn new(position: Fl, orientation: Fl, torsion: Fl) -> Self {
        Self {
            position,
            orientation,
            torsion,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ConfSize {
    pub ligands: Vec<usize>,
    pub flex: Vec<usize>,
}

impl ConfSize {
    pub fn num_degrees_of_freedom(&self) -> usize {
        sum(&self.ligands) + sum(&self.flex) + 6 * self.ligands.len()
    }
}

pub fn torsions_set_to_null(torsions: &mut [Fl]) {
    for torsion in torsions {
        *torsion = 0.0;
    }
}

pub fn torsions_increment(torsions: &mut [Fl], change: &[Fl], factor: Fl) {
    assert_eq!(torsions.len(), change.len());
    for (torsion, delta) in torsions.iter_mut().zip(change) {
        *torsion += normalized_angle(factor * *delta);
        normalize_angle(torsion);
    }
}

pub fn torsions_randomize(torsions: &mut [Fl], generator: &mut Rng64) {
    for torsion in torsions {
        *torsion = random_fl(-PI, PI, generator);
    }
}

pub fn torsions_too_close(torsions1: &[Fl], torsions2: &[Fl], cutoff: Fl) -> bool {
    assert_eq!(torsions1.len(), torsions2.len());
    torsions1
        .iter()
        .zip(torsions2)
        .all(|(a, b)| normalized_angle(*a - *b).abs() <= cutoff)
}

pub fn torsions_generate(
    torsions: &mut [Fl],
    spread: Fl,
    rp: Fl,
    rs: Option<&[Fl]>,
    generator: &mut Rng64,
) {
    if let Some(rs) = rs {
        assert_eq!(rs.len(), torsions.len());
    }
    for (index, torsion) in torsions.iter_mut().enumerate() {
        if let Some(rs) = rs {
            if random_fl(0.0, 1.0, generator) < rp {
                *torsion = rs[index];
                continue;
            }
        }
        *torsion += random_fl(-spread, spread, generator);
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidChange {
    pub position: Vec3,
    pub orientation: Vec3,
}

impl Default for RigidChange {
    fn default() -> Self {
        Self {
            position: ZERO_VEC,
            orientation: ZERO_VEC,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidConf {
    pub position: Vec3,
    pub orientation: Quaternion,
}

impl Default for RigidConf {
    fn default() -> Self {
        Self {
            position: ZERO_VEC,
            orientation: QT_IDENTITY,
        }
    }
}

impl RigidConf {
    pub fn set_to_null(&mut self) {
        self.position = ZERO_VEC;
        self.orientation = QT_IDENTITY;
    }

    pub fn increment(&mut self, change: &RigidChange, factor: Fl) {
        self.position += factor * change.position;
        quaternion_increment(&mut self.orientation, factor * change.orientation);
    }

    pub fn randomize(&mut self, corner1: Vec3, corner2: Vec3, generator: &mut Rng64) {
        self.position = random_in_box(corner1, corner2, generator);
        self.orientation = random_orientation(generator);
    }

    pub fn too_close(&self, other: &Self, position_cutoff: Fl, orientation_cutoff: Fl) -> bool {
        vec_distance_sqr(self.position, other.position) <= position_cutoff * position_cutoff
            && quaternion_difference(self.orientation, other.orientation).norm_sqr()
                <= orientation_cutoff * orientation_cutoff
    }

    pub fn mutate_position(&mut self, spread: Fl, generator: &mut Rng64) {
        self.position += spread * random_inside_sphere(generator);
    }

    pub fn mutate_orientation(&mut self, spread: Fl, generator: &mut Rng64) {
        quaternion_increment(
            &mut self.orientation,
            spread * random_inside_sphere(generator),
        );
    }

    pub fn generate(
        &mut self,
        position_spread: Fl,
        orientation_spread: Fl,
        rp: Fl,
        rs: Option<&RigidConf>,
        generator: &mut Rng64,
    ) {
        if let Some(rs) = rs {
            if random_fl(0.0, 1.0, generator) < rp {
                self.position = rs.position;
            } else {
                self.mutate_position(position_spread, generator);
            }
            if random_fl(0.0, 1.0, generator) < rp {
                self.orientation = rs.orientation;
            } else {
                self.mutate_orientation(orientation_spread, generator);
            }
        } else {
            self.mutate_position(position_spread, generator);
            self.mutate_orientation(orientation_spread, generator);
        }
    }

    pub fn apply(&self, input: &[Vec3], output: &mut [Vec3], begin: usize, end: usize) {
        assert_eq!(input.len(), output.len());
        let matrix = quaternion_to_r3(self.orientation);
        for i in begin..end {
            output[i] = matrix.mul_vec(input[i]) + self.position;
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct LigandChange {
    pub rigid: RigidChange,
    pub torsions: Vec<Fl>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct LigandConf {
    pub rigid: RigidConf,
    pub torsions: Vec<Fl>,
}

impl LigandConf {
    pub fn set_to_null(&mut self) {
        self.rigid.set_to_null();
        torsions_set_to_null(&mut self.torsions);
    }

    pub fn increment(&mut self, change: &LigandChange, factor: Fl) {
        self.rigid.increment(&change.rigid, factor);
        torsions_increment(&mut self.torsions, &change.torsions, factor);
    }

    pub fn randomize(&mut self, corner1: Vec3, corner2: Vec3, generator: &mut Rng64) {
        self.rigid.randomize(corner1, corner2, generator);
        torsions_randomize(&mut self.torsions, generator);
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ResidueChange {
    pub torsions: Vec<Fl>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ResidueConf {
    pub torsions: Vec<Fl>,
}

impl ResidueConf {
    pub fn set_to_null(&mut self) {
        torsions_set_to_null(&mut self.torsions);
    }

    pub fn increment(&mut self, change: &ResidueChange, factor: Fl) {
        torsions_increment(&mut self.torsions, &change.torsions, factor);
    }

    pub fn randomize(&mut self, generator: &mut Rng64) {
        torsions_randomize(&mut self.torsions, generator);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Change {
    pub ligands: Vec<LigandChange>,
    pub flex: Vec<ResidueChange>,
}

impl Change {
    pub fn new(size: &ConfSize) -> Self {
        Self {
            ligands: size
                .ligands
                .iter()
                .map(|count| LigandChange {
                    rigid: RigidChange::default(),
                    torsions: vec![0.0; *count],
                })
                .collect(),
            flex: size
                .flex
                .iter()
                .map(|count| ResidueChange {
                    torsions: vec![0.0; *count],
                })
                .collect(),
        }
    }

    pub fn get(&self, mut index: usize) -> Fl {
        for ligand in &self.ligands {
            if index < 3 {
                return ligand.rigid.position[index];
            }
            index -= 3;
            if index < 3 {
                return ligand.rigid.orientation[index];
            }
            index -= 3;
            if index < ligand.torsions.len() {
                return ligand.torsions[index];
            }
            index -= ligand.torsions.len();
        }
        for residue in &self.flex {
            if index < residue.torsions.len() {
                return residue.torsions[index];
            }
            index -= residue.torsions.len();
        }
        panic!("change index out of range");
    }

    pub fn get_mut(&mut self, mut index: usize) -> &mut Fl {
        for ligand in &mut self.ligands {
            if index < 3 {
                return &mut ligand.rigid.position[index];
            }
            index -= 3;
            if index < 3 {
                return &mut ligand.rigid.orientation[index];
            }
            index -= 3;
            if index < ligand.torsions.len() {
                return &mut ligand.torsions[index];
            }
            index -= ligand.torsions.len();
        }
        for residue in &mut self.flex {
            if index < residue.torsions.len() {
                return &mut residue.torsions[index];
            }
            index -= residue.torsions.len();
        }
        panic!("change index out of range");
    }

    pub fn num_floats(&self) -> usize {
        self.ligands
            .iter()
            .map(|ligand| 6 + ligand.torsions.len())
            .sum::<usize>()
            + self
                .flex
                .iter()
                .map(|residue| residue.torsions.len())
                .sum::<usize>()
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Conf {
    pub ligands: Vec<LigandConf>,
    pub flex: Vec<ResidueConf>,
}

impl Conf {
    pub fn new(size: &ConfSize) -> Self {
        Self {
            ligands: size
                .ligands
                .iter()
                .map(|count| LigandConf {
                    rigid: RigidConf::default(),
                    torsions: vec![0.0; *count],
                })
                .collect(),
            flex: size
                .flex
                .iter()
                .map(|count| ResidueConf {
                    torsions: vec![0.0; *count],
                })
                .collect(),
        }
    }

    pub fn set_to_null(&mut self) {
        for ligand in &mut self.ligands {
            ligand.set_to_null();
        }
        for residue in &mut self.flex {
            residue.set_to_null();
        }
    }

    pub fn increment(&mut self, change: &Change, factor: Fl) {
        for (ligand, change) in self.ligands.iter_mut().zip(&change.ligands) {
            ligand.increment(change, factor);
        }
        for (residue, change) in self.flex.iter_mut().zip(&change.flex) {
            residue.increment(change, factor);
        }
    }

    pub fn internal_too_close(&self, other: &Conf, torsions_cutoff: Fl) -> bool {
        assert_eq!(self.ligands.len(), other.ligands.len());
        self.ligands
            .iter()
            .zip(&other.ligands)
            .all(|(a, b)| torsions_too_close(&a.torsions, &b.torsions, torsions_cutoff))
    }

    pub fn external_too_close(&self, other: &Conf, cutoff: Scale) -> bool {
        assert_eq!(self.ligands.len(), other.ligands.len());
        assert_eq!(self.flex.len(), other.flex.len());
        self.ligands.iter().zip(&other.ligands).all(|(a, b)| {
            a.rigid
                .too_close(&b.rigid, cutoff.position, cutoff.orientation)
        }) && self
            .flex
            .iter()
            .zip(&other.flex)
            .all(|(a, b)| torsions_too_close(&a.torsions, &b.torsions, cutoff.torsion))
    }

    pub fn too_close(&self, other: &Conf, cutoff: Scale) -> bool {
        self.internal_too_close(other, cutoff.torsion) && self.external_too_close(other, cutoff)
    }

    pub fn generate_internal(
        &mut self,
        torsion_spread: Fl,
        rp: Fl,
        rs: Option<&Conf>,
        generator: &mut Rng64,
    ) {
        for (index, ligand) in self.ligands.iter_mut().enumerate() {
            ligand.rigid.position = ZERO_VEC;
            ligand.rigid.orientation = QT_IDENTITY;
            let torsions_rs = rs.map(|conf| conf.ligands[index].torsions.as_slice());
            torsions_generate(
                &mut ligand.torsions,
                torsion_spread,
                rp,
                torsions_rs,
                generator,
            );
        }
    }

    pub fn generate_external(
        &mut self,
        spread: Scale,
        rp: Fl,
        rs: Option<&Conf>,
        generator: &mut Rng64,
    ) {
        for (index, ligand) in self.ligands.iter_mut().enumerate() {
            let rigid_rs = rs.map(|conf| &conf.ligands[index].rigid);
            ligand
                .rigid
                .generate(spread.position, spread.orientation, rp, rigid_rs, generator);
        }
        for (index, residue) in self.flex.iter_mut().enumerate() {
            let torsions_rs = rs.map(|conf| conf.flex[index].torsions.as_slice());
            torsions_generate(
                &mut residue.torsions,
                spread.torsion,
                rp,
                torsions_rs,
                generator,
            );
        }
    }

    pub fn randomize(&mut self, corner1: Vec3, corner2: Vec3, generator: &mut Rng64) {
        for ligand in &mut self.ligands {
            ligand.randomize(corner1, corner2, generator);
        }
        for residue in &mut self.flex {
            residue.randomize(generator);
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputType {
    pub c: Conf,
    pub e: Fl,
    pub lb: Fl,
    pub ub: Fl,
    pub intra: Fl,
    pub inter: Fl,
    pub conf_independent: Fl,
    pub unbound: Fl,
    pub total: Fl,
    pub coords: Vec<Vec3>,
}

impl OutputType {
    pub fn new(c: Conf, e: Fl) -> Self {
        Self {
            c,
            e,
            lb: 0.0,
            ub: 0.0,
            intra: 0.0,
            inter: 0.0,
            conf_independent: 0.0,
            unbound: 0.0,
            total: 0.0,
            coords: Vec::new(),
        }
    }
}

pub type OutputContainer = Vec<OutputType>;

pub fn output_compare(a: &OutputType, b: &OutputType) -> std::cmp::Ordering {
    a.e.partial_cmp(&b.e).unwrap_or(std::cmp::Ordering::Equal)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_fl;

    #[test]
    fn conf_size_counts_degrees_of_freedom() {
        let size = ConfSize {
            ligands: vec![2, 3],
            flex: vec![4],
        };
        assert_eq!(size.num_degrees_of_freedom(), 21);
    }

    #[test]
    fn change_flat_indexing_matches_reference_order() {
        let size = ConfSize {
            ligands: vec![2],
            flex: vec![1],
        };
        let mut change = Change::new(&size);
        *change.get_mut(0) = 1.0;
        *change.get_mut(3) = 2.0;
        *change.get_mut(6) = 3.0;
        *change.get_mut(8) = 4.0;
        assert_eq!(change.num_floats(), 9);
        assert_eq!(change.get(0), 1.0);
        assert_eq!(change.get(3), 2.0);
        assert_eq!(change.get(6), 3.0);
        assert_eq!(change.get(8), 4.0);
    }

    #[test]
    fn torsion_increment_normalizes() {
        let mut torsions = vec![PI - 0.1];
        torsions_increment(&mut torsions, &[1.0], 1.0);
        assert!(torsions[0] < 0.0);
        assert!(eq_fl(torsions[0], -PI + 0.9));
    }
}
