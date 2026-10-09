use crate::common::{Fl, EPSILON_FL, PI};
use crate::conf::Conf;
use crate::quaternion::quaternion_increment;
use crate::random::{random_fl, random_inside_sphere, random_int, Rng64};

pub trait GyrationRadiusModel {
    fn gyration_radius(&self, ligand_index: usize) -> Fl;
}

pub fn count_mutable_entities(conf: &Conf) -> usize {
    conf.ligands
        .iter()
        .map(|ligand| 2 + ligand.torsions.len())
        .sum::<usize>()
        + conf
            .flex
            .iter()
            .map(|residue| residue.torsions.len())
            .sum::<usize>()
}

pub fn mutate_conf<M: GyrationRadiusModel>(
    conf: &mut Conf,
    model: &M,
    amplitude: Fl,
    generator: &mut Rng64,
) {
    let mutable_entities_num = count_mutable_entities(conf);
    if mutable_entities_num == 0 {
        return;
    }
    let mut which = random_int(0, mutable_entities_num as i32 - 1, generator) as usize;

    for (i, ligand) in conf.ligands.iter_mut().enumerate() {
        if which == 0 {
            ligand.rigid.position += amplitude * random_inside_sphere(generator);
            return;
        }
        which -= 1;
        if which == 0 {
            let gr = model.gyration_radius(i);
            if gr > EPSILON_FL {
                let rotation = (amplitude / gr) * random_inside_sphere(generator);
                quaternion_increment(&mut ligand.rigid.orientation, rotation);
            }
            return;
        }
        which -= 1;
        if which < ligand.torsions.len() {
            ligand.torsions[which] = random_fl(-PI, PI, generator);
            return;
        }
        which -= ligand.torsions.len();
    }

    for residue in &mut conf.flex {
        if which < residue.torsions.len() {
            residue.torsions[which] = random_fl(-PI, PI, generator);
            return;
        }
        which -= residue.torsions.len();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conf::{Conf, ConfSize};

    struct FakeModel;

    impl GyrationRadiusModel for FakeModel {
        fn gyration_radius(&self, _ligand_index: usize) -> Fl {
            1.0
        }
    }

    #[test]
    fn counts_mutable_entities_like_reference() {
        let conf = Conf::new(&ConfSize {
            ligands: vec![2],
            flex: vec![3],
        });
        assert_eq!(count_mutable_entities(&conf), 7);
    }

    #[test]
    fn mutation_handles_empty_conf() {
        let mut conf = Conf::default();
        mutate_conf(&mut conf, &FakeModel, 2.0, &mut Rng64::new(1));
        assert!(conf.ligands.is_empty());
    }
}
