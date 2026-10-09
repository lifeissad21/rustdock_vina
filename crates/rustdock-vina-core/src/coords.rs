use crate::common::{vec_distance_sqr, Fl, Vec3, MAX_FL};
use crate::conf::{output_compare, OutputContainer, OutputType};

pub fn rmsd_upper_bound(a: &[Vec3], b: &[Vec3]) -> Fl {
    assert_eq!(a.len(), b.len());
    let acc: Fl = a.iter().zip(b).map(|(x, y)| vec_distance_sqr(*x, *y)).sum();
    if a.is_empty() {
        0.0
    } else {
        (acc / a.len() as Fl).sqrt()
    }
}

pub fn find_closest(a: &[Vec3], b: &[OutputType]) -> (usize, Fl) {
    let mut tmp = (b.len(), MAX_FL);
    for (i, output) in b.iter().enumerate() {
        let res = rmsd_upper_bound(a, &output.coords);
        if i == 0 || res < tmp.1 {
            tmp = (i, res);
        }
    }
    tmp
}

pub fn add_to_output_container(
    out: &mut OutputContainer,
    t: OutputType,
    min_rmsd: Fl,
    max_size: usize,
) {
    let closest_rmsd = find_closest(&t.coords, out);
    if closest_rmsd.0 < out.len() && closest_rmsd.1 < min_rmsd {
        if t.e < out[closest_rmsd.0].e {
            out[closest_rmsd.0] = t;
        }
    } else if out.len() < max_size {
        out.push(t);
    } else if !out.is_empty() && t.e < out.last().expect("non-empty checked").e {
        let last = out.len() - 1;
        out[last] = t;
    }
    out.sort_by(output_compare);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::eq_fl;
    use crate::conf::{Conf, OutputType};

    #[test]
    fn computes_rmsd_upper_bound() {
        let a = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let b = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0)];
        assert!(eq_fl(rmsd_upper_bound(&a, &b), 2.0_f64.sqrt()));
    }

    #[test]
    fn replaces_similar_higher_energy_pose() {
        let mut out = vec![OutputType {
            coords: vec![Vec3::new(0.0, 0.0, 0.0)],
            ..OutputType::new(Conf::default(), 2.0)
        }];
        let better = OutputType {
            coords: vec![Vec3::new(0.1, 0.0, 0.0)],
            ..OutputType::new(Conf::default(), 1.0)
        };
        add_to_output_container(&mut out, better, 0.5, 10);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].e, 1.0);
    }
}
