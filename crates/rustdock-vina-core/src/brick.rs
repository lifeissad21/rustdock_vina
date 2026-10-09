use crate::common::{vec_distance_sqr, Fl, Vec3};

pub fn closest_between(begin: Fl, end: Fl, x: Fl) -> Fl {
    assert!(begin <= end);
    if x <= begin {
        begin
    } else if x >= end {
        end
    } else {
        x
    }
}

pub fn brick_closest(begin: Vec3, end: Vec3, v: Vec3) -> Vec3 {
    Vec3::new(
        closest_between(begin[0], end[0], v[0]),
        closest_between(begin[1], end[1], v[1]),
        closest_between(begin[2], end[2], v[2]),
    )
}

pub fn brick_distance_sqr(begin: Vec3, end: Vec3, v: Vec3) -> Fl {
    vec_distance_sqr(brick_closest(begin, end, v), v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clamps_to_brick_boundary() {
        let begin = Vec3::new(0.0, 0.0, 0.0);
        let end = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(
            brick_closest(begin, end, Vec3::new(-1.0, 1.0, 4.0)),
            Vec3::new(0.0, 1.0, 3.0)
        );
        assert_eq!(
            brick_distance_sqr(begin, end, Vec3::new(-1.0, 1.0, 4.0)),
            2.0
        );
    }
}
