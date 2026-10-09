use crate::common::Fl;

pub fn int_pow<const N: usize>(x: Fl) -> Fl {
    let mut result = 1.0;
    let mut i = 0;
    while i < N {
        result *= x;
        i += 1;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_recursive_template_results() {
        assert_eq!(int_pow::<0>(5.0), 1.0);
        assert_eq!(int_pow::<3>(2.0), 8.0);
    }
}
