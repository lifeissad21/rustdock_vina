#[macro_export]
macro_rules! vina_check {
    ($condition:expr $(,)?) => {
        assert!($condition)
    };
}

#[macro_export]
macro_rules! vina_show {
    ($value:expr $(,)?) => {
        println!("{} = {:?}", stringify!($value), $value)
    };
}

#[macro_export]
macro_rules! vina_eshow {
    ($value:expr $(,)?) => {
        eprintln!("{} = {:?}", stringify!($value), $value)
    };
}

pub fn range(end: usize) -> std::ops::Range<usize> {
    0..end
}

pub fn range_from(begin: usize, end: usize) -> std::ops::Range<usize> {
    begin..end
}

#[cfg(test)]
mod tests {
    #[test]
    fn range_helpers_match_loop_bounds() {
        assert_eq!(super::range(3).collect::<Vec<_>>(), vec![0, 1, 2]);
        assert_eq!(super::range_from(1, 3).collect::<Vec<_>>(), vec![1, 2]);
    }
}
