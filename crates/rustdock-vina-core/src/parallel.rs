pub fn parallel_for<F>(size: usize, mut f: F)
where
    F: FnMut(usize),
{
    for i in 0..size {
        f(i);
    }
}

pub fn parallel_iter<T, F>(values: &mut [T], mut f: F)
where
    F: FnMut(&mut T),
{
    for value in values {
        f(value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn visits_all_indices() {
        let mut values = Vec::new();
        parallel_for(4, |i| values.push(i));
        assert_eq!(values, vec![0, 1, 2, 3]);
    }
}
