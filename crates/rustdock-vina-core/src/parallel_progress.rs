use crate::incrementable::Incrementable;

pub struct ParallelProgress<'a> {
    count: u64,
    value: u64,
    callback: Option<Box<dyn FnMut(f64) + 'a>>,
}

impl<'a> ParallelProgress<'a> {
    pub fn new(callback: Option<Box<dyn FnMut(f64) + 'a>>) -> Self {
        Self {
            count: 0,
            value: 0,
            callback,
        }
    }

    pub fn init(&mut self, n: u64) {
        self.count = n;
        self.value = 0;
    }

    pub fn value(&self) -> u64 {
        self.value
    }
    pub fn advance(&mut self, steps: u64) {
        self.value += steps;
        if let Some(callback) = &mut self.callback {
            if self.count > 0 {
                callback(self.value as f64 / self.count as f64);
            }
        }
    }
}

impl Incrementable for ParallelProgress<'_> {
    fn increment(&mut self) {
        self.advance(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn increments_count() {
        let mut progress = ParallelProgress::new(None);
        progress.init(2);
        progress.increment();
        assert_eq!(progress.value(), 1);
    }
}
