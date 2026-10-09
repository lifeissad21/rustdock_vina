pub trait Incrementable {
    fn increment(&mut self);
}

impl<F> Incrementable for F
where
    F: FnMut(),
{
    fn increment(&mut self) {
        self();
    }
}
