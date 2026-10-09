use crate::atom::Atom;
use crate::common::{cross_product, Fl, Mat3, Vec3, EPSILON_FL, ZERO_VEC};
use crate::conf::{LigandChange, LigandConf, ResidueChange, ResidueConf, RigidChange};
use crate::quaternion::{angle_axis_to_quaternion, quaternion_to_r3, Quaternion, QT_IDENTITY};

#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    origin: Vec3,
    orientation_m: Mat3,
    orientation_q: Quaternion,
}

impl Frame {
    pub fn new(origin: Vec3) -> Self {
        Self {
            origin,
            orientation_m: quaternion_to_r3(QT_IDENTITY),
            orientation_q: QT_IDENTITY,
        }
    }

    pub fn local_to_lab(&self, local_coords: Vec3) -> Vec3 {
        self.origin + self.orientation_m.mul_vec(local_coords)
    }

    pub fn local_to_lab_direction(&self, local_direction: Vec3) -> Vec3 {
        self.orientation_m.mul_vec(local_direction)
    }

    pub fn orientation(&self) -> Quaternion {
        self.orientation_q
    }

    pub fn origin(&self) -> Vec3 {
        self.origin
    }

    fn set_origin(&mut self, origin: Vec3) {
        self.origin = origin;
    }

    fn set_orientation(&mut self, q: Quaternion) {
        self.orientation_q = q;
        self.orientation_m = quaternion_to_r3(q);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomRange {
    pub begin: usize,
    pub end: usize,
}

impl AtomRange {
    pub const fn new(begin: usize, end: usize) -> Self {
        Self { begin, end }
    }

    pub fn len(self) -> usize {
        self.end - self.begin
    }

    pub fn is_empty(self) -> bool {
        self.begin == self.end
    }

    pub fn transform<F>(&mut self, mut f: F)
    where
        F: FnMut(usize) -> usize,
    {
        let diff = self.end - self.begin;
        self.begin = f(self.begin);
        self.end = self.begin + diff;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomFrame {
    pub frame: Frame,
    pub range: AtomRange,
}

impl AtomFrame {
    pub fn new(origin: Vec3, begin: usize, end: usize) -> Self {
        Self {
            frame: Frame::new(origin),
            range: AtomRange::new(begin, end),
        }
    }

    pub fn set_coords(&self, atoms: &[Atom], coords: &mut [Vec3]) {
        assert!(self.range.end <= atoms.len());
        assert!(self.range.end <= coords.len());
        for i in self.range.begin..self.range.end {
            coords[i] = self.frame.local_to_lab(atoms[i].coords);
        }
    }

    pub fn sum_force_and_torque(&self, coords: &[Vec3], forces: &[Vec3]) -> ForceTorque {
        assert!(self.range.end <= coords.len());
        assert!(self.range.end <= forces.len());
        let mut tmp = ForceTorque::default();
        for i in self.range.begin..self.range.end {
            tmp.force += forces[i];
            tmp.torque += cross_product(coords[i] - self.frame.origin(), forces[i]);
        }
        tmp
    }

    pub fn transform<F>(&mut self, f: F)
    where
        F: FnMut(usize) -> usize,
    {
        self.range.transform(f);
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForceTorque {
    pub force: Vec3,
    pub torque: Vec3,
}

impl Default for ForceTorque {
    fn default() -> Self {
        Self {
            force: ZERO_VEC,
            torque: ZERO_VEC,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RigidBody {
    pub atom_frame: AtomFrame,
}

impl RigidBody {
    pub fn new(origin: Vec3, begin: usize, end: usize) -> Self {
        Self {
            atom_frame: AtomFrame::new(origin, begin, end),
        }
    }

    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], c: &crate::conf::RigidConf) {
        self.atom_frame.frame.set_origin(c.position);
        self.atom_frame.frame.set_orientation(c.orientation);
        self.atom_frame.set_coords(atoms, coords);
    }

    pub fn set_derivative(&self, force_torque: ForceTorque, c: &mut RigidChange) {
        c.position = force_torque.force;
        c.orientation = force_torque.torque;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AxisFrame {
    pub atom_frame: AtomFrame,
    pub axis: Vec3,
}

impl AxisFrame {
    pub fn new(origin: Vec3, begin: usize, end: usize, axis_root: Vec3) -> Self {
        let diff = origin - axis_root;
        let norm = diff.norm();
        assert!(norm >= EPSILON_FL);
        Self {
            atom_frame: AtomFrame::new(origin, begin, end),
            axis: (1.0 / norm) * diff,
        }
    }

    pub fn set_derivative(&self, force_torque: ForceTorque, c: &mut Fl) {
        *c = force_torque.torque.dot(self.axis);
    }
}

pub trait TreeNode: Clone {
    fn atom_frame(&self) -> &AtomFrame;
    fn atom_frame_mut(&mut self) -> &mut AtomFrame;
    fn set_conf(
        &mut self,
        parent: &Frame,
        atoms: &[Atom],
        coords: &mut [Vec3],
        torsions: &mut dyn Iterator<Item = Fl>,
    );
    fn set_derivative(&self, force_torque: ForceTorque, c: &mut Fl);

    fn count_torsions(&self, out: &mut usize) {
        *out += 1;
    }

    fn transform<F>(&mut self, f: F)
    where
        F: FnMut(usize) -> usize,
    {
        self.atom_frame_mut().transform(f);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub axis_frame: AxisFrame,
    pub relative_axis: Vec3,
    pub relative_origin: Vec3,
}

impl Segment {
    pub fn new(origin: Vec3, begin: usize, end: usize, axis_root: Vec3, parent: &Frame) -> Self {
        assert_eq!(parent.orientation(), QT_IDENTITY);
        let axis_frame = AxisFrame::new(origin, begin, end, axis_root);
        Self {
            relative_axis: axis_frame.axis,
            relative_origin: origin - parent.origin(),
            axis_frame,
        }
    }
}

impl TreeNode for Segment {
    fn atom_frame(&self) -> &AtomFrame {
        &self.axis_frame.atom_frame
    }

    fn atom_frame_mut(&mut self) -> &mut AtomFrame {
        &mut self.axis_frame.atom_frame
    }

    fn set_conf(
        &mut self,
        parent: &Frame,
        atoms: &[Atom],
        coords: &mut [Vec3],
        torsions: &mut dyn Iterator<Item = Fl>,
    ) {
        let torsion = torsions.next().expect("segment torsion missing");
        let origin = parent.local_to_lab(self.relative_origin);
        self.axis_frame.atom_frame.frame.set_origin(origin);
        self.axis_frame.axis = parent.local_to_lab_direction(self.relative_axis);
        let mut q = angle_axis_to_quaternion(self.axis_frame.axis, torsion) * parent.orientation();
        q.normalize_approx(1e-6);
        self.axis_frame.atom_frame.frame.set_orientation(q);
        self.axis_frame.atom_frame.set_coords(atoms, coords);
    }

    fn set_derivative(&self, force_torque: ForceTorque, c: &mut Fl) {
        self.axis_frame.set_derivative(force_torque, c);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FirstSegment {
    pub axis_frame: AxisFrame,
}

impl FirstSegment {
    pub fn new(origin: Vec3, begin: usize, end: usize, axis_root: Vec3) -> Self {
        Self {
            axis_frame: AxisFrame::new(origin, begin, end, axis_root),
        }
    }

    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], torsion: Fl) {
        self.axis_frame
            .atom_frame
            .frame
            .set_orientation(angle_axis_to_quaternion(self.axis_frame.axis, torsion));
        self.axis_frame.atom_frame.set_coords(atoms, coords);
    }

    pub fn set_derivative(&self, force_torque: ForceTorque, c: &mut Fl) {
        self.axis_frame.set_derivative(force_torque, c);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tree<T: TreeNode> {
    pub node: T,
    pub children: Vec<Tree<Segment>>,
}

impl<T: TreeNode> Tree<T> {
    pub fn new(node: T) -> Self {
        Self {
            node,
            children: Vec::new(),
        }
    }

    pub fn set_conf(
        &mut self,
        parent: &Frame,
        atoms: &[Atom],
        coords: &mut [Vec3],
        torsions: &mut dyn Iterator<Item = Fl>,
    ) {
        self.node.set_conf(parent, atoms, coords, torsions);
        branches_set_conf(
            &mut self.children,
            self.node.atom_frame().frame.clone(),
            atoms,
            coords,
            torsions,
        );
    }

    pub fn derivative(
        &self,
        coords: &[Vec3],
        forces: &[Vec3],
        torsion_changes: &mut dyn Iterator<Item = &mut Fl>,
    ) -> ForceTorque {
        let mut force_torque = self.node.atom_frame().sum_force_and_torque(coords, forces);
        let d = torsion_changes
            .next()
            .expect("segment derivative slot missing");
        branches_derivative(
            &self.children,
            self.node.atom_frame().frame.origin(),
            coords,
            forces,
            &mut force_torque,
            torsion_changes,
        );
        self.node.set_derivative(force_torque, d);
        force_torque
    }
}

pub type Branch = Tree<Segment>;
pub type Branches = Vec<Branch>;

pub fn branches_set_conf(
    branches: &mut [Branch],
    parent: Frame,
    atoms: &[Atom],
    coords: &mut [Vec3],
    torsions: &mut dyn Iterator<Item = Fl>,
) {
    for branch in branches {
        branch.set_conf(&parent, atoms, coords, torsions);
    }
}

pub fn branches_derivative<'a>(
    branches: &[Branch],
    origin: Vec3,
    coords: &[Vec3],
    forces: &[Vec3],
    out: &mut ForceTorque,
    torsion_changes: &mut dyn Iterator<Item = &'a mut Fl>,
) {
    for branch in branches {
        let force_torque = branch.derivative(coords, forces, torsion_changes);
        out.force += force_torque.force;
        let r = branch.node.atom_frame().frame.origin() - origin;
        out.torque += cross_product(r, force_torque.force) + force_torque.torque;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Heterotree<Node> {
    pub node: Node,
    pub children: Branches,
}

impl<Node> Heterotree<Node> {
    pub fn new(node: Node) -> Self {
        Self {
            node,
            children: Vec::new(),
        }
    }
}

impl Heterotree<RigidBody> {
    pub fn set_ligand_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], c: &LigandConf) {
        self.node.set_conf(atoms, coords, &c.rigid);
        let mut torsions = c.torsions.iter().copied();
        branches_set_conf(
            &mut self.children,
            self.node.atom_frame.frame.clone(),
            atoms,
            coords,
            &mut torsions,
        );
        assert!(torsions.next().is_none());
    }

    pub fn derivative_ligand(&self, coords: &[Vec3], forces: &[Vec3], c: &mut LigandChange) {
        let mut force_torque = self.node.atom_frame.sum_force_and_torque(coords, forces);
        let mut torsions = c.torsions.iter_mut();
        branches_derivative(
            &self.children,
            self.node.atom_frame.frame.origin(),
            coords,
            forces,
            &mut force_torque,
            &mut torsions,
        );
        self.node.set_derivative(force_torque, &mut c.rigid);
        assert!(torsions.next().is_none());
    }
}

impl Heterotree<FirstSegment> {
    pub fn set_residue_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], c: &ResidueConf) {
        let mut torsions = c.torsions.iter().copied();
        let first = torsions.next().expect("first segment torsion missing");
        self.node.set_conf(atoms, coords, first);
        branches_set_conf(
            &mut self.children,
            self.node.axis_frame.atom_frame.frame.clone(),
            atoms,
            coords,
            &mut torsions,
        );
        assert!(torsions.next().is_none());
    }

    pub fn derivative_residue(&self, coords: &[Vec3], forces: &[Vec3], c: &mut ResidueChange) {
        let mut force_torque = self
            .node
            .axis_frame
            .atom_frame
            .sum_force_and_torque(coords, forces);
        let mut torsions = c.torsions.iter_mut();
        let first = torsions
            .next()
            .expect("first segment derivative slot missing");
        branches_derivative(
            &self.children,
            self.node.axis_frame.atom_frame.frame.origin(),
            coords,
            forces,
            &mut force_torque,
            &mut torsions,
        );
        self.node.set_derivative(force_torque, first);
        assert!(torsions.next().is_none());
    }
}

pub trait TorsionTree {
    fn count_torsions(&self, out: &mut usize);
}

impl<T: TreeNode> TorsionTree for Tree<T> {
    fn count_torsions(&self, out: &mut usize) {
        self.node.count_torsions(out);
        for child in &self.children {
            child.count_torsions(out);
        }
    }
}

impl TorsionTree for Heterotree<RigidBody> {
    fn count_torsions(&self, out: &mut usize) {
        for child in &self.children {
            child.count_torsions(out);
        }
    }
}

impl TorsionTree for Heterotree<FirstSegment> {
    fn count_torsions(&self, out: &mut usize) {
        *out += 1;
        for child in &self.children {
            child.count_torsions(out);
        }
    }
}

pub fn count_torsions<T: TorsionTree>(tree: &T) -> usize {
    let mut out = 0;
    tree.count_torsions(&mut out);
    out
}

pub type FlexibleBody = Heterotree<RigidBody>;
pub type MainBranch = Heterotree<FirstSegment>;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct VectorMutable<T> {
    pub items: Vec<T>,
}

impl<T> VectorMutable<T> {
    pub fn new(items: Vec<T>) -> Self {
        Self { items }
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

impl VectorMutable<FlexibleBody> {
    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], confs: &[LigandConf]) {
        assert_eq!(self.items.len(), confs.len());
        for (item, conf) in self.items.iter_mut().zip(confs) {
            item.set_ligand_conf(atoms, coords, conf);
        }
    }

    pub fn derivative(&self, coords: &[Vec3], forces: &[Vec3], changes: &mut [LigandChange]) {
        assert_eq!(self.items.len(), changes.len());
        for (item, change) in self.items.iter().zip(changes) {
            item.derivative_ligand(coords, forces, change);
        }
    }
}

impl VectorMutable<MainBranch> {
    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], confs: &[ResidueConf]) {
        assert_eq!(self.items.len(), confs.len());
        for (item, conf) in self.items.iter_mut().zip(confs) {
            item.set_residue_conf(atoms, coords, conf);
        }
    }

    pub fn derivative(&self, coords: &[Vec3], forces: &[Vec3], changes: &mut [ResidueChange]) {
        assert_eq!(self.items.len(), changes.len());
        for (item, change) in self.items.iter().zip(changes) {
            item.derivative_residue(coords, forces, change);
        }
    }
}

impl<T: TorsionTree> VectorMutable<T> {
    pub fn count_torsions(&self) -> Vec<usize> {
        self.items.iter().map(count_torsions).collect()
    }
}

pub trait TransformRanges {
    fn transform_ranges<F>(&mut self, f: &mut F)
    where
        F: FnMut(usize) -> usize;
}

impl<T: TreeNode> TransformRanges for Tree<T> {
    fn transform_ranges<F>(&mut self, f: &mut F)
    where
        F: FnMut(usize) -> usize,
    {
        self.node.atom_frame_mut().transform(&mut *f);
        for child in &mut self.children {
            child.transform_ranges(f);
        }
    }
}

impl TransformRanges for Heterotree<RigidBody> {
    fn transform_ranges<F>(&mut self, f: &mut F)
    where
        F: FnMut(usize) -> usize,
    {
        self.node.atom_frame.transform(&mut *f);
        for child in &mut self.children {
            child.transform_ranges(f);
        }
    }
}

impl TransformRanges for Heterotree<FirstSegment> {
    fn transform_ranges<F>(&mut self, f: &mut F)
    where
        F: FnMut(usize) -> usize,
    {
        self.node.axis_frame.atom_frame.transform(&mut *f);
        for child in &mut self.children {
            child.transform_ranges(f);
        }
    }
}

pub fn transform_ranges<T, F>(tree: &mut T, mut f: F)
where
    T: TransformRanges,
    F: FnMut(usize) -> usize,
{
    tree.transform_ranges(&mut f);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{eq_fl, eq_vec};
    use crate::conf::{LigandConf, RigidConf};

    fn atom(x: Fl, y: Fl, z: Fl) -> Atom {
        Atom {
            coords: Vec3::new(x, y, z),
            ..Atom::default()
        }
    }

    #[test]
    fn atom_range_transform_preserves_length() {
        let mut range = AtomRange::new(2, 7);
        range.transform(|i| i + 10);
        assert_eq!(range, AtomRange::new(12, 17));
    }

    #[test]
    fn rigid_body_set_conf_transforms_atoms() {
        let atoms = vec![atom(1.0, 0.0, 0.0)];
        let mut coords = vec![ZERO_VEC];
        let mut body = RigidBody::new(ZERO_VEC, 0, 1);
        body.set_conf(
            &atoms,
            &mut coords,
            &RigidConf {
                position: Vec3::new(1.0, 2.0, 3.0),
                orientation: QT_IDENTITY,
            },
        );
        assert!(eq_vec(coords[0], Vec3::new(2.0, 2.0, 3.0)));
    }

    #[test]
    fn segment_set_conf_uses_parent_frame_and_torsion() {
        let atoms = vec![atom(0.0, 1.0, 0.0)];
        let mut coords = vec![ZERO_VEC];
        let parent = Frame::new(Vec3::new(10.0, 0.0, 0.0));
        let mut segment = Segment::new(Vec3::new(11.0, 0.0, 0.0), 0, 1, parent.origin(), &parent);
        let mut torsions = [core::f64::consts::FRAC_PI_2].into_iter();
        segment.set_conf(&parent, &atoms, &mut coords, &mut torsions);
        assert!(eq_vec(
            segment.axis_frame.atom_frame.frame.origin(),
            Vec3::new(11.0, 0.0, 0.0)
        ));
        assert!(eq_fl(coords[0][0], 11.0));
        assert!(eq_fl(coords[0][1], 0.0));
        assert!(eq_fl(coords[0][2], 1.0));
    }

    #[test]
    fn ligand_derivative_writes_rigid_change() {
        let mut body = FlexibleBody::new(RigidBody::new(ZERO_VEC, 0, 1));
        let atoms = vec![atom(1.0, 0.0, 0.0)];
        let mut coords = vec![ZERO_VEC];
        body.set_ligand_conf(&atoms, &mut coords, &LigandConf::default());
        let mut change = LigandChange::default();
        body.derivative_ligand(&coords, &[Vec3::new(0.0, 1.0, 0.0)], &mut change);
        assert!(eq_vec(change.rigid.position, Vec3::new(0.0, 1.0, 0.0)));
        assert!(eq_vec(change.rigid.orientation, Vec3::new(0.0, 0.0, 1.0)));
    }

    #[test]
    fn counts_torsions_through_children() {
        let parent = Frame::new(ZERO_VEC);
        let mut body = FlexibleBody::new(RigidBody::new(ZERO_VEC, 0, 0));
        body.children.push(Branch::new(Segment::new(
            Vec3::new(1.0, 0.0, 0.0),
            0,
            0,
            ZERO_VEC,
            &parent,
        )));
        body.children[0].children.push(Branch::new(Segment::new(
            Vec3::new(2.0, 0.0, 0.0),
            0,
            0,
            Vec3::new(1.0, 0.0, 0.0),
            &Frame::new(Vec3::new(1.0, 0.0, 0.0)),
        )));
        assert_eq!(count_torsions(&body), 2);
    }
}
