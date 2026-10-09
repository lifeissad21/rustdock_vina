//! PDBQT parsing and torsion-tree assembly, following reference/src/lib/parse_pdbqt.cpp.
use crate::atom::{Atom, AtomIndex};
use crate::atom_constants::*;
use crate::atom_type::AtomTyping;
use crate::common::{Vec3, EPSILON_FL, ZERO_VEC};
use crate::model::{Context, DistanceType, DistanceTypeMatrix, Ligand, Model, ParsedLine, Residue};
use crate::parse_error::PdbqtParseError as Error;
use crate::tree::{
    AtomRange, Branch, FirstSegment, FlexibleBody, Frame, MainBranch, RigidBody, Segment,
};
use std::collections::HashSet;

#[derive(Clone)]
struct ParsedAtom {
    atom: Atom,
    serial: usize,
    context: usize,
    children: Vec<Node>,
}
#[derive(Default, Clone)]
struct Node {
    atoms: Vec<ParsedAtom>,
    immobile: Option<usize>,
    axis_begin: Option<AtomIndex>,
    axis_end: Option<AtomIndex>,
}

fn error(line: &str, message: &str) -> Error {
    Error::with_line(message, line)
}
fn field<T: std::str::FromStr>(line: &str, start: usize, end: usize) -> Result<T, Error> {
    line.get(start..end)
        .ok_or_else(|| error(line, "This line is too short or is not ASCII."))?
        .trim()
        .parse()
        .map_err(|_| error(line, "Invalid numeric field."))
}
pub fn parse_pdbqt_atom_string(line: &str) -> Result<(usize, Atom), Error> {
    if !line.is_ascii() {
        return Err(error(line, "PDBQT atom records must be ASCII."));
    }
    let serial = field(line, 6, 11)?;
    let coords = Vec3::new(
        field(line, 30, 38)?,
        field(line, 38, 46)?,
        field(line, 46, 54)?,
    );
    let charge = if line.get(68..76).is_some_and(|s| s.trim().is_empty()) {
        0.0
    } else {
        field(line, 68, 76)?
    };
    if !coords.data.iter().all(|v| v.is_finite()) || !f64::is_finite(charge) {
        return Err(error(line, "Coordinates and charges must be finite."));
    }
    let name = line
        .get(77..)
        .ok_or_else(|| error(line, "Missing atom type."))?
        .trim();
    let mut atom = Atom {
        coords,
        ..Atom::default()
    };
    atom.base.charge = charge;
    atom.base.atom_type.ad = string_to_ad_type(name);
    if is_non_ad_metal_name(name) {
        atom.base.atom_type.xs = XS_TYPE_MET_D;
    }
    if !atom.base.atom_type.acceptable_type() {
        return Err(error(line, "Invalid AutoDock atom type (case-sensitive)."));
    }
    Ok((serial, atom))
}
fn ignored(line: &str) -> bool {
    line.trim().is_empty() || line.starts_with("REMARK") || line.starts_with("WARNING")
}
fn atom_line(line: &str) -> bool {
    line.starts_with("ATOM  ") || line.starts_with("HETATM")
}
fn numbers(line: &str, count: usize) -> Result<Vec<usize>, Error> {
    let values = line
        .split_whitespace()
        .skip(1)
        .map(|s| s.parse::<usize>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| error(line, "Invalid atom number or TORSDOF."))?;
    if values.len() != count {
        return Err(error(line, "Wrong number of arguments."));
    }
    Ok(values)
}
struct Parser<'a> {
    lines: Vec<&'a str>,
    pos: usize,
    context: Context,
    serials: HashSet<usize>,
}
impl<'a> Parser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            lines: text.lines().collect(),
            pos: 0,
            context: vec![],
            serials: HashSet::new(),
        }
    }
    fn next(&mut self) -> Option<&'a str> {
        let line = *self.lines.get(self.pos)?;
        self.pos += 1;
        self.context.push(ParsedLine::new(line, None));
        Some(line)
    }
    fn atom(&mut self, line: &str) -> Result<ParsedAtom, Error> {
        let (serial, atom) = parse_pdbqt_atom_string(line)?;
        if !self.serials.insert(serial) {
            return Err(error(line, "Duplicate atom number."));
        }
        Ok(ParsedAtom {
            atom,
            serial,
            context: self.context.len() - 1,
            children: vec![],
        })
    }
    fn branch(&mut self, from: usize, to: usize) -> Result<Node, Error> {
        let mut node = Node::default();
        while let Some(line) = self.next() {
            if ignored(line) {
                continue;
            }
            if atom_line(line) {
                let atom = self.atom(line)?;
                if atom.serial == to {
                    node.immobile = Some(node.atoms.len());
                }
                node.atoms.push(atom);
            } else if line.starts_with("BRANCH ") {
                self.child(&mut node, line)?;
            } else if line.starts_with("ENDBRANCH ") {
                if numbers(line, 2)? != [from, to] {
                    return Err(error(line, "Inconsistent ENDBRANCH."));
                }
                let index = node
                    .immobile
                    .ok_or_else(|| error(line, "Branch endpoint atom is missing."))?;
                if node.atoms[index]
                    .atom
                    .coords
                    .data
                    .iter()
                    .any(|v| !v.is_finite())
                {
                    return Err(error(line, "Invalid branch endpoint."));
                }
                return Ok(node);
            } else {
                return Err(error(line, "Unknown or inappropriate tag in branch."));
            }
        }
        Err(Error::new("Missing ENDBRANCH."))
    }
    fn child(&mut self, node: &mut Node, line: &str) -> Result<(), Error> {
        let ns = numbers(line, 2)?;
        let parent = node
            .atoms
            .iter_mut()
            .find(|a| a.serial == ns[0])
            .ok_or_else(|| error(line, "Parent atom is missing in this branch."))?;
        let child = self.branch(ns[0], ns[1])?;
        let origin = child.atoms[child.immobile.unwrap()].atom.coords;
        if (origin - parent.atom.coords).norm() < EPSILON_FL {
            return Err(error(line, "Rotatable bond has zero length."));
        }
        parent.children.push(child);
        Ok(())
    }
    fn molecule(&mut self, residue: bool) -> Result<(Node, usize), Error> {
        let mut node = Node::default();
        let mut root = false;
        while let Some(line) = self.next() {
            if ignored(line) {
                continue;
            }
            if line == "ROOT" {
                root = true;
                break;
            }
            return Err(error(line, "Expected ROOT."));
        }
        if !root {
            return Err(Error::new("Missing ROOT."));
        }
        let mut ended = false;
        while let Some(line) = self.next() {
            if ignored(line) {
                continue;
            }
            if atom_line(line) {
                node.atoms.push(self.atom(line)?);
            } else if line == "ENDROOT" {
                ended = true;
                break;
            } else {
                return Err(error(line, "Expected atom or ENDROOT."));
            }
        }
        if !ended || node.atoms.is_empty() {
            return Err(Error::new("Missing ENDROOT or no root atoms."));
        }
        let mut torsdof = None;
        while let Some(line) = self.next() {
            if ignored(line) {
                continue;
            }
            if line.starts_with("BRANCH ") {
                self.child(&mut node, line)?;
            } else if !residue && line.starts_with("TORSDOF ") {
                if torsdof.is_some() {
                    return Err(error(line, "Duplicate TORSDOF."));
                }
                torsdof = Some(numbers(line, 1)?[0]);
            } else if residue && line.starts_with("END_RES") {
                return Ok((node, 0));
            } else {
                return Err(error(line, "Unknown or inappropriate ligand/residue tag."));
            }
        }
        if residue {
            return Err(Error::new("Missing END_RES."));
        }
        Ok((node, torsdof.ok_or_else(|| Error::new("Missing TORSDOF."))?))
    }
}

// Fixed endpoints belong to the parent frame; only the remaining branch atoms rotate.
#[derive(Default)]
struct Assembly {
    movable: Vec<(Atom, Vec3)>,
    inflex: Vec<(Atom, Vec3)>,
    bonds: Vec<(AtomIndex, AtomIndex, DistanceType)>,
}
impl Assembly {
    fn insert(
        &mut self,
        atom: &mut ParsedAtom,
        context: &mut Context,
        origin: Vec3,
        inflex: bool,
    ) -> AtomIndex {
        let index = AtomIndex::new(
            if inflex {
                self.inflex.len()
            } else {
                self.movable.len()
            },
            inflex,
        );
        for child in &mut atom.children {
            child.axis_begin = Some(index);
        }
        let lab = atom.atom.coords;
        let mut a = atom.atom.clone();
        a.coords = if inflex { ZERO_VEC } else { lab - origin };
        if inflex {
            self.inflex.push((a, lab));
        } else {
            self.movable.push((a, lab));
            context[atom.context].atom_index = Some(index.i);
        }
        index
    }
    fn immobile(&mut self, node: &mut Node, context: &mut Context, origin: Vec3, inflex: bool) {
        let index = node.immobile.unwrap();
        node.axis_end = Some(self.insert(&mut node.atoms[index], context, origin, inflex));
    }
    fn fill(
        &mut self,
        node: &mut Node,
        context: &mut Context,
        frame: &Frame,
    ) -> (AtomRange, Vec<Branch>) {
        let begin = self.movable.len();
        for (i, atom) in node.atoms.iter_mut().enumerate() {
            if Some(i) != node.immobile {
                self.insert(atom, context, frame.origin(), false);
            }
            for child in &mut atom.children {
                self.immobile(child, context, frame.origin(), false);
            }
        }
        let end = self.movable.len();
        for axis in [node.axis_begin, node.axis_end].into_iter().flatten() {
            for i in begin..end {
                if axis != AtomIndex::new(i, false) {
                    self.bonds
                        .push((axis, AtomIndex::new(i, false), DistanceType::Fixed));
                }
            }
        }
        if let (Some(a), Some(b)) = (node.axis_begin, node.axis_end) {
            self.bonds.push((a, b, DistanceType::Rotor));
        }
        for i in begin..end {
            for j in i + 1..end {
                self.bonds.push((
                    AtomIndex::new(i, false),
                    AtomIndex::new(j, false),
                    DistanceType::Fixed,
                ));
            }
        }
        let mut branches = vec![];
        for atom in &mut node.atoms {
            for child in &mut atom.children {
                let index = child.immobile.unwrap();
                if child.atoms.len() == 1 && child.atoms[index].children.is_empty() {
                    continue;
                }
                let origin = child.atoms[index].atom.coords;
                let mut segment = Segment::new(origin, 0, 0, atom.atom.coords, frame);
                let (range, children) =
                    self.fill(child, context, &segment.axis_frame.atom_frame.frame);
                segment.axis_frame.atom_frame.range = range;
                branches.push(Branch {
                    node: segment,
                    children,
                });
            }
        }
        (AtomRange::new(begin, end), branches)
    }
    fn finish(self, model: &mut Model) {
        let n = self.movable.len();
        for (atom, lab) in self.movable {
            model.push_atom(atom, true);
            *model.coords.last_mut().unwrap() = lab;
        }
        for (atom, lab) in self.inflex {
            model.push_atom(atom, false);
            *model.coords.last_mut().unwrap() = lab;
        }
        let mut mobility = DistanceTypeMatrix::new(model.num_atoms(), DistanceType::Variable);
        for i in n..model.num_atoms() {
            for j in i + 1..model.num_atoms() {
                *mobility.get_mut(i, j) = DistanceType::Fixed;
            }
        }
        for (a, b, t) in self.bonds {
            let i = a.i + if a.in_grid { n } else { 0 };
            let j = b.i + if b.in_grid { n } else { 0 };
            if i != j {
                *mobility.get_mut(i.min(j), i.max(j)) = t;
            }
        }
        model.initialize(&mobility);
    }
}
pub fn parse_ligand_pdbqt_from_string(text: &str, typing: AtomTyping) -> Result<Model, Error> {
    let mut parser = Parser::new(text);
    let (mut node, torsdof) = parser.molecule(false)?;
    let mut body = FlexibleBody::new(RigidBody::new(node.atoms[0].atom.coords, 0, 0));
    let mut assembly = Assembly::default();
    let (range, children) =
        assembly.fill(&mut node, &mut parser.context, &body.node.atom_frame.frame);
    body.node.atom_frame.range = range;
    body.children = children;
    let mut ligand = Ligand::new(body, torsdof);
    ligand.context = parser.context;
    let mut model = Model::new(typing);
    model.ligands.push(ligand);
    assembly.finish(&mut model);
    Ok(model)
}
pub fn parse_ligand_pdbqt_from_file(path: &str, typing: AtomTyping) -> Result<Model, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    parse_ligand_pdbqt_from_string(&text, typing).map_err(|e| e.to_string())
}
pub fn parse_receptor_pdbqt_strings(
    rigid: &str,
    flex: &str,
    typing: AtomTyping,
) -> Result<Model, Error> {
    let mut model = Model::new(typing);
    for line in rigid.lines() {
        if ignored(line) || line.starts_with("TER") || line.starts_with("END") {
            continue;
        }
        if !atom_line(line) {
            return Err(error(line, "Unknown or inappropriate rigid receptor tag."));
        }
        model.push_grid_atom(parse_pdbqt_atom_string(line)?.1);
    }
    let mut parser = Parser::new(flex);
    let mut assembly = Assembly::default();
    while let Some(line) = parser.next() {
        if ignored(line) {
            continue;
        }
        if !line.starts_with("BEGIN_RES") {
            return Err(error(line, "Expected BEGIN_RES."));
        }
        parser.serials.clear();
        let (mut node, _) = parser.molecule(true)?;
        for atom in &mut node.atoms {
            assembly.insert(atom, &mut parser.context, ZERO_VEC, true);
            for child in &mut atom.children {
                assembly.immobile(child, &mut parser.context, ZERO_VEC, true);
            }
        }
        for atom in &mut node.atoms {
            for child in &mut atom.children {
                let index = child.immobile.unwrap();
                if child.atoms.len() == 1 && child.atoms[index].children.is_empty() {
                    continue;
                }
                let mut branch = MainBranch::new(FirstSegment::new(
                    child.atoms[index].atom.coords,
                    0,
                    0,
                    atom.atom.coords,
                ));
                let (range, children) = assembly.fill(
                    child,
                    &mut parser.context,
                    &branch.node.axis_frame.atom_frame.frame,
                );
                branch.node.axis_frame.atom_frame.range = range;
                branch.children = children;
                model.flex.push(Residue::new(branch));
            }
        }
    }
    model.flex_context = parser.context;
    assembly.finish(&mut model);
    Ok(model)
}
pub fn parse_receptor_pdbqt(rigid: &str, flex: &str, typing: AtomTyping) -> Result<Model, String> {
    let read = |p: &str| {
        if p.is_empty() {
            Ok(String::new())
        } else {
            std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}"))
        }
    };
    parse_receptor_pdbqt_strings(&read(rigid)?, &read(flex)?, typing).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_ligand_roundtrips_and_has_rotors() {
        let input =
            include_str!("../../../reference/example/basic_docking/solution/1iep_ligand.pdbqt");
        let mut model = parse_ligand_pdbqt_from_string(input, AtomTyping::Xs).unwrap();
        let coords = model.coords.clone();
        model.set(&model.get_initial_conf());
        assert_eq!(
            model.num_atoms(),
            input.lines().filter(|l| atom_line(l)).count()
        );
        assert!(model.get_size().ligands[0] > 0);
        for (a, b) in coords.iter().zip(&model.coords) {
            assert!((*a - *b).norm() < 1e-10);
        }
        assert!(model
            .atoms
            .iter()
            .any(|a| a.bonds.iter().any(|b| b.rotatable)));
    }
    #[test]
    fn malformed_inputs_return_errors() {
        for text in [
            "",
            "ROOT\nENDROOT\nTORSDOF 0",
            "ROOT\nATOM\nENDROOT\nTORSDOF 0",
            "MODEL 1",
        ] {
            assert!(parse_ligand_pdbqt_from_string(text, AtomTyping::Xs).is_err());
        }
        assert!(parse_pdbqt_atom_string("ATOM").is_err());
    }
}
