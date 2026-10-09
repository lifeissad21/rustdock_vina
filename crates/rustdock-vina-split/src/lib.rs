use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SplitModel {
    pub ligand: Vec<String>,
    pub flex: Vec<String>,
}

pub fn default_prefix(input_name: &str, add: &str) -> String {
    let stem = input_name.strip_suffix(".pdbqt").unwrap_or(input_name);
    format!("{stem}{add}")
}

pub fn parse_multimodel_pdbqt_reader<R: BufRead>(reader: R) -> Result<Vec<SplitModel>, String> {
    let mut models = Vec::new();
    let mut parsing_model = false;
    let mut parsing_ligand = true;

    for (line_index, line) in reader.lines().enumerate() {
        let count = line_index + 1;
        let line = line.map_err(|err| err.to_string())?;
        if line.starts_with("MODEL") {
            if parsing_model || !parsing_ligand {
                return Err(format!("Misplaced MODEL tag at line {count}."));
            }
            models.push(SplitModel::default());
            parsing_model = true;
        } else if line.starts_with("ENDMDL") {
            if !parsing_model || !parsing_ligand {
                return Err(format!("Misplaced ENDMDL tag at line {count}."));
            }
            parsing_model = false;
        } else if line.starts_with("BEGIN_RES") {
            if !parsing_model || !parsing_ligand {
                return Err(format!("Misplaced BEGIN_RES tag at line {count}."));
            }
            parsing_ligand = false;
            models.last_mut().expect("model exists").flex.push(line);
        } else if line.starts_with("END_RES") {
            if !parsing_model || parsing_ligand {
                return Err(format!("Misplaced END_RES tag at line {count}."));
            }
            parsing_ligand = true;
            models.last_mut().expect("model exists").flex.push(line);
        } else {
            if !parsing_model {
                return Err(format!("Input occurs outside MODEL at line {count}."));
            }
            if parsing_ligand {
                models.last_mut().expect("model exists").ligand.push(line);
            } else {
                models.last_mut().expect("model exists").flex.push(line);
            }
        }
    }

    if parsing_model {
        return Err("Missing ENDMDL tag at line EOF.".to_string());
    }
    Ok(models)
}

pub fn parse_multimodel_pdbqt(path: impl AsRef<Path>) -> Result<Vec<SplitModel>, String> {
    let input = File::open(path).map_err(|err| err.to_string())?;
    parse_multimodel_pdbqt_reader(BufReader::new(input))
}

pub fn write_pdbqt(lines: &[String], name: impl AsRef<Path>) -> Result<(), String> {
    if lines.is_empty() {
        return Ok(());
    }
    let mut out = File::create(name).map_err(|err| err.to_string())?;
    for line in lines {
        writeln!(out, "{line}").map_err(|err| err.to_string())?;
    }
    Ok(())
}

pub fn output_names(count: usize, ligand_prefix: &str, flex_prefix: &str) -> Vec<(String, String)> {
    let width = count.to_string().len();
    (1..=count)
        .map(|index| {
            let add = format!("{index:0width$}.pdbqt");
            (
                format!("{ligand_prefix}{add}"),
                format!("{flex_prefix}{add}"),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ligand_and_flex_blocks() {
        let input = "MODEL 1\nLIG\nBEGIN_RES\nFLEX\nEND_RES\nENDMDL\n";
        let models = parse_multimodel_pdbqt_reader(input.as_bytes()).unwrap();
        assert_eq!(models.len(), 1);
        assert_eq!(models[0].ligand, vec!["LIG"]);
        assert_eq!(models[0].flex, vec!["BEGIN_RES", "FLEX", "END_RES"]);
    }

    #[test]
    fn derives_prefix_and_numbered_outputs() {
        assert_eq!(default_prefix("x.pdbqt", "_ligand_"), "x_ligand_");
        assert_eq!(
            output_names(12, "lig_", "flex_")[0],
            ("lig_01.pdbqt".to_string(), "flex_01.pdbqt".to_string())
        );
    }
}
