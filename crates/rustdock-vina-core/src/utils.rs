use std::path::{Path, PathBuf};

use crate::file::{open_input, FileError};

pub fn separator() -> char {
    std::path::MAIN_SEPARATOR
}

pub fn make_path(value: &str) -> PathBuf {
    PathBuf::from(value)
}

pub fn doing(value: &str, verbosity: i32, level: i32) {
    if verbosity > level {
        print!("{value} ... ");
        use std::io::Write;
        let _ = std::io::stdout().flush();
    }
}

pub fn done(verbosity: i32, level: i32) {
    if verbosity > level {
        println!("done.");
    }
}

pub fn default_output(input_name: &str) -> String {
    let stem = input_name.strip_suffix(".pdbqt").unwrap_or(input_name);
    format!("{stem}_out.pdbqt")
}

pub fn default_output_instance(input_name: &str, idx: i32) -> String {
    let stem = input_name.strip_suffix(".pdbqt").unwrap_or(input_name);
    format!("{stem}_instance{idx}_out.pdbqt")
}

pub fn default_output_in_dir(input_name: &str, directory_pathname: &str) -> String {
    Path::new(directory_pathname)
        .join(default_output(input_name))
        .to_string_lossy()
        .into_owned()
}

pub fn default_output_instance_in_dir(
    input_name: &str,
    directory_pathname: &str,
    idx: i32,
) -> String {
    Path::new(directory_pathname)
        .join(default_output_instance(input_name, idx))
        .to_string_lossy()
        .into_owned()
}

pub fn is_directory(directory_pathname: &str) -> bool {
    Path::new(directory_pathname).is_dir()
}

pub fn get_filename(value: &str) -> String {
    Path::new(value)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| value.to_string())
}

pub fn get_file_contents(filename: &str) -> Result<String, FileError> {
    use std::io::Read;

    let mut input = open_input(filename)?;
    let mut contents = String::new();
    input
        .read_to_string(&mut contents)
        .map_err(|source| FileError {
            name: PathBuf::from(filename),
            input: true,
            source,
        })?;
    Ok(contents)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_outputs_match_reference_suffixes() {
        assert_eq!(default_output("ligand.pdbqt"), "ligand_out.pdbqt");
        assert_eq!(
            default_output_instance("ligand.pdbqt", 2),
            "ligand_instance2_out.pdbqt"
        );
    }

    #[test]
    fn extracts_filename() {
        assert_eq!(get_filename("/tmp/a/b.pdbqt"), "b.pdbqt");
    }
}
