use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct FileError {
    pub name: PathBuf,
    pub input: bool,
    pub source: io::Error,
}

impl std::fmt::Display for FileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let direction = if self.input { "reading" } else { "writing" };
        write!(f, "could not open {:?} for {}", self.name, direction)
    }
}

impl std::error::Error for FileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

pub fn open_input(path: impl AsRef<Path>) -> Result<File, FileError> {
    let path = path.as_ref();
    File::open(path).map_err(|source| FileError {
        name: path.to_path_buf(),
        input: true,
        source,
    })
}

pub fn create_output(path: impl AsRef<Path>) -> Result<File, FileError> {
    let path = path.as_ref();
    File::create(path).map_err(|source| FileError {
        name: path.to_path_buf(),
        input: false,
        source,
    })
}
