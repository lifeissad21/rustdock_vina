#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbqtParseError {
    message: String,
    pdbqt_line: Option<String>,
}

impl PdbqtParseError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            pdbqt_line: None,
        }
    }

    pub fn with_line(message: impl Into<String>, pdbqt_line: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            pdbqt_line: Some(pdbqt_line.into()),
        }
    }

    pub fn message(&self) -> &str {
        &self.message
    }

    pub fn pdbqt_line(&self) -> Option<&str> {
        self.pdbqt_line.as_deref()
    }
}

impl std::fmt::Display for PdbqtParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.pdbqt_line {
            Some(line) => write!(
                f,
                "\n\nPDBQT parsing error: {}\n > {}\n",
                self.message, line
            ),
            None => write!(f, "\n\nPDBQT parsing error: {}\n", self.message),
        }
    }
}

impl std::error::Error for PdbqtParseError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_like_reference_error() {
        let error = PdbqtParseError::with_line("bad atom", "ATOM");
        assert_eq!(
            error.to_string(),
            "\n\nPDBQT parsing error: bad atom\n > ATOM\n"
        );
    }
}
