use std::io::{self, IsTerminal, Write};

pub struct Output {
    enabled: bool,
    color: bool,
}
impl Output {
    pub fn new(verbosity: u32) -> Self {
        Self {
            enabled: verbosity > 0,
            color: io::stdout().is_terminal()
                && std::env::var_os("NO_COLOR").is_none()
                && std::env::var("TERM").as_deref() != Ok("dumb"),
        }
    }
    fn accent(&self, text: &str) -> String {
        if self.color {
            format!("\x1b[1;36m{text}\x1b[0m")
        } else {
            text.into()
        }
    }
    pub fn heading(&self, text: &str) {
        if self.enabled {
            println!("\n{}\n{}", self.accent(text), "─".repeat(64));
        }
    }
    pub fn field(&self, label: &str, value: impl std::fmt::Display) {
        if self.enabled {
            println!("  {label:<16} {value}");
        }
    }
    pub fn stage(&self, text: &str) {
        if self.enabled {
            println!("\n  {} {text}", self.accent("›"));
            let _ = io::stdout().flush();
        }
    }
    pub fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        if self.enabled {
            println!("{}", table(headers, rows));
        }
    }
    pub fn poses(&self, rows: &[(f64, f64, f64)]) {
        if !self.enabled {
            return;
        }
        self.heading("Docking results");
        println!("{}", pose_table(rows));
        if let Some(row) = rows.first() {
            self.field("Best affinity", format!("{:.3} kcal/mol", row.0));
        }
        println!("  RMSD bounds compare each pose with the best predicted pose.");
    }
}

fn table(headers: &[&str], rows: &[Vec<String>]) -> String {
    let widths: Vec<_> = headers
        .iter()
        .enumerate()
        .map(|(i, header)| {
            rows.iter()
                .map(|row| row[i].chars().count())
                .max()
                .unwrap_or(0)
                .max(header.chars().count())
        })
        .collect();
    let border = format!(
        "  +{}+",
        widths
            .iter()
            .map(|w| "-".repeat(w + 2))
            .collect::<Vec<_>>()
            .join("+")
    );
    let line = |cells: &[String]| {
        format!(
            "  | {} |",
            cells
                .iter()
                .zip(&widths)
                .map(|(cell, width)| format!("{cell:>width$}"))
                .collect::<Vec<_>>()
                .join(" | ")
        )
    };
    let mut lines = vec![
        border.clone(),
        line(&headers.iter().map(|h| h.to_string()).collect::<Vec<_>>()),
        border.clone(),
    ];
    lines.extend(rows.iter().map(|row| line(row)));
    lines.push(border);
    lines.join("\n")
}

fn pose_table(rows: &[(f64, f64, f64)]) -> String {
    let border = "  +------+---------------------+--------------+--------------+";
    let mut text = format!(
        "{border}\n  | Pose | Affinity (kcal/mol) | RMSD LB (Å)  | RMSD UB (Å)  |\n{border}\n"
    );
    for (i, (energy, lower, upper)) in rows.iter().enumerate() {
        text.push_str(&format!(
            "  | {:>4} | {:>19.3} | {:>12.3} | {:>12.3} |\n",
            i + 1,
            energy,
            lower,
            upper
        ));
    }
    text.push_str(border);
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn benchmark_table_aligns_unicode_and_large_values() {
        let text = table(
            &["Time (s)", "RMSD (Å)"],
            &[
                vec!["123456.789".into(), "—".into()],
                vec!["0.001".into(), "3.000".into()],
            ],
        );
        let lengths: Vec<_> = text.lines().map(|l| l.chars().count()).collect();
        assert!(lengths.iter().all(|n| *n == lengths[0]));
        assert!(!text.contains('\x1b'));
    }
    #[test]
    fn pose_table_labels_units_and_keeps_columns_aligned() {
        let table = pose_table(&[(-6.215, 0.0, 0.0), (-5.100, 1.234, 2.345)]);
        assert!(table.contains("Affinity (kcal/mol)"));
        assert!(table.contains("-6.215"));
        assert!(table.contains("1.234"));
        let lengths: Vec<_> = table.lines().map(|l| l.chars().count()).collect();
        assert!(lengths.iter().all(|n| *n == lengths[0]));
        assert!(!table.contains('\x1b'));
    }
}
