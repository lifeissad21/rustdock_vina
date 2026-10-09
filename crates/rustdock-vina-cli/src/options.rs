use std::collections::BTreeMap;

pub const FLAGS: &[&str] = &[
    "check",
    "help",
    "help_advanced",
    "version",
    "score_only",
    "local_only",
    "randomize_only",
    "no_refine",
    "force_even_voxels",
    "autobox",
];
pub const VALUES: &[&str] = &[
    "reference-vina",
    "engine",
    "trials",
    "backends",
    "benchmark-dir",
    "manifest",
    "limit",
    "metal",
    "lanes",
    "steps",
    "local-steps",
    "config",
    "receptor",
    "flex",
    "ligand",
    "batch",
    "scoring",
    "maps",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "out",
    "dir",
    "write_maps",
    "cpu",
    "seed",
    "exhaustiveness",
    "max_evals",
    "num_modes",
    "min_rmsd",
    "energy_range",
    "spacing",
    "verbosity",
    "unbound_energy",
    "weight_gauss1",
    "weight_gauss2",
    "weight_repulsion",
    "weight_hydrophobic",
    "weight_hydrogen",
    "weight_rot",
    "weight_vinardo_gauss1",
    "weight_vinardo_repulsion",
    "weight_vinardo_hydrophobic",
    "weight_vinardo_hydrogen",
    "weight_vinardo_rot",
    "weight_ad4_vdw",
    "weight_ad4_hb",
    "weight_ad4_elec",
    "weight_ad4_dsolv",
    "weight_ad4_rot",
    "weight_glue",
];
pub type Options = BTreeMap<String, Vec<String>>;
fn insert(options: &mut Options, key: &str, values: Vec<String>) -> Result<(), String> {
    if !FLAGS.contains(&key) && !VALUES.contains(&key) {
        return Err(format!("Unknown option --{key}"));
    }
    if FLAGS.contains(&key)
        && (values.len() != 1 || !["true", "false", "1", "0"].contains(&values[0].as_str()))
    {
        return Err(format!("Invalid boolean --{key}"));
    }
    if options.insert(key.into(), values).is_some() {
        return Err(format!("Duplicate option --{key}"));
    }
    Ok(())
}
pub fn parse(args: &[String]) -> Result<Options, String> {
    let mut options = Options::new();
    let mut i = 0;
    while i < args.len() {
        let arg = args[i]
            .strip_prefix("--")
            .or_else(|| if args[i] == "-h" { Some("help") } else { None })
            .ok_or_else(|| format!("Expected an option, got {}", args[i]))?;
        let (key, inline) = arg
            .split_once('=')
            .map_or((arg, None), |(k, v)| (k, Some(v)));
        i += 1;
        let mut values = Vec::new();
        if FLAGS.contains(&key) {
            if let Some(value) = inline {
                if !["true", "false", "1", "0"].contains(&value) {
                    return Err(format!("Invalid boolean --{key}"));
                }
                values.push(value.into());
            } else {
                values.push("true".into());
            }
        } else {
            if let Some(value) = inline {
                values.push(value.into());
            }
            while i < args.len()
                && !args[i].starts_with("--")
                && args[i] != "-h"
                && (values.is_empty() || key == "ligand" || key == "batch")
            {
                values.push(args[i].clone());
                i += 1;
            }
            if values.is_empty() {
                return Err(format!("Missing value for --{key}"));
            }
        }
        insert(&mut options, key, values)?;
    }
    if let Some(path) = options.get("config").and_then(|v| v.first()).cloned() {
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
        let mut config = Options::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap().trim();
            if line.is_empty() {
                continue;
            }
            let (key, value) = line
                .split_once('=')
                .ok_or_else(|| format!("{path}:{}: expected key = value", i + 1))?;
            let key = key.trim();
            let value = value.trim();
            let values = if key == "ligand" || key == "batch" {
                value.split_whitespace().map(String::from).collect()
            } else {
                vec![value.into()]
            };
            insert(&mut config, key, values)?;
        }
        config.extend(options);
        options = config;
    }
    Ok(options)
}
pub fn value<'a>(o: &'a Options, key: &str, default: &'a str) -> &'a str {
    o.get(key)
        .and_then(|v| v.first())
        .map_or(default, String::as_str)
}
pub fn flag(o: &Options, key: &str) -> bool {
    matches!(value(o, key, "false"), "true" | "1")
}
pub fn number<T: std::str::FromStr>(o: &Options, key: &str, default: &str) -> Result<T, String> {
    value(o, key, default)
        .parse()
        .map_err(|_| format!("Invalid numeric value for --{key}"))
}
