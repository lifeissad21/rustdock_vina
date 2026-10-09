use rustdock_vina_cli::options::{flag, number, parse, value, Options, FLAGS};
use rustdock_vina_cli::output::Output;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path, process::Command, time::Instant};

const KEYS: &[&str] = &[
    "reference-vina",
    "engine",
    "trials",
    "backends",
    "benchmark-dir",
    "manifest",
    "limit",
    "check",
];

fn modes(o: &Options) -> Result<Vec<String>, String> {
    let mut modes = Vec::new();
    for mode in value(o, "backends", value(o, "metal", "vina,off,on")).split(',') {
        let mode = if mode == "cpu" { "off" } else { mode };
        if !["vina", "off", "on", "auto"].contains(&mode) || modes.iter().any(|m| m == mode) {
            return Err(
                "--backends must contain unique vina, off (or cpu), on, auto modes, separated by commas"
                    .into(),
            );
        }
        modes.push(mode.to_owned());
    }
    if modes.contains(&"auto".into())
        && ["lanes", "steps", "local-steps"]
            .iter()
            .any(|k| o.contains_key(*k))
    {
        return Err("GPU overrides require on and cannot be used with auto".into());
    }
    Ok(modes)
}

fn child_args(o: &Options, mode: &str, seed: i32, out: &Path) -> Vec<String> {
    let mut args = Vec::new();
    for (key, values) in o {
        if KEYS.contains(&key.as_str())
            || ["config", "metal", "seed", "out"].contains(&key.as_str())
            || mode != "on" && ["lanes", "steps", "local-steps"].contains(&key.as_str())
        {
            continue;
        }
        if FLAGS.contains(&key.as_str()) {
            if mode == "vina" {
                if flag(o, key) {
                    args.push(format!("--{key}"));
                }
            } else {
                args.push(format!("--{key}={}", values[0]));
            }
        } else {
            args.push(format!("--{key}"));
            args.extend(values.clone());
        }
    }
    if mode != "vina" {
        args.extend(["--metal".into(), mode.into()]);
    }
    args.extend([
        "--seed".into(),
        seed.to_string(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
    ]);
    args
}

fn pose(text: &str) -> Result<(f64, BTreeMap<u32, [f64; 3]>), String> {
    let mut score = None;
    let mut atoms = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with("ENDMDL") {
            break;
        }
        if let Some(rest) = line.strip_prefix("REMARK VINA RESULT:") {
            score = rest
                .split_whitespace()
                .next()
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|s| s.is_finite());
        }
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            if ["H", "HD", "HS"].contains(&line.split_whitespace().last().unwrap_or("")) {
                continue;
            }
            let field = |a, b| line.get(a..b).ok_or("Malformed PDBQT atom");
            let serial = field(6, 11)?
                .trim()
                .parse::<u32>()
                .map_err(|_| "Invalid atom serial")?;
            let mut xyz = [0.0; 3];
            for (i, coordinate) in xyz.iter_mut().enumerate() {
                *coordinate = field(30 + 8 * i, 38 + 8 * i)?
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| "Invalid atom coordinate")?;
                if !coordinate.is_finite() {
                    return Err("Nonfinite atom coordinate".into());
                }
            }
            if atoms.insert(serial, xyz).is_some() {
                return Err("Duplicate atom serial".into());
            }
        }
    }
    if atoms.is_empty() {
        return Err("No heavy atoms in docking output".into());
    }
    Ok((
        score.ok_or("Missing finite affinity in docking output")?,
        atoms,
    ))
}
fn rmsd(a: &BTreeMap<u32, [f64; 3]>, b: &BTreeMap<u32, [f64; 3]>) -> Result<f64, String> {
    if a.is_empty() || a.keys().ne(b.keys()) {
        return Err("Pose heavy atom serials do not match".into());
    }
    Ok((a
        .iter()
        .map(|(serial, xyz)| (0..3).map(|i| (xyz[i] - b[serial][i]).powi(2)).sum::<f64>())
        .sum::<f64>()
        / a.len() as f64)
        .sqrt())
}

fn backend_label(mode: &str) -> &str {
    match mode {
        "vina" => "Official Vina",
        "off" => "Rust CPU",
        "on" => "Rust Metal",
        _ => "Auto",
    }
}
fn metric(value: &Value) -> String {
    value.as_f64().map_or("—".into(), |n| format!("{n:.3}"))
}

pub fn run(args: &[String]) -> Result<(), String> {
    let mut o = parse(args)?;
    if args.is_empty() || flag(&o, "help") {
        println!("Usage: vina-benchmark --receptor FILE --ligand FILE --config FILE [normal docking options]\n  --reference-vina FILE Official AutoDock Vina executable\n                       [../vina-multicore-benchmark/bin/vina]\n  --engine FILE        Docking binary [vina beside this executable]\n  --trials N           Trials per backend [3]; paired seed = seed + trial index\n  --backends MODES     Comma-separated vina, off/cpu, on, auto [vina,off,on]\n  --metal off|on|auto   Select one backend instead of --backends\n  --benchmark-dir DIR  New output directory [benchmark-results]\n  --seed N             First seed [20260717]; must be positive\n  --manifest FILE      Run each pair in a benchmark manifest\n  --limit N            First N manifest cases [all]\n  --check              Validate inputs without docking\nGPU controls apply only to on. CPU and metal off are the same backend.\nOfficial comparisons with a receptor generate maps per run for all backends.\nEach run uses a fresh process. Outputs: poses, logs, results.csv and summary.json.\nNo fallback for on; auto records the actual backend. Existing output directories are never overwritten.");
        return Ok(());
    }
    let trials: u32 = number(&o, "trials", "3")?;
    let seed: i32 = number(&o, "seed", "20260717")?;
    if trials == 0 || seed <= 0 || i64::from(seed) + i64::from(trials) - 1 > i64::from(i32::MAX) {
        return Err("Trials must be positive and seed + trials - 1 must fit a positive 32-bit signed integer".into());
    }
    if o.contains_key("metal") && o.contains_key("backends") {
        return Err("Choose --metal or --backends".into());
    }
    let backends = modes(&o)?;
    let baseline_mode = if backends.iter().any(|m| m == "vina") {
        "vina"
    } else {
        "off"
    };
    let reference = if baseline_mode == "vina" {
        let path = Path::new(value(
            &o,
            "reference-vina",
            "../vina-multicore-benchmark/bin/vina",
        ));
        if !path.is_file() {
            return Err("Official Vina is required: supply --reference-vina /path/to/official/vina, or choose --backends off,on for Rust only".into());
        }
        let path = fs::canonicalize(path).map_err(|e| e.to_string())?;
        let result = Command::new(&path)
            .arg("--version")
            .output()
            .map_err(|e| e.to_string())?;
        let version = format!(
            "{}{}",
            String::from_utf8_lossy(&result.stdout),
            String::from_utf8_lossy(&result.stderr)
        )
        .trim()
        .to_string();
        if !result.status.success()
            || version.contains("RustDock")
            || !version.to_lowercase().contains("vina")
        {
            return Err("--reference-vina must be an official AutoDock Vina executable with working --version".into());
        }
        Some((path, version))
    } else {
        None
    };
    let generated_maps = reference.is_some() && o.contains_key("receptor");
    if generated_maps {
        o.remove("maps");
    }

    for key in [
        "batch",
        "dir",
        "out",
        "write_maps",
        "score_only",
        "local_only",
        "randomize_only",
        "version",
    ] {
        if o.contains_key(key) {
            return Err(format!("--{key} is not supported in a docking benchmark"));
        }
    }
    for (key, max) in [
        ("lanes", 65536),
        ("steps", u32::MAX),
        ("local-steps", u32::MAX),
    ] {
        if o.contains_key(key) {
            let n: u32 = number(&o, key, "")?;
            if n == 0 || n > max {
                return Err(format!("--{key} must be between 1 and {max}"));
            }
        }
    }
    if !backends.contains(&"on".into())
        && ["lanes", "steps", "local-steps"]
            .iter()
            .any(|k| o.contains_key(*k))
    {
        return Err("GPU controls require an on backend".into());
    }
    let display = Output::new(1);
    let dir = Path::new(value(&o, "benchmark-dir", "benchmark-results"));
    if o.contains_key("manifest") {
        let manifest: Value = serde_json::from_str(
            &fs::read_to_string(value(&o, "manifest", "")).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        let cases = manifest["cases"]
            .as_array()
            .ok_or("Manifest requires cases array")?;
        let limit: usize = number(&o, "limit", &cases.len().to_string())?;
        if limit == 0 || limit > cases.len() {
            return Err("Invalid manifest limit".into());
        }
        if !flag(&o, "check") {
            fs::create_dir(dir)
                .map_err(|e| format!("{}: {e}; choose a new benchmark directory", dir.display()))?;
        }
        display.heading("RustDock Vina · dataset benchmark");
        display.field("Cases", limit);
        display.field("Results", dir.display());
        for (i, case) in cases.iter().take(limit).enumerate() {
            display.stage(&format!(
                "Case {} / {limit} · {} · ligand {}",
                i + 1,
                case["target"].as_str().unwrap_or("unknown"),
                case["ligand"].as_str().unwrap_or("unknown")
            ));
            let mut pair = o.clone();
            pair.remove("manifest");
            pair.remove("limit");
            for (key, source) in [
                ("receptor", "receptor"),
                ("ligand", "pdbqt"),
                ("config", "config"),
            ] {
                pair.insert(
                    key.into(),
                    vec![case[source]
                        .as_str()
                        .ok_or("Missing manifest input")?
                        .into()],
                );
            }
            pair.insert(
                "maps".into(),
                vec![format!(
                    "{}/affinity",
                    case["maps"].as_str().ok_or("Missing manifest maps")?
                )],
            );
            pair.insert(
                "benchmark-dir".into(),
                vec![dir
                    .join(format!("case-{:03}", i + 1))
                    .to_string_lossy()
                    .into_owned()],
            );
            pair.insert(
                "seed".into(),
                vec![seed
                    .checked_add(i as i32)
                    .ok_or("Manifest seed overflow")?
                    .to_string()],
            );
            let mut pair_args = Vec::new();
            for (key, values) in pair {
                if FLAGS.contains(&key.as_str()) {
                    pair_args.push(format!("--{key}={}", values[0]));
                } else {
                    pair_args.push(format!("--{key}"));
                    pair_args.extend(values);
                }
            }
            run(&pair_args)?;
        }
        return Ok(());
    }
    if o.contains_key("limit") {
        return Err("--limit requires --manifest".into());
    }
    if o.get("ligand").is_none_or(|v| v.len() != 1) {
        return Err("Benchmark requires exactly one --ligand per receptor/ligand pair".into());
    }
    if !o.contains_key("receptor") && !o.contains_key("maps") {
        return Err("Specify --receptor or --maps".into());
    }
    for key in ["ligand", "receptor", "flex"] {
        if let Some(paths) = o.get(key) {
            for path in paths {
                if !Path::new(path).is_file() {
                    return Err(format!("Missing input: {path}"));
                }
            }
        }
    }
    if let Some(prefix) = o.get("maps").and_then(|v| v.first()) {
        let prefix = Path::new(prefix);
        let parent = prefix
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let name = prefix
            .file_name()
            .ok_or("Invalid maps prefix")?
            .to_string_lossy();
        if !fs::read_dir(parent)
            .map_err(|e| e.to_string())?
            .filter_map(Result::ok)
            .any(|e| {
                let file = e.file_name().to_string_lossy().into_owned();
                file.starts_with(&format!("{name}.")) && file.ends_with(".map")
            })
        {
            return Err(format!("No affinity maps for {}", prefix.display()));
        }
    } else if ![
        "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
    ]
    .iter()
    .all(|k| o.contains_key(*k))
    {
        return Err("Supply a search box via --config or center/size options, or --maps".into());
    }
    if flag(&o, "check") {
        println!("Inputs checked: {}", value(&o, "ligand", ""));
        return Ok(());
    }
    let exe = if o.contains_key("engine") {
        std::path::PathBuf::from(value(&o, "engine", ""))
    } else {
        std::env::current_exe()
            .map_err(|e| e.to_string())?
            .with_file_name(if cfg!(windows) { "vina.exe" } else { "vina" })
    };
    if !exe.is_file() {
        return Err(format!(
            "Docking binary not found: {}; build vina or specify --engine",
            exe.display()
        ));
    }
    if let Some((path, _)) = &reference {
        if fs::canonicalize(&exe).map_err(|e| e.to_string())? == *path {
            return Err("Reference and Rust engine must be different executables".into());
        }
    }
    fs::create_dir(dir)
        .map_err(|e| format!("{}: {e}; choose a new benchmark directory", dir.display()))?;
    let mut records = Vec::new();
    let mut csv = String::from("trial,requested_backend,actual_backend,seed,seconds,best_affinity,baseline_backend,speedup_vs_baseline,score_delta,rmsd\n");
    display.heading("RustDock Vina · docking benchmark");
    display.field("Receptor", value(&o, "receptor", value(&o, "maps", "")));
    display.field("Ligand", value(&o, "ligand", ""));
    display.field("Trials / backend", trials);
    display.field(
        "Backends",
        backends
            .iter()
            .map(|m| backend_label(m))
            .collect::<Vec<_>>()
            .join(", "),
    );
    display.field("Baseline", backend_label(baseline_mode));
    if let Some((path, version)) = &reference {
        display.field("Official Vina", version);
        display.field("Reference binary", path.display());
    }
    if generated_maps {
        display.field("Maps", "Generated per run from the same search box");
    }
    display.field("First seed", seed);
    display.field("CPU threads", value(&o, "cpu", "automatic"));
    display.field("Exhaustiveness", value(&o, "exhaustiveness", "8"));
    display.field("Requested poses", value(&o, "num_modes", "9"));
    for (key, label) in [
        ("lanes", "GPU lanes"),
        ("steps", "GPU steps"),
        ("local-steps", "GPU local steps"),
    ] {
        if o.contains_key(key) {
            display.field(label, value(&o, key, ""));
        }
    }
    display.field("Results", dir.display());
    for trial in 0..trials {
        let seed = seed + trial as i32;
        let mut baseline = None;
        // Run off first so every other backend can be compared within this trial.
        let mut ordered = backends.clone();
        ordered.sort_by_key(|m| m != baseline_mode);
        for mode in ordered {
            let stem = format!("trial-{:03}-{mode}", trial + 1);
            let out = dir.join(format!("{stem}.pdbqt"));
            let arguments = child_args(&o, &mode, seed, &out);
            display.stage(&format!(
                "Trial {} / {trials} · {} · seed {seed} · searching…",
                trial + 1,
                backend_label(&mode)
            ));
            let start = Instant::now();
            let binary = if mode == "vina" {
                &reference.as_ref().unwrap().0
            } else {
                &exe
            };
            let result = Command::new(binary)
                .args(&arguments)
                .env("NO_COLOR", "1")
                .output()
                .map_err(|e| e.to_string())?;
            let seconds = start.elapsed().as_secs_f64();
            let stdout = String::from_utf8_lossy(&result.stdout);
            let stderr = String::from_utf8_lossy(&result.stderr);
            let log = dir.join(format!("{stem}.log"));
            fs::write(&log, format!("{stdout}\n{stderr}")).map_err(|e| e.to_string())?;
            if !result.status.success() {
                return Err(format!(
                    "{mode} trial {} failed ({}); see {}",
                    trial + 1,
                    result.status,
                    log.display()
                ));
            }
            let actual = if mode == "vina" {
                "vina"
            } else if mode == "off" || mode == "auto" && stderr.contains("using CPU:") {
                "off"
            } else {
                "on"
            };
            let (score, atoms) = pose(&fs::read_to_string(&out).map_err(|e| e.to_string())?)?;
            let comparison = baseline
                .as_ref()
                .map(|(s, a, time)| rmsd(a, &atoms).map(|r| (score - s, r, time / seconds)))
                .transpose()?;
            if mode == baseline_mode {
                baseline = Some((score, atoms, seconds));
            }
            let (delta, rmsd, speedup) =
                comparison.map_or((None, None, None), |(d, r, s)| (Some(d), Some(r), Some(s)));
            display.field(
                "Completed",
                format!(
                    "{seconds:.3} s · {score:.3} kcal/mol · {}",
                    backend_label(actual)
                ),
            );
            csv.push_str(&format!(
                "{},{mode},{actual},{seed},{seconds:.6},{score:.3},{baseline_mode},{},{},{}\n",
                trial + 1,
                speedup.map_or(String::new(), |n| n.to_string()),
                delta.map_or(String::new(), |n| n.to_string()),
                rmsd.map_or(String::new(), |n| n.to_string())
            ));
            records.push(json!({"trial":trial+1,"requested_backend":mode,"actual_backend":actual,"seed":seed,"seconds":seconds,"best_affinity":score,"baseline_backend":baseline_mode,"speedup_vs_baseline":speedup,"score_delta":delta,"rmsd":rmsd,"score_delta_vs_off":if baseline_mode == "off" {delta} else {None},"rmsd_vs_off":if baseline_mode == "off" {rmsd} else {None},"binary":binary,"command_args":arguments}));
            fs::write(dir.join("results.csv"), &csv).map_err(|e| e.to_string())?;
            fs::write(
                dir.join("summary.json"),
                serde_json::to_string_pretty(&json!({"settings":o,"reference":reference.as_ref().map(|(path,version)|json!({"path":path,"version":version})),"baseline_backend":baseline_mode,"runs":records}))
                    .map_err(|e| e.to_string())?,
            )
            .map_err(|e| e.to_string())?;
        }
    }
    display.heading("Trial results");
    let rows: Vec<_> = records
        .iter()
        .map(|r| {
            vec![
                r["trial"].to_string(),
                backend_label(r["requested_backend"].as_str().unwrap()).into(),
                backend_label(r["actual_backend"].as_str().unwrap()).into(),
                metric(&r["seconds"]),
                metric(&r["best_affinity"]),
                metric(&r["score_delta"]),
                metric(&r["rmsd"]),
            ]
        })
        .collect();
    display.table(
        &[
            "Trial",
            "Requested",
            "Actual",
            "Time (s)",
            "Affinity (kcal/mol)",
            "Δ vs baseline (kcal/mol)",
            "RMSD (Å)",
        ],
        &rows,
    );
    display.field(
        "Comparisons",
        &format!(
            "Best poses versus paired {}; RMSD without fitting",
            backend_label(baseline_mode)
        ),
    );
    display.heading("Backend averages");
    let mut averages = Vec::new();
    let mut baseline_mean = None;
    for mode in ["vina", "off", "on", "auto"] {
        let runs: Vec<_> = records
            .iter()
            .filter(|r| r["requested_backend"] == mode)
            .collect();
        if runs.is_empty() {
            continue;
        }
        let mean = runs
            .iter()
            .map(|r| r["seconds"].as_f64().unwrap())
            .sum::<f64>()
            / runs.len() as f64;
        let affinity = runs
            .iter()
            .map(|r| r["best_affinity"].as_f64().unwrap())
            .sum::<f64>()
            / runs.len() as f64;
        if mode == baseline_mode {
            baseline_mean = Some(mean);
        }
        let speedup = baseline_mean.map_or("—".into(), |cpu| format!("{:.2}×", cpu / mean));
        averages.push(vec![
            backend_label(mode).into(),
            runs.len().to_string(),
            format!("{mean:.3}"),
            format!("{affinity:.3}"),
            speedup,
        ]);
    }
    display.table(
        &[
            "Backend",
            "Runs",
            "Mean time (s)",
            "Mean affinity (kcal/mol)",
            "Baseline speedup",
        ],
        &averages,
    );
    display.field(
        "Timing includes",
        "Startup, maps, search and backend refinement",
    );
    display.heading("Complete");
    display.field("Measurements", dir.join("results.csv").display());
    display.field("Run details", dir.join("summary.json").display());
    display.field("Poses and logs", dir.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn backend_selection_and_forwarding() {
        let o = parse(
            &[
                "--backends",
                "cpu,on",
                "--lanes",
                "37",
                "--no_refine",
                "--trials",
                "2",
            ]
            .map(String::from),
        )
        .unwrap();
        assert_eq!(modes(&o).unwrap(), vec!["off", "on"]);
        let cpu = child_args(&o, "off", 42, Path::new("pose.pdbqt"));
        assert!(!cpu.contains(&"--lanes".into()));
        assert!(!cpu.contains(&"--trials".into()));
        assert!(cpu.contains(&"--no_refine=true".into()));
        assert!(child_args(&o, "on", 42, Path::new("pose.pdbqt")).contains(&"--lanes".into()));
        for modes_text in ["off,cpu", "invalid", "", "on,on"] {
            let mut o = Options::new();
            o.insert("backends".into(), vec![modes_text.into()]);
            assert!(modes(&o).is_err());
        }
    }
    #[test]
    fn official_arguments_omit_gpu_options_and_use_boost_boolean_switches() {
        let o = parse(
            &[
                "--reference-vina",
                "official",
                "--backends",
                "vina,on",
                "--lanes",
                "37",
                "--no_refine",
                "--force_even_voxels=false",
            ]
            .map(String::from),
        )
        .unwrap();
        let args = child_args(&o, "vina", 42, Path::new("pose.pdbqt"));
        for forbidden in [
            "--metal",
            "--lanes",
            "--reference-vina",
            "--no_refine=true",
            "--force_even_voxels",
        ] {
            assert!(!args.contains(&forbidden.into()));
        }
        assert!(args.contains(&"--no_refine".into()));
        assert!(args.contains(&"--seed".into()));
        assert!(args.contains(&"42".into()));
    }

    #[test]
    fn compares_first_pose_heavy_atoms_and_validates_results() {
        let atom = |x: f64| {
            format!(
                "ATOM      1  C   LIG A   1    {x:8.3}{:8.3}{:8.3}  1.00  0.00     0.000 C\n",
                0., 0.
            )
        };
        let text = format!(
            "REMARK VINA RESULT: -7.123 0 0\n{}ENDMDL\n{}",
            atom(0.),
            atom(9.)
        );
        let (score, a) = pose(&text).unwrap();
        assert_eq!(score, -7.123);
        assert_eq!(a.len(), 1);
        let (_, b) = pose(&format!("REMARK VINA RESULT: -6 0 0\n{}", atom(3.))).unwrap();
        assert_eq!(rmsd(&a, &b).unwrap(), 3.);
        assert!(pose("REMARK VINA RESULT: NaN 0 0").is_err());
        assert!(rmsd(&a, &BTreeMap::new()).is_err());
    }
    #[cfg(unix)]
    #[test]
    fn executes_paired_trials_records_fallback_and_preserves_existing_outputs() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!(
            "vina-benchmark-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&dir).unwrap();
        let engine = dir.join("fake-vina");
        fs::write(&engine, r#"#!/bin/sh
while [ "$#" -gt 0 ]; do
  case "$1" in
    --metal) mode="$2"; shift ;;
    --out) out="$2"; shift ;;
  esac
  shift
done
score=-6.000
x=0.000
if [ "$mode" = on ]; then score=-7.000; x=3.000; fi
if [ "$mode" = auto ]; then echo 'Metal unavailable for this run; using CPU: synthetic fallback' >&2; fi
printf 'REMARK VINA RESULT: %s 0 0\nATOM      1  C   LIG A   1    %8s%8s%8s  1.00  0.00     0.000 C\nENDMDL\n' "$score" "$x" 0.000 0.000 > "$out"
"#).unwrap();
        fs::set_permissions(&engine, fs::Permissions::from_mode(0o755)).unwrap();
        let ligand = dir.join("ligand.pdbqt");
        fs::write(&ligand, "").unwrap();
        let receptor = dir.join("receptor.pdbqt");
        fs::write(&receptor, "").unwrap();
        let config = dir.join("box.txt");
        fs::write(
            &config,
            "center_x = 0\ncenter_y = 0\ncenter_z = 0\nsize_x = 10\nsize_y = 10\nsize_z = 10\n",
        )
        .unwrap();
        let output = dir.join("results");
        let args = vec![
            "--engine".into(),
            engine.to_string_lossy().into_owned(),
            "--ligand".into(),
            ligand.to_string_lossy().into_owned(),
            "--receptor".into(),
            receptor.to_string_lossy().into_owned(),
            "--config".into(),
            config.to_string_lossy().into_owned(),
            "--benchmark-dir".into(),
            output.to_string_lossy().into_owned(),
            "--trials".into(),
            "2".into(),
            "--backends".into(),
            "off,on,auto".into(),
            "--seed".into(),
            "42".into(),
        ];
        run(&args).unwrap();
        let summary: Value =
            serde_json::from_str(&fs::read_to_string(output.join("summary.json")).unwrap())
                .unwrap();
        let runs = summary["runs"].as_array().unwrap();
        assert_eq!(runs.len(), 6);
        for trial in 0..2 {
            for record in &runs[trial * 3..trial * 3 + 3] {
                assert_eq!(record["seed"], 42 + trial);
            }
            let metal = &runs[trial * 3 + 1];
            assert_eq!(metal["score_delta_vs_off"], -1.0);
            assert_eq!(metal["rmsd_vs_off"], 3.0);
            assert_eq!(runs[trial * 3 + 2]["actual_backend"], "off");
        }
        assert_eq!(
            fs::read_to_string(output.join("results.csv"))
                .unwrap()
                .lines()
                .count(),
            7
        );
        assert!(run(&args).unwrap_err().contains("new benchmark directory"));
        let reference = dir.join("official-vina");
        let script = fs::read_to_string(&engine).unwrap().replacen(
            "#!/bin/sh",
            "#!/bin/sh\nif [ \"$1\" = --version ]; then echo 'AutoDock Vina v1.2.7'; exit 0; fi",
            1,
        );
        fs::write(&reference, script).unwrap();
        fs::set_permissions(&reference, fs::Permissions::from_mode(0o755)).unwrap();
        let mut official_args = args.clone();
        let backends = official_args
            .iter()
            .position(|v| v == "--backends")
            .unwrap();
        official_args[backends + 1] = "vina,off,on".into();
        let output_index = official_args
            .iter()
            .position(|v| v == "--benchmark-dir")
            .unwrap();
        let official_output = dir.join("official-results");
        official_args[output_index + 1] = official_output.to_string_lossy().into_owned();
        official_args.extend([
            "--reference-vina".into(),
            reference.to_string_lossy().into_owned(),
            "--maps".into(),
            "unused-maps-prefix".into(),
        ]);
        run(&official_args).unwrap();
        let summary: Value = serde_json::from_str(
            &fs::read_to_string(official_output.join("summary.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(summary["baseline_backend"], "vina");
        assert_eq!(summary["reference"]["version"], "AutoDock Vina v1.2.7");
        let runs = summary["runs"].as_array().unwrap();
        assert_eq!(runs.len(), 6);
        for trial in 0..2 {
            let reference = &runs[trial * 3];
            assert_eq!(reference["actual_backend"], "vina");
            assert!(!reference["command_args"]
                .as_array()
                .unwrap()
                .iter()
                .any(|v| v == "--metal"));
            for record in &runs[trial * 3..trial * 3 + 3] {
                assert_eq!(record["seed"], 42 + trial);
                assert!(!record["command_args"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|v| v == "--maps"));
            }
            assert_eq!(runs[trial * 3 + 2]["score_delta"], -1.0);
            assert_eq!(runs[trial * 3 + 2]["rmsd"], 3.0);
            assert!(runs[trial * 3 + 2]["speedup_vs_baseline"].as_f64().unwrap() > 0.);
            assert!(runs[trial * 3 + 2]["score_delta_vs_off"].is_null());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn rejects_invalid_settings_before_running() {
        for args in [
            vec!["--trials", "0"],
            vec!["--seed", "2147483647", "--trials", "2"],
            vec!["--out", "input.pdbqt"],
            vec!["--metal", "on", "--backends", "off,on"],
        ] {
            assert!(run(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
    }
}

fn main() {
    if let Err(error) = run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("Benchmark error: {error}");
        std::process::exit(1);
    }
}
