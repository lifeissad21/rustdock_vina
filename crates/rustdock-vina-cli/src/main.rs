use rustdock_vina_cli::options::{flag, number, parse, value, Options, VALUES};
use rustdock_vina_core::{common::Vec3, vina::Vina};
use std::path::Path;
mod metal;
use rustdock_vina_cli::output;

fn metal_controls(o: &Options) -> Result<metal::Controls, String> {
    let optional = |key: &str| -> Result<Option<u32>, String> {
        if o.contains_key(key) {
            number(o, key, "").map(Some)
        } else {
            Ok(None)
        }
    };
    let c = metal::Controls {
        lanes: optional("lanes")?,
        steps: optional("steps")?,
        local_steps: optional("local-steps")?,
    };
    c.validate()?;
    if c.is_set() && value(o, "metal", "off") != "on" {
        return Err(
            "Experimental GPU controls require --metal on; they cannot apply to CPU fallback"
                .into(),
        );
    }
    Ok(c)
}
fn weights(o: &Options, scoring: &str) -> Result<Vec<f64>, String> {
    let specs: &[(&str, &str)] = match scoring {
        "vina" => &[
            ("weight_gauss1", "-0.035579"),
            ("weight_gauss2", "-0.005156"),
            ("weight_repulsion", "0.840245"),
            ("weight_hydrophobic", "-0.035069"),
            ("weight_hydrogen", "-0.587439"),
            ("weight_glue", "50"),
            ("weight_rot", "0.05846"),
        ],
        "vinardo" => &[
            ("weight_vinardo_gauss1", "-0.045"),
            ("weight_vinardo_repulsion", "0.8"),
            ("weight_vinardo_hydrophobic", "-0.035"),
            ("weight_vinardo_hydrogen", "-0.600"),
            ("weight_glue", "50"),
            ("weight_vinardo_rot", "0.05846"),
        ],
        _ => &[
            ("weight_ad4_vdw", "0.1662"),
            ("weight_ad4_hb", "0.1209"),
            ("weight_ad4_elec", "0.1406"),
            ("weight_ad4_dsolv", "0.1322"),
            ("weight_glue", "50"),
            ("weight_ad4_rot", "0.2983"),
        ],
    };
    specs
        .iter()
        .map(|(key, default)| number(o, key, default))
        .collect()
}
fn score_output(es: [f64; 8]) {
    for (label, energy) in [
        ("Affinity", es[0]),
        ("Intermolecular energy", es[1] + es[2]),
        ("Internal energy", es[3] + es[4] + es[5]),
        ("Torsional energy", es[6]),
        ("Unbound energy", es[7]),
    ] {
        println!("  {label:<23}: {energy:>10.3} kcal/mol");
    }
}

fn run_case(
    engine: &mut Vina,
    o: &Options,
    ligands: &[String],
    output: &str,
) -> Result<(), String> {
    let started = std::time::Instant::now();
    let verbosity: u32 = number(o, "verbosity", "1")?;
    let display = output::Output::new(verbosity);
    display.heading("RustDock Vina · 0.1.0");
    display.field(
        "Receptor",
        value(o, "receptor", value(o, "maps", "flexible residues")),
    );
    for ligand in ligands {
        display.field("Ligand", ligand);
    }
    display.field("Scoring", value(o, "scoring", "vina"));
    display.field("Seed", engine.seed);
    display.field("CPU threads", engine.cpu);
    if !["score_only", "local_only", "randomize_only"]
        .iter()
        .any(|f| flag(o, f))
    {
        display.field(
            "Acceleration",
            match value(o, "metal", "off") {
                "on" => "Metal GPU",
                "auto" => "Metal GPU / CPU fallback",
                _ => "CPU",
            },
        );
        display.field("Exhaustiveness", value(o, "exhaustiveness", "8"));
        display.field("Requested poses", value(o, "num_modes", "9"));
        for (key, label) in [
            ("lanes", "GPU lanes"),
            ("steps", "GPU steps"),
            ("local-steps", "GPU local steps"),
        ] {
            if o.contains_key(key) {
                display.field(label, value(o, key, ""));
            }
        }
    }
    display.stage("Preparing ligand topology…");
    engine.set_ligands_from_files(ligands)?;
    if !value(o, "maps", "").is_empty() {
        display.stage("Loading affinity maps…");
        engine.load_maps(value(o, "maps", ""))?;
    } else {
        let (center, size) =
            if flag(o, "autobox") || flag(o, "score_only") && !o.contains_key("center_x") {
                engine.grid_dimensions_from_ligand(4.0)?
            } else {
                for key in [
                    "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
                ] {
                    if !o.contains_key(key) {
                        return Err(format!(
                            "Missing --{key}; supply a search box or use --autobox"
                        ));
                    }
                }
                (
                    Vec3::new(
                        number(o, "center_x", "0")?,
                        number(o, "center_y", "0")?,
                        number(o, "center_z", "0")?,
                    ),
                    Vec3::new(
                        number(o, "size_x", "0")?,
                        number(o, "size_y", "0")?,
                        number(o, "size_z", "0")?,
                    ),
                )
            };
        display.field(
            "Box center",
            format!("{:.3}, {:.3}, {:.3} Å", center[0], center[1], center[2]),
        );
        display.field(
            "Box size",
            format!("{:.3} × {:.3} × {:.3} Å", size[0], size[1], size[2]),
        );
        display.stage("Computing affinity maps…");
        engine.compute_vina_maps(
            center,
            size,
            number(o, "spacing", "0.375")?,
            flag(o, "force_even_voxels"),
        )?;
    }
    if o.contains_key("write_maps") {
        engine.write_maps(value(o, "write_maps", ""))?;
    }
    if flag(o, "score_only") {
        let es = if o.contains_key("unbound_energy") {
            engine.score_with_unbound(number(o, "unbound_energy", "0")?)?
        } else {
            engine.score()?
        };
        display.heading("Score summary");
        score_output(es);
    } else if flag(o, "local_only") {
        display.stage("Optimizing ligand pose…");
        let es = engine.optimize(0)?;
        display.heading("Score summary");
        score_output(es);
        engine.write_pose(
            output,
            &format!("REMARK VINA RESULT: {:9.3}     0.000     0.000\n", es[0]),
        )?;
    } else if flag(o, "randomize_only") {
        display.stage("Randomizing ligand pose…");
        engine.randomize(10000)?;
        engine.write_pose(output, "")?;
    } else {
        let modes = number(o, "num_modes", "9")?;
        let exhaustiveness = number(o, "exhaustiveness", "8")?;
        let min_rmsd = number(o, "min_rmsd", "1")?;
        let max_evals = number(o, "max_evals", "0")?;
        let backend = value(o, "metal", "off");
        let mut used_metal = false;
        display.stage("Searching poses and refining candidates…");
        if backend != "off" {
            match metal::dock(
                engine,
                exhaustiveness,
                modes,
                min_rmsd,
                max_evals,
                number(o, "verbosity", "1")?,
                metal_controls(o)?,
            ) {
                Ok(()) => used_metal = true,
                Err(error) if backend == "auto" => {
                    eprintln!("Metal unavailable for this run; using CPU: {error}")
                }
                Err(error) => return Err(error),
            }
        }
        if !used_metal {
            engine.global_search(exhaustiveness, modes, min_rmsd, max_evals)?;
        }
        let range: f64 = number(o, "energy_range", "3")?;
        engine.write_poses(output, modes, range)?;
        let rows: Vec<_> = engine
            .poses
            .iter()
            .take(modes)
            .take_while(|p| p.e <= engine.poses[0].e + range)
            .map(|p| (p.e, p.lb, p.ub))
            .collect();
        display.poses(&rows);
    }
    display.heading("Complete");
    if !flag(o, "score_only") {
        display.field("Saved poses", output);
    }
    if o.contains_key("write_maps") {
        display.field("Saved maps", value(o, "write_maps", ""));
    }
    display.field(
        "Total time",
        format!("{:.2} s", started.elapsed().as_secs_f64()),
    );
    Ok(())
}
fn run(args: &[String]) -> Result<(), String> {
    let o = parse(args)?;
    if [
        "reference-vina",
        "engine",
        "trials",
        "backends",
        "benchmark-dir",
        "manifest",
        "limit",
        "check",
    ]
    .iter()
    .any(|key| o.contains_key(*key))
    {
        return Err("Benchmark options belong to the vina-benchmark executable".into());
    }
    metal_controls(&o)?;
    if !["on", "off", "auto"].contains(&value(&o, "metal", "off")) {
        return Err("--metal must be on, off or auto".into());
    }
    if value(&o, "metal", "off") == "on"
        && ["score_only", "local_only", "randomize_only"]
            .iter()
            .any(|f| flag(&o, f))
    {
        return Err("Metal accelerates docking search; use --metal off for scoring, local optimization or randomization".into());
    }
    if flag(&o, "version") {
        println!("RustDock Vina 0.1.0 (native Rust port)");
        return Ok(());
    }
    if args.is_empty() || flag(&o, "help") || flag(&o, "help_advanced") {
        println!(
            r#"RustDock Vina · native Rust docking

Usage
  vina-benchmark --help       Repeated CPU / Metal comparisons
  vina --receptor receptor.pdbqt --ligand ligand.pdbqt
       --config box.txt --out poses.pdbqt

Inputs
  --receptor FILE             Rigid receptor
  --flex FILE                 Flexible receptor residues
  --ligand FILE...            One or more ligands
  --batch FILE/DIR...         Batch ligands; requires --dir
  --maps PREFIX               Load existing affinity maps
  --config FILE               Read key = value settings

Search
  --scoring vina|vinardo|ad4   Scoring function [vina]
  --metal off|on|auto          CPU, require GPU, or automatic fallback [off]
  --exhaustiveness N          Search effort [8]; higher searches longer
  --num_modes N               Maximum poses to return [9]
  --cpu N                     CPU threads [0 = automatic]
  --seed N                    Random seed [0 = generated]
  --max_evals N               CPU evaluation limit [0 = unlimited]

Experimental GPU controls (require --metal on)
  --lanes N                   Independent GPU searches [automatic; 1–65536]
  --steps N                   Global mutation steps per lane [automatic]
  --local-steps N             BFGS iterations per GPU optimization [automatic]
                              Positive integers; Rust final refinement is unchanged

Search box
  --center_x/y/z N            Box center in Angstrom
  --size_x/y/z N              Box dimensions in Angstrom
  --spacing N                 Grid spacing [0.375]
  --autobox                   Infer box for scoring/local optimization

Actions and output
  --score_only                Score the input pose
  --local_only                Optimize the input pose
  --randomize_only            Randomize the input pose
  --out FILE                  Save predicted poses
  --dir DIR                   Batch output directory
  --min_rmsd N                Minimum pose separation [1 Angstrom]
  --energy_range N            Pose energy window [3 kcal/mol]
  --write_maps PREFIX         Save affinity maps
  --force_even_voxels          Use even grid dimensions
  --no_refine                 Use affinity maps for final scoring
  --verbosity 0|1|2           Output detail [1]
  --help_advanced             Show scoring weights
  --version                   Show version

Terminal colors are automatic; set NO_COLOR to disable them."#
        );
        if flag(&o, "help_advanced") {
            println!(
                "Advanced weights: {}",
                VALUES
                    .iter()
                    .filter(|v| v.starts_with("weight_"))
                    .map(|v| format!("--{v}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        return Ok(());
    }
    if ["score_only", "local_only", "randomize_only"]
        .iter()
        .filter(|f| flag(&o, f))
        .count()
        > 1
    {
        return Err("Choose only one action".into());
    }
    if flag(&o, "autobox") && !flag(&o, "score_only") && !flag(&o, "local_only") {
        return Err("--autobox is only supported for scoring and local optimization".into());
    }
    let verbosity: u32 = number(&o, "verbosity", "1")?;
    if verbosity > 2 {
        return Err("--verbosity must be 0, 1 or 2".into());
    }
    let scoring = value(&o, "scoring", "vina");
    let mut engine = Vina::new(
        scoring,
        number(&o, "cpu", "0")?,
        number(&o, "seed", "0")?,
        flag(&o, "no_refine"),
    )?;
    engine.set_weights(weights(&o, scoring)?)?;
    let rigid = value(&o, "receptor", "");
    let flex = value(&o, "flex", "");
    if !rigid.is_empty() || !flex.is_empty() {
        engine.set_receptor(rigid, flex)?;
    }
    if o.contains_key("batch") {
        if o.contains_key("ligand") || o.contains_key("out") {
            return Err("--batch cannot be combined with --ligand or --out".into());
        }
        let dir = Path::new(value(&o, "dir", ""));
        if value(&o, "dir", "").is_empty() || !dir.is_dir() {
            return Err("Batch mode requires an existing --dir output directory".into());
        }
        let mut paths = vec![];
        for path in &o["batch"] {
            if Path::new(path).is_dir() {
                for entry in std::fs::read_dir(path).map_err(|e| e.to_string())? {
                    let path = entry.map_err(|e| e.to_string())?.path();
                    if path.extension().is_some_and(|s| s == "pdbqt") {
                        paths.push(path.to_string_lossy().into_owned());
                    }
                }
            } else {
                paths.push(path.clone());
            }
        }
        paths.sort();
        if paths.is_empty() {
            return Err("No PDBQT ligands in batch".into());
        }
        let mut names = std::collections::HashSet::new();
        for path in &paths {
            let name = Path::new(path)
                .file_stem()
                .ok_or("Invalid ligand path")?
                .to_string_lossy();
            if !names.insert(name.into_owned()) {
                return Err("Batch ligand basenames collide".into());
            }
        }
        for path in paths {
            engine = Vina::new(scoring, engine.cpu, engine.seed, flag(&o, "no_refine"))?;
            engine.set_weights(weights(&o, scoring)?)?;
            if !rigid.is_empty() || !flex.is_empty() {
                engine.set_receptor(rigid, flex)?;
            }
            let name = Path::new(&path).file_stem().unwrap().to_string_lossy();
            let out = dir.join(format!("{name}_out.pdbqt"));
            run_case(&mut engine, &o, &[path], &out.to_string_lossy())?;
        }
    } else {
        let ligands = o.get("ligand").ok_or("Specify --ligand or --batch")?;
        let default = Path::new(&ligands[0]).with_file_name(format!(
            "{}_out.pdbqt",
            Path::new(&ligands[0])
                .file_stem()
                .ok_or("Invalid ligand path")?
                .to_string_lossy()
        ));
        let output = value(&o, "out", default.to_str().ok_or("Non-UTF8 output path")?);
        if ligands.iter().any(|p| Path::new(p) == Path::new(output))
            || Path::new(rigid) == Path::new(output)
            || !flex.is_empty() && Path::new(flex) == Path::new(output)
        {
            return Err("Output must not overwrite an input file".into());
        }
        run_case(&mut engine, &o, ligands, output)?;
    }
    Ok(())
}
fn main() {
    if let Err(error) = run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("Vina error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn metal_flag_validates_modes_and_actions() {
        for mode in ["on", "off", "auto"] {
            let o = parse(&["--metal".into(), mode.into()]).unwrap();
            assert_eq!(value(&o, "metal", "off"), mode);
        }
        assert!(run(&["--metal".into(), "yes".into()])
            .unwrap_err()
            .contains("on, off or auto"));
        assert!(run(&["--metal".into(), "on".into(), "--score_only".into()])
            .unwrap_err()
            .contains("accelerates docking search"));
    }
    #[test]
    fn experimental_controls_parse_and_reject_cpu_or_invalid_values() {
        let args = [
            "--metal",
            "on",
            "--lanes",
            "1024",
            "--steps",
            "100",
            "--local-steps",
            "10",
        ]
        .map(String::from);
        let c = metal_controls(&parse(&args).unwrap()).unwrap();
        assert_eq!(
            (c.lanes, c.steps, c.local_steps),
            (Some(1024), Some(100), Some(10))
        );
        for args in [
            vec!["--lanes", "100"],
            vec!["--metal", "auto", "--steps", "10"],
            vec!["--metal", "on", "--lanes", "0"],
            vec!["--metal", "on", "--steps", "-1"],
            vec!["--metal", "on", "--local-steps", "0"],
        ] {
            assert!(run(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>()).is_err());
        }
    }
    #[test]
    fn parses_multiple_ligands_and_negative_coordinates() {
        let args = [
            "--ligand",
            "one.pdbqt",
            "two.pdbqt",
            "--center_x",
            "-10.5",
            "--score_only",
        ]
        .map(String::from);
        let o = parse(&args).unwrap();
        assert_eq!(o["ligand"].len(), 2);
        assert_eq!(number::<f64>(&o, "center_x", "0").unwrap(), -10.5);
        assert!(flag(&o, "score_only"));
        assert!(parse(&["--unknown".into(), "x".into()]).is_err());
        assert!(parse(&["--receptor".into()]).is_err());
    }
}
