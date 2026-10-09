//! Stateful local Python bridge. stdout carries only newline-delimited JSON.
use rustdock_vina_core::{common::Vec3, vina::Vina};
use serde_json::{json, Value};
use std::io::{self, BufRead, Write};

fn argument<T: serde_json::value::Index + Copy>(args: &Value, index: T) -> Result<&Value, String> {
    args.get(index).ok_or("Missing argument".into())
}
fn string(args: &Value, index: usize) -> Result<&str, String> {
    argument(args, index)?
        .as_str()
        .ok_or("Expected string".into())
}
fn uint(args: &Value, index: usize) -> Result<usize, String> {
    usize::try_from(
        argument(args, index)?
            .as_u64()
            .ok_or("Expected nonnegative integer")?,
    )
    .map_err(|e| e.to_string())
}
fn float(args: &Value, index: usize) -> Result<f64, String> {
    argument(args, index)?
        .as_f64()
        .filter(|v| v.is_finite())
        .ok_or("Expected finite number".into())
}
fn boolean(args: &Value, index: usize) -> Result<bool, String> {
    argument(args, index)?
        .as_bool()
        .ok_or("Expected boolean".into())
}
fn strings(value: &Value) -> Result<Vec<String>, String> {
    if let Some(s) = value.as_str() {
        return Ok(vec![s.into()]);
    }
    value
        .as_array()
        .ok_or("Expected string or array")?
        .iter()
        .map(|v| v.as_str().map(String::from).ok_or("Expected string".into()))
        .collect()
}
fn dispatch(engine: &mut Option<Vina>, request: Value) -> Result<Value, String> {
    let method = request["method"].as_str().ok_or("Missing method")?;
    let args = &request["args"];
    if method == "new" {
        let vina = Vina::new(
            string(args, 0)?,
            uint(args, 1)?,
            uint(args, 2)? as u64,
            boolean(args, 4)?,
        )?;
        let seed = vina.seed;
        *engine = Some(vina);
        return Ok(json!(seed));
    }
    let v = engine.as_mut().ok_or("Engine has not been initialized")?;
    match method {
        "seed" => return Ok(json!(v.seed)),
        "set_receptor" => v.set_receptor(
            string(args, 0)?,
            args.get(1).and_then(Value::as_str).unwrap_or(""),
        )?,
        "set_ligand_from_file" => v.set_ligands_from_files(&strings(argument(args, 0)?)?)?,
        "set_ligand_from_string" => v.set_ligands_from_strings(&strings(argument(args, 0)?)?)?,
        "set_vina_weights" | "set_vinardo_weights" | "set_ad4_weights" => v.set_weights(
            args.as_array()
                .ok_or("Expected weights array")?
                .iter()
                .map(|n| {
                    n.as_f64()
                        .filter(|x| x.is_finite())
                        .ok_or("Expected finite weight".into())
                })
                .collect::<Result<Vec<_>, String>>()?,
        )?,
        "compute_vina_maps" => v.compute_vina_maps(
            Vec3::new(float(args, 0)?, float(args, 1)?, float(args, 2)?),
            Vec3::new(float(args, 3)?, float(args, 4)?, float(args, 5)?),
            float(args, 6)?,
            boolean(args, 7)?,
        )?,
        "load_maps" => v.load_maps(string(args, 0)?)?,
        "write_maps" => v.write_maps_with_metadata(
            string(args, 0)?,
            string(args, 1)?,
            string(args, 2)?,
            string(args, 3)?,
        )?,
        "write_pose" => v.write_pose(string(args, 0)?, string(args, 1)?)?,
        "write_poses" => v.write_poses(string(args, 0)?, uint(args, 1)?, float(args, 2)?)?,
        "get_poses" => return Ok(json!(v.get_poses(uint(args, 0)?, float(args, 1)?)?)),
        "get_poses_coordinates" => {
            return Ok(json!(
                v.get_poses_coordinates(uint(args, 0)?, float(args, 1)?)?
            ))
        }
        "get_poses_energies" => {
            return Ok(json!(v.get_poses_energies(uint(args, 0)?, float(args, 1)?)?))
        }
        "score" => {
            return Ok(json!(if args.as_array().is_some_and(|a| !a.is_empty()) {
                v.score_with_unbound(float(args, 0)?)?
            } else {
                v.score()?
            }))
        }
        "optimize" => {
            return Ok(json!(v.optimize(
                u32::try_from(uint(args, 0)?).map_err(|e| e.to_string())?
            )?))
        }
        "randomize" => v.randomize(uint(args, 0)?)?,
        "global_search" => v.global_search(
            uint(args, 0)?,
            uint(args, 1)?,
            float(args, 2)?,
            u32::try_from(uint(args, 3)?).map_err(|e| e.to_string())?,
        )?,
        _ => return Err(format!("Unknown method {method}")),
    }
    Ok(Value::Null)
}
fn main() {
    let mut engine = None;
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for line in io::stdin().lock().lines() {
        let result = line
            .map_err(|e| e.to_string())
            .and_then(|s| serde_json::from_str(&s).map_err(|e| e.to_string()))
            .and_then(|request| dispatch(&mut engine, request));
        let response = match result {
            Ok(value) => json!({"result":value}),
            Err(error) => json!({"error":error}),
        };
        if writeln!(out, "{response}")
            .and_then(|_| out.flush())
            .is_err()
        {
            break;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn protocol_reports_invalid_requests_and_uninitialized_state() {
        let mut engine = None;
        assert!(dispatch(&mut engine, json!({"method":"score","args":[]})).is_err());
        assert_eq!(
            dispatch(
                &mut engine,
                json!({"method":"new","args":["vina",1,42,0,false]})
            )
            .unwrap(),
            json!(42)
        );
        assert!(dispatch(&mut engine, json!({"method":"optimize","args":[-1]})).is_err());
        assert!(dispatch(&mut engine, json!({"method":"unknown","args":[]})).is_err());
    }
}
