use std::{env, path::PathBuf, process::Command};
fn main() {
    println!("cargo:rustc-check-cfg=cfg(has_metal)");
    println!("cargo:rerun-if-changed=metal/main.swift");
    println!("cargo:rerun-if-changed=metal/grid_score.metal");
    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("macos")
        || env::var_os("CARGO_FEATURE_METAL").is_none()
    {
        return;
    }
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("rustdock-vina-metal");
    let arch = match env::var("CARGO_CFG_TARGET_ARCH").as_deref() {
        Ok("aarch64") => "arm64",
        Ok("x86_64") => "x86_64",
        _ => return,
    };
    let target = format!("{arch}-apple-macosx13.0");
    let result = Command::new("xcrun")
        .args(["swiftc", "-target", &target, "-module-cache-path"])
        .arg(output.parent().unwrap().join("swift-cache"))
        .args(["-O", "metal/main.swift", "-o"])
        .arg(&output)
        .status();
    if result.is_ok_and(|s| s.success()) {
        println!("cargo:rustc-cfg=has_metal");
    } else {
        println!("cargo:warning=Metal helper could not be built; CPU backend remains available. Install Apple Command Line Tools and rebuild to enable Metal.");
    }
}
