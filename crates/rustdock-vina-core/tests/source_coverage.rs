use std::collections::BTreeSet;
use std::path::Path;

use rustdock_vina_core::porting::{PortStatus, SOURCE_PORTS};

#[test]
fn every_reference_source_file_has_a_rust_counterpart() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest_dir
        .parent()
        .and_then(Path::parent)
        .expect("core crate should live under crates/rustdock-vina-core");
    let reference_src = repo_root.join("reference/src");

    let mut reference_files = BTreeSet::new();
    collect_source_files(&reference_src, repo_root, &mut reference_files);

    let mapped_reference_files = SOURCE_PORTS
        .iter()
        .map(|port| port.reference_path)
        .collect::<BTreeSet<_>>();

    for file in &reference_files {
        assert!(
            mapped_reference_files.contains(file.as_str()),
            "missing Rust mapping for {file}"
        );
    }

    for port in SOURCE_PORTS {
        assert_eq!(
            port.status,
            PortStatus::Ported,
            "Unfinished source port: {}",
            port.reference_path
        );
        let rust = std::fs::read_to_string(repo_root.join(port.rust_path)).unwrap();
        assert!(
            !rust.contains("Port scaffold"),
            "Source mapping points to a placeholder: {}",
            port.rust_path
        );
        assert!(
            repo_root.join(port.reference_path).is_file(),
            "mapped reference file does not exist: {}",
            port.reference_path
        );
        assert!(
            repo_root.join(port.rust_path).is_file(),
            "mapped Rust counterpart does not exist for {}: {}",
            port.reference_path,
            port.rust_path
        );
    }
}

fn collect_source_files(path: &Path, repo_root: &Path, out: &mut BTreeSet<String>) {
    for entry in std::fs::read_dir(path).expect("reference source directory should be readable") {
        let entry = entry.expect("reference source entry should be readable");
        let path = entry.path();
        if path.is_dir() {
            collect_source_files(&path, repo_root, out);
            continue;
        }
        let extension = path.extension().and_then(|value| value.to_str());
        if matches!(extension, Some("cpp" | "h")) {
            let normalized = path
                .strip_prefix(repo_root)
                .expect("source file should be under repo root")
                .to_string_lossy()
                .replace('\\', "/");
            let reference_index = normalized
                .find("reference/src/")
                .expect("source file should be under reference/src");
            out.insert(normalized[reference_index..].to_string());
        }
    }
}
