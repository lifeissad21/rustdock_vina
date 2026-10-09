use rustdock_vina_split::{default_prefix, output_names, parse_multimodel_pdbqt, write_pdbqt};

fn main() {
    let mut args = std::env::args().skip(1);
    let mut input = None;
    let mut ligand_prefix = None;
    let mut flex_prefix = None;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" => input = args.next(),
            "--ligand" => ligand_prefix = args.next(),
            "--flex" => flex_prefix = args.next(),
            "--help" => {
                println!("Usage: vina_split --input <PDBQT> [--ligand <prefix>] [--flex <prefix>]");
                return;
            }
            "--version" => {
                println!("AutoDock Vina PDBQT Split rust-port");
                return;
            }
            _ => {
                eprintln!("Unknown argument: {arg}");
                std::process::exit(1);
            }
        }
    }

    let Some(input) = input else {
        eprintln!("Missing input.");
        std::process::exit(1);
    };

    let ligand_prefix = ligand_prefix.unwrap_or_else(|| default_prefix(&input, "_ligand_"));
    let flex_prefix = flex_prefix.unwrap_or_else(|| default_prefix(&input, "_flex_"));
    let models = match parse_multimodel_pdbqt(&input) {
        Ok(models) => models,
        Err(err) => {
            eprintln!("\n\nPDBQT parsing error: {err}\n");
            std::process::exit(1);
        }
    };

    for (model, (ligand_name, flex_name)) in
        models
            .iter()
            .zip(output_names(models.len(), &ligand_prefix, &flex_prefix))
    {
        if let Err(err) = write_pdbqt(&model.ligand, ligand_name)
            .and_then(|_| write_pdbqt(&model.flex, flex_name))
        {
            eprintln!("{err}");
            std::process::exit(1);
        }
    }
}
