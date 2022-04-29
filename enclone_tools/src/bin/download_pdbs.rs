// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Download gzipped PDB structures in antibody_sets/antibody_antigen_structures into
// antibody_sets/pdbs, if they have not already been downloaded.
//
// These need to be git added.

use io_utils::*;
use pretty_trace::*;
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let pdb_list = include_str!("antibody_antigen_structures");
    for line in pdb_list.lines() {
        if !line.starts_with('#') && line.len() > 0 {
            let x = line.after(" ").split(',').collect::<Vec<&str>>();
            for pdb in x.iter() {
                if !path_exists(&format!("antibody_sets/pdbs/{}.gz", pdb)) {
                    let x = Command::new("wget")
                        .arg("-q")
                        .arg("-O")
                        .arg(&format!("antibody_sets/pdbs/{}.gz", pdb))
                        .arg(&format!("https://files.rcsb.org/download/{}.cif.gz", pdb))
                        .output()
                        .expect(&format!("failed to execute wget"));
                    if x.status.code() != Some(0) {
                        eprintln!("\nDownload of {} failed.\n", pdb);
                        std::process::exit(1);
                    }
                }
            }
        }
    }
}
