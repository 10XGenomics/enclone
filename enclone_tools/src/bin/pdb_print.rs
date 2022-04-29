// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use enclone_tools::pdb::*;
use pretty_trace::*;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let pdb_list = include_str!("antibody_antigen_structures");
    for line in pdb_list.lines() {
        if !line.starts_with('#') && line.len() > 0 {
            let x = line.after(" ").split(',').collect::<Vec<&str>>();
            let pdb = fetch_pdb_structure(&x[0]).unwrap();
            let (mut heavy, mut light) = (None, None);
            for i in 0..pdb.chain_names.len() {
                if pdb.chain_names[i].contains("light")
                    || pdb.chain_names[i].contains("Light")
                    || pdb.chain_names[i].contains("LIGHT")
                    || pdb.chain_names[i].contains("Fab L")
                    || pdb.chain_names[i].contains("L chain")
                    || pdb.chain_names[i].contains("Kappa")
                    || pdb.chain_names[i].contains("kappa")
                {
                    light = Some(i);
                } else if pdb.chain_names[i].contains("heavy")
                    || pdb.chain_names[i].contains("Heavy")
                    || pdb.chain_names[i].contains("HEAVY")
                    || pdb.chain_names[i].contains("H chain")
                    || pdb.chain_names[i].contains("Fab H")
                {
                    heavy = Some(i);
                }
            }
            let s = line.after(" ").replace(",", "=");
            println!(
                "{} {} {}",
                s,
                strme(&pdb.chains[heavy.unwrap()]),
                strme(&pdb.chains[light.unwrap()])
            );
        }
    }
}
