// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Mine STCRDab  = "The Structural T-Cell Receptor Database"
// http://opig.stats.ox.ac.uk/webapps/stcrdab
// to get CDR{1,2,3} truth data.  Data are for mouse and human.
//
// This was run on 7/22/20 and stdout captured as cdr_truth_data.fasta.

use pretty_trace::*;
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();

    // Traverse the species.

    for pass in 1..=2 {
        let species;
        if pass == 1 {
            species = "mouse";
        } else {
            species = "human";
        }

        // Get list of all TCR structures.

        let mut structs = Vec::<(String, String, String)>::new();
        let url = format!(
            "http://opig.stats.ox.ac.uk/webapps/stcrdab/Browser?TCRtype=All&species={}&\
            subgroup=All&MHCtype=All&agtype=All&peptide=&resolution=&r_factor=&aglen=&affinity=All",
            species
        );
        let o = Command::new("curl")
            .arg(url)
            .output()
            .expect("failed to execute opig http");
        let m = String::from_utf8(o.stdout).unwrap();
        for s in m.lines() {
            if s.contains("<td>X-RAY DIFFRACTION</td>") {
                let x = s.between("<tr><td>", "</td>").to_string();
                if s.contains("<br>VA:") {
                    let va = s.between("<br>VA: ", "<br>").to_string();
                    structs.push((x.clone(), "A".to_string(), va));
                }
                if s.contains("<br>VB:") {
                    let vb = s.between("<br>VB: ", "<br>").to_string();
                    structs.push((x.clone(), "B".to_string(), vb));
                }
            }
        }

        // Get the fasta files for each.

        for (x, y, v) in structs.iter() {
            let url = format!(
                "http://opig.stats.ox.ac.uk/webapps/stcrdab/Fasta?\
                chain_type={}&pdb={}&chain={}",
                y, x, v
            );
            let o = Command::new("curl")
                .arg(url)
                .output()
                .expect("failed to execute opig http");
            if o.status.code().unwrap() != 0 {
                println!("\nfailed\n");
                std::process::exit(1);
            }
            let m = String::from_utf8(o.stdout).unwrap();
            if m.contains("Internal Server Error") {
                // This happened a couple of times, not sure exactly why.
                // It is replicable at the same place.
            } else {
                print!("{}", m);
            }
        }
    }
}
