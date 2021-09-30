// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Read V gene predictions, extract features, and line things up.
//
// usage: squeeze_novo feature length
// e.g. squeeze_novo FWR3 39
// optional extra arg: a string that full reference names should match

use amino::*;
use enclone_denovo::vdj_features::*;
use fasta_tools::*;
use pretty_trace::*;
use std::env;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let feature = &args[1];
    let len = args[2].force_usize();
    let mut required = String::new();
    if args.len() == 4 {
        required = args[3].to_string();
    }
    let dir = "/mnt/deck5/david.jaffe/denovo_ref";
    let refs = std::fs::read_dir(&dir).unwrap();
    for f in refs {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let mut g = f.rev_after("/").rev_before(".fasta").to_string();
        if g.contains(&required) {
            if g.contains(':') {
                g = format!("{}:{}", g.before(":"), g.rev_after(":"));
            }
            let x = read_fasta_to_vec_vec_u8(f);
            for i in (0..x.len()).step_by(2) {
                let chain_type = strme(&x[i]).before("V");
                if chain_type == "IGH" {
                    let aa = aa_seq(&x[i + 1], 0);
                    let bb;
                    if feature == "FWR1" {
                        bb = fwr1(&aa, chain_type, false).unwrap();
                    } else if feature == "FWR2" {
                        bb = fwr2(&aa, chain_type, false).unwrap();
                    } else if feature == "FWR3" {
                        bb = fwr3(&aa, chain_type, false).unwrap();
                    } else if feature == "CDR1" {
                        bb = cdr1(&aa, chain_type, false).unwrap();
                    } else if feature == "CDR2" {
                        bb = cdr2(&aa, chain_type, false).unwrap();
                    } else {
                        eprintln!("\nUnrecognized feature.\n");
                        std::process::exit(1);
                    }
                    if bb.len() == len {
                        println!("{}  {}", strme(&bb), g);
                    }
                }
            }
        }
    }
}
