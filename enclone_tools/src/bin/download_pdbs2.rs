// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// download_pdbs2 list-of-structures target_directory
//
// Download pdb files.  The list should have one line per pdb file.  Blank lines and lines starting
// with # are ignored.  If a file has already been downloaded to the location, it will not be
// downloaded again.  Files are gzipped.
//
// This is actually quite fast, only a second or so per pdb file.
//
// This runs between
// 1. antibody_sets/candidate_structures
// 2. group_structures_by_antigen

use io_utils::*;
use pretty_trace::*;
use std::env;
use std::io::BufRead;
use std::process::Command;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (list, dir) = (&args[1], &args[2]);
    let f = open_for_read![&list];
    let mut count = 0;
    for line in f.lines() {
        let pdb = line.unwrap();
        if !pdb.starts_with('#') && pdb.len() > 0 {
            count += 1;
            if !path_exists(&format!("{}/{}.gz", dir, pdb)) {
                println!("downloading {} = {}", count, pdb);
                let x = Command::new("wget")
                    .arg("-q")
                    .arg("-O")
                    .arg(&format!("{}/{}.gz", dir, pdb))
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
