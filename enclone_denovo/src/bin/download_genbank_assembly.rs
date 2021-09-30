// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Download a GenBank assembly.
// usage:
// download_genbank_assembly GCA_id.n
// in the directory where you want the file.
// Creates GCA_id.n.fasta.gz.

use pretty_trace::PrettyTrace;
use std::env;
use std::process::Command;
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let id = &args[1];
    let url = format!("https://www.ncbi.nlm.nih.gov/assembly/{}", id);
    let o = Command::new("curl")
        .arg(url)
        .output()
        .expect("failed to execute curl 1");
    let m = String::from_utf8(o.stdout).unwrap();
    for line in m.lines() {
        if line.contains("FTP directory for GenBank assembly") {
            let ftpdir = line.between("href=\"", "\"");
            let o = Command::new("curl")
                .arg(&format!("{}/", ftpdir))
                .output()
                .expect("failed to execute curl 2");
            let m = String::from_utf8(o.stdout).unwrap();
            let mut ps = Vec::<String>::new();
            for line in m.lines() {
                if line.contains("href=") {
                    let x = line.between("\"", "\"");
                    if x.ends_with(".fna.gz") {
                        let p = format!("{}/{}", ftpdir, x);
                        if p.contains("_rna_") || p.contains("_cds_") {
                            continue;
                        }
                        ps.push(p);
                    }
                }
            }
            if ps.len() != 1 {
                eprintln!("m = {}\n", m);
            }
            assert_eq!(ps.len(), 1);
            println!("running ftp {}", ps[0]);
            let fasta_file = format!("{}.fasta.gz", id);
            let _ = Command::new("wget")
                .arg(&ps[0])
                .arg(&format!("-O {}", fasta_file))
                .output()
                .expect("failed to execute wget");
            // no idea why this is needed
            std::fs::rename(&format!(" {}", fasta_file), &fasta_file).unwrap();
            std::process::exit(0);
        }
    }
    eprintln!("\nFailed to find FTP URL.\n");
    std::process::exit(1);
}
