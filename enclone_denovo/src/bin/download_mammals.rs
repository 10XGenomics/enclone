// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// usage: download_mammals genome_list directory
//
// Download the mammals in the genome_list, that are not already in the given directory.
//
// Creates a binary_vec_vec file.
//
// It would be better if this did not create a temporary fasta file.

use binary_vec_io::binary_write_vec_vec;
use fasta_tools::read_fasta_to_vec_vec_u8;
use pretty_trace::PrettyTrace;
use std::env;
use std::fs::{read_dir, File};
use std::process::Command;
use string_utils::TextUtils;
use vector_utils::bin_member;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let genome_list = &args[1];
    let dir = &args[2];
    let all = read_dir(&dir).unwrap();
    let mut owned = Vec::<String>::new();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        owned.push(f.rev_after("/").before(":").to_string());
    }
    owned.sort();
    let f = std::fs::read_to_string(&genome_list).unwrap();
    for line in f.lines() {
        if !line.contains("GC") {
            continue;
        }
        let acc = line.before(":").to_string();
        if !bin_member(&owned, &acc) {
            println!("downloading {}", line);
            let outname = format!("{}/{}", dir, line);
            let url = format!("https://www.ncbi.nlm.nih.gov/assembly/{}", acc);
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
                    if ps.len() > 1 {
                        for j in 0..ps.len() {
                            println!("[{}] {}", j + 1, ps[j]);
                        }
                    }
                    assert_eq!(ps.len(), 1);
                    println!("running ftp {}", ps[0]);
                    let fasta_file = format!("{}.fasta.gz", acc);
                    let _ = Command::new("wget")
                        .arg(&ps[0])
                        .arg(&format!("-O {}", fasta_file))
                        .output()
                        .expect("failed to execute wget");
                    // no idea why the space is needed
                    std::fs::rename(&format!(" {}", fasta_file), &fasta_file).unwrap();
                    let x = read_fasta_to_vec_vec_u8(&fasta_file);
                    binary_write_vec_vec::<u8>(&mut File::create(&outname).unwrap(), &x).unwrap();
                    let _ = std::fs::remove_file(&fasta_file);
                    break;
                }
            }
        }
    }
}
