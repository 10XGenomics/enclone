// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Temporary code for playing with melting temperatures.

use dna::*;
use pretty_trace::*;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    // enclone BCR=123085 PER_CELL AMINO=cdr3 MIN_CHAINS_EXACT=2 CDR3=CARDRIAGRFGYGMDVW
    //         POUT=stdout PCOLS=barcode PCELL

    // Compute Tm for IGHG primer.

    let s = "TCCTGAGGACTGTAGGACAGC";
    let tg = tm_nearest_neighbor(&s);
    println!("\nIGHG = {} ==> {:.1}°\n", s, tg);

    // Compute Tm for ...

    let pr1 = b"CTACACGACGCTCTTCCGATCT";
    let s = [
        "ATGCGATGTCTCAACA",
        "CAACTAGGTAGAGGAA",
        "CATCCACCAGCGAACA",
        "CCTTTCTAGGACGAAA",
        "CGATCGGCATTGGCGC",
        "CTCGAAACAAGCCGTC",
        "TAAGCGTTCAAAGACA",
        "TGGCTGGAGGATGTAT",
        "TGTATTCAGGTCGGAT",
        "GATCAGTGTCGAGTTT",
    ];
    for i in 0..s.len() {
        let mut res = Vec::<(f64, f64, String)>::new();
        for j in 0..=pr1.len() {
            let x = format!("{}{}", strme(&pr1[pr1.len() - j..pr1.len()]), s[i]);
            let mut xp = format!("{}+{}", strme(&pr1[pr1.len() - j..pr1.len()]), s[i]);
            if j == 0 {
                xp = x.clone();
            }
            let t = tm_nearest_neighbor(&x);
            res.push(((t - tg).abs(), t, xp));
        }
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let t = res[0].1;
        let xp = res[0].2.clone();
        println!("{} = {} ==> {:.1}°", i + 1, xp, t);
    }
    println!("");
}
