// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Temporary code for playing with melting temperatures.

use dna::*;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();
    println!("");
    // enclone BCR=123085 PER_CELL AMINO=cdr3 MIN_CHAINS_EXACT=2 CDR3=CARDRIAGRFGYGMDVW 
    //         POUT=stdout PCOLS=barcode PCELL
    let s = ["ATGCGATGTCTCAACA", "CAACTAGGTAGAGGAA", "CATCCACCAGCGAACA", "CCTTTCTAGGACGAAA",
        "CGATCGGCATTGGCGC", "CTCGAAACAAGCCGTC", "TAAGCGTTCAAAGACA", "TGGCTGGAGGATGTAT",
        "TGTATTCAGGTCGGAT", "GATCAGTGTCGAGTTT"];
    for i in 0..s.len() {
        let t = tm_nearest_neighbor(&s[i]);
        println!( "{} = {} ==> {:.1}°", i+1, s[i], t );
    }
    let s = "TCCTGAGGACTGTAGGACAGC";
    let t = tm_nearest_neighbor(&s);
    println!( "\nIGHG = {} ==> {:.1}°\n", s, t );
}
