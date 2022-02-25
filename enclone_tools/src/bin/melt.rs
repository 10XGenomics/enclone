// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Temporary code for playing with melting temperatures.
//
// Barcodes are from:
// enclone PRE= 153782,153783,153784 PER_CELL SEG="IGHG1|IGHG2|IGHG3|IGHG4"
// Only used IGHG, excluded clonotypes having three chains.

use dna::*;
use pretty_trace::*;
use string_utils::*;
use tables::*;

fn main() {
    PrettyTrace::new().on();

    // Define barcodes.  Data show (id = clonotype.cell, barcode, constant regions),
    // where the numbering scheme for clonotype and cell is for this list.

    let mut source = Vec::<(String, Vec<u8>, String)>::new();
    source.push((
        "1 = 1.1".to_string(),
        b"ATGTGTGAGACTGGGT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "2 = 1.2".to_string(),
        b"CACACCTAGATACACA".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "3 = 1.3".to_string(),
        b"GTGCAGCCATTTCACT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "4 = 1.4".to_string(),
        b"TGACGGCCAATGGAGC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "5 = 2.1".to_string(),
        b"GGGAATGGTAAGGGCT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "6 = 3.1".to_string(),
        b"AGCAGCCAGATGTGGC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "7 = 3.2".to_string(),
        b"TGTTCCGGTTACTGAC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "8 = 4.1".to_string(),
        b"AAAGTAGTCCTCCTAG".to_vec(),
        "IGHG3, IGKC".to_string(),
    ));
    source.push((
        "9 = 5.1".to_string(),
        b"ATCTACTCACGAAACG".to_vec(),
        "IGHG1, IGLC1".to_string(),
    ));
    source.push((
        "10 = 5.2".to_string(),
        b"TTAGGACCAGTACACT".to_vec(),
        "IGHG1, IGLC1".to_string(),
    ));
    source.push((
        "11 = 6.1".to_string(),
        b"TCGTACCGTTCCATGA".to_vec(),
        "IGHG1, IGLC3".to_string(),
    ));
    source.push((
        "12 = 6.2".to_string(),
        b"TGCTGCTGTAAATACG".to_vec(),
        "IGHG1, IGLC3".to_string(),
    ));

    // to delete:
    // enclone BCR=123085 PER_CELL AMINO=cdr3 MIN_CHAINS_EXACT=2 CDR3=CARDRIAGRFGYGMDVW
    //         POUT=stdout PCOLS=barcode PCELL

    // Compute Tm for the relevant constant region primers.

    println!("");
    let s = "TCCTGAGGACTGTAGGACAGC";
    let tg = tm_nearest_neighbor(&s);
    println!("IGHG = {} ==> {:.1}°", s, tg);
    let s = "TAGCTGCTGGCCGC";
    let tl = tm_nearest_neighbor(&s);
    println!("IGLC = {} ==> {:.1}°", s, tl);
    let s = "GCGTTATCCACCTTCCACTGT";
    let tk = tm_nearest_neighbor(&s);
    println!("IGKC = {} ==> {:.1}°", s, tk);
    println!("");
    let tx = 57.5;

    // Sequence of pR1.  This is immediately before the barcode.

    let pr1 = b"CTACACGACGCTCTTCCGATCT";
    let npr = pr1.len();

    // Compute Tm for ...

    /*
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
    */
    let mut rows = Vec::<Vec<String>>::new();
    let row = vec![
        "#".to_string(),
        "const".to_string(),
        "method1".to_string(),
        "t1".to_string(),
        "method2".to_string(),
        "t2".to_string(),
    ];
    rows.push(row);
    rows.push(vec!["\\hline".to_string(); 6]);
    for i in 0..source.len() {
        let b = source[i].1.clone(); // barcode

        // Method 1: (tail of pR1 + barcode).

        let mut res = Vec::<(f64, f64, String)>::new();
        for j in 0..=pr1.len() {
            let x = format!("{}{}", strme(&pr1[npr - j..npr]), strme(&b));
            let mut xp = format!("{}+{}", strme(&pr1[npr - j..npr]), strme(&b));
            if j == 0 {
                xp = x.clone();
            }
            let t = tm_nearest_neighbor(&x);
            res.push(((t - tx).abs(), t, xp));
        }
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let t1 = res[0].1;
        let xp1 = res[0].2.clone();

        // Method 2: (tail of pR1 + first 8 bases of barcode).

        let b = &b[0..8];
        let mut res = Vec::<(f64, f64, String)>::new();
        for j in 0..=pr1.len() {
            let x = format!("{}{}", strme(&pr1[npr - j..npr]), strme(&b));
            let mut xp = format!("{}+{}", strme(&pr1[npr - j..npr]), strme(&b));
            if j == 0 {
                xp = x.clone();
            }
            let t = tm_nearest_neighbor(&x);
            res.push(((t - tx).abs(), t, xp));
        }
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let t2 = res[0].1;
        let xp2 = res[0].2.clone();

        // Add row to table.

        let row = vec![
            source[i].0.clone(),
            source[i].2.clone(),
            xp1,
            format!("{:.1}°", t1),
            xp2,
            format!("{:.1}°", t2),
        ];
        rows.push(row);
    }
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, &b"l|l|l|l|l|l".to_vec(), false, false);
    println!("{}", log);
    for i in 0..rows.len() {
        if i == 1 {
            continue;
        }
        for j in 0..rows[i].len() {
            let mut s = rows[i][j].clone();
            if s.contains(',') {
                s = format!("\"{}\"", s);
            }
            if j > 0 {
                print!(",");
            }
            print!("{}", s);
        }
        println!("");
    }
}
