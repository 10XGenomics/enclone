// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Load known covid antibodies and compute.

use flate2::read::GzDecoder;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::io::Read;
use std::mem::swap;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();

    // Load the transition matrix.  Code copied.

    let f = include_bytes!["../../../enclone_paper/data/mat.575142"].to_vec();
    let f = stringme(&f);
    let mut m = Vec::<Vec<f64>>::new();
    for line in f.lines() {
        let mut s = line.to_string();
        if s.starts_with("    A") {
            continue;
        }
        if s.len() > 2 && s.as_bytes()[0] >= b'A' {
            s = s[2..].to_string();
        }
        let fields = s.split(' ').collect::<Vec<&str>>();
        let mut row = Vec::<f64>::new();
        for i in 0..fields.len() {
            row.push(fields[i].force_f64());
        }
        m.push(row);
    }
    let mut penalty = vec![vec![0.0; 27]; 27];
    let aa = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
    for i1 in 0..aa.len() {
        let c1 = (aa[i1] - b'A') as usize;
        for i2 in 0..aa.len() {
            let c2 = (aa[i2] - b'A') as usize;
            penalty[c1][c2] = m[i1][i2];
        }
    }

    // Load covid antibodies and parse.

    let f = include_bytes!["../../../enclone_paper/data/CoV-AbDab_040322.csv.gz"].to_vec();
    let mut d = GzDecoder::new(&*f);
    let mut s = String::new();
    d.read_to_string(&mut s).unwrap();
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(String, String, usize, String, String)>::new();
    for (i, line) in s.lines().enumerate() {
        let fields = parse_csv(&line);
        if i == 0 {
            for (j, x) in fields.iter().enumerate() {
                tof.insert(x.clone(), j);
            }
        } else {
            let vh = &fields[tof["Heavy V Gene"]];
            let vl = &fields[tof["Light V Gene"]];
            let h = &fields[tof["VH or VHH"]];

            // Get the CDR3H and extend by one amino acid on each end.

            let mut cdrh3 = fields[tof["CDRH3"]].clone();
            cdrh3 = cdrh3.replace(" ", "");
            if cdrh3 == "ND" || !h.contains(&cdrh3) {
                continue;
            }
            let hvec = &h.as_bytes();
            let h3vec = cdrh3.as_bytes();
            let mut pp = Vec::<usize>::new();
            for i in 0..=hvec.len() - h3vec.len() {
                if *h3vec == hvec[i..i + h3vec.len()] {
                    pp.push(i);
                }
            }
            if pp.len() != 1 {
                continue;
            }
            let p = pp[0];
            cdrh3 = stringme(&hvec[p - 1..p + cdrh3.len() + 1]);

            // Save.

            if vh.contains("Human") && vl.contains("Human") {
                let vh = vh.replace(" ", "");
                let vl = vl.replace(" ", "");
                let vl = vl.replace("D", "");
                let vh = vh.before("(").to_string();
                let vl = vl.before("(").to_string();
                let mut source = fields[tof["Sources"]].clone();
                if source.len() > 90 {
                    source.truncate(90);
                }
                data.push((vh, vl, cdrh3.len(), cdrh3, source));
            }
        }
    }
    unique_sort(&mut data);
    let mut to_delete = vec![false; data.len()];
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        if j - i == 1 {
            to_delete[i] = true;
        }
        i = j;
    }
    erase_if(&mut data, &to_delete);

    // Look for covid antibodies that are proximate.

    let mut total = 0;
    let mut same = 0;
    let mut errs = Vec::<(String, String)>::new();
    for i1 in 0..data.len() {
        for i2 in i1 + 1..data.len() {
            let x = &data[i1];
            let y = &data[i2];

            // Require same heavy chain gene and same CDRH3 length.

            if x.0 != y.0 || x.3.len() != y.3.len() {
                continue;
            }

            // Don't consider antibodies with same author, as best we can tell.

            if x.4.contains(" ") && y.4.contains(" ") {
                if x.4.after(" ").contains(" ") && y.4.after(" ").contains(" ") {
                    let t1 = x.4.after(" ").after(" ");
                    let a1 = t1.rev_before(&t1);
                    let t2 = y.4.after(" ").after(" ");
                    let a2 = t2.rev_before(&t2);
                    if a1 == a2 {
                        continue;
                    }
                }
            }

            // Don't consider antibodies with identical CDRH3.

            if x.3 == y.3 {
                continue;
            }

            // Compare.

            let c1 = x.3.as_bytes();
            let c2 = y.3.as_bytes();
            let mut err = 0.0;
            for i in 0..c1.len() {
                let d1 = (c1[i] - b'A') as usize;
                let d2 = (c2[i] - b'A') as usize;
                err += penalty[d1][d2];
            }
            err /= c1.len() as f64;
            if err <= 0.1 {
                if x.1 != y.1 {
                    println!("\n{}, {}, {}", x.0, x.1, x.3);
                    println!("{}, {}, {}", y.0, y.1, y.3);
                    println!("err = {:.2}", err);
                    println!("source1 = {}", x.4);
                    println!("source2 = {}", y.4);
                    let mut a1 = x.1.clone();
                    let mut a2 = y.1.clone();
                    if a1 > a2 {
                        swap(&mut a1, &mut a2);
                    }
                    errs.push((a1, a2));
                }
                total += 1;
                if x.1 == y.1 {
                    same += 1;
                }
            }
        }
    }
    let coherence = 100.0 * same as f64 / total as f64;
    println!("\nlight chain coherence = {:.1}% of {}\n", coherence, total);
    errs.sort();
    let mut freq = Vec::<(u32, (String, String))>::new();
    make_freq(&errs, &mut freq);
    for i in 0..freq.len() {
        println!("[{}] {}, {}", freq[i].0, freq[i].1 .0, freq[i].1 .1);
    }
    println!("");
}
