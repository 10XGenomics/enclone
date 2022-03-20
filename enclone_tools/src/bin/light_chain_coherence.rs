// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Optimize amino acid penalty matrix to increase joining given light chain coherence of at
// least 75%.

// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > per_cell_stuff

use io_utils::*;
use perf_stats::elapsed;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use std::time::Instant;
use string_utils::*;
use tables::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(String, usize, Vec<u8>, String, String, usize)>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("donors_cell"));
            assert!(tof.contains_key("v_name1"));
            assert!(tof.contains_key("v_name2"));
            assert!(tof.contains_key("dref"));
            assert!(tof.contains_key("cdr3_aa1"));
            first = false;
        } else {
            let dref = fields[tof["dref"]].force_usize();
            if dref > 0 {
                data.push((
                    fields[tof["v_name1"]].to_string(),
                    fields[tof["cdr3_aa1"]].len(),
                    fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                    fields[tof["donors_cell"]].to_string(),
                    fields[tof["v_name2"]].to_string(),
                    fields[tof["dref"]].force_usize(),
                ));
            }
        }
    }
    data.sort();

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Replace light chain genes and donors by numbers.

    let mut datax = Vec::<(String, usize, Vec<u8>, usize, usize, usize)>::new();
    let mut lights = Vec::<String>::new();
    let mut donors = Vec::<String>::new();
    for i in 0..data.len() {
        lights.push(data[i].4.clone());
        donors.push(data[i].3.clone());
    }
    unique_sort(&mut lights);
    unique_sort(&mut donors);
    for i in 0..data.len() {
        datax.push((
            data[i].0.clone(),
            data[i].1,
            data[i].2.clone(),
            bin_position(&donors, &data[i].3) as usize,
            bin_position(&lights, &data[i].4) as usize,
            data[i].5,
        ));
    }
    let mut data = datax;

    // Convert CDRH3 to [0,20).

    let aa = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
    for i in 0..data.len() {
        let mut x = Vec::<u8>::new();
        for j in 0..data[i].2.len() {
            let c = data[i].2[j];
            let p = bin_position(&aa, &c) as u8;
            x.push(p);
        }
        data[i].2 = x;
    }

    // Define penalty matrix.

    let mut penalty = vec![vec![1.0; 20]; 20];
    for i in 0..20 {
        penalty[i][i] = 0.0;
    }

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, (usize, usize))>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, (0, 0)));
        i = j;
    }

    // Form buckets.

    let bucket_size = 100_000;
    let mut buckets = Vec::<(Vec<(usize, usize)>, (usize, usize))>::new();
    {
        let mut bucket = Vec::<(usize, usize)>::new();
        for m in 0..bounds.len() {
            let i = bounds[m].0;
            let j = bounds[m].1;
            for k1 in i..j {
                for k2 in k1 + 1..j {
                    // Require different donors.

                    if data[k1].3 == data[k2].3 {
                        continue;
                    }
                    bucket.push((k1, k2));
                    if bucket.len() == bucket_size {
                        buckets.push((bucket.clone(), (0, 0)));
                        bucket.clear();
                    }
                }
            }
        }
        if bucket.len() > 0 {
            buckets.push((bucket, (0, 0)));
        }
    }

    // Loop.

    println!("");
    let mut rand = 0i64;
    let mut best_n = 0;
    let mut canonical_n = 0;
    for count in 1.. {
        // Mutate penalty matrix.

        let penalty_save = penalty.clone();
        if count > 0 {
            let rand1 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand2 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand1)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand3 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand2)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand4 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand3)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand5 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand4)
                .wrapping_add(1_442_695_040_888_963_407);
            rand = rand5;
            let pert = (rand1 % 1_000_000i64) as f32 / 1_000_000.0; // in [0,1)
            let mul = 1.0 + pert;
            let a1 = (rand2 as usize) % 20;
            let a2 = (rand3 as usize) % 20;
            let b1 = (rand4 as usize) % 20;
            let b2 = (rand5 as usize) % 20;
            if a1 == a2 || b1 == b2 {
                continue;
            }
            penalty[a1][a2] *= mul;
            penalty[a2][a1] *= mul;
            penalty[b1][b2] /= mul;
            penalty[b2][b1] /= mul;
            penalty[a1][a2] = format!("{:.3}", penalty[a1][a2]).parse::<f32>().unwrap();
            penalty[a2][a1] = format!("{:.3}", penalty[a2][a1]).parse::<f32>().unwrap();
            penalty[b1][b2] = format!("{:.3}", penalty[b1][b2]).parse::<f32>().unwrap();
            penalty[b2][b1] = format!("{:.3}", penalty[b2][b1]).parse::<f32>().unwrap();
        }
        let mut penaltyx = Vec::<f32>::new();
        for i in 0..20 {
            for j in 0..20 {
                penaltyx.push(penalty[i][j]);
            }
        }

        // Results = for each percent identity, rounded down:
        // 1. count for equal light chain gene names and dref1 > 0 and dref2 > 0
        // 2. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.

        buckets.par_iter_mut().for_each(|res| {
            for m in 0..res.0.len() {
                let k1 = res.0[m].0;
                let k2 = res.0[m].1;

                // Compute stuff.

                let mut err = 0.0;
                for m in 0..data[k1].2.len() {
                    let c1 = data[k1].2[m] as usize;
                    let c2 = data[k2].2[m] as usize;
                    err += penaltyx[20 * c1 + c2];
                }
                err /= data[k1].2.len() as f32;

                // Add to results.

                if err <= 0.1 {
                    let eq_light = data[k1].4 == data[k2].4;
                    if eq_light {
                        res.1 .0 += 1;
                    } else {
                        res.1 .1 += 1;
                    }
                }
            }
        });

        // Sum.

        let mut res = (0, 0);
        for i in 0..buckets.len() {
            res.0 += buckets[i].1 .0;
            res.1 += buckets[i].1 .1;
        }

        // Print.

        let n = res.0 + res.1;
        if count == 1 {
            canonical_n = n;
        }
        let nznz = 100.0 * res.0 as f32 / n as f32;
        if n > best_n && nznz >= 75.0 {
            if count > 1 {
                let nrel = n as f32 / canonical_n as f32;
                print!("count = {count}, nrel = {nrel:.4}, light chain coherence = {nznz:.1}%");
                println!(", used {:.1} minutes", elapsed(&t) / 60.0);
            }
            best_n = n;
        } else {
            penalty = penalty_save;
        }
        if count % 200 == 0 {
            let mut rows = Vec::<Vec<String>>::new();
            let mut row = Vec::<String>::new();
            row.push(String::new());
            for i in 0..20 {
                row.push((aa[i] as char).to_string());
            }
            rows.push(row);
            for i in 0..20 {
                let mut row = Vec::<String>::new();
                row.push((aa[i] as char).to_string());
                for j in 0..20 {
                    row.push(format!("{:.3}", penalty[i][j]));
                }
                rows.push(row);
            }
            let mut log = Vec::<u8>::new();
            print_tabular(&mut log, &rows, 1, Some(vec![b'r'; 21]));
            println!("\n{}", strme(&log));
        }

        // Reset.

        for i in 0..buckets.len() {
            buckets[i].1 = (0, 0);
        }
    }
}
