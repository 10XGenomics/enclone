// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Optimize amino acid penalty matrix to increase joining given light chain coherence ≥ 75%.
//
// This is very delicate piece of this that is easy to break.  Namely, do we correctly handle
// cases such as 1 amino acid different out of 10, which should be 90% identity.  This can be
// broken because of floating point arithmetic.
//
// Usage: light_chain_coherence per_cell_stuff
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1 > per_cell_stuff
//
// or with more fields
//
// Output:
//
// ...
// count = 272920, nrel = 8.2689, light chain coherence = 75.0%, used 1626.1 minutes
//
//    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
// A 0.0 0.1 0.6 0.4 0.4 0.4 0.5 0.4 0.3 0.1 0.2 0.4 0.1 8.0 0.4 0.4 0.2 0.4 4.6 0.7
// C 0.1 0.0 4.9 0.2 1.5 1.0 0.3 0.2 0.2 0.1 0.2 0.2 0.5 4.5 6.5 0.3 0.4 0.2 8.0 0.2
// D 0.6 4.9 0.0 0.2 0.7 1.0 0.9 8.0 6.7 1.1 8.0 0.6 8.0 0.7 1.2 0.9 0.7 0.8 8.0 0.8
// E 0.4 0.2 0.2 0.0 8.0 0.9 0.2 0.2 0.3 0.8 6.8 5.7 0.6 0.4 0.9 0.8 0.7 0.6 8.0 8.0
// F 0.4 1.5 0.7 8.0 0.0 0.7 0.2 0.6 0.9 0.7 8.0 0.6 0.9 8.0 0.9 0.4 0.3 0.7 0.4 0.5
// G 0.4 1.0 1.0 0.9 0.7 0.0 8.0 6.3 0.4 1.1 0.5 0.7 1.0 6.9 0.8 0.6 0.6 0.6 1.0 1.4
// H 0.5 0.3 0.9 0.2 0.2 8.0 0.0 8.0 8.0 0.8 0.3 4.2 1.0 0.6 8.0 8.0 0.3 0.5 5.7 0.3
// I 0.4 0.2 8.0 0.2 0.6 6.3 8.0 0.0 0.1 0.1 0.2 8.0 0.7 0.5 0.8 0.7 0.1 0.1 6.4 0.6
// K 0.3 0.2 6.7 0.3 0.9 0.4 8.0 0.1 0.0 0.1 0.1 0.1 8.0 0.1 0.2 8.0 8.0 0.2 0.1 0.1
// L 0.1 0.1 1.1 0.8 0.7 1.1 0.8 0.1 0.1 0.0 0.4 0.4 8.0 0.7 0.9 0.2 0.1 0.2 0.6 0.5
// M 0.2 0.2 8.0 6.8 8.0 0.5 0.3 0.2 0.1 0.4 0.0 8.0 0.7 0.4 0.8 0.5 0.7 0.3 0.5 0.6
// N 0.4 0.2 0.6 5.7 0.6 0.7 4.2 8.0 0.1 0.4 8.0 0.0 0.5 0.2 0.9 0.5 8.0 0.3 1.0 0.6
// P 0.1 0.5 8.0 0.6 0.9 1.0 1.0 0.7 8.0 8.0 0.7 0.5 0.0 5.8 1.1 8.0 8.0 0.8 8.0 8.0
// Q 8.0 4.5 0.7 0.4 8.0 6.9 0.6 0.5 0.1 0.7 0.4 0.2 5.8 0.0 8.0 0.7 0.4 0.2 0.7 0.5
// R 0.4 6.5 1.2 0.9 0.9 0.8 8.0 0.8 0.2 0.9 0.8 0.9 1.1 8.0 0.0 0.7 0.5 1.0 8.0 1.2
// S 0.4 0.3 0.9 0.8 0.4 0.6 8.0 0.7 8.0 0.2 0.5 0.5 8.0 0.7 0.7 0.0 0.3 0.3 8.0 0.5
// T 0.2 0.4 0.7 0.7 0.3 0.6 0.3 0.1 8.0 0.1 0.7 8.0 8.0 0.4 0.5 0.3 0.0 0.1 8.0 1.0
// V 0.4 0.2 0.8 0.6 0.7 0.6 0.5 0.1 0.2 0.2 0.3 0.3 0.8 0.2 1.0 0.3 0.1 0.0 0.7 0.5
// W 4.6 8.0 8.0 8.0 0.4 1.0 5.7 6.4 0.1 0.6 0.5 1.0 8.0 0.7 8.0 8.0 8.0 0.7 0.0 0.7
// Y 0.7 0.2 0.8 8.0 0.5 1.4 0.3 0.6 0.1 0.5 0.6 0.6 8.0 0.5 1.2 0.5 1.0 0.5 0.7 0.0

use io_utils::*;
use perf_stats::elapsed;
use pretty_trace::PrettyTrace;
use rand_chacha;
use rand_chacha::rand_core::RngCore;
use rand_chacha::rand_core::SeedableRng;
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
    let mut randme = rand_chacha::ChaCha8Rng::seed_from_u64(123456789);
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
    let mut best_n = 0;
    let mut canonical_n = 0;
    let mut changed = false;
    for count in 1.. {
        if count % 10000 == 0 {
            println!(
                "\n(count = {}, time = {:.1} minutes)",
                add_commas(count),
                elapsed(&t) / 60.0
            );
        }

        // Mutate penalty matrix.

        let penalty_save = penalty.clone();
        if count > 1 {
            let rand1 = randme.next_u64();
            let rand2 = randme.next_u64();
            let rand3 = randme.next_u64();
            let rand4 = randme.next_u64();
            let rand5 = randme.next_u64();
            let pert = (rand1 % 1_000_000u64) as f64 / 1_000_000.0; // in [0,1)
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
            penalty[a1][a2] = penalty[a1][a2].min(8.0);
            penalty[a2][a1] = penalty[a2][a1].min(8.0);
            penalty[b1][b2] = penalty[b1][b2].min(8.0);
            penalty[b2][b1] = penalty[b2][b1].min(8.0);
            penalty[a1][a2] = format!("{:.1}", penalty[a1][a2]).parse::<f64>().unwrap();
            penalty[a2][a1] = format!("{:.1}", penalty[a2][a1]).parse::<f64>().unwrap();
            penalty[b1][b2] = format!("{:.1}", penalty[b1][b2]).parse::<f64>().unwrap();
            penalty[b2][b1] = format!("{:.1}", penalty[b2][b1]).parse::<f64>().unwrap();
        }
        let mut penaltyx = Vec::<f64>::new();
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
                err /= data[k1].2.len() as f64;

                // Add to results.

                if 1.0 - err >= 0.9 {
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
        let nznz = 100.0 * res.0 as f64 / n as f64;
        if n > best_n && nznz >= 75.0 {
            changed = true;
            if count > 1 {
                let nrel = n as f64 / canonical_n as f64;
                print!("count = {count}, nrel = {nrel:.4}, light chain coherence = {nznz:.1}%");
                println!(", used {:.1} minutes", elapsed(&t) / 60.0);
            }
            best_n = n;
        } else {
            penalty = penalty_save;
        }
        if count % 500 == 0 && changed {
            changed = false;
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
                    row.push(format!("{:.1}", penalty[i][j]));
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
