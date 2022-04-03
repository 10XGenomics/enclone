// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// insertion_analysis per_cell_stuff

use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::{stringme, strme, TextUtils};
use vector_utils::{make_freq, unique_sort};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load data.

    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(
        String,
        usize,
        Vec<u8>,
        String,
        String,
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        String,
    )>::new();
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
            first = false;
        } else {
            data.push((
                /* 0 */ fields[tof["v_name1"]].to_string(),
                /* 1 */ fields[tof["cdr3_aa1"]].len(),
                /* 2 */ fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                /* 3 */ fields[tof["donors_cell"]].to_string(),
                /* 4 */ fields[tof["v_name2"]].to_string(),
                /* 5 */ fields[tof["dref"]].force_usize(),
                /* 6 */ fields[tof["jun_mat"]].force_usize(),
                /* 7 */ fields[tof["jun_sub"]].force_usize(),
                /* 8 */ fields[tof["hcomp"]].force_usize(),
                /* 9 */ fields[tof["jun_ins"]].force_usize(),
                /* 10 */ fields[tof["datasets_cell"]].force_usize(),
                /* 11 */ fields[tof["d1_name1"]].to_string(),
            ));
        }
    }
    data.sort();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // For naive cells with junction insertion length zero, show the substitution and
    // substitution rate distribution.

    println!("\nsubstitutions for naive cells with junction insertion length 0\n");
    for pass in 0..5 {
        if pass == 0 {
            println!("all donors");
        } else {
            println!("\ndonor {pass}");
        }
        let mut x = Vec::<usize>::new();
        let mut rates = Vec::<f64>::new();
        for k in 0..data.len() {
            let dref = data[k].5;
            let jun_ins = data[k].9;
            if dref == 0 && jun_ins == 0 {
                let donor = &data[k].3;
                if pass == 0 || *donor == format!("d{pass}") {
                    let matches = data[k].6;
                    let subs = data[k].7;
                    x.push(subs);
                    let rate = subs as f64 / (subs + matches) as f64;
                    rates.push(rate);
                }
            }
        }
        x.sort();
        let mut freq = Vec::<(u32, usize)>::new();
        make_freq(&x, &mut freq);
        let total: usize = x.iter().sum();
        let total_rates: f64 = rates.iter().sum();
        if pass == 0 {
            println!(
                "\nmost frequent substitution values for naive cells with junction \
                    insertion length 0 (of {})",
                x.len()
            );
            for i in 0..10 {
                println!(
                    "{} [{:.1}%]",
                    freq[i].1,
                    100.0 * freq[i].0 as f64 / x.len() as f64
                );
            }
            println!("mean = {:.1}", total as f64 / x.len() as f64);
            println!("\nsubstitution rates");
            let mut bins = vec![0; 21];
            for i in 0..rates.len() {
                bins[(20.0 * rates[i]).floor() as usize] += 1;
            }
            for i in 0..bins.len() {
                if bins[i] > 0 {
                    println!(
                        "{}-{}% ==> {:.1}%",
                        5 * i,
                        5 * (i + 1),
                        100.0 * bins[i] as f64 / rates.len() as f64
                    );
                }
            }
            println!("");
        }
        println!(
            "mean substitution rate = {:.2}%",
            100.0 * total_rates / rates.len() as f64
        );
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Compute substitution rates for public naive cells having jun_ins = 0.

    let mut i = 0;
    let mut rates = Vec::<f64>::new();
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        let mut donors = Vec::<String>::new();
        for k in i..j {
            if data[k].5 == 0 {
                donors.push(data[k].3.clone());
            }
        }
        unique_sort(&mut donors);
        if donors.len() > 1 {
            for k in i..j {
                let dref = data[k].5;
                let jun_ins = data[k].9;
                if dref == 0 && jun_ins == 0 {
                    let matches = data[k].6;
                    let subs = data[k].7;
                    let rate = subs as f64 / (subs + matches) as f64;
                    rates.push(rate);
                }
            }
        }
        i = j;
    }
    let total: f64 = rates.iter().sum();
    println!("\nsubstitution rates for public naive cells");
    let mut bins = vec![0; 21];
    for i in 0..rates.len() {
        bins[(20.0 * rates[i]).floor() as usize] += 1;
    }
    for i in 0..bins.len() {
        if bins[i] > 0 {
            println!(
                "{}-{}% ==> {:.1}%",
                5 * i,
                5 * (i + 1),
                100.0 * bins[i] as f64 / rates.len() as f64
            );
        }
    }
    println!(
        "mean substitution rate = {:.1}%",
        100.0 * total as f64 / rates.len() as f64
    );

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // For naive cells with junction insertion length zero, show the D gene distribution.

    let mut x = Vec::<String>::new();
    for k in 0..data.len() {
        let dref = data[k].5;
        let jun_ins = data[k].9;
        if dref == 0 && jun_ins == 0 {
            let dname = &data[k].11;
            x.push(dname.clone());
        }
    }
    x.sort();
    let mut freq = Vec::<(u32, String)>::new();
    make_freq(&x, &mut freq);
    println!(
        "\nmost frequent D genes for naive cells with junction insertion length 0 (of {})",
        x.len()
    );
    for i in 0..10 {
        println!(
            "{} [{:.1}%]",
            freq[i].1,
            100.0 * freq[i].0 as f64 / x.len() as f64
        );
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Compute jun_ins frequency for memory and naive cells.

    let mut ins = vec![vec![0; 1000]; 2];
    let mut total = vec![0; 2];
    for pass in 0..2 {
        for i in 0..data.len() {
            let dref = data[i].5;
            let jun_ins = data[i].9;
            if (pass == 0 && dref > 0) || (pass == 1 && dref == 0) {
                ins[pass][jun_ins] += 1;
                total[pass] += 1;
            }
        }
    }
    println!("\njunction insertion length frequency for memory cells");
    for i in 0..=20 {
        if ins[0][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[0][i] as f64 / total[0] as f64);
        }
    }
    println!("\njunction insertion length frequency for naive cells");
    for i in 0..=20 {
        if ins[1][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[1][i] as f64 / total[1] as f64);
        }
    }

    // Compute mean jun_ins for various classes of cells.

    let mut n_memory = 0;
    let mut n_naive = 0;
    let mut ins_memory = 0;
    let mut ins_naive = 0;
    for i in 0..data.len() {
        let dref = data[i].5;
        let jun_ins = data[i].9;
        if dref == 0 {
            n_naive += 1;
            ins_naive += jun_ins
        } else {
            n_memory += 1;
            ins_memory += jun_ins
        }
    }
    println!(
        "\nmean junction insertion bases for memory = {:.1}",
        ins_memory as f64 / n_memory as f64
    );
    println!(
        "mean junction insertion bases for naive = {:.1}",
        ins_naive as f64 / n_naive as f64
    );
    let mut n_memory = 0;
    let mut n_naive = 0;
    let mut ins_memory = 0;
    let mut ins_naive = 0;
    let mut max_ins_memory = 0;
    let mut max_ins_memory_cdr3 = String::new();
    let mut max_ins_naive = 0;
    let mut i = 0;
    let mut last = String::new();
    let mut ins_mem = vec![0; 100];
    let mut ins = vec![vec![0; 1000]; 2];
    let mut total = vec![0; 2];
    println!("");
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        let mut donors = Vec::<String>::new();
        for k in i..j {
            donors.push(data[k].3.clone());
        }
        unique_sort(&mut donors);
        if donors.len() > 1 {
            for k in i..j {
                let dref = data[k].5;
                let jun_ins = data[k].9;
                if dref == 0 {
                    n_naive += 1;
                    ins_naive += jun_ins;
                    max_ins_naive = max(max_ins_naive, jun_ins);
                    ins[1][jun_ins] += 1;
                    total[1] += 1;
                } else {
                    n_memory += 1;
                    ins_memory += jun_ins;
                    if k == i {
                        ins_mem[jun_ins] += 1;
                    }
                    if jun_ins >= 6 && last != strme(&data[k].2) {
                        println!(
                            "public memory with CDR3H = {} has jun_ins = {}",
                            stringme(&data[k].2),
                            jun_ins
                        );
                        last = stringme(&data[k].2);
                    }
                    if jun_ins > max_ins_memory {
                        max_ins_memory_cdr3 = stringme(&data[k].2);
                    }
                    max_ins_memory = max(max_ins_memory, jun_ins);
                    ins[0][jun_ins] += 1;
                    total[0] += 1;
                }
            }
        }
        i = j;
    }
    println!("\njunction insertion length frequency for public memory cells");
    for i in 0..=20 {
        if ins[0][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[0][i] as f64 / total[0] as f64);
        }
    }
    println!("\njunction insertion length frequency for public naive cells");
    for i in 0..=20 {
        if ins[1][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[1][i] as f64 / total[1] as f64);
        }
    }
    println!(
        "mean junction insertion bases for public memory = {:.1}",
        ins_memory as f64 / n_memory as f64
    );
    println!(
        "max = {} at CDR3H = {}",
        max_ins_memory, max_ins_memory_cdr3
    );
    println!(
        "mean junction insertion bases for public naive = {:.1}",
        ins_naive as f64 / n_naive as f64
    );
    println!("max = {}", max_ins_naive);
    println!("\npublic memory junction insertion length distribution, by cell group");
    for i in 0..ins_mem.len() {
        if ins_mem[i] > 0 {
            println!("{} ==> {}", i, ins_mem[i]);
        }
    }

    // Compute light chain coherence for memory cells as a function of insertion length.

    let mut n = vec![0; 100];
    let mut same = vec![0; 100];
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in k1 + 1..j {
                if data[k1].5 > 0 && data[k2].5 > 0 {
                    if data[k1].3 != data[k2].3 {
                        let ins = data[k1].9;
                        if ins == data[k2].9 {
                            n[ins] += 1;
                            if data[k1].4 == data[k2].4 {
                                same[ins] += 1;
                            }
                        }
                    }
                }
            }
        }
        i = j;
    }
    println!(
        "\nlight chain coherence for public memory cells as function of junction insertion length"
    );
    for i in 0..=4 {
        if n[i] > 0 {
            println!(
                "ins = {}  ==> coherence = {:.1}% ({} of {})",
                i,
                100.0 * same[i] as f64 / n[i] as f64,
                same[i],
                n[i],
            );
        }
    }
    let mut n_big = 0;
    let mut same_big = 0;
    for i in 5..n.len() {
        n_big += n[i];
        same_big += same[i];
    }
    println!(
        "ins >= 5 ==> coherence = {:.1}% ({} of {})",
        100.0 * same_big as f64 / n_big as f64,
        same_big,
        n_big,
    );

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, Vec<Vec<(usize, usize, usize, usize)>>)>::new();
    let mut bounds2 = Vec::<(usize, usize, Vec<usize>, Vec<usize>)>::new();
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_9(&data, i as i32) as usize;
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, vec![vec![(0, 0, 0, 0); 11]; 7]));
        bounds2.push((i, j, vec![0; 100], vec![0; 100]));
        i = j;
    }

    // Compute light chain coherence for memory cells as a function of insertion length, as
    // earlier but now using 90% CDRH3 identity.

    bounds2.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let mut samex = 0;
                for m in 0..data[k1].2.len() {
                    if data[k1].2[m] == data[k2].2[m] {
                        samex += 1;
                    }
                }
                let ident = 100.0 * samex as f64 / data[k1].2.len() as f64;
                let ident = ident.floor() as usize;
                if ident >= 90 {
                    if data[k1].5 > 0 && data[k2].5 > 0 {
                        if data[k1].3 != data[k2].3 {
                            let ins = data[k1].9;
                            if ins == data[k2].9 {
                                res.3[ins] += 1;
                                if data[k1].4 == data[k2].4 {
                                    res.2[ins] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    });
    let mut n = vec![0; 100];
    let mut same = vec![0; 100];
    for i in 0..bounds2.len() {
        for j in 0..100 {
            same[j] += bounds2[i].2[j];
            n[j] += bounds2[i].3[j];
        }
    }
    println!(
        "\nlight chain coherence for public (>= 90%) memory cells as function of junction \
        insertion length"
    );
    for i in 0..=9 {
        if n[i] > 0 {
            println!(
                "ins = {}  ==> coherence = {:.1}% ({} of {})",
                i,
                100.0 * same[i] as f64 / n[i] as f64,
                same[i],
                n[i],
            );
        }
    }
    let mut n_big = 0;
    let mut same_big = 0;
    for i in 10..n.len() {
        n_big += n[i];
        same_big += same[i];
    }
    println!(
        "ins >= 10 ==> coherence = {:.1}% ({} of {})",
        100.0 * same_big as f64 / n_big as f64,
        same_big,
        n_big,
    );
}
