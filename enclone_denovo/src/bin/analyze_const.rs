// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Analyze the file igh*.data.  For now, only use the first entry for a given dataset.
//
// usage: analyze_igh_data x
//        where x is a or d or e or g or m (heavy chains) or k or l (light chains)

use pretty_trace::PrettyTrace;
use std::env;
use string_utils::{abbrev_list, TextUtils};
use vector_utils::{make_freq, sort_sync2};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let chain = &args[1];
    let chainx;
    if chain == "a" || chain == "d" || chain == "e" || chain == "g" || chain == "m" {
        chainx = format!("IGH{}", chain.to_uppercase());
    } else {
        chainx = format!("IG{}C", chain.to_uppercase());
    }

    // Load data.

    let mut data = String::new();
    let data_a = include_str!("../igha.data");
    let data_d = include_str!("../ighd.data");
    let data_e = include_str!("../ighe.data");
    let data_g = include_str!("../ighg.data");
    let data_m = include_str!("../ighm.data");
    let data_k = include_str!("../igkc.data");
    let data_l = include_str!("../iglc.data");
    if chain == "a" {
        data = data_a.to_string();
    }
    if chain == "d" {
        data = data_d.to_string();
    }
    if chain == "e" {
        data = data_e.to_string();
    }
    if chain == "g" {
        data = data_g.to_string();
    }
    if chain == "m" {
        data = data_m.to_string();
    }
    if chain == "k" {
        data = data_k.to_string();
    }
    if chain == "l" {
        data = data_l.to_string();
    }

    // Parse data.

    let mut xs = Vec::<Vec<u8>>::new();
    let mut names = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in data.lines() {
        lines.push(line.to_string());
        let x = line.before(" ").as_bytes().to_vec();
        let name = line.after("  ").to_string();
        if !names.is_empty() && name == names[names.len() - 1] {
            continue;
        }
        xs.push(x);
        names.push(name);
    }

    // Collate bases.

    let n = xs[0].len();
    let mut calls = vec![Vec::<u8>::new(); n];
    for i in 0..xs.len() {
        for j in 0..n {
            let c = xs[i][j];
            if c != b'|' && c != b'.' {
                calls[j].push(c);
            }
        }
    }

    // Define motif and PWM.

    let mut motifs = Vec::<Vec<u8>>::new();
    let mut bits = 0;
    let mut nmotifs = 0;
    const MIN_FRAC: f64 = 0.9;
    let mut pwm = String::new();
    for i in 0..n {
        let k = calls[i].len();
        let mut m = Vec::<u8>::new();
        if k > 0 {
            calls[i].sort_unstable();
            let mut freqs = Vec::<(u32, u8)>::new();
            make_freq(&calls[i], &mut freqs);
            if i > 0 {
                pwm += ":";
            }
            for (u, b) in [b'A', b'C', b'G', b'T'].iter().enumerate() {
                if u > 0 {
                    pwm += ",";
                }
                let mut count = 0;
                for m in 0..freqs.len() {
                    if freqs[m].1 == *b {
                        count = freqs[m].0;
                    }
                }
                pwm += &format!("{}", count);
            }
            if freqs[0].0 as f64 / k as f64 >= MIN_FRAC {
                m.push(freqs[0].1);
                bits += 2;
                nmotifs += 1;
            } else if (freqs[0].0 + freqs[1].0) as f64 / k as f64 >= MIN_FRAC {
                m.push(freqs[0].1);
                m.push(freqs[1].1);
                bits += 1;
                nmotifs += 1;
            }
        }
        motifs.push(m);
    }
    println!("\n{}_PWM = {}", chainx, pwm);
    println!("\n{} calls = {}", chainx, lines.len());
    println!("different genomes = {}", xs.len());
    println!("bits = {}", bits);

    // Print motif.

    for i in 0..n {
        if (i + 1) % 10 == 0 && i + 1 >= 100 {
            print!("{}", (i + 1) / 100);
        } else {
            print!(" ");
        }
    }
    println!();
    for i in 0..n {
        if (i + 1) % 10 == 0 {
            print!("{}", ((i + 1) / 10) % 10);
        } else {
            print!(" ");
        }
    }
    println!();
    for i in 0..n {
        print!("{}", (i + 1) % 10);
    }
    println!();
    for pass in 0..2 {
        for i in 0..n {
            if motifs[i].len() > pass {
                print!("{}", motifs[i][pass] as char);
            } else if xs[0][i] == b'|' {
                print!("|");
            } else {
                print!(" ");
            }
        }
        println!();
    }

    // Find errors relative to motifs.

    let mut pre = 0;
    for i in 0..n {
        if xs[0][i] == b'|' {
            pre = i;
            break;
        }
    }
    println!("\npre = {}", pre);
    const ALLOWED_PRE_ERRS: usize = 3;
    println!("pre errs exceeding {}:", ALLOWED_PRE_ERRS);
    let mut all_errs = Vec::<usize>::new();
    for i in 0..xs.len() {
        let mut errs = 0;
        let mut p = 0;
        for j in 0..n {
            if !motifs[j].is_empty()
                && xs[i][j] != motifs[j][0]
                && (motifs[j].len() == 1 || xs[i][j] != motifs[j][1])
            {
                errs += 1;
                if j < pre {
                    p += 1;
                }
            }
        }
        all_errs.push(errs);
        if p > ALLOWED_PRE_ERRS {
            println!("* genome = {}", names[i]);
        }
    }
    sort_sync2(&mut all_errs, &mut names);
    println!("\nnerrs = {}", abbrev_list(&all_errs));
    println!("top errs:");
    for j in 1..=5 {
        println!(
            "* genome = {}, errs = {}",
            names[names.len() - j],
            all_errs[names.len() - j]
        );
    }

    // Assess power.

    const MIN_THRESH: f64 = 13.5;
    for e in (1..=20).rev() {
        // let z = log10( 2^bits / (nmotifs choose e) ).
        let mut z = bits as f64 * 2.0_f64.log10();
        for j in 1..=e {
            z += (j as f64).log10();
        }
        for j in nmotifs - e + 1..=nmotifs {
            z -= (j as f64).log10();
        }

        if z >= MIN_THRESH {
            println!("\nacceptable errs = {}\n", e);
            break;
        }
    }
}
