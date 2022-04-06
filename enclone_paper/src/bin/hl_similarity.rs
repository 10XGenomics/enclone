// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Compute heavy chain similarity to light chain similarity, considering only memory cells from
// different donors.
//
// hl_similarity per_cell_stuff
//
// Second and third arguments: option SVG file name for plots.
//
// Generate a CSV file as output with fields as follows.
// class: 1 or 2; 2 has the additional restriction that the CDRH3 lengths are the same
// donor1: d1 or d2 or d3 or d4
// donor2: d1 or d2 or d3 or d4
// const1: const region name for first cell
// const2: const region name for first cell
// hd: heavy chain edit distance, excluding leader
// ld: light chain edit distance, excluding leader.

use enclone_tail::plot_points::plot_points;
use io_utils::*;
use pretty_trace::PrettyTrace;
use rand_chacha;
use rand_chacha::rand_core::RngCore;
use rand_chacha::rand_core::SeedableRng;
use std::collections::{HashMap, HashSet};
use std::env;
use std::io::{BufRead, Write};
use string_utils::TextUtils;
use triple_accel::levenshtein;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load data.

    let f = open_for_read![&args[1]];
    let mut svg_file1 = String::new();
    let mut svg_file2 = String::new();
    if args.len() > 3 {
        svg_file1 = args[2].to_string();
        svg_file2 = args[3].to_string();
    }
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
        String,
        usize,
        usize,
        String,
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
                first = false;
            }
        } else {
            data.push((
                /* 0 */ fields[tof["v_name1"]].to_string(),
                /* 1 */ fields[tof["cdr3_aa1"]].len(),
                /* 2 */ fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                /* 3 */ fields[tof["donors_cell"]].to_string(),
                /* 4 */ fields[tof["v_name2"]].to_string(),
                /* 5 */ fields[tof["dref"]].force_usize(),
                /* 6 */ fields[tof["clonotype_ncells"]].force_usize(),
                /* 7 */ fields[tof["const1"]].to_string(),
                /* 8 */ fields[tof["hcomp"]].force_usize(),
                /* 9 */ fields[tof["datasets_cell"]].force_usize(),
                /* 10 */ fields[tof["vj_aa_nl1"]].to_string(),
                /* 11 */ fields[tof["vj_aa_nl2"]].to_string(),
            ));
        }
    }
    data.sort();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Print CSV header.

    if svg_file1.len() == 0 {
        println!("class,donor1,donor2,const1,const2,hd,ld");
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Define groups based on equal CDRH3 amino acid sequences.

    let mut bounds = Vec::<(usize, usize)>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j));
        i = j;
    }

    // Find pairs.

    let mut bucket = Vec::<(usize, usize)>::new();
    for m in 0..bounds.len() {
        let i = bounds[m].0;
        let j = bounds[m].1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let donor1 = &data[k1].3;
                let donor2 = &data[k2].3;
                let dref1 = data[k1].5;
                let dref2 = data[k2].5;
                if donor1 != donor2 && dref1 > 0 && dref2 > 0 {
                    let _hname1 = &data[k1].0;
                    let _hname2 = &data[k2].0;
                    if
                    /* hname1 == hname2 && */
                    data[k1].2 == data[k2].2 {
                        bucket.push((k1, k2));
                    }
                }
            }
        }
    }

    // Add points.

    let mut randme = rand_chacha::ChaCha8Rng::seed_from_u64(123456789);
    let mut points2 = Vec::<(u32, (u8, u8, u8), f32, f32)>::new();
    let mut x = Vec::<f64>::new();
    let mut y = Vec::<f64>::new();

    let mut select = vec![false; bucket.len()];
    let mut count = 0;
    loop {
        if count == bucket.len() / 3 {
            break;
        }
        let m = randme.next_u64() as usize % bucket.len();
        if !select[m] {
            count += 1;
            select[m] = true;
        }
    }
    let mut lds = Vec::<usize>::new();
    for m in 0..bucket.len() {
        let i1 = bucket[m].0;
        let i2 = bucket[m].1;
        let h1 = &data[i1].10.as_bytes();
        let h2 = &data[i2].10.as_bytes();
        let hd = levenshtein(&h1, &h2) as usize;
        let l1 = &data[i1].11.as_bytes();
        let l2 = &data[i2].11.as_bytes();
        let ld = levenshtein(&l1, &l2) as usize;
        x.push(hd as f64);
        y.push(ld as f64);
        lds.push(ld);
        if !select[m] {
            continue;
        }
        let class = 2;
        let donor1 = &data[i1].3;
        let donor2 = &data[i2].3;
        let const1 = &data[i1].7;
        let const2 = &data[i2].7;
        if svg_file1.len() == 0 {
            println!("{class},{donor1},{donor2},{const1},{const2},{hd},{ld}");
        }
        points2.push((0, (255, 0, 0), hd as f32, ld as f32));
    }
    let mut small = 0;
    for i in 0..lds.len() {
        if lds[i] <= 20 {
            small += 1;
        }
    }
    println!(
        "\nsame CDR3H_AA: {:.1}% of cell pairs have ld <= 20\n",
        100.0 * small as f64 / lds.len() as f64
    );

    // Compute R^2.  Not particularly useful.

    let n = x.len();
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;
    let mut sd_x = 0.0;
    for i in 0..n {
        sd_x += (x[i] - mean_x) * (x[i] - mean_x);
    }
    sd_x = sd_x / n as f64;
    sd_x = sd_x.sqrt();
    let mut sd_y = 0.0;
    for i in 0..n {
        sd_y += (y[i] - mean_y) * (y[i] - mean_y);
    }
    sd_y = sd_y / n as f64;
    sd_y = sd_y.sqrt();
    println!("\nx has mean = {mean_x:.1}, sd = {sd_x:.1}");
    println!("y has mean = {mean_y:.1}, sd = {sd_y:.1}");
    let mut cov = 0.0;
    for i in 0..n {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    cov /= n as f64;
    let pearson_correl_coeff = cov / (sd_x * sd_y);
    let r2 = pearson_correl_coeff * pearson_correl_coeff;
    println!("R^2 = {r2:.3}\n");

    // Select pairs of cells having different donors at random.  Compute their heavy and light
    // chain edit distances.

    let mut seen = HashSet::<(usize, usize)>::new();
    let mut points1 = Vec::<(u32, (u8, u8, u8), f32, f32)>::new();
    let mut lds = Vec::<usize>::new();
    while seen.len() < count {
        let i1 = randme.next_u64() as usize % data.len();
        let i2 = randme.next_u64() as usize % data.len();
        if i1 < i2 && !seen.contains(&(i1, i2)) {
            let donor1 = &data[i1].3;
            let donor2 = &data[i2].3;
            let dref1 = data[i1].5;
            let dref2 = data[i2].5;
            if donor1 != donor2 && dref1 > 0 && dref2 > 0 {
                seen.insert((i1, i2));
                let h1 = &data[i1].10.as_bytes();
                let h2 = &data[i2].10.as_bytes();
                let hd = levenshtein(&h1, &h2) as usize;
                let l1 = &data[i1].11.as_bytes();
                let l2 = &data[i2].11.as_bytes();
                let ld = levenshtein(&l1, &l2) as usize;
                lds.push(ld);
                let class = 2;
                let const1 = &data[i1].7;
                let const2 = &data[i2].7;
                if svg_file1.len() == 0 {
                    println!("{class},{donor1},{donor2},{const1},{const2},{hd},{ld}");
                }
                points1.push((0, (0, 0, 255), hd as f32, ld as f32));
            }
        }
    }
    let mut small = 0;
    for i in 0..lds.len() {
        if lds[i] <= 20 {
            small += 1;
        }
    }
    println!(
        "\nunrestricted: {:.1}% of cell pairs have ld <= 20\n",
        100.0 * small as f64 / lds.len() as f64
    );

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Make plots.

    if svg_file1.len() > 0 {
        for pass in 1..=2 {
            let mut points = if pass == 1 {
                points1.clone()
            } else {
                points2.clone()
            };
            points.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut pointsx = Vec::new();
            let mut i = 0;
            while i < points.len() {
                let mut j = i + 1;
                while j < points.len() {
                    if points[j] != points[i] {
                        break;
                    }
                    j += 1;
                }
                let mut p = points[i].clone();
                let n = j - i;
                let r = (n as f64).sqrt();
                p.0 = r.round() as u32;
                pointsx.push(p);
                i = j;
            }
            let mut svg = String::new();
            plot_points(
                &pointsx,
                "heavy chain edit distance",
                "light chain edit distance",
                &mut svg,
                true,
                Some(String::new()),
                Some(0.0),
                Some(90.0),
                Some(0.0),
                Some(90.0),
                None,
            )
            .unwrap();
            let svg_file = if pass == 1 { &svg_file1 } else { &svg_file2 };
            let mut f = open_for_write_new![&svg_file];
            fwrite!(f, "{}", svg);
        }
    }
}
