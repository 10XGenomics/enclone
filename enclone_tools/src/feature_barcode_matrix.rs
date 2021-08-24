// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Input = id (integer).
//
// Create a cell barcode x feature barcode matrix for the most frequent barcodes.
// 1. Extract {cell barcode, umi, feature barcode} for each read.
// 2. Discard feature barcode = GGGGGGGGGGGGGGG.
// 3. Discard data where a cell barcode has only one read.
// 4. For a given {cell barcode, umi}, pick the most frequent feature barcode.
// 5. Report data for only the top 100 feature barcodes.
// 6. List the barcodes in order by frequency.

use enclone_core::defs::get_config;
use enclone_core::*;
use flate2::read::MultiGzDecoder;
use io_utils::*;
use itertools::Itertools;
use mirror_sparse_matrix::*;
use perf_stats::*;
use rayon::prelude::*;
use serde_json::Value;
use std::collections::HashMap;
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

pub struct SequencingDef {
    read_path: String,
    sample_indices: Vec<String>,
    lanes: Vec<usize>,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn feature_barcode_matrix_head1(id: usize) -> SequencingDef {
    println!("getting feature barcode matrix for {}", id);

    // Get configuration.

    let mut config = HashMap::<String, String>::new();
    let mut config_file = String::new();
    for (key, value) in env::vars() {
        if key == "ENCLONE_CONFIG" {
            config_file = value.to_string();
            if config_file.contains(',') {
                config_file = config_file.after(",").to_string();
            }
        }
    }
    let _ = get_config(&config_file, &mut config);

    // Get pipestance path.

    let url = format!("{}/{}", config["ones"], id);
    let m = fetch_url(&url).unwrap();
    if m.contains("502 Bad Gateway") {
        eprintln!(
            "\nWell, this is sad.  The URL \
            {} returned a 502 Bad Gateway \
            message.  Please try again later or ask someone for help.\n",
            url
        );
        std::process::exit(1);
    }
    let v: Value = serde_json::from_str(&m).unwrap();
    let latest = &v["latest_analysis_run"];
    let pipestance = latest["path"].to_string().between("\"", "\"").to_string();

    // Get read path and lanes and sample indices from the invocation file.

    let invocation = format!("{}/_invocation", pipestance);
    if !path_exists(&invocation) {
        eprintln!("\n_invocation does not exist for {}\n", pipestance);
        std::process::exit(1);
    }
    let mut read_path = String::new(); // path to reads
    let mut si = Vec::<String>::new(); // sample indices
    let mut lanes = Vec::<usize>::new(); // lanes
    {
        let mut f = File::open(&invocation).unwrap();
        let mut bytes = Vec::<u8>::new();
        f.read_to_end(&mut bytes).unwrap();
        let mut sample_def = Vec::<u8>::new();
        let mut bracks: isize = 0;
        for i in 0..bytes.len() {
            if bytes[i] == b'[' {
                bracks += 1;
            } else if bytes[i] == b']' {
                bracks -= 1;
            }
            if bracks > 0 {
                sample_def.push(bytes[i]);
            } else if bracks == 0 && !sample_def.is_empty() {
                sample_def.append(&mut b"]\n".to_vec());
                break;
            }
        }
        let mut sample_def = stringme(&sample_def);
        sample_def = sample_def.replace(",\n        }", "\n        }");
        let sample_def = format!(
            "{}{}",
            sample_def.rev_before(","),
            sample_def.rev_after(",")
        );
        let v: Value = serde_json::from_str(&sample_def).unwrap();
        let sample_def = &v.as_array().unwrap();
        for x in sample_def.iter() {
            if x["library_type"] == "Antibody Capture" {
                read_path = x["read_path"].to_string().between("\"", "\"").to_string();
                let y = x["sample_indices"].as_array().unwrap().to_vec();
                for i in 0..y.len() {
                    si.push(y[i].to_string());
                }
                let y = x["lanes"].as_array().unwrap().to_vec();
                for i in 0..y.len() {
                    lanes.push(y[i].to_string().force_usize());
                }
            }
        }
        for i in 0..si.len() {
            si[i] = si[i].between("\"", "\"").to_string();
        }
        if read_path.len() == 0 {
            eprintln!("\nfailed to find read path\n");
            std::process::exit(1);
        }
        if !path_exists(&read_path) {
            eprintln!("\nread path does not exist");
            std::process::exit(1);
        }
    }
    if lanes.len() == 0 {
        eprintln!("\nfailed to find lanes\n");
        std::process::exit(1);
    }
    let attr = fs::metadata(&read_path).unwrap();
    if !attr.is_dir() {
        eprintln!("\nread path is not a directory\n");
        std::process::exit(1);
    }

    let seq_def = SequencingDef {
        read_path: read_path,
        sample_indices: si,
        lanes: lanes,
    };
    seq_def
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn feature_barcode_matrix(id: usize, verbose: bool) -> Result<MirrorSparseMatrix, String> {
    let t = Instant::now();
    let seq_def = feature_barcode_matrix_head1(id);

    // Find the read files.

    let mut read_files = Vec::<String>::new();
    let x = dir_list(&seq_def.read_path);
    if verbose {
        println!("");
    }
    for f in x.iter() {
        for sample_index in seq_def.sample_indices.iter() {
            for lane in seq_def.lanes.iter() {
                if f.contains(&format!("read-RA_si-{}_lane-00{}-", sample_index, lane)) {
                    if verbose {
                        println!("{}", f);
                    }
                    if f.ends_with("__evap") {
                        eprintln!(
                            "\nfound an evaporated read file =\n{}\n",
                            format!("{}/{}", seq_def.read_path, f),
                        );
                        std::process::exit(1);
                    }
                    read_files.push(f.clone());
                }
            }
        }
    }
    if read_files.is_empty() {
        eprintln!("\nreads do not exist\n");
        std::process::exit(1);
    }

    // Report what we found.

    if verbose {
        println!("\nread path = {}", seq_def.read_path);
        println!("lanes = {}", seq_def.lanes.iter().format(","));
        println!("sample indices = {}", seq_def.sample_indices.iter().format(","));
        println!("used {:.1} seconds\n", elapsed(&t));
    }

    // Traverse the reads.

    eprintln!("start parsing reads for {}", id);
    let mut buf = Vec::<(Vec<u8>, Vec<u8>, Vec<u8>)>::new(); // {(barcode, umi, fb)}
    let mut total_reads = 0;
    let mut junk = 0;
    for rf in read_files.iter() {
        let f = format!("{}/{}", seq_def.read_path, rf);
        let gz = MultiGzDecoder::new(File::open(&f).unwrap());
        let b = BufReader::new(gz);

        // Paired reads are in groups of eight lines.  Line 2 is the cell barcode-umi read,
        // and line 6 is the read that contains the feature barcode.

        let mut count = 0;
        let mut barcode = Vec::<u8>::new();
        let mut umi = Vec::<u8>::new();
        let mut fb;
        for line in b.lines() {
            count += 1;
            if count % 8 == 2 || count % 8 == 6 {
                let s = line.unwrap();
                let s = s.as_bytes();
                if count % 8 == 2 {
                    if s.len() < 28 {
                        return Err(format!(
                            "\nencountered read of length {} < 28 in {}\n",
                            s.len(),
                            f
                        ));
                    }
                    assert!(s.len() >= 28);
                    barcode = s[0..16].to_vec();
                    umi = s[16..28].to_vec();
                } else {
                    fb = s[10..25].to_vec();
                    total_reads += 1;
                    if fb == b"GGGGGGGGGGGGGGG" {
                        junk += 1;
                    } else {
                        buf.push((barcode.clone(), umi.clone(), fb.clone()));
                    }
                }
            }
        }
    }
    if verbose {
        println!("there are {} read pairs", buf.len());
        let junk_percent = 100.0 * junk as f64 / total_reads as f64;
        println!("GGGGGGGGGGGGGGG fraction = {:.1}%", junk_percent);
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }

    // Reduce to UMI counts.

    buf.par_sort();
    let mut bfu = Vec::<(Vec<u8>, Vec<u8>, Vec<u8>)>::new(); // {(barcode, fb, umi)}
    let mut singletons = 0;
    let mut i = 0;
    while i < buf.len() {
        let j = next_diff12_3(&buf, i as i32) as usize;
        let mut sing = true;
        if i > 0 && buf[i].0 == buf[i - 1].0 {
            sing = false;
        }
        if i < buf.len() - 1 && buf[i].0 == buf[i + 1].0 {
            sing = false;
        }
        if sing {
            singletons += 1;
            i = j;
            continue;
        }
        let mut bfs = Vec::<String>::new();
        for k in i..j {
            bfs.push(stringme(&buf[k].2));
        }
        let mut uniq = true;
        for k in i + 1..j {
            if buf[k].2 != buf[k - 1].2 {
                uniq = false;
            }
        }
        if uniq {
            bfu.push((buf[i].0.clone(), buf[i].2.clone(), buf[i].1.clone()));
        } else {
            // Sloppy, just take the most frequent feature barcode, ignoring quality scores.
            let mut fs = Vec::<Vec<u8>>::new();
            for k in i..j {
                fs.push(buf[k].2.clone());
            }
            fs.sort();
            let mut freq = Vec::<(u32, Vec<u8>)>::new();
            make_freq(&fs, &mut freq);
            bfu.push((buf[i].0.clone(), freq[0].1.clone(), buf[i].1.clone()));
        }
        i = j;
    }
    bfu.par_sort();
    if verbose {
        let singleton_percent = 100.0 * singletons as f64 / total_reads as f64;
        println!("singleton fraction = {:.1}%", singleton_percent);
    }
    let mut bfn = Vec::<(Vec<u8>, Vec<u8>, usize)>::new(); // {(barcode, fb, numis)}
    let mut i = 0;
    while i < bfu.len() {
        let j = next_diff12_3(&bfu, i as i32) as usize;
        let mut n = 1;
        for k in i + 1..j {
            if bfu[k].2 != bfu[k - 1].2 {
                n += 1;
            }
        }
        bfn.push((bfu[i].0.clone(), bfu[i].1.clone(), n));
        i = j;
    }
    if verbose {
        println!("there are {} uniques", bfu.len());
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }

    // Report common feature barcodes.

    const TOP_FEATURE_BARCODES: usize = 100;
    if verbose {
        println!("common feature barcodes\n");
    }
    let mut fbx = Vec::<Vec<u8>>::new();
    for i in 0..bfu.len() {
        fbx.push(bfu[i].1.clone());
    }
    fbx.par_sort();
    let mut freq = Vec::<(u32, Vec<u8>)>::new();
    make_freq(&fbx, &mut freq);
    if verbose {
        for i in 0..10 {
            println!("{} = {}", strme(&freq[i].1), freq[i].0);
        }
    }
    let mut tops = Vec::<Vec<u8>>::new();
    for i in 0..std::cmp::min(TOP_FEATURE_BARCODES, freq.len()) {
        tops.push(freq[i].1.clone());
    }
    let mut tops_sorted = tops.clone();
    let mut ids = Vec::<usize>::new();
    for i in 0..tops.len() {
        ids.push(i);
    }
    sort_sync2(&mut tops_sorted, &mut ids);
    if verbose {
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }

    // Generate a feature-barcode matrix for the common feature barcodes.

    if verbose {
        println!("making mirror sparse matrix");
    }
    let mut x = Vec::<Vec<(i32, i32)>>::new();
    let mut row_labels = Vec::<String>::new();
    let mut col_labels = Vec::<String>::new();
    for x in tops.iter() {
        col_labels.push(stringme(x));
    }
    let mut i = 0;
    while i < bfn.len() {
        let j = next_diff1_3(&bfn, i as i32) as usize;
        let mut y = Vec::<(i32, i32)>::new();
        for k in i..j {
            let p = bin_position(&tops_sorted, &bfn[k].1);
            if p >= 0 {
                y.push((ids[p as usize] as i32, bfn[k].2 as i32));
            }
        }
        if !y.is_empty() {
            row_labels.push(format!("{}-1", strme(&bfn[i].0)));
            x.push(y);
        }
        i = j;
    }
    let m = MirrorSparseMatrix::build_from_vec(&x, &row_labels, &col_labels);
    if verbose {
        println!("used {:.1} seconds\n", elapsed(&t));
    }
    Ok(m)
}
