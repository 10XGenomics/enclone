// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Input = id (integer).
//
// Create a cell barcode x feature barcode matrix for the most frequent barcodes.
// 1. Extract {cell barcode, umi, feature barcode} for each read.
// 2. For a given {cell barcode, umi}, pick the most frequent feature barcode.
// 3. Report data for only the top 100 feature barcodes.
// 4. List the barcodes in order by frequency.
//
// Also compute several other statistical entities.

use chrono::prelude::*;
use enclone_core::defs::get_config;
use enclone_core::fetch_url;
use flate2::read::MultiGzDecoder;
use io_utils::{dir_list, path_exists};
use itertools::Itertools;
use mirror_sparse_matrix::MirrorSparseMatrix;
use perf_stats::*;
use rayon::prelude::*;
use serde_json::Value;
use stats_utils::*;
use std::collections::HashMap;
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::time::Instant;
use string_utils::{stringme, strme, TextUtils};
use vector_utils::*;

pub struct SequencingDef {
    pub read_path: String,
    pub sample_indices: Vec<String>,
    pub lanes: Vec<usize>,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn feature_barcode_matrix_seq_def(id: usize) -> Option<SequencingDef> {
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
        let mut sample_def = format!(
            "{}{}",
            sample_def.rev_before(","),
            sample_def.rev_after(",")
        );

        // Remove some trailing commas to get proper json.

        let mut chars = Vec::<char>::new();
        for char in sample_def.chars() {
            chars.push(char);
        }
        let mut to_delete = vec![false; chars.len()];
        for i in 0..chars.len() {
            if chars[i] == ',' {
                for j in i + 1..chars.len() {
                    if chars[j] == ' ' || chars[j] == '\n' {
                    } else if chars[j] == ']' {
                        to_delete[i] = true;
                    } else {
                        break;
                    }
                }
            }
        }
        sample_def.clear();
        for i in 0..chars.len() {
            if !to_delete[i] {
                sample_def.push(chars[i]);
            }
        }

        // Convert string to json.

        let v: Result<Value, _> = serde_json::from_str(&sample_def);
        if v.is_err() {
            println!("\nsample_def = $$${}$$$", sample_def);
        }
        let v = v.unwrap();
        let sample_defx = &v.as_array().unwrap();
        for x in sample_defx.iter() {
            if x["library_type"] == "Antibody Capture" {
                read_path = x["read_path"].to_string().between("\"", "\"").to_string();
                let y = x["sample_indices"].as_array().unwrap().to_vec();
                for i in 0..y.len() {
                    // Manually converting one case for backward compatibility.  Could use a
                    // full table if needed.
                    if y[i] == "SI-P01-E1" {
                        si.push("CGCCCGTA".to_string());
                        si.push("GTTTGCCT".to_string());
                        si.push("TAAGTTAG".to_string());
                        si.push("ACGAAAGC".to_string());
                    } else {
                        si.push(y[i].to_string());
                    }
                }
                let y = x["lanes"].as_array().unwrap().to_vec();
                for i in 0..y.len() {
                    lanes.push(y[i].to_string().force_usize());
                }
            }
        }
        for i in 0..si.len() {
            if si[i].contains("\"") {
                si[i] = si[i].between("\"", "\"").to_string();
            }
        }
        if read_path.is_empty() {
            println!("\nfailed to find read path, presumably because there's no antibody data\n");
            println!("(this is probably ok)\n");
            return None;
        }
        if !path_exists(&read_path) {
            eprintln!("\nread path does not exist");
            std::process::exit(1);
        }
    }
    if lanes.is_empty() {
        eprintln!("\nfailed to find lanes\n");
        std::process::exit(1);
    }
    let attr = fs::metadata(&read_path).unwrap();
    if !attr.is_dir() {
        eprintln!("\nread path is not a directory\n");
        std::process::exit(1);
    }

    let seq_def = SequencingDef {
        read_path,
        sample_indices: si,
        lanes,
    };
    Some(seq_def)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This computes a feature barcode matrix and several other statistical entities:
//
// 1 = MirrorSparseMatrix
// • rows = cell barcodes
// • columns = frequent feature barcodes
// • entries = number of UMIs
//
// 2 = u64 = "total UMIs" = total number of pairs (cell barcode, UMI) observed
//
// 3 = Vec<(String, u32, u32)> = for each cell barcode the number of UMIs whose feature barcode
//     is reference, and the number whose feature barcode is nonreference.
//
// 4 = MirrorSparseMatrix
// • rows = cell barcodes
// • columns = frequent feature barcodes
// • entries = number of nondegenerate reads
//
// 5 = u64 = total number of reads (meaning as usual, read pairs)
//
// 6 = Vec<(String, u32, u32)> = for each cell barcode the number of nondegenerate reads whose
//     feature barcode is reference, and the number whose feature barcode is nonreference.
//
// 7 = Vec<(String, u32, u32, u32)> = for each cell barcode the number of degenerate reads,
//     the number of those that are canonical, and the number of those that are semicanonical.

pub fn feature_barcode_matrix(
    seq_def: &SequencingDef,
    id: usize,
    verbosity: usize,
    ref_fb: &Vec<String>,
) -> Result<
    (
        MirrorSparseMatrix,
        u64,
        Vec<(String, u32, u32)>,
        MirrorSparseMatrix,
        u64,
        Vec<(String, u32, u32)>,
        Vec<(String, u32, u32, u32)>,
    ),
    String,
> {
    let t = Instant::now();

    // Find the read files.

    let mut read_files = Vec::<String>::new();
    let x = dir_list(&seq_def.read_path);
    if verbosity > 0 {
        println!();
    }
    for f in x.iter() {
        for sample_index in seq_def.sample_indices.iter() {
            for lane in seq_def.lanes.iter() {
                if f.starts_with(&format!("read-RA_si-{}_lane-00{}-", sample_index, lane)) {
                    if verbosity > 0 {
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
        eprintln!("\nreads do not exist");
        eprintln!("read path =\n{}", seq_def.read_path);
        eprintln!(
            "sample indices = {}",
            seq_def.sample_indices.iter().format(",")
        );
        eprintln!("lanes = {}", seq_def.lanes.iter().format(","));
        eprintln!("filename structure = read-RA_si-<sample index>_lane-00<lane>-\n");
        std::process::exit(1);
    }

    // Report what we found.

    if verbosity > 0 {
        println!("\nread path = {}", seq_def.read_path);
        println!("lanes = {}", seq_def.lanes.iter().format(","));
        println!(
            "sample indices = {}",
            seq_def.sample_indices.iter().format(",")
        );
        println!("used {:.1} seconds\n", elapsed(&t));
    }

    // A read pair is called degenerate if the first ten bases of R2 are GGGGGGGGGG.
    // The following sequence is the 22-base end of the // Illumina Nextera-version of
    // the R2 primer = CTGTCTCTTATACACATCTCCGAGCCCACGAGAC.  We call a read pair canonical if it
    // is degenerate and R1 contains the 22-base sequence, and semicanonical if it does not, but
    // does contain the first ten bases of it.

    let canonical = b"CACATCTCCGAGCCCACGAGAC".to_vec(); // 22 bases

    // Traverse the reads.

    println!("start parsing reads for {}", id);
    let mut buf = Vec::<Vec<u8>>::new(); // {barcode_umi_fb}
                                         // let mut buf = Vec::<(Vec<u8>, Vec<u8>, Vec<u8>)>::new(); // {(barcode, umi, fb)}
    let mut bdcs = Vec::<(String, u32, u32, u32)>::new();
    let (mut ncanonical, mut nsemicanonical) = (0, 0);
    let mut ndegen = 0;
    for pass in 1..=2 {
        let mut degen = Vec::<(Vec<u8>, Vec<u8>)>::new(); // {(barcode, umi)}
        if verbosity > 0 {
            println!("megapass {}", pass);
        }
        for (i, rf) in read_files.iter().enumerate() {
            let local: DateTime<Local> = Local::now();
            let local = format!("{:?}", local);
            let time = local.between("T", ".");
            if verbosity > 0 {
                println!(
                    "- at {}, processing dataset {} of {}; buf has size {} and degen has size {}",
                    time,
                    i + 1,
                    read_files.len(),
                    buf.len(),
                    degen.len()
                );
            }
            let f = format!("{}/{}", seq_def.read_path, rf);
            let gz = MultiGzDecoder::new(File::open(&f).unwrap());
            let b = BufReader::new(gz);

            // Paired reads are in groups of eight lines.  Line 2 is the cell barcode-umi read,
            // and line 6 is the read that contains the feature barcode.

            let mut count = 0;
            let mut barcode = Vec::<u8>::new();
            let mut umi = Vec::<u8>::new();
            let mut read1 = Vec::<u8>::new();
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
                        read1 = s.to_vec();
                    } else {
                        fb = s[10..25].to_vec();
                        let mut degenerate = true;
                        for i in 0..10 {
                            if s[i] != b'G' {
                                degenerate = false;
                                break;
                            }
                        }

                        // Save.

                        if verbosity == 2 {
                            let (mut is_canonical, mut is_semicanonical) = (false, false);
                            if degenerate {
                                for j in 0..=read1.len() - canonical.len() {
                                    if read1[j..j + canonical.len()] == canonical {
                                        is_canonical = true;
                                        break;
                                    }
                                }
                            }
                            if degenerate && !is_canonical {
                                for i in 0..18 {
                                    let mut w = true;
                                    for j in 0..10 {
                                        if read1[i + j] != canonical[j] {
                                            w = false;
                                            break;
                                        }
                                    }
                                    if w {
                                        is_semicanonical = true;
                                        break;
                                    }
                                }
                            }
                            print!(
                                "r: {} {} {} {} {} {}",
                                strme(&barcode),
                                strme(&umi),
                                strme(&s[0..10]),
                                strme(&fb),
                                strme(&s[25..35]),
                                strme(&s[35..55]),
                            );
                            if is_canonical {
                                print!(" = canon");
                            } else if is_semicanonical {
                                print!(" = semi");
                            }
                            println!("");
                        }
                        if degenerate && pass == 1 {
                            degen.push((barcode.clone(), umi.clone()));
                        } else if pass == 2 {
                            let mut x = barcode.clone();
                            x.append(&mut umi.clone());
                            x.append(&mut fb.clone());
                            buf.push(x);
                        }
                    }
                }
            }
        }

        // Build data structure for the degenerate reads.

        if pass == 1 {
            if verbosity > 0 {
                println!("parallel sorting");
            }
            degen.par_sort();
            if verbosity > 0 {
                println!("build data structure for degenerate reads");
            }
            let mut i = 0;
            while i < degen.len() {
                let j = next_diff1_2(&degen, i as i32) as usize;
                let (mut canon, mut semi) = (0, 0);
                for k in i..j {
                    let mut read1 = degen[k].0.clone();
                    read1.append(&mut degen[k].1.clone());
                    let mut is_canonical = false;
                    for j in 0..=read1.len() - canonical.len() {
                        if read1[j..j + canonical.len()] == canonical {
                            is_canonical = true;
                            canon += 1;
                            ncanonical += 1;
                            break;
                        }
                    }
                    if !is_canonical {
                        for i in 0..18 {
                            let mut w = true;
                            for j in 0..10 {
                                if read1[i + j] != canonical[j] {
                                    w = false;
                                    break;
                                }
                            }
                            if w {
                                semi += 1;
                                nsemicanonical += 1;
                                break;
                            }
                        }
                    }
                }
                bdcs.push((stringme(&degen[i].0), (j - i) as u32, canon, semi));
                i = j;
            }
            ndegen = degen.len();
            drop(degen);
        }
    }

    // Print summary stats.

    let total_reads = ndegen + buf.len();
    if verbosity > 0 {
        println!("there are {} read pairs", total_reads);
        println!(
            "of which {:.1}% are degenerate",
            percent_ratio(ndegen, total_reads)
        );
        let canonical_percent = 100.0 * ncanonical as f64 / total_reads as f64;
        println!("canonical fraction = {:.1}%", canonical_percent);
        let semicanonical_percent = 100.0 * nsemicanonical as f64 / total_reads as f64;
        println!("semicanonical fraction = {:.1}%", semicanonical_percent);
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }

    // Find the most common feature barcodes, by read count.

    let mut fbs = Vec::<Vec<u8>>::new();
    for i in 0..buf.len() {
        fbs.push(buf[i][28..43].to_vec());
    }
    fbs.par_sort();
    let mut fb_freq = Vec::<(u32, Vec<u8>)>::new();
    make_freq(&fbs, &mut fb_freq);
    drop(fbs);
    const TOP_FEATURE_BARCODES: usize = 100;
    if fb_freq.len() > TOP_FEATURE_BARCODES {
        fb_freq.truncate(TOP_FEATURE_BARCODES);
    }
    let mut fb_freq_sorted = Vec::<Vec<u8>>::new();
    for i in 0..fb_freq.len() {
        fb_freq_sorted.push(fb_freq[i].1.clone());
    }
    let mut ids_ffs = Vec::<usize>::new();
    for i in 0..fb_freq_sorted.len() {
        ids_ffs.push(i);
    }
    sort_sync2(&mut fb_freq_sorted, &mut ids_ffs);

    // Generate a feature-barcode matrix for the common feature barcodes, by read count.  Also,
    // make a table of all barcodes, and for each, the number of reference and nonreference reads.

    if verbosity > 0 {
        println!("making mirror sparse matrix, by read counts");
    }
    let mut bf = Vec::<(Vec<u8>, Vec<u8>)>::new();
    for i in 0..buf.len() {
        bf.push((buf[i][0..16].to_vec(), buf[i][28..43].to_vec()));
    }
    bf.par_sort();
    let mut brnr = Vec::<(String, u32, u32)>::new();
    let mut x = Vec::<Vec<(i32, i32)>>::new();
    let mut row_labels = Vec::<String>::new();
    let mut col_labels = Vec::<String>::new();
    for z in fb_freq.iter() {
        col_labels.push(stringme(&z.1));
    }
    if verbosity > 0 {
        println!("start making m_reads, used {:.1} seconds", elapsed(&t));
    }
    let mut i = 0;
    while i < bf.len() {
        let j = next_diff1_2(&bf, i as i32) as usize;
        let mut y = Vec::<(i32, i32)>::new();
        let mut k = i;
        while k < j {
            let l = next_diff(&bf, k);
            let p = bin_position(&fb_freq_sorted, &bf[k].1);
            if p >= 0 {
                y.push((ids_ffs[p as usize] as i32, (l - k) as i32));
            }
            k = l;
        }
        if !y.is_empty() {
            row_labels.push(format!("{}-1", strme(&bf[i].0)));
            x.push(y);
        }
        let (mut refx, mut nrefx) = (0, 0);
        for k in i..j {
            if bin_member(ref_fb, &stringme(&bf[k].1)) {
                refx += 1;
            } else {
                nrefx += 1;
            }
        }
        brnr.push((stringme(&bf[i].0), refx, nrefx));
        i = j;
    }
    drop(bf);
    let m_reads = MirrorSparseMatrix::build_from_vec(&x, &row_labels, &col_labels);
    if verbosity > 0 {
        println!("after making m_reads, used {:.1} seconds", elapsed(&t));
    }

    // Sort buf.

    buf.par_sort();

    // Reduce buf to UMI counts.

    let mut bfu = Vec::<Vec<u8>>::new(); // {barcode_fb_umi}
    let mut singletons = 0;
    let mut i = 0;
    while i < buf.len() {
        let mut j = i + 1;
        while j < buf.len() {
            if buf[j][0..28] != buf[i][0..28] {
                break;
            }
            j += 1;
        }
        let mut sing = true;
        if i > 0 && buf[i][0..16] == buf[i - 1][0..16] {
            sing = false;
        }
        if i < buf.len() - 1 && buf[i][0..16] == buf[i + 1][0..16] {
            sing = false;
        }
        if sing {
            singletons += 1;
        }
        let mut bfs = Vec::<String>::new();
        for k in i..j {
            bfs.push(stringme(&buf[k][28..43].to_vec()));
        }
        let mut uniq = true;
        for k in i + 1..j {
            if buf[k][28..43] != buf[k - 1][28..43] {
                uniq = false;
            }
        }
        if uniq {
            let mut x = buf[i][0..16].to_vec();
            x.append(&mut buf[i][28..43].to_vec());
            x.append(&mut buf[i][16..28].to_vec());
            bfu.push(x);
        } else {
            // Sloppy, just take the most frequent feature barcode, ignoring quality scores.
            let mut fs = Vec::<Vec<u8>>::new();
            for k in i..j {
                fs.push(buf[k][28..43].to_vec());
            }
            fs.sort();
            let mut freq = Vec::<(u32, Vec<u8>)>::new();
            make_freq(&fs, &mut freq);
            let mut x = buf[i][0..16].to_vec();
            x.append(&mut freq[0].1.clone());
            x.append(&mut buf[i][16..28].to_vec());
            bfu.push(x);
        }
        i = j;
    }
    drop(buf);
    bfu.par_sort();

    // Make table of all barcodes, and for each, the number of reference and nonreference UMIs.

    let mut brn = Vec::<(String, u32, u32)>::new();
    let mut i = 0;
    while i < bfu.len() {
        let mut j = i + 1;
        while j < bfu.len() {
            if bfu[j][0..16] != bfu[i][0..16] {
                break;
            }
            j += 1;
        }
        let mut refx = 0;
        let mut nrefx = 0;
        for k in i..j {
            if bin_member(ref_fb, &stringme(&bfu[k][16..31])) {
                refx += 1;
            } else {
                nrefx += 1;
            }
        }
        brn.push((stringme(&bfu[i][0..16]), refx, nrefx));
        i = j;
    }

    // Proceed.

    if verbosity > 0 {
        let singleton_percent = 100.0 * singletons as f64 / total_reads as f64;
        println!("singleton fraction = {:.1}%", singleton_percent);
        let mut bc = 0;
        let mut i = 0;
        while i < bfu.len() {
            let mut j = i + 1;
            while j < bfu.len() {
                if bfu[j][0..16] != bfu[i][0..16] {
                    break;
                }
                j += 1;
            }
            bc += 1;
            i = j;
        }
        println!("total barcodes = {}", bc);
    }
    let mut bfn = Vec::<(Vec<u8>, Vec<u8>, usize)>::new(); // {(barcode, fb, numis)}
    let mut i = 0;
    while i < bfu.len() {
        let mut j = i + 1;
        while j < bfu.len() {
            if bfu[j][0..31] != bfu[i][0..31] {
                break;
            }
            j += 1;
        }
        let mut n = 1;
        for k in i + 1..j {
            if bfu[k][31..43] != bfu[k - 1][31..43] {
                n += 1;
            }
        }
        bfn.push((bfu[i][0..16].to_vec(), bfu[i][16..31].to_vec(), n));
        i = j;
    }
    if verbosity > 0 {
        println!("there are {} uniques", bfu.len());
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }
    let mut total_umis = 0;
    for i in 0..bfn.len() {
        total_umis += bfn[i].2 as u64;
    }
    if verbosity > 0 {
        println!("total UMIs = {}\n", total_umis);
    }

    // Report common feature barcodes, by UMI count.

    if verbosity > 0 {
        println!("common feature barcodes and their total UMI counts\n");
    }
    let mut fbx = Vec::<Vec<u8>>::new();
    for i in 0..bfu.len() {
        fbx.push(bfu[i][16..31].to_vec());
    }
    drop(bfu);
    if verbosity > 0 {
        println!("parallel sorting fbx");
    }
    fbx.par_sort();
    let mut freq = Vec::<(u32, Vec<u8>)>::new();
    make_freq(&fbx, &mut freq);
    drop(fbx);
    if verbosity > 0 {
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
    if verbosity > 0 {
        println!("\nused {:.1} seconds\n", elapsed(&t));
    }

    // Generate a feature-barcode matrix for the common feature barcodes, by UMI count.

    if verbosity > 0 {
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
    drop(bfn);
    if verbosity > 0 {
        println!("building last mirror sparse matrix");
    }
    let m = MirrorSparseMatrix::build_from_vec(&x, &row_labels, &col_labels);
    if verbosity > 0 {
        println!("used {:.1} seconds\n", elapsed(&t));
        #[cfg(not(target_os = "windows"))]
        {
            println!(
                "peak mem usage for feature barcode matrix = {:.1} GB\n",
                peak_mem_usage_gb()
            );
        }
    }
    Ok((m, total_umis, brn, m_reads, total_reads as u64, brnr, bdcs))
}
