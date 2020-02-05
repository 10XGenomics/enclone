// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use vdj_ann::*;

use self::annotate::*;
use self::refx::*;
use self::transcript::*;
use crate::defs::*;
use crate::explore::*;
use debruijn::dna_string::*;
use io_utils::*;
use perf_stats::*;
use rayon::prelude::*;
use serde_json::Value;
use std::{collections::HashMap, io::BufReader, time::Instant};
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn json_error(json: Option<&str>) {
    eprint!(
        "\nThere is something wrong with the contig annotations in the Cell Ranger output \
         file"
    );
    if json.is_some() {
        eprint!("\n{}.", json.unwrap());
    } else {
        eprint!(".");
    }
    eprintln!(
        "\n\nHere are possible sources of this problem:\n\n\
         1. If the file was generated using \
         Cell Ranger version < 3.1, please either\nregenerate the file using the \
         current version, or else run this program with the RE option to\n\
         regenerate annotations from scratch, but we warn you that this code \
         is not guaranteed to run\ncorrectly on outdated json files.\n\n\
         2. Make sure you have the correct chain type, TCR or BCR.\n\n\
         3. Make sure you have the correct reference sequence.  See \
         \"enclone help faq\".\n\n\
         4. If none of these apply, please report the problem to \
         enclone@10xgenomics.com.  But please\nfirst rerun with RE to confirm the problem.\n"
    );
    std::process::exit(1);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the json annotations file.
//
// Should try using:
// https://martian-lang.github.io/martian-rust/doc/martian_filetypes/json_file/
// index.html#lazy-readwrite-example.
// Tried at one point but was unable to compile.
//
// Tracking contigs using bc_cdr3_aa, which is a mess, should improve later.
//
// This section requires 3.1.  If you want to avoid that, do something to make tig_start
// and tig_stop always nonnegative.  Or use the RE option.
//
// Computational performance.  It would appear that nearly all the time here is spent in
// two lines:
//
// read_vector_entry_from_json(&mut f) {
// let v: Value = serde_json::from_str(strme(&x)).unwrap();
// (Should retest.)
//
// and simply reading the file lines is several times faster.  So the way we parse the
// files is suboptimal.  If we want to make this faster, one option would be to speed up
// this code.  Another would be to write out a binary version of the json file that contains
// only the information that we need.

pub fn read_json(
    li: usize,
    json: &String,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    reannotate: bool,
) -> Vec<Vec<TigData>> {
    let mut tigs = Vec::<TigData>::new();
    let mut jsonx = json.clone();
    if !path_exists(&json) {
        jsonx = format!("{}.lz4", json);
    }
    if jsonx.contains('/') {
        let p = jsonx.rev_before("/");
        if !path_exists(&p) {
            eprintln!(
                "\nThere should be a directory\n\
                 \"{}\"\n\
                 but it does not exist.  Please check how you have specified the\n\
                 input files to enclone, including the PRE argument.\n",
                p
            );
            std::process::exit(1);
        }
    }
    if !path_exists(&jsonx) {
        eprintln!(
            "\nThe path\n\
             \"{}\"\n\
             does not exist.  Please check how you have specified the\n\
             input files to enclone, including the PRE argument.\n",
            jsonx
        );
        std::process::exit(1);
    }
    let mut f = BufReader::new(open_maybe_compressed(&jsonx));
    // ◼ This loop could be speeded up, see comments above.
    loop {
        match read_vector_entry_from_json(&mut f) {
            None => break,
            Some(x) => {
                let v: Value = serde_json::from_str(strme(&x)).unwrap();
                if !v["productive"].as_bool().unwrap_or(false) {
                    continue;
                }
                if !v["high_confidence"].as_bool().unwrap_or(false) {
                    continue;
                }
                if !v["is_cell"].as_bool().unwrap_or(false) {
                    continue;
                }
                let barcode = &v["barcode"].to_string().between("\"", "\"").to_string();
                let tigname = &v["contig_name"].to_string().between("\"", "\"").to_string();
                let full_seq = &v["sequence"].to_string().between("\"", "\"").to_string();
                let mut left = false;
                let (mut v_ref_id, mut j_ref_id) = (1000000, 0);
                let mut d_ref_id: Option<usize> = None;
                let mut c_ref_id = None;
                let mut chain_type = String::new();
                let mut u_ref_id = None;
                let (mut tig_start, mut tig_stop) = (-1 as isize, -1 as isize);
                let mut v_stop = 0;
                let mut v_stop_ref = 0;
                let mut j_start = 0;
                let mut j_start_ref = 0;
                let mut c_start = None;
                let mut annv = Vec::<(i32, i32, i32, i32, i32)>::new();
                let cdr3_aa: String;
                let cdr3_dna: String;
                let mut cdr3_start: usize;

                // Reannotate.

                if reannotate {
                    let x = DnaString::from_dna_string(&full_seq);
                    let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
                    annotate_seq(&x, &refdata, &mut ann, true, false, true);
                    let mut log = Vec::<u8>::new();
                    if !is_valid(&x, &refdata, &ann, false, &mut log) {
                        continue;
                    }
                    let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
                    get_cdr3_using_ann(&x, &refdata, &ann, &mut cdr3);
                    cdr3_aa = stringme(&cdr3[0].1);
                    cdr3_start = cdr3[0].0;
                    cdr3_dna = x
                        .slice(cdr3_start, cdr3_start + 3 * cdr3_aa.len())
                        .to_string();
                    let mut seen_j = false;
                    for i in 0..ann.len() {
                        let t = ann[i].2 as usize;
                        if refdata.is_u(t) {
                            u_ref_id = Some(t);
                        } else if refdata.is_v(t) && !seen_j {
                            v_ref_id = t;
                            annv.push(ann[i].clone());
                            chain_type = refdata.name[t][0..3].to_string();
                            if chain_type == "IGH".to_string() || chain_type == "TRB".to_string() {
                                left = true;
                            }
                            if ann[i].3 == 0 {
                                tig_start = ann[i].0 as isize;
                                cdr3_start -= tig_start as usize;
                            }
                            v_stop = (ann[i].0 + ann[i].1) as usize;
                            v_stop_ref = (ann[i].3 + ann[i].1) as usize;
                        } else if refdata.is_d(t) {
                            d_ref_id = Some(t);
                        } else if refdata.is_j(t) {
                            j_ref_id = t;
                            tig_stop = (ann[i].0 + ann[i].1) as isize;
                            j_start = ann[i].0 as usize;
                            j_start_ref = ann[i].3 as usize;
                            seen_j = true;
                        } else if refdata.is_c(t) {
                            c_ref_id = Some(t);
                            c_start = Some(ann[i].0 as usize);
                        }
                    }
                    for i in (0..annv.len()).rev() {
                        annv[i].0 -= annv[0].0;
                    }
                } else {
                    // Use annotations from json file.

                    cdr3_aa = v["cdr3"].to_string().between("\"", "\"").to_string();
                    cdr3_dna = v["cdr3_seq"].to_string().between("\"", "\"").to_string();
                    cdr3_start = v["cdr3_start"].as_u64().unwrap() as usize;
                    let ann = v["annotations"].as_array().unwrap();
                    let mut cigarv = String::new(); // cigar for V segment
                    for i in 0..ann.len() {
                        let a = &ann[i];
                        let region_type = &a["feature"]["region_type"];
                        let feature_id = a["feature"]["feature_id"].as_u64().unwrap() as usize;
                        if !to_ref_index.contains_key(&feature_id) {
                            continue;
                        }
                        let feature_id = to_ref_index[&feature_id];
                        let ref_start = a["annotation_match_start"].as_u64().unwrap() as usize;
                        if region_type == "L-REGION+V-REGION" {
                            v_stop = a["contig_match_end"].as_i64().unwrap() as usize;
                            v_stop_ref = a["annotation_match_end"].as_i64().unwrap() as usize;
                        }
                        if region_type == "L-REGION+V-REGION" && ref_start == 0 {
                            let chain = a["feature"]["chain"]
                                .to_string()
                                .between("\"", "\"")
                                .to_string();
                            // if !chain.starts_with("IG") { continue; } // *******************
                            tig_start = a["contig_match_start"].as_i64().unwrap() as isize;
                            cdr3_start -= tig_start as usize;
                            chain_type = chain.clone();
                            if chain == "IGH".to_string() || chain == "TRB".to_string() {
                                left = true;
                            }
                            v_ref_id = feature_id;
                            cigarv = a["cigar"].to_string().between("\"", "\"").to_string();
                        } else {
                            // also check for IG chain?????????????????????????????????????????
                            let ref_stop = a["annotation_match_end"].as_u64().unwrap() as usize;
                            let ref_len = a["annotation_length"].as_u64().unwrap() as usize;
                            if region_type == "J-REGION" && ref_stop == ref_len {
                                tig_stop = a["contig_match_end"].as_i64().unwrap() as isize;
                                j_ref_id = feature_id;
                                j_start = a["contig_match_start"].as_i64().unwrap() as usize;
                                j_start_ref =
                                    a["annotation_match_start"].as_i64().unwrap() as usize;
                            }
                            if region_type == "5'UTR" {
                                u_ref_id = Some(feature_id);
                            }
                            if region_type == "D-REGION" {
                                d_ref_id = Some(feature_id);
                            }
                            if region_type == "C-REGION" {
                                c_ref_id = Some(feature_id);
                                c_start = Some(a["contig_match_start"].as_i64().unwrap() as usize);
                            }
                        }
                    }
                    if v_ref_id == 1000000 {
                        continue;
                    }

                    // Compute annv from cigarv.  We don't compute the mismatch entry.

                    let mut cg = Vec::<Vec<u8>>::new(); // pieces of cigar string
                    let mut piece = Vec::<u8>::new();
                    for c in cigarv.chars() {
                        piece.push(c as u8);
                        if c.is_ascii_alphabetic() {
                            cg.push(piece.clone());
                            piece.clear();
                        }
                    }
                    let t = v_ref_id as i32;
                    let (mut len1, mut len2) = (0, 0);
                    let (mut ins, mut del) = (0, 0);
                    for i in 0..cg.len() {
                        let x = strme(&cg[i][0..cg[i].len() - 1]).force_i32();
                        if cg[i][cg[i].len() - 1] == b'M' {
                            if len1 == 0 {
                                len1 = x;
                            } else if len2 == 0 {
                                len2 = x;
                            } else {
                                // probably can't happen
                                len1 = 0;
                                len2 = 0;
                                break;
                            }
                        }
                        if cg[i][cg[i].len() - 1] == b'I' {
                            ins = x;
                        }
                        if cg[i][cg[i].len() - 1] == b'D' {
                            del = x;
                        }
                    }
                    annv.push((0 as i32, len1, t, 0, 0));
                    if ins > 0 && ins % 3 == 0 && del == 0 && len2 > 0 {
                        let start = (len1 + ins) as i32;
                        annv.push((start, len2, t, len1, 0));
                    } else if del > 0 && del % 3 == 0 && ins == 0 && len2 > 0 {
                        annv.push((len1, len2, t, len1 + del, 0));
                    }
                }

                // Keep going.

                if tig_start < 0 || tig_stop < 0 {
                    eprintme!(tig_start, tig_stop);
                    json_error(Some(&json));
                }
                let (tig_start, tig_stop) = (tig_start as usize, tig_stop as usize);
                let quals0 = v["quals"].to_string();
                let quals0 = quals0.after("\"").as_bytes();
                let mut quals = Vec::<u8>::new();
                let mut slashed = false;
                for i in 0..quals0.len() - 1 {
                    if !slashed && quals0[i] == b'\\'
                    /* && ( i == 0 || quals0[i-1] != b'\\' ) */
                    {
                        slashed = true;
                        continue;
                    }
                    slashed = false;
                    quals.push(quals0[i]);
                }
                assert_eq!(full_seq.len(), quals.len());
                let seq = full_seq[tig_start..tig_stop].to_string();
                for i in 0..quals.len() {
                    quals[i] -= 33 as u8;
                }
                let full_quals = quals.clone();
                let quals = quals[tig_start..tig_stop].to_vec();
                // let cdr3_dna = &v["cdr3_seq"].to_string().between("\"", "\"").to_string();
                let umi_count = v["umi_count"].as_i64().unwrap() as usize;
                let read_count = v["read_count"].as_i64().unwrap() as usize;
                tigs.push(TigData {
                    cdr3_dna: cdr3_dna.to_string(),
                    len: seq.len(),
                    seq: seq.as_bytes().to_vec(),
                    v_start: tig_start,
                    v_stop: v_stop,
                    v_stop_ref: v_stop_ref,
                    j_start: j_start,
                    j_start_ref: j_start_ref,
                    j_stop: tig_stop,
                    c_start: c_start,
                    full_seq: full_seq.as_bytes().to_vec(),
                    v_ref_id: v_ref_id,
                    d_ref_id: d_ref_id,
                    j_ref_id: j_ref_id,
                    c_ref_id: c_ref_id,
                    u_ref_id: u_ref_id,
                    cdr3_aa: cdr3_aa.to_string(),
                    cdr3_start: cdr3_start,
                    quals: quals,
                    full_quals: full_quals,
                    barcode: barcode.to_string(),
                    tigname: tigname.to_string(),
                    left: left,
                    dataset_index: li,
                    umi_count: umi_count,
                    read_count: read_count,
                    chain_type: chain_type.clone(),
                    annv: annv.clone(),
                });
            }
        }
    }
    let mut tig_bc = Vec::<Vec<TigData>>::new();
    let mut r = 0;
    while r < tigs.len() {
        let mut s = r + 1;
        while s < tigs.len() {
            if tigs[s].barcode != tigs[r].barcode {
                break;
            }
            s += 1;
        }
        /*
        let (mut have_left, mut have_right) = (false, false);
        for u in r..s {
            if tigs[u].left {
                have_left = true;
            } else {
                have_right = true;
            }
        }
        */

        // For now we require at most four contigs (but we don't yet merge foursies).

        if
        /* have_left && have_right && */
        s - r <= 4 {
            let mut bc_tigs = Vec::<TigData>::new();
            for u in r..s {
                bc_tigs.push(tigs[u].clone());
            }
            bc_tigs.sort();
            tig_bc.push(bc_tigs);
        }
        r = s;
    }
    tig_bc
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the json annotations files.

pub fn parse_json_annotations_files(
    ctl: &EncloneControl,
    tig_bc: &mut Vec<Vec<TigData>>,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
) {
    let tl = Instant::now();
    // (lena index, contig name, V..J length): (?)
    let mut results = Vec::<(
        usize,
        Vec<(String, usize)>,
        Vec<Vec<TigData>>,
        Vec<Vec<u8>>, // logs
    )>::new();
    for i in 0..ctl.sample_info.dataset_path.len() {
        results.push((
            i,
            Vec::<(String, usize)>::new(),
            Vec::<Vec<TigData>>::new(),
            Vec::<Vec<u8>>::new(),
        ));
    }
    // note: only tracking truncated seq and quals initially
    results.par_iter_mut().for_each(|res| {
        let li = res.0;
        let json = format!(
            "{}/all_contig_annotations.json",
            ctl.sample_info.dataset_path[li]
        );
        let json_lz4 = format!(
            "{}/all_contig_annotations.json.lz4",
            ctl.sample_info.dataset_path[li]
        );
        if !path_exists(&json) && !path_exists(&json_lz4) {
            eprintln!("can't find {} or {}", json, json_lz4);
            std::process::exit(1);
        }
        let tig_bc: Vec<Vec<TigData>> =
            read_json(li, &json, &refdata, &to_ref_index, ctl.gen_opt.reannotate);
        explore(li, &tig_bc, &ctl);
        res.2 = tig_bc;
    });
    for i in 0..results.len() {
        tig_bc.append(&mut results[i].2.clone());
    }
    if ctl.comp {
        println!("used {:.1} seconds loading from json", elapsed(&tl));
    }
}
