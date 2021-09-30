// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use self::annotate::{annotate_seq, get_cdr3_using_ann, print_some_annotations};
use self::refx::RefData;
use self::transcript::is_valid;
use debruijn::dna_string::DnaString;
use enclone_core::defs::{EncloneControl, OriginInfo, TigData};
use io_utils::{open_maybe_compressed, path_exists, read_vector_entry_from_json};
use rayon::prelude::*;
use serde_json::Value;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;
use std::{collections::HashMap, io::BufReader};
use string_utils::{stringme, strme, TextUtils};
use vdj_ann::{annotate, refx, transcript};
use vector_utils::{bin_position, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn json_error(
    json: Option<&str>,
    ctl: &EncloneControl,
    exiting: &AtomicBool,
    msg: &str,
) -> Result<(), String> {
    // The following line prevents error messages from this function from being
    // printed multiple times.
    let mut msgx = String::new();
    if !exiting.swap(true, Ordering::Relaxed) {
        msgx = "\nThere is something wrong with the contig annotations in the cellranger output \
             file"
            .to_string();
        if json.is_some() {
            msgx += &mut format!("\n{}.", json.unwrap());
        } else {
            msgx += ".";
        }
        if ctl.gen_opt.internal_run {
            msgx += &mut format!("\n\npossibly relevant internal data: {}\n", msg);
        }
        if ctl.gen_opt.internal_run {
            msgx += &mut "\n\nATTENTION INTERNAL 10X USERS!\n\
                Quite possibly you are using data from a cellranger run carried out using a \
                version\n\
                between 3.1 and 4.0.  For certain of these versions, it is necessary to add the\n\
                argument CURRENT_REF to your command line.  If that doesn't work, \
                please see below.\n"
                .to_string();
        }
        msgx += &mut "\n\nHere is what you should do:\n\n\
             1. If you used cellranger version ≥ 4.0, the problem is very likely\n\
                that the directory outs/vdj_reference was not retained, so enclone\n\
                didn't see it, and had to guess what the reference sequence was.\n\
                Fix this and everything should be fine.\n\n\
             2. If you used cellranger version 3.1, then you need to add a command-line\n\
                argument REF=<vdj_reference_fasta_file_name>, or if you already did that,\n\
                make sure it is the *same* as that which you gave cellranger.\n\n\
             3. If you used cellranger version < 3.1 (the only other possibility), then\n\
                you have options:\n\
                • rerun cellranger using the current version\n\
                • or provide an argument REF= as above and RE to force reannotation\n\
                • or provide the argument BUILT_IN to use the current reference and force\n  \
                  reannotation (and MOUSE if you used mouse); only works with human and mouse.\n\n\
             Note that one way to get the error is to specify TCR when you meant BCR, or the\n\
             other way.\n\n\
             If you're stuck, please write to us at enclone@10xgenomics.com.\n"
            .to_string();
    }
    Err(msgx)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn parse_vector_entry_from_json(
    x: &Vec<u8>,
    json: &String,
    accept_inconsistent: bool,
    origin_info: &OriginInfo,
    li: usize,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    reannotate: bool,
    ctl: &EncloneControl,
    vdj_cells: &mut Vec<String>,
    gex_cells: &mut Vec<String>,
    gex_cells_specified: &mut bool,
    cr_version: &mut String,
    tigs: &mut Vec<TigData>,
    exiting: &AtomicBool,
) -> Result<(), String> {
    let v: Result<Value, _> = serde_json::from_str(strme(x));
    if v.is_err() {
        return Err(format!(
            "\nInternal error, failed to parse a value from a string.  The string is:\n{}\n",
            strme(x)
        ));
    }
    let v = v.unwrap();
    let barcode = &v["barcode"].to_string().between("\"", "\"").to_string();

    // Get cell status.  Sometime after CR 4.0 was released, and before 4.1 was released,
    // we added new fields is_asm_cell and is_gex_cell to the json file.  The value of
    // is_asm_cell is the original determination of "cell" in the VDJ pipeline, whereas the
    // value of is_gex_cell is that for the GEX pipeline.

    let mut is_cell = v["is_cell"].as_bool().unwrap_or(false);
    let is_asm_cell = v["is_asm_cell"].as_bool().unwrap_or(false);
    if is_asm_cell {
        is_cell = true;
    }

    let is_gex_cell = v["is_gex_cell"].as_bool();
    if is_gex_cell.is_some() {
        *gex_cells_specified = true;
    }
    if is_gex_cell == Some(true) {
        gex_cells.push(barcode.clone());
    }

    if !ctl.gen_opt.ncell && !is_cell {
        return Ok(());
    }
    if is_cell {
        vdj_cells.push(barcode.clone());
    }

    // Proceed.

    if !ctl.gen_opt.reprod && !v["productive"].as_bool().unwrap_or(false) {
        return Ok(());
    }
    if !ctl.gen_opt.reprod && !ctl.gen_opt.ncell && !v["high_confidence"].as_bool().unwrap_or(false)
    {
        return Ok(());
    }
    let tigname = &v["contig_name"].to_string().between("\"", "\"").to_string();
    let full_seq = &v["sequence"].to_string().between("\"", "\"").to_string();
    let mut left = false;
    let (mut v_ref_id, mut j_ref_id) = (1000000, 0);
    let mut d_ref_id: Option<usize> = None;
    let mut c_ref_id = None;
    let mut chain_type = String::new();
    let mut u_ref_id = None;
    let (mut tig_start, mut tig_stop) = (-1_isize, -1_isize);
    let mut v_stop = 0;
    let mut v_stop_ref = 0;
    let mut d_start = None;
    let mut j_start = 0;
    let mut j_start_ref = 0;
    let mut c_start = None;
    let mut annv = Vec::<(i32, i32, i32, i32, i32)>::new();
    let mut cdr3_aa: String;
    let mut cdr3_dna: String;
    let mut cdr3_start: usize;
    if v.get("version").is_some() {
        *cr_version = v["version"].to_string().between("\"", "\"").to_string();
    }

    // Read validated and non-validated UMIs.

    let mut validated_umis = Vec::<String>::new();
    let mut validated_umis_present = false;
    let val = v["validated_umis"].as_array();
    if val.is_some() {
        validated_umis_present = true;
        let val = val.unwrap();
        for i in 0..val.len() {
            validated_umis.push(val[i].to_string().between("\"", "\"").to_string());
        }
    }
    let mut non_validated_umis = Vec::<String>::new();
    let mut non_validated_umis_present = false;
    let non_val = v["non_validated_umis"].as_array();
    if non_val.is_some() {
        non_validated_umis_present = true;
        let non_val = non_val.unwrap();
        for i in 0..non_val.len() {
            non_validated_umis.push(non_val[i].to_string().between("\"", "\"").to_string());
        }
    }
    let mut invalidated_umis = Vec::<String>::new();
    let mut invalidated_umis_present = false;
    let inval = v["invalidated_umis"].as_array();
    if inval.is_some() {
        invalidated_umis_present = true;
        let inval = inval.unwrap();
        for i in 0..inval.len() {
            invalidated_umis.push(inval[i].to_string().between("\"", "\"").to_string());
        }
    }

    // Read fraction_of_reads_for_this_barcode_provided_as_input_to_assembly.

    let mut frac_reads_used = None;
    let f = v["fraction_of_reads_for_this_barcode_provided_as_input_to_assembly"].as_f64();
    if f.is_some() {
        frac_reads_used = Some((f.unwrap() * 1_000_000.0).round() as u32);
    }

    // Reannotate.

    if reannotate || ctl.gen_opt.reprod {
        let x = DnaString::from_dna_string(full_seq);
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&x, refdata, &mut ann, true, false, true);

        // If there are multiple V segment alignments, possibly reduce to just one.

        let mut ann2 = Vec::<(i32, i32, i32, i32, i32)>::new();
        let mut j = 0;
        while j < ann.len() {
            let t = ann[j].2 as usize;
            let mut k = j + 1;
            while k < ann.len() {
                if refdata.segtype[ann[k].2 as usize] != refdata.segtype[t] {
                    break;
                }
                k += 1;
            }
            if refdata.segtype[t] == *"V" && k - j > 1 {
                let mut entries = 1;
                if j < ann.len() - 1
                    && ann[j + 1].2 as usize == t
                    && ((ann[j].0 + ann[j].1 == ann[j + 1].0 && ann[j].3 + ann[j].1 < ann[j + 1].3)
                        || (ann[j].0 + ann[j].1 < ann[j + 1].0
                            && ann[j].3 + ann[j].1 == ann[j + 1].3))
                {
                    entries = 2;
                }
                for l in j..j + entries {
                    ann2.push(ann[l]);
                }
            } else {
                for l in j..k {
                    ann2.push(ann[l]);
                }
            }
            j = k;
        }
        ann = ann2;

        // Proceed.

        if ctl.gen_opt.trace_barcode == *barcode {
            let mut log = Vec::<u8>::new();
            print_some_annotations(refdata, &ann, &mut log, false);
            print!("\n{}", strme(&log));
        }
        let mut log = Vec::<u8>::new();
        if ctl.gen_opt.trace_barcode == *barcode {
            if !is_valid(&x, refdata, &ann, true, &mut log) {
                print!("{}", strme(&log));
                println!("invalid");
                return Ok(());
            }
        } else if !is_valid(&x, refdata, &ann, false, &mut log) {
            return Ok(());
        }
        let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        get_cdr3_using_ann(&x, refdata, &ann, &mut cdr3);
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
                annv.push(ann[i]);
                chain_type = refdata.name[t][0..3].to_string();
                if chain_type == *"IGH" || chain_type == *"TRB" {
                    left = true;
                }
                if ann[i].3 == 0 {
                    tig_start = ann[i].0 as isize;
                    if tig_start > cdr3_start as isize {
                        panic!(
                            "Something is wrong with the CDR3 start for this contig:\n\n{}.",
                            &full_seq
                        );
                    }
                    cdr3_start -= tig_start as usize;
                }
                v_stop = (ann[i].0 + ann[i].1) as usize;
                v_stop_ref = (ann[i].3 + ann[i].1) as usize;
            } else if refdata.is_d(t) {
                d_start = Some(ann[i].0 as usize);
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
            let feature_idx = to_ref_index[&feature_id];
            let ref_start = a["annotation_match_start"].as_u64().unwrap() as usize;
            if region_type == "L-REGION+V-REGION" {
                v_stop = a["contig_match_end"].as_i64().unwrap() as usize;
                v_stop_ref = a["annotation_match_end"].as_i64().unwrap() as usize;
            }
            let gene_name = a["feature"]["gene_name"]
                .to_string()
                .between("\"", "\"")
                .to_string();
            if refdata.name[feature_idx] != gene_name
                && !accept_inconsistent
                && !exiting.swap(true, Ordering::Relaxed)
            {
                return Err(format!(
                    "\nThere is an inconsistency between the reference \
                     file used to create the Cell Ranger output files in\n{}\nand the \
                     reference that enclone is using.\n\nFor example, the feature \
                     numbered {} is\nthe gene {} in one and the gene {} in the other.\n\n\
                     As far as we know, this type of error can only occur with Cell Ranger \
                     versions before 4.0.\n\n\
                     If this is mouse data, please use the argument MOUSE, and that may \
                     solve the problem.\n\n\
                     If this is human or mouse data, and you are OK with using the current \
                     built-in reference that\nenclone has, \
                     you can instead add the argument BUILT_IN to the command line.  This \
                     forces\nrecomputation of annotations and may be somewhat slower.\n\n\
                     A solution that should always work is to supply\n\
                     REF=vdj_reference_fasta_filename as an argument to enclone.\n",
                    json.rev_before("/"),
                    feature_id,
                    gene_name,
                    refdata.name[feature_idx]
                ));
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
                if chain == *"IGH" || chain == *"TRB" {
                    left = true;
                }
                v_ref_id = feature_idx;
                cigarv = a["cigar"].to_string().between("\"", "\"").to_string();
            } else {
                // also check for IG chain?????????????????????????????????????????
                let ref_stop = a["annotation_match_end"].as_u64().unwrap() as usize;
                let ref_len = a["annotation_length"].as_u64().unwrap() as usize;
                if region_type == "J-REGION" && ref_stop == ref_len {
                    tig_stop = a["contig_match_end"].as_i64().unwrap() as isize;
                    j_ref_id = feature_idx;
                    j_start = a["contig_match_start"].as_i64().unwrap() as usize;
                    j_start_ref = a["annotation_match_start"].as_i64().unwrap() as usize;
                }
                if region_type == "5'UTR" {
                    u_ref_id = Some(feature_idx);
                }
                if region_type == "D-REGION" {
                    d_start = Some(a["contig_match_start"].as_i64().unwrap() as usize);
                    d_ref_id = Some(feature_idx);
                }
                if region_type == "C-REGION" {
                    c_ref_id = Some(feature_idx);
                    c_start = Some(a["contig_match_start"].as_i64().unwrap() as usize);
                }
            }
        }
        if v_ref_id == 1000000 {
            return Ok(());
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
        annv.push((0_i32, len1, t, 0, 0));
        if ins > 0 && ins % 3 == 0 && del == 0 && len2 > 0 {
            let start = (len1 + ins) as i32;
            annv.push((start, len2, t, len1, 0));
        } else if del > 0 && del % 3 == 0 && ins == 0 && len2 > 0 {
            annv.push((len1, len2, t, len1 + del, 0));
        }
        let rt = &refdata.refs[v_ref_id as usize];
        if annv.len() == 2 && annv[0].1 as usize > rt.len() {
            let msg = format!("annv[0].1 = {}, rt.len() = {}", annv[0].1, rt.len());
            json_error(None, ctl, exiting, &msg)?;
        }

        // Check to see if the CDR3 sequence has changed.  This could happen if the cellranger
        // version for all_contig_annotations.json used an older version of the CDR3 calculation
        // than is used in the current version of enclone.  This could result in internal
        // inconsistencies, leading to an assert somewhere downstream.

        let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        let x = DnaString::from_dna_string(full_seq);
        get_cdr3_using_ann(&x, refdata, &annv, &mut cdr3);
        if cdr3.is_empty() {
            return Ok(());
        }
        let cdr3_aa_alt = stringme(&cdr3[0].1);
        if cdr3_aa != cdr3_aa_alt {
            // This is particularly pathological and rare:

            if tig_start as usize > cdr3[0].0 {
                return Ok(());
            }

            // Define start.

            cdr3_start = cdr3[0].0 - tig_start as usize;

            // Define cdr3.

            cdr3_aa = cdr3_aa_alt;
            cdr3_dna = x
                .slice(cdr3_start, cdr3_start + 3 * cdr3_aa.len())
                .to_string();
        }
    }

    // Test for two very rare conditions where the CDR3 is busted.  This could be confusing to
    // users if they hit one of these.
    // Case 1: seen on 47680, barcode CGCCAAGTCCATGAAC-1.
    // Case 2: seen on 99640, barcode CAGTAACCATGTCGAT-1.
    // It is not known if these correspond to bugs in cellranger that were subsequently fixed.

    if cdr3_aa.contains('*') {
        return Ok(());
    }
    if cdr3_start + 3 * cdr3_aa.len() > tig_stop as usize - tig_start as usize {
        return Ok(());
    }

    // Keep going.

    if tig_start < 0 || tig_stop < 0 {
        let msg = format!("tig_start = {}, tig_stop = {}", tig_start, tig_stop);
        json_error(Some(json), ctl, exiting, &msg)?;
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
        quals[i] -= 33_u8;
    }
    let full_quals = quals.clone();
    let quals = quals[tig_start..tig_stop].to_vec();
    let umi_count = v["umi_count"].as_i64().unwrap() as usize;
    let read_count = v["read_count"].as_i64().unwrap() as usize;
    let mut origin = None;
    let mut donor = None;
    let mut tag = None;
    if origin_info.origin_for_bc[li].contains_key(&barcode.clone()) {
        origin = Some(origin_info.origin_for_bc[li][&barcode.clone()].clone());
    } else {
        // the way we use s1 here is flaky
        if !origin_info.origin_id[li].is_empty()
            && (origin_info.origin_id[li] != *"s1" || origin_info.origin_for_bc[li].is_empty())
        {
            origin = Some(origin_info.origin_id[li].clone());
        }
    }
    if origin_info.donor_for_bc[li].contains_key(&barcode.clone()) {
        donor = Some(origin_info.donor_for_bc[li][&barcode.clone()].clone());
    } else {
        // the way we use d1 here is flaky
        if !origin_info.origin_id[li].is_empty()
            && (origin_info.donor_id[li] != *"d1" || origin_info.donor_for_bc[li].is_empty())
        {
            donor = Some(origin_info.donor_id[li].clone());
        }
    }
    if origin_info.tag[li].contains_key(&barcode.clone()) {
        tag = Some(origin_info.tag[li][&barcode.clone()].clone());
    }
    let mut origin_index = None;
    let mut donor_index = None;
    let mut tag_index = None;
    if origin.is_some() {
        if origin.is_some() {
            origin_index = Some(bin_position(&origin_info.origin_list, &origin.unwrap()) as usize);
        }
        if donor.is_some() {
            donor_index = Some(bin_position(&origin_info.donor_list, &donor.unwrap()) as usize);
        }
    }
    if tag.is_some() {
        tag_index = Some(bin_position(&origin_info.tag_list, &tag.unwrap()) as usize);
    }
    let mut valu = None;
    if validated_umis_present {
        valu = Some(validated_umis);
    }
    let mut non_valu = None;
    if non_validated_umis_present {
        non_valu = Some(non_validated_umis);
    }
    let mut invalu = None;
    if invalidated_umis_present {
        invalu = Some(invalidated_umis);
    }
    tigs.push(TigData {
        cdr3_dna,
        len: seq.len(),
        v_start: tig_start,
        v_stop,
        v_stop_ref,
        d_start,
        j_start,
        j_start_ref,
        j_stop: tig_stop,
        c_start,
        full_seq: full_seq.as_bytes().to_vec(),
        v_ref_id,
        d_ref_id,
        j_ref_id,
        c_ref_id,
        u_ref_id,
        fr1_start: 0,
        cdr1_start: None,
        fr2_start: None,
        cdr2_start: None,
        fr3_start: None,
        cdr3_aa,
        cdr3_start,
        quals,
        full_quals,
        barcode: barcode.to_string(),
        tigname: tigname.to_string(),
        left,
        dataset_index: li,
        origin_index,
        donor_index,
        tag_index,
        umi_count,
        read_count,
        chain_type,
        annv: annv.clone(),
        validated_umis: valu,
        non_validated_umis: non_valu,
        invalidated_umis: invalu,
        frac_reads_used,
    });
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the JSON annotations file.
//
// In the future could be converted to LazyWrite:
// https://martian-lang.github.io/martian-rust/doc/martian_filetypes/json_file/
// index.html#lazy-readwrite-example.
//
// Tracking contigs using bc_cdr3_aa; could improve later.
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
// this code.  Another would be to write out a binary version of the JSON file that contains
// only the information that we need.

pub fn read_json(
    accept_inconsistent: bool,
    origin_info: &OriginInfo,
    li: usize,
    json: &String,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    reannotate: bool,
    cr_version: &mut String,
    ctl: &EncloneControl,
    mut vdj_cells: &mut Vec<String>,
    mut gex_cells: &mut Vec<String>,
    gex_cells_specified: &mut bool,
) -> Result<Vec<Vec<TigData>>, String> {
    *gex_cells_specified = false;
    let mut tigs = Vec::<TigData>::new();
    let mut jsonx = json.clone();
    if !path_exists(json) {
        jsonx = format!("{}.lz4", json);
    }
    if jsonx.contains('/') {
        let p = jsonx.rev_before("/");
        if !path_exists(p) {
            return Err(format!(
                "\nThere should be a directory\n\
                 \"{}\"\n\
                 but it does not exist.  Please check how you have specified the\n\
                 input files to enclone, including the PRE argument.\n",
                p
            ));
        }
    }
    if !path_exists(&jsonx) {
        return Err(format!(
            "\nThe path\n\
             \"{}\"\n\
             does not exist.  Please check how you have specified the\n\
             input files to enclone, including the PRE argument.\n",
            jsonx
        ));
    }
    let mut f = BufReader::new(open_maybe_compressed(&jsonx));
    // ◼ This loop could be speeded up, see comments above.
    let mut xs = Vec::<Vec<u8>>::new();
    loop {
        match read_vector_entry_from_json(&mut f) {
            None => break,
            Some(x) => {
                xs.push(x);
            }
        }
    }
    let mut results = Vec::<(
        usize,
        Vec<String>,
        Vec<String>,
        bool,
        String,
        Vec<TigData>,
        String,
    )>::new();
    for i in 0..xs.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<String>::new(),
            false,
            String::new(),
            Vec::<TigData>::new(),
            String::new(),
        ));
    }
    let exiting = AtomicBool::new(false);
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let resx = parse_vector_entry_from_json(
            &xs[i],
            json,
            accept_inconsistent,
            origin_info,
            li,
            refdata,
            to_ref_index,
            reannotate,
            ctl,
            &mut res.1,
            &mut res.2,
            &mut res.3,
            &mut res.4,
            &mut res.5,
            &exiting,
        );
        if resx.is_err() {
            res.6 = resx.unwrap_err();
        }
    });
    for i in 0..results.len() {
        if !results[i].6.is_empty() {
            return Err(results[i].6.clone());
        }
    }
    for i in 0..xs.len() {
        vdj_cells.append(&mut results[i].1);
        gex_cells.append(&mut results[i].2);
        if results[i].3 {
            *gex_cells_specified = true;
        }
        if !results[i].4.is_empty() {
            *cr_version = results[i].4.clone();
        }
        tigs.append(&mut results[i].5);
    }
    unique_sort(&mut gex_cells);
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

        // For now we require at most four contigs (but we don't yet merge foursies).

        if s - r <= 4 {
            let mut bc_tigs = Vec::<TigData>::new();
            for u in r..s {
                bc_tigs.push(tigs[u].clone());
            }
            bc_tigs.sort();
            tig_bc.push(bc_tigs);
        }
        r = s;
    }
    unique_sort(&mut vdj_cells);
    Ok(tig_bc)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the JSON annotations file(s).

pub fn parse_json_annotations_files(
    ctl: &EncloneControl,
    tig_bc: &mut Vec<Vec<TigData>>,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    vdj_cells: &mut Vec<Vec<String>>,
    gex_cells: &mut Vec<Vec<String>>,
    gex_cells_specified: &mut Vec<bool>,
) -> Result<(), String> {
    // (origin index, contig name, V..J length): (?)
    let mut results = Vec::<(
        usize,
        Vec<(String, usize)>,
        Vec<Vec<TigData>>,
        Vec<Vec<u8>>, // logs
        String,
        Vec<String>,
        Vec<String>,
        bool,
        String,
    )>::new();
    for i in 0..ctl.origin_info.dataset_path.len() {
        results.push((
            i,
            Vec::<(String, usize)>::new(),
            Vec::<Vec<TigData>>::new(),
            Vec::<Vec<u8>>::new(),
            String::new(),
            Vec::<String>::new(),
            Vec::<String>::new(),
            false,
            String::new(),
        ));
    }
    // Note: only tracking truncated seq and quals initially
    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
    results.par_iter_mut().for_each(|res| {
        let li = res.0;
        let json = format!("{}/{}", ctl.origin_info.dataset_path[li], ann);
        let json_lz4 = format!("{}/{}.lz4", ctl.origin_info.dataset_path[li], ann);
        if !path_exists(&json) && !path_exists(&json_lz4) {
            res.8 = format!("\ncan't find {} or {}\n", json, json_lz4);
            return;
        }
        let resx = read_json(
            ctl.gen_opt.accept_inconsistent,
            &ctl.origin_info,
            li,
            &json,
            refdata,
            to_ref_index,
            ctl.gen_opt.reannotate,
            &mut res.4,
            ctl,
            &mut res.5,
            &mut res.6,
            &mut res.7,
        );
        if resx.is_ok() {
            let tig_bc: Vec<Vec<TigData>> = resx.unwrap();
            res.5.sort();
            res.2 = tig_bc;
        } else {
            res.8 = resx.err().unwrap();
        }
    });
    for i in 0..results.len() {
        if !results[i].8.is_empty() {
            return Err(results[i].8.clone());
        }
    }
    let mut versions = Vec::<String>::new();
    for i in 0..results.len() {
        tig_bc.append(&mut results[i].2.clone());
        // ctl.gen_opt.cr_version = results[i].4.clone();
        if results[i].4.is_empty() {
            versions.push("≤3.1".to_string());
        } else {
            versions.push(results[i].4.clone());
        }
        vdj_cells.push(results[i].5.clone());
        gex_cells.push(results[i].6.clone());
        gex_cells_specified.push(results[i].7);
    }
    /*
    if !ctl.gen_opt.internal_run {
        unique_sort(&mut versions);
        if versions.len() > 1
            && versions != vec!["4.0".to_string(), "4009.52.0-82-g2244c685a".to_string()]
        {
            let args: Vec<String> = env::args().collect();
            return Err(format!(
                "\nYou're using output from multiple Cell Ranger versions = {},\n\
                 which is not allowed.  Your command was:\n{}\n",
                versions.iter().format(", "),
                args.iter().format(","),
            ));
        }
    }
    */
    Ok(())
}
