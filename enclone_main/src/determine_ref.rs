// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Start of code to determine the reference sequence that is to be used.

use enclone_core::defs::*;
use io_utils::*;
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn determine_ref(ctl: &mut EncloneControl, refx: &mut String) {
    // First check for the existence of a json file.

    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
    let mut jsonx = String::new();
    if ctl.origin_info.n() > 0 {
        let json = format!("{}/{}", ctl.origin_info.dataset_path[0], ann);
        let json_lz4 = format!("{}/{}.lz4", ctl.origin_info.dataset_path[0], ann);
        if !path_exists(&json) && !path_exists(&json_lz4) {
            eprintln!(
                "\nUnable to find a VDJ input file: can't find\n{}\nor {}.\n\n\
                There are various possible reasons for this, including an incorrectly \
                specified path, or incorrect\nspecification of PRE, or a partially copied outs \
                directory that does not include \
                all the needed\nfiles, or a mixup between VDJ and GEX path names.\n",
                json, json_lz4
            );
            std::process::exit(1);
        }
        jsonx = json.clone();
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
    }

    // Step 1.  Test to see if CURRENT_REF or BUILT_IN is specified.  Kind of a mess that we
    // have both.

    if ctl.gen_opt.current_ref || ctl.gen_opt.built_in {
        if !ctl.gen_opt.mouse {
            *refx = human_ref();
        } else {
            *refx = mouse_ref();
        }
        if ctl.gen_opt.built_in {
            ctl.gen_opt.reannotate = true;
        }
    }

    // Step 2. Test to see if IMGT is specified.

    if ctl.gen_opt.imgt && ctl.gen_opt.internal_run {
        if !ctl.gen_opt.mouse {
            let imgt =
                "/mnt/opt/refdata_cellranger/vdj/vdj_IMGT_human_20200415-0.0.0/fasta/regions.fa";
            if ctl.gen_opt.descrip {
                println!("using IMGT human reference");
            }
            let f = open_for_read![imgt];
            for line in f.lines() {
                let mut s = line.unwrap();
                if ctl.gen_opt.imgt_fix {
                    // Fix IGHJ6.
                    if s == "ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA"
                        .to_string()
                        || s == "ATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCA"
                            .to_string()
                    {
                        s += "G";
                    }
                }
                *refx += &s;
                *refx += &"\n";
            }
        } else {
            let imgt =
                "/mnt/opt/refdata_cellranger/vdj/vdj_IMGT_mouse_20180723-2.2.0/fasta/regions.fa";
            if ctl.gen_opt.descrip {
                println!("using IMGT mouse reference");
            }
            let f = open_for_read![imgt];
            for line in f.lines() {
                let s = line.unwrap();
                *refx += &s;
                *refx += &"\n";
            }
        }
        ctl.gen_opt.reannotate = true;
    }

    // Step 3.  Test to see if REF is specified.

    if refx.len() == 0 && ctl.gen_opt.refname.len() > 0 {
        if std::path::Path::new(&ctl.gen_opt.refname).is_dir() {
            eprintln!(
                "\nProblem with REF: \"{}\"\nis a directory, not a file.\n",
                ctl.gen_opt.refname
            );
            std::process::exit(1);
        }
        if ctl.gen_opt.descrip {
            println!("using reference = {}", ctl.gen_opt.refname);
        }
        let fx = File::open(&ctl.gen_opt.refname);
        if fx.is_err() {
            eprintln!(
                "\nProblem with REF: unable to read from the file\n\
                 \"{}\".\nPlease check that that path makes sense and that you have read \
                 permission along that path.\n",
                ctl.gen_opt.refname
            );
            std::process::exit(1);
        }
        let f = BufReader::new(fx.unwrap());
        let mut nheader = 0;
        let mut bases = 0;
        let mut na = 0;
        let mut nc = 0;
        let mut ng = 0;
        let mut nt = 0;
        for line in f.lines() {
            let s = line.unwrap();
            *refx += &s;
            *refx += &"\n";
            if s.starts_with('>') {
                nheader += 1;
            } else {
                for c in s.chars() {
                    bases += 1;
                    if c == 'A' || c == 'a' {
                        na += 1;
                    } else if c == 'C' || c == 'c' {
                        nc += 1;
                    } else if c == 'G' || c == 'g' {
                        ng += 1;
                    } else if c == 'T' || c == 't' {
                        nt += 1;
                    }
                }
            }
        }
        if nheader == 0 || bases == 0 || (na + nc + ng + nt) as f64 / (bases as f64) < 0.95 {
            eprintln!("\nProblem with REF: it is not a FASTA file.\n");
            std::process::exit(1);
        }
    }

    // Step 4.  Test for presence of a reference file in the VDJ directories.

    if refx.len() == 0 && ctl.gen_opt.refname.len() == 0 {
        let rpaths = [
            "outs/vdj_reference/fasta/regions.fa",
            "vdj_reference/fasta/regions.fa",
            "regions.fa",
        ];
        let mut refs = Vec::<String>::new();
        for li in 0..ctl.origin_info.n() {
            for r in rpaths.iter() {
                let fasta = format!("{}/{}", ctl.origin_info.dataset_path[li], r);
                if path_exists(&fasta) {
                    refs.push(std::fs::read_to_string(&fasta).unwrap());
                    break;
                }
            }
        }
        if refs.len() > 0 {
            unique_sort(&mut refs);
            if refs.len() > 1 {
                eprintln!(
                    "The VDJ reference sequences that were supplied to Cell Ranger are not \
                    identical with each other.\nAs a consequence, the VDJ output files are not \
                    compatible with each other, so enclone can't run.\nYou have some options as \
                    to how to proceed:\n\
                    1. You can rerun Cell Ranger using the same reference.\n\
                    2. You can select one of the references, and supply that to enclone using the \
                    REF option.\n   You will also need to supply the argument RE to get enclone to \
                    recompute annotations,\n   and that will make it somewhat slower.\n\n"
                );
                std::process::exit(1);
            }
            if ctl.gen_opt.mouse {
                eprintln!(
                    "\nSince the reference sequence is already in the VDJ input directories that\n\
                    you supplied to enclone, it is not necessary to supply the MOUSE argument.\n\
                    Please remove that argument.  Exiting now because of possible unintended\n\
                    consequences.\n"
                );
                std::process::exit(1);
            }
            *refx = refs[0].clone();
        }
    }

    // Step 5.  Attempt to determine the reference that was used by reading far enough into the
    // first json file to find a distinguishing entry.

    if refx.len() == 0 && jsonx.len() > 0 {
        let cr_ver = ["2.0", "3.1", "4.0", "current"];
        let species = ["human", "mouse"];
        let mut refhash = Vec::<(HashMap<usize, (usize, String)>, String)>::new();
        for i in 0..cr_ver.len() {
            for j in 0..species.len() {
                let refx;
                if cr_ver[i] == "2.0" && species[j] == "human" {
                    refx = human_ref_2_0();
                } else if cr_ver[i] == "2.0" && species[j] == "mouse" {
                    continue;
                } else if cr_ver[i] == "3.1" && species[j] == "human" {
                    refx = human_ref_3_1();
                } else if cr_ver[i] == "3.1" && species[j] == "mouse" {
                    refx = mouse_ref_3_1();
                } else if cr_ver[i] == "4.0" && species[j] == "human" {
                    refx = human_ref_4_0();
                } else if cr_ver[i] == "4.0" && species[j] == "mouse" {
                    refx = mouse_ref_4_0();
                } else if cr_ver[i] == "current" && species[j] == "human" {
                    refx = human_ref();
                } else if cr_ver[i] == "current" && species[j] == "mouse" {
                    refx = mouse_ref();
                } else {
                    panic!("Failed match for reference.");
                }
                let mut h = HashMap::<usize, (usize, String)>::new();
                let mut lines = Vec::<String>::new();
                for line in refx.lines() {
                    lines.push(line.to_string());
                }
                for k in 0..lines.len() {
                    if lines[k].starts_with('>') {
                        let id = lines[k].between(">", "|").force_usize();
                        let gene = lines[k].between("|", " ").to_string();
                        let len = lines[k + 1].len();
                        h.insert(id, (len, gene));
                    }
                }
                refhash.push((h, refx));
            }
        }
        let mut to_delete = vec![false; refhash.len()];
        for i1 in 0..refhash.len() {
            for i2 in i1 + 1..refhash.len() {
                if refhash[i1].0 == refhash[i2].0 {
                    to_delete[i2] = true;
                }
            }
        }
        erase_if(&mut refhash, &to_delete);
        let mut f = BufReader::new(open_maybe_compressed(&jsonx));
        'json_entry: loop {
            let x = read_vector_entry_from_json(&mut f);
            if x.is_none() {
                break;
            }
            let v: Value = serde_json::from_str(strme(&x.unwrap())).unwrap();
            let ann = v["annotations"].as_array().unwrap();
            for i in 0..ann.len() {
                let a = &ann[i];
                let id = a["feature"]["feature_id"].as_u64().unwrap() as usize;
                let gene = a["feature"]["gene_name"].to_string();
                let gene = gene.between("\"", "\"").to_string();
                let len = a["annotation_length"].as_u64().unwrap() as usize;
                let mut matches = Vec::<usize>::new();
                for j in 0..refhash.len() {
                    if refhash[j].0.contains_key(&id) && refhash[j].0[&id] == (len, gene.clone()) {
                        matches.push(j);
                    }
                }
                if matches.is_empty() {
                    break 'json_entry;
                } else if matches.solo() {
                    *refx = refhash[matches[0]].1.clone();
                    break 'json_entry;
                }
            }
        }
    }
    if refx.len() == 0 && jsonx.len() > 0 {
        eprintln!(
            "\nenclone was unable to determine the reference sequence that you used.  You \
            have two options:\n\
            1. If you used cellranger version 4.0 or later, copy the vdj_reference directory\n   \
               from there to the outs directory that contains your other enclone input data.\n\
            2. Use the REF argument to specify the name of the reference fasta file.\n",
        );
        std::process::exit(1);
    }
}
