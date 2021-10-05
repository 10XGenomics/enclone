// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Start of code to determine the reference sequence that is to be used.

use enclone_core::defs::EncloneControl;
use io_utils::{open_maybe_compressed, path_exists, read_vector_entry_from_json};
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::{strme, TextUtils};
use vdj_ann::refx::{
    human_ref, human_ref_2_0, human_ref_3_1, human_ref_4_0, mouse_ref, mouse_ref_3_1, mouse_ref_4_0,
};
use vector_utils::{erase_if, unique_sort, VecUtils};

pub fn determine_ref(ctl: &mut EncloneControl, refx: &mut String) -> Result<(), String> {
    // First check for the existence of a json file.

    let ann = "contig_annotations.json";
    let mut jsonx = String::new();
    if ctl.origin_info.n() > 0 {
        let json = format!("{}/{}", ctl.origin_info.dataset_path[0], ann);
        let json_lz4 = format!("{}/{}.lz4", ctl.origin_info.dataset_path[0], ann);
        if !path_exists(&json) && !path_exists(&json_lz4) {
            return Err(format!(
                "\nUnable to find a VDJ input file: can't find\n{}\nor {}.\n\n\
                There are various possible reasons for this, including an incorrectly \
                specified path, or incorrect\nspecification of PRE, or a partially copied outs \
                directory that does not include \
                all the needed\nfiles, or a mixup between VDJ and GEX path names.\n",
                json, json_lz4
            ));
        }
        jsonx = json.clone();
        if !path_exists(&json) {
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
    }

    // Step 1.  Test to see if REF is specified.

    if refx.is_empty() && !ctl.gen_opt.refname.is_empty() {
        if std::path::Path::new(&ctl.gen_opt.refname).is_dir() {
            return Err(format!(
                "\nProblem with REF: \"{}\"\nis a directory, not a file.\n",
                ctl.gen_opt.refname
            ));
        }
        if ctl.gen_opt.descrip {
            println!("using reference = {}", ctl.gen_opt.refname);
        }
        let fx = File::open(&ctl.gen_opt.refname);
        if fx.is_err() {
            return Err(format!(
                "\nProblem with REF: unable to read from the file\n\
                 \"{}\".\nPlease check that that path makes sense and that you have read \
                 permission along that path.\n",
                ctl.gen_opt.refname
            ));
        }
        let f = BufReader::new(fx.unwrap());
        for line in f.lines() {
            let s = line.unwrap();
            *refx += &s;
            *refx += "\n";
        }
    }

    // Step 2.  Test for presence of a reference file in the VDJ directories.

    if refx.is_empty() && ctl.gen_opt.refname.is_empty() {
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
        if !refs.is_empty() {
            unique_sort(&mut refs);
            *refx = refs[0].clone();
        }
    }

    // Step 3.  Attempt to determine the reference that was used by reading far enough into the
    // first json file to find a distinguishing entry.

    if refx.is_empty() && !jsonx.is_empty() {
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
    Ok(())
}
