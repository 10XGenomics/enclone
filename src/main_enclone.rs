// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use vdj_ann::*;

use self::refx::*;
use crate::allele::*;
use crate::defs::*;
use crate::explore::*;
use crate::graph_filter::*;
use crate::info::*;
use crate::join::*;
use crate::load_gex::*;
use crate::loupe::*;
use crate::misc1::*;
use crate::misc2::*;
use crate::misc3::*;
use crate::print_clonotypes::*;
use crate::proc_args2::*;
use crate::proc_args_check::*;
use crate::read_json::*;
use debruijn::dna_string::*;
use equiv::EquivRel;
use io_utils::*;
use perf_stats::*;
use regex::Regex;
use serde_json::Value;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

pub fn main_enclone(args: &Vec<String>) {
    // Set up stuff, read args, etc.

    let tall = Instant::now();
    let mut ctl = EncloneControl::default();
    setup(&mut ctl, &args);

    // Read external data.

    if ctl.gen_opt.ext.len() > 0 {
        let f = open_for_read![&ctl.gen_opt.ext];
        let mut exts = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            ctl.gen_opt.extc.insert(
                (fields[0].to_string(), fields[1].to_string()),
                fields[2].to_string(),
            );
            exts.push(fields[2].to_string());
        }
        ctl.clono_print_opt.lvars.push("ext".to_string());
        exts.sort();
        let mut i = 0;
        while i < exts.len() {
            let j = next_diff(&exts, i);
            ctl.gen_opt.extn.insert(exts[i].clone(), j - i);
            i = j;
        }
    }

    // Get gene expression and antibody counts.  Sanity check variables in cases where that
    // has to occur after loading gex data.  Actually, it could occur after loading only
    // the feature list, which would be better.

    let gex_info = get_gex_info(&mut ctl);
    check_lvars(&ctl, &gex_info);
    check_pcols(&ctl, &gex_info);

    // Find matching features for <regular expression>_g etc.

    ctl.clono_print_opt.lvars_match =
        vec![vec![Vec::<usize>::new(); ctl.clono_print_opt.lvars.len()]; ctl.sample_info.n()];
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for x in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
        }
    }
    for (i, x) in ctl.clono_print_opt.lvars.iter().enumerate() {
        for y in ends.iter() {
            if x.ends_with(y) {
                let p = x.rev_before(y);
                if !p.is_empty() && Regex::new(&p).is_ok() {
                    let mut ok = true;
                    let mut px = false;
                    let b = p.as_bytes();
                    for i in 0..p.len() {
                        if !((b[i] >= b'A' && b[i] <= b'Z')
                            || (b[i] >= b'a' && b[i] <= b'z')
                            || (b[i] >= b'0' && b[i] <= b'9')
                            || b".-_[]()*".contains(&b[i]))
                        {
                            ok = false;
                            break;
                        }
                        if b"[]()*".contains(&b[i]) {
                            px = true;
                        }
                    }
                    if ok && px {
                        let reg = Regex::new(&format!("^{}$", p));
                        for li in 0..ctl.sample_info.n() {
                            for j in 0..gex_info.gex_features[li].len() {
                                let f = &gex_info.gex_features[li][j];
                                let ff = f.split('\t').collect::<Vec<&str>>();
                                let mut ok = false;
                                if ff[2].starts_with("Antibody") {
                                    if y.contains("_ab") {
                                        ok = true;
                                    }
                                } else if ff[2].starts_with("Antigen") {
                                    if y.contains("_ag") {
                                        ok = true;
                                    }
                                } else if ff[2].starts_with("CRISPR") {
                                    if y.contains("_cr") {
                                        ok = true;
                                    }
                                } else if ff[2].starts_with("Custom") {
                                    if y.contains("_cu") {
                                        ok = true;
                                    }
                                } else if y.contains("_g") {
                                    ok = true;
                                }
                                if ok
                                    && (reg.as_ref().unwrap().is_match(&ff[0])
                                        || reg.as_ref().unwrap().is_match(&ff[1]))
                                {
                                    ctl.clono_print_opt.lvars_match[li][i].push(j);
                                }
                            }
                        }
                        let mut matches = false;
                        for li in 0..ctl.sample_info.n() {
                            if !ctl.clono_print_opt.lvars_match[i][li].is_empty() {
                                matches = true;
                            }
                        }
                        if !matches {
                            eprintln!(
                                "\nLead variable {} contains a pattern that matches \
                                no features.\n",
                                x
                            );
                            std::process::exit(1);
                        }
                        break;
                    }
                }
            }
        }
    }

    // Determine the Cell Ranger version that was used.  Really painful.

    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
    let json = format!("{}/{}", ctl.sample_info.dataset_path[0], ann);
    let json_lz4 = format!("{}/{}.lz4", ctl.sample_info.dataset_path[0], ann);
    if !path_exists(&json) && !path_exists(&json_lz4) {
        eprintln!("can't find {} or {}", json, json_lz4);
        std::process::exit(1);
    }
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
    match read_vector_entry_from_json(&mut f) {
        None => {
            eprintln!("\nFailure reading {}.\n", jsonx);
        }
        Some(x) => {
            let v: Value = serde_json::from_str(strme(&x)).unwrap();
            if v.get("version").is_some() {
                ctl.gen_opt.cr_version = v["version"].to_string().between("\"", "\"").to_string();
            }
        }
    }
    if ctl.gen_opt.current_ref || ctl.gen_opt.cellranger {
        ctl.gen_opt.cr_version = "4.0".to_string();
    }

    // Build reference data.

    let tr = Instant::now();
    let mut refdata = RefData::new();
    let mut refx = String::new();
    if ctl.gen_opt.refname.len() > 0 {
        if std::path::Path::new(&ctl.gen_opt.refname).is_dir() {
            eprintln!(
                "\nProblem with REF: \"{}\"\nis a directory, not a file.\n",
                ctl.gen_opt.refname
            );
            std::process::exit(1);
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
            refx += &s;
            refx += &"\n";
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
    } else if ctl.gen_opt.mouse {
        if ctl.gen_opt.cr_version == "".to_string() && !ctl.gen_opt.reannotate {
            refx = mouse_ref_old();
        } else {
            refx = mouse_ref();
        }
    } else {
        if ctl.gen_opt.cr_version == "".to_string() && !ctl.gen_opt.reannotate {
            refx = human_ref_old();
        } else {
            refx = human_ref();
        }
    }
    let ext_refx = String::new();
    let (mut is_tcr, mut is_bcr) = (true, true);
    if ctl.gen_opt.tcr {
        is_bcr = false;
    }
    if ctl.gen_opt.bcr {
        is_tcr = false;
    }
    make_vdj_ref_data_core(&mut refdata, &refx, &ext_refx, is_tcr, is_bcr, None);
    let mut to_ref_index = HashMap::<usize, usize>::new();
    for i in 0..refdata.refs.len() {
        to_ref_index.insert(refdata.id[i] as usize, i);
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds building reference, peak mem = {:.2} GB",
            elapsed(&tr),
            peak_mem_usage_gb()
        );
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds up through building reference",
            elapsed(&tall)
        );
    }

    // Parse the json annotations file.

    let mut tig_bc = Vec::<Vec<TigData>>::new();
    parse_json_annotations_files(&mut ctl, &mut tig_bc, &refdata, &to_ref_index);

    // Search for SHM indels.  Exploratory.

    let tproto = Instant::now();
    search_for_shm_indels(&ctl, &tig_bc);

    // Filter using light --> heavy graph.

    if !ctl.gen_opt.ngraph_filter {
        graph_filter(&mut tig_bc, ctl.gen_opt.graph);
    }

    // Sort tig_bc.

    sort_tig_bc(&ctl, &mut tig_bc, &refdata);

    // Cross filter.

    cross_filter(&ctl, &mut tig_bc);

    // Look for barcode reuse.

    check_for_barcode_reuse(&ctl, &tig_bc);
    if ctl.comp {
        println!("used {:.2} seconds in proto stuff", elapsed(&tproto));
    }

    // Find exact subclonotypes.

    let mut exact_clonotypes = find_exact_subclonotypes(&ctl, &tig_bc, &refdata);

    // Filter out some foursie artifacts.

    if ctl.clono_filt_opt.weak_foursies {
        let t = Instant::now();
        let mut to_delete = vec![false; exact_clonotypes.len()];
        let mut twosies = Vec::<(Vec<u8>, Vec<u8>)>::new();
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            if ex.share.len() == 2 && (ex.share[0].left ^ ex.share[1].left) && ex.ncells() >= 10 {
                twosies.push((ex.share[0].seq.clone(), ex.share[1].seq.clone()));
            }
        }
        unique_sort(&mut twosies);
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            if ex.share.len() == 4 {
                for i1 in 0..4 {
                    for i2 in i1 + 1..4 {
                        if ex.share[i1].left ^ ex.share[i2].left {
                            let p = (ex.share[i1].seq.clone(), ex.share[i2].seq.clone());
                            if bin_member(&twosies, &p) {
                                to_delete[i] = true;
                            }
                        }
                    }
                }
            }
        }
        erase_if(&mut exact_clonotypes, &to_delete);
        if ctl.comp {
            println!("used {:.2} seconds filtering foursies", elapsed(&t));
        }
    }

    // Look for insertions (experimental).

    find_insertions(&ctl, &exact_clonotypes);

    // Build info about clonotypes.  Note that this edits the V reference sequence to perform
    // an indel in some cases.

    let mut info: Vec<CloneInfo> = build_info(&refdata, &ctl, &mut exact_clonotypes);

    // Derive consensus sequences for alternate alleles of V segments.  Then create donor
    // reference sequences for Loupe.

    let alt_refs : Vec<(usize,usize,DnaString)>  // {(donor, ref id, alt seq)}
        = find_alleles( &refdata, &ctl, &exact_clonotypes );
    if ctl.gen_opt.dref_file.len() > 0 {
        let f = File::create(&ctl.gen_opt.dref_file);
        if f.is_err() {
            eprintln!(
                "\nProblem with DONOR_REF_FILE: unable to write to the file\n\
                 \"{}\".\nPlease check that that path makes sense and that you have write \
                 permission for it.\n",
                ctl.gen_opt.dref_file
            );
            std::process::exit(1);
        }
        let mut f = BufWriter::new(f.unwrap());
        let mut count = 0;
        for i in 0..alt_refs.len() {
            let donor = alt_refs[i].0;
            let ref_id = alt_refs[i].1;
            if i > 0 && (donor != alt_refs[i - 1].0 || ref_id != alt_refs[i - 1].1) {
                count = 0;
            }
            let alt_seq = &alt_refs[i].2;
            fwriteln!(
                f,
                ">{}.{}.{}, {}\n{}",
                refdata.id[ref_id],
                donor + 1,
                count + 1,
                refdata.name[ref_id],
                alt_seq.to_string()
            );
            count += 1;
        }
    }
    let tdonor = Instant::now();
    let drefs = make_donor_refs(&alt_refs, &refdata);
    if ctl.comp {
        println!("used {:.2} seconds making donor refs", elapsed(&tdonor));
    }

    // Update reference sequences for V segments by substituting in alt alleles if better.

    sub_alts(&ctl, &alt_refs, &mut info, &mut exact_clonotypes);

    // Form equivalence relation on exact subclonotypes.

    let mut join_info = Vec::<(usize, usize, bool, Vec<u8>)>::new();
    let eq: EquivRel = join_exacts(
        is_bcr,
        &refdata,
        &ctl,
        &exact_clonotypes,
        &info,
        &mut join_info,
    );
    /*
    if ctl.comp {
        if ctl.clono_filt_opt.ncells_low < ctl.clono_filt_opt.ncells_high {
            println!("");
        }
    }
    */

    // Lookup for heavy chain reuse (special purpose experimental option).

    lookup_heavy_chain_reuse(&ctl, &exact_clonotypes, &info, &eq);

    // Find and print clonotypes.

    let torb = Instant::now();
    print_clonotypes(
        &tall,
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &eq,
        &gex_info,
        &join_info,
    );
    if ctl.comp {
        if !ctl.gen_opt.noprint {
            println!("");
        }
        println!("used {:.2} seconds making orbits", elapsed(&torb));
    }

    // Report computational performance.

    if ctl.comp {
        println!(
            "\nused {:.2} seconds total, peak mem = {:.2} GB",
            elapsed(&tall),
            peak_mem_usage_gb()
        );
    }
    println!("");
    // It's not totally clear that the exit below actually saves time.  Would need more testing.
    if !ctl.gen_opt.cellranger {
        std::process::exit(0);
    }
}
