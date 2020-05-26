// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use vdj_ann::*;

use self::refx::*;
use debruijn::dna_string::DnaString;
use enclone::allele::*;
use enclone::explore::*;
use enclone::graph_filter::*;
use enclone::info::*;
use enclone::join::*;
use enclone::load_gex::*;
use enclone::misc1::*;
use enclone::misc2::*;
use enclone::misc3::*;
use enclone::proc_args::*;
use enclone::proc_args2::*;
use enclone::proc_args_check::*;
use enclone::read_json::*;
use enclone_core::defs::*;
use enclone_core::*;
use enclone_help::help1::*;
use enclone_help::help2::*;
use enclone_help::help3::*;
use enclone_help::help4::*;
use enclone_help::help5::*;
use enclone_help::help_utils::*;
use enclone_print::loupe::*;
use enclone_print::print_clonotypes::*;
use equiv::EquivRel;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use regex::Regex;
use serde_json::Value;
use stats_utils::*;
use std::{
    cmp::max,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn binomial_sum(n: usize, k: usize, p: f64) -> f64 {
    assert!(n >= 1);
    assert!(k <= n);
    let mut sum = 0.0;
    let mut choose = 1.0;
    for _ in 0..n {
        choose *= 1.0 - p;
    }
    let q = p / (1.0 - p);
    for i in 0..=k {
        sum += choose;
        choose *= (n - i) as f64;
        choose /= (i + 1) as f64;
        choose *= q;
    }
    sum
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    // Provide help if requested.

    {
        if args.len() == 2 && (args[1] == "version" || args[1] == "--version") {
            println!("{} : {}", env!("CARGO_PKG_VERSION"), version_string());
            std::process::exit(0);
        }
        let mut args = args.clone();
        let mut to_delete = vec![false; args.len()];
        let mut nopager = false;
        let mut plain = false;
        let mut long_help = false;
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                nopager = true;
                to_delete[i] = true;
            } else if args[i] == "HTML" {
                ctl.gen_opt.html = true;
                ctl.gen_opt.html_title = "enclone output".to_string();
                to_delete[i] = true;
            } else if args[i].starts_with("HTML=") {
                ctl.gen_opt.html = true;
                let mut title = args[i].after("HTML=").to_string();
                if title.starts_with("\"") && title.ends_with("\"") {
                    title = title.between("\"", "\"").to_string();
                }
                ctl.gen_opt.html_title = title;
                to_delete[i] = true;
            } else if args[i] == "SVG" {
                ctl.gen_opt.svg = true;
                to_delete[i] = true;
            } else if args[i] == "STABLE_DOC" {
                ctl.gen_opt.stable_doc = true;
                to_delete[i] = true;
            } else if args[i] == "FORCE_EXTERNAL" {
                to_delete[i] = true;
            } else if args[i] == "LONG_HELP" {
                long_help = true;
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PRE=") {
                to_delete[i] = true;
            } else if args[i] == "PLAIN" {
                to_delete[i] = true;
                plain = true;
                unsafe {
                    PLAIN = true;
                }
            }
        }
        if ctl.gen_opt.html && ctl.gen_opt.svg {
            eprintln!("\nBoth HTML and SVG cannot be used at the same time.\n");
            std::process::exit(1);
        }
        erase_if(&mut args, &to_delete);
        if args.len() == 1 || args.contains(&"help".to_string()) {
            PrettyTrace::new().on();
            setup_pager(!nopager);
        }
        let mut help_all = false;
        if args.len() >= 3 && args[1] == "help" && args[2] == "all" {
            unsafe {
                HELP_ALL = true;
            }
            help_all = true;
        }
        let mut h = HelpDesk::new(plain, help_all, long_help, ctl.gen_opt.html);
        help1(&args, &mut h);
        help2(&args, &ctl, &mut h);
        help3(&args, &mut h);
        help4(&args, &mut h);
        help5(&args, &ctl, &mut h);
    }

    // Pretest for some options.

    ctl.pretty = true;
    let mut nopretty = false;
    ctl.gen_opt.h5 = true;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN") {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NOPRETTY") {
            nopretty = true;
        }
        if is_simple_arg(&args[i], "COMP") {
            ctl.comp = true;
        }
        if is_simple_arg(&args[i], "COMP2") {
            ctl.comp = true;
            ctl.comp2 = true;
        }
        if is_simple_arg(&args[i], "CELLRANGER") {
            ctl.gen_opt.cellranger = true;
        }
        if is_simple_arg(&args[i], "NH5") {
            ctl.gen_opt.h5 = false;
        }
    }

    // Test for happening mode and turn on pretty trace.

    if !nopretty {
        let mut happening = 0;
        let mut ctrlc = false;
        for i in 1..args.len() {
            if args[i].starts_with("HAPS=") {
                // should actually test for usize
                happening = args[i].after("HAPS=").force_usize();
            }
            if is_simple_arg(&args[i], "CTRLC") {
                ctrlc = true;
            }
        }
        let thread_message = new_thread_message();
        if happening > 0 {
            PrettyTrace::new()
                .message(&thread_message)
                .profile(happening)
                .whitelist(&vec![
                    "amino",
                    "ansi_escape",
                    "binary_vec_io",
                    "enclone",
                    "equiv",
                    "graph_simple",
                    "io_utils",
                    "marsoc",
                    "mirror_sparse_matrix",
                    "perf_stats",
                    "stats_utils",
                    "stirling_numbers",
                    "string_utils",
                    "tables",
                    "vector_utils",
                ])
                .ctrlc()
                .on();
        } else if ctrlc {
            PrettyTrace::new().message(&thread_message).ctrlc().on();
        } else {
            let exit_message: String;
            if !ctl.gen_opt.cellranger {
                exit_message = format!(
                    "Something has gone badly wrong.  You have probably encountered an internal \
                    error in enclone.\n\n\
                    Please email us at enclone@10xgenomics.com, including the traceback shown\n\
                    above and also the following version information:\n\
                    {} : {}.\n\n\
                    Thank you and have a nice day!",
                    env!("CARGO_PKG_VERSION"),
                    version_string()
                );
            } else {
                exit_message = format!(
                    "Something has gone badly wrong.  You have probably \
                     encountered an internal error\nin cellranger.  \
                     Please email us at support@10xgenomics.com, including the traceback\nshown \
                     above."
                );
            }
            PrettyTrace::new().exit_message(&exit_message).on();
            let mut nopager = false;
            for i in 1..args.len() {
                if args[i] == "NOPAGER" {
                    nopager = true;
                }
            }
            setup_pager(!nopager);
        }
    }

    // Process args (and set defaults for them).

    proc_args(&mut ctl, &args);

    // Dump lenas.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "DUMP_INTERNAL_IDS") {
            let mut x = Vec::<usize>::new();
            for y in ctl.sample_info.dataset_id.iter() {
                x.push(y.force_usize());
            }
            x.sort();
            println!("\n{}\n", x.iter().format(","));
            std::process::exit(0);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone(args: &Vec<String>) {
    // Set up stuff, read args, etc.

    let tall = Instant::now();
    let (mut print_cpu, mut print_cpu_info) = (false, false);
    for i in 1..args.len() {
        if args[i] == "PRINT_CPU" {
            print_cpu = true;
        }
        if args[i] == "PRINT_CPU_INFO" {
            print_cpu_info = true;
        }
    }
    let (mut cpu_all_start, mut cpu_this_start) = (0, 0);
    if print_cpu || print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_start = t.before(" ").force_usize();
            break;
        }
        let f = open_for_read![&format!("/proc/{}/stat", std::process::id())];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            cpu_this_start = fields[13].force_usize();
        }
    }
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

    // Get gene expression and feature barcode counts.  Sanity check variables in cases where that
    // has to occur after loading GEX data.  This could also occur after loading only
    // the feature list, which would be better.

    let gex_info = get_gex_info(&mut ctl);
    check_lvars(&ctl, &gex_info);
    check_pcols(&ctl, &gex_info);

    // Find matching features for <regular expression>_g etc.

    ctl.clono_print_opt.regex_match =
        vec![HashMap::<String, Vec<usize>>::new(); ctl.sample_info.n()];
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let ends1 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g", "_ab", "_ag", "_cr", "_cu", "_g",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    let mut endsz = Vec::<String>::new();
    for (ix, x) in ends0.iter().enumerate() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
            endsz.push(ends1[ix].to_string());
        }
    }
    let mut vars = ctl.clono_print_opt.lvars.clone();
    vars.append(&mut ctl.parseable_opt.pcols.clone());
    unique_sort(&mut vars);
    for x in vars.iter() {
        for (iy, y) in ends.iter().enumerate() {
            let mut xc = x.clone();
            if x.ends_with("_cell") {
                xc = xc.rev_before("_cell").to_string();
            }
            if xc.ends_with(y) {
                let mut p = xc.rev_before(y);
                if p.contains(':') {
                    p = p.after(":");
                }
                let pp = format!("{}{}", p, endsz[iy]);
                if !p.is_empty() && Regex::new(&p).is_ok() {
                    let mut ok = true;
                    let mut px = false;
                    let b = p.as_bytes();
                    for i in 0..p.len() {
                        if !((b[i] >= b'A' && b[i] <= b'Z')
                            || (b[i] >= b'a' && b[i] <= b'z')
                            || (b[i] >= b'0' && b[i] <= b'9')
                            || b".-_[]()|*".contains(&b[i]))
                        {
                            ok = false;
                            break;
                        }
                        if b"[]()|*".contains(&b[i]) {
                            px = true;
                        }
                    }
                    if ok && px {
                        let reg = Regex::new(&format!("^{}$", p));
                        for li in 0..ctl.sample_info.n() {
                            let mut js = Vec::<usize>::new();
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
                                    js.push(j);
                                }
                            }
                            if js.len() > 0 {
                                ctl.clono_print_opt.regex_match[li].insert(pp.clone(), js);
                            }
                        }
                        let mut matches = false;
                        for li in 0..ctl.sample_info.n() {
                            if ctl.clono_print_opt.regex_match[li].contains_key(&pp) {
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

    // Determine the Cell Ranger version that was used.

    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
    let json = format!("{}/{}", ctl.sample_info.dataset_path[0], ann);
    let json_lz4 = format!("{}/{}.lz4", ctl.sample_info.dataset_path[0], ann);
    if !path_exists(&json) && !path_exists(&json_lz4) {
        eprintln!("\ncan't find {} or {}\n", json, json_lz4);
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

    // Find the VDJ reference.

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
            if ctl.gen_opt.descrip {
                println!("using old mouse reference");
            }
            refx = mouse_ref_old();
        } else {
            if ctl.gen_opt.descrip {
                println!("using new mouse reference");
            }
            refx = mouse_ref();
        }
    } else {
        if ctl.gen_opt.imgt && ctl.gen_opt.internal_run {
            let imgt =
                "/mnt/opt/refdata_cellranger/vdj/vdj_IMGT_human_20200415-0.0.0/fasta/regions.fa";
            if ctl.gen_opt.descrip {
                println!("using imgt human reference");
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
                refx += &s;
                refx += &"\n";
            }
            ctl.gen_opt.reannotate = true;
        } else if ctl.gen_opt.cr_version == "".to_string() && !ctl.gen_opt.reannotate {
            if ctl.gen_opt.descrip {
                println!("using old human reference");
            }
            refx = human_ref_old();
        } else {
            if ctl.gen_opt.descrip {
                println!("using new human reference");
            }
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

    /*

    // Remove V sequences not beginning with a start codon and do some tidying.
    // Commented out until proven useful.

    let lines = refx.split('\n').collect::<Vec<&str>>();
    let mut refx2 = String::new();
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with(">") {
                break;
            }
            j += 1;
        }
        let mut seq = String::new();
        for k in i + 1..j {
            seq += &lines[k];
        }
        seq = seq.replace('a', "A");
        seq = seq.replace('c', "C");
        seq = seq.replace('g', "G");
        seq = seq.replace('t', "T");
        let mut ok = true;
        if lines[i].contains("V-REGION") {
            if !seq.starts_with("ATG") {
                ok = false;
            }
        }
        if ok {
            refx2 += &format!("{}\n{}\n", lines[i], seq);
        }
        i = j;
    }

    */
    let refx2 = &refx;

    // Build reference data.

    make_vdj_ref_data_core(&mut refdata, &refx2, &ext_refx, is_tcr, is_bcr, None);
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
    let mut vdj_cells = Vec::<Vec<String>>::new();
    parse_json_annotations_files(
        &mut ctl,
        &mut tig_bc,
        &refdata,
        &to_ref_index,
        &mut vdj_cells,
    );

    // Search for SHM indels.

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

    // For B cells, filter based on UMI counts.  More details in heuristics.html.
    // Find all clonotypes having one cell which has two chains,
    // one heavy and one light.  Get the sum of the chain UMI counts for this cell.
    //
    // For each cell, let umish be the umi count for its heavy chain having the most umis, and
    // similarly define umisl.  Let umitot = umish + umisl.
    //
    // If every cell in a clonotype would have been deleted, first find the exact subclonotype for
    // which the sum of its umitot values is greatest, and then in it, find the cell having
    // highest umitot value.  Protect this cell, so long as it has at least two chains.

    let tumi = Instant::now();
    let mut orbits = Vec::<Vec<i32>>::new();
    let mut reps = Vec::<i32>::new();
    eq.orbit_reps(&mut reps);
    if is_tcr || !ctl.clono_filt_opt.umi_filt {
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            orbits.push(o);
        }
    }
    if is_bcr
        && (ctl.gen_opt.baseline || ctl.clono_filt_opt.umi_filt || ctl.clono_filt_opt.umi_filt_mark)
    {
        let mut umis = vec![Vec::<usize>::new(); ctl.sample_info.n()];
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            if o.solo() {
                let x: &CloneInfo = &info[o[0] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                if ex.ncells() == 1 && ex.share.duo() && ex.share[0].left != ex.share[1].left {
                    umis[ex.clones[0][0].dataset_index]
                        .push(ex.clones[0][0].umi_count + ex.clones[0][1].umi_count);
                }
            }
        }
        let mut nu = vec![0; ctl.sample_info.n()];
        let mut umin = vec![0.0; ctl.sample_info.n()];
        for l in 0..ctl.sample_info.n() {
            umis[l].sort();
            nu[l] = umis[l].len();
            if ctl.gen_opt.baseline {
                println!("\n{} umi counts for dataset {}", nu[l], l + 1);
            }
            if nu[l] > 0 {
                let n10 = umis[l][nu[l] / 10] as f64;
                let n50 = umis[l][nu[l] / 2] as f64;
                umin[l] = n10.min(n50 - (4.0 * n50.sqrt()));
            }
            if nu[l] > 0 && ctl.gen_opt.baseline {
                println!("1% ==> {}", umis[l][umis[l].len() / 100]);
                println!("2% ==> {}", umis[l][umis[l].len() / 50]);
                println!("5% ==> {}", umis[l][umis[l].len() / 20]);
                println!("10% ==> {}", umis[l][umis[l].len() / 10]);
                println!("20% ==> {}", umis[l][umis[l].len() / 5]);
                println!("50% ==> {}", umis[l][umis[l].len() / 2]);
                println!("umin = {:.2}", umin[l]);
            }
        }
        if ctl.clono_filt_opt.umi_filt || ctl.clono_filt_opt.umi_filt_mark {
            const MIN_BASELINE_CELLS: usize = 20;
            for i in 0..reps.len() {
                let mut o = Vec::<i32>::new();
                eq.orbit(reps[i], &mut o);
                let mut ncells = 0;
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &exact_clonotypes[x.clonotype_index];
                    ncells += ex.ncells();
                }
                let mut nbads = 0;
                if ncells >= 2 {
                    let mut to_deletex = vec![false; o.len()];
                    let (mut best_ex, mut best_ex_sum) = (0, 0);
                    let (mut best_cell, mut best_cell_count) = (0, 0);
                    let mut baselined = true;
                    let mut protected = false;
                    for pass in 1..=3 {
                        if pass == 2 {
                            if nbads == 0 {
                                protected = true;
                            } else {
                                let p = 0.1;
                                let bound = 0.01;

                                // Find probability of observing nbads or more events of probability
                                // p in a sample of size ncells, and if that is at least bound,
                                // don't delete any cells (except onesies).

                                if binomial_sum(ncells, ncells - nbads, 1.0 - p) >= bound {
                                    protected = true;
                                }
                            }
                        }
                        for j in 0..o.len() {
                            let x: &CloneInfo = &info[o[j] as usize];
                            let ex = &mut exact_clonotypes[x.clonotype_index];
                            let mut to_delete = vec![false; ex.ncells()];
                            let mut ex_sum = 0;
                            for k in 0..ex.ncells() {
                                let li = ex.clones[k][0].dataset_index;
                                if nu[li] >= MIN_BASELINE_CELLS {
                                    let (mut umish, mut umisl) = (0, 0);
                                    for l in 0..ex.share.len() {
                                        if ex.share[l].left {
                                            umish = max(umish, ex.clones[k][l].umi_count);
                                        } else {
                                            umisl = max(umish, ex.clones[k][l].umi_count);
                                        }
                                    }
                                    let umitot = umish + umisl;
                                    if pass == 1 {
                                        ex_sum += umitot;
                                    }
                                    if pass == 2
                                        && j == best_ex
                                        && umitot > best_cell_count
                                        && ex.share.len() > 1
                                    {
                                        best_cell = k;
                                        best_cell_count = umitot;
                                    }
                                    if (umitot as f64) < umin[li] {
                                        if pass == 1 {
                                            nbads += 1;
                                        } else if pass == 3 && protected {
                                            if ex.share.len() == 1 {
                                                to_delete[k] = true;
                                                if ctl.clono_filt_opt.umi_filt_mark {
                                                    ex.clones[k][0].marked = true;
                                                }
                                            }
                                        } else if pass == 3 {
                                            if !baselined
                                                || (best_ex, best_cell) != (j, k)
                                                || ex.share.len() == 1
                                            {
                                                to_delete[k] = true;
                                                if ctl.clono_filt_opt.umi_filt_mark {
                                                    ex.clones[k][0].marked = true;
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    baselined = false;
                                }
                            }
                            if pass == 1 && ex_sum > best_ex_sum {
                                best_ex = j;
                                best_ex_sum = ex_sum;
                            }
                            if pass == 3 && ctl.clono_filt_opt.umi_filt {
                                erase_if(&mut ex.clones, &to_delete);
                            }
                        }
                    }
                    for j in 0..o.len() {
                        let x: &CloneInfo = &info[o[j] as usize];
                        let ex = &mut exact_clonotypes[x.clonotype_index];
                        if ex.ncells() == 0 {
                            to_deletex[j] = true;
                        }
                    }
                    erase_if(&mut o, &to_deletex);
                }
                if ctl.clono_filt_opt.umi_filt && !o.is_empty() {
                    orbits.push(o.clone());
                }
            }
        }
    }

    // Filter B cells based on UMI count ratios.  For now just marking.  This assumes V..J
    // identity to filter.

    if ctl.clono_filt_opt.umi_ratio_filt || ctl.clono_filt_opt.umi_ratio_filt_mark {
        const MIN_UMI_RATIO: usize = 500;
        let mut orbits2 = Vec::<Vec<i32>>::new();
        'orbit: for i in 0..orbits.len() {
            let mut ncells = 0;
            let mut o = orbits[i].clone();
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                ncells += ex.ncells();
            }
            let mut nbads = 0;
            for pass in 1..=2 {
                if pass == 2 {
                    if nbads == 0 {
                        orbits2.push(o.clone());
                        continue 'orbit;
                    } else {
                        let p = 0.1;
                        let bound = 0.01;

                        // Find probability of observing nbads or more events of probability
                        // p in a sample of size ncells, and if that is at least bound,
                        // don't delete any cells.

                        if binomial_sum(ncells, ncells - nbads, 1.0 - p) >= bound {
                            orbits2.push(o.clone());
                            continue 'orbit;
                        }
                    }
                }
                let mut to_deletex = vec![false; o.len()];
                let mut z = Vec::<(Vec<u8>, usize, usize, usize, usize)>::new();
                let mut to_delete = Vec::<Vec<bool>>::new();
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &mut exact_clonotypes[x.clonotype_index];
                    to_delete.push(vec![false; ex.ncells()]);
                    for k in 0..ex.ncells() {
                        let mut tot = 0;
                        for m in 0..ex.clones[k].len() {
                            tot += ex.clones[k][m].umi_count;
                        }
                        for m in 0..ex.clones[k].len() {
                            z.push((
                                ex.share[m].seq.clone(),
                                ex.clones[k][m].umi_count,
                                j,
                                k,
                                tot,
                            ));
                        }
                    }
                }
                reverse_sort(&mut z);
                let mut j = 0;
                while j < z.len() {
                    let k = next_diff1_5(&z, j as i32) as usize;
                    for l in j..k {
                        if z[j].1 >= MIN_UMI_RATIO * z[l].4 {
                            to_delete[z[l].2][z[l].3] = true;
                        }
                    }
                    j = k;
                }
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &mut exact_clonotypes[x.clonotype_index];
                    for l in 0..ex.ncells() {
                        if to_delete[j][l] {
                            if ctl.clono_filt_opt.umi_ratio_filt_mark {
                                ex.clones[l][0].marked = true;
                            }
                            nbads += 1;
                        }
                    }
                    if pass == 2 && ctl.clono_filt_opt.umi_ratio_filt {
                        erase_if(&mut ex.clones, &to_delete[j]);
                        if ex.ncells() == 0 {
                            to_deletex[j] = true;
                        }
                    }
                }
                if pass == 2 {
                    if ctl.clono_filt_opt.umi_ratio_filt {
                        erase_if(&mut o, &to_deletex);
                        if !o.is_empty() {
                            orbits2.push(o.clone());
                        }
                    }
                }
            }
        }
        if ctl.clono_filt_opt.umi_ratio_filt {
            orbits = orbits2;
        }
    }

    // Remove cells that are not called cells by GEX or feature barcodes.

    if !ctl.clono_filt_opt.ngex {
        let mut orbits2 = Vec::<Vec<i32>>::new();
        for i in 0..orbits.len() {
            let mut o = orbits[i].clone();
            let mut to_deletex = vec![false; o.len()];
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                let mut to_delete = vec![false; ex.ncells()];
                for k in 0..ex.ncells() {
                    let li = ex.clones[k][0].dataset_index;
                    let bc = &ex.clones[k][0].barcode;
                    if ctl.sample_info.gex_path[li].len() > 0 {
                        let gbc = &gex_info.gex_cell_barcodes[li];
                        if !bin_member(&gbc, &bc) {
                            to_delete[k] = true;
                        }
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
                if ex.ncells() == 0 {
                    to_deletex[j] = true;
                }
            }
            erase_if(&mut o, &to_deletex);
            if !o.is_empty() {
                orbits2.push(o.clone());
            }
        }
        orbits = orbits2;
    }
    if ctl.comp {
        println!("used {:.2} seconds umi filtering and such", elapsed(&tumi));
    }

    // Check for disjoint orbits.  This is an incomplete test.

    let mut orbits2 = Vec::<Vec<i32>>::new();
    for i in 0..orbits.len() {
        let o = orbits[i].clone();
        let mut eqx = EquivRel::new(o.len() as i32);
        for i1 in 0..o.len() {
            for i2 in i1 + 1..o.len() {
                if eqx.class_id(i1 as i32) != eqx.class_id(i2 as i32) {
                    let x1: &CloneInfo = &info[o[i1] as usize];
                    let x2: &CloneInfo = &info[o[i2] as usize];
                    if x1.clonotype_index == x2.clonotype_index {
                        eqx.join(i1 as i32, i2 as i32);
                    } else {
                        let ex1 = &exact_clonotypes[x1.clonotype_index];
                        let ex2 = &exact_clonotypes[x2.clonotype_index];
                        'cloop: for m1 in 0..ex1.nchains() {
                            for m2 in 0..ex2.nchains() {
                                if ex1.share[m1].seq_del.len() == ex2.share[m2].seq_del.len()
                                    && ex1.share[m1].cdr3_aa.len() == ex2.share[m2].cdr3_aa.len()
                                {
                                    eqx.join(i1 as i32, i2 as i32);
                                    break 'cloop;
                                }
                            }
                        }
                    }
                }
            }
        }
        if eqx.norbits() == 1 {
            orbits2.push(o.clone());
        } else {
            let mut repsx = Vec::<i32>::new();
            eqx.orbit_reps(&mut repsx);
            for j in 0..repsx.len() {
                let mut ox = Vec::<i32>::new();
                eqx.orbit(repsx[j], &mut ox);
                let mut o2 = Vec::<i32>::new();
                for k in 0..ox.len() {
                    o2.push(o[ox[k] as usize]);
                }
                orbits2.push(o2);
            }
        }
    }
    orbits = orbits2;

    // Find and print clonotypes.

    let torb = Instant::now();
    print_clonotypes(
        &tall,
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &orbits,
        &gex_info,
        &join_info,
        &vdj_cells,
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

    let (mut cpu_all_stop, mut cpu_this_stop) = (0, 0);
    if print_cpu || print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_stop = t.before(" ").force_usize();
            break;
        }
        let f = open_for_read![&format!("/proc/{}/stat", std::process::id())];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            cpu_this_stop = fields[13].force_usize();
        }
        let (this_used, all_used) = (cpu_this_stop - cpu_this_start, cpu_all_stop - cpu_all_start);
        if print_cpu {
            println!("{}", this_used);
        } else {
            println!(
                "used cpu = {} = {:.1}% of total",
                this_used,
                percent_ratio(this_used, all_used)
            );
        }
    }

    println!("");
    // It's not totally clear that the exit below actually saves time.  Would need more testing.
    if !ctl.gen_opt.cellranger {
        std::process::exit(0);
    }
}
