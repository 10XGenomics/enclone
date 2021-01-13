// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use self::refx::*;
use amino::*;
use debruijn::dna_string::DnaString;
use enclone::allele::*;
use enclone::explore::*;
use enclone::fwr3_freqs::*;
use enclone::graph_filter::*;
use enclone::info::*;
use enclone::innate::*;
use enclone::join::*;
use enclone::load_gex::*;
use enclone::misc1::*;
use enclone::misc2::*;
use enclone::misc3::*;
use enclone::proc_args::*;
use enclone::proc_args2::*;
use enclone::proc_args_check::*;
use enclone::read_json::*;
use enclone::secret::*;
use enclone_core::defs::*;
use enclone_core::vdj_features::*;
use enclone_core::*;
use enclone_help::help1::*;
use enclone_help::help2::*;
use enclone_help::help3::*;
use enclone_help::help4::*;
use enclone_help::help5::*;
use enclone_help::help_utils::*;
use enclone_print::loupe::*;
use enclone_print::print_clonotypes::*;
use enclone_print::print_utils4::*;
use enclone_tail::tail::tail_code;
use equiv::EquivRel;
use evalexpr::*;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use regex::Regex;
use serde_json::Value;
use stats_utils::*;
use std::{
    cmp::max,
    collections::HashMap,
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    process::Command,
    thread, time,
    time::Instant,
};
use string_utils::*;
use tables::*;
use vdj_ann::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    let t = Instant::now();
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

    if !nopretty && !ctl.gen_opt.cellranger {
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
            let args: Vec<String> = env::args().collect();
            let exit_message = format!(
                "Something has gone badly wrong.  You have probably encountered an internal \
                error in enclone.\n\n\
                Please email us at enclone@10xgenomics.com, including the traceback shown\n\
                above and also the following version information:\n\
                {} : {}.\n\n\
                Your command was:\n\n{}\n\n\
                🌸 Thank you and have a nice day! 🌸",
                env!("CARGO_PKG_VERSION"),
                version_string(),
                args.iter().format(" "),
            );
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
    ctl.perf_stats(&t, "in first part of setup");

    // Process args (and set defaults for them).

    proc_args(&mut ctl, &args);

    // Dump lenas.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "DUMP_INTERNAL_IDS") {
            let mut x = Vec::<usize>::new();
            for y in ctl.origin_info.dataset_id.iter() {
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
    let (mut comp, mut comp2) = (false, false);
    for i in 1..args.len() {
        if args[i] == "PRINT_CPU" {
            print_cpu = true;
        }
        if args[i] == "PRINT_CPU_INFO" {
            print_cpu_info = true;
        }
        if args[i] == "COMP" {
            comp = true;
        }
        if args[i] == "COMP2" {
            comp2 = true;
        }
    }
    if comp && !comp2 {
        println!("");
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
    ctl.perf_stats(&tall, "before setup");
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
    let twoof = Instant::now();
    check_lvars(&ctl, &gex_info);
    check_pcols(&ctl, &gex_info);

    // Find matching features for <regular expression>_g etc.

    ctl.clono_print_opt.regex_match =
        vec![HashMap::<String, Vec<usize>>::new(); ctl.origin_info.n()];
    let ends0 = [
        "_g", "_ab", "_cr", "_cu", "_g_μ", "_ab_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let ends1 = ["_g", "_ab", "_cr", "_cu", "_g", "_ab", "_cr", "_cu", "_g"];
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
                        for li in 0..ctl.origin_info.n() {
                            let mut js = Vec::<usize>::new();
                            for j in 0..gex_info.gex_features[li].len() {
                                let f = &gex_info.gex_features[li][j];
                                let ff = f.split('\t').collect::<Vec<&str>>();
                                let mut ok = false;
                                if ff[2].starts_with("Antibody") {
                                    if y.contains("_ab") {
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
                        for li in 0..ctl.origin_info.n() {
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
    ctl.perf_stats(&twoof, "doing miscellaneous stuff");

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Start of code to determine the reference sequence that is to be used.  First check for the
    // existence of a json file.

    let tr = Instant::now();
    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
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

    // Step 1.  Test to see if CURRENT_REF or BUILT_IN is specified.  Kind of a mess that we
    // have both.

    let mut refx = String::new();
    if ctl.gen_opt.current_ref || ctl.gen_opt.built_in {
        if !ctl.gen_opt.mouse {
            refx = human_ref();
        } else {
            refx = mouse_ref();
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
                refx += &s;
                refx += &"\n";
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
                refx += &s;
                refx += &"\n";
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
            refx = refs[0].clone();
        }
    }

    // Step 5.  Attempt to determine the reference that was used by reading far enough into the
    // first json file to find a distinguishing entry.

    if refx.len() == 0 {
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
                    refx = refhash[matches[0]].1.clone();
                    break 'json_entry;
                }
            }
        }
    }
    if refx.len() == 0 {
        eprintln!(
            "\nenclone was unable to determine the reference sequence that you used.  You \
            have two options:\n\
            1. If you used cellranger version 4.0 or later, copy the vdj_reference directory\n   \
               from there to the outs directory that contains your other enclone input data.\n\
            2. Use the REF argument to specify the name of the reference fasta file.\n",
        );
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

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

    // Build reference data.

    let refx2 = &refx;
    let mut refdata = RefData::new();
    let ext_refx = String::new();
    let (mut is_tcr, mut is_bcr) = (true, true);
    if ctl.gen_opt.tcr {
        is_bcr = false;
    }
    if ctl.gen_opt.bcr {
        is_tcr = false;
    }
    make_vdj_ref_data_core(&mut refdata, &refx2, &ext_refx, is_tcr, is_bcr, None);
    let mut to_ref_index = HashMap::<usize, usize>::new();
    for i in 0..refdata.refs.len() {
        to_ref_index.insert(refdata.id[i] as usize, i);
    }

    // Flag defective reference sequences.

    let mut log = Vec::<u8>::new();
    let mut count = 0;
    let mut broken = vec![false; refdata.refs.len()];
    for i in 0..refdata.refs.len() {
        // Determine chain type and exclude those other than IGH, IGK, IGL, TRA and TRB.

        let rtype = refdata.rtype[i];
        let chain_type;
        if rtype == 0 {
            chain_type = "IGH";
        } else if rtype == 1 {
            chain_type = "IGK";
        } else if rtype == 2 {
            chain_type = "IGL";
        } else if rtype == 3 {
            chain_type = "TRA";
        } else if rtype == 4 {
            chain_type = "TRB";
        } else {
            continue;
        }

        // Look for problems.

        if refdata.is_c(i) {
            // This is very ugly.  We are exempting mouse IGHG2B because is is in our current
            // reference but has an extra base at the beginning.  See also comments below at
            // TRBV21-1.  Also, we're not actually checking for mouse.

            if refdata.name[i] == "IGHG2B" {
                continue;
            }

            // Continue.

            let seq = refdata.refs[i].to_ascii_vec();
            let aa0 = aa_seq(&seq, 0);
            let aa2 = aa_seq(&seq, 2);
            if aa2.contains(&b'*') && !aa0.contains(&b'*') {
                count += 1;
                broken[i] = true;
                fwriteln!(
                    log,
                    "{}. The following C segment reference sequence appears to have \
                    an extra base at its beginning:\n",
                    count
                );
                fwriteln!(log, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
            }
        } else if refdata.is_v(i) {
            // This is very ugly.  We are exempting human TRBV21-1 because it is in our current
            // reference (twice), but has multiple stop codons.  It should be deleted from the
            // reference, and then we should remove this test.  But probably in the future so as
            // not to inconvenience users.

            if ctl.gen_opt.species != "mouse" && refdata.name[i] == "TRBV21-1" {
                broken[i] = true;
                continue;
            }

            // This is very ugly.  We are exempting human IGHV1-12 because it is in our current
            // human reference, but is only 60 amino acids long.  Also we're not checking for mouse.

            if refdata.name[i] == "IGHV1-12" {
                broken[i] = true;
                continue;
            }

            // Ugly.  Frameshifted.  Also this is mouse and not checking for that.

            if refdata.name[i] == "TRAV23" {
                broken[i] = true;
                continue;
            }

            // Ugly.  Truncated on right.  Human.

            if refdata.name[i] == "IGLV5-48" {
                broken[i] = true;
                continue;
            }

            // Test for broken.

            let seq = refdata.refs[i].to_ascii_vec();
            let aa = aa_seq(&seq, 0);
            let mut reasons = Vec::<String>::new();
            if !aa.starts_with(b"M") {
                reasons.push("does not begin with a start codon".to_string());
            }
            let stops = aa.iter().filter(|&n| *n == b'*').count();
            if stops > 1 {
                reasons.push("has more than one stop codon".to_string());
            }
            if aa.len() < 100 {
                reasons.push("appears truncated (has less than 100 amino acids)".to_string());
            } else if aa.len() < 105 && chain_type == "IGH" {
                reasons.push("appears truncated (has less than 105 amino acids)".to_string());
            }
            if aa.len() >= 30 {
                let mut aap = aa.clone();
                aap.push(b'C');
                if cdr3_score(&aap, &chain_type, false) > 4 + cdr3_score(&aa, &chain_type, false) {
                    reasons.push("appears to need a C to be appended to its right end".to_string());
                }
            }
            if stops > 0 {
                let mut fixable = false;
                const TRIM: usize = 10;
                for j in 0..aa.len() - TRIM {
                    if aa[j] == b'*' {
                        let mut seqx = seq.clone();
                        for _ in 1..=2 {
                            let _ = seqx.remove(3 * j);
                            let aax = aa_seq(&seqx, 0);
                            if !aax.contains(&b'*') {
                                fixable = true;
                            }
                        }
                    }
                }
                if fixable {
                    reasons.push("appears to be frameshifted".to_string());
                }
            }
            if aa.len() >= 31 {
                for del in 1..=2 {
                    let aad = aa_seq(&seq, del);
                    if cdr3_score(&aad, &chain_type, false)
                        > 4 + cdr3_score(&aa, &chain_type, false)
                    {
                        reasons.push("appears to be frameshifted".to_string());
                    }
                }
            }
            if reasons.is_empty() {
                let cs2 = cdr2_start(&aa, &chain_type, false);
                if cs2.is_some() {
                    let fr3 = fr3_start(&aa, &chain_type, false).unwrap();
                    if cs2.unwrap() > fr3 {
                        reasons.push(
                            "appears to be defective, because our computed \
                            CDR2 start exceeds our computed FWR3 start"
                                .to_string(),
                        );
                    }
                }
            }
            if reasons.is_empty() {
                if aa.len() >= 31 {
                    // Pretty crappy frameshift test.  One should see high aa and dna similarity
                    // to other seqs if shifted.  Or use more aas.
                    let score = cdr3_score(&aa, &chain_type, false);
                    let mut frameshift = false;
                    for del in 1..=2 {
                        let aad = aa_seq(&seq, del);
                        if score <= 6 && cdr3_score(&aad, &chain_type, false) >= 3 + score {
                            frameshift = true;
                        }
                    }
                    if frameshift {
                        reasons.push("appears to be frameshifted".to_string());
                    }
                }
            }
            if reasons.is_empty() {
                let r;
                if chain_type == "IGH" {
                    r = 0;
                } else if chain_type == "IGK" {
                    r = 1;
                } else if chain_type == "IGL" {
                    r = 2;
                } else if chain_type == "TRA" {
                    r = 3;
                } else {
                    assert_eq!(chain_type, "TRB");
                    r = 4;
                }
                let freqs = fwr3_freqs();
                let score = score_fwr3(&aa, r, &freqs);
                if score < 8.0 && score4(&aa, r) < 5 {
                    reasons.push("appears to be frameshifted or truncated".to_string());
                }
            }

            // Report results.

            unique_sort(&mut reasons);
            if !reasons.is_empty() {
                let msg = format!(
                    "The following V segment reference sequence {}",
                    reasons.iter().format(", and ")
                );
                count += 1;
                broken[i] = true;
                fwriteln!(log, "{}. {}:\n", count, msg);
                fwriteln!(log, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
            }
        }
    }

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata);

    // Test for okness of sec/mem args.

    let mut vars = ctl.parseable_opt.pcols.clone();
    vars.append(&mut ctl.clono_print_opt.lvars.clone());
    unique_sort(&mut vars);
    ctl.gen_opt.using_secmem =
        bin_member(&vars, &"sec".to_string()) || bin_member(&vars, &"mem".to_string());
    if !ctl.gen_opt.using_secmem
        && ctl.parseable_opt.pout.len() > 0
        && ctl.parseable_opt.pcols.len() == 0
    {
        if ctl.gen_opt.species == "human" || ctl.gen_opt.species == "mouse" {
            if is_bcr {
                let mut have_bam = true;
                for g in ctl.origin_info.gex_path.iter() {
                    if g.len() == 0 {
                        have_bam = false;
                        break;
                    }
                    let bam = format!("{}/possorted_genome_bam.bam", g);
                    if !path_exists(&bam) {
                        have_bam = false;
                        break;
                    }
                }
                if have_bam {
                    let o = Command::new("samtools")
                        .arg("--help")
                        .output()
                        .expect("failed to execute samtools");
                    let status = o.status.code().unwrap();
                    if status == 0 {
                        ctl.gen_opt.using_secmem = true;
                    }
                }
            }
        }
    }
    if bin_member(&vars, &"sec".to_string()) || bin_member(&vars, &"mem".to_string()) {
        if ctl.gen_opt.species != "human" && ctl.gen_opt.species != "mouse" {
            eprintln!("\nThe lvars sec and mem can only be used for data from human and mouse.\n");
            std::process::exit(1);
        }
        if !is_bcr {
            eprintln!("\nThe lvars sec and mem do not make sense for TCR data.\n");
            std::process::exit(1);
        }
        for g in ctl.origin_info.gex_path.iter() {
            if g.len() == 0 {
                eprintln!("\nThe lvars sec and mem can only be used if GEX data are provided.\n");
                std::process::exit(1);
            }
            let bam = format!("{}/possorted_genome_bam.bam", g);
            if !path_exists(&bam) {
                eprintln!(
                    "\nThe lvars sec and mem can only be used if the file\n\
                    pos_sorted_genome_bam.bam is provided.  We did not see it at this path\n\
                    {}.",
                    g
                );
                std::process::exit(1);
            }
        }
        let o = Command::new("samtools")
            .arg("--help")
            .output()
            .expect("failed to execute samtools");
        let status = o.status.code().unwrap();
        if status != 0 {
            eprintln!(
                "\nThe lvars sec and mem can only be used if the samtools\n\
                executable is in your path.\n"
            );
            std::process::exit(1);
        }
    }
    ctl.perf_stats(&tr, "building reference and other things");

    // If sec (secreted) or mem (membrane) lvars have been specified, gather those data.

    if ctl.gen_opt.using_secmem {
        fetch_secmem(&mut ctl);
    }

    // Parse the json annotations file.

    let mut tig_bc = Vec::<Vec<TigData>>::new();
    let mut vdj_cells = Vec::<Vec<String>>::new();
    let mut gex_cells = Vec::<Vec<String>>::new();
    let mut gex_cells_specified = Vec::<bool>::new();
    let tparse = Instant::now();
    parse_json_annotations_files(
        &mut ctl,
        &mut tig_bc,
        &refdata,
        &to_ref_index,
        &mut vdj_cells,
        &mut gex_cells,
        &mut gex_cells_specified,
    );
    ctl.perf_stats(&tparse, "loading from json");

    // Populate features.

    let mut fr1_starts = vec![0; refdata.refs.len()];
    let mut fr2_starts = vec![None; refdata.refs.len()];
    let mut fr3_starts = vec![None; refdata.refs.len()];
    let mut cdr1_starts = vec![None; refdata.refs.len()];
    let mut cdr2_starts = vec![None; refdata.refs.len()];
    let mut fail = false;
    for i in 0..refdata.refs.len() {
        if refdata.is_v(i) {
            if broken[i] && ctl.gen_opt.require_unbroken_ok {
                continue;
            }
            let aa = aa_seq(&refdata.refs[i].to_ascii_vec(), 0);
            let rtype = refdata.rtype[i];
            let chain_type;
            if rtype == 0 {
                chain_type = "IGH";
            } else if rtype == 1 {
                chain_type = "IGK";
            } else if rtype == 2 {
                chain_type = "IGL";
            } else if rtype == 3 {
                chain_type = "TRA";
            } else if rtype == 4 {
                chain_type = "TRB";
            } else {
                continue;
            }
            let fs1 = fr1_start(&aa, &chain_type);
            fr1_starts[i] = 3 * fs1;
            let fs2 = fr2_start(&aa, &chain_type, false);
            if fs2.is_some() {
                fr2_starts[i] = Some(3 * fs2.unwrap());
            } else if ctl.gen_opt.require_unbroken_ok {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR2 start \
                    could not be computed\nfor this reference sequence:"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
            let fs3 = fr3_start(&aa, &chain_type, false);
            if fs3.is_some() {
                fr3_starts[i] = Some(3 * fs3.unwrap());
            } else if ctl.gen_opt.require_unbroken_ok {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR3 start \
                    could not be computed\nfor this reference sequence:"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
            let cs1 = cdr1_start(&aa, &chain_type, false);
            if cs1.is_some() {
                cdr1_starts[i] = Some(3 * cs1.unwrap());
                if fs2.is_some() && cs1.unwrap() > fs2.unwrap() && ctl.gen_opt.require_unbroken_ok {
                    eprintln!(
                        "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR1 start \
                        exceeds the FWR2 start for this reference sequence:"
                    );
                    let seq = refdata.refs[i].to_ascii_vec();
                    eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                    fail = true;
                }
            } else if ctl.gen_opt.require_unbroken_ok {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR1 start \
                    could not be computed\nfor this reference sequence:\n"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
            let cs2 = cdr2_start(&aa, &chain_type, false);
            if cs2.is_some() {
                cdr2_starts[i] = Some(3 * cs2.unwrap());
                if ctl.gen_opt.require_unbroken_ok && fs3.is_some() && cs2.unwrap() > fs3.unwrap() {
                    eprintln!(
                        "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR2 start \
                        exceeds the FWR3 start for this reference sequence:"
                    );
                    let seq = refdata.refs[i].to_ascii_vec();
                    eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                    fail = true;
                }
            } else if ctl.gen_opt.require_unbroken_ok {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR2 start \
                    could not be computed\nfor this reference sequence:"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
            if cs1.is_some() && fs1 > cs1.unwrap() && ctl.gen_opt.require_unbroken_ok {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR1 start \
                    exceeds the CDR1 start for this reference sequence:\n"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
            if cs2.is_some()
                && fs2.is_some()
                && fs2.unwrap() > cs2.unwrap()
                && ctl.gen_opt.require_unbroken_ok
            {
                eprintln!(
                    "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR2 start \
                    exceeds the CDR2 start for this reference sequence:"
                );
                let seq = refdata.refs[i].to_ascii_vec();
                eprintln!(">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
                fail = true;
            }
        }
    }
    if fail {
        std::process::exit(1);
    }

    // Report on broken reference sequences.  This comes after the json loading because possibly
    // the user supplied the wrong reference, so there is no value in criticizing the reference
    // in that case.

    if !log.is_empty() && !ctl.gen_opt.cellranger && !ctl.gen_opt.accept_broken {
        eprintln!(
            "\nSome errors were detected in the reference sequences supplied to enclone.\n\
            Please see comments at end for what you can do about this.\n",
        );
        eprint!("{}", strme(&log));
        eprintln!(
"🌼  Dear user, some defects were detected in the reference sequences supplied to enclone.   🌼\n\
 🌼  Some of these defects may be small.  Generally they are associated with V segments that 🌼\n\
 🌼  are frameshifted or truncated, or with C segments that have an extra base at the        🌼\n\
 🌼  beginning.  We are letting you know about this because they could result in             🌼\n\
 🌼  misannotation.                                                                          🌼\n"
        );

        let mut rows = Vec::<Vec<String>>::new();
        rows.push(vec![
            "You can make enclone ignore these defects by adding the additional argument"
                .to_string(),
        ]);
        rows.push(vec![
            "ACCEPT_BROKEN to the enclone command line.  Or you can obtain the same".to_string(),
        ]);
        rows.push(vec![
            "behavior by defining the environment variable ENCLONE_ACCEPT_BROKEN.".to_string(),
        ]);

        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l".to_vec(), false, true);
        eprintln!("{}", log);

        eprintln!(
        "This is probably OK, but if your sample is human or mouse, you may wish to either:\n\
        • rerun cellranger using the cleaned up reference sequences that come prepackaged with \
          it\n  (noting that your might have used an older, less clean version of that)\n\
        • or add the argument BUILT_IN to enclone, which will force use of the built-in reference\n  \
        sequences.  This will be a bit slower because all the contigs will need to be\n  \
        reannotated.  If you're using mouse, you'll also need to add the argument MOUSE.\n"
        );
        std::process::exit(1);
    }

    if ctl.gen_opt.require_unbroken_ok {
        std::process::exit(0);
    }
    for i in 0..tig_bc.len() {
        for j in 0..tig_bc[i].len() {
            let x = &mut tig_bc[i][j];
            x.fr1_start = fr1_starts[x.v_ref_id];
            x.fr2_start = fr2_starts[x.v_ref_id];
            x.fr3_start = fr3_starts[x.v_ref_id];
            x.cdr1_start = cdr1_starts[x.v_ref_id];
            x.cdr2_start = cdr2_starts[x.v_ref_id];
        }
    }

    // Search for SHM indels.

    let tproto = Instant::now();
    search_for_shm_indels(&ctl, &tig_bc);

    // Filter using light --> heavy graph.

    let mut fate = vec![HashMap::<String, String>::new(); vdj_cells.len()];
    graph_filter(&ctl, &mut tig_bc, ctl.gen_opt.graph, &mut fate);

    // Sort tig_bc.

    sort_tig_bc(&ctl, &mut tig_bc, &refdata);

    // Cross filter.

    cross_filter(&ctl, &mut tig_bc, &mut fate);

    // Look for barcode reuse.

    check_for_barcode_reuse(&ctl, &tig_bc);
    ctl.perf_stats(&tproto, "in proto stuff");

    // Find exact subclonotypes.

    let texact = Instant::now();
    let mut exact_clonotypes = find_exact_subclonotypes(&ctl, &tig_bc, &refdata, &mut fate);
    ctl.perf_stats(&texact, "finding exact subclonotypes");

    // Test for consistency between VDJ cells and GEX cells.  This is designed to work even if
    // NCELL is used.  We take up to 100 VDJ cells having both heavy and light (or TRB and TRA)
    // chains, and having the highest VDJ UMI count total (but using only one cell per exact
    // subclonotype), and find those that are GEX cells.
    //
    // If n cells were taken, and k of those are GEX cells, we require that
    // binomial_sum(n, k, 0.7) >= 0.00002.  For n = 100, this is the same as requiring that
    // k >= 50.  Using a binomial sum threshold allows the stringency of the requirement to be
    // appropriately lower when n is small.  When we tested on 260 libraries, the lowest value
    // observed for k/n was 0.65, and the vast majority of values were 0.9 or higher.
    //
    // This code is inefficient because for every dataset, it searches the entirety of tig_bc, but
    // it doesn't matter much because not much time is spent here.

    let tinc = Instant::now();
    let mut fail = false;
    for li in 0..ctl.origin_info.n() {
        if ctl.origin_info.gex_path[li].len() > 0 && !ctl.gen_opt.allow_inconsistent {
            let vdj = &vdj_cells[li];
            let gex = &gex_info.gex_cell_barcodes[li];
            let (mut heavy, mut light) = (vec![false; vdj.len()], vec![false; vdj.len()]);
            let mut exid = vec![0; vdj.len()];
            let mut inex = vec![false; vdj.len()];
            for i in 0..exact_clonotypes.len() {
                let ex = &exact_clonotypes[i];
                for j in 0..ex.clones.len() {
                    let p = bin_position(&vdj, &ex.clones[j][0].barcode);
                    if p >= 0 {
                        inex[p as usize] = true;
                        exid[p as usize] = i;
                    }
                }
            }
            let mut numi = vec![0; vdj.len()];
            for i in 0..tig_bc.len() {
                if tig_bc[i][0].dataset_index == li {
                    let p = bin_position(&vdj, &tig_bc[i][0].barcode);
                    if p >= 0 {
                        for j in 0..tig_bc[i].len() {
                            numi[p as usize] += tig_bc[i][j].umi_count;
                            if tig_bc[i][j].left {
                                heavy[p as usize] = true;
                            } else {
                                light[p as usize] = true;
                            }
                        }
                    }
                }
            }
            let mut x = Vec::<(usize, bool, usize)>::new();
            for i in 0..vdj.len() {
                if heavy[i] && light[i] {
                    x.push((numi[i], bin_member(&gex, &vdj[i]), i));
                }
            }
            reverse_sort(&mut x);
            let mut used = vec![false; exact_clonotypes.len()];
            let (mut total, mut good) = (0, 0);
            for i in 0..x.len() {
                let m = x[i].2;
                if inex[m] && used[exid[m]] {
                    continue;
                }
                total += 1;
                if x[i].1 {
                    good += 1;
                }
                if inex[m] {
                    used[exid[m]] = true;
                }
                if total == 100 {
                    break;
                }
            }
            if total >= 1 {
                let bino = binomial_sum(total, good, 0.7);
                if bino < 0.00002 {
                    fail = true;
                    eprint!(
                        "\nThe VDJ dataset with path\n{}\nand the GEX dataset with path\n\
                        {}\nshow insufficient sharing of barcodes.  ",
                        ctl.origin_info.dataset_path[li], ctl.origin_info.gex_path[li],
                    );
                    eprintln!(
                        "Of the {} VDJ cells that were tested,\nonly {} were GEX cells.",
                        total, good
                    );
                }
            }
        }
    }
    if fail {
        eprintln!(
            "\nThis test is restricted to VDJ cells having both chain types, uses at most \
            one cell\nper exact subclonotype, and uses up to 100 cells having the highest \
            UMI counts."
        );
        eprintln!(
            "\nThe data suggest a laboratory or informatic mixup.  If you believe \
            that this is not the case,\nyou can force enclone to run by adding \
            the argument ALLOW_INCONSISTENT to the command line.\n"
        );
        std::process::exit(1);
    }
    ctl.perf_stats(&tinc, "testing for inconsistency");

    // Filter out some foursie artifacts.

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
                            for j in 0..ex.clones.len() {
                                fate[ex.clones[j][0].dataset_index].insert(
                                    ex.clones[j][0].barcode.clone(),
                                    "failed FOURSIE_KILL filter".to_string(),
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    if ctl.clono_filt_opt.weak_foursies {
        erase_if(&mut exact_clonotypes, &to_delete);
    }
    ctl.perf_stats(&t, "filtering foursies");

    // Look for insertions (experimental).

    find_insertions(&ctl, &exact_clonotypes);

    // Build info about clonotypes.  Note that this edits the V reference sequence to perform
    // an indel in some cases.

    let tinfo = Instant::now();
    let mut info: Vec<CloneInfo> = build_info(&refdata, &ctl, &mut exact_clonotypes, &mut fate);
    ctl.perf_stats(&tinfo, "building info");

    // Derive consensus sequences for alternate alleles of V segments.  Then create donor
    // reference sequences for Loupe.

    let talt = Instant::now();
    let alt_refs : Vec<(usize,usize,DnaString)>  // {(donor, ref id, alt seq)}
        = find_alleles( &refdata, &ctl, &exact_clonotypes );
    ctl.perf_stats(&talt, "finding alt alleles");
    if ctl.gen_opt.dref_file.len() > 0 {
        let f = File::create(&ctl.gen_opt.dref_file);
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
                ">{}:{}:{}:{} (reference record id : donor name : allele number : gene name)\n{}",
                refdata.id[ref_id],
                ctl.origin_info.donor_id[donor],
                count + 1,
                refdata.name[ref_id],
                alt_seq.to_string()
            );
            count += 1;
        }
    }
    let tdonor = Instant::now();
    let drefs = make_donor_refs(&alt_refs, &refdata);
    ctl.perf_stats(&tdonor, "making donor refs");

    // Update reference sequences for V segments by substituting in alt alleles if better.

    sub_alts(&refdata, &ctl, &alt_refs, &mut info, &mut exact_clonotypes);

    // Form equivalence relation on exact subclonotypes.

    let mut join_info = Vec::<(usize, usize, bool, Vec<u8>)>::new();
    let mut eq: EquivRel = join_exacts(
        is_bcr,
        &refdata,
        &ctl,
        &exact_clonotypes,
        &info,
        &mut join_info,
    );

    // If NWEAK_ONESIES is not specified, disintegrate certain onesie clonotypes into single
    // cell clonotypes.  This requires editing of exact_clonotypes, info, eq and join_info.

    let mut disintegrated = Vec::<bool>::new();
    if ctl.clono_filt_opt.weak_onesies {
        let ncells_total = exact_clonotypes.iter().map(|x| x.ncells()).sum();
        let mut to_info = HashMap::<usize, usize>::new();
        let mut exacts2 = Vec::<ExactClonotype>::new();
        for i in 0..info.len() {
            to_info.insert(info[i].clonotype_index, i);
        }
        let mut to_exact_new = Vec::<Vec<usize>>::new();
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            let mut enew = Vec::<usize>::new();
            if ex.share.len() == 1
                && ex.ncells() > 1
                && ex.ncells() * 1000 < ncells_total
                && to_info.contains_key(&i)
                && eq.orbit_size(to_info[&i] as i32) == 1
            {
                for j in 0..ex.clones.len() {
                    enew.push(exacts2.len());
                    exacts2.push(ExactClonotype {
                        share: ex.share.clone(),
                        clones: vec![ex.clones[j].clone()],
                    });
                    disintegrated.push(true);
                }
            } else {
                enew.push(exacts2.len());
                exacts2.push(exact_clonotypes[i].clone());
                disintegrated.push(false);
            }
            to_exact_new.push(enew);
        }

        let mut join_info2 = Vec::new();
        for i in 0..join_info.len() {
            let (u1, u2) = (join_info[i].0, join_info[i].1);
            for v1 in to_exact_new[u1].iter() {
                for v2 in to_exact_new[u2].iter() {
                    let mut x = join_info[i].clone();
                    x.0 = *v1;
                    x.1 = *v2;
                    join_info2.push(x);
                }
            }
        }
        join_info = join_info2;
        exact_clonotypes = exacts2;
        let mut info2 = Vec::<CloneInfo>::new();
        let mut to_info2 = Vec::<Vec<usize>>::new();
        for i in 0..info.len() {
            let j = info[i].clonotype_index;
            let mut x = Vec::<usize>::new();
            for k in 0..to_exact_new[j].len() {
                info[i].clonotype_index = to_exact_new[j][k];
                info[i].clonotype_id = to_exact_new[j][k];
                x.push(info2.len());
                info2.push(info[i].clone());
            }
            to_info2.push(x);
        }
        info = info2;
        let mut reps = Vec::<i32>::new();
        eq.orbit_reps(&mut reps);
        let mut eq2 = EquivRel::new(info.len() as i32);
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            if o.len() > 1 {
                for j in 0..o.len() - 1 {
                    eq2.join(
                        to_info2[o[j] as usize][0] as i32,
                        to_info2[o[j + 1] as usize][0] as i32,
                    );
                }
            }
        }
        eq = eq2;
    }

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
    if is_tcr {
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            orbits.push(o);
        }
    } else {
        let mut umis = vec![Vec::<usize>::new(); ctl.origin_info.n()];
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
        let mut nu = vec![0; ctl.origin_info.n()];
        let mut umin = vec![0.0; ctl.origin_info.n()];
        for l in 0..ctl.origin_info.n() {
            umis[l].sort();
            nu[l] = umis[l].len();
            if ctl.gen_opt.baseline {
                println!(
                    "\n{} umi counts for dataset {} = {}",
                    nu[l],
                    l + 1,
                    ctl.origin_info.dataset_id[l]
                );
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
        // if ctl.clono_filt_opt.umi_filt || ctl.clono_filt_opt.umi_filt_mark {
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
                        if pass == 3 {
                            for i in 0..ex.clones.len() {
                                if to_delete[i] {
                                    fate[ex.clones[i][0].dataset_index].insert(
                                        ex.clones[i][0].barcode.clone(),
                                        "failed UMI filter".to_string(),
                                    );
                                }
                            }
                            if ctl.clono_filt_opt.umi_filt {
                                erase_if(&mut ex.clones, &to_delete);
                            }
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
            if !o.is_empty() {
                orbits.push(o.clone());
            }
        }
        // }
    }

    // Filter B cells based on UMI count ratios.  This assumes V..J identity to filter.

    if is_bcr {
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

                    if pass == 2 {
                        for i in 0..ex.clones.len() {
                            if to_delete[j][i] {
                                fate[ex.clones[i][0].dataset_index].insert(
                                    ex.clones[i][0].barcode.clone(),
                                    "failed UMI_RATIO filter".to_string(),
                                );
                            }
                        }
                        if ctl.clono_filt_opt.umi_ratio_filt {
                            erase_if(&mut ex.clones, &to_delete[j]);
                            if ex.ncells() == 0 {
                                to_deletex[j] = true;
                            }
                        }
                    }
                }
                if pass == 2 {
                    if ctl.clono_filt_opt.umi_ratio_filt {
                        erase_if(&mut o, &to_deletex);
                    }
                    if !o.is_empty() {
                        orbits2.push(o.clone());
                    }
                }
            }
        }
        orbits = orbits2;
    }

    // Remove cells that are not called cells by GEX or feature barcodes.

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
                if ctl.gen_opt.cellranger {
                    if gex_cells_specified[li] && !bin_member(&gex_cells[li], &bc) {
                        to_delete[k] = true;
                        fate[li].insert(bc.clone(), "failed GEX filter".to_string());
                    }
                } else if ctl.origin_info.gex_path[li].len() > 0 {
                    let gbc = &gex_info.gex_cell_barcodes[li];
                    if !bin_member(&gbc, &bc) {
                        fate[li].insert(bc.clone(), "failed GEX filter".to_string());
                        if !ctl.clono_filt_opt.ngex {
                            to_delete[k] = true;
                        }
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

    // Filter using constraints imposed by FCELL.

    if !ctl.clono_filt_opt.fcell.is_empty() {
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
                    let mut keep = true;
                    for x in ctl.clono_filt_opt.fcell.iter() {
                        let alt = &ctl.origin_info.alt_bc_fields[li];
                        let vars = x.iter_variable_identifiers().collect::<Vec<&str>>();
                        let mut vals = Vec::<String>::new();
                        for m in 0..vars.len() {
                            let mut val = String::new();
                            'uloop: for u in 0..alt.len() {
                                if alt[u].0 == vars[m] {
                                    if alt[u].1.contains_key(&bc.clone()) {
                                        val = alt[u].1[&bc.clone()].clone();
                                        break 'uloop;
                                    }
                                }
                            }
                            vals.push(val);
                        }
                        let mut c = HashMapContext::new();
                        for m in 0..vars.len() {
                            if vals[m].parse::<i64>().is_ok() {
                                c.set_value(
                                    vars[m].into(),
                                    evalexpr::Value::from(vals[m].force_i64()),
                                )
                                .unwrap();
                            } else if vals[m].parse::<f64>().is_ok() {
                                c.set_value(
                                    vars[m].into(),
                                    evalexpr::Value::from(vals[m].force_f64()),
                                )
                                .unwrap();
                            } else {
                                c.set_value(vars[m].into(), vals[m].clone().into()).unwrap();
                            }
                        }
                        let res = x.eval_with_context(&c);
                        let ok = res == Ok(evalexpr::Value::from(true));
                        if !ok {
                            keep = false;
                        }
                    }
                    if !keep {
                        to_delete[k] = true;
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

    // Delete exact subclonotypes that appear to represent doublets.

    if ctl.clono_filt_opt.doublet {
        // Define pure subclonotypes.  To do this we break each clonotype up by chain signature.
        // Note duplication of code with print_clonotypes.rs.  And this is doing some
        // superfluous compute.

        let mut results = Vec::<(usize, Vec<Vec<usize>>)>::new();
        for i in 0..orbits.len() {
            results.push((i, Vec::new()));
        }
        let mut pures = Vec::<Vec<usize>>::new();
        results.par_iter_mut().for_each(|res| {
            let i = res.0;
            let o = orbits[i].clone();
            let mut od = Vec::<(Vec<usize>, usize, i32)>::new();
            for id in o.iter() {
                let x: &CloneInfo = &info[*id as usize];
                od.push((x.origin.clone(), x.clonotype_id, *id));
            }
            od.sort();
            let mut exacts = Vec::<usize>::new();
            let mut cdr3s_len = Vec::<Vec<(String, usize)>>::new();
            let mut js = Vec::<usize>::new();
            let mut j = 0;
            while j < od.len() {
                let k = next_diff12_3(&od, j as i32) as usize;
                let mut z_len = Vec::<(String, usize)>::new();
                for l in j..k {
                    let x: &CloneInfo = &info[od[l].2 as usize];
                    for m in 0..x.cdr3_aa.len() {
                        let mut c = x.chain_types[m].clone();
                        if c.starts_with("TRB") {
                            c = c.replacen("TRB", "TRX", 1);
                        } else if c.starts_with("TRA") {
                            c = c.replacen("TRA", "TRY", 1);
                        }
                        z_len.push((format!("{}:{}", c, x.cdr3_aa[m]), x.lens[m]));
                    }
                }
                unique_sort(&mut z_len);
                cdr3s_len.push(z_len);
                js.push(j);
                exacts.push(od[j].1);
                j = k;
            }
            let mat = define_mat(&ctl, &exact_clonotypes, &cdr3s_len, &js, &od, &info);
            let nexacts = mat[0].len();
            let mut priority = Vec::<Vec<bool>>::new();
            for u in 0..nexacts {
                let mut typex = vec![false; mat.len()];
                for col in 0..mat.len() {
                    if mat[col][u].is_some() {
                        typex[col] = true;
                    }
                }
                priority.push(typex.clone());
            }
            sort_sync2(&mut priority, &mut exacts);
            let mut j = 0;
            while j < priority.len() {
                let k = next_diff(&priority, j);
                let mut p = Vec::<usize>::new();
                for l in j..k {
                    p.push(exacts[l] as usize);
                }
                res.1.push(p);
                j = k;
            }
        });
        for i in 0..results.len() {
            pures.append(&mut results[i].1);
        }

        // Define the number of cells in each pure subclonotype.

        let mut npure = vec![0; pures.len()];
        for j in 0..pures.len() {
            for id in pures[j].iter() {
                npure[j] += exact_clonotypes[*id].ncells();
            }
        }

        // Find the pairs of pure subclonotypes that share identical CDR3 sequences.

        let mut shares = Vec::<(usize, usize)>::new();
        {
            let mut content = Vec::<(String, usize)>::new();
            for j in 0..pures.len() {
                for id in pures[j].iter() {
                    let ex = &exact_clonotypes[*id];
                    for k in 0..ex.share.len() {
                        content.push((ex.share[k].cdr3_dna.clone(), j));
                    }
                }
            }
            unique_sort(&mut content);
            content.sort();
            let mut j = 0;
            while j < content.len() {
                let k = next_diff1_2(&content, j as i32) as usize;
                for l1 in j..k {
                    for l2 in j + 1..k {
                        shares.push((content[l1].1, content[l2].1));
                        shares.push((content[l2].1, content[l1].1));
                    }
                }
                j = k;
            }
            unique_sort(&mut shares);
        }

        // Find triples of pure subclonotypes in which the first two have no share, but both
        // of the first two share with the third.

        const MIN_MULT_DOUBLET: usize = 5;
        let mut trips = Vec::<(usize, usize, usize)>::new();
        {
            let mut us = Vec::<usize>::new();
            let mut vs = Vec::<Vec<usize>>::new();
            let mut j = 0;
            while j < shares.len() {
                let k = next_diff1_2(&shares, j as i32) as usize;
                let u = shares[j].0;
                us.push(u);
                let mut x = Vec::<usize>::new();
                for l in j..k {
                    let v = shares[l].1;
                    if MIN_MULT_DOUBLET * npure[u] <= npure[v] {
                        x.push(v);
                    }
                }
                vs.push(x);
                j = k;
            }
            let mut results = Vec::<(usize, Vec<(usize, usize, usize)>)>::new();
            for i in 0..us.len() {
                results.push((i, Vec::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let i = res.0;
                let u = us[i];
                let vs = &vs[i];
                for l1 in 0..vs.len() {
                    for l2 in l1 + 1..vs.len() {
                        let v1 = vs[l1];
                        let v2 = vs[l2];
                        if !bin_member(&shares, &(v1, v2)) {
                            res.1.push((v1, v2, u));
                        }
                    }
                }
            });
            for j in 0..results.len() {
                trips.append(&mut results[j].1.clone());
            }
        }

        // Delete some of the third members of the triples.

        // const MIN_DIST_DOUBLET: usize = 50;
        let mut to_delete = vec![false; exact_clonotypes.len()];
        for j in 0..trips.len() {
            let (v0, v1, v2) = (trips[j].2, trips[j].0, trips[j].1);
            {
                /*
                let ex1 = &exact_clonotypes[v1];
                let ex2 = &exact_clonotypes[v2];
                for i1 in 0..ex1.share.len() {
                    for i2 in 0..ex2.share.len() {
                        if ex1.share[i1].seq_del.len() == ex2.share[i2].seq_del.len() {
                            let mut diffs = 0;
                            for k in 0..ex1.share[i1].seq_del.len() {
                                if ex1.share[i1].seq_del[k] != ex2.share[i2].seq_del[k] {
                                    diffs += 1;
                                }
                            }
                            if diffs < MIN_DIST_DOUBLET {
                                ok = false;
                            }
                        }
                    }
                }
                */

                let verbose = false;
                if verbose {
                    println!("\n{}, {}, {}", v0, v1, v2);
                    println!("DELETING");
                    for (u, m) in pures[v0].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                    println!("USING");
                    for (u, m) in pures[v1].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                    println!("AND");
                    for (u, m) in pures[v2].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                }
                for m in pures[v0].iter() {
                    to_delete[*m] = true;
                }
            }
        }
        let mut orbits2 = Vec::<Vec<i32>>::new();
        for i in 0..orbits.len() {
            let mut o = orbits[i].clone();

            let mut del2 = vec![false; o.len()];
            for j in 0..o.len() {
                let id = info[o[j] as usize].clonotype_index;
                if to_delete[id] {
                    del2[j] = true;
                }
            }
            erase_if(&mut o, &del2);
            orbits2.push(o);
        }
        orbits = orbits2;
    }

    // Merge onesies where totally unambiguous.  Possibly inefficient and should optimize.

    if ctl.join_alg_opt.merge_onesies {
        // ctl.join_alg_opt.merge_onesies is always true
        let mut eqo = EquivRel::new(orbits.len() as i32);
        let mut to_orbit = vec![None; info.len()];
        for i in 0..orbits.len() {
            for j in 0..orbits[i].len() {
                to_orbit[orbits[i][j] as usize] = Some(i);
            }
        }
        let mut ncells_total = 0;
        for i in 0..exact_clonotypes.len() {
            ncells_total += exact_clonotypes[i].ncells();
        }
        let mut onesies = Vec::<usize>::new();
        for i in 0..info.len() {
            if to_orbit[i].is_some() && info[i].tigs.len() == 1 {
                if !ctl.clono_filt_opt.weak_onesies || !disintegrated[info[i].clonotype_index] {
                    onesies.push(i);
                }
            }
        }
        let mut alltigs2 = Vec::<(Vec<u8>, usize)>::new();
        for i in 0..info.len() {
            if to_orbit[i].is_some() && info[i].tigs.len() >= 2 {
                for j in 0..info[i].tigs.len() {
                    alltigs2.push((info[i].tigs[j].clone(), i));
                }
            }
        }
        alltigs2.sort();
        for x in onesies.iter() {
            let low = lower_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let high = upper_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let mut ms = Vec::<usize>::new();
            for m in low..high {
                if alltigs2[m as usize].0 == info[*x].tigs[0] {
                    ms.push(m as usize);
                }
            }
            let mut ok = ms.len() > 0;
            let mut exacts = Vec::<usize>::new();
            for j in 0..ms.len() {
                if eq.class_id(alltigs2[ms[j]].1 as i32) != eq.class_id(alltigs2[ms[0]].1 as i32) {
                    ok = false;
                }
                let mut o = Vec::<i32>::new();
                eq.orbit(alltigs2[ms[j]].1 as i32, &mut o);
                for z in o.iter() {
                    exacts.push(info[*z as usize].clonotype_index);
                }
            }
            unique_sort(&mut exacts);
            if ctl.join_alg_opt.merge_onesies_ctl {
                let ncells0 = exact_clonotypes[info[*x].clonotype_index].ncells();
                if ncells0 * 10000 < ncells_total {
                    ok = false;
                }
            }
            if ok {
                let orb1 = to_orbit[*x].unwrap();
                let orb2 = to_orbit[alltigs2[ms[0]].1].unwrap();

                // Test for donor mixing.

                if !ctl.clono_filt_opt.donor {
                    let mut donors = vec![Vec::<Option<usize>>::new(); 2];
                    let orbs = [&orb1, &orb2];
                    for (pass, orb) in orbs.iter().enumerate() {
                        for id in orbits[**orb as usize].iter() {
                            let ex = &exact_clonotypes[info[*id as usize].clonotype_id];
                            for i in 0..ex.clones.len() {
                                donors[pass].push(ex.clones[i][0].donor_index);
                            }
                        }
                        unique_sort(&mut donors[pass]);
                    }
                    if donors[0] != donors[1] {
                        continue;
                    }
                }

                // Make join.

                eqo.join(orb1 as i32, orb2 as i32);
            }
        }
        let mut orbits2 = Vec::<Vec<i32>>::new();
        let mut repsx = Vec::<i32>::new();
        eqo.orbit_reps(&mut repsx);
        for i in 0..repsx.len() {
            let mut ox = Vec::<i32>::new();
            eqo.orbit(repsx[i], &mut ox);
            let mut o = Vec::<i32>::new();
            for j in 0..ox.len() {
                o.append(&mut orbits[ox[j] as usize].clone());
            }
            orbits2.push(o);
        }
        orbits = orbits2;
    }

    // Check for disjoint orbits.  The test here should correspond exactly to the test in
    // define_mat, but it probably doesn't.

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
                                    let mut diffs = 0;
                                    let c1 = &ex1.share[m1].cdr3_aa.as_bytes();
                                    let c2 = &ex2.share[m2].cdr3_aa.as_bytes();
                                    for k in 0..c1.len() {
                                        if c1[k] != c2[k] {
                                            diffs += 1;
                                        }
                                    }
                                    if diffs <= MAX_CDR3_DIFFS_TO_JOIN {
                                        let mut tdiffs = 0;
                                        for k in 0..ex1.share[m1].seq_del.len() {
                                            if ex1.share[m1].seq_del[k] != ex2.share[m2].seq_del[k]
                                            {
                                                tdiffs += 1;
                                            }
                                        }
                                        if tdiffs <= ctl.heur.max_diffs {
                                            eqx.join(i1 as i32, i2 as i32);
                                            break 'cloop;
                                        }
                                    }
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

    // Mark VDJ noncells.

    if ctl.clono_filt_opt.non_cell_mark {
        for i in 0..exact_clonotypes.len() {
            let ex = &mut exact_clonotypes[i];
            for j in 0..ex.clones.len() {
                let di = ex.clones[j][0].dataset_index;
                if !bin_member(&vdj_cells[di], &ex.clones[j][0].barcode) {
                    ex.clones[j][0].marked = true;
                }
            }
        }
    }

    ctl.perf_stats(&tumi, "umi filtering and such");

    // Load the GEX and FB data.

    let tdi = Instant::now();
    let mut d_readers = Vec::<Option<hdf5::Reader>>::new();
    let mut ind_readers = Vec::<Option<hdf5::Reader>>::new();
    for li in 0..ctl.origin_info.n() {
        if ctl.origin_info.gex_path[li].len() > 0 && !gex_info.gex_matrices[li].initialized() {
            let x = gex_info.h5_data[li].as_ref();
            if x.is_none() {
                // THIS FAILS SPORADICALLY, OBSERVED MULTIPLE TIMES,
                // CAUSING PUSH TO D_READERS BELOW TO FAIL.
                eprintln!("\nWeird, gex_info.h5_data[li].as_ref() is None.");
                eprintln!("Path = {}.", ctl.origin_info.gex_path[li]);
                let current = env::current_dir().unwrap();
                println!(
                    "The current working directory is {}",
                    current.canonicalize().unwrap().display()
                );
                if path_exists(&ctl.origin_info.gex_path[li]) {
                    eprintln!(
                        "The directory that is supposed to contain \
                        raw_feature_bc_matrix.h5 exists."
                    );
                    let list = dir_list(&ctl.origin_info.gex_path[li]);
                    eprintln!(
                        "This directory is {} and its contents are:",
                        ctl.origin_info.gex_path[li]
                    );
                    for i in 0..list.len() {
                        eprintln!("{}.  {}", i + 1, list[i]);
                    }
                    let h5_path =
                        format!("{}/raw_feature_bc_matrix.h5", ctl.origin_info.gex_path[li]);
                    eprintln!("H5 path = {}.", h5_path);
                    if !path_exists(&h5_path) {
                        eprintln!("H5 path {} does not exist.", h5_path);
                        eprintln!("Retrying a few times to see if it appears.");
                        for _ in 0..5 {
                            eprintln!("Sleeping for 0.1 seconds.");
                            thread::sleep(time::Duration::from_millis(100));
                            if !path_exists(&h5_path) {
                                eprintln!("Now h5 path does not exist.");
                            } else {
                                eprintln!("Now h5 path exists.");
                                break;
                            }
                        }
                        eprintln!("Aborting.\n");
                        std::process::exit(1);
                    } else {
                        eprintln!("h5 path exists.");
                    }
                } else {
                    eprintln!("Path exists.");
                }
                eprintln!("");
            }
            d_readers.push(Some(x.unwrap().as_reader()));
            ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
        } else {
            d_readers.push(None);
            ind_readers.push(None);
        }
    }
    let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
    for li in 0..ctl.origin_info.n() {
        h5_data.push((li, Vec::new(), Vec::new()));
    }
    h5_data.par_iter_mut().for_each(|res| {
        let li = res.0;
        if ctl.origin_info.gex_path[li].len() > 0
            && !gex_info.gex_matrices[li].initialized()
            && ctl.gen_opt.h5_pre
        {
            res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
            res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
        }
    });
    ctl.perf_stats(&tdi, "setting up readers");

    // Find and print clonotypes.  (But we don't actually print them here.)

    let torb = Instant::now();
    let mut pics = Vec::<String>::new();
    let mut exacts = Vec::<Vec<usize>>::new(); // ugly reuse of name
    let mut in_center = Vec::<bool>::new();
    let mut rsi = Vec::<ColInfo>::new(); // ditto
    let mut out_datas = Vec::<Vec<HashMap<String, String>>>::new();
    let mut tests = Vec::<usize>::new();
    let mut controls = Vec::<usize>::new();
    print_clonotypes(
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &orbits,
        &gex_info,
        &vdj_cells,
        &d_readers,
        &ind_readers,
        &h5_data,
        &mut pics,
        &mut exacts,
        &mut in_center,
        &mut rsi,
        &mut out_datas,
        &mut tests,
        &mut controls,
        &mut fate,
    );

    // Process the SUBSET_JSON option.

    if ctl.gen_opt.subset_json.len() > 0 {
        let mut barcodes = Vec::<String>::new();
        for l in 0..exacts.len() {
            for u in 0..exacts[l].len() {
                let ex = &exact_clonotypes[exacts[l][u]];
                for j in 0..ex.clones.len() {
                    barcodes.push(ex.clones[j][0].barcode.clone());
                }
            }
        }
        unique_sort(&mut barcodes);
        let mut g = open_for_write_new![&ctl.gen_opt.subset_json];
        fwriteln!(g, "[");
        for li in 0..ctl.origin_info.dataset_path.len() {
            let json = format!("{}/{}", ctl.origin_info.dataset_path[li], ann);
            let mut jsonx = json.clone();
            if !path_exists(&json) {
                jsonx = format!("{}.lz4", json);
            }
            let mut xs = Vec::<Vec<u8>>::new();
            let mut f = BufReader::new(open_maybe_compressed(&jsonx));
            loop {
                match read_vector_entry_from_json(&mut f) {
                    None => break,
                    Some(x) => {
                        let v: Value = serde_json::from_str(strme(&x)).unwrap();
                        let barcode = &v["barcode"].to_string().between("\"", "\"").to_string();
                        if bin_member(&barcodes, &barcode) {
                            xs.push(x);
                        }
                    }
                }
            }
            for j in 0..xs.len() {
                fwrite!(g, "{}", strme(&xs[j]));
                if j < xs.len() - 1 {
                    fwrite!(g, ",");
                }
                fwriteln!(g, "");
            }
        }
        fwriteln!(g, "]");
    }
    ctl.perf_stats(&torb, "making orbits");

    // Tail code.

    let ttail = Instant::now();
    tail_code(
        &tall,
        &refdata,
        &pics,
        &exacts,
        &in_center,
        &rsi,
        &exact_clonotypes,
        &ctl,
        &mut out_datas,
        &join_info,
        &gex_info,
        &vdj_cells,
        &fate,
        &tests,
        &controls,
        &h5_data,
        &d_readers,
        &ind_readers,
        &drefs,
    );
    ctl.perf_stats(&ttail, "in tail code");

    // Report computational performance.

    let delta;
    unsafe {
        delta = elapsed(&tall) - WALLCLOCK;
    }
    ctl.perf_stats(&tall, "total");
    if ctl.comp {
        println!("used {:.2} seconds unaccounted for", delta);
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

    if !(ctl.gen_opt.noprint && ctl.parseable_opt.pout == "stdout") {
        println!("");
    }
    // It's not totally clear that the exit below actually saves time.  Would need more testing.
    if !ctl.gen_opt.cellranger {
        std::process::exit(0);
    }
}
