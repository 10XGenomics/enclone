// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use self::refx::*;
use crate::determine_ref::*;
use crate::doublets::*;
use crate::filter_umi::*;
use crate::flag_defective::*;
use crate::inconsistent::*;
use crate::merge_onesies::*;
use crate::setup::*;
use crate::split_orbits::*;
use amino::*;
use debruijn::dna_string::DnaString;
use enclone::allele::*;
use enclone::explore::*;
use enclone::graph_filter::*;
use enclone::info::*;
use enclone::innate::*;
use enclone::join::*;
use enclone::load_gex::*;
use enclone::misc1::*;
use enclone::misc2::*;
use enclone::misc3::*;
use enclone::proc_args_check::*;
use enclone::read_json::*;
use enclone::secret::*;
use enclone_core::defs::*;
use enclone_core::vdj_features::*;
use enclone_print::loupe::*;
use enclone_print::print_clonotypes::*;
use enclone_tail::tail::tail_code;
use equiv::EquivRel;
use evalexpr::*;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use regex::Regex;
use serde_json::Value;
use stats_utils::*;
use std::{
    collections::HashMap,
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    process::Command,
    thread, time,
    time::Instant,
};
use stirling_numbers::*;
use string_utils::*;
use tables::*;
use vdj_ann::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone(args: &Vec<String>) {
    // Set up stuff, read args, etc.

    let mut ctl = EncloneControl::default();
    for i in 0..args.len() {
        let arg = &args[i];
        if arg == "PROFILE" {
            ctl.gen_opt.profile = true;
        }
    }
    if ctl.gen_opt.profile {
        let blacklist = [
            "alloc",
            "build",
            "core",
            "core-arch",
            "crossbeam-deque",
            "crossbeam-epoch",
            "debruijn",
            "float-ord",
            "hashbrown",
            "hdf5-rust",
            "hdf5-types",
            "lock_api",
            "lz4",
            "ndarray",
            "parking_lot",
            "parking_lot_core",
            "rayon",
            "rayon-core",
            "regex",
            "regex-syntax",
            "rust-bio",
            "serde",
            "serde_json",
            "std",
            "superslice",
            "unknown",
        ];
        let mut b = Vec::<String>::new();
        for x in blacklist.iter() {
            b.push(x.to_string());
        }
        start_profiling(&b);
    }
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
    ctl.perf_stats(&tall, "before setup");
    setup(&mut ctl, &args);

    // Read external data.

    if ctl.gen_opt.ext.len() > 0 {
        let f = open_userfile_for_read(&ctl.gen_opt.ext);
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
    // has to occur after loading GEX data.  This could also occur after loading only the feature
    // list, which would be better.

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

    // Start of code to determine the reference sequence that is to be used.

    let tr = Instant::now();
    let mut refx = String::new();
    let ann;
    if !ctl.gen_opt.cellranger {
        ann = "all_contig_annotations.json";
    } else {
        ann = "contig_annotations.json";
    }
    determine_ref(&mut ctl, &mut refx);
    if refx.len() == 0 && ctl.origin_info.n() == 0 {
        eprintln!("\nNo data and no TCR or BCR data have been specified.\n");
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

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
    let mut broken = Vec::<bool>::new();
    flag_defective(&ctl, &refdata, &mut log, &mut broken);

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

    // Test for no data.

    if ctl.origin_info.n() == 0 {
        eprintln!("\nNo TCR or BCR data have been specified.\n");
        std::process::exit(1);
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

    // Test for consistency between VDJ cells and GEX cells.

    test_vdj_gex_inconsistent(&ctl, &tig_bc, &exact_clonotypes, &vdj_cells, &gex_info);

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

    // Compute to_bc, which maps (dataset_index, clonotype_id) to {barcodes}.
    // This is intended as a replacement for some old code below.

    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for i in 0..exact_clonotypes.len() {
        for j in 0..exact_clonotypes[i].clones.len() {
            let x = &exact_clonotypes[i].clones[j][0];
            if !to_bc.contains_key(&(x.dataset_index, i)) {
                to_bc.insert((x.dataset_index, i), vec![x.barcode.clone()]);
            } else {
                to_bc
                    .get_mut(&(x.dataset_index, i))
                    .unwrap()
                    .push(x.barcode.clone());
            }
        }
    }

    // Make stirling ratio table.  Not sure that fixing the size of this is safe.

    let sr = stirling2_ratio_table::<f64>(3000);

    // Form equivalence relation on exact subclonotypes.  We also keep the raw joins, consisting
    // of pairs of info indices, that were originally joined.

    let mut join_info = Vec::<(usize, usize, bool, Vec<u8>)>::new();
    let mut raw_joins = Vec::<(i32, i32)>::new();
    let mut eq: EquivRel = join_exacts(
        is_bcr,
        &to_bc,
        &refdata,
        &ctl,
        &exact_clonotypes,
        &info,
        &mut join_info,
        &mut raw_joins,
        &sr,
    );

    // If NWEAK_ONESIES is not specified, disintegrate certain onesie clonotypes into single cell
    // clonotypes.  This requires editing of exact_clonotypes, info, eq, join_info and raw_joins.

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
                let mut origins = Vec::<usize>::new();
                let ex = &exact_clonotypes[info[i].clonotype_index];
                for i in 0..ex.clones.len() {
                    origins.push(ex.clones[i][0].dataset_index);
                }
                unique_sort(&mut origins);
                info[i].origin = origins;
                x.push(info2.len());
                info2.push(info[i].clone());
            }
            to_info2.push(x);
        }
        info = info2;
        let mut raw_joins2 = Vec::<(i32, i32)>::new();
        for i in 0..raw_joins.len() {
            let (j1, j2) = (
                &to_info2[raw_joins[i].0 as usize],
                &to_info2[raw_joins[i].1 as usize],
            );
            raw_joins2.push((j1[0] as i32, j2[0] as i32));
        }
        raw_joins = raw_joins2;
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

    // Update to_bc.

    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for i in 0..exact_clonotypes.len() {
        for j in 0..exact_clonotypes[i].clones.len() {
            let x = &exact_clonotypes[i].clones[j][0];
            if !to_bc.contains_key(&(x.dataset_index, i)) {
                to_bc.insert((x.dataset_index, i), vec![x.barcode.clone()]);
            } else {
                to_bc
                    .get_mut(&(x.dataset_index, i))
                    .unwrap()
                    .push(x.barcode.clone());
            }
        }
    }

    // Restructure raw joins.

    raw_joins.sort();
    let mut raw_joins2 = vec![Vec::<usize>::new(); info.len()];
    for i in 0..raw_joins.len() {
        raw_joins2[raw_joins[i].0 as usize].push(raw_joins[i].1 as usize);
        raw_joins2[raw_joins[i].1 as usize].push(raw_joins[i].0 as usize);
    }
    let raw_joins = raw_joins2;

    // Lock info.

    let info = &info;

    // Lookup for heavy chain reuse (special purpose experimental option).

    lookup_heavy_chain_reuse(&ctl, &exact_clonotypes, &info, &eq);

    // Filter B cells based on UMI counts.

    let tumi = Instant::now();
    let mut orbits = Vec::<Vec<i32>>::new();
    filter_umi(
        &eq,
        &mut orbits,
        &ctl,
        &mut exact_clonotypes,
        &info,
        &mut fate,
    );

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

    delete_doublets(
        &mut orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );

    // Merge onesies where totally unambiguous, then check for disjoint orbits.

    merge_onesies(
        &mut orbits,
        &ctl,
        &exact_clonotypes,
        &info,
        &eq,
        &disintegrated,
    );
    split_orbits(
        &mut orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );

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
        is_bcr,
        &to_bc,
        &sr,
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &orbits,
        &raw_joins,
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

    // Report profiling.

    if ctl.gen_opt.profile {
        let t = Instant::now();
        stop_profiling();
        ctl.perf_stats(&t, "summarizing profiling");
    }

    // Report computational performance.

    let delta;
    unsafe {
        delta = elapsed(&tall) - WALLCLOCK;
    }
    ctl.perf_stats(&tall, "total");
    if ctl.comp {
        println!("used {:.2} seconds unaccounted for", delta);
        println!("peak mem usage = {:.1} MB", peak_mem_usage_gb() * 1000.0);
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
