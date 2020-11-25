// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::proc_args2::*;
use crate::proc_args3::*;
use crate::proc_args_check::*;
use enclone_core::defs::*;
use enclone_core::testlist::*;
use evalexpr::*;
use io_utils::*;
use itertools::Itertools;
use regex::Regex;
use std::fs::{remove_file, File};
use std::io::{BufRead, BufReader};
use std::{env, process::Command, time::Instant};
use string_utils::*;
use tilde_expand::*;
use vector_utils::*;

// Process arguments.

pub fn proc_args(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    // Knobs.

    let targs = Instant::now();
    let heur = ClonotypeHeuristics {
        max_diffs: 50,
        max_degradation: 3,
        ref_v_trim: 15,
        ref_j_trim: 15,
    };
    ctl.heur = heur;

    // Form the combined set of command-line arguments and "command-line" arguments
    // implied by environment variables.

    let mut args = args.clone();
    let mut args2 = Vec::<String>::new();
    args2.push(args[0].clone());
    for (key, value) in env::vars() {
        if key.starts_with("ENCLONE_") {
            args2.push(format!("{}={}", key.after("ENCLONE_"), value));
        }
    }
    for i in 1..args.len() {
        args2.push(args[i].clone());
    }
    args = args2;

    // Test for internal run.

    for (key, value) in env::vars() {
        if (key == "HOST" || key == "HOSTNAME") && value.ends_with(".fuzzplex.com") {
            ctl.gen_opt.internal_run = true;
        }
    }
    for i in 1..args.len() {
        if args[i] == "FORCE_EXTERNAL".to_string() {
            ctl.gen_opt.internal_run = false;
        }
    }
    if ctl.gen_opt.internal_run {
        ctl.gen_opt.pre = vec![
            format!("/mnt/assembly/vdj/current{}", TEST_FILES_VERSION),
            format!("enclone/test/inputs"),
            format!("enclone_main"),
        ];
    } else if !ctl.gen_opt.cellranger {
        let home = dirs::home_dir().unwrap().to_str().unwrap().to_string();
        ctl.gen_opt.pre = vec![
            format!("{}/enclone/datasets", home),
            format!("{}/enclone/datasets2", home),
        ];
    }

    // Process special option SPLIT_COMMAND.

    let mut split = false;
    for i in 1..args.len() {
        if args[i] == "SPLIT_BY_COMMAND" {
            split = true;
        }
    }
    if split {
        let (mut bcr, mut gex) = (Vec::<&str>::new(), Vec::<&str>::new());
        let mut args2 = Vec::<String>::new();
        for i in 1..args.len() {
            if args[i] == "SPLIT_BY_COMMAND" {
            } else if args[i].starts_with("BCR=") {
                bcr = args[i].after("BCR=").split(',').collect::<Vec<&str>>();
            } else if args[i].starts_with("GEX=") {
                gex = args[i].after("GEX=").split(',').collect::<Vec<&str>>();
            } else {
                args2.push(args[i].to_string());
            }
        }
        for i in 0..bcr.len() {
            let mut args = args2.clone();
            args.push(format!("BCR={}", bcr[i]));
            args.push(format!("GEX={}", gex[i]));
            println!("\nenclone {}\n", args.iter().format(" "));
            let o = Command::new("enclone")
                .args(&args)
                .output()
                .expect("failed to execute enclone");
            print!("{}{}", strme(&o.stdout), strme(&o.stderr));
            if o.status.code() != Some(0) {
                println!("FAILED!\n");
                std::process::exit(1);
            }
        }
        std::process::exit(0);
    }

    // Set up general options.

    ctl.gen_opt.h5_pre = true;
    ctl.gen_opt.min_cells_exact = 1;
    ctl.gen_opt.min_chains_exact = 1;
    ctl.gen_opt.exact = None;
    for i in 1..args.len() {
        if args[i].starts_with("PRE=") {
            let pre = args[i].after("PRE=").split(',').collect::<Vec<&str>>();
            ctl.gen_opt.pre.clear();
            for x in pre.iter() {
                ctl.gen_opt.pre.push(x.to_string());
            }
        }
    }
    ctl.gen_opt.full_counts = true;
    ctl.gen_opt.color = "codon".to_string();
    ctl.silent = true;
    ctl.gen_opt.peer_group_dist = "MFL".to_string();
    ctl.gen_opt.color_by_rarity_pc = -1.0;

    // Set up clonotyping control parameters.

    ctl.clono_filt_opt.ncells_low = 1;
    ctl.clono_filt_opt.ncells_high = 1_000_000_000;
    ctl.clono_filt_opt.min_umi = 0;
    ctl.clono_filt_opt.max_chains = 1000000;
    ctl.clono_filt_opt.qual_filter = true;
    ctl.clono_filt_opt.weak_chains = true;
    ctl.clono_filt_opt.weak_onesies = true;
    ctl.clono_filt_opt.weak_foursies = true;
    ctl.clono_filt_opt.bc_dup = true;
    ctl.clono_filt_opt.max_datasets = 1000000000;
    ctl.clono_filt_opt.umi_filt = true;
    ctl.clono_filt_opt.umi_ratio_filt = true;

    ctl.clono_print_opt.amino = vec![
        "cdr3".to_string(),
        "var".to_string(),
        "share".to_string(),
        "donor".to_string(),
    ];
    ctl.clono_print_opt.cvars = vec!["u".to_string(), "const".to_string(), "notes".to_string()];
    ctl.clono_print_opt.lvars = vec!["datasets".to_string(), "n".to_string()];

    ctl.clono_group_opt.min_group = 1;

    ctl.allele_alg_opt.min_mult = 4;
    ctl.allele_alg_opt.min_alt = 4;

    ctl.join_alg_opt.max_score = 1_000_000.0;
    ctl.join_alg_opt.merge_onesies = true; // should just kill this as an option
    ctl.join_alg_opt.merge_onesies_ctl = true;
    ctl.join_alg_opt.max_cdr3_diffs = 10;

    ctl.join_print_opt.pfreq = 1_000_000_000;
    ctl.join_print_opt.quiet = true;

    ctl.parseable_opt.pchains = 4;

    // Pretest for consistency amongst TCR, BCR, GEX and META.  Also preparse GEX.

    let (mut have_tcr, mut have_bcr) = (false, false);
    let mut have_gex = false;
    let mut have_meta = false;
    let mut gex = String::new();
    let mut bc = String::new();
    let mut metas = Vec::<String>::new();
    let mut xcrs = Vec::<String>::new();
    for i in 1..args.len() {
        if args[i].starts_with("BI=") {
            have_bcr = true;
            have_gex = true;
        } else if args[i].starts_with("TCR=") {
            have_tcr = true;
        } else if args[i].starts_with("BCR=") {
            have_bcr = true;
        } else if args[i].starts_with("GEX=") {
            have_gex = true;
        } else if args[i].starts_with("META=") {
            have_meta = true;
        }
        if args[i].starts_with("GEX=") {
            gex = args[i].after("GEX=").to_string();
        }
        if args[i].starts_with("BC=") {
            bc = args[i].after("BC=").to_string();
        }
        if is_simple_arg(&args[i], "MARK_STATS") {
            ctl.gen_opt.mark_stats = true;
        }
        if is_simple_arg(&args[i], "MARK_STATS2") {
            ctl.gen_opt.mark_stats2 = true;
        }
        if is_simple_arg(&args[i], "MARKED_B") {
            ctl.clono_filt_opt.marked_b = true;
        }
    }
    if have_meta && (have_tcr || have_bcr || have_gex || bc.len() > 0) {
        eprintln!("\nIf META is specified, then none of TCR, BCR, GEX or BC can be specified.\n");
        std::process::exit(1);
    }
    if have_tcr && have_bcr {
        eprintln!("\nKindly please do not specify both TCR and BCR.\n");
        std::process::exit(1);
    }
    let mut using_plot = false;

    // Preprocess BI argument.

    if ctl.gen_opt.internal_run {
        for i in 1..args.len() {
            if args[i].starts_with("BI=") {
                let x = args[i].after("BI=").split(',').collect::<Vec<&str>>();
                let mut y = Vec::<String>::new();
                for j in 0..x.len() {
                    if x[j].contains('-') {
                        let (start, stop) = (x[j].before("-"), x[j].after("-"));
                        if !start.parse::<usize>().is_ok()
                            || !stop.parse::<usize>().is_ok()
                            || start.force_usize() > stop.force_usize()
                        {
                            eprintln!("\nIllegal range in BI argument.\n");
                            std::process::exit(1);
                        }
                        let (start, stop) = (start.force_usize(), stop.force_usize());
                        for j in start..=stop {
                            y.push(format!("{}", j));
                        }
                    } else {
                        y.push(x[j].to_string());
                    }
                }
                let mut args2 = Vec::<String>::new();
                for j in 0..i {
                    args2.push(args[j].clone());
                }
                let f = include_str!["enclone.testdata.bcr.gex"];
                let (mut bcrv, mut gexv) = (Vec::<String>::new(), Vec::<String>::new());
                for n in y.iter() {
                    if *n != "m1" {
                        if !n.parse::<usize>().is_ok()
                            || n.force_usize() < 1
                            || n.force_usize() > 12
                        {
                            eprintln!(
                                "\nBI only works for values n with if 1 <= n <= 12, or n = m1.\n"
                            );
                            std::process::exit(1);
                        }
                    } else if y.len() > 1 {
                        eprintln!("\nFor BI, if you specify m1, you can only specify m1.\n");
                        std::process::exit(1);
                    }
                    let mut found = false;
                    for s in f.lines() {
                        if s == format!("DONOR={}", n) {
                            found = true;
                        } else if found && s.starts_with("DONOR=") {
                            break;
                        }
                        if found {
                            if s.starts_with("BCR=") {
                                bcrv.push(s.after("BCR=").to_string());
                            }
                            if s.starts_with("GEX=") {
                                gexv.push(s.after("GEX=").to_string());
                            }
                            if s == "SPECIES=mouse" {
                                args2.push("MOUSE".to_string());
                            }
                        }
                    }
                }
                args2.push(format!("BCR={}", bcrv.iter().format(";")));
                args2.push(format!("GEX={}", gexv.iter().format(";")));
                gex = format!("{}", gexv.iter().format(";"));
                for j in i + 1..args.len() {
                    args2.push(args[j].clone());
                }
                args = args2;
                break;
            }
        }
    }

    // Preprocess NALL and NALL_GEX.

    for i in 1..args.len() {
        if args[i] == "NALL".to_string() || args[i] == "NALL_CELL" || args[i] == "NALL_GEX" {
            let f = [
                "NCELL",
                "NGEX",
                "NCROSS",
                "NUMI",
                "NUMI_RATIO",
                "NGRAPH_FILTER",
                "NQUAL",
                "NWEAK_CHAINS",
                "NWEAK_ONESIES",
                "NFOURSIE_KILL",
                "NWHITEF",
                "NBC_DUP",
                "MIX_DONORS",
                "NIMPROPER",
            ];
            for j in 0..f.len() {
                if f[j] == "NCELL" {
                    if args[i] != "NALL_CELL" {
                        args.push(f[j].to_string());
                    }
                } else if f[j] == "NGEX" {
                    if args[i] != "NALL_GEX" {
                        args.push(f[j].to_string());
                    }
                } else {
                    args.push(f[j].to_string());
                }
            }
            break;
        }
    }

    // Define arguments that set something to true.

    let mut set_true = vec![
        ("ACCEPT_BROKEN", &mut ctl.gen_opt.accept_broken),
        ("ACCEPT_INCONSISTENT", &mut ctl.gen_opt.accept_inconsistent),
        ("ACCEPT_REUSE", &mut ctl.gen_opt.accept_reuse),
        ("ALLOW_INCONSISTENT", &mut ctl.gen_opt.allow_inconsistent),
        ("ANN", &mut ctl.join_print_opt.ann),
        ("ANN0", &mut ctl.join_print_opt.ann0),
        ("BARCODES", &mut ctl.clono_print_opt.barcodes),
        ("BASELINE", &mut ctl.gen_opt.baseline),
        ("BCJOIN", &mut ctl.join_alg_opt.bcjoin),
        ("BUILT_IN", &mut ctl.gen_opt.built_in),
        ("CDIFF", &mut ctl.clono_filt_opt.cdiff),
        ("CHAIN_BRIEF", &mut ctl.clono_print_opt.chain_brief),
        ("COMPLETE", &mut ctl.gen_opt.complete),
        ("CON", &mut ctl.allele_print_opt.con),
        ("CON_CON", &mut ctl.gen_opt.con_con),
        ("CON_TRACE", &mut ctl.allele_print_opt.con_trace),
        ("CURRENT_REF", &mut ctl.gen_opt.current_ref),
        ("DEBUG_TABLE_PRINTING", &mut ctl.debug_table_printing),
        ("DEL", &mut ctl.clono_filt_opt.del),
        ("DESCRIP", &mut ctl.gen_opt.descrip),
        ("EASY", &mut ctl.join_alg_opt.easy),
        ("ECHO", &mut ctl.gen_opt.echo),
        ("EXP", &mut ctl.gen_opt.exp),
        ("FORCE", &mut ctl.force),
        ("FULL_SEQC", &mut ctl.clono_print_opt.full_seqc),
        ("GRAPH", &mut ctl.gen_opt.graph),
        ("GROUP_HEAVY_CDR3", &mut ctl.clono_group_opt.heavy_cdr3_aa),
        ("GROUP_VJ_REFNAME", &mut ctl.clono_group_opt.vj_refname),
        (
            "GROUP_VJ_REFNAME_STRONG",
            &mut ctl.clono_group_opt.vj_refname_strong,
        ),
        ("HAVE_ONESIE", &mut ctl.clono_filt_opt.have_onesie),
        ("HEAVY_CHAIN_REUSE", &mut ctl.gen_opt.heavy_chain_reuse),
        ("IMGT", &mut ctl.gen_opt.imgt),
        ("IMGT_FIX", &mut ctl.gen_opt.imgt_fix),
        ("INDELS", &mut ctl.gen_opt.indels),
        ("INKT", &mut ctl.clono_filt_opt.inkt),
        ("INSERTIONS", &mut ctl.gen_opt.insertions),
        ("JC1", &mut ctl.gen_opt.jc1),
        ("MAIT", &mut ctl.clono_filt_opt.mait),
        ("MARKED", &mut ctl.clono_filt_opt.marked),
        ("MEAN", &mut ctl.clono_print_opt.mean),
        ("MIX_DONORS", &mut ctl.clono_filt_opt.donor),
        ("MOUSE", &mut ctl.gen_opt.mouse),
        ("NCELL", &mut ctl.gen_opt.ncell),
        ("NCROSS", &mut ctl.clono_filt_opt.ncross),
        ("NEWICK", &mut ctl.gen_opt.newick),
        ("NGEX", &mut ctl.clono_filt_opt.ngex),
        ("NGRAPH_FILTER", &mut ctl.gen_opt.ngraph_filter),
        ("NGROUP", &mut ctl.gen_opt.ngroup),
        ("NIMPROPER", &mut ctl.merge_all_impropers),
        ("NON_CELL_MARK", &mut ctl.clono_filt_opt.non_cell_mark),
        ("NOPRINT", &mut ctl.gen_opt.noprint),
        ("NOTE_SIMPLE", &mut ctl.clono_print_opt.note_simple),
        ("NPLAIN", &mut ctl.pretty),
        ("NWHITEF", &mut ctl.gen_opt.nwhitef),
        ("NWARN", &mut ctl.gen_opt.nwarn),
        ("PCELL", &mut ctl.parseable_opt.pbarcode),
        ("PG_READABLE", &mut ctl.gen_opt.peer_group_readable),
        ("PER_CELL", &mut ctl.clono_print_opt.bu),
        ("PROTECT_BADS", &mut ctl.clono_filt_opt.protect_bads),
        ("RE", &mut ctl.gen_opt.reannotate),
        ("REPROD", &mut ctl.gen_opt.reprod),
        ("REQUIRE_UNBROKEN_OK", &mut ctl.gen_opt.require_unbroken_ok),
        ("REUSE", &mut ctl.gen_opt.reuse),
        ("SEQC", &mut ctl.clono_print_opt.seqc),
        ("SHOW_BC", &mut ctl.join_print_opt.show_bc),
        ("STABLE_DOC", &mut ctl.gen_opt.stable_doc),
        ("SUM", &mut ctl.clono_print_opt.sum),
        ("SUMMARY", &mut ctl.gen_opt.summary),
        ("SUMMARY_CLEAN", &mut ctl.gen_opt.summary_clean),
        ("SUMMARY_CSV", &mut ctl.gen_opt.summary_csv),
        ("TOY", &mut ctl.toy),
        ("UMI_FILT_MARK", &mut ctl.clono_filt_opt.umi_filt_mark),
        (
            "UMI_RATIO_FILT_MARK",
            &mut ctl.clono_filt_opt.umi_ratio_filt_mark,
        ),
        ("UTR_CON", &mut ctl.gen_opt.utr_con),
        ("VDUP", &mut ctl.clono_filt_opt.vdup),
        ("WEAK", &mut ctl.gen_opt.weak),
        ("WHITEF", &mut ctl.clono_filt_opt.whitef),
    ];

    // Define arguments that set something to false.

    let mut set_false = vec![
        ("H5_SLICE", &mut ctl.gen_opt.h5_pre),
        ("NBC_DUP", &mut ctl.clono_filt_opt.bc_dup),
        ("NFOURSIE_KILL", &mut ctl.clono_filt_opt.weak_foursies),
        ("NMERGE_ONESIES", &mut ctl.join_alg_opt.merge_onesies_ctl),
        ("NQUAL", &mut ctl.clono_filt_opt.qual_filter),
        ("NSILENT", &mut ctl.silent),
        ("NUMI", &mut ctl.clono_filt_opt.umi_filt),
        ("NUMI_RATIO", &mut ctl.clono_filt_opt.umi_ratio_filt),
        ("NWEAK_CHAINS", &mut ctl.clono_filt_opt.weak_chains),
        ("NWEAK_ONESIES", &mut ctl.clono_filt_opt.weak_onesies),
        ("PRINT_FAILED_JOINS", &mut ctl.join_print_opt.quiet),
    ];

    // Define arguments that set something to a usize.

    let set_usize = [
        ("CHAINS_EXACT", &mut ctl.gen_opt.chains_exact),
        ("MAX_CDR3_DIFFS", &mut ctl.join_alg_opt.max_cdr3_diffs),
        ("MAX_DATASETS", &mut ctl.clono_filt_opt.max_datasets),
        ("MAX_DEGRADATION", &mut ctl.heur.max_degradation),
        ("MAX_DIFFS", &mut ctl.heur.max_diffs),
        ("MIN_ALT", &mut ctl.allele_alg_opt.min_alt),
        ("MIN_CELLS_EXACT", &mut ctl.gen_opt.min_cells_exact),
        ("MIN_CHAINS_EXACT", &mut ctl.gen_opt.min_chains_exact),
        (
            "MIN_DATASET_RATIO",
            &mut ctl.clono_filt_opt.min_dataset_ratio,
        ),
        ("MIN_DATASETS", &mut ctl.clono_filt_opt.min_datasets),
        ("MIN_EXACTS", &mut ctl.clono_filt_opt.min_exacts),
        ("MIN_GROUP", &mut ctl.clono_group_opt.min_group),
        ("MIN_MULT", &mut ctl.allele_alg_opt.min_mult),
        ("MIN_UMI", &mut ctl.clono_filt_opt.min_umi),
        ("PCHAINS", &mut ctl.parseable_opt.pchains),
        ("PFREQ", &mut ctl.join_print_opt.pfreq),
    ];

    // Define arguments that set something to a string.

    let set_string = [
        ("CLUSTAL_AA", &mut ctl.gen_opt.clustal_aa),
        ("CLUSTAL_DNA", &mut ctl.gen_opt.clustal_dna),
        ("EXT", &mut ctl.gen_opt.ext),
        ("PHYLIP_AA", &mut ctl.gen_opt.phylip_aa),
        ("PHYLIP_DNA", &mut ctl.gen_opt.phylip_dna),
        ("POUT", &mut ctl.parseable_opt.pout),
        ("REF", &mut ctl.gen_opt.refname),
        ("TRACE_BARCODE", &mut ctl.gen_opt.trace_barcode),
    ];

    // Define arguments that set something to a string that is an output file name.

    let set_string_writeable = [
        ("BINARY", &mut ctl.gen_opt.binary),
        ("DONOR_REF_FILE", &mut ctl.gen_opt.dref_file),
        ("PEER_GROUP", &mut ctl.gen_opt.peer_group_filename),
        ("PROTO", &mut ctl.gen_opt.proto),
        ("SUBSET_JSON", &mut ctl.gen_opt.subset_json),
    ];

    // Define arguments that set something to a string that is an input file name.

    let set_string_readable = [("PROTO_METADATA", &mut ctl.gen_opt.proto_metadata)];

    // Define arguments that do nothing (because already parsed), and which have no "= value" part.

    let set_nothing_simple = [
        "CELLRANGER",
        "COMP",
        "COMP2",
        "CTRLC",
        "DUMP_INTERNAL_IDS",
        "FORCE_EXTERNAL",
        "LONG_HELP",
        "MARKED_B",
        "MARK_STATS",
        "MARK_STATS2",
        "NALL",
        "NALL_CELL",
        "NALL_GEX",
        "NOPAGER",
        "NOPRETTY",
        "PLAIN",
        "PRINT_CPU",
        "PRINT_CPU_INFO",
        "SVG",
    ];

    // Define arguments that do nothing (because already parsed), and which may have
    // an "= value" part.

    let set_nothing = ["BC", "BI", "EMAIL", "GEX", "HAPS", "HTML", "PRE"];

    // Traverse arguments.

    'args_loop: for i in 1..args.len() {
        let mut arg = args[i].to_string();

        // Replace deprecated option.

        if arg == "KEEP_IMPROPER".to_string() {
            arg = "NIMPROPER".to_string();
        }

        // Strip out certain quoted expressions.

        if arg.contains("=\"") && arg.ends_with("\"") {
            let mut quotes = 0;
            for c in arg.chars() {
                if c == '\"' {
                    quotes += 1;
                }
            }
            if quotes == 2 {
                arg = format!("{}={}", arg.before("="), arg.between("\"", "\""));
            }
        }

        // Check for weird case that might arise if testing code is screwed up.

        if arg.len() == 0 {
            eprintln!(
                "\nYou've passed a null argument to enclone.  Normally that isn't \
                 possible.\nPlease take a detailed look at how you're invoking enclone.\n"
            );
            std::process::exit(1);
        }

        // Process set_true arguments.

        for j in 0..set_true.len() {
            if arg == set_true[j].0.to_string() {
                *(set_true[j].1) = true;
                continue 'args_loop;
            }
        }

        // Process set_false arguments.

        for j in 0..set_false.len() {
            if arg == set_false[j].0.to_string() {
                *(set_false[j].1) = false;
                continue 'args_loop;
            }
        }

        // Process set_usize args.

        for j in 0..set_usize.len() {
            if is_usize_arg(&arg, &set_usize[j].0) {
                *(set_usize[j].1) = arg.after(&format!("{}=", set_usize[j].0)).force_usize();
                continue 'args_loop;
            }
        }

        // Process set_string args.

        for j in 0..set_string.len() {
            if is_string_arg(&arg, &set_string[j].0) {
                *(set_string[j].1) = arg.after(&format!("{}=", set_string[j].0)).to_string();
                continue 'args_loop;
            }
        }

        // Process set_string_writeable args.

        for j in 0..set_string_writeable.len() {
            let var = &set_string_writeable[j].0;
            if is_string_arg(&arg, var) {
                *(set_string_writeable[j].1) = arg.after(&format!("{}=", var)).to_string();
                let val = &set_string_writeable[j].1;
                let f = File::create(&val);
                if f.is_err() {
                    eprintln!(
                        "\nYou've specified an output file\n{}\nthat cannot be written.",
                        val
                    );
                    if val.contains("/") {
                        let dir = val.rev_before("/");
                        let msg;
                        if path_exists(&dir) {
                            msg = "exists";
                        } else {
                            msg = "does not exist";
                        }
                        eprintln!("Note that the path {} {}.", dir, msg);
                    }
                    eprintln!("");
                    std::process::exit(1);
                }
                remove_file(&val).unwrap();
                continue 'args_loop;
            }
        }

        // Process set_string_readable args.

        for j in 0..set_string_readable.len() {
            let var = &set_string_readable[j].0;
            if is_string_arg(&arg, var) {
                let val = arg.after(&format!("{}=", var));
                if val.is_empty() {
                    eprintln!("\nFilename input in {} cannot be empty\n", val);
                    std::process::exit(1);
                }
                *(set_string_readable[j].1) = Some(val.to_string());
                if let Err(e) = File::open(&val) {
                    eprintln!(
                        "\nYou've specified an input file\n{}\nthat cannot be read due to {}\n",
                        val, e
                    );
                    std::process::exit(1);
                }
                continue 'args_loop;
            }
        }

        // Process set_nothing_simple args.

        for j in 0..set_nothing_simple.len() {
            if arg == set_nothing_simple[j].to_string() {
                continue 'args_loop;
            }
        }

        // Process set_nothing args.

        for j in 0..set_nothing.len() {
            if arg == set_nothing[j].to_string() || arg.starts_with(&format!("{}=", set_nothing[j]))
            {
                continue 'args_loop;
            }
        }

        // Process the argument.

        if is_simple_arg(&arg, "SEQ") {
            ctl.join_print_opt.seq = true;

        // Not movable.
        } else if arg.starts_with("PG_DIST=") {
            let dist = arg.after("PG_DIST=");
            if dist != "MFL" {
                eprintln!("\nCurrently the only allowed value for PG_DIST is MFL.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.peer_group_dist = dist.to_string();
        } else if is_simple_arg(&arg, "H5") {
            ctl.gen_opt.force_h5 = true;
        } else if is_simple_arg(&arg, "NH5") {
            ctl.gen_opt.force_h5 = false;
        } else if arg == "LEGEND" {
            ctl.gen_opt.use_legend = true;
        } else if is_usize_arg(&arg, "REQUIRED_FPS") {
            ctl.gen_opt.required_fps = Some(arg.after("REQUIRED_FPS=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_CELLS") {
            ctl.gen_opt.required_cells = Some(arg.after("REQUIRED_CELLS=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_DONORS") {
            ctl.gen_opt.required_donors = Some(arg.after("REQUIRED_DONORS=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_CLONOTYPES") {
            ctl.gen_opt.required_clonotypes = Some(arg.after("REQUIRED_CLONOTYPES=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_TWO_CELL_CLONOTYPES") {
            ctl.gen_opt.required_two_cell_clonotypes =
                Some(arg.after("REQUIRED_TWO_CELL_CLONOTYPES=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_TWO_CHAIN_CLONOTYPES") {
            ctl.gen_opt.required_two_chain_clonotypes =
                Some(arg.after("REQUIRED_TWO_CHAIN_CLONOTYPES=").force_usize());
        } else if is_usize_arg(&arg, "REQUIRED_DATASETS") {
            ctl.gen_opt.required_datasets = Some(arg.after("REQUIRED_DATASETS=").force_usize());
        } else if is_usize_arg(&arg, "EXACT") {
            ctl.gen_opt.exact = Some(arg.after("EXACT=").force_usize());
        } else if is_usize_arg(&arg, "MIN_CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("MIN_CHAINS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CHAINS") {
            ctl.clono_filt_opt.max_chains = arg.after("MAX_CHAINS=").force_usize();
        } else if is_usize_arg(&arg, "MIN_CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("MIN_CELLS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CELLS") {
            ctl.clono_filt_opt.ncells_high = arg.after("MAX_CELLS=").force_usize();
        } else if arg.starts_with("EXFASTA=") {
            ctl.gen_opt.fasta = arg.after("EXFASTA=").to_string();
        } else if arg.starts_with("FASTA=") {
            ctl.gen_opt.fasta_filename = arg.after("FASTA=").to_string();
        } else if arg.starts_with("FASTA_AA=") {
            ctl.gen_opt.fasta_aa_filename = arg.after("FASTA_AA=").to_string();

        // Other.
        } else if arg.starts_with("DIFF_STYLE=") {
            ctl.gen_opt.diff_style = arg.after("=").to_string();
            if ctl.gen_opt.diff_style != "C1" && ctl.gen_opt.diff_style != "C2" {
                eprintln!("\nThe only allowed values for DIFF_STYLE are C1 and C2.\n");
                std::process::exit(1);
            }
        } else if arg.starts_with("COLOR=") {
            ctl.gen_opt.color = arg.after("COLOR=").to_string();
            if ctl.gen_opt.color != "codon".to_string()
                && ctl.gen_opt.color != "property".to_string()
            {
                let mut ok = false;
                if arg.starts_with("COLOR=peer.") {
                    let pc = arg.after("COLOR=peer.");
                    if pc.parse::<f64>().is_ok() {
                        let pc = pc.force_f64();
                        if pc >= 0.0 && pc <= 100.0 {
                            ok = true;
                            ctl.gen_opt.color_by_rarity_pc = pc;
                        }
                    }
                }
                if !ok {
                    eprintln!(
                        "The specified value for COLOR is not allowed.  Please see \
                        \"enclone help color\".\n"
                    );
                    std::process::exit(1);
                }
            }
        } else if arg == "TREE" {
            ctl.gen_opt.tree = ".".to_string();
        } else if arg == "TREE=const" {
            ctl.gen_opt.tree = "const".to_string();
        } else if arg.starts_with("FCELL=") {
            let mut condition = arg.after("FCELL=").to_string();
            let con = condition.as_bytes();
            for i in 0..con.len() {
                if i > 0 && i < con.len() - 1 && con[i] == b'=' {
                    if con[i - 1] != b'=' && con[i + 1] != b'=' {
                        eprintln!(
                            "\nConstraints for FCELL cannot use =.  Please use == instead.\n"
                        );
                        std::process::exit(1);
                    }
                }
            }
            condition = condition.replace("'", "\"");
            let compiled = build_operator_tree(&condition);
            if !compiled.is_ok() {
                eprintln!("\nFCELL usage incorrect.\n");
                std::process::exit(1);
            }
            ctl.clono_filt_opt.fcell.push(compiled.unwrap());
        } else if is_simple_arg(&arg, "FAIL_ONLY=true") {
            ctl.clono_filt_opt.fail_only = true;
        } else if arg.starts_with("LEGEND=") {
            let x = parse_csv(&arg.after("LEGEND="));
            if x.len() == 0 || x.len() % 2 != 0 {
                eprintln!("\nValue of LEGEND doesn't make sense.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.use_legend = true;
            for i in 0..x.len() / 2 {
                ctl.gen_opt
                    .legend
                    .push((x[2 * i].clone(), x[2 * i + 1].clone()));
            }
        } else if arg.starts_with("BARCODE=") {
            let bcs = arg.after("BARCODE=").split(',').collect::<Vec<&str>>();
            let mut x = Vec::<String>::new();
            for j in 0..bcs.len() {
                if !bcs[j].contains('-') {
                    eprintln!(
                        "\nValue for a barcode in BARCODE argument is invalid, must contain -.\n"
                    );
                    std::process::exit(1);
                }
                x.push(bcs[j].to_string());
            }
            ctl.clono_filt_opt.barcode = x;
        } else if arg.starts_with("F=") {
            let filt = arg.after("F=").to_string();
            ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt));
        } else if arg.starts_with("SCAN=") {
            let mut x = arg.after("SCAN=").to_string();
            x = x.replace(" ", "").to_string();
            let x = x.split(',').collect::<Vec<&str>>();
            if x.len() != 3 {
                eprintln!("\nArgument to SCAN must have three components.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.gene_scan_test = Some(LinearCondition::new(&x[0]));
            ctl.gen_opt.gene_scan_control = Some(LinearCondition::new(&x[1]));
            let threshold = LinearCondition::new(&x[2]);
            for i in 0..threshold.var.len() {
                if threshold.var[i] != "t".to_string() && threshold.var[i] != "c".to_string() {
                    eprintln!("\nIllegal variable in threshold for scan.\n");
                    std::process::exit(1);
                }
            }
            ctl.gen_opt.gene_scan_threshold = Some(threshold);
        } else if arg.starts_with("PLOT=") {
            using_plot = true;
            let x = arg.after("PLOT=").split(',').collect::<Vec<&str>>();
            if x.is_empty() {
                eprintln!("\nArgument to PLOT is invalid.\n");
                std::process::exit(1);
            }
            ctl.gen_opt.plot_file = x[0].to_string();
            for j in 1..x.len() {
                if !x[j].contains("->") {
                    eprintln!("\nArgument to PLOT is invalid.\n");
                    std::process::exit(1);
                }
                ctl.gen_opt
                    .origin_color_map
                    .insert(x[j].before("->").to_string(), x[j].after("->").to_string());
            }
        } else if arg.starts_with("PLOT_BY_ISOTYPE=") {
            ctl.gen_opt.plot_by_isotype = true;
            ctl.gen_opt.plot_file = arg.after("PLOT_BY_ISOTYPE=").to_string();
            if ctl.gen_opt.plot_file.is_empty() {
                eprintln!("\nFilename value needs to be supplied to PLOT_BY_ISOTYPE.\n");
                std::process::exit(1);
            }
        } else if arg.starts_with("PLOT_BY_MARK=") {
            ctl.gen_opt.plot_by_mark = true;
            ctl.gen_opt.plot_file = arg.after("PLOT_BY_MARK=").to_string();
            if ctl.gen_opt.plot_file.is_empty() {
                eprintln!("\nFilename value needs to be supplied to PLOT_BY_MARK.\n");
                std::process::exit(1);
            }
        } else if is_simple_arg(&arg, "FAIL_ONLY=false") {
            ctl.clono_filt_opt.fail_only = false;
        } else if is_usize_arg(&arg, "MAX_CORES") {
            let nthreads = arg.after("MAX_CORES=").force_usize();
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build_global();
        } else if arg.starts_with("PCOLS=") {
            ctl.parseable_opt.pcols.clear();
            let p = arg.after("PCOLS=").split(',').collect::<Vec<&str>>();
            for i in 0..p.len() {
                let mut x = p[i].to_string();
                x = x.replace("_sum", "_Σ");
                x = x.replace("_mean", "_μ");
                ctl.parseable_opt.pcols.push(x.to_string());
                ctl.parseable_opt.pcols_sort = ctl.parseable_opt.pcols.clone();
                ctl.parseable_opt.pcols_sortx = ctl.parseable_opt.pcols.clone();
                for j in 0..ctl.parseable_opt.pcols_sortx.len() {
                    if ctl.parseable_opt.pcols_sortx[j].contains(":") {
                        ctl.parseable_opt.pcols_sortx[j] =
                            ctl.parseable_opt.pcols_sortx[j].before(":").to_string();
                    }
                }
                unique_sort(&mut ctl.parseable_opt.pcols_sort);
                unique_sort(&mut ctl.parseable_opt.pcols_sortx);
            }
        } else if arg.starts_with("VJ=") {
            ctl.clono_filt_opt.vj = arg.after("VJ=").as_bytes().to_vec();
            for c in ctl.clono_filt_opt.vj.iter() {
                if !(*c == b'A' || *c == b'C' || *c == b'G' || *c == b'T') {
                    eprintln!("\nIllegal value for VJ, must be over alphabet ACGT.\n");
                    std::process::exit(1);
                }
            }
        } else if arg.starts_with("AMINO=") {
            ctl.clono_print_opt.amino.clear();
            for x in arg.after("AMINO=").split(',').collect::<Vec<&str>>() {
                if x != "" {
                    ctl.clono_print_opt.amino.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.amino.iter() {
                let mut ok = false;
                if *x == "cdr1"
                    || *x == "cdr2"
                    || *x == "cdr3"
                    || *x == "fwr1"
                    || *x == "fwr2"
                    || *x == "fwr3"
                    || *x == "fwr4"
                    || *x == "var"
                    || *x == "share"
                    || *x == "donor"
                    || *x == "donorn"
                {
                    ok = true;
                } else if x.contains('-') {
                    let (start, stop) = (x.before("-"), x.after("-"));
                    if start.parse::<usize>().is_ok() && stop.parse::<usize>().is_ok() {
                        if start.force_usize() <= stop.force_usize() {
                            ok = true;
                        }
                    }
                }
                if !ok {
                    eprintln!(
                        "\nUnrecognized variable {} for AMINO.  Please type \
                         \"enclone help amino\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
        } else if arg.starts_with("CVARS=") {
            ctl.clono_print_opt.cvars.clear();
            for x in arg.after("CVARS=").split(',').collect::<Vec<&str>>() {
                if x.len() > 0 {
                    ctl.clono_print_opt.cvars.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.cvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("CVARSP=") {
            for x in arg.after("CVARSP=").split(',').collect::<Vec<&str>>() {
                if x.len() > 0 {
                    ctl.clono_print_opt.cvars.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.cvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("LVARS=") {
            ctl.clono_print_opt.lvars.clear();
            for x in arg.after("LVARS=").split(',').collect::<Vec<&str>>() {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
            for x in ctl.clono_print_opt.lvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if arg.starts_with("LVARSP=") {
            let lvarsp = arg.after("LVARSP=").split(',').collect::<Vec<&str>>();
            for x in lvarsp {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
            for x in ctl.clono_print_opt.lvars.iter_mut() {
                *x = x.replace("_sum", "_Σ");
                *x = x.replace("_mean", "_μ");
            }
        } else if is_f64_arg(&arg, "MAX_SCORE") {
            ctl.join_alg_opt.max_score = arg.after("MAX_SCORE=").force_f64();
        } else if is_f64_arg(&arg, "MAX_LOG_SCORE") {
            let x = arg.after("MAX_LOG_SCORE=").force_f64();
            ctl.join_alg_opt.max_score = 10.0_f64.powf(x);
        } else if arg.starts_with("CONST_IGH=") {
            let reg = Regex::new(&format!("^{}$", arg.after("CONST_IGH=")));
            if !reg.is_ok() {
                eprintln!(
                    "\nYour CONST_IGH value {} could not be parsed as a regular expression.\n",
                    arg.after("CONST_IGH=")
                );
                std::process::exit(1);
            }
            ctl.gen_opt.const_igh = Some(reg.unwrap());
        } else if arg.starts_with("CONST_IGKL=") {
            let reg = Regex::new(&format!("^{}$", arg.after("CONST_IGKL=")));
            if !reg.is_ok() {
                eprintln!(
                    "\nYour CONST_IGKL value {} could not be parsed as a regular expression.\n",
                    arg.after("CONST_IGKL=")
                );
                std::process::exit(1);
            }
            ctl.gen_opt.const_igkl = Some(reg.unwrap());
        } else if arg.starts_with("CDR3=") {
            let fields = arg.split('|').collect::<Vec<&str>>();
            let mut lev = true;
            for i in 0..fields.len() {
                if !Regex::new(r"[A-Z]+~[0-9]+")
                    .as_ref()
                    .unwrap()
                    .is_match(fields[i])
                {
                    lev = false;
                }
            }
            if lev {
                ctl.clono_filt_opt.cdr3_lev = arg.after("=").to_string();
            } else {
                let reg = Regex::new(&format!("^{}$", arg.after("CDR3=")));
                if !reg.is_ok() {
                    eprintln!(
                        "\nYour CDR3 value {} could not be parsed as a regular expression.\n",
                        arg.after("CDR3=")
                    );
                    std::process::exit(1);
                }
                ctl.clono_filt_opt.cdr3 = Some(reg.unwrap());
            }
        } else if is_usize_arg(&arg, "CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("CHAINS=").force_usize();
            ctl.clono_filt_opt.max_chains = arg.after("CHAINS=").force_usize();
        } else if arg.starts_with("SEG=") {
            let fields = arg.after("SEG=").split('|').collect::<Vec<&str>>();
            let mut y = Vec::<String>::new();
            for x in fields.iter() {
                y.push(x.to_string());
            }
            y.sort();
            ctl.clono_filt_opt.seg.push(y);
        } else if arg.starts_with("SEGN=") {
            let fields = arg.after("SEGN=").split('|').collect::<Vec<&str>>();
            let mut y = Vec::<String>::new();
            for x in fields.iter() {
                if !x.parse::<i32>().is_ok() {
                    eprintln!("\nInvalid argument to SEGN.\n");
                    std::process::exit(1);
                }
                y.push(x.to_string());
            }
            y.sort();
            ctl.clono_filt_opt.segn.push(y);
        } else if is_usize_arg(&arg, "CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("CELLS=").force_usize();
            ctl.clono_filt_opt.ncells_high = ctl.clono_filt_opt.ncells_low;
        } else if arg.starts_with("META=") {
            let f = arg.after("META=");
            metas.push(f.to_string());
        } else if arg.starts_with("TCR=")
            || arg.starts_with("BCR=")
            || (arg.len() > 0 && arg.as_bytes()[0] >= b'0' && arg.as_bytes()[0] <= b'9')
        {
            xcrs.push(arg.to_string());
        } else {
            eprintln!("\nUnrecognized argument {}.\n", arg);
            std::process::exit(1);
        }
    }
    ctl.perf_stats(&targs, "in main args loop");

    // Expand ~ and ~user in output file names.

    let t = Instant::now();
    let mut files = [
        &mut ctl.gen_opt.plot_file,
        &mut ctl.gen_opt.fasta_filename,
        &mut ctl.gen_opt.fasta_aa_filename,
        &mut ctl.gen_opt.dref_file,
        &mut ctl.parseable_opt.pout,
    ];
    for f in files.iter_mut() {
        **f = stringme(&tilde_expand(&f.as_bytes()));
    }

    // Sanity check arguments (and more below).

    if ctl.clono_filt_opt.cdr3.is_some() && ctl.clono_filt_opt.cdr3_lev.len() > 0 {
        eprintln!(
            "\nPlease use the CDR3 argument to specify either a regular expression or a\n\
            Levenshtein distance pattern, but not both.\n"
        );
        std::process::exit(1);
    }
    if ctl.gen_opt.clustal_aa != "".to_string() && ctl.gen_opt.clustal_aa != "stdout".to_string() {
        if !ctl.gen_opt.clustal_aa.ends_with(".tar") {
            eprintln!("\nIf the value of CLUSTAL_AA is not stdout, it must end in .tar.\n");
            std::process::exit(1);
        }
    }
    if ctl.gen_opt.clustal_dna != "".to_string() && ctl.gen_opt.clustal_dna != "stdout".to_string()
    {
        if !ctl.gen_opt.clustal_dna.ends_with(".tar") {
            eprintln!("\nIf the value of CLUSTAL_DNA is not stdout, it must end in .tar.\n");
            std::process::exit(1);
        }
    }
    if ctl.gen_opt.phylip_aa != "".to_string() && ctl.gen_opt.phylip_aa != "stdout".to_string() {
        if !ctl.gen_opt.phylip_aa.ends_with(".tar") {
            eprintln!("\nIf the value of PHYLIP_AA is not stdout, it must end in .tar.\n");
            std::process::exit(1);
        }
    }
    if ctl.gen_opt.phylip_dna != "".to_string() && ctl.gen_opt.phylip_dna != "stdout".to_string() {
        if !ctl.gen_opt.phylip_dna.ends_with(".tar") {
            eprintln!("\nIf the value of PHYLIP_DNA is not stdout, it must end in .tar.\n");
            std::process::exit(1);
        }
    }
    if ctl.clono_filt_opt.umi_filt && ctl.clono_filt_opt.umi_filt_mark {
        eprintln!(
            "\nIf you use UMI_FILT_MARK, you should also use NUMI, to turn off \
            the filter,\nas otherwise nothing will be marked.\n"
        );
        std::process::exit(1);
    }
    if ctl.clono_filt_opt.umi_ratio_filt && ctl.clono_filt_opt.umi_ratio_filt_mark {
        eprintln!(
            "\nIf you use UMI_RATIO_FILT_MARK, you should also use NUMI_RATIO, to turn off \
            the filter,\nas otherwise nothing will be marked.\n"
        );
        std::process::exit(1);
    }
    ctl.perf_stats(&t, "after main args loop 1");

    // Process TCR, BCR and META.

    let t = Instant::now();
    check_cvars(&ctl);
    if metas.len() > 0 {
        let f = &metas[metas.len() - 1];
        let f = get_path_fail(&f, &ctl, "META");
        proc_meta(&f, &mut ctl);
    }
    ctl.perf_stats(&t, "in proc_meta");
    if xcrs.len() > 0 {
        let arg = &xcrs[xcrs.len() - 1];
        proc_xcr(&arg, &gex, &bc, have_gex, &mut ctl);
    }

    // More argument sanity checking.

    let bcr_only = [
        "PEER_GROUP",
        "PG_READABLE",
        "PG_DIST",
        "COLOR=peer",
        "CONST_IGH",
        "CONST_IGL",
    ];
    if !ctl.gen_opt.bcr {
        for i in 1..args.len() {
            let arg = &args[i];
            for x in bcr_only.iter() {
                if arg == x || arg.starts_with(&format!("{}=", x)) {
                    eprintln!("\nThe option {} does not make sense for TCR.\n", x);
                    std::process::exit(1);
                }
            }
        }
    }

    // Proceed.

    let t = Instant::now();
    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);
    for con in ctl.clono_filt_opt.fcell.iter() {
        for var in con.iter_variable_identifiers() {
            if !bin_member(&alt_bcs, &var.to_string()) {
                eprintln!(
                    "\nYou've used a variable {} as part of an FCELL argument that has not\n\
                    been specified using BC or bc (via META).\n",
                    var
                );
                std::process::exit(1);
            }
        }
        for _ in con.iter_function_identifiers() {
            eprintln!("\nSomething is wrong with your FCELL value.\n");
            std::process::exit(1);
        }
    }
    for i in 0..ctl.origin_info.n() {
        let (mut cells_cr, mut rpc_cr) = (None, None);
        if ctl.gen_opt.internal_run {
            let p = &ctl.origin_info.dataset_path[i];
            let mut f = format!("{}/metrics_summary_csv.csv", p);
            if !path_exists(&f) {
                f = format!("{}/metrics_summary.csv", p);
            }
            if path_exists(&f) {
                let f = open_for_read![&f];
                let mut count = 0;
                let (mut cells_field, mut rpc_field) = (None, None);
                for line in f.lines() {
                    count += 1;
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    for (i, x) in fields.iter().enumerate() {
                        if count == 1 {
                            if *x == "Estimated Number of Cells" {
                                cells_field = Some(i);
                            } else if *x == "Mean Read Pairs per Cell" {
                                rpc_field = Some(i);
                            }
                        } else if count == 2 {
                            if Some(i) == cells_field {
                                let mut n = x.to_string();
                                if n.contains("\"") {
                                    n = n.between("\"", "\"").to_string();
                                }
                                n = n.replace(",", "");
                                cells_cr = Some(n.force_usize());
                            } else if Some(i) == rpc_field {
                                let mut n = x.to_string();
                                if n.contains("\"") {
                                    n = n.between("\"", "\"").to_string();
                                }
                                n = n.replace(",", "");
                                rpc_cr = Some(n.force_usize());
                            }
                        }
                    }
                }
            }
        }
        ctl.origin_info.cells_cellranger.push(cells_cr);
        ctl.origin_info
            .mean_read_pairs_per_cell_cellranger
            .push(rpc_cr);
    }
    if ctl.gen_opt.plot_by_isotype {
        if using_plot || ctl.gen_opt.use_legend {
            eprintln!("\nPLOT_BY_ISOTYPE cannot be used with PLOT or LEGEND.\n");
            std::process::exit(1);
        }
        if !ctl.gen_opt.bcr {
            eprintln!("\nPLOT_BY_ISOTYPE can only be used with BCR data.\n");
            std::process::exit(1);
        }
        if ctl.gen_opt.plot_by_mark {
            eprintln!("\nPLOT_BY_ISOTYPE and PLOT_BY_MARK cannot be used together.\n");
            std::process::exit(1);
        }
    }
    if ctl.gen_opt.plot_by_mark {
        if using_plot || ctl.gen_opt.use_legend {
            eprintln!("\nPLOT_BY_MARK cannot be used with PLOT or LEGEND.\n");
            std::process::exit(1);
        }
    }
    if ctl.parseable_opt.pbarcode && ctl.parseable_opt.pout.len() == 0 {
        eprintln!("\nIt does not make sense to specify PCELL unless POUT is also specified.\n");
        std::process::exit(1);
    }
    if ctl.origin_info.n() == 0 {
        eprintln!("\nNo TCR or BCR data have been specified.\n");
        std::process::exit(1);
    }
    let mut donors = Vec::<String>::new();
    let mut origins = Vec::<String>::new();
    let mut tags = Vec::<String>::new();
    let mut origin_for_bc = Vec::<String>::new();
    let mut donor_for_bc = Vec::<String>::new();
    for i in 0..ctl.origin_info.n() {
        for x in ctl.origin_info.origin_for_bc[i].iter() {
            origins.push(x.1.clone());
            origin_for_bc.push(x.1.clone());
        }
        for x in ctl.origin_info.donor_for_bc[i].iter() {
            donors.push(x.1.clone());
            donor_for_bc.push(x.1.clone());
        }
        for x in ctl.origin_info.tag[i].iter() {
            tags.push((x.1).clone());
        }
        donors.push(ctl.origin_info.donor_id[i].clone());
        origins.push(ctl.origin_info.origin_id[i].clone());
    }
    unique_sort(&mut donors);
    unique_sort(&mut origins);
    unique_sort(&mut tags);
    unique_sort(&mut origin_for_bc);
    unique_sort(&mut donor_for_bc);
    ctl.origin_info.donors = donors.len();
    ctl.origin_info.dataset_list = ctl.origin_info.dataset_id.clone();
    unique_sort(&mut ctl.origin_info.dataset_list);
    ctl.origin_info.origin_list = origins.clone();
    ctl.origin_info.donor_list = donors.clone();
    ctl.origin_info.tag_list = tags;
    for i in 0..ctl.origin_info.donor_for_bc.len() {
        if ctl.origin_info.donor_for_bc[i].len() > 0 {
            ctl.clono_filt_opt.donor = true;
        }
    }
    ctl.perf_stats(&t, "after main args loop 2");
    proc_args_tail(&mut ctl, &args);

    // Check for invalid variables in linear conditions.

    for i in 0..ctl.clono_filt_opt.bounds.len() {
        ctl.clono_filt_opt.bounds[i].require_valid_variables(&ctl);
    }
    if ctl.gen_opt.gene_scan_test.is_some() {
        ctl.gen_opt
            .gene_scan_test
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl);
        ctl.gen_opt
            .gene_scan_control
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl);
    }
}
