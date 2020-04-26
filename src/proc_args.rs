// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use crate::proc_args2::*;
use crate::proc_args3::*;
use crate::proc_args_check::*;
use crate::testlist::*;
use perf_stats::*;
use regex::Regex;
use std::{env, time::Instant};
use string_utils::*;
use vector_utils::*;

// Process arguments.

pub fn proc_args(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    // Knobs.

    let heur = ClonotypeHeuristics {
        max_diffs: 50,
        ref_v_trim: 15,
        ref_j_trim: 15,
    };
    ctl.heur = heur;

    // Form the combined set of command-line arguments and "command-line" arguments
    // implied by environment variables.

    let targs = Instant::now();
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
        ctl.gen_opt.pre = vec![format!("/mnt/assembly/vdj/current{}", TEST_FILES_VERSION)];
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
    ctl.silent = true;

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

    ctl.join_print_opt.pfreq = 1_000_000_000;
    ctl.join_print_opt.quiet = true;

    ctl.parseable_opt.pchains = 4;

    ctl.onesie_mult = 10_000;

    // Pretest for consistency amongst TCR, BCR, GEX and META.  Also preparse GEX.

    let mut have_tcr = false;
    let mut have_bcr = false;
    let mut have_gex = false;
    let mut have_meta = false;
    let mut gex = String::new();
    let mut bc = String::new();
    for i in 1..args.len() {
        if args[i].starts_with("TCR=") {
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
    }
    if have_meta && (have_tcr || have_bcr || have_gex || bc.len() > 0) {
        eprintln!("\nIf META is specified, then none of TCR, BCR, GEX or BC can be specified.\n");
        std::process::exit(1);
    }
    if have_tcr && have_bcr {
        eprintln!("\nKindly please do not specify both TCR and BCR.\n");
        std::process::exit(1);
    }

    // Traverse arguments.

    for i in 1..args.len() {
        let mut arg = args[i].to_string();

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

        // Process the argument.

        if is_simple_arg(&arg, "SEQ") {
            ctl.join_print_opt.seq = true;
        } else if is_simple_arg(&arg, "ANN") {
            ctl.join_print_opt.ann = true;
        } else if is_simple_arg(&arg, "ANN0") {
            ctl.join_print_opt.ann0 = true;
        } else if is_simple_arg(&arg, "DUMP_LENAS") {
        } else if is_simple_arg(&arg, "SHOW_BC") {
            ctl.join_print_opt.show_bc = true;
        } else if is_simple_arg(&arg, "PER_CELL") {
            ctl.clono_print_opt.bu = true;
        } else if is_simple_arg(&arg, "COMP") {
        } else if is_simple_arg(&arg, "COMP2") {
        } else if is_simple_arg(&arg, "LONG_HELP") {
        } else if is_simple_arg(&arg, "CON") {
            ctl.allele_print_opt.con = true;
        } else if is_simple_arg(&arg, "CON_TRACE") {
            ctl.allele_print_opt.con_trace = true;
        } else if is_simple_arg(&arg, "EXP") {
            ctl.gen_opt.exp = true;
        } else if is_simple_arg(&arg, "JC1") {
            ctl.gen_opt.jc1 = true;
        } else if is_simple_arg(&arg, "NGROUP") {
            ctl.gen_opt.ngroup = true;
        } else if is_simple_arg(&arg, "STABLE_DOC") {
            ctl.gen_opt.stable_doc = true;
        } else if is_simple_arg(&arg, "IMGT") {
            ctl.gen_opt.imgt = true;
        } else if is_simple_arg(&arg, "IMGT_FIX") {
            ctl.gen_opt.imgt_fix = true;
        } else if arg == "LEGEND" {
            ctl.gen_opt.use_legend = true;
        } else if arg == "HTML" {
        } else if arg == "SVG" {
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
        } else if is_simple_arg(&arg, "H5") {
            ctl.gen_opt.force_h5 = true;
        } else if is_simple_arg(&arg, "CURRENT_REF") {
            ctl.gen_opt.current_ref = true;
        } else if is_simple_arg(&arg, "SUM") {
            ctl.clono_print_opt.sum = true;
        } else if is_simple_arg(&arg, "MEAN") {
            ctl.clono_print_opt.mean = true;
        } else if is_simple_arg(&arg, "NH5") {
            ctl.gen_opt.force_h5 = false;
        } else if is_simple_arg(&arg, "H5_SLICE") {
            ctl.gen_opt.h5_pre = false;
        } else if is_simple_arg(&arg, "DESCRIP") {
            ctl.gen_opt.descrip = true;
        } else if is_simple_arg(&arg, "CTRLC") {
        } else if is_simple_arg(&arg, "FORCE") {
            ctl.force = true;
        } else if is_simple_arg(&arg, "CELLRANGER") {
        } else if is_simple_arg(&arg, "WEAK") {
            ctl.gen_opt.weak = true;
        } else if is_simple_arg(&arg, "REUSE") {
            ctl.gen_opt.reuse = true;
        } else if is_simple_arg(&arg, "FORCE_EXTERNAL") {
        } else if is_simple_arg(&arg, "NWARN") {
            ctl.gen_opt.nwarn = true;
        } else if arg.starts_with("BINARY=") {
            ctl.gen_opt.binary = arg.after("BINARY=").to_string();
        } else if arg.starts_with("PROTO=") {
            ctl.gen_opt.proto = arg.after("PROTO=").to_string();
        } else if is_simple_arg(&arg, "PRINT_FAILED_JOINS") {
            ctl.join_print_opt.quiet = false;
        } else if is_simple_arg(&arg, "NOTE_SIMPLE") {
            ctl.clono_print_opt.note_simple = true;
        } else if is_simple_arg(&arg, "SEQC") {
            ctl.clono_print_opt.seqc = true;
        } else if is_simple_arg(&arg, "FULL_SEQC") {
            ctl.clono_print_opt.full_seqc = true;
        } else if is_simple_arg(&arg, "BARCODES") {
            ctl.clono_print_opt.barcodes = true;
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
        } else if is_simple_arg(&arg, "GRAPH") {
            ctl.gen_opt.graph = true;
        } else if is_simple_arg(&arg, "ACCEPT_INCONSISTENT") {
            ctl.gen_opt.accept_inconsistent = true;
        } else if is_simple_arg(&arg, "NGEX") {
            ctl.clono_filt_opt.ngex = true;
        } else if is_simple_arg(&arg, "NCROSS") {
            ctl.clono_filt_opt.ncross = true;
        } else if is_simple_arg(&arg, "NWEAK_CHAINS") {
            ctl.clono_filt_opt.weak_chains = false;
        } else if is_simple_arg(&arg, "NWEAK_ONESIES") {
            ctl.clono_filt_opt.weak_onesies = false;
        } else if is_simple_arg(&arg, "NFOURSIE_KILL") {
            ctl.clono_filt_opt.weak_foursies = false;
        } else if is_simple_arg(&arg, "NBC_DUP") {
            ctl.clono_filt_opt.bc_dup = false;
        } else if is_simple_arg(&arg, "MIX_DONORS") {
            ctl.clono_filt_opt.donor = true;
        } else if is_simple_arg(&arg, "HAVE_ONESIE") {
            ctl.clono_filt_opt.have_onesie = true;
        } else if is_simple_arg(&arg, "UTR_CON") {
            ctl.gen_opt.utr_con = true;
        } else if is_simple_arg(&arg, "CON_CON") {
            ctl.gen_opt.con_con = true;
        } else if is_simple_arg(&arg, "MOUSE") {
            ctl.gen_opt.mouse = true;
        } else if is_simple_arg(&arg, "SUMMARY") {
            ctl.gen_opt.summary = true;
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
                    .sample_color_map
                    .insert(x[j].before("->").to_string(), x[j].after("->").to_string());
            }
        } else if is_simple_arg(&arg, "SUMMARY_CLEAN") {
            ctl.gen_opt.summary_clean = true;
        } else if arg.starts_with("EMAIL=") {
        } else if arg.starts_with("REF=") {
            ctl.gen_opt.refname = arg.after("REF=").to_string();
        } else if is_simple_arg(&arg, "NSILENT") {
            ctl.silent = false;
        } else if is_simple_arg(&arg, "TOY") {
            ctl.toy = true;
        } else if is_simple_arg(&arg, "RE") {
            ctl.gen_opt.reannotate = true;
        } else if is_simple_arg(&arg, "WHITEF") {
            ctl.clono_filt_opt.whitef = true;
        } else if is_simple_arg(&arg, "PROTECT_BADS") {
            ctl.clono_filt_opt.protect_bads = true;
        } else if is_simple_arg(&arg, "NWHITEF") {
            ctl.gen_opt.nwhitef = true;
        } else if is_simple_arg(&arg, "FAIL_ONLY=true") {
            ctl.clono_filt_opt.fail_only = true;
        } else if is_simple_arg(&arg, "FAIL_ONLY=false") {
            ctl.clono_filt_opt.fail_only = false;
        } else if is_simple_arg(&arg, "CHAIN_BRIEF") {
            ctl.clono_print_opt.chain_brief = true;
        } else if is_simple_arg(&arg, "NGRAPH_FILTER") {
            ctl.gen_opt.ngraph_filter = true;
        } else if is_simple_arg(&arg, "INDELS") {
            ctl.gen_opt.indels = true;
        } else if is_simple_arg(&arg, "INSERTIONS") {
            ctl.gen_opt.insertions = true;
        } else if is_simple_arg(&arg, "DEBUG_TABLE_PRINTING") {
            ctl.debug_table_printing = true;
        } else if is_simple_arg(&arg, "KEEP_IMPROPER") {
            ctl.merge_all_impropers = true;
        } else if is_simple_arg(&arg, "NQUAL") {
            ctl.clono_filt_opt.qual_filter = false;
        } else if is_simple_arg(&arg, "HEAVY_CHAIN_REUSE") {
            ctl.gen_opt.heavy_chain_reuse = true;
        } else if is_simple_arg(&arg, "GROUP_HEAVY_CDR3") {
            ctl.clono_group_opt.heavy_cdr3_aa = true;
        } else if is_simple_arg(&arg, "GROUP_VJ_REFNAME") {
            ctl.clono_group_opt.vj_refname = true;
        } else if is_simple_arg(&arg, "NPLAIN") {
            ctl.pretty = true;
        } else if is_simple_arg(&arg, "NO_REUSE") {
            ctl.gen_opt.no_reuse = true;
        } else if is_simple_arg(&arg, "NOPAGER") {
        } else if is_simple_arg(&arg, "NOPRINT") {
            ctl.gen_opt.noprint = true;
        } else if arg.starts_with("POUT=") {
            ctl.parseable_opt.pout = arg.after("POUT=").to_string();
        } else if is_simple_arg(&arg, "PCELL") {
            ctl.parseable_opt.pbarcode = true;
        } else if arg.starts_with("DONOR_REF_FILE=") {
            ctl.gen_opt.dref_file = arg.after("DONOR_REF_FILE=").to_string();
        } else if arg.starts_with("EXT=") {
            ctl.gen_opt.ext = arg.after("EXT=").to_string();
        } else if arg.starts_with("TRACE_BARCODE=") {
            ctl.gen_opt.trace_barcode = arg.after("TRACE_BARCODE=").to_string();
        } else if is_usize_arg(&arg, "PCHAINS") {
            ctl.parseable_opt.pchains = arg.after("PCHAINS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CORES") {
            let nthreads = arg.after("MAX_CORES=").force_usize();
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build_global();
        } else if is_usize_arg(&arg, "REQUIRED_FPS") {
            ctl.gen_opt.required_fps = Some(arg.after("REQUIRED_FPS=").force_usize());
        } else if arg.starts_with("PCOLS=") {
            ctl.parseable_opt.pcols.clear();
            let p = arg.after("PCOLS=").split(',').collect::<Vec<&str>>();
            for i in 0..p.len() {
                let mut x = p[i].to_string();
                x = x.replace("_sum", "_Σ");
                x = x.replace("_mean", "_μ");
                ctl.parseable_opt.pcols.push(x.to_string());
                ctl.parseable_opt.pcols_sort = ctl.parseable_opt.pcols.clone();
                unique_sort(&mut ctl.parseable_opt.pcols_sort);
            }
        } else if is_simple_arg(&arg, "PLAIN") {
        } else if is_simple_arg(&arg, "NOPRETTY") {
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
                if *x == "cdr3" || *x == "var" || *x == "share" || *x == "donor" || *x == "donorn" {
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
            check_cvars(&ctl);
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
            check_cvars(&ctl);
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
        } else if is_simple_arg(&arg, "CDIFF") {
            ctl.clono_filt_opt.cdiff = true;
        } else if is_simple_arg(&arg, "DEL") {
            ctl.clono_filt_opt.del = true;
        } else if is_usize_arg(&arg, "MIN_DATASETS") {
            ctl.clono_filt_opt.min_datasets = arg.after("MIN_DATASETS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_DATASETS") {
            ctl.clono_filt_opt.max_datasets = arg.after("MAX_DATASETS=").force_usize();
        } else if is_usize_arg(&arg, "MIN_GROUP") {
            ctl.clono_group_opt.min_group = arg.after("MIN_GROUP=").force_usize();
        } else if is_simple_arg(&arg, "BCJOIN") {
            ctl.join_alg_opt.bcjoin = true;
        } else if arg.starts_with("EXFASTA=") {
            ctl.gen_opt.fasta = arg.after("EXFASTA=").to_string();
        } else if arg.starts_with("FASTA=") {
            ctl.gen_opt.fasta_filename = arg.after("FASTA=").to_string();
        } else if arg.starts_with("FASTA_AA=") {
            ctl.gen_opt.fasta_aa_filename = arg.after("FASTA_AA=").to_string();
        } else if arg.starts_with("CDR3=") {
            let reg = Regex::new(&format!("^{}$", arg.after("CDR3=")));
            if !reg.is_ok() {
                eprintln!(
                    "\nYour CDR3 value {} could not be parsed as a regular expression.\n",
                    arg.after("CDR3=")
                );
                std::process::exit(1);
            }
            ctl.clono_filt_opt.cdr3 = Some(reg.unwrap());
        } else if arg.starts_with("GEX=") {
        } else if arg.starts_with("BC=") {
        } else if is_usize_arg(&arg, "MIN_MULT") {
            ctl.allele_alg_opt.min_mult = arg.after("MIN_MULT=").force_usize();
        } else if is_usize_arg(&arg, "MIN_EXACTS") {
            ctl.clono_filt_opt.min_exacts = arg.after("MIN_EXACTS=").force_usize();
        } else if is_simple_arg(&arg, "VDUP") {
            ctl.clono_filt_opt.vdup = true;
        } else if is_usize_arg(&arg, "MIN_CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("MIN_CHAINS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CHAINS") {
            ctl.clono_filt_opt.max_chains = arg.after("MAX_CHAINS=").force_usize();
        } else if is_usize_arg(&arg, "CHAINS") {
            ctl.clono_filt_opt.min_chains = arg.after("CHAINS=").force_usize();
            ctl.clono_filt_opt.max_chains = arg.after("CHAINS=").force_usize();
        } else if arg.starts_with("SEG=") {
            let fields = arg.after("SEG=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                ctl.clono_filt_opt.seg.push(x.to_string());
            }
            ctl.clono_filt_opt.seg.sort();
        } else if arg.starts_with("SEGN=") {
            let fields = arg.after("SEGN=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                if !x.parse::<i32>().is_ok() {
                    eprintln!("\nInvalid argument to SEGN.\n");
                    std::process::exit(1);
                }
                ctl.clono_filt_opt.segn.push(x.to_string());
            }
            ctl.clono_filt_opt.segn.sort();
        } else if is_usize_arg(&arg, "MIN_CELLS_EXACT") {
            ctl.gen_opt.min_cells_exact = arg.after("MIN_CELLS_EXACT=").force_usize();
        } else if is_usize_arg(&arg, "MIN_CHAINS_EXACT") {
            ctl.gen_opt.min_chains_exact = arg.after("MIN_CHAINS_EXACT=").force_usize();
        } else if is_usize_arg(&arg, "CHAINS_EXACT") {
            ctl.gen_opt.chains_exact = arg.after("CHAINS_EXACT=").force_usize();
        } else if is_usize_arg(&arg, "EXACT") {
            ctl.gen_opt.exact = Some(arg.after("EXACT=").force_usize());
        } else if is_usize_arg(&arg, "MIN_UMI") {
            ctl.clono_filt_opt.min_umi = arg.after("MIN_UMI=").force_usize();
        } else if is_usize_arg(&arg, "MIN_ALT") {
            ctl.allele_alg_opt.min_alt = arg.after("MIN_ALT=").force_usize();
        } else if is_usize_arg(&arg, "ONESIE_MULT") {
            ctl.onesie_mult = arg.after("ONESIE_MULT=").force_usize();
        } else if arg.starts_with("PRE=") {
        } else if is_usize_arg(&arg, "MIN_CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("MIN_CELLS=").force_usize();
        } else if is_usize_arg(&arg, "MAX_CELLS") {
            ctl.clono_filt_opt.ncells_high = arg.after("MAX_CELLS=").force_usize();
        } else if is_usize_arg(&arg, "CELLS") {
            ctl.clono_filt_opt.ncells_low = arg.after("CELLS=").force_usize();
            ctl.clono_filt_opt.ncells_high = ctl.clono_filt_opt.ncells_low;
        } else if is_simple_arg(&arg, "EASY") {
            ctl.join_alg_opt.easy = true;
        } else if is_usize_arg(&arg, "PFREQ") {
            ctl.join_print_opt.pfreq = arg.after("PFREQ=").force_usize();
        } else if arg.starts_with("HAPS=") {
            // done above
        } else if arg.starts_with("META=") {
            let f = arg.after("META=");
            proc_meta(&f, &mut ctl);
        } else if arg.starts_with("TCR=")
            || arg.starts_with("BCR=")
            || (arg.len() > 0 && arg.as_bytes()[0] >= b'0' && arg.as_bytes()[0] <= b'9')
        {
            proc_xcr(&arg, &gex, &bc, have_gex, &mut ctl);
        } else {
            eprintln!("\nUnrecognized argument {}.\n", arg);
            std::process::exit(1);
        }
    }
    if ctl.parseable_opt.pbarcode && ctl.parseable_opt.pout.len() == 0 {
        eprintln!("\nIt does not make sense to specify PCELL unless POUT is also specified.\n");
        std::process::exit(1);
    }
    if ctl.sample_info.n() == 0 {
        eprintln!("\nNo TCR or BCR data have been specified.\n");
        std::process::exit(1);
    }
    let mut donors = Vec::<String>::new();
    let mut samples = Vec::<String>::new();
    let mut tags = Vec::<String>::new();
    let mut sample_for_bc = Vec::<String>::new();
    let mut donor_for_bc = Vec::<String>::new();
    for i in 0..ctl.sample_info.n() {
        for x in ctl.sample_info.sample_for_bc[i].iter() {
            samples.push(x.1.clone());
            sample_for_bc.push(x.1.clone());
        }
        for x in ctl.sample_info.donor_for_bc[i].iter() {
            donors.push(x.1.clone());
            donor_for_bc.push(x.1.clone());
        }
        for x in ctl.sample_info.tag[i].iter() {
            tags.push((x.1).clone());
        }
        donors.push(ctl.sample_info.donor_id[i].clone());
        samples.push(ctl.sample_info.sample_id[i].clone());
    }
    unique_sort(&mut donors);
    unique_sort(&mut samples);
    unique_sort(&mut tags);
    unique_sort(&mut sample_for_bc);
    unique_sort(&mut donor_for_bc);
    ctl.sample_info.donors = donors.len();
    ctl.sample_info.dataset_list = ctl.sample_info.dataset_id.clone();
    unique_sort(&mut ctl.sample_info.dataset_list);
    ctl.sample_info.sample_list = samples.clone();
    ctl.sample_info.donor_list = donors.clone();
    ctl.sample_info.tag_list = tags;
    for i in 0..ctl.sample_info.donor_for_bc.len() {
        if ctl.sample_info.donor_for_bc[i].len() > 0 {
            ctl.clono_filt_opt.donor = true;
        }
    }
    if ctl.comp2 {
        println!("\n-- used {:.2} seconds processing args", elapsed(&targs));
    }
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
