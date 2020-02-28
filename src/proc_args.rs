// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use crate::proc_args2::*;
use crate::proc_args3::*;
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

    // Mine environment variables and fetch command line args.

    let targs = Instant::now();
    let mut args = args.clone();
    let mut args2 = Vec::<String>::new();
    args2.push(args[0].clone());
    let mut internal_run = false;
    for (key, value) in env::vars() {
        if key.starts_with("ENCLONE_") {
            args2.push(format!("{}={}", key.after("ENCLONE_"), value));
        } else if (key == "HOST" || key == "HOSTNAME") && value.ends_with(".fuzzplex.com") {
            internal_run = true;
            ctl.gen_opt.pre = "/mnt/assembly/vdj/current14".to_string();
        }
    }
    for i in 1..args.len() {
        args2.push(args[i].clone());
    }
    args = args2;

    // Set up general options.

    ctl.gen_opt.min_cells_exact = 1;
    ctl.gen_opt.min_chains_exact = 1;
    ctl.gen_opt.exact = None;
    for i in 1..args.len() {
        if args[i].starts_with("PRE=") {
            ctl.gen_opt.pre = args[i].after("PRE=").to_string();
        }
    }
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

    ctl.clono_print_opt.amino = vec![
        "cdr3".to_string(),
        "var".to_string(),
        "share".to_string(),
        "donor".to_string(),
    ];
    ctl.clono_print_opt.cvars = vec!["umed".to_string(), "const".to_string(), "notes".to_string()];
    ctl.clono_print_opt.lvars = vec!["datasets".to_string(), "ncells".to_string()];

    ctl.clono_group_opt.min_group = 1;

    ctl.allele_alg_opt.min_mult = 4;
    ctl.allele_alg_opt.min_alt = 4;

    ctl.join_alg_opt.max_score = 1_000_000.0;
    ctl.join_alg_opt.merge_onesies = true; // should just kill this as an option

    ctl.join_print_opt.pfreq = 1_000_000_000;
    ctl.join_print_opt.quiet = true;

    ctl.parseable_opt.pchains = 4;

    ctl.onesie_mult = 10_000;

    let cvars_allowed = vec![
        "var", "umed", "umax", "comp", "utot", "rmed", "const", "white", "cdr3_dna", "ulen",
        "clen", "cdiff", "udiff", "notes", "d_univ", "d_donor",
    ];

    // Pretest for consistency amongst TCR, BCR, GEX and META.  Also preparse GEX.

    let mut have_tcr = false;
    let mut have_bcr = false;
    let mut have_gex = false;
    let mut have_meta = false;
    let mut gex = String::new();
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
    }
    if have_meta && (have_tcr || have_bcr || have_gex) {
        eprintln!("\nIf META is specified, then none of TCR, BCR or GEX can be specified.\n");
        std::process::exit(1);
    }
    if have_tcr && have_bcr {
        eprintln!("\nPlease do not specify both TCR and BCR.\n");
        std::process::exit(1);
    }

    // Traverse arguments.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "SEQ") {
            ctl.join_print_opt.seq = true;
        } else if is_simple_arg(&args[i], "ANN") {
            ctl.join_print_opt.ann = true;
        } else if is_simple_arg(&args[i], "ANN0") {
            ctl.join_print_opt.ann0 = true;
        } else if is_simple_arg(&args[i], "DUMP_LENAS") {
        } else if is_simple_arg(&args[i], "BC") {
            ctl.join_print_opt.show_bc = true;
        } else if is_simple_arg(&args[i], "PER_BC") {
            ctl.clono_print_opt.bu = true;
        } else if is_simple_arg(&args[i], "COMP") {
        } else if is_simple_arg(&args[i], "CON") {
            ctl.allele_print_opt.con = true;
        } else if is_simple_arg(&args[i], "CON_TRACE") {
            ctl.allele_print_opt.con_trace = true;
        } else if is_simple_arg(&args[i], "EXP") {
            ctl.gen_opt.exp = true;
        } else if is_simple_arg(&args[i], "CURRENT_REF") {
            ctl.gen_opt.current_ref = true;
        } else if is_simple_arg(&args[i], "SUM") {
            ctl.clono_print_opt.sum = true;
        } else if is_simple_arg(&args[i], "MEAN") {
            ctl.clono_print_opt.mean = true;
        } else if is_simple_arg(&args[i], "NH5") {
        } else if is_simple_arg(&args[i], "H5_PRE") {
            ctl.gen_opt.h5_pre = true;
        } else if is_simple_arg(&args[i], "DESCRIP") {
            ctl.gen_opt.descrip = true;
        } else if is_simple_arg(&args[i], "CTRLC") {
        } else if is_simple_arg(&args[i], "FORCE") {
            ctl.force = true;
        } else if is_simple_arg(&args[i], "CELLRANGER") {
        } else if is_simple_arg(&args[i], "WEAK") {
            ctl.gen_opt.weak = true;
        } else if is_simple_arg(&args[i], "REUSE") {
            ctl.gen_opt.reuse = true;
        } else if is_simple_arg(&args[i], "NWARN") {
            ctl.gen_opt.nwarn = true;
        } else if args[i].starts_with("BINARY=") {
            ctl.gen_opt.binary = args[i].after("BINARY=").to_string();
        } else if args[i].starts_with("PROTO=") {
            ctl.gen_opt.proto = args[i].after("PROTO=").to_string();
        } else if is_simple_arg(&args[i], "PRINT_FAILED_JOINS") {
            ctl.join_print_opt.quiet = false;
        } else if is_simple_arg(&args[i], "NOTE_SIMPLE") {
            ctl.clono_print_opt.note_simple = true;
        } else if is_simple_arg(&args[i], "SEQC") {
            ctl.clono_print_opt.seqc = true;
        } else if is_simple_arg(&args[i], "FULL_SEQC") {
            ctl.clono_print_opt.full_seqc = true;
        } else if is_simple_arg(&args[i], "BARCODES") {
            ctl.clono_print_opt.barcodes = true;
        } else if args[i].starts_with("BARCODE=") {
            let bcs = args[i].after("BARCODE=").split(',').collect::<Vec<&str>>();
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
        } else if is_simple_arg(&args[i], "GRAPH") {
            ctl.gen_opt.graph = true;
        } else if is_simple_arg(&args[i], "ACCEPT_INCONSISTENT") {
            ctl.gen_opt.accept_inconsistent = true;
        } else if is_simple_arg(&args[i], "NCROSS") {
            ctl.clono_filt_opt.ncross = true;
        } else if is_simple_arg(&args[i], "NWEAK_CHAINS") {
            ctl.clono_filt_opt.weak_chains = false;
        } else if is_simple_arg(&args[i], "NWEAK_ONESIES") {
            ctl.clono_filt_opt.weak_onesies = false;
        } else if is_simple_arg(&args[i], "NFOURSIE_KILL") {
            ctl.clono_filt_opt.weak_foursies = false;
        } else if is_simple_arg(&args[i], "NBC_DUP") {
            ctl.clono_filt_opt.bc_dup = false;
        } else if is_simple_arg(&args[i], "NDONOR") {
            ctl.clono_filt_opt.donor = true;
        } else if is_simple_arg(&args[i], "HAVE_ONESIE") {
            ctl.clono_filt_opt.have_onesie = true;
        } else if is_simple_arg(&args[i], "UTR_CON") {
            ctl.gen_opt.utr_con = true;
        } else if is_simple_arg(&args[i], "CON_CON") {
            ctl.gen_opt.con_con = true;
        } else if is_simple_arg(&args[i], "MOUSE") {
            ctl.gen_opt.mouse = true;
        } else if is_simple_arg(&args[i], "SUMMARY") {
            ctl.gen_opt.summary = true;
        } else if args[i].starts_with("F=") {
            let filt = args[i].after("F=").to_string();
            ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt));
        } else if args[i].starts_with("SCAN=") {
            let mut x = args[i].after("SCAN=").to_string();
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
        } else if args[i].starts_with("PLOT=") {
            let x = args[i].after("PLOT=").split(',').collect::<Vec<&str>>();
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
        } else if is_simple_arg(&args[i], "SUMMARY_CLEAN") {
            ctl.gen_opt.summary_clean = true;
        } else if args[i].starts_with("EMAIL=") {
        } else if args[i].starts_with("REF=") {
            ctl.gen_opt.refname = args[i].after("REF=").to_string();
        } else if is_simple_arg(&args[i], "NSILENT") {
            ctl.silent = false;
        } else if is_simple_arg(&args[i], "TOY") {
            ctl.toy = true;
        } else if is_simple_arg(&args[i], "RE") {
            ctl.gen_opt.reannotate = true;
        } else if is_simple_arg(&args[i], "WHITEF") {
            ctl.clono_filt_opt.whitef = true;
        } else if is_simple_arg(&args[i], "PROTECT_BADS") {
            ctl.clono_filt_opt.protect_bads = true;
        } else if is_simple_arg(&args[i], "NWHITEF") {
            ctl.gen_opt.nwhitef = true;
        } else if is_simple_arg(&args[i], "FAIL_ONLY=true") {
            ctl.clono_filt_opt.fail_only = true;
        } else if is_simple_arg(&args[i], "FAIL_ONLY=false") {
            ctl.clono_filt_opt.fail_only = false;
        } else if is_simple_arg(&args[i], "CHAIN_BRIEF") {
            ctl.clono_print_opt.chain_brief = true;
        } else if is_simple_arg(&args[i], "NGRAPH_FILTER") {
            ctl.gen_opt.ngraph_filter = true;
        } else if is_simple_arg(&args[i], "INDELS") {
            ctl.gen_opt.indels = true;
        } else if is_simple_arg(&args[i], "INSERTIONS") {
            ctl.gen_opt.insertions = true;
        } else if is_simple_arg(&args[i], "DEBUG_TABLE_PRINTING") {
            ctl.debug_table_printing = true;
        } else if is_simple_arg(&args[i], "KEEP_IMPROPER") {
            ctl.merge_all_impropers = true;
        } else if is_simple_arg(&args[i], "NQUAL") {
            ctl.clono_filt_opt.qual_filter = false;
        } else if is_simple_arg(&args[i], "HEAVY_CHAIN_REUSE") {
            ctl.gen_opt.heavy_chain_reuse = true;
        } else if is_simple_arg(&args[i], "GROUP_HEAVY_CDR3") {
            ctl.clono_group_opt.heavy_cdr3_aa = true;
        } else if is_simple_arg(&args[i], "GROUP_VJ_REFNAME") {
            ctl.clono_group_opt.vj_refname = true;
        } else if is_simple_arg(&args[i], "NPLAIN") {
            ctl.pretty = true;
        } else if is_simple_arg(&args[i], "NO_REUSE") {
            ctl.gen_opt.no_reuse = true;
        } else if is_simple_arg(&args[i], "NOPAGER") {
        } else if is_simple_arg(&args[i], "NOPRINT") {
            ctl.gen_opt.noprint = true;
        } else if args[i].starts_with("POUT=") {
            ctl.parseable_opt.pout = args[i].after("POUT=").to_string();
        } else if args[i].starts_with("DONOR_REF_FILE=") {
            ctl.gen_opt.dref_file = args[i].after("DONOR_REF_FILE=").to_string();
        } else if args[i].starts_with("EXT=") {
            ctl.gen_opt.ext = args[i].after("EXT=").to_string();
        } else if is_usize_arg(&args[i], "PCHAINS") {
            ctl.parseable_opt.pchains = args[i].after("PCHAINS=").force_usize();
        } else if is_usize_arg(&args[i], "MAX_THREADS") {
            let nthreads = args[i].after("MAX_THREADS=").force_usize();
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build_global();
        } else if is_usize_arg(&args[i], "REQUIRED_FPS") {
            ctl.gen_opt.required_fps = Some(args[i].after("REQUIRED_FPS=").force_usize());
        } else if args[i].starts_with("PCOLS=") {
            ctl.parseable_opt.pcols.clear();
            for x in args[i].after("PCOLS=").split(',').collect::<Vec<&str>>() {
                ctl.parseable_opt.pcols.push(x.to_string());
                ctl.parseable_opt.pcols_sort = ctl.parseable_opt.pcols.clone();
                unique_sort(&mut ctl.parseable_opt.pcols_sort);
            }
        } else if is_simple_arg(&args[i], "PLAIN") {
        } else if is_simple_arg(&args[i], "NOPRETTY") {
        } else if args[i].starts_with("VJ=") {
            ctl.clono_filt_opt.vj = args[i].after("VJ=").as_bytes().to_vec();
            for c in ctl.clono_filt_opt.vj.iter() {
                if !(*c == b'A' || *c == b'C' || *c == b'G' || *c == b'T') {
                    eprintln!("\nIllegal value for VJ, must be over alphabet ACGT.\n");
                    std::process::exit(1);
                }
            }
        } else if args[i].starts_with("AMINO=") {
            ctl.clono_print_opt.amino.clear();
            for x in args[i].after("AMINO=").split(',').collect::<Vec<&str>>() {
                if x != "" {
                    ctl.clono_print_opt.amino.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.amino.iter() {
                if !(*x == "cdr3"
                    || *x == "var"
                    || *x == "share"
                    || *x == "donor"
                    || *x == "donorn")
                {
                    eprintln!(
                        "\nUnrecognized variable {} for AMINO.  Please type \
                         \"enclone help amino\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
        } else if args[i].starts_with("CVARS=") {
            ctl.clono_print_opt.cvars.clear();
            for x in args[i].after("CVARS=").split(',').collect::<Vec<&str>>() {
                if x.len() > 0 {
                    ctl.clono_print_opt.cvars.push(x.to_string());
                }
            }
            for x in ctl.clono_print_opt.cvars.iter() {
                let mut ok = cvars_allowed.contains(&(*x).as_str());
                if x.starts_with("ndiff")
                    && x.after("ndiff").parse::<usize>().is_ok()
                    && x.after("ndiff").force_usize() >= 1
                {
                    ok = true;
                }
                if !ok {
                    eprintln!(
                        "\nUnrecognized variable {} for CVARS.  Please type \
                         \"enclone help cvars\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
        } else if args[i].starts_with("CVARSP=") {
            let cvarsp = args[i].after("CVARSP=").split(',').collect::<Vec<&str>>();
            for x in cvarsp.iter() {
                let mut ok = cvars_allowed.contains(&x);
                if x.starts_with("ndiff")
                    && x.after("ndiff").parse::<usize>().is_ok()
                    && x.after("ndiff").force_usize() >= 1
                {
                    ok = true;
                }
                if !ok {
                    eprintln!(
                        "\nUnrecognized variable {} for CVARSP.  Please type \
                         \"enclone help cvars\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
            for x in cvarsp {
                ctl.clono_print_opt.cvars.push(x.to_string());
            }
        } else if args[i].starts_with("LVARS=") {
            ctl.clono_print_opt.lvars.clear();
            for x in args[i].after("LVARS=").split(',').collect::<Vec<&str>>() {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
        } else if args[i].starts_with("LVARSP=") {
            let lvarsp = args[i].after("LVARSP=").split(',').collect::<Vec<&str>>();
            for x in lvarsp {
                ctl.clono_print_opt.lvars.push(x.to_string());
            }
        } else if is_f64_arg(&args[i], "MAX_SCORE") {
            ctl.join_alg_opt.max_score = args[i].after("MAX_SCORE=").force_f64();
        } else if is_simple_arg(&args[i], "CDIFF") {
            ctl.clono_filt_opt.cdiff = true;
        } else if is_simple_arg(&args[i], "DEL") {
            ctl.clono_filt_opt.del = true;
        } else if is_usize_arg(&args[i], "MIN_DATASETS") {
            ctl.clono_filt_opt.min_datasets = args[i].after("MIN_DATASETS=").force_usize();
        } else if is_usize_arg(&args[i], "MIN_GROUP") {
            ctl.clono_group_opt.min_group = args[i].after("MIN_GROUP=").force_usize();
        } else if is_simple_arg(&args[i], "BCJOIN") {
            ctl.join_alg_opt.bcjoin = true;
        } else if args[i].starts_with("EXFASTA=") {
            ctl.gen_opt.fasta = args[i].after("EXFASTA=").to_string();
        } else if args[i].starts_with("FASTA=") {
            ctl.gen_opt.fasta_filename = args[i].after("FASTA=").to_string();
        } else if args[i].starts_with("FASTA_AA=") {
            ctl.gen_opt.fasta_aa_filename = args[i].after("FASTA_AA=").to_string();
        } else if args[i].starts_with("CDR3=") {
            let reg = Regex::new(&format!("^{}$", args[i].after("CDR3=")));
            if !reg.is_ok() {
                eprintln!(
                    "\nYour CDR3 value {} could not be parsed as a regular expression.\n",
                    args[i].after("CDR3=")
                );
                std::process::exit(1);
            }
            ctl.clono_filt_opt.cdr3 = Some(reg.unwrap());
        } else if args[i].starts_with("GEX=") {
        } else if is_usize_arg(&args[i], "MIN_MULT") {
            ctl.allele_alg_opt.min_mult = args[i].after("MIN_MULT=").force_usize();
        } else if is_usize_arg(&args[i], "MIN_EXACTS") {
            ctl.clono_filt_opt.min_exacts = args[i].after("MIN_EXACTS=").force_usize();
        } else if is_simple_arg(&args[i], "VDUP") {
            ctl.clono_filt_opt.vdup = true;
        } else if is_usize_arg(&args[i], "MIN_CHAINS") {
            ctl.clono_filt_opt.min_chains = args[i].after("MIN_CHAINS=").force_usize();
        } else if is_usize_arg(&args[i], "MAX_CHAINS") {
            ctl.clono_filt_opt.max_chains = args[i].after("MAX_CHAINS=").force_usize();
        } else if is_usize_arg(&args[i], "CHAINS") {
            ctl.clono_filt_opt.min_chains = args[i].after("CHAINS=").force_usize();
            ctl.clono_filt_opt.max_chains = args[i].after("CHAINS=").force_usize();
        } else if args[i].starts_with("SEG=") {
            let fields = args[i].after("SEG=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                ctl.clono_filt_opt.seg.push(x.to_string());
            }
            ctl.clono_filt_opt.seg.sort();
        } else if args[i].starts_with("SEGN=") {
            let fields = args[i].after("SEGN=").split('|').collect::<Vec<&str>>();
            for x in fields.iter() {
                if !x.parse::<i32>().is_ok() {
                    eprintln!("\nInvalid argument to SEGN.\n");
                    std::process::exit(1);
                }
                ctl.clono_filt_opt.segn.push(x.to_string());
            }
            ctl.clono_filt_opt.segn.sort();
        } else if is_usize_arg(&args[i], "MIN_CELLS_EXACT") {
            ctl.gen_opt.min_cells_exact = args[i].after("MIN_CELLS_EXACT=").force_usize();
        } else if is_usize_arg(&args[i], "MIN_CHAINS_EXACT") {
            ctl.gen_opt.min_chains_exact = args[i].after("MIN_CHAINS_EXACT=").force_usize();
        } else if is_usize_arg(&args[i], "EXACT") {
            ctl.gen_opt.exact = Some(args[i].after("EXACT=").force_usize());
        } else if is_usize_arg(&args[i], "MIN_UMI") {
            ctl.clono_filt_opt.min_umi = args[i].after("MIN_UMI=").force_usize();
        } else if is_usize_arg(&args[i], "MIN_ALT") {
            ctl.allele_alg_opt.min_alt = args[i].after("MIN_ALT=").force_usize();
        } else if is_usize_arg(&args[i], "ONESIE_MULT") {
            ctl.onesie_mult = args[i].after("ONESIE_MULT=").force_usize();
        } else if args[i].starts_with("PRE=") {
        } else if args[i] == "NOPAR" {
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build_global();
        } else if is_usize_arg(&args[i], "MIN_CELLS") {
            ctl.clono_filt_opt.ncells_low = args[i].after("MIN_CELLS=").force_usize();
        } else if is_usize_arg(&args[i], "MAX_CELLS") {
            ctl.clono_filt_opt.ncells_high = args[i].after("MAX_CELLS=").force_usize();
        } else if is_usize_arg(&args[i], "CELLS") {
            ctl.clono_filt_opt.ncells_low = args[i].after("CELLS=").force_usize();
            ctl.clono_filt_opt.ncells_high = ctl.clono_filt_opt.ncells_low;
        } else if is_simple_arg(&args[i], "EASY") {
            ctl.join_alg_opt.easy = true;
        } else if is_usize_arg(&args[i], "PFREQ") {
            ctl.join_print_opt.pfreq = args[i].after("PFREQ=").force_usize();
        } else if args[i].starts_with("HAPS=") {
            // done above
        } else if args[i].starts_with("META=") {
            let f = args[i].after("META=");
            proc_meta(&f, &mut ctl);
        } else if args[i].starts_with("TCR=")
            || args[i].starts_with("BCR=")
            || (args[i].len() > 0 && args[i].as_bytes()[0] >= b'0' && args[i].as_bytes()[0] <= b'9')
        {
            proc_xcr(&args[i], &gex, have_gex, internal_run, &mut ctl);
        } else {
            eprintln!("\nUnrecognized argument {}.\n", args[i]);
            std::process::exit(1);
        }
    }
    if ctl.sample_info.n() == 0 {
        eprintln!("\nNo TCR or BCR data have been specified.\n");
        std::process::exit(1);
    }
    let mut donors = Vec::<String>::new();
    let mut samples = Vec::<String>::new();
    let mut tags = Vec::<String>::new();
    let mut sample_donor = Vec::<(String, String)>::new();
    for i in 0..ctl.sample_info.n() {
        for x in ctl.sample_info.sample_donor[i].iter() {
            donors.push((x.1).1.clone());
            samples.push((x.1).0.clone());
            sample_donor.push(((x.1).0.clone(), (x.1).1.clone()));
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
    unique_sort(&mut sample_donor);
    ctl.sample_info.donors = donors.len();
    ctl.sample_info.donor_list = donors.clone();
    ctl.sample_info.sample_list = samples.clone();
    ctl.sample_info.tag_list = tags;
    let mut sample_donor_list = Vec::<(usize, usize)>::new();
    for i in 0..sample_donor.len() {
        sample_donor_list.push((
            bin_position(&samples, &sample_donor[i].0) as usize,
            bin_position(&donors, &sample_donor[i].1) as usize,
        ));
    }
    unique_sort(&mut sample_donor_list);
    ctl.sample_info.sample_donor_list = sample_donor_list;
    if ctl.comp {
        println!("-- used {:.2} seconds processing args", elapsed(&targs));
    }
    proc_args_tail(&mut ctl, &args, internal_run);

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
