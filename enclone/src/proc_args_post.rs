// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::proc_args2::*;
use crate::proc_args3::*;
use crate::proc_args_check::*;
use enclone_core::defs::*;
use io_utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use string_utils::*;
use tilde_expand::*;
use vector_utils::*;

pub fn proc_args_post(
    mut ctl: &mut EncloneControl,
    args: &Vec<String>,
    metas: &Vec<String>,
    metaxs: &Vec<String>,
    xcrs: &Vec<String>,
    have_gex: bool,
    gex: &String,
    bc: &String,
    using_plot: bool,
) {
    // Process INFO.

    let t = Instant::now();
    if ctl.gen_opt.info.is_some() {
        let f = open_for_read![&ctl.gen_opt.info.as_ref().unwrap()];
        let mut lines = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
        if lines.is_empty() {
            eprintln!(
                "\nThe file {} is empty.\n",
                ctl.gen_opt.info.as_ref().unwrap()
            );
            std::process::exit(1);
        }
        let fields = lines[0].split(',').collect::<Vec<&str>>();
        if !fields.contains(&"vj_seq1") || !fields.contains(&"vj_seq2") {
            eprintln!(
                "\nThe CSV file {} needs to have fields vj_seq1 and vj_seq2.\n",
                ctl.gen_opt.info.as_ref().unwrap()
            );
            std::process::exit(1);
        }
        for i in 0..fields.len() {
            if fields[i] != "vj_seq1" && fields[i] != "vj_seq2" {
                ctl.gen_opt.info_fields.push(fields[i].to_string());
            }
        }
        let mut tags = Vec::<String>::new();
        for i in 1..lines.len() {
            let vals = parse_csv(&lines[i]);
            if vals.len() != fields.len() {
                eprintln!(
                    "\nINFO file line {} has length {} whereas the file has {} fields. \
                    The line is\n{}\n",
                    i + 1,
                    vals.len(),
                    fields.len(),
                    lines[i]
                );
            }
            let (mut vj1, mut vj2) = (String::new(), String::new());
            let mut other = Vec::<String>::new();
            for i in 0..vals.len() {
                if fields[i] == "vj_seq1" {
                    vj1 = vals[i].to_string();
                } else if fields[i] == "vj_seq2" {
                    vj2 = vals[i].to_string();
                } else {
                    other.push(vals[i].to_string());
                }
            }
            let tag = format!("{}_{}", vj1, vj2);
            tags.push(tag.clone());
            ctl.gen_opt.info_data.insert(tag, other);
        }
        tags.sort();
        let mut i = 0;
        while i < tags.len() {
            let j = next_diff(&tags, i);
            if j - i > 1 {
                eprintln!(
                    "\nThe immune receptor sequence pair\n{},\n {}\nappears more than once \
                    in the file {}.\n",
                    tags[i].before("_"),
                    tags[i].after("_"),
                    ctl.gen_opt.info.as_ref().unwrap(),
                );
                std::process::exit(1);
            }
            i = j;
        }
    }

    // Expand ~ and ~user in output file names.

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

    // Sanity check grouping arguments.

    let mut group_styles = 0;
    if ctl.clono_group_opt.heavy_cdr3_aa {
        group_styles += 1;
    }
    if ctl.clono_group_opt.vj_refname {
        group_styles += 1;
    }
    if ctl.clono_group_opt.vj_refname_strong {
        group_styles += 1;
    }
    if ctl.clono_group_opt.asymmetric {
        group_styles += 1;
    }
    if group_styles > 1 {
        eprintln!(
            "\nOnly one of the options\n\
            GROUP_HEAVY_CDR3, GROUP_VJ_REFNAME, GROUP_VJ_REFNAME_STRONG, AGROUP\n\
            may be specified at a time.\n"
        );
        std::process::exit(1);
    }
    if ctl.clono_group_opt.asymmetric {
        if ctl.clono_group_opt.asymmetric_center.len() == 0
            || ctl.clono_group_opt.asymmetric_dist_formula.len() == 0
            || ctl.clono_group_opt.asymmetric_dist_bound.len() == 0
        {
            eprintln!(
                "\nIf the AGROUP option is used to specify asymmetric grouping, then all\n\
                of the options AG_CENTER, AG_DIST_FORMULA and AG_DIST_BOUND must also be \
                specified.\n"
            );
            std::process::exit(1);
        }
    }
    if ctl.clono_group_opt.asymmetric_center.len() > 0
        || ctl.clono_group_opt.asymmetric_dist_formula.len() > 0
        || ctl.clono_group_opt.asymmetric_dist_bound.len() > 0
    {
        if !ctl.clono_group_opt.asymmetric {
            eprintln!(
                "\nIf any of the asymmetric grouping options AG_CENTER or \
                    AG_DIST_FORMULA or\nAG_DIST_BOUND are specified, then the option AGROUP \
                    must also be specified, to turn on asymmetric grouping.\n"
            );
            std::process::exit(1);
        }
    }
    if ctl.clono_group_opt.asymmetric {
        if ctl.clono_group_opt.asymmetric_center != "from_filters"
            && ctl.clono_group_opt.asymmetric_center != "copy_filters"
        {
            eprintln!(
                "\nThe only allowed forms for AG_CENTER are AG_CENTER=from_filters\n\
                and AG_CENTER=copy_filters.\n"
            );
            std::process::exit(1);
        }
        if ctl.clono_group_opt.asymmetric_dist_formula != "cdr3_edit_distance" {
            eprintln!("\nThe only allowed form for AG_DIST_FORMULA is cdr3_edit_distance.\n");
            std::process::exit(1);
        }
        let ok1 = ctl
            .clono_group_opt
            .asymmetric_dist_bound
            .starts_with("top=")
            && ctl
                .clono_group_opt
                .asymmetric_dist_bound
                .after("top=")
                .parse::<usize>()
                .is_ok();
        let ok2 = ctl
            .clono_group_opt
            .asymmetric_dist_bound
            .starts_with("max=")
            && ctl
                .clono_group_opt
                .asymmetric_dist_bound
                .after("max=")
                .parse::<f64>()
                .is_ok();
        if !ok1 && !ok2 {
            eprintln!(
                "\nThe only allowed forms for AG_DIST_BOUND are top=n, where n is an\n\
                integer, and max=d, where d is a number.\n"
            );
            std::process::exit(1);
        }
    }

    // Sanity check other arguments (and more below).

    if ctl.clono_print_opt.conx && ctl.clono_print_opt.conp {
        eprintln!("\nPlease specify at most one of CONX and CONP.\n");
        std::process::exit(1);
    }
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
    if metaxs.len() > 0 {
        let lines0 = metaxs[metaxs.len() - 1].split(';').collect::<Vec<&str>>();
        let mut lines = Vec::<String>::new();
        for i in 0..lines0.len() {
            lines.push(lines0[i].to_string());
        }
        proc_meta_core(&lines, &mut ctl);
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
                let f = open_userfile_for_read(&f);
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
