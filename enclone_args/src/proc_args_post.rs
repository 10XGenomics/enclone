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
) -> Result<(), String> {
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
            return Err(format!(
                "\nThe file {} is empty.\n",
                ctl.gen_opt.info.as_ref().unwrap()
            ));
        }
        let fields = lines[0].split(',').collect::<Vec<&str>>();
        if !fields.contains(&"vj_seq1") || !fields.contains(&"vj_seq2") {
            return Err(format!(
                "\nThe CSV file {} needs to have fields vj_seq1 and vj_seq2.\n",
                ctl.gen_opt.info.as_ref().unwrap()
            ));
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
            if ctl.gen_opt.info_resolve && ctl.gen_opt.info_data.contains_key(&tag) {
                continue;
            }
            tags.push(tag.clone());
            ctl.gen_opt.info_data.insert(tag, other);
        }
        tags.sort();
        let mut i = 0;
        while i < tags.len() {
            let j = next_diff(&tags, i);
            if j - i > 1 {
                return Err(format!(
                    "\nThe immune receptor sequence pair\n{},\n {}\nappears more than once \
                    in the file {}.\n",
                    tags[i].before("_"),
                    tags[i].after("_"),
                    ctl.gen_opt.info.as_ref().unwrap(),
                ));
            }
            i = j;
        }
    }

    // Expand ~ and ~user in output file names.

    let mut files = [
        &mut ctl.plot_opt.plot_file,
        &mut ctl.gen_opt.fasta_filename,
        &mut ctl.gen_opt.fasta_aa_filename,
        &mut ctl.gen_opt.dref_file,
        &mut ctl.parseable_opt.pout,
    ];
    for f in files.iter_mut() {
        **f = stringme(&tilde_expand(&f.as_bytes()));
    }

    // Sanity check grouping arguments.

    if ctl.clono_group_opt.style == "asymmetric" {
        if ctl.clono_group_opt.asymmetric_center.len() == 0
            || ctl.clono_group_opt.asymmetric_dist_formula.len() == 0
            || ctl.clono_group_opt.asymmetric_dist_bound.len() == 0
        {
            return Err(format!(
                "\nIf the AGROUP option is used to specify asymmetric grouping, then all\n\
                of the options AG_CENTER, AG_DIST_FORMULA and AG_DIST_BOUND must also be \
                specified.\n"
            ));
        }
    }
    if ctl.clono_group_opt.asymmetric_center.len() > 0
        || ctl.clono_group_opt.asymmetric_dist_formula.len() > 0
        || ctl.clono_group_opt.asymmetric_dist_bound.len() > 0
    {
        if ctl.clono_group_opt.style == "symmetric" {
            return Err(format!(
                "\nIf any of the asymmetric grouping options AG_CENTER or \
                    AG_DIST_FORMULA or\nAG_DIST_BOUND are specified, then the option AGROUP \
                    must also be specified, to turn on asymmetric grouping.\n"
            ));
        }
    }
    if ctl.clono_group_opt.style == "asymmetric" {
        if ctl.clono_group_opt.asymmetric_center != "from_filters"
            && ctl.clono_group_opt.asymmetric_center != "copy_filters"
        {
            return Err(format!(
                "\nThe only allowed forms for AG_CENTER are AG_CENTER=from_filters\n\
                and AG_CENTER=copy_filters.\n"
            ));
        }
        if ctl.clono_group_opt.asymmetric_dist_formula != "cdr3_edit_distance" {
            return Err(format!(
                "\nThe only allowed form for AG_DIST_FORMULA is cdr3_edit_distance.\n"
            ));
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
            return Err(format!(
                "\nThe only allowed forms for AG_DIST_BOUND are top=n, where n is an\n\
                integer, and max=d, where d is a number.\n"
            ));
        }
    }

    // Sanity check other arguments (and more below).

    if !ctl.parseable_opt.pchains.parse::<usize>().is_ok() && ctl.parseable_opt.pchains != "max" {
        return Err(format!(
            "\nThe only allowed values for PCHAINS are a positive integer and max.\n"
        ));
    }
    if ctl.gen_opt.align_jun_align_consistency && ctl.pretty {
        return Err(format!(
            "\nIf you use ALIGN_JALIGN_CONSISTENCY, you should also use PLAIN.\n"
        ));
    }
    if ctl.gen_opt.gene_scan_exact && ctl.gen_opt.gene_scan_test.is_none() {
        return Err(format!(
            "\nIt doesn't make sense to specify SCAN_EXIT unless SCAN is also specified.\n"
        ));
    }
    if ctl.clono_print_opt.conx && ctl.clono_print_opt.conp {
        return Err(format!("\nPlease specify at most one of CONX and CONP.\n"));
    }
    if ctl.clono_filt_opt.cdr3.is_some() && ctl.clono_filt_opt.cdr3_lev.len() > 0 {
        return Err(format!(
            "\nPlease use the CDR3 argument to specify either a regular expression or a\n\
            Levenshtein distance pattern, but not both.\n"
        ));
    }
    if ctl.gen_opt.clustal_aa != "".to_string() && ctl.gen_opt.clustal_aa != "stdout".to_string() {
        if !ctl.gen_opt.clustal_aa.ends_with(".tar") {
            return Err(format!(
                "\nIf the value of CLUSTAL_AA is not stdout, it must end in .tar.\n"
            ));
        }
    }
    if ctl.gen_opt.clustal_dna != "".to_string() && ctl.gen_opt.clustal_dna != "stdout".to_string()
    {
        if !ctl.gen_opt.clustal_dna.ends_with(".tar") {
            return Err(format!(
                "\nIf the value of CLUSTAL_DNA is not stdout, it must end in .tar.\n"
            ));
        }
    }
    if ctl.gen_opt.phylip_aa != "".to_string() && ctl.gen_opt.phylip_aa != "stdout".to_string() {
        if !ctl.gen_opt.phylip_aa.ends_with(".tar") {
            return Err(format!(
                "\nIf the value of PHYLIP_AA is not stdout, it must end in .tar.\n"
            ));
        }
    }
    if ctl.gen_opt.phylip_dna != "".to_string() && ctl.gen_opt.phylip_dna != "stdout".to_string() {
        if !ctl.gen_opt.phylip_dna.ends_with(".tar") {
            return Err(format!(
                "\nIf the value of PHYLIP_DNA is not stdout, it must end in .tar.\n"
            ));
        }
    }
    if ctl.clono_filt_opt_def.umi_filt && ctl.clono_filt_opt_def.umi_filt_mark {
        return Err(format!(
            "\nIf you use UMI_FILT_MARK, you should also use NUMI, to turn off \
            the filter,\nas otherwise nothing will be marked.\n"
        ));
    }
    if ctl.clono_filt_opt_def.umi_ratio_filt && ctl.clono_filt_opt_def.umi_ratio_filt_mark {
        return Err(format!(
            "\nIf you use UMI_RATIO_FILT_MARK, you should also use NUMI_RATIO, to turn off \
            the filter,\nas otherwise nothing will be marked.\n"
        ));
    }
    ctl.perf_stats(&t, "after main args loop 1");

    // Process TCR, BCR and META.

    let t = Instant::now();
    check_cvars(&ctl)?;
    if metas.len() > 0 {
        let f = &metas[metas.len() - 1];
        let f = get_path_fail(&f, &ctl, "META")?;
        proc_meta(&f, &mut ctl)?;
    }
    if metaxs.len() > 0 {
        let lines0 = metaxs[metaxs.len() - 1].split(';').collect::<Vec<&str>>();
        let mut lines = Vec::<String>::new();
        for i in 0..lines0.len() {
            lines.push(lines0[i].to_string());
        }
        proc_meta_core(&lines, &mut ctl)?;
    }
    ctl.perf_stats(&t, "in proc_meta");
    if xcrs.len() > 0 {
        let arg = &xcrs[xcrs.len() - 1];
        proc_xcr(&arg, &gex, &bc, have_gex, &mut ctl)?;
    }

    // More argument sanity checking.

    let t = Instant::now();
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
                    return Err(format!("\nThe option {} does not make sense for TCR.\n", x));
                }
            }
        }
    }

    // Proceed.

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
    if ctl.plot_opt.plot_by_isotype {
        if using_plot || ctl.plot_opt.use_legend {
            return Err(format!(
                "\nPLOT_BY_ISOTYPE cannot be used with PLOT or LEGEND.\n"
            ));
        }
        if !ctl.gen_opt.bcr {
            return Err(format!(
                "\nPLOT_BY_ISOTYPE can only be used with BCR data.\n"
            ));
        }
        if ctl.plot_opt.plot_by_mark {
            return Err(format!(
                "\nPLOT_BY_ISOTYPE and PLOT_BY_MARK cannot be used together.\n"
            ));
        }
    }
    if ctl.plot_opt.plot_by_mark {
        if using_plot || ctl.plot_opt.use_legend {
            return Err(format!(
                "\nPLOT_BY_MARK cannot be used with PLOT or LEGEND.\n"
            ));
        }
    }
    if ctl.parseable_opt.pbarcode && ctl.parseable_opt.pout.len() == 0 {
        return Err(format!(
            "\nIt does not make sense to specify PCELL unless POUT is also specified.\n"
        ));
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
            ctl.clono_filt_opt_def.donor = true;
        }
    }
    ctl.perf_stats(&t, "after main args loop 2");
    proc_args_tail(&mut ctl, &args)?;

    // Sort chains_to_align.

    unique_sort(&mut ctl.gen_opt.chains_to_align);
    unique_sort(&mut ctl.gen_opt.chains_to_align2);
    unique_sort(&mut ctl.gen_opt.chains_to_jun_align);
    unique_sort(&mut ctl.gen_opt.chains_to_jun_align2);

    // Check for invalid variables in linear conditions.

    for i in 0..ctl.clono_filt_opt.bounds.len() {
        ctl.clono_filt_opt.bounds[i].require_valid_variables(&ctl)?;
    }
    if ctl.gen_opt.gene_scan_test.is_some() {
        ctl.gen_opt
            .gene_scan_test
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl)?;
        ctl.gen_opt
            .gene_scan_control
            .as_ref()
            .unwrap()
            .require_valid_variables(&ctl)?;
    }
    Ok(())
}
