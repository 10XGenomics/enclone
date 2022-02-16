// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Process a special argument, i.e. one that does not fit into a neat bucket.

use crate::proc_args2::{is_f64_arg, is_usize_arg};
use enclone_core::defs::EncloneControl;
use enclone_core::linear_condition::LinearCondition;
use evalexpr::build_operator_tree;
use regex::Regex;
use string_utils::{parse_csv, stringme, TextUtils};
use tilde_expand::tilde_expand;
use vector_utils::unique_sort;

pub fn process_special_arg2(
    arg: &str,
    ctl: &mut EncloneControl,
    metas: &mut Vec<String>,
    metaxs: &mut Vec<String>,
    xcrs: &mut Vec<String>,
    using_plot: &mut bool,
) -> Result<(), String> {
    // Process the argument.

    if arg.starts_with("GROUP=") {
        if ctl.clono_group_opt.style == "asymmetric" {
            return Err(
                "\nSymmetric and asymmetric grouping options cannot both be specified.\n"
                    .to_string(),
            );
        }
        ctl.clono_group_opt.style = "symmetric".to_string();
        let c = arg.after("GROUP=").split(',').collect::<Vec<&str>>();
        for x in c.iter() {
            let x = *x;
            if x == "vj_refname" {
                ctl.clono_group_opt.vj_refname = true;
            } else if x == "vj_heavy_refname" {
                ctl.clono_group_opt.vj_heavy_refname = true;
            } else if x == "vdj_refname" {
                ctl.clono_group_opt.vdj_refname = true;
            } else if x == "vdj_heavy_refname" {
                ctl.clono_group_opt.vdj_heavy_refname = true;
            } else if x == "len" {
                ctl.clono_group_opt.vj_len = true;
            } else if x == "cdr3_len" {
                ctl.clono_group_opt.cdr3_len = true;
            } else if x.starts_with("≥light") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.light_pc = Some(val.force_f64());
            } else if x.starts_with("≥aa_light") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("light>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.light_pc = Some(val.force_f64());
            } else if x.starts_with("aa_light>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("light⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.light_pc = Some(val.force_f64());
            } else if x.starts_with("aa_light⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("heavy≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.heavy_pc = Some(val.force_f64());
            } else if x.starts_with("aa_heavy≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("heavy>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.heavy_pc = Some(val.force_f64());
            } else if x.starts_with("aa_heavy>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("heavy⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.heavy_pc = Some(val.force_f64());
            } else if x.starts_with("aa_heavy⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.aa_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_light≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_light≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_light>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_light>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_light⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_light⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_light in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_light_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_heavy≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_heavy≥") && x.ends_with('%') {
                let val = x.after("≥").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_heavy>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_heavy>=") && x.ends_with('%') {
                let val = x.after(">=").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_heavy⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_heavy_pc = Some(val.force_f64());
            } else if x.starts_with("cdr3_aa_heavy⩾") && x.ends_with('%') {
                let val = x.after("⩾").rev_before("%");
                if val.parse::<f64>().is_err() {
                    return Err("\nIllegal value for cdr3_aa_heavy in GROUP.\n".to_string());
                }
                ctl.clono_group_opt.cdr3_aa_heavy_pc = Some(val.force_f64());
            } else {
                return Err(format!(
                    "\nUnrecognized condition {} in GROUP argument.\n",
                    x
                ));
            }
        }
    } else if arg.starts_with("DIFF_STYLE=") {
        ctl.gen_opt.diff_style = arg.after("=").to_string();
        if ctl.gen_opt.diff_style != "C1" && ctl.gen_opt.diff_style != "C2" {
            return Err("\nThe only allowed values for DIFF_STYLE are C1 and C2.\n".to_string());
        }
    } else if arg.starts_with("COLOR=") {
        ctl.gen_opt.color = arg.after("COLOR=").to_string();
        if ctl.gen_opt.color != *"codon"
            && ctl.gen_opt.color != *"codon-diffs"
            && ctl.gen_opt.color != *"property"
        {
            let mut ok = false;
            if arg.starts_with("COLOR=peer.") {
                let pc = arg.after("COLOR=peer.");
                if pc.parse::<f64>().is_ok() {
                    let pc = pc.force_f64();
                    if (0.0..=100.0).contains(&pc) {
                        ok = true;
                        ctl.gen_opt.color_by_rarity_pc = pc;
                    }
                }
            }
            if !ok {
                return Err(
                    "\nThe specified value for COLOR is not allowed.  Please see \
                    \"enclone help color\".\n"
                        .to_string(),
                );
            }
        }
    } else if arg == "TREE" {
        ctl.gen_opt.tree_on = true;
    } else if arg == "TREE=const" {
        // this is for backward compatibility
        ctl.gen_opt.tree_on = true;
        ctl.gen_opt.tree.push("const1".to_string());
    } else if arg.starts_with("TREE=") {
        ctl.gen_opt.tree_on = true;
        let p = arg.after("TREE=").split(',').collect::<Vec<&str>>();
        for i in 0..p.len() {
            ctl.gen_opt.tree.push(p[i].to_string());
        }
    } else if arg.starts_with("FCELL=") // FCELL retained for backward compatibility
        || arg.starts_with("KEEP_CELL_IF")
    {
        let mut condition;
        if arg.starts_with("FCELL") {
            condition = arg.after("FCELL=").to_string();
        } else {
            condition = arg.after("KEEP_CELL_IF=").to_string();
        }
        let con = condition.as_bytes();
        for i in 0..con.len() {
            if i > 0
                && i < con.len() - 1
                && con[i] == b'='
                && con[i - 1] != b'='
                && con[i - 1] != b'<'
                && con[i - 1] != b'>'
                && con[i + 1] != b'='
            {
                return Err(format!(
                    "\nConstraints for {} cannot use =.  Please use == instead.\n",
                    arg.before("="),
                ));
            }
        }
        condition = condition.replace("'", "\"");
        let compiled = build_operator_tree(&condition);
        if compiled.is_err() {
            return Err(format!("\n{} usage incorrect.\n", arg.before("=")));
        }
        ctl.clono_filt_opt_def.fcell.push(compiled.unwrap());
    } else if arg.starts_with("LEGEND=") {
        let x = parse_csv(arg.after("LEGEND="));
        if x.is_empty() || x.len() % 2 != 0 {
            return Err("\nValue of LEGEND doesn't make sense.\n".to_string());
        }
        ctl.plot_opt.use_legend = true;
        for i in 0..x.len() / 2 {
            ctl.plot_opt
                .legend
                .push((x[2 * i].clone(), x[2 * i + 1].clone()));
        }
    } else if arg.starts_with("BARCODE=") {
        let bcs = arg.after("BARCODE=").split(',').collect::<Vec<&str>>();
        let mut x = Vec::<String>::new();
        for j in 0..bcs.len() {
            if !bcs[j].contains('-') {
                return Err(
                    "\nValue for a barcode in BARCODE argument is invalid, must contain -.\n"
                        .to_string(),
                );
            }
            x.push(bcs[j].to_string());
        }
        ctl.clono_filt_opt.barcode = x;
    } else if arg.starts_with("F=") {
        // deprecated but retained for backward compatibility
        let filt = arg.after("F=").to_string();
        ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt)?);
        ctl.clono_filt_opt.bound_type.push("mean".to_string());
    } else if arg.starts_with("KEEP_CLONO_IF_CELL_MEAN=") {
        let filt = arg.after("KEEP_CLONO_IF_CELL_MEAN=").to_string();
        ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt)?);
        ctl.clono_filt_opt.bound_type.push("mean".to_string());
    } else if arg.starts_with("KEEP_CLONO_IF_CELL_MIN=") {
        let filt = arg.after("KEEP_CLONO_IF_CELL_MIN=").to_string();
        ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt)?);
        ctl.clono_filt_opt.bound_type.push("min".to_string());
    } else if arg.starts_with("KEEP_CLONO_IF_CELL_MAX=") {
        let filt = arg.after("KEEP_CLONO_IF_CELL_MAX=").to_string();
        ctl.clono_filt_opt.bounds.push(LinearCondition::new(&filt)?);
        ctl.clono_filt_opt.bound_type.push("max".to_string());
    } else if arg.starts_with("SCAN=") {
        let mut x = arg.after("SCAN=").to_string();
        x = x.replace(" ", "");
        let x = x.split(',').collect::<Vec<&str>>();
        if x.len() != 3 {
            return Err("\nArgument to SCAN must have three components.\n".to_string());
        }
        ctl.gen_opt.gene_scan_test = Some(LinearCondition::new(x[0])?);
        ctl.gen_opt.gene_scan_control = Some(LinearCondition::new(x[1])?);
        let threshold = LinearCondition::new(x[2])?;
        for i in 0..threshold.var.len() {
            if threshold.var[i] != *"t" && threshold.var[i] != *"c" {
                return Err("\nIllegal variable in threshold for scan.\n".to_string());
            }
        }
        ctl.gen_opt.gene_scan_threshold = Some(threshold);
    } else if arg.starts_with("PLOT=") {
        *using_plot = true;
        let x = arg.after("PLOT=").split(',').collect::<Vec<&str>>();
        if x.is_empty() {
            return Err("\nArgument to PLOT is invalid.\n".to_string());
        }
        ctl.plot_opt.plot_file = x[0].to_string();
        for j in 1..x.len() {
            if !x[j].contains("->") {
                return Err("\nArgument to PLOT is invalid.\n".to_string());
            }
            ctl.gen_opt
                .origin_color_map
                .insert(x[j].before("->").to_string(), x[j].after("->").to_string());
        }
    } else if arg.starts_with("PLOT2=") {
        *using_plot = true;
        let x = arg.after("PLOT2=").split(',').collect::<Vec<&str>>();
        if x.is_empty() {
            return Err("\nArgument to PLOT is invalid.\n".to_string());
        }
        if x.len() % 2 != 1 {
            return Err("\nArgument to PLOT is invalid.\n".to_string());
        }
        ctl.plot_opt.plot_file = x[0].to_string();
        for j in (1..x.len()).step_by(2) {
            let condition = x[j].to_string();
            let color = x[j + 1].to_string();
            if !condition.contains('=') {
                return Err("\nArgument to PLOT is invalid.\n".to_string());
            }
            ctl.plot_opt.plot_conditions.push(condition);
            ctl.plot_opt.plot_colors.push(color);
        }
    } else if arg.starts_with("PLOT_BY_ISOTYPE=") {
        ctl.plot_opt.plot_by_isotype = true;
        ctl.plot_opt.plot_file = arg.after("PLOT_BY_ISOTYPE=").to_string();
        if ctl.plot_opt.plot_file.is_empty() {
            return Err("\nFilename value needs to be supplied to PLOT_BY_ISOTYPE.\n".to_string());
        }
    } else if arg.starts_with("PLOT_BY_ISOTYPE_COLOR=") {
        if arg.after("PLOT_BY_ISOTYPE_COLOR=").is_empty() {
            return Err(
                "\nA value needs to be specified for the PLOT_BY_ISOTYPE_COLOR \
                argument.\n"
                    .to_string(),
            );
        }
        let fields = arg
            .after("PLOT_BY_ISOTYPE_COLOR=")
            .split(',')
            .collect::<Vec<&str>>();
        for i in 0..fields.len() {
            ctl.plot_opt
                .plot_by_isotype_color
                .push(fields[i].to_string());
        }
    } else if arg.starts_with("PLOT_BY_MARK=") {
        ctl.plot_opt.plot_by_mark = true;
        ctl.plot_opt.plot_file = arg.after("PLOT_BY_MARK=").to_string();
        if ctl.plot_opt.plot_file.is_empty() {
            return Err("\nFilename value needs to be supplied to PLOT_BY_MARK.\n".to_string());
        }
    } else if is_usize_arg(arg, "MAX_CORES")? {
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
                if ctl.parseable_opt.pcols_sortx[j].contains(':') {
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
                return Err("\nIllegal value for VJ, must be over alphabet ACGT.\n".to_string());
            }
        }
    } else if arg.starts_with("AMINO=") {
        ctl.clono_print_opt.amino.clear();
        for x in arg.after("AMINO=").split(',').collect::<Vec<&str>>() {
            if !x.is_empty() {
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
                if start.parse::<usize>().is_ok()
                    && stop.parse::<usize>().is_ok()
                    && start.force_usize() <= stop.force_usize()
                {
                    ok = true;
                }
            }
            if !ok {
                return Err(format!(
                    "\nUnrecognized variable {} for AMINO.  Please type \
                     \"enclone help amino\".\n",
                    x
                ));
            }
        }
    } else if arg.starts_with("CVARS=") {
        ctl.clono_print_opt.cvars.clear();
        for x in arg.after("CVARS=").split(',').collect::<Vec<&str>>() {
            if !x.is_empty() {
                ctl.clono_print_opt.cvars.push(x.to_string());
            }
        }
        for x in ctl.clono_print_opt.cvars.iter_mut() {
            *x = x.replace("_sum", "_Σ");
            *x = x.replace("_mean", "_μ");
        }
    } else if arg.starts_with("CVARSP=") {
        for x in arg.after("CVARSP=").split(',').collect::<Vec<&str>>() {
            if !x.is_empty() {
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
    } else if arg.starts_with("DVARS=") {
        ctl.gen_opt.dvars.clear();
        for x in arg.after("DVARS=").split(',').collect::<Vec<&str>>() {
            ctl.gen_opt.dvars.push(x.to_string());
        }
    } else if arg.starts_with("GVARS=") {
        ctl.gen_opt.gvars.clear();
        for x in arg.after("GVARS=").split(',').collect::<Vec<&str>>() {
            ctl.gen_opt.gvars.push(x.to_string());
        }
    } else if is_f64_arg(arg, "MAX_SCORE")? {
        ctl.join_alg_opt.max_score = arg.after("MAX_SCORE=").force_f64();
    } else if is_f64_arg(arg, "MAX_LOG_SCORE")? {
        let x = arg.after("MAX_LOG_SCORE=").force_f64();
        ctl.join_alg_opt.max_score = 10.0_f64.powf(x);
    } else if arg.starts_with("CONST_IGH=") {
        let reg = Regex::new(&format!("^{}$", arg.after("CONST_IGH=")));
        if reg.is_err() {
            return Err(format!(
                "\nYour CONST_IGH value {} could not be parsed as a regular expression.\n",
                arg.after("CONST_IGH=")
            ));
        }
        ctl.clono_filt_opt.const_igh = Some(reg.unwrap());
    } else if arg.starts_with("CONST_IGKL=") {
        let reg = Regex::new(&format!("^{}$", arg.after("CONST_IGKL=")));
        if reg.is_err() {
            return Err(format!(
                "\nYour CONST_IGKL value {} could not be parsed as a regular expression.\n",
                arg.after("CONST_IGKL=")
            ));
        }
        ctl.clono_filt_opt.const_igkl = Some(reg.unwrap());
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
            let mut ok = true;
            for i in 0..fields.len() {
                if !fields[i].contains('~') {
                    ok = false;
                } else {
                    let _f1 = fields[i].before("~");
                    let f2 = fields[i].after("~");
                    if !f2.parse::<usize>().is_ok() {
                        ok = false;
                    }
                }
            }
            if !ok {
                return Err(format!(
                    "\nLooks like your CDR3 value {} is trying to be an Levenshtein distance\n\
                        pattern, but it is not.\n",
                    arg.after("CDR3=")
                ));
            }
            ctl.clono_filt_opt.cdr3_lev = arg.after("=").to_string();
        } else {
            let reg = Regex::new(&format!("^{}$", arg.after("CDR3=")));
            if reg.is_err() {
                return Err(format!(
                    "\nYour CDR3 value {} could not be parsed as a regular expression.\n",
                    arg.after("CDR3=")
                ));
            }
            ctl.clono_filt_opt.cdr3 = Some(reg.unwrap());
        }
    } else if is_usize_arg(arg, "CHAINS")? {
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
            if x.parse::<i32>().is_err() {
                return Err("\nInvalid argument to SEGN.\n".to_string());
            }
            y.push(x.to_string());
        }
        y.sort();
        ctl.clono_filt_opt.segn.push(y);
    } else if arg.starts_with("NSEG=") {
        let fields = arg.after("NSEG=").split('|').collect::<Vec<&str>>();
        let mut y = Vec::<String>::new();
        for x in fields.iter() {
            y.push(x.to_string());
        }
        y.sort();
        ctl.clono_filt_opt.nseg.push(y);
    } else if arg.starts_with("NSEGN=") {
        let fields = arg.after("NSEGN=").split('|').collect::<Vec<&str>>();
        let mut y = Vec::<String>::new();
        for x in fields.iter() {
            if x.parse::<i32>().is_err() {
                return Err("\nInvalid argument to NSEGN.\n".to_string());
            }
            y.push(x.to_string());
        }
        y.sort();
        ctl.clono_filt_opt.nsegn.push(y);
    } else if is_usize_arg(arg, "CELLS")? {
        ctl.clono_filt_opt.ncells_low = arg.after("CELLS=").force_usize();
        ctl.clono_filt_opt.ncells_high = ctl.clono_filt_opt.ncells_low;
    } else if arg.starts_with("META=") {
        let f = arg.after("META=");
        let f = stringme(&tilde_expand(f.as_bytes()));
        metas.push(f);
    } else if arg.starts_with("METAX=") {
        let f = arg.after("METAX=");
        let f = f.chars().filter(|c| !c.is_whitespace()).collect();
        metaxs.push(f);
    } else if arg.starts_with("TCR=")
        || arg.starts_with("BCR=")
        || arg.starts_with("TCRGD=")
        || (!arg.is_empty() && arg.as_bytes()[0] >= b'0' && arg.as_bytes()[0] <= b'9')
    {
        xcrs.push(arg.to_string());
    } else if arg != "--help" {
        return Err(format!("\nUnrecognized argument {}.\n", arg));
    }
    Ok(())
}
