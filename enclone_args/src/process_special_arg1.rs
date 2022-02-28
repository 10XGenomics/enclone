// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Process a special argument, i.e. one that does not fit into a neat bucket.

use crate::proc_args::test_writeable;
use crate::proc_args2::{is_simple_arg, is_usize_arg};
use enclone_core::cell_color::*;
use enclone_core::defs::EncloneControl;
use enclone_vars::encode_arith;
use evalexpr::build_operator_tree;
use expr_tools::test_functions_in_node;
use io_utils::path_exists;
use itertools::Itertools;
use std::fs::{remove_file, File};
use string_utils::{stringme, TextUtils};
use tilde_expand::tilde_expand;
use vector_utils::{unique_sort, VecUtils};

pub fn process_special_arg1(
    arg: &str,
    ctl: &mut EncloneControl,
    _metas: &mut Vec<String>,
    _metaxs: &mut Vec<String>,
    _xcrs: &mut Vec<String>,
    _using_plot: &mut bool,
) -> Result<bool, String> {
    // Process the argument.

    if is_simple_arg(arg, "SEQ")? {
        ctl.join_print_opt.seq = true;

    // Not movable.
    } else if arg.starts_with("PG_DIST=") {
        let dist = arg.after("PG_DIST=");
        if dist != "MFL" {
            return Err("\nCurrently the only allowed value for PG_DIST is MFL.\n".to_string());
        }
        ctl.gen_opt.peer_group_dist = dist.to_string();
    } else if is_simple_arg(arg, "H5")? {
        ctl.gen_opt.force_h5 = true;
    } else if is_simple_arg(arg, "NH5")? {
        ctl.gen_opt.force_h5 = false;
    } else if arg == "LEGEND" {
        ctl.plot_opt.use_legend = true;
    } else if arg == "MAX_HEAVIES=1" {
        ctl.gen_opt.max_heavies = 1;
    } else if arg.starts_with("ALIGN_2ND") {
        let n = arg.after("ALIGN_2ND");
        if n.parse::<usize>().is_err() || n.force_usize() == 0 {
            return Err(format!("\nArgument {} is not properly specified.\n", arg));
        }
        ctl.gen_opt.chains_to_align2.push(n.force_usize());
    } else if arg.starts_with("ALIGN") {
        let n = arg.after("ALIGN");
        if n.parse::<usize>().is_err() || n.force_usize() == 0 {
            return Err(format!("\nArgument {} is not properly specified.\n", arg));
        }
        ctl.gen_opt.chains_to_align.push(n.force_usize());
    } else if arg.starts_with("JALIGN_2ND") {
        let n = arg.after("JALIGN_2ND");
        if n.parse::<usize>().is_err() || n.force_usize() == 0 {
            return Err(format!("\nArgument {} is not properly specified.\n", arg));
        }
        ctl.gen_opt.chains_to_jun_align2.push(n.force_usize());
    } else if arg.starts_with("ALL_BC=") || arg.starts_with("ALL_BCH=") {
        let parts;
        if arg.starts_with("ALL_BC=") {
            parts = arg.after("ALL_BC=").split(',').collect::<Vec<&str>>();
        } else {
            parts = arg.after("ALL_BCH=").split(',').collect::<Vec<&str>>();
            ctl.gen_opt.all_bc_human = true;
        }
        if parts.is_empty() || parts[0].len() == 0 {
            return Err(
                "\nFor ALL_BC/ALL_BCH, at a minimum, a filename must be provided.\n".to_string(),
            );
        }
        if ctl.gen_opt.all_bc_filename.len() > 0 {
            return Err("\nThe argument ALL_BC/ALL_BCH may only be used once.\n".to_string());
        }
        ctl.gen_opt.all_bc_filename = parts[0].to_string();
        test_writeable(&ctl.gen_opt.all_bc_filename, ctl.gen_opt.evil_eye)?;
        for i in 1..parts.len() {
            ctl.gen_opt.all_bc_fields.push(parts[i].to_string());
        }
        ctl.gen_opt.all_bc_fields_orig = ctl.gen_opt.all_bc_fields.clone();
    } else if arg.starts_with("JOIN_BASIC=") {
        let val = arg.after("JOIN_BASIC=");
        if !val.parse::<f64>().is_ok() || val.force_f64() < 0.0 || val.force_f64() > 100.0 {
            return Err(
                "\nArgument to JOIN_BASIC needs to be a number between 0 and 100.\n".to_string(),
            );
        }
        ctl.join_alg_opt.basic = Some(val.force_f64());
    } else if arg.starts_with("JOIN_BASIC_H=") {
        let val = arg.after("JOIN_BASIC_H=");
        if !val.parse::<f64>().is_ok() || val.force_f64() < 0.0 || val.force_f64() > 100.0 {
            return Err(
                "\nArgument to JOIN_BASIC_H needs to be a number between 0 and 100.\n".to_string(),
            );
        }
        ctl.join_alg_opt.basic_h = Some(val.force_f64());
    } else if arg.starts_with("JOIN_CDR3_IDENT=") {
        let val = arg.after("JOIN_CDR3_IDENT=");
        if !val.parse::<f64>().is_ok() || val.force_f64() < 0.0 || val.force_f64() > 100.0 {
            return Err(
                "\nArgument to JOIN_CDR3_IDENT needs to be a number between 0 and 100.\n"
                    .to_string(),
            );
        }
        ctl.join_alg_opt.join_cdr3_ident = val.force_f64();
    } else if arg.starts_with("FWR1_CDR12_DELTA=") {
        let val = arg.after("FWR1_CDR12_DELTA=");
        if !val.parse::<f64>().is_ok() || val.force_f64() < 0.0 || val.force_f64() > 100.0 {
            return Err(
                "\nArgument to FWR1_CDR12_DELTA needs to be a number between 0 and 100.\n"
                    .to_string(),
            );
        }
        ctl.join_alg_opt.fwr1_cdr12_delta = val.force_f64();
    } else if arg.starts_with("HONEY=") {
        let mut parts = Vec::<Vec<String>>::new();
        {
            let subparts = arg.after("HONEY=").split(',').collect::<Vec<&str>>();
            if subparts.is_empty() || !subparts[0].contains('=') {
                return Err("\nSyntax for HONEY=... is incorrect.\n".to_string());
            }
            let mut part = Vec::<String>::new();
            for i in 0..subparts.len() {
                if subparts[i].contains('=') && !part.is_empty() {
                    parts.push(part.clone());
                    part.clear();
                }
                part.push(subparts[i].to_string());
            }
            if !part.is_empty() {
                parts.push(part);
            }
        }
        ctl.plot_opt.use_legend = true;
        let mut out_count = 0;
        let mut legend_count = 0;
        let mut color_count = 0;
        let (mut min, mut max) = (None, None);
        let (mut var, mut display_var) = (String::new(), String::new());
        let mut schema = String::new();
        for p in parts.iter() {
            let mut p = p.clone();
            let part_name = p[0].before("=").to_string();
            p[0] = p[0].after("=").to_string();
            let err = format!(
                "\nUnrecognized {} specification {}.\n",
                part_name,
                p.iter().format(",")
            );
            if part_name == "out" {
                if p.len() > 2 {
                    return Err(err);
                }
                let filename = p[0].clone();
                if p.len() == 2 && p[1].parse::<usize>().is_ok() {
                    ctl.plot_opt.png_width = Some(p[1].force_usize());
                    if !filename.ends_with(".png") {
                        return Err("\nWidth specification for the HONEY argument only \
                            makes sense if the filename ends with .png.\n"
                            .to_string());
                    }
                }
                if filename != "stdout"
                    && filename != "stdout.png"
                    && filename != "gui"
                    && !filename.ends_with(".svg")
                    && !filename.ends_with(".png")
                {
                    return Err(
                        "\nHONEY out filename needs to end with .svg or .png.\n".to_string()
                    );
                }
                ctl.plot_opt.plot_file = filename;
                out_count += 1;
            } else if part_name == "legend" {
                if p.solo() && p[0] == "none" {
                    ctl.plot_opt.use_legend = false;
                    legend_count += 1;
                } else {
                    return Err(err);
                }
            } else if part_name == "color" {
                color_count += 1;
                if p.len() == 1 && p[0] == "dataset" {
                    schema = "dataset".to_string();
                    let v = ColorByDataset {};
                    let cc = CellColor::ByDataset(v);
                    ctl.plot_opt.cell_color = cc;
                } else {
                    if p[0] != "var" || p.len() < 2 {
                        return Err(err);
                    }
                    schema = "variable".to_string();
                    var = p[1].to_string();
                    display_var = var.clone();
                    if var.contains(':') {
                        display_var = var.before(":").to_string();
                        var = var.after(":").to_string();
                    }
                    if p.len() >= 3 && !p[2].is_empty() && p[2] != "turbo" {
                        return Err(err);
                    }
                    if p.len() >= 4 {
                        let scale = &p[3..];
                        if !scale.is_empty() && scale[0] != "minmax" {
                            return Err(err);
                        }
                        if scale.len() >= 2 {
                            if scale[1].parse::<f64>().is_err() {
                                return Err(err);
                            }
                            min = Some(scale[1].force_f64());
                        }
                        if scale.len() >= 3 {
                            if scale[2].parse::<f64>().is_err() {
                                return Err(err);
                            }
                            max = Some(scale[2].force_f64());
                        }
                        if min.is_some() && max.is_some() && min >= max {
                            return Err(err);
                        }
                    }
                }
            } else {
                return Err(format!("\nUnrecognized specification {}=....\n", part_name));
            }
        }
        if out_count == 0 {
            return Err("\nHONEY=... must specify out=....\n".to_string());
        }
        if out_count > 1 {
            return Err("\nHONEY=... must specify out=... only once.\n".to_string());
        }
        if legend_count > 1 {
            return Err("\nHONEY=... may specify legend=... only once.\n".to_string());
        }
        if color_count == 0 {
            return Err("\nHONEY=... must specify color=....\n".to_string());
        }
        if color_count > 1 {
            return Err("\nHONEY=... must specify color=... only once.\n".to_string());
        }
        if schema == "dataset" {
            let v = ColorByDataset {};
            let cc = CellColor::ByDataset(v);
            ctl.plot_opt.cell_color = cc;
        } else if schema == "variable" {
            let v = ColorByVariableValue {
                var,
                display_var,
                min,
                max,
            };
            let cc = CellColor::ByVariableValue(v);
            ctl.plot_opt.cell_color = cc;
        }
    } else if arg.starts_with("VAR_DEF=") {
        let val = arg.after("VAR_DEF=");
        if !val.contains(':') {
            return Err(format!("\nCould not find : in {}.\n", arg));
        }
        let name = val.before(":");
        let expr = val.after(":");
        let eval = encode_arith(expr);
        let compiled = build_operator_tree(&eval);
        if compiled.is_err() {
            return Err(format!(
                "\nUnable to represent \"{}\" as a valid expression.  You might \
                check the following:\n\
                • arithmetic operators + - * / must have a blank on both sides\n\
                • parentheses must be balanced\n",
                expr,
            ));
        }
        let compiled = compiled.unwrap();
        let res = test_functions_in_node(&compiled);
        if res.is_err() {
            let err = res.as_ref().err().unwrap();
            return Err(format!(
                "\n{}\nYou might check the following:\n\
                • arithmetic operators + - * / must have a blank on both sides\n",
                err,
            ));
        }
        ctl.gen_opt
            .var_def
            .push((name.to_string(), eval, compiled, expr.to_string()));
    } else if arg.starts_with("MIN_DONORS") {
        let n = arg.after("MIN_DONORS=");
        if n.parse::<usize>().is_err() || n.force_usize() == 0 {
            return Err(format!("\nArgument {} is not properly specified.\n", arg));
        }
        let n = n.force_usize();
        ctl.clono_filt_opt.min_donors = n;
        if n >= 2 {
            ctl.clono_filt_opt_def.donor = true;
        }
    } else if arg.starts_with("JALIGN") {
        let n = arg.after("JALIGN");
        if n.parse::<usize>().is_err() || n.force_usize() == 0 {
            return Err(format!("\nArgument {} is not properly specified.\n", arg));
        }
        ctl.gen_opt.chains_to_jun_align.push(n.force_usize());
    } else if arg.starts_with("FB_SHOW=") {
        let fields = arg.after("FB_SHOW=").split(',').collect::<Vec<&str>>();
        let mut found_k = false;
        let mut ok = true;
        for i in 0..fields.len() {
            if fields[i].parse::<usize>().is_ok() {
                if found_k {
                    return Err("\nFB_SHOW argument contains more than one integer.\n".to_string());
                }
                found_k = true;
            } else {
                if fields[i].len() != 15 {
                    ok = false;
                }
                for c in fields[i].chars() {
                    if c != 'A' && c != 'C' && c != 'G' && c != 'T' {
                        ok = false;
                    }
                }
            }
        }
        if !ok {
            return Err("\nFB_SHOW argument must be a comma-separated list \
                containing at most one nonnegative integer and zero or more DNA \
                sequences of length 15 (in the alphabet A,C,G,T).\n"
                .to_string());
        }
        ctl.gen_opt.fb_show = arg.after("FB_SHOW=").to_string();
    } else if arg.starts_with("SIM_MAT_PLOT=") {
        let fields = arg.after("SIM_MAT_PLOT=").split(',').collect::<Vec<&str>>();
        if fields.len() < 2 {
            return Err(
                "\nSIM_MAT_PLOT requires at least two comma-separated arguments.\n".to_string(),
            );
        }
        let mut val = fields[0].to_string();
        val = stringme(&tilde_expand(val.as_bytes()));
        ctl.plot_opt.sim_mat_plot_file = val.clone();
        if val != "stdout" && val != "stdouth" && val != "gui" {
            let f = File::create(&val);
            if f.is_err() {
                let mut emsg = format!(
                    "\nYou've specified an output file\n{}\nthat cannot be written.\n",
                    val
                );
                if val.contains('/') {
                    let dir = val.rev_before("/");
                    let msg;
                    if path_exists(dir) {
                        msg = "exists";
                    } else {
                        msg = "does not exist";
                    }
                    emsg += &mut format!("Note that the path {} {}.\n", dir, msg);
                }
                return Err(emsg);
            }
            remove_file(&val).unwrap_or_else(|_| panic!("could not remove file {}", val));
        }
        ctl.plot_opt.sim_mat_plot_vars.clear();
        for j in 1..fields.len() {
            ctl.plot_opt.sim_mat_plot_vars.push(fields[j].to_string());
        }
    } else if arg.starts_with("G=") {
        let mut x = Vec::<usize>::new();
        if arg != "G=all" {
            let s = arg.after("G=").split(',').collect::<Vec<&str>>();
            let mut ok = false;
            for i in 0..s.len() {
                if s[i].parse::<usize>().is_ok() {
                    let n = s[i].force_usize();
                    if n >= 1 {
                        x.push(n);
                        ok = true;
                    }
                } else if s[i].contains('-') {
                    let (a, b) = (s[i].before("-"), s[i].after("-"));
                    if a.parse::<usize>().is_ok() && b.parse::<usize>().is_ok() {
                        let (a, b) = (a.force_usize(), b.force_usize());
                        if 1 <= a && a <= b {
                            for j in a..=b {
                                x.push(j);
                            }
                            ok = true;
                        }
                    }
                }
                if !ok {
                    return Err(
                        "\nArgument to G= must be a comma separated list of positive integers or \
                            hyphenated rangers of positive integers or all.\n"
                            .to_string(),
                    );
                }
            }
            unique_sort(&mut x);
        }
        ctl.gen_opt.group_post_filter = Some(x);
    } else if arg.starts_with("PLOTXY_EXACT=") {
        let fields = arg.after("PLOTXY_EXACT=").split(',').collect::<Vec<&str>>();
        if fields.len() != 3 && fields.len() != 4 {
            return Err(
                "\nPLOTXY_EXACT requires three or four comma-separated arguments.\n".to_string(),
            );
        }
        if fields.len() == 4 && fields[3] != "sym" {
            return Err(
                "\nIf four arguments are supplied to PLOTXY_EXACT, then the fourth argument \
                    must be sym.\n"
                    .to_string(),
            );
        }
        ctl.plot_opt.plot_xy_sym = fields.len() == 4;
        if fields[0].is_empty() || fields[1].is_empty() || fields[2].is_empty() {
            return Err("\nArguments to PLOTXY_EXACT must be non-null.\n".to_string());
        }
        let mut xvar = fields[0].to_string();
        let mut yvar = fields[1].to_string();
        if xvar.starts_with("log10(") && xvar.ends_with(')') {
            xvar = xvar.between("log10(", ")").to_string();
            ctl.plot_opt.plot_xy_x_log10 = true;
        }
        if yvar.starts_with("log10(") && yvar.ends_with(')') {
            yvar = yvar.between("log10(", ")").to_string();
            ctl.plot_opt.plot_xy_y_log10 = true;
        }
        ctl.plot_opt.plot_xy_xvar = xvar;
        ctl.plot_opt.plot_xy_yvar = yvar;
        let mut val = fields[2].to_string();
        val = stringme(&tilde_expand(val.as_bytes()));
        ctl.plot_opt.plot_xy_filename = val.clone();
        if val != "stdout" && val != "stdouth" && val != "gui" {
            let f = File::create(&val);
            if f.is_err() {
                let mut emsg = format!(
                    "\nYou've specified an output file\n{}\nthat cannot be written.\n",
                    val
                );
                if val.contains('/') {
                    let dir = val.rev_before("/");
                    let msg;
                    if path_exists(dir) {
                        msg = "exists";
                    } else {
                        msg = "does not exist";
                    }
                    emsg += &mut format!("Note that the path {} {}.\n", dir, msg);
                }
                return Err(emsg);
            }
            remove_file(&val).unwrap_or_else(|_| panic!("could not remove file {}", val));
        }
    } else if is_usize_arg(arg, "REQUIRED_FPS")? {
        ctl.gen_opt.required_fps = Some(arg.after("REQUIRED_FPS=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_CELLS")? {
        ctl.gen_opt.required_cells = Some(arg.after("REQUIRED_CELLS=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_DONORS")? {
        ctl.gen_opt.required_donors = Some(arg.after("REQUIRED_DONORS=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_CLONOTYPES")? {
        ctl.gen_opt.required_clonotypes = Some(arg.after("REQUIRED_CLONOTYPES=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_TWO_CELL_CLONOTYPES")? {
        ctl.gen_opt.required_two_cell_clonotypes =
            Some(arg.after("REQUIRED_TWO_CELL_CLONOTYPES=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_TWO_CHAIN_CLONOTYPES")? {
        ctl.gen_opt.required_two_chain_clonotypes =
            Some(arg.after("REQUIRED_TWO_CHAIN_CLONOTYPES=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_THREE_CHAIN_CLONOTYPES")? {
        ctl.gen_opt.required_three_chain_clonotypes =
            Some(arg.after("REQUIRED_THREE_CHAIN_CLONOTYPES=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_FOUR_CHAIN_CLONOTYPES")? {
        ctl.gen_opt.required_four_chain_clonotypes =
            Some(arg.after("REQUIRED_FOUR_CHAIN_CLONOTYPES=").force_usize());
    } else if is_usize_arg(arg, "REQUIRED_DATASETS")? {
        ctl.gen_opt.required_datasets = Some(arg.after("REQUIRED_DATASETS=").force_usize());
    } else if is_usize_arg(arg, "EXACT")? {
        ctl.gen_opt.exact = Some(arg.after("EXACT=").force_usize());
    } else if is_usize_arg(arg, "MIN_CHAINS")? {
        ctl.clono_filt_opt.min_chains = arg.after("MIN_CHAINS=").force_usize();
    } else if is_usize_arg(arg, "MAX_CHAINS")? {
        ctl.clono_filt_opt.max_chains = arg.after("MAX_CHAINS=").force_usize();
    } else if is_usize_arg(arg, "MIN_CELLS")? {
        ctl.clono_filt_opt.ncells_low = arg.after("MIN_CELLS=").force_usize();
    } else if is_usize_arg(arg, "MAX_CELLS")? {
        ctl.clono_filt_opt.ncells_high = arg.after("MAX_CELLS=").force_usize();
    } else if arg.starts_with("EXFASTA=") {
        ctl.gen_opt.fasta = arg.after("EXFASTA=").to_string();
    } else if arg.starts_with("FASTA=") {
        ctl.gen_opt.fasta_filename = arg.after("FASTA=").to_string();
    } else if arg.starts_with("FASTA_AA=") {
        ctl.gen_opt.fasta_aa_filename = arg.after("FASTA_AA=").to_string();

    // Other.
    } else if arg == "AGROUP" {
        if ctl.clono_group_opt.style == "symmetric" {
            return Err(
                "\nSymmetric and asymmetric grouping options cannot both be specified.\n"
                    .to_string(),
            );
        }
        ctl.clono_group_opt.style = "asymmetric".to_string();
    } else if arg == "GROUP_VJ_REFNAME" {
        ctl.clono_group_opt.style = "symmetric".to_string();
        ctl.clono_group_opt.vj_refname = true;
    } else if arg == "GROUP_VJ_REFNAME_HEAVY" {
        ctl.clono_group_opt.style = "symmetric".to_string();
        ctl.clono_group_opt.vj_heavy_refname = true;
    } else if arg == "GROUP_VDJ_REFNAME_HEAVY" {
        ctl.clono_group_opt.style = "symmetric".to_string();
        ctl.clono_group_opt.vdj_heavy_refname = true;
    } else if arg == "GROUP_VJ_REFNAME_STRONG" {
        ctl.clono_group_opt.style = "symmetric".to_string();
        ctl.clono_group_opt.vj_refname = true;
        ctl.clono_group_opt.vj_len = true;
        ctl.clono_group_opt.cdr3_len = true;
    } else {
        return Ok(false);
    }
    Ok(true)
}
