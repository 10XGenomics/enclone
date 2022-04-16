// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use self::refx::{make_vdj_ref_data_core, RefData};
use crate::blacklist::profiling_blacklist;
use crate::determine_ref::determine_ref;
use crate::sec_mem::test_sec_mem;
use crate::setup::{critical_args, setup};
use crate::stop::main_enclone_stop;
use enclone::innate::species;
use enclone::secret::fetch_secmem;
use enclone_args::load_gex::get_gex_info;
use enclone_args::proc_args2::is_simple_arg;
use enclone_args::proc_args_check::{
    check_gvars, check_lvars, check_one_lvar, check_pcols, get_known_features,
};
use enclone_core::cell_color::CellColor;
use enclone_core::defs::EncloneControl;
use enclone_core::enclone_structs::*;
use enclone_core::version_string;
use enclone_stuff::start::*;
use enclone_stuff::vars::match_vars;
use enclone_vars::decode_arith;
use expr_tools::vars_of_node;
use io_utils::{open_for_read, open_userfile_for_read, path_exists};
use itertools::Itertools;
#[cfg(not(target_os = "windows"))]
use pretty_trace::start_profiling;
use std::{collections::HashMap, env, fs, fs::read_to_string, io::BufRead, time::Instant};
use string_utils::TextUtils;
use vdj_ann::refx;
use vector_utils::{bin_member, next_diff, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone(args: &Vec<String>) -> Result<EncloneState, String> {
    let setup = main_enclone_setup(args)?;
    if setup.tall.is_none() {
        return Ok(EncloneState::default());
    }
    let inter = main_enclone_start(setup)?;
    if inter.setup.tall.is_none() {
        return Ok(EncloneState::default());
    }
    main_enclone_stop(inter)
}

pub fn main_enclone_setup(args: &Vec<String>) -> Result<EncloneSetup, String> {
    let tall = Instant::now();

    // Test for enclone --check.

    if args.len() == 2 && args[1] == "--check" {
        let version1 = env!("CARGO_PKG_VERSION");
        let home = dirs::home_dir().unwrap().to_str().unwrap().to_string();
        let version_file = format!("{}/enclone/version", home);
        if !path_exists(&version_file) {
            return Err("\nError: the file ~/enclone/version does not exist.\n\
                Please visit bit.ly/enclone_install_issues.\n"
                .to_string());
        }
        let mut version2 = read_to_string(&version_file).unwrap();
        if !version2.starts_with('v') || !version2.ends_with('\n') {
            return Err(format!(
                "\nThe file ~/enclone/version appears to be damaged.\n\
                Its content is \"{}\".\n\
                Please visit bit.ly/enclone_install_issues.\n",
                version2,
            ));
        }
        version2 = version2.between("v", "\n").to_string();
        if version2 != version1 {
            return Err(format!(
                "\nError: enclone sees version {} but you downloaded version {}.\n\
                Please visit bit.ly/enclone_install_issues.\n",
                version1, version2
            ));
        }
        println!("\nCheck complete: it appears that your install of enclone was successful!\n");
        print!("Your version is: ");
        println!("{} : {}.\n", env!("CARGO_PKG_VERSION"), version_string());
        return Ok(EncloneSetup::default());
    }

    // Set up stuff, read args, etc.

    let args_orig = args.clone();
    let mut ctl = EncloneControl::default();
    let args = critical_args(args, &mut ctl)?;
    ctl.start_time = Some(tall);
    for i in 0..args.len() {
        let arg = &args[i];
        if arg == "PROFILE" {
            ctl.gen_opt.profile = true;
        }
    }
    if ctl.gen_opt.profile {
        start_profiling(&profiling_blacklist());
    }
    let (mut comp, mut comp2) = (false, false);
    for i in 1..args.len() {
        if args[i] == "PRINT_CPU" {
            ctl.gen_opt.print_cpu = true;
        }
        if args[i] == "PRINT_CPU_INFO" {
            ctl.gen_opt.print_cpu_info = true;
        }
        if args[i] == "COMP" || args[i] == "COMPE" {
            comp = true;
        }
        if args[i] == "COMPE" {
            ctl.perf_opt.comp_enforce = true;
        }
        if args[i] == "COMP2" {
            comp2 = true;
        }
    }
    if comp && !comp2 {
        println!();
    }
    ctl.gen_opt.cpu_all_start = 0;
    ctl.gen_opt.cpu_this_start = 0;
    if ctl.gen_opt.print_cpu || ctl.gen_opt.print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        if let Some(line) = f.lines().next() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            ctl.gen_opt.cpu_all_start = t.before(" ").force_usize();
        }
        let f = open_for_read![&format!("/proc/{}/stat", std::process::id())];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            ctl.gen_opt.cpu_this_start = fields[13].force_usize();
        }
    }
    if args_orig.len() == 2 && (args_orig[1] == "version" || args_orig[1] == "--version") {
        println!("{} : {}", env!("CARGO_PKG_VERSION"), version_string());
        return Ok(EncloneSetup::default());
    }
    if ctl.gen_opt.evil_eye {
        println!("calling perf_stats, before setup");
    }
    ctl.perf_stats(&tall, "before setup");
    let mut argsx = Vec::<String>::new();
    setup(&mut ctl, &args, &mut argsx, &args_orig)?;
    if ctl.gen_opt.split {
        return Ok(EncloneSetup::default());
    }
    let mut argsy = Vec::<String>::new();
    for i in 0..args_orig.len() {
        if args_orig[i] != "HTML"
            && args_orig[i] != "STABLE_DOC"
            && args_orig[i] != "NOPAGER"
            && args_orig[i] != "FORCE_EXTERNAL"
            && args_orig[i] != "NO_KILL"
            && !args_orig[i].starts_with("PRE=")
            && !args_orig[i].starts_with("PREPOST=")
            && !args_orig[i].starts_with("MAX_CORES=")
        {
            argsy.push(args_orig[i].clone());
        }
    }
    if argsy.len() == 1 || (argsy.len() > 1 && (argsy[1] == "help" || argsy[1] == "--help")) {
        return Ok(EncloneSetup::default());
    }

    // Dump internal ids.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "DUMP_INTERNAL_IDS")? {
            let mut x = Vec::<usize>::new();
            for y in ctl.origin_info.dataset_id.iter() {
                x.push(y.force_usize());
            }
            x.sort_unstable();
            println!("\n{}\n", x.iter().format(","));
            return Ok(EncloneSetup::default());
        }
    }

    // Read external data.

    if !ctl.gen_opt.ext.is_empty() {
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

    let gex_info = get_gex_info(&mut ctl)?;
    check_lvars(&ctl, &gex_info)?;
    let twoof = Instant::now();
    check_gvars(&ctl)?;
    check_pcols(
        &ctl,
        &gex_info,
        &ctl.parseable_opt.pcols,
        ctl.parseable_opt.pbarcode,
    )?;
    check_pcols(
        &ctl,
        &gex_info,
        &ctl.gen_opt.tree,
        ctl.parseable_opt.pbarcode,
    )?;
    if !ctl.plot_opt.plot_xy_filename.is_empty() {
        check_pcols(
            &ctl,
            &gex_info,
            &vec![
                ctl.plot_opt.plot_xy_xvar.clone(),
                ctl.plot_opt.plot_xy_yvar.clone(),
            ],
            ctl.parseable_opt.pbarcode,
        )?;
    }
    match ctl.plot_opt.cell_color {
        CellColor::ByVariableValue(ref x) => {
            check_pcols(&ctl, &gex_info, &vec![x.var.clone()], true)?;
        }
        CellColor::ByCategoricalVariableValue(ref x) => {
            check_pcols(&ctl, &gex_info, &x.vars, true)?;
        }
        _ => {}
    };
    let mut bound_vars = Vec::<String>::new();
    for bi in 0..ctl.clono_filt_opt.bounds.len() {
        let x = &ctl.clono_filt_opt.bounds[bi];
        for i in 0..x.n() {
            bound_vars.push(x.var[i].clone());
        }
    }
    unique_sort(&mut bound_vars);
    check_pcols(&ctl, &gex_info, &bound_vars, ctl.parseable_opt.pbarcode)?;
    check_pcols(
        &ctl,
        &gex_info,
        &ctl.plot_opt.sim_mat_plot_vars,
        ctl.parseable_opt.pbarcode,
    )?;
    let mut var_def_vars = Vec::<String>::new();
    for i in 0..ctl.gen_opt.var_def.len() {
        let n = &ctl.gen_opt.var_def[i].2;
        let vars = vars_of_node(&n);
        for v in vars.iter() {
            let w = decode_arith(&*v);
            var_def_vars.push(w);
        }
    }
    check_pcols(&ctl, &gex_info, &var_def_vars, ctl.parseable_opt.pbarcode)?;

    // Check ALL_BC arguments.

    if ctl.gen_opt.all_bc_filename.len() > 0 {
        if gex_info.gex_barcodes.is_empty() {
            return Err(format!("\nYou can't use ALL_BC with VDJ data alone.\n"));
        }
        for li in 0..ctl.origin_info.n() {
            if !gex_info.gex_matrices[li].initialized() {
                return Err(format!(
                    "\nALL_BC only works if feature_barcode_matrix.bin has been \
                    generated.\nPlease type \"enclone help input\" for more information.\n"
                ));
            }
        }
        let known_features = get_known_features(&gex_info)?;
        let extras = ["gex", "type", "clust", "cell"];
        for i in 0..ctl.gen_opt.all_bc_fields.len() {
            let var = &ctl.gen_opt.all_bc_fields[i];
            let mut ok = false;
            if bin_member(&known_features, &var) {
                ok = true;
            }
            let mut nd_var = false;
            for j in 0..extras.len() {
                if *var == extras[j] {
                    if *var != "cell" {
                        check_one_lvar(
                            var,
                            &ctl,
                            &gex_info,
                            &mut nd_var,
                            &Vec::<String>::new(),
                            true,
                        )?;
                    }
                    ok = true;
                }
            }
            if !ok {
                return Err(format!("\nIllegal variable {} in ALL_BC.\n", var));
            }
        }
    }
    ctl.perf_stats(&twoof, "checking pcols");

    // Check DVARS.

    let tfcell = Instant::now();
    if !ctl.gen_opt.dvars.is_empty() {
        let known_features = get_known_features(&gex_info)?;
        for j in 0..ctl.gen_opt.dvars.len() {
            let mut var = ctl.gen_opt.dvars[j].clone();
            if var.contains(':') {
                var = var.after(":").to_string();
            }
            let mut found = false;
            for k in 0..gex_info.json_metrics.len() {
                if gex_info.json_metrics[k].contains_key(&var.to_string()) {
                    found = true;
                }
            }
            if !found {
                let feature;
                if var.ends_with("_cellular_r") {
                    feature = var.before("_cellular_r").to_string();
                } else if var.ends_with("_cellular_u") {
                    feature = var.before("_cellular_u").to_string();
                } else {
                    return Err(format!("\nUnknown DVAR = {}.\n", var));
                }
                if !bin_member(&known_features, &feature) {
                    return Err(format!(
                        "\nIn DVAR = {}, the feature {} is unknown.\n",
                        var, feature,
                    ));
                }
            }
        }
    }

    // Test fcell.

    if !ctl.clono_filt_opt_def.fcell.is_empty() {
        let mut alt_bcs = Vec::<String>::new();
        for li in 0..ctl.origin_info.alt_bc_fields.len() {
            for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
                alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
            }
        }
        unique_sort(&mut alt_bcs);
        let mut test2 = Vec::<String>::new();
        for con in ctl.clono_filt_opt_def.fcell.iter() {
            for var in con.iter_variable_identifiers() {
                if !bin_member(&alt_bcs, &var.to_string()) {
                    test2.push(var.to_string());
                }
            }
            if let Some(_) = con.iter_function_identifiers().next() {
                return Err(
                    "\nSomething is wrong with your KEEP_CELL_IF or FCELL value.\n".to_string(),
                );
            }
        }
        if !test2.is_empty() {
            let known_features = get_known_features(&gex_info)?; // note duplicated computation
            for var in test2.iter() {
                if !bin_member(&known_features, var) {
                    return Err(format!(
                        "\nYou've used an illegal variable {} as part of an KEEP_CELL_IF or FCELL constraint.\n",
                        var
                    ));
                }
            }
        }
    }
    ctl.perf_stats(&tfcell, "checking fcell");

    // Find matching features for <regular expression>_g etc.

    match_vars(&mut ctl, &gex_info)?;

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
    determine_ref(&mut ctl, &mut refx)?;
    if refx.is_empty() && ctl.origin_info.n() == 0 {
        return Err("\nNo data and no TCR or BCR data have been specified.\n".to_string());
    }
    ctl.perf_stats(&tr, "starting reference");

    // Build reference data.

    let tr = Instant::now();
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
    make_vdj_ref_data_core(&mut refdata, refx2, &ext_refx, is_tcr, is_bcr, None);
    let mut to_ref_index = HashMap::<usize, usize>::new();
    for i in 0..refdata.refs.len() {
        to_ref_index.insert(refdata.id[i] as usize, i);
    }

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata);

    // Process for sec (secreted) or mem (membrane) if specified.

    test_sec_mem(&mut ctl)?;
    if ctl.gen_opt.using_secmem {
        fetch_secmem(&mut ctl)?;
    }
    ctl.perf_stats(&tr, "building reference and other things");

    // Get VDJ data paths.

    for li in 0..ctl.origin_info.dataset_path.len() {
        let json = format!("{}/{}", ctl.origin_info.dataset_path[li], ann);
        let json_lz4 = format!("{}/{}.lz4", ctl.origin_info.dataset_path[li], ann);
        if !path_exists(&json) && !path_exists(&json_lz4) {
            return Err(format!("\ncan't find {} or {}\n", json, json_lz4));
        } else if path_exists(&json) {
            ctl.pathlist.push(json);
        } else {
            ctl.pathlist.push(json_lz4);
        }
    }

    // Get last modified info for pathlist.

    for i in 0..ctl.pathlist.len() {
        let metadata = fs::metadata(&ctl.pathlist[i]);
        if metadata.is_err() {
            return Err(format!(
                "\nUnable to get file metadata for {}.\n",
                ctl.pathlist[i],
            ));
        }
        let modified = metadata.unwrap().modified();
        if modified.is_err() {
            return Err(format!(
                "\nUnable to determine modification date of {}.\n",
                ctl.pathlist[i],
            ));
        } else {
            ctl.last_modified.push(modified.unwrap());
        }
    }

    // Return.

    Ok(EncloneSetup {
        ctl,
        refdata,
        ann: ann.to_string(),
        gex_info,
        tall: Some(tall),
        is_bcr,
        to_ref_index,
    })
}
