// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use self::refx::*;
use crate::blacklist::*;
use crate::determine_ref::*;
use crate::disintegrate::*;
use crate::fcell::*;
use crate::filter_umi::*;
use crate::flag_defective::*;
use crate::inconsistent::*;
use crate::populate_features::*;
use crate::sec_mem::*;
use crate::setup::*;
use crate::some_filters::*;
use crate::stop::*;
use crate::vars::*;
use debruijn::dna_string::DnaString;
use enclone::allele::*;
use enclone::graph_filter::*;
use enclone::info::*;
use enclone::innate::*;
use enclone::join::*;
use enclone::misc1::*;
use enclone::misc2::*;
use enclone::misc3::*;
use enclone::secret::*;
use enclone_args::load_gex::*;
use enclone_args::proc_args2::*;
use enclone_args::proc_args_check::*;
use enclone_args::read_json::*;
use enclone_core::defs::*;
use enclone_core::*;
use enclone_print::loupe::*;
use equiv::EquivRel;
use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;
use std::{
    collections::HashMap,
    env, fs,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    time::Instant,
};
use stirling_numbers::*;
use string_utils::*;
use vdj_ann::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Clone, Debug, Default)]
pub struct MainEncloneOutput {
    pub pics: Vec<String>, // clonotype tables
    pub last_widths: Vec<usize>,
    pub svgs: Vec<String>, // SVG objects
    pub noprint: bool,
    pub noprintx: bool,
    pub html: bool,
    pub ngroup: bool,
    pub pretty: bool,
}

#[derive(Default)]
pub struct EncloneState {
    pub inter: EncloneIntermediates,
    pub outs: MainEncloneOutput,
}

#[derive(Default)]
pub struct EncloneSetup {
    pub ctl: EncloneControl,
    pub ann: String,
    pub gex_info: GexInfo,
    pub tall: Option<Instant>,
    pub refdata: RefData,
    pub is_bcr: bool,
    pub to_ref_index: HashMap<usize, usize>,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone(args: &Vec<String>) -> Result<EncloneState, String> {
    let setup = main_enclone_setup(&args)?;
    if setup.tall.is_none() {
        return Ok(EncloneState::default());
    }
    let inter = main_enclone_start(setup)?;
    if inter.setup.tall.is_none() {
        return Ok(EncloneState::default());
    }
    Ok(main_enclone_stop(inter)?)
}

pub fn main_enclone_setup(args: &Vec<String>) -> Result<EncloneSetup, String> {
    let tall = Instant::now();
    let args_orig = args.clone();
    let args = process_source(&args)?;

    // Set up stuff, read args, etc.

    let mut ctl = EncloneControl::default();
    ctl.start_time = Some(tall.clone());
    for i in 0..args.len() {
        let arg = &args[i];
        if arg == "PROFILE" {
            ctl.gen_opt.profile = true;
        }
        if arg == "EVIL_EYE" {
            ctl.gen_opt.evil_eye = true;
        }
    }
    if ctl.gen_opt.evil_eye {
        println!("the evil eye is on");
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
        println!("");
    }
    ctl.gen_opt.cpu_all_start = 0;
    ctl.gen_opt.cpu_this_start = 0;
    if ctl.gen_opt.print_cpu || ctl.gen_opt.print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            ctl.gen_opt.cpu_all_start = t.before(" ").force_usize();
            break;
        }
        let f = open_for_read![&format!("/proc/{}/stat", std::process::id())];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            ctl.gen_opt.cpu_this_start = fields[13].force_usize();
        }
    }
    if args.len() == 2 && (args[1] == "version" || args[1] == "--version") {
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
    if argsx.len() == 1 || (argsx.len() > 1 && (argsx[1] == "help" || argsx[1] == "--help")) {
        return Ok(EncloneSetup::default());
    }

    // Dump internal ids.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "DUMP_INTERNAL_IDS")? {
            let mut x = Vec::<usize>::new();
            for y in ctl.origin_info.dataset_id.iter() {
                x.push(y.force_usize());
            }
            x.sort();
            println!("\n{}\n", x.iter().format(","));
            return Ok(EncloneSetup::default());
        }
    }

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

    let gex_info = get_gex_info(&mut ctl)?;
    check_lvars(&ctl, &gex_info)?;
    let twoof = Instant::now();
    check_gvars(&ctl)?;
    check_pcols(&ctl, &gex_info, &ctl.parseable_opt.pcols)?;
    check_pcols(&ctl, &gex_info, &ctl.gen_opt.tree)?;
    if ctl.plot_opt.plot_xy_filename.len() > 0 {
        check_pcols(
            &ctl,
            &gex_info,
            &vec![
                ctl.plot_opt.plot_xy_xvar.clone(),
                ctl.plot_opt.plot_xy_yvar.clone(),
            ],
        )?;
    }
    let mut bound_vars = Vec::<String>::new();
    for bi in 0..ctl.clono_filt_opt.bounds.len() {
        let x = &ctl.clono_filt_opt.bounds[bi];
        for i in 0..x.n() {
            bound_vars.push(x.var[i].clone());
        }
    }
    unique_sort(&mut bound_vars);
    check_pcols(&ctl, &gex_info, &bound_vars)?;
    check_pcols(&ctl, &gex_info, &ctl.plot_opt.sim_mat_plot_vars)?;
    ctl.perf_stats(&twoof, "checking pcols");

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
    if refx.len() == 0 && ctl.origin_info.n() == 0 {
        return Err(format!(
            "\nNo data and no TCR or BCR data have been specified.\n"
        ));
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
    make_vdj_ref_data_core(&mut refdata, &refx2, &ext_refx, is_tcr, is_bcr, None);
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
        } else {
            if path_exists(&json) {
                ctl.pathlist.push(json);
            } else {
                ctl.pathlist.push(json_lz4);
            }
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
        ctl: ctl,
        refdata: refdata,
        ann: ann.to_string(),
        gex_info: gex_info,
        tall: Some(tall),
        is_bcr: is_bcr,
        to_ref_index: to_ref_index,
    })
}

pub fn main_enclone_start(setup: EncloneSetup) -> Result<EncloneIntermediates, String> {
    let tr = Instant::now();
    let ctl = &setup.ctl;
    let gex_info = &setup.gex_info;
    let refdata = &setup.refdata;
    let is_bcr = setup.is_bcr;
    let to_ref_index = &setup.to_ref_index;

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Flag defective reference sequences.

    let mut log = Vec::<u8>::new();
    let mut broken = Vec::<bool>::new();
    flag_defective(&ctl, &refdata, &mut log, &mut broken);
    ctl.perf_stats(&tr, "flagging defective references");

    // Parse the json annotations file.

    let tparse = Instant::now();
    let mut tig_bc = Vec::<Vec<TigData>>::new();
    let mut vdj_cells = Vec::<Vec<String>>::new();
    let mut gex_cells = Vec::<Vec<String>>::new();
    let mut gex_cells_specified = Vec::<bool>::new();
    parse_json_annotations_files(
        &ctl,
        &mut tig_bc,
        &refdata,
        &to_ref_index,
        &mut vdj_cells,
        &mut gex_cells,
        &mut gex_cells_specified,
    )?;
    ctl.perf_stats(&tparse, "loading from json");

    // Populate features.

    let tpop = Instant::now();
    let mut fr1_starts = Vec::<usize>::new();
    let mut fr2_starts = Vec::<Option<usize>>::new();
    let mut fr3_starts = Vec::<Option<usize>>::new();
    let mut cdr1_starts = Vec::<Option<usize>>::new();
    let mut cdr2_starts = Vec::<Option<usize>>::new();
    populate_features(
        &ctl,
        &refdata,
        &broken,
        &mut fr1_starts,
        &mut fr2_starts,
        &mut fr3_starts,
        &mut cdr1_starts,
        &mut cdr2_starts,
        &mut log,
    )?;
    if ctl.gen_opt.require_unbroken_ok {
        return Ok(EncloneIntermediates::default());
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
    ctl.perf_stats(&tpop, "populating features");

    // Test for no data.

    let tproto = Instant::now();
    if ctl.origin_info.n() == 0 {
        return Err(format!("\nNo TCR or BCR data have been specified.\n"));
    }

    // Search for SHM indels.

    search_for_shm_indels(&ctl, &tig_bc);
    if ctl.gen_opt.indels {
        return Ok(EncloneIntermediates::default());
    }

    // Record fate of non-cells.

    let mut fate = vec![HashMap::<String, String>::new(); vdj_cells.len()];
    if ctl.gen_opt.ncell {
        for i in 0..tig_bc.len() {
            let bc = &tig_bc[i][0].barcode;
            let li = tig_bc[i][0].dataset_index;
            if !bin_member(&vdj_cells[li], bc) {
                fate[li].insert(bc.clone(), "fails CELL filter".to_string());
            }
        }
    }

    // Filter using light --> heavy graph.

    graph_filter(&ctl, &mut tig_bc, ctl.gen_opt.graph, &mut fate);

    // Sort tig_bc.

    sort_tig_bc(&ctl, &mut tig_bc, &refdata);

    // Cross filter.

    cross_filter(&ctl, &mut tig_bc, &mut fate);

    // Look for barcode reuse.

    check_for_barcode_reuse(&ctl, &tig_bc)?;
    ctl.perf_stats(&tproto, "in proto stuff");

    // Find exact subclonotypes.

    let mut exact_clonotypes = find_exact_subclonotypes(&ctl, &tig_bc, &refdata, &mut fate);
    if ctl.gen_opt.utr_con || ctl.gen_opt.con_con {
        return Ok(EncloneIntermediates::default());
    }

    // Test for consistency between VDJ cells and GEX cells.

    test_vdj_gex_inconsistent(&ctl, &tig_bc, &exact_clonotypes, &vdj_cells, &gex_info)?;

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
    if ctl.clono_filt_opt_def.weak_foursies {
        erase_if(&mut exact_clonotypes, &to_delete);
    }

    // Filter if MAX_HEAVIES = 1 set.

    if ctl.gen_opt.max_heavies == 1 {
        let mut to_delete = vec![false; exact_clonotypes.len()];
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            let mut heavies = 0;
            for j in 0..ex.share.len() {
                if ex.share[j].left {
                    heavies += 1;
                }
            }
            if heavies > 1 {
                to_delete[i] = true;
            }
        }
        erase_if(&mut exact_clonotypes, &to_delete);
    }
    ctl.perf_stats(&t, "filtering foursies");

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

    let tbc = Instant::now();
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
    ctl.perf_stats(&tbc, "computing to_bc");

    // Make stirling ratio table.  Not sure that fixing the size of this is safe.

    let tsr = Instant::now();
    let sr = stirling2_ratio_table::<f64>(3000);
    ctl.perf_stats(&tsr, "computing stirling number table");

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
    disintegrate_onesies(
        &ctl,
        &mut disintegrated,
        &mut eq,
        &mut exact_clonotypes,
        &mut info,
        &mut join_info,
        &mut raw_joins,
    );

    // Update to_bc.

    let txxx = Instant::now();
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
    if ctl.gen_opt.heavy_chain_reuse {
        return Ok(EncloneIntermediates::default());
    }
    ctl.perf_stats(&txxx, "in some odds and ends");

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
                        if !ctl.clono_filt_opt_def.ngex {
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

    filter_by_fcell(&ctl, &mut orbits, &info, &mut exact_clonotypes);
    ctl.perf_stats(&tumi, "umi filtering and such");

    // Run some filters.

    some_filters(
        &mut orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
        &eq,
        &disintegrated,
        &mut fate,
    );

    // Mark VDJ noncells.

    let tmark = Instant::now();
    if ctl.clono_filt_opt_def.non_cell_mark {
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
    ctl.perf_stats(&tmark, "marking vdj noncells");
    Ok(EncloneIntermediates {
        setup: setup,
        ex: EncloneExacts {
            to_bc: to_bc,
            exact_clonotypes: exact_clonotypes,
            raw_joins: raw_joins,
            info: info.to_vec(),
            orbits: orbits,
            vdj_cells: vdj_cells,
            join_info: join_info,
            drefs: drefs,
            sr: sr,
            fate: fate,
            is_bcr: is_bcr,
        },
    })
}
