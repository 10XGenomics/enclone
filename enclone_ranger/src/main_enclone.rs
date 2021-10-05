// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use self::refx::{make_vdj_ref_data_core, RefData};
use crate::setup::{critical_args, setup};
use crate::stop::main_enclone_stop;
use enclone::innate::species;
use enclone_args::load_gex::get_gex_info;
use enclone_core::defs::EncloneControl;
use enclone_core::enclone_structs::*;
use enclone_stuff::start::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use vdj_ann::refx;

pub fn main_enclone(args: &Vec<String>) -> Result<(), String> {
    let setup = main_enclone_setup(args)?;
    let inter = main_enclone_start(setup)?;
    main_enclone_stop(inter)
}

pub fn main_enclone_setup(args: &Vec<String>) -> Result<EncloneSetup, String> {
    let tall = Instant::now();

    // Set up stuff, read args, etc.

    let args_orig = args.clone();
    let mut ctl = EncloneControl::default();
    let args = critical_args(args, &mut ctl)?;
    ctl.start_time = Some(tall);
    ctl.gen_opt.cpu_all_start = 0;
    ctl.gen_opt.cpu_this_start = 0;
    let mut argsx = Vec::<String>::new();
    setup(&mut ctl, &args, &mut argsx)?;
    let mut argsy = Vec::<String>::new();
    for i in 0..args_orig.len() {
        if args_orig[i] != "HTML"
            && args_orig[i] != "NOPAGER"
            && args_orig[i] != "FORCE_EXTERNAL"
            && !args_orig[i].starts_with("PRE=")
            && !args_orig[i].starts_with("MAX_CORES=")
        {
            argsy.push(args_orig[i].clone());
        }
    }
    if argsy.len() == 1 || (argsy.len() > 1 && (argsy[1] == "help" || argsy[1] == "--help")) {
        return Ok(EncloneSetup::default());
    }

    // Get gene expression and feature barcode counts.  

    let gex_info = get_gex_info(&mut ctl)?;

    // Determine the reference sequence that is to be used.

    let mut refx = String::new();
    let ann = "contig_annotations.json";
    let fx = File::open(&ctl.gen_opt.refname);
    let f = BufReader::new(fx.unwrap());
    for line in f.lines() {
        let s = line.unwrap();
        refx += &s;
        refx += "\n";
    }

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
    make_vdj_ref_data_core(&mut refdata, refx2, &ext_refx, is_tcr, is_bcr, None);
    let mut to_ref_index = HashMap::<usize, usize>::new();
    for i in 0..refdata.refs.len() {
        to_ref_index.insert(refdata.id[i] as usize, i);
    }

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata);

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
