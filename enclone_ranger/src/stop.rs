// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::ColInfo;
use enclone_core::enclone_structs::*;
use enclone_print::print_clonotypes::print_clonotypes;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn main_enclone_stop_ranger(mut inter: EncloneIntermediates) -> Result<(), String> {
    // Unpack inputs.

    let to_bc = &inter.ex.to_bc;
    let exact_clonotypes = &inter.ex.exact_clonotypes;
    let raw_joins = &inter.ex.raw_joins;
    let info = &inter.ex.info;
    let orbits = &inter.ex.orbits;
    let vdj_cells = &inter.ex.vdj_cells;
    let refdata = &inter.setup.refdata;
    let drefs = &inter.ex.drefs;
    let gex_info = &inter.setup.gex_info;
    let sr = &inter.ex.sr;
    let fate = &mut inter.ex.fate;
    let ctl = &inter.setup.ctl;
    let is_bcr = inter.ex.is_bcr;

    // Load the GEX and FB data.  This is quite horrible: the code and computation are duplicated
    // verbatim in fcell.rs.

    let mut d_readers = Vec::<Option<hdf5::Reader>>::new();
    let mut ind_readers = Vec::<Option<hdf5::Reader>>::new();
    for li in 0..ctl.origin_info.n() {
        if !ctl.origin_info.gex_path[li].is_empty() && !gex_info.gex_matrices[li].initialized() {
            let x = gex_info.h5_data[li].as_ref();
            d_readers.push(Some(x.unwrap().as_reader()));
            ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
        } else {
            d_readers.push(None);
            ind_readers.push(None);
        }
    }
    let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
    for li in 0..ctl.origin_info.n() {
        h5_data.push((li, Vec::new(), Vec::new()));
    }
    h5_data.par_iter_mut().for_each(|res| {
        let li = res.0;
        if !ctl.origin_info.gex_path[li].is_empty()
            && !gex_info.gex_matrices[li].initialized()
            && ctl.gen_opt.h5_pre
        {
            res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
            res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
        }
    });

    // Find and print clonotypes.  (But we don't actually print them here.)

    let mut pics = Vec::<String>::new();
    let mut exacts = Vec::<Vec<usize>>::new(); // ugly reuse of name
    let mut in_center = Vec::<bool>::new();
    let mut rsi = Vec::<ColInfo>::new(); // ditto
    let mut out_datas = Vec::<Vec<HashMap<String, String>>>::new();
    let mut tests = Vec::<usize>::new();
    let mut controls = Vec::<usize>::new();
    print_clonotypes(
        is_bcr,
        to_bc,
        sr,
        refdata,
        drefs,
        ctl,
        exact_clonotypes,
        info,
        orbits,
        raw_joins,
        gex_info,
        vdj_cells,
        &d_readers,
        &ind_readers,
        &h5_data,
        &mut pics,
        &mut exacts,
        &mut in_center,
        &mut rsi,
        &mut out_datas,
        &mut tests,
        &mut controls,
        fate,
    )?;
    Ok(())
}
