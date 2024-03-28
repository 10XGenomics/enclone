// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::opt_d_val::make_opt_d_val;
use crate::subset::subset_json;
use enclone_core::defs::ColInfo;
use enclone_core::enclone_structs::*;
use enclone_print::print_clonotypes::{print_clonotypes, PrintClonotypesResult};
use enclone_tail::grouper::grouper;
use enclone_tail::tail::tail_code;
use io_utils::{dir_list, open_for_read, path_exists};
use rayon::prelude::*;
use stats_utils::percent_ratio;
use std::{collections::HashMap, env, io::BufRead, thread, time};
use string_utils::TextUtils;
use vector_utils::*;

use hdf5::Reader;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone_stop(mut inter: EncloneIntermediates) -> Result<EncloneState, String> {
    // Unpack inputs.

    let to_bc = &inter.ex.to_bc;
    let exact_clonotypes = &inter.ex.exact_clonotypes;
    let raw_joins = &inter.ex.raw_joins;
    let info = &inter.ex.info;
    let orbits = &inter.ex.orbits;
    let vdj_cells = &inter.ex.vdj_cells;
    let refdata = &inter.setup.refdata;
    let join_info = &inter.ex.join_info;
    let drefs = &inter.ex.drefs;
    let gex_info = &inter.setup.gex_info;
    let sr = &inter.ex.sr;
    let ann = &inter.setup.ann;
    let fate = &mut inter.ex.fate;
    let ctl = &inter.setup.ctl;
    let is_bcr = inter.ex.is_bcr;
    let tall = &inter.setup.tall.unwrap();
    let allele_data = &inter.ex.allele_data;

    // Load the GEX and FB data.  This is quite horrible: the code and computation are duplicated
    // verbatim in fcell.rs.

    let mut d_readers = Vec::<Option<Reader>>::new();
    let mut ind_readers = Vec::<Option<Reader>>::new();
    for li in 0..ctl.origin_info.n() {
        if !ctl.origin_info.gex_path[li].is_empty() {
            let x = gex_info.h5_data[li].as_ref();
            if x.is_none() {
                // THIS FAILS SPORADICALLY, OBSERVED MULTIPLE TIMES,
                // CAUSING PUSH TO D_READERS BELOW TO FAIL.
                eprintln!("\nWeird, gex_info.h5_data[li].as_ref() is None.");
                eprintln!("Path = {}.", ctl.origin_info.gex_path[li]);
                let current = env::current_dir().unwrap();
                println!(
                    "The current working directory is {}",
                    current.canonicalize().unwrap().display()
                );
                if path_exists(&ctl.origin_info.gex_path[li]) {
                    eprintln!(
                        "The directory that is supposed to contain \
                        raw_feature_bc_matrix.h5 exists."
                    );
                    let list = dir_list(&ctl.origin_info.gex_path[li]);
                    eprintln!(
                        "This directory is {} and its contents are:",
                        ctl.origin_info.gex_path[li]
                    );
                    for i in 0..list.len() {
                        eprintln!("{}.  {}", i + 1, list[i]);
                    }
                    let h5_path =
                        format!("{}/raw_feature_bc_matrix.h5", ctl.origin_info.gex_path[li]);
                    eprintln!("H5 path = {}.", h5_path);
                    if !path_exists(&h5_path) {
                        let mut msg = format!("H5 path {} does not exist.\n", h5_path);
                        msg += "Retrying a few times to see if it appears.\n";
                        for _ in 0..5 {
                            msg += "Sleeping for 0.1 seconds.";
                            thread::sleep(time::Duration::from_millis(100));
                            if !path_exists(&h5_path) {
                                msg += "Now h5 path does not exist.\n";
                            } else {
                                msg += "Now h5 path exists.\n";
                                break;
                            }
                        }
                        msg += "Aborting.\n";
                        return Err(msg);
                    } else {
                        println!("h5 path exists.");
                    }
                } else {
                    println!("Path exists.");
                }
                println!();
            }
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
        if !ctl.origin_info.gex_path[li].is_empty() && ctl.gen_opt.h5_pre {
            res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
            res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
        }
    });

    // Find and print clonotypes.  (But we don't actually print them here.)
    if !ctl.gen_opt.trace_barcode.is_empty() {
        for u in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[u];
            for j in 0..ex.clones.len() {
                if ex.clones[j][0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    let PrintClonotypesResult {
        mut pics,
        mut exacts,
        in_center,
        mut rsi,
        mut out_datas,
        gene_scan_result,
    } = print_clonotypes(
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
        fate,
        allele_data,
    )?;

    // Gather some data for gene scan.
    let (mut tests, mut controls) = (vec![], vec![]);
    if ctl.gen_opt.gene_scan.is_some() {
        if !ctl.gen_opt.gene_scan_exact {
            for (i, in_sets) in gene_scan_result.iter().enumerate() {
                for in_set in in_sets {
                    if in_set.test {
                        tests.push(i);
                    }
                    if in_set.control {
                        controls.push(i);
                    }
                }
            }
        } else {
            for (in_sets, e) in gene_scan_result.iter().zip(exacts.iter()) {
                for (&ej, in_set) in e.iter().zip(in_sets) {
                    if in_set.test {
                        tests.push(ej);
                    }
                    if in_set.control {
                        controls.push(ej);
                    }
                }
            }
        }
    }

    // Process the SUBSET_JSON option.

    subset_json(ctl, exact_clonotypes, &exacts, ann)?;

    // Assign a D segment to each "left" column in a clonotype (if we need this information).
    // The assignments are to exact subclonotypes, and might differ across a clonotype, even
    // though the true values have to be the same.  This is also true for V and J segments,
    // although they are less likely to vary.

    let mut opt_d_val = Vec::<(usize, Vec<Vec<Vec<usize>>>)>::new();
    make_opt_d_val(
        ctl,
        exact_clonotypes,
        &exacts,
        &rsi,
        refdata,
        drefs,
        &mut opt_d_val,
    );

    // Group clonotypes.

    let mut groups = grouper(
        refdata,
        &exacts,
        &in_center,
        exact_clonotypes,
        ctl,
        &rsi,
        &opt_d_val,
        drefs,
    );

    // Remove clonotypes that are not in groups.

    let mut to_delete = vec![true; exacts.len()];
    for i in 0..groups.len() {
        for j in 0..groups[i].len() {
            to_delete[groups[i][j].0 as usize] = false;
        }
    }
    let mut count = 0;
    let mut to_new = HashMap::<usize, usize>::new();
    for i in 0..exacts.len() {
        if !to_delete[i] {
            to_new.insert(i, count);
            count += 1;
        }
    }
    for i in 0..groups.len() {
        for j in 0..groups[i].len() {
            groups[i][j].0 = to_new[&(groups[i][j].0 as usize)] as i32;
        }
    }
    erase_if(&mut exacts, &to_delete);
    erase_if(&mut out_datas, &to_delete);
    erase_if(&mut rsi, &to_delete);
    erase_if(&mut pics, &to_delete);

    // Tail code.

    let mut svgs = Vec::<String>::new();
    let mut group_pics = Vec::<String>::new();
    let mut last_widths = Vec::<u32>::new();
    let mut summary = String::new();
    tail_code(
        tall,
        refdata,
        &pics,
        &mut group_pics,
        &mut last_widths,
        &exacts,
        &rsi,
        exact_clonotypes,
        ctl,
        &mut out_datas,
        join_info,
        gex_info,
        vdj_cells,
        fate,
        &tests,
        &controls,
        &h5_data,
        &d_readers,
        &ind_readers,
        drefs,
        &groups,
        &opt_d_val,
        &mut svgs,
        &mut summary,
    )?;

    let (mut cpu_all_stop, mut cpu_this_stop) = (0, 0);
    if ctl.gen_opt.print_cpu || ctl.gen_opt.print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        if let Some(line) = f.lines().next() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_stop = t.before(" ").force_usize();
        }
        let f = open_for_read![&format!("/proc/{}/stat", std::process::id())];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(' ').collect::<Vec<&str>>();
            cpu_this_stop = fields[13].force_usize();
        }
        let (this_used, all_used) = (
            cpu_this_stop - ctl.gen_opt.cpu_this_start,
            cpu_all_stop - ctl.gen_opt.cpu_all_start,
        );
        if ctl.gen_opt.print_cpu {
            println!("{}", this_used);
        } else {
            println!(
                "used cpu = {} = {:.1}% of total",
                this_used,
                percent_ratio(this_used, all_used)
            );
        }
    }

    if !(ctl.gen_opt.noprint && ctl.parseable_opt.pout == "stdout") && !ctl.gen_opt.no_newline {
        println!();
    }
    let outs = MainEncloneOutput {
        pics: group_pics,
        last_widths,
        svgs,
        summary,
        metrics: gex_info.metrics.clone(),
        dataset_names: ctl.origin_info.dataset_id.clone(),
        parseable_stdouth: ctl.parseable_opt.pout == "stdouth",
        noprint: ctl.gen_opt.noprint,
        noprintx: ctl.gen_opt.noprintx,
        html: ctl.gen_opt.html,
        ngroup: ctl.clono_group_opt.ngroup,
        pretty: ctl.pretty,
    };
    Ok(EncloneState { inter, outs })
}
