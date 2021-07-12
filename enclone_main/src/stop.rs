// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::main_enclone::*;
use crate::opt_d_val::*;
use crate::subset::*;
use enclone_core::defs::*;
use enclone_print::print_clonotypes::*;
use enclone_proto::types::DonorReferenceItem;
use enclone_tail::grouper::*;
use enclone_tail::tail::tail_code;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use stats_utils::*;
use std::{
    collections::HashMap,
    env,
    fs::File,
    io::{BufRead, BufReader},
    thread, time,
    time::Instant,
};
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
pub struct EncloneIntermediates {
    pub setup: EncloneSetup,
    pub ex: EncloneExacts,
}

#[derive(Default, Clone)]
pub struct EncloneExacts {
    pub to_bc: HashMap<(usize, usize), Vec<String>>,
    pub exact_clonotypes: Vec<ExactClonotype>,
    pub raw_joins: Vec<Vec<usize>>,
    pub info: Vec<CloneInfo>,
    pub orbits: Vec<Vec<i32>>,
    pub vdj_cells: Vec<Vec<String>>,
    pub join_info: Vec<(usize, usize, bool, Vec<u8>)>,
    pub drefs: Vec<DonorReferenceItem>,
    pub sr: Vec<Vec<f64>>,
    pub fate: Vec<HashMap<String, String>>, // GETS MODIFIED SUBSEQUENTLY
    pub is_bcr: bool,
}

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
    let mut fate = &mut inter.ex.fate;
    let ctl = &inter.setup.ctl;
    let is_bcr = inter.ex.is_bcr;
    let tall = &inter.setup.tall.unwrap();

    // Load the GEX and FB data.

    let tdi = Instant::now();
    let mut d_readers = Vec::<Option<hdf5::Reader>>::new();
    let mut ind_readers = Vec::<Option<hdf5::Reader>>::new();
    for li in 0..ctl.origin_info.n() {
        if ctl.origin_info.gex_path[li].len() > 0 && !gex_info.gex_matrices[li].initialized() {
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
                println!("");
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
        if ctl.origin_info.gex_path[li].len() > 0
            && !gex_info.gex_matrices[li].initialized()
            && ctl.gen_opt.h5_pre
        {
            res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
            res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
        }
    });
    ctl.perf_stats(&tdi, "setting up readers");

    // Find and print clonotypes.  (But we don't actually print them here.)

    let torb = Instant::now();
    let mut pics = Vec::<String>::new();
    let mut exacts = Vec::<Vec<usize>>::new(); // ugly reuse of name
    let mut in_center = Vec::<bool>::new();
    let mut rsi = Vec::<ColInfo>::new(); // ditto
    let mut out_datas = Vec::<Vec<HashMap<String, String>>>::new();
    let mut tests = Vec::<usize>::new();
    let mut controls = Vec::<usize>::new();
    print_clonotypes(
        is_bcr,
        &to_bc,
        &sr,
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &orbits,
        &raw_joins,
        &gex_info,
        &vdj_cells,
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
        &mut fate,
    )?;

    // Lock data structures so they can't be changed accidentally.

    let ctl = ctl;
    let refdata = refdata;
    let exact_clonotypes = exact_clonotypes;
    let exacts = exacts;
    let pics = pics;

    // Process the SUBSET_JSON option.

    subset_json(&ctl, &exact_clonotypes, &exacts, &ann);
    ctl.perf_stats(&torb, "making orbits");

    // Assign a D segment to each "left" column in a clonotype (if we need this information).
    // The assignments are to exact subclonotypes, and might differ across a clonotype, even
    // though the true values have to be the same.  This is also true for V and J segments,
    // although they are less likely to vary.

    let mut opt_d_val = Vec::<(usize, Vec<Vec<Vec<usize>>>)>::new();
    make_opt_d_val(
        &ctl,
        &exact_clonotypes,
        &exacts,
        &rsi,
        &refdata,
        &drefs,
        &mut opt_d_val,
    );

    // Group clonotypes.
    let t = Instant::now();
    let groups = grouper(
        &refdata,
        &exacts,
        &in_center,
        &exact_clonotypes,
        &ctl,
        &rsi,
        &opt_d_val,
    );
    ctl.perf_stats(&t, "in grouper");
    // Process TOY_COM option.

    if ctl.gen_opt.toy_com {
        println!(
            "\nHello, enclone is now in server mode.  Hopefully you have already started\n\
            enclone_client in a separate terminal window, before starting enclone, because\n\
            otherwise the system won't work.  The client should now show a prompt.\n"
        );
        /*
        enclone_server(&ctl, &refdata, &exacts, &exact_clonotypes, &groups, &pics)
            .await
            .unwrap();
        */
        return Ok(EncloneState::default());
    }

    // Tail code.

    let mut svgs = Vec::<String>::new();
    let mut group_pics = Vec::<String>::new();
    let mut last_widths = Vec::<usize>::new();
    let mut summary = String::new();
    tail_code(
        &tall,
        &refdata,
        &pics,
        &mut group_pics,
        &mut last_widths,
        &exacts,
        &rsi,
        &exact_clonotypes,
        &ctl,
        &mut out_datas,
        &join_info,
        &gex_info,
        &vdj_cells,
        &fate,
        &tests,
        &controls,
        &h5_data,
        &d_readers,
        &ind_readers,
        &drefs,
        &groups,
        &opt_d_val,
        &mut svgs,
        &mut summary,
    )?;

    // Report profiling.

    if ctl.gen_opt.profile {
        let t = Instant::now();
        stop_profiling();
        ctl.perf_stats(&t, "summarizing profiling");
    }

    // Report computational performance.

    let delta;
    unsafe {
        delta = elapsed(&tall) - WALLCLOCK;
    }
    let deltas = format!("{:.2}", delta);
    ctl.perf_stats(&tall, "total");
    if ctl.perf_opt.comp {
        println!("used {} seconds unaccounted for", deltas);
        println!("peak mem usage = {:.1} MB", peak_mem_usage_gb() * 1000.0);
    }
    if ctl.perf_opt.comp_enforce {
        if deltas.force_f64() > 0.03 {
            return Err(format!(
                "\nUnaccounted time = {} seconds, but COMPE option required that it \
                be at most 0.03.\n\n\
                Note that this may fail for a small fraction of runs, even though \
                nothing is wrong.\n",
                deltas
            ));
        }
    }
    let (mut cpu_all_stop, mut cpu_this_stop) = (0, 0);
    if ctl.gen_opt.print_cpu || ctl.gen_opt.print_cpu_info {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_stop = t.before(" ").force_usize();
            break;
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

    if !(ctl.gen_opt.noprint && ctl.parseable_opt.pout == "stdout") {
        println!("");
    }
    let outs = MainEncloneOutput {
        pics: group_pics,
        last_widths: last_widths,
        svgs: svgs,
        summary: summary,
        noprint: ctl.gen_opt.noprint,
        noprintx: ctl.gen_opt.noprintx,
        html: ctl.gen_opt.html,
        ngroup: ctl.clono_group_opt.ngroup,
        pretty: ctl.pretty,
    };
    Ok(EncloneState {
        inter: inter,
        outs: outs,
    })
}
