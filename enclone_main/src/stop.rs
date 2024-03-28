// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::opt_d_val::make_opt_d_val;
use crate::subset::subset_json;

use enclone_core::enclone_structs::*;
use enclone_process::process_clonotypes::process_clonotypes;
use enclone_tail::grouper::grouper;
use enclone_tail::print_clonotypes::{EncloneOrbitProcessor, PrintClonotypesResult};
use enclone_tail::tail::tail_code;
use io_utils::open_for_read;
use stats_utils::percent_ratio;
use std::{collections::HashMap, io::BufRead};
use string_utils::TextUtils;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone_stop(
    setup: &EncloneSetup,
    exacts: &EncloneExacts,
    fate: Vec<BarcodeFates>,
) -> Result<(), String> {
    // Unpack inputs.

    let exact_clonotypes = &exacts.exact_clonotypes;
    let vdj_cells = &exacts.vdj_cells;
    let refdata = &setup.refdata;
    let join_info = &exacts.join_info;
    let drefs = &exacts.drefs;
    let gex_info = &setup.gex_info;
    let ann = &setup.ann;
    let ctl = &setup.ctl;
    let tall = &setup.tall.unwrap();

    let gex_readers = setup.create_gex_readers();

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

    let mut proc = EncloneOrbitProcessor::new(setup, &exacts.vdj_cells);

    process_clonotypes(setup, exacts, &gex_readers, &fate, &mut proc)?;

    let PrintClonotypesResult {
        mut pics,
        mut exacts,
        in_center,
        mut rsi,
        mut out_datas,
        gene_scan_result,
    } = proc.result;

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
        &fate,
        &tests,
        &controls,
        &gex_readers,
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
    Ok(())
}
