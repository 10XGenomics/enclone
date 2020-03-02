// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use vdj_ann::*;

use self::refx::*;
use crate::allele::*;
use crate::defs::*;
use crate::explore::*;
use crate::graph_filter::*;
use crate::info::*;
use crate::join::*;
use crate::load_gex::*;
use crate::loupe::*;
use crate::misc1::*;
use crate::misc2::*;
use crate::misc3::*;
use crate::print_clonotypes::*;
use crate::proc_args2::*;
use crate::read_json::*;
use debruijn::dna_string::*;
use equiv::EquivRel;
use io_utils::*;
use perf_stats::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    time::Instant,
};
use vector_utils::*;

pub fn main_enclone(args: &Vec<String>) {
    // Set up stuff, read args, etc.

    let tall = Instant::now();
    let mut ctl = EncloneControl::default();
    setup(&mut ctl, &args);

    // Read external data.

    if ctl.gen_opt.ext.len() > 0 {
        let f = open_for_read![&ctl.gen_opt.ext];
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

    // Get gene expression and antibody counts.

    let gex_info = get_gex_info(&mut ctl);

    // Build reference data.

    let tr = Instant::now();
    let mut refdata = RefData::new();
    let mut refx = String::new();
    if ctl.gen_opt.refname.len() > 0 {
        let fx = File::open(&ctl.gen_opt.refname);
        if fx.is_err() {
            eprintln!(
                "\nProblem with REF: unable to read from the file\n\
                 \"{}\".\nPlease check that that path makes sense and that you have read \
                 permission along that path.\n",
                ctl.gen_opt.refname
            );
            std::process::exit(1);
        }
        let f = BufReader::new(fx.unwrap());
        let mut nheader = 0;
        let mut bases = 0;
        let mut na = 0;
        let mut nc = 0;
        let mut ng = 0;
        let mut nt = 0;
        for line in f.lines() {
            let s = line.unwrap();
            refx += &s;
            refx += &"\n";
            if s.starts_with('>') {
                nheader += 1;
            } else {
                for c in s.chars() {
                    bases += 1;
                    if c == 'A' || c == 'a' {
                        na += 1;
                    } else if c == 'C' || c == 'c' {
                        nc += 1;
                    } else if c == 'G' || c == 'g' {
                        ng += 1;
                    } else if c == 'T' || c == 't' {
                        nt += 1;
                    }
                }
            }
        }
        if nheader == 0 || bases == 0 || (na + nc + ng + nt) as f64 / (bases as f64) < 0.95 {
            eprintln!("\nProblem with REF: it is not a FASTA file.\n");
            std::process::exit(1);
        }
    } else if ctl.gen_opt.mouse {
        refx = mouse_ref();
    } else {
        refx = human_ref();
    }
    let ext_refx = String::new();
    let (mut is_tcr, mut is_bcr) = (true, true);
    if ctl.gen_opt.tcr {
        is_bcr = false;
    }
    if ctl.gen_opt.bcr {
        is_tcr = false;
    }
    make_vdj_ref_data_core(&mut refdata, &refx, &ext_refx, is_tcr, is_bcr, None);
    let mut to_ref_index = HashMap::<usize, usize>::new();
    for i in 0..refdata.refs.len() {
        to_ref_index.insert(refdata.id[i] as usize, i);
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds building reference, peak mem = {:.2} GB",
            elapsed(&tr),
            peak_mem_usage_gb()
        );
    }
    if ctl.comp {
        println!(
            "used {:.2} seconds up through building reference",
            elapsed(&tall)
        );
    }

    // Parse the json annotations file.

    let mut tig_bc = Vec::<Vec<TigData>>::new();
    parse_json_annotations_files(&ctl, &mut tig_bc, &refdata, &to_ref_index);

    // Search for SHM indels.  Exploratory.

    let tproto = Instant::now();
    search_for_shm_indels(&ctl, &tig_bc);

    // Filter using light --> heavy graph.

    if !ctl.gen_opt.ngraph_filter {
        graph_filter(&mut tig_bc, ctl.gen_opt.graph);
    }

    // Sort tig_bc.

    sort_tig_bc(&mut tig_bc, &refdata);

    // Cross filter.

    cross_filter(&ctl, &mut tig_bc);

    // Look for barcode reuse.

    check_for_barcode_reuse(&ctl, &tig_bc);
    if ctl.comp {
        println!("used {:.2} seconds in proto stuff", elapsed(&tproto));
    }

    // Find exact subclonotypes.

    let mut exact_clonotypes = find_exact_subclonotypes(&ctl, &tig_bc, &refdata);

    // Filter out some foursie artifacts.

    if ctl.clono_filt_opt.weak_foursies {
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
                            }
                        }
                    }
                }
            }
        }
        erase_if(&mut exact_clonotypes, &to_delete);
        if ctl.comp {
            println!("used {:.2} seconds filtering foursies", elapsed(&t));
        }
    }

    // Look for insertions (experimental).

    find_insertions(&ctl, &exact_clonotypes);

    // Build info about clonotypes.  Note that this edits the V reference sequence to perform
    // an indel in some cases.

    let mut info: Vec<CloneInfo> = build_info(&refdata, &ctl, &mut exact_clonotypes);

    // Derive consensus sequences for alternate alleles of V segments.  Then create donor
    // reference sequences for Loupe.

    let alt_refs : Vec<(usize,usize,DnaString)>  // {(donor, ref id, alt seq)}
        = find_alleles( &refdata, &ctl, &exact_clonotypes );
    if ctl.gen_opt.dref_file.len() > 0 {
        let f = File::create(&ctl.gen_opt.dref_file);
        if f.is_err() {
            eprintln!(
                "\nProblem with DONOR_REF_FILE: unable to write to the file\n\
                 \"{}\".\nPlease check that that path makes sense and that you have write \
                 permission for it.\n",
                ctl.gen_opt.dref_file
            );
            std::process::exit(1);
        }
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
                ">{}.{}.{}, {}\n{}",
                refdata.id[ref_id],
                donor + 1,
                count + 1,
                refdata.name[ref_id],
                alt_seq.to_string()
            );
            count += 1;
        }
    }
    let tdonor = Instant::now();
    let drefs = make_donor_refs(&alt_refs, &refdata);
    if ctl.comp {
        println!("\nused {:.2} seconds making donor refs", elapsed(&tdonor));
    }

    // Update reference sequences for V segments by substituting in alt alleles if better.

    sub_alts(&ctl, &alt_refs, &mut info, &mut exact_clonotypes);

    // Form equivalence relation on exact subclonotypes.

    let eq: EquivRel = join_exacts(is_bcr, &refdata, &ctl, &exact_clonotypes, &info);
    if ctl.comp {
        if ctl.clono_filt_opt.ncells_low < ctl.clono_filt_opt.ncells_high {
            println!("");
        }
    }

    // Lookup for heavy chain reuse (special purpose experimental option).

    lookup_heavy_chain_reuse(&ctl, &exact_clonotypes, &info, &eq);

    // Find and print clonotypes.

    let torb = Instant::now();
    print_clonotypes(
        &refdata,
        &drefs,
        &ctl,
        &exact_clonotypes,
        &info,
        &eq,
        &gex_info,
    );
    if ctl.comp {
        println!("\nused {:.2} seconds making orbits", elapsed(&torb));
    }

    // Report total cells.

    let mut ncells = 0;
    for i in 0..info.len() {
        ncells += exact_clonotypes[info[i].clonotype_index].ncells();
    }
    if !ctl.silent {
        if ctl.clono_filt_opt.ncells_low < ctl.clono_filt_opt.ncells_high {
            println!("");
        }
        println!("total cells in clonotypes = {}", ncells);
    }

    // Report computational performance.

    if ctl.comp {
        println!(
            "{:.2} seconds total, peak mem = {:.2} GB",
            elapsed(&tall),
            peak_mem_usage_gb()
        );
    }
    println!("");
    // It's not totally clear that the exit below actually saves time.  Would need more testing.
    if !ctl.gen_opt.cellranger {
        std::process::exit(0);
    }
}
