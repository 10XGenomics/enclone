// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::define_mat::*;
use crate::filter::*;
use crate::finish_table::*;
use crate::gene_scan::*;
use crate::loupe::*;
use crate::print_utils1::*;
use crate::print_utils2::*;
use crate::print_utils3::*;
use crate::print_utils4::*;
use crate::print_utils5::*;
use enclone_core::allowed_vars::*;
use enclone_core::defs::*;
use enclone_core::mammalian_fixed_len::*;
use enclone_proto::types::*;
use equiv::EquivRel;
use rayon::prelude::*;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// Print clonotypes.  A key challenge here is to define the columns that represent shared
// chains.  This is given below by the code that forms an equivalence relation on the CDR3_AAs.
//
// This code carries out a second function, which is to filter out exact subclonotypes in orbits
// that appear to be junk.  Exactly how these should be reflected in files is TBD.
//
// Some inputs for this section:
// refdata                = reference sequence data
// ctl                    = control parameters
// exact_clonotypes       = the vector of all exact subclonotypes
// info                   = vector of clonotype info
// eq                     = equivalence relation on info

pub fn print_clonotypes(
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    orbits: &Vec<Vec<i32>>,
    raw_joins: &Vec<Vec<usize>>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    d_readers: &Vec<Option<hdf5::Reader>>,
    ind_readers: &Vec<Option<hdf5::Reader>>,
    h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
    pics: &mut Vec<String>,
    exacts: &mut Vec<Vec<usize>>,
    in_center: &mut Vec<bool>,
    rsi: &mut Vec<ColInfo>,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    tests: &mut Vec<usize>,
    controls: &mut Vec<usize>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    let lvars = &ctl.clono_print_opt.lvars;

    // Compute extra args.

    let extra_args = extra_args(&ctl);

    // Define parseable output columns.  The entire machinery for parseable output is controlled
    // by macros that begin with "speak".

    let mut parseable_fields = Vec::<String>::new();
    set_speakers(&ctl, &mut parseable_fields);
    let pcols_sort = &ctl.parseable_opt.pcols_sort;

    // Identify certain extra parseable variables.  These arise from parameterizable cvars.

    let mut extra_parseables = Vec::<String>::new();
    {
        let mut exclusions = ctl.clono_print_opt.cvars.clone();
        for v in CVARS_ALLOWED.iter() {
            exclusions.push(v.to_string());
        }
        for v in CVARS_ALLOWED_PCELL.iter() {
            exclusions.push(v.to_string());
        }
        unique_sort(&mut exclusions);
        for x in pcols_sort.iter() {
            let mut chars = Vec::<char>::new();
            for c in x.chars() {
                chars.push(c);
            }
            let n = chars.len();
            let mut trim = 0;
            for i in (0..n).rev() {
                if !chars[i].is_digit(10) {
                    break;
                }
                trim += 1;
            }
            if trim > 0 {
                let mut v = String::new();
                for i in 0..n - trim {
                    v.push(chars[i]);
                }
                if !bin_member(&exclusions, &v) {
                    extra_parseables.push(v);
                }
            }
        }
        unique_sort(&mut extra_parseables);
    }

    // Compute all_vars.

    let rsi_vars = &ctl.clono_print_opt.cvars;
    let mut all_vars = rsi_vars.clone();
    for j in 0..CVARS_ALLOWED.len() {
        let var = &CVARS_ALLOWED[j];
        if !rsi_vars.contains(&var.to_string()) {
            all_vars.push(var.to_string());
        }
    }
    for j in 0..CVARS_ALLOWED_PCELL.len() {
        let var = &CVARS_ALLOWED_PCELL[j];
        if !rsi_vars.contains(&var.to_string()) {
            all_vars.push(var.to_string());
        }
    }
    all_vars.append(&mut extra_parseables.clone());
    for x in extra_args.iter() {
        if !rsi_vars.contains(&x) {
            all_vars.push(x.clone());
        }
    }

    // Compute total cells.

    let mut total_cells = 0;
    for i in 0..exact_clonotypes.len() {
        total_cells += exact_clonotypes[i].ncells();
    }

    // Test for presence of GEX/FB data.

    let mut have_gex = false;
    for i in 0..ctl.origin_info.gex_path.len() {
        if ctl.origin_info.gex_path[i].len() > 0 {
            have_gex = true;
        }
    }

    // Gather alt_bcs_fields.

    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);

    // Compute number of vdj cells that are gex.

    let mut n_vdj_gex = Vec::<usize>::new();
    for li in 0..ctl.origin_info.n() {
        let mut n = 0;
        for y in gex_info.pca[li].iter() {
            if bin_member(&vdj_cells[li], &y.0) {
                n += 1;
            }
        }
        n_vdj_gex.push(n);
    }

    // Compute peer groups.

    let peer_groups = mammalian_fixed_len_peer_groups(&refdata);

    // Traverse the orbits.

    // 0: index in reps
    // 1: vector of clonotype pictures
    // 2: vector of some clonotype info
    //    [parallel to 1]
    // next to last three entries = whitelist contam, denominator for that, low gex count
    // added out_datas (used to be next to last three, now one more)
    let mut results = Vec::<(
        usize,
        Vec<String>,
        Vec<(Vec<usize>, ColInfo)>,
        usize,
        usize,
        usize,
        Vec<Clonotype>,
        Vec<Vec<HashMap<String, String>>>,
        isize,
        Vec<bool>,
        Vec<bool>,
        Vec<(usize, String, String)>,
        Vec<bool>,
    )>::new();
    for i in 0..orbits.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<(Vec<usize>, ColInfo)>::new(),
            0,
            0,
            0,
            Vec::<Clonotype>::new(),
            Vec::<Vec<HashMap<String, String>>>::new(),
            0,
            Vec::<bool>::new(),
            Vec::<bool>::new(),
            Vec::new(),
            Vec::new(),
        ));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let o = &orbits[i];
        let mut od = Vec::<(Vec<usize>, usize, i32)>::new();
        for id in o.iter() {
            let x: &CloneInfo = &info[*id as usize];
            od.push((x.origin.clone(), x.clonotype_id, *id));
        }
        od.sort();

        // Reconstruct the participating clones.  This is needed because most exact subclonotypes
        // having more than two chains have been split up.
        //
        // Capture these data into parallel data structures, one per exact subclonotype:
        // exacts: the exact subclonotype ids
        // mults:  number of cells [redundant, might remove]

        let mut exacts = Vec::<usize>::new();
        let mut mults = Vec::<usize>::new();
        let mut j = 0;
        let loupe_clonotypes = &mut res.6;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            let mut mult = 0 as usize;
            for l in j..k {
                let x: &CloneInfo = &info[od[l].2 as usize];
                let m = x.clonotype_index;
                mult = exact_clonotypes[m].clones.len();
            }
            mults.push(mult);
            exacts.push(od[j].1);
            j = k;
        }

        // There are two passes.  On the first pass we only identify the exact subclonotypes that
        // are junk.  On the second pass we remove those and then print the orbit.

        let mut bads = vec![false; exacts.len()];
        for pass in 1..=2 {
            // Delete weak exact subclonotypes.

            if pass == 2 && !ctl.clono_filt_opt.protect_bads {
                erase_if(&mut mults, &bads);
                erase_if(&mut exacts, &bads);
            }

            // Sort exact subclonotypes.

            let mat = define_mat(
                is_bcr,
                &to_bc,
                &sr,
                &ctl,
                &exact_clonotypes,
                &exacts,
                &od,
                &info,
                &raw_joins,
            );
            let mut priority = Vec::<(Vec<bool>, usize, usize)>::new();
            for u in 0..exacts.len() {
                let mut typex = vec![false; mat.len()];
                for col in 0..mat.len() {
                    if mat[col][u].is_some() {
                        typex[col] = true;
                    }
                }
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                let mut utot0 = 0;
                let mid = mat[0][u];
                if mid.is_some() {
                    let mid = mid.unwrap();
                    let ex = &exact_clonotypes[clonotype_id];
                    for j in 0..ex.clones.len() {
                        utot0 += ex.clones[j][mid].umi_count;
                    }
                }
                priority.push((typex.clone(), ex.ncells(), utot0));
            }
            let permutation = permutation::sort(&priority[..]);
            exacts = permutation.apply_slice(&exacts[..]);
            mults = permutation.apply_slice(&mults[..]);
            exacts.reverse();
            mults.reverse();

            // Define a matrix mat[col][ex] which is the column of the exact subclonotype
            // corresponding to the given column col of the clonotype, which may or may not be
            // defined.  Then define other information associated to each chain.  These are
            // reference sequence identifiers, CDR3 start positions, and the like.

            let nexacts = exacts.len();
            let mat = define_mat(
                is_bcr,
                &to_bc,
                &sr,
                &ctl,
                &exact_clonotypes,
                &exacts,
                &od,
                &info,
                &raw_joins,
            );
            let cols = mat.len();
            let mut rsi = define_column_info(&ctl, &exacts, &exact_clonotypes, &mat, &refdata);
            rsi.mat = mat;
            let mat = &rsi.mat;

            // Let n be the total number of cells in this pass.

            let n: usize = mults.iter().sum();

            // Filter.

            let mut in_center = true;
            if pass == 2
                && !survives_filter(&exacts, &rsi, &ctl, &exact_clonotypes, &refdata, &gex_info)
            {
                if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                    in_center = false;
                } else {
                    continue;
                }
            }

            // Generate Loupe data.

            if (ctl.gen_opt.binary.len() > 0 || ctl.gen_opt.proto.len() > 0) && pass == 2 {
                loupe_clonotypes.push(make_loupe_clonotype(
                    &exact_clonotypes,
                    &exacts,
                    &rsi,
                    &refdata,
                    &dref,
                ));
            }

            // Set up for parseable output.

            let mut out_data = Vec::<HashMap<String, String>>::new();

            // Print the orbit.
            // ◼ An assumption of this code is that a productive pair does not have two contigs
            // ◼ having identical CDR3_AA sequences.  At present this is not enforced by the
            // ◼ assembly stage, so the assumption is violated.  To work around this there are
            // ◼ some unsavory workarounds below.

            let mut mlog = Vec::<u8>::new();
            if n >= ctl.clono_filt_opt.ncells_low
                || ctl.clono_group_opt.asymmetric_center == "from_filters"
            {
                // Start to generate parseable output.

                if pass == 2 {
                    start_gen(
                        &ctl,
                        &exacts,
                        &exact_clonotypes,
                        &mut out_data,
                        &mut mlog,
                        &extra_args,
                    );
                }

                // Find variant positions.  And some other things.

                let mut vars = Vec::<Vec<usize>>::new();
                let mut vars_amino = Vec::<Vec<usize>>::new();
                let mut shares_amino = Vec::<Vec<usize>>::new();
                vars_and_shares(
                    pass,
                    &ctl,
                    &exacts,
                    &exact_clonotypes,
                    &rsi,
                    &refdata,
                    &dref,
                    &mut vars,
                    &mut vars_amino,
                    &mut shares_amino,
                    &mut out_data,
                );

                // Mark some weak exact subclonotypes for deletion.

                if pass == 1 {
                    delete_weaks(
                        &ctl,
                        &exacts,
                        &mults,
                        &exact_clonotypes,
                        total_cells,
                        &mat,
                        &refdata,
                        &vars,
                        &mut bads,
                        &mut res.11,
                    );
                }

                // Done unless on second pass.  Unless there are bounds or COMPLETE specified.

                if pass == 1 && ctl.clono_filt_opt.bounds.len() == 0 && !ctl.gen_opt.complete {
                    continue;
                }

                // Define amino acid positions to show.

                let show_aa = build_show_aa(
                    &ctl,
                    &rsi,
                    &vars_amino,
                    &shares_amino,
                    &refdata,
                    &dref,
                    &exacts,
                    &exact_clonotypes,
                );

                // Define field types corresponding to the amino acid positions to show.

                let mut field_types = vec![Vec::new(); cols];
                for cx in 0..cols {
                    let mut ft = vec![0 as u8; show_aa[cx].len()];
                    let cs1 = rsi.cdr1_starts[cx];
                    let cs2 = rsi.cdr2_starts[cx];
                    let cs3 = rsi.cdr3_starts[cx];
                    let n3 = rsi.cdr3_lens[cx];
                    let fs1 = rsi.fr1_starts[cx];
                    let fs2 = rsi.fr2_starts[cx];
                    let fs3 = rsi.fr3_starts[cx];
                    let show_cdr1 = cs1.is_some()
                        && fs2.is_some()
                        && cs1.unwrap() <= fs2.unwrap()
                        && ctl.clono_print_opt.amino.contains(&"cdr1".to_string());
                    let show_cdr2 = cs2.is_some()
                        && fs3.is_some()
                        && cs2.unwrap() <= fs3.unwrap()
                        && ctl.clono_print_opt.amino.contains(&"cdr2".to_string());
                    let show_cdr3 = ctl.clono_print_opt.amino.contains(&"cdr3".to_string());
                    let show_fwr1 = cs1.is_some()
                        && rsi.fr1_starts[cx] <= cs1.unwrap()
                        && ctl.clono_print_opt.amino.contains(&"fwr1".to_string());
                    let show_fwr2 = fs2.is_some()
                        && cs2.is_some()
                        && fs2.unwrap() <= cs2.unwrap()
                        && ctl.clono_print_opt.amino.contains(&"fwr2".to_string());
                    let show_fwr3 = fs3.is_some()
                        && fs3.unwrap() <= rsi.cdr3_starts[cx]
                        && ctl.clono_print_opt.amino.contains(&"fwr3".to_string());
                    let show_fwr4 = ctl.clono_print_opt.amino.contains(&"fwr4".to_string());
                    for (j, p) in show_aa[cx].iter().enumerate() {
                        if show_cdr1 && *p >= cs1.unwrap() / 3 && *p < fs2.unwrap() / 3 {
                            ft[j] = 1;
                        } else if show_cdr2 && *p >= cs2.unwrap() / 3 && *p < fs3.unwrap() / 3 {
                            ft[j] = 2;
                        } else if show_cdr3 && *p >= cs3 / 3 && *p < cs3 / 3 + n3 {
                            ft[j] = 3;
                        } else if show_fwr1 && *p >= fs1 / 3 && *p < cs1.unwrap() / 3 {
                            ft[j] = 4;
                        } else if show_fwr2 && *p >= fs2.unwrap() / 3 && *p < cs2.unwrap() / 3 {
                            ft[j] = 5;
                        } else if show_fwr3
                            && *p >= fs3.unwrap() / 3
                            && *p < rsi.cdr3_starts[cx] / 3
                        {
                            ft[j] = 6;
                        } else if show_fwr4 && *p >= cs3 / 3 + n3 {
                            ft[j] = 7;
                        }
                    }
                    field_types[cx] = ft;
                }

                // Build varmat.

                let mut varmat = vec![vec![vec![b'-']; cols]; nexacts];
                for col in 0..cols {
                    for u in 0..nexacts {
                        let m = mat[col][u];
                        if m.is_some() {
                            let mut v = Vec::<u8>::new();
                            let seq = rsi.seqss[col][u].clone();
                            for p in vars[col].iter() {
                                v.push(seq[*p]);
                            }
                            varmat[u][col] = v;
                        }
                    }
                }

                // Find the fields associated to nd<k> if used.

                let mut lvarsc = lvars.clone();
                let mut nd_fields = Vec::<String>::new();
                for (i, x) in lvars.iter().enumerate() {
                    if x.starts_with("nd")
                        && x.after("nd").parse::<usize>().is_ok()
                        && x.after("nd").force_usize() >= 1
                    {
                        lvarsc.clear();
                        for m in 0..i {
                            lvarsc.push(lvars[m].clone());
                        }
                        let k = x.after("nd").force_usize();
                        let mut n = vec![0 as usize; ctl.origin_info.n()];
                        for u in 0..nexacts {
                            let ex = &exact_clonotypes[exacts[u]];
                            for l in 0..ex.ncells() {
                                n[ex.clones[l][0].dataset_index] += 1;
                            }
                        }
                        let mut datasets = ctl.origin_info.dataset_id.clone();
                        // does not work for unknown reason, so "manually" replaced
                        // sort_sync2(&mut n, &mut datasets);
                        let permutation = permutation::sort(&n[..]);
                        n = permutation.apply_slice(&n[..]);
                        datasets = permutation.apply_slice(&datasets[..]);
                        n.reverse();
                        datasets.reverse();
                        for l in 0..n.len() {
                            if n[l] == 0 {
                                n.truncate(l);
                                datasets.truncate(l);
                                break;
                            }
                        }
                        for l in 0..k {
                            if l >= n.len() {
                                break;
                            }
                            nd_fields.push(format!("n_{}", datasets[l].clone()));
                            lvarsc.push(format!("n_{}", datasets[l].clone()));
                        }
                        if n.len() > k {
                            nd_fields.push("n_other".to_string());
                            lvarsc.push("n_other".to_string());
                        }
                        for m in i + 1..lvars.len() {
                            lvarsc.push(lvars[m].clone());
                        }
                    }
                }
                let lvars = lvarsc.clone();

                // Now build table content.

                let mut sr = Vec::<(Vec<String>, Vec<Vec<String>>, Vec<Vec<u8>>, usize)>::new();
                let mut groups = HashMap::<usize, Vec<usize>>::new();
                for i in 0..lvars.len() {
                    if lvars[i].starts_with('g') && lvars[i].after("g").parse::<usize>().is_ok() {
                        let d = lvars[i].after("g").force_usize();
                        if groups.contains_key(&d) {
                            continue;
                        }
                        let mut e: EquivRel = EquivRel::new(nexacts as i32);
                        for u1 in 0..nexacts {
                            let ex1 = &exact_clonotypes[exacts[u1]];
                            for u2 in u1 + 1..nexacts {
                                if e.class_id(u1 as i32) == e.class_id(u2 as i32) {
                                    continue;
                                }
                                let ex2 = &exact_clonotypes[exacts[u2]];
                                let mut diffs = 0;
                                'comp: for cx in 0..cols {
                                    let m1 = mat[cx][u1];
                                    if m1.is_none() {
                                        continue;
                                    }
                                    let m2 = mat[cx][u2];
                                    if m2.is_none() {
                                        continue;
                                    }
                                    let (m1, m2) = (m1.unwrap(), m2.unwrap());
                                    for p in vars[cx].iter() {
                                        if ex1.share[m1].seq_del[*p] != ex2.share[m2].seq_del[*p] {
                                            diffs += 1;
                                            if diffs > d {
                                                break 'comp;
                                            }
                                        }
                                    }
                                }
                                if diffs <= d {
                                    e.join(u1 as i32, u2 as i32);
                                }
                            }
                        }
                        let mut c = Vec::<usize>::new();
                        let mut reps = Vec::<i32>::new();
                        e.orbit_reps(&mut reps);
                        for u in 0..nexacts {
                            c.push(bin_position(&reps, &e.class_id(u as i32)) as usize);
                        }
                        groups.insert(d, c);
                    }
                }

                // Set up to record stats that assign a value to each cell for a given variable.

                let mut stats = Vec::<(String, Vec<f64>)>::new();

                // Compute some stats;

                let mut cred = Vec::<Vec<String>>::new();
                let mut pe = Vec::<Vec<String>>::new();
                let mut ppe = Vec::<Vec<String>>::new();
                let mut npe = Vec::<Vec<String>>::new();
                compute_some_stats(
                    &ctl,
                    &lvars,
                    &exacts,
                    &exact_clonotypes,
                    &gex_info,
                    &vdj_cells,
                    &n_vdj_gex,
                    &mut cred,
                    &mut pe,
                    &mut ppe,
                    &mut npe,
                );

                // Precompute for near and far.

                let mut fp = vec![Vec::<usize>::new(); varmat.len()]; // footprints
                for i in 0..varmat.len() {
                    for j in 0..varmat[i].len() {
                        if varmat[i][j] != vec![b'-'] {
                            fp[i].push(j);
                        }
                    }
                }

                // Build rows.

                let mut cell_count = 0;
                for u in 0..nexacts {
                    let mut typex = vec![false; cols];
                    let mut row = Vec::<String>::new();
                    let mut gex_low = 0;
                    let mut cx = Vec::<Vec<String>>::new();
                    for col in 0..cols {
                        cx.push(vec![String::new(); rsi.cvars[col].len()]);
                    }
                    let clonotype_id = exacts[u];
                    let ex = &exact_clonotypes[clonotype_id];
                    let mut d_all = vec![Vec::<u32>::new(); ex.clones.len()];
                    let mut ind_all = vec![Vec::<u32>::new(); ex.clones.len()];
                    row_fill(
                        pass,
                        u,
                        &ctl,
                        &exacts,
                        &mults,
                        &exact_clonotypes,
                        &gex_info,
                        &refdata,
                        &varmat,
                        &fp,
                        &vars_amino,
                        &show_aa,
                        &field_types,
                        &mut bads,
                        &mut gex_low,
                        &mut row,
                        &mut out_data,
                        &mut cx,
                        &mut d_all,
                        &mut ind_all,
                        &rsi,
                        &dref,
                        &groups,
                        &d_readers,
                        &ind_readers,
                        &h5_data,
                        &mut stats,
                        &vdj_cells,
                        &n_vdj_gex,
                        &lvars,
                        &nd_fields,
                        &peer_groups,
                        &extra_args,
                        &all_vars,
                        &fate,
                    );
                    let mut bli = Vec::<(String, usize, usize)>::new();
                    for l in 0..ex.clones.len() {
                        bli.push((
                            ex.clones[l][0].barcode.clone(),
                            ex.clones[l][0].dataset_index,
                            l,
                        ));
                    }
                    // WHY ARE WE SORTING HERE?
                    bli.sort();
                    for col in 0..cols {
                        if mat[col][u].is_some() {
                            typex[col] = true;
                        }
                    }
                    for r in 0..cx.len() {
                        for s in 0..cx[r].len() {
                            row.push(cx[r][s].clone());
                        }
                    }
                    res.5 = gex_low;

                    // Compute per-cell entries.

                    let mut subrows = Vec::<Vec<String>>::new();
                    compute_bu(
                        u,
                        cell_count,
                        &exacts,
                        &lvars,
                        &ctl,
                        &bli,
                        &ex,
                        &exact_clonotypes,
                        &mut row,
                        &mut subrows,
                        &varmat,
                        have_gex,
                        &gex_info,
                        &rsi,
                        &mut sr,
                        &fate,
                        &nd_fields,
                        &alt_bcs,
                        &cred,
                        &pe,
                        &ppe,
                        &npe,
                        &d_all,
                        &ind_all,
                        &mat,
                    );
                    cell_count += ex.clones.len();
                }
                let mut rord = Vec::<usize>::new(); // note that this is now superfluous
                for j in 0..sr.len() {
                    rord.push(j);
                }

                // Combine stats for the same variable.  This is needed because each exact
                // subclonotype contributes.  Note that we don't care about the order of the
                // values here (other than stability) because what we're going to do with them is
                // compute the mean or max.

                stats.sort_by(|a, b| a.0.cmp(&b.0));
                let stats_orig = stats.clone();
                let mut stats2 = Vec::<(String, Vec<f64>)>::new();
                let mut i = 0;
                while i < stats.len() {
                    let mut j = i + 1;
                    while j < stats.len() {
                        if stats[j].0 != stats[i].0 {
                            break;
                        }
                        j += 1;
                    }
                    let mut all = Vec::<f64>::new();
                    for k in i..j {
                        all.append(&mut stats[k].1.clone());
                    }
                    stats2.push((stats[i].0.clone(), all));
                    i = j;
                }
                stats = stats2;

                // Traverse the bounds and apply them.
                // Notes:
                // 1. This seems to run during both pass 1 and 2, and should only run
                //    during pass 1.
                // 2. The results of this can be counterintuitive, because the filtering is
                //    applied during pass 1, when there could be cells in the clonotype, that
                //    are removed by other filters.

                for bi in 0..ctl.clono_filt_opt.bounds.len() {
                    let x = &ctl.clono_filt_opt.bounds[bi];
                    let mut means = Vec::<f64>::new();
                    let mut maxs = Vec::<f64>::new();
                    // traverse the coefficients on the left hand side (each having a variable)
                    let mut fail = false;
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new(); // the stats for the variable
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                break;
                            }
                        }
                        let mut mean = 0.0;
                        let mut max = -1000_000_000.0_f64;
                        let mut count = 0;
                        for j in 0..vals.len() {
                            if !vals[j].is_nan() {
                                mean += vals[j];
                                max = max.max(vals[j]);
                                count += 1;
                            }
                        }
                        if count == 0 {
                            fail = true;
                        } else {
                            mean /= count as f64;
                            means.push(mean);
                            maxs.push(max);
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "mean" && (fail || !x.satisfied(&means))
                    {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for u in 0..nexacts {
                                bads[u] = true;
                            }
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "max" && (fail || !x.satisfied(&maxs)) {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for u in 0..nexacts {
                                bads[u] = true;
                            }
                        }
                    }
                }

                // Process COMPLETE.

                if ctl.gen_opt.complete {
                    let mut used = vec![false; cols];
                    for u in 0..nexacts {
                        if !bads[u] {
                            for m in 0..cols {
                                if mat[m][u].is_some() {
                                    used[m] = true;
                                }
                            }
                        }
                    }
                    for u in 0..nexacts {
                        for m in 0..cols {
                            if used[m] && mat[m][u].is_none() {
                                bads[u] = true;
                            }
                        }
                    }
                }

                // Done unless on second pass.

                if pass == 1 {
                    continue;
                }

                // See if we're in the test and control sets for gene scan.

                gene_scan_test(
                    &ctl,
                    &stats,
                    &stats_orig,
                    nexacts,
                    n,
                    &mut res.9,
                    &mut res.10,
                );

                // Make the table.

                let mut logz = String::new();
                finish_table(
                    n,
                    &ctl,
                    &exacts,
                    &exact_clonotypes,
                    &rsi,
                    &vars,
                    &show_aa,
                    &field_types,
                    &lvars,
                    &refdata,
                    &dref,
                    &peer_groups,
                    &mut mlog,
                    &mut logz,
                    &stats,
                    &mut sr,
                    &extra_args,
                    &pcols_sort,
                    &mut out_data,
                    &rord,
                    pass,
                );

                // Save.

                res.1.push(logz);
                res.2.push((exacts.clone(), rsi.clone()));
                res.12.push(in_center);
                for u in 0..exacts.len() {
                    res.8 += exact_clonotypes[exacts[u]].ncells() as isize;
                }
            }
            if pass == 2 {
                res.7.push(out_data);
            }
        }
    });

    for i in 0..results.len() {
        for j in 0..results[i].11.len() {
            fate[results[i].11[j].0].insert(results[i].11[j].1.clone(), results[i].11[j].2.clone());
        }
    }

    // Sort results in descending order by number of cells.

    results.sort_by_key(|x| -x.8);

    // Write loupe output.

    let mut all_loupe_clonotypes = Vec::<Clonotype>::new();
    for i in 0..results.len() {
        all_loupe_clonotypes.append(&mut results[i].6);
    }
    loupe_out(&ctl, all_loupe_clonotypes, &refdata, &dref);

    // Set up to group and print clonotypes.

    for i in 0..orbits.len() {
        for j in 0..results[i].1.len() {
            pics.push(results[i].1[j].clone());
            exacts.push(results[i].2[j].0.clone());
            rsi.push(results[i].2[j].1.clone());
            in_center.push(results[i].12[j]);
        }
        out_datas.append(&mut results[i].7);
    }

    // Gather some data for gene scan.

    if ctl.gen_opt.gene_scan_test.is_some() && !ctl.gen_opt.gene_scan_exact {
        for i in 0..orbits.len() {
            for j in 0..results[i].1.len() {
                if results[i].9[j] {
                    tests.push(i);
                }
                if results[i].10[j] {
                    controls.push(i);
                }
            }
        }
    }
    if ctl.gen_opt.gene_scan_test.is_some() && ctl.gen_opt.gene_scan_exact {
        for i in 0..exacts.len() {
            for j in 0..exacts[i].len() {
                if results[i].9[j] {
                    tests.push(exacts[i][j]);
                }
                if results[i].10[j] {
                    controls.push(exacts[i][j]);
                }
            }
        }
    }
}
