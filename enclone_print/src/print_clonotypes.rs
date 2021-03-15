// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::define_mat::*;
use crate::filter::*;
use crate::loupe::*;
use crate::print_utils1::*;
use crate::print_utils2::*;
use crate::print_utils3::*;
use crate::print_utils4::*;
use crate::print_utils5::*;
use amino::*;
use ansi_escape::*;
use enclone_core::defs::*;
use enclone_core::mammalian_fixed_len::*;
use enclone_core::print_tools::*;
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
    let mut tree_args = ctl.gen_opt.tree.clone();
    unique_sort(&mut tree_args);

    // Compute total cells.

    let mut total_cells = 0;
    for i in 0..exact_clonotypes.len() {
        total_cells += exact_clonotypes[i].ncells();
    }

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
                        &refdata,
                        &rsi,
                        &mut out_data,
                        &mut mlog,
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
                        &extra_parseables,
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

                // Combine stats for the same variable.  This is needed because each dataset
                // contributes.  Note that we don't care about the order of the values here
                // (other than stability) because what we're going to do with them is compute the
                // mean or max.

                stats.sort_by(|a, b| a.0.cmp(&b.0));
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

                for bi in 0..ctl.clono_filt_opt.bounds.len() {
                    let x = &ctl.clono_filt_opt.bounds[bi];
                    let mut means = Vec::<f64>::new();
                    let mut maxs = Vec::<f64>::new();
                    // traverse the coefficients on the left hand side (each having a variable)
                    for i in 0..x.n() {
                        let mut found = false;
                        let mut vals = Vec::<f64>::new(); // the stats for the variable
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            for u in 0..nexacts {
                                let clonotype_id = exacts[u];
                                let ex = &exact_clonotypes[clonotype_id];
                                for _ in 0..ex.clones.len() {
                                    vals.push(0.0);
                                }
                            }

                            /*
                            eprintln!(
                                "\nFailed to find the variable {} used in a \
                                 bound.  Please see \"enclone help filter\".\n",
                                x.var[i]
                            );
                            std::process::exit(1);
                            */
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
                        mean /= count as f64;
                        means.push(mean);
                        maxs.push(max);
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "mean" && !x.satisfied(&means) {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for u in 0..nexacts {
                                bads[u] = true;
                            }
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "max" && !x.satisfied(&maxs) {
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

                // See if we're in the test and control sets for gene scan.
                // uses: ctl, stats

                if ctl.gen_opt.gene_scan_test.is_some() {
                    let x = ctl.gen_opt.gene_scan_test.clone().unwrap();
                    let mut means = Vec::<f64>::new();
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new();
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                // found = true;
                                break;
                            }
                        }
                        let mut mean = 0.0;
                        for j in 0..vals.len() {
                            mean += vals[j];
                        }
                        mean /= n as f64;
                        means.push(mean);
                    }
                    res.9.push(x.satisfied(&means));
                    let x = ctl.gen_opt.gene_scan_control.clone().unwrap();
                    let mut means = Vec::<f64>::new();
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new();
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                break;
                            }
                        }
                        let mut mean = 0.0;
                        for j in 0..vals.len() {
                            mean += vals[j];
                        }
                        mean /= n as f64;
                        means.push(mean);
                    }
                    res.10.push(x.satisfied(&means));
                }

                // Done unless on second pass.

                if pass == 1 {
                    continue;
                }

                // Fill in exact_subclonotype_id, reorder.

                if ctl.parseable_opt.pout.len() > 0 || !ctl.gen_opt.tree.is_empty() {
                    for u in 0..nexacts {
                        macro_rules! speak {
                            ($u:expr, $var:expr, $val:expr) => {
                                if pass == 2
                                    && (ctl.parseable_opt.pout.len() > 0
                                        || !ctl.gen_opt.tree.is_empty())
                                {
                                    if pcols_sort.is_empty()
                                        || bin_member(&pcols_sort, &$var.to_string())
                                        || bin_member(&tree_args, &$var.to_string())
                                    {
                                        out_data[$u].insert($var.to_string(), $val);
                                    }
                                }
                            };
                        }
                        speak![rord[u], "exact_subclonotype_id", format!("{}", u + 1)];
                    }
                    let mut out_data2 = vec![HashMap::<String, String>::new(); nexacts];
                    for v in 0..nexacts {
                        out_data2[v] = out_data[rord[v]].clone();
                    }
                    out_data = out_data2;
                }

                // Add header text to mlog.

                add_header_text(&ctl, &exacts, &exact_clonotypes, &rord, &mat, &mut mlog);

                // Build table stuff.

                let mut row1 = Vec::<String>::new();
                let mut justify = Vec::<u8>::new();
                let mut rows = Vec::<Vec<String>>::new();
                let mut drows = Vec::<Vec<String>>::new();
                build_table_stuff(
                    &ctl,
                    &exacts,
                    &exact_clonotypes,
                    &rsi,
                    &vars,
                    &show_aa,
                    &field_types,
                    &mut row1,
                    &mut justify,
                    &mut drows,
                    &mut rows,
                    &lvars,
                );

                // Insert universal and donor reference rows.

                insert_reference_rows(
                    &ctl,
                    &rsi,
                    &show_aa,
                    &field_types,
                    &refdata,
                    &dref,
                    &row1,
                    &mut drows,
                    &mut rows,
                    &exacts,
                    &exact_clonotypes,
                    &peer_groups,
                );

                // Insert consensus row.

                if ctl.clono_print_opt.conx || ctl.clono_print_opt.conp {
                    let style;
                    if ctl.clono_print_opt.conx {
                        style = "x";
                    } else {
                        style = "p";
                    }
                    let mut row = Vec::<String>::new();
                    row.push("consensus".to_string());
                    for _ in 1..row1.len() {
                        row.push("\\ext".to_string());
                    }
                    let classes = aa_classes();
                    for col in 0..rsi.mat.len() {
                        for m in 0..rsi.cvars[col].len() {
                            if rsi.cvars[col][m] == "amino".to_string() {
                                let mut xdots = String::new();
                                for k in 0..show_aa[col].len() {
                                    if k > 0 && field_types[col][k] != field_types[col][k - 1] {
                                        xdots.push(' ');
                                    }
                                    let p = show_aa[col][k];
                                    let mut codons = Vec::<Vec<u8>>::new();
                                    for u in 0..nexacts {
                                        if mat[col][u].is_some() {
                                            let seq_amino = rsi.seqss_amino[col][u].clone();
                                            if 3 * p + 3 <= seq_amino.len() {
                                                codons.push(seq_amino[3 * p..3 * p + 3].to_vec());
                                            }
                                        }
                                    }
                                    unique_sort(&mut codons);
                                    let mut gap = false;
                                    for x in codons.iter() {
                                        if x.contains(&b'-') {
                                            gap = true;
                                        }
                                    }
                                    if codons.solo() && gap {
                                        xdots += &"g";
                                    } else if codons.solo() {
                                        let codon = &codons[0];
                                        let aa = codon_to_aa(&codon);
                                        let mut log = Vec::<u8>::new();
                                        emit_codon_color_escape(&codon, &mut log);
                                        log.push(aa);
                                        emit_end_escape(&mut log);
                                        xdots += &strme(&log);
                                    } else if gap {
                                        xdots += &"X";
                                    } else {
                                        let mut aas = Vec::<u8>::new();
                                        for x in codons.iter() {
                                            aas.push(codon_to_aa(x));
                                        }
                                        unique_sort(&mut aas);
                                        if aas.solo() {
                                            xdots.push(aas[0] as char);
                                        } else if style == "x" {
                                            xdots += &"X";
                                        } else {
                                            for m in classes.iter() {
                                                if meet_size(&aas, &m.1) == aas.len() {
                                                    xdots.push(m.0 as char);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                row.push(xdots);
                            } else {
                                row.push("".to_string());
                            }
                        }
                    }
                    rows.push(row);
                }

                // Insert horizontal line.

                if !drows.is_empty() {
                    let mut width = 1 + lvars.len();
                    for col in 0..cols {
                        width += rsi.cvars[col].len();
                    }
                    rows.push(vec!["\\hline".to_string(); width]);
                }

                // Build the diff row.

                build_diff_row(
                    &ctl,
                    &rsi,
                    &mut rows,
                    &mut drows,
                    &row1,
                    nexacts,
                    &field_types,
                    &show_aa,
                );

                // Finish building table content.

                for j in 0..sr.len() {
                    sr[j].0[0] = format!("{}", j + 1); // row number (#)
                    rows.push(sr[j].0.clone());
                    rows.append(&mut sr[j].1.clone());
                }

                // Add sum and mean rows.

                if ctl.clono_print_opt.sum {
                    let mut row = Vec::<String>::new();
                    row.push("Σ".to_string());
                    for i in 0..lvars.len() {
                        let mut x = lvars[i].clone();
                        if x.contains(':') {
                            x = x.before(":").to_string();
                        }
                        let mut found = false;
                        let mut total = 0.0;
                        for j in 0..stats.len() {
                            if stats[j].0 == x {
                                found = true;
                                for k in 0..stats[j].1.len() {
                                    total += stats[j].1[k];
                                }
                            }
                        }
                        if !found {
                            row.push(String::new());
                        } else {
                            if !lvars[i].ends_with("_%") {
                                row.push(format!("{}", total.round() as usize));
                            } else {
                                row.push(format!("{:.2}", total));
                            }
                        }
                    }
                    // This is necessary but should not be:
                    for cx in 0..cols {
                        for _ in 0..rsi.cvars[cx].len() {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }
                if ctl.clono_print_opt.mean {
                    let mut row = Vec::<String>::new();
                    row.push("μ".to_string());
                    for i in 0..lvars.len() {
                        let mut x = lvars[i].clone();
                        if x.contains(':') {
                            x = x.before(":").to_string();
                        }
                        let mut found = false;
                        let mut total = 0.0;
                        for j in 0..stats.len() {
                            if stats[j].0 == x {
                                found = true;
                                for k in 0..stats[j].1.len() {
                                    total += stats[j].1[k];
                                }
                            }
                        }
                        let mean = total / n as f64;
                        if !found {
                            row.push(String::new());
                        } else {
                            if !lvars[i].ends_with("_%") {
                                row.push(format!("{:.1}", mean));
                            } else {
                                row.push(format!("{:.2}", mean));
                            }
                        }
                    }
                    // This is necessary but should not be:
                    for cx in 0..cols {
                        for _ in 0..rsi.cvars[cx].len() {
                            row.push(String::new());
                        }
                    }
                    rows.push(row);
                }

                // Make table.

                for i in 0..rows.len() {
                    for j in 0..rows[i].len() {
                        rows[i][j] = rows[i][j].replace("|TRX", "TRB");
                        rows[i][j] = rows[i][j].replace("|TRY", "TRA");
                    }
                }
                for cx in 0..cols {
                    justify.push(b'|');
                    for m in 0..rsi.cvars[cx].len() {
                        justify.push(justification(&rsi.cvars[cx][m]));
                    }
                }
                let mut logz = String::new();
                make_table(&ctl, &mut rows, &justify, &mlog, &mut logz);

                // Add phylogeny.

                if ctl.toy {
                    let mut vrefs = Vec::<Vec<u8>>::new();
                    let mut jrefs = Vec::<Vec<u8>>::new();
                    for cx in 0..cols {
                        let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
                        for u in 0..nexacts {
                            let m = rsi.mat[cx][u];
                            if m.is_some() {
                                let m = m.unwrap();
                                jref = exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
                            }
                            let vseq1 = refdata.refs[rsi.vids[cx]].to_ascii_vec();
                            if rsi.vpids[cx].is_some() {
                                vref = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
                            } else {
                                vref = vseq1.clone();
                            }
                        }
                        vrefs.push(vref);
                        jrefs.push(jref);
                    }
                    for u1 in 0..nexacts {
                        let ex1 = &exact_clonotypes[exacts[u1]];
                        for u2 in u1 + 1..nexacts {
                            let ex2 = &exact_clonotypes[exacts[u2]];
                            let (mut d1, mut d2) = (0, 0);
                            let mut d = 0;
                            for cx in 0..cols {
                                let (m1, m2) = (rsi.mat[cx][u1], rsi.mat[cx][u2]);
                                if m1.is_none() || m2.is_none() {
                                    continue;
                                }
                                let (m1, m2) = (m1.unwrap(), m2.unwrap());
                                let (s1, s2) = (&ex1.share[m1].seq_del, &ex2.share[m2].seq_del);
                                let n = s1.len();
                                let (vref, jref) = (&vrefs[cx], &jrefs[cx]);
                                for j in 0..vars[cx].len() {
                                    let p = vars[cx][j];
                                    if s1[p] != s2[p] {
                                        if p < vref.len() - ctl.heur.ref_v_trim {
                                            if s1[p] == vref[p] {
                                                d1 += 1;
                                            } else if s2[p] == vref[p] {
                                                d2 += 1;
                                            }
                                        } else if p >= n - (jref.len() - ctl.heur.ref_j_trim) {
                                            if s1[p] == jref[jref.len() - (n - p)] {
                                                d1 += 1;
                                            } else if s2[p] == jref[jref.len() - (n - p)] {
                                                d2 += 1;
                                            }
                                        } else {
                                            d += 1;
                                        }
                                    }
                                }
                            }
                            if (d1 == 0) ^ (d2 == 0) {
                                if d1 == 0 {
                                    logz += &format!("{} ==> {}", u1 + 1, u2 + 1);
                                } else {
                                    logz += &format!("{} ==> {}", u2 + 1, u1 + 1);
                                }
                                let s = format!(
                                    "; u1 = {}, u2 = {}, d1 = {}, d2 = {}, d = {}\n",
                                    u1 + 1,
                                    u2 + 1,
                                    d1,
                                    d2,
                                    d
                                );
                                logz += &s;
                            }
                        }
                    }
                }

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

    // Gather some data for gene scan.

    if ctl.gen_opt.gene_scan_test.is_some() {
        let mut count = 0;
        for i in 0..orbits.len() {
            for j in 0..results[i].1.len() {
                if results[i].9[j] {
                    tests.push(count);
                }
                if results[i].10[j] {
                    controls.push(count);
                }
                count += 1;
            }
        }
    }

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
}
