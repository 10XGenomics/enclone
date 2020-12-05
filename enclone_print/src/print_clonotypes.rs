// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::filter::*;
use crate::loupe::*;
use crate::print_utils1::*;
use crate::print_utils2::*;
use crate::print_utils3::*;
use crate::print_utils4::*;
use crate::print_utils5::*;
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
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    orbits: &Vec<Vec<i32>>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    d_readers: &Vec<Option<hdf5::Reader>>,
    ind_readers: &Vec<Option<hdf5::Reader>>,
    h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
    pics: &mut Vec<String>,
    exacts: &mut Vec<Vec<usize>>,
    rsi: &mut Vec<ColInfo>,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    tests: &mut Vec<usize>,
    controls: &mut Vec<usize>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    // Make an abbreviation.

    let lvars = &ctl.clono_print_opt.lvars;

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
    // added out_datas
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
        // cdr3s:  sorted list of (chain_type:cdr3)
        // js:     indices into od (BE VERY CARFUL ABOUT USING THIS)

        let mut exacts = Vec::<usize>::new();
        let mut mults = Vec::<usize>::new();
        let mut cdr3s = Vec::<Vec<String>>::new();
        let mut cdr3s_len = Vec::<Vec<(String, usize)>>::new();
        let mut js = Vec::<usize>::new();
        let mut j = 0;
        let loupe_clonotypes = &mut res.6;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            let mut mult = 0 as usize;
            let mut z = Vec::<String>::new();
            let mut z_len = Vec::<(String, usize)>::new();
            for l in j..k {
                let x: &CloneInfo = &info[od[l].2 as usize];
                let m = x.clonotype_index;
                mult = exact_clonotypes[m].clones.len();
                for m in 0..x.cdr3_aa.len() {
                    // Do something EXTREMELY ugly.  To force TRB columns to come before TRA
                    // columns, rename then TRX and TRY.  This is reversed at the very end.

                    let mut c = x.chain_types[m].clone();
                    if c.starts_with("TRB") {
                        c = c.replacen("TRB", "TRX", 1);
                    } else if c.starts_with("TRA") {
                        c = c.replacen("TRA", "TRY", 1);
                    }
                    z.push(format!("{}:{}", c, x.cdr3_aa[m]));
                    z_len.push((format!("{}:{}", c, x.cdr3_aa[m]), x.lens[m]));
                }
            }
            unique_sort(&mut z);
            unique_sort(&mut z_len);
            cdr3s.push(z);
            cdr3s_len.push(z_len);
            js.push(j);
            let mut x = Vec::<usize>::new();
            for l in j..k {
                x.push(l);
            }
            mults.push(mult);
            exacts.push(od[j].1);
            j = k;
        }

        // There are two passes.  On the first pass we only identify the exact subclonotypes that
        // are junk.  On the second pass we remove those and then print the orbit.

        let mut bads = vec![false; cdr3s.len()];
        for pass in 1..=2 {
            // Delete weak exact subclonotypes.

            if pass == 2 && !ctl.clono_filt_opt.protect_bads {
                erase_if(&mut cdr3s, &bads);
                erase_if(&mut cdr3s_len, &bads);
                erase_if(&mut js, &bads);
                erase_if(&mut mults, &bads);
                erase_if(&mut exacts, &bads);
            }

            // Sort exact subclonotypes.

            let mat = define_mat(&ctl, &exact_clonotypes, &cdr3s_len, &js, &od, &info);
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
            cdr3s = permutation.apply_slice(&cdr3s[..]);
            cdr3s_len = permutation.apply_slice(&cdr3s_len[..]);
            js = permutation.apply_slice(&js[..]);
            exacts.reverse();
            mults.reverse();
            cdr3s.reverse();
            cdr3s_len.reverse();
            js.reverse();

            // Define a matrix mat[col][ex] which is the column of the exact subclonotype
            // corresponding to the given column col of the clonotype, which may or may not be
            // defined.  Then define other information associated to each chain.  These are
            // reference sequence identifiers, CDR3 start positions, and the like.

            let nexacts = exacts.len();
            let mat = define_mat(&ctl, &exact_clonotypes, &cdr3s_len, &js, &od, &info);
            let cols = mat.len();
            let mut rsi = define_column_info(&ctl, &exacts, &exact_clonotypes, &mat, &refdata);
            rsi.mat = mat;
            let mat = &rsi.mat;

            // Let n be the total number of cells in this pass.

            let n: usize = mults.iter().sum();

            // Filter.

            if pass == 2
                && !survives_filter(&exacts, &rsi, &ctl, &exact_clonotypes, &refdata, &gex_info)
            {
                continue;
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
            if n >= ctl.clono_filt_opt.ncells_low {
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

                // Compute "cred" stats (credibility/# of neighboring cells that are also
                // B cells).

                let mut cred = vec![Vec::<String>::new(); lvars.len()];
                for k in 0..lvars.len() {
                    if lvars[k] == "cred".to_string() {
                        for u in 0..nexacts {
                            let clonotype_id = exacts[u];
                            let ex = &exact_clonotypes[clonotype_id];
                            for l in 0..ex.clones.len() {
                                let bc = &ex.clones[l][0].barcode;
                                let li = ex.clones[l][0].dataset_index;
                                if gex_info.pca[li].contains_key(&bc.clone()) {
                                    let mut creds = 0;
                                    let mut z = Vec::<(f64, String)>::new();
                                    let x = &gex_info.pca[li][&bc.clone()];
                                    for y in gex_info.pca[li].iter() {
                                        let mut dist2 = 0.0;
                                        for m in 0..x.len() {
                                            dist2 += (y.1[m] - x[m]) * (y.1[m] - x[m]);
                                        }
                                        z.push((dist2, y.0.clone()));
                                    }
                                    z.sort_by(|a, b| a.partial_cmp(b).unwrap());
                                    let top = n_vdj_gex[li];
                                    for i in 0..top {
                                        if bin_member(&vdj_cells[li], &z[i].1) {
                                            creds += 1;
                                        }
                                    }
                                    let pc = 100.0 * creds as f64 / top as f64;
                                    cred[k].push(format!("{:.1}", pc));
                                } else {
                                    cred[k].push("".to_string());
                                }
                            }
                        }
                    }
                }

                // Compute pe (PCA distance).

                let mut pe = vec![Vec::<String>::new(); lvars.len()];
                for k in 0..lvars.len() {
                    if lvars[k].starts_with("pe") {
                        let n = lvars[k].after("pe").force_usize();
                        let mut bcs = Vec::<String>::new();
                        let mut count = 0;
                        let mut to_index = Vec::<usize>::new();
                        for u in 0..nexacts {
                            let clonotype_id = exacts[u];
                            let ex = &exact_clonotypes[clonotype_id];
                            for l in 0..ex.clones.len() {
                                let bc = &ex.clones[l][0].barcode;
                                let li = ex.clones[l][0].dataset_index;
                                if gex_info.pca[li].contains_key(&bc.clone()) {
                                    bcs.push(bc.to_string());
                                    to_index.push(count);
                                }
                                count += 1;
                            }
                        }
                        let mut e: EquivRel = EquivRel::new(bcs.len() as i32);
                        let li = 0; // BEWARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        let mut mat = vec![Vec::<f64>::new(); bcs.len()];
                        for i in 0..bcs.len() {
                            mat[i] = gex_info.pca[li][&bcs[i].clone()].clone();
                        }
                        for i1 in 0..bcs.len() {
                            for i2 in i1 + 1..bcs.len() {
                                if e.class_id(i1 as i32) != e.class_id(i2 as i32) {
                                    let mut d = 0.0;
                                    for j in 0..mat[i1].len() {
                                        d += (mat[i1][j] - mat[i2][j]) * (mat[i1][j] - mat[i2][j]);
                                    }
                                    d = d.sqrt();
                                    if d <= n as f64 {
                                        e.join(i1 as i32, i2 as i32);
                                    }
                                }
                            }
                        }
                        pe[k] = vec![String::new(); count];
                        let mut ids = Vec::<i32>::new();
                        for i in 0..bcs.len() {
                            ids.push(e.class_id(i as i32));
                        }
                        unique_sort(&mut ids);
                        let mut reps = Vec::<i32>::new();
                        e.orbit_reps(&mut reps);
                        reps.sort();
                        for i in 0..bcs.len() {
                            pe[k][to_index[i]] =
                                format!("{}", bin_position(&ids, &e.class_id(i as i32)));
                        }
                    }
                }

                // Compute ppe (PCA distance).

                let mut ppe = vec![Vec::<String>::new(); lvars.len()];
                for k in 0..lvars.len() {
                    if lvars[k].starts_with("ppe") {
                        let n = lvars[k].after("ppe").force_usize();
                        let mut bcs = Vec::<String>::new();
                        let mut count = 0;
                        let mut to_index = Vec::<usize>::new();
                        for u in 0..nexacts {
                            let clonotype_id = exacts[u];
                            let ex = &exact_clonotypes[clonotype_id];
                            for l in 0..ex.clones.len() {
                                let bc = &ex.clones[l][0].barcode;
                                let li = ex.clones[l][0].dataset_index;
                                if gex_info.pca[li].contains_key(&bc.clone()) {
                                    bcs.push(bc.to_string());
                                    to_index.push(count);
                                }
                                count += 1;
                            }
                        }
                        let li = 0; // BEWARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        let mut mat = vec![Vec::<f64>::new(); bcs.len()];
                        for i in 0..bcs.len() {
                            mat[i] = gex_info.pca[li][&bcs[i].clone()].clone();
                        }
                        let mut matg = Vec::<Vec<f64>>::new();
                        for i in gex_info.pca[li].iter() {
                            matg.push(i.1.to_vec());
                        }
                        let mut x = vec![0; bcs.len()];
                        for i1 in 0..mat.len() {
                            for i2 in 0..matg.len() {
                                let m1 = &mat[i1];
                                let m2 = &matg[i2];
                                let mut d = 0.0;
                                for j in 0..m1.len() {
                                    d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                                }
                                d = d.sqrt();
                                if d <= n as f64 {
                                    x[i1] += 1;
                                }
                            }
                        }
                        let mut y = vec![0; bcs.len()];
                        for i1 in 0..mat.len() {
                            for i2 in 0..mat.len() {
                                let m1 = &mat[i1];
                                let m2 = &mat[i2];
                                let mut d = 0.0;
                                for j in 0..m1.len() {
                                    d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                                }
                                d = d.sqrt();
                                if d <= n as f64 {
                                    y[i1] += 1;
                                }
                            }
                        }
                        ppe[k] = vec![String::new(); count];
                        for i in 0..bcs.len() {
                            ppe[k][to_index[i]] =
                                format!("{:.1}", 100.0 * y[i] as f64 / x[i] as f64);
                        }
                    }
                }

                // Compute npe (PCA distance).

                let mut npe = vec![Vec::<String>::new(); lvars.len()];
                for k in 0..lvars.len() {
                    if lvars[k].starts_with("npe") {
                        let n = lvars[k].after("npe").force_usize();
                        let mut bcs = Vec::<String>::new();
                        let mut count = 0;
                        let mut to_index = Vec::<usize>::new();
                        for u in 0..nexacts {
                            let clonotype_id = exacts[u];
                            let ex = &exact_clonotypes[clonotype_id];
                            for l in 0..ex.clones.len() {
                                let bc = &ex.clones[l][0].barcode;
                                let li = ex.clones[l][0].dataset_index;
                                if gex_info.pca[li].contains_key(&bc.clone()) {
                                    bcs.push(bc.to_string());
                                    to_index.push(count);
                                }
                                count += 1;
                            }
                        }
                        let li = 0; // BEWARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        let mut mat = vec![Vec::<f64>::new(); bcs.len()];
                        for i in 0..bcs.len() {
                            mat[i] = gex_info.pca[li][&bcs[i].clone()].clone();
                        }
                        let mut y = vec![0; bcs.len()];
                        for i1 in 0..mat.len() {
                            for i2 in 0..mat.len() {
                                let m1 = &mat[i1];
                                let m2 = &mat[i2];
                                let mut d = 0.0;
                                for j in 0..m1.len() {
                                    d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                                }
                                d = d.sqrt();
                                if d <= n as f64 {
                                    y[i1] += 1;
                                }
                            }
                        }
                        npe[k] = vec![String::new(); count];
                        for i in 0..bcs.len() {
                            npe[k][to_index[i]] = format!("{}", y[i]);
                        }
                    }
                }

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

                    // Very bad computation because of embedded binary search.

                    let mut subrows = Vec::<Vec<String>>::new();
                    if ctl.clono_print_opt.bu {
                        for bcl in bli.iter() {
                            let mut row = Vec::<String>::new();
                            let bc = &bcl.0;
                            let li = bcl.1;
                            let di = ex.clones[bcl.2][0].dataset_index;
                            row.push(format!("$  {}", bc.clone()));
                            let ex = &exact_clonotypes[exacts[u]];
                            for k in 0..lvars.len() {
                                let nr = row.len();
                                let mut filled = false;
                                for l in 0..ctl.origin_info.n() {
                                    if lvars[k] == format!("n_{}", ctl.origin_info.dataset_id[l]) {
                                        let mut n = 0;
                                        if di == l {
                                            n = 1;
                                        }
                                        row.push(format!("{}", n));
                                        filled = true;
                                    }
                                }
                                if filled {
                                } else if lvars[k] == "n_b".to_string() {
                                    let mut n = 0;
                                    let li = ex.clones[bcl.2][0].dataset_index;
                                    if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                        if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                                            n = 1;
                                        }
                                    }
                                    row.push(format!("{}", n));
                                } else if lvars[k] == "filter".to_string() {
                                    let mut f = String::new();
                                    if fate[li].contains_key(&bc.clone()) {
                                        f = fate[li][&bc.clone()].clone();
                                        f = f.between(" ", " ").to_string();
                                    }
                                    row.push(f);
                                } else if lvars[k] == "n_other".to_string() {
                                    let mut n = 0;
                                    let di = ex.clones[bcl.2][0].dataset_index;
                                    let f = format!("n_{}", ctl.origin_info.dataset_id[di]);
                                    let mut found = false;
                                    for i in 0..nd_fields.len() {
                                        if f == nd_fields[i] {
                                            found = true;
                                        }
                                    }
                                    if !found {
                                        n = 1;
                                    }
                                    row.push(format!("{}", n));
                                } else if lvars[k] == "sec".to_string() {
                                    let mut n = 0;
                                    if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                                        n = ctl.origin_info.secmem[li][&bc.clone()].0;
                                    }
                                    row.push(format!("{}", n));
                                } else if lvars[k] == "mem".to_string() {
                                    let mut n = 0;
                                    if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                                        n = ctl.origin_info.secmem[li][&bc.clone()].1;
                                    }
                                    row.push(format!("{}", n));
                                } else if bin_member(&alt_bcs, &lvars[k]) {
                                    let mut val = String::new();
                                    let alt = &ctl.origin_info.alt_bc_fields[li];
                                    for j in 0..alt.len() {
                                        if alt[j].0 == lvars[k] {
                                            if alt[j].1.contains_key(&bc.clone()) {
                                                val = alt[j].1[&bc.clone()].clone();
                                            }
                                        }
                                    }
                                    row.push(val);
                                } else if lvars[k] == "datasets".to_string() {
                                    row.push(format!("{}", ctl.origin_info.dataset_id[li].clone()));
                                } else if lvars[k] == "clust".to_string() && have_gex {
                                    let mut cid = 0;
                                    if gex_info.cluster[li].contains_key(&bc.clone()) {
                                        cid = gex_info.cluster[li][&bc.clone()];
                                    }
                                    row.push(format!("{}", cid));
                                } else if lvars[k].starts_with("pe") && have_gex {
                                    row.push(format!("{}", pe[k][cell_count + bcl.2]));
                                } else if lvars[k].starts_with("npe") && have_gex {
                                    row.push(format!("{}", npe[k][cell_count + bcl.2]));
                                } else if lvars[k].starts_with("ppe") && have_gex {
                                    row.push(format!("{}", ppe[k][cell_count + bcl.2]));
                                } else if lvars[k] == "cred".to_string() && have_gex {
                                    row.push(format!("{}", cred[k][cell_count + bcl.2]));
                                } else if lvars[k] == "type".to_string() && have_gex {
                                    let mut cell_type = "".to_string();
                                    if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                        cell_type = gex_info.cell_type[li][&bc.clone()].clone();
                                    }
                                    row.push(cell_type);
                                } else if lvars[k] == "n_gex".to_string() && have_gex {
                                    let mut n_gex = 0;
                                    if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                                        n_gex = 1;
                                    }
                                    row.push(format!("{}", n_gex));
                                } else if lvars[k] == "mark".to_string() {
                                    let mut mark = String::new();
                                    if ex.clones[bcl.2][0].marked {
                                        mark = "x".to_string();
                                    }
                                    row.push(mark);
                                } else if lvars[k] == "entropy".to_string() && have_gex {
                                    // NOTE DUPLICATION WITH CODE BELOW.
                                    let mut gex_count = 0;
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        let mut raw_count = 0;
                                        if gex_info.gex_matrices[li].initialized() {
                                            let row = gex_info.gex_matrices[li].row(p as usize);
                                            for j in 0..row.len() {
                                                let f = row[j].0;
                                                let n = row[j].1;
                                                if gex_info.is_gex[li][f] {
                                                    raw_count += n;
                                                }
                                            }
                                        } else {
                                            let l = bcl.2;
                                            for j in 0..d_all[l].len() {
                                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                                    raw_count += d_all[l][j] as usize;
                                                }
                                            }
                                        }
                                        gex_count = raw_count;
                                    }
                                    let mut entropy = 0.0;
                                    if p >= 0 {
                                        if gex_info.gex_matrices[li].initialized() {
                                            let row = gex_info.gex_matrices[li].row(p as usize);
                                            for j in 0..row.len() {
                                                let f = row[j].0;
                                                let n = row[j].1;
                                                if gex_info.is_gex[li][f] {
                                                    let q = n as f64 / gex_count as f64;
                                                    entropy -= q * q.log2();
                                                }
                                            }
                                        } else {
                                            let l = bcl.2;
                                            for j in 0..d_all[l].len() {
                                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                                    let n = d_all[l][j] as usize;
                                                    let q = n as f64 / gex_count as f64;
                                                    entropy -= q * q.log2();
                                                }
                                            }
                                        }
                                    }
                                    row.push(format!("{:.2}", entropy));
                                } else if have_gex {
                                    // this calc isn't needed except in _% case below
                                    // TODO: ELIMINATE UNNEEDED CALC
                                    let mut gex_count = 0.0;
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        let mut raw_count = 0 as f64;
                                        if gex_info.gex_matrices[li].initialized() {
                                            let row = gex_info.gex_matrices[li].row(p as usize);
                                            for j in 0..row.len() {
                                                let f = row[j].0;
                                                let n = row[j].1;
                                                if gex_info.is_gex[li][f] {
                                                    raw_count += n as f64;
                                                }
                                            }
                                        } else {
                                            let l = bcl.2;
                                            for j in 0..d_all[l].len() {
                                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                                    raw_count += d_all[l][j] as f64;
                                                }
                                            }
                                        }
                                        if !ctl.gen_opt.full_counts {
                                            gex_count = raw_count * gex_info.gex_mults[li];
                                        } else {
                                            gex_count = raw_count;
                                        }
                                    }
                                    if lvars[k] == "gex".to_string() {
                                        row.push(format!("{}", gex_count.round()));
                                    } else {
                                        let mut y = lvars[k].clone();
                                        if y.contains(':') {
                                            y = y.after(":").to_string();
                                        }
                                        let y0 = y.clone();
                                        let suffixes = ["_min", "_max", "_μ", "_Σ", "_cell", "_%"];
                                        for s in suffixes.iter() {
                                            if y.ends_with(s) {
                                                y = y.rev_before(&s).to_string();
                                                break;
                                            }
                                        }
                                        let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                        let mut computed = false;
                                        let mut count = 0.0;
                                        let l = bcl.2;
                                        if p >= 0 {
                                            let mut ux = Vec::<usize>::new();
                                            if ctl.clono_print_opt.regex_match[li].contains_key(&y)
                                            {
                                                ux =
                                                    ctl.clono_print_opt.regex_match[li][&y].clone();
                                            }
                                            if ux.len() > 0 {
                                                computed = true;
                                                for fid in ux.iter() {
                                                    let counti = get_gex_matrix_entry(
                                                        &ctl, &gex_info, *fid, &d_all, &ind_all,
                                                        li, l, p as usize, &y,
                                                    );
                                                    count += counti;
                                                }
                                            } else if gex_info.feature_id[li].contains_key(&y) {
                                                computed = true;
                                                let fid = gex_info.feature_id[li][&y];
                                                count = get_gex_matrix_entry(
                                                    &ctl, &gex_info, fid, &d_all, &ind_all, li, l,
                                                    p as usize, &y,
                                                );
                                            }
                                        }
                                        if computed {
                                            // note unneeded calculation above in certain cases
                                            // TODO: ELIMINATE!
                                            if y0.ends_with("_min") {
                                            } else if y0.ends_with("_max") {
                                            } else if y0.ends_with("_μ") {
                                            } else if y0.ends_with("_Σ") {
                                            } else if y0.ends_with("_%") {
                                                row.push(format!(
                                                    "{:.2}",
                                                    (100.0 * count) / gex_count
                                                ));
                                            } else {
                                                row.push(format!("{}", count.round()));
                                            }
                                        }
                                    }
                                }
                                if row.len() == nr {
                                    row.push("".to_string());
                                }
                            }
                            let mut ncall = 0;
                            for k in 0..cols {
                                ncall += rsi.cvars[k].len();
                            }
                            let mut cx = vec!["".to_string(); ncall];
                            let mut cp = 0;
                            for col in 0..cols {
                                let m = mat[col][u];
                                if m.is_some() {
                                    let m = m.unwrap();
                                    for p in 0..rsi.cvars[col].len() {
                                        if rsi.cvars[col][p] == "u".to_string() {
                                            let numi = ex.clones[bcl.2][m].umi_count;
                                            cx[cp + p] = format!("{}", numi);
                                        } else if rsi.cvars[col][p] == "r".to_string() {
                                            let r = ex.clones[bcl.2][m].read_count;
                                            cx[cp + p] = format!("{}", r);
                                        }
                                    }
                                }
                                cp += rsi.cvars[col].len();
                            }
                            row.append(&mut cx);
                            subrows.push(row);
                        }
                    }
                    sr.push((row, subrows, varmat[u].clone(), u));
                    cell_count += ex.clones.len();
                }
                let mut rord = Vec::<usize>::new(); // note that this is now superfluous
                for j in 0..sr.len() {
                    rord.push(j);
                }

                // Apply bounds.  Before sorting we check for non-numbers because otherwise you'll
                // get an inscrutable traceback.

                for i in 0..stats.len() {
                    for j in 0..stats[i].1.len() {
                        if !stats[i].1[j].is_finite() {
                            panic!(
                                "About to sort but there's a non-finite value, which would \
                                cause the sort to fail.  This is a bug."
                            );
                        }
                    }
                }
                stats.sort_by(|a, b| a.partial_cmp(b).unwrap());
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
                for i in 0..ctl.clono_filt_opt.bounds.len() {
                    let x = &ctl.clono_filt_opt.bounds[i];
                    let mut means = Vec::<f64>::new();
                    let mut maxs = Vec::<f64>::new();
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new();
                        let mut found = false;
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            eprintln!(
                                "\nFailed to find the variable {} used in a \
                                 bound.  Please see \"enclone help filter\".\n",
                                x.var[i]
                            );
                            std::process::exit(1);
                        }
                        let mut mean = 0.0;
                        let mut max = -1000_000_000.0_f64;
                        for j in 0..vals.len() {
                            mean += vals[j];
                            max = max.max(vals[j]);
                        }
                        mean /= n as f64;
                        means.push(mean);
                        maxs.push(max);
                    }
                    if ctl.clono_filt_opt.bound_type[i] == "mean" && !x.satisfied(&means) {
                        for u in 0..nexacts {
                            bads[u] = true;
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[i] == "max" && !x.satisfied(&maxs) {
                        for u in 0..nexacts {
                            bads[u] = true;
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

                if ctl.parseable_opt.pout.len() > 0 {
                    for u in 0..nexacts {
                        macro_rules! speak {
                            ($u:expr, $var:expr, $val:expr) => {
                                if pass == 2 && ctl.parseable_opt.pout.len() > 0 {
                                    if pcols_sort.is_empty()
                                        || bin_member(&pcols_sort, &$var.to_string())
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

                // Insert horizontal line.

                if !drows.is_empty() {
                    let mut width = 1 + lvars.len();
                    for col in 0..cols {
                        width += rsi.cvars[col].len();
                    }
                    rows.push(vec!["\\hline".to_string(); width]);
                }

                // Build the diff row.

                let diff_pos = rows.len();
                if !drows.is_empty() {
                    let mut row = row1.clone();
                    for col in 0..cols {
                        for m in 0..rsi.cvars[col].len() {
                            if rsi.cvars[col][m] == "amino".to_string() {
                                let mut xdots = String::new();
                                for k in 0..show_aa[col].len() {
                                    if k > 0 && field_types[col][k] != field_types[col][k - 1] {
                                        xdots.push(' ');
                                    }
                                    let p = show_aa[col][k];
                                    let q = 3 * p;
                                    let leader = q < rsi.fr1_starts[col];
                                    let mut cdr = false;
                                    if rsi.cdr1_starts[col].is_some()
                                        && q >= rsi.cdr1_starts[col].unwrap()
                                        && q < rsi.fr2_starts[col].unwrap()
                                    {
                                        cdr = true;
                                    }
                                    if q >= rsi.cdr2_starts[col].unwrap()
                                        && q < rsi.fr3_starts[col].unwrap()
                                    {
                                        cdr = true;
                                    }
                                    if q >= rsi.cdr3_starts[col]
                                        && q < rsi.cdr3_starts[col] + 3 * rsi.cdr3_lens[col]
                                    {
                                        cdr = true;
                                    }
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
                                    if codons.len() > 1 {
                                        if cdr {
                                            if ctl.gen_opt.diff_style == "C1".to_string() {
                                                xdots.push('C');
                                            } else if ctl.gen_opt.diff_style == "C2".to_string() {
                                                xdots.push('');
                                                xdots.push('[');
                                                xdots.push('0');
                                                xdots.push('1');
                                                xdots.push('m');
                                                xdots.push('');
                                                xdots.push('[');
                                                xdots.push('3');
                                                xdots.push('1');
                                                xdots.push('m');
                                                xdots.push('◼');
                                                xdots.push('');
                                                xdots.push('[');
                                                xdots.push('0');
                                                xdots.push('1');
                                                xdots.push('m');
                                                xdots.push('');
                                                xdots.push('[');
                                                xdots.push('3');
                                                xdots.push('0');
                                                xdots.push('m');
                                            } else {
                                                xdots.push('x');
                                            }
                                        } else if !leader {
                                            if ctl.gen_opt.diff_style == "C1".to_string() {
                                                xdots.push('F');
                                            } else if ctl.gen_opt.diff_style == "C2".to_string() {
                                                xdots.push('▮');
                                            } else {
                                                xdots.push('x');
                                            }
                                        } else {
                                            if ctl.gen_opt.diff_style == "C1".to_string() {
                                                xdots.push('L');
                                            } else if ctl.gen_opt.diff_style == "C2".to_string() {
                                                xdots.push('▮');
                                            } else {
                                                xdots.push('x');
                                            }
                                        }
                                    } else {
                                        xdots.push('.');
                                    }
                                }
                                row.push(xdots);
                            } else {
                                row.push(rsi.cvars[col][m].clone());
                            }
                        }
                        for i in 0..row.len() {
                            row[i] = format!("[01m{}[0m", row[i]);
                        }
                    }
                    rows.push(row);
                } else {
                    for i in 0..row1.len() {
                        rows[diff_pos - 1][i] = row1[i].clone();
                    }
                }

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
        }
        out_datas.append(&mut results[i].7);
    }
}
