// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::defs::*;
use crate::filter::*;
use crate::group::*;
use crate::loupe::*;
use crate::plot::*;
use crate::print_utils1::*;
use crate::print_utils2::*;
use crate::print_utils3::*;
use crate::print_utils4::*;
use crate::print_utils5::*;
use crate::types::*;
use equiv::EquivRel;
use io_utils::*;
use ndarray::s;
use rayon::prelude::*;
use stats_utils::*;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
use string_utils::*;
use tables::*;
use vdj_ann::refx::*;
use vector_utils::*;

// TO DO: provide alternate coloring by amino acid properties:
// 1. Aliphatic: A, G, I, L, P, V
// 2. Aromatic: F, W, Y
// 3. Acidic: D, E
// 4. Basic: R, H, K
// 5. Hydroxylic: S, T
// 6. Sulfurous: C, M
// 7. Amidic: N, Q

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
    tall: &Instant,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    eq: &EquivRel,
    gex_info: &GexInfo,
    join_info: &Vec<(usize, usize, bool, Vec<u8>)>,
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

    // Compute orbits.

    let mut reps = Vec::<i32>::new();
    eq.orbit_reps(&mut reps);

    // Test for presence of gex data.

    let mut have_gex = false;
    for i in 0..ctl.sample_info.gex_path.len() {
        if ctl.sample_info.gex_path[i].len() > 0 {
            have_gex = true;
        }
    }

    // Load the GEX data.

    let mut d_readers = Vec::<Option<h5::Reader>>::new();
    let mut ind_readers = Vec::<Option<h5::Reader>>::new();
    if ctl.gen_opt.h5 {
        for li in 0..ctl.sample_info.n() {
            if ctl.sample_info.gex_path[li].len() > 0 {
                d_readers.push(Some(gex_info.h5_data[li].as_ref().unwrap().as_reader()));
                ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
            } else {
                d_readers.push(None);
                ind_readers.push(None);
            }
        }
    }
    let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
    if ctl.gen_opt.h5 && ctl.gen_opt.h5_pre {
        for li in 0..ctl.sample_info.n() {
            h5_data.push((li, Vec::new(), Vec::new()));
        }
        h5_data.par_iter_mut().for_each(|res| {
            let li = res.0;
            if ctl.sample_info.gex_path[li].len() > 0 {
                res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
                res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
            }
        });
    }

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
        Vec<(Vec<usize>, Vec<Vec<Option<usize>>>)>,
        usize,
        usize,
        usize,
        Vec<Clonotype>,
        Vec<Vec<HashMap<String, String>>>,
        isize,
        Vec<bool>,
        Vec<bool>,
    )>::new();
    for i in 0..reps.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<(Vec<usize>, Vec<Vec<Option<usize>>>)>::new(),
            0,
            0,
            0,
            Vec::<Clonotype>::new(),
            Vec::<Vec<HashMap<String, String>>>::new(),
            0,
            Vec::<bool>::new(),
            Vec::<bool>::new(),
        ));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        eq.orbit(reps[i], &mut o);
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
        let mut js = Vec::<usize>::new();
        let mut j = 0;
        let loupe_clonotypes = &mut res.6;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            let mut mult = 0 as usize;
            let mut z = Vec::<String>::new();
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
                }
            }
            unique_sort(&mut z);
            cdr3s.push(z);
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
                erase_if(&mut js, &bads);
                erase_if(&mut mults, &bads);
                erase_if(&mut exacts, &bads);
            }

            // Sort exact subclonotypes.

            let mat = define_mat(&ctl, &exact_clonotypes, &cdr3s, &js, &od, &info);
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
            js = permutation.apply_slice(&js[..]);
            exacts.reverse();
            mults.reverse();
            cdr3s.reverse();
            js.reverse();

            // Define a matrix mat[col][ex] which is the column of the exact subclonotype
            // corresponding to the given column col of the clonotype, which may or may not be
            // defined.  Then define other information associated to each chain.  These are
            // reference sequence identifiers, CDR3 start positions, and the like.

            let nexacts = exacts.len();
            let mat = define_mat(&ctl, &exact_clonotypes, &cdr3s, &js, &od, &info);
            let cols = mat.len();
            let mut rsi = define_column_info(&ctl, &exacts, &exact_clonotypes, &mat, &refdata);
            rsi.mat = mat;
            let mat = &rsi.mat;

            // Let n be the total number of cells in this pass.

            let n: usize = mults.iter().sum();

            // Filter.

            if pass == 2 && !survives_filter(&exacts, &rsi, &ctl, &exact_clonotypes, &refdata) {
                continue;
            }

            // Generate loupe data.

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
                        &vars,
                        &mut bads,
                    );
                }

                // Done unless on second pass.  Unless there are bounds.

                if pass == 1 && ctl.clono_filt_opt.bounds.len() == 0 {
                    continue;
                }

                // Define amino acid positions to show.

                let show_aa =
                    build_show_aa(&ctl, &rsi, &vars_amino, &shares_amino, &refdata, &dref);

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

                // Now build table content.
                // sr: now pretty pointless

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

                // Build rows.

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
                        &vars_amino,
                        &show_aa,
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
                    );
                    let mut bli = Vec::<(String, usize, usize)>::new();
                    for l in 0..ex.clones.len() {
                        bli.push((
                            ex.clones[l][0].barcode.clone(),
                            ex.clones[l][0].dataset_index,
                            l,
                        ));
                    }
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
                        for (kb, bcl) in bli.iter().enumerate() {
                            let mut row = Vec::<String>::new();
                            let bc = &bcl.0;
                            let li = bcl.1;
                            row.push(format!("$  {}", bc.clone()));
                            for k in 0..lvars.len() {
                                let nr = row.len();
                                if lvars[k] == "datasets".to_string() {
                                    row.push(format!("{}", ctl.sample_info.dataset_id[li].clone()));
                                } else if lvars[k] == "n_gex".to_string() && have_gex {
                                    let mut n_gex = 0;
                                    if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                                        n_gex = 1;
                                    }
                                    row.push(format!("{}", n_gex));
                                } else if lvars[k] == "entropy".to_string() && have_gex {
                                    // NOTE DUPLICATION WITH CODE BELOW.
                                    let mut gex_count = 0;
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        let mut raw_count = 0;
                                        if !ctl.gen_opt.h5 {
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
                                        if !ctl.gen_opt.h5 {
                                            let row = gex_info.gex_matrices[li].row(p as usize);
                                            for j in 0..row.len() {
                                                let f = row[j].0;
                                                let n = row[j].1;
                                                if gex_info.is_gex[li][f] {
                                                    let q = n as f64 / gex_count as f64;
                                                    entropy += q * q.log2();
                                                }
                                            }
                                        } else {
                                            let l = bcl.2;
                                            for j in 0..d_all[l].len() {
                                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                                    let n = d_all[l][j] as usize;
                                                    let q = n as f64 / gex_count as f64;
                                                    entropy += q * q.log2();
                                                }
                                            }
                                        }
                                    }
                                    row.push(format!("{:.2}", entropy));
                                } else if lvars[k] == "gex".to_string() && have_gex {
                                    let mut gex_count = 0;
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        let mut raw_count = 0 as f64;
                                        if !ctl.gen_opt.h5 {
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
                                            gex_count = (raw_count * gex_info.gex_mults[li]).round()
                                                as usize;
                                        } else {
                                            gex_count = raw_count.round() as usize;
                                        }
                                    }
                                    row.push(format!("{}", gex_count));
                                } else {
                                    let mut y = lvars[k].clone();
                                    if y.contains(':') {
                                        y = y.after(":").to_string();
                                    }
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        if gex_info.feature_id[li].contains_key(&y) {
                                            let fid = gex_info.feature_id[li][&y];
                                            let mut raw_count = 0 as f64;
                                            if !ctl.gen_opt.h5 {
                                                raw_count = gex_info.gex_matrices[li]
                                                    .value(p as usize, fid)
                                                    as f64;
                                            } else {
                                                for j in 0..d_all[kb].len() {
                                                    if ind_all[kb][j] == fid as u32 {
                                                        raw_count = d_all[kb][j] as f64;
                                                        break;
                                                    }
                                                }
                                            }
                                            let mult: f64;
                                            if y.ends_with("_g") {
                                                mult = gex_info.gex_mults[li];
                                            } else {
                                                mult = gex_info.fb_mults[li];
                                            }
                                            let count;
                                            if !ctl.gen_opt.full_counts {
                                                count = (raw_count * mult).round() as f64;
                                            } else {
                                                count = raw_count.round() as f64;
                                            }
                                            row.push(format!("{}", count));
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
                            let ex = &exact_clonotypes[exacts[u]];
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
                }
                let mut rord = Vec::<usize>::new(); // note that this is now superfluous
                for j in 0..sr.len() {
                    rord.push(j);
                }

                // Apply bounds.

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
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new();
                        // let mut found = false;
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                // found = true;
                                break;
                            }
                        }
                        /*
                        if !found {
                            eprintln!(
                                "\nFailed to find the variable {} used in a \
                                 bound.  Please see \"enclone help filter\".\n",
                                x.var[i]
                            );
                            std::process::exit(1);
                        }
                        */
                        let mut mean = 0.0;
                        for j in 0..vals.len() {
                            mean += vals[j];
                        }
                        mean /= n as f64;
                        means.push(mean);
                    }
                    if !x.satisfied(&means) {
                        for u in 0..nexacts {
                            bads[u] = true;
                        }
                    }
                }

                // See if we're in the test and control sets for gene scan.

                if ctl.gen_opt.gene_scan_test.is_some() {
                    let x = ctl.gen_opt.gene_scan_test.clone().unwrap();
                    let mut means = Vec::<f64>::new();
                    for i in 0..x.n() {
                        let mut vals = Vec::<f64>::new();
                        // let mut found = false;
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                // found = true;
                                break;
                            }
                        }
                        /*
                        if !found {
                            eprintln!(
                                "\nFailed to find the variable {} used in a \
                                 bound.  Please see \"enclone help filter\".\n",
                                x.var[i]
                            );
                            std::process::exit(1);
                        }
                        */
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
                        // let mut found = false;
                        for j in 0..stats.len() {
                            if stats[j].0 == x.var[i] {
                                vals.append(&mut stats[j].1.clone());
                                // found = true;
                                break;
                            }
                        }
                        /*
                        if !found {
                            eprintln!(
                                "\nFailed to find the variable {} used in a \
                                 bound.  Please see \"enclone help filter\".\n",
                                x.var[i]
                            );
                            std::process::exit(1);
                        }
                        */
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
                    &mut row1,
                    &mut justify,
                    &mut drows,
                    &mut rows,
                );

                // Insert universal and donor reference rows.

                insert_reference_rows(
                    &ctl, &rsi, &show_aa, &refdata, &dref, &row1, &mut drows, &mut rows,
                );

                // Insert horizontal line.

                let mut width = 1 + lvars.len();
                for col in 0..cols {
                    width += rsi.cvars[col].len();
                }
                rows.push(vec!["\\hline".to_string(); width]);

                // Insert placeholder for dots row.

                let cvars = &ctl.clono_print_opt.cvars;
                let diff_pos = rows.len();
                if !ctl.clono_print_opt.amino.is_empty() || cvars.contains(&"var".to_string()) {
                    let row = Vec::<String>::new();
                    rows.push(row);
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
                            row.push(format!("{}", total.round() as usize));
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
                            row.push(format!("{:.1}", mean));
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

                // Make the diff row.

                make_diff_row(&ctl, &rsi, cols, diff_pos, &drows, &mut row1, &mut rows);

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
                        if rsi.cvars[cx][m] == "amino".to_string()
                            || rsi.cvars[cx][m] == "var".to_string()
                            || rsi.cvars[cx][m] == "const".to_string()
                            || rsi.cvars[cx][m] == "cdr3_dna".to_string()
                            || rsi.cvars[cx][m] == "cdiff".to_string()
                            || rsi.cvars[cx][m] == "notes".to_string()
                            || rsi.cvars[cx][m] == "edit".to_string()
                        {
                            justify.push(b'l');
                        } else {
                            justify.push(b'r');
                        }
                    }
                }
                let mut logz = String::new();
                make_table(&ctl, &mut rows, &justify, &mlog, &mut logz);

                // Add phyologeny.

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
                res.2.push((exacts.clone(), mat.clone()));
                for u in 0..exacts.len() {
                    res.8 += exact_clonotypes[exacts[u]].ncells() as isize;
                }
            }
            if pass == 2 {
                res.7.push(out_data);
            }
        }
    });

    // Sort results in descending order by number of cells.

    results.sort_by_key(|x| -x.8);

    // Write loupe output.

    let mut all_loupe_clonotypes = Vec::<Clonotype>::new();
    for i in 0..results.len() {
        all_loupe_clonotypes.append(&mut results[i].6);
    }
    loupe_out(&ctl, all_loupe_clonotypes, &refdata, &dref);

    // Group and print clonotypes.  For now, limited functionality.

    let mut pics = Vec::<String>::new();
    let mut exacts = Vec::<Vec<usize>>::new(); // ugly reuse of name
    let mut mat = Vec::<Vec<Vec<Option<usize>>>>::new(); // ditto
    let mut out_datas = Vec::<Vec<HashMap<String, String>>>::new();
    for i in 0..reps.len() {
        for j in 0..results[i].1.len() {
            pics.push(results[i].1[j].clone());
            exacts.push(results[i].2[j].0.clone());
            mat.push(results[i].2[j].1.clone());
        }
        out_datas.append(&mut results[i].7);
    }
    group_and_print_clonotypes(
        &tall,
        &refdata,
        &pics,
        &exacts,
        &mat,
        &exact_clonotypes,
        &ctl,
        &parseable_fields,
        &mut out_datas,
        &join_info,
    );

    // Do gene scan.

    if ctl.gen_opt.gene_scan_test.is_some() {
        println!("\nFEATURE SCAN\n");
        let mut tests = Vec::<usize>::new();
        let mut controls = Vec::<usize>::new();
        let mut count = 0;
        for i in 0..reps.len() {
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
        let mut test_cells = 0;
        for i in tests.iter() {
            for u in exacts[*i].iter() {
                test_cells += exact_clonotypes[*u].ncells();
            }
        }
        println!(
            "{} clonotypes containing {} cells in test set",
            tests.len(),
            test_cells
        );
        let mut control_cells = 0;
        for i in controls.iter() {
            for u in exacts[*i].iter() {
                control_cells += exact_clonotypes[*u].ncells();
            }
        }
        println!(
            "{} clonotypes containing {} cells in control set\n",
            controls.len(),
            control_cells
        );
        if tests.len() == 0 {
            eprintln!("Gene scan failed, no test clonotypes.\n");
            std::process::exit(1);
        }
        if controls.len() == 0 {
            eprintln!("Gene scan failed, no control clonotypes.\n");
            std::process::exit(1);
        }
        println!("enriched features\n");
        let mut results = Vec::<(usize, Vec<u8>, f64, f64, f64)>::new();
        let nf = gex_info.gex_features[0].len();
        for fid in 0..nf {
            results.push((fid, Vec::<u8>::new(), 0.0, 0.0, 0.0));
        }
        results.par_iter_mut().for_each(|res| {
            let fid = res.0;
            // NOT SURE THIS IS BACKWARD COMPATIBLE!
            let gene = gex_info.gex_features[0][fid]
                .after("\t")
                .after("\t")
                .contains("Gene");
            let mut test_values = Vec::<f64>::new();
            let mut control_values = Vec::<f64>::new();
            for pass in 1..=2 {
                let tc;
                let vals;
                if pass == 1 {
                    tc = &tests;
                    vals = &mut test_values;
                } else {
                    tc = &controls;
                    vals = &mut control_values;
                }
                for j in 0..tc.len() {
                    for m in 0..exacts[tc[j]].len() {
                        let ex = &exact_clonotypes[exacts[tc[j]][m]];
                        for l in 0..ex.clones.len() {
                            let li = ex.clones[l][0].dataset_index;
                            let bc = ex.clones[l][0].barcode.clone();
                            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                            if p >= 0 {
                                let mut raw_count = 0 as f64;
                                if !ctl.gen_opt.h5 {
                                    raw_count =
                                        gex_info.gex_matrices[li].value(p as usize, fid) as f64;
                                } else {
                                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                                    // p+1 OK?
                                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize;
                                    let d: Vec<u32>;
                                    let ind: Vec<u32>;
                                    if ctl.gen_opt.h5_pre {
                                        d = h5_data[li].1[z1..z2].to_vec();
                                        ind = h5_data[li].2[z1..z2].to_vec();
                                    } else {
                                        d = d_readers[li]
                                            .as_ref()
                                            .unwrap()
                                            .read_slice(&s![z1..z2])
                                            .unwrap()
                                            .to_vec();
                                        ind = ind_readers[li]
                                            .as_ref()
                                            .unwrap()
                                            .read_slice(&s![z1..z2])
                                            .unwrap()
                                            .to_vec();
                                    }
                                    for j in 0..d.len() {
                                        if ind[j] == fid as u32 {
                                            raw_count = d[j] as f64;
                                            break;
                                        }
                                    }
                                }
                                let mult: f64;
                                if gene {
                                    mult = gex_info.gex_mults[li];
                                } else {
                                    mult = gex_info.fb_mults[li];
                                }
                                if !ctl.gen_opt.full_counts {
                                    vals.push(raw_count * mult);
                                } else {
                                    vals.push(raw_count);
                                }
                            }
                        }
                    }
                }
            }
            let mut test_mean = 0.0;
            for i in 0..test_values.len() {
                test_mean += test_values[i];
            }
            test_mean /= test_values.len() as f64;
            let mut control_mean = 0.0;
            for i in 0..control_values.len() {
                control_mean += control_values[i];
            }
            control_mean /= control_values.len() as f64;
            let mut vals = Vec::<f64>::new();
            let threshold = ctl.gen_opt.gene_scan_threshold.clone().unwrap();
            for i in 0..threshold.var.len() {
                if threshold.var[i] == "t".to_string() {
                    vals.push(test_mean);
                } else {
                    vals.push(control_mean);
                }
            }
            if threshold.satisfied(&vals) {
                fwrite!(res.1, "{}", gex_info.gex_features[0][fid]);
                res.2 = test_mean;
                res.3 = control_mean;
                res.4 = test_mean / control_mean;
            }
        });
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "id".to_string(),
            "name".to_string(),
            "library_type".to_string(),
            "test".to_string(),
            "control".to_string(),
            "enrichment".to_string(),
        ];
        rows.push(row);
        for fid in 0..nf {
            if results[fid].1.len() > 0 {
                let stuff = strme(&results[fid].1);
                let fields = stuff.split('\t').collect::<Vec<&str>>();
                let mut row = Vec::<String>::new();
                row.push(fields[0].to_string());
                row.push(fields[1].to_string());
                row.push(fields[2].to_string());
                row.push(format!("{:.2}", results[fid].2));
                row.push(format!("{:.2}", results[fid].3));
                row.push(format!("{:.2}", results[fid].4));
                rows.push(row);
            }
        }
        let mut log = Vec::<u8>::new();
        print_tabular(&mut log, &rows, 2, Some(b"lllrrr".to_vec()));
        print!("{}", strme(&log));
    }

    // Plot clonotypes.

    plot_clonotypes(&ctl, &exacts, &exact_clonotypes);

    // Tally low gene expression count.
    // WARNING: THIS MAY ONLY WORK IF YOU RUN WITH CLONES=1 AND NO OTHER FILTERS.
    // And probably you should run on only one dataset at a time.
    // And this probably doesn't belong inside print_clonotypes.

    if gex_info.gex_features.len() > 0 && !ctl.silent {
        let mut ncells = 0;
        for i in 0..exact_clonotypes.len() {
            ncells += exact_clonotypes[i].clones.len();
        }
        let mut bads = 0;
        for i in 0..results.len() {
            bads += results[i].5;
        }
        let bad_rate = percent_ratio(bads, ncells);
        println!(
            "fraction of cells having gex count < 100 = {:.2}%",
            bad_rate
        );
    }
}
