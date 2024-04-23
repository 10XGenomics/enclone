// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

mod build_table_stuff;
mod finish_table;
mod gene_scan;
mod print_utils1;
mod print_utils2;
mod print_utils3;
mod print_utils4;
mod print_utils5;
mod proc_cvar_auto;
mod proc_lvar2;
mod proc_lvar_auto;

use self::finish_table::{finish_table, Sr};
use self::gene_scan::{gene_scan_test, InSet};
use self::print_utils1::{compute_field_types, extra_args, start_gen};
use self::print_utils2::row_fill;
use self::print_utils3::{consensus_codon_cdr3, get_extra_parseables, process_complete};
use self::print_utils4::{build_show_aa, compute_bu, compute_some_stats, SomeStats};
use self::print_utils5::vars_and_shares;
use crate::mammalian_fixed_len::mammalian_fixed_len_peer_groups;
use enclone_args::proc_args_check::involves_gex_fb;
use enclone_core::allowed_vars::{CVARS_ALLOWED, CVARS_ALLOWED_PCELL, LVARS_ALLOWED};
use enclone_core::defs::ColInfo;
use enclone_core::enclone_structs::{BarcodeFates, EncloneExacts, EncloneSetup, GexReaders};
use enclone_process::process_clonotypes::OrbitProcessor;
use equiv::EquivRel;
use itertools::izip;

use std::collections::{HashMap, HashSet};

use string_utils::TextUtils;

use vector_utils::{bin_member, bin_position, unique_sort};

type Stats = Vec<Vec<(String, Vec<String>)>>;

#[derive(Default)]
pub struct TraverseResult {
    pic: String,
    exacts: Vec<usize>,
    rsi: ColInfo,
    in_center: bool,
    out_data: Vec<HashMap<String, String>>,
    gene_scan_membership: Vec<InSet>,
}

pub struct EncloneOrbitProcessor {
    pub result: PrintClonotypesResult,
    n_vdj_gex: Vec<usize>,
    peer_groups: Vec<Vec<(usize, u8, u32)>>,
    alt_bcs: Vec<String>,
    have_gex: bool,
    need_gex: bool,
    all_vars: Vec<String>,
    extra_args: Vec<String>,
}

impl EncloneOrbitProcessor {
    pub fn new(setup: &EncloneSetup, vdj_cells: &[Vec<String>]) -> Self {
        let EncloneSetup {
            ctl,
            ann: _,
            gex_info,
            tall: _,
            refdata,
        } = setup;

        let lvars = &ctl.clono_print_opt.lvars;

        // Compute extra args.

        let extra_args = extra_args(ctl);

        // Determine if any lvars need gex info.

        let need_gex = {
            lvars.iter().map(String::as_str).any(involves_gex_fb)
                || {
                    if ctl.parseable_opt.pout.is_empty() {
                        false
                    } else if ctl.parseable_opt.pcols.is_empty() {
                        LVARS_ALLOWED.iter().any(|var| involves_gex_fb(var))
                    } else {
                        ctl.parseable_opt
                            .pcols
                            .iter()
                            .map(String::as_str)
                            .any(involves_gex_fb)
                    }
                }
                || extra_args.iter().map(String::as_str).any(involves_gex_fb)
        };

        // Compute all_vars.

        let rsi_vars = &ctl.clono_print_opt.cvars;
        let mut all_vars = rsi_vars.clone();
        for var in CVARS_ALLOWED {
            if !rsi_vars.iter().any(|v| v == var) {
                all_vars.push(var.to_string());
            }
        }
        for var in CVARS_ALLOWED_PCELL {
            if !rsi_vars.iter().any(|v| v == var) {
                all_vars.push(var.to_string());
            }
        }
        // Identify certain extra parseable variables.  These arise from parameterizable cvars.
        all_vars.extend(
            get_extra_parseables(ctl, &ctl.parseable_opt.pcols_sort)
                .into_iter()
                .map(String::from),
        );
        for x in &extra_args {
            if !rsi_vars.contains(x) {
                all_vars.push(x.clone());
            }
        }

        // Test for presence of GEX/FB data.

        let mut have_gex = false;
        for i in 0..ctl.origin_info.gex_path.len() {
            if !ctl.origin_info.gex_path[i].is_empty() {
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

        // Compute peer groups.

        let peer_groups = mammalian_fixed_len_peer_groups(refdata);

        // Compute number of vdj cells that are gex.

        let mut n_vdj_gex = Vec::<usize>::new();
        for (gex, vdj) in gex_info
            .pca
            .iter()
            .zip(vdj_cells.iter())
            .take(ctl.origin_info.n())
        {
            let mut n = 0;
            for y in gex {
                if bin_member(vdj, y.0) {
                    n += 1;
                }
            }
            n_vdj_gex.push(n);
        }

        Self {
            result: Default::default(),
            n_vdj_gex,
            peer_groups,
            alt_bcs,
            have_gex,
            need_gex,
            all_vars,
            extra_args,
        }
    }
}

impl OrbitProcessor<TraverseResult, Stats> for &mut EncloneOrbitProcessor {
    fn filter(
        &self,
        setup: &EncloneSetup,
        enclone_exacts: &EncloneExacts,
        gex_readers: &[Option<GexReaders<'_>>],
        fate: &[BarcodeFates],
        exacts: &[usize],
        mults: &[usize],
        n: usize,
        rsi: &ColInfo,
        bads: &mut [bool],
        in_center: bool,
    ) -> Result<Stats, String> {
        let ctl = &setup.ctl;
        if ctl.clono_filt_opt.bounds.is_empty()
            && !ctl.gen_opt.complete
            && ctl.gen_opt.var_def.is_empty()
        {
            return Ok(Default::default());
        }
        let mut stats_pass1 = Vec::new();
        let res = process_orbit_tail_enclone_only(
            1,
            setup,
            enclone_exacts,
            gex_readers,
            fate,
            &self.all_vars,
            &self.alt_bcs,
            &self.extra_args,
            self.need_gex,
            self.have_gex,
            exacts,
            mults,
            n,
            rsi,
            &self.n_vdj_gex,
            &self.peer_groups,
            bads,
            &mut stats_pass1,
            in_center,
        )?;
        assert!(res.is_none());
        Ok(stats_pass1)
    }

    fn finalize(
        &self,
        setup: &EncloneSetup,
        enclone_exacts: &EncloneExacts,
        gex_readers: &[Option<GexReaders<'_>>],
        fate: &[BarcodeFates],
        exacts: &[usize],
        mults: &[usize],
        n: usize,
        rsi: &ColInfo,
        bads: &mut [bool],
        in_center: bool,
        mut stats_pass1: Stats,
    ) -> Result<Option<TraverseResult>, String> {
        process_orbit_tail_enclone_only(
            2,
            setup,
            enclone_exacts,
            gex_readers,
            fate,
            &self.all_vars,
            &self.alt_bcs,
            &self.extra_args,
            self.need_gex,
            self.have_gex,
            exacts,
            mults,
            n,
            rsi,
            &self.n_vdj_gex,
            &self.peer_groups,
            bads,
            &mut stats_pass1,
            in_center,
        )
    }

    fn collect(&mut self, result: Option<TraverseResult>) {
        if let Some(data) = result {
            self.result.pics.push(data.pic);
            self.result.exacts.push(data.exacts);
            self.result.rsi.push(data.rsi);
            self.result.in_center.push(data.in_center);
            self.result.out_datas.push(data.out_data);
            self.result.gene_scan_result.push(data.gene_scan_membership);
        } else {
            self.result.out_datas.push(Vec::new());
            self.result.gene_scan_result.push(Vec::new());
        }
    }
}

#[derive(Default)]
pub struct PrintClonotypesResult {
    pub pics: Vec<String>,
    pub exacts: Vec<Vec<usize>>,
    pub in_center: Vec<bool>,
    pub rsi: Vec<ColInfo>,
    pub out_datas: Vec<Vec<HashMap<String, String>>>,
    pub gene_scan_result: Vec<Vec<InSet>>,
}

fn process_orbit_tail_enclone_only(
    pass: usize,
    setup: &EncloneSetup,
    enclone_exacts: &EncloneExacts,
    gex_readers: &[Option<GexReaders<'_>>],
    fate: &[BarcodeFates],
    all_vars: &[String],
    alt_bcs: &[String],
    extra_args: &[String],
    need_gex: bool,
    have_gex: bool,
    exacts: &[usize],
    mults: &[usize],
    n: usize,
    rsi: &ColInfo,
    n_vdj_gex: &[usize],
    peer_groups: &[Vec<(usize, u8, u32)>],
    bads: &mut [bool],
    stats_pass1: &mut Vec<Vec<(String, Vec<String>)>>,
    mut in_center: bool,
) -> Result<Option<TraverseResult>, String> {
    let EncloneSetup {
        ctl,
        ann: _,
        gex_info,
        tall: _,
        refdata,
    } = setup;
    let EncloneExacts {
        to_bc: _,
        exact_clonotypes,
        raw_joins: _,
        info: _,
        orbits: _,
        vdj_cells,
        join_info: _,
        drefs: dref,
        sr: _,
        allele_data,
    } = enclone_exacts;
    let nexacts = exacts.len();
    let mat = &rsi.mat;
    let cols = mat.len();
    let lvars = &ctl.clono_print_opt.lvars;

    let mut mlog = Vec::<u8>::new();
    let mut out_data = Vec::new();

    // Start to generate parseable output.

    // Print the orbit.
    // ◼ An assumption of this code is that a productive pair does not have two contigs
    // ◼ having identical CDR3_AA sequences.  At present this is not enforced by the
    // ◼ assembly stage, so the assumption is violated.  To work around this there are
    // ◼ some unsavory workarounds below.

    if pass == 2 {
        start_gen(
            ctl,
            exacts,
            exact_clonotypes,
            &mut out_data,
            &mut mlog,
            extra_args,
        );
    }

    // Find variant positions.  And some other things.

    let mut vars = Vec::<Vec<usize>>::new();
    let mut vars_amino = Vec::<Vec<usize>>::new();
    let mut shares_amino = Vec::<Vec<usize>>::new();
    let mut ref_diff_pos = Vec::<Vec<Vec<usize>>>::new();
    vars_and_shares(
        pass,
        ctl,
        exacts,
        exact_clonotypes,
        rsi,
        refdata,
        dref,
        &mut vars,
        &mut vars_amino,
        &mut shares_amino,
        &mut ref_diff_pos,
        &mut out_data,
    );

    // Define amino acid positions to show.

    let show_aa = build_show_aa(
        ctl,
        rsi,
        &vars_amino,
        &shares_amino,
        refdata,
        dref,
        exacts,
        exact_clonotypes,
    );

    // Define field types corresponding to the amino acid positions to show.
    let field_types = compute_field_types(ctl, rsi, &show_aa);

    // Build varmat matrix of size (nexacts, cols).
    let mut varmat = vec![vec![vec![b'-']; cols]; nexacts];
    for (col, (mat_slice, seqss_slice, vars_slice)) in
        izip!(mat, &rsi.seqss, &vars).take(cols).enumerate()
    {
        for (varmat_u, m, seq) in izip!(&mut varmat, mat_slice, seqss_slice) {
            varmat_u[col] = if m.is_some() {
                vars_slice
                    .iter()
                    .map(|&p| *seq.get(p).unwrap_or(&b'?'))
                    .collect()
            } else {
                vec![b'-']
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
            lvarsc.extend(lvars.iter().take(i).cloned());
            let k = x.after("nd").force_usize();
            let mut n = vec![0_usize; ctl.origin_info.n()];
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
            for (l, ds) in datasets.iter().take(k).enumerate() {
                if l >= n.len() {
                    break;
                }
                nd_fields.push(format!("n_{}", ds.as_str()));
                lvarsc.push(format!("n_{}", ds.as_str()));
            }
            if n.len() > k {
                nd_fields.push("n_other".to_string());
                lvarsc.push("n_other".to_string());
            }
            lvarsc.extend(lvars.iter().skip(i + 1).cloned());
            break;
        }
    }
    let lvars = lvarsc.clone();
    let mut lvarsh = HashSet::<String>::new();
    for x in &lvars {
        lvarsh.insert(x.to_string());
    }

    // Now build table content.

    let mut sr = Vec::<Sr>::new();
    let mut groups = HashMap::<usize, Vec<usize>>::new();
    for lvar in &lvars {
        if let Some(Ok(d)) = lvar.strip_prefix('g').map(str::parse::<usize>) {
            if groups.contains_key(&d) {
                continue;
            }
            let mut e: EquivRel = EquivRel::new(nexacts as i32);
            for (u1, &e1) in exacts.iter().take(nexacts).enumerate() {
                let ex1 = &exact_clonotypes[e1];
                for (u2, &e2) in exacts.iter().enumerate().take(nexacts).skip(u1 + 1) {
                    if e.class_id(u1 as i32) == e.class_id(u2 as i32) {
                        continue;
                    }
                    let ex2 = &exact_clonotypes[e2];
                    let mut diffs = 0;
                    'comp: for (mm, vars) in mat.iter().zip(vars.iter()).take(cols) {
                        if let (Some(m1), Some(m2)) = (mm[u1], mm[u2]) {
                            let (s1, s2) = (&ex1.share[m1].seq_del, &ex2.share[m2].seq_del);
                            for &p in vars {
                                if s1[p] != s2[p] {
                                    diffs += 1;
                                    if diffs > d {
                                        break 'comp;
                                    }
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

    let mut stats = Vec::<(String, Vec<String>)>::new();

    // Compute some stats;

    let SomeStats { cred, pe, ppe, npe } = compute_some_stats(
        ctl,
        &lvars,
        exacts,
        exact_clonotypes,
        gex_info,
        vdj_cells,
        n_vdj_gex,
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

    // Form CDR3 consensus sequences.

    let mut cdr3_con = Vec::<Vec<u8>>::new();
    if ctl.gen_opt.color == "codon-diffs" {
        cdr3_con = consensus_codon_cdr3(rsi, exacts, exact_clonotypes);
    }

    // Build rows.

    let mut cell_count = 0;
    for u in 0..nexacts {
        let mut typex = vec![false; cols];
        let mut row = Vec::<String>::new();
        let mut cx = Vec::<Vec<String>>::new();
        for col in 0..cols {
            cx.push(vec![String::new(); rsi.cvars[col].len()]);
        }
        let clonotype_id = exacts[u];
        let ex = &exact_clonotypes[clonotype_id];
        let mut d_all = vec![Vec::<u32>::new(); ex.clones.len()];
        let mut ind_all = vec![Vec::<u32>::new(); ex.clones.len()];
        let mut these_stats = Vec::<(String, Vec<String>)>::new();
        row_fill(
            pass,
            u,
            ctl,
            exacts,
            mults,
            exact_clonotypes,
            gex_info,
            refdata,
            &varmat,
            &fp,
            &vars_amino,
            &show_aa,
            &ref_diff_pos,
            &field_types,
            bads,
            &mut row,
            &mut out_data,
            &mut cx,
            &mut d_all,
            &mut ind_all,
            rsi,
            dref,
            &groups,
            gex_readers,
            &mut these_stats,
            stats_pass1,
            vdj_cells,
            n_vdj_gex,
            &lvars,
            &lvarsh,
            &nd_fields,
            peer_groups,
            extra_args,
            all_vars,
            need_gex,
            fate,
            &cdr3_con,
            allele_data,
        )?;
        stats.append(&mut these_stats.clone());
        if pass == 1 {
            stats_pass1.push(these_stats.clone());
        }
        these_stats.sort_by(|a, b| a.0.cmp(&b.0));
        let mut bli = ex
            .clones
            .iter()
            .enumerate()
            .map(|(l, clone)| (clone[0].barcode.clone(), clone[0].dataset_index, l))
            .collect::<Vec<_>>();
        // WHY ARE WE SORTING HERE?
        bli.sort();
        for col in 0..cols {
            if mat[col][u].is_some() {
                typex[col] = true;
            }
        }
        for mut cxr in cx {
            row.append(&mut cxr);
        }

        // Compute per-cell entries.

        if pass == 2 {
            let mut subrows = Vec::<Vec<String>>::new();
            compute_bu(
                u,
                cell_count,
                exacts,
                &lvars,
                ctl,
                &bli,
                ex,
                exact_clonotypes,
                &mut row,
                &mut subrows,
                have_gex,
                gex_info,
                rsi,
                &mut sr,
                fate,
                &nd_fields,
                alt_bcs,
                &cred,
                &pe,
                &ppe,
                &npe,
                &d_all,
                &ind_all,
                mat,
                &these_stats,
                refdata,
            );
        }
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
    let mut stats2 = Vec::<(String, Vec<String>)>::new();
    let mut i = 0;
    while i < stats.len() {
        let mut j = i + 1;
        while j < stats.len() {
            if stats[j].0 != stats[i].0 {
                break;
            }
            j += 1;
        }
        let all = stats[i..j]
            .iter()
            .flat_map(|s| s.1.iter().cloned())
            .collect();
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
        let mut mins = Vec::<f64>::new();
        let mut maxs = Vec::<f64>::new();
        // traverse the coefficients on the left hand side (each having a variable)
        let mut fail = false;
        for vi in x.var.iter().take(x.n()) {
            let mut vals = Vec::<f64>::new(); // the stats for the variable
            for stat in &stats {
                if stat.0 == *vi {
                    for sk in &stat.1 {
                        if let Ok(sk) = sk.parse::<f64>() {
                            vals.push(sk);
                        }
                    }
                    break;
                }
            }
            let mut min = 1_000_000_000.0_f64;
            let mut mean = 0.0;
            let mut max = -1_000_000_000.0_f64;
            let mut count = 0;
            for val in vals {
                if !val.is_nan() {
                    min = min.min(val);
                    mean += val;
                    max = max.max(val);
                    count += 1;
                }
            }
            if count == 0 {
                fail = true;
            } else {
                mins.push(min);
                mean /= count as f64;
                means.push(mean);
                maxs.push(max);
            }
        }
        if ctl.clono_filt_opt.bound_type[bi] == "mean" && (fail || !x.satisfied(&means)) {
            if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                in_center = false;
            } else {
                for b in bads.iter_mut().take(nexacts) {
                    *b = true;
                }
            }
        }
        if ctl.clono_filt_opt.bound_type[bi] == "min" && (fail || !x.satisfied(&mins)) {
            if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                in_center = false;
            } else {
                for b in bads.iter_mut().take(nexacts) {
                    *b = true;
                }
            }
        }
        if ctl.clono_filt_opt.bound_type[bi] == "max" && (fail || !x.satisfied(&maxs)) {
            if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                in_center = false;
            } else {
                for b in bads.iter_mut().take(nexacts) {
                    *b = true;
                }
            }
        }
    }

    // Process COMPLETE.

    process_complete(ctl, nexacts, bads, mat);

    // Done unless on second pass.

    if pass == 1 {
        return Ok(None);
    }

    // See if we're in the test and control sets for gene scan.
    let gene_scan_membership = ctl
        .gen_opt
        .gene_scan
        .as_ref()
        .map(|gene_scan_opts| {
            gene_scan_test(
                gene_scan_opts,
                ctl.gen_opt.gene_scan_exact,
                &stats,
                &stats_orig,
                nexacts,
                n,
            )
        })
        .unwrap_or_default();

    // Make the table.

    let clonotype_pic = finish_table(
        n,
        ctl,
        exacts,
        exact_clonotypes,
        rsi,
        &vars,
        &show_aa,
        &field_types,
        &lvars,
        refdata,
        dref,
        peer_groups,
        &mut mlog,
        &stats,
        sr,
        extra_args,
        &mut out_data,
        &rord,
        pass,
        &cdr3_con,
    );

    Ok(Some(TraverseResult {
        pic: clonotype_pic,
        exacts: exacts.to_vec(),
        rsi: rsi.clone(),
        in_center,
        out_data,
        gene_scan_membership,
    }))
}
