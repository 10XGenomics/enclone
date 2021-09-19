// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::doublets::*;
use crate::merge_onesies::*;
use crate::split_orbits::*;
use crate::weak_chains::*;
use enclone_core::defs::*;
use equiv::EquivRel;
use std::collections::HashMap;
use std::time::Instant;

pub fn some_filters(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    eq: &EquivRel,
    disintegrated: &Vec<bool>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    // Delete exact subclonotypes that appear to represent doublets.

    let tdoublet = Instant::now();
    delete_doublets(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tdoublet, "doublet filtering");

    // Merge onesies where totally unambiguous.

    let tmerge = Instant::now();
    merge_onesies(orbits, &ctl, &exact_clonotypes, &info, &eq, &disintegrated);
    ctl.perf_stats(&tmerge, "merging onesies");

    // Check for disjoint orbits.

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 1");

    // Test for weak chains.

    let tweak = Instant::now();
    weak_chains(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
        fate,
    );
    ctl.perf_stats(&tweak, "weak chain filtering");

    // Check for disjoint orbits (again).

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 2");
}
