// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Potential join structure.

pub struct PotentialJoin {
    pub k1: usize,
    pub k2: usize,
    pub nrefs: usize,
    pub cd: isize,
    pub diffs: usize,
    pub bcs1: Vec<String>,
    pub bcs2: Vec<String>,
    pub shares: Vec<isize>,
    pub indeps: Vec<isize>,
    pub shares_details: Vec<Vec<usize>>,
    pub share_pos_v: Vec<Vec<usize>>,
    pub share_pos_j: Vec<Vec<usize>>,
    pub score: f64,
    pub err: bool,
    pub p1: f64,
    pub mult: f64,
}
