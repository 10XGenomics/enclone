// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use amino::nucleotide_to_aminoacid_sequence;
use std::collections::HashMap;
use string_utils::TextUtils;
use superslice::Ext;
use vdj_ann::refx::RefData;
use vdj_ann::vdj_features::{cdr1_start, cdr2_start, cdr3_start, fr1_start, fr2_start, fr3_start};

// {chain, feature, len, {{(count, amino_acid)}}}

#[allow(clippy::type_complexity)]
pub fn mammalian_fixed_len() -> Vec<(&'static str, &'static str, usize, Vec<Vec<(u32, u8)>>)> {
    const X: &str = include_str!("mammalian_fixed_len.table");
    X.lines()
        .map(|line| {
            let mut y = line.splitn(5, ',');
            (
                y.next().unwrap(),
                y.next().unwrap(),
                y.next().unwrap().force_usize(),
                y.next()
                    .unwrap()
                    .split('+')
                    .map(|zi| {
                        zi.split('/')
                            .map(|wj| {
                                (
                                    wj.before(":").force_usize() as u32,
                                    wj.after(":").as_bytes()[0],
                                )
                            })
                            .collect()
                    })
                    .collect(),
            )
        })
        .collect()
}

// Calculate peer groups for each V segment reference sequence.

pub fn mammalian_fixed_len_peer_groups(refdata: &RefData) -> Vec<Vec<(usize, u8, u32)>> {
    let m = mammalian_fixed_len();
    let mut start = HashMap::<(&str, &str), usize>::new();
    let mut stop = HashMap::<(&str, &str), usize>::new();
    for chain in ["IGH", "IGK", "IGL", "TRA", "TRB"] {
        for feature in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"] {
            let low = m.lower_bound_by_key(&(chain, feature), |(a, b, _c, _d)| (a, b));
            let high = m.upper_bound_by_key(&(chain, feature), |(a, b, _c, _d)| (a, b));
            start.insert((chain, feature), low);
            stop.insert((chain, feature), high);
        }
    }
    let mut pg = vec![Vec::<(usize, u8, u32)>::new(); refdata.refs.len()];
    for (i, pgi) in pg.iter_mut().enumerate() {
        if refdata.is_v(i) {
            let aa = nucleotide_to_aminoacid_sequence(&refdata.refs[i].to_ascii_vec(), 0);
            let rtype = refdata.rtype[i];
            let chain_type = if rtype == 0 {
                "IGH"
            } else if rtype == 1 {
                "IGK"
            } else if rtype == 2 {
                "IGL"
            } else if rtype == 3 {
                "TRA"
            } else if rtype == 4 {
                "TRB"
            } else {
                continue;
            };
            let fs1 = fr1_start(&aa, chain_type);
            let cs1 = cdr1_start(&aa, chain_type, false);
            let fs2 = fr2_start(&aa, chain_type, false);
            let cs2 = cdr2_start(&aa, chain_type, false);
            let fs3 = fr3_start(&aa, chain_type, false);
            let cs3 = cdr3_start(&aa, chain_type, false);
            if let Some(cs1) = cs1 {
                if fs1 < cs1 {
                    let x1 = start[&(chain_type, "fwr1")];
                    let x2 = stop[&(chain_type, "fwr1")];
                    let len = cs1 - fs1;
                    for mx in &m[x1..x2] {
                        if mx.2 == len {
                            for (j, mj) in mx.3.iter().take(len).enumerate() {
                                pgi.extend(mj.iter().map(|mjk| (fs1 + j, mjk.1, mjk.0)));
                            }
                        }
                    }
                }
            }
            if let (Some(cs1), Some(fs2)) = (cs1, fs2) {
                if cs1 < fs2 {
                    let x1 = start[&(chain_type, "cdr1")];
                    let x2 = stop[&(chain_type, "cdr1")];
                    let len = fs2 - cs1;
                    for mx in &m[x1..x2] {
                        if mx.2 == len {
                            for (j, mj) in mx.3.iter().take(len).enumerate() {
                                pgi.extend(mj.iter().map(|mjk| (cs1 + j, mjk.1, mjk.0)));
                            }
                        }
                    }
                }
            }
            if let (Some(fs2), Some(cs2)) = (fs2, cs2) {
                if fs2 < cs2 {
                    let x1 = start[&(chain_type, "fwr2")];
                    let x2 = stop[&(chain_type, "fwr2")];
                    let len = cs2 - fs2;
                    for mx in &m[x1..x2] {
                        if mx.2 == len {
                            for (j, mj) in mx.3.iter().take(len).enumerate() {
                                pgi.extend(mj.iter().map(|mjk| (fs2 + j, mjk.1, mjk.0)));
                            }
                        }
                    }
                }
            }
            if let (Some(cs2), Some(fs3)) = (cs2, fs3) {
                if cs2 < fs3 {
                    let x1 = start[&(chain_type, "cdr2")];
                    let x2 = stop[&(chain_type, "cdr2")];
                    let len = fs3 - cs2;
                    for mx in &m[x1..x2] {
                        if mx.2 == len {
                            for (j, mj) in mx.3.iter().take(len).enumerate() {
                                pgi.extend(mj.iter().map(|mjk| (cs2 + j, mjk.1, mjk.0)));
                            }
                        }
                    }
                }
            }
            if let Some(fs3) = fs3 {
                if fs3 < cs3 {
                    let x1 = start[&(chain_type, "fwr3")];
                    let x2 = stop[&(chain_type, "fwr3")];
                    let len = cs3 - fs3;
                    for mx in &m[x1..x2] {
                        if mx.2 == len {
                            for (j, mj) in mx.3.iter().take(len).enumerate() {
                                pgi.extend(mj.iter().map(|mjk| (fs3 + j, mjk.1, mjk.0)));
                            }
                        }
                    }
                }
            }
        }
    }
    pg
}
