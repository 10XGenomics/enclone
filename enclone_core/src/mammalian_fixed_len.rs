// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::vdj_features::{cdr1_start, cdr2_start, cdr3_start, fr1_start, fr2_start, fr3_start};
use amino::aa_seq;
use std::collections::HashMap;
use string_utils::TextUtils;
use superslice::Ext;
use vdj_ann::refx::RefData;

// {chain, feature, len, {{(count, amino_acid)}}}

pub fn mammalian_fixed_len() -> Vec<(String, String, usize, Vec<Vec<(u32, u8)>>)> {
    let mut p = Vec::<(String, String, usize, Vec<Vec<(u32, u8)>>)>::new();
    let x = include_str!("mammalian_fixed_len.table");
    for line in x.lines() {
        let y = line.split(',').collect::<Vec<&str>>();
        let mut pp = Vec::<Vec<(u32, u8)>>::new();
        let z = y[3].split('+').collect::<Vec<&str>>();
        for i in 0..z.len() {
            let mut ppp = Vec::<(u32, u8)>::new();
            let w = z[i].split('/').collect::<Vec<&str>>();
            for j in 0..w.len() {
                ppp.push((
                    w[j].before(":").force_usize() as u32,
                    w[j].after(":").as_bytes()[0] as u8,
                ));
            }
            pp.push(ppp);
        }
        p.push((y[0].to_string(), y[1].to_string(), y[2].force_usize(), pp));
    }
    p
}

// Calculate peer groups for each V segment reference sequence.

pub fn mammalian_fixed_len_peer_groups(refdata: &RefData) -> Vec<Vec<(usize, u8, u32)>> {
    let mut pg = vec![Vec::<(usize, u8, u32)>::new(); refdata.refs.len()];
    let m = mammalian_fixed_len();
    let mut start = HashMap::<(String, String), usize>::new();
    let mut stop = HashMap::<(String, String), usize>::new();
    for chain in ["IGH", "IGK", "IGL", "TRA", "TRB"].iter() {
        for feature in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"].iter() {
            let low = m.lower_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            let high = m.upper_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            start.insert((chain.to_string(), feature.to_string()), low);
            stop.insert((chain.to_string(), feature.to_string()), high);
        }
    }
    for i in 0..refdata.refs.len() {
        if refdata.is_v(i) {
            let aa = aa_seq(&refdata.refs[i].to_ascii_vec(), 0);
            let rtype = refdata.rtype[i];
            let chain_type;
            if rtype == 0 {
                chain_type = "IGH";
            } else if rtype == 1 {
                chain_type = "IGK";
            } else if rtype == 2 {
                chain_type = "IGL";
            } else if rtype == 3 {
                chain_type = "TRA";
            } else if rtype == 4 {
                chain_type = "TRB";
            } else {
                continue;
            }
            let fs1 = fr1_start(&aa, chain_type);
            let cs1 = cdr1_start(&aa, chain_type, false);
            let fs2 = fr2_start(&aa, chain_type, false);
            let cs2 = cdr2_start(&aa, chain_type, false);
            let fs3 = fr3_start(&aa, chain_type, false);
            let cs3 = cdr3_start(&aa, chain_type, false);
            if cs1.is_some() && fs1 < cs1.unwrap() {
                let cs1 = cs1.unwrap();
                let x1 = start[&(chain_type.to_string(), "fwr1".to_string())];
                let x2 = stop[&(chain_type.to_string(), "fwr1".to_string())];
                let len = cs1 - fs1;
                for x in x1..x2 {
                    if m[x].2 == len {
                        for j in 0..len {
                            for k in 0..m[x].3[j].len() {
                                pg[i].push((fs1 + j, m[x].3[j][k].1, m[x].3[j][k].0));
                            }
                        }
                    }
                }
            }
            if cs1.is_some() && fs2.is_some() && cs1.unwrap() < fs2.unwrap() {
                let cs1 = cs1.unwrap();
                let fs2 = fs2.unwrap();
                let x1 = start[&(chain_type.to_string(), "cdr1".to_string())];
                let x2 = stop[&(chain_type.to_string(), "cdr1".to_string())];
                let len = fs2 - cs1;
                for x in x1..x2 {
                    if m[x].2 == len {
                        for j in 0..len {
                            for k in 0..m[x].3[j].len() {
                                pg[i].push((cs1 + j, m[x].3[j][k].1, m[x].3[j][k].0));
                            }
                        }
                    }
                }
            }
            if fs2.is_some() && cs2.is_some() && fs2.unwrap() < cs2.unwrap() {
                let fs2 = fs2.unwrap();
                let cs2 = cs2.unwrap();
                let x1 = start[&(chain_type.to_string(), "fwr2".to_string())];
                let x2 = stop[&(chain_type.to_string(), "fwr2".to_string())];
                let len = cs2 - fs2;
                for x in x1..x2 {
                    if m[x].2 == len {
                        for j in 0..len {
                            for k in 0..m[x].3[j].len() {
                                pg[i].push((fs2 + j, m[x].3[j][k].1, m[x].3[j][k].0));
                            }
                        }
                    }
                }
            }
            if cs2.is_some() && fs3.is_some() && cs2.unwrap() < fs3.unwrap() {
                let cs2 = cs2.unwrap();
                let fs3 = fs3.unwrap();
                let x1 = start[&(chain_type.to_string(), "cdr2".to_string())];
                let x2 = stop[&(chain_type.to_string(), "cdr2".to_string())];
                let len = fs3 - cs2;
                for x in x1..x2 {
                    if m[x].2 == len {
                        for j in 0..len {
                            for k in 0..m[x].3[j].len() {
                                pg[i].push((cs2 + j, m[x].3[j][k].1, m[x].3[j][k].0));
                            }
                        }
                    }
                }
            }
            if fs3.is_some() && fs3.unwrap() < cs3 {
                let fs3 = fs3.unwrap();
                let x1 = start[&(chain_type.to_string(), "fwr3".to_string())];
                let x2 = stop[&(chain_type.to_string(), "fwr3".to_string())];
                let len = cs3 - fs3;
                for x in x1..x2 {
                    if m[x].2 == len {
                        for j in 0..len {
                            for k in 0..m[x].3[j].len() {
                                pg[i].push((fs3 + j, m[x].3[j][k].1, m[x].3[j][k].0));
                            }
                        }
                    }
                }
            }
        }
    }
    pg
}
