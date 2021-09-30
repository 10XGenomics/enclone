// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use string_utils::TextUtils;

pub fn mammalian_pwms() -> Vec<(String, String, usize, Vec<Vec<(u32, u8)>>)> {
    let mut p = Vec::<(String, String, usize, Vec<Vec<(u32, u8)>>)>::new();
    let x = include_str!("mammalian_pwms.table");
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
