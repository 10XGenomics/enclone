// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use vector_utils::*;

// Look for insertions (initial exploration).

pub fn find_insertions(ctl: &EncloneControl, exact_clonotypes: &Vec<ExactClonotype>) {
    if ctl.gen_opt.insertions {
        println!("CDR3s associated with possible SHM insertions");
        let mut z = Vec::<(String, usize, isize)>::new(); // {(cdr3_aa, v_ref_id, delta)}
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            for j in 0..ex.share.len() {
                let sh = &ex.share[j];
                if sh.annv.len() == 2 && sh.annv[1].0 > sh.annv[0].0 + sh.annv[0].1 {
                    let ins = sh.annv[1].0 - sh.annv[0].0 - sh.annv[0].1;
                    z.push((sh.cdr3_aa.clone(), sh.v_ref_id, ins as isize));
                } else if sh.annv.len() == 1 {
                    z.push((sh.cdr3_aa.clone(), sh.v_ref_id, 0));
                }
            }
        }
        unique_sort(&mut z);
        z.sort();
        let mut i = 0;
        while i < z.len() {
            let j = next_diff12_3(&z, i as i32) as usize;
            if j - i > 1 {
                let mut have_zero = false;
                for k in i..j {
                    if z[k].2 == 0 {
                        have_zero = true;
                    }
                }
                if have_zero {
                    for k in i..j {
                        if z[k].2 > 0 {
                            println!("{} ==> {}", z[k].0, z[k].2);
                        }
                    }
                }
            }
            i = j;
        }
        println!("");
        std::process::exit(0);
    }
}
