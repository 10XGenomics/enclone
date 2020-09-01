// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;
use vector_utils::*;

// Exploratory code, in which we examine the data, so as to help decide on the
// best definition of exact subclonotype.  Super antiquated.

pub fn explore(li: usize, tig_bc: &Vec<Vec<TigData>>, ctl: &EncloneControl) {
    if ctl.gen_opt.exp {
        let mut logs = Vec::<Vec<u8>>::new();
        let mut r = 0;
        while r < tig_bc.len() {
            let mut s = r + 1;
            while s < tig_bc.len() {
                if tig_bc[s].len() != tig_bc[r].len() {
                    break;
                }
                let mut ok = true;
                for m in 0..tig_bc[r].len() {
                    if tig_bc[s][m].cdr3_dna != tig_bc[r][m].cdr3_dna
                        || tig_bc[s][m].len != tig_bc[r][m].len
                    {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    break;
                }
                s += 1;
            }
            if s > r {
                if s - r > 1 {
                    // if only one barcode, should save but not analyze
                    let mut log = Vec::<u8>::new();
                    fwriteln!(log, "lena = {}", ctl.origin_info.dataset_id[li]);
                    fwriteln!(log, "there are {} barcodes", s - r);

                    // printme!( r, s ); // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                    /*
                    for t in r..s { // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                        println!( "bc {}", tig_bc[t][0].barcode ); // YYYYYYYYYYYYYYYYYYYYY
                    } // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                    */

                    fwriteln!(log, "there are {} chains", tig_bc[r].len());
                    for m in 0..tig_bc[r].len() {
                        fwrite!(
                            log,
                            "{} = {} = {} = {}",
                            m + 1,
                            tig_bc[r][m].cdr3_aa,
                            tig_bc[r][m].cdr3_dna,
                            tig_bc[r][m].len
                        );
                        if m < tig_bc[r].len() - 1 {
                            fwriteln!(log, " ");
                        }
                    }
                    fwriteln!(log, "");
                    let mut bc_print = Vec::<usize>::new();
                    // let mut different_columns = 0;
                    // let mut clear_cut = true;
                    let mut printme = false;
                    let mut bcs = Vec::<(usize, usize)>::new();
                    for m in 0..tig_bc[r].len() {
                        // go through the chains
                        for p in 0..tig_bc[r][m].len {
                            // go through positions

                            // Find {(base,qual,bc)}.

                            let mut bqb = Vec::<(u8, u8, usize)>::new();
                            for u in r..s {
                                bqb.push((tig_bc[u][m].seq[p], tig_bc[u][m].quals[p], u - r));
                            }
                            bqb.sort();

                            // Ignore case where only one base seen.

                            let mut nbases = 1;
                            for j in 1..bqb.len() {
                                if bqb[j].0 != bqb[j - 1].0 {
                                    nbases += 1;
                                }
                            }
                            if nbases == 1 {
                                continue;
                            }

                            // Gather frequency info.

                            let mut freq = Vec::<(usize, u8, Vec<(u8, usize)>)>::new();
                            let mut u = 0;
                            while u < bqb.len() {
                                let mut v = u + 1;
                                while v < bqb.len() {
                                    if bqb[v].0 != bqb[u].0 {
                                        break;
                                    }
                                    v += 1;
                                }
                                let mut y = Vec::<(u8, usize)>::new();
                                for t in u..v {
                                    y.push((bqb[t].1, bqb[t].2));
                                }
                                freq.push((y.len(), bqb[u].0, y));
                                u = v;
                            }
                            reverse_sort(&mut freq);
                            /*
                            if freq.len() > 1 {
                                different_columns += 1;
                            }
                            */

                            // A column is declared "clear cut" if there is at most one
                            // base having only non-Q60 support.

                            /*
                            let mut non_q60s = 0;
                            for m in 0..freq.len() {
                                let mut have_q60 = false;
                                for j in 0..freq[m].2.len() {
                                    if freq[m].2[j].0 >= 60 {
                                        have_q60 = true;
                                    }
                                }
                                if !have_q60 {
                                    non_q60s += 1;
                                }
                            }
                            if non_q60s > 1 {
                                clear_cut = false;
                            }
                            */

                            // Print.

                            fwrite!(log, "chain {}, pos {}:", m + 1, p + 1);
                            for f in 0..freq.len() {
                                if freq[f].0 > 1 && f == 0 {
                                    fwrite!(log, " {}[{}]", freq[f].1 as char, freq[f].0);
                                } else {
                                    fwrite!(log, " {}[{}", freq[f].1 as char, freq[f].0);
                                    for j in 0..freq[f].2.len() {
                                        let bc_id = freq[f].2[j].1;
                                        bcs.push((bc_id, m));
                                        bc_print.push(bc_id);
                                        fwrite!(log, "; q = {}, bc = {}", freq[f].2[j].0, bc_id);
                                    }
                                    fwrite!(log, "]");
                                }
                                let mut have_q60 = false;
                                for j in 0..freq[f].2.len() {
                                    if freq[f].2[j].0 >= 60 {
                                        have_q60 = true;
                                    }
                                }
                                if !have_q60 {
                                    fwrite!(log, " = WEAK");
                                }
                            }
                            fwriteln!(log, "");
                        }
                        bcs.sort();
                        let mut bcs_count = Vec::<(usize, usize)>::new();
                        let mut i = 0;
                        while i < bcs.len() {
                            let j = next_diff(&bcs, i);
                            bcs_count.push((bcs[i].0, j - i));
                            i = j;
                        }
                        let mut i = 0;
                        while i < bcs_count.len() {
                            let j = next_diff1_2(&bcs_count, i as i32) as usize;
                            if j - i == 1 && bcs_count[i].1 >= 10 {
                                printme = true;
                            }
                            i = j;
                        }
                    }
                    // if different_columns > 5 || !clear_cut {
                    if printme || ctl.gen_opt.weak {
                        unique_sort(&mut bc_print);
                        for j in 0..bc_print.len() {
                            let bc_id = bc_print[j];
                            fwriteln!(log, "bc {} = {}", bc_id, tig_bc[bc_id + r][0].barcode);
                        }

                        logs.push(log);
                    }
                }
            }
            r = s;
        }
        for i in 0..logs.len() {
            println!("\ncase {}", i + 1);
            print!("{}", strme(&logs[i]));
        }
        std::process::exit(0);
    }
}

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
