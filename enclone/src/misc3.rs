// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use enclone_core::defs::*;
use io_utils::*;
use itertools::*;
use std::cmp::{max, min, Ordering};
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn sort_tig_bc(ctl: &EncloneControl, tig_bc: &mut Vec<Vec<TigData>>, refdata: &RefData) {
    tig_bc.sort_by(|x, y| -> Ordering {
        for i in 0..x.len() {
            // Order by number of chains.

            if i >= y.len() {
                return Ordering::Greater;
            }

            // Order by cdr3_dna.

            if x[i].cdr3_dna < y[i].cdr3_dna {
                return Ordering::Less;
            } else if x[i].cdr3_dna > y[i].cdr3_dna {
                return Ordering::Greater;

            // Order by chain length.
            } else if x[i].len < y[i].len {
                return Ordering::Less;
            } else if x[i].len > y[i].len {
                return Ordering::Greater;

            // Order by chain sequence.
            } else if x[i].seq() < y[i].seq() {
                return Ordering::Less;
            } else if x[i].seq() > y[i].seq() {
                return Ordering::Greater;
            }

            // Working around a bug here and below.  For TCR, there are two TRBC1 records in the,
            // reference, having a SNP after our primer, and we appear to pick one of
            // them at random.  Also not sure this fully respects the sort order.
            // And of course a customer could have the same feature in their reference.

            let (cid1, cid2) = (x[i].c_ref_id, y[i].c_ref_id);
            if cid1.is_none() && cid2.is_some() {
                return Ordering::Less;
            } else if cid2.is_none() && cid1.is_some() {
                return Ordering::Greater;

            // Order by constant region name.
            } else if cid1.is_some()
                && cid2.is_some()
                && refdata.name[cid1.unwrap()] < refdata.name[cid2.unwrap()]
            {
                return Ordering::Less;
            } else if cid1.is_some()
                && cid2.is_some()
                && refdata.name[cid1.unwrap()] > refdata.name[cid2.unwrap()]
            {
                return Ordering::Greater;

            // Order by JC delta.
            } else if x[i].c_start.is_some()
                && y[i].c_start.is_some()
                && x[i].c_start.unwrap() + y[i].j_stop < y[i].c_start.unwrap() + x[i].j_stop
            {
                return Ordering::Less;
            } else if x[i].c_start.is_some()
                && y[i].c_start.is_some()
                && x[i].c_start.unwrap() + y[i].j_stop > y[i].c_start.unwrap() + x[i].j_stop
            {
                return Ordering::Greater;

            // Order by donor if MIX_DONORS option used.
            } else if !ctl.clono_filt_opt.donor && x[i].donor_index < y[i].donor_index {
                return Ordering::Less;
            } else if !ctl.clono_filt_opt.donor && x[i].donor_index > y[i].donor_index {
                return Ordering::Greater;
            }
        }
        if x.len() < y.len() {
            return Ordering::Less;
        }
        return Ordering::Equal;
    });
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Exploratory code to understand exact subclonotype consensus formation.
//
// 1. Not clear if we want to do this for exact subclonotypes or for clonotypes.
// 2. Or should we be computing by UTR?
// 3. And should we set this aside as not needing to be solved now? ***************
//
// - might be nice to show forward instead of reverse
// - don't show the most trivial case with one UTR and all agree
// - find code simplifications.

pub fn study_consensus(
    _count: &mut usize,
    ctl: &EncloneControl,
    share: &Vec<TigData1>,
    clones: &Vec<Vec<TigData0>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
) {
    if ctl.gen_opt.utr_con {
        for z in 0..clones[0].len() {
            let mut log = Vec::<u8>::new();
            let (mut utr_ids, mut v_ids) = (Vec::<usize>::new(), Vec::<usize>::new());
            for _ in 0..clones.len() {
                if share[z].u_ref_id.is_some() {
                    utr_ids.push(share[z].u_ref_id.unwrap());
                }
                v_ids.push(share[z].v_ref_id);
            }
            unique_sort(&mut utr_ids);
            unique_sort(&mut v_ids);
            fwriteln!(
                log,
                "\n[{}] rev lefts for exact subclonotype {}, chain {}, \
                 vs = {}, utrs = {}, cdr3 = {}\n",
                *_count + 1,
                exact_clonotypes.len(),
                z,
                v_ids.iter().format(","),
                utr_ids.iter().format(","),
                share[z].cdr3_aa
            );
            let _len = share[z].seq.len();
            let mut lefts = Vec::<Vec<u8>>::new();
            for m in 0..clones.len() {
                let start = clones[m][z].v_start;
                let mut x = clones[m][z].full_seq[0..start].to_vec();
                x.reverse();
                lefts.push(x.to_vec());
            }
            let mut rutrs = Vec::<Vec<u8>>::new();
            for i in 0..utr_ids.len() {
                let mut x = refdata.refs[utr_ids[i]].to_string().as_bytes().to_vec();
                x.reverse();
                rutrs.push(x.to_vec());
            }
            let mut minlen = 1_000_000;
            let mut maxlen = 0;
            for i in 0..lefts.len() {
                minlen = min(minlen, lefts[i].len());
                maxlen = max(maxlen, lefts[i].len());
            }
            for i in 0..rutrs.len() {
                minlen = min(minlen, rutrs[i].len());
                maxlen = max(maxlen, rutrs[i].len());
            }
            let mut dots = Vec::<u8>::new();
            let mut diffs = 0;
            for j in 0..maxlen {
                let mut bases = Vec::<u8>::new();
                for i in 0..lefts.len() {
                    if j < lefts[i].len() {
                        bases.push(lefts[i][j]);
                    }
                }
                for i in 0..rutrs.len() {
                    if j < rutrs[i].len() {
                        bases.push(rutrs[i][j]);
                    }
                }
                let mut diff = false;
                for i in 1..bases.len() {
                    if bases[i] != bases[0] {
                        diff = true;
                    }
                }
                if diff {
                    diffs += 1;
                    dots.push(b'x');
                } else {
                    dots.push(b'.');
                }
            }
            fwriteln!(log, "     {}", strme(&dots));
            for i in 0..rutrs.len() {
                fwriteln!(log, " U = {}", strme(&rutrs[i]));
            }
            for i in 0..lefts.len() {
                if i + 1 <= 9 {
                    fwrite!(log, " ");
                }
                fwriteln!(log, "{} = {}", i + 1, strme(&lefts[i]));
            }
            if !(minlen == maxlen && diffs == 0 && utr_ids.len() == 1) {
                print!("{}", strme(&log));
                *_count += 1;
            }
        }
    }
    if ctl.gen_opt.con_con && clones.len() > 0 {
        // ???????????????????????????????????????
        // NOTE TRUNCATED TO 120 BASES!
        const SHOW: usize = 120;
        for z in 0..clones[0].len() {
            let mut log = Vec::<u8>::new();
            let mut c_ref_ids = Vec::<Option<usize>>::new();
            c_ref_ids.push(share[z].c_ref_id);
            unique_sort(&mut c_ref_ids);
            fwriteln!(
                log,
                "\n[{}] rights for exact subclonotype {}, chain {}, cs = {:?}\n",
                *_count + 1,
                exact_clonotypes.len(),
                z,
                c_ref_ids.iter().format(",")
            );
            let _len = share[z].seq.len();
            let mut rights = Vec::<Vec<u8>>::new();
            let mut bcs = Vec::<String>::new();
            for m in 0..clones.len() {
                let start = clones[m][z].j_stop;
                let mut x = clones[m][z].full_seq[start..].to_vec();
                if x.len() > SHOW {
                    x.truncate(SHOW);
                }
                rights.push(x.to_vec());
                bcs.push(clones[m][0].barcode.clone());
            }
            let mut rconst = Vec::<Vec<u8>>::new();
            for i in 0..c_ref_ids.len() {
                let cid = c_ref_ids[i];
                if cid.is_none() {
                    continue;
                }
                let mut x = refdata.refs[cid.unwrap()].to_string().as_bytes().to_vec();
                if x.len() > SHOW {
                    x.truncate(SHOW);
                }
                /*
                // WARNING!  TO INVESTIGATE, AND NOT NECESSARILY VALID FOR MOUSE!!!!!!!!!!!!!!!
                let n = refdata.name[cid.unwrap()].after("IG");
                if n == "HM" || n == "HA1" || n == "HA2" ||  n == "HG1" || n == "HG2"
                    || n == "HG4"{
                    x.remove(0);
                }
                */
                rconst.push(x.to_vec());
            }
            let mut minlen = 1_000_000;
            let mut maxlen = 0;
            for i in 0..rights.len() {
                minlen = min(minlen, rights[i].len());
                maxlen = max(maxlen, rights[i].len());
            }
            for i in 0..rights.len() {
                minlen = min(minlen, rights[i].len());
                maxlen = max(maxlen, rights[i].len());
            }
            let mut dots = Vec::<u8>::new();
            let mut diffs = 0;
            for j in 0..maxlen {
                let mut bases = Vec::<u8>::new();
                for i in 0..rights.len() {
                    if j < rights[i].len() {
                        bases.push(rights[i][j]);
                    }
                }
                for i in 0..rconst.len() {
                    if j < rconst[i].len() {
                        bases.push(rconst[i][j]);
                    }
                }
                let mut diff = false;
                for i in 1..bases.len() {
                    if bases[i] != bases[0] {
                        diff = true;
                    }
                }
                if diff {
                    diffs += 1;
                    dots.push(b'x');
                } else {
                    dots.push(b'.');
                }
            }
            if dots.len() > SHOW {
                dots.truncate(SHOW);
            }
            fwriteln!(log, "     {}, diffs = {}", strme(&dots), diffs);
            for i in 0..rconst.len() {
                let cid = c_ref_ids[i];
                if cid.is_some() {
                    fwriteln!(
                        log,
                        " C = {}, {}|{}",
                        strme(&rconst[i]),
                        refdata.id[cid.unwrap()],
                        refdata.name[cid.unwrap()]
                    );
                }
            }
            for i in 0..rights.len() {
                if i + 1 <= 9 {
                    fwrite!(log, " ");
                }
                fwriteln!(log, "{} = {} = {}", i + 1, strme(&rights[i]), bcs[i]);
            }
            // if !( minlen == maxlen && diffs == 0 && utr_ids.len() == 1 ) {
            print!("{}", strme(&log));
            *_count += 1;
            // }
        }
    }
}
