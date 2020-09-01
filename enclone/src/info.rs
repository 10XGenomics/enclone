// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file provides the single function build_info.

use vdj_ann::*;

use self::refx::*;
use crate::read_json::*;
use amino::*;
use ansi_escape::*;
use debruijn::{dna_string::*, Mer};
use enclone_core::defs::*;
use enclone_core::print_tools::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::AtomicBool;
use string_utils::*;
use vector_utils::*;

pub fn build_info(
    refdata: &RefData,
    ctl: &EncloneControl,
    exact_clonotypes: &mut Vec<ExactClonotype>,
    fate: &mut Vec<HashMap<String, String>>,
) -> Vec<CloneInfo> {
    // Build info about clonotypes.  We create a data structure info.
    // An entry in info is a clonotype having appropriate properties.
    //
    // Much of the information in a CloneInfo object is redundant.  So we could probably
    // improve both time and space computational performance by reducing that redundancy.

    let exiting = AtomicBool::new(false);
    let mut info = Vec::<CloneInfo>::new();
    let mut results = Vec::<(
        usize,
        Vec<CloneInfo>,
        ExactClonotype,
        Vec<(usize, String, String)>,
    )>::new();
    for i in 0..exact_clonotypes.len() {
        results.push((i, Vec::new(), exact_clonotypes[i].clone(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut lens = Vec::<usize>::new();
        let mut tigs = Vec::<Vec<u8>>::new();
        let mut tigs_amino = Vec::<Vec<u8>>::new();
        let mut aa_mod_indel = Vec::<Vec<u8>>::new();
        let mut tigs_ins = Vec::<Vec<(usize, Vec<u8>)>>::new();
        let mut tigsp = Vec::<DnaString>::new();
        let mut has_del = Vec::<bool>::new();
        let mut orig_tigs = Vec::<DnaString>::new();
        let (mut vs, mut js) = (Vec::<DnaString>::new(), Vec::<DnaString>::new());
        let mut vsids = Vec::<usize>::new();
        let mut jsids = Vec::<usize>::new();
        let mut cdr3s = Vec::<String>::new();
        let mut cdr3_aa = Vec::<String>::new();
        let mut chain_types = Vec::<String>::new();
        let mut vs_notes = Vec::<String>::new();
        let mut vs_notesx = Vec::<String>::new();
        let p = &mut res.2;
        for j in 0..p.share.len() {
            let x = &mut p.share[j];
            tigsp.push(DnaString::from_acgt_bytes(&x.seq));
            // INCORRECT, TO DO SOMETHING ABOUT LATER:
            orig_tigs.push(DnaString::from_acgt_bytes(&x.full_seq));
            let jid = x.j_ref_id;
            js.push(refdata.refs[jid].clone());

            // If there is a deletion in a V segment region, edit the contig sequence,
            // inserting hyphens where the deletion is, and if there is an insertion, delete it.

            let vid = x.v_ref_id;
            let jid = x.j_ref_id;
            let mut annv = x.annv.clone();
            vsids.push(vid);
            jsids.push(jid);
            let mut vsnx = String::new();
            // DELETION
            if annv.len() == 2 && annv[1].0 == annv[0].0 + annv[0].1 {
                let mut t = Vec::<u8>::new();
                let (mut del_start, mut del_stop) = (annv[0].1, annv[1].3);
                for i in 0..del_start {
                    t.push(x.seq[i as usize]);
                }
                for _ in del_start..del_stop {
                    t.push(b'-');
                }
                for i in (annv[1].0 as usize)..x.seq.len() {
                    t.push(x.seq[i as usize]);
                }
                lens.push(t.len());
                tigs.push(t.clone());
                if del_start % 3 != 0 {
                    // Bad solution here, should pick optimal choice.
                    let offset = del_start % 3 - 3;
                    del_start -= offset;
                    del_stop -= offset;
                    t.clear();
                    for i in 0..del_start {
                        t.push(x.seq[i as usize]);
                    }
                    for _ in del_start..del_stop {
                        t.push(b'-');
                    }
                    for i in ((annv[1].0 - offset) as usize)..x.seq.len() {
                        t.push(x.seq[i as usize]);
                    }
                }
                annv[0].1 += (del_stop - del_start) + annv[1].1;
                annv.truncate(1);
                tigs_amino.push(t.clone());
                let mut aa = Vec::<u8>::new();
                for p in (0..=t.len() - 3).step_by(3) {
                    if t[p] == b'-' {
                        aa.push(b'-');
                    } else {
                        aa.push(codon_to_aa(&t[p..p + 3]));
                    }
                }

                aa_mod_indel.push(aa);
                tigs_ins.push(Vec::new());
                has_del.push(true);
            // INSERTION
            } else if annv.len() == 2 && annv[1].3 == annv[0].3 + annv[0].1 {
                let ins_len = (annv[1].0 - annv[0].0 - annv[0].1) as usize;
                let ins_pos = (annv[0].0 + annv[0].1) as usize;
                let mut t = Vec::<u8>::new();
                let mut nt = Vec::<u8>::new();
                for i in 0..x.seq.len() {
                    if i < ins_pos || i >= ins_pos + ins_len {
                        t.push(x.seq[i]);
                    } else {
                        nt.push(x.seq[i]);
                    }
                }
                has_del.push(true); // DOES NOT MAKE SENSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                lens.push(t.len());
                tigs.push(t.clone());
                let ins = vec![(ins_pos, nt)];
                tigs_ins.push(ins);
                tigs_amino.push(t);

                // Optimize to compute entry in aa_mod_indel and the inserted aa sequence.

                let aa_full = aa_seq(&x.seq, 0);
                let ref_aa = aa_seq(&refdata.refs[vid].to_ascii_vec(), 0);
                let ins_len_aa = ins_len / 3;
                const EXT: usize = 10;
                let ins_pos_low;
                if ins_pos / 3 < EXT {
                    ins_pos_low = 0;
                } else {
                    ins_pos_low = ins_pos / 3 - EXT;
                }
                let mut ins_pos_high =
                    std::cmp::min(ins_pos / 3 + EXT, aa_full.len() - ins_len_aa + 1);
                ins_pos_high = std::cmp::min(ins_pos_high, ref_aa.len() - ins_len_aa + 1);
                let mut mis = Vec::<(usize, usize, Vec<u8>)>::new();
                for j in ins_pos_low..ins_pos_high {
                    let mut y = Vec::<u8>::new();
                    for k in 0..aa_full.len() {
                        if k < j || k >= j + ins_len_aa {
                            y.push(aa_full[k].clone());
                        }
                    }
                    let mut m = 0;
                    for l in 0..ref_aa.len() {
                        if l < y.len() && ref_aa[l] != y[l] {
                            m += 1;
                        }
                    }
                    mis.push((m, j, y.clone()));
                }
                mis.sort();
                aa_mod_indel.push(mis[0].2.clone());
                let ins_aa_pos = mis[0].1.clone();
                let mut aax = Vec::<u8>::new();
                let b = 3 * ins_aa_pos;
                for p in 0..ins_len_aa {
                    emit_codon_color_escape(&x.seq[b + 3 * p..b + 3 * p + 3], &mut aax);
                    let aa = codon_to_aa(&x.seq[b + 3 * p..b + 3 * p + 3]);
                    aax.push(aa);
                    emit_end_escape(&mut aax);
                }
                vsnx = format!("ins = {} at {}", strme(&aax), ins_aa_pos);

                // Finish up ann.

                annv[0].1 += annv[1].1;
                annv.truncate(1);
            } else {
                has_del.push(false);
                lens.push(x.seq.len());
                tigs.push(x.seq.clone());
                tigs_amino.push(x.seq.clone());
                aa_mod_indel.push(aa_seq(&x.seq, 0));
                tigs_ins.push(Vec::new());
            }

            // Save reference V segment.  However in the case where there is a
            // single indel between the contig and the reference sequence, edit the
            // reference V segment accordingly.

            let rt = &refdata.refs[vid as usize];
            if x.annv.len() == 2 {
                if x.annv[0].1 as usize > rt.len() {
                    let msg = format!("x.annv[0].1 = {}, rt.len() = {}", x.annv[0].1, rt.len());
                    json_error(None, &ctl, &exiting, &msg);
                }
                let mut r = rt.slice(0, x.annv[0].1 as usize).to_owned();
                // deletion
                if x.annv[1].0 == x.annv[0].0 + x.annv[0].1 {
                    // DEAD CODE
                    for m in x.annv[1].3 as usize..rt.len() {
                        r.push(rt.get(m));
                    }
                    vs.push(r.clone());
                    vs_notes.push(format!(
                        "has deletion of {} bases relative to reference",
                        x.annv[1].3 - x.annv[0].1
                    ));
                    vs_notesx.push("".to_string());
                // insertion
                } else if x.annv[1].3 == x.annv[0].3 + x.annv[0].1 {
                    /*
                    for m in x.annv[0].0 + x.annv[0].1..x.annv[1].0 {
                        if m as usize >= x.seq.len() {
                            eprintln!( "\nannotation problem with {}", strme(&x.seq) );
                        }
                        r.push( *x.seq.get(m as usize).unwrap() );
                    }
                    */
                    for m in x.annv[1].3 as usize..rt.len() {
                        r.push(rt.get(m));
                    }
                    vs.push(r.clone());
                    vs_notes.push("".to_string());
                } else {
                    // maybe can't happen
                    vs.push(rt.clone());
                    // At one point there was a bug in which the following line was missing.
                    // This caused a traceback on "enclone 123085 RE".  It is interesting because
                    // the traceback did not get back to the main program, even with
                    // "enclone 123085 RE NOPRETTY".
                    vs_notes.push("".to_string());
                    vsnx = "".to_string();
                }
            } else {
                vs.push(rt.clone());
                vs_notes.push(String::new());
                vsnx = "".to_string();
            }
            cdr3s.push(x.cdr3_dna.clone());
            cdr3_aa.push(x.cdr3_aa.clone());
            chain_types.push(x.chain_type.clone());

            // Add to notes if there's a J/C delta. This likely represents an error.

            let z = &p.clones[0][j];
            if z.c_start.is_some() {
                let delta = z.c_start.unwrap() as isize - z.j_stop as isize;
                if delta != 0 {
                    if vsnx.len() > 0 {
                        vsnx += "; ";
                    }
                    if delta > 0 {
                        vsnx += &mut format!("gap from J stop to C start = {}", delta);
                    } else {
                        if delta != -1 || ctl.gen_opt.jc1 {
                            vsnx += &mut format!("J and C segs overlap by {}", -delta);
                        }
                    }
                }
            }
            vs_notesx.push(vsnx);

            // Modify the exact subclonotype to fill in some members.
            // This is the only place where build_info modifies the exact subclonotype.

            x.seq_del = tigs[tigs.len() - 1].clone();
            x.seq_del_amino = tigs_amino[tigs_amino.len() - 1].clone();
            x.aa_mod_indel = aa_mod_indel[aa_mod_indel.len() - 1].clone();
            x.ins = tigs_ins[tigs_ins.len() - 1].clone();
            x.vs = vs[vs.len() - 1].clone();
            x.vs_notesx = vs_notesx[vs_notesx.len() - 1].clone();
            x.js = js[js.len() - 1].clone();
        }
        let mut origin = Vec::<usize>::new();
        for j in 0..exact_clonotypes[i].clones.len() {
            origin.push(exact_clonotypes[i].clones[j][0].dataset_index);
        }
        unique_sort(&mut origin);
        let shares = &exact_clonotypes[i].share;
        let mut placed = false;
        for i1 in 0..shares.len() {
            if shares[i1].left {
                for i2 in 0..shares.len() {
                    if !shares[i2].left {
                        placed = true;
                        let lensx = [lens[i1], lens[i2]].to_vec();
                        let tigsx = [tigs[i1].clone(), tigs[i2].clone()].to_vec();
                        let tigs_aminox = [tigs_amino[i1].clone(), tigs_amino[i2].clone()].to_vec();
                        let tigspx = [tigsp[i1].clone(), tigsp[i2].clone()].to_vec();
                        let has_delx = [has_del[i1], has_del[i2]].to_vec();
                        let orig_tigsx = [orig_tigs[i1].clone(), orig_tigs[i2].clone()].to_vec();
                        let vsx = [vs[i1].clone(), vs[i2].clone()].to_vec();
                        let jsx = [js[i1].clone(), js[i2].clone()].to_vec();
                        let cdr3sx = [cdr3s[i1].clone(), cdr3s[i2].clone()].to_vec();
                        let cdr3_aax = [cdr3_aa[i1].clone(), cdr3_aa[i2].clone()].to_vec();
                        let chain_typesx =
                            [chain_types[i1].clone(), chain_types[i2].clone()].to_vec();
                        let vsidsx = [vsids[i1], vsids[i2]].to_vec();
                        let jsidsx = [jsids[i1], jsids[i2]].to_vec();
                        let exact_cols = vec![i1, i2];
                        res.1.push(CloneInfo {
                            lens: lensx,
                            tigs: tigsx,
                            tigs_amino: tigs_aminox,
                            tigsp: tigspx,
                            has_del: has_delx,
                            orig_tigs: orig_tigsx,
                            clonotype_id: i,
                            exact_cols: exact_cols,
                            clonotype_index: i, // CLEARLY UNNEEDED
                            origin: origin.clone(),
                            vs: vsx.clone(),
                            dref: vec![None; vsx.len()],
                            vs_notesx: vs_notesx.clone(),
                            js: jsx,
                            vsids: vsidsx,
                            jsids: jsidsx,
                            cdr3s: cdr3sx,
                            cdr3_aa: cdr3_aax,
                            chain_types: chain_typesx,
                        });
                    }
                }
            }
        }

        // Incorporate improper cells if they are onesies.  Note that we're dropping the
        // improper cells having two or more chains.

        if !placed && shares.len() > 1 && !ctl.merge_all_impropers {
            let ex = &exact_clonotypes[i];
            for j in 0..ex.clones.len() {
                res.3.push((
                    ex.clones[j][0].dataset_index,
                    ex.clones[j][0].barcode.clone(),
                    "failed IMPROPER filter".to_string(),
                ));
            }
        }
        if !placed && (shares.len() == 1 || ctl.merge_all_impropers) {
            let mut exact_cols = Vec::<usize>::new();
            for i in 0..tigs.len() {
                exact_cols.push(i);
            }
            res.1.push(CloneInfo {
                lens: lens,
                tigs: tigs,
                tigs_amino: tigs_amino,
                tigsp: tigsp,
                has_del: has_del,
                orig_tigs: orig_tigs,
                clonotype_id: i,
                exact_cols: exact_cols,
                clonotype_index: i, // CLEARLY UNNEEDED
                origin: origin.clone(),
                vs: vs.clone(),
                dref: vec![None; vs.len()],
                vs_notesx: vs_notesx,
                js: js,
                vsids: vsids,
                jsids: jsids,
                cdr3s: cdr3s,
                cdr3_aa: cdr3_aa,
                chain_types: chain_types,
            });
        }
    });

    // Cumulate info.  This is single threaded and could probably be speeded up.

    for i in 0..results.len() {
        info.append(&mut results[i].1);
        exact_clonotypes[i] = results[i].2.clone();
        for j in 0..results[i].3.len() {
            fate[results[i].3[j].0].insert(results[i].3[j].1.clone(), results[i].3[j].2.clone());
        }
    }

    // Sort info.

    info.par_sort();

    // Done.

    info
}
