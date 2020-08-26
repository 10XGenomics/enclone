// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::print_utils1::*;
use amino::*;
use ansi_escape::*;
use enclone_core::defs::*;
use enclone_core::print_tools::*;
use enclone_proto::types::*;
use io_utils::*;
use itertools::Itertools;
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define reference segment identifiers, one per chain.  For the V segment, there is
// also an optional donor reference sequence alignment.  (For now.  We might do the
// same thing later for J segments.)
// ◼ 1. We vote using the number of cells, whereas a better way would be to test all
// ◼    the alternatives to find the best match.
// ◼ 2. Defining a constant region identifier for the entire clonotype is
// ◼    biologically dubious.
// ◼ 3. Maybe we only need to do this for pass 2.

pub fn define_column_info(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    mat: &Vec<Vec<Option<usize>>>,
    refdata: &RefData,
) -> ColInfo {
    let cols = mat.len();

    // Define cvars.

    let mut cvars = Vec::<Vec<String>>::new();
    for cx in 0..cols {
        let mut have_notes = false;
        for u in 0..exacts.len() {
            let ex = &exact_clonotypes[exacts[u]];
            let m = mat[cx][u];
            if m.is_some() {
                let m = m.unwrap();
                let ex = &ex.share[m];
                if ex.vs_notesx.len() > 0 {
                    have_notes = true;
                }
            }
        }
        let mut cv = Vec::<String>::new();
        for i in 0..ctl.clono_print_opt.cvars.len() {
            if ctl.clono_print_opt.cvars[i] == "notes" && !have_notes {
                continue;
            }
            cv.push(ctl.clono_print_opt.cvars[i].to_string());
        }
        cvars.push(cv);
    }

    // Compute CDR3 starts, etc.

    let mut fr1_starts = Vec::<usize>::new();
    let mut fr2_starts = Vec::<usize>::new();
    let mut fr3_starts = Vec::<usize>::new();
    let mut cdr1_starts = Vec::<usize>::new();
    let mut cdr2_starts = Vec::<usize>::new();
    let mut cdr3_starts = Vec::<usize>::new();
    let mut cdr3_lens = Vec::<usize>::new();
    let mut seq_lens = Vec::<usize>::new();
    let mut seq_del_lens = Vec::<usize>::new();
    for cx in 0..cols {
        for u in 0..exacts.len() {
            let ex = &exact_clonotypes[exacts[u]];
            let m = mat[cx][u];
            if m.is_some() {
                let m = m.unwrap();
                let exm = &ex.share[m];
                cdr3_lens.push(exm.cdr3_aa.len());
                seq_lens.push(exm.seq.len());
                seq_del_lens.push(exm.seq_del.len());

                // The logic below with testing i < start while incrementing start seems fishy.

                let mut start = exm.cdr1_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                cdr1_starts.push(start);
                let mut start = exm.cdr2_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                cdr2_starts.push(start);
                let mut start = exm.fr1_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                fr1_starts.push(start);
                let mut start = exm.fr2_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                fr2_starts.push(start);
                let mut start = exm.fr3_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                fr3_starts.push(start);
                let mut start = exm.cdr3_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                cdr3_starts.push(start);
                break;
            }
        }
    }

    // Compute reference info.

    let mut uids = vec![None; cols];
    let mut vids = vec![0; cols];
    let mut vpids = vec![None; cols];
    let mut vpids_d = vec![None; cols];
    let mut vpids_a = vec![None; cols];
    let mut dids = vec![None; cols];
    let mut jids = vec![0; cols];
    let mut cids = vec![None; cols];
    for col in 0..cols {
        let mut u = Vec::<usize>::new();
        let mut v = Vec::<usize>::new();
        let mut vp = Vec::<(usize, Option<usize>, Option<usize>, Option<usize>)>::new();
        let mut d = Vec::<usize>::new();
        let mut j = Vec::<usize>::new();
        let mut c = Vec::<usize>::new();
        for e in 0..exacts.len() {
            let clonotype_id = exacts[e];
            let ex = &exact_clonotypes[clonotype_id];
            let m = mat[col][e];
            if m.is_some() {
                let x = &ex.share[m.unwrap()];
                if x.u_ref_id.is_some() {
                    for _ in 0..ex.ncells() {
                        u.push(x.u_ref_id.unwrap());
                    }
                }
                for _ in 0..ex.ncells() {
                    v.push(x.v_ref_id);
                    vp.push((
                        x.v_ref_id,
                        x.v_ref_id_donor,
                        x.v_ref_id_donor_donor,
                        x.v_ref_id_donor_alt_id,
                    ));
                }
                if x.d_ref_id.is_some() {
                    for _ in 0..ex.ncells() {
                        d.push(x.d_ref_id.unwrap());
                    }
                }
                for _ in 0..ex.ncells() {
                    j.push(x.j_ref_id);
                }
                if x.c_ref_id.is_some() {
                    for _ in 0..ex.ncells() {
                        c.push(x.c_ref_id.unwrap());
                    }
                }
            }
        }
        u.sort();
        v.sort();
        vp.sort();
        d.sort();
        j.sort();
        c.sort();
        let mut uf = Vec::<(u32, usize)>::new();
        make_freq(&u, &mut uf);
        if !uf.is_empty() {
            uids[col] = Some(uf[0].1);
        }
        let mut vf = Vec::<(u32, usize)>::new();
        make_freq(&v, &mut vf);
        vids[col] = vf[0].1;
        let mut to_delete = vec![false; vp.len()];
        for i in 0..vp.len() {
            if vp[i].0 != vids[col] {
                to_delete[i] = true;
            }
        }
        erase_if(&mut vp, &to_delete);
        let mut vpf = Vec::<(u32, (usize, Option<usize>, Option<usize>, Option<usize>))>::new();
        make_freq(&vp, &mut vpf);
        vpids[col] = (vpf[0].1).1;
        vpids_d[col] = (vpf[0].1).2;
        vpids_a[col] = (vpf[0].1).3;
        let mut df = Vec::<(u32, usize)>::new();
        make_freq(&d, &mut df);
        if !df.is_empty() {
            dids[col] = Some(df[0].1);
        }
        let mut jf = Vec::<(u32, usize)>::new();
        make_freq(&j, &mut jf);
        jids[col] = jf[0].1;
        let mut cf = Vec::<(u32, usize)>::new();
        make_freq(&c, &mut cf);
        if !cf.is_empty() {
            cids[col] = Some(cf[0].1);
        }
    }

    // Compute seqss and seqss_amino.

    let mut seqss = Vec::<Vec<Vec<u8>>>::new();
    let mut seqss_amino = Vec::<Vec<Vec<u8>>>::new();
    let nexacts = exacts.len();
    for cx in 0..cols {
        let mut seqs = Vec::<Vec<u8>>::new();
        let mut seqs_amino = Vec::<Vec<u8>>::new();
        for u in 0..nexacts {
            let m = mat[cx][u];
            if m.is_some() {
                let m = m.unwrap();
                seqs.push(exact_clonotypes[exacts[u]].share[m].seq_del.clone());
                seqs_amino.push(exact_clonotypes[exacts[u]].share[m].seq_del_amino.clone());
            } else {
                seqs.push(Vec::<u8>::new());
                seqs_amino.push(Vec::<u8>::new());
            }
        }
        seqss.push(seqs.clone());
        seqss_amino.push(seqs_amino.clone());
    }

    // Show segment names.  We used ◼ as a separator character, but that does not render well
    // as a fixed-width character in Google Docs.  So we changed it to ◆.

    let mut chain_descrip = vec![String::new(); cols];
    for cx in 0..cols {
        let vid = vids[cx];
        let mut vdescrip = format!("{}", refdata.id[vid]);
        if vpids[cx].is_some() {
            vdescrip = format!(
                "{}.{}.{}",
                vdescrip,
                vpids_d[cx].unwrap() + 1,
                vpids_a[cx].unwrap() + 1
            );
        }
        chain_descrip[cx] = format!("{}|{}", vdescrip, refdata.name[vid]);
        let did = dids[cx];
        if did.is_some() {
            let did = did.unwrap();
            chain_descrip[cx] += &format!(" ◆ {}|{}", refdata.id[did], refdata.name[did]);
        }
        let jid = jids[cx];
        chain_descrip[cx] += &format!(" ◆ {}|{}", refdata.id[jid], refdata.name[jid]);
    }

    // Return.

    ColInfo {
        uids: uids,
        vids: vids,
        vpids: vpids,
        dids: dids,
        jids: jids,
        cids: cids,
        fr1_starts: fr1_starts,
        fr2_starts: fr2_starts,
        fr3_starts: fr3_starts,
        cdr1_starts: cdr1_starts,
        cdr2_starts: cdr2_starts,
        cdr3_starts: cdr3_starts,
        cdr3_lens: cdr3_lens,
        seq_lens: seq_lens,
        seq_del_lens: seq_del_lens,
        seqss: seqss,
        seqss_amino: seqss_amino,
        chain_descrip: chain_descrip,
        mat: Vec::<Vec<Option<usize>>>::new(),
        cvars: cvars,
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Add header text to mlog.

pub fn add_header_text(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    rord: &Vec<usize>,
    mat: &Vec<Vec<Option<usize>>>,
    mut mlog: &mut Vec<u8>,
) {
    let nexacts = exacts.len();
    let cols = mat.len();
    for cx in 0..cols {
        let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
        for u in 0..nexacts {
            let m = mat[cx][u];
            if m.is_some() {
                let m = m.unwrap();
                vref = exact_clonotypes[exacts[u]].share[m].vs.to_ascii_vec();
                jref = exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
            }
        }
        let mut seqs = Vec::<Vec<u8>>::new();
        let mut full_seqs = Vec::<Vec<u8>>::new();
        for u in 0..nexacts {
            let ex = &exact_clonotypes[exacts[rord[u]]];
            let m = mat[cx][rord[u]];
            if m.is_some() {
                let m = m.unwrap();
                seqs.push(ex.share[m].seq_del.clone());
                full_seqs.push(ex.share[m].full_seq.clone());
            } else {
                seqs.push(Vec::<u8>::new());
                full_seqs.push(Vec::<u8>::new());
            }
        }
        let mut simple = false;
        if ctl.clono_print_opt.note_simple && vref.len() + jref.len() >= seqs[0].len() {
            let n = seqs[0].len() - jref.len();
            let mut vj = vref[0..n].to_vec();
            vj.append(&mut jref.clone());
            if vj == seqs[0] {
                simple = true;
            }
        }
        if ctl.clono_print_opt.seqc || ctl.clono_print_opt.full_seqc || simple {
            fwriteln!(&mut mlog, "CHAIN {}", cx + 1);
        }
        if ctl.clono_print_opt.seqc {
            fwriteln!(&mut mlog, "• {}", strme(&seqs[0]));
        }
        if ctl.clono_print_opt.full_seqc {
            fwriteln!(&mut mlog, "• {}", strme(&full_seqs[0]));
        }
        if simple {
            fwriteln!(&mut mlog, "• This chain is simple.");
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Insert universal and donor reference rows.
// Possibly buggy WRT reference indels that we insert.

pub fn insert_reference_rows(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    show_aa: &Vec<Vec<usize>>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    row1: &Vec<String>,
    drows: &mut Vec<Vec<String>>,
    rows: &mut Vec<Vec<String>>,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
) {
    let cols = rsi.seq_del_lens.len();
    if drows.len() >= 1 {
        for pass in 1..=2 {
            let mut row = Vec::<String>::new();
            if pass == 1 {
                row.push("reference".to_string());
            } else {
                row.push("donor ref".to_string());
            }
            for _ in 1..row1.len() {
                row.push("\\ext".to_string());
            }
            for cz in 0..cols {
                let mut refseq = Vec::<u8>::new();
                let mut vlen: usize;
                let vseq: Vec<u8>;
                if pass == 1 {
                    vlen = refdata.refs[rsi.vids[cz]].len();
                    vseq = refdata.refs[rsi.vids[cz]].to_ascii_vec();
                } else if rsi.vpids[cz].is_none() {
                    vlen = refdata.refs[rsi.vids[cz]].len();
                    vseq = refdata.refs[rsi.vids[cz]].to_ascii_vec();
                } else {
                    vlen = dref[rsi.vpids[cz].unwrap()].nt_sequence.len();
                    vseq = dref[rsi.vpids[cz].unwrap()].nt_sequence.clone();
                }
                let mut trim = ctl.heur.ref_v_trim;
                vlen -= trim;
                let mut jlen = refdata.refs[rsi.jids[cz]].len() - trim;
                let jseq = refdata.refs[rsi.jids[cz]].to_ascii_vec();
                let mut gap = rsi.seq_del_lens[cz] as isize - vlen as isize - jlen as isize;

                if gap < -2 * (trim as isize) {
                    let mut bcs = Vec::<String>::new();
                    for u in 0..exacts.len() {
                        let ex = &exact_clonotypes[exacts[u]];
                        for i in 0..ex.clones.len() {
                            bcs.push(ex.clones[i][0].barcode.clone());
                        }
                    }
                    bcs.sort();
                    eprintln!("\ncz = {}", cz);
                    eprintln!("pass = {}", pass);
                    eprintln!("seq_del.len() = {}", rsi.seq_del_lens[cz]);
                    eprintln!("vlen = {}", vlen);
                    eprintln!("jlen = {}", jlen);
                    eprintln!("gap = seq_del.len() - vlen - jlen");
                    panic!(
                        "Something is wrong because gap is {}, which is negative.\n\
                        This is happening for the clonotype with these barcodes:\n{}.",
                        gap,
                        bcs.iter().format(",")
                    );
                }

                if gap < 0 {
                    let mut ptrim = (-gap) / 2;
                    if (-gap) % 2 == 1 {
                        ptrim += 1;
                    }
                    vlen += ptrim as usize;
                    jlen += ptrim as usize;
                    gap += 2 * ptrim;
                    trim -= ptrim as usize;
                }

                let gap = gap as usize;
                for j in 0..vlen {
                    refseq.push(vseq[j]);
                }
                for _ in 0..gap {
                    refseq.push(b'-');
                }
                for j in 0..jlen {
                    refseq.push(jseq[j + trim]);
                }
                let mut refx = String::new();
                let n1 = (rsi.fr2_starts[cz] - rsi.cdr1_starts[cz]) / 3;
                let n2 = (rsi.fr3_starts[cz] - rsi.cdr2_starts[cz]) / 3;
                let cs3 = rsi.cdr3_starts[cz] / 3;
                let n3 = rsi.cdr3_lens[cz];
                let amino = &ctl.clono_print_opt.amino;
                let show1 = amino.contains(&"cdr1".to_string())
                    && rsi.cdr1_starts[cz] <= rsi.fr2_starts[cz];
                let show2 = amino.contains(&"cdr2".to_string())
                    && rsi.cdr2_starts[cz] <= rsi.fr3_starts[cz];
                for k in 0..show_aa[cz].len() {
                    let p = show_aa[cz][k];
                    if show1 && k > 0 && p == rsi.cdr1_starts[cz] / 3 {
                        refx += " ";
                    }
                    if show2
                        && k > 0
                        && p == rsi.cdr2_starts[cz] / 3
                        && (!show1 || show_aa[cz][k - 1] != rsi.cdr1_starts[cz] / 3 + n1 - 1)
                    {
                        refx += " ";
                    }
                    if k > 0 && p == cs3 {
                        refx += " ";
                    }
                    if 3 * p + 3 > refseq.len() || refseq[3 * p..3 * p + 3].contains(&b'-') {
                        refx += "◦";
                    } else {
                        let mut log = Vec::<u8>::new();
                        let aa = codon_to_aa(&refseq[3 * p..3 * p + 3]);
                        if ctl.gen_opt.color == "codon".to_string() {
                            emit_codon_color_escape(&refseq[3 * p..3 * p + 3], &mut log);
                            log.push(aa);
                            emit_end_escape(&mut log);
                        } else {
                            color_by_property(&vec![aa], &mut log);
                        }
                        refx += strme(&log);
                    }

                    if k < show_aa[cz].len() - 1 {
                        if show1 && p == rsi.cdr1_starts[cz] / 3 + n1 - 1 {
                            refx += " ";
                        }
                        if show2 && p == rsi.cdr2_starts[cz] / 3 + n2 - 1 {
                            refx += " ";
                        }
                        if p == cs3 + n3 - 1 {
                            refx += " ";
                        }
                    }
                }
                row.push(refx);
                for _ in 1..rsi.cvars[cz].len() {
                    row.push("".to_string());
                }
            }
            rows.push(row);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn build_table_stuff(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    rsi: &ColInfo,
    vars: &Vec<Vec<usize>>,
    show_aa: &Vec<Vec<usize>>,
    row1: &mut Vec<String>,
    justify: &mut Vec<u8>,
    drows: &mut Vec<Vec<String>>,
    rows: &mut Vec<Vec<String>>,
    lvars: &Vec<String>,
) {
    // Build lead header row and justification to match.

    let cols = rsi.vids.len();
    let nexacts = exacts.len();
    if !ctl.clono_print_opt.bu {
        row1.push("#".to_string());
    } else {
        row1.push("#  barcode".to_string());
    }
    justify.push(b'l');
    for i in 0..lvars.len() {
        let mut x = lvars[i].to_string();
        if x.contains(':') {
            x = x.before(":").to_string();
        }
        row1.push(x.clone());
        justify.push(justification(&x));
    }

    // Insert main chain row.  Then insert chain info row if we're using CHAIN_SPLIT.

    let mut row = vec!["".to_string(); row1.len()];
    for j in 0..cols {
        if rsi.chain_descrip[j].contains(&"IGH".to_string())
            || rsi.chain_descrip[j].contains(&"TRB".to_string())
        {
            row.push(bold(&format!("CHAIN {}", j + 1)));
        } else {
            row.push(format!("CHAIN {}", j + 1));
        }
        for _ in 1..rsi.cvars[j].len() {
            row.push("\\ext".to_string());
        }
    }
    rows.push(row);
    let mut row = vec!["".to_string(); row1.len()];
    for j in 0..cols {
        if rsi.chain_descrip[j].contains(&"IGH".to_string())
            || rsi.chain_descrip[j].contains(&"TRB".to_string())
        {
            row.push(bold(&format!("{}", rsi.chain_descrip[j])));
        } else {
            row.push(format!("{}", rsi.chain_descrip[j]));
        }
        for _ in 1..rsi.cvars[j].len() {
            row.push("\\ext".to_string());
        }
    }
    rows.push(row);

    // Insert divider row (horizontal line across the chains).

    let mut row = vec!["".to_string(); lvars.len() + 1];
    let mut ncall = 0;
    for j in 0..cols {
        ncall += rsi.cvars[j].len();
    }
    row.append(&mut vec!["\\hline".to_string(); ncall]);
    rows.push(row);

    // Insert position rows.

    *drows = insert_position_rows(&ctl, &rsi, &show_aa, &vars, &row1);
    let mut drows2 = drows.clone();
    rows.append(&mut drows2);

    // Insert main per-chain header row.

    let mut row = vec!["".to_string(); row1.len()];
    for cx in 0..cols {
        let show = &show_aa[cx];
        for j in 0..rsi.cvars[cx].len() {
            if rsi.cvars[cx][j] != "amino".to_string() {
                if drows.is_empty() {
                    row.push(rsi.cvars[cx][j].to_string());
                } else {
                    row.push("".to_string());
                }
            } else {
                for u in 0..nexacts {
                    let m = rsi.mat[cx][u];
                    if m.is_some() {
                        let m = m.unwrap();
                        let mut s = String::new();
                        let amino = &ctl.clono_print_opt.amino;
                        let printing_cdr1 = amino.contains(&"cdr1".to_string())
                            && rsi.cdr1_starts[cx] <= rsi.fr2_starts[cx];
                        let printing_cdr2 = amino.contains(&"cdr2".to_string())
                            && rsi.cdr2_starts[cx] <= rsi.fr3_starts[cx];
                        let printing_cdr3 = amino.contains(&"cdr3".to_string());
                        if printing_cdr1 {
                            let cs1 = rsi.cdr1_starts[cx] / 3;
                            let lead = show.iter().position(|x| *x == cs1).unwrap();
                            s += &" ".repeat(lead);
                            if lead > 0 {
                                s += " ";
                            }
                            let mut n = 0;
                            let x = &exact_clonotypes[exacts[u]].share[m];
                            if x.cdr1_start <= x.fr2_start {
                                n = (x.fr2_start - x.cdr1_start) / 3;
                            }
                            if n >= 4 {
                                let left = (n - 3) / 2;
                                let right = n - left - 4;
                                s += &"═".repeat(left);
                                s += "CDR1";
                                s += &"═".repeat(right);
                            } else if n == 3 {
                                s += "CR1";
                            } else if n == 2 {
                                s += "C1";
                            } else if n == 1 {
                                s += "1";
                            }
                        }
                        if printing_cdr2 {
                            let lead;
                            let cs2 = rsi.cdr2_starts[cx] / 3;
                            let pos2 = show.iter().position(|x| *x == cs2).unwrap();
                            if !printing_cdr1 {
                                lead = pos2;
                            } else {
                                let cs1 = rsi.cdr1_starts[cx] / 3;
                                let pos1 = show.iter().position(|x| *x == cs1).unwrap();
                                lead = pos2 - pos1 - (rsi.fr2_starts[cx] - rsi.cdr1_starts[cx]) / 3;
                            }
                            s += &" ".repeat(lead);
                            if lead > 0 {
                                s += " ";
                            }
                            if printing_cdr1 {
                                s += " ";
                            }
                            let mut n = 0;
                            if rsi.cdr2_starts[cx] <= rsi.fr3_starts[cx] {
                                n = (rsi.fr3_starts[cx] - rsi.cdr2_starts[cx]) / 3;
                            }
                            if n >= 4 {
                                let left = (n - 3) / 2;
                                let right = n - left - 4;
                                s += &"═".repeat(left);
                                s += "CDR2";
                                s += &"═".repeat(right);
                            } else if n == 3 {
                                s += "CR2";
                            } else if n == 2 {
                                s += "C2";
                            } else if n == 1 {
                                s += "2";
                            }
                        }
                        if printing_cdr3 {
                            let lead;
                            let cs3 = rsi.cdr3_starts[cx] / 3;
                            let pos3 = show.iter().position(|x| *x == cs3).unwrap();
                            if !printing_cdr1 && !printing_cdr2 {
                                lead = pos3;
                            } else if printing_cdr2 {
                                let cs2 = rsi.cdr2_starts[cx] / 3;
                                let pos2 = show.iter().position(|x| *x == cs2).unwrap();
                                lead = 1 + pos3
                                    - pos2
                                    - (rsi.fr3_starts[cx] - rsi.cdr2_starts[cx]) / 3;
                            } else {
                                let cs1 = rsi.cdr1_starts[cx] / 3;
                                let pos1 = show.iter().position(|x| *x == cs1).unwrap();
                                lead = 1 + pos3
                                    - pos1
                                    - (rsi.fr2_starts[cx] - rsi.cdr1_starts[cx]) / 3;
                            }
                            s += &" ".repeat(lead);
                            if lead > 0 {
                                s += " ";
                            }
                            let n = exact_clonotypes[exacts[u]].share[m].cdr3_aa.len();
                            let left = (n - 3) / 2;
                            let right = n - left - 4;
                            s += &"═".repeat(left);
                            s += "CDR3";
                            s += &"═".repeat(right);
                        }
                        row.push(s);
                        break;
                    }
                }
            }
        }
    }
    rows.push(row);
}
