// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use amino::*;
use enclone_core::defs::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::*;

pub fn generate_fasta(
    i: usize,
    j: usize,
    oo: usize,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    refdata: &RefData,
    logx: &mut Vec<u8>,
    fout: &mut Box<dyn std::io::Write>,
    faaout: &mut Box<dyn std::io::Write>,
) {
    // Generate FASTA output.

    if !ctl.gen_opt.fasta_filename.is_empty() {
        for (k, u) in exacts[oo].iter().enumerate() {
            for m in 0..rsi[oo].mat.len() {
                if rsi[oo].mat[m][k].is_some() {
                    let r = rsi[oo].mat[m][k].unwrap();
                    let ex = &exact_clonotypes[*u];
                    if ctl.gen_opt.fasta_filename != *"stdout" {
                        fwriteln!(
                            fout,
                            ">group{}.clonotype{}.exact{}.chain{}",
                            i + 1,
                            j + 1,
                            k + 1,
                            m + 1
                        );
                    } else {
                        fwriteln!(
                            logx,
                            ">group{}.clonotype{}.exact{}.chain{}",
                            i + 1,
                            j + 1,
                            k + 1,
                            m + 1
                        );
                    }
                    let mut seq = ex.share[r].seq.clone();
                    let mut cid = ex.share[r].c_ref_id;
                    if cid.is_none() {
                        for l in 0..exacts[oo].len() {
                            if rsi[oo].mat[m][l].is_some() {
                                let r2 = rsi[oo].mat[m][l].unwrap();
                                let ex2 = &exact_clonotypes[exacts[oo][l]];
                                let cid2 = ex2.share[r2].c_ref_id;
                                if cid2.is_some() {
                                    cid = cid2;
                                    break;
                                }
                            }
                        }
                    }
                    if cid.is_some() {
                        let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                        seq.append(&mut cseq);
                        if ctl.gen_opt.fasta_filename != *"stdout" {
                            fwriteln!(fout, "{}", strme(&seq));
                        } else {
                            fwriteln!(logx, "{}", strme(&seq));
                        }
                    }
                }
            }
        }
    }

    // Generate fasta amino acid output.

    if !ctl.gen_opt.fasta_aa_filename.is_empty() {
        for (k, u) in exacts[oo].iter().enumerate() {
            for m in 0..rsi[oo].mat.len() {
                if rsi[oo].mat[m][k].is_some() {
                    let r = rsi[oo].mat[m][k].unwrap();
                    let ex = &exact_clonotypes[*u];
                    if ctl.gen_opt.fasta_aa_filename != *"stdout" {
                        fwriteln!(
                            faaout,
                            ">group{}.clonotype{}.exact{}.chain{}",
                            i + 1,
                            j + 1,
                            k + 1,
                            m + 1
                        );
                    } else {
                        fwriteln!(
                            logx,
                            ">group{}.clonotype{}.exact{}.chain{}",
                            i + 1,
                            j + 1,
                            k + 1,
                            m + 1
                        );
                    }
                    let mut seq = ex.share[r].seq.clone();
                    let mut cid = ex.share[r].c_ref_id;
                    if cid.is_none() {
                        for l in 0..exacts[oo].len() {
                            if rsi[oo].mat[m][l].is_some() {
                                let r2 = rsi[oo].mat[m][l].unwrap();
                                let ex2 = &exact_clonotypes[exacts[oo][l]];
                                let cid2 = ex2.share[r2].c_ref_id;
                                if cid2.is_some() {
                                    cid = cid2;
                                    break;
                                }
                            }
                        }
                    }
                    if cid.is_some() {
                        let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                        seq.append(&mut cseq);
                        if ctl.gen_opt.fasta_aa_filename != *"stdout" {
                            fwriteln!(faaout, "{}", strme(&aa_seq(&seq, 0)));
                        } else {
                            fwriteln!(logx, "{}", strme(&aa_seq(&seq, 0)));
                        }
                    }
                }
            }
        }
    }
}
