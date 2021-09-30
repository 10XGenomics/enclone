// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Generate sequential PHYLIP output.  See:
// 1. http://www.atgc-montpellier.fr/phyml/usersguide.php?type=phylip
// 2. http://evolution.genetics.washington.edu/phylip/doc/sequence.html.
// We don't fold lines because it may not be necessary.  See giant value for W;
// code left in place in case it turns out that folding is needed.  This will involve
// a bit more than lowering W.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use io_utils::{fwrite, fwriteln};
use std::fs::File;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};
use string_utils::strme;
use tar::{Builder, Header};

pub fn print_phylip(
    i: usize,
    j: usize,
    oo: usize,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    logx: &mut Vec<u8>,
    phylip_aa: &mut Option<Builder<File>>,
    phylip_dna: &mut Option<Builder<File>>,
) {
    if !ctl.gen_opt.phylip_aa.is_empty() {
        let stdout = ctl.gen_opt.phylip_aa == *"stdout";
        let mut data = Vec::<u8>::new();
        let mut nbases = 0;
        for m in 0..rsi[oo].mat.len() {
            nbases += rsi[oo].seq_del_lens[m];
        }
        if stdout {
            fwriteln!(logx, "");
            fwriteln!(logx, "{} {}", exacts[oo].len(), nbases / 3);
        } else {
            fwriteln!(data, "{} {}", exacts[oo].len(), nbases / 3);
        }
        let mut aa = Vec::<Vec<u8>>::new();
        let mut names = Vec::<String>::new();
        for (k, u) in exacts[oo].iter().enumerate() {
            let ex = &exact_clonotypes[*u];
            let mut seq = Vec::<u8>::new();
            for m in 0..rsi[oo].mat.len() {
                if rsi[oo].mat[m][k].is_none() {
                    seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m] / 3]);
                } else {
                    let r = rsi[oo].mat[m][k].unwrap();
                    let mut z = ex.share[r].aa_mod_indel.clone();
                    seq.append(&mut z);
                }
            }
            names.push(format!("{}", k + 1,));
            aa.append(&mut vec![seq; 1]);
        }
        const W: usize = 10000;
        const PAD: usize = 4;
        let mut name_width = 0;
        for i in 0..names.len() {
            name_width = std::cmp::max(name_width, names[i].len());
        }
        for start in (0..aa[0].len()).step_by(W) {
            if start > 0 {
                if stdout {
                    fwriteln!(logx, "");
                } else {
                    fwriteln!(data, "");
                }
            }
            let stop = std::cmp::min(start + W, aa[0].len());
            for i in 0..aa.len() {
                if stdout {
                    fwrite!(logx, "{}", names[i]);
                    fwrite!(
                        logx,
                        "{}",
                        strme(&vec![b' '; name_width + PAD - names[i].len()])
                    );
                    fwriteln!(logx, "{}", strme(&aa[i][start..stop]));
                } else {
                    fwrite!(data, "{}", names[i]);
                    fwrite!(
                        data,
                        "{}",
                        strme(&vec![b' '; name_width + PAD - names[i].len()])
                    );
                    fwriteln!(data, "{}", strme(&aa[i][start..stop]));
                }
            }
            if stdout {
                fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
            } else {
                fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
            }
        }
        if !stdout {
            let mut header = Header::new_gnu();
            header.set_size(data.len() as u64);
            header.set_cksum();
            header.set_mode(0o0644);
            let now = SystemTime::now();
            header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
            let filename = format!("{}.{}", i + 1, j + 1);
            phylip_aa
                .as_mut()
                .unwrap()
                .append_data(&mut header, &filename, &data[..])
                .unwrap();
        }
    }
    if !ctl.gen_opt.phylip_dna.is_empty() {
        let stdout = ctl.gen_opt.phylip_dna == *"stdout";
        let mut data = Vec::<u8>::new();
        let mut nbases = 0;
        for m in 0..rsi[oo].mat.len() {
            nbases += rsi[oo].seq_del_lens[m];
        }
        if stdout {
            fwriteln!(logx, "");
            fwriteln!(logx, "{} {}", exacts[oo].len(), nbases);
        } else {
            fwriteln!(data, "{} {}", exacts[oo].len(), nbases);
        }
        let mut dna = Vec::<Vec<u8>>::new();
        let mut names = Vec::<String>::new();
        for (k, u) in exacts[oo].iter().enumerate() {
            let ex = &exact_clonotypes[*u];
            let mut seq = Vec::<u8>::new();
            for m in 0..rsi[oo].mat.len() {
                if rsi[oo].mat[m][k].is_none() {
                    seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m]]);
                } else {
                    let r = rsi[oo].mat[m][k].unwrap();
                    let mut s = ex.share[r].seq_del_amino.clone();
                    seq.append(&mut s);
                }
            }
            names.push(format!("{}", k + 1,));
            dna.append(&mut vec![seq; 1]);
        }
        const W: usize = 10000;
        const PAD: usize = 4;
        let mut name_width = 0;
        for i in 0..names.len() {
            name_width = std::cmp::max(name_width, names[i].len());
        }
        for start in (0..dna[0].len()).step_by(W) {
            if start > 0 {
                if stdout {
                    fwriteln!(logx, "");
                } else {
                    fwriteln!(data, "");
                }
            }
            let stop = std::cmp::min(start + W, dna[0].len());
            for i in 0..dna.len() {
                if stdout {
                    fwrite!(logx, "{}", names[i]);
                    fwrite!(
                        logx,
                        "{}",
                        strme(&vec![b' '; name_width + PAD - names[i].len()])
                    );
                    fwriteln!(logx, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                } else {
                    fwrite!(data, "{}", names[i]);
                    fwrite!(
                        data,
                        "{}",
                        strme(&vec![b' '; name_width + PAD - names[i].len()])
                    );
                    fwriteln!(data, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                }
            }
            if stdout {
                fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
            } else {
                fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
            }
        }
        if !stdout {
            let mut header = Header::new_gnu();
            header.set_size(data.len() as u64);
            header.set_cksum();
            header.set_mode(0o0644);
            let now = SystemTime::now();
            header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
            let filename = format!("{}.{}", i + 1, j + 1);
            phylip_dna
                .as_mut()
                .unwrap()
                .append_data(&mut header, &filename, &data[..])
                .unwrap();
        }
    }
}
