// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Generate clustal output.  See:
// 1. http://meme-suite.org/doc/clustalw-format.html
// 2. https://www.ebi.ac.uk/seqdb/confluence/display/THD/Help+-+Clustal+Omega+FAQ
//    at "What do the consensus symbols mean in the alignment?".

use enclone_core::defs::*;
use io_utils::*;
use std::fs::File;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};
use string_utils::*;
use tar::{Builder, Header};
use vector_utils::*;

pub fn print_clustal(
    i: usize,
    j: usize,
    oo: usize,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    logx: &mut Vec<u8>,
    clustal_aa: &mut Option<Builder<File>>,
    clustal_dna: &mut Option<Builder<File>>,
) {
    if ctl.gen_opt.clustal_aa.len() > 0 {
        let stdout = ctl.gen_opt.clustal_aa == "stdout".to_string();
        let mut data = Vec::<u8>::new();
        if stdout {
            fwriteln!(logx, "");
            fwriteln!(logx, "CLUSTALW\n");
        } else {
            fwriteln!(data, "CLUSTALW\n");
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
            names.push(format!("{}.{}.{}", i + 1, j + 1, k + 1,));
            aa.append(&mut vec![seq; 1]);
        }
        const W: usize = 60;
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
                    fwriteln!(logx, "{}  {}", strme(&aa[i][start..stop]), stop - start);
                } else {
                    fwrite!(data, "{}", names[i]);
                    fwrite!(
                        data,
                        "{}",
                        strme(&vec![b' '; name_width + PAD - names[i].len()])
                    );
                    fwriteln!(data, "{}  {}", strme(&aa[i][start..stop]), stop - start);
                }
            }
            if stdout {
                fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
            } else {
                fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
            }
            for p in start..stop {
                let mut res = Vec::<u8>::new();
                for i in 0..aa.len() {
                    res.push(aa[i][p]);
                }
                unique_sort(&mut res);
                if res.solo() {
                    if stdout {
                        fwrite!(logx, "*");
                    } else {
                        fwrite!(data, "*");
                    }
                } else {
                    let mut con = false;
                    'pass: for pass in 1..=2 {
                        let x: Vec<&[u8]>;
                        // Conservative mutations
                        if pass == 1 {
                            x = vec![
                                b"STA", b"NEQK", b"NHQK", b"NDEQ", b"QHRK", b"MILV", b"MILF",
                                b"HY", b"FYW",
                            ];
                        } else {
                            // Semi-conservative mutations
                            x = vec![
                                b"CSA", b"ATV", b"SAG", b"STNK", b"STPA", b"SGND", b"SNDEQK",
                                b"NDEQHK", b"NEQHRK", b"FVLIM", b"HFY",
                            ];
                        }
                        for y in x.iter() {
                            let mut sub = true;
                            for c in res.iter() {
                                if !y.contains(c) {
                                    sub = false;
                                    break;
                                }
                            }
                            if sub {
                                let sym;
                                if pass == 1 {
                                    sym = ":";
                                } else {
                                    sym = ".";
                                }
                                if stdout {
                                    fwrite!(logx, "{}", sym);
                                } else {
                                    fwrite!(data, "{}", sym);
                                }
                                con = true;
                                break 'pass;
                            }
                        }
                    }
                    if !con {
                        if stdout {
                            fwrite!(logx, " ");
                        } else {
                            fwrite!(data, " ");
                        }
                    }
                }
            }
            if stdout {
                fwriteln!(logx, "");
            } else {
                fwriteln!(data, "");
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
            clustal_aa
                .as_mut()
                .unwrap()
                .append_data(&mut header, &filename, &data[..])
                .unwrap();
        }
    }
    if ctl.gen_opt.clustal_dna.len() > 0 {
        let stdout = ctl.gen_opt.clustal_dna == "stdout".to_string();
        let mut data = Vec::<u8>::new();
        if stdout {
            fwriteln!(logx, "");
            fwriteln!(logx, "CLUSTALW\n");
        } else {
            fwriteln!(data, "CLUSTALW\n");
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
            names.push(format!("{}.{}.{}", i + 1, j + 1, k + 1,));
            dna.append(&mut vec![seq; 1]);
        }
        const W: usize = 60;
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
            for p in start..stop {
                let mut res = Vec::<u8>::new();
                for i in 0..dna.len() {
                    res.push(dna[i][p]);
                }
                unique_sort(&mut res);
                if res.solo() {
                    if stdout {
                        fwrite!(logx, "*");
                    } else {
                        fwrite!(data, "*");
                    }
                } else {
                    if stdout {
                        fwrite!(logx, " ");
                    } else {
                        fwrite!(data, " ");
                    }
                }
            }
            if stdout {
                fwriteln!(logx, "");
            } else {
                fwriteln!(data, "");
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
            clustal_dna
                .as_mut()
                .unwrap()
                .append_data(&mut header, &filename, &data[..])
                .unwrap();
        }
    }
}
