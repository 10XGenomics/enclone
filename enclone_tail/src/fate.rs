// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::barcode_fate::BarcodeFate;
use enclone_core::defs::EncloneControl;
use enclone_core::enclone_structs::BarcodeFates;
use io_utils::*;

use std::io::Write;
use tables::*;
use vector_utils::*;

pub fn print_fate(ctl: &EncloneControl, fate: &[BarcodeFates], logx: &mut Vec<u8>) {
    // Print barcode fate.

    fwriteln!(logx, "2. barcode fate");
    let mut fates = Vec::<String>::new();
    for i in 0..fate.len() {
        for f in fate[i].iter() {
            match f.1 {
                BarcodeFate::WeakChains => {
                    if !ctl.clono_filt_opt_def.weak_chains {
                        continue;
                    }
                }
                BarcodeFate::FoursieKill => {
                    if !ctl.clono_filt_opt_def.weak_foursies {
                        continue;
                    }
                }
                BarcodeFate::NotGexCell => {
                    if ctl.clono_filt_opt_def.ngex {
                        continue;
                    }
                }
                BarcodeFate::Qual => {
                    if !ctl.clono_filt_opt.qual_filter {
                        continue;
                    }
                }
                BarcodeFate::Umi => {
                    if !ctl.cr_opt.umi_filt {
                        continue;
                    }
                }
                BarcodeFate::UmiRatio => {
                    if !ctl.cr_opt.umi_ratio_filt {
                        continue;
                    }
                }
                BarcodeFate::GelBeadContamination => {
                    if ctl.gen_opt.nwhitef {
                        continue;
                    }
                }
                BarcodeFate::DuplicatedBarcode => {
                    if !ctl.clono_filt_opt_def.bc_dup {
                        continue;
                    }
                }
                BarcodeFate::Cross => {
                    if ctl.clono_filt_opt_def.ncross {
                        continue;
                    }
                }
                BarcodeFate::Improper => {
                    if ctl.merge_all_impropers {
                        continue;
                    }
                }
                BarcodeFate::GraphFilter => {
                    if ctl.cr_opt.ngraph_filter {
                        continue;
                    }
                }
                BarcodeFate::Doublet => {
                    if !ctl.clono_filt_opt_def.doublet {
                        continue;
                    }
                }
                BarcodeFate::Signature => {
                    if !ctl.clono_filt_opt_def.signature {
                        continue;
                    }
                }
                BarcodeFate::NotAsmCell => {
                    if ctl.gen_opt.ncell {
                        continue;
                    }
                }
                BarcodeFate::NonProductive => {}
            }

            fates.push(format!("failed {} filter", f.1.label()));
        }
    }
    fates.sort();
    let mut freq = Vec::<(u32, String)>::new();
    make_freq(&fates, &mut freq);
    let mut rows = Vec::<Vec<String>>::new();
    rows.push(vec!["barcodes".to_string(), "why deleted".to_string()]);
    rows.push(vec!["\\hline".to_string(); 2]);
    for i in 0..freq.len() {
        rows.push(vec![format!("{}", freq[i].0), freq[i].1.clone()]);
    }
    rows.push(vec![format!("{}", fates.len()), "total".to_string()]);
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, b"r|l".as_ref(), false, false);
    log.truncate(log.len() - 1);
    log = log.replace('\n', "\n   ");
    fwrite!(logx, "   {}\n", log);
}
