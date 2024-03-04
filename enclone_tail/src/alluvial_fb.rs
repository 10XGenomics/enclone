// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Make alluvial tables for feature barcode data.  We determine cellular using vdj_cells,
// which is not the only way of doing it.

use enclone_core::defs::EncloneControl;
use enclone_core::stringulate::*;
use itertools::Itertools;
use tables::print_tabular_vbox;

pub fn description_table(ctl: &EncloneControl, logx: &mut Vec<u8>) {
    if ctl.gen_opt.vis_dump {
        let mut need = false;
        let n = ctl.origin_info.n();
        for i in 0..n {
            if ctl.origin_info.descrips[i].len() > 0
                && ctl.origin_info.descrips[i] != ctl.origin_info.dataset_id[i]
            {
                need = true;
            }
        }
        if need {
            let mut rows = Vec::<Vec<String>>::new();
            let mut csv_rows = Vec::<Vec<String>>::new();
            rows.push(vec!["id".to_string(), "description".to_string()]);
            for i in 0..n {
                rows.push(vec![
                    ctl.origin_info.dataset_id[i].clone(),
                    ctl.origin_info.descrips[i].clone(),
                ]);
                csv_rows.push(vec![
                    ctl.origin_info.dataset_id[i].clone(),
                    ctl.origin_info.descrips[i].clone(),
                ]);
            }
            let mut display_text = String::new();
            print_tabular_vbox(&mut display_text, &rows, 0, &b"l|l".to_vec(), false, false);
            let mut spreadsheet_text = String::new();
            for r in csv_rows.iter() {
                spreadsheet_text += &mut format!("{}\n", r.iter().format("\t "));
            }
            let f = DescriptionTable {
                display_text: display_text.clone(),
                spreadsheet_text: spreadsheet_text.clone(),
            };
            logx.append(&mut f.to_string().as_bytes().to_vec());
        }
    }
}
