// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use io_utils::*;
use std::process::Command;
use vector_utils::*;

pub fn test_sec_mem(ctl: &mut EncloneControl) -> Result<(), String> {
    let is_bcr = !ctl.gen_opt.tcr;

    // Test for okness of sec/mem args.

    let mut vars = ctl.parseable_opt.pcols.clone();
    vars.append(&mut ctl.clono_print_opt.lvars.clone());
    unique_sort(&mut vars);
    ctl.gen_opt.using_secmem =
        bin_member(&vars, &"sec".to_string()) || bin_member(&vars, &"mem".to_string());
    if !ctl.gen_opt.using_secmem
        && ctl.parseable_opt.pout.len() > 0
        && ctl.parseable_opt.pcols.len() == 0
    {
        if ctl.gen_opt.species == "human" || ctl.gen_opt.species == "mouse" {
            if is_bcr {
                let mut have_bam = true;
                for g in ctl.origin_info.gex_path.iter() {
                    if g.len() == 0 {
                        have_bam = false;
                        break;
                    }
                    let bam = format!("{}/possorted_genome_bam.bam", g);
                    if !path_exists(&bam) {
                        have_bam = false;
                        break;
                    }
                }
                if have_bam {
                    let o = Command::new("samtools")
                        .arg("--help")
                        .output()
                        .expect("failed to execute samtools");
                    let status = o.status.code().unwrap();
                    if status == 0 {
                        ctl.gen_opt.using_secmem = true;
                    }
                }
            }
        }
    }
    if bin_member(&vars, &"sec".to_string()) || bin_member(&vars, &"mem".to_string()) {
        if ctl.gen_opt.species != "human" && ctl.gen_opt.species != "mouse" {
            return Err(format!(
                "\nThe lvars sec and mem can only be used for data from human and mouse.\n"
            ));
        }
        if !is_bcr {
            return Err(format!(
                "\nThe lvars sec and mem do not make sense for TCR data.\n"
            ));
        }
        for g in ctl.origin_info.gex_path.iter() {
            if g.len() == 0 {
                return Err(format!(
                    "\nThe lvars sec and mem can only be used if GEX data are provided.\n"
                ));
            }
            let bam = format!("{}/possorted_genome_bam.bam", g);
            if !path_exists(&bam) {
                return Err(format!(
                    "\nThe lvars sec and mem can only be used if the file\n\
                    pos_sorted_genome_bam.bam is provided.  We did not see it at this path\n\
                    {}.",
                    g
                ));
            }
        }
        let o = Command::new("samtools")
            .arg("--help")
            .output()
            .expect("failed to execute samtools");
        let status = o.status.code().unwrap();
        if status != 0 {
            return Err(format!(
                "\nThe lvars sec and mem can only be used if the samtools\n\
                executable is in your path.\n"
            ));
        }
    }
    Ok(())
}
