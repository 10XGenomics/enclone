// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::USING_PAGER;
use enclone_args::proc_args::proc_args;
use enclone_core::defs::EncloneControl;
use std::sync::atomic::Ordering::SeqCst;
use vector_utils::erase_if;

pub fn setup(
    mut ctl: &mut EncloneControl,
    args: &Vec<String>,
    argsx: &mut Vec<String>,
) -> Result<(), String> {
    ctl.gen_opt.nopager = true;
    {
        let mut args = args.clone();
        let mut to_delete = vec![false; args.len()];
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                to_delete[i] = true;
            } else if args[i] == "FORCE_EXTERNAL" {
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PRE=") {
                to_delete[i] = true;
            }
        }
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
    }
    ctl.pretty = true;
    ctl.gen_opt.h5 = true;
    USING_PAGER.store(false, SeqCst);

    // Process args (and set defaults for them).

    proc_args(&mut ctl, args)?;
    Ok(())
}
