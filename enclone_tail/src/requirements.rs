// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test requirements.

use enclone_core::defs::{EncloneControl, ExactClonotype};

pub fn test_requirements(
    pics: &Vec<String>,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    nclono2: usize,
    two_chain: usize,
    three_chain: usize,
) -> Result<(), String> {
    // Test for required number of false positives.

    let mut msg = String::new();
    if ctl.gen_opt.required_fps.is_some() {
        let mut fps = 0;
        for i in 0..pics.len() {
            if pics[i].contains("WARNING:") {
                fps += 1;
            }
        }
        if fps != ctl.gen_opt.required_fps.unwrap() {
            msg += &mut format!(
                "\nA \"false positive\" is a clonotype that contains cells from multiple\n\
                 donors.  You invoked enclone with the argument REQUIRED_FPS={}, but we found\n\
                 {} false positives, so the requirement is not met.\n",
                ctl.gen_opt.required_fps.unwrap(),
                fps
            );
        }
    }

    // Test for required number of cells.

    if ctl.gen_opt.required_cells.is_some() {
        let mut ncells = 0;
        for i in 0..pics.len() {
            for x in exacts[i].iter() {
                ncells += exact_clonotypes[*x].ncells();
            }
        }
        if ctl.gen_opt.required_cells.unwrap() != ncells {
            msg += &mut format!(
                "\nThe required number of cells is {}, but you actually have {}.\n",
                ctl.gen_opt.required_cells.unwrap(),
                ncells,
            );
        }
    }

    // Test for required number of clonotypes.

    let nclono = exacts.len();
    if ctl.gen_opt.required_clonotypes.is_some()
        && ctl.gen_opt.required_clonotypes.unwrap() != nclono
    {
        msg += &mut format!(
            "\nThe required number of clonotypes is {}, but you actually have {}.\n",
            ctl.gen_opt.required_clonotypes.unwrap(),
            nclono,
        );
    }

    // Test for required number of donors

    if ctl.gen_opt.required_donors.is_some() {
        let ndonors = ctl.origin_info.donors;
        if ctl.gen_opt.required_donors.unwrap() != ndonors {
            msg += &format!(
                "\nThe required number of donors is {}, but you actually have {}.\n",
                ctl.gen_opt.required_donors.unwrap(),
                ndonors,
            );
        }
    }

    // Test for required number of >=2 cell clonotypes.

    if ctl.gen_opt.required_two_cell_clonotypes.is_some()
        && ctl.gen_opt.required_two_cell_clonotypes.unwrap() != nclono2
    {
        msg += &mut format!(
            "\nThe required number of two-cell clonotypes is {}, but you actually have {}.\n",
            ctl.gen_opt.required_two_cell_clonotypes.unwrap(),
            nclono2,
        );
    }

    // Test for required number of two chain clonotypes.

    if ctl.gen_opt.required_two_chain_clonotypes.is_some()
        && ctl.gen_opt.required_two_chain_clonotypes.unwrap() != two_chain
    {
        msg += &mut format!(
            "\nThe required number of two-chain clonotypes is {}, but you actually have {}.\n",
            ctl.gen_opt.required_two_chain_clonotypes.unwrap(),
            two_chain,
        );
    }

    // Test for required number of three chain clonotypes.

    if ctl.gen_opt.required_three_chain_clonotypes.is_some()
        && ctl.gen_opt.required_three_chain_clonotypes.unwrap() != three_chain
    {
        msg += &mut format!(
            "\nThe required number of three-chain clonotypes is {}, but you actually have {}.\n",
            ctl.gen_opt.required_three_chain_clonotypes.unwrap(),
            three_chain,
        );
    }

    // Test for required number of datasets

    if ctl.gen_opt.required_datasets.is_some() {
        let ndatasets = ctl.origin_info.n();
        if ctl.gen_opt.required_datasets.unwrap() != ndatasets {
            msg += &mut format!(
                "\nThe required number of datasets is {}, but you actually have {}.\n",
                ctl.gen_opt.required_datasets.unwrap(),
                ndatasets,
            );
        }
    }
    if !msg.is_empty() {
        return Err(msg);
    }
    Ok(())
}
