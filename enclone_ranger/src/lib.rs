// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod blacklist;
pub mod determine_ref;
pub mod flag_defective;
pub mod inconsistent;
pub mod main_enclone;
pub mod opt_d_val;
pub mod sec_mem;
pub mod setup;
pub mod stop;
pub mod subset;
pub mod vars;

use std::sync::atomic::AtomicBool;

pub static USING_PAGER: AtomicBool = AtomicBool::new(false);
