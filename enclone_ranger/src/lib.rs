// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod determine_ref;
pub mod main_enclone;
pub mod setup;
pub mod stop;
pub mod vars;

use std::sync::atomic::AtomicBool;

pub static USING_PAGER: AtomicBool = AtomicBool::new(false);
