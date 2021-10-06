// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod main_enclone;
pub mod stop;

use std::sync::atomic::AtomicBool;

pub static USING_PAGER: AtomicBool = AtomicBool::new(false);
