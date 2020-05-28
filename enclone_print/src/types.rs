// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

//!
//! This crate defines the data structure that would represent the clonotypes
//! computed by enclone.
//!

include!(concat!(env!("OUT_DIR"), "/enclone.types.rs"));

impl From<&bio::alignment::Alignment> for Alignment {
    fn from(al: &bio::alignment::Alignment) -> Self {
        Alignment {
            ref_start: al.ystart as u32,
            cigar: al.cigar(false),
        }
    }
}

// TODO: Donor names?
// TODO: Aggr metadata structure
