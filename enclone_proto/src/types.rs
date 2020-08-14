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

impl ClonotypeChain {
    pub fn cdr3_nt(&self) -> &[u8] {
        &self.nt_sequence[self.cdr3_start as usize..self.cdr3_end as usize]
    }
    pub fn cdr3_nt_string(&self) -> String {
        std::str::from_utf8(self.cdr3_nt()).unwrap().to_string()
    }
    pub fn cdr3_aa(&self) -> &[u8] {
        let start = (self.cdr3_start - self.v_start) / 3;
        let end = (self.cdr3_end - self.v_start) / 3;
        &self.aa_sequence[start as usize..end as usize]
    }
    pub fn cdr3_aa_string(&self) -> String {
        std::str::from_utf8(self.cdr3_aa()).unwrap().to_string()
    }
}

impl ExactSubClonotypeChain {
    pub fn cdr3_nt(&self) -> &[u8] {
        &self.nt_sequence[self.cdr3_start as usize..self.cdr3_end as usize]
    }
    pub fn cdr3_nt_string(&self) -> String {
        std::str::from_utf8(self.cdr3_nt()).unwrap().to_string()
    }
    pub fn cdr3_aa(&self) -> &[u8] {
        let start = (self.cdr3_start - self.v_start) / 3;
        let end = (self.cdr3_end - self.v_start) / 3;
        &self.aa_sequence[start as usize..end as usize]
    }
    pub fn cdr3_aa_string(&self) -> String {
        std::str::from_utf8(self.cdr3_aa()).unwrap().to_string()
    }
}
