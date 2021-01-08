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
    pub fn fwr1_nt(&self) -> Option<Vec<u8>> {
        if self.fwr1_start.is_none() || self.cdr1_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.fwr1_start.unwrap() as usize..self.cdr1_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn fwr1_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr1_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr1_aa(&self) -> Option<Vec<u8>> {
        if self.fwr1_start.is_none() || self.cdr1_start.is_none() {
            return None;
        }
        let start = (self.fwr1_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr1_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr1_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr1_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr1_nt(&self) -> Option<Vec<u8>> {
        if self.cdr1_start.is_none() || self.fwr2_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.cdr1_start.unwrap() as usize..self.fwr2_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn cdr1_nt_string(&self) -> Option<String> {
        String::from_utf8(self.cdr1_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr1_aa(&self) -> Option<Vec<u8>> {
        if self.cdr1_start.is_none() || self.fwr2_start.is_none() {
            return None;
        }
        let start = (self.cdr1_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.fwr2_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn cdr1_aa_string(&self) -> Option<String> {
        String::from_utf8(self.cdr1_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr2_nt(&self) -> Option<Vec<u8>> {
        if self.fwr2_start.is_none() || self.cdr2_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.fwr2_start.unwrap() as usize..self.cdr2_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn fwr2_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr2_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr2_aa(&self) -> Option<Vec<u8>> {
        if self.fwr2_start.is_none() || self.cdr2_start.is_none() {
            return None;
        }
        let start = (self.fwr2_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr2_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr2_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr2_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr2_nt(&self) -> Option<Vec<u8>> {
        if self.cdr2_start.is_none() || self.fwr3_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.cdr2_start.unwrap() as usize..self.fwr3_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn cdr2_nt_string(&self) -> Option<String> {
        String::from_utf8(self.cdr2_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr2_aa(&self) -> Option<Vec<u8>> {
        if self.cdr2_start.is_none() || self.fwr3_start.is_none() {
            return None;
        }
        let start = (self.cdr2_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.fwr3_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn cdr2_aa_string(&self) -> Option<String> {
        String::from_utf8(self.cdr2_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr3_nt(&self) -> Option<Vec<u8>> {
        if self.fwr3_start.is_none() {
            return None;
        }
        Some(self.nt_sequence[self.fwr3_start.unwrap() as usize..self.cdr3_start as usize].to_vec())
    }
    pub fn fwr3_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr3_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr3_aa(&self) -> Option<Vec<u8>> {
        if self.fwr3_start.is_none() {
            return None;
        }
        let start = (self.fwr3_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr3_start - self.v_start) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr3_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr3_aa().as_deref().unwrap().to_vec()).ok()
    }
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
    pub fn fwr4_nt(&self) -> Option<Vec<u8>> {
        if self.fwr4_end.is_none() {
            return None;
        }
        Some(self.nt_sequence[self.cdr3_end as usize..self.fwr4_end.unwrap() as usize].to_vec())
    }
    pub fn fwr4_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr4_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr4_aa(&self) -> Option<Vec<u8>> {
        if self.fwr4_end.is_none() {
            return None;
        }
        let start = (self.cdr3_end as usize - self.v_start as usize) / 3;
        let end = (self.fwr4_end.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr4_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr4_aa().as_deref().unwrap().to_vec()).ok()
    }
}

impl ExactSubClonotypeChain {
    pub fn fwr1_nt(&self) -> Option<Vec<u8>> {
        if self.fwr1_start.is_none() || self.cdr1_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.fwr1_start.unwrap() as usize..self.cdr1_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn fwr1_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr1_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr1_aa(&self) -> Option<Vec<u8>> {
        if self.fwr1_start.is_none() || self.cdr1_start.is_none() {
            return None;
        }
        let start = (self.fwr1_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr1_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr1_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr1_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr1_nt(&self) -> Option<Vec<u8>> {
        if self.cdr1_start.is_none() || self.fwr2_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.cdr1_start.unwrap() as usize..self.fwr2_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn cdr1_nt_string(&self) -> Option<String> {
        String::from_utf8(self.cdr1_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr1_aa(&self) -> Option<Vec<u8>> {
        if self.cdr1_start.is_none() || self.fwr2_start.is_none() {
            return None;
        }
        let start = (self.cdr1_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.fwr2_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn cdr1_aa_string(&self) -> Option<String> {
        String::from_utf8(self.cdr1_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr2_nt(&self) -> Option<Vec<u8>> {
        if self.fwr2_start.is_none() || self.cdr2_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.fwr2_start.unwrap() as usize..self.cdr2_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn fwr2_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr2_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr2_aa(&self) -> Option<Vec<u8>> {
        if self.fwr2_start.is_none() || self.cdr2_start.is_none() {
            return None;
        }
        let start = (self.fwr2_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr2_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr2_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr2_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr2_nt(&self) -> Option<Vec<u8>> {
        if self.cdr2_start.is_none() || self.fwr3_start.is_none() {
            return None;
        }
        Some(
            self.nt_sequence[self.cdr2_start.unwrap() as usize..self.fwr3_start.unwrap() as usize]
                .to_vec(),
        )
    }
    pub fn cdr2_nt_string(&self) -> Option<String> {
        String::from_utf8(self.cdr2_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn cdr2_aa(&self) -> Option<Vec<u8>> {
        if self.cdr2_start.is_none() || self.fwr3_start.is_none() {
            return None;
        }
        let start = (self.cdr2_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.fwr3_start.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn cdr2_aa_string(&self) -> Option<String> {
        String::from_utf8(self.cdr2_aa().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr3_nt(&self) -> Option<Vec<u8>> {
        if self.fwr3_start.is_none() {
            return None;
        }
        Some(self.nt_sequence[self.fwr3_start.unwrap() as usize..self.cdr3_start as usize].to_vec())
    }
    pub fn fwr3_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr3_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr3_aa(&self) -> Option<Vec<u8>> {
        if self.fwr3_start.is_none() {
            return None;
        }
        let start = (self.fwr3_start.unwrap() as usize - self.v_start as usize) / 3;
        let end = (self.cdr3_start - self.v_start) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr3_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr3_aa().as_deref().unwrap().to_vec()).ok()
    }
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
    pub fn fwr4_nt(&self) -> Option<Vec<u8>> {
        if self.fwr4_end.is_none() {
            return None;
        }
        Some(self.nt_sequence[self.cdr3_end as usize..self.fwr4_end.unwrap() as usize].to_vec())
    }
    pub fn fwr4_nt_string(&self) -> Option<String> {
        String::from_utf8(self.fwr4_nt().as_deref().unwrap().to_vec()).ok()
    }
    pub fn fwr4_aa(&self) -> Option<Vec<u8>> {
        if self.fwr4_end.is_none() {
            return None;
        }
        let start = (self.cdr3_end as usize - self.v_start as usize) / 3;
        let end = (self.fwr4_end.unwrap() as usize - self.v_start as usize) / 3;
        Some(self.aa_sequence[start as usize..end as usize].to_vec())
    }
    pub fn fwr4_aa_string(&self) -> Option<String> {
        String::from_utf8(self.fwr4_aa().as_deref().unwrap().to_vec()).ok()
    }
}
