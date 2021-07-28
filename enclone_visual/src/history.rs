// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Storage of enclone visual history, and functions to save and restore.  These have not been
// optimized.

use crate::*;
use serde::{Deserialize, Serialize};

#[derive(Default, Deserialize, Serialize)]
pub struct EncloneVisualHistory {
    //
    // more or less uniqued history:
    //
    pub svg_hist_uniq: Vec<String>,     // each entry is an SVG
    pub summary_hist_uniq: Vec<String>, // each entry is a summary
    pub input1_hist_uniq: Vec<String>,  // each entry is the originating command 1
    pub input2_hist_uniq: Vec<String>,  // each entry is the originating command 2
    pub translated_input_hist_uniq: Vec<String>, // each entry is the translated originating command
    pub displayed_tables_hist_uniq: Vec<String>, // each entry is the tables that are displayed
    pub table_comp_hist_uniq: Vec<Vec<u8>>, // each entry is the compressed list of all tables
    pub last_widths_hist_uniq: Vec<Vec<u32>>,
    //
    // parallel vectors, with one entry for each command entered in the text box:
    //
    pub svg_history: Vec<u32>,              // each entry is an SVG
    pub summary_history: Vec<u32>,          // each entry is a summary
    pub input1_history: Vec<u32>,           // each entry is the originating command 1
    pub input2_history: Vec<u32>,           // each entry is the originating command 2
    pub translated_input_history: Vec<u32>, // each entry is the translated originating command
    pub displayed_tables_history: Vec<u32>, // each entry is the tables that are displayed
    pub table_comp_history: Vec<u32>,       // each entry is the compressed list of all tables
    pub last_widths_history: Vec<u32>,      // last widths for clonotype pictures
    pub is_blank: Vec<bool>,                // set if the SVG is empty
    //
    // index of "current" position in those vectors, plus one:
    //
    pub history_index: u32,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn u32_bytes(x: usize) -> Vec<u8> {
    (x as u32).to_le_bytes().to_vec()
}

pub fn u32_from_bytes(x: &[u8]) -> u32 {
    u32::from_le_bytes([x[0], x[1], x[2], x[3]])
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_string(x: &Vec<String>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        bytes.append(&mut x[i].as_bytes().to_vec());
    }
    bytes
}

pub fn restore_vec_string(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<String>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![String::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        let s = String::from_utf8(x[*pos..*pos + k].to_vec());
        if s.is_err() {
            return Err(());
        }
        *pos += k;
        y[j] = s.unwrap();
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u8(x: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        bytes.append(&mut x[i].clone());
    }
    bytes
}

pub fn restore_vec_vec_u8(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<Vec<u8>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u8>::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        y[j] = x[*pos..*pos + k].to_vec();
        *pos += k;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u32(x: &Vec<Vec<u32>>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut u32_bytes(x[i].len()));
        for j in 0..x[i].len() {
            bytes.append(&mut x[i][j].to_le_bytes().to_vec());
        }
    }
    bytes
}

pub fn restore_vec_vec_u32(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<Vec<u32>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u32>::new(); n];
    for j in 0..n {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + 4 * k > x.len() {
            return Err(());
        }
        for _ in 0..k {
            y[j].push(u32_from_bytes(&x[*pos..*pos + 4]));
            *pos += 4;
        }
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_bool(x: &Vec<bool>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.push(if x[i] { 1 } else { 0 });
    }
    bytes
}

pub fn restore_vec_bool(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<bool>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + n > x.len() {
        return Err(());
    }
    let mut y = vec![false; n];
    for j in 0..n {
        y[j] = if *pos == 1 { true } else { false };
        *pos += 1;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_u32(x: &Vec<u32>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.append(&mut u32_bytes(x.len()));
    for i in 0..x.len() {
        bytes.append(&mut x[i].to_le_bytes().to_vec());
    }
    bytes
}

pub fn restore_vec_u32(x: &Vec<u8>, pos: &mut usize) -> Result<Vec<u32>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + 4 * n > x.len() {
        return Err(());
    }
    let mut y = vec![0; n];
    for j in 0..n {
        y[j] = u32_from_bytes(&x[*pos..*pos + 4]);
        *pos += 4;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_u32(x: u32) -> Vec<u8> {
    x.to_le_bytes().to_vec()
}

pub fn restore_u32(x: &Vec<u8>, pos: &mut usize) -> Result<u32, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let y = u32_from_bytes(&x[*pos..*pos + 4]);
    *pos += 4;
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

const ENCLONE_VISUAL_HISTORY_VERSION: usize = 1;

impl EncloneVisualHistory {
    //
    // requirements for saving and restoring history
    //
    // 1. File structure must not change, except when we explicitly change it.
    // 2. Must be able to read old files.
    // 3. Must be versioned.
    // 4. Must have text header.
    // 5. Must store reasonably efficiently.
    // 6. Must be able to extract command list without reading entire file.
    // 7. File structure must be guaranteed the same across environments (so, e.g., no usize).
    //
    // structure description
    //
    // 1. text header = 40 bytes = "enclone visual history file version ***\n",
    //    where *** is a positive integer, padded on the right with blanks
    // 2. data for each of the fields in EncloneVisualHistory,
    //    with translated_input_history and translated_input_hist_uniq stored first.

    pub fn save_as_bytes(&self) -> Vec<u8> {
        let mut bytes = format!(
            "enclone visual history file version{:<4}\n",
            ENCLONE_VISUAL_HISTORY_VERSION
        )
        .as_bytes()
        .to_vec();
        if ENCLONE_VISUAL_HISTORY_VERSION == 1 {
            bytes.append(&mut save_vec_u32(&self.translated_input_history));
            bytes.append(&mut save_vec_string(&self.translated_input_hist_uniq));
            bytes.append(&mut save_vec_string(&self.svg_hist_uniq));
            bytes.append(&mut save_vec_string(&self.summary_hist_uniq));
            bytes.append(&mut save_vec_string(&self.input1_hist_uniq));
            bytes.append(&mut save_vec_string(&self.input2_hist_uniq));
            bytes.append(&mut save_vec_string(&self.displayed_tables_hist_uniq));
            bytes.append(&mut save_vec_vec_u8(&self.table_comp_hist_uniq));
            bytes.append(&mut save_vec_vec_u32(&self.last_widths_hist_uniq));
            bytes.append(&mut save_vec_u32(&self.svg_history));
            bytes.append(&mut save_vec_u32(&self.summary_history));
            bytes.append(&mut save_vec_u32(&self.input1_history));
            bytes.append(&mut save_vec_u32(&self.input2_history));
            bytes.append(&mut save_vec_u32(&self.displayed_tables_history));
            bytes.append(&mut save_vec_u32(&self.table_comp_history));
            bytes.append(&mut save_vec_u32(&self.last_widths_history));
            bytes.append(&mut save_vec_bool(&self.is_blank));
            bytes.append(&mut save_u32(self.history_index));
        }
        bytes
    }

    pub fn restore_from_bytes(bytes: &[u8]) -> Self {
        serde_json::from_str(&strme(&bytes)).unwrap()
    }
}
