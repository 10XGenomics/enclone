// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Storage of enclone visual history, and functions to save and restore.

use crate::packing::*;
use io_utils::*;
use std::fs::File;
use std::io::{BufReader, Read, Write};

#[derive(Default, PartialEq)]
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

const ENCLONE_VISUAL_HISTORY_VERSION: usize = 1;
const HEADER_LENGTH: usize = 40;

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
    // 2. total bytes in file (u32)
    // 3. total bytes up through translated_hist_uniq (u32)
    // 4. data for each of the fields in EncloneVisualHistory,
    //    with translated_input_history and translated_input_hist_uniq stored first.

    pub fn save_as_bytes(&self) -> Vec<u8> {
        let mut bytes = format!(
            "enclone visual history file version{:<4}\n",
            ENCLONE_VISUAL_HISTORY_VERSION
        )
        .as_bytes()
        .to_vec();
        if ENCLONE_VISUAL_HISTORY_VERSION == 1 {
            bytes.append(&mut vec![0 as u8; 8]);
            bytes.append(&mut save_vec_u32(&self.translated_input_history));
            bytes.append(&mut save_vec_string(&self.translated_input_hist_uniq));
            let b = u32_bytes(bytes.len());
            for i in 0..4 {
                bytes[HEADER_LENGTH + 4 + i] = b[i];
            }
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
            let b = u32_bytes(bytes.len());
            for i in 0..4 {
                bytes[HEADER_LENGTH + i] = b[i];
            }
        }
        bytes
    }

    pub fn restore_from_bytes(bytes: &Vec<u8>) -> Result<Self, ()> {
        if bytes.len() < HEADER_LENGTH {
            return Err(());
        }
        let expected_header = format!(
            "enclone visual history file version{:<4}\n",
            ENCLONE_VISUAL_HISTORY_VERSION
        )
        .as_bytes()
        .to_vec();
        if bytes[0..HEADER_LENGTH].to_vec() != expected_header {
            return Err(());
        }
        let mut pos = HEADER_LENGTH + 8;
        let translated_input_history = restore_vec_u32(&bytes, &mut pos)?;
        let translated_input_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let svg_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let summary_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let input1_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let input2_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let displayed_tables_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let table_comp_hist_uniq = restore_vec_vec_u8(&bytes, &mut pos)?;
        let last_widths_hist_uniq = restore_vec_vec_u32(&bytes, &mut pos)?;
        let svg_history = restore_vec_u32(&bytes, &mut pos)?;
        let summary_history = restore_vec_u32(&bytes, &mut pos)?;
        let input1_history = restore_vec_u32(&bytes, &mut pos)?;
        let input2_history = restore_vec_u32(&bytes, &mut pos)?;
        let displayed_tables_history = restore_vec_u32(&bytes, &mut pos)?;
        let table_comp_history = restore_vec_u32(&bytes, &mut pos)?;
        let last_widths_history = restore_vec_u32(&bytes, &mut pos)?;
        let is_blank = restore_vec_bool(&bytes, &mut pos)?;
        let history_index = restore_u32(&bytes, &mut pos)?;
        Ok(EncloneVisualHistory {
            translated_input_history: translated_input_history,
            translated_input_hist_uniq: translated_input_hist_uniq,
            svg_hist_uniq: svg_hist_uniq,
            summary_hist_uniq: summary_hist_uniq,
            input1_hist_uniq: input1_hist_uniq,
            input2_hist_uniq: input2_hist_uniq,
            displayed_tables_hist_uniq: displayed_tables_hist_uniq,
            table_comp_hist_uniq: table_comp_hist_uniq,
            last_widths_hist_uniq: last_widths_hist_uniq,
            svg_history: svg_history,
            summary_history: summary_history,
            input1_history: input1_history,
            input2_history: input2_history,
            displayed_tables_history: displayed_tables_history,
            table_comp_history: table_comp_history,
            last_widths_history: last_widths_history,
            is_blank: is_blank,
            history_index: history_index,
        })
    }

    pub fn save_restore_works(&self) -> bool {
        let bytes = self.save_as_bytes();
        let new = EncloneVisualHistory::restore_from_bytes(&bytes);
        if new.is_err() {
            return false;
        }
        *self == new.unwrap()
    }
}

pub fn read_command_list(filename: &str) -> Result<Vec<String>, ()> {
    let n;
    {
        let mut f = open_for_read![&filename];
        let mut buf = vec![0 as u8; HEADER_LENGTH + 8];
        let res = f.read_exact(&mut buf);
        if res.is_err() {
            return Err(());
        }
        n = u32_from_bytes(&buf[HEADER_LENGTH + 4..HEADER_LENGTH + 8]);
    }
    let mut bytes = vec![0 as u8; n as usize];
    let mut f = open_for_read![&filename];
    let res = f.read_exact(&mut bytes);
    if res.is_err() {
        return Err(());
    }
    let mut pos = HEADER_LENGTH + 8;
    let translated_input_history = restore_vec_u32(&bytes, &mut pos)?;
    let translated_input_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
    let mut commands = Vec::<String>::new();
    for i in translated_input_history.iter() {
        let i = *i as usize;
        if i >= translated_input_history.len() {
            return Err(());
        }
        commands.push(translated_input_hist_uniq[i].clone());
    }
    Ok(commands)
}

pub fn read_enclone_visual_history(filename: &str) -> Result<EncloneVisualHistory, ()> {
    let n;
    {
        let mut f = open_for_read![&filename];
        let mut buf = vec![0 as u8; HEADER_LENGTH + 4];
        let res = f.read_exact(&mut buf);
        if res.is_err() {
            return Err(());
        }
        n = u32_from_bytes(&buf[HEADER_LENGTH..HEADER_LENGTH + 4]);
    }
    let mut bytes = vec![0 as u8; n as usize];
    let mut f = open_for_read![&filename];
    let res = f.read_exact(&mut bytes);
    if res.is_err() {
        return Err(());
    }
    EncloneVisualHistory::restore_from_bytes(&bytes)
}

pub fn write_enclone_visual_history(evh: &EncloneVisualHistory, filename: &str) -> Result<(), ()> {
    let f = File::create(&filename);
    if f.is_err() {
        return Err(());
    }
    let bytes = evh.save_as_bytes();
    let res = f.unwrap().write_all(&bytes);
    if res.is_err() {
        return Err(());
    }
    Ok(())
}
