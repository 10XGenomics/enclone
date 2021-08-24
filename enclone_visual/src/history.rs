// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Storage of enclone visual history, and functions to save and restore.

use crate::packing::*;
use io_utils::*;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Error, ErrorKind, Read, Seek, SeekFrom, Write};

#[derive(Default, PartialEq, Clone)]
pub struct EncloneVisualHistory {
    //
    // more or less uniqued history:
    //
    pub svg_hist_uniq: Vec<String>,       // each entry is an SVG
    pub summary_hist_uniq: Vec<String>,   // each entry is a summary
    pub input1_hist_uniq: Vec<String>,    // each entry is the originating command 1
    pub input2_hist_uniq: Vec<String>,    // each entry is the originating command 2
    pub narrative_hist_uniq: Vec<String>, // each entry is the narrative
    pub translated_input_hist_uniq: Vec<String>, // each entry is the translated originating command
    pub displayed_tables_hist_uniq: Vec<String>, // each entry is the tables that are displayed
    pub table_comp_hist_uniq: Vec<Vec<u8>>, // each entry is the compressed list of all tables
    pub last_widths_hist_uniq: Vec<Vec<u32>>,
    pub descrip_hist_uniq: Vec<String>, // descriptions (not used yet)
    //
    // parallel vectors, with one entry for each command entered in the text box:
    //
    pub svg_history: Vec<u32>,              // each entry is an SVG
    pub summary_history: Vec<u32>,          // each entry is a summary
    pub input1_history: Vec<u32>,           // each entry is the originating command 1
    pub input2_history: Vec<u32>,           // each entry is the originating command 2
    pub narrative_history: Vec<u32>,        // each entry is the narrative
    pub translated_input_history: Vec<u32>, // each entry is the translated originating command
    pub displayed_tables_history: Vec<u32>, // each entry is the tables that are displayed
    pub table_comp_history: Vec<u32>,       // each entry is the compressed list of all tables
    pub last_widths_history: Vec<u32>,      // last widths for clonotype pictures
    pub is_blank: Vec<bool>,                // set if the SVG is empty
    pub descrip_history: Vec<u32>,          // each entry is the description (not used yet)
    //
    // index of "current" position in those vectors, plus one:
    //
    pub history_index: u32, // this is saved but we ignore it, which seems better
    //
    // name of this session and narrative
    //
    pub name_value: String,
    pub orig_name_value: String,
    pub narrative: String, // not used yet
    //
    // origin of this session (if shared)
    //
    pub origin: String,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// clean: remove unused elements of hist_uniq vectors

pub fn clean_history<T: Clone>(hist_uniq: &mut Vec<T>, history: &mut Vec<u32>) {
    let mut used = vec![false; hist_uniq.len()];
    for i in 0..history.len() {
        used[history[i] as usize] = true;
    }
    let mut to_new = vec![0; hist_uniq.len()];
    let mut j = 0;
    for i in 0..hist_uniq.len() {
        if used[i] {
            if i != j {
                hist_uniq[j] = hist_uniq[i].clone();
            }
            to_new[i] = j;
            j += 1;
        }
    }
    hist_uniq.truncate(j);
    for i in 0..history.len() {
        history[i] = to_new[history[i] as usize] as u32;
    }
}

impl EncloneVisualHistory {
    pub fn clean_history(&mut self) {
        clean_history(&mut self.svg_hist_uniq, &mut self.svg_history);
        clean_history(&mut self.summary_hist_uniq, &mut self.summary_history);
        clean_history(&mut self.input1_hist_uniq, &mut self.input1_history);
        clean_history(&mut self.input2_hist_uniq, &mut self.input2_history);
        clean_history(&mut self.narrative_hist_uniq, &mut self.narrative_history);
        clean_history(
            &mut self.translated_input_hist_uniq,
            &mut self.translated_input_history,
        );
        clean_history(
            &mut self.displayed_tables_hist_uniq,
            &mut self.displayed_tables_history,
        );
        clean_history(&mut self.table_comp_hist_uniq, &mut self.table_comp_history);
        clean_history(
            &mut self.last_widths_hist_uniq,
            &mut self.last_widths_history,
        );
        clean_history(&mut self.descrip_hist_uniq, &mut self.descrip_history);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

const ENCLONE_VISUAL_HISTORY_VERSION: usize = 1;
const HEADER_LENGTH: usize = 40;
const NAME_BYTES: usize = 160;

impl EncloneVisualHistory {
    //
    // requirements for saving and restoring history
    //
    // 1. File structure must not change, except when we explicitly change it.
    // 2. Must be able to read old files.
    // 3. Must be versioned.
    // 4. Must have text header.
    // 5. Must store reasonably efficiently.
    // 6. Must be able to extract command list and name without reading entire file.
    // 7. Name can be changed without rewriting entire file.
    // 8. File structure must be guaranteed the same across environments (so, e.g., no usize).
    //
    // structure description
    //
    // 1. text header = 40 bytes = "enclone visual history file version ***\n",
    //    where *** is a positive integer, padded on the right with blanks
    // 2. total bytes in file (u32)
    // 3. total bytes in narrative (u32)
    // 4. total bytes up through translated_hist_uniq (u32)
    // 5. data for each of the fields in EncloneVisualHistory, with translated_input_history and
    //    translated_input_hist_uniq stored first, and narrative stored last (as raw bytes).

    pub fn save_as_bytes(&self) -> Vec<u8> {
        let mut bytes = format!(
            "enclone visual history file version{:<4}\n",
            ENCLONE_VISUAL_HISTORY_VERSION
        )
        .as_bytes()
        .to_vec();
        if ENCLONE_VISUAL_HISTORY_VERSION == 1 {
            bytes.append(&mut vec![0 as u8; 12]);
            let mut name_bytes = vec![0 as u8; NAME_BYTES];
            for i in 0..std::cmp::min(NAME_BYTES, self.name_value.as_bytes().len()) {
                name_bytes[i] = self.name_value.as_bytes()[i];
            }
            bytes.append(&mut name_bytes);
            bytes.append(&mut save_vec_u32(&self.translated_input_history));
            bytes.append(&mut save_vec_string(&self.translated_input_hist_uniq));
            bytes.append(&mut save_string(&self.origin));
            let b = u32_bytes(bytes.len());
            for i in 0..4 {
                bytes[HEADER_LENGTH + 4 + i] = b[i];
            }
            bytes.append(&mut save_vec_string_comp(&self.svg_hist_uniq));
            bytes.append(&mut save_vec_string(&self.summary_hist_uniq));
            bytes.append(&mut save_vec_string(&self.input1_hist_uniq));
            bytes.append(&mut save_vec_string(&self.input2_hist_uniq));
            bytes.append(&mut save_vec_string(&self.narrative_hist_uniq));
            bytes.append(&mut save_vec_string_comp(&self.displayed_tables_hist_uniq));
            bytes.append(&mut save_vec_vec_u8(&self.table_comp_hist_uniq));
            bytes.append(&mut save_vec_vec_u32(&self.last_widths_hist_uniq));
            bytes.append(&mut save_vec_string(&self.descrip_hist_uniq));
            bytes.append(&mut save_vec_u32(&self.svg_history));
            bytes.append(&mut save_vec_u32(&self.summary_history));
            bytes.append(&mut save_vec_u32(&self.input1_history));
            bytes.append(&mut save_vec_u32(&self.input2_history));
            bytes.append(&mut save_vec_u32(&self.narrative_history));
            bytes.append(&mut save_vec_u32(&self.displayed_tables_history));
            bytes.append(&mut save_vec_u32(&self.table_comp_history));
            bytes.append(&mut save_vec_u32(&self.last_widths_history));
            bytes.append(&mut save_vec_bool(&self.is_blank));
            bytes.append(&mut save_vec_u32(&self.descrip_history));
            bytes.append(&mut save_u32(self.history_index));
            bytes.append(&mut self.narrative.as_bytes().to_vec());
            let b = u32_bytes(bytes.len());
            for i in 0..4 {
                bytes[HEADER_LENGTH + i] = b[i];
            }
            let b = u32_bytes(self.narrative.len());
            for i in 0..4 {
                bytes[HEADER_LENGTH + 8 + i] = b[i];
            }
        }
        bytes
    }

    pub fn restore_from_bytes(bytes: &Vec<u8>) -> Result<Self, ()> {
        if bytes.len() < HEADER_LENGTH + 12 {
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
        let narrative_len = u32_from_bytes(&bytes[HEADER_LENGTH + 8..HEADER_LENGTH + 12]) as usize;
        let mut pos = HEADER_LENGTH + 12;
        let mut name_value_bytes = Vec::<u8>::new();
        for i in 0..NAME_BYTES {
            if bytes[pos + i] == 0 {
                break;
            }
            name_value_bytes.push(bytes[pos + i]);
        }
        let res = String::from_utf8(name_value_bytes);
        if res.is_err() {
            return Err(());
        }
        let name_value = res.unwrap();
        pos += NAME_BYTES;
        let translated_input_history = restore_vec_u32(&bytes, &mut pos)?;
        let translated_input_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let origin = restore_string(&bytes, &mut pos)?;
        let svg_hist_uniq = restore_vec_string_comp(&bytes, &mut pos)?;
        let summary_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let input1_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let input2_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let narrative_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let displayed_tables_hist_uniq = restore_vec_string_comp(&bytes, &mut pos)?;
        let table_comp_hist_uniq = restore_vec_vec_u8(&bytes, &mut pos)?;
        let last_widths_hist_uniq = restore_vec_vec_u32(&bytes, &mut pos)?;
        let descrip_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
        let svg_history = restore_vec_u32(&bytes, &mut pos)?;
        let summary_history = restore_vec_u32(&bytes, &mut pos)?;
        let input1_history = restore_vec_u32(&bytes, &mut pos)?;
        let input2_history = restore_vec_u32(&bytes, &mut pos)?;
        let narrative_history = restore_vec_u32(&bytes, &mut pos)?;
        let displayed_tables_history = restore_vec_u32(&bytes, &mut pos)?;
        let table_comp_history = restore_vec_u32(&bytes, &mut pos)?;
        let last_widths_history = restore_vec_u32(&bytes, &mut pos)?;
        let is_blank = restore_vec_bool(&bytes, &mut pos)?;
        let descrip_history = restore_vec_u32(&bytes, &mut pos)?;
        let history_index = restore_u32(&bytes, &mut pos)?;
        if pos + narrative_len != bytes.len() {
            return Err(());
        }
        let narrative = String::from_utf8(bytes[pos..].to_vec());
        if narrative.is_err() {
            return Err(());
        }
        Ok(EncloneVisualHistory {
            translated_input_history: translated_input_history,
            translated_input_hist_uniq: translated_input_hist_uniq,
            svg_hist_uniq: svg_hist_uniq,
            summary_hist_uniq: summary_hist_uniq,
            input1_hist_uniq: input1_hist_uniq,
            input2_hist_uniq: input2_hist_uniq,
            narrative_hist_uniq: narrative_hist_uniq,
            displayed_tables_hist_uniq: displayed_tables_hist_uniq,
            table_comp_hist_uniq: table_comp_hist_uniq,
            last_widths_hist_uniq: last_widths_hist_uniq,
            descrip_hist_uniq: descrip_hist_uniq,
            svg_history: svg_history,
            summary_history: summary_history,
            input1_history: input1_history,
            input2_history: input2_history,
            narrative_history: narrative_history,
            displayed_tables_history: displayed_tables_history,
            table_comp_history: table_comp_history,
            last_widths_history: last_widths_history,
            is_blank: is_blank,
            descrip_history: descrip_history,
            history_index: history_index,
            name_value: name_value.clone(),
            orig_name_value: name_value,
            origin: origin,
            narrative: narrative.unwrap(),
        })
    }
}

pub fn rewrite_name(filename: &str, name: &str) -> Result<(), std::io::Error> {
    let mut f = OpenOptions::new().write(true).read(true).open(&filename)?;
    f.seek(SeekFrom::Start((HEADER_LENGTH + 12) as u64))?;
    let mut name_bytes = vec![0 as u8; NAME_BYTES];
    for i in 0..std::cmp::min(NAME_BYTES, name.as_bytes().len()) {
        name_bytes[i] = name.as_bytes()[i];
    }
    f.write_all(&name_bytes)?;
    Ok(())
}

pub fn rewrite_narrative(filename: &str, narrative: &str) -> Result<(), std::io::Error> {
    let mut f = OpenOptions::new().write(true).read(true).open(&filename)?;
    let mut buf = vec![0 as u8; HEADER_LENGTH + 12];
    let res = f.read_exact(&mut buf);
    if res.is_err() {
        return Err(Error::new(ErrorKind::Other, "file appears truncated"));
    }
    let total = u32_from_bytes(&buf[HEADER_LENGTH..HEADER_LENGTH + 4]) as usize;
    let narrative_length = u32_from_bytes(&buf[HEADER_LENGTH + 8..HEADER_LENGTH + 12]) as usize;
    if narrative_length > total {
        return Err(Error::new(
            ErrorKind::Other,
            "narrative length exceeds total length",
        ));
    }
    let total_less_narrative = (total - narrative_length) as u64;
    f.seek(SeekFrom::Start(total_less_narrative))?;
    write!(f, "{}", narrative)?;
    let new_total = total + narrative.len() - narrative_length;
    if new_total < total {
        f.set_len(new_total as u64)?;
    }
    let new_narrative_length = narrative.len();
    f.seek(SeekFrom::Start(0))?;
    let x = u32_bytes(new_total);
    for i in 0..4 {
        buf[HEADER_LENGTH + i] = x[i];
    }
    let x = u32_bytes(new_narrative_length);
    for i in 0..4 {
        buf[HEADER_LENGTH + 8 + i] = x[i];
    }
    f.write_all(&buf)?;
    Ok(())
}

pub fn read_metadata(filename: &str) -> Result<(Vec<String>, String, String, String), ()> {
    let total;
    let n;
    let narrative_length;
    let narrative;
    {
        let mut f = open_for_read![&filename];
        let mut buf = vec![0 as u8; HEADER_LENGTH + 12];
        let res = f.read_exact(&mut buf);
        if res.is_err() {
            return Err(());
        }
        total = u32_from_bytes(&buf[HEADER_LENGTH..HEADER_LENGTH + 4]) as usize;
        n = u32_from_bytes(&buf[HEADER_LENGTH + 4..HEADER_LENGTH + 8]) as usize;
        narrative_length = u32_from_bytes(&buf[HEADER_LENGTH + 8..HEADER_LENGTH + 12]) as usize;
        if narrative_length > total {
            return Err(());
        }
        let total_less_narrative = (total - narrative_length) as u64;
        let res = f.seek(SeekFrom::Start(total_less_narrative));
        if res.is_err() {
            return Err(());
        }
        let mut buf = vec![0 as u8; narrative_length];
        let res = f.read_exact(&mut buf);
        if res.is_err() {
            return Err(());
        }
        let s = String::from_utf8(buf);
        if s.is_err() {
            return Err(());
        }
        narrative = s.unwrap();
    }
    let mut bytes = vec![0 as u8; n];
    let mut f = open_for_read![&filename];
    let res = f.read_exact(&mut bytes);
    if res.is_err() {
        return Err(());
    }
    let mut pos = HEADER_LENGTH + 12;
    let mut name_value_bytes = Vec::<u8>::new();
    for i in 0..NAME_BYTES {
        if bytes[pos + i] == 0 {
            break;
        }
        name_value_bytes.push(bytes[pos + i]);
    }
    let name = String::from_utf8(name_value_bytes);
    if name.is_err() {
        return Err(());
    }
    let name = name.unwrap();
    pos += NAME_BYTES;
    let translated_input_history = restore_vec_u32(&bytes, &mut pos)?;
    let translated_input_hist_uniq = restore_vec_string(&bytes, &mut pos)?;
    let mut commands = Vec::<String>::new();
    for i in translated_input_history.iter() {
        let i = *i as usize;
        if i >= translated_input_hist_uniq.len() {
            return Err(());
        }
        commands.push(translated_input_hist_uniq[i].clone());
    }
    let origin = restore_string(&bytes, &mut pos)?;
    Ok((commands, name, origin, narrative))
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
    let mut evh = evh.clone();
    evh.clean_history();
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Testing functions.

impl EncloneVisualHistory {
    pub fn save_restore_works(&self) -> bool {
        let bytes = self.save_as_bytes();
        let new = EncloneVisualHistory::restore_from_bytes(&bytes);
        if new.is_err() {
            return false;
        }
        *self == new.unwrap()
    }
}

pub fn test_evh_read_write(evh: &EncloneVisualHistory, filename: &str) {
    write_enclone_visual_history(&evh, &filename).unwrap();
    let mut evh = evh.clone();
    evh.clean_history();
    let evh2 = read_enclone_visual_history(&filename).unwrap();
    if evh != evh2 {
        eprintln!("");
        if evh.svg_hist_uniq != evh2.svg_hist_uniq {
            eprintln!("svg_hist_uniq changed");
            if evh.svg_hist_uniq.len() != evh2.svg_hist_uniq.len() {
                eprintln!(
                    "length changed from {} to {}",
                    evh.svg_hist_uniq.len(),
                    evh2.svg_hist_uniq.len()
                );
            }
        }
        if evh.history_index != evh2.history_index {
            eprintln!("history index changed");
        }
        if evh.orig_name_value != evh2.orig_name_value {
            eprintln!("orig name changed");
        }
        if evh.name_value != evh2.name_value {
            eprintln!("name changed");
        }
        if evh.narrative != evh2.narrative {
            eprintln!("narrative changed");
        }
        if evh.origin != evh2.origin {
            eprintln!("origin changed");
        }
        panic!("evh != evh2 for {}", filename);
    }
    for name in ["gerbilspit the great".to_string(), "shorter".to_string()].iter() {
        rewrite_name(&filename, &*name).unwrap();
        let evh2 = read_enclone_visual_history(&filename).unwrap();
        assert_eq!(evh2.name_value, *name);
    }
}
