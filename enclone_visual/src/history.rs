// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

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
    pub last_widths_hist_uniq: Vec<Vec<usize>>,
    //
    // parallel vectors, with one entry for each command entered in the text box:
    //
    pub svg_history: Vec<usize>,              // each entry is an SVG
    pub summary_history: Vec<usize>,          // each entry is a summary
    pub input1_history: Vec<usize>,           // each entry is the originating command 1
    pub input2_history: Vec<usize>,           // each entry is the originating command 2
    pub translated_input_history: Vec<usize>, // each entry is the translated originating command
    pub displayed_tables_history: Vec<usize>, // each entry is the tables that are displayed
    pub table_comp_history: Vec<usize>,       // each entry is the compressed list of all tables
    pub last_widths_history: Vec<usize>,
    pub is_blank: Vec<bool>, // set if the SVG is empty
    //
    // index of "current" position in those vectors, plus one:
    //
    pub history_index: usize,
}

impl EncloneVisualHistory {
    pub fn save_as_bytes(&self) -> Vec<u8> {
        serde_json::to_string(&self).unwrap().as_bytes().to_vec()
    }

    pub fn restore_from_bytes(bytes: &[u8]) -> Self {
        serde_json::from_str(&strme(&bytes)).unwrap()
    }
}
