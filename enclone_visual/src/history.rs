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
    // 2. nbytes for truncated file including just the first two fields (u32)
    // 3. nbytes for total file (u32)
    // 2. number of fields
    // 3. sequence of fields of the form
    //    (a) nbytes for member name (u8)
    //    (b) member name bytes
    //    (c) nbytes for data (u32)
    //    (d) data.
    // And the first two fields are translated_input_history, translated_input_hist_uniq.
    //
    // supported subtypes
    // - u32
    // - Vec<u32>
    // - Vec<bool>
    // - Vec<String>
    // - Vec<Vec<u8>>
    // - Vec<Vec<u32>>

    pub fn save_as_bytes(&self) -> Vec<u8> {
        serde_json::to_string(&self).unwrap().as_bytes().to_vec()
    }

    pub fn restore_from_bytes(bytes: &[u8]) -> Self {
        serde_json::from_str(&strme(&bytes)).unwrap()
    }
}
