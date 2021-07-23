// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use canvas_view::CanvasView;
use iced::{button, scrollable, text_input, Color};
// use iced::Subscription;
// use iced_native::{window, Event};
use std::collections::HashMap;
use std::time::Instant;

#[derive(PartialEq)]
pub enum ComputeState {
    WaitingForRequest,
    Thinking,
}

impl Default for ComputeState {
    fn default() -> ComputeState {
        WaitingForRequest
    }
}

use ComputeState::*;

#[derive(Default)]
pub struct EncloneVisual {
    pub scroll: scrollable::State,
    pub input1: text_input::State,
    pub input2: text_input::State,
    pub input1_value: String,
    pub input2_value: String,
    pub input_value: String,
    pub translated_input_value: String,
    pub output_value: String,
    pub svg_value: String,
    pub png_value: Vec<u8>,
    pub summary_value: String,
    pub table_comp_value: Vec<u8>,
    pub last_widths_value: Vec<usize>,
    pub submit_button_text: String,
    // pub should_exit: bool,
    pub compute_state: ComputeState,
    pub copy_image_button_color: Color,
    pub canvas_view: CanvasView,
    pub cookbook: HashMap<String, String>,
    pub window_id: usize,
    pub start_command: Option<Instant>,
    pub help_mode: bool,
    pub cookbook_mode: bool,
    pub summary_mode: bool,
    //
    // current tables: suboptimal, as it would be better to keep some sort of vector of compressed
    // strings (allowing for compression to extend across the vector); see also
    // table_comp_hist_uniq[i], which is compressed, but doesn't allow access to individual entries
    //
    pub current_tables: Vec<String>,
    //
    // button states:
    //
    pub button: button::State,
    pub back_button: button::State,
    pub forward_button: button::State,
    pub del_button: button::State,
    pub exec_button: button::State,
    pub summary_button: button::State,
    pub open_state_cookbook: button::State,
    pub exit_state: button::State,
    pub copy_image_button: button::State,
    pub command_copy_button: button::State,
    pub null_button1: button::State,
    pub null_button2: button::State,
    pub null_button3: button::State,
    pub null_button: button::State,
    pub state_pos_button_null: button::State,
    pub clear_button: button::State,
    pub open_state: button::State,
    pub help_button: button::State,
    pub cookbook_button: button::State,
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
    //
    // current window dimensions
    //
    pub width: u32,
    pub height: u32,
}

#[derive(Default)]
pub struct ModalState {
    pub cancel_state: button::State,
    pub cookbook: bool,
}

impl EncloneVisual {
    pub fn hi(&self) -> usize {
        self.history_index - 1
    }
    pub fn state_count(&self) -> usize {
        self.svg_history.len()
    }
    pub fn svg_current(&self) -> String {
        return self.svg_hist_uniq[self.svg_history[self.hi()]].clone();
    }
    pub fn summary_current(&self) -> String {
        return self.summary_hist_uniq[self.summary_history[self.hi()]].clone();
    }
    pub fn input1_current(&self) -> String {
        return self.input1_hist_uniq[self.input1_history[self.hi()]].clone();
    }
    pub fn input2_current(&self) -> String {
        return self.input2_hist_uniq[self.input2_history[self.hi()]].clone();
    }
    pub fn translated_input_current(&self) -> String {
        return self.translated_input_hist_uniq[self.translated_input_history[self.hi()]].clone();
    }
    pub fn displayed_tables_current(&self) -> String {
        return self.displayed_tables_hist_uniq[self.displayed_tables_history[self.hi()]].clone();
    }
    pub fn table_comp_current(&self) -> Vec<u8> {
        return self.table_comp_hist_uniq[self.table_comp_history[self.hi()]].clone();
    }
    pub fn last_widths_current(&self) -> Vec<usize> {
        return self.last_widths_hist_uniq[self.last_widths_history[self.hi()]].clone();
    }
    pub fn is_blank_current(&self) -> bool {
        return self.is_blank[self.hi()];
    }
    pub fn sanity_check(&self) {
        let n = self.svg_history.len();
        assert_eq!(n, self.summary_history.len());
        assert_eq!(n, self.input1_history.len());
        assert_eq!(n, self.input2_history.len());
        assert_eq!(n, self.translated_input_history.len());
        assert_eq!(n, self.displayed_tables_history.len());
        assert_eq!(n, self.table_comp_history.len());
        assert_eq!(n, self.last_widths_history.len());
        assert_eq!(n, self.is_blank.len());
    }
}
