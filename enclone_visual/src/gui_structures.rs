// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::history::*;
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
    pub last_widths_value: Vec<u32>,
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
    pub console_mode: bool,
    pub archive_mode: bool,
    pub save_on_exit: bool,
    pub enabled: bool,
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
    pub console_open_button: button::State,
    pub console_close_button: button::State,
    pub save_on_exit_button: button::State,
    pub archive_open_button: button::State,
    pub archive_close_button: button::State,
    //
    // history
    //
    pub h: EncloneVisualHistory,
    //
    // current window dimensions
    //
    pub width: u32,
    pub height: u32,
}

impl EncloneVisual {
    pub fn hi(&self) -> usize {
        self.h.history_index as usize - 1
    }
    pub fn state_count(&self) -> usize {
        self.h.svg_history.len()
    }
    pub fn svg_current(&self) -> String {
        return self.h.svg_hist_uniq[self.h.svg_history[self.hi()] as usize].clone();
    }
    pub fn summary_current(&self) -> String {
        return self.h.summary_hist_uniq[self.h.summary_history[self.hi()] as usize].clone();
    }
    pub fn input1_current(&self) -> String {
        return self.h.input1_hist_uniq[self.h.input1_history[self.hi()] as usize].clone();
    }
    pub fn input2_current(&self) -> String {
        return self.h.input2_hist_uniq[self.h.input2_history[self.hi()] as usize].clone();
    }
    pub fn translated_input_current(&self) -> String {
        return self.h.translated_input_hist_uniq
            [self.h.translated_input_history[self.hi()] as usize]
            .clone();
    }
    pub fn displayed_tables_current(&self) -> String {
        return self.h.displayed_tables_hist_uniq
            [self.h.displayed_tables_history[self.hi()] as usize]
            .clone();
    }
    pub fn table_comp_current(&self) -> Vec<u8> {
        return self.h.table_comp_hist_uniq[self.h.table_comp_history[self.hi()] as usize].clone();
    }
    pub fn last_widths_current(&self) -> Vec<u32> {
        return self.h.last_widths_hist_uniq[self.h.last_widths_history[self.hi()] as usize]
            .clone();
    }
    pub fn is_blank_current(&self) -> bool {
        return self.h.is_blank[self.hi()];
    }
    pub fn sanity_check(&self) {
        let n = self.h.svg_history.len();
        assert_eq!(n, self.h.summary_history.len());
        assert_eq!(n, self.h.input1_history.len());
        assert_eq!(n, self.h.input2_history.len());
        assert_eq!(n, self.h.translated_input_history.len());
        assert_eq!(n, self.h.displayed_tables_history.len());
        assert_eq!(n, self.h.table_comp_history.len());
        assert_eq!(n, self.h.last_widths_history.len());
        assert_eq!(n, self.h.is_blank.len());
    }
}
