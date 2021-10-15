// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::history::*;
use crate::messages::*;
use crate::*;
use canvas_view::CanvasView;
use chrono::prelude::*;
use enclone_core::packing::*;
use flate2::read::GzDecoder;
use iced::{button, scrollable, text_input, Color};
// use iced::Subscription;
// use iced_native::{window, Event};
use std::collections::HashMap;
use std::io::Read;
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

pub fn convert_bytes_to_string(bytes: &[u8]) -> String {
    base64::encode(&bytes)
}

pub fn convert_string_to_bytes(s: &str) -> Vec<u8> {
    base64::decode(&s).unwrap()
}

pub struct Summary {
    pub summary: String,
    pub dataset_names: Vec<String>,
    pub metrics: Vec<Vec<String>>,
    pub metric_selected: Vec<bool>,
    pub metrics_condensed: bool,
}

impl Summary {
    pub fn pack(&self) -> String {
        let mut bytes = Vec::<u8>::new();
        bytes.append(&mut save_string(&self.summary));
        bytes.append(&mut save_vec_string(&self.dataset_names));
        bytes.append(&mut save_vec_vec_string(&self.metrics));
        bytes.append(&mut save_vec_bool(&self.metric_selected));
        bytes.append(&mut save_bool(self.metrics_condensed));
        convert_bytes_to_string(&bytes)
    }

    pub fn unpack(s: &str) -> Self {
        let bytes = convert_string_to_bytes(&s);
        let mut pos = 0;
        let summary = restore_string(&bytes, &mut pos).unwrap();
        let dataset_names = restore_vec_string(&bytes, &mut pos).unwrap();
        let metrics = restore_vec_vec_string(&bytes, &mut pos).unwrap();
        let metric_selected = restore_vec_bool(&bytes, &mut pos).unwrap();
        let metrics_condensed = restore_bool(&bytes, &mut pos).unwrap();
        Summary {
            summary: summary,
            dataset_names: dataset_names,
            metrics: metrics,
            metric_selected: metric_selected,
            metrics_condensed: metrics_condensed,
        }
    }
}

use ComputeState::*;

#[derive(Default)]
pub struct EncloneVisual {
    pub modified: bool,
    pub scroll: scrollable::State,
    pub summary_scroll: scrollable::State,
    pub clonotypes_scroll: scrollable::State,
    pub input1: text_input::State,
    pub input2: text_input::State,
    pub input1_value: String,
    pub input2_value: String,
    pub narrative_value: String,
    pub input_value: String,
    pub translated_input_value: String,
    pub output_value: String,
    pub svg_value: String,
    pub png_value: Vec<u8>,
    pub summary_value: String,
    pub table_comp_value: Vec<u8>,
    pub last_widths_value: Vec<u32>,
    pub descrip_value: String,
    pub submit_button_text: String,
    // pub should_exit: bool,
    pub compute_state: ComputeState,
    pub copy_image_button_color: Color,
    pub snapshot_button_color: Color,
    pub graphic_snapshot_button_color: Color,
    pub summary_snapshot_button_color: Color,
    pub sanity_button_color: Color,
    pub clonotypes_copy_button_color: Color,
    pub tooltip_toggle_button_color: Color,
    pub descrips_copy_button_color: Color,
    pub canvas_view: CanvasView,
    pub cookbook: HashMap<String, String>,
    pub window_id: usize,
    pub start_command: Option<Instant>,
    pub help_mode: bool,
    pub cookbook_mode: bool,
    pub summary_mode: bool,
    pub console_mode: bool,
    pub archive_mode: bool,
    pub clonotypes_mode: bool,
    pub graphic_mode: bool,
    pub save: bool,
    pub save_in_progress: bool,
    pub save_on_exit: bool,
    pub shares: Vec<Share>,
    pub visual: String,
    pub meta_pos: usize,
    pub metric_selected: Vec<bool>,
    pub metrics_condensed: bool,
    pub snapshot_start: Option<Instant>,
    pub graphic_snapshot_start: Option<Instant>,
    pub summary_snapshot_start: Option<Instant>,
    pub sanity_check_start: Option<Instant>,
    pub alluvial_tables_for_spreadsheet: String,
    pub alluvial_tables_copy_button_color: Color,
    pub alluvial_reads_tables_for_spreadsheet: String,
    pub alluvial_reads_tables_copy_button_color: Color,
    pub descrips_for_spreadsheet: String,
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
    pub save_button: button::State,
    pub save_on_exit_button: button::State,
    pub archive_open_button: button::State,
    pub archive_close_button: button::State,
    pub archive_refresh_button: button::State,
    pub open_archive_doc_button: button::State,
    pub close_archive_doc_button: button::State,
    pub open_alluvial_reads_doc_button: button::State,
    pub close_alluvial_reads_doc_button: button::State,
    pub archive_name_change_button: Vec<button::State>,
    pub archive_narrative_button: Vec<button::State>,
    pub copy_archive_narrative_button: Vec<button::State>,
    pub copy_cookbook_narrative_button: Vec<button::State>,
    pub narrative_button: button::State,
    pub copy_narrative_button: button::State,
    pub cookbook_narrative_button: Vec<button::State>,
    pub summary_copy_button: button::State,
    pub copy_summary_button_color: Color,
    pub metric_button: Vec<button::State>,
    pub condense_metrics_button: button::State,
    pub copy_selected_metrics_button: button::State,
    pub copy_selected_metrics_button_color: Color,
    pub snapshot_button: button::State,
    pub graphic_snapshot_button: button::State,
    pub recompute_button: button::State,
    pub sanity_button: button::State,
    pub clonotypes_open_button: button::State,
    pub clonotypes_close_button: button::State,
    pub graphic_open_button: button::State,
    pub graphic_close_button: button::State,
    pub clonotypes_copy_button: button::State,
    pub tooltip_toggle_button: button::State,
    pub alluvial_tables_copy_button: button::State,
    pub alluvial_reads_tables_copy_button: button::State,
    pub descrips_copy_button: button::State,
    //
    // more
    //
    pub this_meta: Vec<Message>,
    pub save_name: String,
    pub cookbooks: Vec<Vec<u8>>,
    //
    // history
    //
    pub h: EncloneVisualHistory,
    //
    // archive information and logic
    //
    pub archive_dir: Option<String>,
    pub archive_list: Vec<String>,
    pub expand_archive_entry: Vec<bool>,
    pub expand_cookbook_entry: Vec<bool>,
    pub restore_requested: Vec<bool>,
    pub restore_cookbook_requested: Vec<bool>,
    pub delete_requested: Vec<bool>,
    pub deleted: Vec<bool>,
    pub restore_msg: Vec<String>,
    pub restore_cookbook_msg: Vec<String>,
    pub just_restored: bool,
    pub archived_command_list: Vec<Option<Vec<String>>>,
    pub cookbook_command_list: Vec<Option<Vec<String>>>,
    pub archive_name: Vec<text_input::State>,
    pub archive_name_value: Vec<String>,
    pub name: text_input::State,
    pub name_change_requested: bool,
    pub orig_archive_name: Vec<String>,
    pub archive_share_requested: Vec<bool>,
    pub archive_origin: Vec<String>,
    pub archive_narrative: Vec<String>,
    pub archive_doc_open: bool,
    pub alluvial_reads_doc_open: bool,
    pub share_start: Option<Instant>,
    pub archive_name_change_button_color: Vec<Color>,
    pub copy_archive_narrative_button_color: Vec<Color>,
    pub copy_cookbook_narrative_button_color: Vec<Color>,
    pub copy_narrative_button_color: Color,
    pub cookbook_name: Vec<String>,
    pub cookbook_narrative: Vec<String>,
    //
    // users for sharing
    //
    pub sharing_enabled: bool,
    pub user: Vec<text_input::State>,
    pub user_value: Vec<String>,
    pub user_selected: Vec<bool>,
    pub user_valid: Vec<bool>,
    pub do_share: bool,
    pub do_share_complete: bool,
    pub archive_refresh_button_color: Color,
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
    pub fn narrative_current(&self) -> String {
        return self.h.narrative_hist_uniq[self.h.narrative_history[self.hi()] as usize].clone();
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
    pub fn descrip_current(&self) -> String {
        return self.h.descrip_hist_uniq[self.h.descrip_history[self.hi()] as usize].clone();
    }
    pub fn is_blank_current(&self) -> bool {
        return self.h.is_blank[self.hi()];
    }
    pub fn sanity_check(&self) {
        let n = self.h.svg_history.len();
        assert_eq!(n, self.h.summary_history.len());
        assert_eq!(n, self.h.input1_history.len());
        assert_eq!(n, self.h.input2_history.len());
        assert_eq!(n, self.h.narrative_history.len());
        assert_eq!(n, self.h.translated_input_history.len());
        assert_eq!(n, self.h.displayed_tables_history.len());
        assert_eq!(n, self.h.table_comp_history.len());
        assert_eq!(n, self.h.last_widths_history.len());
        assert_eq!(n, self.h.descrip_history.len());
        assert_eq!(n, self.h.is_blank.len());
    }
    pub fn update_to_current(&mut self) {
        if self.h.history_index == 0 {
            self.input1_value.clear();
            self.input2_value.clear();
            self.svg_value.clear();
            self.png_value.clear();
            self.submit_button_text.clear();
            self.summary_value.clear();
            self.output_value.clear();
            self.table_comp_value.clear();
            self.last_widths_value.clear();
            self.translated_input_value.clear();
            self.current_tables.clear();
            self.descrip_value.clear();
        } else {
            let x = self.svg_current();
            self.svg_value = self.svg_current();
            self.post_svg(&x);
            self.summary_value = self.summary_current();
            self.output_value = self.displayed_tables_current();
            self.table_comp_value = self.table_comp_current();
            self.last_widths_value = self.last_widths_current();
            self.descrip_value = self.descrip_current();
            self.input1_value = self.input1_current();
            self.input2_value = self.input2_current();
            self.translated_input_value = self.translated_input_current();
            if self.table_comp_value.len() > 0 {
                let mut gunzipped = Vec::<u8>::new();
                let mut d = GzDecoder::new(&*self.table_comp_value);
                d.read_to_end(&mut gunzipped).unwrap();
                self.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
            } else {
                self.current_tables.clear();
            }
        }
    }
    pub fn save_as(&mut self, filename: &str) {
        let path = format!("{}/{}", self.archive_dir.as_ref().unwrap(), filename);
        let res = write_enclone_visual_history(&self.h, &path);
        if res.is_err() {
            xprintln!(
                "Was Unable to write history to the file {}, \
                so Save failed.\n",
                path
            );
            std::process::exit(1);
        }
        self.archive_list.insert(0, filename.to_string());
        self.restore_requested.insert(0, false);
        self.delete_requested.insert(0, false);
        self.deleted.insert(0, false);
        self.expand_archive_entry.insert(0, false);
        self.restore_msg.insert(0, String::new());
        self.archived_command_list.insert(0, None);
        self.archive_name
            .insert(0, iced::text_input::State::default());
        self.archive_name_value.insert(0, String::new());
        self.archive_name_change_button_color
            .insert(0, Color::from_rgb(0.0, 0.0, 0.0));
        self.copy_archive_narrative_button_color
            .insert(0, Color::from_rgb(0.0, 0.0, 0.0));
        self.archive_name_change_button
            .insert(0, button::State::default());
        self.archive_narrative_button
            .insert(0, button::State::default());
        self.copy_archive_narrative_button
            .insert(0, button::State::default());
        self.archive_share_requested.insert(0, false);
        self.archive_origin.insert(0, String::new());
        self.archive_narrative.insert(0, String::new());
        self.orig_archive_name.insert(0, String::new());
    }
    pub fn save(&mut self) {
        let mut now = format!("{:?}", Local::now());
        now = now.replace("T", "___");
        now = now.before(".").to_string();
        self.save_as(&now);
    }
}
