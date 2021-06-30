// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use canvas_view::CanvasView;
use iced::{button, scrollable, text_input, Color};
// use iced::Subscription;
use iced_aw::modal;
// use iced_native::{window, Event};
use std::collections::HashMap;

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
    pub input: text_input::State,
    pub input_value: String,
    pub translated_input_value: String,
    pub output_value: String,
    pub svg_value: String,
    pub png_value: Vec<u8>,
    pub button: button::State,
    pub back_button: button::State,
    pub forward_button: button::State,
    pub exec_button: button::State,
    pub submit_button_text: String,
    pub open_state: button::State,
    pub open_state_cookbook: button::State,
    pub exit_state: button::State,
    pub modal_state_help: modal::State<ModalState>,
    // pub should_exit: bool,
    pub compute_state: ComputeState,
    pub copy_image_button: button::State,
    pub copy_image_button_color: Color,
    pub canvas_view: CanvasView,
    pub command_copy_button: button::State,
    pub null_button1: button::State,
    pub null_button2: button::State,
    pub null_button3: button::State,
    pub null_button: button::State,
    pub cookbook: HashMap<String, String>,
    pub clear_button: button::State,
    pub window_id: usize,

    // parallel vectors:
    pub svg_history: Vec<String>,
    pub command_history: Vec<String>,
    pub is_blank: Vec<bool>,

    // index of "current" position in those vectors:
    pub history_index: usize,
}

#[derive(Default)]
pub struct ModalState {
    pub cancel_state: button::State,
    pub cookbook: bool,
}
