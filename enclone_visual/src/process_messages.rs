// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_clipboard;
use crate::messages::*;
use crate::testsuite::TESTS;
use crate::*;
use enclone_core::combine_group_pics::*;
use flate2::read::GzDecoder;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use itertools::Itertools;
use std::io::Read;
use std::time::{Duration, Instant};
use vector_utils::*;

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::GroupClicked(_message) => {
                let group_id = GROUP_ID.load(SeqCst);
                self.input_value = format!("{}", group_id);
                GROUP_ID_CLICKED_ON.store(false, SeqCst);
                Command::perform(noop(), Message::SubmitButtonPressed)
            }
            Message::SubmitButtonPressed(_) => {
                let mut group_spec = true;
                let mut group_ids = Vec::<usize>::new();
                let s = &self.input_value.split(',').collect::<Vec<&str>>();
                for i in 0..s.len() {
                    let mut ok = false;
                    if s[i].parse::<usize>().is_ok() {
                        let n = s[i].force_usize();
                        if n >= 1 {
                            group_ids.push(n);
                            ok = true;
                        }
                    } else if s[i].contains("-") {
                        let (a, b) = (s[i].before("-"), s[i].after("-"));
                        if a.parse::<usize>().is_ok() && b.parse::<usize>().is_ok() {
                            let (a, b) = (a.force_usize(), b.force_usize());
                            if 1 <= a && a <= b {
                                for j in a..=b {
                                    group_ids.push(j);
                                }
                                ok = true;
                            }
                        }
                    }
                    if !ok {
                        group_spec = false;
                    }
                }
                if group_ids.is_empty() {
                    group_spec = false;
                }
                unique_sort(&mut group_ids);
                if self.compute_state != WaitingForRequest {
                    Command::none()
                } else {
                    if group_spec {
                        self.translated_input_value = self.input_value.clone();
                        let mut reply_text;
                        let new = self.translated_input_current();
                        let args = new.split(' ').collect::<Vec<&str>>();
                        if self.input_history.is_empty() {
                            reply_text = "Group identifier can only be supplied if another \
                                command has already been run."
                                .to_string();
                        } else if group_ids[group_ids.len() - 1] > self.current_tables.len() {
                            reply_text = "Group identifier is too large.".to_string();
                        } else {
                            let mut group_pics = Vec::<String>::new();
                            let mut last_widths = Vec::<usize>::new();
                            for x in group_ids.iter() {
                                group_pics.push(self.current_tables[*x - 1].clone());
                                last_widths.push(self.last_widths_value[*x - 1]);
                            }
                            reply_text = combine_group_pics(
                                &group_pics,
                                &last_widths,
                                args.contains(&"NOPRINT"),
                                true, // .noprintx
                                args.contains(&"HTML"),
                                args.contains(&"NGROUP"),
                                false, // .pretty
                            );
                            reply_text += "\n \n \n"; // papering over truncation bug in display
                        }
                        self.input_history.push(self.input_hist_uniq.len());
                        self.input_hist_uniq.push(self.input_value.clone());
                        let mut args2 = Vec::<String>::new();
                        for x in args.iter() {
                            if x.len() > 0 && !x.starts_with("G=") {
                                args2.push(x.to_string());
                            }
                        }
                        args2.push(format!("G={}", self.input_value));
                        self.translated_input_history
                            .push(self.translated_input_hist_uniq.len());
                        self.translated_input_hist_uniq
                            .push(args2.iter().format(" ").to_string());
                        self.svg_history.push(*self.svg_history.last().unwrap());
                        self.summary_history
                            .push(*self.summary_history.last().unwrap());
                        self.displayed_tables_history
                            .push(self.displayed_tables_hist_uniq.len());
                        self.displayed_tables_hist_uniq.push(reply_text.to_string());
                        self.table_comp_history
                            .push(*self.table_comp_history.last().unwrap());
                        self.last_widths_history
                            .push(*self.last_widths_history.last().unwrap());
                        self.is_blank.push(self.is_blank[self.is_blank.len() - 1]);
                        self.history_index = self.input_history.len();
                        self.output_value = reply_text.to_string();
                        if !TEST_MODE.load(SeqCst) {
                            Command::none()
                        } else {
                            Command::perform(noop(), Message::Capture)
                        }
                    } else {
                        self.compute_state = Thinking;
                        self.start_command = Some(Instant::now());
                        // The following sleep is needed for button text to consistenly update.
                        thread::sleep(Duration::from_millis(20));
                        if self.input_value.starts_with('#')
                            && self.cookbook.contains_key(&self.input_value)
                        {
                            self.translated_input_value = self.cookbook[&self.input_value].clone();
                        } else {
                            self.translated_input_value = self.input_value.clone();
                        }
                        USER_REQUEST.lock().unwrap().clear();
                        USER_REQUEST
                            .lock()
                            .unwrap()
                            .push(self.translated_input_value.clone());
                        PROCESSING_REQUEST.store(true, SeqCst);
                        Command::perform(compute(), Message::ComputationDone)
                    }
                }
            }

            Message::BackButtonPressed(_) => {
                self.history_index -= 1;
                let x = self.svg_current();
                self.post_svg(&x);
                self.summary_value = self.summary_current();
                self.output_value = self.displayed_tables_current();
                self.table_comp_value = self.table_comp_current();
                self.last_widths_value = self.last_widths_current();
                self.input_value = self.input_current();
                self.translated_input_value = self.translated_input_current();
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if self.table_comp_value.len() > 0 {
                    let mut gunzipped = Vec::<u8>::new();
                    let mut d = GzDecoder::new(&*self.table_comp_value);
                    d.read_to_end(&mut gunzipped).unwrap();
                    self.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
                } else {
                    self.current_tables.clear();
                }
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop(), Message::Capture)
                }
            }

            Message::ForwardButtonPressed(_) => {
                self.history_index += 1;
                let x = self.svg_current();
                self.post_svg(&x);
                self.summary_value = self.summary_current();
                self.output_value = self.displayed_tables_current();
                self.table_comp_value = self.table_comp_current();
                self.last_widths_value = self.last_widths_current();
                self.input_value = self.input_current();
                self.translated_input_value = self.translated_input_current();
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if self.table_comp_value.len() > 0 {
                    let mut gunzipped = Vec::<u8>::new();
                    let mut d = GzDecoder::new(&*self.table_comp_value);
                    d.read_to_end(&mut gunzipped).unwrap();
                    self.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
                } else {
                    self.current_tables.clear();
                }
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop(), Message::Capture)
                }
            }

            Message::RunTests(_) => {
                let count = COUNT.load(SeqCst);
                if count == 0 {
                    self.window_id = get_window_id();
                }
                if count < TESTS.len() {
                    if TESTS[count].0.len() > 0 {
                        self.input_value = TESTS[count].0.to_string();
                    }
                } else {
                    std::process::exit(0);
                }
                COUNT.store(COUNT.load(SeqCst) + 1, SeqCst);
                Command::perform(noop(), TESTS[count].1)
            }

            Message::OpenModalHelp => {
                COOKBOOK.store(false, SeqCst);
                SUMMARY.store(false, SeqCst);
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::OpenModalCookbook => {
                COOKBOOK.store(true, SeqCst);
                SUMMARY.store(false, SeqCst);
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::OpenModalSummary => {
                COOKBOOK.store(false, SeqCst);
                SUMMARY.store(true, SeqCst);
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::CloseModalHelp => {
                self.modal_state_help.show(false);
                Command::none()
            }

            Message::Exit => {
                if true {
                    std::process::exit(0);
                }
                Command::none()
            }

            Message::CancelButtonPressed => {
                self.modal_state_help.show(false);
                Command::none()
            }

            Message::InputChanged(ref value) => {
                self.input_value = value.to_string();
                Command::none()
            }

            Message::ClearButtonPressed => {
                self.input_value.clear();
                Command::none()
            }

            Message::ComputationDone(_) => {
                let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                if reply_text.contains("enclone failed") {
                    reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
                }
                if reply_text.len() == 0 {
                    reply_text = "Looks like you used the NOPRINT option, and there are no \
                        clonotypes to see."
                        .to_string();
                }

                // Start storing values.

                let reply_table_comp = SERVER_REPLY_TABLE_COMP.lock().unwrap()[0].clone();
                self.table_comp_value = reply_table_comp.clone();
                let len = self.table_comp_hist_uniq.len();
                if len > 0 && self.table_comp_hist_uniq[len - 1] == reply_table_comp {
                    self.table_comp_history.push(len - 1);
                } else {
                    self.table_comp_history.push(len);
                    self.table_comp_hist_uniq.push(reply_table_comp.clone());
                    if self.table_comp_value.len() > 0 {
                        let mut gunzipped = Vec::<u8>::new();
                        let mut d = GzDecoder::new(&*reply_table_comp);
                        d.read_to_end(&mut gunzipped).unwrap();
                        self.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
                    } else {
                        self.current_tables.clear();
                    }
                }

                // Keep going.

                reply_text += "\n \n \n"; // papering over truncation bug in display
                let reply_summary = SERVER_REPLY_SUMMARY.lock().unwrap()[0].clone();
                let reply_last_widths = SERVER_REPLY_LAST_WIDTHS.lock().unwrap()[0].clone();
                let mut reply_svg = String::new();
                let mut blank = false;
                if SERVER_REPLY_SVG.lock().unwrap().len() > 0 {
                    reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                    if reply_svg.len() == 0 {
                        reply_svg = blank_svg();
                        blank = true;
                    }
                }

                // Continue storing values.
                //
                // We want to push as little as possible onto the hist_uniq vectors,
                // and we want to do this as rapidly as possible.  The code here is not
                // optimal, for two reasons:
                // 1. We only compare to the last entry.
                // 2. We make comparisons in cases where we should already know the answer.

                let len = self.last_widths_hist_uniq.len();
                if len > 0 && self.last_widths_hist_uniq[len - 1] == reply_last_widths {
                    self.last_widths_history.push(len - 1);
                } else {
                    self.last_widths_history.push(len);
                    self.last_widths_hist_uniq.push(reply_last_widths.clone());
                }
                let len = self.svg_hist_uniq.len();
                if len > 0 && self.svg_hist_uniq[len - 1] == reply_svg {
                    self.svg_history.push(len - 1);
                } else {
                    self.svg_history.push(len);
                    self.svg_hist_uniq.push(reply_svg.clone());
                }
                let len = self.summary_hist_uniq.len();
                if len > 0 && self.summary_hist_uniq[len - 1] == reply_summary {
                    self.summary_history.push(len - 1);
                } else {
                    self.summary_history.push(len);
                    self.summary_hist_uniq.push(reply_summary.clone());
                }
                let len = self.displayed_tables_hist_uniq.len();
                if len > 0 && self.displayed_tables_hist_uniq[len - 1] == reply_text {
                    self.displayed_tables_history.push(len - 1);
                } else {
                    self.displayed_tables_history.push(len);
                    self.displayed_tables_hist_uniq.push(reply_text.clone());
                }
                let len = self.input_hist_uniq.len();
                if len > 0 && self.input_hist_uniq[len - 1] == self.input_value {
                    self.input_history.push(len - 1);
                } else {
                    self.input_history.push(len);
                    self.input_hist_uniq.push(self.input_value.clone());
                }
                let len = self.translated_input_hist_uniq.len();
                if len > 0
                    && self.translated_input_hist_uniq[len - 1] == self.translated_input_value
                {
                    self.translated_input_history.push(len - 1);
                } else {
                    self.translated_input_history.push(len);
                    self.translated_input_hist_uniq
                        .push(self.translated_input_value.clone());
                }
                self.is_blank.push(blank);
                self.history_index = self.input_history.len();
                self.output_value = reply_text.to_string();
                self.svg_value = reply_svg.to_string();
                self.summary_value = reply_summary.to_string();
                self.last_widths_value = reply_last_widths.clone();
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if self.svg_value.len() > 0 {
                    self.post_svg(&reply_svg);
                }
                self.compute_state = WaitingForRequest;
                eprintln!(
                    "total time to run command = {:.1} seconds\n",
                    elapsed(&self.start_command.unwrap())
                );
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    self.sanity_check();
                    Command::perform(noop(), Message::Capture)
                }
            }

            Message::Capture(_) => {
                let count = COUNT.load(SeqCst);
                if count >= 1 {
                    capture(count, self.window_id);
                }
                Command::perform(noop(), Message::RunTests)
            }

            // Catch exit (when the upper left red button is pushed) and store DONE to make
            // the server thread exit gracefully.  Otherwise you will get a an error message
            // and a traceback.
            /*
            Message::EventOccurred(ref event) => {
                if let Event::Window(window::Event::CloseRequested) = event {
                    DONE.store(true, SeqCst);
                    thread::sleep(Duration::from_millis(50));
                    self.should_exit = true;
                }
                Command::none()
            }
            */
            Message::GraphicsCopyButtonPressed => {
                self.copy_image_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(
                    flash_copy_image_button(),
                    Message::GraphicsCopyButtonFlashed,
                )
            }

            Message::GraphicsCopyButtonFlashed(_) => {
                // Convert to PNG and copy to clipboard, and flash the button for the maximum of
                // the conversion time and MIN_FLASH_SECONDS.
                const MIN_FLASH_SECONDS: f64 = 0.4;
                let t = Instant::now();
                if self.png_value.is_empty() {
                    self.png_value = convert_svg_to_png(&self.svg_value.as_bytes());
                }
                copy_png_bytes_to_clipboard(&self.png_value);
                let used = elapsed(&t);
                let extra = MIN_FLASH_SECONDS - used;
                if extra > 0.0 {
                    thread::sleep(Duration::from_millis((extra * 1000.0).round() as u64));
                }
                self.copy_image_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CommandCopyButtonPressed => {
                copy_bytes_to_clipboard(&self.translated_input_current().as_bytes());
                Command::none()
            }

            Message::DoNothing => Command::none(),
        }
    }
}
