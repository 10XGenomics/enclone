// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_mac_clipboard;
use crate::messages::*;
use crate::testsuite::TESTS;
use crate::*;
use flate2::read::GzDecoder;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use std::io::Read;
use std::time::{Duration, Instant};

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::SubmitButtonPressed(_) => {
                if self.compute_state != WaitingForRequest {
                        Command::none()
                } else {
                    if self.input_value.parse::<usize>().is_ok() {
                        self.translated_input_value = self.input_value.clone();
                        let id = self.input_value.force_usize();
                        let mut reply_text;
                        if self.command_history.is_empty() {
                            reply_text = "Group identifier can only be supplied if another \
                                command has already been run.".to_string();
                        } else if id == 0 {
                            reply_text = "Group identifiers start at 1".to_string();
                        } else if id > self.current_tables.len() {
                            reply_text = "Group identifier is too large.".to_string();
                        } else {
                            reply_text = self.current_tables[id - 1].clone();
                            reply_text += "\n \n \n"; // papering over truncation bug in display
                        }
                        self.svg_history.push(self.svg_history[self.svg_history.len() - 1]);
                        self.summary_history.push(self.summary_history[self.summary_history.len() - 1]);
                        self.displayed_tables_history.push(self.displayed_tables_history[self.displayed_tables_history.len() - 1]);
                        self.table_comp_history.push(self.table_comp_history[self.table_comp_history.len() - 1]);
                        self.command_history.push(self.command_hist_uniq.len());
                        self.command_hist_uniq.push(self.input_value.clone());
                        self.is_blank.push(self.is_blank[self.is_blank.len() - 1]);
                        self.history_index += 1;
                        self.output_value = reply_text.to_string();
                        if !TEST_MODE.load(SeqCst) {
                            Command::none()
                        } else {
                            let count = COUNT.load(SeqCst);
                            if count > 1 {
                                capture(count, self.window_id);
                            }
                            Command::perform(noop(), Message::RunTests)
                        }
                    } else {
                        self.compute_state = Thinking;
                        self.start_command = Some(Instant::now());
                        // The following sleep is needed to get the button text to consistenly 
                        // update.
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
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    let count = COUNT.load(SeqCst);
                    if count > 1 {
                        capture(count, self.window_id);
                    }
                    Command::perform(noop(), Message::RunTests)
                }
            }

            Message::ForwardButtonPressed(_) => {
                self.history_index += 1;
                let x = self.svg_current();
                self.post_svg(&x);
                self.summary_value = self.summary_current();
                self.output_value = self.displayed_tables_current();
                self.table_comp_value = self.table_comp_current();
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    let count = COUNT.load(SeqCst);
                    if count > 1 {
                        capture(count, self.window_id);
                    }
                    Command::perform(noop(), Message::RunTests)
                }
            }

            Message::RunTests(_) => {
                let mut count = COUNT.load(SeqCst);
                if count == 0 {
                    self.window_id = get_window_id();
                }
                if count < TESTS.len() {
                    if TESTS[count].0.len() > 0 {
                        self.input_value = TESTS[count].0.to_string();
                    }
                } else {
                    count += 1;
                    capture(count, self.window_id);
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

            // same as above except for first line
            Message::ExecuteButtonPressed => {
                self.input_value = self.command_current();
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    // The following sleep is needed to get the button text to consistenly update.
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
                } else {
                    Command::none()
                }
            }

            Message::ComputationDone(_) => {
                let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                if reply_text.contains("enclone failed") {
                    reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
                }
                reply_text += "\n \n \n"; // papering over truncation bug in display
                let reply_summary = SERVER_REPLY_SUMMARY.lock().unwrap()[0].clone();
                let reply_table_comp = SERVER_REPLY_TABLE_COMP.lock().unwrap()[0].clone();
                let mut reply_svg = String::new();
                let mut blank = false;
                if SERVER_REPLY_SVG.lock().unwrap().len() > 0 {
                    reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                    if reply_svg.len() == 0 {
                        reply_svg = blank_svg();
                        blank = true;
                    }
                }

                // Store values.
                //
                // We want to push as little as possible onto the hist_uniq vectors,
                // and we want to do this as rapidly as possible.  The code here is not
                // optimal, for two reasons:
                // 1. We only compare to the last entry.
                // 2. We make comparisons in cases where we should already know the answer.

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
                let len = self.table_comp_hist_uniq.len();
                if len > 0 && self.table_comp_hist_uniq[len - 1] == reply_table_comp {
                    self.table_comp_history.push(len - 1);
                } else {
                    self.table_comp_history.push(len);
                    self.table_comp_hist_uniq.push(reply_table_comp.clone());
                    let mut gunzipped = Vec::<u8>::new();
                    let mut d = GzDecoder::new(&*reply_table_comp);
                    d.read_to_end(&mut gunzipped).unwrap();
                    self.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
                }
                let len = self.command_hist_uniq.len();
                if len > 0 && self.command_hist_uniq[len - 1] == self.translated_input_value {
                    self.command_history.push(len - 1);
                } else {
                    self.command_history.push(len);
                    self.command_hist_uniq
                        .push(self.translated_input_value.clone());
                }
                self.is_blank.push(blank);
                self.history_index += 1;
                self.output_value = reply_text.to_string();
                self.svg_value = reply_svg.to_string();
                self.summary_value = reply_summary.to_string();
                self.table_comp_value = reply_table_comp.clone();
                SUMMARY_CONTENTS.lock().unwrap().clear();
                SUMMARY_CONTENTS
                    .lock()
                    .unwrap()
                    .push(self.summary_value.clone());
                if self.svg_value.len() > 0 {
                    self.post_svg(&reply_svg);
                }
                self.compute_state = WaitingForRequest;
                println!(
                    "total time to run command = {:.1} seconds\n",
                    elapsed(&self.start_command.unwrap())
                );
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    let count = COUNT.load(SeqCst);
                    if count > 1 {
                        capture(count, self.window_id);
                    }
                    Command::perform(noop(), Message::RunTests)
                }
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
                copy_png_bytes_to_mac_clipboard(&self.png_value);
                let used = elapsed(&t);
                let extra = MIN_FLASH_SECONDS - used;
                if extra > 0.0 {
                    thread::sleep(Duration::from_millis((extra * 1000.0).round() as u64));
                }
                self.copy_image_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CommandCopyButtonPressed => {
                copy_bytes_to_mac_clipboard(&self.command_current().as_bytes());
                Command::none()
            }

            Message::DoNothing => Command::none(),
        }
    }
}
