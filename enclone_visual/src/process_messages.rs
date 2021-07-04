// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_mac_clipboard;
use crate::messages::*;
use crate::testsuite::TESTS;
use crate::*;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use std::time::{Duration, Instant};

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::SubmitButtonPressed(_) => {
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    self.start_command = Some(Instant::now());
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

            Message::BackButtonPressed(_) => {
                self.history_index -= 1;
                let x = self.svg_history[self.history_index - 1].clone();
                self.post_svg(&x);
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
                let x = self.svg_history[self.history_index - 1].clone();
                self.post_svg(&x);
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
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::OpenModalCookbook => {
                COOKBOOK.store(true, SeqCst);
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
                self.input_value = self.command_history[self.history_index - 1].clone();
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
                reply_text += "\n \n \n"; // papering over truncation bug
                let mut reply_svg = String::new();
                if SERVER_REPLY_SVG.lock().unwrap().len() > 0 {
                    reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                    let mut blank = false;
                    if reply_svg.len() == 0 {
                        reply_svg = blank_svg();
                        blank = true;
                    }
                    if reply_svg.len() > 0 && self.input_value.parse::<usize>().is_err() {
                        self.svg_history.push(reply_svg.clone());
                        self.history_index += 1;
                        self.command_history
                            .push(self.translated_input_value.clone());
                        self.is_blank.push(blank);
                    }
                }
                self.output_value = reply_text.to_string();
                self.svg_value = reply_svg.to_string();
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
                copy_bytes_to_mac_clipboard(
                    &self.command_history[self.history_index - 1].as_bytes(),
                );
                Command::none()
            }

            Message::DoNothing => Command::none(),
        }
    }
}
