// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_clipboard;
use crate::history::*;
use crate::messages::*;
use crate::testsuite::TESTS;
use crate::*;
use chrono::prelude::*;
use enclone_core::combine_group_pics::*;
use flate2::read::GzDecoder;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use io_utils::*;
use itertools::Itertools;
use std::env;
use std::io::Read;
use std::time::{Duration, Instant};
use vector_utils::*;

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::Restore(_) => Command::none(),
            Message::Resize(width, height) => {
                self.width = width;
                self.height = height;
                Command::none()
            }
            Message::GroupClicked(_message) => {
                let group_id = GROUP_ID.load(SeqCst);
                self.input_value = format!("{}", group_id);
                self.input1_value = format!("{}", group_id);
                self.input2_value.clear();
                GROUP_ID_CLICKED_ON.store(false, SeqCst);
                Command::perform(noop0(), Message::SubmitButtonPressed)
            }
            Message::SubmitButtonPressed(_) => {
                let mut group_spec = true;
                let mut group_ids = Vec::<usize>::new();
                let s = self.input_value.split(',').collect::<Vec<&str>>();
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
                        if self.h.input1_history.is_empty() && self.h.input2_history.is_empty() {
                            reply_text = "Group identifier can only be supplied if another \
                                command has already been run."
                                .to_string();
                        } else if group_ids[group_ids.len() - 1] > self.current_tables.len() {
                            reply_text = "Group identifier is too large.".to_string();
                        } else {
                            let mut group_pics = Vec::<String>::new();
                            let mut last_widths = Vec::<u32>::new();
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
                        let mut args2 = Vec::<String>::new();
                        for x in args.iter() {
                            if x.len() > 0 && !x.starts_with("G=") {
                                args2.push(x.to_string());
                            }
                        }
                        args2.push(format!("G={}", self.translated_input_value));
                        self.output_value = reply_text.to_string();
                        let hi = self.h.history_index;
                        self.h
                            .input1_history
                            .insert(hi as usize, self.h.input1_hist_uniq.len() as u32);
                        self.h.input1_hist_uniq.push(self.input1_value.clone());
                        self.h
                            .input2_history
                            .insert(hi as usize, self.h.input2_hist_uniq.len() as u32);
                        self.h.input2_hist_uniq.push(self.input2_value.clone());
                        self.h
                            .translated_input_history
                            .insert(hi as usize, self.h.translated_input_hist_uniq.len() as u32);
                        self.h
                            .translated_input_hist_uniq
                            .push(args2.iter().format(" ").to_string());
                        self.h
                            .svg_history
                            .insert(hi as usize, self.h.svg_history[(hi - 1) as usize]);
                        self.h
                            .summary_history
                            .insert(hi as usize, self.h.summary_history[(hi - 1) as usize]);
                        self.h
                            .displayed_tables_history
                            .insert(hi as usize, self.h.displayed_tables_hist_uniq.len() as u32);
                        self.h
                            .displayed_tables_hist_uniq
                            .push(reply_text.to_string());
                        self.h
                            .table_comp_history
                            .insert(hi as usize, self.h.table_comp_history[(hi - 1) as usize]);
                        self.h
                            .last_widths_history
                            .insert(hi as usize, self.h.last_widths_history[(hi - 1) as usize]);
                        self.h.is_blank.insert(hi as usize, self.is_blank_current());
                        self.h.history_index += 1;
                        if !TEST_MODE.load(SeqCst) {
                            Command::none()
                        } else {
                            Command::perform(noop0(), Message::Capture)
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
                self.h.history_index -= 1;
                let x = self.svg_current();
                self.post_svg(&x);
                self.summary_value = self.summary_current();
                self.output_value = self.displayed_tables_current();
                self.table_comp_value = self.table_comp_current();
                self.last_widths_value = self.last_widths_current();
                self.input1_value = self.input1_current();
                self.input2_value = self.input2_current();
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
                    Command::perform(noop0(), Message::Capture)
                }
            }

            Message::DelButtonPressed(_) => {
                let h = self.h.history_index - 1;
                self.h.svg_history.remove(h as usize);
                self.h.summary_history.remove(h as usize);
                self.h.input1_history.remove(h as usize);
                self.h.input2_history.remove(h as usize);
                self.h.translated_input_history.remove(h as usize);
                self.h.displayed_tables_history.remove(h as usize);
                self.h.table_comp_history.remove(h as usize);
                self.h.last_widths_history.remove(h as usize);
                self.h.is_blank.remove(h as usize);
                if self.state_count() == 0 {
                    self.h.history_index -= 1;
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
                } else {
                    if h > 0 {
                        self.h.history_index -= 1;
                    }
                    let x = self.svg_current();
                    self.post_svg(&x);
                    self.summary_value = self.summary_current();
                    self.output_value = self.displayed_tables_current();
                    self.table_comp_value = self.table_comp_current();
                    self.last_widths_value = self.last_widths_current();
                    self.input1_value = self.input1_current();
                    self.input2_value = self.input2_current();
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
                }
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop0(), Message::Capture)
                }
            }

            Message::ForwardButtonPressed(_) => {
                self.h.history_index += 1;
                let x = self.svg_current();
                self.post_svg(&x);
                self.summary_value = self.summary_current();
                self.output_value = self.displayed_tables_current();
                self.table_comp_value = self.table_comp_current();
                self.last_widths_value = self.last_widths_current();
                self.input1_value = self.input1_current();
                self.input2_value = self.input2_current();
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
                    Command::perform(noop0(), Message::Capture)
                }
            }

            Message::RunTests(_) => {
                let count = COUNT.load(SeqCst);
                if count == 0 {
                    self.window_id = get_window_id();
                }
                if count < TESTS.len() {
                    if TESTS[count].0.len() > 0 {
                        if TESTS[count].0.contains(" + ") {
                            self.input1_value = TESTS[count].0.before(" + ").to_string();
                            self.input2_value = TESTS[count].0.after(" + ").to_string();
                            self.input_value =
                                format!("{} {}", self.input1_value, self.input2_value);
                        } else {
                            self.input_value = TESTS[count].0.to_string();
                            self.input1_value = TESTS[count].0.to_string();
                            self.input2_value.clear();
                        }
                    }
                } else {
                    std::process::exit(0);
                }
                COUNT.store(COUNT.load(SeqCst) + 1, SeqCst);
                Command::perform(noop(), TESTS[count].1)
            }

            Message::HelpOpen(_) => {
                self.help_mode = true;
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

            Message::HelpClose(_) => {
                self.help_mode = false;
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

            Message::CookbookOpen => {
                self.cookbook_mode = true;
                Command::none()
            }

            Message::CookbookClose => {
                self.cookbook_mode = false;
                Command::none()
            }

            Message::SummaryOpen(_) => {
                self.summary_mode = true;
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

            Message::SummaryClose(_) => {
                self.summary_mode = false;
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

            Message::ConsoleOpen => {
                self.console_mode = true;
                Command::none()
            }

            Message::ConsoleClose => {
                self.console_mode = false;
                Command::none()
            }

            Message::ArchiveOpen => {
                self.archive_mode = true;
                Command::none()
            }

            Message::ArchiveClose => {
                self.archive_mode = false;
                Command::none()
            }

            Message::SaveOnExit => {
                self.save_on_exit = !self.save_on_exit;
                Command::none()
            }

            Message::Exit => {
                if true {
                    if self.save_on_exit {
                        let mut home = String::new();
                        for (key, value) in env::vars() {
                            if key == "HOME" {
                                home = value.clone();
                            }
                        }
                        if home.len() == 0 {
                            xprintln!(
                                "Weird, unable to determine your home directory, \
                                so Save on Exit failed.\n"
                            );
                            std::process::exit(1);
                        }
                        let enclone = format!("{}/enclone", home);
                        if !path_exists(&enclone) {
                            xprintln!(
                                "You do not have a directory ~/enclone, \
                                so Save on Exit failed.\n"
                            );
                            std::process::exit(1);
                        }
                        let dir = format!("{}/visual_history", enclone);
                        if !path_exists(&dir) {
                            let res = std::fs::create_dir(&dir);
                            if res.is_err() {
                                xprintln!(
                                    "Unable to create the directory \
                                    ~/enclone/visual_history, so Save on Exit Failed.\n"
                                );
                                std::process::exit(1);
                            }
                        }
                        let mut now = format!("{:?}", Local::now());
                        now = now.replace("T", "___");
                        now = now.before(".").to_string();
                        let filename = format!("{}/{}", dir, now);
                        let res = write_enclone_visual_history(&self.h, &filename);
                        if res.is_err() {
                            xprintln!(
                                "Was Unable to write history to the file {}, \
                                so Save on Exit Failed.\n",
                                filename
                            );
                            std::process::exit(1);
                        }
                    }
                    std::process::exit(0);
                }
                Command::none()
            }

            Message::InputChanged1(ref value) => {
                self.input1_value = value.to_string();
                self.input_value = self.input1_value.clone();
                if self.input1_value.len() > 0 && self.input2_value.len() > 0 {
                    self.input_value += " ";
                }
                self.input_value += &mut self.input2_value.clone();
                Command::none()
            }

            Message::InputChanged2(ref value) => {
                self.input2_value = value.to_string();
                self.input_value = self.input1_value.clone();
                if self.input1_value.len() > 0 && self.input2_value.len() > 0 {
                    self.input_value += " ";
                }
                self.input_value += &mut self.input2_value.clone();
                Command::none()
            }

            Message::ClearButtonPressed => {
                self.input1_value.clear();
                self.input2_value.clear();
                Command::none()
            }

            Message::ComputationDone(_) => {
                let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                if reply_text.contains("enclone failed") {
                    reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
                }
                if reply_text.len() == 0 {
                    if self.translated_input_value.contains(" NOPRINT") {
                        reply_text = "You used the NOPRINT option, so there are no \
                            clonotypes to see."
                            .to_string();
                    } else {
                        reply_text = "There are no clonotypes.  Please have a look at the summary."
                            .to_string();
                    }
                }

                // Start storing values.

                let reply_table_comp = SERVER_REPLY_TABLE_COMP.lock().unwrap()[0].clone();
                self.table_comp_value = reply_table_comp.clone();
                let hi = self.h.history_index;
                let len = self.h.table_comp_hist_uniq.len();
                if len > 0 && self.h.table_comp_hist_uniq[len - 1] == reply_table_comp {
                    self.h
                        .table_comp_history
                        .insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.table_comp_history.insert(hi as usize, len as u32);
                    self.h.table_comp_hist_uniq.push(reply_table_comp.clone());
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

                let len = self.h.last_widths_hist_uniq.len();
                if len > 0 && self.h.last_widths_hist_uniq[len - 1] == reply_last_widths {
                    self.h
                        .last_widths_history
                        .insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.last_widths_history.insert(hi as usize, len as u32);
                    self.h.last_widths_hist_uniq.push(reply_last_widths.clone());
                }
                let len = self.h.svg_hist_uniq.len();
                if len > 0 && self.h.svg_hist_uniq[len - 1] == reply_svg {
                    self.h.svg_history.insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.svg_history.insert(hi as usize, len as u32);
                    self.h.svg_hist_uniq.push(reply_svg.clone());
                }
                let len = self.h.summary_hist_uniq.len();
                if len > 0 && self.h.summary_hist_uniq[len - 1] == reply_summary {
                    self.h.summary_history.insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.summary_history.insert(hi as usize, len as u32);
                    self.h.summary_hist_uniq.push(reply_summary.clone());
                }
                let len = self.h.displayed_tables_hist_uniq.len();
                if len > 0 && self.h.displayed_tables_hist_uniq[len - 1] == reply_text {
                    self.h
                        .displayed_tables_history
                        .insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h
                        .displayed_tables_history
                        .insert(hi as usize, len as u32);
                    self.h.displayed_tables_hist_uniq.push(reply_text.clone());
                }
                let len = self.h.input1_hist_uniq.len();
                if len > 0 && self.h.input1_hist_uniq[len - 1] == self.input1_value {
                    self.h.input1_history.insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.input1_history.insert(hi as usize, len as u32);
                    self.h.input1_hist_uniq.push(self.input1_value.clone());
                }
                let len = self.h.input2_hist_uniq.len();
                if len > 0 && self.h.input2_hist_uniq[len - 1] == self.input2_value {
                    self.h.input2_history.insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.input2_history.insert(hi as usize, len as u32);
                    self.h.input2_hist_uniq.push(self.input2_value.clone());
                }
                let len = self.h.translated_input_hist_uniq.len();
                if len > 0
                    && self.h.translated_input_hist_uniq[len - 1] == self.translated_input_value
                {
                    self.h
                        .translated_input_history
                        .insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h
                        .translated_input_history
                        .insert(hi as usize, len as u32);
                    self.h
                        .translated_input_hist_uniq
                        .push(self.translated_input_value.clone());
                }
                self.h.is_blank.insert(hi as usize, blank);
                self.h.history_index += 1;
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
                xprintln!(
                    "total time to run command = {:.1} seconds",
                    elapsed(&self.start_command.unwrap())
                );
                let maxrss_self;
                unsafe {
                    let mut rusage: libc::rusage = std::mem::zeroed();
                    let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
                    assert_eq!(retval, 0);
                    maxrss_self = rusage.ru_maxrss;
                }
                let peak_mem_mb = maxrss_self as f64 / ((1024 * 1024) as f64);
                xprintln!(
                    "all time peak mem of this process is {:.1} MB\n",
                    peak_mem_mb
                );
                if VERBOSE.load(SeqCst) {
                    let mb = (1024 * 1024) as f64;
                    let mut total_svg = 0;
                    for x in self.h.svg_hist_uniq.iter() {
                        total_svg += x.len();
                    }
                    xprintln!("stored svgs = {:.1} MB", total_svg as f64 / mb);
                    let mut total_tables = 0;
                    for x in self.h.table_comp_hist_uniq.iter() {
                        total_tables += x.len();
                    }
                    xprintln!("stored tables = {:.1} MB", total_tables as f64 / mb);
                    xprintln!("");
                }

                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    self.sanity_check();
                    assert!(self.h.save_restore_works());
                    Command::perform(noop0(), Message::Capture)
                }
            }

            Message::Capture(_) => {
                let verbose = false;
                if verbose {
                    xprintln!("\ncapturing, input history:");
                    for i in 0..self.h.input1_history.len() {
                        let mark = if i + 1 == self.h.history_index as usize {
                            "*"
                        } else {
                            ""
                        };
                        xprintln!(
                            "[{}] {} {} {}",
                            i + 1,
                            self.h.input1_hist_uniq[self.h.input1_history[i] as usize],
                            self.h.input2_hist_uniq[self.h.input2_history[i] as usize],
                            mark
                        );
                    }
                }
                let count = COUNT.load(SeqCst);
                if count >= 1 {
                    capture(count, self.window_id);
                }
                Command::perform(noop0(), Message::RunTests)
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
