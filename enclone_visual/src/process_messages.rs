// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_clipboard;
use crate::history::*;
use crate::messages::*;
use crate::proc1::*;
use crate::proc2::*;
use crate::summary::*;
use crate::testsuite::TESTS;
use crate::*;
use chrono::prelude::*;
use iced::{Color, Command};
use itertools::Itertools;
use std::time::{Duration, Instant};

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        MESSAGE_HISTORY
            .lock()
            .unwrap()
            .push(format!("{:?}", message));
        match message {
            Message::ArchiveSnapshot => {
                self.archive_snapshot_start = Some(Instant::now());
                self.archive_snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteArchiveSnapshot)
            }

            Message::CompleteArchiveSnapshot(_) => {
                snapshot(&self.archive_snapshot_start);
                self.archive_snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::ClonotypesSnapshot => {
                self.clonotypes_snapshot_start = Some(Instant::now());
                self.clonotypes_snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteClonotypesSnapshot)
            }

            Message::CompleteClonotypesSnapshot(_) => {
                snapshot(&self.clonotypes_snapshot_start);
                self.clonotypes_snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CommandSnapshot => {
                self.command_snapshot_start = Some(Instant::now());
                self.command_snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteCommandSnapshot)
            }

            Message::CompleteCommandSnapshot(_) => {
                snapshot(&self.command_snapshot_start);
                self.command_snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CommandOpen(_) => {
                self.command_mode = true;
                Command::none()
            }

            Message::CommandClose => {
                self.command_mode = false;
                Command::none()
            }

            Message::GraphicHelp => {
                if self.graphic_help_title == "Help" {
                    self.graphic_help_title = "Close help".to_string();
                } else {
                    self.graphic_help_title = "Help".to_string();
                }
                Command::none()
            }

            Message::GraphicPng => {
                self.summary_png_start = Some(Instant::now());
                self.png_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteGraphicPng)
            }

            Message::CompleteGraphicPng(_) => {
                if self.graphic_png_title == "PNG" {
                    if self.png_value.len() == 0 {
                        self.png_value = convert_svg_to_png(&self.svg_value.as_bytes(), 2000);
                    }
                }
                const MIN_SLEEP: f64 = 0.4;
                let used = elapsed(&self.summary_png_start.unwrap());
                if used < MIN_SLEEP {
                    let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
                    thread::sleep(Duration::from_millis(ms));
                }
                if self.graphic_png_title == "PNG" {
                    self.graphic_png_title = "SVG".to_string();
                } else {
                    self.graphic_png_title = "PNG".to_string();
                }
                self.png_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopyDescrips => {
                self.descrips_copy_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.descrips_for_spreadsheet.as_bytes());
                Command::perform(noop1(), Message::CompleteCopyDescrips)
            }

            Message::CompleteCopyDescrips(_) => {
                self.descrips_copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::OpenAlluvialReadsDoc => {
                self.alluvial_reads_doc_open = true;
                Command::none()
            }

            Message::CloseAlluvialReadsDoc => {
                self.alluvial_reads_doc_open = false;
                Command::none()
            }

            Message::SetSummaryScrollablePos(p) => {
                self.summary_scroll.snap_to(p);
                Command::none()
            }

            Message::SummarySnapshot => {
                self.summary_snapshot_start = Some(Instant::now());
                self.summary_snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteSummarySnapshot)
            }

            Message::CompleteSummarySnapshot(_) => {
                snapshot(&self.summary_snapshot_start);
                self.summary_snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopyAlluvialReadsTables => {
                self.alluvial_reads_tables_copy_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.alluvial_reads_tables_for_spreadsheet.as_bytes());
                Command::perform(noop1(), Message::CompleteCopyAlluvialReadsTables)
            }

            Message::CompleteCopyAlluvialReadsTables(_) => {
                self.alluvial_reads_tables_copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopyAlluvialTables => {
                self.alluvial_tables_copy_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.alluvial_tables_for_spreadsheet.as_bytes());
                Command::perform(noop1(), Message::CompleteCopyAlluvialTables)
            }

            Message::CompleteCopyAlluvialTables(_) => {
                self.alluvial_tables_copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::TooltipToggle => {
                self.tooltip_toggle_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                let pos = (TOOLTIP_POS.load(SeqCst) + 1) % 4;
                TOOLTIP_POS.store(pos, SeqCst);
                Command::perform(noop1(), Message::CompleteTooltipToggle)
            }

            Message::CompleteTooltipToggle(_) => {
                self.tooltip_toggle_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::GraphicSnapshot => {
                self.graphic_snapshot_start = Some(Instant::now());
                self.graphic_snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteGraphicSnapshot)
            }

            Message::CompleteGraphicSnapshot(_) => {
                snapshot(&self.graphic_snapshot_start);
                self.graphic_snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::GraphicOpen(_) => {
                self.graphic_mode = true;
                GRAPHIC_MODE.store(true, SeqCst);
                Command::none()
            }

            Message::GraphicClose => {
                self.graphic_mode = false;
                self.graphic_help_mode = true;
                self.graphic_help_title = "Help".to_string();
                GRAPHIC_MODE.store(false, SeqCst);
                self.graphic_png_title = "PNG".to_string();
                Command::none()
            }

            Message::ClonotypesCopy => {
                self.clonotypes_copy_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                let mut s = String::new();
                for line in self.output_value.lines() {
                    let t = line.replace(" ", "");
                    if t.len() > 0 {
                        s += &mut format!("{}\n", line);
                    }
                }
                copy_bytes_to_clipboard(&s.as_bytes());
                Command::perform(noop1(), Message::CompleteClonotypesCopy)
            }

            Message::CompleteClonotypesCopy(_) => {
                self.clonotypes_copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::SanityCheck => {
                self.sanity_check_start = Some(Instant::now());
                self.sanity_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteSanityCheck)
            }

            Message::CompleteSanityCheck(_) => {
                self.sanity_check();
                let used = elapsed(&self.sanity_check_start.unwrap());
                const MIN_SLEEP: f64 = 0.4;
                if used < MIN_SLEEP {
                    let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
                    thread::sleep(Duration::from_millis(ms));
                }
                self.sanity_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::Sleep(ms) => {
                thread::sleep(Duration::from_millis(ms));
                Command::none()
            }

            Message::CopyLastNarrative => {
                let index = self.h.narrative_history[(self.h.history_index - 2) as usize];
                let last = self.h.narrative_hist_uniq[index as usize].clone();
                let len = self.h.narrative_hist_uniq.len();
                self.h.narrative_hist_uniq.push(last);
                self.h.narrative_history[(self.h.history_index - 1) as usize] = len as u32;
                Command::none()
            }

            Message::Recompute => {
                let n = self.state_count();
                if n == 0 {
                    return Command::none();
                }
                let mut messages = Vec::<Message>::new();
                messages.push(Message::ConsoleClose);
                let k = self.hi();
                for _ in 0..k {
                    messages.push(Message::BackButtonPressed(Ok(())));
                }
                for i in 0..n {
                    messages.push(Message::SubmitButtonPressed(Ok(())));
                    messages.push(Message::CopyLastNarrative);
                    messages.push(Message::BackButtonPressed(Ok(())));
                    messages.push(Message::DelButtonPressed(Ok(())));
                    if i < n - 1 {
                        messages.push(Message::ForwardButtonPressed(Ok(())));
                        if i > 0 {
                            messages.push(Message::ForwardButtonPressed(Ok(())));
                        }
                    }
                }
                if n >= 2 {
                    for _ in 0..n - 2 {
                        messages.push(Message::BackButtonPressed(Ok(())));
                    }
                }
                messages.push(Message::ConsoleOpen);
                self.meta_pos = 0;
                self.this_meta = messages;
                META_TESTING.store(true, SeqCst);
                PSEUDO_META.store(true, SeqCst);
                Command::perform(noop0(), Message::Meta)
            }

            Message::Snapshot => {
                self.snapshot_start = Some(Instant::now());
                self.snapshot_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop0(), Message::CompleteSnapshot)
            }

            Message::CompleteSnapshot(_) => {
                snapshot(&self.snapshot_start);
                self.snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopySelectedMetrics => {
                self.copy_selected_metrics_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                let show = &self.metric_selected;
                copy_bytes_to_clipboard(
                    &expand_summary_as_csv(&self.summary_current(), &show).as_bytes(),
                );
                Command::perform(noop1(), Message::CompleteCopySelectedMetrics)
            }

            Message::CompleteCopySelectedMetrics(_) => {
                self.copy_selected_metrics_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CondenseMetrics => {
                self.metrics_condensed = !self.metrics_condensed;
                Command::none()
            }

            Message::MetricButton(index) => {
                self.metric_selected[index] = !self.metric_selected[index];
                Command::none()
            }

            Message::WaitCommand(_) => {
                while PROCESSING_REQUEST.load(SeqCst) {
                    thread::sleep(Duration::from_millis(10));
                }
                Command::perform(noop1(), Message::Meta)
            }

            Message::CopySummary => {
                self.copy_summary_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                let show;
                if !self.metrics_condensed {
                    show = vec![true; self.metric_selected.len()];
                } else {
                    show = self.metric_selected.clone();
                }
                copy_bytes_to_clipboard(
                    &expand_summary(&self.summary_current(), true, &show).as_bytes(),
                );
                Command::perform(noop1(), Message::CompleteCopySummary)
            }

            Message::CompleteCopySummary(_) => {
                self.copy_summary_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopyNarrative => {
                self.copy_narrative_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.narrative_current().as_bytes());
                Command::perform(noop1(), Message::CompleteCopyNarrative)
            }

            Message::CompleteCopyNarrative(_) => {
                self.copy_narrative_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopyArchiveNarrative(i) => {
                self.copy_archive_narrative_button_color[i] = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.archive_narrative[i].as_bytes());
                Command::perform(noop1(), Message::CompleteCopyArchiveNarrative)
            }

            Message::CompleteCopyArchiveNarrative(_) => {
                for i in 0..self.copy_archive_narrative_button_color.len() {
                    self.copy_archive_narrative_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
                }
                Command::none()
            }

            Message::CopyCookbookNarrative(i) => {
                self.copy_cookbook_narrative_button_color[i] = Color::from_rgb(1.0, 0.0, 0.0);
                copy_bytes_to_clipboard(&self.cookbook_narrative[i].as_bytes());
                Command::perform(noop1(), Message::CompleteCopyCookbookNarrative)
            }

            Message::CompleteCopyCookbookNarrative(_) => {
                for i in 0..self.copy_cookbook_narrative_button_color.len() {
                    self.copy_cookbook_narrative_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
                }
                Command::none()
            }

            Message::Snap(x) => {
                capture(&x, self.window_id);
                Command::none()
            }

            Message::SetName(x) => {
                self.save_name = x.to_string();
                Command::none()
            }

            Message::Meta(_) => do_meta(self),

            Message::CompleteMeta(_) => {
                capture(&self.save_name, self.window_id);
                if self.meta_pos == self.this_meta.len() {
                    if PSEUDO_META.load(SeqCst) {
                        PSEUDO_META.store(false, SeqCst);
                        META_TESTING.store(false, SeqCst);
                        return Command::none();
                    }
                    if PLAYBACK.load(SeqCst) {
                        println!("message history:\n");
                        let messages = compressed_message_history();
                        for i in 0..messages.len() {
                            println!("[{}] {}", i + 1, messages[i]);
                        }
                    }
                    println!("");
                    std::process::exit(0);
                }
                Command::perform(noop0(), Message::Meta)
            }

            Message::NullMeta(_) => Command::perform(noop0(), Message::Meta),

            Message::Narrative => {
                self.modified = true;
                let copy = get_clipboard_content();
                if copy.is_some() {
                    let copy = copy.unwrap();
                    self.narrative_value = copy.clone();
                    let len = self.h.narrative_hist_uniq.len();
                    self.h.narrative_hist_uniq.push(copy);
                    self.h.narrative_history[(self.h.history_index - 1) as usize] = len as u32;
                }
                Command::none()
            }

            Message::ArchiveNarrative(i) => {
                self.modified = true;
                let copy = get_clipboard_content();
                if copy.is_some() {
                    let copy = copy.unwrap();
                    self.archive_narrative[i] = copy.clone();
                    let filename = format!(
                        "{}/{}",
                        self.archive_dir.as_ref().unwrap(),
                        &self.archive_list[i]
                    );
                    let res = rewrite_narrative(&filename, &copy);
                    if res.is_err() {
                        xprintln!(
                            "\nSomething went wrong changing the narrative of\n{}\n\
                            Possibly the file has been corrupted.\n",
                            filename,
                        );
                        std::process::exit(1);
                    }
                }
                Command::none()
            }

            Message::OpenArchiveDoc => {
                self.archive_doc_open = true;
                Command::none()
            }

            Message::CloseArchiveDoc => {
                self.archive_doc_open = false;
                Command::none()
            }

            Message::DoShare(check_val) => do_share_button_pressed(self, check_val),

            Message::CompleteDoShare(_) => {
                self.do_share_complete = true;
                Command::none()
            }

            Message::Restore(check_val, index) => {
                if !self.just_restored && !self.delete_requested[index] {
                    self.restore_requested[index] = check_val;
                    self.restore_msg[index] =
                        "Restore scheduled!  Now click Dismiss or Save and dismiss at top."
                            .to_string();
                    self.just_restored = true;
                    self.modified = false;
                }
                Command::none()
            }

            Message::RestoreCookbook(check_val, index) => {
                if !self.just_restored {
                    self.restore_cookbook_requested[index] = check_val;
                    self.restore_cookbook_msg[index] =
                        "Restore scheduled!  Now click Dismiss or Save and dismiss at top."
                            .to_string();
                    self.just_restored = true;
                    self.modified = false;
                }
                Command::none()
            }

            Message::Save => {
                self.save_in_progress = true;
                self.save("");
                Command::perform(noop1(), Message::CompleteSave)
            }

            Message::SaveAs(x) => {
                self.save_in_progress = true;
                self.save_as(&x, "");
                Command::perform(noop1(), Message::CompleteSave)
            }

            Message::CompleteSave(_) => {
                self.save_in_progress = false;
                Command::none()
            }

            Message::Exit => {
                if true {
                    if self.sharing_enabled {
                        let share_bytes = unsafe { self.shares.align_to::<u8>().1.to_vec() };
                        let share_file = format!("{}/shares", self.visual);
                        std::fs::write(&share_file, &share_bytes).unwrap();
                    }
                    if self.save_on_exit {
                        let mut now = format!("{:?}", Local::now());
                        now = now.replace("T", "___");
                        now = now.before(".").to_string();
                        let filename = format!("{}/{}", self.archive_dir.as_ref().unwrap(), now);
                        let res = write_enclone_visual_history(&self.h, &filename);
                        if res.is_err() {
                            xprintln!(
                                "Was Unable to write history to the file {}, \
                                so Save on Exit failed.\n",
                                filename
                            );
                            std::process::exit(1);
                        }
                    }
                    if PLAYBACK.load(SeqCst) {
                        println!("message history:\n");
                        let messages = compressed_message_history();
                        for i in 0..messages.len() {
                            println!("[{}] {}", i + 1, messages[i]);
                        }
                    }
                    println!("");
                    std::process::exit(0);
                }
                Command::none()
            }

            Message::UserName(x, index) => {
                self.user_value[index] = x.to_string();
                Command::none()
            }

            Message::UserSelected(check_val, index) => {
                self.user_selected[index] = check_val;
                if check_val {
                    let user = &self.user_value[index];
                    USER_NAME.lock().unwrap().clear();
                    USER_NAME.lock().unwrap().push(user.to_string());
                    TESTING_USER_NAME.store(true, SeqCst);
                    loop {
                        thread::sleep(Duration::from_millis(10));
                        if !TESTING_USER_NAME.load(SeqCst) {
                            self.user_valid[index] = USER_NAME_VALID.load(SeqCst);
                            break;
                        }
                    }
                    if self.user_valid[index] {
                        if self.user_value[self.user_value.len() - 1].len() > 0 {
                            self.user.push(iced::text_input::State::new());
                            self.user_value.push(String::new());
                            self.user_selected.push(false);
                            self.user_valid.push(false);
                        }
                    }
                }
                Command::none()
            }

            Message::ArchiveShare(check_val, index) => do_archive_share(self, check_val, index),

            Message::NameChange(check_val) => {
                self.name_change_requested = check_val;
                if !check_val {
                    self.h.name_value = self.h.orig_name_value.clone();
                }
                Command::none()
            }

            Message::Name(x) => {
                self.h.name_value = x.to_string();
                Command::none()
            }

            Message::ArchiveNameChange(i) => {
                self.archive_name_change_button_color[i] = Color::from_rgb(1.0, 0.0, 0.0);
                if self.archive_name_value[i] != self.orig_archive_name[i] {
                    let filename = format!(
                        "{}/{}",
                        self.archive_dir.as_ref().unwrap(),
                        &self.archive_list[i]
                    );
                    let res = rewrite_name(&filename, &self.archive_name_value[i]);
                    if res.is_err() {
                        xprintln!(
                            "Something went wrong changing the name of\n{}\n\
                            Possibly the file has been corrupted.\n",
                            filename,
                        );
                        std::process::exit(1);
                    }
                    self.orig_archive_name[i] = self.archive_name_value[i].clone();
                }
                Command::perform(noop1(), Message::CompleteArchiveNameChange)
            }

            Message::CompleteArchiveNameChange(_) => {
                for i in 0..self.archive_name_change_button_color.len() {
                    self.archive_name_change_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
                }
                Command::none()
            }

            Message::ArchiveName(x, index) => {
                let mut y = String::new();
                for (i, char) in x.chars().enumerate() {
                    if i == 30 {
                        break;
                    }
                    y.push(char);
                }
                self.archive_name_value[index] = y;
                Command::none()
            }

            Message::DeleteArchiveEntry(check_val, index) => {
                if !self.just_restored {
                    self.delete_requested[index] = check_val;
                    if check_val {
                        self.expand_archive_entry[index] = false;
                        self.restore_msg[index] = "Will be deleted upon refresh or dismissal of \
                            this page.  Before then, you can change your mind \
                            and unclick!"
                            .to_string();
                    } else {
                        self.restore_msg[index].clear();
                    }
                }
                Command::none()
            }

            Message::ExpandArchiveEntry(check_val, index) => {
                if !self.delete_requested[index] {
                    self.expand_archive_entry[index] = check_val;
                }
                Command::none()
            }

            Message::ExpandCookbookEntry(check_val, index) => {
                self.expand_cookbook_entry[index] = check_val;
                Command::none()
            }

            Message::Resize(width, height) => {
                self.width = width;
                CURRENT_WIDTH.store(width as usize, SeqCst);
                self.height = height;
                CURRENT_HEIGHT.store(height as usize, SeqCst);
                Command::none()
            }

            Message::GroupClicked(_message) => {
                if GROUP_ID_CLICKED_ON.load(SeqCst) {
                    GROUP_ID_CLICKED_ON.store(false, SeqCst);
                    if TOOLTIP_TEXT.lock().unwrap().is_empty() {
                        // This case happens rarely, and we don't know why, but haven't
                        // studied the problem.
                        return Command::none();
                    }
                    self.modified = true;
                    let group_id = GROUP_ID.load(SeqCst);
                    self.input_value = format!("{}", group_id);
                    self.input1_value = format!("{}", group_id);
                    self.input2_value.clear();
                    for i in 0..self.inputn_value.len() {
                        self.inputn_value[i].clear();
                    }
                    let tt = TOOLTIP_TEXT.lock().unwrap()[0].clone();
                    copy_bytes_to_clipboard(&tt.as_bytes());
                    Command::perform(noop0(), Message::SubmitButtonPressed)
                } else {
                    Command::none()
                }
            }

            Message::SubmitButtonPressed(_) => do_submit_button_pressed(self),

            Message::DelButtonPressed(_) => do_del_button_pressed(self),

            Message::BackButtonPressed(_) => {
                self.h.history_index -= 1;
                self.update_to_current();
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop0(), Message::Capture)
                }
            }

            Message::ForwardButtonPressed(_) => {
                self.h.history_index += 1;
                self.update_to_current();
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
                        for i in 0..self.inputn_value.len() {
                            self.inputn_value[i].clear();
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

            Message::ClonotypesOpen(_) => {
                self.clonotypes_mode = true;
                Command::none()
            }

            Message::ClonotypesClose => {
                self.clonotypes_mode = false;
                Command::none()
            }

            Message::SummaryOpen(_) => {
                self.summary_mode = true;
                let summaryx = unpack_summary(&self.summary_value);
                self.metric_selected = summaryx.metric_selected.clone();
                self.metrics_condensed = summaryx.metrics_condensed;
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

            Message::SummaryClose(_) => {
                self.summary_mode = false;
                let mut summaryx = unpack_summary(&self.summary_value);
                summaryx.metric_selected = self.metric_selected.clone();
                summaryx.metrics_condensed = self.metrics_condensed;
                self.summary_value = summaryx.pack();
                self.h.summary_hist_uniq
                    [self.h.summary_history[self.h.history_index as usize - 1] as usize] =
                    self.summary_value.clone();
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

            Message::ArchiveRefresh => {
                self.share_start = Some(Instant::now());
                self.archive_refresh_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop(), Message::ArchiveRefreshComplete)
            }

            Message::ArchiveRefreshComplete(_) => do_archive_refresh_complete(self),

            Message::ArchiveClose => do_archive_close(self, false),

            Message::ArchiveSaveClose => do_archive_close(self, true),

            Message::ArchiveOpen(_) => do_archive_open(self),

            Message::SaveOnExit => {
                self.save_on_exit = !self.save_on_exit;
                Command::none()
            }

            Message::InputChanged1(ref value) => {
                self.modified = true;
                self.input1_value = value.to_string();
                self.input_value = self.input1_value.clone();
                if self.input1_value.len() > 0 && self.input2_value.len() > 0 {
                    self.input_value += " ";
                }
                self.input_value += &mut self.input2_value.clone();
                Command::none()
            }

            Message::InputChanged2(ref value) => {
                self.modified = true;
                self.input2_value = value.to_string();
                self.input_value = self.input1_value.clone();
                if self.input1_value.len() > 0 && self.input2_value.len() > 0 {
                    self.input_value += " ";
                }
                self.input_value += &mut self.input2_value.clone();
                Command::none()
            }

            Message::InputChangedN(ref value, i) => {
                self.modified = true;
                self.inputn_value[i] = value.to_string();
                let mut values = Vec::<String>::new();
                if self.input1_value.len() > 0 {
                    values.push(self.input1_value.clone());
                }
                if self.input2_value.len() > 0 {
                    values.push(self.input2_value.clone());
                }
                for j in 0..self.inputn_value.len() {
                    if self.inputn_value[j].len() > 0 {
                        values.push(self.inputn_value[j].clone());
                    }
                }
                self.input_value = format!("{}", values.iter().format(" "));
                Command::none()
            }

            Message::ClearButtonPressed => {
                self.modified = true;
                self.input1_value.clear();
                self.input2_value.clear();
                for j in 0..self.inputn_value.len() {
                    self.inputn_value[j].clear();
                }
                Command::none()
            }

            Message::ComputationDone(_) => do_computation_done(self),

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
                        // note not showing inputn
                    }
                }
                let count = COUNT.load(SeqCst);
                if count >= 1 {
                    capture(&TESTS[count - 1].2, self.window_id);
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
                let mut width = 2000;
                let copy = get_clipboard_content();
                if copy.is_some() {
                    let copy = copy.unwrap();
                    if copy.parse::<usize>().is_ok() {
                        let w = copy.force_usize();
                        if w >= 1000 && w <= 4000 {
                            width = w as u32;
                        }
                    }
                }
                let png = convert_svg_to_png(&self.svg_value.as_bytes(), width);
                copy_png_bytes_to_clipboard(&png);
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
