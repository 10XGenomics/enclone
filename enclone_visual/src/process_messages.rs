// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_clipboard;
use crate::history::*;
use crate::messages::*;
use crate::proc1::*;
use crate::share::*;
use crate::summary::*;
use crate::testsuite::TESTS;
use crate::*;
use chrono::prelude::*;
use iced::{Clipboard, Color, Command};
use io_utils::*;
use std::fs::{remove_file, File};
use std::io::Read;
use std::time::{Duration, Instant};

impl EncloneVisual {
    pub fn process_message(
        &mut self,
        message: Message,
        clipboard: &mut Clipboard,
    ) -> Command<Message> {
        MESSAGE_HISTORY
            .lock()
            .unwrap()
            .push(format!("{:?}", message));
        match message {
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
                let filename = "/tmp/enclone_visual_snapshot.png";
                capture_as_file(&filename, get_window_id());
                let mut bytes = Vec::<u8>::new();
                {
                    let mut f = File::open(&filename).unwrap();
                    f.read_to_end(&mut bytes).unwrap();
                }
                remove_file(&filename).unwrap();
                copy_png_bytes_to_clipboard(&bytes);
                const MIN_SLEEP: f64 = 0.4;
                let used = elapsed(&self.snapshot_start.unwrap());
                if used < MIN_SLEEP {
                    let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
                    thread::sleep(Duration::from_millis(ms));
                }
                self.snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::CopySelectedMetrics => {
                self.copy_selected_metrics_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                let show = &self.metric_selected;
                copy_bytes_to_clipboard(
                    &expand_summary(&self.summary_current(), false, &show).as_bytes(),
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

            Message::Meta(_) => {
                if self.meta_pos == self.this_meta.len() {
                    if PSEUDO_META.load(SeqCst) {
                        PSEUDO_META.store(false, SeqCst);
                        META_TESTING.store(false, SeqCst);
                        return Command::none();
                    }
                    std::process::exit(0);
                }
                let mut done = false;
                let mut null = false;
                let mut submit = false;
                let mut wait = false;
                for i in self.meta_pos..self.this_meta.len() {
                    if i == 0 {
                        self.window_id = get_window_id();
                    }

                    match self.this_meta[i] {
                        Message::SubmitButtonPressed(_) => {
                            self.meta_pos = i + 1;
                            submit = true;
                            break;
                        }
                        Message::WaitCommand(_) => {
                            self.meta_pos = i + 1;
                            wait = true;
                            break;
                        }
                        _ => {}
                    }

                    self.update(self.this_meta[i].clone(), clipboard);
                    match self.this_meta[i] {
                        Message::SetName(_) => {
                            self.meta_pos = i + 1;
                            done = true;
                            break;
                        }
                        Message::WaitCommand(_) => {
                            self.meta_pos = i + 1;
                            null = true;
                            break;
                        }
                        _ => {}
                    }
                    if i == self.this_meta.len() - 1 {
                        self.meta_pos = i + 1;
                        done = true;
                    }
                }
                if submit {
                    Command::perform(noop0(), Message::SubmitButtonPressed)
                } else if wait {
                    Command::perform(noop0(), Message::WaitCommand)
                } else if null {
                    Command::perform(noop0(), Message::NullMeta)
                } else if done {
                    Command::perform(noop0(), Message::CompleteMeta)
                } else {
                    Command::none()
                }
            }

            Message::CompleteMeta(_) => {
                capture(&self.save_name, self.window_id);
                if self.meta_pos == self.this_meta.len() {
                    if PSEUDO_META.load(SeqCst) {
                        PSEUDO_META.store(false, SeqCst);
                        META_TESTING.store(false, SeqCst);
                        return Command::none();
                    }
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
                    let mut index = index;
                    self.restore_requested[index] = check_val;
                    if self.modified && self.state_count() > 0 {
                        self.save();
                        index += 1;
                    }
                    let filename = format!(
                        "{}/{}",
                        self.archive_dir.as_ref().unwrap(),
                        self.archive_list[index]
                    );
                    let res = read_enclone_visual_history(&filename);
                    if res.is_ok() {
                        self.h = res.unwrap();
                        // Ignore history index and instead rewind.
                        if self.h.history_index > 1 {
                            self.h.history_index = 1;
                        }
                        self.update_to_current();
                        self.restore_msg[index] =
                            "Restored!  Now click Dismiss at top.".to_string();
                        self.just_restored = true;
                        self.modified = false;
                    } else {
                        self.restore_msg[index] = format!(
                            "Oh dear, restoration of the file {} \
                            failed.",
                            filename
                        );
                    }
                }
                Command::none()
            }

            Message::RestoreCookbook(check_val, index) => {
                if !self.just_restored {
                    self.restore_cookbook_requested[index] = check_val;
                    if self.modified {
                        self.save();
                    }
                    let res = EncloneVisualHistory::restore_from_bytes(&self.cookbooks[index]);
                    self.h = res.unwrap();
                    // Ignore history index and instead rewind.
                    if self.h.history_index > 1 {
                        self.h.history_index = 1;
                    }
                    self.update_to_current();
                    self.restore_cookbook_msg[index] =
                        "Restored!  Now click Dismiss at top.".to_string();
                    self.just_restored = true;
                    self.modified = false;
                }
                Command::none()
            }

            Message::Save => {
                self.save_in_progress = true;
                self.save();
                Command::perform(noop1(), Message::CompleteSave)
            }

            Message::SaveAs(x) => {
                self.save_in_progress = true;
                self.save_as(&x);
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
                self.h.summary_hist_uniq[self.h.history_index as usize - 1] =
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

            Message::ArchiveClose => {
                for i in 0..self.archive_name_value.len() {
                    self.archive_name_value[i] = self.orig_archive_name[i].clone();
                }
                self.archive_mode = false;
                self.do_share = false;
                self.do_share_complete = false;
                self.user.clear();
                self.user_value.clear();
                self.user_selected.clear();
                self.user_valid.clear();
                for i in 0..self.archive_share_requested.len() {
                    self.archive_share_requested[i] = false;
                }
                for i in 0..self.expand_archive_entry.len() {
                    self.expand_archive_entry[i] = false;
                }
                for i in 0..self.cookbooks.len() {
                    self.expand_cookbook_entry[i] = false;
                    self.restore_cookbook_requested[i] = false;
                }
                for i in 0..self.restore_msg.len() {
                    self.restore_msg[i].clear();
                    self.restore_requested[i] = false;
                    if self.delete_requested[i] {
                        let filename = format!(
                            "{}/{}",
                            self.archive_dir.as_ref().unwrap(),
                            self.archive_list[i]
                        );
                        if path_exists(&filename) {
                            std::fs::remove_file(&filename).unwrap();
                        }
                        self.deleted[i] = true;
                    }
                }
                for i in 0..self.restore_cookbook_msg.len() {
                    self.restore_cookbook_msg[i].clear();
                }
                self.just_restored = false;
                Command::none()
            }

            Message::ArchiveOpen(_) => {
                self.archive_mode = true;
                update_shares(self);
                let n = self.archive_name.len();
                for i in 0..n {
                    self.archive_name_change_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
                    self.copy_archive_narrative_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
                    // This is a dorky way of causing loading of command lists, etc. from disk
                    // occurs just once per session, and only if the archive button is pushed.
                    if self.archived_command_list[i].is_none() {
                        let x = &self.archive_list[i];
                        let path = format!("{}/{}", self.archive_dir.as_ref().unwrap(), x);
                        let res = read_metadata(&path);
                        if res.is_err() {
                            panic!(
                                "Unable to read the history file at\n{}\n\
                                This could either be a bug in enclone or it could be that \
                                the file is corrupted.\n",
                                path,
                            );
                        }
                        let (command_list, name, origin, narrative) = res.unwrap();
                        self.archived_command_list[i] = Some(command_list);
                        self.archive_name_value[i] = name;
                        self.archive_origin[i] = origin;
                        self.archive_narrative[i] = narrative;
                    }
                }
                self.orig_archive_name = self.archive_name_value.clone();
                self.h.orig_name_value = self.h.name_value.clone();
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop1(), Message::Capture)
                }
            }

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

            Message::ClearButtonPressed => {
                self.modified = true;
                self.input1_value.clear();
                self.input2_value.clear();
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
