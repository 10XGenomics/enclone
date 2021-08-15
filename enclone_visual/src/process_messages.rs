// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::copy_bytes_to_clipboard;
use crate::history::*;
use crate::messages::*;
use crate::proc1::*;
use crate::share::*;
use crate::testsuite::TESTS;
use crate::*;
use chrono::prelude::*;
use clipboard::{ClipboardContext, ClipboardProvider};
use flate2::read::GzDecoder;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use io_utils::*;
use std::fs::File;
use std::io::Read;
use std::time::{Duration, Instant};
use vector_utils::*;

impl EncloneVisual {
    pub fn process_message(&mut self, message: Message) -> Command<Message> {
        MESSAGE_HISTORY
            .lock()
            .unwrap()
            .push(format!("{:?}", message));
        match message {
            Message::Narrative => {
                self.modified = true;
                let ctx: Result<ClipboardContext, _> = ClipboardProvider::new();
                if ctx.is_err() {
                    xprintln!("\nSomething went wrong accessing clipboard.");
                    xprintln!("This is weird so please ask for help.");
                    std::process::exit(1);
                }
                let mut ctx = ctx.unwrap();
                let copy = ctx.get_contents();
                if copy.is_err() {
                    xprintln!("\nSomething went wrong copying from clipboard.");
                    xprintln!("This is weird so please ask for help.");
                    std::process::exit(1);
                }
                let copy = format!("{}", ctx.get_contents().unwrap());
                self.narrative_value = copy.clone();
                let len = self.h.narrative_hist_uniq.len();
                self.h.narrative_hist_uniq.push(copy);
                self.h.narrative_history[(self.h.history_index - 1) as usize] = len as u32;
                Command::none()
            }

            Message::ArchiveNarrative(i) => {
                self.modified = true;
                let ctx: Result<ClipboardContext, _> = ClipboardProvider::new();
                if ctx.is_err() {
                    xprintln!("\nSomething went wrong accessing clipboard.");
                    xprintln!("This is weird so please ask for help.");
                    std::process::exit(1);
                }
                let mut ctx = ctx.unwrap();
                let copy = ctx.get_contents();
                if copy.is_err() {
                    xprintln!("\nSomething went wrong copying from clipboard.");
                    xprintln!("This is weird so please ask for help.");
                    std::process::exit(1);
                }
                let copy = format!("{}", ctx.get_contents().unwrap());
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

            Message::DoShare(check_val) => {
                self.do_share = check_val;
                if !check_val {
                    self.do_share_complete = false;
                } else {
                    let mut recipients = Vec::<String>::new();
                    for i in 0..self.user_value.len() {
                        if self.user_valid[i] {
                            recipients.push(self.user_value[i].clone());
                        }
                    }
                    let mut index = 0;
                    for i in 0..self.archive_share_requested.len() {
                        if self.archive_share_requested[i] {
                            index = i;
                        }
                    }
                    let path = format!(
                        "{}/{}",
                        self.archive_dir.as_ref().unwrap(),
                        self.archive_list[index]
                    );
                    if !path_exists(&path) {
                        xprintln!("could not find path for archive file\n");
                        std::process::exit(1);
                    }
                    let mut content = Vec::<u8>::new();
                    let f = File::open(&path);
                    if f.is_err() {
                        xprintln!("could not open archive file\n");
                        std::process::exit(1);
                    }
                    let mut f = f.unwrap();
                    let res = f.read_to_end(&mut content);
                    if res.is_err() {
                        xprintln!("could not read archive file\n");
                        std::process::exit(1);
                    }
                    SHARE_CONTENT.lock().unwrap().clear();
                    SHARE_CONTENT.lock().unwrap().push(content);
                    SHARE_RECIPIENTS.lock().unwrap().clear();
                    let days = Utc::now().num_days_from_ce();
                    for i in 0..recipients.len() {
                        SHARE_RECIPIENTS.lock().unwrap().push(recipients[i].clone());
                        let mut user_name = [0 as u8; 32];
                        for j in 0..recipients[i].len() {
                            user_name[j] = recipients[i].as_bytes()[j];
                        }
                        self.shares.push(Share {
                            days_since_ce: days,
                            user_id: user_name,
                        });
                    }
                    SENDING_SHARE.store(true, SeqCst);
                }
                Command::perform(compute_share(), Message::CompleteDoShare)
            }

            Message::CompleteDoShare(_) => {
                self.do_share_complete = true;
                Command::none()
            }

            Message::Restore(check_val, index) => {
                if !self.just_restored && !self.delete_requested[index] {
                    let mut index = index;
                    self.restore_requested[index] = check_val;
                    if self.modified {
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
                        self.restore_msg[index] = "Oh dear, restoration failed.".to_string();
                    }
                }
                Command::none()
            }

            Message::Save => {
                self.save_in_progress = true;
                self.save();
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
                        let dir;
                        if VISUAL_HISTORY_DIR.lock().unwrap().len() > 0 {
                            dir = VISUAL_HISTORY_DIR.lock().unwrap()[0].clone();
                        } else {
                            dir = format!("{}/history", self.visual);
                        }
                        let mut now = format!("{:?}", Local::now());
                        now = now.replace("T", "___");
                        now = now.before(".").to_string();
                        let filename = format!("{}/{}", dir, now);
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

            Message::ArchiveShare(check_val, index) => {
                if !check_val {
                    self.do_share = false;
                    self.do_share_complete = false;
                }
                let mut already_sharing = false;
                for i in 0..self.archive_share_requested.len() {
                    if i != index && self.archive_share_requested[i] {
                        already_sharing = true;
                    }
                }
                if !already_sharing {
                    self.archive_share_requested[index] = check_val;
                    if !check_val {
                        self.user.clear();
                        self.user_value.clear();
                        self.user_selected.clear();
                        self.user_valid.clear();
                    } else {
                        let mut names = Vec::<String>::new();
                        for i in 0..self.shares.len() {
                            let mut j = 0;
                            while j < self.shares[i].user_id.len() {
                                if self.shares[i].user_id[j] == 0 {
                                    break;
                                }
                                j += 1;
                            }
                            names.push(stringme(&self.shares[i].user_id[0..j]));
                        }
                        names.sort();
                        let mut freq = Vec::<(u32, String)>::new();
                        make_freq(&names, &mut freq);
                        const MAX_USERS_TO_SHOW: usize = 6;
                        let show = std::cmp::min(MAX_USERS_TO_SHOW, freq.len());
                        for i in 0..show {
                            self.user.push(iced::text_input::State::new());
                            self.user_value.push(freq[i].1.clone());
                            self.user_selected.push(false);
                            self.user_valid.push(false);
                        }
                        self.user.push(iced::text_input::State::new());
                        self.user_value.push(String::new());
                        self.user_selected.push(false);
                        self.user_valid.push(false);
                    }
                }
                Command::none()
            }

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
                if !self.delete_requested[index] {
                    self.delete_requested[index] = check_val;
                    self.expand_archive_entry[index] = false;
                    let filename = format!(
                        "{}/{}",
                        self.archive_dir.as_ref().unwrap(),
                        self.archive_list[index]
                    );
                    if path_exists(&filename) {
                        std::fs::remove_file(&filename).unwrap();
                    }
                    self.restore_msg[index] = "Deleted.".to_string();
                }
                Command::none()
            }

            Message::ExpandArchiveEntry(check_val, index) => {
                if !self.delete_requested[index] {
                    self.expand_archive_entry[index] = check_val;
                }
                Command::none()
            }

            Message::Resize(width, height) => {
                self.width = width;
                self.height = height;
                Command::none()
            }

            Message::GroupClicked(_message) => {
                self.modified = true;
                let group_id = GROUP_ID.load(SeqCst);
                self.input_value = format!("{}", group_id);
                self.input1_value = format!("{}", group_id);
                self.input2_value.clear();
                GROUP_ID_CLICKED_ON.store(false, SeqCst);
                Command::perform(noop0(), Message::SubmitButtonPressed)
            }

            Message::SubmitButtonPressed(_) => submit_button_pressed(self),

            Message::DelButtonPressed(_) => {
                self.modified = true;
                let h = self.h.history_index - 1;
                self.h.svg_history.remove(h as usize);
                self.h.summary_history.remove(h as usize);
                self.h.input1_history.remove(h as usize);
                self.h.input2_history.remove(h as usize);
                self.h.narrative_history.remove(h as usize);
                self.h.translated_input_history.remove(h as usize);
                self.h.displayed_tables_history.remove(h as usize);
                self.h.table_comp_history.remove(h as usize);
                self.h.last_widths_history.remove(h as usize);
                self.h.is_blank.remove(h as usize);
                self.h.descrip_history.remove(h as usize);
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
                    self.descrip_value.clear();
                    self.translated_input_value.clear();
                    self.current_tables.clear();
                } else {
                    if h > 0 {
                        self.h.history_index -= 1;
                    }
                    self.update_to_current();
                }
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    Command::perform(noop0(), Message::Capture)
                }
            }

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

            Message::UpdateShares => {
                self.share_start = Some(Instant::now());
                self.receive_shares_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                Command::perform(noop(), Message::UpdateSharesComplete)
            }

            Message::UpdateSharesComplete(_) => {
                update_shares(self);
                let n = self.archive_name.len();
                for i in 0..n {
                    // This is a dorky way of causing loading of command lists, etc. from disk
                    // occurs just once per session, and only if the archive button is pushed.
                    if self.archived_command_list[i].is_none() {
                        let x = &self.archive_list[i];
                        let path = format!("{}/{}", self.archive_dir.as_ref().unwrap(), x);
                        let (command_list, name, origin, narrative) = read_metadata(&path).unwrap();
                        self.archived_command_list[i] = Some(command_list);
                        self.archive_name_value[i] = name;
                        self.archive_origin[i] = origin;
                        self.archive_narrative[i] = narrative;
                    }
                }
                self.orig_archive_name = self.archive_name_value.clone();
                self.h.orig_name_value = self.h.name_value.clone();
                self.receive_shares_button_color = Color::from_rgb(0.0, 0.0, 0.0);

                // Sleep so that total time for updating of shares is at least 0.4 seconds.  This
                // keeps the "Receive shares" button red for at least that amount of time.

                const MIN_SLEEP: f64 = 0.4;
                let used = elapsed(&self.share_start.unwrap());
                if used < MIN_SLEEP {
                    let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
                    thread::sleep(Duration::from_millis(ms));
                }

                Command::none()
            }

            Message::ArchiveOpen(_) => {
                self.archive_mode = true;
                if self.sharing_enabled {
                    update_shares(self);
                }
                let n = self.archive_name.len();
                for i in 0..n {
                    self.archive_name_change_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
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
                for i in 0..self.restore_msg.len() {
                    self.restore_msg[i].clear();
                    self.restore_requested[i] = false;
                    if self.delete_requested[i] {
                        self.deleted[i] = true;
                    }
                }
                self.just_restored = false;
                Command::none()
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
                self.narrative_value.clear();
                let len = self.h.narrative_hist_uniq.len();
                if len > 0 && self.h.narrative_hist_uniq[len - 1] == self.narrative_value {
                    self.h
                        .narrative_history
                        .insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.narrative_history.insert(hi as usize, len as u32);
                    self.h
                        .narrative_hist_uniq
                        .push(self.narrative_value.clone());
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
                let len = self.h.descrip_hist_uniq.len();
                if len > 0 && self.h.descrip_hist_uniq[len - 1] == self.descrip_value {
                    self.h.descrip_history.insert(hi as usize, (len - 1) as u32);
                } else {
                    self.h.descrip_history.insert(hi as usize, len as u32);
                    self.h.descrip_hist_uniq.push(self.descrip_value.clone());
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
                    test_evh_read_write(&self.h, "/tmp/evh_test");
                    std::fs::remove_file("/tmp/evh_test").unwrap();
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
