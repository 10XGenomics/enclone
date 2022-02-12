// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::gui_structures::ComputeState::WaitingForRequest;
use crate::history::*;
use crate::messages::*;
use crate::share::*;
use crate::*;
use iced::{Color, Command};
use io_utils::*;
use std::time::Duration;
use vector_utils::*;

pub fn do_archive_close(slf: &mut EncloneVisual, save: bool) -> Command<Message> {
    let mut index = None;
    for i in 0..slf.restore_requested.len() {
        if slf.restore_requested[i] {
            index = Some(i);
        }
    }
    if index.is_some() {
        let mut index = index.unwrap();
        if save {
            slf.save("(saved upon restore)");
            index += 1;
        }
        let filename = format!(
            "{}/{}",
            slf.archive_dir.as_ref().unwrap(),
            slf.archive_list[index]
        );
        let res = read_enclone_visual_history(&filename);
        if res.is_ok() {
            slf.h = res.unwrap();
            // Ignore history index and instead rewind.
            if slf.h.history_index > 1 {
                slf.h.history_index = 1;
            }
            slf.update_to_current();
        } else {
            slf.restore_msg[index] = format!(
                "Oh dear, restoration of the file {} \
                failed.",
                filename
            );
        }
    }
    let mut index = None;
    for i in 0..slf.restore_cookbook_requested.len() {
        if slf.restore_cookbook_requested[i] {
            index = Some(i);
        }
    }
    if index.is_some() {
        let index = index.unwrap();
        if save {
            slf.save("(saved upon restore)");
        }
        let res = EncloneVisualHistory::restore_from_bytes(&slf.cookbooks[index]);
        slf.h = res.unwrap();
        // Ignore history index and instead rewind.
        if slf.h.history_index > 1 {
            slf.h.history_index = 1;
        }
        slf.update_to_current();
    }
    for i in 0..slf.archive_name_value.len() {
        slf.archive_name_value[i] = slf.orig_archive_name[i].clone();
    }
    slf.archive_mode = false;
    slf.do_share = false;
    slf.do_share_complete = false;
    slf.user.clear();
    slf.user_value.clear();
    slf.user_selected.clear();
    slf.user_valid.clear();
    for i in 0..slf.archive_share_requested.len() {
        slf.archive_share_requested[i] = false;
    }
    for i in 0..slf.expand_archive_entry.len() {
        slf.expand_archive_entry[i] = false;
    }
    for i in 0..slf.cookbooks.len() {
        slf.expand_cookbook_entry[i] = false;
        slf.restore_cookbook_requested[i] = false;
    }
    for i in 0..slf.restore_msg.len() {
        slf.restore_msg[i].clear();
        slf.restore_requested[i] = false;
        if slf.delete_requested[i] {
            let filename = format!(
                "{}/{}",
                slf.archive_dir.as_ref().unwrap(),
                slf.archive_list[i]
            );
            if path_exists(&filename) {
                std::fs::remove_file(&filename).unwrap();
            }
            slf.deleted[i] = true;
        }
    }
    for i in 0..slf.restore_cookbook_msg.len() {
        slf.restore_cookbook_msg[i].clear();
    }
    slf.just_restored = false;
    Command::none()
}

pub fn do_meta(slf: &mut EncloneVisual) -> Command<Message> {
    if slf.meta_pos == slf.this_meta.len() {
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
    for i in slf.meta_pos..slf.this_meta.len() {
        if i == 0 {
            slf.window_id = get_window_id();
        }

        match slf.this_meta[i] {
            Message::SubmitButtonPressed(_) => {
                slf.meta_pos = i + 1;
                submit = true;
                break;
            }
            Message::WaitCommand(_) => {
                slf.meta_pos = i + 1;
                wait = true;
                break;
            }
            _ => {}
        }

        // slf.update(slf.this_meta[i].clone(), clipboard);
        slf.update(slf.this_meta[i].clone());
        match slf.this_meta[i] {
            Message::SetName(_) => {
                slf.meta_pos = i + 1;
                done = true;
                break;
            }
            Message::WaitCommand(_) => {
                slf.meta_pos = i + 1;
                null = true;
                break;
            }
            _ => {}
        }
        if i == slf.this_meta.len() - 1 {
            slf.meta_pos = i + 1;
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

pub fn do_del_button_pressed(slf: &mut EncloneVisual) -> Command<Message> {
    if slf.compute_state != WaitingForRequest {
        Command::none()
    } else {
        slf.modified = true;
        let h = slf.h.history_index - 1;
        slf.h.svg_history.remove(h as usize);
        slf.h.summary_history.remove(h as usize);
        slf.h.input1_history.remove(h as usize);
        slf.h.input2_history.remove(h as usize);
        slf.h.inputn_history.remove(h as usize);
        slf.h.narrative_history.remove(h as usize);
        slf.h.translated_input_history.remove(h as usize);
        slf.h.displayed_tables_history.remove(h as usize);
        slf.h.table_comp_history.remove(h as usize);
        slf.h.last_widths_history.remove(h as usize);
        slf.h.is_blank.remove(h as usize);
        slf.h.descrip_history.remove(h as usize);
        if slf.state_count() == 0 {
            slf.h.history_index -= 1;
            slf.input1_value.clear();
            slf.input2_value.clear();
            slf.inputn_value.clear();
            slf.svg_value.clear();
            slf.png_value.clear();
            slf.submit_button_text.clear();
            slf.summary_value.clear();
            slf.output_value.clear();
            slf.table_comp_value.clear();
            slf.last_widths_value.clear();
            slf.descrip_value.clear();
            slf.translated_input_value.clear();
            slf.current_tables.clear();
        } else {
            if h > 0 {
                slf.h.history_index -= 1;
            }
            slf.update_to_current();
        }
        if !TEST_MODE.load(SeqCst) {
            Command::none()
        } else {
            Command::perform(noop0(), Message::Capture)
        }
    }
}

pub fn do_archive_share(
    slf: &mut EncloneVisual,
    check_val: bool,
    index: usize,
) -> Command<Message> {
    if !check_val {
        slf.do_share = false;
        slf.do_share_complete = false;
    }
    let mut already_sharing = false;
    for i in 0..slf.archive_share_requested.len() {
        if i != index && slf.archive_share_requested[i] {
            already_sharing = true;
        }
    }
    if !already_sharing {
        slf.archive_share_requested[index] = check_val;
        if !check_val {
            slf.user.clear();
            slf.user_value.clear();
            slf.user_selected.clear();
            slf.user_valid.clear();
        } else {
            let mut names = Vec::<String>::new();
            for i in 0..slf.shares.len() {
                let mut j = 0;
                while j < slf.shares[i].user_id.len() {
                    if slf.shares[i].user_id[j] == 0 {
                        break;
                    }
                    j += 1;
                }
                names.push(stringme(&slf.shares[i].user_id[0..j]));
            }
            names.sort();
            let mut freq = Vec::<(u32, String)>::new();
            make_freq(&names, &mut freq);
            const MAX_USERS_TO_SHOW: usize = 6;
            let show = std::cmp::min(MAX_USERS_TO_SHOW, freq.len());
            for i in 0..show {
                slf.user.push(iced::text_input::State::new());
                slf.user_value.push(freq[i].1.clone());
                slf.user_selected.push(false);
                slf.user_valid.push(false);
            }
            slf.user.push(iced::text_input::State::new());
            slf.user_value.push(String::new());
            slf.user_selected.push(false);
            slf.user_valid.push(false);
        }
    }
    Command::none()
}

pub fn do_archive_refresh_complete(slf: &mut EncloneVisual) -> Command<Message> {
    update_shares(slf);
    let n = slf.archive_name.len();
    for i in 0..n {
        if slf.delete_requested[i] {
            let filename = format!(
                "{}/{}",
                slf.archive_dir.as_ref().unwrap(),
                slf.archive_list[i]
            );
            if path_exists(&filename) {
                std::fs::remove_file(&filename).unwrap();
            }
        }
    }
    slf.orig_archive_name = slf.archive_name_value.clone();
    slf.h.orig_name_value = slf.h.name_value.clone();
    slf.archive_refresh_button_color = Color::from_rgb(0.0, 0.0, 0.0);

    // Sleep so that total time for updating of shares is at least 0.4 seconds.  This
    // keeps the Refresh button red for at least that amount of time.

    const MIN_SLEEP: f64 = 0.4;
    let used = elapsed(&slf.share_start.unwrap());
    if used < MIN_SLEEP {
        let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
        thread::sleep(Duration::from_millis(ms));
    }
    for i in 0..n {
        if slf.archived_command_list[i].is_none() {
            let x = &slf.archive_list[i];
            let path = format!("{}/{}", slf.archive_dir.as_ref().unwrap(), x);
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
            slf.archived_command_list[i] = Some(command_list);
            slf.archive_name_value[i] = name;
            slf.archive_origin[i] = origin;
            slf.archive_narrative[i] = narrative;
        }
    }
    slf.do_share = false;
    slf.do_share_complete = false;
    slf.user.clear();
    slf.user_value.clear();
    slf.user_selected.clear();
    slf.user_valid.clear();
    for i in 0..slf.archive_share_requested.len() {
        slf.archive_share_requested[i] = false;
    }
    for i in 0..slf.expand_archive_entry.len() {
        slf.expand_archive_entry[i] = false;
    }
    for i in 0..slf.cookbooks.len() {
        slf.expand_cookbook_entry[i] = false;
        slf.restore_cookbook_requested[i] = false;
    }
    for i in 0..slf.restore_msg.len() {
        slf.restore_msg[i].clear();
        slf.restore_requested[i] = false;
        if slf.delete_requested[i] {
            slf.deleted[i] = true;
        }
    }
    for i in 0..slf.restore_cookbook_msg.len() {
        slf.restore_cookbook_msg[i].clear();
    }
    slf.just_restored = false;
    for i in 0..n {
        slf.archive_name_value[i] = slf.orig_archive_name[i].clone();
        slf.archive_name_change_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
        slf.copy_archive_narrative_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
    }
    if !TEST_MODE.load(SeqCst) {
        Command::none()
    } else {
        Command::perform(noop1(), Message::Capture)
    }
}

pub fn do_archive_open(slf: &mut EncloneVisual) -> Command<Message> {
    slf.archive_mode = true;
    update_shares(slf);
    let n = slf.archive_name.len();
    for i in 0..n {
        slf.archive_name_change_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
        slf.copy_archive_narrative_button_color[i] = Color::from_rgb(0.0, 0.0, 0.0);
        // This is a dorky way of causing loading of command lists, etc. from disk
        // occurs just once per session, and only if the archive button is pushed.
        if slf.archived_command_list[i].is_none() {
            let x = &slf.archive_list[i];
            let path = format!("{}/{}", slf.archive_dir.as_ref().unwrap(), x);
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
            slf.archived_command_list[i] = Some(command_list);
            slf.archive_name_value[i] = name;
            slf.archive_origin[i] = origin;
            slf.archive_narrative[i] = narrative;
        }
    }
    slf.orig_archive_name = slf.archive_name_value.clone();
    slf.h.orig_name_value = slf.h.name_value.clone();
    if !TEST_MODE.load(SeqCst) {
        Command::none()
    } else {
        Command::perform(noop1(), Message::Capture)
    }
}
