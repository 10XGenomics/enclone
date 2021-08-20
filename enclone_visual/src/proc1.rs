// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::messages::*;
use crate::*;
use chrono::prelude::*;
use enclone_core::combine_group_pics::*;
use gui_structures::ComputeState::*;
use iced::Command;
use io_utils::*;
use itertools::Itertools;
use std::fs::File;
use std::io::Read;
use std::time::{Duration, Instant};
use vector_utils::*;

pub fn do_share_button_pressed(slf: &mut EncloneVisual, check_val: bool) -> Command<Message> {
    slf.do_share = check_val;
    if !check_val {
        slf.do_share_complete = false;
    } else {
        let mut recipients = Vec::<String>::new();
        for i in 0..slf.user_value.len() {
            if slf.user_valid[i] {
                recipients.push(slf.user_value[i].clone());
            }
        }
        let mut index = 0;
        for i in 0..slf.archive_share_requested.len() {
            if slf.archive_share_requested[i] {
                index = i;
            }
        }
        let path = format!(
            "{}/{}",
            slf.archive_dir.as_ref().unwrap(),
            slf.archive_list[index]
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
            slf.shares.push(Share {
                days_since_ce: days,
                user_id: user_name,
            });
        }
        SENDING_SHARE.store(true, SeqCst);
    }
    Command::perform(compute_share(), Message::CompleteDoShare)
}

pub fn submit_button_pressed(slf: &mut EncloneVisual) -> Command<Message> {
    slf.modified = true;
    let mut group_spec = true;
    let mut group_ids = Vec::<usize>::new();
    let s = slf.input_value.split(',').collect::<Vec<&str>>();
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
    if slf.compute_state != WaitingForRequest {
        Command::none()
    } else {
        if group_spec {
            slf.translated_input_value = slf.input_value.clone();
            let mut reply_text;
            let new = slf.translated_input_current();
            let args = new.split(' ').collect::<Vec<&str>>();
            if slf.h.input1_history.is_empty() && slf.h.input2_history.is_empty() {
                reply_text = "Group identifier can only be supplied if another \
                    command has already been run."
                    .to_string();
            } else if group_ids[group_ids.len() - 1] > slf.current_tables.len() {
                reply_text = "Group identifier is too large.".to_string();
            } else {
                let mut group_pics = Vec::<String>::new();
                let mut last_widths = Vec::<u32>::new();
                for x in group_ids.iter() {
                    group_pics.push(slf.current_tables[*x - 1].clone());
                    last_widths.push(slf.last_widths_value[*x - 1]);
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
            args2.push(format!("G={}", slf.translated_input_value));
            slf.output_value = reply_text.to_string();
            let hi = slf.h.history_index;
            slf.h
                .input1_history
                .insert(hi as usize, slf.h.input1_hist_uniq.len() as u32);
            slf.h.input1_hist_uniq.push(slf.input1_value.clone());
            slf.h
                .input2_history
                .insert(hi as usize, slf.h.input2_hist_uniq.len() as u32);
            slf.h.input2_hist_uniq.push(slf.input2_value.clone());
            slf.h
                .narrative_history
                .insert(hi as usize, slf.h.narrative_hist_uniq.len() as u32);
            slf.h.narrative_hist_uniq.push(slf.narrative_value.clone());
            slf.h
                .translated_input_history
                .insert(hi as usize, slf.h.translated_input_hist_uniq.len() as u32);
            slf.h
                .translated_input_hist_uniq
                .push(args2.iter().format(" ").to_string());
            slf.h
                .svg_history
                .insert(hi as usize, slf.h.svg_history[(hi - 1) as usize]);
            slf.h
                .summary_history
                .insert(hi as usize, slf.h.summary_history[(hi - 1) as usize]);
            slf.h
                .displayed_tables_history
                .insert(hi as usize, slf.h.displayed_tables_hist_uniq.len() as u32);
            slf.h
                .displayed_tables_hist_uniq
                .push(reply_text.to_string());
            slf.h
                .table_comp_history
                .insert(hi as usize, slf.h.table_comp_history[(hi - 1) as usize]);
            slf.h
                .last_widths_history
                .insert(hi as usize, slf.h.last_widths_history[(hi - 1) as usize]);
            slf.h.is_blank.insert(hi as usize, slf.is_blank_current());
            slf.h
                .descrip_history
                .insert(hi as usize, slf.h.descrip_history[(hi - 1) as usize]);
            slf.h.history_index += 1;
            if !TEST_MODE.load(SeqCst) {
                Command::none()
            } else {
                Command::perform(noop0(), Message::Capture)
            }
        } else {
            slf.compute_state = Thinking;
            slf.start_command = Some(Instant::now());
            // The following sleep is needed for button text to consistenly update.
            thread::sleep(Duration::from_millis(20));
            if slf.input_value.starts_with('#') && slf.cookbook.contains_key(&slf.input_value) {
                slf.translated_input_value = slf.cookbook[&slf.input_value].clone();
            } else {
                slf.translated_input_value = slf.input_value.clone();
            }
            USER_REQUEST.lock().unwrap().clear();
            USER_REQUEST
                .lock()
                .unwrap()
                .push(slf.translated_input_value.clone());
            PROCESSING_REQUEST.store(true, SeqCst);
            Command::perform(compute(), Message::ComputationDone)
        }
    }
}
