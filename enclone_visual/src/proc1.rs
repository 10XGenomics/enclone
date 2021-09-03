// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::history::*;
use crate::messages::*;
use crate::share::*;
use crate::*;
use chrono::prelude::*;
use enclone_core::combine_group_pics::*;
use flate2::read::GzDecoder;
use gui_structures::ComputeState::*;
use iced::{Color, Command};
use io_utils::*;
use itertools::Itertools;
use std::fs::File;
use std::io::Read;
use std::time::{Duration, Instant};
use vector_utils::*;

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

pub fn do_computation_done(slf: &mut EncloneVisual) -> Command<Message> {
    let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
    if reply_text.contains("enclone failed") {
        reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
    }
    if reply_text.len() == 0 {
        if slf.translated_input_value.contains(" NOPRINT") {
            reply_text = "You used the NOPRINT option, so there are no \
                clonotypes to see."
                .to_string();
        } else {
            reply_text = "There are no clonotypes.  Please have a look at the summary.".to_string();
        }
    }

    // Start storing values.

    let reply_table_comp = SERVER_REPLY_TABLE_COMP.lock().unwrap()[0].clone();
    slf.table_comp_value = reply_table_comp.clone();
    let hi = slf.h.history_index;
    let len = slf.h.table_comp_hist_uniq.len();
    if len > 0 && slf.h.table_comp_hist_uniq[len - 1] == reply_table_comp {
        slf.h
            .table_comp_history
            .insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.table_comp_history.insert(hi as usize, len as u32);
        slf.h.table_comp_hist_uniq.push(reply_table_comp.clone());
        if slf.table_comp_value.len() > 0 {
            let mut gunzipped = Vec::<u8>::new();
            let mut d = GzDecoder::new(&*reply_table_comp);
            d.read_to_end(&mut gunzipped).unwrap();
            slf.current_tables = serde_json::from_str(&strme(&gunzipped)).unwrap();
        } else {
            slf.current_tables.clear();
        }
    }

    // Get the summary, and stuff the dataset names and metrics into it.  The reason for this
    // grotesque operation was to avoid updating the history data structure.

    let reply_summary = SERVER_REPLY_SUMMARY_PLUS.lock().unwrap()[0].clone();

    // Keep going.

    reply_text += "\n \n \n"; // papering over truncation bug in display
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

    let len = slf.h.last_widths_hist_uniq.len();
    if len > 0 && slf.h.last_widths_hist_uniq[len - 1] == reply_last_widths {
        slf.h
            .last_widths_history
            .insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.last_widths_history.insert(hi as usize, len as u32);
        slf.h.last_widths_hist_uniq.push(reply_last_widths.clone());
    }
    let len = slf.h.svg_hist_uniq.len();
    if len > 0 && slf.h.svg_hist_uniq[len - 1] == reply_svg {
        slf.h.svg_history.insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.svg_history.insert(hi as usize, len as u32);
        slf.h.svg_hist_uniq.push(reply_svg.clone());
    }
    let len = slf.h.summary_hist_uniq.len();
    if len > 0 && slf.h.summary_hist_uniq[len - 1] == reply_summary {
        slf.h.summary_history.insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.summary_history.insert(hi as usize, len as u32);
        slf.h.summary_hist_uniq.push(reply_summary.clone());
    }
    slf.narrative_value.clear();
    let len = slf.h.narrative_hist_uniq.len();
    if len > 0 && slf.h.narrative_hist_uniq[len - 1] == slf.narrative_value {
        slf.h
            .narrative_history
            .insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.narrative_history.insert(hi as usize, len as u32);
        slf.h.narrative_hist_uniq.push(slf.narrative_value.clone());
    }
    let len = slf.h.displayed_tables_hist_uniq.len();
    if len > 0 && slf.h.displayed_tables_hist_uniq[len - 1] == reply_text {
        slf.h
            .displayed_tables_history
            .insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h
            .displayed_tables_history
            .insert(hi as usize, len as u32);
        slf.h.displayed_tables_hist_uniq.push(reply_text.clone());
    }
    let len = slf.h.input1_hist_uniq.len();
    if len > 0 && slf.h.input1_hist_uniq[len - 1] == slf.input1_value {
        slf.h.input1_history.insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.input1_history.insert(hi as usize, len as u32);
        slf.h.input1_hist_uniq.push(slf.input1_value.clone());
    }
    let len = slf.h.input2_hist_uniq.len();
    if len > 0 && slf.h.input2_hist_uniq[len - 1] == slf.input2_value {
        slf.h.input2_history.insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.input2_history.insert(hi as usize, len as u32);
        slf.h.input2_hist_uniq.push(slf.input2_value.clone());
    }
    let len = slf.h.descrip_hist_uniq.len();
    if len > 0 && slf.h.descrip_hist_uniq[len - 1] == slf.descrip_value {
        slf.h.descrip_history.insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h.descrip_history.insert(hi as usize, len as u32);
        slf.h.descrip_hist_uniq.push(slf.descrip_value.clone());
    }
    let len = slf.h.translated_input_hist_uniq.len();
    if len > 0 && slf.h.translated_input_hist_uniq[len - 1] == slf.translated_input_value {
        slf.h
            .translated_input_history
            .insert(hi as usize, (len - 1) as u32);
    } else {
        slf.h
            .translated_input_history
            .insert(hi as usize, len as u32);
        slf.h
            .translated_input_hist_uniq
            .push(slf.translated_input_value.clone());
    }
    slf.h.is_blank.insert(hi as usize, blank);
    slf.h.history_index += 1;
    slf.output_value = reply_text.to_string();
    slf.svg_value = reply_svg.to_string();
    slf.summary_value = reply_summary.to_string();
    slf.last_widths_value = reply_last_widths.clone();
    if slf.svg_value.len() > 0 {
        slf.post_svg(&reply_svg);
    }
    slf.compute_state = WaitingForRequest;
    xprintln!(
        "total time to run command = {:.1} seconds",
        elapsed(&slf.start_command.unwrap())
    );
    let maxrss_slf;
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        assert_eq!(retval, 0);
        maxrss_slf = rusage.ru_maxrss;
    }
    let peak_mem_mb = maxrss_slf as f64 / ((1024 * 1024) as f64);
    xprintln!(
        "all time peak mem of this process is {:.1} MB\n",
        peak_mem_mb
    );
    if VERBOSE.load(SeqCst) {
        let mb = (1024 * 1024) as f64;
        let mut total_svg = 0;
        for x in slf.h.svg_hist_uniq.iter() {
            total_svg += x.len();
        }
        xprintln!("stored svgs = {:.1} MB", total_svg as f64 / mb);
        let mut total_tables = 0;
        for x in slf.h.table_comp_hist_uniq.iter() {
            total_tables += x.len();
        }
        xprintln!("stored tables = {:.1} MB", total_tables as f64 / mb);
        xprintln!("");
    }

    if META_TESTING.load(SeqCst) {
        Command::perform(noop0(), Message::Meta)
    } else if !TEST_MODE.load(SeqCst) {
        Command::none()
    } else {
        slf.sanity_check();
        assert!(slf.h.save_restore_works());
        test_evh_read_write(&slf.h, "/tmp/evh_test");
        std::fs::remove_file("/tmp/evh_test").unwrap();
        Command::perform(noop0(), Message::Capture)
    }
}

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
    slf.input_value = slf.input1_value.clone();
    if slf.input1_value.len() > 0 && slf.input2_value.len() > 0 {
        slf.input_value += " ";
    }
    slf.input_value += &mut slf.input2_value.clone();
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
