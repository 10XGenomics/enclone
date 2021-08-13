// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::history::*;
use crate::*;
use chrono::prelude::*;
use std::time::Instant;

pub fn update_shares(slf: &mut gui_structures::EncloneVisual) {
    let t = Instant::now();

    // Import shares.

    GET_MY_SHARES.store(true, SeqCst);
    while GET_MY_SHARES.load(SeqCst) {
        thread::sleep(Duration::from_millis(10));
    }
    let k = RECEIVED_SHARES_CONTENT.lock().unwrap().len();
    let mut new_filenames = Vec::<String>::new();
    for i in 0..k {
        let bytes = &RECEIVED_SHARES_CONTENT.lock().unwrap()[i];
        let origin = RECEIVED_SHARES_MESSAGES.lock().unwrap()[i].clone();
        let mut evh = EncloneVisualHistory::restore_from_bytes(&bytes).unwrap();
        evh.origin = origin;
        let dir;
        if VISUAL_HISTORY_DIR.lock().unwrap().len() > 0 {
            dir = VISUAL_HISTORY_DIR.lock().unwrap()[0].clone();
        } else {
            dir = format!("{}/history", slf.visual);
        }
        let mut now = format!("{:?}", Local::now());
        now = now.replace("T", "___");
        now = now.before(".").to_string();
        let filename = format!("{}.{}", now, i + 1);
        let path = format!("{}/{}", dir, filename);
        let res = write_enclone_visual_history(&evh, &path);
        if res.is_err() {
            xprintln!(
                "Was Unable to write history to the file {}, \
                so Save on Exit failed.\n",
                path
            );
            std::process::exit(1);
        }
        new_filenames.push(filename);
    }
    prepend_to_vec(&mut slf.archive_list, &new_filenames);
    prepend_to_vec(&mut slf.restore_requested, &vec![false; k]);
    prepend_to_vec(&mut slf.delete_requested, &vec![false; k]);
    prepend_to_vec(&mut slf.deleted, &vec![false; k]);
    prepend_to_vec(&mut slf.expand_archive_entry, &vec![false; k]);
    prepend_to_vec(&mut slf.restore_msg, &vec![String::new(); k]);
    prepend_to_vec(&mut slf.archived_command_list, &vec![None; k]);
    prepend_to_vec(
        &mut slf.archive_name,
        &vec![iced::text_input::State::default(); k],
    );
    prepend_to_vec(&mut slf.archive_name_value, &vec![String::new(); k]);
    prepend_to_vec(&mut slf.archive_name_change_requested, &vec![false; k]);
    prepend_to_vec(&mut slf.archive_share_requested, &vec![false; k]);
    prepend_to_vec(&mut slf.archive_origin, &vec![String::new(); k]);

    // Delete remote shares.

    RELEASE_MY_SHARES.store(true, SeqCst);
    while RELEASE_MY_SHARES.load(SeqCst) {
        thread::sleep(Duration::from_millis(10));
    }

    // Sleep so that total time for the update_shares function is at least 0.4 seconds.  This 
    // keeps the "Receive shares" button red for at least that amount of time.

    const MIN_SLEEP: f64 = 0.4;
    let used = elapsed(&t);
    if used < MIN_SLEEP {
        let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
        thread::sleep(Duration::from_millis(ms));
    }
}
