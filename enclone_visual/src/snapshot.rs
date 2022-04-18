// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.

use crate::capture_as_file;
use crate::copy_image_to_clipboard::*;
use crate::get_window_id;
use perf_stats::*;
use std::env;
use std::fs::{remove_file, File};
use std::io::Read;
use std::thread;
use std::time::{Duration, Instant};

// Copy window image to clipboard.  If the environment variable ENCLONE_VIS_SNAPSHOT is defined,
// also save to that file.

pub fn snapshot(start: &Option<Instant>) {
    let mut filename = "/tmp/enclone_visual_snapshot.png".to_string();
    let mut snapshot = false;
    for (key, value) in env::vars() {
        if key == "ENCLONE_VIS_SNAPSHOT" {
            snapshot = true;
            filename = value.to_string();
        }
    }
    capture_as_file(&filename, get_window_id());
    let mut bytes = Vec::<u8>::new();
    {
        let mut f = File::open(&filename).unwrap();
        f.read_to_end(&mut bytes).unwrap();
    }
    if !snapshot {
        remove_file(&filename).unwrap();
    }
    copy_png_bytes_to_clipboard(&bytes);
    const MIN_SLEEP: f64 = 0.4;
    let used = elapsed(&start.unwrap());
    if used < MIN_SLEEP {
        let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
        thread::sleep(Duration::from_millis(ms));
    }
}
