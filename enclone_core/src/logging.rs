// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use lazy_static::lazy_static;
use std::fs::OpenOptions;
use std::io::Write;
use std::sync::Mutex;

lazy_static! {
    pub static ref SERVER_LOGFILE: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

pub fn logme(s: &str) {
    if SERVER_LOGFILE.lock().unwrap().len() > 0 {
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&SERVER_LOGFILE.lock().unwrap()[0])
            .unwrap();
        writeln!(file, "{}", s).unwrap();
    }
}
