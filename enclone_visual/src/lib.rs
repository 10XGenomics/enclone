// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod enclone_client;
pub mod enclone_server;

pub mod proto {
    tonic::include_proto!("enclone");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Global variables for client.

use lazy_static::lazy_static;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Mutex;

pub static REMOTE: AtomicBool = AtomicBool::new(false);
pub static USING_SETUP: AtomicBool = AtomicBool::new(false);
pub static CLEANED_UP: AtomicBool = AtomicBool::new(false);

pub static REMOTE_SERVER_ID: AtomicUsize = AtomicUsize::new(0);
pub static SETUP_PID: AtomicUsize = AtomicUsize::new(0);

pub static PROCESSING_REQUEST: AtomicBool = AtomicBool::new(false);

pub static DONE: AtomicBool = AtomicBool::new(false);

lazy_static! {
    pub static ref VERSION: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref USER_REQUEST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_TEXT: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_SVG: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref CONFIG_FILE: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}
