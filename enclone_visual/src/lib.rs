// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use failure::Error;
use lazy_static::lazy_static;
use libc::SIGINT;
use nix::sys::signal::{kill, Signal, SIGINT as SIGINT_nix};
use nix::sys::signal::{sigaction, SaFlags, SigAction, SigHandler, SigSet};
use nix::unistd::Pid;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Mutex;

pub mod enclone_client;
pub mod enclone_server;

pub mod proto {
    tonic::include_proto!("enclone");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Global variables for client.

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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn _truncate(s: &str) -> String {
    const MAX_LINES: usize = 10;
    let mut t = String::new();
    let mut extra = 0;
    for (i, line) in s.lines().enumerate() {
        if i < MAX_LINES {
            t += &mut format!("{}\n", line);
        } else {
            extra += 1;
        }
    }
    if extra > 0 {
        t += &mut format!("(+ {} more lines)", extra);
    }
    t
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Cleanup code to make sure processes are killed.

pub fn cleanup() {
    if !CLEANED_UP.load(SeqCst) {
        CLEANED_UP.store(true, SeqCst);
        if REMOTE.load(SeqCst) {
            if USING_SETUP.load(SeqCst) {
                kill(Pid::from_raw(SETUP_PID.load(SeqCst) as i32), SIGINT_nix).unwrap();
            }
            if HOST.lock().unwrap().len() > 0 {
                let host = &HOST.lock().unwrap()[0];
                let _ = Command::new("ssh")
                    .arg(&host)
                    .arg("kill")
                    .arg("-9")
                    .arg(&format!("{}", REMOTE_SERVER_ID.load(SeqCst)))
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn();
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Redirect SIGINT interrupts to the function "handler".  There may be issues with reliablity,
// since a CTRL-C could happen at any point, including in the memory manager.

pub fn install_signal_handler() -> Result<(), Error> {
    let handler = SigHandler::Handler(handler);
    let action = SigAction::new(handler, SaFlags::SA_RESTART, SigSet::empty());
    unsafe {
        sigaction(Signal::SIGINT, &action)?;
    }
    Ok(())
}

extern "C" fn handler(sig: i32) {
    if sig == SIGINT {
        cleanup();
        std::process::exit(0);
    }
}

pub extern "C" fn exit_handler() {
    cleanup();
}
