// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::copy_image_to_clipboard::*;
use crate::gui_structures::EncloneVisual;
use crate::testsuite::TESTS;
use convert_svg_to_png::*;
use failure::Error;
use iced::{Application, Font, Settings};
use lazy_static::lazy_static;
use libc::SIGINT;
use nix::sys::signal::{kill, Signal, SIGINT as SIGINT_nix};
use nix::sys::signal::{sigaction, SaFlags, SigAction, SigHandler, SigSet};
use nix::unistd::Pid;
use perf_stats::*;
use std::cmp::max;
use std::collections::HashMap;
use std::convert::TryInto;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Mutex;
use std::thread;
use std::time::{Duration, Instant};
use string_utils::*;
use svg_to_geometry::*;
use tables::*;

pub mod canvas_view;
pub mod compare_images;
pub mod convert_svg_to_png;
pub mod copy_image_to_clipboard;
pub mod enclone_client;
pub mod enclone_server;
pub mod geometry;
pub mod gui;
pub mod gui_structures;
pub mod messages;
pub mod process_messages;
pub mod style;
pub mod svg_to_geometry;
pub mod testsuite;
pub mod update_restart;

type MsgFn = fn(Result<(), String>) -> messages::Message;

async fn noop() -> Result<(), String> {
    // Increasing this time to 2000ms will prevent the screen from going dark on initialization
    // in test mode.
    thread::sleep(Duration::from_millis(100));
    Ok(())
}

async fn compute() -> Result<(), String> {
    let t = Instant::now();
    while PROCESSING_REQUEST.load(SeqCst) {
        thread::sleep(Duration::from_millis(10));
    }
    print!(
        "time used processing command = {:.1} seconds, ",
        elapsed(&t)
    );
    Ok(())
}

async fn flash_copy_image_button() -> Result<(), String> {
    Ok(())
}

const DEJAVU_BOLD: Font = Font::External {
    name: "DEJAVU_BOLD",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono-Bold.ttf"),
};

impl EncloneVisual {
    pub fn post_svg(&mut self, svg: &str) {
        self.png_value.clear();
        let geometry = svg_to_geometry(&svg, false);
        let mut using_geometry = false;
        if geometry.is_some() {
            let mut ok = true;
            for i in 0..geometry.as_ref().unwrap().len() {
                match &geometry.as_ref().unwrap()[i] {
                    crate::geometry::Geometry::Text(ttt) => {
                        if ttt.rotate != [0.0; 3] {
                            ok = false;
                        }
                    }
                    _ => {}
                }
            }
            if ok {
                using_geometry = true;
                self.canvas_view.state.geometry_value = geometry;
            }
        } else if VERBOSE.load(SeqCst) {
            println!("translation from svg to geometries failed");
        }
        if !using_geometry {
            self.canvas_view.state.geometry_value = None;
            self.png_value = convert_svg_to_png(&svg.as_bytes());
        }
    }
}

pub fn blank_svg() -> String {
    r###"<svg version="1.1" baseProfile="full" width="400" height="400"
xmlns="http://www.w3.org/2000/svg">
<rect x="0" y="0" width="400" height="400" style="fill:white" />
</svg>
"###
    .to_string()
}

pub mod proto {
    tonic::include_proto!("enclone");
}

pub async fn launch_gui() -> iced::Result {
    let mut settings = Settings::default();
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = (1100 as u32, 1060 as u32); // reasonable minimum size
    settings.window = window_settings;
    settings.exit_on_close_request = false;
    EncloneVisual::run(settings)
}

pub fn capture(count: usize, window_id: usize) {
    if TESTS[count - 1].2.len() > 0 {
        thread::sleep(Duration::from_millis(20));
        let o = std::process::Command::new("screencapture")
            .arg("-x")
            .arg(&format!("-l{}", window_id))
            .arg(&format!(
                "enclone_visual/outputs/{}.png",
                TESTS[count - 1].2
            ))
            .output()
            .expect("failed to execute screencapture");
        if o.status.code() != Some(0) {
            eprintln!("\nCall to screencapture failed.");
            eprintln!("stderr =\n{}\n", strme(&o.stderr));
            std::process::exit(1);
        }
    }
}

pub fn get_window_id() -> usize {
    let o = Command::new("GetWindowID")
        .arg("enclone")
        .arg("--list")
        .output()
        .expect("failed to execute GetWindowID");
    if o.status.code() != Some(0) {
        eprintln!("\nCall to GetWindowID failed.\n");
        std::process::exit(1);
    }
    let mut m = String::from_utf8(o.stdout).unwrap();
    m = m.replace("\n", "");
    if !m.contains("id=") || m.after("id=").parse::<usize>().is_err() {
        eprintln!("\nGetWindowId could not find id\n");
        std::process::exit(1);
    }
    m.after("id=").force_usize()
}

pub fn fold(line: &str, max_line: usize) -> Vec<String> {
    let mut pieces = Vec::<String>::new();
    let words = line.split(' ').collect::<Vec<&str>>();
    let mut current = String::new();
    let mut i = 0;
    while i < words.len() {
        if current.len() > 0 && current.len() + 1 + words[i].len() > max_line {
            pieces.push(current.clone());
            current.clear();
            i -= 1;
        } else if words[i].len() >= max_line {
            let mut w = words[i].as_bytes().to_vec();
            loop {
                let n = std::cmp::min(max_line, w.len());
                let sub = stringme(&w[0..n]);
                if n < w.len() {
                    pieces.push(sub);
                    w = w[n..w.len()].to_vec();
                } else {
                    current = stringme(&w);
                    break;
                }
            }
        } else if current.len() == 0 {
            current += &mut words[i].clone();
        } else {
            current += &mut format!(" {}", words[i]);
        }
        i += 1;
    }
    if current.len() > 0 {
        pieces.push(current);
    }
    pieces
}

pub fn format_cookbook() -> String {
    let c;
    if INTERNAL.load(SeqCst) {
        c = include_str!["cookbook.internal"];
    } else {
        c = include_str!["cookbook.external"];
    }
    let mut rows = Vec::<Vec<String>>::new();
    let row = vec![
        "tag".to_string(),
        "command".to_string(),
        "action".to_string(),
    ];
    rows.push(row);
    let mut row = Vec::<String>::new();
    for line in c.lines() {
        if line.len() > 0 {
            row.push(line.to_string());
            if row.len() == 3 {
                if rows.len() > 0 {
                    rows.push(vec!["\\hline".to_string(); 3]);
                }
                rows.push(row.clone());
                row.clear();
            }
        }
    }
    let mut rows2 = Vec::<Vec<String>>::new();
    for i in 0..rows.len() {
        let m1 = fold(&rows[i][1], 60);
        let m2 = fold(&rows[i][2], 60);
        if m1.len() == 1 && m2.len() == 1 {
            rows2.push(rows[i].clone());
        } else {
            for j in 0..max(m1.len(), m2.len()) {
                let mut row = Vec::<String>::new();
                if j == 0 {
                    row.push(rows[i][0].clone());
                } else {
                    row.push("".to_string());
                }
                if j < m1.len() {
                    row.push(m1[j].clone());
                } else {
                    row.push("".to_string());
                }
                if j < m2.len() {
                    row.push(m2[j].clone());
                } else {
                    row.push("".to_string());
                }
                rows2.push(row);
            }
        }
    }
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows2, 0, &b"l|l|l".to_vec(), false, true);
    log
}

pub fn parse_cookbook() -> HashMap<String, String> {
    let c;
    if INTERNAL.load(SeqCst) {
        c = include_str!["cookbook.internal"];
    } else {
        c = include_str!["cookbook.external"];
    }
    let mut lines = Vec::<String>::new();
    for line in c.lines() {
        if line.len() > 0 {
            lines.push(line.to_string());
        }
    }
    let mut h = HashMap::<String, String>::new();
    for i in (0..lines.len()).step_by(3) {
        h.insert(lines[i].clone(), lines[i + 1].clone());
    }
    h
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Global variables for client.

pub static REMOTE: AtomicBool = AtomicBool::new(false);
pub static USING_SETUP: AtomicBool = AtomicBool::new(false);
pub static CLEANED_UP: AtomicBool = AtomicBool::new(false);
pub static VERBOSE: AtomicBool = AtomicBool::new(false);
pub static COOKBOOK: AtomicBool = AtomicBool::new(false);
pub static SUMMARY: AtomicBool = AtomicBool::new(false);
pub static INTERNAL: AtomicBool = AtomicBool::new(false);
pub static TEST_MODE: AtomicBool = AtomicBool::new(false);
pub static PROCESSING_REQUEST: AtomicBool = AtomicBool::new(false);
pub static DONE: AtomicBool = AtomicBool::new(false);
pub static GROUP_ID_CLICKED_ON: AtomicBool = AtomicBool::new(false);

pub static REMOTE_SERVER_ID: AtomicUsize = AtomicUsize::new(0);
pub static SERVER_PROCESS_PID: AtomicUsize = AtomicUsize::new(0);
pub static SETUP_PID: AtomicUsize = AtomicUsize::new(0);
pub static COUNT: AtomicUsize = AtomicUsize::new(0);
pub static GROUP_ID: AtomicUsize = AtomicUsize::new(0);

lazy_static! {
    pub static ref VERSION: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref USER_REQUEST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_TEXT: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_SVG: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_SUMMARY: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SERVER_REPLY_TABLE_COMP: Mutex<Vec<Vec<u8>>> = Mutex::new(Vec::<Vec<u8>>::new());
    pub static ref SERVER_REPLY_LAST_WIDTHS: Mutex<Vec<Vec<usize>>> =
        Mutex::new(Vec::<Vec<usize>>::new());
    pub static ref CONFIG_FILE: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref COOKBOOK_CONTENTS: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref SUMMARY_CONTENTS: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
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
        } else {
            kill(
                Pid::from_raw(SERVER_PROCESS_PID.load(SeqCst) as i32),
                SIGINT_nix,
            )
            .unwrap();
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is copied from the palaver crate 0.3.0.  We copied it because it had not been updated
// recently and was causing crate duplication.

// Count the number of threads of the current process. Uses
// [`/proc/self/stat`](http://man7.org/linux/man-pages/man5/proc.5.html):`num_threads` on Linux,
// [`task_threads`](http://web.mit.edu/darwin/src/modules/xnu/osfmk/man/task_threads.html)
// on macOS.

pub fn thread_count() -> usize {
    #[cfg(any(target_os = "android", target_os = "linux"))]
    {
        // This pulls in the procfs crate.  Not really sure this is necessary.

        procfs::process::Process::myself()
            .unwrap()
            .stat
            .num_threads
            .try_into()
            .unwrap()
    }
    #[cfg(any(target_os = "macos", target_os = "ios"))]
    {
        use mach::{
            kern_return::{kern_return_t, KERN_SUCCESS},
            mach_types::thread_act_array_t,
            message::mach_msg_type_number_t,
            task::task_threads,
            traps::mach_task_self,
            vm_types::{vm_address_t, vm_map_t, vm_size_t},
        };
        use std::{mem, ptr};
        extern "C" {
            pub fn vm_deallocate(
                target_task: vm_map_t,
                address: vm_address_t,
                size: vm_size_t,
            ) -> kern_return_t;
        }

        let this_task = unsafe { mach_task_self() };

        let mut thread_list: thread_act_array_t = ptr::null_mut();
        let mut thread_count: mach_msg_type_number_t = 0;
        let kret = unsafe { task_threads(this_task, &mut thread_list, &mut thread_count) };
        assert_eq!(kret, KERN_SUCCESS);
        let thread_count: usize = thread_count.try_into().unwrap();

        for i in 0..thread_count {
            let kret = unsafe {
                mach::mach_port::mach_port_deallocate(
                    this_task,
                    *thread_list.offset(i.try_into().unwrap()),
                )
            };
            assert_eq!(kret, KERN_SUCCESS);
        }
        let kret = unsafe {
            vm_deallocate(
                this_task,
                thread_list as usize,
                mem::size_of_val(&*thread_list) * thread_count,
            )
        };
        assert_eq!(kret, KERN_SUCCESS);
        thread_count
    }
    #[cfg(not(any(
        target_os = "android",
        target_os = "linux",
        target_os = "macos",
        target_os = "ios"
    )))]
    unimplemented!()
}
