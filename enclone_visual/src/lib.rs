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
pub mod gui;

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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is copied from the palaver crate 0.3.0.  We copied it because it had not been updated
// recently and was causing crate duplication.

// Count the number of threads of the current process. Uses
// [`/proc/self/stat`](http://man7.org/linux/man-pages/man5/proc.5.html):`num_threads` on Linux,
// [`task_threads`](http://web.mit.edu/darwin/src/modules/xnu/osfmk/man/task_threads.html)
// on macOS.

use std::convert::TryInto;

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
