// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Note: if enclone is run from the command line, and fails, it will still return exit status
// zero.  As far as we know, in all other cases where it is not run from the command line, it
// returns exit status zero.

use enclone_main::main_enclone::main_enclone;
use enclone_main::USING_PAGER;
use enclone_visual::enclone_client::enclone_client;
use enclone_visual::enclone_server::enclone_server;
use nix::sys::signal::{kill, SIGINT};
use nix::unistd::getppid;
use nix::unistd::Pid;
use pretty_trace::*;
use std::env;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t = Instant::now();
    PrettyTrace::new().on();
    let mut args: Vec<String> = env::args().collect();
    let mut no_kill = false;
    for i in 1..args.len() {
        if args[i] == "NO_KILL" {
            no_kill = true;
        }
    }

    // Client run of enclone.

    for i in 0..args.len() {
        if args[i] == "VIS" || args[i].starts_with("VIS=") {
            enclone_client(&t).await?;
        }
    }

    // Standard run of enclone.

    if args.len() < 2 || args[1] != "SERVER" {
        let res = main_enclone(&mut args).await;
        if res.is_err() {
            // TURNED OFF BECAUSE WE GOT EXIT STATUS ZERO SOMETIMES WHEN WE USED THROUGH COMMAND.
            //
            // If there was an error and we had used the pager, then std::process::exit(1) will
            // result in exit status 0 if enclone was invoked from a terminal window, and
            // probably not otherwise.  To get nonzero exit status, we instead kill the parent
            // process, which is less.  That's a little surprising, but that's how it works:
            // it is the parent that is less.
            //
            // The little big of sleep seems to prevent printing of an extra newline, but this
            // is flaky.
            //
            // The kill makes the screen flash.  This is pretty horrible.

            eprintln!("{}", res.unwrap_err());
            if !no_kill && USING_PAGER.load(SeqCst) && 0 == 1 {
                thread::sleep(Duration::from_millis(10));
                let ppid = getppid();
                kill(Pid::from_raw(i32::from(ppid)), SIGINT).unwrap();
                thread::sleep(Duration::from_millis(10));
            } else {
                std::process::exit(1);
            }
        }
        std::process::exit(0);
    }

    // Server run of enclone.

    enclone_server().await?;
    Ok(())
}
