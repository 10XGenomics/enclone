// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// TO DO
//
// 1.  Make sure server dies.
//
// 2.  Need a regression test.
//
// 3.  Trim features and duplicated crates; reduce binary size.
//
// 4.  Vertical placement of legend in PLOT_BY_ISOTYPE is not great.
//
// 5.  Handle the case where tsh hasn't been started
//     server says this: david.jaffe@tsh-jump.txgmesh.net: Permission denied (publickey).
//     kex_exchange_identification: Connection closed by remote host
//
// 6.  What do we do about the truncation of printing to fifty clonotypes?
//
// 7.  tooltip
//
// 8.  the wraparound problem
//
// 9.  Make sure that client and server are the same version.
//
// 10. Handle the case where button is pushed twice, etc.
//
// 11. Have text-only mode for testing and development.
//
// 12. Need auto-update for binary (at least for Mac).
//
// 13. Clonotype printouts need to show not just the pictures, but also the group info
//     and any additional stuff like trees.
//
// 14. Improve data shown for tooltip.
//
// 15. Refactor to keep tokio etc. out of cellranger.
//
// 16. Make code work for 0, 2 or more SVG files, see the code svgs[0] which is broken in general.
//
// WAITING ON ICED
//
// 1.  Can't cut and paste text from the GUI window, except for the text input box.
//     Looks like this is https://github.com/hecrj/iced/issues/36.
//
// 2.  Pretty, not plain.
//     Enabling e.g. multicolor text is on the iced roadmap.
//
// 3.  Text in SVG objects does not show up.
//     Known regression = https://github.com/hecrj/iced/issues/870.
//
// 4.  Place the scrollbar on the left side of the scrollable window.
//     Asked on zulip chat if this is possible.
//
// 5.  Make the TextInput gracefully handle long inputs.
//     This is https://github.com/hecrj/iced/issues/320.  Fix targetted for 1.0.
//
// 6.  Flaky behavior scrolling large texts.
//     Issue filed: https://github.com/hecrj/iced/issues/897.
//
// NICE TO HAVE
//
// 1.  Add font size slider for clonotype tables.
//
// 2.  Can carriage return be used instead of pushing a button?
//
// 3.  Bold some text on the help page.
//     Probably not possible at this time, could ask on zulip chat.
//
// 4.  Use native ssh calls rather than ssh commands, as this might be faster and more robust.
//
// 5.  Enable plotting after initially starting GUI.
//
// 6.  Put plots in a scrollable window.
//
// 7.  Add scale slider for plots.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is a GUI client, work in progress.
//
// This starts enclone in SERVER mode, which can be either local or remote.
//
// For now accepts a single argument, which is VIS=x where x is a "configuration name".  It then
// looks for an environment variable ENCLONE_VIS_x, and parses that into blank-separated
// arguments, which may be:
//
// argument                  interpretation
//
// REMOTE_HOST=...           name of remote host for password-free ssh
// REMOTE_IP=...             IP number of remote host
// REMOTE_SETUP=...          command to be forked to use port through firewall, may include $port
// REMOVE_BIN=...            directory on remote host containing the enclone executable

use crate::proto::{analyzer_client::AnalyzerClient, ClonotypeRequest, EncloneRequest};
use crate::*;
use enclone_core::parse_bsv;
use gui::launch_gui;
use itertools::Itertools;
use libc::atexit;
use nix::sys::signal::{kill, SIGINT as SIGINT_nix};
use nix::unistd::Pid;
use perf_stats::*;
use std::collections::HashMap;
use std::env;
use std::io::Read;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn enclone_client(t: &Instant) -> Result<(), Box<dyn std::error::Error>> {
    //
    // Set up to catch CTRL-C events.  Parse arguments.

    let _ = install_signal_handler();
    let args: Vec<String> = env::args().collect();
    let mut verbose = false;
    let mut monitor_threads = false;
    let mut config_name = String::new();
    for i in 1..args.len() {
        let arg = &args[i];
        if arg == "VERBOSE" {
            verbose = true;
        } else if arg == "MONITOR_THREADS" {
            monitor_threads = true;
        } else if arg.starts_with("VIS=") {
            config_name = arg.after("VIS=").to_string();
        } else if arg != "VIS" {
            eprintln!(
                "\nCurrently the only allowed arguments are VIS, VIS=x where x is a\n\
                configuration name, VERBOSE, and MONITOR_THREADS.\n"
            );
            std::process::exit(1);
        }
    }

    // Set enclone visual version.

    let version = "0.0000000000000000000000000000001";
    VERSION.lock().unwrap().push(version.to_string());

    // Monitor threads.

    if monitor_threads {
        tokio::spawn(async move {
            loop {
                thread::sleep(Duration::from_millis(5000));
                println!("[using {} threads", thread_count());
            }
        });
    }

    // Announce.

    let version_float = format!("1e-{}", -version.force_f64().log10());
    if !verbose {
        println!(
            "\nHi! You are using enclone visual {} = {}.\n\n\
            If you get an error \
            message or a window does not pop up, and it's not clear what to do,\nplease \
            rerun the command with the added argument VERBOSE, and then ask for help.",
            version, version_float,
        );
    }

    // Get config file name if defined.

    for (key, value) in env::vars() {
        if key == "ENCLONE_CONFIG" {
            CONFIG_FILE.lock().unwrap().push(value.to_string());
        }
    }

    // Get configuration.

    let mut configuration = None;
    if config_name.len() > 0 {
        let env_var = format!("ENCLONE_VIS_{}", config_name);
        for (key, value) in env::vars() {
            if key == env_var {
                configuration = Some(value.clone());
            }
        }
        if configuration.is_none() {
            eprintln!(
                "\nYou specified the configuration name {}, but the environment variable {} \
                is not defined.\n",
                config_name, env_var,
            );
            std::process::exit(1);
        } else if verbose {
            println!(
                "\nusing configuration\n▓{}▓",
                configuration.as_ref().unwrap()
            );
        }
    }
    let mut config = HashMap::<String, String>::new();
    if configuration.is_some() {
        let configuration = configuration.unwrap();
        let x = parse_bsv(&configuration);
        for arg in x.iter() {
            if !arg.contains("=") {
                eprintln!(
                    "\nYour configuration has an argument {} that does not contain =.\n",
                    arg
                );
                std::process::exit(1);
            }
            config.insert(arg.before("=").to_string(), arg.after("=").to_string());
        }
    }

    // Determine if the server is remote.

    let remote = config.contains_key("REMOTE_HOST")
        || config.contains_key("REMOTE_IP")
        || config.contains_key("REMOTE_BIN");
    if remote {
        if !config.contains_key("REMOTE_HOST")
            || !config.contains_key("REMOTE_IP")
            || !config.contains_key("REMOTE_BIN")
        {
            eprintln!(
                "\nTo use a remote host, please specify all of REMOTE_HOST, \
                REMOTE_IP, and REMOTE_BIN.\n"
            );
            eprintln!("Here is what is specified:");
            for (key, value) in config.iter() {
                eprintln!("{}={}", key, value);
            }
            std::process::exit(1);
        } else {
            REMOTE.store(true, SeqCst);
        }
    }

    // If server is remote, see if we can ssh to it.  Otherwise, busted.

    if remote {
        let t = Instant::now();
        let host = config["REMOTE_HOST"].clone();
        let o = Command::new("ssh")
            .arg(&host)
            .arg("-n")
            .arg("echo")
            .output()
            .expect("failed to execute initial ssh");
        println!("\ninitial test ssh took {:.1} seconds", elapsed(&t));
        if o.status.code() != Some(0) {
            let m = String::from_utf8(o.stderr).unwrap();
            println!("\ntest ssh failed with error message =\n{}", m);
            println!(
                "Attempt to ssh to {} as specified by REMOTE_HOST failed.",
                host
            );
            println!("Here are two possible explanations:");
            println!("1. You have the wrong REMOTE_HOST.");
            println!("2. You first need to do something to enable crossing a firewall.");
            println!("   If so, ask one of your colleagues how to do this.\n");
            std::process::exit(1);
        }
    }

    // Set exit handler to force cleanup at end of process.

    unsafe {
        atexit(exit_handler);
    }

    // Loop through random ports until we get one that works.

    loop {
        let start = SystemTime::now();
        let since_the_epoch = start
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards");
        let nanos = since_the_epoch.subsec_nanos() as u64;
        let port: u16 = (nanos % 65536) as u16;
        if port < 1024 {
            continue;
        }

        // Attempt to fork the server.

        let server_process;
        let mut local_host = "127.0.0.1".to_string();
        println!("\ntrying random port {}", port);
        if remote {
            let host = config["REMOTE_HOST"].clone();
            HOST.lock().unwrap().push(host.clone());
            let ip = &config["REMOTE_IP"];
            let bin = &config["REMOTE_BIN"];
            if verbose {
                println!(
                    "\nstarting remote server using\nssh {} {}/enclone {}:{} SERVER",
                    host, bin, ip, port
                );
            }
            server_process = Command::new("ssh")
                .arg("-n")
                .arg(&host)
                .arg(&format!("{}/enclone", bin))
                .arg("SERVER")
                .arg(&format!("{}:{}", ip, port))
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
            local_host = "localhost".to_string();
        } else {
            let ip = "127.0.0.1";
            server_process = Command::new("enclone")
                .arg("SERVER")
                .arg(&format!("{}:{}", ip, port))
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
        }
        if !server_process.is_ok() {
            eprintln!(
                "\nfailed to launch server, err =\n{}.\n",
                server_process.unwrap_err()
            );
            std::process::exit(1);
        }
        let mut server_process = server_process.unwrap();
        let server_process_id = server_process.id();

        // Wait until server has printed something.

        let mut buffer = [0; 50];
        let server_stdout = server_process.stdout.as_mut().unwrap();
        let tread = Instant::now();
        server_stdout.read(&mut buffer).unwrap();
        println!(
            "time spent waiting to read bytes from server = {:.1} seconds",
            elapsed(&tread)
        );
        // seems to not be needed
        // thread::sleep(Duration::from_millis(100));

        // Look at stderr.

        let mut ebuffer = [0; 200];
        let server_stderr = server_process.stderr.as_mut().unwrap();
        server_stderr.read(&mut ebuffer).unwrap();
        let emsg = strme(&ebuffer);
        if emsg.len() > 0 {
            if emsg.contains("already in use") {
                println!("oops, that port is in use, trying a different one");
                continue;
            }
            if verbose {
                println!("\nserver says this:\n{}", emsg);
            }
            if emsg.contains("tput: No value") {
                println!(
                    "\nThat's an odd message that we've observed once and don't \
                    understand.  When it happened,\nit was associated with setting the \
                    environment variable PS1 on the remote server.  If it happens\nto you, \
                    please contact us at enclone@10xgenomics.com.\n"
                );
            }
        }

        // Get server process id, possibly remote.

        let mut remote_id = None;
        if remote {
            if emsg.contains("I am process ") && emsg.after("I am process ").contains(".") {
                let id = emsg.between("I am process ", ".");
                if id.parse::<usize>().is_ok() {
                    remote_id = Some(id.force_usize());
                    REMOTE_SERVER_ID.store(remote_id.unwrap(), SeqCst);
                }
            }
            if remote_id.is_none() {
                eprintln!("\nUnable to determine remote process id.\n");
                std::process::exit(1);
            }
        }

        // Form local URL.

        let url = format!("http://{}:{}", local_host, port);

        // Fork remote setup command if needed.

        if config.contains_key("REMOTE_SETUP") {
            let mut setup = config["REMOTE_SETUP"].clone();
            if setup.starts_with("\"") && setup.ends_with("\"") {
                setup = setup.after("\"").rev_before("\"").to_string();
            }
            setup = setup.replace("$port", &format!("{}", port));
            let argsp = setup.split(' ').collect::<Vec<&str>>();
            let args = argsp[1..].to_vec();
            if verbose {
                println!("\nrunning setup command = {}", argsp.iter().format(" "));
            }
            let setup_process = Command::new(argsp[0])
                .args(args)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
            if !setup_process.is_ok() {
                eprintln!(
                    "\nfailed to launch setup, err =\n{}.\n",
                    setup_process.unwrap_err()
                );
                kill(Pid::from_raw(server_process_id as i32), SIGINT_nix).unwrap();
                cleanup();
                std::process::exit(1);
            }
            let setup_process = setup_process.unwrap();
            let setup_process_id = setup_process.id();
            SETUP_PID.store(setup_process_id as usize, SeqCst);
            USING_SETUP.store(true, SeqCst);
            // Reducing sleep time below to 500 ms causes frequent failures.
            thread::sleep(Duration::from_millis(1000));
        }

        // Connect to client.

        if verbose {
            println!("connecting to {}", url);
        }
        let mut client = AnalyzerClient::connect(url.clone()).await;
        if client.is_err() {
            // If connection failed, sleep and try again.  This happens maybe 10% of the time.

            println!("connection attempt failed, waiting one second and will try again");
            thread::sleep(Duration::from_millis(1000));
            client = AnalyzerClient::connect(url).await;

            // Test for second failure.

            if client.is_err() {
                eprintln!("\nconnection failed with error\n{:?}\n", client);
                cleanup();
                std::process::exit(1);
            } else {
                println!("excellent, it's OK now");
            }
        }
        println!("connected");
        println!("time since startup = {:.1} seconds\n", elapsed(&t));
        let mut client = client.unwrap();

        // Process commands via the server in the background.

        tokio::spawn(async move {
            loop {
                thread::sleep(Duration::from_millis(10));
                if DONE.load(SeqCst) {
                    break;
                }
                if PROCESSING_REQUEST.load(SeqCst) {
                    let input = USER_REQUEST.lock().unwrap()[0].clone();
                    let mut line = input.to_string();
                    let output;
                    let mut svg_output = String::new();
                    if line == "q" {
                        cleanup();
                        std::process::exit(0);
                    }
                    if line == "d" {
                        line = "BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui".to_string();
                    }
                    if line.parse::<usize>().is_ok() {
                        let n = line.force_usize();
                        if n == 0 {
                            output = "clonotype numbers start at 1\n".to_string();
                        } else {
                            let request = tonic::Request::new(ClonotypeRequest {
                                clonotype_number: (n - 1) as u32,
                            });
                            let response = client.get_clonotype(request).await;
                            if response.is_err() {
                                output = "clonotype number is too large\n".to_string();
                            } else {
                                let response = response.unwrap();
                                let r = response.into_inner();
                                output = r.table.clone();
                            }
                        }
                    } else {
                        if CONFIG_FILE.lock().unwrap().len() > 0 {
                            line += &mut format!(" CONFIG={}", CONFIG_FILE.lock().unwrap()[0]);
                        }
                        let request = tonic::Request::new(EncloneRequest { args: line });
                        let response = client.enclone(request).await;
                        if response.is_err() {
                            let left = r###"message: "\n"###;
                            let right = r###"\n""###;
                            let mut err = format!("{:?}", response);
                            if err.contains(&left) && err.after(&left).contains(&right) {
                                err = err.between(&left, &right).to_string();
                            }
                            err = err.replace("\\n", "\n");
                            output =
                                format!("\nThe enclone server is unhappy.  It says:\n\n{}", err);

                            let mut ebuffer = [0; 10000];
                            let server_stderr = server_process.stderr.as_mut().unwrap();
                            server_stderr.read(&mut ebuffer).unwrap();
                            let emsg = strme(&ebuffer);
                            println!("server error =\n{}\n", emsg);
                        } else {
                            let response = response.unwrap();
                            let r = response.into_inner();
                            svg_output = r.plot.clone();
                            output = format!("\n\n{}", r.table);
                        }
                        SERVER_REPLY_SVG.lock().unwrap().clear();
                        SERVER_REPLY_SVG.lock().unwrap().push(svg_output);
                    }
                    SERVER_REPLY_TEXT.lock().unwrap().clear();
                    SERVER_REPLY_TEXT.lock().unwrap().push(output.clone());
                    PROCESSING_REQUEST.store(false, SeqCst);
                }
            }
        });

        // Launch GUI.

        launch_gui();
        cleanup();
        return Ok(());
    }
}
