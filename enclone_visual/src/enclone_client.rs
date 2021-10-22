// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

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
//
// Alternatively, you can define an environment variable ENCLONE_CONFIG=filename or
// ENCLONE_CONFIG=filehost:filename, in which case that file is fetched, and lines are found
// that look like
// vis.x.variable=value
// and in such cases, we use variable=value is the source of variable definitions.
//
// enclone VIS -- run the serve locally
//
// The special argument SERVER_DEBUG causes the server to print debuggin information.  However
// you will only see this if you run the server locally using enclone VIS.

use crate::client_requests::*;
use crate::launch_gui;
use crate::proto::analyzer_client::AnalyzerClient;
use crate::update_restart::*;
use crate::*;
use enclone_core::parse_bsv;
use enclone_core::prepare_for_apocalypse::*;
use enclone_core::REMOTE_HOST;
use enclone_version::*;
use io_utils::*;
use itertools::Itertools;
use libc::atexit;
use nix::sys::signal::{kill, SIGINT as SIGINT_nix};
use nix::unistd::Pid;
use perf_stats::peak_mem_usage_gb;
use std::env;
use std::io::Read;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn enclone_client(t: &Instant) -> Result<(), Box<dyn std::error::Error>> {
    //
    // Fail if not running on a Mac.

    #[cfg(not(any(target_os = "macos", target_os = "ios")))]
    {
        xprintln!(
            "\nenclone visual only runs on a Mac at present.  Please let us know if you\n\
            are interested in running it under Linux or Windows.\n"
        );
        if 0 == 0 {
            std::process::exit(1);
        }
    }

    // Set up to catch CTRL-C events.  Parse arguments.

    let _ = install_signal_handler();
    let args: Vec<String> = env::args().collect();
    let mut verbose = false;
    let mut monitor_threads = false;
    let mut config_name = String::new();
    let mut fixed_port = None;
    for i in 1..args.len() {
        let arg = &args[i];

        // General options.

        if arg.starts_with("VIS=") {
            config_name = arg.after("VIS=").to_string();
        } else if arg == "VIS" {
        } else if arg == "VERBOSE" {
            verbose = true;

        // Special testing options.
        } else if arg == "MONITOR_THREADS" {
            monitor_threads = true;
        } else if arg.starts_with("PORT=") {
            fixed_port = Some(arg.after("PORT=").parse::<u16>().unwrap());
            assert!(fixed_port.unwrap() >= 1024);
        } else if arg == "TEST" {
            TEST_MODE.store(true, SeqCst);
        } else if arg.starts_with("VISUAL_DIR=") {
            let dir = arg.after("VISUAL_DIR=").to_string();
            VISUAL_DIR.lock().unwrap().push(dir);
        } else if arg == "PLAYBACK" {
            PLAYBACK.store(true, SeqCst);
        } else if arg == "FAIL_ON_ERROR" {
            FAIL_ON_ERROR.store(true, SeqCst);
        } else if arg.starts_with("SERVER_LOGFILE=") {
            // don't use tilde because we don't expand it
            let server_logfile = arg.after("SERVER_LOGFILE=").to_string();
            SERVER_LOGFILE.lock().unwrap().push(server_logfile);
        } else if arg.starts_with("META=") {
            META_TESTING.store(true, SeqCst);
            META.store(arg.after("META=").force_usize() - 1, SeqCst);
        } else if arg.starts_with("EXEC=") {
            EXEC.lock().unwrap().push(arg.after("EXEC=").to_string());
        } else {
            xprintln!(
                "\nCurrently the only allowed arguments are VIS, VIS=x where x is a\n\
                configuration name and VERBOSE, and some special testing options.\n"
            );
            std::process::exit(1);
        }
    }
    for (key, _value) in env::vars() {
        if key == "ENCLONE_VERBOSE" {
            verbose = true;
        }
    }
    if verbose {
        VERBOSE.store(true, SeqCst);
    }

    // Set enclone visual version.

    let version = "0.00000000000001";
    VERSION.lock().unwrap().push(version.to_string());

    // Monitor threads.

    if monitor_threads {
        tokio::spawn(async move {
            loop {
                thread::sleep(Duration::from_millis(5000));
                xprintln!("[using {} threads", thread_count());
            }
        });
    }

    // Announce.

    let version_float = format!("1e-{}", -version.force_f64().log10());
    if !verbose {
        xprintln!(
            "\nHi! You are using enclone visual {} = {}.\n\n\
            If you get an error \
            message or a window does not pop up, and it's not clear what to do,\nplease \
            rerun the command with the added argument VERBOSE, and then ask for help.",
            version,
            version_float,
        );
    }

    // Get configuration from environment variable if defined.

    let mut config = HashMap::<String, String>::new();
    let mut configuration = None;
    let mut found = false;
    if config_name.len() > 0 {
        let env_var = format!("ENCLONE_VIS_{}", config_name);
        for (key, value) in env::vars() {
            if key == env_var {
                configuration = Some(value.clone());
            }
        }
        if configuration.is_some() {
            let configuration = configuration.unwrap();
            found = true;
            if verbose {
                xprintln!("\nusing configuration\n▓{}▓", configuration);
            }
            let x = parse_bsv(&configuration);
            for arg in x.iter() {
                if !arg.contains("=") {
                    xprintln!(
                        "\nYour configuration has an argument {} that does not contain =.\n",
                        arg
                    );
                    std::process::exit(1);
                }
                config.insert(arg.before("=").to_string(), arg.after("=").to_string());
            }
        }
    }

    // Get config file if defined.  The config file is specified by an environment variable
    // of the form ENCLONE_CONFIG=filename or ENCLONE_CONFIG=filehost:filename.
    //
    // Note that even if ENCLONE_VIS is specified, we still look for ENCLONE_CONFIG.

    let mut filehost = String::new();
    let mut filehost_used = false;
    let mut auto_update = false;
    let mut internal = false;
    let mut config_file_contents = String::new();
    let mut ssh_cat_time = 0.0;
    for (key, value) in env::vars() {
        if key == "ENCLONE_INTERNAL" {
            internal = true;
        }

        if key == "ENCLONE_CONFIG" {
            CONFIG_FILE.lock().unwrap().push(value.to_string());
            let hf = value.to_string();
            let mut filename = hf.clone();
            if hf.contains(":") {
                filehost = hf.before(":").to_string();
                filename = hf.after(":").to_string();
            }
            if filehost.len() > 0 {
                REMOTE_HOST.lock().unwrap().push(filehost.clone());
            }

            // If filename exists on this server, use it, and ignore filehost.

            if path_exists(&filename) {
                config_file_contents = std::fs::read_to_string(&filename).unwrap();

            // Otherwise, fetch the file from the host.
            } else if filehost.len() > 0 && config_name.len() > 0 {
                filehost_used = true;
                let t = Instant::now();
                let o = Command::new("ssh")
                    .arg(&filehost)
                    .arg("cat")
                    .arg(&filename)
                    .output()
                    .expect("failed to execute ssh cat");
                ssh_cat_time = elapsed(&t);
                xprintln!("\nssh cat to {} took {:.1} seconds", filehost, ssh_cat_time);
                if o.status.code() != Some(0) {
                    let m = String::from_utf8(o.stderr).unwrap();
                    xprintln!("\ntest ssh failed with error message =\n{}", m);
                    xprintln!(
                        "Attempt to ssh to {} as specified by environment variable \
                        ENCLONE_CONFIG failed.",
                        filehost,
                    );
                    xprintln!("Here are two possible explanations:");
                    xprintln!("1. The host is wrong.");
                    xprintln!("2. You are not connected to the internet.");
                    xprintln!(
                        "3. You first need to do something to enable crossing \
                              a firewall.  If so, ask a colleague.\n"
                    );
                    std::process::exit(1);
                }
                config_file_contents = strme(&o.stdout).to_string();
            }
        }
    }

    // Get configuration.

    for line in config_file_contents.lines() {
        if line == "visual_auto_update=true" {
            auto_update = true;
        }
        if line == "internal=true" {
            internal = true;
        }
        if line.starts_with("cookbooks=") {
            let cb = line.after("cookbooks=").split(',').collect::<Vec<&str>>();
            for i in 0..cb.len() {
                COOKBOOK_DIRS.lock().unwrap().push(cb[i].to_string());
            }
        }
        if config_name.len() > 0 {
            let prefix = format!("vis.{}.", config_name);
            if line.starts_with(&prefix) {
                let def = line.after(&prefix);
                if def.contains("=") {
                    config.insert(def.before("=").to_string(), def.after("=").to_string());
                    found = true;
                }
            }
        }
    }
    if config_name.len() > 0 && !found {
        xprintln!(
            "\nYou specified the configuration name {}, but the content of that configuration \
               was not found.\n",
            config_name,
        );
        if CONFIG_FILE.lock().unwrap().len() > 0 {
            xprintln!(
                "The value of ENCLONE_CONFIG is {}.\n",
                CONFIG_FILE.lock().unwrap()[0]
            );
        } else {
            xprintln!("The environment variable ENCLONE_CONFIG is not defined.\n");
        }
        xprintln!("Here is what's in your configuration file:\n");
        for line in config_file_contents.lines() {
            xprintln!("{}", line);
        }
        xprintln!("");
        std::process::exit(1);
    }
    if internal {
        INTERNAL.store(true, SeqCst);
    }

    // Set up proper tracebacks.

    let mut bug_reports = "enclone@10xgenomics.com".to_string();
    for i in 1..args.len() {
        if args[i] == "BUG_REPORTS" {
            bug_reports = "".to_string();
        } else if args[i].starts_with("BUG_REPORTS=") {
            bug_reports = args[i].after("BUG_REPORTS=").to_string();
        }
    }
    for (key, value) in env::vars() {
        if key == "ENCLONE_BUG_REPORTS" {
            bug_reports = value.to_string();
        }
    }
    prepare_for_apocalypse(&args, internal, &bug_reports);
    BUG_REPORTS.lock().unwrap().push(bug_reports);

    // Save remote share.

    if config.contains_key("REMOTE_SHARE") {
        REMOTE_SHARE
            .lock()
            .unwrap()
            .push(config["REMOTE_SHARE"].clone());
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
            xprintln!(
                "\nTo use a remote host, please specify all of REMOTE_HOST, \
                REMOTE_IP, and REMOTE_BIN.\n"
            );
            xprintln!("Here is what is specified:");
            for (key, value) in config.iter() {
                xprintln!("{}={}", key, value);
            }
            std::process::exit(1);
        } else {
            REMOTE.store(true, SeqCst);
        }
    }

    // If server is remote, see if we can ssh to it.  Otherwise, busted.

    if remote {
        let host = config["REMOTE_HOST"].clone();
        if !filehost_used || filehost != host {
            let t = Instant::now();
            let o = Command::new("ssh")
                .arg(&host)
                .arg("-n")
                .arg("echo")
                .output()
                .expect("failed to execute initial ssh");
            xprintln!(
                "\ninitial test ssh to {} took {:.1} seconds",
                host,
                elapsed(&t)
            );
            if o.status.code() != Some(0) {
                let m = String::from_utf8(o.stderr).unwrap();
                xprintln!("\ntest ssh to {} failed with error message =\n{}", host, m);
                xprintln!(
                    "Attempt to ssh to {} as specified by REMOTE_HOST failed.",
                    host
                );
                xprintln!("Here are two possible explanations:");
                xprintln!("1. You have the wrong REMOTE_HOST.");
                xprintln!("2. You first need to do something to enable crossing a firewall.");
                xprintln!("   If so, ask one of your colleagues how to do this.\n");
                std::process::exit(1);
            }
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
        let mut port: u16 = (nanos % 65536) as u16;
        if fixed_port.is_some() {
            port = fixed_port.unwrap();
        }
        if port < 1024 {
            continue;
        }

        // Attempt to fork the server.

        let server_process;
        let mut local_host = "127.0.0.1".to_string();
        xprintln!("\ntrying random port {}", port);
        if remote {
            let host = config["REMOTE_HOST"].clone();
            HOST.lock().unwrap().push(host.clone());
            let ip = &config["REMOTE_IP"];
            let bin = &config["REMOTE_BIN"];
            if verbose {
                xprintln!(
                    "\nstarting remote server using\nssh {} {}/enclone {}:{} SERVER",
                    host,
                    bin,
                    ip,
                    port
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
            let t = Instant::now();
            server_process = Command::new("enclone")
                .arg("SERVER")
                .arg(&format!("{}:{}", ip, port))
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
            xprintln!("used {:.1} seconds launching local server", elapsed(&t));
        }
        if verbose {
            xprintln!("server forked");
        }
        if !server_process.is_ok() {
            xprintln!(
                "\nfailed to launch server, err =\n{}.\n",
                server_process.as_ref().unwrap_err()
            );
            std::process::exit(1);
        }
        let mut server_process = server_process.unwrap();
        let server_process_id = server_process.id();
        SERVER_PROCESS_PID.store(server_process_id as usize, SeqCst);

        // Wait until server has printed something.

        let mut buffer = [0; 50];
        let server_stdout = server_process.stdout.as_mut().unwrap();
        let tread = Instant::now();
        if verbose {
            xprintln!("waiting for server response");
        }
        pub static READ_DONE: AtomicBool = AtomicBool::new(false);
        thread::spawn(move || {
            // At one time we used a sleep time of 5 seconds, but some people have slower
            // connections and this resulted in failures.  The heuristic in the next
            // line is the solution, but perhaps there is a better way to handle this.
            let sleep_time = 5.0_f64.max(2.0 * ssh_cat_time);
            thread::sleep(Duration::from_millis((sleep_time * 1000.0).round() as u64));
            if !READ_DONE.load(SeqCst) {
                xprintln!("darn, we seem to have hit a bad port, so restarting");
                restart_enclone();
            }
        });
        server_stdout.read(&mut buffer).unwrap();
        READ_DONE.store(true, SeqCst);
        if verbose {
            xprintln!(
                "time spent waiting to read bytes from server = {:.1} seconds",
                elapsed(&tread)
            );
        }

        // Look at stderr.  We read exactly 200 bytes.  By design, this is enough to know that the
        // server succeeded and enough to contain the information that the server is passing to
        // the client.

        let mut ebuffer = [0; 200];
        let server_stderr = server_process.stderr.as_mut().unwrap();
        let tread = Instant::now();
        server_stderr.read_exact(&mut ebuffer).unwrap();
        if verbose {
            xprintln!(
                "used {:.1} seconds reading from server stderr",
                elapsed(&tread)
            );
        }
        let emsg = strme(&ebuffer);
        if emsg.len() > 0 {
            if emsg.contains("already in use") {
                xprintln!("oops, that port is in use, trying a different one");
                continue;
            }
            if verbose {
                xprintln!("\nserver says this:\n{}", emsg);
            }
            if emsg.contains("tput: No value") {
                xprintln!(
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
                xprintln!("\nUnable to determine remote process id.\n");
                xprintln!("message = {}", emsg);
                std::process::exit(1);
            }
            let remote_version;
            if emsg.contains("enclone version = ")
                && emsg.after("enclone version = ").contains("\n")
            {
                remote_version = emsg.between("enclone version = ", "\n").to_string();
            } else {
                xprint!("\nUnable to determine remote enclone version.\n");
                std::process::exit(1);
            }
            let remote_version_string;
            if emsg.contains("version string = ") && emsg.after("version string = ").contains("\n")
            {
                remote_version_string = emsg.between("version string = ", "\n").to_string();
            } else {
                xprint!("\nUnable to determine remote version string.\n");
                std::process::exit(1);
            }
            let local_version = env!("CARGO_PKG_VERSION");
            if local_version != remote_version {
                xprintln!("\nremote enclone version = {}", remote_version);
                xprintln!("local enclone version = {}", local_version);
                xprintln!("\nYour enclone version is not up to date.");
                if auto_update {
                    update_enclone();
                    xprintln!("Done, restarting!\n");
                    restart_enclone();
                } else {
                    xprintln!(
                        "Please update, following \
                        the instructions at bit.ly/enclone, then restart.  Thank you!\n"
                    );
                    std::process::exit(1);
                }
            }

            // Check for identity of local and remote enclone versions.  This is complicated
            // because in general, the version_string() function does not return the current
            // value.  So we only test for version identity if we're in the directory where
            // enclone as compiled.  And that same condition is partially tested on the remote.

            let current_dir = std::env::current_dir()?;
            let current_dir = current_dir.display();
            let current_executable = std::env::current_exe()?;
            let current_executable = current_executable.display();
            println!("current dir = {}", current_dir);
            println!("current executable = {}", current_executable);
            if format!("{}", current_executable) == format!("{}/target/debug/enclone", current_dir)
            {
                let local_version_string = current_version_string();
                let local = format!("{} : {}", local_version, local_version_string);
                let remote = format!("{} : {}", remote_version, remote_version_string);
                let mut xlocal = local.clone();
                xlocal = xlocal.replace(": macos :", "");
                xlocal = xlocal.replace(": linux :", "");
                let mut xremote = remote.clone();
                xremote = xremote.replace(": macos :", "");
                xremote = xremote.replace(": linux :", "");
                if xlocal != xremote && verbose {
                    xprintln!("\n🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴");
                    xprintln!("----------------------------------------------------------------------------------");
                    xprintln!("😱 WARNING: INCOMPATIBLE SERVER/CLIENT VERSIONS DETECTED!!! 😱");
                    xprintln!("local version = {}", local);
                    xprintln!("remote version = {}", remote);
                    xprintln!("THIS CAN CAUSE VERY BAD AND INSCRUTABLE THINGS TO HAPPEN!");
                    xprintln!("😱 PROCEED AT RISK.................😱");
                    xprintln!("----------------------------------------------------------------------------------");
                    xprintln!("🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴🔴");
                } else if verbose {
                    xprintln!("local version = {}", local);
                    xprintln!("remote version = {}", remote);
                }
            }
        }

        // Form local URL.

        let url = format!("http://{}:{}", local_host, port);

        // Fork remote setup command if needed.

        if config.contains_key("REMOTE_SETUP") {
            let tremote = Instant::now();
            let mut setup = config["REMOTE_SETUP"].clone();
            if setup.starts_with("\"") && setup.ends_with("\"") {
                setup = setup.after("\"").rev_before("\"").to_string();
            }
            setup = setup.replace("$port", &format!("{}", port));
            let argsp = setup.split(' ').collect::<Vec<&str>>();
            let args = argsp[1..].to_vec();
            if verbose {
                xprintln!("\nrunning setup command = {}", argsp.iter().format(" "));
            }
            let setup_process = Command::new(argsp[0])
                .args(args)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
            if !setup_process.is_ok() {
                xprintln!(
                    "\nfailed to launch setup, err =\n{}.\n",
                    setup_process.as_ref().unwrap_err()
                );
                kill(Pid::from_raw(server_process_id as i32), SIGINT_nix).unwrap();
                cleanup();
                std::process::exit(1);
            }
            let setup_process = setup_process.unwrap();
            let setup_process_id = setup_process.id();
            SETUP_PID.store(setup_process_id as usize, SeqCst);
            USING_SETUP.store(true, SeqCst);
            if verbose {
                xprintln!("used {:.1} seconds connecting to remote", elapsed(&tremote));
            }
        }

        // Connect to client.

        if verbose {
            xprintln!("connecting to {}", url);
        }
        let tconnect = Instant::now();
        const MAX_CONNECT_MS: u64 = 10_000;
        const CONNECT_WAIT_MS: u64 = 250;
        use tonic::transport::Channel;
        let mut client: Result<AnalyzerClient<Channel>, tonic::transport::Error>;
        let mut wait_time = 0;
        loop {
            thread::sleep(Duration::from_millis(CONNECT_WAIT_MS));
            wait_time += CONNECT_WAIT_MS;
            client = AnalyzerClient::connect(url.clone()).await;
            if client.is_ok() {
                break;
            }
            if wait_time >= MAX_CONNECT_MS {
                xprintln!("\nconnection failed with error\n{:?}\n", client);
                if !verbose {
                    xprintln!("Please retry after adding the VERBOSE argument to your command.\n");
                }
                xprintln!("Please report this problem.  It is possible that the maximum");
                xprintln!("connection time used by enclone visual needs to be increased.\n");
                cleanup();
                std::process::exit(1);
            }
        }
        if verbose {
            xprintln!("used {:.1} seconds connecting", elapsed(&tconnect));
        }
        xprintln!("connected");
        xprintln!(
            "time since startup = {:.1} seconds, peak mem = {:.1} GB\n",
            elapsed(&t),
            peak_mem_usage_gb()
        );
        let mut client = client.unwrap();

        // Process commands via the server in the background.

        tokio::spawn(async move {
            process_requests(&mut client, &mut server_process, verbose).await;
        });

        // Launch GUI.

        launch_gui().await?;
        cleanup();
        return Ok(());
    }
}
