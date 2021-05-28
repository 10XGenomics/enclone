// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// TO DO
//
// 1.  Make sure server dies.
//
// 2.  Need a regression test.
//
// 3.  Doesn't properly handle connection refused.
//
// 4.  Vertical placement of legend in PLOT_BY_ISOTYPE is not great.
//
// 5.  Handle the case where tsh hasn't been started
//     server says this: david.jaffe@tsh-jump.txgmesh.net: Permission denied (publickey).
//     kex_exchange_identification: Connection closed by remote host
//
// 6.  Add local server capability.
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
// 13. Trim features and duplicated crates; reduce binary size.
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
// NICE TO HAVE
//
// 1.  Make font a little darker.
//
// 2.  Can carriage return be used instead of pushing a button?


// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is a GUI client, work in progress.
//
// This starts enclone in SERVER mode, which can be either local or remote.
//
// For now accepts a single argument, which is COM=x where x is a "configuration name".  It then
// looks for an environment variable ENCLONE_COM_x, and parses that into blank-separated
// arguments, which may be:
//
// argument                  interpretation
//
// REMOTE_HOST=...           name of remote host for password-free ssh
// REMOTE_IP=...             IP number of remote host
// REMOTE_SETUP=...          command to be forked to use port through firewall, may include $port
// REMOVE_BIN=...            directory on remote host containing the enclone executable

use crate::proto::{analyzer_client::AnalyzerClient, ClonotypeRequest, EncloneRequest};
use enclone_core::parse_bsv;
use failure::Error;
use iced::svg::Handle;
use iced::Length::Units;
use iced::{
    button, scrollable, text_input, Align, Button, Color, Column, /* Container, */ Element,
    Font, HorizontalAlignment, Length, Row, Rule, Sandbox, Scrollable, Settings, Svg, Text,
    TextInput, VerticalAlignment,
};
use iced_aw::{modal, Card, Modal};
use itertools::Itertools;
use lazy_static::lazy_static;
use libc::{atexit, SIGINT};
use nix::sys::signal::{kill, Signal, SIGINT as SIGINT_nix};
use nix::sys::signal::{sigaction, SaFlags, SigAction, SigHandler, SigSet};
use nix::unistd::Pid;
use perf_stats::*;
use std::collections::HashMap;
use std::env;
use std::io::Read;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Mutex;
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use string_utils::*;

const DEJAVU: Font = Font::External {
    name: "DEJAVU",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono.ttf"),
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Global variables.

static REMOTE: AtomicBool = AtomicBool::new(false);
static USING_SETUP: AtomicBool = AtomicBool::new(false);
static CLEANED_UP: AtomicBool = AtomicBool::new(false);

static REMOTE_SERVER_ID: AtomicUsize = AtomicUsize::new(0);
static SETUP_PID: AtomicUsize = AtomicUsize::new(0);

static PROCESSING_REQUEST: AtomicBool = AtomicBool::new(false);

static DONE: AtomicBool = AtomicBool::new(false);

lazy_static! {
    static ref HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    static ref USER_REQUEST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    static ref SERVER_REPLY_TEXT: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    static ref SERVER_REPLY_SVG: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    static ref CONFIG_FILE: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn _truncate(s: &str) -> String {
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

fn cleanup() {
    if !CLEANED_UP.load(SeqCst) {
        CLEANED_UP.store(true, SeqCst);
        if REMOTE.load(SeqCst) {
            if USING_SETUP.load(SeqCst) {
                kill(Pid::from_raw(SETUP_PID.load(SeqCst) as i32), SIGINT_nix).unwrap();
            }
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

// Redirect SIGINT interrupts to the function "handler".  There may be issues with reliablity,
// since a CTRL-C could happen at any point, including in the memory manager.

fn install_signal_handler() -> Result<(), Error> {
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

extern "C" fn exit_handler() {
    cleanup();
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn enclone_client(t: &Instant) -> Result<(), Box<dyn std::error::Error>> {
    //
    // Set up to catch CTRL-C events.  Parse arguments.

    let _ = install_signal_handler();
    let args: Vec<String> = env::args().collect();
    let mut verbose = false;
    let mut config_name = String::new();
    for i in 1..args.len() {
        let arg = &args[i];
        if arg == "VERBOSE" {
            verbose = true;
        } else if arg.starts_with("COM=") {
            config_name = arg.after("COM=").to_string();
        } else {
            eprintln!(
                "\nCurrently the only allowed arguments are COM=x where x is a \
                configuration name, and VERBOSE.\n"
            );
            std::process::exit(1);
        }
    }
    unsafe {
        atexit(exit_handler);
    }

    // Announce.

    if !verbose {
        println!(
            "\nHi! You are using the experimental enclone GUI client.  If you get an error \
            message or a\nGUI window does not pop up, please rerun the command with the added \
            argument VERBOSE,\nand then ask for help."
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
        let env_var = format!("ENCLONE_COM_{}", config_name);
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

    // Loop through random ports until we get one that works.

    loop {
        let start = SystemTime::now();
        let since_the_epoch = start
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards");
        let nanos = since_the_epoch.subsec_nanos() as u64;
        let port: u16 = (nanos % 65536) as u16;

        // let port: u16 = rng.gen();

        if port < 1024 {
            continue;
        }

        // Attempt to fork the server.

        let mut remote = false;
        let server_process;
        let mut local_host = "127.0.0.1".to_string();
        if config.contains_key("REMOTE_HOST")
            || config.contains_key("REMOTE_IP")
            || config.contains_key("REMOTE_BIN")
        {
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
            }
            remote = true;
            REMOTE.store(true, SeqCst);
        }
        if config.contains_key("REMOTE_SETUP") {
            if !remote {
                eprintln!("\nYou specified REMOTE_SETUP but not REMOTE_HOST.\n");
                std::process::exit(1);
            }
        }
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
            "time spent waiting to read bytes from remote server = {:.1} seconds",
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

        // Form local URL.

        let url = format!("http://{}:{}", local_host, port);

        // Connect to client.

        if verbose {
            println!("connecting to {}", url);
        }
        let client = AnalyzerClient::connect(url).await;
        if client.is_err() {
            eprintln!("\nconnection failed with error\n{:?}\n", client);
            cleanup();
            std::process::exit(1);
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
                        let request = tonic::Request::new(ClonotypeRequest {
                            clonotype_number: n as u32,
                        });
                        let response = client.get_clonotype(request).await;
                        if response.is_err() {
                            eprintln!("\nclonotype request failed\n");
                            std::process::exit(1);
                        }
                        let response = response.unwrap();
                        let r = response.into_inner();
                        output = r.table.clone();
                    } else {
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

        // Launch GUI.  We set the size to a reasonable minimum.

        let mut settings = Settings::default();
        let mut window_settings = iced::window::Settings::default();
        window_settings.size = (1100 as u32, 1060 as u32);
        settings.window = window_settings;
        let _ = EncloneVisual::run(settings);
        cleanup();
        return Ok(());
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
struct EncloneVisual {
    scroll: scrollable::State,
    input: text_input::State,
    input_value: String,
    output_value: String,
    svg_value: String,
    button: button::State,

    open_state: button::State,
    modal_state: modal::State<ModalState>,
    last_message: Option<Message>,
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
    OpenModal,
    CloseModal,
    CancelButtonPressed,
}

#[derive(Default)]
struct ModalState {
    cancel_state: button::State,
}

impl Sandbox for EncloneVisual {
    type Message = Message;

    fn new() -> Self {
        EncloneVisual::default()
    }

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::OpenModal => self.modal_state.show(true),
            Message::CloseModal => self.modal_state.show(false),
            Message::CancelButtonPressed => self.modal_state.show(false),
            Message::InputChanged(ref value) => self.input_value = value.to_string(),
            Message::ButtonPressed => {
                // Need to figure what to do if we are already processing a request, for example
                // if the user pushes the button twice or enters a second command and pushes the
                // button before the first one has completed.  For now, do nothing.

                if !PROCESSING_REQUEST.load(SeqCst) {
                    let t = Instant::now();
                    USER_REQUEST.lock().unwrap().clear();
                    let mut plus = format!("{}", self.input_value.clone());
                    if CONFIG_FILE.lock().unwrap().len() > 0 {
                        let config_file = CONFIG_FILE.lock().unwrap()[0].clone();
                        plus = format!("{} CONFIG={}", plus, config_file);
                    }
                    USER_REQUEST.lock().unwrap().push(plus);
                    PROCESSING_REQUEST.store(true, SeqCst);
                    while PROCESSING_REQUEST.load(SeqCst) {
                        thread::sleep(Duration::from_millis(10));
                    }
                    let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                    if reply_text.contains("enclone failed") {
                        reply_text =
                            format!("enclone failed{}", reply_text.after("enclone failed"));
                    }
                    let reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                    self.output_value = reply_text.to_string();
                    self.svg_value = reply_svg.to_string();
                    println!(
                        "time used processing command = {:.1} seconds\n",
                        elapsed(&t)
                    );
                }
            }
        }
        self.last_message = Some(message)
    }

    fn view(&mut self) -> Element<Message> {
        let text_input = TextInput::new(
            &mut self.input,
            "",
            &self.input_value,
            Message::InputChanged,
        )
        .padding(10)
        .size(20);

        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);

        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .style(style::Squeak)
            .push(Text::new(&self.output_value).font(DEJAVU).size(13));

        // Display the SVG.

        let svg = Svg::new(Handle::from_memory(self.svg_value.as_bytes().to_vec()))
            .width(Units(400))
            .height(Units(400));

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .max_width(1500) // this governs the max window width upon manual resizing
            .push(Row::new().spacing(10).align_items(Align::Center).push(
                Button::new(&mut self.open_state, Text::new("Help")).on_press(Message::OpenModal),
            ))
            .push(Row::new().spacing(10).push(text_input).push(button))
            .push(Row::new().spacing(10).push(svg))
            .push(Rule::horizontal(10))
            .push(
                Row::new()
                    .height(Length::Units(1000)) // Height of scrollable window, maybe??
                    .align_items(Align::Center)
                    .push(scrollable),
            );

        use iced_aw::style::{
            card::{Style, StyleSheet},
            colors,
        };

        #[derive(Clone, Copy)]
        pub struct Gerbil;

        impl StyleSheet for Gerbil {
            fn active(&self) -> Style {
                Style {
                    background: iced::Background::Color(Color::from_rgb(0.9, 1.0, 0.9)),
                    border_width: 0.0,
                    border_color: iced::Color::from_rgb(1.0, 1.0, 1.0),
                    head_background: iced::Background::Color(Color::from_rgb(0.9, 1.0, 0.9)),
                    head_text_color: colors::WHITE,
                    close_color: colors::WHITE,
                    ..Style::default()
                }
            }
        }

        let style = Gerbil;

        Modal::new(&mut self.modal_state, content, move |state| {
            Card::new(
                Text::new(""),
                Text::new(
                    "Welcome to enclone visual 0.000...0001!\n\n\
                        To use it, type in the box \
                       (see below)\nand then push the Submit button.  Here are the things \
                       that you can type:\n\n\
                        • an enclone command, without the enclone part\n\
                        • an clonotype id (number)\n\
                        • d, for a demo, same as BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui\n\
                        • q to quit\n\n\
                     \n\
                     Major limitations of this version:\n\
                     1. There is no color in the clonotype tables.\n\
                     2. Text in plots does not show up.\n\
                     3. Cutting and pasting from clonotype tables doesn't work.",
                )
                .height(Units(400))
                .vertical_alignment(VerticalAlignment::Center),
            )
            .style(style)
            .foot(
                Row::new().spacing(10).push(
                    Button::new(
                        &mut state.cancel_state,
                        Text::new("Dismiss").horizontal_alignment(HorizontalAlignment::Left),
                    )
                    // .width(Length::Fill)
                    .on_press(Message::CancelButtonPressed),
                ),
            )
            .width(Units(1100))
            .height(Units(1060))
            .on_close(Message::CloseModal)
            .into()
        })
        .backdrop(Message::CloseModal)
        .on_esc(Message::CloseModal)
        .into()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

mod style {

    pub struct Squeak;

    use iced::{scrollable, Color};

    impl scrollable::StyleSheet for Squeak {
        fn active(&self) -> scrollable::Scrollbar {
            scrollable::Scrollbar {
                background: Color::from_rgb(0.75, 0.75, 0.75).into(),
                border_radius: 2.0,
                border_width: 0.0,
                border_color: Color::TRANSPARENT,
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    border_radius: 2.0,
                    border_width: 0.0,
                    border_color: Color::TRANSPARENT,
                },
            }
        }

        fn hovered(&self) -> scrollable::Scrollbar {
            let active = self.active();
            scrollable::Scrollbar {
                background: Color {
                    a: 0.5,
                    ..Color::from_rgb(0.0, 0.0, 0.0)
                }
                .into(),
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    ..active.scroller
                },
                ..active
            }
        }

        fn dragging(&self) -> scrollable::Scrollbar {
            let hovered = self.hovered();
            scrollable::Scrollbar {
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    ..hovered.scroller
                },
                ..hovered
            }
        }
    }
}
