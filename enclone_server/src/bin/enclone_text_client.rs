// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a text client, which is useless except for debugging and experimentation.
//
// This starts enclone_server, which can be either local or remote.
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
// REMOVE_BIN=...            directory on remote host containing the enclone_server executable

use async_trait::async_trait;
use enclone_core::parse_bsv;
use enclone_server::proto::{analyzer_client::AnalyzerClient, ClonotypeRequest, EncloneRequest};
use iced::{
    button, scrollable, text_input, Align, Button, Column, Container, Element, Length, Row, Rule,
    Sandbox, Scrollable, Settings, Text, TextInput,
};
use itertools::Itertools;
use lazy_static::lazy_static;
use nix::sys::signal::{kill, SIGINT};
use nix::unistd::Pid;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::io::{self, BufRead, Read, Write};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::atomic::Ordering::SeqCst;
use std::sync::Mutex;
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use string_utils::*;
use tonic::transport::Channel;

type Com = Option<AnalyzerClient<Channel>>;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn truncate(s: &str) -> String {
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

static REMOTE: AtomicBool = AtomicBool::new(false);
static REMOTE_SERVER_ID: AtomicUsize = AtomicUsize::new(0);

lazy_static! {
    static ref HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

    // let sep = HAPPENING.lock().unwrap().sep;
    // HAPPENING.lock().unwrap().on = happening.on;

fn cleanup() {
    if REMOTE.load(SeqCst) {
        let host = &HOST.lock().unwrap()[0];
        Command::new("ssh")
            .arg(&host)
            .arg("kill")
            .arg("-9")
            .arg(&format!("{}", REMOTE_SERVER_ID.load(SeqCst)))
            .output()
            .expect("failed to execute ssh to kill");
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Set up initial communications with server.

async fn initialize_com() -> Com {

    // Get configuration.

    let mut configuration = None;
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 {
        let arg = args[1].to_string();
        if !arg.starts_with("COM=") {
            eprintln!(
                "\nCurrently the only allowed argument is COM=x where x is a \
                configuration name.\n"
            );
            std::process::exit(1);
        }
        let config_name = arg.after("COM=");
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

    // let mut rng = rand::thread_rng();
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
        let mut host = String::new();
        if remote {
            host = config["REMOTE_HOST"].clone();
            HOST.lock().unwrap().push(host.clone());
            let ip = &config["REMOTE_IP"];
            let bin = &config["REMOTE_BIN"];
            println!(
                "\nstarting remote server using\nssh {} {}/enclone_server {}:{}",
                host, bin, ip, port
            );
            server_process = Command::new("ssh")
                .arg("-n")
                .arg(&host)
                .arg(&format!("{}/enclone_server", bin))
                .arg(&format!("{}:{}", ip, port))
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn();
            local_host = "localhost".to_string();
        } else {
            let ip = "127.0.0.1";
            server_process = Command::new("enclone_server")
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
        let server_process = server_process.unwrap();
        let server_process_id = server_process.id();

        // Wait until server has printed something.

        let mut buffer = [0; 50];
        let mut server_stdout = server_process.stdout.unwrap();
        server_stdout.read(&mut buffer).unwrap();
        thread::sleep(Duration::from_millis(100));

        // Look at stderr.

        let mut ebuffer = [0; 200];
        let mut server_stderr = server_process.stderr.unwrap();
        server_stderr.read(&mut ebuffer).unwrap();
        let emsg = strme(&ebuffer);
        if emsg.len() > 0 {
            if emsg.contains("already in use") {
                println!("oops, that port is in use, trying a different one");
                continue;
            }
            eprintln!("\nserver says this:\n{}", emsg);
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
            eprintln!("\nrunning setup command = {}", argsp.iter().format(" "));
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
                kill(Pid::from_raw(server_process_id as i32), SIGINT).unwrap();
                cleanup();
                std::process::exit(1);
            }
            let _setup_process = setup_process.unwrap();
            thread::sleep(Duration::from_millis(2000));
        }

        // Form local URL.

        let url = format!("http://{}:{}", local_host, port);

        // Connect to client.

        println!("connecting to {}", url);
        let client = AnalyzerClient::connect(url).await;
        if client.is_err() {
            eprintln!("\nconnection failed with error\n{:?}\n", client);
            cleanup();
            std::process::exit(1);
        }
        println!("connected");
        return Some(client.unwrap());
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

async fn process_command(input: &str, com: &mut Com) -> String {
    let mut line = input.to_string();
    let mut output = String::new();
    if line == "q" {
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
        let response = com.as_mut().unwrap().get_clonotype(request).await;
        if response.is_err() {
            eprintln!("\nclonotype request failed\n");
            std::process::exit(1);
        }
        let response = response.unwrap();
        let r = response.into_inner();
        output = format!("\ntable = {}", truncate(&r.table));
    } else {
        let request = tonic::Request::new(EncloneRequest { args: line });
        let response = com.as_mut().unwrap().enclone(request).await;
        if response.is_err() {
            let left = r###"message: "\n"###;
            let right = r###"\n""###;
            let mut err = format!("{:?}", response);
            if err.contains(&left) && err.after(&left).contains(&right) {
                err = err.between(&left, &right).to_string();
            }
            output = format!("\nThe server is unhappy.  It says:\n{}", err);
        } else {
            let response = response.unwrap();
            let r = response.into_inner();
            output = format!("\nargs = {}", r.args);
            output += &format!("\nplot = {}", truncate(&r.plot));
            output += &format!("\ntable = {}", truncate(&r.table));
        }
    }
    output
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
struct Styling {
    scroll: scrollable::State,
    input: text_input::State,
    input_value: String,
    submitted_input_value: String,
    button: button::State,
    client: Com,
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
}

#[async_trait]
impl Sandbox for Styling {
    type Message = Message;

    async fn new() -> Self {
        let mut sss = Styling::default();
        sss.client = initialize_com().await; // the await is new, maybe not right
        sss
    }

            // Kill server.

            // kill(Pid::from_raw(server_process_id as i32), SIGINT).unwrap();
            // cleanup();


    fn title(&self) -> String {
        String::from("Styling - Iced")
    }

    async fn update(&mut self, message: Message) {
        match message {
            Message::InputChanged(value) => self.input_value = value,
            Message::ButtonPressed => {
                let mut line = self.input_value.clone();
                let mut output = String::new();
                let output = process_command(&line, &mut self.client);
                self.submitted_input_value = output.await; // await is new, maybe wrong
            }
        }
    }

    fn view(&mut self) -> Element<Message> {
        let text_input = TextInput::new(
            &mut self.input,
            "Enter command...",
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
            .push(Text::new(&self.submitted_input_value));

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .max_width(600)
            .push(Row::new().spacing(10).push(text_input).push(button))
            .push(
                Row::new()
                    .spacing(10)
                    .height(Length::Units(100))
                    .align_items(Align::Center)
                    .push(scrollable)
                    .push(Rule::vertical(38)),
            );

        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .center_x()
            .center_y()
            .into()
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    PrettyTrace::new().on();
    Styling::run(Settings::default());
    Ok(())
}
