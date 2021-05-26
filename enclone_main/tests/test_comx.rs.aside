// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Partially test enclone server/client interaction.  This uses the actual server and the
// (useless) text client.
//
// To test this alone, use "cargo test --test test_comx -- --nocapture".
//
// ./test runs this and many other tests.

use io_utils::*;
use nix::sys::signal::{kill, SIGINT};
use nix::unistd::Pid;
use pretty_trace::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use std::process::{Command, Stdio};
use std::thread;
use std::time;
use string_utils::*;

#[test]
fn test_comx() {
    PrettyTrace::new().on();

    // Define commmands to be tested.

    let requests = ["BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui", "7", "q"];

    // Define sleep time in milliseconds.

    let pause = 100;

    // Start the server.

    let server_process = match Command::new("enclone_server")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn server: {}", why),
        Ok(server_process) => server_process,
    };

    let server_process_id = server_process.id();

    // Wait until server has printed something.

    let mut buffer = [0; 50];
    let mut server_stdout = server_process.stdout.unwrap();
    server_stdout.read(&mut buffer).unwrap();
    thread::sleep(time::Duration::from_millis(pause));

    // Look at stderr.

    let mut ebuffer = [0; 200];
    let mut server_stderr = server_process.stderr.unwrap();
    server_stderr.read(&mut ebuffer).unwrap();
    let emsg = strme(&ebuffer);
    if emsg.contains("already in use") {
        eprintln!("server failed, error = {}", emsg);
        std::process::exit(1);
    }

    // Fork client.

    let client_process = match Command::new("enclone_text_client")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn client: {}", why),
        Ok(client_process) => client_process,
    };

    let client_process_id = client_process.id();

    // Wait until client has printed something.

    let mut client_buffer = [0; 224]; // much larger than needed
    let mut client_stdout = client_process.stdout.unwrap();
    client_stdout.read(&mut client_buffer).unwrap();

    // Send enclone_client the requests.

    let mut client_input = client_process.stdin.unwrap();
    for x in requests.iter() {
        let msg = format!("{}\n", x);
        client_input.write_all(msg.as_bytes()).unwrap();
        thread::sleep(time::Duration::from_millis(pause));
    }

    // Read the rest of what client printed.

    let mut rest = String::new();
    client_stdout.read_to_string(&mut rest).unwrap();

    // Form total message from client.

    let total = format!("{}{}", strme(&client_buffer), rest);
    let save_output = false;
    if save_output {
        let mut f = open_for_write_new!["test_output"];
        fwrite!(f, "{}", total);
    }

    // Kill processes.  If this step is omitted, then the port is in some sense "hung" after
    // the test completes.

    kill(Pid::from_raw(server_process_id as i32), SIGINT).unwrap();
    kill(Pid::from_raw(client_process_id as i32), SIGINT).unwrap();

    // Verify that client output is correct.

    let total_control = include_str!["test_output"];
    if total != total_control {
        println!("expected output =\n{}", total);
        println!("actual output =\n{}", total_control);
    }
    assert_eq!(total, total_control);
}
