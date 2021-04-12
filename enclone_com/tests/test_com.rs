// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test enclone server/client interaction.
//
// To test this alone, use "cargo test --test test_com -- --nocapture".
//
// ./test runs this and many other tests.

use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use std::process::{Command, Stdio};
use std::thread;
use std::time;
use string_utils::*;

#[test]
fn test_com() {
    PrettyTrace::new().on();

    // Define server command.

    let args = ["BCR=123085", "TOY_COM"];

    // Define requests that would be given to enclone_client.

    let requests = ["7", "100", "PLOT_BY_ISOTYPE", "q"];

    // Define sleep time in milliseconds.

    let pause = 100;

    // Fork client.

    let client_process = match Command::new("enclone_client")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn client: {}", why),
        Ok(client_process) => client_process,
    };

    // Wait until client has printed something.

    let mut client_buffer = [0; 50];
    let mut client_stdout = client_process.stdout.unwrap();
    client_stdout.read(&mut client_buffer).unwrap();

    // Start the server.

    let server_process = match Command::new("enclone")
        .args(&args)
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn server: {}", why),
        Ok(server_process) => server_process,
    };

    // Wait until server has printed something.

    let mut buffer = [0; 50];
    let mut server_stdout = server_process.stdout.unwrap();
    server_stdout.read(&mut buffer).unwrap();
    thread::sleep(time::Duration::from_millis(pause));

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

    // Verify that client output is correct.

    let total_control = include_str!["test_output"];
    assert_eq!(total, total_control);
}
