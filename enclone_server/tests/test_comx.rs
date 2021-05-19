// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Partially test enclone server/client interaction.  This uses the actual server and the
// (useless) text client.
//
// To test this alone, use "cargo test --test test_comx -- --nocapture".
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
    println!("starting test"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Define commmands to be tested.

    let requests = ["BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui", "7", "q"];

    // Define sleep time in milliseconds.

    // XXXXXXXXXXXXXXXXX should be 100 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let pause = 5000;

    // Start the server.

    println!("starting server"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let server_process = match Command::new("enclone_server")
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn server: {}", why),
        Ok(server_process) => server_process,
    };
    println!("done"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Wait until server has printed something.

    let mut buffer = [0; 50];
    let mut server_stdout = server_process.stdout.unwrap();
    server_stdout.read(&mut buffer).unwrap();
    thread::sleep(time::Duration::from_millis(pause));
    println!("server stdout = {}", strme(&buffer)); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Fork client.

    let client_process = match Command::new("enclone_text_client")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn client: {}", why),
        Ok(client_process) => client_process,
    };

    // Wait until client has printed something.

    let mut client_buffer = [0; 224];
    let mut client_stdout = client_process.stdout.unwrap();
    client_stdout.read(&mut client_buffer).unwrap();
    println!("client has printed: {}", strme(&client_buffer)); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Send enclone_client the requests.

    let mut client_input = client_process.stdin.unwrap();
    for x in requests.iter() {
        let msg = format!("{}\n", x);
        print!("msg = {}", msg); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        client_input.write_all(msg.as_bytes()).unwrap();
        println!("sleeping"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        thread::sleep(time::Duration::from_millis(pause));
    }

    // Read the rest of what client printed.

    println!("reading rest of what client printed"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let mut rest = String::new();
    client_stdout.read_to_string(&mut rest).unwrap();

    // Form total message from client.

    let total = format!("{}{}", strme(&client_buffer), rest);
    let save_output = false;
    if save_output {
        let mut f = open_for_write_new!["test_output"];
        fwrite!(f, "{}", total);
    }
    println!("output =\n{}\n", total); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Verify that client output is correct.

    /*
    let total_control = include_str!["test_output"];
    assert_eq!(total, total_control);
    */
}
