// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone client.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use enclone_com::typed_com::*;
use pretty_trace::*;
use std::error::Error;
use std::io::{self, BufRead, Write};
use string_utils::*;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    PrettyTrace::new().on();
    println!(
        "\nHello, I am the goofball enclone client.  Now that you've started me, you need\n\
        to start enclone in a separate window, with the argument TOY_COM and enough other\n\
        arguments to run enclone (e.g. providing a dataset).\n"
    );

    // Initiate communication.

    let addr = "127.0.0.1:8080";
    let listener = TcpListener::bind(&addr).await?;
    let (mut socket, _) = listener.accept().await?;
    let mut buf = vec![0; 1024];

    // Loop forever.

    loop {
        // Ask for a request from the user.

        print!("clonotype number or PLOT_BY_ISOTYPE or q to quit? ");
        std::io::stdout().flush().unwrap();
        let stdin = io::stdin();
        let line = stdin.lock().lines().next().unwrap().unwrap();
        if line == "q" {
            println!("");
            std::process::exit(0);
        }
        if line != "PLOT_BY_ISOTYPE" && !line.parse::<usize>().is_ok() {
            println!("That doesn't make sense.\n");
            continue;
        }

        // Send the request to the server.

        let id = line.as_bytes();
        let msg = pack_object(ENCLONE_SUBCHANNEL, "request", &id);
        socket
            .write_all(&msg)
            .await
            .expect("client failed to write data to socket");

        // Get the response.

        buf.clear();
        loop {
            let mut bufbit = vec![0; 1024];
            let n = socket
                .read(&mut bufbit)
                .await
                .expect("client failed to read data from socket");
            buf.append(&mut bufbit[0..n].to_vec());
            if n < 1024 {
                break;
            }
        }
        let n = buf.len();

        // Unpack the response.

        let mut id = 0_u64;
        let mut type_name = Vec::<u8>::new();
        let mut body = Vec::<u8>::new();
        unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);
        if !String::from_utf8(body.clone()).is_ok() {
            println!("received body but it's not UTF-8");
            std::process::exit(1);
        }

        // Check for error.

        if type_name == b"error" {
            print!("{}", strme(&body));
            continue;
        }

        // Check for the expected response.

        if type_name == b"colored-text" {
            println!("\n{}", strme(&body));
            continue;
        }
        if type_name == b"svg" {
            println!("\n{}", strme(&body));
            continue;
        }

        // Otherwise it's an internal error.

        println!(
            "response from server has type {}, which is unknown",
            strme(&type_name)
        );
        std::process::exit(1);
    }
}
