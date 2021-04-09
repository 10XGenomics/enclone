// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone client.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use enclone_com::typed_com::*;
use pretty_trace::*;
use std::env;
use std::error::Error;
use std::io::{self, BufRead};
use string_utils::*;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;

use std::thread;
use std::time::Duration;

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    PrettyTrace::new().on();
    // Initiate communication.

    let addr = env::args()
        .nth(1)
        .unwrap_or_else(|| "127.0.0.1:8080".to_string());
    let listener = TcpListener::bind(&addr).await?;
    // println!("Listening on: {}", addr);

    // Not sure what this outer loop is doing.

    loop {
        let (mut socket, _) = listener.accept().await?;
        tokio::spawn(async move {
            let mut buf = vec![0; 1024]; // not sure why we can't just use Vec::<u8>::new()

            // Loop forever.

            println!("");
            loop {
                // Request a clonotype number from the user.

                print!("clonotype number? ");
                use std::io::Write;
                std::io::stdout().flush().unwrap();
                let stdin = io::stdin();
                let line = stdin.lock().lines().next().unwrap().unwrap();
                if !line.parse::<usize>().is_ok() {
                    println!("That doesn't make sense.\n");
                    continue;
                }
                println!("");

                // Send the clonotype number to the server.

                let id = line.as_bytes();
                let msg = pack_object(ENCLONE_SUBCHANNEL, "request", &id);
                socket
                    .write_all(&msg)
                    .await
                    .expect("client failed to write data to socket");

                // Get the response.

                thread::sleep(std::time::Duration::from_millis(1000));

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


                    // Ping server.

                    /*
                    let ping = [0_u8];
                    socket
                        .write_all(&ping)
                        .await
                        .expect("client failed to write data to socket");
                    */


                }
                let n = buf.len();

                /*
                let n = socket
                    .read(&mut buf)
                    .await
                    .expect("client failed to read data from socket");
                */

                // Unpack the response.

                let mut id = 0_u64;
                let mut type_name = Vec::<u8>::new();
                let mut body = Vec::<u8>::new();
                unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);
                if !String::from_utf8(body.clone()).is_ok() {
                    println!("received body but it's not UTF-8");
                    println!("first byte of body is {}", body[0]);
                    std::process::exit(1);
                }

                // Check for error.

                if type_name == b"error" {
                    print!("{}", strme(&body));
                    continue;
                }

                // Check for the expected response.

                if type_name == b"colored-text" {
                    println!("{}", strme(&body));
                    continue;
                }

                // Otherwise it's an internal error.

                println!(
                    "response from server has type {}, which is unknown",
                    strme(&type_name)
                );
                std::process::exit(1);
            }
        });
    }
}
