// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone server.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use crate::typed_com::*;
use enclone_core::*;
use std::env;
use std::error::Error;
use string_utils::*;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;

pub async fn enclone_server() -> Result<(), Box<dyn Error>> {
    // Initiate communication.

    let addr = env::args()
        .nth(1)
        .unwrap_or_else(|| "127.0.0.1:8080".to_string());
    let listener = TcpListener::bind(&addr).await?;
    println!("Listening on: {}", addr);

    // Not sure what this outer loop is doing.

    loop {
        let (mut socket, _) = listener.accept().await?;
        tokio::spawn(async move {
            let mut buf = vec![0; 1024]; // not sure why we can't just use Vec::<u8>::new()

            // Loop forever.

            loop {
                // Wait for message from client.

                let n = socket
                    .read(&mut buf)
                    .await
                    .expect("failed to read data from socket");

                // Unpack the message.

                let mut id = 0_u64;
                let mut type_name = Vec::<u8>::new();
                let mut body = Vec::<u8>::new();
                unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);

                // Get the response.

                let n = socket
                    .read(&mut buf)
                    .await
                    .expect("failed to read data from socket");

                // Unpack the response.

                let mut id = 0_u64;
                let mut type_name = Vec::<u8>::new();
                let mut body = Vec::<u8>::new();
                unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);

                // Check for request.

                if type_name == b"request" {
                    if !strme(&body).parse::<usize>().is_ok() {
                        println!(
                            "received the request \"{}\", which doesn't make sense",
                            strme(&body)
                        );
                        continue;
                    }

                    let msg;

                    {
                        let pics = PICS.lock().unwrap();
                        let id = strme(&body).force_usize();
                        if id >= pics.len() {
                            println!("clonotype id is too large");
                            continue;
                        }

                        // Send back the clonotype picture.

                        msg = pack_object(ENCLONE_SUBCHANNEL, "colored-text", &pics[id].as_bytes());
                    }

                    socket
                        .write_all(&msg)
                        .await
                        .expect("server failed to write data to socket");
                    continue;
                }

                // Otherwise it's an internal error.

                println!(
                    "client has sent type {}, which is unknown",
                    strme(&type_name)
                );
                std::process::exit(1);
            }
        });
    }
}
