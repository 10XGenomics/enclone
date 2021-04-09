// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone server.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use crate::typed_com::*;
use enclone_core::*;
use string_utils::*;
use tokio::io::AsyncWriteExt;

use std::{error::Error, net::SocketAddr};
use tokio::net::TcpStream;

use std::thread;

use std::cmp::min;

pub async fn enclone_server() -> Result<(), Box<dyn Error>> {
    // Fixed address for now.

    let addr_raw = "127.0.0.1:8080";

    // Create socket.

    println!("entering enclone_server"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let addr = addr_raw.parse::<SocketAddr>()?;

    let mut stream = TcpStream::connect(addr).await.unwrap();

    // Loop forever.

    let mut buf = vec![0; 1024]; // not sure why we can't just use Vec::<u8>::new()
    loop {
        // Wait for message from client.

        let n;
        loop {
            stream.readable().await?;
            let result = stream.try_read(&mut buf);
            if result.is_ok() {
                n = result.unwrap();
                break;
            }
            thread::sleep(std::time::Duration::from_millis(100));
        }

        // Unpack the message.

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

                let pic_bytes = pics[id].as_bytes().to_vec();
                msg = pack_object(ENCLONE_SUBCHANNEL, "colored-text", &pic_bytes);
                println!("packed message has length {}", msg.len());
            }
            let mut start = 0;
            while start < msg.len() {
                let stop = min(start + 1024, msg.len());
                let n = stream.write(&msg[start..stop]).await.unwrap();
                start += n;
            }
            continue;
        }
        println!(
            "client has sent type {}, which is unknown",
            strme(&type_name)
        );
        std::process::exit(1);
    }
}
