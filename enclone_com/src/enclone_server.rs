// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone server.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use crate::typed_com::*;
use enclone_core::*;
use string_utils::*;
use tokio::io::AsyncWriteExt;

use tokio::io;

use bytes::Bytes;
use futures::{future, Sink, SinkExt, Stream, StreamExt};
use std::{error::Error, net::SocketAddr};
use tokio::net::TcpStream;
use tokio_util::codec::{BytesCodec, FramedRead, FramedWrite};

use std::thread;

use std::cmp::min;

pub async fn connect(
    addr: &SocketAddr,
    mut stdin: impl Stream<Item = Result<Bytes, io::Error>> + Unpin,
    mut stdout: impl Sink<Bytes, Error = io::Error> + Unpin,
) -> Result<(), Box<dyn Error>> {
    let stream = TcpStream::connect(addr).await;
    if stream.is_err() {
        println!("The connect function in enclone server failed on calling TcpStream::connect.");
        println!("error = {:?}", stream.err());
        println!(
            "\nThe most likely explanation for this is that you failed to start enclone_client\n\
            before starting enclone.\n"
        );
        std::process::exit(1);
    }
    let mut stream = stream.unwrap();
    let (r, w) = stream.split();
    let mut sink = FramedWrite::new(w, BytesCodec::new());
    // filter map Result<BytesMut, Error> stream into just a Bytes stream to match stdout Sink
    // on the event of an Error, log the error and end the stream
    let mut stream = FramedRead::new(r, BytesCodec::new())
        .filter_map(|i| match i {
            //BytesMut into Bytes
            Ok(i) => future::ready(Some(i.freeze())),
            Err(e) => {
                println!("failed to read from socket; error={}", e);
                future::ready(None)
            }
        })
        .map(Ok);

    match future::join(sink.send_all(&mut stdin), stdout.send_all(&mut stream)).await {
        (Err(e), _) | (_, Err(e)) => Err(e.into()),
        _ => Ok(()),
    }
}

pub async fn enclone_server() -> Result<(), Box<dyn Error>> {
    // Fixed address for now.

    let addr_raw = "127.0.0.1:8080";

    // Create socket.

    println!("entering enclone_server"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let addr = addr_raw.parse::<SocketAddr>()?;

    let mut stream = TcpStream::connect(addr).await.unwrap();

    // Initiate communication.

    // let listener = TcpListener::bind(&addr).await?;
    println!("Listening on: {}", addr);

    let mut buf = vec![0; 1024]; // not sure why we can't just use Vec::<u8>::new()

    // Loop forever.

    println!("start loop"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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

        println!("received message of length {}", n); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        // Unpack the message.

        let mut id = 0_u64;
        let mut type_name = Vec::<u8>::new();
        let mut body = Vec::<u8>::new();
        unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);
        println!("message unpacked"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        println!("type name = {}", strme(&type_name));

        // Check for request.

        if type_name == b"request" {
            println!("Request received.");
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

                println!("packing clonotype picture");

                let mut chars = Vec::<char>::new();
                for char in pics[id].chars() {
                    chars.push(char);
                }

                let mut s = String::new();
                for i in 0..min(chars.len(), 10000) {
                    s.push(chars[i]);
                }
                let pic_bytes = s.as_bytes().to_vec();

                if !String::from_utf8(pic_bytes.clone()).is_ok() {
                    println!("pic_bytes is not UTF-8");
                    std::process::exit(1);
                }

                msg = pack_object(ENCLONE_SUBCHANNEL, "colored-text", &pic_bytes);
                println!("packed message has length {}", msg.len());
            }

            println!("writing message to stream");

            // let _result = stream.write(&msg).await?;

            let mut start = 0;
            while start < msg.len() {
                let stop = min(start + 1024, msg.len());
                let n = stream.write(&msg[start..stop]).await.unwrap();
                println!("sent {} bytes", n);
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
