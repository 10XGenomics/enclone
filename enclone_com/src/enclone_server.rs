// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone server.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use crate::typed_com::*;
use enclone_core::*;
use string_utils::*;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;

use tokio::io;

use bytes::Bytes;
use futures::{future, Sink, SinkExt, Stream, StreamExt};
use std::{error::Error, net::SocketAddr};
use tokio::net::TcpStream;
use tokio_util::codec::{BytesCodec, FramedRead, FramedWrite};

pub async fn connect(
    addr: &SocketAddr,
    mut stdin: impl Stream<Item = Result<Bytes, io::Error>> + Unpin,
    mut stdout: impl Sink<Bytes, Error = io::Error> + Unpin,
) -> Result<(), Box<dyn Error>> {
    let stream = TcpStream::connect(addr).await;
    println!("back from TcpStream"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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

    println!("at future::join"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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




    /*
    let stdin = FramedRead::new(io::stdin(), BytesCodec::new());
    let stdin = stdin.map(|i| i.map(|bytes| bytes.freeze()));
    let stdout = FramedWrite::new(io::stdout(), BytesCodec::new());
    connect(&addr, stdin, stdout).await?;
    */

    // let stream = TcpStream::connect(addr).await;




    println!("now connected"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Initiate communication.

    let listener = TcpListener::bind(&addr).await?;
    println!("Listening on: {}", addr);

    // Not sure what this outer loop is doing.

    loop {
        let (mut socket, _) = listener.accept().await?;
        tokio::spawn(async move {
            let mut buf = vec![0; 1024]; // not sure why we can't just use Vec::<u8>::new()

            // Loop forever.

            println!("start loop"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            loop {
                // Wait for message from client.

                println!("waiting for message from client"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                let n = socket
                    .read(&mut buf)
                    .await
                    .expect("failed to read data from socket");
                println!("received message of length {}", n); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                // Unpack the message.

                let mut id = 0_u64;
                let mut type_name = Vec::<u8>::new();
                let mut body = Vec::<u8>::new();
                unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);
                println!("message unpacked"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                /*

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

                */

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
