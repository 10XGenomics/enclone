// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Toy enclone server.  Not useful, except pedagogically.
//
// based on code in the rust tokio crate

use crate::typed_com::*;
use enclone_core::defs::*;
use enclone_tail::plot::*;
use std::cmp::min;
use std::error::Error;
use string_utils::*;
use tokio::io::AsyncWriteExt;
use tokio::net::TcpStream;
use vdj_ann::refx::RefData;

pub async fn enclone_server(
    ctl: &EncloneControl,
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    pics: &Vec<String>,
) -> Result<(), Box<dyn Error>> {
    // Fixed address for now.

    let addr = "127.0.0.1:8080";

    // Create socket.

    let stream = TcpStream::connect(addr).await;
    if !stream.is_ok() {
        println!(
            "Connection failed.  The most likely explanation is that you did not start \
            enclone_client before\nyou started enclone.  However, in case that's not the problem, \
            here is the error message:\n{:?}\n",
            stream.err()
        );
        std::process::exit(1);
    }
    let mut stream = stream.unwrap();

    // Loop forever.

    let mut buf = vec![0; 1024];
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
        }

        // Unpack the message.

        let mut id = 0_u64;
        let mut type_name = Vec::<u8>::new();
        let mut body = Vec::<u8>::new();
        unpack_message(&mut buf[0..n], &mut id, &mut type_name, &mut body);

        // Process request.

        if type_name == b"request" {
            let bs = strme(&body);
            if bs != "PLOT_BY_ISOTYPE" && !bs.parse::<usize>().is_ok() {
                println!(
                    "received the request \"{}\", which doesn't make sense",
                    strme(&body)
                );
                continue;
            }
            let msg;
            if bs == "PLOT_BY_ISOTYPE" {
                let mut svg = String::new();
                let mut plot_opt = PlotOpt::default();
                plot_opt.plot_by_isotype = true;
                plot_opt.plot_file = "stdout".to_string();
                plot_clonotypes(
                    &ctl,
                    &plot_opt,
                    &refdata,
                    &exacts,
                    &exact_clonotypes,
                    &mut svg,
                );
                let svg_bytes = svg.as_bytes().to_vec();
                msg = pack_object(ENCLONE_SUBCHANNEL, "svg", &svg_bytes);
            } else {
                {
                    let id = strme(&body).force_usize();
                    if id >= pics.len() {
                        println!("clonotype id is too large");
                        continue;
                    }

                    // Send back the clonotype picture.

                    let pic_bytes = pics[id].as_bytes().to_vec();
                    msg = pack_object(ENCLONE_SUBCHANNEL, "colored-text", &pic_bytes);
                }
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
