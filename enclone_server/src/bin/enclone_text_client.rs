// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a text client, which is useless except for debugging and experimentation.
// Accepts one argument, the IP address of the server that enclone_server is running on.

use enclone_server::proto::{analyzer_client::AnalyzerClient, ClonotypeRequest, EncloneRequest};
use pretty_trace::*;
use std::env;
use std::io::{self, BufRead, Write};
use string_utils::*;

fn truncate(s: &str) -> String {
    const MAX_LINES: usize = 10;
    let mut t = String::new();
    let mut extra = 0;
    for (i, line) in s.lines().enumerate() {
        if i < MAX_LINES {
            t += &mut format!("{}\n", line);
        } else {
            extra += 1;
        }
    }
    if extra > 0 {
        t += &mut format!("(+ {} more lines)", extra);
    }
    t
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut url = "http://127.0.0.1:7000".to_string();
    if args.len() > 1 {
        url = format!("http://{}", args[1]);
    }
    let mut client = AnalyzerClient::connect(url).await?;

    loop {
        print!(
            "\nenter one of the following:\n\
            • an enclone command, without the enclone part\n\
            • an clonotype id (number)\n\
            • q to quit\n\n% "
        );
        std::io::stdout().flush().unwrap();
        let stdin = io::stdin();
        let line = stdin.lock().lines().next().unwrap().unwrap();
        if line == "q" {
            println!("");
            break;
        }

        if line.parse::<usize>().is_ok() {
            let n = line.force_usize();
            let request = tonic::Request::new(ClonotypeRequest {
                clonotype_number: n as u32,
            });
            let response = client.get_clonotype(request).await?;
            let r = response.into_inner();
            println!("\ntable = {}", truncate(&r.table));
        } else {
            let request = tonic::Request::new(EncloneRequest { args: line });
            let response = client.enclone(request).await?;
            let r = response.into_inner();
            println!("\nargs = {}", r.args);
            println!("\nplot = {}", truncate(&r.plot));
            println!("\ntable = {}", truncate(&r.table));
        }
    }

    Ok(())
}
