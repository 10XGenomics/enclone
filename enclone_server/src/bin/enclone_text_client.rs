// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a text client, which is useless except for debugging and experimentation.

use enclone_server::proto::{
    analyzer_client::AnalyzerClient,
    EncloneRequest,
};
use std::io::{self, BufRead, Write};

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
    let mut client = AnalyzerClient::connect("http://127.0.0.1:7000").await?;

    loop {
        print!("\nenter an enclone command, without the enclone part, or type q to quit\n% ");
        std::io::stdout().flush().unwrap();
        let stdin = io::stdin();
        let line = stdin.lock().lines().next().unwrap().unwrap();
        if line == "q" {
            println!("");
            break;
        }

        let request = tonic::Request::new(EncloneRequest {
            // args: "BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui".into(),
            args: line,
        });
    
        let response = client.enclone(request).await?;
    
        let r = response.into_inner();
        println!("\nargs = {}", r.args);
        println!("\nplot = {}", truncate(&r.plot));
        println!("\ntable = {}", truncate(&r.table));
    }
    
    Ok(())
}
