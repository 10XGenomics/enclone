// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::proto::{
    analyzer_client::AnalyzerClient,
    analyzer_server::{Analyzer, AnalyzerServer},
    ClonotypeRequest, ClonotypeResponse, EncloneRequest, EncloneResponse, Unit,
};
use enclone_core::combine_group_pics::*;
use enclone_core::parse_bsv;
use enclone_main::main_enclone::{main_enclone, MainEncloneOutput};
use itertools::Itertools;
use log::{error, warn};
use pretty_trace::*;
use std::env;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use tokio::net::TcpListener;
use tokio_stream::wrappers::TcpListenerStream;
use tonic::{transport::Server, Code, Request, Response, Status};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub struct EncloneAnalyzer {
    enclone_command: Arc<Mutex<String>>,
    enclone_output: Arc<Mutex<MainEncloneOutput>>, // Caches result from main_enclone.
}

#[tonic::async_trait]
impl Analyzer for EncloneAnalyzer {
    async fn ping(&self, _request: Request<Unit>) -> Result<Response<Unit>, Status> {
        Ok(Response::new(Unit {}))
    }

    async fn enclone(
        &self,
        request: Request<EncloneRequest>,
    ) -> Result<Response<EncloneResponse>, Status> {
        // TODO: Actually parse the arguments etc
        let req: EncloneRequest = request.into_inner();

        // Override the output file
        let fields = parse_bsv(&req.args);
        let mut args = Vec::<String>::new();
        let mut server_debug = false;
        for j in 0..fields.len() {
            if fields[j].len() > 0 {
                if fields[j] == "SERVER_DEBUG" {
                    server_debug = true;
                } else {
                    args.push(fields[j].to_string());
                }
            }
        }
        args.push("NOPRINTX".to_string());
        args.push("NOPAGER".to_string());
        args.push("PLAIN".to_string()); // until colored text can be rendered
        eprintln!("Running enclone:\n  {}", args.join(" "));
        let result = main_enclone(&args).await;
        if result.is_err() {
            let err_msg = format!("{}", result.unwrap_err().to_string());
            let mut msg = format!("enclone failed, here is the error message:\n{}\n", err_msg);
            if server_debug {
                msg += &mut format!(
                    "The arguments provided to the server were\n{}.\n",
                    args.iter().format(" ")
                );
            }
            let response = EncloneResponse {
                args: req.args,
                plot: String::new(),
                table: msg,
            };
            return Ok(Response::new(response));
        }
        let output = result.unwrap();
        eprintln!("Enclone done, updating in-memory cache");
        // Update stored command
        {
            let mut enclone_command = self.enclone_command.lock().unwrap();
            *enclone_command = req.args.clone();
        }
        // Update stored result
        let response;
        {
            let mut enclone_output = self.enclone_output.lock().unwrap();
            *enclone_output = output;
            let mut table = enclone_output.pics.clone();
            let mut widths = enclone_output.last_widths.clone();
            if table.len() > 100 {
                table.truncate(100);
                widths.truncate(100);
            }
            let table_string = combine_group_pics(
                &table,
                &widths,
                enclone_output.noprint,
                enclone_output.noprintx,
                enclone_output.html,
                enclone_output.ngroup,
                enclone_output.pretty,
            );
            let mut plot = String::new();
            if enclone_output.svgs.len() > 0 {
                plot = enclone_output.svgs[0].clone();
            }
            response = EncloneResponse {
                args: req.args,
                plot: plot,
                table: table_string,
            };
            if server_debug {
                println!("sending response as follows:");
                println!("args = {}", response.args);
                println!("plot = {}", response.plot);
                println!("table = {}", response.table);
            }
        }
        if server_debug {
            println!("returning response");
        }
        Ok(Response::new(response))
    }

    async fn get_clonotype(
        &self,
        request: Request<ClonotypeRequest>,
    ) -> Result<Response<ClonotypeResponse>, Status> {
        let req: ClonotypeRequest = request.into_inner();
        let id = req.clonotype_number as usize;
        let enclone_output = self.enclone_output.lock().unwrap();
        if id >= enclone_output.pics.len() {
            return Err(Status::new(Code::Internal, "group id too large"));
        }

        // Send back the clonotype picture.
        let table = &enclone_output.pics[id];
        Ok(Response::new(ClonotypeResponse {
            table: table.to_string(),
        }))
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn enclone_server() -> Result<(), Box<dyn std::error::Error>> {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut ip_port = "127.0.0.1:7000".to_string();
    if args.len() > 2 {
        ip_port = args[2].clone();
    }

    // Force exit after 24 hours.  Although the client is supposed to kill the server, sometimes
    // this doesn't work for users (for unclear reasons).  This is the fallback.  A more
    // sophisticated version would wait for a specified period of time after the last activity,
    // and also communicate with the client before exiting.

    tokio::spawn(async move {
        std::thread::sleep(Duration::from_secs(60 * 60 * 24));
        std::process::exit(0);
    });

    // Start server.

    let addr = ip_port;
    let enclone_command = Arc::new(Mutex::new("".to_string()));
    let enclone_output = Arc::new(Mutex::new(MainEncloneOutput::default()));
    let analyzer = EncloneAnalyzer {
        enclone_command: Arc::clone(&enclone_command),
        enclone_output: Arc::clone(&enclone_output),
    };

    let listener = TcpListener::bind(addr).await?;
    let local_addr = listener.local_addr()?;

    // Thread waits to print PORT for client until we can connect to our own endpoints.

    tokio::spawn(async move {
        let dest = format!("http://{}", local_addr);
        let tick = Instant::now();
        loop {
            tokio::time::sleep(Duration::from_secs(1)).await;
            match AnalyzerClient::connect(dest.clone()).await {
                Ok(mut client) => match client.ping(Unit {}).await {
                    Ok(_) => {
                        println!("For debugging:");
                        println!(
                            "  grpcurl -plaintext -import-path ./enclone \
                             -proto ./enclone/server.proto 127.0.0.1:{}",
                            local_addr.port()
                        );
                        println!("To run the client (in another terminal window):");
                        println!("  cd enclone_client; yarn start");
                        return;
                    }
                    Err(e) => warn!("failed to ping, ({:?}), reattempting in 1s", e),
                },
                Err(e) => warn!("failed to connect ({:?}), reattempting in 1s", e),
            }
            if (Instant::now() - tick).as_secs() > 11 {
                error!("Failed to initialize gRPC service, exiting...");
                std::process::exit(1);
            }
        }
    });

    eprintln!("I am process {}.", std::process::id());
    eprintln!("enclone version = {}", env!("CARGO_PKG_VERSION"));
    eprintln!("Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    eprintln!("Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    eprintln!("Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    eprintln!("Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    Server::builder()
        .add_service(AnalyzerServer::new(analyzer))
        .serve_with_incoming(TcpListenerStream::new(listener))
        .await?;

    Ok(())
}
