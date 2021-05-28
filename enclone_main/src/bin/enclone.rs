// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_main::enclone_client::enclone_client;
use enclone_main::main_enclone::{main_enclone, MainEncloneOutput};
use enclone_main::proto::{
    analyzer_client::AnalyzerClient,
    analyzer_server::{Analyzer, AnalyzerServer},
    ClonotypeRequest, ClonotypeResponse, EncloneRequest, EncloneResponse, Unit,
};
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
        let fields = &req.args.split(' ').collect::<Vec<&str>>();
        let mut args = vec!["enclone".to_string()];
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
        args.push("NOPRINT".to_string());
        args.push("NOPAGER".to_string());
        args.push("PLAIN".to_string()); // until colored text can be rendered
        args.push("INTERNAL".to_string());
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
            table.truncate(100);
            response = EncloneResponse {
                args: req.args,
                plot: enclone_output.svgs[0].clone(),
                table: table.join("\n"),
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
            return Err(Status::new(Code::Internal, "clonotype id too large"));
        }

        // Send back the clonotype picture.
        let table = &enclone_output.pics[id];
        Ok(Response::new(ClonotypeResponse {
            table: table.to_string(),
        }))
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t = Instant::now();
    PrettyTrace::new().on();
    let mut args: Vec<String> = env::args().collect();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Client run of enclone.

    for i in 0..args.len() {
        if args[i].starts_with("COM=") {
            enclone_client(&t).await?;
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Standard run of enclone.

    if args.len() < 2 || args[1] != "SERVER" {
        let res = main_enclone(&mut args).await;
        if res.is_err() {
            eprintln!("{}", res.unwrap_err());
            std::process::exit(1);
        }
        std::process::exit(0);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Server run of enclone.

    let mut ip_port = "127.0.0.1:7000".to_string();
    if args.len() > 2 {
        ip_port = args[2].clone();
    }

    // Start server

    let addr = ip_port;
    let enclone_command = Arc::new(Mutex::new("".to_string()));
    let enclone_output = Arc::new(Mutex::new(MainEncloneOutput::default()));
    let analyzer = EncloneAnalyzer {
        enclone_command: Arc::clone(&enclone_command),
        enclone_output: Arc::clone(&enclone_output),
    };

    let listener = TcpListener::bind(addr).await?;
    let local_addr = listener.local_addr()?;

    // thread waits to print PORT for client until we can connect to our own endpoints

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
