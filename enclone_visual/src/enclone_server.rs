// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::proto::{
    analyzer_client::AnalyzerClient,
    analyzer_server::{Analyzer, AnalyzerServer},
    ClonotypeRequest, ClonotypeResponse, EncloneRequest, EncloneResponse, Unit,
};
use enclone_core::combine_group_pics::*;
use enclone_core::parse_bsv;
use enclone_main::main_enclone::*;
use enclone_main::stop::*;
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
    enclone_state: Arc<Mutex<EncloneState>>, // caches enclone state
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
        let setup = main_enclone_setup(&args);
        if setup.is_err() {
            let err_msg = format!("{}", setup.err().unwrap());
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
        let setup = setup.unwrap();
        if setup.tall.is_none() {
            let response = EncloneResponse {
                args: req.args,
                plot: String::new(),
                table: String::new(),
            };
            return Ok(Response::new(response));
        }

        // Check for change to setup that could change intermediates.  We are very conservative
        // about this, and only allow changes to:
        // * start_time
        // * clono_filt_opt
        // * plot_opt.
        // More exceptions could be added.

        let mut changed = false;
        {
            // last_setup must be scoped or enclone VIS will mysteriously fail
            let last_setup = &self.enclone_state.lock().unwrap().inter.setup;
            if setup.ctl.perf_opt != last_setup.ctl.perf_opt {
                changed = true;
            }
            if setup.ctl.gen_opt != last_setup.ctl.gen_opt {
                changed = true;
            }
            if setup.ctl.pretty != last_setup.ctl.pretty {
                changed = true;
            }
            if setup.ctl.silent != last_setup.ctl.silent {
                changed = true;
            }
            if setup.ctl.force != last_setup.ctl.force {
                changed = true;
            }
            if setup.ctl.debug_table_printing != last_setup.ctl.debug_table_printing {
                changed = true;
            }
            if setup.ctl.merge_all_impropers != last_setup.ctl.merge_all_impropers {
                changed = true;
            }
            if setup.ctl.heur != last_setup.ctl.heur {
                changed = true;
            }
            if setup.ctl.origin_info != last_setup.ctl.origin_info {
                changed = true;
            }
            if setup.ctl.clono_filt_opt_def != last_setup.ctl.clono_filt_opt_def {
                changed = true;
            }
            if setup.ctl.allele_alg_opt != last_setup.ctl.allele_alg_opt {
                changed = true;
            }
            if setup.ctl.allele_print_opt != last_setup.ctl.allele_print_opt {
                changed = true;
            }
            if setup.ctl.join_alg_opt != last_setup.ctl.join_alg_opt {
                changed = true;
            }
            if setup.ctl.clono_print_opt != last_setup.ctl.clono_print_opt {
                changed = true;
            }
            if setup.ctl.clono_group_opt != last_setup.ctl.clono_group_opt {
                changed = true;
            }
            if setup.ctl.parseable_opt != last_setup.ctl.parseable_opt {
                changed = true;
            }
            if setup.ctl.pathlist != last_setup.ctl.pathlist {
                changed = true;
            }
            if setup.ctl.last_modified != last_setup.ctl.last_modified {
                changed = true;
            }
        }

        // Now proceed with the computation.

        let result;
        if !changed {
            result = main_enclone_stop(EncloneIntermediates {
                setup: setup,
                ex: self.enclone_state.lock().unwrap().inter.ex.clone(),
            });
        } else {
            let inter = main_enclone_start(setup);
            if inter.is_err() {
                let err_msg = format!("{}", inter.err().unwrap());
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
            let inter = inter.unwrap();
            if inter.setup.tall.is_none() {
                let response = EncloneResponse {
                    args: req.args,
                    plot: String::new(),
                    table: String::new(),
                };
                return Ok(Response::new(response));
            }
            result = main_enclone_stop(inter);
        }
        if result.is_err() {
            let err_msg = format!("{}", result.err().unwrap());
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
            let mut enclone_state = self.enclone_state.lock().unwrap();
            *enclone_state = output;
            let mut table = enclone_state.outs.pics.clone();
            let mut widths = enclone_state.outs.last_widths.clone();
            if table.len() > 100 {
                table.truncate(100);
                widths.truncate(100);
            }
            let table_string = combine_group_pics(
                &table,
                &widths,
                enclone_state.outs.noprint,
                enclone_state.outs.noprintx,
                enclone_state.outs.html,
                enclone_state.outs.ngroup,
                enclone_state.outs.pretty,
            );
            let mut plot = String::new();
            if enclone_state.outs.svgs.len() > 0 {
                plot = enclone_state.outs.svgs[0].clone();
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
        let enclone_state = self.enclone_state.lock().unwrap();
        if id >= enclone_state.outs.pics.len() {
            return Err(Status::new(Code::Internal, "group id too large"));
        }

        // Send back the clonotype picture.
        let table = &enclone_state.outs.pics[id];
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
    let enclone_state = Arc::new(Mutex::new(EncloneState::default()));
    let analyzer = EncloneAnalyzer {
        enclone_command: Arc::clone(&enclone_command),
        enclone_state: Arc::clone(&enclone_state),
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
                        println!("using port {}", local_addr.port());
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
