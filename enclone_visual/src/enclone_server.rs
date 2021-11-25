// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::proto::{
    analyzer_client::AnalyzerClient,
    analyzer_server::{Analyzer, AnalyzerServer},
    *,
};
use crate::*;
use chrono::prelude::*;
use enclone_core::combine_group_pics::*;
use enclone_core::enclone_structs::*;
use enclone_core::logging::*;
use enclone_core::parse_bsv;
use enclone_core::version_string;
use enclone_main::main_enclone::*;
use enclone_main::stop::*;
use enclone_stuff::start::main_enclone_start;
// use enclone_version::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use io_utils::*;
use itertools::Itertools;
use log::{error, warn};
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::os::unix::fs::PermissionsExt;
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

        let mut fields = parse_bsv(&req.args);
        for j in 0..fields.len() {
            fields[j] = fields[j].replace("\"", "");
        }
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
        let mut g_specified = false;
        for j in 0..args.len() {
            if args[j].starts_with("G=") {
                g_specified = true;
            }
        }
        args.push("SUMMARY".to_string());
        args.push("NOPRINTX".to_string());
        args.push("NOPAGER".to_string());
        args.push("PLAIN".to_string()); // until colored text can be rendered
        args.push("VISUAL".to_string());
        if req.server_logfile.is_some() {
            if enclone_core::logging::SERVER_LOGFILE
                .lock()
                .unwrap()
                .is_empty()
            {
                enclone_core::logging::SERVER_LOGFILE
                    .lock()
                    .unwrap()
                    .push(req.server_logfile.as_ref().unwrap().clone());
            }
        }
        logme(&format!("Running enclone:\n  {}", args.join(" ")));
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
                summary: String::new(),
                metrics: Vec::<String>::new(),
                dataset_names: Vec::<String>::new(),
                table_comp: Vec::<u8>::new(),
                last_widths: Vec::<u32>::new(),
            };
            return Ok(Response::new(response));
        }
        let setup = setup.unwrap();
        if setup.tall.is_none() {
            let response = EncloneResponse {
                args: req.args,
                plot: String::new(),
                table: String::new(),
                summary: String::new(),
                metrics: Vec::<String>::new(),
                dataset_names: Vec::<String>::new(),
                table_comp: Vec::<u8>::new(),
                last_widths: Vec::<u32>::new(),
            };
            return Ok(Response::new(response));
        }

        // Check for change to setup that could change intermediates.  We are very conservative
        // about this, and only allow changes to:
        // * start_time
        // * clono_filt_opt or clono_print_opt
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
                    summary: String::new(),
                    metrics: Vec::<String>::new(),
                    dataset_names: Vec::<String>::new(),
                    table_comp: Vec::<u8>::new(),
                    last_widths: Vec::<u32>::new(),
                };
                return Ok(Response::new(response));
            }
            let inter = inter.unwrap();
            if inter.setup.tall.is_none() {
                let response = EncloneResponse {
                    args: req.args,
                    plot: String::new(),
                    table: String::new(),
                    summary: String::new(),
                    metrics: Vec::<String>::new(),
                    dataset_names: Vec::<String>::new(),
                    table_comp: Vec::<u8>::new(),
                    last_widths: Vec::<u32>::new(),
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
                summary: String::new(),
                metrics: Vec::<String>::new(),
                dataset_names: Vec::<String>::new(),
                table_comp: Vec::<u8>::new(),
                last_widths: Vec::<u32>::new(),
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
            let widths = enclone_state.outs.last_widths.clone();
            if !g_specified {
                const MAX_TABLE: usize = 50;
                if table.len() > MAX_TABLE {
                    table.truncate(MAX_TABLE);
                }
            }
            let mut last_widths = Vec::<u32>::new();
            for i in 0..widths.len() {
                last_widths.push(widths[i] as u32);
            }
            let table_string = combine_group_pics(
                &table,
                &widths,
                enclone_state.outs.parseable_stdouth,
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
            let full_table = enclone_state.outs.pics.clone();
            let serialized = serde_json::to_string(&full_table)
                .unwrap()
                .as_bytes()
                .to_vec();
            let mut e = GzEncoder::new(Vec::new(), Compression::default());
            let _ = e.write_all(&serialized);
            let gzipped = e.finish().unwrap();
            logme(&format!("plot=\n{}", plot));
            response = EncloneResponse {
                args: req.args,
                plot: plot,
                table: table_string,
                summary: enclone_state.outs.summary.clone(),
                metrics: enclone_state.outs.metrics.clone(),
                dataset_names: enclone_state.outs.dataset_names.clone(),
                table_comp: gzipped,
                last_widths: last_widths,
            };
            if server_debug {
                println!("sending response as follows:");
                println!("args = {}", response.args);
                println!("plot = {}", response.plot);
                println!("table = {}", response.table);
                println!("summary = {}", response.summary);
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

    async fn test_user_name(
        &self,
        request: Request<UserNameRequest>,
    ) -> Result<Response<UserNameResponse>, Status> {
        let req: UserNameRequest = request.into_inner();
        let valid = is_user_name_valid(&req.user_name);
        Ok(Response::new(UserNameResponse { value: valid }))
    }

    async fn share_session(
        &self,
        request: Request<SendShareRequest>,
    ) -> Result<Response<SendShareResponse>, Status> {
        let req: SendShareRequest = request.into_inner();
        for recip in req.recipients.iter() {
            let mut bytes = req.content.clone();
            let rbytes = &recip.as_bytes();
            for i in 0..bytes.len() {
                bytes[i] = bytes[i].wrapping_add(rbytes[i % rbytes.len()]);
            }
            let rdir = format!("{}/{}", req.share_dir, recip);
            let perms = std::fs::Permissions::from_mode(0o777);

            // Create directory if needed.

            let dir_exists = path_exists(&rdir);
            if !dir_exists {
                let res = std::fs::create_dir(&rdir);
                if res.is_err() {
                    return Err(Status::new(
                        Code::Internal,
                        "unable to create share directory",
                    ));
                }
            }

            // Set permissions to allow group and world write on the directory.  Note that if the
            // directory already existed, this may not work.  In such cases group and world write
            // should already be enabled, because some other user will have executed this same
            // code to create the directory.

            let res = std::fs::set_permissions(&rdir, perms.clone());
            if !dir_exists && res.is_err() {
                return Err(Status::new(
                    Code::Internal,
                    format!("unable to set permissions on share directory {}", rdir),
                ));
            }

            // Now write the file.

            let mut now = format!("{:?}", Local::now());
            now = now.replace("T", "___");
            now = now.before(".").to_string();
            let filename = format!("{}/{}_{}", rdir, now, req.sender);
            let res = std::fs::write(&filename, &bytes);
            if res.is_err() {
                return Err(Status::new(Code::Internal, "unable to write share file"));
            }
            let res = std::fs::set_permissions(&filename, perms);
            if res.is_err() {
                return Err(Status::new(
                    Code::Internal,
                    format!("unable to set permissions on share file {}", filename),
                ));
            }
        }
        Ok(Response::new(SendShareResponse { ok: true }))
    }

    async fn get_my_shares(
        &self,
        request: Request<GetMySharesRequest>,
    ) -> Result<Response<GetMySharesResponse>, Status> {
        let req: GetMySharesRequest = request.into_inner();
        let dir = &req.share_dir;
        let me_only = req.me_only;
        let me = users::get_current_username();
        if me.is_none() {
            return Err(Status::new(Code::Internal, "unable to determine user name"));
        }
        let me = me.unwrap();
        let me = me.to_string_lossy();
        if !path_exists(&dir) {
            return Err(Status::new(
                Code::Internal,
                "share directory does not exist",
            ));
        }
        let rdir = format!("{}/{}", dir, me);
        if !path_exists(&rdir) {
            let res = std::fs::create_dir(&rdir);
            if res.is_err() {
                return Err(Status::new(
                    Code::Internal,
                    format!("unable to create my share directory {}", rdir),
                ));
            }
        }
        let all = dir_list(&rdir);
        let n = all.len();
        let mut content = vec![Vec::<u8>::new(); n];
        let mut messages = vec![String::new(); n];
        let mut filenames = vec![String::new(); n];
        let rbytes = &me.as_bytes();
        for i in 0..n {
            let filename = format!("{}/{}", rdir, all[i]);
            let f = File::open(&filename);
            if f.is_err() {
                return Err(Status::new(Code::Internal, "unable to open share file"));
            }
            let mut f = f.unwrap();
            let mut bytes = Vec::<u8>::new();
            let res = f.read_to_end(&mut bytes);
            if res.is_err() {
                return Err(Status::new(Code::Internal, "unable to read share file"));
            }
            for i in 0..bytes.len() {
                bytes[i] = bytes[i].wrapping_sub(rbytes[i % rbytes.len()]);
            }
            content[i] = bytes;
            filenames[i] = all[i].clone();
            if !all[i].contains("_") {
                return Err(Status::new(Code::Internal, "malformed file name"));
            }
            let sender = all[i].rev_after("_");
            if me_only && sender != me {
                continue;
            }
            let when = all[i].rev_before("_");
            if !when.contains("___") {
                return Err(Status::new(
                    Code::Internal,
                    format!("ill-formed file name {}", all[i]),
                ));
            }
            let (date, time) = (when.before("___"), when.after("___"));
            let msg = format!("session shared by {} on {} at {}", sender, date, time);
            messages[i] = msg;
        }
        Ok(Response::new(GetMySharesResponse {
            content: content,
            messages: messages,
            filenames: filenames,
        }))
    }

    async fn get_my_cookbooks(
        &self,
        request: Request<GetMyCookbooksRequest>,
    ) -> Result<Response<GetMyCookbooksResponse>, Status> {
        let req: GetMyCookbooksRequest = request.into_inner();
        let dirs = &req.cookbook_dirs;
        let mut cookbooks = Vec::<Vec<u8>>::new();
        for dir in dirs.iter() {
            if path_exists(&*dir) {
                let x = std::fs::read_dir(&*dir);
                if x.is_ok() {
                    let x = x.unwrap();
                    for f in x {
                        let s: String = f.unwrap().file_name().into_string().unwrap();
                        let cb = format!("{}/{}", dir, s);
                        let mut f = File::open(&cb).unwrap();
                        let mut bytes = Vec::<u8>::new();
                        f.read_to_end(&mut bytes).unwrap();
                        cookbooks.push(bytes);
                    }
                }
            }
        }
        Ok(Response::new(GetMyCookbooksResponse {
            cookbooks: cookbooks,
        }))
    }

    async fn release_my_shares(
        &self,
        request: Request<ReleaseMySharesRequest>,
    ) -> Result<Response<ReleaseMySharesResponse>, Status> {
        let req: ReleaseMySharesRequest = request.into_inner();
        let me = users::get_current_username();
        if me.is_none() {
            return Err(Status::new(Code::Internal, "unable to determine user name"));
        }
        let me = me.unwrap();
        let me = me.to_string_lossy();
        for i in 0..req.filenames.len() {
            let path = format!("{}/{}/{}", req.share_dir, me, req.filenames[i]);
            if path_exists(&path) {
                let res = std::fs::remove_file(&path);
                if res.is_err() {
                    return Err(Status::new(Code::Internal, "unable to remove file"));
                }
            } else {
                return Err(Status::new(
                    Code::Internal,
                    format!("file to be removed does not exist: {}", path),
                ));
            }
        }
        Ok(Response::new(ReleaseMySharesResponse { ok: true }))
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

    // Some some info.

    let current_dir = std::env::current_dir()?;
    let current_dir = current_dir.display();
    let current_executable = std::env::current_exe()?;
    let current_executable = current_executable.display();

    // Get version.
    //
    // This is what we had before.  It is confusing and unsound:
    // 1. current_version_string() can fail when it invokes the git command
    // 2. it seems like we set version, and then sometimes set version to the same thing.
    //
    // Therefore, this code is commented out.  But perhaps something in it needs to be salvaged.
    //
    // let mut version = current_version_string();
    // if format!("{}", current_executable) == format!("{}/target/debug/enclone", current_dir) {
    //     version = current_version_string();
    // }
    //
    // Subsequently we changed this to
    //
    // let version = env!("CARGO_PKG_VERSION");
    //
    // but that wasn't right either, because enclone_client.rs would read that back and do
    // the wrong thing.
    //
    // So finally we changed it to what is below.  Note that this suffers from the defect that
    // it is not properly updated.

    let version = version_string();

    // Announce.

    let mut emsg = format!("I am process {}.\n", std::process::id());
    emsg += &mut format!("enclone version = {}\n", env!("CARGO_PKG_VERSION"));
    emsg += &mut format!("version string = {}\n", version);
    emsg += &mut format!(
        "Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    );
    emsg += &mut format!(
        "Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    );
    emsg += &mut format!(
        "Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    );
    emsg += &mut format!(
        "Welcome!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    );
    eprint!("{}", emsg);

    println!("current dir = {}", current_dir);
    println!("current executable = {}", current_executable);
    Server::builder()
        .add_service(AnalyzerServer::new(analyzer))
        .serve_with_incoming(TcpListenerStream::new(listener))
        .await?;

    Ok(())
}
