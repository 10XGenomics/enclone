// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Process commands via the server in the background.

use crate::proto::{analyzer_client::AnalyzerClient, *};
use crate::summary::*;
use crate::*;
use std::io::Read;
use std::process::Child;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::Duration;
use tonic::transport::Channel;

pub async fn process_requests(
    client: &mut AnalyzerClient<Channel>,
    server_process: &mut Child,
    verbose: bool,
) {
    loop {
        thread::sleep(Duration::from_millis(10));
        if DONE.load(SeqCst) {
            cleanup();
            std::process::exit(0);
        }

        if GET_MY_COOKBOOKS.load(SeqCst) {
            let cookbook_dirs = COOKBOOK_DIRS.lock().unwrap().clone();
            let request = tonic::Request::new(GetMyCookbooksRequest {
                cookbook_dirs: cookbook_dirs,
            });
            let response = client.get_my_cookbooks(request).await;
            if response.is_err() {
                eprintln!("Attempt to get cookbooks failed.");
                let err = format!("{:?}", response);
                eprintln!("err = {}\n", err);
                eprintln!("Please ask for help!\n");
                std::process::exit(1);
            }
            let res = response.unwrap().into_inner();
            for i in 0..res.cookbooks.len() {
                REMOTE_COOKBOOKS
                    .lock()
                    .unwrap()
                    .push(res.cookbooks[i].clone());
            }
            GET_MY_COOKBOOKS.store(false, SeqCst);
        }

        if SENDING_SHARE.load(SeqCst) {
            let share_dir = REMOTE_SHARE.lock().unwrap()[0].clone();
            let sender = users::get_current_username();
            let sender = sender.unwrap().to_string_lossy().to_string();
            let content = SHARE_CONTENT.lock().unwrap()[0].clone();
            let nrecip = SHARE_RECIPIENTS.lock().unwrap().len();
            let mut recipients = Vec::<String>::new();
            for i in 0..nrecip {
                recipients.push(SHARE_RECIPIENTS.lock().unwrap()[i].clone());
            }
            let request = tonic::Request::new(SendShareRequest {
                share_dir: share_dir,
                content: content,
                sender: sender,
                recipients: recipients,
            });
            let response = client.share_session(request).await;
            if response.is_err() {
                eprintln!("Attempt to share session failed.");
                let err = format!("{:?}", response);
                eprintln!("err = {}\n", err);
                eprintln!("Please ask for help!\n");
                std::process::exit(1);
            }
            SENDING_SHARE.store(false, SeqCst);
        }

        if GET_MY_SHARES.load(SeqCst) {
            let share_dir = REMOTE_SHARE.lock().unwrap()[0].clone();
            let request = tonic::Request::new(GetMySharesRequest {
                share_dir: share_dir,
                me_only: META_TESTING.load(SeqCst),
            });
            let response = client.get_my_shares(request).await;
            if response.is_err() {
                eprintln!("Attempt to retrieve shared sessions failed.");
                let err = format!("{:?}", response);
                eprintln!("err = {}\n", err);
                eprintln!("Please ask for help!\n");
                std::process::exit(1);
            }
            let res = response.unwrap().into_inner();
            let n = res.content.len();
            for i in 0..n {
                RECEIVED_SHARES_CONTENT
                    .lock()
                    .unwrap()
                    .push(res.content[i].clone());
                RECEIVED_SHARES_MESSAGES
                    .lock()
                    .unwrap()
                    .push(res.messages[i].clone());
                RECEIVED_SHARES_FILENAMES
                    .lock()
                    .unwrap()
                    .push(res.filenames[i].clone());
            }
            GET_MY_SHARES.store(false, SeqCst);
        }

        if RELEASE_MY_SHARES.load(SeqCst) {
            let share_dir = REMOTE_SHARE.lock().unwrap()[0].clone();
            let n = RECEIVED_SHARES_FILENAMES.lock().unwrap().len();
            let mut filenames = Vec::<String>::new();
            for i in 0..n {
                filenames.push(RECEIVED_SHARES_FILENAMES.lock().unwrap()[i].clone());
            }
            let request = tonic::Request::new(ReleaseMySharesRequest {
                share_dir: share_dir,
                filenames: filenames,
            });
            let response = client.release_my_shares(request).await;
            if response.is_err() {
                eprintln!("Attempt to release shared sessions failed.");
                let err = format!("{:?}", response);
                eprintln!("err = {}\n", err);
                eprintln!("Please ask for help!\n");
                std::process::exit(1);
            }
            RECEIVED_SHARES_CONTENT.lock().unwrap().clear();
            RECEIVED_SHARES_MESSAGES.lock().unwrap().clear();
            RECEIVED_SHARES_FILENAMES.lock().unwrap().clear();
            RELEASE_MY_SHARES.store(false, SeqCst);
        }

        if TESTING_USER_NAME.load(SeqCst) {
            let user_name = USER_NAME.lock().unwrap()[0].clone();
            let request = tonic::Request::new(UserNameRequest {
                user_name: user_name,
            });
            let response = client.test_user_name(request).await;
            if response.is_err() {
                eprintln!("\nWeird, validity test for user name failed.");
                let err = format!("{:?}", response);
                eprintln!("err = {}\n", err);
                std::process::exit(1);
            }
            let valid = response.unwrap().into_inner().value;
            USER_NAME_VALID.store(valid, SeqCst);
            TESTING_USER_NAME.store(false, SeqCst);
        }

        if PROCESSING_REQUEST.load(SeqCst) {
            let input = USER_REQUEST.lock().unwrap()[0].clone();
            let mut line = input.to_string();
            if verbose {
                xprintln!("processing command {}", line);
            }
            let output;
            let mut svg_output = String::new();
            let mut summary = String::new();
            let mut metrics = Vec::<String>::new();
            let mut dataset_names = Vec::<String>::new();
            let mut table_comp = Vec::<u8>::new();
            let mut last_widths = Vec::<u32>::new();
            if line != "enclone" && !line.starts_with("enclone ") {
                if line.starts_with("#") {
                    output = "unknown tag".to_string();
                } else {
                    output =
                        "an actual enclone command needs to start with \"enclone\"".to_string();
                    if FAIL_ON_ERROR.load(SeqCst) {
                        eprintln!("\nFAIL_ON_ERROR: command must start with enclone\n");
                        std::process::exit(1);
                    }
                }
            } else {
                if CONFIG_FILE.lock().unwrap().len() > 0 {
                    line += &mut format!(" CONFIG={}", CONFIG_FILE.lock().unwrap()[0]);
                }
                let mut server_logfile = None;
                if SERVER_LOGFILE.lock().unwrap().len() > 0 {
                    server_logfile = Some(SERVER_LOGFILE.lock().unwrap()[0].clone());
                }
                let request = tonic::Request::new(EncloneRequest {
                    args: line,
                    server_logfile: server_logfile,
                });
                let response = client.enclone(request).await;
                if response.is_err() {
                    let left = r###"message: "\n"###;
                    let right = r###"\n""###;
                    let mut err = format!("{:?}", response);
                    if err.contains(&left) && err.after(&left).contains(&right) {
                        err = err.between(&left, &right).to_string();
                    }
                    err = err.replace("\\n", "\n");
                    let crash = err.contains("transport error: connection error: broken pipe");
                    if crash {
                        output = format!(
                            "\nIt would appear the the enclone server \
                            crashed.\nPlease look in the terminal window for a traceback \
                            and report it.\n"
                        );
                    } else {
                        let msg = format!("\nThe enclone server is unhappy.  It says:\n\n{}", err);
                        output = msg.clone();
                        xprintln!("{}", msg);
                    }
                    let mut ebuffer = [0; 10000];
                    let server_stderr = server_process.stderr.as_mut().unwrap();
                    server_stderr.read(&mut ebuffer).unwrap();
                    let emsg = strme(&ebuffer);
                    xprint!("server error =\n{}", emsg);
                } else {
                    let response = response.unwrap();
                    let r = response.into_inner();
                    svg_output = r.plot.clone();
                    metrics = r.metrics.clone();
                    dataset_names = r.dataset_names.clone();
                    output = format!("{}", r.table);
                    table_comp = r.table_comp.clone();
                    for x in r.last_widths.iter() {
                        last_widths.push(*x as u32);
                    }
                    summary = r.summary.clone();
                }
            }
            SERVER_REPLY_SVG.lock().unwrap().clear();
            SERVER_REPLY_SVG.lock().unwrap().push(svg_output);
            SERVER_REPLY_SUMMARY.lock().unwrap().clear();
            SERVER_REPLY_SUMMARY.lock().unwrap().push(summary);
            SERVER_REPLY_METRICS.lock().unwrap().clear();
            let n = metrics.len();
            for i in 0..n {
                SERVER_REPLY_METRICS
                    .lock()
                    .unwrap()
                    .push(metrics[i].clone());
            }
            SERVER_REPLY_DATASET_NAMES.lock().unwrap().clear();
            let n = dataset_names.len();
            for i in 0..n {
                SERVER_REPLY_DATASET_NAMES
                    .lock()
                    .unwrap()
                    .push(dataset_names[i].clone());
            }
            let summary = form_summary_from_server_response().pack();
            SERVER_REPLY_SUMMARY_PLUS.lock().unwrap().clear();
            SERVER_REPLY_SUMMARY_PLUS.lock().unwrap().push(summary);
            SERVER_REPLY_TABLE_COMP.lock().unwrap().clear();
            SERVER_REPLY_TABLE_COMP.lock().unwrap().push(table_comp);
            SERVER_REPLY_LAST_WIDTHS.lock().unwrap().clear();
            SERVER_REPLY_LAST_WIDTHS.lock().unwrap().push(last_widths);
            SERVER_REPLY_TEXT.lock().unwrap().clear();
            SERVER_REPLY_TEXT.lock().unwrap().push(output.clone());
            PROCESSING_REQUEST.store(false, SeqCst);
        }
    }
}
