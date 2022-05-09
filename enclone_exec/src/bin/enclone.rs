// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Note: if enclone is run from the command line, and fails, it will still return exit status
// zero.  As far as we know, in all other cases where it is not run from the command line, it
// returns exit status zero.

use chrono::prelude::*;
use enclone_core::combine_group_pics::combine_group_pics;
use enclone_main::main_enclone::main_enclone;
use enclone_main::USING_PAGER;
#[cfg(feature = "enclone_visual")]
use enclone_visual::enclone_client::enclone_client;
#[cfg(feature = "enclone_visual")]
use enclone_visual::enclone_server::enclone_server;
#[cfg(feature = "enclone_visual")]
use enclone_visual::gui_structures::Summary;
#[cfg(feature = "enclone_visual")]
use enclone_visual::history::write_enclone_visual_history;
#[cfg(feature = "enclone_visual")]
use enclone_visual::history::EncloneVisualHistory;
use flate2::write::GzEncoder;
use flate2::Compression;
use io_utils::*;
use itertools::Itertools;
#[cfg(not(target_os = "windows"))]
use nix::sys::signal::{kill, SIGINT};
#[cfg(not(target_os = "windows"))]
use nix::unistd::getppid;
#[cfg(not(target_os = "windows"))]
use nix::unistd::Pid;
use pretty_trace::*;
use std::env;
use std::io::Write;
use std::process::Command;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::Duration;
use string_utils::*;
use vector_utils::unique_sort;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    PrettyTrace::new().on();
    let mut args: Vec<String> = env::args().collect();
    let mut no_kill = false;
    let mut update = false;
    for i in 1..args.len() {
        if args[i] == "NO_KILL" {
            no_kill = true;
        } else if args[i] == "UPDATE" {
            update = true;
        }
    }

    // Update mode.

    if update {
        if args.len() != 2 {
            eprintln!(
                "\nYou've specified UPDATE, but in update mode, we expect only a single \
                argument.\nIf you'd like to update, please type just enclone UPDATE.\n"
            );
            std::process::exit(1);
        }
        let mut home = String::new();
        for (key, value) in env::vars() {
            if key == "HOME" {
                home = value.clone();
            }
        }
        if home.len() == 0 {
            eprintln!("Weird, unable to determine your home directory.\n");
            std::process::exit(1);
        }
        let datasets = format!("{}/enclone/datasets", home);
        if !path_exists(&datasets) || !std::path::Path::new(&datasets).is_dir() {
            eprintln!(
                "\nSomething odd has happened.  There should be a directory ~/enclone and \
                inside that, a directory datasets.\n"
            );
            std::process::exit(1);
        }
        let list = dir_list(&datasets);
        let size;
        if list.len() <= 3 {
            size = "small";
        } else if list.len() <= 40 {
            size = "medium";
        } else {
            size = "large";
        }
        println!("updating enclone using size = {}", size);
        let o = Command::new("bash")
            .arg("-c")
            .arg(&format!(
                "curl -sSf -L bit.ly/enclone_install | bash -s {}",
                size
            ))
            .output()
            .expect("failed to execute curl");
        print!("{}{}", strme(&o.stdout), strme(&o.stderr));
        std::process::exit(0);
    }

    // Client run of enclone.
    #[cfg(feature = "enclone_visual")]
    {
        let t = std::time::Instant::now();
        for i in 0..args.len() {
            if args[i] == "VIS" || args[i].starts_with("VIS=") {
                enclone_client(&t).await?;
            }
        }
    }

    // Standard run of enclone.

    if args.len() < 2 || args[1] != "SERVER" {
        let res = main_enclone(&mut args);

        // Test for error.

        if res.is_err() {
            // TURNED OFF BECAUSE WE GOT EXIT STATUS ZERO SOMETIMES WHEN WE USED THROUGH COMMAND.
            //
            // If there was an error and we had used the pager, then std::process::exit(1) will
            // result in exit status 0 if enclone was invoked from a terminal window, and
            // probably not otherwise.  To get nonzero exit status, we instead kill the parent
            // process, which is less.  That's a little surprising, but that's how it works:
            // it is the parent that is less.
            //
            // The little bit of sleep seems to prevent printing of an extra newline, but this
            // is flaky.
            //
            // The kill makes the screen flash.  This is pretty horrible.

            eprintln!("{}", res.as_ref().err().unwrap());
            if !no_kill && USING_PAGER.load(SeqCst) && 0 == 1 {
                thread::sleep(Duration::from_millis(10));
                #[cfg(not(target_os = "windows"))]
                {
                    let ppid = getppid();
                    kill(Pid::from_raw(i32::from(ppid)), SIGINT).unwrap();
                }
                thread::sleep(Duration::from_millis(10));
            } else {
                std::process::exit(1);
            }
        }

        // Process VIS_DUMP.

        #[cfg(feature = "enclone_visual")]
        {
            let res = res.unwrap();
            let ctl = &res.inter.setup.ctl;
            if ctl.gen_opt.vis_dump {
                let mut evh = EncloneVisualHistory::default();
                let mut svg = String::new();
                if !res.outs.svgs.is_empty() {
                    svg = res.outs.svgs.last().unwrap().clone();
                }
                evh.svg_hist_uniq.push(svg.clone());
                evh.svg_history.push(0);

                let dataset_names = &res.outs.dataset_names;
                let metrics0 = &res.outs.metrics;
                let mut metrics = Vec::<Vec<String>>::new();
                let n = metrics0.len();
                for i in 0..n {
                    let m = &metrics0[i];
                    let mut lines = Vec::<String>::new();
                    for line in m.lines() {
                        lines.push(line.to_string());
                    }
                    metrics.push(lines);
                }
                let mut all_metric_names = Vec::<String>::new();
                for i in 0..metrics.len() {
                    for j in 0..metrics[i].len() {
                        let s = parse_csv(&metrics[i][j]);
                        let m = format!("{},{}", s[0], s[1]);
                        all_metric_names.push(m);
                    }
                }
                unique_sort(&mut all_metric_names);
                let nm = all_metric_names.len();
                let s = Summary {
                    summary: res.outs.summary.clone(),
                    dataset_names: dataset_names.to_vec(),
                    metrics: metrics,
                    metric_selected: vec![false; nm],
                    metrics_condensed: false,
                };
                let summary = s.pack();

                evh.summary_hist_uniq.push(summary);
                evh.summary_history.push(0);
                let mut args2 = Vec::<String>::new();
                for i in 0..args.len() {
                    if args[i] == "VIS_DUMP" {
                        continue;
                    }
                    if args[i].starts_with("SESSION_NAME=") {
                        continue;
                    }
                    if args[i].starts_with("SESSION_NARRATIVE=") {
                        continue;
                    }
                    if args[i].starts_with("STATE_NARRATIVE=") {
                        continue;
                    }
                    args2.push(args[i].clone());
                }
                let command = format!("{}", args2.iter().format(" "));
                evh.input1_hist_uniq.push(command.clone());
                evh.input1_history.push(0);
                evh.input2_hist_uniq.push(String::new());
                evh.input2_history.push(0);
                evh.inputn_hist_uniq.push(vec![String::new(); 8]);
                evh.inputn_history.push(0);
                evh.translated_input_hist_uniq.push(command.clone());
                evh.translated_input_history.push(0);
                evh.history_index = 1;
                evh.is_blank.push(svg.is_empty());
                let mut table = res.outs.pics.clone();
                let widths = res.outs.last_widths.clone();
                const MAX_TABLE: usize = 50;
                if table.len() > MAX_TABLE {
                    table.truncate(MAX_TABLE);
                }
                let mut reply_text = combine_group_pics(
                    &table,
                    &widths,
                    res.outs.parseable_stdouth,
                    res.outs.noprint,
                    res.outs.noprintx,
                    res.outs.html,
                    res.outs.ngroup,
                    res.outs.pretty,
                );
                if reply_text.contains("enclone failed") {
                    reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
                }
                if reply_text.len() == 0 {
                    if command.contains(" NOPRINT") {
                        reply_text = "You used the NOPRINT option, so there are no \
                            clonotypes to see."
                            .to_string();
                    } else {
                        reply_text = "There are no clonotypes.  Please have a look at the summary."
                            .to_string();
                    }
                }
                let mut reply_text_clean = String::new();
                let mut chars = Vec::<char>::new();
                for c in reply_text.chars() {
                    chars.push(c);
                }
                let mut escaped = false;
                for l in 0..chars.len() {
                    if chars[l] == '' {
                        escaped = true;
                    }
                    if escaped {
                        if chars[l] == 'm' {
                            escaped = false;
                        }
                        continue;
                    }
                    reply_text_clean.push(chars[l]);
                }
                reply_text = reply_text_clean;
                reply_text += "\n \n \n"; // papering over truncation bug in display
                evh.displayed_tables_hist_uniq.push(reply_text);
                evh.displayed_tables_history.push(0);
                evh.last_widths_hist_uniq.push(widths);
                evh.last_widths_history.push(0);
                let full_table = res.outs.pics.clone();
                let serialized = serde_json::to_string(&full_table)
                    .unwrap()
                    .as_bytes()
                    .to_vec();
                let mut e = GzEncoder::new(Vec::new(), Compression::default());
                let _ = e.write_all(&serialized);
                let gzipped = e.finish().unwrap();
                evh.table_comp_hist_uniq.push(gzipped);
                evh.table_comp_history.push(0);
                evh.name_value = ctl.gen_opt.session_name.clone();
                evh.orig_name_value = ctl.gen_opt.session_name.clone();
                evh.narrative_hist_uniq
                    .push(ctl.gen_opt.state_narrative.clone());
                evh.narrative_history.push(0);
                evh.descrip_hist_uniq.push(String::new());
                evh.descrip_history.push(0);
                evh.narrative = ctl.gen_opt.session_narrative.clone();
                let mut now = format!("{:?}", Local::now());
                now = now.replace("T", "___");
                now = now.before(".").to_string();
                now = now.replace(":", "-");
                let home = home::home_dir();
                if home.is_none() {
                    eprintln!(
                        "Unable to determine home directory.  This is unexpected and \
                        pathological.\nPlease report this problem!\n"
                    );
                    std::process::exit(1);
                }
                let enclone = format!("{}/enclone", home.unwrap().display());
                let history = format!("{}/visual/history", enclone);
                if !path_exists(&history) {
                    let res = std::fs::create_dir_all(&history);
                    if res.is_err() {
                        eprintln!(
                            "Unable to create the directory ~/enclone/visual/history.  \
                            This is odd and unexpected.\nPlease report this problem!\n"
                        );
                        std::process::exit(1);
                    }
                }
                let filename = format!("{}/{}", history, now);
                let res = write_enclone_visual_history(&evh, &filename);
                if res.is_err() {
                    eprintln!("\nWas unable to write history to the file {}.\n", filename);
                    std::process::exit(1);
                }
            }
        }

        // Done.

        std::process::exit(0);
    }

    // Server run of enclone.
    #[cfg(feature = "enclone_visual")]
    {
        enclone_server().await?;
    }
    Ok(())
}
