// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Run pgen in parallel on some simulated sequences.
//
// Usage:
// pgen_sim_compute dir source
// where dir is a directory containing a CSV file source, with fields
// junction_dna,junction_aa,heavy_v_gene,heavy_j_gene
// but no header line.
// 
// This assumes:
// 1. You have olga-compute_pgen set up in your environment.
// 2. You have a working qsub.

use io_utils::*;
use pretty_trace::*;
use std::process::Command;
use std::env;
use std::io::{BufRead, Write};
use string_utils::strme;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (dir, source) = (&args[1], &args[2]);
    let scripts = format!("{}/scripts", dir);
    let ins = format!("{}/ins", dir);
    let (outs, errs) = (format!("{}/outs", dir), format!("{}/errs", dir));
    let results = format!("{}/results", dir);
    for d in [&scripts, &ins, &outs, &errs, &results].iter() {
        if !path_exists(&d) {
            std::fs::create_dir(&d).unwrap();
        }
    }
    let inputs = format!("{}/{}", dir, source);
    let mut lines = Vec::<String>::new();
    let f = open_for_read![&inputs];
    for line in f.lines() {
        let s = line.unwrap();
        lines.push(s);
    }
    let n = lines.len();
    const BATCH: usize = 1000;
    let mut count = 1;
    for start in (0..n).step_by(BATCH) {
        let stop = std::cmp::min(start + BATCH, n);
        let source = format!("{ins}/{count}.csv");
        let mut f = open_for_write_new![&source];
        for j in start..stop {
            fwriteln!(f, "{}", lines[j]);
        }
        let script = format!("{scripts}/{count}.sh");
        let mut g = open_for_write_new![&script];
        fwriteln!(g, 
            "#!/usr/bin/env bash\n\
            #$ -N \"pgen run{count}\"\n\
            #$ -pe threads 1\n\
            #$ -l mem_free=1G\n\
            #$ -o {outs}/{count}.out\n\
            #$ -o {errs}/{count}.err\n\
            #$ -l h_rt=48:00:00\n\
            #$ -S \"/usr/bin/env bash\"\n\
            #$ -V\n\
            #$ -cwd\n\
            olga-compute_pgen --delimiter ',' --comment_delimiter=# --humanIGH \
                -i {source} --seq_in 0 -o {results}/pgen.{count}"
        );
        let new = Command::new("qsub")
            .arg(&script)
            .output()
            .expect(&format!("failed to execute qsub"));
        if new.status.code() != Some(0) {
            eprint!(
                "\nsomething went wrong running qsub, stderr =\n{}",
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        count += 1;
    }
}
