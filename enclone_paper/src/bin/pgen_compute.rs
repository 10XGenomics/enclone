// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Run pgen in parallel on some sequences.
//
// Usage:
// 1. pgen_compute dir source
// where dir is a directory containing a CSV file source, with fields
// junction_dna,junction_aa,heavy_v_gene,heavy_j_gene
// but no header line.
// 2. wait until qsubbed jobs finish (otherwise next step will fail gracefully)
// 3. pgen_compute dir source MERGE
// 4. output is dir/source.out.
//
// This assumes:
// 1. You have olga-compute_pgen set up in your environment.
// 2. You have a working qsub.

use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, Write};
use std::process::{Command, Stdio};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (dir, source) = (&args[1], &args[2]);
    let mut merge = false;
    if args.len() >= 4 {
        if args[3] == "MERGE" {
            merge = true;
        } else {
            eprintln!("Illegal argument.");
            std::process::exit(1);
        }
    }
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
    if merge {
        let mut lines = Vec::<String>::new();
        let mut count = 1;
        for start in (0..n).step_by(BATCH) {
            let stop = std::cmp::min(start + BATCH, n);
            let len = stop - start;
            let res = format!("{results}/pgen.{count}");
            if !path_exists(&res) {
                println!("{count} not done");
                std::process::exit(0);
            }
            let f = open_for_read![&res];
            let mut k = 0;
            for line in f.lines() {
                let s = line.unwrap();
                lines.push(s);
                k += 1;
            }
            if k != len {
                println!("{count} not done, has only {k} lines of {len}");
                std::process::exit(0);
            }
            count += 1;
        }
        let outfile = format!("{inputs}.out");
        {
            let mut f = open_for_write_new![&outfile];
            for line in lines.iter() {
                fwriteln!(f, "{}", line);
            }
        }
        std::process::exit(0);
    }
    let mut count = 1;
    for start in (0..n).step_by(BATCH) {
        let stop = std::cmp::min(start + BATCH, n);
        let source = format!("{ins}/{count}.csv");
        let mut f = open_for_write_new![&source];
        for j in start..stop {
            fwriteln!(f, "{}", lines[j]);
        }
        let script = format!("{scripts}/{count}.sh");
        let results_file = format!("{results}/pgen.{count}");
        if path_exists(&results_file) {
            std::fs::remove_file(&results_file).unwrap();
        }
        {
            let mut g = open_for_write_new![&script];
            fwriteln!(
                g,
                "#!/usr/bin/env bash\n\
            #$ -N \"pgen run{count}\"\n\
            #$ -pe threads 1\n\
            #$ -l mem_free=1G\n\
            #$ -o {outs}/{count}.out\n\
            #$ -e {errs}/{count}.err\n\
            #$ -l h_rt=48:00:00\n\
            #$ -S \"/usr/bin/env bash\"\n\
            #$ -V\n\
            #$ -cwd\n\
            olga-compute_pgen --delimiter ',' --comment_delimiter=# --humanIGH \
                -i {source} --seq_in 0 -o {results_file}"
            );
        }
        let script = File::open(&script).unwrap();
        let child = Command::new("qsub")
            .stdin(script)
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let output = child.wait_with_output().unwrap();
        assert!(output.status.success());
        count += 1;
    }
}
