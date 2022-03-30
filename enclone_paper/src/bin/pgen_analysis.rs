// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Analyze pgen data.
//
// Usage:
// pgen_analysis per_cell_stuff per_cell_stuff_pgen per_cell_stuff_pgen_sim
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=datasets_cell,donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,
//         const1,hcomp,jun_ins,d1_name1
//         > per_cell_stuff

use io_utils::*;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;
use tables::*;
use vector_utils::unique_sort;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load simulated heavy chain probabilities.

    let f = open_for_read![&args[3]];
    let mut first = true;
    let mut hlike_sim = Vec::<f64>::new();
    let mut tof = HashMap::<String, usize>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("haa_pgen"));
            first = false;
        } else {
            hlike_sim.push(-fields[tof["haa_pgen"]].force_f64().log10());
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load heavy chain probabilities.

    let f = open_for_read![&args[2]];
    let mut first = true;
    let mut hlike = Vec::<f64>::new();
    let mut tof = HashMap::<String, usize>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("haa_pgen"));
            first = false;
        } else {
            hlike.push(-fields[tof["haa_pgen"]].force_f64().log10());
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load main data.

    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(
        String,
        usize,
        Vec<u8>,
        String,
        String,
        usize,
        usize,
        String,
        usize,
        usize,
        String,
        f64,
    )>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("datasets_cell"));
            assert!(tof.contains_key("donors_cell"));
            assert!(tof.contains_key("v_name1"));
            assert!(tof.contains_key("v_name2"));
            assert!(tof.contains_key("dref"));
            assert!(tof.contains_key("cdr3_aa1"));
            assert!(tof.contains_key("clonotype_ncells"));
            assert!(tof.contains_key("const1"));
            assert!(tof.contains_key("hcomp"));
            assert!(tof.contains_key("jun_ins"));
            assert!(tof.contains_key("d1_name1"));
            first = false;
        } else {
            data.push((
                /* 0 */ fields[tof["v_name1"]].to_string(),
                /* 1 */ fields[tof["cdr3_aa1"]].len(),
                /* 2 */ fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                /* 3 */ fields[tof["donors_cell"]].to_string(),
                /* 4 */ fields[tof["v_name2"]].to_string(),
                /* 5 */ fields[tof["dref"]].force_usize(),
                /* 6 */ fields[tof["clonotype_ncells"]].force_usize(),
                /* 7 */ fields[tof["const1"]].to_string(),
                /* 8 */ fields[tof["hcomp"]].force_usize(),
                /* 9 */ fields[tof["datasets_cell"]].force_usize(),
                /* 10 */ fields[tof["d1_name1"]].to_string(),
                /* 11 */ hlike[data.len()],
            ));
        }
    }
    data.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Show distribution of hlike for naive cells, and for simulated, and for public naive cells.

    let mut rows = Vec::<Vec<String>>::new();
    let row = vec![
        "hlike".to_string(),
        "naive%".to_string(),
        "sim%".to_string(),
        "public naive%".to_string(),
    ];
    rows.push(row);
    let mut nhb = vec![0; 100];
    let mut inf = 0;
    let mut nhb_sim = vec![0; 100];
    let mut inf_sim = 0;
    let mut total = 0;
    for i in 0..data.len() {
        if data[i].5 == 0 {
            total += 1;
            let n = data[i].11.round() as usize;
            if n < 50 {
                nhb[n] += 1;
            } else {
                inf += 1;
            }
            let n = hlike_sim[total].round() as usize;
            if n < 50 {
                nhb_sim[n] += 1;
            } else {
                inf_sim += 1;
            }
        }
    }

    let mut x = Vec::<f64>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        let mut donors = Vec::<String>::new();
        for k in i..j {
            donors.push(data[k].3.clone());
        }
        unique_sort(&mut donors);
        if donors.len() > 1 {
            for k in i..j {
                let dref = data[k].5;
                if dref == 0 {
                    x.push(data[k].11);
                }
            }
        }
        i = j;
    }

    let mut nhb_pub = vec![0; 100];
    let mut inf_pub = 0;
    for i in 0..x.len() {
        let n = x[i].round() as usize;
        if n < 50 {
            nhb_pub[n] += 1;
        } else {
            inf_pub += 1;
        }
    }

    for i in 0..100 {
        if nhb[i] > 0 || nhb_pub[i] > 0 {
            let mut row = Vec::<String>::new();
            row.push(format!("{}", i));
            row.push(format!("{:.3}", 100.0 * nhb[i] as f64 / total as f64));
            row.push(format!("{:.3}", 100.0 * nhb_sim[i] as f64 / total as f64));
            row.push(format!("{:.3}", 100.0 * nhb_pub[i] as f64 / x.len() as f64));
            rows.push(row);
        }
    }
    let mut row = Vec::<String>::new();
    row.push(">= 50".to_string());
    row.push(format!("{:.3}", 100.0 * inf as f64 / total as f64));
    row.push(format!("{:.3}", 100.0 * inf_sim as f64 / total as f64));
    row.push(format!("{:.3}", 100.0 * inf_pub as f64 / x.len() as f64));
    rows.push(row);

    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 0, &b"r|r|r|r".to_vec(), false, false);
    println!("\n{}", log);
}
