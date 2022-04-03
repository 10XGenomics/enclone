// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Splice junctions from olga.
// An optional second argument is the number of sequences to process.

use bio_edit::alignment::AlignmentOperation::{Del, Ins, Match, Subst};
use debruijn::dna_string::DnaString;
use enclone_core::align_to_vdj_ref::align_to_vdj_ref;
use enclone_core::defs::Junction;
use enclone_core::opt_d::{jflank, opt_d};
use enclone_tail::align_n::print_vis_align;
use io_utils::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::cmp::min;
use std::env;
use std::io::{BufRead, Write};
use string_utils::strme;
use string_utils::TextUtils;
use vdj_ann::annotate::{annotate_seq, get_cdr3_using_ann};
use vdj_ann::refx::make_vdj_ref_data_core;
use vdj_ann::refx::RefData;
use vdj_ann_ref::human_ref;
use vector_utils::make_freq;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut n = 0;
    if args.len() >= 3 {
        n = args[2].force_usize();
    }

    // Get the human reference used by enclone.

    let mut refnames = Vec::<String>::new();
    let mut refs = Vec::<Vec<u8>>::new();
    let mut utr = Vec::<bool>::new();
    let href = human_ref();
    {
        for (i, line) in href.lines().enumerate() {
            if i % 2 == 0 {
                let n = line.between("|", " ").to_string();
                utr.push(line.contains("5'UTR"));
                refnames.push(n);
            } else {
                refs.push(line.as_bytes().to_vec());
            }
        }
    }

    // Add extra genes.  For IGHV3-NL1, we took the sequence in olga.seq, and found the leader
    // in GenBank AC245166.2:
    // ATGGAGAAATAGAGAGACTGAGTGTGAGTGAACAT
    // GAGTGAGAAAAACTGGATTTGTGTGGCATTTTCTGATAACGGTGTCCTTCTGTTTGCAGGTGTCCAGTGT.
    // However, this leader has stop codons.  Therefore we put in a FAKE leader, that for
    // IGHV3-11.

    refnames.push("IGHV3-43D".to_string());
    utr.push(false);
    refs.push(
        b"ATGGAGTTTGGACTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCT\
        GGGGGAGTCGTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTG\
        GGTCCGTCAAGCTCCGGGGAAGGGTCTGGAGTGGGTCTCTCTTATTAGTTGGGATGGTGGTAGCACCTACTATGCAGACTCTGTGA\
        AGGGTCGATTCACCATCTCCAGAGACAACAGCAAAAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACCGCCTTG\
        TATTACTGTGCAAAAGATA"
            .to_vec(),
    );
    refnames.push("IGHV3-30-3".to_string());
    utr.push(false);
    refs.push(
        b"ATGGAGTTTGGGCTGAGCTGGGTTTTCCTCGTTGCTCTTTTAAGAGGTGTCCAGTGTCAGGTGCAGCTGGTGGAGTCT\
        GGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTG\
        GGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGCAATAAATACTACGCAGACTCCGTGA\
        AGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTG\
        TATTACTGTGCGAGAGA"
            .to_vec(),
    );
    refnames.push("IGHV3-NL1".to_string());
    refs.push(
        b"ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTATAAAAGGTG\
        TCCAGTGTCAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCCTGAGACT\
        CTCCTGTGCAGCGTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTCT\
        CAGTTATTTATAGCGGTGGTAGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAAC\
        ACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAAAGA"
            .to_vec(),
    );
    utr.push(false);
    refnames.push("IGHV4-30-2".to_string());
    refs.push(
        b"ATGAAACACCTGTGGTTCTTCCTCCTGCTGGTGGCAGCTCCCAGATGGGTCCTGTCCCAGCTGCAGCTGCAGGAGTCC\
        GGCTCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTCCTG\
        GAGCTGGATCCGGCAGCCACCAGGGAAGGGCCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCC\
        TCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCAGACACGGCT\
        GTGTATTACTGTGCGAGACA"
            .to_vec(),
    );
    utr.push(false);
    let nref = refs.len();
    let mut refdata = RefData::new();
    let ext_ref = String::new();
    make_vdj_ref_data_core(&mut refdata, &href, &ext_ref, false, true, None);

    // Load the input file.

    let mut jun = Vec::<Vec<u8>>::new();
    let mut hv = Vec::<String>::new();
    let mut hj = Vec::<String>::new();
    let mut cdr3 = Vec::<String>::new();
    let mut cdr3_dna = Vec::<String>::new();
    let f = open_for_read![&args[1]];
    for line in f.lines() {
        let s = line.unwrap();
        let fields = s.split(',').collect::<Vec<&str>>();
        jun.push(fields[0].as_bytes().to_vec());
        hv.push(fields[2].to_string());
        hj.push(fields[3].to_string());
        cdr3.push(fields[1].to_string());
        cdr3_dna.push(fields[0].to_string());
    }
    if n == 0 {
        n = jun.len();
    }

    // Process entries.

    #[derive(Default)]
    struct Data {
        out: Vec<u8>,
        fail: bool,
        drefname: String,
        d: String,
        jun_ins: usize,
        subs: usize,
        rate: f64,
        cdrh3_len: usize,
    }
    let mut results = Vec::<(usize, Data)>::new();
    for i in 0..n {
        results.push((i, Data::default()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut log = Vec::<u8>::new();
        fwriteln!(
            log,
            "\n-------------------------------------------------------------------------\
            --------------------------"
        );
        fwriteln!(log, "\nsequence {}", i + 1);
        let mut vseq = Vec::<u8>::new();
        let mut jseq = Vec::<u8>::new();
        for j in 0..nref {
            if refnames[j] == hv[i] && !utr[j] {
                vseq = refs[j].clone();
                break;
            }
        }
        if vseq.is_empty() {
            eprintln!(
                "\nEntry {}, unable to find V gene named {}.\n",
                i + 1,
                hv[i]
            );
            std::process::exit(1);
        }
        for j in 0..nref {
            if refnames[j] == hj[i] {
                jseq = refs[j].clone();
                break;
            }
        }
        if jseq.is_empty() {
            eprintln!(
                "\nEntry {}, unable to find J gene named {}.\n",
                i + 1,
                hj[i]
            );
            std::process::exit(1);
        }
        let junv = &jun[i][0..2];
        fwriteln!(log, "\nV gene = {}", hv[i]);
        let jrefname = &hj[i];
        fwriteln!(log, "J gene = {}", jrefname);
        // fwriteln!(log, "\nlooking for {} in {}", strme(&junv), strme(&vseq));
        let mut vstart = None;
        for k in (0..=vseq.len() - junv.len()).rev() {
            if vseq[k..].starts_with(&junv) {
                if k % 3 != 0 {
                    continue;
                }
                vstart = Some(k);
                // fwriteln!(log, "vstart for {} found at pos -{}", i + 1, vseq.len() - k);
                break;
            }
        }
        if vstart.is_none() {
            eprintln!("\nfailed to find vstart for entry {}\n", i + 1);
            std::process::exit(1);
        }
        let mut jtail = jseq[jseq.len() - 31..].to_vec(); // concatenate to junction
        let mut seq = vseq[0..vstart.unwrap()].to_vec();
        seq.append(&mut jun[i].clone());
        seq.append(&mut jtail);
        fwriteln!(log, "\ncdr3[{}] = {}\n", i + 1, cdr3[i]);
        fwriteln!(log, "seq[{}] = {}\n", i + 1, strme(&seq));
        let x = DnaString::from_dna_string(&strme(&seq));
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&x, &refdata, &mut ann, true, false, true);
        let mut annv = Vec::<(i32, i32, i32, i32, i32)>::new();
        for i in 0..ann.len() {
            let t = ann[i].2 as usize;
            if refdata.is_v(t) {
                annv.push(ann[i]);
            }
        }
        let mut cdr3x = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        get_cdr3_using_ann(&x, &refdata, &ann, &mut cdr3x);
        if cdr3x.len() != 1 || strme(&cdr3x[0].1) != cdr3[i] {
            fwriteln!(log, "failed to find unique or matching CDR3\n");
            fwriteln!(log, "found {} CDR3s\n", cdr3x.len());
            for i in 0..cdr3x.len() {
                fwriteln!(log, "CDR3 = {}", strme(&cdr3x[i].1));
            }
            res.1 = Data {
                out: log,
                fail: true,
                drefname: String::new(),
                d: String::new(),
                jun_ins: 0,
                subs: 0,
                rate: 0.0,
                cdrh3_len: 0,
            };
        } else {
            fwriteln!(log, "CDR3 = {}", strme(&cdr3x[0].1));

            // Analyze the junction, following hcomp.rs.

            let mut vref = vseq.clone();
            let cdr3_start = cdr3x[0].0;
            let vstart = cdr3_start - 2;
            vref = vref[vstart..vref.len()].to_vec();
            let mut concat = vref.clone();
            let mut drefx = Vec::<u8>::new();
            let mut d2ref = Vec::<u8>::new();
            let mut drefname = String::new();
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            let jscore_match = 20;
            let jscore_mismatch = -20;
            let jscore_gap_open = -120;
            let jscore_gap_extend = -20;
            let jscore_bits_multiplier = 2.2;
            let mut v_ref_id = None;
            let mut j_ref_id = None;
            for i in 0..ann.len() {
                let t = ann[i].2 as usize;
                if refdata.is_v(t) {
                    v_ref_id = Some(t);
                } else if refdata.is_j(t) {
                    j_ref_id = Some(t);
                }
            }
            let v_ref_id = v_ref_id.unwrap();
            if j_ref_id.is_none() {
                res.1 = Data {
                    out: log,
                    fail: true,
                    drefname: String::new(),
                    d: String::new(),
                    jun_ins: 0,
                    subs: 0,
                    rate: 0.0,
                    cdrh3_len: 0,
                };
            } else {
                let j_ref_id = j_ref_id.unwrap();
                fwriteln!(log, "\nO:{},{},{},{}\n", cdr3_dna[i], cdr3[i], hv[i], hj[i],);
                fwriteln!(
                    log,
                    "Z:{},{},{},{}\n",
                    cdr3_dna[i],
                    cdr3[i],
                    refdata.name[v_ref_id],
                    refdata.name[j_ref_id],
                );
                let j_ref_id = j_ref_id;
                opt_d(
                    v_ref_id,
                    j_ref_id,
                    &seq,
                    &ann,
                    &strme(&cdr3x[0].1),
                    &refdata,
                    &Vec::new(),
                    &mut scores,
                    &mut ds,
                    jscore_match,
                    jscore_mismatch,
                    jscore_gap_open,
                    jscore_gap_extend,
                    jscore_bits_multiplier,
                    None,
                );
                let mut opt = Vec::new();
                if !ds.is_empty() {
                    opt = ds[0].clone();
                }
                for j in 0..opt.len() {
                    let d = opt[j];
                    if j == 0 {
                        drefx = refdata.refs[d].to_ascii_vec();
                    } else {
                        d2ref = refdata.refs[d].to_ascii_vec();
                    }
                    if j > 0 {
                        drefname += ":";
                    }
                    drefname += &mut refdata.name[d].clone();
                }
                concat.append(&mut drefx.clone());
                concat.append(&mut d2ref.clone());
                let mut jref = refdata.refs[j_ref_id].to_ascii_vec();
                let jend = jflank(&seq, &jref);
                let mut seq_start = vstart as isize;
                // probably not exactly right
                if annv.len() > 1 {
                    let q1 = annv[0].0 + ann[0].1;
                    let q2 = annv[1].0;
                    seq_start += q2 as isize - q1 as isize;
                }
                let mut seq_end = seq.len() - (jref.len() - jend);
                if seq_start as usize > seq_end {
                    seq_start = vstart as isize;
                }
                if seq_end <= seq_start as usize {
                    seq_end = seq.len();
                }
                seq = seq[seq_start as usize..seq_end].to_vec();
                jref = jref[0..jend].to_vec();
                concat.append(&mut jref.clone());
                let (ops, _score) = align_to_vdj_ref(
                    &seq,
                    &vref,
                    &drefx,
                    &d2ref,
                    &jref,
                    &drefname,
                    true,
                    jscore_match,
                    jscore_mismatch,
                    jscore_gap_open,
                    jscore_gap_extend,
                    jscore_bits_multiplier,
                );
                let mut tigpos = 0;
                let mut hcomp = 0;
                let mut jun_ins = 0;
                let mut indels = Vec::<(usize, isize)>::new();
                let mut ins_start = 0;
                let mut del_len = 0;
                let mut matches = 0;
                let mut mismatches = 0;
                for i in 0..ops.len() {
                    if ops[i] == Subst {
                        if tigpos >= 2 {
                            mismatches += 1;
                        }
                        hcomp += 1;
                        tigpos += 1;
                    } else if ops[i] == Match {
                        if tigpos >= 2 {
                            matches += 1;
                        }
                        tigpos += 1;
                    } else if ops[i] == Ins {
                        hcomp += 1;
                        jun_ins += 1;
                        if i == 0 || ops[i - 1] != Ins {
                            ins_start = tigpos;
                        }
                        tigpos += 1;
                        if i == ops.len() - 1 || ops[i + 1] != Ins {
                            let ins_len = tigpos - ins_start;
                            indels.push((ins_start, ins_len as isize));
                        }
                    } else if ops[i] == Del {
                        if i == 0 || ops[i - 1] != Del {
                            hcomp += 1;
                        }
                        del_len += 1;
                        if i == ops.len() - 1 || ops[i + 1] != Del {
                            indels.push((tigpos, -(del_len as isize)));
                            del_len = 0;
                        }
                    }
                }
                let _jun = Junction {
                    hcomp: hcomp,
                    matches: matches,
                    mismatches: mismatches,
                    jun_ins: jun_ins,
                    d: ds[0].clone(),
                    vstart: vstart,
                    indels: indels,
                };
                fwriteln!(log, "D = {}", drefname);
                fwriteln!(log, "junction insertion = {}", jun_ins);
                fwriteln!(log, "junction mismatches = {}", mismatches);
                let mut d = String::new();
                for j in 0..ds[0].len() {
                    if j > 0 {
                        d += ":";
                    }
                    d += &mut refdata.name[ds[0][j]].clone();
                }
                let rate = mismatches as f64 / (matches + mismatches) as f64;
                let width = 100;
                let mut drefname = drefname.clone();
                if drefname == *"" {
                    drefname = "none".to_string();
                }
                let rank = "1ST";
                let add = format!("  •  D = {} = {}", rank, drefname);
                let frame = seq_start as usize % 3;
                let pretty = true;
                print_vis_align(
                    &seq,
                    &concat,
                    1,
                    i,
                    &vref,
                    &drefx,
                    &d2ref,
                    &jref,
                    &jrefname,
                    true,
                    &mut log,
                    width,
                    &add,
                    frame,
                    true,
                    jscore_match,
                    jscore_mismatch,
                    jscore_gap_open,
                    jscore_gap_extend,
                    jscore_bits_multiplier,
                    pretty,
                );
                res.1 = Data {
                    out: log,
                    fail: false,
                    drefname: drefname,
                    d: d,
                    jun_ins: jun_ins,
                    subs: mismatches,
                    rate: rate,
                    cdrh3_len: cdr3[i].len(),
                };
            }
        }
    });

    // Tally results and report.

    let mut drefnames = Vec::<String>::new();
    let mut ds_all = Vec::<String>::new();
    let mut subs = Vec::<usize>::new();
    let mut jun_ins = Vec::<usize>::new();
    let mut rates = Vec::<f64>::new();
    let mut fails = 0;
    let mut dd = 0;
    let mut cdrh3_lens = Vec::<usize>::new();
    for i in 0..results.len() {
        print!("{}", strme(&results[i].1.out));
        if results[i].1.fail {
            fails += 1;
        } else {
            let drefname = results[i].1.drefname.clone();
            if drefname.contains(":") {
                dd += 1;
            }
            drefnames.push(drefname);
            ds_all.push(results[i].1.d.clone());
            cdrh3_lens.push(results[i].1.cdrh3_len);
            jun_ins.push(results[i].1.jun_ins);
            if results[i].1.jun_ins == 0 {
                subs.push(results[i].1.subs);
                rates.push(results[i].1.rate);
            }
        }
    }
    println!("\nThere were {} fails.\n", fails);
    println!("CDRH3 length distribution");
    let mut bins = vec![0; 100];
    let mut total = 0;
    for k in 0..cdrh3_lens.len() {
        let len = cdrh3_lens[k];
        bins[len / 5] += 1;
        total += 1;
    }
    for i in 0..bins.len() {
        if bins[i] > 0 {
            println!(
                "{}-{} ==> {:.1}%",
                5 * i,
                5 * (i + 1),
                100.0 * bins[i] as f64 / total as f64
            );
        }
    }
    let mut total = 0;
    for i in 0..jun_ins.len() {
        total += jun_ins[i];
    }
    println!(
        "mean junction insertion length = {:.1}\n",
        total as f64 / jun_ins.len() as f64
    );
    if ds_all.len() > 0 {
        ds_all.sort();
        let mut freq = Vec::<(u32, String)>::new();
        make_freq(&ds_all, &mut freq);
        println!(
            "most frequent D genes for naive cells with junction insertion length 0 (of {})",
            ds_all.len()
        );
        for i in 0..min(10, freq.len()) {
            println!(
                "{} [{:.1}%]",
                freq[i].1,
                100.0 * freq[i].0 as f64 / ds_all.len() as f64
            );
        }
    }
    if subs.len() > 0 {
        subs.sort();
        let mut freq = Vec::<(u32, usize)>::new();
        make_freq(&subs, &mut freq);
        println!(
            "\nmost frequent substitution values for naive cells with junction insertion length 0 (of {})",
            subs.len()
        );
        for i in 0..min(10, freq.len()) {
            println!(
                "{} [{:.1}%]",
                freq[i].1,
                100.0 * freq[i].0 as f64 / subs.len() as f64
            );
        }
        let mut total = 0;
        for i in 0..subs.len() {
            total += subs[i];
        }
        println!("mean = {:.1}", total as f64 / subs.len() as f64);
        let mut bins = vec![0; 21];
        for i in 0..rates.len() {
            bins[(20.0 * rates[i]).floor() as usize] += 1;
        }
        let mut total = 0.0;
        for i in 0..rates.len() {
            total += rates[i];
        }
        println!("\nsubstitution rates");
        for i in 0..bins.len() {
            if bins[i] > 0 {
                println!(
                    "{}-{}% ==> {:.1}%",
                    5 * i,
                    5 * (i + 1),
                    100.0 * bins[i] as f64 / rates.len() as f64
                );
            }
        }
        println!(
            "mean substitution rate = {:.1}%",
            100.0 * total as f64 / rates.len() as f64
        );
        println!(
            "\nDD fraction = {:.2}%\n",
            100.0 * dd as f64 / drefnames.len() as f64
        );
    }
}
