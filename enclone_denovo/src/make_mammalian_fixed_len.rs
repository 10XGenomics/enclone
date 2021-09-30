// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// From IMGT mammalian reference sequences, find the
// (per chain, per feature, per length, per position) amino acid distribution.

use crate::vdj_features::{cdr1, cdr2, cdr3_score, fwr1, fwr2, fwr3};
use amino::aa_seq;
use debruijn::dna_string::DnaString;
use fasta_tools::read_fasta_into_vec_dna_string_plus_headers;
use std::fs::read_dir;
use string_utils::TextUtils;
use vector_utils::make_freq;

pub fn make_mammalian_fixed_len() -> Vec<(String, String, usize, Vec<Vec<(u32, u8)>>)> {
    // Set up to track calls.

    let mut calls = Vec::<(String, String, usize, usize, u8)>::new();

    // Define reference sequence data.

    let mut refs_all = Vec::<Vec<DnaString>>::new();
    let mut headers_all = Vec::<Vec<String>>::new();

    // Gather the reference sequences.

    let all = read_dir(&"v_segments").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        if !f.contains("-mkvdjref") {
            continue;
        }
        let species = f.between("/", "-mkvdjref");

        // Skip fish.

        if species == "Salmo_salar" || species == "Oncorhynchus_mykiss" || species == "Danio_rerio"
        {
            continue;
        }

        // Skip birds.

        if species == "Gallus_gallus" {
            continue;
        }

        let mut refs = Vec::<DnaString>::new();
        let mut headers = Vec::<String>::new();
        read_fasta_into_vec_dna_string_plus_headers(&f.to_string(), &mut refs, &mut headers);

        // Skip broken references: Canus_lupus_familiaris and a bunch of others.

        if headers.is_empty() {
            continue;
        }
        refs_all.push(refs);
        headers_all.push(headers);
    }

    // Go through all the species.

    for m in 0..refs_all.len() {
        let refs = &refs_all[m];
        let headers = &headers_all[m];

        // Go through the reference sequences for the species.

        for i in 0..refs.len() {
            // Only look at V genes.

            if !headers[i].contains("V-REGION") {
                continue;
            }
            if headers[i].contains("TRD") || headers[i].contains("TRG") {
                continue;
            }
            let seq = refs[i].to_ascii_vec();
            let aa = aa_seq(&seq, 0);
            let chain_type = headers[i].after("REGION|").between("|", "|");

            // Exclude junk for the non-10x references.

            if aa.len() < 100 {
                continue;
            }
            let stop = aa.contains(&b'*');
            if stop {
                continue;
            }
            if aa[0] != b'M' {
                continue;
            }
            if aa.len() >= 30 {
                let mut aap = aa.clone();
                aap.push(b'C');
                if cdr3_score(&aap, chain_type, false) > 4 + cdr3_score(&aa, chain_type, false) {
                    continue;
                }
            }
            /*
            let mut fixable = false;
            const TRIM: usize = 10;
            for j in 0..aa.len() - TRIM {
                if aa[j] == b'*' {
                    let mut seqx = seq.clone();
                    for _ in 1..=2 {
                        let _ = seqx.remove(3 * j);
                        let aax = aa_seq(&seqx, 0);
                        if !aax.contains(&b'*') {
                            fixable = true;
                        }
                    }
                }
            }
            if fixable {
                continue;
            }
            */
            if aa.len() < 105 && chain_type == "IGH" {
                continue;
            }
            if aa.len() >= 31 {
                // Pretty crappy frameshift test.  One should see high aa and dna similarity
                // to other seqs if shifted.  Or use more aas.
                let score = cdr3_score(&aa, chain_type, false);
                let mut frameshift = false;
                for del in 1..=2 {
                    let aad = aa_seq(&seq, del);
                    if score <= 6 && cdr3_score(&aad, chain_type, false) >= 3 + score {
                        // println!("frameshift = {} = {}", species, headers[i].before(" "));
                        // use io_utils::*;
                        // printme!(cdr3_score(&aa, &chain_type, false));
                        frameshift = true;
                    }
                }
                if frameshift {
                    continue;
                }
            }

            // Gather calls.

            let cdr1 = cdr1(&aa, chain_type, false);
            let cdr2 = cdr2(&aa, chain_type, false);
            let fwr1 = fwr1(&aa, chain_type, false);
            let fwr2 = fwr2(&aa, chain_type, false);
            let fwr3 = fwr3(&aa, chain_type, false);
            if cdr1.is_some() {
                let cdr1 = cdr1.unwrap();
                for i in 0..cdr1.len() {
                    let ct = chain_type.to_string();
                    calls.push((ct, "cdr1".to_string(), cdr1.len(), i, cdr1[i]));
                }
            }
            if cdr2.is_some() {
                let cdr2 = cdr2.unwrap();
                for i in 0..cdr2.len() {
                    let ct = chain_type.to_string();
                    calls.push((ct, "cdr2".to_string(), cdr2.len(), i, cdr2[i]));
                }
            }
            if fwr1.is_some() {
                let fwr1 = fwr1.unwrap();
                for i in 0..fwr1.len() {
                    let ct = chain_type.to_string();
                    calls.push((ct, "fwr1".to_string(), fwr1.len(), i, fwr1[i]));
                }
            }
            if fwr2.is_some() {
                let fwr2 = fwr2.unwrap();
                for i in 0..fwr2.len() {
                    let ct = chain_type.to_string();
                    calls.push((ct, "fwr2".to_string(), fwr2.len(), i, fwr2[i]));
                }
            }
            if fwr3.is_some() {
                let fwr3 = fwr3.unwrap();
                for i in 0..fwr3.len() {
                    let ct = chain_type.to_string();
                    calls.push((ct, "fwr3".to_string(), fwr3.len(), i, fwr3[i]));
                }
            }
        }
    }

    // Organize calls.

    let mut c2 = Vec::<(String, String, usize, Vec<(u32, u8)>)>::new();
    calls.sort();
    let mut i = 0;
    while i < calls.len() {
        let mut j = i + 1;
        while j < calls.len() {
            if calls[j].0 != calls[i].0 || calls[j].1 != calls[i].1 || calls[j].2 != calls[i].2 {
                break;
            }
            if calls[j].3 != calls[i].3 {
                break;
            }
            j += 1;
        }
        let mut c = Vec::<u8>::new();
        for k in i..j {
            c.push(calls[k].4);
        }
        let mut freq = Vec::<(u32, u8)>::new();
        make_freq(&c, &mut freq);
        c2.push((calls[i].0.clone(), calls[i].1.clone(), calls[i].2, freq));
        i = j;
    }
    let mut x = Vec::<(String, String, usize, Vec<Vec<(u32, u8)>>)>::new();
    let mut i = 0;
    while i < c2.len() {
        let mut j = i + 1;
        while j < c2.len() {
            if c2[j].0 != c2[i].0 || c2[j].1 != c2[i].1 || c2[j].2 != c2[i].2 {
                break;
            }
            j += 1;
        }
        let mut y = Vec::<Vec<(u32, u8)>>::new();
        for k in i..j {
            y.push(c2[k].3.clone());
        }
        x.push((c2[i].0.clone(), c2[i].1.clone(), c2[i].2, y));
        i = j;
    }
    x
}
