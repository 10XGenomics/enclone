// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build PWMs for IMGT mammalian reference sequences.  For each pair (chain, feature) there
// are one or more matrices.  There are several steps.
//
// 1. For the following lengths, and all sequences (chain, feature, length), directly form a PWM.
//
//        IGH        IGK           IGL
// fwr1   25         23            22
// fwr2   19         15            15
// fwr3   39         32            34
// cdr1   7,11       12,17         14
// cdr2   8          7             7,11
//
// 2. Reference sequences longer than the max length are ignored.
//
// 3. Other reference sequences are aligned to the smallest length that is at least that length,
// and then the PWM is updated.  We ignore sequences that have insertions.

use crate::vdj_features::*;
use amino::*;
use bio::alignment::pairwise::banded::*;
use bio::alignment::AlignmentOperation::Del;
use bio::alignment::AlignmentOperation::Ins;
use debruijn::dna_string::*;
use fasta_tools::*;
use std::cmp::max;
use std::collections::HashMap;
use std::fs::read_dir;
use string_utils::*;
use superslice::Ext;
use vector_utils::*;

pub fn make_mammalian_pwms() -> Vec<(String, String, usize, Vec<Vec<(u32, u8)>>)> {
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

    // Define initial lengths, per the table documented at the top of this file.

    let mut starts = HashMap::<(String, String), Vec<usize>>::new();
    starts.insert(("IGH".to_string(), "fwr1".to_string()), vec![25]);
    starts.insert(("IGK".to_string(), "fwr1".to_string()), vec![23]);
    starts.insert(("IGL".to_string(), "fwr1".to_string()), vec![22]);
    starts.insert(("IGH".to_string(), "fwr2".to_string()), vec![19]);
    starts.insert(("IGK".to_string(), "fwr2".to_string()), vec![15]);
    starts.insert(("IGL".to_string(), "fwr2".to_string()), vec![15]);
    starts.insert(("IGH".to_string(), "fwr3".to_string()), vec![39]);
    starts.insert(("IGK".to_string(), "fwr3".to_string()), vec![32]);
    starts.insert(("IGL".to_string(), "fwr3".to_string()), vec![34]);
    starts.insert(("IGH".to_string(), "cdr1".to_string()), vec![7, 11]);
    starts.insert(("IGK".to_string(), "cdr1".to_string()), vec![12, 17]);
    starts.insert(("IGL".to_string(), "cdr1".to_string()), vec![14]);
    starts.insert(("IGH".to_string(), "cdr2".to_string()), vec![8]);
    starts.insert(("IGK".to_string(), "cdr2".to_string()), vec![7]);
    starts.insert(("IGL".to_string(), "cdr2".to_string()), vec![7, 11]);

    // Set up to track feature sequences.  The structure is:
    // (chain type, feature, length, sequence).

    let mut calls = Vec::<(String, String, usize, Vec<u8>)>::new();

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
            if cdr1.is_none() || cdr2.is_none() {
                continue;
            }
            if fwr1.is_none() || fwr2.is_none() || fwr3.is_none() {
                continue;
            }
            let cdr1 = cdr1.unwrap();
            let cdr2 = cdr2.unwrap();
            let fwr1 = fwr1.unwrap();
            let fwr2 = fwr2.unwrap();
            let fwr3 = fwr3.unwrap();
            if cdr1.is_empty() || cdr2.is_empty() {
                continue;
            }
            if fwr1.is_empty() || fwr2.is_empty() || fwr3.is_empty() {
                continue;
            }
            calls.push((
                chain_type.to_string(),
                "cdr1".to_string(),
                cdr1.len(),
                cdr1.clone(),
            ));
            calls.push((
                chain_type.to_string(),
                "cdr2".to_string(),
                cdr2.len(),
                cdr2.clone(),
            ));
            calls.push((
                chain_type.to_string(),
                "fwr1".to_string(),
                fwr1.len(),
                fwr1.clone(),
            ));
            calls.push((
                chain_type.to_string(),
                "fwr2".to_string(),
                fwr2.len(),
                fwr2.clone(),
            ));
            calls.push((
                chain_type.to_string(),
                "fwr3".to_string(),
                fwr3.len(),
                fwr3.clone(),
            ));
        }
    }
    calls.sort();

    // Form the PWMs.

    let mut pwms = Vec::<(String, String, usize, Vec<Vec<(u32, u8)>>)>::new();
    for s in starts.iter() {
        let chain_type = &s.0 .0;
        let feature = &s.0 .1;
        for (i, len) in s.1.iter().enumerate() {
            // Form an initial PWM.

            let mut pwm = Vec::<Vec<(u32, u8)>>::new();
            let m1 =
                calls.lower_bound_by_key(&(chain_type, feature, len), |(a, b, c, _d)| (a, b, c));
            let m2 =
                calls.upper_bound_by_key(&(chain_type, feature, len), |(a, b, c, _d)| (a, b, c));
            for j in 0..*len {
                let mut c = Vec::<u8>::new();
                for m in m1..m2 {
                    c.push(calls[m].3[j]);
                }
                c.sort_unstable();
                let mut freq = Vec::<(u32, u8)>::new();
                make_freq(&c, &mut freq);
                pwm.push(freq);
            }

            // Find the data that are to be merged into the initial PWM.

            let mut llen = 0;
            if i > 0 {
                llen = s.1[i - 1];
            }
            let hlen = *len - 1;
            let m1 =
                calls.lower_bound_by_key(&(chain_type, feature, &llen), |(a, b, c, _d)| (a, b, c));
            let m2 =
                calls.upper_bound_by_key(&(chain_type, feature, &hlen), |(a, b, c, _d)| (a, b, c));

            // Merge.

            'm_loop: for m in m1..m2 {
                let x1 = &calls[m].3;

                // Define scoring scheme.  This is very ugly, because the second argument b is
                // treated as the column number in the position weight matrix.

                let mut n = 0_i32;
                for i in 0..pwm[0].len() {
                    n += pwm[0][i].0 as i32;
                }
                let score = |a: u8, b: u8| {
                    let mut dot_count = 0;
                    for x in pwm[b as usize].iter() {
                        if x.1 == b'.' {
                            dot_count = x.0 as i32;
                        }
                    }
                    for x in pwm[b as usize].iter() {
                        if a == x.1 {
                            return x.0 as i32 - dot_count;
                        }
                    }
                    -n
                };
                let (gap_open, gap_extend) = (-1 * n as i32, -n as i32);

                // Set up the aligner.

                let mut x2 = Vec::<u8>::new();
                for i in 0..pwm.len() {
                    x2.push(i as u8);
                }
                let (n1, n2) = (x1.len() as u32, x2.len() as u32);
                assert!(n1 <= n2);
                let l = max(1, (n2 - n1) / 2);
                let corners = vec![(0, l), (n1 - 1, n1 - 1 + l)];
                let path = vec![0, 1];
                let mut aligner = Aligner::new(gap_open, gap_extend, &score, 1, l as usize);

                // Align.

                let align = aligner.custom_with_match_path(x1, &x2, &corners, &path);
                let ops = &align.operations;

                // Punt if there was an insertion.

                for i in 0..ops.len() {
                    if ops[i] == Ins {
                        continue 'm_loop;
                    }
                }

                // Update the PWM.

                assert_eq!(ops.len(), x2.len());
                let mut n = 0;
                for i in 0..ops.len() {
                    let c;
                    if ops[i] == Del {
                        c = b'.';
                    } else {
                        c = x1[n];
                        n += 1;
                    }
                    let mut found = false;
                    for j in 0..pwm[i].len() {
                        if pwm[i][j].1 == c {
                            pwm[i][j].0 += 1;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        pwm[i].push((1_u32, c));
                    }
                }
            }
            pwms.push((chain_type.clone(), feature.clone(), *len, pwm));
        }
    }
    pwms.sort();
    pwms
}
