// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Compute FWR3 for all IMGT mammalian reference sequences, and report out a per-chain
// amino acid frequency table.

use crate::vdj_features::{cdr3_score, cdr3_start, fr3_start};
use amino::aa_seq;
use debruijn::dna_string::DnaString;
use fasta_tools::read_fasta_into_vec_dna_string_plus_headers;
use std::fs::read_dir;
use string_utils::TextUtils;
// use vdj_ann::refx::*;
use vector_utils::make_freq;

pub fn make_fwr3_freqs() -> Vec<Vec<Vec<(u32, u8)>>> {
    // Define constants.

    const MOTIF: usize = 30;

    // Set up to track calls.

    let mut calls = vec![vec![Vec::<u8>::new(); MOTIF]; 5];

    // Define reference sequence data.

    let mut refs_all = Vec::<Vec<DnaString>>::new();
    let mut headers_all = Vec::<Vec<String>>::new();
    let mut spx = Vec::<String>::new();

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
        spx.push(species.to_string());
    }
    /*
    for pass in 1..=2 {
        let refx;
        let species;
        if pass == 1 {
            refx = human_ref();
            species = "human";
        } else {
            refx = mouse_ref();
            species = "mouse";
        }
        let mut refs = Vec::<DnaString>::new();
        let mut headers = Vec::<String>::new();
        read_fasta_contents_into_vec_dna_string_plus_headers(&refx, &mut refs, &mut headers);
        let mut refs2 = Vec::<DnaString>::new();
        let mut headers2 = Vec::<String>::new();
        for i in 0..refs.len() {
            if headers[i].contains("36|IGHV1-12") {
                continue;
            }
            refs2.push(refs[i].clone());
            headers2.push(headers[i].clone());
        }
        refs_all.push(refs2);
        headers_all.push(headers2);
        spx.push(species.to_string());
    }
    */

    // Go through all the species.

    for m in 0..refs_all.len() {
        let refs = &refs_all[m];
        let headers = &headers_all[m];
        let species = &spx[m];

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

            // Exclude junk for the non-10x references.  (MOOT NOW.)

            if species != "human" && species != "mouse" {
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
                if aa.len() >= MOTIF {
                    let mut aap = aa.clone();
                    aap.push(b'C');
                    if cdr3_score(&aap, chain_type, false) > 4 + cdr3_score(&aa, chain_type, false)
                    {
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
            }

            // Gather calls.

            let fr3_start = fr3_start(&aa, chain_type, false).unwrap();
            let cdr3_start = cdr3_start(&aa, chain_type, false);
            let f3 = &aa[fr3_start..cdr3_start];
            let r;
            if chain_type == "IGH" {
                r = 0;
            } else if chain_type == "IGK" {
                r = 1;
            } else if chain_type == "IGL" {
                r = 2;
            } else if chain_type == "TRA" {
                r = 3;
            } else {
                assert_eq!(chain_type, "TRB");
                r = 4;
            }

            // Save data.

            let n = f3.len();
            if n >= MOTIF {
                for j in 0..MOTIF {
                    calls[r][j].push(f3[n - j - 1]);
                }
            }
        }
    }

    // Organize the calls.

    let mut freqs = vec![vec![Vec::<(u32, u8)>::new(); MOTIF]; 5];
    for r in 0..5 {
        for j in 0..MOTIF {
            calls[r][j].sort_unstable();
            make_freq(&calls[r][j], &mut freqs[r][j]);
        }
    }

    // Done.

    freqs
}
