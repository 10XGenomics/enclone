// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Given two PDB antibody/antigen structures, score their similarity.
//
// Run this from the top level of the repo.  This some fragile hardcoding here to identify the
// chains.  This is also currently hardcoded for sars-cov2.
//
// Usage: pdb2 structure-name1 structure-name2 [show]
// or     pdb2 all [show] [nice] [equal] [hiv] [equiv] [seq]
//
// equal: include equal antibodies
// hiv: include hiv antigens
// equiv: show equivalence classes under seq_dist = 0
// seq: show full chain sequences
//
// The second form does an all-vs-all comparison using the entries in
// antibody_sets/antibody_antigen_structures.
//
// Some challenging examples:

// 7B3O                 IGHV3-53    CAR---DVAD----AFDIW          CQQLNSYPPFTF  dist = 16.98
// 7JMP                 IGHV3-53    CARAHVDTAMVESGAFDIW          CCSYAGSSTWVF
//
// 6XC2                 IGHV3-53    CARDLDVYG-LDVW               CQQLNSYPPKFTF dist = 2.57
// 6XCM=6XCN            IGHV3-53    CARGEGWELPYDYW               CSSYEGSNNFVVF
//
//                                     *******  * *                  *   **                   12
// 7BZ5                 IGHV3-53    CAREA-----YGMDVW             CQQLNSYPPYTF  dist = 2.56
// 7KFV                 IGHV3-53    CARGDVSGYRYGLDYW             CQQLISYP-GTF
//
//                                     ******** *                 ** *  *****                 17
// 7CDI                 IGHV3-53    CARDLVVY-GMDVW               CQQY-GSS--PTF dist = 1.57
// 6XCM=6XCN            IGHV3-53    CARGEGWELPYDYW               CSSYEGSNNFVVF

use amino::*;
use debruijn::dna_string::*;
use edit_distance::*;
use enclone_tools::pdb::*;
use equiv::*;
use fasta_tools::*;
use io_utils::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::cmp::min;
use std::env;
use std::io::Write;
use string_utils::*;
use vdj_ann::vdj_features::*;
use vdj_ann_ref::*;
use vector_utils::*;

pub fn cdr3_start_longer(aa: &[u8], _chain_type: &str, _verbose: bool) -> usize {
    let motif = [b"LQPEDSAVYYC", b"VEASQTGTYFC", b"ATSGQASLYLC"];
    let nm = motif[0].len();
    let reach = 400;
    let mut scores = Vec::<(usize, usize)>::new();
    for j in aa.len() as isize - nm as isize - reach as isize..=aa.len() as isize - nm as isize {
        if j < 0 {
            continue;
        }
        let j = j as usize;
        let mut score = 0;
        for k in 0..nm {
            for l in 0..motif.len() {
                if aa[j + k] == motif[l][k] {
                    score += 1;
                    if aa[j + k] == b'Q' {
                        break;
                    }
                }
            }
        }
        scores.push((score, j + nm));
    }
    reverse_sort(&mut scores);
    scores[0].1
}

pub fn fetch_atoms(pdb: &PdbStructure, chain: usize, seq: &[u8]) -> Vec<[f32; 3]> {
    for i in 0..pdb.chains[chain].len() {
        if pdb.chains[chain][i..].starts_with(&seq) {
            let mut u = Vec::<[f32; 3]>::new();
            for j in 0..pdb.atoms.len() {
                if pdb.atoms[j].chain as usize == chain {
                    if (pdb.atoms[j].chain_pos as usize) >= i
                        && (pdb.atoms[j].chain_pos as usize) < i + seq.len()
                    {
                        let x_s = pdb.atoms[j].x;
                        let y_s = pdb.atoms[j].y;
                        let z_s = pdb.atoms[j].z;
                        u.push([x_s, y_s, z_s]);
                    }
                }
            }
            return u;
        }
    }
    eprintln!("fetch_atoms failed");
    std::process::exit(1);
}

fn main() {
    PrettyTrace::new().on();
    if !path_exists("antibody_sets/pdbs") {
        eprintln!("\nYou need to run pdb2.\n");
        std::process::exit(1);
    }
    let args: Vec<String> = env::args().collect();
    let mut pdbs = Vec::<PdbStructure>::new();
    let mut codes = Vec::<String>::new();
    let mut antigens = Vec::<String>::new();
    let mut show = false;
    let mut equal = false;
    let mut hiv = false;
    let mut equiv = false;
    // let mut nice = false;
    let mut seq = false;
    let mut start = 3;
    if args[1] == "all" {
        start = 2;
    }
    for i in start..args.len() {
        if args[i] == "show" {
            show = true;
        } else if args[i] == "nice" {
            // nice = true;
        } else if args[i] == "equal" {
            equal = true;
        } else if args[i] == "hiv" {
            hiv = true;
        } else if args[i] == "equiv" {
            equiv = true;
        } else if args[i] == "seq" {
            seq = true;
        }
    }
    if args[1] == "all" {
        let pdb_list = include_str!("antibody_antigen_structures");
        let mut results = Vec::<(String, PdbStructure)>::new();
        for line in pdb_list.lines() {
            if !line.starts_with('#') && line.len() > 0 {
                if !line.contains("hiv") || hiv {
                    let mut x = line.after(" ").split(',').collect::<Vec<&str>>();
                    if !equal {
                        x.truncate(1);
                    }
                    for p in x.iter() {
                        results.push((p.to_string(), PdbStructure::default()));
                        codes.push(p.to_string());
                        antigens.push(line.before(" ").to_string());
                    }
                }
            }
        }
        results.par_iter_mut().for_each(|res| {
            res.1 = fetch_pdb_structure(&res.0).unwrap();
        });
        for i in 0..results.len() {
            pdbs.push(results[i].1.clone());
        }
    } else {
        for i in 1..=2 {
            pdbs.push(fetch_pdb_structure(&args[i]).unwrap());
            codes.push(args[i].to_string());
            antigens.push("".to_string());
        }
    }

    // Compute amino acids sequences for human reference IGHV segments (FWR1..FWR3 inclusive).

    let mut vs = Vec::<Vec<u8>>::new();
    let mut vs_names = Vec::<String>::new();
    let mut human_count = 0;
    for (j, r) in [&human_ref(), &mouse_ref()].iter().enumerate() {
        if j == 1 {
            human_count = vs.len();
        }
        let mut basesp = Vec::<DnaString>::new();
        let mut headers = Vec::<String>::new();
        read_fasta_contents_into_vec_dna_string_plus_headers(r, &mut basesp, &mut headers);
        for i in 0..basesp.len() {
            let h = &headers[i];
            if h.contains("V-REGION") {
                let mut aa = aa_seq(&basesp[i].to_ascii_vec(), 0);
                let f1 = fr1_start(&aa, "IGH");
                let c3 = cdr3_start(&aa, "IGH", false);
                aa = aa[f1..c3].to_vec();
                vs.push(aa);
                vs_names.push(h.between("|", " ").to_string());
            }
        }
    }

    // Compute distances.  Note duplication with code below.

    let mut distances = vec![(0, vec![Vec::<(Vec<u8>, f32)>::new(); 6]); pdbs.len()];
    for z in 0..pdbs.len() {
        distances[z].0 = z;
    }
    distances.par_iter_mut().for_each(|res| {
        let z = res.0;
        let pdb = &pdbs[z];
        let (mut heavy, mut light) = (None, None);
        for i in 0..pdb.chain_names.len() {
            if pdb.chain_names[i].contains("light")
                || pdb.chain_names[i].contains("Light")
                || pdb.chain_names[i].contains("LIGHT")
                || pdb.chain_names[i].contains("Fab L")
                || pdb.chain_names[i].contains("L chain")
                || pdb.chain_names[i].contains("Kappa")
                || pdb.chain_names[i].contains("kappa")
            {
                light = Some(i as u16);
            } else if pdb.chain_names[i].contains("heavy")
                || pdb.chain_names[i].contains("Heavy")
                || pdb.chain_names[i].contains("HEAVY")
                || pdb.chain_names[i].contains("H chain")
                || pdb.chain_names[i].contains("Fab H")
            {
                heavy = Some(i as u16);
            }
        }
        assert!(heavy.is_some());
        assert!(light.is_some());
        let (light, heavy) = (light.unwrap() as usize, heavy.unwrap() as usize);
        assert!(light != heavy);
        let mut spike = 0;
        for i in 0..pdb.chain_names.len() {
            if i != light && i != heavy {
                spike = i;
            }
        }

        // Find the CDRs.

        let heavy_fwr4_stop = ig_j_gene(&pdb.chains[heavy]).unwrap().1;
        let light_chain_type = ig_j_gene(&pdb.chains[light]).unwrap().0;
        let light_fwr4_stop = ig_j_gene(&pdb.chains[light]).unwrap().1;
        let heavy_thru_cdr3 = &pdb.chains[heavy][0..heavy_fwr4_stop - 10];
        let light_thru_cdr3 = &pdb.chains[light][0..light_fwr4_stop - 9];
        let heavy_cdr3_stop = heavy_thru_cdr3.len();
        let light_cdr3_stop = light_thru_cdr3.len();
        let heavy_cdr3_start = cdr3_start_longer(&heavy_thru_cdr3, "IGH", false) - 1;
        let light_cdr3_start = cdr3_start_longer(&light_thru_cdr3, &light_chain_type, false) - 1;
        let heavy_trim_stop = min(heavy_cdr3_start + 10, heavy_thru_cdr3.len());
        let heavy_trim = format!(
            "MXXXXXXXXXXXXXXXXXXX{}",
            strme(&heavy_thru_cdr3[0..heavy_trim_stop])
        )
        .as_bytes()
        .to_vec();
        let light_trim_stop = min(light_cdr3_start + 10, light_thru_cdr3.len());
        let light_trim = format!(
            "MXXXXXXXXXXXXXXXXXXX{}",
            strme(&light_thru_cdr3[0..light_trim_stop])
        )
        .as_bytes()
        .to_vec();
        let heavy_cdr1_start = cdr1_start(&heavy_trim, "IGH", false).unwrap() - 20;
        let heavy_cdr1_stop = fr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
        let light_cdr1_start = cdr1_start(&light_trim, &light_chain_type, false);
        if heavy_cdr1_start < heavy_cdr1_stop && light_cdr1_start.is_some() {
            let light_cdr1_start = light_cdr1_start.unwrap() - 20;
            let light_cdr1_stop = fr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
            let heavy_cdr2_start = cdr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
            let heavy_cdr2_stop = fr3_start(&heavy_trim, "IGH", false).unwrap() - 20;
            let light_cdr2_start = cdr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
            let light_cdr2_stop = fr3_start(&light_trim, &light_chain_type, false).unwrap() - 20;

            // Get atoms.

            let heavy_cdr1_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr1_start, heavy_cdr1_stop);
            let heavy_cdr2_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr2_start, heavy_cdr2_stop);
            let heavy_cdr3_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr3_start, heavy_cdr3_stop);
            let light_cdr1_atoms = pdb.fetch_atoms_range(light, light_cdr1_start, light_cdr1_stop);
            let light_cdr2_atoms = pdb.fetch_atoms_range(light, light_cdr2_start, light_cdr2_stop);
            let light_cdr3_atoms = pdb.fetch_atoms_range(light, light_cdr3_start, light_cdr3_stop);

            // Gather distances.

            for i in 0..pdb.atoms.len() {
                if pdb.atoms[i].chain as usize == spike {
                    let x_s = pdb.atoms[i].x;
                    let y_s = pdb.atoms[i].y;
                    let z_s = pdb.atoms[i].z;
                    let pos = pdb.atoms[i].chain_pos as usize;
                    let sp = &pdb.chains[spike];
                    const FLANK: usize = 3;
                    if pos >= FLANK && pos + FLANK < sp.len() {
                        let word = &sp[pos - FLANK..=pos + FLANK];
                        for (a, atoms) in [
                            &heavy_cdr1_atoms,
                            &heavy_cdr2_atoms,
                            &heavy_cdr3_atoms,
                            &light_cdr1_atoms,
                            &light_cdr2_atoms,
                            &light_cdr3_atoms,
                        ]
                        .iter()
                        .enumerate()
                        {
                            let mut min_dist = 1_000_000_000.0_f32;
                            for atom in atoms.iter() {
                                let s1 = x_s - atom[0];
                                let s2 = y_s - atom[1];
                                let s3 = z_s - atom[2];
                                let dist = (s1 * s1 + s2 * s2 + s3 * s3).sqrt();
                                min_dist = min_dist.min(dist);
                            }
                            let mut found = false;

                            // This computation is inefficient and unnecessary.

                            for k in 0..res.1[a].len() {
                                if res.1[a][k].0 == word {
                                    res.1[a][k].1 = res.1[a][k].1.min(min_dist);
                                    found = true;
                                    break;
                                }
                            }
                            if !found {
                                res.1[a].push((word.to_vec(), min_dist));
                            }
                        }
                    }
                }
            }
        }
    });

    // Proceed.

    const MAX_DIST: f32 = 5.0;
    let cdrs = ["CDR1H", "CDR2H", "CDR3H", "CDR1L", "CDR2L", "CDR3L"];
    let mut results = Vec::<(usize, Vec<f32>, Vec<f32>, Vec<Vec<u8>>, Vec<(usize, usize)>)>::new();
    for z1 in 0..pdbs.len() {
        results.push((
            z1,
            Vec::new(),
            Vec::new(),
            Vec::<Vec<u8>>::new(),
            Vec::new(),
        ));
    }
    results.par_iter_mut().for_each(|res| {
        let z1 = res.0;
        for z2 in z1 + 1..pdbs.len() {
            if antigens[z2] != antigens[z1] {
                continue;
            }
            let mut heavyx = vec![String::new(); 2];
            let mut lightx = vec![String::new(); 2];
            let mut heavyy = Vec::new();
            let mut hfwr1 = vec![String::new(); 2];
            let mut hcdr1 = vec![String::new(); 2];
            let mut lcdr1 = vec![String::new(); 2];
            let mut hcdr2 = vec![String::new(); 2];
            let mut lcdr2 = vec![String::new(); 2];
            let mut hfwr3 = vec![String::new(); 2];
            let mut lfwr3 = vec![String::new(); 2];
            let mut hcdr3 = vec![String::new(); 2];
            let mut lcdr3 = vec![String::new(); 2];
            let mut hfwr4 = vec![String::new(); 2];
            let mut lfwr4 = vec![String::new(); 2];
            let mut spike_seq = Vec::<Vec<u8>>::new();
            let mut hvrefs = Vec::new();
            let mut hvrefnames = Vec::new();
            let mut hvspecies = Vec::new();
            let mut lvrefs = Vec::new();
            let mut lvrefnames = Vec::new();
            let mut lvspecies = Vec::new();
            for pi in 0..2 {
                let pdb: &PdbStructure;
                let _code: &str;
                if pi == 0 {
                    pdb = &pdbs[z1];
                    _code = &codes[z1];
                } else {
                    pdb = &pdbs[z2];
                    _code = &codes[z2];
                }
                if args[1] != "all" {
                    println!("");
                    for i in 0..pdb.chain_names.len() {
                        println!(
                            "{} = {} (len = {})",
                            i + 1,
                            pdb.chain_names[i],
                            pdb.chains[i].len()
                        );
                    }
                }
                assert_eq!(pdb.chain_names.len(), 3);
                let (mut heavy, mut light) = (None, None);
                for i in 0..pdb.chain_names.len() {
                    if pdb.chain_names[i].contains("light")
                        || pdb.chain_names[i].contains("Light")
                        || pdb.chain_names[i].contains("LIGHT")
                        || pdb.chain_names[i].contains("Fab L")
                        || pdb.chain_names[i].contains("L chain")
                        || pdb.chain_names[i].contains("Kappa")
                        || pdb.chain_names[i].contains("kappa")
                    {
                        light = Some(i as u16);
                    } else if pdb.chain_names[i].contains("heavy")
                        || pdb.chain_names[i].contains("Heavy")
                        || pdb.chain_names[i].contains("HEAVY")
                        || pdb.chain_names[i].contains("H chain")
                        || pdb.chain_names[i].contains("Fab H")
                    {
                        heavy = Some(i as u16);
                    }
                }
                let (light, heavy) = (light.unwrap() as usize, heavy.unwrap() as usize);
                heavyx[pi] = stringme(&pdb.chains[heavy]);
                lightx[pi] = stringme(&pdb.chains[light]);

                // Find some features.

                let heavy_fwr4_stop = ig_j_gene(&pdb.chains[heavy]).unwrap().1;
                let light_chain_type = ig_j_gene(&pdb.chains[light]).unwrap().0;
                let light_fwr4_stop = ig_j_gene(&pdb.chains[light]).unwrap().1;
                let heavy_thru_cdr3 = &pdb.chains[heavy][0..heavy_fwr4_stop - 10];
                let light_thru_cdr3 = &pdb.chains[light][0..light_fwr4_stop - 9];
                let heavy_cdr3_stop = heavy_thru_cdr3.len();
                let light_cdr3_stop = light_thru_cdr3.len();
                let heavy_cdr3_start = cdr3_start_longer(&heavy_thru_cdr3, "IGH", false) - 1;
                let light_cdr3_start =
                    cdr3_start_longer(&light_thru_cdr3, &light_chain_type, false) - 1;

                // Find best V gene match for heavy chain.  First we truncate to fwr1..fwr3.

                let mut ref_id = 0;
                let mut aa = pdb.chains[heavy].clone();
                aa.truncate(heavy_cdr3_start + 1);
                heavyy.push(stringme(&aa));
                let f1 = fr1_start(&aa, "IGH");
                aa = aa[f1..].to_vec();
                let mut r = 1000000;
                for j in 0..vs.len() {
                    let dist = edit_distance(&strme(&aa), &strme(&vs[j]));
                    if dist < r {
                        r = dist;
                        ref_id = j;
                    }
                }
                hvrefs.push(ref_id);
                hvrefnames.push(vs_names[ref_id].clone());
                if ref_id < human_count {
                    hvspecies.push("human".to_string());
                } else {
                    hvspecies.push("mouse".to_string());
                }
                let mut ref_id = 0;
                let mut aa = pdb.chains[light].clone();
                aa.truncate(light_cdr3_start + 1);
                let f1 = fr1_start(&aa, &light_chain_type);
                aa = aa[f1..].to_vec();
                let mut r = 1000000;
                for j in 0..vs.len() {
                    let dist = edit_distance(&strme(&aa), &strme(&vs[j]));
                    if dist < r {
                        r = dist;
                        ref_id = j;
                    }
                }
                lvrefs.push(ref_id);
                lvrefnames.push(vs_names[ref_id].clone());
                if ref_id < human_count {
                    lvspecies.push("human".to_string());
                } else {
                    lvspecies.push("mouse".to_string());
                }

                // Find the CDRs.

                let heavy_trim_stop = min(heavy_cdr3_start + 10, heavy_thru_cdr3.len());
                let heavy_trim = format!(
                    "MXXXXXXXXXXXXXXXXXXX{}",
                    strme(&heavy_thru_cdr3[0..heavy_trim_stop])
                )
                .as_bytes()
                .to_vec();
                let light_trim_stop = min(light_cdr3_start + 10, light_thru_cdr3.len());
                let light_trim = format!(
                    "MXXXXXXXXXXXXXXXXXXX{}",
                    strme(&light_thru_cdr3[0..light_trim_stop])
                )
                .as_bytes()
                .to_vec();
                let heavy_cdr1_start = cdr1_start(&heavy_trim, "IGH", false).unwrap() - 20;
                let heavy_cdr1_stop = fr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
                let light_cdr1_start =
                    cdr1_start(&light_trim, &light_chain_type, false).unwrap() - 20;
                let light_cdr1_stop =
                    fr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
                if heavy_cdr1_start >= heavy_cdr1_stop {
                    continue;
                }

                let heavy_fwr1_start = fr1_start(&pdb.chains[heavy], "IGH");
                let heavy_fwr1;
                if heavy_fwr1_start > heavy_cdr1_start {
                    heavy_fwr1 = Vec::<u8>::new();
                } else {
                    heavy_fwr1 = pdb.chains[heavy][heavy_fwr1_start..heavy_cdr1_start].to_vec();
                }

                let heavy_cdr1 = &pdb.chains[heavy][heavy_cdr1_start..heavy_cdr1_stop];
                let light_cdr1 = &pdb.chains[light][light_cdr1_start..light_cdr1_stop];
                let heavy_cdr2_start = cdr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
                let heavy_cdr2_stop = fr3_start(&heavy_trim, "IGH", false).unwrap() - 20;
                let light_cdr2_start =
                    cdr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
                let light_cdr2_stop =
                    fr3_start(&light_trim, &light_chain_type, false).unwrap() - 20;

                if heavy_cdr2_start > heavy_cdr2_stop {
                    continue;
                }

                let heavy_cdr2 = &pdb.chains[heavy][heavy_cdr2_start..heavy_cdr2_stop];
                let light_cdr2 = &pdb.chains[light][light_cdr2_start..light_cdr2_stop];
                let heavy_cdr3 = &pdb.chains[heavy][heavy_cdr3_start..heavy_cdr3_stop];
                let light_cdr3 = &pdb.chains[light][light_cdr3_start..light_cdr3_stop];
                let heavy_fwr3_start = fr3_start(&heavy_trim, "IGH", false).unwrap() - 20;
                let light_fwr3_start =
                    fr3_start(&light_trim, &light_chain_type, false).unwrap() - 20;
                let heavy_fwr3 = &pdb.chains[heavy][heavy_fwr3_start..heavy_cdr3_start];
                let light_fwr3 = &pdb.chains[light][light_fwr3_start..light_cdr3_start];
                let heavy_fwr4 = &pdb.chains[heavy][heavy_fwr4_stop - 10..heavy_fwr4_stop];
                let light_fwr4 = &pdb.chains[heavy][light_fwr4_stop - 9..light_fwr4_stop];
                if args[1] != "all" {
                    println!("heavy cdr2 = IGH : {}", strme(&heavy_cdr2));
                    println!("light cdr2 = {} : {}", light_chain_type, strme(&light_cdr2));
                    println!("heavy fwr3 = IGH : {}", strme(&heavy_fwr3));
                    println!("heavy cdr3 = IGH : {}", strme(&heavy_cdr3));
                    println!("light fwr3 = {} : {}", light_chain_type, strme(&light_fwr3));
                    println!("light cdr3 = {} : {}", light_chain_type, strme(&light_cdr3));
                    println!("heavy fwr4 = IGH : {}", strme(&heavy_fwr4));
                    println!("light fwr4 = {} : {}", light_chain_type, strme(&light_fwr4));
                }
                let mut spike = 0;
                for j in 0..3 {
                    if j != light && j != heavy {
                        spike = j;
                    }
                }
                spike_seq.push(pdb.chains[spike].clone());
                hfwr1[pi] = stringme(&heavy_fwr1);
                hcdr1[pi] = stringme(&heavy_cdr1);
                lcdr1[pi] = stringme(&light_cdr1);
                hcdr2[pi] = stringme(&heavy_cdr2);
                lcdr2[pi] = stringme(&light_cdr2);
                hfwr3[pi] = stringme(&heavy_fwr3);
                lfwr3[pi] = stringme(&light_fwr3);
                hcdr3[pi] = stringme(&heavy_cdr3);
                lcdr3[pi] = stringme(&light_cdr3);
                hfwr4[pi] = stringme(&heavy_fwr4);
                lfwr4[pi] = stringme(&light_fwr4);
            }

            // Compute edit distance.

            /*
            let ed =
                  2.0 * (edit_distance(&hcdr3[0], &hcdr3[1])
                + edit_distance(&lcdr3[0], &lcdr3[1])) as f32
                  // 10.0 same
                + 10.0 * (lfwr3[0].len() as isize - lfwr3[1].len() as isize).abs() as f32
                  // 1.0 worse, 3.0 worse
                + 2.0 * ((hcdr3[0].len() as isize - hcdr3[1].len() as isize).abs() as f32)
                  // 2.0 same, 4.0 same
                + 3.0 * (hfwr3[0].len() as isize - hfwr3[1].len() as isize).abs() as f32
                  // 1.0 worse, 3.0 same
                + 2.0 * (edit_distance(&hcdr2[0], &hcdr2[1])
                + edit_distance(&lcdr2[0], &lcdr2[1])) as f32
                  // 0.5 same, 2.0 worse
                + 1.0 * (edit_distance(&hfwr3[0], &hfwr3[1])
                + edit_distance(&lfwr3[0], &lfwr3[1])) as f32;
            */

            let ed = edit_distance(&hcdr3[0], &hcdr3[1]) + edit_distance(&lcdr3[0], &lcdr3[1]);

            let fed = edit_distance(&hfwr3[0], &hfwr3[1]) + edit_distance(&lfwr3[0], &lfwr3[1]);
            /*
            let hed = edit_distance(&hcdr3[0], &hcdr3[1])
                + 2 * ((hcdr3[0].len() as isize - hcdr3[1].len() as isize).abs() as usize);
            let xed = edit_distance(&hfwr4[0], &hfwr4[1]) + edit_distance(&lfwr4[0], &lfwr4[1]);
            let c1d = edit_distance(&hcdr1[0], &hcdr1[1]) + edit_distance(&lcdr1[0], &lcdr1[1]);
            let c2d = edit_distance(&hcdr2[0], &hcdr2[1]) + edit_distance(&lcdr2[0], &lcdr2[1]);
            */

            /*
            let gap_pen = -10;
            let x1 = &hcdr3[0].as_bytes();
            let x2 = &hcdr3[1].as_bytes();
            let mut aligner = Aligner::with_capacity(x1.len(), x2.len(), 0, gap_pen, &pam120);
            let mut ed = -aligner.global(&x1, &x2).score;
            let x1 = &lcdr3[0].as_bytes();
            let x2 = &lcdr3[1].as_bytes();
            let mut aligner = Aligner::with_capacity(x1.len(), x2.len(), 0, gap_pen, &pam120);
            ed -= aligner.global(&x1, &x2).score;
            */

            /*
            let score = |a: u8, b: u8| {
                if a == b {
                    0i32
                } else if (a == b'L' && b == b'I') || (a == b'I' && b == b'L') {
                    -100i32
                } else if (a == b'L' && b == b'V') || (a == b'V' && b == b'L') {
                    -100i32
                } else if (a == b'I' && b == b'V') || (a == b'V' && b == b'I') {
                    -100i32
                } else {
                    -100i32
                }
            };
            let x1 = &hcdr3[0].as_bytes();
            let x2 = &hcdr3[1].as_bytes();
            let mut aligner = Aligner::with_capacity(x1.len(), x2.len(), 0, -110, &score);
            let alignment = aligner.global(x1, x2);
            let mut ed = -alignment.score;
            let x1 = &lcdr3[0].as_bytes();
            let x2 = &lcdr3[1].as_bytes();
            let mut aligner = Aligner::with_capacity(x1.len(), x2.len(), 0, -110, &score);
            let alignment = aligner.global(x1, x2);
            ed -= alignment.score;
            */

            // Compare distances.

            if args[1] != "all" && show {
                println!("");
            }
            let mut sum = 0.0_f32;
            let mut n = 0;
            let mut lines = Vec::<(usize, String)>::new();
            for i in 0..cdrs.len() {
                for k1 in 0..distances[z1].1[i].len() {
                    for k2 in 0..distances[z2].1[i].len() {
                        if distances[z1].1[i][k1].0 == distances[z2].1[i][k2].0 {
                            let d1 = distances[z1].1[i][k1].1;
                            let d2 = distances[z2].1[i][k2].1;
                            if d1 <= MAX_DIST || d2 <= MAX_DIST {
                                if d1 < 1_000_000_000.0 && d2 < 1_000_000_000.0 {
                                    n += 1;
                                    let d = d1 - d2;
                                    sum += d * d;
                                    if args[1] != "all" && show {
                                        let x = &distances[z1].1[i][k1].0;
                                        let mut pos = 0;
                                        for j in 0..spike_seq[0].len() {
                                            if spike_seq[0][j..].starts_with(&x) {
                                                pos = j + 3;
                                            }
                                        }
                                        let mut log = Vec::<u8>::new();
                                        fwriteln!(
                                            log,
                                            "{} at distances {:.1} and {:.1} from {}",
                                            pos,
                                            d1,
                                            d2,
                                            cdrs[i]
                                        );
                                        lines.push((pos, stringme(&log)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            lines.sort();
            for j in 0..lines.len() {
                print!("{}", lines[j].1);
            }
            let dist = sum.sqrt() / (n as f32).sqrt();

            // Compute seq_dist.

            let mut seq_dist = 0;

            if edit_distance(&hcdr3[0], &hcdr3[1]) >= 10 {
                seq_dist += 10;
            }
            const MAX_REF_DIFFS: usize = 0;
            let ref_diffs = edit_distance(&strme(&vs[hvrefs[0]]), &strme(&vs[hvrefs[1]]));
            let hdiffs = edit_distance(&heavyy[0], &heavyy[1]);
            if ref_diffs > MAX_REF_DIFFS && hdiffs > 1 && edit_distance(&hcdr3[0], &hcdr3[1]) >= 1 {
                seq_dist += 10;
            }

            /*
            let mut exception = false;
            let trim = 2;
            if hcdr3[0].len() == hcdr3[1].len() && hcdr3[0].len() >= 2*trim {
                let mut max_match = 0;
                let x1 = &hcdr3[0].as_bytes();
                let x2 = &hcdr3[1].as_bytes();
                let mut j = trim;
                while j < x1.len() - trim {
                    if x1[j] != x2[j] {
                        j += 1;
                    } else {
                        let mut k = j;
                        while k < x1.len() && x1[k] == x2[k] {
                            k += 1;
                        }
                        max_match = max(max_match, k - j);
                        j = k;
                    }
                    if max_match >= 5 {
                        exception = true;
                    }
                }
                if exception && edit_distance(&hcdr3[0], &hcdr3[1]) < 10 {
                    seq_dist = 0;
                }
            }
            */

            /*
            if ed > 15 {
                seq_dist += ed;
            }
            if edit_distance(&hfwr1[0], &hfwr1[1]) + edit_distance(&hfwr3[0], &hfwr3[1]) > 5 {
                seq_dist += 10;
            }
            */

            if seq_dist == 0 {
                res.4.push((z1, z2));
            }

            // Flag.

            const STRUCT_DIST_CUTOFF: f32 = 3.50;
            let flag = (dist < STRUCT_DIST_CUTOFF) == (seq_dist > 0);

            // Save.

            if !dist.is_nan() {
                // nan case very rare
                res.1.push(seq_dist as f32);
                res.2.push(dist);
                let mut log = Vec::<u8>::new();
                if args[1] != "all" {
                    fwriteln!(log, "");
                }
                if antigens[z1].len() > 0 {
                    fwrite!(log, "{} ", antigens[z1]);
                }
                fwrite!(
                    log,
                    "{} {} {} {}",
                    hvspecies[0],
                    lvspecies[0],
                    codes[z1],
                    codes[z2]
                );
                fwrite!(log, ", {}/{}", hvrefnames[0], hvrefnames[1]);
                fwrite!(log, ", {}/{}", lvrefnames[0], lvrefnames[1]);
                fwrite!(log, ", struct d = {:.2}", dist);
                fwrite!(log, ", seq d = {}", seq_dist);
                fwrite!(log, ", ed = {}", ed);
                // fwrite!(log, ", cdr3 d = {}", ed);
                fwrite!(log, ", fwr3 d = {}", fed);
                // fwrite!(log, ", fwr4 d = {}", xed);
                // fwrite!(log, ", cdr1 d = {}", c1d);
                // fwrite!(log, ", cdr2 d = {}", c2d);
                /*
                fwrite!(log, ", hcdr3 ld = {}",
                    (hcdr3[0].len() as isize -hcdr3[1].len() as isize ).abs());
                fwrite!(log, ", hed = {}", hed);
                */
                if flag {
                    fwrite!(log, " [FLAG]");
                }
                fwriteln!(log, "");
                fwriteln!(log, "FWR1  {}", hfwr1[0]);
                fwriteln!(log, "FWR1  {}", hfwr1[1]);
                if hfwr1[0].len() == hfwr1[1].len() {
                    fwrite!(log, "      ");
                    let f1 = &hfwr1[0].as_bytes();
                    let f2 = &hfwr1[1].as_bytes();
                    for j in 0..f1.len() {
                        if f1[j] != f2[j] {
                            fwrite!(log, "*");
                        } else {
                            fwrite!(log, " ");
                        }
                    }
                    fwriteln!(log, "");
                }
                fwriteln!(log, "CDR3  {}  {}", hcdr3[0], lcdr3[0]);
                fwriteln!(log, "CDR3  {}  {}", hcdr3[1], lcdr3[1]);
                fwriteln!(log, "FWR3  {}  {}", hfwr3[0], lfwr3[0]);
                fwriteln!(log, "FWR3  {}  {}", hfwr3[1], lfwr3[1]);
                if hfwr3[0].len() == hfwr3[1].len() {
                    fwrite!(log, "      ");
                    let f1 = &hfwr3[0].as_bytes();
                    let f2 = &hfwr3[1].as_bytes();
                    for j in 0..f1.len() {
                        if f1[j] != f2[j] {
                            fwrite!(log, "*");
                        } else {
                            fwrite!(log, " ");
                        }
                    }
                    if lfwr3[0].len() == lfwr3[1].len() {
                        fwrite!(log, "  ");
                        let f1 = &lfwr3[0].as_bytes();
                        let f2 = &lfwr3[1].as_bytes();
                        for j in 0..f1.len() {
                            if f1[j] != f2[j] {
                                fwrite!(log, "*");
                            } else {
                                fwrite!(log, " ");
                            }
                        }
                    }
                    fwriteln!(log, "");
                }
                if (flag && dist <= STRUCT_DIST_CUTOFF) || seq {
                    fwriteln!(log, "heavy1 = {}", heavyx[0]);
                    fwriteln!(log, "light1 = {}", lightx[0]);
                    fwriteln!(log, "heavy2 = {}", heavyx[1]);
                    fwriteln!(log, "light2 = {}", lightx[1]);
                }
                fwriteln!(log, "");
                res.3.push(log);
            }
        }
    });
    let mut flags = 0;
    let mut dists_logs = Vec::new();
    for i in 0..results.len() {
        for j in 0..results[i].2.len() {
            if strme(&results[i].3[j]).contains("FLAG") {
                flags += 1;
            }
            dists_logs.push((results[i].2[j], results[i].3[j].clone()));
        }
    }
    if args[1] != "all" || show {
        dists_logs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for i in 0..dists_logs.len() {
            print!("{}", strme(&dists_logs[i].1));
        }
        if args[1] != "all" {
            std::process::exit(0);
        }
    }
    if equiv {
        let mut joins = Vec::<(usize, usize)>::new();
        for i in 0..results.len() {
            joins.append(&mut results[i].4.clone());
        }
        let mut e = EquivRel::new(pdbs.len() as i32);
        for i in 0..joins.len() {
            e.join(joins[i].0 as i32, joins[i].1 as i32);
        }
        let mut reps = Vec::<i32>::new();
        e.orbit_reps(&mut reps);
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            e.orbit(reps[i], &mut o);
            if o.len() >= 2 {
                print!("{} orbit:", antigens[o[0] as usize]);
                for j in 0..o.len() {
                    print!(" {}", codes[o[j] as usize]);
                }
                println!("");
            }
        }
    }

    // Start analyzing the results.

    let mut x = Vec::<(f32, f32)>::new(); // {(seq dist, struct dist)}
    let mut struct_dists = Vec::<f32>::new();
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            x.push((results[i].1[j], results[i].2[j]));
            struct_dists.push(results[i].2[j]);
        }
    }
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let all = x.len();
    println!("number of distance values = {}", all);
    const EVAL_PC: f64 = 10.0;
    let using = (x.len() as f64 * EVAL_PC / 100.0).round() as usize;

    // Find the fraction of entries that we're using that have zero sequence distance.

    let mut using_zero = 0;
    for i in 0..using {
        if x[i].0 > 0.0 {
            break;
        }
        using_zero += 1;
    }
    let using_zero_pc = 100.0 * using_zero as f64 / using as f64;
    println!(
        "percent of top sequence distances that are zero = {:.1}",
        using_zero_pc
    );

    // Select entries that we're using.  This is tricky because of the last group of entries that
    // have identical sequence distances, but for which the structure distances are ordered.  If
    // we're not careful, we'll bias towards lower structure distances in that group.

    let mut y = Vec::<f32>::new();
    let mut i = 0;
    let mut count = 0;
    'select: while i < all {
        let mut j = i + 1;
        while j < all {
            if x[j].0 != x[i].0 {
                break;
            }
            j += 1;
        }
        let mut mean_struct_dist = 0.0;
        for k in i..j {
            mean_struct_dist += x[k].1;
        }
        mean_struct_dist /= (j - i) as f32;
        let mut z = Vec::<(f32, f32)>::new();
        for k in i..j {
            z.push(((x[k].1 - mean_struct_dist).abs(), x[k].1));
        }
        z.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for k in 0..z.len() {
            y.push(z[k].1);
            count += 1;
            if count == using {
                break 'select;
            }
        }
        i = j;
    }

    // Compute mean rank.

    struct_dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut ranks = Vec::<usize>::new();
    for i in 0..y.len() {
        let mut j = 0;
        while j < struct_dists.len() {
            let mut k = j + 1;
            while k < struct_dists.len() && struct_dists[k] == struct_dists[j] {
                k += 1;
            }
            if y[i] == struct_dists[j] {
                ranks.push((j + k - 1) / 2);
                break;
            }
            j = k;
        }
    }
    assert_eq!(ranks.len(), using);
    let mut mean_rank = 0.0;
    for i in 0..ranks.len() {
        mean_rank += ranks[i] as f64;
    }
    mean_rank /= using as f64;

    // Compute normalized mean rank.  Informal explanation: Let's say we're using the top 10% of
    // antibody pairs, as ordered by sequence distance, and that there are n of those.  Each of
    // those has a rank in the ordered list of all structure distances, and the mean rank is the
    // mean of those ranks.  The best possible mean rank is n/2.  The "random mean rank" is N/2,
    // where N is the total number of antibody pairs.  Then the normalized mean rank is the mean
    // rank scaled so that the best possible value is 0 and the random value is 100.  (One could
    // do worse than that.)

    let optimal_mean_rank = using as f64 / 2.0;
    let random_rank = all as f64 / 2.0;
    let mut normalized_mean_rank = mean_rank - optimal_mean_rank;
    normalized_mean_rank = 100.0 * normalized_mean_rank / (random_rank - optimal_mean_rank);
    println!("normalized mean rank = {:.1}", normalized_mean_rank);
    println!("total flags = {}", flags);
}
