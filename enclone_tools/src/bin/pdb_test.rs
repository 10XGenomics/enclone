// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Run this from the top level of the repo.

use enclone_tools::pdb::*;
use pretty_trace::*;
use std::env;
use string_utils::*;
use vdj_ann::vdj_features::*;
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

pub fn mean_dist(a: &Vec<[f32; 3]>, b: &Vec<[f32; 3]>) -> f32 {
    let count = a.len() * b.len();
    let mut sum = 0.0;
    for i in 0..a.len() {
        for j in 0..b.len() {
            let m1 = a[i][0] - b[j][0];
            let m2 = a[i][1] - b[j][1];
            let m3 = a[i][2] - b[j][2];
            sum += (m1 * m1 + m2 * m2 + m3 * m3).sqrt();
        }
    }
    sum / count as f32
}

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let code = &args[1];
    let pdb = fetch_pdb_structure(&code).unwrap();
    for i in 0..pdb.chain_names.len() {
        println!(
            "{} = {} (len = {})",
            i + 1,
            pdb.chain_names[i],
            pdb.chains[i].len()
        );
    }
    println!("{} atoms", pdb.atoms.len());
    assert_eq!(pdb.chain_names.len(), 3);
    let mut spike = None;
    let mut heavy = None;
    let mut light = None;
    for i in 0..pdb.chain_names.len() {
        if pdb.chain_names[i].contains("spike")
            || pdb.chain_names[i].contains("Spike")
            || pdb.chain_names[i].contains("Receptor Binding Domain")
            || pdb.chain_names[i].contains("receptor binding domain")
            || pdb.chain_names[i].contains("RBD")
            || pdb.chain_names[i].contains("glycoprotein")
        {
            spike = Some(i as u16);
        } else if pdb.chain_names[i].contains("light")
            || pdb.chain_names[i].contains("Light")
            || pdb.chain_names[i].contains("Fab L")
            || pdb.chain_names[i].contains("Kappa")
        {
            light = Some(i as u16);
        } else if pdb.chain_names[i].contains("heavy")
            || pdb.chain_names[i].contains("Heavy")
            || pdb.chain_names[i].contains("Fab H")
        {
            heavy = Some(i as u16);
        }
    }
    assert!(spike.is_some());
    assert!(heavy.is_some());
    assert!(light.is_some());
    let spike = spike.unwrap() as usize;
    let light = light.unwrap() as usize;
    let heavy = heavy.unwrap() as usize;
    // println!("spike = {}", strme(&pdb.chains[spike]));

    // Find EGFN in the spike.

    let mut egfn = false;
    let mut noffset = 0;
    for i in 0..pdb.chains[spike].len() {
        if pdb.chains[spike][i..].starts_with(b"EGFN") {
            egfn = true;
            println!("EGFN at {}", i);
            assert!(i <= 600);
            noffset = (600 - i) as u16;
        }
    }
    if !egfn {
        eprintln!(
            "\nFailed to find EGFN in spike:\n{}.\n",
            strme(&pdb.chains[spike])
        );
        std::process::exit(1);
    }

    // Find SASFSTF in the spike.

    for i in 0..pdb.chains[spike].len() {
        if pdb.chains[spike][i..].starts_with(b"SASFSTF") {
            println!("SASFSTF at {}", i);
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
    let heavy_trim = format!(
        "MXXXXXXXXXXXXXXXXXXX{}",
        strme(&heavy_thru_cdr3[0..heavy_cdr3_start + 10])
    )
    .as_bytes()
    .to_vec();
    let light_trim = format!(
        "MXXXXXXXXXXXXXXXXXXX{}",
        strme(&light_thru_cdr3[0..light_cdr3_start + 10])
    )
    .as_bytes()
    .to_vec();
    let heavy_cdr1_start = cdr1_start(&heavy_trim, "IGH", false).unwrap() - 20;
    let heavy_cdr1_stop = fr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
    let light_cdr1_start = cdr1_start(&light_trim, &light_chain_type, false).unwrap() - 20;
    let light_cdr1_stop = fr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
    let heavy_cdr2_start = cdr2_start(&heavy_trim, "IGH", false).unwrap() - 20;
    let heavy_cdr2_stop = fr3_start(&heavy_trim, "IGH", false).unwrap() - 20;
    let light_cdr2_start = cdr2_start(&light_trim, &light_chain_type, false).unwrap() - 20;
    let light_cdr2_stop = fr3_start(&light_trim, &light_chain_type, false).unwrap() - 20;

    println!(
        "heavy cdr1 = {}",
        strme(&pdb.chains[heavy][heavy_cdr1_start..heavy_cdr1_stop])
    );
    println!(
        "heavy cdr2 = {}",
        strme(&pdb.chains[heavy][heavy_cdr2_start..heavy_cdr2_stop])
    );
    let heavy_cdr3 = &pdb.chains[heavy][heavy_cdr3_start..heavy_cdr3_stop];
    println!("heavy cdr3 = {}", strme(&heavy_cdr3));
    println!("light chain type = {}", light_chain_type);
    println!(
        "light cdr1 = {}",
        strme(&pdb.chains[light][light_cdr1_start..light_cdr1_stop])
    );
    println!(
        "light cdr2 = {}",
        strme(&pdb.chains[light][light_cdr2_start..light_cdr2_stop])
    );
    let light_cdr3 = &pdb.chains[light][light_cdr3_start..light_cdr3_stop];
    println!("light cdr3 = {}", strme(&light_cdr3));

    // Get atoms.

    let heavy_cdr1_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr1_start, heavy_cdr1_stop);
    let heavy_cdr2_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr2_start, heavy_cdr2_stop);
    let heavy_cdr3_atoms = pdb.fetch_atoms_range(heavy, heavy_cdr3_start, heavy_cdr3_stop);
    let light_cdr1_atoms = pdb.fetch_atoms_range(light, light_cdr1_start, light_cdr1_stop);
    let light_cdr2_atoms = pdb.fetch_atoms_range(light, light_cdr2_start, light_cdr2_stop);
    let light_cdr3_atoms = pdb.fetch_atoms_range(light, light_cdr3_start, light_cdr3_stop);

    // Next.

    let cdrs = ["CDR1H", "CDR2H", "CDR3H", "CDR1L", "CDR2L", "CDR3L"];
    let mut nears = Vec::<(usize, usize, f32)>::new();
    for i in 0..pdb.atoms.len() {
        if pdb.atoms[i].chain as usize == spike {
            let x_s = pdb.atoms[i].x;
            let y_s = pdb.atoms[i].y;
            let z_s = pdb.atoms[i].z;
            let pos = pdb.atoms[i].chain_pos;
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
                if min_dist <= 5.0 {
                    nears.push(((pos + noffset) as usize, a, min_dist));
                }
            }
        }
    }
    nears.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut i = 0;
    while i < nears.len() {
        let mut j = i + 1;
        while j < nears.len() {
            if nears[j].0 != nears[i].0 || nears[j].1 != nears[i].1 {
                break;
            }
            j += 1;
        }
        println!(
            "{} at distance {:.1} from {}",
            nears[i].0, nears[i].2, cdrs[nears[i].1]
        );
        i = j;
    }

    // Proceed.

    let mut sp = Vec::<u16>::new();
    for i in 0..pdb.atoms.len() {
        if pdb.atoms[i].chain as usize == spike {
            let x_s = pdb.atoms[i].x;
            let y_s = pdb.atoms[i].y;
            let z_s = pdb.atoms[i].z;
            let mut hdist = 1_000_000_000.0_f32;
            let mut ldist = 1_000_000_000.0_f32;
            for c in [light, heavy].iter() {
                let mut dist = 1_000_000_000.0_f32;
                for j in 0..pdb.atoms.len() {
                    if pdb.atoms[j].chain as usize == *c {
                        if *c == heavy {
                            if pdb.atoms[j].chain_pos < heavy_cdr3_start as u16
                                || pdb.atoms[j].chain_pos > heavy_cdr3_stop as u16
                            {
                                continue;
                            }
                        }
                        if *c == light {
                            if pdb.atoms[j].chain_pos < light_cdr3_start as u16
                                || pdb.atoms[j].chain_pos > light_cdr3_stop as u16
                            {
                                continue;
                            }
                        }
                        let x_c = pdb.atoms[j].x;
                        let y_c = pdb.atoms[j].y;
                        let z_c = pdb.atoms[j].z;
                        let u1 = x_s - x_c;
                        let u2 = y_s - y_c;
                        let u3 = z_s - z_c;
                        let d = (u1 * u1 + u2 * u2 + u3 * u3).sqrt();
                        dist = dist.min(d);
                    }
                }
                if *c == heavy {
                    hdist = dist;
                } else {
                    ldist = dist;
                }
            }
            if hdist <= 8.0 && ldist <= 8.0 {
                sp.push(pdb.atoms[i].chain_pos);
            }
        }
    }
    unique_sort(&mut sp);
    for i in 0..sp.len() {
        if i > 0 {
            print!(",");
        }
        print!(
            "{}({})",
            sp[i] + noffset,
            pdb.chains[spike][sp[i] as usize] as char
        );
    }
    println!("");

    // let s1 = fetch_atoms(&pdb, spike, b"EGFN");
    // let s2 = fetch_atoms(&pdb, spike, b"SASFSTF");
    let l = fetch_atoms(&pdb, light, &light_cdr3);
    let h = fetch_atoms(&pdb, heavy, &heavy_cdr3);
    println!(
        "mean distance from heavy CDR3 to light CDR3 = {:.1}",
        mean_dist(&h, &l)
    );
}
