// Copyright (c) 2020 10x Genomics, Inc. All rights reserved.

// This file contains functions that take as input an amino acid reference sequence for a V
// segment, along with its chain type (IGH, IGK, IGL, TRA or TRB), and attempt to find features in
// the sequence.
//
// These functions can fail, in particular for reference sequences that are damaged and possibly
// for other sequences.
//
// The same functions can also be applied to a protein sequence, however the sequence needs to be
// modified to add a fake leader sequence on the left (we use MXXXXXXXXXXXXXXXXXXXX), and to
// truncate on the right to trim a bit beyond the start of the CDR3.

// THIS IS A COPY FROM THE ENCLONE REPO TO FACILITATE EXPERIMENTATION!

use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// framework one region, which is immmediately after the signal peptide.  Chain type is one of
// IGH, IGK, IGL, TRA or TRB, and is not used at the moment.
//
// We find the best scoring location for this motif, starting in the first 50 amino acids
// (positions start at one).
//
// pos   weight  values
// 1     150     Q, E, D, G, K
// 2      50     V, I, Q, A
// 4     100     L, V, M
// 6     250     Q, E
// 22    250     C
// 23    250     C
//
// If the starting amino acid is C, we add one to the start position.

pub fn fr1_start(aa: &Vec<u8>, chain_type: &str) -> usize {
    // Define PWM.

    let mut pwm = Vec::<Vec<(usize, u8)>>::new();

    // #1

    pwm.push(vec![
        (150, b'Q'),
        (150, b'E'),
        (150, b'D'),
        (150, b'G'),
        (150, b'K'),
    ]);

    // #2

    pwm.push(vec![(50, b'V'), (50, b'I'), (50, b'Q'), (50, b'A')]);

    // #3

    pwm.push(vec![]);

    // #4

    pwm.push(vec![(100, b'L'), (100, b'V'), (100, b'M')]);

    // #5

    pwm.push(vec![]);

    // #6

    pwm.push(vec![(250, b'Q'), (250, b'E')]);

    for _ in 0..22 - 6 - 1 {
        pwm.push(vec![]);
    }

    // #22

    if chain_type == "IGH" {
        pwm.push(vec![(500, b'C')]);
    } else {
        pwm.push(vec![(250, b'C')]);
    }

    // #23

    pwm.push(vec![(250, b'C')]);

    // Score positions.

    let mut score_pos = Vec::<(usize, usize)>::new();
    for j in 0..=aa.len() - pwm.len() {
        if j > 40 || (chain_type == "IGL" && j > 25) {
            break;
        }
        let mut score = 0;
        for p in 0..pwm.len() {
            for l in 0..pwm[p].len() {
                if pwm[p][l].1 == aa[j + p] {
                    score += pwm[p][l].0;
                }
            }
        }
        score_pos.push((score, j));
    }
    reverse_sort(&mut score_pos);
    let mut p = score_pos[0].1;
    if aa[p] == b'C' {
        p += 1;
    }
    p
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// CDR1 region.  Note that there is more than one convention regarding the CDR1 start, and these
// conventions appear to differ by fixed offsets.  The convention used here is for IMGT.
// Chain type is one of IGH, IGK, IGL, TRA or TRB.

pub fn cdr1_start(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<usize> {
    // Define PWM for eight amino acids.

    let mut pwm = Vec::<Vec<(usize, u8)>>::new();

    let z = 19;

    // #19

    pwm.push(vec![(50, b'V')]);

    // #20

    pwm.push(vec![(30, b'T')]);

    // #21

    pwm.push(vec![(200, b'L'), (200, b'I'), (200, b'V'), (200, b'M')]);

    // #22

    pwm.push(vec![(80, b'S'), (80, b'T'), (80, b'R')]);

    // #23

    pwm.push(vec![(250, b'C')]);

    // #24

    pwm.push(vec![]);

    // #25

    pwm.push(vec![]);

    // #26

    pwm.push(vec![(100, b'S'), (100, b'I'), (100, b'D')]);

    // Score positions.

    let mut score_pos = Vec::<(usize, usize)>::new();
    let fr1 = fr1_start(&aa, &chain_type);
    for j in 0..=aa.len() - pwm.len() {
        if j + pwm.len() > fr1 + 27 {
            continue;
        }
        if j < 7 + z - 1 {
            continue;
        }
        let mut score = 0;
        for p in 0..pwm.len() {
            for l in 0..pwm[p].len() {
                if pwm[p][l].1 == aa[j + p] {
                    score += pwm[p][l].0;
                }
            }
        }
        if verbose {
            let seq = &aa[j..j + pwm.len()];
            println!("j = {}, seq = {}, score = {}", j, strme(&seq), score);
        }
        score_pos.push((score, j));
    }
    reverse_sort(&mut score_pos);
    let mut add = 3;
    if chain_type.starts_with("TR") || chain_type == "IGH" {
        add = 6;
    }
    if score_pos.is_empty() {
        return None;
    }
    Some(score_pos[0].1 + add + 2)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// FR2 region.  Chain type is one of IGH, IGK, IGL, TRA or TRB'

pub fn fr2_start(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<usize> {
    // Define PWM for six amino acids.

    let mut pwm = Vec::<Vec<(usize, u8)>>::new();

    let z = 39;

    // #39

    pwm.push(vec![(50, b'L'), (50, b'M'), (50, b'V'), (50, b'F')]);

    // #40

    if chain_type == "IGL" {
        pwm.push(vec![(250, b'Y')]);
    } else {
        pwm.push(vec![]);
    }

    // #41

    pwm.push(vec![(250, b'W')]);

    // #42

    pwm.push(vec![(150, b'Y')]);

    // #43

    pwm.push(vec![(100, b'R')]);

    // #44

    pwm.push(vec![(250, b'Q')]);

    // #45-46

    pwm.push(vec![]);
    pwm.push(vec![]);

    // #47

    pwm.push(vec![(110, b'G')]);

    // #48

    pwm.push(vec![(60, b'K'), (60, b'Q')]);

    // #49

    pwm.push(vec![(40, b'G'), (40, b'K'), (40, b'A')]);

    // Score positions.

    let mut score_pos = Vec::<(usize, usize)>::new();
    for j in 0..=aa.len() - pwm.len() {
        let start = 2;
        let mut stop = 35;
        if chain_type == "IGH" {
            stop = 24;
        }
        if j > stop + z - 1 || j < start + z - 1 {
            continue;
        }
        let mut score = 0;
        for p in 0..pwm.len() {
            for l in 0..pwm[p].len() {
                if pwm[p][l].1 == aa[j + p] {
                    score += pwm[p][l].0;
                }
            }
        }
        if verbose {
            println!("j = {}, score = {}", j, score);
        }
        score_pos.push((score, j));
    }
    reverse_sort(&mut score_pos);
    let mut add = 0;
    if chain_type != "IGH" {
        add = 1;
    }
    if chain_type == "IGK" || chain_type == "IGL" {
        add += 2;
    }
    if score_pos.is_empty() {
        None
    } else {
        Some(score_pos[0].1 + add - 1)
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// CDR2 region.  Chain type is one of IGH, IGK, IGL, TRA or TRB.

pub fn cdr2_start(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<usize> {
    let s2 = fr2_start(&aa, chain_type, false);
    if s2.is_none() {
        return None;
    }
    let s2 = s2.unwrap();
    let mut add = 0 as isize;
    if chain_type == "IGH" {
        // Six amino acids preceeding the CDR2 start.

        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // #50

        pwm.push(vec![(80, b'L')]);

        // #51

        pwm.push(vec![(80, b'E')]);

        // #52

        pwm.push(vec![(80, b'W')]);

        // #53

        pwm.push(vec![(40, b'V'), (40, b'M'), (40, b'I'), (40, b'L')]);

        // #54

        pwm.push(vec![(40, b'G'), (40, b'S'), (40, b'A')]);

        // #55

        pwm.push(vec![]);

        // Score positions.

        let mut score_pos = Vec::<(usize, usize)>::new();
        for j in 0..=aa.len() - pwm.len() {
            if j < s2 + 8 || j > s2 + 13 {
                continue;
            }
            let mut score = 0;
            for p in 0..pwm.len() {
                for l in 0..pwm[p].len() {
                    if pwm[p][l].1 == aa[j + p] {
                        score += pwm[p][l].0;
                    }
                }
            }
            if verbose {
                println!("j = {}, score = {}", j, score);
            }
            score_pos.push((score, j));
        }
        reverse_sort(&mut score_pos);
        Some(score_pos[0].1 + 7)
    } else if chain_type == "TRA" {
        // Six amino acids preceeding the CDR2 start.

        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // #50

        pwm.push(vec![(15, b'P'), (15, b'L')]);

        // #51

        pwm.push(vec![
            (15, b'Q'),
            (15, b'V'),
            (15, b'E'),
            (15, b'T'),
            (15, b'I'),
        ]);

        // #52

        pwm.push(vec![(20, b'L'), (20, b'F')]);

        // #53

        pwm.push(vec![(35, b'L')]);

        // #54

        pwm.push(vec![(15, b'L'), (15, b'I')]);

        // #55

        pwm.push(vec![]);

        // Score positions.

        let mut score_pos = Vec::<(usize, usize)>::new();
        for j in 0..=aa.len() - pwm.len() {
            if j < s2 + 10 || j > s2 + 12 {
                continue;
            }
            let mut score = 0;
            for p in 0..pwm.len() {
                for l in 0..pwm[p].len() {
                    if pwm[p][l].1 == aa[j + p] {
                        score += pwm[p][l].0;
                    }
                }
            }
            // if verbose {
            //     println!("j = {}, score = {}", j, score);
            // }
            score_pos.push((score, j));
        }
        reverse_sort(&mut score_pos);
        if score_pos.is_empty() {
            return None;
        } else {
            return Some(score_pos[0].1 + 6);
        }
    } else {
        if chain_type == "IGK" || chain_type == "IGL" {
            add = -2;
        }
        Some(s2 + (17 + add) as usize)
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn cdr3_start(aa: &Vec<u8>, _chain_type: &str, _verbose: bool) -> usize {
    let motif = [b"LQPEDSAVYYC", b"VEASQTGTYFC", b"ATSGQASLYLC"];
    let nm = motif[0].len();
    let reach = 18;
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

pub fn cdr3_score(aa: &Vec<u8>, _chain_type: &str, _verbose: bool) -> usize {
    let motif = [b"LQPEDSAVYYC", b"VEASQTGTYFC", b"ATSGQASLYLC"];
    let nm = motif[0].len();
    let reach = 18;
    let mut scores = Vec::<(usize, usize)>::new();
    for j in aa.len() - nm - reach..=aa.len() - nm {
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
    scores[0].0
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// FR3 region.

pub fn fr3_start(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<usize> {
    // First find the start of the CDR3.

    let cdr3_start = cdr3_start(&aa, &chain_type, verbose);
    /*
    use string_utils::*;
    println!(
        "bases before cdr3 = {}",
        strme(&aa[cdr3_start - 6..cdr3_start])
    );
    */

    // Do IGK and IGL.

    if chain_type == "IGK" || chain_type == "IGL" {
        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // x1

        pwm.push(vec![(100, b'G')]);

        // x2

        pwm.push(vec![]);

        // x3

        pwm.push(vec![(100, b'P')]);

        // x4

        pwm.push(vec![]);

        // x5

        pwm.push(vec![(100, b'R')]);

        // x6

        pwm.push(vec![(100, b'F')]);

        // x7

        pwm.push(vec![]);

        // x8

        pwm.push(vec![(100, b'G')]);

        // Score positions.

        let mut score_pos = Vec::<(usize, usize)>::new();
        for j in cdr3_start - 35..=cdr3_start - 28 {
            // changed to 39
            let mut score = 0;
            for p in 0..pwm.len() {
                for l in 0..pwm[p].len() {
                    if pwm[p][l].1 == aa[j + p] {
                        score += pwm[p][l].0;
                    }
                }
            }
            if verbose {
                println!("j = {}, score = {}", j, score);
            }
            score_pos.push((score, j));
        }
        reverse_sort(&mut score_pos);
        if !score_pos.is_empty() {
            return Some(score_pos[0].1);
        } else {
            return None;
        }

    // Do IGH.
    } else if chain_type == "IGH" {
        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // #2

        pwm.push(vec![(600, b'Y'), (600, b'N')]);

        // #3

        pwm.push(vec![(500, b'Y')]);

        // #4-6

        pwm.push(vec![(400, b'A'), (400, b'N')]);
        pwm.push(vec![]);
        pwm.push(vec![]);

        // #7

        pwm.push(vec![(850, b'F'), (850, b'L')]);

        // #8

        pwm.push(vec![(800, b'K'), (800, b'Q'), (800, b'R')]);

        // #9

        pwm.push(vec![]);

        // #10

        pwm.push(vec![(1000, b'R'), (1000, b'K')]);

        // #11

        pwm.push(vec![(700, b'F'), (700, b'V'), (700, b'A'), (700, b'L')]);

        // Score positions.

        let mut score_pos = Vec::<(usize, isize)>::new();
        for j in cdr3_start as isize - 42 + 2..=cdr3_start as isize - 32 + 2 {
            if j >= 0 {
                let mut score = 0;
                for p in 0..pwm.len() {
                    for l in 0..pwm[p].len() {
                        if pwm[p][l].1 == aa[j as usize + p] {
                            score += pwm[p][l].0;
                        }
                    }
                }
                // use string_utils::*;
                // println!("score of {} = {}", strme(&aa[j..j + pwm.len()]), score);
                if verbose {
                    println!("j = {}, score = {}", j, score);
                }
                score_pos.push((score, -(j as isize)));
            }
        }
        reverse_sort(&mut score_pos);
        if !score_pos.is_empty() {
            if -score_pos[0].1 >= 1 {
                return Some((-score_pos[0].1) as usize - 1);
            } else {
                return None;
            }
        } else {
            return None;
        }

    // Do TRA.
    } else if chain_type == "TRA" {
        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // #m0

        pwm.push(vec![(50, b'E'), (50, b'V'), (50, b'N'), (50, b'K')]);

        // #m1

        pwm.push(vec![(50, b'T'), (50, b'K'), (50, b'A'), (50, b'E')]);

        // #m2

        pwm.push(vec![(50, b'E'), (50, b'S')]);

        // #m3

        pwm.push(vec![(50, b'N'), (50, b'D'), (50, b'S')]);

        // #m4

        pwm.push(vec![(50, b'G'), (50, b'N')]);

        // #m5

        pwm.push(vec![(80, b'R'), (80, b'G'), (80, b'M')]);

        // #m6

        pwm.push(vec![(50, b'F'), (50, b'Y'), (50, b'A'), (50, b'I')]);

        // #m7

        pwm.push(vec![(50, b'S'), (50, b'T')]);

        // #m8

        pwm.push(vec![(50, b'A'), (50, b'V')]);

        // #m9

        pwm.push(vec![(50, b'T'), (50, b'E')]);

        // #m10

        pwm.push(vec![]);

        // #m11

        pwm.push(vec![(50, b'N'), (50, b'D')]);

        // Score positions.

        let mut score_pos = Vec::<(usize, usize)>::new();
        for j in cdr3_start - 36..=cdr3_start - 33 {
            let mut score = 0;
            for p in 0..pwm.len() {
                for l in 0..pwm[p].len() {
                    if pwm[p][l].1 == aa[j + p] {
                        score += pwm[p][l].0;
                    }
                }
            }
            // println!("score of {} = {}", strme(&aa[j..j + pwm.len()]), score);
            if verbose {
                println!("start = {}, score = {}", j + 1, score);
            }
            score_pos.push((score, j));
        }
        reverse_sort(&mut score_pos);
        if !score_pos.is_empty() {
            return Some(score_pos[0].1 + 1);
        } else {
            return None;
        }

    // Do TRB.
    } else {
        let mut pwm = Vec::<Vec<(usize, u8)>>::new();

        // #m1

        pwm.push(vec![]);

        // #m2

        pwm.push(vec![]);

        // #m3

        pwm.push(vec![(50, b'K'), (50, b'E'), (50, b'D')]);

        // #m4

        pwm.push(vec![(200, b'G'), (200, b'S'), (200, b'Q')]);

        // #m5

        pwm.push(vec![(200, b'D'), (200, b'E'), (200, b'G'), (200, b'S')]);

        // #m6

        pwm.push(vec![(200, b'I'), (200, b'V'), (200, b'L'), (200, b'M')]);

        // #m7

        pwm.push(vec![(100, b'P'), (100, b'S')]);

        // Score positions.

        let mut score_pos = Vec::<(usize, usize)>::new();
        for j in cdr3_start - 38..=cdr3_start - 35 {
            let mut score = 0;
            for p in 0..pwm.len() {
                for l in 0..pwm[p].len() {
                    if pwm[p][l].1 == aa[j + p] {
                        score += pwm[p][l].0;
                    }
                }
            }
            // println!("score of {} = {}", strme(&aa[j..j + pwm.len()]), score);
            if verbose {
                println!("j = {}, score = {}", j, score);
            }
            score_pos.push((score, j));
        }
        reverse_sort(&mut score_pos);
        if !score_pos.is_empty() {
            return Some(score_pos[0].1);
        } else {
            return None;
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, return the CDR1 sequence.
// Chain type is one of IGH, IGK, IGL, TRA or TRB, and is not used at the moment.

pub fn cdr1(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<Vec<u8>> {
    let fr2 = fr2_start(&aa, chain_type, verbose);
    if fr2.is_none() {
        return None;
    }
    Some(aa[cdr1_start(&aa, chain_type, verbose).unwrap()..fr2.unwrap()].to_vec())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, return the CDR2 sequence.
// Chain type is one of IGH, IGK, IGL, TRA or TRB.

pub fn cdr2(aa: &Vec<u8>, chain_type: &str, verbose: bool) -> Option<Vec<u8>> {
    let start = cdr2_start(&aa, chain_type, verbose);
    if start.is_none() {
        return None;
    }
    let start = start.unwrap();
    let stop = fr3_start(&aa, chain_type, false);
    if stop.is_none() {
        return None;
    }
    let stop = stop.unwrap();
    if start > stop {
        panic!(
            "Error in cdr2(...): cdr2_start = {} exceeds fr3_start = {}",
            start, stop
        );
    }
    Some(aa[start..stop].to_vec())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn score_fwr3(aa: &[u8], r: usize, freqs: &Vec<Vec<Vec<(u32, u8)>>>) -> f64 {
    let chain_type;
    if r == 0 {
        chain_type = "IGH";
    } else if r == 1 {
        chain_type = "IGK";
    } else if r == 2 {
        chain_type = "IGL";
    } else if r == 3 {
        chain_type = "TRA";
    } else {
        chain_type = "TRB";
    }
    let cdr3 = cdr3_start(&aa.to_vec(), chain_type, false);
    let motif = freqs[0].len();
    let mut score = 0.0;
    for j in 0..motif {
        let x = aa[cdr3 - j - 1];
        let mut m = 0;
        let mut total = 0;
        for k in 0..freqs[r][j].len() {
            let count = freqs[r][j][k].0;
            let y = freqs[r][j][k].1;
            total += count;
            if y == x {
                m += count;
            }
        }
        score += m as f64 / total as f64;
    }
    score
}

pub fn score4(aa: &[u8], r: usize) -> usize {
    let chain_type;
    if r == 0 {
        chain_type = "IGH";
    } else if r == 1 {
        chain_type = "IGK";
    } else if r == 2 {
        chain_type = "IGL";
    } else if r == 3 {
        chain_type = "TRA";
    } else {
        chain_type = "TRB";
    }
    let cdr3 = cdr3_start(&aa.to_vec(), chain_type, false);
    let n = aa.len();
    assert!(n >= 22);
    let mut score = 0;
    let x = aa[cdr3 - 4];
    if x == b'V' || x == b'T' || x == b'L' {
        score += 1;
    }
    let x = aa[cdr3 - 3];
    if x == b'Y' {
        score += 1;
    }
    let x = aa[cdr3 - 2];
    if x == b'Y' || x == b'F' || x == b'L' {
        score += 1;
    }
    let x = aa[cdr3 - 1];
    if x == b'C' {
        score += 3;
    }
    score
}
