// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Flag defective reference sequences.

use amino::*;
use enclone_core::defs::*;
use enclone_core::vdj_features::*;
use io_utils::*;
use itertools::Itertools;
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::RefData;
use vector_utils::*;

pub fn flag_defective(
    ctl: &EncloneControl,
    refdata: &RefData,
    log: &mut Vec<u8>,
    broken: &mut Vec<bool>,
) {
    *log = Vec::<u8>::new();
    let mut count = 0;
    *broken = vec![false; refdata.refs.len()];

    // Compute freqs.

    let mut freqs = Vec::<Vec<Vec<(u32, u8)>>>::new();
    let f = include_str!["../../enclone/src/fwr3_freqs.data"];
    let mut x = Vec::<(u32, u8)>::new();
    let mut y = Vec::<Vec<(u32, u8)>>::new();
    let mut rlast = 0;
    let mut ilast = 0;
    for line in f.lines() {
        let fields = line.split(',').collect::<Vec<&str>>();
        let r = fields[0].force_usize();
        let i = fields[1].force_usize();
        let count = fields[3].force_usize() as u32;
        let res = fields[4].as_bytes()[0];
        if r != rlast || i != ilast {
            y.push(x.clone());
            x.clear();
            if r != rlast {
                freqs.push(y.clone());
                y.clear();
            }
        }
        x.push((count, res));
        rlast = r;
        ilast = i;
    }
    y.push(x);
    freqs.push(y);

    // Study the reference.

    for i in 0..refdata.refs.len() {
        // Determine chain type and exclude those other than IGH, IGK, IGL, TRA and TRB.

        let rtype = refdata.rtype[i];
        let chain_type;
        if rtype == 0 {
            chain_type = "IGH";
        } else if rtype == 1 {
            chain_type = "IGK";
        } else if rtype == 2 {
            chain_type = "IGL";
        } else if rtype == 3 {
            chain_type = "TRA";
        } else if rtype == 4 {
            chain_type = "TRB";
        } else {
            continue;
        }

        // Look for problems.

        if refdata.is_c(i) {
            // This is very ugly.  We are exempting mouse IGHG2B because is is in our current
            // reference but has an extra base at the beginning.  See also comments below at
            // TRBV21-1.  Also, we're not actually checking for mouse.

            if refdata.name[i] == "IGHG2B" {
                continue;
            }

            // Continue.

            let seq = refdata.refs[i].to_ascii_vec();
            let aa0 = aa_seq(&seq, 0);
            let aa2 = aa_seq(&seq, 2);
            if aa2.contains(&b'*') && !aa0.contains(&b'*') {
                count += 1;
                broken[i] = true;
                fwriteln!(
                    log,
                    "{}. The following C segment reference sequence appears to have \
                    an extra base at its beginning:\n",
                    count
                );
                fwriteln!(log, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
            }
        } else if refdata.is_v(i) {
            // This is very ugly.  We are exempting human TRBV21-1 because it is in our current
            // reference (twice), but has multiple stop codons.  It should be deleted from the
            // reference, and then we should remove this test.  But probably in the future so as
            // not to inconvenience users.

            if ctl.gen_opt.species != "mouse" && refdata.name[i] == "TRBV21-1" {
                broken[i] = true;
                continue;
            }

            // This is very ugly.  We are exempting human IGHV1-12 because it is in our current
            // human reference, but is only 60 amino acids long.  Also we're not checking for mouse.

            if refdata.name[i] == "IGHV1-12" {
                broken[i] = true;
                continue;
            }

            // Ugly.  Frameshifted.  Also this is mouse and not checking for that.

            if refdata.name[i] == "TRAV23" {
                broken[i] = true;
                continue;
            }

            // Ugly.  Truncated on right.  Human.

            if refdata.name[i] == "IGLV5-48" {
                broken[i] = true;
                continue;
            }

            // Test for broken.

            let seq = refdata.refs[i].to_ascii_vec();
            let aa = aa_seq(&seq, 0);
            let mut reasons = Vec::<String>::new();
            if !aa.starts_with(b"M") {
                reasons.push("does not begin with a start codon".to_string());
            }
            let stops = aa.iter().filter(|&n| *n == b'*').count();
            if stops > 1 {
                reasons.push("has more than one stop codon".to_string());
            }
            if aa.len() < 100 {
                reasons.push("appears truncated (has less than 100 amino acids)".to_string());
            } else if aa.len() < 105 && chain_type == "IGH" {
                reasons.push("appears truncated (has less than 105 amino acids)".to_string());
            }
            if aa.len() >= 30 {
                let mut aap = aa.clone();
                aap.push(b'C');
                if cdr3_score(&aap, &chain_type, false) > 4 + cdr3_score(&aa, &chain_type, false) {
                    reasons.push("appears to need a C to be appended to its right end".to_string());
                }
            }
            if stops > 0 {
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
                    reasons.push("appears to be frameshifted".to_string());
                }
            }
            if aa.len() >= 31 {
                for del in 1..=2 {
                    let aad = aa_seq(&seq, del);
                    if cdr3_score(&aad, &chain_type, false)
                        > 4 + cdr3_score(&aa, &chain_type, false)
                    {
                        reasons.push("appears to be frameshifted".to_string());
                    }
                }
            }
            if reasons.is_empty() {
                let cs2 = cdr2_start(&aa, &chain_type, false);
                if cs2.is_some() {
                    let fr3 = fr3_start(&aa, &chain_type, false).unwrap();
                    if cs2.unwrap() > fr3 {
                        reasons.push(
                            "appears to be defective, because our computed \
                            CDR2 start exceeds our computed FWR3 start"
                                .to_string(),
                        );
                    }
                }
            }
            if reasons.is_empty() {
                if aa.len() >= 31 {
                    // Pretty crappy frameshift test.  One should see high aa and dna similarity
                    // to other seqs if shifted.  Or use more aas.
                    let score = cdr3_score(&aa, &chain_type, false);
                    let mut frameshift = false;
                    for del in 1..=2 {
                        let aad = aa_seq(&seq, del);
                        if score <= 6 && cdr3_score(&aad, &chain_type, false) >= 3 + score {
                            frameshift = true;
                        }
                    }
                    if frameshift {
                        reasons.push("appears to be frameshifted".to_string());
                    }
                }
            }
            if reasons.is_empty() {
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
                let score = score_fwr3(&aa, r, &freqs);
                if score < 8.0 && score4(&aa, r) < 5 {
                    reasons.push("appears to be frameshifted or truncated".to_string());
                }
            }

            // Report results.

            unique_sort(&mut reasons);
            if !reasons.is_empty() {
                let msg = format!(
                    "The following V segment reference sequence {}",
                    reasons.iter().format(", and ")
                );
                count += 1;
                broken[i] = true;
                fwriteln!(log, "{}. {}:\n", count, msg);
                fwriteln!(log, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq));
            }
        }
    }
}
