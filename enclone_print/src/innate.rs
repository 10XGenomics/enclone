// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Functions relating to the identification if iNKT and MAIT cells.

// species: return "human" or "mouse" or "unknown", based on a 60-base perfect match between
// the TRAC sequence in the provided reference sequences, and the internally provided reference
// sequences for human and mouse.

// use enclone_core::defs::*;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn species(refdata: &RefData) -> String {
    let mut my_trac = Vec::<Vec<u8>>::new();
    for i in 0..refdata.refs.len() {
        if refdata.name[i].starts_with("TRAC") {
            my_trac.push(refdata.refs[i].to_ascii_vec());
        }
    }
    const K: usize = 60;
    let mut kmers = Vec::<Vec<u8>>::new();
    for i in 0..my_trac.len() {
        for j in 0..=my_trac[i].len() - K {
            kmers.push(my_trac[i][j..j + K].to_vec());
        }
    }
    unique_sort(&mut kmers);
    let mut counts = Vec::<(usize, String)>::new();
    for pass in 1..=2 {
        let mut count = 0;
        let species;
        let refx;
        if pass == 1 {
            refx = human_ref();
            species = "human".to_string();
        } else {
            refx = mouse_ref();
            species = "mouse".to_string();
        }
        let mut trac = Vec::<u8>::new();
        let mut in_trac = false;
        for line in refx.lines() {
            if line.starts_with(">") && line.contains("|TRAC") {
                in_trac = true;
                continue;
            } else if line.starts_with(">") {
                break;
            }
            if in_trac {
                trac.append(&mut line.as_bytes().to_vec());
            }
        }
        for i in 0..=trac.len() - K {
            let kmer = trac[i..i + K].to_vec();
            if bin_member(&kmers, &kmer) {
                count += 1;
            }
        }
        counts.push((count, species));
    }
    reverse_sort(&mut counts);
    if counts[0] == counts[1] {
        return "unknown".to_string();
    } else {
        return counts[0].1.clone();
    }
}

// innate_cdr3: for the given species and given class (iNKT or MAIT), return the list of CDR3_AA
// sequences that are known to occur for that class.  These are defined in files in this directory.

pub fn innate_cdr3(species: &str, class: &str) -> Vec<String> {
    assert!(class == "iNKT" || class == "MAIT");
    let mut json = String::new();
    if species == "human" && class == "iNKT" {
        json = include_str!["human_iNKT_CDR3.json"].to_string();
    } else if species == "human" && class == "MAIT" {
        json = include_str!["human_MAIT_CDR3.json"].to_string();
    }
    let mut cdr3 = Vec::<String>::new();
    for line in json.lines() {
        if line.contains("\"cdr3\": ") {
            cdr3.push(line.after("\"cdr3\": ").between("\"", "\"").to_string());
        }
    }
    cdr3
}

// mark_innate: for each exact subclonotype, return "iNKT" or "MAIT" or "".

/*

pub fn mark_innate(refdata: &RefData, ex: &Vec<ExactClonotype>) -> Vec<String> {
    let species = species(&refdata);
    let inkt_cdr3 = innate_cdr3(&species, "iNKT");
    iet mait_cdr3 = innate_cdr3(&species, "MAIT");
    let innate = vec![String::new(); ex.len()];
    for i in 0..ex.len() {
        let (mut have_mait_tra, mut have_mait_trb) = (false, false);
        let (mut have_mait_tra_cdr3, mut have_mait_trb_cdr3) = (false, false);
        let mut vname = refdata.name[ex.share[j].v_ref_id].clone();
        if vname.contain('*') {
            vname = vname.before("*").to_string();
        }
        for j in 0..ex.share.len() {
            let mut vname = &refdata.name[ex.share[j].v_ref_id];
            if vname.contain('*') {
                vname = vname.before("*");
            }
            let mut jname = &refdata.name[ex.share[j].j_ref_id];
            if jname.contain('*') {
                jname = jname.before("*");
            }
            if species == "human".to_string() {
                if vname == "TRAV10" && jname == "TRAJ18" {
                    have_inkt_tra = true;
                } else if vname == "TRBV25-1" {
                    have_inkt_trb = true;
                } else if vname == "TRAV1-2"
                    && (jname == "TRAJ33" || jname == "TRAJ20" || jname == "TRAJ12") {
                    have_mait_tra = true;
                } else if vname == "TRBV20" || vname == TRBV6" {
                    have_mait_trb = true;
                }
            } else if species == "mouse".to_string() {
                if vname == "TRAV1" && jname == "TRAJ33" {
                    have_mait_tra = true;
                } else if vname == "TRBV19" || vname == "TRBV13" {
                    have_mait_trb = true;
                } else if vname == "TRAV11" && jname == "TRAJ18" {
                    have inkt_tra = true;
                } else if vname == "TRBV13-2" || vname == "TRBV1" || vname == "TRBV29" {
                    have_inkt_trb = true;
                }
            }
            if ex.share[j].left {
                if bin_member(&inkt_cdr3, &ex.share[j].cdr3_aa) {
                    have_inkt_trb_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &ex.share[j].cdr3_aa) {
                    have_mait_trb_cdr3 = true;
                }
            } else {
                if bin_member(&inkt_cdr3, &ex.share[j].cdr3_aa) {
                    have_inkt_tra_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &ex.share[j].cdr3_aa) {
                    have_mait_tra_cdr3 = true;
                }
            }
        }
        let have_inkt = have_inkt_tra && have_inkt_trb;
        let have_mait = have_mait_tra && have_mait_trb;

*/

/*

// iNKT and MAIT cell types
enum InvariantTCellType {
    INKT = 0; // invariant natural killer T cells
    MAIT = 1; // mucosal associated invariant T cells
}
// The amount of evidence we have about a single chain on whether the cell is iNKT/MAIT
enum InvariantTCellEvidence {
    GeneMatch = 0; // Matches the expected V/J genes
    JunctionMatch = 1; // Matches the V/J genes and the CDR3
}
// iNKT/MAIT annotation details
message InvariantTCellAnnotation {
    required InvariantTCellType cell_type = 1;
    // What level of evidence do we have in the alpha chain. None indicates either the chain is
    // missing or there is no match to the expected genes.
    optional InvariantTCellEvidence alpha_chain_evidence = 2;
    // What level of evidence do we have in the beta chain. None indicates either the chain is
    // missing or there is no match to the expected genes.
    optional InvariantTCellEvidence beta_chain_evidence = 3;
    // Total score for the evidence across both chains with a score of 1 for GeneMatch
    // and 2 for JunctionMatch
    required uint32 evidence_score = 4;
}

*/
