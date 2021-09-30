// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Functions relating to the identification if iNKT and MAIT cells.

// species: return "human" or "mouse" or "unknown", based on a 60-base perfect match between
// the TRAC sequence in the provided reference sequences, and the internally provided reference
// sequences for human and mouse.

use enclone_core::defs::ExactClonotype;
use string_utils::TextUtils;
use vdj_ann::refx::{human_ref, mouse_ref, RefData};
use vector_utils::{bin_member, reverse_sort, unique_sort};

pub fn species(refdata: &RefData) -> String {
    let mut my_trac = Vec::<Vec<u8>>::new();
    for i in 0..refdata.refs.len() {
        if refdata.name[i].starts_with("TRAC") || refdata.name[i].starts_with("IGHM") {
            my_trac.push(refdata.refs[i].to_ascii_vec());
        }
    }
    const K: usize = 60;
    let mut kmers = Vec::<Vec<u8>>::new();
    for i in 0..my_trac.len() {
        if my_trac[i].len() >= K {
            for j in 0..=my_trac[i].len() - K {
                kmers.push(my_trac[i][j..j + K].to_vec());
            }
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
            if line.starts_with('>') && (line.contains("|TRAC") || line.contains("|IGHM")) {
                in_trac = true;
                continue;
            }
            if in_trac {
                trac.append(&mut line.as_bytes().to_vec());
                trac.push(b' ');
                in_trac = false;
            }
        }
        if trac.len() < K {
            return "unknown".to_string();
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
    if counts[0].0 == counts[1].0 {
        "unknown".to_string()
    } else {
        counts[0].1.clone()
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
    unique_sort(&mut cdr3);
    cdr3
}

// mark_innate: for each exact subclonotype, fill in iNKT and MAIT fields.

pub fn mark_innate(refdata: &RefData, ex: &mut Vec<ExactClonotype>) {
    let species = species(refdata);
    let inkt_cdr3 = innate_cdr3(&species, "iNKT");
    let mait_cdr3 = innate_cdr3(&species, "MAIT");
    for i in 0..ex.len() {
        let (mut have_mait_tra, mut have_mait_trb) = (false, false);
        let (mut have_mait_tra_cdr3, mut have_mait_trb_cdr3) = (false, false);
        let (mut have_inkt_tra, mut have_inkt_trb) = (false, false);
        let (mut have_inkt_tra_cdr3, mut have_inkt_trb_cdr3) = (false, false);
        for j in 0..ex[i].share.len() {
            let mut vname = refdata.name[ex[i].share[j].v_ref_id].clone();
            if vname.contains('*') {
                vname = vname.before("*").to_string();
            }
            let mut jname = refdata.name[ex[i].share[j].j_ref_id].clone();
            if jname.contains('*') {
                jname = jname.before("*").to_string();
            }
            if species == *"human" {
                if vname == "TRAV10" && jname == "TRAJ18" {
                    have_inkt_tra = true;
                }
                if vname == "TRBV25-1" {
                    have_inkt_trb = true;
                }
                if vname == "TRAV1-2"
                    && (jname == "TRAJ33" || jname == "TRAJ20" || jname == "TRAJ12")
                {
                    have_mait_tra = true;
                }
                if vname.starts_with("TRBV20") || vname.starts_with("TRBV6") {
                    have_mait_trb = true;
                }
            } else if species == *"mouse" {
                if vname == "TRAV1" && jname == "TRAJ33" {
                    have_mait_tra = true;
                }
                if vname == "TRBV19"
                    || vname == "TRBV13-1"
                    || vname == "TRBV13-2"
                    || vname == "TRBV13-3"
                {
                    have_mait_trb = true;
                }
                if (vname == "TRAV11" || vname == "TRAV11D") && jname == "TRAJ18" {
                    have_inkt_tra = true;
                }
                if vname == "TRBV13-2" || vname == "TRBV1" || vname == "TRBV29" {
                    have_inkt_trb = true;
                }
            }
            if ex[i].share[j].left {
                if bin_member(&inkt_cdr3, &ex[i].share[j].cdr3_aa) {
                    have_inkt_trb_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &ex[i].share[j].cdr3_aa) {
                    have_mait_trb_cdr3 = true;
                }
            } else {
                if bin_member(&inkt_cdr3, &ex[i].share[j].cdr3_aa) {
                    have_inkt_tra_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &ex[i].share[j].cdr3_aa) {
                    have_mait_tra_cdr3 = true;
                }
            }
        }
        for j in 0..ex[i].share.len() {
            ex[i].share[j].inkt_alpha_chain_gene_match = have_inkt_tra;
            ex[i].share[j].inkt_alpha_chain_junction_match = have_inkt_tra_cdr3;
            ex[i].share[j].inkt_beta_chain_gene_match = have_inkt_trb;
            ex[i].share[j].inkt_beta_chain_junction_match = have_inkt_trb_cdr3;
            ex[i].share[j].mait_alpha_chain_gene_match = have_mait_tra;
            ex[i].share[j].mait_alpha_chain_junction_match = have_mait_tra_cdr3;
            ex[i].share[j].mait_beta_chain_gene_match = have_mait_trb;
            ex[i].share[j].mait_beta_chain_junction_match = have_mait_trb_cdr3;
        }
    }
}
