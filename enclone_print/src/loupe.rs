// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// This set of functions writes a protobuf data structure that
// Loupe uses to access clonotype data.

use enclone_proto::proto_io::write_proto;
use enclone_proto::PROTO_VERSION;
use vdj_ann::*;

use self::refx::*;
use amino::*;
use bio::alignment::pairwise::*;

use debruijn::dna_string::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use io_utils::*;
use vector_utils::*;

// Export donor reference/inferred alt allele sequences
pub fn make_donor_refs(
    alt_refs: &Vec<(usize, usize, DnaString)>,
    refdata: &RefData,
) -> Vec<DonorReferenceItem> {
    let mut drefs = Vec::<DonorReferenceItem>::new();
    let mut i = 0;
    while i < alt_refs.len() {
        let j = next_diff12_3(&alt_refs, i as i32) as usize;
        for k in i..j {
            let x = &alt_refs[k];
            let donor_id = x.0;
            let ref_id = x.1;
            let alt = x.2.to_ascii_vec();
            let alt_name = format!(
                "{}, donor {}, alt allele {}",
                refdata.name[ref_id].clone(),
                donor_id,
                k - i + 1
            );
            let refx = refdata.refs[ref_id].to_ascii_vec();
            let mut cigar = "".to_string();
            let mut matches = 0;
            for p in 0..refx.len() {
                if alt[p] == refx[p] {
                    matches += 1;
                } else {
                    if matches > 0 {
                        cigar.push_str(&format!("{}=", matches));
                    }
                    cigar.push_str("1X");
                    matches = 0;
                }
            }
            if matches > 0 {
                cigar.push_str(&format!("{}=", matches));
            }
            drefs.push(DonorReferenceItem {
                universal_idx: ref_id as u32,
                donor_idx: donor_id as u32,
                display_name: alt_name,
                region: Region::V.into(),
                nt_sequence: alt,
                universal_aln: Alignment {
                    ref_start: 0,
                    cigar: cigar,
                },
            });
        }
        i = j;
    }
    drefs
}

fn amino_acid(seq: &[u8], start: usize) -> Vec<u8> {
    seq[start..]
        .chunks_exact(3)
        .map(|codon| codon_to_aa(&codon))
        .collect()
}

pub fn make_loupe_clonotype(
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<usize>,
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
) -> Clonotype {
    // Define concatenated universal and donor reference sequences.

    let mat = &rsi.mat;
    let cols = rsi.uids.len();
    let mut concatu = vec![Vec::<u8>::new(); cols];
    let mut vstartu = vec![0; cols]; // Index of the start of the V region in concatu
    let mut concatd = vec![Vec::<u8>::new(); cols];
    let mut vstartd = vec![0; cols]; // Index of the start of the V region in concatd

    for cx in 0..cols {
        if rsi.uids[cx].is_some() {
            let mut x = refdata.refs[rsi.uids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append(&mut x);
            concatd[cx].append(&mut x);
        }
        let mut x = refdata.refs[rsi.vids[cx]].to_ascii_vec();
        vstartu[cx] = concatu[cx].len();
        vstartd[cx] = concatd[cx].len();
        concatu[cx].append(&mut x);
        if rsi.vpids[cx].is_none() {
            concatd[cx].append(&mut x);
        } else {
            let mut y = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
            concatd[cx].append(&mut y);
        }
        if rsi.dids[cx].is_some() {
            let mut x = refdata.refs[rsi.dids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append(&mut x);
            concatd[cx].append(&mut x);
        }
        let mut x = refdata.refs[rsi.jids[cx]].to_ascii_vec();
        concatu[cx].append(&mut x);
        concatd[cx].append(&mut x);
        if rsi.cids[cx].is_some() {
            let mut x = refdata.refs[rsi.cids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append(&mut x);
            concatd[cx].append(&mut x);
        }
    }

    // Define ClonotypeChains.

    let mut xchains = Vec::<ClonotypeChain>::new();
    for cx in 0..cols {
        let mut m0 = 0;
        let mut u0 = 0;
        for z in 0..mat[cx].len() {
            if mat[cx][z].is_some() {
                u0 = z;
                m0 = mat[cx][z].unwrap();
                break;
            }
        }
        let ex = &exact_clonotypes[exacts[u0]];
        let nt_sequence = ex.share[m0].full_seq.clone();
        let chain_type = ex.share[m0].chain_type.clone();

        let donor_v_idx = rsi.vpids[cx];
        let donor_j_idx = None;
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let universal_reference = concatu[cx].clone();
        let donor_reference = concatd[cx].clone();
        let mut aligner = Aligner::new(-6, -1, &score);
        let al = aligner.semiglobal(&nt_sequence, &universal_reference);
        let universal_reference_aln = Alignment::from(&al);
        let al = aligner.semiglobal(&nt_sequence, &donor_reference);
        let donor_reference_aln = Alignment::from(&al);

        let aa_sequence = amino_acid(&nt_sequence, ex.share[m0].v_start);
        let aa_sequence_universal = amino_acid(&universal_reference, vstartu[cx]);
        let aa_sequence_donor = amino_acid(&donor_reference, vstartd[cx]);

        let v_start = ex.share[m0].v_start;
        let fwr1_start = Some((v_start + ex.share[m0].fr1_start) as u32);
        let mut cdr1_start = None;
        if ex.share[m0].cdr1_start.is_some() {
            cdr1_start = Some((v_start + ex.share[m0].cdr1_start.unwrap()) as u32);
        }
        let mut fwr2_start = None;
        if ex.share[m0].fr2_start.is_some() {
            fwr2_start = Some((v_start + ex.share[m0].fr2_start.unwrap()) as u32);
        }
        let mut cdr2_start = None;
        if ex.share[m0].cdr2_start.is_some() {
            cdr2_start = Some((v_start + ex.share[m0].cdr2_start.unwrap()) as u32);
        }
        let mut fwr3_start = None;
        if ex.share[m0].fr3_start.is_some() {
            fwr3_start = Some((v_start + ex.share[m0].fr3_start.unwrap()) as u32);
        }
        let fwr4_end = Some(ex.share[m0].full_seq.len() as u32);

        xchains.push(ClonotypeChain {
            nt_sequence: nt_sequence,
            aa_sequence: aa_sequence,
            u_idx: rsi.uids[cx].map(|idx| idx as u32),
            v_idx: rsi.vids[cx] as u32,
            d_idx: rsi.dids[cx].map(|idx| idx as u32),
            j_idx: rsi.jids[cx] as u32,
            c_idx: rsi.cids[cx].map(|idx| idx as u32),
            donor_v_idx: donor_v_idx.map(|idx| idx as u32),
            donor_j_idx: donor_j_idx,
            universal_reference: universal_reference,
            universal_reference_aln: universal_reference_aln,
            aa_sequence_universal,
            donor_reference: donor_reference,
            donor_reference_aln: donor_reference_aln,
            aa_sequence_donor,
            v_start: ex.share[m0].v_start as u32,
            v_end: ex.share[m0].v_stop as u32,
            v_end_ref: ex.share[m0].v_stop_ref as u32,
            j_start: ex.share[m0].j_start as u32,
            j_start_ref: ex.share[m0].j_start_ref as u32,
            j_end: ex.share[m0].j_stop as u32,
            fwr1_start: fwr1_start,
            cdr1_start: cdr1_start,
            fwr2_start: fwr2_start,
            cdr2_start: cdr2_start,
            fwr3_start: fwr3_start,
            cdr3_start: ex.share[m0].v_start as u32 + ex.share[m0].cdr3_start as u32,
            cdr3_end: ex.share[m0].v_start as u32
                + (ex.share[m0].cdr3_start + 3 * ex.share[m0].cdr3_aa.len()) as u32,
            fwr4_end: fwr4_end,
            chain_type: chain_type,
        });
    }

    // Define an EClonotype.

    let mut ecl = Vec::<ExactSubClonotype>::new();
    for j in 0..exacts.len() {
        let mut chains = Vec::<Option<ExactSubClonotypeChain>>::new();
        let ex = &exact_clonotypes[exacts[j]];
        for cx in 0..cols {
            let m = mat[cx][j];
            if !m.is_some() {
                chains.push(None);
                continue;
            }
            let m = m.unwrap();
            let nt_sequence = ex.share[m].full_seq.clone();
            let v_start = ex.share[m].v_start;
            let mut aa_sequence = Vec::<u8>::new();
            let mut p = v_start;
            while p + 3 <= nt_sequence.len() {
                aa_sequence.push(codon_to_aa(&nt_sequence[p..p + 3]));
                p += 3;
            }
            let j_end = ex.share[m].j_stop;
            let c_region_idx = rsi.cids[cx];
            let fwr1_start = Some((v_start + ex.share[m].fr1_start) as u32);
            let mut cdr1_start = None;
            if ex.share[m].cdr1_start.is_some() {
                cdr1_start = Some((v_start + ex.share[m].cdr1_start.unwrap()) as u32);
            }
            let mut fwr2_start = None;
            if ex.share[m].fr2_start.is_some() {
                fwr2_start = Some((v_start + ex.share[m].fr2_start.unwrap()) as u32);
            }
            let mut cdr2_start = None;
            if ex.share[m].cdr2_start.is_some() {
                cdr2_start = Some((v_start + ex.share[m].cdr2_start.unwrap()) as u32);
            }
            let mut fwr3_start = None;
            if ex.share[m].fr3_start.is_some() {
                fwr3_start = Some((v_start + ex.share[m].fr3_start.unwrap()) as u32);
            }
            let cdr3_start = v_start + ex.share[m].cdr3_start;
            let cdr3_end = cdr3_start + ex.share[m].cdr3_dna.len();
            let fwr4_end = Some(ex.share[m].full_seq.len() as u32);
            let mut umi_counts = Vec::<u32>::new();
            let mut read_counts = Vec::<u32>::new();
            let mut contig_ids = Vec::<String>::new();
            for l in 0..ex.clones.len() {
                umi_counts.push(ex.clones[l][m].umi_count as u32);
                read_counts.push(ex.clones[l][m].read_count as u32);
                contig_ids.push(ex.clones[l][m].tigname.clone());
            }

            // Build alignments.  Don't know that these are sensible
            // alignment parameters.

            let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
            let mut aligner = Aligner::new(-6, -1, &score);
            let al = aligner.semiglobal(&nt_sequence, &xchains[cx].nt_sequence);
            let clonotype_consensus_aln = Alignment::from(&al);
            let al = aligner.semiglobal(&nt_sequence, &concatu[cx]);
            let universal_reference_aln = Alignment::from(&al);
            let al = aligner.semiglobal(&nt_sequence, &concatd[cx]);
            let donor_reference_aln = Alignment::from(&al);

            // Finally, define the ExactClonotypeChain.

            chains.push(Some(ExactSubClonotypeChain {
                nt_sequence: nt_sequence,
                aa_sequence: aa_sequence,
                v_start: v_start as u32,
                j_end: j_end as u32,
                c_region_idx: c_region_idx.map(|idx| idx as u32),
                fwr1_start: fwr1_start,
                cdr1_start: cdr1_start,
                fwr2_start: fwr2_start,
                cdr2_start: cdr2_start,
                fwr3_start: fwr3_start,
                cdr3_start: cdr3_start as u32,
                cdr3_end: cdr3_end as u32,
                fwr4_end: fwr4_end,
                umi_counts: umi_counts,
                read_counts: read_counts,
                contig_ids: contig_ids,
                clonotype_consensus_aln: clonotype_consensus_aln,
                donor_reference_aln: donor_reference_aln,
                universal_reference_aln: universal_reference_aln,
            }));
        }
        let mut cell_barcodes = Vec::<String>::new();
        for l in 0..ex.clones.len() {
            cell_barcodes.push(ex.clones[l][0].barcode.clone());
        }
        let inkt_evidence = InvariantTCellAnnotation {
            alpha_chain_gene_match: ex.share[0].inkt_alpha_chain_gene_match,
            alpha_chain_junction_match: ex.share[0].inkt_alpha_chain_junction_match,
            beta_chain_gene_match: ex.share[0].inkt_beta_chain_gene_match,
            beta_chain_junction_match: ex.share[0].inkt_beta_chain_junction_match,
        };
        let mait_evidence = InvariantTCellAnnotation {
            alpha_chain_gene_match: ex.share[0].mait_alpha_chain_gene_match,
            alpha_chain_junction_match: ex.share[0].mait_alpha_chain_junction_match,
            beta_chain_gene_match: ex.share[0].mait_beta_chain_gene_match,
            beta_chain_junction_match: ex.share[0].mait_beta_chain_junction_match,
        };
        ecl.push(ExactSubClonotype {
            chains: chains
                .into_iter()
                .enumerate()
                .filter_map(|(index, opt_chain)| {
                    opt_chain.map(|chain| ExactSubClonotypeChainInfo {
                        index: index as u32,
                        chain,
                    })
                })
                .collect(),
            cell_barcodes: cell_barcodes,
            inkt_evidence: inkt_evidence,
            mait_evidence: mait_evidence,
        });
    }

    // Build Clonotype.

    let mut n = 0;
    for i in 0..ecl.len() {
        n += ecl[i].cell_barcodes.len();
    }
    Clonotype {
        chains: xchains,
        exact_clonotypes: ecl,
        frequency: n as u32,
    }
}

pub fn loupe_out(
    ctl: &EncloneControl,
    all_loupe_clonotypes: Vec<Clonotype>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
) {
    if ctl.gen_opt.binary.len() > 0 || ctl.gen_opt.proto.len() > 0 {
        let mut uref = Vec::new();
        for i in 0..refdata.refs.len() {
            uref.push(UniversalReferenceItem {
                ref_idx: refdata.id[i] as u32,
                display_name: refdata.name[i].clone(),
                region: match refdata.segtype[i].as_str() {
                    "U" => Region::U.into(),
                    "V" => Region::V.into(),
                    "D" => Region::D.into(),
                    "J" => Region::J.into(),
                    "C" => Region::C.into(),
                    _ => unreachable!(),
                },
                nt_sequence: refdata.refs[i].to_ascii_vec(),
            });
        }
        let metadata = match &ctl.gen_opt.proto_metadata {
            Some(fname) => serde_json::from_reader(
                std::fs::File::open(fname).expect(&format!("Error while reading {}", fname)),
            )
            .expect(&format!("Unable to deserialize Metadata from {}", fname)),
            None => Metadata::default(),
        };
        let enclone_outputs = EncloneOutputs {
            version: PROTO_VERSION.into(),
            metadata,
            num_clonotypes: all_loupe_clonotypes.len() as u32,
            clonotypes: all_loupe_clonotypes,
            universal_reference: UniversalReference { items: uref },
            donor_reference: DonorReference {
                items: dref.to_vec(),
            },
        };
        if ctl.gen_opt.binary.len() > 0 {
            write_obj(&enclone_outputs, &ctl.gen_opt.binary);
        }
        if ctl.gen_opt.proto.len() > 0 {
            write_proto(enclone_outputs, &ctl.gen_opt.proto).unwrap();
        }
    }
}
