// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

use vdj_ann::*;

use amino::*;
use bio::alignment::AlignmentOperation::*;
use bio::alignment::pairwise::*;
use debruijn::dna_string::*;
use defs::*;
use io_utils::*;
use self::refx::*;
use types::*;
use vector_utils::*;

pub fn make_donor_refs( alt_refs: &Vec<(usize,usize,DnaString)>, refdata: &RefData ) 
    -> Vec<DonorReferenceItem> {
    let mut drefs = Vec::<DonorReferenceItem>::new();
    // if ctl.gen_opt.loupe.len() > 0 {
        let mut i = 0;
        while i < alt_refs.len() {
            let j = next_diff12_3(&alt_refs, i as i32) as usize;
            for k in i..j {
                let x = &alt_refs[k];
                let donor_id = x.0;
                let ref_id = x.1;
                let alt = x.2.to_ascii_vec();
                let alt_name = format!( "{}, donor {}, alt allele {}", 
                    refdata.name[ref_id].clone(), donor_id, k-i+1 );
                let refx = refdata.refs[ref_id].to_ascii_vec();
                let mut cigar = Vec::<CigarOps>::new();
                let mut matches = 0;
                for p in 0..refx.len() {
                    if alt[p] == refx[p] {
                        matches += 1;
                    } else {
                        if matches > 0 {
                            cigar.push( CigarOps::Match(matches as u16) );
                        }
                        cigar.push( CigarOps::Mismatch(1 as u16) );
                        matches = 0;
                    }
                }
                if matches > 0 {
                    cigar.push( CigarOps::Match(matches as u16) );
                }
                drefs.push( DonorReferenceItem {
                    universal_idx: ref_id,
                    donor_idx: donor_id,
                    display_name: alt_name,
                    region: Region::V,
                    nt_sequence: alt,
                    universal_aln: Alignment { ref_start: 0, cigar: cigar },
                } );
            }
            i = j;
        }
    // }
    drefs
}

pub fn convert_alignment_bio_to_types( al : &bio::alignment::Alignment ) -> Alignment {
    let mut cig = Vec::<CigarOps>::new();
    let mut i = 0;
    while i < al.operations.len() {
        let j = next_diff(&al.operations, i);
        let mult = (j-i) as u16;
        match al.operations[i] {
            Match => { cig.push( CigarOps::Match(mult) ); },
            Subst => { cig.push( CigarOps::Mismatch(mult) ); },
            Del => { cig.push( CigarOps::Deletion(mult) ); },
            Ins => { cig.push( CigarOps::Insertion(mult) ); },
            Xclip(n) => {
                for _ in i..j {
                    cig.push( CigarOps::SoftClip(n as u16) );
                }
            }
            // THIS SEEMS VERY DANGEROUS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            _ => { panic!(); }
        };
        i = j;
    }
    Alignment {
        ref_start: al.ystart,
        cigar: cig,
    }
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
    let mut concatu = vec![ Vec::<u8>::new(); cols ];
    let mut concatd = vec![ Vec::<u8>::new(); cols ];
    for cx in 0..cols {
        if rsi.uids[cx].is_some() {
            let mut x = refdata.refs[rsi.uids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append( &mut x );
            concatd[cx].append( &mut x );
        }
        let mut x = refdata.refs[rsi.vids[cx]].to_ascii_vec();
        concatu[cx].append( &mut x );
        if rsi.vpids[cx].is_none() {
            concatd[cx].append( &mut x );
        } else {
            let mut y = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
            concatd[cx].append( &mut y );
        }
        if rsi.dids[cx].is_some() {
            let mut x = refdata.refs[rsi.dids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append( &mut x );
            concatd[cx].append( &mut x );
        }
        let mut x = refdata.refs[rsi.jids[cx]].to_ascii_vec();
        concatu[cx].append( &mut x );
        concatd[cx].append( &mut x );
        if rsi.cids[cx].is_some() {
            let mut x = refdata.refs[rsi.cids[cx].unwrap()].to_ascii_vec();
            concatu[cx].append( &mut x );
            concatd[cx].append( &mut x );
        }
    }

    // Define ClonotypeChains.

    let mut xchains = Vec::<ClonotypeChain>::new();
    for cx in 0..cols {

        // Find the first exact subclonotype that has an entry in the given column.  Over
        // 99% of the time, this will be the first exact subclonotype.

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

        // Compute the ClonotypeChain from this exact subclonotype.

        let nt_sequence = ex.share[m0].full_seq.clone();
        let u_idx = rsi.uids[cx];
        let v_idx = rsi.vids[cx];
        let d_idx = rsi.dids[cx];
        let j_idx = rsi.jids[cx];
        let donor_v_idx = rsi.vpids[cx];
        let donor_j_idx = None;
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = Aligner::new(-6, -1, &score);
        let mut donor_reference_aln = None;
        if rsi.vpids[cx].is_some() {
            let al = aligner.semiglobal(&nt_sequence, &concatd[cx]);
            donor_reference_aln = Some(convert_alignment_bio_to_types(&al));
        }
        let al = aligner.semiglobal(&nt_sequence, &concatu[cx]);
        let universal_reference_aln = convert_alignment_bio_to_types(&al);
        xchains.push( ClonotypeChain {
            nt_sequence: nt_sequence,
            u_idx: u_idx,
            v_idx: v_idx,
            d_idx: d_idx,
            j_idx: j_idx,
            donor_v_idx: donor_v_idx,
            donor_j_idx: donor_j_idx,
            donor_reference_aln: donor_reference_aln,
            universal_reference_aln: universal_reference_aln,
        } );
    }

    // Define an EClonotype.

    let mut ecl = Vec::<EClonotype>::new();
    for j in 0..exacts.len() {
        let mut chains = Vec::<Option<ExactClonotypeChain>>::new();
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
                aa_sequence.push( codon_to_aa(&nt_sequence[p.. p + 3]) );
                p += 3;
            }
            let j_end = ex.share[m].j_stop;
            let c_region_idx = rsi.cids[cx];
            let cdr3_start = v_start + ex.share[m].cdr3_start;
            let cdr3_end = cdr3_start + ex.share[m].cdr3_dna.len();
            let mut umi_counts =  Vec::<u32>::new();
            let mut contig_ids = Vec::<String>::new();
            if mat[cx][j].is_some() {
                let m = mat[cx][j].unwrap();
                for l in 0..ex.clones.len() {
                    umi_counts.push( ex.clones[l][m].umi_count as u32 );
                    contig_ids.push( ex.clones[l][m].tigname.clone() );
                }
            }

            // Build alignments.  Don't know that these are sensible
            // alignment parameters.

            let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
            let mut aligner = Aligner::new(-6, -1, &score);
            let al = aligner.semiglobal(&ex.share[m].seq, 
                &xchains[cx].nt_sequence);
            let clonotype_consensus_aln = convert_alignment_bio_to_types(&al);
            let al = aligner.semiglobal(&ex.share[m].seq, &concatu[cx]);
            let universal_reference_aln = convert_alignment_bio_to_types(&al);
            let al = aligner.semiglobal(&ex.share[m].seq, &concatd[cx]);
            let donor_reference_aln;
            if rsi.vpids[cx].is_none() {
                donor_reference_aln = None;
            } else {
                donor_reference_aln = Some(convert_alignment_bio_to_types(&al));
            }

            // Finally, define the ExactClonotypeChain.

            chains.push( Some( ExactClonotypeChain {
                nt_sequence: nt_sequence,
                aa_sequence: aa_sequence,
                v_start: v_start,
                j_end: j_end,
                c_region_idx: c_region_idx,
                cdr3_start: cdr3_start,
                cdr3_end: cdr3_end,
                umi_counts: umi_counts,
                contig_ids: contig_ids,
                clonotype_consensus_aln: clonotype_consensus_aln,
                donor_reference_aln: donor_reference_aln,
                universal_reference_aln: universal_reference_aln,
            } ) );
        }
        let mut cell_barcodes = Vec::<String>::new();
        for l in 0..ex.clones.len() {
            cell_barcodes.push(ex.clones[l][0].barcode.clone());
        }
        cell_barcodes.sort(); // not sure this makes sense
        ecl.push( EClonotype{ chains: chains, cell_barcodes: cell_barcodes } );
    }

    // Build Clonotype.

    let mut n = 0;
    for i in 0..ecl.len() {
        n += ecl[i].cell_barcodes.len();
    }
    Clonotype { chains: xchains, exact_clonotypes: ecl, frequency: n as u32 }
}

pub fn loupe_out( ctl: &EncloneControl, all_loupe_clonotypes: Vec<Clonotype>,
    refdata: &RefData, dref: &Vec<DonorReferenceItem> ) {
    if ctl.gen_opt.loupe.len() > 0 {
        let mut uref = UniversalReference::new();
        for i in 0..refdata.refs.len() {
            uref.push( UniversalReferenceItem {
                ref_idx: refdata.id[i] as u32,
                display_name: refdata.name[i].clone(),
                region: match refdata.segtype[i].as_str() {
                    "U" => Region::U,
                    "V" => Region::V,
                    "D" => Region::D,
                    "J" => Region::J,
                    "C" => Region::C,
                    _ => Region::U, // nonsense!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                },
                nt_sequence: refdata.refs[i].to_ascii_vec(),
            } );
        }
        let enclone_outputs = EncloneOutputs {
            clonotypes: all_loupe_clonotypes,
            universal_reference: uref,
            donor_reference: dref.to_vec(),
        };
        write_obj(&enclone_outputs, &ctl.gen_opt.loupe);
    }
}
