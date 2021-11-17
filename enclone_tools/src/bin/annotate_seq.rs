// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

// annotate_seq species sequence [option_1 ... option_n]
//
// where species is human or mouse, sequence is a DNA sequence or a fasta filename
// and the following options are allowed:
// RC     reverse complement the sequence
// IMGT   use IMGT reference
// V      use verbose logging
// TCR    use only TCR reference
// BCR    use only BCR reference
// EXT    use both the plain and extended references, and do both fw and rc
// D      for IGH and TRB, look for the best aligning D segments and report them.
//        (Mostly this is confusing.  It is very hard to say why one D segment alignment is
//         better than another.)
// REF=fasta_file_name = use this reference
//
// Annotate the given DNA sequences.

extern crate debruijn;
extern crate fasta_tools;
#[macro_use]
extern crate io_utils;
extern crate pretty_trace;
extern crate string_utils;

extern crate vdj_ann;
use vdj_ann::*;

use annotate::*;
use bio_edit::alignment::pairwise::*;
use bio_edit::alignment::AlignmentOperation::*;
use debruijn::{dna_string::*, *};
use fasta_tools::*;
use itertools::Itertools;
use pretty_trace::*;
use refx::*;
use std::{cmp::min, env, io::Write};
use string_utils::*;
use vdj_ann::transcript::*;
use vdj_ann_ref::make_vdj_ref_data;
use vector_utils::*;

fn main() {
    // Set up and parse args.

    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let species = args[1].to_string();
    let mut is_dna = true;
    for c in args[2].clone().into_bytes() {
        if c != b'A' && c != b'C' && c != b'G' && c != b'T' {
            is_dna = false;
        }
    }
    let mut seqs = Vec::<DnaString>::new();
    let mut headers = Vec::<String>::new();
    if is_dna {
        seqs.push(DnaString::from_dna_string(&args[2]));
        headers.push(String::new());
    } else {
        read_fasta_into_vec_dna_string_plus_headers(&args[2], &mut seqs, &mut headers);
    }
    let mut imgt = false;
    let mut verbose = false;
    let (mut is_tcr, mut is_bcr) = (true, true);
    let mut ext = false;
    let mut d = false;
    let mut rc = false;
    let mut ref_fasta = String::new();
    for j in 3..args.len() {
        if args[j] == "IMGT" {
            imgt = true;
        }
        if args[j] == "V" {
            verbose = true;
        }
        if args[j] == "TCR" {
            is_bcr = false;
        }
        if args[j] == "BCR" {
            is_tcr = false;
        }
        if args[j] == "EXT" {
            ext = true;
        }
        if args[j] == "D" {
            d = true;
        }
        if args[j] == "RC" {
            rc = true;
        }
        if args[j].starts_with("REF=") {
            ref_fasta = args[j].after("REF=").to_string();
        }
    }

    // Make reference data.

    let mut refdata = RefData::new();
    if ref_fasta.len() > 0 {
        let refx = std::fs::read_to_string(&ref_fasta).unwrap();
        let ext_ref = String::new();
        make_vdj_ref_data_core(&mut refdata, &refx, &ext_ref, true, true, None);
    } else {
        make_vdj_ref_data(&mut refdata, imgt, &species, false, is_tcr, is_bcr);
    }

    // Traverse the sequences.

    println!("");
    for i in 0..seqs.len() {
        if !is_dna {
            println!("{}\n", headers[i]);
        }
        let mut seq = seqs[i].clone();
        if rc {
            seq = seq.rc();
        }
        let mut log = Vec::<u8>::new();

        // Annotate using plain reference.

        if ext {
            fwriteln!(log, "\nFW ANNOTATION VERSUS PLAIN REFERENCE\n");
        }
        print_annotations(&seq, &refdata, &mut log, false, true, verbose);
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&seq, &refdata, &mut ann, true, false, true);
        print_cdr3_using_ann(&seq, &refdata, &ann, &mut log);
        print_start_codon_positions(&seq, &mut log);
        if is_valid(&seq, &refdata, &ann, true, &mut log) {
            fwriteln!(log, "VALID");
        }

        // Handle d.  The code here mirrors that in the function row_fill (for "comp") in enclone.
        //
        // Except that -- and this is a cause for concern -- we use the end of the J rather than
        // the end of the CDR3 as the alignment stop point.
        //
        // Also separated edits by centered dots.

        if d {
            let (mut v, mut j) = (None, None);
            let (mut vstart, mut jstop) = (0, 0);
            for i in 0..ann.len() {
                let t = ann[i].2 as usize;
                if refdata.rtype[t] == 0 as i32 || refdata.rtype[t] == 4 as i32 {
                    if refdata.is_v(t) {
                        v = Some(t);
                    } else if refdata.is_j(t) {
                        j = Some(t);
                        jstop = (ann[i].0 + ann[i].1) as usize;
                    }
                }
            }
            for i in 0..ann.len() {
                let t = ann[i].2 as usize;
                if refdata.is_v(t) {
                    vstart = ann[i].0 as usize;
                    break;
                }
            }
            let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
            get_cdr3_using_ann(&seq, &refdata, &ann, &mut cdr3);
            if v.is_some() && j.is_some() && cdr3.len() == 1 {
                let tig = seq.to_ascii_vec();
                let tig = &tig[vstart..jstop];
                let (v, j) = (v.unwrap(), j.unwrap());
                let mut results = Vec::<(usize, usize, String)>::new();
                let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                let mut aligner = Aligner::new(-6, -1, &score);
                let start = cdr3[0].0 - vstart;
                let stop = jstop - vstart;
                for z in 0..refdata.ds.len() {
                    let d = refdata.ds[z];
                    if refdata.name[v][0..2] == refdata.name[d][0..2] {
                        let mut concat = refdata.refs[v].to_ascii_vec();
                        concat.append(&mut refdata.refs[d].to_ascii_vec());
                        concat.append(&mut refdata.refs[j].to_ascii_vec());

                        // Align the VDJ sequence on the contig to the reference concatenation.

                        let al = aligner.semiglobal(&tig, &concat);
                        let mut m = 0;
                        let mut pos = al.xstart;
                        let mut rpos = (al.ystart as isize) - (refdata.refs[v].len() as isize);
                        let mut count = 0;
                        let mut edits = Vec::<String>::new();
                        while m < al.operations.len() {
                            let n = next_diff(&al.operations, m);
                            match al.operations[m] {
                                Match => {
                                    pos += 1;
                                    rpos += 1;
                                }
                                Subst => {
                                    if pos >= start && pos < stop {
                                        count += 1;
                                        edits.push(format!("S{}", rpos));
                                    }
                                    pos += 1;
                                    rpos += 1;
                                }
                                Del => {
                                    if pos >= start && pos < stop {
                                        count += 1;
                                        edits.push(format!("D{}:{}", rpos, n - m));
                                    }
                                    pos += n - m;
                                    m = n - 1;
                                }
                                Ins => {
                                    if pos >= start && pos < stop {
                                        count += 1;
                                        edits.push(format!("I{}:{}", rpos, n - m));
                                    }
                                    rpos += (n - m) as isize;
                                    m = n - 1;
                                }
                                _ => {}
                            };
                            m += 1;
                        }
                        let edits = format!("{}", edits.iter().format("•"));
                        results.push((count, d, edits));
                    }
                }
                results.sort();
                fwriteln!(log, "");
                for j in 0..min(5, results.len()) {
                    let comp = results[j].0;
                    let t = results[j].1;
                    let edits = results[j].2.clone();
                    fwriteln!(
                        log,
                        "[{}] {} = {} = {}, score = {}, edit = {}",
                        j + 1,
                        refdata.id[t],
                        refdata.name[t],
                        refdata.refs[t].to_string(),
                        comp,
                        edits
                    );
                }
            }
        }

        // Handle ext.

        if ext {
            let seq_rc = seq.rc();
            fwriteln!(log, "\nRC ANNOTATION VERSUS PLAIN REFERENCE\n");
            print_annotations(&seq_rc, &refdata, &mut log, false, false, verbose);

            // Annotate using the extended reference.

            let mut refdatax = RefData::new();
            make_vdj_ref_data(&mut refdatax, imgt, &species, true, is_tcr, is_bcr);
            fwriteln!(log, "\nFW ANNOTATION VERSUS EXTENDED REFERENCE\n");
            print_annotations(&seq, &refdatax, &mut log, false, false, verbose);
            fwriteln!(log, "\nRC ANNOTATION VERSUS EXTENDED REFERENCE\n");
            print_annotations(&seq_rc, &refdatax, &mut log, false, false, verbose);
        }

        // Print.

        print!("{}\n", stringme(&log));
    }
}
