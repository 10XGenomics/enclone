// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Usage:
//
// denovo species region <optional args>
//
// species: may be human or mouse or dog or any string having a unique match to a record
// in a directory having genomes, see genomes_list.
// or all to do everything
//
// region is V or D or J or C, or any combination of those, e.g. VDJC
// For V, the 5'-UTR TAG is included as a separate record.
//
// optional extras for V segments:
// - trans: show transition
// - dna: print bases (however this also perturbs results)
// - all: print all
// - pwm: print PWM for start site
// - fasta: print fasta instead of other stuff and print all
//          CAN PRINT MULTIPLE V FASTA SEQS PER FWR3 START SITE.
//          FOR NOW PRINTS BOTH FW AND RC FOR D GENES.
// - aa: print amino acid fasta (only implemented for J genes)
// - fasta_store: same as fasta but store in somewhere/denovo_ref
// - 10x: use 10x style for fasta header lines
//
// Currently we use the reference to disambiguate groups which differ only in how far left the
// leader extends.  In such cases we do not have a good basis for choosing and do not believe
// that the reference was chosen on such a basis either.
//
// optional extra arg = amb: don't disambiguate such cases
//
// This code is not fully functional at present because it refers to directories that don't
// in general exist.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// LIST OF CONDITIONS/HEURISTICS
//
// Two numbers are shown for each test (unless turning it off won't work at all).  The first number
// is the number of genes that are missed if the test is turned off, whereas the second is the
// number of fake genes that are added.  Since these numbers were computed before turning off some
// tests, they may no longer be accurate.
//
// 1.  required  FWR3 must end with YYC or YFC or YHC or FYC or YIC
// 2.  required  pick single best stopping point for FWR3
// 3.  required  seach for exon one start between -175 and 520 bases behind the ORF
// 4.  required  require that exon one is between 19 and 73 bases long
// 5.  crashes   require that the sequence following exon one is GT or GC
// 6.  required  intron has size between 76 and 443 bases
// 7.  0   5     exon two start score has to be at least -9.0
// 8.  crashes   the intron sequence preceeding exon two must be CAG or TAG or AAG
// 9.  1   1     the last base of exon 1 may not be C
// 10. crashes   the length of leader..FWR3 in amino acids must be between 102 and 125
// 11. crashes   the leader must have length between 17 and 28 amino acids
// 12. 2   7     score1 of the leader must be at least 4
// 13. 0   4     if the second intron base is C, then the score1 must be at least 6
// 14. required  all features must be computable
// 15. crashes   FWR1 must have length between 20 and 26 amino acids
// 16. epsilon   for IGL, the FWR1 must have length at most 22
// 17. crashes   score1 must be at least 1.8
// 18. 0   1     muscore >= 12.0
// 19. 0  40     nlogp2 >= 58.0
//               (following conditions compare s1 to s2)
// 20. 1  10     l2(s1) >= 10 * l2(s2) & l1(s1) >= l1(s2) & post(s1) != C = kill s2
// 21. 0   2     l2(s1) >= l2(s2) & l1(s1) >= l1(s2) & post(s1) = T & post(s2) = C = kill s2
// 22. 2   2     intron(s2) >= 150 + intron(s1) & l2(s1) >= 10 * l2(s2)
//               & start_s(s1) >= start_s(s2) = kill s2
// 23. 0  10     have proline at reverse pos 5 or 6 on the leader is better than not having one
// 24. 0  18     reverse leader score test 2
// 25. 0  31     exon one same start test
// 26. 0   5     recombination site filtering

/*
PREBUILD: ran this:
% cd ensembl/release-94/fasta/mus_musculus/dna
% shrink_fasta Mus_musculus.GRCm38.dna.toplevel.fa.gz Mus_musculus.GRCm38.dna.toplevel.trunc_1000.fa
% cd ensembl/release-94/fasta/homo_sapiens/dna
% shrink_fasta Homo_sapiens.GRCh38.dna.toplevel.fa.gz Homo_sapiens.GRCh38.dna.toplevel.trunc_1000.fa
% convert_fasta release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.trunc_1000.fa ensembl/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.trunc_1000.fa.binary
% convert_fasta ensembl/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.trunc_1000.fa ensembl/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.trunc_1000.fa.binary
*/

/*

Notes on mouse edge cases:

[7.1] 131 ==> (17) %T% MYRSTLELLFHLFSSSS 🌸 DVLLTQTPLFLPVSLGDQASISCSSSQSLVHSNGNSYLEWHLQKSGQSLQLLIYEVSKRHSGVPDRFSGSGSGTDFTLKISRVEPEDLGVYYC, rev_ldr_score=2.41, start_score=-16.31 nlogp=35.88, nlogp2=24.07, ct = IGK (93) IGK
==> short leader, tail matches IGKV1-115*01 but that has a stop codon

[9.1] 218 ==> (20) %T% MRCLAEFLRLLVLWIPGATG 🌸 DIVMTQAAPSVPANPGESVSISCRSSKSLLHSSGNTYLYWFLQRPGQSPQLLIYYISNLASGVPDRFSGSGSGTDFTLRISRVEAEDVGVYYC, rev_ldr_score=2.94, start_score=-11.70 nlogp=32.47, nlogp2=24.07, ct = IGK (93) IGK
==> does not seem to be real but can't figure out what's wrong with it

[5.1] 123 ==> (20) %T% MTAPTQLLGLLLLWIIGTRC 🌸 DIRMIQTPASLSGSLGESVTITCQASQDIGKSLLWYQQKTGNPPKILIYTTSNLADGISSRVSGSGSGTQFFLKFSSLKPEDTATYYC, rev_ldr_score=3.01, start_score=-14.58 nlogp=49.41, nlogp2=36.38, ct = IGK (88) IGK
==> appears in the data exactly, but always frameshifted; could something be wrong on the far end?

[3.1] 122 ==> (20) %T% MSVPTQLLALLLLWLTDARC 🌸 DIQVTQSPASLSAPVGESVSITCKASEEIYSALNWYQQKPGKSPQLLIYYATSLGDDVPSRFSGSKSGTQYSLKISSLQPEDLATYYC, rev_ldr_score=2.84, start_score=-14.73 muscore=-3.38 nlogp=36.88, nlogp2=23.38, ct = IGK (88) IGK
==> equals IMGT IGKV12-47*01, but not in our data; annotated as pseudogene

[1.1] 114 ==> (20) %T% MRAPALFLGILLLWFPGARC 🌸 DIQMTQSPSSMSASLRERVSLTCQASQGINGNLHWFQQKSGGTLKHLIYSTSNLDSGVPSRFSGSGSGSDYSLTISSLESEDFAVYYC, rev_ldr_score=3.21, start_score=-13.83 muscore=-1.00 nlogp=36.19, nlogp2=25.57, ct = IGK (88) IGK
==> the same except for left extension, not in data, IGK but not close to a reference

[5.1] 131 ==> (17) %T% MRTPAPFLGLLLFCARC 🌸 DVLMTQSPSSLSASLGERVSLTCQASQGISNNLNWYQQTPGKAPRLLIYDASKLEDGVPSRFSGTGYRTDFNFTISSLEEEDVATYFC, rev_ldr_score=2.11, start_score=-12.23 muscore=-3.78 nlogp=38.02, nlogp2=28.76, ct = IGK (88) IGK
==> matches IMGT IGKV11-106*01, listed as pseudogene, except for a two-base indel and a different
leader; in our data, except for the leader; seq in data has leader length 14, which is very low,
and fails GT|GC test after exon 1; note that the leader that's shown is short too

[14.2] 144 ==> (20) %T% MRIRGQLLGLLVLWITGVLC 🌸 DIQMTQSSSFLSASLGDHLTINCRASKDINKYFAWVQQKPRKAPRMLIHFASTLLPGVPEKFSGSGSGTDFSLTIRNIESEDIAMYYC, rev_ldr_score=3.60, start_score=-12.14 muscore=2.25 nlogp=55.62, nlogp2=40.46, ct = IGK HIGH (88) IGK
==> somewhat high nlogp2 score

*/

// Notes on human pseudogenes.
//
// [1] IMGT >546|IGKV1/OR2-0*01 googles nonfunctional, but we find
// MRAPTQLLGLLVLWLPGARC
// 🌸 DIQMTQSPSSLSASVGDRVTITCRASQGISNNLNWYQQKPGKTPKLLIYAASSLQSGIPSRFSDSGSGADYTLTIRSLQPEDFATYYC
// which differs from that at two bases (and one amino acid).
//
// The IMGT sequence probably does not occur in any 10x data.  The command
// enclone BI=1-12 IMGT RE SEGN=546 ACCEPT_BROKEN ALLOW_INCONSISTENT
// yields four cells but the sequences are very mutated, suggesting misalignment.
//
// Evidence that the gene might be defective is scant:
// 1. Its FWR3 has some rare amino acids
// (a) D in position 8 occurs only once in all IMGT mammalian data
// (b) R in position 20 occurs only twice.
// 2. Its initial base sequence ATGAGGGCCCCC seems like it might hairpin.  Yet it does initiate 16
// of the 5795 mammalian reference sequences.  (Possibly these are all pseudogenes; hard to know.)
//
// The gene blasts to about the right place on chr12 (and elsewhere).
//
// [2] IMGT 278|IGHV1/OR15-9*01 is
// MGWTWRILFLVVIAAGAQS 🌸
// QVQLMQSGAEVKKPGASVRISCKASGYTFTSYCMHWVCQAHAQGLEWMGLVCPSDGSTSYAQKFQGRVTITRDTSMGTAYMELSSLRSEDTAMYYC
// and the feature lengths are totally standard.  Haven't looked deeper than that.
// Doesn't appear in BI=6-12.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

use amino::{aa_seq, codon_to_aa};
use binary_vec_io::binary_read_vec_vec;
use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::AlignmentOperation::Del;
use bio::alignment::AlignmentOperation::Ins;
use debruijn::dna_string::DnaString;
use enclone_denovo::const_ighd::ighd_score2;
use enclone_denovo::make_fwr3_freqs::make_fwr3_freqs;
use enclone_denovo::mammalian_fixed_len::mammalian_fixed_len;
use enclone_denovo::mammalian_pwms::mammalian_pwms;
use enclone_denovo::vdj_features::{
    cdr1, cdr2, cdr3_start, fr1_start, fwr1, fwr2, fwr3, score_fwr3_at_end,
};
use fasta_tools::read_fasta_contents_into_vec_dna_string_plus_headers;
use io_utils::{fwrite, fwriteln, open_for_write_new};
use itertools::Itertools;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::{add_commas, stringme, strme, TextUtils};
use superslice::Ext;
use vdj_ann_ref::{human_ref, mouse_ref};
use vector_utils::{bin_member, bin_position, erase_if, reverse_sort, sort_sync2, unique_sort};

// copied from tenkit2/pack_dna.rs:

pub fn reverse_complement(x: &mut Vec<u8>) {
    x.reverse();
    for j in 0..x.len() {
        x[j] = match x[j] {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => x[j],
        }
    }
}

fn make_upper_acgtn(x: &mut u8) {
    if *x == b'A' {
    } else if *x == b'a' {
        *x = b'A';
    } else if *x == b'C' {
    } else if *x == b'c' {
        *x = b'C';
    } else if *x == b'G' {
    } else if *x == b'g' {
        *x = b'G';
    } else if *x == b'T' {
    } else if *x == b't' {
        *x = b'T';
    } else {
        *x = b'N';
    }
}

fn _make_upper_acgt(x: &mut u8) {
    if *x == b'A' {
    } else if *x == b'a' {
        *x = b'A';
    } else if *x == b'C' {
    } else if *x == b'c' {
        *x = b'C';
    } else if *x == b'G' {
    } else if *x == b'g' {
        *x = b'G';
    } else {
        *x = b'T';
    }
}

// score exon 2 start by base frequency in human and mouse data

fn exon2_start_score(x: &[u8]) -> f64 {
    let counts = [
        [41, 16, 342, 29],
        [33, 54, 56, 285],
        [84, 1, 255, 88],
        [1, 199, 2, 226],
        [32, 355, 8, 33],
        [145, 205, 25, 53],
        [179, 27, 116, 106],
        [87, 100, 174, 67],
        [15, 0, 118, 295],
        [3, 206, 212, 7],
    ];
    let mut p = 1.0;
    for i in 0..counts.len() {
        let total: usize = counts[i].iter().sum();
        let n;
        if x[i] == b'A' {
            n = counts[i][0];
        } else if x[i] == b'C' {
            n = counts[i][1];
        } else if x[i] == b'G' {
            n = counts[i][2];
        } else {
            n = counts[i][3];
        }
        p *= n as f64 / total as f64;
    }
    p.log10()
}

// lscore1 of an IG leader segment: the max length of a substring consisting of two or more
// letters from the alphabet VFAILT, except that such substrings are joined over single exceptions,
// and scored as the total number of letters from the alphabet in the string.
// Example
//       12345-678     (score = 8)
// MDWIWRILFLVGAATGAHS.
//
// Now capped at 7.

pub fn lscore1(x: &[u8]) -> usize {
    let all = [b'V', b'F', b'A', b'I', b'L', b'M', b'T'];
    let mut good = vec![false; x.len()];
    for i in 0..x.len() {
        for j in 0..all.len() {
            if x[i] == all[j] {
                good[i] = true;
                break;
            }
        }
    }
    let mut m = 0;
    for i in 0..x.len() - 1 {
        if !good[i] || !good[i + 1] {
            continue;
        }
        let mut j = i + 2;
        let mut count = 2;
        while j < x.len() {
            if good[j] {
                j += 1;
                count += 1;
            } else if j + 2 < good.len() && good[j + 1] && good[j + 2] {
                j += 3;
                count += 2;
            } else {
                break;
            }
        }
        m = max(m, count);
    }
    std::cmp::min(m, 7)
}

// lscore2 of an IG leader segment: number of times the length has been seen in reference Vs.

pub fn lscore2(x: &[u8]) -> usize {
    let n = x.len();
    if n == 19 {
        1256
    } else if n == 20 {
        321
    } else if n == 22 {
        142
    } else if n == 18 {
        48
    } else if n == 17 {
        14
    } else if n == 21 {
        9
    } else if n == 23 {
        6
    } else if n == 26 {
        4
    } else if n == 24 {
        3
    } else if n == 28 {
        2
    } else if n == 25 {
        1
    } else if n == 14 {
        1
    } else if n == 13 {
        1
    } else {
        0
    }
}

fn main() {
    PrettyTrace::new().on();

    // Get the species and load the genome reference.  We read the entire fasta file into a
    // Vec<Vec<u8>>, with records alternating between headers and bases.  This may not be
    // correctly handled in subsequent code.

    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("\nusage = denovo genome region (with additional optional args)\n");
        std::process::exit(1);
    }
    let species = &args[1];
    let region_type = &args[2];
    let mut use_v = false;
    let mut use_d = false;
    let mut use_j = false;
    let mut use_c = false;
    for c in region_type.chars() {
        if c == 'V' {
            use_v = true;
        } else if c == 'D' {
            use_d = true;
        } else if c == 'J' {
            use_j = true;
        } else if c == 'C' {
            use_c = true;
        } else {
            eprintln!("\nIllegal character {} in region type.\n", c);
            std::process::exit(1);
        }
    }
    let genomes0 = std::fs::read_dir("somewhere/genomes").unwrap();
    let mut genomes = Vec::<String>::new();
    for f in genomes0 {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        genomes.push(f.to_string());
    }

    // Launch on multiple genomes.

    if species == "all" {
        let (mut todo, mut names) = (Vec::<String>::new(), Vec::<String>::new());
        if species == "all" {
            todo = vec!["human".to_string(), "mouse".to_string(), "dog".to_string()];
            names = todo.clone();
            for i in 0..genomes.len() {
                todo.push(genomes[i].rev_after("/").rev_before(":").to_string());
                names.push(genomes[i].rev_after("/").to_string());
            }
        }
        println!();
        for i in 0..todo.len() {
            println!("running denovo {}", names[i]);
            let o = Command::new("denovo")
                .arg(&todo[i])
                .args(&args[2..])
                .output()
                .expect("failed to execute denovo");
            if o.status.code().unwrap() != 0 {
                println!("\nOOPS FAILED!");
                std::process::exit(1);
            }
            let m = String::from_utf8(o.stdout).unwrap();
            println!("{}", m);
        }
        std::process::exit(0);
    }

    // Continuing now on a single genome.

    let mut fasta_file = String::new();
    let order;
    if species == "human" {
        order = "Primates".to_string();
    } else if species == "mouse" {
        order = "Rodentia".to_string();
    } else if species == "dog" {
        order = "Carnivora".to_string();
    } else {
        let mut fns = Vec::<String>::new();
        for f in genomes.iter() {
            let sp = species.clone();
            if f.contains(&sp) {
                fns.push(f.clone());
            }
        }
        if fns.is_empty() {
            eprintln!("\nCan't find your species.\n");
            std::process::exit(1);
        }
        if fns.len() > 1 {
            eprintln!("\nThere are multiple matches for your species:\n");
            for i in 0..fns.len() {
                eprintln!("{}", fns[i]);
            }
            eprintln!();
            std::process::exit(1);
        }
        fasta_file = fns[0].clone();
        order = fasta_file.after(":").between(":", ":").to_string();
    }
    let id_name;
    if fasta_file.is_empty() {
        id_name = species.clone();
    } else {
        id_name = fasta_file.rev_after("/").rev_before(".").to_string();
    }
    let mut fasta_log = Vec::<u8>::new();
    let fasta_out_dir = "somewhere/denovo_ref";

    // Get the reference.

    let mut refx = Vec::<Vec<u8>>::new();
    // let t = Instant::now();
    if species == "human" || species == "mouse" {
        let root = "ensembl/release-94/fasta";
        let xref;
        if species == "human" {
            xref = format!(
                "{}/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.trunc_1000.fa.binary",
                root
            );
        } else {
            xref = format!(
                "{}/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.trunc_1000.fa.binary",
                root
            );
        }
        let mut f = std::fs::File::open(&xref).unwrap();
        binary_read_vec_vec(&mut f, &mut refx).unwrap();
    } else if species == "dog" {
        fasta_file = "somewhere/genomes/GCA_000002285.4:some_mam:Carnivora:Canis_lupus_familiaris:Dog.vecvec_u8".to_string();
        let mut f = std::fs::File::open(&fasta_file).unwrap();
        binary_read_vec_vec(&mut f, &mut refx).unwrap();
    } else {
        let mut f = std::fs::File::open(&fasta_file).unwrap();
        binary_read_vec_vec(&mut f, &mut refx).unwrap();
    }
    // println!("used {:.2} seconds loading genome", elapsed(&t));

    // Parse other arguments.

    let mut show_transition = false;
    let mut print_bases = false;
    let mut print_all = false;
    let mut amb = false;
    let mut pwm = false;
    let mut print_fasta = false;
    let mut store_fasta = false;
    let mut print_aa = false;
    let mut tenx = false;
    for i in 3..args.len() {
        if args[i] == "trans" {
            show_transition = true;
        } else if args[i] == "dna" {
            print_bases = true;
        } else if args[i] == "all" {
            print_all = true;
        } else if args[i] == "pwm" {
            pwm = true;
        } else if args[i] == "amb" {
            amb = true;
        } else if args[i] == "fasta" {
            print_fasta = true;
        } else if args[i] == "fasta_store" {
            store_fasta = true;
        } else if args[i] == "10x" {
            tenx = true;
        } else if args[i] == "aa" {
            print_aa = true;
        } else {
            eprintln!("\nIllegal argument {}.\n", args[i]);
            std::process::exit(1);
        }
    }
    let mut fasta_count = 0;

    // Empirical constants.

    const WALK_BACK: usize = 520; // distance to walk back from ORF start to find start codon
                                  // LOW1 could probably be quite a bit higher
    const LOW1: usize = 19; // minimum length of exon 1 in bases (from start codon)
    const HIGH1: usize = 73; // maximum length of exon 1 in bases (from start codon)
    const LOW_INTRON: usize = 76; // minimum length of intron in bases
    const HIGH_INTRON: usize = 443; // maximum length of intron in bases
    const MIN_LEN: usize = 102; // minimum length from start codon to end of FWR3
    const MAX_LEN: usize = 125; // maximum length from start codon to end of FWR3
    const MIN_SCORE: f64 = 6.75; // minimum CDR3 score
    const MIN_FWR1: usize = 20; // minimum FWR1 length
    const MAX_FWR1: usize = 26; // maximum FWR1 length
    const MIN_LEADER: usize = 17; // minimum leader length
    const MAX_LEADER: usize = 28; // maximum leader length
    const TAG_LEN: usize = 40; // length of upstream region used as tag
    const MOTIF1: &[u8; 7] = b"CACAGTG"; // start motif for initial upstream tags
    const MOTIF2: &[u8; 9] = b"ACAAAAACC"; // mid motif for upstream tags

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Order genome records by length, to speed up subsequent parallel loops.  This has a huge
    // effect.  Weirdly it is this order and not reverse order that works.  Note that we do not
    // track the index changes.

    {
        let mut refx2 = Vec::<Vec<u8>>::new();
        let mut lens = Vec::<(usize, usize)>::new();
        for i in (0..refx.len()).step_by(2) {
            lens.push((refx[i + 1].len(), i));
        }
        lens.sort_unstable();
        for i in 0..lens.len() {
            refx2.push(refx[lens[i].1].clone());
            refx2.push(refx[lens[i].1 + 1].clone());
        }
        refx = refx2;
    }

    // Upper case genome bases.

    refx.par_iter_mut().for_each(|res| {
        if res[0] != b'>' {
            for j in 0..res.len() {
                if res[j] == b'a' {
                    res[j] = b'A';
                } else if res[j] == b'c' {
                    res[j] = b'C';
                } else if res[j] == b'g' {
                    res[j] = b'G';
                } else if res[j] == b't' {
                    res[j] = b'T';
                }
            }
        }
    });

    // Break long genome records into shorter overlapping pieces.  Also drop headers.

    let mut refy = Vec::<Vec<u8>>::new();
    let mut originy = Vec::<(usize, usize)>::new();
    for k in (0..refx.len()).step_by(2) {
        let i = k + 1;
        if refx[i].len() <= MAX_TIG {
            refy.push(refx[i].clone());
            originy.push((k, 0));
        } else {
            let mut stops = Vec::<usize>::new();
            for j in (0..refx[i].len()).step_by(MAX_TIG) {
                let mut start = 0;
                if j > 0 {
                    start = stops[j / MAX_TIG - 1] - TIG_OVERLAP;
                }
                let stop = min(start + MAX_TIG, refx[i].len());
                stops.push(stop);
                refy.push(refx[i][start..stop].to_vec());
                originy.push((k, start));
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // BEGIN ANALYSIS OF C REGIONS
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Constant region analysis limitations.
    // 1. There is no regression testing in place, either to prior or published results.
    // 2. Gene naming (e.g. IGHM versus IGHD) is entirely theoretical, based on consistency
    //    with genomes for which the answer is known, and might be wrong.
    // 3. We only report the first part of CH1, so we don't see the biology of what comes
    //    after, and for species having functional genes missing CH1 (e.g. camelid IGHD), we
    //    don't see the gene at all.
    // 4. For the 200 mammal assemblies, we are not necessarily extracting full power from the
    //    underlying data.  We might be better off with a graph assembly that did not drop small
    //    contigs, as we may be losing both sequence and order information.
    // 5. The alignment algorithm has some flaws:
    //    (a) It only allows deletions, and not insertions.  This likely gets us the wrong answer
    //        e.g. for bamboo rate, possibly also pronghorn (one allele) and hedgehog.
    //    (b) It requires the deletions occur at mod three boundaries.
    //    (c) It does not use dynamic programming and so is likely slow.
    // 6. For the PWMs, we change zero counts to two (which is arbitrary), and don't change one
    //    counts.
    // 7. There is likely more power to be gotten by having motifs/PWMs clade by clade.
    // 8. For Artiodactyla, if IGHM and IGHD are not on the same reference contig, we call them
    //    both IGHM.  Probable Examples include chevrotain and reindeer.  Possible approaches
    //    include looking at a longer pre sequence and seeing if there's a signal that
    //    distinguishes the two cases, perhaps in conjunction with clade-specific information.
    //    See also point 4.
    // 9. Many bats appear to be missing IGKC.  To reexamine using isoseq data.
    // 10. Abbreviate contig ids and show contig lengths.
    // 11. Motif finder should averge over sequences for a given genome rather than takine the
    //     first one.
    // 12. Motif finding should be iterated.

    // Read constant region motifs.

    let chains = ["IGHA", "IGHD", "IGHE", "IGHG", "IGHM", "IGKC", "IGLC"];
    let mut chain_motifs = vec![Vec::<Vec<u8>>::new(); chains.len()];
    let mut pre_max_errs = vec![0; chains.len()];
    let mut max_errs = vec![0; chains.len()];
    let mut pwmx = vec![Vec::<f32>::new(); chains.len()];
    let mut ins_start = vec![0; chains.len()];
    const IGH_CONST_DEL_PENALTY: f32 = 1.0;
    {
        let motifs = include_str!("../const_motifs");
        for line in motifs.lines() {
            if line.starts_with("IG") {
                let mut chain = None;
                for i in 0..chains.len() {
                    if line[0..4] == chains[i][..] {
                        chain = Some(i);
                    }
                }
                if chain.is_none() {
                    continue;
                }
                let chain = chain.unwrap();
                let ch = &line[0..4];
                let f = format!("{} = ", ch);
                if line.starts_with(&f) {
                    chain_motifs[chain].push(line.after(&f).to_string().as_bytes().to_vec());
                }
                let f = format!("{}_PRE_MAX = ", ch);
                if line.starts_with(&f) {
                    pre_max_errs[chain] = line.after(&f).force_usize();
                }
                let f = format!("{}_MAX = ", ch);
                if line.starts_with(&f) {
                    max_errs[chain] = line.after(&f).force_usize();
                }
                let f = format!("{}_PWM = ", ch);
                if line.starts_with(&f) {
                    let fields = line.after(&f).split(':').collect::<Vec<&str>>();
                    for i in 0..fields.len() {
                        let mut c = Vec::<usize>::new();
                        let counts = fields[i].split(',').collect::<Vec<&str>>();
                        for j in 0..counts.len() {
                            c.push(counts[j].force_usize());
                        }
                        for j in 0..c.len() {
                            if c[j] == 0 {
                                c[j] = 2;
                            }
                        }
                        let s: usize = c.iter().sum();
                        for j in 0..c.len() {
                            pwmx[chain].push(-(c[j] as f32 / s as f32).log10());
                        }
                        pwmx[chain].push(IGH_CONST_DEL_PENALTY);
                    }
                }
            }
        }
        for i in 0..chains.len() {
            if !chain_motifs[i].is_empty() {
                let mut sep = 0;
                for j in 0..chain_motifs[i][0].len() {
                    if chain_motifs[i][0][j] == b'|' {
                        sep += 1;
                        if sep == 2 {
                            break;
                        }
                    } else {
                        ins_start[i] += 1;
                    }
                }
            }
        }
    }

    // Define a data structure that can be used to efficiently determine if a base matches
    // the chain motif.  It has a "one" for "wrong" bases.

    let n = pwmx[0].len();
    let mut chain_motifs_penalty = vec![vec![vec![1_u8; 4]; n]; chains.len()];
    for i in 0..chains.len() {
        for k in 0..chain_motifs[i][0].len() {
            if chain_motifs[i][0][k] == b' ' {
                for m in 0..4 {
                    chain_motifs_penalty[i][k][m] = 0;
                }
            }
        }
        for j in 0..2 {
            for k in 0..chain_motifs[i][j].len() {
                let c = chain_motifs[i][j][k];
                if c == b'A' {
                    chain_motifs_penalty[i][k][0] = 0;
                } else if c == b'C' {
                    chain_motifs_penalty[i][k][1] = 0;
                } else if c == b'G' {
                    chain_motifs_penalty[i][k][2] = 0;
                } else if c == b'T' {
                    chain_motifs_penalty[i][k][3] = 0;
                }
            }
        }
    }

    // Define the maximum allowed deletion, relative to the PWM, measured in amino acids.  Note
    // confusing deletion/insertion language.

    const MAX_DEL_LEN: usize = 6;

    // Score.

    #[derive(Clone, PartialOrd, PartialEq)]
    struct Chit {
        pub tig: usize,
        pub pass: usize,
        pub start: usize,
        pub score: f32,
        pub log: Vec<u8>,
        pub bases: Vec<u8>,
        pub ins_pos: usize,
        pub ins_len: usize,
        pub chain: usize,
        pub region: String,
    }
    let mut results = Vec::<(usize, Vec<Chit>)>::new();
    if use_c {
        for i in 0..refy.len() {
            results.push((i, Vec::<Chit>::new()));
        }
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut r = refy[i].clone();
        for j in 0..r.len() {
            make_upper_acgtn(&mut r[j]);
        }
        for pass in 0..2 {
            if pass == 1 {
                reverse_complement(&mut r);
            }
            let or;
            if pass == 0 {
                or = b'+';
            } else {
                or = b'-';
            }
            let mut s = r.clone();
            for j in 0..r.len() {
                if r[j] == b'A' {
                    s[j] = 0;
                } else if r[j] == b'C' {
                    s[j] = 1;
                } else if r[j] == b'G' {
                    s[j] = 2;
                } else {
                    s[j] = 3;
                }
            }
            const MIN_AA_IGH: usize = 60;
            const MAX_MIS: i32 = 2;
            for frame in 0..3 {
                let aa = aa_seq_safe(&r, frame);
                let mut finds = Vec::<(f32, usize, Vec<u8>, usize, usize, usize)>::new();
                let mut j = 0;
                if aa.len() < MIN_AA_IGH {
                    continue;
                }
                while j <= aa.len() - MIN_AA_IGH {
                    if aa[j] == b'*' {
                        j += 1;
                        continue;
                    }
                    let mut k = j + 1;
                    while k < aa.len() && aa[k] != b'*' {
                        k += 1;
                    }
                    if k - j >= MIN_AA_IGH {
                        for m in j..=k - MIN_AA_IGH {
                            let start = 3 * m + frame;

                            // Traverse each of the chain types.

                            for chain in 0..chains.len() {
                                if start < ins_start[chain] {
                                    continue;
                                }

                                // Screen the bases before the gene start for a match to the motif.

                                let mut mis = 0;
                                for u in 0..12 {
                                    mis += chain_motifs_penalty[chain][u]
                                        [s[start - 14 + u] as usize]
                                        as i32;
                                    if mis > MAX_MIS {
                                        break;
                                    }
                                }
                                if mis > MAX_MIS {
                                    continue;
                                }

                                // Now score the position.

                                let width = pwmx[0].len() / 5;
                                let bb = &s[start - ins_start[chain]..start + width];
                                let (mut ins_pos, mut ins_len) = (0, 0);
                                let score = ighd_score2(
                                    bb,
                                    &pwmx[chain],
                                    MAX_DEL_LEN,
                                    ins_start[chain],
                                    &mut ins_pos,
                                    &mut ins_len,
                                );
                                if score <= 70.0 {
                                    let pre = strme(&r[start - ins_start[chain]..start - 2]);
                                    let mid = strme(&r[start - 2..start]);
                                    let mut beast = r[start..start + 120].to_vec();
                                    for _ in 0..ins_len {
                                        beast.insert(ins_pos - ins_start[chain], b'.');
                                    }
                                    beast.truncate(beast.len() - ins_len);
                                    let acc;
                                    let mut name;
                                    if species == "human" || species == "mouse" || species == "dog"
                                    {
                                        acc = "           ".to_string();
                                        name = species.clone();
                                    } else {
                                        let f = fasta_file.rev_after("/");
                                        acc = f.between("_", ":").to_string();
                                        name = f
                                            .after(":")
                                            .after(":")
                                            .after(":")
                                            .between(":", ".")
                                            .to_string();
                                    }
                                    name = name.replace("_", " ");
                                    let tig = originy[i].0;
                                    let best = format!(
                                        "{}{}.{}  {}  {}|{}|{}  {:.1}  {}  {}",
                                        or as char,
                                        tig,
                                        add_commas(start),
                                        chains[chain],
                                        pre,
                                        mid,
                                        strme(&beast),
                                        score,
                                        acc,
                                        name
                                    );
                                    finds.push((
                                        score,
                                        m,
                                        best.as_bytes().to_vec(),
                                        chain,
                                        ins_pos,
                                        ins_len,
                                    ));
                                }
                            }
                        }
                    }
                    j = k;
                }

                // Save results.

                for m in 0..finds.len() {
                    let mut log = Vec::<u8>::new();
                    let tig = originy[i].0;
                    let start = originy[i].1 + frame + 3 * finds[m].1;
                    fwriteln!(log, "{}", strme(&finds[m].2));
                    let rstart = frame + 3 * finds[m].1;
                    let chain = finds[m].3;
                    let ins_pos = finds[m].4;
                    let ins_len = finds[m].5;
                    res.1.push(Chit {
                        tig,
                        pass,
                        start,
                        score: finds[m].0,
                        log,
                        bases: r[rstart - 14..rstart + 120].to_vec(),
                        ins_pos,
                        ins_len,
                        chain,
                        region: chains[chain].to_string(),
                    });
                }
            }
        }
    });
    let mut chits = Vec::<Chit>::new();
    for i in 0..results.len() {
        chits.append(&mut results[i].1.clone());
    }
    chits.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // For Artiodactyla, correct the "second" IGHM to be IGHD.

    const MAX_IGHM_IGHD_DIFF: usize = 15_000;
    if order == "Artiodactyla" {
        for i in 1..chits.len() {
            if chits[i].region == "IGHM"
                && chits[i - 1].region == "IGHM"
                && chits[i].tig == chits[i - 1].tig
                && chits[i].pass == chits[i - 1].pass
                && chits[i].start - chits[i - 1].start < MAX_IGHM_IGHD_DIFF
            {
                chits[i].region = "IGHD".to_string();
                chits[i].log = stringme(&chits[i].log)
                    .replace(" IGHM", " IGHD")
                    .as_bytes()
                    .to_vec();
            }
        }
    }

    // Remove inferior hits.  They must differ significantly in score and either:
    // (a) have nearly the same start position or
    // (b) have the same gene and differ signifantly in their last eight pre bases.

    const MIN_IGH_DIFF: f32 = 5.0;
    let mut to_delete = vec![false; chits.len()];
    for i1 in 0..chits.len() {
        for i2 in 0..chits.len() {
            if chits[i2].score >= MIN_IGH_DIFF + chits[i1].score {
                if i2 > i1
                    && chits[i2].tig == chits[i1].tig
                    && chits[i2].pass == chits[i1].pass
                    && chits[i2].start - chits[i1].start <= 3 * MAX_DEL_LEN
                {
                    to_delete[i2] = true;
                }
                if chits[i2].region == chits[i1].region {
                    let mut pre_diffs = 0;
                    for j in 4..12 {
                        if chits[i1].bases[j] != chits[i2].bases[j] {
                            pre_diffs += 1;
                        }
                    }
                    if pre_diffs >= 3 {
                        to_delete[i2] = true;
                    }
                }
            }
        }
    }
    erase_if(&mut chits, &to_delete);

    // Print.

    if !print_fasta && !store_fasta {
        let mut last_eq = 0;
        for i in 0..chits.len() {
            let eq = chits[i].log.iter().position(|&r| r == b'I').unwrap();
            last_eq = max(last_eq, eq);
        }
        for i in 0..chits.len() {
            let eq = chits[i].log.iter().position(|&r| r == b'I').unwrap();
            let add = last_eq - eq;
            for _ in 0..add {
                chits[i].log.insert(eq, b' ');
            }
        }
        for i in 0..chits.len() {
            print!("{}", strme(&chits[i].log));
        }
    } else {
        let mut count = vec![0; chains.len()];
        for i in 0..chits.len() {
            for j in 0..chains.len() {
                if chains[j] == chits[i].region {
                    count[j] += 1;
                    if !tenx {
                        fwriteln!(fasta_log, ">{}{}", chits[i].region, count[j]);
                    } else {
                        fasta_count += 1;
                        let gene = format!("{}{}", chits[i].region, count[j]);
                        fwriteln!(
                            fasta_log,
                            ">{}|{} enclone|{}|C-REGION|IG|{}|{}|00",
                            fasta_count,
                            gene,
                            gene,
                            strme(&gene.as_bytes()[0..3]),
                            strme(&gene.as_bytes()[3..]),
                        );
                    }
                    fwriteln!(fasta_log, "{}", strme(&chits[i].bases[12..]));
                }
            }
        }
    }

    // horse references
    // - https://pubmed.ncbi.nlm.nih.gov/15322185
    // - https://www.jimmunol.org/content/jimmunol/173/5/3230.full.pdf

    pub fn aa_seq_safe(x: &Vec<u8>, start: usize) -> Vec<u8> {
        let mut a = Vec::<u8>::new();
        if x.len() >= 3 {
            for j in (start..x.len() - 3 + 1).step_by(3) {
                if x[j] == b'-' && x[j + 1] == b'-' && x[j + 2] == b'-' {
                    a.push(b'-');
                } else if x[j] == b'N' || x[j + 1] == b'N' || x[j + 2] == b'N' {
                    a.push(b'*');
                } else {
                    a.push(codon_to_aa(&x[j..j + 3]));
                }
            }
        }
        a
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // PWM for start codon.

    let start_pwm = [
        [66, 20, 36, 18],
        [22, 33, 67, 18],
        [22, 10, 57, 51],
        [47, 66, 13, 14],
        [27, 73, 9, 31],
        [30, 38, 19, 53],
        [17, 107, 12, 4],
        [126, 1, 12, 1],
        [14, 92, 25, 9],
        [6, 121, 6, 7],
        [140, 0, 0, 0],
        [0, 0, 0, 140],
        [0, 0, 140, 0],
        [20, 2, 115, 3],
        [85, 38, 13, 4],
        [22, 73, 41, 4],
        [29, 23, 10, 78],
        [9, 18, 49, 64],
        [9, 35, 73, 23],
        [55, 34, 36, 15],
        [4, 62, 43, 31],
        [11, 34, 58, 37],
        [4, 61, 41, 34],
    ];

    // Define some pseudogenes.  Truncated after FWR3.

    let mut pseudos = Vec::new();
    pseudos.push((
        "human",
        "IGKV1/OR2-0*01",
        b"MRAPTQLLGLLVLWLPGARC\
        DIQMTQSPSSLSASVGDRVTITCRASQGISNNLNWYQQKPGKTPKLLIYAAPSLQSGIPSRFSDSGSGADYTLTIRSLQPEDFATYYC"
            .to_vec(),
    ));
    pseudos.push(("human", "GHV1/OR15-9*01",
        b"MGWTWRILFLVVIAAGAQS\
        QVQLMQSGAEVKKPGASVRISCKASGYTFTSYCMHWVCQAHAQGLEWMGLVCPSDGSTSYAQKFQGRVTITRDTSMGTAYMELSSLRSEDTAMYYC"
        .to_vec()));

    // Get the reference sequences, extract IG V seqs, truncate at end of FWR3, and build index.

    let mut to_ref = HashMap::<Vec<u8>, String>::new();
    if species == "human" || species == "mouse" {
        let r;
        if species == "human" {
            r = human_ref();
        } else {
            r = mouse_ref();
        }
        let mut zref = Vec::<DnaString>::new();
        let mut zheaders = Vec::<String>::new();
        read_fasta_contents_into_vec_dna_string_plus_headers(&r, &mut zref, &mut zheaders);
        for i in 0..zref.len() {
            if !zheaders[i].contains("V-REGION") {
                continue;
            }
            if !zheaders[i].contains("|IGH")
                && !zheaders[i].contains("|IGK")
                && !zheaders[i].contains("|IGL")
            {
                continue;
            }
            let chain_type;
            if zheaders[i].contains("|IGH") {
                chain_type = "IGH";
            } else if zheaders[i].contains("|IGK") {
                chain_type = "IGK";
            } else {
                chain_type = "IGL";
            }
            let aa = aa_seq(&zref[i].to_ascii_vec(), 0);

            // Skip short transcripts.  There is one, one of the two IGHV1-12 sequences, which has
            // length 60 amino acids, and should be deleted from the reference rather than
            // excluded here.

            if aa.len() < 100 {
                continue;
            }

            let cdr3 = cdr3_start(&aa, chain_type, false);
            to_ref.insert(aa[0..cdr3].to_vec(), zheaders[i].before(" ").to_string());
        }
    }

    // Make freqs.

    let freqs = make_fwr3_freqs();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Define fwr1 IG freqs.  These are frequencies for the first ten amino acids in FWR1 for
    // IG sequences, obtained by running fwr3_matrix.

    #[allow(dead_code)]
    let mut freq1 = Vec::<Vec<(usize, u8)>>::new();
    freq1.push(vec![
        (817, b'Q'),
        (494, b'E'),
        (308, b'D'),
        (93, b'S'),
        (51, b'A'),
        (15, b'N'),
        (7, b'R'),
        (5, b'V'),
        (5, b'G'),
        (4, b'T'),
        (4, b'L'),
        (3, b'H'),
        (1, b'W'),
        (1, b'P'),
        (1, b'K'),
    ]);
    freq1.push(vec![
        (881, b'V'),
        (415, b'I'),
        (144, b'S'),
        (113, b'A'),
        (58, b'Y'),
        (47, b'P'),
        (33, b'E'),
        (30, b'T'),
        (19, b'L'),
        (16, b'Q'),
        (14, b'M'),
        (13, b'N'),
        (13, b'F'),
        (5, b'G'),
        (3, b'R'),
        (2, b'H'),
        (2, b'D'),
        (1, b'K'),
    ]);
    freq1.push(vec![
        (884, b'Q'),
        (553, b'V'),
        (98, b'K'),
        (81, b'T'),
        (53, b'E'),
        (35, b'A'),
        (27, b'L'),
        (19, b'G'),
        (18, b'M'),
        (10, b'R'),
        (6, b'N'),
        (5, b'W'),
        (5, b'H'),
        (4, b'S'),
        (4, b'I'),
        (2, b'Y'),
        (2, b'P'),
        (2, b'D'),
        (1, b'C'),
    ]);
    freq1.push(vec![
        (1395, b'L'),
        (283, b'M'),
        (78, b'V'),
        (16, b'K'),
        (10, b'I'),
        (8, b'E'),
        (6, b'Q'),
        (5, b'P'),
        (4, b'F'),
        (2, b'T'),
        (1, b'W'),
        (1, b'G'),
    ]);
    freq1.push(vec![
        (813, b'T'),
        (417, b'V'),
        (274, b'Q'),
        (157, b'K'),
        (46, b'E'),
        (22, b'L'),
        (20, b'R'),
        (17, b'I'),
        (14, b'N'),
        (8, b'S'),
        (8, b'M'),
        (6, b'A'),
        (3, b'Y'),
        (3, b'D'),
        (1, b'H'),
    ]);
    freq1.push(vec![
        (1142, b'Q'),
        (622, b'E'),
        (26, b'S'),
        (7, b'K'),
        (6, b'V'),
        (3, b'P'),
        (1, b'R'),
        (1, b'I'),
        (1, b'F'),
    ]);
    freq1.push(vec![
        (1199, b'S'),
        (325, b'P'),
        (151, b'T'),
        (53, b'E'),
        (21, b'L'),
        (14, b'A'),
        (10, b'W'),
        (6, b'G'),
        (6, b'D'),
        (4, b'V'),
        (4, b'R'),
        (4, b'C'),
        (3, b'K'),
        (3, b'F'),
        (2, b'Y'),
        (1, b'Q'),
        (1, b'N'),
        (1, b'M'),
        (1, b'H'),
    ]);
    freq1.push(vec![
        (932, b'G'),
        (621, b'P'),
        (85, b'S'),
        (70, b'A'),
        (28, b'T'),
        (15, b'E'),
        (12, b'H'),
        (10, b'Q'),
        (7, b'D'),
        (5, b'R'),
        (4, b'V'),
        (4, b'N'),
        (4, b'L'),
        (4, b'F'),
        (3, b'I'),
        (2, b'W'),
        (2, b'K'),
        (1, b'C'),
    ]);
    freq1.push(vec![
        (490, b'S'),
        (406, b'G'),
        (358, b'P'),
        (357, b'A'),
        (98, b'L'),
        (21, b'K'),
        (18, b'E'),
        (17, b'D'),
        (16, b'T'),
        (12, b'F'),
        (6, b'R'),
        (5, b'V'),
        (2, b'H'),
        (1, b'Y'),
        (1, b'Q'),
    ]);
    freq1.push(vec![
        (527, b'G'),
        (334, b'S'),
        (285, b'V'),
        (257, b'E'),
        (106, b'L'),
        (58, b'T'),
        (55, b'A'),
        (49, b'I'),
        (44, b'D'),
        (41, b'F'),
        (12, b'M'),
        (8, b'Y'),
        (8, b'R'),
        (6, b'P'),
        (6, b'N'),
        (5, b'K'),
        (5, b'H'),
        (1, b'Q'),
        (1, b'C'),
    ]);

    // Add reverse leader frequencies.

    let mut ldrf = Vec::<Vec<(usize, u8)>>::new();
    ldrf.push(vec![
        (602, b'S'),
        (484, b'C'),
        (326, b'G'),
        (325, b'A'),
        (17, b'V'),
        (14, b'T'),
        (10, b'F'),
        (6, b'I'),
        (6, b'D'),
        (5, b'Y'),
        (5, b'R'),
        (5, b'E'),
        (2, b'P'),
        (1, b'L'),
        (1, b'H'),
    ]);
    ldrf.push(vec![
        (426, b'Q'),
        (343, b'L'),
        (206, b'R'),
        (191, b'W'),
        (183, b'H'),
        (111, b'S'),
        (71, b'V'),
        (64, b'T'),
        (52, b'C'),
        (33, b'N'),
        (32, b'D'),
        (19, b'Y'),
        (17, b'K'),
        (16, b'I'),
        (11, b'F'),
        (11, b'E'),
        (8, b'G'),
        (6, b'A'),
        (5, b'M'),
        (4, b'P'),
    ]);
    ldrf.push(vec![
        (829, b'V'),
        (518, b'S'),
        (244, b'A'),
        (86, b'T'),
        (74, b'I'),
        (16, b'G'),
        (14, b'F'),
        (7, b'M'),
        (7, b'C'),
        (4, b'Y'),
        (4, b'L'),
        (4, b'H'),
        (1, b'R'),
        (1, b'Q'),
    ]);
    ldrf.push(vec![
        (1310, b'G'),
        (109, b'W'),
        (92, b'C'),
        (69, b'D'),
        (49, b'S'),
        (49, b'A'),
        (33, b'M'),
        (31, b'Y'),
        (24, b'V'),
        (15, b'L'),
        (11, b'E'),
        (6, b'T'),
        (6, b'R'),
        (2, b'N'),
        (2, b'I'),
        (1, b'Q'),
    ]);
    ldrf.push(vec![
        (442, b'T'),
        (301, b'K'),
        (283, b'P'),
        (251, b'S'),
        (130, b'A'),
        (114, b'Q'),
        (108, b'R'),
        (57, b'I'),
        (33, b'E'),
        (19, b'N'),
        (14, b'G'),
        (13, b'V'),
        (13, b'H'),
        (11, b'L'),
        (6, b'W'),
        (4, b'Y'),
        (4, b'F'),
        (3, b'C'),
        (2, b'D'),
        (1, b'M'),
    ]);
    ldrf.push(vec![
        (512, b'L'),
        (273, b'C'),
        (249, b'P'),
        (196, b'V'),
        (136, b'A'),
        (119, b'I'),
        (111, b'T'),
        (97, b'F'),
        (48, b'G'),
        (21, b'Y'),
        (20, b'S'),
        (13, b'D'),
        (4, b'W'),
        (4, b'M'),
        (2, b'Q'),
        (1, b'R'),
        (1, b'K'),
        (1, b'H'),
        (1, b'E'),
    ]);
    ldrf.push(vec![
        (359, b'W'),
        (280, b'L'),
        (215, b'I'),
        (215, b'A'),
        (213, b'V'),
        (143, b'H'),
        (105, b'F'),
        (55, b'T'),
        (51, b'S'),
        (42, b'C'),
        (40, b'G'),
        (38, b'Q'),
        (25, b'Y'),
        (10, b'P'),
        (7, b'R'),
        (4, b'M'),
        (4, b'E'),
        (2, b'N'),
        (1, b'K'),
    ]);
    ldrf.push(vec![
        (691, b'A'),
        (458, b'L'),
        (268, b'T'),
        (181, b'S'),
        (57, b'V'),
        (52, b'F'),
        (45, b'I'),
        (16, b'G'),
        (12, b'C'),
        (10, b'M'),
        (5, b'W'),
        (4, b'P'),
        (3, b'R'),
        (2, b'N'),
        (2, b'H'),
        (1, b'Y'),
        (1, b'K'),
        (1, b'E'),
    ]);
    ldrf.push(vec![
        (842, b'L'),
        (645, b'V'),
        (109, b'M'),
        (60, b'F'),
        (52, b'S'),
        (42, b'A'),
        (24, b'I'),
        (24, b'G'),
        (3, b'T'),
        (3, b'C'),
        (2, b'R'),
        (1, b'W'),
        (1, b'Q'),
        (1, b'P'),
    ]);
    ldrf.push(vec![
        (1518, b'L'),
        (120, b'I'),
        (63, b'V'),
        (50, b'F'),
        (15, b'T'),
        (12, b'M'),
        (9, b'G'),
        (8, b'P'),
        (3, b'S'),
        (3, b'Q'),
        (3, b'C'),
        (2, b'R'),
        (1, b'Y'),
        (1, b'W'),
        (1, b'A'),
    ]);

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // BEGIN ANALYSIS OF J REGIONS
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Diagram of J gene structure:
    //
    // left    . . . . . . . . J = between 38 and 75 bases . . . . . .  right
    // flank                           012345678901234567890123456789   flank
    // ------  -------------------------------------------------------  -----
    // ...GTG                          last ten AAs very constrained G  GTAAG
    //                                 this is FWR4                  C    G
    // recomb                 CDR3 tail
    // site                 amino acids
    //                      constrained

    // We tested many primate assemblies for presence of IGHJ4 (numbering as in human), as shown
    // in the table below.  The table suggests that presence of this gene depends on the technology
    // (wetware plus software) used to generate the assemblies, and that for this specific purpose,
    // some technologies are better than others.

    /*
    ---------------------------------------------------------------------------------------------------
    SEQUENCE                                         ACCESSION        SAMPLE         TECHNOLOGY
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GRCh38           human (IGHJ4)  finished BACs
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GCF_008122165.1  gorilla        PB RSII, Ill
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GCA_000306695.2  CHM1           Illumina
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GCA_009914755.1  CHM13          ONT, PB RSII, Ill
    ACTACTTTGAATACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GCA_004024825.1  douc           Ill 2x250 PCR-free
    ACTACTTTGATTACTGGGGCCAGGGGACCCTGGTCACCGTCTCCTCAG GCA_004024785.1  spider monkey  Ill 2x250 PCR-free
    ACTACTTTGATTACTGGGGCCAGGGGACCCTGGTCACCGTCTCCTCAG GCA_004027835.1  howler monkey  Ill 2x250 PCR-free
    ACTACTTTGACTACTGGGGCCAGGGAGTCCTGGTCACCGTCTCCTCAG GCA_004027335.1  Patas monkey   Ill 2x250 PCR-free
    ACTACTTTGAGTACTGGGGCAAAGGGACCCTGGTCACCGTCTCCTCAG GCA_004027275.1  brown lemur    Ill 2x250 PCR-free
    ACTACTTTGAATACTGGGGCCAGGGGACCCTGGTCACCGTCTCCTCAG GCA_004024885.1  tamarin        Ill 2x250 PCR-free
    ACTACTTTGATTACTGGGGCCAGGGGACCCTGGTCACTGTCTCCTCAG GCA_004026645.1  saki           Ill 2x250 PCR-free
    ACTACTTTGAATACTGGGGCCAGGGGACCCAGGTCACCGTCTCCTCAG GCA_004027715.1  titi           Ill 2x250 PCR-free
    ACTACTTTGAATACTGGGGCCAGGGGACCCTGGTCACCGTCTCCTCAG GCA_004027755.1  capuchin       Ill 2x250 PCR-free
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG GCA_000325845.1  chimp          Sanger (2003)
                                           (missing) GCF_002880755.1  chimp          Sanger+kitchen sink
    ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG unpublished      NA12878        Ill 2x250 PCR-free
                                           (missing) GCA_002077035.3  NA12878        PB RSII
                                           (missing) GCA_002022845.1  NA12878        10x
                                           (missing) GCA_001524155.4  NA19240        PB
                                           (missing) GCA_002022975.1  NA19240        10x
    ---------------------------------------------------------------------------------------------------
    */

    // Shown below is the motif for the first five genome bases after a IG J gene, based on human
    // and mouse data.

    let post_j = [
        [b'C', b'G'].to_vec(),
        [b'T'].to_vec(),
        [b'A', b'G'].to_vec(),
        [b'A'].to_vec(),
        [b'G'].to_vec(),
    ];

    // Shown below are motifs for the 3' end of J genes, by chain type.  These are based on
    // human and mouse genes in the 10x reference sequences, after removing the following
    // sequences, which are not present in 10x data (or at least in the large set that was checked):
    //
    // >human-310|IGLJ2 ENST00000390322|IGLJ2|5'UTR|IG|IGL|None|00
    // GTGTGGGGGCCATGTGGACTCCCTCATGAGCAG
    // >human-312|IGLJ3 ENST00000390324|IGLJ3|5'UTR|IG|IGL|None|00
    // GTGTGGGGGCCATGTGGACTCCCTC
    // >human-314|IGLJ4 ENST00000390326|IGLJ4|J-REGION|IG|IGL|None|00
    // GTATTTGGTGGAGGAACCCAGCTGATCATTTTA
    // >human-315|IGLJ5 ENST00000390327|IGLJ5|J-REGION|IG|IGL|None|00
    // CAGAGAGGGTTTTTGTATGAGCCTGTGTCACAGCACTGGGTGTTTGGTGAGGGGACCGAGCTGACCGTCCTA
    // >mouse-342|IGLJ4 ENSMUST00000198313|IGLJ4|J-REGION|IG|IGL|None|00
    // TTGGGTGTTCGGAGGTGGAACCAGATTGACTGTCCTAGATGA
    //
    // The data below simply show all bases that are observed.  In some cases one might continue
    // on, adding more positions.

    let iglj_rev_motif = [
        [b'G'].to_vec(),                   // 1
        [b'A', b'C'].to_vec(),             // 2
        [b'T'].to_vec(),                   // 3
        [b'C'].to_vec(),                   // 4
        [b'C'].to_vec(),                   // 5
        [b'T'].to_vec(),                   // 6
        [b'G'].to_vec(),                   // 7
        [b'C', b'T'].to_vec(),             // 8
        [b'C'].to_vec(),                   // 9
        [b'A'].to_vec(),                   // 10
        [b'C', b'G'].to_vec(),             // 11
        [b'T'].to_vec(),                   // 12
        [b'C', b'G'].to_vec(),             // 13
        [b'A', b'G'].to_vec(),             // 14
        [b'A'].to_vec(),                   // 15
        [b'A', b'C'].to_vec(),             // 16
        [b'C'].to_vec(),                   // 17
        [b'C'].to_vec(),                   // 18
        [b'A'].to_vec(),                   // 19
        [b'A', b'C', b'G'].to_vec(),       // 20
        [b'G'].to_vec(),                   // 21
        [b'G'].to_vec(),                   // 22
        [b'A', b'T'].to_vec(),             // 23
        [b'C', b'G'].to_vec(),             // 24
        [b'A', b'G'].to_vec(),             // 25
        [b'A', b'C', b'T'].to_vec(),       // 26
        [b'G'].to_vec(),                   // 27
        [b'G'].to_vec(),                   // 28
        [b'C'].to_vec(),                   // 29
        [b'T'].to_vec(),                   // 30
        [b'T'].to_vec(),                   // 31
        [b'A', b'C', b'G', b'T'].to_vec(), // 32
        [b'T'].to_vec(),                   // 33
        [b'A', b'G'].to_vec(),             // 34
        [b'T', b'G'].to_vec(),             // 35
    ];

    let igkj_rev_motif = [
        [b'C'].to_vec(),                   // 1
        [b'A'].to_vec(),                   // 2
        [b'A'].to_vec(),                   // 3
        [b'A'].to_vec(),                   // 4
        [b'A', b'C', b'G', b'T'].to_vec(), // 5
        [b'T'].to_vec(),                   // 6
        [b'A', b'C'].to_vec(),             // 7
        [b'A', b'G', b'T'].to_vec(),       // 8
        [b'A'].to_vec(),                   // 9
        [b'G'].to_vec(),                   // 10
        [b'G'].to_vec(),                   // 11
        [b'T'].to_vec(),                   // 12
        [b'C', b'G', b'T'].to_vec(),       // 13
        [b'A', b'G'].to_vec(),             // 14
        [b'A', b'G'].to_vec(),             // 15
        [b'A', b'C'].to_vec(),             // 16
        [b'A', b'C'].to_vec(),             // 17
        [b'C'].to_vec(),                   // 18
        [b'A'].to_vec(),                   // 19
        [b'C', b'G'].to_vec(),             // 20
        [b'G'].to_vec(),                   // 21
        [b'G'].to_vec(),                   // 22
        [b'A', b'G', b'T'].to_vec(),       // 23
        [b'A', b'C', b'G'].to_vec(),       // 24
        [b'C', b'G', b'T'].to_vec(),       // 25
        [b'A', b'C', b'T'].to_vec(),       // 26
        [b'G'].to_vec(),                   // 27
        [b'A', b'G'].to_vec(),             // 28
        [b'C', b'T'].to_vec(),             // 29
        [b'T'].to_vec(),                   // 30
        [b'T'].to_vec(),                   // 31
        [b'A', b'C', b'G', b'T'].to_vec(), // 32
        [b'C', b'G'].to_vec(),             // 33
        [b'A'].to_vec(),                   // 34
        [b'C', b'G'].to_vec(),             // 35
    ];

    let ighj_rev_motif = [
        [b'G'].to_vec(),             // 1
        [b'A'].to_vec(),             // 2
        [b'C'].to_vec(),             // 3
        [b'G', b'T'].to_vec(),       // 4
        [b'C', b'T'].to_vec(),       // 5
        [b'C'].to_vec(),             // 6
        [b'T'].to_vec(),             // 7
        [b'C'].to_vec(),             // 8
        [b'T'].to_vec(),             // 9
        [b'G'].to_vec(),             // 10
        [b'A', b'C', b'T'].to_vec(), // 11
        [b'C'].to_vec(),             // 12
        [b'A'].to_vec(),             // 13
        [b'C'].to_vec(),             // 14
        [b'T'].to_vec(),             // 15
        [b'C', b'G'].to_vec(),       // 16
        [b'A', b'G', b'T'].to_vec(), // 17
        [b'C', b'T'].to_vec(),       // 18
        [b'A', b'C', b'T'].to_vec(), // 19
        [b'A', b'C', b'T'].to_vec(), // 20
        [b'C'].to_vec(),             // 21
        [b'A'].to_vec(),             // 22
        [b'A', b'C', b'G'].to_vec(), // 23
        [b'G'].to_vec(),             // 24
        [b'G'].to_vec(),             // 25
        [b'A', b'G', b'T'].to_vec(), // 26
        [b'A', b'C', b'G'].to_vec(), // 27
        [b'A', b'C'].to_vec(),       // 28
        [b'C', b'T'].to_vec(),       // 29
        [b'G'].to_vec(),             // 30
        [b'G'].to_vec(),             // 31
        [b'G'].to_vec(),             // 32
        [b'G'].to_vec(),             // 33
        [b'T'].to_vec(),             // 34
        [b'C'].to_vec(),             // 35
    ];

    // Define patterns that extend these motifs and are set up for fast computation.

    fn to01(z: &Vec<u8>) -> [u8; 4] {
        let mut x = [0; 4];
        for t in z.iter() {
            if *t == b'A' {
                x[0] = 1;
            } else if *t == b'C' {
                x[1] = 1;
            } else if *t == b'G' {
                x[2] = 1;
            } else {
                x[3] = 1;
            }
        }
        x
    }
    let mut igkj_rev = Vec::<[u8; 4]>::new();
    for i in (0..post_j.len()).rev() {
        igkj_rev.push(to01(&post_j[i]));
    }
    for i in 0..igkj_rev_motif.len() {
        igkj_rev.push(to01(&igkj_rev_motif[i]));
    }
    let mut iglj_rev = Vec::<[u8; 4]>::new();
    for i in (0..post_j.len()).rev() {
        iglj_rev.push(to01(&post_j[i]));
    }
    for i in 0..iglj_rev_motif.len() {
        iglj_rev.push(to01(&iglj_rev_motif[i]));
    }
    let mut ighj_rev = Vec::<[u8; 4]>::new();
    for i in (0..post_j.len()).rev() {
        ighj_rev.push(to01(&post_j[i]));
    }
    for i in 0..ighj_rev_motif.len() {
        ighj_rev.push(to01(&ighj_rev_motif[i]));
    }

    // Heuristics for J gene detection.

    const MAX_J_DIFFS: usize = 2;
    const MIN_J: isize = 38;
    const MAX_J: isize = 75;

    // Search for J region matches.  For mouse and human, by design, 0 works.

    let mut j_pat_len = igkj_rev.len();
    j_pat_len = max(j_pat_len, iglj_rev.len());
    j_pat_len = max(j_pat_len, ighj_rev.len());
    #[derive(Clone)]
    struct Jhit {
        pub tig: usize,
        pub pass: usize,
        pub stop: usize,
        pub rtype: String,
        pub jmatch: Vec<u8>,
        pub seq: Vec<u8>,
        pub lscore1: usize,
        pub lscore2: usize,
        pub errs: usize,
    }
    let mut results = Vec::<(usize, Vec<Jhit>)>::new();
    if use_j {
        for i in 0..refy.len() {
            results.push((i, Vec::<Jhit>::new()));
        }
    }
    let t = Instant::now();
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut r = refy[i].clone();
        let mut rr = vec![0_u8; r.len()];
        for j in 0..r.len() {
            if r[j] == b'A' {
                rr[j] = 0;
            } else if r[j] == b'C' {
                rr[j] = 1;
            } else if r[j] == b'G' {
                rr[j] = 2;
            } else {
                rr[j] = 3;
            }
        }
        for pass in 0..2 {
            if pass == 1 {
                rr.reverse();
                for j in 0..rr.len() {
                    rr[j] = 3 - rr[j];
                }
                reverse_complement(&mut r);
            }
            let js = ["IGKJ", "IGLJ", "IGHJ"];
            for pos in (j_pat_len..=rr.len()).rev() {
                for (u, ig) in [&igkj_rev, &iglj_rev, &ighj_rev].iter().enumerate() {
                    let mut errs = 0;
                    for z in 0..ig.len() {
                        errs += 1 - ig[z][rr[pos - 1 - z] as usize] as usize;
                        if errs > MAX_J_DIFFS {
                            break;
                        }
                    }
                    if errs <= MAX_J_DIFFS {
                        // Look for the 5' end, which is a recombination site.

                        let jstop = pos - post_j.len();
                        let mut top = Vec::<(usize, usize)>::new();
                        for jstart in jstop as isize - MAX_J..=jstop as isize - MIN_J {
                            if jstart < 0 {
                                continue;
                            }
                            let motif;
                            if u == 0 || u == 2 {
                                // IGKJ or IGHJ
                                motif = b"GGTTTTTGT.......................CACTGTG".to_vec();
                            } else {
                                // IGLJ
                                motif = b"GGTTTTTGT............CACTGTG".to_vec();
                            }
                            let mut score = 0;
                            if jstart as usize >= motif.len() {
                                for m in 0..motif.len() {
                                    if motif[motif.len() - m - 1] == r[jstart as usize - 1 - m] {
                                        score += 1;
                                    }
                                }
                            }
                            top.push((score, jstart as usize));
                        }
                        if top.is_empty() {
                            continue;
                        }
                        reverse_sort(&mut top);
                        if !print_fasta && !store_fasta && !print_aa {
                            println!("{}/{}, {}/{}", top[0].0, top[0].1, top[1].0, top[1].1);
                        }
                        let jstart = top[0].1;

                        // Record the J gene.

                        res.1.push(Jhit {
                            tig: i,
                            pass,
                            stop: pos,
                            rtype: js[u].to_string(),
                            jmatch: r[pos - ig.len()..pos].to_vec(),
                            seq: r[jstart..jstop].to_vec(),
                            lscore1: top[0].0,
                            lscore2: top[1].0,
                            errs,
                        });
                    }
                }
            }
        }
    });
    let mut jhits = Vec::<Jhit>::new();
    for i in 0..results.len() {
        jhits.append(&mut results[i].1.clone());
    }

    // Define J region truth data.

    let true_j_human = [
        [
            "49|IGJH1",
            "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
        ],
        // 1. commented out because adding C at the beginning causes a slightly better match on the RSS,
        // and then it equals 51:
        // adding a C at the beginning because it results in a slightly better match on the RSS,
        // this is sort of in a gray area
        // ["50|IGHJ2", "TACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG"],
        [
            "51|IGHJ2",
            "CTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG",
        ],
        // 2. commented out because adding T at the beginning causes a much better match on the RSS,
        // and then it equals 53:
        // ["52|IGHJ3", "GATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG"],
        [
            "53|IGHJ3",
            "TGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG",
        ],
        // 3. commented out because adding AC at the beginning causes a much better match on the RSS,
        // and then it equals 55:
        // ["54|IGHJ4", "TACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"],
        [
            "55|IGHJ4",
            "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        ],
        // 4. commented out because adding AC at the beginning causes a much better match on the RSS,
        // and then it equals 57:
        // ["56|IGHJ5", "AACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"],
        [
            "57|IGHJ5",
            "ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        ],
        // 5. commented out because adding AT at the beginning causes a much better match on the RSS,
        // and then it equals 59:
        // ["58|IGHJ6", "TACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG"],
        [
            "59|IGHJ6",
            "ATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG",
        ],
        ["213|IGKJ1", "GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC"],
        ["214|IGKJ2", "TGTGCAGTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC"],
        ["215|IGKJ3", "ATTCACTTTCGGCCCTGGGACCAAAGTGGATATCAAAC"],
        // 6. added G at beginning because then there's a much better match on the RSS:
        ["216|IGKJ4", "GCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAAC"],
        ["217|IGKJ5", "GATCACCTTCGGCCAAGGGACACGACTGGAGATTAAAC"],
        ["309|IGLJ1", "TTATGTCTTCGGAACTGGGACCAAGGTCACCGTCCTAG"],
        // 7. not found in data so excluded but should check for possible correction:
        // ["310|IGLJ2", "GTGTGGGGGCCATGTGGACTCCCTCATGAGCAG"],
        ["311|IGLJ2", "TGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG"],
        // 8. not found in data so excluded but should check for possible correction:
        // ["312|IGLJ3", "GTGTGGGGGCCATGTGGACTCCCTC"],
        ["313|IGLJ3", "TTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTAG"],
        // 9. not found in data so excluded but should check for possible correction:
        // ["314|IGLJ4", "GTATTTGGTGGAGGAACCCAGCTGATCATTTTA"],
        // 10. not found in data so excluded but should check for possible correction:
        // ["315|IGLJ5", "CAGAGAGGGTTTTTGTATGAGCCTGTGTCACAGCACTGGGTGTTTGGTGAGGGGACCGAGCTGACCGTCCTA"],
        // 11. deleted 32 bases at the beginning because then there a much better match on the RSS:
        ["316|IGLJ6", "TAATGTGTTCGGCAGTGGCACCAAGGTGACCGTCCTCG"],
        // 12. deleted 8 bases at the beginning because then there's a much better match on RSS:
        ["317|IGLJ7", "TGCTGTGTTCGGAGGAGGCACCCAGCTGACCGTCCTCG"],
        // may be OK but excluded from analysis because not present in GRCh38:
        // ["737|IGHJ6", "ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG"],
    ];

    let true_j_mouse = [
        [
            "28|IGHJ1",
            "CTACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG",
        ],
        [
            "29|IGHJ2",
            "ACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACAGTCTCCTCAG",
        ],
        [
            "30|IGHJ3",
            "CCTGGTTTGCTTACTGGGGCCAAGGGACTCTGGTCACTGTCTCTGCAG",
        ],
        [
            "31|IGHJ4",
            "ATTACTATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCAG",
        ],
        ["181|IGKJ1", "GTGGACGTTCGGTGGAGGCACCAAGCTGGAAATCAAAC"],
        ["182|IGKJ2", "TGTACACGTTCGGAGGGGGGACCAAGCTGGAAATAAAAC"],
        ["183|IGKJ3", "AATCACATTCAGTGATGGGACCAGACTGGAAATAAAAC"],
        ["184|IGKJ4", "ATTCACGTTCGGCTCGGGGACAAAGTTGGAAATAAAAC"],
        ["185|IGKJ5", "GCTCACGTTCGGTGCTGGGACCAAGCTGGAGCTGAAAC"],
        ["338|IGLJ1", "CTGGGTGTTCGGTGGAGGAACCAAACTGACTGTCCTAG"],
        ["339|IGLJ2", "TTATGTTTTCGGCGGTGGAACCAAGGTCACTGTCCTAG"],
        ["340|IGLJ3", "GTTTATTTTCGGCAGTGGAACCAAGGTCACTGTCCTAG"],
        // excluded because labeled a pseudogene
        // ["341|IGLJ3PSEUDOGENE", "AGGTTCTTTTTCCTCAAATGGCCTATTGTATGCAGGAG"],
        // excluded because not found in 10x data:
        // ["342|IGLJ4", "TTGGGTGTTCGGAGGTGGAACCAGATTGACTGTCCTAGATGA"],
    ];

    // Analyze J results.

    if use_j {
        if print_fasta || store_fasta || print_aa {
            let (mut h, mut k, mut l) = (0, 0, 0);
            for i in 0..jhits.len() {
                let count;
                if jhits[i].rtype == "IGHJ" {
                    h += 1;
                    count = h;
                } else if jhits[i].rtype == "IGKJ" {
                    k += 1;
                    count = k;
                } else {
                    l += 1;
                    count = l;
                }
                let gene = format!("{}{}", jhits[i].rtype, count);
                let seq = &jhits[i].seq;
                if print_fasta && !tenx {
                    fwriteln!(fasta_log, ">{}{}", jhits[i].rtype, count);
                } else if print_fasta {
                    fasta_count += 1;
                    fwriteln!(
                        fasta_log,
                        ">{}|{} enclone|{}|J-REGION|IG|{}|None|00",
                        fasta_count,
                        gene,
                        gene,
                        strme(&gene.as_bytes()[0..3]),
                    );
                } else if print_aa {
                    println!(">{}", gene);
                    let n = (seq.len() - 1) % 3;
                    println!("{}", strme(&aa_seq(seq, n)));
                }
                if print_fasta {
                    fwriteln!(fasta_log, "{}", strme(seq));
                }
            }
        } else {
            let mut trues = Vec::<Vec<u8>>::new();
            if species == "human" {
                for i in 0..true_j_human.len() {
                    trues.push(true_j_human[i][1].as_bytes().to_vec());
                }
            }
            if species == "mouse" {
                for i in 0..true_j_mouse.len() {
                    trues.push(true_j_mouse[i][1].as_bytes().to_vec());
                }
            }
            let mut ids = Vec::<usize>::new();
            for i in 0..trues.len() {
                ids.push(i);
            }
            sort_sync2(&mut trues, &mut ids);
            let mut found = vec![false; trues.len()];
            println!("\nfound {} J hits", jhits.len());
            for i in 0..jhits.len() {
                let mut x = jhits[i].jmatch.clone();
                x.reverse();
                let or;
                if jhits[i].pass == 0 {
                    or = b'+';
                } else {
                    or = b'-';
                }
                let p = bin_position(&trues, &jhits[i].seq);
                if p >= 0 {
                    found[ids[p as usize]] = true;
                }
                print!(
                    "{} {}{}.{} {} {} {}, {} jseq = {} errs = {}",
                    jhits[i].rtype,
                    or as char,
                    jhits[i].tig,
                    jhits[i].stop,
                    strme(&x[0..5]),
                    strme(&x[5..40]),
                    jhits[i].lscore1,
                    jhits[i].lscore2,
                    strme(&jhits[i].seq),
                    jhits[i].errs
                );
                if p < 0 {
                    print!(" NOT FOUND");
                }
                println!();
            }
            println!();
            if species == "human" {
                for i in 0..true_j_human.len() {
                    if !found[i] {
                        println!("missing {} = {}", true_j_human[i][0], true_j_human[i][1]);
                    }
                }
            } else if species == "mouse" {
                for i in 0..true_j_mouse.len() {
                    if !found[i] {
                        println!("missing {} = {}", true_j_mouse[i][0], true_j_mouse[i][1]);
                    }
                }
            }
            println!("\nused {:.2} seconds\n", elapsed(&t));
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // BEGIN ANALYSIS OF D REGIONS
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // True human D genes.  Mostly not tested for existence in data.
    // Excluded IGHD1-14 because not found in BI=6-12.

    let mut true_d_human = [
        b"CTAACTGGGGA".to_vec(),
        b"TGACTACAGTAACTAC".to_vec(),
        b"TGACTACGGTGACTAC".to_vec(),
        b"GGTACAACTGGAACGAC".to_vec(),
        b"GGTATAACTGGAACGAC".to_vec(),
        b"GGTATAACTGGAACTAC".to_vec(),
        b"GAGTATAGCAGCTCGTCC".to_vec(),
        b"GGGTATAGCAGCGGCTAC".to_vec(),
        b"TGACTACGGTGGTAACTCC".to_vec(),
        b"GGTATAGTGGGAGCTACTAC".to_vec(),
        b"GTAGAGATGGCTACAATTAC".to_vec(),
        b"GTGGATACAGCTATGGTTAC".to_vec(),
        b"GGGTATAGCAGCAGCTGGTAC".to_vec(),
        b"GGGTATAGCAGTGGCTGGTAC".to_vec(),
        b"GTGGATATAGTGGCTACGATTAC".to_vec(),
        b"AGCATATTGTGGTGGTGACTGCTATTCC".to_vec(),
        b"AGGATATTGTACTAATGGTGTATGCTATACC".to_vec(),
        b"AGGATATTGTAGTAGTACCAGCTGCTATGCC".to_vec(),
        b"AGGATATTGTAGTGGTGGTAGCTGCTACTCC".to_vec(),
        b"GTATTACGATATTTTGACTGGTTATTATAAC".to_vec(),
        b"GTATTACGATTTTTGGAGTGGTTATTATACC".to_vec(),
        b"GTATTACTATGATAGTAGTGGTTATTACTAC".to_vec(),
        b"GTATTACTATGGTTCGGGGAGTTATTATAAC".to_vec(),
        b"GTATTATGATTACGTTTGGGGGAGTTATCGTTATACC".to_vec(),
    ]
    .to_vec();
    let mut true_d_human_pseudo = [
        b"GGTATAACTGGAACAAC".to_vec(),
        b"TGACTATGGTGCTAACTAC".to_vec(),
        b"GTGGATATAGTGTCTACGATTAC".to_vec(),
        b"AGAATATTGTAATAGTACTACTTTCTATGCC".to_vec(),
        b"GTATTATGATTTTTGGACTGGTTATTATACC".to_vec(),
    ]
    .to_vec();
    unique_sort(&mut true_d_human);
    unique_sort(&mut true_d_human_pseudo);
    let mut to_delete = vec![false; true_d_human.len()];
    for i in 0..true_d_human.len() {
        if !to_delete[i] {
            let mut x = true_d_human[i].clone();
            reverse_complement(&mut x);
            let j = bin_position(&true_d_human, &x);
            if j >= 0 {
                to_delete[j as usize] = true;
            }
        }
    }
    erase_if(&mut true_d_human, &to_delete);

    // True mouse D genes.
    // Excluded IGHD5-7 and IGHD5-8 because we don't find them in 10x mouse data.

    let mut true_d_mouse = [
        b"GAATACCTAC".to_vec(),
        b"GACTACCTAC".to_vec(),
        b"CTAACTGGGAC".to_vec(),
        b"AGACAGCTCAGGCTAC".to_vec(),
        b"CCTACTATAGTAACTAC".to_vec(),
        b"GGCACAGCTCGGGCTAC".to_vec(),
        b"TCTACTATGATTACGAC".to_vec(),
        b"TCTACTATGGTAACTAC".to_vec(),
        b"TCTACTATGGTTACGAC".to_vec(),
        b"TCTATGATGGTTACTAC".to_vec(),
        b"TTTATTACTACGGTAGTAGCTAC".to_vec(),
    ]
    .to_vec();
    unique_sort(&mut true_d_mouse);
    let mut to_delete = vec![false; true_d_mouse.len()];
    for i in 0..true_d_mouse.len() {
        if !to_delete[i] {
            let mut x = true_d_mouse[i].clone();
            reverse_complement(&mut x);
            let j = bin_position(&true_d_mouse, &x);
            if j >= 0 {
                to_delete[j as usize] = true;
            }
        }
    }
    erase_if(&mut true_d_mouse, &to_delete);
    let true_d_mouse_pseudo = Vec::<Vec<u8>>::new();

    // D region test code:
    // 1. PWM for recombination motif flanking both sides of D from human and mouse, using right
    // version (left is the reverse complement).  Conservation continues much further.
    // 2. Heuristic for match.
    // 3. Bounds on observed lengths of D in human and mouse.

    let dright = [
        [0, 160, 2, 0],   // 1
        [160, 0, 2, 0],   // 2
        [0, 157, 0, 5],   // 3
        [116, 10, 7, 29], // 4
        [6, 1, 150, 5],   // 5
        [0, 10, 0, 152],  // 6
        [21, 12, 129, 0], // 7
        [101, 2, 35, 24], // 8
        [0, 75, 43, 44],  // 9
        [119, 0, 23, 20], // 10
        [26, 101, 26, 9], // 11
        [100, 8, 35, 19], // 12
        [36, 27, 78, 21], // 13; tried to nullify
        [51, 61, 50, 0],  // 14
        [1, 117, 39, 5],  // 15
        [30, 100, 0, 32], // 16
        [0, 0, 0, 0],     // 17; [27, 57, 24, 54]
        [17, 69, 26, 50], // 18
        [60, 19, 39, 44], // 19
        [55, 23, 27, 57], // 20; tried to nullify
        [4, 145, 0, 13],  // 21
        [116, 38, 4, 4],  // 22
        [93, 18, 47, 4],  // 23
        [0, 0, 0, 0],     // 24; [137, 12, 1, 12],
        [160, 2, 0, 0],   // 25
        [124, 2, 2, 34],  // 26
        [3, 129, 18, 12], // 27
        [0, 127, 7, 28],  // 28
    ];
    let mut dright_use = Vec::new();
    for i in 0..dright.len() {
        let mut sum = 0;
        for j in 0..4 {
            sum += dright[i][j];
        }
        dright_use.push(sum > 0);
    }
    pub fn d_right_flank_score(
        x: &[u8],
        dright: &[[usize; 4]],
        dright_use: &Vec<bool>,
        total: usize,
    ) -> f64 {
        const MIN_FRAC: f64 = 0.02;
        let mut p = 1.0;
        for i in 0..dright.len() {
            if dright_use[i] {
                let m = dright[i][x[i] as usize];
                let mut q = m as f64 / total as f64;
                q = q.max(MIN_FRAC);
                p *= q;
            }
        }
        p
    }
    let mut total_d = 0;
    for i in 0..4 {
        total_d += dright[0][i];
    }
    fn to0123(x: &[u8]) -> Vec<u8> {
        let mut y = vec![0_u8; x.len()];
        for i in 0..x.len() {
            if x[i] == b'A' {
                y[i] = 0;
            } else if x[i] == b'C' {
                y[i] = 1;
            } else if x[i] == b'G' {
                y[i] = 2;
            } else {
                y[i] = 3;
            }
        }
        y
    }
    let mut min_p = d_right_flank_score(
        &to0123(b"CACGGCCTGGCACCCCCTGACAATAACCACACCTGGAACT"),
        &dright,
        &dright_use,
        total_d,
    );
    // the following two from IGHD3-16, which is real
    min_p = min_p.min(d_right_flank_score(
        &to0123(b"CACAGCATCACACGGTCCATCAGAAACCCATGCCACAGCC"),
        &dright,
        &dright_use,
        total_d,
    ));
    min_p = min_p.min(d_right_flank_score(
        &to0123(b"CACAGTGACACAGACCTCACTTCAAACCTACCATCTGGCC"),
        &dright,
        &dright_use,
        total_d,
    ));
    min_p *= 0.1;
    const D_LOW: usize = 10;
    const D_HIGH: usize = 37;
    use perf_stats::elapsed;
    use std::time::Instant;
    let t = Instant::now();
    const MAX_TIG: usize = 10_000_000;
    const TIG_OVERLAP: usize = 100_000;
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    if use_d {
        for i in 0..refy.len() {
            results.push((i, Vec::<Vec<u8>>::new()));
        }
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut r = refy[i].clone();
        let mut rr = vec![0_u8; r.len()];
        for j in 0..r.len() {
            if r[j] == b'A' {
                rr[j] = 0;
            } else if r[j] == b'C' {
                rr[j] = 1;
            } else if r[j] == b'G' {
                rr[j] = 2;
            } else {
                rr[j] = 3;
            }
        }
        let mut matches = Vec::<(usize, usize)>::new();
        for pass in 0..2 {
            if pass == 1 {
                rr.reverse();
                for j in 0..rr.len() {
                    rr[j] = 3 - rr[j];
                }
                reverse_complement(&mut r);
            }
            if rr.len() < dright.len() {
                continue;
            }
            for j in 0..=rr.len() - dright.len() {
                if d_right_flank_score(&rr[j..], &dright, &dright_use, total_d) >= min_p {
                    let pos;
                    if pass == 0 {
                        pos = j;
                    } else {
                        pos = rr.len() - j;
                    }
                    matches.push((pos, 1 - pass));
                }
            }
        }
        rr.reverse();
        for j in 0..rr.len() {
            rr[j] = 3 - rr[j];
        }
        reverse_complement(&mut r);
        matches.sort_unstable();
        for z in 0..matches.len() {
            // not sure if we're forcing orientation here
            if matches[z].1 == 0 {
                for j in z + 1..matches.len() {
                    let start = matches[z].0;
                    let stop = matches[j].0;
                    let dlen = stop - start;
                    if dlen > D_HIGH + 1 {
                        break;
                    }
                    if matches[j].1 == 1 && dlen >= D_LOW - 1 {
                        if !print_fasta && !store_fasta {
                            let mut s = r[start..stop].to_vec();
                            reverse_complement(&mut s);
                            print!(
                                "{} = {}', {}.{}-{} of {}",
                                strme(&r[start..stop]),
                                strme(&s),
                                i,
                                start,
                                stop,
                                r.len()
                            );
                        }
                        let mut x = r[start..stop].to_vec();
                        if !print_fasta && !store_fasta {
                            let mut good = false;
                            if species == "human" {
                                good = bin_member(&true_d_human, &x);
                            } else if species == "mouse" {
                                good = bin_member(&true_d_mouse, &x);
                            }
                            if !good {
                                reverse_complement(&mut x);
                                if species == "human" {
                                    good = bin_member(&true_d_human, &x);
                                } else if species == "mouse" {
                                    good = bin_member(&true_d_mouse, &x);
                                }
                            }
                            if good {
                                print!(" GOOD");
                            }
                            println!();
                        }
                        res.1.push(r[start..stop].to_vec());
                    }
                }
            }
        }
    });
    if use_d {
        if print_fasta || store_fasta {
            let mut count = 0;
            for i in 0..results.len() {
                for j in 0..results[i].1.len() {
                    let mut x = results[i].1[j].clone();
                    count += 1;
                    if !tenx {
                        fwriteln!(fasta_log, ">IGHD{}_FW", count);
                    } else {
                        fasta_count += 1;
                        let gene = format!("IGHD{}_FW", count);
                        fwriteln!(
                            fasta_log,
                            ">{}|{} enclone|{}|D-REGION|IG|{}|None|00",
                            fasta_count,
                            gene,
                            gene,
                            strme(&gene.as_bytes()[0..3]),
                        );
                    }
                    fwriteln!(fasta_log, "{}", strme(&x));
                    reverse_complement(&mut x);
                    if !tenx {
                        fwriteln!(fasta_log, ">IGHD{}_RC", count);
                    } else {
                        fasta_count += 1;
                        let gene = format!("IGHD{}_RC", count);
                        fwriteln!(
                            fasta_log,
                            ">{}|{} enclone|{}|D-REGION|IG|{}|None|00",
                            fasta_count,
                            gene,
                            gene,
                            strme(&gene.as_bytes()[0..3]),
                        );
                    }
                    fwriteln!(fasta_log, "{}", strme(&x));
                }
            }
        } else {
            let mut all = Vec::<Vec<u8>>::new();
            for i in 0..results.len() {
                all.append(&mut results[i].1.clone());
            }
            unique_sort(&mut all);
            let mut to_delete = vec![false; all.len()];
            for i in 0..all.len() {
                if !to_delete[i] {
                    let mut x = all[i].clone();
                    reverse_complement(&mut x);
                    let j = bin_position(&all, &x);
                    if j >= 0 {
                        to_delete[j as usize] = true;
                    }
                }
            }
            erase_if(&mut all, &to_delete);
            println!();
            let mut hits = 0;
            let mut phits = 0;
            let mut true_d = Vec::<Vec<u8>>::new();
            let mut true_d_pseudo = Vec::<Vec<u8>>::new();
            if species == "human" {
                true_d = true_d_human.clone();
                true_d_pseudo = true_d_human_pseudo.clone();
            } else if species == "mouse" {
                true_d = true_d_mouse.clone();
                true_d_pseudo = true_d_mouse_pseudo;
            }
            for i in 0..true_d.len() {
                let mut x = true_d[i].clone();
                if bin_member(&all, &x) {
                    hits += 1;
                } else {
                    reverse_complement(&mut x);
                    if bin_member(&all, &x) {
                        hits += 1;
                    } else {
                        println!("missed {}", strme(&true_d[i]));
                    }
                }
            }
            for i in 0..true_d_pseudo.len() {
                let mut x = true_d_pseudo[i].clone();
                if bin_member(&all, &x) {
                    phits += 1;
                } else {
                    reverse_complement(&mut x);
                    if bin_member(&all, &x) {
                        phits += 1;
                    } else {
                        println!("missed pseudo {}", strme(&true_d_human_pseudo[i]));
                    }
                }
            }
            println!(
                "found {}; {} of {} trues; {} of {} pseudos",
                all.len(),
                hits,
                true_d.len(),
                phits,
                true_d_pseudo.len()
            );
            println!("used {:.2} seconds\n", elapsed(&t));
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Get mammalian fixed len data.

    let mmm = mammalian_fixed_len();
    let mut mstart = HashMap::<(String, String), usize>::new();
    let mut mstop = HashMap::<(String, String), usize>::new();
    for chain in ["IGH", "IGK", "IGL", "TRA", "TRB"].iter() {
        for feature in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"].iter() {
            let low = mmm.lower_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            let high = mmm.upper_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            mstart.insert((chain.to_string(), feature.to_string()), low);
            mstop.insert((chain.to_string(), feature.to_string()), high);
        }
    }

    // Get mammalian pwm data.

    let mmp = mammalian_pwms();
    let mut mpstart = HashMap::<(String, String), usize>::new();
    let mut mpstop = HashMap::<(String, String), usize>::new();
    for chain in ["IGH", "IGK", "IGL", "TRA", "TRB"].iter() {
        for feature in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"].iter() {
            let low = mmp.lower_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            let high = mmp.upper_bound_by_key(&(*chain, *feature), |(a, b, _c, _d)| (a, b));
            mpstart.insert((chain.to_string(), feature.to_string()), low);
            mpstop.insert((chain.to_string(), feature.to_string()), high);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // LOOK FOR V GENES
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // We parallelize over refy, which has overlapping segments.  Haven't checked to see if
    // this might result in artifacts.

    let mut results = Vec::<(
        usize,                  // 0: index
        Vec<String>,            // 1: log
        Vec<String>,            // 2: defunct
        Vec<(String, Vec<u8>)>, // 3: special for pwm option
        Vec<Vec<u8>>,           // 4: upstream tag (synced to log)
        Vec<String>,            // 5: reference matches (synced to log)
        Vec<Vec<Vec<u8>>>,      // 6: bases through end of FWR3 (synced to log)
        Vec<Vec<String>>,       // 7: chain type (synced to log)
        Vec<Vec<Vec<u8>>>,      // 8: utr bases tag (synced to log)
    )>::new();
    if use_v {
        for i in 0..2 * refy.len() {
            results.push((
                i,
                Vec::<String>::new(),
                Vec::<String>::new(),
                Vec::<(String, Vec<u8>)>::new(),
                Vec::<Vec<u8>>::new(),
                Vec::<String>::new(),
                Vec::<Vec<Vec<u8>>>::new(),
                Vec::<Vec<String>>::new(),
                Vec::<Vec<Vec<u8>>>::new(),
            ));
        }
    }
    results.par_iter_mut().for_each(|res| {
        let p = res.0;
        let i = p / 2;
        let pass = p % 2;
        let mut r = refy[i].clone();
        if pass == 1 {
            reverse_complement(&mut r);
        }
        for j in 0..r.len() {
            if r[j] != b'A' && r[j] != b'C' && r[j] != b'G' && r[j] != b'T' {
                r[j] = b'A';
            }
        }
        for frame in 0..3 {
            let aa = aa_seq(&r, frame);
            let mut j = 0;
            while j < aa.len() {
                if aa[j] == b'*' {
                    j += 1;
                    continue;
                }
                let mut k = j;
                let mut good = false;
                while k < aa.len() && aa[k] != b'*' {
                    if aa[k..].starts_with(b"YYC")
                        || aa[k..].starts_with(b"YFC")
                        || aa[k..].starts_with(b"YHC")
                        || aa[k..].starts_with(b"FYC")
                        || aa[k..].starts_with(b"YIC")
                    {
                        good = true;
                    }
                    k += 1;
                }

                if !good {
                    j = k;
                    continue;
                }

                // Find stopping points k2 for a putative FWR3 region.

                let mut k2s = Vec::<usize>::new();
                let mut k2 = k;
                while k2 >= j + 3 {
                    if aa[j..k2].ends_with(b"YYC")
                        || aa[j..k2].ends_with(b"YFC")
                        || aa[j..k2].ends_with(b"YHC")
                        || aa[j..k2].ends_with(b"FYC")
                        || aa[j..k2].ends_with(b"YIC")
                    {
                        k2s.push(k2);
                    }
                    k2 -= 1;
                }
                if k2s.is_empty() {
                    j = k;
                    continue;
                }

                // Pick best stopping point.

                let mut max_score = 0.0;
                let mut best_k2 = k2s[0];
                for z in 0..k2s.len() {
                    let k2 = k2s[z];
                    let x = &aa[j..k2];
                    let score = score_fwr3_at_end(x, 0, &freqs);
                    if score > max_score {
                        max_score = score;
                        best_k2 = k2;
                    }
                }
                let k2 = best_k2;
                let score = max_score;

                // Test for RAG site.

                const MAX_CDR3_PART: usize = 44;
                let first = 3 * k2 + frame;
                const UPSTREAM: usize = TAG_LEN + MAX_CDR3_PART;
                if first + UPSTREAM > r.len() {
                    j = k;
                    continue;
                }
                let upstream = r[first..first + UPSTREAM].to_vec();

                // Proceed.

                // println!("ORF = {}, score = {}", strme(&aa[j..k2]), score);
                if score < MIN_SCORE {
                    j = k;
                    continue;
                }

                // Set up to gather info: a leader sequence plus a message plus a bunch of stuff.

                let mut info = Vec::<(
                    Vec<u8>,
                    String,
                    u8,
                    usize,
                    Vec<u8>,
                    usize,
                    usize,
                    f64,
                    u8,
                    Vec<u8>,
                    Vec<u8>,
                    String,
                    Vec<u8>,
                )>::new();

                // Look for possible exon 1 sequences.

                let mut scount = 0;
                let mut log = String::new();
                let mut mlp = Vec::<(isize, usize, u8)>::new();
                for m in -175..WALK_BACK as isize {
                    // start = start of exon 1 in bases
                    if m > 3 * j as isize + frame as isize {
                        continue;
                    }
                    let start = (3 * j as isize + frame as isize - m) as usize;
                    if start >= r.len() {
                        continue;
                    }
                    if !r[start..].starts_with(b"ATG") {
                        continue;
                    }
                    // l = length of exon 1 in bases (from start codon)
                    for l in (LOW1..=HIGH1).step_by(3) {
                        if start + l + 1 >= r.len() {
                            break;
                        }
                        let bb = aa_seq(&r[start..start + l].to_vec(), 0);
                        if bb.contains(&b'*') {
                            continue;
                        }

                        // Check that the sequence after exon 1 is GT or
                        // GC.  In human and mouse we observe GT, and rarely GC.

                        if r[start + l] != b'G' {
                            continue;
                        }
                        if r[start + l + 1] != b'T' && r[start + l + 1] != b'C' {
                            continue;
                        }

                        // Save.

                        mlp.push((m, l, r[start + l + 1]));
                    }
                }

                // Look for possible exon 2 sequences.

                for pp in 0..mlp.len() {
                    let m = mlp[pp].0;
                    let l = mlp[pp].1; // length of exon 1 in bases
                    let post = mlp[pp].2; // T and rarely C
                                          // start = start of exon 1 in bases
                    let start = (3 * j as isize + frame as isize - m) as usize;
                    if start < 10 {
                        continue;
                    }
                    let bb = aa_seq(&r[start..start + l].to_vec(), 0);
                    for mx in LOW_INTRON..=HIGH_INTRON {
                        let start2 = start + l + mx;
                        if start2 > r.len() || r.len() - start2 < 10 {
                            continue;
                        }
                        let e2ss = exon2_start_score(&r[start2..]);
                        if e2ss < -9.0 {
                            continue;
                        }

                        // In human and mouse, we observe that the sequence preceeding exon 2 is
                        // CAG, with TAG and AAG as rare exceptions (and very rare other cases).

                        if !r[start2 - 3..].starts_with(b"CAG")
                            && !r[start2 - 3..].starts_with(b"TAG")
                            && !r[start2 - 3..].starts_with(b"AAG")
                        {
                            continue;
                        }
                        let stop2 = 3 * k2 + frame;
                        if start2 + 2 >= stop2 {
                            continue;
                        }

                        // Don't allow last base of exon 1 to be C, because it does not happen
                        // in the human and mouse reference sequences.

                        if r[start + l - 1] == b'C' {
                            continue;
                        }

                        // Keep going.

                        let cc = aa_seq(&r[start2 + 2..stop2].to_vec(), 0);
                        let mid = vec![r[start + l - 1], r[start2], r[start2 + 1]];
                        let mut s = String::new();
                        if show_transition {
                            s = format!(
                                "{}({}) ..{}.. <{}>({}){}",
                                strme(&bb),
                                r[start + l - 1] as char,
                                codon_to_aa(&mid) as char,
                                strme(&r[start2 - 10..start2]),
                                strme(&r[start2..start2 + 2]),
                                strme(&cc),
                            );
                        }
                        let mut full = bb.clone();
                        full.push(codon_to_aa(&mid));
                        full.append(&mut cc.clone());
                        if full.len() < MIN_LEN || full.len() > MAX_LEN || full.contains(&b'*') {
                            continue;
                        }
                        #[derive(Clone, PartialOrd, PartialEq)]
                        struct ChainData {
                            pub score: f64,
                            pub flen: usize,
                            pub f1_start: usize,
                            pub ldrx: Vec<u8>,
                            pub ct: String,
                            pub start_score: f64,
                            pub s: String,
                        }
                        let mut chain_datas = Vec::<ChainData>::new();
                        for ct in ["IGH", "IGK", "IGL"].iter() {
                            let f1_start = fr1_start(&full, ct);

                            // Test leader.

                            let ldr = &full[0..f1_start];
                            if ldr.len() < MIN_LEADER || ldr.len() > MAX_LEADER {
                                continue;
                            }
                            if lscore1(ldr) < 4 {
                                continue;
                            }
                            if post == b'C' && lscore1(ldr) < 6 {
                                continue;
                            }
                            if *ct == "IGL" && ldr.len() < 17 {
                                continue;
                            }

                            // If the leader has length 20 amino acids, and is preceeded by amino
                            // acids MX, where X is any amino acid, punt, because the longer
                            // (22 amino acid) form is the correct one.  Ditto for a MXY.
                            // And just M.

                            fn is_stop(x: &[u8], j: usize) -> bool {
                                x[j..j + 3].to_vec() == b"TAA".to_vec()
                                    || x[j..j + 3].to_vec() == b"TAG".to_vec()
                                    || x[j..j + 3].to_vec() == b"TGA".to_vec()
                            }
                            if f1_start == 19
                                && start >= 9
                                && r[start - 9..start - 6].to_vec() == b"ATG".to_vec()
                                && !is_stop(&r, start - 6)
                                && !is_stop(&r, start - 3)
                            {
                                continue;
                            }
                            if f1_start == 19
                                && start >= 3
                                && r[start - 3..start].to_vec() == b"ATG".to_vec()
                            {
                                continue;
                            }

                            if f1_start == 26 {
                                let mut m = false;
                                for k in 1..=7 {
                                    if full[k] == b'M' {
                                        m = true;
                                    }
                                }
                                if !m {
                                    continue;
                                }
                            }

                            // Keep going.

                            let flen = full.len() - f1_start;
                            let f1 = fwr1(&full, ct, false);
                            if f1.is_some() {
                                let f1 = f1.unwrap();
                                let c1 = cdr1(&full, ct, false);
                                if f1.len() >= MIN_FWR1
                                    && f1.len() <= MAX_FWR1
                                    && c1.is_some()
                                    && fwr2(&full, ct, false).is_some()
                                    && cdr2(&full, ct, false).is_some()
                                    && fwr3(&full, ct, false).is_some()
                                {
                                    let f2 = fwr2(&full, ct, false).unwrap();
                                    let f3 = fwr3(&full, ct, false).unwrap();
                                    if *ct == "IGL" && f1.len() < 20 {
                                        continue;
                                    }
                                    if *ct == "IGL" && f1.len() > 22 {
                                        continue;
                                    }
                                    let c1 = cdr1(&full, ct, false).unwrap();
                                    let c2 = cdr2(&full, ct, false).unwrap();
                                    let mut score1 = 0.0;
                                    for my in 0..freq1.len() {
                                        let mut total = 0;
                                        let mut hit = 0;
                                        for p in 0..freq1[my].len() {
                                            total += freq1[my][p].0;
                                            if freq1[my][p].1 == f1[my] {
                                                hit = freq1[my][p].0;
                                            }
                                        }
                                        score1 += hit as f64 / total as f64;
                                    }

                                    if score1 >= 1.8 {
                                        // Compute peer group negative log probability two.
                                        // We exclude cdr1 and cdr2 because for human IGLV9-49,
                                        // we observe that the CDR2 sequence appears totally novel,
                                        // relative to other CDR2 sequences.  It also has length
                                        // 12, which is one longer than expected.  The gene
                                        // is probably functional because we seem to observe
                                        // class switching and SHM in our data.

                                        let mut nlogp2 = 0.0;
                                        for f in ["fwr1", "fwr2", "fwr3"].iter() {
                                            let x1 = mpstart[&(ct.to_string(), f.to_string())];
                                            let x2 = mpstop[&(ct.to_string(), f.to_string())];
                                            let z;
                                            if *f == "fwr1" {
                                                z = &f1;
                                            } else if *f == "fwr2" {
                                                z = &f2;
                                            } else if *f == "fwr3" {
                                                z = &f3;
                                            } else if *f == "cdr1" {
                                                z = &c1;
                                            } else {
                                                z = &c2;
                                            }
                                            let len = z.len();
                                            if mmp[x2 - 1].2 < len {
                                                nlogp2 += 20.0;
                                                continue;
                                            }
                                            for x in x1..x2 {
                                                if mmp[x].2 >= len {
                                                    let x1 = &z;
                                                    let pwm = &mmp[x].3;

                                                    // Define scoring scheme.  This is very ugly,
                                                    // because the second argument b is treated as
                                                    // the column number in the position weight
                                                    // matrix.

                                                    let mut n = 0_i32;
                                                    for i in 0..pwm[0].len() {
                                                        n += pwm[0][i].0 as i32;
                                                    }
                                                    let score = |a: u8, b: u8| {
                                                        let mut dot_count = 0;
                                                        for x in pwm[b as usize].iter() {
                                                            if x.1 == b'.' {
                                                                dot_count = x.0 as i32;
                                                            }
                                                        }
                                                        for x in pwm[b as usize].iter() {
                                                            if a == x.1 {
                                                                return x.0 as i32 - dot_count;
                                                            }
                                                        }
                                                        -n
                                                    };
                                                    let (gap_open, gap_extend) =
                                                        (-1 * n as i32, -n as i32);

                                                    // Set up the aligner.

                                                    let mut x2 = Vec::<u8>::new();
                                                    for i in 0..pwm.len() {
                                                        x2.push(i as u8);
                                                    }
                                                    let (n1, n2) =
                                                        (x1.len() as u32, x2.len() as u32);
                                                    assert!(n1 <= n2);
                                                    let l = max(1, (n2 - n1) / 2);
                                                    let corners =
                                                        vec![(0, l), (n1 - 1, n1 - 1 + l)];
                                                    let path = vec![0, 1];
                                                    let mut aligner = Aligner::new(
                                                        gap_open, gap_extend, &score, 1, l as usize,
                                                    );

                                                    // Align.

                                                    let align = aligner.custom_with_match_path(
                                                        x1, &x2, &corners, &path,
                                                    );
                                                    let ops = &align.operations;

                                                    // Evaluate.  It's not totally clear that
                                                    // the scoring of insertions and deletions
                                                    // makes sense.

                                                    let (mut pwm_pos, mut tig_pos) = (0, 0);
                                                    for q in 0..ops.len() {
                                                        if ops[q] == Ins {
                                                            let hit = 1.0;
                                                            let mut total = 0;
                                                            if pwm_pos >= mmp[x].3.len() {
                                                                continue;
                                                            }
                                                            for z in mmp[x].3[pwm_pos].iter() {
                                                                total += z.0;
                                                            }
                                                            nlogp2 -= (hit / total as f64).log10();
                                                            tig_pos += 1;
                                                        } else if ops[q] == Del {
                                                            let mut hit = 0;
                                                            let mut total = 0;
                                                            for r in mmp[x].3[pwm_pos].iter() {
                                                                if r.1 == b'.' {
                                                                    hit += r.0;
                                                                }
                                                                total += r.0;
                                                            }
                                                            hit = max(1, hit);
                                                            nlogp2 -=
                                                                (hit as f64 / total as f64).log10();
                                                            pwm_pos += 1;
                                                        } else {
                                                            let mut hit = 0;
                                                            let mut total = 0;
                                                            for r in mmp[x].3[pwm_pos].iter() {
                                                                if r.1 == z[tig_pos] {
                                                                    hit += r.0;
                                                                }
                                                                total += r.0;
                                                            }
                                                            hit = max(1, hit);
                                                            nlogp2 -=
                                                                (hit as f64 / total as f64).log10();
                                                            tig_pos += 1;
                                                            pwm_pos += 1;
                                                        }
                                                    }
                                                    break;
                                                }
                                            }
                                        }

                                        // Compute reverse leader score.

                                        let mut rlscore = 0.0;
                                        for my in 0..ldrf.len() {
                                            let mut total = 0;
                                            let mut hit = 0;
                                            for p in 0..ldrf[my].len() {
                                                total += ldrf[my][p].0;
                                                if ldrf[my][p].1 == ldr[ldr.len() - my - 1] {
                                                    hit = ldrf[my][p].0;
                                                }
                                            }
                                            rlscore += hit as f64 / total as f64;
                                        }
                                        let mut start_score = 1.0;
                                        for m in 0..23 {
                                            let pos = start - 10 + m;
                                            let total: usize = start_pwm[m].iter().sum();
                                            let mut hit = 0;
                                            if r[pos] == b'A' {
                                                hit += start_pwm[m][0];
                                            } else if r[pos] == b'C' {
                                                hit += start_pwm[m][1];
                                            } else if r[pos] == b'G' {
                                                hit += start_pwm[m][2];
                                            } else {
                                                hit += start_pwm[m][3];
                                            }
                                            start_score *= hit as f64 / total as f64;
                                        }
                                        let muscore = 4.0 * rlscore + start_score.log10();
                                        if muscore < -12.0 {
                                            continue;
                                        }
                                        if nlogp2 >= 58.0 {
                                            continue;
                                        }
                                        let mut s = String::new();
                                        if !show_transition {
                                            s = format!(
                                                "%{}% {} 🌸 {}, rev_ldr_score={:.2}, \
                                                start_score={:.2} \
                                                muscore={:.2} e2ss={:.2} nlogp2={:.2}, ct = {}",
                                                post as char,
                                                strme(&full[0..f1_start]),
                                                strme(&full[f1_start..]),
                                                rlscore,
                                                start_score.log10(),
                                                muscore,
                                                e2ss,
                                                nlogp2,
                                                ct
                                            );
                                        }
                                        chain_datas.push(ChainData {
                                            score: nlogp2,
                                            start_score,
                                            flen,
                                            f1_start,
                                            ldrx: ldr.to_vec(),
                                            ct: ct.to_string(),
                                            s,
                                        });
                                    }
                                }
                            }
                        }
                        if chain_datas.is_empty() {
                            continue;
                        }
                        chain_datas.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        const MAX_SCORE_DIFF: f64 = 0.0; // ONLY KEEPING TOP ENTRY!
                        let mut to_delete = vec![false; chain_datas.len()];
                        for m in 1..chain_datas.len() {
                            if chain_datas[m].score >= chain_datas[0].score + MAX_SCORE_DIFF {
                                to_delete[m] = true;
                            }
                        }
                        erase_if(&mut chain_datas, &to_delete);
                        let max_score = chain_datas.last().unwrap().score;
                        let start_score = chain_datas[0].start_score; // not quite right
                        for m in 0..chain_datas.len() {
                            s += &chain_datas[m].s.clone();
                        }
                        let ldrx = chain_datas[0].ldrx.clone();
                        if max_score >= 40.0 {
                            // flag high scores
                            s += " HIGH";
                        }
                        let mut flens = Vec::<usize>::new();
                        let mut f1_starts = Vec::<usize>::new();
                        let mut cts = Vec::<String>::new();
                        for m in 0..chain_datas.len() {
                            flens.push(chain_datas[m].flen);
                            f1_starts.push(chain_datas[m].f1_start);
                            cts.push(chain_datas[m].ct.clone());
                        }
                        unique_sort(&mut flens);
                        unique_sort(&mut f1_starts);
                        let mut msg = format!(" ({})", flens.iter().format(","));
                        const MAX_MIS: usize = 3;
                        let mut is_match = false;
                        if to_ref.contains_key(&full) {
                            msg += &format!(" 🔴 {}", to_ref[&full].clone());
                            let true_chain = &to_ref[&full].after("|")[0..3];
                            if !s.contains(&format!("ct = {}", true_chain)) {
                                msg += " WRONG CHAIN";
                            }
                            if pwm {
                                let mut x = Vec::<u8>::new();
                                x.append(&mut r[start - 10..start + 13].to_vec());
                                res.3.push((to_ref[&full].clone(), x));
                            }
                        } else {
                            for x in to_ref.iter() {
                                if x.0.len() == full.len() {
                                    let mut diffs = 0;
                                    for z in 0..full.len() {
                                        if x.0[z] != full[z] {
                                            diffs += 1;
                                        }
                                    }
                                    if diffs <= MAX_MIS {
                                        is_match = true;
                                    }
                                }
                            }
                        }
                        scount += 1;
                        let mut logp = format!(
                            "[{}] {} ==> ({}) {}{}",
                            scount,
                            mx,
                            f1_starts.iter().format(","),
                            s,
                            msg
                        );
                        logp += &format!(" {}", cts.iter().format(","));
                        let mut pseudo = false;
                        for p in pseudos.iter() {
                            if p.2.len() == full.len() {
                                let mut diffs = 0;
                                for z in 0..full.len() {
                                    if p.2[z] != full[z] {
                                        diffs += 1;
                                    }
                                }
                                if diffs <= MAX_MIS {
                                    pseudo = true;
                                }
                            }
                        }
                        if is_match {
                            logp += " MATCH";
                        } else if pseudo {
                            logp += " PSEUDO";
                        }
                        let start1 = start;
                        let stop1 = start + l;
                        let mut bases = r[start1..stop1].to_vec();
                        bases.append(&mut r[start2..stop2].to_vec());
                        if start1 < 80 {
                            continue;
                        }
                        let utr_tag = r[start1 - 80..start1].to_vec();
                        if print_bases {
                            logp += " ";
                            logp += strme(&bases);
                        }
                        logp += "\n";
                        info.push((
                            ldrx,
                            logp,
                            post,
                            mx,
                            full.clone(),
                            start2,
                            start,
                            start_score,
                            r[start + l - 1],
                            upstream.clone(),
                            bases,
                            cts[0].clone(),
                            utr_tag,
                        ));
                    }
                }

                // Competitively screen, part 1.  We do not apply the screening to "subset"
                // cases M...M..., although possibly we should.  The problem is that we don't
                // know if the short form or the long form is correct.  Conservation would suggest
                // that the short form is likely to be correct (as it would match the expected
                // length), but we just don't know.  Perhaps CHIP-seq data could settle this

                let mut to_delete = vec![false; info.len()];
                for i1 in 0..info.len() {
                    for i2 in 0..info.len() {
                        let ldr1 = &info[i1].0;
                        let ldr2 = &info[i2].0;
                        if !ldr2.ends_with(ldr1) {
                            let post1 = info[i1].2;
                            let post2 = info[i2].2;
                            if lscore2(ldr1) >= 10 * lscore2(ldr2)
                                && lscore1(ldr1) >= lscore1(ldr2)
                                && post1 != b'C'
                            {
                                to_delete[i2] = true;
                            }
                            if lscore2(ldr1) >= lscore2(ldr2)
                                && lscore1(ldr1) >= lscore1(ldr2)
                                && post1 == b'T'
                                && post2 == b'C'
                            {
                                to_delete[i2] = true;
                            }
                        }
                    }
                }
                erase_if(&mut info, &to_delete);

                // If we've gotten this far, favor a much shorter intron.

                let mut to_delete = vec![false; info.len()];
                for i1 in 0..info.len() {
                    for i2 in 0..info.len() {
                        let ldr1 = &info[i1].0;
                        let ldr2 = &info[i2].0;
                        let start_score1 = info[i1].7;
                        let start_score2 = info[i2].7;
                        if info[i2].3 >= 150 + info[i1].3
                            && lscore2(ldr1) * 10 > lscore2(ldr2)
                            && start_score1 >= start_score2
                        {
                            to_delete[i2] = true;
                        }
                    }
                }
                erase_if(&mut info, &to_delete);

                // Having a proline at reverse position 5 or 6 is better than not having one there.

                let mut to_delete = vec![false; info.len()];
                for i1 in 0..info.len() {
                    for i2 in 0..info.len() {
                        let ldr1 = &info[i1].0;
                        let ldr2 = &info[i2].0;
                        let p1 = ldr1[ldr1.len() - 5] == b'P' || ldr1[ldr1.len() - 6] == b'P';
                        let p2 = ldr2[ldr2.len() - 5] == b'P' || ldr2[ldr2.len() - 6] == b'P';
                        if p1 && !p2 && lscore2(ldr2) < 100 * lscore2(ldr1) {
                            to_delete[i2] = true;
                        }
                    }
                }
                erase_if(&mut info, &to_delete);

                // Favor higher reverse leader score.  Ugly code duplication.

                let mut to_delete = vec![false; info.len()];
                for i1 in 0..info.len() {
                    for i2 in 0..info.len() {
                        let ldr1 = &info[i1].0;
                        let ldr2 = &info[i2].0;
                        let mut rlscore1 = 0.0;
                        for my in 0..ldrf.len() {
                            let mut total = 0;
                            let mut hit = 0;
                            for p in 0..ldrf[my].len() {
                                total += ldrf[my][p].0;
                                if ldrf[my][p].1 == ldr1[ldr1.len() - my - 1] {
                                    hit = ldrf[my][p].0;
                                }
                            }
                            rlscore1 += hit as f64 / total as f64;
                        }
                        let mut rlscore2 = 0.0;
                        for my in 0..ldrf.len() {
                            let mut total = 0;
                            let mut hit = 0;
                            for p in 0..ldrf[my].len() {
                                total += ldrf[my][p].0;
                                if ldrf[my][p].1 == ldr2[ldr2.len() - my - 1] {
                                    hit = ldrf[my][p].0;
                                }
                            }
                            rlscore2 += hit as f64 / total as f64;
                        }
                        if rlscore1 > rlscore2 {
                            // Test for better start score.
                            if info[i1].7 > info[i2].7 {
                                to_delete[i2] = true;
                            }
                        }
                    }
                }
                erase_if(&mut info, &to_delete);

                // If there are two exon 1 sequences having the same start, and one ends in G,
                // that's better.

                let mut to_delete = vec![false; info.len()];
                for i1 in 0..info.len() {
                    for i2 in 0..info.len() {
                        if info[i1].6 == info[i2].6 && info[i1].8 == b'G' && info[i2].8 != b'G' {
                            to_delete[i2] = true;
                        }
                    }
                }
                erase_if(&mut info, &to_delete);

                // Finally, we use the reference to disambiguate in cases where the only difference
                // is how far left the leader goes.  See discussion at beginning.  If
                // there is no reference hit, pick based on leader length frequency.

                if !amb {
                    let mut to_delete = vec![false; info.len()];
                    for i1 in 0..info.len() {
                        let full1 = &info[i1].4;
                        let ref1 = to_ref.contains_key(full1);
                        for i2 in 0..info.len() {
                            let full2 = &info[i2].4;
                            let ref2 = to_ref.contains_key(full2);
                            let ldr1 = &info[i1].0;
                            let ldr2 = &info[i2].0;
                            if (full1.ends_with(full2) || full2.ends_with(full1))
                                && ((ref1 && !ref2)
                                    || (!ref1 && !ref2 && lscore2(ldr1) > lscore2(ldr2)))
                            {
                                to_delete[i2] = true;
                            }
                        }
                    }
                    erase_if(&mut info, &to_delete);
                    let mut to_delete = vec![false; info.len()];
                    for i1 in 0..info.len() {
                        let full1 = &info[i1].4;
                        if info[i1].1.contains(" MATCH") {
                            let mut ext = true;
                            for i2 in 0..info.len() {
                                let full2 = &info[i2].4;
                                if !full1.ends_with(full2) && !full2.ends_with(full1) {
                                    ext = false;
                                }
                            }
                            if ext {
                                for i2 in 0..info.len() {
                                    if i2 != i1 {
                                        to_delete[i2] = true;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    erase_if(&mut info, &to_delete);
                }

                // Save.

                if !info.is_empty() {
                    // Record true positives.

                    let mut refs = Vec::<String>::new();
                    for m in 0..info.len() {
                        let full = &info[m].4;
                        if to_ref.contains_key(full) {
                            refs.push(to_ref[full].clone());
                        }
                    }
                    unique_sort(&mut refs);
                    res.5.push(format!("{}", refs.iter().format(",")));

                    // Record upstream and bases and chain types.

                    res.4.push(info[0].9.clone());
                    let mut bases = Vec::<Vec<u8>>::new();
                    let mut utr_tag = Vec::<Vec<u8>>::new();
                    for m in 0..info.len() {
                        bases.push(info[m].10.clone());
                        utr_tag.push(info[m].12.clone());
                    }
                    res.6.push(bases);
                    res.8.push(utr_tag);
                    let mut cts = Vec::<String>::new();
                    for m in 0..info.len() {
                        cts.push(info[m].11.clone());
                    }
                    res.7.push(cts);

                    // Save and advance.

                    for m in 0..info.len() {
                        log += &info[m].1;
                    }
                    res.1.push(log);
                }

                // Advance.

                j = k;
            }
        }
    });

    // Merge results.

    let mut all = Vec::new();
    let mut upstream = Vec::<Vec<u8>>::new();
    let mut bases = Vec::<Vec<Vec<u8>>>::new();
    let mut utr_tag = Vec::<Vec<Vec<u8>>>::new();
    let mut cts = Vec::<Vec<String>>::new();
    let mut refs = Vec::new();
    let mut start_motifs = Vec::<(String, Vec<u8>)>::new();
    for i in 0..results.len() {
        all.append(&mut results[i].1.clone());
        upstream.append(&mut results[i].4.clone());
        bases.append(&mut results[i].6.clone());
        utr_tag.append(&mut results[i].8.clone());
        cts.append(&mut results[i].7.clone());
        start_motifs.append(&mut results[i].3.clone());
        refs.append(&mut results[i].5.clone());
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // ANALYZE RECOMBINATION SITE
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Analyze upstream.  Here we are looking for recombination signal sequences, see
    // https://en.wikipedia.org/wiki/Recombination_signal_sequences.  However, these sequences
    // are not as conserved as the literature might suggest.  To locate them, we score using a
    // combination of proximity to the known sequences, and proximity to sequences similar to
    // a given one.

    let show_tags = false;
    if show_tags {
        println!("\nbegin upstream analysis, phase 1");
    }
    let mut tags = Vec::<Vec<u8>>::new();
    let mut tag_origin = Vec::<(usize, usize)>::new();
    for i in 0..upstream.len() {
        for j in 0..=upstream[i].len() - TAG_LEN {
            if upstream[i][j..].starts_with(&MOTIF1.to_vec()) {
                tag_origin.push((i, j));
                tags.push(upstream[i][j..j + TAG_LEN].to_vec());
            }
        }
    }
    if show_tags {
        println!("phase 1, found {} upstream tags", tags.len());
        println!("\nbegin upstream analysis, phase 2");
    }
    let mut tags2 = Vec::<Vec<u8>>::new();
    let mut extra = Vec::<Vec<u8>>::new();
    let mut tag_origin2 = Vec::<(usize, usize)>::new();

    fn cmp_tags(t1: &[u8], t2: &[u8]) -> usize {
        let mut x = 0;
        for k in 0..t1.len() - 1 {
            if t1[k] == t2[k] && t1[k + 1] == t2[k + 1] {
                x += 1;
            }
        }
        x
    }
    for i in 0..upstream.len() {
        let mut best = 0;
        let mut best_score = 0;
        for j in 0..=upstream[i].len() - TAG_LEN {
            let mut score = 0;
            for k in 0..tags.len() {
                score += cmp_tags(&upstream[i][j..j + TAG_LEN], &tags[k]);
            }
            if score > best_score {
                best = j;
                best_score = score;
            }
        }
        tags2.push(upstream[i][best..best + TAG_LEN].to_vec());
        extra.push(upstream[i][0..best].to_vec());
        tag_origin2.push((i, best));
    }
    if show_tags {
        println!("phase 2, found {} upstream tags\n", tags2.len());
    }
    let mut friends = vec![Vec::<usize>::new(); tags2.len()];
    for i in 0..tags2.len() {
        for j in 0..tags2.len() {
            if j == i {
                continue;
            }
            let mut x = cmp_tags(&tags2[i], &tags2[j]);
            let t1 = &tags2[i];
            let mut s = 0;
            for u in 0..MOTIF1.len() {
                if t1[u] == MOTIF1[u] {
                    s += 1;
                }
            }
            let mut s1 = s;
            for u in 0..MOTIF2.len() {
                if t1[MOTIF1.len() + 12 + u] == MOTIF2[u] {
                    s1 += 1;
                }
            }
            let mut s2 = s;
            for u in 0..MOTIF2.len() {
                if t1[MOTIF1.len() + 23 + u] == MOTIF2[u] {
                    s2 += 1;
                }
            }
            x += max(s1, s2);
            friends[i].push(x);
        }
        reverse_sort(&mut friends[i]);
    }
    let mut to_delete = vec![false; all.len()];
    for i in 0..tags2.len() {
        // known bads go up to 23, known goods go down to 27
        let good = friends[i][3] >= 25;
        if !good {
            to_delete[i] = true;
        }
        if show_tags {
            print!(
                "{}  {}  {}",
                strme(&tags2[i]),
                friends[i][0..10].iter().format(", "),
                refs[i]
            );
            if !good {
                print!(" JUNK");
            }
            println!();
        }
    }
    erase_if(&mut all, &to_delete);
    erase_if(&mut upstream, &to_delete);
    erase_if(&mut bases, &to_delete);
    erase_if(&mut utr_tag, &to_delete);
    erase_if(&mut extra, &to_delete);
    erase_if(&mut cts, &to_delete);
    erase_if(&mut refs, &to_delete);

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    //
    // ASSESS V GENE RESULTS
    //
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Gather reference hits.

    let mut finds = Vec::<String>::new();
    for i in 0..refs.len() {
        let r = refs[i].split(',').collect::<Vec<&str>>();
        for j in 0..r.len() {
            finds.push(r[j].to_string());
        }
    }

    // Collate.  Note that with the dna option on, the unique_sort option does something slightly
    // different, because it sees the DNA sequence, and hence results will be perturbed.

    unique_sort(&mut all);

    // Print if "fasta" argument specified.

    let (mut hcount, mut kcount, mut lcount) = (0, 0, 0);
    if print_fasta || store_fasta || print_aa {
        for i in 0..bases.len() {
            for j in 0..bases[i].len() {
                for pass in 0..2 {
                    let ct = &cts[i][j];
                    let count;
                    if ct == "IGH" {
                        if pass == 0 {
                            hcount += 1;
                        }
                        count = hcount;
                    } else if ct == "IGK" {
                        if pass == 0 {
                            kcount += 1;
                        }
                        count = kcount;
                    } else {
                        if pass == 0 {
                            lcount += 1;
                        }
                        count = lcount;
                    }
                    if !tenx {
                        fwrite!(fasta_log, ">{}V{}", ct, count);
                        if pass == 0 {
                            fwrite!(fasta_log, "-5'UTR");
                        }
                        fwriteln!(fasta_log, "");
                    } else {
                        fasta_count += 1;
                        let gene = format!("{}V{}", ct, count);
                        if pass == 1 {
                            fwriteln!(
                                fasta_log,
                                ">{}|{} enclone|{}|L-REGION+V-REGION|IG|{}|None|00",
                                fasta_count,
                                gene,
                                gene,
                                strme(&gene.as_bytes()[0..3]),
                            );
                        } else {
                            fwriteln!(
                                fasta_log,
                                ">{}|{} enclone|{}|5'UTR|IG|{}|None|00",
                                fasta_count,
                                gene,
                                gene,
                                strme(&gene.as_bytes()[0..3]),
                            );
                        }
                    }
                    if pass == 0 {
                        let b = utr_tag[i][j].clone();
                        fwriteln!(fasta_log, "{}", strme(&b));
                    } else {
                        let mut b = bases[i][j].clone();
                        b.append(&mut extra[i].clone());
                        fwriteln!(fasta_log, "{}", strme(&b));
                    }
                }
            }
        }
    } else {
        // Print in other cases.

        let mut mcount = 0;
        let mut nonsimples = 0;
        let mut unannotated = 0;
        let mut wrongs = 0;
        let mut ambs = 0;
        for k in 0..all.len() {
            let mut m = all[k].clone();
            let n = m.matches("==").count();
            let annotated = m.contains('🔴');
            let wrong = m.contains("WRONG CHAIN");
            if wrong {
                wrongs += 1;
            }
            let amb = m.matches('🌸').count() > n;
            if amb {
                ambs += 1;
            }
            let goods = m.matches('🔴').count()
                + m.matches(" MATCH").count()
                + m.matches(" PSEUDO").count();
            if goods == n && !wrong && (!annotated || !print_all) && !amb {
                continue;
            }
            if wrong || amb {
                mcount += 1;
                m = m.replace("[", &format!("[{}.", mcount));
                print!("\n{}", m);
                // print!(" upstream={}", strme(&upstream[k]));
                println!();
            } else if n == 1 && annotated && !print_all {
            } else if n > 0 {
                nonsimples += n;
                if annotated {
                    nonsimples -= 1;
                }
                mcount += 1;
                m = m.replace("[", &format!("[{}.", mcount));
                print!("\n{}", m);
                // print!(" upstream={}", strme(&upstream[k]));
                println!();
                if !annotated {
                    unannotated += 1;
                }
            }
        }

        // Find missing genes.

        unique_sort(&mut finds);
        println!();
        let mut missing = 0;
        let mut found = 0;
        for x in to_ref.iter() {
            if !bin_member(&finds, x.1) {
                // Exceptions because not seen in BI=6-12.

                if species == "human"
                    && (*x.1 == "251|IGKV1D-37".to_string()
                        || *x.1 == "285|IGKV3-7".to_string()
                        || *x.1 == "367|IGLV3-32".to_string()
                        || *x.1 == "382|IGLV5-48".to_string()
                        || x.1.ends_with("IGLV11-55"))
                {
                    continue;
                }

                // Mouse exceptions because not seen in BCR="1022446-1022449,1022518-1022521;77990;70838;1022418,1022419-1022421,1022434-1022437,1022506-1022509,1022490-1022493;1022426-1022428,1022442-1022445,1022498-1022501,1022514-1022516;1022414-1022417,1022430-1022433,1022502-1022505,1022486-1022489;1022422-1022425,1022438-1022441,1022494-1022497,1022510,1022512,1022513;1032106,1032107,1032110,1032111,1032098,1032099,1032102,1032103,1033023,1033024,1033027,1033028;1023684-1023691,1022722-1022729".

                if species == "mouse" {
                    if (*x.1).ends_with("IGHV1-62-1") {
                        continue;
                    }
                    if (*x.1).ends_with("IGHV7-2") {
                        continue;
                    }
                    if (*x.1).ends_with("IGKV4-60") {
                        continue;
                    }
                }

                // Exceptions because not in the reference.

                if species == "human" && *x.1 == "736|IGHV1-8".to_string() {
                    continue;
                }
                if species == "mouse" {
                    if (*x.1).ends_with("IGHV1-unknown1") {
                        continue;
                    }
                    if (*x.1).ends_with("IGHV12-1") {
                        continue;
                    }
                    if (*x.1).ends_with("IGHV8-9") {
                        continue;
                    }
                }

                // Exception.  The FWR3 ends with G instead of C, but data show C.

                if species == "mouse" && (*x.1).ends_with("IGHV8-2") {
                    continue;
                }

                // Declare missing.

                println!("failed to find {}", x.1);
                missing += 1;
            } else {
                found += 1;
            }
        }

        // Tally start motif info.

        if pwm {
            println!("let pwm = [");
            unique_sort(&mut start_motifs);
            for i in 0..start_motifs[0].1.len() {
                let mut a = 0;
                let mut c = 0;
                let mut g = 0;
                let mut t = 0;
                for j in 0..start_motifs.len() {
                    let b = start_motifs[j].1[i];
                    if b == b'A' {
                        a += 1;
                    } else if b == b'C' {
                        c += 1;
                    } else if b == b'G' {
                        g += 1;
                    } else if b == b'T' {
                        t += 1;
                    }
                }
                println!("    [{}, {}, {}, {}],", a, c, g, t);
            }
            println!("];\n");
        }

        // Print stats.

        println!("total missing genes = {}", missing);
        println!("total found genes = {}", found);
        println!("total unannotated = {}", unannotated);
        println!("total nonsimples = {}", nonsimples);
        println!("total wrongs = {}", wrongs);
        println!("total ambs = {}\n", ambs);

        // Test for regression.

        if species == "human" && (missing > 0 || nonsimples > 0 || wrongs > 0 || ambs > 0) {
            println!("REGRESSED!\n");
            std::process::exit(1);
        }
        if species == "mouse" && (missing > 0 || nonsimples > 9 || wrongs > 0 || ambs > 0) {
            println!("REGRESSED!\n");
            std::process::exit(1);
        }
    }

    // Output fasta for all genes.

    if print_fasta {
        print!("{}", strme(&fasta_log));
    }
    if store_fasta {
        let outname = format!("{}.fasta", id_name);
        let mut f = open_for_write_new![&format!("{}/{}", fasta_out_dir, outname)];
        fwrite!(f, "{}", strme(&fasta_log));
    }
}
