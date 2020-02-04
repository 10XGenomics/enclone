// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::*;
use h5::Dataset;
use mirror_sparse_matrix::*;
use regex::Regex;
use std::cmp::max;
use std::collections::HashMap;
use vector_utils::*;

// Clonotyping algorithm heuristics.

#[derive(Default)]
pub struct ClonotypeHeuristics {
    pub max_diffs: usize,
    pub ref_v_trim: usize,
    pub ref_j_trim: usize,
}

// Sample info data structure.

#[derive(Default)]
pub struct SampleInfo {
    // parallel vectors
    pub descrips: Vec<String>,     // map dataset index to dataset long name
    pub dataset_path: Vec<String>, // map dataset index to vdj path
    pub gex_path: Vec<String>,     // map dataset index to gex path
    pub dataset_id: Vec<String>,   // map dataset index to dataset short name
    pub donor_index: Vec<usize>,   // map dataset index to donor index
    pub donor_id: Vec<String>,     // map dataset index to donor short name
    pub sample_id: Vec<String>,    // map dataset id to sample short name
    // other
    pub dataset_list: Vec<Vec<usize>>, // map donor index to list of dataset indices
    pub donors: usize,                 // number of donors
    pub name_list: HashMap<String, Vec<usize>>, // map short name to list of dataset indices
}

impl SampleInfo {
    // number of datasets
    pub fn n(&self) -> usize {
        self.dataset_path.len()
    }
}

// Miscellaneous general options.

#[derive(Default)]
pub struct GeneralOpt {
    pub pre: String,
    pub insertions: bool,
    pub indels: bool,
    pub reannotate: bool,
    pub heavy_chain_reuse: bool,
    pub ngraph_filter: bool,
    pub graph: bool,
    pub utr_con: bool,
    pub con_con: bool,
    pub nwhitef: bool,
    pub weak: bool,
    pub tcr: bool,
    pub bcr: bool,
    pub exp: bool,
    pub reuse: bool,
    pub fasta: String,
    pub fasta_filename: String,
    pub min_cells_exact: usize,
    pub min_chains_exact: usize,
    pub exact: Option<usize>,
    pub binary: String,
    pub proto: String,
    pub h5: bool,
    pub no_reuse: bool,
    pub descrip: bool,
    pub ext: String,
    pub extc: HashMap<(String, String), String>,
    pub extn: HashMap<String, usize>,
    pub dref_file: String,
    pub mouse: bool,
    pub refname: String,
    pub noprint: bool,
    pub required_fps: Option<usize>,
    pub cellranger: bool,
    pub summary: bool,
}

// Allele finding algorithmic options.

#[derive(Default)]
pub struct AlleleAlgOpt {
    pub min_mult: usize,
    pub min_alt: usize,
}

// Allele finding print options.

#[derive(Default)]
pub struct AllelePrintOpt {
    pub con: bool,       // print alternate consensus sequences
    pub con_trace: bool, // tracing for con
}

// Join printing options.

#[derive(Default)]
pub struct JoinPrintOpt {
    pub seq: bool,     // print sequences of contigs, before truncation to V..J
    pub ann: bool,     // print annotations of contigs
    pub ann0: bool,    // print annotations of contigs, after truncation to V..J
    pub show_bc: bool, // show barcodes
    pub quiet: bool,   // don't print join events
    pub pfreq: usize,  // show data for 1/n joins even if correct
}

// Join algorithmic options.

#[derive(Default)]
pub struct JoinAlgOpt {
    pub max_score: f64,      // max score for join
    pub easy: bool,          // make joins even if core condition violated
    pub merge_onesies: bool, // create and merge onesies where completely unambiguous
    pub bcjoin: bool,        // join only by barcode identity
}

// Clonotype filtering options.

#[derive(Default)]
pub struct ClonoFiltOpt {
    pub ncells_low: usize,   // only show clonotypes with at least this many cells
    pub ncells_high: usize,  // only show clonotypes with at most this many cells
    pub min_umi: usize,      // only show clonotypes with at least this many UMIs in some contig
    pub min_datasets: usize, // only show clonotypes involving at least this many datasets
    pub min_chains: usize,   // only show clonotypes with at least this many chains
    pub max_chains: usize,   // only show clonotypes with at most this many chains
    pub ncross: bool,        // turn off cross filtering,
    pub cdr3: Option<Regex>, // only show clonotypes having one of these CDR3_AA sequences
    pub whitef: bool,        // only show clonotypes exhibiting whitelist contamination
    pub protect_bads: bool,  // protect bads from deletion
    pub fail_only: bool,     // only print fails
    pub seg: Vec<String>,    // only show clonotypes using one of these VDJ segment names
    pub segn: Vec<String>,   // only show clonotypes using one of these VDJ segment numbers
    pub min_exacts: usize,   // only show clonotypes having at least this many exact subclonotypes
    pub vj: Vec<u8>,         // only show clonotypes having exactly this full length V..J sequence
    pub vdup: bool,          // only show clonotypes having a same V segment in two chains
    pub have_onesie: bool,   // only show clonotypes including a onesie exact subclonotype
    pub cdiff: bool,         // only show clonotypes having a constant region difference
    pub del: bool,           // only show clonotypes exhibiting a deletion
    pub qual_filter: bool,   // filter out exact subclonotypes having a weak base
    pub weak_chains: bool,   // filter weak chains from clonotypes
    pub weak_onesies: bool,  // filter weak onesies
    pub weak_foursies: bool, // filter weak foursies
    pub bc_dup: bool,        // filter duplicated barcodes within an exact subclonotype
}

// Clonotype printing options.

#[derive(Default)]
pub struct ClonoPrintOpt {
    pub bu: bool,           // print barcodes and UMI counts
    pub seqc: bool,         // print V..J sequence for each chain if constant across clonotype
    pub full_seqc: bool,    // print contig sequence for each chain if constant across clonotype
    pub barcodes: bool,     // print the list of barcodes
    pub note_simple: bool,  // note if V..J is simple
    pub amino: Vec<String>, // categories for amino acid columns (per-chain per-exact subclonotype)
    pub cvars: Vec<String>, // per-chain per-exact-clonotype columns
    pub lvars: Vec<String>, // per-exact-clonotype ('lead') columns
    pub chain_brief: bool,  // show abbreviated chain headers
}

// Clonotype grouping options.

#[derive(Default)]
pub struct ClonoGroupOpt {
    pub heavy_cdr3_aa: bool, // group by perfect identity of cdr3_aa IGH or TRB
    pub min_group: usize,    // minimum number of clonotypes in group to print
}

// Parseable output options.

#[derive(Default)]
pub struct ParseableOpt {
    pub pout: String,            // name of parseable output file
    pub pchains: usize,          // number of chains to show in parseable output
    pub pcols: Vec<String>,      // column names to show in parseable output
    pub pcols_sort: Vec<String>, // sorted column names to show in parseable output
}

// Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
// run of enclone.

#[derive(Default)]
pub struct EncloneControl {
    pub gen_opt: GeneralOpt,              // miscellaneous general options
    pub pretty: bool,                     // use escape characters to enhance view
    pub silent: bool,                     // turn off extra logging
    pub force: bool,                      // make joins even if redundant
    pub comp: bool,                       // print computational performance stats
    pub debug_table_printing: bool,       // turn on debugging for table printing
    pub onesie_mult: usize,               // see main.rs
    pub merge_all_impropers: bool,        // merge all improper exact subclonotypes
    pub heur: ClonotypeHeuristics,        // algorithmic heuristics
    pub sample_info: SampleInfo,          // sample info
    pub allele_alg_opt: AlleleAlgOpt,     // algorithmic options for allele finding
    pub allele_print_opt: AllelePrintOpt, // print options for allele finding
    pub join_alg_opt: JoinAlgOpt,         // algorithmic options for join
    pub join_print_opt: JoinPrintOpt,     // printing options for join operations
    pub clono_filt_opt: ClonoFiltOpt,     // filtering options for clonotypes
    pub clono_print_opt: ClonoPrintOpt,   // printing options for clonotypes
    pub clono_group_opt: ClonoGroupOpt,   // grouping options for clonotypes
    pub parseable_opt: ParseableOpt,      // parseable output options
    pub toy: bool,                        // toy with phylogeny
}

// Set up data structure to track clonotype data.  A TigData is for one contig;
// a Vec<TigData> is for one barcode, and an ExactClonotype is for an exact subclonotype.

#[derive(Eq, Ord, PartialEq, PartialOrd, Default, Clone)] // not sure these are all needed
pub struct TigData {
    pub cdr3_dna: String,                     // CDR3 DNA sequence
    pub len: usize,                           // length of V..J sequence
    pub seq: Vec<u8>,                         // V..J contig subsequence
    pub v_start: usize,                       // start of V on full contig sequence
    pub v_stop: usize,                        // stop of aligned V on full contig sequence
    pub v_stop_ref: usize,                    // stop of aligned V on reference V
    pub j_start: usize,                       // start of aligned J on full contig sequence
    pub j_start_ref: usize,                   // start of aligned J on reference J
    pub j_stop: usize,                        // stop of J on full contig sequence
    pub c_start: Option<usize>,               // start of C on full contig sequence
    pub full_seq: Vec<u8>,                    // full contig sequence
    pub u_ref_id: Option<usize>,              // index of 5'-UTR in ref file if found
    pub v_ref_id: usize,                      // index of V segment reference sequence in ref file
    pub d_ref_id: Option<usize>,              // index of D segment reference sequence in ref file
    pub j_ref_id: usize,                      // index of J segment reference sequence in ref file
    pub c_ref_id: Option<usize>,              // index of C segment reference sequence in ref file
    pub cdr3_aa: String,                      // CDR3 amino acid sequence
    pub cdr3_start: usize,                    // start position in bases of CDR3 on V..J
    pub quals: Vec<u8>,                       // quality scores, truncated to V..J
    pub full_quals: Vec<u8>,                  // quality scores
    pub barcode: String,                      // barcode
    pub tigname: String,                      // name of contig
    pub left: bool,                           // true if this is IGH or TRA
    pub lena_index: usize,                    // index of lena id
    pub umi_count: usize,                     // number of UMIs supporting contig
    pub read_count: usize,                    // number of reads supporting contig
    pub chain_type: String,                   // e.g. IGH
    pub annv: Vec<(i32, i32, i32, i32, i32)>, // V annotation (one or two entries), for V..J
}

// The ExactClonotype data structure stores information that could be exhibited as a
// Vec<Vec<Vec<TigData>>>, but it avoids repetition of identical data.
//
// TigData0: data for each cell
// TigData1: shared data

pub struct TigData0 {
    pub quals: Vec<u8>,                       // quality scores, truncated to V..J
    pub v_start: usize,                       // start of V on full contig sequence
    pub j_stop: usize,                        // stop of J on full contig sequence
    pub c_start: Option<usize>,               // start of C on full contig sequence
    pub full_seq: Vec<u8>,                    // full contig sequence
    pub barcode: String,                      // barcode
    pub tigname: String,                      // name of contig
    pub lena_index: usize,                    // index of lena id
    pub umi_count: usize,                     // number of UMIs supporting contig
    pub read_count: usize,                    // number of reads supporting contig
}

pub struct TigData1 {
    pub cdr3_dna: String,                     // CDR3 DNA sequence
    pub seq: Vec<u8>,                         // V..J contig subsequence
    pub seq_del: Vec<u8>,                     // V..J, possibly with mod 3 del
    pub seq_del_amino: Vec<u8>,               // V..J, possibly with mod 3 del at mod 3 start
    pub full_seq: Vec<u8>,                    // full contig sequence (consensus)
    pub v_start: usize,                       // start of V on full contig sequence
    pub v_stop: usize,                        // stop of aligned V on full contig sequence
    pub v_stop_ref: usize,                    // stop of aligned V on reference V
    pub j_start: usize,                       // start of aligned J on full contig sequence
    pub j_start_ref: usize,                   // start of aligned J on reference J
    pub j_stop: usize,                        // stop of J on full contig sequence
    pub u_ref_id: Option<usize>,              // index of 5'-UTR in ref file if found
    pub v_ref_id: usize,                      // index of V segment reference sequence in ref file
    pub v_ref_id_donor: Option<usize>,        // optional index into alt_refs
    pub v_ref_id_donor_donor: Option<usize>,  // donor id for v_ref_id_donor
    pub v_ref_id_donor_alt_id: Option<usize>, // alt ref id for donor id for v_ref_id_donor
    pub d_ref_id: Option<usize>,              // index of D segment reference sequence in ref file
    pub j_ref_id: usize,                      // index of J segment reference sequence in ref file
    pub c_ref_id: Option<usize>,              // index of C segment reference sequence in ref file
    pub cdr3_aa: String,                      // CDR3 amino acid sequence
    pub cdr3_start: usize,                    // start position in bases of CDR3 on V..J
    pub left: bool,                           // true if this is IGH or TRA
    pub chain_type: String,                   // e.g. IGH
    pub annv: Vec<(i32, i32, i32, i32, i32)>, // V annotation (one or two entries), for V..J
    pub vs: DnaString,                        // reference V segment (possibly donor allele)
    pub vs_notesx: String, // notes on reference V segment (probably to be replaced)
    pub js: DnaString,     // reference J segment
}

pub struct ExactClonotype {
    pub share: Vec<TigData1>,       // clone info that is shared
    pub clones: Vec<Vec<TigData0>>, // clone info, excluding shared stuff
}

impl ExactClonotype {
    pub fn ncells(&self) -> usize {
        self.clones.len()
    }
    pub fn max_umi_count(&self) -> usize {
        let mut m = 0;
        for i in 0..self.clones.len() {
            for j in 0..self.clones[i].len() {
                m = max(m, self.clones[i][j].umi_count);
            }
        }
        m
    }
    pub fn lena_indices(&self) -> Vec<usize> {
        let mut x = Vec::<usize>::new();
        for i in 0..self.clones.len() {
            x.push(self.clones[i][0].lena_index);
        }
        unique_sort(&mut x);
        x
    }
}

// Define clonotype info data structure.  The fact that we have multiple data structures
// encapsulating the same info is legacy and to be cleaned up.
//
// The vectors in a CloneInfo object mostly have length two.  The exceptions are in
// improper clones (not having chains of both types).

#[derive(Eq, Ord, PartialEq, PartialOrd, Default, Clone)] // not sure we need all these
pub struct CloneInfo {
    pub lens: Vec<usize>,   // V..J contig lengths (will sort by this)
    pub tigs: Vec<Vec<u8>>, // contigs, truncated to V..J (with possible - chars inserted)
    // note only need tigs in has_del case, could improve
    pub tigs_amino: Vec<Vec<u8>>, // same as tigs, but deletion shifted to mod 3 position
    // if there is one (rare, so wasteful, should be Option)
    pub tigsp: Vec<DnaString>, // contigs, truncated to V..J, packed (doesn't show - chars)
    pub has_del: Vec<bool>,    // if - chars inserted to represent deletion
    pub orig_tigs: Vec<DnaString>, // untruncated contigs
    pub clonotype_id: usize,   // index into exact_clonotypes
    pub exact_cols: Vec<usize>, // the columns of the exact_clonotype that were extracted (used?)
    pub clonotype_index: usize, // index into vector of all clonotypes (across samples)
    pub origin: Vec<usize>,    // sample indices
    pub vs: Vec<DnaString>,    // reference V segments (possibly donor allele)
    pub dref: Vec<Option<usize>>, // indices into alt_refs
    pub vs_notesx: Vec<String>, // notes on reference V segments (probably to be replaced)
    pub js: Vec<DnaString>,    // reference J segments
    pub vsids: Vec<usize>,     // ids of V segments
    pub jsids: Vec<usize>,     // ids of J segments
    pub cdr3s: Vec<String>,    // cdr3 nucleotide seqs
    pub cdr3_aa: Vec<String>,  // cdr3 amino acid seqs
    pub chain_types: Vec<String>, // chain types
}

// Gene expression stuff.

#[derive(Default)]
pub struct GexInfo {
    pub gex_features: Vec<Vec<String>>,
    pub gex_barcodes: Vec<Vec<String>>,
    pub gex_matrices: Vec<MirrorSparseMatrix>,
    pub gex_mults: Vec<f64>,
    pub fb_mults: Vec<f64>,
    pub h5_data: Vec<Option<Dataset>>,
    pub h5_indices: Vec<Option<Dataset>>,
    pub h5_indptr: Vec<Vec<u32>>,
    pub is_gex: Vec<Vec<bool>>,
    pub feature_id: Vec<HashMap<String, usize>>,
}

// Every entry in a ColInfo is a vector whose number of entries is the number of chains
// in a clonotype.

#[derive(Default)]
pub struct ColInfo {
    pub uids: Vec<Option<usize>>,
    pub vids: Vec<usize>,
    pub vpids: Vec<Option<usize>>,
    pub dids: Vec<Option<usize>>,
    pub jids: Vec<usize>,
    pub cids: Vec<Option<usize>>,
    pub cdr3_starts: Vec<usize>,
    pub cdr3_lens: Vec<usize>,
    pub seq_lens: Vec<usize>, // not sure we should be computing or using this
    pub seq_del_lens: Vec<usize>,
    pub seqss: Vec<Vec<Vec<u8>>>,
    pub seqss_amino: Vec<Vec<Vec<u8>>>,
    pub chain_descrip: Vec<String>,
    pub mat: Vec<Vec<Option<usize>>>,
    pub cvars: Vec<Vec<String>>,
}
