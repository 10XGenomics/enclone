// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::*;
use hdf5::Dataset;
use mirror_sparse_matrix::*;
use regex::Regex;
use std::cmp::max;
use std::collections::HashMap;
use string_utils::*;
use vector_utils::*;

// Field (variable) names.

pub const LVARS_ALLOWED: [&str; 18] = [
    "datasets",
    "samples",
    "donors",
    "n",
    "gex",
    "gex_min",
    "gex_max",
    "gex_μ",
    "gex_Σ",
    "gex_cell",
    "n_gex_cell",
    "n_gex",
    "clust",
    "type",
    "entropy",
    "near",
    "far",
    "ext",
];

pub const CVARS_ALLOWED: [&str; 24] = [
    "var", "u", "u_min", "u_max", "u_Σ", "u_μ", "comp", "edit", "r", "r_min", "r_max", "r_Σ",
    "r_μ", "const", "white", "cdr3_dna", "ulen", "vjlen", "clen", "cdiff", "udiff", "notes",
    "d_univ", "d_donor",
];

pub const CVARS_ALLOWED_PCELL: [&str; 2] = ["u_cell", "r_cell"];

pub const PLVARS_ALLOWED: [&str; 7] = [
    "group_id",
    "group_ncells",
    "clonotype_id",
    "clonotype_ncells",
    "nchains",
    "exact_subclonotype_id",
    "barcodes",
];

pub const PCVARS_ALLOWED: [&str; 19] = [
    "v_name",
    "d_name",
    "j_name",
    "v_id",
    "d_id",
    "j_id",
    "var_indices_dna",
    "var_indices_aa",
    "share_indices_dna",
    "share_indices_aa",
    "v_start",
    "const_id",
    "utr_id",
    "utr_name",
    "cdr3_start",
    "cdr3_aa",
    "seq",
    "vj_seq",
    "var_aa",
];

// Clonotyping algorithm heuristics.

#[derive(Default)]
pub struct ClonotypeHeuristics {
    pub max_diffs: usize,
    pub ref_v_trim: usize,
    pub ref_j_trim: usize,
}

#[derive(Clone)]
pub struct LinearCondition {
    pub coeff: Vec<f64>,  // left hand side (lhs) coefficients
    pub var: Vec<String>, // left hand side variables (parallel to coefficients)
    pub rhs: f64,         // right hand side; sum of lhs must exceed rhs
    pub sense: String,    // le, ge, lt, gt
}

impl LinearCondition {
    pub fn n(&self) -> usize {
        self.coeff.len()
    }

    pub fn new(x: &str) -> LinearCondition {
        let y = x.replace(" ", "");
        let lhs: String;
        let mut rhs: String;
        let sense: String;
        if y.contains(">=") {
            lhs = y.before(">=").to_string();
            rhs = y.after(">=").to_string();
            sense = "ge".to_string();
        } else if y.contains('≥') {
            lhs = y.before("≥").to_string();
            rhs = y.after("≥").to_string();
            sense = "ge".to_string();
        } else if y.contains("<=") {
            lhs = y.before("<=").to_string();
            rhs = y.after("<=").to_string();
            sense = "le".to_string();
        } else if y.contains('≤') {
            lhs = y.before("≤").to_string();
            rhs = y.after("≤").to_string();
            sense = "le".to_string();
        } else if y.contains('<') {
            lhs = y.before("<").to_string();
            rhs = y.after("<").to_string();
            sense = "lt".to_string();
        } else if y.contains('>') {
            lhs = y.before(">").to_string();
            rhs = y.after(">").to_string();
            sense = "gt".to_string();
        } else {
            eprintln!(
                "\nImproperly formatted condition, no inequality symbol, \
                 please type \"enclone help display\": {}.\n",
                x
            );
            std::process::exit(1);
        }
        if !rhs.contains('.') {
            rhs += ".0";
        }
        if !rhs.parse::<f64>().is_ok() {
            eprintln!(
                "\nImproperly formatted condition, right-hand side invalid: {}.\n",
                x
            );
            std::process::exit(1);
        }
        let rhs = rhs.force_f64();
        let mut parts = Vec::<String>::new();
        let mut last = 0;
        let lhsx = lhs.as_bytes();
        let mut parens = 0 as isize;
        for i in 0..lhsx.len() {
            if i > 0 && parens == 0 && (lhsx[i] == b'+' || lhsx[i] == b'-') {
                if lhsx[last] != b'+' {
                    parts.push(stringme(&lhsx[last..i]));
                } else {
                    parts.push(stringme(&lhsx[last + 1..i]));
                }
                last = i;
            }
            if lhsx[i] == b'(' {
                parens += 1;
            } else if lhsx[i] == b')' {
                parens -= 1;
            }
        }
        let mut coeff = Vec::<f64>::new();
        let mut var = Vec::<String>::new();
        parts.push(lhs[last..].to_string());
        for i in 0..parts.len() {
            parts[i] = parts[i].replace("(", "");
            parts[i] = parts[i].replace(")", "");
            if parts[i].contains('*') {
                let mut coeffi = parts[i].before("*").to_string();
                let vari = parts[i].after("*");
                if !coeffi.contains('.') {
                    coeffi += ".0";
                }
                if !coeffi.parse::<f64>().is_ok() {
                    eprintln!(
                        "\nImproperly formatted condition, coefficient {} is invalid: {}.\n",
                        coeffi, x
                    );
                    std::process::exit(1);
                }
                coeff.push(coeffi.force_f64());
                var.push(vari.to_string());
            } else {
                let mut coeffi = 1.0;
                let mut start = 0;
                if parts[i].starts_with('-') {
                    coeffi = -1.0;
                    start = 1;
                }
                coeff.push(coeffi);
                var.push(parts[i][start..].to_string());
            }
        }
        LinearCondition {
            coeff: coeff,
            var: var,
            rhs: rhs,
            sense: sense,
        }
    }

    pub fn satisfied(&self, val: &Vec<f64>) -> bool {
        let mut lhs = 0.0;
        for i in 0..self.coeff.len() {
            lhs += self.coeff[i] * val[i];
        }
        if self.sense == "lt".to_string() {
            return lhs < self.rhs;
        } else if self.sense == "gt".to_string() {
            return lhs > self.rhs;
        } else if self.sense == "le".to_string() {
            return lhs <= self.rhs;
        } else {
            return lhs >= self.rhs;
        }
    }

    pub fn require_valid_variables(&self, ctl: &EncloneControl) {
        let lvars = &ctl.clono_print_opt.lvars;
        let mut lvars0 = Vec::<String>::new();
        let exclude = vec![
            "datasets",
            "donors",
            "near",
            "far",
            "n_gex_cell",
            "n_gex",
            "clust",
            "type",
            "gex",
            "gex_min",
            "gex_max",
            "gex_mean",
            "gex_sum",
            "entropy",
            "ext",
        ];
        for j in 0..lvars.len() {
            let mut ok = true;
            for m in 0..exclude.len() {
                if lvars[j] == exclude[m] {
                    ok = false;
                }
            }
            if lvars[j].starts_with("g") && lvars[j].after("g").parse::<usize>().is_ok() {
                ok = false;
            }
            if ok {
                lvars0.push(lvars[j].clone());
            }
        }
        unique_sort(&mut lvars0);
        for i in 0..self.var.len() {
            if !bin_member(&lvars0, &self.var[i]) {
                eprintln!(
                    "\nFound invalid variable {} in linear condition.\n",
                    self.var[i]
                );
                std::process::exit(1);
            }
        }
    }
}

// Sample info data structure.

#[derive(Default)]
pub struct SampleInfo {
    // parallel vectors
    pub descrips: Vec<String>,     // map dataset index to dataset long name
    pub dataset_path: Vec<String>, // map dataset index to vdj path
    pub gex_path: Vec<String>,     // map dataset index to gex path
    pub dataset_id: Vec<String>,   // map dataset index to dataset short name
    pub donor_id: Vec<String>,     // map dataset index to donor short name
    pub sample_id: Vec<String>,    // map dataset id to sample short name
    pub color: Vec<String>,        // map dataset to color
    // other
    pub dataset_list: Vec<String>, // unique-sorted list of dataset short names
    pub sample_list: Vec<String>,  // unique-sorted list of sample short names
    pub donor_list: Vec<String>,   // unique-sorted list of donor short names
    pub tag_list: Vec<String>,     // unique-sorted list of tag short names
    pub sample_donor_list: Vec<(usize, usize)>, // unique-sorted list of (sample, donor) indices
    pub donors: usize,             // number of donors
    // map dataset index to map of barcode to (sample,donor):
    pub sample_donor: Vec<HashMap<String, (String, String)>>,
    // map dataset index to map of barcode to tag:
    pub tag: Vec<HashMap<String, String>>,
    // map dataset index to map of barcode to color:
    pub barcode_color: Vec<HashMap<String, String>>,
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
    pub fasta_aa_filename: String,
    pub min_cells_exact: usize,
    pub min_chains_exact: usize,
    pub exact: Option<usize>,
    pub binary: String,
    pub proto: String,
    pub h5: bool,
    pub h5_pre: bool,
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
    pub summary_clean: bool,
    pub cr_version: String,
    pub nwarn: bool,
    pub gene_scan_test: Option<LinearCondition>,
    pub gene_scan_control: Option<LinearCondition>,
    pub gene_scan_threshold: Option<LinearCondition>,
    pub plot_file: String,
    pub sample_color_map: HashMap<String, String>,
    pub use_legend: bool,
    pub legend: Vec<(String, String)>,
    pub accept_inconsistent: bool, // TEMPORARY!
    pub current_ref: bool,         // TEMPORARY!
    pub internal_run: bool,
    pub force_h5: bool,
    pub full_counts: bool,
    pub html: bool,
    pub stable_doc: bool,
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
    pub max_datasets: usize, // only show clonotypes involving at most this many datasets
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
    pub donor: bool,         // allow cells from different donors to be placed in the same clonotype
    pub bounds: Vec<LinearCondition>, // bounds on certain variables
    pub barcode: Vec<String>, // requires one of these barcodes
}

// Clonotype printing options.

#[derive(Default)]
pub struct ClonoPrintOpt {
    pub bu: bool,                          // print barcodes and UMI counts
    pub seqc: bool, // print V..J sequence for each chain if constant across clonotype
    pub full_seqc: bool, // print contig sequence for each chain if constant across clonotype
    pub barcodes: bool, // print the list of barcodes
    pub note_simple: bool, // note if V..J is simple
    pub amino: Vec<String>, // categories for amino acid columns (per-chain per-exact subclonotype)
    pub cvars: Vec<String>, // per-chain per-exact-clonotype columns
    pub lvars: Vec<String>, // per-exact-clonotype ('lead') columns
    pub lvars_match: Vec<Vec<Vec<usize>>>, // matching features for <regular expression>_g etc.
    pub chain_brief: bool, // show abbreviated chain headers
    pub sum: bool,  // print sum row
    pub mean: bool, // print mean row
}

// Clonotype grouping options.

#[derive(Default)]
pub struct ClonoGroupOpt {
    pub heavy_cdr3_aa: bool, // group by perfect identity of cdr3_aa IGH or TRB
    pub vj_refname: bool,    // group by having the same VJ reference names
    pub min_group: usize,    // minimum number of clonotypes in group to print
}

// Parseable output options.

#[derive(Default)]
pub struct ParseableOpt {
    pub pout: String,            // name of parseable output file
    pub pchains: usize,          // number of chains to show in parseable output
    pub pcols: Vec<String>,      // column names to show in parseable output
    pub pcols_sort: Vec<String>, // sorted column names to show in parseable output
    pub pbarcode: bool,          // generate output per barcode rather than per exact subclonotype
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
    pub comp2: bool,                      // print more detailed computational performance stats
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
    pub dataset_index: usize,                 // index of dataset
    pub sample_index: Option<usize>,          // index of sample
    pub donor_index: Option<usize>,           // index of donor
    pub tag_index: Option<usize>,             // index of tag
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

#[derive(Clone)]
pub struct TigData0 {
    pub quals: Vec<u8>,              // quality scores, truncated to V..J
    pub v_start: usize,              // start of V on full contig sequence
    pub j_stop: usize,               // stop of J on full contig sequence
    pub c_start: Option<usize>,      // start of C on full contig sequence
    pub full_seq: Vec<u8>,           // full contig sequence
    pub barcode: String,             // barcode
    pub tigname: String,             // name of contig
    pub dataset_index: usize,        // index of dataset
    pub sample_index: Option<usize>, // index of sample
    pub donor_index: Option<usize>,  // index of donor
    pub tag_index: Option<usize>,    // index of tag
    pub umi_count: usize,            // number of UMIs supporting contig
    pub read_count: usize,           // number of reads supporting contig
}

#[derive(Clone)]
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

#[derive(Clone)]
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
    pub fn dataset_indices(&self) -> Vec<usize> {
        let mut x = Vec::<usize>::new();
        for i in 0..self.clones.len() {
            x.push(self.clones[i][0].dataset_index);
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
    pub clonotype_index: usize, // index into vector of all exact subclonotypes (across samples)
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
    pub gex_cell_barcodes: Vec<Vec<String>>,
    pub cluster: Vec<HashMap<String, usize>>,
    pub cell_type: Vec<HashMap<String, String>>,
    pub gex_mults: Vec<f64>,
    pub fb_mults: Vec<f64>,
    pub h5_data: Vec<Option<Dataset>>,
    pub h5_indices: Vec<Option<Dataset>>,
    pub h5_indptr: Vec<Vec<u32>>,
    pub is_gex: Vec<Vec<bool>>,
    pub feature_id: Vec<HashMap<String, usize>>,
    pub have_gex: bool,
    pub have_fb: bool,
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
