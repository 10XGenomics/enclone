// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::*;
use hdf5::Dataset;
use mirror_sparse_matrix::*;
use perf_stats::*;
use regex::Regex;
use std::cmp::max;
use std::collections::HashMap;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

// Field (variable) names.
// Lead variables for exact subclonotypes and cells.
pub const LVARS_ALLOWED: [&str; 29] = [
    "datasets",
    "origins",
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
    "n_b",
    "clust",
    "cred",
    "cred_cell",
    "type",
    "entropy",
    "entropy_cell",
    "near",
    "far",
    "dref",
    "dref_aa",
    "ext",
    "mark",
    "inkt",
    "mait",
    "sec",
    "mem",
];

// Chain variables that can be used for contigs and chains

pub const CVARS_ALLOWED: [&str; 46] = [
    "var", "u", "u_min", "u_max", "u_Σ", "u_μ", "comp", "edit", "r", "r_min", "r_max", "r_Σ",
    "r_μ", "const", "white", "cdr1_dna", "cdr2_dna", "cdr3_dna", "cdr1_len", "cdr2_len",
    "cdr3_len", "cdr1_aa", "cdr2_aa", "cdr3_aa", "fwr1_dna", "fwr2_dna", "fwr3_dna", "fwr4_dna",
    "fwr1_len", "fwr2_len", "fwr3_len", "fwr4_len", "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa",
    "ulen", "vjlen", "clen", "cdiff", "udiff", "notes", "d_univ", "d_donor", "aa%", "dna%",
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

pub const PCVARS_ALLOWED: [&str; 22] = [
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
    "vj_seq_nl",
    "vj_aa",
    "vj_aa_nl",
    "var_aa",
];

// Clonotyping algorithm heuristics.

#[derive(Default)]
pub struct ClonotypeHeuristics {
    pub max_diffs: usize,
    pub max_degradation: usize,
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
            "dref",
            "dref_aa",
            "n_gex_cell",
            "n_gex",
            "n_b",
            "clust",
            "cred",
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

// Origin info data structure.

#[derive(Default)]
pub struct OriginInfo {
    // parallel vectors
    pub descrips: Vec<String>,     // map dataset index to dataset long name
    pub dataset_path: Vec<String>, // map dataset index to vdj path
    pub gex_path: Vec<String>,     // map dataset index to gex path
    pub dataset_id: Vec<String>,   // map dataset index to dataset short name
    pub donor_id: Vec<String>,     // map dataset index to donor short name
    pub origin_id: Vec<String>,    // map dataset id to origin (sample) short name
    pub color: Vec<String>,        // map dataset to color
    // other
    pub dataset_list: Vec<String>, // unique-sorted list of dataset short names
    pub origin_list: Vec<String>,  // unique-sorted list of origin (sample) short names
    pub donor_list: Vec<String>,   // unique-sorted list of donor short names
    pub tag_list: Vec<String>,     // unique-sorted list of tag short names
    pub donors: usize,             // number of donors
    // map dataset index to map of barcode to origin:
    pub origin_for_bc: Vec<HashMap<String, String>>,
    // map dataset index to map of barcode to donor:
    pub donor_for_bc: Vec<HashMap<String, String>>,
    // map dataset index to map of barcode to tag:
    pub tag: Vec<HashMap<String, String>>,
    // map dataset index to map of barcode to color:
    pub barcode_color: Vec<HashMap<String, String>>,
    pub alt_bc_fields: Vec<Vec<(String, HashMap<String, String>)>>,
    pub cells_cellranger: Vec<Option<usize>>,
    pub mean_read_pairs_per_cell_cellranger: Vec<Option<usize>>,
    // map dataset index to a map of barcode to (secreted, membrane) UMI counts
    pub secmem: Vec<HashMap<String, (usize, usize)>>,
}

impl OriginInfo {
    // number of datasets
    pub fn n(&self) -> usize {
        self.dataset_path.len()
    }
}

// Miscellaneous general options.

#[derive(Default)]
pub struct GeneralOpt {
    pub pre: Vec<String>,
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
    pub clustal_aa: String,
    pub clustal_dna: String,
    pub phylip_aa: String,
    pub phylip_dna: String,
    pub min_cells_exact: usize,
    pub min_chains_exact: usize,
    pub chains_exact: usize,
    pub complete: bool,
    pub exact: Option<usize>,
    pub binary: String,
    pub proto: String,
    // Optional path to a json file containing metadata
    pub proto_metadata: Option<String>,
    pub h5: bool,
    pub h5_pre: bool,
    pub accept_reuse: bool,
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
    pub summary_csv: bool,
    pub cr_version: String,
    pub nwarn: bool,
    pub gene_scan_test: Option<LinearCondition>,
    pub gene_scan_control: Option<LinearCondition>,
    pub gene_scan_threshold: Option<LinearCondition>,
    pub plot_file: String,
    pub plot_by_isotype: bool,
    pub plot_by_mark: bool,
    pub origin_color_map: HashMap<String, String>,
    pub use_legend: bool,
    pub legend: Vec<(String, String)>,
    pub accept_inconsistent: bool, // TEMPORARY!
    pub current_ref: bool,         // TEMPORARY!
    pub internal_run: bool,
    pub force_h5: bool,
    pub full_counts: bool,
    pub html: bool,
    pub html_title: String,
    pub svg: bool,
    pub stable_doc: bool,
    pub imgt: bool,
    pub imgt_fix: bool,
    pub ngroup: bool,
    pub jc1: bool,
    pub trace_barcode: String,
    pub ncell: bool,
    pub baseline: bool,
    pub echo: bool,
    pub mark_stats: bool,
    pub mark_stats2: bool,
    pub print_cpu: bool,
    pub print_cpu_info: bool,
    pub newick: bool,
    pub tree: String,
    pub allow_inconsistent: bool,
    pub color: String,
    pub species: String, // human or mouse or unknown, determined from the reference sequence
    pub using_secmem: bool,
    pub const_igh: Option<Regex>,
    pub const_igkl: Option<Regex>,
    pub diff_style: String,
}

// Allele-finding algorithmic options.

#[derive(Default)]
pub struct AlleleAlgOpt {
    pub min_mult: usize,
    pub min_alt: usize,
}

// Allele-finding print options.

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
    pub max_score: f64,          // max score for join
    pub easy: bool,              // make joins even if core condition violated
    pub merge_onesies: bool,     // create and merge onesies where completely unambiguous
    pub merge_onesies_ctl: bool, // restriction on onesie merger
    pub bcjoin: bool,            // join only by barcode identity
    pub max_cdr3_diffs: usize,
}

// Clonotype filtering options.
// These fall into 2 categories: 1) on by default and 2) user-specified.

#[derive(Default)]
pub struct ClonoFiltOpt {
    pub ncells_low: usize,   // only show clonotypes with at least this many cells
    pub ncells_high: usize,  // only show clonotypes with at most this many cells
    pub min_umi: usize,      // only show clonotypes with at least this many UMIs in some contig
    pub min_datasets: usize, // only show clonotypes involving at least this many datasets
    pub max_datasets: usize, // only show clonotypes involving at most this many datasets
    pub min_dataset_ratio: usize, // see "enclone help filter"
    pub min_chains: usize,   // only show clonotypes with at least this many chains
    pub max_chains: usize,   // only show clonotypes with at most this many chains
    pub ngex: bool,          // turn off gex filtering,
    pub ncross: bool,        // turn off cross filtering,
    pub cdr3: Option<Regex>, // only show clonotypes whose CDR3_AA matches regular expression
    pub cdr3_lev: String,    // only show clonotypes whose CDR3_AA matches Levenshtein dist pattern
    pub whitef: bool,        // only show clonotypes exhibiting whitelist contamination
    pub protect_bads: bool,  // protect bads from deletion
    pub fail_only: bool,     // only print fails
    pub seg: Vec<Vec<String>>, // only show clonotypes using one of these VDJ segment names
    pub segn: Vec<Vec<String>>, // only show clonotypes using one of these VDJ segment numbers
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
    pub umi_filt: bool,      // umi count filter
    pub umi_filt_mark: bool, // umi count filter (but only mark)
    pub non_cell_mark: bool,
    pub marked: bool,                 // only print clonotypes having a mark
    pub marked_b: bool, // only print clonotypes having a mark and which are typed as B cells
    pub umi_ratio_filt: bool, // umi ratio filter
    pub umi_ratio_filt_mark: bool, // umi ratio filter (but only mark)
    pub fcell: Vec<(String, String)>, // constaints from FCELL
    pub inkt: bool,
    pub mait: bool,
}

// Clonotype printing options.

#[derive(Default)]
pub struct ClonoPrintOpt {
    pub bu: bool,                                      // print barcodes and UMI counts
    pub seqc: bool, // print V..J sequence for each chain if constant across clonotype
    pub full_seqc: bool, // print contig sequence for each chain if constant across clonotype
    pub barcodes: bool, // print the list of barcodes
    pub note_simple: bool, // note if V..J is simple
    pub amino: Vec<String>, // categories for amino acid columns (per-chain per-exact subclonotype)
    pub cvars: Vec<String>, // per-chain per-exact-clonotype columns
    pub lvars: Vec<String>, // per-exact-clonotype ('lead') columns
    pub regex_match: Vec<HashMap<String, Vec<usize>>>, // matching features for <regex>_g etc.
    pub chain_brief: bool, // show abbreviated chain headers
    pub sum: bool,  // print sum row
    pub mean: bool, // print mean row
}

// Clonotype grouping options.

#[derive(Default)]
pub struct ClonoGroupOpt {
    pub heavy_cdr3_aa: bool, // group by perfect identity of cdr3_aa IGH or TRB
    pub vj_refname: bool,    // group by having the same VJ reference names
    pub vj_refname_strong: bool, // group by having the same VJ reference names, but stronger
    pub min_group: usize,    // minimum number of clonotypes in group to print
}

// Parseable output options.

#[derive(Default)]
pub struct ParseableOpt {
    pub pout: String,             // name of parseable output file
    pub pchains: usize,           // number of chains to show in parseable output
    pub pcols: Vec<String>,       // column names to show in parseable output
    pub pcols_sort: Vec<String>,  // sorted column names to show in parseable output
    pub pcols_sortx: Vec<String>, // same but before colon if present
    pub pbarcode: bool,           // generate output per barcode rather than per exact subclonotype
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
    pub merge_all_impropers: bool,        // merge all improper exact subclonotypes
    pub heur: ClonotypeHeuristics,        // algorithmic heuristics
    pub origin_info: OriginInfo,          // origin (sample) info
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

pub static mut WALLCLOCK: f64 = 0.0;

impl EncloneControl {
    pub fn perf_stats(&self, t: &Instant, msg: &str) {
        let used = elapsed(&t);
        if self.comp {
            println!(
                "used {:.2} seconds {}, peak mem = {:.2} GB",
                used,
                msg,
                peak_mem_usage_gb()
            );
        }
        unsafe {
            WALLCLOCK += used;
        }
    }
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
    pub fr1_start: usize,                     // start position in bases of FWR1 on V..J
    pub cdr1_start: usize,                    // start position in bases of CDR1 on V..J
    pub fr2_start: usize,                     // start position in bases of FWR2 on V..J
    pub cdr2_start: usize,                    // start position in bases of CDR2 on V..J
    pub fr3_start: usize,                     // start position in bases of FWR3 on V..J
    pub cdr3_aa: String,                      // CDR3 amino acid sequence
    pub cdr3_start: usize,                    // start position in bases of CDR3 on V..J
    pub quals: Vec<u8>,                       // quality scores, truncated to V..J
    pub full_quals: Vec<u8>,                  // quality scores
    pub barcode: String,                      // barcode
    pub tigname: String,                      // name of contig
    pub left: bool,                           // true if this is IGH or TRB
    pub dataset_index: usize,                 // index of dataset
    pub origin_index: Option<usize>,          // index of origin (sample)
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
    pub origin_index: Option<usize>, // index of origin (sample)
    pub donor_index: Option<usize>,  // index of donor
    pub tag_index: Option<usize>,    // index of tag
    pub umi_count: usize,            // number of UMIs supporting contig
    pub read_count: usize,           // number of reads supporting contig
    pub marked: bool,                // if marked for possible deletion
}

#[derive(Clone)]
pub struct TigData1 {
    pub cdr3_dna: String,           // CDR3 DNA sequence
    pub seq: Vec<u8>,               // V..J contig subsequence
    pub seq_del: Vec<u8>,           // V..J, possibly with mod 3 del
    pub seq_del_amino: Vec<u8>,     // V..J, possibly with mod 3 del at mod 3 start
    pub aa_mod_indel: Vec<u8>,      // amino acid sequence, after removing indel if present
    pub ins: Vec<(usize, Vec<u8>)>, // insertions in V..J (currently at most one)
    // **before** the given position
    pub full_seq: Vec<u8>,                   // full contig sequence (consensus)
    pub v_start: usize,                      // start of V on full contig sequence
    pub v_stop: usize,                       // stop of aligned V on full contig sequence
    pub v_stop_ref: usize,                   // stop of aligned V on reference V
    pub j_start: usize,                      // start of aligned J on full contig sequence
    pub j_start_ref: usize,                  // start of aligned J on reference J
    pub j_stop: usize,                       // stop of J on full contig sequence
    pub u_ref_id: Option<usize>,             // index of 5'-UTR in ref file if found
    pub v_ref_id: usize,                     // index of V segment reference sequence in ref file
    pub v_ref_id_donor: Option<usize>,       // optional index into alt_refs
    pub v_ref_id_donor_donor: Option<usize>, // donor id for v_ref_id_donor
    pub v_ref_id_donor_alt_id: Option<usize>, // alt ref id for donor id for v_ref_id_donor
    pub d_ref_id: Option<usize>,             // index of D segment reference sequence in ref file
    pub j_ref_id: usize,                     // index of J segment reference sequence in ref file
    pub c_ref_id: Option<usize>,             // index of C segment reference sequence in ref file
    pub fr1_start: usize,                    // start position in bases of FWR1 on V..J
    pub cdr1_start: usize,                   // start position in bases of CDR1 on V..J
    pub fr2_start: usize,                    // start position in bases of FWR2 on V..J
    pub cdr2_start: usize,                   // start position in bases of CDR2 on V..J
    pub fr3_start: usize,                    // start position in bases of FWR3 on V..J
    pub cdr3_aa: String,                     // CDR3 amino acid sequence
    pub cdr3_start: usize,                   // start position in bases of CDR3 on V..J
    pub left: bool,                          // true if this is IGH or TRB
    pub chain_type: String,                  // e.g. IGH
    pub annv: Vec<(i32, i32, i32, i32, i32)>, // V annotation (one or two entries), for V..J
    pub vs: DnaString,                       // reference V segment (possibly donor allele)
    pub vs_notesx: String, // notes on reference V segment (probably to be replaced)
    pub js: DnaString,     // reference J segment
    pub inkt_alpha_chain_gene_match: bool,
    pub inkt_alpha_chain_junction_match: bool,
    pub inkt_beta_chain_gene_match: bool,
    pub inkt_beta_chain_junction_match: bool,
    pub mait_alpha_chain_gene_match: bool,
    pub mait_alpha_chain_junction_match: bool,
    pub mait_beta_chain_gene_match: bool,
    pub mait_beta_chain_junction_match: bool,
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
    pub fn nchains(&self) -> usize {
        self.share.len()
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
    pub clonotype_index: usize, // index into vector of all exact subclonotypes (across origins)
    pub origin: Vec<usize>,    // origin indices
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

// Gene expression and feature barcoding stuff.

#[derive(Default)]
pub struct GexInfo {
    pub gex_features: Vec<Vec<String>>,
    pub gex_barcodes: Vec<Vec<String>>,
    pub gex_matrices: Vec<MirrorSparseMatrix>,
    pub gex_cell_barcodes: Vec<Vec<String>>,
    pub cluster: Vec<HashMap<String, usize>>,
    pub cell_type: Vec<HashMap<String, String>>,
    pub cell_type_specified: Vec<bool>,
    pub pca: Vec<HashMap<String, Vec<f64>>>,
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

#[derive(Clone, Default)]
pub struct ColInfo {
    pub uids: Vec<Option<usize>>,
    pub vids: Vec<usize>,
    pub vpids: Vec<Option<usize>>,
    pub dids: Vec<Option<usize>>,
    pub jids: Vec<usize>,
    pub cids: Vec<Option<usize>>,
    pub fr1_starts: Vec<usize>,
    pub cdr1_starts: Vec<usize>,
    pub fr2_starts: Vec<usize>,
    pub cdr2_starts: Vec<usize>,
    pub fr3_starts: Vec<usize>,
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn justification(x: &str) -> u8 {
    if x == "amino"
        || x == "var"
        || x == "const"
        || (x.ends_with("_aa") && x != "dref_aa")
        || x.ends_with("_dna")
        || x.ends_with("_name")
        || x == "cdiff"
        || x == "notes"
        || x == "edit"
        || x == "datasets"
        || x == "donors"
        || x == "ext"
        || x == "barcode"
        || x == "barcodes"
        || x.starts_with("vj_seq")
        || x.starts_with("vj_seq_nl")
        || x.starts_with("vj_aa_nl")
        || x.starts_with("seq")
        || x.starts_with("q")
        || x.starts_with("var_indices")
        || x.starts_with("share_indices")
        || x.ends_with("_barcode")
        || x.ends_with("_barcodes")
    {
        return b'l';
    } else {
        return b'r';
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define the set "parseable_fields" of fields that could occur in parseable output.
//
// The overlap with code in proc_args_check.rs is not nice.

pub fn set_speakers(ctl: &EncloneControl, parseable_fields: &mut Vec<String>) {
    // Make some abbreviations.

    let lvars = &ctl.clono_print_opt.lvars;

    // Define parseable output columns.  The entire machinery for parseable output is controlled
    // by macros that begin with "speak".

    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    macro_rules! speaker {
        ($var:expr) => {
            if ctl.parseable_opt.pcols.is_empty() || bin_member(&pcols_sort, &$var.to_string()) {
                parseable_fields.push($var.to_string());
            }
        };
    }
    let mut have_gex = false;
    for i in 0..ctl.origin_info.gex_path.len() {
        if ctl.origin_info.gex_path[i].len() > 0 {
            have_gex = true;
        }
    }
    let mut all_lvars = lvars.clone();
    for i in 0..LVARS_ALLOWED.len() {
        let x = &LVARS_ALLOWED[i];
        if !have_gex {
            if *x == "gex".to_string()
                || x.starts_with("gex_")
                || x.ends_with("_g")
                || x.ends_with("_g_μ")
                || *x == "n_gex_cell".to_string()
                || *x == "n_gex".to_string()
                || *x == "n_b".to_string()
                || *x == "clust".to_string()
                || *x == "type".to_string()
                || *x == "entropy".to_string()
                || *x == "cred".to_string()
                || *x == "cred_cell".to_string()
            {
                continue;
            }
        }
        if !lvars.contains(&x.to_string()) {
            all_lvars.push(x.to_string());
        }
    }
    for x in all_lvars.iter() {
        if (*x == "sec" || *x == "mem") && !ctl.gen_opt.using_secmem {
            continue;
        }
        speaker!(x);
    }

    // Define chain variables for parseable output.

    macro_rules! speakerc {
        ($col:expr, $var:expr) => {
            let varc = format!("{}{}", $var, $col + 1);
            if ctl.parseable_opt.pcols.is_empty() || bin_member(&pcols_sort, &varc) {
                parseable_fields.push(format!("{}{}", $var, $col + 1));
            }
        };
    }
    for col in 0..ctl.parseable_opt.pchains {
        for x in CVARS_ALLOWED.iter() {
            speakerc!(col, x);
        }
        if ctl.parseable_opt.pbarcode {
            for x in CVARS_ALLOWED_PCELL.iter() {
                speakerc!(col, x);
            }
        }
        for x in &[
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
            "seq",
            "vj_seq",
            "vj_seq_nl",
            "vj_aa",
            "vj_aa_nl",
            "var_aa",
        ] {
            speakerc!(col, x);
        }
        for i in 0..pcols_sort.len() {
            if pcols_sort[i].starts_with('q') && pcols_sort[i].ends_with(&format!("_{}", col + 1)) {
                let x = pcols_sort[i].after("q").rev_before("_");
                if x.parse::<usize>().is_ok() {
                    parseable_fields.push(pcols_sort[i].clone());
                }
            }
        }
    }

    // Define more lead variables for parseable output.

    speaker!("group_id");
    speaker!("group_ncells");
    speaker!("clonotype_id");
    speaker!("clonotype_ncells");
    speaker!("nchains");
    speaker!("exact_subclonotype_id");
    speaker!("barcodes");
    for x in ctl.origin_info.dataset_list.iter() {
        if x.len() > 0 {
            speaker!(&format!("{}_barcodes", x));
        }
    }
    if ctl.parseable_opt.pbarcode {
        speaker!("barcode");
        for x in ctl.origin_info.dataset_list.iter() {
            speaker!(&format!("{}_barcode", x));
        }
    }
}
