// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::cell_color::CellColor;
use crate::linear_condition::LinearCondition;
use debruijn::dna_string::DnaString;
use evalexpr::Node;
use hdf5::Dataset;
use io_utils::{open_for_read, path_exists};
use mirror_sparse_matrix::MirrorSparseMatrix;
use perf_stats::{elapsed, peak_mem_usage_gb};
use regex::Regex;
use std::cmp::max;
use std::collections::HashMap;

use std::io::BufRead;
use std::sync::atomic::AtomicBool;
use std::time::{Instant, SystemTime};
use string_utils::TextUtils;
use vector_utils::unique_sort;

pub static FAILED: AtomicBool = AtomicBool::new(false);

pub const HELP_PAGES: [&str; 21] = [
    "all",
    "amino",
    "color",
    "command",
    "cvars",
    "display",
    "example1",
    "example2",
    "faq",
    "filter",
    "glossary",
    "how",
    "indels",
    "input",
    "input_tech",
    "lvars",
    "main",
    "parseable",
    "quick",
    "setup",
    "special",
];

pub const MAX_CDR3_DIFFS_TO_JOIN: usize = 5;

// Clonotyping algorithm heuristics.

#[derive(Default, PartialEq)]
pub struct ClonotypeHeuristics {
    pub max_diffs: usize,
    pub max_degradation: usize,
    pub ref_v_trim: usize,
    pub ref_j_trim: usize,
}

// Origin info data structure.

#[derive(Default, PartialEq)]
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

#[derive(Default, PartialEq)]
pub struct GeneralOpt {
    pub pre: Vec<String>,
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
    pub tcrgd: bool,
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
    pub noprintx: bool,
    pub required_fps: Option<usize>,
    pub required_cells: Option<usize>,
    pub required_clonotypes: Option<usize>,
    pub required_donors: Option<usize>,
    pub required_two_cell_clonotypes: Option<usize>,
    pub required_two_chain_clonotypes: Option<usize>,
    pub required_three_chain_clonotypes: Option<usize>,
    pub required_four_chain_clonotypes: Option<usize>,
    pub required_datasets: Option<usize>,
    pub cellranger: bool,
    pub summary: bool,
    pub summary_clean: bool,
    pub summary_csv: bool,
    pub cr_version: String,
    pub nwarn: bool,
    pub gene_scan_test: Option<LinearCondition>,
    pub gene_scan_control: Option<LinearCondition>,
    pub gene_scan_threshold: Option<LinearCondition>,
    pub gene_scan_exact: bool,
    pub clonotype_group_names: Option<String>,
    pub origin_color_map: HashMap<String, String>,
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
    pub tree_on: bool,
    pub tree: Vec<String>,
    pub allow_inconsistent: bool,
    pub color: String,
    pub color_by_rarity_pc: f64,
    pub species: String, // human or mouse or unknown, determined from the reference sequence
    pub using_secmem: bool,
    pub diff_style: String,
    pub accept_broken: bool,
    pub require_unbroken_ok: bool,
    pub built_in: bool,
    pub reprod: bool,
    pub peer_group_filename: String,
    pub peer_group_dist: String,
    pub peer_group_readable: bool,
    pub subset_json: String,
    pub fold_headers: bool,
    pub no_uncap_sim: bool,
    pub profile: bool,
    pub nopager: bool,
    pub info: Option<String>,
    pub info_fields: Vec<String>,
    pub info_data: HashMap<String, Vec<String>>,
    pub info_resolve: bool,
    pub internal_data_dir: String,
    pub row_fill_verbose: bool,
    pub config_file: String,
    pub config: HashMap<String, String>,
    pub top_genes: bool,
    pub toy_com: bool,
    pub chains_to_align: Vec<usize>,
    pub chains_to_align2: Vec<usize>,
    pub chains_to_jun_align: Vec<usize>,
    pub chains_to_jun_align2: Vec<usize>,
    pub align_jun_align_consistency: bool,
    pub dvars: Vec<String>, // per dataset variables
    pub gvars: Vec<String>, // per run variables
    pub jscore_match: i32,
    pub jscore_mismatch: i32,
    pub jscore_bits_multiplier: f64,
    pub jscore_gap_open: i32,
    pub jscore_gap_extend: i32,
    pub split: bool,
    pub max_heavies: usize,
    pub cpu_all_start: usize,
    pub cpu_this_start: usize,
    pub evil_eye: bool, // extra printing to try to trace hangs
    pub toy: bool,      // toy with phylogeny
    pub group_post_filter: Option<Vec<usize>>,
    pub no_newline: bool,
    pub fb_show: String,
    pub var_def: Vec<(String, String, Node, String)>, // {(variable, value, compiled value, expr)}
    pub nospaces: bool,
    pub subsample: f64,
    pub all_bc_filename: String,
    pub all_bc_human: bool,
    pub all_bc_fields: Vec<String>,
    pub all_bc_fields_orig: Vec<String>,
    pub gamma_delta: bool,
    pub pre_eval: bool,
}

// Some plot options.  Note that plot options are not allowed to affect intermediate computation.

#[derive(Clone, Default)]
pub struct PlotOpt {
    pub cell_color: CellColor,
    pub plot_xy_filename: String,
    pub plot_xy_xvar: String,
    pub plot_xy_yvar: String,
    pub plot_xy_x_log10: bool,
    pub plot_xy_y_log10: bool,
    pub plot_xy_sym: bool,
    pub plot_conditions: Vec<String>,
    pub plot_colors: Vec<String>,
    pub plot_file: String,
    pub plot_by_isotype: bool,
    pub plot_by_isotype_nolegend: bool,
    pub plot_by_isotype_color: Vec<String>,
    pub plot_by_mark: bool,
    pub plot_quad: bool,
    pub use_legend: bool,
    pub legend: Vec<(String, String)>,
    pub sim_mat_plot_file: String,
    pub sim_mat_plot_vars: Vec<String>,
    pub honey_in: Option<String>,
    pub honey_out: String,
    pub split_plot_by_origin: bool,
    pub png_width: Option<usize>,
}

// Allele-finding algorithmic options.

#[derive(Default, PartialEq)]
pub struct AlleleAlgOpt {
    pub min_mult: usize,
    pub min_alt: usize,
}

// Allele-finding print options.

#[derive(Default, PartialEq)]
pub struct AllelePrintOpt {
    pub con: bool,       // print alternate consensus sequences
    pub con_trace: bool, // tracing for con
}

// Join printing options.

#[derive(Default, PartialEq)]
pub struct JoinPrintOpt {
    pub seq: bool,     // print sequences of contigs, before truncation to V..J
    pub ann: bool,     // print annotations of contigs
    pub ann0: bool,    // print annotations of contigs, after truncation to V..J
    pub show_bc: bool, // show barcodes
    pub quiet: bool,   // do not print join events
    pub pfreq: usize,  // show data for 1/n joins even if correct
}

// Join algorithmic options.

#[derive(Default, PartialEq)]
pub struct JoinAlgOpt {
    pub max_score: f64,          // max score for join
    pub easy: bool,              // make joins even if core condition violated
    pub merge_onesies: bool,     // create and merge onesies where completely unambiguous
    pub merge_onesies_ctl: bool, // restriction on onesie merger
    pub bcjoin: bool,            // join only by barcode identity
    pub max_cdr3_diffs: usize,
    pub cdr3_mult: f64, // multiplier for checking CDR3 SHM concentration
    pub old_mult: bool,
    pub mult_pow: f64,
    pub old_light: bool,
    pub basic_h: Option<f64>,
    pub basic: Option<f64>,
    pub basicx: bool,
    pub join_full_diff: bool,
    pub join_cdr3_ident: f64,
    pub cdr3_normal_len: usize,
    pub auto_share: usize,
}

// Clonotype filtering options.
// These fall into 2 categories: 1) on by default and 2) user-specified.
// Note that ClonoFiltOpt options are not allowed to affect intermediate computation.

#[derive(Default, PartialEq)]
pub struct ClonoFiltOptDefault {
    pub marked_b: bool, // only print clonotypes having a mark and which are typed as B cells
    pub donor: bool,    // allow cells from different donors to be placed in the same clonotype
    pub weak_foursies: bool, // filter weak foursies
    pub ngex: bool,     // turn off gex filtering,
    pub non_cell_mark: bool,
    pub weak_onesies: bool,        // filter weak onesies
    pub doublet: bool,             // filter putative doublets
    pub fcell: Vec<Node>,          // constraints from FCELL
    pub umi_filt: bool,            // umi count filter
    pub umi_filt_mark: bool,       // umi count filter (but only mark)
    pub umi_ratio_filt: bool,      // umi ratio filter
    pub umi_ratio_filt_mark: bool, // umi ratio filter (but only mark)
    pub weak_chains: bool,         // filter weak chains from clonotypes
    pub whitef: bool,              // only show clonotypes exhibiting whitelist contamination
    pub ncross: bool,              // turn off cross filtering,
    pub bc_dup: bool,              // filter duplicated barcodes within an exact subclonotype
    pub signature: bool,           // signature filtering
    pub nmax: bool,                // turn off max contigs filter
}

#[derive(Default)]
pub struct ClonoFiltOpt {
    pub ncells_low: usize,   // only show clonotypes with at least this many cells
    pub ncells_high: usize,  // only show clonotypes with at most this many cells
    pub min_umi: usize,      // only show clonotypes with at least this many UMIs in some contig
    pub min_datasets: usize, // only show clonotypes involving at least this many datasets
    pub max_datasets: usize, // only show clonotypes involving at most this many datasets
    pub min_dataset_ratio: usize, // see "enclone help filter"
    pub min_donors: usize,   // only show clonotypes having at least this many donors
    pub min_origins: usize,  // only show clonotypes involving at least this many origins
    pub min_chains: usize,   // only show clonotypes with at least this many chains
    pub max_chains: usize,   // only show clonotypes with at most this many chains
    pub cdr3: Option<Regex>, // only show clonotypes whose CDR3_AA matches regular expression
    pub cdr3_lev: String,    // only show clonotypes whose CDR3_AA matches Levenshtein dist pattern
    pub protect_bads: bool,  // protect bads from deletion
    pub fail_only: bool,     // only print fails
    pub seg: Vec<Vec<String>>, // only show clonotypes using one of these VDJ segment names
    pub segn: Vec<Vec<String>>, // only show clonotypes using one of these VDJ segment numbers
    pub nseg: Vec<Vec<String>>, // do not show clonotypes using one of these VDJ segment names
    pub nsegn: Vec<Vec<String>>, // do not show clonotypes using one of these VDJ segment numbers
    pub min_exacts: usize,   // only show clonotypes having at least this many exact subclonotypes
    pub max_exacts: usize,
    pub vj: Vec<u8>, // only show clonotypes having exactly this full length V..J sequence
    pub vdup: bool,  // only show clonotypes having a same V segment in two chains
    pub have_onesie: bool, // only show clonotypes including a onesie exact subclonotype
    pub cdiff: bool, // only show clonotypes having a constant region difference
    pub del: bool,   // only show clonotypes exhibiting a deletion
    pub qual_filter: bool, // filter out exact subclonotypes having a weak base
    pub bounds: Vec<LinearCondition>, // bounds on certain variables
    pub bound_type: Vec<String>, // types of those bounds
    pub barcode: Vec<String>, // requires one of these barcodes
    pub marked: bool, // only print clonotypes having a mark
    pub inkt: bool,
    pub mait: bool,
    pub d_inconsistent: bool,
    pub d_none: bool,
    pub d_second: bool,
    pub const_igh: Option<Regex>,
    pub const_igkl: Option<Regex>,
}

// Clonotype printing options.

#[derive(Default, PartialEq)]
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
    pub conx: bool,
    pub conp: bool,
}

// Clonotype grouping options.

#[derive(Default, PartialEq)]
pub struct ClonoGroupOpt {
    // SYMMETRIC AND ASYMMETRIC
    pub ngroup: bool,            // do not print group headers
    pub min_group: usize,        // minimum number of clonotypes in group to print
    pub min_group_donors: usize, // minimum number of donors in a group to print
    pub style: String,           // symmetric or unsymmetric or unspecified
    // SYMMETRIC GROUPING CONTROLS
    pub vj_refname: bool,        // group by having the same VJ reference names
    pub vj_heavy_refname: bool,  // group by having the same heavy VJ reference names
    pub vdj_refname: bool,       // group by having the same VDJ reference names
    pub vdj_heavy_refname: bool, // group by having the same heavy VDJ reference names
    pub vj_len: bool,            // group by V..J of same length
    pub cdr3_len: bool,          // group by CDR3 of same length
    pub cdr3_aa_heavy_pc: Option<f64>, // group if CDR3 aa identity >= given percent on heavy chain
    pub cdr3_aa_light_pc: Option<f64>, // group if CDR3 aa identity >= given percent on light chain
    pub aa_heavy_pc: Option<f64>, // group if amino acid identity >= given percent on heavy chain
    pub aa_light_pc: Option<f64>, // group if amino acid identity >= given percent on light chain
    // ASYMMETRIC GROUPING CONTROLS
    pub asymmetric_center: String, // definition of center for asymmetric grouping
    pub asymmetric_dist_formula: String, // definition of distance formula for asymmetric grouping
    pub asymmetric_dist_bound: String, // definition of distance bound for asymmetric grouping
    // DEPRECATED
    pub vj_refname_strong: bool, // group by having the same VJ reference names, but stronger
}

// Parseable output options.

#[derive(Default, PartialEq)]
pub struct ParseableOpt {
    pub pout: String,             // name of parseable output file
    pub pchains: String,          // number of chains to show in parseable output
    pub pcols: Vec<String>,       // column names to show in parseable output
    pub pcols_sort: Vec<String>,  // sorted column names to show in parseable output
    pub pcols_sortx: Vec<String>, // same but before colon if present
    pub pbarcode: bool,           // generate output per barcode rather than per exact subclonotype
}

// Computational performance options.

#[derive(Default, PartialEq)]
pub struct PerfOpt {
    pub comp: bool,         // print computational performance stats
    pub comp2: bool,        // print more detailed computational performance stats
    pub unaccounted: bool,  // show unaccounted time at each step
    pub comp_enforce: bool, // comp plus enforce no unaccounted time
}

// Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
// run of enclone.  If you add something to this, be sure to update the "changed" section in
// enclone_server.rs, if needed.

#[derive(Default)]
pub struct EncloneControl {
    pub visual_mode: bool,                       // running as enclone visual
    pub perf_opt: PerfOpt,                       // computational performance options
    pub start_time: Option<Instant>,             // enclone start time
    pub gen_opt: GeneralOpt,                     // miscellaneous general options
    pub plot_opt: PlotOpt,                       // plot options
    pub pretty: bool,                            // use escape characters to enhance view
    pub silent: bool,                            // turn off extra logging
    pub force: bool,                             // make joins even if redundant
    pub debug_table_printing: bool,              // turn on debugging for table printing
    pub merge_all_impropers: bool,               // merge all improper exact subclonotypes
    pub heur: ClonotypeHeuristics,               // algorithmic heuristics
    pub origin_info: OriginInfo,                 // origin (sample) info
    pub allele_alg_opt: AlleleAlgOpt,            // algorithmic options for allele finding
    pub allele_print_opt: AllelePrintOpt,        // print options for allele finding
    pub join_alg_opt: JoinAlgOpt,                // algorithmic options for join
    pub join_print_opt: JoinPrintOpt,            // printing options for join operations
    pub clono_filt_opt_def: ClonoFiltOptDefault, // default filtering options for clonotypes
    pub clono_filt_opt: ClonoFiltOpt,            // filtering options for clonotypes
    pub clono_print_opt: ClonoPrintOpt,          // printing options for clonotypes
    pub clono_group_opt: ClonoGroupOpt,          // grouping options for clonotypes
    pub parseable_opt: ParseableOpt,             // parseable output options
    pub pathlist: Vec<String>,                   // list of input files
    pub last_modified: Vec<SystemTime>,          // last modified for pathlist
}

pub static mut WALLCLOCK: f64 = 0.0;
pub static mut LAST_IPEAK: f64 = -0.0;

impl EncloneControl {
    pub fn perf_stats(&self, t: &Instant, msg: &str) {
        let used = elapsed(t);
        let t2 = Instant::now();
        let mut usedx = String::new();
        if self.perf_opt.comp {
            let peak = peak_mem_usage_gb();
            let ipeak = (100.0 * peak).round();
            let peak_mem = format!("peak mem = {:.2} GB", peak);
            usedx = format!("{:.2}", used);
            let mut ipeak_changed = false;
            unsafe {
                if ipeak != LAST_IPEAK {
                    ipeak_changed = true;
                    LAST_IPEAK = ipeak;
                }
            }
            if usedx != "0.00" || ipeak_changed {
                println!("used {} seconds {}, {}", usedx, msg, peak_mem);
            }
        }

        // Check for time used in the above computation, which could otherwise introduce a
        // discrepancy into the time accounting stats.  Surprisingly, the time spent in that
        // section can be nontrivial.

        let used2 = elapsed(&t2);
        let used2x = format!("{:.2}", used2);
        if self.perf_opt.comp && used2x != "0.00" {
            println!("used {} seconds computing perf stats for {}", used2x, msg);
        }

        // Update total time used.

        unsafe {
            WALLCLOCK += used + used2;
        }

        // Report unaccounted time.

        if self.perf_opt.comp && self.perf_opt.unaccounted && msg != "total" {
            let delta;
            unsafe {
                delta = elapsed(&self.start_time.unwrap()) - WALLCLOCK;
            }
            let deltas = format!("{:.2}", delta);
            if deltas != "0.00" {
                if usedx == "0.00" {
                    println!("used 0.00 seconds {}", msg);
                }
                println!("used {} seconds unaccounted for", deltas);
            }
        }
    }
}

// Set up data structure to track clonotype data.  A TigData is for one contig;
// a Vec<TigData> is for one barcode, and an ExactClonotype is for an exact subclonotype.

#[derive(Eq, Ord, PartialEq, PartialOrd, Default, Clone)] // not sure these are all needed
pub struct TigData {
    pub cdr3_dna: String,                        // CDR3 DNA sequence
    pub len: usize,                              // length of V..J sequence
    pub v_start: usize,                          // start of V on full contig sequence
    pub v_stop: usize,                           // stop of aligned V on full contig sequence
    pub v_stop_ref: usize,                       // stop of aligned V on reference V
    pub d_start: Option<usize>,                  // start of aligned D on full contig sequence
    pub j_start: usize,                          // start of aligned J on full contig sequence
    pub j_start_ref: usize,                      // start of aligned J on reference J
    pub j_stop: usize,                           // stop of J on full contig sequence
    pub c_start: Option<usize>,                  // start of C on full contig sequence
    pub full_seq: Vec<u8>,                       // full contig sequence
    pub u_ref_id: Option<usize>,                 // index of 5'-UTR in ref file if found
    pub v_ref_id: usize, // index of V segment reference sequence in ref file
    pub d_ref_id: Option<usize>, // index of D segment reference sequence in ref file
    pub j_ref_id: usize, // index of J segment reference sequence in ref file
    pub c_ref_id: Option<usize>, // index of C segment reference sequence in ref file
    pub fr1_start: usize, // start position in bases of FWR1 on V..J
    pub cdr1_start: Option<usize>, // start position in bases of CDR1 on V..J
    pub fr2_start: Option<usize>, // start position in bases of FWR2 on V..J
    pub cdr2_start: Option<usize>, // start position in bases of CDR2 on V..J
    pub fr3_start: Option<usize>, // start position in bases of FWR3 on V..J
    pub cdr3_aa: String, // CDR3 amino acid sequence
    pub cdr3_start: usize, // start position in bases of CDR3 on V..J
    pub quals: Vec<u8>,  // quality scores, truncated to V..J
    pub full_quals: Vec<u8>, // quality scores
    pub barcode: String, // barcode
    pub tigname: String, // name of contig
    pub left: bool,      // true if this is IGH or TRB (or TRD in gamma/delta mode)
    pub dataset_index: usize, // index of dataset
    pub origin_index: Option<usize>, // index of origin (sample)
    pub donor_index: Option<usize>, // index of donor
    pub tag_index: Option<usize>, // index of tag
    pub umi_count: usize, // number of UMIs supporting contig
    pub read_count: usize, // number of reads supporting contig
    pub chain_type: String, // e.g. IGH
    pub annv: Vec<(i32, i32, i32, i32, i32)>, // V annotation (one or two entries), for V..J
    pub validated_umis: Option<Vec<String>>, // validated UMIs
    pub non_validated_umis: Option<Vec<String>>, // non-validated UMIs
    pub invalidated_umis: Option<Vec<String>>, // invalidated UMIs
    pub frac_reads_used: Option<u32>, // fraction of reads passed to assembly stage in CR
}

impl TigData {
    pub fn seq(&self) -> &[u8] {
        &self.full_seq[self.v_start..self.j_stop]
    }
}

// The ExactClonotype data structure stores information that could be exhibited as a
// Vec<Vec<Vec<TigData>>>, but it avoids repetition of identical data.
//
// TigData0: data for each cell
// TigData1: shared data

#[derive(Clone)]
pub struct TigData0 {
    pub quals: Vec<u8>,                          // quality scores, truncated to V..J
    pub v_start: usize,                          // start of V on full contig sequence
    pub j_stop: usize,                           // stop of J on full contig sequence
    pub c_start: Option<usize>,                  // start of C on full contig sequence
    pub full_seq: Vec<u8>,                       // full contig sequence
    pub barcode: String,                         // barcode
    pub tigname: String,                         // name of contig
    pub dataset_index: usize,                    // index of dataset
    pub origin_index: Option<usize>,             // index of origin (sample)
    pub donor_index: Option<usize>,              // index of donor
    pub tag_index: Option<usize>,                // index of tag
    pub umi_count: usize,                        // number of UMIs supporting contig
    pub read_count: usize,                       // number of reads supporting contig
    pub marked: bool,                            // if marked for possible deletion
    pub validated_umis: Option<Vec<String>>,     // validated UMIs
    pub non_validated_umis: Option<Vec<String>>, // non-validated UMIs
    pub invalidated_umis: Option<Vec<String>>,   // invalidated UMIs
    pub frac_reads_used: Option<u32>,            // fraction of reads passed to assembly stage in CR
}

#[derive(Clone)]
pub struct TigData1 {
    pub cdr3_dna: String,           // CDR3 DNA sequence
    pub seq: Vec<u8>,               // V..J contig subsequence
    pub seq_del: Vec<u8>,           // V..J, possibly with mod 3 del
    pub seq_del_amino: Vec<u8>,     // V..J, possibly with mod 3 del at mod 3 start
    pub aa_mod_indel: Vec<u8>,      // amino acid sequence, after removing indel if present
    pub ins: Vec<(usize, Vec<u8>)>, // insertions in V..J (currently at most one) = {(pos, seq)}
    // **before** the given position
    pub full_seq: Vec<u8>,                   // full contig sequence (consensus)
    pub v_start: usize,                      // start of V on full contig sequence
    pub v_stop: usize,                       // stop of aligned V on full contig sequence
    pub v_stop_ref: usize,                   // stop of aligned V on reference V
    pub d_start: Option<usize>,              // start of aligned D on full contig sequence
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
    pub cdr1_start: Option<usize>,           // start position in bases of CDR1 on V..J
    pub fr2_start: Option<usize>,            // start position in bases of FWR2 on V..J
    pub cdr2_start: Option<usize>,           // start position in bases of CDR2 on V..J
    pub fr3_start: Option<usize>,            // start position in bases of FWR3 on V..J
    pub cdr3_aa: String,                     // CDR3 amino acid sequence
    pub cdr3_start: usize,                   // start position in bases of CDR3 on V..J
    pub left: bool,         // true if this is IGH or TRB (or TRD in gamma/delta mode)
    pub chain_type: String, // e.g. IGH
    pub annv: Vec<(i32, i32, i32, i32, i32)>, // V annotation (one or two entries), for V..J
    pub vs: DnaString,      // reference V segment (possibly donor allele)
    pub vs_notesx: String,  // notes on reference V segment (probably to be replaced)
    pub js: DnaString,      // reference J segment
    pub inkt_alpha_chain_gene_match: bool,
    pub inkt_alpha_chain_junction_match: bool,
    pub inkt_beta_chain_gene_match: bool,
    pub inkt_beta_chain_junction_match: bool,
    pub mait_alpha_chain_gene_match: bool,
    pub mait_alpha_chain_junction_match: bool,
    pub mait_beta_chain_gene_match: bool,
    pub mait_beta_chain_junction_match: bool,
}

impl TigData1 {
    pub fn ins_len(&self) -> usize {
        let mut x = 0;
        for j in 0..self.ins.len() {
            x += self.ins[j].1.len();
        }
        x
    }
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
// encapsulating the same info is legacy and should be cleaned up.
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
    pub fb_top_matrices: Vec<MirrorSparseMatrix>,
    pub fb_top_barcodes: Vec<Vec<String>>,
    pub fb_top_reads_matrices: Vec<MirrorSparseMatrix>,
    pub fb_top_reads_barcodes: Vec<Vec<String>>,
    pub fb_total_umis: Vec<u64>,
    pub fb_total_reads: Vec<u64>,
    pub fb_brn: Vec<Vec<(String, u32, u32)>>,
    pub fb_brnr: Vec<Vec<(String, u32, u32)>>,
    pub fb_bdcs: Vec<Vec<(String, u32, u32, u32)>>,
    pub feature_refs: Vec<String>,
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
    pub feature_metrics: Vec<HashMap<(String, String), String>>,
    pub json_metrics: Vec<HashMap<String, f64>>,
    pub metrics: Vec<String>,
}

// Every entry in a ColInfo is a vector whose number of entries is the number of chains
// in a clonotype.

#[derive(Clone, Default)]
pub struct ColInfo {
    pub left: Vec<bool>,
    pub uids: Vec<Option<usize>>,
    pub vids: Vec<usize>,
    pub vpids: Vec<Option<usize>>,
    pub dids: Vec<Option<usize>>,
    pub jids: Vec<usize>,
    pub cids: Vec<Option<usize>>,
    pub fr1_starts: Vec<usize>,
    pub cdr1_starts: Vec<Option<usize>>,
    pub fr2_starts: Vec<Option<usize>>,
    pub cdr2_starts: Vec<Option<usize>>,
    pub fr3_starts: Vec<Option<usize>>,
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
        || x.ends_with("_indices")
        || x == "cdiff"
        || x == "notes"
        || x == "edit"
        || x == "datasets"
        || x == "donors"
        || x == "origins"
        || x == "ext"
        || x == "barcode"
        || x == "barcodes"
        || x == "filter"
        || x.starts_with("vj_seq")
        || x.starts_with("vj_seq_nl")
        || x.starts_with("vj_aa_nl")
        || x.starts_with("seq")
        || x.starts_with('q')
        || x.ends_with("_barcode")
        || x.ends_with("_barcodes")
        || (x.starts_with("cdr") && !x.ends_with("len"))
        || (x.starts_with("fwr") && !x.ends_with("len"))
        || x.starts_with("d1_name")
        || x.starts_with("d2_name")
        || x.starts_with("fb") && !x.ends_with("_n")
        || x == "cigar"
        || x.contains("valumis")
        || x.contains("valbcumis")
        || x == "nbc"
    {
        b'l'
    } else {
        b'r'
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// The POUT separator character used to be a semicolon, but because semicolons could appear in the
// fields, that was broken (and there is a test for the associated problem).  We substituted a
// character (the bell character) that would not be generated by enclone and should not be allowed
// in input, although we do not enforce that currently.

pub const POUT_SEP: &str = "\x07";

// Potential join structure.

#[derive(Default)]
pub struct PotentialJoin {
    pub k1: usize,
    pub k2: usize,
    pub nrefs: usize,
    pub cd: isize,
    pub diffs: usize,
    pub bcs1: Vec<String>,
    pub bcs2: Vec<String>,
    pub shares: Vec<isize>,
    pub indeps: Vec<isize>,
    pub shares_details: Vec<Vec<usize>>,
    pub share_pos_v: Vec<Vec<usize>>,
    pub share_pos_j: Vec<Vec<usize>>,
    pub score: f64,
    pub err: bool,
    pub p1: f64,
    pub mult: f64,
    pub k: isize,
    pub d: isize,
    pub n: usize,
}

pub fn get_config(config_file: &str, config: &mut HashMap<String, String>) -> bool {
    if !config_file.is_empty() {
        let mut cf = config_file.to_string();
        if cf.contains(':') {
            cf = cf.after(":").to_string();
        }
        if path_exists(&cf) {
            let f = open_for_read![&cf];
            for line in f.lines() {
                let s = line.unwrap();
                config.insert(s.before("=").to_string(), s.after("=").to_string());
            }
            return true;
        }
    }
    false
}
