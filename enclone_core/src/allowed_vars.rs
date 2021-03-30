// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Field (variable) names.
// Lead variables for exact subclonotypes and cells.
pub const LVARS_ALLOWED: [&str; 31] = [
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
    "filter",
    "nchains",
];

// Chain variables that can be used for contigs and chains

pub const CVARS_ALLOWED: [&str; 84] = [
    "var",
    "u",
    "u_min",
    "u_max",
    "u_Σ",
    "u_μ",
    "comp",
    "edit",
    "r",
    "r_min",
    "r_max",
    "r_Σ",
    "r_μ",
    "const",
    "white",
    "cdr1_dna",
    "cdr1_dna_ref",
    "cdr2_dna",
    "cdr2_dna_ref",
    "cdr3_dna",
    "cdr1_len",
    "cdr2_len",
    "cdr3_len",
    "cdr1_aa",
    "cdr1_aa_north",
    "cdr1_aa_ref",
    "cdr2_aa",
    "cdr2_aa_north",
    "cdr2_aa_ref",
    "cdr3_aa",
    "cdr3_aa_north",
    "cdr3_aa_conx",
    "cdr3_aa_conp",
    "fwr1_dna",
    "fwr1_dna_ref",
    "fwr2_dna",
    "fwr2_dna_ref",
    "fwr3_dna",
    "fwr3_dna_ref",
    "fwr4_dna",
    "fwr4_dna_ref",
    "fwr1_len",
    "fwr2_len",
    "fwr3_len",
    "fwr4_len",
    "fwr1_aa",
    "fwr1_aa_ref",
    "fwr2_aa",
    "fwr2_aa_ref",
    "fwr3_aa",
    "fwr3_aa_ref",
    "fwr4_aa",
    "fwr4_aa_ref",
    "ulen",
    "vjlen",
    "clen",
    "cdiff",
    "udiff",
    "notes",
    "d_univ",
    "d_donor",
    "aa%",
    "dna%",
    "nval",
    "nnval",
    "valumis",
    "nvalumis",
    "ivalumis",
    "valbcumis",
    "nvalbcumis",
    "ivalbcumis",
    "d_frame",
    "d_start",
    "v_name",
    "d_name",
    "j_name",
    "v_id",
    "d_id",
    "j_id",
    "const_id",
    "utr_id",
    "utr_name",
    "cdr3_start",
    "v_start",
];

pub const CVARS_ALLOWED_PCELL: [&str; 2] = ["u_cell", "r_cell"];

pub const PLVARS_ALLOWED: [&str; 6] = [
    "group_id",
    "group_ncells",
    "clonotype_id",
    "clonotype_ncells",
    "exact_subclonotype_id",
    "barcodes",
];

pub const PCVARS_ALLOWED: [&str; 11] = [
    "var_indices_dna",
    "var_indices_aa",
    "share_indices_dna",
    "share_indices_aa",
    "cdr3_aa",
    "seq",
    "vj_seq",
    "vj_seq_nl",
    "vj_aa",
    "vj_aa_nl",
    "var_aa",
];
