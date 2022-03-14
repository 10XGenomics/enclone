// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test for help request.

use crate::help_utils::{explain_alt_versions, HelpDesk};
use enclone_core::defs::EncloneControl;
use enclone_core::testlist::EXAMPLES;
use itertools::Itertools;
use string_utils::strme;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help2(args: &Vec<String>, _ctl: &EncloneControl, h: &mut HelpDesk) -> Result<(), String> {
    // Set up.

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide example1 help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "example1") || h.help_all {
        h.begin_doc("example1")?;
        h.print("\nShown below is the output of the command:\n")?;
        h.print(&format!("\n\\bold{{enclone {}}}\n", EXAMPLES[0]))?;
        if !h.plain {
            h.print_plain(&format!("{}", include_str!("example1")))?;
        } else {
            let s = include_str!("example1").as_bytes();
            let mut x = Vec::<u8>::new();
            let mut escaped = false;
            for l in 0..s.len() {
                if s[l] == b'' {
                    escaped = true;
                }
                if escaped {
                    if s[l] == b'm' {
                        escaped = false;
                    }
                    continue;
                }
                x.push(s[l]);
            }
            h.print_plain(&strme(&x).to_string())?;
        }
        h.print(
            "This shows an invocation of enclone that takes one dataset as input \
             and exhibits\nall clonotypes for which some chain has the given CDR3 sequence.\n\n\
             \
             What you see here is a compressed view of the entire information encoded in the\n\
             full length transcripts of the 13 cells comprising this clonotype: every base!\n\
             There is a lot to explain about the compression, so please read carefully.\n\n\
             \
             • Clonotypes are grouped.  Here we see just one group having one clonotype in it.\n\
             • This clonotype has three exact subclonotypes in it, the first of which has 10 \
             cells.\n\
             • This clonotype has two chains.  The reference segments for them are shown at \
             the top.\n\
             • The notation 740.1.2 says that this V reference sequence is an alternate allele\n  \
             derived from the universal reference sequence (contig in the reference file)\n  \
             numbered 181, that is from donor 1 (\"740.1\") and is alternate allele 2 for that \
             donor.\n\
             • Sometimes chains are missing from exact subclonotypes.\n\
             • Amino acids are assigned different colors depending on which codon they represent.\n\
             • Numbered columns show the state of particular amino acids, e.g. the first column \
             is for amino\n  acid 20 in chain 1 (where 0 is the start codon).  The numbers read \
             vertically, downward!\n\
             • Universal ref: state for the contig in the reference file.\n\
             • Donor ref: state for the inferred donor germline sequence.\n\
             • ◦s are \"holes\" in the recombined region where the reference doesn't make sense.\n\
             • The \"dot and x\" line has xs where there's a difference *within* the clonotype.\n\
             • Amino acids are shown if they differ from the universal reference or are in \
             the CDR3.\n\
             • u = median UMI count for a chain in the exact subclonotype.\n\
             • const = const region name for a chain in the exact subclonotype.\n\n",
        )?;
        h.print(
            "The view you see here is configurable: see the documentation at \
             \\bold{enclone help lvars} and \\bold{enclone help cvars}.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide example2 help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "example2") || h.help_all {
        h.begin_doc("example2")?;
        h.print("\nShown below is the output of the command:\n")?;

        // Remove H5.

        let ex2_args = EXAMPLES[1].split(' ').collect::<Vec<&str>>();
        let mut ex2_args2 = Vec::<&str>::new();
        for i in 0..ex2_args.len() {
            if ex2_args[i] != "H5" {
                ex2_args2.push(ex2_args[i]);
            }
        }

        // Proceed.

        h.print(&format!(
            "\n\\bold{{enclone {}}}\n",
            ex2_args2.iter().format(" ")
        ))?;
        if !h.plain {
            h.print_plain_unchecked(include_str!("example2"));
        } else {
            let s = include_str!("example2").as_bytes();
            let mut x = Vec::<u8>::new();
            let mut escaped = false;
            for l in 0..s.len() {
                if s[l] == b'' {
                    escaped = true;
                }
                if escaped {
                    if s[l] == b'm' {
                        escaped = false;
                    }
                    continue;
                }
                x.push(s[l]);
            }
            h.print_plain_unchecked(&strme(&x).to_string());
        }
        h.print(
            "This shows an invocation of enclone that takes VDJ, and gene expression \
             data as input, and exhibits all clonotypes for which some chain has the \
             given CDR3 sequence.  As well the command requests UMI (molecule) counts for one \
             hand-selected gene.  You can use any gene(s) you like and any \
             antibodies for which you have feature barcodes.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide input help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "input") || h.help_all {
        h.begin_doc("input")?;
        h.print(
            "\nenclone has \\boldred{two} mechanisms for specifying input datasets: either \
             directly on the command line or via a supplementary metadata file. Only one mechanism \
             may be used at a time.\n\n\
             In both cases, you will need to provide paths to directories where the outputs of \
             the Cell Ranger pipeline may be found.  enclone uses only some of the pipeline \
             output files, so it is enough that those files are present in given directory, and \
             the particular files that are needed may be found by typing \
             \\bold{enclone help input_tech}.\n\n",
        )?;
        h.print_with_box(
            "If you use the argument \\bold{PRE=p} then \\bold{p/} will be prepended to all \
             pipeline paths.  A comma-separated list is also allowed \\bold{PRE=p1,...,pn}, in \
             which case these directories are searched from left to right, until one works, and \
             if all fail, the path is used without prepending anything.  Lastly, \
             (see \\bold{enclone help command}), you can avoid putting \
             \\bold{PRE} on the command line by setting the environment variable \
             \\bold{ENCLONE_PRE} to the desired value.  The default value for \\bold{PRE} \
             is\n\\bold{~/enclone/datasets_me,~/enclone/datasets,~/enclone/datasets2}.  There is \
             also an argument \\bold{PREPOST=x} that causes \\bold{/x} to be appended to all \
             entries in \\bold{PRE}.",
            true,
        )?;
        h.print(
            "Both input forms involve abbreviated names (discussed below), which should be as \
             short as possible, as longer abbreviations will increase the width of the clonotype \
             displays.\n\n",
        )?;
        h.print_with_box(
            "enclone can use gene expression and feature barcode data, as represented by a feature \
             matrix.  Cell Ranger stores this matrix in an hdf5 file, which while generally very \
             efficient, is not optimized for interactive use.  Therefore enclone provides an \
             alternate file structure, which speeds up enclone overall by up to \\boldred{50%}.  \
             To use this, add the argument \\bold{NH5} to the enclone command line.  This will \
             work so long as you have write permission on input directories.  The first time you \
             run enclone (using given inputs), an alternate file feature_barcode_matrix.bin will \
             be written; then subsequent invocations will be faster.  Once the file has been \
             created, it will always be used, regardless of whether \\bold{NH5} is used.  \
             However, we may occasionally change the format of the alternate file.  If do that, \
             then if you have previously generated the file, then it will be rewritten when \
             you invoke enclone for that dataset.  \
             Like with other enclone command-line options, if you want \\bold{NH5} on all the \
             time, you can set the environment variable \\bold{ENCLONE_NH5}.",
            true
        )?;
        h.print(
            "\\boldred{█ 1 █} To point directly at input files on the command line, use e.g.\n\
             \\bold{TCR=/home/jdoe/runs/dataset345}\n\
             or likewise for \\bold{BCR}.  A more complicated syntax is allowed in which commas, \
             colons and semicolons act as delimiters.  Commas go between datasets from the \
             same origin, colons between datasets from the same donor, and semicolons separate \
             donors.  If semicolons are used, the value must be quoted.\n\n",
        )?;
        h.print(
            "enclone uses the distinction between datasets, origins and donors in the following \
             ways:\n\
             1. If two datasets come from the same origin, then enclone can filter to remove \
             certain artifacts, unless you specify the option \\bold{NCROSS}.\n\
             See also illusory clonotype expansion page at \\green{bit.ly/enclone}.\n\
             2. If two cells came from different donors, then enclone will not put them in the \
             same clonotype, unless you specify the option \\bold{MIX_DONORS}.\n\
             More information may be found at `enclone help special`.  In addition, this is \
             enclone's way of keeping datasets organized and affects the output of fields like \
             origin, etc.\n\n",
        )?;

        h.print_with_box(
            "\\bold{Naming.}  Using this input system, each dataset is assigned an abbreviated \
             name, which is \
             everything after the final slash in the directory name (e.g. \\bold{dataset345} in the \
             above example), or the entire name if there is no slash; \
             origins and donors are assigned identifiers s1,... and d1,..., respectively; \
             numbering of origins restarts with each new donor.  \\bold{To specify origins}\n\
             \\bold{and donors, use the second input form, and see in particular} \
             \\green{abbr:path}\\bold{.}",
            true,
        )?;
        h.print(
            "Examples:\n\
             \\bold{TCR=p1,p2}   -- input data from two libraries from the same origin\n\
             \\bold{TCR=p1,p2:q} -- input data as above plus another from a different origin \
             from the same donor\n\
             \\bold{TCR=\"a;b\"}   -- input one library from each of two donors.\n\n",
        )?;
        h.print(
            "Matching gene expression and/or feature barcode data may also be supplied using \
             an argument \\bold{GEX=...}, whose right side must have the exact same structure \
             as the \\bold{TCR} or \\bold{BCR} argument.  Specification of both \
             \\bold{TCR} and \\bold{BCR} is not allowed.  If both BCR and GEX data are in the \
             same directory (from a multi run), and single argument \\bold{BCR_GEX=...} may be \
             used, and similarly one may use \\bold{TCR_GEX}.\n\n",
        )?;
        h.print(
            "In addition, barcode-level data may be specified using \\bold{BC=...}, whose right \
             side is a list of paths having the same structure as the \\bold{TCR} or \\bold{BCR} \
             argument.  Each such path must be for a CSV or TSV file, which must include the field \
             \\bold{barcode}, may include special fields \\bold{origin}, \\bold{donor}, \
             \\bold{tag} and \\bold{color}, and may also include arbitrary other fields.  The \
             \\bold{origin} and \\bold{donor} fields allow a particular origin and donor to be \
             associated to a given barcode.  A use case for this is genetic demultiplexing.  The \
             \\bold{tag} field is intended to be used with tag demultiplexing.  The \\bold{color} \
             field is used by the \\bold{PLOT} option.  All other fields are treated as lead \
             variables, but values are only displayed in \\bold{PER_CELL} mode, or for parseable \
             output using \\bold{PCELL}.  These fields should not include existing lead variable \
             names.  Use of \\bold{BC} automatically turns on the \\bold{MIX_DONORS} option.\n\n",
        )?;
        h.print(
            "Alternatively, an argument \\bold{BC_JOINT=filename} may be specified, where the \
            filename is a CSV or TSV file like those for \\bold{BC=...}, but with an additional \
            field \\bold{dataset}, whose value is an abbreviated dataset name, and which enables \
            the information to be split up to mirror the \
            specification of \\bold{TCR} or \\bold{BCR}.\n\n",
        )?;
        h.print(
            "The argument \\bold{BC=...} or equivalently \\bold{BC_JOINT=filename} may be used \
            on conjunction with\n\\bold{KEEP_CELL_IF=...} (see \"enclone help special\") to \
            restrict the barcodes used by enclone to a specified set.\n\n",
        )?;
        h.print("\\boldred{█ 2 █} To specify a metadata file, use the command line argument\n")?;
        h.print("\\bold{META=filename}\n")?;
        h.print(
            "This file should be a CSV (comma-separated values) file, with one line per cell \
             group.  After the first line, blank lines and lines starting with # are ignored.  \
             There must be a \
             field tcr or bcr, and some other fields are allowed:\n",
        )?;
        h.doc3("\\bold{field}", "\\bold{default}", "\\bold{meaning}");
        h.ldoc3(
            "tcr",
            "(required!)",
            "path to dataset, or abbr:path, where abbr is an abbreviated",
        );
        h.doc3(
            "or bcr",
            "",
            "name for the dataset; exactly one of tcr or bcr must be used",
        );
        h.ldoc3(
            "gex",
            "null",
            "path to GEX dataset, which may include or consist entirely",
        );
        h.doc3("", "", "of FB data");
        h.ldoc3("origin", "s1", "abbreviated name of origin");

        h.ldoc3("donor", "d1", "abbreviated name of donor");
        h.ldoc3pr(
            "color",
            "null",
            "color to associate to this dataset (for \\bold{PLOT} option)",
        );
        h.ldoc3pr("bc", "null", "name of CSV file as in the \\bold{BC} option");
        h.print_tab3()?;
        h.print(
            "\nMultiple \\bold{META} arguments are cumulative and we also allow \
            \\bold{META} to be a comma-separated list of filenames.  In both cases the \
            \\bold{META} files must have identical header lines.  \
            In addition, metadata maybe fully specified on the command line via \
            \\bold{METAX=\"l1;...;ln\"} where the \\bold{li} are the lines that you would \
            otherwise put in the \\bold{META} file.",
        )?;

        h.print(
            "\n\\boldred{█ 3 █} enclone can also read an ancillary CSV file that \
            specifies arbitrary fields that are associated to particular immune receptor \
            sequences.  This is done using \\bold{INFO=path.csv}; \
            The CSV file must have fields \\bold{vj_seq1}, \
            specifying the full heavy or TRB sequence from the beginning of the V segment to the \
            end of the J segment, and \\bold{vj_seq2}, for the light or TRA chain.  The other \
            fields are then made accessible as lvars (see \"enclone help lvars\"), which are \
            populated for any exact subclonotype having exactly two chains (heavy/light or \
            TRB/TRA) that match the data in the CSV file.  By default, one cannot have two \
            lines for the same antibody, however a separate argument\n\\bold{INFO_RESOLVE} may \
            be used to \"pick the first one\".",
        )?;

        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide cvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "cvars") || h.help_all {
        h.begin_doc("cvars")?;

        // Header.

        h.print(
            "\n\\bold{per-chain column options}: These options define per-chain variables, \
             which correspond to columns that appear once for each chain in each clonotype, and \
             have one entry for each exact subclonotype.  Please note that for medians of \
             integers, we actually report the \"rounded median\", the result of rounding the \
             true median up to the nearest integer, so that e.g. 6.5 is rounded up to 7.\n\n",
        )?;
        h.print(
            "See also \"enclone help lvars\" and the inventory of all variables at
            https://10xgenomics.github.io/enclone/pages/auto/inventory.html.\n\n",
        )?;
        h.print(
            "Per-column variables are specified using\n\
             \\bold{CVARS=x1,...,xn}\n\
             where each xi is one of:\n\n",
        )?;

        // Main table entries.

        h.doc(
            "var",
            "bases at positions in chain that vary across the clonotype",
        );
        h.ldocpr(
            "u",
            "\\red{●} VDJ UMI count for each exact subclonotype, median across cells",
        );
        h.docpr(
            "r",
            "\\red{●} VDJ read count for each exact subclonotype, median across cells",
        );
        h.ldoc(
            "edit",
            "a string that defines the edit of the reference V(D)J concatenation versus",
        );
        h.doc2("the contig, from the beginning of the CDR3 to the end of the J segment;");
        h.doc2("this uses a coordinate system in which 0 is the first base of the J ref");
        h.doc2("segment (or the first base of the D ref segment for IGH and TRB); for");
        h.doc2("example D-4:4 denotes the deletion of the last 4 bases of the V segment, ");
        h.doc2("I0:2 denotes an insertion of 2 bases after the V");
        h.doc2("and I0:2•S5 denotes that plus a substitution at position 5; in computing");
        h.doc2("\"edit\", for IGH and TRB, we always test every possible D segment,");
        h.doc2("regardless of whether one is annotated, and pick the best one; for this");
        h.doc2("reason, \"edit\" may be slow");
        h.doc(
            "comp",
            "a measure of CDR3 complexity, which is the total number of S, D and I",
        );
        h.doc2("symbols in \"edit\" as defined above");
        h.doc(
            "cigar",
            "the CIGAR string that defines the edit of the V..J contig sequence versus",
        );
        h.doc2("the universal reference V(D)J concatenation");

        h.ldoc(
            "cdr*_aa",
            "the CDR*_AA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_aa_L_R_ext",
            "the CDR*_AA sequence, with L amino acids added on the left and R amino acids",
        );
        h.doc2("added on the right; either may be negative, denoting trimming instead");
        h.doc2("of extension");
        h.doc(
            "cdr*_aa_north",
            "the CDR*_AA sequence for BCR defined by North B et al. (2011), A new",
        );
        h.doc(
            "",
            "clustering of antibody CDR loop conformations, J Mol Biol 406, 228-256.",
        );
        h.doc("", "cdr1_aa_north = cdr1_aa_3_3_ext for heavy chains");
        h.doc("", "cdr1_aa_north = cdr1_aa for light chains");
        h.doc("", "cdr2_aa_north = cdr2_aa_2_3_ext for heavy chains");
        h.doc("", "cdr2_aa_north = cdr2_aa_1_0_ext for light chains");
        h.doc("", "cdr3_aa_north = cdr3_aa_-1_-1_ext");
        h.doc(
            "cdr*_aa_ref",
            "cdr*_aa, for the universal reference sequence (but not for cdr3)",
        );
        h.doc(
            "cdr*_len",
            "number of amino acids in the CDR* sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_dna",
            "the CDR*_DNA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_dna_ref",
            "same, for the universal reference sequence (but not for cdr3)",
        );
        h.ldoc(
            "cdr3_aa_conx",
            "consensus for CDR3 across the clonotype, showing X for each variant residue",
        );
        h.doc(
            "cdr3_aa_conp",
            "consensus for CDR3 across the clonotype, showing a property symbol whenever",
        );
        h.doc2("two different amino acids are observed, per the following table:");
        h.doc2("--------------------------------------------------------------------");
        h.doc2("asparagine or aspartic acid   B   DN");
        h.doc2("glutamine or glutamic acid    Z   EQ");
        h.doc2("leucine or isoleucine         J   IL");
        h.doc2("negatively charged            -   DE");
        h.doc2("positively charged            +   KRH");
        h.doc2("aliphatic (non-aromatic)      Ψ   VILM");
        h.doc2("small                         π   PGAS");
        h.doc2("aromatic                      Ω   FWYH");
        h.doc2("hydrophobic                   Φ   VILFWYM");
        h.doc2("hydrophilic                   ζ   STHNQEDKR");
        h.doc2("any                           X   ADEFGHIKLMNPQRSTVWY");
        h.doc2("--------------------------------------------------------------------");
        h.doc2("The table is searched top to bottom until a matching class is found.");
        h.doc2("In the special case where every amino acid is shown as a gap (-),");
        h.doc2("a \"g\" is printed.");
        h.ldoc(
            "fwr*_aa",
            "the FWR*_AA sequence, or \"unknown\" if not computed",
        );
        h.doc("fwr*_aa_ref", "same, for the universal reference sequence");
        h.doc(
            "fwr*_len",
            "number of amino acids in the FWR* sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "fwr*_dna",
            "the FWR*_DNA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "fwr*_dna_ref",
            "same, for the universal reference sequences",
        );
        h.doc(
            "",
            "For all of these, * is 1 or 2 or 3 (or 4, for the fwr variables).",
        );
        h.doc(
            "",
            "For CDR1 and CDR2, please see \"enclone help amino\" and the page on",
        );
        h.docpr("", "\\green{bit.ly/enclone} on V(D)J features.");
        h.ldoc("v_start", "start of V segment on full DNA sequence");
        h.doc(
            "d_start",
            "start of D segment on full DNA sequence (or null)",
        );
        h.doc(
            "cdr3_start",
            "base position start of CDR3 sequence on full contig",
        );
        h.doc(
            "d_frame",
            "reading frame of D segment, either 0 or 1 or 2 (or null)",
        );
        h.ldoc(
            "aa%",
            "amino acid percent identity with donor reference, outside junction region",
        );
        h.doc(
            "dna%",
            "nucleotide percent identity with donor reference, outside junction region",
        );

        h.ldoc(
            "v_name_orig",
            "name of V region originally assigned (per cell);",
        );
        h.doc2("values below are clonotype consensuses");
        h.doc("utr_name", "name of 5'-UTR region");
        h.doc("v_name", "name of V region");
        h.doc("d_name", "name of D region (or null)");
        h.doc("j_name", "name of J region");
        h.doc("const", "name of constant region");
        h.ldoc("utr_id", "id of 5'-UTR region");
        h.doc("v_id", "id of V region");
        h.doc("d_id", "id of D region (or null)");
        h.doc("j_id", "id of J region");
        h.doc("const_id", "id of constant region (or null, if not known)");
        h.doc2("(these are the numbers after \">\" in the VDJ reference file)");
        h.ldoc(
            "allele",
            "numerical identifier of the computed donor reference allele",
        );
        h.doc2("for this exact subclonotype");
        h.doc(
            "allele_d",
            "variant bases in the allele for this exact subclonotype,",
        );
        h.doc2("and a list of all the possibilities for this");
        h.ldoc("d1_name", "name of optimal D gene, or none");
        h.doc("d2_name", "name of second best D gene, or none");
        h.doc("d1_score", "score for optimal D gene");
        h.doc("d2_score", "score for second best D gene");
        h.doc(
            "d_delta",
            "score difference between first and second best D gene",
        );
        h.doc("d_Δ", "same");
        h.doc2("These are recomputed from scratch and ignore the given assignment.");
        h.doc2("Note that in many cases D gene assignments are essentially random, as");
        h.doc2("it is often not possible to know the true D gene assignment.");
        h.doc2("If the value is \"null\" it means that having no D gene at all scores better");
        h.ldoc(
            "vjlen",
            "number of bases from the start of the V region to the end of the J region",
        );
        h.doc2("Please note that D gene assignments are frequently \"random\" -- it is not");
        h.doc2("possible to know the actual D gene that was assigned.");
        h.doc(
            "clen",
            "length of observed constant region (usually truncated at primer start)",
        );
        h.doc("ulen", "length of observed 5'-UTR sequence;");
        h.doc(
            "",
            "note however that what report is just the start of the V segment",
        );
        h.doc(
            "",
            "on the contig, and thus the length may include junk before the UTR",
        );
        h.doc(
            "cdiff",
            "differences with universal reference constant region, shown in the",
        );
        h.doc(
            "",
            "abbreviated form e.g. 22T (ref changed to T at base 22) or 22T+10",
        );
        h.doc(
            "",
            "(same but contig has 10 additional bases beyond end of ref C region",
        );
        h.doc(
            "",
            "At most five differences are shown, and if there are more, ... is appended.",
        );
        h.doc("udiff", "like cdiff, but for the 5'-UTR");
        h.ldoc("q<n>_", "comma-separated list of the quality");
        h.doc2("scores at zero-based position n, numbered starting at the");
        h.doc2("beginning of the V segment, for each cell in the exact subclonotype");
        h.ldoc(
            "notes",
            "optional note if there is an insertion or the end of J does not exactly abut",
        );
        h.doc(
            "",
            "the beginning of C; elided if empty; also single base overlaps between",
        );
        h.docpr(
            "",
            "J and C are not shown unless you use the special option \\bold{JC1}; we do this",
        );
        h.doc(
            "",
            "because with some VDJ references, one nearly always has such an overlap",
        );
        h.ldoc(
            "ndiff<n>vj",
            "number of base differences within V..J between this exact subclonotype and",
        );
        h.doc("", "exact subclonotype n");
        h.doc(
            "d_univ",
            "distance from universal reference, more specifically,",
        );
        h.doc(
            "",
            "number of base differences within V..J between this exact",
        );
        h.doc(
            "",
            "clonotype and universal reference, exclusive of indels, the last 15",
        );
        h.doc("", "bases of the V and the first 15 bases of the J");
        h.doc("d_donor", "distance from donor reference,");
        h.doc("", "as above but computed using donor reference");

        // The rest.

        h.print_tab2()?;
        h.print("\n")?;
        explain_alt_versions(h)?;
        h.print(
            "\nAt least one variable must be listed.  The default is \\bold{u,const,notes}.  \
             \\bold{CVARSP}: same as \\bold{CVARS} but appends.\n\n",
        )?;
        h.end_doc();
    }
    Ok(())
}
