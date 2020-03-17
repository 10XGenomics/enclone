// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::*;
use ansi_escape::*;
use string_utils::*;
use tables::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help4(args: &Vec<String>) {
    // Set up.

    let mut args = args.clone();
    let mut rows = Vec::<Vec<String>>::new();

    macro_rules! doc {
        ($n1:expr, $n2:expr) => {
            rows.push(vec![$n1.to_string(), $n2.to_string()]);
        };
    }
    macro_rules! docpr {
        ($n1:expr, $n2:expr) => {
            rows.push(vec![print_to($n1), print_to($n2)]);
        };
    }
    macro_rules! ldoc {
        ($n1:expr, $n2:expr) => {
            rows.push(vec!["\\hline".to_string(); 2]);
            rows.push(vec![$n1.to_string(), $n2.to_string()]);
        };
    }
    macro_rules! ldocpr {
        ($n1:expr, $n2:expr) => {
            rows.push(vec!["\\hline".to_string(); 2]);
            rows.push(vec![print_to($n1), print_to($n2)]);
        };
    }
    let mut plain = false;
    for i in 0..args.len() {
        if args[i] == "PLAIN" {
            args.remove(i);
            plain = true;
            unsafe {
                PLAIN = true;
            }
            break;
        }
    }
    if args.len() == 1 || (args.len() >= 2 && args[1] == "help") {
        let mut to_delete = vec![false; args.len()];
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                to_delete[i] = true;
            }
        }
        erase_if(&mut args, &to_delete);
    }
    /*
    macro_rules! doc_red {
        ($n1:expr, $n2:expr) => {
            if !plain {
                let r1 = format!( "[01;31m{}[0m", $n1 );
                let r2 = format!( "[01;31m{}[0m", $n2 );
                rows.push( vec![ r1, r2 ] );
            } else {
        };
    }
    macro_rules! ldoc_red {
        ($n1:expr, $n2:expr) => {
            rows.push( vec![ "\\hline".to_string(); 2 ] );
            if !plain {
                let r1 = format!( "[01;31m{}[0m", $n1 );
                let r2 = format!( "[01;31m{}[0m", $n2 );
                rows.push( vec![ r1, r2 ] );
            } else {
                rows.push( vec![ $n1.to_string(), $n2.to_string() ] );
            }
        };
    }
    */
    macro_rules! bold {
        () => {
            if !plain {
                let mut log = Vec::<u8>::new();
                emit_bold_escape(&mut log);
                print!("{}", strme(&log));
            }
        };
    }
    macro_rules! end_escape {
        () => {
            if !plain {
                let mut log = Vec::<u8>::new();
                emit_end_escape(&mut log);
                print!("{}", strme(&log));
            }
        };
    }
    let mut help_all = false;
    unsafe {
        if HELP_ALL {
            help_all = true;
        }
    }
    macro_rules! begin_doc {
        ($x:expr) => {
            rows.clear();
            if help_all {
                banner($x, plain);
            }
        };
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide special filtering help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "special") || help_all {
        begin_doc!("special");
        print("\n\\bold{special filtering options}\n\n");
        print(
            "This page documents some options, most of which allow noise \
             filters to be turned off, and which normally should not be invoked.  The last \
             two options can be used to simplify the view of a clonotype.\n\n",
        );
        doc!(
            "NCROSS",
            "If you specify that two or more libraries arose from the same sample (i.e."
        );
        doc!(
            "",
            "from the same tube of cells), then the default behavior of enclone is to"
        );
        doc!(
            "",
            "\"cross filter\" so as to remove expanded exact subclonotypes that are"
        );
        doc!(
            "",
            "present in one library but not another, in a fashion that would be highly"
        );
        doc!(
            "",
            "improbable, assuming random draws of cells from the tube.  These are"
        );
        doc!(
            "",
            "believed to arise when a plasma or plasmablast cell breaks up during during"
        );
        doc!(
            "",
            "or after pipetting from the tube, and the resulting fragments seed GEMs,"
        );
        doc!(
            "",
            "yielding expanded 'fake' clonotypes that are residues of real single plasma"
        );
        doc!(
            "",
            "cells.  The NCROSS options turns off this filter, which could be useful so"
        );
        doc!(
            "",
            "long as you interpret the restored clonotypes as representing what are"
        );
        doc!(
            "",
            "probably single cells.  There may also be other situations where the filter"
        );
        doc!(
            "",
            "should be turned off, and in particular the filter can do weird things if"
        );
        doc!(
            "",
            "inputs are somehow mis-specified to enclone.  Note that for purposes of"
        );
        doc!("", "this option, enclone defines a sample by the pair");
        doc!("", "(sample name, donor name).");
        ldoc!(
            "NGRAPH_FILTER",
            "By default, enclone filters to remove exact subclonotypes that by virtue of"
        );
        doc!(
            "",
            "their relationship to other exact subclonotypes, appear to arise from"
        );
        doc!(
            "",
            "background mRNA or a phenotypically similar phenomenon.  The"
        );
        doc!("", "NGRAPH_FILTER option turns off this filtering.");
        ldoc!(
            "NQUAL",
            "By default, enclone filters out exact subclonotypes having a base in V..J"
        );
        doc!(
            "",
            "that looks like it might be wrong.  More specifically, enclone finds bases"
        );
        doc!(
            "",
            "which are not Q60 for a barcode, not Q40 for two barcodes, are not"
        );
        doc!(
            "",
            "supported by other exact subclonotypes, are variant within the clonotype,"
        );
        doc!(
            "",
            "and which disagree with the donor reference.  NQUAL turns this off."
        );
        ldoc!(
            "NWEAK_CHAINS",
            "By default, enclone filters chains from clonotypes that are"
        );
        doc!(
            "",
            "weak and appear to be artifacts, perhaps arising from a stray mRNA molecule"
        );
        doc!(
            "",
            "that floated into a GEM.  The NWEAK_CHAINS option turns off this filter."
        );
        ldoc!(
            "NWEAK_ONESIES",
            "By default, enclone filters out onesie clonotypes having a single exact"
        );
        doc!(
            "",
            "subclonotype, and that are light chain or TRA, and whose number of cells is"
        );
        doc!(
            "",
            "less than 0.1% of the total number of cells.  NWEAK_ONESIES turns this off."
        );
        ldoc!(
            "NFOURSIE_KILL",
            "By default, if enclone finds a foursie exact subclonotype that"
        );
        doc!(
            "",
            "contains a twosie exact subclonotype having at least ten cells, it kills"
        );
        doc!(
            "",
            "the foursie exact subclonotype, no matter how many cells it has.  The"
        );
        doc!(
            "",
            "foursies that are killed are believed to be rare oddball artifacts arising"
        );
        doc!(
            "",
            "from repeated cell doublets or GEMs that contain two cells and multiple gel"
        );
        doc!(
            "",
            "beads.  The argument NFOURSIE_KILL turns off this filtering."
        );
        ldoc!(
            "NWHITEF",
            "By default, enclone filters out rare artifacts arising from \
             contamination"
        );
        doc!(
            "",
            "of oligos on gel beads.  The NWHITEF option turns off this filter."
        );
        ldoc!(
            "NBC_DUP",
            "By default, enclone filters out duplicated barcodes within an exact."
        );
        doc!(
            "",
            "subclonotype.  The NBC_DUP option turns off this filter."
        );
        ldoc!(
            "MIX_DONORS",
            "By default, enclone will prevent cells from different donors from being"
        );
        doc!(
            "",
            "placed in the same clonotype.  The MIX_DONORS option turns off this"
        );
        doc!(
            "",
            "behavior, thus allowing cells from different donors to be placed in the"
        );
        doc!(
            "",
            "same clonotype.  The main use of this option is for specificity testing, in"
        );

        doc!(
            "",
            "which data from different donors are deliberately combined in an attempt"
        );
        doc!(
            "",
            "to find errors.  Use of the bc field for META input specification"
        );
        doc!("", "automatically turns on this option.");
        ldoc!(
            "KEEP_IMPROPER",
            "An exact subclonotype is improper if it does not have one chain"
        );
        doc!(
            "",
            "of each type.  This option causes all improper exact subclonotypes to be"
        );
        doc!(
            "",
            "retained, although they may be removed by other filters."
        );
        ldoc!(
            "MIN_CHAINS_EXACT=n",
            "Delete any exact subclonotype having less than n chains.  You can use this"
        );
        doc!(
            "",
            "to \"purify\" a clonotype so as to display only exact subclonotypes having"
        );
        doc!("", "all their chains.");
        doc!(
            "MIN_CELLS_EXACT=n",
            "Delete any exact subclonotype having less than n cells.  You might want"
        );
        doc!(
            "",
            "to use this if you have a very large and complex expanded clonotype,"
        );
        doc!("", "for which you would like to see a simplified view.");
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Function that provides an explanation used for both enclone help lvars and
    // enclone help cvars.

    fn explain_alt_versions() {
        print!(
            "{}",
            gray_left_bar(&print_to(
                "\\red{◉} These variables have some alternate versions, \
                 as shown in the table below:\n\n"
            ))
        );
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "variable".to_string(),
            "semantics".to_string(),
            "visual".to_string(),
            "visual".to_string(),
            "parseable".to_string(),
            "parseable".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "".to_string(),
            "".to_string(),
            "".to_string(),
            "(one cell)".to_string(),
            "".to_string(),
            "(one cell)".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 6];
        rows.push(row);
        let row = vec![
            "x".to_string(),
            "median over cells".to_string(),
            "yes".to_string(),
            "this cell".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_mean".to_string(),
            "mean over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_μ".to_string(),
            "(same as above)".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_sum".to_string(),
            "sum over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_Σ".to_string(),
            "(same as above)".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_min".to_string(),
            "min over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_max".to_string(),
            "max over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_cell".to_string(),
            "this cell".to_string(),
            "no".to_string(),
            "no".to_string(),
            "no".to_string(),
            "this cell".to_string(),
        ];
        rows.push(row);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l|l|l|l|l".to_vec(), false, false);
        print!("{}", gray_left_bar(&log));
        print!( "{}", gray_left_bar(&print_to(
            "Some explanation is required.  If you use enclone without certain options, you \
             get the \"visual\" column.\n\
             • Add the option \\bold{PER_CELL} \
             (see \"enclone help display\") and then you get visual output with extra lines for \
             each cell within an exact subclonotype, and each of those extra lines is described by \
             the \"visual (one cell)\" column.\n\
             • If you generate parseable output (see \"enclone help parseable\"), then you get \
             the \"parseable\" column for that output, unless you specify \\bold{PCELL}, \
             and then you get the last column.\n\
             • For the forms with μ and Σ, the Greek letters are only used in column headings for \
             visual output (to save space), and optionally, in names of fields on the command \
             line.\n\
             \\green{▶} If you try out these features, you'll see exactly what happens! \
             \\green{◀}\n"
        )));
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide lvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "lvars") || help_all {
        begin_doc!("lvars");
        println!("");
        bold!();
        println!("lead column options\n");
        end_escape!();
        println!(
            "These options define lead variables, which correspond to columns that \
             appear once in each\nclonotype, on the left side, and have one entry for each \
             exact subclonotype row.\n"
        );
        print(
            "Lead variables are specified using \\bold{LVARS=x1,...,xn} \
             where each xi is one of:\n\n",
        );
        doc!("datasets", "dataset identifiers");
        doc!("samples", "sample identifiers");
        doc!("donors", "donor identifiers");
        ldoc!("n", "number of cells");
        doc!(
            "n_<name>",
            "number of cells associated to the given name, which can be a dataset"
        );
        doc!(
            "",
            "or sample or donor or tag short name; may name only one such category"
        );
        ldoc!(
            "near",
            "Hamming distance of V..J DNA sequence to nearest neighbor"
        );
        doc!(
            "far",
            "Hamming distance of V..J DNA sequence to farthest neighbor"
        );
        doc!(
            "",
            "both compare to cells having chains in the same columns of the clonotype,"
        );
        doc!(
            "",
            "with - shown if there is no other exact subclonotype to compare to"
        );
        ldoc!(
            "g<d>",
            "Here d is a nonnegative integer.  Then all the exact subclonotypes are"
        );
        doc!(
            "",
            "grouped according to the Hamming distance of their V..J sequences.  Those"
        );
        doc!(
            "",
            "within distance d are defined to be in the same group, and this is"
        );
        doc!(
            "",
            "extended transitively.  The group identifier 1, 2, ... is shown.  The"
        );
        doc!(
            "",
            "ordering of these identifiers is arbitrary.  This option is best applied"
        );
        doc!(
            "",
            "to cases where all exact subclonotypes have a complete set of chains."
        );
        ldocpr!("gex", "\\red{◉} median gene expression UMI count");
        docpr!("n_gex", "\\blue{◉} number of cells reported by GEX");
        // nonpublic for now as we don't know if this is useful
        /*
        doc!(
            "entropy",
            "Shannon entropy of GEX UMI counts (median across cells)"
        );
        */
        ldocpr!(
            "<gene>_g",
            "\\red{◉} all five feature types: look for a declared feature of the \
             given type"
        );
        doc!(
            "<antibody>_ab",
            "with the given id or name; report the median UMI count for it; we allow"
        );
        doc!(
            "<antigen>_ag",
            "the form e.g. <abbr>:<gene>_g where abbr is an abbreviation to be shown"
        );
        doc!("<crispr>_cr", "");
        doc!("<custom>_cu", "");
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        print!("{}", log);
        print(
            "For gene expression and feature barcode stats, such data must be provided \
             as input to enclone.\n\n",
        );
        explain_alt_versions();
        print(
            "\n\\blue{◉} Similar to the above but simpler: n_gex is just a count of cells, \
             visual (one cell) shows 0 or 1, n_gex_cell is defined for parseable (one cell), \
             and the x_mean etc. forms do not apply.\n\n",
        );
        print(
            "The default is \\bold{datasets,n}, except that datasets is suppressed if \
             there is only one dataset.\n\n",
        );
        print("\\bold{LVARSP=x1,...,xn} is like \\bold{LVARS} but appends to the list.\n\n");
        print(
            "Note: gene expression counts are normalized to 20,000 read pairs per cell, and \
             feature barcode counts are normalized to 5,000 read pairs per cell.  The normalized \
             counts are rounded to the nearest integer.  For this normalization, \
             we simply scale the counts, rather than subsample reads.  If you want to turn off \
             the normalization, add the argument \\bold{FULL_COUNTS} to the command line.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide cvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "cvars") || help_all {
        begin_doc!("cvars");

        // Header.

        print(
            "\n\\bold{per-chain column options}: These options define per-chain variables, \
             which correspond to columns that appear once for each chain in each clonotype, and \
             have one entry for each exact subclonotype.\n\n",
        );
        print(
            "Per-column variables are specified using\n\
             \\bold{CVARS=x1,...,xn}\n\
             where each xi is one of:\n\n",
        );

        // Main table entries.

        doc!(
            "var",
            "bases at positions in chain that vary across the clonotype"
        );
        ldocpr!(
            "u",
            "\\red{◉} VDJ UMI count for each exact subclonotype, median across cells"
        );
        docpr!(
            "r",
            "\\red{◉} VDJ read count for each exact subclonotype, median across cells"
        );
        ldoc!("const", "constant region name");
        ldoc!(
            "edit",
            "a string that partially defines the edit of the reference V(D)J concatenation"
        );
        doc!(
            "",
            "that gives rise to the observed CDR3; this uses a coordinate system in which"
        );
        doc!(
            "",
            "0 is the first base of the J ref segment (or the first base of the D ref"
        );
        doc!(
            "",
            "segment for IGH and TRB); for example D-4:4 denotes the deletion of the last"
        );
        doc!(
            "",
            "4 bases of the V segment, I0:2 denotes an insertion of 2 bases after the V"
        );
        doc!(
            "",
            "and I0:2;S5 denotes that plus a substitution at position 5; in computing"
        );
        doc!(
            "",
            "\"edit\", for IGH and TRB, we always test every possible D segment,"
        );
        doc!(
            "",
            "regardless of whether one is annotated, and pick the best one; for this"
        );
        doc!("", "reason, \"edit\" may be slow");
        doc!(
            "comp",
            "a measure of CDR3 complexity, which is the total number of S, D and I"
        );
        doc!("", "symbols in \"edit\" as defined above");
        ldoc!("cdr3_dna", "the CDR3_DNA sequence");
        ldoc!(
            "vjlen",
            "number of bases from the start of the V region to the end of the J region"
        );
        doc!(
            "clen",
            "length of observed constant region (usually truncated at primer start)"
        );
        doc!("ulen", "length of observed 5'-UTR sequence");
        doc!(
            "cdiff",
            "differences with universal reference constant region, shown in the"
        );
        doc!(
            "",
            "abbreviated form e.g. 22T (ref changed to T at base 22) or 22T+10"
        );
        doc!(
            "",
            "(same but contig has 10 additional bases beyond end of ref C region"
        );
        doc!(
            "",
            "At most five differences are shown, and if there are more, ... is appended."
        );
        doc!("udiff", "like cdiff, but for the 5'-UTR");
        ldoc!(
            "notes",
            "optional note if there is an insertion or the end of J does not exactly abut"
        );
        doc!("", "the beginning of C; elided if empty");
        ldoc!(
            "ndiff<n>",
            "number of base differences within V..J between this exact subclonotype and"
        );
        doc!("", "exact subclonotype n");
        doc!(
            "d_univ",
            "distance from universal reference, more specifically,"
        );
        doc!(
            "",
            "number of base differences within V..J between this exact"
        );
        doc!(
            "",
            "clonotype and universal reference, exclusive of indels, the last 15"
        );
        doc!("", "bases of the V and the first 15 bases of the J");
        doc!("d_donor", "distance from donor reference,");
        doc!("", "as above but computed using donor reference");

        // The rest.

        let mut log = String::new();
        // was ... rows.clone()
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        explain_alt_versions();
        print(
            "\nAt least one variable must be listed.  The default is \\bold{u,const,notes}.  \
             \\bold{CVARSP}: same as \\bold{CVARS} but appends.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide amino help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "amino") || help_all {
        begin_doc!("amino");
        print(
            "\nThere is a complex per-chain column to the left of other \
             per-chain columns, defined by\n\
             \\bold{AMINO=x1,...,xn}: display amino acid columns for the given categories, \
             in one combined ordered group, where each xi is one of:\n\n",
        );
        doc!("cdr3", "CDR3 sequence");
        ldoc!("var", "positions in chain that vary across the clonotype");
        doc!(
            "share",
            "positions in chain that differ consistently from the donor reference"
        );
        ldoc!(
            "donor",
            "positions in chain where the donor reference differs from the universal \
             reference"
        );
        ldoc!(
            "donorn",
            "positions in chain where the donor reference differs nonsynonymously"
        );
        doc!("", "from the universal reference");
        ldoc!(
            "a-b",
            "amino acids numbered a through b (zero-based, inclusive)"
        );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        print(
            "Note that we compute positions in base space, and then divide by three to get \
             positions in amino acid space.  Thus it can happen that a position in amino acid \
             space is shown for both \\bold{var} and \\bold{share}.\n\n",
        );
        print(
            "The default value for \\bold{AMINO} is \\bold{cdr3,var,share,donor}.  \
             Note that we only report amino acids that are strictly within V..J, \
             thus specifically excluding the codon bridging J and C.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide display help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "display") || help_all {
        begin_doc!("display");
        print("\n\\bold{other options that control clonotype display}\n\n");
        doc!(
            "PER_CELL",
            "expand out each exact clonotype line, showing one line per cell,"
        );
        doc!(
            "",
            "for each such line, displaying the barcode name, the number of UMIs assigned,"
        );
        doc!(
            "",
            "and the gene expression UMI count, if applicable, under gex_med"
        );
        ldoc!(
            "BARCODES",
            "print list of all barcodes of the cells in each clonotype, in a"
        );
        doc!(
            "",
            "single line near the top of the printout for a given clonotype"
        );
        ldoc!(
            "SEQC",
            "print V..J sequence for each chain in the first exact subclonotype, near"
        );
        doc!("", "the top of the printout for a given clonotype");
        ldoc!(
            "FULL_SEQC",
            "print full sequence for each chain in the first exact subclonotype,"
        );
        doc!("", "near the top of the printout for a given clonotype");
        ldoc!("SUM", "print sum row for each clonotype");
        doc!("MEAN", "print mean row for each clonotype");
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        print(
            "\\bold{options that control clonotype grouping}\n\n\
             we plan to add grouping capability in a future version of enclone, and for now, as a \
             placeholder, we have the following \"toy\" options:\n\n",
        );
        rows.clear();

        doc!(
            "GROUP_HEAVY_CDR3",
            "group by perfect identity of CDR3 amino acid sequence \
             of IGH or TRB"
        );
        doc!(
            "GROUP_VJ_REFNAME",
            "group by sharing identical V and J reference gene names,"
        );
        doc!("", "but ignores foursies and moresies");
        ldoc!(
            "MIN_GROUP",
            "minimum number of clonotypes in group to print (default = 1)"
        );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        if !help_all {
            std::process::exit(0);
        }
    }
}
