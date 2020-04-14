// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::*;
use tables::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help4(args: &Vec<String>, mut h: &mut HelpDesk) {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide special filtering help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "special") || h.help_all {
        h.begin_doc("special");
        h.print("\n\\bold{special filtering options}\n\n");
        h.print(
            "This page documents some options, most of which allow noise \
             filters to be turned off, and which normally should not be invoked.  The last \
             two options can be used to simplify the view of a clonotype.\n\n",
        );
        h.doc(
            "NGEX",
            "If gene expression and/or feature barcode data are provided, if a barcode",
        );
        h.doc(
            "",
            "is called a cell by the VDJ part of the Cell Ranger pipeline, but not",
        );
        h.doc(
            "",
            "called a cell by the gene expression and/or feature barcode part, then the",
        );
        h.doc(
            "",
            "default behavior of enclone is to remove such cells from clonotypes.  This",
        );
        h.doc("", "option disables that behavior.");
        h.ldoc(
            "NCROSS",
            "If you specify that two or more libraries arose from the same sample (i.e.",
        );
        h.doc(
            "",
            "from the same tube of cells), then the default behavior of enclone is to",
        );
        h.doc(
            "",
            "\"cross filter\" so as to remove expanded exact subclonotypes that are",
        );
        h.doc(
            "",
            "present in one library but not another, in a fashion that would be highly",
        );
        h.doc(
            "",
            "improbable, assuming random draws of cells from the tube.  These are",
        );
        h.doc(
            "",
            "believed to arise when a plasma or plasmablast cell breaks up during during",
        );
        h.doc(
            "",
            "or after pipetting from the tube, and the resulting fragments seed GEMs,",
        );
        h.doc(
            "",
            "yielding expanded 'fake' clonotypes that are residues of real single plasma",
        );
        h.doc(
            "",
            "cells.  The NCROSS options turns off this filter, which could be useful so",
        );
        h.doc(
            "",
            "long as you interpret the restored clonotypes as representing what are",
        );
        h.doc(
            "",
            "probably single cells.  There may also be other situations where the filter",
        );
        h.doc(
            "",
            "should be turned off, and in particular the filter can do weird things if",
        );
        h.doc(
            "",
            "inputs are somehow mis-specified to enclone.  Note that for purposes of",
        );
        h.doc("", "this option, enclone defines a sample by the pair");
        h.doc("", "(sample name, donor name).");
        h.ldoc(
            "NGRAPH_FILTER",
            "By default, enclone filters to remove exact subclonotypes that by virtue of",
        );
        h.doc(
            "",
            "their relationship to other exact subclonotypes, appear to arise from",
        );
        h.doc(
            "",
            "background mRNA or a phenotypically similar phenomenon.  The",
        );
        h.doc("", "NGRAPH_FILTER option turns off this filtering.");
        h.ldoc(
            "NQUAL",
            "By default, enclone filters out exact subclonotypes having a base in V..J",
        );
        h.doc(
            "",
            "that looks like it might be wrong.  More specifically, enclone finds bases",
        );
        h.doc(
            "",
            "which are not Q60 for a barcode, not Q40 for two barcodes, are not",
        );
        h.doc(
            "",
            "supported by other exact subclonotypes, are variant within the clonotype,",
        );
        h.doc(
            "",
            "and which disagree with the donor reference.  NQUAL turns this off.",
        );
        h.ldoc(
            "NWEAK_CHAINS",
            "By default, enclone filters chains from clonotypes that are",
        );
        h.doc(
            "",
            "weak and appear to be artifacts, perhaps arising from a stray mRNA molecule",
        );
        h.doc(
            "",
            "that floated into a GEM.  The NWEAK_CHAINS option turns off this filter.",
        );
        h.ldoc(
            "NWEAK_ONESIES",
            "By default, enclone filters out onesie clonotypes having a single exact",
        );
        h.doc(
            "",
            "subclonotype, and that are light chain or TRA, and whose number of cells is",
        );
        h.doc(
            "",
            "less than 0.1% of the total number of cells.  NWEAK_ONESIES turns this off.",
        );
        h.ldoc(
            "NFOURSIE_KILL",
            "By default, if enclone finds a foursie exact subclonotype that",
        );
        h.doc(
            "",
            "contains a twosie exact subclonotype having at least ten cells, it kills",
        );
        h.doc(
            "",
            "the foursie exact subclonotype, no matter how many cells it has.  The",
        );
        h.doc(
            "",
            "foursies that are killed are believed to be rare oddball artifacts arising",
        );
        h.doc(
            "",
            "from repeated cell doublets or GEMs that contain two cells and multiple gel",
        );
        h.doc(
            "",
            "beads.  The argument NFOURSIE_KILL turns off this filtering.",
        );
        h.ldoc(
            "NWHITEF",
            "By default, enclone filters out rare artifacts arising from \
             contamination",
        );
        h.doc(
            "",
            "of oligos on gel beads.  The NWHITEF option turns off this filter.",
        );
        h.ldoc(
            "NBC_DUP",
            "By default, enclone filters out duplicated barcodes within an exact.",
        );
        h.doc(
            "",
            "subclonotype.  The NBC_DUP option turns off this filter.",
        );
        h.ldoc(
            "MIX_DONORS",
            "By default, enclone will prevent cells from different donors from being",
        );
        h.doc(
            "",
            "placed in the same clonotype.  The MIX_DONORS option turns off this",
        );
        h.doc(
            "",
            "behavior, thus allowing cells from different donors to be placed in the",
        );
        h.doc(
            "",
            "same clonotype.  The main use of this option is for specificity testing, in",
        );

        h.doc(
            "",
            "which data from different donors are deliberately combined in an attempt",
        );
        h.doc(
            "",
            "to find errors.  Use of the bc field for META input specification",
        );
        h.doc("", "automatically turns on this option.");
        h.ldoc(
            "KEEP_IMPROPER",
            "An exact subclonotype is improper if it does not have one chain",
        );
        h.doc(
            "",
            "of each type.  This option causes all improper exact subclonotypes to be",
        );
        h.doc(
            "",
            "retained, although they may be removed by other filters.",
        );
        h.ldoc(
            "MIN_CHAINS_EXACT=n",
            "Delete any exact subclonotype having less than n chains.  You can use this",
        );
        h.doc(
            "",
            "to \"purify\" a clonotype so as to display only exact subclonotypes having",
        );
        h.doc("", "all their chains.");
        h.doc(
            "MIN_CELLS_EXACT=n",
            "Delete any exact subclonotype having less than n cells.  You might want",
        );
        h.doc(
            "",
            "to use this if you have a very large and complex expanded clonotype,",
        );
        h.doc("", "for which you would like to see a simplified view.");
        h.print_tab2();
        h.print("\n");
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Function that provides an explanation used for both enclone help lvars and
    // enclone help cvars.

    fn explain_alt_versions(h: &mut HelpDesk) {
        h.print(&format!(
            "{}",
            gray_left_bar(&print_to(
                "\\red{●} These variables have some alternate versions, \
                 as shown in the table below:\n\n"
            ))
        ));
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
            "x_%".to_string(),
            "% of total GEX (genes only)".to_string(),
            "yes".to_string(),
            "this cell".to_string(),
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
        h.print_plain(&format!("{}", gray_left_bar(&log)));
        h.print_plain(&format!(
            "{}",
            gray_left_bar(&print_to(
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
            ))
        ));
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide lvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "lvars") || h.help_all {
        h.begin_doc("lvars");
        h.print("\n\\bold{lead column options}\n\n");
        h.print(
            "These options define lead variables, which correspond to columns that \
             appear once in each clonotype, on the left side, and have one entry for each \
             exact subclonotype row.\n\n",
        );
        h.print(
            "Lead variables are specified using \\bold{LVARS=x1,...,xn} \
             where each xi is one of:\n\n",
        );
        h.doc("datasets", "dataset identifiers");
        h.doc("samples", "sample identifiers");
        h.doc("donors", "donor identifiers");
        h.ldoc("n", "number of cells");
        h.doc(
            "n_<name>",
            "number of cells associated to the given name, which can be a dataset",
        );
        h.doc(
            "",
            "or sample or donor or tag short name; may name only one such category",
        );
        h.ldoc(
            "near",
            "Hamming distance of V..J DNA sequence to nearest neighbor",
        );
        h.doc(
            "far",
            "Hamming distance of V..J DNA sequence to farthest neighbor",
        );
        h.doc(
            "",
            "both compare to cells having chains in the same columns of the clonotype,",
        );
        h.doc(
            "",
            "with - shown if there is no other exact subclonotype to compare to",
        );
        h.ldoc(
            "g<d>",
            "Here d is a nonnegative integer.  Then all the exact subclonotypes are",
        );
        h.doc(
            "",
            "grouped according to the Hamming distance of their V..J sequences.  Those",
        );
        h.doc(
            "",
            "within distance d are defined to be in the same group, and this is",
        );
        h.doc(
            "",
            "extended transitively.  The group identifier 1, 2, ... is shown.  The",
        );
        h.doc(
            "",
            "ordering of these identifiers is arbitrary.  This option is best applied",
        );
        h.doc(
            "",
            "to cases where all exact subclonotypes have a complete set of chains.",
        );
        h.ldocpr("gex", "\\red{●} median gene expression UMI count");
        h.docpr("n_gex", "\\blue{●} number of cells reported by GEX");
        // nonpublic for now as we don't know if this is useful
        /*
        h.doc(
            "entropy",
            "Shannon entropy of GEX UMI counts (median across cells)"
        );
        */
        h.ldocpr(
            "<gene>_g",
            "\\red{●} all five feature types: look for a declared feature of the \
             given type",
        );
        h.doc(
            "<antibody>_ab",
            "with the given id or name; report the median UMI count for it; we allow",
        );
        h.doc(
            "<antigen>_ag",
            "the form e.g. <abbr>:<gene>_g where abbr is an abbreviation to be shown;",
        );
        h.doc(
            "<crispr>_cr",
            "we also allow <regular expression>_g where g can be replaced by ab, ag, cr",
        );
        h.doc(
            "<custom>_cu",
            "or cu; this represents a sum of UMI counts across the matching features. ●",
        );
        h.ldoc(
            "cred",
            "Short for credibility.  It is a measure of the extent to which cells",
        );
        h.doc(
            "",
            "having gene expression similar to a given putative B cell are themselves",
        );
        h.doc(
            "",
            "B cells.  (Or similarly for T cells.)  For the actual definition, let n",
        );
        h.doc(
            "",
            "be the number of VDJ cells that are also GEX cells.  For a given cell,",
        );
        h.doc(
            "",
            "find the n GEX cells that are closest to it in PCA space, and report the",
        );
        h.doc(
            "",
            "percent of those that are also VDJ cells.  For multiple datasets, it would",
        );
        h.doc(
            "",
            "be better to \"aggr\" the data, however that is not currently supported",
        );
        h.doc(
            "",
            "The computation is also inefficient, so let us know if it's causing",
        );
        h.doc(
            "",
            "problems for you.  And cred makes much better sense for samples that",
        );
        h.doc(
            "",
            "consist of mixed cell types, rather than consisting of pure B or T cells.",
        );
        h.print_tab2();
        h.print(
            "For gene expression and feature barcode stats, such data must be provided \
             as input to enclone.\n\n",
        );
        h.print(
            "● Example: IG.*_g matches all genes that begin with IG.  The regular expression must \
             be in the alphabet A-Za-z0-9+_-.[]()|* and is only interpreted as a regular expression \
             if it contains a character in []()|*.  See \"enclone help filter\" \
             for more information about regular expressions.\n\n",
        );
        explain_alt_versions(&mut h);
        h.print(
            "\n\\blue{●} Similar to the above but simpler: n_gex is just a count of cells, \
             visual (one cell) shows 0 or 1, n_gex_cell is defined for parseable (one cell), \
             and the x_mean etc. forms do not apply.\n\n",
        );
        h.print(
            "The default is \\bold{datasets,n}, except that datasets is suppressed if \
             there is only one dataset.\n\n",
        );
        h.print("\\bold{LVARSP=x1,...,xn} is like \\bold{LVARS} but appends to the list.\n\n");
        h.print(
            "Note: gene expression counts are normalized to 20,000 read pairs per cell, and \
             feature barcode counts are normalized to 5,000 read pairs per cell.  The normalized \
             counts are rounded to the nearest integer.  For this normalization, \
             we simply scale the counts, rather than subsample reads.  If you want to turn off \
             the normalization, add the argument \\bold{FULL_COUNTS} to the command line.\n\n",
        );
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide cvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "cvars") || h.help_all {
        h.begin_doc("cvars");

        // Header.

        h.print(
            "\n\\bold{per-chain column options}: These options define per-chain variables, \
             which correspond to columns that appear once for each chain in each clonotype, and \
             have one entry for each exact subclonotype.\n\n",
        );
        h.print(
            "Per-column variables are specified using\n\
             \\bold{CVARS=x1,...,xn}\n\
             where each xi is one of:\n\n",
        );

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
        h.ldoc("const", "constant region name");
        h.ldoc(
            "edit",
            "a string that partially defines the edit of the reference V(D)J concatenation",
        );
        h.doc(
            "",
            "that gives rise to the observed CDR3; this uses a coordinate system in which",
        );
        h.doc(
            "",
            "0 is the first base of the J ref segment (or the first base of the D ref",
        );
        h.doc(
            "",
            "segment for IGH and TRB); for example D-4:4 denotes the deletion of the last",
        );
        h.doc(
            "",
            "4 bases of the V segment, I0:2 denotes an insertion of 2 bases after the V",
        );
        h.doc(
            "",
            "and I0:2;S5 denotes that plus a substitution at position 5; in computing",
        );
        h.doc(
            "",
            "\"edit\", for IGH and TRB, we always test every possible D segment,",
        );
        h.doc(
            "",
            "regardless of whether one is annotated, and pick the best one; for this",
        );
        h.doc("", "reason, \"edit\" may be slow");
        h.doc(
            "comp",
            "a measure of CDR3 complexity, which is the total number of S, D and I",
        );
        h.doc("", "symbols in \"edit\" as defined above");
        h.ldoc("cdr3_dna", "the CDR3_DNA sequence");
        h.ldoc(
            "vjlen",
            "number of bases from the start of the V region to the end of the J region",
        );
        h.doc(
            "clen",
            "length of observed constant region (usually truncated at primer start)",
        );
        h.doc("ulen", "length of observed 5'-UTR sequence");
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
        h.ldoc(
            "notes",
            "optional note if there is an insertion or the end of J does not exactly abut",
        );
        h.doc("", "the beginning of C; elided if empty");
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

        h.print_tab2();
        h.print("\n");
        explain_alt_versions(&mut h);
        h.print(
            "\nAt least one variable must be listed.  The default is \\bold{u,const,notes}.  \
             \\bold{CVARSP}: same as \\bold{CVARS} but appends.\n\n",
        );
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide amino help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "amino") || h.help_all {
        h.begin_doc("amino");
        h.print(
            "\nThere is a complex per-chain column to the left of other \
             per-chain columns, defined by\n\
             \\bold{AMINO=x1,...,xn}: display amino acid columns for the given categories, \
             in one combined ordered group, where each xi is one of:\n\n",
        );
        h.doc("cdr3", "CDR3 sequence");
        h.ldoc("var", "positions in chain that vary across the clonotype");
        h.doc(
            "share",
            "positions in chain that differ consistently from the donor reference",
        );
        h.ldoc(
            "donor",
            "positions in chain where the donor reference differs from the universal \
             reference",
        );
        h.ldoc(
            "donorn",
            "positions in chain where the donor reference differs nonsynonymously",
        );
        h.doc("", "from the universal reference");
        h.ldoc(
            "a-b",
            "amino acids numbered a through b (zero-based, inclusive)",
        );
        h.print_tab2();
        h.print("\n");
        h.print(
            "Note that we compute positions in base space, and then divide by three to get \
             positions in amino acid space.  Thus it can happen that a position in amino acid \
             space is shown for both \\bold{var} and \\bold{share}.\n\n",
        );
        h.print(
            "The default value for \\bold{AMINO} is \\bold{cdr3,var,share,donor}.  \
             Note that we only report amino acids that are strictly within V..J, \
             thus specifically excluding the codon bridging J and C.\n\n",
        );
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide display help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "display") || h.help_all {
        h.begin_doc("display");
        h.print("\n\\bold{other options that control clonotype display}\n\n");
        h.doc(
            "PER_CELL",
            "expand out each exact clonotype line, showing one line per cell,",
        );
        h.doc(
            "",
            "for each such line, displaying the barcode name, the number of UMIs assigned,",
        );
        h.doc(
            "",
            "and the gene expression UMI count, if applicable, under gex_med",
        );
        h.ldoc(
            "BARCODES",
            "print list of all barcodes of the cells in each clonotype, in a",
        );
        h.doc(
            "",
            "single line near the top of the printout for a given clonotype",
        );
        h.ldoc(
            "SEQC",
            "print V..J sequence for each chain in the first exact subclonotype, near",
        );
        h.doc("", "the top of the printout for a given clonotype");
        h.ldoc(
            "FULL_SEQC",
            "print full sequence for each chain in the first exact subclonotype,",
        );
        h.doc("", "near the top of the printout for a given clonotype");
        h.ldoc("SUM", "print sum row for each clonotype");
        h.doc("MEAN", "print mean row for each clonotype");
        h.print_tab2();
        h.print("\n");
        h.print(
            "\\bold{options that control clonotype grouping}\n\n\
             we plan to add grouping capability in a future version of enclone, and for now, as a \
             placeholder, we have the following \"toy\" options:\n\n",
        );
        h.rows.clear();

        h.doc(
            "GROUP_HEAVY_CDR3",
            "group by perfect identity of CDR3 amino acid sequence \
             of IGH or TRB",
        );
        h.doc(
            "GROUP_VJ_REFNAME",
            "group by sharing identical V and J reference gene names,",
        );
        h.doc("", "but ignores foursies and moresies");
        h.ldoc(
            "MIN_GROUP",
            "minimum number of clonotypes in group to print (default = 1)",
        );
        h.print_tab2();
        h.print("\n");
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }
}
