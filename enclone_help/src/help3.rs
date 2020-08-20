// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::*;
use tables::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help3(args: &Vec<String>, h: &mut HelpDesk) {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide input_tech help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "input_tech") || h.help_all {
        h.begin_doc("input_tech");
        h.print("\n\\bold{information about providing input to enclone (technical notes)}\n\n");
        h.print(
            "enclone only uses certain files, which are all in the outs subdirectory of \
             a Cell Ranger pipeline directory:\n\n",
        );
        h.docpr("\\bold{file}", "\\bold{pipeline}");
        h.ldoc("all_contig_annotations.json", "VDJ");
        h.ldoc("vdj_reference/fasta/regions.fa", "VDJ");
        h.ldoc("metrics_summary.csv", "GEX");
        h.ldoc("raw_feature_bc_matrix.h5", "GEX");
        h.ldoc("analysis/clustering/graphclust/clusters.csv", "GEX");
        h.ldoc("analysis/pca/10_components/projection.csv", "GEX");
        h.print_tab2();
        h.print(
            "\nThe first file is required, and the second should be supplied if Cell Ranger \
            version 4.0 or greater was used.  The others are required, in the indicated \
            structure, if GEX or META/gex arguments are provided.  The exact files \
            that are used could be changed in the future.\n\n",
        );
        h.print(
            "Note that the VDJ outs directories must be from Cell Ranger version \
             \\boldred{≥ 3.1}.  There \
             is a workaround for earlier versions (which you will be informed of if you try), but \
             it is much slower and the results may not be as good.\n\n",
        );
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide parseable output help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "parseable") || h.help_all {
        h.begin_doc("parseable");
        h.print("\n");
        h.print("\\bold{parseable output}\n");
        h.print(
            "\nThe standard output of enclone is designed to be read by humans, but is not \
             readily parseable by computers.  We supplement this with parseable output that can \
             be easily read by computers.\n\n\
             \
             The default behavior for this is to generate a CSV file having \"every possible\" \
             field (over a hundred).  We also provide an option to print only selected fields, \
             and some options which enable inspection, short of generating a separate CSV file.\n\n\
             \
             Parseable output is targeted primarily at R and Python users, because of the ease of \
             wrangling CSV files with these languages.\n\n",
        );
        h.print_with_box(
            "Parseable output is invoked by using the argument\n\
             \\bold{POUT=filename}\n\
             specifying the name of the file that is to be written to.\n\
             \
             [47m [0m The filename \"stdout\" may be used for a preview; in that case \
             parseable output is generated\n\
             [47m [0m separately for each clonotype and the two output types \
             are integrated.  There is also\n\
             [47m [0m \"stdouth\", which is similar, but uses spaces instead \
             of commas, and lines things up in columns.\n\
             \
             By default, we show four chains for each clonotype, regardless of how many chains it\n\
             has, filling in with null entries.  One may instead specify n chains using the \
             argument\n\
             \\bold{PCHAINS=n}\n\
             The parseable output fields may be specified using\n\
             \\bold{PCOLS=x1,...,xn}\n\
             where each xi is one of the field names shown below.\n\\boldred{This option reduces} \
             \\boldred{run time and memory usage, and prevents voluminous output.  Please use it!}",
            true,
        );
        h.print(
            "Over time additional fields may be added and the order of fields may \
             change.\n\n",
        );
        h.print(
            "There is an alternate parseable output mode in which one line is emitted for each \
             cell, rather then each exact subclonotype.  This mode is enabled by adding the \
             argument \\bold{PCELL} to the command line.  Each exact subclonotype then yields a \
             sequence of output lines that are identical except as noted below.\n\n",
        );
        h.print(
            "If you want to completely suppress the generation of visual clonotypes, add \
             \\bold{NOPRINT} to the enclone command line.\n\n",
        );
        h.print_with_box(
            "\\bold{FASTA output.}  This is a separate feature.  \
             To generate nucleotide FASTA output for each chain in each exact subclonotype, \
             use the argument \\bold{FASTA=filename}.  The special case \\bold{stdout} will \
             cause the FASTA records to be shown as part of standard output.  The FASTA records \
             that are generated are of the form V(D)JC, where V is the full V segment (including \
             the leader) and C is the full constant region, copied verbatim from the reference.  \
             If a particular chain in a particular exact subclonotype is not assigned a constant \
             region, then we use the constant region that was assigned to the clonotype.  If no \
             constant region at all was assigned, then the FASTA record is omitted.  \
             Similarly, \\bold{FASTA_AA=filename} may be used to generate a matching amino acid \
             FASTA file.",
            true,
        );
        h.print(
            "\\boldred{───────────────────────}\n\
             \\boldred{parseable output fields}\n\
             \\boldred{───────────────────────}\n\n",
        );
        h.rows.clear();
        h.print("\\bold{1. per clonotype group fields}\n\n");
        h.doc("group_id", "identifier of clonotype group - 0,1, ...");
        h.ldoc("group_ncells", "total number of cells in the group");
        h.print_tab2();
        h.print("\n");

        h.rows.clear();
        h.print("\\bold{2. per clonotype fields}\n\n");
        h.doc(
            "clonotype_id",
            "identifier of clonotype within the clonotype group = 0, 1, ...",
        );
        h.ldoc("clonotype_ncells", "total number of cells in the clonotype");
        h.ldoc("nchains", "total number of chains in the clonotype");
        h.print_tab2();
        h.print("\n");

        h.rows.clear();
        h.print(
            "\\bold{3. per chain fields, where <i> is 1,2,... (see above)\n\
             each of these has the same value for each exact clonotype}\n\n",
        );
        h.doc("v_name<i>", "name of V segment");
        h.doc("d_name<i>", "name of D segment (or null)");
        h.doc("j_name<i>", "name of J segment");
        h.ldoc("v_id<i>", "id of V segment");
        h.doc("d_id<i>", "id of D segment (or null)");
        h.doc("j_id<i>", "id of J segment");
        h.ldoc(
            "var_indices_dna<i>",
            "DNA positions in chain that vary across the clonotype",
        );
        h.doc(
            "var_indices_aa<i>",
            "amino acid positions in chain that vary across the clonotype",
        );
        h.doc(
            "share_indices_dna<i>",
            "DNA positions in chain that are constant across the \
             clonotype,",
        );
        h.doc("", "but differ from the donor ref");
        h.doc(
            "share_indices_aa<i>",
            "amino acid positions in chain that are constant across the \
             clonotype,",
        );
        h.doc("", "all of these are comma-separated lists");
        h.doc("", "but differ from the donor ref");
        h.print_tab2();
        h.print("\n");

        h.rows.clear();
        h.print("\\bold{4. per exact subclonotype fields}\n\n");
        h.doc(
            "exact_subclonotype_id",
            "identifer of exact subclonotype = 1, 2, ...",
        );
        h.ldoc(
            "barcodes",
            "comma-separated list of barcodes for the exact subclonotype",
        );
        h.doc(
            "<dataset>_barcodes",
            "like \"barcodes\", but restricted to the dataset with the given name",
        );
        h.doc("barcode", "if PCELL is specified, barcode for one cell");
        h.doc(
            "<dataset>_barcode",
            "if PCELL is specified, barcode for one cell, or null, if the barcode is",
        );
        h.doc("", "not from the given dataset");
        h.ldoc(
            "In addition, every lead variable may be specified as a field.  \
             See \"enclone help lvars\".",
            "\\ext",
        );
        h.print_tab2();
        h.print("\n");

        h.rows.clear();
        h.print(
            "\\bold{5. per chain, per exact subclonotype fields, where <i> is 1,2,... \
             (see above)}\n\n",
        );
        h.print("[all apply to chain i of a particular exact clonotype]\n\n");
        h.doc("vj_seq<i>", "DNA sequence of V..J");
        h.doc("seq<i>", "full DNA sequence");
        h.doc(
            "q<n>_<i>",
            "special option to display a comma-separated list of the quality",
        );
        h.doc(
            "",
            "scores for chain i, at zero-based position n, numbered starting at the",
        );
        h.doc(
            "",
            "beginning of the V segment, for each cell in the exact subclonotype",
        );
        h.ldoc("v_start<i>", "start of V segment on full DNA sequence");
        h.ldoc(
            "const_id<i>",
            "numerical identifier of constant region (or null, if not known)",
        );
        h.ldoc(
            "utr_id<i>",
            "numerical identifier of 5'-UTR region (or null, if not known)",
        );
        h.doc(
            "utr_name<i>",
            "name of 5'-UTR region (or null, if not known)",
        );
        h.ldoc(
            "cdr3_start<i>",
            "base position start of CDR3 sequence on full contig",
        );
        h.doc("cdr3_aa<i>", "amino acid sequence of CDR3");
        h.ldoc(
            "var_aa<i>",
            "amino acids that vary across the clonotype (synonymous changes included)",
        );
        h.ldoc(
            "In addition, every chain variable, after suffixing by <i>, may be used as a field.",
            "\\ext",
        );
        h.doc("See \"enclone help cvars\".", "\\ext");
        h.print_tab2();
        h.print("\n");
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide filter help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "filter") || h.help_all {
        h.begin_doc("filter");

        // intro

        h.print("\n\\bold{clonotype filtering options}\n\n");
        h.print(
            "enclone provides filtering by cell, by exact subclonotype, and by clonotype.  This \
            page describes filtering by clonotype.  \
            These options cause only certain clonotypes to be printed.  See also \
            \"enclone help special\", which describes other filtering options.  This page \
            also described scanning for feature enrichment.\n\n",
        );

        // doc *CELLS

        h.doc(
            "MIN_CELLS=n",
            "only show clonotypes having at least n cells",
        );
        h.doc("MAX_CELLS=n", "only show clonotypes having at most n cells");
        h.doc("CELLS=n", "only show clonotypes having exactly n cells");

        // doc MIN_UMIS

        h.ldoc(
            "MIN_UMIS=n",
            "only show clonotypes having ≳ n UMIs on some chain on some cell",
        );

        // doc *CHAINS

        h.ldoc(
            "MIN_CHAINS=n",
            "only show clonotypes having at least n chains",
        );
        h.doc(
            "MAX_CHAINS=n",
            "only show clonotypes having at most n chains",
        );
        h.doc("CHAINS=n", "only show clonotypes having exactly n chains");

        // doc CDR3

        h.ldoc(
            "CDR3=<pattern>",
            "only show clonotypes having a CDR3 amino acid seq that matches",
        );
        h.doc("", "the given pattern*, from beginning to end");

        // doc SEG and SEGN

        h.ldoc(
            "SEG=\"s_1|...|s_n\"",
            "only show clonotypes using one of the given reference segment names",
        );
        h.doc(
            "SEGN=\"s_1|...|s_n\"",
            "only show clonotypes using one of the given reference segment numbers",
        );
        h.doc(
            "",
            "both: looks for V, D, J and C segments; double quote only",
        );
        h.doc("", "needed if n > 1");
        h.doc(
            "",
            "For both SEG and SEGN, multiple instances are allowed, and their",
        );
        h.doc("", "effects are cumulative.");

        // doc MIN_EXACTS

        h.ldoc(
            "MIN_EXACTS=n",
            "only show clonotypes having at least n exact subclonotypes",
        );

        // doc VJ

        h.ldoc(
            "VJ=seq",
            "only show clonotypes using exactly the given V..J sequence",
        );
        h.doc("", "(string in alphabet ACGT)");

        // doc MIN_DATASETS and MAX_DATASETS and MIN_DATASET_RATIO

        h.ldoc(
            "MIN_DATASETS=n",
            "only show clonotypes containing cells from at least n datasets",
        );
        h.doc(
            "MAX_DATASETS=n",
            "only show clonotypes containing cells from at most n datasets",
        );
        h.doc(
            "MIN_DATASET_RATIO=n",
            "only show clonotypes having at least n cells and for which the ratio",
        );
        h.doc(
            "",
            "of the number of cells in the must abundant dataset to the next most",
        );
        h.doc("", "abundant one is at least n");

        // doc CDIFF

        h.ldoc(
            "CDIFF",
            "only show clonotypes having a difference in constant region with the",
        );
        h.doc("", "universal reference");

        // doc DEL

        h.ldoc("DEL", "only show clonotypes exhibiting a deletion");

        // doc BARCODE

        h.ldoc(
            "BARCODE=bc1,...,bcn",
            "only show clonotypes that use one of the given barcodes; note that such",
        );
        h.doc(
            "",
            "clonotypes will typically contain cells that are \\bold{not} in your",
        );
        h.doc(
            "",
            "list; if you want to fully restrict to a list of barcodes you can use",
        );
        h.doc(
            "",
            "the \\bold{FCELL} option, please see \"enclone help special\"",
        );
        // doc INKT and MAIT

        h.ldoc(
            "INKT",
            "only show clonotypes for which some exact subclonotype is annotated as",
        );
        h.docpr(
            "",
            "having some iNKT evidence, see \\green{bit.ly/enclone} for details",
        );
        h.ldoc(
            "MAIT",
            "only show clonotypes for which some exact subclonotype is annotated as",
        );
        h.docpr(
            "",
            "having some MAIT evidence, see \\green{bit.ly/enclone} for details",
        );

        // print main table

        h.print_tab2();
        h.print("\n");

        // footnote for CDR3

        h.print("* Examples of how to specify CDR3:\n\n");
        h.print(
            "Two pattern types are allowed: either regular expressions, or \"Levenshtein \
            distance patterns\", as exhibited by examples below.\n\n",
        );
        let mut rows = Vec::<Vec<String>>::new();
        rows.push(vec![
            "CDR3=CARPKSDYIIDAFDIW".to_string(),
            "have exactly this sequence as a CDR3".to_string(),
        ]);
        rows.push(vec![
            "CDR3=\"CARPKSDYIIDAFDIW|CQVWDSSSDHPYVF\"".to_string(),
            "have at least one of these sequences as a CDR3".to_string(),
        ]);
        rows.push(vec![
            "CDR3=\".*DYIID.*\"".to_string(),
            "have a CDR3 that contains DYIID inside it".to_string(),
        ]);
        rows.push(vec![
            "CDR3=\"CQTWGTGIRVF~3\"".to_string(),
            "CDR3s within Levenshtein distance 3 of CQTWGTGIRVF".to_string(),
        ]);
        rows.push(vec![
            "CDR3=\"CQTWGTGIRVF~3|CQVWDSSSDHPYVF~2\"".to_string(),
            "CDR3s within Levenshtein distance 3 of CQTWGTGIRVF".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "or Levenshtein distance 2 of CQVWDSSSDHPYVF".to_string(),
        ]);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        h.print(&format!("{}\n", log));
        h.print(
            "Note that double quotes should be used if the pattern \
             contains characters other than letters.\n\n",
        );
        h.print(
            "A gentle introduction to regular expressions may be found at\n\
             https://en.wikipedia.org/wiki/Regular_expression#Basic_concepts, and a precise\n\
             specification for the regular expression version used by enclone may be found at\n\
             https://docs.rs/regex.\n\n",
        );

        // linear conditions

        h.print(
            "\\bold{linear conditions}\n\n\
             enclone understands linear conditions of the form\n\
             \\bold{c1*v1 ± ... ± cn*vn > d}\n\
             where each ci is a constant, \"ci*\" may be omitted, each vi is a variable, \
             and d is a constant.  Blank spaces are ignored.  The > sign may be replaced by \
             >= or ≥ or < or <= or ≤.  \
             Each vi is a lead variable (see \"\\bold{enclone help lvars}\") that \
             represents a \
             origin/donor/tag count or gene/feature barcode UMI count.  In evaluating the \
             condition, each vi is \
             replaced by the \\bold{mean} of its values across all cells in the clonotype.  \
             Because the minus sign - doubles as a hyphen and is used in some feature names, we \
             allow parentheses around variable names to prevent erroneous parsing, like this \
             \\bold{(IGHV3-7_g) >= 1}.\n\n",
        );

        // bounds

        h.print(
            "\\bold{filtering by linear conditions}\n\n\
             enclone has the capability to filter by bounding certain lead variables, using \
             the command-line argument:\n\
             \\bold{F=\"L\"}\n\
             where L is a linear condition (as defined above).  Currently this is limited to \
             the case where the lead variables have been selected using \\bold{LVARS} or \
             \\bold{LVARSP}!  Multiple bounds may be imposed by using\n\
             multiple instances of \\bold{F=...} .\n\n",
        );

        // feature scanning

        h.print(
            "\\bold{feature scanning}\n\n\
            If gene expression and/or feature barcode data have been generated, \
            enclone can scan all features to find those that are enriched \
            in certain clonotypes relative to certain other clonotypes.  This feature is turned \
            on using the command line argument\n\
            \\bold{SCAN=\"test,control,threshold\"}\n\
            where each of \\bold{test}, \\bold{control} and \\bold{threshold} are linear \
            conditions as defined above.  Blank spaces are ignored.  The \\bold{test} condition \
            defines the \"test clonotypes\" and the \\bold{control} condition defines the \
            \"control clonotypes\".  Currently, the lead variables in \\bold{test} and \
            \\bold{control} must be specified by\n\
            \\bold{LVARS} or \\bold{LVARSP}!  \
            The \\bold{threshold} condition is special: it may use \
            only the variables \"t\" and \"c\" that represent the raw UMI count for \
            a particular gene or feature, for the test (t) or control (c) clonotypes.  \
            To get a meaningful result, you should specify \\bold{MIN_CELLS} appropriately \
            and manually examine the test and control clonotypes to make sure that they make \
            sense.\n\n\
            \
            \\bold{an example}\n\nSuppose that your data are comprised of two origins with datasets
            named pre and post, representing time points relative to some event.  Then\n\
            \\bold{SCAN=\"n_post - 10*n_pre >= 0, n_pre - 0.5*n_post >= 0, t - 2*c >= 0.1\"}\n\
            would define the test clonotypes to be those satisfying \
            n_post >= 10*n_pre (so having far more post cells then pre cells), \
            the control clonotypes to be those satisfying n_pre >= 0.5*n_post (so having lots of \
            pre cells), and thresholding on t >= 2*c * 0.1, so that the feature must \
            have a bit more than twice as many UMIs in the test than the control.  The 0.1 \
            is there to exclude noise from features having very low UMI counts.\n\n\
            \
            Feature scanning is not a proper statistical test.  It is a tool for generating a list \
            of feature candidates that may then be examined in more detail by rerunning \
            enclone using some of the detected features as lead variables (appropriately \
            suffixed).  Ultimately the power of the scan is determined by having \"enough\" \
            cells in both the test and control sets, and in having those sets cleanly defined.\n\n\
            Currently feature scanning requires that each dataset have identical features.\n\n",
        );

        // done

        h.end_doc();
    }
}
