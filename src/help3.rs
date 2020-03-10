// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::*;
use ansi_escape::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;
use tables::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help3(args: &Vec<String>) {
    // Set up.

    let mut args = args.clone();
    let mut rows = Vec::<Vec<String>>::new();
    macro_rules! doc {
        ($n1:expr, $n2:expr) => {
            rows.push(vec![$n1.to_string(), $n2.to_string()]);
        };
    }
    macro_rules! ldoc {
        ($n1:expr, $n2:expr) => {
            rows.push(vec!["\\hline".to_string(); 2]);
            rows.push(vec![$n1.to_string(), $n2.to_string()]);
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

    // Provide input_tech help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "input_tech") || help_all {
        begin_doc!("input_tech");
        println!("");
        bold!();
        println!("information about providing input to enclone (technical notes)\n");
        end_escape!();
        print(
            "enclone only uses certain files, which are all in the outs subdirectory of \
             a Cell Ranger pipeline directory:\n\n",
        );
        let mut log1 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log1);
        }
        log1.append(&mut b"file".to_vec());
        if !plain {
            emit_end_escape(&mut log1);
        }
        let s1 = stringme(&log1);
        let mut log2 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log2);
        }
        log2.append(&mut b"pipeline".to_vec());
        if !plain {
            emit_end_escape(&mut log2);
        }
        let s2 = stringme(&log2);
        doc!(&s1, &s2);
        ldoc!("all_contig_annotations.json", "VDJ");
        ldoc!("metrics_summary_json.json", "GEX");
        ldoc!("raw_gene_bc_matrices_h5.h5", "GEX");
        ldoc!("raw_feature_bc_matrix/barcodes.tsv.gz", "GEX");
        doc!("raw_feature_bc_matrix/features.tsv.gz", "GEX");
        doc!("filtered_feature_bc_matrix/barcodes.tsv.gz", "GEX");
        print_tab2(&rows);
        println!("\nThe exact files that are used could be changed in the future.\n");
        // end_escape!();
        print(
            "Note that you must use the output of Cell Ranger version \\boldred{≥ 3.1}.  There \
             is a workaround for earlier versions (which you will be informed of if you try), but \
             it is much slower and the results may not be as good.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide parseable output help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "parseable") || help_all {
        begin_doc!("parseable");
        println!("");
        bold!();
        println!("parseable output");
        end_escape!();
        println!(
            "\nThe standard output of enclone is designed to be read by humans, but is not\n\
             readily parseable by computers.  We supplement this with parseable output that can\n\
             be read by computers, but is highly verbose and not intended to be read by humans.\n\n\
             This output is in the form of a CSV file having an exhaustive list of around\n\
             140 fields by default.\n\n\
             It is also possible to specify any subset of these fields, and there are a few other\n\
             choices, which we describe.\n\n\
             Parseable output is targeted primarily at R and Python users, because of the ease of\n\
             wrangling CSV files with these languages.\n"
        );
        print_with_box(
            "Parseable output is invoked by using the argument\n\
             \\bold{POUT=filename}\n\
             specifying the name of the file that is to be written to.\n\
             [The filename \"stdout\" may be used for a preview; in that case the parseable \
             output is\n\
             generated separately for each clonotype and the two output types are integrated.]\n\
             By default, we show four chains for each clonotype, regardless of how many chains it\n\
             has, filling in with null entries.  One may instead specify n chains using the \
             argument\n\
             \\bold{PCHAINS=n}\n\
             The parseable output fields may be specified using\n\
             \\bold{PCOLS=x1,...,xn}\n\
             where each xi is one of the field names shown below.",
            true,
        );
        print(
            "Over time additional fields may be added and the order of fields may \
             change.\n\n",
        );
        print(
            "There is an alternate parseable output mode in which one line is emitted for each \
             cell, rather then each exact subclonotype.  This mode is enabled by adding the \
             argument \\bold{PCELL} to the command line.  Each exact subclonotype then yields a \
             sequence of output lines that are identical except as noted below.\n\n",
        );
        print(
            "If you want to completely suppress the generation of visual clonotypes, add \
             \\bold{NOPRINT} to the enclone command line.\n\n",
        );
        print_with_box(
            "\\bold{FASTA output.}  This is a separate feature. \
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
        let mut log = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log);
            emit_red_escape(&mut log);
        }
        fwriteln!(log, "───────────────────────");
        if !plain {
            emit_bold_escape(&mut log);
            emit_red_escape(&mut log);
        }
        log.append(&mut b"parseable output fields\n".to_vec());
        if !plain {
            emit_bold_escape(&mut log);
            emit_red_escape(&mut log);
        }
        fwrite!(log, "───────────────────────");
        if !plain {
            emit_end_escape(&mut log);
        }
        println!("{}\n", strme(&log));
        rows.clear();
        bold!();
        println!("1. per clonotype group fields\n");
        end_escape!();
        doc!("group_id", "identifier of clonotype group - 0,1, ...");
        ldoc!("group_ncells", "total number of cells in the group");
        print_tab2(&rows);
        println!("");

        rows.clear();
        bold!();
        println!("2. per clonotype fields\n");
        end_escape!();
        doc!(
            "clonotype_id",
            "identifier of clonotype within the clonotype group = 0, 1, ..."
        );
        ldoc!("clonotype_ncells", "total number of cells in the clonotype");
        ldoc!("nchains", "total number of chains in the clonotype");
        print_tab2(&rows);
        println!("");

        rows.clear();
        print(
            "\\bold{3. per chain fields, where <i> is 1,2,... (see above)\n\
             each of these has the same value for each exact clonotype}\n\n",
        );
        doc!("v_name<i>", "name of V segment");
        doc!("d_name<i>", "name of D segment (or null)");
        doc!("j_name<i>", "name of J segment");
        ldoc!("v_id<i>", "id of V segment");
        doc!("d_id<i>", "id of D segment (or null)");
        doc!("j_id<i>", "id of J segment");
        ldoc!(
            "var_indices_dna<i>",
            "DNA positions in chain that vary across the clonotype"
        );
        doc!(
            "var_indices_aa<i>",
            "amino acid positions in chain that vary across the clonotype"
        );
        doc!(
            "share_indices_dna<i>",
            "DNA positions in chain that are constant across the \
             clonotype,"
        );
        doc!("", "but differ from the donor ref");
        doc!(
            "share_indices_aa<i>",
            "amino acid positions in chain that are constant across the \
             clonotype,"
        );
        doc!("", "all of these are comma-separated lists");
        doc!("", "but differ from the donor ref");
        print_tab2(&rows);
        println!("");

        rows.clear();
        print("\\bold{4. per exact subclonotype fields}\n\n");
        doc!(
            "exact_subclonotype_id",
            "identifer of exact subclonotype = 1, 2, ..."
        );
        ldoc!(
            "barcodes",
            "comma-separated list of barcodes for the exact subclonotype"
        );
        doc!(
            "<dataset>_barcodes",
            "like \"barcodes\", but restricted to the dataset with the given name"
        );
        doc!("barcode", "if PCELL is specified, barcode for one cell");
        doc!(
            "<dataset>_barcode",
            "if PCELL is specified, barcode for one cell, or null, if the barcode is"
        );
        doc!("", "not from the given dataset");
        ldoc!(
            "In addition, every lead variable may be specified as a field.  \
             See \"enclone help lvars\".",
            "\\ext"
        );
        doc!(
            "However, to use such a field here, the lead variable must be specified using \
             the LVARS",
            "\\ext"
        );
        doc!("or LVARSP options (or be in the default set).", "\\ext");
        print_tab2(&rows);
        println!("");

        rows.clear();
        print(
            "\\bold{5. per chain, per exact subclonotype fields, where <i> is 1,2,... \
             (see above)}\n\n",
        );
        println!("[all apply to chain i of a particular exact clonotype]\n");
        doc!("vj_seq<i>", "DNA sequence of V..J");
        doc!("seq<i>", "full DNA sequence");
        doc!(
            "q<n>_<i>",
            "special option to display a comma-separated list of the quality"
        );
        doc!(
            "",
            "scores for chain i, at zero-based position n, numbered starting at the"
        );
        doc!(
            "",
            "beginning of the V segment, for each cell in the exact subclonotype"
        );
        ldoc!("v_start<i>", "start of V segment on full DNA sequence");
        ldoc!(
            "const_id<i>",
            "numerical identifier of constant region (or null, if not known)"
        );
        ldoc!(
            "utr_id<i>",
            "numerical identifier of 5'-UTR region (or null, if not known)"
        );
        doc!(
            "utr_name<i>",
            "name of 5'-UTR region (or null, if not known)"
        );
        ldoc!(
            "cdr3_start<i>",
            "base position start of CDR3 sequence on full contig"
        );
        doc!("cdr3_aa<i>", "amino acid sequence of CDR3");
        ldoc!(
            "var_aa<i>",
            "amino acids that vary across the clonotype (synonymous changes included)"
        );
        ldoc!(
            "In addition, every chain variable, after suffixing by <i>, may be used as a field.",
            "\\ext"
        );
        doc!(
            "See \"enclone help cvars\".  However to use such a field here, the chain variable \
             must be ",
            "\\ext"
        );
        doc!(
            "specified using the CVARS or CVARSP options (or be in the default set).",
            "\\ext"
        );
        print_tab2(&rows);
        println!("");
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide filter help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "filter") || help_all {
        begin_doc!("filter");
        println!("");

        // intro

        bold!();
        println!("clonotype filtering options\n");
        end_escape!();
        println!("these options cause only certain clonotypes to be printed\n");

        // doc *CELLS

        doc!(
            "MIN_CELLS=n",
            "only show clonotypes having at least n cells"
        );
        doc!("MAX_CELLS=n", "only show clonotypes having at most n cells");
        doc!("CELLS=n", "only show clonotypes having exactly n cells");

        // doc MIN_UMIS

        ldoc!(
            "MIN_UMIS=n",
            "only show clonotypes having ≳ n UMIs on some chain on some cell"
        );

        // doc *CHAINS

        ldoc!(
            "MIN_CHAINS=n",
            "only show clonotypes having at least n chains"
        );
        doc!(
            "MAX_CHAINS=n",
            "only show clonotypes having at most n chains"
        );
        doc!("CHAINS=n", "only show clonotypes having exactly n chains");

        // doc CDR3

        ldoc!(
            "CDR3=<pattern>",
            "only show clonotypes having a CDR3 amino acid seq that matches"
        );
        doc!(
            "",
            "the given pattern (regular expression)*, from beginning to end"
        );

        // doc SEG and SEGN

        ldoc!(
            "SEG=\"s_1|...|s_n\"",
            "only show clonotypes using one of the given VDJ segment names"
        );
        doc!("", "(double quotes only needed if n > 1)");
        doc!(
            "SEGN=\"s_1|...|s_n\"",
            "only show clonotypes using one of the given VDJ segment numbers"
        );
        doc!("", "(double quotes only needed if n > 1)");

        // doc MIN_EXACTS

        ldoc!(
            "MIN_EXACTS=n",
            "only show clonotypes having at least n exact subclonotypes"
        );

        // doc VJ

        ldoc!(
            "VJ=seq",
            "only show clonotypes using exactly the given V..J sequence"
        );
        doc!("", "(string in alphabet ACGT)");

        // doc MIN_DATASETS

        ldoc!(
            "MIN_DATASETS=n",
            "only show clonotypes containing cells from at least n datasets"
        );

        // doc CDIFF

        ldoc!(
            "CDIFF",
            "only show clonotypes having a difference in constant region with the"
        );
        doc!("", "universal reference");

        // doc DEL

        ldoc!("DEL", "only show clonotypes exhibiting a deletion");

        // doc BARCODE

        ldoc!(
            "BARCODE=bc1,...,bcn",
            "only show clonotypes that use one of the given barcodes"
        );

        // print main table

        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);

        // footnote for CDR3

        println!("* Examples of how to specify CDR3:\n");
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
        log.clear();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false, false);
        println!("{}", log);
        println!(
            "Note that double quotes should be used if the pattern \
             contains characters other than letters.\n"
        );
        println!(
            "A gentle introduction to regular expressions may be found at\n\
             https://en.wikipedia.org/wiki/Regular_expression#Basic_concepts, and a precise\n\
             specification for the regular expression version used by enclone may be found at\n\
             https://docs.rs/regex.\n"
        );

        // linear conditions

        print(
            "\\bold{linear conditions}\n\n\
             enclone understands linear conditions of the form\n\
             \\bold{c1*v1 ± ... ± cn*vn > d}\n\
             where each ci is a constant, \"ci*\" may be omitted, each vi is a variable, \
             and d is a constant.  Blank spaces are ignored.  The > sign may be replaced by \
             >= or ≥ or < or <= or ≤.  \
             Each vi is a lead variable (see \"\\bold{enclone help lvars}\") that \
             represents a \
             sample/donor/tag count or gene/feature barcode UMI count.  In evaluating the \
             condition, each vi is \
             replaced by the \\bold{mean} of its values across all cells in the clonotype.  \
             Because the minus sign - doubles as a hyphen and is used in some feature names, we \
             allow parentheses around variable names to prevent erroneous parsing, like this \
             \\bold{(IGHV3-7_g) >= 1}.\n\n",
        );

        // bounds

        print(
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

        print(
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
            only the variables \"t\" and \"c\" that represent the normalized UMI count for \
            a particular gene or feature, for the test (t) or control (c) clonotypes.  \
            To get a meaningful result, you should specify \\bold{MIN_CELLS} appropriately \
            and manually examine the test and control clonotypes to make sure that they make \
            sense.\n\n\
            \
            \\bold{an example}\n\nSuppose that your data are comprised of two samples named pre \
            and post, representing time points relative to some event.  Then\n\
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
            Currently feature scanning requires that each dataset have identical features.\n\n"
        );

        // done

        if !help_all {
            std::process::exit(0);
        }
    }
}
