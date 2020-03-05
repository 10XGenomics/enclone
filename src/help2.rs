// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Test for help request.

use crate::help_utils::*;
use ansi_escape::*;
use std::env;
use string_utils::*;
use tables::*;
use vector_utils::*;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help2(args: &Vec<String>) {
    // Set up.

    let mut args = args.clone();
    let mut rows = Vec::<Vec<String>>::new();
    /*
    macro_rules! doc {
        ($n1:expr, $n2:expr) => {
            rows.push( vec![ $n1.to_string(), $n2.to_string() ] );
        };
    }
    macro_rules! ldoc {
        ($n1:expr, $n2:expr) => {
            rows.push( vec![ "\\hline".to_string(); 2 ] );
            rows.push( vec![ $n1.to_string(), $n2.to_string() ] );
        };
    }
    */
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

    // Provide example1 help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "example1") || help_all {
        begin_doc!("example1");
        println!(
            "\nSuppose you have placed the datasets that enclone comes with in the\n\
             directory /users/jdoe/enclone_data.  Then you can run this command:"
        );
        println!("\n% enclone PRE=/users/jdoe/enclone_data BCR=123089 CDR3=CARRYFGVVADAFDIW");
        if !plain {
            print!("{}", include_str!("example1"));
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
            print!("{}", strme(&x));
        }
        println!(
            "{}",
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
             • The notation 181.1.1 says that this V reference sequence is an alternate allele\n  \
             derived from the universal reference sequence (contig in the reference file)\n  \
             numbered 181, that is from donor 1 (\"181.1\") and is alternate allele 1 for that \
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
             • umed = median UMI count for a chain in the exact subclonotype.\n\
             • const = const region name for a chain in the exact subclonotype.\n"
        );
        print(
            "The view you see here is configurable: see the documentation at \
             \\bold{enclone help lvars} and \\bold{enclone help cvars}.  See also \
             \\bold{enclone help command} for how to remove the \\bold{PRE} part of the \
             command.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide example2 help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "example2") || help_all {
        begin_doc!("example2");
        println!(
            "\nSuppose you have placed the datasets that enclone comes with in the\n\
             directory /users/jdoe/enclone_data.  Then you can run this command:"
        );
        println!(
            "\n% enclone PRE=/users/jdoe/enclone_data BCR=123085 GEX=126106 \
             LVARSP=gex_med,IGHV2-5_g,CD4_a CDR3=CALMGTYCSGDNCYSWFDPW"
        );
        if !plain {
            print!("{}", include_str!("example2"));
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
            print!("{}", strme(&x));
        }
        print(
            "This shows an invocation of enclone that takes VDJ, gene expression and feature \
             barcode data as input, and exhibits all clonotypes for which some chain has the \
             given CDR3 sequence.  As well the command requests UMI (molecule) counts for one \
             hand-selected gene and one antibody.  You can use any gene(s) you like and any \
             antibodies for which you have feature barcodes.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide support help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "support") || help_all {
        begin_doc!("support");
        print(
            "\n\\red{enclone (beta) is provided as an open-source tool for use by the community.  \
             Although we cannot guarantee full support for the software, please email us at \
             enclone@10xgenomics.com if you have problems, questions or comments (see below).  \
             If you prefer you may submit a GitHub issue.}\n\n\
             \\blue{Please note that syntax and features in enclone will change over time.  See}\n\
             \\green{https://github.com/10XDev/enclone/blob/master/HISTORY.md} \
             \\blue{for the history of what was changed}\n\
             \\blue{and when.  We will try not to break} \
             \\blue{things, but when we first introduce a feature, it may}\n\
             \\blue{not be just right, and we may have to to perturb it later.}\n\n\
             \\bold{Critical things we want to hear from you}\n\n\
             1. If you have trouble understanding the documentation.\n\n\
             2. If enclone crashes.  We always need to see the output you got.  Often we will \
             need data to reproduce the problem.  Please also send this version information:\n",
        );
        println!("{} = {}.\n", env!("CARGO_PKG_VERSION"), VERSION_STRING);
        print(
            "3. If you're sure that enclone made a mistake.  Usually an actionable mistake is \
             exhibited via a single clonotype or two, along with your explanation as to why it's \
             wrong!  And we may need the enclone input \
             files to reproduce the problem.\n\n\
             4. Your ideas about how to make enclone better.\n\n\
             \\bold{Things where we may have trouble}\n\n\
             1. If you can't get your terminal window to work or need help working from the \
             command line.  Please tell us, but probably \
             you'll need help from someone who can interact with you on site.\n\n\
             2. If you get different results by two methods or observe suspicious statistics.  \
             Feel free to tell us, but we are much more likely to be able to respond if you have \
             specific data as in point 3 above.\n\n\
             3. Details about the algorithms.  We can answer some questions but may need to \
             point you to the source code to read yourself.\n\n\
             4. If for whatever reason, we get stuck.  We may not be able to fix every problem, \
             even if it's our fault.  And sometimes the right solution will take time.  We'll do \
             our best!\n\n",
        );
        print(
            "\\red{Please communicate with us by writing to} \
             \\boldred{enclone@10xgenomics.com.}\n\n\
             We use only this email address for enclone questions.  Support from 10x Genomics \
             for enclone needs to go through this point of access (or the next one).\n\n",
        );
        print(
            "If you are GitHub-savvy, we also welcome GitHub issues if applicable and you \
             prefer that route.  If you can make the enclone code better, go for it!  Make sure \
             the internal tests are successful, submit a \
             pull request, and if we can use your code (no guarantees), we'll add you as an \
             author on the site.\n\n",
        );
        print(
            "\\blue{Before writing to us, please check our faq by typing} \
             \\boldblue{enclone help faq}\\blue{.}\n\n",
        );
        print("\\red{Happy encloning!}\n\n");
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide plot help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "plot") || help_all {
        begin_doc!("plot");
        println!("");
        bold!();
        println!("plotting clonotypes\n");
        end_escape!();
        print(
            "enclone can create a \"honeycomb\" plot showing each clonotype as a cluster of \
             dots, one per cell.  You can see an example at \
             \\green{https://github.com/10XDev/enclone/blob/master/README.md#honeycomb}.\n\n\
             To generate such a plot, use the argument\n\
             \\bold{PLOT=\"filename,sample1->color1,...,samplen->colorn\"}\n\
             which creates an svg file of the given name, and assigns the given colors to the \
             given samples.  Unspecified samples will be black.\n\n\
             \
             The colors should be valid colors \
             for use in an svg file.  They can be named colors like red or blue (see \
             \\green{https://www.w3.org/TR/SVG11/types.html#ColorKeywords} for a full list) \
             or a hex specification like #00FFFF for aqua.  The full color description for svg \
             is at\n\
             \\green{https://www.w3.org/TR/SVGColor12}.\n\n\
             \
             Each cell is shown as a small \
             disk having the given color, and each clonotype is shown as a cluster of these small \
             disks, which are positioned at random.  We suggest using the \
             \\bold{MIN_CELLS} option \
             (see \"enclone help filter\") so that tiny clonotypes do not dominate.  The filename \
             argument may be \"stdout\".  Note that plotting is potentially slow.\n\n\
             To add a legend to the graph, add the argument \\bold{LEGEND} to your command \
             line.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide input help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "input") || help_all {
        begin_doc!("input");
        print(
            "\nenclone has \\boldred{two} mechanisms for specifying input datasets: either \
             directly on the command line or via a supplementary metadata file. Only one mechanism \
             may be used at a time.\n\n\
             In both cases, you will need to provide paths to directories where the outputs of \
             the Cell Ranger pipeline may be found.  enclone uses only some of the pipeline \
             output files, so it is enough that those files are present in given directory, and \
             the particular files that are needed may be found by typing \
             \\bold{enclone help input_tech}.  \
             If you use the argument \\bold{PRE=p} then \\bold{p/} will be prepended to all \
             pipeline paths.  Moreover (see \\bold{enclone help command}), you can avoid putting \
             \\bold{PRE} on the command line by setting the environment variable \
             \\bold{ENCLONE_PRE} to \\bold{p}.\n\n",
        );
        print(
            "Both input forms involve abbreviated names (discussed below), which should be as \
             short as possible, as longer abbreviations will increase the width of the clonotype \
             displays.\n\n",
        );
        print(
            "┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n\
             ┃ enclone can use gene expression and feature barcode data, as represented by a feature matrix.  ┃\n\
             ┃ Cell Ranger stores this matrix in an hdf5 file, which while generally very efficient, is not   ┃\n\
             ┃ optimized for interactive use.  Therefore enclone provides an alternate file structure, which  ┃\n\
             ┃ speeds up enclone overall by up to \\boldred{50%}.  To use this, add the argument \\bold{NH5} to the enclone      ┃\n\
             ┃ command line.  This will work so long as you have write permission on input directories.  The  ┃\n\
             ┃ first time you run enclone (using given inputs), the alternate file will be written; then      ┃\n\
             ┃ subsequent invocations will be faster.   Like with other enclone command-line options, if you  ┃\n\
             ┃ you want it on all the time, you can set the environment variable \\bold{ENCLONE_NH5}.                 ┃\n\
             ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛\n\n",
        );
        print(
            "\\boldred{█ 1 █} To point directly at input files on the command line, use e.g.\n\
             \\bold{TCR=/home/jdoe/runs/sample345}\n\
             or likewise for \\bold{BCR}.  A more complicated syntax is allowed in which commas, \
             colons and semicolons act as delimiters.  Commas go between datasets from the \
             same sample, colons between datasets from the same donor, and semicolons separate \
             donors.  If semicolons are used, the value must be quoted.\n\n",
        );
        print(
            "Using this input system, each dataset is assigned an abbreviated name, which is \
             everything \
             after the final slash in the directory name (e.g. \\bold{sample345} in the above \
             example), or the entire name if there is no slash; \
             samples and donors are assigned identifers s1,... and d1,..., respectively.\n\n",
        );
        print(
            "Examples:\n\
             \\bold{TCR=p1,p2}   -- input data from two libraries from the same sample\n\
             \\bold{TCR=p1,p2:q} -- input data as above plus another from a different sample \
             from the same donor\n\
             \\bold{TCR=\"a;b\"}   -- input one library from each of two donors.\n\n",
        );
        print(
            "Matching gene expression and/or feature barcode data may also be supplied using \
             an argument \\bold{GEX=...}, whose right side must have the exact same structure \
             as the \\bold{TCR} or \\bold{BCR} argument.  Specification of both \
             \\bold{TCR} and \\bold{BCR} is not allowed.\n\n",
        );
        print("\\boldred{█ 2 █} To specify a metadata file, use the command line argument\n");
        bold!();
        println!("META=filename");
        end_escape!();
        print(
            "This file should be a CSV (comma-separated values) file, with one line per cell \
             group.  After the first line, lines starting with # are ignored.  There must be a \
             field tcr or bcr, and some other fields are allowed:\n",
        );
        let mut log1 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log1);
        }
        log1.append(&mut b"field".to_vec());
        if !plain {
            emit_end_escape(&mut log1);
        }
        let s1 = stringme(&log1);
        let mut log2 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log2);
        }
        log2.append(&mut b"default".to_vec());
        if !plain {
            emit_end_escape(&mut log2);
        }
        let s2 = stringme(&log2);

        let mut log3 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log3);
        }
        log3.append(&mut b"meaning".to_vec());
        if !plain {
            emit_end_escape(&mut log3);
        }
        let s3 = stringme(&log3);
        rows.push(vec![s1, s2, s3]);
        rows.push(vec!["\\hline".to_string(); 3]);
        rows.push(vec![
            "tcr".to_string(),
            "(required!)".to_string(),
            "path to dataset, or abbr:path, where abbr is an abbreviated".to_string(),
        ]);
        rows.push(vec![
            "or bcr".to_string(),
            "".to_string(),
            "name for the dataset; exactly one of tcr or bcr must be used".to_string(),
        ]);
        rows.push(vec![
            "gex".to_string(),
            "null".to_string(),
            "path to GEX dataset, which may include or consist entirely".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "of FB data".to_string(),
        ]);
        rows.push(vec!["\\hline".to_string(); 3]);
        rows.push(vec![
            "sample".to_string(),
            "s1".to_string(),
            "abbreviated name of sample".to_string(),
        ]);
        rows.push(vec!["\\hline".to_string(); 3]);
        rows.push(vec![
            "donor".to_string(),
            "d1".to_string(),
            "abbreviated name of donor".to_string(),
        ]);
        rows.push(vec!["\\hline".to_string(); 3]);
        rows.push(vec![
            "bc".to_string(),
            "null".to_string(),
            "name of CSV file with header \"barcode,sample,donor\" that".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "assigns a sample and donor name to each barcode; if sample and/or".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "donor are also specified, then those are treated as default values".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "to be used in case a particular barcode is not specified by bc;".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "a fourth field \"tag\" is allowed and may be arbitrarily specified".to_string(),
        ]);
        rows.push(vec![
            "".to_string(),
            "".to_string(),
            "Use of bc automatically turns on the MIX_DONORS option".to_string(),
        ]);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 3, &b"l|l|l".to_vec(), false);
        println!("{}", log);
        if !help_all {
            std::process::exit(0);
        }
    }
}
