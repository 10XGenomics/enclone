// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test for help request.

use crate::help_utils::*;
use enclone_core::defs::*;
use enclone_core::testlist::*;
use itertools::Itertools;
use string_utils::*;

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
            h.print_plain(&format!("{}", strme(&x)))?;
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
                ex2_args2.push(ex2_args[i].clone());
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
            h.print_plain_unchecked(&format!("{}", strme(&x)));
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
             is\n\\bold{~/enclone/datasets,~/enclone/datasets2}.",
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
             origins and donors are assigned identifers s1,... and d1,..., respectively; \
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
        h.print("\\boldred{█ 2 █} To specify a metadata file, use the command line argument\n")?;
        h.print("\\bold{META=filename}\n")?;
        h.print(
            "This file should be a CSV (comma-separated values) file, with one line per cell \
             group.  After the first line, lines starting with # are ignored.  There must be a \
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
            "\nIn addition, metadata maybe fully specified on the command line via \
            \\bold{METAX=\"l1;...;ln\"} where the \\bold{li} are the lines that you would \
            otherwise put in the \\bold{META} file.\n",
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
    Ok(())
}
