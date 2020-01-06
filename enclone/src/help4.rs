// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use ansi_escape::*;
use help_utils::*;
use std::env;
use string_utils::*;
use tables::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help4() {

    // Set up.

    let mut args: Vec<String> = env::args().collect();
    let mut rows = Vec::<Vec<String>>::new();
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
                print!( "{}", strme(&log) );
            }
        };
    }
    macro_rules! end_escape {
        () => {
            if !plain {
                let mut log = Vec::<u8>::new();
                emit_end_escape(&mut log);
                print!( "{}", strme(&log) );
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

    if ( args.len() == 3 && args[1] == "help" && args[2] == "special" ) || help_all {
        begin_doc!("special");
        print( "\n\\bold{special filtering options}\n\n" );
        print( "This page documents some options, most of which allow noise \
            filters to be turned off, and which normally should not be invoked.  The last \
            two options can be used to simplify the view of a clonotype.\n\n" );
        doc!( "NCROSS", 
            "If you specify that two or more libraries arose from the same sample (i.e." );
        doc!( "", "from the same tube of cells), then the default behavior of enclone is to" );
        doc!( "", "\"cross filter\" so as to remove expanded exact subclonotypes that are" );
        doc!( "", "present in one library but not another, in a fashion that would be highly" );
        doc!( "", "improbable, assuming random draws of cells from the tube.  These are" );
        doc!( "", "believed to arise when a plasma or plasmablast cell breaks up during during" );
        doc!( "", "or after pipetting from the tube, and the resulting fragments seed GEMS," );
        doc!( "", "yielding 'fake' cells.  The NCROSS option turns off this filtering step." );
        ldoc!( "NGRAPH_FILTER", 
            "By default, enclone filters to remove exact subclonotypes that by virtue of" );
        doc!( "", "their relationship to other exact subclonotypes, appear to arise from" );
        doc!( "", "background mRNA or a phenotypically similar phenomenon.  The" );
        doc!( "", "NGRAPH_FILTER option turns off this filtering." );
        ldoc!( "NQUAL",
            "By default, enclone filters out exact subclonotypes having a base in V..J" );
        doc!( "", "that looks like it might be wrong.  More specifically, enclone finds bases" );
        doc!( "", "which are not Q60 for a barcode, not Q40 for two barcodes, are not" );
        doc!( "", "supported by other exact subclonotypes, are variant within the clonotype," );
        doc!( "", "and which disagree with the donor reference.  NQUAL turns this off." );
        ldoc!( "NWEAK_CHAINS", "By default, enclone filters chains from clonotypes that are" );
        doc!( "", "weak and appear to be artifacts, perhaps arising from a stray mRNA molecule" );
        doc!( "", "that floated into a GEM.  The NWEAK_CHAINS option turns off this filter." );
        ldoc!( "NWEAK_ONESIES", 
            "By default, enclone filters out onesie clonotypes having a single exact" );
        doc!( "", "subclonotype, and that are light chain or TRA, and whose number of cells is" );
        doc!( "", "less than 0.1% of the total number of cells.  NWEAK_ONESIES turns this off." );
        ldoc!( "NFOURSIE_KILL", "By default, if enclone finds a foursie exact subclonotype that" );
        doc!( "", "contains a twosie exact subclonotype having at least ten cells, it kills" );
        doc!( "", "the foursie exact subclonotype, no matter how many cells it has.  The" );
        doc!( "", "foursies that are killed are believed to be rare oddball artifacts arising" );
        doc!( "", "from repeated cell doublets or GEMs that contain two cells and multiple gel" );
        doc!( "", "beads.  The argument NFOURSIE_KILL turns off this filtering." );
        ldoc!( "NWHITEF", "By default, enclone filters out rare artifacts arising from \
                contamination" );
        doc!( "", "of oligos on gel beads.  The NWHITEF option turns off this filter." );
        ldoc!( "KEEP_IMPROPER", "An exact subclonotype is improper if it does not have one chain" );
        doc!( "", "of each type.  This option causes all improper exact subclonotypes to be" );
        doc!( "", "retained, although they may be removed by other filters." );
        ldoc!( "MIN_CHAINS_EXACT=n",
            "Delete any exact subclonotype having less than n chains.  You can use this" );
        doc!( "", "to \"purify\" a clonotype so as to display only exact subclonotypes having" );
        doc!( "", "all their chains." );
        doc!( "MIN_CELLS_EXACT=n", 
            "Delete any exact subclonotype having less than n cells.  You might want" );
        doc!( "", "to use this if you have a very large and complex expanded clonotype," );
        doc!( "", "for which you would like to see a simplified view." );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide lvars help.

    if ( args.len() == 3 && args[1] == "help" && args[2] == "lvars" ) || help_all {
        begin_doc!("lvars");
        println!("");
        bold!();
        println!( "lead column options\n" );
        end_escape!();
        println!( "These options define lead variables, which correspond to columns that\n\
            appear once in each clonotype, on the left side, and have one entry for each\n\
            exact subclonotype row.\n" );
        print( "Lead variables are specified using\n\
                \\bold{LVARS=x1,...,xn}\n\
                where each xi is one of:\n\n" );
        doc!( "datasets", "dataset identifiers" );
        doc!( "donors", "donor identifiers" );
        ldoc!( "ncells", "number of cells" );
        doc!( "n_<name>", 
            "number of cells associated to the given name, which can be a dataset short" );
        doc!( "", 
            "name, or a sample short name, or a donor short name; it may not name more than" );
        doc!( "", "one such category" );
        ldoc!( "gex_med", "median gene expression UMI count" );
        doc!( "gex_max", "max gene expression UMI count" );
        doc!( "", "(requires that gene expression data are provided as input)" );
        ldoc!( "near", "Hamming distance of V..J DNA sequence to nearest neighbor" );
        doc!( "far", "Hamming distance of V..J DNA sequence to farthest neighbor" );
        doc!( "", 
            "both only compare to cells having chains in the same columns of the clonotype;" );
        doc!( "", "with - shown if there is no other exact subclonotype to compare to" );
        ldoc!( "g<d>", 
            "Here d is a nonnegative integer.  Then all the exact subclonotypes are grouped" );
        doc!( "", "according to the Hamming distance of their V..J sequences.  Those within" );
        doc!( "", "distance d are defined to be in the same group, and this is extended" );
        doc!( "", "transitively.  The group identifier 1, 2, ... is shown.  The ordering of" );
        doc!( "", "these identifiers is arbitrary.  This option is best applied to cases where" );
        doc!( "", "all exact subclonotypes have a complete set of chains." );
        ldoc!( "<antibody>_a", "assuming that feature barcode data has been provided," );
        doc!( "", "look for a feature line that starts with the given name, and" );
        doc!( "", "then has a tab; report the mean UMI count value" );
        doc!( "<gene name>_g", "assuming that gene expression data has been provided," );
        doc!( "", "look for a feature line that has the given name in the second" );
        doc!( "", "tab-delimited column; report the mean UMI count value" );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        print( "The default is \\bold{datasets,ncells}, except that datasets is suppressed if \
            there is only one dataset.\n\n" );
        print( "\\bold{LVARSP=x1,...,xn} is like \\bold{LVARS} but appends to the default \
            list.\n\n" );
        print( "Note: gene expression counts are normalized to 20,000 read pairs per cell, and \
            feature barcode counts are normalized to 5,000 read pairs per cell.  The normalized \
            counts are rounded to the nearest integer.\n\n" );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide cvars help.

    if ( args.len() == 3 && args[1] == "help" && args[2] == "cvars" ) || help_all {
        begin_doc!("cvars");

        // Header.

        print( "\n\\bold{per-chain column options}: These options define per-chain variables, \
            which correspond to columns that appear once for each chain in each clonotype, and \
            have one entry for each exact subclonotype.\n\n" );
        print( "Per-column variables are specified using\n\
                \\bold{CVARS=x1,...,xn}\n\
                where each xi is one of:\n\n" );

        // Main table entries.

        doc!( "var", "bases at positions in chain that vary across the clonotype" );
        ldoc!( "umed", "median VDJ umi count for each exact subclonotype" );
        doc!( "umax", "max VDJ umi count for each exact subclonotype" );
        doc!( "utot", "total VDJ umi count for each exact subclonotype" );
        ldoc!( "rmed", "median VDJ read count for each exact subclonotype" );
        ldoc!( "const", "constant region name" );
        ldoc!( "comp", "a measure of CDR3 complexity (currently computed inefficiently)" );
        ldoc!( "cdr3_dna", "the CDR3_DNA sequence" );
        ldoc!( "clen", "length of observed constant sequence (usually truncated at primer start)" );
        doc!( "ulen", "length of observed 5'-UTR sequence" );
        doc!( "cdiff", "differences with universal reference constant region, shown in the" );
        doc!( "", "abbreviated form e.g. 22T (ref changed to T at base 22) or 22T+10" );
        doc!( "", "(same but contig has 10 additional bases beyond end of ref C region" );
        doc!( "", "At most five differences are shown, and if there are more, ... is appended." );
        doc!( "udiff", "like cdiff, but for the 5'-UTR" );
        ldoc!( "notes", "optional note if there is an insertion, elided if empty" );
        ldoc!( "ndiff<n>", "number of base differences within V..J with exact subclonotype n" );

        // The rest.

        let mut log = String::new();
        // was ... rows.clone()
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        print( "At least one variable must be listed.  The default is \\bold{umed,const,notes}.  \
            \\bold{CVARSP}: same as \\bold{CVARS} but appends.\n\n" );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide amino help.

    if ( args.len() == 3 && args[1] == "help" && args[2] == "amino" ) || help_all {
        begin_doc!("amino");
        print( "\nThere is a complex per-chain column to the left of other \
            per-chain columns, defined by\n\
            \\bold{AMINO=x1,...,xn}: display amino acid columns for the given categories, \
            in one combined ordered group, where each xi is one of:\n\n" );
        doc!( "cdr3", "CDR3 sequence" );
        ldoc!( "var", "positions in chain that vary across the clonotype" );
        doc!( "share", "positions in chain that differ consistently from the donor reference" );
        ldoc!( "donor", "positions in chain where the donor reference differs from the universal \
            reference" );
        ldoc!( "donorn", "positions in chain where the donor reference differs nonsynonymously" );
        doc!( "", "from the universal reference" );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        print( "Note that we compute positions in base space, and then divide by three to get \
            positions in amino acid space.  Thus it can happen that a position in amino acid \
            space is shown for both \\bold{var} and \\bold{share}.\n\n" );
        print( "The default value for \\bold{AMINO} is \\bold{cdr3,var,share,donor}.  \
            Note that we only report amino acids that are strictly within V..J, \
            thus specifically excluding the codon bridging J and C.\n\n" );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide display help.

    if ( args.len() == 3 && args[1] == "help" && args[2] == "display" ) || help_all {
        begin_doc!("display");
        print( "\n\\bold{other options that control clonotype display}\n\n" );
        doc!( "PER_BC", "expand out each exact clonotype line, showing one line per barcode," );
        doc!( "", "for each such line, displaying the barcode name, the number of UMIs assigned," );
        doc!( "", "and the gene expression UMI count, if applicable, under gex_med" );
        ldoc!( "BARCODES", "print list of all barcodes of the cells in each clonotype, in a" );
        doc!( "", "single line near the top of the printout for a given clonotype" );
        ldoc!( "SEQC", "print V..J sequence for each chain in the first exact subclonotype, near" );
        doc!( "", "the top of the printout for a given clonotype" );
        ldoc!( "FULL_SEQC", "print full sequence for each chain in the first exact subclonotype," );
        doc!( "", "near the top of the printout for a given clonotype" );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        print( "\\bold{options that control clonotype grouping}\n\n\
            we plan to add grouping capability in a future version of enclone, and for now, as a \
            placeholder, we have the following \"toy\" options:\n\n" );
        rows.clear();
        doc!( "GROUP_HEAVY_CDR3", "group by perfect identity of CDR3 amino acid sequence \
            of IGH or TRB" );
        ldoc!( "MIN_GROUP", "minimum number of clonotypes in group to print (default = 1)" );
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
        println!( "{}", log );
        if !help_all {
            std::process::exit(0);
        }
    }
}
