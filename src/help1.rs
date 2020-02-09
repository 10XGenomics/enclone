// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Test for help request.

use crate::help_utils::*;
use crate::misc1::*;
use ansi_escape::*;
use pretty_trace::*;
use std::env;
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help1(args: &Vec<String>) {
    // Set up.

    let mut args = args.clone();
    if args.len() == 1 || (args.len() >= 2 && args[1] == "help") {
        PrettyTrace::new().on();
        let mut nopager = false;
        let mut to_delete = vec![false; args.len()];
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                nopager = true;
                to_delete[i] = true;
            }
        }
        erase_if(&mut args, &to_delete);
        setup_pager(!nopager);
    }
    let mut help_all = false;
    if args.len() >= 3 && args[1] == "help" && args[2] == "all" {
        unsafe {
            HELP_ALL = true;
        }
        help_all = true;
    }
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
    macro_rules! doc_red {
        ($n1:expr, $n2:expr) => {
            if !plain {
                let r1 = format!("[01;31m{}[0m", $n1);
                let r2 = format!("[01;31m{}[0m", $n2);
                rows.push(vec![r1, r2]);
            } else {
                rows.push(vec![$n1.to_string(), $n2.to_string()]);
            }
        };
    }
    macro_rules! ldoc_red {
        ($n1:expr, $n2:expr) => {
            rows.push(vec!["\\hline".to_string(); 2]);
            if !plain {
                let r1 = format!("[01;31m{}[0m", $n1);
                let r2 = format!("[01;31m{}[0m", $n2);
                rows.push(vec![r1, r2]);
            } else {
                rows.push(vec![$n1.to_string(), $n2.to_string()]);
            }
        };
    }
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
    macro_rules! begin_doc {
        ($x:expr) => {
            rows.clear();
            if help_all {
                banner($x, plain);
            }
        };
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Setup test.

    if (args.len() == 2 && args[1] == "help")
        || (args.len() == 2 && args[1] == "--help")
        || help_all
    {
        if help_all {
            println!("");
        }
        begin_doc!("");
        print(
            "\nWelcome to enclone!\n\n\
             The purpose of this first page is to help you make sure that you're set up properly\n\
             to run enclone.  PLEASE READ!\n\n\
             (for the main help page, please type instead: enclone help main)\n\n\n\
             Here we go through several setup tests.  If you have any problem that you can't\n\
             resolve, please email us at enclone@10xgenomics.com.\n\n\n\
             1. Are you using a fixed width font?\n\
             Look at this:\n\
             A FAT BROWN CAT JUMPED OVER THE WALL\n\
             ||||||||||||||||||||||||||||||||||||\n\
             Do those two lines end at the same position?  If not, you need to switch your \
             font.\n\n\n\
             2. Is your terminal window wide enough to see the help pages?\n\
             Your terminal needs to be at least 100 columns wide.  Look at this:\n\
             01234567890123456789012345678901234567890123456789\
             01234567890123456789012345678901234567890123456789\n\
             Does it appear as a single line?  If not, please widen your window.\n\n\n\
             3. Can your terminal display box characters?\n\
             Look at this:\n\
             ┌────────┬─────────┐\n\
             │banana  │  peel   │\n\
             ├────────┼─────────┤\n\
             │oops    │  slipped│\n\
             └────────┴─────────┘\n\
             Do you see a neat rectangle composed of four rectangles with words inside them?\n\
             If not, something is wrong with your terminal!\n\n\
             4. Can your terminal correctly display ANSI escape sequences?\n\
             The following word should be \\bold{bold}.  \
             The following word should be \\blue{blue}.\n\
             If that doesn't make sense, or is messed up, something is wrong, and you have \
             two options:\n\
             (a) seek help to fix your terminal window\n\
             (b) turn off escape sequences by adding PLAIN to every enclone command, or set\n\
             the environment variable ENCLONE_PLAIN.\n\
             But that should be only a last resort.\n\n\
             5. Can your terminal correctly display unicode characters?\n\
             Do you see a centered dot here • ?\n\
             If not, your terminal has a problem!\n\n\
             6. Does this entire help page appear at once in your terminal window?\n\
             If not, please increase the number of rows in your window to 56.\n\n\n\
             If you go through all those tests and everything worked, you should be \
             good to go!\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // Provide main help.

    if args.len() == 1 || (args.len() == 3 && args[1] == "help" && args[2] == "main") || help_all {
        begin_doc!("main");
        print!("\nThis is version {} (beta) of ", env!("CARGO_PKG_VERSION"));
        print_enclone(plain);
        print!(".  The mission of ");
        print_enclone(plain);
        println!(" is to:\n");
        bold!();
        println!("  Find and display the clonotypes within single cell VDJ datasets:");
        println!("  groups of cells having the same fully rearranged common ancestor.\n");
        end_escape!();
        print(
            "This help page catalogs all the enclone help pages.  We strongly \
             recommend studying at least those in \\red{red} below.  \
             Pages fit in 100 wide x 56 high \
             windows, except those labeled \"long\" or \"wide\".\n\n",
        );
        print(
            "\\boldblue{enclone is part of the 10x Genomics immune profiling tools, including \
             Cell Ranger and Loupe,}\n\\boldblue{which enclone is integrated with.}  enclone \
             uses output from Cell Ranger version \\boldred{≥ 3.1.}\n\n",
        );
        let mut log1 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log1);
        }
        log1.append(&mut b"command".to_vec());
        if !plain {
            emit_end_escape(&mut log1);
        }
        let s1 = stringme(&log1);
        let mut log2 = Vec::<u8>::new();
        if !plain {
            emit_bold_escape(&mut log2);
        }
        log2.append(&mut b"what it provides".to_vec());
        if !plain {
            emit_end_escape(&mut log2);
        }
        let s2 = stringme(&log2);
        doc!(&s1, &s2);
        ldoc_red!("enclone help", "help to test for correct setup");
        doc_red!("enclone", "what you see here: guide to all the doc");
        ldoc_red!("enclone help quick", "quick guide to getting started");
        ldoc_red!("enclone help how", "how enclone works (long)");
        ldoc_red!(
            "enclone help command",
            "info about enclone command line argument processing"
        );
        ldoc_red!(
            "enclone help glossary",
            "glossary of terms used by enclone, and conventions"
        );
        ldoc_red!("enclone help example1", "explanation of an example");
        doc_red!(
            "enclone help example2",
            "example showing gene expression \
             and feature barcodes (wide)"
        );
        ldoc_red!(
            "enclone help support",
            "how we can help, enclone@10xgenomics.com"
        );
        ldoc!("enclone help input", "how to provide input to enclone");
        doc!(
            "enclone help input_tech",
            "how to provide input to enclone (technical notes)"
        );
        ldoc!("enclone help parseable", "parseable output (long)");
        ldoc!("enclone help filter", "clonotype filtering options");
        doc!("enclone help special", "special filtering options (long)");
        ldoc!("enclone help lvars", "lead column options");
        doc!("enclone help cvars", "per chain column options");
        doc!(
            "enclone help amino",
            "per chain column options for amino acids"
        );
        doc!("enclone help display", "other clonotype display options");
        ldoc!("enclone help indels", "insertion and deletion handling");
        ldoc!(
            "enclone help color",
            "how enclone uses color, and related things"
        );
        doc!(
            "enclone help ideas",
            "ideas for features that might be implemented"
        );
        doc!("enclone help faq", "frequently asked questions (long)");
        doc!(
            "enclone help all",
            "concatenation of all the help pages (long, wide);"
        );
        doc!("", "you can use this to search all the help pages");
        print_tab2(&rows);
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide quick help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "quick") || help_all {
        begin_doc!("quick");
        println!("");
        bold!();
        println!("quick guide to getting started\n");
        end_escape!();
        print(
            "Just type this:\n\n\
             \\bold{enclone BCR=p}\n\n\
             where \\bold{p} is the path to your Cell Ranger VDJ directory.\n\n\
             Substitute \\bold{TCR} if that's what you've got.\n\n\
             This will show you all the clonotypes, in descending order by number of cells.\n\n\
             You'll need to make your window wide enough so that lines are not folded.  This \
             depends on the dataset.\n\n\
             Only one page of output is shown at a time.  To navigate within the full output, \
             use the space bar to go forward and the b key to go backward.\n\n\
             See \\bold{enclone help example1} for a detailed guide to how to read the enclone \
             output.  A few key things you should know:\n\n\
             1. You'll see numbers near the top.  These are amino acid position numbers, and\n   \
             they read downwards.  Numbering starts at the start codon, numbered zero.\n\n\
             2. Each numbered line represents an exact subclonotype: cells having identical \
             V(D)J transcripts.\n\n\
             3. By default, you'll see data in amino acid space.  Only \"interesting\" amino acids \
             are shown.\n\n\
             Please read on to learn more!\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide how help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "how") || help_all {
        // Start.

        begin_doc!("how");
        println!("");
        bold!();
        println!("information about how enclone works\n");
        end_escape!();
        print(
            "The goal of enclone is to find and display the clonotypes within single cell \
             VDJ datasets: groups of cells having the same fully rearranged common ancestor.\n\n\
             \
             enclone provides the foundation for fully understanding each cell's antigen \
             affinity and the evolutionary relationship between cells within a sample.  \
             This starts with, for each cell, \
             \\bold{the full length sequence of all its VDJ receptor chains}.  Such data may be \
             obtained using the 10x Genomics immune profiling platform.\n\n\
             \
             For this, there are fundamental challenges:\n\n",
        );

        // Print challenges.

        print_with_box(
            "1. It is extremely easy to get false positives: the incorrect \
             appearance that two cells have a common ancestor.\n\n\
             \
             2. Because of somatic hypermutation in B cells, it can be difficult to know that \
             two B cells share a common ancestor.\n\n\
             \
             3. There is always some background noise, e.g. from ambient mRNA.  When building \
             large clonotypes, this noise tends to pile up, yielding ectopic chains, i.e. chains \
             within a clonotype that are artifacts and do not represent true biology.",
        );

        // Print boxed algorithm.

        print(
            "To address these challenges, the enclone algorithm has a number of steps, which we \
             outline:\n\n",
        );
        print(
            "\\boldred{1}. enclone gets its information from the file all_contig_annotations.json \
             that is \
             produced by Cell Ranger.  Only productive contigs are used.  Each has an annotated \
             V and J segment.  The V segment alignment may have a single indel whose length is \
             divisible by three, and in that case, the V reference sequence is edited either to \
             delete or insert sequence.  In the insertion case, the bases are taken from the \
             contig.  These indels are noted in the enclone output.\n\n\
             \
             \\boldred{2}. enclone groups cells into exact subclonotypes, provided that they have \
             the same \
             number of chains, identical V..J sequences, identical C segment assignments, \
             and the same distance between the J stop and the C start (which is usually zero).\n\n\
             \
             \\boldred{3}. For samples from a given donor, enclone derives \"donor reference \
             sequences\" \
             for the V chains present in the donor's genome.  This is powerful, even though \
             based on imperfect information.  V segments vary in their expression frequency and \
             thus the more cells which are present, the more complete the information will be.  It \
             is also not possible to accurately determine the last ~15 bases in a V chain from \
             transcript data.\n\n\
             \
             \\boldred{4}. Pairs of exact subclonotypes are considered for joining, as described \
             below.  This \
             process only considers exact subclonotypes have two or three chains.  There is some \
             separate joining for the case of one chain.  Exact subclonotypes having four chains \
             are not joined at present.  These cases are clearly harder because these exact \
             subclonotypes are highly enriched for cell doublets, which we discard if we can \
             identify as such.\n\n\
             \
             \\boldred{5}. For each pair of exact subclonotypes, and for each pair of chains in \
             each of the \
             two exact subclonotypes, for which V..J has the same length for the corresponding \
             chains, and the CDR3 segments have the same length for the corresponding chains, \
             enclone considers joining the exact subclonotypes into the same clonotype.\n\n\
             \
             \\boldred{6}. To proceed, as a minimum requirement, there must be at most 50 total \
             mismatches \
             between the two exact subclonotypes, within the given two V..J segments.\n\n\
             \
             \\boldred{7}. enclone next finds shared differences betweens exact subclonotypes, \
             that is, for \
             two exact subclonotypes, common mutations from the reference sequence, using the \
             donor reference for the V segments and the universal reference for the J segments. \
             We also compute the number of independent mutations.\n\n\
             \
             \\boldred{8}. Two cells sharing sufficiently many shared differences and \
             sufficiently few \
             CDR3 differences are deemed to be in the same clonotype.  This depends on heuristics \
             which are too detailed to describe on this page and will be described elsewhere.\n\n\
             \
             \\boldred{9}. Spurious chains are filtered out based on frequency and \
             connections.\n\n",
        );

        // Finish.

        print( "We are actively working to improve the algorithm.  To test the performance of the \
            current version, we combined data from 443 BCR libraries from 30 donors, which yielded \
            \\boldred{9573} clonotypes having at least two cells each, of \
            which \\boldred{15 (0.16%)} \
            contained data from multiple donors.  These are errors.\n\n" );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide command line help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "command") || help_all {
        begin_doc!("command");
        println!("");
        bold!();
        println!("information about enclone command-line argument processing\n");
        end_escape!();
        print("\\bold{1. Order of processing}\n\n");
        print(
            "• Before processing its command line, enclone first checks for environment\n\
             variables of the form \\bold{ENCLONE_<x>}.  These are converted into command-line \
             arguments.\n\
             • For example, setting the environment variable \\bold{ENCLONE_PRE} to \
             \\bold{/Users/me/enclone_data} \
             is equivalent to providing the command-line argument \
             \\bold{PRE=/Users/me/enclone_data}.\n\
             • After checking environment variables, arguments on the command line are read from \
             left to right; if an argument name is repeated, only the \
             rightmost value is used.\n\n",
        );
        print("\\bold{2. Color}\n\n");
        print_enclone(plain);
        print(
            " uses ANSI escape codes for color and bolding, frivolously, for emphasis, \
             and more\nimportantly for amino acids, to represent different codons.\n\n\
             \\boldred{PLEASE READ THIS:}\n\n\
             You can turn off escape codes by adding \\bold{PLAIN} to any command.  Use this if \
             you want to peruse output using a text editor which does not grok the escape \
             codes.  However some things will not make sense without color.\n\n",
        );
        print("\\bold{3. Paging}\n\n");
        print("• enclone automatically pipes its output to \\bold{less -R -F -X}.\n");
        print(
            "• The effect of this will be that you'll see only the first screen of output.  \
             You can then use the spacebar to go forward, b to go backward, and q to quit.  \
             The \\bold{-R} option causes escape characters to be correctly displayed, the \
             \\bold{-F} option causes an automatic exit if output fits on a single screen, and \
             the \\bold{-X} option prevents output from being sent to the \"alternate screen\" \
             under certain platform/version combinations.\n",
        );
        print("• Type \\bold{man less} if you need more information.\n");
        print(
            "• If for whatever reason you need to turn off output paging, add the argument \
             \\bold{NOPAGER} to the enclone command.\n\n",
        );
        if !help_all {
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide glossary help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "glossary") || help_all {
        begin_doc!("glossary");
        println!("");

        // intro

        bold!();
        println!("glossary of terms used by enclone\n");
        end_escape!();

        // doc V..J

        doc!(
            "V..J",
            "the full sequence of a V(D)J transcript, from the beginning of the V"
        );
        doc!(
            "",
            "segment to the end of the J segment; this sequence begins with a stop codon"
        );
        doc!("", "and ends with a partial codon (its first base)");

        // doc clonotype

        ldoc!(
            "clonotype",
            "all the cells descended from a single fully rearranged T or B cell"
        );
        doc!("", "(approximated computationally)");

        // doc exact subclonotype

        let x1 = "exact subclonotype".to_string();
        let mut w2 = b"all cells having identical transcripts".to_vec();
        emit_bold_escape(&mut w2);
        emit_red_escape(&mut w2);
        w2.append(&mut " ○".as_bytes().to_vec());
        emit_end_escape(&mut w2);
        let x2 = stringme(&w2);
        rows.push(vec![x1, x2]);
        doc!("", "(every clonotype is a union of exact subclonotypes)");

        // doc clone

        doc!(
            "clone",
            "a cell in a clonotype, or in an exact subclonotype"
        );

        // doc onesie etc.

        ldoc!(
            "onesie",
            "a clonotype or exact subclonotype having exactly one chain"
        );
        doc!(
            "twosie",
            "a clonotype or exact subclonotype having exactly two chains"
        );
        doc!(
            "threesie",
            "a clonotype or exact subclonotype having exactly three chains;"
        );
        doc!(
            "",
            "these frequently represent true biological events, arising from expression"
        );
        doc!("", "of both alleles");
        doc!(
            "foursie",
            "a clonotype or exact subclonotype having exactly four chains;"
        );
        doc!("", "these very rarely represent true biological events");
        doc!("moresie", "a clonotype having more than four chains;");
        doc!(
            "",
            "these sad clonotypes do not represent true biological events"
        );

        // doc donor etc.

        ldoc!("donor", "an individual from whom samples are obtained");
        doc!(
            "sample",
            "a tube of cells from a donor, from a particular tissue at a"
        );
        doc!(
            "",
            "particular point in time, and possibly enriched for particular cells"
        );
        doc!(
            "cell group",
            "an aliquot from a sample, presumed to be a random draw"
        );
        doc!(
            "dataset",
            "all sequencing data obtained from a particular library type"
        );
        doc!(
            "",
            "(e.g. TCR or BCR or GEX or FB), from one cell group, processed by running"
        );
        doc!("", "through the Cell Ranger pipeline");

        // print main table

        print_tab2(&rows);

        // print footnote

        print(
            "\\boldred{○} The exact requirements for being in the same exact subclonotype are \
             that cells:\n\
             • have the same number of productive contigs identified\n\
             • that these have identical bases within V..J\n\
             • that they are assigned the same constant region reference sequences\n\
             • and that the difference between the V stop and the C start is the same\n  \
             (noting that this difference is nearly always zero).\n\
             Note that we allow mutations within the 5'-UTR and constant regions.\n\n",
        );

        // conventions

        bold!();
        println!("conventions\n");
        end_escape!();
        println!(
            "• When we refer to \"V segments\", we always include the leader segment.\n\
             • Zero or one?  We number exact subclonotypes as 1, 2, ... and likewise with\n\
             chains within a clonotype, however DNA and amino-acid positions are numbered starting \
             at zero.\n"
        );

        // done

        if !help_all {
            std::process::exit(0);
        }
    }
}
