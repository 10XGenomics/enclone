// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test for help request.

use crate::help_utils::HelpDesk;

pub fn help1(args: &Vec<String>, h: &mut HelpDesk) -> Result<(), String> {
    // Provide main help.

    if args.len() == 1 || (args.len() == 3 && args[1] == "help" && args[2] == "main") || h.help_all
    {
        if h.help_all {
            h.print("\n")?;
        }
        h.begin_doc("")?;
        h.print("\nThe mission of ")?;
        h.print_enclone()?;
        h.print(" is to:\n\n")?;
        h.print("\\bold{  Find and display the clonotypes within single cell VDJ datasets:}\n")?;
        h.print("\\bold{  groups of cells having the same fully rearranged common ancestor.}\n\n")?;
        h.print(
            "\\boldblue{enclone is part of the 10x Genomics immune profiling tools, including \
             Cell Ranger and Loupe.  enclone uses output from Cell Ranger version ≥ 3.1.}\n\n\
             The complete enclone documentation is at \\green{bit.ly/enclone}.  This page \
             catalogs the subset of those pages that are directly accessible from the enclone \
             command line.  These pages can be viewed in a 100 wide x 56 high window, except for \
             those labeled \"long\" or \"wide\".\n\n",
        )?;
        h.docpr("\\bold{command}", "\\bold{what it provides}");
        h.ldoc("enclone help", "help to test for correct setup");
        h.doc("enclone", "what you see here: guide to all the doc");
        h.ldoc("enclone help quick", "quick guide to getting started");
        h.doc("enclone help how", "how enclone works (long)");
        h.doc(
            "enclone help command",
            "info about enclone command line argument processing",
        );
        h.ldoc(
            "enclone help glossary",
            "glossary of terms used by enclone, and conventions",
        );
        h.ldoc("enclone help example1", "explanation of an example");
        h.doc(
            "enclone help example2",
            "example showing gene expression \
             and feature barcodes (wide)",
        );
        h.ldoc(
            "enclone help input",
            "how to provide input to enclone (long)",
        );
        h.doc(
            "enclone help input_tech",
            "how to provide input to enclone (technical notes)",
        );
        h.doc("enclone help parseable", "parseable output (long)");
        h.ldoc(
            "enclone help filter",
            "clonotype filtering options, scanning for feature enrichment (long)",
        );
        h.doc("enclone help special", "special filtering options (long)");
        h.ldoc("enclone help lvars", "lead column options (long)");
        h.doc("enclone help cvars", "per chain column options (long)");
        h.doc(
            "enclone help amino",
            "per chain column options for amino acids",
        );
        h.doc(
            "enclone help display",
            "other clonotype display options (long)",
        );
        h.ldoc("enclone help indels", "insertion and deletion handling");
        h.ldoc(
            "enclone help color",
            "how enclone uses color, and related things",
        );
        h.doc("enclone help faq", "frequently asked questions (long)");
        h.ldoc_greenish(
            "enclone help all",
            "concatenation of all the help pages (long, wide)",
        );
        h.doc_greenish("", "███ USE THIS TO SEARCH ALL THE HELP PAGES! ███");
        h.print_tab2()?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Setup test.

    if (args.len() == 2 && args[1] == "help")
        || (args.len() == 2 && args[1] == "--help")
        || h.help_all
    {
        h.begin_doc("setup")?;
        h.print(
            "\n\nWelcome to enclone!\n\n\
             The purpose of this first page is to help you make sure that you're set up properly\n\
             to run enclone.  PLEASE READ!\n\n\
             (for the main help page, please type instead: enclone)\n\n\
             Here we go through several setup tests.  If you have any problem that you can't\n\
             resolve, please email us at enclone@10xgenomics.com.\n\n\n\
             1. Are you using a fixed width font?\n\
             Look at this:\n\
             A FAT BROWN CAT JUMPED OVER THE WALL\n\
             ||||||||||||||||||||||||||||||||||||\n\
             Do those two lines end at the same position?  If not, you need to switch your \
             font.\n\n\
             2. Is your terminal window wide enough to see the help pages?\n\
             Your terminal needs to be at least 100 columns wide.  Look at this:\n\
             01234567890123456789012345678901234567890123456789\
             01234567890123456789012345678901234567890123456789\n\
             Does it appear as a single line?  If not, please widen your window.\n\n\
             3. Can your terminal display box characters?\n\
             Look at this:\n\
             ┌────────┬─────────┐\n\
             │banana  │  peel   │\n\
             ├────────┼─────────┤\n\
             │oops    │  slipped│\n\
             └────────┴─────────┘\n\
             Do you see a neat rectangle composed of four rectangles with words inside them?  \
             Are the vertical lines contiguous?  \
             If not, something is wrong with your terminal!  You may need to change the terminal \
             font.  We use Menlo Regular at 13pt.  However, you may still observe small vertical\n\
             gaps between characters in some instances, depending on your computer and \
             apparently\nresulting from bugs in font rendition.\n\n\
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
             good to go!\n\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide quick help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "quick") || h.help_all {
        h.begin_doc("quick")?;
        h.print("\n")?;
        h.print("\\bold{quick guide to getting started}\n\n")?;
        h.print(
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
             Please read on to learn more!\n\n\
             \\bold{navigation in enclone}\n\n\
             enclone automatically sends its output through the program \"less\".  This allows you \
             to navigate within the output, using the following keys \
             (and many more, not shown, and which you don't need to know):\n\
             • space: causes output to page forward\n\
             • b: causes output to page backward\n\
             • /string: finds instances of \"string\" in the output\n\
             • n: having done the previous, jump to the next instance\n\
             • q: quit, to return to the command line.\n\n\
             When enclone uses less, it passes the argument -R, which causes certain characters \
             to be hidden, namely escape codes that color or bold text.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide how help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "how") || h.help_all {
        // Start.

        h.begin_doc("how")?;
        h.print("\n")?;
        h.print("\\bold{information about how enclone works}\n\n")?;
        h.print(
            "The goal of enclone is to find and display the clonotypes within single cell \
             VDJ datasets: groups of cells having the same fully rearranged common ancestor.\n\n\
             \
             enclone provides the foundation for fully understanding each cell's antigen \
             affinity and the evolutionary relationship between cells within one or more datasets.  \
             This starts with, for each cell, \
             \\bold{the full length sequence of all its VDJ receptor chains}.  Such data may be \
             obtained using the 10x Genomics immune profiling platform.\n\n\
             See also the heuristics page at \\green{bit.ly/enclone}.\n\n\
             \
             For this, there are fundamental challenges:\n\n",
        )?;

        // Print challenges.

        h.print_with_box(
            "1. It is extremely easy to get false positives: the incorrect \
             appearance that two cells have a common ancestor.\n\n\
             \
             2. Because of somatic hypermutation in B cells, it can be difficult to know that \
             two B cells share a common ancestor.\n\n\
             \
             3. There is always some background noise, e.g. from ambient mRNA.  When building \
             large clonotypes, this noise tends to pile up, yielding ectopic chains, i.e. chains \
             within a clonotype that are artifacts and do not represent true biology.",
            false,
        )?;

        // Print boxed algorithm.

        h.print(
            "To address these challenges, the enclone algorithm has several steps, which we \
             outline:\n\n",
        )?;
        h.print(
            "\\boldred{1}.  Input data.  \
             enclone gets its information from the file all_contig_annotations.json that is \
             produced by Cell Ranger.  Only productive contigs are used.  Each has an annotated \
             V and J segment.  The V segment alignment may have a single indel whose length is \
             divisible by three, and in that case, the V reference sequence is edited either to \
             delete or insert sequence.  In the insertion case, the bases are taken from the \
             contig.  These indels are noted in the enclone output.\n\n\
             \
             \\boldred{2}.  Exact subclonotypes.  \
             enclone groups cells into exact subclonotypes, provided that they have the same \
             number of chains, identical V..J sequences, identical C segment assignments, \
             and the same distance between the J stop and the C start (which is usually zero).\n\n\
             \
             \\boldred{3}.  Finding the germline sequences.  For datasets from a given donor, \
             enclone derives \"donor reference sequences\" for the V chains present in the donor's \
             genome.  This is powerful, even though based on imperfect information.  V segments \
             vary in their expression frequency and thus the more cells which are present, the \
             more complete the information will be.  It is also not possible to accurately \
             determine the terminal bases in a V chain from transcript data alone because these \
             bases mutate during recombination and because of non-templated nucleotide addition.\n\n\
             \
             The idea for how this is done is roughly the following: for each V segment, we choose \
             one cell from each clonotype (although these have not actually been computed yet, so \
             it's an approximation).  Next for each position on the V segment, excluding the last \
             15 bases, we determine the distribution of bases that occur within these selected \
             cells.  We only consider those positions where a non-reference base occurs at least \
             four times and is at least 25% of the total.  Then each cell has a footprint relative \
             to these positions; we require that these footprints satisfy similar evidence \
             criteria.  Each such non-reference footprint then defines an \"alternate allele\".  \
             We do not restrict the number of alternate alleles because they may arise from \
             duplicated gene copies.\n\n\
             \
             A similar approach was attempted for J segments but at the time of testing did not \
             appear to enhance clonotyping specificity.  This could be revisited later and might \
             be of interest even if it does not improve specificity.\n\n\
             \
             \\boldred{4}.  What joins are tested.  \
             Pairs of exact subclonotypes are considered for joining, as described below.  This \
             process only considers exact subclonotypes have two or three chains.  There is some \
             separate joining for the case of one chain.  Exact subclonotypes having four chains \
             are not joined at present.  These cases are clearly harder because these exact \
             subclonotypes are highly enriched for cell doublets, which we discard if we can \
             identify as such.\n\n\
             \
             \\boldred{5}.  Initial grouping.  \
             For each pair of exact subclonotypes, and for each pair of chains in each of the \
             two exact subclonotypes, for which V..J has the same length for the corresponding \
             chains, and the CDR3 segments have the same length for the corresponding chains, \
             enclone considers joining the exact subclonotypes into the same clonotype.\n\n\
             \
             \\boldred{6}.  Error bounding.  \
             To proceed, as a minimum requirement, there must be at most \\bold{55} total \
             mismatches between the two exact subclonotypes, within the given two V..J segments.\n\
             This can be changed by setting \\bold{MAX_DIFFS=n} on the command line.  (Note
             that for CellRanger version 5.0, the value is instead \\bold{50}.)\n\n\
             \
             \\boldred{7}.  Shared mutations.  \
             enclone next finds shared mutations betweens exact subclonotypes, that is, for \
             two exact subclonotypes, common mutations from the reference sequence, using the \
             donor reference for the V segments and the universal reference for the J segments.  \
             Shared mutations are supposed to be somatic hypermutations, that would be evidence \
             of common ancestry.  By using the donor reference sequences, most shared germline \
             mutations are excluded, and this is critical for the algorithm's success.\n\n\
             \
             \\boldred{8}.  Are there enough shared mutations?  \
             We find the probability p that “the shared mutations occur by chance”.  More \
             specifically, given d shared mutations, and k total mutations (across the two cells), \
             we compute the probability p that a sample with replacement of k items from a set \
             whose size is the total number of bases in the V..J segments, yields at most k – d \
             distinct elements.  The probability is an approximation, for the method please see\n\
             \\green{https://docs.rs/stirling_numbers/0.1.0/stirling_numbers}.\n\n\
             \
             \\boldred{9}.  Are there too many CDR3 mutations?  \
             Next, let N be \"the number of DNA sequences that differ from the given CDR3 \
             sequences by at most the number of observed differences\".  More specifically, if \
             cd is the number of differences between the given CDR3 nucleotide sequences, and n \
             is the total length in nucleotides of the CDR3 sequences (for the two chains), we \
             compute the total number N of strings of length n that are obtainable by perturbing \
             a given string of length n, which is\nsum( choose(n,m), m = 0..=cd) ).  We also \
             require that cd is at most 15 (and this bound is adjustable via the command-line \
             argument \\bold{MAX_CDR3_DIFFS}).  (The value for Cell Ranger 5.0 is 10.)\n\n\
             \
             \\boldred{10}.  Key join criteria.  \
             Two cells sharing sufficiently many shared differences and sufficiently few \
             CDR3 differences are deemed to be in the same clonotype.  That is, The lower p is, \
             and the lower N is, the more likely it is that the shared mutations represent bona \
             fide shared ancestry.  Accordingly, the smaller p*N is, the more likely it is that \
             two cells lie in the same true clonotype.  To join two cells into the same \
             clonotype, we require that the bound p*n ≤ C is satisfied, where C is the \
             constant 500,000.  (The value for Cell Ranger 5.0 is 1,000,000.) \
             The value may be adjusted using the command-line argument \\bold{MAX_SCORE}, or the \
             log10 of this,\n\
             \\bold{MAX_LOG_SCORE}.  This constant was arrived at by empirically balancing \
             sensitivity and specificity across a large collection of datasets.  See results \
             described at \\green{bit.ly/enclone}.\n\n\
             \
             \\boldred{11}.  Other join criteria.  We do not join two clonotypes which were \
             assigned different reference sequences unless those reference sequences differ by \
             at most \\bold{2} positions.  This value can be controlled using the \
             command-line argument \\bold{MAX_DEGRADATION}.  (Note that for Cell Ranger 5.0, \
             the value is instead \\bold{3}.)  There is an additional restriction \
             imposed when creating two-cell clonotypes: we require that that \
             cd ≤ d, where cd is the number of CDR3 differences and d is the number of shared \
             mutations, as above.  This filter may be turned off \
             using the command-line argument \\bold{EASY}.\n\n\
             \
             \\boldred{12}.  Junk.  \
             Spurious chains are filtered out based on frequency and connections. See \
             \"enclone help special\" for a description of the filters.\n\n",
        )?;
        h.print(
            "We are actively working to improve the algorithm.  Test results for the current \
             version may be found at \\green{bit.ly/enclone}.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide command line help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "command") || h.help_all {
        h.begin_doc("command")?;
        h.print("\n")?;
        h.print("\\bold{information about enclone command-line argument processing}\n\n")?;
        h.print("\\bold{1. Order of processing}\n\n")?;
        h.print(
            "• Before processing its command line, enclone first checks for environment\n\
             variables of the form \\bold{ENCLONE_<x>}.  These are converted into command-line \
             arguments.  You can set any command-line argument this way.  The reason why you might \
             want to use this feature is if you find yourself using the same \
             command-line option over and over, and it is more convenient to set it once as \
             an environment variable.\n\
             • For example, setting the environment variable \\bold{ENCLONE_PRE} to \
             \\bold{/Users/me/enclone_data} \
             is equivalent to providing the command-line argument \
             \\bold{PRE=/Users/me/enclone_data}.\n\
             • After checking environment variables, arguments on the command line are read from \
             left to right; if an argument name is repeated, only the \
             rightmost value is used, except as noted specifically in the documentation.\n\n",
        )?;
        h.print("\\bold{2. Importing arguments}\n\n")?;
        h.print(
            "Extra arguments can be imported on the command line using \\bold{SOURCE=filename}.  \
            The file may have newlines, and more than one SOURCE command may be used.  Any \
            line starting with # is treated as a comment.\n\n",
        )?;
        h.print("\\bold{3. Color}\n\n")?;
        h.print_enclone()?;
        h.print(
            " uses ANSI escape codes for color and bolding, frivolously, for emphasis, \
             and more\nimportantly for amino acids, to represent different codons.  This is \
             done automatically but you can turn it off....\n\n\
             \\boldred{PLEASE READ THIS:}\n\n\
             You can turn off escape codes by adding \\bold{PLAIN} to any command.  Use this if \
             you want to peruse output using a text editor which does not grok the escape \
             codes.  However some things will not make sense without color.\n\n",
        )?;
        h.print("\\bold{4. Paging}\n\n")?;
        h.print("• enclone automatically pipes its output to \\bold{less -R -F -X}.\n")?;
        h.print(
            "• The effect of this will be that you'll see only the first screen of output.  \
             You can then use the spacebar to go forward, b to go backward, and q to quit.  \
             The \\bold{-R} option causes escape characters to be correctly displayed, the \
             \\bold{-F} option causes an automatic exit if output fits on a single screen, and \
             the \\bold{-X} option prevents output from being sent to the \"alternate screen\" \
             under certain platform/version combinations.\n",
        )?;
        h.print("• Type \\bold{man less} if you need more information.\n")?;
        h.print(
            "• If for whatever reason you need to turn off output paging, add the argument \
             \\bold{NOPAGER} to the enclone command.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide glossary help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "glossary") || h.help_all {
        h.begin_doc("glossary")?;
        h.print("\n")?;

        // intro

        h.print("\\bold{glossary of terms used by enclone}\n\n")?;

        // doc V..J

        h.doc(
            "V..J",
            "the full sequence of a V(D)J transcript, from the beginning of the V",
        );
        h.doc2("segment to the end of the J segment; this sequence begins with a stop codon");
        h.doc2("and ends with a partial codon (its first base)");

        // doc CDR3

        h.docf2(
            "CDR3",
            "The terms CDR3 and junction are commonly mistaken and often used \
            interchangeably.  In enclone's nomenclature, \"CDR3\" actually refers to the \
            junction (the CDR3 loop plus the canonical C and W/F at the N and C termini \
            respectively).",
            60,
        )?;

        // doc clonotype

        h.ldoc(
            "clonotype",
            "all the cells descended from a single fully rearranged T or B cell",
        );
        h.doc("", "(approximated computationally)");

        // doc exact subclonotype

        h.docpr(
            "exact subclonotype",
            "all cells having identical transcripts \\boldred{○}",
        );
        h.doc2("(every clonotype is a union of exact subclonotypes)");

        // doc clone

        h.doc(
            "clone",
            "a cell in a clonotype, or in an exact subclonotype",
        );

        // doc onesie etc.

        h.ldoc(
            "onesie",
            "a clonotype or exact subclonotype having exactly one chain",
        );
        h.doc(
            "twosie",
            "a clonotype or exact subclonotype having exactly two chains",
        );
        h.doc(
            "threesie",
            "a clonotype or exact subclonotype having exactly three chains;",
        );
        h.doc(
            "",
            "these frequently represent true biological events, arising from expression",
        );
        h.doc2("of both alleles");
        h.doc(
            "foursie",
            "a clonotype or exact subclonotype having exactly four chains;",
        );
        h.doc("", "these very rarely represent true biological events");
        h.doc("moresie", "a clonotype having more than four chains;");
        h.doc2("these sad clonotypes do not represent true biological events");

        // doc donor etc.

        h.ldoc(
            "donor",
            "an individual from whom datasets of an origin are obtained",
        );
        h.doc(
            "origin",
            "a tube of cells from a donor, from a particular tissue at a",
        );
        h.doc2("particular point in time, and possibly enriched for particular cells");
        h.doc(
            "cell group",
            "an aliquot from an origin, presumed to be a random draw",
        );
        h.doc(
            "dataset",
            "all sequencing data obtained from a particular library type",
        );
        h.doc2("(e.g. TCR or BCR or GEX or FB), from one cell group, processed by running");
        h.doc2("through the Cell Ranger pipeline");

        // print main table

        h.print_tab2()?;
        h.print("\n")?;

        // print footnote

        h.print(
            "\\boldred{○} The exact requirements for being in the same exact subclonotype are \
             that cells:\n\
             • have the same number of productive contigs identified\n\
             • that these have identical bases within V..J\n\
             • that they are assigned the same constant region reference sequences\n\
             • and that the difference between the V stop and the C start is the same\n  \
             (noting that this difference is nearly always zero).\n\
             Note that we allow mutations within the 5'-UTR and constant regions.\n\n",
        )?;

        // conventions

        h.print("\\bold{conventions}\n\n")?;
        h.print(
            "• When we refer to \"V segments\", we always include the leader segment.\n\
             • Zero or one?  We number exact subclonotypes as 1, 2, ... and likewise with\n\
             chains within a clonotype, however DNA and amino-acid positions are numbered starting \
             at zero.\n\n",
        )?;

        // done

        h.end_doc();
    }
    Ok(())
}
