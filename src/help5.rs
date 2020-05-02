// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::defs::*;
use crate::help_utils::*;
use ansi_escape::*;
use string_utils::*;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help5(args: &Vec<String>, ctl: &EncloneControl, h: &mut HelpDesk) {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide indels help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "indels") || h.help_all {
        h.begin_doc("indels");
        h.print("\n\\bold{handling of insertions and deletions}\n\n");
        h.print(
            "enclone can recognize and display a single insertion or deletion in a contig \
             relative to the reference, so long as its length is divisible by three, is relatively \
             short, and occurs within the V segment, not too close to its right end.\n\n\
             These indels could be germline, however most such events are already captured in a \
             reference sequence.  Currently the donor reference code in enclone does not recognize \
             indels.\n\n\
             SHM deletions are rare, and SHM insertions are even more rare.\n\n\
             Deletions are displayed using hyphens (-).  If you use the \\bold{var} option for \
             \\bold{cvars}, the hyphens will be displayed in base space, where they are initially \
             observed.  For the \\bold{AMINO} option, the deletion is first shifted by up to two \
             bases, so that the deletion starts at a base position that is divisible by three.  \
             Then the deleted amino acids are shown as hyphens.\n\n\
             Insertions are shown only in amino acid space, in a special per-chain column called \
             \\bold{notes} that \
             appears if there is an insertion.  Colored amino acids are shown for the insertion, \
             and the position of the insertion is shown.  The position is the position of the \
             amino acid after which the insertion appears, where the first amino acid (start \
             codon) is numbered 0.\n\n",
        );
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide ideas help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "ideas") || h.help_all {
        h.begin_doc("ideas");
        h.print("\n\\bold{features that might be implemented in enclone}\n\n");
        h.doc("speed", "make enclone faster");
        h.ldoc(
            "CDRn",
            "make CDR1 and CDR2 viewable in the same way that CDR3 is now",
        );
        h.ldoc(
            "distance grouping",
            "provide an option to group clonotypes by distance",
        );
        h.ldoc("cloning", "package V..J region into a cloning vector");
        h.ldoc(
            "phylogeny",
            "generate a phylogeny for the exact clonotypes within a clonotype",
        );

        h.ldoc("windows", "make enclone work on windows computers");
        h.print_tab2();
        h.print(
            "\nPlease let us know if you are interested in these features, or if there are \
             other features that you would like us to implement!\n\n",
        );
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide color help.
    //
    // Here, and in substitute_enclone_color in plot.rs, we swap the order of colors, placing the
    // last three before the first three.  This is because the last three seem to make a better
    // three-color palette.

    if (args.len() == 3 && args[1] == "help" && args[2] == "color") || h.help_all {
        h.begin_doc("color");
        h.print("\nHere is the color palette that enclone uses for amino acids:\n\n");
        let mut pal = String::new();
        for i in 0..6 {
            let s = best_color_order(i);
            let mut log = Vec::<u8>::new();
            if !h.plain {
                print_color(s, &mut log);
                pal += &stringme(&log);
            }
            pal.push('█');
            let mut log = Vec::<u8>::new();
            if !h.plain {
                emit_end_escape(&mut log);
                pal += &stringme(&log);
            }
            if i < 6 {
                pal.push(' ');
            }
        }
        h.print_plain(&format!("{}\n", pal));
        h.print(
            "\nWhen enclone shows amino acids, it colors each codon differently, via \
             the following scheme:\n\n",
        );
        h.print_plain(&format!("{}\n\n", colored_codon_table(h.plain)));
        h.print(
            "Colored amino acids enable the compact display of all the information in a \
             clonotype.\n\n",
        );
        h.print(
            "The coloring is done using special characters, called ANSI escape characters.  \
             Color is used occasionally elsewhere by enclone, and there is also some  \
             bolding, accomplished using the same mechanism.\n\n\
             Correct display of colors and bolding depends on having a terminal window \
             that is properly set up.  As far as we know, this may always be the case, \
             but it is possible that there are exceptions.  In addition, in general, text \
             editors do not correctly interpret escape characters.\n\n\
             For both of these reasons, you may wish to turn off the \"special effects\", \
             either some or all of the time.  You can do this by adding the argument\n",
        );
        h.print("\\bold{PLAIN}\n");
        h.print("to any enclone command.\n\n");
        h.print(
            "We know of two methods to get enclone output into another document, along \
             with colors:\n\
             1. Take a screenshot.\n\
             2. Open a new terminal window, type the enclone command, and then convert the \
             terminal window into a pdf.  See \\bold{enclone help faq} for related \
             instructions.\n\n",
        );
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide faq help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "faq") || h.help_all {
        h.begin_doc("faq");
        h.print("\n\\boldred{Frequently Asked Questions}\n\n");
        h.print(
            "We're sorry you're having difficulty!  Please see the answers below, check out \
             the other help guides, and if you're still stuck, write to us at \
             enclone@10xgenomics.com.\n\n",
        );

        h.print("\\boldblue{1. Why is my enclone output garbled?}\n\n");
        h.print(
            "We can think of two possibilities:\n\n\
            A. The escape characters that enclone emits for color and bolding are not getting\n\
            translated.  You have some options:\n\
            (a) Turn off escape character generation by adding PLAIN to your enclone commands.\n\
            This will work but you'll lose some information.\n\
            (b) If your terminal window is not translating escape characters, ask someone\n\
            with appropriate expertise to help you.  We have not observed this phenomenon,\n\
            but it should be fixable.\n\
            (c) If you're trying to view enclone output, with escape characters, using an editor,\n\
            that's probably not going to work well.\n\n\
            B. Perhaps enclone is emitting very wide lines.  Here are things you can do about \
            this:\n\
            (a) Make your terminal window wider or reduce the font size.\n\
            (b) Identify the field that is very wide and use the column controls to remove that\n\
            field.  See the help for lvars and cvars.  For example,\n",
        );
        h.print(
            "\\bold{AMINO=cdr3}\n\
             may help, or even\n\
             \\bold{AMINO=}\n\n",
        );

        h.print("\\boldblue{2. Can I convert the enclone visual output into other forms?}\n\n");
        h.print(
            "Yes, there are choices:\n\
             \\bold{A}. On a Mac, you can screenshot from a terminal window.\n\
             \\bold{B}. Add the argument \\bold{HTML} to the enclone command line.  Then the \
             output will be presented as html, with title \"enclone output\".  If you want to \
             set the title, use \\bold{HTML=\"...\"}.\n\
             \\bold{C}. You can then convert the html to pdf.  The best way on a Mac is to open \
             Safari, which is the best browser for this particular purpose, \
             select the file where you've saved the html, and then export as pdf.  Do not convert \
             to pdf via printing, which produces a less readable file, and also distorts colors.  \
             (We do not know why the colors are distorted.)\n\
             \\bold{D}. If you want to put enclone output in a Google Doc, you can do it via \
             approach \\bold{A}, although then you won't be able to select text \
             within the copied region.  \
             Alternatively, if you open the html file in a browser, you can then select \
             text (including clonotype box text) and paste into a Google Doc.  It will be pretty \
             ugly, but will capture color and correctly render the box structure, provided that \
             you use an appropriate fixed-width font for that part of the Doc.  We found that \
             Courier New works, with line spacing set to 0.88.  You may have to reduce the font \
             size.\n\n",
        );

        h.print("\\boldblue{3. Why is enclone slow for me?}\n\n");
        h.print(
            "On a single VDJ dataset, it typically runs for us in a few seconds, on a Mac or Linux \
             server.  Runs where we combine several hundred datasets execute in a couple minutes \
             (on a server).  Your mileage could vary, and we are interested in cases where \
             it is underperforming.  Let us know.  We are aware of several things that could be \
             done to speed up enclone.\n\n",
        );

        h.print("\\boldblue{4. What platforms does enclone run on?}\n\n");
        h.print(
            "1. Linux/x86-64 (that's most servers)\n\
               2. Mac.\n\n\
             However, we have not and cannot test every possible configuration of these \
             platforms.  Please let us know if you encounter problems!\n\n",
        );

        h.print("\\boldblue{5. How can I print out all the donor reference sequences?}\n\n");
        h.print(
            "Add the argument \\bold{DONOR_REF_FILE=filename} to your enclone command, \
             and fasta for the donor reference sequences will be dumped there.\n\n",
        );

        h.print("\\boldblue{6. How does enclone know what VDJ reference sequences I'm using?}\n\n");
        h.print(
            "It does not!  It assumes that you have the \\bold{human} reference sequences that \
             shipped with the latest version of Cell Ranger.  If instead your sequences are mouse, \
             then you can specify that by adding the argument \\bold{MOUSE} to your command \
             line.  If you are simply using another reference sequence, please specify that \
             by adding \\bold{REF=f} to your command line, where \\bold{f} is the name of your \
             VDJ reference fasta file.\n\n",
        );

        h.print("\\boldblue{7. Can I provide data from more than one donor?}\n\n");
        h.print(
            "Yes.  Type \\bold{enclone help input} for details.  The default behavior of \
             enclone is to prevent cells from different donors from being placed in the same \
             clonotype.  The \\bold{MIX_DONORS} option may be used to turn off this behavior.  If \
             you employ this option, then clonotypes containing cells from more than one donor \
             will be flagged as errors, unless you use the \\bold{NWARN} option to turn off those \
             warnings.  The primary reason for allowing entry of data from multiple \
             donors is to allow estimation of enclone's error rate.\n\n",
        );

        h.print("\\boldblue{8. What are some command line argument values quoted?}\n\n");
        h.print(
            "Command line argument values that contain any of these characters ;|* need to \
             be quoted like so\n\
             \\bold{TCR=\"a;b\"}\n\
             to prevent the shell from interpreting them for a purpose completely unrelated to \
             enclone.  This is a trap, because forgetting to add the quotes can result in \
             nonsensical and confusing behavior!\n\n",
        );

        h.print("\\boldblue{9. If enclone fails, does it return nonzero exit status?}\n\n");
        h.print(
            "Yes, unless output of enclone is going to a terminal.  In that case, you'll always \
             get zero.\n\n",
        );

        h.print("\\boldblue{10. Could a cell be missing from an enclone clonotype?}\n\n");
        h.print(
            "Yes, some cells are deliberately deleted.  The cell might have been deleted by \
             one of the filters described in \\bold{enclone help special}, and which you can \
             turn off.  We also delete cells for which more than four chains were found.\n\n",
        );

        h.print("\\boldblue{11. Can enclone print summary stats?}\n\n");
        h.print(
            "Yes, if you add the option \\bold{SUMMARY}, then some summary stats will be \
             printed.  If you only want to see the summary stats, then also add the option \
             \\bold{NOPRINT}.\n\n",
        );

        h.print("\\boldblue{12. What is the notes column?}\n\n");
        h.print(
            "The notes column appears if one of two relatively rare events occurs:\n\n\
             1. An insertion is detected in a chain sequence, relative to the reference.\n\n\
             2. The end of the J segment on a chain sequence does not exactly coincide with\n   \
             the beginning of the C segment.\n\
             The latter could correspond to one of several phenomena:\n\
             a. A transcript has an insertion between its J and C segments.\n   \
             This can happen.  See e.g. Behlke MA, Loh DY.\n   \
             Alternative splicing of murine T-cell receptor beta-chain transcripts.\n   \
             Nature 322(1986), 379-382.\n\
             b. There is an error in a reference sequence segment.\n   \
             We have tried to eliminate all such errors from the built-in references for\n   \
             human and mouse.\n\
             c. A cell produced a nonstandard transcript and also standard ones, and the\n   \
             Cell Ranger pipeline just happened to pick a nonstandard one.\n\
             d. There was a technical artifact and the sequence does not actually represent\n   \
             an mRNA molecule.\n\n\
             Overlaps of length exactly one between J and C segments are not shown unless you \
             specify the option \\bold{JC1}.  The reason for this is that certain reference \
             sequences (notably those from IMGT and those supplied with Cell Ranger 3.1) often \
             have an extra base at the beginning of their C segments, resulting \
             in annoying overlap notes for a large fraction of clonotypes.\n\n",
        );

        h.print("\\boldblue{13. Can I cap the number of threads used by enclone?}\n\n");
        h.print(
            "You can use the command-line argument \\bold{MAX_CORES=n} to cap the number of \
             cores used in parallel loops.  The number of threads used is typically one \
             higher.\n\n",
        );

        h.print("\\boldblue{14. Does enclone work under Windows?}\n\n");
        h.print(
            "No.  There are nontrivial technical problems with getting this to work.  If you're \
             sufficiently curious, see the notes in the source code file misc1.rs.  Please let us \
             know if you're interested in support for Windows.\n\n",
        );

        h.print("\\boldblue{15. Can I use enclone if I have only gene expression data?}\n\n");
        h.print(
            "Possibly.  In some cases this works very well, but in other cases it does not.  \
            Success depends on dataset characteristics that have not been carefully investigated.  \
            To attempt this, you need to invoke Cell Ranger on the GEX dataset as if \
            it was a VDJ dataset, and you need to specify to Cell Ranger that the run is to be \
            treated as BCR or TCR.  Two separate invocations can be used to get both.  Note also \
            that Cell Ranger has been only minimally tested for this configuration and that this \
            is not an officially supported Cell Ranger configuration.\n\n",
        );

        h.print("\\boldblue{16. How can I cite enclone?}\n\n");
        h.print("This version of enclone has been provided under a non-disclosure agreement,\n");
        h.print(
            "however once enclone has officially launched, you will be able to cite this \
             version as:\n",
        );
        let mut log = Vec::<u8>::new();
        emit_green_escape(&mut log);
        h.print(&format!("{}", strme(&log)));
        if !ctl.gen_opt.stable_doc {
            h.print(&format!(
                "10x Genomics, https://github.com/10XGenomics/enclone,\nversion {}.\n",
                VERSION_STRING
            ));
        } else {
            h.print(
                "10x Genomics, https://github.com/10XGenomics/enclone,\n\
                    (your enclone version information will be printed here).\n",
            );
        }
        let mut log = Vec::<u8>::new();
        emit_end_escape(&mut log);
        h.print(&format!("{}", strme(&log)));
        h.print(
            "At some point subsequent to that, there will be a white paper to which you can refer, \
            in addition to a DOI minted at Zenodo.  In the spirit of reproducibility, you should \
            provide the arguments that you used when you ran enclone and indicate the version of \
            Cell Ranger that you used to generate the input data.\n\n",
        );

        h.print("\\boldblue{17. Can I print the enclone version?}\n\n");
        h.print("Yes, type \"enclone version\".\n\n");

        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide developer help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "developer") || h.help_all {
        h.begin_doc("developer");
        h.print("\n\\bold{a few options for developers}\n\n");
        h.print(
            "For instructions on how to compile, please see\n\
             \\green{bit.ly/enclone}.\n\n",
        );
        h.doc(
            "COMP",
            "report computational performance stats; use this with NOPRINT if you",
        );
        h.doc(
            "",
            "only want to see the computational performance stats, and with NOPAGER if you",
        );
        h.doc("", "want output to be unbuffered");
        h.doc(
            "COMP2",
            "like COMP, but adds more detailed lines that are prefixed with --",
        );
        h.ldoc(
            "LONG_HELP",
            "allow long lines in help pages, which will otherwise trigger an assert",
        );
        h.ldoc(
            "CTRLC",
            "upon CTRL-C, emit a traceback and then exit; can be used as a primitive",
        );
        h.doc(
            "",
            "but easy profiling method, to know what the code is doing if it seems to be",
        );
        h.doc("", "very slow");
        h.ldoc(
            "HAPS=n",
            "Interrupt code n times, at one second intervals, get a traceback, and then tally",
        );
        h.doc(
            "",
            "the tracebacks.  This only works if the n tracebacks can be obtained before",
        );
        h.doc(
            "",
            "enclone terminates.  Interrupts that occur in the allocator are ignored, and",
        );
        h.doc(
            "",
            "in some cases, this accounts for most interrupts, resulting in confusing",
        );
        h.doc(
            "",
            "output.  In such cases, consider using CTRLC or a more sophisticated tool",
        );
        h.doc(
            "",
            "like perf.  Also HAPS only reports on the master thread, so to get useful",
        );
        h.doc(
            "",
            "information, you probably need to change an instance in the code of",
        );
        h.doc(
            "",
            "par_iter_mut to iter_mut, to turn off parallelization for a strategically",
        );
        h.doc("", "selected section.");
        h.print_tab2();
        h.print("\n");
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Finish enclone help all.

    if h.help_all {
        h.dump();
        std::process::exit(0);

    // Catch unrecognized help requests.
    } else if args.len() >= 2 {
        let mut x = args[1].clone();
        x.make_ascii_lowercase();
        if x.contains("help") {
            println!("\nYour help request doesn't match one known to enclone.\n");
            println!("Please type \"enclone\" to see the help options.\n");
            std::process::exit(1);
        }
    }
}
