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
    // Set up.

    /*
    let mut args = args.clone();
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
    let mut h.help_all = false;
    unsafe {
        if HELP_ALL {
            h.help_all = true;
        }
    }
    */

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
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
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
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
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
            pal.push('▓');
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
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
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

        h.print("\\boldblue{2. How can I print the entire enclone documentation?}\n\n");
        h.print(
            "We don't know how in general, but the following works for us from a mac:\n\n\
             A. open a new terminal window\n\
             B. make it 111 characters wide; font should be fixed width and roughly 12pt\n\
             C. type \"enclone help all NOPAGER\"\n\
             D. type command-A to select all\n\
             E. type option-command-P to print selection\n\
             F. click the PDF button in the lower left (drop down menu)\n\
             G. click \"Open in Preview\"\n\
             H. then print (or save the pdf, if you prefer).\n\n",
        );

        h.print("\\boldblue{3. Why is enclone slow for me?}\n\n");
        h.print(
            "On a single VDJ dataset, it typically runs for us in a few seconds, on a Mac or Linux \
             server.  Runs where we combine several hundred datasets execute in a couple minutes \
             (on a server).  Your mileage could vary, and we are interested in cases where \
             it is underperforming.  Let us know.  We are aware of several things that could be \
             done to speed up enclone.\n\n",
        );

        h.print(
            "\\boldblue{4. How does enclone fit into the 10x Genomics software ecosystem?}\n\n",
        );
        h.print(
            "There are several parts to the answer:\n\
             • enclone is a standalone executable that by default produces human-readable output.\n\
             • You can also run enclone to produce parseable output \
             (see \\bold{enclone help parseable}), \
             and that output can be digested using code that you write (for example, in R).\n\
             • When you run Cell Ranger to process 10x single cell immune profiling data, it in \
             effect calls enclone with a special option that yields only an output file for \
             the 10x visualization tool Loupe.\n\
             • Clonotypes may then be viewed using Loupe.  The view of a clonotype provided by \
             Loupe is different than the view provided by enclone.  Loupe shows a continuous \
             expanse of bases across each chain, which you can scroll across, rather than the \
             compressed view of \"important\" bases or amino acids that enclone shows.\n\n",
        );

        h.print("\\boldblue{5. What platforms does enclone run on?}\n\n");
        h.print(
            "1. linux/x86-64 (that's most servers)\n\
             2. mac.\n\n\
             However, we have not and cannot test every possible configuration of these \
             platforms.  Please let us know if you encounter problems!\n\n",
        );

        h.print("\\boldblue{6. How can I print out all the donor reference sequences?}\n\n");
        h.print(
            "Add the argument \\bold{DONOR_REF_FILE=filename} to your enclone command, \
             and fasta for the donor reference sequences will be dumped there.\n\n",
        );

        h.print("\\boldblue{7. How does enclone know what VDJ reference sequences I'm using?}\n\n");
        h.print(
            "It does not!  It assumes that you have the \\bold{human} reference sequences that \
             shipped with the latest version of Cell Ranger.  If instead your sequences are mouse, \
             then you can specify that by adding the argument \\bold{MOUSE} to your command \
             line.  If you are simply using another reference sequence, please specify that \
             by adding \\bold{REF=f} to your command line, where \\bold{f} is the name of your \
             VDJ reference fasta file.\n\n",
        );

        h.print("\\boldblue{8. Can I provide data from more than one donor?}\n\n");
        h.print(
            "Yes.  Type \\bold{enclone help input} for details.  The default behavior of \
             enclone is to prevent cells from different donors from being placed in the same \
             clonotype.  The \\bold{MIX_DONORS} option may be used to turn off this behavior.  If \
             you employ this option, then clonotypes containing cells from more than one donor \
             will be flagged as errors, unless you use the \\bold{NWARN} option to turn off those \
             warnings.  The primary reason for allowing entry of data from multiple \
             donors is to allow estimation of enclone's error rate.\n\n",
        );

        h.print("\\boldblue{9. What are some command line argument values quoted?}\n\n");
        h.print(
            "Command line argument values that contain any of these characters ;|* need to \
             be quoted like so\n\
             \\bold{TCR=\"a;b\"}\n\
             to prevent the shell from interpreting them for a purpose completely unrelated to \
             enclone.  This is a trap, because forgetting to add the quotes can result in \
             nonsensical and confusing behavior!\n\n",
        );

        h.print("\\boldblue{10. If enclone fails, does it return nonzero exit status?}\n\n");
        h.print(
            "Yes, unless output of enclone is going to a terminal.  In that case, you'll always \
             get zero.\n\n",
        );

        h.print("\\boldblue{11. Could a cell be missing from an enclone clonotype?}\n\n");
        h.print(
            "Yes, some cells are deliberately deleted.  The cell might have been deleted by \
             one of the filters described in \\bold{enclone help special}, and which you can \
             turn off.  We also delete cells for which more than four chains were found.\n\n",
        );

        h.print("\\boldblue{12. Can enclone print summary stats?}\n\n");
        h.print(
            "Yes, if you add the option \\bold{SUMMARY}, then some summary stats will be \
             printed.  If you only want to see the summary stats, then also add the option \
             \\bold{NOPRINT}.\n\n",
        );

        h.print("\\boldblue{13. What is the notes column?}\n\n");
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
             If you supply enclone with output from Cell Ranger 3.1, \
             and J and C segments overlap by exactly one, this will not be noted.  \
             The reason for this is that many of the reference sequences supplied with \
             Cell Ranger 3.1 had an extra base at the beginning of their C segments, resulting \
             in annoying overlap notes for a large fraction of clonotypes.\n\n",
        );

        h.print("\\boldblue{14. Can I cap the number of threads used by enclone?}\n\n");
        h.print(
            "You can use the command-line argument \\bold{MAX_CORES=n} to cap the number of \
             cores used in parallel loops.  The number of threads used is typically one \
             higher.\n\n",
        );

        h.print("\\boldblue{15. Does enclone work under Windows?}\n\n");
        h.print(
            "No.  There are nontrivial technical problems with getting this to work.  If you're \
             sufficiently curious, see the notes in the source code file misc1.rs.  Please let us \
             know if you're interested in support for Windows.\n\n",
        );

        h.print("\\boldblue{16. Can I use enclone if I have only gene expression data?}\n\n");
        h.print(
            "Possibly.  In some cases this works very well, but in other cases it does not.  \
            Success depends on dataset characteristics that have not been carefully investigated.  \
            To attempt this, you need to invoke Cell Ranger on the GEX dataset as if \
            it was a VDJ dataset, and you need to specify to Cell Ranger that the run is to be \
            treated as BCR or TCR.  Two separate invocations can be used to get both.  Note also \
            that Cell Ranger has been only minimally tested for this configuration and that this \
            is not an officially supported Cell Ranger configuration.\n\n",
        );

        h.print("\\boldblue{17. How can I cite enclone?}\n\n");
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
                "10x Genomics, https://github.com/10XGenomics/enclone, version {}.\n",
                VERSION_STRING.before(",")
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

        h.print("\\boldblue{18. Can enclone output html?}\n\n");
        h.print(
            "Yes, just add the argument \\bold{HTML} to the command line.  Currently this does \
            not work with help.\n\n",
        );

        h.dump();
        std::process::exit(0);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide developer help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "developer") || h.help_all {
        h.begin_doc("developer");
        h.print("\n\\bold{a few options for developers}\n\n");
        h.print(
            "For instructions on how to compile, please see\n\
             \\green{https://github.com/10XDev/enclone/blob/master/COMPILE.md}.\n\n",
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
        if !h.help_all {
            h.dump();
            std::process::exit(0);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Finish enclone help all.

    if h.help_all {
        h.dump();
    }

    // Catch unrecognized help requests.

    if args.len() >= 2 {
        let mut x = args[1].clone();
        x.make_ascii_lowercase();
        if x.contains("help") {
            println!("\nYour help request doesn't match one known to enclone.\n");
            println!("Please type \"enclone\" to see the help options.\n");
            std::process::exit(1);
        }
    }
}
