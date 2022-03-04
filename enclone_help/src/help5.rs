// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::{colored_codon_table, HelpDesk};
use ansi_escape::{best_color_order, emit_end_escape, emit_green_escape, print_color};
use enclone_core::defs::EncloneControl;
use enclone_core::print_tools::color_by_property;
use enclone_core::version_string;
use io_utils::fwrite;
use std::io::Write;
use string_utils::{stringme, strme};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help5(args: &Vec<String>, ctl: &EncloneControl, h: &mut HelpDesk) -> Result<(), String> {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide display help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "display") || h.help_all {
        h.begin_doc("display")?;
        h.print("\n\\bold{other options that control clonotype display}\n\n")?;
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
        h.ldoc(
            "SUM",
            "print sum row for each clonotype (sum is across cells)",
        );
        h.doc(
            "MEAN",
            "print mean row for each clonotype (mean is across cells)",
        );

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "DIFF_STYLE=C1",
            "instead of showing an x for each amino acid column containing a difference, \
             show a C if the column lies within a complementarity-determining region, \
             and F if it lies in a framework region, and an L if it lies in the leader",
            75,
        )?;
        h.doc(
            "DIFF_STYLE=C2",
            "instead of showing an x for each amino acid column containing a difference,",
        );
        h.docpr(
            "",
            "show a \\boldred{◼} if the column lies within a complementarity-determining region,",
        );
        h.docpr("", "and otherwise show a \\bold{▮}.");
        h.ldoc(
            "CONX",
            "add an additional row to each clonotype table, showing the amino acid",
        );
        h.doc(
            "",
            "consensus across the clonotype, with X for each variant residue",
        );
        h.doc(
            "CONP",
            "add an additional row to each clonotype table, showing the amino acid",
        );
        h.doc(
            "",
            "consensus across the clonotype, with a property symbol whenever two different",
        );
        h.doc("", "amino acids are observed, see \"enclone help cvars\"");
        h.ldoc(
            "ALIGN<n>",
            "exhibit a visual alignment for chain n (for each exact subclonotype) to the",
        );
        h.doc2("donor V(D)J reference, picking the best D for heavy chains / TRB");
        h.doc2("Multiple values of n may be specified using multiple arguments.");
        h.doc(
            "ALIGN_2ND<n>",
            "same as ALIGN<n> but use second best D segment",
        );
        h.doc(
            "JALIGN<n>",
            "same as ALIGN<n> but only show the region from 15 bases before the end of the",
        );
        h.doc2("V segment to 35 bases into the J segment");
        h.doc(
            "JALIGN_2ND<n>",
            "same as JALIGN<n> but use second best D segment",
        );
        h.print_tab2()?;
        h.print("\n")?;

        h.print(
            "\\red{━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━}\n\n",
        )?;

        h.print(
            "\\bold{options that control clonotype grouping}\n\n\
             By default, enclone organizes clonotypes into groups, and each group contains \
             just one clonotype!  \
             We offer some options to do actual grouping, with the intention of reflecting \
             functional (antigen-binding) differences, but with many caveats because this is a \
             hard problem.\n\n\
             \
             These options are experimental.  There are many natural extensions that we have \
             not implemented.\n\n\
             \
             enclone has two types of grouping: symmetric and asymmetric.  Symmetric grouping \
             creates nonoverlapping groups, whereas asymmetric grouping creates groups that may \
             overlap.\n\n\
             \
             To turn on symmetric grouping, one uses a command of the form\n\
             \\bold{GROUP=c1,...,cn}\n\
             where each \\bold{ci} is a condition.  Two clonotypes are placed in the same group \
             if all the conditions are satisfied, and that grouping is extended transitively.\n\n\
             In what follows, \\bold{heavy chain} means IGH or TRB, and \\bold{light chain} means \
                IGK or IGL or TRA.  \n\n\
             Here are the conditions:\n\n",
        )?;
        h.rows.clear();

        h.doc(
            "vj_refname",
            "V segments have the same reference sequence name,",
        );
        h.doc2("and likewise for J segments");

        h.doc(
            "v_heavy_refname",
            "heavy chain V segments have the same reference sequence name",
        );
        h.doc(
            "vj_heavy_refname",
            "heavy chain V segments have the same reference sequence name,",
        );
        h.doc2("and likewise for J segments");
        h.doc2("(only applied to heavy chains)");
        h.doc(
            "vdj_refname",
            "V segments have the same reference sequence name,",
        );
        h.doc2("and likewise for D segments, computed from scratch, and J segments");
        h.doc(
            "vdj_heavy_refname",
            "V segments have the same reference sequence name,",
        );
        h.doc2("and likewise for D segments, computed from scratch, and J segments");
        h.doc2("(only applied to heavy chains)");

        h.ldoc(
            "len",
            "the lengths of V..J are the same (after correction for indels)",
        );
        h.doc("cdr3_len", "CDR3 sequences have the same length");
        h.doc(
            "cdr3_heavy_len",
            "heavy chain CDR3 sequences have the same length",
        );
        h.doc(
            "cdr3_light_len",
            "light chain CDR3 sequences have the same length",
        );
        h.ldoc(
            "cdr3_heavy≥n%",
            "nucleotide identity on heavy chain CDR3 sequences is at least n%",
        );
        h.doc(
            "cdr3_light≥n%",
            "nucleotide identity on light chain CDR3 sequences is at least n%",
        );
        h.doc(
            "cdr3_aa_heavy≥n%",
            "amino acid identity on heavy chain CDR3 sequences is at least n%",
        );
        h.doc(
            "cdr3_aa_light≥n%",
            "amino acid identity on light chain CDR3 sequences is at least n%",
        );
        h.doc2("(note that use of either of these options without at least one of the");
        h.doc2("earlier options may be slow)");
        h.doc2("(in both cases, we also recognize >= (with quoting) and ⩾)");
        h.ldoc(
            "heavy≥n%",
            "nucleotide identity on heavy chain V..J sequences is at least n%",
        );
        h.doc(
            "light≥n%",
            "nucleotide identity on light chain V..J sequences is at least n%",
        );

        h.doc(
            "aa_heavy≥n%",
            "amino acid identity on heavy chain V..J sequences is at least n%",
        );
        h.doc(
            "aa_light≥n%",
            "amino acid identity on light chain V..J sequences is at least n%",
        );
        h.doc2("(note that use of either of these options without at least one of the");
        h.doc2("earlier options may be very slow)");
        h.doc2("(in both cases, we also recognize >= (with quoting) and ⩾)");
        h.print_tab2()?;
        h.print("\n")?;

        h.print(
            "To instead turn on asymmetric grouping, one uses the \\bold{AGROUP} option.  To use \
            this, it is in addition necessary to define \"center clonotypes\", \
            a \"distance formula\", and a \"distance bound\".  Each group will then consist \
            of the center clonotype (which comes first), followed by, in order by distance \
            (relative to the formula), all those clonotypes that satisfy the distance bound \
            (with ties broken arbitrarily).  For each clonotype in a group, we print its \
            distance from the first clonotype, and this is also available as a parseable variable \
            \\bold{dist_center}.\n\n",
        )?;

        h.print(
            "\\bold{Center clonotypes.}  These are in principle any set of clonotypes.  \
            For now we allow two options:\n\
            \\bold{AG_CENTER=from_filters}\n\
            which causes all the filters described at \"enclone help filters\" to NOT filter \
            clonotypes in the usual way, but instead filter to define the center, and\n\
            \\bold{AG_CENTER=copy_filters}\n\
            which effectively does nothing -- it just says that filters apply to all clonotypes, \
            whether in the center or not.\n\n",
        )?;

        h.print_with_box(
            "Please note that asymmetric grouping is very time consuming, and run time is \
            roughly a linear function of (number of center clonotypes) * (number of clonotypes).  \
            So it is advisable to restrict the number of center clonotypes.",
            false,
        )?;

        h.print(
            "\\bold{Distance formula.}  This could in principle be any function that takes as \
            input two clonotypes and returns a number.  For now we allow only:\n\
            \\bold{AG_DIST_FORMULA=cdr3_edit_distance}\n\
            which is the \"Levenshtein CDR3 edit distance between two clonotypes\".  This is the \
            minimum, over all pairs of exact subclonotypes, one from each of the two clonotypes, \
            of the edit distance between two exact subclonotypes, which is the sum of the edit \
            distances between the heavy chains and between the light chains.\n\n",
        )?;

        h.print(
            "Technical note.  This is the explanation for the case where there are two chains \
            of different types.  Here is the explanation for the \"non-standard\" cases.  \
            We take the sum, over all \
            pairs of heavy chains, one from each of the two exact subclonotypes, of the \
            edit distance between the CDR3 sequences for the heavy chains, plus the same for \
            light chains.  Exact subclonotypes that lack a heavy or a light chain are ignored \
            by this computation.  Also the distance between two clonotypes is declared infinite \
            if one of them lacks a heavy chain or one of them lacks a light chain.\n\n",
        )?;

        h.print(
            "\\bold{Distance bound.}  For now we allow the following two forms:\n\
            \\bold{AG_DIST_BOUND=top=n}\n\
            which returns the top n clonotypes (plus the center), and\n\
            \\bold{AG_DIST_BOUND=max=d}\n\
            which returns all clonotypes having distance ≤ d from the center clonotype.\n\n",
        )?;

        h.print(
            "In addition, there are the following grouping options, for both the symmetric \
            and asymmetric cases:\n\n",
        )?;
        h.rows.clear();
        h.doc(
            "MIN_GROUP",
            "minimum number of clonotypes in group to print (default = 1)",
        );
        h.doc(
            "MIN_GROUP_DONORS",
            "minimum number of donors for a group to be printed (default = 1)",
        );
        h.doc(
            "GROUP_CDR3H_LEN_VAR",
            "only print groups having variable heavy chain CDR3 length",
        );
        h.doc(
            "GROUP_CDR3=x",
            "only print groups containing the CDR3 amino acid sequence x",
        );
        h.doc(
            "GROUP_DONOR=d",
            "only print groups containing a cell from the given donor;",
        );
        h.doc2("multiple instances may be used to jointly restrict");
        h.ldoc("NGROUP", "don't display grouping messages");
        h.print_tab2()?;
        h.print("\n")?;

        h.print(
            "\\red{━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━}\n\n",
        )?;

        h.print(
            "\\bold{options that display dataset variables}\n\n\
             enclone has some variables that are computed for each dataset, and whose values may \
             by printed as a table in the summary, and not otherwise used (currently).  \
             These may be specified using\n\
             \\bold{DVARS=var1,...,varn}.  The dataset-level variables that are supported \
             currently are:\n\
             <feature>_cellular_r\n\
             <feature>_cellular_u\n\
             which are, respectively, the percentage of reads [UMIs] for the given feature that \
             are in cells that were called by the cellranger pipeline.  A feature is e.g. \
             IGHG1_g etc. as discussed at \"enclone help lvars\".  To compute the metrics, the \
             cellranger output file per_feature_metrics.csv is read.  In addition, one may also \
             use numeric values defined in the file metrics_summary_json.json, but this file is in \
             general not available.  To get it, it may be necessary to rerun the cellranger \
             pipeline using --vdrmode=disable and then copy the json file to outs.  Finally, \
             variable names may be prefaced with abbreviation:, and in such cases, it is the \
             abbreviation that is displayed in the table.",
        )?;

        h.print(
            "\n\n\\red{━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━}\n\n",
        )?;

        h.print(
            "\\bold{options that control global variables}\n\n\
             enclone has some global variables that can be computed, with values printed in the \
             summary, and not otherwise used (currently).  These may be specified using \
             \\bold{GVARS=var1,...,varn}.  The global variables that are supported currently are:\n\
             d_inconsistent_%\n\
             d_inconsistent_n\n\
             Please see https://10xgenomics.github.io/enclone/pages/auto/d_genes.html for more \
             information.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide indels help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "indels") || h.help_all {
        h.begin_doc("indels")?;
        h.print("\n\\bold{handling of insertions and deletions}\n\n")?;
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
             and the position of the insertion is shown.  The notation e.g.\n\
             ins = TFT at 46\n\
             means that TFT is inserted after the first 46 amino acids.  Since the first amino \
             acid (which is a start codon) is numbered 0, the insertion actually occurs after \
             the amino acid numbered 45.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide color help.
    //
    // Here, and in substitute_enclone_color in plot.rs, we swap the order of colors, so that the
    // first three are as given, because they seem to make a better three-color palette.

    if (args.len() == 3 && args[1] == "help" && args[2] == "color") || h.help_all {
        h.begin_doc("color")?;
        h.print("\nHere is the color palette that enclone uses for amino acids:\n\n")?;
        let mut pal = String::new();
        for i in 0..7 {
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
            if i < 7 {
                pal.push(' ');
            }
        }
        h.print_plain(&format!("{}\n", pal))?;
        h.print(
            "\nWhen enclone shows amino acids, it uses one of three coloring schemes.  The first \
             scheme (the default, or using the argument \\bold{COLOR=codon}), colors amino \
             acids by codon, according to the following scheme:\n\n",
        )?;
        h.print_plain(&format!("{}\n\n", colored_codon_table(h.plain)))?;
        h.print(
            "Colored amino acids enable the compact display of all the information in a \
             clonotype.\n\n",
        )?;
        h.print(
            "The second scheme, \\bold{COLOR=codon-diffs}, is the same as the first, except that \
             some amino acids are \"grayed out\".  An amino acid is highlighted (not grayed out) \
             if (a) its codon differs from the universal reference or (b) it is in a CDR3 and the \
             codon is shared by half or less of the exact subclonotypes having the given \
             chain.  You may wish to use this with the \\bold{CONX} or \\bold{CONP} option.\n\n",
        )?;
        h.print(
            "The third scheme for coloring amino acids, \\bold{COLOR=property}, colors amino \
             acids by their properties, according to the following scheme:\n\n",
        )?;
        {
            let mut log = Vec::<u8>::new();
            if !h.plain {
                fwrite!(log, "1. Aliphatic: ");
                color_by_property(b"A G I L P V\n", &mut log);
                fwrite!(log, "2. Aromatic: ");
                color_by_property(b"F W Y\n", &mut log);
                fwrite!(log, "3. Acidic: ");
                color_by_property(b"D E\n", &mut log);
                fwrite!(log, "4. Basic: ");
                color_by_property(b"R H K\n", &mut log);
                fwrite!(log, "5. Hydroxylic: ");
                color_by_property(b"S T\n", &mut log);
                fwrite!(log, "6. Sulfurous: ");
                color_by_property(b"C M\n", &mut log);
                fwrite!(log, "7. Amidic: ");
                color_by_property(b"N Q\n", &mut log);
                h.print_plain(&format!("{}\n", stringme(&log)))?;
            } else {
                h.print(
                    "1. Aliphatic: A G I L P V\n\
                    2. Aromatic: F W Y\n\
                    3. Acidic: D E\n\
                    4. Basic: R H K\n\
                    5. Hydroxylic: S T\n\
                    6. Sulfurous: C M\n\
                    7. Amidic: N Q\n\n",
                )?;
            }
        }
        /*
        h.print(
            "The third scheme for coloring amino acids is for BCR, and has the form \
             \\bold{COLOR=peer.p} \
             where p is a decimal number (at most 100), representing a percentage.  For example, \
             one could have \\bold{COLOR=peer.1}, which would cause amino acids having \
             general frequency zero to be colored red, those having \
             general frequency ≤ 1% to be colored light blue, \
             and all others colored black.  For more information about this, please see\n\
             can just navigate to there from \\green{bit.ly/enclone}.)\n\n",
        )?;
        */
        h.print(
            "In all cases, \
             the coloring is done using special characters, called ANSI escape characters.  \
             Color is used occasionally elsewhere by enclone, and there is also some \
             bolding, accomplished using the same mechanism.\n\n\
             Correct display of colors and bolding depends on having a terminal window \
             that is properly set up.  As far as we know, this may always be the case, \
             but it is possible that there are exceptions.  In addition, in general, text \
             editors do not correctly interpret escape characters.\n\n\
             For both of these reasons, you may wish to turn off the \"special effects\", \
             either some or all of the time.  You can do this by adding the argument\n",
        )?;
        h.print("\\bold{PLAIN}\n")?;
        h.print("to any enclone command.\n\n")?;
        h.print(
            "We know of two methods to get enclone output into another document, along \
             with colors:\n\
             1. Take a screenshot.\n\
             2. Open a new terminal window, type the enclone command, and then convert the \
             terminal window into a pdf.  See \\bold{enclone help faq} for related \
             instructions.\n\n",
        )?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide faq help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "faq") || h.help_all {
        h.begin_doc("faq")?;
        h.print("\n\\boldred{Frequently Asked Questions}\n\n")?;
        h.print(
            "We're sorry you're having difficulty!  Please see the answers below, check out \
             the other help guides, and if you're still stuck, write to us at \
             enclone@10xgenomics.com.\n\n",
        )?;

        h.print("\\boldblue{1. Why is my enclone output garbled?}\n\n")?;
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
        )?;
        h.print(
            "\\bold{AMINO=cdr3}\n\
             may help, or even\n\
             \\bold{AMINO=}\n",
        )?;
        h.print("These options may also help: \\bold{CVARS=u FOLD_HEADERS}.\n\n")?;

        h.print("\\boldblue{2. Can I convert the enclone visual output into other forms?}\n\n")?;
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
        )?;

        h.print("\\boldblue{3. Why is enclone slow for me?}\n\n")?;
        h.print(
            "On a single VDJ dataset, it typically runs for us in a few seconds, on a Mac or Linux \
             server.  Runs where we combine several hundred datasets execute in a couple minutes \
             (on a server).  Your mileage could vary, and we are interested in cases where \
             it is underperforming.  Let us know.  We are aware of several things that could be \
             done to speed up enclone.\n\n",
        )?;

        h.print(
            "\\boldblue{4. How does enclone fit into the 10x Genomics software ecosystem?}\n\n",
        )?;
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
        )?;

        h.print("\\boldblue{5. What platforms does enclone run on?}\n\n")?;
        h.print(
            "1. Linux/x86-64 (that's most servers)\n\
               2. Mac.\n\n\
             However, we have not and cannot test every possible configuration of these \
             platforms.  Please let us know if you encounter problems!\n\n",
        )?;

        h.print("\\boldblue{6. How can I print out all the donor reference sequences?}\n\n")?;
        h.print(
            "Add the argument \\bold{DONOR_REF_FILE=filename} to your enclone command, \
             and fasta for the donor reference sequences will be dumped there.\n\n",
        )?;

        h.print(
            "\\boldblue{7. How does enclone know what VDJ reference sequences I'm using?}\n\n",
        )?;
        h.print(
            "If you used Cell Ranger version 4.0 or greater, then the VDJ reference file was \
             included in the outs directory, and so enclone knows the reference sequence from \
             that.\n\n\
             For outs from older Cell Ranger versions, enclone has to guess which VDJ \
             reference sequences were used, and may or may not do so correctly.  As part of this, \
             if you have mouse data from older Cell Ranger versions, you need to supply the \
             argument \\bold{MOUSE} on the command line.\n\n\
             It is also possible to set the reference sequence directly by adding \
             by adding \\bold{REF=f} to your command line, where \\bold{f} is the name of your \
             VDJ reference fasta file, but if that is different than the reference \
             supplied to Cell Ranger, then you will have to add the additional argument \
             \\bold{RE} to recompute annotations, and that will slow down enclone somewhat.\n\n",
        )?;

        h.print("\\boldblue{8. Can I provide data from more than one donor?}\n\n")?;
        h.print(
            "Yes.  Type \\bold{enclone help input} for details.  The default behavior of \
             enclone is to prevent cells from different donors from being placed in the same \
             clonotype.  The \\bold{MIX_DONORS} option may be used to turn off this behavior.  If \
             you employ this option, then clonotypes containing cells from more than one donor \
             will be flagged as errors, unless you use the \\bold{NWARN} option to turn off those \
             warnings.  The primary reason for allowing entry of data from multiple \
             donors is to allow estimation of enclone's error rate.\n\n",
        )?;

        h.print("\\boldblue{9. What are some command line argument values quoted?}\n\n")?;
        h.print(
            "Command line argument values that contain any of these characters ;|* need to \
             be quoted like so\n\
             \\bold{TCR=\"a;b\"}\n\
             to prevent the shell from interpreting them for a purpose completely unrelated to \
             enclone.  This is a trap, because forgetting to add the quotes can result in \
             nonsensical and confusing behavior!\n\n",
        )?;

        h.print("\\boldblue{10. If enclone fails, does it return nonzero exit status?}\n\n")?;
        h.print(
            "Yes, unless output of enclone is going to a terminal.  In that case, you'll always \
             get zero.\n\n",
        )?;

        h.print("\\boldblue{11. Could a cell be missing from an enclone clonotype?}\n\n")?;
        h.print(
            "Yes, some cells are deliberately deleted.  The cell might have been deleted by \
             one of the filters described in \\bold{enclone help special}, and which you can \
             turn off.  We also delete cells for which more than four chains were found.\n\n",
        )?;

        h.print("\\boldblue{12. Can enclone print summary stats?}\n\n")?;
        h.print(
            "Yes, if you add the option \\bold{SUMMARY}, then some summary stats will be \
             printed.  If you wish to suppress visual output, then also add the option \
             \\bold{NOPRINT}.\n\n",
        )?;

        h.print("\\boldblue{13. What is the notes column?}\n\n")?;
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
        )?;

        h.print("\\boldblue{14. Can I cap the number of threads used by enclone?}\n\n")?;
        h.print(
            "You can use the command-line argument \\bold{MAX_CORES=n} to cap the number of \
             cores used in parallel loops.  The number of threads used is typically one \
             higher.\n\n",
        )?;

        h.print("\\boldblue{15. Can I use enclone if I have only gene expression data?}\n\n")?;
        h.print(
            "Possibly.  In some cases this works very well, but in other cases it does not.  \
            Success depends on dataset characteristics that have not been carefully investigated.  \
            To attempt this, you need to invoke Cell Ranger on the GEX dataset as if \
            it was a VDJ dataset, and you need to specify to Cell Ranger that the run is to be \
            treated as BCR or TCR.  Two separate invocations can be used to get both.  Note also \
            that Cell Ranger has been only minimally tested for this configuration and that this \
            is not an officially supported Cell Ranger configuration.\n\n",
        )?;

        h.print("\\boldblue{16. How can I cite enclone?}\n\n")?;
        let mut log = Vec::<u8>::new();
        emit_green_escape(&mut log);
        h.print(&strme(&log).to_string())?;
        if !ctl.gen_opt.stable_doc {
            h.print(&format!(
                "10x Genomics, https://github.com/10XGenomics/enclone,\nversion {}.\n",
                version_string()
            ))?;
        } else {
            h.print(
                "10x Genomics, https://github.com/10XGenomics/enclone,\n\
                    (your enclone version information will be printed here).\n",
            )?;
        }
        let mut log = Vec::<u8>::new();
        emit_end_escape(&mut log);
        h.print(&strme(&log).to_string())?;
        h.print(
            "At some point in the future, there will be a white paper to which you can refer, \
            in addition to a DOI minted at Zenodo.  In the spirit of reproducibility, you should \
            provide the arguments that you used when you ran enclone and indicate the version of \
            Cell Ranger that you used to generate the input data.\n\n",
        )?;

        h.print("\\boldblue{17. Can I print the enclone version?}\n\n")?;
        h.print("Yes, type \"enclone version\".\n\n")?;

        h.print("\\boldblue{18. Can enclone ingest multiple datasets from the same library?}\n\n")?;
        h.print(
            "If enclone detects significant (≥ 25%) barcode reuse between datasets, it will exit.  \
            This behavior can be overridden using the argument \\bold{ACCEPT_REUSE}.\n\n"
        )?;

        h.print("\\boldblue{19. Can I turn off all the filters used in joining clonotypes?}\n\n")?;
        h.print(
            "Pretty much.  You can run with the following arguments:\n\
            MAX_CDR3_DIFFS=100\n\
            MAX_LOG_SCORE=100\n\
            EASY\n\
            MAX_DIFFS=200\n\
            MAX_DEGRADATION=150,\n\
            however this will in general be very slow and not produce useful results.  Depending \
            on what your goal is, you may find it helpful to use some of these arguments, and \
            with lower values.  You can see the meaning of the arguments and their default values \
            by typing \"enclone help how\".\n\n",
        )?;

        h.print("\\boldblue{20. How can I send the developers an example?}\n\n")?;
        h.print(
            "Use filters to select a clonotype or clonotypes of interest.  \
            Then you can cut and paste enclone output into an email.  If you want the example to \
            be reproducible by us, add the argument\n\
            \\bold{SUBSET_JSON=filename} to the command line, which will create a json file \
            containing only data for the barcodes in the clonotype.  Then send us the file \
            as an email attachment.  This only works for VDJ data, and we do not have a \
            parallel mechanism for gene expression and antibody data.  Please note also that \
            running enclone on the barcodes from a single clonotype will not necessarily \
            reproduce the results you observed, because the clonotyping algorithm uses all the \
            data, even if only some clonotypes are selected.\n\n",
        )?;

        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Finish enclone help all.

    if h.help_all {
        h.dump();

    // Catch unrecognized help requests.
    } else if args.len() >= 2 && !h.ok {
        let mut x = args[1].clone();
        x.make_ascii_lowercase();
        if x.contains("help") {
            return Err("\nYour help request doesn't match one known to enclone.\n\
                Please type \"enclone\" to see the help options.\n"
                .to_string());
        }
    }
    Ok(())
}
