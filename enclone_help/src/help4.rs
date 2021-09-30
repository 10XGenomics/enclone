// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Test for help request, under development.

use crate::help_utils::*;
use tables::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn help4(args: &Vec<String>, mut h: &mut HelpDesk) -> Result<(), String> {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide special filtering help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "special") || h.help_all {
        h.begin_doc("special")?;
        h.print("\n\\bold{special filtering options}\n\n")?;
        h.print(
            "This page documents some options, most of which allow noise \
             filters to be turned off, and which normally should not be invoked.  Some of \
             these options delete barcodes, and a summary of this action is included in the \
             SUMMARY option.  See also the lead variable \"filter\", see \
             \"enclone help lvars\".  At the bottom of this page we provide some other options \
             that are not noise filters.\n\n",
        )?;

        h.docf2(
            "NALL",
            "Turn off all the noise filters shown below.  This may yield quite a mess.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NCELL",
            "Use contigs found by Cell Ranger even if they were not in a called cell, \
            or not called high confidence.",
            55,
        )?;

        h.doc(
            "NALL_CELL",
            "Turn off all the noise filters except for the cell filter.",
        );

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NGEX",
            "If gene expression and/or feature barcode data are provided, if a barcode \
            is called a cell by the VDJ part of the Cell Ranger pipeline, but not \
            called a cell by the gene expression and/or feature barcode part, then the \
            default behavior of enclone is to remove such cells from clonotypes.  This \
            option disables that behavior.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NCROSS",
            "If you specify that two or more libraries arose from the same origin (i.e. \
            cells from the same tube or tissue), then by default enclone will \
            \"cross filter\" so as to remove expanded exact subclonotypes that are \
            present in one library but not another, in a fashion that would be highly \
            improbable, assuming random draws of cells from the tube.  These are \
            believed to arise when a plasma or plasmablast cell breaks up during during \
            or after pipetting from the tube, and the resulting fragments seed GEMs, \
            yielding expanded 'fake' clonotypes that are residues of real single plasma \
            cells.  The NCROSS options turns off this filter, which could be useful so \
            long as you interpret the restored clonotypes as representing what are \
            probably single cells.  There may also be other situations where the filter \
            should be turned off, and in particular the filter can do weird things if \
            inputs are somehow mis-specified to enclone.  Note that for purposes of \
            this option, enclone defines an origin by the pair \
            (origin name, donor name).",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NUMI",
            "Filter out B cells based on low BCR UMI counts.  The heuristics",
            65,
        )?;

        h.docpr(
            "",
            "for this are described on the enclone site at \\green{bit.ly/enclone}.",
        );
        h.doc(
            "NUMI_RATIO",
            "Filter out B cells based on low BCR UMI counts relative to another",
        );
        h.doc2("cell in a given clonotype.  The heuristics for this");
        h.docpr(
            "",
            "are described on the enclone site at \\green{bit.ly/enclone}.",
        );

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NGRAPH_FILTER",
            "By default, enclone filters to remove exact subclonotypes that by virtue of \
            their relationship to other exact subclonotypes, appear to arise from \
            background mRNA or a phenotypically similar phenomenon.  The \
            NGRAPH_FILTER option turns off this filtering.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NQUAL",
            "By default, enclone filters out exact subclonotypes having a base in V..J \
            that looks like it might be wrong.  More specifically, enclone finds bases \
            which are not Q60 for a barcode, not Q40 for two barcodes, are not \
            supported by other exact subclonotypes, are variant within the clonotype, \
            and which disagree with the donor reference.  NQUAL turns this off.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NWEAK_CHAINS",
            "By default, enclone filters chains from clonotypes that are \
            weak and appear to be artifacts, perhaps arising from a stray mRNA molecule \
            that floated into a GEM.  The NWEAK_CHAINS option turns off this filter.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NWEAK_ONESIES",
            "By default, enclone disintegrates certain untrusted clonotypes into single cell \
            clonotypes.  The untrusted clonotypes are onesies \
            that are light chain or TRA and whose number of cells is less than 0.1% of the total \
            number of cells.  This operation reduces the likelihood of creating clonotypes \
            containing cells that arose from different recombination events.  NWEAK_ONESIES turns \
            this operation off.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NMERGE_ONESIES",
            "enclone merges certain onesie clonotypes into clonotypes having two or more chains.  \
            By default, this merger is prevented if the number of cells in the onesie is less \
            than 0.01% of the total number of cells.  NMERGE_ONESIES causes these merges to \
            happen anyway.  The naming of this option is confusing.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NFOURSIE_KILL",
            "Under certain circumstances, enclone will delete foursie exact subclonotypes.  \
            Please see 10xgenomics.github.io/enclone/pages/auto/default_filters.html.  \
            The foursies that are killed are believed to be artifacts arising \
            from repeated cell doublets or GEMs that contain two cells and multiple gel \
            beads.  The argument NFOURSIE_KILL turns off this filtering.",
            62,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NDOUBLET",
            "Under certain circumstances, enclone will delete exact subclonotypes that appear \
            to represent doublets.  \
            Please see 10xgenomics.github.io/enclone/pages/auto/default_filters.html.  \
            The argument NDOUBLET turns off this filtering.",
            65,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NSIG",
            "Under certain circumstances, enclone will delete exact subclonotypes that appear \
            to be contaminants, based on their chain signature.  \
            Please see 10xgenomics.github.io/enclone/pages/auto/default_filters.html.  \
            The argument NSIG turns off this filtering.",
            65,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NWHITEF",
            "By default, enclone filters out rare artifacts arising from contamination \
            of oligos on gel beads.  The NWHITEF option turns off this filter.",
            55,
        )?;

        h.ldoc(
            "NBC_DUP",
            "By default, enclone filters out duplicated barcodes within an exact",
        );
        h.doc(
            "",
            "subclonotype.  The NBC_DUP option turns off this filter.",
        );

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "MIX_DONORS",
            "By default, enclone will prevent cells from different donors from being \
            placed in the same clonotype.  The MIX_DONORS option turns off this \
            behavior, thus allowing cells from different donors to be placed in the \
            same clonotype.  The main use of this option is for specificity testing, \
            in which data from different donors are deliberately combined in an attempt \
            to find errors.  Use of the bc field for META input specification \
            automatically turns on this option.",
            55,
        )?;

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "NIMPROPER",
            "enclone filters out exact subclonotypes having more than one chain, but all of the \
            same type.  For example, the filter removes all exact subclonotypes having two TRA \
            chains and no other chains.  The NIMPROPER option turns off this filter.",
            55,
        )?;

        // Documentation section.

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "MIN_CHAINS_EXACT=n",
            "Delete any exact subclonotype having less than n chains.  You can use this \
            to \"purify\" a clonotype so as to display only exact subclonotypes having \
            all their chains.",
            55,
        )?;

        h.doc(
            "CHAINS_EXACT=n",
            "Delete any exact subclonotype not having exactly n chains.",
        );
        h.doc(
            "MIN_CELLS_EXACT=n",
            "Delete any exact subclonotype having less than n cells.  You might",
        );
        h.doc2("want to use this if you have a very large and complex expanded.");
        h.doc2("clonotype.");
        h.doc(
            "COMPLETE",
            "delete any exact subclonotype that has less chains than the",
        );
        h.doc2("clonotype for which you would like to see a simplified view.");
        h.docf2(
            "CONST_IGH=\"<pattern>\"",
            "for BCR, keep only exact subclonotypes having a heavy chain whose constant region \
            gene name matches the given pattern (meaning regular expression, see \
            \"enclone help filter\")",
            55,
        )?;
        h.docf2(
            "CONST_IGKL=\"<pattern>\"",
            "for BCR, keep only exact subclonotypes having a light chain whose constant region \
            gene name matches the given pattern (meaning regular expression, see \
            \"enclone help filter\")",
            55,
        )?;
        h.doc(
            "MAX_HEAVIES=1",
            "ignore any cell having more than one IGH or TRB chain",
        );

        // Documentation section.

        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "KEEP_CELL_IF=constraint",
            "Let \"constraint\" be any constraint involving arithmetic and \
            boolean operators, and variables that are specified as fields using the BC option \
            (or equivalently, using bc, via META), see \"enclone help input\", or feature \
            variables: <gene>_g or <antibody>_ab or <crispr>_cr or <custom>_cu, as described at \
            \"enclone help lvars\" (but without regular expressions, as these would conflict \
            with arithmetic operators).  This \
            option filters out all barcodes that do not satisfy the given constraint.  \
            Note that for purposes of testing the constraint, if the value for a \
            particular barcode has not been specified, then its value is \
            taken to be null.  Also multiple instances of KEEP_CELL_IF may be used to impose \
            multiple filters.  See the examples below, and be very careful about syntax, \
            which should match the given examples exactly.  In particular,",
            55,
        )?;
        h.doc2("• use == for equality, and not =");
        h.doc2("• put string values in single quotes");
        h.doc2("• put the entire expression in double quotes.");
        h.doc2("");
        h.doc(
            "",
            "As a toy example, suppose you had a CSV file f having five lines:",
        );
        h.doc2("barcode,nice,rank");
        h.doc2("AGCATACTCAGAGGTG-1,true,3");
        h.doc2("CGTGAGCGTATATGGA-1,true,7");
        h.doc2("CGTTAGAAGGAGTAGA-1,false,99");
        h.doc2("CGTTAGAAGGAGTAGA-1,dunno,43");
        h.doc2("then the command");
        h.doc2("enclone BCR=123085 BC=f KEEP_CELL_IF=\"nice == 'true'\"");
        h.doc2("would cause enclone to use only the first two barcodes shown in");
        h.doc2("the file, and the command");
        h.doc2("enclone BCR=123085 BC=f KEEP_CELL_IF=\"nice == 'true' && rank <= 5\"");
        h.doc2("would cause only the first barcode to be used.");
        h.doc2("");
        h.doc2("See also KEEP_CLONO_IF_CELL_MEAN=... and");
        h.doc2("KEEP_CLONO_IF_CELL_MAX=... at \"enclone help filter\".");

        // Done.

        h.print_tab2()?;
        h.print("\n")?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Function that provides an explanation used for both enclone help lvars and
    // enclone help cvars.

    fn explain_alt_versions(h: &mut HelpDesk) -> Result<(), String> {
        h.print(&gray_left_bar(&print_to(
            "\\red{●} These variables have some alternate versions, \
                 as shown in the table below.\n\n",
        )))?;
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "variable".to_string(),
            "semantics".to_string(),
            "visual".to_string(),
            "visual".to_string(),
            "parseable".to_string(),
            "parseable".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "".to_string(),
            "".to_string(),
            "".to_string(),
            "(one cell)".to_string(),
            "".to_string(),
            "(one cell)".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 6];
        rows.push(row);
        let row = vec![
            "x".to_string(),
            "median over cells".to_string(),
            "yes".to_string(),
            "this cell".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_mean".to_string(),
            "mean over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_μ".to_string(),
            "(same as above)".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_sum".to_string(),
            "sum over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_Σ".to_string(),
            "(same as above)".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_min".to_string(),
            "min over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_max".to_string(),
            "max over cells".to_string(),
            "yes".to_string(),
            "null".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_%".to_string(),
            "% of total GEX (genes only)".to_string(),
            "yes".to_string(),
            "this cell".to_string(),
            "yes".to_string(),
            "yes".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "x_cell".to_string(),
            "this cell".to_string(),
            "no".to_string(),
            "no".to_string(),
            "no".to_string(),
            "this cell".to_string(),
        ];
        rows.push(row);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|l|l|l|l|l".to_vec(), false, false);
        h.print_plain(&gray_left_bar(&log))?;
        h.print_plain(&gray_left_bar(&print_to(
            "Some explanation is required.  If you use enclone without certain options, you \
             get the \"visual\" column.\n\
             • Add the option \\bold{PER_CELL} \
             (see \"enclone help display\") and then you get visual output with extra lines for \
             each cell within an exact subclonotype, and each of those extra lines is described by \
             the \"visual (one cell)\" column.\n\
             • If you generate parseable output (see \"enclone help parseable\"), then you get \
             the \"parseable\" column for that output, unless you specify \\bold{PCELL}, \
             and then you get the last column.\n\
             • For the forms with μ and Σ, the Greek letters are only used in column headings for \
             visual output (to save space), and optionally, in names of fields on the command \
             line.\n\
             \\green{▶} If you try out these features, you'll see exactly what happens! \
             \\green{◀}\n",
        )))?;
        Ok(())
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide lvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "lvars") || h.help_all {
        h.begin_doc("lvars")?;
        h.print("\n\\bold{lead column options}\n\n")?;
        h.print(
            "These options define lead variables, which are variables that are computed for each \
             exact subclonotype, and if using the \\bold{PER_CELL} option, also computed for each \
             cell.  In addition, lead variables can be used for parseable output.\n\n\
             Lead variables appear in columns that \
             appear once in each clonotype, on the left side, and have one entry for each \
             exact subclonotype row.\n\n\
             Note that for medians of integers, \
             we actually report the \"rounded median\", the result of rounding the \
             true median up to the nearest integer, so that e.g. 6.5 is rounded up to 7.\n\n",
        )?;
        h.print(
            "See also \"enclone help cvars\" and the inventory of all variables at
            https://10xgenomics.github.io/enclone/pages/auto/inventory.html.\n\n",
        )?;
        h.print(
            "Lead variables are specified using \\bold{LVARS=x1,...,xn} \
             where each xi is one of:\n\n",
        )?;
        h.doc("nchains", "total number of chains in the clonotype");
        h.doc(
            "nchains_present",
            "number of chains present in an exact subclonotype",
        );
        h.ldoc("datasets", "dataset identifiers");
        h.doc("origin", "origin identifiers");
        h.doc("donors", "donor identifiers");
        h.ldoc("n", "number of cells");
        h.doc(
            "n_<name>",
            "number of cells associated to the given name, which can be a dataset",
        );
        h.doc2("or origin or donor or tag short name; may name only one such category");
        h.doc("clonotype_ncells", "total number of cells in the clonotype");
        h.ldoc(
            "nd<k>",
            "For k a positive integer, this creates k+1 fields, that are specific to each",
        );
        h.doc2("clonotype.  The first field is n_<d1>, where d1 is the name of the dataset");
        h.doc2("having the most cells in the clonotype.  If k ≥ 2, then you'll get a");
        h.doc2("\"runner-up\" field n_<d2>, etc.  Finally you get a field n_other, however");
        h.doc2("fields will be elided if they represent no cells.  Use a variable of this");
        h.doc2("type at most once.");
        h.ldoc(
            "near",
            "Hamming distance of V..J DNA sequence to nearest neighbor",
        );
        h.doc(
            "far",
            "Hamming distance of V..J DNA sequence to farthest neighbor",
        );
        h.doc2("both compare to cells having chains in the same columns of the clonotype,");
        h.doc2("with - shown if there is no other exact subclonotype to compare to");
        h.doc(
            "dref",
            "Hamming distance of V..J DNA sequence to donor reference, excluding",
        );
        h.doc2("region of recombination");
        h.doc(
            "dref_aa",
            "Hamming distance of V..J amino acid sequence to donor reference, excluding",
        );
        h.doc2("region of recombination");
        h.ldoc(
            "count_<reg>",
            "Number of matches of the V..J amino acid sequences of all chains to the given",
        );
        h.doc2("regular expression, which is treated as a subset match, so for example,");
        h.doc2("count_CAR would count the total number of occurrences of the string CAR in");
        h.doc2("all the chains.  Please see \"enclone help filter\" for a discussion");
        h.doc2("about regular expressions.  We also allow the form abbr:count_<regex>,");
        h.doc2("where abbr is an abbreviation that will appear as the field label.");
        h.doc(
            "count_<f>_<reg>",
            "Supposing that f is in {cdr1,..,cdr3,fwr1,..,fwr4,cdr,fwr}, this is similar",
        );
        h.doc2("to the above but restricted to motifs lying entirely within");
        h.doc2("a given feature or feature set.");
        h.ldoc(
            "inkt",
            "A string showing the extent to which the T cells in an exact subclonotype",
        );
        h.doc2("have evidence for being an iNKT cell.  The most evidence is denoted 𝝰gj𝝱gj,");
        h.doc2("representing both gene name and junction sequence (CDR3) requirements for");
        h.docpr(
            "",
            "both chains.  See \\green{bit.ly/enclone} for details on the requirements.",
        );
        h.doc("mait", "Same as with inkt but for MAIT cells instead.");
        h.ldoc(
            "g<d>",
            "Here d is a nonnegative integer.  Then all the exact subclonotypes are",
        );
        h.doc2("grouped according to the Hamming distance of their V..J sequences.  Those");
        h.doc2("within distance d are defined to be in the same group, and this is");
        h.doc2("extended transitively.  The group identifier 1, 2, ... is shown.  The");
        h.doc2("ordering of these identifiers is arbitrary.  This option is best applied");
        h.doc2("to cases where all exact subclonotypes have a complete set of chains.");
        h.ldocpr("gex", "\\red{●} median gene expression UMI count");
        h.docpr(
            "n_gex",
            "\\blue{●} number of cells found by cellranger using GEX or Ab data",
        );
        // nonpublic for now as we don't know if this is useful
        /*
        h.doc(
            "entropy",
            "Shannon entropy of GEX UMI counts (median across cells)"
        );
        */
        h.ldocpr(
            "<gene>_g",
            "\\red{●} all four feature types: look for a declared feature of the \
             given type",
        );
        h.doc(
            "<antibody>_ab",
            "with the given id or name; report the median UMI count for it; we allow",
        );
        h.doc(
            "<crispr>_cr",
            "we also allow <regular expression>_g where g can be replaced by ab, ag, cr",
        );
        h.doc(
            "<custom>_cu",
            "or cu; this represents a sum of UMI counts across the matching features. ●",
        );

        // sec and mem: deprecated because not enough signal

        /*
        h.ldoc(
            "sec",
            "for human or mouse BCR, number of GEX UMIs that are characterized as secreted",
        );
        h.doc(
            "mem",
            "for human or mouse BCR, number of GEX UMIs that are characterized as membrane",
        );
        h.doc2(
            "For both of these, the algorithm looks for reads that are aligned through the",
        );
        h.doc2(
            "right end of a constant region CH3 exon, and then read into a CH3-CHS or",
        );
        h.doc2(
            "CH4-CHS exon, in the secreted case, or a M, M1 or M2 exon, in the membrane case.",
        );
        h.doc2(
            "This choice is determined by sequence tables in the code, and we cannot be",
        );
        h.doc2("absolutely certain that these tables are complete.");
        h.docpr(
            "",
            "\\bold{These fields require the presence of the files possorted_genome_bam.bam}",
        );
        h.docpr("", "\\bold{and possorted_genome_bam.bam.bai.}");
        h.docpr(
            "",
            "\\bold{These fields also require that you have samtools in your path.}",
        );
        h.doc2("Note that these counts tend to be low.");
        h.docpr(
            "",
            "\\boldred{PLEASE NOTE: THIS IS EXPERIMENTAL AND UNLIKELY TO BE FULLY CORRECT.}",
        );
        */

        h.ldoc(
            "cred",
            "Short for credibility.  It is a measure of the extent to which cells",
        );
        h.doc2("having gene expression similar to a given putative B cell are themselves");
        h.doc2("B cells.  (Or similarly for T cells.)  For the actual definition, let n");
        h.doc2("be the number of VDJ cells that are also GEX cells.  For a given cell,");
        h.doc2("find the n GEX cells that are closest to it in PCA space, and report the");
        h.doc2("percent of those that are also VDJ cells.  For multiple datasets, it would");
        h.doc2("be better to \"aggr\" the data, however that is not currently supported");
        h.doc2("The computation is also inefficient, so let us know if it's causing");
        h.doc2("problems for you.  And cred makes much better sense for datasets that");
        h.doc2("consist of mixed cell types, rather than consisting of pure B or T cells.");
        h.rows.push(vec!["\\hline".to_string(); 2]);
        h.docf2(
            "filter",
            "See \"enclone help special\".  Use with PER_CELL.  If you turn off some default \
            filters (or all default filters, e.g. with NALL_CELL), and this cell would have been \
            deleted by one of the default filters, then this will show the name of \
            the last filter that would have been applied to delete the cell.  Note that there \
            are complex interactions between filters, so the actual effect with all default \
            filters on may be significantly different.  Note also that use of NALL_CELL will \
            typically result in peculiar artifacts, so this should only be used as an \
            exploratory tool.",
            75,
        )?;
        h.print_tab2()?;
        h.print(
            "For gene expression and feature barcode stats, such data must be provided \
             as input to enclone.\n\n",
        )?;
        h.print(
            "● Example: IG.*_g matches all genes that begin with IG, and TR(A|B).*_g matches \
             all genes that begin with TRA or TRB.  Double quotes as in \\bold{LVARS=\"...\"} \
             may be needed.  The regular expression must \
             be in the alphabet A-Za-z0-9+_-.[]()|* and is only interpreted as a regular \
             expression if it contains a character in []()|*.  \
             See \"enclone help filter\" \
             for more information about regular expressions.\n\n",
        )?;
        explain_alt_versions(&mut h)?;
        h.print(
            "\n\\blue{●} Similar to the above but simpler: n_gex is just a count of cells, \
             visual (one cell) shows 0 or 1, n_gex_cell is defined for parseable (one cell), \
             and the x_mean etc. forms do not apply.\n\n",
        )?;
        h.print(
            "The default is \\bold{datasets,n}, except that datasets is suppressed if \
             there is only one dataset.\n\n",
        )?;
        h.print("\\bold{LVARSP=x1,...,xn} is like \\bold{LVARS} but appends to the list.\n\n")?;
        h.end_doc();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Provide cvars help.

    if (args.len() == 3 && args[1] == "help" && args[2] == "cvars") || h.help_all {
        h.begin_doc("cvars")?;

        // Header.

        h.print(
            "\n\\bold{per-chain column options}: These options define per-chain variables, \
             which correspond to columns that appear once for each chain in each clonotype, and \
             have one entry for each exact subclonotype.  Please note that for medians of \
             integers, we actually report the \"rounded median\", the result of rounding the \
             true median up to the nearest integer, so that e.g. 6.5 is rounded up to 7.\n\n",
        )?;
        h.print(
            "See also \"enclone help lvars\" and the inventory of all variables at
            https://10xgenomics.github.io/enclone/pages/auto/inventory.html.\n\n",
        )?;
        h.print(
            "Per-column variables are specified using\n\
             \\bold{CVARS=x1,...,xn}\n\
             where each xi is one of:\n\n",
        )?;

        // Main table entries.

        h.doc(
            "var",
            "bases at positions in chain that vary across the clonotype",
        );
        h.ldocpr(
            "u",
            "\\red{●} VDJ UMI count for each exact subclonotype, median across cells",
        );
        h.docpr(
            "r",
            "\\red{●} VDJ read count for each exact subclonotype, median across cells",
        );
        h.ldoc(
            "edit",
            "a string that defines the edit of the reference V(D)J concatenation versus",
        );
        h.doc2("the contig, from the beginning of the CDR3 to the end of the J segment;");
        h.doc2("this uses a coordinate system in which 0 is the first base of the J ref");
        h.doc2("segment (or the first base of the D ref segment for IGH and TRB); for");
        h.doc2("example D-4:4 denotes the deletion of the last 4 bases of the V segment, ");
        h.doc2("I0:2 denotes an insertion of 2 bases after the V");
        h.doc2("and I0:2•S5 denotes that plus a substitution at position 5; in computing");
        h.doc2("\"edit\", for IGH and TRB, we always test every possible D segment,");
        h.doc2("regardless of whether one is annotated, and pick the best one; for this");
        h.doc2("reason, \"edit\" may be slow");
        h.doc(
            "comp",
            "a measure of CDR3 complexity, which is the total number of S, D and I",
        );
        h.doc2("symbols in \"edit\" as defined above");
        h.doc(
            "cigar",
            "the CIGAR string that defines the edit of the V..J contig sequence versus",
        );
        h.doc2("the universal reference V(D)J concatenation");

        h.ldoc(
            "cdr*_aa",
            "the CDR*_AA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_aa_L_R_ext",
            "the CDR*_AA sequence, with L amino acids added on the left and R amino acids",
        );
        h.doc2("added on the right; either may be negative, denoting trimming instead");
        h.doc2("of extension");
        h.doc(
            "cdr*_aa_north",
            "the CDR*_AA sequence for BCR defined by North B et al. (2011), A new",
        );
        h.doc(
            "",
            "clustering of antibody CDR loop conformations, J Mol Biol 406, 228-256.",
        );
        h.doc("", "cdr1_aa_north = cdr1_aa_3_3_ext for heavy chains");
        h.doc("", "cdr1_aa_north = cdr1_aa for light chains");
        h.doc("", "cdr2_aa_north = cdr2_aa_2_3_ext for heavy chains");
        h.doc("", "cdr2_aa_north = cdr2_aa_1_0_ext for light chains");
        h.doc("", "cdr3_aa_north = cdr3_aa_-1_-1_ext");
        h.doc(
            "cdr*_aa_ref",
            "cdr*_aa, for the universal reference sequence (but not for cdr3)",
        );
        h.doc(
            "cdr*_len",
            "number of amino acids in the CDR* sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_dna",
            "the CDR*_DNA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "cdr*_dna_ref",
            "same, for the universal reference sequence (but not for cdr3)",
        );
        h.ldoc(
            "cdr3_aa_conx",
            "consensus for CDR3 across the clonotype, showing X for each variant residue",
        );
        h.doc(
            "cdr3_aa_conp",
            "consensus for CDR3 across the clonotype, showing a property symbol whenever",
        );
        h.doc2("two different amino acids are observed, per the following table:");
        h.doc2("--------------------------------------------------------------------");
        h.doc2("asparagine or aspartic acid   B   DN");
        h.doc2("glutamine or glutamic acid    Z   EQ");
        h.doc2("leucine or isoleucine         J   IL");
        h.doc2("negatively charged            -   DE");
        h.doc2("positively charged            +   KRH");
        h.doc2("aliphatic (non-aromatic)      Ψ   VILM");
        h.doc2("small                         π   PGAS");
        h.doc2("aromatic                      Ω   FWYH");
        h.doc2("hydrophobic                   Φ   VILFWYM");
        h.doc2("hydrophilic                   ζ   STHNQEDKR");
        h.doc2("any                           X   ADEFGHIKLMNPQRSTVWY");
        h.doc2("--------------------------------------------------------------------");
        h.doc2("The table is searched top to bottom until a matching class is found.");
        h.doc2("In the special case where every amino acid is shown as a gap (-),");
        h.doc2("a \"g\" is printed.");
        h.ldoc(
            "fwr*_aa",
            "the FWR*_AA sequence, or \"unknown\" if not computed",
        );
        h.doc("fwr*_aa_ref", "same, for the universal reference sequence");
        h.doc(
            "fwr*_len",
            "number of amino acids in the FWR* sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "fwr*_dna",
            "the FWR*_DNA sequence, or \"unknown\" if not computed",
        );
        h.doc(
            "fwr*_dna_ref",
            "same, for the universal reference sequences",
        );
        h.doc(
            "",
            "For all of these, * is 1 or 2 or 3 (or 4, for the fwr variables).",
        );
        h.doc(
            "",
            "For CDR1 and CDR2, please see \"enclone help amino\" and the page on",
        );
        h.docpr("", "\\green{bit.ly/enclone} on V(D)J features.");
        h.ldoc("v_start", "start of V segment on full DNA sequence");
        h.doc(
            "d_start",
            "start of D segment on full DNA sequence (or null)",
        );
        h.doc(
            "cdr3_start",
            "base position start of CDR3 sequence on full contig",
        );
        h.doc(
            "d_frame",
            "reading frame of D segment, either 0 or 1 or 2 (or null)",
        );
        h.ldoc(
            "aa%",
            "amino acid percent identity with donor reference, outside junction region",
        );
        h.doc(
            "dna%",
            "nucleotide percent identity with donor reference, outside junction region",
        );

        h.ldoc("utr_name", "name of 5'-UTR region");
        h.doc("v_name", "name of V region");
        h.doc("d_name", "name of D region (or null)");
        h.doc("j_name", "name of J region");
        h.doc("const", "name of constant region");
        h.ldoc("utr_id", "id of 5'-UTR region");
        h.doc("v_id", "id of V region");
        h.doc("d_id", "id of D region (or null)");
        h.doc("j_id", "id of J region");
        h.doc("const_id", "id of constant region (or null, if not known)");
        h.doc2("(these are the numbers after \">\" in the VDJ reference file)");
        h.ldoc("d1_name", "name of optimal D gene, or none");
        h.doc("d2_name", "name of second best D gene, or none");
        h.doc("d1_score", "score for optimal D gene");
        h.doc("d2_score", "score for second best D gene");
        h.doc(
            "d_delta",
            "score difference between first and second best D gene",
        );
        h.doc("d_Δ", "same");
        h.doc2("These are recomputed from scratch and ignore the given assignment.");
        h.doc2("Note that in many cases D gene assignments are essentially random, as");
        h.doc2("it is often not possible to know the true D gene assignment.");
        h.doc2("If the value is \"null\" it means that having no D gene at all scores better");
        h.ldoc(
            "vjlen",
            "number of bases from the start of the V region to the end of the J region",
        );
        h.doc2("Please note that D gene assignments are frequently \"random\" -- it is not");
        h.doc2("possible to know the actual D gene that was assigned.");
        h.doc(
            "clen",
            "length of observed constant region (usually truncated at primer start)",
        );
        h.doc("ulen", "length of observed 5'-UTR sequence;");
        h.doc(
            "",
            "note however that what report is just the start of the V segment",
        );
        h.doc(
            "",
            "on the contig, and thus the length may include junk before the UTR",
        );
        h.doc(
            "cdiff",
            "differences with universal reference constant region, shown in the",
        );
        h.doc(
            "",
            "abbreviated form e.g. 22T (ref changed to T at base 22) or 22T+10",
        );
        h.doc(
            "",
            "(same but contig has 10 additional bases beyond end of ref C region",
        );
        h.doc(
            "",
            "At most five differences are shown, and if there are more, ... is appended.",
        );
        h.doc("udiff", "like cdiff, but for the 5'-UTR");
        h.ldoc("q<n>_", "comma-separated list of the quality");
        h.doc2("scores at zero-based position n, numbered starting at the");
        h.doc2("beginning of the V segment, for each cell in the exact subclonotype");
        h.ldoc(
            "notes",
            "optional note if there is an insertion or the end of J does not exactly abut",
        );
        h.doc(
            "",
            "the beginning of C; elided if empty; also single base overlaps between",
        );
        h.docpr(
            "",
            "J and C are not shown unless you use the special option \\bold{JC1}; we do this",
        );
        h.doc(
            "",
            "because with some VDJ references, one nearly always has such an overlap",
        );
        h.ldoc(
            "ndiff<n>vj",
            "number of base differences within V..J between this exact subclonotype and",
        );
        h.doc("", "exact subclonotype n");
        h.doc(
            "d_univ",
            "distance from universal reference, more specifically,",
        );
        h.doc(
            "",
            "number of base differences within V..J between this exact",
        );
        h.doc(
            "",
            "clonotype and universal reference, exclusive of indels, the last 15",
        );
        h.doc("", "bases of the V and the first 15 bases of the J");
        h.doc("d_donor", "distance from donor reference,");
        h.doc("", "as above but computed using donor reference");

        // The rest.

        h.print_tab2()?;
        h.print("\n")?;
        explain_alt_versions(&mut h)?;
        h.print(
            "\nAt least one variable must be listed.  The default is \\bold{u,const,notes}.  \
             \\bold{CVARSP}: same as \\bold{CVARS} but appends.\n\n",
        )?;
        h.end_doc();
    }
    Ok(())
}
