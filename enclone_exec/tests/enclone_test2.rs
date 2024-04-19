// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// There are three categories of tests here:
// 1. basic tests (feature = basic), runs without additional data requirements
// 2. nonbasic tests, requires extended dataset distributed with enclone
// 3. speed test (feature = cpu), requires non-public datasets.

use ansi_escape::*;
use anyhow::Error;
use enclone_core::*;
use enclone_proto::proto_io::{read_proto, ClonotypeIter};
use enclone_proto::types::EncloneOutputs;
use enclone_testlist::main_testlist::*;
use enclone_testlist::*;
use enclone_tools::html::*;
use enclone_tools::run_test::*;
use flate2::read::GzDecoder;
use io_utils::*;
use itertools::Itertools;

use rayon::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};
use stats_utils::*;
use std::cmp::min;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::{metadata, read_dir, read_to_string, remove_dir_all, remove_file, File};
use std::io;
use std::io::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::process::{Command, Stdio};
use std::thread;
use std::time;
use std::time::{Duration, Instant};
use string_utils::*;
use vector_utils::*;

const LOUPE_OUT_FILENAME: &str = "testx/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 10. The following is a single test, containing many subtests, each of which is a regression test
// for a given enclone command line.
//
// If you ever need to change the output of all tests, use the main program
// update_all_main_tests.rs in enclone/src/bin.  Note that there is some duplicated code there.
//
// If this code is modified to also include the DTESTS, the test time is almost exactly the same.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone() {
    run_tests(env!("CARGO_BIN_EXE_enclone"), "test", 0, &TESTS);
}

#[cfg(not(feature = "cpu"))]
#[test]
/// Tests that are affected by the D region alignment algorithm.
fn test_enclone_d() {
    run_tests(
        env!("CARGO_BIN_EXE_enclone"),
        "dtest",
        40,
        &[
            (
                1,
                "test ALIGN_2ND<n>",
                r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW ALIGN_2ND1 CVARS=d1_name"###,
            ),
            (
                2,
                "test JALIGN_2ND<n>",
                r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW JALIGN_2ND1 CVARS=d2_name"###,
            ),
            (
                3,
                "test ALIGN_JALIGN_CONSISTENCY",
                r###"BCR=123085 CELLS=1 CHAINS=2 ALIGN1 JALIGN1 ALIGN_JALIGN_CONSISTENCY AMINO=cdr3
             PLAIN NOPAGER EXPECT_OK"###,
            ),
            (
                4,
                "test D_INCONSISTENT, and lock number of inconsistencies",
                r###"BCR=123085 D_INCONSISTENT CVARS=d1_name COMPLETE NGROUP"###,
            ),
            (
                5,
                "the JALIGN1 in this example had a boundary location that was off by one",
                r###"BCR=165807 JALIGN1 AMINO=cdr3 CVARS=d1_score,d2_score CDR3=CAKEYYDFWSGYSDVRGVIPNIDYW"###,
            ),
            (
                6,
                "the JALIGN1 in this example had a boundary location that was off by one",
                r###"BCR=123085 CELLS=2 JALIGN1 AMINO=cdr3 CVARS=d1_name CDR3=CAKAGPTESGYYVWYFDLW"###,
            ),
            (
                7,
                "test d_inconsistent_{%,n}",
                r###"BCR=123085 GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT SUMMARY SUMMARY_CLEAN"###,
            ),
            (
                8,
                "test ALIGN<n>",
                r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW ALIGN1 CVARS=d1_name"###,
            ),
            (
                9,
                "test ALIGN<n> and JALIGN<n>, case where there's a D segment",
                r###"BCR=85333 ALIGN1 JALIGN1 CDR3=CARGYDFWSGYLVGNWAGDYYYYMDVW"###,
            ),
            (
                10,
                "test ALIGN<n> and JALIGN<n>, case where there is no D segment",
                r###"BCR=85333 ALIGN1 JALIGN1 CDR3=CAKGKGFRNYYYYMDVW"###,
            ),
            (
                11,
                "test d1 etc.",
                r###"BCR=123085 CVARS=d1_name,d2_name,d_Δ,d_delta AMINO=cdr3 CDR3=CARVRDILTGDYGMDVW"###,
            ),
            (
                12,
                "test GROUP_VDJ_REFNAME_HEAVY (deprecated but supported)",
                r###"BCR=86237 GROUP_VDJ_REFNAME_HEAVY CDR3="CAKAVAGKAVAGGWDYW|CAKVSTGIAVAGPGDYW" COMPLETE"###,
            ),
            (
                13,
                "test GROUP_VJ_REFNAME_HEAVY (deprecated but supported)",
                r###"BCR=86237 GROUP_VJ_REFNAME_HEAVY CDR3="CARGVLWFGELGAFDIW|CARAGLGVVLAARGAFDIW""###,
            ),
            (
                14,
                "test placement of indel, needed shifting right",
                r###"BCR=123085 CELLS=1 CHAINS=2 AMINO=cdr3 JALIGN2 CDR3=CAKDKSRPPTHYYGSGSYYSRILDNW"###,
            ),
            (
                15,
                "test placement of indel, needed shifting left",
                r###"BCR=123085 CELLS=1 CHAINS=2 AMINO=cdr3 JALIGN2 CDR3=CARMAQFYSGSGTYYIGPYYFEYW"###,
            ),
        ],
    );
}

#[cfg(not(feature = "cpu"))]
#[test]
/// Tests that are affected by the grouping algorithm.
fn test_grouping() {
    run_tests(
        env!("CARGO_BIN_EXE_enclone"),
        "gtest",
        0,
        &[
            (
                1,
                "test 5/8 for newline correctness (this grouping option deprecated but supported)",
                r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN SET_IN_STONE"###,
            ),
            (
                2,
                "test 6/8 for newline correctness (this grouping option deprecated but supported)",
                r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN NGROUP SET_IN_STONE"###,
            ),
            (
                3,
                "test 7/8 for newline correctness (this grouping option deprecated but supported)",
                r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML SET_IN_STONE"###,
            ),
            (
                4,
                "test 8/8 for newline correctness (this grouping option deprecated but supported)",
                r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML NGROUP SET_IN_STONE"###,
            ),
            (
                5,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_refname,cdr3_aa_heavy≥80%,cdr3_aa_light≥80% CVARS=cdr3_len
             AMINO=cdr3 CDR3="CARHLQWELP.*W""###,
            ),
            (
                6,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_refname,len,cdr3_len MIN_GROUP=2 MIN_CHAINS=2 CDR3="CQQSY.*TLATF"
             CVARS=cdr3_len"###,
            ),
            (
                7,
                "test of GROUP",
                r###"BCR=123085 GROUP=cdr3_aa_heavy≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
             CDR3=CARPKSDYIIDAFDIW"###,
            ),
            (
                8,
                "test of GROUP",
                r###"BCR=123085 GROUP=cdr3_aa_light≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
             CDR3=CQTWGTGPWVF"###,
            ),
            (
                9,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_refname,aa_heavy≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
             CDR3=CARVPYYYDRSYYYYGMDVW"###,
            ),
            (
                10,
                "test of AGROUP",
                r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CARHSYSSGWYDEWDYW
             AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=top=2"###,
            ),
            (
                11,
                "test of AGROUP",
                r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CAKDGGEHYYDSSGYYASYYFDYW 
             AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=14"###,
            ),
            (
                12,
                "test of AGROUP",
                r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CAKDGGEHYYDSSGYYASYYFDYW 
             AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=13"###,
            ),
            (
                13,
                "test of AGROUP",
                r###"BCR=123085 AGROUP AG_CENTER=copy_filters MIN_CELLS=2 MAX_CELLS=2
             AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=3 MIN_GROUP=2"###,
            ),
            (
                14,
                "test symmetric grouping stats",
                r###"BCR=123085 GROUP=vj_refname,cdr3_aa_heavy≥80%,cdr3_aa_light≥80% NOPRINT
             SUMMARY SUMMARY_CLEAN"###,
            ),
            (
                15,
                "test of GROUP",
                r###"BCR=123085 GROUP=cdr3_heavy≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
             CDR3=CARPKSDYIIDAFDIW"###,
            ),
            (
                16,
                "test of GROUP",
                r###"BCR=123085 GROUP="cdr3_light>=100%" MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
             CDR3=CQTWGTGPWVF"###,
            ),
            (
                17,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_refname,heavy≥96.6% MIN_GROUP=2 MIN_CHAINS=2 
             CDR3="CARVIVGPKKLEGRLYSSSLHFDCW|CARVIVGPEKQEGRLYSSSLHFDYW" POUT=stdout PCOLS=vj_seq1"###,
            ),
            (
                18,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_heavy_refname,cdr3_heavy_len,cdr3_heavy≥80% LVARS=n,donors,dref 
             CVARS=const,cdr3_len AMINO=cdr3 CHAINS=2 MIN_GROUP=2
             CDR3="CARDLHGYDPYGMDVW|CARELRHYDTYGMDVW""###,
            ),
            (
                19,
                "test of GROUP",
                r###"BCR=123085 GROUP=vj_refname,cdr3_light_len LVARS=n,donors,dref CVARS=const,cdr3_len 
             AMINO=cdr3 CHAINS=2 MIN_GROUP=2 CDR3="CARESAVAGDMDVW|CARDYGDYRWWVDGMDVW""###,
            ),
            (
                20,
                "test of group",
                r###"BCR=123085 GROUP=vj_refname GROUP_CDR3=CACFGRIGVVVRAAHYW"###,
            ),
            (
                21,
                "test of group, asserted at one time",
                r###"BCR=123085 GROUP=vj_refname,cdr3_len MIN_GROUP=3 HONEY=out=stdout,color=var,u1 
             EXPECT_OK"###,
            ),
        ],
    );
}

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
/// Test using datasets that are either in the extended public dataset collection,
/// or which require samtools.
fn test_extended() {
    run_tests(
        env!("CARGO_BIN_EXE_enclone"),
        "ext_test",
        0,
        &[
            (1, "Make sure that POUT works on full dataset. If we experience failures on other PUBLIC ids, we can add them to this list.",
            r###"BCR=86237 RE POUT=/dev/null NOPRINT EXPECT_OK NO_PRE"###),
            (2, "tests nd2",
            r###"BCR=47199,47200,47212 AMINO=cdr3 NCROSS LVARS=nd2 CDR3=CVKGKSGSFWYYFENW
             NO_PRE"###),
            (3, "test sec and mem [requires samtools]",
            r###"BCR=123085 GEX=123217 LVARSP=sec,mem CDR3=CVKDRVTGTITELDYW"###),
            (4, "crashed at one point",
            r###"BCR=128037,128040 GEX=127798,127801 LVARSP=pe1 NOPRINT EXPECT_OK NO_PRE"###),
            //
            (5, "this added because it got better when a bug in bads detection was fixed",
            r###"TCR=163914 CDR3=CASRLGGEETQYF NO_PRE"###),
            (6, "Test PCHAINS=max.  For this we need a clonotype having at least five chains, and the \
                question is whether the header line represents cvars for all the chains.  The output of
                this is expected to change whenever variables are added.",
            r###"BCR=123085,123089,124547 NWEAK_CHAINS NDOUBLET MIN_CHAINS=5 POUT=stdout PCHAINS=max
             NOPRINT RE NO_PRE"###),
            (7, "test MIN_GROUP_DONORS",
            r###"BCR="40953;43899" MIX_DONORS MIN_GROUP=2 NO_PRE
             GROUP="cdr3_len,cdr3_aa_heavy>=85%,cdr3_aa_light>=85%,vj_refname" MIN_GROUP_DONORS=2"###),
            (8, "this asserted at one point",
            r###"BUILT_IN GROUP=vj_refname,cdr3_aa_heavy≥90% MIN_CHAINS_EXACT=2 MIN_GROUP=2 
             KEEP_CLONO_IF_CELL_MEAN="cdr3_len1>=18" BCR=1018096-1018098 JALIGN1 NO_PRE
             EXPECT_OK"###),
            (9, "this clonotype included a junk chain before we made a change, and test /outs",
            r###"TCR=163911/outs CDR3=CAPSAGDKIIF AMINO=donor NO_PRE"###),
            (10, "test case where digit rows are just barely present",
            r###"TCR=163911 CDR3=CASSLVQPSTDTQYF AMINO=donor NO_PRE"###),
            (11, "this added because it got better when a noise filter was added, also tests u_max",
            r###"TCR=163914 CDR3=CASSLVQPSTDTQYF CVARSP=u_max NO_PRE"###),
            (12, "this added because it got better when a noise filter was added; also test FASTA",
            r###"TCR=163914 CDR3=CAFRGGSYIPTF FASTA=stdout NO_PRE"###),
        ],
    );
}

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
/// Crash tests.  These are tests to make sure that certain options do not result in a crash, even
/// when run on a large and complex dataset.  The options are in groups because not all are
/// compatible with each other.  The datasets are defined by a single fixed list, to be enlarged
/// over time based on discovery of pathologies in particular PUBLIC datasets.
/// All run with certain shared options.
fn test_crash() {
    let crash_tests: Vec<_> = [
        (
            1,
            "CONP SEQC SUM MEAN BARCODES DIFF_STYLE=C1 GROUP_VJ_REFNAME",
        ),
        (
            2,
            "CONX FULL_SEQC DIFF_STYLE=C2 POUT=stdout PCOLS=count_CAR",
        ),
        (
            3,
            "AMINO=fwr1,cdr1,fwr2,cdr2,fwr3,cdr3,fwr4 CVARS=d1_name,d2_name,d_delta,d_Δ,cigar",
        ),
        (
            4,
            "PLOT_BY_ISOTYPE=stdout MIN_CELLS=3 GROUP_VJ_REFNAME_HEAVY ALIGN1 JALIGN1",
        ),
        (
            5,
            "GROUP_VDJ_REFNAME_HEAVY GVARS=d_inconsistent_%,d_inconsistent_n",
        ),
        (6, "GROUP=vj_refname,cdr3_aa_heavy≥90%,cdr3_aa_light≥90%"),
    ]
    .iter()
    .map(|(num, crash_set)| {
        (
            num,
            "crash test",
            format!(
            "BCR=\"45977;123085;testx/inputs/flaky\" {crash_set} NOPRINT BUILT_IN EXPECT_OK NO_PRE",
        ),
        )
    })
    .collect();
    run_tests(
        env!("CARGO_BIN_EXE_enclone"),
        "crash_test",
        40,
        &crash_tests
            .iter()
            .map(|(num, comment, args)| (**num, *comment, args.as_str()))
            .collect_vec(),
    );
}

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
/// Regression tests for internal features.
fn test_internal() {
    run_tests(
        env!("CARGO_BIN_EXE_enclone"),
        "internal_test",
        40,
        &[
            (1, "gave wrong result",
            r###"123085 CDR3=CARDRIAGRFGYGMDVW"###),
            (2, "test human + IMGT; note that specifying by number forces BCR+TCR reference checks",
            r###"123085 REQUIRE_UNBROKEN_OK IMGT ACCEPT_BROKEN EXPECT_NULL"###),
            (3, "this crashed; it is not exactly an internal feature test but uses an internal feature (IMGT) to exhibit the phenomenon",
            r###"BCR=123085 IMGT RE ACCEPT_BROKEN POUT=stdout PCELL BARCODE=AGCAGCCCATTAGGCT-1 EXPECT_OK"###),
        ],
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 14.

// NOT BASIC

// Test site for broken links and spellcheck.
//
// Two approaches for checking broken links left in place for now, to delete one, and the
// corresponding crate from Cargo.toml.
//
// This looks for
// ▓<a href="..."▓
// ▓http:...[, '")}<#\n]▓
// ▓https:...[, '")}<#\n]▓
// ▓<img src="..."▓.
// (These also test termination by ". ".)
// SHOULD also look for at least:
// ▓ href="..."▓
// ▓ href='...'▓.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_for_broken_links_and_spellcheck() {
    use std::time::Duration;

    // Set up link exceptions.  These are links that have been observed to break periodically.
    // The web.archive.org ones are probably just too slow, and we should allow for that rather
    // than have it on the unreliable list.  The "period" version is because of a parsing bug.
    // See also the test in ./test.

    let unreliable_links = include_str!("../../pages/unreliable_links")
        .split('\n')
        .collect::<Vec<&str>>();

    // Set up list of archived links.  These are broken a lot and we have archived versions, so
    // we don't test these at all.  If we determine that they're permanently broken, we should
    // do something different.

    let archived_links = [
        "http://www.bioinf.org.uk/abs/info.html#martinnum",
        "http://opig.stats.ox.ac.uk/webapps/stcrdab",
        "http://www.imgt.org",
    ];

    // Set up dictionary exceptions.  We should rewrite the code to avoid looking in certain
    // places and reduce the dictionary exceptions accordingly.

    let extra_words =
        "abybank actgtgcgagag actgtgcgagagc adefghiklmnpqrstvwy amazonaws anarci autoremove \
        barcode barcodes barcoding bcn \
        bioinf biorxiv cdiff cellranger chmod clen clono clonotype clonotypes \
        clonotyping codebase colorn contig contigs cqvwdsssdhpyvf cred crispr cshlp \
        csv ctrlc cvar cvars datalayer dejavusansmono dotplot \
        dref dyiid enclone exe executables false fcell \
        fixedtextbox foursie foursies frameshifted frameshifts frontiersin fwr fwyh ganesh \
        genomics germline ggctttgactactgg gggctttgactactgg github githubusercontent google \
        googletagmanager grok gz hcomp html \
        hypermutation hypermutations igblast igh ighd igk igl ighm igkc imgt immunoglobulins \
        indel indels inkt intradonor ireceptor \
        jsdelivr json krh levenshtein lgc linux loh lvar lvars \
        macbook mait metadata minmax mkdir \
        moresies multiomic nall ncbi nchains ncross ndoublet newick nimproper \
        nopager noprint nospaces nqual nseg nsegn nsig nwhitef oligos onesie onesies osx parseable \
        pbmc pcell pcols pdb pgas phad phylip png \
        plasmablast powershell preinstalled prepends pwm pwms recombinants redownloads \
        researchsquare rustup samtools screenshot segn \
        sloooooooow spacebar stackexchange standalone stcrdab stdout sthnqedkr subclonotype \
        subclonotypes sudo svg tattgtagtggtggtagct tctgtgcgagata tctgtgcgagat tctgtgcgagata \
        testlist thresholding timeline timepoint \
        tracebacks trb tsv \
        tttctgtgcgaga tttctgtgcgagat \
        twosie ubuntu udiff umi umis underperforming unicode untarring \
        vddj vdj vdjc vilella vilfwym vilm vjlen wallclock website wget whitef whitelist wikimedia \
        wikipedia workaround workflow xcode xf xhtml xkcd \
        xxxxxxxxxxx xxxxxxxxxxxxxxxxxxxxxxx xy yvar zenodo zx";
    let extra_words = extra_words.split(' ').collect::<Vec<&str>>();

    // Set up dictionary.

    let dictionary0 = read_to_string("../enclone-data/english_wordlist").unwrap();
    let dictionary0 = dictionary0.split('\n').collect::<Vec<&str>>();
    let mut dictionary = Vec::<String>::new();
    for w in dictionary0.iter() {
        let mut x = w.to_string();
        x.make_ascii_lowercase();
        dictionary.push(x);
    }
    for w in extra_words {
        dictionary.push(w.to_string());
    }
    unique_sort(&mut dictionary);

    // Find html pages on site.

    let mut htmls = vec!["../index.html".to_string()];
    let pages = read_dir("../pages").unwrap();
    for page in pages {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(page.to_string());
        }
    }
    let auto = read_dir("../pages/auto").unwrap();
    for page in auto {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(page.to_string());
        }
    }

    // Hardcoded exceptions to link testing, because of slowness.

    let mut tested = HashSet::<String>::new();
    tested.insert("https://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd".to_string());
    tested.insert("http://www.w3.org/1999/xhtml".to_string());

    // Hardcode exception for funny svn URL.

    tested.insert("https://github.com/10XGenomics/enclone/trunk".to_string());

    // Test each html.

    let mut dict_fail = false;
    for x in htmls.iter() {
        let mut bads = HashSet::<String>::new();
        let f = open_for_read![x];
        let depth = x.matches('/').count();
        let mut line_no = 0;
        for line in f.lines() {
            line_no += 1;
            let mut s = line.unwrap();

            // Test spelling.  Case insensitive.

            let mut s0 = s.replace(',', " ");
            s0 = s0.replace('.', " ");
            s0 = s0.replace(';', " ");
            let words = s0.split(' ').collect::<Vec<&str>>();
            for i in 0..words.len() {
                let mut ok = true;
                let w = words[i].to_string();
                for c in w.chars() {
                    if !c.is_ascii_alphabetic() {
                        ok = false;
                    }
                }
                if w.is_empty() || !ok {
                    continue;
                }
                let mut wl = w.clone();
                wl.make_ascii_lowercase();
                if !bin_member(&dictionary, &wl.to_string()) && !bads.contains(&wl.to_string()) {
                    let mut wu = w.clone();
                    wu.make_ascii_uppercase();
                    if w != wu || w.len() < 20 {
                        // arbitrary long uppercase strings allowed
                        bads.insert(wl.to_string());
                        eprintln!(
                            "\nthe word \"{}\" in file {} isn't in the dictionary",
                            wl, x
                        );
                        dict_fail = true;
                    }
                }
            }

            // Check links.

            if cfg!(feature = "linkless") {
                continue;
            }
            let mut links = Vec::<String>::new();
            let mut chars = Vec::<char>::new();
            for c in s.chars() {
                chars.push(c);
            }
            let mut i = 0;
            let terminators = [',', ' ', '\'', '"', ')', '}', '<', '#', '\n'];
            while i < chars.len() {
                let http = chars[i..].starts_with(&['h', 't', 't', 'p', ':']);
                let https = chars[i..].starts_with(&['h', 't', 't', 'p', 's', ':']);
                if http || https {
                    for j in i + 5..chars.len() {
                        if terminators.contains(&chars[j])
                            || (chars[j] == '.' && j < chars.len() - 1 && chars[j + 1] == ' ')
                        {
                            let mut link = String::new();
                            for k in i..j {
                                link.push(chars[k]);
                            }
                            if !tested.contains(&link.to_string()) {
                                links.push(link.clone());
                                tested.insert(link.to_string());
                            }
                            i = j - 1;
                            break;
                        }
                    }
                }
                i += 1;
            }
            let s2 = s.clone();
            while s.contains("<a href=\"") {
                let link = s.between("<a href=\"", "\"");
                if tested.contains(&link.to_string()) {
                    s = s.after("<a href=\"").to_string();
                    continue;
                }
                tested.insert(link.to_string());

                // Allow mailto to enclone.

                if link == "mailto:enclone@10xgenomics.com" {
                    s = s.after("<a href=\"").to_string();
                    continue;
                }

                // Otherwise if not http..., assume it's a file path.

                if !link.starts_with("http") {
                    let mut link = link.to_string();
                    if link.contains('#') {
                        link = link.before("#").to_string();
                    }
                    let mut z = link.clone();
                    for _ in 0..depth - 1 {
                        if !z.starts_with("../") {
                            eprintln!("something wrong with file {} on page {}", link, x);
                            panic!("failed");
                        }
                        z = z.after("../").to_string();
                    }
                    z = format!("../{}", z);
                    if !path_exists(&z) {
                        eprintln!("failed to find file {} on page {}", link, x);
                        panic!("failed");
                    }
                    s = s.after("<a href=\"").to_string();
                    continue;
                }

                // And finally do http....

                links.push(link.to_string());
                s = s.after("<a href=\"").to_string();
            }
            s = s2;
            while s.contains("<img src=\"") {
                let path = s.between("<img src=\"", "\"");
                if tested.contains(&path.to_string()) {
                    s = s.after("<img src=\"").to_string();
                    continue;
                }
                tested.insert(path.to_string());
                let path = path.to_string();
                let mut z = path.clone();
                for _ in 0..depth - 1 {
                    if !path.starts_with("../") {
                        eprintln!("something wrong with file {} on page {}", path, x);
                        panic!("failed");
                    }
                    z = z.after("../").to_string();
                }
                z = format!("../{}", z);
                if !path_exists(&z) {
                    eprintln!("failed to find file {} on page {}", path, x);
                    panic!("failed");
                }
                s = s.after("<img src=\"").to_string();
            }
            'links: for link in links {
                // Temporary workaround.

                if link == "https://10xgenomics.github.io/enclone/install.sh" {
                    continue;
                }

                // Test for known unreliable links.

                let mut unreliable = false;
                for l in unreliable_links.iter() {
                    if *l == link {
                        unreliable = true;
                    }
                }
                if unreliable {
                    continue;
                }
                let mut archived = false;
                for l in archived_links.iter() {
                    if *l == link {
                        archived = true;
                    }
                }
                if archived {
                    continue;
                }

                // Test for some links that don't exist yet, but will exist once page is live.

                for h in htmls.iter() {
                    if link == format!("https://10xgenomics.github.io/enclone/{}", h.after("../")) {
                        continue 'links;
                    }
                }

                // eprintln!("checking link \"{}\"", link);

                // Approach 1 to testing if link works.  This seemed to hang once in spite of
                // the timeout.

                use attohttpc::*;
                const LINK_RETRIES: usize = 5;
                for i in 0..LINK_RETRIES {
                    if i > 0 {
                        thread::sleep(time::Duration::from_millis(100));
                        eprintln!("retrying link {}, attempt {}", link, i);
                    }
                    let req = attohttpc::get(link.clone()).read_timeout(Duration::new(10, 0));
                    let response = req.send();
                    if response.is_err() {
                        eprintln!(
                            "\ncould not read link {} on page {} line {}\n",
                            link, x, line_no
                        );
                        if i == LINK_RETRIES - 1 {
                            panic!("failed");
                        }
                    } else {
                        let response = response.unwrap();
                        if response.is_success() {
                            break;
                        }
                        eprintln!(
                            "\ncould not read link {} on page {} line {}\n",
                            link, x, line_no
                        );
                        if i == LINK_RETRIES - 1 {
                            panic!("failed");
                        }
                    }
                }

                // Approach 2 to testing if link works.  This may not have a timeout and does
                // not auto retry like approach 1.  Also may not compile anymore.

                /*
                use reqwest::StatusCode;
                let req = reqwest::blocking::get(link);
                if req.is_err() {
                    eprintln!("\ncould not read link {} on page {}\n", link, x);
                    panic!("failed");
                }
                if req.unwrap().status() == StatusCode::NOT_FOUND {
                    eprintln!("\ncould not read link {} on page {}\n", link, x);
                    panic!("failed");
                }
                */
            }
        }
    }
    if dict_fail {
        eprintln!();
        panic!("failed");
    }
}
