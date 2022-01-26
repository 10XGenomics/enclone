// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// There are three categories of tests here:
// 1. basic tests (feature = basic), runs without additional data requirements
// 2. nonbasic tests, requires extended dataset distributed with enclone
// 3. speed test (feature = cpu), requires non-public datasets.

use ansi_escape::*;
use anyhow::Error;
use enclone_core::main_testlist::*;
use enclone_core::testlist::*;
use enclone_core::*;
use enclone_proto::proto_io::{read_proto, ClonotypeIter};
use enclone_proto::types::EncloneOutputs;
use enclone_tools::html::*;
use enclone_tools::run_test::*;
use flate2::read::GzDecoder;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
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
    PrettyTrace::new().on();
    let t = Instant::now();
    println!("running tests using {}", env!("CARGO_BIN_EXE_enclone"));
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        results.push((i, false, String::new()));
    }
    let this = include_str!("../../enclone_core/src/main_testlist.rs");
    let mut tracking = false;
    let mut comments = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in this.lines() {
        if line.starts_with("pub const TESTS: ") {
            tracking = true;
            continue;
        }
        if tracking {
            if line == "];" {
                break;
            }
            lines.push(line.to_string());
        }
    }
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with("    // ") {
                let c = lines[j].after("    // ").as_bytes()[0];
                if c >= b'0' && c <= b'9' {
                    break;
                }
            }
            j += 1;
        }
        comments.push(lines[i..j].iter().format("\n").to_string());
        i = j;
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = TESTS[it].to_string();
        if test.split(' ').collect::<Vec<&str>>().contains(&"NFORCE") {
            eprintln!(
                "\nOn main tests, NFORCE is not allowed, because it can cause \
                failure in the GitHub Actions tests.\n"
            );
            std::process::exit(1);
        }
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            &comments[it],
            &test,
            "test",
            &mut res.1,
            &mut res.2,
            &mut out,
            0,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ntotal time for {} enclone subtests = {:.2} seconds\n",
        TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 10a.  unaccounted time test

#[cfg(not(feature = "basic"))]
#[cfg(feature = "cpu")]
#[test]
fn test_accounting() {
    PrettyTrace::new().on();
    let t = Instant::now();
    println!("running tests using {}", env!("CARGO_BIN_EXE_enclone"));
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..UNAC_TESTS.len() {
        results.push((i, false, String::new()));
    }
    let this = include_str!("../../enclone_core/src/testlist.rs");
    let mut tracking = false;
    let mut comments = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in this.lines() {
        if line.starts_with("pub const UNAC_TESTS: ") {
            tracking = true;
            continue;
        }
        if tracking {
            if line == "];" {
                break;
            }
            lines.push(line.to_string());
        }
    }
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with("    // ") {
                let c = lines[j].after("    // ").as_bytes()[0];
                if c >= b'0' && c <= b'9' {
                    break;
                }
            }
            j += 1;
        }
        comments.push(lines[i..j].iter().format("\n").to_string());
        i = j;
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = UNAC_TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            &comments[it],
            &test,
            "unac_test",
            &mut res.1,
            &mut res.2,
            &mut out,
            0,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ntotal time for {} enclone unaccounted subtests = {:.2} seconds\n",
        UNAC_TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 10b.  tests that are affected by D region assignment

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_d() {
    PrettyTrace::new().on();
    let t = Instant::now();
    println!("running tests using {}", env!("CARGO_BIN_EXE_enclone"));
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..DTESTS.len() {
        results.push((i, false, String::new()));
    }
    let this = include_str!("../../enclone_core/src/testlist.rs");
    let mut tracking = false;
    let mut comments = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in this.lines() {
        if line.starts_with("pub const DTESTS: ") {
            tracking = true;
            continue;
        }
        if tracking {
            if line == "];" {
                break;
            }
            lines.push(line.to_string());
        }
    }
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with("    // ") {
                let c = lines[j].after("    // ").as_bytes()[0];
                if c >= b'0' && c <= b'9' {
                    break;
                }
            }
            j += 1;
        }
        comments.push(lines[i..j].iter().format("\n").to_string());
        i = j;
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = DTESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            &comments[it],
            &test,
            "dtest",
            &mut res.1,
            &mut res.2,
            &mut out,
            40,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ntotal time for {} enclone d gene subtests = {:.2} seconds\n",
        DTESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 10c.  tests that are affected grouping

#[cfg(not(feature = "cpu"))]
#[test]
fn test_grouping() {
    PrettyTrace::new().on();
    let t = Instant::now();
    println!("running tests using {}", env!("CARGO_BIN_EXE_enclone"));
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..GTESTS.len() {
        results.push((i, false, String::new()));
    }
    let this = include_str!("../../enclone_core/src/testlist.rs");
    let mut tracking = false;
    let mut comments = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in this.lines() {
        if line.starts_with("pub const GTESTS: ") {
            tracking = true;
            continue;
        }
        if tracking {
            if line == "];" {
                break;
            }
            lines.push(line.to_string());
        }
    }
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with("    // ") {
                let c = lines[j].after("    // ").as_bytes()[0];
                if c >= b'0' && c <= b'9' {
                    break;
                }
            }
            j += 1;
        }
        comments.push(lines[i..j].iter().format("\n").to_string());
        i = j;
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = GTESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            &comments[it],
            &test,
            "gtest",
            &mut res.1,
            &mut res.2,
            &mut out,
            0,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ntotal time for {} enclone grouping subtests = {:.2} seconds\n",
        GTESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 11.

// NOT BASIC

// Regression tests using the extended public dataset collection.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_extended() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..EXTENDED_TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = EXTENDED_TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "ext_test",
            &mut res.1,
            &mut res.2,
            &mut out,
            0,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\nextended tests total time for {} enclone subtests = {:.2} seconds\n",
        EXTENDED_TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 12.

// NOT BASIC

// Crash tests.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_crash() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mut crash_tests = Vec::<String>::new();
    for i in 0..CRASH_SETS.len() {
        crash_tests.push(format!("{} {} {}", CRASH_DATA, CRASH_SETS[i], CRASH_OPTS));
    }
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..crash_tests.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = &crash_tests[it];
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "crash_test",
            &mut res.1,
            &mut res.2,
            &mut out,
            40,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ncrash tests total time for {} enclone subtests = {:.2} seconds\n",
        crash_tests.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 13.

// NOT BASIC

// Regression tests for internal features.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_internal() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..INTERNAL_TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = INTERNAL_TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "internal_test",
            &mut res.1,
            &mut res.2,
            &mut out,
            40,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ninternal tests total time for {} enclone subtests = {:.2} seconds\n",
        INTERNAL_TESTS.len(),
        elapsed(&t)
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
        "abybank actgtgcgagag actgtgcgagagc adefghiklmnpqrstvwy amazonaws anarci barcode \
        barcodes barcoding bcn \
        bioinf cdiff cellranger chmod clen clono clonotype clonotypes \
        clonotyping codebase colorn contig contigs cqvwdsssdhpyvf cred crispr cshlp \
        csv ctrlc cvar cvars datalayer dejavusansmono dotplot \
        dref dyiid enclone executables false fcell \
        fixedtextbox foursie foursies frameshifted frameshifts frontiersin fwr fwyh ganesh \
        genomics germline ggctttgactactgg gggctttgactactgg github githubusercontent google \
        googletagmanager grok gz html \
        hypermutation hypermutations igblast igh ighd igk igl ighm igkc imgt immunoglobulins \
        indel indels inkt intradonor jsdelivr json krh levenshtein lgc linux loh lvar lvars \
        macbook mait metadata minmax mkdir \
        moresies multiomic nall ncbi nchains ncross ndoublet newick nimproper \
        nopager noprint nospaces nqual nseg nsegn nsig nwhitef oligos onesie onesies parseable \
        pbmc pcell pcols pdb pgas phad phylip png \
        plasmablast preinstalled prepends pwm pwms recombinants redownloads \
        researchsquare samtools screenshot segn \
        sloooooooow spacebar stackexchange standalone stcrdab stdout sthnqedkr subclonotype \
        subclonotypes svg tattgtagtggtggtagct tctgtgcgagata tctgtgcgagat tctgtgcgagata testlist \
        thresholding timepoint \
        tracebacks trb tsv \
        tttctgtgcgaga \
        tttctgtgcgagat \
        twosie ubuntu udiff umi umis underperforming unicode untarring \
        vddj vdj vdjc vilella vilfwym vilm vjlen website wget whitef whitelist wikimedia \
        wikipedia workaround workflow xf xhtml xkcd \
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
            htmls.push(format!("{}", page));
        }
    }
    let auto = read_dir("../pages/auto").unwrap();
    for page in auto {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(format!("{}", page));
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
            let terminators = vec![',', ' ', '\'', '"', ')', '}', '<', '#', '\n'];
            while i < chars.len() {
                let http = chars[i..].starts_with(&vec!['h', 't', 't', 'p', ':']);
                let https = chars[i..].starts_with(&vec!['h', 't', 't', 'p', 's', ':']);
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
                            std::process::exit(1);
                        }
                        z = z.after("../").to_string();
                    }
                    z = format!("../{}", z);
                    if !path_exists(&z) {
                        eprintln!("failed to find file {} on page {}", link, x);
                        std::process::exit(1);
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
                        std::process::exit(1);
                    }
                    z = z.after("../").to_string();
                }
                z = format!("../{}", z);
                if !path_exists(&z) {
                    eprintln!("failed to find file {} on page {}", path, x);
                    std::process::exit(1);
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
                            std::process::exit(1);
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
                            std::process::exit(1);
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
                    std::process::exit(1);
                }
                if req.unwrap().status() == StatusCode::NOT_FOUND {
                    eprintln!("\ncould not read link {} on page {}\n", link, x);
                    std::process::exit(1);
                }
                */
            }
        }
    }
    if dict_fail {
        eprintln!("");
        std::process::exit(1);
    }
}
