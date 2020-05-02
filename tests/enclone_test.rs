// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// somewhat misplaced comments:
//
// Run enclone on a test case and verify that the output is identical to what
// was gotten before.  If the output is different, look at it
// and decide if the change is justified, and if so update the output file.
//
// To test just this test, use:
//
// cargo test -p enclone enclone -- --nocapture

use ansi_escape::*;
use enclone::html::insert_html;
use enclone::misc3::parse_bsv;
use enclone::proto_io::read_proto;
use enclone::testlist::*;
use enclone::types::EncloneOutputs;
use failure::Error;
use flate2::read::GzDecoder;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashSet;
use std::fs::{read_dir, read_to_string, remove_file, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::process::{Command, Stdio};
use std::thread;
use std::time;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

const LOUPE_OUT_FILENAME: &str = "test/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that files are rustfmt'ed.

#[test]
fn test_formatting() {
    let mut rs = Vec::<String>::new();
    let src = read_dir("src").unwrap();
    for x in src {
        let x = x.unwrap().path();
        let x = x.to_str().unwrap();
        if x.ends_with(".rs") {
            rs.push(x.to_string());
        }
    }
    rs.push("build.rs".to_string());
    rs.push("tests/enclone_test.rs".to_string());
    use itertools::Itertools;
    let new = Command::new("rustfmt")
        .arg("--check")
        .args(&rs)
        .output()
        .expect(&format!("failed to execute test_formatting"));
    if new.status.code().unwrap() != 0 {
        eprintln!("\nYou need to run rustfmt.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// The following is a single test, containing many subtests, each of which is a regression test
// for a given enclone command line.
//
// If you ever need to change the output of all tests, use the main program
// update_all_main_tests.rs in enclone/src/bin.  Note that there is some duplicated code there.

#[test]
fn test_enclone() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let mut test = TESTS[it].to_string();
        test = test.replace("\n", "");
        for _ in 0..3 {
            test = test.replace("  ", " ");
        }
        let mut expect_null = false;
        let mut expect_fail = false;
        let mut expect_ok = false;
        let mut set_in_stone = false;
        let mut no_pre = false;
        if test.contains(" EXPECT_NULL") {
            test = test.replace(" EXPECT_NULL", "");
            expect_null = true;
        }
        if test.contains(" EXPECT_FAIL") {
            test = test.replace(" EXPECT_FAIL", "");
            expect_fail = true;
        }
        if test.contains(" EXPECT_OK") {
            test = test.replace(" EXPECT_OK", "");
            expect_ok = true;
        }
        if test.contains(" SET_IN_STONE") {
            test = test.replace(" SET_IN_STONE", "");
            set_in_stone = true;
        }
        if test.contains(" NO_PRE") {
            test = test.replace(" NO_PRE", "");
            no_pre = true;
        }
        test = test.replace("{TEST_FILES_VERSION}", &format!("{}", TEST_FILES_VERSION));
        let mut log = Vec::<u8>::new();
        let out_file = format!("test/inputs/outputs/enclone_test{}_output", it + 1);
        let mut pre_arg = format!("PRE=test/inputs/version{}", TEST_FILES_VERSION);
        if no_pre {
            pre_arg = String::new();
        }
        if !path_exists(&out_file) && !expect_fail && !expect_ok {
            fwriteln!(log, "\nYou need to create the output file {}.\n", out_file);
            fwriteln!(
                log,
                "Do this by executing the following command from \
                 cellranger/lib/rust/enclone:\n"
            );
            emit_bold_escape(&mut log);
            fwriteln!(
                log,
                "enclone {} {} > test/inputs/outputs/enclone_test{}_output; \
                 git add test/inputs/outputs/enclone_test{}_output\n",
                pre_arg,
                test,
                it + 1,
                it + 1
            );
            emit_end_escape(&mut log);
            res.2 = stringme(&log);
        } else {
            let mut old = String::new();
            if !expect_fail && !expect_ok {
                old = read_to_string(&out_file).unwrap();
            }
            let args = parse_bsv(&test);

            // Form the command and execute it.

            let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
            let mut new = new.arg(&args[0]);
            if !no_pre {
                new = new.arg(&pre_arg);
            }
            for i in 1..args.len() {
                new = new.arg(&args[i]);
            }
            // dubious use of expect:
            let new = new
                .arg("FORCE_EXTERNAL")
                // Cap number of cores at 24.  Surprisingly, for testing on a 64-core
                // server, this significantly reduces wallclock.  And substituting either
                // 16 or 64 is slower.  Slower at the time of testing!  As we add tests or
                // change the algorithms, this may change.
                .arg("MAX_CORES=24")
                .output()
                .expect(&format!("failed to execute enclone for test{}", it + 1));
            let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
            let new2 = stringme(&new.stdout);

            // Process tests that were supposed to fail or supposed to succeed.

            if expect_fail || expect_ok {
                res.1 = false;
                if new.status.code().is_none() {
                    fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                    fwriteln!(
                        log,
                        "Something really funky happened, status code unavailable.\n"
                    );
                } else {
                    let status = new.status.code().unwrap();
                    if expect_fail {
                        if status == 0 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have failed, but instead \
                                 succeeded.\n"
                            );
                        } else if status != 1 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have failed with exit status 1,\n\
                                 but instead failed with exit status {}.\n",
                                status
                            );
                        } else {
                            res.1 = true;
                        }
                    } else {
                        if status != 0 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have succeeded, but instead \
                                 failed.\n"
                            );
                        } else {
                            res.1 = true;
                        }
                    }
                }
                res.2 = stringme(&log);

            // Process tests that yield the expected stdout.
            } else if old == new2 {
                res.1 = true;
                if old.len() <= 1 && !expect_null {
                    fwriteln!(
                        log,
                        "\nWarning: old output for subtest {} has {} bytes.\n",
                        it + 1,
                        old.len()
                    );
                }
                if new.stderr.len() > 0 {
                    fwriteln!(log, "Command for subtest {} failed.\n", it + 1);
                    fwriteln!(log, "stderr has {} bytes:", new.stderr.len());
                    fwrite!(log, "{}", strme(&new.stderr));
                    res.1 = false;
                }
                res.2 = stringme(&log);

            // Process tests that yield unexpected stdout.
            } else {
                fwriteln!(log, "\nSubtest {}: old and new differ", it + 1);
                fwriteln!(
                    log,
                    "old has u8 length {} and new has u8 length {}",
                    old.len(),
                    new2.len()
                );
                let mut oldc = Vec::<char>::new();
                let mut newc = Vec::<char>::new();
                for c in old.chars() {
                    oldc.push(c);
                }
                for c in new2.chars() {
                    newc.push(c);
                }
                fwriteln!(
                    log,
                    "old has char length {} and new has char length {}",
                    oldc.len(),
                    newc.len()
                );
                for i in 0..min(oldc.len(), newc.len()) {
                    if oldc[i] != newc[i] {
                        fwriteln!(
                            log,
                            "the first difference is at character {}: old = \"{}\", \
                             new = \"{}\"\n",
                            i,
                            oldc[i],
                            newc[i]
                        );
                        break;
                    }
                }
                fwrite!(log, "old:\n{}", old);
                fwrite!(log, "new:\n{}", new2);
                if new_err.len() != 1 || new_err[0].len() != 0 {
                    fwriteln!(log, "stderr has {} lines:", new_err.len());
                    for i in 0..new_err.len() {
                        fwriteln!(log, "{}", new_err[i]);
                    }
                }
                // let f = format!(
                //     "test/inputs/version{}/{}/outs/all_contig_annotations.json.lz4",
                //         version, args[0].after("=") );
                // if !path_exists(&f) {
                //     println!( "Perhaps you forgot to lz4 compress the json file.\n" );
                //     std::process::exit(1);
                // }
                // println!( "The size of {} is {} bytes.", f, fs::metadata(&f).unwrap().len() );

                fwriteln!(
                    log,
                    "enclone subtest {} failed.  If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     cellranger/lib/rust/enclone (essential!):\n",
                    it + 1
                );
                if set_in_stone {
                    fwriteln!(
                        log,
                        "🔴 However, the output of this test was not supposed to have changed.\n\
                         🔴 Please be extremely careful if you change it.\n",
                    );
                }
                emit_bold_escape(&mut log);
                fwriteln!(
                    log,
                    "enclone {} {} \
                     > test/inputs/outputs/enclone_test{}_output\n",
                    pre_arg,
                    test,
                    it + 1
                );
                emit_end_escape(&mut log);
                fwrite!(log, "and then committing the changed file.  ");
                fwriteln!(
                    log,
                    "You can then retest using:\n\n\
                     cargo test -p enclone enclone  -- --nocapture"
                );
                if new2.len() > 0 {
                    fwriteln!(log, "");
                    res.2 = stringme(&log);
                } else if old != new2 {
                    fwriteln!(log, "old != new");
                    res.2 = stringme(&log);
                }
            }
        }
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

// Test site for broken links and spellcheck.
//
// Two approaches for checking broken links left in place for now, to delete one, and the
// corresponding crate from Cargo.toml.
//
// This looks for
// ▓<a href="..."▓
// ▓http:...[, '")}<#\n]▓
// ▓https:...[, '")}<#\n]▓.
// (These also test termination by ". ".)
// Should also look for at least:
// ▓ href="..."▓
// ▓ href='...'▓
// ▓ src="..."▓.

#[cfg(not(feature = "basic"))]
#[test]
fn test_for_broken_links_and_spellcheck() {
    extern crate attohttpc;
    use std::time::Duration;

    // Set up dictionary exceptions.

    let extra_words = "barcode barcoding clonotype clonotypes clonotyping codebase contig contigs \
        csv cvars enclone executables genomics germline github grok hypermutation hypermutations \
        indel indels linux loh lvars metadata onesie parseable pbmc spacebar subclonotype \
        subclonotypes svg umi umis underperforming vdj website workflow zenodo";
    let extra_words = extra_words.split(' ').collect::<Vec<&str>>();

    // Set up dictionary.

    let dictionary0 = read_to_string("src/english_wordlist").unwrap();
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

    let mut htmls = vec!["index.html".to_string()];
    let pages = read_dir("pages").unwrap();
    for page in pages {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(format!("{}", page));
        }
    }
    let auto = read_dir("pages/auto").unwrap();
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

    for x in htmls {
        let f = open_for_read![x];
        let depth = x.matches('/').count();
        for line in f.lines() {
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
                if !bin_member(&dictionary, &wl.to_string()) {
                    eprintln!(
                        "\nthe word \"{}\" in file {} isn't in the dictionary\n",
                        w, x
                    );
                    std::process::exit(1);
                }
            }

            // Check links.

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
                    for _ in 0..depth {
                        if !z.starts_with("../") {
                            eprintln!("something wrong with file {} on page {}", link, x);
                            std::process::exit(1);
                        }
                        z = z.after("../").to_string();
                    }
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
            for link in links {
                eprintln!("checking link \"{}\"", link);

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
                        eprintln!("\ncould not read link {} on page {}\n", link, x);
                        if i == LINK_RETRIES - 1 {
                            std::process::exit(1);
                        }
                    } else {
                        let response = response.unwrap();
                        if response.is_success() {
                            break;
                        }
                        eprintln!("\ncould not read link {} on page {}\n", link, x);
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
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test site examples to make sure they are what they claim to be, and that the
// merged html files are correct.

#[cfg(not(feature = "basic"))]
#[test]
fn test_site_examples() {
    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let test = SITE_EXAMPLES[i].1;
        let in_stuff = read_to_string(&format!("pages/auto/{}.html", example_name)).unwrap();
        let args = parse_bsv(&test);
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&args)
            .output()
            .expect(&format!("failed to execute test_site_examples"));
        let out_stuff = stringme(&new.stdout);
        if in_stuff != out_stuff {
            eprintln!("\nThe output for site example {} has changed.\n", i + 1);
            eprintln!(
                "Possibly this could be because you're running \"cargo t\" in an \
                environment without the\n\
                extended dataset collection.  Possibly you should run \
                \"cargo test basic -- --nocapture\" instead.\n"
            );
            std::process::exit(1);
        }
    }

    insert_html("pages/index.html.src", "test/outputs/index.html");
    insert_html("pages/expanded.html.src", "test/outputs/expanded.html");

    if read_to_string("index.html").unwrap() != read_to_string("test/outputs/index.html").unwrap() {
        eprintln!("\nContent of index.html has changed.\n");
        std::process::exit(1);
    }
    if read_to_string("pages/auto/expanded.html").unwrap()
        != read_to_string("test/outputs/expanded.html").unwrap()
    {
        eprintln!("\nContent of expanded.html has changed.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that examples are what we claim they are.

#[test]
fn test_enclone_examples() {
    PrettyTrace::new().on();
    for t in 0..EXAMPLES.len() {
        let testn = format!("{}", EXAMPLES[t]);
        let out_file = format!("src/example{}", t + 1);
        let old = read_to_string(&out_file).unwrap();
        let args = testn.split(' ').collect::<Vec<&str>>();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
        for i in 0..args.len() {
            new = new.arg(&args[i]);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_enclone_examples"));
        let new2 = stringme(&new.stdout);
        if old != new2 {
            eprintln!(
                "\nenclone_test_examples: the file example{} is not up to date\n",
                t + 1
            );
            eprintln!("old output =\n{}", old);
            eprintln!("new output =\n{}\n", new2);
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that references to the dataset version in README.md are current.

#[test]
fn test_version_number_in_readme() {
    PrettyTrace::new().on();
    let readme = read_to_string("README.md").unwrap();
    let fields = readme.split('/').collect::<Vec<&str>>();
    for x in fields {
        if x.starts_with("version") {
            let y = x.after("version");
            if y.parse::<usize>().is_ok() {
                let v = y.force_usize();
                assert_eq!(v, TEST_FILES_VERSION as usize);
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that the DejaVuSansMono definition in enclone.css has not changed.  We put this here
// because that definition has to be manually tested, and we don't want it accidentally changed
// and broken.  This is really gross, but it's not clear how to do it better.
//
// Absolutely hideous implementation to verify that
// cat pages/enclone.css | head -36 = "2474276863 1467".
//
// Only works with high probability.

#[test]
fn test_dejavu() {
    PrettyTrace::new().on();
    let mut cat_output_child = Command::new("cat")
        .arg("pages/enclone.css")
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    if let Some(cat_output) = cat_output_child.stdout.take() {
        let mut head_output_child = Command::new("head")
            .arg("-36")
            .stdin(cat_output)
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        cat_output_child.wait().unwrap();
        if let Some(head_output) = head_output_child.stdout.take() {
            let cksum_output_child = Command::new("cksum")
                .stdin(head_output)
                .stdout(Stdio::piped())
                .spawn()
                .unwrap();
            let cksum_stdout = cksum_output_child.wait_with_output().unwrap();
            head_output_child.wait().unwrap();
            let cksum = String::from_utf8(cksum_stdout.stdout).unwrap();
            println!("cksum = {}", cksum);
            assert!(cksum == "2474276863 1467\n".to_string());
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that enclone help all HTML works (without STABLE_DOC).

#[test]
fn test_help_no_stable() {
    PrettyTrace::new().on();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg("help");
    new = new.arg("all");
    new = new.arg("HTML");
    new = new.arg("NOPAGER");
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute test_help_output"));
    if new.status.code() != Some(0) {
        eprintln!("Attempt to run enclone help all without STABLE_DOC failed.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that help output hasn't changed.

#[test]
fn test_help_output() {
    PrettyTrace::new().on();
    let pages = vec![
        "setup",
        "main",
        "quick",
        "how",
        "command",
        "glossary",
        "example1",
        "example2",
        "support",
        "input",
        "input_tech",
        "parseable",
        "plot",
        "filter",
        "special",
        "lvars",
        "cvars",
        "amino",
        "display",
        "indels",
        "color",
        "ideas",
        "faq",
        "developer",
        "all",
    ];
    for p in pages {
        let mut command = format!("enclone help {}", p);
        if p == "setup" {
            command = "enclone help".to_string();
        } else if p == "main" {
            command = "enclone".to_string();
        }
        let out_file = format!("pages/auto/help.{}.html", p);
        let old = read_to_string(&out_file).unwrap();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg("HTML");
        if p == "setup" {
            new = new.arg("help");
        } else if p == "main" {
        } else {
            new = new.arg("help");
            new = new.arg(p);
        }
        new = new.arg("STABLE_DOC");
        new = new.arg("NOPAGER");
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_help_output"));
        if new.status.code() != Some(0) {
            eprintln!("Attempt to run {} failed.\n", command);
            std::process::exit(1);
        }
        let new2 = stringme(&new.stdout);
        if old != new2 {
            eprintme!(old.len(), new2.len());
            eprintln!(
                "\nHelp test failed on {}.\n\
                 You need to update help output by typing \"./build_help\", \
                    assuming that the change is expected.\n",
                p
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that PREBUILD works.

#[test]
fn test_enclone_prebuild() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mb = format!(
        "test/inputs/version{}/123749/outs/raw_feature_bc_matrix/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    );
    if path_exists(&mb) {
        remove_file(&mb).unwrap();
    }

    // First pass: run with NH5.

    let test_id = 48;
    let it = test_id - 1;
    let testn = format!("{} NH5", TESTS[it]);
    let out_file = format!("test/inputs/outputs/enclone_test{}_output", test_id);
    let old = read_to_string(&out_file).unwrap();
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
    for i in 0..args.len() {
        new = new.arg(&args[i]);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute test_enclone_prebuild"));
    // let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nenclone_test_prebuild: first pass output has changed.\n\
             You may want to add more info to this failure message.\n"
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }
    if !path_exists(&format!(
        "test/inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!("\nenclone_test_prebuild: did not create feature_barcode_matrix.bin.");
    }

    // Second pass: run without PREBUILD but using the feature_barcode_matrix.bin that the first
    // pass created.

    let testn = TESTS[it];
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
    for i in 0..args.len() {
        new = new.arg(&args[i]);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute enclone_test_prebuild"));
    // let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nenclone_test_prebuild: second pass output has changed.\n\
             You may want to add more info to this failure message.\n\
             And don't forget to remove feature_barcode_matrix.bin.\n"
        );
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }

    // Clean up: delete feature_barcode_matrix.bin.

    std::fs::remove_file(&format!(
        "test/inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    ))
    .unwrap();
    println!("\nused {:.2} seconds in enclone_test_prebuild", elapsed(&t));
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This test runs enclone for two test inputs, with LOUPE output
// turned on. It will then read both the bincode and proto file created
// and asserts that we get the same data structure either way.
//
// It also tests to make sure that the LOUPE output is unchanged.  If it changed for a good
// reason, update the output file.  Otherwise perhaps something has gone wrong!

#[test]
fn test_proto_write() -> Result<(), Error> {
    let tests = vec!["BCR=123085", "TCR=101287"];
    let pre_arg = format!("PRE=test/inputs/version{}", TEST_FILES_VERSION);
    let binary_arg = format!("BINARY={}.bin", LOUPE_OUT_FILENAME);
    let proto_arg = format!("PROTO={}.proto", LOUPE_OUT_FILENAME);
    for t in tests.iter() {
        // FIXME: It would be nicer to use the enclone API here
        std::process::Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&[&pre_arg, *t, &binary_arg, &proto_arg])
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));

        // Test to make sure output is unchanged.

        let oldx = format!("test/inputs/{}.binary.gz", t.after("="));
        let newx = format!("{}.bin", LOUPE_OUT_FILENAME);
        let mut f = File::open(&oldx)?;
        let mut oldbufgz = Vec::<u8>::new();
        f.read_to_end(&mut oldbufgz)?;
        let mut gz = GzDecoder::new(&oldbufgz[..]);
        let mut oldbuf = Vec::<u8>::new();
        gz.read_to_end(&mut oldbuf)?;
        let mut f = File::open(&newx)?;
        let mut newbuf = Vec::<u8>::new();
        f.read_to_end(&mut newbuf)?;
        if oldbuf != newbuf {
            eprintln!(
                "\nThe binary output of enclone on {} has changed.  If this is expected,\n\
                please regenerate the file {}\n",
                t, oldx
            );
            std::process::exit(1);
        }

        // Test to make sure proto and bin are consistent.

        let outputs_proto = read_proto(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(format!("{}.bin", LOUPE_OUT_FILENAME));
        std::fs::remove_file(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        std::fs::remove_file(format!("{}.bin", LOUPE_OUT_FILENAME))?;
        assert!(outputs_proto == outputs_bin);
    }

    Ok(())
}
