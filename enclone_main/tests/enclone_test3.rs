// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

use ansi_escape::*;
use enclone::html::*;
use enclone_core::defs::*;
use enclone_core::run_test::*;
use enclone_core::testlist::*;
use enclone_core::*;
use enclone_proto::proto_io::{read_proto, ClonotypeIter};
use enclone_proto::types::EncloneOutputs;
use failure::Error;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use sha2::Digest;
use std::cmp::min;
use std::env;
use std::fs::{read_to_string, remove_dir_all, remove_file, File};
use std::io::{BufWriter, Read, Write};
use std::process::{Command, Stdio};
use std::time::Instant;
use string_utils::*;

const LOUPE_OUT_FILENAME: &str = "testx/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 15.

// NOT BASIC

// Test site examples to make sure they are what they claim to be, and that the
// merged html files are correct.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_site_examples() {
    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let test = SITE_EXAMPLES[i].1;
        let in_file = format!("../{}", example_name);
        let in_stuff = read_to_string(&in_file).expect(&format!("couldn't find {}", in_file));
        let args = parse_bsv(&test);
        if args.contains(&"GEX=123217".to_string()) && !args.contains(&"H5".to_string()) {
            panic!("Oops please fix this, to prevent sporadic failures.");
        }
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&args)
            .output()
            .expect(&format!("failed to execute test_site_examples"));
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_site_examples: example {} failed to execute, stderr =\n{}",
                i + 1,
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        let out_stuff = strme(&new.stdout);
        if in_stuff != out_stuff {
            eprintln!("\nThe output for site example {} has changed.\n", i + 1);
            eprintln!("stderr:\n{}", strme(&new.stderr));
            let old_lines = in_stuff.split('\n').collect::<Vec<&str>>();
            let new_lines = out_stuff.split('\n').collect::<Vec<&str>>();
            eprintln!(
                "old stdout has {} lines; new stdout has {} lines",
                old_lines.len(),
                new_lines.len(),
            );
            for i in 0..min(old_lines.len(), new_lines.len()) {
                if old_lines[i] != new_lines[i] {
                    eprintln!(
                        "first different stdout line is line {}\nold = {}\nnew = {}",
                        i + 1,
                        old_lines[i],
                        new_lines[i],
                    );
                    break;
                }
            }
            let save = format!("testx/outputs/{}", example_name.rev_after("/"));
            {
                let mut f = open_for_write_new![&save];
                fwrite!(f, "{}", out_stuff);
            }
            let mut in_filex = in_file.clone();
            if in_filex.starts_with("../") {
                in_filex = in_filex.after("../").to_string();
            }
            eprintln!("\nPlease diff {} enclone_main/{}.", in_filex, save);
            eprintln!(
                "\nPossibly this could be because you're running \"cargo t\" in an \
                environment without the\n\
                extended dataset collection.  Possibly you should run \
                \"cd enclone; cargo test basic -- --nocapture\" instead.\n\n\
                Otherwise, if you're satisfied with the new output, you can update using\n\n\
                enclone {} > {}.\n",
                args.iter().format(" "),
                example_name
            );
            std::process::exit(1);
        }
    }

    insert_html(
        "../pages/index.html.src",
        "testx/outputs/index.html",
        true,
        0,
    );
    insert_html(
        "../pages/expanded.html.src",
        "testx/outputs/expanded.html",
        true,
        2,
    );
    let new_index = read_to_string("testx/outputs/index.html").unwrap();
    if read_to_string("../index.html").unwrap() != new_index {
        eprintln!("\nContent of index.html has changed.");
        {
            let mut f = open_for_write_new!["testx/outputs/index.html.new"];
            fwrite!(f, "{}", new_index);
        }
        eprintln!("Please diff index.html enclone_main/testx/outputs/index.html.new.\n");
        std::process::exit(1);
    }
    /*
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != edit_html(&read_to_string("testx/outputs/expanded.html").unwrap())
    */
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != read_to_string("testx/outputs/expanded.html").unwrap()
    {
        eprintln!("\nContent of expanded.html has changed.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 16. Test that examples are what we claim they are.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_examples() {
    PrettyTrace::new().on();
    for t in 0..EXAMPLES.len() {
        let testn = format!("{}", EXAMPLES[t]);
        let out_file = format!("../enclone_help/src/example{}", t + 1);
        let old = read_to_string(&out_file).unwrap();
        let args = testn.split(' ').collect::<Vec<&str>>();
        if args.contains(&"GEX=123217") && !args.contains(&"H5") {
            panic!("Oops please fix this, to prevent sporadic failures.");
        }
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg(format!(
            "PRE=../enclone-data/big_inputs/version{}",
            TEST_FILES_VERSION
        ));
        for i in 0..args.len() {
            new = new.arg(&args[i]);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_enclone_examples"));
        let new2 = stringme(&new.stdout);
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_test_examples: example{} failed to execute, stderr =\n{}",
                t + 1,
                strme(&new.stderr),
            );
            eprintln!("If it's not clear what is happening, make sure you've run ./build.\n");
            std::process::exit(1);
        }
        if old != new2 {
            eprintln!(
                "\nenclone_test_examples: the file enclone_help/src/example{} is not up to date\n",
                t + 1
            );
            eprintln!("old output =\n{}", old);
            eprintln!("new output =\n{}\n", new2);
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 17. Test that references to the dataset version in README.md are current.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_version_number_in_readme() {
    PrettyTrace::new().on();
    let readme = read_to_string("../README.md").unwrap();
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

// 18. Test that the DejaVuSansMono definition in enclone_css_v2.css has not changed.  We put this
// here because that definition has to be manually tested, and we don't want it accidentally
// changed and broken.  This is really gross, but it's not clear how to do it better.
//
// Absolutely hideous implementation to verify that
// cat ../pages/enclone_css_v2.css | head -36 = "2474276863 1467".
//
// Only works with high probability.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_dejavu() {
    PrettyTrace::new().on();
    let mut cat_output_child = Command::new("cat")
        .arg("../pages/enclone_css_v2.css")
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
            // println!("cksum = {}", cksum);
            assert!(cksum == "2474276863 1467\n".to_string());
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 19. Test that help output hasn't changed.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_output() {
    PrettyTrace::new().on();
    for p in HELP_PAGES.iter() {
        let mut command = format!("enclone help {}", p);
        if *p == "setup" {
            command = "enclone help".to_string();
        } else if *p == "main" {
            command = "enclone".to_string();
        }
        let out_file = format!("../pages/auto/help.{}.html", p);
        let old = read_to_string(&out_file).unwrap();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg("HTML");
        if *p == "setup" {
            new = new.arg("help");
        } else if *p == "main" {
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
        let new2 = edit_html(&stringme(&new.stdout));
        if old != new2 {
            eprintme!(old.len(), new2.len());
            eprintln!(
                "\nHelp test failed on {}.\n\
                 You need to update help output by typing \"./build\", \
                    assuming that the change is expected.\n",
                p
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 20. Test that enclone help all HTML works (without STABLE_DOC).

#[cfg(not(feature = "cpu"))]
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

// 21. Test that PREBUILD works.  Because this creates and then deletes the .bin file, it
// would play havoc with any test that runs with GEX=123217, unless it also has H5, resulting
// in sporadic (rare) test failures.  So don't do that.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_prebuild() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mb = format!(
        "../enclone-data/big_inputs/version{}/123217/outs/raw_feature_bc_matrix/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    );
    if path_exists(&mb) {
        remove_file(&mb).unwrap();
    }

    // First pass: run with NH5.

    let test_id = 48;
    let it = test_id - 1;
    let testn = format!("{} NH5", TESTS[it]);
    let out_file = format!("testx/inputs/outputs/enclone_test{}_output", test_id);
    let old = read_to_string(&out_file).unwrap();
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
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
            If you are OK with the new output, it should work to type:\n\
            enclone {} > enclone_main/{}\n",
            testn, out_file
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        eprintln!("new stderr = \n{}\n", strme(&new.stderr));
        std::process::exit(1);
    }
    if !path_exists(&format!(
        "../enclone-data/big_inputs/version{}/123217/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!("\nenclone_test_prebuild: did not create feature_barcode_matrix.bin.");
    }

    // Second pass: run without PREBUILD but using the feature_barcode_matrix.bin that the first
    // pass created.

    let testn = TESTS[it];
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
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
        "../enclone-data/big_inputs/version{}/123217/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    ))
    .unwrap();
    println!("\nused {:.2} seconds in enclone_test_prebuild", elapsed(&t));
}

fn check_enclone_outs_consistency(enclone_outs: &EncloneOutputs) {
    let uref_items = &enclone_outs.universal_reference.items;
    for cl in &enclone_outs.clonotypes {
        for cl_chain in &cl.chains {
            assert!(uref_items.len() > cl_chain.v_idx as usize);
            assert!(uref_items.len() > cl_chain.j_idx as usize);
            assert_eq!(uref_items[cl_chain.v_idx as usize].region, 1);
            assert_eq!(uref_items[cl_chain.j_idx as usize].region, 3);
            if let Some(d_idx) = cl_chain.d_idx {
                assert!(uref_items.len() > d_idx as usize);
                assert_eq!(uref_items[d_idx as usize].region, 2);
            }
            if let Some(u_idx) = cl_chain.u_idx {
                assert!(uref_items.len() > u_idx as usize);
                assert_eq!(uref_items[u_idx as usize].region, 0);
            }
            for ex_cl in &cl.exact_clonotypes {
                for ex_cl_chain in ex_cl.chains.iter().map(|info| &info.chain) {
                    if let Some(c_region_idx) = ex_cl_chain.c_region_idx {
                        assert!(uref_items.len() > c_region_idx as usize);
                        assert_eq!(uref_items[c_region_idx as usize].region, 4);
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 22. This test runs enclone for two test inputs, with LOUPE output
// turned on. It will then read both the bincode and proto file created
// and asserts that we get the same data structure either way.
//
// It also tests to make sure that the LOUPE output is unchanged.  If it changed for a good
// reason, update the output file.  Otherwise perhaps something has gone wrong!

#[cfg(not(feature = "cpu"))]
#[test]
fn test_proto_write() -> Result<(), Error> {
    let tests = vec!["BCR=123085", "TCR=101287"];
    let pre_arg = format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    );

    let bin_file = format!("{}.binary", LOUPE_OUT_FILENAME);
    let proto_file = format!("{}.proto", LOUPE_OUT_FILENAME);

    let binary_arg = format!("BINARY={}", bin_file);
    let proto_arg = format!("PROTO={}", proto_file);
    for t in tests.iter() {
        // FIXME: It would be nicer to use the enclone API here
        std::process::Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&[&pre_arg, *t, &binary_arg, &proto_arg])
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));

        // Test to make sure proto and bin are consistent.

        let outputs_proto = read_proto(&proto_file)?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(&bin_file);
        if outputs_proto != outputs_bin {
            eprintln!("\noutputs_proto is not equal to outputs_bin\n");
            std::process::exit(1);
        }

        // Test to make sure that the clonotype iterator works

        let clonotypes: Vec<_> = ClonotypeIter::from_file(&proto_file).unwrap().collect();
        assert!(clonotypes == outputs_proto.clonotypes);

        // Check consistency

        check_enclone_outs_consistency(&outputs_proto);

        // Test to make sure output is unchanged.

        let oldx = format!("testx/inputs/{}.binary.sha256", t.after("="));
        let mut fold = std::fs::File::open(&oldx)?;
        let mut cksum_old = Vec::<u8>::new();
        fold.read_to_end(&mut cksum_old)?;
        let mut fnew = std::fs::File::open(&bin_file)?;
        let mut cksum_new = Vec::<u8>::new();
        fnew.read_to_end(&mut cksum_new)?;
        let cksum_new = format!("{:x}", sha2::Sha256::digest(&cksum_new));
        std::fs::remove_file(&proto_file)?;
        std::fs::remove_file(&bin_file)?;
        if strme(&cksum_old) != cksum_new {
            eprintln!(
                "\nThe binary output of enclone on {} has changed.  If this is expected,\n\
                please run the command\n\
                echo -n {} > enclone_main/testx/inputs/{}.binary.sha256",
                t,
                &cksum_new,
                t.after("=")
            );
            eprintln!(
                "\nIf you can't figure out what happened, first turn off this test and see\n\
                if any other tests failed.  If not, try the following:\n\n\
                1. build old code\n\
                2. enclone BCR=123085 NOPRINT PROTO=~/old\n\
                3. <build new code>\n\
                4. enclone BCR=123085 NOPRINT PROTO=~/new\n\
                5. cmp ~/old ~/new\n\
                6. in the quasi-readable output, locate the barcode at the first point after the\n\
                point where the two files differ (or conceivably before)\n\
                7. analyze what is happening with that barcode\n"
            );
            std::process::exit(1);
        }
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 23. Test for no change to the output of the command that appears in the enclone annotated
// example on the landing page.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_annotated_example() {
    PrettyTrace::new().on();
    let it = 0;
    let test = "BCR=123085 CDR3=CTRDRDLRGATDAFDIW";
    let mut out = String::new();
    let mut log = String::new();
    let mut ok = false;
    run_test(
        env!("CARGO_BIN_EXE_enclone"),
        it,
        "",
        &test,
        "annotated_example_test",
        &mut ok,
        &mut log,
        &mut out,
    );
    print!("{}", log);
    if !ok {
        let mut log = Vec::<u8>::new();
        emit_red_escape(&mut log);
        fwriteln!(
            log,
            "Oh no: the results for the annotated example on the landing page have \
            changed.  Assuming that\nthe change is intentional, to fix this you \
            need to do two things:\n\
            1. Update the test output.\n\
            2. Manually update the annotated example output.\n\
            Because the second item is such a big pain, we stopped doing it, but you should\n\
            check to make sure that it has not changed too much from what is shown."
        );
        emit_end_escape(&mut log);
        eprintln!("{}", strme(&log));
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 24. Test SUBSET_JSON option.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_subset_json() {
    // Note need create_dir_all because testx/outputs may not exist for GitHub Actions.
    std::fs::create_dir_all("testx/outputs/woof").unwrap();
    let test = r###"BCR=123085,123089 CDR3=CARVGSFLSSSWHPRDYYYYGMDVW LVARS=n SUBSET_JSON=testx/outputs/woof/all_contig_annotations.json"###;
    let pre_arg = format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    );
    let args = parse_bsv(&test);
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .arg(&pre_arg)
        .args(&args)
        .output()
        .expect(&format!("failed to execute test_subset_json 1"));
    if new.status.code() != Some(0) {
        eprint!(
            "\nsubset json test 1: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let o1 = new.stdout;
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .arg("BCR=testx/outputs/woof")
        .output()
        .expect(&format!("failed to execute test_subset_json 2"));
    if new.status.code() != Some(0) {
        eprint!(
            "\nsubset json test 2: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let o2 = new.stdout;
    if o1 != o2 {
        eprintln!("\nSubset json test failed: outputs are unequal.\n");
        eprintln!("output 1:\n{}\n", strme(&o1));
        eprintln!("output 2:\n{}\n", strme(&o2));
        std::process::exit(1);
    }
    let _ = remove_dir_all("testx/outputs/woof");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 25. Test that ordinary cell-exact variables have a _cell version.  By ordinary we mean
// that the variable is described in variables.html.src by a name in the alphabet a-zA-Z0-9_,
// although that's not exactly how we test.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_cell_exact() {
    let vars = include_str!["../../pages/variables.html.src"];
    let mut lines = Vec::<String>::new();
    for line in vars.lines() {
        lines.push(line.to_string());
    }
    let mut fails = Vec::<String>::new();
    for i in 0..lines.len() {
        // WARNING: NOTE THAT FOOTNOTE 3 TEST IS HARDCODED.
        if lines[i].contains("<td> cell-exact") && !lines[i - 2].contains(" 3 ") {
            let var = lines[i - 3].between("<code> ", " </code>");
            if !var.contains("&") {
                let lvar = lines[i - 1].contains("lvar");
                let varp;
                if lvar {
                    varp = format!("{}_cell", var);
                } else {
                    varp = format!("{}_cell1", var);
                }
                let pre_arg = format!(
                    "PRE=../enclone-data/big_inputs/version{}",
                    TEST_FILES_VERSION
                );
                let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
                    .arg(&pre_arg)
                    .arg("BCR=86237")
                    .arg("GEX=85679")
                    .arg("NH5")
                    .arg("POUT=stdout")
                    .arg(&format!("PCOLS={}", varp))
                    .arg("PCELL")
                    .output()
                    .expect(&format!("failed to execute test_cel_exact"));
                if new.status.code() != Some(0) {
                    fails.push(var.to_string());
                }
            }
        }
    }
    if !fails.is_empty() {
        let mut msg = String::new();
        for i in 0..fails.len() {
            msg += &format!(
                "{} is cell-exact but does not have a _cell version\n",
                fails[i]
            );
        }
        eprintln!("\n{}", msg);
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 26. Test running with just reference.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_ref_only() {
    let test = "REF=../enclone-data/big_inputs/version15/█≈ΠΠΠ≈█/outs/\
        vdj_reference/fasta/regions.fa";
    let args = parse_bsv(&test);
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .args(&args)
        .output()
        .expect(&format!("failed to execute test_subset_json 1"));
    if new.status.code() == Some(0) {
        eprint!("\ntest_ref_only: enclone command should not have succeeded.\n");
        std::process::exit(1);
    }
    if !strme(&new.stderr).contains("No TCR or BCR data have been specified.") {
        eprintln!("\ntest_ref_only: unexpected error message\n");
        std::process::exit(1);
    }
}
