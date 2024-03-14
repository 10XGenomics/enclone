// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

use ansi_escape::*;
use anyhow::Error;
use enclone_core::defs::*;
use enclone_core::*;
use enclone_proto::proto_io::{read_proto, ClonotypeIter};
use enclone_proto::types::EncloneOutputs;
use enclone_testlist::main_testlist::*;
use enclone_testlist::*;
use enclone_tools::html::*;
use enclone_tools::run_test::*;
use io_utils::*;
use itertools::Itertools;

use rayon::prelude::*;
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
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..SITE_EXAMPLES.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let example_name = SITE_EXAMPLES[i].0;
        let test = SITE_EXAMPLES[i].1;
        let in_file = format!("../{}", example_name);
        let mut f = File::open(&in_file).expect(&format!("couldn't find {}", in_file));
        let mut in_stuff = Vec::<u8>::new();
        f.read_to_end(&mut in_stuff).unwrap();
        let mut args = parse_bsv(&test);
        args.push("NO_KILL");
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&args)
            .output()
            .expect(&format!("failed to execute test_site_examples"));
        if new.status.code() != Some(0) {
            res.2 = format!(
                "\nenclone_site_examples: example {} failed to execute, stderr =\n{}",
                i + 1,
                strme(&new.stderr),
            );
            res.1 = true;
        }
        let out_stuff = new.stdout.to_vec();
        if in_stuff != out_stuff {
            res.1 = true;
            res.2 += &mut format!("\nThe output for site example {} has changed.\n\n", i + 1);
            res.2 += &mut format!("stderr:\n{}\n", strme(&new.stderr));
            if String::from_utf8(in_stuff.clone()).is_err()
                && String::from_utf8(out_stuff.clone()).is_err()
            {
                res.2 += &mut format!("Both the old and new files are binary.\n\n");
            } else if String::from_utf8(in_stuff.clone()).is_err() {
                res.2 += &mut format!("The new file is binary but the old is not.\n\n");
            } else if String::from_utf8(out_stuff.clone()).is_err() {
                res.2 += &mut format!("The old file is binary but the new is not.\n\n");
            } else {
                let old_lines = strme(&in_stuff).split('\n').collect::<Vec<&str>>();
                let new_lines = strme(&out_stuff).split('\n').collect::<Vec<&str>>();
                res.2 += &mut format!(
                    "old stdout has {} lines; new stdout has {} lines\n",
                    old_lines.len(),
                    new_lines.len(),
                );
                for i in 0..min(old_lines.len(), new_lines.len()) {
                    if old_lines[i] != new_lines[i] {
                        res.2 += &mut format!(
                            "first different stdout line is line {}\nold = {}\nnew = {}\n",
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
                    fwrite!(f, "{}", strme(&out_stuff));
                }
                let mut in_filex = in_file.clone();
                if in_filex.starts_with("../") {
                    in_filex = in_filex.after("../").to_string();
                }
                res.2 += &mut format!("\nPlease diff {} enclone_exec/{}.\n", in_filex, save);
                res.2 += &mut format!(
                    "\nPossibly this could be because you're running \"cargo t\" in an \
                    environment without the\n\
                    extended dataset collection.  Possibly you should run \
                    \"cd enclone; cargo test basic -- --nocapture\" instead.\n\n\
                    Otherwise, if you're satisfied with the new output, you can update using\n\n\
                    enclone {} > {}.\n\n",
                    args.iter().format(" "),
                    example_name
                );
            }
        }
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if results[i].1 {
            panic!("failed");
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
        eprintln!("Please diff index.html enclone_exec/testx/outputs/index.html.new.\n");
        panic!("failed");
    }
    /*
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != edit_html(&read_to_string("testx/outputs/expanded.html").unwrap())
    */
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != read_to_string("testx/outputs/expanded.html").unwrap()
    {
        eprintln!("\nContent of expanded.html has changed.\n");
        panic!("failed");
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 16. Test that examples are what we claim they are.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_examples() {
    let examples = [
        (1, r###"BCR=123089 CDR3=CARRYFGVVADAFDIW"###),
        (
            2,
            r###"BCR=123085 GEX=123217 LVARSP=gex,IGHV2-5_g_μ CDR3=CALMGTYCSGDNCYSWFDPW"###,
        ),
    ];
    for (num, args) in examples {
        let out_file = format!("../enclone_help/src/example{}", num);
        let old = read_to_string(&out_file).unwrap();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg(format!(
            "PRE=../enclone-data/big_inputs/version{}",
            TEST_FILES_VERSION
        ));
        for arg in args.split(' ') {
            new = new.arg(arg);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .arg("NO_KILL")
            .output()
            .expect(&format!("failed to execute test_enclone_examples"));
        let new2 = stringme(&new.stdout);
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_test_examples: example{} failed to execute, stderr =\n{}",
                num,
                strme(&new.stderr),
            );
            eprintln!("If it's not clear what is happening, make sure you've run ./build.\n");
            panic!("failed");
        }
        if old != new2 {
            eprintln!(
                "\nenclone_test_examples: the file enclone_help/src/example{} is not up to date\n",
                num
            );
            eprintln!("old output =\n{}", old);
            eprintln!("new output =\n{}\n", new2);
            panic!("failed");
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 17. Test that references to the dataset version in README.md are current.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_version_number_in_readme() {
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
    for p in HELP_PAGES.iter() {
        let mut command = format!("enclone help {}", p);
        if *p == "setup" {
            command = "enclone help".to_string();
        } else if *p == "main" {
            command = "enclone".to_string();
        }
        let out_file = format!("../pages/auto/help.{}.html", p);
        let old = read_to_string(&out_file).unwrap();
        let mut args = Vec::<String>::new();
        if *p == "setup" {
            args.push("help".to_string());
        } else if *p == "main" {
        } else {
            args.push("help".to_string());
            args.push(p.to_string());
        }
        args.push("HTML".to_string());
        args.push("STABLE_DOC".to_string());
        args.push("NOPAGER".to_string());
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let new = new.args(&args);
        let new = new
            .arg("FORCE_EXTERNAL")
            .arg("NO_KILL")
            .output()
            .expect(&format!("failed to execute test_help_output"));
        if new.status.code() != Some(0) {
            eprintln!("Attempt to run {} failed.\n", command);
            eprintln!("stderr = {}", strme(&new.stderr));
            panic!("failed");
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
            panic!("failed");
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 20. Test that enclone help all HTML works (without STABLE_DOC).

#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_no_stable() {
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg("help");
    new = new.arg("all");
    new = new.arg("HTML");
    new = new.arg("NOPAGER");
    let new = new
        .arg("FORCE_EXTERNAL")
        .arg("NO_KILL")
        .output()
        .expect(&format!("failed to execute test_help_output"));
    if new.status.code() != Some(0) {
        eprintln!("Attempt to run enclone help all without STABLE_DOC failed.\n");
        panic!("failed");
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 21. Test that PREBUILD works.  This test used to exercise creation of the
// binary matrix format, since removed.  It's not entirely clear what purpose this
// test is serving now, beyond yet another regression test.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_prebuild() {
    let t = Instant::now();
    let mb = format!(
        "../enclone-data/big_inputs/version{}/123217/outs/raw_feature_bc_matrix/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    );
    if path_exists(&mb) {
        remove_file(&mb).unwrap();
    }

    let (test_num, comment, args) = TESTS[47];
    assert_eq!(48, test_num);
    let out_file = format!("testx/inputs/outputs/enclone_test{}_output", test_num);
    let old = read_to_string(&out_file).unwrap();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
    for arg in args.split(' ') {
        new = new.arg(arg);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .arg("NO_KILL")
        .output()
        .expect(&format!("failed to execute test_enclone_prebuild"));
    // let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nenclone_test_prebuild: first pass output has changed.\n\
            If you are OK with the new output, it should work to type:\n\
            enclone {} > enclone_exec/{}\n",
            args, out_file
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        eprintln!("new stderr = \n{}\n", strme(&new.stderr));
        panic!("failed");
    }

    // Second pass: run without PREBUILD
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
    for arg in args.split(' ') {
        new = new.arg(arg);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .arg("NO_KILL")
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
        panic!("failed");
    }

    println!(
        "\nused {:.2} seconds in enclone_test_prebuild",
        t.elapsed().as_secs_f64()
    );
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
            .arg("NO_KILL")
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));

        // Test to make sure proto and bin are consistent.

        let outputs_proto = read_proto(&proto_file)?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(&bin_file);
        if outputs_proto != outputs_bin {
            eprintln!("\noutputs_proto is not equal to outputs_bin\n");
            panic!("failed");
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
                echo -n {} > enclone_exec/testx/inputs/{}.binary.sha256",
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
            panic!("failed");
        }
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 23. Test for no change to the output of the command that appears in the enclone annotated
// example on the landing page.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_annotated_example() -> Result<(), String> {
    run_test(
        env!("CARGO_BIN_EXE_enclone"),
        "annotated_example_test",
        0,
        1,
        "",
        "BCR=123085 CDR3=CTRDRDLRGATDAFDIW",
    )
    .map_err(|res| {
        format!(
            "{}\n\nOh no: the results for the annotated example on the landing page have \
        changed. Assuming that the change is intentional, to fix this you \
        need to do two things:\n\
        1. Update the test output.\n\
        2. Manually update the annotated example output.\n\
        Because the second item is such a big pain, we stopped doing it, but you should \
        check to make sure that it has not changed too much from what is shown.",
            res.log
        )
    })?;
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 24. Test SUBSET_JSON option.
//
// This uses two examples, that were broken for different reasons.

// not BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_subset_json() {
    for pass in 1..=2 {
        // Note need create_dir_all because testx/outputs may not exist for GitHub Actions.
        std::fs::create_dir_all("testx/outputs/woof").unwrap();
        let test;
        if pass == 1 {
            test = r###"BCR=123085,123089 CDR3=CARVGSFLSSSWHPRDYYYYGMDVW LVARS=n SUBSET_JSON=testx/outputs/woof/all_contig_annotations.json"###;
        } else {
            test = r###"BCR=140703,140704 BUILT_IN CDR3=CARGGALYSSSASYYYYYYGMDVW LVARS=n SUBSET_JSON=testx/outputs/woof/all_contig_annotations.json NOPRINT POUT=stdout PCOLS=n,nchains"###;
        }
        /*
        let pre_arg = format!(
            "PRE=../enclone-data/big_inputs/version{}",
            TEST_FILES_VERSION
        );
        */
        let args = parse_bsv(&test);
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            // .arg(&pre_arg)
            .args(&args)
            .arg("NO_KILL")
            .output()
            .expect(&format!(
                "failed to execute test_subset_json 1, pass = {}",
                pass
            ));
        if new.status.code() != Some(0) {
            eprint!(
                "\nsubset json test 1, pass = {}: failed to execute, stderr =\n{}",
                pass,
                strme(&new.stderr),
            );
            panic!("failed");
        }
        let o1 = new.stdout;
        let mut args = vec!["BCR=testx/outputs/woof".to_string()];
        if pass == 2 {
            args.push("BUILT_IN".to_string());
            args.push("NOPRINT".to_string());
            args.push("POUT=stdout".to_string());
            args.push("PCOLS=n,nchains".to_string());
        }
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&args)
            .arg("NO_KILL")
            .output()
            .expect(&format!(
                "failed to execute test_subset_json 2, pass = {}",
                pass
            ));
        if new.status.code() != Some(0) {
            eprint!(
                "\nsubset json test 2, pass = {}: failed to execute, stderr =\n{}",
                pass,
                strme(&new.stderr),
            );
            panic!("failed");
        }
        let o2 = new.stdout;
        if o1 != o2 {
            eprintln!(
                "\nSubset json test failed: outputs are unequal, pass = {}.\n",
                pass
            );
            eprintln!("output 1:\n{}\n", strme(&o1));
            eprintln!("output 2:\n{}\n", strme(&o2));
            panic!("failed");
        }
        let _ = remove_dir_all("testx/outputs/woof");
    }
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
                    .arg("POUT=stdout")
                    .arg(&format!("PCOLS={}", varp))
                    .arg("PCELL")
                    .arg("NO_KILL")
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
        panic!("failed");
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
        .arg("NO_KILL")
        .output()
        .expect(&format!("failed to execute test_subset_json 1"));
    if new.status.code() == Some(0) {
        eprint!("\ntest_ref_only: enclone command should not have succeeded.\n");
        panic!("failed");
    }
    if !strme(&new.stderr).contains("No TCR or BCR data have been specified.") {
        eprintln!("\ntest_ref_only: unexpected error message\n");
        panic!("failed");
    }
}
