// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// somewhat misplaced comments:
//
// Run enclone on a test case and verify that the output is identical to what
// was gotten before.  If the output is different, look at it
// and decide if the change is justified, and if so update the output file.
//
// This test only runs if you use cargo test --release.  The test for
// not debug_assertions is a proxy for that.
//
// To test just this test, use:
//
// cargo test --release -p enclone enclone -- --nocapture

use ansi_escape::*;
use enclone::proto_io::read_proto;
use enclone::types::EncloneOutputs;
use failure::Error;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::cmp::min;
use std::fs::read_to_string;
use std::io::Write;
use std::process::Command;
use std::time::Instant;
use string_utils::*;

const TEST_FILES_VERSION: u8 = 14;
const LOUPE_OUT_FILENAME: &str = "test/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[cfg(debug_assertions)]
#[test]
fn test_enclone_fail() {
    println!("\n\"cargo test\" deliberately fails here because without running in release mode,");
    println!("the test in enclone would be too slow.  We could simply elide the test, but");
    println!("then you wouldn't know that you're missing an important test.  If you really want");
    println!("to run all tests except this test in debug(dev) mode, please use");
    println!("\"cargo test --all --exclude enclone\"");
    println!("However please also note that even with the extra test, \"cargo test --release\"");
    println!("will be faster then the above.\n");
    assert!(0 == 1);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define a bunch of tests.

const TESTS: [&str; 34] = [
    // 1. tests variant base after CDR3, parseable output
    "BCR=123089 CDR3=CVRDRQYYFDYW POUT=stdout \
     PCOLS=exact_subclonotype_id,ncells,v_name1,v_name2,nchains,var_indices_aa1,barcodes",
    // 2. tests many donor ref differences, test comp and var and donorn
    "BCR=123089 CDR3=CARRYFGVVADAFDIW CVARSP=comp,var AMINO=cdr3,var,share,donorn",
    // 3. tests motif in CDR3, CHAINS, utot, flipped args in CVARS, on tiny dataset
    "BCR=85333 CDR3=\"CAA.*\" CHAINS=2 CVARS=const,utot",
    // 4. tests gex and antibody, FULL_SEQC, ulen, udiff, on tiny dataset
    "BCR=86237 GEX=85679 LVARSP=gex_med,CD19_ab,CD25_ab,IGLV3-1_g,RPS27_g CELLS=3 FULL_SEQC \
     CVARSP=ulen,udiff",
    // 5. tests TCR and correct grouping of onesies on AGBT Donor 2 dataset
    "TCR=101287 MIN_CELLS=100",
    // 6. tests AMINO=
    "BCR=86237 CELLS=3 AMINO= CVARS=umed,rmed,cdr3_dna",
    // 7. tests SHM deletion
    "BCR=123085 CVARSP=var,clen,cdiff CDR3=CAREPLYYDFWSAYFDYW LVARSP=near,far",
    // 8. this clonotype included a junk chain before we made a change
    "TCR=163911 CDR3=CAPSAGDKIIF AMINO=donor",
    // 9. tests PER_BC
    "BCR=85333 CDR3=CAKGDRTGYSYGGGIFDYW PER_BC",
    // 10. tests multiple datasets and also LVARS=ncells,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    "BCR=\"123085;123089\" CDR3=CVKDRVTGTITELDYW LVARS=ncells,donors,datasets AMINO=share MIX_DONORS",
    // 11. tests META
    "META=test/inputs/meta_test11 CDR3=CARSFFGDTAMVMFQAFDPW LVARSP=donors,gex_med",
    // 12. this added because it got better when a noise filter was added, also tests umax
    "TCR=163914 CDR3=CASSLVQPSTDTQYF CVARSP=umax",
    // 13. this added because it got better when a noise filter was added; also test FASTA
    "TCR=163914 CDR3=CAFRGGSYIPTF FASTA=stdout",
    // 14. this added because it got better when a bug in bads detection was fixed
    "TCR=163914 CDR3=CASRLGGEETQYF",
    // 15. tests insertion and AMINO range
    "BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113",
    // BCR=123085 CDR3=CARHPAPNYGFWSGYYKTDNWFDPW ==> alt example if we need to dump 86233
    // 16. tests number of cells broken out by dataset
    "BCR=123085,123089 LVARS=ncells,n_123085,n_123089 CDR3=CTRDRDLRGATDAFDIW",
    // 17. tests gex with PER_BC and tests n_gex
    // See also enclone_test_prebuild below, that tests nearly the same thing,
    // and tests versus the same output file.
    "BCR=86237 GEX=85679 LVARSP=gex_max,gex_med,n_gex,CD19_ab CELLS=3 PER_BC",
    // 18. makes sure cross filtering is isn't applied to two samples from same donor
    "BCR=123085:123089 CDR3=CVRDEGGARPNKWNYEGAFDIW",
    // 19. there was a bug that caused twosie to be deleted, and there was foursie junk
    "BCR=123085 CDR3=CARRYFGVVADAFDIW",
    // 20. example affected by whitelist (gel bead oligo contamination) filtering
    "BCR=52177 AMINO=cdr3 PER_BC CDR3=CATWDDSLSGPNWVF",
    // 21. test MIN_CHAINS_EXACT
    "BCR=123089 CDR3=CGTWHSNSKPNWVF MIN_CHAINS_EXACT=3",
    // 22. there was a false positive clonotype
    "BCR=\"165807;165808\" FAIL_ONLY=true EXPECT_NULL",
    // 23. here we were generating a fake alternate allele
    "BCR=83808 CDR3=CAREGRGMVTTNPFDYW MIN_CELLS_EXACT=30",
    // 24. an example that uses IGHE
    "BCR=52177 CDR3=CSTGWGLDFDFWSGYYTAGYHW",
    // 25. add mouse B6 example that had messed up constant regions
    "TCR=74396 MOUSE CVARSP=cdiff CDR3=CASSDAGDTQYF",
    // 26. tests multiple datasets and also LVARS=ncells,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    // Here we make sure that non-specification of MIX_DONORS works.
    "BCR=\"123085;123089\" CDR3=CVKDRVTGTITELDYW",
    // 27. tests SUMMARY and NOPRINT
    "BCR=123085 SUMMARY SUMMARY_CLEAN NOPRINT",
    // 28. tests BARCODE option
    "BCR=165807 BARCODE=CCCATACGTGATGATA-1,TCTATTGAGCTGAAAT-1",
    // 29. tests parenthesized variable in F, SUM and MEAN
    "BCR=86237 GEX=85679 LVARSP=IGHV3-7_g F=\"(IGHV3-7_g)>=4.5\" MIN_CHAINS=2 SUM MEAN",
    // 30. tests d_univ and d_donor
    "BCR=123085 CVARSP=d_univ,d_donor CDR3=CVKDRVTGTITELDYW",
    // 31. tests Cell Ranger 3.1 output
    "BCR=../3.1/123085 CDR3=CVKDRVTGTITELDYW",
    // 32. tests Cell Ranger 2.0 output and RE
    "BCR=../2.0/124550 CDR3=CAREPLYYDFWSAYFDYW RE",
    // 33. tests SCAN
    "BCR=123085 GEX=123201 LVARSP=IGHV1-69D_g MIN_CELLS=10 SCAN=\"(IGHV1-69D_g)>=100,(IGHV1-69D_g)<=1,t-10*c>=0.1\" NOPRINT",
    // 34. tests honeycomb plot
    // (This yields a lot of output so will be annoying to debug if something changes.)
    "BCR=123085:123089 MIN_CELLS=50 PLOT=\"stdout,s1->blue,s2->red\" NOPRINT",
];

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// The following is a single test, containing many subtests, each of which is a regression test
// for a given enclone command line.

#[cfg(not(debug_assertions))]
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
        let mut expect_null = false;
        if test.contains(" EXPECT_NULL") {
            test = test.replace(" EXPECT_NULL", "");
            expect_null = true;
        }
        let testn = test.replace("\"", "");
        let mut log = Vec::<u8>::new();
        let out_file = format!("test/inputs/enclone_test{}_output", it + 1);
        if !path_exists(&out_file) {
            fwriteln!(log, "\nYou need to create the output file {}.\n", out_file);
            fwriteln!(
                log,
                "Do this by executing the following command from \
                 cellranger/lib/rust/enclone:\n"
            );
            emit_bold_escape(&mut log);
            fwriteln!(
                log,
                "enclone PRE=test/inputs/version{} {} \
                 > test/inputs/enclone_test{}_output\n",
                TEST_FILES_VERSION,
                test,
                it + 1
            );
            emit_end_escape(&mut log);
            fwriteln!(log, "and then adding/committing the new file.");
            res.2 = stringme(&log);
        } else {
            let old = read_to_string(&out_file).unwrap();
            let args = testn.split(' ').collect::<Vec<&str>>();
            let mut new = Command::new("target/release/enclone");
            let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
            for i in 0..args.len() {
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
            if old == new2 {
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
                fwriteln!(log, "stderr has {} lines:", new_err.len());
                for i in 0..new_err.len() {
                    fwriteln!(log, "{}", new_err[i]);
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
                emit_bold_escape(&mut log);
                fwriteln!(
                    log,
                    "enclone PRE=test/inputs/version{} {} \
                     > test/inputs/enclone_test{}_output\n",
                    TEST_FILES_VERSION,
                    test,
                    it + 1
                );
                emit_end_escape(&mut log);
                fwrite!(log, "and then committing the changed file.  ");
                fwriteln!(
                    log,
                    "You can then retest using:\n\n\
                     cargo test --release -p enclone enclone  -- --nocapture"
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

// Test that PREBUILD works.  This reuses the output of test 17 above, so if you change that,
// you also have to change this.
//
// WARNING: if this test in interrupted, then you could accidentally be left with matrix.bin,
// and you would need to delete that file.

#[cfg(not(debug_assertions))]
#[test]
fn test_enclone_prebuild() {
    PrettyTrace::new().on();
    let t = Instant::now();

    // See if we're in a broken state.

    if path_exists(&format!(
        "test/inputs/version{}/85679/outs/raw_feature_bc_matrix/matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!(
            "\nenclone_test_prebuild: the file matrix.bin already exists.\n\
             Perhaps a previous run of this test was interrupted.  Please delete the file."
        );
    }

    // First pass: run with NH5.

    let test_id = 17;
    let it = test_id - 1;
    let testn = format!("{} NH5", TESTS[it]);
    let out_file = format!("test/inputs/enclone_test{}_output", test_id);
    let old = read_to_string(&out_file).unwrap();
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new("target/release/enclone");
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
            "\nenclone_test_prebuild: first pass output has changed.\n\
             You may want to add more info to this failure message.\n\
             And don't forget to remove matrix.bin.\n"
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }
    if !path_exists(&format!(
        "test/inputs/version{}/85679/outs/raw_feature_bc_matrix/matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!("\nenclone_test_prebuild: did not create matrix.bin.");
    }

    // Second pass: run without PREBUILD but using the matrix.bin that the first pass created.

    let testn = TESTS[it];
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new("target/release/enclone");
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
             And don't forget to remove matrix.bin.\n"
        );
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }

    // Clean up: delete matrix.bin.

    std::fs::remove_file(&format!(
        "test/inputs/version{}/85679/outs/raw_feature_bc_matrix/matrix.bin",
        TEST_FILES_VERSION
    ))
    .unwrap();
    println!("\nused {:.2} seconds in enclone_test_prebuild", elapsed(&t));
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This test runs enclone for a few test inputs, with LOUPE output
// turned on. It will then read both the bincode and proto file created
// and asserts that we get the same data structure either way.

#[cfg(not(debug_assertions))]
#[test]
fn test_proto_write() -> Result<(), Error> {
    let tests = vec!["BCR=123085", "TCR=101287"];
    let pre_arg = format!("PRE=test/inputs/version{}", TEST_FILES_VERSION);
    let binary_arg = format!("BINARY={}.bin", LOUPE_OUT_FILENAME);
    let proto_arg = format!("PROTO={}.proto", LOUPE_OUT_FILENAME);
    for t in tests.iter() {
        // FIXME: It would be nicer to use the enclone API here
        std::process::Command::new("target/release/enclone")
            .args(&[&pre_arg, *t, &binary_arg, &proto_arg])
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));
        let outputs_proto = read_proto(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(format!("{}.bin", LOUPE_OUT_FILENAME));
        std::fs::remove_file(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        std::fs::remove_file(format!("{}.bin", LOUPE_OUT_FILENAME))?;
        assert!(outputs_proto == outputs_bin);
    }

    Ok(())
}
