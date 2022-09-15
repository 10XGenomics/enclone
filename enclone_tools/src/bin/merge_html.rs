// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Build html files by inserting other html files.
//
// If supplied the single argument BUILD, also rebuild from source.

use enclone_core::parse_bsv;
use enclone_tools::html::{edit_html, insert_html};
use io_utils::{fwrite, open_for_read, open_for_write_new};
use itertools::Itertools;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::env;
use std::fs::{read_dir, File};
use std::io::{BufRead, Write};
use std::process::Command;
use string_utils::{strme, TextUtils};
use vector_utils::{bin_member, unique_sort};

// List of examples on site.

pub const SITE_EXAMPLES: [(&str, &str); 28] = [
    // 1.
    // Do not use NH5 because the bin file is too big for git.
    (
        "pages/auto/clonotype_with_gex.html",
        "BCR=123085 CDR3=CTRDRDLRGATDAFDIW GEX=123217 H5 LVARSP=gex,IGHV3-49_g NUMI \
         HTML=\"enclone example with gex\"",
    ),
    // 2.
    (
        "pages/auto/illusory1.html",
        "BCR=128037,128040 NCROSS CDR3=CARGGTTTYFISW NGROUP NUMI NUMI_RATIO \
         HTML=\"illusory clonotype expansion 1\"",
    ),
    // 3.
    (
        "pages/auto/illusory2.html",
        "BCR=128037,128040 CDR3=CARGGTTTYFISW NGROUP NUMI NUMI_RATIO \
      HTML=\"illusory clonotype expansion 2\"",
    ),
    // 4.
    (
        "pages/auto/illusory3.html",
        "BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW NGROUP NUMI NUMI_RATIO \
         HTML=\"illusory clonotype expansion 3\"",
    ),
    // 5.
    (
        "pages/auto/illusory4.html",
        "BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW PER_CELL LVARSP=gex,cred MIN_CHAINS_EXACT=2 NUMI \
         NUMI_RATIO NGROUP HTML=\"illusory clonotype expansion 4\"",
    ),
    // 6.
    (
        "pages/auto/illusory5.html",
        "BCR=128040 GEX=127801 BC=testx/inputs/128024_cells.csv \
         CDR3=CARGGTTTYFISW PER_CELL NUMI NUMI_RATIO \
         LVARSP=gex,cred,T CHAINS_EXACT=2 NGROUP HTML=\"illusory clonotype expansion 5\"",
    ),
    // 7.
    (
        "img/samples.svg",
        "BCR=123085:123089 MIN_CELLS=10 PLOT=\"stdout,s1->blue,s2->red\" NOPRINT \
         LEGEND=blue,123085,red,123089",
    ),
    // 8.
    (
        "img/iso.svg",
        "BCR=123085,123089 MIN_CELLS=5 MIN_CHAINS_EXACT=2 NOPRINT PLOT_BY_ISOTYPE=stdout",
    ),
    // 9.
    (
        "pages/auto/tree_example.html",
        "BCR=123085 TREE COMPLETE CDR3=CARDQNFDESSGYDAFDIW LVARSP=dref HTML",
    ),
    // 10.
    (
        "pages/auto/mait_example.html",
        "TCR=101287 LVARSP=mait CDR3=CSAGQGDTEAFF HTML",
    ),
    // 11.
    (
        "pages/auto/foursie1.html",
        "BCR=123085 CDR3=CARRYFGVVADAFDIW NFOURSIE_KILL HTML",
    ),
    // 12.
    (
        "pages/auto/foursie2.html",
        "BCR=123085 CDR3=CARRYFGVVADAFDIW HTML",
    ),
    // 13.
    (
        "img/quad_hive.svg",
        "BCR=123085:123089 PLOT=\"stdout,s1->blue,s2->red\" QUAD_HIVE NOPRINT",
    ),
    // 14.
    (
        "img/two_genes.svg",
        "BCR=123085 GEX=123217 NOPRINT PLOTXY_EXACT=HLA-A_g,CD74_g,stdout H5",
    ),
    // 15.
    (
        "pages/auto/variable_demo.html",
        "BCR=123085 CDR3=CALMGTYCSGDNCYSWFDPW PER_CELL POUT=stdouth PCELL PCOLS=barcode,u1,u_cell1 \
         HTML=\"variable demo\"",
    ),
    // 16.
    (
        "pages/auto/variable_demo2.html",
        "BCR=123085 CDR3=CALMGTYCSGDNCYSWFDPW POUT=stdouth PCOLS=barcodes,u1 \
         HTML=\"variable demo2\"",
    ),
    // 17.
    (
        "pages/auto/d_gene_example1.html",
        "BCR=123085 CVARS=d1_name,d2_name,d_Δ CDR3=CTRDRDLRGATDAFDIW \
         HTML=\"D gene example1\"",
    ),
    // 18.
    (
        "pages/auto/d_gene_example1b.html",
        "BCR=123085 CVARS=d1_name,d2_name,d_Δ CDR3=CAREGGVGVVTATDWYFDLW COMPLETE \
         HTML=\"D gene example1b\"",
    ),
    // 19.
    (
        "pages/auto/d_gene_example2.html",
        "BCR=123085 GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT SUMMARY SUMMARY_CLEAN \
         HTML=\"D gene example2\"",
    ),
    // 20.
    (
        "pages/auto/align_example.html",
        "BCR=123085 ALIGN1 CDR3=CARYIVVVVAATINVGWFDPW CVARSP=d1_name \
         HTML=\"ALIGN example\"",
    ),
    // 21.
    (
        "pages/auto/jun_align_example.html",
        "BCR=123085 JALIGN1 CDR3=CARYIVVVVAATINVGWFDPW CVARSP=d1_name \
         HTML=\"JALIGN example\"",
    ),
    // 22.
    (
        "pages/auto/vddj.html",
        "BCR=165808 JALIGN1 CDR3=CARAYDILTGYYERGYSYGWGFDYW \
         HTML=\"VDDJ example\"",
    ),
    // 23.
    (
        "img/sim_mat_plot.svg",
        "BCR=123085 GEX=123217 SIM_MAT_PLOT=stdout,CDKN1A_g,CDKN1B_g,RBX1_g,IGLC1_g,IGLV3-21_g \
         NOPRINT H5",
    ),
    // 24.
    (
        "img/twin_plot.svg",
        "BCR=123085:123089 PLOT_BY_ISOTYPE=stdout SPLIT_PLOT_BY_ORIGIN NOPRINT",
    ),
    // 25.
    (
        "img/var.png",
        "BCR=123085 MIN_CELLS=10 HONEY=out=stdout.png,color=var,u_cell1 NOPRINT NO_NEWLINE",
    ),
    // 26.
    (
        "img/by_dataset.svg",
        "BCR=123085,123089,124547 MIN_CELLS=5 HONEY=out=stdout,color=dataset NOPRINT",
    ),
    // 27.
    (
        "pages/auto/var_def.html",
        r###"BCR=86237 GEX=85679 VAR_DEF="sum:CD19_ab + CD25_ab" LVARSP=CD19_ab,CD25_ab,sum CDR3=CARSFFGDTAMVMFQAFDPW PER_CELL FOLD_HEADERS HTML"###,
    ),
    // 28.
    (
        "img/cat_var.svg",
        "BCR=123085 HONEY=out=stdout,color=catvar,v_name1+v_name2,maxcat:10 NOPRINT CHAINS_EXACT=2",
    ),

// Notes on how to add to the above SITE_EXAMPLES:
//
// Be very careful: there are svg and html examples above.  Mimic one or the other.
//
// 1. cargo b
// 2. merge_html BUILD
// 3. ./build

];

fn main() {
    PrettyTrace::new().on();

    // Build from source.

    let args: Vec<String> = env::args().collect();
    if args.len() == 2 && args[1] == "BUILD" {
        let mut results = Vec::<(usize, Vec<u8>)>::new();
        for i in 0..SITE_EXAMPLES.len() {
            results.push((i, Vec::new()));
        }
        results.par_iter_mut().for_each(|r| {
            let i = r.0;
            let test = SITE_EXAMPLES[i].1;
            let args = parse_bsv(test);
            let new = Command::new("target/debug/enclone")
                .args(&args)
                .arg("MAX_CORES=24")
                .output()
                .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
            if new.status.code() != Some(0) {
                eprint!(
                    "\nenclone_site_examples: example {} failed to execute, stderr =\n{}",
                    i + 1,
                    strme(&new.stderr),
                );
                std::process::exit(1);
            }
            r.1 = new.stdout;
        });
        for i in 0..SITE_EXAMPLES.len() {
            // Move file to the site location.

            let example_name = SITE_EXAMPLES[i].0;
            let out_file = example_name.to_string();
            let mut file = File::create(&out_file).unwrap();
            file.write_all(&results[i].1).unwrap();
        }
    }
    let mut site_ex = Vec::<String>::new();
    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let out_file = example_name.to_string();
        site_ex.push(out_file);
    }
    unique_sort(&mut site_ex);

    // Apply insert_html.

    let all = read_dir("pages").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        if f.ends_with(".html.src") {
            let mut level = 2;
            let mut target = format!("pages/auto/{}.html", f.between("/", "."));
            if target == *"pages/auto/index.html" {
                target = "index.html".to_string();
                level = 0;
            }
            insert_html(f, &target, false, level);
        }
    }

    // This is ugly.  Edit the html pages.

    let all = read_dir("pages/auto").unwrap();
    let mut allx = vec!["index.html".to_string()];
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap().to_string();
        if f.ends_with(".html") && !bin_member(&site_ex, &f) {
            allx.push(f);
        }
    }
    for f in allx {
        let mut edited = false;
        let mut lines = Vec::<String>::new();
        {
            let h = open_for_read![&f.to_string()];
            for line in h.lines() {
                let s = line.unwrap();
                lines.push(s.clone());
                if s.contains("googletag") {
                    edited = true;
                }
            }
        }
        if !edited {
            let x = format!("{}\n", lines.iter().format("\n"));
            let mut h = open_for_write_new![&f.to_string()];
            fwrite!(h, "{}", edit_html(&x));
        }
    }
}
