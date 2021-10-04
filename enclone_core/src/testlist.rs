// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Information about enclone tests.

pub fn enclone_testdata() -> String {
    include_str!["enclone.testdata"].to_string()
}
pub fn enclone_testdata_public_bcr_human() -> String {
    include_str!["testdata.public.bcr.human"].to_string()
}
pub fn enclone_testdata_public_tcr_human() -> String {
    include_str!["testdata.public.tcr.human"].to_string()
}
pub fn enclone_testdata_public_tcr_mouse() -> String {
    include_str!["testdata.public.tcr.mouse"].to_string()
}
pub fn enclone_testdata_public_gex_human() -> String {
    include_str!["testdata.public.gex.human"].to_string()
}

pub const TEST_FILES_VERSION: u8 = 15;

// Unaccounted time test.  This is separated out so that we can avoid running it in parallel with
// other tests, to reduce the sporadic failure rate.  Also this is off in GitHub Actions.

pub const UNAC_TESTS: [&str; 1] = [
    // 1. enforce no unaccounted time
    r###"BCR=123085 COMPE UNACCOUNTED NOPRINT EXPECT_OK"###,
];

// Tests that are affected by the D region alignment algorithm.  All such tests should go here.

pub const DTESTS: [&str; 15] = [
    // 1. test ALIGN_2ND<n>
    r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW ALIGN_2ND1 CVARS=d1_name"###,
    // 2. test JALIGN_2ND<n>
    r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW JALIGN_2ND1 CVARS=d2_name"###,
    // 3. test ALIGN_JALIGN_CONSISTENCY
    r###"BCR=123085 CELLS=1 CHAINS=2 ALIGN1 JALIGN1 ALIGN_JALIGN_CONSISTENCY AMINO=cdr3
         PLAIN NOPAGER EXPECT_OK"###,
    // 4. test D_INCONSISTENT, and lock number of inconsistencies
    r###"BCR=123085 D_INCONSISTENT CVARS=d1_name COMPLETE NGROUP"###,
    // 5. the JALIGN1 in this example had a boundary location that was off by one
    r###"BCR=165807 JALIGN1 AMINO=cdr3 CVARS=d1_score,d2_score CDR3=CAKEYYDFWSGYSDVRGVIPNIDYW"###,
    // 6. the JALIGN1 in this example had a boundary location that was off by one
    r###"BCR=123085 CELLS=2 JALIGN1 AMINO=cdr3 CVARS=d1_name CDR3=CAKAGPTESGYYVWYFDLW"###,
    // 7. test d_inconsistent_{%,n}
    r###"BCR=123085 GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT SUMMARY SUMMARY_CLEAN"###,
    // 8. test ALIGN<n>
    r###"BCR=123085 CDR3=CKVMLYDSRGSDYYYVMDVW ALIGN1 CVARS=d1_name"###,
    // 9. test ALIGN<n> and JALIGN<n>, case where there's a D segment
    r###"BCR=85333 ALIGN1 JALIGN1 CDR3=CARGYDFWSGYLVGNWAGDYYYYMDVW"###,
    // 10. test ALIGN<n> and JALIGN<n>, case where there is no D segment
    r###"BCR=85333 ALIGN1 JALIGN1 CDR3=CAKGKGFRNYYYYMDVW"###,
    // 11. test d1 etc.
    r###"BCR=123085 CVARS=d1_name,d2_name,d_Δ,d_delta AMINO=cdr3 CDR3=CARVRDILTGDYGMDVW"###,
    // 12. test GROUP_VDJ_REFNAME_HEAVY (deprecated but supported)
    r###"BCR=86237 GROUP_VDJ_REFNAME_HEAVY CDR3="CAKAVAGKAVAGGWDYW|CAKVSTGIAVAGPGDYW" COMPLETE"###,
    // 13. test GROUP_VJ_REFNAME_HEAVY (deprecated but supported)
    r###"BCR=86237 GROUP_VJ_REFNAME_HEAVY CDR3="CARGVLWFGELGAFDIW|CARAGLGVVLAARGAFDIW""###,
    // 14. test placement of indel, needed shifting right
    r###"BCR=123085 CELLS=1 CHAINS=2 AMINO=cdr3 JALIGN2 CDR3=CAKDKSRPPTHYYGSGSYYSRILDNW"###,
    // 15. test placement of indel, needed shifting left
    r###"BCR=123085 CELLS=1 CHAINS=2 AMINO=cdr3 JALIGN2 CDR3=CARMAQFYSGSGTYYIGPYYFEYW"###,
];

// Tests that are affected by the grouping algorithm.  All such tests should go here, if not
// already in DTESTS.

pub const GTESTS: [&str; 13] = [
    // 1. test 5/8 for newline correctness (this grouping option deprecated but supported)
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN SET_IN_STONE"###,
    // 2. test 6/8 for newline correctness (this grouping option deprecated but supported)
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN NGROUP SET_IN_STONE"###,
    // 3. test 7/8 for newline correctness (this grouping option deprecated but supported)
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML SET_IN_STONE"###,
    // 4. test 8/8 for newline correctness (this grouping option deprecated but supported)
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML NGROUP SET_IN_STONE"###,
    // 5. test of GROUP
    r###"BCR=123085 GROUP=vj_refname,cdr3_aa_heavy≥80%,cdr3_aa_light≥80% CVARS=cdr3_len
         AMINO=cdr3 CDR3="CARHLQWELP.*W""###,
    // 6. test of GROUP
    r###"BCR=123085 GROUP=vj_refname,len,cdr3_len MIN_GROUP=2 MIN_CHAINS=2 CDR3="CQQSY.*TLATF"
         CVARS=cdr3_len"###,
    // 7. test of GROUP
    r###"BCR=123085 GROUP=cdr3_aa_heavy≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
         CDR3=CARPKSDYIIDAFDIW"###,
    // 8. test of GROUP
    r###"BCR=123085 GROUP=cdr3_aa_light≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
         CDR3=CQTWGTGPWVF"###,
    // 9. test of GROUP
    r###"BCR=123085 GROUP=vj_refname,aa_heavy≥100% MIN_GROUP=2 MIN_CHAINS=2 CVARS=cdr3_len
         CDR3=CARVPYYYDRSYYYYGMDVW"###,
    // 10. test of AGROUP
    r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CARHSYSSGWYDEWDYW
         AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=top=2"###,
    // 11. test of AGROUP
    r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CAKDGGEHYYDSSGYYASYYFDYW 
         AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=14"###,
    // 12. test of AGROUP
    r###"BCR=123085 AGROUP AG_CENTER=from_filters CDR3=CAKDGGEHYYDSSGYYASYYFDYW 
         AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=13"###,
    // 13. test of AGROUP
    r###"BCR=123085 AGROUP AG_CENTER=copy_filters MIN_CELLS=2 MAX_CELLS=2
         AG_DIST_FORMULA=cdr3_edit_distance AG_DIST_BOUND=max=3 MIN_GROUP=2"###,
];

// Crash tests.  These are tests to make sure that certain options do not result in a crash, even
// when run on a large and complex dataset.  The options are in groups because not all are
// compatible with each other.  The datasets are defined by a single fixed list, to be enlarged
// over time based on discovery of pathologies in particular datasets.  In general these datasets
// are not public.  All run with certain shared options.

pub const CRASH_DATA: &str = "BCR=\"45987;123085;testx/inputs/flaky\"";
pub const CRASH_OPTS: &str = "NOPRINT BUILT_IN EXPECT_OK NO_PRE NFORCE";
pub const CRASH_SETS: [&str; 6] = [
    /* 1 */ "CONP SEQC SUM MEAN BARCODES DIFF_STYLE=C1 GROUP_VJ_REFNAME",
    //
    /* 2 */ "CONX FULL_SEQC DIFF_STYLE=C2 POUT=stdout PCOLS=count_CAR",
    //
    /* 3 */
    "AMINO=fwr1,cdr1,fwr2,cdr2,fwr3,cdr3,fwr4 CVARS=d1_name,d2_name,d_delta,d_Δ,cigar",
    //
    /* 4 */
    "PLOT_BY_ISOTYPE=stdout MIN_CELLS=3 GROUP_VJ_REFNAME_HEAVY ALIGN1 JALIGN1",
    //
    /* 5 */
    "GROUP_VDJ_REFNAME_HEAVY GVARS=d_inconsistent_%,d_inconsistent_n",
    //
    /* 6 */
    "GROUP=vj_refname,cdr3_aa_heavy≥90%,cdr3_aa_light≥90%",
];

// Test using datasets that are either in the extended public dataset collection, or which are
// not publicly avaiable, or which require samtools.

pub const EXTENDED_TESTS: [&str; 36] = [
    // 1. test that used to crash on a particular barcode; this also gave the wrong
    // answer for an insertion until it was fixed
    r###"BCR=40955 NCELL BARCODE=GCGCAGTCAAAGTGCG-1 AMINO=cdr3 NO_PRE NFORCE"###,
    // 2. tests nd2
    r###"BCR=47199,47200,47212 AMINO=cdr3 NCROSS LVARS=nd2 CDR3=CVKGKSGSFWYYFENW
         NO_PRE NFORCE"###,
    // 3. test sec and mem [requires samtools]
    r###"BCR=123085 GEX=123217 LVARSP=sec,mem CDR3=CVKDRVTGTITELDYW H5"###,
    // 4. test MOUSE + IMGT; note that specifying by number forces BCR+TCR reference checks
    r###"70838 MOUSE NOPRINT SUMMARY SUMMARY_CLEAN IMGT ACCEPT_BROKEN NO_PRE NFORCE"###,
    // 5. this crashed (and didn't check if this is in extended public dataset collection)
    r###"BCR=83809 CDR3=CARVSLGYCSGGSCNSNYYFDYW NO_PRE NFORCE"###,
    // 6. this crashed (and didn't check if this is in extended public dataset collection)
    r###"BCR=47680 BARCODE=CGCCAAGTCCATGAAC-1 NO_PRE NFORCE"###,
    // 7. this crashed (and didn't check if this is in extended public dataset collection)
    r###"BCR=99640 BARCODE=CAGTAACCATGTCGAT-1 NO_PRE NFORCE"###,
    // 8. test MOUSE BCR + our reference (this crashed) -- LOOKS REDUNDANT NOW
    r###"BCR=70838 MOUSE NOPRINT NO_PRE NFORCE EXPECT_NULL"###,
    // 9. this clonotype included a junk chain before we made a change, and test "/outs"
    r###"TCR=163911/outs CDR3=CAPSAGDKIIF AMINO=donor NO_PRE NFORCE"###,
    // 10. test case where digit rows are just barely present
    r###"TCR=163911 CDR3=CASSLVQPSTDTQYF AMINO=donor NO_PRE NFORCE"###,
    // 11. this added because it got better when a noise filter was added, also tests u_max
    r###"TCR=163914 CDR3=CASSLVQPSTDTQYF CVARSP=u_max NO_PRE NFORCE"###,
    // 12. this added because it got better when a noise filter was added; also test FASTA
    r###"TCR=163914 CDR3=CAFRGGSYIPTF FASTA=stdout NO_PRE NFORCE"###,
    // 13. this added because it got better when a bug in bads detection was fixed
    r###"TCR=163914 CDR3=CASRLGGEETQYF NO_PRE NFORCE"###,
    // 14. this crashed before a bug was fixed
    r###"BCR=1021341 NCELL CDR3=CQQANSYPLTF SEG=IGHV1-69D NO_PRE NFORCE"###,
    // 15. test that LVARSP=gex fails on Ab-only data
    r###"BCR=1031851 GEX=1031779 NGEX LVARSP=gex EXPECT_FAIL NO_PRE NFORCE"###,
    // 16. test Ab-only data
    r###"BCR=1031851 GEX=1031779 NGEX LVARSP=n_gex,CD19_ab
         CDR3="CARDELDILTGYNIPTFGGCVYW|CAHHGSARYSSSWHAAPGPYYFDYW" BUILT_IN NO_PRE NFORCE"###,
    // 17. test for very long (120 amino acid) CDR3
    // Note that this long CDR3 is likely part of a nonproductive chain.  The test is here because
    // there may be long productive CDR3 sequences in data from other species, although we do not
    // have such data.
    r###"BCR=1020665 BUILT_IN REPROD CVARSP=cdr3_len CDR3=CARDGGGQPFDLW AMINO= NO_PRE NFORCE"###,
    // 18. an example that triggered an internal inconsistency test, which we subsequently removed;
    // there are three chains and the middle one was the problem
    r###"TCR=48602 BARCODE=CCAGCGAAGTGTTGAA-1 REPROD NO_PRE NFORCE"###,
    // 19. Make sure that POUT works on full dataset.
    // If we experience failures on other lena ids, we can add them to this list.
    r###"BCR="86213;86237" RE POUT=/dev/null NOPRINT EXPECT_OK NO_PRE NFORCE"###,
    // 20. Make sure that FP join output includes join error details.
    // If somehow we fix the FP join occurring here, another one should be substituted.
    r###"BCR="131036;140707" ANN SHOW_BC FAIL_ONLY=true PRINT_FAILED_JOINS MIX_DONORS
         NO_PRE NFORCE"###,
    // 21. the result of this changed when sub_alts was changed
    r###"BCR="40086;132888" SEG=IGHV3-43 MIX_DONORS MAX_DIFFS=80 CDR3=CVKGDWGSAFDIW
         NO_PRE NFORCE"###,
    // 22. clonotype that was two clonotypes before raising MAX_DIFFS to 60
    r###"BCR=1084461-1084462 CDR3=CAKEFGNGGFDTFDIW NO_PRE NFORCE"###,
    // 23. test BCR_GEX and GD_BC
    r###"BCR_GEX=1089851 GD_BC=1089848 NOPRINT NO_PRE NFORCE EXPECT_OK"###,
    // 24. This used to appear as a four-chain clonotype, and is now split.
    r###"BCR=123085,123090 BUILT_IN BARCODE=AAAGTAGCAAGCCATT-1,ATGGGAGTCCATGAGT-1 NO_PRE NFORCE"###,
    // 25. Test a tweak to the weak chains filter.  This should have two chains.
    r###"BCR=174957 CDR3=CARPRGYCSGGSCFPFASW BUILT_IN NO_PRE NFORCE"###,
    // 26. crashed at one point
    r###"BCR=123085,123086 GEX=123749,123750 LVARSP=pe1 BUILT_IN NOPRINT EXPECT_OK NO_PRE
         NFORCE"###,
    // 27. parseable value for fwr4_aa was wrong
    r###"BCR=1117070 AMINO=fwr4 CDR3=CAKDVNGYSSGWAFENW POUT=stdout PCOLS=fwr4_aa1 NO_PRE NFORCE"###,
    // 28. conp value was truncated
    r###"BCR=1117069 CONP CDR3=CVRDPPEELELFDYW NO_PRE NFORCE"###,
    // 29. Test PCHAINS=max.  For this we need a clonotype having at least five chains, and the
    // question is whether the header line represents cvars for all the chains.  The output of
    // this is expected to change whenever variables are added.
    r###"BCR=140696,140697,140701,140704 MIN_CHAINS=5 BUILT_IN AMINO= FOLD_HEADERS LVARS=
         POUT=stdout PCHAINS=max NOPRINT NO_PRE NFORCE"###,
    // 30. test on PD multi pipestance; failed before bug fix
    r###"BCR_GEX=1084461 NOPRINT EXPECT_OK NO_PRE NFORCE"###,
    // 31. previously this yielded a disconnected clonotype
    r###"BUILT_IN BCR=140699,140705-140706 AMINO=cdr3 CDR3="CAKDRQAGGIGEVDDW|CARDRVPGGIGEVDYW"
         NO_PRE NFORCE"###,
    // 32. test fb variables
    r###"BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX LVARSP=fb2,fb2_n,Ag_APC-C0956_ab PER_CELL
         AMINO=cdr3 CVARS= FOLD_HEADERS POUT=stdouth PCOLS=fb2,fb2_n,fb2_n_cell PCELL 
         CDR3=CAKLLVALHYW NO_PRE NFORCE"###,
    // 33. test DVARS
    r###"TCR_GEX=1175300-1175301 DVARS=Ag_PE-C0951_ab_cellular_u,Ag_PE-C0951_ab_cellular_r
         NOPRINT SUMMARY SUMMARY_CLEAN NFORCE"###,
    // 34. test signature filtering ON
    // This is super annoying.  It is the only case we could find where signature filtering has
    // an effect.  However, other changes will also affect this.  See the next test and make sure
    // that the results are different from it.
    r###"BCR=83808-83809 BUILT_IN NOPRINT SUMMARY SUMMARY_CLEAN NO_PRE NFORCE"###,
    // 35. test signature filtering OFF
    // This is super annoying.  It is the only case we could find where signature filtering has
    // an effect.  However, other changes will also affect this.  See the previous test and make
    // sure that the results are different from it.
    r###"BCR=83808-83809 NSIG BUILT_IN NOPRINT SUMMARY SUMMARY_CLEAN NO_PRE NFORCE"###,
    // 36. test MIN_GROUP_DONORS
    r###"BCR="40953;43899" MIX_DONORS MIN_GROUP=2 NFORCE
         GROUP="cdr3_len,cdr3_aa_heavy>=85%,cdr3_aa_light>=85%,vj_refname" MIN_GROUP_DONORS=2"###,
];

// Tests of internal features.

pub const INTERNAL_TESTS: [&str; 4] = [
    // 1. gave wrong result
    r###"123085 CDR3=CARDRIAGRFGYGMDVW NFORCE"###,
    // 2. test human + IMGT; note that specifying by number forces BCR+TCR reference checks
    r###"123085 REQUIRE_UNBROKEN_OK IMGT ACCEPT_BROKEN EXPECT_NULL"###,
    // 3. test mouse + IMGT; note that specifying by number forces BCR+TCR reference checks
    r###"70838 REQUIRE_UNBROKEN_OK IMGT ACCEPT_BROKEN MOUSE NO_PRE NFORCE EXPECT_NULL"###,
    // 4. this crashed; it is not exactly an internal feature test but uses an internal feature
    // (IMGT) to exhibit the phenomenon
    r###"BCR=123085 IMGT RE ACCEPT_BROKEN POUT=stdout PCELL BARCODE=AGCAGCCCATTAGGCT-1
         EXPECT_OK"###,
];

// List of examples in documentation.

pub const EXAMPLES: [&str; 2] = [
    // 1.
    r###"BCR=123089 CDR3=CARRYFGVVADAFDIW"###,
    // 2.
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 H5 LVARSP=gex,IGHV2-5_g_μ CDR3=CALMGTYCSGDNCYSWFDPW"###,
];

// List of examples on site.

pub const SITE_EXAMPLES: [(&str, &str); 25] = [
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
        "BCR=123085:123089 PLOT_BY_ISOTYPE=stdout SPLIT_PLOT_BY_ORIGIN NOPRINT H5",
    ),
    // 25.
    (
        "img/var.png",
        "BCR=123085 MIN_CELLS=10 HONEY=out=stdout.png,color=var,u_cell1 NOPRINT NO_NEWLINE",
    ),

// Notes on how to add to the above SITE_EXAMPLES:
//
// Be very careful: there are svg and html examples above.  Mimic one or the other.
//
// 1. cargo b
// 2. merge_html BUILD
// 3. ./build

];
