// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

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

pub const TESTS: [&str; 174] = [
    // 1. tests variant base after CDR3, parseable output
    r###"BCR=123089 CDR3=CVRDRQYYFDYW POUT=stdout
     PCOLS=exact_subclonotype_id,n,v_name1,v_name2,nchains,var_indices_aa1,barcodes"###,
    // 2. tests many donor ref differences, test comp, edit and var and donorn
    r###"BCR=123089 CDR3=CARRYFGVVADAFDIW CVARSP=comp,edit,var AMINO=cdr3,var,share,donorn"###,
    // 3. tests motif in CDR3, CHAINS, u_sum, ulen, flipped args in CVARS, on tiny dataset
    r###"BCR=85333 CDR3="CAA.*" CHAINS=2 CVARS=const,u_sum,ulen"###,
    // 4. tests gex and antibody, FULL_SEQC, ulen, udiff, on tiny dataset
    r###"BCR=86237 GEX=85679 LVARSP=gex,CD19_ab_μ,CD25_ab_μ,IGLV3-1_g_μ,IGLV3-1_g_%,RPS27_g_μ
     CELLS=3 FULL_SEQC SUM MEAN
     CVARSP=ulen,udiff"###,
    // 5. tests TCR and correct grouping of onesies on AGBT Donor 2 dataset
    r###"TCR=101287 MIN_CELLS=100"###,
    // 6. tests AMINO= and vjlen and other things
    r###"BCR=86237 CELLS=3 AMINO= CVARS=u,r,cdr3_dna,cdr3_len,vjlen"###,
    // 7. tests SHM deletion
    r###"BCR=123085 CVARSP=var,clen,cdiff CDR3=CAREPLYYDFWSAYFDYW LVARSP=near,far"###,
    // 8. DUPLICATE; TO REPLACE WITH A NEW TEST
    r###"BCR=123085 CVARSP=var,clen,cdiff CDR3=CAREPLYYDFWSAYFDYW LVARSP=near,far"###,
    // 9. tests PER_CELL and unicode
    r###"BCR=█≈ΠΠΠ≈█ CDR3=CAKGDRTGYSYGGGIFDYW PER_CELL"###,
    // 10. tests multiple datasets and also LVARS=n,origins,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    r###"BCR="123085;123089" CDR3=CVKDRVTGTITELDYW LVARS=n,origins,donors,datasets AMINO=share
     MIX_DONORS"###,
    // 11. tests META, and CONST_IGH + META, which was broken at one point
    r###"META=testx/inputs/test11_meta CDR3=CARSFFGDTAMVMFQAFDPW LVARSP=donors,gex
     CONST_IGH=IGHD"###,
    // 12. test colon lvar in KEEP_CLONO_IF_CELL_MEAN= and test for parsing error at +
    r###"BCR=86237 GEX=85679 LVARSP=g37:IGHV3-7_g_μ KEEP_CLONO_IF_CELL_MEAN="n + g37 >= 5.5"
        MIN_CHAINS=2 NH5"###,
    // 13. DUPLICATE; TO REPLACE WITH A NEW TEST.
    r###"BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113"###,
    // 14. DUPLICATE; TO REPLACE WITH A NEW TEST.
    r###"BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113"###,
    // 15. tests insertion and AMINO range; also this incorrectly reported an insertion before
    // it was fixed
    r###"BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113"###,
    // 16. tests number of cells broken out by dataset
    r###"BCR=123085,123089 LVARS=n,n_123085,n_123089 CDR3=CTRDRDLRGATDAFDIW"###,
    // 17. tests gex with PER_CELL and tests n_gex
    // See also enclone_test_prebuild below, that tests nearly the same thing,
    // and tests versus the same output file.
    r###"BCR=86237 GEX=85679 LVARSP=gex_max,gex,n_gex,CD19_ab_μ CELLS=3 PER_CELL NH5"###,
    // 18. makes sure cross filtering isn't applied to two origins from same donor
    r###"BCR=123085:123089 CDR3=CVRDEGGARPNKWNYEGAFDIW"###,
    // 19. there was a bug that caused a twosie to be deleted, and there was foursie junk
    // There were also some cells that were lost due to a bug in graph filtering.
    r###"BCR=123085 CDR3=CARRYFGVVADAFDIW"###,
    // 20. example affected by whitelist (gel bead oligo contamination) filtering, and test u_Σ
    r###"BCR=52177 AMINO=cdr3 PER_CELL CDR3=CATWDDSLSGPNWVF CVARSP=u_Σ"###,
    // 21. test MIN_CHAINS_EXACT
    r###"BCR=123089 CDR3=CGTWHSNSKPNWVF MIN_CHAINS_EXACT=3"###,
    // 22. there was a false positive clonotype
    r###"BCR="165807;165808" FAIL_ONLY=true EXPECT_NULL"###,
    // 23. here we were generating a fake alternate allele
    r###"BCR=83808 CDR3=CAREGRGMVTTNPFDYW MIN_CELLS_EXACT=30"###,
    // 24. an example that uses IGHE, and test NGROUP
    r###"BCR=52177 CDR3=CSTGWGLDFDFWSGYYTAGYHW NGROUP"###,
    // 25. add mouse B6 example that had messed up constant regions
    r###"TCR=74396 MOUSE CVARSP=cdiff CDR3=CASSDAGDTQYF"###,
    // 26. tests multiple datasets and also LVARS=n,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    // Here we make sure that non-specification of MIX_DONORS works.
    r###"BCR="123085;123089" CDR3=CVKDRVTGTITELDYW"###,
    // 27. tests SUMMARY and NOPRINT
    r###"BCR=123085 SUMMARY SUMMARY_CLEAN NOPRINT"###,
    // 28. tests BARCODE option
    r###"BCR=165807 BARCODE=CCCATACGTGATGATA-1,TCTATTGAGCTGAAAT-1"###,
    // 29. tests KEEP_CLONO_IF_CELL_MAX and parenthesized variable in it, SUM and MEAN
    r###"BCR=123085 GEX=123217 LVARSP=IGHV3-7_g,IGHV3-7_g_μ
        KEEP_CLONO_IF_CELL_MAX="(IGHV3-7_g_μ)>=10000.0" MIN_CHAINS=2 SUM MEAN"###,
    // 30. tests d_univ and d_donor
    r###"BCR=123085 CVARSP=d_univ,d_donor CDR3=CVKDRVTGTITELDYW"###,
    // 31. tests Cell Ranger 3.1 output
    r###"BCR=../3.1/123085 CDR3=CVKDRVTGTITELDYW ACCEPT_BROKEN"###,
    // 32. tests Cell Ranger 2.0 output and RE
    r###"BCR=../2.0/124550 CDR3=CAREPLYYDFWSAYFDYW RE ACCEPT_BROKEN"###,
    // 33. tests SCAN
    r###"BCR=123085 GEX=123217 LVARSP=IGHV1-69D_g_μ MIN_CELLS=10 NGEX
     SCAN="(IGHV1-69D_g_μ)>=100,(IGHV1-69D_g_μ)<=1,t-10*c>=0.1" NOPRINT H5"###,
    // 34. tests honeycomb plot
    // (This yields a lot of output so will be annoying to debug if something changes.)
    r###"BCR=123085:123089 MIN_CELLS=10 PLOT="stdout,s1->red,s2->blue" NOPRINT
     LEGEND=red,"cell from 123085",blue,"cell from 123089""###,
    // 35. tests barcode-by-barcode specification of colors, and tests LEGEND=
    // Note that the specification of PRE overrides our usual specification.
    // (This yields a lot of output so will be annoying to debug if something changes.)
    r###"PRE=../enclone-data/big_inputs/version{TEST_FILES_VERSION},. META=testx/inputs/test35_meta MIN_CELLS=10 MIN_CHAINS_EXACT=2 NOPRINT PLOT=stdout NO_PRE
     LEGEND=red,IGHG1,green,IGHG3,blue,IGHA1,orange,IGHM,black,unassigned"###,
    // 36. tests PCELL and u_Σ in PCOLS (both forms)
    r###"BCR=85333 CDR3=CARDGMTTVTTTAYYGMDVW POUT=stdout PCELL CVARSP=u_Σ
        PCOLS=barcode,const1,const2,u_Σ1,u_sum1"###,
    // 37. tests parseable output of barcodes for a given dataset
    r###"BCR=123085,123089 POUT=stdout PCOLS=123085_barcodes,123089_barcodes
     CDR3=CAVTIFGVRTALPYYYALDVW"###,
    // 38. tests parseable output of barcodes for a given dataset, using PCELL
    r###"BCR=123085,123089 POUT=stdout PCOLS=123085_barcode,123089_barcode PCELL
     CDR3=CAVTIFGVRTALPYYYALDVW"###,
    // 39. tests u and r fields in parseable output, and tests stdouth
    r###"BCR=85333 POUT=stdouth PCOLS=barcode,u1,u_cell1,r2,r_cell2 PCELL PER_CELL CVARSP=r
        CDR3=CAADGGGDQYYYMDVW"###,
    // 40. indel was wrong
    r###"BCR=86237 GEX=85679 LVARSP=IGHV3-7_g_μ F="(IGHV3-7_g_μ)>=4.5" MIN_CHAINS=2 SUM MEAN
        NH5"###,
    // 41. test case for gex_cell
    r###"BCR=86237 GEX=85679 CDR3=CAKAVAGKAVAGGWDYW POUT=stdouth PCOLS=gex_cell PCELL NH5"###,
    // 42. test case that should fail because gex_cell doesn't make sense without gex data
    r###"BCR=85333 CDR3=CQQRSNWPLYTF POUT=stdouth PCOLS=gex_cell PCELL PER_CELL EXPECT_FAIL"###,
    // 43. test case that should fail because _cell variables can't be used in LVARS
    r###"BCR=86237 GEX=85679 CDR3=CAKAVAGKAVAGGWDYW LVARS=gex_cell EXPECT_FAIL"###,
    // 44. test _cell
    r###"BCR=86237 GEX=85679 LVARSP=gex,RPS27_g_μ CELLS=3 POUT=stdouth
        PCOLS=barcode,gex_cell,CD19_ab,CD19_ab_cell NH5 PCELL"###,
    // 45. test ndiff...
    r###"BCR=123085 CVARSP=ndiff1vj,ndiff2vj CDR3=CARDQNFDESSGYDAFDIW"###,
    // 46. test u_μ, u_min, r_μ, r_min and r_max
    r###"BCR=85333 CVARSP=u_μ,u_min,u_max,r,r_μ,r_min,r_max AMINO=cdr3 CDR3=CAADGGGDQYYYMDVW
        POUT=stdouth PCOLS=u_μ1,u_min1,u_max1,r2,r_μ2,r_min2,r_max2"###,
    // 47. this should fail
    r###"BCR=85333 CDR3=CAREEYYYDSSGDAFDIW LVARSP=gex_mean EXPECT_FAIL"###,
    // 48. test gex_mean and gex_Σ and NGEX
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 LVARSP=gex_mean,gex_Σ CDR3=CASRKSGNYIIYW NGEX H5"###,
    // 49. test HTML
    r###"BCR=85333 CDR3=CAAWDDSLNGWVF CHAINS=1 POUT=stdouth PCOLS=barcodes,n FASTA=stdout
        FASTA_AA=stdout HTML=CAAWDDSLNGWVF"###,
    // 50. make sure this doesn't fail
    r###"NOPAGER EXPECT_OK"###,
    // 51. make sure this fails gracefully
    r###"BCR=123085 PLOT=/nonexistent/broken.svg NOPRINT MIN_CELLS=50 EXPECT_FAIL"###,
    // 52. add test for some gene patterns
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 CDR3=CARPKSDYIIDAFDIW MIN_CELLS=10
        LVARSP="(IGHV5-51|IGLV1-47)_g_%,IGH.*_g_%,IG(K|L).*_g_%""###,
    // 53. add test for _% with PER_CELL
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 LVARSP="gex,n_gex,JCHAIN_g_%,IG%:IG.*_g_%" CVARS=u_μ,const
        MIN_CHAINS_EXACT=2 CDR3=CAREGGVGVVTATDWYFDLW PER_CELL"###,
    // 54. make sure this fails gracefully
    r###"BCR=86237 GEX=85679 LVARSP=GERBULXXX123_g_% EXPECT_FAIL"###,
    // 55. test cred
    r###"BCR=86237 GEX=85679 LVARSP=cred PCELL PER_CELL POUT=stdouth PCOLS=cred_cell
        CDR3=CARSFFGDTAMVMFQAFDPW"###,
    // 56. test SVG
    r###"BCR=85333 CDR3=CARDPRGWGVELLYYMDVW SVG NGROUP"###,
    // 57. test 1/8 for newline correctness
    r###"BCR=85333 CDR3="CLLSYSGARVF|CQSADSSGTYKVF" AMINO= PLAIN SET_IN_STONE"###,
    // 58. test 2/8 for newline correctness
    r###"BCR=85333 CDR3="CLLSYSGARVF|CQSADSSGTYKVF" AMINO= PLAIN NGROUP SET_IN_STONE"###,
    // 59. test 3/8 for newline correctness
    r###"BCR=85333 CDR3="CLLSYSGARVF|CQSADSSGTYKVF" AMINO= PLAIN HTML SET_IN_STONE"###,
    // 60. test 4/8 for newline correctness
    r###"BCR=85333 CDR3="CLLSYSGARVF|CQSADSSGTYKVF" AMINO= PLAIN NGROUP HTML SET_IN_STONE"###,
    // 61. test 5/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN SET_IN_STONE"###,
    // 62. test 6/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN NGROUP SET_IN_STONE"###,
    // 63. test 7/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML SET_IN_STONE"###,
    // 64. test 8/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML NGROUP SET_IN_STONE"###,
    // 65. test NCELL
    r###"BCR=86237 NCELL CDR3=CAKTATTLGGYYSHGLDVW MIN_CELLS=2"###,
    // 66. test BC in combination with PER_CELL and PCELL
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T PCELL
        POUT=stdouth PCOLS=barcode,T CDR3=CAKAGPTESGYYVWYFDLW MIN_CELLS=2"###,
    // 67. expect fail if garbage PRE
    r###"PRE=garbage_gerbil_stuff BCR=86237 CELLS=3 EXPECT_FAIL NO_PRE"###,
    // 68. a test of PRE
    r###"PRE=mumbo_jumbo,../enclone-data/big_inputs/version{TEST_FILES_VERSION} BCR=86237 NO_PRE
        CDR3=CARENHPVEYCSSTSCYKAYYYGMDVW"###,
    // 69. another test of pre
    r###"PRE=mumbo_jumbo BCR=../enclone-data/big_inputs/version{TEST_FILES_VERSION}/86237 NO_PRE
        CDR3=CARENHPVEYCSSTSCYKAYYYGMDVW"###,
    // 70. another test of META
    r###"META=mumbo_jumbo EXPECT_FAIL"###,
    // 71. another test of META
    r###"PRE=../enclone-data/big_inputs/version{TEST_FILES_VERSION},testx/inputs META=test11_meta
        CDR3=CARSFFGDTAMVMFQAFDPW LVARSP=donors,gex NO_PRE"###,
    // 72. test SUMMARY_CSV
    r###"BCR=86237 NOPRINT SUMMARY_CSV"###,
    // 73. DUPLICATE; TO REPLACE WITH A NEW TEST
    r###"BCR=86237 NOPRINT SUMMARY_CSV"###,
    // 74. this changed after a bug was fixed; the RE can probably be dropped later when we
    // rerun all the datasets
    r###"BCR=123085 RE CDR3=CARGYEDFTMKYGMDVW POUT=stdouth PCOLS=utr_id2"###,
    // 75. this changed after a bug in RE was fixed, and this is in fact testing RE
    r###"BCR=123085 CDR3=CQQSYSTPRTF RE"###,
    // 76. test PLOT_BY_ISOTYPE
    r###"BCR=123085 MIN_CELLS=10 PLOT_BY_ISOTYPE=stdout NOPRINT MIN_CHAINS_EXACT=2"###,
    // 77. DUPLICATE, SHOULD BE DELETED
    r###"BCR=86237 POUT=/dev/null NOPRINT EXPECT_OK"###,
    // 78. make sure that POUT with PCELL works on full dataset
    r###"BCR=86237 POUT=stdout PCELL EXPECT_OK"###,
    // 79. make sure that POUT works on full dataset with gex
    r###"BCR=86237 GEX=85679 POUT=stdout NGEX NCELL EXPECT_OK"###,
    // 80. make sure that POUT with PCELL works on full dataset with gex
    r###"BCR=86237 GEX=85679 POUT=stdout PCELL NGEX NCELL EXPECT_OK"###,
    // 81. IG:IG.*_g_%_cell and variants in parseable output
    r###"BCR=86237 GEX=85679 CDR3=CARSFFGDTAMVMFQAFDPW POUT=stdouth PCELL
        PCOLS="barcode,IG:IG.*_g_%_cell,IG.*_g_%_cell,IGN:IG.*_g_%,IG.*_g_%""###,
    // 82. test entropy
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 LVARSP=entropy PER_CELL POUT=stdouth PCELL
        PCOLS=barcode,entropy,entropy_cell CDR3=CARAQRHDFWGGYYHYGMDVW"###,
    // 83. test COMPLETE and dref
    r###"BCR=86237 CDR3=CARSFFGDTAMVMFQAFDPW COMPLETE LVARSP=dref"###,
    // 84. test CLUSTAL_AA
    r###"BCR=123085 CDR3=CAADRQLWSRSPGDYIYYGMQVW CLUSTAL_AA=stdout"###,
    // 85. test NALL
    r###"BCR=86237 NALL CDR3=CARAPEDTSRWPQYNYSGLDVW SEG=IGKV3-15"###,
    // 86. test CLUSTAL_DNA
    r###"BCR=86237 CDR3=CARSFFGDTAMVMFQAFDPW CLUSTAL_DNA=stdout"###,
    // 87. test PHYLIP_AA and COLOR=codon
    r###"BCR=123085 CDR3=CAADRQLWSRSPGDYIYYGMQVW PHYLIP_AA=stdout COLOR=codon"###,
    // 88. test PHYLIP_DNA and COLOR=default
    r###"BCR=123085 CDR3=CAADRQLWSRSPGDYIYYGMQVW PHYLIP_DNA=stdout COLOR=property"###,
    // 89. test TREE and NEWICK
    r###"BCR=123085 COMPLETE TREE NEWICK CDR3=CARDLGGRYYGSKDPW"###,
    // 90. test KEEP_CELL_IF with non-null value
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T
        CDR3=CARGYEDFTMKYGMDVW KEEP_CELL_IF="keeper == 'yes'""###,
    // 91. test FCELL with null value
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T
        CDR3=CARGYEDFTMKYGMDVW FCELL="keeper == ''""###,
    // 92. test NALL_CELL
    r###"BCR=123085 NALL_CELL CDR3=CQKYDSAPLTF MIN_CELLS=20"###,
    // 93. test MIN_DATASET_RATIO
    r###"BCR=123085,123089 MIN_DATASET_RATIO=8 LVARSP=nd2"###,
    // 94. test use of SEG twice
    r###"BCR=123085 SEG=IGHV5-51 SEG=IGKV1D-39"###,
    // 95. test TREE=const
    r###"BCR=123085 TREE=const CDR3=CARPKSDYIIDAFDIW MIN_CELLS=2"###,
    // 96. test MAX_LOG_SCORE
    r###"BCR=123085 CDR3=CARDQNFDESSGYDAFDIW MAX_LOG_SCORE=0.0"###,
    // 97. test MAX_CDR3_DIFFS
    r###"BCR=123085 CDR3=CARESVVGLLPIFDYW MAX_CDR3_DIFFS=1"###,
    // 98. test reduced stringency D alignment
    // (RE can be removed once cellranger rerun)
    r###"TCR=101287 CDR3=CASSPAGTSGKVWGTDTQYF RE"###,
    // 99. test mait (redundant with mait_example.html below, so could delete)
    r###"TCR=101287 LVARSP=mait CDR3=CSAGQGDTEAFF"###,
    // 100. test inkt and INKT
    r###"TCR=101287 LVARSP=inkt INKT MIN_CELLS=2"###,
    // 101. test MAIT
    r###"TCR=101287 LVARSP=mait MAIT MIN_CELLS=50"###,
    // 102. test BINARY with unwriteable path
    r###"BCR=123085 BINARY=/gerbilspam/bumblebee EXPECT_FAIL"###,
    // 103. test POUT without PCOLS (somewhat annoying, because easily triggered to change)
    r###"BCR=85333 POUT=stdout CDR3=CQSADSSGTYKVF"###,
    // 104. test EASY
    r###"BCR=123085 CDR3="CARVIVGPKKLEGRLYSSSLHFDCW|CARVIVGPEKQEGRLYSSSLHFDYW" EASY
        MAX_LOG_SCORE=100"###,
    // 105. test MAX_DEGRADATION and MAX_DIFFS
    r###"BCR=123085,123089 MAX_LOG_SCORE=100 MAX_DEGRADATION=150 MAX_DIFFS=200
        MAX_CDR3_DIFFS=100 CDR3=CVRILGRALTVRVYFYYGIDVW"###,
    // 106. test for failed interaction between POUT and COMPLETE (crashed at one point)
    r###"BCR=123085 CDR3=CAKANQLLYGGRQYYYGMDVW COMPLETE POUT=stdout
        PCOLS=clonotype_id,exact_subclonotype_id,n,d_donor1,d_donor2"###,
    // 107. part 1 of test for weak onesies filter
    r###"TCR=101287 CDR3=CASSQVAGAGQPQHF"###,
    // 108. part 2 of test for weak onesies filter
    r###"TCR=101287 CDR3=CASSQVAGAGQPQHF NWEAK_ONESIES"###,
    // 109. test Levenshtein distance pattern
    r###"BCR=123085 CDR3="CAKDKVPRRSSWSVFDYYGMDVW~9|CAVTIFGVRTALPYYYALDVW~9" NGROUP"###,
    // 110. test dref_aa
    r###"BCR=123085 LVARSP=dref,dref_aa CDR3=CAREKGIGSSGWDWGAFDIW"###,
    // 111. test for fail if F used with unsupported variable
    r###"BCR=123085 LVARSP=near F="near>=0" EXPECT_FAIL"###,
    // 112. test 1 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr1,cdr2,cdr3"###,
    // 113. test 2 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr1,cdr3"###,
    // 114. test 3 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr2,cdr3"###,
    // 115. test 4 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr1,cdr2"###,
    // 116. test 5 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr1"###,
    // 117. test 6 of 6 for cdr1/cdr2 in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr2"###,
    // 118. test cdr1_aa and cdr2_aa
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=cdr1_aa,cdr2_aa AMINO=cdr1"###,
    // 119. test cdr3_aa
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=cdr3_aa AMINO=cdr3"###,
    // 120. test cdr1_dna and cdr2_dna
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARS=cdr1_dna,cdr2_dna AMINO="###,
    // 121. test cdr1_len and cdr2_len
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARS=cdr1_len,cdr2_len AMINO="###,
    // 122. test insertion in CDR1
    r###"BCR=123089 CDR3=CARARPYSSGWSLDAFDIW AMINO=cdr1,cdr3 CVARSP=cdr1_aa"###,
    // 123. test fwr1_dna and fwr2_dna
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=fwr1_dna,fwr2_dna AMINO=cdr3"###,
    // 124. test fwr3_dna
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=fwr3_dna AMINO=cdr3"###,
    // 125. test fwr1_aa and fwr2_aa and fwr3_aa
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=fwr1_aa,fwr2_aa,fwr3_aa AMINO=cdr3"###,
    // 126. test fwr1_len and fwr2_len and fwr3_len
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARSP=fwr1_len,fwr2_len,fwr3_len AMINO=cdr3"###,
    // 127. test 1/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr1,cdr1"###,
    // 128. test 2/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr1,fwr2"###,
    // 129. test 3/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr2,cdr2"###,
    // 130. test 4/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr2,fwr3"###,
    // 131. test 5/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr3,cdr3"###,
    // 132. test 6/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr1,fwr2"###,
    // 133. test 7/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr2,fwr3"###,
    // 134. test 8/8 for fwr* in AMINO
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,fwr2,cdr2,fwr3"###,
    // 135. test CONST_IGH
    r###"BCR=123085 CDR3=CARPKSDYIIDAFDIW SEG=IGLV1-47 CONST_IGH="IGHG.""###,
    // 136. test CONST_IGKL
    r###"BCR=123085 CDR3=CARPKSDYIIDAFDIW SEG=IGLV1-47 CONST_IGKL=IGLC3"###,
    // 137. test 1/2 of fwr4
    r###"BCR=85333 CDR3=CARDLRVEGFDYW AMINO=var,share,donor,cdr3,fwr4 CVARS=fwr4_aa"###,
    // 138. test 2/2 of fwr4
    r###"BCR=85333 CDR3=CARDLRVEGFDYW CVARS=fwr4_dna,fwr4_len"###,
    // 139. test cvar vj_seq_nl
    r###"BCR=85333 CHAINS=1 CDR3=CAAWDDSLNGWVF POUT=stdout PCOLS=vj_seq_nl1"###,
    // 140. test cvar vj_aa_nl
    r###"BCR=85333 CHAINS=1 CDR3=CAAWDDSLNGWVF POUT=stdout PCOLS=vj_aa_nl1"###,
    // 141. test cvar aa%
    r###"BCR=85333 CDR3=CAKGDRTGYSYGGGIFDYW CVARS=aa%,dna%"###,
    // 142. test 1/3 of DIFF_STYLE
    r###"BCR=123085 CDR3=CARVRDILTGDYGMDVW DIFF_STYLE=C1"###,
    // 143. test 2/3 of DIFF_STYLE
    r###"BCR=123085 CDR3=CARVRDILTGDYGMDVW DIFF_STYLE=C2"###,
    // 144. test 3/3 of DIFF_STYLE
    r###"BCR=123085 CDR3=CAREPLYYDFWSAYFDYW DIFF_STYLE=C1"###,
    // 145. test the lead variable "filter"
    r###"BCR=123085 NALL LVARSP=filter PER_CELL CHAINS=2 CDR3=CQQSYSTPPYTF SEG=IGKV1D-39
        SEG=IGLV3-21"###,
    // 146. test BUILT_IN
    r###"BCR=../2.0/124550 CDR3=CAREPLYYDFWSAYFDYW BUILT_IN"###,
    // 147. test NALL_GEX
    r###"BCR=86237 GEX=85679 NALL_GEX LVARSP=n_gex,filter PER_CELL BARCODE=CTTGGCTGTTAAGACA-1"###,
    // 148. test that LVARSP=n_gex fails if only BCR provided
    r###"BCR=1031851 LVARSP=n_gex EXPECT_FAIL"###,
    // 149. test FCELL with complex expression
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=keeper,rank
        FCELL="keeper == 'no' && rank > 10""###,
    // 150. test FCELL with a more complex expression
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=keeper,rank
        FCELL="(keeper == 'no' && rank > 10) || keeper == 'maybe'""###,
    // 151. test PEER_GROUP
    r###"BCR=85333 CDR3=CAKGRYSSPQYYFDYW PEER_GROUP=stdout"###,
    // 152. test PEER_GROUP with PG_READABLE
    r###"BCR=85333 CDR3=CAKGRYSSPQYYFDYW PEER_GROUP=stdout PG_READABLE"###,
    // 153. test d_start and d_frame
    r###"BCR=86237 CDR3=CARGHPNYDYVWGSYRYRAYYFDYW POUT=stdouth
        PCOLS=d_start1,d_frame1,d_start2,d_frame2"###,
    // 154. test POUT=stdout with NOPRINT
    r###"BCR=85333 CDR3="CARTSNRGIVATIFRAFDIW|CARDPRGWGVELLYYMDVW" NOPRINT POUT=stdout
        PCOLS=cdr3_aa1"###,
    // 155. test count_<regex> and F for that
    r###"BCR=123085 LVARSP="z:count_CAKTG" F="z > 0""###,
    // 156. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr1_aa1,fwr1_aa_ref1 AMINO=fwr1"###,
    // 157. test ref variables
    r##"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr1_dna1,fwr1_dna_ref1 AMINO=fwr1"##,
    // 158. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr2_aa1,fwr2_aa_ref1 AMINO=fwr2"###,
    // 159. test ref variables
    r##"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr2_dna1,fwr2_dna_ref1 AMINO=fwr2"##,
    // 160. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr3_aa1,fwr3_aa_ref1 AMINO=fwr3"###,
    // 161. test ref variables
    r##"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr3_dna1,fwr3_dna_ref1 AMINO=fwr3"##,
    // 162. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=fwr4_aa1,fwr4_aa_ref1 AMINO=fwr4"###,
    // 163. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=cdr1_aa2,cdr1_aa_ref2 AMINO=cdr1"###,
    // 164. test ref variables
    r##"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=cdr1_dna2,cdr1_dna_ref2 AMINO=cdr1"##,
    // 165. test ref variables
    r###"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=cdr2_aa2,cdr2_aa_ref2 AMINO=cdr2"###,
    // 166. test ref variables
    r##"BCR=123085 CDR3=CAREVEQWLERNTLDYW POUT=stdouth PCOLS=cdr2_dna2,cdr2_dna_ref2 AMINO=cdr2"##,
    // 167. Test that for TCR, the number of two-chain clonotypes does not change.  It is probably
    // OK for it to change a little bit, but a big change would be indicative of a problem.  At
    // one point we had a release with such a problem and this test is here to prevent that from
    // happening again.
    r###"TCR=101287 NOPRINT REPROD REQUIRED_TWO_CHAIN_CLONOTYPES=849 EXPECT_OK"###,
    // 168. Test POUT without PCELL, where a per-barcode variable is converted into a
    // comma-separated list.
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv POUT=stdout PCOLS=rank
        CDR3=CAKAGPTESGYYVWYFDLW MIN_CELLS=2"###,
    // 169. this crashed at one point because the heavy chain CDR3 computed by cellranger was
    // different than the current one, resulting in an inconsistency
    r###"BCR=85333 CDR3=CQQYNSYSYTF CVARSP=fwr3_aa_ref"###,
    // 170. doublet filter, before
    r###"BCR=123085 CDR3=CAREGGVGVVTATDWYFDLW NDOUBLET"###,
    // 171. doublet filter, after
    r###"BCR=123085 CDR3=CAREGGVGVVTATDWYFDLW"###,
    // 172. this crashed at one point
    r###"META=testx/inputs/test11_meta LVARSP=CD56_ab NOPRINT EXPECT_OK"###,
    // 173. test MIN_UMIS
    r###"BCR=85333 MIN_UMIS=100"###,
    // 174. test METAX, and also the origins printed by this was wrong at one point
    r###"METAX="bcr,origin,donor;toast:86237,c,d;zip:123085,a,b" LVARSP=origins,donors
        POUT=stdouth PCOLS=origins,donors CDR3=CARSFFGDTAMVMFQAFDPW"###,
];

// Test using datasets that are either in the extended public dataset collection, or which are
// not publicly avaiable, or which require samtools.

pub const EXTENDED_TESTS: [&str; 22] = [
    // 1. test that used to crash on a particular barcode; this also gave the wrong
    // answer for an insertion until it was fixed
    r###"BCR=40955 NCELL BARCODE=GCGCAGTCAAAGTGCG-1 AMINO=cdr3 NO_PRE NFORCE"###,
    // 2. tests nd2
    r###"BCR=47199,47200,47212 AMINO=cdr3 NCROSS LVARS=nd2 CDR3=CVKGKSGSFWYYFENW
     NO_PRE NFORCE"###,
    // 3. test sec and mem [requires samtools]
    r###"BCR=123085 GEX=123217 LVARSP=sec,mem CDR3=CVKDRVTGTITELDYW"###,
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
    // If we experience failurs on other lena ids, we can add them to this list.
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
    r###"BCR=123085 GEX=123217 LVARSP=gex,IGHV2-5_g_μ CDR3=CALMGTYCSGDNCYSWFDPW"###,
];

// List of examples on site.

pub const SITE_EXAMPLES: [(&str, &str); 12] = [
    // 1.
    // Do not use NH5 because the bin file is too big for git.
    (
        "pages/auto/clonotype_with_gex.html",
        "BCR=123085 CDR3=CTRDRDLRGATDAFDIW GEX=123217 LVARSP=gex,IGHV3-49_g NUMI \
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
];
