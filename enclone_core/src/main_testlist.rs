// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

pub const TESTS: [&str; 294] = [
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
    // 8. test KEEP_CELL_IF with >= and <=
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=rank
         KEEP_CELL_IF="rank >= 2 && rank <= 3""###,
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
    // 13. check TSV file with BC
    r###"BCR=123085 BC=testx/inputs/123077_cells.tsv PER_CELL LVARSP=T CDR3=CARGYEDFTMKYGMDVW"###,
    // 14. test cdr3_aa_conp
    r###"BCR=123085 CVARSP=cdr3_aa_conp CDR3=CAKTGDLELRYFDWDMDVW"###,
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
    // 22. DUPLICATE, TO REMOVE
    r###"BCR=123089 CDR3=CGTWHSNSKPNWVF MIN_CHAINS_EXACT=3"###,
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
    // 29. tests KEEP_CLONO_IF_CELL_MAX and parenthesized variable in it, SUM and MEAN, use of ≥
    r###"BCR=123085 GEX=123217 LVARSP=IGHV3-7_g,IGHV3-7_g_μ
         KEEP_CLONO_IF_CELL_MAX="(IGHV3-7_g_μ)≥10000.0" MIN_CHAINS=2 SUM MEAN H5"###,
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
    r###"PRE=../enclone-data/big_inputs/version{TEST_FILES_VERSION},.
         META=testx/inputs/test35_meta MIN_CELLS=10 MIN_CHAINS_EXACT=2 NOPRINT PLOT=stdout NO_PRE
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
    // Note that F is deprecated, equals KEEP_CLONO_IF_CELL_MEAN.  Also test ⩾.
    r###"BCR=86237 GEX=85679 LVARSP=IGHV3-7_g_μ F="(IGHV3-7_g_μ)⩾4.5" MIN_CHAINS=2 SUM MEAN
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
    r###"BCR=123085 GEX=123217 CDR3=CARPKSDYIIDAFDIW MIN_CELLS=10 H5
         LVARSP="(IGHV5-51|IGLV1-47)_g_%,IGH.*_g_%,IG(K|L).*_g_%""###,
    // 53. add test for _% with PER_CELL
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 LVARSP="gex,n_gex,JCHAIN_g_%,IG%:IG.*_g_%" CVARS=u_μ,const
         MIN_CHAINS_EXACT=2 CDR3=CAREGGVGVVTATDWYFDLW PER_CELL H5"###,
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
    // 61. test that enclone doesn't crash on CS multi 6.1 pipeline
    r###"BCR_GEX=tiny_multi_CS_6.1 ALLOW_INCONSISTENT EXPECT_OK"###,
    // 62. make sure color from BC can be used as lead variable was broken)
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=color
         AMINO=cdr3 BARCODE=CATATGGTCAGTTGAC-1"###,
    // 63. DUPLICATE, TO REPLACE
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML SET_IN_STONE"###,
    // 64. DUPLICATE, TO REPLACE
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML NGROUP SET_IN_STONE"###,
    // 65. test NCELL
    r###"BCR=86237 NCELL CDR3=CAKTATTLGGYYSHGLDVW MIN_CELLS=2"###,
    // 66. test BC in combination with PER_CELL and PCELL
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T PCELL
         POUT=stdouth PCOLS=barcode,T CDR3=CAKAGPTESGYYVWYFDLW MIN_CELLS=2 H5"###,
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
    // 73. test cdr3_aa_conx
    r###"BCR=123085 CVARSP=cdr3_aa_conx CDR3=CAKTGDLELRYFDWDMDVW"###,
    // 74. this changed after a bug was fixed; the RE can probably be dropped later when we
    // rerun all the datasets
    r###"BCR=123085 RE CDR3=CARGYEDFTMKYGMDVW POUT=stdouth PCOLS=utr_id2"###,
    // 75. this changed after a bug in RE was fixed, and this is in fact testing RE
    r###"BCR=123085 CDR3=CQQSYSTPRTF RE"###,
    // 76. test PLOT_BY_ISOTYPE
    r###"BCR=123085 MIN_CELLS=10 PLOT_BY_ISOTYPE=stdout NOPRINT MIN_CHAINS_EXACT=2"###,
    // 77. test PLOT_BY_ISOTYPE_COLOR
    r###"BCR=123085 MIN_CELLS=10 PLOT_BY_ISOTYPE=stdout NOPRINT MIN_CHAINS_EXACT=2
         PLOT_BY_ISOTYPE_COLOR=red,green,blue,yellow,black,orange,turquoise,pink,gray,purple"###,
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
         PCOLS=barcode,entropy,entropy_cell CDR3=CARAQRHDFWGGYYHYGMDVW H5"###,
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
    r###"BCR=123085 GEX=123217 H5 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T
         CDR3=CARGYEDFTMKYGMDVW KEEP_CELL_IF="keeper == 'yes'""###,
    // 91. test FCELL with null value
    // Do not use NH5 because the bin file is too big for git.
    r###"BCR=123085 GEX=123217 H5 BC=testx/inputs/123077_cells.csv PER_CELL LVARSP=gex,cred,T
         CDR3=CARGYEDFTMKYGMDVW FCELL="keeper == ''""###,
    // 92. test NALL_CELL
    r###"BCR=123085 NALL_CELL CDR3=CQKYDSAPLTF MIN_CELLS=20"###,
    // 93. test MIN_DATASET_RATIO
    r###"BCR=123085,123089 MIN_DATASET_RATIO=6 LVARSP=nd2"###,
    // 94. test use of SEG twice
    r###"BCR=123085 SEG=IGHV5-51 SEG=IGKV1D-39"###,
    // 95. test TREE=const
    r###"BCR=123085 TREE=const CDR3=CARPKSDYIIDAFDIW MIN_CELLS=2"###,
    // 96. test MAX_LOG_SCORE
    r###"BCR=123085 CDR3=CARDQNFDESSGYDAFDIW MAX_LOG_SCORE=0.0"###,
    // 97. Test MAX_CDR3_DIFFS.  This is also an instance where an exact subclonotype has
    // two chains with indentical CDR3s, and this is the right answer, until and unless we change
    // cellranger to somehow not emit two such chains.
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
    // 111. test for fail if F used with unsupported variable (but now supported)
    // Note that F is deprecated, equals KEEP_CLONO_IF_CELL_MEAN.
    r###"BCR=123085 LVARSP=near F="near>=0" EXPECT_OK"###,
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
    // 122. test insertion in CDR1 and test cdr3_start when there is an insertion
    r###"BCR=123089 CDR3=CARARPYSSGWSLDAFDIW AMINO=cdr1,cdr3 CVARSP=cdr1_aa
         POUT=stdout PCOLS=cdr3_start1"###,
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
    // Note that F is deprecated, equals KEEP_CLONO_IF_CELL_MEAN.
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
    // 175. test some variables
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth PCOLS=cdr1_aa1,cdr1_aa_1_2_ext1"###,
    // 176. test some variables
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth PCOLS=cdr2_aa1,cdr2_aa_1_2_ext1"###,
    // 177. test some variables
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth PCOLS=cdr3_aa1,cdr3_aa_1_2_ext1"###,
    // 178. test an ndiff variable as a parseable variable
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth PCOLS=ndiff1vj1"###,
    // 179. test cdr1_aa_north etc.
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth
         PCOLS=cdr1_aa_north1,cdr1_aa_north2,cdr2_aa_north1,cdr2_aa_north2,cdr3_aa_north1,cdr3_aa_north2"###,
    // 180. test some count vars
    r###"BCR=85333 LVARS=all:count_C,c:count_cdr_C,c1:count_cdr1_C,c3:count_cdr3_C,f:count_fwr_C,f1:count_fwr1_C CDR3=CARDKEGLSGYAVERAFDYW"###,
    // 181. test some count vars
    r###"BCR=85333 LVARS=f2:count_fwr2_C CDR3=CVKDIRESSGPLLSHSFDLW"###,
    // 182. test some count vars
    r###"BCR=85333 LVARS=f3:count_fwr3_C CDR3=CARGGFSHAFDIW"###,
    // 183. test some count vars
    r###"BCR=123085 LVARS=f4:count_fwr4_V CDR3=CTRDRDLRGATDAFDIW"###,
    // 184. test some count vars
    r###"BCR=123085 LVARS=c2:count_cdr2_C CDR3=CARQQDVYTRSWYFDYW CELLS=1"###,
    // 185. test SUPPRESS_ISOTYPE_LEGEND
    r###"BCR=123085 MIN_CELLS=10 PLOT_BY_ISOTYPE=stdout NOPRINT MIN_CHAINS_EXACT=2
         SUPPRESS_ISOTYPE_LEGEND"###,
    // 186. test LVAR= (with no value)
    r###"BCR=123085 CDR3=CAREPLYYDFWSAYFDYW LVARS="###,
    // 187. test FOLD_HEADERS
    r###"BCR=123085 CDR3=CAREADYCSGGSCYFSDW FOLD_HEADERS AMINO=cdr3 CVARS=u"###,
    // 188. test for correct handling of COMPLETE + r_cell1 (asserted at one point)
    r###"BCR=85333 POUT=stdout PCOLS=r_cell1 COMPLETE PCELL CDR3=CARGQEGSGWYRPWDYW"###,
    // 189. test CONP
    r###"BCR=123085 CONP CDR3=CVKRASGSAFTAPYW"###,
    // 190. test CONX
    r###"BCR=123085 CONX CDR3=CVKRASGSAFTAPYW"###,
    // 191. test CONP when there's a gap
    r###"BCR=123085 CONP CDR3=CALGGYTWFDPW"###,
    // 192. test INFO
    r###"BCR=123085 CDR3=CAREGGVGVVTATDWYFDLW INFO=testx/inputs/123085_info.csv LVARSP=funny"###,
    // 193. check that this fails gracefully
    r###"NOPRINT EXPECT_FAIL"###,
    // 194. this crashed at one point
    r###"BCR=86237 GEX=85679 LVARSP=g37:IGHV3-7_g_μ NH5 POUT=stdout PCOLS=g37 EXPECT_OK"###,
    // 195. failed at one point
    r###"BCR=86237 GEX=85679 LVARSP=woof:IGHV3-7_g_μ NH5 POUT=stdout PCOLS=woof EXPECT_OK"###,
    // 196. test TREE=n
    r###"BCR=123085 COMPLETE TREE=n CDR3=CARDLGGRYYGSKDPW"###,
    // 197. failed at one point
    r###"BCR=123085 INFO=testx/inputs/123085_info.csv LVARSP=funny EXPECT_OK"###,
    // 198. test TREE=n,cdr2_aa1
    r###"BCR=123085 AMINO=cdr3 CDR3=CAVTIFGVRTALPYYYALDVW TREE=n,cdr2_aa1"###,
    // 199. test KEEP_CLONO_IF_CELL_MEAN with INFO
    r###"BCR=123085 INFO=testx/inputs/123085_info.csv LVARSP=moo
         KEEP_CLONO_IF_CELL_MEAN="moo>0""###,
    // 200. test SCAN_EXACT
    r###"BCR=123085 GEX=123217 LVARSP=IGHV1-69D_g_μ,IGHV3-64D_g_μ MIN_CELLS=10 
         SCAN="(IGHV1-69D_g_μ)>=1800,(IGHV3-64D_g_μ)>=100,t-10*c>=5.0" NOPRINT H5 SCAN_EXACT"###,
    // 201. test SOURCE
    r###"SOURCE=testx/inputs/123085_args AMINO=cdr2,cdr3"###,
    // 202. DUPLICATE TO REPLACE
    r###"SOURCE=testx/inputs/123085_args AMINO=cdr2,cdr3 EXPECT_OK"###,
    // 203. test plotting with using the BC option to set color
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv PLOT=stdout NOPRINT"###,
    //
    // TESTS WITH PER_CELL AND PCELL
    //
    // 204. test INFO with PER_CELL and PCELL
    r###"BCR=123085 CDR3=CAREGGVGVVTATDWYFDLW INFO=testx/inputs/123085_info.csv POUT=stdout
         PCOLS=moo LVARS=moo PCELL PER_CELL"###,
    // 205. test g<d> with PER_CELL and PCELL
    r###"BCR=123085 GEX=123217 AMINO=cdr3 LVARS=g15 CDR3=CARVRDILTGDYGMDVW POUT=stdout PCOLS=g15
         PCELL PER_CELL H5"###,
    // 206. test origins with PER_CELL and PCELL
    r###"BCR=123085:123089 AMINO= CDR3=CTRAGFLSYQLLSYYYYGMDVW FOLD_HEADERS POUT=stdout PCELL
         PER_CELL PCOLS=origins LVARSP=origins"###,
    // 207. test datasets with PER_CELL and PCELL
    r###"BCR=123085:123089 CELLS=5 AMINO= CDR3=CTRAGFLSYQLLSYYYYGMDVW FOLD_HEADERS POUT=stdout
         PCELL PER_CELL PCOLS=datasets LVARSP=datasets"###,
    // 208. test donors with PER_CELL and PCELL
    r###"BCR="123085;123089" AMINO= CDR3=CTRAGFLSYQLLSYYYYGMDVW FOLD_HEADERS POUT=stdout PCELL
         PER_CELL PCOLS=donors LVARS=donors MIX_DONORS CHAINS=2"###,
    // 209. test n with PER_CELL and PCELL
    r###"BCR=123085 AMINO=cdr3 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdout PCELL PER_CELL PCOLS=n"###,
    // 210. test filter with PER_CELL and PCELL
    r###"BCR=123085 AMINO=cdr3 FOLD_HEADERS POUT=stdout PCELL PER_CELL PCOLS=filter LVARSP=filter
         NALL_CELL CDR3=CAKHQRGGGRQNYYYGMDVW"###,
    // 211. test inkt with PER_CELL and PCELL
    r###"TCR=101287 INKT MIN_CELLS=2 AMINO=cdr3 FOLD_HEADERS POUT=stdout PCELL PER_CELL
         PCOLS=inkt LVARSP=inkt"###,
    // 212. test mait with PER_CELL and PCELL
    r###"TCR=101287 AMINO=cdr3 FOLD_HEADERS POUT=stdout PCELL PER_CELL PCOLS=mait LVARSP=mait
         CDR3=CSAGQGDTEAFF"###,
    // 213. test cred with PER_CELL and PCELL
    r###"BCR=123085 GEX=123217 AMINO=cdr3 LVARS=cred CVARS=u POUT=stdout PCOLS=cred,cred_cell
         PCELL PER_CELL H5 CDR3=CARDPEDIVLMVYAMGGNYGMDVW"###,
    // 214. test n_<name> with PER_CELL and PCELL
    r###"BCR=123085:123089 AMINO=cdr3 FOLD_HEADERS POUT=stdout PCELL PER_CELL PCOLS=n_s1
         LVARS=datasets,n_s1 CDR3=CARDLFVLVPAAITYYYGMDVW CVARS=u"###,
    // 215. test n_gex with PER_CELL and PCELL
    r###"BCR=123085 GEX=123217 AMINO=cdr3 LVARS=n_gex POUT=stdout PCOLS=n_gex,n_gex_cell PCELL 
         PER_CELL H5 CDR3=CAKDKVPRRSSWSVFDYYGMDVW"###,
    // 216. test near with PER_CELL and PCELL
    r###"BCR=123085 AMINO=cdr3 POUT=stdout PCELL PER_CELL LVARSP=near PCOLS=near CVARS=u
         CDR3=CARHLQWELPYW"###,
    // 217. test far with PER_CELL and PCELL
    r###"BCR=123085 AMINO=cdr3 POUT=stdout PCELL PER_CELL LVARSP=far PCOLS=far CVARS=u
         CDR3=CARHLQWELPYW"###,
    // 218. test dref with PER_CELL and PCELL
    r###"BCR=123085 AMINO=cdr3 POUT=stdout PCELL PER_CELL LVARSP=dref PCOLS=dref CVARS=u
         CDR3=CSRVFGNSTYYSSRVGGYW"###,
    // 219. test count_cdr_C with PER_CELL and PCELL
    r###"BCR=85333 LVARSP=count_cdr_C CDR3=CARDKEGLSGYAVERAFDYW POUT=stdout PCELL PER_CELL
         PCOLS=count_cdr_C CVARS=u"###,
    // 220. test cdr3_aa_conp with PER_CELL and PCELL
    r###"BCR=123085 AMINO= CDR3=CARHLQWELPYW FOLD_HEADERS POUT=stdout PCELL PER_CELL
         PCOLS=cdr3_aa_conp2 CVARS=cdr3_aa_conp"###,
    // 221. test RPS27_g with PER_CELL and PCELL
    r###"BCR=123085 GEX=123217 AMINO=cdr3 POUT=stdout PCOLS=RPS27_g,RPS27_g_cell PCELL
         PER_CELL H5 CDR3=CAREVEQWLERNTLDYW LVARSP=RPS27_g"###,
    //
    // OTHER TESTS
    //
    // 222. test for busted reference
    r###"BCR=85333 REF=testx/inputs/busted_regions.fa EXPECT_FAIL"###,
    // 223. test {v,d,j}_name and _id
    r###"BCR=86237 CDR3=CARGHPNYDYVWGSYRYRAYYFDYW POUT=stdouth
         PCOLS=v_name1,d_name1,j_name1,v_id1,d_id1,j_id1"###,
    // 224. test const_id and utr_name
    r###"BCR=86237 CDR3=CARSFFGDTAMVMFQAFDPW POUT=stdouth
         PCOLS=const_id1,utr_name1"###,
    // 225. test q<n>_
    r###"BCR=123085 CDR3=CANFGRGGDVAFDIW CVARS=q10_"###,
    //
    // MORE TESTS OF PER_CELL AND PCELL
    //
    // 226. test RPS27_g_mean with PER_CELL and PCELL
    r###"BCR=123085 GEX=123217 AMINO=cdr3 POUT=stdout PCOLS=RPS27_g_mean PCELL PER_CELL H5
         CDR3=CAREVEQWLERNTLDYW LVARSP=RPS27_g_mean CVARS=u"###,
    // 227. test datasets, donors, origins with PER_CELL and PCELL
    r###"BCR=123085:123089 CELLS=5 AMINO= CDR3=CTRAGFLSYQLLSYYYYGMDVW FOLD_HEADERS POUT=stdout
         PCELL PER_CELL PCOLS=datasets,datasets_cell,origins,origins_cell,donors,donors_cell
         LVARSP=origins,donors"###,
    // 228. test clonotype_ncells with PER_CELL and PCELL
    r###"BCR=123085 AMINO= CDR3=CARHLQWELPYW FOLD_HEADERS POUT=stdout PCELL PER_CELL
         PCOLS=clonotype_ncells LVARSP=clonotype_ncells"###,
    //
    // OTHER TESTS
    //
    // 229. test KEEP_CLONO_IF_CELL_MAX with comp
    r###"BCR=123085 CVARSP=comp KEEP_CLONO_IF_CELL_MAX="comp1 >= 18" AMINO=cdr3"###,
    // 230. not really clear what this is doing, but don't delete, as it used to represent
    // strange behavior
    r###"BCR=123085 CDR3=CTRDRDLRGATDAFDIW"###,
    // 231. test ≤
    r###"BCR=86237 KEEP_CLONO_IF_CELL_MEAN="u2≤150" NOPRINT SUMMARY SUMMARY_CLEAN"###,
    // 232. test nonsense variable in linear constraint
    r###"BCR=86237 KEEP_CLONO_IF_CELL_MAX="gexzz > 8000" EXPECT_FAIL H5"###,
    // 233. test use of two linear constraints
    r###"BCR=123085 GEX=123217
         KEEP_CLONO_IF_CELL_MAX="gex > 8000" KEEP_CLONO_IF_CELL_MAX="gex < 8200" H5"###,
    // 234. test tooltip comments; this is via a testing-only filename option gui_stdout
    r###"BCR=123085 MIN_CELLS=10 PLOT_BY_ISOTYPE=gui_stdout NOPRINT MIN_CHAINS_EXACT=2"###,
    // 235. test that v_name etc. do not appear in parseable output if chain is absent
    r###"BCR=123085 CDR3=CVRGLRTW PCOLS=barcodes,v_name1,j_name1,v_id1,j_id1 POUT=stdouth"###,
    // 236. test MAX_HEAVIES=1
    r###"BCR=123085 CDR3=CASPVPYYYDSSGYPYW MAX_HEAVIES=1 EXPECT_NULL"###,
    // 237. test enclone --help
    r###"--help NO_PRE EXPECT_OK"###,
    // 238. test cigar
    r###"BCR=123085 AMINO=cdr3 POUT=stdout PCOLS=cigar2 CDR3=CTRSSTTPRDPTMIVVAYYYYGMDVW"###,
    // 239. test group post filtering
    r###"BCR=123085 G=2,100-101 NGROUP"###,
    // 240. test sym option for PLOTXY_EXACT
    r###"BCR=123085 PLOTXY_EXACT=u1,u2,stdout,sym NOPRINT"###,
    // 241. test KEEP_CELL_IF on gex var
    r###"BCR=123085 GEX=123217 LVARSP=IGHM_g KEEP_CELL_IF="IGHM_g>=10" CDR3=CARRYFGVVADAFDIW H5"###,
    // 242. test nchains_present
    r###"BCR=86237 LVARSP=nchains_present CDR3=CARSFFGDTAMVMFQAFDPW"###,
    // 243. this crashed at one point
    r###"BCR=123085 GROUP="cdr3_aa_heavy>=85%,vj_refname" MIN_GROUP=2 PLOT=/dev/null
         NOPRINT EXPECT_OK"###,
    // 244. test for very long (120 amino acid) CDR3
    // Note that this long CDR3 is likely part of a nonproductive chain.  The test is here because
    // there may be long productive CDR3 sequences in data from other species, although we do not
    // have such data.  This is from 1020665.
    r###"BCR=testx/inputs/flaky BUILT_IN REPROD CVARSP=cdr3_len CDR3=CARDGGGQPFDLW AMINO="###,
    // 245. Test a tweak to the weak chains filter.  This should have two chains.  From 174957.
    r###"BCR=testx/inputs/flaky2 CDR3=CARPRGYCSGGSCFPFASW BUILT_IN"###,
    // 246. test that used to crash on a particular barcode; this also gave the wrong
    // answer for an insertion until it was fixed
    r###"BCR=testx/inputs/flaky3 NCELL CDR3=CARNWRYCTSVSCQHREYFYYMDVW AMINO=cdr3"###,
    // 247. this crashed
    r###"BCR=testx/inputs/flaky4"###,
    // 248. this crashed
    r###"BCR=testx/inputs/flaky5"###,
    // 249. an example that triggered an internal inconsistency test, which we subsequently removed;
    // there are three chains and the middle one was the problem
    r###"TCR=testx/inputs/flaky6 BARCODE=CCAGCGAAGTGTTGAA-1 REPROD EXPECT_OK"###,
    // 250. test MOUSE + IMGT; note that specifying by number forces BCR+TCR reference checks
    // Added FORCE_EXTERNAL because couldn't reproduce the result.  Don't understand.
    r###"74396 MOUSE NOPRINT SUMMARY SUMMARY_CLEAN IMGT ACCEPT_BROKEN FORCE_EXTERNAL"###,
    // 251. test mouse + IMGT; note that specifying by number forces BCR+TCR reference checks
    r###"74396 MOUSE REQUIRE_UNBROKEN_OK IMGT ACCEPT_BROKEN EXPECT_NULL"###,
    // 252. this exhibits what happens when signature filtering is ON, see next
    // this was the only example we could find
    // based on 83808-83809, derived using modified version of minimal_fail, and also shrink_json
    r###"BCR=testx/inputs/flaky7 BUILT_IN REPROD REQUIRED_TWO_CHAIN_CLONOTYPES=1
         REQUIRED_THREE_CHAIN_CLONOTYPES=0 NOPRINT EXPECT_OK"###,
    // 253. this exhibits what happens when signature filtering is OFF, see previous
    // this was the only example we could find
    // based on 83808-83809, derived using modified version of minimal_fail, and also shrink_json
    r###"BCR=testx/inputs/flaky7 BUILT_IN REPROD NSIG REQUIRED_TWO_CHAIN_CLONOTYPES=0
         REQUIRED_THREE_CHAIN_CLONOTYPES=1 NOPRINT EXPECT_OK"###,
    // 254. parseable value for fwr4_aa was wrong, from 1117070
    r###"BCR=testx/inputs/flaky AMINO=fwr4 CDR3=CAKDVNGYSSGWAFENW POUT=stdout PCOLS=fwr4_aa1
         BUILT_IN"###,
    // 255. conp value was truncated, from 1117069
    r###"BCR=testx/inputs/flaky CONP CDR3=CVRDPPEELELFDYW BUILT_IN"###,
    // 256. Make sure that FP join output includes join error details.
    // If somehow we fix the FP join occurring here, another one should be substituted.
    // This is from BCR="131036;140707".
    r###"PRE=testx/inputs BCR="flaky8a;flaky8b" ANN SHOW_BC MIN_DONORS=2
         PRINT_FAILED_JOINS BUILT_IN NO_PRE"###,
    // 257. clonotype that was two clonotypes before raising MAX_DIFFS to 60, from 1084461-1084462
    r###"BCR=testx/inputs/flaky CDR3=CAKEFGNGGFDTFDIW BUILT_IN AMINO=cdr3"###,
    // 258. This used to appear as a four-chain clonotype, and is now split.  From 123085,123090.
    r###"BCR=testx/inputs/flaky9 BUILT_IN REQUIRED_FOUR_CHAIN_CLONOTYPES=0 EXPECT_OK"###,
    // 259. this crashed at one point, from 83809
    r###"BCR=testx/inputs/flaky10 EXPECT_OK"###,
    // 260. the result of this changed when sub_alts was changed, from 40086;132888
    r###"BCR=testx/inputs/flaky11 MAX_DIFFS=80 CDR3=CVKGDWGSAFDIW BUILT_IN"###,
    // 261. previously this yielded a disconnected clonotype, from 140699,140705-140706
    r###"BCR=testx/inputs/flaky12 AMINO=cdr3 CDR3="CAKDRQAGGIGEVDDW|CARDRVPGGIGEVDYW" BUILT_IN"###,
    // 262. test NSEG
    r###"BCR=86237 SEG=IGHV4-59 NSEG="IGHJ3|IGHJ4|IGHJ6""###,
    // 263. test NSEGN
    r###"BCR=86237 SEG=IGHV4-34 NSEGN="51|54|55|57|321""###,
    // 264. test MIN_ORIGINS
    r###"BCR=123085:123089 MAX_CELLS=2 SEG=IGHV3-49 MIN_ORIGINS=2"###,
    // 265. test DVARS
    // The output is a bit flaky because we imported some but not all of the special files
    // for 85679.
    r###"BCR=86237 GEX=85679 DVARS=CD19_ab_cellular_u,CD19_ab_cellular_r
         NOPRINT SUMMARY SUMMARY_CLEAN"###,
    // 266. a test of VAR_DEF
    r###"BCR=86237 GEX=85679 VAR_DEF="mu:CD19_ab + CD25_ab" LVARSP=gex,CD19_ab,CD25_ab,mu
         CDR3=CARSFFGDTAMVMFQAFDPW FOLD_HEADERS PER_CELL AMINO="###,
    // 267. a test of VAR_DEF
    r###"BCR=86237 GEX=85679 VAR_DEF=x19:CD19_ab VAR_DEF=x25:CD25_ab VAR_DEF="mu:x19 + x25"
         LVARSP=gex,CD19_ab,CD25_ab,mu CDR3=CARSFFGDTAMVMFQAFDPW FOLD_HEADERS PER_CELL AMINO="###,
    // 268. a test of VAR_DEF
    r###"BCR=86237 GEX=85679 VAR_DEF="pink:PINK1-AS_g" LVARSP=pink CDR3=CARSFFGDTAMVMFQAFDPW
         FOLD_HEADERS PER_CELL AMINO="###,
    // 269. test fb variables
    r###"BCR=86237 GEX=85679 ALLOW_INCONSISTENT NGEX LVARSP=fb1,fb1_n PER_CELL AMINO=cdr3 CVARS=             FOLD_HEADERS POUT=stdouth PCOLS=fb2,fb2_n,fb2_n_cell PCELL CDR3=CARSFFGDTAMVMFQAFDPW"###,
    // 270. test NOSPACES
    r###"BCR=123085 CDR3=CTRDRDLRGATDAFDIW AMINO=cdr3,fwr4 NOSPACES CONX"###,
    // 271. test for weird path bug
    r###"BCR_GEX=tiny_multi_PD_broken EXPECT_OK"###,
    // 272. a test for validated UMI variables
    r###"BCR=tiny_multi_PD CVARS=u,nval,nnval,nival BARCODE=AAAGCAAGTGGCTCCA-1 AMINO= PER_CELL
         POUT=stdouth PCOLS=nval1,nval2,nval3,valumis3,valbcumis2"###,
    // 273. a test for validated UMI variables
    r###"BCR=tiny_multi_PD CVARS=u,nval,nnval,nival AMINO= PER_CELL POUT=stdouth
         PCOLS=ivalumis1,ivalbcumis1,nvalbcumis2 BARCODE=TACCTTAAGAGCCCAA-1"###,
    // 274. at one point this printed bell characters
    r###"CVARS=u,nval,nnval,nival AMINO= PER_CELL POUT=stdouth PCOLS=nval1,nval2 BCR=123085
         BARCODE=ACAGCCGAGATAGGAG-1"###,
    // 275. test _ext var with negative extensions
    r###"BCR=123085 CDR3=CAKDKVPRRSSWSVFDYYGMDVW POUT=stdouth
         PCOLS=cdr3_aa1,cdr3_aa_-1_-2_ext1"###,
    // 276. this failed at one time
    r###"BCR=40970_subset NCELL NOPRINT EXPECT_OK"###,
    // 277. test nbc
    r###"BCR=85333 CDR3=CARDGMTTVTTTAYYGMDVW LVARSP=nbc PER_CELL CVARS= FOLD_HEADERS"###,
    // 278. test some count_fwr variables
    r###"BCR=85333 LVARS=count_fwr1_C,count_fwr_C CDR3=CARGGFSHAFDIW AMINO=cdr3"###,
    // 279. barcode having five contigs
    r###"BCR=123085 NALL CDR3=CAKKHYRYYDSSGYNPLGYYYYGMDVW CVARS=u AMINO= FOLD_HEADERS"###,
    // 280. this asserted at one point
    r###"BCR=86237 CDR3=CARRGPRFRPRFLRGRKTGNWFDPW CVARS= AMINO= EXPECT_FAIL"###,
    // 281. this yielded the wrong aa_nl_2 value
    r###"BCR=123085 CDR3=CARHPAPNYGFWSGYYKTDNWFDPW POUT=stdout PCOLS=vj_aa_nl2 CVARS=u,notes
         AMINO=fwr1"###,
    // 282. fwr3_aa1 was wrong
    r###"BCR=123085 CDR3=CALGGYTWFDPW POUT=stdout PCOLS=fwr3_aa1"###,
    // 283. another test of VAR_DEF
    r###"BCR=86237 GEX=85679 VAR_DEF="mu:CD19_ab + CD25_ab" LVARSP=gex,CD19_ab,CD25_ab,mu
         CDR3=CARSFFGDTAMVMFQAFDPW FOLD_HEADERS PER_CELL AMINO= POUT=stdout PCOLS=mu"###,
    // 284. another test of VAR_DEF
    r###"BCR=86237 GEX=85679 VAR_DEF="mu:CD19_ab + CD25_ab" LVARSP=mu CDR3=CARSFFGDTAMVMFQAFDPW
         PER_CELL POUT=stdout PCOLS=mu PCELL FOLD_HEADERS"###,
    // 285. another test of VAR_DEF
    r###"BCR=123085 VAR_DEF=x:u1 LVARSP=x CDR3=CAKDGYSSSWYVVDW SEG=IGHV3-30"###,
    // 286. this asserted at one point
    r###"BUILT_IN BCR=testx/inputs/flaky2/outs/,testx/inputs/flaky3/outs/ EXPECT_OK"###,
    // 287. test gamma delta data (pos control)
    r###"TCRGD=testx/inputs/gamma_delta1 GAMMA_DELTA MOUSE BUILT_IN"###,
    // 288. test gamma delta data without GAMMA_DELTA tag in regular TCR pipe (neg control)
    r###"TCR=testx/inputs/gamma_delta1 MOUSE BUILT_IN REQUIRED_CLONOTYPES=0 EXPECT_NULL"###,
    // 289. test dref_max
    r###"BCR=86237 LVARSP=dref_max CDR3=CARAPEDTSRWPQYNYSGLDVW AMINO=cdr3"###,
    // 290. test v_name_orig
    r###"BCR=123089 CVARS=v_name_orig PCELL POUT=stdout PCOLS=v_name_orig_cell2 PER_CELL
         CDR3=CARDRIDDSSGYYYAYYYGMDVW"###,
    // 291. test SPLIT_PLOT_BY_DATASET
    r###"BCR=123085,123089 PLOT_BY_ISOTYPE=stdout SPLIT_PLOT_BY_DATASET NOPRINT"###,
    // 292. this asserted
    r###"BCR=85333 CDR3=”CAKGDRTGYSYGGGIFDYW~3” NOPRINT SUMMARY EXPECT_FAIL"###,
    // 293. test BC var in color by variable
    r###"BCR=123085 BC=testx/inputs/123077_cells.csv KEEP_CELL_IF="rank >= 1"
         HONEY="out=stdout,color=var,rank" NOPRINT"###,
    // 294. this asserted
    r###"BCR=86237 HONEY=out=stdout,color=var,cdr3_aa1 NOPRINT EXPECT_FAIL"###,
];
