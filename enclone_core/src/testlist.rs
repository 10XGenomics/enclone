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
    r###"BCR=123085 GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT"###,
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

pub const TESTS: [&str; 237] = [
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
    // 230. weird but correct result, becaused filtering is applied simultaneously with other
    // filters
    r###"BCR=123085 KEEP_CLONO_IF_CELL_MAX="nchains > 2" CDR3=CTRDRDLRGATDAFDIW"###,
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
    "AMINO=fwr1,cdr1,fwr2,cdr2,fwr3,cdr3,fwr4 CVARS=d1_name,d2_name,d_delta,d_Δ",
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

pub const EXTENDED_TESTS: [&str; 30] = [
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

pub const SITE_EXAMPLES: [(&str, &str); 22] = [
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
        "BCR=123085 GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT \
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

// Notes on how to add to the above SITE_EXAMPLES:
// 1. cargo b
// 2. merge_html BUILD
// 3. ./build

];
