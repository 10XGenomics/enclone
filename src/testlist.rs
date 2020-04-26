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

pub const TEST_FILES_VERSION: u8 = 14;

pub const TESTS: [&str; 64] = [
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
    // 6. tests AMINO= and vjlen
    r###"BCR=86237 CELLS=3 AMINO= CVARS=u,r,cdr3_dna,vjlen"###,
    // 7. tests SHM deletion
    r###"BCR=123085 CVARSP=var,clen,cdiff CDR3=CAREPLYYDFWSAYFDYW LVARSP=near,far"###,
    // 8. this clonotype included a junk chain before we made a change, and test "/outs"
    r###"TCR=163911/outs CDR3=CAPSAGDKIIF AMINO=donor"###,
    // 9. tests PER_CELL
    r###"BCR=85333 CDR3=CAKGDRTGYSYGGGIFDYW PER_CELL"###,
    // 10. tests multiple datasets and also LVARS=n,samples,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    r###"BCR="123085;123089" CDR3=CVKDRVTGTITELDYW LVARS=n,samples,donors,datasets AMINO=share
     MIX_DONORS"###,
    // 11. tests META
    r###"META=test/inputs/test11_meta CDR3=CARSFFGDTAMVMFQAFDPW LVARSP=donors,gex"###,
    // 12. this added because it got better when a noise filter was added, also tests u_max
    r###"TCR=163914 CDR3=CASSLVQPSTDTQYF CVARSP=u_max"###,
    // 13. this added because it got better when a noise filter was added; also test FASTA
    r###"TCR=163914 CDR3=CAFRGGSYIPTF FASTA=stdout"###,
    // 14. this added because it got better when a bug in bads detection was fixed
    r###"TCR=163914 CDR3=CASRLGGEETQYF"###,
    // 15. tests insertion and AMINO range
    r###"BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113"###,
    // 16. tests number of cells broken out by dataset
    r###"BCR=123085,123089 LVARS=n,n_123085,n_123089 CDR3=CTRDRDLRGATDAFDIW"###,
    // 17. tests gex with PER_CELL and tests n_gex
    // See also enclone_test_prebuild below, that tests nearly the same thing,
    // and tests versus the same output file.
    r###"BCR=86237 GEX=85679 LVARSP=gex_max,gex,n_gex,CD19_ab_μ CELLS=3 PER_CELL NH5"###,
    // 18. makes sure cross filtering is isn't applied to two samples from same donor
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
    // 29. tests parenthesized variable in F, SUM and MEAN
    r###"BCR=86237 GEX=85679 LVARSP=IGHV3-7_g_μ F="(IGHV3-7_g_μ)>=4.5" MIN_CHAINS=2 SUM MEAN
     NH5"###,
    // 30. tests d_univ and d_donor
    r###"BCR=123085 CVARSP=d_univ,d_donor CDR3=CVKDRVTGTITELDYW"###,
    // 31. tests Cell Ranger 3.1 output
    r###"BCR=../3.1/123085 CDR3=CVKDRVTGTITELDYW"###,
    // 32. tests Cell Ranger 2.0 output and RE
    r###"BCR=../2.0/124550 CDR3=CAREPLYYDFWSAYFDYW RE"###,
    // 33. tests SCAN
    r###"BCR=123085 GEX=123749 LVARSP=IGHV1-69D_g_μ MIN_CELLS=10 NGEX
     SCAN="(IGHV1-69D_g_μ)>=100,(IGHV1-69D_g_μ)<=1,t-10*c>=0.1" NOPRINT H5"###,
    // 34. tests honeycomb plot
    // (This yields a lot of output so will be annoying to debug if something changes.)
    r###"BCR=123085:123089 MIN_CELLS=10 PLOT="stdout,s1->red,s2->blue" NOPRINT
     LEGEND=red,"cell from 123085",blue,"cell from 123089""###,
    // 35. tests barcode-by-barcode specification of colors, and tests LEGEND=
    // Note that the specification of PRE overrides our usual specification.
    // (This yields a lot of output so will be annoying to debug if something changes.)
    r###"PRE= META=test/inputs/test35_meta MIN_CELLS=10 MIN_CHAINS_EXACT=2 NOPRINT PLOT=stdout
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
    // 40. test case where digit rows are just barely present
    r###"TCR=163911 CDR3=CASSLVQPSTDTQYF AMINO=donor"###,
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
    r###"BCR=123085 GEX=123749 LVARSP=gex_mean,gex_Σ CDR3=CASRKSGNYIIYW NGEX H5"###,
    // 49. test HTML
    r###"BCR=85333 CDR3=CAAWDDSLNGWVF CHAINS=1 POUT=stdouth PCOLS=barcodes,n FASTA=stdout
        FASTA_AA=stdout HTML"###,
    // 50. make sure this doesn't fail
    r###"NOPAGER EXPECT_OK"###,
    // 51. make sure this fails gracefully
    r###"BCR=123085 PLOT=/nonexistent/broken.svg NOPRINT MIN_CELLS=50 EXPECT_FAIL"###,
    // 52. add test for some gene patterns
    r###"BCR=123085 GEX=123749 H5 CDR3=CARPKSDYIIDAFDIW MIN_CELLS=10
        LVARSP="(IGHV5-51|IGLV1-47)_g_%,IGH.*_g_%,IG(K|L).*_g_%""###,
    // 53. add test for _% with PER_CELL
    r###"BCR=123085 GEX=123749 LVARSP="gex,n_gex,JCHAIN_g_%,IG%:IG.*_g_%" CVARS=u_μ,const
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
    // 61. test 5/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN SET_IN_STONE"###,
    // 62. test 6/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN NGROUP SET_IN_STONE"###,
    // 63. test 7/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML SET_IN_STONE"###,
    // 64. test 8/8 for newline correctness
    r###"BCR=85333 GROUP_VJ_REFNAME MIN_GROUP=2 AMINO= PLAIN HTML NGROUP SET_IN_STONE"###,
];

// List of examples in documentation.

pub const EXAMPLES: [&str; 2] = [
    // 1.
    r###"BCR=123089 CDR3=CARRYFGVVADAFDIW"###,
    // 2.
    r###"BCR=123085 GEX=123749 LVARSP=gex,IGHV2-5_g_μ,CD4_ab_μ CDR3=CALMGTYCSGDNCYSWFDPW"###,
];

// List of examples on site.

pub const SITE_EXAMPLES: [(&str, &str); 6] = [
    // 1.
    (
        "clonotype_with_gex",
        "BCR=123085 CDR3=CQQRSNWPPSITF GEX=123749 LVARSP=gex,IGHV3-49_g,CD19_ab",
    ),
    // 2.
    (
        "illusory1",
        "BCR=128037,128040 NCROSS CDR3=CARGGTTTYFISW NGROUP",
    ),
    // 3.
    ("illusory2", "BCR=128037,128040 CDR3=CARGGTTTYFISW NGROUP"),
    // 4.
    (
        "illusory3",
        "BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW NGROUP",
    ),
    // 5.
    (
        "illusory4",
        "BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW PER_CELL LVARSP=gex,cred MIN_CHAINS_EXACT=2 \
         NGROUP",
    ),
    // 6.
    (
        "illusory5",
        "BCR=128040 GEX=127801 BC=128024_cells.csv CDR3=CARGGTTTYFISW PER_CELL \
         LVARSP=gex,cred,T CHAINS_EXACT=2 NGROUP",
    ),
];
