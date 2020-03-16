// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

pub const TEST_FILES_VERSION: u8 = 14;

// Main list of tests for enclone.

pub const TESTS: [&str; 39] = [
    // 1. tests variant base after CDR3, parseable output
    r###"BCR=123089 CDR3=CVRDRQYYFDYW POUT=stdout
     PCOLS=exact_subclonotype_id,ncells,v_name1,v_name2,nchains,var_indices_aa1,barcodes"###,
    // 2. tests many donor ref differences, test comp and var and donorn
    r###"BCR=123089 CDR3=CARRYFGVVADAFDIW CVARSP=comp,var AMINO=cdr3,var,share,donorn"###,
    // 3. tests motif in CDR3, CHAINS, u_sum, flipped args in CVARS, on tiny dataset
    r###"BCR=85333 CDR3="CAA.*" CHAINS=2 CVARS=const,u_sum"###,
    // 4. tests gex and antibody, FULL_SEQC, ulen, udiff, on tiny dataset
    r###"BCR=86237 GEX=85679 LVARSP=gex_med,CD19_ab_μ,CD25_ab_μ,IGLV3-1_g_μ,RPS27_g_μ
     CELLS=3 FULL_SEQC
     CVARSP=ulen,udiff"###,
    // 5. tests TCR and correct grouping of onesies on AGBT Donor 2 dataset
    r###"TCR=101287 MIN_CELLS=100"###,
    // 6. tests AMINO=
    r###"BCR=86237 CELLS=3 AMINO= CVARS=u_med,r_med,cdr3_dna"###,
    // 7. tests SHM deletion
    r###"BCR=123085 CVARSP=var,clen,cdiff CDR3=CAREPLYYDFWSAYFDYW LVARSP=near,far"###,
    // 8. this clonotype included a junk chain before we made a change
    r###"TCR=163911 CDR3=CAPSAGDKIIF AMINO=donor"###,
    // 9. tests PER_CELL
    r###"BCR=85333 CDR3=CAKGDRTGYSYGGGIFDYW PER_CELL"###,
    // 10. tests multiple datasets and also LVARS=ncells,samples,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    r###"BCR="123085;123089" CDR3=CVKDRVTGTITELDYW LVARS=ncells,samples,donors,datasets AMINO=share
     MIX_DONORS"###,
    // 11. tests META
    r###"META=test/inputs/test11_meta CDR3=CARSFFGDTAMVMFQAFDPW LVARSP=donors,gex_med"###,
    // 12. this added because it got better when a noise filter was added, also tests u_max
    r###"TCR=163914 CDR3=CASSLVQPSTDTQYF CVARSP=u_max"###,
    // 13. this added because it got better when a noise filter was added; also test FASTA
    r###"TCR=163914 CDR3=CAFRGGSYIPTF FASTA=stdout"###,
    // 14. this added because it got better when a bug in bads detection was fixed
    r###"TCR=163914 CDR3=CASRLGGEETQYF"###,
    // 15. tests insertion and AMINO range
    r###"BCR=86233 CDR3=CARGLVVVYAIFDYW CVARS=notes AMINO=cdr3,105-113"###,
    // 16. tests number of cells broken out by dataset
    r###"BCR=123085,123089 LVARS=ncells,n_123085,n_123089 CDR3=CTRDRDLRGATDAFDIW"###,
    // 17. tests gex with PER_CELL and tests n_gex
    // See also enclone_test_prebuild below, that tests nearly the same thing,
    // and tests versus the same output file.
    r###"BCR=86237 GEX=85679 LVARSP=gex_max,gex,n_gex,CD19_ab_μ CELLS=3 PER_CELL"###,
    // 18. makes sure cross filtering is isn't applied to two samples from same donor
    r###"BCR=123085:123089 CDR3=CVRDEGGARPNKWNYEGAFDIW"###,
    // 19. there was a bug that caused twosie to be deleted, and there was foursie junk
    r###"BCR=123085 CDR3=CARRYFGVVADAFDIW"###,
    // 20. example affected by whitelist (gel bead oligo contamination) filtering, and test u_Σ
    r###"BCR=52177 AMINO=cdr3 PER_CELL CDR3=CATWDDSLSGPNWVF CVARSP=u_Σ"###,
    // 21. test MIN_CHAINS_EXACT
    r###"BCR=123089 CDR3=CGTWHSNSKPNWVF MIN_CHAINS_EXACT=3"###,
    // 22. there was a false positive clonotype
    r###"BCR="165807;165808" FAIL_ONLY=true EXPECT_NULL"###,
    // 23. here we were generating a fake alternate allele
    r###"BCR=83808 CDR3=CAREGRGMVTTNPFDYW MIN_CELLS_EXACT=30"###,
    // 24. an example that uses IGHE
    r###"BCR=52177 CDR3=CSTGWGLDFDFWSGYYTAGYHW"###,
    // 25. add mouse B6 example that had messed up constant regions
    r###"TCR=74396 MOUSE CVARSP=cdiff CDR3=CASSDAGDTQYF"###,
    // 26. tests multiple datasets and also LVARS=ncells,donors,datasets, and share
    // Note that we have deliberately "faked" two donors.  In reality there is one.
    // Here we make sure that non-specification of MIX_DONORS works.
    r###"BCR="123085;123089" CDR3=CVKDRVTGTITELDYW"###,
    // 27. tests SUMMARY and NOPRINT
    r###"BCR=123085 SUMMARY SUMMARY_CLEAN NOPRINT"###,
    // 28. tests BARCODE option
    r###"BCR=165807 BARCODE=CCCATACGTGATGATA-1,TCTATTGAGCTGAAAT-1"###,
    // 29. tests parenthesized variable in F, SUM and MEAN
    r###"BCR=86237 GEX=85679 LVARSP=IGHV3-7_g_μ F="(IGHV3-7_g_μ)>=4.5" MIN_CHAINS=2 SUM MEAN"###,
    // 30. tests d_univ and d_donor
    r###"BCR=123085 CVARSP=d_univ,d_donor CDR3=CVKDRVTGTITELDYW"###,
    // 31. tests Cell Ranger 3.1 output
    r###"BCR=../3.1/123085 CDR3=CVKDRVTGTITELDYW"###,
    // 32. tests Cell Ranger 2.0 output and RE
    r###"BCR=../2.0/124550 CDR3=CAREPLYYDFWSAYFDYW RE"###,
    // 33. tests SCAN
    r###"BCR=123085 GEX=123201 LVARSP=IGHV1-69D_g_μ MIN_CELLS=10
     SCAN="(IGHV1-69D_g_μ)>=100,(IGHV1-69D_g_μ)<=1,t-10*c>=0.1" NOPRINT"###,
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
    // 39. tests u and r fields in parseable output
    r###"BCR=85333 POUT=stdout PCOLS=barcode,u_med1,u1,r_med2,r2 PCELL CVARSP=r_med
        CDR3=CAADGGGDQYYYMDVW"###,
];
