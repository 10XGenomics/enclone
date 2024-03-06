// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Temporary code for playing with melting temperatures.
//
// Barcodes are from:
// enclone PRE= 153782,153783,153784 PER_CELL SEG="IGHG1|IGHG2|IGHG3|IGHG4"
// Only used IGHG, excluded clonotypes having three chains.

use pretty_trace::*;
use string_utils::*;
use tables::*;

fn main() {
    PrettyTrace::new().on();

    // Define barcodes.  Data show (id = clonotype.cell, barcode, constant regions),
    // where the numbering scheme for clonotype and cell is for this list.

    let mut source = Vec::<(String, Vec<u8>, String)>::new();
    source.push((
        "1 = 1.1".to_string(),
        b"ATGTGTGAGACTGGGT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "2 = 1.2".to_string(),
        b"CACACCTAGATACACA".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "3 = 1.3".to_string(),
        b"GTGCAGCCATTTCACT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "4 = 1.4".to_string(),
        b"TGACGGCCAATGGAGC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "5 = 2.1".to_string(),
        b"GGGAATGGTAAGGGCT".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "6 = 3.1".to_string(),
        b"AGCAGCCAGATGTGGC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "7 = 3.2".to_string(),
        b"TGTTCCGGTTACTGAC".to_vec(),
        "IGHG1, IGKC".to_string(),
    ));
    source.push((
        "8 = 4.1".to_string(),
        b"AAAGTAGTCCTCCTAG".to_vec(),
        "IGHG3, IGKC".to_string(),
    ));
    source.push((
        "9 = 5.1".to_string(),
        b"ATCTACTCACGAAACG".to_vec(),
        "IGHG1, IGLC1".to_string(),
    ));
    source.push((
        "10 = 5.2".to_string(),
        b"TTAGGACCAGTACACT".to_vec(),
        "IGHG1, IGLC1".to_string(),
    ));
    source.push((
        "11 = 6.1".to_string(),
        b"TCGTACCGTTCCATGA".to_vec(),
        "IGHG1, IGLC3".to_string(),
    ));
    source.push((
        "12 = 6.2".to_string(),
        b"TGCTGCTGTAAATACG".to_vec(),
        "IGHG1, IGLC3".to_string(),
    ));

    // to delete:
    // enclone BCR=123085 PER_CELL AMINO=cdr3 MIN_CHAINS_EXACT=2 CDR3=CARDRIAGRFGYGMDVW
    //         POUT=stdout PCOLS=barcode PCELL

    // Compute Tm for the relevant constant region primers.

    println!("");
    let s = "TCCTGAGGACTGTAGGACAGC";
    let tg = tm_nearest_neighbor(&s);
    println!("IGHG = {} ==> {:.1}°", s, tg);
    let s = "TAGCTGCTGGCCGC";
    let tl = tm_nearest_neighbor(&s);
    println!("IGLC = {} ==> {:.1}°", s, tl);
    let s = "GCGTTATCCACCTTCCACTGT";
    let tk = tm_nearest_neighbor(&s);
    println!("IGKC = {} ==> {:.1}°", s, tk);
    println!("");
    let tx = 57.5;

    // Sequence of pR1.  This is immediately before the barcode.

    let pr1 = b"CTACACGACGCTCTTCCGATCT";
    let npr = pr1.len();

    // Compute Tm for ...

    /*
    let s = [
        "ATGCGATGTCTCAACA",
        "CAACTAGGTAGAGGAA",
        "CATCCACCAGCGAACA",
        "CCTTTCTAGGACGAAA",
        "CGATCGGCATTGGCGC",
        "CTCGAAACAAGCCGTC",
        "TAAGCGTTCAAAGACA",
        "TGGCTGGAGGATGTAT",
        "TGTATTCAGGTCGGAT",
        "GATCAGTGTCGAGTTT",
    ];
    */
    let mut rows = Vec::<Vec<String>>::new();
    let row = vec![
        "#".to_string(),
        "const".to_string(),
        "method1".to_string(),
        "t1".to_string(),
        "method2".to_string(),
        "t2".to_string(),
    ];
    rows.push(row);
    rows.push(vec!["\\hline".to_string(); 6]);
    for i in 0..source.len() {
        let b = source[i].1.clone(); // barcode

        // Method 1: (tail of pR1 + barcode).

        let mut res = Vec::<(f64, f64, String)>::new();
        for j in 0..=pr1.len() {
            let x = format!("{}{}", strme(&pr1[npr - j..npr]), strme(&b));
            let mut xp = format!("{}+{}", strme(&pr1[npr - j..npr]), strme(&b));
            if j == 0 {
                xp = x.clone();
            }
            let t = tm_nearest_neighbor(&x);
            res.push(((t - tx).abs(), t, xp));
        }
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let t1 = res[0].1;
        let xp1 = res[0].2.clone();

        // Method 2: (tail of pR1 + first 8 bases of barcode).

        let b = &b[0..8];
        let mut res = Vec::<(f64, f64, String)>::new();
        for j in 0..=pr1.len() {
            let x = format!("{}{}", strme(&pr1[npr - j..npr]), strme(&b));
            let mut xp = format!("{}+{}", strme(&pr1[npr - j..npr]), strme(&b));
            if j == 0 {
                xp = x.clone();
            }
            let t = tm_nearest_neighbor(&x);
            res.push(((t - tx).abs(), t, xp));
        }
        res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let t2 = res[0].1;
        let xp2 = res[0].2.clone();

        // Add row to table.

        let row = vec![
            source[i].0.clone(),
            source[i].2.clone(),
            xp1,
            format!("{:.1}°", t1),
            xp2,
            format!("{:.1}°", t2),
        ];
        rows.push(row);
    }
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, &b"l|l|l|l|l|l".to_vec(), false, false);
    println!("{}", log);
    for i in 0..rows.len() {
        if i == 1 {
            continue;
        }
        for j in 0..rows[i].len() {
            let mut s = rows[i][j].clone();
            if s.contains(',') {
                s = format!("\"{}\"", s);
            }
            if j > 0 {
                print!(",");
            }
            print!("{}", s);
        }
        println!("");
    }
}

// This file provides a function tm_nearest_neighbor.  All the code in this file is a verbatim
// translation to rust of C++ code in the BroadCRD codebase (copyright 2006), as
// dna/DNAHybridization.{cc,h}.  It is conceivable that mistakes were introduced in translation.
// That code in turn was a verbatim translation of calculations in published work as cited below.

// tm_nearest_neighbor: compute melting temperature of a given DNA sequence s
// at molarity s_mol (default 0.25 microM), where Na+ concentration is na_mol
// (default 50 mM), using a nearest-neighbor model, as used by
//
// http://www.idtdna.com/ANALYZER/Applications/OligoAnalyzer
//
// and described in
//
// [1] Allawi, H. T. and SantaLucia J. Jr.  Thermodynamics and NMR of internal
// G.T mismatches in DNA.  Biochemistry 36 (1997), 10581-10594.
//
// [2] Owczarzy, et al.  Effects of sodium ions on DNA duplex oligomers: improved
// predictions of melting temperatures.  Biochemistry 43 (2004), 3537-3554.
//
// [3] Patricia M. McTigue, Raymond J. Peterson, Jason D. Kuhn.  Sequence-dependent
// thermodynamic parameters for locked nucleic acid (LNA)-DNA duplex formation.
// Biochemistry 43 (2004), 5388-5405.
// (See also www.celadonlabs.com/Software/ModChem/default.aspx.)
//
// This method should be quite good for short sequences.
//
// The results are close to (but not identical to) the results given by the IDT
// OligoAnalyzer.  I don't know what the discrepancy is due to.
//
// We allow some nucleotides to be locked, as specified by the vector "locked".
// The following assumptions are enforced:
// (a) no two locked bases in a row;
// (b) no locked base at the beginning or end, or adjacent to those positions.
// These are consistent with the assumptions in [3].
//
// Alternatively, if S contains + symbols, they will be interpreted as instructions
// to "lock the following nucleotide", and the vector "locked" will be generated
// accordingly.

fn tm_nearest_neighbor(s: &str) -> f64 {
    let locked = Vec::<bool>::new();
    tm_nearest_neighbor_full(s, 0.00000025, 0.05, &locked)
}

fn tm_nearest_neighbor_full(s: &str, s_mol: f64, na_mol: f64, locked: &[bool]) -> f64 {
    // Allow for + symbols.

    if s.contains('+') {
        assert_eq!(locked.len(), 0);
        assert!(!s.contains("++"));
        assert!(!s.ends_with('+'));
        let mut sx = String::new();
        let schars: Vec<char> = s.chars().collect();
        let mut lockedx = Vec::<bool>::with_capacity(schars.len());
        let mut i = 0;
        while i < schars.len() {
            if schars[i] != '+' {
                sx.push(schars[i]);
                lockedx.push(false);
            } else {
                sx.push(schars[i + 1]);
                lockedx.push(true);
                i += 1;
            }
        }
        return tm_nearest_neighbor_full(&sx, s_mol, na_mol, &lockedx);
    }

    // Validate sequence.

    verify_dna(s);

    // Compute thermodynamic sums.

    let mut dh_sum = 0.0;
    let mut ds_sum = 0.0;
    let mut dg_sum = 0.0;
    thermodynamic_sums_dna(s, &mut dh_sum, &mut ds_sum, &mut dg_sum, true, true, locked);

    // Compute melting temperature based on nearest-neighbor model.

    let ideal_gas_const = 1.987; // calories per Kelvin per mole
    let kelvin_to_celsius = 273.15;
    let mut temp = 1000.0 * dh_sum / (ds_sum + ideal_gas_const * s_mol.ln()) - kelvin_to_celsius;

    // Correct for Na concentration, following [2].

    let mut sx = Vec::<char>::new();
    for c in s.chars() {
        sx.push(c);
    }
    let mut gc = 0;
    for i in 0..sx.len() {
        if sx[i] == 'G' || sx[i] == 'C' {
            gc += 1;
        }
    }
    let gc_fract = gc as f64 / sx.len() as f64;
    let ln_na_mol = na_mol.ln();

    temp = -kelvin_to_celsius
        + 1.0
            / (1.0 / (temp + kelvin_to_celsius)
                + (4.29 * gc_fract - 3.95) * 0.00001 * ln_na_mol
                + 9.40 * 0.000001 * ln_na_mol * ln_na_mol);
    temp
}

fn verify_dna(s: &str) {
    for c in s.chars() {
        assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T');
    }
}

// get_thermodynamic_parameters_dna.
//
// Return nearest neighbor thermodynamic parameters, from table 1 in reference [1],
// Allawi and SantaLucia (1997).
//
// dH = enthalpy (kcal/mole) = dH^\circ
// dS = entropy = dS^\circ
// dG = Gibb's free energy at 37 degrees = dG^\circ_37
// dH, dS, dG are all {A,C,G,T} x {A,C,G,T} matrices.
// (d means delta)
//
// Although there are dG entries in the table, I've computed them here directly
// from dH and dS.  I checked a few entries and they agree with what's in the table.
// dG = dH - T dS (theoretical)
// dG = dH - (37 + 273.15) dS / 1000 (used).
// I'm not sure quite why the factor of 1000 is needed.

fn get_thermodynamic_parameters_dna(
    dh: &mut Vec<Vec<f64>>,
    ds: &mut Vec<Vec<f64>>,
    dg: &mut Vec<Vec<f64>>,
    dh_g_or_c_init: &mut f64,
    dh_a_or_t_init: &mut f64,
    ds_g_or_c_init: &mut f64,
    ds_a_or_t_init: &mut f64,
    dg_g_or_c_init: &mut f64,
    dg_a_or_t_init: &mut f64,
    dh_symmetry_correction: &mut f64,
    ds_symmetry_correction: &mut f64,
    dg_symmetry_correction: &mut f64,
) {
    let a = 0;
    let c = 1;
    let g = 2;
    let t = 3;

    *dh = vec![vec![0.0; 4]; 4];
    *ds = vec![vec![0.0; 4]; 4];
    *dg = vec![vec![0.0; 4]; 4];

    dh[a][a] = -7.9;
    dh[t][t] = -7.9;
    dh[a][t] = -7.2;
    dh[t][a] = -7.2;
    dh[c][a] = -8.5;
    dh[t][g] = -8.5;
    dh[g][t] = -8.4;
    dh[a][c] = -8.4;
    dh[c][t] = -7.8;
    dh[a][g] = -7.8;
    dh[g][a] = -8.2;
    dh[t][c] = -8.2;
    dh[c][g] = -10.6;
    dh[g][c] = -9.8;
    dh[g][g] = -8.0;
    dh[c][c] = -8.0;
    ds[a][a] = -22.2;
    ds[t][t] = -22.2;
    ds[a][t] = -20.4;
    ds[t][a] = -21.3;
    ds[c][a] = -22.7;
    ds[t][g] = -22.7;
    ds[g][t] = -22.4;
    ds[a][c] = -22.4;
    ds[c][t] = -21.0;
    ds[a][g] = -21.0;
    ds[g][a] = -22.2;
    ds[t][c] = -22.2;
    ds[c][g] = -27.2;
    ds[g][c] = -24.4;
    ds[g][g] = -19.9;
    ds[c][c] = -19.9;

    /*
    // here are the parameters from santaLucia and hicks (2004).  they are almost the same.
    dh[a][a] = dh[t][t] = -7.6;
    dh[a][t] = -7.2;
    dh[t][a] = -7.2;
    dh[c][a] = dh[t][g] = -8.5;
    dh[g][t] = dh[a][c] = -8.4;
    dh[c][t] = dh[a][g] = -7.8;
    dh[g][a] = dh[t][c] = -8.2;
    dh[c][g] = -10.6;
    dh[g][c] = -9.8;
    dh[g][g] = dh[c][c] = -8.0;
    ds[a][a] = ds[t][t] = -21.3;
    ds[a][t] = -20.4;
    ds[t][a] = -21.3;
    ds[c][a] = ds[t][g] = -22.7;
    ds[g][t] = ds[a][c] = -22.4;
    ds[c][t] = ds[a][g] = -21.0;
    ds[g][a] = ds[t][c] = -22.2;
    ds[c][g] = -27.2;
    ds[g][c] = -24.4;
    ds[g][g] = ds[c][c] = -19.9;
    */

    let temp = 37.0 + 273.15;
    for i in 0..4 {
        for j in 0..4 {
            dg[i][j] = dh[i][j] - temp * ds[i][j] / 1000.0;
        }
    }
    *dh_g_or_c_init = 0.1;
    *dh_a_or_t_init = 2.3;
    *ds_g_or_c_init = -2.8;
    *ds_a_or_t_init = 4.1;
    *dg_g_or_c_init = *dh_g_or_c_init - temp * *ds_g_or_c_init / 1000.0;
    *dg_a_or_t_init = *dh_a_or_t_init - temp * *ds_a_or_t_init / 1000.0;
    *dh_symmetry_correction = 0.0;
    *ds_symmetry_correction = -1.4;
    *dg_symmetry_correction = *dh_symmetry_correction - temp * *ds_symmetry_correction / 1000.0;
}

// get_locked_thermodynamic_parameters_dna.  Return the 32 nearest-neighbor
// thermodynamic parameter set for LNA incorporation, from the left half of Table 4
// in reference [3].

fn get_locked_thermodynamic_parameters_dna(
    ddh_left: &mut Vec<Vec<f64>>,
    dds_left: &mut Vec<Vec<f64>>,
    ddg_left: &mut Vec<Vec<f64>>,
    ddh_right: &mut Vec<Vec<f64>>,
    dds_right: &mut Vec<Vec<f64>>,
    ddg_right: &mut Vec<Vec<f64>>,
) {
    let a = 0;
    let c = 1;
    let g = 2;
    let t = 3;
    *ddh_left = vec![vec![0.0; 4]; 4];
    *dds_left = vec![vec![0.0; 4]; 4];
    *ddg_left = vec![vec![0.0; 4]; 4];
    *ddh_right = vec![vec![0.0; 4]; 4];
    *dds_right = vec![vec![0.0; 4]; 4];
    *ddg_right = vec![vec![0.0; 4]; 4];
    ddh_left[a][a] = 0.707;
    ddh_left[a][c] = 1.131;
    ddh_left[a][g] = 0.264;
    ddh_left[a][t] = 2.282;
    ddh_left[c][a] = 1.049;
    ddh_left[c][c] = 2.096;
    ddh_left[c][g] = 0.785;
    ddh_left[c][t] = 0.708;
    ddh_left[g][a] = 3.162;
    ddh_left[g][c] = -0.360;
    ddh_left[g][g] = -2.844;
    ddh_left[g][t] = -0.212;
    ddh_left[t][a] = -0.046;
    ddh_left[t][c] = 1.893;
    ddh_left[t][g] = -1.540;
    ddh_left[t][t] = 1.528;
    dds_left[a][a] = 2.477;
    dds_left[a][c] = 4.064;
    dds_left[a][g] = 2.613;
    dds_left[a][t] = 7.457;
    dds_left[c][a] = 4.320;
    dds_left[c][c] = 7.996;
    dds_left[c][g] = 3.709;
    dds_left[c][t] = 4.175;
    dds_left[g][a] = 10.544;
    dds_left[g][c] = -0.251;
    dds_left[g][g] = -6.680;
    dds_left[g][t] = 0.073;
    dds_left[t][a] = 1.562;
    dds_left[t][c] = 6.685;
    dds_left[t][g] = -3.044;
    dds_left[t][t] = 5.298;
    ddg_left[a][a] = -0.092;
    ddg_left[a][c] = -0.122;
    ddg_left[a][g] = -0.561;
    ddg_left[a][t] = -0.007;
    ddg_left[c][a] = -0.270;
    ddg_left[c][c] = -0.457;
    ddg_left[c][g] = -0.332;
    ddg_left[c][t] = -0.666;
    ddg_left[g][a] = -0.072;
    ddg_left[g][c] = -0.414;
    ddg_left[g][g] = -0.700;
    ddg_left[g][t] = -0.194;
    ddg_left[t][a] = -0.563;
    ddg_left[t][c] = -0.208;
    ddg_left[t][g] = -0.548;
    ddg_left[t][t] = -0.130;
    ddh_right[a][a] = 0.992;
    ddh_right[a][c] = 2.890;
    ddh_right[a][g] = -1.200;
    ddh_right[a][t] = 1.816;
    ddh_right[c][a] = 1.358;
    ddh_right[c][c] = 2.063;
    ddh_right[c][g] = -0.276;
    ddh_right[c][t] = -1.671;
    ddh_right[g][a] = 0.444;
    ddh_right[g][c] = -0.925;
    ddh_right[g][g] = -0.943;
    ddh_right[g][t] = -0.635;
    ddh_right[t][a] = 1.591;
    ddh_right[t][c] = 0.609;
    ddh_right[t][g] = 2.165;
    ddh_right[t][t] = 2.326;
    dds_right[a][a] = 4.065;
    dds_right[a][c] = 10.576;
    dds_right[a][g] = -1.826;
    dds_right[a][t] = 6.863;
    dds_right[c][a] = 4.367;
    dds_right[c][c] = 7.565;
    dds_right[c][g] = -0.718;
    dds_right[c][t] = -4.070;
    dds_right[g][a] = 2.898;
    dds_right[g][c] = -1.111;
    dds_right[g][g] = -0.933;
    dds_right[g][t] = -0.342;
    dds_right[t][a] = 5.281;
    dds_right[t][c] = 3.169;
    dds_right[t][g] = 7.163;
    dds_right[t][t] = 8.051;
    ddg_right[a][a] = -0.396;
    ddg_right[a][c] = -0.390;
    ddg_right[a][g] = -0.603;
    ddg_right[a][t] = -0.309;
    ddg_right[c][a] = 0.046;
    ddg_right[c][c] = -0.404;
    ddg_right[c][g] = -0.003;
    ddg_right[c][t] = -0.409;
    ddg_right[g][a] = -0.437;
    ddg_right[g][c] = -0.535;
    ddg_right[g][g] = -0.666;
    ddg_right[g][t] = -0.520;
    ddg_right[t][a] = 0.004;
    ddg_right[t][c] = -0.396;
    ddg_right[t][g] = -0.106;
    ddg_right[t][t] = -0.212;
}

// thermodynamic_sums_dna.  Compute dh_sum, ds_sum, dg_sum.  We allow some nucleotides
// to be locked, as specified by the variable "locked".  The following assumptions
// are enforced:
// (a) no two locked bases in a row;
// (b) no locked base at the beginning or end, or adjacent to those positions.
// These are consistent with the assumptions in [3].

fn thermodynamic_sums_dna(
    s: &str,
    dh_sum: &mut f64,
    ds_sum: &mut f64,
    dg_sum: &mut f64,
    include_symmetry_correction: bool,
    include_initiation_terms: bool,
    locked: &[bool],
) {
    // defaults for last: true, true, empty

    let mut dh = Vec::<Vec<f64>>::new();
    let mut ds = Vec::<Vec<f64>>::new();
    let mut dg = Vec::<Vec<f64>>::new();

    let mut dh_g_or_c_init = 0.0;
    let mut dh_a_or_t_init = 0.0;
    let mut ds_g_or_c_init = 0.0;
    let mut ds_a_or_t_init = 0.0;
    let mut dg_g_or_c_init = 0.0;
    let mut dg_a_or_t_init = 0.0;
    let mut dh_symmetry_correction = 0.0;
    let mut ds_symmetry_correction = 0.0;
    let mut dg_symmetry_correction = 0.0;

    get_thermodynamic_parameters_dna(
        &mut dh,
        &mut ds,
        &mut dg,
        &mut dh_g_or_c_init,
        &mut dh_a_or_t_init,
        &mut ds_g_or_c_init,
        &mut ds_a_or_t_init,
        &mut dg_g_or_c_init,
        &mut dg_a_or_t_init,
        &mut dh_symmetry_correction,
        &mut ds_symmetry_correction,
        &mut dg_symmetry_correction,
    );
    *dh_sum = 0.0;
    *ds_sum = 0.0;
    *dg_sum = 0.0;
    if include_symmetry_correction {
        *dh_sum += dh_symmetry_correction;
        *ds_sum += ds_symmetry_correction;
        *dg_sum += dg_symmetry_correction;
    }
    let mut sx = Vec::<char>::new();
    for c in s.chars() {
        sx.push(c);
    }
    if include_initiation_terms {
        if sx[0] == 'A' || sx[0] == 'T' {
            *dh_sum += dh_a_or_t_init;
            *ds_sum += ds_a_or_t_init;
            *dg_sum += dg_a_or_t_init;
        } else {
            *dh_sum += dh_g_or_c_init;
            *ds_sum += ds_g_or_c_init;
            *dg_sum += dg_g_or_c_init;
        }
        if sx[sx.len() - 1] == 'A' || sx[sx.len() - 1] == 'T' {
            *dh_sum += dh_a_or_t_init;
            *ds_sum += ds_a_or_t_init;
            *dg_sum += dg_a_or_t_init;
        } else {
            *dh_sum += dh_g_or_c_init;
            *ds_sum += ds_g_or_c_init;
            *dg_sum += dg_g_or_c_init;
        }
    }
    let mut sx = Vec::<char>::new();
    for c in s.chars() {
        sx.push(c);
    }
    let mut b = Vec::<usize>::new();
    for i in 0..sx.len() {
        if sx[i] == 'A' {
            b.push(0);
        } else if sx[i] == 'C' {
            b.push(1);
        } else if sx[i] == 'G' {
            b.push(2);
        } else {
            b.push(3);
        }
    }
    for i in 0..sx.len() - 1 {
        *dh_sum += dh[b[i]][b[i + 1]];
        *ds_sum += ds[b[i]][b[i + 1]];
        *dg_sum += dg[b[i]][b[i + 1]];
    }

    // Handle locked bases.

    let mut have_lock = false;
    for i in 0..locked.len() {
        if locked[i] {
            have_lock = true;
        }
    }
    if have_lock {
        for i in 1..locked.len() {
            if locked[i] {
                assert!(!locked[i - 1]);
            }
        }
        assert!(locked.len() >= 5);
        assert!(!locked[0] && !locked[1]);
        assert!(!locked[locked.len() - 1]);
        assert!(!locked[locked.len() - 2]);
        let mut ddh_left = Vec::<Vec<f64>>::new();
        let mut dds_left = Vec::<Vec<f64>>::new();
        let mut ddg_left = Vec::<Vec<f64>>::new();
        let mut ddh_right = Vec::<Vec<f64>>::new();
        let mut dds_right = Vec::<Vec<f64>>::new();
        let mut ddg_right = Vec::<Vec<f64>>::new();
        get_locked_thermodynamic_parameters_dna(
            &mut ddh_left,
            &mut dds_left,
            &mut ddg_left,
            &mut ddh_right,
            &mut dds_right,
            &mut ddg_right,
        );
        for i in 0..locked.len() {
            if locked[i] {
                *dh_sum += ddh_left[b[i]][b[i + 1]] + ddh_right[b[i - 1]][b[i]];
                *ds_sum += dds_left[b[i]][b[i + 1]] + dds_right[b[i - 1]][b[i]];
                *dg_sum += ddg_left[b[i]][b[i + 1]] + ddg_right[b[i - 1]][b[i]];
            }
        }
    }
}
