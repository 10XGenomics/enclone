// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Tools to work with PDB files.

use amino::*;
use flate2::read::MultiGzDecoder;
use io_utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;

#[derive(Default, Clone)]
pub struct PdbStructure {
    pub chain_names: Vec<String>,
    pub chains: Vec<Vec<u8>>,
    pub atoms: Vec<PdbAtom>,
}

#[derive(Default, Clone)]
pub struct PdbAtom {
    pub atom_name: [u8; 2], // e.g. C for carbon etc.; second byte zero unless two letter atom
    pub amino_acid: u8,     // e.g. L for leucine etc.
    pub chain: u16,         // zero-based chain number
    pub chain_pos: u16,     // zero-based index of amino acid on chain
    pub x: f32,             // x coordinate of atom
    pub y: f32,             // y coordinate of atom
    pub z: f32,             // y coordinate of atom
}

impl PdbStructure {
    pub fn fetch_atoms_range(&self, chain: usize, start: usize, stop: usize) -> Vec<[f32; 3]> {
        let mut u = Vec::<[f32; 3]>::new();
        for j in 0..self.atoms.len() {
            if self.atoms[j].chain as usize == chain {
                if (self.atoms[j].chain_pos as usize) >= start
                    && (self.atoms[j].chain_pos as usize) < stop
                {
                    let x_s = self.atoms[j].x;
                    let y_s = self.atoms[j].y;
                    let z_s = self.atoms[j].z;
                    u.push([x_s, y_s, z_s]);
                }
            }
        }
        return u;
    }
}

pub fn fetch_pdb_structure(name: &str) -> Result<PdbStructure, ()> {
    fetch_pdb_structure_gen(&name, "antibody_sets/pdbs")
}

pub fn fetch_pdb_structure_gen(name: &str, dir: &str) -> Result<PdbStructure, ()> {
    let mut lines = Vec::<String>::new();
    let opath = format!("{}/{}.gz", dir, name);
    if path_exists(&opath) {
        let gz = MultiGzDecoder::new(File::open(&opath).unwrap());
        let b = BufReader::new(gz);
        for line in b.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
    } else {
        let x = Command::new("wget")
            .arg("-q")
            .arg("-O")
            .arg("-")
            .arg(&format!("https://files.rcsb.org/download/{}.cif", name))
            .output()
            .expect(&format!("failed to execute wget"));
        for line in strme(&x.stdout).lines() {
            lines.push(line.to_string());
        }
    }
    if lines.is_empty() {
        eprintln!("");
        eprintme!(name, dir);
        eprintln!("No lines read.\n");
        std::process::exit(1);
    }
    let mut chain_names = Vec::<String>::new();
    let mut chains = Vec::<Vec<u8>>::new();
    let mut ms = Vec::<PdbAtom>::new();
    let mut i = 0;
    while i < lines.len() {
        // Parse the chain names.  This is madness.

        if lines[i].starts_with("_entity.details") {
            let mut elines = Vec::<String>::new();
            let mut j = i + 1;
            while !lines[j].starts_with('#') {
                elines.push(lines[j].to_string());
                j += 1;
            }
            let mtypes = ["polymer", "non-polymer", "branched", "water", "?"];
            let mut count = 1;
            let mut starts = Vec::<usize>::new();
            for m in 0..elines.len() {
                let mut s = elines[m].clone();
                while s.contains("  ") {
                    s = s.replace("  ", " ");
                }
                if s.starts_with(&format!("{} ", count)) {
                    if s.after(" ").contains(" ") {
                        let mtype = s.between(" ", " ");
                        let mut known = false;
                        for u in 0..mtypes.len() {
                            if mtype == mtypes[u] {
                                known = true;
                            }
                        }
                        if !known {
                            eprintln!("unknown type {}", mtype);
                            std::process::exit(1);
                        }
                    }
                    starts.push(m);
                    count += 1;
                }
            }
            starts.push(elines.len());
            let mut flines = Vec::<String>::new();
            for z in 0..starts.len() - 1 {
                let mut s = String::new();
                for k in starts[z]..starts[z + 1] {
                    s += &elines[k];
                }
                flines.push(s);
            }
            for m in 0..flines.len() {
                let mut s = flines[m].replace("  ", " ");
                while s.contains("  ") {
                    s = s.replace("  ", " ");
                }
                if s.between(" ", " ") == "polymer" {
                    let mut t = s.after("polymer").to_string();
                    if t.contains("'") && t.after("'").contains("'") {
                        t = t.between("'", "'").to_string();
                    }
                    chain_names.push(t);
                }
            }
            i = j + 1;

        // Parse the chain sequences.
        } else if lines[i].starts_with("_entity_poly.pdbx_target_identifier") {
            let mut j = i + 1;
            while j < lines.len() && !lines[j].starts_with('#') {
                if lines[j].contains("'polypeptide(L)'") {
                    j += 1;
                    if lines[j].starts_with("#") {
                        break;
                    }
                    if !lines[j].contains(";") {
                        return Err(());
                    }
                    let mut s = lines[j].after(";").as_bytes().to_vec();
                    while j < lines.len() {
                        j += 1;
                        if lines[j].starts_with(";") {
                            break;
                        }
                        s.append(&mut lines[j].as_bytes().to_vec());
                    }
                    chains.push(s);
                } else {
                    j += 1;
                }
            }
            i = j;
        } else if lines[i].starts_with("ATOM ") {
            let fields = lines[i].split_ascii_whitespace().collect::<Vec<&str>>();
            let atom0 = fields[2].as_bytes();
            let mut atom = [0 as u8; 2];
            atom[0] = atom0[0];
            if atom0.len() == 2 {
                atom[1] = atom0[1];
            }
            assert!(atom0.len() <= 2);

            // Check for defective amino acid.  Possibly these cases could be rescued.

            if fields[5].len() != 3 || fields[5] == "UNK" {
                return Err(());
            }

            // Create atom.

            let aa = aa3_to_aa(&fields[5].as_bytes());
            let pos_on_entity = fields[8].force_usize() as u16;
            let entity = fields[7].force_usize() as u16;
            let xcoord = fields[10].force_f64() as f32;
            let ycoord = fields[11].force_f64() as f32;
            let zcoord = fields[12].force_f64() as f32;
            let m = PdbAtom {
                atom_name: atom,
                amino_acid: aa,
                chain: entity - 1,
                chain_pos: pos_on_entity - 1,
                x: xcoord,
                y: ycoord,
                z: zcoord,
            };
            ms.push(m);
            i += 1;
        } else {
            i += 1;
        }
    }
    if chains.len() != chain_names.len() {
        Err(())
    } else {
        for i in 0..chains.len() {
            loop {
                let mut found = false;
                for j in 0..chains[i].len() - 4 {
                    if chains[i][j] == b'(' && chains[i][j + 4] == b')' {
                        let s = stringme(&chains[i]);
                        let t = s.replace(&strme(&chains[i][j..j + 4]), "X");
                        chains[i] = t.as_bytes().to_vec();
                        found = true;
                        break;
                    }
                }
                if !found {
                    break;
                }
            }
            if chains[i].contains(&b'(') {
                eprintln!("\nProblem with {} in {}.\n", strme(&chains[i]), name);
                std::process::exit(1);
            }
        }
        Ok(PdbStructure {
            chain_names: chain_names,
            chains: chains,
            atoms: ms,
        })
    }
}
