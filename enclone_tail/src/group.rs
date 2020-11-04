// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.
//
// To keep compilation time down, this crate should not reach into the enclone crate.

use crate::display_tree::*;
use crate::neighbor::*;
use crate::newick::*;
use amino::*;
use ansi_escape::ansi_to_html::*;
use ansi_escape::*;
use enclone_core::defs::*;
use enclone_core::mammalian_fixed_len::*;
use enclone_core::print_tools::*;
use enclone_proto::types::*;
use equiv::EquivRel;
use io_utils::*;
use itertools::*;
use perf_stats::*;
use stats_utils::*;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;
use std::io::*;
use std::mem::swap;
use std::path::Path;
use std::time::{Instant, SystemTime, UNIX_EPOCH};
use string_utils::*;
use tables::*;
use tar::{Builder, Header};
use vdj_ann::refx::*;
use vector_utils::*;

pub fn group_and_print_clonotypes(
    tall: &Instant,
    refdata: &RefData,
    pics: &Vec<String>,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    join_info: &Vec<(usize, usize, bool, Vec<u8>)>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &Vec<HashMap<String, String>>,
    dref: &Vec<DonorReferenceItem>,
) {
    // Build index to join info.

    let mut to_join_info = vec![Vec::<usize>::new(); exact_clonotypes.len()];
    for i in 0..join_info.len() {
        to_join_info[join_info[i].0].push(i);
        to_join_info[join_info[i].1].push(i);
    }

    // Set up for parseable output.

    let mut parseable_fields = Vec::<String>::new();
    set_speakers(&ctl, &mut parseable_fields);
    #[allow(bare_trait_objects)]
    let mut pout = match ctl.parseable_opt.pout.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        "stdouth" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.parseable_opt.pout);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    let mut pcols = ctl.parseable_opt.pcols.clone();
    for i in 0..pcols.len() {
        pcols[i] = pcols[i].replace("_Σ", "_sum");
        pcols[i] = pcols[i].replace("_μ", "_mean");
    }
    if pcols.is_empty() {
        pcols = parseable_fields.clone();
    }
    let mut pcols2 = Vec::<String>::new();
    for i in 0..pcols.len() {
        if pcols[i].contains(":") {
            pcols2.push(pcols[i].before(":").to_string());
        } else {
            pcols2.push(pcols[i].clone());
        }
    }
    pcols = pcols2;
    if !ctl.parseable_opt.pout.is_empty()
        && ctl.parseable_opt.pout != "stdout".to_string()
        && ctl.parseable_opt.pout != "stdouth".to_string()
    {
        fwriteln!(pout, "{}", pcols.iter().format(","));
    }

    // Set up for fasta output.

    #[allow(bare_trait_objects)]
    let mut fout = match ctl.gen_opt.fasta_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.fasta_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    #[allow(bare_trait_objects)]
    let mut faaout = match ctl.gen_opt.fasta_aa_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.fasta_aa_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };

    // Set up for clustal output.

    let (mut clustal_aa, mut clustal_dna) = (None, None);
    if ctl.gen_opt.clustal_aa.len() > 0 && ctl.gen_opt.clustal_aa != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.clustal_aa).unwrap();
        clustal_aa = Some(Builder::new(file));
    }
    if ctl.gen_opt.clustal_dna.len() > 0 && ctl.gen_opt.clustal_dna != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.clustal_dna).unwrap();
        clustal_dna = Some(Builder::new(file));
    }

    // Set up for phylip output.

    let (mut phylip_aa, mut phylip_dna) = (None, None);
    if ctl.gen_opt.phylip_aa.len() > 0 && ctl.gen_opt.phylip_aa != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.phylip_aa).unwrap();
        phylip_aa = Some(Builder::new(file));
    }
    if ctl.gen_opt.phylip_dna.len() > 0 && ctl.gen_opt.phylip_dna != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.phylip_dna).unwrap();
        phylip_dna = Some(Builder::new(file));
    }

    // Set up for peer group output.

    #[allow(bare_trait_objects)]
    let mut pgout = match ctl.gen_opt.peer_group_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.peer_group_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    if ctl.gen_opt.peer_group_filename.len() > 0 {
        if ctl.gen_opt.peer_group_filename != "stdout".to_string() {
            if !ctl.gen_opt.peer_group_readable {
                fwriteln!(pgout, "group,clonotype,chain,pos,amino_acid,count");
            } else {
                fwriteln!(pgout, "group,clonotype,chain,pos,distribution");
            }
        }
    }

    // Group clonotypes and make output.

    let mut last_width = 0;
    let mut e: EquivRel = EquivRel::new(pics.len() as i32);
    if ctl.clono_group_opt.heavy_cdr3_aa {
        let mut all = Vec::<(String, usize)>::new();
        for i in 0..pics.len() {
            for x in exacts[i].iter() {
                for m in 0..exact_clonotypes[*x].share.len() {
                    let y = &exact_clonotypes[*x].share[m];
                    if y.left {
                        all.push((y.cdr3_aa.clone(), i));
                    }
                }
            }
        }
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                e.join(all[i].1 as i32, all[k].1 as i32);
            }
            i = j;
        }
    }
    if ctl.clono_group_opt.vj_refname || ctl.clono_group_opt.vj_refname_strong {
        let mut all = Vec::<(Vec<String>, usize)>::new();
        for i in 0..pics.len() {
            let ex = &exact_clonotypes[exacts[i][0]];
            let mut s = Vec::<String>::new();
            for j in 0..ex.share.len() {
                s.push(refdata.name[ex.share[j].v_ref_id].clone());
                s.push(refdata.name[ex.share[j].j_ref_id].clone());
            }
            s.sort();
            all.push((s, i));
        }
        // Note duplication with above code.
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                let m1 = all[i].1;
                let m2 = all[k].1;
                if ctl.clono_group_opt.vj_refname_strong {
                    let ex1 = &exact_clonotypes[exacts[m1][0]];
                    let ex2 = &exact_clonotypes[exacts[m2][0]];
                    let mut lens1 = Vec::<(usize, usize)>::new();
                    let mut lens2 = Vec::<(usize, usize)>::new();
                    for j in 0..ex1.share.len() {
                        lens1.push((ex1.share[j].seq_del.len(), ex1.share[j].cdr3_aa.len()));
                    }
                    for j in 0..ex2.share.len() {
                        lens2.push((ex2.share[j].seq_del.len(), ex2.share[j].cdr3_aa.len()));
                    }
                    lens1.sort();
                    lens2.sort();
                    if lens1 != lens2 {
                        continue;
                    }
                }
                e.join(m1 as i32, m2 as i32);
            }
            i = j;
        }
    }
    let mut groups = 0;
    let mut greps = Vec::<i32>::new();
    e.orbit_reps(&mut greps);

    // Sort so that larger groups (as measured by cells) come first.

    let mut grepsn = Vec::<(usize, usize)>::new();
    for i in 0..greps.len() {
        let mut o = Vec::<i32>::new();
        e.orbit(greps[i], &mut o);
        if o.len() < ctl.clono_group_opt.min_group {
            continue;
        }
        let mut n = 0;
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                n += exact_clonotypes[s[k]].clones.len();
            }
        }
        grepsn.push((n, i));
    }
    reverse_sort(&mut grepsn);

    // Echo command.

    let mut logx = Vec::<u8>::new();
    if ctl.gen_opt.echo {
        let args: Vec<String> = env::args().collect();
        fwriteln!(logx, "\n{}", args.iter().format(" "));
        if ctl.gen_opt.html {
            fwriteln!(logx, "");
        }
    }

    // Now print clonotypes.

    for z in 0..grepsn.len() {
        let i = grepsn[z].1;
        let n = grepsn[z].0;
        let mut o = Vec::<i32>::new();
        e.orbit(greps[i], &mut o);
        groups += 1;

        // Generate human readable output.  Getting the newlines right is tricky, so
        // they're marked.

        if !ctl.gen_opt.noprint {
            if !ctl.gen_opt.html && !ctl.gen_opt.ngroup {
                fwriteln!(logx, ""); // NEWLINE 1
            }

            // If we just printed a clonotype box, output a bar.

            if last_width > 0 {
                if ctl.gen_opt.ngroup || ctl.gen_opt.html {
                    fwriteln!(logx, ""); // NEWLINE 2
                }
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_eight_bit_color_escape(&mut log, 44);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwrite!(logx, "╺{}╸", "━".repeat(last_width - 2));
                if !ctl.gen_opt.ngroup {
                    fwriteln!(logx, ""); // NEWLINE 3
                }
                fwriteln!(logx, ""); // NEWLINE 4
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    fwrite!(logx, "{}", strme(&log));
                }
            }

            // If NGROUP is not on, output a GROUP line, including a newline at the end.

            if !ctl.gen_opt.ngroup {
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_bold_escape(&mut log);
                    emit_eight_bit_color_escape(&mut log, 27);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwrite!(
                    logx,
                    "[{}] GROUP = {} CLONOTYPES = {} CELLS",
                    groups,
                    o.len(),
                    n
                );
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwriteln!(logx, ""); // NEWLINE 5
            }
        }
        let mut group_ncells = 0;
        for j in 0..o.len() {
            let oo = o[j] as usize;
            for l in 0..exacts[oo].len() {
                group_ncells += exact_clonotypes[exacts[oo][l]].ncells();
            }
        }
        for j in 0..o.len() {
            let oo = o[j] as usize;
            if !ctl.gen_opt.noprint {
                if z > 0 || j > 0 || !(ctl.gen_opt.html && ctl.gen_opt.ngroup) {
                    fwrite!(logx, "\n"); // NEWLINE 6
                }
                if ctl.gen_opt.svg {
                    const FONT_SIZE: usize = 15;
                    let s = format!("[{}.{}] {}", groups, j + 1, pics[oo]);

                    // Generate svg.  This does not generate the shortest possible string.  One
                    // thing that could be done is to use only one text tag and instead use
                    // relative positions in the tspan tags to avoid repeating the font family,
                    // etc.  But there are probably other economizations.
                    //
                    // The other thing is that the aspect ratio is just a little bit off.

                    fwrite!(
                        logx,
                        "{}",
                        convert_text_with_ansi_escapes_to_svg(&s, "Menlo", FONT_SIZE)
                    );
                } else {
                    fwrite!(logx, "[{}.{}] {}", groups, j + 1, pics[oo]);
                }
            }
            let x = &pics[oo];
            let mut y = Vec::<char>::new();
            for c in x.chars() {
                y.push(c);
            }
            y.reverse();
            let mut m = 2;
            while m < y.len() {
                if y[m] == '\n' {
                    break;
                }
                m += 1;
            }
            last_width = m - 1;

            // Print join info.

            let mut ji = Vec::<usize>::new();
            for u in exacts[oo].iter() {
                ji.append(&mut to_join_info[*u].clone());
            }
            unique_sort(&mut ji);
            for i in 0..ji.len() {
                fwriteln!(logx, "{}", strme(&join_info[ji[i]].3));
            }

            // Generate clustal output.  See:
            // 1. http://meme-suite.org/doc/clustalw-format.html
            // 2. https://www.ebi.ac.uk/seqdb/confluence/display/THD/Help+-+Clustal+Omega+FAQ
            //    at "What do the consensus symbols mean in the alignment?".

            if ctl.gen_opt.clustal_aa.len() > 0 {
                let stdout = ctl.gen_opt.clustal_aa == "stdout".to_string();
                let mut data = Vec::<u8>::new();
                if stdout {
                    fwriteln!(logx, "");
                    fwriteln!(logx, "CLUSTALW\n");
                } else {
                    fwriteln!(data, "CLUSTALW\n");
                }
                let mut aa = Vec::<Vec<u8>>::new();
                let mut names = Vec::<String>::new();
                for (k, u) in exacts[oo].iter().enumerate() {
                    let ex = &exact_clonotypes[*u];
                    let mut seq = Vec::<u8>::new();
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_none() {
                            seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m] / 3]);
                        } else {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let mut z = ex.share[r].aa_mod_indel.clone();
                            seq.append(&mut z);
                        }
                    }
                    names.push(format!("{}.{}.{}", groups, j + 1, k + 1,));
                    aa.append(&mut vec![seq; 1]);
                }
                const W: usize = 60;
                const PAD: usize = 4;
                let mut name_width = 0;
                for i in 0..names.len() {
                    name_width = std::cmp::max(name_width, names[i].len());
                }
                for start in (0..aa[0].len()).step_by(W) {
                    if start > 0 {
                        if stdout {
                            fwriteln!(logx, "");
                        } else {
                            fwriteln!(data, "");
                        }
                    }
                    let stop = std::cmp::min(start + W, aa[0].len());
                    for i in 0..aa.len() {
                        if stdout {
                            fwrite!(logx, "{}", names[i]);
                            fwrite!(
                                logx,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(logx, "{}  {}", strme(&aa[i][start..stop]), stop - start);
                        } else {
                            fwrite!(data, "{}", names[i]);
                            fwrite!(
                                data,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(data, "{}  {}", strme(&aa[i][start..stop]), stop - start);
                        }
                    }
                    if stdout {
                        fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
                    } else {
                        fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
                    }
                    for p in start..stop {
                        let mut res = Vec::<u8>::new();
                        for i in 0..aa.len() {
                            res.push(aa[i][p]);
                        }
                        unique_sort(&mut res);
                        if res.solo() {
                            if stdout {
                                fwrite!(logx, "*");
                            } else {
                                fwrite!(data, "*");
                            }
                        } else {
                            let mut con = false;
                            'pass: for pass in 1..=2 {
                                let x: Vec<&[u8]>;
                                // Conservative mutations
                                if pass == 1 {
                                    x = vec![
                                        b"STA", b"NEQK", b"NHQK", b"NDEQ", b"QHRK", b"MILV",
                                        b"MILF", b"HY", b"FYW",
                                    ];
                                } else {
                                    // Semi-conservative mutations
                                    x = vec![
                                        b"CSA", b"ATV", b"SAG", b"STNK", b"STPA", b"SGND",
                                        b"SNDEQK", b"NDEQHK", b"NEQHRK", b"FVLIM", b"HFY",
                                    ];
                                }
                                for y in x.iter() {
                                    let mut sub = true;
                                    for c in res.iter() {
                                        if !y.contains(c) {
                                            sub = false;
                                            break;
                                        }
                                    }
                                    if sub {
                                        let sym;
                                        if pass == 1 {
                                            sym = ":";
                                        } else {
                                            sym = ".";
                                        }
                                        if stdout {
                                            fwrite!(logx, "{}", sym);
                                        } else {
                                            fwrite!(data, "{}", sym);
                                        }
                                        con = true;
                                        break 'pass;
                                    }
                                }
                            }
                            if !con {
                                if stdout {
                                    fwrite!(logx, " ");
                                } else {
                                    fwrite!(data, " ");
                                }
                            }
                        }
                    }
                    if stdout {
                        fwriteln!(logx, "");
                    } else {
                        fwriteln!(data, "");
                    }
                }
                if !stdout {
                    let mut header = Header::new_gnu();
                    header.set_size(data.len() as u64);
                    header.set_cksum();
                    header.set_mode(0o0644);
                    let now = SystemTime::now();
                    header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
                    let filename = format!("{}.{}", groups, j + 1);
                    clustal_aa
                        .as_mut()
                        .unwrap()
                        .append_data(&mut header, &filename, &data[..])
                        .unwrap();
                }
            }
            if ctl.gen_opt.clustal_dna.len() > 0 {
                let stdout = ctl.gen_opt.clustal_dna == "stdout".to_string();
                let mut data = Vec::<u8>::new();
                if stdout {
                    fwriteln!(logx, "");
                    fwriteln!(logx, "CLUSTALW\n");
                } else {
                    fwriteln!(data, "CLUSTALW\n");
                }
                let mut dna = Vec::<Vec<u8>>::new();
                let mut names = Vec::<String>::new();
                for (k, u) in exacts[oo].iter().enumerate() {
                    let ex = &exact_clonotypes[*u];
                    let mut seq = Vec::<u8>::new();
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_none() {
                            seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m]]);
                        } else {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let mut s = ex.share[r].seq_del_amino.clone();
                            seq.append(&mut s);
                        }
                    }
                    names.push(format!("{}.{}.{}", groups, j + 1, k + 1,));
                    dna.append(&mut vec![seq; 1]);
                }
                const W: usize = 60;
                const PAD: usize = 4;
                let mut name_width = 0;
                for i in 0..names.len() {
                    name_width = std::cmp::max(name_width, names[i].len());
                }
                for start in (0..dna[0].len()).step_by(W) {
                    if start > 0 {
                        if stdout {
                            fwriteln!(logx, "");
                        } else {
                            fwriteln!(data, "");
                        }
                    }
                    let stop = std::cmp::min(start + W, dna[0].len());
                    for i in 0..dna.len() {
                        if stdout {
                            fwrite!(logx, "{}", names[i]);
                            fwrite!(
                                logx,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(logx, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                        } else {
                            fwrite!(data, "{}", names[i]);
                            fwrite!(
                                data,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(data, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                        }
                    }
                    if stdout {
                        fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
                    } else {
                        fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
                    }
                    for p in start..stop {
                        let mut res = Vec::<u8>::new();
                        for i in 0..dna.len() {
                            res.push(dna[i][p]);
                        }
                        unique_sort(&mut res);
                        if res.solo() {
                            if stdout {
                                fwrite!(logx, "*");
                            } else {
                                fwrite!(data, "*");
                            }
                        } else {
                            if stdout {
                                fwrite!(logx, " ");
                            } else {
                                fwrite!(data, " ");
                            }
                        }
                    }
                    if stdout {
                        fwriteln!(logx, "");
                    } else {
                        fwriteln!(data, "");
                    }
                }
                if !stdout {
                    let mut header = Header::new_gnu();
                    header.set_size(data.len() as u64);
                    header.set_cksum();
                    header.set_mode(0o0644);
                    let now = SystemTime::now();
                    header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
                    let filename = format!("{}.{}", groups, j + 1);
                    clustal_dna
                        .as_mut()
                        .unwrap()
                        .append_data(&mut header, &filename, &data[..])
                        .unwrap();
                }
            }

            // Generate sequential PHYLIP output.  See:
            // 1. http://www.atgc-montpellier.fr/phyml/usersguide.php?type=phylip
            // 2. http://evolution.genetics.washington.edu/phylip/doc/sequence.html.
            // We don't fold lines because it may not be necessary.  See giant value for W;
            // code left in place in case it turns out that folding is needed.  This will involve
            // a bit more than lowering W.

            if ctl.gen_opt.phylip_aa.len() > 0 {
                let stdout = ctl.gen_opt.phylip_aa == "stdout".to_string();
                let mut data = Vec::<u8>::new();
                let mut nbases = 0;
                for m in 0..rsi[oo].mat.len() {
                    nbases += rsi[oo].seq_del_lens[m];
                }
                if stdout {
                    fwriteln!(logx, "");
                    fwriteln!(logx, "{} {}", exacts[oo].len(), nbases / 3);
                } else {
                    fwriteln!(data, "{} {}", exacts[oo].len(), nbases / 3);
                }
                let mut aa = Vec::<Vec<u8>>::new();
                let mut names = Vec::<String>::new();
                for (k, u) in exacts[oo].iter().enumerate() {
                    let ex = &exact_clonotypes[*u];
                    let mut seq = Vec::<u8>::new();
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_none() {
                            seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m] / 3]);
                        } else {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let mut z = ex.share[r].aa_mod_indel.clone();
                            seq.append(&mut z);
                        }
                    }
                    names.push(format!("{}", k + 1,));
                    aa.append(&mut vec![seq; 1]);
                }
                const W: usize = 10000;
                const PAD: usize = 4;
                let mut name_width = 0;
                for i in 0..names.len() {
                    name_width = std::cmp::max(name_width, names[i].len());
                }
                for start in (0..aa[0].len()).step_by(W) {
                    if start > 0 {
                        if stdout {
                            fwriteln!(logx, "");
                        } else {
                            fwriteln!(data, "");
                        }
                    }
                    let stop = std::cmp::min(start + W, aa[0].len());
                    for i in 0..aa.len() {
                        if stdout {
                            fwrite!(logx, "{}", names[i]);
                            fwrite!(
                                logx,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(logx, "{}", strme(&aa[i][start..stop]));
                        } else {
                            fwrite!(data, "{}", names[i]);
                            fwrite!(
                                data,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(data, "{}", strme(&aa[i][start..stop]));
                        }
                    }
                    if stdout {
                        fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
                    } else {
                        fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
                    }
                }
                if !stdout {
                    let mut header = Header::new_gnu();
                    header.set_size(data.len() as u64);
                    header.set_cksum();
                    header.set_mode(0o0644);
                    let now = SystemTime::now();
                    header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
                    let filename = format!("{}.{}", groups, j + 1);
                    phylip_aa
                        .as_mut()
                        .unwrap()
                        .append_data(&mut header, &filename, &data[..])
                        .unwrap();
                }
            }
            if ctl.gen_opt.phylip_dna.len() > 0 {
                let stdout = ctl.gen_opt.phylip_dna == "stdout".to_string();
                let mut data = Vec::<u8>::new();
                let mut nbases = 0;
                for m in 0..rsi[oo].mat.len() {
                    nbases += rsi[oo].seq_del_lens[m];
                }
                if stdout {
                    fwriteln!(logx, "");
                    fwriteln!(logx, "{} {}", exacts[oo].len(), nbases);
                } else {
                    fwriteln!(data, "{} {}", exacts[oo].len(), nbases);
                }
                let mut dna = Vec::<Vec<u8>>::new();
                let mut names = Vec::<String>::new();
                for (k, u) in exacts[oo].iter().enumerate() {
                    let ex = &exact_clonotypes[*u];
                    let mut seq = Vec::<u8>::new();
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_none() {
                            seq.append(&mut vec![b'-'; rsi[oo].seq_del_lens[m]]);
                        } else {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let mut s = ex.share[r].seq_del_amino.clone();
                            seq.append(&mut s);
                        }
                    }
                    names.push(format!("{}", k + 1,));
                    dna.append(&mut vec![seq; 1]);
                }
                const W: usize = 10000;
                const PAD: usize = 4;
                let mut name_width = 0;
                for i in 0..names.len() {
                    name_width = std::cmp::max(name_width, names[i].len());
                }
                for start in (0..dna[0].len()).step_by(W) {
                    if start > 0 {
                        if stdout {
                            fwriteln!(logx, "");
                        } else {
                            fwriteln!(data, "");
                        }
                    }
                    let stop = std::cmp::min(start + W, dna[0].len());
                    for i in 0..dna.len() {
                        if stdout {
                            fwrite!(logx, "{}", names[i]);
                            fwrite!(
                                logx,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(logx, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                        } else {
                            fwrite!(data, "{}", names[i]);
                            fwrite!(
                                data,
                                "{}",
                                strme(&vec![b' '; name_width + PAD - names[i].len()])
                            );
                            fwriteln!(data, "{}  {}", strme(&dna[i][start..stop]), stop - start);
                        }
                    }
                    if stdout {
                        fwrite!(logx, "{}", strme(&vec![b' '; name_width + PAD]));
                    } else {
                        fwrite!(data, "{}", strme(&vec![b' '; name_width + PAD]));
                    }
                }
                if !stdout {
                    let mut header = Header::new_gnu();
                    header.set_size(data.len() as u64);
                    header.set_cksum();
                    header.set_mode(0o0644);
                    let now = SystemTime::now();
                    header.set_mtime(now.duration_since(UNIX_EPOCH).unwrap().as_secs());
                    let filename = format!("{}.{}", groups, j + 1);
                    phylip_dna
                        .as_mut()
                        .unwrap()
                        .append_data(&mut header, &filename, &data[..])
                        .unwrap();
                }
            }

            // Generate experimental tree output (options NEWICK0 and TREE).

            if ctl.gen_opt.newick || ctl.gen_opt.tree != "".to_string() {
                // Compute the n x n distance matrix for the exact subclonotypes.

                let n = exacts[oo].len();
                let cols = rsi[oo].mat.len();
                let mut dist = vec![vec![0; n]; n];
                for i1 in 0..n {
                    for i2 in 0..n {
                        let ex1 = &exact_clonotypes[exacts[oo][i1]];
                        let ex2 = &exact_clonotypes[exacts[oo][i2]];
                        for m in 0..cols {
                            if rsi[oo].mat[m][i1].is_some() && rsi[oo].mat[m][i2].is_some() {
                                let r1 = rsi[oo].mat[m][i1].unwrap();
                                let r2 = rsi[oo].mat[m][i2].unwrap();
                                let seq1 = &ex1.share[r1].seq_del_amino;
                                let seq2 = &ex2.share[r2].seq_del_amino;
                                for j in 0..seq1.len() {
                                    if seq1[j] != seq2[j] {
                                        dist[i1][i2] += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // Add a zeroeth entry for a "root subclonotype" which is defined to have the
                // donor reference away from the recombination region, and is undefined within it.
                // Define its distance to actual exact subclonotypes by only computing away from
                // the recombination region.  This yields an (n+1) x (n+1) matrix.

                let mut droot = vec![0; n];
                for i in 0..n {
                    let ex = &exact_clonotypes[exacts[oo][i]];
                    for m in 0..cols {
                        if rsi[oo].mat[m][i].is_some() {
                            let r = rsi[oo].mat[m][i].unwrap();
                            let seq = &ex.share[r].seq_del_amino;
                            let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                            if rsi[oo].vpids[m].is_some() {
                                vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                            }
                            let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                            let z = seq.len();
                            for p in 0..z {
                                let b = seq[p];
                                if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                                    droot[i] += 1;
                                }
                                if p >= z - (jref.len() - ctl.heur.ref_j_trim)
                                    && b != jref[jref.len() - (z - p)]
                                {
                                    droot[i] += 1;
                                }
                            }
                        }
                    }
                }
                let mut distp = vec![vec![0.0; n + 1]; n + 1];
                for i1 in 0..n {
                    for i2 in 0..n {
                        distp[i1 + 1][i2 + 1] = dist[i1][i2] as f64;
                    }
                }
                for i in 0..n {
                    distp[i + 1][0] = droot[i] as f64;
                    distp[0][i + 1] = droot[i] as f64;
                }

                // Generate the neighborhood joining tree associated to these data.

                let mut tree = neighbor_joining(&distp);
                let mut nvert = 0;
                for i in 0..tree.len() {
                    nvert = max(nvert, tree[i].0 + 1);
                    nvert = max(nvert, tree[i].1 + 1);
                }

                // Use the root to direct the edges.

                let r = 0;
                let mut index = vec![Vec::<usize>::new(); nvert];
                for i in 0..tree.len() {
                    index[tree[i].0].push(i);
                    index[tree[i].1].push(i);
                }
                let mut rooted = vec![false; nvert];
                rooted[r] = true;
                let mut roots = vec![r];
                for i in 0..nvert {
                    let v = roots[i];
                    for j in index[v].iter() {
                        let e = &mut tree[*j];

                        if e.1 == v && !rooted[e.0] {
                            swap(&mut e.0, &mut e.1);
                        }
                        if e.0 == v && !rooted[e.1] {
                            rooted[e.1] = true;
                            roots.push(e.1);
                        }
                    }
                }

                // Output in Newick format.

                if ctl.gen_opt.newick {
                    let mut vnames = Vec::<String>::new();
                    for i in 0..=n {
                        vnames.push(format!("{}", i));
                    }
                    let mut edges = Vec::<(usize, usize, String)>::new();
                    for i in 0..tree.len() {
                        edges.push((tree[i].0, tree[i].1, format!("{:.2}", tree[i].2)));
                    }
                    for i in n + 1..nvert {
                        vnames.push(format!("I{}", i - n));
                    }
                    let nw = newick(&vnames, 0, &edges);
                    fwriteln!(logx, "\n{}", nw);
                }

                // Output as visual tree.

                if ctl.gen_opt.tree != "".to_string() {
                    let mut edges = Vec::<(usize, usize, f64)>::new();
                    let mut nvert = 0;
                    for i in 0..tree.len() {
                        edges.push((tree[i].0, tree[i].1, tree[i].2));
                        nvert = max(nvert, tree[i].0 + 1);
                        nvert = max(nvert, tree[i].1 + 1);
                    }

                    // Make edge names.

                    let mut vnames = Vec::<String>::new();
                    for i in 0..nvert {
                        let mut len = 0.0;
                        for j in 0..edges.len() {
                            if edges[j].1 == i {
                                len = edges[j].2;
                            }
                        }
                        let mut c = String::new();
                        if i > 0 && i <= n && ctl.gen_opt.tree == "const".to_string() {
                            let ex = &exact_clonotypes[exacts[oo][i - 1]];
                            let mut h = Vec::<String>::new();
                            for m in 0..ex.share.len() {
                                if ex.share[m].left {
                                    if ex.share[m].c_ref_id.is_none() {
                                        h.push("?".to_string());
                                    } else {
                                        h.push(refdata.name[ex.share[m].c_ref_id.unwrap()].clone());
                                    }
                                }
                            }
                            c = format!(",{}", h.iter().format("+"));
                        }
                        if i == 0 {
                            vnames.push("•".to_string());
                        } else if i <= n {
                            if ctl.pretty {
                                vnames.push(format!("[01m[31m{}[0m [{:.2}{}]", i, len, c));
                            } else {
                                vnames.push(format!("{} [{:.2}{}]", i, len, c));
                            }
                        } else {
                            vnames.push(format!("• [{:.2}{}]", len, c));
                        }
                    }

                    // Display the tree.

                    let nw = display_tree(&vnames, &edges, 0, 100);
                    fwrite!(logx, "\n{}", nw);
                }
            }

            // Generate peer group output.

            if ctl.gen_opt.peer_group_filename.len() > 0 {
                let pg = mammalian_fixed_len_peer_groups(&refdata);
                if !ctl.gen_opt.peer_group_readable {
                    if ctl.gen_opt.peer_group_filename == "stdout".to_string() {
                        fwriteln!(logx, "group,clonotype,chain,pos,amino_acid,count");
                    }
                } else {
                    if ctl.gen_opt.peer_group_filename == "stdout".to_string() {
                        fwriteln!(logx, "group,clonotype,chain,pos,distribution");
                    }
                }
                let chain_types = ["IGH", "IGK", "IGL", "TRA", "TRB"];
                for q in 0..rsi[oo].mat.len() {
                    let id = rsi[oo].vids[q];
                    if ctl.gen_opt.peer_group_filename != "stdout".to_string() {
                        if !ctl.gen_opt.peer_group_readable {
                            for y in pg[id].iter() {
                                fwriteln!(
                                    pgout,
                                    "{},{},{},{},{},{},{}",
                                    groups,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    y.0,
                                    y.1 as char,
                                    y.2
                                );
                            }
                        } else {
                            let mut k = 0;
                            while k < pg[id].len() {
                                let l = next_diff1_3(&pg[id], k as i32) as usize;
                                let mut s = Vec::<String>::new();
                                for m in k..l {
                                    s.push(format!("{}={}", pg[id][m].1 as char, pg[id][m].2));
                                }
                                fwriteln!(
                                    pgout,
                                    "{},{},{},{},{},{}",
                                    groups,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    pg[id][k].0,
                                    s.iter().format(":")
                                );
                                k = l;
                            }
                        }
                    } else {
                        if !ctl.gen_opt.peer_group_readable {
                            for y in pg[id].iter() {
                                fwriteln!(
                                    logx,
                                    "{},{},{},{},{},{},{}",
                                    groups,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    y.0,
                                    y.1 as char,
                                    y.2
                                );
                            }
                        } else {
                            let mut k = 0;
                            while k < pg[id].len() {
                                let l = next_diff1_3(&pg[id], k as i32) as usize;
                                let mut s = Vec::<String>::new();
                                for m in k..l {
                                    s.push(format!("{}={}", pg[id][m].1 as char, pg[id][m].2));
                                }
                                fwriteln!(
                                    logx,
                                    "{},{},{},{},{},{}",
                                    groups,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    pg[id][k].0,
                                    s.iter().format(":")
                                );
                                k = l;
                            }
                        }
                    }
                }
            }

            // Generate FASTA output.

            if ctl.gen_opt.fasta_filename.len() > 0 {
                for (k, u) in exacts[oo].iter().enumerate() {
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_some() {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let ex = &exact_clonotypes[*u];
                            if ctl.gen_opt.fasta_filename != "stdout".to_string() {
                                fwriteln!(
                                    fout,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            } else {
                                fwriteln!(
                                    logx,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            }
                            let mut seq = ex.share[r].seq.clone();
                            let mut cid = ex.share[r].c_ref_id;
                            if cid.is_none() {
                                for l in 0..exacts[oo].len() {
                                    if rsi[oo].mat[m][l].is_some() {
                                        let r2 = rsi[oo].mat[m][l].unwrap();
                                        let ex2 = &exact_clonotypes[exacts[oo][l]];
                                        let cid2 = ex2.share[r2].c_ref_id;
                                        if cid2.is_some() {
                                            cid = cid2;
                                            break;
                                        }
                                    }
                                }
                            }
                            if cid.is_some() {
                                let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                                seq.append(&mut cseq);
                                if ctl.gen_opt.fasta_filename != "stdout".to_string() {
                                    fwriteln!(fout, "{}", strme(&seq));
                                } else {
                                    fwriteln!(logx, "{}", strme(&seq));
                                }
                            }
                        }
                    }
                }
            }

            // Generate fasta amino acid output.

            if ctl.gen_opt.fasta_aa_filename.len() > 0 {
                for (k, u) in exacts[oo].iter().enumerate() {
                    for m in 0..rsi[oo].mat.len() {
                        if rsi[oo].mat[m][k].is_some() {
                            let r = rsi[oo].mat[m][k].unwrap();
                            let ex = &exact_clonotypes[*u];
                            if ctl.gen_opt.fasta_aa_filename != "stdout".to_string() {
                                fwriteln!(
                                    faaout,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            } else {
                                fwriteln!(
                                    logx,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            }
                            let mut seq = ex.share[r].seq.clone();
                            let mut cid = ex.share[r].c_ref_id;
                            if cid.is_none() {
                                for l in 0..exacts[oo].len() {
                                    if rsi[oo].mat[m][l].is_some() {
                                        let r2 = rsi[oo].mat[m][l].unwrap();
                                        let ex2 = &exact_clonotypes[exacts[oo][l]];
                                        let cid2 = ex2.share[r2].c_ref_id;
                                        if cid2.is_some() {
                                            cid = cid2;
                                            break;
                                        }
                                    }
                                }
                            }
                            if cid.is_some() {
                                let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                                seq.append(&mut cseq);
                                if ctl.gen_opt.fasta_aa_filename != "stdout".to_string() {
                                    fwriteln!(faaout, "{}", strme(&aa_seq(&seq, 0)));
                                } else {
                                    fwriteln!(logx, "{}", strme(&aa_seq(&seq, 0)));
                                }
                            }
                        }
                    }
                }
            }

            // Generate parseable output.

            if ctl.parseable_opt.pout.len() > 0 {
                let mut rows = Vec::<Vec<String>>::new();
                for m in 0..out_datas[oo].len() {
                    out_datas[oo][m].insert("group_id".to_string(), format!("{}", groups));
                    out_datas[oo][m]
                        .insert("group_ncells".to_string(), format!("{}", group_ncells));
                    out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
                }
                if ctl.parseable_opt.pout == "stdout".to_string() {
                    fwriteln!(logx, "{}", pcols.iter().format(","));
                }
                if ctl.parseable_opt.pout == "stdouth".to_string() {
                    rows.push(pcols.clone());
                }
                let x = &out_datas[oo];
                for (u, y) in x.iter().enumerate() {
                    if !ctl.parseable_opt.pbarcode {
                        if ctl.parseable_opt.pout != "stdouth".to_string() {
                            for (i, c) in pcols.iter().enumerate() {
                                if i > 0 {
                                    if ctl.parseable_opt.pout != "stdout".to_string() {
                                        fwrite!(pout, ",");
                                    } else {
                                        fwrite!(logx, ",");
                                    }
                                }
                                if y.contains_key(c) {
                                    let val = &y[c];
                                    if !val.contains(',') {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "{}", val);
                                        } else {
                                            fwrite!(logx, "{}", val);
                                        }
                                    } else {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "\"{}\"", val);
                                        } else {
                                            fwrite!(logx, "\"{}\"", val);
                                        }
                                    }
                                } else {
                                    if ctl.parseable_opt.pout != "stdout".to_string() {
                                        fwrite!(pout, "");
                                    } else {
                                        fwrite!(logx, "");
                                    }
                                }
                            }
                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                fwriteln!(pout, "");
                            } else {
                                fwriteln!(logx, "");
                            }
                        } else {
                            let mut row = Vec::<String>::new();
                            for c in pcols.iter() {
                                if y.contains_key(c) {
                                    let val = &y[c];
                                    row.push(val.clone());
                                } else {
                                    row.push("".to_string());
                                }
                            }
                            rows.push(row);
                        }
                    } else {
                        let ex = &exact_clonotypes[exacts[oo][u]];
                        let n = ex.ncells();
                        if ctl.parseable_opt.pout != "stdouth".to_string() {
                            // Traverse the cells in the exact subclonotype.
                            for m in 0..n {
                                // Traverse the parseable fields to be displayed.
                                for (i, c) in pcols.iter().enumerate() {
                                    if i > 0 {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, ",");
                                        } else {
                                            fwrite!(logx, ",");
                                        }
                                    }
                                    // Test for whether the out_data contain the field.
                                    if y.contains_key(c) {
                                        let mut id = 0;
                                        let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                        if vals.len() > 1 {
                                            id = m;
                                        }
                                        if id >= vals.len() {
                                            panic!(
                                                "id >= vals.len() where id = {} and vals.len() \
                                                = {},\nparseable variable = {}, barcodes include \
                                                {}, n = {}, y[c] = {}",
                                                id,
                                                vals.len(),
                                                c,
                                                ex.clones[0][0].barcode,
                                                n,
                                                y[c],
                                            );
                                        }
                                        let val = vals[id];
                                        if !val.contains(',') {
                                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                                fwrite!(pout, "{}", val);
                                            } else {
                                                fwrite!(logx, "{}", val);
                                            }
                                        } else {
                                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                                fwrite!(pout, "\"{}\"", val);
                                            } else {
                                                fwrite!(logx, "\"{}\"", val);
                                            }
                                        }
                                    } else {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "");
                                        } else {
                                            fwrite!(logx, "");
                                        }
                                    }
                                }
                                if ctl.parseable_opt.pout != "stdout".to_string() {
                                    fwriteln!(pout, "");
                                } else {
                                    fwriteln!(logx, "");
                                }
                            }
                        } else {
                            for m in 0..n {
                                let mut row = Vec::<String>::new();
                                for c in pcols.iter() {
                                    if y.contains_key(c) {
                                        let mut id = 0;
                                        let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                        if vals.len() > 1 {
                                            id = m;
                                        }
                                        let val = vals[id];
                                        row.push(val.to_string());
                                    } else {
                                        row.push("".to_string());
                                    }
                                }
                                rows.push(row);
                            }
                        }
                    }
                }
                if ctl.parseable_opt.pout == "stdouth".to_string() {
                    let mut log = Vec::<u8>::new();
                    let mut justify = Vec::<u8>::new();
                    for x in rows[0].iter() {
                        justify.push(justification(&x));
                    }
                    print_tabular(&mut log, &rows, 2, Some(justify));
                    fwrite!(logx, "{}", strme(&log));
                }
            }
        }
    }

    // Finish CLUSTAL.

    if clustal_aa.is_some() {
        clustal_aa.unwrap().finish().unwrap();
    }
    if clustal_dna.is_some() {
        clustal_dna.unwrap().finish().unwrap();
    }

    // Compute two umi stats.

    let nclono = exacts.len();
    let mut umish = Vec::<usize>::new();
    let mut umisl = Vec::<usize>::new();
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for k in 0..ex.clones.len() {
                for l in 0..ex.share.len() {
                    if ex.share[l].left {
                        umish.push(ex.clones[k][l].umi_count);
                    } else {
                        umisl.push(ex.clones[k][l].umi_count);
                    }
                }
            }
        }
    }
    umish.sort();
    umisl.sort();
    let (mut middleh, mut denomh) = (0, 0);
    for j in umish.len() / 3..(2 * umish.len()) / 3 {
        middleh += umish[j];
        denomh += 1;
    }
    let mut middle_mean_umish = 0.0;
    if denomh > 0 {
        middle_mean_umish = (middleh as f64) / (denomh as f64);
    }
    let (mut middlel, mut denoml) = (0, 0);
    for j in umisl.len() / 3..(2 * umisl.len()) / 3 {
        middlel += umisl[j];
        denoml += 1;
    }
    let mut middle_mean_umisl = 0.0;
    if denoml > 0 {
        middle_mean_umisl = (middlel as f64) / (denoml as f64);
    }

    // Compute n1 and n2 and n23 and n4.

    let mut n1 = 0;
    let mut n2 = 0;
    let mut n23 = 0;
    let mut n4 = 0;
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            if ex.nchains() == 1 {
                n1 += ex.ncells();
            }
            if ex.nchains() == 2 {
                n2 += ex.ncells();
            }
            if ex.nchains() == 2 || ex.nchains() == 3 {
                n23 += ex.ncells();
            }
            if ex.nchains() == 4 {
                n4 += ex.ncells();
            }
        }
    }

    // Print summary stats.

    let mut ncells = 0;
    let mut nclono2 = 0;
    let mut ncc = Vec::<(usize, usize)>::new();
    let mut sd = Vec::<(Option<usize>, Option<usize>)>::new();
    for i in 0..nclono {
        let mut n = 0;
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            n += ex.ncells();
            for k in 0..ex.clones.len() {
                let x = &ex.clones[k][0];
                sd.push((x.origin_index, x.donor_index));
            }
        }
        if n >= 2 {
            nclono2 += 1;
        }
        ncells += n;
        ncc.push((rsi[i].mat.len(), n));
    }
    sd.sort();
    let mut sdx = Vec::<(Option<usize>, Option<usize>, usize)>::new();
    let mut i = 0;
    while i < sd.len() {
        let j = next_diff(&sd, i);
        sdx.push((sd[i].0, sd[i].1, j - i));
        i = j;
    }
    if ctl.gen_opt.summary {
        fwriteln!(logx, "\nSUMMARY STATISTICS");
        fwriteln!(logx, "1. overall");
        fwriteln!(logx, "   • number of datasets = {}", ctl.origin_info.n());
        fwriteln!(logx, "   • number of donors = {}", ctl.origin_info.donors);
        let mut vcells = 0;
        for i in 0..vdj_cells.len() {
            vcells += vdj_cells[i].len();
        }
        fwriteln!(logx, "   • original number of cells = {}", vcells);

        // Print mean reads per cell if known.

        let mut known = true;
        for i in 0..ctl.origin_info.n() {
            if ctl.origin_info.cells_cellranger[i].is_none() {
                known = false;
            } else if ctl.origin_info.mean_read_pairs_per_cell_cellranger[i].is_none() {
                known = false;
            }
        }
        if known {
            let (mut cells, mut read_pairs) = (0, 0);
            for i in 0..ctl.origin_info.n() {
                let c = ctl.origin_info.cells_cellranger[i].unwrap();
                let rpc = ctl.origin_info.mean_read_pairs_per_cell_cellranger[i].unwrap();
                cells += c;
                read_pairs += cells * rpc;
            }
            let rpc = ((read_pairs as f64) / (cells as f64)).round();
            fwriteln!(logx, "   • read pairs per cell = {}", rpc);
        }

        // Print computational performance stats.

        if !ctl.gen_opt.summary_clean {
            fwriteln!(
                logx,
                "   • total elapsed time = {:.1} seconds",
                elapsed(&tall)
            );
            #[cfg(not(target_os = "macos"))]
            fwriteln!(logx, "   • peak memory = {:.1} GB", peak_mem_usage_gb());
        }

        // Compute marking stats.

        let (mut nmarked, mut nmarked_good, mut ndubious) = (0, 0, 0);
        let (mut nfake, mut nfake_marked, mut ngood, mut ngood_marked) = (0, 0, 0, 0);
        if ctl.gen_opt.mark_stats || ctl.gen_opt.mark_stats2 {
            for i in 0..nclono {
                let mut datasets = Vec::<usize>::new();
                let mut ncells = 0;
                for j in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][j]];
                    ncells += ex.ncells();
                    for l in 0..ex.ncells() {
                        datasets.push(ex.clones[l][0].dataset_index);
                    }
                }
                datasets.sort();
                let mut freq = Vec::<(u32, usize)>::new();
                make_freq(&datasets, &mut freq);
                let mut fake = false;
                let mut di = -1;
                if freq.len() == 1 || freq[0].0 >= 10 * freq[1].0 {
                    di = freq[0].1 as isize;
                }
                if ncells >= 10 && di >= 0 {
                    nfake += ncells - 1;
                    fake = true;
                }
                if ncells >= 2 {
                    for j in 0..exacts[i].len() {
                        let ex = &exact_clonotypes[exacts[i][j]];

                        // Determine if cell is called a B cell.

                        for l in 0..ex.ncells() {
                            let mut b = false;
                            let li = ex.clones[l][0].dataset_index;
                            let bc = &ex.clones[l][0].barcode;
                            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                                    b = true;
                                }
                            }

                            // Record accordingly.

                            if ex.clones[l][0].dataset_index as isize == di || !b {
                                ndubious += 1;
                            }
                            if !fake
                                && ncells >= 10
                                && b
                                && (ex.share.len() == 2 || ex.share.len() == 3)
                            {
                                ngood += 1;
                                if ex.clones[l][0].marked {
                                    ngood_marked += 1;
                                }
                            }
                        }
                    }
                }
                if fake {
                    ngood += 1;
                }
                let mut fake_marks = 0;
                for j in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][j]];
                    for l in 0..ex.ncells() {
                        if ex.clones[l][0].marked {
                            nmarked += 1;
                            if fake {
                                nfake_marked += 1;
                                fake_marks += 1;
                            }
                            let chains_ok = ex.nchains() >= 2 && ex.nchains() <= 3;
                            let mut b = false;
                            let li = ex.clones[l][0].dataset_index;
                            let bc = &ex.clones[l][0].barcode;
                            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                                    b = true;
                                }
                            }
                            if chains_ok && freq.len() >= 2 && b {
                                nmarked_good += 1;
                            }
                        }
                    }
                }
                if fake_marks == ncells {
                    nfake_marked -= 1;
                }
            }
        }

        // Print barcode fate.

        fwriteln!(logx, "2. barcode fate");
        let mut fates = Vec::<String>::new();
        for i in 0..fate.len() {
            for f in fate[i].iter() {
                if f.1.contains(" GEX ") && ctl.clono_filt_opt.ngex {
                    continue;
                }
                if f.1.contains(" CROSS ") && ctl.clono_filt_opt.ncross {
                    continue;
                }
                if f.1.contains(" UMI ") && !ctl.clono_filt_opt.umi_filt {
                    continue;
                }
                if f.1.contains(" UMI_RATIO ") && !ctl.clono_filt_opt.umi_ratio_filt {
                    continue;
                }
                if f.1.contains(" GRAPH_FILTER ") && ctl.gen_opt.ngraph_filter {
                    continue;
                }
                if f.1.contains(" QUAL") && !ctl.clono_filt_opt.qual_filter {
                    continue;
                }
                if f.1.contains(" WEAK_CHAINS ") && !ctl.clono_filt_opt.weak_chains {
                    continue;
                }
                if f.1.contains(" FOURSIE_KILL ") && !ctl.clono_filt_opt.weak_foursies {
                    continue;
                }
                if f.1.contains(" WHITEF ") && ctl.gen_opt.nwhitef {
                    continue;
                }
                if f.1.contains(" BC_DUP ") && !ctl.clono_filt_opt.bc_dup {
                    continue;
                }
                if f.1.contains(" IMPROPER ") && ctl.merge_all_impropers {
                    continue;
                }
                fates.push(f.1.clone());
            }
        }
        fates.sort();
        let mut freq = Vec::<(u32, String)>::new();
        make_freq(&fates, &mut freq);
        let mut rows = Vec::<Vec<String>>::new();
        rows.push(vec!["barcodes".to_string(), "why deleted".to_string()]);
        rows.push(vec!["\\hline".to_string(); 2]);
        for i in 0..freq.len() {
            rows.push(vec![format!("{}", freq[i].0), freq[i].1.clone()]);
        }
        rows.push(vec![format!("{}", fates.len()), "total".to_string()]);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"r|l".to_vec(), false, false);
        log.truncate(log.len() - 1);
        log = log.replace("\n", "\n   ");
        fwrite!(logx, "   {}\n", log);

        // Print other stats.

        fwriteln!(logx, "3. for the selected clonotypes");

        // Print summary table for chains / clonotypes / cells.

        ncc.sort();
        let mut i = 0;
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "chains".to_string(),
            "clonotypes with this".to_string(),
            "cells in these".to_string(),
            "%".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "".to_string(),
            "number of chains".to_string(),
            "clonotypes".to_string(),
            "".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 4];
        rows.push(row);
        while i < ncc.len() {
            let j = next_diff1_2(&ncc, i as i32) as usize;
            let nchains_this = ncc[i].0;
            let nclono_this = j - i;
            let mut ncells_this = 0;
            for k in i..j {
                ncells_this += ncc[k].1;
            }
            let row = vec![
                format!("{}", nchains_this),
                format!("{}", nclono_this),
                format!("{}", ncells_this),
                format!("{:.1}", percent_ratio(ncells_this, ncells)),
            ];
            rows.push(row);
            i = j;
        }
        let row = vec![
            "total".to_string(),
            format!("{}", nclono),
            format!("{}", ncells),
            "100.0".to_string(),
        ];
        rows.push(row);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|r|r|r".to_vec(), false, false);
        log = log.replace("\n", "\n   ");
        fwrite!(logx, "   {}", log);

        // Print other cell/clonotype stats.

        fwriteln!(
            logx,
            "• number of clonotypes having at least two cells = {}",
            nclono2
        );
        fwriteln!(logx, "   • number of cells having 1 chain = {}", n1);
        fwriteln!(logx, "   • number of cells having 2 or 3 chains = {}", n23);
        let mut doublet_rate = 0.0;
        if n2 > 0 || n4 > 0 {
            doublet_rate = n4 as f64 / (n2 + n4) as f64;
        }
        let celltype;
        if ctl.gen_opt.bcr {
            celltype = "B";
        } else {
            celltype = "T";
        }
        fwrite!(
            logx,
            "   • estimated {}-{} doublet rate = {:.1}% = {}/{}",
            celltype,
            celltype,
            100.0 * doublet_rate,
            n4,
            n2 + n4
        );
        fwriteln!(logx, " = cells with 4 chains / cells with 2 or 4 chains");

        // Print UMI stats.

        let hchain;
        let lchain;
        if ctl.gen_opt.bcr {
            hchain = "heavy chain";
            lchain = "light chain";
        } else {
            hchain = "TRB";
            lchain = "TRA";
        }
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts ({}) = {:.2}",
            hchain,
            middle_mean_umish,
        );
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts ({}) = {:.2}",
            lchain,
            middle_mean_umisl,
        );

        // Print marking stats.

        if ctl.gen_opt.mark_stats {
            fwriteln!(logx, "   --------------------------------");
            fwriteln!(logx, "   • number of dubious cells = {}", ndubious);
            fwriteln!(logx, "   • number of marked cells = {}", nmarked);
            fwriteln!(logx, "   • number of good marked cells = {}", nmarked_good);
        }
        if ctl.gen_opt.mark_stats2 {
            fwriteln!(logx, "   --------------------------------");
            fwriteln!(
                logx,
                "   • number of fake expanded clonotype cells = {}",
                nfake
            );
            fwriteln!(
                logx,
                "   • number of these that are marked = {}",
                nfake_marked
            );
            fwriteln!(logx, "   • residual = {}", nfake - nfake_marked);
            fwriteln!(
                logx,
                "   • number of good expanded clonotype cells = {}",
                ngood
            );
            fwriteln!(
                logx,
                "   • number of these that are marked = {}",
                ngood_marked
            );
        }

        // Print origin (sample)/donor table, but only if there is more than one.

        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "origin".to_string(),
            "donor".to_string(),
            "cells".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 3];
        rows.push(row);
        for i in 0..sdx.len() {
            let mut row = Vec::<String>::new();
            if sdx[i].0.is_some() {
                row.push(format!(
                    "{}",
                    ctl.origin_info.origin_list[sdx[i].0.unwrap()]
                ));
            } else {
                row.push("?".to_string());
            }
            if sdx[i].1.is_some() {
                row.push(format!("{}", ctl.origin_info.donor_list[sdx[i].1.unwrap()]));
            } else {
                row.push("?".to_string());
            }
            row.push(format!("{}", sdx[i].2));
            rows.push(row);
        }
        if rows.len() > 3 {
            let mut log = String::new();
            print_tabular_vbox(&mut log, &rows, 2, &b"llr".to_vec(), false, false);
            log = log.replace("\n", "\n   ");
            fwrite!(logx, "   {}", log);
        }

        // Print gene expression cluster analysis, but only if there is one library and the
        // cell types file was provided (which would be for PD cellranger runs with GEX data).

        if ctl.origin_info.n() == 1 && gex_info.cell_type_specified[0] {
            let cell_type = &gex_info.cell_type[0];
            let cluster = &gex_info.cluster[0];
            if cluster.len() > 0 {
                let mut nclust = 0;
                for x in cluster.iter() {
                    nclust = max(*x.1, nclust);
                }
                let mut ncells = vec![0; nclust];
                let mut ncells_type = vec![0; nclust];
                let mut ncells_vdj = vec![0; nclust];
                let mut ncells_vdj_max = vec![0; nclust];
                for x in cluster.iter() {
                    ncells[x.1 - 1] += 1;
                    let bc = &x.0;
                    let typex = &cell_type[bc.clone()];
                    if ctl.gen_opt.bcr && typex.starts_with('B') {
                        ncells_type[x.1 - 1] += 1;
                    } else if ctl.gen_opt.tcr && typex.starts_with('T') {
                        ncells_type[x.1 - 1] += 1;
                    }
                }
                for i in 0..nclono {
                    let mut cs = Vec::<usize>::new();
                    for j in 0..exacts[i].len() {
                        let ex = &exact_clonotypes[exacts[i][j]];
                        for k in 0..ex.ncells() {
                            let bc = &ex.clones[k][0].barcode;
                            if cluster.contains_key(&bc.clone()) {
                                cs.push(cluster[&bc.clone()] - 1);
                            }
                        }
                    }
                    cs.sort();
                    let mut freq = Vec::<(u32, usize)>::new();
                    make_freq(&cs, &mut freq);
                    for x in freq.iter() {
                        ncells_vdj[x.1] += x.0;
                        ncells_vdj_max[x.1] = max(x.0, ncells_vdj_max[x.1]);
                    }
                }
                let mut rows = Vec::<Vec<String>>::new();
                let bt;
                if ctl.gen_opt.bcr {
                    bt = "B";
                } else {
                    bt = "T";
                }
                rows.push(vec![
                    "gex cluster".to_string(),
                    "cells".to_string(),
                    format!("%{}", bt),
                    "%in clonotypes".to_string(),
                    "max % in one clonotype".to_string(),
                ]);
                rows.push(vec!["\\hline".to_string(); 5]);
                for i in 0..nclust {
                    let n = ncells[i];
                    let mut row = Vec::<String>::new();
                    row.push(format!("{}", i + 1));
                    row.push(format!("{}", n));
                    row.push(format!("{:.1}", percent_ratio(ncells_type[i], n)));
                    row.push(format!("{:.1}", percent_ratio(ncells_vdj[i] as usize, n)));
                    row.push(format!(
                        "{:.1}",
                        percent_ratio(ncells_vdj_max[i] as usize, n)
                    ));
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 2, &b"r|r|r|r|r".to_vec(), false, false);
                log = log.replace("\n", "\n   ");
                fwrite!(logx, "   {}", log);
            }
        }
    }

    // Print summary csv stats.

    if ctl.gen_opt.summary_csv {
        println!("\nmiddle_mean_umis_heavy,middle_mean_umis_light,n_twothreesie");
        println!("{:.2},{:.2},{}", middle_mean_umish, middle_mean_umisl, n23);
    }

    // Print to stdout.

    if !ctl.gen_opt.html {
        print!("{}", compress_ansi_escapes(&strme(&logx)));
    } else {
        // Note that we do not link to the css file, because it is less fragile then including
        // the font face information directly.  In particular, the css file could be accidentally
        // deleted or renamed, which would break previously generated user html files.  This
        // actually happened!
        let s = convert_text_with_ansi_escapes_to_html(
            strme(&logx),
            "", // source
            &ctl.gen_opt.html_title,
            &format!("<style type=\"text/css\">\n{}</style>", font_face_in_css()),
            "DejaVuSansMono",
            14,
        );
        print!("{}", s);
    }

    // Test for required number of false positives.

    if ctl.gen_opt.required_fps.is_some() {
        let mut fps = 0;
        for i in 0..pics.len() {
            if pics[i].contains("WARNING:") {
                fps += 1;
            }
        }
        if fps != ctl.gen_opt.required_fps.unwrap() {
            eprintln!(
                "\nA \"false positive\" is a clonotype that contains cells from multiple\n\
                 donors.  You invoked enclone with the argument REQUIRED_FPS={}, but we found\n\
                 {} false positives, so the requirement is not met.\n",
                ctl.gen_opt.required_fps.unwrap(),
                fps
            );
            std::process::exit(1);
        }
    }

    // Test for required number of cells.

    if ctl.gen_opt.required_cells.is_some() {
        let mut ncells = 0;
        for i in 0..pics.len() {
            for x in exacts[i].iter() {
                ncells += exact_clonotypes[*x].ncells();
            }
        }
        if ctl.gen_opt.required_cells.unwrap() != ncells {
            eprintln!(
                "\nThe required number of cells is {}, but you actually have {}.\n",
                ctl.gen_opt.required_cells.unwrap(),
                ncells,
            );
            std::process::exit(1);
        }
    }

    // Test for required number of clonotypes.

    if ctl.gen_opt.required_clonotypes.is_some() {
        if ctl.gen_opt.required_clonotypes.unwrap() != nclono {
            eprintln!(
                "\nThe required number of clonotypes is {}, but you actually have {}.\n",
                ctl.gen_opt.required_clonotypes.unwrap(),
                nclono,
            );
            std::process::exit(1);
        }
    }

    // Test for required number of donors

    if ctl.gen_opt.required_donors.is_some() {
        let ndonors = ctl.origin_info.donors;
        if ctl.gen_opt.required_donors.unwrap() != ndonors {
            eprintln!(
                "\nThe required number of donors is {}, but you actually have {}.\n",
                ctl.gen_opt.required_donors.unwrap(),
                ndonors,
            );
            std::process::exit(1);
        }
    }

    // Test for required number of >=2 cell clonotypes.

    if ctl.gen_opt.required_two_cell_clonotypes.is_some() {
        if ctl.gen_opt.required_two_cell_clonotypes.unwrap() != nclono2 {
            eprintln!(
                "\nThe required number of two-cell clonotypes is {}, but you actually have {}.\n",
                ctl.gen_opt.required_two_cell_clonotypes.unwrap(),
                nclono2,
            );
            std::process::exit(1);
        }
    }

    // Test for required number of datasets

    if ctl.gen_opt.required_datasets.is_some() {
        let ndatasets = ctl.origin_info.n();
        if ctl.gen_opt.required_datasets.unwrap() != ndatasets {
            eprintln!(
                "\nThe required number of datasets is {}, but you actually have {}.\n",
                ctl.gen_opt.required_datasets.unwrap(),
                ndatasets,
            );
            std::process::exit(1);
        }
    }
}
