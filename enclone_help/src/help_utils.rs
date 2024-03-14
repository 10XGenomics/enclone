// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::ansi_to_html::convert_text_with_ansi_escapes_to_html;
use ansi_escape::{
    emit_blue_escape, emit_bold_escape, emit_end_escape, emit_green_escape, emit_red_escape,
    print_color,
};
use enclone_core::print_tools::emit_codon_color_escape;
use io_utils::fwrite;
use std::io::Write;
use string_utils::{stringme, strme};
use tables::print_tabular_vbox;

pub fn font_face_in_css() -> String {
    let f = include_str!["enclone_css_v2.css"];
    let mut x = String::new();
    let mut in_font_face = false;
    let mut count = 0;
    for line in f.lines() {
        if line.starts_with("@font-face") {
            in_font_face = true;
            count += 1;
        }
        if in_font_face {
            x += &format!("{}\n", line);
        }
        if line == "}" {
            in_font_face = false;
        }
    }
    assert_eq!(count, 2); // because there are two fonts: regular and bold
    x
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
pub struct HelpDesk {
    pub plain: bool,
    pub help_all: bool,
    pub long_help: bool,
    pub html: bool,
    pub rows: Vec<Vec<String>>,
    pub log: Vec<u8>,
    pub title: String,
    pub ok: bool,
}

impl HelpDesk {
    pub fn new(plain: bool, help_all: bool, long_help: bool, html: bool) -> HelpDesk {
        HelpDesk {
            plain,
            help_all,
            long_help,
            html,
            rows: Vec::<Vec<String>>::new(),
            log: Vec::<u8>::new(),
            title: String::new(),
            ok: false,
        }
    }
    pub fn doc(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec![x1.to_string(), x2.to_string()]);
    }
    pub fn doc2(&mut self, x2: &str) {
        self.rows.push(vec!["".to_string(), x2.to_string()]);
    }

    // docf2: like doc, but fold x2 to n2 chars.

    pub fn docf2(&mut self, x1: &str, x2: &str, n2: usize) -> Result<(), String> {
        let mut c2 = Vec::<char>::new();
        for m in x2.chars() {
            c2.push(m);
        }
        let mut y2 = Vec::<String>::new();
        let mut start = 0;
        let mut i = 0;
        while i < c2.len() {
            if c2[i] == ' ' && i - start > n2 {
                i -= 1;
                while i > 0 && c2[i] != ' ' {
                    i -= 1;
                }
                let mut s = String::new();
                for j in start..i {
                    s.push(c2[j]);
                }
                if s.is_empty() {
                    return Err(format!(
                        "\nError in docf2, x2 = \n\n{}\n\nThere is probably a string in there that \
                        exceeds your argument {}.\n",
                        x2,
                        n2,
                    ));
                }
                y2.push(s);
                start = i + 1;
            }
            i += 1;
        }
        let mut s = String::new();
        for j in start..c2.len() {
            s.push(c2[j]);
        }
        y2.push(s);
        for i in 0..y2.len() {
            if i == 0 {
                self.doc(x1, &y2[i]);
            } else {
                self.doc("", &y2[i]);
            }
        }
        Ok(())
    }

    pub fn ldoc(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec!["\\hline".to_string(); 2]);
        self.rows.push(vec![x1.to_string(), x2.to_string()]);
    }
    pub fn doc3(&mut self, x1: &str, x2: &str, x3: &str) {
        self.rows.push(vec![
            self.print_to(x1),
            self.print_to(x2),
            self.print_to(x3),
        ]);
    }
    pub fn ldoc3(&mut self, x1: &str, x2: &str, x3: &str) {
        self.rows.push(vec!["\\hline".to_string(); 3]);
        self.rows
            .push(vec![x1.to_string(), x2.to_string(), x3.to_string()]);
    }
    pub fn docpr(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec![self.print_to(x1), self.print_to(x2)]);
    }
    pub fn ldocpr(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec!["\\hline".to_string(); 2]);
        self.rows.push(vec![self.print_to(x1), self.print_to(x2)]);
    }
    pub fn ldoc3pr(&mut self, x1: &str, x2: &str, x3: &str) {
        self.rows.push(vec!["\\hline".to_string(); 3]);
        self.rows.push(vec![
            self.print_to(x1),
            self.print_to(x2),
            self.print_to(x3),
        ]);
    }
    pub fn doc_red(&mut self, x1: &str, x2: &str) {
        if !self.plain {
            let r1 = format!("[31m{}[0m", x1);
            let r2 = format!("[31m{}[0m", x2);
            self.rows.push(vec![r1, r2]);
        } else {
            self.rows.push(vec![x1.to_string(), x2.to_string()]);
        }
    }
    pub fn ldoc_red(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec!["\\hline".to_string(); 2]);
        if !self.plain {
            let r1 = format!("[31m{}[0m", x1);
            let r2 = format!("[31m{}[0m", x2);
            self.rows.push(vec![r1, r2]);
        } else {
            self.rows.push(vec![x1.to_string(), x2.to_string()]);
        }
    }
    pub fn doc_greenish(&mut self, x1: &str, x2: &str) {
        if !self.plain {
            let r1 = format!("[38;5;36m{}[0m", x1);
            let r2 = format!("[38;5;36m{}[0m", x2);
            self.rows.push(vec![r1, r2]);
        } else {
            self.rows.push(vec![x1.to_string(), x2.to_string()]);
        }
    }
    pub fn ldoc_greenish(&mut self, x1: &str, x2: &str) {
        self.rows.push(vec!["\\hline".to_string(); 2]);
        if !self.plain {
            let r1 = format!("[38;5;36m{}[0m", x1);
            let r2 = format!("[38;5;36m{}[0m", x2);
            self.rows.push(vec![r1, r2]);
        } else {
            self.rows.push(vec![x1.to_string(), x2.to_string()]);
        }
    }
    pub fn print_enclone(&mut self) -> Result<(), String> {
        if self.plain {
            self.print("enclone")?;
        } else {
            let mut log = Vec::<u8>::new();
            print_color(3, &mut log);
            log.push(b'e');
            emit_end_escape(&mut log);
            print_color(1, &mut log);
            log.push(b'n');
            emit_end_escape(&mut log);
            print_color(2, &mut log);
            log.push(b'c');
            emit_end_escape(&mut log);
            print_color(0, &mut log);
            log.push(b'l');
            emit_end_escape(&mut log);
            print_color(4, &mut log);
            log.push(b'o');
            emit_end_escape(&mut log);
            print_color(5, &mut log);
            log.push(b'n');
            emit_end_escape(&mut log);
            print_color(1, &mut log);
            log.push(b'e');
            emit_end_escape(&mut log);
            self.print(strme(&log))?;
        }
        Ok(())
    }
    pub fn print_tab2(&mut self) -> Result<(), String> {
        let mut log = String::new();
        print_tabular_vbox(&mut log, &self.rows, 2, b"l|l".as_ref(), false, false);
        self.print_plain(&log.to_string())?;
        Ok(())
    }
    pub fn print_tab3(&mut self) -> Result<(), String> {
        let mut log = String::new();
        print_tabular_vbox(&mut log, &self.rows, 2, b"l|l|l".as_ref(), false, false);
        self.print_plain(&log.to_string())?;
        Ok(())
    }
    pub fn begin_doc(&mut self, title: &str) -> Result<(), String> {
        self.title = format!("enclone help {}", title);
        self.rows.clear();
        if self.help_all {
            let mut log = Vec::<u8>::new();
            if !self.plain {
                emit_blue_escape(&mut log);
            }
            self.print_plain(strme(&log))?;
            for _ in 1..100 {
                self.print_plain("▓")?;
            }
            self.print_plain(strme(&log))?;
            if title.is_empty() {
                self.print_plain(
                    "\nenclone main help page (what you get by typing \
                    \"enclone\")\n",
                )?;
            } else if title == "setup" {
                self.print_plain(
                    "\nenclone setup page (for one time use, what you get by typing \
                    \"enclone help\")\n",
                )?;
            } else {
                self.print_plain(&format!("\nenclone help {}\n", title))?;
            }
            let mut log = Vec::<u8>::new();
            if !self.plain {
                emit_blue_escape(&mut log);
            }
            self.print_plain(strme(&log))?;
            for _ in 1..100 {
                self.print_plain("▓")?;
            }
            let mut log = Vec::<u8>::new();
            if !self.plain {
                emit_end_escape(&mut log);
            }
            self.print_plain(&format!("{}\n", strme(&log)))?;
        }
        Ok(())
    }
    pub fn end_doc(&mut self) {
        if !self.help_all {
            self.dump();
        }
        self.ok = true;
    }
    pub fn print_with_box(&mut self, x: &str, bold_box: bool) -> Result<(), String> {
        let y = self.print_to(x);
        let mut rows = Vec::<Vec<String>>::new();
        let lines = y.split('\n').collect::<Vec<&str>>();
        for z in lines {
            rows.push(vec![z.to_string()]);
        }
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, b"l".as_ref(), false, bold_box);
        self.print_plain(&format!("{}\n", log))?;
        Ok(())
    }
    pub fn print(&mut self, x: &str) -> Result<(), String> {
        self.print_plain(&self.print_to(x))?;
        Ok(())
    }
    pub fn print_plain_unchecked(&mut self, x: &str) {
        fwrite!(self.log, "{}", &x);
    }
    pub fn print_plain(&mut self, x: &str) -> Result<(), String> {
        if !self.long_help {
            let mut count = 0;
            let mut escaped = false;
            let mut line = String::new();
            for c in x.chars() {
                line.push(c);
                if c == '\n' {
                    count = 0;
                    line.clear();
                } else if escaped {
                    if c == 'm' {
                        escaped = false;
                    }
                } else if c == '' {
                    escaped = true;
                } else {
                    count += 1;
                    if count > 100 {
                        return Err(format!(
                            "\nHelp line is too long:\n\n{}\n\
                            \nTry running with LONG_HELP to locate the problem.\n",
                            line
                        ));
                    }
                }
            }
        }
        fwrite!(self.log, "{}", &x);
        Ok(())
    }
    pub fn dump(&self) {
        if !self.html {
            print!("{}", strme(&self.log));
        } else {
            // Note that we do not link to the css file, because it is less fragile then including
            // the font face information directly.  In particular, the css file could be
            // accidentally deleted or renamed, which would break previously generated user html
            // files.  This actually happened!
            let s = convert_text_with_ansi_escapes_to_html(
                strme(&self.log),
                "", // source
                &self.title,
                &format!("<style type=\"text/css\">\n{}</style>", font_face_in_css()),
                "DejaVuSansMono",
                14,
            );
            print!("{}", s);
        }
    }
}

// This encodes the color codes for each possible codon of a given amino acid
// that could be found in a BCR or TCR sequence.

pub fn colored_codon_table(plainx: bool) -> String {
    let plain = b"\
        Alanine        A  GCT GCC GCA GCG\n\
        Arginine       R  CGT CGC CGA CGG AGA AGG\n\
        Asparagine     N  AAT AAC\n\
        Aspartic Acid  D  GAT GAC\n\
        Cysteine       C  TGT TGC\n\
        Glutamine      Q  CAA CAG\n\
        Glutamic Acid  E  GAA GAG\n\
        Glycine        G  GGT GGC GGA GGG\n\
        Histidine      H  CAT CAC\n\
        Isoleucine     I  ATT ATC ATA\n\
        Leucine        L  TTA TTG CTT CTC CTA CTG\n\
        Lysine         K  AAA AAG\n\
        Methionine     M  ATG\n\
        Phenylalanine  F  TTT TTC\n\
        Proline        P  CCT CCC CCA CCG\n\
        Serine         S  TCT TCC TCA TCG AGT AGC\n\
        Threonine      T  ACT ACC ACA ACG\n\
        Tryptophan     W  TGG\n\
        Tyrosine       Y  TAT TAC\n\
        Valine         V  GTT GTC GTA GTG";
    let mut colored = Vec::<u8>::new();
    let mut p = 0;
    while p < plain.len() {
        if (plain[p] as char).is_uppercase() && (plain[p + 1] as char).is_uppercase() {
            let mut log = Vec::<u8>::new();
            if !plainx {
                emit_codon_color_escape(&plain[p..p + 3], &mut log);
            }
            for i in 0..3 {
                log.push(plain[p + i]);
            }
            if !plainx {
                emit_end_escape(&mut log);
            }
            colored.append(&mut log);
            p += 3;
        } else {
            colored.push(plain[p]);
            p += 1;
        }
    }
    stringme(&colored)
}

impl HelpDesk {
    /// Print a string, making the following conversions, the first three of which are governed
    /// by the state of self.plain:
    /// • Change \bold{x} into a bolded string by issuing appropriate escape characters.
    /// • Change \red{x} into a red string by issuing appropriate escape characters.
    /// • Change \boldred{x} into a bold red string by issuing appropriate escape characters.
    /// • Fold at 99 characters.
    fn print_to(&self, x: &str) -> String {
        let mut y = Vec::<char>::new();
        for c in x.chars() {
            y.push(c);
        }
        let mut s = String::new();
        let mut i = 0;
        while i < y.len() {
            if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', '{']) {
                let mut j = i + 6;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_bold_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 6..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                }
            } else if y[i..].starts_with(&['\\', 'r', 'e', 'd', '{']) {
                let mut j = i + 5;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_red_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 5..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                } else {
                    i += 1;
                }
            } else if y[i..].starts_with(&['\\', 'b', 'l', 'u', 'e', '{']) {
                let mut j = i + 6;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_blue_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 6..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                } else {
                    i += 1;
                }
            } else if y[i..].starts_with(&['\\', 'g', 'r', 'e', 'e', 'n', '{']) {
                let mut j = i + 7;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_green_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 7..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                } else {
                    i += 1;
                }
            } else if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', 'r', 'e', 'd', '{']) {
                let mut j = i + 9;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_bold_escape(&mut log);
                        emit_red_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 9..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                } else {
                    i += 1;
                }
            } else if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', 'b', 'l', 'u', 'e', '{']) {
                let mut j = i + 10;
                while j < y.len() {
                    if y[j] == '}' {
                        break;
                    }
                    j += 1;
                }
                if j < y.len() {
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_bold_escape(&mut log);
                        emit_blue_escape(&mut log);
                    }

                    s += strme(&log);
                    for k in i + 10..j {
                        s.push(y[k]);
                    }
                    let mut log = Vec::<u8>::new();

                    if !self.plain {
                        emit_end_escape(&mut log);
                    }

                    s += strme(&log);
                    i = j + 1;
                } else {
                    i += 1;
                }
            } else {
                s.push(y[i]);
                i += 1;
            }
        }
        let mut x = Vec::<char>::new();
        for c in s.chars() {
            x.push(c);
        }
        let mut printed = 0;
        let mut escaped = false;
        let mut y = Vec::<char>::new();
        let mut i = 0;
        while i < x.len() {
            if x[i] == '' {
                escaped = true;
            }
            if escaped {
                if x[i] == 'm' {
                    escaped = false;
                }
                y.push(x[i]);
                i += 1;
                continue;
            }
            if x[i] == ' ' {
                let mut j = i + 1;
                while j < x.len() {
                    if x[j] == ' ' || x[j] == '\n' || x[j] == '' {
                        break;
                    }
                    j += 1;
                }
                if printed + j - i >= 100 - 1 {
                    y.push('\n');
                    printed = 0;
                    i += 1;
                    continue;
                }
            }
            y.push(x[i]);
            printed += 1;
            if x[i] == '\n' {
                printed = 0;
            }
            i += 1;
        }
        let mut ans = String::new();
        for i in 0..y.len() {
            ans.push(y[i]);
        }
        ans
    }
}

pub fn print_tab2(rows: &Vec<Vec<String>>) {
    let mut log = String::new();
    print_tabular_vbox(&mut log, rows, 2, b"l|l".as_ref(), false, false);
    print!("{}", log);
}

// Given a string, preface every line in in by a gray left bar.

pub fn gray_left_bar(s: &str, plain: bool) -> String {
    let mut gray = "[47m [0m ".to_string();
    if plain {
        gray = "┃ ".to_string();
    }

    let mut x = Vec::<char>::new();
    for c in s.chars() {
        x.push(c);
    }
    let mut t = gray.to_string();
    for i in 0..x.len() - 1 {
        t.push(x[i]);
        if x[i] == '\n' {
            t += &gray;
        }
    }
    t.push(x[x.len() - 1]);
    t
}

// Function that provides an explanation used for both enclone help lvars and
// enclone help cvars.
pub fn explain_alt_versions(h: &mut HelpDesk) -> Result<(), String> {
    h.print(&gray_left_bar(
        &h.print_to(
            "\\red{●} These variables have some alternate versions, \
             as shown in the table below.\n\n",
        ),
        h.plain,
    ))?;
    let mut rows = Vec::<Vec<String>>::new();
    let row = vec![
        "variable".to_string(),
        "semantics".to_string(),
        "visual".to_string(),
        "visual".to_string(),
        "parseable".to_string(),
        "parseable".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "".to_string(),
        "".to_string(),
        "".to_string(),
        "(one cell)".to_string(),
        "".to_string(),
        "(one cell)".to_string(),
    ];
    rows.push(row);
    let row = vec!["\\hline".to_string(); 6];
    rows.push(row);
    let row = vec![
        "x".to_string(),
        "median over cells".to_string(),
        "yes".to_string(),
        "this cell".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_mean".to_string(),
        "mean over cells".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_μ".to_string(),
        "(same as above)".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_sum".to_string(),
        "sum over cells".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_Σ".to_string(),
        "(same as above)".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_min".to_string(),
        "min over cells".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_max".to_string(),
        "max over cells".to_string(),
        "yes".to_string(),
        "null".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_%".to_string(),
        "% of total GEX (genes only)".to_string(),
        "yes".to_string(),
        "this cell".to_string(),
        "yes".to_string(),
        "yes".to_string(),
    ];
    rows.push(row);
    let row = vec![
        "x_cell".to_string(),
        "this cell".to_string(),
        "no".to_string(),
        "no".to_string(),
        "no".to_string(),
        "this cell".to_string(),
    ];
    rows.push(row);
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, b"l|l|l|l|l|l".as_ref(), false, false);
    h.print_plain(&gray_left_bar(&log, h.plain))?;
    h.print_plain(&gray_left_bar(
        &h.print_to(
            "Some explanation is required.  If you use enclone without certain options, you \
         get the \"visual\" column.\n\
         • Add the option \\bold{PER_CELL} \
         (see \"enclone help display\") and then you get visual output with extra lines for \
         each cell within an exact subclonotype, and each of those extra lines is described by \
         the \"visual (one cell)\" column.\n\
         • If you generate parseable output (see \"enclone help parseable\"), then you get \
         the \"parseable\" column for that output, unless you specify \\bold{PCELL}, \
         and then you get the last column.\n\
         • For the forms with μ and Σ, the Greek letters are only used in column headings for \
         visual output (to save space), and optionally, in names of fields on the command \
         line.\n\
         \\green{▶} If you try out these features, you'll see exactly what happens! \
         \\green{◀}\n",
        ),
        h.plain,
    ))?;
    Ok(())
}
