// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Utility for inserting html files.  It also changes all instances of #enclone to
// a preset format for that.

use io_utils::*;
use stats_utils::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn gtag() -> String {
    "<!-- Global site tag (gtag.js) - Google Analytics -->\n\
    <script async src=\"https://www.googletagmanager.com/gtag/js?id=UA-58278925-3\"></script>\n\
    <script>\n\
      window.dataLayer = window.dataLayer || [];\n\
      function gtag(){{dataLayer.push(arguments);}}\n\
      gtag('js', new Date());\n\
      gtag('config', 'UA-58278925-3');\n\
    </script>\n"
        .to_string()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Insert google tag and banner.

pub fn edit_html(html: &str) -> String {
    let mut lines2 = Vec::<String>::new();
    for line in html.lines() {
        if line == "</head>" {
            let g = gtag();
            for x in g.lines() {
                lines2.push(x.to_string());
            }
        }
        lines2.push(line.to_string());
        if line == "<body>" {
            lines2.push("".to_string());
            lines2.push("<br>".to_string());
            lines2.push(
                "<img src=\"../../img/enclone_banner.png\" \
                alt=\"enclone banner\" title=\"enclone banner\" width=100% />"
                    .to_string(),
            );
        }
    }
    let mut x = String::new();
    for i in 0..lines2.len() {
        x += &format!("{}\n", lines2[i]);
    }
    x
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn html_header(level: usize, title: &str, extra_head: &str) -> String {
    assert!(level == 0 || level == 2);
    let ltext;
    if level == 0 {
        ltext = "pages";
    } else {
        ltext = "..";
    }
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\
        <!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \n\
        \"https://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n\
        <!--  -->\n\
        <html xmlns=\"http://www.w3.org/1999/xhtml\">\n\
        <head>\n\
        <meta http-equiv=\"Content-Type\" content=\"application/xml+xhtml; charset=UTF-8\"/>\n\
        <title>{}</title>\n\
        <link rel=\"stylesheet\" type=\"text/css\" href=\"{}/enclone_css_v2.css\">\n\
        {}
        {}
        </head>\n",
        title,
        ltext,
        gtag(),
        extra_head,
    )
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn insert_html(in_file: &str, out_file: &str, up: bool, level: usize) {
    // Extract requirements for the big test from the file enclone.test.

    let x = include_str!("enclone.test");
    let mut required_fps = None;
    let mut required_cells = None;
    let mut required_donors = None;
    let mut required_two_cell_clonotypes = None;
    let mut required_datasets = None;
    let mut required_clonotypes = None;
    for line in x.lines() {
        if line.contains("REQUIRED_FPS=") {
            required_fps = Some(line.between("=", " ").force_usize());
        } else if line.contains("REQUIRED_CELLS=") {
            required_cells = Some(line.between("=", " ").force_usize());
        } else if line.contains("REQUIRED_CLONOTYPES=") {
            required_clonotypes = Some(line.between("=", " ").force_usize());
        } else if line.contains("REQUIRED_DONORS=") {
            required_donors = Some(line.between("=", " ").force_usize());
        } else if line.contains("REQUIRED_TWO_CELL_CLONOTYPES=") {
            required_two_cell_clonotypes = Some(line.between("=", " ").force_usize());
        } else if line.contains("REQUIRED_DATASETS=") {
            required_datasets = Some(line.between("=", " ").force_usize());
        }
    }
    let required_fps = required_fps.unwrap();
    let required_cells = required_cells.unwrap();
    let required_clonotypes = required_clonotypes.unwrap();
    let required_donors = required_donors.unwrap();
    let required_two_cell_clonotypes = required_two_cell_clonotypes.unwrap();
    let required_datasets = required_datasets.unwrap();
    let required_fp_percent = format!(
        "{:.2}%",
        percent_ratio(required_fps, required_two_cell_clonotypes)
    );

    // Keep going.

    const ENCLONE_FORMATTED: &str =
        "<span style=\"color:rgb(120,123,175);font-weight:900\">enclone</span>";
    let pwd = env::current_dir().unwrap();
    let pwd = pwd.to_str().unwrap();
    let mut title = String::new();
    let mut extra_head = String::new();
    {
        let f = BufReader::new(File::open(&in_file).expect(&format!(
            "In directory {}, could not open file \"{}\"",
            pwd, &in_file
        )));
        let mut in_head = false;
        for line in f.lines() {
            let s = line.unwrap();
            if s == "</head>" {
                break;
            }
            if in_head {
                extra_head += &format!("{}\n", s);
            }
            if s.contains("<title>") {
                title = s.between("<title>", "</title>").to_string();
            }
            if s == "<head>" {
                in_head = true;
            }
        }
    }
    let f = open_for_read![&in_file];
    let mut g = open_for_write_new![&out_file];
    fwrite!(g, "{}", html_header(level, &title, &extra_head));
    fwriteln!(
        g,
        "
        <! ––\n
        💩 💩 💩 🔴 🔨 🔨 🔨 🔨 🔨 🔨 🔴 💩 💩 💩\n
        PUT DOWN YOUR HAMMER.
        THIS IS AN AUTO-GENERATED FILE.  PLEASE DO NOT EDIT IT.
        THANK YOU FOR YOUR COOPERATION,\n
        SINCERELY,
        THE BENEVOLENT OVERLORDS\n
        💩 💩 💩 🔴 🔨 🔨 🔨 🔨 🔨 🔨 🔴 💩 💩 💩\n
        ––>"
    );
    let mut in_head = false;
    for line in f.lines() {
        let mut s = line.unwrap();
        if s == "<head>" {
            in_head = true;
        }
        if s == "</head>" {
            in_head = false;
            continue;
        }
        if in_head {
            continue;
        }
        if s.starts_with("<title>") {
        } else if s.starts_with("#include ") {
            let mut f = format!("../{}", s.after("#include "));
            if !up {
                f = format!("{}", s.after("#include "));
            }
            let h = open_for_read![&f];
            let mut started = false;
            let mut count = 0;
            for line in h.lines() {
                count += 1;
                let t = line.unwrap();
                if t == "<body>" {
                    started = true;
                } else if t == "</body>" {
                    break;
                } else if started && !t.contains("enclone_banner") {
                    fwriteln!(g, "{}", t);
                }
            }
            if count == 0 {
                eprintln!("\nThe file {} is empty.\n", f);
                std::process::exit(1);
            }
        } else {
            s = s.replace("#enclone", ENCLONE_FORMATTED);
            s = s.replace("#required_fps", &format!("{}", required_fps));
            s = s.replace("#required_cells", &add_commas(required_cells));
            s = s.replace("#required_clonotypes", &add_commas(required_clonotypes));
            s = s.replace("#required_donors", &format!("{}", required_donors));
            s = s.replace("#required_datasets", &format!("{}", required_datasets));
            s = s.replace(
                "#required_two_cell_clonotypes",
                &add_commas(required_two_cell_clonotypes),
            );
            s = s.replace("#required_fp_percent", &format!("{}", &required_fp_percent));
            fwriteln!(g, "{}", s);
        }
    }
}
