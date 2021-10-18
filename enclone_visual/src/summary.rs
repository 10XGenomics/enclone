// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::gui_structures::*;
use crate::style::{ButtonBoxStyle1, ButtonBoxStyle2};
use crate::*;
use enclone_core::stringulate::*;
use iced::Length::Units;
use iced::{Button, Color, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text};
use itertools::{izip, Itertools};
use messages::Message;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This section contains the packing and unpacking of the summary.  The main reason for packing
// various pieces of information into the summary was to avoid invalidating session files that
// people might have written.  We could have done this by raising the output file version, but
// that would have been messy.
//
// There are actually two versions of the packing.  The first existed for about a week.

// pub fn pack_summary() -> String {
//     let mut reply_summary = SERVER_REPLY_SUMMARY.lock().unwrap()[0].clone();
//     let n = SERVER_REPLY_DATASET_NAMES.lock().unwrap().len();
//     for i in 0..n {
//         reply_summary += &mut format!("$$${}", SERVER_REPLY_DATASET_NAMES.lock().unwrap()[i]);
//     }
//     let n = SERVER_REPLY_METRICS.lock().unwrap().len();
//     for i in 0..n {
//         reply_summary += &mut format!("###{}", SERVER_REPLY_METRICS.lock().unwrap()[i]);
//     }
//     reply_summary
// }

pub fn unpack_summary(s: &str) -> Summary {
    // Handle the case of the format that existed only briefly.  This is flaky because it's
    // conceivable that the string $$$ would appear in the second format.

    if s.contains("$$$") {
        let mut dataset_names = Vec::<String>::new();
        let mut metrics = Vec::<Vec<String>>::new();
        let p1 = s.split("$$$").collect::<Vec<&str>>();
        let mut p1x = Vec::<String>::new();
        for i in 0..p1.len() {
            p1x.push(p1[i].to_string());
        }
        let summary = p1x[0].to_string();
        for i in 1..p1x.len() {
            if p1x[i].contains("###") {
                dataset_names.push(p1x[i].before("###").to_string());
                let p2 = p1x[i].split("###").collect::<Vec<&str>>();
                for j in 1..p2.len() {
                    let mut lines = Vec::<String>::new();
                    for line in p2[j].lines() {
                        lines.push(line.to_string());
                    }
                    metrics.push(lines);
                }
            } else {
                dataset_names.push(p1x[i].to_string());
            }
        }
        let mut all_metric_names = Vec::<String>::new();
        for i in 0..metrics.len() {
            for j in 0..metrics[i].len() {
                let s = parse_csv(&metrics[i][j]);
                let m = format!("{},{}", s[0], s[1]);
                all_metric_names.push(m);
            }
        }
        unique_sort(&mut all_metric_names);
        let nm = all_metric_names.len();
        Summary {
            summary: summary,
            dataset_names: dataset_names,
            metrics: metrics,
            metric_selected: vec![false; nm],
            metrics_condensed: false,
        }

    // Handle the current case.
    } else {
        Summary::unpack(&s)
    }
}

pub fn form_summary_from_server_response() -> Summary {
    let summary = SERVER_REPLY_SUMMARY.lock().unwrap()[0].clone();
    let mut dataset_names = Vec::<String>::new();
    let n = SERVER_REPLY_DATASET_NAMES.lock().unwrap().len();
    for i in 0..n {
        dataset_names.push(SERVER_REPLY_DATASET_NAMES.lock().unwrap()[i].clone());
    }
    let mut metrics = Vec::<Vec<String>>::new();
    let n = SERVER_REPLY_METRICS.lock().unwrap().len();
    for i in 0..n {
        let m = SERVER_REPLY_METRICS.lock().unwrap()[i].clone();
        let mut lines = Vec::<String>::new();
        for line in m.lines() {
            lines.push(line.to_string());
        }
        metrics.push(lines);
    }
    let mut all_metric_names = Vec::<String>::new();
    for i in 0..metrics.len() {
        for j in 0..metrics[i].len() {
            let s = parse_csv(&metrics[i][j]);
            let m = format!("{},{}", s[0], s[1]);
            all_metric_names.push(m);
        }
    }
    unique_sort(&mut all_metric_names);
    let nm = all_metric_names.len();
    Summary {
        summary: summary,
        dataset_names: dataset_names,
        metrics: metrics,
        metric_selected: vec![false; nm],
        metrics_condensed: false,
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Clone, Default)]
pub struct Metric {
    pub name: String,
    pub values: Vec<String>,
}

pub fn get_metrics(metricsx: &Vec<Vec<String>>, nd: usize) -> Vec<Metric> {
    let mut all_metric_names = Vec::<String>::new();
    for i in 0..metricsx.len() {
        for j in 0..metricsx[i].len() {
            let s = parse_csv(&metricsx[i][j]);
            let m = format!("{},{}", s[0], s[1]);
            all_metric_names.push(m);
        }
    }
    unique_sort(&mut all_metric_names);
    let nm = all_metric_names.len();
    let mut values = vec![vec![String::new(); nd]; nm];
    for i in 0..nd {
        for j in 0..metricsx[i].len() {
            let s = parse_csv(&metricsx[i][j]);
            let value = s[2].clone();
            let m = format!("{},{}", s[0], s[1]);
            let p = bin_position(&all_metric_names, &m) as usize;
            values[p][i] = value;
        }
    }
    let mut metrics = vec![Metric::default(); nm];
    for i in 0..nm {
        metrics[i] = Metric {
            name: all_metric_names[i].clone(),
            values: values[i].clone(),
        }
    }
    metrics
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn expand_summary(summary: &str, all: bool, show: &Vec<bool>) -> String {
    let summaryx = unpack_summary(&summary);
    let mut summary = String::new();
    if all {
        summary = summaryx.summary.clone();
    }
    let n = summaryx.metrics.len();
    if n > 0 && n == summaryx.dataset_names.len() {
        let dataset_names = summaryx.dataset_names.clone();
        let nd = dataset_names.len();
        let metricsx = summaryx.metrics.clone();
        let metrics = get_metrics(&metricsx, nd);
        let nm = metrics.len();
        let mut categories = Vec::<String>::new();
        for i in 0..nm {
            categories.push(metrics[i].name.before(",").to_string());
        }
        unique_sort(&mut categories);
        for cat in categories.iter() {
            let catc = format!("{},", cat);
            let mut have_some = false;
            for i in 0..nm {
                if show[i] {
                    if metrics[i].name.starts_with(&catc) {
                        have_some = true;
                    }
                }
            }
            if have_some {
                let upcat = cat.to_ascii_uppercase();
                let mut rows = Vec::<Vec<String>>::new();
                let mut row = vec!["metric".to_string()];
                row.append(&mut dataset_names.clone());
                rows.push(row);
                for i in 0..nm {
                    if show[i] {
                        if metrics[i].name.starts_with(&catc) {
                            let mut row = vec![metrics[i].name.after(&catc).to_string()];
                            for j in 0..nd {
                                row.push(metrics[i].values[j].clone());
                            }
                            rows.push(vec!["\\hline".to_string(); nd + 1]);
                            rows.push(row);
                        }
                    }
                }
                let mut log = String::new();
                let mut just = vec![b'l'];
                for _ in 0..nd {
                    just.push(b'|');
                    just.push(b'r');
                }
                print_tabular_vbox(&mut log, &rows, 0, &just, false, false);
                summary += &mut format!("\n{} METRICS BY DATASET\n", upcat);
                summary += &mut log;
            }
        }
    }
    summary
}

pub fn expand_summary_as_csv(summary: &str, show: &Vec<bool>) -> String {
    let summaryx = unpack_summary(&summary);
    let mut summary = String::new();
    let n = summaryx.metrics.len();
    let mut first = true;
    if n > 0 && n == summaryx.dataset_names.len() {
        let dataset_names = summaryx.dataset_names.clone();
        let nd = dataset_names.len();
        let metricsx = summaryx.metrics.clone();
        let metrics = get_metrics(&metricsx, nd);
        let nm = metrics.len();
        let mut categories = Vec::<String>::new();
        for i in 0..nm {
            categories.push(metrics[i].name.before(",").to_string());
        }
        unique_sort(&mut categories);
        for cat in categories.iter() {
            let catc = format!("{},", cat);
            let mut have_some = false;
            for i in 0..nm {
                if show[i] {
                    if metrics[i].name.starts_with(&catc) {
                        have_some = true;
                    }
                }
            }
            if have_some {
                let mut rows = Vec::<Vec<String>>::new();
                if first {
                    let mut row = vec!["class".to_string(), "metric".to_string()];
                    row.append(&mut dataset_names.clone());
                    rows.push(row);
                    first = false;
                }
                for i in 0..nm {
                    if show[i] {
                        if metrics[i].name.starts_with(&catc) {
                            let mut row =
                                vec![cat.clone(), metrics[i].name.after(&catc).to_string()];
                            for j in 0..nd {
                                let m = metrics[i].values[j].replace(",", "");
                                row.push(m);
                            }
                            rows.push(row);
                        }
                    }
                }
                for i in 0..rows.len() {
                    summary += &mut format!("{}\n", rows[i].iter().format("\t "));
                }
            }
        }
    }
    summary
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn max_line_val(s: &str) -> usize {
    let mut m = 0;
    for line in s.lines() {
        let mut nchars = 0;
        for _ in line.chars() {
            nchars += 1;
        }
        m = std::cmp::max(m, nchars);
    }
    m
}

pub fn appropriate_font_size(s: &str, w: u32, max_size: usize) -> usize {
    let mut font_size = 20;
    let max_line = max_line_val(&s);
    const FUDGE: f32 = 175.0;
    let width = (max_line * font_size) as f32 * DEJAVU_WIDTH_OVER_HEIGHT + FUDGE;
    let iwidth = width.ceil() as u32;
    if iwidth > w {
        let fs = w as f32 / width * (font_size as f32);
        font_size = fs.floor() as usize;
    }
    if font_size > max_size {
        font_size = max_size;
    }
    font_size
}

pub fn summary(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let width = (slf.width - 65) as u16; // for text, so scrollbar is not on top of text
    let summary_title = Text::new(&format!("Summary")).size(30);

    // Expand summary.

    let summaryx = unpack_summary(&slf.summary_value);
    let summary = summaryx.summary.clone();
    let n = summaryx.metrics.len();

    // Unpack the summary.  For now we are assuming that it consists of one or two objects,
    // because that's all we have implemented.

    let hets = unpack_to_het_string(&summary);

    // Determine initial font size.

    let mut sum =
        "The summary is empty.  Usually this happens if the command failed.\n".to_string();
    if !hets.is_empty() {
        sum = hets[0].content.clone();
    }
    let mut max_line = max_line_val(&sum);
    let mut font_size = appropriate_font_size(&sum, slf.width, 20);
    const FUDGE: f32 = 175.0;
    let orig_font_size = font_size;

    // Put first part of summary into a scrollable.

    let summary = format!("{}\n \n", sum);
    let mut summary_scrollable = Scrollable::new(&mut slf.summary_scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        .push(
            Text::new(&format!("{}", summary))
                .font(DEJAVU_BOLD)
                .size(orig_font_size as u16),
        );

    // If there is a DescriptionTable, add that.

    let mut descrips = None;
    for j in 1..hets.len() {
        if hets[j].name == "DescriptionTable" {
            descrips = Some(j);
        }
    }
    if descrips.is_some() {
        let j = descrips.unwrap();
        let table = DescriptionTable::from_string(&hets[j].content);
        slf.descrips_for_spreadsheet = table.spreadsheet_text.clone();
        let tables_font_size = appropriate_font_size(&table.display_text, slf.width, 16);
        summary_scrollable = summary_scrollable
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)))
            .push(
                Text::new("Dataset descriptions")
                    .size(25)
                    .color(Color::from_rgb(0.9, 0.0, 0.9)),
            )
            .push(Space::with_height(Units(8)))
            .push(
                Button::new(
                    &mut slf.descrips_copy_button,
                    Text::new("Copy").color(slf.descrips_copy_button_color),
                )
                .on_press(Message::CopyDescrips),
            )
            .push(Space::with_height(Units(8)))
            .push(
                Text::new(&format!("{}", table.display_text))
                    .font(DEJAVU_BOLD)
                    .size(tables_font_size as u16),
            );
    }

    // If there is a FeatureBarcodeAlluvialReadsTableSet, add that.

    let mut alluv = None;
    for j in 1..hets.len() {
        if hets[j].name == "FeatureBarcodeAlluvialReadsTableSet" {
            alluv = Some(j);
        }
    }
    if alluv.is_some() {
        let j = alluv.unwrap();
        let tables = FeatureBarcodeAlluvialReadsTableSet::from_string(&hets[j].content);
        slf.alluvial_reads_tables_for_spreadsheet.clear();
        let mut tables_text = String::new();
        for i in 0..tables.s.len() {
            tables_text += &mut format!(
                "\nfeature barcode read distribution for {}\n{}",
                tables.s[i].id, tables.s[i].display_text
            );
            slf.alluvial_reads_tables_for_spreadsheet += &mut tables.s[i].spreadsheet_text.clone();
        }
        tables_text += "\n \n";
        let tables_font_size = appropriate_font_size(&tables_text, slf.width, 16);
        summary_scrollable = summary_scrollable
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)))
            .push(
                Text::new("Feature barcode read count alluvial tables")
                    .size(25)
                    .color(Color::from_rgb(0.9, 0.0, 0.9)),
            )
            .push(Space::with_height(Units(8)));
        if slf.alluvial_reads_doc_open {
            summary_scrollable = summary_scrollable
                .push(
                    Button::new(
                        &mut slf.close_alluvial_reads_doc_button,
                        Text::new("Hide documentation"),
                    )
                    .on_press(Message::CloseAlluvialReadsDoc),
                )
                .push(Space::with_height(Units(8)))
                .push(
                    Text::new(
                        "For each dataset, we show a table that classifies its feature barcode \
                     reads.\n\n\
                     \
                     These reads are classified as cellular, if their cell barcode was \
                     identified as a cell by the Cell Ranger VDJ pipeline, else noncellular.\n\n\
                     \
                     Each of these categories is subdivided into:\n\
                     • degenerate: R2 starts with at least ten Gs\n\
                     • reference: nondegenerate + feature barcode is in the reference\n\
                     • nonreference: otherwise.\n\n\
                     \
                     In the noncellular degenerate category, canonical reads are \
                     those whose R1 contains CACATCTCCGAGCCCACGAGAC.  This is the end of the \
                     Illumina Nextera version of the R2 primer = \
                     CTGTCTCTTATACACATCTCCGAGCCCACGAGAC.  \
                     If R1 contains the first ten bases = CACATCTCCG, and the read is not \
                     canonical, we call it semicanonical.\n\n\
                     \
                     The number of barcodes shown can be controlled using the extra argument\n\
                     FB_SHOW=k\n\
                     where k is the maximum number of feature barcodes shown for reference (and \
                     likewise for nonreference).  The default value for k is 3.",
                    )
                    .width(Units(width)),
                )
                .push(Space::with_height(Units(8)))
                .push(
                    Text::new(
                        "All the tables can be copied at once, for inclusion in \
                     a spreadsheet, by pushing the button below.  This copies the numbers in the \
                     last column, but not the numbers in the earlier columns.",
                    )
                    .width(Units(width)),
                );
        } else {
            summary_scrollable = summary_scrollable.push(
                Button::new(
                    &mut slf.open_alluvial_reads_doc_button,
                    Text::new("Expand documentation"),
                )
                .on_press(Message::OpenAlluvialReadsDoc),
            );
        }
        summary_scrollable = summary_scrollable
            .push(Space::with_height(Units(8)))
            .push(
                Button::new(
                    &mut slf.alluvial_reads_tables_copy_button,
                    Text::new("Copy").color(slf.alluvial_reads_tables_copy_button_color),
                )
                .on_press(Message::CopyAlluvialReadsTables),
            )
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)))
            .push(
                Text::new(&format!("{}", tables_text))
                    .font(DEJAVU_BOLD)
                    .size(tables_font_size as u16),
            );
    }

    // If there is a FeatureBarcodeAlluvialTableSet, add that.

    let mut alluv = None;
    for j in 1..hets.len() {
        if hets[j].name == "FeatureBarcodeAlluvialTableSet" {
            alluv = Some(j);
        }
    }
    if alluv.is_some() {
        let j = alluv.unwrap();
        let tables = FeatureBarcodeAlluvialTableSet::from_string(&hets[j].content);
        slf.alluvial_tables_for_spreadsheet.clear();
        let mut tables_text = String::new();
        for i in 0..tables.s.len() {
            tables_text += &mut format!(
                "\nfeature barcode UMI distribution for {}\n{}",
                tables.s[i].id, tables.s[i].display_text
            );
            slf.alluvial_tables_for_spreadsheet += &mut tables.s[i].spreadsheet_text.clone();
        }
        tables_text += "\n \n";
        let tables_font_size = appropriate_font_size(&tables_text, slf.width, 20);
        summary_scrollable = summary_scrollable
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)))
            .push(
                Text::new("Feature barcode UMI count alluvial tables")
                    .size(25)
                    .color(Color::from_rgb(0.9, 0.0, 0.9)),
            )
            .push(Space::with_height(Units(8)))
            .push(Text::new(
                "These tables are similar to the tables for reads.  See the documentation there.  \
                 UMIs for degenerate reads are not included in these tables.  \
                 Note that the FB_SHOW option can also be used here.\n\n\
                 All the tables can be copied at once, in a form suitable for inclusion in \
                 a spreadsheet, by pushing the button below.  This copies the numbers in the \
                 last column, but not the numbers in the earlier columns.",
            ))
            .push(Space::with_height(Units(8)))
            .push(
                Button::new(
                    &mut slf.alluvial_tables_copy_button,
                    Text::new("Copy").color(slf.alluvial_tables_copy_button_color),
                )
                .on_press(Message::CopyAlluvialTables),
            )
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)))
            .push(
                Text::new(&format!("{}", tables_text))
                    .font(DEJAVU_BOLD)
                    .size(tables_font_size as u16),
            );
    }

    // Suppose we have dataset level metrics.

    let mut button_text_row = Row::new();
    if n > 0 && n == summaryx.dataset_names.len() {
        let dataset_names = summaryx.dataset_names.clone();
        let nd = dataset_names.len();
        let metricsx = summaryx.metrics.clone();
        let metrics = get_metrics(&metricsx, nd);
        let nm = metrics.len();
        if slf.metric_button.len() != nm {
            slf.metric_button = vec![iced::button::State::default(); nm];
        }
        let mut categories = Vec::<String>::new();
        for i in 0..nm {
            categories.push(metrics[i].name.before(",").to_string());
        }
        unique_sort(&mut categories);

        // Make text for metrics.

        let mut text = String::new();
        for cat in categories.iter() {
            let mut have_some = false;
            let catc = format!("{},", cat);
            if !slf.metrics_condensed {
                have_some = true;
            } else {
                for i in 0..nm {
                    if metrics[i].name.starts_with(&catc) {
                        if slf.metric_selected[i] {
                            have_some = true;
                        }
                    }
                }
            }
            if have_some {
                let upcat = cat.to_ascii_uppercase();
                let mut rows = Vec::<Vec<String>>::new();
                let mut row = vec!["metric".to_string()];
                row.append(&mut dataset_names.clone());
                rows.push(row);
                for i in 0..nm {
                    if metrics[i].name.starts_with(&catc) {
                        if slf.metrics_condensed && !slf.metric_selected[i] {
                            continue;
                        }
                        let mut row = vec![metrics[i].name.after(&catc).to_string()];
                        for j in 0..nd {
                            row.push(metrics[i].values[j].clone());
                        }
                        rows.push(vec!["\\hline".to_string(); nd + 1]);
                        rows.push(row);
                    }
                }
                let mut log = String::new();
                let mut just = vec![b'l'];
                for _ in 0..nd {
                    just.push(b'|');
                    just.push(b'r');
                }
                print_tabular_vbox(&mut log, &rows, 0, &just, false, false);
                text += &mut format!("\n{} METRICS BY DATASET\n", upcat);
                text += &mut log;
            }
        }
        text = format!("{}\n \n", text);

        // Update font size.

        for line in text.lines() {
            let mut nchars = 0;
            for _ in line.chars() {
                nchars += 1;
            }
            max_line = std::cmp::max(max_line, nchars);
        }
        let width = (max_line * font_size) as f32 * DEJAVU_WIDTH_OVER_HEIGHT + FUDGE;
        let iwidth = width.ceil() as u32;
        if iwidth > slf.width {
            let fs = slf.width as f32 / width * (font_size as f32);
            font_size = fs.floor() as usize;
        }

        // Make text column for metrics.

        let font_size = font_size as u16;
        let text_column = Column::new().push(Text::new(&text).font(DEJAVU_BOLD).size(font_size));

        // Make button column for metrics.

        let mut button_column = Column::new();
        for (i, y) in izip!(0..nm, slf.metric_button.iter_mut()) {
            if i == 0 || metrics[i].name.before(",") != metrics[i - 1].name.before(",") {
                button_column = button_column.push(Space::with_height(Units(4 * font_size)));
            }
            if i > 0 && metrics[i].name.before(",") != metrics[i - 1].name.before(",") {
                button_column = button_column.push(Space::with_height(Units(font_size)));
            }
            let mut button = Button::new(
                y,
                Text::new("")
                    .height(Units(font_size))
                    .width(Units(font_size)),
            );
            if !slf.metric_selected[i] {
                button = button.style(ButtonBoxStyle1);
            } else {
                button = button.style(ButtonBoxStyle2);
            }
            button = button.padding(0).on_press(Message::MetricButton(i));
            button_column = button_column
                .push(Space::with_height(Units(font_size)))
                .push(button);
        }

        // Put together buttons and text.

        if slf.metrics_condensed {
            button_text_row = Row::new().push(text_column);
        } else {
            button_text_row = Row::new()
                .push(button_column)
                .push(Space::with_width(Units(font_size / 4)))
                .push(text_column);
        }
    }

    // Build final structure.

    let summary_snapshot_button = Button::new(
        &mut slf.graphic_snapshot_button,
        Text::new("Snapshot").color(slf.summary_snapshot_button_color),
    )
    .on_press(Message::SummarySnapshot);
    let summary_copy_button = Button::new(
        &mut slf.summary_copy_button,
        Text::new("Copy").color(slf.copy_summary_button_color),
    )
    .on_press(Message::CopySummary);
    let summary_close_button = Button::new(&mut slf.summary_button, Text::new("Dismiss"))
        .on_press(Message::SummaryClose(Ok(())));
    let top_bar = Row::new()
        .push(summary_title)
        .push(Space::with_width(Length::Fill))
        .push(summary_snapshot_button)
        .push(Space::with_width(Units(8)))
        .push(summary_copy_button)
        .push(Space::with_width(Units(8)))
        .push(summary_close_button);
    let mut have_metrics = false;
    for m in summaryx.metrics.iter() {
        if m.len() > 0 {
            have_metrics = true;
        }
    }
    if have_metrics && n == summaryx.dataset_names.len() {
        summary_scrollable = summary_scrollable
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)));
        if !slf.metrics_condensed {
            summary_scrollable = summary_scrollable
                .push(
                    Text::new("Dataset level metrics")
                        .size(25)
                        .color(Color::from_rgb(0.9, 0.0, 0.9)),
                )
                .push(Space::with_height(Units(8)))
                .push(Text::new(
                    "Metrics below can be selectively displayed by clicking on boxes, \
                    and then pushing the button below.",
                ))
                .push(Space::with_height(Units(4)))
                .push(Text::new(
                    "The display choices made here are \
                    saveable, but cannot be recapitulated using an enclone command.",
                ))
                .push(Space::with_height(Units(4)))
                .push(Text::new(
                    "The copy selected metrics button may be used to copy the selected \
                    metrics to the clipboard, in a form that can be pasted into a spreadsheet.",
                ))
                .push(Space::with_height(Units(12)));
        }
        let text = if slf.metrics_condensed {
            "Show all metrics".to_string()
        } else {
            "Show selected metrics".to_string()
        };
        summary_scrollable = summary_scrollable.push(
            Button::new(&mut slf.condense_metrics_button, Text::new(&text))
                .on_press(Message::CondenseMetrics),
        );
        summary_scrollable = summary_scrollable.push(Space::with_height(Units(8)));
        summary_scrollable = summary_scrollable.push(
            Button::new(
                &mut slf.copy_selected_metrics_button,
                Text::new("Copy selected metrics").color(slf.copy_selected_metrics_button_color),
            )
            .on_press(Message::CopySelectedMetrics),
        );
        summary_scrollable = summary_scrollable
            .push(Space::with_height(Units(8)))
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(button_text_row);
    }
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(summary_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
