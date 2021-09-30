// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::gui_structures::*;
use crate::style::{ButtonBoxStyle1, ButtonBoxStyle2};
use crate::*;
use iced::Length::Units;
use iced::{Button, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text};
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

pub fn summary(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let summary_title = Text::new(&format!("Summary")).size(30);

    // Expand summary.

    let summaryx = unpack_summary(&slf.summary_value);
    let summary = summaryx.summary.clone();
    let n = summaryx.metrics.len();

    // Determine initial font size.

    let mut font_size = 20;
    let mut max_line = 0;
    for line in summary.lines() {
        let mut nchars = 0;
        for _ in line.chars() {
            nchars += 1;
        }
        max_line = std::cmp::max(max_line, nchars);
    }
    const FUDGE: f32 = 175.0;
    let width = (max_line * font_size) as f32 * DEJAVU_WIDTH_OVER_HEIGHT + FUDGE;
    let iwidth = width.ceil() as u32;
    if iwidth > slf.width {
        let fs = slf.width as f32 / width * (font_size as f32);
        font_size = fs.floor() as usize;
    }
    let orig_font_size = font_size;

    // Suppose we have dataset level metrics.

    let mut button_text_row = Row::new();
    if n > 0 && n == summaryx.dataset_names.len() {
        let dataset_names = summaryx.dataset_names.clone();
        let nd = dataset_names.len();
        let metricsx = summaryx.metrics.clone();
        let metrics = get_metrics(&metricsx, nd);
        let nm = metrics.len();
        if slf.metric_button.is_empty() {
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

    let summary = format!("{}\n \n", summary);
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
        .push(summary_copy_button)
        .push(Space::with_width(Units(8)))
        .push(summary_close_button);
    let mut summary_scrollable = Scrollable::new(&mut slf.scroll)
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
    if n > 0 && n == summaryx.dataset_names.len() {
        summary_scrollable = summary_scrollable
            .push(Rule::horizontal(10).style(style::RuleStyle2))
            .push(Space::with_height(Units(8)));
        if !slf.metrics_condensed {
            summary_scrollable = summary_scrollable
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
                .push(Space::with_height(Units(8)));
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
