// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::Length::Units;
use iced::{Button, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text};
use messages::Message;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This section contains the grotesque packing and unpacking of the summary.  We did this to
// avoid changing history (which would require introducing a new history version)

pub fn pack_summary() -> String {
    let mut reply_summary = SERVER_REPLY_SUMMARY.lock().unwrap()[0].clone();
    let n = SERVER_REPLY_DATASET_NAMES.lock().unwrap().len();
    for i in 0..n {
        reply_summary += &mut format!("$$${}", SERVER_REPLY_DATASET_NAMES.lock().unwrap()[i]);
    }
    let n = SERVER_REPLY_METRICS.lock().unwrap().len();
    for i in 0..n {
        reply_summary += &mut format!("###{}", SERVER_REPLY_METRICS.lock().unwrap()[i]);
    }
    reply_summary
}

// #[derive(Default, PartialEq, Clone)]
pub struct SummaryStuff {
    pub summary: String,
    pub dataset_names: Vec<String>,
    pub metrics: Vec<Vec<String>>,
}

impl SummaryStuff {
    pub fn unpack_summary(s: &str) -> Self {
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
        SummaryStuff {
            summary: summary,
            dataset_names: dataset_names,
            metrics: metrics,
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn expand_summary(summary: &str) -> String {
    let mut summary = summary.to_string();
    let summary_stuff = SummaryStuff::unpack_summary(&summary);
    let n = summary_stuff.metrics.len();
    if n > 0 && n == summary_stuff.dataset_names.len() {
        summary = summary_stuff.summary.clone();
        let dataset_names = summary_stuff.dataset_names.clone();
        let metrics = summary_stuff.metrics.clone();
        let mut all_metric_names = Vec::<String>::new();
        for i in 0..metrics.len() {
            for j in 0..metrics[i].len() {
                let s = parse_csv(&metrics[i][j]);
                let m = format!("{},{}", s[0], s[1]);
                all_metric_names.push(m);
            }
        }
        unique_sort(&mut all_metric_names);
        let nd = dataset_names.len();
        let nm = all_metric_names.len();
        let mut values = vec![vec![String::new(); nm]; nd];
        for i in 0..nd {
            for j in 0..metrics[i].len() {
                let s = parse_csv(&metrics[i][j]);
                let value = s[2].clone();
                let m = format!("{},{}", s[0], s[1]);
                let p = bin_position(&all_metric_names, &m) as usize;
                values[i][p] = value;
            }
        }
        let mut categories = Vec::<String>::new();
        for i in 0..nm {
            categories.push(all_metric_names[i].before(",").to_string());
        }
        unique_sort(&mut categories);
        for cat in categories.iter() {
            let catc = format!("{},", cat);
            let upcat = cat.to_ascii_uppercase();
            let mut rows = Vec::<Vec<String>>::new();
            let mut row = vec!["metric".to_string()];
            row.append(&mut dataset_names.clone());
            rows.push(row);
            for i in 0..nm {
                if all_metric_names[i].starts_with(&catc) {
                    let mut row = vec![all_metric_names[i].clone().after(&catc).to_string()];
                    for j in 0..nd {
                        row.push(values[j][i].clone());
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
            summary += &mut format!("\n{} METRICS BY DATASET\n", upcat);
            summary += &mut log;
        }
    }
    summary
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn summary(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let summary_title = Text::new(&format!("Summary")).size(30);
    let summary = expand_summary(&slf.summary_value);
    let summary = format!("{}\n \n", summary);
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
    let summary_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        .push(
            Text::new(&format!("{}", summary))
                .font(DEJAVU_BOLD)
                .size(font_size as u16),
        );
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
