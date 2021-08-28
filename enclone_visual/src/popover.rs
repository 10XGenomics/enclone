// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::Length::Units;
use iced::{Button, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text};
use itertools::Itertools;
use messages::Message;
use vector_utils::*;

pub fn console(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let console_title = Text::new(&format!("Console")).size(30);
    let mut console = String::new();
    let n = CONSOLE.lock().unwrap().len();
    for i in 0..n {
        console += &mut format!("{}\n", CONSOLE.lock().unwrap()[i]);
    }
    console += " \n\n\n";
    let console_close_button = Button::new(&mut slf.console_close_button, Text::new("Dismiss"))
        .on_press(Message::ConsoleClose);
    let top_bar = Row::new()
        .push(console_title)
        .push(Space::with_width(Length::Fill))
        .push(console_close_button);
    let font_size = 15;
    let console_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        .push(
            Text::new(&format!("{}", console))
                .font(DEJAVU_BOLD)
                .size(font_size as u16),
        );
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(console_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

pub fn summary(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let summary_title = Text::new(&format!("Summary")).size(30);
    let mut summary = SUMMARY_CONTENTS.lock().unwrap()[0].clone();

    // Disect the summary, unscrewing the grotesque way that information is packed in it.

    let mut dataset_names = Vec::<String>::new();
    let mut metrics = Vec::<Vec<String>>::new();
    if summary.contains("$$$") {
        let p1 = summary.split("$$$").collect::<Vec<&str>>();
        let mut p1x = Vec::<String>::new();
        for i in 0..p1.len() {
            p1x.push(p1[i].to_string());
        }
        summary = p1x[0].to_string();
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
        if dataset_names.len() == metrics.len() {
            let mut all_metric_names = Vec::<String>::new();
            for i in 0..metrics.len() {
                for j in 1..metrics[i].len() {
                    let mut s = parse_csv(&metrics[i][j]);
                    s.truncate(s.len() - 1);
                    let m = format!("{}", s.iter().format(","));
                    all_metric_names.push(m);
                }
            }
            unique_sort(&mut all_metric_names);
            let nd = dataset_names.len();
            let nm = all_metric_names.len();
            let mut values = vec![vec![String::new(); nm]; nd];
            for i in 0..nd {
                for j in 1..metrics[i].len() {
                    let mut s = parse_csv(&metrics[i][j]);
                    let value = s[s.len() - 1].clone();
                    s.truncate(s.len() - 1);
                    let m = format!("{}", s.iter().format(","));
                    let p = bin_position(&all_metric_names, &m) as usize;
                    values[i][p] = value;
                }
            }
            let mut rows = Vec::<Vec<String>>::new();
            let mut row = vec!["metric".to_string()];
            row.append(&mut dataset_names.clone());
            rows.push(row);
            for i in 0..nm {
                let mut row = vec![all_metric_names[i].clone()];
                for j in 0..nd {
                    row.push(values[j][i].clone());
                }
                rows.push(vec!["\\hline".to_string(); nd + 1]);
                rows.push(row);
            }
            let mut log = String::new();
            let mut just = vec![b'l'];
            for _ in 0..nd {
                just.push(b'|');
                just.push(b'r');
            }
            print_tabular_vbox(&mut log, &rows, 0, &just, false, false);
            summary += "\n\n";
            summary += &mut log;
        }
    }

    // Proceed.

    let summary = format!("{}\n \n", summary);
    let font_size = 20;
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

pub fn cookbook(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let cookbook_title = Text::new(&format!("Cookbook")).size(30);
    let preamble = "Type the tag into the input box to run the given command.\n\n";
    let cookbook_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        .push(
            Text::new(&format!(
                "{}{}",
                preamble,
                COOKBOOK_CONTENTS.lock().unwrap()[0]
            ))
            .font(DEJAVU_BOLD)
            .size(14),
        );
    let cookbook_close_button =
        Button::new(&mut slf.open_state, Text::new("Vanish!")).on_press(Message::CookbookClose);
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(cookbook_title)
        .push(cookbook_scrollable)
        .push(cookbook_close_button);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
