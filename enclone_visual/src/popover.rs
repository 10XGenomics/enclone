// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::Length::Units;
use iced::{
    Button, Checkbox, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text,
};
use io_utils::*;
use messages::Message;
use std::env;
use std::fs::metadata;

pub fn archive(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let archive_title = Text::new(&format!("Archive")).size(30);
    let archive_close_button = Button::new(&mut slf.archive_close_button, Text::new("Dismiss"))
        .on_press(Message::ArchiveClose);
    let top_bar = Row::new()
        .push(archive_title)
        .push(Space::with_width(Length::Fill))
        .push(archive_close_button);
    let text1 = Text::new(
        "For a given enclone visual session, if you click the Save On Exit \
        box (which will make the text red), then when you later push the Exit button, your session \
        will be saved.  Pushing repeatedly toggles the state.",
    );
    let text2 = Text::new(
        "You can restore a previously saved session by clicking on one of the \
        boxes below.",
    );
    let mut hist_dir;
    let mut hist = Vec::<String>::new();
    for (key, value) in env::vars() {
        if key == "HOME" {
            let home = value.clone();
            let enclone = format!("{}/enclone", home);
            hist_dir = format!("{}/visual_history", enclone);
            if path_exists(&hist_dir) {
                if metadata(&hist_dir).unwrap().is_dir() {
                    hist = dir_list(&hist_dir);
                }
            }
        }
    }
    hist.reverse();
    let mut archive_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle);
    for (i, x) in hist.iter().enumerate() {
        let row = Row::new()
            .push(
                Text::new(&format!(
                    "{:<3} {}    {}    ",
                    i + 1,
                    x.before("___"),
                    x.after("___")
                ))
                .font(DEJAVU),
            )
            .push(Checkbox::new(slf.enabled, "", Message::Restore));
        if i > 0 {
            archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
        }
        archive_scrollable = archive_scrollable.push(row);
    }
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(text1)
        .push(text2)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(archive_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

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
    let summary = SUMMARY_CONTENTS.lock().unwrap()[0].clone();
    let nlines = summary.chars().filter(|&n| n == '\n').count();
    let font_size = (15 * nlines) / 38;
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
    let summary_close_button = Button::new(&mut slf.open_state, Text::new("Vanish!"))
        .on_press(Message::SummaryClose(Ok(())));
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(summary_title)
        .push(summary_scrollable)
        .push(summary_close_button);
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
