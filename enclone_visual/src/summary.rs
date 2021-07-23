// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::{Button, Column, Container, Element, Length, Scrollable, Text};
use messages::Message;

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
    let summary_close_button =
        Button::new(&mut slf.open_state, Text::new("Vanish!")).on_press(Message::SummaryClose);
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
