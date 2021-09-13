// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::{Button, Column, Container, Element, Length, Row, Rule, Scrollable, Space, Text};
use messages::Message;

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
    let console_text = "The current console output is shown below.";
    let recompute_text = "The recompute button causes each state within the given session to be \
        recomputed.  The purpose of this is to update the session to reflect the results of the \
        current enclone visual code, in the event that it has changed.  The recompute button will \
        cause the screen to flash as it recapitulates each computation.  When done, it will return \
        to this page.\n\n\
        Typically the way you would use this is to go to archive page, restore a session, then \
        come here and push Recompute, then exit this page.  If the recomputed session looks right, \
        then you would push Save, and finally go to the archive page and delete the original \
        version.";
    let recompute_button =
        Button::new(&mut slf.recompute_button, Text::new("Recompute")).on_press(Message::Recompute);
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
        .push(Text::new(recompute_text))
        .push(recompute_button)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Text::new(console_text))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(console_scrollable);
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
