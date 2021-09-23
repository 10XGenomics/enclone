// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::Length::Units;
use iced::{Button, Column, Container, Element, Image, Length, Row, Rule, Scrollable, Space, Text};
use messages::Message;

pub fn graphic(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let graphic_title = Text::new(&format!("Graphic")).size(30);
    let graphic_snapshot_button = Button::new(
        &mut slf.graphic_snapshot_button,
        Text::new("Snapshot").color(slf.graphic_snapshot_button_color),
    )
    .on_press(Message::GraphicSnapshot);
    let graphic_close_button = Button::new(&mut slf.graphic_close_button, Text::new("Dismiss"))
        .on_press(Message::GraphicClose);
    let top_bar = Row::new()
        .push(graphic_title)
        .push(Space::with_width(Length::Fill))
        .push(graphic_snapshot_button)
        .push(Space::with_width(Units(8)))
        .push(graphic_close_button);
    let svg_height = CURRENT_HEIGHT.load(SeqCst) as u16;
    let have_canvas = slf.canvas_view.state.geometry_value.is_some();
    let mut graphic_row = Row::new().spacing(10);
    if slf.svg_value.len() > 0 {
        if have_canvas {
            graphic_row = graphic_row
                .push(
                    slf.canvas_view
                        .view()
                        .map(move |message| Message::GroupClicked(message)),
                )
                .height(Units(svg_height));
        } else {
            let svg_as_png = Image::new(iced::image::Handle::from_memory(slf.png_value.clone()))
                .height(Units(svg_height));
            graphic_row = graphic_row.push(svg_as_png);
        }
    }
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(graphic_row);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

pub fn clonotypes(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let clonotypes_title = Text::new(&format!("Clonotypes")).size(30);

    let copy_button = Button::new(
        &mut slf.clonotypes_copy_button,
        Text::new("Copy text").color(slf.clonotypes_copy_button_color),
    )
    .on_press(Message::ClonotypesCopy);

    let clonotypes_close_button =
        Button::new(&mut slf.clonotypes_close_button, Text::new("Dismiss"))
            .on_press(Message::ClonotypesClose);
    const CLONOTYPE_FONT_SIZE: u16 = 13;
    let font_width = CLONOTYPE_FONT_SIZE as f32 * DEJAVU_WIDTH_OVER_HEIGHT;
    let available = slf.width - (3 * SPACING + SCROLLBAR_WIDTH) as u32;
    let nchars = (available as f32 / font_width).round() as usize;
    let mut trunc = String::new();
    let failed = slf.output_value.contains("enclone failed");
    for line in slf.output_value.lines() {
        for (i, c) in line.chars().enumerate() {
            if i == nchars && !failed {
                break;
            }
            trunc.push(c);
        }
        trunc.push('\n');
    }
    let top_bar = Row::new()
        .push(clonotypes_title)
        .push(Space::with_width(Length::Fill))
        .push(copy_button)
        .push(Space::with_width(Units(8)))
        .push(clonotypes_close_button);
    let clonotypes_scrollable = Scrollable::new(&mut slf.clonotypes_scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        .push(
            Text::new(&trunc)
                .font(DEJAVU_BOLD)
                .size(CLONOTYPE_FONT_SIZE),
        );
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(clonotypes_scrollable);
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
    let intro_text = "This console is mostly for developers.  It shows some logging and provides \
        some tools for tracing bugs and maintaining cookbooks.";
    let sanity_text = "The sanity check button is for debugging.  It causes some sanity checks \
        to be run on the current enclone visual state.";
    let sanity_button = Button::new(
        &mut slf.sanity_button,
        Text::new("Sanity check").color(slf.sanity_button_color),
    )
    .on_press(Message::SanityCheck);
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
        .push(Text::new(intro_text))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Text::new(sanity_text))
        .push(sanity_button)
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
