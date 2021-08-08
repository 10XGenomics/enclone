// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Be very careful about editing the archive function as it is delicate and under-tested.

use crate::*;
use iced::Length::Units;
use iced::{
    Align, Button, Checkbox, Color, Column, Container, Element, Length, Row, Rule, Scrollable,
    Space, Text, TextInput,
};
use io_utils::*;
use messages::Message;

pub fn archive(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let archive_title = Text::new(&format!("Archive")).size(30);
    let archive_close_button = Button::new(&mut slf.archive_close_button, Text::new("Dismiss"))
        .on_press(Message::ArchiveClose);
    let top_bar = Row::new()
        .push(archive_title)
        .push(Space::with_width(Length::Fill))
        .push(archive_close_button);
    let text0 =
        Text::new("enclone visual can save sessions to the directory ~/enclone/visual_history.");
    let text1 = Text::new(
        "For a given enclone visual session:\n\
         • click the Save On Exit box on the main screen; it will turn red \
           (pushing again toggles state)\n\
         • when you later push the Exit button, your session will be saved.",
    );
    let text2 = Text::new("Display the commands in a saved session by clicking on the expand box.");
    let text3 = Text::new(
        "Restore a previously saved session by clicking on the restore box.  \
            This will delete your current session!",
    );
    let text4 = Text::new("Delete a previously saved session by clicking on the delete box.");
    let text5 = Text::new(
        "Name of this or a previous session is displayed:\n\
         • change it using the rectangular box and then check to its right\n\
         • unchecking restores the previous name\n\
         • long names are allowed but incompletely displayed.",
    );
    let labels = Text::new("#   date          time        expand     restore    delete    name")
        .font(DEJAVU);
    let mut archive_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle);
    let row = Row::new()
        .align_items(Align::Center)
        .push(Text::new("0   today         now").font(DEJAVU))
        .push(Space::with_width(Units(424)))
        .push(TextInput::new(&mut slf.name, "", &slf.h.name_value, Message::Name).padding(2))
        .push(Space::with_width(Units(8)))
        .push(Checkbox::new(
            slf.name_change_requested,
            "",
            Message::NameChange,
        ));
    archive_scrollable = archive_scrollable.push(row);
    archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
    let mut count = 0;
    for (i, y) in slf.archive_name.iter_mut().enumerate() {
        let x = &slf.archive_list[i];
        let path = format!("{}/{}", slf.archive_dir.as_ref().unwrap(), x);
        let mut make_row = x.contains("___") && !slf.deleted[i];
        if !path_exists(&path) && !slf.delete_requested[i] {
            make_row = false;
        }
        if make_row {
            let mut row = Row::new()
                .align_items(Align::Center)
                .push(
                    Text::new(&format!(
                        "{:<3} {}    {}    ",
                        count + 1,
                        x.before("___"),
                        x.after("___")
                    ))
                    .font(DEJAVU),
                )
                .push(Checkbox::new(
                    slf.expand_archive_entry[i],
                    "",
                    move |x: bool| Message::ExpandArchiveEntry(x, i),
                ))
                .push(Space::with_width(Units(80)))
                .push(Checkbox::new(
                    slf.restore_requested[i],
                    "",
                    move |x: bool| Message::Restore(x, i),
                ))
                .push(Space::with_width(Units(78)))
                .push(Checkbox::new(
                    slf.delete_requested[i],
                    "",
                    move |x: bool| Message::DeleteArchiveEntry(x, i),
                ));
            row = row.push(Space::with_width(Units(68)));
            row = row.push(
                TextInput::new(y, "", &slf.archive_name_value[i], move |x: String| {
                    Message::ArchiveName(x, i)
                })
                .padding(2),
            );
            row = row.push(Space::with_width(Units(8)));
            row = row.push(Checkbox::new(
                slf.archive_name_change_requested[i],
                "",
                move |x: bool| Message::ArchiveNameChange(x, i),
            ));
            /*
            if slf.restore_msg[i].len() > 0 {
                row = row.push(Space::with_width(Units(70)));
                row = row.push(Text::new(&format!("{}", slf.restore_msg[i])).font(DEJAVU));
            }
            */
            if i > 0 {
                archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
            }
            archive_scrollable = archive_scrollable.push(row);
            if slf.expand_archive_entry[i] {
                let clist = &slf.archived_command_list[i].as_ref().unwrap();
                if !clist.is_empty() {
                    archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
                }
                for (j, y) in clist.iter().enumerate() {
                    // Determine max_line.  Sloppy, as we don't actually do the correct math.
                    const MAX_LINE: usize = 110;
                    let max_line = MAX_LINE as f64 * slf.width as f64 / INITIAL_WIDTH as f64;
                    let max_line = max_line.round() as usize;
                    let lines = fold(&y, max_line);
                    for k in 0..lines.len() {
                        let s = if k == 0 {
                            format!("     {} = {}", j + 1, lines[k])
                        } else {
                            let len = format!("     {} = ", j + 1).len();
                            let b = stringme(&vec![b' '; len]);
                            format!("{}{}", b, lines[k])
                        };
                        let row = Row::new().push(
                            Text::new(&s)
                                .font(DEJAVU)
                                .color(Color::from_rgb(0.0, 0.0, 0.4))
                                .size(16),
                        );
                        archive_scrollable = archive_scrollable.push(row);
                        if j < clist.len() - 1 || k < lines.len() - 1 {
                            archive_scrollable =
                                archive_scrollable.push(Space::with_height(Units(4)));
                        }
                    }
                }
                if !clist.is_empty() {
                    archive_scrollable = archive_scrollable.push(Space::with_height(Units(3)));
                }
            }
            count += 1;
        }
    }
    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(text0)
        .push(text1)
        .push(text2)
        .push(text3)
        .push(text4)
        .push(text5)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(labels)
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
