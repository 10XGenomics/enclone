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
        Text::new("enclone visual can save sessions to the directory ~/enclone/visual/history.");
    let text1 = Text::new(
        "For a given enclone visual session:\n\
         • click the Save On Exit box on the main screen; it will turn red \
           (pushing again toggles state)\n\
         • when you later push the Exit button, your session will be saved.",
    );
    let text2 =
        Text::new("▒ expand - Display the commands in a saved session by checking the expand box.");
    let text3 = Text::new(
        "▒ restore - Restore a previously saved session by checking the restore box.  \
            This deletes your current session!",
    );
    let text4 =
        Text::new("▒ delete - Delete a previously saved session by checking the delete box.");
    let text5 = Text::new(
        "▒ share - Share with other users by checking the share box.  You will be prompted for \
            their names.\n\
         Conversely another user may share with you, and if so you will be given the option \
            to accept the share, below here.",
    );
    let text6 = Text::new(
        "▒ name - Name of this or a previous session is displayed:\n\
         • change it using the rectangular box and then check to its right\n\
         • unchecking restores the previous name\n\
         • long names are allowed but incompletely displayed.",
    );
    let labels =
        Text::new("#   date        time     expand     restore    delete    share    name")
            .font(DEJAVU);

    fn share_col() -> Column<'static, Message> {
        let c = Color::from_rgb(0.4, 0.1, 0.2);
        Column::new()
        .push(Space::with_height(Units(8)))
        .push(
            Text::new(
                "           Enter user names (usually first.last), one per line, \
            and check the box to the right.",
            )
            .size(16)
            .color(c),
        )
        .push(Space::with_height(Units(4)))
        .push(
            Text::new("           A new line will appear once you do so.")
                .size(16)
                .color(c),
        )
        .push(Space::with_height(Units(4)))
        .push(
            Text::new(
                "           Lines may be prepopulated based on your recent shares; \
            you can check them.",
            )
            .size(16)
            .color(c),
        )
        .push(Space::with_height(Units(4)))
        .push(
            Text::new("           Push Dismiss at the top when you're done, or uncheck share to cancel.")
                .size(16)
                .color(c),
        )
    }

    let mut archive_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle);
    let row = Row::new()
        .align_items(Align::Center)
        .push(Text::new("0   today       now         ").font(DEJAVU))
        .push(Space::with_width(Units(301)))
        .push(Checkbox::new(slf.share_requested, "", Message::Share))
        .push(Space::with_width(Units(58)))
        .push(TextInput::new(&mut slf.name, "", &slf.h.name_value, Message::Name).padding(2))
        .push(Space::with_width(Units(8)))
        .push(Checkbox::new(
            slf.name_change_requested,
            "",
            Message::NameChange,
        ));
    archive_scrollable = archive_scrollable.push(row);
    if slf.share_requested {
        archive_scrollable = archive_scrollable.push(share_col());
    }
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
            let mut row = Row::new().align_items(Align::Center);
            let date = x.before("___").to_string();
            let date1 = stringme(&date.as_bytes()[2..4]);
            let date2 = date.after("-").to_string();
            let mut date = format!("{}-{}", date2, date1);
            let mut time = x.after("___").rev_before(":").to_string();
            if TEST_MODE.load(SeqCst) {
                date = stringme(&vec![b' '; 8]);
                time = stringme(&vec![b' '; 5]);
            }
            row = row.push(
                Text::new(&format!("{:<3} {}    {}    ", count + 1, date, time,)).font(DEJAVU),
            );
            row = row
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
                .push(Space::with_width(Units(79)))
                .push(Checkbox::new(
                    slf.delete_requested[i],
                    "",
                    move |x: bool| Message::DeleteArchiveEntry(x, i),
                ))
                .push(Space::with_width(Units(68)))
                .push(Checkbox::new(
                    slf.archive_share_requested[i],
                    "",
                    move |x: bool| Message::ArchiveShare(x, i),
                ));
            row = row.push(Space::with_width(Units(58)));
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
            if i > 0 {
                archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
            }
            archive_scrollable = archive_scrollable.push(row);
            if slf.restore_msg[i].len() > 0 {
                archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
                let mut row = Row::new();
                row = row.push(Space::with_width(Units(43)));
                row = row.push(
                    Text::new(&format!("{}", slf.restore_msg[i]))
                        .font(DEJAVU)
                        .color(Color::from_rgb(1.0, 0.0, 0.0))
                        .size(16),
                );
                archive_scrollable = archive_scrollable.push(row);
            }

            if slf.archive_share_requested[i] {
                archive_scrollable = archive_scrollable.push(share_col());
            }

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
        .push(text6)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(labels)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(archive_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
