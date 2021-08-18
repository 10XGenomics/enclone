// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Be very careful about editing the archive function as it is delicate and under-tested.

use crate::*;
use iced::Length::Units;
use iced::{
    Align, Button, Checkbox, Color, Column, Container, Element, Length, Row, Rule, Scrollable,
    Space, Text, TextInput,
};
use io_utils::*;
use itertools::izip;
use messages::Message;

pub fn archive(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let archive_title = Text::new(&format!("Archive")).size(30);
    let archive_close_button = Button::new(&mut slf.archive_close_button, Text::new("Dismiss"))
        .on_press(Message::ArchiveClose);
    let open_archive_doc_button = Button::new(
        &mut slf.open_archive_doc_button,
        Text::new("Expand documentation"),
    )
    .on_press(Message::OpenArchiveDoc);
    let close_archive_doc_button = Button::new(
        &mut slf.close_archive_doc_button,
        Text::new("Hide documentation"),
    )
    .on_press(Message::CloseArchiveDoc);
    let receive_shares_button = Button::new(
        &mut slf.receive_shares_button,
        Text::new("Receive shares").color(slf.receive_shares_button_color),
    )
    .on_press(Message::UpdateShares);
    let mut top_bar = Row::new()
        .push(archive_title)
        .push(Space::with_width(Length::Fill));
    if slf.sharing_enabled {
        top_bar = top_bar.push(receive_shares_button);
        top_bar = top_bar.push(Space::with_width(Units(8)));
    }
    top_bar = top_bar.push(archive_close_button);

    let text0 =
        Text::new("enclone visual can save sessions to the directory ~/enclone/visual/history.");
    let text1 = Text::new(
        "For a given enclone visual session:\n\
         • click the Save box on the main screen to immediately save your session\n\
         • or click the On Exit box on the main screen; it will turn red \
           (pushing again toggles state)\n\
         • when you later push the Exit button, your session will be saved.",
    );
    let text2 =
        Text::new("▒ expand - Display the commands in a saved session by checking the expand box.");
    let text3 = Text::new(
        "▒ restore - Restore a saved session by checking the restore box.  \
            This automatically saves your current session, if started.",
    );
    let text4 = Text::new("▒ delete - Delete a saved session by checking the delete box.");
    let text5 = Text::new(
        "▒ share - Share with other users by checking the share box.  You will be prompted for \
            their names.\n\
         Conversely another user may share with you.  Use the top button to receive shares.",
    );
    let text6 = Text::new(
        "▒ name - Name of a previous session is displayed:\n\
         • change it using the rectangular box and then click Apply\n\
         • up to 30 characters are allowed.",
    );
    let text7 = Text::new(
        "▒ Narrative is displayed in a box:\n\
         • if there is no narrative, a tiny box is seen\n\
         • clicking on the box will change the narrative to whatever is on your clipboard.",
    );
    let mut help_col = Column::new();
    if slf.archive_doc_open {
        help_col = help_col
            .spacing(SPACING)
            .push(text1)
            .push(text2)
            .push(text3)
            .push(text4);
        if slf.sharing_enabled {
            help_col = help_col.push(text5);
        }
        help_col = help_col.push(text6);
        help_col = help_col.push(text7);
        help_col = help_col.push(close_archive_doc_button);
    } else {
        help_col = help_col.push(open_archive_doc_button);
    }

    let mut labels = "#   date        time     expand    restore   delete".to_string();
    if slf.sharing_enabled {
        labels += "    share";
    }
    labels += "    name";
    let labels = Text::new(&labels).font(DEJAVU);

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
    }

    let mut archive_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle);

    archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
    let mut sharing = false;
    for i in 0..slf.archive_share_requested.len() {
        if slf.archive_share_requested[i] {
            sharing = true;
        }
    }
    let mut share_body = Column::new();
    if sharing {
        share_body = share_body.push(share_col());
        for (j, u) in slf.user.iter_mut().enumerate() {
            share_body = share_body.push(Space::with_height(Units(8)));
            let mut row = Row::new()
                .push(Text::new("    ").font(DEJAVU))
                .push(
                    TextInput::new(u, "", &slf.user_value[j], move |x: String| {
                        Message::UserName(x, j)
                    })
                    .width(Units(345))
                    .font(DEJAVU)
                    .padding(2),
                )
                .push(Space::with_width(Units(8)))
                .push(Checkbox::new(slf.user_selected[j], "", move |x: bool| {
                    Message::UserSelected(x, j)
                }))
                .push(Space::with_width(Units(4)));
            if !slf.user_selected[j] {
                row = row.push(Text::new("       ").font(DEJAVU));
            } else {
                if slf.user_valid[j] {
                    row = row.push(Text::new("valid  ").font(DEJAVU));
                } else {
                    row = row.push(Text::new("invalid").font(DEJAVU));
                }
            }
            share_body = share_body.push(row);
        }
        let mut valids = 0;
        for x in slf.user_valid.iter() {
            if *x {
                valids += 1;
            }
        }
        if valids > 0 {
            let row = Row::new()
                .align_items(Align::Center)
                .push(Text::new("           check to complete share").size(16))
                .push(Space::with_width(Units(14)))
                .push(Checkbox::new(slf.do_share, "", Message::DoShare))
                .push(Text::new("or uncheck share to cancel").size(16));
            share_body = share_body.push(Space::with_height(Units(8)));
            share_body = share_body.push(row);
            if slf.do_share_complete {
                share_body = share_body
                    .push(Space::with_height(Units(8)))
                    .push(
                        Text::new(
                            "           Done, your session has been shared!  \
                        You can uncheck the share box to make the share information vanish.",
                        )
                        .size(16),
                    )
                    .push(Space::with_height(Units(8)))
                    .push(
                        Text::new(
                            "           You have to uncheck the share box to do \
                        another share.",
                        )
                        .size(16),
                    );
            }
        }
    }
    let mut share_body = Some(share_body);

    let mut count = 0;
    for (i, y, q, narbut) in izip!(
        0..slf.archive_name.len(),
        slf.archive_name.iter_mut(),
        slf.archive_name_change_button.iter_mut(),
        slf.archive_narrative_button.iter_mut()
    ) {
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
                .push(Space::with_width(Units(70)))
                .push(Checkbox::new(
                    slf.restore_requested[i],
                    "",
                    move |x: bool| Message::Restore(x, i),
                ))
                .push(Space::with_width(Units(69)))
                .push(Checkbox::new(
                    slf.delete_requested[i],
                    "",
                    move |x: bool| Message::DeleteArchiveEntry(x, i),
                ))
                .push(Space::with_width(Units(68)));
            if slf.sharing_enabled {
                row = row.push(Checkbox::new(
                    slf.archive_share_requested[i],
                    "",
                    move |x: bool| Message::ArchiveShare(x, i),
                ));
                row = row.push(Space::with_width(Units(58)));
            }
            row = row.push(
                TextInput::new(y, "", &slf.archive_name_value[i], move |x: String| {
                    Message::ArchiveName(x, i)
                })
                .width(Units(300))
                .max_width(300)
                .padding(2),
            );
            row = row.push(Space::with_width(Units(8)));
            row = row.push(
                Button::new(
                    q,
                    Text::new("Apply").color(slf.archive_name_change_button_color[i]),
                )
                .on_press(Message::ArchiveNameChange(i)),
            );
            if i > 0 {
                archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
            }
            archive_scrollable = archive_scrollable.push(row);

            if slf.archive_origin[i].len() > 0 {
                archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
                archive_scrollable = archive_scrollable.push(
                    Text::new(&format!("     {}", &slf.archive_origin[i]))
                        .font(DEJAVU)
                        .size(16),
                );
            }

            // Show narrative.

            const MAX_LINE: usize = 120;
            let mut log = String::new();
            let mut rows = Vec::<Vec<String>>::new();
            let folds = fold(&slf.archive_narrative[i], MAX_LINE);
            for i in 0..folds.len() {
                rows.push(vec![folds[i].clone()]);
            }
            for i in 0..rows.len() {
                if i > 0 {
                    log += "\n";
                }
                log += &mut rows[i][0].clone();
            }
            let row = Row::new()
                .push(Text::new("    ").font(DEJAVU))
                .push(Button::new(narbut, Text::new(&log)).on_press(Message::ArchiveNarrative(i)));
            archive_scrollable = archive_scrollable.push(Space::with_height(Units(8)));
            archive_scrollable = archive_scrollable.push(row);

            // Show if share requested.

            if slf.archive_share_requested[i] {
                archive_scrollable = archive_scrollable.push(share_body.take().unwrap());
            }

            // Show restore message.

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

            // Show the commands if expansion requested.

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
        .push(help_col)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(labels)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(archive_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
