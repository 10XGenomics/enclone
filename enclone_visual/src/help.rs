// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
use iced::Length::Units;
use iced::{
    Alignment, Button, Column, Container, Element, Image, Length, Row, Rule, Scrollable, Space,
    Text,
};
use messages::Message;

pub fn help(slf: &mut gui_structures::EncloneVisual) -> Element<Message> {
    let version = VERSION.lock().unwrap()[0].clone();
    let help_title = Text::new(&format!("Help")).font(LIBERATION_SANS).size(30);
    let help_close_button = Button::new(
        &mut slf.open_state,
        Text::new("Dismiss").font(LIBERATION_SANS),
    )
    .on_press(Message::HelpClose(Ok(())));
    let top_bar = Row::new()
        .push(help_title)
        .push(Space::with_width(Length::Fill))
        .push(help_close_button);
    let max_width = Units((slf.width - 60) as u16);
    let max_width2 = Units((slf.width - 250) as u16);
    let png_input_region = include_bytes!("../images/input_region.png").to_vec();
    let input_region =
        Image::new(iced::image::Handle::from_memory(png_input_region)).width(max_width);
    let png_history_region = include_bytes!("../images/history_region.png").to_vec();
    let history_region =
        Image::new(iced::image::Handle::from_memory(png_history_region)).height(Units(240));
    let png_right_region = include_bytes!("../images/right_region.png").to_vec();
    let right_region =
        Image::new(iced::image::Handle::from_memory(png_right_region)).height(Units(115));
    let png_middle_region = include_bytes!("../images/middle_region.png").to_vec();
    let middle_region = Image::new(iced::image::Handle::from_memory(png_middle_region))
        .height(Units(300))
        .width(Units(290));
    let png_top_region = include_bytes!("../images/top_region.png").to_vec();
    let top_region = Image::new(iced::image::Handle::from_memory(png_top_region)).height(Units(84));
    let help_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        //
        // Intro.
        //
        .push(Text::new("Introduction").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(&format!("Welcome to enclone visual version {}!", version,))
                .font(LIBERATION_SANS),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "enclone visual is a semi-graphical \
            version of enclone.  You can find out more about enclone at the site",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(Text::new("bit.ly/enclone.").font(DEJAVU_BOLD))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "enclone visual simultaneously displays the text and graphical \
            output that enclone can produce.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "This is alpha software: there are many more bugs than there \
            are in enclone.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        //
        // Top.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Top of the page").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Row::new()
                .push(top_region)
                .push(Space::with_width(Units(15)))
                .push(
                    Column::new()
                        .push(
                            Text::new(
                                "Here are two buttons that appear in the upper left \
                    corner of the screen:",
                            )
                            .font(LIBERATION_SANS),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "1.  Exit, to leave enclone visual.  Note that the \
                    circular red button in the extreme upper left corner is busted.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width),
                        )
                        .push(Text::new("2.  Help, to get to this page.").font(LIBERATION_SANS)),
                ),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Row::new()
                .align_items(Alignment::Center)
                .push(
                    Column::new()
                        .push(Space::with_height(Units(10)))
                        .push(
                            Text::new(
                                "Here are six buttons that appear in the upper right \
                             corner of the screen:",
                            )
                            .font(LIBERATION_SANS),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "1.  Snapshot, to copy a screenshot of the entire window \
                                 to the clipboard.  This button is only present on Macs.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "2.  Console, to show what's in the terminal window.  \
                                 At the moment this is mostly of interest to developers.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        )
                        .push(
                            Text::new("3.  Save, to cause the session to be saved.")
                                .font(LIBERATION_SANS)
                                .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "4.  On Exit, to cause the session to be saved \
                                 when the Exit button is pushed.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "5.  Tooltip, to toggle the position of the tooltip text window.  \
                                 Clicking causes clockwise rotation of the position amongst the \
                                 four corners of the graphics window.  One position may yield a \
                                 better overall view than another.  If you click, you'll see the \
                                 button flash, but no other change, until you mouse over a cell.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "6.  Archive.  Opens a page to allow restoration or sharing of a \
                                 previous session.",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        )
                        .push(Space::with_height(Units(10)))
                        .push(
                            Text::new(
                                "More information about saving and restoring session may be \
                                obtained by \
                                pushing the Archive button on the main page.  The archive page \
                                also provides access to cookbooks, which you should work through!",
                            )
                            .font(LIBERATION_SANS)
                            .width(max_width2),
                        ),
                )
                .push(Space::with_width(Length::Fill))
                .push(right_region)
                .push(Space::with_width(Units(20))),
        )
        //
        // Layout.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Overall layout").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("There are input boxes near the top (described next).")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Once you've typed your first command, the screen will \
            split into two main parts:",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("1.  A graphics subwindow, which may or may not be populated.")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "2.  A text subwindow, which typically has clonotypes in it.  \
                To make scrolling manageable, only the first fifty clonotypes are displayed.  \
                Please see the special commands section for how to see a different set of \
                clonotypes.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        //
        // Input.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Entering input").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(input_region)
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Above, you can see two boxes.  You can type a command into \
            these.  The reason for having two boxes is that it allows for longer \
            commands: you can split a command between the two boxes.  Pushing the Cmd button \
            (see below) will provide even more boxes.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Except for special cases (see below), every command begins with \
            the word enclone.  You can see examples \
            in the cookbooks on the Archive page.  You can learn about \
            enclone commands in general by going to the site bit.ly/enclone.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("Once you've entered your command, push the Submit button.").width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Then wait for the command to finish!  Some buttons are disabled during \
                this period.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        //
        // Special commands.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Special commands").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "You can type a number into the text box, where the number \
            is the number of a clonotype group.  Things like this",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(10)))
        .push(Text::new("1,7,10-15").font(DEJAVU_BOLD).size(20))
        .push(Space::with_height(Units(10)))
        .push(Text::new("also work.").font(LIBERATION_SANS))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "If you've displayed a honeycomb plot (see cookbooks for examples), \
            then positioning your mouse over a cell will cause a \"tooltip\" box to appear that \
            provides some information about that cell.  See also the Tooltip button, that \
            controls the position of this box.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Clicking on a cell is the same as typing its number into \
            the input box, thus displaying the corresponding clonotype table.  In addition, \
            clicking on a cell will copy the tooltip text to your clipboard.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("Group ids are converted into a special enclone argument")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        .push(Space::with_height(Units(10)))
        .push(Text::new("G=...").font(DEJAVU_BOLD).size(20))
        .push(Space::with_height(Units(10)))
        .push(
            Text::new("that can also be supplied to enclone.  In addition, G=all works.")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        //
        // History.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("History, AKA the time machine")
                .font(LIBERATION_SANS)
                .size(24),
        )
        .push(Space::with_height(Units(15)))
        .push(
            Row::new()
                .push(
                    Column::new()
                        .push(Space::with_height(Units(5)))
                        .push(
                            Text::new(
                                "enclone visual remembers your previous commands and \
                    their outputs.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "On the right, you can see boxes, that will appear on \
                    the right of your screen once you've entered your first command.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Initially, some of the boxes will be blank, meaning \
                    that they don't make sense yet and won't do anything.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(Text::new("There are four boxes:").font(LIBERATION_SANS))
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "• The number at the top is the index of the current \
                    state.  This is not a button that can be pushed.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(
                            Text::new(
                                "• Push the up arrow to go back to the previous state, \
                    meaning the last command that you typed.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(
                            Text::new("• Push the down arrow to go forward to the next state.")
                                .font(LIBERATION_SANS)
                                .width(Units((slf.width - 120) as u16)),
                        )
                        .push(
                            Text::new(
                                "• Push the Del button to delete the current state, and go \
                    backward, if that makes sense.",
                            )
                            .font(LIBERATION_SANS),
                        )
                        .width(Units((slf.width - 120) as u16)),
                )
                .push(Space::with_height(Units(20)))
                .push(history_region),
        )
        //
        // The middle boxes.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("The middle boxes").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(15)))
        .push(
            Row::new()
                .push(
                    Column::new()
                        .push(Space::with_height(Units(5)))
                        .push(
                            Text::new(
                                "Just to the left of the history boxes are more, \
                    samples of which you can see on the right.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "The top box is the translated command.  It is the same \
                    as the command you typed, unless you used a special command.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Below it there is a button to copy the graphics image to \
                    your clipboard, assuming that you have a graphics image.  The default width \
                    of the image is 2000 pixels.  To obtain lower or higher resolution, copy a \
                    number between 1000 and 4000 onto your clipboard, before pushing Copy image.  \
                    Then that number will be used as the width.  Note also that direct image copy \
                    from the screen can be a good alternative.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Next to it is a button to copy the command to your \
                    clipboard.  This copied command can be reentered in enclone visual, \
                    or supplied to \"regular\" enclone, so long as you change instances of gui \
                    to actual file names.  You also need to add double quotes \
                    around arguments containing certain characters including < and >.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Below that is a button to open a window displaying just the \
                                graphic, and one for just the clonotypes, \
                                a button to display the summary stats for your enclone command, \
                                and a button for entering long commands.  Each of these four \
                                buttons opens a separate page.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Finally near the bottom is the narrative box.  This is text for \
                    a given state that you can set by first copying text to your clipboard \
                    (outside enclone visual), and then clicking on the narrative box to copy \
                    the text into it.  Conversely, there is a Copy button below the narrative box \
                    that will copy the narrative text onto your clipboard.",
                            )
                            .font(LIBERATION_SANS)
                            .width(Units((slf.width - 350) as u16)),
                        ),
                )
                .push(Space::with_height(Units(20)))
                .push(middle_region),
        )
        //
        // Smarts.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Smarts").font(LIBERATION_SANS).size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "If you run a command, and then run a similar command after it, \
            enclone visual may elide some calculations from the previous command, so as \
            to respond faster.  This capability is not pushed as far as it could be.",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        //
        // Limitations.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("Limitations, AKA big bugs")
                .font(LIBERATION_SANS)
                .size(24),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "There are two main limitations of the current version of enclone \
            visual:",
            )
            .font(LIBERATION_SANS)
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("1.  The clonotype tables are black and white.")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        .push(
            Text::new("2.  You can't use the mouse to copy text, except from text input boxes.")
                .font(LIBERATION_SANS)
                .width(max_width),
        )
        //
        // Bottom.
        //
        .push(Space::with_height(Units(20)));

    let content = Column::new()
        .spacing(SPACING)
        .padding(20)
        .push(top_bar)
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(help_scrollable);
    Container::new(content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
