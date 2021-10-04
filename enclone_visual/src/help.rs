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
    let version_float = format!("1e-{}", -version.force_f64().log10());
    let help_title = Text::new(&format!("Help")).size(30);
    let help_close_button =
        Button::new(&mut slf.open_state, Text::new("Dismiss")).on_press(Message::HelpClose(Ok(())));
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
    let top_region =
        Image::new(iced::image::Handle::from_memory(png_top_region)).height(Units(120));
    let help_scrollable = Scrollable::new(&mut slf.scroll)
        .width(Length::Fill)
        .height(Length::Fill)
        .scrollbar_width(SCROLLBAR_WIDTH)
        .scroller_width(12)
        .style(style::ScrollableStyle)
        //
        // Intro.
        //
        .push(Text::new("Introduction").size(24))
        .push(Space::with_height(Units(20)))
        .push(Text::new(&format!(
            "Welcome to enclone visual version {} = {}!",
            version, version_float,
        )))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "enclone visual is a semi-graphical \
            version of enclone.  You can find out more about enclone at the site",
            )
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
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "This is alpha software: there are many more bugs than there \
            are in enclone.",
            )
            .width(max_width),
        )
        //
        // Top.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Top of the page").size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Row::new()
                .push(top_region)
                .push(Space::with_width(Units(15)))
                .push(
                    Column::new()
                        .push(Space::with_height(Units(10)))
                        .push(Text::new(
                            "Here are three buttons that appear in the upper left \
                    corner of the screen:",
                        ))
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "1.  Exit, to leave enclone visual.  Note that the \
                    circular red button in the extreme upper left corner is busted.",
                            )
                            .width(max_width),
                        )
                        .push(Text::new("2.  Help, to get to this page."))
                        .push(Text::new(
                            "3.  Cookbook, to show some sample commands.  And \
                            see below, at Archive.",
                        )),
                ),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Row::new()
                .align_items(Alignment::Center)
                .push(
                    Column::new()
                        .push(Space::with_height(Units(10)))
                        .push(Text::new(
                            "Here are six buttons that appear in the upper right \
                             corner of the screen:",
                        ))
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "1.  Snapshot, to copy a screenshot of the entire window \
                                 to the clipboard.",
                            )
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "2.  Console, to show what's in the terminal window.  \
                                 At the moment this is mostly of interest to developers.",
                            )
                            .width(max_width2),
                        )
                        .push(
                            Text::new("3.  Save, to cause the session to be saved.")
                                .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "4.  On Exit, to cause the session to be saved \
                                 when the Exit button is pushed.",
                            )
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
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "6.  Archive.  Opens a page to allow restoration or sharing of a \
                                 previous session.",
                            )
                            .width(max_width2),
                        )
                        .push(
                            Text::new(
                                "More information for the last two buttons may be obtained by \
                                pushing the Archive button on the main page.",
                            )
                            .width(max_width2),
                        )
                        .push(
                            Text::new("The archive page also gives you access to more cookbooks!")
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
        .push(Text::new("Overall layout").size(24))
        .push(Space::with_height(Units(20)))
        .push(Text::new("There are input boxes near the top (described next).").width(max_width))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Once you've typed your first command, the screen will \
            split into two main parts:",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("1.  A graphics subwindow, which may or may not be populated.")
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
            .width(max_width),
        )
        //
        // Input.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Entering input").size(24))
        .push(Space::with_height(Units(20)))
        .push(input_region)
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Above, you can see two boxes.  You can type a command into \
            these.  The reason for having two boxes is that it allows for longer \
            commands: you can split a command between the two boxes.  If you have a particularly \
            long command, you can also widen the enclone visual window.",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Except for special cases (see below), every command begins with \
            the word enclone.  You can see examples by pushing the Cookbook button on the \
            main screen, and from the cookbooks at the Archive page.  You can learn about \
            enclone commands in general by going to the site bit.ly/enclone.",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new("Once you've entered your command, push the Submit button.").width(max_width),
        )
        //
        // Special commands.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Special commands").size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "In the cookbook, you'll find abbreviations for commands, \
            called tags, for example #1.  You can type these into the input box.",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "You can also type a number into the text box, where the number \
            is the number of a clonotype group.  Things like this",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(10)))
        .push(Text::new("1,7,10-15").font(DEJAVU_BOLD).size(20))
        .push(Space::with_height(Units(10)))
        .push(Text::new("also work."))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "If you've displayed a honeycomb plot (see cookbook for examples), \
            then positioning your mouse over a cell will cause a \"tooltip\" box to appear that \
            provides some information about that cell.  See also the Tooltip button, that \
            controls the position of this box.",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "Clicking on a cell is the same as typing its number into \
            the input box, thus displaying the corresponding clonotype table.  In addition, \
            clicking on a cell will copy the tooltip text to your clipboard.",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(Text::new("Group ids are converted into a special enclone argument").width(max_width))
        .push(Space::with_height(Units(10)))
        .push(Text::new("G=...").font(DEJAVU_BOLD).size(20))
        .push(Space::with_height(Units(10)))
        .push(
            Text::new("that can also be supplied to enclone.  In addition, G=all works.")
                .width(max_width),
        )
        //
        // History.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("History, AKA the time machine").size(24))
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
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "On the right, you can see boxes, that will appear on \
                    the right of your screen once you've entered your first command.",
                            )
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Initially, some of the boxes will be blank, meaning \
                    that they don't make sense yet and won't do anything.",
                            )
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(Text::new("There are four boxes:"))
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "• The number at the top is the index of the current \
                    state.  This is not for pushing.",
                            )
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(
                            Text::new(
                                "• Push the up arrow to go back to the previous state, \
                    meaning the last command that you typed.",
                            )
                            .width(Units((slf.width - 120) as u16)),
                        )
                        .push(
                            Text::new("• Push the down arrow to go forward to the next state.")
                                .width(Units((slf.width - 120) as u16)),
                        )
                        .push(Text::new(
                            "• Push the Del button to delete the current state, and go \
                    backward, if that makes sense.",
                        ))
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
        .push(Text::new("The middle boxes").size(24))
        .push(Space::with_height(Units(15)))
        .push(
            Row::new()
                .push(
                    Column::new()
                        .push(Space::with_height(Units(5)))
                        .push(
                            Text::new(
                                "Just to the left of the history boxes are some more, \
                    samples of which you can see on the right.",
                            )
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "The top box is the translated command.  It is the same \
                    as the command you typed, unless you used a special command.",
                            )
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
                            .width(Units((slf.width - 350) as u16)),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(
                            Text::new(
                                "Below that is a button to open a window displaying just the \
                                graphic, and one for just the clonotypes, and a button to display \
                                the summary stats for your enclone command.",
                            )
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
        .push(Text::new("Smarts").size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "If you run a command, and then run a similar command after it, \
            enclone visual may elide some calculations from the previous command, so as \
            to respond faster.  This capability is not pushed as far as it could be.",
            )
            .width(max_width),
        )
        //
        // Limitations.
        //
        .push(Space::with_height(Units(20)))
        .push(Rule::horizontal(10).style(style::RuleStyle2))
        .push(Space::with_height(Units(20)))
        .push(Text::new("Limitations, AKA big bugs").size(24))
        .push(Space::with_height(Units(20)))
        .push(
            Text::new(
                "There are two main limitations of the current version of enclone \
            visual:",
            )
            .width(max_width),
        )
        .push(Space::with_height(Units(20)))
        .push(Text::new("1.  The clonotype tables are black and white.").width(max_width))
        .push(
            Text::new(
                "2.  You can't use the mouse to copy text from the graphics \
            window or the text window.",
            )
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
