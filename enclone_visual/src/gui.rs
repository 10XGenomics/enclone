// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use gui_structures::ComputeState::*;
use gui_structures::*;
use iced::Length::Units;
use iced::{
    Align, Application, Button, Clipboard, Color, Column, Command, Container, Element, Image,
    Length, Row, Rule, Scrollable, Space, Subscription, Text, TextInput,
};
// use iced::Subscription;
// use iced_native::{window, Event};
use iced_native::{event, subscription, window, Event};
use messages::Message;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::Duration;

fn handle_resize(width: u32, height: u32) -> Option<Message> {
    Some(Message::Resize(width, height))
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

impl Application for EncloneVisual {
    type Executor = iced::executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (EncloneVisual, Command<Message>) {
        COOKBOOK_CONTENTS.lock().unwrap().push(format_cookbook());
        let mut x = EncloneVisual::default();
        x.submit_button_text = "Submit".to_string();
        x.compute_state = WaitingForRequest;
        x.copy_image_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.cookbook = parse_cookbook();
        x.width = INITIAL_WIDTH;
        x.height = INITIAL_HEIGHT;
        if !TEST_MODE.load(SeqCst) {
            (x, Command::none())
        } else {
            thread::sleep(Duration::from_millis(1000));
            (x, Command::perform(noop(), Message::RunTests))
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message, _clipboard: &mut Clipboard) -> Command<Message> {
        self.process_message(message)
    }

    /*
    fn should_exit(&self) -> bool {
        self.should_exit
    }

    fn subscription(&self) -> Subscription<Message> {
        iced_native::subscription::events().map(Message::EventOccurred)
    }
    */

    // The subscription detects window resize events and uses those to reset
    // self.width and self.height.

    fn subscription(&self) -> Subscription<Message> {
        subscription::events_with(|event, status| {
            if let event::Status::Captured = status {
                return None;
            }
            match event {
                Event::Window(window::Event::Resized { width, height }) => {
                    handle_resize(width, height)
                }
                _ => None,
            }
        })
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    fn view(&mut self) -> Element<Message> {
        const SCROLLBAR_WIDTH: u16 = 12;
        const SPACING: u16 = 20;
        let version = VERSION.lock().unwrap()[0].clone();
        let version_float = format!("1e-{}", -version.force_f64().log10());

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Handle the summary case.

        if self.summary_mode {
            let summary_title = Text::new(&format!("Summary")).size(30);
            let summary = SUMMARY_CONTENTS.lock().unwrap()[0].clone();
            let nlines = summary.chars().filter(|&n| n == '\n').count();
            let font_size = (15 * nlines) / 38;
            let summary_scrollable = Scrollable::new(&mut self.scroll)
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
            let summary_close_button = Button::new(&mut self.open_state, Text::new("Vanish!"))
                .on_press(Message::SummaryClose);
            let content = Column::new()
                .spacing(SPACING)
                .padding(20)
                .push(summary_title)
                .push(summary_scrollable)
                .push(summary_close_button);
            return Container::new(content)
                .width(Length::Fill)
                .height(Length::Fill)
                .into();
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Handle the help case.

        if self.help_mode {
            let help_title = Text::new(&format!("Help")).size(30);
            let help_close_button = Button::new(&mut self.open_state, Text::new("Dismiss"))
                .on_press(Message::HelpClose);
            let top_bar = Row::new()
                .push(help_title)
                .push(Space::with_width(Length::Fill))
                .push(help_close_button);
            let max_width = Units((self.width - 60) as u16);
            let png_input_region = include_bytes!("../images/input_region.png").to_vec();
            let input_region =
                Image::new(iced::image::Handle::from_memory(png_input_region)).width(max_width);
            let png_history_region = include_bytes!("../images/history_region.png").to_vec();
            let history_region =
                Image::new(iced::image::Handle::from_memory(png_history_region)).height(Units(240));
            let png_middle_region = include_bytes!("../images/middle_region.png").to_vec();
            let middle_region = Image::new(iced::image::Handle::from_memory(png_middle_region))
                .height(Units(300))
                .width(Units(290));
            let png_top_region = include_bytes!("../images/top_region.png").to_vec();
            let top_region =
                Image::new(iced::image::Handle::from_memory(png_top_region)).height(Units(120));
            let help_scrollable = Scrollable::new(&mut self.scroll)
                .width(Length::Fill)
                .height(Length::Fill)
                .scrollbar_width(SCROLLBAR_WIDTH)
                .scroller_width(12)
                .style(style::ScrollableStyle)
                //
                // Intro.
                //
                .push(Space::with_height(Units(10)))
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
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                            corner of the screen.",
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
                                .push(Text::new("3.  Cookbook, to show some sample commands.")),
                        ),
                )
                //
                // Layout.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
                .push(Space::with_height(Units(20)))
                .push(Text::new("Overall layout").size(24))
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new("There are input boxes near the top (described next).")
                        .width(max_width),
                )
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
                .push(
                    Text::new("2.  A text subwindow, which typically has clonotypes in it.")
                        .width(max_width),
                )
                //
                // Input.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
                .push(Space::with_height(Units(20)))
                .push(Text::new("Entering input").size(24))
                .push(Space::with_height(Units(20)))
                .push(input_region)
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new(
                        "Above, you can see two boxes.  You can type a command into \
                    these.  The reason for having two boxes is that it allows for longer \
                    commands: you can split a command between the two boxes",
                    )
                    .width(max_width),
                )
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new(
                        "Except for special cases (see below), every command begins with \
                    the word enclone.  You can see examples by pushing the Cookbook button on the \
                    main screen.  You can learn about enclone commands in general by going to \
                    the site bit.ly/enclone.",
                    )
                    .width(max_width),
                )
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new("Once you've entered your command, push the Submit button.")
                        .width(max_width),
                )
                //
                // Special commands.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                .push(Text::new("also work"))
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new(
                        "If you've displayed a honeycomb plot (see cookbook for examples), \
                    then positioning your mouse over a cell will cause a box to appear that \
                    provides some information about that cell.",
                    )
                    .width(max_width),
                )
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new(
                        "And clicking on a cell is the same as typing its number into \
                    the input box!",
                    )
                    .width(max_width),
                )
                .push(Space::with_height(Units(20)))
                .push(
                    Text::new("Group ids are converted into a special enclone argument")
                        .width(max_width),
                )
                .push(Space::with_height(Units(10)))
                .push(Text::new("G=...").font(DEJAVU_BOLD).size(20))
                .push(Space::with_height(Units(10)))
                .push(Text::new("that can also be supplied to enclone.").width(max_width))
                //
                // History.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                                    .width(max_width),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "On the right, you can see boxes, that will appear on \
                            the right of your screen once you've entered your first command.",
                                    )
                                    .width(max_width),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "Initially, some of the boxes will be blank, meaning \
                            that they don't make sense yet and won't do anything.",
                                    )
                                    .width(max_width),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(Text::new("There are four boxes:"))
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "• The number at the top is the index of the current \
                            state.  This is not for pushing.",
                                    )
                                    .width(max_width),
                                )
                                .push(
                                    Text::new(
                                        "• Push the up arrow to go back to the previous state, \
                            meaning the last command that you typed.",
                                    )
                                    .width(max_width),
                                )
                                .push(
                                    Text::new(
                                        "• Push the down arrow to go forward to the next state.",
                                    )
                                    .width(max_width),
                                )
                                .push(Text::new(
                                    "• Push the Del button to delete the current state, and go \
                            backward, if that makes sense.",
                                ))
                                .width(max_width),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(history_region),
                )
                //
                // The middle boxes.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                                    .width(Units((self.width - 350) as u16)),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "The top box is the translated command.  It is the same \
                            as the command you typed, unless you used a special command.",
                                    )
                                    .width(Units((self.width - 350) as u16)),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "Below it is a button to copy the command to your \
                            clipboard.  This copied command can be reentered in enclone visual, \
                            or supplied to \"regular\" enclone",
                                    )
                                    .width(Units((self.width - 350) as u16)),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "Next there is a button to copy the graphics image to \
                            your clipboard, assuming that you have a graphics image.",
                                    )
                                    .width(Units((self.width - 350) as u16)),
                                )
                                .push(Space::with_height(Units(20)))
                                .push(
                                    Text::new(
                                        "And at the bottom is a button to display the summary \
                            stats for your enclone command.",
                                    )
                                    .width(Units((self.width - 350) as u16)),
                                ),
                        )
                        .push(Space::with_height(Units(20)))
                        .push(middle_region),
                )
                //
                // Smarts.
                //
                .push(Space::with_height(Units(20)))
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                .push(Rule::horizontal(10).style(style::RuleStyle))
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
                .push(Rule::horizontal(10).style(style::RuleStyle))
                .push(help_scrollable);
            return Container::new(content)
                .width(Length::Fill)
                .height(Length::Fill)
                .into();
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Handle the cookbook case.

        if self.cookbook_mode {
            let cookbook_title = Text::new(&format!("Cookbook")).size(30);
            let preamble = "Type the tag into the input box to run the given command.\n\n";
            let cookbook_scrollable = Scrollable::new(&mut self.scroll)
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
            let cookbook_close_button = Button::new(&mut self.open_state, Text::new("Vanish!"))
                .on_press(Message::CookbookClose);
            let content = Column::new()
                .spacing(SPACING)
                .padding(20)
                .push(cookbook_title)
                .push(cookbook_scrollable)
                .push(cookbook_close_button);
            return Container::new(content)
                .width(Length::Fill)
                .height(Length::Fill)
                .into();
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Now do everything else.

        let text_input1 = TextInput::new(
            &mut self.input1,
            "",
            &self.input1_value,
            Message::InputChanged1,
        )
        .padding(7)
        .font(DEJAVU_BOLD)
        .size(16);
        let text_input2 = TextInput::new(
            &mut self.input2,
            "",
            &self.input2_value,
            Message::InputChanged2,
        )
        .padding(7)
        .font(DEJAVU_BOLD)
        .size(16);
        let text_input_column = Column::new()
            .spacing(5)
            .width(iced::Length::Fill)
            .push(text_input1)
            .push(text_input2);

        let button = Button::new(
            &mut self.button,
            Text::new(if self.compute_state == WaitingForRequest {
                "Submit"
            } else {
                "thinking"
            }),
        )
        .padding(10)
        .on_press(Message::SubmitButtonPressed(Ok(())));
        let clear_button = Button::new(&mut self.clear_button, Text::new("Clear"))
            .padding(10)
            .on_press(Message::ClearButtonPressed);

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Define the button complex that is the "control panel".

        let mut command_complex = Row::new().spacing(10);
        {
            const FB_BUTTON_FONT_SIZE: u16 = 45;
            let back_button = Button::new(
                &mut self.back_button,
                Text::new("⇧").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::BackButtonPressed(Ok(())));

            let forward_button = Button::new(
                &mut self.forward_button,
                Text::new("⇩").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::ForwardButtonPressed(Ok(())));

            const COPY_BUTTON_FONT_SIZE: u16 = 15;
            let del_button = Button::new(
                &mut self.del_button,
                Text::new("Del ").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DelButtonPressed(Ok(())));

            let null_button1 = Button::new(
                &mut self.null_button1,
                Text::new(" ").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);

            let null_button2 = Button::new(
                &mut self.null_button2,
                Text::new(" ").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);

            let copy_image_button = Button::new(
                &mut self.copy_image_button,
                Text::new("Copy image")
                    .size(COPY_BUTTON_FONT_SIZE)
                    .color(self.copy_image_button_color),
            )
            .on_press(Message::GraphicsCopyButtonPressed);

            let null_copy_image_button = Button::new(
                &mut self.null_button3,
                Text::new("          ").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::GraphicsCopyButtonPressed);

            let mut state_pos = format!("{}", self.history_index);
            if state_pos.len() == 1 {
                state_pos += "    ";
            } else if state_pos.len() == 2 {
                state_pos += "  ";
            }
            let state_pos_button = Button::new(
                &mut self.state_pos_button_null,
                Text::new(&state_pos).size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);

            let mut button_column2 = Column::new().spacing(8);
            button_column2 = button_column2.push(state_pos_button);
            if self.history_index > 1 {
                button_column2 = button_column2.push(back_button);
            } else {
                button_column2 = button_column2.push(null_button1);
            }
            if self.history_index < self.svg_history.len() {
                button_column2 = button_column2.push(forward_button);
            } else {
                button_column2 = button_column2.push(null_button2);
            }
            button_column2 = button_column2.push(del_button);

            // Add command box.

            const MAX_LINE: usize = 35;
            let mut log = String::new();
            if self.history_index >= 1 {
                let cmd = self.translated_input_hist_uniq
                    [self.translated_input_history[self.history_index - 1]]
                    .clone();
                let mut rows = Vec::<Vec<String>>::new();
                let folds = fold(&cmd, MAX_LINE);
                for i in 0..folds.len() {
                    rows.push(vec![folds[i].clone()]);
                }
                for i in 0..rows.len() {
                    if i > 0 {
                        log += "\n";
                    }
                    log += &mut rows[i][0].clone();
                }
            }

            // Create summary button.

            let summary_button = Button::new(
                &mut self.summary_button,
                Text::new("Summary").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::SummaryOpen);

            // Build the command column.

            let mut col = Column::new().spacing(8).align_items(Align::End);
            col = col
                .push(
                    Button::new(
                        &mut self.null_button,
                        Text::new(&log).font(DEJAVU_BOLD).size(12),
                    )
                    .on_press(Message::DoNothing),
                )
                .push(
                    Button::new(
                        &mut self.command_copy_button,
                        Text::new("Copy command").size(COPY_BUTTON_FONT_SIZE),
                    )
                    .on_press(Message::CommandCopyButtonPressed),
                );
            if self.history_index >= 1 && !self.is_blank[self.history_index - 1] {
                col = col.push(copy_image_button);
            } else {
                col = col.push(null_copy_image_button);
            }
            col = col.push(summary_button);

            // Add the command column to the row.

            command_complex = command_complex.push(col);

            // Add up and down arrows.

            command_complex = command_complex.push(button_column2);
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Build the scrollable for clonotypes.  We truncate lines to prevent wrapping.

        const CLONOTYPE_FONT_SIZE: u16 = 13;
        let font_width = CLONOTYPE_FONT_SIZE as f32 * 0.5175;
        let available = self.width - (3 * SPACING + SCROLLBAR_WIDTH) as u32;
        let nchars = (available as f32 / font_width).round() as usize;
        let mut trunc = String::new();
        for line in self.output_value.lines() {
            for (i, c) in line.chars().enumerate() {
                if i == nchars {
                    break;
                }
                trunc.push(c);
            }
            trunc.push('\n');
        }
        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(SCROLLBAR_WIDTH)
            .scroller_width(12)
            .style(style::ScrollableStyle)
            .push(
                Text::new(&trunc)
                    .font(DEJAVU_BOLD)
                    .size(CLONOTYPE_FONT_SIZE),
            );

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Fix the height of the SVG.  This needs to be set so that there is enough room for
        // the clonotype tables.  We do not set the width because it's the height that we need
        // to control.

        // let svg_height = if !self.is_blank_current() {
        let mut blank = false;
        if self.history_index > 0 {
            blank = self.is_blank[self.history_index - 1];
        }
        let svg_height = if !blank { SVG_HEIGHT } else { SVG_NULL_HEIGHT };

        // Display the SVG.

        let png_banner = include_bytes!("../../img/enclone_banner.png").to_vec();
        let banner = Image::new(iced::image::Handle::from_memory(png_banner)).width(Units(500));

        let have_canvas = self.canvas_view.state.geometry_value.is_some();
        let mut graphic_row = Row::new().spacing(10);
        if self.svg_value.len() > 0 {
            // Show the graphic.

            if have_canvas {
                graphic_row = graphic_row
                    .push(
                        self.canvas_view
                            .view()
                            .map(move |message| Message::GroupClicked(message)),
                    )
                    .height(Units(svg_height));
            } else {
                let svg_as_png =
                    Image::new(iced::image::Handle::from_memory(self.png_value.clone()))
                        .height(Units(svg_height));
                graphic_row = graphic_row.push(svg_as_png);
            }

            // Insert space to push the graphic to the left and the command column to the right.

            if !have_canvas {
                graphic_row = graphic_row.push(Space::with_width(Length::Fill));
            }

            graphic_row = graphic_row.push(command_complex);
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Put it all together.

        let left_buttons = Column::new()
            .spacing(8)
            .push(Button::new(&mut self.exit_state, Text::new("Exit")).on_press(Message::Exit))
            .push(Button::new(&mut self.open_state, Text::new("Help")).on_press(Message::HelpOpen))
            .push(
                Button::new(&mut self.open_state_cookbook, Text::new("Cookbook"))
                    .on_press(Message::CookbookOpen),
            );
        let mut content = Column::new()
            .spacing(SPACING)
            .padding(20)
            .max_width(1500) // this governs the max window width upon manual resizing
            .push(
                Row::new()
                    .spacing(100)
                    .align_items(Align::Center)
                    .push(left_buttons)
                    .push(banner),
            )
            .push(
                Row::new()
                    .spacing(10)
                    .align_items(Align::Center)
                    .push(text_input_column)
                    .push(button)
                    .push(clear_button),
            )
            // .push(Row::new().spacing(10).push(svg))
            .push(graphic_row);
        if self.svg_history.len() > 0 {
            content = content.push(Rule::horizontal(10).style(style::RuleStyle));
        }
        content = content.push(
            Row::new()
                .height(Length::Units(1000)) // Height of scrollable window, maybe??
                .align_items(Align::Center)
                .push(scrollable),
        );
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}
