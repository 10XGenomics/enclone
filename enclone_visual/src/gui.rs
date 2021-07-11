// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use gui_structures::ComputeState::*;
use gui_structures::*;
use iced::svg::Handle;
use iced::Length::Units;
use iced::{
    Align, Application, Button, Clipboard, Color, Column, Command, Element, HorizontalAlignment,
    Image, Length, Row, Rule, Scrollable, Space, Svg, Text, TextInput, VerticalAlignment,
};
// use iced::Subscription;
use iced_aw::{Card, Modal};
// use iced_native::{window, Event};
use messages::Message;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::Duration;

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
        if !TEST_MODE.load(SeqCst) {
            (x, Command::none())
        } else {
            thread::sleep(Duration::from_millis(1000));
            (x, Command::perform(noop(), Message::RunTests))
        }
    }

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

    fn view(&mut self) -> Element<Message> {
        let text_input = TextInput::new(
            &mut self.input,
            "",
            &self.input_value,
            Message::InputChanged,
        )
        .padding(10)
        .font(DEJAVU_BOLD)
        .size(16);

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

            const COPY_BUTTON_FONT_SIZE: u16 = 15;
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

            let mut button_column2 = Column::new().spacing(8);
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

            // Add command box.

            const MAX_LINE: usize = 45;
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

            // Create summary buttons.

            let summary_button = Button::new(
                &mut self.summary_button,
                Text::new("Summary").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::OpenModalSummary);

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

        // Build the scrollable for clonotypes.

        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .style(style::ScrollableStyle)
            .push(Text::new(&self.output_value).font(DEJAVU_BOLD).size(13));

        // Fix the height of the SVG.  This needs to be set so that there is enough room for
        // the clonotype tables.  We do not set the width because it's the height that we need
        // to control.

        const SVG_HEIGHT: u16 = 400;

        // Display the SVG.
        //
        // WARNING!  When we changed the width and height to 400, the performance of scrolling
        // in the clonotype table window gradually degraded, becoming less and less responsive.
        // After a couple minutes, the app crashed, with thirty threads running.

        let svg = Svg::new(Handle::from_memory(self.svg_value.as_bytes().to_vec()))
            .height(Units(SVG_HEIGHT));
        let _svg = &svg; // to temporarily prevent warning

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
                            .map(move |_message| Message::SubmitButtonPressed(Ok(()))),
                    )
                    .height(Units(SVG_HEIGHT));
            } else {
                let svg_as_png =
                    Image::new(iced::image::Handle::from_memory(self.png_value.clone()))
                        .height(Units(SVG_HEIGHT));
                graphic_row = graphic_row.push(svg_as_png);
            }

            // Insert space to push the graphic to the left and the command column to the right.

            if !have_canvas {
                graphic_row = graphic_row.push(Space::with_width(Length::Fill));
            }

            graphic_row = graphic_row.push(command_complex);
        }

        // Put it all together.

        let left_buttons = Column::new()
            .spacing(8)
            .push(Button::new(&mut self.exit_state, Text::new("Exit")).on_press(Message::Exit))
            .push(
                Button::new(&mut self.open_state, Text::new("Help"))
                    .on_press(Message::OpenModalHelp),
            )
            .push(
                Button::new(&mut self.open_state_cookbook, Text::new("Cookbook"))
                    .on_press(Message::OpenModalCookbook),
            );

        let mut content = Column::new()
            .spacing(20)
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
                    .push(text_input)
                    .push(button)
                    .push(clear_button),
            )
            // .push(Row::new().spacing(10).push(svg))
            .push(graphic_row);
        if !self.input_history.is_empty() {
            content = content.push(Rule::horizontal(10).style(style::RuleStyle));
        }
        content = content.push(
            Row::new()
                .height(Length::Units(1000)) // Height of scrollable window, maybe??
                .align_items(Align::Center)
                .push(scrollable),
        );

        let version = VERSION.lock().unwrap()[0].clone();
        let version_float = format!("1e-{}", -version.force_f64().log10());
        Modal::new(&mut self.modal_state_help, content, move |state| {
            Card::new(
                Text::new(""),
                if !COOKBOOK.load(SeqCst) && !SUMMARY.load(SeqCst) {
                    Text::new(&format!(
                        "Welcome to enclone visual {} = {}!\n\n{}",
                        version,
                        version_float,
                        include_str!["help.txt"],
                    ))
                } else if SUMMARY.load(SeqCst) {
                    let summary = SUMMARY_CONTENTS.lock().unwrap()[0].clone();
                    let nlines = summary.chars().filter(|&n| n == '\n').count();
                    let font_size = (15 * nlines) / 38;
                    Text::new(&format!("{}", summary))
                        .font(DEJAVU_BOLD)
                        .size(font_size as u16)
                } else {
                    let preamble = "Type the tag into the input box to run the given command.\n\n";
                    Text::new(&format!(
                        "{}{}",
                        preamble,
                        COOKBOOK_CONTENTS.lock().unwrap()[0]
                    ))
                    .font(DEJAVU_BOLD)
                    .size(14)
                }
                .height(Units(800))
                .vertical_alignment(VerticalAlignment::Center),
            )
            .style(style::Help)
            .foot(
                Row::new().spacing(10).push(
                    Button::new(
                        &mut state.cancel_state,
                        Text::new("Vanish!").horizontal_alignment(HorizontalAlignment::Left),
                    )
                    .on_press(Message::CancelButtonPressed),
                ),
            )
            .width(Units(1100))
            .height(Units(1060))
            .on_close(Message::CloseModalHelp)
            .into()
        })
        .backdrop(Message::CloseModalHelp)
        .on_esc(Message::CloseModalHelp)
        .into()
    }
}
