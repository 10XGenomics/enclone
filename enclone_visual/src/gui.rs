// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use iced::svg::Handle;
use iced::Length::Units;
use iced::{
    button, scrollable, text_input, Align, Button, Color, Column, Element, Font,
    HorizontalAlignment, Image, Length, Row, Rule, Sandbox, Scrollable, Settings, Svg, Text,
    TextInput, VerticalAlignment,
};
use iced_aw::{modal, Card, Modal};
use perf_stats::*;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant};
use string_utils::*;

const DEJAVU: Font = Font::External {
    name: "DEJAVU",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono-Bold.ttf"),
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn launch_gui() {
    let mut settings = Settings::default();
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = (1100 as u32, 1060 as u32); // reasonable minimum size
    settings.window = window_settings;
    let _ = EncloneVisual::run(settings);
}

#[derive(Default)]
struct EncloneVisual {
    scroll: scrollable::State,
    input: text_input::State,
    input_value: String,
    output_value: String,
    svg_value: String,
    button: button::State,

    open_state: button::State,
    modal_state: modal::State<ModalState>,
    last_message: Option<Message>,
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
    OpenModal,
    CloseModal,
    CancelButtonPressed,
}

#[derive(Default)]
struct ModalState {
    cancel_state: button::State,
}

impl Sandbox for EncloneVisual {
    type Message = Message;

    fn new() -> Self {
        EncloneVisual::default()
    }

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::OpenModal => self.modal_state.show(true),
            Message::CloseModal => self.modal_state.show(false),
            Message::CancelButtonPressed => self.modal_state.show(false),
            Message::InputChanged(ref value) => self.input_value = value.to_string(),
            Message::ButtonPressed => {
                // Need to figure what to do if we are already processing a request, for example
                // if the user pushes the button twice or enters a second command and pushes the
                // button before the first one has completed.  For now, do nothing.

                if !PROCESSING_REQUEST.load(SeqCst) {
                    let t = Instant::now();
                    USER_REQUEST.lock().unwrap().clear();
                    USER_REQUEST.lock().unwrap().push(self.input_value.clone());
                    PROCESSING_REQUEST.store(true, SeqCst);
                    while PROCESSING_REQUEST.load(SeqCst) {
                        thread::sleep(Duration::from_millis(10));
                    }
                    let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                    if reply_text.contains("enclone failed") {
                        reply_text =
                            format!("enclone failed{}", reply_text.after("enclone failed"));
                    }
                    reply_text += "\n \n"; // papering over truncation bug
                    let reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                    self.output_value = reply_text.to_string();
                    self.svg_value = reply_svg.to_string();
                    println!(
                        "time used processing command = {:.1} seconds\n",
                        elapsed(&t)
                    );
                }
            }
        }
        self.last_message = Some(message)
    }

    fn view(&mut self) -> Element<Message> {
        let text_input = TextInput::new(
            &mut self.input,
            "",
            &self.input_value,
            Message::InputChanged,
        )
        .padding(10)
        .size(14);

        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);

        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .style(style::Squeak)
            .push(Text::new(&self.output_value).font(DEJAVU).size(13));

        // Display the SVG.
        //
        // WARNING!  When we changed the width and height to 400, the performance of scolling
        // in the clonotype table window gradually degraded, becoming less and less responsive.
        // After a couple minutes, the app crashed, with thirty threads running.

        let svg = Svg::new(Handle::from_memory(self.svg_value.as_bytes().to_vec()))
            .width(Units(300))
            .height(Units(300));

        let png = include_bytes!("../../img/enclone_banner.png").to_vec();
        let banner = Image::new(iced::image::Handle::from_memory(png)).width(Units(500));

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .max_width(1500) // this governs the max window width upon manual resizing
            .push(
                Row::new()
                    .spacing(230)
                    .align_items(Align::Center)
                    .push(
                        Button::new(&mut self.open_state, Text::new("Help"))
                            .on_press(Message::OpenModal),
                    )
                    .push(banner),
            )
            .push(Row::new().spacing(10).push(text_input).push(button))
            .push(Row::new().spacing(10).push(svg))
            .push(Rule::horizontal(10).style(style::RuleStyle))
            .push(
                Row::new()
                    .height(Length::Units(1000)) // Height of scrollable window, maybe??
                    .align_items(Align::Center)
                    .push(scrollable),
            );

        use iced_aw::style::{
            card::{Style, StyleSheet},
            colors,
        };

        #[derive(Clone, Copy)]
        pub struct Gerbil;

        impl StyleSheet for Gerbil {
            fn active(&self) -> Style {
                Style {
                    background: iced::Background::Color(Color::from_rgb(0.9, 1.0, 0.9)),
                    border_width: 0.0,
                    border_color: iced::Color::from_rgb(1.0, 1.0, 1.0),
                    head_background: iced::Background::Color(Color::from_rgb(0.9, 1.0, 0.9)),
                    head_text_color: colors::WHITE,
                    close_color: colors::WHITE,
                    ..Style::default()
                }
            }
        }

        let style = Gerbil;

        let version = VERSION.lock().unwrap()[0].clone();
        let version_float = format!("1e-{}", -version.force_f64().log10());
        Modal::new(&mut self.modal_state, content, move |state| {
            Card::new(
                Text::new(""),
                Text::new(&format!(
                    "Welcome to enclone visual {} = {}!\n\n\
                     Please type bit.ly/enclone in a browser to learn more about enclone.\n\n\
                     To use enclone visual, type in the box \
                     (see below)\nand then push the Submit button.  Here are the things \
                     that you can type:\n\n\
                     • an enclone command, without the enclone part\n\
                     • an clonotype id (number)\n\
                     • d, for a demo, same as BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui\n\
                     • q to quit\n\n\
                     Major limitations of this version:\n\
                     1. There is no color in the clonotype tables.\n\
                     2. Text in plots does not show up.\n\
                     3. Cutting and pasting from clonotype tables doesn't work.\n\
                     4. Long commands are hard to work with in the input box.\n\
                     5. Very wide clonotype tables wrap, making them unintelligible, and \
                     only solvable by window resizing, and sometimes not that.",
                    version, version_float,
                ))
                .height(Units(450))
                .vertical_alignment(VerticalAlignment::Center),
            )
            .style(style)
            .foot(
                Row::new().spacing(10).push(
                    Button::new(
                        &mut state.cancel_state,
                        Text::new("Dismiss").horizontal_alignment(HorizontalAlignment::Left),
                    )
                    // .width(Length::Fill)
                    .on_press(Message::CancelButtonPressed),
                ),
            )
            .width(Units(1100))
            .height(Units(1060))
            .on_close(Message::CloseModal)
            .into()
        })
        .backdrop(Message::CloseModal)
        .on_esc(Message::CloseModal)
        .into()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

mod style {

    pub struct RuleStyle;

    impl iced::rule::StyleSheet for RuleStyle {
        fn style(&self) -> iced::rule::Style {
            iced::rule::Style {
                color: Color::from_rgb(0.0, 1.0, 1.0),
                width: 3,
                radius: 1.0,
                fill_mode: iced::rule::FillMode::Percent(100.0),
            }
        }
    }

    pub struct Squeak;

    use iced::{scrollable, Color};

    impl scrollable::StyleSheet for Squeak {
        fn active(&self) -> scrollable::Scrollbar {
            scrollable::Scrollbar {
                background: Color::from_rgb(0.75, 0.75, 0.75).into(),
                border_radius: 2.0,
                border_width: 0.0,
                border_color: Color::TRANSPARENT,
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    border_radius: 2.0,
                    border_width: 0.0,
                    border_color: Color::TRANSPARENT,
                },
            }
        }

        fn hovered(&self) -> scrollable::Scrollbar {
            let active = self.active();
            scrollable::Scrollbar {
                background: Color {
                    a: 0.5,
                    ..Color::from_rgb(0.0, 0.0, 0.0)
                }
                .into(),
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    ..active.scroller
                },
                ..active
            }
        }

        fn dragging(&self) -> scrollable::Scrollbar {
            let hovered = self.hovered();
            scrollable::Scrollbar {
                scroller: scrollable::Scroller {
                    color: Color::from_rgb(0.0, 0.0, 0.0),
                    ..hovered.scroller
                },
                ..hovered
            }
        }
    }
}
