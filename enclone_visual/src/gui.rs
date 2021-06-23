// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::convert_svg_to_png::*;
use crate::copy_image_to_clipboard::*;
use crate::svg_to_geometry::*;
use crate::*;
use canvas_view::CanvasView;
use iced::svg::Handle;
use iced::Length::Units;
use iced::{
    button, scrollable, text_input, Align, Application, Button, Clipboard, Color, Column, Command,
    Element, Font, HorizontalAlignment, Image, Length, Row, Rule, Scrollable, Settings,
    Subscription, Svg, Text, TextInput, VerticalAlignment,
};
use iced_aw::{modal, Card, Modal};
use iced_native::{window, Event};
use perf_stats::*;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant};
use string_utils::*;

const DEJAVU_BOLD: Font = Font::External {
    name: "DEJAVU_BOLD",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono-Bold.ttf"),
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn launch_gui() -> iced::Result {
    let mut settings = Settings::default();
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = (1100 as u32, 1060 as u32); // reasonable minimum size
    settings.window = window_settings;
    settings.exit_on_close_request = false;
    EncloneVisual::run(settings)
}

#[derive(PartialEq)]
enum ComputeState {
    WaitingForRequest,
    Thinking,
}

impl Default for ComputeState {
    fn default() -> ComputeState {
        WaitingForRequest
    }
}

use ComputeState::*;

#[derive(Default)]
struct EncloneVisual {
    scroll: scrollable::State,
    input: text_input::State,
    input_value: String,
    output_value: String,
    svg_value: String,
    png_value: Vec<u8>,
    button: button::State,
    submit_button_text: String,
    open_state: button::State,
    modal_state: modal::State<ModalState>,
    should_exit: bool,
    compute_state: ComputeState,
    copy_button: button::State,
    copy_button_color: Color,
    canvas_view: CanvasView,
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
    OpenModal,
    CloseModal,
    CancelButtonPressed,
    ComputationDone(Result<(), String>),
    EventOccurred(iced_native::Event),
    CopyButtonPressed,
    CopyButtonFlashed(Result<(), String>),
    // CopyButtonFlashed,
}

#[derive(Default)]
struct ModalState {
    cancel_state: button::State,
}

impl Application for EncloneVisual {
    type Executor = iced::executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (EncloneVisual, Command<Message>) {
        let mut x = EncloneVisual::default();
        x.submit_button_text = "Submit".to_string();
        x.compute_state = WaitingForRequest;
        x.copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        (x, Command::none())
    }

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message, _clipboard: &mut Clipboard) -> Command<Message> {
        match message {
            Message::OpenModal => {
                self.modal_state.show(true);
                Command::none()
            }
            Message::CloseModal => {
                self.modal_state.show(false);
                Command::none()
            }
            Message::CancelButtonPressed => {
                self.modal_state.show(false);
                Command::none()
            }
            Message::InputChanged(ref value) => {
                self.input_value = value.to_string();
                Command::none()
            }
            Message::ButtonPressed => {
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    // The following sleep is needed to get the button text to consistenly update.
                    thread::sleep(Duration::from_millis(20));
                    USER_REQUEST.lock().unwrap().clear();
                    USER_REQUEST.lock().unwrap().push(self.input_value.clone());
                    PROCESSING_REQUEST.store(true, SeqCst);
                    Command::perform(compute(), Message::ComputationDone)
                } else {
                    Command::none()
                }
            }
            Message::ComputationDone(_) => {
                let mut reply_text = SERVER_REPLY_TEXT.lock().unwrap()[0].clone();
                if reply_text.contains("enclone failed") {
                    reply_text = format!("enclone failed{}", reply_text.after("enclone failed"));
                }
                reply_text += "\n \n \n"; // papering over truncation bug
                let mut reply_svg = String::new();
                if SERVER_REPLY_SVG.lock().unwrap().len() > 0 {
                    reply_svg = SERVER_REPLY_SVG.lock().unwrap()[0].clone();
                }
                self.output_value = reply_text.to_string();
                self.svg_value = reply_svg.to_string();
                if self.svg_value.len() > 0 {
                    self.png_value = convert_svg_to_png(&reply_svg.as_bytes());
                    let geometry = svg_to_geometry(&reply_svg, false);
                    if geometry.is_some() {
                        let mut ok = true;
                        for i in 0..geometry.as_ref().unwrap().len() {
                            match &geometry.as_ref().unwrap()[i] {
                                crate::geometry::Geometry::Text(ttt) => {
                                    if ttt.rotate != [0.0; 3] {
                                        ok = false;
                                    }
                                }
                                _ => {}
                            }
                        }
                        if ok {
                            self.canvas_view.state.geometry_value = geometry;
                        } else {
                            self.canvas_view.state.geometry_value = None;
                        }
                    } else {
                        if VERBOSE.load(SeqCst) {
                            println!("translation from svg to geometries failed");
                        }
                        self.canvas_view.state.geometry_value = None;
                    }
                }
                self.compute_state = WaitingForRequest;
                Command::none()
            }

            // Catch exit (when the upper left red button is pushed) and store DONE to make
            // the server thread exit gracefully.  Otherwise you will get a an error message
            // and a traceback.
            Message::EventOccurred(ref event) => {
                if let Event::Window(window::Event::CloseRequested) = event {
                    DONE.store(true, SeqCst);
                    thread::sleep(Duration::from_millis(50));
                    self.should_exit = true;
                }
                Command::none()
            }

            Message::CopyButtonPressed => {
                self.copy_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_png_bytes_to_mac_clipboard(&self.png_value);
                Command::perform(flash_copy_button(), Message::CopyButtonFlashed)
            }

            Message::CopyButtonFlashed(_) => {
                self.copy_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }
        }
    }

    fn should_exit(&self) -> bool {
        self.should_exit
    }

    fn subscription(&self) -> Subscription<Message> {
        iced_native::subscription::events().map(Message::EventOccurred)
    }

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
        .on_press(Message::ButtonPressed);

        let copy_button = Button::new(
            &mut self.copy_button,
            Text::new("Copy").color(self.copy_button_color),
        )
        .padding(10)
        .on_press(Message::CopyButtonPressed);

        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .style(style::Squeak)
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

        let png = include_bytes!("../../img/enclone_banner.png").to_vec();
        let banner = Image::new(iced::image::Handle::from_memory(png)).width(Units(500));

        let svg_as_png = Image::new(iced::image::Handle::from_memory(self.png_value.clone()))
            .height(Units(SVG_HEIGHT));

        let mut graphic_row = Row::new().spacing(10);
        if self.png_value.len() > 0 {
            if self.canvas_view.state.geometry_value.is_some() {
                graphic_row = graphic_row
                    .push(
                        self.canvas_view
                            .view()
                            .map(move |_message| Message::ButtonPressed),
                    )
                    .height(Units(SVG_HEIGHT));
            } else {
                graphic_row = graphic_row.push(svg_as_png);
            }
            graphic_row = graphic_row.push(copy_button);
        }

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
            // .push(Row::new().spacing(10).push(svg))
            .push(graphic_row)
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
                     • an enclone command\n\
                     • an group id (number)\n\
                     • d, for a demo, same as enclone BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui\n\
                     • q to quit\n\n\
                     Major limitations of this version:\n\
                     1. There is no color in the clonotype tables.\n\
                     2. Cutting and pasting from clonotype tables doesn't work.\n\
                     3. Long commands are hard to work with in the input box.\n\
                     4. Very wide clonotype tables wrap, making them unintelligible, and \
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

async fn compute() -> Result<(), String> {
    let t = Instant::now();
    while PROCESSING_REQUEST.load(SeqCst) {
        thread::sleep(Duration::from_millis(10));
    }
    println!(
        "time used processing command = {:.1} seconds\n",
        elapsed(&t)
    );
    Ok(())
}

async fn flash_copy_button() -> Result<(), String> {
    thread::sleep(Duration::from_millis(400));
    Ok(())
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
