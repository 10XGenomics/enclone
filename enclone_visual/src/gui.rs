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
    Element, Font, HorizontalAlignment, Image, Length, Row, Rule, Scrollable, Settings, Svg, Text,
    TextInput, VerticalAlignment,
};
// use iced::Subscription;
use iced_aw::{modal, Card, Modal};
// use iced_native::{window, Event};
use perf_stats::*;
use std::collections::HashMap;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant};

const DEJAVU_BOLD: Font = Font::External {
    name: "DEJAVU_BOLD",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono-Bold.ttf"),
};

fn blank_svg() -> String {
    r###"<svg version="1.1" baseProfile="full" width="400" height="400"
xmlns="http://www.w3.org/2000/svg">
<rect x="0" y="0" width="400" height="400" style="fill:white" />
</svg>
"###
    .to_string()
}

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
    translated_input_value: String,
    output_value: String,
    svg_value: String,
    png_value: Vec<u8>,
    button: button::State,
    back_button: button::State,
    forward_button: button::State,
    exec_button: button::State,
    submit_button_text: String,
    open_state: button::State,
    open_state_cookbook: button::State,
    exit_state: button::State,
    modal_state_help: modal::State<ModalState>,
    // should_exit: bool,
    compute_state: ComputeState,
    copy_image_button: button::State,
    copy_image_button_color: Color,
    canvas_view: CanvasView,
    command_copy_button: button::State,
    null_button1: button::State,
    null_button2: button::State,
    null_button3: button::State,
    null_button: button::State,
    cookbook: HashMap<String, String>,
    clear_button: button::State,

    // parallel vectors:
    svg_history: Vec<String>,
    command_history: Vec<String>,
    is_blank: Vec<bool>,

    // index of "current" position in those vectors:
    history_index: usize,
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
    BackButtonPressed,
    ForwardButtonPressed,
    ExecuteButtonPressed,
    OpenModalHelp,
    CloseModalHelp,
    OpenModalCookbook,
    CancelButtonPressed,
    ComputationDone(Result<(), String>),
    // EventOccurred(iced_native::Event),
    GraphicsCopyButtonPressed,
    GraphicsCopyButtonFlashed(Result<(), String>),
    CommandCopyButtonPressed,
    DoNothing,
    Exit,
    ClearButtonPressed,
}

#[derive(Default)]
pub struct ModalState {
    cancel_state: button::State,
    pub cookbook: bool,
}

impl EncloneVisual {
    pub fn post_svg(&mut self, svg: &str) {
        self.png_value = convert_svg_to_png(&svg.as_bytes());
        let geometry = svg_to_geometry(&svg, false);
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
}

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
        (x, Command::none())
    }

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message, _clipboard: &mut Clipboard) -> Command<Message> {
        match message {
            Message::OpenModalHelp => {
                COOKBOOK.store(false, SeqCst);
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::OpenModalCookbook => {
                COOKBOOK.store(true, SeqCst);
                self.modal_state_help.show(true);
                Command::none()
            }

            Message::CloseModalHelp => {
                self.modal_state_help.show(false);
                Command::none()
            }

            Message::Exit => {
                if true {
                    std::process::exit(0);
                }
                Command::none()
            }

            Message::CancelButtonPressed => {
                self.modal_state_help.show(false);
                Command::none()
            }

            Message::InputChanged(ref value) => {
                self.input_value = value.to_string();
                Command::none()
            }

            Message::ClearButtonPressed => {
                self.input_value.clear();
                Command::none()
            }

            Message::ButtonPressed => {
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    // The following sleep is needed to get the button text to consistenly update.
                    thread::sleep(Duration::from_millis(20));
                    if self.input_value.starts_with('#')
                        && self.cookbook.contains_key(&self.input_value)
                    {
                        self.translated_input_value = self.cookbook[&self.input_value].clone();
                    } else {
                        self.translated_input_value = self.input_value.clone();
                    }
                    USER_REQUEST.lock().unwrap().clear();
                    USER_REQUEST
                        .lock()
                        .unwrap()
                        .push(self.translated_input_value.clone());
                    PROCESSING_REQUEST.store(true, SeqCst);
                    Command::perform(compute(), Message::ComputationDone)
                } else {
                    Command::none()
                }
            }

            // same as above except for first line
            Message::ExecuteButtonPressed => {
                self.input_value = self.command_history[self.history_index - 1].clone();
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    // The following sleep is needed to get the button text to consistenly update.
                    thread::sleep(Duration::from_millis(20));
                    if self.input_value.starts_with('#')
                        && self.cookbook.contains_key(&self.input_value)
                    {
                        self.translated_input_value = self.cookbook[&self.input_value].clone();
                    } else {
                        self.translated_input_value = self.input_value.clone();
                    }
                    USER_REQUEST.lock().unwrap().clear();
                    USER_REQUEST
                        .lock()
                        .unwrap()
                        .push(self.translated_input_value.clone());
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
                    let mut blank = false;
                    if reply_svg.len() == 0 {
                        reply_svg = blank_svg();
                        blank = true;
                    }
                    if reply_svg.len() > 0 && self.input_value.parse::<usize>().is_err() {
                        self.svg_history.push(reply_svg.clone());
                        self.history_index += 1;
                        self.command_history
                            .push(self.translated_input_value.clone());
                        self.is_blank.push(blank);
                    }
                }
                self.output_value = reply_text.to_string();
                self.svg_value = reply_svg.to_string();
                if self.svg_value.len() > 0 {
                    self.post_svg(&reply_svg);
                }
                self.compute_state = WaitingForRequest;
                Command::none()
            }

            // Catch exit (when the upper left red button is pushed) and store DONE to make
            // the server thread exit gracefully.  Otherwise you will get a an error message
            // and a traceback.
            /*
            Message::EventOccurred(ref event) => {
                if let Event::Window(window::Event::CloseRequested) = event {
                    DONE.store(true, SeqCst);
                    thread::sleep(Duration::from_millis(50));
                    self.should_exit = true;
                }
                Command::none()
            }
            */
            Message::GraphicsCopyButtonPressed => {
                self.copy_image_button_color = Color::from_rgb(1.0, 0.0, 0.0);
                copy_png_bytes_to_mac_clipboard(&self.png_value);
                Command::perform(
                    flash_copy_image_button(),
                    Message::GraphicsCopyButtonFlashed,
                )
            }

            Message::GraphicsCopyButtonFlashed(_) => {
                self.copy_image_button_color = Color::from_rgb(0.0, 0.0, 0.0);
                Command::none()
            }

            Message::BackButtonPressed => {
                self.history_index -= 1;
                let x = self.svg_history[self.history_index - 1].clone();
                self.post_svg(&x);
                Command::none()
            }

            Message::ForwardButtonPressed => {
                self.history_index += 1;
                let x = self.svg_history[self.history_index - 1].clone();
                self.post_svg(&x);
                Command::none()
            }

            Message::CommandCopyButtonPressed => {
                copy_bytes_to_mac_clipboard(
                    &self.command_history[self.history_index - 1].as_bytes(),
                );
                Command::none()
            }

            Message::DoNothing => Command::none(),
        }
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
        .on_press(Message::ButtonPressed);

        let clear_button = Button::new(&mut self.clear_button, Text::new("Clear"))
            .padding(10)
            .on_press(Message::ClearButtonPressed);

        const FB_BUTTON_FONT_SIZE: u16 = 45;
        let back_button = Button::new(
            &mut self.back_button,
            Text::new("⇧").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
        )
        .on_press(Message::BackButtonPressed);

        let forward_button = Button::new(
            &mut self.forward_button,
            Text::new("⇩").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
        )
        .on_press(Message::ForwardButtonPressed);

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

        let mut graphic_row = Row::new().spacing(10);
        if self.png_value.len() > 0 {
            // Show the graphic.

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

            // Add button column.

            // graphic_row = graphic_row.push(button_column);

            // Add command box.

            const MAX_LINE: usize = 45;
            let cmd = &self.command_history[self.history_index - 1];
            let mut rows = Vec::<Vec<String>>::new();
            {
                let words = cmd.split(' ').collect::<Vec<&str>>();
                let mut current = String::new();
                let mut i = 0;
                while i < words.len() {
                    if current.len() > 0 && current.len() + 1 + words[i].len() > MAX_LINE {
                        rows.push(vec![current.clone()]);
                        current.clear();
                        i -= 1;
                    } else if words[i].len() >= MAX_LINE {
                        let mut w = words[i].as_bytes().to_vec();
                        loop {
                            let n = std::cmp::min(MAX_LINE, w.len());
                            let sub = stringme(&w[0..n]);
                            if n < w.len() {
                                rows.push(vec![sub]);
                                w = w[n..w.len()].to_vec();
                            } else {
                                current = stringme(&w);
                                break;
                            }
                        }
                    } else if current.len() == 0 {
                        current += &mut words[i].clone();
                    } else {
                        current += &mut format!(" {}", words[i]);
                    }
                    i += 1;
                }
                if current.len() > 0 {
                    rows.push(vec![current]);
                }
            }
            let mut log = String::new();
            for i in 0..rows.len() {
                if i > 0 {
                    log += "\n";
                }
                log += &mut rows[i][0].clone();
            }

            let exec_button = Button::new(
                &mut self.exec_button,
                Text::new("Execute command").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::ExecuteButtonPressed);

            let mut col = Column::new()
                .spacing(8)
                .align_items(Align::End)
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
            if !self.is_blank[self.history_index - 1] {
                col = col.push(copy_image_button);
            } else {
                col = col.push(null_copy_image_button);
            }
            col = col.push(exec_button);

            graphic_row = graphic_row.push(col);

            // Add up and down arrows.

            graphic_row = graphic_row.push(button_column2);
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
        if !self.command_history.is_empty() {
            content = content.push(Rule::horizontal(10).style(style::RuleStyle));
        }
        content = content.push(
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
        Modal::new(&mut self.modal_state_help, content, move |state| {
            Card::new(
                Text::new(""),
                if !COOKBOOK.load(SeqCst) {
                    Text::new(&format!(
                        "Welcome to enclone visual {} = {}!\n\n\
                         Please type bit.ly/enclone in a browser to learn more about enclone.\n\n\
                         To use enclone visual, type in the box \
                         (see below)\nand then push the Submit button.  Here are the things \
                         that you can type:\n\n\
                         • an enclone command\n\
                         • an group id (number)\n\
                         • a recipe tag (see cookbook)\n\
                         • q to quit\n\n\
                         Some limitations of this version:\n\
                         1. There is no color in the clonotype tables.\n\
                         2. Cutting and pasting from clonotype tables doesn't work.\n\
                         3. Long commands are hard to work with in the input box.\n\
                         4. Very wide clonotype tables wrap, making them unintelligible, and \
                         only solvable by window resizing, and sometimes not that.",
                        version, version_float,
                    ))
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
            .on_close(Message::CloseModalHelp)
            .into()
        })
        .backdrop(Message::CloseModalHelp)
        .on_esc(Message::CloseModalHelp)
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

async fn flash_copy_image_button() -> Result<(), String> {
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
