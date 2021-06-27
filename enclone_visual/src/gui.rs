// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::convert_svg_to_png::*;
use crate::copy_image_to_clipboard::*;
use crate::svg_to_geometry::*;
use crate::*;
use gui_structures::*;
use gui_structures::ComputeState::*;
use iced::svg::Handle;
use iced::Length::Units;
use iced::{
    Align, Application, Button, Clipboard, Color, Column, Command,
    Element, Font, HorizontalAlignment, Image, Length, Row, Rule, Scrollable, Settings, Space,
    Svg, Text, TextInput, VerticalAlignment,
};
// use iced::Subscription;
use iced_aw::{Card, Modal};
// use iced_native::{window, Event};
use messages::Message;
use perf_stats::*;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::{Duration, Instant};

const DEJAVU_BOLD: Font = Font::External {
    name: "DEJAVU_BOLD",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono-Bold.ttf"),
};

pub static COUNT: AtomicUsize = AtomicUsize::new(0);

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub async fn launch_gui() -> iced::Result {
    let mut settings = Settings::default();
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = (1100 as u32, 1060 as u32); // reasonable minimum size
    settings.window = window_settings;
    settings.exit_on_close_request = false;
    EncloneVisual::run(settings)
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
        match message {

            // Identical to ButtonPressed, except for function signature.
            // This code duplication is horrible and should be fixed.

            Message::ButtonPressedX(_) => {
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

            Message::RunTests(_) => {
                // thread::sleep(Duration::from_millis(1000));
                self.input_value = "#1".to_string();
                self.view();
                let count = COUNT.load(SeqCst);
                if count == 0 {
                    self.window_id = get_window_id();
                    self.input_value = "#1".to_string();
                } else if count == 1 {
                    self.input_value = "#2".to_string();
                } else if count == 2 {
                    self.input_value = "#3".to_string();
                } else {
                    let count = count + 1;
                    capture(count, self.window_id);
                    std::process::exit(0);
                }
                COUNT.store(COUNT.load(SeqCst) + 1, SeqCst);
                Command::perform(noop(), Message::ButtonPressedX)
            }

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
                if !TEST_MODE.load(SeqCst) {
                    Command::none()
                } else {
                    let count = COUNT.load(SeqCst);
                    if count > 1 {
                        capture(count, self.window_id);
                    }
                    Command::perform(noop(), Message::RunTests)
                }
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






        // Define the button complex that is the "control panel".

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
            let cmd = &self.command_history[self.history_index - 1];
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

        // Create execute button.

        let exec_button = Button::new(
            &mut self.exec_button,
            Text::new("Execute command").size(COPY_BUTTON_FONT_SIZE),
        )
        .on_press(Message::ExecuteButtonPressed);

        // Build the command column.

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
        if self.history_index >= 1 && !self.is_blank[self.history_index - 1] {
            col = col.push(copy_image_button);
        } else {
            col = col.push(null_copy_image_button);
        }
        col = col.push(exec_button);





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

        let png = include_bytes!("../../img/enclone_banner.png").to_vec();
        let banner = Image::new(iced::image::Handle::from_memory(png)).width(Units(500));

        let svg_as_png = Image::new(iced::image::Handle::from_memory(self.png_value.clone()))
            .height(Units(SVG_HEIGHT));

        let have_canvas = self.canvas_view.state.geometry_value.is_some();
        let mut graphic_row = Row::new().spacing(10);
        if self.png_value.len() > 0 {
            // Show the graphic.

            if have_canvas {
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

            // Insert space to push the graphic to the left and the command column to the right.

            if !have_canvas {
                graphic_row = graphic_row.push(Space::with_width(Length::Fill));
            }

            let mut command_complex = Row::new().spacing(10);

            // Add the command column to the row.

            command_complex = command_complex.push(col);

            // Add up and down arrows.

            command_complex = command_complex.push(button_column2);

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
        if !self.command_history.is_empty() {
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
                if !COOKBOOK.load(SeqCst) {
                    Text::new(&format!(
                        "Welcome to enclone visual {} = {}!\n\n{}",
                        version, version_float, include_str!["help.txt"], 
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
                .height(Units(600))
                .vertical_alignment(VerticalAlignment::Center),
            )
            .style(style::Help)
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

async fn noop() -> Result<(), String> {
    // Increasing this time to 2000ms will prevent the screen from going dark on initialization
    // in test mode.
    thread::sleep(Duration::from_millis(100));
    Ok(())
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
