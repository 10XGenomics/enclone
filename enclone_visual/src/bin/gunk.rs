// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    button, text_input, Application, Button, Clipboard, Column, Command, 
    Container, Element, Length, Row, Settings, Subscription, Text, TextInput,
};
use iced_native::{window, Event};
use std::thread;
use std::time::Duration;

fn main() -> iced::Result {
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
    input: text_input::State,
    input_value: String,
    button: button::State,
    submit_button_text: String,
    should_exit: bool,
    compute_state: ComputeState,
    
}

#[derive(Debug, Clone)]
struct Gerbil {
}

#[derive(Debug, Clone)]
enum Message {
    InputChanged(String),
    ButtonPressed,
    ComputationDone(Result<Gerbil, String>),
    EventOccurred(iced_native::Event),
}

impl Application for EncloneVisual {
    type Executor = iced::executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (EncloneVisual, Command<Message>) {
        let mut x = EncloneVisual::default();
        x.submit_button_text = "Submit".to_string();
        x.compute_state = WaitingForRequest;
        (x, Command::none())
    }

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(
        &mut self, 
        message: Message,
       _clipboard: &mut Clipboard,
    ) -> Command<Message> {
        match message {
            Message::InputChanged(ref value) => {
                self.input_value = value.to_string();
                Command::none()
            }
            Message::ButtonPressed => {
                println!("pushed");
                if self.compute_state == WaitingForRequest {
                    self.compute_state = Thinking;
                    Command::perform(compute(), Message::ComputationDone)
                } else {
                    Command::none()
                }
            }
            Message::ComputationDone(_) => {
                self.compute_state = WaitingForRequest;
                Command::none()
            }

            // Catch exit (when the upper left red button is pushed) and store DONE to make
            // the server thread exit gracefully.  Otherwise you will get a an error message
            // and a traceback.

            Message::EventOccurred(ref event) => {
                if let Event::Window(window::Event::CloseRequested) = event {
                    thread::sleep(Duration::from_millis(50));
                    self.should_exit = true;
                }
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
        .size(14);

        let button = Button::new(&mut self.button, 
            Text::new(if self.compute_state == WaitingForRequest { "Submit" } else { "thinking" }))
            .padding(10)
            .on_press(Message::ButtonPressed);

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .max_width(1500)
            .push(Row::new().spacing(10).push(text_input).push(button));

        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .center_x()
            .center_y()
            .into()

    }
}

async fn compute() -> Result<Gerbil, String> {
    thread::sleep(Duration::from_millis(2000));
    Ok(Gerbil{})
}

/*
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
}
*/
