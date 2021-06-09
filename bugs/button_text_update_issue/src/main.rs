// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    button, Application, Button, Clipboard, Column, Command, Container, Element, Length, Row,
    Settings, Subscription, Text,
};
use iced_native::{window, Event};

fn main() -> iced::Result {
    StrangeThing::run(Settings::default())
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
struct StrangeThing {
    button: button::State,
    should_exit: bool,
    compute_state: ComputeState,
}

#[derive(Debug, Clone)]
enum Message {
    ButtonPressed,
    ComputationDone(Result<(), String>),
    EventOccurred(iced_native::Event),
}

impl Application for StrangeThing {
    type Executor = iced::executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (StrangeThing, Command<Message>) {
        let mut x = StrangeThing::default();
        x.compute_state = WaitingForRequest;
        (x, Command::none())
    }

    fn title(&self) -> String {
        String::from("crazy")
    }

    fn update(&mut self, message: Message, _clipboard: &mut Clipboard) -> Command<Message> {
        match message {
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
            Message::EventOccurred(ref event) => {
                if let Event::Window(window::Event::CloseRequested) = event {
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
        let button = Button::new(
            &mut self.button,
            Text::new(if self.compute_state == WaitingForRequest {
                "Submit"
            } else {
                "thinking"
            }),
        )
        .on_press(Message::ButtonPressed);
        let content = Column::new().push(Row::new().spacing(10).push(button));
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .center_x()
            .center_y()
            .into()
    }
}

async fn compute() -> Result<(), String> {
    std::thread::sleep(std::time::Duration::from_millis(3000));
    Ok(())
}
