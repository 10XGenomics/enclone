// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    button, Button, Column, Container, Element, Length, Sandbox, Settings, Text,
};
use engine::Engine;
use std::time::{SystemTime, UNIX_EPOCH};

pub fn main() -> iced::Result {
    Circles::run(Settings::default())
}

#[derive(Default)]
struct Circles {
    engine: Engine,
    button: button::State,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    ButtonPressed,
}

impl Sandbox for Circles {
    type Message = Message;

    fn new() -> Self {
        Circles::default()
    }

    fn title(&self) -> String {
        String::from("Circles")
    }

    fn update(
        &mut self,
        message: Message,
    ) {
        match message {
            Message::ButtonPressed => {
                self.engine.state.button_pressed = true;
                let start = SystemTime::now();
                let since_the_epoch = start
                    .duration_since(UNIX_EPOCH)
                    .expect("Time went backwards");
                let nanos = since_the_epoch.subsec_nanos() as u64;
                let radius = (nanos % 99) as f32;
                self.engine.state.radius = radius;
            }
        }
    }

    fn view(&mut self) -> Element<Message> {
        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);
        let content = Column::new()
            .push(button)
            .push(
                self.engine
                    .view()
                    .map(move |_message| Message::ButtonPressed),
            );
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

mod engine {
    use iced::{
        canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
        Color, Element, Length, Rectangle,
    };
    use std::time::{SystemTime, UNIX_EPOCH};

    #[derive(Default)]
    pub struct State {
        pub button_pressed: bool,
        pub radius: f32,
    }

    pub struct Engine {
        pub state: State,
        pub last_radius: f32,
    }

    #[derive(Debug, Clone)]
    pub enum Message {
    }

    impl Default for Engine {
        fn default() -> Self {
            Self { 
                state: State::default(),
                last_radius: 0.0,
            }
        }
    }

    impl Engine {

        pub fn update(&mut self, message: Message) {
            match message {
            }
        }

        pub fn view<'a>(&'a mut self) -> Element<'a, Message> {
            Canvas::new(self)
                .width(Length::Fill)
                .height(Length::Fill)
                .into()
        }
    }

    impl<'a> canvas::Program<Message> for Engine {
        fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
            let mut frame = Frame::new(bounds.size());
            /*
            let start = SystemTime::now();
            let since_the_epoch = start
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards");
            let nanos = since_the_epoch.subsec_nanos() as u64;
            let radius = (nanos % 99) as f32;
            */
            let radius = self.state.radius;
            // self.last_radius = radius;
            let circle = Path::circle(frame.center(), radius);
            frame.fill(&circle, Color::BLACK);
            vec![frame.into_geometry()]
        }
    }
}
