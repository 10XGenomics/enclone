// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    button, Button, Column, Container, Element, Length, Sandbox, Settings, Text,
};

use grid::Grid;

pub fn main() -> iced::Result {
    Clock::run(Settings::default())
}

#[derive(Default)]
struct Clock {
    grid: Grid,
    button: button::State,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    ButtonPressed,
}

impl Sandbox for Clock {
    type Message = Message;

    fn new() -> Self {
        Clock::default()
    }

    fn title(&self) -> String {
        String::from("Clock - Iced")
    }

    fn update(
        &mut self,
        message: Message,
    ) {
        match message {
            Message::ButtonPressed => {
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
                self.grid
                    .view()
                    .map(move |_message| Message::ButtonPressed),
            );
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

mod grid {
    use iced::{
        canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
        Color, Element, Length, Rectangle,
    };
    use std::time::{SystemTime, UNIX_EPOCH};

    #[derive(Default)]
    struct State {
    }

    pub struct Grid {
        state: State,
    }

    #[derive(Debug, Clone)]
    pub enum Message {
    }

    impl Default for Grid {
        fn default() -> Self {
            Self { 
                state: State::default(),
            }
        }
    }

    impl Grid {

        /*
        pub fn update(&mut self, message: Message) {
            match message {
            }
        }
        */

        pub fn view<'a>(&'a mut self) -> Element<'a, Message> {
            Canvas::new(self)
                .width(Length::Fill)
                .height(Length::Fill)
                .into()
        }
    }

    impl<'a> canvas::Program<Message> for Grid {
        fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
            let mut frame = Frame::new(bounds.size());
            let start = SystemTime::now();
            let since_the_epoch = start
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards");
            let nanos = since_the_epoch.subsec_nanos() as u64;
            let radius = (nanos % 99) as f32;
            let circle = Path::circle(frame.center(), radius);
            frame.fill(&circle, Color::BLACK);
            vec![frame.into_geometry()]
        }
    }
}
