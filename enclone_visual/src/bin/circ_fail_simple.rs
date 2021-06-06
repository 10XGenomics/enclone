// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use engine::Engine;
use iced::{
    button, Button, Column, Container, Element, Length, Sandbox, scrollable,
    Scrollable, Settings, Text
};

pub fn main() -> iced::Result {
    Circles::run(Settings::default())
}

#[derive(Default)]
struct Circles {
    engine: Engine,
    button: button::State,
    scroll: scrollable::State,
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

    fn update(&mut self, message: Message) {
        match message {
            Message::ButtonPressed => {
            }
        }
    }

    fn view(&mut self) -> Element<Message> {
        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);
        let engine = 
            self.engine
                .view()
                .map(move |_message| Message::ButtonPressed);
        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .push(engine);
        let content = Column::new().push(button).push(scrollable);
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

mod engine {
    use iced::{
        canvas::{self, Canvas, Cursor, Frame, Geometry, Path}, Color, Element, Length, Rectangle,
    };

    #[derive(Default)]
    pub struct State {
        pub radius: f32,
    }

    pub struct Engine {
        pub state: State,
    }

    #[derive(Debug, Clone)]
    pub enum Message {}

    impl Default for Engine {
        fn default() -> Self {
            Self {
                state: State {
                    radius: 10.0,
                }
            }
        }
    }

    impl Engine {
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
            let circle = Path::circle(frame.center(), 20.0);
            frame.fill(&circle, Color::from_rgb(0.5, 0.5, 1.0));
            vec![frame.into_geometry()]
        }
    }
}
