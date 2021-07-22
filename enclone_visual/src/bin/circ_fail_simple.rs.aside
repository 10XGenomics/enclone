// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is also in enclone/bugs/canvas_in_scrollable, where we can use the latest version of iced.

use iced::canvas::{self, Canvas, Cursor, Frame, Geometry, Path};
use iced::{
    scrollable, Color, Column, Container, Element, Length, Rectangle, Sandbox, Scrollable, Settings,
};

pub fn main() -> iced::Result {
    Circles::run(Settings::default())
}

#[derive(Default)]
struct Circles {
    engine: Engine,
    scroll: scrollable::State,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    Something,
}

impl Sandbox for Circles {
    type Message = Message;

    fn new() -> Self {
        Circles::default()
    }

    fn title(&self) -> String {
        String::from("Circles")
    }

    fn update(&mut self, _message: Message) {}

    fn view(&mut self) -> Element<Message> {
        let engine = self.engine.view().map(move |_message| Message::Something);
        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            .push(engine);
        let content = Column::new().push(scrollable);
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

#[derive(Default)]
pub struct Engine {}

pub enum EngineMessage {}

impl Engine {
    pub fn view<'a>(&'a mut self) -> Element<'a, EngineMessage> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

impl<'a> canvas::Program<EngineMessage> for Engine {
    fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
        let mut frame = Frame::new(bounds.size());
        let circle = Path::circle(frame.center(), 20.0);
        frame.fill(&circle, Color::BLACK);
        vec![frame.into_geometry()]
    }
}
