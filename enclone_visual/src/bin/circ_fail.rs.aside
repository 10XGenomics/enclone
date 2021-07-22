// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use engine::Engine;
use iced::{
    button, scrollable, Button, Column, Container, Element, Length, Sandbox, Scrollable, Settings,
    Text,
};
use std::time::{SystemTime, UNIX_EPOCH};

pub fn rotate(r: i64) -> i64 {
    6_364_136_223_846_793_005i64
        .wrapping_mul(r)
        .wrapping_add(1_442_695_040_888_963_407)
}

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
                self.engine.state.button_pressed = true;
                let start = SystemTime::now();
                let since_the_epoch = start
                    .duration_since(UNIX_EPOCH)
                    .expect("Time went backwards");
                let nanos = since_the_epoch.subsec_nanos() as u64;
                let radius = (nanos % 99) as f32;
                self.engine.state.radius = radius;
                self.engine.state.rand = rotate(self.engine.state.rand);
            }
        }
    }

    fn view(&mut self) -> Element<Message> {
        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);
        let engine = self
            .engine
            .view()
            .map(move |_message| Message::ButtonPressed);

        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(12)
            .scroller_width(12)
            // .style(style::Squeak)
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
    use crate::rotate;
    use iced::{
        canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
        Color, Element, Length, Rectangle,
    };
    use iced_native::Point;
    use iced_native::Vector;

    #[derive(Default)]
    pub struct State {
        pub button_pressed: bool,
        pub radius: f32,
        pub rand: i64,
    }

    pub struct Engine {
        pub state: State,
    }

    #[derive(Debug, Clone)]
    pub enum Message {}

    impl Default for Engine {
        fn default() -> Self {
            Self {
                state: State::default(),
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
        fn draw(&self, bounds: Rectangle, cursor: Cursor) -> Vec<Geometry> {
            let pos = cursor.position_in(&bounds);
            let mut frame = Frame::new(bounds.size());
            let radius = self.state.radius;
            let center = frame.center();
            if pos.is_some() {
                let xdiff = pos.unwrap().x - center.x;
                let ydiff = pos.unwrap().y - center.y;
                let dist = (xdiff * xdiff + ydiff * ydiff).sqrt();
                if dist <= radius {
                    frame.translate(Vector { x: 100.0, y: 100.0 });
                    frame.fill_text(format!(
                        "delta = ({:.1}, {:.1}); in circle one at distance {:.1} <= {:.1}",
                        pos.unwrap().x - center.x,
                        pos.unwrap().y - center.y,
                        dist,
                        radius
                    ));
                    frame.translate(Vector {
                        x: -100.0,
                        y: -100.0,
                    });
                }
            }
            let circle1 = Path::circle(center, radius);
            frame.fill(&circle1, Color::from_rgb(0.5, 0.5, 1.0));
            let circlex = Path::circle(center, 2.0);
            frame.fill(&circlex, Color::from_rgb(0.0, 0.0, 0.0));

            let mut r = self.state.rand;
            for _ in 0..10000 {
                r = rotate(r);
                let x = r % 200;
                r = rotate(r);
                let y = r % 200;
                r = rotate(r);
                let rad = r % 10;
                let circle2 = Path::circle(
                    Point {
                        x: center.x + x as f32,
                        y: center.y + y as f32,
                    },
                    rad as f32,
                );
                r = rotate(r);
                let c1 = 0.7 + (r % 1000) as f32 / 3000.0;
                r = rotate(r);
                let c2 = 0.7 + (r % 1000) as f32 / 3000.0;
                r = rotate(r);
                let c3 = 0.7 + (r % 1000) as f32 / 3000.0;
                frame.fill(&circle2, Color::from_rgb(c1, c2, c3));
            }

            vec![frame.into_geometry()]
        }
    }
}
