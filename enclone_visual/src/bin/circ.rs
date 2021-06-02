use iced::{
    // canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
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

        /*
        let canvas = Canvas::new(self)
            .width(Length::Units(400))
            .height(Length::Units(400));

        Container::new(canvas)
            .width(Length::Fill)
            .height(Length::Fill)
            .padding(20)
            .center_x()
            .center_y()
            .into()
        */

        let content = Column::new()
            .push(button)
            .push(
                self.grid
                    .view()
                    .map(move |message| Message::ButtonPressed),
            );

        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            // .style(style::Container)
            .into()

    }
}

mod grid {
    use iced::{
        // canvas::event::{self, Event},
        canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
        Color, Element, Length, Rectangle,
    };

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
            // Self::from_preset(Preset::default())
        }
    }

    impl Grid {

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

        /*
        fn visible_region(&self, size: Size) -> Region {
            let width = size.width / self.scaling;
            let height = size.height / self.scaling;

            Region {
                x: -self.translation.x - width / 2.0,
                y: -self.translation.y - height / 2.0,
                width,
                height,
            }
        }
        */
    }

    impl<'a> canvas::Program<Message> for Grid {
        fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
            let mut frame = Frame::new(bounds.size());
            let radius = 100.0;
            let circle = Path::circle(frame.center(), radius);
            frame.fill(&circle, Color::BLACK);
            vec![frame.into_geometry()]
        }
    }
}




/*
impl canvas::Program<Message> for Clock {
    fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
        let mut frame = Frame::new(bounds.size());
        let radius = 100.0;
        let circle = Path::circle(frame.center(), radius);
        frame.fill(&circle, Color::BLACK);
        vec![frame.into_geometry()]
    }
}
*/
