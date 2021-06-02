use iced::{
    canvas::{self, Canvas, Cursor, Frame, Geometry, Path},
    button, Button, Color, Container, Element, Length, Rectangle, Sandbox, Settings, Text,
};

pub fn main() -> iced::Result {
    Clock::run(Settings::default())
}

#[derive(Default)]
struct Clock {
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
    }
}

impl canvas::Program<Message> for Clock {
    fn draw(&self, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
        let mut frame = Frame::new(bounds.size());
        let radius = 100.0;
        let circle = Path::circle(frame.center(), radius);
        frame.fill(&circle, Color::BLACK);
        vec![frame.into_geometry()]
    }
}
