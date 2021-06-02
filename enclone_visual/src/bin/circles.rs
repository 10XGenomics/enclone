//! Circles
use iced::{
    canvas::{self, Cursor, Path, Stroke},
    executor, time, Application, Canvas, Clipboard, Color, Command,
    Element, Length, Point, Rectangle, Settings, Size, Subscription, Vector,
};

use std::time::Instant;

pub fn main() -> iced::Result {
    SolarSystem::run(Settings {
        antialiasing: true,
        ..Settings::default()
    })
}

struct SolarSystem {
    state: State,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    Tick(Instant),
}

impl Application for SolarSystem {
    type Executor = executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (Self, Command<Message>) {
        (
            SolarSystem {
                state: State::new(),
            },
            Command::none(),
        )
    }

    fn title(&self) -> String {
        String::from("Solar system - Iced")
    }

    fn update(
        &mut self,
        message: Message,
        _clipboard: &mut Clipboard,
    ) -> Command<Message> {
        match message {
            Message::Tick(instant) => {
                self.state.update(instant);
            }
        }

        Command::none()
    }

    fn subscription(&self) -> Subscription<Message> {
        time::every(std::time::Duration::from_millis(10))
            .map(|instant| Message::Tick(instant))
    }

    fn view(&mut self) -> Element<Message> {

        /*
        let button = Button::new(&mut self.button, Text::new("Submit"))
            .padding(10)
            .on_press(Message::ButtonPressed);

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .push(button)
            .push(
                Canvas::new(&mut self.state)
                    .width(Length::Fill)
                    .height(Length::Fill)
                    .into()
            );
        */




        Canvas::new(&mut self.state)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

#[derive(Debug)]
struct State {
    space_cache: canvas::Cache,
    system_cache: canvas::Cache,
    cursor_position: Point,
    start: Instant,
    now: Instant,
}

impl State {
    const SUN_RADIUS: f32 = 70.0;
    const ORBIT_RADIUS: f32 = 150.0;
    const EARTH_RADIUS: f32 = 12.0;
    const MOON_RADIUS: f32 = 4.0;
    const MOON_DISTANCE: f32 = 28.0;

    pub fn new() -> State {
        let now = Instant::now();
        // let (width, height) = window::Settings::default().size;

        State {
            space_cache: Default::default(),
            system_cache: Default::default(),
            cursor_position: Point::ORIGIN,
            start: now,
            now,
        }
    }

    pub fn update(&mut self, now: Instant) {
        self.now = now;
        self.system_cache.clear();
    }

}

impl<Message> canvas::Program<Message> for State {
    fn draw(
        &self,
        bounds: Rectangle,
        _cursor: Cursor,
    ) -> Vec<canvas::Geometry> {
        use std::f32::consts::PI;

        let background = self.space_cache.draw(bounds.size(), |frame| {
            let space = Path::rectangle(Point::new(0.0, 0.0), frame.size());
            frame.fill(&space, Color::BLACK);
            frame.translate(frame.center() - Point::ORIGIN);
        });

        let system = self.system_cache.draw(bounds.size(), |frame| {
            let center = frame.center();

            let sun = Path::circle(center, Self::SUN_RADIUS);
            let orbit = Path::circle(center, Self::ORBIT_RADIUS);

            frame.fill(&sun, Color::from_rgb8(0xF9, 0xD7, 0x1C));
            frame.stroke(
                &orbit,
                Stroke {
                    width: 1.0,
                    color: Color::from_rgba8(0, 153, 255, 0.1),
                    ..Stroke::default()
                },
            );

            let elapsed = self.now - self.start;
            let rotation = (2.0 * PI / 60.0) * elapsed.as_secs() as f32
                + (2.0 * PI / 60_000.0) * elapsed.subsec_millis() as f32;

            frame.with_save(|frame| {
                frame.translate(Vector::new(center.x, center.y));
                frame.rotate(rotation);
                frame.translate(Vector::new(Self::ORBIT_RADIUS, 0.0));

                let earth = Path::circle(Point::ORIGIN, Self::EARTH_RADIUS);
                let shadow = Path::rectangle(
                    Point::new(0.0, -Self::EARTH_RADIUS),
                    Size::new(
                        Self::EARTH_RADIUS * 4.0,
                        Self::EARTH_RADIUS * 2.0,
                    ),
                );

                frame.fill(&earth, Color::from_rgb8(0x6B, 0x93, 0xD6));

                frame.with_save(|frame| {
                    frame.rotate(rotation * 10.0);
                    frame.translate(Vector::new(0.0, Self::MOON_DISTANCE));

                    let moon = Path::circle(Point::ORIGIN, Self::MOON_RADIUS);
                    frame.fill(&moon, Color::WHITE);
                });

                frame.fill(
                    &shadow,
                    Color {
                        a: 0.7,
                        ..Color::BLACK
                    },
                );
            });
        });

        vec![background, system]
    }
}
