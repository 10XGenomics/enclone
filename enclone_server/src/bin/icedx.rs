// mod style;

use iced::Font;

use iced::{
    button, scrollable, Button, Column, Container, Element, Length, Radio, Row, Rule, Sandbox,
    Scrollable, Settings, Space, Text,
};

mod style {

    use iced::{container, radio, rule, scrollable};

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Theme {
        Light,
        Dark,
    }

    impl Theme {
        pub const ALL: [Theme; 2] = [Theme::Light, Theme::Dark];
    }

    impl Default for Theme {
        fn default() -> Theme {
            Theme::Light
        }
    }

    impl From<Theme> for Box<dyn container::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default(),
                Theme::Dark => dark::Container.into(),
            }
        }
    }

    impl From<Theme> for Box<dyn radio::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default(),
                Theme::Dark => dark::Radio.into(),
            }
        }
    }

    impl From<Theme> for Box<dyn scrollable::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default(),
                Theme::Dark => dark::Scrollable.into(),
            }
        }
    }

    impl From<Theme> for Box<dyn rule::StyleSheet> {
        fn from(theme: Theme) -> Self {
            match theme {
                Theme::Light => Default::default(),
                Theme::Dark => dark::Rule.into(),
            }
        }
    }

    mod dark {
        use iced::{container, radio, rule, scrollable, Color};

        const BACKGROUND: Color = Color::from_rgb(
            0x36 as f32 / 255.0,
            0x39 as f32 / 255.0,
            0x3F as f32 / 255.0,
        );

        const SURFACE: Color = Color::from_rgb(
            0x40 as f32 / 255.0,
            0x44 as f32 / 255.0,
            0x4B as f32 / 255.0,
        );

        const ACCENT: Color = Color::from_rgb(
            0x6F as f32 / 255.0,
            0xFF as f32 / 255.0,
            0xE9 as f32 / 255.0,
        );

        const ACTIVE: Color = Color::from_rgb(
            0x72 as f32 / 255.0,
            0x89 as f32 / 255.0,
            0xDA as f32 / 255.0,
        );

        const SCROLLBAR: Color = Color::from_rgb(
            0x2E as f32 / 255.0,
            0x33 as f32 / 255.0,
            0x38 as f32 / 255.0,
        );

        const SCROLLER: Color = Color::from_rgb(
            0x20 as f32 / 255.0,
            0x22 as f32 / 255.0,
            0x25 as f32 / 255.0,
        );

        pub struct Container;

        impl container::StyleSheet for Container {
            fn style(&self) -> container::Style {
                container::Style {
                    background: Color {
                        a: 0.99,
                        ..BACKGROUND
                    }
                    .into(),
                    text_color: Color::WHITE.into(),
                    ..container::Style::default()
                }
            }
        }

        pub struct Radio;

        impl radio::StyleSheet for Radio {
            fn active(&self) -> radio::Style {
                radio::Style {
                    background: SURFACE.into(),
                    dot_color: ACTIVE,
                    border_width: 1.0,
                    border_color: ACTIVE,
                }
            }

            fn hovered(&self) -> radio::Style {
                radio::Style {
                    background: Color { a: 0.5, ..SURFACE }.into(),
                    ..self.active()
                }
            }
        }

        pub struct Scrollable;

        impl scrollable::StyleSheet for Scrollable {
            fn active(&self) -> scrollable::Scrollbar {
                scrollable::Scrollbar {
                    background: Color {
                        a: 0.8,
                        ..SCROLLBAR
                    }
                    .into(),
                    border_radius: 2.0,
                    border_width: 0.0,
                    border_color: Color::TRANSPARENT,
                    scroller: scrollable::Scroller {
                        color: Color { a: 0.7, ..SCROLLER },
                        border_radius: 2.0,
                        border_width: 0.0,
                        border_color: Color::TRANSPARENT,
                    },
                }
            }

            fn hovered(&self) -> scrollable::Scrollbar {
                let active = self.active();

                scrollable::Scrollbar {
                    background: SCROLLBAR.into(),
                    scroller: scrollable::Scroller {
                        color: SCROLLER,
                        ..active.scroller
                    },
                    ..active
                }
            }

            fn dragging(&self) -> scrollable::Scrollbar {
                let hovered = self.hovered();

                scrollable::Scrollbar {
                    scroller: scrollable::Scroller {
                        color: ACCENT,
                        ..hovered.scroller
                    },
                    ..hovered
                }
            }
        }

        pub struct Rule;

        impl rule::StyleSheet for Rule {
            fn style(&self) -> rule::Style {
                rule::Style {
                    color: SURFACE,
                    width: 2,
                    radius: 1.0,
                    fill_mode: rule::FillMode::Percent(30.0),
                }
            }
        }
    }
}

const CQ_MONO: Font = Font::External {
    name: "CQ_MONO",
    bytes: include_bytes!("../../../fonts/DejaVuLGCSansMono.ttf"),
};

pub fn main() -> iced::Result {
    ScrollableDemo::run(Settings::default())
}

struct ScrollableDemo {
    theme: style::Theme,
    variants: Vec<Variant>,
}

#[derive(Debug, Clone)]
enum Message {
    ThemeChanged(style::Theme),
}

impl Sandbox for ScrollableDemo {
    type Message = Message;

    fn new() -> Self {
        ScrollableDemo {
            theme: Default::default(),
            variants: Variant::all(),
        }
    }

    fn title(&self) -> String {
        String::from("Scrollable - Iced")
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::ThemeChanged(theme) => self.theme = theme,
        }
    }

    fn view(&mut self) -> Element<Message> {
        let ScrollableDemo {
            theme, variants, ..
        } = self;

        let choose_theme = style::Theme::ALL.iter().fold(
            Column::new().spacing(10).push(Text::new("Choose a theme:")),
            |column, option| {
                column.push(
                    Radio::new(
                        *option,
                        &format!("{:?}", option),
                        Some(*theme),
                        Message::ThemeChanged,
                    )
                    .style(*theme),
                )
            },
        );

        let scrollable_row = Row::with_children(
            variants
                .iter_mut()
                .map(|variant| {
                    let mut scrollable = Scrollable::new(&mut variant.scrollable)
                        .padding(10)
                        .spacing(10)
                        .width(Length::Fill)
                        .height(Length::Fill)
                        .style(*theme)
                        .push(Text::new(variant.title));

                    if let Some(scrollbar_width) = variant.scrollbar_width {
                        scrollable = scrollable
                            .scrollbar_width(scrollbar_width)
                            .push(Text::new(format!("scrollbar_width: {:?}", scrollbar_width)));
                    }

                    if let Some(scrollbar_margin) = variant.scrollbar_margin {
                        scrollable = scrollable
                            .scrollbar_margin(scrollbar_margin)
                            .push(Text::new(format!(
                                "scrollbar_margin: {:?}",
                                scrollbar_margin
                            )));
                    }

                    if let Some(scroller_width) = variant.scroller_width {
                        scrollable = scrollable
                            .scroller_width(scroller_width)
                            .push(Text::new(format!("scroller_width: {:?}", scroller_width)));
                    }

                    scrollable = scrollable
                        .push(Space::with_height(Length::Units(100)))
                        .push(
                            Text::new(
                                "Some content that should wrap within the \
                            scrollable. Let's output a lot of short words, so \
                            that we'll make sure to see how wrapping works \
                            with these scrollbars.\n\n\
\n\
┌───────────┬────────────────────────────────────────────┬──────────────────────────┐\n\
│           │  CHAIN 1                                   │  CHAIN 2                 │\n\
│           │  144.1.2|IGHV3-49 ◆ 737|IGHJ6              │  265|IGKV2-28 ◆ 213|IGKJ1│\n\
│           ├────────────────────────────────────────────┼──────────────────────────┤\n\
│           │     11111111111111111111111111             │  11111111111             │\n\
│           │  35 11112222222222333333333344             │  11111111222             │\n\
│           │  15 67890123456789012345678901             │  23456789012             │\n\
│           │     ═══════════CDR3═══════════             │  ════CDR3═══             │\n\
│reference  │  QV ◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦GMDVW             │  CMQ◦◦◦◦◦◦◦◦             │\n\
│donor ref  │  QF ◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦GMDVW             │  CMQ◦◦◦◦◦◦◦◦             │\n\
├───────────┼────────────────────────────────────────────┼──────────────────────────┤\n\
│#  n       │  .. ..........................   u  const  │  ...........   u  const  │\n\
│1  1       │  KF CTRSSTTPRDPTMIVVAYYYYGMDVW  17  IGHG1  │  CMQALQTSWTF  25  IGKC   │\n\
└───────────┴────────────────────────────────────────────┴──────────────────────────┘\n",
                            )
                            .font(CQ_MONO),
                        )
                        .push(Space::with_height(Length::Units(1200)))
                        .push(Text::new("Middle"))
                        .push(Space::with_height(Length::Units(1200)))
                        .push(
                            Button::new(&mut variant.button, Text::new("I am a button"))
                                .width(Length::Fill)
                                .padding(10),
                        )
                        .push(Text::new("The End."));

                    Container::new(scrollable)
                        .width(Length::Fill)
                        .height(Length::Fill)
                        .style(*theme)
                        .into()
                })
                .collect(),
        )
        .spacing(20)
        .width(Length::Fill)
        .height(Length::Fill);

        let content = Column::new()
            .spacing(20)
            .padding(20)
            .push(choose_theme)
            .push(Rule::horizontal(20).style(self.theme))
            .push(scrollable_row);

        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .center_x()
            .center_y()
            .style(self.theme)
            .into()
    }
}

/// A version of a scrollable
struct Variant {
    title: &'static str,
    scrollable: scrollable::State,
    button: button::State,
    scrollbar_width: Option<u16>,
    scrollbar_margin: Option<u16>,
    scroller_width: Option<u16>,
}

impl Variant {
    pub fn all() -> Vec<Self> {
        vec![
            Self {
                title: "Default Scrollbar",
                scrollable: scrollable::State::new(),
                button: button::State::new(),
                scrollbar_width: None,
                scrollbar_margin: None,
                scroller_width: None,
            },
            /*
            Self {
                title: "Slimmed & Margin",
                scrollable: scrollable::State::new(),
                button: button::State::new(),
                scrollbar_width: Some(4),
                scrollbar_margin: Some(3),
                scroller_width: Some(4),
            },
            Self {
                title: "Wide Scroller",
                scrollable: scrollable::State::new(),
                button: button::State::new(),
                scrollbar_width: Some(4),
                scrollbar_margin: None,
                scroller_width: Some(10),
            },
            Self {
                title: "Narrow Scroller",
                scrollable: scrollable::State::new(),
                button: button::State::new(),
                scrollbar_width: Some(10),
                scrollbar_margin: None,
                scroller_width: Some(4),
            },
            */
        ]
    }
}
