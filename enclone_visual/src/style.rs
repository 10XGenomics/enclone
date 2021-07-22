// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// This file defines styles for some elements in the GUI.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub struct RuleStyle;

impl iced::rule::StyleSheet for RuleStyle {
    fn style(&self) -> iced::rule::Style {
        iced::rule::Style {
            color: Color::from_rgb(0.0, 1.0, 1.0),
            width: 3,
            radius: 1.0,
            fill_mode: iced::rule::FillMode::Percent(100.0),
        }
    }
}

pub struct RuleStyle2;

impl iced::rule::StyleSheet for RuleStyle2 {
    fn style(&self) -> iced::rule::Style {
        iced::rule::Style {
            color: Color::from_rgb(0.8, 0.7, 0.7),
            width: 4,
            radius: 1.0,
            fill_mode: iced::rule::FillMode::Percent(100.0),
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub struct ScrollableStyle;

use iced::{scrollable, Color};

impl scrollable::StyleSheet for ScrollableStyle {
    fn active(&self) -> scrollable::Scrollbar {
        scrollable::Scrollbar {
            background: Color::from_rgb(0.75, 0.75, 0.75).into(),
            border_radius: 2.0,
            border_width: 0.0,
            border_color: Color::TRANSPARENT,
            scroller: scrollable::Scroller {
                color: Color::from_rgb(0.0, 0.0, 0.0),
                border_radius: 2.0,
                border_width: 0.0,
                border_color: Color::TRANSPARENT,
            },
        }
    }

    fn hovered(&self) -> scrollable::Scrollbar {
        let active = self.active();
        scrollable::Scrollbar {
            background: Color {
                a: 0.5,
                ..Color::from_rgb(0.0, 0.0, 0.0)
            }
            .into(),
            scroller: scrollable::Scroller {
                color: Color::from_rgb(0.0, 0.0, 0.0),
                ..active.scroller
            },
            ..active
        }
    }

    fn dragging(&self) -> scrollable::Scrollbar {
        let hovered = self.hovered();
        scrollable::Scrollbar {
            scroller: scrollable::Scroller {
                color: Color::from_rgb(0.0, 0.0, 0.0),
                ..hovered.scroller
            },
            ..hovered
        }
    }
}
