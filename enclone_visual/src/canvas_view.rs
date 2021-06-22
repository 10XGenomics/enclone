// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    canvas::{self, Canvas, Cursor, Frame, Geometry, Path, Stroke},
    Color, Element, Length, Rectangle, Size,
};
use iced_native::{Point, Vector};

#[derive(Default)]
pub struct State {
    pub geometry_value: Option<Vec<crate::geometry::Geometry>>,
}

pub struct CanvasView {
    pub state: State,
}

#[derive(Debug, Clone)]
pub enum Message {}

impl Default for CanvasView {
    fn default() -> Self {
        Self {
            state: State::default(),
        }
    }
}

impl CanvasView {
    pub fn view<'a>(&'a mut self) -> Element<'a, Message> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

fn to_color(c: &crate::geometry::Color) -> Color {
    Color {
        r: c.r as f32 / 255.0,
        g: c.g as f32 / 255.0,
        b: c.b as f32 / 255.0,
        a: c.t as f32 / 255.0,
    }
}

impl<'a> canvas::Program<Message> for CanvasView {
    fn draw(&self, bounds: Rectangle, cursor: Cursor) -> Vec<Geometry> {
        let mut frame = Frame::new(bounds.size());
        if self.state.geometry_value.is_some() {
            let g = self.state.geometry_value.as_ref().unwrap();
            for i in 0..g.len() {
                match &g[i] {
                    crate::geometry::Geometry::Rectangle(rect) => {
                        let r = Path::rectangle(
                            Point {
                                x: rect.p.x,
                                y: rect.p.y,
                            },
                            Size::new(rect.width, rect.height),
                        );
                        frame.fill(&r, to_color(&rect.fill_color));
                        let c = to_color(&rect.stroke_color);
                        frame.stroke(
                            &r,
                            Stroke::default()
                                .with_color(c)
                                .with_width(rect.stroke_width),
                        );
                    }
                    crate::geometry::Geometry::PolySegment(segs) => {
                        for i in 0..segs.p.len() - 1 {
                            let p = Path::line(
                                Point {
                                    x: segs.p[i].x,
                                    y: segs.p[i].y,
                                },
                                Point {
                                    x: segs.p[i + 1].x,
                                    y: segs.p[i + 1].y,
                                },
                            );
                            let c = to_color(&segs.c);
                            frame.stroke(&p, Stroke::default().with_color(c).with_width(segs.w));
                        }
                    }
                    crate::geometry::Geometry::Segment(seg) => {
                        let p = Path::line(
                            Point {
                                x: seg.p1.x,
                                y: seg.p1.y,
                            },
                            Point {
                                x: seg.p2.x,
                                y: seg.p2.y,
                            },
                        );
                        let c = to_color(&seg.c);
                        frame.stroke(&p, Stroke::default().with_color(c).with_width(seg.w));
                    }
                    crate::geometry::Geometry::CircleWithTooltip(circ) => {
                        let circle = Path::circle(
                            Point {
                                x: circ.p.x,
                                y: circ.p.y,
                            },
                            circ.r,
                        );
                        let c = &circ.c;
                        frame.fill(&circle, to_color(c));
                    }
                    crate::geometry::Geometry::Circle(circ) => {
                        let circle = Path::circle(
                            Point {
                                x: circ.p.x,
                                y: circ.p.y,
                            },
                            circ.r,
                        );
                        let c = &circ.c;
                        frame.fill(&circle, to_color(c));
                    }
                    _ => {}
                };
            }
            let pos = cursor.position_in(&bounds);
            if pos.is_some() {
                for i in 0..g.len() {
                    match &g[i] {
                        crate::geometry::Geometry::CircleWithTooltip(circ) => {
                            let xdiff = pos.unwrap().x - circ.p.x;
                            let ydiff = pos.unwrap().y - circ.p.y;
                            let dist = (xdiff * xdiff + ydiff * ydiff).sqrt();
                            if dist <= circ.r {
                                let mut stext = circ.t.clone();
                                stext = stext.replace(",", "\n");
                                frame.translate(Vector { x: 400.0, y: 10.0 });
                                let text = canvas::Text {
                                    content: stext,
                                    size: 25.0,
                                    ..canvas::Text::default()
                                };
                                frame.fill_text(text);
                                frame.translate(Vector {
                                    x: -400.0,
                                    y: -10.0,
                                });
                                break;
                            }
                        }
                        _ => {}
                    };
                }
            }
        }
        vec![frame.into_geometry()]
    }
}
