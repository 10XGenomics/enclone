// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use iced::{
    canvas::{self, Canvas, Cursor, Frame, Geometry, Path, Stroke},
    Color, Element, Length, Rectangle, Size,
};
use iced_native::{Font, Point, Vector};
use string_utils::*;
use tables::*;

// not bold
const DEJAVU: Font = Font::External {
    name: "DEJAVU",
    bytes: include_bytes!("../../fonts/DejaVuLGCSansMono.ttf"),
};

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
                                let stext = circ.t.clone();
                                let xs = stext.split(',').collect::<Vec<&str>>();
                                let mut rows = Vec::<Vec<String>>::new();
                                for i in 0..xs.len() {
                                    if i > 0 {
                                        rows.push(vec!["\\hline".to_string(); 2]);
                                    }
                                    let mut row = Vec::<String>::new();
                                    row.push(xs[i].before("=").to_string());
                                    row.push(xs[i].after("=").to_string());
                                    rows.push(row);
                                }
                                let mut log = String::new();
                                print_tabular_vbox(&mut log, &rows, 1, &b"l|r".to_vec(), 
                                    false, true);
                                frame.translate(Vector { x: 400.0, y: 10.0 });
                                let text = canvas::Text {
                                    content: log,
                                    size: 20.0,
                                    font: DEJAVU,
                                    color: Color::from_rgb(0.8, 0.4, 0.6),
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
