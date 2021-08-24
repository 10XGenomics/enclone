// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::GROUP_ID;
use crate::GROUP_ID_CLICKED_ON;
use crate::*;
use iced::canvas::event::{self, Event};
use iced::{
    canvas::{self, Canvas, Cursor, Frame, Geometry, Path, Stroke, Text},
    mouse, Color, Element, HorizontalAlignment, Length, Rectangle, Size, VerticalAlignment,
};
use iced_native::{Font, Point, Vector};
use lazy_static::lazy_static;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering::SeqCst;
use std::sync::Mutex;

const LIBERATION_SANS: Font = Font::External {
    name: "LIBERATION_SANS",
    bytes: include_bytes!("../../fonts/LiberationSans-Regular.ttf"),
};

lazy_static! {
    pub static ref IN_GEOMETRIES: Mutex<Vec<crate::geometry::Geometry>> =
        Mutex::new(Vec::<crate::geometry::Geometry>::new());
    pub static ref OUT_GEOMETRIES: Mutex<Vec<Geometry>> = Mutex::new(Vec::<Geometry>::new());
    pub static ref OUT_GEOMETRIES_TOOLTIP: Mutex<Vec<Geometry>> =
        Mutex::new(Vec::<Geometry>::new());
    pub static ref POS: Mutex<Point> = Mutex::new(Point::default());
}
pub static POS_IS_SOME: AtomicBool = AtomicBool::new(false);

#[derive(Default)]
pub struct State {
    pub geometry_value: Option<Vec<crate::geometry::Geometry>>,
}

pub struct CanvasView {
    pub state: State,
}

#[derive(Debug, Clone)]
pub enum Message {
    DoNothing,
    GroupClick,
}

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

impl CanvasView {
    fn dimensions(&self) -> (f32, f32) {
        let g = self.state.geometry_value.as_ref().unwrap();
        let mut height = 0.0 as f32;
        let mut width = 0.0 as f32;
        for i in 0..g.len() {
            match &g[i] {
                crate::geometry::Geometry::Text(o) => {
                    height = height.max(o.p.y);
                    // not right: need to add text length
                    width = width.max(o.p.x);
                }
                crate::geometry::Geometry::Rectangle(rect) => {
                    height = height.max(rect.p.y + rect.height);
                    width = width.max(rect.p.x + rect.width);
                }
                crate::geometry::Geometry::PolySegment(segs) => {
                    for i in 0..segs.p.len() - 1 {
                        height = height.max(segs.p[i].y);
                        width = width.max(segs.p[i].x);
                    }
                }
                crate::geometry::Geometry::Segment(seg) => {
                    height = height.max(seg.p1.y);
                    height = height.max(seg.p2.y);
                    width = width.max(seg.p1.x);
                    width = width.max(seg.p2.x);
                }
                crate::geometry::Geometry::CircleWithTooltip(circ) => {
                    height = height.max(circ.p.y + circ.r);
                    width = width.max(circ.p.x + circ.r);
                }
                crate::geometry::Geometry::Circle(circ) => {
                    height = height.max(circ.p.y + circ.r);
                    width = width.max(circ.p.x + circ.r);
                }
            };
        }
        (width, height)
    }
}

const MAX_WIDTH: f32 = 770.0;

fn get_scale(width: f32, height: f32, empty: bool) -> f32 {
    let mut max_height = SVG_HEIGHT as f32;
    if empty {
        max_height = SVG_NULL_HEIGHT as f32;
    }
    max_height -= 5.0;
    let scale_x = MAX_WIDTH / width;
    let scale_y = max_height / height;
    let mut scale = scale_y;
    if scale_x < scale_y {
        scale = scale_x;
    }
    scale
}

impl<'a> canvas::Program<Message> for CanvasView {
    fn update(
        &mut self,
        event: Event,
        bounds: Rectangle,
        cursor: Cursor,
    ) -> (event::Status, Option<Message>) {
        let _cursor_position = if let Some(position) = cursor.position_in(&bounds) {
            position
        } else {
            return (event::Status::Ignored, None);
        };
        match event {
            Event::Mouse(mouse_event) => match mouse_event {
                mouse::Event::ButtonPressed(button) => {
                    let message = match button {
                        mouse::Button::Left => {
                            let g = self.state.geometry_value.as_ref().unwrap();
                            let (width, height) = self.dimensions();
                            let scale = get_scale(width, height, g.len() == 1);
                            let mut group_id = None;
                            let pos = cursor.position_in(&bounds);
                            for i in 0..g.len() {
                                match &g[i] {
                                    crate::geometry::Geometry::CircleWithTooltip(circ) => {
                                        let xdiff = pos.unwrap().x - circ.p.x * scale;
                                        let ydiff = pos.unwrap().y - circ.p.y * scale;
                                        let dist = (xdiff * xdiff + ydiff * ydiff).sqrt();
                                        if dist <= circ.r {
                                            let stext = circ.t.clone();
                                            let xs = stext.split(',').collect::<Vec<&str>>();
                                            for j in 0..xs.len() {
                                                if xs[j].starts_with("group_id=") {
                                                    group_id = Some(xs[j].after("=").force_usize());
                                                }
                                            }
                                            let group_id = group_id.unwrap();
                                            GROUP_ID_CLICKED_ON.store(true, SeqCst);
                                            GROUP_ID.store(group_id, SeqCst);
                                            break;
                                        }
                                    }
                                    _ => {}
                                }
                            }
                            if group_id.is_some() {
                                Some(Message::GroupClick)
                            } else {
                                None
                            }
                        }
                        _ => None,
                    };
                    (event::Status::Captured, message)
                }
                _ => (event::Status::Captured, None),
            },
            _ => (event::Status::Captured, None),
        }
    }

    fn draw(&self, bounds: Rectangle, cursor: Cursor) -> Vec<Geometry> {
        // Suppose there is no geometry.

        if self.state.geometry_value.is_none() {
            POS_IS_SOME.store(false, SeqCst);
            IN_GEOMETRIES.lock().unwrap().clear();
            return Vec::new();
        }

        // Suppose we have geometries.

        let pos = cursor.position_in(&bounds);

        // Use cached geometries if possible.
        //
        // Annoyingly, can't do this:
        // if pos_changed = POS.lock().unwrap() != pos {
        // or
        // if IN_GEOMETRIES.lock().unwrap() == self.state.geometry_value.as_ref().unwrap() {

        let mut pos_same = true;
        let last_pos_is_some = POS_IS_SOME.load(SeqCst);
        if !last_pos_is_some && pos.is_none() {
        } else if last_pos_is_some && pos.is_none() {
            pos_same = false;
        } else if !last_pos_is_some && pos.is_some() {
            pos_same = false;
        } else if POS.lock().unwrap().x != pos.unwrap().x {
            pos_same = false;
        } else if POS.lock().unwrap().y != pos.unwrap().y {
            pos_same = false;
        }
        let mut geom_same = true;
        let in_geometries = &self.state.geometry_value.as_ref().unwrap();
        if IN_GEOMETRIES.lock().unwrap().len() != in_geometries.len() {
            geom_same = false;
        } else {
            for i in 0..in_geometries.len() {
                if IN_GEOMETRIES.lock().unwrap()[i] != in_geometries[i] {
                    geom_same = false;
                    break;
                }
            }
        }
        if pos_same && geom_same {
            let mut v = OUT_GEOMETRIES.lock().unwrap().clone();
            v.append(&mut OUT_GEOMETRIES_TOOLTIP.lock().unwrap().clone());
            return v;
        }
        if !geom_same {
            IN_GEOMETRIES.lock().unwrap().clear();
            IN_GEOMETRIES
                .lock()
                .unwrap()
                .append(&mut self.state.geometry_value.as_ref().unwrap().clone());
        }
        if pos.is_none() {
            POS_IS_SOME.store(false, SeqCst);
        } else {
            POS.lock().unwrap().x = pos.unwrap().x;
            POS.lock().unwrap().y = pos.unwrap().y;
        }

        // Compute width and height and scale.
        //
        // for now not scaling stroke width, not sure what is optimal
        // scaling seems to be needed only because .height(SVG_WIDTH) doesn't work on a canvas
        // should file bug
        // This could be 400.0 but we would need to take account of text width
        // in computing the max.

        let g = self.state.geometry_value.as_ref().unwrap();
        let (width, height) = self.dimensions();
        let scale = get_scale(width, height, g.len() == 1);

        // Rebuild geometries if needed.

        let mut v;
        if geom_same {
            v = OUT_GEOMETRIES.lock().unwrap().clone();
        } else {
            let mut frame = Frame::new(bounds.size());
            for i in 0..g.len() {
                match &g[i] {
                    crate::geometry::Geometry::Text(o) => {
                        // rotate not implemented because not a feature yet in iced
                        let x = Text {
                            content: o.t.clone(),
                            size: o.font_size * scale,
                            color: to_color(&o.c),
                            position: Point {
                                x: o.p.x * scale,
                                // font bit is compensation for vertical alignment issue
                                // don't understand why / 4.0 makes sense
                                y: o.p.y * scale + o.font_size * scale / 4.0,
                            },
                            font: match o.font.as_str() {
                                "DejaVuSansMono" => DEJAVU,
                                "Arial" => LIBERATION_SANS,
                                _ => LIBERATION_SANS,
                            },
                            // Center doesn't seem to work, should report bug
                            // nor does Top
                            vertical_alignment: VerticalAlignment::Bottom,
                            horizontal_alignment: match o.halign {
                                crate::geometry::HorizontalAlignment::Left => {
                                    HorizontalAlignment::Left
                                }
                                crate::geometry::HorizontalAlignment::Center => {
                                    HorizontalAlignment::Center
                                }
                                crate::geometry::HorizontalAlignment::Right => {
                                    HorizontalAlignment::Right
                                }
                            },
                        };
                        frame.fill_text(x);
                    }
                    crate::geometry::Geometry::Rectangle(rect) => {
                        let r = Path::rectangle(
                            Point {
                                x: rect.p.x * scale,
                                y: rect.p.y * scale,
                            },
                            Size::new(rect.width * scale, rect.height * scale),
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
                                    x: segs.p[i].x * scale,
                                    y: segs.p[i].y * scale,
                                },
                                Point {
                                    x: segs.p[i + 1].x * scale,
                                    y: segs.p[i + 1].y * scale,
                                },
                            );
                            let c = to_color(&segs.c);
                            frame.stroke(&p, Stroke::default().with_color(c).with_width(segs.w));
                        }
                    }
                    crate::geometry::Geometry::Segment(seg) => {
                        let p = Path::line(
                            Point {
                                x: seg.p1.x * scale,
                                y: seg.p1.y * scale,
                            },
                            Point {
                                x: seg.p2.x * scale,
                                y: seg.p2.y * scale,
                            },
                        );
                        let c = to_color(&seg.c);
                        frame.stroke(&p, Stroke::default().with_color(c).with_width(seg.w));
                    }
                    crate::geometry::Geometry::CircleWithTooltip(circ) => {
                        let circle = Path::circle(
                            Point {
                                x: circ.p.x * scale,
                                y: circ.p.y * scale,
                            },
                            circ.r * scale,
                        );
                        let c = &circ.c;
                        frame.fill(&circle, to_color(c));
                    }
                    crate::geometry::Geometry::Circle(circ) => {
                        let circle = Path::circle(
                            Point {
                                x: circ.p.x * scale,
                                y: circ.p.y * scale,
                            },
                            circ.r * scale,
                        );
                        let c = &circ.c;
                        frame.fill(&circle, to_color(c));
                    }
                };
            }
            v = vec![frame.into_geometry()];
            OUT_GEOMETRIES.lock().unwrap().clear();
            OUT_GEOMETRIES.lock().unwrap().append(&mut v.clone());
        }
        if pos_same && pos.is_some() {
            v.append(&mut OUT_GEOMETRIES_TOOLTIP.lock().unwrap().clone());
        }
        if !pos_same && pos.is_some() {
            let mut frame = Frame::new(bounds.size());
            for i in 0..g.len() {
                match &g[i] {
                    crate::geometry::Geometry::CircleWithTooltip(circ) => {
                        let xdiff = pos.unwrap().x - circ.p.x * scale;
                        let ydiff = pos.unwrap().y - circ.p.y * scale;
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
                            print_tabular_vbox(&mut log, &rows, 0, &b"l|r".to_vec(), false, true);
                            let xpos = 15.0 + width * scale;
                            frame.translate(Vector { x: xpos, y: 0.0 });

                            TOOLTIP_TEXT.lock().unwrap().clear();
                            TOOLTIP_TEXT.lock().unwrap().push(log.clone());

                            let mut logp = String::new();
                            for char in log.chars() {
                                if char == '\n' {
                                    logp.push(char);
                                } else {
                                    logp.push('█');
                                }
                            }

                            // We put a layer of black below the tooltip text, which is going to
                            // be white.  There are two approaches.  First, if the tooltip box lies
                            // strictly within the canvas (and not to the right of it), we display
                            // a black rectangle.  Otherwise, we contruct the layer out of box
                            // characters.  This is not fully satisfactory because there are small
                            // gaps between them.

                            let tooltip_font_size: f32 = 13.5;
                            let mut width_in_chars = 0;
                            let mut height_in_chars = 0;
                            for line in log.lines() {
                                let mut nchars = 0;
                                for _char in line.chars() {
                                    nchars += 1;
                                }
                                width_in_chars = std::cmp::max(width_in_chars, nchars);
                                height_in_chars += 1;
                            }
                            let box_width = width_in_chars as f32
                                * tooltip_font_size
                                * DEJAVU_WIDTH_OVER_HEIGHT;
                            let box_height = height_in_chars as f32 * tooltip_font_size;
                            if xpos + box_width <= MAX_WIDTH {
                                frame.fill_rectangle(
                                    Point { x: 0.0, y: 0.0 },
                                    Size {
                                        width: box_width,
                                        height: box_height,
                                    },
                                    iced::canvas::Fill::from(Color::BLACK),
                                );
                            } else {
                                let text = canvas::Text {
                                    content: logp,
                                    size: tooltip_font_size,
                                    font: DEJAVU,
                                    color: Color::from_rgb(0.0, 0.0, 0.0),
                                    ..canvas::Text::default()
                                };
                                frame.fill_text(text);
                            }

                            // Now display the actual text in the tooltip box.

                            let text = canvas::Text {
                                content: log,
                                size: tooltip_font_size,
                                font: DEJAVU,
                                color: Color::from_rgb(1.0, 1.0, 1.0),
                                ..canvas::Text::default()
                            };
                            frame.fill_text(text);

                            frame.translate(Vector { x: -xpos, y: -10.0 });
                            break;
                        }
                    }
                    _ => {}
                };
            }
            let mut w = vec![frame.into_geometry()];
            OUT_GEOMETRIES_TOOLTIP.lock().unwrap().clear();
            OUT_GEOMETRIES_TOOLTIP
                .lock()
                .unwrap()
                .append(&mut w.clone());
            v.append(&mut w);
        }
        v
    }
}
