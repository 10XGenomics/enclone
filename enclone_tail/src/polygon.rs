// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Code for working with polygons.

use std::f64::consts::PI;

#[derive(Clone, Copy)]
pub struct Point {
    x: f64,
    y: f64,
}

// Polygon structure.  The vertices represent a clockwise traversal and do not repeat the first
// vertex at the end.  The polygon is not assumed to be convex.  We assume that the edges do not
// cross each other, but that is not tested.

#[derive(Default)]
pub struct Polygon {
    pub v: Vec<Point>,
}

impl Polygon {
    // Generate an svg path for the polygon.

    pub fn to_svg(&self) -> String {
        let mut svg = format!("M {:.4} {:.4}\n", self.v[0].x, self.v[0].y);
        for i in 0..=self.v.len() {
            let i = i % self.v.len();
            svg += &format!("L {:.4} {:.4}\n", self.v[i].x, self.v[i].y);
        }
        svg
    }

    // Generate an svg representation of a smooth path that passes through the vertices of the
    // polygon.  The code was derived trivially from src/chart.rs in lorikeet-dash crate
    // version 0.1.0, whose license is MIT/Apache-2.0.

    pub fn catmull_bezier_svg(&self) -> String {
        let mut svg = format!("M {:.4} {:.4}\n", self.v[0].x, self.v[0].y);
        let last = self.v.len() - 1;
        for i in 0..=last + 1 {
            let p0 = if i == 0 { self.v[last] } else { self.v[i - 1] };
            let p1 = self.v[i % self.v.len()];
            let p2 = self.v[(i + 1) % self.v.len()];
            let p3 = self.v[(i + 2) % self.v.len()];
            let c1 = Point {
                x: (-p0.x + 6.0 * p1.x + p2.x) / 6.0,
                y: (-p0.y + 6.0 * p1.y + p2.y) / 6.0,
            };
            let c2 = Point {
                x: (p1.x + 6.0 * p2.x - p3.x) / 6.0,
                y: (p1.y + 6.0 * p2.y - p3.y) / 6.0,
            };
            let end = p2;
            svg += &format!(
                "C {:.4} {:.4}, {:.4} {:.4}, {:.4} {:.4}\n",
                c1.x, c1.y, c2.x, c2.y, end.x, end.y
            );
        }
        return svg;
    }

    // Determine if a point lies inside the polygon (including the boundary).  This function is
    // O(n) where n is the number of vertices of the polygon.  Note that with some precompute for
    // the polygon, this could probably be made O(ln(n)), at least for non-pathological cases.
    // Reference: http://geomalgorithms.com/a03-_inclusion.html (not used).

    pub fn inside(&self, p: Point) -> bool {
        let (x, y) = (p.x, p.y);

        // Find the closest point below (x,y) on the vertical line X = x that lies on an edge
        // of the polygon.

        let mut inside = false;
        let mut d = None;
        for i1 in 0..self.v.len() {
            let i2 = (i1 + 1) % self.v.len();
            let (x1, y1) = (self.v[i1].x, self.v[i1].y);
            let (x2, y2) = (self.v[i2].x, self.v[i2].y);
            if x1 == x2 {
                if x1 == x {
                    if (y >= y1 && y <= y2) || (y >= y2 && y <= y1) {
                        return true; // special case: on a vertical boundary edge
                    }
                }
            } else {
                let m = (y1 - y2) / (x1 - x2);
                let b = y2 - m * x2;
                let yy = m * x + b;
                if (yy >= y1 && yy <= y2) || (yy >= y2 && yy <= y1) {
                    if yy == y && ((x >= x1 && x <= x2) || (x >= x2 && x <= x1)) {
                        return true; // special case: on a non-vertical boundary edge
                    }
                    if yy >= y {
                        if d.is_none() || (d.is_some() && d.unwrap() > yy - y) {
                            d = Some(yy - y);
                            inside = x2 < x1;
                        }
                    }
                }
            }
        }
        inside
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find an n-vertex polygon that encloses circles c = {(r, x, y)} with separation of at least
// d from them.

pub fn enclosing_polygon(c: &Vec<(f64, f64, f64)>, d: f64, n: usize) -> Polygon {
    // First find the center of mass (mx, my);

    let (mut m, mut mx, mut my) = (0.0, 0.0, 0.0);
    for i in 0..c.len() {
        let w = c[i].0 * c[i].0;
        m += w;
        mx += w * c[i].1;
        my += w * c[i].2;
    }
    mx /= m;
    my /= m;

    // Initialize the distances of the polygonal vertices from the center of mass.

    let mut pr = vec![0.0_f64; n];

    // Traverse the polygonal vertices.

    for i in 0..n {
        let theta1 = 2.0 * PI * (i as f64) / (n as f64);
        let theta2 = 2.0 * PI * ((i + 1) as f64) / (n as f64);
        for j in 0..c.len() {
            let (r, x, y) = (c[j].0, c[j].1, c[j].2);

            // Find the maximum distance of a point on the circle from the center of mass, with the
            // angle from the center of mass confined between theta1 and theta2.  However, instead
            // of solving the math problem, do this for a bunch of points on the circle.  This is
            // lazy.

            let np = 100;
            for k in 0..np {
                let rho = 2.0 * PI * (k as f64) / (np as f64);
                let px = rho.cos() * r + x;
                let py = rho.sin() * r + y;
                let (dx, dy) = (px - mx, py - my);
                let mut theta = dy.atan2(dx);
                if theta < 0.0 {
                    theta += 2.0 * PI;
                }
                if theta1 <= theta && theta <= theta2 {
                    pr[i] = pr[i].max((dx * dx + dy * dy).sqrt() + d);
                    let mut ip = i + 1;
                    if ip == n {
                        ip = 0;
                    }
                    pr[ip] = pr[ip].max((dx * dx + dy * dy).sqrt() + d);
                }
            }
        }
    }

    // Form the polygon.

    let mut p = Polygon::default();
    for i in 0..n {
        let theta = 2.0 * PI * (i as f64) / (n as f64);
        p.v.push(Point {
            x: mx + pr[i] * theta.cos(),
            y: my + pr[i] * theta.sin(),
        });
    }
    p
}
