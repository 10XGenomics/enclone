// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Code for working with polygons.  Possibly rename to planar_objects.rs.

use core::mem::swap;
use float_ord::FloatOrd;
use std::f64::consts::PI;
use superslice::Ext;
use vector_utils::unique_sort;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
//
// POINTS AND SEGMENTS
//
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Point in two-space.

#[derive(Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn origin_dist(&self) -> f64 {
        ((self.x * self.x) + (self.y * self.y)).sqrt()
    }

    pub fn dist(&self, p: Point) -> f64 {
        let (dx, dy) = (self.x - p.x, self.y - p.y);
        (dx * dx + dy * dy).sqrt()
    }

    pub fn dist_to_segment(&self, s: Segment) -> f64 {
        let mut d = self.dist(s.p1).min(self.dist(s.p2));
        if s.p1.x == s.p2.x {
            if (self.y >= s.p1.y && self.y <= s.p2.y) || (self.y >= s.p2.y && self.y <= s.p1.y) {
                d = d.min((self.x - s.p1.x).abs());
            }
        } else {
            // The line L extending the segment is Y = mX + b.
            let m = (s.p2.y - s.p1.y) / (s.p2.x - s.p1.x);
            let b = s.p1.y - m * s.p1.x;
            let x;
            let y;
            // Find the intersection (x, y) of the line M orthogonal to L and passing through self.
            // First suppose L is horizontal.
            if m == 0.0 {
                // The line M is given by X = self.x.
                x = self.x;
                y = m * x + b;
            // Now suppose L is not horizontal.
            } else {
                // The line M is given by Y = nX + c.
                let n = -1.0 / m;
                let c = self.y - n * self.x;
                x = (c - b) / (m - n);
                y = m * x + b;
            }
            // Determine if the intersection lies on the segment.
            if (x >= s.p1.x && x <= s.p2.x) || (x >= s.p2.x && x <= s.p1.x) {
                d = d.min(self.dist(Point { x, y }));
            }
        }
        d
    }
}

// Line segment in two-space.

pub struct Segment {
    pub p1: Point,
    pub p2: Point,
}

impl Segment {
    pub fn mid(&self) -> Point {
        Point {
            x: (self.p1.x + self.p2.x) / 2.0,
            y: (self.p1.y + self.p2.y) / 2.0,
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
//
// INTERVALS
//
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Interval in one-space.

#[derive(Clone, Default)]
pub struct Interval {
    pub x1: f64,
    pub x2: f64,
}

impl Interval {
    pub fn touches_interval(&self, m: Interval) -> bool {
        // Does [p.y - d, p.y + d] overlap [y1, y2]?
        if m.x1 <= self.x1 && m.x2 >= self.x1 {
            return true;
        } else if m.x1 <= self.x2 && m.x2 >= self.x2 {
            return true;
        } else if m.x1 >= self.x1 && m.x2 <= self.x2 {
            return true;
        }
        false
    }
}

// Interval vector with fast lookup of which closed intervals a number lies in.  We allow
// degenerate intervals [x, x].

#[derive(Clone, Default)]
pub struct IntervalVec {
    // PUBLIC
    pub is: Vec<Interval>,
    // PRIVATE
    ends: Vec<FloatOrd<f64>>,
    locs: Vec<Vec<usize>>,
    locsi: Vec<Vec<usize>>,
}

impl IntervalVec {
    // precompute is O(n * ln(n)).  Call this before calling get_containers.

    pub fn precompute(&mut self) {
        self.ends.clear();
        self.locs.clear();
        for i in 0..self.is.len() {
            assert!(self.is[i].x1 <= self.is[i].x2);
            self.ends.push(FloatOrd(self.is[i].x1));
            self.ends.push(FloatOrd(self.is[i].x2));
        }
        unique_sort(&mut self.ends);
        self.locs = vec![Vec::new(); self.ends.len() + 2];
        self.locsi = vec![Vec::new(); self.ends.len() + 2];
        for i in 0..self.is.len() {
            let p1 = self.ends.lower_bound(&FloatOrd(self.is[i].x1));
            let p2 = self.ends.lower_bound(&FloatOrd(self.is[i].x2));
            for j in p1..=p2 {
                self.locs[j + 1].push(i);
            }
            for j in p1..p2 {
                self.locsi[j + 1].push(i);
            }
        }
    }

    // Find the indices of the closed intervals that contain a given number.  You need to have
    // called precompute first.  This function is O(ln(n) + c), where n is the total number of
    // intervals in the vector and c is the number of intervals that contain the number.

    pub fn get_containers(&self, x: f64) -> &Vec<usize> {
        let x = FloatOrd(x);
        if self.is.is_empty() {
            return &self.locs[0];
        } else if x > *self.ends.last().unwrap() {
            return &self.locs[0];
        } else if x == *self.ends.last().unwrap() {
            return self.locs.last().unwrap();
        }
        let m = self.ends.upper_bound(&x);
        if m == 0 || self.ends[m - 1] == x {
            &self.locs[m]
        } else {
            &self.locsi[m]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_trace::PrettyTrace;
    #[test]
    fn test_interval_vec() {
        PrettyTrace::new().on();
        let mut v = IntervalVec::default();
        v.is.push(Interval { x1: 1.1, x2: 1.2 });
        v.is.push(Interval { x1: 1.0, x2: 2.0 });
        v.is.push(Interval { x1: 0.0, x2: 1.0 });
        v.is.push(Interval {
            x1: 10.0,
            x2: 10.01,
        });
        v.is.push(Interval { x1: 3.0, x2: 3.0 });
        v.is.push(Interval { x1: 1.0, x2: 1.0 });
        v.precompute();
        for x in [
            -1.0, 0.0, 0.5, 1.0, 1.01, 1.15, 1.8, 2.0, 3.1, 10.0, 10.05, 10.5,
        ]
        .iter()
        {
            let mut y = Vec::<usize>::new();
            for i in 0..v.is.len() {
                if *x >= v.is[i].x1 && *x <= v.is[i].x2 {
                    y.push(i);
                }
            }
            assert!(*v.get_containers(*x) == y);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
//
// POLYGONS
//
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Polygon structure.  The vertices represent a clockwise traversal and do not repeat the first
// vertex at the end.  The polygon is not assumed to be convex.  We assume that the edges do not
// cross each other, but that is not tested.

#[derive(Clone, Default)]
pub struct Polygon {
    // PUBLIC
    pub v: Vec<Point>,
    // PRIVATE
    xranges: IntervalVec,
}

impl Polygon {
    // Precompute xranges.

    pub fn precompute(&mut self) {
        for i1 in 0..self.v.len() {
            let i2 = (i1 + 1) % self.v.len();
            let mut x1 = self.v[i1].x;
            let mut x2 = self.v[i2].x;
            if x1 > x2 {
                swap(&mut x1, &mut x2);
            }
            self.xranges.is.push(Interval { x1, x2 });
        }
        self.xranges.precompute();
    }

    // Enlarge the polygon by a specified distance d.  This moves each vertex d farther away
    // from the center of mass, as defined by the vertices.  Not sure that makes sense.

    pub fn enlarge(&mut self, d: f64) {
        let (mut cx, mut cy) = (0.0, 0.0);
        for i in 0..self.v.len() {
            cx += self.v[i].x;
            cy += self.v[i].y;
        }
        cx /= self.v.len() as f64;
        cy /= self.v.len() as f64;
        let c = Point { x: cx, y: cy };
        for i in 0..self.v.len() {
            let x = self.v[i].x;
            let y = self.v[i].y;
            let dc = self.v[i].dist(c);
            let r = (dc + d) / dc;
            self.v[i].x = cx + (x - cx) * r;
            self.v[i].y = cy + (y - cy) * r;
        }
    }

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
        svg
    }

    // Generate an svg representation of a smooth path that passes through the vertices of the
    // polygon, and whose points lie within distance d of the polygon.  This is based on
    // catmull_bezier_svg.  The code is ridiculously inefficient and yields a suboptimal solution.

    pub fn catmull_bezier_bounded_svg(&self, d: f64) -> String {
        let mut p = self.clone();
        'loop_loop: loop {
            let last = p.v.len() - 1;
            for i in 0..=last + 1 {
                let mut is = Vec::<usize>::new();
                if i == 0 {
                    is.push(last);
                } else {
                    is.push(i - 1);
                }
                is.push(i % p.v.len());
                is.push((i + 1) % p.v.len());
                is.push((i + 2) % p.v.len());
                let (p0, p1, p2, p3) = (p.v[is[0]], p.v[is[1]], p.v[is[2]], p.v[is[3]]);

                // Define the control points.

                let c1 = Point {
                    x: (-p0.x + 6.0 * p1.x + p2.x) / 6.0,
                    y: (-p0.y + 6.0 * p1.y + p2.y) / 6.0,
                };
                let c2 = Point {
                    x: (p1.x + 6.0 * p2.x - p3.x) / 6.0,
                    y: (p1.y + 6.0 * p2.y - p3.y) / 6.0,
                };

                // Determine if the control points are farther than d from the three edges defined
                // by p0 ==> p1 ==> p2 ==> p3.  This is what makes the solution suboptimal.
                // Optimally, we would compute the distance of the farthest point on the actual
                // Bezier curve from the entire polygon (although that may not be worth the
                // effort).

                let mut dp1 = c1.dist_to_segment(Segment { p1: p0, p2: p1 });
                dp1 = dp1.min(c1.dist_to_segment(Segment { p1, p2 }));
                dp1 = dp1.min(c1.dist_to_segment(Segment { p1: p2, p2: p3 }));
                let mut dp2 = c2.dist_to_segment(Segment { p1: p0, p2: p1 });
                dp2 = dp2.min(c2.dist_to_segment(Segment { p1, p2 }));
                dp2 = dp2.min(c2.dist_to_segment(Segment { p1: p2, p2: p3 }));
                if dp1 > d || dp2 > d {
                    // Add three points to the polygon and start over.

                    let mut p_new = Polygon::default();
                    for j in 0..p.v.len() {
                        p_new.v.push(p.v[j]);
                        let mut k = None;
                        if j == is[0] {
                            k = Some(0);
                        } else if j == is[1] {
                            k = Some(1);
                        } else if j == is[2] {
                            k = Some(2);
                        }
                        if k.is_some() {
                            let k = k.unwrap();
                            // Add point halfway between pj and p[is[k+1]].
                            p_new.v.push(
                                Segment {
                                    p1: p.v[j],
                                    p2: p.v[is[k + 1]],
                                }
                                .mid(),
                            );
                        }
                    }
                    p = p_new;
                    continue 'loop_loop;
                }
            }
            break;
        }
        p.catmull_bezier_svg()
    }

    // Determine if a point lies inside the polygon (including the boundary).  This function is
    // O(n) where n is the number of vertices of the polygon.  See also inside_log.
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
                if x1 == x && ((y >= y1 && y <= y2) || (y >= y2 && y <= y1)) {
                    return true; // special case: on a vertical boundary edge
                }
            } else {
                let m = (y1 - y2) / (x1 - x2);
                let b = y2 - m * x2;
                let yy = m * x + b;
                if (yy >= y1 && yy <= y2) || (yy >= y2 && yy <= y1) {
                    if yy == y && ((x >= x1 && x <= x2) || (x >= x2 && x <= x1)) {
                        return true; // special case: on a non-vertical boundary edge
                    }
                    if yy >= y && (d.is_none() || (d.is_some() && d.unwrap() > yy - y)) {
                        d = Some(yy - y);
                        inside = x2 < x1;
                    }
                }
            }
        }
        inside
    }

    // Determine if a point lies inside the polygon (including the boundary), in O(ln(n)) time,
    // where n is the number of vertices of the polygon, except in pathological cases.  You need
    // to call precompute first.  See also inside.

    pub fn inside_log(&self, p: Point) -> bool {
        let (x, y) = (p.x, p.y);
        let hits = self.xranges.get_containers(x);
        let mut inside = false;
        let mut d = None;
        for j in 0..hits.len() {
            let i1 = hits[j];
            let i2 = (i1 + 1) % self.v.len();
            let (x1, y1) = (self.v[i1].x, self.v[i1].y);
            let (x2, y2) = (self.v[i2].x, self.v[i2].y);
            if x1 == x2 {
                if x1 == x && ((y >= y1 && y <= y2) || (y >= y2 && y <= y1)) {
                    return true; // special case: on a vertical boundary edge
                }
            } else {
                let m = (y1 - y2) / (x1 - x2);
                let b = y2 - m * x2;
                let yy = m * x + b;
                if (yy >= y1 && yy <= y2) || (yy >= y2 && yy <= y1) {
                    if yy == y && ((x >= x1 && x <= x2) || (x >= x2 && x <= x1)) {
                        return true; // special case: on a non-vertical boundary edge
                    }
                    if yy >= y && (d.is_none() || (d.is_some() && d.unwrap() > yy - y)) {
                        d = Some(yy - y);
                        inside = x2 < x1;
                    }
                }
            }
        }
        inside
    }

    // Determine if the disk centered at p with radius r touches the polygon.  This assumes that
    // you've run precompute on the polygon.

    pub fn touches_disk(&self, p: Point, r: f64) -> bool {
        // Case 1: the point is inside the polygon.

        if self.inside_log(p) {
            return true;
        }

        // Case 2: the polygon is inside the disk.

        let mut enclosed = true;
        for i in 0..self.v.len() {
            if p.dist(self.v[i]) > r {
                enclosed = false;
                break;
            }
        }
        if enclosed {
            return true;
        }

        // Case 3: the disk boundary (circle) touches an edge.

        for i1 in 0..self.v.len() {
            let i2 = (i1 + 1) % self.v.len();
            let (p1, p2) = (self.v[i1], self.v[i2]);
            if p1.x == p2.x {
                if (p1.x - p.x).abs() <= r {
                    let dx = p1.x - p.x;
                    let d = (r * r - dx * dx).sqrt();
                    let (mut y1, mut y2) = (p1.y, p2.y);
                    if y1 > y2 {
                        swap(&mut y1, &mut y2);
                    }
                    let i1 = Interval {
                        x1: p.y - d,
                        x2: p.y + d,
                    };
                    let i2 = Interval { x1: y1, x2: y2 };
                    if i1.touches_interval(i2) {
                        return true;
                    }
                }
            } else {
                // Let Y = mX + b be the line from p1 to p2.
                let m = (p2.y - p1.y) / (p2.x - p1.x);
                let b = p1.y - m * p1.x;
                // Y = mX + b (from p1 to p2)  ==?overlaps?== (X - p.x)^2 + (Y - p.y)^2 = r^2
                // ==> solve (X - p.x)^2 + (mX + b - p.y)^2 = r^2
                // ==> X^2 - 2p.x*X + p.x^2 + m^2*X^2 + (b - p.y)^2 + 2 * mX * (b - p.y) = r^2
                // ==> (1 + m^2)X^2 + (2m * (b - p.y) - 2p.x)*X = r^2 - p.x^2 - (b - p.y)^2
                let c = 1.0 + m * m;
                let mut d = m * (b - p.y) - p.x;
                let mut u = r * r - p.x * p.x - (b - p.y) * (b - p.y);
                // ==> c*X^2 + 2*d*X = u
                d /= c;
                u /= c;
                // ==> X^2 + 2*d*X = u
                // ==> (X + d)^2 = u + d^2
                u += d * d;
                // ==> (X + d)^2 = u
                if u >= 0.0 {
                    let z1 = -u.sqrt() - d;
                    let z2 = u.sqrt() - d;
                    let (mut x1, mut x2) = (p1.x, p2.x);
                    if x1 > x2 {
                        swap(&mut x1, &mut x2);
                    }
                    let i1 = Interval { x1: z1, x2: z2 };
                    let i2 = Interval { x1, x2 };
                    if i1.touches_interval(i2) {
                        return true;
                    }
                }
            }
        }
        false
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
//
// ENCLOSING POLYGON
//
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
//
// WHICH CIRCLE
//
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a collection of circles, determine which if any circles contain a given point.  Requires
// precompute.

pub fn precompute_for_circle_containment(circles: &Vec<(f64, f64, f64)>) -> IntervalVec {
    let mut x = IntervalVec::default();
    x.is.resize(circles.len(), Interval::default());
    for (i, c) in circles.iter().enumerate() {
        x.is[i].x1 = c.0 - c.2;
        x.is[i].x2 = c.0 + c.2;
    }
    x.precompute();
    x
}

pub fn circles_containing_point(
    p: Point,
    circles: &Vec<(f64, f64, f64)>,
    circles_precompute: &IntervalVec,
) -> Vec<usize> {
    let xhome = circles_precompute.get_containers(p.x);
    let mut c = Vec::<usize>::new();
    for i in xhome.iter() {
        let (x, y, r) = &circles[*i];
        let (xdiff, ydiff) = (p.x - x, p.y - y);
        if xdiff * xdiff + ydiff * ydiff <= r * r {
            c.push(*i);
        }
    }
    c
}
