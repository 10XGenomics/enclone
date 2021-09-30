// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::polygon::{Point, Polygon};
use rayon::prelude::*;

// Pack circles of given radii, which should be in descending order.  Return centers for the
// circles.  There is probably a literature on this, and this is probably a very crappy algorithm.
//
// Should google "d3 circle packing algorithm".  The code is open source.
//
// quad: force into first quadrant rather than anywhere
//
// Blacklisted polygons are avoided.

pub fn pack_circles(r: &Vec<f64>, blacklist: &Vec<Polygon>, quad: bool) -> Vec<(f64, f64)> {
    // Set up to track the centers of the circles.

    let mut c = Vec::<(f64, f64)>::new();
    if r.is_empty() {
        return c;
    }

    // Set up to track the x extents of the circles, which are intervals, and which we break into
    // subintervals equal to the extent of the smallest circle (and a possible fractional part).
    // We track the starting point of each subinterval, and the index of the circle.
    //
    // The overall code would presumably be faster if we did not store the data in a sorted
    // vector that needs to be sorted after each center is added.  However it's not actually
    // clear that much time is spent in the sort.

    let mut ints = Vec::<(f64, usize)>::new();
    let radius0 = *r.last().unwrap();

    // Define function to add a center.

    fn push_center(
        x: f64,
        y: f64,
        index: usize,
        radius: f64,
        radius0: f64,
        c: &mut Vec<(f64, f64)>,
        ints: &mut Vec<(f64, usize)>,
    ) {
        c.push((x, y));

        // need to cover [x-radius..x+radius] using segments of length radius0

        let mut pos = x - radius;
        while pos < x + radius {
            ints.push((pos, index));
            pos += radius0;
        }
        ints.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    // Define data for the first circle.

    if blacklist.is_empty() {
        if !quad {
            push_center(0.0, 0.0, 0, r[0], radius0, &mut c, &mut ints);
        } else {
            push_center(r[0], r[0], 0, r[0], radius0, &mut c, &mut ints);
        }
    }

    // Compute the maximum distance "bigr" from the origin.

    let mut bigr = if blacklist.is_empty() { r[0] } else { 0.0 };
    for p in blacklist.iter() {
        for x in p.v.iter() {
            bigr = bigr.max(x.origin_dist());
        }
    }

    // Proceed.

    let mut rand = 0i64;
    // We use a ridiculously large sample.  Reducing it to 10,000 noticeably reduces symmetry.
    // Presumably as the number of clusters increases, the sample would need to be increased
    // (ideally) to increase symmetry.
    const SAMPLE: usize = 100000;
    const MUL: f64 = 1.5;
    let mut centers = vec![(0.0, 0.0); SAMPLE];
    let start = if blacklist.is_empty() { 1 } else { 0 };
    let mut max_dist;
    if !quad {
        max_dist = 0.0;
    } else {
        max_dist = r[0];
    }
    for i in start..r.len() {
        let mut found = false;
        let mut best_r1 = 0.0;
        let mut best_r2 = 0.0;
        let mut best_val = 0.0;

        // Loop until we find a placement for the circle.

        loop {
            for z in 0..SAMPLE {
                // Get a random point in [-1,+1] x [-1,+1].  Using a hand-rolled random number
                // generator (from the internet) for speed and reproducibility, although there
                // might be something better in the standard packages.

                let rand1 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand)
                    .wrapping_add(1_442_695_040_888_963_407);
                let rand2 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand1)
                    .wrapping_add(1_442_695_040_888_963_407);
                rand = rand2;
                let mut r1;
                let mut r2;
                if !quad {
                    r1 = (2.0 * (rand1 % 1_000_000i64) as f64 / 1_000_000.0) - 1.0;
                    r2 = (2.0 * (rand2 % 1_000_000i64) as f64 / 1_000_000.0) - 1.0;
                } else {
                    r1 = ((rand1 % 1_000_000i64) as f64 / 1_000_000.0 + 1.0) / 2.0;
                    r2 = ((rand2 % 1_000_000i64) as f64 / 1_000_000.0 + 1.0) / 2.0;
                }

                // Make it bigger.

                if !quad {
                    r1 *= (bigr + r[i]) * MUL;
                    r2 *= (bigr + r[i]) * MUL;
                } else {
                    r1 = r1 * (bigr + r[i]) * MUL + r[i];
                    r2 = r2 * (bigr + r[i]) * MUL + r[i];
                }
                const RELAXATION_FACTOR: f64 = 0.8;
                if (r1 * r1 + r2 * r2).sqrt() < RELAXATION_FACTOR * max_dist {
                    continue;
                }
                centers[z] = (r1, r2);
            }
            let mut results = Vec::<(usize, f64, f64, bool)>::new();
            for z in 0..SAMPLE {
                results.push((z, 0.0, 0.0, false));
            }
            results.par_iter_mut().for_each(|res| {
                let z = res.0;
                let r1 = centers[z].0;
                let r2 = centers[z].1;

                // Test for overlap with a blacklisted polygon.

                let mut ok = true;
                let m = Point { x: r1, y: r2 };
                for p in blacklist.iter() {
                    if p.touches_disk(m, r[i]) {
                        ok = false;
                        break;
                    }
                }

                // See if circle at (r1,r2) with radius r[i] overlaps any of the existing circles.
                // This involves a binary search, which reduces the complexity of the algorithm.

                if ok {
                    let low = ints.binary_search_by(|v| {
                        v.partial_cmp(&(r1 - r[i] - radius0, 0))
                            .expect("Comparison failed.")
                    });
                    let low = match low {
                        Ok(i) => i,
                        Err(i) => i,
                    };
                    let high = ints.binary_search_by(|v| {
                        v.partial_cmp(&(r1 + r[i], 0)).expect("Comparison failed.")
                    });
                    let high = match high {
                        Ok(i) => i,
                        Err(i) => i,
                    };
                    for m in low..high {
                        if m > low && ints[m].1 == ints[m - 1].1 {
                            continue;
                        }
                        let k = ints[m].1;
                        let d = (c[k].0 - r1) * (c[k].0 - r1) + (c[k].1 - r2) * (c[k].1 - r2);
                        if d < (r[i] + r[k]) * (r[i] + r[k]) {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    res.1 = r1;
                    res.2 = r2;
                    res.3 = true;
                }
            });
            for z in 0..SAMPLE {
                if results[z].3 {
                    let (r1, r2) = (results[z].1, results[z].2);
                    let val;
                    if blacklist.is_empty() || i == 0 {
                        val = r1 * r1 + r2 * r2;
                    } else {
                        val = (r1 - c[0].0) * (r1 - c[0].0) + (r2 - c[0].1) * (r2 - c[0].1);
                    }
                    if !found || val < best_val {
                        best_r1 = r1;
                        best_r2 = r2;
                        best_val = val;
                    }
                    found = true;
                }
            }
            if found {
                break;
            }
        }
        max_dist = max_dist.max((best_r1 * best_r1 + best_r2 * best_r2).sqrt());
        push_center(best_r1, best_r2, i, r[i], radius0, &mut c, &mut ints);
        bigr = bigr.max(r[i] + (c[i].0 * c[i].0 + c[i].1 * c[i].1).sqrt());
    }
    c
}
