// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Make group colors.

pub fn make_group_colors(ngroups: usize) -> Vec<String> {
    let mut group_color = Vec::<String>::new();

    // Build colors.

    let sum = 40; // a+b+c, where color is rgb(255-a, 255-b, 255-c); smaller is closer to white
    let mut rand = 0i64;
    let mut points = Vec::<(f64, f64, f64)>::new();
    while points.len() < 10_000 {
        let rand1 = 6_364_136_223_846_793_005i64
            .wrapping_mul(rand)
            .wrapping_add(1_442_695_040_888_963_407);
        let rand2 = 6_364_136_223_846_793_005i64
            .wrapping_mul(rand1)
            .wrapping_add(1_442_695_040_888_963_407);
        let rand3 = 6_364_136_223_846_793_005i64
            .wrapping_mul(rand2)
            .wrapping_add(1_442_695_040_888_963_407);
        rand = rand3;
        let mut r1 = (rand1 % 1_000_000i64) as f64 / 1_000_000.0;
        let mut r2 = (rand2 % 1_000_000i64) as f64 / 1_000_000.0;
        let mut r3 = (rand3 % 1_000_000i64) as f64 / 1_000_000.0;
        r1 = (r1 + 1.0) / 2.0;
        r2 = (r2 + 1.0) / 2.0;
        r3 = (r3 + 1.0) / 2.0;
        let r = r1 + r2 + r3;
        if r > 0.0 {
            r1 /= r;
            r2 /= r;
            r3 /= r;
            points.push((r1, r2, r3));
        }
    }

    // Prepopulate colors with what appear to be an optimal sequence.  This is as measured by
    // distance from each other.  These colors are not optimal relative to color blindness.

    let mut fracs = vec![
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
        (0.0, 1.0 / 3.0, 2.0 / 3.0),
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.0),
        (2.0 / 3.0, 0.0, 1.0 / 3.0),
        (0.0, 2.0 / 3.0, 1.0 / 3.0),
        (1.0 / 3.0, 0.0, 2.0 / 3.0),
        (5.0 / 9.0, 2.0 / 9.0, 2.0 / 9.0),
        (4.0 / 9.0, 1.0 / 9.0, 4.0 / 9.0),
        (7.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0),
        (1.0 / 9.0, 4.0 / 9.0, 4.0 / 9.0),
        (2.0 / 9.0, 5.0 / 9.0, 2.0 / 9.0),
        (4.0 / 9.0, 4.0 / 9.0, 1.0 / 9.0),
        (2.0 / 9.0, 2.0 / 9.0, 5.0 / 9.0),
        (1.0 / 9.0, 1.0 / 9.0, 7.0 / 9.0),
        (1.0 / 9.0, 7.0 / 9.0, 1.0 / 9.0),
        (1.0 / 9.0, 2.0 / 9.0, 6.0 / 9.0),
        (3.0 / 9.0, 1.0 / 9.0, 5.0 / 9.0),
        (3.0 / 9.0, 2.0 / 9.0, 4.0 / 9.0),
        (1.0 / 9.0, 8.0 / 9.0, 0.0 / 9.0),
        (3.0 / 9.0, 4.0 / 9.0, 2.0 / 9.0),
        (5.0 / 9.0, 3.0 / 9.0, 1.0 / 9.0),
    ];
    if fracs.len() >= ngroups {
        fracs.truncate(ngroups);

    // Add more colors if needed.
    } else {
        for _ in fracs.len()..ngroups {
            let mut max_dist = 0.0;
            let mut best = (0.0, 0.0, 0.0);
            for p in points.iter() {
                let mut min_dist = 100.0_f64;
                for q in fracs.iter() {
                    let d = (p.0 - q.0) * (p.0 - q.0)
                        + (p.1 - q.1) * (p.1 - q.1)
                        + (p.2 - q.2) * (p.2 - q.2);
                    min_dist = min_dist.min(d);
                }
                if min_dist > max_dist {
                    max_dist = min_dist;
                    best = *p;
                }
            }
            fracs.push(best);
        }
    }
    for i in 0..fracs.len() {
        let r = 255 - (fracs[i].0 * sum as f64).round() as usize;
        let g = 255 - (fracs[i].1 * sum as f64).round() as usize;
        let b = 255 - (fracs[i].2 * sum as f64).round() as usize;
        group_color.push(format!("rgb({},{},{})", r, g, b));
    }
    group_color
}
