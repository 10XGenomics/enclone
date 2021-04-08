// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// For radius r and n = 0, 1, ..., consider a counterclockwise spiral of lattice-packed disks of
// radius r, starting at the origin and going first to the right.  Return the coordinates of the
// center of the nth disk.  See this picture:
// https://www.researchgate.net/profile/Guorui_Li4/publication/220270050/figure/fig1/
//         AS:393993713143808@1470946829076/The-hexagonal-coordinate-system.png
// There is no attempt at efficiency.

pub fn hex_coord(n: usize, r: f64) -> (f64, f64) {
    // Special case.
    if n == 0 {
        return (0.0, 0.0);
    }
    // If the hexagons are numbered 0, 1, ... outward, which hexagon "hid" are we on and
    // which position "hpos" on that are we at?
    let mut hid = 1;
    let mut k = 6;
    let mut hpos = n - 1;
    loop {
        if hpos < k {
            break;
        }
        hpos -= k;
        hid += 1;
        k += 6;
    }
    // Find coordinates.
    let c = r * 3.0f64.sqrt() / 2.0; // center to center distance, divided by 2
    let mut x = hid as f64 * 2.0 * c;
    let mut y = 0.0;
    let mut p = hpos;
    if p > 0 {
        // Traverse the six faces, as far as we have to go.
        for _ in 0..hid {
            x -= c;
            y += 1.5;
            p -= 1;
            if p == 0 {
                break;
            }
        }
        if p > 0 {
            for _ in 0..hid {
                x -= 2.0 * c;
                p -= 1;
                if p == 0 {
                    break;
                }
            }
            if p > 0 {
                for _ in 0..hid {
                    x -= c;
                    y -= 1.5;
                    p -= 1;
                    if p == 0 {
                        break;
                    }
                }
                if p > 0 {
                    for _ in 0..hid {
                        x += c;
                        y -= 1.5;
                        p -= 1;
                        if p == 0 {
                            break;
                        }
                    }
                    if p > 0 {
                        for _ in 0..hid {
                            x += 2.0 * c;
                            p -= 1;
                            if p == 0 {
                                break;
                            }
                        }
                        if p > 0 {
                            for _ in 0..hid - 1 {
                                x += c;
                                y += 1.5;
                                p -= 1;
                                if p == 0 {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    x *= 2.0 / 3.0f64.sqrt();
    y *= 2.0 / 3.0f64.sqrt();
    (x, y)
}
