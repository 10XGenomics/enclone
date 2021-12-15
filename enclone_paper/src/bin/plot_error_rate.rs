// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_paper::plot_points_gen::plot_points_gen;

// ================================================================================================

// LEGEND
// Clonotyping accuracy.  The probability that the clonotyping algorithm places two unrelated
// cells in the same clonotype is roughly 10^-9.  This probability varies stochastically as input
// data vary and increases with the number of cells.  Our observed value is 0.9 x 10^-9
// for 1.3M cells.
// * Two B cells are called unrelated if they arose from different fully recombined ancestors.
// * The probability that two unrelated cells are co-clonotyped
//   was estimated by clonotyping a combined dataset, containing cells from multiple
//   individuals, and determining the probability that two cells from different individuals are
//   co-clonotyped.
// * Red point: combined dataset consisting of 1,325,190 cells from 38 individuals.
// * Orange points: combined dataset was cellwise subsampled at 10%, 20%, ..., 90%,
//                  mean of 20 replicates is shown.
// * Blue points: each point represents all data from some of the individuals.

// ================================================================================================

// Special plot_points.rs ==> copy file into place.

// ================================================================================================

fn main() {

    // Define constants.


    // Set up.

    let mut points = Vec::<(u32, (u8, u8, u8), f32, f32)>::new();
    let large = 5;
    let small = 2;

    // Plot individual data points for cellwise subsampling.

    let orange = (255, 153, 51);
    let cells = 1325190.0 as f32;
    #[rustfmt::skip]
    let srx = [
        (0.1, [1.77, 0.16, 0.48, 0.00, 0.00, 0.32, 0.16, 0.32, 0.00, 0.48, 
               0.00, 0.32, 0.48, 0.16, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00]),
        (0.2, [0.29, 0.74, 0.12, 0.73, 0.37, 0.41, 0.37, 0.45, 0.24, 0.41, 
               0.85, 0.45, 0.36, 0.45, 0.08, 0.28, 0.61, 0.41, 0.24, 0.61]),
        (0.3, [0.71, 0.46, 0.68, 0.59, 0.64, 0.42, 0.18, 0.29, 0.57, 0.29, 
               0.11, 0.71, 0.42, 0.82, 0.49, 0.51, 0.18, 0.33, 0.46, 0.15]),
        (0.4, [1.03, 0.92, 0.31, 0.30, 0.46, 0.48, 0.41, 0.41, 0.56, 0.87, 
               0.59, 0.54, 0.41, 0.27, 0.52, 0.46, 1.26, 0.80, 0.50, 0.80]),
        (0.5, [0.44, 0.63, 0.67, 0.21, 0.64, 0.74, 0.39, 0.60, 0.36, 0.36, 
               0.46, 0.84, 0.85, 0.45, 0.52, 0.69, 0.62, 0.79, 0.48, 0.67]),
        (0.6, [0.97, 0.48, 0.86, 0.78, 0.60, 0.76, 0.90, 0.45, 0.59, 0.70, 
               0.74, 0.85, 0.78, 0.52, 0.84, 0.59, 0.74, 0.82, 0.73, 0.60]),
        (0.7, [0.87, 0.84, 0.61, 0.83, 0.71, 0.99, 0.80, 0.70, 0.93, 0.57, 
               1.14, 0.40, 1.02, 0.80, 0.84, 0.95, 0.83, 0.60, 0.90, 0.77]),
        (0.8, [1.13, 0.83, 0.76, 0.88, 0.74, 0.83, 0.83, 0.82, 0.84, 0.89, 
               0.74, 0.95, 0.83, 0.93, 0.75, 0.88, 0.80, 0.64, 0.78, 0.76]),
        (0.9, [0.93, 0.82, 0.97, 0.79, 0.86, 0.97, 0.85, 0.96, 0.99, 0.77, 
               0.70, 0.66, 0.88, 0.76, 0.97, 0.91, 0.87, 0.86, 0.80, 0.88]),
    ];
    for i in 0..srx.len() {
        for j in 0..srx[i].1.len() {
            points.push((small, orange, srx[i].0 as f32 * cells, srx[i].1[j]));
        }
    }

    // Plot points for dataset-level sampling.

    let blue = (0, 0, 255);
    for i in 0..srx.len() {
        let x = srx[i].0 as f32 * cells;
        let mut y = 0.0;
        for j in 0..srx[i].1.len() {
            y += srx[i].1[j];
        }
        y /= srx[i].1.len() as f32;
        points.push((large, blue, x, y));
    }

    // Plot point for all the data.

    let red = (255, 0, 0);
    let points3 = vec![(1325190.0, 0.90)];
    for i in 0..points3.len() {
        points.push((large, red, points3[i].0, points3[i].1));
    }

    // Plot mean points for cellwise sampling.

    let sr = [
        (0.1, 0.25),
        (0.2, 0.42),
        (0.3, 0.45),
        (0.4, 0.60),
        (0.5, 0.57),
        (0.6, 0.72),
        (0.7, 0.81),
        (0.8, 0.76), // preliminary
    ];
    let mut points4 = Vec::<(f32, f32)>::new();
    for i in 0..sr.len() {
        points4.push((sr[i].0 as f32 * cells, sr[i].1));
    }
    for i in 0..points4.len() {
        points.push((large, orange, points4[i].0, points4[i].1));
    }

    // Make the plot.

    let mut svg = String::new();
    plot_points_gen(
        &points,
        "number of cells",
        "p(two unrelated cells are co-clonotyped) x 10^9",
        &mut svg,
        false,
        Some("clonotyping accuracy".to_string()),
        Some(0.0),
        Some(1500000.0),
        Some(0.0),
        Some(2.0),
        Some(35),
    )
    .unwrap();
    print!("{}", svg);
}
