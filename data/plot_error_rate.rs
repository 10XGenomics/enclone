
use enclone_tail::plot_points::plot_points;

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

// Temporarily changed plot_points.rs to 
// * set xlow to 0 and xhigh to 1500000
// * set ylow to 0 and yhigh to 1.4.
// * set title to "clonotyping accuracy"
// * increase margin from 25 to 35
// *     let point_color = RGBColor(255,153,51);
// * added:
/*
    let point_color2 = BLUE;
    let points2 = vec![
        (652537.0, 0.66), 
        (672653.0, 1.23),
        (405168.0, 0.75),
        (245629.0, 0.43),
        (301997.0, 0.68),
    ];
    chart
        .draw_series(PointSeries::of_element(
            points2.clone(),
            point_size,
            &point_color2,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();
    let point_color3 = RED;
    let points3 = vec![(1325190.0, 0.90)];
    chart
        .draw_series(PointSeries::of_element(
            points3.clone(),
            point_size,
            &point_color3,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();
*/

// ================================================================================================

fn main() {

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

    let cells = 1325190.0 as f32;

    let mut points = Vec::<(f32, f32)>::new();
    for i in 0..sr.len() {
        points.push((sr[i].0 as f32 * cells, sr[i].1));
    }

    let mut svg = String::new();

    plot_points(
        &points, 
        "number of cells",
        "p(two unrelated cells are co-clonotyped) x 10^9",
        &mut svg,
        false
    ).unwrap();

    print!("{}", svg);

}
