// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Print dataset table.
//
// dataset_table

use enclone_core::test_def::test_donor_id;
use pretty_trace::PrettyTrace;

const NAIVE: [usize; 40] = [
    1279049, 1279053, 1279057, 1279061, 1279065, 1279069, 1279073, 1279077, 1287144, 1287145,
    1287146, 1287147, 1287152, 1287153, 1287154, 1287155, 1287160, 1287161, 1287162, 1287163,
    1287168, 1287169, 1287170, 1287171, 1287176, 1287177, 1287178, 1287179, 1287184, 1287185,
    1287186, 1287187, 1287192, 1287193, 1287194, 1287195, 1287200, 1287201, 1287202, 1287203,
];
const PLASMABLAST: [usize; 6] = [1279052, 1279060, 1279068, 1279072, 1279076, 1279080];
const SWITCHED: [usize; 24] = [
    1279051, 1279055, 1279059, 1279063, 1279067, 1279071, 1279075, 1279079, 1287150, 1287151,
    1287158, 1287159, 1287166, 1287167, 1287174, 1287175, 1287182, 1287183, 1287190, 1287191,
    1287198, 1287199, 1287206, 1287207,
];
const UNSWITCHED: [usize; 24] = [
    1279050, 1279054, 1279058, 1279062, 1279066, 1279070, 1279074, 1279078, 1287148, 1287149,
    1287156, 1287157, 1287164, 1287165, 1287172, 1287173, 1287180, 1287181, 1287188, 1287189,
    1287196, 1287197, 1287204, 1287205,
];

fn main() {
    PrettyTrace::new().on();
    let mut all = Vec::<usize>::new();
    all.append(&mut NAIVE.to_vec());
    all.append(&mut UNSWITCHED.to_vec());
    all.append(&mut SWITCHED.to_vec());
    all.append(&mut PLASMABLAST.to_vec());
    all.sort();
    println!("dataset,donor,flow_class");
    for i in 0..all.len() {
        let dataset = all[i];
        for pass in 1..=2 {
            if pass == 1 && dataset.to_string().starts_with("128") {
                continue;
            }
            if pass == 2 && dataset.to_string().starts_with("127") {
                continue;
            }
            let donor_id = test_donor_id(dataset);
            if NAIVE.contains(&dataset) {
                println!("{dataset},{donor_id},naive");
            } else if UNSWITCHED.contains(&dataset) {
                if pass == 1 || donor_id == 1 {
                    println!("{dataset},{donor_id},unswitched");
                } else {
                    println!("{dataset},{donor_id},unswitched + naive");
                }
            } else if SWITCHED.contains(&dataset) {
                if pass == 1 {
                    println!("{dataset},{donor_id},switched");
                } else {
                    println!("{dataset},{donor_id},switched + naive");
                }
            } else if PLASMABLAST.contains(&dataset) {
                println!("{dataset},{donor_id},plasmablast");
            }
        }
    }
}
