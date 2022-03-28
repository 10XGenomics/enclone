// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Define enclone paper test sets.

use crate::expand_integer_ranges;
use string_utils::TextUtils;

const TEST1: &str = "1279053,1279061,1287192-1287195,1287200-1287203:\
        1279050,1279058,1287196-1287197,1287204-1287205:\
        1279051,1279059,1287198-1287199,1287206-1287207:1279052,1279060";

const TEST2: &str = "1279049,1279057,1287176-1287179,1287184-1287187:\
        1279054,1279062,1287180-1287181,1287188-1287189:\
        1279055,1279063,1287182-1287183,1287190-1287191";

const TEST3: &str = "1279065,1279073,1287144-1287147,1287152-1287155:\
        1279066,1279074,1287156-1287157,1287148-1287149:\
        1279067,1279075,1287150-1287151,1287158-1287159:1279068,1279076";

const TEST4: &str = "1279069,1279077,1287160-1287163,1287168-1287171:\
        1279070,1279078,1287164-1287165,1287172-1287173:\
        1279071,1279079,1287166-1287167,1287174-1287175:1279072,1279080";

pub fn replace_at_test(x: &mut String) {
    *x = x.replace("@test1", &TEST1);
    *x = x.replace("@test2", &TEST2);
    *x = x.replace("@test3", &TEST3);
    *x = x.replace("@test4", &TEST4);
    *x = x.replace("@test", &format!("{};{};{};{}", TEST1, TEST2, TEST3, TEST4));
    *x = x.replace("@training", "1-3,5-9,11-12,14-16,18-43");
}

pub fn test_donor_id(x: usize) -> usize {
    let test1 = expand_integer_ranges(&TEST1.replace(":", ","));
    let test1 = test1.split(',').collect::<Vec<&str>>();
    for t in test1.iter() {
        if t.force_usize() == x {
            return 1;
        }
    }
    let test2 = expand_integer_ranges(&TEST2.replace(":", ","));
    let test2 = test2.split(',').collect::<Vec<&str>>();
    for t in test2.iter() {
        if t.force_usize() == x {
            return 2;
        }
    }
    let test3 = expand_integer_ranges(&TEST3.replace(":", ","));
    let test3 = test3.split(',').collect::<Vec<&str>>();
    for t in test3.iter() {
        if t.force_usize() == x {
            return 3;
        }
    }
    let test4 = expand_integer_ranges(&TEST4.replace(":", ","));
    let test4 = test4.split(',').collect::<Vec<&str>>();
    for t in test4.iter() {
        if t.force_usize() == x {
            return 4;
        }
    }
    panic!("unknown test donor id");
}
