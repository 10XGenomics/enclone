// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is the start of a variable registry for enclone.
//
// The intention is that all variables would be registered here.  Or perhaps we could instead
// have a distributed system.
//
// String fields are to be replaced by structures.

pub struct Variable {
    pub name: String,
    pub scope: String,
    pub prereqs: Vec<String>,
    pub value_type: String,
    pub function: String,
    pub doc: String,
}

#[rustfmt::skip]
pub fn variable_registry() -> Vec<Variable> {
    let mut reg = Vec::<Variable>::new();
    
    // <FeatureName>_cellular_u

    reg.push( Variable {
        name:       "<FeatureName>_cellular_u".to_string(),
        scope:      "dataset".to_string(),
        prereqs:    vec!["per_feature_metrics.csv".to_string()],
        value_type: "float[0,100].precision(1)".to_string(),
        function:   "in ***.rs".to_string(),
        doc:        "For a given feature, the percent of UMIs that are identified by the \
                     cellranger pipeline as lying in a cell.".to_string(),
    });
    
    // <FeatureName>_cellular_r

    reg.push( Variable {
        name:       "<FeatureName>_cellular_r".to_string(),
        scope:      "dataset".to_string(),
        prereqs:    vec!["per_feature_metrics.csv".to_string()],
        value_type: "float[0,100].precision(1)".to_string(),
        function:   "in ***.rs".to_string(),
        doc:        "For a given feature, the percent of reads that are identified by the \
                     cellranger pipeline as lying in a cell.".to_string(),
    });

    reg
}
