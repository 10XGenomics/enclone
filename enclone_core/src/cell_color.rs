// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Define the coloring scheme for cells.
//
// color_spec:
// - const,color                              all cells are assigned the given color (default/black)
// - iso,color1,...,colorn                    by isotype; uses default colors if none provided
// - var,name,minmax,min,max                  by values of given variable
// - catvar,varlist,n                         by categorical variables
// - dataset                                  by dataset, using the color field in META
// - origin,origin1,color1,...,originn,colorn by origin, using the given assignment
// - bc                                       by color to barcode assignment, via BC or META/bc
// - mark                                     (internal)
//
// Related, doesn't really belong here:
//
// HONEY=file:color-spec:legend-spec
//                       none
//                       origin,blue,123085,red,123089
//
// The default for :color-spec is const,black.  The :legend-spec would usually be omitted.
//
// Multiple HONEY commands can be supplied (but only one for enclone visual).

use std::collections::HashMap;

#[derive(Clone)]
pub struct ColorByIsotype {
    pub color: Vec<String>,
    pub show_legend: bool,
}

// For coloring by variable value, the value of a variable is first truncated to the range
// [min, max], then scaled to [0, 1], then converted to 0,...,255, and then converted to a color
// using the turbo color scheme.  If the value of a variable is unspecified for a given cell,
// or not convertible to a number, black is used.

#[derive(Clone)]
pub struct ColorByVariableValue {
    pub var: String,
    pub display_var: String,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

#[derive(Clone)]
pub struct ColorByCategoricalVariableValue {
    pub vars: Vec<String>,
    pub maxcat: usize,
}

#[derive(Clone)]
pub struct ColorByDataset {}

#[derive(Clone)]
pub struct ColorBySample {
    pub by_meta: bool,
    pub specification: HashMap<String, String>,
}

#[derive(Clone)]
pub struct ColorByBarcodeSpecification {}

#[derive(Clone)]
pub enum CellColor {
    Unspecified,
    ByIsotype(ColorByIsotype),
    ByVariableValue(ColorByVariableValue),
    ByCategoricalVariableValue(ColorByCategoricalVariableValue),
    BySample(ColorBySample),
    ByBarcodeSpecification(ColorByBarcodeSpecification),
    ByDataset(ColorByDataset),
}

impl Default for CellColor {
    fn default() -> Self {
        CellColor::Unspecified
    }
}
