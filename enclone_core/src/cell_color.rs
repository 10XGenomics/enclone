// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Define the coloring scheme for cells.
//
// CELL_COLOR
// =isotype,color1,...,colorn                   by isotype; uses default colors if none provided
// =var,name,min,max                            by values of given variable
// =sample_from_meta                            by sample, using the color field in META
// =sample,sample1->color1,...,samplen->colorn  by sample, using the given assignment
// =bc                                          by color to barcode assignment, via BC or META/bc
//
// For the default (Unspecified) schema, all cells are black.

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
    pub min: Option<f64>,
    pub max: Option<f64>,
}

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
    BySample(ColorBySample),
    ByBarcodeSpecification(ColorByBarcodeSpecification),
}

impl Default for CellColor {
    fn default() -> Self {
        CellColor::Unspecified
    }
}
