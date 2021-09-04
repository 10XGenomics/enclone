// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Define the coloring scheme for cells.
//
// CELL_COLOR
// =isotype,color1,...,colorn                   by isotype; uses default colors if none provided
// =variable,name,min,max                       by values of given variable
// =sample_from_meta                            by sample, using the color field in META
// =sample,sample1->color1,...,samplen->colorn  by sample, using the given assignment
// =bc                                          by color to barcode assignment, via BC or META/bc

use std::collections::HashMap;

#[derive(Clone)]
pub struct ColorByIsotype {
    pub color: Vec<String>,
    pub show_legend: bool,
}

#[derive(Clone)]
pub struct ColorByVariableValue {
    pub var: String,
    pub min: f64,
    pub max: f64,
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
