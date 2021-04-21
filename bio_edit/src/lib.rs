// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

extern crate approx;

extern crate custom_derive;

#[macro_use]
extern crate lazy_static;

extern crate newtype_derive;

#[macro_use]
extern crate serde_derive;

extern crate strum_macros;

extern crate getset;

#[cfg(feature = "phylogeny")]
#[macro_use]
extern crate pest_derive;

pub mod alignment;
pub mod alphabets;
pub mod utils;
