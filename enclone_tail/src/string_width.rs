// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use vector_utils::*;

// Estimate the width in pixels of an Arial string at a given font size.
// This uses a hardcoded table of widths of Arial 1000 point characters.  The table
// is incomplete, and we use a fixed value for all other characters.

pub fn arial_width(s: &str, font_size: usize) -> f64 {
    const DEFAULT_WIDTH: usize = 1000;
    let mut len = 0;
    for c in s.chars() {
        let p = bin_position1_2(&ARIAL_1000_WIDTH_TABLE, &c);
        if p < 0 {
            len += DEFAULT_WIDTH;
        } else {
            len += ARIAL_1000_WIDTH_TABLE[p as usize].1;
        }
    }
    len as f64 * font_size as f64 / 1000.0
}

// Note that the following table must be sorted.

const ARIAL_1000_WIDTH_TABLE: [(char, usize); 66] = [
    (' ', 278),
    ('-', 334),
    ('.', 278),
    ('0', 557),
    ('1', 557),
    ('2', 557),
    ('3', 557),
    ('4', 557),
    ('5', 557),
    ('6', 557),
    ('7', 557),
    ('8', 557),
    ('9', 557),
    ('A', 667),
    ('B', 667),
    ('C', 723),
    ('D', 723),
    ('E', 667),
    ('F', 611),
    ('G', 778),
    ('H', 723),
    ('I', 278),
    ('J', 500),
    ('K', 667),
    ('L', 557),
    ('M', 834),
    ('N', 723),
    ('O', 778),
    ('P', 667),
    ('Q', 778),
    ('R', 723),
    ('S', 667),
    ('T', 611),
    ('U', 723),
    ('V', 667),
    ('W', 944),
    ('X', 667),
    ('Y', 667),
    ('Z', 611),
    ('_', 557),
    ('a', 557),
    ('b', 557),
    ('c', 500),
    ('d', 557),
    ('e', 557),
    ('f', 278),
    ('g', 557),
    ('h', 557),
    ('i', 223),
    ('j', 223),
    ('k', 500),
    ('l', 223),
    ('m', 834),
    ('n', 557),
    ('o', 557),
    ('p', 557),
    ('q', 557),
    ('r', 334),
    ('s', 500),
    ('t', 278),
    ('u', 557),
    ('v', 500),
    ('w', 723),
    ('x', 500),
    ('y', 500),
    ('z', 500),
];

// Computed using this html code from the internet.  For each character,
// edit the code, refresh the html window, and click the button.
// Obviously one write write a code that did all characters and generated
// the rust code directly with one button push.

/*
<!DOCTYPE html>
<html>
<head>
    <title>Calculate the text width with JavaScript</title>
</head>
<body>
    <p><span class="output"></span></p>
    <button onclick="getTextWidth()">Calculate text width</button>
    <script type="text/javascript">
        function getTextWidth() {

            inputText = "A";
            font = "1000px arial";

            canvas = document.createElement("canvas");
            context = canvas.getContext("2d");
            context.font = font;
            width = context.measureText(inputText).width;
            formattedWidth = Math.ceil(width) + "px";
            document.querySelector('.output').textContent
                        = formattedWidth;
        }
    </script>
</body>
</html>
*/
