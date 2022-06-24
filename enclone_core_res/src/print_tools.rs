// Extract the @font-face content from the current css file.

pub fn font_face_in_css() -> String {
    let f = include_str!["enclone_css_v2.css"];
    let mut x = String::new();
    let mut in_font_face = false;
    let mut count = 0;
    for line in f.lines() {
        if line.starts_with("@font-face") {
            in_font_face = true;
            count += 1;
        }
        if in_font_face {
            x += &format!("{}\n", line);
        }
        if line == "}" {
            in_font_face = false;
        }
    }
    assert_eq!(count, 2); // because there are two fonts: regular and bold
    x
}

//
