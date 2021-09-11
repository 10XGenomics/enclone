// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::apocalypse::*;
use crate::archive::*;
use crate::help::*;
use crate::history::*;
use crate::popover::*;
use crate::summary::*;
use crate::testsuite::*;
use crate::*;
use gui_structures::ComputeState::*;
use gui_structures::*;
use iced::Length::Units;
use iced::{
    Align, Application, Button, Clipboard, Color, Column, Command, Container, Element, Image,
    Length, Row, Rule, Scrollable, Space, Subscription, Text, TextInput,
};
// use iced::Subscription;
// use iced_native::{window, Event};
use iced_native::{event, subscription, window, Event};
use io_utils::*;
use messages::Message;
use std::env;
use std::fs::{create_dir, create_dir_all, metadata, remove_file, File};
use std::io::{Read, Write};
use std::sync::atomic::Ordering::SeqCst;
use std::thread;
use std::time::Duration;

fn handle_resize(width: u32, height: u32) -> Option<Message> {
    Some(Message::Resize(width, height))
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

impl Application for EncloneVisual {
    type Executor = iced::executor::Default;
    type Message = Message;
    type Flags = ();

    fn new(_flags: ()) -> (EncloneVisual, Command<Message>) {
        prepare_for_apocalypse_visual();
        COOKBOOK_CONTENTS.lock().unwrap().push(format_cookbook());
        let mut x = EncloneVisual::default();
        x.submit_button_text = "Submit".to_string();
        x.compute_state = WaitingForRequest;
        x.copy_image_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.snapshot_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.archive_refresh_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.copy_selected_metrics_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.cookbook = parse_cookbook();
        x.width = INITIAL_WIDTH;
        CURRENT_WIDTH.store(INITIAL_WIDTH as usize, SeqCst);
        CURRENT_WIDTH_LAST_SEEN.store(INITIAL_WIDTH as usize, SeqCst);
        x.height = INITIAL_HEIGHT;
        let mut home = String::new();
        for (key, value) in env::vars() {
            if key == "HOME" {
                home = value.clone();
            }
        }
        if home.len() == 0 {
            eprintln!(
                "Unable to determine home directory.  This is unexpected and \
                pathological.\nPlease report this problem!\n"
            );
            std::process::exit(1);
        }
        let enclone = format!("{}/enclone", home);
        if !path_exists(&enclone) {
            eprintln!(
                "Oddly, you do not have a directory ~/enclone.  Normally this would be\n\
                created by following the installation instructions at bit.ly/enclone.  Please do \
                that or at least create the directory.\n"
            );
            std::process::exit(1);
        }
        x.visual = format!("{}/visual", enclone);
        if VISUAL_DIR.lock().unwrap().len() > 0 {
            x.visual = VISUAL_DIR.lock().unwrap()[0].clone();
        }
        let history = format!("{}/history", x.visual);
        if !path_exists(&history) {
            let res = create_dir_all(&history);
            if res.is_err() {
                eprintln!(
                    "Unable to create the directory ~/enclone/visual/history.  This is odd and \
                    unexpected.\nPlease report this problem!\n"
                );
                std::process::exit(1);
            }
        }
        x.sharing_enabled = REMOTE_SHARE.lock().unwrap().len() > 0;
        x.archive_dir = Some(history.clone());

        // Read shares.  If the file is corrupted, silently ignore it.

        if x.sharing_enabled {
            let shares = format!("{}/shares", x.visual);
            if path_exists(&shares) {
                let share_size = std::fs::metadata(&shares).unwrap().len() as usize;
                let n = std::mem::size_of::<Share>();
                if share_size % n == 0 {
                    let mut bytes = Vec::<u8>::new();
                    let mut f = File::open(&shares).unwrap();
                    f.read_to_end(&mut bytes).unwrap();
                    assert_eq!(bytes.len(), share_size);
                    unsafe {
                        x.shares = bytes.align_to::<Share>().1.to_vec();
                    }
                }
            }
        }

        // Keep going.

        if x.archive_dir.is_some() {
            let arch_dir = &x.archive_dir.as_ref().unwrap();
            if path_exists(&arch_dir) {
                if metadata(&arch_dir).unwrap().is_dir() {
                    x.archive_list = dir_list(&arch_dir);
                }
            }
        }
        x.archive_list.reverse();
        let n = x.archive_list.len();
        x.restore_requested = vec![false; n];
        x.delete_requested = vec![false; n];
        x.deleted = vec![false; n];
        x.expand_archive_entry = vec![false; n];
        x.restore_msg = vec![String::new(); n];
        x.archived_command_list = vec![None; n];
        x.archive_name = vec![iced::text_input::State::default(); n];
        x.archive_name_value = vec![String::new(); n];
        x.archive_name_change_button_color = vec![Color::from_rgb(0.0, 0.0, 0.0); n];
        x.copy_archive_narrative_button_color = vec![Color::from_rgb(0.0, 0.0, 0.0); n];
        x.archive_name_change_button = vec![iced::button::State::default(); n];
        x.archive_narrative_button = vec![iced::button::State::default(); n];
        x.copy_archive_narrative_button = vec![iced::button::State::default(); n];
        x.archive_share_requested = vec![false; n];
        x.archive_origin = vec![String::new(); n];
        x.archive_narrative = vec![String::new(); n];
        x.copy_narrative_button_color = Color::from_rgb(0.0, 0.0, 0.0);
        x.copy_summary_button_color = Color::from_rgb(0.0, 0.0, 0.0);

        // Fetch cookbooks.  For remote cookbooks, we keep a local copy, which may or may not
        // be current.

        let cookbook_dir = include_dir::include_dir!("src/cookbooks");
        for cookbook in cookbook_dir.find("*.cb").unwrap() {
            let f = cookbook_dir.get_file(cookbook.path()).unwrap();
            x.cookbooks.push(f.contents().to_vec());
        }
        GET_MY_COOKBOOKS.store(true, SeqCst);
        while GET_MY_COOKBOOKS.load(SeqCst) {
            thread::sleep(Duration::from_millis(10));
        }
        if !META_TESTING.load(SeqCst) && !TEST_MODE.load(SeqCst) {
            let local_remote = format!("{}/remote_cookbooks", x.visual);
            if REMOTE.load(SeqCst) {
                if path_exists(&local_remote) {
                    let list = dir_list(&local_remote);
                    for f in list.iter() {
                        remove_file(&format!("{}/{}", local_remote, f)).unwrap();
                    }
                } else {
                    create_dir(&local_remote).unwrap();
                }
                let n = REMOTE_COOKBOOKS.lock().unwrap().len();
                for i in 0..n {
                    let content = &REMOTE_COOKBOOKS.lock().unwrap()[i];
                    x.cookbooks.push(content.clone());
                    let f = format!("{}/{}", local_remote, i + 1);
                    let mut file = File::create(&f).unwrap();
                    file.write_all(&content).unwrap();
                }
            } else if path_exists(&local_remote) {
                let list = dir_list(&local_remote);
                for f in list.iter() {
                    let mut h = File::open(&format!("{}/{}", local_remote, f)).unwrap();
                    let mut content = Vec::<u8>::new();
                    h.read_to_end(&mut content).unwrap();
                    x.cookbooks.push(content);
                }
            }
        }
        let nc = x.cookbooks.len();
        x.expand_cookbook_entry = vec![false; nc];
        x.restore_cookbook_requested = vec![false; nc];
        x.cookbook_narrative_button = vec![iced::button::State::default(); nc];
        x.copy_cookbook_narrative_button = vec![iced::button::State::default(); nc];
        x.copy_cookbook_narrative_button_color = vec![Color::from_rgb(0.0, 0.0, 0.0); nc];
        x.restore_cookbook_msg = vec![String::new(); nc];
        for i in 0..nc {
            let evh = EncloneVisualHistory::restore_from_bytes(&x.cookbooks[i]).unwrap();
            x.cookbook_name.push(evh.name_value.clone());
            let mut cc = Vec::<String>::new();
            for j in 0..evh.translated_input_history.len() {
                cc.push(
                    evh.translated_input_hist_uniq[evh.translated_input_history[j] as usize]
                        .clone(),
                );
            }
            x.cookbook_command_list.push(Some(cc));
            x.cookbook_narrative.push(evh.narrative.clone());
        }

        // Handle test and meta modes.

        if !TEST_MODE.load(SeqCst) && !META_TESTING.load(SeqCst) {
            (x, Command::none())
        } else if !META_TESTING.load(SeqCst) {
            thread::sleep(Duration::from_millis(1000));
            (x, Command::perform(noop(), Message::RunTests))
        } else {
            let id = META.load(SeqCst); // id of meta test
            x.this_meta = metatests()[id].clone();
            x.meta_pos = 0;
            (x, Command::perform(noop0(), Message::Meta))
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    fn title(&self) -> String {
        String::from("EncloneVisual")
    }

    fn update(&mut self, message: Message, clipboard: &mut Clipboard) -> Command<Message> {
        self.process_message(message, clipboard)
    }

    /*
    fn should_exit(&self) -> bool {
        self.should_exit
    }

    fn subscription(&self) -> Subscription<Message> {
        iced_native::subscription::events().map(Message::EventOccurred)
    }
    */

    // The subscription detects window resize events and uses those to reset
    // self.width and self.height.

    fn subscription(&self) -> Subscription<Message> {
        subscription::events_with(|event, status| {
            if let event::Status::Captured = status {
                return None;
            }
            match event {
                Event::Window(window::Event::Resized { width, height }) => {
                    handle_resize(width, height)
                }
                _ => None,
            }
        })
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    fn view(&mut self) -> Element<Message> {
        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Handle popovers.

        if self.summary_mode {
            return summary(self);
        }
        if self.console_mode {
            return console(self);
        }
        if self.cookbook_mode {
            return cookbook(self);
        }
        if self.archive_mode {
            return archive(self);
        }
        if self.help_mode {
            return help(self);
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Now do everything else.

        let text_input1 = TextInput::new(
            &mut self.input1,
            "",
            &self.input1_value,
            Message::InputChanged1,
        )
        .padding(7)
        .font(DEJAVU_BOLD)
        .size(16);
        let text_input2 = TextInput::new(
            &mut self.input2,
            "",
            &self.input2_value,
            Message::InputChanged2,
        )
        .padding(7)
        .font(DEJAVU_BOLD)
        .size(16);
        let text_input_column = Column::new()
            .spacing(5)
            .width(iced::Length::Fill)
            .push(text_input1)
            .push(text_input2);

        let button = Button::new(
            &mut self.button,
            Text::new(if self.compute_state == WaitingForRequest {
                "Submit"
            } else {
                "thinking"
            }),
        )
        .padding(10)
        .on_press(Message::SubmitButtonPressed(Ok(())));
        let clear_button = Button::new(&mut self.clear_button, Text::new("Clear"))
            .padding(10)
            .on_press(Message::ClearButtonPressed);

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Define the button complex that is the "control panel".

        let mut blank = false;
        if self.h.history_index > 0 {
            blank = self.h.is_blank[self.h.history_index as usize - 1];
        }
        let command_complex_height;
        let mut command_complex = Row::new().spacing(10);
        {
            const FB_BUTTON_FONT_SIZE: u16 = 45;
            let back_button = Button::new(
                &mut self.back_button,
                Text::new("⇧").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::BackButtonPressed(Ok(())));

            let forward_button = Button::new(
                &mut self.forward_button,
                Text::new("⇩").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::ForwardButtonPressed(Ok(())));

            const COPY_BUTTON_FONT_SIZE: u16 = 15;
            let del_button = Button::new(
                &mut self.del_button,
                Text::new("Del ").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DelButtonPressed(Ok(())));

            let null_button1 = Button::new(
                &mut self.null_button1,
                Text::new(" ").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);

            let null_button2 = Button::new(
                &mut self.null_button2,
                Text::new(" ").font(DEJAVU_BOLD).size(FB_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);
            let copy_image_button = Button::new(
                &mut self.copy_image_button,
                Text::new("Copy image")
                    .size(COPY_BUTTON_FONT_SIZE)
                    .color(self.copy_image_button_color),
            )
            .on_press(Message::GraphicsCopyButtonPressed);
            let null_copy_image_button = Button::new(
                &mut self.null_button3,
                Text::new("          ").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::GraphicsCopyButtonPressed);
            let mut state_pos = format!("{}", self.h.history_index);
            if state_pos.len() == 1 {
                state_pos += "    ";
            } else if state_pos.len() == 2 {
                state_pos += "  ";
            }
            let state_pos_button = Button::new(
                &mut self.state_pos_button_null,
                Text::new(&state_pos).size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::DoNothing);

            let mut button_column2 = Column::new().spacing(8);
            button_column2 = button_column2.push(state_pos_button);
            if self.h.history_index > 1 {
                button_column2 = button_column2.push(back_button);
            } else {
                button_column2 = button_column2.push(null_button1);
            }
            if (self.h.history_index as usize) < self.h.svg_history.len() {
                button_column2 = button_column2.push(forward_button);
            } else {
                button_column2 = button_column2.push(null_button2);
            }
            button_column2 = button_column2.push(del_button);

            // Add command box.

            const MAX_LINE: usize = 35;
            let mut log_lines = 1;
            let mut log = String::new();
            if self.h.history_index >= 1 {
                let cmd = self.h.translated_input_hist_uniq
                    [self.h.translated_input_history[self.h.history_index as usize - 1] as usize]
                    .clone();
                let mut rows = Vec::<Vec<String>>::new();
                let folds = fold(&cmd, MAX_LINE);
                log_lines = folds.len();
                for i in 0..folds.len() {
                    rows.push(vec![folds[i].clone()]);
                }
                for i in 0..rows.len() {
                    if i > 0 {
                        log += "\n";
                    }
                    log += &mut rows[i][0].clone();
                }
            }

            // Create summary button.

            let summary_button = Button::new(
                &mut self.summary_button,
                Text::new("Summary").size(COPY_BUTTON_FONT_SIZE),
            )
            .on_press(Message::SummaryOpen(Ok(())));

            // Create narrative button.

            let narrative_width;
            if !blank {
                narrative_width = MAX_LINE;
            } else {
                narrative_width = 80;
            }
            let mut logx = String::new();
            let mut logx_lines = 1;
            let mut have_narrative = false;
            if self.h.history_index >= 1 {
                let mut cmd = self.h.narrative_hist_uniq
                    [self.h.narrative_history[self.h.history_index as usize - 1] as usize]
                    .clone();
                if cmd.len() == 0 {
                    cmd = "Narrative: click to paste in clipboard".to_string();
                } else {
                    have_narrative = true;
                }
                let mut rows = Vec::<Vec<String>>::new();
                let folds = fold(&cmd, narrative_width);
                logx_lines = folds.len();
                for i in 0..folds.len() {
                    rows.push(vec![folds[i].clone()]);
                }
                for i in 0..rows.len() {
                    if i > 0 {
                        logx += "\n";
                    }
                    logx += &mut rows[i][0].clone();
                }
            }
            let narrative_button = Button::new(
                &mut self.narrative_button,
                Text::new(&logx)
                    .font(DEJAVU_BOLD)
                    .size(12)
                    .color(Color::from_rgb(1.0, 0.0, 0.5)),
            )
            .on_press(Message::Narrative);

            let copy_narrative_button = Button::new(
                &mut self.copy_narrative_button,
                Text::new("Copy")
                    .size(COPY_BUTTON_FONT_SIZE)
                    .color(self.copy_narrative_button_color),
            )
            .on_press(Message::CopyNarrative);

            // Build the command column.

            let mut row = Row::new().spacing(8);
            if self.h.history_index >= 1 && !self.h.is_blank[self.h.history_index as usize - 1] {
                row = row.push(copy_image_button);
            } else {
                row = row.push(null_copy_image_button);
            }
            row = row.push(
                Button::new(
                    &mut self.command_copy_button,
                    Text::new("Copy command").size(COPY_BUTTON_FONT_SIZE),
                )
                .on_press(Message::CommandCopyButtonPressed),
            );
            let mut col = Column::new().spacing(8).align_items(Align::End);
            const SMALL_FONT: u16 = 12;
            command_complex_height = ((1 + log_lines + logx_lines) * SMALL_FONT as usize)
                + (4 * 8)
                + (2 * COPY_BUTTON_FONT_SIZE as usize);
            col = col.push(
                Button::new(
                    &mut self.null_button,
                    Text::new(&log).font(DEJAVU_BOLD).size(SMALL_FONT),
                )
                .on_press(Message::DoNothing),
            );
            col = col.push(row);
            col = col.push(summary_button);
            col = col.push(narrative_button);
            if have_narrative {
                col = col.push(copy_narrative_button);
            }

            // Add the command column to the row.

            command_complex = command_complex.push(col);

            // Add up and down arrows.

            command_complex = command_complex.push(button_column2);
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Build the scrollable for clonotypes.  We truncate lines to prevent wrapping.

        const CLONOTYPE_FONT_SIZE: u16 = 13;
        let font_width = CLONOTYPE_FONT_SIZE as f32 * DEJAVU_WIDTH_OVER_HEIGHT;
        let available = self.width - (3 * SPACING + SCROLLBAR_WIDTH) as u32;
        let nchars = (available as f32 / font_width).round() as usize;
        let mut trunc = String::new();
        let failed = self.output_value.contains("enclone failed");
        for line in self.output_value.lines() {
            for (i, c) in line.chars().enumerate() {
                if i == nchars && !failed {
                    break;
                }
                trunc.push(c);
            }
            trunc.push('\n');
        }
        let scrollable = Scrollable::new(&mut self.scroll)
            .width(Length::Fill)
            .height(Length::Units(100))
            .scrollbar_width(SCROLLBAR_WIDTH)
            .scroller_width(12)
            .style(style::ScrollableStyle)
            .push(
                Text::new(&trunc)
                    .font(DEJAVU_BOLD)
                    .size(CLONOTYPE_FONT_SIZE),
            );

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Fix the height of the SVG.  This needs to be set so that there is enough room for
        // the clonotype tables.  We do not set the width because it's the height that we need
        // to control.

        let mut svg_height = if !blank { SVG_HEIGHT } else { SVG_NULL_HEIGHT };
        // 60 is a fudge factor:
        svg_height = std::cmp::max(svg_height, command_complex_height as u16 + 60);

        // Display the SVG.

        let png_banner = include_bytes!("../../img/enclone_banner.png").to_vec();
        let banner = Image::new(iced::image::Handle::from_memory(png_banner)).width(Units(500));

        let have_canvas = self.canvas_view.state.geometry_value.is_some();
        let mut graphic_row = Row::new().spacing(10);
        if self.svg_value.len() > 0 {
            // Show the graphic.

            if have_canvas {
                graphic_row = graphic_row
                    .push(
                        self.canvas_view
                            .view()
                            .map(move |message| Message::GroupClicked(message)),
                    )
                    .height(Units(svg_height));
            } else {
                let svg_as_png =
                    Image::new(iced::image::Handle::from_memory(self.png_value.clone()))
                        .height(Units(svg_height));
                graphic_row = graphic_row.push(svg_as_png);
            }

            // Insert space to push the graphic to the left and the command column to the right.

            if !have_canvas {
                graphic_row = graphic_row.push(Space::with_width(Length::Fill));
            }
            graphic_row = graphic_row.push(command_complex);
        }

        // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

        // Put it all together.

        let left_buttons = Column::new()
            .spacing(8)
            .push(Button::new(&mut self.exit_state, Text::new("Exit")).on_press(Message::Exit))
            .push(
                Button::new(&mut self.open_state, Text::new("Help"))
                    .on_press(Message::HelpOpen(Ok(()))),
            )
            .push(
                Button::new(&mut self.open_state_cookbook, Text::new("Cookbook"))
                    .on_press(Message::CookbookOpen),
            );
        let console_button = Button::new(&mut self.console_open_button, Text::new("Console"))
            .on_press(Message::ConsoleOpen);
        let snapshot_button = Button::new(
            &mut self.snapshot_button,
            Text::new("Snapshot").color(self.snapshot_button_color),
        )
        .on_press(Message::Snapshot);
        let console_row = Row::new()
            .spacing(8)
            .push(snapshot_button)
            .push(console_button);
        let mut save_text = Text::new("Save");
        if self.save_in_progress {
            save_text = save_text.color(Color::from_rgb(1.0, 0.0, 0.0));
        }
        let save_button = Button::new(&mut self.save_button, save_text).on_press(Message::Save);
        let mut save_on_exit_text = Text::new("On Exit").width(Units(66));
        if self.save_on_exit {
            save_on_exit_text = save_on_exit_text.color(Color::from_rgb(1.0, 0.0, 0.0));
        }
        let save_on_exit_button = Button::new(&mut self.save_on_exit_button, save_on_exit_text)
            .on_press(Message::SaveOnExit);
        let save_row = Row::new()
            .spacing(8)
            .push(save_button)
            .push(save_on_exit_button);
        let archive_button = Button::new(
            &mut self.archive_open_button,
            Text::new("Archive").width(Units(66)),
        )
        .on_press(Message::ArchiveOpen(Ok(())));
        let mut top_row = Row::new()
            .align_items(Align::Center)
            .push(left_buttons)
            .push(Space::with_width(Length::Fill))
            .push(banner)
            .push(Space::with_width(Length::Fill));
        let right_col = Column::new()
            .align_items(Align::End)
            .spacing(8)
            .push(console_row)
            .push(save_row)
            .push(archive_button);
        top_row = top_row.push(right_col);
        let mut content = Column::new()
            .spacing(SPACING)
            .padding(20)
            .push(top_row)
            .push(
                Row::new()
                    .spacing(10)
                    .align_items(Align::Center)
                    .push(text_input_column)
                    .push(button)
                    .push(clear_button),
            )
            // .push(Row::new().spacing(10).push(svg))
            .push(graphic_row);
        if self.h.svg_history.len() > 0 {
            content = content.push(Rule::horizontal(10).style(style::RuleStyle));
        }
        content = content.push(
            Row::new()
                .height(Length::Units(1000)) // Height of scrollable window, maybe??
                .align_items(Align::Center)
                .push(scrollable),
        );
        Container::new(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}
