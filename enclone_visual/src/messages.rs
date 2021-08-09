// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#[derive(Debug, Clone)]
pub enum Message {
    InputChanged1(String),
    InputChanged2(String),
    SubmitButtonPressed(Result<(), String>),
    BackButtonPressed(Result<(), String>),
    ForwardButtonPressed(Result<(), String>),
    DelButtonPressed(Result<(), String>),
    ComputationDone(Result<(), String>),
    // EventOccurred(iced_native::Event),
    GraphicsCopyButtonPressed,
    GraphicsCopyButtonFlashed(Result<(), String>),
    CommandCopyButtonPressed,
    DoNothing,
    Exit,
    ClearButtonPressed,
    RunTests(Result<(), String>),
    Capture(Result<(), String>),
    GroupClicked(crate::canvas_view::Message),
    Resize(u32, u32),
    HelpOpen(Result<(), String>),
    HelpClose(Result<(), String>),
    CookbookOpen,
    CookbookClose,
    SummaryOpen(Result<(), String>),
    SummaryClose(Result<(), String>),
    ConsoleOpen,
    ConsoleClose,
    ArchiveOpen(Result<(), String>),
    ArchiveClose,
    SaveOnExit,
    Restore(bool, usize),
    ExpandArchiveEntry(bool, usize),
    DeleteArchiveEntry(bool, usize),
    ArchiveName(String, usize),
    ArchiveNameChange(bool, usize),
    Name(String),
    NameChange(bool),
}
