// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#[derive(Debug, Clone)]
pub enum Message {
    InputChanged(String),
    SubmitButtonPressed(Result<(), String>),
    BackButtonPressed(Result<(), String>),
    ForwardButtonPressed(Result<(), String>),
    OpenModalSummary,
    OpenModalHelp,
    CloseModalHelp,
    OpenModalCookbook,
    CancelButtonPressed,
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
}
