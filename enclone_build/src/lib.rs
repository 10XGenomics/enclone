use chrono::{TimeZone, Utc};
use itertools::Itertools;

use std::env;
use string_utils::TextUtils;

use backtrace::Backtrace;

use std::{panic, sync::atomic::AtomicBool, sync::atomic::Ordering::SeqCst, thread};

static PANICKING: AtomicBool = AtomicBool::new(false);

/// Set up panic handling.
/// This function ensures that we always collect a stack trace, regardless of
/// whether or not an env var was set. It also limits concurrent panics to a
/// single thread, whichever panics first, to ensure that if more than one thread
/// panics at once, we don't end up with interleaved stack trace messages.
pub fn set_panic_handler(args: &[String]) {
    let now = Utc::now().naive_utc().timestamp();
    let build_date = version_string().after(":").between(": ", " :").to_string();
    let build_datetime = format!("{} 00:00:00", build_date);
    let then = Utc
        .datetime_from_str(&build_datetime, "%Y-%m-%d %H:%M:%S")
        .unwrap()
        .timestamp();
    let days_since_build = (now - then) / (60 * 60 * 24);
    let mut elapsed_message = String::new();
    if days_since_build > 30 {
        elapsed_message = format!(
            "Your build is {} days old.  You might want to check \
            to see if there is a newer build now.\n\n",
            days_since_build
        );
    }

    let trailer = format!(
        "You have probably encountered an internal \
            error in enclone.\n\n\
            Please email us at enclone@10xgenomics.com, including the traceback shown\n\
            above and also the following version information:\n\
            {} : {}.\n\n\
            Your command was:\n\n{}\n\n\
            {}\
            🌸 Thank you for reporting this bug and have a nice day! 🌸",
        env!("CARGO_PKG_VERSION"),
        version_string(),
        args.iter().format(" "),
        elapsed_message,
    );

    assert_eq!(thread::current().name().unwrap_or("<none>"), "main");
    let _ = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        // Prevent multiple threads from issuing tracebacks.
        if PANICKING.load(SeqCst) {
            return;
        }
        PANICKING.store(true, SeqCst);

        // Get backtrace.
        let backtrace = Backtrace::new();

        let thread = thread::current();
        let thread = thread.name().unwrap_or("unnamed");
        match info.location() {
            Some(location) => eprintln!(
                "thread '{}' panicked at {}:{}",
                thread,
                location.file(),
                location.line()
            ),
            None => eprintln!("thread '{thread}' panicked "),
        };

        let msg = match info.payload().downcast_ref::<&'static str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => s.as_str(),
                None => "Box<Any>",
            },
        };
        eprintln!("{msg}\n{backtrace:?}\n{trailer}");

        std::process::exit(101);
    }));
}

const VERSION_STRING: &str = env!("VERSION_STRING");

// WARNING: the version string will not be up to date unless enclone_build/build.rs is touched
// before compiling.  Use current_version_string() to always get the current version.

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}
