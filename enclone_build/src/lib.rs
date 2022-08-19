pub mod prepare_for_apocalypse;

const VERSION_STRING: &str = env!("VERSION_STRING");

// WARNING: the version string will not be up to date unless enclone_build/build.rs is touched
// before compiling.  Use current_version_string() to always get the current version.

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}
