//! Error type and the assertion helper used throughout the crate.
//!
//! Every check in this crate has the form `ensure(condition, code, message)?`:
//! a side-effect-free Boolean test that, when it fails, returns an error to
//! the caller. Errors carry no heap data, so they can be created after the
//! allocator has been locked.

use core::fmt;

/// Class of a failure. It decides the exit code of the program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Code {
    /// The command line is wrong.
    Usage,
    /// The equilibrium file is missing, unreadable or unusable.
    Equilibrium,
    /// The coil file is missing, malformed or unusable.
    CoilFile,
    /// A run setting is missing or invalid.
    Settings,
    /// The computed configuration does not satisfy the schema.
    Schema,
    /// Coils and equilibrium do not match and no override was given.
    Inconsistent,
    /// Reading or writing a file failed.
    Io,
    /// A computation received or produced a value outside its valid range.
    Numeric,
    /// A size exceeds the fixed limits of the program.
    Limit,
}

/// A failed check: what kind, a fixed description, and where in the input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Error {
    /// Class of the failure.
    pub code: Code,
    /// What went wrong.
    pub message: &'static str,
    /// Line of the input file the failure refers to; zero if none.
    pub line: usize,
}

impl Error {
    /// An error without a position in an input file.
    #[must_use]
    pub const fn new(code: Code, message: &'static str) -> Self {
        Self { code, message, line: 0 }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line > 0 {
            write!(out, "line {}: {}", self.line, self.message)
        } else {
            write!(out, "{}", self.message)
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(_source: std::io::Error) -> Self {
        Self::new(Code::Io, "writing output failed")
    }
}

impl From<fmt::Error> for Error {
    fn from(_source: fmt::Error) -> Self {
        Self::new(Code::Io, "formatting output failed")
    }
}

impl From<core::num::TryFromIntError> for Error {
    fn from(_source: core::num::TryFromIntError) -> Self {
        Self::new(Code::Limit, "a count does not fit its number type")
    }
}

/// The assertion of this crate: `Ok` if `condition` holds, an error if not.
///
/// # Errors
/// Returns an [`Error`] with the given code and message when `condition` is
/// false.
pub const fn ensure(condition: bool, code: Code, message: &'static str) -> Result<(), Error> {
    if condition {
        Ok(())
    } else {
        Err(Error::new(code, message))
    }
}

/// Like [`ensure`], for a check that refers to a line of an input file.
///
/// # Errors
/// Returns an [`Error`] with the given code, message and line when
/// `condition` is false.
pub const fn ensure_at(
    condition: bool,
    code: Code,
    message: &'static str,
    line: usize,
) -> Result<(), Error> {
    if condition {
        Ok(())
    } else {
        Err(Error { code, message, line })
    }
}
