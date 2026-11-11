//! Command-line entry point.

use std::process::ExitCode;
use stellarator_configprocessor::cli::{execute, prepare, report_failure, EXIT_INPUT, EXIT_OK};

/// Counts memory requests and refuses them once start-up is complete.
#[global_allocator]
static GUARD: alloc_guard::Guard = alloc_guard::Guard;

fn main() -> ExitCode {
    let mut out = std::io::stdout();
    let mut err = std::io::stderr();
    let status = match prepare(std::env::args().skip(1), &mut out) {
        Ok(Some(mut session)) => {
            // Start-up is complete: from here on no memory may be requested.
            alloc_guard::lock();
            execute(&mut session, &mut out, &mut err)
        }
        Ok(None) => EXIT_OK,
        Err(error) => report_failure(&mut err, &error, None).unwrap_or(EXIT_INPUT),
    };
    ExitCode::from(status)
}
