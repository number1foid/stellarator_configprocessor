//! Shared helpers of the test suite.

use std::io::Write as _;
use std::path::PathBuf;
use stellarator_configprocessor::equilibrium::{read_wout, Equilibrium};
use stellarator_configprocessor::error::Error;
use stellarator_configprocessor::makegrid::{read_makegrid, CoilSet, GroupFactors};

/// No group factors.
pub(crate) const NO_FACTORS: GroupFactors =
    GroupFactors { count: 0, group: [0.0; 32], factor: [1.0; 32] };

/// Fail the test unless `value` is within `tolerance` (relative) of `expected`.
pub(crate) fn close(value: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = expected.abs().max(1.0e-300);
    assert!(
        (value - expected).abs() <= tolerance * scale,
        "{what}: got {value:e}, expected {expected:e} (relative tolerance {tolerance:e})"
    );
}

/// Fail the test unless `value` is within `tolerance` (absolute) of `expected`.
pub(crate) fn near(value: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (value - expected).abs() <= tolerance,
        "{what}: got {value:e}, expected {expected:e} (absolute tolerance {tolerance:e})"
    );
}

/// Path of a file in `data/local`, which holds the public reference inputs
/// fetched by `scripts/fetch_reference_data.py`.
pub(crate) fn data_file(name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data");
    path.push("local");
    path.push(name);
    path
}

/// The reference machine `wout`/`coils` pair, or `None` (with a note on the
/// test output) if the files have not been fetched.
pub(crate) fn reference_case(
    wout: &str,
    coils: &str,
) -> Result<Option<(Equilibrium, CoilSet)>, Error> {
    let (Ok(bytes), Ok(text)) =
        (std::fs::read(data_file(wout)), std::fs::read_to_string(data_file(coils)))
    else {
        writeln!(
            std::io::stderr(),
            "skipped: {wout} or {coils} missing; run scripts/fetch_reference_data.py"
        )?;
        return Ok(None);
    };
    Ok(Some((read_wout(&bytes)?, read_makegrid(&text, &NO_FACTORS)?)))
}
