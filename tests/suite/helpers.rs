//! Shared helpers of the test suite.

use crate::synthetic::{circular_coils_text, torus_wout};
use std::io::Write as _;
use std::path::PathBuf;
use stellarator_configprocessor::equilibrium::{read_wout, Equilibrium};
use stellarator_configprocessor::error::{Code, Error};
use stellarator_configprocessor::makegrid::{read_makegrid, CoilSet, GroupFactors};
use stellarator_configprocessor::num::MU0;
use stellarator_configprocessor::pipeline::{new_workspace, run, Outcome, Settings, Workspace};
use stellarator_configprocessor::schema::KEYS;

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

/// The error of a result that must have failed.
pub(crate) fn failure<T>(result: Result<T, Error>, what: &str) -> Error {
    let error = result.err();
    assert!(error.is_some(), "{what}: an error was expected but the call succeeded");
    error.unwrap_or(Error::new(Code::Numeric, "no error"))
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

/// A circular torus (`R = 10 m`, `a = 2 m`, five field periods) inside 20
/// circular coils of radius 4 m carrying 1 MA each. The equilibrium's field
/// on axis is that of the ideal toroidal solenoid.
pub(crate) fn torus_case() -> Result<(Equilibrium, CoilSet), Error> {
    let field = ideal_toroid_field(20.0, 1.0e6, 10.0);
    let eq = read_wout(&torus_wout(10.0, 2.0, 5, field, 11)?)?;
    let set = read_makegrid(&circular_coils_text(20, 10.0, 4.0, 1.0e6, 120, 5)?, &NO_FACTORS)?;
    Ok((eq, set))
}

/// Field of an ideal toroidal solenoid at major radius `r` [T].
pub(crate) fn ideal_toroid_field(n_coils: f64, current: f64, r: f64) -> f64 {
    MU0 * n_coils * current / (core::f64::consts::TAU * r)
}

/// Run the whole calculation.
pub(crate) fn compute(
    eq: &mut Equilibrium,
    set: &mut CoilSet,
    settings: &Settings,
) -> Result<(Outcome, Workspace), Error> {
    let mut work = new_workspace(eq, set, settings)?;
    let outcome = run(eq, set, settings, &mut work, &mut std::io::sink())?;
    Ok((outcome, work))
}

/// The value stored under `name` in the configuration of `outcome`.
pub(crate) fn value(outcome: &Outcome, name: &str) -> f64 {
    let index = KEYS.iter().position(|key| key.name == name);
    assert!(index.is_some(), "the schema has no key `{name}`");
    outcome.config.values[index.unwrap_or(0)]
}
