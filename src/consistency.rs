//! Check that a coil set and an equilibrium describe the same machine.
//!
//! The coil field is evaluated along the equilibrium's magnetic axis. For a
//! matching pair, the axis links all coil current (Ampere's law), the field
//! is tangent to the axis, and its strength equals the equilibrium's.

use crate::coils::SampledCoils;
use crate::equilibrium::{
    axis_field_strength, cartesian, interpolate, surface_point, Coefficients, Equilibrium,
};
use crate::error::{ensure, Code, Error};
use crate::field::filament_field;
use crate::limits::AXIS_POINTS;
use crate::num::{dot, norm, real, scale, sub, Vec3, MU0};
use core::f64::consts::TAU;

const BAD: Code = Code::Numeric;

/// Limits for the consistency tests. All are relative.
#[derive(Debug, Clone, Copy)]
pub struct Tolerances {
    /// Deviation from one of the linked fraction of the coil current.
    pub linkage: f64,
    /// Mean perpendicular field on the axis divided by the mean field.
    pub axis_misalignment: f64,
    /// Relative difference of the mean field on axis.
    pub field: f64,
}

/// The default tolerances.
pub const DEFAULT_TOLERANCES: Tolerances =
    Tolerances { linkage: 0.01, axis_misalignment: 0.02, field: 0.03 };

/// Coil field sampled along the equilibrium's magnetic axis; a reusable buffer.
#[derive(Debug, Clone)]
pub struct AxisField {
    /// Axis points `[m]`.
    pub points: Vec<Vec3>,
    /// Axis line elements `[m]`.
    pub line: Vec<Vec3>,
    /// Coil field at the points `[T]`.
    pub field: Vec<Vec3>,
    /// One flag per coil, all set: every coil contributes.
    pub all_coils: Vec<bool>,
}

/// Number of consistency tests.
pub const N_TESTS: usize = 5;
/// Test: the coil file declares the same number of field periods.
pub const TEST_PERIODS: usize = 0;
/// Test: the number of coils is a multiple of the number of field periods.
pub const TEST_DIVISIBILITY: usize = 1;
/// Test: the axis links all of the coil current.
pub const TEST_LINKAGE: usize = 2;
/// Test: the coil field is tangent to the axis.
pub const TEST_ALIGNMENT: usize = 3;
/// Test: the field strength on the axis agrees.
pub const TEST_FIELD: usize = 4;

/// Outcome of the consistency tests.
#[derive(Debug, Clone, Copy)]
pub struct Consistency {
    /// Field periods of the equilibrium.
    pub nfp_equilibrium: usize,
    /// Field periods declared in the coil file.
    pub nfp_coil_file: usize,
    /// Number of coils.
    pub n_coils: usize,
    /// Current linked by the axis divided by the total coil current.
    pub linkage_ratio: f64,
    /// Mean perpendicular field on the axis divided by the mean field.
    pub axis_misalignment: f64,
    /// Arc-length mean of the coil field strength on the axis `[T]`.
    pub b_axis_coils: f64,
    /// Arc-length mean of the equilibrium's field strength on the axis `[T]`.
    pub b_axis_equilibrium: f64,
    /// Arc-length mean of the toroidal component of the coil field `[T]`.
    pub b_axis_toroidal: f64,
    /// `b_axis_coils / b_axis_equilibrium`.
    pub field_ratio: f64,
    /// For each test, whether it failed; indexed by the `TEST_` constants.
    pub failed: [bool; N_TESTS],
    /// Whether every test passed.
    pub passed: bool,
}

/// A buffer for the axis field of `n_coils` coils.
///
/// # Errors
/// Returns an error if there are no coils.
pub fn new_axis_field(n_coils: usize) -> Result<AxisField, Error> {
    ensure(n_coils >= 1, BAD, "at least one coil is needed")?;
    ensure(AXIS_POINTS >= 16, BAD, "the axis needs at least 16 points")?;
    Ok(AxisField {
        points: vec![[0.0; 3]; AXIS_POINTS],
        line: vec![[0.0; 3]; AXIS_POINTS],
        field: vec![[0.0; 3]; AXIS_POINTS],
        all_coils: vec![true; n_coils],
    })
}

/// Evaluate the coil field along the full magnetic axis.
///
/// # Errors
/// Returns an error if the buffer does not match the coil set or the axis
/// or the field cannot be evaluated.
pub fn axis_field(
    eq: &Equilibrium,
    coils: &SampledCoils,
    surface: &mut Coefficients,
    axis: &mut AxisField,
) -> Result<(), Error> {
    ensure(eq.nfp >= 1, Code::Equilibrium, "the equilibrium has no field periods")?;
    ensure(axis.all_coils.len() == coils.n_coils, BAD, "the axis buffer does not match the coils")?;
    ensure(axis.points.len() == AXIS_POINTS, BAD, "the axis buffer has the wrong size")?;
    interpolate(eq, 0.0, surface)?;
    for (index, point) in axis.points.iter_mut().enumerate() {
        let phi = TAU * real(index)? / real(AXIS_POINTS)?;
        *point = cartesian(&surface_point(eq, surface, 0.0, phi)?, phi);
    }
    // Central differences on the closed axis give second-order line elements.
    for index in 0..AXIS_POINTS {
        let ahead = axis.points[(index + 1) % AXIS_POINTS];
        let behind = axis.points[(index + AXIS_POINTS - 1) % AXIS_POINTS];
        axis.line[index] = scale(sub(ahead, behind), 0.5);
        axis.field[index] = filament_field(axis.points[index], coils, &axis.all_coils)?;
    }
    Ok(())
}

/// Sums along the axis from which the test quantities follow.
#[derive(Debug, Clone, Copy, Default)]
struct AxisSums {
    arc: f64,
    strength: f64,
    toroidal: f64,
    circulation: f64,
    perpendicular: f64,
    equilibrium: f64,
}

fn axis_sums(eq: &Equilibrium, axis: &AxisField) -> Result<AxisSums, Error> {
    ensure(axis.line.len() == AXIS_POINTS, BAD, "the axis buffer has the wrong size")?;
    ensure(axis.field.len() == AXIS_POINTS, BAD, "the axis buffer has the wrong size")?;
    let mut sums = AxisSums::default();
    for index in 0..AXIS_POINTS {
        let phi = TAU * real(index)? / real(AXIS_POINTS)?;
        let (field, line) = (axis.field[index], axis.line[index]);
        let arc = norm(line);
        ensure(arc > 0.0, BAD, "the magnetic axis has no length")?;
        let parallel = dot(field, line) / arc;
        let across = sub(field, scale(line, parallel / arc));
        sums.arc += arc;
        sums.strength += norm(field) * arc;
        sums.toroidal += (field[1] * phi.cos() - field[0] * phi.sin()) * arc;
        sums.circulation += dot(field, line);
        sums.perpendicular += norm(across) * arc;
        sums.equilibrium += axis_field_strength(eq, phi)?.abs() * arc;
    }
    Ok(sums)
}

/// Run the consistency tests. A failed test is recorded, not returned as an
/// error; see [`enforce`].
///
/// `check_field = false` skips the field-strength comparison. It is used
/// when the coil currents are rescaled to match the equilibrium anyway.
///
/// # Errors
/// Returns an error if the axis field has not been evaluated or a tolerance
/// is not positive.
pub fn check_consistency(
    eq: &Equilibrium,
    coils: &SampledCoils,
    periods_in_file: usize,
    tolerances: &Tolerances,
    axis: &AxisField,
    check_field: bool,
) -> Result<Consistency, Error> {
    ensure(
        eq.nfp >= 1 && periods_in_file >= 1,
        BAD,
        "the number of field periods must be at least 1",
    )?;
    ensure(coils.currents.len() == coils.n_coils, BAD, "the coil buffers are damaged")?;
    ensure(
        tolerances.linkage > 0.0 && tolerances.axis_misalignment > 0.0 && tolerances.field > 0.0,
        Code::Settings,
        "tolerances must be positive",
    )?;
    let sums = axis_sums(eq, axis)?;
    let mut total_current = 0.0;
    for current in &coils.currents {
        total_current += current.abs();
    }
    ensure(total_current > 0.0, Code::CoilFile, "all coil currents are zero")?;
    let linkage_ratio = sums.circulation.abs() / (MU0 * total_current);
    let b_axis_coils = sums.strength / sums.arc;
    let b_axis_equilibrium = sums.equilibrium / sums.arc;
    let axis_misalignment =
        if b_axis_coils > 0.0 { sums.perpendicular / sums.strength } else { f64::NAN };
    let field_ratio =
        if b_axis_equilibrium > 0.0 { b_axis_coils / b_axis_equilibrium } else { f64::NAN };
    let mut failed = [false; N_TESTS];
    failed[TEST_PERIODS] = periods_in_file != eq.nfp;
    failed[TEST_DIVISIBILITY] = coils.n_coils % eq.nfp != 0;
    failed[TEST_LINKAGE] = (linkage_ratio - 1.0).abs() > tolerances.linkage;
    failed[TEST_ALIGNMENT] =
        axis_misalignment.is_nan() || axis_misalignment > tolerances.axis_misalignment;
    failed[TEST_FIELD] =
        check_field && (field_ratio.is_nan() || (field_ratio - 1.0).abs() > tolerances.field);
    Ok(Consistency {
        nfp_equilibrium: eq.nfp,
        nfp_coil_file: periods_in_file,
        n_coils: coils.n_coils,
        linkage_ratio,
        axis_misalignment,
        b_axis_coils,
        b_axis_equilibrium,
        b_axis_toroidal: sums.toroidal / sums.arc,
        field_ratio,
        failed,
        passed: failed == [false; N_TESTS],
    })
}

/// Stop if the check failed and no override is set.
///
/// # Errors
/// Returns an error of class [`Code::Inconsistent`] if a test failed and
/// `allow_inconsistent` is false, or if the result is not finite.
pub fn enforce(result: &Consistency, allow_inconsistent: bool) -> Result<(), Error> {
    ensure(result.b_axis_coils.is_finite(), BAD, "the coil field on the axis is not finite")?;
    ensure(
        result.passed || allow_inconsistent,
        Code::Inconsistent,
        "coils and equilibrium appear inconsistent",
    )
}
