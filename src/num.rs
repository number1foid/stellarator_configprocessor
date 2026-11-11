//! Constants, three-component vectors and conversions between counts and
//! real numbers.
//!
//! The conversions avoid `as` casts: a count becomes a real number through
//! `u32`, which is exact, and a real number becomes a count by bisection.

use crate::error::{ensure, Code, Error};
use crate::limits::BISECTION_STEPS;

/// Vacuum permeability [T m / A]. The classical value is used so that
/// results match VMEC and the `0.2` factor in the PROCESS peak-field formula.
pub const MU0: f64 = 4.0e-7 * core::f64::consts::PI;
/// Version of the PROCESS `stella_conf.json` schema this program writes.
pub const SCHEMA_VERSION: &str = "process-3.4";
/// Program name and version stored in the output files.
pub const GENERATOR: &str = "stellarator_configprocessor 0.2.0";

/// A point or vector in Cartesian coordinates.
pub type Vec3 = [f64; 3];

/// Sum of two vectors.
#[must_use]
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Difference of two vectors.
#[must_use]
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// A vector multiplied by a number.
#[must_use]
pub fn scale(a: Vec3, factor: f64) -> Vec3 {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

/// Scalar product.
#[must_use]
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Vector product.
#[must_use]
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Length of a vector.
#[must_use]
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// Whether all three components are finite.
#[must_use]
pub fn finite_vector(a: Vec3) -> bool {
    a[0].is_finite() && a[1].is_finite() && a[2].is_finite()
}

/// A count as a real number, exactly.
///
/// # Errors
/// Returns an error if the count does not fit in 32 bits.
pub fn real(count: usize) -> Result<f64, Error> {
    ensure(count <= 0x7FFF_FFFF, Code::Limit, "a count is too large to be used as a real number")?;
    let small = u32::try_from(count)?;
    Ok(f64::from(small))
}

/// The largest count not above `value`, but at most `limit`.
///
/// # Errors
/// Returns an error if `value` is negative or not finite, or `limit` is too
/// large to be compared with a real number.
pub fn floor_count(value: f64, limit: usize) -> Result<usize, Error> {
    ensure(
        value.is_finite() && value >= 0.0,
        Code::Numeric,
        "a count was computed from a negative or non-finite number",
    )?;
    ensure(real(limit)? >= 0.0, Code::Limit, "limit is not a count")?;
    let mut low = 0_usize;
    let mut high = limit;
    for _ in 0..BISECTION_STEPS {
        if low < high {
            let middle = low + (high - low).div_ceil(2);
            if real(middle)? <= value {
                low = middle;
            } else {
                high = middle - 1;
            }
        }
    }
    Ok(low)
}

/// The smallest count not below `value`.
///
/// # Errors
/// Returns an error if `value` is negative, not finite or above `limit`.
pub fn ceil_count(value: f64, limit: usize) -> Result<usize, Error> {
    ensure(value.is_finite(), Code::Numeric, "a count was computed from a non-finite number")?;
    let below = floor_count(value, limit)?;
    let count = if real(below)? < value { below + 1 } else { below };
    ensure(count <= limit, Code::Limit, "a computed count exceeds its limit")?;
    Ok(count)
}

/// Whether every value of a slice is finite.
#[must_use]
pub fn all_finite(values: &[f64]) -> bool {
    let mut finite = true;
    for value in values {
        finite = finite && value.is_finite();
    }
    finite
}
