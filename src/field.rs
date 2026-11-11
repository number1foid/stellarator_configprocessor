//! Biot-Savart field of closed polygons of straight current segments.
//!
//! A closed filament is represented as a polygon. The field of each straight
//! segment is evaluated with the exact closed-form expression, so the only
//! discretisation error is the polygon approximation of the curve.

use crate::coils::SampledCoils;
use crate::error::{ensure, Code, Error};
use crate::limits::MAX_INDUCTANCE_SAMPLES;
use crate::num::{add, cross, dot, finite_vector, norm, scale, sub, Vec3, MU0};

const BAD: Code = Code::Numeric;

/// Add to `total` the field `[T]` at `point` of a closed polygon through
/// `vertices` that carries `current` `[A]` in the order of the vertices.
///
/// A point lying exactly on a segment receives no contribution from that
/// segment (the field of a line current is undefined there).
///
/// # Errors
/// Returns an error if the polygon has fewer than two or too many vertices,
/// or the current, the point or the result is not finite.
pub fn add_polygon_field(
    point: Vec3,
    vertices: &[Vec3],
    current: f64,
    total: &mut Vec3,
) -> Result<(), Error> {
    ensure(vertices.len() >= 2, BAD, "a polygon needs at least two vertices")?;
    ensure(
        vertices.len() <= MAX_INDUCTANCE_SAMPLES,
        Code::Limit,
        "a polygon has too many vertices",
    )?;
    ensure(
        current.is_finite() && finite_vector(point),
        BAD,
        "a field point or current is not finite",
    )?;
    let mut sum = [0.0; 3];
    let mut from_start = sub(point, vertices[vertices.len() - 1]);
    let mut start_distance = norm(from_start);
    for vertex in vertices {
        let from_end = sub(point, *vertex);
        let end_distance = norm(from_end);
        let product = start_distance * end_distance;
        let denominator = product * (product + dot(from_start, from_end));
        // The denominator vanishes only for points on the segment itself.
        if denominator > 1.0e-14 * product * product + 1.0e-300 {
            let factor = (start_distance + end_distance) / denominator;
            sum = add(sum, scale(cross(from_start, from_end), factor));
        }
        from_start = from_end;
        start_distance = end_distance;
    }
    *total = add(*total, scale(sum, current * MU0 / (4.0 * core::f64::consts::PI)));
    ensure(finite_vector(*total), BAD, "a computed field is not finite")?;
    Ok(())
}

/// Field `[T]` at `point` of the centre filaments of the coils selected by
/// `include`, which holds one flag per coil.
///
/// # Errors
/// Returns an error if `include` does not hold one flag per coil, or a coil
/// field cannot be evaluated.
pub fn filament_field(point: Vec3, coils: &SampledCoils, include: &[bool]) -> Result<Vec3, Error> {
    ensure(coils.n_points >= 2, BAD, "a coil needs at least two samples")?;
    ensure(include.len() == coils.n_coils, BAD, "one flag per coil is needed")?;
    ensure(
        coils.points.len() == coils.n_coils * coils.n_points,
        BAD,
        "the coil buffers are damaged",
    )?;
    let mut total = [0.0; 3];
    for (coil, selected) in include.iter().enumerate() {
        if *selected {
            let vertices = &coils.points[coil * coils.n_points..(coil + 1) * coils.n_points];
            add_polygon_field(point, vertices, coils.currents[coil], &mut total)?;
        }
    }
    Ok(total)
}
