//! Closed space curves: periodic cubic spline, resampling, curvature.
//!
//! The spline parameter is the cumulative chord length, so it runs from 0 to
//! the period. The curve is twice continuously differentiable.

use crate::error::{ensure, Code, Error};
use crate::limits::{BISECTION_STEPS, MAX_FILE_POINTS, MAX_INDUCTANCE_SAMPLES};
use crate::num::{add, cross, finite_vector, norm, real, scale, sub, Vec3};

const BAD: Code = Code::Numeric;

/// Reusable buffers for fitting and sampling a spline.
#[derive(Debug, Clone)]
pub struct SplineWork {
    /// Chord length from each knot to the next, cyclic.
    chord: Vec<f64>,
    /// Parameter of each knot, with one extra entry for the period.
    knots: Vec<f64>,
    /// Second derivative at each knot.
    second: Vec<Vec3>,
    /// Scratch: modified upper diagonal of the linear system.
    upper: Vec<f64>,
    /// Scratch: solution of the system for the corner correction.
    corner: Vec<f64>,
    /// Cumulative arc length on the dense parameter grid.
    arc: Vec<f64>,
}

/// Position and derivatives of the spline at one parameter value.
#[derive(Debug, Clone, Copy)]
struct Sample {
    position: Vec3,
    first: Vec3,
    second: Vec3,
}

/// Number of dense intervals used to measure arc length.
fn dense_count(n_points: usize, n_samples: usize) -> usize {
    (16 * n_points).max(8 * n_samples)
}

/// Buffers for curves of up to `max_points` points sampled at up to
/// `max_samples` points.
///
/// # Errors
/// Returns an error if a size is zero or above the limits.
pub fn new_spline_work(max_points: usize, max_samples: usize) -> Result<SplineWork, Error> {
    ensure(max_points >= 3, BAD, "a closed curve needs at least 3 points")?;
    ensure(
        max_points <= MAX_FILE_POINTS && max_samples <= MAX_INDUCTANCE_SAMPLES,
        Code::Limit,
        "a curve has more points than this program supports",
    )?;
    Ok(SplineWork {
        chord: vec![0.0; max_points],
        knots: vec![0.0; max_points + 1],
        second: vec![[0.0; 3]; max_points],
        upper: vec![0.0; max_points],
        corner: vec![0.0; max_points],
        arc: vec![0.0; dense_count(max_points, max_samples) + 1],
    })
}

/// Largest index `i` below `values.len() - 1` with `values[i] <= target`.
fn interval_index(values: &[f64], target: f64) -> usize {
    let mut low = 0_usize;
    let mut high = values.len().saturating_sub(2);
    for _ in 0..BISECTION_STEPS {
        if low < high {
            let middle = low + (high - low).div_ceil(2);
            if values[middle] <= target {
                low = middle;
            } else {
                high = middle - 1;
            }
        }
    }
    low
}

/// Chord lengths, knots and the right-hand side of the spline equations.
fn prepare_system(points: &[Vec3], work: &mut SplineWork) -> Result<(), Error> {
    let n = points.len();
    ensure(n >= 3, BAD, "a closed curve needs at least 3 distinct points")?;
    ensure(work.chord.len() >= n, Code::Limit, "the spline buffers are too small")?;
    work.knots[0] = 0.0;
    for index in 0..n {
        let length = norm(sub(points[(index + 1) % n], points[index]));
        ensure(
            length > 0.0 && length.is_finite(),
            BAD,
            "a curve has repeated or non-finite points",
        )?;
        work.chord[index] = length;
        work.knots[index + 1] = work.knots[index] + length;
    }
    for index in 0..n {
        let before = (index + n - 1) % n;
        let ahead = scale(sub(points[(index + 1) % n], points[index]), 1.0 / work.chord[index]);
        let behind = scale(sub(points[index], points[before]), 1.0 / work.chord[before]);
        work.second[index] = scale(sub(ahead, behind), 6.0);
    }
    Ok(())
}

/// Solve the cyclic tridiagonal system for the second derivatives.
///
/// The cyclic system is reduced to a tridiagonal one plus a rank-one
/// correction (Sherman-Morrison). Both right-hand sides are eliminated in
/// one sweep.
fn solve_system(n: usize, work: &mut SplineWork) -> Result<(), Error> {
    ensure(n >= 3 && work.chord.len() >= n, BAD, "the spline system is too small")?;
    let last = work.chord[n - 1];
    let gamma = -2.0 * (last + work.chord[0]);
    let mut pivot = 2.0 * (last + work.chord[0]) - gamma;
    work.upper[0] = work.chord[0] / pivot;
    work.second[0] = scale(work.second[0], 1.0 / pivot);
    work.corner[0] = gamma / pivot;
    for index in 1..n {
        let lower = work.chord[index - 1];
        let mut diagonal = 2.0 * (lower + work.chord[index]);
        let mut rhs_corner = 0.0;
        if index == n - 1 {
            diagonal -= last * last / gamma;
            rhs_corner = last;
        }
        pivot = diagonal - lower * work.upper[index - 1];
        ensure(pivot.abs() > 0.0, BAD, "the spline system is singular")?;
        work.upper[index] = work.chord[index] / pivot;
        let reduced = sub(work.second[index], scale(work.second[index - 1], lower));
        work.second[index] = scale(reduced, 1.0 / pivot);
        work.corner[index] = (rhs_corner - lower * work.corner[index - 1]) / pivot;
    }
    for index in (0..n - 1).rev() {
        let upper = work.upper[index];
        work.second[index] = sub(work.second[index], scale(work.second[index + 1], upper));
        work.corner[index] -= upper * work.corner[index + 1];
    }
    let divisor = 1.0 + work.corner[0] + last * work.corner[n - 1] / gamma;
    ensure(divisor.abs() > 0.0, BAD, "the spline system is singular")?;
    let numerator = add(work.second[0], scale(work.second[n - 1], last / gamma));
    for index in 0..n {
        let correction = scale(numerator, work.corner[index] / divisor);
        work.second[index] = sub(work.second[index], correction);
    }
    Ok(())
}

/// Evaluate the fitted spline at parameter `t` in `[0, period]`.
fn evaluate(points: &[Vec3], work: &SplineWork, t: f64) -> Sample {
    let index = interval_index(&work.knots[..=points.len()], t);
    let next = (index + 1) % points.len();
    let width = work.chord[index];
    let near = (work.knots[index + 1] - t) / width;
    let far = 1.0 - near;
    let (y0, y1) = (points[index], points[next]);
    let (m0, m1) = (work.second[index], work.second[next]);
    let bend = add(scale(m0, near * near * near - near), scale(m1, far * far * far - far));
    let position = add(add(scale(y0, near), scale(y1, far)), scale(bend, width * width / 6.0));
    let slope = scale(sub(y1, y0), 1.0 / width);
    let twist = sub(
        scale(m1, (3.0 * far * far - 1.0) / 6.0 * width),
        scale(m0, (3.0 * near * near - 1.0) / 6.0 * width),
    );
    Sample { position, first: add(slope, twist), second: add(scale(m0, near), scale(m1, far)) }
}

/// Cumulative arc length on a dense, uniform parameter grid. Returns the
/// number of intervals.
fn measure_arc(points: &[Vec3], work: &mut SplineWork, n_samples: usize) -> Result<usize, Error> {
    ensure(n_samples >= 8, BAD, "need at least 8 samples per curve")?;
    let n = points.len();
    let dense = dense_count(n, n_samples);
    ensure(work.arc.len() > dense, Code::Limit, "the spline buffers are too small")?;
    let step = work.knots[n] / real(dense)?;
    ensure(step > 0.0 && step.is_finite(), BAD, "a curve has no length")?;
    work.arc[0] = 0.0;
    for index in 0..dense {
        let middle = (real(index)? + 0.5) * step;
        let speed = norm(evaluate(points, work, middle).first);
        work.arc[index + 1] = work.arc[index] + speed * step;
    }
    Ok(dense)
}

/// Fit a periodic spline through `points` and sample it at points uniformly
/// spaced in arc length. Returns the length of the curve `[m]`.
///
/// The three output slices must have the same length, which is the number of
/// samples: positions, unit tangents and curvature [1/m].
///
/// # Errors
/// Returns an error if the curve has fewer than three points or repeated
/// points, fewer than eight samples are requested, the buffers are too
/// small, or the result is not finite.
pub fn resample(
    points: &[Vec3],
    work: &mut SplineWork,
    positions: &mut [Vec3],
    tangents: &mut [Vec3],
    curvature: &mut [f64],
) -> Result<f64, Error> {
    let limit = MAX_FILE_POINTS.max(MAX_INDUCTANCE_SAMPLES);
    ensure(points.len() <= limit, Code::Limit, "a curve has too many points")?;
    let count = positions.len();
    ensure(count >= 8, BAD, "need at least 8 samples per curve")?;
    ensure(
        tangents.len() == count && curvature.len() == count,
        BAD,
        "sample buffers differ in size",
    )?;
    prepare_system(points, work)?;
    solve_system(points.len(), work)?;
    let dense = measure_arc(points, work, count)?;
    let length = work.arc[dense];
    let step = work.knots[points.len()] / real(dense)?;
    for index in 0..count {
        let target = real(index)? * length / real(count)?;
        let cell = interval_index(&work.arc[..=dense], target);
        let span = work.arc[cell + 1] - work.arc[cell];
        let t = (real(cell)? + (target - work.arc[cell]) / span) * step;
        let sample = evaluate(points, work, t);
        let speed = norm(sample.first);
        ensure(
            speed > 0.0 && finite_vector(sample.position),
            BAD,
            "a sampled curve point is not finite",
        )?;
        positions[index] = sample.position;
        tangents[index] = scale(sample.first, 1.0 / speed);
        curvature[index] = norm(cross(sample.first, sample.second)) / (speed * speed * speed);
    }
    ensure(length.is_finite() && length > 0.0, BAD, "a curve has no length")?;
    Ok(length)
}
