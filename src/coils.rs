//! Coils sampled at the same number of points, uniform in arc length.

use crate::curve::{resample, SplineWork};
use crate::error::{ensure, Code, Error};
use crate::limits::{MAX_COILS, MAX_SAMPLES, MIN_SAMPLES};
use crate::makegrid::CoilSet;
use crate::num::{add, all_finite, real, scale, Vec3};

const BAD: Code = Code::Numeric;

/// All coils, each sampled at `n_points` points. Arrays over points are
/// stored coil after coil: sample `j` of coil `i` is at `i * n_points + j`.
#[derive(Debug, Clone)]
pub struct SampledCoils {
    /// Number of coils.
    pub n_coils: usize,
    /// Number of samples per coil.
    pub n_points: usize,
    /// Sample positions [m].
    pub points: Vec<Vec3>,
    /// Unit tangents at the samples.
    pub tangents: Vec<Vec3>,
    /// Curvature at the samples [1/m].
    pub curvature: Vec<f64>,
    /// Length of each coil [m].
    pub lengths: Vec<f64>,
    /// Current of each coil [A].
    pub currents: Vec<f64>,
    /// Mean position of the samples of each coil [m].
    pub centroids: Vec<Vec3>,
}

/// Buffers for `n_coils` coils sampled at `n_points` points each.
///
/// # Errors
/// Returns an error if a size is outside the limits of this program.
pub fn new_sampled_coils(n_coils: usize, n_points: usize) -> Result<SampledCoils, Error> {
    ensure(
        (1..=MAX_COILS).contains(&n_coils),
        Code::Limit,
        "the number of coils is outside the limits",
    )?;
    ensure(
        (MIN_SAMPLES..=MAX_SAMPLES).contains(&n_points),
        Code::Settings,
        "the number of samples per coil must be between 32 and 1024",
    )?;
    let total = n_coils * n_points;
    Ok(SampledCoils {
        n_coils,
        n_points,
        points: vec![[0.0; 3]; total],
        tangents: vec![[0.0; 3]; total],
        curvature: vec![0.0; total],
        lengths: vec![0.0; n_coils],
        currents: vec![0.0; n_coils],
        centroids: vec![[0.0; 3]; n_coils],
    })
}

/// The samples of one coil.
///
/// # Errors
/// Returns an error if `coil` is not a valid coil index.
pub fn coil_points(coils: &SampledCoils, coil: usize) -> Result<&[Vec3], Error> {
    ensure(coil < coils.n_coils, BAD, "a coil index is out of range")?;
    ensure(
        coils.points.len() == coils.n_coils * coils.n_points,
        BAD,
        "the coil buffers are damaged",
    )?;
    Ok(&coils.points[coil * coils.n_points..(coil + 1) * coils.n_points])
}

/// Arc length represented by each sample of one coil [m].
///
/// # Errors
/// Returns an error if `coil` is not a valid coil index or its length is
/// not positive.
pub fn arc_step(coils: &SampledCoils, coil: usize) -> Result<f64, Error> {
    ensure(coil < coils.n_coils, BAD, "a coil index is out of range")?;
    let step = coils.lengths[coil] / real(coils.n_points)?;
    ensure(step > 0.0 && step.is_finite(), BAD, "a coil has no length")?;
    Ok(step)
}

/// Sample every coil of `set` into `out`.
///
/// # Errors
/// Returns an error if `out` was created for a different number of coils,
/// a coil cannot be sampled, or a current is not finite.
pub fn sample_coils(
    set: &CoilSet,
    work: &mut SplineWork,
    out: &mut SampledCoils,
) -> Result<(), Error> {
    let n = out.n_points;
    ensure(set.currents.len() == out.n_coils, BAD, "the coil buffers do not match the coil set")?;
    ensure(set.offsets.len() == out.n_coils + 1, BAD, "the coil set is damaged")?;
    ensure(all_finite(&set.currents), BAD, "a coil current is not finite")?;
    for coil in 0..out.n_coils {
        let source = &set.points[set.offsets[coil]..set.offsets[coil + 1]];
        let span = coil * n..(coil + 1) * n;
        out.lengths[coil] = resample(
            source,
            work,
            &mut out.points[span.clone()],
            &mut out.tangents[span.clone()],
            &mut out.curvature[span.clone()],
        )?;
        let mut sum = [0.0; 3];
        for point in &out.points[span] {
            sum = add(sum, *point);
        }
        out.centroids[coil] = scale(sum, 1.0 / real(n)?);
        out.currents[coil] = set.currents[coil];
    }
    Ok(())
}
