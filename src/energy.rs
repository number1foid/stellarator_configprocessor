//! Magnetic energy and inductance of the coil set.
//!
//! The energy is `W = 1/2 * sum of I_i * Phi_i`, with `Phi_i` the flux linked
//! by coil `i`. The flux from other coils follows from Neumann's formula for
//! the centre filaments. The flux of a coil through itself diverges for a
//! filament; it is regularised with the geometric mean distance (GMD) of the
//! rectangular winding pack from itself.

use crate::coilgeom::SymmetryClasses;
use crate::coils::{arc_step, coil_points, SampledCoils};
use crate::curve::{resample, SplineWork};
use crate::error::{ensure, Code, Error};
use crate::limits::MAX_INDUCTANCE_SAMPLES;
use crate::num::{ceil_count, dot, real, sub, Vec3, MU0};
use crate::winding_pack::WindingPack;
use core::f64::consts::PI;

const BAD: Code = Code::Numeric;
/// Samples per geometric mean distance along the filament.
const RESOLUTION: f64 = 4.0;

/// Buffers for one coil resampled for the self-inductance integral.
#[derive(Debug, Clone)]
pub struct InductanceWork {
    /// Sample positions `[m]`.
    points: Vec<Vec3>,
    /// Unit tangents.
    tangents: Vec<Vec3>,
    /// Curvature; computed by the sampler but not used here.
    curvature: Vec<f64>,
}

/// Energy of the coil set.
#[derive(Debug, Clone, Copy)]
pub struct EnergyResult {
    /// Total magnetic energy `[J]`.
    pub energy: f64,
    /// Current per coil used for the effective inductance `[A]`.
    pub reference_current: f64,
    /// `2 W / I^2` with `I` the reference current per coil `[H]`.
    pub effective_inductance: f64,
}

/// Buffers for the self-inductance integral.
///
/// # Errors
/// Returns an error if the fixed sample limit is too small.
pub fn new_inductance_work() -> Result<InductanceWork, Error> {
    ensure(MAX_INDUCTANCE_SAMPLES >= 64, BAD, "the inductance sample limit is too small")?;
    ensure(
        MAX_INDUCTANCE_SAMPLES <= 1 << 20,
        Code::Limit,
        "the inductance sample limit is too large",
    )?;
    Ok(InductanceWork {
        points: vec![[0.0; 3]; MAX_INDUCTANCE_SAMPLES],
        tangents: vec![[0.0; 3]; MAX_INDUCTANCE_SAMPLES],
        curvature: vec![0.0; MAX_INDUCTANCE_SAMPLES],
    })
}

/// Square of the geometric mean distance of a rectangle from itself [m^2].
///
/// Closed form of `exp(2 <ln d>)` for two points uniformly distributed in a
/// rectangle with sides `a` and `b`. For a square, `GMD = 0.44705 a`.
///
/// # Errors
/// Returns an error if the winding pack has no area or the result is not
/// positive.
pub fn gmd_squared(wp: &WindingPack) -> Result<f64, Error> {
    let (a, b) = (wp.radial, wp.toroidal);
    ensure(a > 0.0 && b > 0.0, BAD, "winding pack dimensions must be positive")?;
    let (a2, b2) = (a * a, b * b);
    let exponent = 4.0 * b / (3.0 * a) * (a / b).atan()
        + 4.0 * a / (3.0 * b) * (b / a).atan()
        + b2 / (6.0 * a2) * (b / a).ln()
        + a2 / (6.0 * b2) * (a / b).ln()
        - (a2 * a2 - 6.0 * a2 * b2 + b2 * b2) / (6.0 * a2 * b2) * (a / b + b / a).ln();
    let squared = a * b * (exponent - 25.0 / 6.0).exp();
    ensure(
        squared > 0.0 && squared.is_finite(),
        BAD,
        "the geometric mean distance is not positive",
    )?;
    Ok(squared)
}

/// Self-inductance of one closed filament with a rectangular winding pack `[H]`.
///
/// The filament is resampled so that the spacing is a quarter of the GMD,
/// which resolves the regularised kernel.
///
/// # Errors
/// Returns an error if the curve cannot be sampled, or the winding pack is
/// so small relative to the coil that more samples are needed than allowed.
pub fn self_inductance(
    points: &[Vec3],
    wp: &WindingPack,
    spline: &mut SplineWork,
    work: &mut InductanceWork,
) -> Result<f64, Error> {
    ensure(points.len() >= 3, BAD, "a coil needs at least three samples")?;
    let softening = gmd_squared(wp)?;
    let coarse = points.len().max(16);
    ensure(coarse <= MAX_INDUCTANCE_SAMPLES, Code::Limit, "a coil has too many samples")?;
    let length = resample(
        points,
        spline,
        &mut work.points[..coarse],
        &mut work.tangents[..coarse],
        &mut work.curvature[..coarse],
    )?;
    let wanted = RESOLUTION * length / softening.sqrt();
    ensure(
        wanted < real(MAX_INDUCTANCE_SAMPLES)?,
        Code::Limit,
        "the winding pack is too small relative to the coil length for the inductance integral",
    )?;
    let count = ceil_count(wanted, MAX_INDUCTANCE_SAMPLES)?.max(64);
    let length = resample(
        points,
        spline,
        &mut work.points[..count],
        &mut work.tangents[..count],
        &mut work.curvature[..count],
    )?;
    let mut total = 0.0;
    for first in 0..count {
        for second in 0..count {
            let gap = sub(work.points[first], work.points[second]);
            total += dot(work.tangents[first], work.tangents[second])
                / (dot(gap, gap) + softening).sqrt();
        }
    }
    let step = length / real(count)?;
    ensure(total.is_finite() && total > 0.0, BAD, "the self-inductance is not positive")?;
    Ok(MU0 / (4.0 * PI) * total * step * step)
}

/// Flux through coil `coil` produced by all other coils `[Wb]`.
///
/// # Errors
/// Returns an error if `coil` is not a valid index or the result is not
/// finite.
pub fn flux_from_others(coils: &SampledCoils, coil: usize) -> Result<f64, Error> {
    ensure(coils.tangents.len() == coils.points.len(), BAD, "the coil buffers are damaged")?;
    ensure(coil < coils.n_coils, BAD, "a coil index is out of range")?;
    let n = coils.n_points;
    let own = coil_points(coils, coil)?;
    let mut total = 0.0;
    for other in 0..coils.n_coils {
        if other != coil {
            // Tangent times line element times current of every source sample.
            let strength = arc_step(coils, other)? * coils.currents[other];
            let mut linked = 0.0;
            for (index, point) in own.iter().enumerate() {
                let tangent = coils.tangents[coil * n + index];
                for source in other * n..(other + 1) * n {
                    let gap = sub(*point, coils.points[source]);
                    linked += dot(tangent, coils.tangents[source]) / dot(gap, gap).sqrt();
                }
            }
            total += linked * strength;
        }
    }
    ensure(total.is_finite(), BAD, "the flux linked by a coil is not finite")?;
    Ok(MU0 / (4.0 * PI) * total * arc_step(coils, coil)?)
}

/// Total magnetic energy of the coil set.
///
/// Coils in the same symmetry class contribute equally, so only one
/// representative per class is evaluated. The self-inductance of each
/// representative is written to `self_inductances`.
///
/// # Errors
/// Returns an error if the classes do not cover the coil set, a buffer is
/// too small, or the energy is not positive.
pub fn magnetic_energy(
    coils: &SampledCoils,
    wp: &WindingPack,
    classes: &SymmetryClasses,
    spline: &mut SplineWork,
    work: &mut InductanceWork,
    self_inductances: &mut [f64],
) -> Result<EnergyResult, Error> {
    ensure(coils.currents.len() == coils.n_coils, BAD, "the coil buffers are damaged")?;
    ensure(classes.n_classes >= 1, BAD, "at least one symmetry class is needed")?;
    ensure(
        self_inductances.len() >= classes.n_classes,
        Code::Limit,
        "the class buffer is too small",
    )?;
    let (mut energy, mut members, mut current_sum) = (0.0, 0, 0.0);
    for (class, slot) in self_inductances.iter_mut().enumerate().take(classes.n_classes) {
        let coil = classes.representatives[class];
        let current = coils.currents[coil];
        let own = self_inductance(coil_points(coils, coil)?, wp, spline, work)?;
        let flux = own * current + flux_from_others(coils, coil)?;
        energy += 0.5 * real(classes.multiplicity[class])? * current * flux;
        members += classes.multiplicity[class];
        *slot = own;
    }
    for current in &coils.currents {
        current_sum += current.abs();
    }
    ensure(members == coils.n_coils, BAD, "symmetry classes do not cover the coil set")?;
    ensure(energy > 0.0 && energy.is_finite(), BAD, "the magnetic energy is not positive")?;
    let reference_current = current_sum / real(coils.n_coils)?;
    let squared = reference_current * reference_current;
    ensure(squared > 0.0, BAD, "the reference current is zero")?;
    Ok(EnergyResult { energy, reference_current, effective_inductance: 2.0 * energy / squared })
}
