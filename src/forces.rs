//! Electromagnetic forces on the coils.
//!
//! The force per unit length of the coil centre line is the Lorentz force
//! density `j x B` integrated over the winding pack cross-section:
//!
//! ```text
//! F' = I * sum over the cross-section of  w * (dP/ds) x B
//! ```
//!
//! where `P` is a point of the cross-section, `dP/ds` its derivative along
//! the centre line and `w` the quadrature weight. `dP/ds` is the local
//! current direction times the ratio of the local to the centre-line length:
//! on the outside of a bend the conductor is longer than its centre line.
//! That factor matters for the coil's own field, which changes sign across
//! the cross-section, and is needed to obtain the correct hoop force.

use crate::coilgeom::SymmetryClasses;
use crate::coils::{arc_step, SampledCoils};
use crate::error::{ensure, Code, Error};
use crate::num::{add, cross, dot, norm, real, scale, sub, Vec3};
use crate::winding_pack::{Frames, PackLayout, PackWork};

const BAD: Code = Code::Numeric;

/// Forces on the coils, from the representative of each symmetry class.
#[derive(Debug, Clone, Copy)]
pub struct ForceResult {
    /// Largest force per unit length [N/m].
    pub max_line_force: f64,
    /// Largest cross-section-averaged force density [N/m^3].
    pub max_force_density: f64,
    /// Largest component along the in-plane normal, per area [N/m^3].
    pub max_radial_force_density: f64,
    /// Largest component along the toroidal normal, per area [N/m^3].
    pub max_lateral_force_density: f64,
    /// Largest net force towards the torus axis on one coil [N].
    pub centering_max: f64,
    /// Smallest net force towards the torus axis on one coil [N].
    pub centering_min: f64,
    /// Mean net force towards the torus axis per coil [N].
    pub centering_avg: f64,
}

/// Force per unit length at one sample, from the field at its cross-section.
fn line_force(
    layout: &PackLayout,
    work: &PackWork,
    base: usize,
    stations: (usize, usize),
    scale_factor: f64,
) -> Vec3 {
    let stride = layout.n_toroidal + 1;
    let corners = (layout.n_radial + 1) * stride;
    let (station, n) = stations;
    let here = base + station * corners;
    let ahead = base + ((station + 1) % n) * corners;
    let behind = base + ((station + n - 1) % n) * corners;
    let mut force = [0.0; 3];
    for a in 0..=layout.n_radial {
        for b in 0..=layout.n_toroidal {
            let corner = a * stride + b;
            let along =
                sub(work.samples[ahead + corner].point, work.samples[behind + corner].point);
            let lorentz = cross(along, work.samples[here + corner].field);
            force = add(force, scale(lorentz, layout.weight[a][b] * scale_factor));
        }
    }
    force
}

/// Forces from the winding pack field of the class representatives.
///
/// `work` must hold the field of every representative at every sample, as
/// computed for `layout`. The net force towards the torus axis on each
/// representative is written to `centering`.
///
/// # Errors
/// Returns an error if the buffers do not match or a force is not finite.
pub fn coil_forces(
    coils: &SampledCoils,
    frames: &Frames,
    classes: &SymmetryClasses,
    layout: &PackLayout,
    work: &PackWork,
    centering: &mut [f64],
) -> Result<ForceResult, Error> {
    let n = coils.n_points;
    let corners = (layout.n_radial + 1) * (layout.n_toroidal + 1);
    ensure(
        work.samples.len() >= classes.n_classes * n * corners,
        BAD,
        "forces need the field at every sample of each coil",
    )?;
    ensure(centering.len() >= classes.n_classes, Code::Limit, "the class buffer is too small")?;
    let (mut line, mut radial, mut lateral) = (0.0_f64, 0.0_f64, 0.0_f64);
    let (mut highest, mut lowest, mut weighted, mut members) =
        (f64::NEG_INFINITY, f64::INFINITY, 0.0, 0);
    for (class, slot) in centering.iter_mut().enumerate().take(classes.n_classes) {
        let coil = classes.representatives[class];
        let step = arc_step(coils, coil)?;
        let mut net = [0.0; 3];
        for station in 0..n {
            let index = coil * n + station;
            let factor = coils.currents[coil] / (2.0 * step);
            let force = line_force(layout, work, class * n * corners, (station, n), factor);
            line = line.max(norm(force));
            radial = radial.max(dot(force, frames.radial[index]).abs());
            lateral = lateral.max(dot(force, frames.toroidal[index]).abs());
            net = add(net, scale(force, step));
        }
        let centroid = coils.centroids[coil];
        let outward =
            (net[0] * centroid[0] + net[1] * centroid[1]) / centroid[0].hypot(centroid[1]);
        ensure(outward.is_finite(), BAD, "a coil force is not finite")?;
        *slot = -outward;
        highest = highest.max(-outward);
        lowest = lowest.min(-outward);
        weighted -= outward * real(classes.multiplicity[class])?;
        members += classes.multiplicity[class];
    }
    ensure(members >= 1 && line.is_finite(), BAD, "the coil forces are not finite")?;
    let area = layout.wp.area;
    Ok(ForceResult {
        max_line_force: line,
        max_force_density: line / area,
        max_radial_force_density: radial / area,
        max_lateral_force_density: lateral / area,
        centering_max: highest,
        centering_min: lowest,
        centering_avg: weighted / real(members)?,
    })
}
