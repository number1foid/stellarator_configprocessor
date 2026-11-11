//! Magnetic field inside coils with a finite rectangular winding pack.
//!
//! The field of an ideal filament diverges at the filament, so the peak
//! field on the conductor and the forces need a finite cross-section. The
//! model follows Lion et al. 2021, section 3.7: a rectangular winding pack
//! (WP) carrying a uniform current density, with one side normal along the
//! toroidal direction.
//!
//! Numerical scheme:
//!
//! * The WP of a coil is replaced by `n_r x n_t` filaments at the centres of
//!   equal cells, each carrying an equal share of the current.
//! * The field is evaluated on the `(n_r + 1) x (n_t + 1)` cell corners,
//!   which are the points farthest from the filaments and include the WP
//!   surface.
//! * Local correction. Close to the conductor the difference between a
//!   uniform current density and discrete filaments is dominated by the
//!   two-dimensional field of a straight bar. That difference is known in
//!   closed form and is added, so that only the smooth, curvature-related
//!   part of the field relies on the filament discretisation.
//! * Coils far from the evaluated coil are represented by their centre
//!   filament.

use crate::coilgeom::SymmetryClasses;
use crate::coils::SampledCoils;
use crate::error::{ensure, Code, Error};
use crate::field::{add_polygon_field, filament_field};
use crate::limits::{MAX_CORNERS, MAX_SCAN, MAX_SUB};
use crate::num::{add, cross, dot, is_finite, norm, real, scale, sub, Vec3, MU0};
use core::f64::consts::TAU;

const BAD: Code = Code::Numeric;
/// Largest number of filament points held in memory.
const MAX_BUNDLE_POINTS: usize = 16_000_000;
/// Coils closer than this many WP diagonals get the full filament bundle.
const NEAR_FACTOR: f64 = 3.0;

/// Rectangular winding pack cross-section.
#[derive(Debug, Clone, Copy)]
pub struct WindingPack {
    /// Width in the coil's in-plane normal direction [m].
    pub radial: f64,
    /// Width in the toroidal direction [m].
    pub toroidal: f64,
    /// Cross-section area [m^2].
    pub area: f64,
    /// Length of the diagonal [m].
    pub diagonal: f64,
}

/// Unit normals of the winding pack at every sample point; a reusable buffer.
#[derive(Debug, Clone)]
pub struct Frames {
    /// In-plane normal, pointing away from the coil centroid.
    pub radial: Vec<Vec3>,
    /// The toroidal direction made perpendicular to the filament tangent.
    pub toroidal: Vec<Vec3>,
}

/// Filament and corner positions of one winding pack size.
#[derive(Debug, Clone, Copy)]
pub struct PackLayout {
    /// The winding pack.
    pub wp: WindingPack,
    /// Filaments across the radial width.
    pub n_radial: usize,
    /// Filaments across the toroidal width.
    pub n_toroidal: usize,
    /// Radial offsets of the filaments [m].
    pub centre_u: [f64; MAX_SUB],
    /// Toroidal offsets of the filaments [m].
    pub centre_v: [f64; MAX_SUB],
    /// Radial offsets of the field points [m].
    pub corner_u: [f64; MAX_CORNERS],
    /// Toroidal offsets of the field points [m].
    pub corner_v: [f64; MAX_CORNERS],
    /// Uniform bar minus filament bundle, radial kernel, per field point.
    pub delta_u: [[f64; MAX_CORNERS]; MAX_CORNERS],
    /// Uniform bar minus filament bundle, toroidal kernel, per field point.
    pub delta_v: [[f64; MAX_CORNERS]; MAX_CORNERS],
    /// Quadrature weight of each field point; the weights sum to one.
    pub weight: [[f64; MAX_CORNERS]; MAX_CORNERS],
}

/// A field point in a winding pack and the field there.
#[derive(Debug, Clone, Copy, Default)]
pub struct FieldSample {
    /// Position [m].
    pub point: Vec3,
    /// Field [T].
    pub field: Vec3,
}

/// Buffers for the field in the winding pack of the class representatives.
#[derive(Debug, Clone)]
pub struct PackWork {
    /// Filament points of every coil for the current winding pack size.
    pub bundle: Vec<Vec3>,
    /// Centre-filament distance from each representative to every coil.
    pub distances: Vec<f64>,
    /// Per coil: modelled with the full bundle (the target and its neighbours).
    pub close: Vec<bool>,
    /// Per coil: modelled by its centre filament.
    pub far: Vec<bool>,
    /// Field points of the reference winding pack and the field there, per
    /// class and station.
    pub samples: Vec<FieldSample>,
    /// Largest field strength in the cross-section, per class and station [T].
    pub station_peak: Vec<f64>,
}

/// Fit of the peak field against winding pack size, in PROCESS's form.
#[derive(Debug, Clone, Copy)]
pub struct PeakFieldFit {
    /// Fit parameter.
    pub a1: f64,
    /// Fit parameter of the winding pack size term.
    pub a2: f64,
    /// Number of winding pack sizes used.
    pub count: usize,
    /// Radial width of each size [m].
    pub widths: [f64; MAX_SCAN],
    /// Computed peak field of each size [T].
    pub peaks: [f64; MAX_SCAN],
    /// Value of the fit at each size [T].
    pub fitted: [f64; MAX_SCAN],
    /// Largest relative deviation of the fit from the computed peak fields.
    pub max_relative_error: f64,
}

/// A winding pack with the given widths.
///
/// # Errors
/// Returns an error if a width is not positive and finite.
pub fn winding_pack(radial: f64, toroidal: f64) -> Result<WindingPack, Error> {
    ensure(
        radial.is_finite() && radial > 0.0,
        Code::Settings,
        "wp_radial must be a positive number",
    )?;
    ensure(
        toroidal.is_finite() && toroidal > 0.0,
        Code::Settings,
        "wp_toroidal must be a positive number",
    )?;
    Ok(WindingPack { radial, toroidal, area: radial * toroidal, diagonal: radial.hypot(toroidal) })
}

/// Buffers for the frames of `n_coils` coils with `n_points` samples each.
///
/// # Errors
/// Returns an error if a size is zero.
pub fn new_frames(n_coils: usize, n_points: usize) -> Result<Frames, Error> {
    ensure(n_coils >= 1, BAD, "at least one coil is needed")?;
    ensure(n_points >= 3, BAD, "at least three samples per coil are needed")?;
    Ok(Frames {
        radial: vec![[0.0; 3]; n_coils * n_points],
        toroidal: vec![[0.0; 3]; n_coils * n_points],
    })
}

/// Compute the unit normals of the winding pack at every sample point.
///
/// # Errors
/// Returns an error if the buffer does not match the coils, or a coil runs
/// nearly parallel to the toroidal direction.
pub fn coil_frames(coils: &SampledCoils, frames: &mut Frames) -> Result<(), Error> {
    let n = coils.n_points;
    ensure(
        frames.radial.len() == coils.points.len(),
        BAD,
        "the frame buffer does not match the coils",
    )?;
    ensure(
        frames.toroidal.len() == coils.points.len(),
        BAD,
        "the frame buffer does not match the coils",
    )?;
    for coil in 0..coils.n_coils {
        let mut outward = 0.0;
        for index in coil * n..(coil + 1) * n {
            let (point, tangent) = (coils.points[index], coils.tangents[index]);
            let phi = point[1].atan2(point[0]);
            let e_phi = [-phi.sin(), phi.cos(), 0.0];
            let across = sub(e_phi, scale(tangent, dot(e_phi, tangent)));
            let length = norm(across);
            ensure(
                length >= 0.2,
                Code::CoilFile,
                "a coil runs nearly parallel to the toroidal direction; only modular coils are supported",
            )?;
            frames.toroidal[index] = scale(across, 1.0 / length);
            frames.radial[index] = cross(tangent, frames.toroidal[index]);
            outward += dot(frames.radial[index], sub(point, coils.centroids[coil]));
        }
        if outward < 0.0 {
            for normal in &mut frames.radial[coil * n..(coil + 1) * n] {
                *normal = scale(*normal, -1.0);
            }
        }
    }
    Ok(())
}

/// Antiderivative of `x / (x^2 + y^2)` in `x` and `y`, without the `-y`
/// term, which cancels in the double difference.
fn bar_primitive(x: f64, y: f64) -> f64 {
    let squared = x * x + y * y;
    let log_term = if squared > 0.0 { 0.5 * y * squared.ln() } else { 0.0 };
    let atan_term = if x == 0.0 { 0.0 } else { x * (y / x).atan() };
    log_term + atan_term
}

/// Field integrals of a uniform straight rectangular bar.
///
/// For a bar occupying `|u'| <= width_u / 2`, `|v'| <= width_v / 2`, returns
/// the integrals over the cross-section of `(u - u') / rho^2` and
/// `(v - v') / rho^2` with `rho^2 = (u - u')^2 + (v - v')^2`. The field is
/// `mu0 J / (2 pi) * (t x n_u * F_u + t x n_v * F_v)`.
///
/// # Errors
/// Returns an error if a width is not positive or a position is not finite.
pub fn bar_kernel(u: f64, v: f64, width_u: f64, width_v: f64) -> Result<(f64, f64), Error> {
    ensure(width_u > 0.0 && width_v > 0.0, BAD, "bar widths must be positive")?;
    ensure(u.is_finite() && v.is_finite(), BAD, "a bar position is not finite")?;
    let (x1, x2) = (u - 0.5 * width_u, u + 0.5 * width_u);
    let (y1, y2) = (v - 0.5 * width_v, v + 0.5 * width_v);
    let f_u = bar_primitive(x2, y2) - bar_primitive(x1, y2) - bar_primitive(x2, y1)
        + bar_primitive(x1, y1);
    let f_v = bar_primitive(y2, x2) - bar_primitive(y1, x2) - bar_primitive(y2, x1)
        + bar_primitive(y1, x1);
    Ok((f_u, f_v))
}

/// Filament and field-point positions, quadrature weights and the local
/// correction for a winding pack with `n_radial x n_toroidal` filaments.
///
/// # Errors
/// Returns an error if a filament count is zero or above the limit.
pub fn pack_layout(
    wp: &WindingPack,
    n_radial: usize,
    n_toroidal: usize,
) -> Result<PackLayout, Error> {
    ensure(
        n_radial >= 1 && n_toroidal >= 1,
        Code::Settings,
        "wp_filaments entries must be at least 1",
    )?;
    ensure(
        n_radial <= MAX_SUB && n_toroidal <= MAX_SUB,
        Code::Settings,
        "wp_filaments entries must be at most 12",
    )?;
    let zero = [[0.0; MAX_CORNERS]; MAX_CORNERS];
    let mut layout = PackLayout {
        wp: *wp,
        n_radial,
        n_toroidal,
        centre_u: [0.0; MAX_SUB],
        centre_v: [0.0; MAX_SUB],
        corner_u: [0.0; MAX_CORNERS],
        corner_v: [0.0; MAX_CORNERS],
        delta_u: zero,
        delta_v: zero,
        weight: zero,
    };
    let (nr, nt) = (real(n_radial)?, real(n_toroidal)?);
    for index in 0..=n_radial {
        layout.corner_u[index] = wp.radial * (real(index)? / nr - 0.5);
    }
    for index in 0..n_radial {
        layout.centre_u[index] = wp.radial * ((real(index)? + 0.5) / nr - 0.5);
    }
    for index in 0..=n_toroidal {
        layout.corner_v[index] = wp.toroidal * (real(index)? / nt - 0.5);
    }
    for index in 0..n_toroidal {
        layout.centre_v[index] = wp.toroidal * ((real(index)? + 0.5) / nt - 0.5);
    }
    fill_correction(&mut layout)?;
    Ok(layout)
}

/// Quadrature weights and the closed-form difference between a uniform bar
/// and the filament bundle at every field point.
fn fill_correction(layout: &mut PackLayout) -> Result<(), Error> {
    let (n_radial, n_toroidal) = (layout.n_radial, layout.n_toroidal);
    ensure(n_radial <= MAX_SUB && n_toroidal <= MAX_SUB, BAD, "too many filaments")?;
    let cell = layout.wp.area / (real(n_radial)? * real(n_toroidal)?);
    ensure(cell > 0.0, BAD, "the winding pack has no area")?;
    for a in 0..=n_radial {
        let edge_a = if a == 0 || a == n_radial { 0.5 } else { 1.0 };
        for b in 0..=n_toroidal {
            let edge_b = if b == 0 || b == n_toroidal { 0.5 } else { 1.0 };
            let (u, v) = (layout.corner_u[a], layout.corner_v[b]);
            let (mut f_u, mut f_v) = bar_kernel(u, v, layout.wp.radial, layout.wp.toroidal)?;
            for centre_u in &layout.centre_u[..n_radial] {
                for centre_v in &layout.centre_v[..n_toroidal] {
                    let (du, dv) = (u - centre_u, v - centre_v);
                    f_u -= cell * du / (du * du + dv * dv);
                    f_v -= cell * dv / (du * du + dv * dv);
                }
            }
            layout.delta_u[a][b] = f_u;
            layout.delta_v[a][b] = f_v;
            layout.weight[a][b] = edge_a * edge_b / (real(n_radial)? * real(n_toroidal)?);
        }
    }
    Ok(())
}

/// Buffers for coils with `n_points` samples and up to `n_sub` filaments
/// across each side of the winding pack.
///
/// # Errors
/// Returns an error if the sizes are zero or need more memory than allowed.
pub fn new_pack_work(
    n_coils: usize,
    n_points: usize,
    n_sub: (usize, usize),
) -> Result<PackWork, Error> {
    ensure(n_coils >= 1 && n_points >= 3, BAD, "at least one coil with three samples is needed")?;
    ensure(
        (1..=MAX_SUB).contains(&n_sub.0) && (1..=MAX_SUB).contains(&n_sub.1),
        Code::Settings,
        "wp_filaments entries must be between 1 and 12",
    )?;
    let bundle = n_coils * n_points * n_sub.0 * n_sub.1;
    let corners = n_coils * n_points * (n_sub.0 + 1) * (n_sub.1 + 1);
    ensure(
        bundle <= MAX_BUNDLE_POINTS && corners <= MAX_BUNDLE_POINTS,
        Code::Limit,
        "too many coils, samples and winding pack filaments for the available memory",
    )?;
    Ok(PackWork {
        bundle: vec![[0.0; 3]; bundle],
        distances: vec![0.0; n_coils * n_coils],
        close: vec![false; n_coils],
        far: vec![false; n_coils],
        samples: vec![FieldSample::default(); corners],
        station_peak: vec![0.0; n_coils * n_points],
    })
}

/// Compute the filament points of every coil for the given layout.
///
/// # Errors
/// Returns an error if a buffer does not match the coils and the layout.
pub fn build_bundles(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    bundle: &mut [Vec3],
) -> Result<(), Error> {
    let n = coils.n_points;
    let filaments = layout.n_radial * layout.n_toroidal;
    ensure(
        bundle.len() >= coils.n_coils * filaments * n,
        Code::Limit,
        "the filament buffer is too small",
    )?;
    ensure(
        frames.radial.len() == coils.points.len(),
        BAD,
        "the frame buffer does not match the coils",
    )?;
    for coil in 0..coils.n_coils {
        for a in 0..layout.n_radial {
            for b in 0..layout.n_toroidal {
                let start = (coil * filaments + a * layout.n_toroidal + b) * n;
                for station in 0..n {
                    let index = coil * n + station;
                    let shift = add(
                        scale(frames.radial[index], layout.centre_u[a]),
                        scale(frames.toroidal[index], layout.centre_v[b]),
                    );
                    bundle[start + station] = add(coils.points[index], shift);
                }
            }
        }
    }
    Ok(())
}

/// Mark the coils modelled with the full bundle (`close`) and by one
/// filament (`far`) for target coil `coil`, whose distances to all coils
/// are given.
///
/// # Errors
/// Returns an error if the slices do not hold one entry per coil.
pub fn select_close(
    distances: &[f64],
    coil: usize,
    reach: f64,
    close: &mut [bool],
    far: &mut [bool],
) -> Result<(), Error> {
    ensure(coil < distances.len(), BAD, "a coil index is out of range")?;
    ensure(
        close.len() == distances.len() && far.len() == distances.len(),
        BAD,
        "one flag per coil is needed",
    )?;
    ensure(reach > 0.0, BAD, "the neighbour distance must be positive")?;
    for (index, distance) in distances.iter().enumerate() {
        close[index] = index == coil || *distance < reach;
        far[index] = !close[index];
    }
    Ok(())
}

/// Position of field point `(a, b)` of the cross-section at sample `index`.
fn corner_point(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    index: usize,
    a: usize,
    b: usize,
) -> Vec3 {
    let shift = add(
        scale(frames.radial[index], layout.corner_u[a]),
        scale(frames.toroidal[index], layout.corner_v[b]),
    );
    add(coils.points[index], shift)
}

/// Field of the filament bundles of the coils flagged in `close`.
fn bundle_field(
    point: Vec3,
    coils: &SampledCoils,
    filaments: usize,
    bundle: &[Vec3],
    close: &[bool],
) -> Result<Vec3, Error> {
    ensure(close.len() == coils.n_coils, BAD, "one flag per coil is needed")?;
    ensure(filaments >= 1, BAD, "at least one filament is needed")?;
    let n = coils.n_points;
    let mut total = [0.0; 3];
    for (coil, flagged) in close.iter().enumerate() {
        if *flagged {
            let share = coils.currents[coil] / real(filaments)?;
            for filament in 0..filaments {
                let start = (coil * filaments + filament) * n;
                add_polygon_field(point, &bundle[start..start + n], share, &mut total)?;
            }
        }
    }
    Ok(total)
}

/// Field of distant coils at the four outer corners of a cross-section.
///
/// That field varies slowly over the cross-section, so it is evaluated at
/// the outer corners only and interpolated bilinearly.
fn outer_fields(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    index: usize,
    far: &[bool],
) -> Result<[Vec3; 4], Error> {
    ensure(index < coils.points.len(), BAD, "a sample index is out of range")?;
    ensure(far.len() == coils.n_coils, BAD, "one flag per coil is needed")?;
    let (nr, nt) = (layout.n_radial, layout.n_toroidal);
    Ok([
        filament_field(corner_point(coils, frames, layout, index, 0, 0), coils, far)?,
        filament_field(corner_point(coils, frames, layout, index, nr, 0), coils, far)?,
        filament_field(corner_point(coils, frames, layout, index, 0, nt), coils, far)?,
        filament_field(corner_point(coils, frames, layout, index, nr, nt), coils, far)?,
    ])
}

/// Field at the field points of the cross-section at sample `index` of a
/// coil whose neighbours are flagged in `close` and `far`. The points and
/// the field are written to `out`; the largest field strength is returned.
///
/// # Errors
/// Returns an error if a buffer is too small or a field is not finite.
pub fn station_field(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    bundle: &[Vec3],
    masks: (&[bool], &[bool]),
    index: usize,
    out: &mut [FieldSample],
) -> Result<f64, Error> {
    let (nr, nt) = (layout.n_radial, layout.n_toroidal);
    let stride = nt + 1;
    ensure(out.len() >= (nr + 1) * stride, BAD, "the station buffer is too small")?;
    let outer = outer_fields(coils, frames, layout, index, masks.1)?;
    let coil = index / coils.n_points;
    let tangent = coils.tangents[index];
    let density = MU0 * coils.currents[coil] / layout.wp.area / TAU;
    let (lift_u, lift_v) =
        (cross(tangent, frames.radial[index]), cross(tangent, frames.toroidal[index]));
    let mut peak = 0.0_f64;
    for a in 0..=nr {
        let wu = real(a)? / real(nr)?;
        for b in 0..=nt {
            let wv = real(b)? / real(nt)?;
            let point = corner_point(coils, frames, layout, index, a, b);
            let mut field = bundle_field(point, coils, nr * nt, bundle, masks.0)?;
            let low = add(scale(outer[0], 1.0 - wu), scale(outer[1], wu));
            let high = add(scale(outer[2], 1.0 - wu), scale(outer[3], wu));
            field = add(field, add(scale(low, 1.0 - wv), scale(high, wv)));
            let local =
                add(scale(lift_u, layout.delta_u[a][b]), scale(lift_v, layout.delta_v[a][b]));
            field = add(field, scale(local, density));
            ensure(is_finite(field), BAD, "the field in the winding pack is not finite")?;
            out[a * stride + b] = FieldSample { point, field };
            peak = peak.max(norm(field));
        }
    }
    Ok(peak)
}

/// Field in the reference winding pack of every class representative, at
/// every sample. Fills `work.samples` and `work.station_peak`
/// and returns the largest field strength found [T].
///
/// # Errors
/// Returns an error if the buffers do not match, or a field cannot be
/// evaluated.
pub fn reference_pack_field(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    classes: &SymmetryClasses,
    work: &mut PackWork,
) -> Result<f64, Error> {
    let n = coils.n_points;
    let corners = (layout.n_radial + 1) * (layout.n_toroidal + 1);
    ensure(
        work.samples.len() >= classes.n_classes * n * corners,
        Code::Limit,
        "the pack buffer is too small",
    )?;
    ensure(classes.n_classes <= coils.n_coils, BAD, "more classes than coils")?;
    build_bundles(coils, frames, layout, &mut work.bundle)?;
    let reach = NEAR_FACTOR * layout.wp.diagonal;
    let mut peak = 0.0_f64;
    for class in 0..classes.n_classes {
        let coil = classes.representatives[class];
        let row = &work.distances[class * coils.n_coils..(class + 1) * coils.n_coils];
        select_close(row, coil, reach, &mut work.close, &mut work.far)?;
        for station in 0..n {
            let span = (class * n + station) * corners..(class * n + station + 1) * corners;
            let station_peak = station_field(
                coils,
                frames,
                layout,
                &work.bundle,
                (&work.close[..], &work.far[..]),
                coil * n + station,
                &mut work.samples[span],
            )?;
            work.station_peak[class * n + station] = station_peak;
            peak = peak.max(station_peak);
        }
    }
    Ok(peak)
}

/// Peak field for another winding pack size, evaluated only at the stations
/// where the reference pack came within 80 % of its peak, plus the strongest
/// station of each class. `work.station_peak` must hold the reference result.
///
/// # Errors
/// Returns an error if the buffers do not match, or a field cannot be
/// evaluated.
pub fn scan_pack_peak(
    coils: &SampledCoils,
    frames: &Frames,
    layout: &PackLayout,
    classes: &SymmetryClasses,
    work: &mut PackWork,
    reference_peak: f64,
) -> Result<f64, Error> {
    let n = coils.n_points;
    ensure(reference_peak > 0.0, BAD, "the reference peak field must be positive")?;
    ensure(
        work.station_peak.len() >= classes.n_classes * n,
        Code::Limit,
        "the pack buffer is too small",
    )?;
    build_bundles(coils, frames, layout, &mut work.bundle)?;
    let mut samples = [FieldSample::default(); MAX_CORNERS * MAX_CORNERS];
    let mut peak = 0.0_f64;
    for class in 0..classes.n_classes {
        let coil = classes.representatives[class];
        let row = &work.distances[class * coils.n_coils..(class + 1) * coils.n_coils];
        select_close(row, coil, NEAR_FACTOR * layout.wp.diagonal, &mut work.close, &mut work.far)?;
        let mut strongest = 0.0_f64;
        for value in &work.station_peak[class * n..(class + 1) * n] {
            strongest = strongest.max(*value);
        }
        for station in 0..n {
            let reference = work.station_peak[class * n + station];
            if reference >= 0.8 * reference_peak || reference >= strongest {
                let masks = (&work.close[..], &work.far[..]);
                let out = &mut samples[..];
                let index = coil * n + station;
                peak = peak.max(station_field(
                    coils,
                    frames,
                    layout,
                    &work.bundle,
                    masks,
                    index,
                    out,
                )?);
            }
        }
    }
    Ok(peak)
}

/// Signed distance from `point` to the cross-section of coil `other` nearest
/// to it: negative inside the winding pack.
fn section_distance(
    point: Vec3,
    coils: &SampledCoils,
    frames: &Frames,
    wp: &WindingPack,
    other: usize,
) -> f64 {
    let n = coils.n_points;
    let (mut nearest, mut smallest) = (other * n, f64::INFINITY);
    for index in other * n..(other + 1) * n {
        let offset = sub(point, coils.points[index]);
        if dot(offset, offset) < smallest {
            smallest = dot(offset, offset);
            nearest = index;
        }
    }
    let offset = sub(point, coils.points[nearest]);
    let du = dot(offset, frames.radial[nearest]).abs() - 0.5 * wp.radial;
    let dv = dot(offset, frames.toroidal[nearest]).abs() - 0.5 * wp.toroidal;
    du.max(0.0).hypot(dv.max(0.0)) + du.max(dv).min(0.0)
}

/// Smallest gap between the winding packs of different coils [m].
///
/// For every corner of the winding pack of each class representative, the
/// nearest cross-section of every nearby coil is located and the distance
/// from the corner to that rectangle is measured. A value of zero or less
/// means two winding packs touch or overlap. Because only corners are
/// tested, the result is slightly optimistic for packs that meet edge to
/// edge. `distances` holds the centre-filament distances from each
/// representative to every coil.
///
/// # Errors
/// Returns an error if `distances` does not hold one row per class.
pub fn pack_clearance(
    coils: &SampledCoils,
    frames: &Frames,
    wp: &WindingPack,
    classes: &SymmetryClasses,
    distances: &[f64],
) -> Result<f64, Error> {
    let n = coils.n_points;
    ensure(
        distances.len() >= classes.n_classes * coils.n_coils,
        BAD,
        "one distance row per class is needed",
    )?;
    ensure(
        frames.radial.len() == coils.points.len(),
        BAD,
        "the frame buffer does not match the coils",
    )?;
    let mut clearance = f64::INFINITY;
    for class in 0..classes.n_classes {
        let coil = classes.representatives[class];
        let row = &distances[class * coils.n_coils..(class + 1) * coils.n_coils];
        let mut closest = f64::INFINITY;
        for distance in row {
            closest = closest.min(*distance);
        }
        // Nearby coils, and always the nearest one.
        let reach = (2.0 * wp.diagonal).max(1.01 * closest);
        for (other, distance) in row.iter().enumerate() {
            if *distance <= reach {
                for index in coil * n..(coil + 1) * n {
                    for corner in 0..4 {
                        let side_u = if corner < 2 { -0.5 } else { 0.5 };
                        let side_v = if corner % 2 == 0 { -0.5 } else { 0.5 };
                        let shift = add(
                            scale(frames.radial[index], side_u * wp.radial),
                            scale(frames.toroidal[index], side_v * wp.toroidal),
                        );
                        let point = add(coils.points[index], shift);
                        clearance =
                            clearance.min(section_distance(point, coils, frames, wp, other));
                    }
                }
            }
        }
    }
    Ok(clearance)
}

/// PROCESS's peak-field formula (`bmax_from_awp`) [T].
#[must_use]
pub fn process_bmax(scale_factor: f64, r_major: f64, a1: f64, a2: f64, radial_width: f64) -> f64 {
    scale_factor * (a1 + a2 * r_major / radial_width)
}

/// Least-squares fit of `a1` and `a2` in PROCESS's peak-field formula
/// `B = 0.2 * I[MA] * N / (R - a) * (a1 + a2 * R / width)`.
///
/// `sizes` holds the radial width and the computed peak field of each
/// winding pack size; `count` of them are used.
///
/// # Errors
/// Returns an error if fewer than two sizes are given, a width or a peak
/// field is not positive, or the coil radii are not usable.
pub fn fit_peak_field(
    sizes: (&[f64; MAX_SCAN], &[f64; MAX_SCAN]),
    count: usize,
    current_ma: f64,
    n_coils: usize,
    radii: (f64, f64),
) -> Result<PeakFieldFit, Error> {
    let (widths, peaks) = sizes;
    let (r_major, r_minor) = radii;
    ensure(
        (2..=MAX_SCAN).contains(&count),
        Code::Settings,
        "at least two winding pack sizes are needed for the fit",
    )?;
    ensure(current_ma > 0.0 && n_coils >= 1, BAD, "the fit needs current and at least one coil")?;
    ensure(
        r_major > r_minor && r_minor > 0.0,
        BAD,
        "the coil major radius must exceed its minor radius",
    )?;
    let factor = 0.2 * current_ma * real(n_coils)? / (r_major - r_minor);
    let (mut x_mean, mut y_mean) = (0.0, 0.0);
    for index in 0..count {
        ensure(
            widths[index] > 0.0 && peaks[index] > 0.0,
            BAD,
            "widths and peak fields must be positive",
        )?;
        x_mean += r_major / widths[index] / real(count)?;
        y_mean += peaks[index] / factor / real(count)?;
    }
    let (mut covariance, mut variance) = (0.0, 0.0);
    for index in 0..count {
        let x = r_major / widths[index] - x_mean;
        covariance += x * (peaks[index] / factor - y_mean);
        variance += x * x;
    }
    ensure(variance > 0.0, BAD, "the winding pack sizes for the fit must differ")?;
    let a2 = covariance / variance;
    let a1 = y_mean - a2 * x_mean;
    let mut fit = PeakFieldFit {
        a1,
        a2,
        count,
        widths: *widths,
        peaks: *peaks,
        fitted: [0.0; MAX_SCAN],
        max_relative_error: 0.0,
    };
    for index in 0..count {
        fit.fitted[index] = process_bmax(factor, r_major, a1, a2, widths[index]);
        fit.max_relative_error =
            fit.max_relative_error.max((fit.fitted[index] / peaks[index] - 1.0).abs());
    }
    Ok(fit)
}
