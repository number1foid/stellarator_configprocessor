//! Run every calculation and assemble the PROCESS configuration.
//!
//! All memory is allocated by [`new_workspace`]. [`run`] works in those
//! buffers and allocates nothing.

use crate::coilgeom::{
    coil_distances, coil_measures, fitted_slope, lofted_surface_area, new_loft_work,
    new_symmetry_classes, plasma_coil_distance, port_width, select_near_period, symmetry_classes,
    CoilMeasures, LoftWork, SymmetryClasses,
};
use crate::coils::{new_sampled_coils, sample_coils, SampledCoils};
use crate::consistency::{
    axis_field, check_consistency, enforce, new_axis_field, AxisField, Consistency, Tolerances,
    DEFAULT_TOLERANCES,
};
use crate::curve::{new_spline_work, SplineWork};
use crate::energy::{magnetic_energy, new_inductance_work, EnergyResult, InductanceWork};
use crate::equilibrium::{
    new_coefficients, scale_equilibrium, Coefficients, Equilibrium, Reported,
};
use crate::error::{ensure, Code, Error};
use crate::forces::{coil_forces, ForceResult};
use crate::limits::{
    FIT_SURFACES, MAX_INDUCTANCE_SAMPLES, MAX_SAMPLES, MAX_SCAN, MAX_SURFACE, MIN_SAMPLES,
};
use crate::makegrid::{scale_coils, CoilSet, FileNotes};
use crate::num::real;
use crate::plasma::{surface_geometry, SurfaceGeometry};
use crate::schema::{set_value, validate, Config, EMPTY_CONFIG};
use crate::winding_pack::{
    coil_frames, fit_peak_field, new_frames, new_pack_work, pack_clearance, pack_layout,
    reference_pack_field, scan_pack_peak, winding_pack, Frames, PackLayout, PackWork, PeakFieldFit,
    WindingPack,
};
use std::io::Write;

/// Number of kinds of warning.
pub const N_WARNINGS: usize = 10;
/// Warning: coils with zero current were ignored.
pub const WARN_ZERO_CURRENT: usize = 0;
/// Warning: a closing point in the coil file has a non-zero current.
pub const WARN_CLOSING_CURRENT: usize = 1;
/// Warning: the current varies along a coil in the coil file.
pub const WARN_VARYING_CURRENT: usize = 2;
/// Warning: a coil in the coil file is not closed.
pub const WARN_NOT_CLOSED: usize = 3;
/// Warning: the equilibrium is not stellarator-symmetric.
pub const WARN_ASYMMETRIC: usize = 4;
/// Warning: the consistency check failed and was overridden.
pub const WARN_OVERRIDDEN: usize = 5;
/// Warning: half the winding pack width exceeds the plasma-coil distance.
pub const WARN_PACK_EXCEEDS_DISTANCE: usize = 6;
/// Warning: coil currents differ between coils.
pub const WARN_CURRENTS_DIFFER: usize = 7;
/// Warning: winding packs of neighbouring coils overlap.
pub const WARN_OVERLAP: usize = 8;
/// Warning: winding pack sizes were left out of the peak-field fit.
pub const WARN_SIZES_SKIPPED: usize = 9;

/// Factor from SI units to mega-units (MA, MN).
const MEGA: f64 = 1.0e-6;
/// Flux surfaces used for the derivative of the plasma-coil distance.
const FIT_FLUX: [f64; FIT_SURFACES] = [1.0, 0.95, 0.9, 0.85];

/// Run settings. The first four have no default on purpose.
#[derive(Debug, Clone, Copy)]
pub struct Settings {
    /// Winding pack width in the coil's in-plane normal direction `[m]`.
    pub wp_radial: f64,
    /// Winding pack width in the toroidal direction `[m]`.
    pub wp_toroidal: f64,
    /// Effective helical ripple, passed through to PROCESS (not computed).
    pub epseff: f64,
    /// Neutron wall load peaking factor, passed through (not computed).
    pub neutron_peakfactor: f64,
    /// Rescale the machine to this plasma major radius `[m]`.
    pub rmajor_ref: Option<f64>,
    /// Rescale the coil currents to give this mean field on axis `[T]`.
    pub bt_ref: Option<f64>,
    /// Scale the coil currents to the equilibrium's field before checking.
    pub match_field: bool,
    /// Continue when coils and equilibrium appear inconsistent.
    pub allow_inconsistent: bool,
    /// Limits of the consistency tests.
    pub tolerances: Tolerances,
    /// Samples per coil.
    pub n_points: usize,
    /// Filaments across the winding pack (radial, toroidal).
    pub wp_filaments: (usize, usize),
    /// Number of winding pack sizes used for the a1, a2 fit.
    pub n_scan: usize,
    /// Those sizes, relative to the reference. The default only shrinks the
    /// pack: larger packs often touch the neighbouring coils, and PROCESS at
    /// reactor scale works with a pack that is smaller relative to the coil
    /// than the reference.
    pub scan_factors: [f64; MAX_SCAN],
    /// Angle resolution of the plasma surface for distance calculations.
    pub surface_resolution: usize,
}

/// Every buffer the calculation needs.
#[derive(Debug)]
pub struct Workspace {
    /// Spline fitting and sampling.
    pub spline: SplineWork,
    /// The sampled coils.
    pub coils: SampledCoils,
    /// Fourier coefficients of one flux surface.
    pub surface: Coefficients,
    /// Coil field along the magnetic axis.
    pub axis: AxisField,
    /// Symmetry classes of the coils.
    pub classes: SymmetryClasses,
    /// Winding pack orientation at every sample.
    pub frames: Frames,
    /// Field in the winding pack.
    pub pack: PackWork,
    /// Self-inductance integral.
    pub inductance: InductanceWork,
    /// Per coil: can be closest to the first field period.
    pub near: Vec<bool>,
    /// Coil indices in toroidal order.
    pub order: Vec<usize>,
    /// Surface lofted between neighbouring coils.
    pub loft: LoftWork,
    /// Self-inductance of each class representative `[H]`.
    pub class_self_inductance: Vec<f64>,
    /// Net force towards the torus axis on each class representative `[N]`.
    pub class_centering: Vec<f64>,
    /// The latest consistency check, kept for the error report.
    pub last_check: Option<Consistency>,
}

/// Factors applied to the inputs to reach the reference point.
#[derive(Debug, Clone, Copy)]
pub struct Scaling {
    /// Factor on all lengths.
    pub length_factor: f64,
    /// Factor on the field.
    pub field_factor: f64,
    /// Factor on the currents that makes the coil field match the equilibrium.
    pub current_factor: f64,
    /// Factor on the currents in total. Field ~ current / length, so
    /// currents follow both other factors.
    pub total_current_factor: f64,
}

/// Coil geometry values shared by several stages.
#[derive(Debug, Clone, Copy)]
pub struct CoilGeometry {
    /// Overall dimensions of the coil set.
    pub measures: CoilMeasures,
    /// Minimum distance between the centre filaments of two coils `[m]`.
    pub d_min: f64,
    /// Largest distance between the centroids of neighbouring coils `[m]`.
    pub port_width: f64,
    /// Area of the surface spanned by straight lines between neighbouring
    /// coils [m^2].
    pub lofted_surface: f64,
    /// Effective minor radius of the flux surfaces used for the fit `[m]`.
    pub fit_radii: [f64; FIT_SURFACES],
    /// Minimum plasma-coil distance for those surfaces `[m]`.
    pub fit_distances: [f64; FIT_SURFACES],
    /// Derivative of the distance with respect to the minor radius.
    pub slope: f64,
}

/// Everything a run produces, apart from the per-class arrays that stay in
/// the workspace.
#[derive(Debug, Clone, Copy)]
pub struct Outcome {
    /// The values PROCESS reads.
    pub config: Config,
    /// Consistency of the inputs as given.
    pub consistency: Consistency,
    /// Factors applied to the inputs.
    pub scaling: Scaling,
    /// Geometry of the last closed flux surface.
    pub plasma: SurfaceGeometry,
    /// Scalars as reported by the equilibrium file, after scaling.
    pub reported: Reported,
    /// Coil geometry.
    pub geometry: CoilGeometry,
    /// Mean toroidal component of the coil field on the axis `[T]`.
    pub b_axis_toroidal: f64,
    /// Smallest gap between the winding packs of different coils `[m]`.
    pub clearance: f64,
    /// Peak field in the reference winding pack `[T]`.
    pub peak_field: f64,
    /// Fit of the peak field against winding pack size.
    pub fit: PeakFieldFit,
    /// Magnetic energy.
    pub energy: EnergyResult,
    /// Forces on the coils.
    pub forces: ForceResult,
    /// Irregularities met while reading the coil file.
    pub notes: FileNotes,
    /// For each kind of warning, whether it applies; see the `WARN_` constants.
    pub warnings: [bool; N_WARNINGS],
    /// For each scan factor, whether it was left out of the fit.
    pub skipped: [bool; MAX_SCAN],
}

/// Settings with the default resolution, tolerances and scan sizes.
#[must_use]
pub const fn default_settings(wp: (f64, f64), epseff: f64, neutron_peakfactor: f64) -> Settings {
    Settings {
        wp_radial: wp.0,
        wp_toroidal: wp.1,
        epseff,
        neutron_peakfactor,
        rmajor_ref: None,
        bt_ref: None,
        match_field: false,
        allow_inconsistent: false,
        tolerances: DEFAULT_TOLERANCES,
        n_points: 160,
        wp_filaments: (4, 4),
        n_scan: 5,
        scan_factors: [0.5, 0.6, 0.7, 0.85, 1.0, 0.0, 0.0, 0.0],
        surface_resolution: 96,
    }
}

/// Check the run settings.
///
/// # Errors
/// Returns an error naming the first setting that is missing or invalid.
pub fn validate_settings(settings: &Settings) -> Result<(), Error> {
    let bad = Code::Settings;
    ensure(
        settings.wp_radial.is_finite() && settings.wp_radial > 0.0,
        bad,
        "wp_radial must be a positive number",
    )?;
    ensure(
        settings.wp_toroidal.is_finite() && settings.wp_toroidal > 0.0,
        bad,
        "wp_toroidal must be a positive number",
    )?;
    ensure(
        settings.epseff.is_finite() && settings.epseff >= 0.0,
        bad,
        "epseff must be a non-negative number",
    )?;
    ensure(
        settings.neutron_peakfactor.is_finite() && settings.neutron_peakfactor > 0.0,
        bad,
        "neutron_peakfactor must be a positive number",
    )?;
    let radius = settings.rmajor_ref.unwrap_or(1.0);
    let field = settings.bt_ref.unwrap_or(1.0);
    ensure(radius.is_finite() && radius > 0.0, bad, "rmajor_ref must be positive")?;
    ensure(field.is_finite() && field > 0.0, bad, "bt_ref must be positive")?;
    ensure(
        (MIN_SAMPLES..=MAX_SAMPLES).contains(&settings.n_points),
        bad,
        "points must be between 32 and 1024",
    )?;
    ensure(
        (8..=MAX_SURFACE).contains(&settings.surface_resolution),
        bad,
        "the surface resolution must be 8 to 512",
    )?;
    ensure(
        (2..=MAX_SCAN).contains(&settings.n_scan),
        bad,
        "wp_scan_factors must contain 1.0 and at least one other size",
    )?;
    let mut reference = 0;
    for factor in &settings.scan_factors[..settings.n_scan] {
        ensure(factor.is_finite() && *factor > 0.0, bad, "wp_scan_factors must be positive")?;
        reference += usize::from((factor - 1.0).abs() < 1.0e-12);
    }
    ensure(reference == 1, bad, "wp_scan_factors must contain 1.0 exactly once")?;
    Ok(())
}

/// Allocate every buffer needed to process the given inputs.
///
/// # Errors
/// Returns an error if the settings are invalid or a size exceeds the limits.
pub fn new_workspace(
    eq: &Equilibrium,
    set: &CoilSet,
    settings: &Settings,
) -> Result<Workspace, Error> {
    ensure(eq.n_modes >= 1 && !set.currents.is_empty(), Code::Numeric, "the inputs are empty")?;
    validate_settings(settings)?;
    let n_coils = set.currents.len();
    ensure(set.offsets.len() == n_coils + 1, Code::CoilFile, "the coil set is damaged")?;
    let mut longest = settings.n_points;
    for coil in 0..n_coils {
        longest = longest.max(set.offsets[coil + 1] - set.offsets[coil]);
    }
    Ok(Workspace {
        spline: new_spline_work(longest, MAX_INDUCTANCE_SAMPLES)?,
        coils: new_sampled_coils(n_coils, settings.n_points)?,
        surface: new_coefficients(eq.n_modes)?,
        axis: new_axis_field(n_coils)?,
        classes: new_symmetry_classes(n_coils)?,
        frames: new_frames(n_coils, settings.n_points)?,
        pack: new_pack_work(n_coils, settings.n_points, settings.wp_filaments)?,
        inductance: new_inductance_work()?,
        near: vec![false; n_coils],
        order: vec![0; n_coils],
        loft: new_loft_work(n_coils, settings.n_points)?,
        class_self_inductance: vec![0.0; n_coils],
        class_centering: vec![0.0; n_coils],
        last_check: None,
    })
}

/// Length, field and current factors that bring the machine to the
/// reference point.
fn scaling_factors(
    eq: &Equilibrium,
    check: &Consistency,
    settings: &Settings,
    surface: &mut Coefficients,
) -> Result<Scaling, Error> {
    ensure(
        check.b_axis_coils > 0.0,
        Code::Settings,
        "the coils give no field on the axis; cannot rescale",
    )?;
    let current_factor =
        if settings.match_field { check.b_axis_equilibrium / check.b_axis_coils } else { 1.0 };
    ensure(
        current_factor.is_finite() && current_factor > 0.0,
        Code::Settings,
        "the equilibrium has no field on the axis; cannot match the field",
    )?;
    let now = check.b_axis_coils * current_factor;
    let field_factor = settings.bt_ref.unwrap_or(now) / now;
    let radius = surface_geometry(eq, 1.0, surface)?.major_radius;
    let length_factor = settings.rmajor_ref.unwrap_or(radius) / radius;
    Ok(Scaling {
        length_factor,
        field_factor,
        current_factor,
        total_current_factor: current_factor * field_factor * length_factor,
    })
}

/// Sample the coils, check them against the equilibrium and rescale both.
/// Returns the check of the inputs as given, the check after scaling and
/// the scaling.
fn check_inputs<W: Write>(
    eq: &mut Equilibrium,
    set: &mut CoilSet,
    settings: &Settings,
    work: &mut Workspace,
    log: &mut W,
) -> Result<(Consistency, Consistency, Scaling), Error> {
    ensure(set.periods >= 1, Code::CoilFile, "the coil file declares no field periods")?;
    writeln!(log, "... sampling coils")?;
    sample_coils(set, &mut work.spline, &mut work.coils)?;
    writeln!(log, "... checking that coils and equilibrium match")?;
    axis_field(eq, &work.coils, &mut work.surface, &mut work.axis)?;
    let compare = !settings.match_field;
    let before =
        check_consistency(eq, &work.coils, set.periods, &settings.tolerances, &work.axis, compare)?;
    work.last_check = Some(before);
    enforce(&before, settings.allow_inconsistent)?;
    let scaling = scaling_factors(eq, &before, settings, &mut work.surface)?;
    let mut after = before;
    let changed = (scaling.length_factor - 1.0).abs() > 0.0
        || (scaling.total_current_factor - 1.0).abs() > 0.0;
    if changed {
        scale_equilibrium(eq, scaling.length_factor, scaling.field_factor)?;
        scale_coils(set, scaling.length_factor, scaling.total_current_factor)?;
        sample_coils(set, &mut work.spline, &mut work.coils)?;
        axis_field(eq, &work.coils, &mut work.surface, &mut work.axis)?;
        after = check_consistency(
            eq,
            &work.coils,
            set.periods,
            &settings.tolerances,
            &work.axis,
            false,
        )?;
    }
    ensure(
        work.coils.n_coils % eq.nfp == 0,
        Code::Settings,
        "the coils cannot be divided into the field periods of the equilibrium",
    )?;
    Ok((before, after, scaling))
}

/// Symmetry classes, coil dimensions and distances.
fn coil_geometry(
    eq: &Equilibrium,
    settings: &Settings,
    work: &mut Workspace,
) -> Result<CoilGeometry, Error> {
    let count = work.coils.n_coils;
    ensure(
        work.pack.distances.len() >= count * count,
        Code::Limit,
        "the distance buffer is too small",
    )?;
    symmetry_classes(&work.coils, eq.nfp, !eq.lasym, &mut work.classes)?;
    let mut d_min = f64::INFINITY;
    for class in 0..work.classes.n_classes {
        let row = &mut work.pack.distances[class * count..(class + 1) * count];
        coil_distances(&work.coils, work.classes.representatives[class], row)?;
        for distance in row.iter() {
            d_min = d_min.min(*distance);
        }
    }
    select_near_period(&work.coils, eq.nfp, &mut work.near)?;
    let (mut fit_radii, mut fit_distances) = ([0.0; FIT_SURFACES], [0.0; FIT_SURFACES]);
    for (index, flux) in FIT_FLUX.iter().enumerate() {
        fit_radii[index] = surface_geometry(eq, *flux, &mut work.surface)?.minor_radius;
        fit_distances[index] = plasma_coil_distance(
            eq,
            &work.coils,
            *flux,
            settings.surface_resolution,
            &mut work.surface,
            &work.near,
        )?;
    }
    ensure(d_min > 0.0, Code::CoilFile, "two coils touch or cross")?;
    let port = port_width(&work.coils, &mut work.order)?;
    Ok(CoilGeometry {
        measures: coil_measures(&work.coils)?,
        d_min,
        port_width: port,
        lofted_surface: lofted_surface_area(&work.coils, &work.order, &mut work.loft)?,
        fit_radii,
        fit_distances,
        slope: fitted_slope(&fit_radii, &fit_distances)?,
    })
}

/// Peak field for several winding pack sizes and the fit of a1 and a2.
fn peak_field_fit(
    settings: &Settings,
    layout: &PackLayout,
    reference_peak: f64,
    measures: &CoilMeasures,
    work: &mut Workspace,
    skipped: &mut [bool; MAX_SCAN],
) -> Result<PeakFieldFit, Error> {
    let (nr, nt) = settings.wp_filaments;
    ensure(settings.n_scan <= MAX_SCAN, Code::Settings, "too many winding pack sizes")?;
    let (mut widths, mut peaks, mut count) = ([0.0; MAX_SCAN], [0.0; MAX_SCAN], 0);
    for (index, factor) in settings.scan_factors[..settings.n_scan].iter().enumerate() {
        let trial = winding_pack(layout.wp.radial * factor, layout.wp.toroidal * factor)?;
        // A trial pack must stay clear of its neighbours by one filament cell,
        // otherwise the neighbour's discrete filaments distort the field.
        let cell = (trial.radial / real(nr)?).max(trial.toroidal / real(nt)?);
        let reference = (factor - 1.0).abs() < 1.0e-12;
        let gap =
            pack_clearance(&work.coils, &work.frames, &trial, &work.classes, &work.pack.distances)?;
        skipped[index] = !reference && gap < cell;
        if !skipped[index] {
            let scaled = pack_layout(&trial, nr, nt)?;
            widths[count] = trial.radial;
            peaks[count] = if reference {
                reference_peak
            } else {
                scan_pack_peak(
                    &work.coils,
                    &work.frames,
                    &scaled,
                    &work.classes,
                    &mut work.pack,
                    reference_peak,
                )?
            };
            count += 1;
        }
    }
    ensure(
        count >= 2,
        Code::Settings,
        "fewer than two usable winding pack sizes for the a1, a2 fit; use smaller sizes",
    )?;
    let mut current_sum = 0.0;
    for current in &work.coils.currents {
        current_sum += current.abs();
    }
    let current_ma = current_sum / real(work.coils.n_coils)? / 1.0e6;
    fit_peak_field(
        (&widths, &peaks),
        count,
        current_ma,
        work.coils.n_coils,
        (measures.r_major, measures.r_minor),
    )
}

/// Decide which warnings apply.
fn collect_warnings(
    eq: &Equilibrium,
    wp: &WindingPack,
    coils: &SampledCoils,
    outcome: &mut Outcome,
) -> Result<(), Error> {
    ensure(coils.n_coils >= 1, Code::Numeric, "at least one coil is needed")?;
    let (mut weakest, mut strongest) = (f64::INFINITY, 0.0_f64);
    for current in &coils.currents {
        weakest = weakest.min(current.abs());
        strongest = strongest.max(current.abs());
    }
    ensure(strongest > 0.0, Code::CoilFile, "no coil carries current")?;
    let mut any_skipped = false;
    for flag in &outcome.skipped {
        any_skipped = any_skipped || *flag;
    }
    let notes = outcome.notes;
    outcome.warnings[WARN_ZERO_CURRENT] = notes.zero_current > 0;
    outcome.warnings[WARN_CLOSING_CURRENT] = notes.closing_current > 0;
    outcome.warnings[WARN_VARYING_CURRENT] = notes.varying_current > 0;
    outcome.warnings[WARN_NOT_CLOSED] = notes.not_closed > 0;
    outcome.warnings[WARN_ASYMMETRIC] = eq.lasym;
    outcome.warnings[WARN_OVERRIDDEN] = !outcome.consistency.passed;
    outcome.warnings[WARN_PACK_EXCEEDS_DISTANCE] =
        0.5 * wp.radial >= outcome.geometry.fit_distances[0];
    outcome.warnings[WARN_CURRENTS_DIFFER] = strongest - weakest > 1.0e-6 * strongest;
    outcome.warnings[WARN_OVERLAP] = outcome.clearance <= 0.0;
    outcome.warnings[WARN_SIZES_SKIPPED] = any_skipped;
    Ok(())
}

/// Fill the values PROCESS reads, keyed as in the file, and validate them.
fn fill_config(
    eq: &Equilibrium,
    settings: &Settings,
    wp: &WindingPack,
    work: &Workspace,
    checked: &Consistency,
    outcome: &mut Outcome,
) -> Result<(), Error> {
    ensure(eq.nfp >= 1, Code::Equilibrium, "the equilibrium has no field periods")?;
    ensure(wp.toroidal > 0.0, Code::Settings, "wp_toroidal must be a positive number")?;
    let (plasma, coil, forces) = (outcome.plasma, outcome.geometry, outcome.forces);
    let config = &mut outcome.config;
    set_value(config, "symmetry", real(eq.nfp)?)?;
    set_value(config, "coilspermodule", real(work.coils.n_coils / eq.nfp)?)?;
    set_value(config, "rmajor_ref", plasma.major_radius)?;
    set_value(config, "rminor_ref", plasma.minor_radius)?;
    set_value(config, "aspect_ref", plasma.aspect_ratio)?;
    set_value(config, "bt_ref", checked.b_axis_coils)?;
    set_value(config, "vol_plasma", plasma.volume)?;
    set_value(config, "plasma_surface", plasma.area)?;
    set_value(config, "coil_rmajor", coil.measures.r_major)?;
    set_value(config, "coil_rminor", coil.measures.r_minor)?;
    set_value(config, "coillength", coil.measures.total_length)?;
    set_value(config, "coilsurface", coil.measures.surface_area)?;
    set_value(config, "maximal_coil_height", coil.measures.height)?;
    set_value(config, "dmin", coil.d_min)?;
    set_value(config, "max_portsize_width", coil.port_width)?;
    set_value(config, "min_bend_radius", coil.measures.min_bend_radius)?;
    set_value(config, "min_plasma_coil_distance", coil.fit_distances[0])?;
    set_value(config, "derivative_min_LCFS_coils_dist", coil.slope)?;
    set_value(config, "WP_area", wp.area)?;
    set_value(config, "WP_ratio", wp.radial / wp.toroidal)?;
    set_value(config, "i0", outcome.energy.reference_current * MEGA)?;
    set_value(config, "WP_bmax", outcome.peak_field)?;
    set_value(config, "a1", outcome.fit.a1)?;
    set_value(config, "a2", outcome.fit.a2)?;
    set_value(config, "inductance", outcome.energy.effective_inductance)?;
    set_value(config, "max_force_density", forces.max_force_density * MEGA)?;
    set_value(config, "max_force_density_MNm", forces.max_line_force * MEGA)?;
    set_value(config, "max_lateral_force_density", forces.max_lateral_force_density * MEGA)?;
    set_value(config, "max_radial_force_density", forces.max_radial_force_density * MEGA)?;
    set_value(config, "centering_force_max_MN", forces.centering_max * MEGA)?;
    set_value(config, "centering_force_min_MN", forces.centering_min * MEGA)?;
    set_value(config, "centering_force_avg_MN", forces.centering_avg * MEGA)?;
    set_value(config, "epseff", settings.epseff)?;
    set_value(config, "neutron_peakfactor", settings.neutron_peakfactor)?;
    validate(config)
}

/// Results of the electromagnetic stages.
#[derive(Debug, Clone, Copy)]
struct Electromagnetics {
    clearance: f64,
    peak_field: f64,
    forces: ForceResult,
    fit: PeakFieldFit,
    energy: EnergyResult,
    skipped: [bool; MAX_SCAN],
}

/// Field in the winding pack, forces, peak-field fit and magnetic energy.
fn electromagnetics<W: Write>(
    settings: &Settings,
    wp: &WindingPack,
    measures: &CoilMeasures,
    work: &mut Workspace,
    log: &mut W,
) -> Result<Electromagnetics, Error> {
    ensure(wp.area > 0.0, Code::Settings, "the winding pack has no area")?;
    ensure(work.classes.n_classes >= 1, Code::Numeric, "no symmetry class was found")?;
    writeln!(log, "... field in the winding pack")?;
    coil_frames(&work.coils, &mut work.frames)?;
    let layout = pack_layout(wp, settings.wp_filaments.0, settings.wp_filaments.1)?;
    let distances = &work.pack.distances;
    let clearance = pack_clearance(&work.coils, &work.frames, wp, &work.classes, distances)?;
    let peak_field =
        reference_pack_field(&work.coils, &work.frames, &layout, &work.classes, &mut work.pack)?;
    writeln!(log, "... forces")?;
    let centering = &mut work.class_centering;
    let forces =
        coil_forces(&work.coils, &work.frames, &work.classes, &layout, &work.pack, centering)?;
    writeln!(log, "... peak field for several winding pack sizes")?;
    let mut skipped = [false; MAX_SCAN];
    let fit = peak_field_fit(settings, &layout, peak_field, measures, work, &mut skipped)?;
    writeln!(log, "... magnetic energy")?;
    let energy = magnetic_energy(
        &work.coils,
        wp,
        &work.classes,
        &mut work.spline,
        &mut work.inductance,
        &mut work.class_self_inductance,
    )?;
    Ok(Electromagnetics { clearance, peak_field, forces, fit, energy, skipped })
}

/// Compute the configuration for one equilibrium and coil set.
///
/// The equilibrium and the coil set are rescaled in place to the reference
/// point. Progress messages are written to `log`.
///
/// # Errors
/// Returns an error if the settings are invalid, coils and equilibrium are
/// inconsistent and no override is set, a calculation fails, or the
/// resulting configuration does not satisfy the schema.
pub fn run<W: Write>(
    eq: &mut Equilibrium,
    set: &mut CoilSet,
    settings: &Settings,
    work: &mut Workspace,
    log: &mut W,
) -> Result<Outcome, Error> {
    validate_settings(settings)?;
    let matching = work.coils.n_coils == set.currents.len();
    ensure(matching, Code::Numeric, "the workspace does not match the coil set")?;
    let modes = eq.n_modes == work.surface.rc.len();
    ensure(modes, Code::Numeric, "the workspace does not match the equilibrium")?;
    let (consistency, checked, scaling) = check_inputs(eq, set, settings, work, log)?;
    let wp = winding_pack(settings.wp_radial, settings.wp_toroidal)?;
    writeln!(log, "... plasma geometry")?;
    let plasma = surface_geometry(eq, 1.0, &mut work.surface)?;
    writeln!(log, "... coil geometry")?;
    let geometry = coil_geometry(eq, settings, work)?;
    let results = electromagnetics(settings, &wp, &geometry.measures, work, log)?;
    let mut outcome = Outcome {
        config: EMPTY_CONFIG,
        consistency,
        scaling,
        plasma,
        reported: eq.reported,
        geometry,
        b_axis_toroidal: checked.b_axis_toroidal,
        clearance: results.clearance,
        peak_field: results.peak_field,
        fit: results.fit,
        energy: results.energy,
        forces: results.forces,
        notes: set.notes,
        warnings: [false; N_WARNINGS],
        skipped: results.skipped,
    };
    collect_warnings(eq, &wp, &work.coils, &mut outcome)?;
    fill_config(eq, settings, &wp, work, &checked, &mut outcome)?;
    Ok(outcome)
}
