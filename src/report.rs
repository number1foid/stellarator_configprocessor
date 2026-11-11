//! Human-readable and JSON run reports.

use crate::consistency::{
    Consistency, Tolerances, N_TESTS, TEST_ALIGNMENT, TEST_DIVISIBILITY, TEST_FIELD, TEST_LINKAGE,
    TEST_PERIODS,
};
use crate::error::{ensure, Code, Error};
use crate::limits::FIT_SURFACES;
use crate::num::{GENERATOR, SCHEMA_VERSION};
use crate::pipeline::{
    Outcome, Settings, Workspace, N_WARNINGS, WARN_ASYMMETRIC, WARN_CLOSING_CURRENT,
    WARN_CURRENTS_DIFFER, WARN_NOT_CLOSED, WARN_OVERLAP, WARN_OVERRIDDEN,
    WARN_PACK_EXCEEDS_DISTANCE, WARN_SIZES_SKIPPED, WARN_VARYING_CURRENT, WARN_ZERO_CURRENT,
};
use crate::schema::{write_json_real, write_json_text, Kind, KEYS};
use std::io::Write;

const BAD: Code = Code::Numeric;

/// Write a real number with six significant digits.
///
/// # Errors
/// Returns an error if writing fails.
pub fn write_real<W: Write>(out: &mut W, value: f64) -> Result<(), Error> {
    let size = value.abs();
    if !value.is_finite() || size == 0.0 {
        write!(out, "{value}")?;
    } else if (1.0e-4..1.0e6).contains(&size) {
        let (mut decimals, mut bound) = (9_usize, 1.0e-3);
        for _ in 0..9 {
            if size >= bound {
                decimals -= 1;
                bound *= 10.0;
            }
        }
        write!(out, "{value:.decimals$}")?;
    } else {
        write!(out, "{value:.5e}")?;
    }
    Ok(())
}

/// Write the description of failed consistency test `index`.
///
/// # Errors
/// Returns an error if `index` is not a test or writing fails.
pub fn failure_text<W: Write>(
    out: &mut W,
    check: &Consistency,
    tolerances: &Tolerances,
    index: usize,
) -> Result<(), Error> {
    ensure(index < N_TESTS, BAD, "unknown consistency test")?;
    if index == TEST_PERIODS {
        write!(
            out,
            "field periods differ: coil file declares {}, equilibrium has {}",
            check.nfp_coil_file, check.nfp_equilibrium
        )?;
    } else if index == TEST_DIVISIBILITY {
        write!(
            out,
            "{} coils cannot be divided into {} field periods",
            check.n_coils, check.nfp_equilibrium
        )?;
    } else if index == TEST_LINKAGE {
        write!(
            out,
            "the magnetic axis links {:.4} of the total coil current (expected 1): some coils do \
             not enclose the axis, or coil currents have opposing directions",
            check.linkage_ratio
        )?;
    } else if index == TEST_ALIGNMENT {
        write!(
            out,
            "coil field is not tangent to the equilibrium's magnetic axis: mean perpendicular \
             fraction {:.4} exceeds {}; the axis position or shape does not match the coils",
            check.axis_misalignment, tolerances.axis_misalignment
        )?;
    } else if index == TEST_FIELD {
        write!(
            out,
            "field on axis differs: coils give {:.4} T, equilibrium has {:.4} T (ratio {:.4}, \
             tolerance {})",
            check.b_axis_coils, check.b_axis_equilibrium, check.field_ratio, tolerances.field
        )?;
    }
    Ok(())
}

/// Text of the warnings that come from the coil file.
fn file_warning_text<W: Write>(out: &mut W, outcome: &Outcome, index: usize) -> Result<(), Error> {
    let notes = outcome.notes;
    ensure(index <= WARN_NOT_CLOSED, BAD, "not a coil file warning")?;
    if index == WARN_ZERO_CURRENT {
        write!(out, "{} coil(s) with zero current were ignored", notes.zero_current)?;
    } else if index == WARN_CLOSING_CURRENT {
        write!(
            out,
            "{} coil(s) have a closing point with non-zero current",
            notes.closing_current
        )?;
    } else if index == WARN_VARYING_CURRENT {
        write!(
            out,
            "{} coil(s) have a current that varies along the filament; the first value is used",
            notes.varying_current
        )?;
    } else {
        write!(
            out,
            "{} coil(s) are not closed in the coil file; each is closed with a straight segment",
            notes.not_closed
        )?;
    }
    Ok(())
}

/// Write the text of warning `index`.
///
/// # Errors
/// Returns an error if `index` is not a warning or writing fails.
pub fn warning_text<W: Write>(
    out: &mut W,
    settings: &Settings,
    outcome: &Outcome,
    index: usize,
) -> Result<(), Error> {
    ensure(index < N_WARNINGS, BAD, "unknown warning")?;
    if index <= WARN_NOT_CLOSED {
        file_warning_text(out, outcome, index)?;
    } else if index == WARN_ASYMMETRIC {
        write!(out, "the equilibrium is not stellarator-symmetric; coils are grouped by field-period rotation only")?;
    } else if index == WARN_OVERRIDDEN {
        write!(out, "the consistency check failed and was overridden")?;
    } else if index == WARN_PACK_EXCEEDS_DISTANCE {
        write!(out, "half the winding pack radial width exceeds the plasma-coil distance")?;
    } else if index == WARN_CURRENTS_DIFFER {
        write!(out, "coil currents differ between coils; PROCESS scales a single current, taken here as the mean absolute current")?;
    } else if index == WARN_OVERLAP {
        write!(
            out,
            "winding packs of neighbouring coils overlap by {:.4} m (minimum distance between \
             filaments: {:.4} m); the peak field and forces are not meaningful. Use smaller \
             winding pack dimensions.",
            -outcome.clearance, outcome.geometry.d_min
        )?;
    } else if index == WARN_SIZES_SKIPPED {
        write!(out, "winding pack sizes left out of the a1, a2 fit because they would touch neighbouring coils: factors")?;
        for (factor, skipped) in settings.scan_factors.iter().zip(outcome.skipped.iter()) {
            if *skipped {
                write!(out, " {factor}")?;
            }
        }
    }
    Ok(())
}

/// Largest field strength in the winding pack of the representative of
/// class `class` [T].
fn class_peak(work: &Workspace, class: usize) -> Result<f64, Error> {
    let n = work.coils.n_points;
    ensure(class < work.classes.n_classes, BAD, "a class index is out of range")?;
    ensure(work.pack.station_peak.len() >= (class + 1) * n, BAD, "the pack buffer is too small")?;
    let mut peak = 0.0_f64;
    for value in &work.pack.station_peak[class * n..(class + 1) * n] {
        peak = peak.max(*value);
    }
    Ok(peak)
}

fn text_consistency<W: Write>(
    out: &mut W,
    settings: &Settings,
    outcome: &Outcome,
) -> Result<(), Error> {
    let check = outcome.consistency;
    ensure(check.n_coils >= 1, BAD, "the consistency result is empty")?;
    writeln!(out, "\nConsistency of coils and equilibrium")?;
    writeln!(
        out,
        "  field periods: equilibrium {}, coil file {}, {} coils",
        check.nfp_equilibrium, check.nfp_coil_file, check.n_coils
    )?;
    write!(out, "  current linked by the magnetic axis / total current: ")?;
    write_real(out, check.linkage_ratio)?;
    write!(out, "\n  field perpendicular to the axis / field: ")?;
    write_real(out, check.axis_misalignment)?;
    write!(out, "\n  field on axis: coils ")?;
    write_real(out, check.b_axis_coils)?;
    write!(out, " T, equilibrium ")?;
    write_real(out, check.b_axis_equilibrium)?;
    write!(out, " T (ratio ")?;
    write_real(out, check.field_ratio)?;
    writeln!(out, ")\n  result: {}", if check.passed { "passed" } else { "FAILED (overridden)" })?;
    for (index, failed) in check.failed.iter().enumerate() {
        if *failed {
            write!(out, "  - ")?;
            failure_text(out, &check, &settings.tolerances, index)?;
            writeln!(out)?;
        }
    }
    let scaling = outcome.scaling;
    if (scaling.length_factor - 1.0).abs() > 0.0 || (scaling.total_current_factor - 1.0).abs() > 0.0
    {
        write!(out, "\nScaling applied to the inputs\n  lengths x ")?;
        write_real(out, scaling.length_factor)?;
        write!(out, ", field x ")?;
        write_real(out, scaling.field_factor)?;
        write!(out, ", currents x ")?;
        write_real(out, scaling.total_current_factor)?;
        writeln!(out)?;
    }
    Ok(())
}

fn text_values<W: Write>(out: &mut W, name: &str, outcome: &Outcome) -> Result<(), Error> {
    ensure(!name.is_empty(), BAD, "the configuration needs a name")?;
    writeln!(out, "\nValues written to the configuration file")?;
    for (index, entry) in KEYS.iter().enumerate() {
        ensure(
            entry.kind == Kind::Text || outcome.config.assigned[index],
            BAD,
            "a value is missing",
        )?;
        write!(out, "  {:<30}  ", entry.name)?;
        if entry.kind == Kind::Text {
            write!(out, "{name}")?;
        } else if entry.kind == Kind::Count {
            write!(out, "{:.0}", outcome.config.values[index])?;
        } else {
            write_real(out, outcome.config.values[index])?;
        }
        writeln!(out, " {}  ({})", entry.unit, entry.method)?;
    }
    Ok(())
}

fn text_diagnostics<W: Write>(
    out: &mut W,
    outcome: &Outcome,
    work: &Workspace,
) -> Result<(), Error> {
    ensure(work.classes.n_classes >= 1, BAD, "no symmetry class was found")?;
    write!(out, "\nDiagnostics\n  magnetic energy: ")?;
    write_real(out, outcome.energy.energy * 1.0e-9)?;
    writeln!(out, " GJ\n  coil classes (one representative coil each):")?;
    for class in 0..work.classes.n_classes {
        let coil = work.classes.representatives[class];
        write!(out, "    coil {:>3} x{:<3} current ", coil, work.classes.multiplicity[class])?;
        write_real(out, work.coils.currents[coil] * 1.0e-6)?;
        write!(out, " MA, length ")?;
        write_real(out, work.coils.lengths[coil])?;
        write!(out, " m, peak field ")?;
        write_real(out, class_peak(work, class)?)?;
        write!(out, " T, centering force ")?;
        write_real(out, work.class_centering[class] * 1.0e-6)?;
        writeln!(out, " MN")?;
    }
    writeln!(out, "  peak-field fit (radial width m: computed T / fitted T):")?;
    for index in 0..outcome.fit.count {
        write!(out, "    ")?;
        write_real(out, outcome.fit.widths[index])?;
        write!(out, ": ")?;
        write_real(out, outcome.fit.peaks[index])?;
        write!(out, " / ")?;
        write_real(out, outcome.fit.fitted[index])?;
        writeln!(out)?;
    }
    write!(out, "    largest relative deviation of the fit: ")?;
    write_real(out, outcome.fit.max_relative_error)?;
    writeln!(out)?;
    Ok(())
}

/// Plain-text report: every value, its unit and how it was obtained.
///
/// `inputs` holds the paths of the equilibrium and the coil file and the
/// configuration name.
///
/// # Errors
/// Returns an error if the outcome is incomplete or writing fails.
pub fn write_text_report<W: Write>(
    out: &mut W,
    inputs: (&str, &str, &str),
    settings: &Settings,
    outcome: &Outcome,
    work: &Workspace,
    runtime: f64,
) -> Result<(), Error> {
    ensure(runtime.is_finite() && runtime >= 0.0, BAD, "the run time is not valid")?;
    writeln!(out, "{GENERATOR}  (schema {SCHEMA_VERSION})")?;
    writeln!(out, "equilibrium: {}\ncoils:       {}", inputs.0, inputs.1)?;
    text_consistency(out, settings, outcome)?;
    text_values(out, inputs.2, outcome)?;
    text_diagnostics(out, outcome, work)?;
    let mut heading = true;
    for (index, active) in outcome.warnings.iter().enumerate() {
        if *active {
            if heading {
                writeln!(out, "\nWarnings")?;
                heading = false;
            }
            write!(out, "  - ")?;
            warning_text(out, settings, outcome, index)?;
            writeln!(out)?;
        }
    }
    writeln!(out, "\nruntime: {runtime:.2} s")?;
    Ok(())
}

/// Write `"name": value` with the given indentation and a trailing comma.
fn json_real<W: Write>(out: &mut W, indent: &str, name: &str, value: f64) -> Result<(), Error> {
    ensure(!name.is_empty(), BAD, "a JSON entry needs a name")?;
    write!(out, "{indent}\"{name}\": ")?;
    write_json_real(out, value)?;
    writeln!(out, ",")?;
    Ok(())
}

/// Write `"name": [values]` with the given indentation, without a comma.
fn json_array<W: Write>(
    out: &mut W,
    indent: &str,
    name: &str,
    values: &[f64],
) -> Result<(), Error> {
    ensure(!name.is_empty(), BAD, "a JSON entry needs a name")?;
    write!(out, "{indent}\"{name}\": [")?;
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            write!(out, ", ")?;
        }
        write_json_real(out, *value)?;
    }
    write!(out, "]")?;
    Ok(())
}

fn json_settings<W: Write>(out: &mut W, name: &str, settings: &Settings) -> Result<(), Error> {
    ensure(settings.n_scan <= settings.scan_factors.len(), BAD, "too many scan factors")?;
    write!(out, "  \"settings\": {{\n    \"name\": ")?;
    write_json_text(out, name)?;
    writeln!(out, ",")?;
    json_real(out, "    ", "wp_radial", settings.wp_radial)?;
    json_real(out, "    ", "wp_toroidal", settings.wp_toroidal)?;
    json_real(out, "    ", "epseff", settings.epseff)?;
    json_real(out, "    ", "neutron_peakfactor", settings.neutron_peakfactor)?;
    json_real(out, "    ", "rmajor_ref", settings.rmajor_ref.unwrap_or(f64::NAN))?;
    json_real(out, "    ", "bt_ref", settings.bt_ref.unwrap_or(f64::NAN))?;
    writeln!(out, "    \"match_field\": {},", settings.match_field)?;
    writeln!(out, "    \"allow_inconsistent\": {},", settings.allow_inconsistent)?;
    writeln!(out, "    \"n_points\": {},", settings.n_points)?;
    writeln!(
        out,
        "    \"wp_filaments\": [{}, {}],",
        settings.wp_filaments.0, settings.wp_filaments.1
    )?;
    json_array(out, "    ", "wp_scan_factors", &settings.scan_factors[..settings.n_scan])?;
    writeln!(out, ",\n    \"surface_resolution\": {},", settings.surface_resolution)?;
    json_real(out, "    ", "tolerance_linkage", settings.tolerances.linkage)?;
    json_real(out, "    ", "tolerance_axis_misalignment", settings.tolerances.axis_misalignment)?;
    write!(out, "    \"tolerance_field\": ")?;
    write_json_real(out, settings.tolerances.field)?;
    writeln!(out, "\n  }},")?;
    Ok(())
}

fn json_values<W: Write>(out: &mut W, name: &str, outcome: &Outcome) -> Result<(), Error> {
    ensure(!name.is_empty(), BAD, "the configuration needs a name")?;
    writeln!(out, "  \"values\": [")?;
    for (index, entry) in KEYS.iter().enumerate() {
        write!(out, "    {{\"key\": \"{}\", \"value\": ", entry.name)?;
        if entry.kind == Kind::Text {
            write_json_text(out, name)?;
        } else if entry.kind == Kind::Count {
            write!(out, "{:.0}", outcome.config.values[index])?;
        } else {
            write_json_real(out, outcome.config.values[index])?;
        }
        write!(
            out,
            ", \"unit\": \"{}\", \"description\": \"{}\", ",
            entry.unit, entry.description
        )?;
        writeln!(
            out,
            "\"method\": \"{}\"}}{}",
            entry.method,
            if index + 1 < KEYS.len() { "," } else { "" }
        )?;
    }
    writeln!(out, "  ],")?;
    Ok(())
}

fn json_consistency<W: Write>(
    out: &mut W,
    settings: &Settings,
    outcome: &Outcome,
) -> Result<(), Error> {
    let check = outcome.consistency;
    ensure(check.n_coils >= 1, BAD, "the consistency result is empty")?;
    writeln!(out, "  \"consistency\": {{")?;
    writeln!(out, "    \"nfp_equilibrium\": {},", check.nfp_equilibrium)?;
    writeln!(out, "    \"nfp_coil_file\": {},", check.nfp_coil_file)?;
    writeln!(out, "    \"n_coils\": {},", check.n_coils)?;
    json_real(out, "    ", "linkage_ratio", check.linkage_ratio)?;
    json_real(out, "    ", "axis_misalignment", check.axis_misalignment)?;
    json_real(out, "    ", "b_axis_coils_T", check.b_axis_coils)?;
    json_real(out, "    ", "b_axis_equilibrium_T", check.b_axis_equilibrium)?;
    json_real(out, "    ", "field_ratio_coils_over_equilibrium", check.field_ratio)?;
    write!(out, "    \"passed\": {},\n    \"failures\": [", check.passed)?;
    let mut first = true;
    for (index, failed) in check.failed.iter().enumerate() {
        if *failed {
            write!(out, "{}\"", if first { "" } else { ", " })?;
            failure_text(out, &check, &settings.tolerances, index)?;
            write!(out, "\"")?;
            first = false;
        }
    }
    writeln!(out, "]\n  }},\n  \"scaling\": {{")?;
    json_real(out, "    ", "length_factor", outcome.scaling.length_factor)?;
    json_real(out, "    ", "field_factor", outcome.scaling.field_factor)?;
    json_real(out, "    ", "current_factor_to_match_equilibrium", outcome.scaling.current_factor)?;
    write!(out, "    \"total_current_factor\": ")?;
    write_json_real(out, outcome.scaling.total_current_factor)?;
    writeln!(out, "\n  }},")?;
    Ok(())
}

fn json_classes<W: Write>(out: &mut W, work: &Workspace) -> Result<(), Error> {
    let count = work.classes.n_classes;
    ensure(
        count >= 1 && count <= work.class_centering.len(),
        BAD,
        "the class buffers are damaged",
    )?;
    writeln!(out, "    \"coil_classes\": [")?;
    for class in 0..count {
        let coil = work.classes.representatives[class];
        write!(
            out,
            "      {{\"representative_coil\": {coil}, \"count\": {}, \"current_MA\": ",
            work.classes.multiplicity[class]
        )?;
        write_json_real(out, work.coils.currents[coil] * 1.0e-6)?;
        write!(out, ", \"length_m\": ")?;
        write_json_real(out, work.coils.lengths[coil])?;
        write!(out, ", \"self_inductance_H\": ")?;
        write_json_real(out, work.class_self_inductance[class])?;
        write!(out, ", \"peak_field_T\": ")?;
        write_json_real(out, class_peak(work, class)?)?;
        write!(out, ", \"centering_force_MN\": ")?;
        write_json_real(out, work.class_centering[class] * 1.0e-6)?;
        writeln!(out, "}}{}", if class + 1 < count { "," } else { "" })?;
    }
    writeln!(out, "    ],")?;
    Ok(())
}

fn json_reported<W: Write>(out: &mut W, outcome: &Outcome) -> Result<(), Error> {
    let reported = outcome.reported;
    ensure(outcome.plasma.volume > 0.0, BAD, "the plasma geometry is missing")?;
    writeln!(out, "    \"plasma_as_reported_by_equilibrium\": {{")?;
    json_real(out, "      ", "Rmajor_p", reported.rmajor_p)?;
    json_real(out, "      ", "Aminor_p", reported.aminor_p)?;
    json_real(out, "      ", "aspect", reported.aspect)?;
    json_real(out, "      ", "volume_p", reported.volume_p)?;
    json_real(out, "      ", "b0", reported.b0)?;
    json_real(out, "      ", "volavgB", reported.volavg_b)?;
    json_real(out, "      ", "betatotal", reported.betatotal)?;
    json_real(out, "      ", "rbtor", reported.rbtor)?;
    write!(out, "      \"ier_flag\": ")?;
    write_json_real(out, reported.ier_flag)?;
    writeln!(
        out,
        "
    }},"
    )?;
    Ok(())
}

fn json_diagnostics<W: Write>(
    out: &mut W,
    outcome: &Outcome,
    work: &Workspace,
) -> Result<(), Error> {
    let fit = outcome.fit;
    ensure(fit.count <= fit.widths.len(), BAD, "the peak-field fit is damaged")?;
    ensure(outcome.geometry.fit_radii.len() == FIT_SURFACES, BAD, "the distance fit is damaged")?;
    writeln!(out, "  \"diagnostics\": {{")?;
    json_reported(out, outcome)?;
    json_real(out, "    ", "field_on_axis_mean_toroidal_component_T", outcome.b_axis_toroidal)?;
    json_real(out, "    ", "magnetic_energy_GJ", outcome.energy.energy * 1.0e-9)?;
    json_classes(out, work)?;
    writeln!(out, "    \"peak_field_fit\": {{")?;
    json_array(out, "      ", "radial_width_m", &fit.widths[..fit.count])?;
    writeln!(out, ",")?;
    json_array(out, "      ", "peak_field_T", &fit.peaks[..fit.count])?;
    writeln!(out, ",")?;
    json_array(out, "      ", "fitted_T", &fit.fitted[..fit.count])?;
    write!(out, ",\n      \"max_relative_error\": ")?;
    write_json_real(out, fit.max_relative_error)?;
    writeln!(out, "\n    }},\n    \"plasma_coil_distance_fit\": {{")?;
    json_array(out, "      ", "minor_radius_m", &outcome.geometry.fit_radii)?;
    writeln!(out, ",")?;
    json_array(out, "      ", "distance_m", &outcome.geometry.fit_distances)?;
    writeln!(out, "\n    }},")?;
    json_real(out, "    ", "winding_pack_clearance_m", outcome.clearance)?;
    write!(out, "    \"mean_coil_spacing_at_coil_rmajor_m\": ")?;
    let spacing = core::f64::consts::TAU * outcome.geometry.measures.r_major;
    write_json_real(out, spacing / f64::from(u32::try_from(work.coils.n_coils)?))?;
    writeln!(out, "\n  }},")?;
    Ok(())
}

/// JSON report: settings, every value with its method, the consistency
/// check, diagnostics and warnings.
///
/// `inputs` holds the paths of the equilibrium and the coil file and the
/// configuration name.
///
/// # Errors
/// Returns an error if the outcome is incomplete or writing fails.
pub fn write_json_report<W: Write>(
    out: &mut W,
    inputs: (&str, &str, &str),
    settings: &Settings,
    outcome: &Outcome,
    work: &Workspace,
    runtime: f64,
) -> Result<(), Error> {
    ensure(runtime.is_finite() && runtime >= 0.0, BAD, "the run time is not valid")?;
    writeln!(
        out,
        "{{\n  \"generator\": \"{GENERATOR}\",\n  \"schema_version\": \"{SCHEMA_VERSION}\","
    )?;
    write!(out, "  \"inputs\": {{\"equilibrium\": ")?;
    write_json_text(out, inputs.0)?;
    write!(out, ", \"coils\": ")?;
    write_json_text(out, inputs.1)?;
    writeln!(out, "}},")?;
    json_settings(out, inputs.2, settings)?;
    json_values(out, inputs.2, outcome)?;
    json_consistency(out, settings, outcome)?;
    json_diagnostics(out, outcome, work)?;
    write!(out, "  \"warnings\": [")?;
    let mut first = true;
    for (index, active) in outcome.warnings.iter().enumerate() {
        if *active {
            write!(out, "{}\"", if first { "" } else { ", " })?;
            warning_text(out, settings, outcome, index)?;
            write!(out, "\"")?;
            first = false;
        }
    }
    writeln!(out, "],\n  \"runtime_s\": {runtime:.2}\n}}")?;
    Ok(())
}
