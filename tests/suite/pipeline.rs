//! The whole calculation: an analytic torus and the reference machines.

use crate::helpers::{
    close, compute, failure, ideal_toroid_field, near, reference_case, torus_case, value,
};
use crate::reference::{
    HSX, HSX_ENERGY, HSX_LOFTED, W7X, W7X_ENERGY, W7X_LOFTED, W7X_SCALED, W7X_SCALED_ENERGY,
    W7X_SCALED_LOFTED,
};
use core::f64::consts::TAU;
use stellarator_configprocessor::consistency::{TEST_FIELD, TEST_PERIODS};
use stellarator_configprocessor::error::{Code, Error};
use stellarator_configprocessor::makegrid::scale_coils;
use stellarator_configprocessor::num::MU0;
use stellarator_configprocessor::pipeline::{
    default_settings, validate_settings, Outcome, Settings, WARN_OVERLAP, WARN_OVERRIDDEN,
    WARN_SIZES_SKIPPED,
};
use stellarator_configprocessor::report::{write_json_report, write_text_report};

fn torus_settings() -> Settings {
    let mut settings = default_settings((0.5, 0.5), 0.01, 1.6);
    settings.n_points = 96;
    settings
}

/// Compare every configuration value with the Python reference.
fn check_reference(outcome: &Outcome, expected: &[(&str, f64)], lofted: f64, energy: f64) {
    for (name, reference) in expected {
        close(value(outcome, name), *reference, 1.0e-8, name);
    }
    close(outcome.geometry.lofted_surface, lofted, 1.0e-8, "lofted coil surface");
    close(outcome.energy.energy * 1.0e-9, energy, 1.0e-8, "magnetic energy");
    assert_eq!(outcome.warnings, [false; 10], "unexpected warnings");
}

#[test]
fn torus_with_circular_coils() -> Result<(), Error> {
    let (mut eq, mut set) = torus_case()?;
    let (outcome, work) = compute(&mut eq, &mut set, &torus_settings())?;
    assert!(outcome.consistency.passed);
    close(outcome.consistency.linkage_ratio, 1.0, 1.0e-3, "linked current");
    close(value(&outcome, "symmetry"), 5.0, 0.0, "symmetry");
    close(value(&outcome, "coilspermodule"), 4.0, 0.0, "coils per module");
    close(value(&outcome, "rmajor_ref"), 10.0, 1.0e-10, "plasma major radius");
    close(value(&outcome, "rminor_ref"), 2.0, 1.0e-10, "plasma minor radius");
    close(value(&outcome, "vol_plasma"), TAU * core::f64::consts::PI * 40.0, 1.0e-10, "volume");
    close(value(&outcome, "coil_rmajor"), 10.0, 1.0e-6, "coil major radius");
    close(value(&outcome, "coil_rminor"), 4.0, 1.0e-6, "coil minor radius");
    close(value(&outcome, "coillength"), 20.0 * TAU * 4.0, 1.0e-6, "coil length");
    close(value(&outcome, "maximal_coil_height"), 8.0, 1.0e-4, "coil height");
    close(value(&outcome, "min_bend_radius"), 4.0, 2.0e-3, "bend radius");
    close(value(&outcome, "min_plasma_coil_distance"), 2.0, 1.0e-3, "plasma-coil distance");
    close(value(&outcome, "derivative_min_LCFS_coils_dist"), -1.0, 1.0e-3, "distance derivative");
    close(value(&outcome, "i0"), 1.0, 1.0e-12, "current per coil");
    close(
        value(&outcome, "bt_ref"),
        ideal_toroid_field(20.0, 1.0e6, 10.0),
        1.0e-2,
        "field on axis",
    );
    close(value(&outcome, "WP_area"), 0.25, 1.0e-12, "winding pack area");
    close(value(&outcome, "WP_ratio"), 1.0, 1.0e-12, "winding pack ratio");
    // The innermost point of each coil is at R = 6 m, 1.885 m from its neighbour.
    close(value(&outcome, "dmin"), 2.0 * 6.0 * (TAU / 40.0).sin(), 1.0e-3, "coil distance");
    close(value(&outcome, "max_portsize_width"), 2.0 * 10.0 * (TAU / 40.0).sin(), 1.0e-6, "port");
    // Ideal toroidal solenoid: L = mu0 N^2 (R - sqrt(R^2 - a^2)).
    let ideal = MU0 * 400.0 * (10.0 - 84.0_f64.sqrt());
    // Twenty separate coils store more energy than the ideal current sheet.
    close(value(&outcome, "inductance"), ideal, 0.25, "inductance");
    assert!(value(&outcome, "inductance") > ideal, "less energy than the ideal solenoid");
    // With five field periods and the stellarator-symmetry flip, the coils
    // at 9 and 63 degrees form one class and those at 27 and 45 the other.
    assert_eq!(work.classes.n_classes, 2);
    assert_eq!(work.classes.multiplicity[..2], [10, 10]);
    assert!(value(&outcome, "WP_bmax") > value(&outcome, "bt_ref"));
    assert!(value(&outcome, "centering_force_min_MN") > 0.0, "the coils are pulled inwards");
    let (highest, lowest) =
        (value(&outcome, "centering_force_max_MN"), value(&outcome, "centering_force_min_MN"));
    close(highest, lowest, 1.0e-9, "equal centering force on all coils");
    assert!(outcome.fit.max_relative_error < 0.02, "poor peak-field fit");
    Ok(())
}

#[test]
fn rescaling_follows_the_scaling_laws() -> Result<(), Error> {
    let (mut eq, mut set) = torus_case()?;
    let (base, _) = compute(&mut eq, &mut set, &torus_settings())?;
    let (mut eq, mut set) = torus_case()?;
    let mut settings = torus_settings();
    settings.rmajor_ref = Some(20.0);
    settings.bt_ref = Some(2.0 * value(&base, "bt_ref"));
    settings.wp_radial = 1.0;
    settings.wp_toroidal = 1.0;
    let (scaled, _) = compute(&mut eq, &mut set, &settings)?;
    close(scaled.scaling.length_factor, 2.0, 1.0e-10, "length factor");
    close(scaled.scaling.total_current_factor, 4.0, 1.0e-10, "current factor");
    close(value(&scaled, "rmajor_ref"), 20.0, 1.0e-10, "major radius");
    close(value(&scaled, "bt_ref"), 2.0 * value(&base, "bt_ref"), 1.0e-9, "field");
    close(value(&scaled, "i0"), 4.0 * value(&base, "i0"), 1.0e-9, "current");
    close(value(&scaled, "coillength"), 2.0 * value(&base, "coillength"), 1.0e-9, "length");
    close(value(&scaled, "vol_plasma"), 8.0 * value(&base, "vol_plasma"), 1.0e-9, "volume");
    close(value(&scaled, "WP_bmax"), 2.0 * value(&base, "WP_bmax"), 1.0e-8, "peak field");
    close(value(&scaled, "inductance"), 2.0 * value(&base, "inductance"), 1.0e-8, "inductance");
    close(value(&scaled, "a1"), value(&base, "a1"), 1.0e-7, "a1 is dimensionless");
    close(value(&scaled, "a2"), value(&base, "a2"), 1.0e-7, "a2 is dimensionless");
    Ok(())
}

#[test]
fn inconsistent_inputs_stop_the_run() -> Result<(), Error> {
    let (mut eq, mut set) = torus_case()?;
    scale_coils(&mut set, 1.0, 1.2)?;
    let error = failure(compute(&mut eq, &mut set, &torus_settings()), "20 % more current");
    assert_eq!(error.code, Code::Inconsistent);

    let mut settings = torus_settings();
    settings.match_field = true;
    let (matched, _) = compute(&mut eq, &mut set, &settings)?;
    assert!(matched.consistency.passed, "matching the field removes the difference");
    close(matched.scaling.current_factor, 1.0 / 1.2, 1.0e-2, "current factor");

    let (mut eq, mut set) = torus_case()?;
    scale_coils(&mut set, 1.0, 1.2)?;
    set.periods = 4;
    settings = torus_settings();
    settings.allow_inconsistent = true;
    let (forced, _) = compute(&mut eq, &mut set, &settings)?;
    assert!(!forced.consistency.passed && forced.warnings[WARN_OVERRIDDEN]);
    assert!(forced.consistency.failed[TEST_FIELD] && forced.consistency.failed[TEST_PERIODS]);
    near(forced.consistency.field_ratio, 1.2, 2.0e-2, "field ratio");
    Ok(())
}

#[test]
fn oversized_winding_packs_are_reported() -> Result<(), Error> {
    let (mut eq, mut set) = torus_case()?;
    let mut settings = torus_settings();
    settings.wp_toroidal = 2.2;
    let (outcome, _) = compute(&mut eq, &mut set, &settings)?;
    assert!(outcome.clearance < 0.0 && outcome.warnings[WARN_OVERLAP]);
    assert!(outcome.warnings[WARN_SIZES_SKIPPED] && outcome.skipped[3] && !outcome.skipped[0]);
    Ok(())
}

#[test]
fn invalid_settings_are_rejected() {
    let good = torus_settings();
    assert!(validate_settings(&good).is_ok());
    let mut cases = [good; 7];
    cases[0].wp_radial = 0.0;
    cases[1].epseff = -0.1;
    cases[2].neutron_peakfactor = f64::NAN;
    cases[3].n_points = 8;
    cases[4].bt_ref = Some(-1.0);
    cases[5].scan_factors = [0.5, 0.6, 0.7, 0.85, 0.9, 0.0, 0.0, 0.0];
    cases[6].n_scan = 1;
    for settings in &cases {
        assert_eq!(failure(validate_settings(settings), "settings").code, Code::Settings);
    }
}

#[test]
fn reports_are_complete() -> Result<(), Error> {
    let (mut eq, mut set) = torus_case()?;
    let settings = torus_settings();
    let (outcome, work) = compute(&mut eq, &mut set, &settings)?;
    let inputs = ("wout \"torus\".nc", "coils.torus", "torus");
    let mut text = Vec::new();
    write_text_report(&mut text, inputs, &settings, &outcome, &work, 1.5)?;
    let text = String::from_utf8_lossy(&text);
    for part in ["result: passed", "coilspermodule", "peak-field fit", "runtime: 1.50 s"] {
        assert!(text.contains(part), "the text report lacks `{part}`");
    }
    let mut json = Vec::new();
    write_json_report(&mut json, inputs, &settings, &outcome, &work, 1.5)?;
    let json = String::from_utf8_lossy(&json);
    for part in ["\"wout \\\"torus\\\".nc\"", "\"passed\": true", "\"warnings\": []", "\"a2\""] {
        assert!(json.contains(part), "the JSON report lacks `{part}`");
    }
    let (opened, closed) = (json.matches('{').count(), json.matches('}').count());
    assert_eq!(opened, closed, "unbalanced braces in the JSON report");
    assert_eq!(json.matches('[').count(), json.matches(']').count());
    assert!(!json.contains("NaN") && !json.contains("inf"), "non-finite number in the report");
    Ok(())
}

#[test]
fn w7x_reproduces_the_python_results() -> Result<(), Error> {
    let Some((mut eq, mut set)) = reference_case("wout_w7x_standard.nc", "coils.w7x_standard")?
    else {
        return Ok(());
    };
    let mut settings = default_settings((0.18, 0.15), 0.0095, 1.6);
    settings.bt_ref = Some(3.0);
    let (outcome, _) = compute(&mut eq, &mut set, &settings)?;
    check_reference(&outcome, &W7X, W7X_LOFTED, W7X_ENERGY);
    close(value(&outcome, "bt_ref"), 3.0, 1.0e-12, "requested field");
    Ok(())
}

#[test]
fn w7x_rescaled_reproduces_the_python_results() -> Result<(), Error> {
    let Some((mut eq, mut set)) = reference_case("wout_w7x_standard.nc", "coils.w7x_standard")?
    else {
        return Ok(());
    };
    let mut settings = default_settings((0.6, 0.5), 0.0095, 1.6);
    settings.bt_ref = Some(5.5);
    settings.rmajor_ref = Some(20.0);
    let (outcome, _) = compute(&mut eq, &mut set, &settings)?;
    check_reference(&outcome, &W7X_SCALED, W7X_SCALED_LOFTED, W7X_SCALED_ENERGY);
    Ok(())
}

#[test]
fn hsx_reproduces_the_python_results() -> Result<(), Error> {
    let Some((mut eq, mut set)) = reference_case("wout_hsx_qhs.nc", "coils.hsx_qhs")? else {
        return Ok(());
    };
    assert_eq!((eq.nfp, set.currents.len()), (4, 48));
    let settings = default_settings((0.06, 0.05), 0.005, 1.6);
    let (outcome, work) = compute(&mut eq, &mut set, &settings)?;
    check_reference(&outcome, &HSX, HSX_LOFTED, HSX_ENERGY);
    assert_eq!(work.classes.n_classes, 6);
    Ok(())
}
