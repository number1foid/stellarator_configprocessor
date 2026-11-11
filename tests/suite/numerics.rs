//! Field, winding pack, energy and force building blocks against analytic
//! results.

use crate::helpers::{close, failure, near};
use core::f64::consts::{PI, TAU};
use stellarator_configprocessor::curve::{new_spline_work, resample};
use stellarator_configprocessor::energy::{gmd_squared, new_inductance_work, self_inductance};
use stellarator_configprocessor::error::{Code, Error};
use stellarator_configprocessor::field::add_polygon_field;
use stellarator_configprocessor::num::{ceil_count, floor_count, real, Vec3, MU0};
use stellarator_configprocessor::report::write_real;
use stellarator_configprocessor::winding_pack::{
    bar_kernel, fit_peak_field, pack_layout, winding_pack,
};

/// A circle of the given radius in the plane `z = 0`.
fn circle(radius: f64, count: usize) -> Result<Vec<Vec3>, Error> {
    let mut points = Vec::new();
    for index in 0..count {
        let angle = TAU * real(index)? / real(count)?;
        points.push([radius * angle.cos(), radius * angle.sin(), 0.0]);
    }
    Ok(points)
}

#[test]
fn circular_loop_on_axis() -> Result<(), Error> {
    let (radius, current) = (1.5, 2.0e5);
    let loop_points = circle(radius, 400)?;
    for height in [0.0, 0.4, 2.0] {
        let mut field = [0.0; 3];
        add_polygon_field([0.0, 0.0, height], &loop_points, current, &mut field)?;
        let exact =
            MU0 * current * radius * radius / (2.0 * (radius * radius + height * height).powf(1.5));
        close(field[2], exact, 1.0e-4, "axial field of a loop");
        near(field[0], 0.0, 1.0e-12, "transverse field of a loop");
        near(field[1], 0.0, 1.0e-12, "transverse field of a loop");
    }
    Ok(())
}

#[test]
fn long_straight_wire() -> Result<(), Error> {
    // A very long, thin rectangle: near one long side the other three sides
    // contribute next to nothing.
    let (current, distance, half) = (1.0e4, 0.3, 1.0e4);
    let wire = [[0.0, 0.0, -half], [0.0, 0.0, half], [1.0e7, 0.0, half], [1.0e7, 0.0, -half]];
    let mut field = [0.0; 3];
    add_polygon_field([distance, 0.0, 0.0], &wire, current, &mut field)?;
    // The two short sides add the field of two half-infinite wires at
    // distance `half`.
    let ideal = MU0 * current / (TAU * distance);
    close(field[1], ideal * (1.0 + distance / half), 1.0e-6, "field of a straight wire");
    near(field[0], 0.0, 1.0e-12, "field along the offset");
    Ok(())
}

#[test]
fn field_point_on_the_conductor_is_skipped() -> Result<(), Error> {
    let square = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
    let mut field = [0.0; 3];
    add_polygon_field([0.5, 0.0, 0.0], &square, 1.0, &mut field)?;
    assert!(field[2].is_finite() && field[2] > 0.0, "unexpected field {field:?}");
    let error = failure(add_polygon_field([f64::NAN, 0.0, 0.0], &square, 1.0, &mut field), "nan");
    assert_eq!(error.code, Code::Numeric);
    assert!(add_polygon_field([0.0; 3], &square[..1], 1.0, &mut field).is_err());
    Ok(())
}

#[test]
fn spline_through_a_circle() -> Result<(), Error> {
    let points = circle(2.0, 60)?;
    let mut work = new_spline_work(60, 128)?;
    let (mut positions, mut tangents, mut curvature) =
        (vec![[0.0; 3]; 128], vec![[0.0; 3]; 128], vec![0.0; 128]);
    let length = resample(&points, &mut work, &mut positions, &mut tangents, &mut curvature)?;
    close(length, TAU * 2.0, 1.0e-5, "circumference");
    for (position, value) in positions.iter().zip(curvature.iter()) {
        close(position[0].hypot(position[1]), 2.0, 1.0e-5, "radius of a sample");
        close(*value, 0.5, 2.0e-3, "curvature");
    }
    let first_gap = (positions[1][0] - positions[0][0]).hypot(positions[1][1] - positions[0][1]);
    close(first_gap, length / 128.0, 1.0e-3, "sample spacing");
    assert!(
        resample(&points[..2], &mut work, &mut positions, &mut tangents, &mut curvature).is_err()
    );
    Ok(())
}

#[test]
fn geometric_mean_distance_of_a_square() -> Result<(), Error> {
    let square = winding_pack(0.2, 0.2)?;
    close(gmd_squared(&square)?.sqrt(), 0.447_05 * 0.2, 1.0e-4, "GMD of a square");
    let swapped = (gmd_squared(&winding_pack(0.3, 0.1)?)?, gmd_squared(&winding_pack(0.1, 0.3)?)?);
    close(swapped.0, swapped.1, 1.0e-12, "GMD of a rotated rectangle");
    assert!(winding_pack(0.0, 0.1).is_err());
    assert!(winding_pack(0.1, f64::NAN).is_err());
    Ok(())
}

#[test]
fn self_inductance_of_a_thin_ring() -> Result<(), Error> {
    // L = mu0 R (ln(8 R / GMD) - 2) for a ring that is thin compared with
    // its radius.
    let (radius, side) = (5.0, 0.05);
    let wp = winding_pack(side, side)?;
    let ring = circle(radius, 200)?;
    let mut spline = new_spline_work(200, 16_384)?;
    let mut work = new_inductance_work()?;
    let computed = self_inductance(&ring, &wp, &mut spline, &mut work)?;
    let gmd = gmd_squared(&wp)?.sqrt();
    let exact = MU0 * radius * ((8.0 * radius / gmd).ln() - 2.0);
    close(computed, exact, 2.0e-3, "self-inductance of a ring");
    Ok(())
}

#[test]
fn bar_kernel_limits() -> Result<(), Error> {
    let (f_u, f_v) = bar_kernel(0.0, 0.0, 0.2, 0.1)?;
    near(f_u, 0.0, 1.0e-14, "radial kernel at the centre");
    near(f_v, 0.0, 1.0e-14, "toroidal kernel at the centre");
    // Far away the bar acts like a filament carrying its whole area.
    let (far_u, far_v) = bar_kernel(30.0, 40.0, 0.2, 0.1)?;
    close(far_u, 0.02 * 30.0 / 2500.0, 1.0e-5, "radial kernel far away");
    close(far_v, 0.02 * 40.0 / 2500.0, 1.0e-5, "toroidal kernel far away");
    assert!(bar_kernel(0.0, 0.0, 0.0, 0.1).is_err());
    Ok(())
}

#[test]
fn pack_layout_weights_sum_to_one() -> Result<(), Error> {
    let layout = pack_layout(&winding_pack(0.3, 0.2)?, 4, 3)?;
    let mut total = 0.0;
    for row in &layout.weight {
        for weight in row {
            total += weight;
        }
    }
    close(total, 1.0, 1.0e-12, "sum of the quadrature weights");
    close(layout.corner_u[0], -0.15, 1.0e-12, "first radial corner");
    close(layout.corner_u[4], 0.15, 1.0e-12, "last radial corner");
    close(layout.centre_v[0], -0.1 + 0.1 / 3.0, 1.0e-12, "first toroidal centre");
    assert!(pack_layout(&winding_pack(0.3, 0.2)?, 0, 3).is_err());
    assert!(pack_layout(&winding_pack(0.3, 0.2)?, 4, 13).is_err());
    Ok(())
}

#[test]
fn peak_field_fit_recovers_its_parameters() -> Result<(), Error> {
    let (a1, a2, r_major, r_minor, current_ma, n_coils) = (0.65, 0.03, 10.0, 4.0, 2.0, 20);
    let factor = 0.2 * current_ma * 20.0 / (r_major - r_minor);
    let mut widths = [0.0; 8];
    let mut peaks = [0.0; 8];
    for (index, width) in [0.2, 0.3, 0.4, 0.5].iter().enumerate() {
        widths[index] = *width;
        peaks[index] = factor * (a1 + a2 * r_major / width);
    }
    let fit = fit_peak_field((&widths, &peaks), 4, current_ma, n_coils, (r_major, r_minor))?;
    close(fit.a1, a1, 1.0e-10, "a1");
    close(fit.a2, a2, 1.0e-10, "a2");
    near(fit.max_relative_error, 0.0, 1.0e-12, "fit error");
    assert!(fit_peak_field((&widths, &peaks), 1, current_ma, n_coils, (r_major, r_minor)).is_err());
    assert!(fit_peak_field((&widths, &peaks), 4, current_ma, n_coils, (4.0, 10.0)).is_err());
    Ok(())
}

#[test]
fn counts_from_real_numbers() -> Result<(), Error> {
    assert_eq!(floor_count(0.0, 100)?, 0);
    assert_eq!(floor_count(7.999, 100)?, 7);
    assert_eq!(floor_count(8.0, 100)?, 8);
    assert_eq!(floor_count(1.0e9, 100)?, 100);
    assert_eq!(ceil_count(7.001, 100)?, 8);
    assert_eq!(ceil_count(7.0, 100)?, 7);
    assert!(floor_count(-1.0, 100).is_err());
    assert!(floor_count(f64::NAN, 100).is_err());
    assert!(ceil_count(100.5, 100).is_err());
    close(real(123_456)?, 123_456.0, 0.0, "count as a real number");
    assert!(real(usize::MAX).is_err());
    Ok(())
}

#[test]
fn six_significant_digits() -> Result<(), Error> {
    let cases = [
        (PI, "3.14159"),
        (-0.012_345_678, "-0.0123457"),
        (123_456.7, "123457"),
        (0.0, "0"),
        (1.5e-7, "1.50000e-7"),
        (2.5e9, "2.50000e9"),
    ];
    for (number, expected) in cases {
        let mut text = Vec::new();
        write_real(&mut text, number)?;
        assert_eq!(String::from_utf8_lossy(&text), expected);
    }
    Ok(())
}
