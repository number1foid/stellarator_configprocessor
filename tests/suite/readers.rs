//! Input readers, sampling, plasma geometry and the consistency check.

use crate::helpers::{close, near, reference_case, NO_FACTORS};
use crate::synthetic::{circular_coils_text, torus_wout};
use core::f64::consts::{PI, TAU};
use stellarator_configprocessor::coilgeom::coil_measures;
use stellarator_configprocessor::coils::{new_sampled_coils, sample_coils};
use stellarator_configprocessor::consistency::{
    axis_field, check_consistency, new_axis_field, DEFAULT_TOLERANCES,
};
use stellarator_configprocessor::curve::new_spline_work;
use stellarator_configprocessor::equilibrium::{new_coefficients, read_wout};
use stellarator_configprocessor::error::{Code, Error};
use stellarator_configprocessor::makegrid::read_makegrid;
use stellarator_configprocessor::plasma::surface_geometry;

#[test]
fn torus_geometry_is_exact() -> Result<(), Error> {
    let eq = read_wout(&torus_wout(10.0, 2.0, 5, 5.0, 11)?)?;
    assert_eq!((eq.nfp, eq.n_radial, eq.n_modes, eq.lasym), (5, 11, 2, false));
    let mut surface = new_coefficients(eq.n_modes)?;
    let geometry = surface_geometry(&eq, 1.0, &mut surface)?;
    close(geometry.volume, TAU * PI * 10.0 * 4.0, 1.0e-12, "volume");
    close(geometry.area, 4.0 * PI * PI * 10.0 * 2.0, 1.0e-12, "area");
    close(geometry.minor_radius, 2.0, 1.0e-12, "minor radius");
    close(geometry.major_radius, 10.0, 1.0e-12, "major radius");
    let inner = surface_geometry(&eq, 0.5, &mut surface)?;
    close(inner.minor_radius, 2.0 * 0.5_f64.sqrt(), 1.0e-12, "minor radius at s = 0.5");
    Ok(())
}

#[test]
fn wout_reader_rejects_other_files() -> Result<(), Error> {
    let Err(error) = read_wout(b"\x89HDF\r\n\x1a\n and more") else {
        return Err(Error::new(Code::Numeric, "an HDF5 file was accepted"));
    };
    assert_eq!(error.code, Code::Equilibrium);
    assert!(error.message.contains("NetCDF-4"));
    let mut damaged = torus_wout(10.0, 2.0, 5, 5.0, 11)?;
    damaged.truncate(200);
    assert!(read_wout(&damaged).is_err());
    assert!(read_wout(b"not a file").is_err());
    Ok(())
}

#[test]
fn circular_coils_are_read_and_sampled() -> Result<(), Error> {
    let text = circular_coils_text(20, 10.0, 4.0, 1.0e6, 120, 5)?;
    let set = read_makegrid(&text, &NO_FACTORS)?;
    assert_eq!((set.periods, set.currents.len(), set.points.len()), (5, 20, 20 * 120));
    let mut work = new_spline_work(120, 160)?;
    let mut coils = new_sampled_coils(20, 160)?;
    sample_coils(&set, &mut work, &mut coils)?;
    for coil in 0..20 {
        close(coils.lengths[coil], TAU * 4.0, 1.0e-6, "coil length");
    }
    for curvature in &coils.curvature {
        close(*curvature, 0.25, 1.0e-3, "curvature");
    }
    let measures = coil_measures(&coils)?;
    close(measures.r_major, 10.0, 1.0e-6, "coil major radius");
    close(measures.r_minor, 4.0, 1.0e-6, "coil minor radius");
    close(measures.height, 8.0, 1.0e-5, "coil height");
    close(measures.surface_area, TAU * 10.0 * TAU * 4.0, 1.0e-6, "coil surface");
    Ok(())
}

#[test]
fn w7x_matches_the_python_reference() -> Result<(), Error> {
    let Some((eq, set)) = reference_case("wout_w7x_standard.nc", "coils.w7x_standard")? else {
        return Ok(());
    };
    assert_eq!((eq.nfp, set.periods, set.currents.len()), (5, 5, 50));
    let mut surface = new_coefficients(eq.n_modes)?;
    let geometry = surface_geometry(&eq, 1.0, &mut surface)?;
    close(geometry.major_radius, 5.511_939_539_133_372_5, 1.0e-10, "rmajor");
    close(geometry.minor_radius, 0.512_692_195_535_189_6, 1.0e-10, "rminor");
    close(geometry.volume, 28.598_786_068_665_905, 1.0e-10, "volume");
    close(geometry.area, 136.662_192_597_402_22, 1.0e-10, "area");
    close(geometry.volume, eq.reported.volume_p, 1.0e-9, "volume against VMEC");
    let mut work = new_spline_work(200, 160)?;
    let mut coils = new_sampled_coils(50, 160)?;
    sample_coils(&set, &mut work, &mut coils)?;
    close(coils.lengths[0], 8.557_318_034_722_847, 1.0e-9, "length of coil 0");
    let mut axis = new_axis_field(50)?;
    axis_field(&eq, &coils, &mut surface, &mut axis)?;
    let check = check_consistency(&eq, &coils, set.periods, &DEFAULT_TOLERANCES, &axis, true)?;
    assert!(check.passed);
    close(check.linkage_ratio, 0.999_957_060_969_174_8, 1.0e-9, "linkage");
    close(check.b_axis_coils, 2.784_718_094_363_533_7, 1.0e-9, "field on axis");
    close(check.b_axis_equilibrium, 2.784_731_050_149_441_8, 1.0e-10, "equilibrium field");
    near(check.axis_misalignment, 0.000_252_230_699_990_212_27, 1.0e-9, "misalignment");
    Ok(())
}
