//! Builders for synthetic inputs with known analytic properties.

use core::f64::consts::TAU;
use core::fmt::Write as _;
use stellarator_configprocessor::error::Error;
use stellarator_configprocessor::num::{floor_count, real};

/// One variable of a synthetic `NetCDF` classic file.
struct Entry {
    name: &'static str,
    /// Dimension identifiers; empty for a scalar.
    dims: Vec<u32>,
    /// Whether the values are stored as 32-bit integers instead of doubles.
    integer: bool,
    values: Vec<f64>,
}

fn put_u32(out: &mut Vec<u8>, value: usize) -> Result<(), Error> {
    out.extend_from_slice(&u32::try_from(value)?.to_be_bytes());
    Ok(())
}

fn put_name(out: &mut Vec<u8>, name: &str) -> Result<(), Error> {
    put_u32(out, name.len())?;
    out.extend_from_slice(name.as_bytes());
    for _ in name.len()..name.len().div_ceil(4) * 4 {
        out.push(0);
    }
    Ok(())
}

/// Bytes of a `NetCDF` classic (CDF-1) file with the given dimensions and
/// variables.
fn netcdf_bytes(dims: &[(&str, usize)], entries: &[Entry]) -> Result<Vec<u8>, Error> {
    let mut header_size = 4 + 4 + 8 + 8 + 8;
    for (name, _) in dims {
        header_size += 4 + name.len().div_ceil(4) * 4 + 4;
    }
    for entry in entries {
        header_size += 4 + entry.name.len().div_ceil(4) * 4 + 4 + 4 * entry.dims.len() + 8 + 12;
    }
    let mut out = Vec::new();
    out.extend_from_slice(b"CDF\x01");
    put_u32(&mut out, 0)?;
    put_u32(&mut out, 10)?;
    put_u32(&mut out, dims.len())?;
    for (name, length) in dims {
        put_name(&mut out, name)?;
        put_u32(&mut out, *length)?;
    }
    put_u32(&mut out, 0)?;
    put_u32(&mut out, 0)?;
    put_u32(&mut out, 11)?;
    put_u32(&mut out, entries.len())?;
    let mut begin = header_size;
    for entry in entries {
        let width = if entry.integer { 4 } else { 8 };
        put_name(&mut out, entry.name)?;
        put_u32(&mut out, entry.dims.len())?;
        for id in &entry.dims {
            out.extend_from_slice(&id.to_be_bytes());
        }
        put_u32(&mut out, 0)?;
        put_u32(&mut out, 0)?;
        put_u32(&mut out, if entry.integer { 4 } else { 6 })?;
        put_u32(&mut out, entry.values.len() * width)?;
        put_u32(&mut out, begin)?;
        begin += entry.values.len() * width;
    }
    for entry in entries {
        for value in &entry.values {
            if entry.integer {
                let whole = i32::try_from(floor_count(*value, 1000)?)?;
                out.extend_from_slice(&whole.to_be_bytes());
            } else {
                out.extend_from_slice(&value.to_be_bytes());
            }
        }
    }
    Ok(out)
}

fn scalar(name: &'static str, value: f64, integer: bool) -> Entry {
    Entry { name, dims: Vec::new(), integer, values: vec![value] }
}

/// A `wout` file for a circular torus: `R = r0 + a*sqrt(s)*cos(theta)`.
pub(crate) fn torus_wout(
    r0: f64,
    a: f64,
    nfp: usize,
    b0: f64,
    ns: usize,
) -> Result<Vec<u8>, Error> {
    let (mut r_cos, mut z_sin, mut b_cos) = (Vec::new(), Vec::new(), Vec::new());
    for index in 0..ns {
        let s = real(index)? / real(ns - 1)?;
        r_cos.extend_from_slice(&[r0, a * s.sqrt()]);
        z_sin.extend_from_slice(&[0.0, a * s.sqrt()]);
        b_cos.push(if index == 0 { 0.0 } else { b0 });
    }
    let entries = [
        scalar("nfp", real(nfp)?, true),
        scalar("lasym__logical__", 0.0, true),
        scalar("Rmajor_p", r0, false),
        scalar("Aminor_p", a, false),
        scalar("volume_p", TAU * core::f64::consts::PI * r0 * a * a, false),
        scalar("b0", b0, false),
        Entry { name: "xm", dims: vec![1], integer: false, values: vec![0.0, 1.0] },
        Entry { name: "xn", dims: vec![1], integer: false, values: vec![0.0, 0.0] },
        Entry { name: "xm_nyq", dims: vec![2], integer: false, values: vec![0.0] },
        Entry { name: "xn_nyq", dims: vec![2], integer: false, values: vec![0.0] },
        Entry { name: "rmnc", dims: vec![0, 1], integer: false, values: r_cos },
        Entry { name: "zmns", dims: vec![0, 1], integer: false, values: z_sin },
        Entry { name: "bmnc", dims: vec![0, 2], integer: false, values: b_cos },
        Entry { name: "iotaf", dims: vec![0], integer: false, values: vec![0.9; ns] },
    ];
    netcdf_bytes(&[("radius", ns), ("mn_mode", 2), ("mn_mode_nyq", 1)], &entries)
}

/// MAKEGRID text for planar circular coils of radius `a_coil` centred on the
/// circle `R = r0`, evenly spaced in toroidal angle.
pub(crate) fn circular_coils_text(
    n_coils: usize,
    r0: f64,
    a_coil: f64,
    current: f64,
    n_points: usize,
    periods: usize,
) -> Result<String, Error> {
    let mut text = String::new();
    writeln!(text, "periods {periods}\nbegin filament\nmirror NIL")?;
    for coil in 0..n_coils {
        let phi = TAU * (real(coil)? + 0.5) / real(n_coils)?;
        for index in 0..=n_points {
            let theta = TAU * real(index % n_points)? / real(n_points)?;
            let radius = r0 + a_coil * theta.cos();
            let (x, y, z) = (radius * phi.cos(), radius * phi.sin(), a_coil * theta.sin());
            if index < n_points {
                writeln!(text, "{x:.15E} {y:.15E} {z:.15E} {current:.15E}")?;
            } else {
                writeln!(text, "{x:.15E} {y:.15E} {z:.15E} 0.0 1 TF")?;
            }
        }
    }
    writeln!(text, "end")?;
    Ok(text)
}
