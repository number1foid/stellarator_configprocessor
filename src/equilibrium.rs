//! VMEC `wout` equilibrium: reading the file and evaluating flux surfaces.
//!
//! Conventions (VMEC):
//!
//! * `R(s, theta, phi) = sum r_cos * cos(m*theta - xn*phi)`
//! * `Z(s, theta, phi) = sum z_sin * sin(m*theta - xn*phi)`
//! * `xn` already contains the number of field periods, `xn = n * nfp`.
//! * `phi` is the cylindrical toroidal angle.
//! * Radial index 0 is the magnetic axis, the last index is the boundary.
//!   Geometry coefficients live on the full radial grid, which is uniform in
//!   normalised toroidal flux `s`.

use crate::error::{ensure, Code, Error};
use crate::limits::{MAX_MODES, MAX_PERIODS, MAX_RADIAL};
use crate::netcdf::{find_variable, parse_header, read_values, Header, Variable};
use crate::num::{all_finite, floor_count, real, Vec3};

const BAD: Code = Code::Equilibrium;

/// Scalars as reported in the file, kept for cross-checks and the report.
/// A value the file does not contain is not-a-number.
#[derive(Debug, Clone, Copy)]
pub struct Reported {
    /// Major radius [m].
    pub rmajor_p: f64,
    /// Minor radius [m].
    pub aminor_p: f64,
    /// Aspect ratio.
    pub aspect: f64,
    /// Plasma volume [m^3].
    pub volume_p: f64,
    /// Field on axis [T].
    pub b0: f64,
    /// Volume-averaged field [T].
    pub volavg_b: f64,
    /// Total beta.
    pub betatotal: f64,
    /// `R * B_toroidal` at the edge [T m].
    pub rbtor: f64,
    /// VMEC exit flag.
    pub ier_flag: f64,
}

/// Flux-surface geometry and axis field of one equilibrium.
#[derive(Debug, Clone)]
pub struct Equilibrium {
    /// Number of field periods.
    pub nfp: usize,
    /// Whether the equilibrium lacks stellarator symmetry.
    pub lasym: bool,
    /// Number of radial grid points.
    pub n_radial: usize,
    /// Number of Fourier modes of the geometry.
    pub n_modes: usize,
    /// Poloidal mode numbers, one per mode.
    pub xm: Vec<f64>,
    /// Toroidal mode numbers including `nfp`, one per mode.
    pub xn: Vec<f64>,
    /// `R` cosine coefficients, radial index times `n_modes` plus mode [m].
    pub r_cos: Vec<f64>,
    /// `Z` sine coefficients, same layout [m].
    pub z_sin: Vec<f64>,
    /// `R` sine coefficients; empty unless `lasym` [m].
    pub r_sin: Vec<f64>,
    /// `Z` cosine coefficients; empty unless `lasym` [m].
    pub z_cos: Vec<f64>,
    /// Toroidal mode numbers of the field strength on the axis.
    pub axis_xn: Vec<f64>,
    /// Cosine coefficients of the field strength on the axis [T].
    pub axis_cos: Vec<f64>,
    /// Sine coefficients of the field strength on the axis [T].
    pub axis_sin: Vec<f64>,
    /// Scalars as reported in the file.
    pub reported: Reported,
}

/// Fourier coefficients of one flux surface; a reusable buffer.
#[derive(Debug, Clone)]
pub struct Coefficients {
    /// `R` cosine coefficients [m].
    pub rc: Vec<f64>,
    /// `Z` sine coefficients [m].
    pub zs: Vec<f64>,
    /// `R` sine coefficients; unused unless `lasym` [m].
    pub rs: Vec<f64>,
    /// `Z` cosine coefficients; unused unless `lasym` [m].
    pub zc: Vec<f64>,
}

/// Cylindrical coordinates of a surface point and their angle derivatives.
#[derive(Debug, Clone, Copy, Default)]
pub struct SurfacePoint {
    /// Cylindrical radius [m].
    pub r: f64,
    /// Height [m].
    pub z: f64,
    /// `dR/dtheta`.
    pub r_t: f64,
    /// `dR/dphi`.
    pub r_p: f64,
    /// `dZ/dtheta`.
    pub z_t: f64,
    /// `dZ/dphi`.
    pub z_p: f64,
}

/// Field strength on the axis: mode numbers and coefficients.
#[derive(Debug)]
struct AxisModes {
    xn: Vec<f64>,
    cos: Vec<f64>,
    sin: Vec<f64>,
}

fn read_array(
    bytes: &[u8],
    header: &Header,
    name: &str,
    missing: &'static str,
) -> Result<(Vec<f64>, Variable), Error> {
    let found = find_variable(bytes, header, name)?;
    ensure(found.is_some(), BAD, missing)?;
    let variable = found.ok_or(Error::new(BAD, missing))?;
    ensure(variable.count > 0, BAD, "a wout variable is empty")?;
    let mut values = vec![0.0; variable.count];
    read_values(bytes, &variable, &mut values)?;
    Ok((values, variable))
}

fn read_scalar(bytes: &[u8], header: &Header, name: &str) -> Result<f64, Error> {
    let mut value = [f64::NAN];
    if let Some(variable) = find_variable(bytes, header, name)? {
        ensure(variable.count == 1, BAD, "a wout scalar has more than one value")?;
        ensure(variable.rank <= 1, BAD, "a wout scalar is an array")?;
        read_values(bytes, &variable, &mut value)?;
    }
    Ok(value[0])
}

fn read_reported(bytes: &[u8], header: &Header) -> Result<Reported, Error> {
    let reported = Reported {
        rmajor_p: read_scalar(bytes, header, "Rmajor_p")?,
        aminor_p: read_scalar(bytes, header, "Aminor_p")?,
        aspect: read_scalar(bytes, header, "aspect")?,
        volume_p: read_scalar(bytes, header, "volume_p")?,
        b0: read_scalar(bytes, header, "b0")?,
        volavg_b: read_scalar(bytes, header, "volavgB")?,
        betatotal: read_scalar(bytes, header, "betatotal")?,
        rbtor: read_scalar(bytes, header, "rbtor")?,
        ier_flag: read_scalar(bytes, header, "ier_flag")?,
    };
    ensure(
        reported.rmajor_p.is_nan() || reported.rmajor_p > 0.0,
        BAD,
        "wout file reports a major radius that is not positive",
    )?;
    ensure(
        reported.aminor_p.is_nan() || reported.aminor_p > 0.0,
        BAD,
        "wout file reports a minor radius that is not positive",
    )?;
    Ok(reported)
}

/// Field strength on the axis: mode numbers, cosine and sine coefficients.
///
/// The field strength lives on the half grid: index 0 is unused, index `j`
/// is centred between full-grid points `j - 1` and `j`. Only `m = 0` modes
/// are non-zero on the axis, so those are extrapolated linearly to `s = 0`.
fn read_axis_field(bytes: &[u8], header: &Header, lasym: bool) -> Result<AxisModes, Error> {
    let (xm, _) = read_array(bytes, header, "xm_nyq", "wout file has no variable 'xm_nyq'")?;
    let (xn, _) = read_array(bytes, header, "xn_nyq", "wout file has no variable 'xn_nyq'")?;
    let (b_cos, shape) = read_array(bytes, header, "bmnc", "wout file has no variable 'bmnc'")?;
    let modes = xm.len();
    ensure(
        xn.len() == modes && shape.rank == 2 && shape.shape[1] == modes,
        BAD,
        "inconsistent array shapes of the field strength in wout file",
    )?;
    ensure(shape.shape[0] >= 3, BAD, "wout file has fewer than 3 radial points")?;
    let has_sine = lasym && find_variable(bytes, header, "bmns")?.is_some();
    let b_sin = if has_sine {
        read_array(bytes, header, "bmns", "wout file has no variable 'bmns'")?.0
    } else {
        Vec::new()
    };
    ensure(!has_sine || b_sin.len() == b_cos.len(), BAD, "'bmns' has the wrong size")?;
    let mut axis = AxisModes { xn: Vec::new(), cos: Vec::new(), sin: Vec::new() };
    for mode in 0..modes {
        if xm[mode] == 0.0 {
            axis.xn.push(xn[mode]);
            axis.cos.push(1.5 * b_cos[modes + mode] - 0.5 * b_cos[2 * modes + mode]);
            axis.sin.push(if has_sine {
                1.5 * b_sin[modes + mode] - 0.5 * b_sin[2 * modes + mode]
            } else {
                0.0
            });
        }
    }
    ensure(all_finite(&axis.cos), BAD, "non-finite field strength in wout file")?;
    Ok(axis)
}

fn read_asymmetric(
    bytes: &[u8],
    header: &Header,
    lasym: bool,
    size: usize,
) -> Result<(Vec<f64>, Vec<f64>), Error> {
    ensure(size > 0, BAD, "wout file has no geometry coefficients")?;
    if !lasym {
        return Ok((Vec::new(), Vec::new()));
    }
    let (r_sin, _) = read_array(bytes, header, "rmns", "wout file has no variable 'rmns'")?;
    let (z_cos, _) = read_array(bytes, header, "zmnc", "wout file has no variable 'zmnc'")?;
    ensure(
        r_sin.len() == size && z_cos.len() == size,
        BAD,
        "inconsistent array shapes in wout file",
    )?;
    ensure(
        all_finite(&r_sin) && all_finite(&z_cos),
        BAD,
        "non-finite geometry coefficients in wout file",
    )?;
    Ok((r_sin, z_cos))
}

/// Read a VMEC `wout` file held in memory.
///
/// # Errors
/// Returns an error if the bytes are not a `NetCDF` classic file, a required
/// variable is missing, array shapes are inconsistent, or values are not
/// finite.
pub fn read_wout(bytes: &[u8]) -> Result<Equilibrium, Error> {
    let header = parse_header(bytes)?;
    let periods = read_scalar(bytes, &header, "nfp")?;
    ensure(periods >= 1.0, BAD, "wout file has no valid number of field periods")?;
    let nfp = floor_count(periods, MAX_PERIODS)?;
    let lasym = read_scalar(bytes, &header, "lasym__logical__")? > 0.0;
    let (xm, _) = read_array(bytes, &header, "xm", "wout file has no variable 'xm'")?;
    let (xn, _) = read_array(bytes, &header, "xn", "wout file has no variable 'xn'")?;
    let (r_cos, shape) = read_array(bytes, &header, "rmnc", "wout file has no variable 'rmnc'")?;
    let (z_sin, _) = read_array(bytes, &header, "zmns", "wout file has no variable 'zmns'")?;
    let n_modes = xm.len();
    ensure(
        shape.rank == 2 && shape.shape[1] == n_modes && xn.len() == n_modes,
        BAD,
        "inconsistent array shapes in wout file",
    )?;
    ensure(z_sin.len() == r_cos.len(), BAD, "inconsistent array shapes in wout file")?;
    let n_radial = shape.shape[0];
    ensure(n_radial >= 3, BAD, "wout file has fewer than 3 radial points")?;
    ensure(
        n_radial <= MAX_RADIAL && n_modes <= MAX_MODES,
        Code::Limit,
        "the equilibrium has more radial points or modes than this program supports",
    )?;
    ensure(
        all_finite(&r_cos) && all_finite(&z_sin) && all_finite(&xm) && all_finite(&xn),
        BAD,
        "non-finite geometry coefficients in wout file",
    )?;
    let (r_sin, z_cos) = read_asymmetric(bytes, &header, lasym, r_cos.len())?;
    let axis = read_axis_field(bytes, &header, lasym)?;
    Ok(Equilibrium {
        nfp,
        lasym,
        n_radial,
        n_modes,
        xm,
        xn,
        r_cos,
        z_sin,
        r_sin,
        z_cos,
        axis_xn: axis.xn,
        axis_cos: axis.cos,
        axis_sin: axis.sin,
        reported: read_reported(bytes, &header)?,
    })
}

/// A coefficient buffer for an equilibrium with `n_modes` modes.
///
/// # Errors
/// Returns an error if `n_modes` is zero or above the limit.
pub fn new_coefficients(n_modes: usize) -> Result<Coefficients, Error> {
    ensure(n_modes > 0, BAD, "the equilibrium has no Fourier modes")?;
    ensure(n_modes <= MAX_MODES, Code::Limit, "too many Fourier modes")?;
    Ok(Coefficients {
        rc: vec![0.0; n_modes],
        zs: vec![0.0; n_modes],
        rs: vec![0.0; n_modes],
        zc: vec![0.0; n_modes],
    })
}

/// Fourier coefficients at normalised flux `s`, by linear interpolation
/// between the two neighbouring grid surfaces.
///
/// # Errors
/// Returns an error if `s` is outside `[0, 1]` or `out` has the wrong size.
pub fn interpolate(eq: &Equilibrium, s: f64, out: &mut Coefficients) -> Result<(), Error> {
    ensure((0.0..=1.0).contains(&s), Code::Numeric, "s must be in [0, 1]")?;
    ensure(
        out.rc.len() == eq.n_modes && out.zs.len() == eq.n_modes,
        Code::Numeric,
        "the coefficient buffer does not match the equilibrium",
    )?;
    ensure(eq.n_radial >= 2, BAD, "the equilibrium has fewer than 2 surfaces")?;
    let position = s * real(eq.n_radial - 1)?;
    let inner = floor_count(position, eq.n_radial - 2)?;
    let weight = position - real(inner)?;
    let (low, high) = (inner * eq.n_modes, (inner + 1) * eq.n_modes);
    for mode in 0..eq.n_modes {
        out.rc[mode] = (1.0 - weight) * eq.r_cos[low + mode] + weight * eq.r_cos[high + mode];
        out.zs[mode] = (1.0 - weight) * eq.z_sin[low + mode] + weight * eq.z_sin[high + mode];
        if eq.lasym {
            out.rs[mode] = (1.0 - weight) * eq.r_sin[low + mode] + weight * eq.r_sin[high + mode];
            out.zc[mode] = (1.0 - weight) * eq.z_cos[low + mode] + weight * eq.z_cos[high + mode];
        }
    }
    Ok(())
}

/// Evaluate a flux surface, given by its coefficients, at one pair of angles.
///
/// # Errors
/// Returns an error if an angle is not finite, the buffer does not match the
/// equilibrium, or the result is not finite.
pub fn surface_point(
    eq: &Equilibrium,
    surface: &Coefficients,
    theta: f64,
    phi: f64,
) -> Result<SurfacePoint, Error> {
    ensure(theta.is_finite() && phi.is_finite(), Code::Numeric, "surface angles must be finite")?;
    ensure(
        surface.rc.len() == eq.n_modes && eq.xm.len() == eq.n_modes,
        Code::Numeric,
        "the coefficient buffer does not match the equilibrium",
    )?;
    let mut point = SurfacePoint::default();
    for mode in 0..eq.n_modes {
        let (m, n) = (eq.xm[mode], eq.xn[mode]);
        let (sin, cos) = (m * theta - n * phi).sin_cos();
        let (rc, zs) = (surface.rc[mode], surface.zs[mode]);
        point.r += rc * cos;
        point.z += zs * sin;
        point.r_t -= rc * m * sin;
        point.r_p += rc * n * sin;
        point.z_t += zs * m * cos;
        point.z_p -= zs * n * cos;
        if eq.lasym {
            let (rs, zc) = (surface.rs[mode], surface.zc[mode]);
            point.r += rs * sin;
            point.z += zc * cos;
            point.r_t += rs * m * cos;
            point.r_p -= rs * n * cos;
            point.z_t -= zc * m * sin;
            point.z_p += zc * n * sin;
        }
    }
    ensure(
        point.r.is_finite() && point.z.is_finite(),
        Code::Numeric,
        "a flux surface point is not finite",
    )?;
    Ok(point)
}

/// Cartesian position of a surface point at toroidal angle `phi`.
#[must_use]
pub fn cartesian(point: &SurfacePoint, phi: f64) -> Vec3 {
    [point.r * phi.cos(), point.r * phi.sin(), point.z]
}

/// Field strength on the magnetic axis [T] at toroidal angle `phi`.
///
/// # Errors
/// Returns an error if `phi` or the result is not finite.
pub fn axis_field_strength(eq: &Equilibrium, phi: f64) -> Result<f64, Error> {
    ensure(phi.is_finite(), Code::Numeric, "the toroidal angle must be finite")?;
    let mut strength = 0.0;
    for mode in 0..eq.axis_xn.len() {
        let (sin, cos) = (-eq.axis_xn[mode] * phi).sin_cos();
        strength += eq.axis_cos[mode] * cos + eq.axis_sin[mode] * sin;
    }
    ensure(strength.is_finite(), Code::Numeric, "the axis field is not finite")?;
    Ok(strength)
}

/// Multiply all lengths and all fields of the equilibrium by the factors.
///
/// # Errors
/// Returns an error if a factor is not positive and finite.
pub fn scale_equilibrium(eq: &mut Equilibrium, length: f64, field: f64) -> Result<(), Error> {
    ensure(
        length.is_finite() && length > 0.0,
        Code::Numeric,
        "the length factor must be positive",
    )?;
    ensure(field.is_finite() && field > 0.0, Code::Numeric, "the field factor must be positive")?;
    for value in eq.r_cos.iter_mut().chain(eq.z_sin.iter_mut()) {
        *value *= length;
    }
    for value in eq.r_sin.iter_mut().chain(eq.z_cos.iter_mut()) {
        *value *= length;
    }
    for value in eq.axis_cos.iter_mut().chain(eq.axis_sin.iter_mut()) {
        *value *= field;
    }
    eq.reported.rmajor_p *= length;
    eq.reported.aminor_p *= length;
    eq.reported.volume_p *= length * length * length;
    eq.reported.b0 *= field;
    eq.reported.volavg_b *= field;
    Ok(())
}
