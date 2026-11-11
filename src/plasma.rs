//! Integral properties of a flux surface.
//!
//! * Volume: `V = | integral of (1/2) R^2 dZ/dtheta dtheta dphi |` (divergence
//!   theorem; equivalent to eq. (3) of Lion et al. 2021).
//! * Area: `S = integral of |r_theta x r_phi| dtheta dphi` with
//!   `|r_theta x r_phi|^2 = R^2 (R_t^2 + Z_t^2) + (R_t Z_p - R_p Z_t)^2`
//!   (eq. (4) of the paper).
//! * Cross-section at fixed `phi`: `| integral of R dZ/dtheta dtheta |`.

use crate::equilibrium::{interpolate, surface_point, Coefficients, Equilibrium};
use crate::error::{ensure, Code, Error};
use crate::limits::MAX_SURFACE;
use crate::num::{floor_count, real};
use core::f64::consts::{PI, TAU};

const BAD: Code = Code::Numeric;

/// Integral properties of one flux surface.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceGeometry {
    /// Normalised toroidal flux of the surface.
    pub s: f64,
    /// Enclosed volume [m^3].
    pub volume: f64,
    /// Surface area [m^2].
    pub area: f64,
    /// Toroidally averaged poloidal cross-section [m^2].
    pub cross_section: f64,
    /// Effective minor radius [m]: radius of a circle with the mean cross-section.
    pub minor_radius: f64,
    /// Effective major radius [m], defined so that `V = 2 pi^2 R a^2`.
    pub major_radius: f64,
    /// Major radius divided by minor radius.
    pub aspect_ratio: f64,
}

/// Number of poloidal and toroidal angles that resolve the mode content of
/// the equilibrium over one field period.
///
/// # Errors
/// Returns an error if the mode numbers are not finite or too high for the
/// largest supported grid.
pub fn angle_counts(eq: &Equilibrium) -> Result<(usize, usize), Error> {
    ensure(eq.nfp >= 1, Code::Equilibrium, "the equilibrium has no field periods")?;
    let (mut poloidal, mut toroidal) = (0.0_f64, 0.0_f64);
    for (m, n) in eq.xm.iter().zip(eq.xn.iter()) {
        poloidal = poloidal.max(m.abs());
        toroidal = toroidal.max(n.abs());
    }
    let limit = (MAX_SURFACE - 2) / 6;
    let m_max = floor_count(poloidal, limit + 1)?;
    let n_max = floor_count(toroidal, (limit + 1) * eq.nfp)? / eq.nfp;
    ensure(
        m_max <= limit && n_max <= limit,
        Code::Limit,
        "the equilibrium has higher mode numbers than this program supports",
    )?;
    Ok(((6 * m_max + 2).max(64), (6 * n_max + 2).max(64)))
}

/// Volume, area and mean cross-section of the flux surface at `s`.
///
/// # Errors
/// Returns an error if `s` is outside `[0, 1]`, the surface cannot be
/// evaluated, or it encloses no volume.
pub fn surface_geometry(
    eq: &Equilibrium,
    s: f64,
    surface: &mut Coefficients,
) -> Result<SurfaceGeometry, Error> {
    let (n_theta, n_phi) = angle_counts(eq)?;
    interpolate(eq, s, surface)?;
    let d_theta = TAU / real(n_theta)?;
    let d_phi = TAU / real(eq.nfp)? / real(n_phi)?;
    let (mut volume_sum, mut area_sum, mut cross_sum) = (0.0, 0.0, 0.0);
    for toroidal in 0..n_phi {
        let phi = real(toroidal)? * d_phi;
        let mut column = 0.0;
        for poloidal in 0..n_theta {
            let p = surface_point(eq, surface, real(poloidal)? * d_theta, phi)?;
            let twist = p.r_t * p.z_p - p.r_p * p.z_t;
            volume_sum += 0.5 * p.r * p.r * p.z_t;
            area_sum += (p.r * p.r * (p.r_t * p.r_t + p.z_t * p.z_t) + twist * twist).sqrt();
            column += p.r * p.z_t;
        }
        cross_sum += column.abs();
    }
    let periods = real(eq.nfp)?;
    let volume = volume_sum.abs() * d_theta * d_phi * periods;
    let cross_section = cross_sum / real(n_phi)? * d_theta;
    ensure(volume > 0.0 && volume.is_finite(), BAD, "a flux surface encloses no volume")?;
    ensure(cross_section > 0.0, BAD, "a flux surface has no cross-section")?;
    let minor_radius = (cross_section / PI).sqrt();
    let major_radius = volume / (TAU * cross_section);
    Ok(SurfaceGeometry {
        s,
        volume,
        area: area_sum * d_theta * d_phi * periods,
        cross_section,
        minor_radius,
        major_radius,
        aspect_ratio: major_radius / minor_radius,
    })
}
