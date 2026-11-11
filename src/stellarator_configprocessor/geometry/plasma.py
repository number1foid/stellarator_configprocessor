"""Plasma geometry from flux-surface Fourier coefficients.

All integrals use the trapezoidal rule on a uniform angle grid over one field
period, which converges spectrally for these periodic integrands.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..io.equilibrium import Equilibrium


@dataclass(frozen=True)
class SurfaceGeometry:
    """Integral properties of one flux surface."""

    s: float  # normalised toroidal flux of the surface
    volume: float  # enclosed volume [m^3]
    area: float  # surface area [m^2]
    cross_section: float  # toroidally averaged poloidal cross-section [m^2]

    @property
    def minor_radius(self) -> float:
        """Effective minor radius [m]: radius of a circle with the mean cross-section."""
        return float(np.sqrt(self.cross_section / np.pi))

    @property
    def major_radius(self) -> float:
        """Effective major radius [m], defined so that ``V = 2 pi^2 R a^2``."""
        return self.volume / (2.0 * np.pi * self.cross_section)

    @property
    def aspect_ratio(self) -> float:
        return self.major_radius / self.minor_radius


def angle_grid(eq: Equilibrium, n_theta: int | None = None, n_phi: int | None = None):
    """Uniform angles covering one field period, fine enough for the mode content."""
    m_max = int(np.max(np.abs(eq.xm)))
    n_max = int(np.max(np.abs(eq.xn)) // eq.nfp)
    n_theta = n_theta or max(64, 6 * m_max + 2)
    n_phi = n_phi or max(64, 6 * n_max + 2)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    phi = np.linspace(0.0, 2.0 * np.pi / eq.nfp, n_phi, endpoint=False)
    return theta, phi


def surface_geometry(eq: Equilibrium, s: float = 1.0, n_theta: int | None = None,
                     n_phi: int | None = None) -> SurfaceGeometry:
    """Volume, area and mean cross-section of the flux surface at ``s``.

    * Volume: ``V = | integral of (1/2) R^2 dZ/dtheta dtheta dphi |`` (divergence
      theorem; equivalent to eq. (3) of Lion et al. 2021).
    * Area: ``S = integral of |r_theta x r_phi| dtheta dphi`` with
      ``|r_theta x r_phi|^2 = R^2 (R_t^2 + Z_t^2) + (R_t Z_p - R_p Z_t)^2``
      (eq. (4) of the paper).
    * Cross-section at fixed ``phi``: ``| integral of R dZ/dtheta dtheta |``.
    """
    theta, phi = angle_grid(eq, n_theta, n_phi)
    g = eq.surface(s, theta, phi)
    d_theta = 2.0 * np.pi / len(theta)
    d_phi = 2.0 * np.pi / eq.nfp / len(phi)

    volume = abs(float(np.sum(0.5 * g["R"] ** 2 * g["Z_t"]))) * d_theta * d_phi * eq.nfp
    jac = np.sqrt(g["R"] ** 2 * (g["R_t"] ** 2 + g["Z_t"] ** 2)
                  + (g["R_t"] * g["Z_p"] - g["R_p"] * g["Z_t"]) ** 2)
    area = float(np.sum(jac)) * d_theta * d_phi * eq.nfp
    cross = float(np.mean(np.abs(np.sum(g["R"] * g["Z_t"], axis=0)))) * d_theta
    return SurfaceGeometry(s=s, volume=volume, area=area, cross_section=cross)
