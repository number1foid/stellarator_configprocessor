"""Magnetic energy and inductance of a coil set.

The energy is ``W = 1/2 * sum_ij M_ij I_i I_j`` with

* mutual inductances from Neumann's formula for the centre filaments,
  ``M_ij = mu0 / (4 pi) * double integral of dl_i . dl_j / |r_i - r_j|``;
* self-inductances from the same integral with the distance regularised by
  the geometric mean distance (GMD) of the rectangular winding pack,
  ``|r - r'|^2 -> |r - r'|^2 + GMD^2``. This is the classical result for a
  conductor with uniform current density and reproduces
  ``L = mu0 R (ln(8 R / GMD) - 2)`` for a circular loop.

PROCESS uses one number, the effective inductance ``L = 2 W / I^2`` with
``I`` the current per coil.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._checks import coil_index, points_array, positive, require
from .constants import MU0
from .geometry.coils import SymmetryClasses
from .geometry.curves import sample_closed_curve
from .io.coils import SampledCoils
from .winding_pack import WindingPack


def gmd_squared(wp: WindingPack) -> float:
    """Square of the geometric mean distance of a rectangle from itself [m^2].

    Closed form of ``exp(2 <ln d>)`` for two points uniformly distributed in
    a rectangle with sides ``a`` and ``b``. For a square, ``GMD = 0.44705 a``.
    """
    a, b = wp.radial, wp.toroidal
    k = (4.0 * b / (3.0 * a) * np.arctan(a / b)
         + 4.0 * a / (3.0 * b) * np.arctan(b / a)
         + b**2 / (6.0 * a**2) * np.log(b / a)
         + a**2 / (6.0 * b**2) * np.log(a / b)
         - (a**4 - 6.0 * a**2 * b**2 + b**4) / (6.0 * a**2 * b**2) * np.log(a / b + b / a))
    return float(a * b * np.exp(-25.0 / 6.0 + k))


def self_inductance(points: np.ndarray, wp: WindingPack, resolution: float = 4.0) -> float:
    """Self-inductance of one closed filament with a rectangular winding pack [H].

    The filament is resampled so that the spacing is ``1 / resolution`` of the
    GMD, which resolves the regularised kernel.
    """
    points = points_array(points, "points")
    require(points.ndim == 2 and len(points) >= 3, "a coil needs at least 3 points [n, 3]")
    resolution = positive(resolution, "resolution")
    eps2 = gmd_squared(wp)
    coarse = sample_closed_curve(points, max(16, len(points)))
    n = int(max(64, np.ceil(resolution * coarse.length / np.sqrt(eps2))))
    curve = sample_closed_curve(points, n)
    total = 0.0
    block = max(1, 2_000_000 // n)
    for i in range(0, n, block):
        d = curve.points[i:i + block, None, :] - curve.points[None, :, :]
        dist = np.sqrt(np.einsum("abk,abk->ab", d, d) + eps2)
        total += float(np.sum((curve.tangents[i:i + block] @ curve.tangents.T) / dist))
    return MU0 / (4.0 * np.pi) * total * curve.ds**2


def flux_from_others(coils: SampledCoils, i: int) -> float:
    """Flux through coil ``i`` produced by all other coils [Wb]."""
    i = coil_index(i, coils.n_coils)
    others = np.array([j for j in range(coils.n_coils) if j != i], dtype=int)
    if len(others) == 0:
        return 0.0
    pts_i, tan_i = coils.points[i], coils.tangents[i]
    pts_o = coils.points[others].reshape(-1, 3)
    # Tangent times line element times current of every source sample.
    weight = (coils.tangents[others]
              * (coils.ds[others] * coils.currents[others])[:, None, None]).reshape(-1, 3)
    total = 0.0
    block = max(1, 2_000_000 // len(pts_o))
    for k in range(0, len(pts_i), block):
        d = pts_i[k:k + block, None, :] - pts_o[None, :, :]
        dist = np.sqrt(np.einsum("abk,abk->ab", d, d))
        total += float(np.sum((tan_i[k:k + block] @ weight.T) / dist))
    return MU0 / (4.0 * np.pi) * total * coils.ds[i]


@dataclass(frozen=True)
class EnergyResult:
    energy: float  # total magnetic energy [J]
    reference_current: float  # current per coil used for the effective inductance [A]
    self_inductance: np.ndarray  # [nclass] self-inductance of each representative [H]
    flux: np.ndarray  # [nclass] total flux linked by each representative [Wb]

    @property
    def effective_inductance(self) -> float:
        """``2 W / I^2`` with ``I`` the reference current per coil [H]."""
        require(self.reference_current != 0.0, "the reference current is zero")
        return 2.0 * self.energy / self.reference_current**2


def magnetic_energy(coils: SampledCoils, wp: WindingPack,
                    classes: SymmetryClasses) -> EnergyResult:
    """Total magnetic energy of the coil set.

    Coils in the same symmetry class contribute equally, so only one
    representative per class is evaluated.
    """
    reps = classes.representatives
    require(len(reps) > 0 and len(reps) == len(classes.multiplicity),
            "symmetry classes need one multiplicity per representative")
    require(int(np.sum(classes.multiplicity)) == coils.n_coils,
            "symmetry classes do not cover the coil set")
    self_l = np.array([self_inductance(coils.points[int(i)], wp) for i in reps])
    flux = np.array([self_l[c] * coils.currents[int(i)] + flux_from_others(coils, int(i))
                     for c, i in enumerate(reps)])
    energy = 0.5 * float(np.sum(classes.multiplicity * coils.currents[reps] * flux))
    reference = float(np.mean(np.abs(coils.currents)))
    return EnergyResult(energy=energy, reference_current=reference,
                        self_inductance=self_l, flux=flux)
