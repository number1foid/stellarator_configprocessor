"""Electromagnetic forces on the coils.

The force per unit length of the coil centre line is the Lorentz force
density ``j x B`` integrated over the winding pack cross-section:

    F' = I * sum over the cross-section of  w * (dP/ds) x B

where ``P`` is a point of the cross-section, ``dP/ds`` its derivative along
the centre line and ``w`` the quadrature weight. ``dP/ds`` is the local
current direction times the ratio of the local to the centre-line length:
on the outside of a bend the conductor is longer than its centre line. That
factor matters for the coil's own field, which changes sign across the
cross-section, and is needed to obtain the correct hoop force.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry.coils import SymmetryClasses
from .io.coils import SampledCoils
from .winding_pack import PackField


@dataclass(frozen=True)
class ForceResult:
    """Forces on the representative coil of each symmetry class."""

    line_force: tuple  # per class: [n, 3] force per unit length [N/m]
    radial: tuple  # per class: [n] component along the in-plane normal [N/m]
    lateral: tuple  # per class: [n] component along the toroidal normal [N/m]
    net_force: np.ndarray  # [nclass, 3] total force on the coil [N]
    centering: np.ndarray  # [nclass] net force towards the torus axis [N]
    multiplicity: np.ndarray  # [nclass]
    wp_area: float

    @property
    def max_line_force(self) -> float:
        """Largest force per unit length [N/m]."""
        return float(max(np.linalg.norm(f, axis=1).max() for f in self.line_force))

    @property
    def max_force_density(self) -> float:
        """Largest cross-section-averaged force density [N/m^3]."""
        return self.max_line_force / self.wp_area

    @property
    def max_radial_force_density(self) -> float:
        return float(max(np.abs(f).max() for f in self.radial)) / self.wp_area

    @property
    def max_lateral_force_density(self) -> float:
        return float(max(np.abs(f).max() for f in self.lateral)) / self.wp_area

    @property
    def centering_max(self) -> float:
        return float(self.centering.max())

    @property
    def centering_min(self) -> float:
        return float(self.centering.min())

    @property
    def centering_avg(self) -> float:
        return float(np.sum(self.centering * self.multiplicity) / np.sum(self.multiplicity))


def coil_forces(coils: SampledCoils, pack: PackField, frames, classes: SymmetryClasses
                ) -> ForceResult:
    """Forces from the winding pack field of the class representatives.

    ``pack`` must hold the field of every representative at every sample.
    """
    if not np.array_equal(pack.targets, classes.representatives):
        raise ValueError("the winding pack field must cover the class representatives")
    n_r, n_t = frames
    line, radial, lateral, net, centering = [], [], [], [], []
    for c, i in enumerate(classes.representatives):
        if len(pack.stations[c]) != coils.n_points:
            raise ValueError("forces need the field at every sample of each coil")
        pts = pack.points[c]
        dp_ds = (np.roll(pts, -1, axis=0) - np.roll(pts, 1, axis=0)) / (2.0 * coils.ds[i])
        f = coils.currents[i] * np.einsum("ab,sabk->sk", pack.weights,
                                          np.cross(dp_ds, pack.b[c]))
        line.append(f)
        radial.append(np.einsum("sk,sk->s", f, n_r[i]))
        lateral.append(np.einsum("sk,sk->s", f, n_t[i]))
        total = f.sum(axis=0) * coils.ds[i]
        net.append(total)
        centroid = coils.points[i].mean(axis=0)
        e_r = np.array([centroid[0], centroid[1], 0.0]) / np.hypot(centroid[0], centroid[1])
        centering.append(-float(total @ e_r))
    return ForceResult(line_force=tuple(line), radial=tuple(radial), lateral=tuple(lateral),
                       net_force=np.array(net), centering=np.array(centering),
                       multiplicity=classes.multiplicity.copy(), wp_area=pack.wp.area)
