"""Check that a coil set and an equilibrium describe the same machine.

Four tests:

1. **Field periods.** The number declared in the coil file, if any, equals
   the equilibrium's, and the coil count is a multiple of it.
2. **Linkage.** By Ampere's law the line integral of the coil field along the
   magnetic axis equals ``mu0`` times the current linked by the axis. Every
   coil of a modular set must link the axis in the same sense, so the
   integral must equal ``mu0 * sum(|I|)``.
3. **Axis position.** On the true magnetic axis the field is tangent to the
   axis. The coil field at the equilibrium's axis is tested for this; a
   displaced or differently shaped axis shows up as a perpendicular
   component.
4. **Field on axis.** The mean field strength of the coils along the axis is
   compared with the equilibrium's.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._checks import require
from .constants import MU0
from .field import biot_savart
from .io.coils import SampledCoils
from .io.equilibrium import Equilibrium


class ConsistencyError(RuntimeError):
    """Coils and equilibrium are inconsistent and no override was given."""


@dataclass(frozen=True)
class Tolerances:
    """Limits for the consistency tests. All are relative."""

    linkage: float = 0.01
    #: Mean perpendicular field on the axis divided by the mean field.
    axis_misalignment: float = 0.02
    #: Relative difference of the mean field on axis.
    field: float = 0.03


@dataclass(frozen=True)
class AxisField:
    """Coil field sampled along the equilibrium's magnetic axis."""

    phi: np.ndarray  # [n] toroidal angle of each sample
    points: np.ndarray  # [n, 3]
    dl: np.ndarray  # [n, 3] axis line elements
    b: np.ndarray  # [n, 3] coil field [T]

    @property
    def arc(self) -> np.ndarray:
        return np.linalg.norm(self.dl, axis=1)

    @property
    def mean_strength(self) -> float:
        """Arc-length average of ``|B|`` along the axis [T]."""
        return float(np.sum(np.linalg.norm(self.b, axis=1) * self.arc) / np.sum(self.arc))

    @property
    def mean_toroidal(self) -> float:
        """Arc-length average of the cylindrical toroidal component [T]."""
        e_phi = np.column_stack([-np.sin(self.phi), np.cos(self.phi), np.zeros_like(self.phi)])
        return float(np.sum(np.einsum("ij,ij->i", self.b, e_phi) * self.arc) / np.sum(self.arc))

    @property
    def circulation(self) -> float:
        """Line integral of ``B . dl`` around the axis [T m]."""
        return float(np.sum(np.einsum("ij,ij->i", self.b, self.dl)))


def axis_field(eq: Equilibrium, coils: SampledCoils, n_phi: int = 720) -> AxisField:
    """Evaluate the coil field along the full magnetic axis."""
    require(n_phi >= 16, "the axis needs at least 16 points")
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    pts = eq.axis_xyz(phi)
    # Central differences on the closed axis give second-order line elements.
    dl = 0.5 * (np.roll(pts, -1, axis=0) - np.roll(pts, 1, axis=0))
    start, end, cur = coils.segments()
    return AxisField(phi=phi, points=pts, dl=dl, b=biot_savart(pts, start, end, cur))


@dataclass(frozen=True)
class ConsistencyResult:
    nfp_equilibrium: int
    nfp_coil_file: int | None
    n_coils: int
    linkage_ratio: float
    axis_misalignment: float
    b_axis_coils: float
    b_axis_equilibrium: float
    field_ratio: float
    failures: tuple[str, ...] = field(default_factory=tuple)

    @property
    def ok(self) -> bool:
        return not self.failures

    def as_dict(self) -> dict:
        return {
            "nfp_equilibrium": self.nfp_equilibrium,
            "nfp_coil_file": self.nfp_coil_file,
            "n_coils": self.n_coils,
            "linkage_ratio": self.linkage_ratio,
            "axis_misalignment": self.axis_misalignment,
            "b_axis_coils_T": self.b_axis_coils,
            "b_axis_equilibrium_T": self.b_axis_equilibrium,
            "field_ratio_coils_over_equilibrium": self.field_ratio,
            "passed": self.ok,
            "failures": list(self.failures),
        }


def check_consistency(
    eq: Equilibrium,
    coils: SampledCoils,
    periods_in_file: int | None = None,
    tolerances: Tolerances | None = None,
    axis: AxisField | None = None,
    check_field: bool = True,
) -> ConsistencyResult:
    """Run the consistency tests. Never raises; see :func:`enforce`.

    ``check_field=False`` skips the field-strength comparison. It is used
    when the coil currents are rescaled to match the equilibrium anyway.
    """
    tolerances = tolerances or Tolerances()
    failures: list[str] = []
    n = coils.n_coils

    if periods_in_file is not None and periods_in_file != eq.nfp:
        failures.append(f"field periods differ: coil file declares {periods_in_file}, "
                        f"equilibrium has {eq.nfp}")
    if n % eq.nfp != 0:
        failures.append(f"{n} coils cannot be divided into {eq.nfp} field periods")

    axis = axis if axis is not None else axis_field(eq, coils)
    total_current = float(np.sum(np.abs(coils.currents)))
    if total_current == 0.0:
        failures.append("all coil currents are zero")
        linkage = float("nan")
    else:
        linkage = abs(axis.circulation) / (MU0 * total_current)
        if abs(linkage - 1.0) > tolerances.linkage:
            failures.append(
                f"the magnetic axis links {linkage:.4f} of the total coil current "
                "(expected 1): some coils do not enclose the axis, or coil currents "
                "have opposing directions")

    strength = np.linalg.norm(axis.b, axis=1)
    arc = axis.arc
    tangent = axis.dl / arc[:, None]
    parallel = np.einsum("ij,ij->i", axis.b, tangent)
    perp = np.linalg.norm(axis.b - parallel[:, None] * tangent, axis=1)
    b_coils = axis.mean_strength
    misalignment = (float(np.sum(perp * arc) / np.sum(strength * arc)) if b_coils > 0
                    else float("nan"))
    if not misalignment <= tolerances.axis_misalignment:
        failures.append(
            f"coil field is not tangent to the equilibrium's magnetic axis: mean "
            f"perpendicular fraction {misalignment:.4f} exceeds {tolerances.axis_misalignment}; "
            "the axis position or shape does not match the coils")

    b_eq = float(np.sum(np.abs(eq.axis_field_strength(axis.phi)) * arc) / np.sum(arc))
    ratio = b_coils / b_eq if b_eq > 0 else float("nan")
    if check_field and not abs(ratio - 1.0) <= tolerances.field:
        failures.append(
            f"field on axis differs: coils give {b_coils:.4f} T, equilibrium has "
            f"{b_eq:.4f} T (ratio {ratio:.4f}, tolerance {tolerances.field})")

    return ConsistencyResult(
        nfp_equilibrium=eq.nfp, nfp_coil_file=periods_in_file, n_coils=n,
        linkage_ratio=linkage, axis_misalignment=misalignment, b_axis_coils=b_coils,
        b_axis_equilibrium=b_eq, field_ratio=ratio, failures=tuple(failures),
    )


def enforce(result: ConsistencyResult, allow_inconsistent: bool = False) -> None:
    """Raise :class:`ConsistencyError` if the check failed and no override is set."""
    if result.ok or allow_inconsistent:
        return
    raise ConsistencyError(
        "coils and equilibrium appear inconsistent:\n  - "
        + "\n  - ".join(result.failures)
        + "\nUse the override option to continue anyway.")
