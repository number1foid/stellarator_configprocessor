"""Coil filaments: MAKEGRID ``coils.*`` reader, symmetry expansion, sampling.

MAKEGRID file layout (after an optional Fortran namelist header)::

    periods 5
    begin filament
    mirror NIL
    x y z current
    ...
    x y z 0 group_id group_name      <- closes the coil
    R 0 z current group_id name      <- one-line coil: circle around the z axis
    end

The file lists every coil of the full torus. The current of a coil is the
value in the fourth column, optionally multiplied by a per-group factor (the
``extcur`` array of VMEC).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from ..geometry.curves import CurveError, SampledCurve, clean_closed_points, sample_closed_curve


class CoilFileError(ValueError):
    """The coil file is missing, malformed or not usable."""


@dataclass(frozen=True)
class Coil:
    """One closed filament."""

    points: np.ndarray  # [n, 3] distinct points; the curve closes implicitly
    current: float  # [A] (ampere-turns)
    group: int = 0
    name: str = ""

    def scaled(self, length_factor: float = 1.0, current_factor: float = 1.0) -> "Coil":
        return replace(self, points=self.points * length_factor,
                       current=self.current * current_factor)


@dataclass(frozen=True)
class CoilSet:
    """All coils of the full torus."""

    coils: tuple[Coil, ...]
    periods: int | None = None  # number of field periods declared in the file
    source: str = ""
    warnings: tuple[str, ...] = field(default_factory=tuple)

    def __len__(self) -> int:
        return len(self.coils)

    @property
    def currents(self) -> np.ndarray:
        return np.array([c.current for c in self.coils])

    def scaled(self, length_factor: float = 1.0, current_factor: float = 1.0) -> "CoilSet":
        return replace(self, coils=tuple(c.scaled(length_factor, current_factor)
                                         for c in self.coils))

    def without_zero_current(self) -> "CoilSet":
        return replace(self, coils=tuple(c for c in self.coils if c.current != 0.0))

    def sampled(self, n: int) -> "SampledCoils":
        """Sample every coil at ``n`` points uniform in arc length."""
        curves = [sample_closed_curve(c.points, n) for c in self.coils]
        return SampledCoils(
            points=np.stack([c.points for c in curves]),
            tangents=np.stack([c.tangents for c in curves]),
            curvature=np.stack([c.curvature for c in curves]),
            lengths=np.array([c.length for c in curves]),
            currents=self.currents,
        )


@dataclass(frozen=True)
class SampledCoils:
    """Coils sampled with the same number of points, as stacked arrays."""

    points: np.ndarray  # [ncoil, n, 3]
    tangents: np.ndarray  # [ncoil, n, 3]
    curvature: np.ndarray  # [ncoil, n]
    lengths: np.ndarray  # [ncoil]
    currents: np.ndarray  # [ncoil]

    @property
    def n_coils(self) -> int:
        return self.points.shape[0]

    @property
    def n_points(self) -> int:
        return self.points.shape[1]

    @property
    def ds(self) -> np.ndarray:
        """Arc length per sample of each coil, shape ``[ncoil]``."""
        return self.lengths / self.n_points

    def segments(self, coil_indices=None):
        """Straight segments joining consecutive samples.

        Returns ``(start, end, current)`` with shapes ``[nseg, 3]``,
        ``[nseg, 3]`` and ``[nseg]``.
        """
        idx = np.arange(self.n_coils) if coil_indices is None else np.asarray(coil_indices)
        start = self.points[idx]
        end = np.roll(start, -1, axis=1)
        cur = np.repeat(self.currents[idx], self.n_points)
        return start.reshape(-1, 3), end.reshape(-1, 3), cur


def _to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def read_makegrid(
    path: str | Path,
    group_currents: dict[int, float] | None = None,
    circle_points: int = 128,
) -> CoilSet:
    """Read a MAKEGRID ``coils.*`` file.

    ``group_currents`` maps a coil group id to a factor applied to the
    currents of that group (the VMEC ``extcur`` values). Groups that are not
    listed keep the value from the file.
    """
    path = Path(path)
    if not path.is_file():
        raise CoilFileError(f"coil file not found: {path}")
    lines = path.read_text(errors="replace").splitlines()

    periods = None
    start = None
    for i, line in enumerate(lines):
        tokens = line.split()
        if tokens and tokens[0].lower() == "periods":
            try:
                periods = int(tokens[1])
            except (IndexError, ValueError) as exc:
                raise CoilFileError(f"{path}:{i + 1}: bad 'periods' line") from exc
            start = i + 1
            break
    if start is None:
        raise CoilFileError(f"{path}: no 'periods' line found; not a MAKEGRID coils file")

    coils: list[Coil] = []
    warnings: list[str] = []
    pending: list[list[float]] = []
    pending_currents: list[float] = []
    factor = group_currents or {}

    for lineno, line in enumerate(lines[start:], start=start + 1):
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0].lower()
        if key == "end":
            break
        if key in ("begin", "mirror"):
            continue
        if len(tokens) < 4:
            raise CoilFileError(f"{path}:{lineno}: expected 'x y z current', got: {line.strip()}")
        try:
            x, y, z, cur = (_to_float(t) for t in tokens[:4])
        except ValueError as exc:
            raise CoilFileError(f"{path}:{lineno}: cannot parse numbers: {line.strip()}") from exc

        if len(tokens) < 5:
            pending.append([x, y, z])
            pending_currents.append(cur)
            continue

        # A line with a group id and a name ends a coil.
        try:
            group = int(float(tokens[4]))
        except ValueError as exc:
            raise CoilFileError(f"{path}:{lineno}: bad group id: {line.strip()}") from exc
        name = " ".join(tokens[5:])
        scale = factor.get(group, 1.0)

        if not pending:
            # One-line coil: a circle around the z axis.
            if cur == 0.0:
                continue
            radius = float(np.hypot(x, y))
            ang = np.linspace(0.0, 2.0 * np.pi, circle_points, endpoint=False)
            pts = np.column_stack([radius * np.cos(ang), radius * np.sin(ang), np.full_like(ang, z)])
            coils.append(Coil(points=pts, current=cur * scale, group=group, name=name))
            continue

        pending.append([x, y, z])
        currents = np.array(pending_currents)
        if cur != 0.0:
            warnings.append(f"line {lineno}: closing point of coil '{name}' has non-zero current")
        if np.ptp(currents) > 1e-9 * max(1.0, float(np.max(np.abs(currents)))):
            warnings.append(f"coil ending at line {lineno} ('{name}') has varying current; "
                            "the first value is used")
        pts = np.array(pending)
        closure = float(np.linalg.norm(pts[-1] - pts[0]))
        extent = float(np.max(np.ptp(pts, axis=0)))
        if closure > 1e-6 * extent:
            warnings.append(f"coil ending at line {lineno} ('{name}') is not closed "
                            f"(gap {closure:.3g} m); it is closed with a straight segment")
        try:
            pts = clean_closed_points(pts)
        except CurveError as exc:
            raise CoilFileError(f"{path}: coil ending at line {lineno}: {exc}") from exc
        coils.append(Coil(points=pts, current=float(currents[0]) * scale, group=group, name=name))
        pending, pending_currents = [], []

    if pending:
        raise CoilFileError(f"{path}: file ends inside a coil (no closing line with group id)")
    if not coils:
        raise CoilFileError(f"{path}: no coils found")
    return CoilSet(coils=tuple(coils), periods=periods, source=str(path), warnings=tuple(warnings))


def write_makegrid(path: str | Path, coilset: CoilSet, periods: int) -> None:
    """Write coils in MAKEGRID format. Coils with zero current cannot be stored."""
    out = [f"periods {periods}", "begin filament", "mirror NIL"]
    for coil in coilset.coils:
        if coil.current == 0.0:
            raise CoilFileError("MAKEGRID cannot represent a coil with zero current")
        for p in coil.points:
            out.append(f"{p[0]: .15E} {p[1]: .15E} {p[2]: .15E} {coil.current: .15E}")
        p = coil.points[0]
        name = coil.name or "coil"
        out.append(f"{p[0]: .15E} {p[1]: .15E} {p[2]: .15E} {0.0: .15E} {coil.group} {name}")
    out.append("end")
    Path(path).write_text("\n".join(out) + "\n")


def expand_symmetry(base: CoilSet, nfp: int, stellsym: bool = True) -> CoilSet:
    """Build the full torus from the coils of one half period (or one period).

    Each coil is rotated into every field period. With stellarator symmetry
    each coil is also mirrored, ``(x, y, z) -> (x, -y, -z)``, with the current
    reversed, which keeps the field stellarator-symmetric.
    """
    coils = []
    for k in range(nfp):
        ang = 2.0 * np.pi * k / nfp
        rot = np.array([[np.cos(ang), -np.sin(ang), 0.0],
                        [np.sin(ang), np.cos(ang), 0.0],
                        [0.0, 0.0, 1.0]])
        for flip in ([False, True] if stellsym else [False]):
            for coil in base.coils:
                pts = coil.points
                cur = coil.current
                if flip:
                    pts = pts * np.array([1.0, -1.0, -1.0])
                    cur = -cur
                coils.append(replace(coil, points=pts @ rot.T, current=cur))
    return replace(base, coils=tuple(coils), periods=nfp)


__all__ = ["Coil", "CoilSet", "SampledCoils", "SampledCurve", "CoilFileError",
           "read_makegrid", "write_makegrid", "expand_symmetry"]
