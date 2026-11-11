"""Geometry of a modular coil set and its distance to the plasma.

All functions take :class:`~stellarator_configprocessor.io.coils.SampledCoils`,
in which every coil is sampled at the same number of points, uniform in arc
length.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .._checks import coil_index, points_array, require
from ..io.coils import SampledCoils
from ..io.equilibrium import Equilibrium
from .plasma import angle_grid, surface_geometry

_MAX_PAIRS = 2_000_000


# --------------------------------------------------------------------- basics
def coil_radii(coils: SampledCoils) -> tuple[float, float, float]:
    """Mean major radius, mean minor radius and mean height of the coil set.

    * major radius: arc-length-weighted mean cylindrical radius of all
      filament points;
    * minor radius: arc-length-weighted mean distance of the points from the
      circle ``R = major radius``, ``z = mean z``.

    Both are exact for circular coils centred on that circle.
    """
    weights = np.repeat(coils.ds[:, None], coils.n_points, axis=1)
    r = np.hypot(coils.points[..., 0], coils.points[..., 1])
    z = coils.points[..., 2]
    total = float(np.sum(weights))
    r_major = float(np.sum(weights * r)) / total
    z_mean = float(np.sum(weights * z)) / total
    r_minor = float(np.sum(weights * np.hypot(r - r_major, z - z_mean))) / total
    return r_major, r_minor, z_mean


def total_length(coils: SampledCoils) -> float:
    """Sum of the lengths of all coil centre filaments [m]."""
    return float(np.sum(coils.lengths))


def coil_height(coils: SampledCoils) -> float:
    """Full vertical extent of the coil set [m]."""
    z = coils.points[..., 2]
    return float(z.max() - z.min())


def min_bend_radius(coils: SampledCoils) -> float:
    """Smallest radius of curvature of any filament [m]."""
    return float(1.0 / np.max(coils.curvature))


# ------------------------------------------------------------------ distances
def point_segment_distance(points: np.ndarray, start: np.ndarray, end: np.ndarray):
    """Shortest distance from each point to a set of straight segments.

    Returns ``(distance [npts], segment index [npts])``.
    """
    points = points_array(points, "points").reshape(-1, 3)
    start = points_array(start, "start").reshape(-1, 3)
    end = points_array(end, "end").reshape(-1, 3)
    require(len(start) == len(end) and len(start) > 0,
            "start and end must hold the same, non-zero number of segments")
    seg = end - start
    seg_len2 = np.einsum("sk,sk->s", seg, seg)
    dist = np.empty(len(points))
    index = np.empty(len(points), dtype=int)
    step = max(1, _MAX_PAIRS // max(1, len(start)))
    for i in range(0, len(points), step):
        rel = points[i:i + step, None, :] - start[None, :, :]
        t = np.clip(np.einsum("psk,sk->ps", rel, seg) / seg_len2[None, :], 0.0, 1.0)
        gap = rel - t[:, :, None] * seg[None, :, :]
        d2 = np.einsum("psk,psk->ps", gap, gap)
        index[i:i + step] = np.argmin(d2, axis=1)
        dist[i:i + step] = np.sqrt(np.min(d2, axis=1))
    return dist, index


def coil_distances(coils: SampledCoils, i: int) -> np.ndarray:
    """Shortest distance from coil ``i`` to every coil, shape ``[ncoil]``.

    The entry for coil ``i`` itself is infinite.
    """
    i = coil_index(i, coils.n_coils)
    out = np.full(coils.n_coils, np.inf)
    others = np.array([j for j in range(coils.n_coils) if j != i])
    if len(others) == 0:
        return out
    start = coils.points[others]
    end = np.roll(start, -1, axis=1)
    n = coils.n_points
    rel_start, rel_end = start.reshape(-1, 3), end.reshape(-1, 3)
    # Distance from each point of coil i to each other coil separately.
    points = coils.points[i]
    seg = rel_end - rel_start
    seg_len2 = np.einsum("sk,sk->s", seg, seg)
    best = np.full(len(others), np.inf)
    step = max(1, _MAX_PAIRS // len(rel_start))
    for k in range(0, len(points), step):
        rel = points[k:k + step, None, :] - rel_start[None, :, :]
        t = np.clip(np.einsum("psk,sk->ps", rel, seg) / seg_len2[None, :], 0.0, 1.0)
        gap = rel - t[:, :, None] * seg[None, :, :]
        d2 = np.einsum("psk,psk->ps", gap, gap).reshape(-1, len(others), n)
        best = np.minimum(best, np.sqrt(d2.min(axis=(0, 2))))
    out[others] = best
    return out


# ------------------------------------------------------------------- symmetry
@dataclass(frozen=True)
class SymmetryClasses:
    """Coils grouped into sets that are images of one another."""

    representatives: np.ndarray  # [nclass] index of one coil per class
    multiplicity: np.ndarray  # [nclass] number of coils in each class
    class_of: np.ndarray  # [ncoil] class index of each coil


def symmetry_classes(coils: SampledCoils, nfp: int, stellsym: bool = True,
                     rtol: float = 1e-4) -> SymmetryClasses:
    """Group coils that map onto each other under the device symmetry.

    Two coils are equivalent if one is the other rotated by a multiple of
    ``2 pi / nfp``, optionally combined with the stellarator-symmetry flip
    ``(R, phi, z) -> (R, -phi, -z)``. Equivalence is tested on the centroid,
    the length and the absolute current. If the coil set does not have the
    symmetry, every coil forms its own class and results stay correct.
    """
    require(int(nfp) >= 1, f"nfp must be at least 1, got {nfp!r}")
    require(rtol > 0.0, "rtol must be positive")
    centroid = coils.points.mean(axis=1)
    r = np.hypot(centroid[:, 0], centroid[:, 1])
    phi = np.arctan2(centroid[:, 1], centroid[:, 0])
    z = centroid[:, 2]
    scale = float(np.max(r))
    period = 2.0 * np.pi / nfp
    cur = np.abs(coils.currents)
    cur_scale = float(np.max(cur)) or 1.0

    def same_angle(a, b):
        d = np.mod(a - b + 0.5 * period, period) - 0.5 * period
        return np.abs(d) < rtol

    class_of = np.full(coils.n_coils, -1)
    found: list[int] = []
    for i in range(coils.n_coils):
        if class_of[i] >= 0:
            continue
        c = len(found)
        found.append(i)
        same = (np.abs(r - r[i]) < rtol * scale) \
            & (np.abs(coils.lengths - coils.lengths[i]) < rtol * coils.lengths[i]) \
            & (np.abs(cur - cur[i]) < rtol * cur_scale)
        rotated = same & same_angle(phi, phi[i]) & (np.abs(z - z[i]) < rtol * scale)
        match = rotated
        if stellsym:
            match = match | (same & same_angle(phi, -phi[i]) & (np.abs(z + z[i]) < rtol * scale))
        class_of[match & (class_of < 0)] = c
    reps = np.array(found)
    mult = np.bincount(class_of, minlength=len(reps))
    return SymmetryClasses(representatives=reps, multiplicity=mult, class_of=class_of)


def min_coil_coil_distance(coils: SampledCoils, classes: SymmetryClasses | None = None) -> float:
    """Minimum distance between the centre filaments of two different coils [m]."""
    reps = range(coils.n_coils) if classes is None else classes.representatives
    return float(min(np.min(coil_distances(coils, int(i))) for i in reps))


# ----------------------------------------------------------- toroidal lofting
def toroidal_order(coils: SampledCoils) -> np.ndarray:
    """Coil indices sorted by the toroidal angle of the coil centroid."""
    centroid = coils.points.mean(axis=1)
    return np.argsort(np.arctan2(centroid[:, 1], centroid[:, 0]), kind="stable")


def poloidally_aligned(coils: SampledCoils, n: int | None = None) -> np.ndarray:
    """Points of every coil at the same poloidal angles, shape ``[ncoil, n, 3]``.

    The poloidal angle of a point is measured around the coil's own centroid
    in the (R, z) plane. Index ``j`` of every coil then lies at the same
    poloidal angle, increasing in the same sense. If the angle is not
    monotonic along a coil, that coil is aligned by arc length instead,
    starting from its outermost point.
    """
    n = n or coils.n_points
    target = 2.0 * np.pi * np.arange(n) / n
    out = np.empty((coils.n_coils, n, 3))
    for i in range(coils.n_coils):
        pts = coils.points[i]
        r = np.hypot(pts[:, 0], pts[:, 1])
        r_c, z_c = r.mean(), pts[:, 2].mean()
        theta = np.unwrap(np.arctan2(pts[:, 2] - z_c, r - r_c))
        if theta[-1] < theta[0]:  # make the angle increase along the index
            pts, r = pts[::-1], r[::-1]
            theta = np.unwrap(np.arctan2(pts[:, 2] - z_c, r - r_c))
        if np.all(np.diff(theta) > 0) and theta[-1] - theta[0] < 2.0 * np.pi:
            ext_theta = np.concatenate([theta - 2.0 * np.pi, theta, theta + 2.0 * np.pi])
            ext_pts = np.concatenate([pts, pts, pts])
            base = np.mod(target - theta[0], 2.0 * np.pi) + theta[0]
            out[i] = np.column_stack([np.interp(base, ext_theta, ext_pts[:, k]) for k in range(3)])
        else:
            start = int(np.argmax(np.hypot(pts[:, 0], pts[:, 1])))
            rolled = np.roll(pts, -start, axis=0)
            idx = np.linspace(0, len(rolled), n, endpoint=False)
            lo = np.floor(idx).astype(int)
            w = (idx - lo)[:, None]
            out[i] = (1 - w) * rolled[lo % len(rolled)] + w * rolled[(lo + 1) % len(rolled)]
    return out


def coil_surface_area(coils: SampledCoils) -> float:
    """Mean area of the surfaces of revolution generated by the coils [m^2].

    Each coil, rotated about the vertical axis, sweeps a toroidal surface of
    area ``2 pi * integral of R dl``. The mean over the coils equals
    ``2 pi * coil_rmajor * mean coil length``. Equivalently, every coil length
    element is given a toroidal width of ``2 pi R / N``.
    """
    r = np.hypot(coils.points[..., 0], coils.points[..., 1])
    return float(2.0 * np.pi * np.mean(r.sum(axis=1) * coils.ds))


def lofted_surface_area(coils: SampledCoils) -> float:
    """Area of the surface spanned by straight lines between neighbouring coils [m^2].

    Neighbouring coils (in toroidal order) are joined between points at the
    same poloidal angle; the strips are triangulated and their areas summed.
    For non-planar coils this is smaller than :func:`coil_surface_area`,
    because the coils are not perpendicular to the toroidal direction.
    """
    aligned = poloidally_aligned(coils)[toroidal_order(coils)]
    nxt_coil = np.roll(aligned, -1, axis=0)
    a = aligned
    b = np.roll(aligned, -1, axis=1)
    c = nxt_coil
    d = np.roll(nxt_coil, -1, axis=1)
    tri1 = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=-1)
    tri2 = 0.5 * np.linalg.norm(np.cross(d - b, c - b), axis=-1)
    return float(np.sum(tri1 + tri2))


def port_width(coils: SampledCoils) -> float:
    """Largest centre-to-centre distance between neighbouring coils [m].

    Measured between the centroids of coils that are adjacent in toroidal
    order. For evenly spaced coils this is the coil spacing at the mean coil
    radius. The widths of the coils themselves are not subtracted; PROCESS
    applies factors of 0.4 and 0.8 to obtain port sizes.
    """
    centroid = coils.points.mean(axis=1)[toroidal_order(coils)]
    if len(centroid) < 2:
        return 0.0
    return float(np.linalg.norm(np.roll(centroid, -1, axis=0) - centroid, axis=1).max())


# ------------------------------------------------------------- plasma to coil
def _coils_near_period(coils: SampledCoils, nfp: int) -> np.ndarray:
    """Indices of coils that can be closest to the first field period."""
    if nfp == 1:
        return np.arange(coils.n_coils)
    period = 2.0 * np.pi / nfp
    centroid = coils.points.mean(axis=1)
    phi = np.mod(np.arctan2(centroid[:, 1], centroid[:, 0]) + 0.5 * period, 2.0 * np.pi)
    keep = phi < 2.0 * period  # the period itself plus half a period on each side
    return np.nonzero(keep)[0] if np.any(keep) else np.arange(coils.n_coils)


def plasma_coil_distance(eq: Equilibrium, coils: SampledCoils, s: float = 1.0,
                         n_theta: int = 96, n_phi: int = 96) -> float:
    """Minimum distance between the flux surface at ``s`` and any coil filament [m]."""
    require(0.0 < s <= 1.0, f"s must lie in (0, 1], got {s!r}")
    require(n_theta >= 8 and n_phi >= 8, "the surface needs at least 8 points per angle")
    theta, phi = angle_grid(eq, n_theta, n_phi)
    surface = eq.surface_xyz(s, theta, phi).reshape(-1, 3)
    idx = _coils_near_period(coils, eq.nfp)
    start = coils.points[idx]
    end = np.roll(start, -1, axis=1)
    dist, _ = point_segment_distance(surface, start.reshape(-1, 3), end.reshape(-1, 3))
    return float(dist.min())


def plasma_coil_distance_derivative(eq: Equilibrium, coils: SampledCoils,
                                    s_values=(1.0, 0.95, 0.9, 0.85),
                                    n_theta: int = 96, n_phi: int = 96):
    """Change of the minimum plasma-coil distance with the plasma minor radius.

    The distance is evaluated for several flux surfaces near the boundary and
    fitted linearly against the effective minor radius of each surface.
    Returns ``(slope, minor radii, distances)``. The slope is -1 when all
    surfaces are concentric circles.
    """
    require(len(set(s_values)) >= 2, "the fit needs at least two distinct surfaces")
    radii = np.array([surface_geometry(eq, s).minor_radius for s in s_values])
    dists = np.array([plasma_coil_distance(eq, coils, s, n_theta, n_phi) for s in s_values])
    slope = float(np.polyfit(radii, dists, 1)[0])
    return slope, radii, dists
