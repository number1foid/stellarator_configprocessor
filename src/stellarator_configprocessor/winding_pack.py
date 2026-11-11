"""Magnetic field inside coils with a finite rectangular winding pack.

The field of an ideal filament diverges at the filament, so the peak field on
the conductor and the forces need a finite cross-section. The model follows
Lion et al. 2021, section 3.7: a rectangular winding pack (WP) carrying a
uniform current density, with one side normal along the toroidal direction.

Numerical scheme
----------------
* The WP of a coil is replaced by ``n_r x n_t`` filaments at the centres of
  equal cells, each carrying an equal share of the current.
* The field is evaluated on the ``(n_r + 1) x (n_t + 1)`` cell corners, which
  are the points farthest from the filaments and include the WP surface.
* **Local correction.** Close to the conductor the difference between a
  uniform current density and discrete filaments is dominated by the
  two-dimensional field of a straight bar. That difference is known in closed
  form and is added, so that only the smooth, curvature-related part of the
  field relies on the filament discretisation.
* Coils far from the evaluated coil are represented by their centre filament.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MU0
from .field import biot_savart
from .geometry.coils import coil_distances
from .io.coils import SampledCoils


@dataclass(frozen=True)
class WindingPack:
    """Rectangular winding pack cross-section."""

    radial: float  # width in the coil's in-plane normal direction [m]
    toroidal: float  # width in the toroidal direction [m]

    def __post_init__(self):
        if not (self.radial > 0 and self.toroidal > 0):
            raise ValueError("winding pack dimensions must be positive")

    @property
    def area(self) -> float:
        return self.radial * self.toroidal

    @property
    def ratio(self) -> float:
        """Radial width divided by toroidal width."""
        return self.radial / self.toroidal

    @property
    def diagonal(self) -> float:
        return float(np.hypot(self.radial, self.toroidal))

    def scaled(self, factor: float) -> "WindingPack":
        return WindingPack(self.radial * factor, self.toroidal * factor)


def coil_frames(coils: SampledCoils) -> tuple[np.ndarray, np.ndarray]:
    """Unit normals of the winding pack at every sample point.

    Returns ``(n_r, n_t)``, both ``[ncoil, n, 3]``:

    * ``n_t``: the cylindrical toroidal direction made perpendicular to the
      filament tangent;
    * ``n_r``: perpendicular to the tangent and to ``n_t``, pointing away
      from the coil centroid (outwards from the plasma).
    """
    pts, tan = coils.points, coils.tangents
    phi = np.arctan2(pts[..., 1], pts[..., 0])
    e_phi = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=-1)
    n_t = e_phi - np.sum(e_phi * tan, axis=-1, keepdims=True) * tan
    norm = np.linalg.norm(n_t, axis=-1, keepdims=True)
    if np.any(norm < 0.2):
        raise ValueError(
            "a coil runs nearly parallel to the toroidal direction; the winding pack "
            "orientation is undefined. Only modular (poloidally closed) coils are supported.")
    n_t = n_t / norm
    n_r = np.cross(tan, n_t)
    outward = pts - pts.mean(axis=1, keepdims=True)
    sign = np.sign(np.sum(n_r * outward, axis=(1, 2)))
    sign[sign == 0] = 1.0
    return n_r * sign[:, None, None], n_t


# ------------------------------------------------- two-dimensional bar fields
def bar_kernel(u, v, width_u: float, width_v: float):
    """Field integrals of a uniform straight rectangular bar.

    For a bar occupying ``|u'| <= width_u / 2``, ``|v'| <= width_v / 2``,
    returns the integrals over the cross-section of ``(u - u') / rho^2`` and
    ``(v - v') / rho^2`` with ``rho^2 = (u - u')^2 + (v - v')^2``.
    The field is ``mu0 J / (2 pi) * (t x n_u * F_u + t x n_v * F_v)``.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    def g(x, y):
        # Antiderivative of x / (x^2 + y^2) in x and y, without the -y term,
        # which cancels in the double difference.
        r2 = x * x + y * y
        with np.errstate(divide="ignore", invalid="ignore"):
            log_term = np.where(r2 > 0, 0.5 * y * np.log(np.where(r2 > 0, r2, 1.0)), 0.0)
            atan_term = np.where(x != 0, x * np.arctan(y / np.where(x != 0, x, 1.0)), 0.0)
        return log_term + atan_term

    x1, x2 = u - 0.5 * width_u, u + 0.5 * width_u
    y1, y2 = v - 0.5 * width_v, v + 0.5 * width_v
    f_u = g(x2, y2) - g(x1, y2) - g(x2, y1) + g(x1, y1)
    f_v = g(y2, x2) - g(y1, x2) - g(y2, x1) + g(y1, x1)
    return f_u, f_v


def filament_kernel(u, v, fil_u, fil_v, cell_area: float):
    """The same integrals for a set of straight filaments, each of weight ``cell_area``."""
    du = np.asarray(u, dtype=float)[..., None] - np.asarray(fil_u, dtype=float)
    dv = np.asarray(v, dtype=float)[..., None] - np.asarray(fil_v, dtype=float)
    rho2 = du * du + dv * dv
    return cell_area * np.sum(du / rho2, axis=-1), cell_area * np.sum(dv / rho2, axis=-1)


# -------------------------------------------------------------- field in WP
@dataclass(frozen=True)
class PackField:
    """Field on the cell corners of the winding pack of selected coils."""

    targets: np.ndarray  # [T] coil index of each evaluated coil
    stations: tuple  # T arrays: sample indices evaluated along each coil
    points: tuple  # T arrays [nstation, n_r + 1, n_t + 1, 3]: evaluation points [m]
    b: tuple  # T arrays [nstation, n_r + 1, n_t + 1, 3] in tesla
    weights: np.ndarray  # [n_r + 1, n_t + 1] quadrature weights, sum = 1
    wp: WindingPack

    def peak(self) -> float:
        """Largest field strength found [T]."""
        return float(max(np.linalg.norm(b, axis=-1).max() for b in self.b))

    def station_peaks(self) -> tuple:
        """Largest field strength in the cross-section at each station."""
        return tuple(np.linalg.norm(b, axis=-1).max(axis=(1, 2)) for b in self.b)

    def mean_field(self) -> tuple:
        """Cross-section average of the field at each station, ``[nstation, 3]``."""
        return tuple(np.einsum("ab,sabk->sk", self.weights, b) for b in self.b)


def _offsets(width: float, n: int):
    centres = width * (-0.5 + (np.arange(n) + 0.5) / n)
    corners = width * (-0.5 + np.arange(n + 1) / n)
    return centres, corners


def winding_pack_field(
    coils: SampledCoils,
    wp: WindingPack,
    targets,
    frames: tuple[np.ndarray, np.ndarray] | None = None,
    n_sub: tuple[int, int] = (4, 4),
    near_factor: float = 3.0,
    stations: tuple | None = None,
    local_correction: bool = True,
) -> PackField:
    """Field inside the winding pack of the coils listed in ``targets``.

    Parameters
    ----------
    n_sub
        Number of filaments across the radial and toroidal width.
    near_factor
        Coils closer to the target than ``near_factor`` WP diagonals are
        modelled with the full filament bundle, the rest by one filament.
    stations
        Optional per-target arrays of sample indices to evaluate; all samples
        if omitted.
    """
    n_r, n_t = frames if frames is not None else coil_frames(coils)
    nr, nt = n_sub
    if nr < 1 or nt < 1:
        raise ValueError("n_sub entries must be at least 1")
    cu, ku = _offsets(wp.radial, nr)
    cv, kv = _offsets(wp.toroidal, nt)
    targets = np.atleast_1d(np.asarray(targets, dtype=int))
    all_stations = np.arange(coils.n_points)

    def bundle(k):
        """Segments of the filament bundle of coil k."""
        pts = (coils.points[k][None, None]
               + cu[:, None, None, None] * n_r[k][None, None]
               + cv[None, :, None, None] * n_t[k][None, None])  # [nr, nt, n, 3]
        start = pts.reshape(-1, 3)
        end = np.roll(pts, -1, axis=2).reshape(-1, 3)
        cur = np.full(len(start), coils.currents[k] / (nr * nt))
        return start, end, cur

    # Closed-form difference between a uniform bar and the filament bundle.
    uu, vv = np.meshgrid(ku, kv, indexing="ij")
    bar_u, bar_v = bar_kernel(uu, vv, wp.radial, wp.toroidal)
    fu, fv = np.meshgrid(cu, cv, indexing="ij")
    fil_u, fil_v = filament_kernel(uu, vv, fu.ravel(), fv.ravel(), wp.area / (nr * nt))
    delta_u, delta_v = bar_u - fil_u, bar_v - fil_v

    wr = np.ones(nr + 1)
    wr[[0, -1]] = 0.5
    wt = np.ones(nt + 1)
    wt[[0, -1]] = 0.5
    weights = np.outer(wr / nr, wt / nt)

    fields, used_stations, used_points = [], [], []
    for pos, i in enumerate(targets):
        st = all_stations if stations is None else np.asarray(stations[pos], dtype=int)
        used_stations.append(st)
        if len(st) == 0:
            fields.append(np.zeros((0, nr + 1, nt + 1, 3)))
            used_points.append(np.zeros((0, nr + 1, nt + 1, 3)))
            continue
        base = coils.points[i][st]
        points = (base[:, None, None, :]
                  + ku[None, :, None, None] * n_r[i][st][:, None, None, :]
                  + kv[None, None, :, None] * n_t[i][st][:, None, None, :])

        dist = coil_distances(coils, int(i))
        near = np.nonzero(dist < near_factor * wp.diagonal)[0]
        far = np.array([k for k in range(coils.n_coils) if k != i and k not in set(near)],
                       dtype=int)
        parts = [bundle(int(i))] + [bundle(int(k)) for k in near]
        start = np.concatenate([p[0] for p in parts])
        end = np.concatenate([p[1] for p in parts])
        cur = np.concatenate([p[2] for p in parts])
        b = biot_savart(points.reshape(-1, 3), start, end, cur).reshape(points.shape)
        if len(far):
            # The field of distant coils varies slowly over the cross-section:
            # evaluate it at the four outer corners and interpolate bilinearly.
            outer = points[:, [0, -1]][:, :, [0, -1]]
            b4 = biot_savart(outer.reshape(-1, 3), *coils.segments(far)).reshape(outer.shape)
            wu = ((ku - ku[0]) / (ku[-1] - ku[0]))[None, :, None, None]
            wv = ((kv - kv[0]) / (kv[-1] - kv[0]))[None, None, :, None]
            b = b + ((1 - wu) * (1 - wv) * b4[:, 0, 0][:, None, None, :]
                     + wu * (1 - wv) * b4[:, 1, 0][:, None, None, :]
                     + (1 - wu) * wv * b4[:, 0, 1][:, None, None, :]
                     + wu * wv * b4[:, 1, 1][:, None, None, :])

        if local_correction:
            tan = coils.tangents[i][st]
            t_x_nr = np.cross(tan, n_r[i][st])
            t_x_nt = np.cross(tan, n_t[i][st])
            density = coils.currents[i] / wp.area
            b = b + MU0 * density / (2.0 * np.pi) * (
                t_x_nr[:, None, None, :] * delta_u[None, :, :, None]
                + t_x_nt[:, None, None, :] * delta_v[None, :, :, None])
        fields.append(b)
        used_points.append(points)

    return PackField(targets=targets, stations=tuple(used_stations), points=tuple(used_points),
                     b=tuple(fields), weights=weights, wp=wp)


def pack_clearance(coils: SampledCoils, wp: WindingPack, targets,
                   frames: tuple[np.ndarray, np.ndarray] | None = None) -> float:
    """Smallest gap between the winding packs of different coils [m].

    For every corner of the winding pack of each target coil, the nearest
    cross-section of every other coil is located and the distance from the
    corner to that rectangle is measured. A value of zero or less means two
    winding packs touch or overlap. Because only corners are tested, the
    result is slightly optimistic for packs that meet edge to edge.
    """
    n_r, n_t = frames if frames is not None else coil_frames(coils)
    corners = np.array([[su * 0.5 * wp.radial, sv * 0.5 * wp.toroidal]
                        for su in (-1.0, 1.0) for sv in (-1.0, 1.0)])
    clearance = np.inf
    for i in np.atleast_1d(np.asarray(targets, dtype=int)):
        pts = (coils.points[i][:, None, :]
               + corners[None, :, 0, None] * n_r[i][:, None, :]
               + corners[None, :, 1, None] * n_t[i][:, None, :]).reshape(-1, 3)
        dist = coil_distances(coils, int(i))
        # Nearby coils, and always the nearest one.
        for k in np.nonzero(dist <= max(2.0 * wp.diagonal, 1.01 * dist.min()))[0]:
            rel = pts[:, None, :] - coils.points[k][None, :, :]
            nearest = np.argmin(np.einsum("psk,psk->ps", rel, rel), axis=1)
            d = pts - coils.points[k][nearest]
            du = np.abs(np.einsum("pk,pk->p", d, n_r[k][nearest])) - 0.5 * wp.radial
            dv = np.abs(np.einsum("pk,pk->p", d, n_t[k][nearest])) - 0.5 * wp.toroidal
            # Signed distance to the rectangle: negative inside.
            outside = np.hypot(np.maximum(du, 0.0), np.maximum(dv, 0.0))
            inside = np.minimum(np.maximum(du, dv), 0.0)
            clearance = min(clearance, float(np.min(outside + inside)))
    return clearance


# ------------------------------------------------------------- peak-field fit
@dataclass(frozen=True)
class PeakFieldFit:
    """Fit of the peak field against winding pack size, in PROCESS's form."""

    a1: float
    a2: float
    radial_widths: np.ndarray  # [nsize] radial WP width of each sample [m]
    peak_fields: np.ndarray  # [nsize] computed peak field [T]
    fitted: np.ndarray  # [nsize] value of the fit at the samples [T]

    @property
    def max_relative_error(self) -> float:
        return float(np.max(np.abs(self.fitted / self.peak_fields - 1.0)))


def process_bmax(current_ma: float, n_coils: int, r_major: float, r_minor: float,
                 a1: float, a2: float, radial_width):
    """PROCESS's peak-field formula (``bmax_from_awp``) [T]."""
    return (0.2 * current_ma * n_coils / (r_major - r_minor)
            * (a1 + a2 * r_major / np.asarray(radial_width, dtype=float)))


def fit_peak_field(radial_widths, peak_fields, current_ma: float, n_coils: int,
                   r_major: float, r_minor: float) -> PeakFieldFit:
    """Least-squares fit of ``a1`` and ``a2`` in :func:`process_bmax`."""
    radial_widths = np.asarray(radial_widths, dtype=float)
    peak_fields = np.asarray(peak_fields, dtype=float)
    if len(radial_widths) < 2:
        raise ValueError("at least two winding pack sizes are needed for the fit")
    scale = 0.2 * current_ma * n_coils / (r_major - r_minor)
    design = np.column_stack([np.ones_like(radial_widths), r_major / radial_widths])
    (a1, a2), *_ = np.linalg.lstsq(design, peak_fields / scale, rcond=None)
    fitted = process_bmax(current_ma, n_coils, r_major, r_minor, a1, a2, radial_widths)
    return PeakFieldFit(a1=float(a1), a2=float(a2), radial_widths=radial_widths,
                        peak_fields=peak_fields, fitted=fitted)
