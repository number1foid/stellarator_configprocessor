"""Closed space curves: periodic cubic spline, resampling, curvature."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class CurveError(ValueError):
    """A curve cannot be built from the given points."""


def clean_closed_points(points: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Return the points of a closed curve without repeated points.

    Removes consecutive duplicates and a final point that repeats the first.
    ``tol`` is relative to the largest coordinate extent.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise CurveError(f"points must have shape [n, 3], got {pts.shape}")
    scale = float(np.max(np.ptp(pts, axis=0))) if len(pts) else 0.0
    if scale == 0.0:
        raise CurveError("curve has no extent")
    eps = tol * scale
    keep = np.ones(len(pts), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(pts, axis=0), axis=1) > eps
    pts = pts[keep]
    if len(pts) > 1 and np.linalg.norm(pts[-1] - pts[0]) <= eps:
        pts = pts[:-1]
    if len(pts) < 3:
        raise CurveError("a closed curve needs at least 3 distinct points")
    return pts


class PeriodicSpline:
    """Periodic cubic spline through the points of a closed curve.

    The parameter is the cumulative chord length, so it runs from 0 to
    ``period``. The curve is twice continuously differentiable.
    """

    def __init__(self, points: np.ndarray):
        pts = clean_closed_points(points)
        n = len(pts)
        nxt = np.roll(pts, -1, axis=0)
        h = np.linalg.norm(nxt - pts, axis=1)  # h[i] = t[i+1] - t[i], cyclic
        self._pts = pts
        self._h = h
        self._t = np.concatenate([[0.0], np.cumsum(h)])  # n + 1 knots
        self.period = float(self._t[-1])

        hp = np.roll(h, 1)  # h[i-1]
        a = np.zeros((n, n))
        idx = np.arange(n)
        a[idx, idx] = 2.0 * (hp + h)
        a[idx, (idx - 1) % n] += hp
        a[idx, (idx + 1) % n] += h
        slope = (nxt - pts) / h[:, None]
        rhs = 6.0 * (slope - np.roll(slope, 1, axis=0))
        self._m = np.linalg.solve(a, rhs)  # second derivatives at the knots

    def __call__(self, t: np.ndarray, derivative: int = 0) -> np.ndarray:
        """Evaluate the position or a derivative (1 or 2) at parameters ``t``."""
        t = np.mod(np.asarray(t, dtype=float), self.period)
        n = len(self._pts)
        i = np.clip(np.searchsorted(self._t, t, side="right") - 1, 0, n - 1)
        j = (i + 1) % n
        h = self._h[i][:, None]
        a = ((self._t[i + 1] - t) / self._h[i])[:, None]
        b = 1.0 - a
        y0, y1 = self._pts[i], self._pts[j]
        m0, m1 = self._m[i], self._m[j]
        if derivative == 0:
            return a * y0 + b * y1 + ((a**3 - a) * m0 + (b**3 - b) * m1) * h**2 / 6.0
        if derivative == 1:
            return (y1 - y0) / h - (3 * a**2 - 1) / 6.0 * h * m0 + (3 * b**2 - 1) / 6.0 * h * m1
        if derivative == 2:
            return a * m0 + b * m1
        raise ValueError("derivative must be 0, 1 or 2")


@dataclass(frozen=True)
class SampledCurve:
    """A closed curve sampled at points uniformly spaced in arc length."""

    points: np.ndarray  # [n, 3]
    tangents: np.ndarray  # [n, 3] unit tangents
    curvature: np.ndarray  # [n] in 1/m
    length: float

    @property
    def ds(self) -> float:
        """Arc length represented by each sample."""
        return self.length / len(self.points)


def sample_closed_curve(points: np.ndarray, n: int, oversample: int = 16) -> SampledCurve:
    """Fit a periodic spline and sample it at ``n`` points uniform in arc length."""
    if n < 8:
        raise CurveError("need at least 8 samples per curve")
    spline = PeriodicSpline(points)
    n_dense = max(oversample * len(clean_closed_points(points)), 8 * n)
    t_dense = np.linspace(0.0, spline.period, n_dense + 1)
    speed = np.linalg.norm(spline(t_dense[:-1] + 0.5 * np.diff(t_dense), 1), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(speed * np.diff(t_dense))])
    length = float(arc[-1])
    t = np.interp(np.arange(n) * length / n, arc, t_dense)
    d1 = spline(t, 1)
    d2 = spline(t, 2)
    speed = np.linalg.norm(d1, axis=1)
    curvature = np.linalg.norm(np.cross(d1, d2), axis=1) / speed**3
    return SampledCurve(
        points=spline(t), tangents=d1 / speed[:, None], curvature=curvature, length=length
    )
