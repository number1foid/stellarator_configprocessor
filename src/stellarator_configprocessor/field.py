"""Biot-Savart field of straight current segments.

A closed filament is represented as a polygon. The field of each straight
segment is evaluated with the exact closed-form expression, so the only
discretisation error is the polygon approximation of the curve.
"""

from __future__ import annotations

import numpy as np

from .constants import MU0

#: Upper bound on point-segment pairs evaluated in one vectorised block.
_MAX_PAIRS = 1_500_000


def biot_savart(points: np.ndarray, start: np.ndarray, end: np.ndarray,
                current: np.ndarray) -> np.ndarray:
    """Magnetic field [T] at ``points`` from straight segments.

    Parameters
    ----------
    points : array ``[npts, 3]``
        Evaluation points [m].
    start, end : arrays ``[nseg, 3]``
        Segment end points [m]. Current flows from ``start`` to ``end``.
    current : array ``[nseg]`` or scalar
        Current in each segment [A].

    A point lying exactly on a segment receives no contribution from that
    segment (the field of a line current is undefined there).
    """
    points = np.ascontiguousarray(points, dtype=float).reshape(-1, 3)
    start = np.ascontiguousarray(start, dtype=float).reshape(-1, 3)
    end = np.ascontiguousarray(end, dtype=float).reshape(-1, 3)
    current = np.broadcast_to(np.asarray(current, dtype=float), (len(start),))
    out = np.empty_like(points)
    step = max(1, _MAX_PAIRS // max(1, len(start)))
    for i in range(0, len(points), step):
        out[i:i + step] = _block(points[i:i + step], start, end, current)
    return out


def _block(x, p1, p2, cur):
    r1 = x[:, None, :] - p1[None, :, :]
    r2 = x[:, None, :] - p2[None, :, :]
    n1 = np.sqrt(np.einsum("psk,psk->ps", r1, r1))
    n2 = np.sqrt(np.einsum("psk,psk->ps", r2, r2))
    dot = np.einsum("psk,psk->ps", r1, r2)
    denom = n1 * n2 * (n1 * n2 + dot)
    # denom -> 0 only for points on the segment itself.
    tiny = 1e-14 * (n1 * n2) ** 2 + 1e-300
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(denom > tiny, cur[None, :] * (n1 + n2) / denom, 0.0)
    b = np.empty_like(x)
    b[:, 0] = np.einsum("ps,ps->p", r1[:, :, 1] * r2[:, :, 2] - r1[:, :, 2] * r2[:, :, 1], fac)
    b[:, 1] = np.einsum("ps,ps->p", r1[:, :, 2] * r2[:, :, 0] - r1[:, :, 0] * r2[:, :, 2], fac)
    b[:, 2] = np.einsum("ps,ps->p", r1[:, :, 0] * r2[:, :, 1] - r1[:, :, 1] * r2[:, :, 0], fac)
    return b * (MU0 / (4.0 * np.pi))
