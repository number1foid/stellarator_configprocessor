"""Biot-Savart field of straight current segments.

A closed filament is represented as a polygon. The field of each straight
segment is evaluated with the exact closed-form expression, so the only
discretisation error is the polygon approximation of the curve.

The kernel is the single place where the field is computed, so it can be
replaced by a compiled implementation without touching the callers.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .constants import MU0

#: Point-segment pairs evaluated per vectorised block. Small blocks keep the
#: temporaries in the CPU cache, which is about three times faster.
_MAX_PAIRS = 40_000

#: Worker threads. NumPy releases the interpreter lock inside array operations.
_THREADS = max(1, min(8, (os.cpu_count() or 2) // 2))


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
    current = np.ascontiguousarray(np.broadcast_to(np.asarray(current, dtype=float),
                                                   (len(start),)))
    out = np.zeros_like(points)
    if len(points) == 0 or len(start) == 0:
        return out
    seg = tuple(np.ascontiguousarray(start[:, k]) for k in range(3)) \
        + tuple(np.ascontiguousarray(end[:, k]) for k in range(3))
    step = max(1, _MAX_PAIRS // len(start))
    blocks = [slice(i, i + step) for i in range(0, len(points), step)]

    def work(sl):
        out[sl] = _block(points[sl], seg, current)

    if _THREADS > 1 and len(blocks) > 1:
        with ThreadPoolExecutor(_THREADS) as pool:
            list(pool.map(work, blocks))
    else:
        for sl in blocks:
            work(sl)
    return out


def _block(x, seg, cur):
    ax, ay, az, bx, by, bz = seg
    px, py, pz = x[:, 0, None], x[:, 1, None], x[:, 2, None]
    x1, y1, z1 = px - ax, py - ay, pz - az
    x2, y2, z2 = px - bx, py - by, pz - bz
    n1 = np.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    n2 = np.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    n12 = n1 * n2
    denom = n12 * (n12 + x1 * x2 + y1 * y2 + z1 * z2)
    # denom -> 0 only for points on the segment itself.
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(denom > 1e-14 * n12 * n12 + 1e-300, (n1 + n2) / denom, 0.0)
    fac *= cur
    b = np.empty_like(x)
    b[:, 0] = np.einsum("ps,ps->p", y1 * z2 - z1 * y2, fac)
    b[:, 1] = np.einsum("ps,ps->p", z1 * x2 - x1 * z2, fac)
    b[:, 2] = np.einsum("ps,ps->p", x1 * y2 - y1 * x2, fac)
    return b * (MU0 / (4.0 * np.pi))
