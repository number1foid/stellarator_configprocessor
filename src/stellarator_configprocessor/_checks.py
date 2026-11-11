"""Argument checks shared by the numerical routines.

The checks raise ``ValueError`` rather than use ``assert``, so they stay
active when Python runs with ``-O``.
"""

from __future__ import annotations

import numpy as np


def require(condition: bool, message: str) -> None:
    """Raise ``ValueError(message)`` unless ``condition`` holds."""
    if not condition:
        raise ValueError(message)


def points_array(value, name: str) -> np.ndarray:
    """``value`` as a finite float array whose last axis has length 3."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError(f"{name} must have shape [..., 3], got {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def coil_index(i: int, n_coils: int) -> int:
    """``i`` as a valid coil index."""
    index = int(i)
    if not 0 <= index < n_coils:
        raise ValueError(f"coil index {i} is outside 0..{n_coils - 1}")
    return index


def positive(value: float, name: str) -> float:
    """``value`` as a finite float greater than zero."""
    number = float(value)
    if not (np.isfinite(number) and number > 0.0):
        raise ValueError(f"{name} must be a positive number, got {value!r}")
    return number
