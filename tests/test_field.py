import numpy as np

from stellarator_configprocessor.constants import MU0
from stellarator_configprocessor.field import biot_savart


def loop(radius, n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pts = np.column_stack([radius * np.cos(t), radius * np.sin(t), np.zeros(n)])
    return pts, np.roll(pts, -1, axis=0)


def test_circular_loop_on_axis():
    radius, current = 1.5, 2.0e5
    start, end = loop(radius, 400)
    z = np.array([0.0, 0.4, 2.0])
    pts = np.column_stack([np.zeros(3), np.zeros(3), z])
    b = biot_savart(pts, start, end, current)
    exact = MU0 * current * radius**2 / (2 * (radius**2 + z**2) ** 1.5)
    np.testing.assert_allclose(b[:, 2], exact, rtol=1e-4)
    np.testing.assert_allclose(b[:, :2], 0, atol=1e-12)


def test_long_straight_wire():
    current, d = 1.0e4, 0.3
    start = np.array([[0.0, 0.0, -1.0e4]])
    end = np.array([[0.0, 0.0, 1.0e4]])
    b = biot_savart(np.array([[d, 0.0, 0.0]]), start, end, current)
    np.testing.assert_allclose(b[0], [0, MU0 * current / (2 * np.pi * d), 0], rtol=1e-6, atol=1e-15)


def test_finite_segment_matches_formula():
    # Field at distance d from the middle of a segment of half-length h.
    current, d, h = 3.0e3, 0.5, 0.7
    b = biot_savart(np.array([[d, 0, 0]]), np.array([[0, 0, -h]]), np.array([[0, 0, h]]), current)
    exact = MU0 * current / (4 * np.pi * d) * 2 * h / np.hypot(d, h)
    np.testing.assert_allclose(b[0, 1], exact, rtol=1e-12)


def test_point_on_segment_gives_no_contribution():
    b = biot_savart(np.array([[0, 0, 0.2]]), np.array([[0, 0, 0.0]]), np.array([[0, 0, 1.0]]), 1.0)
    assert np.all(np.isfinite(b)) and np.all(b == 0)


def test_blocking_does_not_change_result(monkeypatch):
    import stellarator_configprocessor.field as field
    rng = np.random.default_rng(0)
    start, end = loop(1.0, 50)
    pts = rng.uniform(-3, 3, (40, 3))
    full = biot_savart(pts, start, end, 1.0)
    monkeypatch.setattr(field, "_MAX_PAIRS", 120)
    np.testing.assert_allclose(biot_savart(pts, start, end, 1.0), full, rtol=1e-13)
