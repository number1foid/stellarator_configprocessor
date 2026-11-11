import numpy as np
import pytest

from helpers import W7X_COILS, W7X_WOUT, circular_tf_coils, ideal_toroid_field
from stellarator_configprocessor.constants import MU0
from stellarator_configprocessor.geometry import coils as cg
from stellarator_configprocessor.io.coils import Coil, CoilSet, read_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout
from stellarator_configprocessor.winding_pack import (
    WindingPack, bar_kernel, coil_frames, filament_kernel, fit_peak_field, process_bmax,
    winding_pack_field)


def test_winding_pack_properties():
    wp = WindingPack(0.6, 0.5)
    assert wp.area == pytest.approx(0.3) and wp.ratio == pytest.approx(1.2)
    assert wp.scaled(2).area == pytest.approx(1.2) and wp.scaled(2).ratio == pytest.approx(1.2)
    with pytest.raises(ValueError):
        WindingPack(0.0, 1.0)


def brute_force_bar(u, v, a, b, n=600):
    xs = a * (-0.5 + (np.arange(n) + 0.5) / n)
    ys = b * (-0.5 + (np.arange(n) + 0.5) / n)
    du = u - xs[:, None]
    dv = v - ys[None, :]
    rho2 = du**2 + dv**2
    cell = a * b / n**2
    return np.sum(du / rho2) * cell, np.sum(dv / rho2) * cell


@pytest.mark.parametrize("u, v", [(0.0, 0.0), (0.2, 0.0), (0.2, 0.15), (0.0, 0.15),
                                  (0.07, -0.04), (0.5, 0.3), (-0.2, 0.0)])
def test_bar_kernel_matches_numerical_integration(u, v):
    a, b = 0.4, 0.3
    f_u, f_v = bar_kernel(u, v, a, b)
    ref_u, ref_v = brute_force_bar(u, v, a, b)
    scale = max(abs(ref_u), abs(ref_v), 1e-3)
    assert float(f_u) == pytest.approx(ref_u, abs=3e-3 * scale + 1e-6)
    assert float(f_v) == pytest.approx(ref_v, abs=3e-3 * scale + 1e-6)


def test_bar_kernel_far_field_is_a_line_current():
    f_u, f_v = bar_kernel(30.0, 40.0, 0.4, 0.3)
    area, rho2 = 0.12, 2500.0
    assert float(f_u) == pytest.approx(area * 30.0 / rho2, rel=1e-4)
    assert float(f_v) == pytest.approx(area * 40.0 / rho2, rel=1e-4)


def test_filament_kernel_of_one_filament():
    f_u, f_v = filament_kernel(np.array([3.0]), np.array([4.0]), [0.0], [0.0], 2.0)
    np.testing.assert_allclose([f_u[0], f_v[0]], [2 * 3 / 25, 2 * 4 / 25])


def test_frames_are_orthonormal_and_oriented():
    sc = circular_tf_coils(n_coils=6, r0=10.0, a_coil=3.0).sampled(64)
    n_r, n_t = coil_frames(sc)
    tan = sc.tangents
    for a, b in [(n_r, n_t), (n_r, tan), (n_t, tan)]:
        np.testing.assert_allclose(np.sum(a * b, axis=-1), 0, atol=1e-9)
    np.testing.assert_allclose(np.linalg.norm(n_r, axis=-1), 1)
    np.testing.assert_allclose(np.linalg.norm(n_t, axis=-1), 1)
    outward = sc.points - sc.points.mean(axis=1, keepdims=True)
    np.testing.assert_allclose(np.sum(n_r * outward, axis=-1), 3.0, rtol=1e-3)
    phi = np.arctan2(sc.points[..., 1], sc.points[..., 0])
    e_phi = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=-1)
    np.testing.assert_allclose(np.abs(np.sum(n_t * e_phi, axis=-1)), 1, atol=1e-6)


def test_toroidally_running_coil_is_rejected():
    t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    ring = np.column_stack([5 * np.cos(t), 5 * np.sin(t), np.zeros_like(t)])
    with pytest.raises(ValueError, match="modular"):
        coil_frames(CoilSet(coils=(Coil(points=ring, current=1.0),)).sampled(64))


@pytest.fixture(scope="module")
def big_loop():
    """A loop so large that its conductor is locally a straight bar."""
    radius = 50.0
    t = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
    pts = np.column_stack([110 + radius * np.cos(t), np.zeros_like(t), radius * np.sin(t)])
    return CoilSet(coils=(Coil(points=pts, current=1.0e6),)).sampled(2000)


def test_local_field_matches_straight_bar(big_loop):
    wp = WindingPack(0.4, 0.3)
    stations = (np.array([0, 500]),)
    pf = winding_pack_field(big_loop, wp, [0], n_sub=(4, 4), stations=stations)
    ku = wp.radial * (-0.5 + np.arange(5) / 4)
    kv = wp.toroidal * (-0.5 + np.arange(5) / 4)
    uu, vv = np.meshgrid(ku, kv, indexing="ij")
    f_u, f_v = bar_kernel(uu, vv, wp.radial, wp.toroidal)
    exact = MU0 * 1.0e6 / wp.area / (2 * np.pi) * np.hypot(f_u, f_v)
    got = np.linalg.norm(pf.b[0], axis=-1)
    # The remaining difference is the field of the rest of the loop.
    assert np.max(np.abs(got - exact[None])) < 0.02 * exact.max()
    assert pf.peak() == pytest.approx(exact.max(), rel=0.02)


def test_local_correction_removes_discretisation_dependence(big_loop):
    wp = WindingPack(0.4, 0.3)
    stations = (np.array([0]),)
    coarse = winding_pack_field(big_loop, wp, [0], n_sub=(2, 2), stations=stations).peak()
    fine = winding_pack_field(big_loop, wp, [0], n_sub=(6, 6), stations=stations).peak()
    assert coarse == pytest.approx(fine, rel=1e-3)
    raw = winding_pack_field(big_loop, wp, [0], n_sub=(2, 2), stations=stations,
                             local_correction=False).peak()
    assert abs(raw / fine - 1) > 0.02  # without the correction it is not converged


def test_weights_and_mean_field_of_uniform_field(big_loop):
    pf = winding_pack_field(big_loop, WindingPack(0.4, 0.3), [0], n_sub=(3, 5),
                            stations=(np.array([0]),))
    assert pf.weights.shape == (4, 6) and pf.weights.sum() == pytest.approx(1.0)
    assert pf.mean_field()[0].shape == (1, 3)


def test_station_subset_matches_full_evaluation():
    sc = circular_tf_coils(n_coils=8, r0=10.0, a_coil=3.0).sampled(48)
    wp = WindingPack(0.5, 0.4)
    full = winding_pack_field(sc, wp, [0, 3], n_sub=(2, 2))
    sub = winding_pack_field(sc, wp, [0, 3], n_sub=(2, 2),
                             stations=(np.array([5, 20]), np.array([0])))
    np.testing.assert_allclose(sub.b[0], full.b[0][[5, 20]], rtol=1e-12)
    np.testing.assert_allclose(sub.b[1], full.b[1][[0]], rtol=1e-12)


def test_tf_coil_peak_field_is_on_the_inboard_side():
    n, r0, ac, cur = 16, 10.0, 4.0, 5.0e6
    sc = circular_tf_coils(n_coils=n, r0=r0, a_coil=ac, current=cur).sampled(96)
    wp = WindingPack(0.8, 0.6)
    pf = winding_pack_field(sc, wp, [0], n_sub=(4, 4))
    strength = pf.station_peaks()[0]
    station = int(np.argmax(strength))
    r_peak = np.hypot(*sc.points[0][station][:2])
    assert r_peak < r0 - 0.9 * ac  # inboard leg
    # Above the ideal toroid field at the inboard filament radius, by the ripple
    # and self-field enhancement, but of the same order.
    ideal = ideal_toroid_field(n, cur, r0 - ac)
    assert 1.0 < pf.peak() / ideal < 1.6


def test_fit_recovers_process_parameters():
    widths = np.array([0.3, 0.45, 0.6, 0.75, 0.9])
    args = dict(current_ma=13.0, n_coils=50, r_major=22.0, r_minor=4.7)
    peaks = process_bmax(a1=0.7, a2=0.03, radial_width=widths, **args)
    fit = fit_peak_field(widths, peaks, **args)
    assert fit.a1 == pytest.approx(0.7) and fit.a2 == pytest.approx(0.03)
    assert fit.max_relative_error < 1e-12
    with pytest.raises(ValueError):
        fit_peak_field(widths[:1], peaks[:1], **args)


def test_process_formula_reproduces_example_file():
    # Values from the example file shipped with PROCESS (Helias 5).
    radial = np.sqrt(0.53333333 * 1.2)
    b = process_bmax(12.98310033, 50, 22.22954927, 4.69015768, 0.0254274, 0.05608523, radial)
    assert float(b) == pytest.approx(11.49413032, rel=0.03)


@pytest.mark.w7x
@pytest.mark.slow
def test_w7x_peak_field_matches_process_builtin():
    eq = read_wout(W7X_WOUT)
    sc = read_makegrid(W7X_COILS).sampled(160)
    classes = cg.symmetry_classes(sc, eq.nfp)
    pf = winding_pack_field(sc, WindingPack(0.18, 0.15), classes.representatives, n_sub=(3, 3))
    # Built-in W7X-50: 6.3 T on the conductor at 3.0 T on axis. The real coils
    # at 15 kA x 108 turns give 2.785 T on axis.
    assert pf.peak() * 3.0 / 2.7847 == pytest.approx(6.3, rel=0.03)
