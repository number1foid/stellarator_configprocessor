import numpy as np
import pytest

from helpers import W7X_COILS, W7X_WOUT, circular_tf_coils, write_torus_wout
from stellarator_configprocessor.geometry import coils as cg
from stellarator_configprocessor.io.coils import Coil, CoilSet, read_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout

N, R0, AC = 20, 10.0, 4.0


@pytest.fixture(scope="module")
def tf():
    return circular_tf_coils(n_coils=N, r0=R0, a_coil=AC, n_points=120).sampled(240)


def test_radii_length_height_bend_radius(tf):
    r_major, r_minor, z_mean = cg.coil_radii(tf)
    assert r_major == pytest.approx(R0, rel=1e-4)
    assert r_minor == pytest.approx(AC, rel=1e-4)
    assert z_mean == pytest.approx(0.0, abs=1e-9)
    assert cg.total_length(tf) == pytest.approx(N * 2 * np.pi * AC, rel=1e-4)
    assert cg.coil_height(tf) == pytest.approx(2 * AC, rel=1e-4)
    assert cg.min_bend_radius(tf) == pytest.approx(AC, rel=5e-3)


def test_point_segment_distance():
    start = np.array([[0.0, 0, 0], [5, 0, 0]])
    end = np.array([[1.0, 0, 0], [5, 0, 2]])
    pts = np.array([[0.5, 2, 0], [-1, 0, 0], [5, 1, 3], [3, 0, 0]])
    dist, idx = cg.point_segment_distance(pts, start, end)
    np.testing.assert_allclose(dist, [2.0, 1.0, np.sqrt(2), 2.0])
    np.testing.assert_array_equal(idx[:3], [0, 0, 1])


def test_min_coil_distance_is_inboard_gap(tf):
    expected = 2 * (R0 - AC) * np.sin(np.pi / N)
    assert cg.min_coil_coil_distance(tf) == pytest.approx(expected, rel=2e-3)
    classes = cg.symmetry_classes(tf, nfp=5)
    assert cg.min_coil_coil_distance(tf, classes) == pytest.approx(expected, rel=2e-3)


def test_symmetry_classes(tf):
    classes = cg.symmetry_classes(tf, nfp=5, stellsym=True)
    assert len(classes.representatives) == 2
    np.testing.assert_array_equal(classes.multiplicity, [10, 10])
    assert np.all(classes.class_of >= 0)
    rotations_only = cg.symmetry_classes(tf, nfp=5, stellsym=False)
    np.testing.assert_array_equal(rotations_only.multiplicity, [5, 5, 5, 5])


def test_unsymmetric_set_gives_one_class_per_coil():
    cs = circular_tf_coils(n_coils=4, r0=R0, a_coil=AC)
    coils = list(cs.coils)
    coils[1] = Coil(points=coils[1].points * 1.05, current=coils[1].current)
    classes = cg.symmetry_classes(CoilSet(coils=tuple(coils)).sampled(64), nfp=4)
    assert len(classes.representatives) == 2
    assert sorted(classes.multiplicity) == [1, 3]


def test_coil_surface_of_circular_coils_is_the_torus(tf):
    assert cg.coil_surface_area(tf) == pytest.approx(4 * np.pi**2 * R0 * AC, rel=1e-4)
    r_major, _, _ = cg.coil_radii(tf)
    assert cg.coil_surface_area(tf) == pytest.approx(
        2 * np.pi * r_major * tf.lengths.mean(), rel=1e-12)


def test_lofted_surface_is_polygonal_torus(tf):
    # Straight strips join neighbouring coils. At poloidal angle theta a strip
    # has chord width 2 (R0 + a cos theta) sin(pi/N), and the coil element makes
    # an angle with the chord whose sine is sqrt(1 - sin^2 theta sin^2(pi/N)).
    theta = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    half = np.pi / N
    integrand = (2 * (R0 + AC * np.cos(theta)) * np.sin(half)
                 * np.sqrt(1 - np.sin(theta) ** 2 * np.sin(half) ** 2))
    expected = N * AC * integrand.mean() * 2 * np.pi
    assert cg.lofted_surface_area(tf) == pytest.approx(expected, rel=1e-3)
    assert cg.lofted_surface_area(tf) < cg.coil_surface_area(tf)


def test_port_width_is_coil_spacing(tf):
    assert cg.port_width(tf) == pytest.approx(2 * R0 * np.sin(np.pi / N), rel=1e-4)


def test_alignment_is_independent_of_start_point_and_direction(tf):
    pts = tf.points[0]
    shuffled = np.roll(pts[::-1], 17, axis=0)
    other = CoilSet(coils=(Coil(points=shuffled, current=1.0),)).sampled(240)
    a = cg.poloidally_aligned(tf)[0]
    b = cg.poloidally_aligned(other)[0]
    np.testing.assert_allclose(a, b, atol=2e-3)


def test_plasma_coil_distance_and_derivative(tmp_path, tf):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=R0, a=2.0, ns=21))
    assert cg.plasma_coil_distance(eq, tf) == pytest.approx(AC - 2.0, rel=2e-3)
    slope, radii, dists = cg.plasma_coil_distance_derivative(eq, tf)
    assert slope == pytest.approx(-1.0, abs=5e-3)
    assert np.all(np.diff(radii) < 0) and np.all(np.diff(dists) > 0)


@pytest.mark.w7x
def test_w7x_geometry_is_close_to_process_builtin():
    """The real coils are close to the FOCUS-generated W7X-50 set in PROCESS."""
    eq = read_wout(W7X_WOUT)
    sc = read_makegrid(W7X_COILS).sampled(160)
    classes = cg.symmetry_classes(sc, eq.nfp)
    np.testing.assert_array_equal(classes.multiplicity, [10] * 5)
    r_major, r_minor, _ = cg.coil_radii(sc)
    assert r_major == pytest.approx(5.62, rel=0.01)
    assert r_minor == pytest.approx(1.18, rel=0.05)
    assert cg.total_length(sc) == pytest.approx(420.67, rel=0.03)
    assert cg.coil_height(sc) == pytest.approx(3.1, rel=0.02)
    assert cg.min_coil_coil_distance(sc, classes) == pytest.approx(0.28, rel=0.03)
    assert cg.min_bend_radius(sc) == pytest.approx(0.39, rel=0.05)
    assert cg.plasma_coil_distance(eq, sc) == pytest.approx(0.39, rel=0.03)
    assert cg.coil_surface_area(sc) == pytest.approx(299.85, rel=0.01)
    slope, _, _ = cg.plasma_coil_distance_derivative(eq, sc)
    assert -1.2 < slope < -0.2
