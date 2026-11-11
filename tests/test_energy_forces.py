import numpy as np
import pytest

from helpers import circular_tf_coils
from stellarator_configprocessor.constants import MU0
from stellarator_configprocessor.energy import (
    flux_from_others, gmd_squared, magnetic_energy, self_inductance)
from stellarator_configprocessor.forces import coil_forces
from stellarator_configprocessor.geometry import coils as cg
from stellarator_configprocessor.io.coils import Coil, CoilSet
from stellarator_configprocessor.winding_pack import WindingPack, coil_frames, winding_pack_field


# ------------------------------------------------------------------ helpers
def elliptic_ke(k2):
    """Complete elliptic integrals K and E of parameter m = k^2 (AGM method)."""
    k2 = np.asarray(k2, dtype=float)
    a, b = np.ones_like(k2), np.sqrt(1.0 - k2)
    c2_sum = 0.5 * k2
    power = 0.5
    for _ in range(12):
        a, b, c = 0.5 * (a + b), np.sqrt(a * b), 0.5 * (a - b)
        power *= 2.0
        c2_sum = c2_sum + power * c * c
    big_k = np.pi / (2.0 * a)
    return big_k, big_k * (1.0 - c2_sum)


def loop_field(r, z, radius, current):
    """Exact (B_R, B_z) of a circular filament in the plane z = 0."""
    alpha2 = radius**2 + r**2 + z**2 - 2 * radius * r
    beta2 = radius**2 + r**2 + z**2 + 2 * radius * r
    big_k, big_e = elliptic_ke(1.0 - alpha2 / beta2)
    c = MU0 * current / np.pi
    beta = np.sqrt(beta2)
    b_z = c / (2 * alpha2 * beta) * ((radius**2 - r**2 - z**2) * big_e + alpha2 * big_k)
    with np.errstate(divide="ignore", invalid="ignore"):
        b_r = np.where(r > 0, c * z / (2 * alpha2 * beta * r)
                       * ((radius**2 + r**2 + z**2) * big_e - alpha2 * big_k), 0.0)
    return b_r, b_z


def ring(radius, n=256, centre_x=0.0, plane="xy"):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if plane == "xy":
        return np.column_stack([radius * np.cos(t), radius * np.sin(t), np.zeros(n)])
    return np.column_stack([centre_x + radius * np.cos(t), np.zeros(n), radius * np.sin(t)])


def test_elliptic_integrals():
    big_k, big_e = elliptic_ke(np.array([0.0, 0.5]))
    np.testing.assert_allclose(big_k, [np.pi / 2, 1.8540746773013719], rtol=1e-12)
    np.testing.assert_allclose(big_e, [np.pi / 2, 1.3506438810476755], rtol=1e-12)


# ---------------------------------------------------------------------- GMD
def test_gmd_of_square_is_classical_value():
    assert np.sqrt(gmd_squared(WindingPack(1.0, 1.0))) == pytest.approx(0.44705, abs=1e-5)


def test_gmd_matches_sampling_for_a_rectangle():
    rng = np.random.default_rng(0)
    a, b = 0.6, 0.2
    p = rng.uniform(0, 1, (2_000_000, 4)) * [a, b, a, b]
    mean_log = np.mean(np.log(np.hypot(p[:, 0] - p[:, 2], p[:, 1] - p[:, 3])))
    assert np.sqrt(gmd_squared(WindingPack(a, b))) == pytest.approx(np.exp(mean_log), rel=2e-3)
    assert gmd_squared(WindingPack(a, b)) == pytest.approx(gmd_squared(WindingPack(b, a)))


# --------------------------------------------------------------- inductance
@pytest.mark.parametrize("radius, width", [(1.0, 0.02), (3.0, 0.1)])
def test_self_inductance_of_circular_loop(radius, width):
    wp = WindingPack(width, width)
    analytic = MU0 * radius * (np.log(8 * radius / np.sqrt(gmd_squared(wp))) - 2.0)
    assert self_inductance(ring(radius), wp) == pytest.approx(analytic, rel=2e-3)


def test_mutual_inductance_of_coaxial_loops():
    r1, r2, d = 1.0, 0.7, 0.5
    a = Coil(points=ring(r1), current=1.0)
    b = Coil(points=ring(r2) + [0, 0, d], current=2.5)
    sc = CoilSet(coils=(a, b)).sampled(400)
    k2 = 4 * r1 * r2 / ((r1 + r2) ** 2 + d**2)
    big_k, big_e = elliptic_ke(k2)
    k = np.sqrt(k2)
    mutual = MU0 * np.sqrt(r1 * r2) * ((2 / k - k) * big_k - 2 / k * big_e)
    assert flux_from_others(sc, 0) == pytest.approx(mutual * 2.5, rel=1e-4)
    assert flux_from_others(sc, 1) == pytest.approx(mutual * 1.0, rel=1e-4)


@pytest.mark.slow
def test_energy_from_inductance_matches_volume_integral_of_field():
    """Independent check: W = 1/2 L I^2 against the integral of B^2 / (2 mu0)."""
    radius, width, current = 1.0, 0.2, 1.0e5
    wp = WindingPack(width, width)
    energy_inductance = 0.5 * self_inductance(ring(radius), wp) * current**2

    # Uniform current density in the square cross-section as n x n exact loops.
    n = 10
    h = width / n
    offsets = width * (-0.5 + (np.arange(n) + 0.5) / n)
    # Grid on the cell corners of the sub-loop lattice, so no grid point sits
    # on a sub-loop. Upper half space only; the z = 0 row gets half weight.
    box = 6.0
    r_grid = h * np.arange(1, int(round(box / h)) + 1)
    z_grid = h * np.arange(0, int(round(box / h)) + 1)
    rr, zz = np.meshgrid(r_grid, z_grid, indexing="ij")
    b_r = np.zeros_like(rr)
    b_z = np.zeros_like(rr)
    for dr in offsets:
        for dz in offsets:
            lr, lz = loop_field(rr, zz - dz, radius + dr, current / n**2)
            b_r += lr
            b_z += lz
    weight = np.ones_like(z_grid)
    weight[0] = 0.5
    density = (b_r**2 + b_z**2) / (2 * MU0) * 2 * np.pi * rr
    energy_field = 2.0 * np.sum(density * weight[None, :]) * h * h
    # Dipole field energy outside a sphere of the box size.
    moment = current * np.pi * radius**2
    energy_field += MU0 * moment**2 / (12 * np.pi * box**3)
    assert energy_inductance == pytest.approx(energy_field, rel=0.02)


def test_energy_of_coil_set_uses_symmetry_consistently():
    cs = circular_tf_coils(n_coils=12, r0=10.0, a_coil=3.0, current=2.0e6)
    sc = cs.sampled(96)
    wp = WindingPack(0.5, 0.4)
    full = magnetic_energy(sc, wp, cg.symmetry_classes(sc, nfp=1, stellsym=False))
    reduced = magnetic_energy(sc, wp, cg.symmetry_classes(sc, nfp=3, stellsym=True))
    assert len(full.self_inductance) == 12 and len(reduced.self_inductance) == 2
    assert reduced.energy == pytest.approx(full.energy, rel=1e-9)
    assert reduced.effective_inductance == pytest.approx(2 * full.energy / 2.0e6**2)


def test_energy_approaches_ideal_toroid():
    n, r0, ac, cur = 60, 10.0, 3.0, 1.0e6
    sc = circular_tf_coils(n_coils=n, r0=r0, a_coil=ac, current=cur).sampled(96)
    res = magnetic_energy(sc, WindingPack(0.5, 0.4), cg.symmetry_classes(sc, nfp=n // 2))
    ideal = MU0 * n**2 * (r0 - np.sqrt(r0**2 - ac**2))  # circular-section toroid
    assert res.effective_inductance == pytest.approx(ideal, rel=0.05)


def test_energy_scales_with_current_squared_and_length():
    cs = circular_tf_coils(n_coils=8, r0=10.0, a_coil=3.0, current=1.0e6)
    wp = WindingPack(0.5, 0.4)

    def result(coilset, pack):
        sc = coilset.sampled(64)
        return magnetic_energy(sc, pack, cg.symmetry_classes(sc, nfp=4))

    base = result(cs, wp)
    assert result(cs.scaled(current_factor=3.0), wp).energy == pytest.approx(9 * base.energy)
    bigger = result(cs.scaled(length_factor=2.0), wp.scaled(2.0))
    assert bigger.effective_inductance == pytest.approx(2 * base.effective_inductance, rel=1e-3)


# ------------------------------------------------------------------- forces
def test_hoop_force_on_circular_loop():
    radius, width, current = 2.0, 0.1, 1.0e6
    sc = CoilSet(coils=(Coil(points=ring(radius, 400, centre_x=10.0, plane="xz"),
                             current=current),)).sampled(400)
    wp = WindingPack(width, width)
    classes = cg.symmetry_classes(sc, nfp=1, stellsym=False)
    frames = coil_frames(sc)
    pack = winding_pack_field(sc, wp, classes.representatives, frames=frames, n_sub=(4, 4))
    res = coil_forces(sc, pack, frames, classes)
    # dW/dR at constant current: F' = mu0 I^2 / (4 pi R) (ln(8R/GMD) - 1), outwards.
    analytic = MU0 * current**2 / (4 * np.pi * radius) * (
        np.log(8 * radius / np.sqrt(gmd_squared(wp))) - 1.0)
    np.testing.assert_allclose(res.radial[0], analytic, rtol=0.02)
    np.testing.assert_allclose(res.lateral[0], 0.0, atol=1e-6 * analytic)
    assert np.linalg.norm(res.net_force[0]) < 1e-6 * analytic * 2 * np.pi * radius
    assert res.max_line_force == pytest.approx(analytic, rel=0.02)
    assert res.max_force_density == pytest.approx(analytic / wp.area, rel=0.02)


def test_centering_force_equals_energy_derivative():
    """Virtual work: N * F_centering = -dW/dR0 at constant current.

    The force comes from the Biot-Savart field in the winding pack, the
    energy from Neumann's inductance formula: two independent code paths.
    """
    n, r0, ac, cur = 24, 10.0, 3.0, 1.0e6
    wp = WindingPack(0.4, 0.3)

    def setup(radius):
        sc = circular_tf_coils(n_coils=n, r0=radius, a_coil=ac, current=cur).sampled(128)
        return sc, cg.symmetry_classes(sc, nfp=n // 2)

    sc, classes = setup(r0)
    frames = coil_frames(sc)
    pack = winding_pack_field(sc, wp, classes.representatives, frames=frames, n_sub=(3, 3))
    res = coil_forces(sc, pack, frames, classes)

    dr = 0.01
    w_plus = magnetic_energy(setup(r0 + dr)[0], wp, classes).energy
    w_minus = magnetic_energy(setup(r0 - dr)[0], wp, classes).energy
    # At constant current the outward force on all coils is +dW/dR0.
    centering_from_energy = -(w_plus - w_minus) / (2 * dr) / n
    assert res.centering_avg == pytest.approx(centering_from_energy, rel=0.02)
    assert res.centering_max == pytest.approx(res.centering_min, rel=1e-6)
    # Below the ideal-toroid value, which assumes a continuous current sheet.
    ideal = 0.5 * MU0 * n * cur**2 * (r0 / np.sqrt(r0**2 - ac**2) - 1.0)
    assert 0.5 * ideal < res.centering_avg < ideal
    # Forces lie in the coil plane for an axisymmetric arrangement.
    assert res.max_lateral_force_density < 1e-3 * res.max_radial_force_density


def test_forces_require_full_field():
    sc = circular_tf_coils(n_coils=4, r0=10.0, a_coil=3.0).sampled(32)
    classes = cg.symmetry_classes(sc, nfp=4)
    frames = coil_frames(sc)
    wp = WindingPack(0.4, 0.3)
    partial = winding_pack_field(sc, wp, classes.representatives, frames=frames, n_sub=(2, 2),
                                 stations=(np.arange(4),))
    with pytest.raises(ValueError, match="every sample"):
        coil_forces(sc, partial, frames, classes)
    wrong = winding_pack_field(sc, wp, [1], frames=frames, n_sub=(2, 2))
    with pytest.raises(ValueError, match="representatives"):
        coil_forces(sc, wrong, frames, classes)
