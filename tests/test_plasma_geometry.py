import numpy as np
import pytest

from helpers import W7X_WOUT, write_torus_wout
from stellarator_configprocessor.geometry.plasma import surface_geometry
from stellarator_configprocessor.io.equilibrium import read_wout


@pytest.mark.parametrize("r0, a", [(10.0, 2.0), (3.0, 1.2)])
def test_circular_torus_matches_analytic(tmp_path, r0, a):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=r0, a=a))
    g = surface_geometry(eq)
    assert g.volume == pytest.approx(2 * np.pi**2 * r0 * a**2, rel=1e-12)
    assert g.area == pytest.approx(4 * np.pi**2 * r0 * a, rel=1e-12)
    assert g.cross_section == pytest.approx(np.pi * a**2, rel=1e-12)
    assert g.minor_radius == pytest.approx(a, rel=1e-12)
    assert g.major_radius == pytest.approx(r0, rel=1e-12)
    assert g.aspect_ratio == pytest.approx(r0 / a, rel=1e-12)


def test_inner_surface_scales_with_flux(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=10.0, a=2.0))
    g = surface_geometry(eq, s=0.5)
    # a(s) = a * sqrt(s): the volume is linear in s for this torus.
    assert g.volume == pytest.approx(0.5 * 2 * np.pi**2 * 10.0 * 4.0, rel=1e-12)
    assert g.minor_radius == pytest.approx(2.0 * np.sqrt(0.5), rel=1e-12)


def test_scaling_laws(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc"))
    g1, g2 = surface_geometry(eq), surface_geometry(eq.scaled(length_factor=1.7))
    assert g2.volume == pytest.approx(1.7**3 * g1.volume)
    assert g2.area == pytest.approx(1.7**2 * g1.area)
    assert g2.minor_radius == pytest.approx(1.7 * g1.minor_radius)


@pytest.mark.w7x
def test_w7x_matches_values_reported_by_vmec():
    eq = read_wout(W7X_WOUT)
    g = surface_geometry(eq)
    assert g.volume == pytest.approx(eq.reported["volume_p"], rel=1e-3)
    assert g.minor_radius == pytest.approx(eq.reported["Aminor_p"], rel=1e-3)
    assert g.major_radius == pytest.approx(eq.reported["Rmajor_p"], rel=1e-3)
    # Resolution independence.
    fine = surface_geometry(eq, n_theta=256, n_phi=256)
    assert fine.volume == pytest.approx(g.volume, rel=1e-9)
    assert fine.area == pytest.approx(g.area, rel=1e-7)
