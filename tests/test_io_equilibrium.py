import numpy as np
import pytest

from helpers import write_torus_wout
from stellarator_configprocessor.io.equilibrium import EquilibriumError, read_wout


def test_reads_torus(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "wout_t.nc", r0=10.0, a=2.0, nfp=5, b0=5.0))
    assert eq.nfp == 5 and not eq.lasym and eq.ns == 11
    theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    phi = np.array([0.0, 0.7])
    surf = eq.surface(1.0, theta, phi)
    np.testing.assert_allclose(surf["R"][:, 0], 10 + 2 * np.cos(theta), atol=1e-12)
    np.testing.assert_allclose(surf["Z"][:, 1], 2 * np.sin(theta), atol=1e-12)
    np.testing.assert_allclose(surf["R_t"][:, 0], -2 * np.sin(theta), atol=1e-12)
    np.testing.assert_allclose(surf["Z_t"][:, 0], 2 * np.cos(theta), atol=1e-12)
    np.testing.assert_allclose(surf["R_p"], 0, atol=1e-12)


def test_surface_at_grid_point_is_exact(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", a=2.0, ns=11))
    surf = eq.surface(0.5, np.array([0.0]), np.array([0.0]))
    np.testing.assert_allclose(surf["R"][0, 0], 10 + 2 * np.sqrt(0.5))


def test_axis_and_axis_field(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=7.0, b0=3.0))
    phi = np.array([0.0, np.pi / 2])
    np.testing.assert_allclose(eq.axis_xyz(phi), [[7, 0, 0], [0, 7, 0]], atol=1e-12)
    np.testing.assert_allclose(eq.axis_field_strength(phi), 3.0)


def test_scaled(tmp_path):
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=10.0, a=2.0, b0=5.0))
    big = eq.scaled(length_factor=2.0, field_factor=0.5)
    np.testing.assert_allclose(big.axis_xyz(np.array([0.0])), [[20, 0, 0]], atol=1e-12)
    np.testing.assert_allclose(big.axis_field_strength(np.array([0.0])), 2.5)
    assert big.reported["volume_p"] == pytest.approx(8 * eq.reported["volume_p"])
    assert big.reported["b0"] == pytest.approx(2.5)


def test_missing_file_and_bad_file(tmp_path):
    with pytest.raises(EquilibriumError, match="not found"):
        read_wout(tmp_path / "nope.nc")
    bad = tmp_path / "bad.nc"
    bad.write_text("not netcdf")
    with pytest.raises(EquilibriumError, match="cannot open"):
        read_wout(bad)


def test_missing_variable_is_named(tmp_path):
    import netCDF4
    path = tmp_path / "empty.nc"
    with netCDF4.Dataset(str(path), "w") as ds:
        ds.createVariable("nfp", "i4")[...] = 5
    with pytest.raises(EquilibriumError, match="'xm'"):
        read_wout(path)
