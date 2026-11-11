"""Second real machine: HSX (4 field periods, 48 coils, 1 T)."""

import pytest

from helpers import HSX_COILS, HSX_WOUT
from stellarator_configprocessor.geometry import coils as cg
from stellarator_configprocessor.io.coils import read_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout
from stellarator_configprocessor.pipeline import Settings, run
from stellarator_configprocessor.winding_pack import WindingPack, pack_clearance

pytestmark = pytest.mark.hsx


@pytest.fixture(scope="module")
def result():
    return run(read_wout(HSX_WOUT), read_makegrid(HSX_COILS),
               Settings(wp_radial=0.06, wp_toroidal=0.05, epseff=0.005, neutron_peakfactor=1.6,
                        name="HSX", wp_filaments=(3, 3)))


@pytest.mark.slow
def test_hsx_is_consistent_and_symmetric(result):
    check = result.report["consistency"]
    assert check["passed"]
    assert check["nfp_equilibrium"] == 4 and check["n_coils"] == 48
    assert check["field_ratio_coils_over_equilibrium"] == pytest.approx(1.0, abs=2e-3)
    assert check["axis_misalignment"] < 5e-3
    classes = result.report["diagnostics"]["coil_classes"]
    assert [k["count"] for k in classes] == [8] * 6
    assert result.report["warnings"] == []


@pytest.mark.slow
def test_hsx_values(result):
    c = result.config
    eq = read_wout(HSX_WOUT)
    assert c["symmetry"] == 4 and c["coilspermodule"] == 12
    assert c["rmajor_ref"] == pytest.approx(eq.reported["Rmajor_p"], rel=1e-6)
    assert c["rminor_ref"] == pytest.approx(eq.reported["Aminor_p"], rel=1e-6)
    assert c["vol_plasma"] == pytest.approx(eq.reported["volume_p"], rel=1e-6)
    assert c["bt_ref"] == pytest.approx(1.0034, rel=2e-3)
    assert c["i0"] == pytest.approx(0.15007255, rel=1e-6)
    # The plasma lies inside the coils and the coils do not touch.
    assert 0 < c["min_plasma_coil_distance"] < c["coil_rminor"] - c["rminor_ref"]
    assert c["dmin"] > 0.05
    # Dimensionless results in the range of PROCESS's built-in machines.
    assert 0.5 < c["a1"] < 1.0
    assert 0.02 < c["a2"] < 0.045
    assert 1.5 < c["WP_bmax"] / c["bt_ref"] < 3.6
    assert -1.2 < c["derivative_min_LCFS_coils_dist"] < -0.2
    assert result.report["diagnostics"]["peak_field_fit"]["max_relative_error"] < 0.01
    assert c["max_force_density"] == pytest.approx(c["max_force_density_MNm"] / c["WP_area"])


def test_hsx_trial_pack_collision_is_detected():
    """The case that exposed the overlap defect: packs meet diagonally."""
    sc = read_makegrid(HSX_COILS).sampled(160)
    classes = cg.symmetry_classes(sc, 4)
    reference = WindingPack(0.06, 0.05)
    d_min = cg.min_coil_coil_distance(sc, classes)
    assert pack_clearance(sc, reference, classes.representatives) > 0.01
    big = reference.scaled(1.5)
    # The toroidal width alone is below the filament distance ...
    assert big.toroidal < d_min
    # ... yet the packs do not keep one filament cell of clearance.
    assert pack_clearance(sc, big, classes.representatives) < big.radial / 4
