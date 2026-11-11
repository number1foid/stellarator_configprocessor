import numpy as np
import pytest

from helpers import W7X_COILS, W7X_WOUT, circular_tf_coils, ideal_toroid_field, write_torus_wout
from stellarator_configprocessor.consistency import (
    ConsistencyError, Tolerances, axis_field, check_consistency, enforce)
from stellarator_configprocessor.io.coils import Coil, CoilSet, read_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout

N, R0, CUR = 20, 10.0, 1.0e6
B_IDEAL = ideal_toroid_field(N, CUR, R0)


@pytest.fixture(scope="module")
def coils():
    return circular_tf_coils(n_coils=N, r0=R0, current=CUR).sampled(120)


def torus(tmp_path, **kw):
    kw.setdefault("r0", R0)
    kw.setdefault("b0", B_IDEAL)
    return read_wout(write_torus_wout(tmp_path / "w.nc", **kw))


def test_axis_field_of_toroidal_coils(tmp_path, coils):
    ax = axis_field(torus(tmp_path), coils)
    assert ax.mean_strength == pytest.approx(B_IDEAL, rel=2e-3)
    assert abs(ax.mean_toroidal) == pytest.approx(ax.mean_strength, rel=1e-4)


def test_consistent_pair_passes(tmp_path, coils):
    res = check_consistency(torus(tmp_path), coils, periods_in_file=5)
    assert res.ok, res.failures
    assert res.linkage_ratio == pytest.approx(1.0, abs=1e-3)
    # Discrete coils make the true axis wiggle radially between coils.
    assert res.axis_misalignment < 5e-3
    assert res.field_ratio == pytest.approx(1.0, abs=2e-3)
    enforce(res)


def test_wrong_periods_fail(tmp_path, coils):
    res = check_consistency(torus(tmp_path), coils, periods_in_file=4)
    assert any("field periods differ" in f for f in res.failures)
    res = check_consistency(torus(tmp_path, nfp=3), coils)
    assert any("cannot be divided" in f for f in res.failures)


def test_field_mismatch_fails_and_can_be_skipped(tmp_path, coils):
    eq = torus(tmp_path, b0=1.08 * B_IDEAL)
    res = check_consistency(eq, coils)
    assert len(res.failures) == 1 and "field on axis differs" in res.failures[0]
    assert check_consistency(eq, coils, check_field=False).ok


def test_axis_outside_coils_fails_linkage(tmp_path, coils):
    res = check_consistency(torus(tmp_path, r0=R0 + 6.0), coils)
    assert any("links" in f for f in res.failures)


def test_tilted_coils_fail_alignment(tmp_path):
    # Tilt every coil about its radial axis: the field at the plasma axis
    # gains a vertical component, which a planar toroidal axis cannot follow.
    base = circular_tf_coils(n_coils=N, r0=R0, current=CUR)
    ang = 0.15
    tilted = []
    for c in base.coils:
        phi = np.arctan2(c.points[:, 1].mean(), c.points[:, 0].mean())
        e_r = np.array([np.cos(phi), np.sin(phi), 0.0])
        d = c.points - R0 * e_r
        rot = (d * np.cos(ang) + np.cross(e_r, d) * np.sin(ang)
               + np.outer(d @ e_r, e_r) * (1 - np.cos(ang)))
        tilted.append(Coil(points=R0 * e_r + rot, current=c.current))
    res = check_consistency(torus(tmp_path), CoilSet(coils=tuple(tilted)).sampled(120),
                            check_field=False)
    assert any("not tangent" in f for f in res.failures)


def test_enforce_raises_unless_overridden(tmp_path, coils):
    res = check_consistency(torus(tmp_path, b0=2 * B_IDEAL), coils)
    with pytest.raises(ConsistencyError, match="field on axis differs"):
        enforce(res)
    enforce(res, allow_inconsistent=True)


def test_tolerances_are_adjustable(tmp_path, coils):
    eq = torus(tmp_path, b0=1.08 * B_IDEAL)
    assert check_consistency(eq, coils, tolerances=Tolerances(field=0.1)).ok


@pytest.mark.w7x
def test_w7x_coils_match_w7x_equilibrium():
    eq = read_wout(W7X_WOUT)
    cs = read_makegrid(W7X_COILS)
    res = check_consistency(eq, cs.sampled(160), cs.periods)
    assert res.ok, res.failures
    assert res.n_coils == 50 and res.nfp_equilibrium == 5
    assert res.field_ratio == pytest.approx(1.0, abs=1e-3)
    assert res.axis_misalignment < 2e-3
    assert res.linkage_ratio == pytest.approx(1.0, abs=1e-3)
