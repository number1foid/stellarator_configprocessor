"""Cross-check against SIMSOPT, an independent implementation.

SIMSOPT is used only here. The tests are skipped when it is not installed.
"""

import numpy as np
import pytest

from helpers import W7X_COILS, W7X_WOUT
from stellarator_configprocessor.field import biot_savart
from stellarator_configprocessor.io.coils import read_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout

simsopt = pytest.importorskip("simsopt")
pytestmark = pytest.mark.w7x


@pytest.fixture(scope="module")
def simsopt_w7x():
    """SIMSOPT's W7-X coils: base curves, base currents and all coils."""
    from simsopt import configs
    from simsopt.field import BiotSavart, coils_via_symmetries

    if hasattr(configs, "get_data"):
        base_curves, base_currents, *_ = configs.get_data("w7x")
    else:
        base_curves, base_currents, _ = configs.get_w7x_data()
    # The first five curves are the non-planar coils; the planar coils carry
    # no current in the standard configuration.
    coils = coils_via_symmetries(base_curves[:5], base_currents[:5], 5, True)
    return base_curves[:5], BiotSavart(coils)


def test_field_on_axis_matches_simsopt(simsopt_w7x):
    _, bs = simsopt_w7x
    eq = read_wout(W7X_WOUT)
    points = np.ascontiguousarray(eq.axis_xyz(np.linspace(0, 2 * np.pi / 5, 25)))
    bs.set_points(points)
    reference = np.asarray(bs.B())
    sc = read_makegrid(W7X_COILS).sampled(400)
    ours = biot_savart(points, *sc.segments())
    np.testing.assert_allclose(np.linalg.norm(ours, axis=1),
                               np.linalg.norm(reference, axis=1), rtol=1e-3)
    # Direction: allow for an overall sign convention difference.
    cos = np.sum(ours * reference, axis=1) / (
        np.linalg.norm(ours, axis=1) * np.linalg.norm(reference, axis=1))
    assert np.all(np.abs(cos) > 1 - 1e-5)
    assert np.all(np.sign(cos) == np.sign(cos[0]))


def test_coil_lengths_match_simsopt(simsopt_w7x):
    from simsopt.geo import CurveLength

    base_curves, _ = simsopt_w7x
    sc = read_makegrid(W7X_COILS).sampled(400)
    reference = sorted(float(CurveLength(c).J()) for c in base_curves)
    ours = sorted(set(np.round(sc.lengths, 6)))
    assert len(ours) == 5
    np.testing.assert_allclose(ours, reference, rtol=1e-4)
