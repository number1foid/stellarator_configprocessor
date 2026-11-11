"""Argument checks at the entry of the numerical routines."""
import numpy as np
import pytest

from stellarator_configprocessor.energy import self_inductance
from stellarator_configprocessor.field import biot_savart
from stellarator_configprocessor.geometry.coils import (
    coil_distances, point_segment_distance, symmetry_classes)
from stellarator_configprocessor.io.coils import SampledCoils
from stellarator_configprocessor.winding_pack import (
    WindingPack, bar_kernel, fit_peak_field, winding_pack_field)

from helpers import circular_tf_coils


@pytest.fixture(scope="module")
def coils():
    return circular_tf_coils(n_coils=6).sampled(32)


def test_biot_savart_rejects_wrong_shapes():
    seg = np.zeros((4, 3))
    with pytest.raises(ValueError, match="shape"):
        biot_savart(np.zeros((5, 2)), seg, seg + 1.0, 1.0)
    with pytest.raises(ValueError, match="same number"):
        biot_savart(np.zeros((5, 3)), seg, np.ones((3, 3)), 1.0)
    with pytest.raises(ValueError, match="one value per segment"):
        biot_savart(np.zeros((5, 3)), seg, seg + 1.0, np.ones(2))


def test_biot_savart_rejects_non_finite_input():
    seg = np.zeros((1, 3))
    with pytest.raises(ValueError, match="non-finite"):
        biot_savart(np.array([[np.nan, 0.0, 0.0]]), seg, seg + 1.0, 1.0)
    with pytest.raises(ValueError, match="non-finite"):
        biot_savart(np.ones((1, 3)), seg, seg + 1.0, np.inf)


def test_point_segment_distance_needs_segments():
    with pytest.raises(ValueError, match="non-zero number"):
        point_segment_distance(np.zeros((2, 3)), np.zeros((0, 3)), np.zeros((0, 3)))


def test_coil_index_is_checked(coils):
    with pytest.raises(ValueError, match="coil index"):
        coil_distances(coils, coils.n_coils)
    with pytest.raises(ValueError, match="coil indices"):
        winding_pack_field(coils, WindingPack(0.1, 0.1), [coils.n_coils])


def test_symmetry_classes_needs_a_period_count(coils):
    with pytest.raises(ValueError, match="nfp"):
        symmetry_classes(coils, 0)


def test_sampled_coils_shapes_are_checked(coils):
    with pytest.raises(ValueError, match="shape"):
        SampledCoils(points=coils.points[0], tangents=coils.tangents[0],
                     curvature=coils.curvature[0], lengths=coils.lengths,
                     currents=coils.currents)
    with pytest.raises(ValueError, match="one value per coil"):
        SampledCoils(points=coils.points, tangents=coils.tangents,
                     curvature=coils.curvature, lengths=coils.lengths,
                     currents=coils.currents[:-1])


def test_kernels_and_fit_reject_bad_sizes():
    with pytest.raises(ValueError, match="width_u"):
        bar_kernel(0.0, 0.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="resolution"):
        self_inductance(np.eye(3), WindingPack(0.1, 0.1), resolution=0.0)
    with pytest.raises(ValueError, match="major radius"):
        fit_peak_field([0.1, 0.2], [5.0, 4.0], 1.0, 10, r_major=1.0, r_minor=2.0)
