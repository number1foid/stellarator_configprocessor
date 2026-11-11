import numpy as np
import pytest

from stellarator_configprocessor.io.coils import (
    Coil, CoilFileError, CoilSet, expand_symmetry, read_makegrid, write_makegrid)

SQUARE = """! header comment
&MGRID_NLI
   MGRID_EXT = 'test'
/
periods 2
begin filament
mirror NIL
 1.0  0.0 -1.0  2.0D3
 1.0  0.0  1.0  2.0D3
 3.0  0.0  1.0  2.0D3
 3.0  0.0 -1.0  2.0D3
 1.0  0.0 -1.0  0.0 1 Mod
 2.5  0.0  0.4  50.0 2 PF
end
"""


def test_reads_filament_and_circular_coil(tmp_path):
    path = tmp_path / "coils.test"
    path.write_text(SQUARE)
    cs = read_makegrid(path)
    assert cs.periods == 2
    assert len(cs) == 2
    square, ring = cs.coils
    assert square.points.shape == (4, 3)  # closing point removed
    assert square.current == 2000.0 and square.group == 1 and square.name == "Mod"
    assert ring.current == 50.0 and ring.group == 2
    np.testing.assert_allclose(np.hypot(ring.points[:, 0], ring.points[:, 1]), 2.5)
    np.testing.assert_allclose(ring.points[:, 2], 0.4)
    assert cs.warnings == ()


def test_group_currents_scale_each_group(tmp_path):
    path = tmp_path / "coils.test"
    path.write_text(SQUARE)
    cs = read_makegrid(path, group_currents={1: 3.0, 2: 0.5})
    np.testing.assert_allclose(cs.currents, [6000.0, 25.0])


def test_open_coil_is_reported(tmp_path):
    path = tmp_path / "coils.open"
    path.write_text(SQUARE.replace(" 1.0  0.0 -1.0  0.0 1 Mod", " 1.2  0.0 -1.0  0.0 1 Mod"))
    cs = read_makegrid(path)
    assert any("not closed" in w for w in cs.warnings)


@pytest.mark.parametrize("text, message", [
    ("just text\n", "no 'periods' line"),
    ("periods 2\nbegin filament\n1 0 0 5\n2 0 0 5\nend\n", "ends inside a coil"),
    ("periods 2\nbegin filament\n1 0 x 5\nend\n", "cannot parse"),
    ("periods 2\nbegin filament\nend\n", "no coils found"),
])
def test_malformed_files_raise(tmp_path, text, message):
    path = tmp_path / "coils.bad"
    path.write_text(text)
    with pytest.raises(CoilFileError, match=message):
        read_makegrid(path)


def test_missing_file_raises(tmp_path):
    with pytest.raises(CoilFileError, match="not found"):
        read_makegrid(tmp_path / "nope")


def test_write_read_round_trip(tmp_path):
    rng = np.random.default_rng(1)
    t = np.linspace(0, 2 * np.pi, 30, endpoint=False)
    pts = np.column_stack([3 + np.cos(t), 0.1 * rng.standard_normal(30), np.sin(t)])
    cs = CoilSet(coils=(Coil(points=pts, current=-1.5e6, group=3, name="A"),))
    path = tmp_path / "coils.rt"
    write_makegrid(path, cs, periods=4)
    back = read_makegrid(path)
    assert back.periods == 4
    np.testing.assert_allclose(back.coils[0].points, pts, rtol=0, atol=1e-13)
    assert back.coils[0].current == -1.5e6 and back.coils[0].group == 3


def test_zero_current_coil_cannot_be_written(tmp_path):
    cs = CoilSet(coils=(Coil(points=np.eye(3), current=0.0),))
    with pytest.raises(CoilFileError, match="zero current"):
        write_makegrid(tmp_path / "c", cs, periods=1)


def test_expand_symmetry_counts_and_images():
    t = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    phi = 0.2
    r = 5 + np.cos(t)
    pts = np.column_stack([r * np.cos(phi), r * np.sin(phi), np.sin(t) + 0.3])
    full = expand_symmetry(CoilSet(coils=(Coil(points=pts, current=2.0),)), nfp=3)
    assert len(full) == 6 and full.periods == 3
    # Mirror image: (x, -y, -z) with reversed current.
    np.testing.assert_allclose(full.coils[1].points, pts * [1, -1, -1])
    assert full.coils[1].current == -2.0
    # Rotated image keeps R and z.
    np.testing.assert_allclose(np.hypot(*full.coils[2].points[:, :2].T), r)
    np.testing.assert_allclose(full.coils[2].points[:, 2], pts[:, 2])


def test_sampled_coils_are_uniform_in_arc_length():
    t = np.linspace(0, 2 * np.pi, 40, endpoint=False)
    t = t + 0.3 * np.sin(t)  # non-uniform input spacing
    pts = np.column_stack([2 * np.cos(t), 2 * np.sin(t), np.zeros_like(t)])
    sc = CoilSet(coils=(Coil(points=pts, current=1.0),)).sampled(64)
    assert sc.points.shape == (1, 64, 3)
    np.testing.assert_allclose(sc.lengths[0], 4 * np.pi, rtol=1e-4)
    seg = np.linalg.norm(np.roll(sc.points[0], -1, axis=0) - sc.points[0], axis=1)
    assert np.ptp(seg) / seg.mean() < 1e-3
    np.testing.assert_allclose(sc.curvature[0], 0.5, rtol=5e-3)
    np.testing.assert_allclose(np.linalg.norm(sc.tangents[0], axis=1), 1.0)
