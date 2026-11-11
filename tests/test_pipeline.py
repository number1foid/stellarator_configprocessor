import json

import numpy as np
import pytest

from helpers import W7X_COILS, W7X_WOUT, circular_tf_coils, ideal_toroid_field, write_torus_wout
from stellarator_configprocessor import cli
from stellarator_configprocessor.consistency import ConsistencyError
from stellarator_configprocessor.constants import MU0
from stellarator_configprocessor.io.coils import Coil, CoilSet, read_makegrid, write_makegrid
from stellarator_configprocessor.io.equilibrium import read_wout
from stellarator_configprocessor.pipeline import Settings, SettingsError, run
from stellarator_configprocessor.schema import KEY_NAMES, validate

N, R0, A_PLASMA, A_COIL, CUR = 20, 10.0, 2.0, 4.0, 1.0e6
B_IDEAL = ideal_toroid_field(N, CUR, R0)
FAST = dict(n_points=64, wp_filaments=(2, 2), surface_resolution=48,
            wp_scan_factors=(0.7, 1.0, 1.3))


def settings(**kw):
    base = dict(wp_radial=0.6, wp_toroidal=0.5, epseff=0.01, neutron_peakfactor=1.6, **FAST)
    base.update(kw)
    return Settings(**base)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("pipe")
    wout = write_torus_wout(tmp / "wout_torus.nc", r0=R0, a=A_PLASMA, nfp=5, b0=B_IDEAL, ns=21)
    coils = tmp / "coils.torus"
    write_makegrid(coils, circular_tf_coils(N, R0, A_COIL, CUR), periods=5)
    return wout, coils


@pytest.fixture(scope="module")
def base(inputs):
    wout, coils = inputs
    return run(read_wout(wout), read_makegrid(coils), settings())


def test_config_is_complete_and_valid(base):
    assert list(base.config) == list(KEY_NAMES)
    validate(base.config)
    json.dumps(base.report)  # the report must be serialisable


def test_torus_values_match_analytic(base):
    c = base.config
    assert c["symmetry"] == 5 and c["coilspermodule"] == 4
    assert c["rmajor_ref"] == pytest.approx(R0, rel=1e-9)
    assert c["rminor_ref"] == pytest.approx(A_PLASMA, rel=1e-9)
    assert c["aspect_ref"] == pytest.approx(R0 / A_PLASMA, rel=1e-9)
    assert c["vol_plasma"] == pytest.approx(2 * np.pi**2 * R0 * A_PLASMA**2, rel=1e-9)
    assert c["plasma_surface"] == pytest.approx(4 * np.pi**2 * R0 * A_PLASMA, rel=1e-9)
    assert c["bt_ref"] == pytest.approx(B_IDEAL, rel=2e-3)
    assert c["coil_rmajor"] == pytest.approx(R0, rel=1e-3)
    assert c["coil_rminor"] == pytest.approx(A_COIL, rel=1e-3)
    assert c["coillength"] == pytest.approx(N * 2 * np.pi * A_COIL, rel=1e-3)
    assert c["coilsurface"] == pytest.approx(4 * np.pi**2 * R0 * A_COIL, rel=1e-3)
    assert c["maximal_coil_height"] == pytest.approx(2 * A_COIL, rel=2e-3)
    assert c["dmin"] == pytest.approx(2 * (R0 - A_COIL) * np.sin(np.pi / N), rel=5e-3)
    assert c["max_portsize_width"] == pytest.approx(2 * R0 * np.sin(np.pi / N), rel=1e-3)
    assert c["min_bend_radius"] == pytest.approx(A_COIL, rel=0.01)
    assert c["min_plasma_coil_distance"] == pytest.approx(A_COIL - A_PLASMA, rel=5e-3)
    assert c["derivative_min_LCFS_coils_dist"] == pytest.approx(-1.0, abs=0.01)
    assert c["WP_area"] == pytest.approx(0.3) and c["WP_ratio"] == pytest.approx(1.2)
    assert c["i0"] == pytest.approx(1.0)
    assert c["epseff"] == 0.01 and c["neutron_peakfactor"] == 1.6


def test_torus_field_energy_and_forces_are_plausible(base):
    c = base.config
    ideal_inboard = ideal_toroid_field(N, CUR, R0 - A_COIL)
    assert 1.0 < c["WP_bmax"] / ideal_inboard < 1.6
    ideal_inductance = MU0 * N**2 * (R0 - np.sqrt(R0**2 - A_COIL**2))
    # 20 discrete coils store more energy near the conductors than a current sheet.
    assert 1.0 < c["inductance"] / ideal_inductance < 1.25
    assert c["max_force_density"] == pytest.approx(c["max_force_density_MNm"] / c["WP_area"])
    assert c["centering_force_max_MN"] == pytest.approx(c["centering_force_avg_MN"], rel=1e-6)
    assert c["centering_force_min_MN"] > 0
    assert c["max_lateral_force_density"] < 1e-3 * c["max_radial_force_density"]
    # The fit reproduces the reference peak field through PROCESS's formula.
    process_bmax = (0.2 * c["i0"] * N / (c["coil_rmajor"] - c["coil_rminor"])
                    * (c["a1"] + c["a2"] * c["coil_rmajor"] / 0.6))
    assert process_bmax == pytest.approx(c["WP_bmax"], rel=0.02)
    assert base.report["consistency"]["passed"]
    assert base.report["warnings"] == []


def test_rescaling_to_reference_radius_and_field(inputs, base):
    wout, coils = inputs
    res = run(read_wout(wout), read_makegrid(coils),
              settings(rmajor_ref=2 * R0, bt_ref=2 * base.config["bt_ref"],
                       wp_radial=1.2, wp_toroidal=1.0))
    c, b = res.config, base.config
    assert c["rmajor_ref"] == pytest.approx(2 * b["rmajor_ref"])
    assert c["bt_ref"] == pytest.approx(2 * b["bt_ref"], rel=1e-6)
    assert c["vol_plasma"] == pytest.approx(8 * b["vol_plasma"])
    assert c["coillength"] == pytest.approx(2 * b["coillength"], rel=1e-6)
    assert c["dmin"] == pytest.approx(2 * b["dmin"], rel=1e-6)
    # Field ~ I / L: doubling both the size and the field needs four times the current.
    assert c["i0"] == pytest.approx(4 * b["i0"], rel=1e-6)
    # A geometrically similar machine: fit parameters are unchanged,
    # inductance scales with length, peak field with the axis field.
    assert c["a1"] == pytest.approx(b["a1"], rel=1e-3)
    assert c["a2"] == pytest.approx(b["a2"], rel=1e-3)
    assert c["inductance"] == pytest.approx(2 * b["inductance"], rel=1e-3)
    assert c["WP_bmax"] == pytest.approx(2 * b["WP_bmax"], rel=1e-3)
    # Force per length ~ I * B.
    assert c["max_force_density_MNm"] == pytest.approx(8 * b["max_force_density_MNm"], rel=1e-3)
    assert res.report["scaling"]["length_factor"] == pytest.approx(2.0)


def test_inconsistent_field_stops_unless_overridden(tmp_path, inputs):
    _, coils = inputs
    eq = read_wout(write_torus_wout(tmp_path / "w.nc", r0=R0, a=A_PLASMA, b0=1.5 * B_IDEAL))
    with pytest.raises(ConsistencyError, match="field on axis differs"):
        run(eq, read_makegrid(coils), settings())
    forced = run(eq, read_makegrid(coils), settings(allow_inconsistent=True))
    assert any("overridden" in w for w in forced.report["warnings"])
    assert forced.config["bt_ref"] == pytest.approx(B_IDEAL, rel=2e-3)
    matched = run(eq, read_makegrid(coils), settings(match_field=True))
    assert matched.report["consistency"]["passed"]
    assert matched.config["bt_ref"] == pytest.approx(1.5 * B_IDEAL, rel=1e-6)
    assert matched.config["i0"] == pytest.approx(1.5, rel=2e-3)


def test_zero_current_coils_are_ignored_with_warning(inputs, base):
    wout, coils = inputs
    cs = read_makegrid(coils)
    t = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    extra = Coil(points=np.column_stack([R0 + 6 * np.cos(t), np.zeros(32), 6 * np.sin(t)]),
                 current=0.0)
    res = run(read_wout(wout), CoilSet(coils=cs.coils + (extra,), periods=5), settings())
    assert any("zero current" in w for w in res.report["warnings"])
    assert res.config["coillength"] == pytest.approx(base.config["coillength"])


def test_overlapping_winding_pack_is_reported(inputs):
    wout, coils = inputs
    res = run(read_wout(wout), read_makegrid(coils),
              settings(wp_radial=2.4, wp_toroidal=2.0, wp_scan_factors=(0.5, 0.7, 1.0)))
    assert any("overlap" in w for w in res.report["warnings"])


@pytest.mark.parametrize("kw, message", [
    (dict(wp_radial=0.0), "wp_radial"),
    (dict(wp_toroidal=-1.0), "wp_toroidal"),
    (dict(epseff=-0.1), "epseff"),
    (dict(neutron_peakfactor=0.0), "neutron_peakfactor"),
    (dict(bt_ref=-3.0), "bt_ref"),
    (dict(n_points=8), "n_points"),
    (dict(wp_scan_factors=(0.8, 1.2)), "wp_scan_factors"),
])
def test_invalid_settings(kw, message):
    with pytest.raises(SettingsError, match=message):
        settings(**kw)


def test_required_settings_have_no_default():
    with pytest.raises(TypeError):
        Settings(wp_radial=0.6, wp_toroidal=0.5, epseff=0.01)  # no neutron_peakfactor


# ------------------------------------------------------------------------ CLI
CLI_FAST = ["--points", "64", "--wp-filaments", "2", "2"]
REQUIRED = ["--wp-radial", "0.6", "--wp-toroidal", "0.5", "--epseff", "0.01",
            "--neutron-peakfactor", "1.6"]


def test_cli_writes_config_and_reports(tmp_path, inputs, capsys):
    wout, coils = inputs
    out = tmp_path / "stella_conf.json"
    code = cli.main([str(wout), str(coils), *REQUIRED, *CLI_FAST, "-o", str(out),
                     "--report", str(tmp_path / "report.json"), "--name", "torus"])
    assert code == 0
    data = json.loads(out.read_text())
    assert data["name"] == "torus" and data["symmetry"] == 5
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["consistency"]["passed"] and len(report["values"]) == len(KEY_NAMES)
    text = capsys.readouterr().out
    assert "Values written to the configuration file" in text and "WP_bmax" in text


def test_cli_text_report_and_process_naming(tmp_path, inputs):
    wout, coils = inputs
    in_dat = tmp_path / "myrun.IN.DAT"
    code = cli.main([str(wout), str(coils), *REQUIRED, *CLI_FAST, "-q",
                     "--process-input", str(in_dat), "--report", str(tmp_path / "report.txt")])
    assert code == 0
    assert (tmp_path / "myrun.stella_conf.json").is_file()
    assert "Consistency of coils and equilibrium" in (tmp_path / "report.txt").read_text()


def test_cli_exit_codes(tmp_path, inputs, capsys):
    wout, coils = inputs
    assert cli.main([str(tmp_path / "missing.nc"), str(coils), *REQUIRED, "-q"]) == 2
    assert "not found" in capsys.readouterr().err
    high_field = write_torus_wout(tmp_path / "w.nc", r0=R0, a=A_PLASMA, b0=2 * B_IDEAL)
    out = tmp_path / "x.json"
    assert cli.main([str(high_field), str(coils), *REQUIRED, *CLI_FAST, "-q", "-o", str(out)]) == 3
    assert "--allow-inconsistent" in capsys.readouterr().err
    assert not out.exists()
    assert cli.main([str(high_field), str(coils), *REQUIRED, *CLI_FAST, "-q", "-o", str(out),
                     "--match-field"]) == 0
    assert out.is_file()


def test_cli_requires_conductor_dimensions(inputs, capsys):
    wout, coils = inputs
    with pytest.raises(SystemExit) as exc:
        cli.main([str(wout), str(coils), "--epseff", "0.01", "--neutron-peakfactor", "1.6"])
    assert exc.value.code == 2
    assert "--wp-radial" in capsys.readouterr().err


def test_cli_group_current(tmp_path, inputs):
    wout, coils = inputs
    out = tmp_path / "g.json"
    code = cli.main([str(wout), str(coils), *REQUIRED, *CLI_FAST, "-q", "-o", str(out),
                     "--group-current", "1=2.0", "--match-field"])
    assert code == 0
    # Doubling the file currents and matching the field leaves i0 unchanged.
    assert json.loads(out.read_text())["i0"] == pytest.approx(1.0, rel=2e-3)


# ------------------------------------------------------------------------ W7-X
@pytest.mark.w7x
@pytest.mark.slow
def test_w7x_reproduces_process_builtin_w7x50():
    """As-built W7-X coils at 3 T against PROCESS's built-in W7X-50 values."""
    res = run(read_wout(W7X_WOUT), read_makegrid(W7X_COILS),
              Settings(wp_radial=0.18, wp_toroidal=0.15, epseff=0.0095, neutron_peakfactor=1.6,
                       bt_ref=3.0, name="W7-X", wp_filaments=(3, 3)))
    c = res.config
    assert res.report["consistency"]["passed"]
    assert c["symmetry"] == 5 and c["coilspermodule"] == 10
    builtin = {  # (value, relative tolerance)
        "i0": (1.745, 0.01), "WP_bmax": (6.3, 0.02), "a1": (0.66, 0.05),
        "inductance": (412.4e-6, 0.03), "coil_rmajor": (5.62, 0.01),
        "coil_rminor": (1.18, 0.05), "coillength": (420.67, 0.03),
        "coilsurface": (299.85, 0.01), "maximal_coil_height": (3.1, 0.02),
        "dmin": (0.28, 0.03), "min_bend_radius": (0.39, 0.05),
        "min_plasma_coil_distance": (0.39, 0.03), "centering_force_max_MN": (2.99, 0.05),
        "centering_force_min_MN": (-1.29, 0.05), "centering_force_avg_MN": (1.61, 0.05),
        "max_radial_force_density": (148.0, 0.03), "rmajor_ref": (5.5, 0.01),
    }
    for key, (value, tol) in builtin.items():
        assert c[key] == pytest.approx(value, rel=tol), key
    # Plasma values against VMEC for the same equilibrium.
    eq = read_wout(W7X_WOUT)
    assert c["vol_plasma"] == pytest.approx(eq.reported["volume_p"], rel=1e-6)
    assert c["rminor_ref"] == pytest.approx(eq.reported["Aminor_p"], rel=1e-6)
    assert res.report["diagnostics"]["peak_field_fit"]["max_relative_error"] < 0.01
