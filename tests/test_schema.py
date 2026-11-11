import json

import pytest

from stellarator_configprocessor import SCHEMA_VERSION
from stellarator_configprocessor.schema import (
    KEY_NAMES, KEYS, SchemaError, process_file_name, validate, write_config)

#: Keys of the example file shipped with PROCESS that PROCESS actually reads.
PROCESS_EXAMPLE_KEYS = {
    "name", "min_plasma_coil_distance", "derivative_min_LCFS_coils_dist", "coilspermodule",
    "coil_rmajor", "coil_rminor", "aspect_ref", "bt_ref", "WP_area", "WP_bmax", "i0", "a1",
    "a2", "dmin", "inductance", "coilsurface", "coillength", "max_portsize_width",
    "maximal_coil_height", "WP_ratio", "max_force_density_MNm", "max_force_density",
    "min_bend_radius", "max_lateral_force_density", "max_radial_force_density",
    "centering_force_max_MN", "centering_force_min_MN", "centering_force_avg_MN", "symmetry",
    "rmajor_ref", "rminor_ref", "vol_plasma", "plasma_surface", "epseff", "neutron_peakfactor",
}

#: Field names of PROCESS's StellaratorConfigData (v3.4.2), without the prefix.
PROCESS_DATACLASS_FIELDS = {
    "name", "symmetry", "coilspermodule", "rmajor_ref", "rminor_ref", "coil_rmajor",
    "coil_rminor", "aspect_ref", "bt_ref", "wp_area", "wp_bmax", "i0", "a1", "a2", "dmin",
    "inductance", "coilsurface", "coillength", "max_portsize_width", "maximal_coil_height",
    "min_plasma_coil_distance", "derivative_min_lcfs_coils_dist", "vol_plasma",
    "plasma_surface", "wp_ratio", "max_force_density", "max_force_density_mnm",
    "min_bend_radius", "epseff", "max_lateral_force_density", "max_radial_force_density",
    "centering_force_max_mn", "centering_force_min_mn", "centering_force_avg_mn",
    "neutron_peakfactor",
}


def valid_config():
    config = {}
    for key in KEYS:
        config[key.name] = {str: "test", int: 5, float: 1.5}[key.kind]
    return config


def test_schema_covers_exactly_what_process_reads():
    assert set(KEY_NAMES) == PROCESS_EXAMPLE_KEYS
    assert {k.lower() for k in KEY_NAMES} == PROCESS_DATACLASS_FIELDS
    assert len(KEY_NAMES) == len(set(KEY_NAMES)) == 35


def test_valid_config_passes():
    validate(valid_config())


@pytest.mark.parametrize("change, message", [
    (lambda c: c.pop("i0"), "missing keys: i0"),
    (lambda c: c.update(extra=1.0), "unknown keys: extra"),
    (lambda c: c.update(WP_bmax=float("nan")), "WP_bmax is not finite"),
    (lambda c: c.update(dmin=-0.1), "dmin must be positive"),
    (lambda c: c.update(symmetry=5.0), "symmetry must be an integer"),
    (lambda c: c.update(a1="x"), "a1 must be a number"),
    (lambda c: c.update(name=3), "name must be a string"),
    (lambda c: c.update(a2=True), "a2 must be a number"),
])
def test_invalid_configs_are_rejected(change, message):
    config = valid_config()
    change(config)
    with pytest.raises(SchemaError, match=message):
        validate(config)


def test_signed_values_are_allowed():
    config = valid_config()
    config.update(centering_force_min_MN=-2.0, derivative_min_LCFS_coils_dist=-0.6, a1=-0.1,
                  max_lateral_force_density=0.0)
    validate(config)


def test_write_config(tmp_path):
    path = write_config(tmp_path / "stella_conf.json", valid_config())
    data = json.loads(path.read_text())
    assert list(data)[: len(KEY_NAMES)] == list(KEY_NAMES)
    assert data["_schema_version"] == SCHEMA_VERSION
    assert data["_generator"].startswith("stellarator_configprocessor")
    assert data["symmetry"] == 5 and isinstance(data["symmetry"], int)


def test_write_rejects_invalid(tmp_path):
    config = valid_config()
    del config["bt_ref"]
    with pytest.raises(SchemaError):
        write_config(tmp_path / "x.json", config)
    assert not (tmp_path / "x.json").exists()


def test_process_file_name(tmp_path):
    assert process_file_name(tmp_path / "run1.IN.DAT") == tmp_path / "run1.stella_conf.json"
    assert process_file_name("IN.DAT").name == "stella_conf.json"
    with pytest.raises(SchemaError):
        process_file_name("run1.txt")
