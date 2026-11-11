"""The ``stella_conf.json`` schema read by PROCESS, and its writer.

See ``docs/schema.md`` for the derivation of every entry from the PROCESS
source. Key spelling follows the example file shipped with PROCESS; PROCESS
itself lower-cases keys on reading.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from . import SCHEMA_VERSION, __version__


class SchemaError(ValueError):
    """A configuration does not satisfy the schema."""


@dataclass(frozen=True)
class Key:
    name: str
    kind: type  # str, int or float
    unit: str
    description: str
    method: str
    positive: bool = False


KEYS: tuple[Key, ...] = (
    Key("name", str, "", "Configuration name", "run setting"),
    Key("symmetry", int, "", "Number of field periods", "from the equilibrium", True),
    Key("coilspermodule", int, "", "Coils per field period",
        "coil count divided by the number of field periods", True),
    Key("rmajor_ref", float, "m", "Reference plasma major radius",
        "volume / (2 pi x mean cross-section) of the last closed flux surface", True),
    Key("rminor_ref", float, "m", "Reference plasma minor radius",
        "radius of a circle with the mean cross-section", True),
    Key("aspect_ref", float, "", "Reference aspect ratio", "rmajor_ref / rminor_ref", True),
    Key("bt_ref", float, "T", "Reference field on axis",
        "arc-length mean of |B| of the coils along the magnetic axis", True),
    Key("vol_plasma", float, "m^3", "Plasma volume",
        "integral over the last closed flux surface (Lion 2021, eq. 3)", True),
    Key("plasma_surface", float, "m^2", "Plasma surface area",
        "integral over the last closed flux surface (Lion 2021, eq. 4)", True),
    Key("coil_rmajor", float, "m", "Mean coil major radius",
        "arc-length-weighted mean cylindrical radius of the filaments", True),
    Key("coil_rminor", float, "m", "Mean coil minor radius",
        "arc-length-weighted mean distance from the circle R = coil_rmajor", True),
    Key("coillength", float, "m", "Total length of all coils",
        "sum of filament lengths", True),
    Key("coilsurface", float, "m^2", "Area of the toroidal surface carrying the coils",
        "mean surface of revolution of the filaments: 2 pi x coil_rmajor x mean coil length",
        True),
    Key("maximal_coil_height", float, "m", "Full height of the coil set",
        "max(z) - min(z) of the filaments", True),
    Key("dmin", float, "m", "Minimum distance between two coil filaments",
        "shortest distance between any two different filaments", True),
    Key("max_portsize_width", float, "m", "Toroidal port width",
        "largest distance between the centroids of neighbouring coils", True),
    Key("min_bend_radius", float, "m", "Minimum bending radius",
        "1 / max curvature of a periodic spline through each filament (Lion 2021, eq. 47)", True),
    Key("min_plasma_coil_distance", float, "m", "Minimum plasma-coil distance",
        "shortest distance between the last closed flux surface and any filament", True),
    Key("derivative_min_LCFS_coils_dist", float, "",
        "Derivative of the plasma-coil distance with respect to the minor radius",
        "linear fit over flux surfaces near the boundary (Lion 2021, eq. 58)"),
    Key("WP_area", float, "m^2", "Reference winding pack area",
        "radial x toroidal width (run settings)", True),
    Key("WP_ratio", float, "", "Winding pack radial / toroidal width", "run settings", True),
    Key("i0", float, "MA", "Current per coil for bt_ref at rmajor_ref",
        "mean absolute coil current", True),
    Key("WP_bmax", float, "T", "Peak field in the reference winding pack",
        "maximum of |B| over the winding pack of all coils (Lion 2021, eq. 38)", True),
    Key("a1", float, "", "Peak-field fit parameter",
        "least-squares fit of the PROCESS formula over several winding pack sizes"),
    Key("a2", float, "", "Peak-field fit parameter (winding pack size term)",
        "least-squares fit of the PROCESS formula over several winding pack sizes"),
    Key("inductance", float, "H", "Effective inductance of the coil set",
        "2 x magnetic energy / i0^2", True),
    Key("max_force_density", float, "MN/m^3", "Maximum mean force density in the winding pack",
        "max_force_density_MNm / WP_area", True),
    Key("max_force_density_MNm", float, "MN/m", "Maximum force per unit length",
        "maximum of |F'| over all coils (Lion 2021, eq. 40)", True),
    Key("max_lateral_force_density", float, "MN/m^3", "Maximum lateral force density",
        "maximum toroidal component of F' / WP_area"),
    Key("max_radial_force_density", float, "MN/m^3", "Maximum radial force density",
        "maximum in-plane normal component of F' / WP_area", True),
    Key("centering_force_max_MN", float, "MN", "Largest centering force on a coil",
        "net coil force towards the torus axis"),
    Key("centering_force_min_MN", float, "MN", "Smallest centering force on a coil",
        "net coil force towards the torus axis; negative is outwards"),
    Key("centering_force_avg_MN", float, "MN", "Mean centering force per coil",
        "average over all coils"),
    Key("epseff", float, "", "Effective helical ripple", "run setting, not computed"),
    Key("neutron_peakfactor", float, "", "Neutron wall load peaking factor",
        "run setting, not computed", True),
)

KEY_NAMES = tuple(k.name for k in KEYS)

#: Extra entries written for traceability. PROCESS ignores unknown keys.
META_KEYS = ("_schema_version", "_generator")


def validate(config: dict) -> None:
    """Check that ``config`` has exactly the schema keys with valid values."""
    problems = []
    missing = [k for k in KEY_NAMES if k not in config]
    extra = [k for k in config if k not in KEY_NAMES and k not in META_KEYS]
    if missing:
        problems.append("missing keys: " + ", ".join(missing))
    if extra:
        problems.append("unknown keys: " + ", ".join(extra))
    for key in KEYS:
        if key.name not in config:
            continue
        value = config[key.name]
        if key.kind is str:
            if not isinstance(value, str):
                problems.append(f"{key.name} must be a string")
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            problems.append(f"{key.name} must be a number, got {type(value).__name__}")
            continue
        if key.kind is int and (not isinstance(value, int)):
            problems.append(f"{key.name} must be an integer")
        if not math.isfinite(value):
            problems.append(f"{key.name} is not finite")
        elif key.positive and not value > 0:
            problems.append(f"{key.name} must be positive, got {value}")
    if problems:
        raise SchemaError("invalid configuration:\n  - " + "\n  - ".join(problems))


def with_metadata(config: dict) -> dict:
    out = {k: config[k] for k in KEY_NAMES}
    out["_schema_version"] = SCHEMA_VERSION
    out["_generator"] = f"stellarator_configprocessor {__version__}"
    return out


def write_config(path: str | Path, config: dict) -> Path:
    """Validate and write the configuration file PROCESS reads with ``istell = 6``."""
    validate(config)
    path = Path(path)
    path.write_text(json.dumps(with_metadata(config), indent=2) + "\n")
    return path


def process_file_name(input_file: str | Path) -> Path:
    """Name PROCESS expects for the configuration of a given input file.

    PROCESS opens ``<prefix>stella_conf.json`` where the prefix is the input
    file name without ``IN.DAT``: ``run1.IN.DAT`` -> ``run1.stella_conf.json``.
    """
    input_file = Path(input_file)
    name = input_file.name
    if not name.endswith("IN.DAT"):
        raise SchemaError(f"a PROCESS input file name ends in 'IN.DAT', got '{name}'")
    return input_file.with_name(name[: -len("IN.DAT")] + "stella_conf.json")
