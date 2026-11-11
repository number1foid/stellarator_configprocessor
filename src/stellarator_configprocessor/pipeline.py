"""Run every calculation and assemble the PROCESS configuration."""

from __future__ import annotations

import time
from dataclasses import dataclass, field, replace

import numpy as np

from . import SCHEMA_VERSION, __version__
from .consistency import Tolerances, axis_field, check_consistency, enforce
from .energy import magnetic_energy
from .forces import coil_forces
from .geometry import coils as cg
from .geometry.plasma import surface_geometry
from .io.coils import CoilSet
from .io.equilibrium import Equilibrium
from .schema import KEYS, validate
from .winding_pack import WindingPack, coil_frames, fit_peak_field, winding_pack_field


class SettingsError(ValueError):
    """Run settings are missing or invalid."""


@dataclass(frozen=True)
class Settings:
    """Run settings. The first four have no default on purpose."""

    #: Winding pack width in the coil's in-plane normal direction [m].
    wp_radial: float
    #: Winding pack width in the toroidal direction [m].
    wp_toroidal: float
    #: Effective helical ripple, passed through to PROCESS (not computed).
    epseff: float
    #: Neutron wall load peaking factor, passed through (not computed).
    neutron_peakfactor: float
    name: str = "stellarator_configprocessor"
    #: Rescale the machine to this plasma major radius [m].
    rmajor_ref: float | None = None
    #: Rescale the coil currents to give this mean field on axis [T].
    bt_ref: float | None = None
    #: Scale the coil currents to the equilibrium's field before checking.
    match_field: bool = False
    #: Continue when coils and equilibrium appear inconsistent.
    allow_inconsistent: bool = False
    tolerances: Tolerances = field(default_factory=Tolerances)
    #: Samples per coil.
    n_points: int = 160
    #: Filaments across the winding pack (radial, toroidal).
    wp_filaments: tuple[int, int] = (4, 4)
    #: Winding pack sizes, relative to the reference, used for the a1, a2 fit.
    wp_scan_factors: tuple[float, ...] = (0.6, 0.8, 1.0, 1.25, 1.5)
    #: Angle resolution of the plasma surface for distance calculations.
    surface_resolution: int = 96

    def __post_init__(self):
        for name in ("wp_radial", "wp_toroidal", "neutron_peakfactor"):
            value = getattr(self, name)
            if not (isinstance(value, (int, float)) and np.isfinite(value) and value > 0):
                raise SettingsError(f"{name} must be a positive number, got {value!r}")
        if not (isinstance(self.epseff, (int, float)) and np.isfinite(self.epseff)
                and self.epseff >= 0):
            raise SettingsError(f"epseff must be a non-negative number, got {self.epseff!r}")
        for name in ("rmajor_ref", "bt_ref"):
            value = getattr(self, name)
            if value is not None and not value > 0:
                raise SettingsError(f"{name} must be positive, got {value!r}")
        if self.n_points < 32:
            raise SettingsError("n_points must be at least 32")
        if 1.0 not in self.wp_scan_factors or len(self.wp_scan_factors) < 2:
            raise SettingsError("wp_scan_factors must contain 1.0 and at least one other size")


@dataclass(frozen=True)
class Result:
    #: The values PROCESS reads, keyed as in the file.
    config: dict
    #: Everything else: methods, diagnostics, warnings, settings.
    report: dict


def run(eq: Equilibrium, coilset: CoilSet, settings: Settings, log=None) -> Result:
    """Compute the configuration for one equilibrium and coil set."""
    log = log or (lambda message: None)
    t_start = time.time()
    warnings: list[str] = list(coilset.warnings)

    # ---------------------------------------------------------------- inputs
    n_all = len(coilset)
    coilset = coilset.without_zero_current()
    if len(coilset) < n_all:
        warnings.append(f"{n_all - len(coilset)} coil(s) with zero current were ignored")
    if len(coilset) == 0:
        raise SettingsError("no coil carries current")
    if eq.lasym:
        warnings.append("the equilibrium is not stellarator-symmetric; coils are grouped "
                        "by field-period rotation only")

    log("sampling coils")
    coils = coilset.sampled(settings.n_points)

    # ----------------------------------------------------------- consistency
    log("checking that coils and equilibrium match")
    axis = axis_field(eq, coils)
    check = check_consistency(eq, coils, coilset.periods, settings.tolerances, axis=axis,
                              check_field=not settings.match_field)
    enforce(check, settings.allow_inconsistent)
    if not check.ok:
        warnings.extend("consistency check overridden: " + f for f in check.failures)

    # --------------------------------------------------------------- scaling
    current_factor = 1.0
    if settings.match_field:
        current_factor = check.b_axis_equilibrium / check.b_axis_coils
    b_now = check.b_axis_coils * current_factor
    field_factor = settings.bt_ref / b_now if settings.bt_ref else 1.0
    plasma_input = surface_geometry(eq)
    length_factor = settings.rmajor_ref / plasma_input.major_radius if settings.rmajor_ref else 1.0
    scaling = {
        "length_factor": length_factor,
        "field_factor": field_factor,
        "current_factor_to_match_equilibrium": current_factor,
        # Field ~ current / length, so currents follow both factors.
        "total_current_factor": current_factor * field_factor * length_factor,
    }
    if length_factor != 1.0 or scaling["total_current_factor"] != 1.0:
        eq = eq.scaled(length_factor, field_factor)
        coilset = coilset.scaled(length_factor, scaling["total_current_factor"])
        coils = coilset.sampled(settings.n_points)
        axis = axis_field(eq, coils)

    wp = WindingPack(settings.wp_radial, settings.wp_toroidal)

    # ---------------------------------------------------------------- plasma
    log("plasma geometry")
    plasma = surface_geometry(eq)

    # ------------------------------------------------------------ coil shape
    log("coil geometry")
    classes = cg.symmetry_classes(coils, eq.nfp, stellsym=not eq.lasym)
    reps = classes.representatives
    r_major, r_minor, _ = cg.coil_radii(coils)
    d_min = cg.min_coil_coil_distance(coils, classes)
    res = settings.surface_resolution
    d_pc = cg.plasma_coil_distance(eq, coils, 1.0, res, res)
    slope, fit_radii, fit_dists = cg.plasma_coil_distance_derivative(
        eq, coils, n_theta=res, n_phi=res)
    if wp.toroidal >= d_min:
        warnings.append(
            f"the winding pack toroidal width ({wp.toroidal:.4g} m) is not smaller than the "
            f"minimum coil-coil distance ({d_min:.4g} m): neighbouring winding packs overlap")
    if 0.5 * wp.radial >= d_pc:
        warnings.append("half the winding pack radial width exceeds the plasma-coil distance")
    if len(coils.currents) and np.ptp(np.abs(coils.currents)) > 1e-6 * np.abs(coils.currents).max():
        warnings.append("coil currents differ between coils; PROCESS scales a single "
                        "current, taken here as the mean absolute current")

    # ------------------------------------------------- field in winding pack
    log("field in the winding pack")
    frames = coil_frames(coils)
    pack = winding_pack_field(coils, wp, reps, frames=frames, n_sub=settings.wp_filaments)
    b_max = pack.peak()

    log("forces")
    forces = coil_forces(coils, pack, frames, classes)

    # ------------------------------------------------------- peak-field fit
    log("peak field for several winding pack sizes")
    station_peaks = pack.station_peaks()
    candidates = tuple(
        np.union1d(np.nonzero(p >= 0.8 * b_max)[0], [int(np.argmax(p))]) for p in station_peaks)
    widths, peaks, skipped = [], [], []
    for factor in settings.wp_scan_factors:
        if factor == 1.0:
            widths.append(wp.radial)
            peaks.append(b_max)
            continue
        scaled = wp.scaled(factor)
        if scaled.toroidal >= 0.98 * d_min:
            skipped.append(factor)
            continue
        scan = winding_pack_field(coils, scaled, reps, frames=frames,
                                  n_sub=settings.wp_filaments, stations=candidates)
        widths.append(scaled.radial)
        peaks.append(scan.peak())
    if skipped:
        warnings.append("winding pack sizes skipped in the a1, a2 fit because they would "
                        f"overlap neighbouring coils: factors {skipped}")
    if len(widths) < 2:
        raise SettingsError("fewer than two usable winding pack sizes for the a1, a2 fit; "
                            "add smaller wp_scan_factors")
    i0_ma = float(np.mean(np.abs(coils.currents))) / 1.0e6
    fit = fit_peak_field(widths, peaks, i0_ma, coils.n_coils, r_major, r_minor)

    # ---------------------------------------------------------------- energy
    log("magnetic energy")
    energy = magnetic_energy(coils, wp, classes)

    # --------------------------------------------------------------- results
    if coils.n_coils % eq.nfp != 0:
        raise SettingsError(f"{coils.n_coils} coils cannot be divided into {eq.nfp} field periods")
    mn = 1.0e-6
    config = {
        "name": settings.name,
        "symmetry": int(eq.nfp),
        "coilspermodule": int(coils.n_coils // eq.nfp),
        "rmajor_ref": plasma.major_radius,
        "rminor_ref": plasma.minor_radius,
        "aspect_ref": plasma.aspect_ratio,
        "bt_ref": axis.mean_strength,
        "vol_plasma": plasma.volume,
        "plasma_surface": plasma.area,
        "coil_rmajor": r_major,
        "coil_rminor": r_minor,
        "coillength": cg.total_length(coils),
        "coilsurface": cg.coil_surface_area(coils),
        "maximal_coil_height": cg.coil_height(coils),
        "dmin": d_min,
        "max_portsize_width": cg.port_width(coils),
        "min_bend_radius": cg.min_bend_radius(coils),
        "min_plasma_coil_distance": d_pc,
        "derivative_min_LCFS_coils_dist": slope,
        "WP_area": wp.area,
        "WP_ratio": wp.ratio,
        "i0": i0_ma,
        "WP_bmax": b_max,
        "a1": fit.a1,
        "a2": fit.a2,
        "inductance": energy.effective_inductance,
        "max_force_density": forces.max_force_density * mn,
        "max_force_density_MNm": forces.max_line_force * mn,
        "max_lateral_force_density": forces.max_lateral_force_density * mn,
        "max_radial_force_density": forces.max_radial_force_density * mn,
        "centering_force_max_MN": forces.centering_max * mn,
        "centering_force_min_MN": forces.centering_min * mn,
        "centering_force_avg_MN": forces.centering_avg * mn,
        "epseff": float(settings.epseff),
        "neutron_peakfactor": float(settings.neutron_peakfactor),
    }
    config = {k: (v if isinstance(v, (str, int)) else float(v)) for k, v in config.items()}
    validate(config)

    report = {
        "generator": f"stellarator_configprocessor {__version__}",
        "schema_version": SCHEMA_VERSION,
        "inputs": {"equilibrium": eq.source, "coils": coilset.source},
        "settings": {
            **{k: v for k, v in replace(settings).__dict__.items() if k != "tolerances"},
            "tolerances": dict(settings.tolerances.__dict__),
        },
        "values": [
            {"key": k.name, "value": config[k.name], "unit": k.unit,
             "description": k.description, "method": k.method}
            for k in KEYS
        ],
        "consistency": check.as_dict(),
        "scaling": scaling,
        "diagnostics": {
            "plasma_as_reported_by_equilibrium": dict(eq.reported),
            "field_on_axis_mean_toroidal_component_T": axis.mean_toroidal,
            "magnetic_energy_GJ": energy.energy * 1.0e-9,
            "coil_classes": [
                {"representative_coil": int(i), "count": int(m),
                 "current_MA": float(coils.currents[i]) * 1.0e-6,
                 "length_m": float(coils.lengths[i]),
                 "self_inductance_H": float(energy.self_inductance[c]),
                 "peak_field_T": float(station_peaks[c].max()),
                 "centering_force_MN": float(forces.centering[c]) * mn}
                for c, (i, m) in enumerate(zip(reps, classes.multiplicity))
            ],
            "peak_field_fit": {
                "radial_width_m": [float(x) for x in fit.radial_widths],
                "peak_field_T": [float(x) for x in fit.peak_fields],
                "fitted_T": [float(x) for x in fit.fitted],
                "max_relative_error": fit.max_relative_error,
            },
            "plasma_coil_distance_fit": {
                "minor_radius_m": [float(x) for x in fit_radii],
                "distance_m": [float(x) for x in fit_dists],
            },
            "coil_surface_lofted_between_filaments_m2": cg.lofted_surface_area(coils),
            "mean_coil_spacing_at_coil_rmajor_m": 2.0 * np.pi * r_major / coils.n_coils,
        },
        "warnings": warnings,
        "runtime_s": round(time.time() - t_start, 2),
    }
    return Result(config=config, report=report)
