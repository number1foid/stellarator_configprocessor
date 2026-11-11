"""Run every calculation and assemble the PROCESS configuration."""

from __future__ import annotations

import time
from dataclasses import dataclass, field, replace

import numpy as np

from . import SCHEMA_VERSION, __version__
from .consistency import (
    AxisField, ConsistencyResult, Tolerances, axis_field, check_consistency, enforce)
from .energy import EnergyResult, magnetic_energy
from .forces import ForceResult, coil_forces
from .geometry import coils as cg
from .geometry.plasma import SurfaceGeometry, surface_geometry
from .io.coils import CoilSet, SampledCoils
from .io.equilibrium import Equilibrium
from .schema import KEYS, validate
from .winding_pack import (
    PackField, PeakFieldFit, WindingPack, coil_frames, fit_peak_field, pack_clearance,
    winding_pack_field)


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
    #: The default only shrinks the pack: larger packs often touch the
    #: neighbouring coils, and PROCESS at reactor scale works with a pack that
    #: is smaller relative to the coil than the reference.
    wp_scan_factors: tuple[float, ...] = (0.5, 0.6, 0.7, 0.85, 1.0)
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


@dataclass(frozen=True)
class _CoilGeometry:
    """Coil shape values shared by several later stages."""

    classes: cg.SymmetryClasses
    r_major: float
    r_minor: float
    d_min: float
    d_pc: float
    slope: float
    fit_radii: np.ndarray
    fit_dists: np.ndarray


@dataclass(frozen=True)
class _Stages:
    """Everything computed by :func:`run` before the output is assembled."""

    eq: Equilibrium
    coilset: CoilSet
    coils: SampledCoils
    axis: AxisField
    check: ConsistencyResult
    scaling: dict
    plasma: SurfaceGeometry
    geometry: _CoilGeometry
    wp: WindingPack
    pack: PackField
    clearance: float
    forces: ForceResult
    fit: PeakFieldFit
    energy: EnergyResult


def _silent(message: str) -> None:
    """Default progress callback: discard the message."""


def _mean_current_ma(coils: SampledCoils) -> float:
    return float(np.mean(np.abs(coils.currents))) / 1.0e6


def _current_carrying(coilset: CoilSet, eq: Equilibrium, warnings: list[str]) -> CoilSet:
    """Drop coils without current and note what the inputs imply."""
    n_all = len(coilset)
    coilset = coilset.without_zero_current()
    if len(coilset) < n_all:
        warnings.append(f"{n_all - len(coilset)} coil(s) with zero current were ignored")
    if len(coilset) == 0:
        raise SettingsError("no coil carries current")
    if eq.lasym:
        warnings.append("the equilibrium is not stellarator-symmetric; coils are grouped "
                        "by field-period rotation only")
    return coilset


def _scaling_factors(eq: Equilibrium, check: ConsistencyResult, settings: Settings) -> dict:
    """Length, field and current factors that bring the machine to the reference."""
    current_factor = 1.0
    if settings.match_field:
        if check.b_axis_coils == 0.0:
            raise SettingsError("the coils give no field on the axis; cannot match the field")
        current_factor = check.b_axis_equilibrium / check.b_axis_coils
    b_now = check.b_axis_coils * current_factor
    if settings.bt_ref and b_now == 0.0:
        raise SettingsError("the coils give no field on the axis; cannot rescale to bt_ref")
    field_factor = settings.bt_ref / b_now if settings.bt_ref else 1.0
    major_radius = surface_geometry(eq).major_radius
    length_factor = settings.rmajor_ref / major_radius if settings.rmajor_ref else 1.0
    return {
        "length_factor": length_factor,
        "field_factor": field_factor,
        "current_factor_to_match_equilibrium": current_factor,
        # Field ~ current / length, so currents follow both factors.
        "total_current_factor": current_factor * field_factor * length_factor,
    }


def _coil_geometry(eq: Equilibrium, coils: SampledCoils, wp: WindingPack, settings: Settings,
                   warnings: list[str]) -> _CoilGeometry:
    classes = cg.symmetry_classes(coils, eq.nfp, stellsym=not eq.lasym)
    r_major, r_minor, _ = cg.coil_radii(coils)
    res = settings.surface_resolution
    d_pc = cg.plasma_coil_distance(eq, coils, 1.0, res, res)
    slope, fit_radii, fit_dists = cg.plasma_coil_distance_derivative(
        eq, coils, n_theta=res, n_phi=res)
    if 0.5 * wp.radial >= d_pc:
        warnings.append("half the winding pack radial width exceeds the plasma-coil distance")
    magnitude = np.abs(coils.currents)
    if np.ptp(magnitude) > 1e-6 * magnitude.max():
        warnings.append("coil currents differ between coils; PROCESS scales a single "
                        "current, taken here as the mean absolute current")
    return _CoilGeometry(classes=classes, r_major=r_major, r_minor=r_minor,
                         d_min=cg.min_coil_coil_distance(coils, classes), d_pc=d_pc,
                         slope=slope, fit_radii=fit_radii, fit_dists=fit_dists)


def _peak_field_fit(coils: SampledCoils, wp: WindingPack, pack: PackField, frames,
                    geometry: _CoilGeometry, settings: Settings,
                    warnings: list[str]) -> PeakFieldFit:
    """Peak field for several winding pack sizes and the fit of a1 and a2."""
    reps = geometry.classes.representatives
    b_max = pack.peak()
    candidates = tuple(
        np.union1d(np.nonzero(p >= 0.8 * b_max)[0], [int(np.argmax(p))])
        for p in pack.station_peaks())
    widths, peaks, skipped = [], [], []
    for factor in settings.wp_scan_factors:
        if factor == 1.0:
            widths.append(wp.radial)
            peaks.append(b_max)
            continue
        scaled = wp.scaled(factor)
        # A trial pack must stay clear of its neighbours by one filament cell,
        # otherwise the neighbour's discrete filaments distort the field.
        cell = max(scaled.radial / settings.wp_filaments[0],
                   scaled.toroidal / settings.wp_filaments[1])
        if pack_clearance(coils, scaled, reps, frames) < cell:
            skipped.append(factor)
            continue
        scan = winding_pack_field(coils, scaled, reps, frames=frames,
                                  n_sub=settings.wp_filaments, stations=candidates)
        widths.append(scaled.radial)
        peaks.append(scan.peak())
    if skipped:
        warnings.append("winding pack sizes left out of the a1, a2 fit because they would "
                        f"touch neighbouring coils: factors {skipped}")
    if len(widths) < 2:
        raise SettingsError("fewer than two usable winding pack sizes for the a1, a2 fit; "
                            "add smaller wp_scan_factors")
    return fit_peak_field(widths, peaks, _mean_current_ma(coils), coils.n_coils,
                          geometry.r_major, geometry.r_minor)


def _config(settings: Settings, s: _Stages) -> dict:
    """The values PROCESS reads, keyed as in the file."""
    mn = 1.0e-6
    g, coils, forces = s.geometry, s.coils, s.forces
    values: dict[str, str | int | float] = {
        "name": settings.name,
        "symmetry": int(s.eq.nfp),
        "coilspermodule": int(coils.n_coils // s.eq.nfp),
        "rmajor_ref": s.plasma.major_radius,
        "rminor_ref": s.plasma.minor_radius,
        "aspect_ref": s.plasma.aspect_ratio,
        "bt_ref": s.axis.mean_strength,
        "vol_plasma": s.plasma.volume,
        "plasma_surface": s.plasma.area,
        "coil_rmajor": g.r_major,
        "coil_rminor": g.r_minor,
        "coillength": cg.total_length(coils),
        "coilsurface": cg.coil_surface_area(coils),
        "maximal_coil_height": cg.coil_height(coils),
        "dmin": g.d_min,
        "max_portsize_width": cg.port_width(coils),
        "min_bend_radius": cg.min_bend_radius(coils),
        "min_plasma_coil_distance": g.d_pc,
        "derivative_min_LCFS_coils_dist": g.slope,
        "WP_area": s.wp.area,
        "WP_ratio": s.wp.ratio,
        "i0": _mean_current_ma(coils),
        "WP_bmax": s.pack.peak(),
        "a1": s.fit.a1,
        "a2": s.fit.a2,
        "inductance": s.energy.effective_inductance,
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
    config = {k: (v if isinstance(v, (str, int)) else float(v)) for k, v in values.items()}
    validate(config)
    return config


def _diagnostics(s: _Stages) -> dict:
    """Intermediate values worth keeping in the run report."""
    g, coils, fit = s.geometry, s.coils, s.fit
    station_peaks = s.pack.station_peaks()
    return {
        "plasma_as_reported_by_equilibrium": dict(s.eq.reported),
        "field_on_axis_mean_toroidal_component_T": s.axis.mean_toroidal,
        "magnetic_energy_GJ": s.energy.energy * 1.0e-9,
        "coil_classes": [
            {"representative_coil": int(i), "count": int(m),
             "current_MA": float(coils.currents[i]) * 1.0e-6,
             "length_m": float(coils.lengths[i]),
             "self_inductance_H": float(s.energy.self_inductance[c]),
             "peak_field_T": float(station_peaks[c].max()),
             "centering_force_MN": float(s.forces.centering[c]) * 1.0e-6}
            for c, (i, m) in enumerate(
                zip(g.classes.representatives, g.classes.multiplicity, strict=True))
        ],
        "peak_field_fit": {
            "radial_width_m": [float(x) for x in fit.radial_widths],
            "peak_field_T": [float(x) for x in fit.peak_fields],
            "fitted_T": [float(x) for x in fit.fitted],
            "max_relative_error": fit.max_relative_error,
        },
        "plasma_coil_distance_fit": {
            "minor_radius_m": [float(x) for x in g.fit_radii],
            "distance_m": [float(x) for x in g.fit_dists],
        },
        "coil_surface_lofted_between_filaments_m2": cg.lofted_surface_area(coils),
        "winding_pack_clearance_m": s.clearance,
        "mean_coil_spacing_at_coil_rmajor_m": 2.0 * np.pi * g.r_major / coils.n_coils,
    }


def _report(settings: Settings, s: _Stages, config: dict, warnings: list[str],
            runtime: float) -> dict:
    """Everything that is not read by PROCESS: methods, diagnostics, warnings."""
    return {
        "generator": f"stellarator_configprocessor {__version__}",
        "schema_version": SCHEMA_VERSION,
        "inputs": {"equilibrium": s.eq.source, "coils": s.coilset.source},
        "settings": {
            **{k: v for k, v in replace(settings).__dict__.items() if k != "tolerances"},
            "tolerances": dict(settings.tolerances.__dict__),
        },
        "values": [
            {"key": k.name, "value": config[k.name], "unit": k.unit,
             "description": k.description, "method": k.method}
            for k in KEYS
        ],
        "consistency": s.check.as_dict(),
        "scaling": s.scaling,
        "diagnostics": _diagnostics(s),
        "warnings": warnings,
        "runtime_s": round(runtime, 2),
    }


def _checked_inputs(eq: Equilibrium, coilset: CoilSet, settings: Settings,
                    warnings: list[str], log):
    """Sample the coils, check them against the equilibrium and rescale both."""
    coilset = _current_carrying(coilset, eq, warnings)
    log("sampling coils")
    coils = coilset.sampled(settings.n_points)

    log("checking that coils and equilibrium match")
    axis = axis_field(eq, coils)
    check = check_consistency(eq, coils, coilset.periods, settings.tolerances, axis=axis,
                              check_field=not settings.match_field)
    enforce(check, settings.allow_inconsistent)
    if not check.ok:
        warnings.extend("consistency check overridden: " + f for f in check.failures)

    scaling = _scaling_factors(eq, check, settings)
    if scaling["length_factor"] != 1.0 or scaling["total_current_factor"] != 1.0:
        eq = eq.scaled(scaling["length_factor"], scaling["field_factor"])
        coilset = coilset.scaled(scaling["length_factor"], scaling["total_current_factor"])
        coils = coilset.sampled(settings.n_points)
        axis = axis_field(eq, coils)
    if coils.n_coils % eq.nfp != 0:
        raise SettingsError(f"{coils.n_coils} coils cannot be divided into {eq.nfp} field periods")
    return eq, coilset, coils, axis, check, scaling


def run(eq: Equilibrium, coilset: CoilSet, settings: Settings, log=None) -> Result:
    """Compute the configuration for one equilibrium and coil set."""
    if not isinstance(eq, Equilibrium) or not isinstance(coilset, CoilSet):
        raise SettingsError("run needs an Equilibrium and a CoilSet")
    if not isinstance(settings, Settings):
        raise SettingsError("run needs a Settings object")
    log = log or _silent
    t_start = time.time()
    warnings: list[str] = list(coilset.warnings)

    eq, coilset, coils, axis, check, scaling = _checked_inputs(eq, coilset, settings, warnings, log)
    wp = WindingPack(settings.wp_radial, settings.wp_toroidal)

    log("plasma geometry")
    plasma = surface_geometry(eq)

    log("coil geometry")
    geometry = _coil_geometry(eq, coils, wp, settings, warnings)
    reps = geometry.classes.representatives

    log("field in the winding pack")
    frames = coil_frames(coils)
    clearance = pack_clearance(coils, wp, reps, frames)
    if clearance <= 0.0:
        warnings.append(
            f"winding packs of neighbouring coils overlap by {-clearance:.4g} m "
            f"(minimum distance between filaments: {geometry.d_min:.4g} m); the peak field "
            "and forces are not meaningful. Use smaller winding pack dimensions.")
    pack = winding_pack_field(coils, wp, reps, frames=frames, n_sub=settings.wp_filaments)

    log("forces")
    forces = coil_forces(coils, pack, frames, geometry.classes)

    log("peak field for several winding pack sizes")
    fit = _peak_field_fit(coils, wp, pack, frames, geometry, settings, warnings)

    log("magnetic energy")
    energy = magnetic_energy(coils, wp, geometry.classes)

    stages = _Stages(eq=eq, coilset=coilset, coils=coils, axis=axis, check=check,
                     scaling=scaling, plasma=plasma, geometry=geometry, wp=wp, pack=pack,
                     clearance=clearance, forces=forces, fit=fit, energy=energy)
    config = _config(settings, stages)
    report = _report(settings, stages, config, warnings, time.time() - t_start)
    return Result(config=config, report=report)
