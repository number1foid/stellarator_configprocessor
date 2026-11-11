"""Human-readable and JSON run reports."""

from __future__ import annotations

import json
from pathlib import Path


def _fmt(value) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_text(report: dict) -> str:
    """Plain-text report: every value, its unit and how it was obtained."""
    lines = [
        f"{report['generator']}  (schema {report['schema_version']})",
        f"equilibrium: {report['inputs']['equilibrium']}",
        f"coils:       {report['inputs']['coils']}",
        "",
        "Consistency of coils and equilibrium",
    ]
    c = report["consistency"]
    lines += [
        f"  field periods: equilibrium {c['nfp_equilibrium']}, coil file {c['nfp_coil_file']}, "
        f"{c['n_coils']} coils",
        f"  current linked by the magnetic axis / total current: {_fmt(c['linkage_ratio'])}",
        f"  field perpendicular to the axis / field: {_fmt(c['axis_misalignment'])}",
        f"  field on axis: coils {_fmt(c['b_axis_coils_T'])} T, equilibrium "
        f"{_fmt(c['b_axis_equilibrium_T'])} T (ratio "
        f"{_fmt(c['field_ratio_coils_over_equilibrium'])})",
        f"  result: {'passed' if c['passed'] else 'FAILED (overridden)'}",
    ]
    s = report["scaling"]
    if s["length_factor"] != 1.0 or s["total_current_factor"] != 1.0:
        lines += ["", "Scaling applied to the inputs",
                  f"  lengths x {_fmt(s['length_factor'])}, field x {_fmt(s['field_factor'])}, "
                  f"currents x {_fmt(s['total_current_factor'])}"]

    lines += ["", "Values written to the configuration file"]
    width = max(len(v["key"]) for v in report["values"])
    for v in report["values"]:
        unit = f" {v['unit']}" if v["unit"] else ""
        lines.append(f"  {v['key']:<{width}}  {_fmt(v['value']) + unit:<18}  {v['method']}")

    d = report["diagnostics"]
    lines += ["", "Diagnostics",
              f"  magnetic energy: {_fmt(d['magnetic_energy_GJ'])} GJ",
              "  coil classes (one representative coil each):"]
    for k in d["coil_classes"]:
        lines.append(
            f"    coil {k['representative_coil']:>3} x{k['count']:<3} current "
            f"{_fmt(k['current_MA'])} MA, length {_fmt(k['length_m'])} m, peak field "
            f"{_fmt(k['peak_field_T'])} T, centering force {_fmt(k['centering_force_MN'])} MN")
    f = d["peak_field_fit"]
    lines.append("  peak-field fit (radial width m: computed T / fitted T):")
    for w, b, fit in zip(f["radial_width_m"], f["peak_field_T"], f["fitted_T"], strict=True):
        lines.append(f"    {_fmt(w)}: {_fmt(b)} / {_fmt(fit)}")
    lines.append(f"    largest relative deviation of the fit: {_fmt(f['max_relative_error'])}")

    if report["warnings"]:
        lines += ["", "Warnings"] + [f"  - {w}" for w in report["warnings"]]
    lines += ["", f"runtime: {report['runtime_s']} s"]
    return "\n".join(lines) + "\n"


def write_report(path: str | Path, report: dict) -> Path:
    """Write the report as JSON (``.json``) or plain text (any other suffix)."""
    path = Path(path)
    if path.suffix.lower() == ".json":
        path.write_text(json.dumps(report, indent=2, default=str) + "\n")
    else:
        path.write_text(format_text(report))
    return path
