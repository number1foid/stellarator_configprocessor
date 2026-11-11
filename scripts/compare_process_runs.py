"""Compare two PROCESS runs that differ only in the stellarator configuration.

Typical use: one run with a built-in configuration (for example ``istell = 5``)
and one with a file written by this tool (``istell = 6``).

``prepare`` turns a PROCESS stellarator input into an evaluation-only input:
PROCESS then evaluates its models once at the given major radius and field,
solving only the power balance for the confinement factor, so every output
difference traces back to the configuration values.

Usage::

    python scripts/compare_process_runs.py prepare base.IN.DAT out.IN.DAT --istell 5
    python scripts/compare_process_runs.py compare a.MFILE.DAT b.MFILE.DAT
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

#: Outputs that depend on the stellarator configuration, in report order.
KEY_OUTPUTS = [
    "n_tf_coils", "rmajor", "rminor", "aspect", "vol_plasma", "a_plasma_surface",
    "b_plasma_toroidal_on_axis", "coil_r", "coil_a", "c_tf_total/n_tf_coils",
    "b_tf_inboard_peak_symmetric", "e_tf_magnetic_stored_total_gj", "inductance",
    "len_tf_coil", "dr_tf_wp_with_insulation", "dx_tf_wp_primary_toroidal", "ap", "j_tf_wp",
    "toroidalgap", "coilcoilgap", "available_radial_space", "required_radial_space",
    "min_bend_radius", "max_force_density", "max_force_density_Mnm",
    "max_lateral_force_density", "max_radial_force_density", "centering_force_max_MN",
    "centering_force_min_MN", "centering_force_avg_MN", "v_tf_coil_dump_quench_kv",
    "vporttmax", "hporttmax", "m_tf_coils_total", "aintmass", "hfact",
    "p_plant_electric_net_mw", "coe",
]


def prepare(source: Path, target: Path, istell: int, rmajor: float | None,
            bt: float | None) -> None:
    """Write an evaluation-only copy of a PROCESS stellarator input file."""
    out = []
    for line in source.read_text(encoding="utf-8").splitlines():
        stripped = line.strip().lower()
        key = stripped.split("=")[0].strip() if "=" in stripped else ""
        if key == "icc" and not re.match(r"icc\s*=\s*2\b", stripped):
            continue  # keep only the global power balance
        if key == "ixc" and not re.match(r"ixc\s*=\s*10\b", stripped):
            continue  # keep only the confinement factor
        text = line
        if key == "n_equality_constraints":
            text = "n_equality_constraints = 1"
        elif key == "i_process_run_mode":
            text = "i_process_run_mode = -2"
        elif key == "istell":
            text = f"istell = {istell}"
        elif key == "rmajor" and rmajor is not None:
            text = f"rmajor = {rmajor}"
        elif key == "b_plasma_toroidal_on_axis" and bt is not None:
            text = f"b_plasma_toroidal_on_axis = {bt}"
        out.append(text)
    target.write_text("\n".join(out) + "\n", encoding="utf-8")


def read_mfile(path: Path) -> dict:
    """Variable name -> value for every numeric entry of a PROCESS MFILE."""
    values = {}
    # The variable name is the last parenthesised token before the value.
    pattern = re.compile(r".*\(([^()\s]+)\)_*\s+(\S+)")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if not match:
            continue
        try:
            values[match.group(1)] = float(match.group(2))
        except ValueError:
            continue
    return values


def compare(file_a: Path, file_b: Path, threshold: float) -> str:
    a, b = read_mfile(file_a), read_mfile(file_b)

    def row(name):
        va, vb = a.get(name), b.get(name)
        if va is None or vb is None:
            return f"  {name:<34} {'missing':>14} {'missing':>14}"
        rel = (vb - va) / abs(va) if va != 0 else float("inf") if vb != 0 else 0.0
        return f"  {name:<34} {va:>14.6g} {vb:>14.6g} {rel:>+9.2%}"

    lines = [f"A: {file_a}", f"B: {file_b}", "",
             f"  {'variable':<34} {'A':>14} {'B':>14} {'B vs A':>9}"]
    lines += [row(name) for name in KEY_OUTPUTS if name in a or name in b]
    others = []
    for name in sorted(set(a) & set(b)):
        if name in KEY_OUTPUTS:
            continue
        va, vb = a[name], b[name]
        scale = max(abs(va), abs(vb))
        if scale > 0 and abs(vb - va) / scale > threshold:
            others.append(name)
    lines += ["", f"{len(others)} other outputs differ by more than {threshold:.0%} "
                  f"(of {len(set(a) & set(b))} common numeric outputs)."]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="write an evaluation-only input file")
    p.add_argument("source", type=Path)
    p.add_argument("target", type=Path)
    p.add_argument("--istell", type=int, required=True)
    p.add_argument("--rmajor", type=float)
    p.add_argument("--bt", type=float)
    c = sub.add_parser("compare", help="compare two MFILE.DAT files")
    c.add_argument("mfile_a", type=Path)
    c.add_argument("mfile_b", type=Path)
    c.add_argument("--threshold", type=float, default=0.05)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.source, args.target, args.istell, args.rmajor, args.bt)
        print(f"wrote {args.target}")
    else:
        print(compare(args.mfile_a, args.mfile_b, args.threshold))


if __name__ == "__main__":
    main()
