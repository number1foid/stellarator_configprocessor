"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import SCHEMA_VERSION, __version__
from .consistency import ConsistencyError, Tolerances
from .io.coils import CoilFileError, read_makegrid
from .io.equilibrium import EquilibriumError, read_wout
from .pipeline import Settings, SettingsError, run
from .report import format_text, write_report
from .schema import SchemaError, process_file_name, write_config

EXIT_OK, EXIT_INPUT, EXIT_INCONSISTENT = 0, 2, 3


def _group_current(text: str) -> tuple[int, float]:
    try:
        group, factor = text.split("=", 1)
        return int(group), float(factor)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"expected GROUP=FACTOR, for example 1=15000, got '{text}'") from exc


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stellarator-configprocessor",
        description="Compute the stellarator configuration file (stella_conf.json) that the "
                    "PROCESS systems code reads with istell = 6, from a VMEC equilibrium and "
                    "a MAKEGRID coil file.",
        epilog=f"Writes schema {SCHEMA_VERSION}. See docs/user_guide.md.")
    p.add_argument("wout", type=Path, help="VMEC wout file (NetCDF)")
    p.add_argument("coils", type=Path, help="MAKEGRID coils file listing all coils")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    req = p.add_argument_group("required settings (no defaults)")
    req.add_argument("--wp-radial", type=float, required=True, metavar="M",
                     help="winding pack width in the coil's in-plane normal direction [m]")
    req.add_argument("--wp-toroidal", type=float, required=True, metavar="M",
                     help="winding pack width in the toroidal direction [m]")
    req.add_argument("--epseff", type=float, required=True,
                     help="effective helical ripple (not computed; passed to PROCESS)")
    req.add_argument("--neutron-peakfactor", type=float, required=True,
                     help="neutron wall load peaking factor (not computed; passed to PROCESS)")

    out = p.add_argument_group("output")
    target = out.add_mutually_exclusive_group()
    target.add_argument("-o", "--output", type=Path, default=None,
                        help="configuration file to write (default: stella_conf.json)")
    target.add_argument("--process-input", type=Path, metavar="IN.DAT",
                        help="PROCESS input file; the configuration is written next to it "
                             "under the name PROCESS expects")
    out.add_argument("--report", type=Path,
                     help="report file; .json gives JSON, anything else plain text")
    out.add_argument("--name", default="stellarator_configprocessor",
                     help="configuration name stored in the file")
    out.add_argument("-q", "--quiet", action="store_true", help="print nothing on success")

    ref = p.add_argument_group("reference point")
    ref.add_argument("--rmajor-ref", type=float, metavar="M",
                     help="rescale the machine to this plasma major radius [m]")
    ref.add_argument("--bt-ref", type=float, metavar="T",
                     help="rescale coil currents to this mean field on axis [T]")
    ref.add_argument("--group-current", type=_group_current, action="append", default=[],
                     metavar="GROUP=FACTOR",
                     help="multiply the currents of a coil group by FACTOR (VMEC extcur); "
                          "repeatable")

    chk = p.add_argument_group("consistency check")
    chk.add_argument("--match-field", action="store_true",
                     help="scale all coil currents so the field on axis equals the "
                          "equilibrium's, instead of checking it")
    chk.add_argument("--allow-inconsistent", action="store_true",
                     help="continue even if coils and equilibrium appear inconsistent")
    defaults = Tolerances()
    chk.add_argument("--tol-field", type=float, default=defaults.field,
                     help=f"relative tolerance on the field on axis (default {defaults.field})")
    chk.add_argument("--tol-axis", type=float, default=defaults.axis_misalignment,
                     help="tolerance on the field component perpendicular to the axis "
                          f"(default {defaults.axis_misalignment})")
    chk.add_argument("--tol-linkage", type=float, default=defaults.linkage,
                     help=f"tolerance on the linked current (default {defaults.linkage})")

    num = p.add_argument_group("resolution")
    num.add_argument("--points", type=int, default=160, help="samples per coil (default 160)")
    num.add_argument("--wp-filaments", type=int, nargs=2, default=(4, 4), metavar=("NR", "NT"),
                     help="filaments across the winding pack (default 4 4)")
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    def say(message: str) -> None:
        if not args.quiet:
            print(message, flush=True)

    try:
        settings = Settings(
            wp_radial=args.wp_radial, wp_toroidal=args.wp_toroidal, epseff=args.epseff,
            neutron_peakfactor=args.neutron_peakfactor, name=args.name,
            rmajor_ref=args.rmajor_ref, bt_ref=args.bt_ref, match_field=args.match_field,
            allow_inconsistent=args.allow_inconsistent,
            tolerances=Tolerances(linkage=args.tol_linkage, axis_misalignment=args.tol_axis,
                                  field=args.tol_field),
            n_points=args.points, wp_filaments=tuple(args.wp_filaments))
        if args.process_input is not None:
            output = process_file_name(args.process_input)
        else:
            output = args.output or Path("stella_conf.json")
        eq = read_wout(args.wout)
        coilset = read_makegrid(args.coils, group_currents=dict(args.group_current) or None)
        result = run(eq, coilset, settings, log=lambda m: say(f"... {m}"))
        write_config(output, result.config)
    except ConsistencyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("Use --allow-inconsistent to continue anyway, or --match-field if only the "
              "current scale differs.", file=sys.stderr)
        return EXIT_INCONSISTENT
    except (EquilibriumError, CoilFileError, SettingsError, SchemaError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_INPUT

    say("")
    say(format_text(result.report))
    say(f"wrote {output}")
    if args.report:
        write_report(args.report, result.report)
        say(f"wrote {args.report}")
    for warning in result.report["warnings"]:
        print(f"warning: {warning}", file=sys.stderr)
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
