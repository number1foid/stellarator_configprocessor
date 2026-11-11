# stellarator_configprocessor

An open pre-processor for the stellarator model of the
[PROCESS](https://github.com/ukaea/process) fusion power plant systems code.

PROCESS can model any modular stellarator if it is given a configuration file
(`stella_conf.json`) of pre-calculated "effective parameters" for that design.
The tool that PROCESS's documentation names for producing this file,
pre-sPROCESS, is not publicly available. This project computes the same file
from two standard inputs:

- a plasma equilibrium: VMEC `wout` file (NetCDF);
- a coil set: MAKEGRID `coils.*` file with filament curves and currents.

The method follows Lion et al., *Nuclear Fusion* 61 (2021) 126021.

## Status

Development version, not yet released. Written against PROCESS 3.4.2.

- All 35 values PROCESS reads are written. Two of them (`epseff`,
  `neutron_peakfactor`) are passed through from the command line, not
  computed.
- PROCESS 3.4.2 accepts the generated file with `istell = 6`.
- For the as-built W7-X coils, the computed values agree with PROCESS's
  built-in W7-X set to within a few percent, and the differences are
  explained in [docs/w7x_comparison.md](docs/w7x_comparison.md).

## Install

Python 3.10 or newer. Dependencies: NumPy and netCDF4.

```
pip install -e .
```

## Use

```
stellarator-configprocessor wout_mydesign.nc coils.mydesign \
    --wp-radial 0.6 --wp-toroidal 0.5 \
    --epseff 0.015 --neutron-peakfactor 1.6 \
    --process-input myrun.IN.DAT --report myrun_report.txt
```

This writes `myrun.stella_conf.json` next to `myrun.IN.DAT`, which is where
PROCESS looks for it, and a report listing every value, its unit and how it
was obtained.

The four settings shown have no defaults: the results depend on the winding
pack dimensions, and the last two cannot be computed by this tool.

See [docs/user_guide.md](docs/user_guide.md) for all options and a worked
example from files to a PROCESS run.

## Documentation

| Document | Content |
|---|---|
| [docs/user_guide.md](docs/user_guide.md) | Options, inputs, worked W7-X example, Python interface, limitations |
| [docs/schema.md](docs/schema.md) | Every key of `stella_conf.json`: unit, how PROCESS uses it, how it is computed; reference values |
| [docs/w7x_comparison.md](docs/w7x_comparison.md) | PROCESS run with the generated W7-X file against the built-in W7-X set |
| [IMPLEMENTATION_OUTLINE.md](IMPLEMENTATION_OUTLINE.md) | Project plan and design decisions |

## Verification

```
pip install -e ".[test]"
python scripts/fetch_w7x_data.py    # optional: public W7-X files for the reference tests
pytest
```

| Quantity | Checked against | Agreement |
|---|---|---|
| Plasma volume, area, radii | circular torus (analytic); VMEC's own values for W7-X | 1e-12; 1e-10 |
| Field of a filament | circular loop and straight wire (analytic) | 1e-4 |
| Field on axis of W7-X coils | the W7-X vacuum equilibrium | 5e-6 |
| Field inside the conductor | straight rectangular bar (closed form) | 2%, limited by loop curvature |
| Self-inductance | circular loop formula; volume integral of B^2 | 0.2%; 2% |
| Mutual inductance | coaxial loops (elliptic integrals) | 1e-4 |
| Hoop force | circular loop formula | 2% |
| Net coil force | derivative of the magnetic energy (virtual work) | 0.6% |
| Full W7-X result | PROCESS built-in W7X-50 values | see docs/w7x_comparison.md |

Tests that need the W7-X files are skipped when the files are absent. The
SIMSOPT cross-check is skipped when SIMSOPT is not installed.

## Limitations

- Modular coils only: every coil must close poloidally around the plasma.
- PROCESS scales one coil current. Coil sets with unequal currents are
  reduced to their mean absolute current, with a warning.
- Coils with zero current are ignored.
- `epseff` and `neutron_peakfactor` are inputs, not results.
- Input formats: VMEC `wout` and MAKEGRID only.

## Licence

MIT. See [LICENSE](LICENSE).
