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
| [docs/user_guide.md](docs/user_guide.md) | Options, inputs, worked W7-X example, Python interface, usage limits |
| [docs/schema.md](docs/schema.md) | Every key of `stella_conf.json`: unit, how PROCESS uses it, how it is computed; reference values |
| [docs/w7x_comparison.md](docs/w7x_comparison.md) | PROCESS run with the generated W7-X file against the built-in W7-X set |
| [docs/hsx_check.md](docs/hsx_check.md) | Second machine (HSX): generality check and what it exposed |
| [IMPLEMENTATION_OUTLINE.md](IMPLEMENTATION_OUTLINE.md) | Project plan and design decisions |

## Verification

```
pip install -e ".[test]"
python scripts/fetch_reference_data.py    # optional: public W7-X and HSX files
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
| Field on axis, coil lengths (W7-X) | SIMSOPT, in CI | 1e-3; 1e-4 |
| Full W7-X result | PROCESS built-in W7X-50 values | see docs/w7x_comparison.md |
| Full HSX result | consistency with its equilibrium; plausibility | see docs/hsx_check.md |

Tests that need the W7-X or HSX files are skipped when the files are absent. The
SIMSOPT cross-check is skipped when SIMSOPT is not installed.

## Known issues

### Not yet verified

- **No independent check of peak field, energy and forces on a real machine.**
  These rest on analytic test cases and on agreement with PROCESS's built-in
  W7-X values, some of which are rounded or placeholders. The SIMSOPT
  cross-check covers only the filament field and the coil lengths.
- **Two unexplained differences from the built-in W7-X values.** `a2` is 25%
  above the built-in value and `max_lateral_force_density` 10% above. The fit
  reproduces the fields this tool computes to 0.06%; which side is right is
  not established.
- **Three definitions are informed guesses.** `coil_rminor`, `coilsurface`
  and `max_portsize_width` have no published formula. The definitions used
  were chosen to match the available reference values; only the authors of
  the original tool can confirm them. See [docs/schema.md](docs/schema.md).
- **No complete optimised PROCESS design.** PROCESS reads and evaluates the
  generated files. A reactor-scale optimisation found no feasible plant for
  W7-X, with either the generated or the built-in configuration, and no
  optimisation has been run for another machine.
- **Paths without a real-data test.** Coil sets with unequal currents,
  equilibria without stellarator symmetry and finite-pressure equilibria are
  handled in the code but have only been run on synthetic cases or not at all.
- **Accuracy on HSX is unknown.** The second machine shows the tool is
  general and its results plausible; there are no reference values for it.

### Not done

- **`epseff` and `neutron_peakfactor` are not computed.** They are required
  inputs, passed through to the file. Until they are computed, the tool does
  not fully replace pre-sPROCESS.
- **PROCESS is not part of the automated tests.** The schema is tied to
  PROCESS 3.4.2 by the documents only, so a schema change upstream would not
  be caught automatically.
- **The implementation outline is partly out of date.** Its "not confirmed"
  list and its Phase 2 step on boundary coefficients predate what the
  PROCESS source showed; [docs/schema.md](docs/schema.md) is authoritative.

### By design

- **Modular coils only.** Every coil must close poloidally around the plasma.
- **One coil current.** PROCESS scales a single current, so coil sets with
  unequal currents are reduced to their mean absolute current, with a
  warning. Coils with zero current are ignored.
- **Input formats.** VMEC `wout` and MAKEGRID only.

### In PROCESS itself

Found while reading and running PROCESS 3.4.2; details in
[docs/schema.md](docs/schema.md), section 6.

- The documentation swaps `istell = 4` and `istell = 5`: in the code, 4 is
  the 30-coil W7-X set and 5 the 50-coil set.
- A key missing from the configuration file is silently read as zero.
- The centering forces are scaled with the coil-length ratio inverted
  relative to eq. (43) of Lion et al. 2021.
- On native Windows a run fails without the workarounds in
  `scripts/run_process_windows.py`.

## Licence

MIT. See [LICENSE](LICENSE).
