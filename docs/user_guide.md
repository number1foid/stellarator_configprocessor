# User guide

## 1. What the tool does

It reads a plasma equilibrium and a coil set, computes the 35 values PROCESS
needs to model that stellarator, and writes them to the JSON file PROCESS
reads when `istell = 6`. It also writes a report of every value and how it was
obtained. The values themselves are defined in [schema.md](schema.md).

## 2. Install

```
pip install -e .
```

Python 3.10 or newer, NumPy and netCDF4. The command
`stellarator-configprocessor` is then available.

## 3. Inputs

### Equilibrium

A VMEC `wout` file in NetCDF format, from VMEC2000 or VMEC++. The last closed
flux surface defines the plasma; inner surfaces are used for one derivative.

### Coils

A MAKEGRID `coils.*` file.

- It must list **every coil of the full torus**, not only one field period.
- Each coil must close poloidally around the plasma (modular coils).
- The fourth column is the coil current in ampere-turns. If the file holds
  unit currents, give the real currents per coil group with
  `--group-current GROUP=FACTOR` (the VMEC `extcur` values), or use
  `--match-field`.
- Coils with zero current are ignored.

### Required settings

These four have no default.

| Option | Meaning |
|---|---|
| `--wp-radial M` | Width of the winding pack in the coil's in-plane normal direction, in metres |
| `--wp-toroidal M` | Width of the winding pack in the toroidal direction, in metres |
| `--epseff X` | Effective helical ripple of the configuration |
| `--neutron-peakfactor X` | Peak over mean neutron wall load |

The winding pack is the current-carrying rectangle of one coil. The peak
field, the forces and the inductance are computed for it and depend on it.
PROCESS later resizes the winding pack itself and rescales those values, so
choose dimensions that are realistic for the reference machine and smaller
than the smallest gap between coils (`dmin` in the report).

`epseff` and `neutron_peakfactor` are not computed. `epseff` comes from a
neoclassical code such as NEO; PROCESS's built-in configurations use 0.015.
`neutron_peakfactor` comes from a wall-load calculation; the built-in
configurations use 1.6.

## 4. The reference point

By default the reference point is the machine as described by the input
files. Two options rescale it:

| Option | Effect |
|---|---|
| `--rmajor-ref M` | Scale all lengths so the plasma major radius is `M` metres. Currents are scaled to keep the field. |
| `--bt-ref T` | Scale the coil currents so the mean field on the magnetic axis is `T` tesla. |

The winding pack dimensions always refer to the rescaled machine.

PROCESS scales the configuration to any size and field itself, so rescaling is
only needed to state the reference point in convenient numbers.

## 5. The consistency check

Before computing anything the tool checks that coils and equilibrium describe
the same machine:

1. the number of field periods agrees and divides the coil count;
2. the magnetic axis links the full coil current (every coil encloses it);
3. the coil field is tangent to the equilibrium's magnetic axis;
4. the mean field on axis from the coils equals the equilibrium's.

If a test fails the tool stops with exit code 3 and says which one.

| Option | Use |
|---|---|
| `--match-field` | The coil file has the right shape but another current scale. All currents are scaled to the equilibrium's field and test 4 is skipped. |
| `--allow-inconsistent` | Continue regardless. The failures are recorded as warnings in the report. |
| `--tol-field`, `--tol-axis`, `--tol-linkage` | Change the limits (defaults 0.03, 0.02, 0.01). |

A finite-pressure equilibrium differs from the vacuum field of its coils by a
few percent on axis, which can exceed the default limits.

Exit codes: 0 success, 2 unusable input or settings, 3 inconsistent inputs.

## 6. Output

### Configuration file

| Option | Result |
|---|---|
| `-o FILE` | Write to `FILE` (default `stella_conf.json`) |
| `--process-input RUN.IN.DAT` | Write to `RUN.stella_conf.json` next to that input file |

**PROCESS does not look for a file called `stella_conf.json` unless the input
file is called `IN.DAT`.** For `myrun.IN.DAT` it opens
`myrun.stella_conf.json` in the same folder. `--process-input` applies that
rule. PROCESS reports no error if the file is found but a key is missing; it
uses zero. This tool always writes all keys.

The file also holds two entries starting with an underscore, the schema
version and the generating program. PROCESS ignores them.

### Report

`--report FILE` writes a report: JSON if the name ends in `.json`, plain text
otherwise. It contains the settings, the consistency check, every value with
unit and method, per-coil diagnostics, the peak-field fit and any warnings.
The same text is printed unless `-q` is given.

## 7. Worked example: W7-X

The files are public but not part of this repository. One script downloads
them and builds the coil file:

```
python scripts/fetch_w7x_data.py
```

This creates, in `data/local/`:

- `wout_w7x_standard.nc`: vacuum equilibrium of the W7-X standard configuration;
- `coils.w7x_standard`: the 50 as-built non-planar coils, 15 kA x 108 turns each;
- `epseff_w7x_standard.txt`: effective ripple profile of that configuration.

Run the tool, stating the reference point at 3 T:

```
stellarator-configprocessor data/local/wout_w7x_standard.nc data/local/coils.w7x_standard \
    --wp-radial 0.18 --wp-toroidal 0.15 \
    --epseff 0.0095 --neutron-peakfactor 1.6 \
    --bt-ref 3.0 --name "W7-X standard" \
    --process-input runs/w7x.IN.DAT --report runs/w7x_report.txt
```

It takes about one minute. The consistency check prints:

```
  field on axis: coils 2.78472 T, equilibrium 2.78473 T (ratio 0.999995)
  result: passed
```

and the main results are:

| Key | Value |
|---|---|
| `rmajor_ref`, `rminor_ref` | 5.512 m, 0.5127 m |
| `vol_plasma`, `plasma_surface` | 28.60 m^3, 136.7 m^2 |
| `i0` | 1.745 MA |
| `WP_bmax` | 6.29 T |
| `a1`, `a2` | 0.655, 0.0308 |
| `inductance` | 4.18e-4 H |
| `dmin`, `min_plasma_coil_distance` | 0.279 m, 0.394 m |
| `centering_force_max_MN` | 2.91 MN |

Then run PROCESS on an input file `runs/w7x.IN.DAT` that sets `istell = 6`:

```
process -i runs/w7x.IN.DAT
```

A stellarator input file to start from is
`tests/regression/input_files/stellarator_helias.IN.DAT` in the PROCESS
repository. [w7x_comparison.md](w7x_comparison.md) shows the result of such a
run against PROCESS's built-in W7-X configuration.

On native Windows, PROCESS 3.4.2 needs three workarounds that are unrelated to
this tool. `python scripts/run_process_windows.py -i runs/w7x.IN.DAT` applies
them; the script explains each.

## 8. Resolution and run time

| Option | Default | Effect |
|---|---|---|
| `--points N` | 160 | Samples per coil. More samples resolve tighter bends. |
| `--wp-filaments NR NT` | 4 4 | Filaments across the winding pack. The W7-X peak field changes by 0.1% between 4 x 4 and 6 x 6. |

Run time grows with the product of the two and with the number of distinct
coil shapes. Coils that are symmetric images of each other are computed once.

## 9. Python interface

```python
from stellarator_configprocessor.io.equilibrium import read_wout
from stellarator_configprocessor.io.coils import read_makegrid
from stellarator_configprocessor.pipeline import Settings, run
from stellarator_configprocessor.schema import write_config
from stellarator_configprocessor.report import write_report

eq = read_wout("wout_mydesign.nc")
coils = read_makegrid("coils.mydesign")
settings = Settings(wp_radial=0.6, wp_toroidal=0.5, epseff=0.015, neutron_peakfactor=1.6)

result = run(eq, coils, settings)
write_config("myrun.stella_conf.json", result.config)
write_report("myrun_report.json", result.report)
```

`result.config` is a dictionary with the 35 keys. Each calculation is also
available on its own, for example
`stellarator_configprocessor.geometry.plasma.surface_geometry`,
`stellarator_configprocessor.winding_pack.winding_pack_field` and
`stellarator_configprocessor.energy.magnetic_energy`.

If a coil file holds only the coils of one half period,
`stellarator_configprocessor.io.coils.expand_symmetry` builds the full set.

## 10. Limitations and known issues

- **Modular coils only.** Helical or toroidally continuous coils are rejected,
  because the winding pack orientation is undefined for them.
- **One current.** PROCESS scales a single coil current. If coil currents
  differ, `i0` is their mean absolute value and a warning is issued. Planar
  coils with zero current, as in the W7-X standard configuration, are ignored.
- **`max_portsize_width`, `coilsurface`, `coil_rminor`** have no published
  definition. The definitions used are stated in [schema.md](schema.md).
- **`derivative_min_LCFS_coils_dist`** is the geometric derivative. PROCESS's
  radial-build formula applies it in a way that changes the available space
  by `(1 - derivative)` per unit of minor radius; see open point 4 in
  [schema.md](schema.md). This matters only if the aspect ratio is an
  iteration variable.
- **Winding pack wider than the coil gap.** If `--wp-toroidal` is not smaller
  than `dmin`, neighbouring winding packs overlap. The tool warns and the
  peak field and forces are then not meaningful.
