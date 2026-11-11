# User guide

## 1. What the tool does

It reads a plasma equilibrium and a coil set, computes the 35 values PROCESS
needs to model that stellarator, and writes them to the JSON file PROCESS
reads when `istell = 6`. It also writes a report of every value and how it was
obtained. The values themselves are defined in [schema.md](schema.md).

## 2. Install

The program is written in Rust and has no other dependencies. Install the
Rust toolchain from <https://rustup.rs> (version 1.82 or newer), then build:

```
cargo build --release
```

The program is then `target/release/stellarator-configprocessor`
(`stellarator-configprocessor.exe` on Windows). `cargo install --path .`
copies it to Cargo's `bin` folder, so that it can be run from anywhere.
`stellarator-configprocessor --help` lists all options.

## 3. Inputs

### Equilibrium

A VMEC `wout` file in the NetCDF classic format (CDF-1 or CDF-2, "64-bit
offset"), which is what VMEC2000 writes by default. The last closed flux
surface defines the plasma; inner surfaces are used for one derivative.

Files in the HDF5-based NetCDF-4 format are rejected with a message. Convert
them with the NetCDF tools: `nccopy -k classic wout_in.nc wout_out.nc`.

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

The output files are opened before the calculation starts. If the run fails,
an output file that already existed is left as it was; one that did not exist
is left behind empty.

## 7. Worked example: W7-X

The files are public but not part of this repository. One script (Python
with NumPy) downloads them and builds the coil file:

```
python scripts/fetch_reference_data.py
```

This creates, in `data/local/`:

- `wout_w7x_standard.nc`: vacuum equilibrium of the W7-X standard configuration;
- `coils.w7x_standard`: the 50 as-built non-planar coils, 15 kA x 108 turns each;
- `epseff_w7x_standard.txt`: effective ripple profile of that configuration;
- `wout_hsx_qhs.nc` and `coils.hsx_qhs`: a second machine, see [hsx_check.md](hsx_check.md).

Run the tool, stating the reference point at 3 T:

```
stellarator-configprocessor data/local/wout_w7x_standard.nc data/local/coils.w7x_standard \
    --wp-radial 0.18 --wp-toroidal 0.15 \
    --epseff 0.0095 --neutron-peakfactor 1.6 \
    --bt-ref 3.0 --name "W7-X standard" \
    --process-input runs/w7x.IN.DAT --report runs/w7x_report.txt
```

It takes about ten seconds. The consistency check prints:

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
| `a1`, `a2` | 0.641, 0.0311 |
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
The W7-X example takes about 9 seconds on one processor core; the program
uses a single core.

All sizes have fixed upper bounds, listed in `src/limits.rs`: for example 512
coils, 1024 samples per coil and 12 filaments across each side of the winding
pack. An input beyond a bound is rejected with a message.

## 9. Use as a library

The crate can be used from other Rust code. The command-line program does
no more than this:

```rust
use stellarator_configprocessor::equilibrium::read_wout;
use stellarator_configprocessor::error::Error;
use stellarator_configprocessor::makegrid::{read_makegrid, GroupFactors};
use stellarator_configprocessor::pipeline::{default_settings, new_workspace, run};
use stellarator_configprocessor::schema::write_config;

fn main() -> Result<(), Error> {
    let bytes = std::fs::read("wout_mydesign.nc")?;
    let text = std::fs::read_to_string("coils.mydesign")?;
    let no_factors = GroupFactors { count: 0, group: [0.0; 32], factor: [1.0; 32] };

    let mut eq = read_wout(&bytes)?;
    let mut coils = read_makegrid(&text, &no_factors)?;
    // Winding pack (radial, toroidal) [m], epseff, neutron peaking factor.
    let settings = default_settings((0.6, 0.5), 0.015, 1.6);

    // All memory is allocated here; `run` allocates nothing.
    let mut work = new_workspace(&eq, &coils, &settings)?;
    let outcome = run(&mut eq, &mut coils, &settings, &mut work, &mut std::io::sink())?;

    let mut file = std::fs::File::create("myrun.stella_conf.json")?;
    write_config(&mut file, &outcome.config, "my design")
}
```

`outcome.config` holds the 35 values in the order of `schema::KEYS`. Each
calculation is also available on its own, for example
`plasma::surface_geometry`, `winding_pack::reference_pack_field` and
`energy::magnetic_energy`. `cargo doc --open` shows the documentation of
every function.

The coil file must list the coils of the full torus. The Python script
`scripts/fetch_reference_data.py` shows how a full set is built from the
coils of one half period.

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
- **Winding pack too large for the coil gap.** If the winding packs of
  neighbouring coils overlap, the tool warns and states the overlap; the peak
  field and forces are then not meaningful. The `a1`, `a2` fit uses packs from
  half the reference size up to the reference size.
- **NetCDF classic only.** See section 3.
- **Coding rules.** The code follows the NASA/JPL "Power of Ten" rules; what
  that means for users (fixed size limits, no memory allocation during the
  calculation, one processor core) is described in
  [power_of_ten.md](power_of_ten.md).
