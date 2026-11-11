# W7-X check: PROCESS with the generated file against the built-in W7X-50

This is the Phase 5 acceptance check: PROCESS must accept the generated
configuration file, and a W7-X run must land close to PROCESS's built-in
W7-X result, with every difference explained.

## Result

- PROCESS 3.4.2 reads the generated file with `istell = 6` and completes.
- At the same design point, the coil outputs agree with the built-in W7X-50
  run (`istell = 5`) to within about 6%, apart from four groups of
  differences that are each traced to a cause below.

## Setup

| | |
|---|---|
| PROCESS | 3.4.2, installed from `github.com/ukaea/process` at tag `v3.4.2` |
| Input file | PROCESS's own `tests/regression/input_files/stellarator_helias.IN.DAT` |
| Mode | Evaluation only (`i_process_run_mode = -2`): models are evaluated once at the given size and field; only the power balance is solved, for the confinement factor |
| Design point | `rmajor = 20 m`, `b_plasma_toroidal_on_axis = 5.5 T` |
| Run A | `istell = 5`, built-in W7X-50 |
| Run B | `istell = 6`, file from this tool: as-built W7-X coils, standard-configuration vacuum equilibrium, `bt_ref = 3.0 T`, winding pack 0.18 m x 0.15 m, `epseff = 0.0095`, `neutron_peakfactor = 1.6` |

Evaluation mode is used because the optimisation run of the same input finds
no feasible 1 GW plant for either W7-X configuration (run A or run B). The
end state of a failed optimisation is not a meaningful comparison.

## Outputs

| PROCESS output | A: built-in | B: this tool | B vs A | Cause |
|---|---|---|---|---|
| Number of coils | 50 | 50 | 0% | |
| Coil major radius (m) | 20.44 | 20.38 | -0.3% | |
| Coil minor radius (m) | 4.291 | 4.152 | -3.2% | definition of `coil_rminor` |
| Current per coil (MA) | 11.63 | 11.61 | -0.2% | |
| Stored magnetic energy (GJ) | 101.5 | 102.3 | +0.8% | |
| Inductance (H) | 1.500e-3 | 1.517e-3 | +1.2% | |
| Mean coil length (m) | 30.59 | 30.90 | +1.0% | |
| Minimum coil distance (m) | 1.018 | 1.013 | -0.5% | |
| Minimum bending radius (m) | 1.510 | 1.542 | +2.1% | |
| Available radial space (m) | 1.410 | 1.429 | +1.3% | |
| Centering force, max (MN) | 10.26 | 10.59 | +3.3% | |
| Centering force, min (MN) | -4.43 | -4.68 | +5.8% | |
| Centering force, mean (MN) | 5.52 | 5.65 | +2.3% | |
| Lateral force density (MN/m^3) | 171.9 | 166.9 | -2.9% | |
| Intercoil structure mass (kg) | 4.18e6 | 4.07e6 | -2.6% | |
| Peak field on conductor (T) | 11.79 | 12.49 | +6.0% | 1 |
| Winding pack area (m^2) | 0.228 | 0.274 | +20% | 1 |
| Winding pack current density (A/m^2) | 5.10e7 | 4.24e7 | -17% | 1 |
| Coil-coil gap left (m) | 0.442 | 0.395 | -11% | 1 |
| Quench voltage (kV) | 0.725 | 0.880 | +21% | 1 |
| Total coil mass (kg) | 3.60e6 | 4.16e6 | +16% | 1 |
| Radial force density (MN/m^3) | 218.5 | 191.3 | -12% | 1 |
| Mean force density (MN/m^3) | 369.1 | 199.9 | -46% | 2 |
| Force per unit length (MN/m) | 1223 | 54.8 | -96% | 2 |
| Vertical port width (m) | 0.436 | 1.115 | +155% | 3 |
| Horizontal port width (m) | 0.873 | 2.230 | +155% | 3 |
| Plasma minor radius (m) | 1.786 | 1.860 | +4.2% | 4 |
| Aspect ratio | 11.2 | 10.75 | -4.0% | 4 |
| Confinement factor needed | 1.206 | 1.145 | -5.0% | 4 |
| Net electric power (MW) | 779 | 876 | +12% | 4 |
| Cost of electricity | 113.4 | 103.3 | -8.9% | 4 |

## Causes of the larger differences

1. **Size dependence of the peak field (`a2`).** Both configurations give the
   same peak field at the reference winding pack (6.29 T against 6.3 T). They
   differ in how the peak field grows as the winding pack shrinks relative to
   the coil: `a2` is 0.0311 here and 0.025 built in. At reactor scale PROCESS
   works at a smaller relative winding pack than the reference, so it finds a
   6% higher peak field and compensates with a 20% larger winding pack. Coil
   mass, current density, quench voltage, the remaining coil gap and the
   radial force density follow from that. The fit used here reproduces the
   peak fields computed for five winding pack sizes to 0.06%. The built-in
   pair, applied to the smallest of those sizes, gives a field 13% below the
   computed one.
2. **Built-in force values are not from one calculation.** The built-in
   `max_force_density_MNm` is 98.0 for all five configurations, including the
   W7-X-sized ones, where the computed value is 4.1 MN/m. The built-in
   `max_force_density` (250) is not the ratio of that value to the winding
   pack area either. The values written by this tool satisfy
   `max_force_density = max_force_density_MNm / WP_area`, as the example file
   shipped with PROCESS does.
3. **Port width definition.** No formula is published. This tool writes the
   largest centre-to-centre distance between neighbouring coils (0.77 m),
   which matches the example file's convention. The built-in value is 0.3 m.
4. **Different reference equilibrium.** The built-in plasma (26.4 m^3, aspect
   ratio 11.2) is smaller than the standard-configuration equilibrium used
   here (28.6 m^3, aspect ratio 10.75). At the same major radius the plasma is
   8% larger, which raises the fusion power and lowers the cost of
   electricity. The plasma values written here agree with those VMEC reports
   for the same file to 1e-10.

## Reproducing

```
python scripts/fetch_reference_data.py

stellarator-configprocessor data/local/wout_w7x_standard.nc data/local/coils.w7x_standard \
    --wp-radial 0.18 --wp-toroidal 0.15 --epseff 0.0095 --neutron-peakfactor 1.6 \
    --bt-ref 3.0 --name "W7-X standard (as built)" -o runs/tool/w7x.stella_conf.json

python scripts/compare_process_runs.py prepare stellarator_helias.IN.DAT runs/builtin/w7x50.IN.DAT --istell 5
python scripts/compare_process_runs.py prepare stellarator_helias.IN.DAT runs/tool/w7x.IN.DAT --istell 6

process -i runs/builtin/w7x50.IN.DAT
process -i runs/tool/w7x.IN.DAT

python scripts/compare_process_runs.py compare runs/builtin/w7x50.MFILE.DAT runs/tool/w7x.MFILE.DAT
```

On native Windows replace `process` with
`python scripts/run_process_windows.py`; see the notes in that script.

## Limits of this check

- The built-in values are the only reference, and several of them are rounded
  or placeholders. Agreement shows that the definitions match; it does not
  measure accuracy. Accuracy rests on the analytic tests in `tests/`.
- SIMSOPT cross-checks the field on axis and the coil lengths in CI. The
  peak field in the conductor, the energy and the forces have no cross-check
  against a second code.
- A second machine, HSX, is covered in [hsx_check.md](hsx_check.md).
- `epseff = 0.0095` is the largest core value of the effective-ripple profile
  published with the equilibrium; `neutron_peakfactor = 1.6` is the built-in
  value. Neither is computed by this tool.
