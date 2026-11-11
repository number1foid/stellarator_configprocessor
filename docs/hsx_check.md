# Second configuration: HSX

W7-X was the only real machine used during development. This check runs the
tool on a different stellarator to find anything that was tuned to W7-X by
accident.

## Why HSX

| | W7-X | HSX |
|---|---|---|
| Field periods | 5 | 4 |
| Coils | 50 (5 shapes) | 48 (6 shapes) |
| Major radius | 5.5 m | 1.2 m |
| Field on axis | 2.8 T | 1.0 T |
| Symmetry type | quasi-isodynamic | quasi-helical |

Both a vacuum equilibrium and the coil filaments are public.
`scripts/fetch_reference_data.py` downloads them and builds
`data/local/coils.hsx_qhs` (48 main coils, 150.07 kA-turns each, from
SIMSOPT's `HSX.dat`) and `data/local/wout_hsx_qhs.nc` (free-boundary vacuum
equilibrium of the quasi-helically symmetric configuration).

## Result

```
stellarator-configprocessor data/local/wout_hsx_qhs.nc data/local/coils.hsx_qhs \
    --wp-radial 0.06 --wp-toroidal 0.05 --epseff 0.005 --neutron-peakfactor 1.6 \
    --name "HSX QHS" -o hsx.stella_conf.json
```

- **The tool runs unchanged** and finishes in about one minute.
- **Consistency check passes.** Field on axis: coils 1.00344 T, equilibrium
  1.00371 T (0.03% apart). The coil field is tangent to the equilibrium's axis
  to 0.16%, and the axis links 99.99% of the coil current.
- **Symmetry is detected**: six classes of eight coils.
- **Plasma values agree with VMEC** for the same file: major radius 1.2113 m,
  minor radius 0.12149 m, volume 0.35293 m^3.
- **PROCESS 3.4.2 accepts the file** with `istell = 6` and evaluates it.

| Key | HSX | W7-X at 3 T | PROCESS built-in range |
|---|---|---|---|
| `i0` (MA) | 0.1501 | 1.745 | |
| `WP_bmax` (T) | 1.906 | 6.292 | |
| `WP_bmax / bt_ref` | 1.90 | 2.10 | 2.04 to 3.53 |
| `a1` | 0.602 | 0.641 | 0.56 to 0.98 |
| `a2` | 0.0297 | 0.0311 | 0.025 to 0.041 |
| largest error of the `a1`, `a2` fit | 0.37% | 0.06% | |
| `derivative_min_LCFS_coils_dist` | -0.54 | -0.44 | -0.61 in the example file |
| `dmin` (m) | 0.0929 | 0.2792 | |
| `min_plasma_coil_distance` (m) | 0.1429 | 0.3937 | |
| `min_bend_radius` (m) | 0.0822 | 0.3956 | |
| `inductance` (H) | 1.017e-4 | 4.182e-4 | |

The dimensionless results of two unrelated machines fall in the same narrow
range, which is also the range of PROCESS's built-in configurations.

## What the check exposed

One defect, now fixed.

**Trial winding packs that run into the neighbouring coil.** For the `a1`,
`a2` fit the tool evaluates the peak field for several winding pack sizes. It
used to reject a size only if its toroidal width exceeded the smallest
distance between filaments. HSX coils approach each other diagonally, so the
largest trial pack passed that test while its corner sat inside the
neighbour's pack. The computed peak field jumped from 1.74 T to 2.66 T and
the fit error was 23%.

Changes:

- `winding_pack.pack_clearance` measures the actual gap between the
  rectangular packs of neighbouring coils.
- A trial size is left out of the fit unless it clears its neighbours by one
  filament cell. An overlap of the reference pack itself is reported as a
  warning with the size of the overlap.
- The default trial sizes now only shrink the pack (0.5 to 1.0 times the
  reference). Larger packs are the ones that collide, and PROCESS at reactor
  scale works with a pack that is smaller relative to the coil.

W7-X was affected too, without visible symptoms: its two largest trial sizes
left less than one cell of clearance. With the new default sizes its `a1`
changed from 0.655 to 0.641 and `a2` from 0.0308 to 0.0311.

## What it did not test

- **Accuracy.** PROCESS has no built-in HSX values, and no published
  effective parameters for HSX were found. This check shows the tool is
  general and its results are plausible; it does not measure their error.
- **Unequal coil currents.** HSX's auxiliary coils carry no current in this
  configuration, like W7-X's planar coils.
- **Equilibria without stellarator symmetry**, and finite-pressure equilibria.
- `epseff = 0.005` is a placeholder of the right order for this
  configuration, not a computed value.
