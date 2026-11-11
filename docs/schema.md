# `stella_conf.json` schema (Phase 0 specification)

This document defines every value PROCESS reads from the stellarator
configuration file, how PROCESS uses it, and how this project computes it.
It is derived from the PROCESS source, not from its documentation.

- **PROCESS version read:** release `v3.4.2` (commit `c0ae5b2`, 2026-08-14).
  The schema is identical on `main` at commit `60eb219` (2026-11-10).
- **Schema version written by this tool:** `process-3.4`.
- **Method reference:** Lion et al., *Nucl. Fusion* 61 (2021) 126021 ("the paper").

## 1. How PROCESS loads the file

Source: `process/models/stellarator/preset_config.py`, `load_stellarator_config`.

- The file is read only when `istell = 6`.
- **File name:** `<output_prefix>stella_conf.json`. For an input file
  `myrun.IN.DAT` the prefix is `myrun.`, so PROCESS opens
  `myrun.stella_conf.json` next to the input file. The plain name
  `stella_conf.json` applies only when the input file is named `IN.DAT`.
- The file is one flat JSON object.
- **Keys are case-insensitive.** Each key is lower-cased and prefixed with
  `stella_config_`, then stored on the `StellaratorConfigData` dataclass
  (`process/data_structure/stellarator_configuration.py`).
- **Unknown keys are ignored without a message.**
- **Missing keys keep the dataclass default, which is `0.0`** (or `""`, `0`).
  PROCESS does not report a missing key. A missing key therefore fails
  silently or as a division by zero later. This tool always writes every key.
- This tool writes keys in the same letter case as the example file shipped
  with PROCESS (`tests/regression/input_files/stellarator_helias.stella_conf.json`).

## 2. Scaling factors PROCESS derives

These appear in the "PROCESS use" column below
(`stellarator.py`, `st_new_config`).

| Symbol | Definition |
|---|---|
| `f_R` | `rmajor / rmajor_ref` |
| `f_a` | `rminor / rminor_ref` |
| `f_B` | `b_plasma_toroidal_on_axis / bt_ref` |
| `f_N` | `n_tf_coils / (coilspermodule * symmetry)`; always 1, because PROCESS overwrites `n_tf_coils` |
| `f_I` | `coil current / i0` |
| `R_c`, `a_c` | scaled coil radii: `R_c = coil_rmajor * f_R`, `a_c = coil_rminor * f_R / f_coil_aspect` |

## 3. Keys

"Reference point" means the machine size and field described by the input
equilibrium and coils, after any rescaling requested in the run settings.

### 3.1 Identification and reference point

| Key | Type, unit | PROCESS use | Definition and calculation |
|---|---|---|---|
| `name` | string | Label only. | Run setting. |
| `symmetry` | int | `n_tf_coils = coilspermodule * symmetry`. The PROCESS docstring says "number of coils"; the code uses it as the number of field periods. | `nfp` from the equilibrium; checked against the coil set. |
| `coilspermodule` | int | See `symmetry`. | Total coil count divided by `symmetry`. Must divide exactly. |
| `rmajor_ref` | float, m | `f_R`. | Major radius of the equilibrium (`Rmajor_p`). |
| `rminor_ref` | float, m | `f_a`; coil shape factor; radial build; neoclassics. | Effective minor radius (`Aminor_p`), equal to `sqrt(V / (2 pi^2 R))`. |
| `aspect_ref` | float | Default aspect ratio when it is not an iteration variable. | `rmajor_ref / rminor_ref`. |
| `bt_ref` | float, T | `f_B`. | Mean field strength along the magnetic axis produced by the coils, paper eq. (35): `(1/l) * integral |B| ds`. On the axis `B` is tangent to the axis, so this equals the mean of the field along the axis. |

### 3.2 Plasma geometry

| Key | Type, unit | PROCESS use | Definition and calculation |
|---|---|---|---|
| `vol_plasma` | float, m^3 | `V = vol_plasma * f_R * f_a^2`. | Volume inside the last closed flux surface, paper eq. (3), from boundary Fourier coefficients. |
| `plasma_surface` | float, m^2 | `S = plasma_surface * f_R * f_a`. | Area of the last closed flux surface, paper eq. (4). |

### 3.3 Coil geometry

| Key | Type, unit | PROCESS use | Definition and calculation |
|---|---|---|---|
| `coil_rmajor` | float, m | `R_c`; peak-field fit; toroidal gap scaling. | Arc-length-weighted mean cylindrical radius of all coil filaments. Exact for circular coils centred at `R`. |
| `coil_rminor` | float, m | `a_c`; scaling of length, height, inductance, surface. | Arc-length-weighted mean distance of filament points from the circle `R = coil_rmajor`, `z = z_mean`. Exact for circular coils. The paper calls it the "average minor coil radius" without a formula, so this definition is a choice. |
| `coillength` | float, m | `len_tf_coil = coillength * (a_c / coil_rminor) / n_tf_coils`. | Sum of the lengths of all coil centre filaments. |
| `coilsurface` | float, m^2 | Cryostat area and inter-coil structure area. | Area of the toroidal surface that passes through all coil filaments (lofted between neighbouring coils at equal poloidal angle). Interpretation: PROCESS calls it "coil surface" and uses it as the area of shells covering the coils. |
| `maximal_coil_height` | float, m | `z_tf_inside_half = 0.5 * height * a_c / coil_rminor`. | `max(z) - min(z)` over all coil filaments. |
| `dmin` | float, m | `toroidalgap = dmin * (R_c - a_c) / (coil_rmajor - coil_rminor)`; constraint 82 requires it to exceed the coil toroidal width. | Minimum distance between the centre filaments of any two different coils, paper sec. 3.10. |
| `max_portsize_width` | float, m | Port widths: `0.4 * w * f_R` (vertical), `0.8 * w * f_R` (horizontal). | Largest gap between two neighbouring coil filaments: maximum over poloidal angle and coil pair of the shortest distance from one filament to its neighbour. Interpretation: the paper defines a port area, PROCESS uses a toroidal width. |
| `min_bend_radius` | float, m | `r_min = min_bend_radius * f_R / (1 - d_wp / (2 a_c))`. | `1 / kappa_max`, paper eq. (47), from a periodic spline through each filament. |
| `min_plasma_coil_distance` | float, m | Coil shape factor `(d + rminor_ref) / coil_rminor`; radial build (constraint 83). | Minimum distance between the last closed flux surface and any coil centre filament. |
| `derivative_min_LCFS_coils_dist` | float | Radial build: `available = a_c * f_shape - a + deriv * (a - f_R * rminor_ref)`. | `d(d_pc) / d(a)`, paper eq. (58): change of the minimum coil distance with the effective minor radius of the flux surface, by finite difference over flux surfaces near the boundary. Close to -1. |

### 3.4 Current, field and winding pack

The winding pack (WP) is a rectangle of radial width `d_r` and toroidal width
`d_t`, uniformly filled with current. Both dimensions are required run
settings with no default.

| Key | Type, unit | PROCESS use | Definition and calculation |
|---|---|---|---|
| `WP_area` | float, m^2 | Force-density scaling `WP_area / a_wp`. | `d_r * d_t`. |
| `WP_ratio` | float | WP shape: `a_wp = d_r^2 / WP_ratio`. | `d_r / d_t`. |
| `i0` | float, MA | `I = f_B * i0 * f_R / f_N`. | Current per coil (ampere-turns) that produces `bt_ref` at `rmajor_ref`. PROCESS assumes one current for all coils, so this is the mean absolute coil current. |
| `WP_bmax` | float, T | Force scaling `B_max / WP_bmax`. | Maximum field strength anywhere in or on the reference WP, paper eq. (38). |
| `a1`, `a2` | float | `B_max = 0.2 * I[MA] * N / (R_c - a_c) * (a1 + a2 * R_c / d_r)` (`coils.py`, `bmax_from_awp`). | Least-squares fit of `B_max` over several WP sizes at fixed `WP_ratio`. `0.2 = mu0 / (2 pi) * 1e6`, so `a1 = 1` is an ideal toroid. The code fits against `R_c / d_r`; the paper's eq. (39) uses `R / sqrt(A_wp)` and names the parameters `a0`, `a1`. This tool follows the code. |
| `inductance` | float, H | `W = 0.5 * L * I^2` with `I` the current per coil; scaled as `L * (a_c / coil_rminor)^2 / f_R`. | Effective inductance `L = 2 W / i0^2`, where `W` is the magnetic energy of the whole coil set. For equal currents this is the sum of all entries of the inductance matrix. Self-inductances use the WP cross-section. |

### 3.5 Forces

All values are at the reference point. `F'` is the force per unit length on a
coil, `F' = I t x <B>`, with `<B>` the field averaged over the WP
cross-section and `t` the filament tangent.

| Key | Type, unit | PROCESS use | Definition and calculation |
|---|---|---|---|
| `max_force_density_MNm` | float, MN/m | Scaled with `f_I * B_max / WP_bmax`. | Maximum of `|F'|` over all coils and positions. |
| `max_force_density` | float, MN/m^3 | Scaled with `f_I * (B_max / WP_bmax) * (WP_area / a_wp)`; gives the WP stress for constraint 32. | `max_force_density_MNm / WP_area`. The example file satisfies this relation exactly. |
| `max_radial_force_density` | float, MN/m^3 | Output only. | Maximum of `|F' . n_r| / WP_area`, with `n_r` the in-plane normal of the coil. PROCESS labels the unit MN/m in one place; the scaling and the example values show it is MN/m^3. |
| `max_lateral_force_density` | float, MN/m^3 | Output only. | Maximum of `|F' . n_t| / WP_area`, with `n_t` the toroidal (out-of-plane) normal. |
| `centering_force_max_MN` | float, MN | Output only. | Largest net force on one coil, projected on the direction from the coil centroid towards the torus axis. |
| `centering_force_min_MN` | float, MN | Output only. | Smallest such force. Negative means outwards. |
| `centering_force_avg_MN` | float, MN | Output only. | Mean over all coils. |

WP orientation follows the paper: `n_t` is the cylindrical toroidal direction
made perpendicular to the tangent, and `n_r = t x n_t`.

### 3.6 Values this tool does not compute

| Key | Type, unit | PROCESS use | Handling |
|---|---|---|---|
| `epseff` | float | Neoclassical flux check (`neoclassics.py`). | Needs a neoclassical code. Required run setting, passed through. The built-in configurations all use 0.015. |
| `neutron_peakfactor` | float | First-wall peaking factor `wallpf`. | Needs a wall model, paper eq. (32) to (34). Required run setting, passed through. The built-in configurations all use 1.6. |

### 3.7 Keys in the example file that PROCESS ignores

`number_nu_star`, `D11_star_mono_input`, `nu_star_mono_input`. They have no
field on the dataclass, so the loader drops them. This is the neoclassical
transport path the documentation calls "not yet accessible". This tool does
not write them.

## 4. Reference data

### 4.1 Built-in configurations

From `preset_config.py`. Note the mapping: **`istell = 4` is the 30-coil set
and `istell = 5` is the 50-coil set.** The PROCESS documentation states the
opposite.

| Key | Helias 5 (`1`) | Helias 4 (`2`) | Helias 3 (`3`) | W7X-30 (`4`) | W7X-50 (`5`) |
|---|---|---|---|---|---|
| `symmetry` | 5 | 4 | 3 | 5 | 5 |
| `coilspermodule` | 10 | 10 | 10 | 6 | 10 |
| `rmajor_ref` | 22.2 | 17.6 | 13.86 | 5.50 | 5.5 |
| `rminor_ref` | 1.80 | 2.0 | 2.18 | 0.49 | 0.49 |
| `aspect_ref` | 12.33 | 8.8 | 6.36 | 11.2 | 11.2 |
| `coil_rmajor` | 22.44 | 18.39 | 14.53 | 5.62 | 5.62 |
| `coil_rminor` | 4.76 | 4.94 | 6.12 | 1.36 | 1.18 |
| `bt_ref` | 5.6 | 5.6 | 5.6 | 3.0 | 3.0 |
| `WP_area` | 0.48 | 0.48 | 0.48 | 0.027 | 0.027 |
| `WP_ratio` | 1.2 | 1.3 | 1.3 | 1.2 | 1.2 |
| `WP_bmax` | 11.44 | 11.51 | 12.346 | 10.6 | 6.3 |
| `i0` | 13.06 | 13.146 | 14.23 | 2.9 | 1.745 |
| `a1` | 0.688 | 0.676 | 0.56 | 0.98 | 0.66 |
| `a2` | 0.025 | 0.029 | 0.030 | 0.041 | 0.025 |
| `vol_plasma` | 1422.63 | 1380.0 | 1300.8 | 26.4 | 26.4 |
| `plasma_surface` | 1960.0 | 1900.0 | 1600.0 | 128.3 | 128.3 |
| `dmin` | 0.84 | 1.08 | 1.145 | 0.21 | 0.28 |
| `max_portsize_width` | 2.12 | 3.24 | 3.24 | 0.5 | 0.3 |
| `maximal_coil_height` | 12.7 | 13.34 | 17.74 | 3.6 | 3.1 |
| `coilsurface` | 4817.7 | 4100.0 | 4240.0 | 370.0 | 299.85 |
| `coillength` | 1680.0 | 1435.07 | 1287.3 | 303.4 | 420.67 |
| `inductance` | 1655.76e-6 | 1290.4e-6 | 1250.7e-6 | 252.7e-6 | 412.4e-6 |
| `max_force_density` | 120.0 | 120.0 | 120.0 | 350.0 | 250.0 |
| `max_force_density_MNm` | 98.0 | 98.0 | 98.0 | 98.0 | 98.0 |
| `max_lateral_force_density` | 92.4 | 87.9 | 96.6 | 271.1 | 116.4 |
| `max_radial_force_density` | 113.5 | 109.9 | 130.5 | 305.2 | 148.0 |
| `centering_force_max_MN` | 189.5 | 226.0 | 428.1 | 7.95 | 2.99 |
| `centering_force_min_MN` | -55.7 | -35.3 | -70.3 | -2.15 | -1.29 |
| `centering_force_avg_MN` | 93.0 | 125.8 | 240.9 | 3.46 | 1.61 |
| `min_plasma_coil_distance` | 1.9 | 1.7 | 1.78 | 0.45 | 0.39 |
| `derivative_min_LCFS_coils_dist` | -1.0 | -1.0 | -1.0 | -1.0 | -1.0 |
| `min_bend_radius` | 1.0 | 0.86 | 1.145 | 0.186 | 0.39 |
| `neutron_peakfactor` | 1.6 | 1.6 | 1.6 | 1.6 | 1.6 |
| `epseff` | 0.015 | 0.015 | 0.015 | 0.015 | 0.015 |

Observations on the built-in values:

- Several entries are rounded or placeholder values (`max_force_density_MNm`
  is 98.0 in all five; `derivative_min_LCFS_coils_dist` is "approximated").
- The built-in `max_force_density` does not equal
  `max_force_density_MNm / WP_area` (98 / 0.48 = 204, listed 120). The example
  file in section 4.2 does satisfy it. The built-in values are therefore not a
  self-consistent output of one pre-processing run.
- The built-in `a1`, `a2` differ in character from the example file
  (0.688, 0.025 against 0.025, 0.056 for Helias 5). Both reproduce
  `WP_bmax` to within 10% through the PROCESS formula.

### 4.2 Example file shipped with PROCESS

`tests/regression/input_files/stellarator_helias.stella_conf.json`, with
`"name": "presprocess"`. It is the only public output of the original tool
found. It is a Helias 5 configuration.

| Key | Value | Key | Value |
|---|---|---|---|
| `symmetry` | 5 | `coilspermodule` | 10 |
| `rmajor_ref` | 22.19309491 | `rminor_ref` | 1.80206932 |
| `aspect_ref` | 12.31533919 | `bt_ref` | 5.6 |
| `coil_rmajor` | 22.22954927 | `coil_rminor` | 4.69015768 |
| `WP_area` | 0.53333333 | `WP_ratio` | 1.2 |
| `WP_bmax` | 11.49413032 | `i0` | 12.98310033 |
| `a1` | 0.0254274 | `a2` | 0.05608523 |
| `vol_plasma` | 1422.62552585 | `plasma_surface` | 1960.01361974 |
| `dmin` | 0.96033081 | `max_portsize_width` | 2.78496289 |
| `maximal_coil_height` | 12.23652989 | `min_bend_radius` | 1.54759778 |
| `coilsurface` | 4740.28625669 | `coillength` | 1698.73434821 |
| `inductance` | 0.00137435 | `min_plasma_coil_distance` | 1.91756031 |
| `derivative_min_LCFS_coils_dist` | -0.61346266 | `max_force_density` | 120.68504161 |
| `max_force_density_MNm` | 64.36535552 | `max_lateral_force_density` | 92.35315199 |
| `max_radial_force_density` | 113.63068984 | `centering_force_max_MN` | 189.52429886 |
| `centering_force_min_MN` | -55.44479209 | `centering_force_avg_MN` | 93.02700909 |
| `epseff` | 0.01464553 | `neutron_peakfactor` | 1.55835632 |

### 4.3 What the W7-X built-in values represent

Section 4.3 of the paper states that the 30- and 50-coil sets were generated
with the coil optimisation code FOCUS to target a W7-X-like equilibrium. They
are not the as-built W7-X coils. The real W7-X has 50 non-planar coils, so its
closest built-in reference is **`istell = 5` (W7X-50)**, not `istell = 4`.
The Phase 5 comparison uses `istell = 5`, and an approximate match is the
most that can be expected.

## 5. Accuracy targets

Tolerances for tests against analytic and cross-code references.

| Quantity group | Target |
|---|---|
| Plasma volume, surface, radii | 0.1% |
| Coil length, height, radii, distances | 0.5% |
| Field on axis, `i0` | 0.5% |
| Inductance and energy | 2% |
| Peak field in the WP, `a1`, `a2` fit | 3% |
| Forces | 5% |

## 6. Open points

1. `coil_rminor`, `coilsurface` and `max_portsize_width` have no published
   formula. The definitions above are choices and are named as such in the
   run report.
2. `epseff` and `neutron_peakfactor` are passed through, not computed.
3. A coil set with unequal coil currents (for example W7-X with its planar
   coils energised) cannot be represented exactly, because PROCESS scales one
   current. The report lists the individual currents.
