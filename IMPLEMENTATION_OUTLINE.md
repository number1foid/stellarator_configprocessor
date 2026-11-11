# Stellarator pre-processor for PROCESS: implementation outline

## 1. Purpose

PROCESS, the UKAEA fusion power plant systems code, can model a stellarator only if it is given a configuration file (`stella_conf.json`) of pre-calculated "effective parameters" for that specific design. The tool meant to produce this file, pre-sPROCESS, is not publicly available as far as could be found.

This project, `stellarator_configprocessor`, builds an open replacement: a command-line program that reads a stellarator's plasma equilibrium and coil set and writes the configuration file PROCESS reads with `istell=6`.

- Git repository: https://github.com/number1foid/stellarator_configprocessor
- Local repository: `D:\NuclearEng\fusion\stellerator_configprocessor`

## 2. Inputs and output

| | What | Format |
|---|---|---|
| Input | Plasma equilibrium | VMEC `wout` file (NetCDF), for example from VMEC++ |
| Input | Coil set | Filament curves and currents (MAKEGRID `coils.*` file, or SIMSOPT curve objects) |
| Input | Run settings | Reference major radius, reference field, conductor cross-section dimensions (required, no default), resolution options |
| Output | Configuration file | `stella_conf.json`, in the schema PROCESS expects |
| Output | Report | Plain-text or JSON log of every computed value and how it was obtained |

## 3. What is known and what is not

**Confirmed from public sources**

- The method is defined in Lion et al., *Nuclear Fusion* 61 (2021) 126021, open access (CC BY 4.0). Its inputs are a reference equilibrium and the associated coil filaments.
- PROCESS is MIT-licensed, so the code that reads `stella_conf.json` can be inspected.
- The PROCESS documentation says the file carries inductances, peak field values and coil force quantities from reference calculations, and that plasma geometry comes from Fourier coefficients of the last closed flux surface.
- PROCESS scales the fixed configuration only in overall size. Coil number and coil shape are fixed by the file.
- SIMSOPT ships the real W7-X coil set (`get_data("w7x")`) and W7-X equilibrium inputs. VMEC++ ships a W7-X equilibrium input.

**Not confirmed**

- The exact list of keys, units and conventions in `stella_conf.json`. The documentation does not list them. They must be read from the PROCESS source (step 1 below).
- The equilibrium and coil files that produced the built-in PROCESS values (`istell=1` to `5`). No HELIAS files were found in the SIMSOPT or VMEC++ repositories.
- Whether the neoclassical transport coefficients can be used at all. The documentation says this path is "not yet accessible".

## 4. Architecture

```
equilibrium file ──► equilibrium reader ──┐
                                          ├──► parameter calculators ──► schema writer ──► stella_conf.json
coil file ─────────► coil reader ─────────┘            │
                                                       └──► report
```

| Module | Responsibility |
|---|---|
| `io.equilibrium` | Read the `wout` file; expose boundary Fourier coefficients, volume, rotational transform, field on axis |
| `io.coils` | Read coil filaments and currents; apply field-period and stellarator symmetry to build the full set |
| `geometry.plasma` | Plasma volume, surface area, effective minor radius, aspect ratio, boundary coefficients |
| `geometry.coils` | Coil lengths, minimum coil-to-coil distance, minimum coil-to-plasma distance, enclosed areas |
| `field` | Biot–Savart evaluation; field on axis; peak field on the conductor |
| `energy` | Stored magnetic energy and inductances |
| `forces` | Force densities on the coils |
| `schema` | Map computed values to the exact keys and units PROCESS expects; validate; write JSON |
| `cli` | Argument parsing, run orchestration, report output |

Each calculator is a pure function of the reader outputs, so it can be tested in isolation.

## 5. Build phases

### Phase 0: Specification (no code)

1. Read the PROCESS stellarator source and list every key read from `stella_conf.json`, with type, unit, and where it is used.
2. Match each key to its definition in the Lion et al. paper.
3. Extract the built-in values for `istell=1` to `5` from PROCESS as reference data.
4. Write the result as a schema document in this repository.

*Done when:* every key has a written definition, a unit, and a planned calculation method.

### Phase 1: Readers

1. Equilibrium reader for `wout` files.
2. Coil reader for the chosen coil format(s), including symmetry expansion.
3. Consistency check that the coils and equilibrium describe the same machine (field periods, axis position, field on axis).

*Done when:* the W7-X equilibrium and coils from SIMSOPT load and pass the consistency check.

### Phase 2: Plasma geometry

1. Volume, surface area, effective minor radius, aspect ratio.
2. Boundary Fourier coefficients in the convention PROCESS uses.

*Done when:* results for a circular torus match the analytic formulas, and W7-X values match those reported by VMEC++.

### Phase 3: Coil geometry and fields

1. Coil lengths and minimum distances.
2. Biot–Savart field at arbitrary points; field on axis.
3. Peak field on the conductor. This needs a finite conductor cross-section, because the field of an ideal filament is infinite at the filament itself.

*Done when:* a single circular coil reproduces the analytic on-axis field, and W7-X field on axis agrees with the equilibrium.

### Phase 4: Energy, inductance, forces

1. Stored magnetic energy and self and mutual inductances.
2. Force densities along each coil.

*Done when:* inductance of a circular loop matches the analytic formula, and energy from inductances agrees with energy from a volume integral of the field.

### Phase 5: Schema writer and PROCESS integration

1. Write `stella_conf.json`.
2. Run PROCESS with `istell=6` on the generated file.
3. Compare against the built-in `istell=5` (W7X-50) run. Phase 0 found that `istell=4` is a 30-coil set and that both built-in W7-X sets are FOCUS-generated, not the as-built coils; see `docs/schema.md` section 4.3.

*Done when:* PROCESS accepts the file and the W7-X run lands close to the built-in result, with every difference explained.

### Phase 6: Command-line tool, documentation, release

1. Command-line interface and report output.
2. User guide with one worked example from files to PROCESS run.
3. Automated tests in continuous integration.

### Phase 7 (optional extensions)

- Generate a family of configuration files across aspect ratio or coil number.
- Neoclassical transport coefficients, if PROCESS makes that path usable.
- Additional equilibrium or coil formats.

## 6. Verification strategy

| Level | Method | Needs non-public data? |
|---|---|---|
| Unit | Analytic cases: circular torus, circular coil, straight wire | No |
| Cross-code | Compare plasma quantities with VMEC++ output and fields with SIMSOPT | No |
| Reference | Reproduce the PROCESS built-in W7-X values from the public W7-X files | No, but expect an approximate match |
| Exact | Reproduce built-in HELIAS values | Yes: the original inputs would have to be requested from the authors |

## 7. Risks

| Risk | Effect | Mitigation |
|---|---|---|
| PROCESS version 3 is declared unstable; the schema may change | Output stops loading | Pin a PROCESS release; tag each output with a schema version |
| Conventions differ between codes (sign of the toroidal angle, units, normalisation) | Silently wrong values | Consistency checks in Phase 1; cross-code tests |
| Peak field depends on the assumed conductor cross-section | Result is sensitive to an arbitrary input | Make the cross-section a required, reported input with no default |
| Public W7-X files differ from those behind the built-in values | No exact reference match | Rely on analytic and cross-code checks for correctness |
| The original authors release pre-sPROCESS | Duplicate effort | None: risk accepted. Contact is deferred until there is a working result (decision 12) |

## 8. Sources

- Lion et al. 2021, open-access PDF: https://scientific-publications.ukaea.uk/wp-content/uploads/LION_2021_NUCL-_FUSION_61_126021.PDF
- PROCESS stellarator documentation: https://ukaea.github.io/PROCESS/fusion-devices/stellarator/
- PROCESS repository: https://github.com/ukaea/process
- SIMSOPT repository: https://github.com/hiddenSymmetries/simsopt
- VMEC++ repository: https://github.com/proximafusion/vmecpp

---

## 9. Creative or arbitrary decisions

These choices have no single correct answer and shape the long-term direction of the software. Each lists the options, a suggestion, and a space for your decision. Decisions 1 to 6 should be made before Phase 1.

### 1. Name

The project and package name. It should not be "pre-sPROCESS", which belongs to the original tool.

**Your decision: "stellarator_configprocessor"**

### 2. Relationship to PROCESS

- **A. Standalone repository** that outputs a file PROCESS reads. Full control, independent release schedule.
- **B. Contribution into the PROCESS repository.** Maximum reach, but bound to their review process, code style and roadmap.
- **C. Standalone first, offer upstream later.**

*Suggestion:* C.

**Your decision: "Option A. Option C is similar to what I have in mind, but I can't say that I will 100% offer upstream"**

The repository is https://github.com/number1foid/stellarator_configprocessor, with the local copy at `D:\NuclearEng\fusion\stellerator_configprocessor`.

### 3. License

- **MIT:** matches PROCESS and VMEC++; anyone, including companies, can use it freely.
- **LGPL or GPL:** changes must be shared back; may deter commercial users and complicates upstreaming into PROCESS.
- **Apache 2.0:** like MIT with an explicit patent grant.

*Suggestion:* MIT.

**Your decision: "MIT"**

### 4. Language

- **A. Pure Python with NumPy.** Simplest; matches PROCESS; fast enough if the field calculations are vectorised.
- **B. Python with a compiled core (C++ or Rust).** Faster for fine coil discretisation; harder to build and distribute.
- **C. Python with JAX.** Makes every output differentiable with respect to the inputs, which allows gradient-based optimisation later; steeper learning curve.

*Suggestion:* A, structured so the field kernel can be swapped later.

**Your decision: "A"**

### 5. Dependence on SIMSOPT and VMEC++

- **A. Hard dependency on SIMSOPT.** Reuses its coil classes, Biot–Savart and file readers; least code to write; heavier installation and tied to its release cycle.
- **B. Own readers and own Biot–Savart, with SIMSOPT optional.** More code and more to verify; lightweight and self-contained.
- **C. Own code, with SIMSOPT used only in tests as a cross-check.**

*Suggestion:* C. Independent code is what makes the cross-check meaningful.

**Your decision: "C"**

### 6. Scope of the first release

- **A. One configuration in, one file out.** Fills the documented gap.
- **B. Also generate families of configurations** across aspect ratio or coil number.
- **C. Also include neoclassical transport coefficients.**

*Suggestion:* A, with B as the first extension.

**Your decision: "A, we indeed do B later"**

### 7. Supported input formats at launch

- Equilibrium: VMEC `wout` only, or also DESC output?
- Coils: MAKEGRID `coils.*` only, or also SIMSOPT JSON and FOCUS files?

*Suggestion:* `wout` and MAKEGRID first; add others on request.

**Your decision: "wout and MAKEGRID for now"**

### 8. Conductor model for peak field

- **A. Rectangular cross-section of user-given dimensions,** filled with many filaments.
- **B. Single filament, field sampled at a fixed offset distance.** Fast but crude.
- **C. Both, selectable.**

Also decide the default dimensions, since the result depends on them.

*Suggestion:* A, with dimensions as a required input and no silent default.

**Your decision: "A. The dimensions are a required input, with no default."**

### 9. Handling PROCESS schema changes

- **A. Support only the latest PROCESS release.**
- **B. Support several PROCESS versions** through a `--process-version` option.
- **C. Define an own stable intermediate format** and convert to each PROCESS schema.

*Suggestion:* A at first; move to C if the schema changes more than once.

**Your decision: "A, C may be added later"**

### 10. Accuracy target

What counts as "matching" a reference value: for example 1%, 0.1%, or a per-quantity tolerance. This sets the default resolution and the test thresholds.

*Suggestion:* per-quantity tolerances, tighter for geometry than for peak field and forces.

**Your decision: "Per-quantity tolerances, tighter for geometry than for peak field and forces."**

### 11. Behaviour on suspicious input

When the coils and equilibrium appear inconsistent (for example field on axis differs by several percent):

- **A. Stop with an error.**
- **B. Warn and continue.**
- **C. Stop by default, with an override flag.**

*Suggestion:* C.

**Your decision: "C"**

### 12. Contact with the original authors and PROCESS maintainers

- **A. Contact them before starting,** to ask whether pre-sPROCESS will be released and whether the reference inputs can be shared.
- **B. Build first, contact once there is a working result.**
- **C. No contact.**

*Suggestion:* A. It could save the whole effort or supply the missing HELIAS inputs.

**Your decision: "B. The risk that pre-sPROCESS is released in the meantime is accepted."**

### 13. Development in public

- **A. Public repository from the first commit.**
- **B. Private until the W7-X check passes.**

**Your decision: "Private until the W7-X check passes"**

### 14. Long-term direction

What should this become if it works?

- **A. A small, finished utility** that does one job.
- **B. A general bridge between stellarator design codes and systems codes,** supporting systems codes other than PROCESS.
- **C. Part of a design optimisation loop,** where a stellarator optimiser calls it to score candidate designs on power plant cost.

This affects decisions 4, 5 and 9.

**Your decision: "Will be decided when the current scope is finished."**