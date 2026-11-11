![A stellarator during assembly, surrounded by scaffolding](docs/images/stellarator_assembly.webp)

# stellarator_configprocessor

[Nederlands](#nederlands) | [English](#english)

---

<a id="nederlands"></a>

# Nederlands

Een open preprocessor voor het stellaratormodel van
[PROCESS](https://github.com/ukaea/process), de systeemcode voor
fusiecentrales.

PROCESS kan elke modulaire stellarator modelleren als het een
configuratiebestand (`stella_conf.json`) krijgt met vooraf berekende
"effectieve parameters" voor dat ontwerp. Het programma dat de documentatie
van PROCESS noemt om dit bestand te maken, pre-sPROCESS, is niet openbaar
beschikbaar. Dit project berekent hetzelfde bestand uit twee gangbare
invoerbestanden:

- een plasma-evenwicht: VMEC-`wout`-bestand (NetCDF, klassiek formaat);
- een spoelenset: MAKEGRID-bestand `coils.*` met filamentkrommen en stromen.

De methode volgt Lion et al., *Nuclear Fusion* 61 (2021) 126021.

**Deze tak (`rust`) bevat de Rust-versie**, geschreven volgens de tien regels
van "The Power of Ten" van NASA/JPL voor veiligheidskritische code. De
oorspronkelijke Python-versie staat op de tak `main`. Beide geven dezelfde
resultaten.

## Status

Ontwikkelversie, nog niet uitgebracht. Geschreven voor PROCESS 3.4.2.

- De Rust-versie geeft voor W7-X, HSX en een herschaalde W7-X dezelfde
  waarden als de Python-versie, tot op beter dan 1e-10, en is ongeveer zes
  keer zo snel (9 s tegenover 51 s voor W7-X).
- De code voldoet aan de tien regels; een geautomatiseerde controle bewaakt
  dat. Zie [docs/power_of_ten.md](docs/power_of_ten.md) voor wat dat precies
  inhoudt, hoe het wordt gecontroleerd en waar een regel voor Rust is
  uitgelegd.

- Alle 35 waarden die PROCESS leest, worden geschreven. Twee daarvan
  (`epseff`, `neutron_peakfactor`) worden doorgegeven vanaf de opdrachtregel
  en niet berekend.
- PROCESS 3.4.2 accepteert het gegenereerde bestand met `istell = 6`.
- Voor de W7-X-spoelen zoals gebouwd komen de berekende waarden binnen enkele
  procenten overeen met de ingebouwde W7-X-set van PROCESS. De verschillen
  worden verklaard in [docs/w7x_comparison.md](docs/w7x_comparison.md).

## Installatie

Installeer de Rust-toolchain via <https://rustup.rs> (versie 1.82 of nieuwer).
Er zijn geen andere afhankelijkheden.

```
cargo build --release
```

Het programma staat daarna in `target/release/stellarator-configprocessor`
(`.exe` op Windows). Met `cargo install --path .` komt het in de `bin`-map
van Cargo, zodat het overal te starten is.

## Gebruik

```
stellarator-configprocessor wout_mydesign.nc coils.mydesign \
    --wp-radial 0.6 --wp-toroidal 0.5 \
    --epseff 0.015 --neutron-peakfactor 1.6 \
    --process-input myrun.IN.DAT --report myrun_report.txt
```

Dit schrijft `myrun.stella_conf.json` naast `myrun.IN.DAT`, de plek waar
PROCESS het zoekt, en een rapport met elke waarde, de eenheid en de manier
waarop ze is verkregen.

De vier getoonde instellingen hebben geen standaardwaarde: de resultaten
hangen af van de afmetingen van het wikkelpakket, en de laatste twee kan dit
programma niet berekenen.

Zie [docs/user_guide.md](docs/user_guide.md) voor alle opties en een
uitgewerkt voorbeeld van bestanden tot een PROCESS-run.

## Documentatie

De documenten zelf zijn in het Engels.

| Document | Inhoud |
|---|---|
| [docs/user_guide.md](docs/user_guide.md) | Opties, invoer, uitgewerkt W7-X-voorbeeld, gebruik als bibliotheek, gebruiksbeperkingen |
| [docs/power_of_ten.md](docs/power_of_ten.md) | Hoe de code aan elk van de tien regels van NASA/JPL voldoet en hoe dat wordt gecontroleerd |
| [docs/schema.md](docs/schema.md) | Elke sleutel van `stella_conf.json`: eenheid, hoe PROCESS haar gebruikt, hoe ze wordt berekend; referentiewaarden |
| [docs/w7x_comparison.md](docs/w7x_comparison.md) | PROCESS-run met het gegenereerde W7-X-bestand tegenover de ingebouwde W7-X-set |
| [docs/hsx_check.md](docs/hsx_check.md) | Tweede machine (HSX): controle op algemene toepasbaarheid en wat die aan het licht bracht |

## Verificatie

```
python scripts/fetch_reference_data.py    # optioneel: openbare W7-X- en HSX-bestanden (Python met NumPy)
cargo fmt --all -- --check
cargo clippy --workspace --all-targets
cargo test --workspace --release
```

`cargo test` voert 42 tests uit: de regelcontrole, de test dat tijdens de
berekening geen geheugen wordt aangevraagd, en de functionele tests.

| Grootheid | Gecontroleerd tegen | Overeenkomst |
|---|---|---|
| Alle uitvoerwaarden voor W7-X, HSX en een herschaalde W7-X | de Python-versie | 1e-10 (test eist 1e-8) |
| Plasmavolume, -oppervlak, stralen | cirkelvormige torus (analytisch); de waarden van VMEC zelf voor W7-X | 1e-12; 1e-9 |
| Veld van een filament | cirkelvormige lus en rechte draad (analytisch) | 1e-4; 1e-6 |
| Veld op de as van de W7-X-spoelen | het vacuümevenwicht van W7-X | 5e-6 |
| Zelfinductie | formule voor een dunne cirkelvormige ring | 0,2% |
| Kern van de rechthoekige staaf | limiet ver weg en symmetrie in het midden | 1e-5 |
| Fit van a1 en a2 | teruggevonden uit zelf opgegeven waarden | 1e-10 |
| Torus met 20 cirkelvormige spoelen | analytische geometrie, ideale toroïdale spoel | 1e-10 tot 1% |
| Herschaling | schaalwetten voor lengte, veld, stroom en inductie | 1e-8 |
| Volledig W7-X-resultaat | ingebouwde W7X-50-waarden van PROCESS | zie docs/w7x_comparison.md |
| Volledig HSX-resultaat | consistentie met het evenwicht; plausibiliteit | zie docs/hsx_check.md |

Tests die de W7-X- of HSX-bestanden nodig hebben, melden dat ze zijn
overgeslagen als die bestanden ontbreken.

De Python-versie op `main` heeft daarnaast analytische tests voor de
wederzijdse inductie, de hoepelkracht, virtuele arbeid en het veld binnen een
staaf, en een kruiscontrole met SIMSOPT. Die zijn niet naar Rust overgezet;
voor de Rust-versie gelden ze via de overeenkomst met de Python-resultaten.

## Bekende problemen

### Nog niet geverifieerd

- **Geen onafhankelijke controle van piekveld, energie en krachten op een
  echte machine.** Deze berusten op analytische testgevallen en op
  overeenkomst met de ingebouwde W7-X-waarden van PROCESS, waarvan sommige
  zijn afgerond of plaatshouders zijn. De kruiscontrole met SIMSOPT dekt
  alleen het filamentveld en de spoellengtes.
- **Twee onverklaarde verschillen met de ingebouwde W7-X-waarden.** `a2` ligt
  25% boven de ingebouwde waarde en `max_lateral_force_density` 10% erboven.
  De fit reproduceert de velden die dit programma berekent tot op 0,06%;
  welke kant gelijk heeft, staat niet vast.
- **Drie definities zijn beredeneerde aannames.** Voor `coil_rminor`,
  `coilsurface` en `max_portsize_width` is geen formule gepubliceerd. De
  gebruikte definities zijn gekozen omdat ze passen bij de beschikbare
  referentiewaarden; alleen de auteurs van het oorspronkelijke programma
  kunnen ze bevestigen. Zie [docs/schema.md](docs/schema.md).
- **Geen volledig geoptimaliseerd PROCESS-ontwerp.** PROCESS leest en
  evalueert de gegenereerde bestanden. Een optimalisatie op reactorschaal
  vond geen haalbare centrale voor W7-X, noch met de gegenereerde, noch met
  de ingebouwde configuratie, en voor een andere machine is geen
  optimalisatie uitgevoerd.
- **Paden zonder test op echte gegevens.** Spoelensets met ongelijke stromen,
  evenwichten zonder stellaratorsymmetrie en evenwichten met eindige druk
  worden door de code afgehandeld, maar zijn alleen op synthetische gevallen
  of helemaal niet uitgevoerd.
- **De nauwkeurigheid voor HSX is onbekend.** De tweede machine laat zien dat
  het programma algemeen toepasbaar is en plausibele resultaten geeft; er
  zijn geen referentiewaarden voor.

### Niet gedaan

- **`epseff` en `neutron_peakfactor` worden niet berekend.** Het zijn
  verplichte invoerwaarden die worden doorgegeven aan het bestand. Zolang ze
  niet worden berekend, vervangt het programma pre-sPROCESS niet volledig.
- **PROCESS maakt geen deel uit van de geautomatiseerde tests.** Het schema
  is alleen via de documenten aan PROCESS 3.4.2 gekoppeld, zodat een
  schemawijziging in PROCESS niet automatisch wordt opgemerkt.

### Alleen in de Rust-versie

- **Alleen het klassieke NetCDF-formaat.** `wout`-bestanden in NetCDF-4
  (HDF5) worden geweigerd; zet ze om met `nccopy -k classic`. VMEC2000
  schrijft standaard het klassieke formaat.
- **Vaste bovengrenzen.** Bijvoorbeeld 512 spoelen, 1024 punten per spoel en
  12 filamenten per zijde van het wikkelpakket. Invoer daarboven wordt
  geweigerd. Zie `src/limits.rs`.
- **Eén processorkern.** Threads vragen geheugen aan en zijn daarom niet
  gebruikt.
- **Leeg uitvoerbestand na een mislukte run.** De uitvoerbestanden worden
  vóór de berekening geopend. Een bestaand bestand blijft bij een fout
  ongewijzigd; een nieuw bestand blijft leeg achter.
- **Niet alle analytische tests zijn overgezet**, zie Verificatie.
- **De dagelijkse statische analyse draait nog niet.** GitHub start geplande
  runs alleen vanaf de standaardtak; op deze tak draait de analyse bij elke
  push.
- **De regels gelden voor het programma zelf**, niet voor de tests, de
  Python-hulpscripts en de standaardbibliotheek van Rust.

### Bewuste keuzes

- **Alleen modulaire spoelen.** Elke spoel moet poloïdaal om het plasma
  sluiten.
- **Eén spoelstroom.** PROCESS schaalt één stroom, dus spoelensets met
  ongelijke stromen worden teruggebracht tot hun gemiddelde absolute stroom,
  met een waarschuwing. Spoelen zonder stroom worden genegeerd.
- **Invoerformaten.** Alleen VMEC-`wout` en MAKEGRID.
- **Geen externe bibliotheken.** Geen enkele bestaande Rust-bibliotheek is
  volgens deze regels geschreven, dus de lezers en schrijvers zijn eigen code.

### In PROCESS zelf

Gevonden bij het lezen en uitvoeren van PROCESS 3.4.2; details in
[docs/schema.md](docs/schema.md), sectie 6.

- De documentatie verwisselt `istell = 4` en `istell = 5`: in de code is 4
  de W7-X-set met 30 spoelen en 5 de set met 50 spoelen.
- Een sleutel die in het configuratiebestand ontbreekt, wordt stilzwijgend
  als nul gelezen.
- De centrerende krachten worden geschaald met de verhouding van de
  spoellengtes omgekeerd ten opzichte van vgl. (43) van Lion et al. 2021.
- Op Windows zelf mislukt een run zonder de noodoplossingen in
  `scripts/run_process_windows.py`.

## Licentie

MIT. Zie [LICENSE](LICENSE).

---

<a id="english"></a>

# English

An open pre-processor for the stellarator model of the
[PROCESS](https://github.com/ukaea/process) fusion power plant systems code.

PROCESS can model any modular stellarator if it is given a configuration file
(`stella_conf.json`) of pre-calculated "effective parameters" for that design.
The tool that PROCESS's documentation names for producing this file,
pre-sPROCESS, is not publicly available. This project computes the same file
from two standard inputs:

- a plasma equilibrium: VMEC `wout` file (NetCDF, classic format);
- a coil set: MAKEGRID `coils.*` file with filament curves and currents.

The method follows Lion et al., *Nuclear Fusion* 61 (2021) 126021.

**This branch (`rust`) holds the Rust version**, written to the ten rules of
NASA/JPL's "The Power of Ten" for safety-critical code. The original Python
version is on the `main` branch. Both give the same results.

## Status

Development version, not yet released. Written against PROCESS 3.4.2.

- For W7-X, HSX and a rescaled W7-X the Rust version gives the same values
  as the Python version, to better than 1e-10, and is about six times faster
  (9 s against 51 s for W7-X).
- The code meets the ten rules, and an automated check keeps it that way.
  See [docs/power_of_ten.md](docs/power_of_ten.md) for what that means
  exactly, how it is checked, and where a rule had to be interpreted for
  Rust.

- All 35 values PROCESS reads are written. Two of them (`epseff`,
  `neutron_peakfactor`) are passed through from the command line, not
  computed.
- PROCESS 3.4.2 accepts the generated file with `istell = 6`.
- For the as-built W7-X coils, the computed values agree with PROCESS's
  built-in W7-X set to within a few percent, and the differences are
  explained in [docs/w7x_comparison.md](docs/w7x_comparison.md).

## Install

Install the Rust toolchain from <https://rustup.rs> (version 1.82 or newer).
There are no other dependencies.

```
cargo build --release
```

The program is then `target/release/stellarator-configprocessor` (`.exe` on
Windows). `cargo install --path .` copies it to Cargo's `bin` folder, so that
it can be run from anywhere.

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
| [docs/user_guide.md](docs/user_guide.md) | Options, inputs, worked W7-X example, use as a library, usage limits |
| [docs/power_of_ten.md](docs/power_of_ten.md) | How the code meets each of the ten NASA/JPL rules and how that is checked |
| [docs/schema.md](docs/schema.md) | Every key of `stella_conf.json`: unit, how PROCESS uses it, how it is computed; reference values |
| [docs/w7x_comparison.md](docs/w7x_comparison.md) | PROCESS run with the generated W7-X file against the built-in W7-X set |
| [docs/hsx_check.md](docs/hsx_check.md) | Second machine (HSX): generality check and what it exposed |

## Verification

```
python scripts/fetch_reference_data.py    # optional: public W7-X and HSX files (Python with NumPy)
cargo fmt --all -- --check
cargo clippy --workspace --all-targets
cargo test --workspace --release
```

`cargo test` runs 42 tests: the rule checker, the test that no memory is
requested during the calculation, and the functional tests.

| Quantity | Checked against | Agreement |
|---|---|---|
| All output values for W7-X, HSX and a rescaled W7-X | the Python version | 1e-10 (the test requires 1e-8) |
| Plasma volume, area, radii | circular torus (analytic); VMEC's own values for W7-X | 1e-12; 1e-9 |
| Field of a filament | circular loop and straight wire (analytic) | 1e-4; 1e-6 |
| Field on axis of W7-X coils | the W7-X vacuum equilibrium | 5e-6 |
| Self-inductance | formula for a thin circular ring | 0.2% |
| Kernel of the rectangular bar | far-field limit and symmetry at the centre | 1e-5 |
| Fit of a1 and a2 | recovered from prescribed values | 1e-10 |
| Torus with 20 circular coils | analytic geometry, ideal toroidal solenoid | 1e-10 to 1% |
| Rescaling | scaling laws for length, field, current and inductance | 1e-8 |
| Full W7-X result | PROCESS built-in W7X-50 values | see docs/w7x_comparison.md |
| Full HSX result | consistency with its equilibrium; plausibility | see docs/hsx_check.md |

Tests that need the W7-X or HSX files report that they were skipped when the
files are absent.

The Python version on `main` also has analytic tests for the mutual
inductance, the hoop force, virtual work and the field inside a bar, and a
cross-check against SIMSOPT. These have not been ported to Rust; they apply
to the Rust version through its agreement with the Python results.

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

### In the Rust version only

- **NetCDF classic format only.** `wout` files in NetCDF-4 (HDF5) are
  rejected; convert them with `nccopy -k classic`. VMEC2000 writes the
  classic format by default.
- **Fixed upper bounds.** For example 512 coils, 1024 samples per coil and 12
  filaments across each side of the winding pack. Inputs beyond a bound are
  rejected. See `src/limits.rs`.
- **One processor core.** Threads request memory and are therefore not used.
- **Empty output file after a failed run.** The output files are opened
  before the calculation. After an error an existing file is unchanged; a new
  file is left behind empty.
- **Not all analytic tests have been ported**, see Verification.
- **The daily static analysis does not run yet.** GitHub starts scheduled
  runs only from the default branch; on this branch the analysis runs on
  every push.
- **The rules cover the program itself**, not the tests, the Python helper
  scripts or Rust's standard library.

### By design

- **Modular coils only.** Every coil must close poloidally around the plasma.
- **One coil current.** PROCESS scales a single current, so coil sets with
  unequal currents are reduced to their mean absolute current, with a
  warning. Coils with zero current are ignored.
- **Input formats.** VMEC `wout` and MAKEGRID only.
- **No external libraries.** No existing Rust library is written to these
  rules, so the readers and writers are own code.

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
