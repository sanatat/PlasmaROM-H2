# PlasmaROM-H2 v12.0 — audited numerical research prototype

**Changes relative to v11:** [Release Notes v12](RELEASE_NOTES_V12.md) and [Verification v12](VERIFICATION_V12.md).

Four verified software-level corrections: (1) `CO` consumes 0.5 mol O₂ per mol CO in the simplified NH₃/H₂/CO equivalence-ratio calculation; unknown species fail closed, (2) duplicate `(case_id, observable)` independent-measurement rows are rejected rather than double-counted, (3) signed superelastic electron-energy exchange is supported and positive losses and negative gains must **each** meet the retained-energy target, and (4) GitHub Actions' real-chemistry workflow now triggers on changes throughout the model, test and example trees and on pull requests.

**Evidence boundary:** numerical regression tests passed locally. Cantera and real two-term EEDF were **not run** here. These fixes do not constitute validated plasma-assisted combustion, published P-DRGEP, experimental confirmation or DNS/LES.

```bash
python -m pip install -e '.[test]'
python -m pytest -q tests/test_v12_regressions.py
python -m examples.run_energy_graph_v11 # synthetic, not physically validated
# Internet-connected system for optional real integration:
bash run_verified_pipeline.sh
```

**Energy-sign interpretation:** in `energy_aware_screen`, positive `energy_losses_eV` signifies energy transferred *out of* electrons; negative entries denote energy *gained by* electrons via superelastic exchange. The two retention quotas are enforced separately; the historical `retained_energy_fraction` now means retention of total **absolute** exchange, not the cancellation-prone net transfer. This remains a synthetic graph-screening method, not plasma P-DRGEP.

**Experimental-replicate policy:** each independent measurement must use a unique `(case_id, observable)`. Repeated experimental replicates must have separate `case_id` values and matching model predictions for the distinct experimental conditions/replicate records. Apparent duplicates are rejected to prevent accidental double-weighting; this is a data-integrity rule, not scientific verification that measurements were independent.

---

# Archived v11 documentation — executable PAC research prototype

**IMPORTANT:** The core numerical code has been tested locally. Optional integration with Cantera 3.2 real NH3/H2 chemistry and its **actual Boltzmann two-term solver** is implemented, but **not yet locally executed** due to network installation restrictions. This is **not** a completed or experimentally validated plasma-assisted-combustion solver, P-DRGEP implementation, 3-D DNS or LES. Do not claim these uncompleted accomplishments in an application.

## Fastest path with internet on Linux Mint

```bash
bash run_verified_pipeline.sh
```

The script installs Cantera and requires the real 3.2 numerical runs to succeed; it stops on error. Evidence is saved in `results_v11/ammonia_h2/`, `results_v11/cantera_boltzmann/`, and `results_v11/transport/`. A successful code run is **not experimental validation**.

For offline verified numerical tasks:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
pytest -q tests/test_v11_extensions.py
python -m plasmarom transport-1d --out results_v11/transport
python examples/run_energy_graph_v11.py
```

## Scientific scope and provenance

- **NH3/H2 numerical chemistry**: `ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml', 'baseline')`, literature-derived example model. Different from plasma-assisted combustion.
- **True numerical EEDF**: built-in Cantera 3.2 `Boltzmann-two-term` method for **N2/O2** with the separately distributed Phelps collision set. No H2/NH3 EEDF or feedback coupling claim.
- **Physics-inspired graph screen**: keeps high-weight kinetic paths and integrated energy channels in a synthetic graph; it is explicitly not published P-DRGEP.
- **Transport**: conservative 1-D constant-density transport and heat release proxy, **not** compressible reacting CFD.
- **Experimental comparison**: schema in `examples/data/`; no measured values are bundled or invented.

See [gap-by-gap scientific audit](GAP_STATUS_V11.md), [release notes](RELEASE_NOTES_V11.md), and [verification](VERIFICATION_V11.md).

---

# PlasmaROM-H2 v0.10.0 — pinned numerical chemistry verification prototype

**Scope:** Independent educational numerical-computing and chemistry mechanism-screening code. **Not a validated plasma-assisted-combustion solver**, not published plasma P-DRGEP, not a completed NH₃/H₂ ignition study, and not DNS/LES. All bundled cross sections and EEDF examples remain **SYNTHETIC**.

## Version 10 — new correctness and reproducibility controls

1. **Gross-flux Third-Body correction.** Net and gross chemical DRGEP use the same collider-aware participation matrix. Explicit collider names and efficiency-listed species count as participants; generic M with a nonzero default efficiency conservatively includes all species. This can legitimately result in **NO REDUCTION** rather than silently discarding a kinetically relevant collider.
2. **OH trace-species integrity.** Alongside the old normalised OH trajectory error, v10 records maximum absolute OH mass-fraction difference and a relative error with a fixed **1×10⁻¹² mass-fraction floor** at active points. A 5000-fold OH false-positive signal is no longer hidden by the older `1e-7` normalisation floor. This deliberately conservative metric may reject many candidate mechanisms and is **not experimental validation**.
3. **Independent Training candidate replay.** The accepted skeletal mechanism must reproduce the **reference and candidate** curves at **all training and holdout conditions**. Training error metrics and the worst training error are saved and checked independently, not just training ignition flags.
4. **Pinned evaluation protocol.** Two committed H₂/O₂ numerical plans in `protocols/` specify conditions, observation horizon, snapshots, thresholds, OH floor, and error tolerances. Their SHA-256 digests are pinned in `plasmarom/protocol.py`. Changing a plan requires a code review/release, rather than silently changing a threshold in the report. This is a *repository-level reproducibility control*, not independent third-party preregistration.
5. **Fail-closed acceptance.** `assess-skeletal` demands a pinned v10 plan. Unpinned exploratory outputs, missing Cantera, unreplayed calculations, fake chemistry files, inconsistent case errors, and non-reduced mechanisms **cannot receive scientific acceptance**. Normal code tests and smoke CI are separate from the manually triggered scientific-acceptance gate.

### Run numerical tests on Linux Mint

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q --ignore=tests/test_v6_scaling.py
python -m pytest -q tests/test_v6_scaling.py
```

### Optional chemistry run (requires separately installed Cantera 3.2+)

```bash
python -m pip install -e '.[test,cantera]'
python -m pytest -q tests/test_cantera_optional.py -rs
python -m plasmarom skeletal --mechanism h2o2.yaml \
  --protocol v10_h2o2_default --out results_skeletal_v10
python -m plasmarom assess-skeletal --out results_skeletal_v10
```

Or use the fixed composition-sweep plan:

```bash
python -m plasmarom skeletal --mechanism h2o2.yaml \
  --protocol v10_h2o2_composition \
  --cases-file examples/cases/h2_composition_sweep.json \
  --out results_composition_v10
python -m plasmarom assess-skeletal --out results_composition_v10
```

**The last command is expected to exit with an error unless a genuine chemical mechanism is reduced and passes all held-out tests.** Successful installation, unit tests, a written YAML, or a successful smoke run alone are *not* evidence of an accepted chemical reduction. In particular, the small `h2o2.yaml` system may have no useful valid reduction under a conservative collider-aware screen.

### Verification in this release

Locally executed without Cantera: **140 numerical/decision tests passed and 5 optional Cantera tests skipped**, excluding the separate scaling module; **2 scaling tests passed separately**. See `VERIFICATION_V10.md` and `results_v10/` for test logs. Local Cantera installation failed due to an unreachable Python package index. **No accepted chemical reduction, PAC physical result, experimental comparison or HPC speedup is claimed.**

For a technical breakdown, see [`RELEASE_NOTES_V10.md`](RELEASE_NOTES_V10.md).

---

## Archived version 9 documentation (historical only)

# PlasmaROM-H2 v0.9 — historical documentation

**Scope:** Independent educational project for numerical ODEs/PDEs, synthetic state-resolved H₂ dynamics, cautious chemistry reduction and reproducibility checks.

> **NOT** a validated plasma-assisted-combustion solver. No complete state-to-state combustion mechanism, experimental plasma benchmarking, validated P-DRGEP, physical 3-D DNS/LES, or demonstrated HPC speedup. **All bundled cross-section and E/N demonstration tables are SYNTHETIC**, not measured or physically derived. The code was developed with AI assistance; users should personally execute, inspect and understand it before describing their own contributions.

## Version 9 — independent training replay, species identities, and mixture protocol

The version 9 **scientific acceptance** gate now replays both the original **training reference** and held-out reference/candidate reactor trajectories when Cantera is available. It checks each training ignition report against freshly calculated flags and two explicit time clocks, and compares all Cantera species thermo/transport input and NASA coefficients rather than composition alone. The reduced phase must still match each source reaction's kinetic identity.

A new per-case composition field expands tests beyond fixed T/P. Cases `[T_K, P_atm, "H2:2,O2:1,N2:4"]` and legacy `[T_K, P_atm]` are supported. Examples include `examples/cases/h2_composition_sweep.json`:

```bash
python -m plasmarom skeletal --mechanism h2o2.yaml \
  --cases-file examples/cases/h2_composition_sweep.json \
  --out results_skeletal_v9
python -m plasmarom assess-skeletal --out results_skeletal_v9
```

**Offline status:** 131 passed, 0 failed, **5 optional Cantera tests skipped**. The Python package index was unreachable; there is **no locally confirmed real chemical mechanism reduction or physically validated plasma discharge**. The strict assessment command is expected to fail if no accepted reduced mechanism exists; no pass was fabricated. Numerical verification can be reproduced independently:

```bash
python examples/reproduce_v9.py
```

See [`RELEASE_NOTES_V9.md`](RELEASE_NOTES_V9.md), [`VERIFICATION_V9.md`](VERIFICATION_V9.md), and [`results_v9/REPRODUCE_V9.json`](results_v9/REPRODUCE_V9.json). Version 8 and earlier sections are historical records.

## Version 8: corrected microreversibility, ignition clocks and chemical audit

**Locally executed: 112 passed, 0 failed, 3 skipped Cantera tests** in the complete pytest suite (34.73 seconds); core and scaling tests also passed separately. Cantera could not be installed due to network restrictions, so physical chemistry reduction remains **unverified**.

- Consistent shared-grid detailed-balance rates now apply to both three-state and energy-coupled models, including non-degenerate Maxwellian limit and configurable level degeneracy.
- Ignition reporting separates the interpolated 20 K threshold time from the peak dT/dt time and records sampling resolution and right-censoring status.
- Independent Cantera acceptance checks reaction equations, all rate parameters, Third-Body data and reversibility, not just species/reaction counts; it also pins the source SHA-256.
- Independent acceptance **replays** each holdout reactor and compares new numerical error metrics. Corrupted inputs or failed replay reject the result.
- New adversarial tests exercise forged reaction data, source changes and mismatched replay outputs.

See [`RELEASE_NOTES_V8.md`](RELEASE_NOTES_V8.md), [`VERIFICATION_V8.md`](VERIFICATION_V8.md) and [`results_v8/REPORT_V8.md`](results_v8/REPORT_V8.md). The older v7 report is preserved below as historical context.

## Version 7: verified scientific-claim safeguards

**Local tests: 103 passed, 0 failed, 3 skipped (Cantera not installed).**

This release **fixes the three reproducible v6 integrity bugs**:

1. **Invalid YAML and self-reported results can no longer pass on a matching checksum alone.** Accepted artifacts require a Cantera-loadable source and skeletal mechanism with independently verified phase counts, source-species subset, and matching species composition. If Cantera cannot run, acceptance is **rejected**, never certified.
2. **Duplicate/overlapping training and holdout conditions are rejected**, both for the synthetic energy-channel workflow and the optional Cantera reduction. Near-equal float duplicates also fail.
3. **Scientific summary statistics are recomputed** from labeled per-case evidence; false ignition flags, impossible counts, negative errors or understated maximum errors fail.

If the provisional chemical trajectory check passes but the independent acceptance check fails, `HOLDOUT_PASS` is removed and the supposed holdout-pass YAML is deleted. Old v6 summary JSON files without training-case evidence cannot be retroactively marked as valid.

See [`RELEASE_NOTES_V7.md`](RELEASE_NOTES_V7.md) and [`VERIFICATION_V7.md`](VERIFICATION_V7.md) for exact tests, results, and remaining gaps. The core physics and validation limitations listed below **still apply**.

## Version 6 changes

| Previous gap | v6 corrective action | Evidence / limits |
|---|---|---|
| Unobserved ignition silently counted as low error | Right-censored ignition statuses, classification-mismatch penalties, and mandatory observed reference ignition for acceptance | Unit tests; **Cantera not run locally** |
| A CI green check could imply false reduction success | Separate *Cantera smoke* from fail-closed `assess-skeletal`; a manual scientific acceptance workflow fails unless actual held-out evidence exists | Unit-tested acceptance; CI workflow not executed in this local environment |
| Important reversible chemistry neglected in selection | Combine chemical net-rate DRGEP with conservative gross bidirectional-flux graph for selecting species | Synthetic graph and mock-reaction regression tests; NOT specialized plasma P-DRGEP |
| Old output YAML could be confused with a fresh valid model | Remove stale files before analysis, require reloadable Cantera YAML and hashed holdout-pass artifact | Code implemented; optional integration unverified without Cantera |
| EEDF had no E/N input mode | Import multiple externally solved EEDFs, interpolate *within* a documented E/N grid, validate normalization and cross-section coverage, hash provenance | Verified on **synthetic** external-table example; NOT a Boltzmann solver |
| A single parallel wall-time was presented | Multiple local fixed-workload runs, serial equality, medians and efficiency | Task-level scaling only; not HPC scaling |
| ML fit could be favored after looking at held-out cases | Disjoint training / validation / held-out sets, model selected by validation with baseline-QSSA fallback | Verified on toy dynamics, NOT physical ML-based lumping |
| 1D flame lacked grid sensitivity check | Independent coarse/fine Cantera solves and speed-change tolerance gate | Tests with mocked solver only; physical Cantera result not verified locally |

## Historical v6 outputs (retained for transparency)

- **Historical v6 result: 77 successful, 0 failed, 2 skipped.** Historical test log: [`results_v6/pytest_v6.txt`](results_v6/pytest_v6.txt).
- `results_v6/field_dependent_rate.json`: synthetic rate coefficient for invented 100 Td mapping; not a published collision coefficient.
- `results_v6/ml_validation_fallback.json`: validation-chosen model, with untouched held-out conditions.
- `results_v6/fixed_workload_scaling.json`: local 1 vs 2 CPU process timing; no speedup was demonstrated.
- Full discussion: [`results_v6/REPORT_V6.md`](results_v6/REPORT_V6.md) and [`RELEASE_NOTES_V6.md`](RELEASE_NOTES_V6.md).

## Run in Linux Mint / Ubuntu

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python examples/reproduce_v6.py
python -m plasmarom verify --out results_verify
```

Reproducible external-EEDF *interface demo* (all example data synthetic):

```bash
python -m plasmarom eedf-field \
  --eedf-manifest examples/data/SYNTHETIC_EoverN_manifest.csv \
  --sigma01 examples/data/SYNTHETIC_H2_v0_v1_cross_section.csv \
  --source01 'SYNTHETIC invented cross section' \
  --field-Td 100 --out results_field
```

Use a real EEDF only after independently solving for your gas composition and E/N (e.g. with BOLSIG+) and verifying energy-probability-density units. Do **not** ingest raw electron-energy probability functions (EEPFs) without a correctly documented conversion. **Field interpolation is not a discharge solution.**

## Optional real Cantera chemistry (requires Cantera installation)

```bash
python -m pip install -e '.[test,cantera]'
python -m pytest -q tests/test_cantera_optional.py
python -m plasmarom skeletal --mechanism h2o2.yaml --out cantera_attempt
python -m plasmarom assess-skeletal --out cantera_attempt
python -m plasmarom flame-grid --mechanism h2o2.yaml --out flame_mesh_study
```

`assess-skeletal` **intentionally exits nonzero** if Cantera is absent or if there is no validated reduction with observed ignition, acceptable held-out error, source/artifact hashes, unchanged reaction identities and independently replayed reactor curves. The small `h2o2.yaml` mechanism may not permit a useful skeletal reduction: **a failed acceptance is a valid negative scientific result.** `flame-grid` checks mesh sensitivity, not experimental agreement. The reference 1D solver comes from Cantera; it is not a solver we coded from scratch.

CI distinction: `.github/workflows/ci.yml` runs numerical tests and a Cantera *smoke* job but makes **no chemical accuracy claim**; `.github/workflows/scientific-acceptance.yml` is an explicitly dispatched **failing** gate unless actual independent evidence passes. Neither GitHub workflow was run locally during release creation.

## Research gaps not closed

1. **Real electron-impact data** and associated citation/license/provenance, including excitation and ionization energy channels.
2. **Self-consistent EEDF/Boltzmann and discharge solution**: field, charge balance, electron number density evolution, chemistry feedback and gas temperature are not jointly solved.
3. **Published plasma P-DRGEP** on a large state-resolved mechanism, with real energy partition channels and held-out combustion conditions.
4. **Combustion evidence**: Cantera integration needs actual runtime, measured ignition/temperature/radical error, and accepted reduced YAML; none are present in this ZIP.
5. **Reactive CFD/DNS/LES**: no 3-D reacting flow, turbulence benchmark, variable-density plasma flame or validated DNS/LES.
6. **Experimental and HPC validation**: no calibrated PAC experiment, strong/weak solver scaling on a cluster, or measured large-network speedup.

These are research-grade development areas, not bugs that can be legitimately marked as solved by adding more toy tests.

## Data references and rights

- [Cantera 3.2 kinetics API](https://cantera.org/3.2/python/kinetics.html) and [mechanism-reduction examples](https://cantera.org/3.2/examples/python/kinetics/mechanism_reduction.html)
- [BOLSIG+ electron Boltzmann solver](https://www.bolsig.laplace.univ-tlse.fr/download.html)
- [LXCat data and contributor-specific terms](https://master.lxcat.net/instructions/how_reference.php)
- [NIST molecular constants](https://cccbdb.nist.gov/expdiatomicsx.asp)

**Do not redistribute third-party LXCat collision tables without checking the specific contributor's terms.** Real data have not been bundled or misrepresented.

## Appropriate CV wording

*Independent, AI-assisted computational research prototype; implemented and tested numerical verification, graph-based chemical-screening safeguards, strict held-out validation, and synthetic EEDF/kinetic workflows.* Document only your personally reviewed and reproducible contributions.
