# PlasmaROM-H2 v0.5 — state-resolved kinetics & mechanism-reduction prototype

**Scientific status:** Independently executable, educational numerical research prototype. **Not** a calibrated plasma-assisted combustion solver or a validated implementation of the VUB group's P-DRGEP. Plasma cross sections in the sample datasets are deliberately **SYNTHETIC**. The physical H₂ vibrational energies come from NIST. This repository was developed with AI assistance; anyone submitting it as their own work should independently execute, review, understand and accurately describe their contributions.

## What v5 adds

1. **Power-coupled electron energy:** a positive Maxwellian electron-temperature state governed by an ODE with prescribed input power; H₂(v=0,1,2) populations, vibrational energy and neutral heating obey a closed numerical energy balance. Unlike v4, electron Te is **not fixed**. Density and gas temperature remain prescribed, so there is no electric-field or discharge model.
2. **Cross-section coverage checks:** `maxwell_rate()` rejects insufficient energy-domain coverage based on an upper Maxwellian energy-tail fraction. This is a *diagnostic* distribution-tail bound, **not** a guaranteed relative error bound on rate coefficients when the cross section beyond the table is unknown.
3. **Non-Maxwellian external EEDF adapter:** reads `energy_eV,probability_per_eV`, checks normalization and coverage, evaluates `<sigma v>`, and records data sources and SHA256 hashes. **Does not solve a Boltzmann EEDF** or take raw EEPF format.
4. **Executable synthetic energy-channel reduction:** ranks two excitation links by integrated gross excitation power, disables a link, reruns the coupled system, and evaluates trajectory error on unseen power settings. This is **not** published P-DRGEP, nor a reduction of a large mechanism.
5. **Reversible-flux diagnostic:** separate DRGEP-adjacent importance option avoids gross reaction-flux cancellation near equilibrium (diagnostic only).
6. **Cantera-related integration hooks:** optional reactor/skeletal submechanism (from v4), optional **1D chemical H₂/air flame** via `FreeFlame` (new, *not run here*), and a strict CLI failure option if skeletal selection doesn't pass holdout.
7. **Local task-parallel demonstration:** check result equivalence and time several independent runs across CPU processes. This is *not* an HPC scalability result; the experiment may show slowdowns.

## Reproducible results

The current verified evidence is in **[`results_v5/REPORT_V5.md`](results_v5/REPORT_V5.md)**, with source JSON files and PNG charts. `results_v4/` preserves historical baseline output, not the v5 benchmark. At release creation, **61 tests passed, 2 Cantera-dependent tests were skipped**, zero failed (see instructions to reproduce the precise state on your machine).

![Synthetic three-level H2 and dynamic electron Te](results_v5/electron_energy.png)

## Setup (Linux Mint / Ubuntu)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
```

Run the actual v5 examples:

```bash
python -m plasmarom electron-energy --out results_v5
python -m plasmarom channel-reduce --out results_v5
python -m plasmarom verify --out results_v5
python -m plasmarom parallel-check --parallel-workers 2 --out results_v5
python -m plasmarom report-v5 --out results_v5
```

`parallel-check` compares **independent parameter-sweep tasks**, not parallel acceleration of the ODE solver. Do not assume a speedup.

## EEDF / cross-section inputs

External data must have **clear provenance**; actual data validity is not certified by this software.

- EEDF input: `energy_eV,probability_per_eV` where the energy probability density approximately integrates to one (up to 0.5% tolerance). It is not an EEPF file. The EEDF should be computed independently using a suitable Boltzmann solver or derived from a justified assumption.
- Cross section input: `energy_eV,sigma_m2`. The energy axis starts at zero, is strictly increasing, and covers both the process threshold and the important EEDF tail.
- The integration uses an energy probability density and `v(E)=sqrt(2eE/me)`.
- Files named `examples/data/SYNTHETIC_*` are **invented mathematical demo inputs, NOT measured cross sections**.

Sample, reproducible CLI:

```bash
python -m plasmarom eedf-rate \
 --eedf examples/data/SYNTHETIC_Maxwellian_EEDF_Te2eV.csv \
 --sigma01 examples/data/SYNTHETIC_H2_v0_v1_cross_section.csv \
 --source-eedf 'synthetic Maxwellian demonstration' \
 --source01 'synthetic excitation demonstration' \
 --out results_v5
```

For a real dataset, use traceable references and verify units, reaction thresholds, applicability and licensing first. A downloaded raw BOLSIG+ EEPF is **not** necessarily a normalized energy probability density; conversion/normalisation must be independently verified.

## Cantera 3.2 integration — not runtime-verified in this release

```bash
python -m pip install -e '.[test,cantera]'
python -m pytest -q tests/test_cantera_optional.py
python -m plasmarom cantera --out results_cantera
python -m plasmarom skeletal --mechanism h2o2.yaml --out results_skeletal
python -m plasmarom skeletal --mechanism gri30.yaml --out results_gri30_h2 --require-holdout-pass
python -m plasmarom flame-1d --mechanism h2o2.yaml --out results_flame
```

`--require-holdout-pass` **exits with error** unless the actual reduced mechanism passes held-out reactor-trajectory errors. Without it, the code may legitimately produce `NO REDUCTION` or `FAILED_HOLDOUT`. Don't call a candidate YAML *validated* if it does not pass. The small H₂/O₂ mechanism may not admit meaningful reduction. `gri30.yaml` is methane-oriented and not a dedicated hydrogen reference mechanism. Flame modelling uses **Cantera's solver**, not a newly implemented CFD solver; it includes chemical energy and species transport but **no plasma**. Grid, literature and experimental validation are still required.

Cantera is **not present in the runtime that built this release** because external package downloads were unavailable, so there are no measured reduced-species counts, validated ignition delays, flame speeds, or chemical simulation plots in this ZIP.

## What is still missing scientifically?

| Research requirement | v5 status |
|---|---|
| Numerical-method verification | **Executed**: periodic FV spatial convergence, conservation, stiff ODE tests |
| Electron-temperature energy equation | **Executed**: 0D fixed-density Maxwellian closure and explicit input power |
| Imported EEDF from external solver | **Input+rate code executed** using **synthetic** PDF; no solved EEDF |
| State-to-state H₂ vibrational kinetics | Three levels, real energies, **synthetic** cross sections |
| Energy-aware reaction pruning | Small synthetic model passes specified holdout; **not** P-DRGEP |
| Real chemical Cantera kinetics | Code and optional tests only: **not executed here** |
| Physical reactive CFD / 1D flame | Optional Cantera FreeFlame code only: **not executed here** |
| 3D DNS/LES and plasma-discharge physics | **Not implemented** |
| Experimental PAC validation | **Not implemented**; no trustworthy measurements bundled |
| HPC strong/weak scaling | **Not implemented**; local multi-process smoke only |

## Primary references

- [NIST H₂ diatomic experimental constants](https://cccbdb.nist.gov/expdiatomicsx.asp)
- [Cantera 3.2 submechanism construction](https://cantera.org/3.2/examples/python/kinetics/extract_submechanism.html)
- [Cantera 3.2 mechanism reduction](https://cantera.org/3.2/examples/python/kinetics/mechanism_reduction.html)
- [Cantera 3.2 FreeFlame numerical example](https://cantera.org/stable/examples/python/onedim/adiabatic_flame.html)
- [pyMARS DRGEP discussion](https://niemeyer-research-group.github.io/pyMARS/theory.html)
- [LXCat cross-section data access](https://us.lxcat.net/)

All code and sample tests are provided to demonstrate methods, **not** to claim completed research validation. This README is a scientific disclosure of the limitations, not merely an installation guide.
