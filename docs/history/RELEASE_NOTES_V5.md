# Release v0.5 — numerical verification and plasma-energy gap reduction

## Gaps narrowed

- Changes electron temperature from a fixed external parameter (v4) to an ODE state tied to electron energy, prescribed power input, vibrational excitation/deexcitation, elastic bath transfer and explicit energy accounting.
- Replaces silent short cross-section integration with Maxwellian energy-tail coverage rejection; cannot prove the unknown high-energy cross-section tail is negligible.
- Adds external EEDF rate integration and input normalization/domain/provenance checks. The sample CSV data are synthetic.
- Adds a small executable energy-channel reduction demonstration with separate training/held-out power conditions and explicit error criteria.
- Adds gross forward/reverse rate screening as a **diagnostic** alongside chemical DRGEP; not a P-DRGEP implementation.
- Adds optional Cantera 1D hydrogen/air chemical flame integration using the vendor's `FreeFlame` solver and a strict holdout acceptance CLI option for skeletal reaction reduction.
- Adds local process-pool parameter-sweep check. Speedups are measured rather than assumed; one run is not HPC scaling.
- Keeps v4 regression tests and new falsification-focused tests for poor EEDF tails, invalid distributions, strong/no power, and serial/parallel identity.

## Open gaps

Actual electron-impact cross-section source datasets, E/N-dependent Boltzmann EEDF, ionization/charge/field balance, specialized P-DRGEP, measured PAC experiments, full variable-density reacting flow, DNS/LES, HPC strong/weak scaling, and runtime validation of Cantera remain open. Do **not** claim these are solved.

## Reproducibility caveat

All numerical figures in `results_v5/` were generated on synthetic plasma data; optional Cantera is absent in this release environment. If Cantera runs on another machine, inspect and attach its JSON outputs separately. Do not mix statuses.
