# PlasmaROM-H2 v0.5 — computed verification report

**Status:** educational, independently executable numerical portfolio prototype. NOT validated plasma-assisted combustion.

## Coupled electron-energy closure (synthetic cross sections)

- Initial / final electron Te: 2 / 0.294608 eV
- Deposited power energy: 0.6 J/m³
- Max relative total-energy residual: 7.089e-09
- Max population conservation residual: 9.992e-16
- Electron Te is solved as an ODE with constant externally prescribed input power and Maxwellian approximation; fixed ne and gas bath. **No field, E/N, Boltzmann solver, or charge balance.**

## Energy-aware reaction channel pruning (synthetic)

- Result: **HELDOUT_PASS**
- Retained: v0<->v1
- Removed: v1<->v2
- Gross excitation energy retained: 99.598297%
- Maximum held-out normalised trajectory error: 0.0022679 (limit 0.05)
- **This is not published plasma P-DRGEP.** There are only two transitions and no real collision dataset.

## External EEDF CSV interface

- Test data: explicitly **SYNTHETIC** Maxwellian EEDF and invented excitation cross section
- Computed cross-section rate: 6.264567e-15 m³/s
- EEDF normalisation integral: 0.99907562
- CSV SHA256 hashes and provenance are saved in `external_eedf_rate.json`. External energy distribution must already be computed; the software **does not solve EEDF**.

## Periodic numerical verification

- Upwind advection observed order: 0.9797 (≈ first order)
- Centered diffusion observed order: 1.9909 (≈ second order)
- These are manufactured-PDE correctness tests, not flame validation.

## Local serial versus two-process check

- Serial walltime: 1.579 s
- Parallel walltime: 2.029 s
- Observed speedup: 0.7781× (may be below one owing to overhead)
- Max difference between results: 0
- One uncalibrated run of independent tasks: **NOT an HPC strong/weak scaling experiment**.

## Cannot be marked verified in this environment

- `Cantera >=3.2` was unavailable; its reactor/skeletal mechanism and 1D H₂/air flame integrations were **not executed**.
- Without a licensed real electron-impact dataset and externally computed EEDF, physical plasma coefficients remain **unverified**.
- No plasma-specific P-DRGEP, full electron Boltzmann EEDF solver, 3D DNS/LES, calibrated reacting-flow CFD, or experimental PAC validation.
- Do not invent a real reduced mechanism, ignition delay, flame speed, measured error or speedup.
