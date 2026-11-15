# PlasmaROM-H2 v0.4 release notes — 2026-11-14

## Fixed and extended

- Carried forward v0.3 finite-volume conservative advection and consistent B=F1+F2 QSSA representation.
- Added manufactured periodic scalar advection/diffusion verification on three grids with measured convergence rates.
- Added a 0D 3-level H2 vibrational population balance using NIST vibrational constants and Maxwellian `<sigma v>` evaluation from table data. Added electron-transfer + vibrational + neutral-quenching energy balance.
- Added a dimensionless sensible-energy proxy coupled to conservative 2D species transport. This is **not** thermochemical CFD.
- Added bounded fitted QSSA-lag regression, explicitly separated training and held-out toy cases.
- Added explicit integrated gross electron-energy channel screening, expressly **not** P-DRGEP.
- Added *optional* Cantera-based construction of a skeletal chemical mechanism with conservative third-body handling, training-only selection, withheld-condition evaluation, and explicit NO REDUCTION / FAILED HOLDOUT statuses.
- Added independent data-validation tests and CI commands; core test suite: **44 passed, 2 Cantera-dependent tests skipped** on release authoring environment.
- Added generated evidence in `results_v4/`: reproducible JSON, Markdown and figures.

## Correct interpretation

The 0D three-level input energies are spectroscopy-based, but **the bundled electron-impact cross sections are synthetic**. Code handles third-body reaction participation conservatively but **Cantera runtime was unavailable**, so construction and reduction metrics could **not** be verified here. The repo does **not** contain validated PAC, a verified P-DRGEP, experimental cross sections, EEDF Boltzmann calculation, real flame CFD, DNS/LES, or HPC scaling.

## Suggested next scientifically defensible upgrades

1. Run Cantera 3.2 in a supported Linux environment and inspect generated `skeletal_summary.json`; repair runtime/API failures if they occur.
2. Obtain traceably licensed H2 electron-impact data from LXCat or an independently documented database; review the applicable EEDF conditions, uncertainties, detailed-balance assumptions, and compare with literature.
3. Add E/N-dependent electron kinetics using an appropriate Boltzmann solver (e.g. BOLSIG+) and a calibrated full kinetic mechanism before attempting true plasma P-DRGEP.
4. Validate one genuine hydrogen ignition/flame case against a published experimental dataset with uncertainty and holdout analysis.
5. Only then consider chemically reacting variable-density CFD and scientifically meaningful parallel strong/weak scaling.
