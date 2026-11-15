# Verification V8 — executed in the local environment

- Source: `PlasmaROM-H2-v7_GitHub_Ready.zip` extracted and modified into v8.
- Python modules tested with installed numpy / scipy / pytest; **Cantera unavailable**.
- Full `python -m pytest -q -rs --disable-warnings`: **112 passed, 3 skipped, 0 failed** in 34.73 seconds.
- `python -m pytest -q -k "not scaling"`: 110 passed, 3 skipped, 2 deselected; 0 failed.
- `python -m pytest -q tests/test_v6_scaling.py`: 2 passed, 0 failed.
- Union: 112 passed, 3 skipped; no locally verified Cantera mechanism.
- `python -m pytest -q tests/test_v8_scientific_hardening.py tests/test_v7_evidence_integrity.py`: 35 passed.

## Two ignition clocks (constructed verification curve, not measured combustion)

- Peak dT/dt numerical time: 2.0 s.
- Interpolated 20 K temperature-threshold crossing: 3.857142857 s.
- The difference is intentional and now explicit in each case report.

## Reversible Maxwellian rate quadrature

- `k_reverse / k_forward = (g_lower/g_upper) exp(DeltaE/Te)` using shared threshold-aware quadrature.
- On synthetic rate tables, maximum ratio relative error was measured at machine precision.
- Energy conservation with the updated rates is still tested, but the bath/electron model is a prescribed closure, not a discharge solver.

## Fail-closed chemistry review

`assess-skeletal` performs source SHA-256 verification, Cantera phase reload, species stoichiometry, complete reaction identity multiset, and independent replay of all reported holdout conditions. The optional Cantera tests skipped on this environment; a physical pass is NOT claimed.

The separate `.github/workflows/scientific-acceptance.yml` is designed to run after GitHub publication and will not issue success if an acceptable reduction/ignition-evidence bundle is absent.
