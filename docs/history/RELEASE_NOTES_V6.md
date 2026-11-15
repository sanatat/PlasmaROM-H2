# Release notes — PlasmaROM-H2 v0.6

Version 6 focuses on scientifically meaningful safeguards rather than inflating test counts.

## Implemented and independently executed

- False-positive ignition protection with right-censoring and resolved numerical temperature-rise criteria.
- Candidate-selection hybrid maximum of net DRGEP and bidirectional gross-flux importance.
- Fail-closed chemical reduction acceptance with per-condition errors, ignition evidence, YAML hash and actual species-count reduction.
- Optional Cantera candidate YAML round-trip and stale-artifact cleanup (code written, **not executed offline**).
- E/N-dependent imported EEDF grid and no-extrapolation interpolation (synthetic tests only).
- QSSA-vs-lag ML baseline fallback selected on disjoint validation conditions.
- Repeated fixed-workload local process-pool benchmark with explicit negative speedup result.
- Optional two-grid Cantera flame speed sensitivity analysis (mocked tests only).
- Isolated manually triggered scientific acceptance CI workflow; routine smoke job is not a passing research-validation claim.

## Tests

77 passed; 2 skipped; zero failed in 29.54s on local Python 3.13.5. Skips correspond to Cantera 3.2, unavailable for installation in this offline environment.

## Remaining explicit open questions

Neither a real E/N-derived EEDF, accepted real chemical skeletal reduction, real plasma P-DRGEP, calibrated ignition dataset, nor physical 3D reacting CFD has been verified. Any AI or physics claim beyond demonstrated toy numerical methods should be treated as a future research goal.
