# PlasmaROM-H2 v0.7.0 — evidence-integrity and leakage patch

## Status

**103 passed, 3 skipped, 0 failed** in a clean local Python run. The three skipped tests require Cantera (not installed in this environment). This is **not** a validated plasma combustion model or a demonstrated successful chemical reduction.

## Fixed — acceptance is now fail-closed

- The JSON acceptance gate **requires an actual Cantera-loadable source and skeletal YAML**. A hash match alone is not sufficient. Missing Cantera or malformed YAML means rejection rather than provisional success.
- The gate compares phase counts to JSON declarations; checks that reduced species are a subset of source species and that elemental species compositions are unchanged.
- A real reduction in both species and reaction counts is required for this strict gate; smaller candidate lists or merely producing a YAML file are not evidence of accepted reduction.
- All holdout trajectory errors are reaggregated from per-case values; a claimed maximum smaller than a case's actual value is rejected.
- Per-case ignition flags must agree with `ignition_error_evaluable`, classifications and reported counts. Ignition errors must be nonnegative and finite. No observed ignition => no positive scientific acceptance.
- Training ignition evidence is stored per training condition and recomputed; aggregate training counts are checked against those records.
- Training and holdout conditions, including near-identical floating-point duplicates, must be distinct within and across splits. Corresponding condition labels must match metric rows.
- The `HOLDOUT_PASS` status is revoked and the holdout-pass YAML removed if independent scientific acceptance rejects the provisional candidate. A failed acceptance remains documented in `scientific_acceptance_failures`.

## Fixed — data leakage

- `evaluate_channel_pruning` now rejects repeated or near-equal training/holdout powers and within-split duplicates **before running simulations**.
- `run_skeletal` now rejects repeated or near-equal training/holdout `(T, P)` pairs **before importing optional Cantera**, enabling deterministic regression tests even without that dependency.

## Added tests

- Forged self-consistent-looking JSON + hash-matching invalid YAML must fail.
- Contradictory ignition states, negative ignition errors, contradictory count totals, understated aggregate error and mislabeled cases must fail.
- Spoofed source/reduced mechanism counts and species identities must fail (mock metadata in unit tests).
- Overlapping and duplicated conditions in both simulation tracks must fail.
- A provisional `HOLDOUT_PASS` rejected by a final gate is downgraded.
- An optional real-Cantera integration test checks rejection of hash-matching but invalid YAML. This test is **skipped locally**, pending a Cantera environment.

## Compatibility note

Version 6 summary JSON files deliberately **do not qualify** for version 7 acceptance if they lack per-training-case ignition observations or valid Cantera phase evidence. Rerun `skeletal` with v7 after installing Cantera; do not manually invent the missing evidence.

## Current research limitations — NOT FIXED

- No real electron-impact cross-section dataset bundled; sample collision and EEDF tables are explicitly synthetic.
- No self-consistent Boltzmann EEDF solution, charge evolution, ionization mechanism, or coupled electric-field solver.
- No independently reproduced published P-DRGEP on a large plasma state-resolved mechanism.
- No successfully validated Cantera chemical reduction on this machine; no experimental flame/ignition reference data.
- No complete reactive 3D CFD, DNS/LES, or demonstrated HPC speedup.

**Scope of v7:** evidence integrity, regression testing, scientific-claim hygiene. Acceptance checks consistency and a real Cantera phase, but do not independently rerun every kinetic trajectory or prove experimental accuracy.
