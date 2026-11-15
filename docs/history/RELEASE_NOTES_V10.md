# PlasmaROM-H2 0.10.0 — code fixes and scientific-claim controls

Release date: 2026-11-15.

## Fixed and verified with local numerical tests

- **Third-body collider omission (gross-flux):** `_reaction_participation_matrix` now supplies identical chemistry participation to both net DRGEP and gross bidirectional-flux calculations, including explicit collider efficiencies and generic M. A mock third-body case verifies that gross connectivity changes from 0 to 1 when a collider is included. Conservative generic M behavior is intentional: this code prioritizes avoiding false kinetic omissions, potentially at the expense of reduction.
- **Trace OH normalized-error blind spot:** `_curve_error` reports `OH_Y_max_abs_difference`, `OH_Y_activity_threshold`, `OH_Y_active_points`, and `OH_Y_trace_relative_max_error` in addition to existing normalized error. The candidate with OH=5e-9 against a 1e-12 reference exceeds the trace relative-error criterion, even when the older 1e-7 normalization gives only ~0.05. This is a numerical diagnostic with a stated 1e-12 mass-fraction floor.
- **Training candidate and reference replay:** a selected candidate stores training-case metrics; scientific acceptance compares fresh reference and skeletal Cantera results for all training and held-out cases and checks training aggregates. The code cannot claim training accuracy based solely on reference ignition labels.
- **Locked evaluation:** `v10_h2o2_default` and `v10_h2o2_composition` are SHA-256-pinned, versioned, non-overlapping evaluation plans. Parameters cannot be changed through CLI while retaining an accepted v10 protocol identifier. The strict CLI / report gate requires a pinned v10 result, a validated chemical file and independent replay. The plan is project-versioned; it was not preregistered with a third party.
- **CI distinction:** the normal workflow tests numerical behavior and a Cantera integration smoke; the separately dispatched `scientific-acceptance.yml` uses only pinned protocol choices and intentionally fails if genuine held-out chemistry evidence does not pass.
- **Metadata:** package and citation versions updated to 0.10.0, with JSON plans packaged for standard Python builds.

## Locally verified results

- Non-scaling suite: 140 passed, 5 skipped (Cantera unavailable), 0 failed.
- Separate scaling module: 2 passed, 0 failed.
- Synthetic / mocked v10 regression cases: 11 passed, 0 failed (subset of the non-scaling suite).
- The complete all-at-once pytest run stalled at the multiprocessing module; it is **not** represented as a successful complete run. Both pieces were tested separately.

## Still unresolved (do not claim as completed)

- No installed Cantera and no locally observed successful reduction of a real mechanism; optional chemistry behavior is not integrated-test verified here.
- No independently computed Boltzmann EEDF, real collision tables from LXCat with permissions/provenance, validated plasma discharge and charge dynamics, or physically validated state-resolved kinetics.
- No implementation validated as plasma-specific P-DRGEP, no published NH3/H2 plasma chemistry benchmark, no experimental flame or ignition comparison.
- No fully coupled reacting-flow CFD, 3D DNS/LES, or demonstrated HPC speedup.

The project is a numerical-computing prototype with reproducible software verification, not a validated research solver. AI-assisted code should be read, reproduced and understood by anyone who represents it as their work.
