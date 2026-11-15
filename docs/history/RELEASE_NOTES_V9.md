# PlasmaROM-H2 v0.9.0 — Release notes

Release date: 2026-11-15. Scope: numerical verification and chemistry-acceptance audit.

## Corrected

1. **Independent training replay:** every *reference training reactor* is now re-simulated from the source Cantera phase during acceptance. The observed ignition classification, interpolated threshold-crossing time, peak dT/dt time and largest reporting time step must agree with the saved per-case record. Tampered, missing and failed replay evidence is rejected. Independent held-out reference/candidate replay from v8 remains active.
2. **Complete species-identity evidence:** the accepted reduced phase must preserve the source species' full Cantera species input (including thermo and transport declarations), composition, charge, size and thermo coefficients / valid-temperature range. Changing a NASA coefficient produces a different fingerprint. Reaction identities, kinetic constants, colliders, full reference hash and skeletal hash are still checked.
3. **Composition-aware training and holdout:** scenarios can explicitly include a third member `(T_K, P_atm, "H2:2,O2:1,N2:4")`. Duplicates are identified by normalized mole fractions, even if a mixture is reordered or uniformly rescaled. Varying H2/O2 stoichiometry and diluent is supported; input species are forced into the candidate roster and source-phase availability is checked. New JSON scenario protocol and CLI `--cases-file` added; the manually dispatched strict CI workflow also accepts a case file and reactor horizon. This does **not** demonstrate generalization until real chemical results are run.
4. **Consistency and reproducibility:** `pyproject.toml`, `CITATION.cff` and CLI/docs now identify 0.9.0. Includes an explicitly labeled version-specific *observed numerical environment pin list*, not a full lock of transitive dependencies.
5. **Automated tests:** adversarial examples for forged training ignition records, missing training timestamps, changed thermo coefficients, malformed mixture strings, composition-alias leakage, and failure-to-load data. New optional Cantera integration tests cover phase YAML round-trip and off-stoichiometric reactor input.

## Locally observed tests

- **131 passed, 0 failed, 5 skipped** in verbose full-suite pytest log (`results_v9/PYTEST_V9.log`). Five are specifically optional Cantera tests; Cantera could not be installed here because the Python package index was unreachable.
- An earlier non-scaling run completed **129 passed, 5 skipped** and the two scaling tests also passed separately. The full verbose run reached the same 131-pass/5-skip summary.
- `examples/reproduce_v9.py` reproduces manufactured-solution numerical checks and composition-aware case validation, not physical chemistry.

## Still not claimed / still missing

- **No confirmed executable reduced H2/O2 or NH3/H2 mechanism in this environment.** CI can run real Cantera tests on a connected Linux GitHub runner; a green unit test is not proof that the *separate* strict chemical-acceptance gate passed.
- No independently solved Boltzmann EEDF, experimental electron-impact cross sections, full discharge/charge-balance physics, published P-DRGEP, fully reacting CFD/DNS/LES or scaled HPC benchmarking.
- No measured PAC dataset or independent experimental combustion validation; no fuel-consumption or plasma-assisted ignition claim should be made from synthetic tables.

This code was developed with AI assistance and needs personal execution, inspection and understanding before representing it as one's research contributions.
