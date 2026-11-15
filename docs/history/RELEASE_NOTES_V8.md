# PlasmaROM-H2 v8 — Scientific correctness and evidence hardening

Release date: 2026-11-15. Based on the v7 project, retaining synthetic inputs and all historical artifacts.

## Fixed and verified locally

1. **Numerically consistent detailed balance.** Introduced `maxwell_reversible_rates()` in `plasmarom/plasma_states.py` using a common threshold-aware quadrature. Excitation/de-excitation rates now satisfy the stated non-degenerate Maxwellian detailed-balance relationship to floating-point precision. Reused in the fixed-Te state model and in the power-driven electron energy ODE. **Not** a Boltzmann EEDF solver; all included cross-sections remain synthetic.
2. **Unambiguous ignition reporting.** `_ignition_assessment()` returns separately the interpolated 20 K threshold crossing time, the maximum dT/dt time (the historical `ignition_time_s` proxy), slope magnitude, resolution, and censoring status. Training and holdout summaries now disclose their numerical clocks. These definitions can disagree; neither is experimental ignition validation.
3. **Reaction identity verification.** The independent acceptance gate compares detailed Cantera reaction signatures: reaction equations, reactants/products, reversible/duplicate flags, orders, type, rate input_data, third-body collider, default efficiency and explicit efficiencies. It requires the entire reduced multiset to be a subset of the source multiset; a match on counts and hashes alone cannot pass.
4. **Pinned source mechanism.** Acceptance independently resolves and hashes the source mechanism YAML and compares it with the originating run; changes to the reference invalidate the report.
5. **Independent reactor replay.** The acceptance gate reruns each holdout reactor with source and skeletal YAML, same recorded time horizon and sample count, and compares the recomputed quantitative metrics with the report. Missing or conflicting results fail closed. Tests use metadata/replay mocks only to prove gate logic; physical Cantera runs still require Cantera.
6. **New adversarial regression tests.** Check changes to Arrhenius/third-body parameters, corrupted source hash, forged yet self-consistent aggregate data, replay integration failure, coarse-grid detailed balance and differing ignition clocks.
7. **CI and package metadata.** Project version advanced to 0.8.0. Both generic and optional Cantera workflows include the v8 regressions. The scientific-acceptance workflow remains a deliberately strict gate, which can fail when there is no acceptable reduced chemistry.

## Local validation

- **Full suite** `python -m pytest -q -rs --disable-warnings`: **112 passed, 3 skipped, 0 failed** in 34.73 s.

- `pytest -q -k 'not scaling'`: **110 passed, 3 skipped, 2 deselected**, ~28.95 s.
- `pytest -q tests/test_v6_scaling.py`: **2 passed**, ~7.27 s.
- Across these complementary runs: **112 unique passed, 3 Cantera integration tests skipped, 0 failures**.
- `pytest -q tests/test_v8_scientific_hardening.py tests/test_v7_evidence_integrity.py`: **35 passed** (subset of the above).
- Reproducible numeric outputs: `results_v8/VERIFIED_NUMERICAL_RESULTS.json`.

## Not verified, not claimed

Cantera cannot be installed in the current environment due to unavailable package-index network resolution. The code calls the documented Cantera 3.2 API, but successful execution and accepted mechanism reduction have **not** been demonstrated locally. Physics-grade P-DRGEP, Boltzmann EEDF, measured collision data, NH3/H2 chemistry, reacting CFD, DNS/LES and experimental validation are **not provided**. A passing numeric unit test suite does not close those gaps.

The exact published rate parameter equality check may be conservative for Cantera mechanisms containing custom rate extensions or user metadata; these must be inspected on an installed Cantera runtime. Failing closed in that case is safer than silently accepting altered chemistry.
