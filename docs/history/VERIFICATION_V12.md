# PlasmaROM-H2 v12 — verification log

Date: 2026-11-15.

**Actually executed on extracted v11 source modified to v12:**

- `python -m pytest -q tests/test_v12_regressions.py tests/test_v11_extensions.py` → **20 passed**.
- `python -m pytest -q tests/test_v12_regressions.py tests/test_v11_extensions.py tests/test_v10_integrity.py tests/test_v9_integrity_and_mixtures.py tests/test_v8_scientific_hardening.py tests/test_v7_evidence_integrity.py` → **85 passed**.
- Legacy numerical and core tests in 10 modules → **57 passed**.
- `python -m pytest -q --ignore=tests/test_v5_improvements.py --ignore=tests/test_v6_scaling.py` → **143 passed, 6 skipped** (skips are optional package integration, not evidence of success).
- `python -m pytest -q tests/test_v5_improvements.py -k 'not plasma_energy_report_round_trip and not v5_pipeline_end_to_end'` → **16 passed, 1 deselected**.

The above suites **overlap**; do not add the counts together without deduplicating. The complete all-tests run was interrupted by an execution time limit, and the v6 multiprocessing scaling module was not rerun in this verification. The singled-out `test_plasma_energy_report_round_trip` was also not rerun in the v12 verification. Complete logs of the two non-overlapping runs are in `results_v12/PYTEST_V12_FAST.log` and `results_v12/PYTEST_V12_V5_SUBSET.log`. Accordingly this is **not** a claim that all tests passed.

Cantera is not installed in this environment, so no real Cantera chemistry or Boltzmann solver execution is evidenced. GitHub Actions status cannot be certified by local unit tests.

**Synthetic regression assertions:** CO example φ=1.4; repeated measurement keys raise ValueError; a superelastic energy-gain channel remains selected and signed; workflow path trigger coverage is verified textually.

**Scientific limitation:** no data or experimental validation of plasma-assisted NH₃/H₂ combustion, plasma-specific P-DRGEP, DNS/LES, or HPC acceleration.
