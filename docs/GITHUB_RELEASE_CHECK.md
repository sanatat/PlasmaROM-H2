# GitHub publication check (2026-11-15)

The extracted source prepared for publication was tested with:

```bash
python -m pytest -q -m 'not slow'
```

**Observed locally:** `160 passed, 6 skipped, 2 deselected in 44.23s` (exit code 0).

- The six skipped tests rely on unavailable optional Cantera integration; they are **not successes**.
- The two deselected tests in `tests/test_v6_scaling.py` exercise local multiprocessing and are available in the manually triggered scaling workflow. They were **not** part of this check.
- The command tests code and numerical invariants, **not** physical plasma-assisted combustion or experiment.
- See [`evidence/PYTEST_GITHUB_FAST.log`](evidence/PYTEST_GITHUB_FAST.log) for the actual test output. `tests/test_v12_regressions.py` and 5 related scientific hardening suites were also run together: 85 passed.

The ZIP was rebuilt with a fresh SHA-256 manifest after README/CI documentation changes. Repository links, syntax, archive integrity and manifest verification were checked. See the release chat for the final archive-verification outcome.
