# Verification log — PlasmaROM-H2 v0.7.0

Executed on 2026-11-14 in the supplied environment.

```text
$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest -q -rs
s.ss.................................................................... [ 67%]
..................................                                       [100%]
SKIPPED [1] tests/test_cantera_optional.py:7: could not import 'cantera': No module named 'cantera'
SKIPPED [1] tests/test_cantera_optional.py:28: could not import 'cantera': No module named 'cantera'
SKIPPED [1] tests/test_cantera_optional.py:46: could not import 'cantera': No module named 'cantera'
103 passed, 3 skipped in 27.92s
```

The three skipped tests are **not** counted as successes and no real Cantera reduction, reaction trajectory or flame is claimed.

## Falsification tests

| Attack / defect | Expected outcome | Locally verified |
|---|---|---|
| Fake YAML with matching SHA-256 | Rejected: missing Cantera or unreadable mechanism | Yes |
| Holdout `candidate_ignition_observed=False` but flag says classifications match | Rejected | Yes |
| Negative ignition proxy error | Rejected | Yes |
| 99 claimed ignitions despite one scenario | Rejected | Yes |
| Claimed max error 0.001 with case error 0.03 | Rejected | Yes |
| Different condition label vs measured row | Rejected | Yes |
| Near-duplicate training and holdout excitation powers | `ValueError` | Yes |
| Duplicate or overlapping Cantera training/holdout conditions | `ValueError` without Cantera import | Yes |
| Candidate YAML altered after hashing | Rejected | Yes |
| Mocked inconsistent chemistry counts / species subsets | Rejected | Yes (mocked phase metadata) |
| Provisional `HOLDOUT_PASS` with rejected final gate | Downgraded; validated YAML removed | Yes (mocked final gate) |

Tests that mock Cantera metadata verify the **decision logic only**, not the validity of a real chemical mechanism. The optional integration test will check YAML parsing in an environment with Cantera installed.

## Running with real Cantera

```bash
python -m pip install -e '.[test,cantera]'
python -m pytest -q tests/test_cantera_optional.py -rs
python -m plasmarom skeletal --mechanism h2o2.yaml --out results_skeletal
python -m plasmarom assess-skeletal --out results_skeletal
```

It is scientifically possible that the small supplied mechanism admits **no** acceptable skeleton under strict criteria. `assess-skeletal` is designed to exit nonzero in that case; do not override or mislabel it as a success.
