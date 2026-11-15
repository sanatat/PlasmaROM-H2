# Reproducibility and testing

## Local requirements

Python >=3.10; on Linux/macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
python -m pytest -q tests/test_v12_regressions.py
```

Other regression tests can be run with `python -m pytest -q -m 'not slow'`. The multiprocessing/scaling module is marked `slow` to avoid expensive/occasionally environment-dependent execution in the fast default CI. Run it explicitly with `python -m pytest -q tests/test_v6_scaling.py` only when the environment supports local child processes. This is not HPC validation.

## External solver integrations

```bash
python -m pip install -e '.[test,cantera]'
python -m plasmarom ammonia-h2 --out results_local/ammonia
python -m plasmarom boltzmann --out results_local/boltzmann
```

External Cantera mechanism/example data must be available. Failing to import Cantera, load a mechanism, or satisfy a pinned acceptance tolerance is a recorded **failure or non-acceptance**; never silently convert it into a success. A manual `Scientific numerical-acceptance gate` GitHub workflow is supplied, and a failed gate does not invalidate the separate synthetic unit tests.

## Evidence provenance

Historical `results_v*` folders document snapshots from different versions. They are not interchangeable with v12 output. The observed v12 verification is recorded in [`history/VERIFICATION_V12.md`](history/VERIFICATION_V12.md), including skipped and interrupted tests. `SHA256SUMS.txt` enumerates distributed files, and may be checked using `sha256sum -c SHA256SUMS.txt` from repository root (Linux).

## Data labels

Files whose names begin with `SYNTHETIC_` are not measurements. The `examples/data/*SCHEMA.csv` files are headers/templates, not validated experiments. Use the experiment adapter only with externally verified measurements and appropriate independent predictions.

The GitHub-preparation test run is recorded in [`GITHUB_RELEASE_CHECK.md`](GITHUB_RELEASE_CHECK.md).
