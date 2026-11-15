# Verification and evidence scope — v9

## Local test matrix

| Check | Result | What it does **not** prove |
|---|---|---|
| Numerical + audit suite | 131 passed, 0 failed | Independent experimental correctness |
| Optional Cantera tests | 5 skipped (not installed) | Executable chemistry or valid reduced mechanism |
| Composition-separated tests | Passed; normalized mixture signatures reject identical compositions | Generalization on real ignition data |
| Tampered Training records | Rejected in mocked Cantera-control regression | Actual Cantera training replay in this offline environment |
| Modified thermo/NASA metadata | Rejected in mocked mechanism-metadata regression | Actual Cantera YAML round-trip in this offline environment |
| Independent 0D training and holdout replay | Implemented, fail-closed on exception | Experimental validation or mechanism-reduction success |
| Local task parallelism | Two scaling tests passed | HPC solver speedup |

The full test summary has been obtained from the verbose pytest log and the separate numerical run. Locally **no** optional Cantera test executed. `REPRODUCE_V9.json` has only manufactured numerical verification and scenario parsing.

## How to run physical-chemistry integration in Linux Mint

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,cantera]'
python -m pytest -q tests/test_cantera_optional.py -rs
python -m plasmarom skeletal --mechanism h2o2.yaml \
    --cases-file examples/cases/h2_composition_sweep.json \
    --reactor-t-end 0.006 --snapshots 80 --out results_skeletal_v9
python -m plasmarom assess-skeletal --out results_skeletal_v9
```

`assess-skeletal` is **expected to exit nonzero** if no mechanism passes all independent checks; this is a valid negative result, not a reason to label a model validated. The physical initial time window and ignition thresholds require independent suitability checks before scientific interpretation. Longer or more appropriate time horizons may be necessary.

## Known limitations

The focus is evidence integrity: unit-tested source/target identity, explicit case separation, and independent replay. State-resolved H2 plasma remains synthetic; no plasma P-DRGEP, real cross-section data, benchmarked NH3/H2 flames, validated reactive CFD, or HPC strong/weak scaling is supplied.
