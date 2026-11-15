# PlasmaROM-H2 v11 — independent verification and scope

Date: 2026-11-15

## Verification actually run in this environment

- Core pytest suite excluding potentially blocking multiprocessing module: **149 passed, 6 skipped**, no failures; see `results_v11/PYTEST_CORE_V11.log`. Skips include Cantera-dependent tests; the current environment has **no Cantera**.
- Multiprocessing module: **2 passed earlier in a separate invocation during this v11 build**, 7.60 s. A later repetition of this module timed out in the tool environment, so its repeat was **not** reported as another pass.
- New v11 tests: **9 passed** before and after code changes, included in core suite above.
- Executed 1D transport: `results_v11/transport/reactive_transport_1d.json`, `.npz`, `.png`. Global mass drift 1.1102230246251565e-16; global total enthalpy drift 0.0 J/kg; max local A+B-1 2.220446049250313e-16. This is not CFD or DNS/LES.
- Executed graph/energy channel screening: `results_v11/SYNTHETIC_energy_graph.json`, explicitly synthetic. It is **not** published plasma P-DRGEP.

## Not executed, not validated

- `cantera==3.2.0` installation failed: DNS/index unavailable. All Cantera real chemistry and plasma-EEDF examples are **implemented but unexecuted** here.
- `bolos` is not installed. Its bridge is **unexecuted**.
- No real NH3/H2 reactor result, plasma-ignition experiment, verified mechanism reduction, physical 3D DNS/LES or HPC-speedup claim.
- No independent experimental measurement file was available to populate the validation input; comparison code was tested on invalid-data cases only.

## Reproduction on internet-connected Linux

```bash
bash run_verified_pipeline.sh
```

If the pipeline passes, this proves certain numerical integrations run, not that plasma-assisted combustion is experimentally validated. CI logs and raw outputs are the relevant evidence.

## Research sources and attribution

- Cantera real two-term Boltzmann EEDF: https://cantera.org/stable/examples/python/thermo/plasma-eedf.html
- Phelps cross-section example: https://github.com/Cantera/cantera-example-data/blob/main/air-plasma-Phelps.yaml
- Alzueta-2023 NH3/H2: https://github.com/Cantera/cantera-example-data/blob/main/ammonia-CO-H2-Alzueta-2023.yaml
- Reference electron collision input and BOLOS: https://github.com/aluque/bolos

**The scientific gaps are not all complete.** `GAP_STATUS_V11.md` has the full ledger.
