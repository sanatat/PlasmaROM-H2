# PlasmaROM-H2 v11.0 — execution-first extensions

## Added
- `plasmarom/ammonia_h2.py`: numerical constant-pressure NH3/H2 reactor with Alzueta-2023 Cantera example mechanism, explicit right-censoring of ignition.
- `plasmarom/cantera_boltzmann.py`: calls **Cantera 3.2's actual two-term Boltzmann EEDF solver** using Phelps N2/O2 cross sections from the separately installed Cantera example data. This is neither H2/NH3 plasma nor a dynamic discharge model.
- `plasmarom/bolos_bridge.py`: optional external BOLOS BOLSIG+ collision input with provenance and source hash; does not rename native f0 to an energy probability PDF.
- `plasmarom/energy_graph.py`: graph + integrated electron-energy channel screening. **Not published P-DRGEP**.
- `plasmarom/reactive_transport_1d.py`: executed finite-volume 1D surrogate A→B with conservative total mass and global enthalpy, plus an analytic homogeneous-reaction test. **Not DNS/LES or real combustion**.
- `plasmarom/experiment_compare.py`: audits imported independent experimental CSV plus uncertainty, observables, units and provenance; performs numerical comparison only when actual evidence is supplied.
- GitHub Actions workflow that **requires** Cantera import, NH3/H2 simulation and real EEDF execution, with artifacts retained on success/failure.
- `run_verified_pipeline.sh` for internet-connected Linux, optional integration test suite, and `GAP_STATUS_V11.md`.

## Local verification

The environment lacks Cantera and BOLOS. Python package installation failed due to unavailable DNS/network; no real Cantera/Boltzmann/NH3 numerical output is claimed. The scalar-transport solver and synthetic energy graph *were* executed and checked locally.

For current test counts and numeric results see `VERIFICATION_V11.md`. Older test logs reflect earlier releases.
