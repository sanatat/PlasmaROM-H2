# PlasmaROM-H2 v12.0 — corrective release

Four defects from the v11 review were addressed in actual source code. Historical v11 outputs are retained and should not be mistaken for v12 experimental data.

## 1. NH₃/H₂/CO equivalence ratio

`plasmarom/ammonia_h2.py::check_initial_mixture`: the stoichiometric oxygen requirement is `0.75 NH3 + 0.50 H2 + 0.50 CO`; denominator = supplied O2. Example `NH3=H2=CO=1, O2=1.25` now gives φ=1.40 (formerly 1.00). Unsupported species are rejected instead of silently omitted. The function is a *limited stoichiometric utility*, not a general equilibrium or chemistry validation routine.

## 2. Experimental duplicate measurement guard

`plasmarom/experiment_compare.py::evaluate_experiment`: repeated `(case_id, observable)` entries are rejected with a descriptive error; no silently inflated `n_measurements` or RMSE. For genuinely independent replicate observations, provide separate case IDs and predictions. The software does not certify source independence.

## 3. Superelastic energy channels

`plasmarom/energy_graph.py::energy_aware_screen`: signed `energy_losses_eV` is now supported (negative = superelastic energy gain by electrons). Integrated positive loss and negative gain magnitudes are tracked separately. **Both** nonzero directions must individually retain at least `min_energy_fraction`; a negative contribution cannot cancel a positive contribution to deceive the selection criterion. New fields are `integrated_energy_loss_eV_units`, `integrated_energy_gain_eV_units`, `retained_loss_fraction`, `retained_gain_fraction`, and `net_energy_exchange_eV_units`. The legacy field `retained_energy_fraction` is an absolute-exchange retention fraction. This code is NOT a published/validated P-DRGEP.

## 4. GitHub CI trigger coverage

`.github/workflows/v12-real-chemistry.yml` replaces the narrowly scoped v11 workflow. It triggers on pushes and pull requests changing `plasmarom/**`, `tests/**`, `examples/**`, `protocols/**`, dependency manifests or the workflow. It includes the v12 regression suite. **Local GitHub-hosted execution not verified** and Cantera-dependent steps remain unexecuted in this environment.

## Validation and constraints

See `VERIFICATION_V12.md`. No Cantera, external electron collision data, experimentally validated ammonia plasma, published P-DRGEP, DNS/LES or HPC speed-up is claimed. The optional Cantera reactor example has not been locally executed in v12.
