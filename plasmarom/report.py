"""Generate an evidence-based run report entirely from the JSON artifacts."""
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy


def write_report(out, root=None):
    out=Path(out)
    kinetic=json.loads((out/'kinetics_results.json').read_text())
    spatial=json.loads((out/'flow_results.json').read_text())
    full=spatial[0];reduced=spatial[1];strong=spatial[-1]
    timing=kinetic['runtimes']
    ratio=timing['full_seconds_median']/timing['QSSA_seconds_median']
    worst=max(case['max_abs_P'] for case in kinetic['validation_grid'])
    graph=kinetic['toy_drgep_and_energy_target_screening']
    ml=kinetic['learned_nonnegative_closure']
    text=f'''# PlasmaROM-H2 v0.3 — computed demonstration report

Generated **{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}** from the actual output files, not manually transcribed results.

**Status:** verified surrogate numerical prototype, **not** calibrated plasma chemistry or combustion DNS/LES.

## Environment
- Python: {platform.python_version()}
- NumPy: {np.__version__}; SciPy: {scipy.__version__}
- Array sizes for displayed run: {full['n']}² (fine) and {spatial[2]['n']}² (coarse)

## 0D surrogate chemical kinetics
- State count: full **5**, QSSA **4** (20% fewer toy variables). Reduced B means F1+F2, never F1 alone.
- QSSA rates now use F1=B/(1+r) and F2=rB/(1+r), where r=k12/(k21+k2r). Full/reduced comparisons use corresponding [F0,B,R,P] observables.
- Nominal maximum absolute converted-pool error: **{kinetic['max_abs_error_F0_B_R_P'][3]:.7g}**.
- Maximum converted-pool error over nine original condition combinations: **{worst:.7g}**.
- Maximum omitted F2 fraction in the reference full model: **{kinetic['max_omitted_F2_population_full']:.7g}**; the reduced model's four variables are effective/lumped pools, not identical full physical fractions.
- Gross full-balance residual: **{kinetic['mass_residual_full']:.3e}**; reduced **{kinetic['mass_residual_QSSA']:.3e}**.
- Median execution times: full **{timing['full_seconds_median']:.4f} s**, QSSA **{timing['QSSA_seconds_median']:.4f} s** (ratio **{ratio:.2f}×** on this machine, not a general speedup claim).
- PCA is only an explained-variance diagnostic; NOT a reduced kinetic integrator.

## Chemical DRGEP and energy-target screening
- Correct chemical DRGEP direct-interaction coefficients and maximum-product path propagation are implemented and unit-tested.
- Toy network chemical importance: {json.dumps(graph['chemical_importance'])}
- Chemical-only candidate: {', '.join(graph['chemical_selected'])}.
- Artificial source-channel energy fractions (weights **1 and 3 toy units**): {json.dumps(graph['toy_energy_fractions'])}.
- Energy-guarded candidate: {', '.join(graph['energy_aware_candidate'])}.
- **These illustrative energy targets are not the Bellemans P-DRGEP algorithm.** No electron-impact cross sections, EEDF, E/N, physical energy branching, or validation of a skeletal plasma mechanism is included.

## Data-fitted, positivity-preserving reduced closure
- Fitted two **nonnegative** coefficients using SciPy NNLS on six reference parameter combinations.
- Weights: {ml['weights']}; training residual L2: {ml['training_residual_l2']:.5g}.
- Three HELD-OUT test conditions (not used for fitting): {json.dumps(ml['heldout_cases'])}.
- The constrained transfer flux preserves the toy material sum. The learned closure may be WORSE than analytic QSSA outside the training regime; data-driven fitting is **not** a guarantee of generalization.
- This is not a PINN, physics-informed deep learning, or a plasma kinetics prediction.

## 2D velocity and reactive scalar transport
- Scalar transport uses conservative face fluxes from a discrete divergence-free streamfunction; each scalar has a global conservation test under pure advection/diffusion.
- Peak spectral velocity divergence: **{full['peak_divergence']:.3e}**.
- Taylor–Green velocity L2 discrepancy: **{full['taylor_green_l2']:.3e}**.
- Minimum surrogate pool: **{full['min_pool']:.3e}**; maximum local sum-of-pools drift **{full['sum_deviation']:.3e}**.
- Converted-pool MAE of same-grid QSSA against full: **{reduced['converted_pool_MAE_vs_full_fine']:.5g}**.
- Coarse-grid converted-pool MAE vs sampled reference without illustrative eddy closure: **{spatial[2]['converted_pool_MAE_vs_full_fine_downsampled']:.5g}**, with closure: **{spatial[3]['converted_pool_MAE_vs_full_fine_downsampled']:.5g}**. These values DO NOT validate LES or DNS.
- Temporal halving on 16²: product MAE **{spatial[-2]['product_MAE']:.3e}**. Single pair is sensitivity only; it is not a measured convergence order.
- **Stress test** pulse strength {strong['pulse_strength']:.0f}, n={strong['n']}, dt={strong['dt']}: min pool **{strong['min_pool']:.3e}**, conservation deviation **{strong['sum_deviation']:.3e}**, max chemical subcycles per half-step **{strong['max_substeps']}**.

## Optional physically dimensioned chemistry
`python -m plasmarom cantera --out results/cantera_case` invokes Cantera's bundled `h2o2.yaml` and generates a constant-pressure adiabatic H2/O2/N2 reactor data table, H/O atom-balance diagnostics, enthalpy-balance diagnostics, and chemical DRGEP species ranking. Cantera installation was unavailable in this environment: this path was **not run or verified locally**. Interface signatures were checked against the Cantera 3.2 documentation, which is not a substitute for executing the optional test. It is a combustion CHEMISTRY demo, not PAC; the selected species are merely CANDIDATES, not a validated reduced mechanism. Requires Cantera >=3.2.

## Scope and limitations
The 2D flow is incompressible and is not a flame; temperature, density, pressure-wave dynamics and coupled thermal-energy equation are absent from the CFD portion. The synthetic state pools are not literal energy levels derived from spectroscopic measurements. The LES-style stress experiment does not represent turbulent DNS/LES validation. The project has no experimental plasma benchmark or scientifically valid ignition-delay prediction.

## Figures (relative to this report)
![Toy kinetics](kinetics.png)
![Channel fluxes](pathway_flux.png)
![2D scalar fields](reactive_flow.png)
![PCA](pca.png)

## Verification
Run `python -m pytest -q` and `python -m plasmarom run --quick` before using this as evidence in an application. Do not state that this report establishes scientific validation of plasma-assisted combustion.
'''
    target=out/'TECHNICAL_REPORT.md'
    target.write_text(text,encoding='utf-8')
    return target
