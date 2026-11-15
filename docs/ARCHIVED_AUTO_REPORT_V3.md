# PlasmaROM-H2 v0.3 — computed demonstration report

Generated **2026-11-14 09:52 UTC** from the actual output files, not manually transcribed results.

**Status:** verified surrogate numerical prototype, **not** calibrated plasma chemistry or combustion DNS/LES.

## Environment
- Python: 3.13.5
- NumPy: 2.3.5; SciPy: 1.17.0
- Array sizes for displayed run: 32² (fine) and 16² (coarse)

## 0D surrogate chemical kinetics
- State count: full **5**, QSSA **4** (20% fewer toy variables). Reduced B means F1+F2, never F1 alone.
- QSSA rates now use F1=B/(1+r) and F2=rB/(1+r), where r=k12/(k21+k2r). Full/reduced comparisons use corresponding [F0,B,R,P] observables.
- Nominal maximum absolute converted-pool error: **0.002718223**.
- Maximum converted-pool error over nine original condition combinations: **0.006618566**.
- Maximum omitted F2 fraction in the reference full model: **0.0540825**; the reduced model's four variables are effective/lumped pools, not identical full physical fractions.
- Gross full-balance residual: **2.414e-12**; reduced **1.669e-11**.
- Median execution times: full **0.0694 s**, QSSA **0.0548 s** (ratio **1.27×** on this machine, not a general speedup claim).
- PCA is only an explained-variance diagnostic; NOT a reduced kinetic integrator.

## Chemical DRGEP and energy-target screening
- Correct chemical DRGEP direct-interaction coefficients and maximum-product path propagation are implemented and unit-tested.
- Toy network chemical importance: {"H2(v0) proxy": 1.0, "H2(v1) proxy": 0.9919892755784794, "H2(v2) proxy": 0.8918527186163895, "Reactive pool": 1.0, "Converted pool": 1.0}
- Chemical-only candidate: H2(v0) proxy, H2(v1) proxy, Reactive pool, Converted pool.
- Artificial source-channel energy fractions (weights **1 and 3 toy units**): {"F0->F1": 0.4211590714917215, "F1->F2": 0.5788409285082785}.
- Energy-guarded candidate: H2(v0) proxy, H2(v1) proxy, H2(v2) proxy, Reactive pool, Converted pool.
- **These illustrative energy targets are not the Bellemans P-DRGEP algorithm.** No electron-impact cross sections, EEDF, E/N, physical energy branching, or validation of a skeletal plasma mechanism is included.

## Data-fitted, positivity-preserving reduced closure
- Fitted two **nonnegative** coefficients using SciPy NNLS on six reference parameter combinations.
- Weights: [0.24086337076530043, 0.0]; training residual L2: 0.45486.
- Three HELD-OUT test conditions (not used for fitting): [{"pulse_strength": 2.5, "pulse_width": 0.12, "learned_max_abs_product_error": 0.003173047082287557, "qssa_max_abs_product_error": 0.004291513270753833, "learned_total_pool_residual": 5.8966165283891314e-12, "learned_min_pool": 0.0}, {"pulse_strength": 6.0, "pulse_width": 0.34, "learned_max_abs_product_error": 0.086379405465144, "qssa_max_abs_product_error": 0.014409603762347345, "learned_total_pool_residual": 4.169997680492088e-12, "learned_min_pool": 0.0}, {"pulse_strength": 0.8, "pulse_width": 0.42, "learned_max_abs_product_error": 0.007892791194471338, "qssa_max_abs_product_error": 0.0024616604323209335, "learned_total_pool_residual": 5.4569682106375694e-12, "learned_min_pool": 0.0}].
- The constrained transfer flux preserves the toy material sum. The learned closure may be WORSE than analytic QSSA outside the training regime; data-driven fitting is **not** a guarantee of generalization.
- This is not a PINN, physics-informed deep learning, or a plasma kinetics prediction.

## 2D velocity and reactive scalar transport
- Scalar transport uses conservative face fluxes from a discrete divergence-free streamfunction; each scalar has a global conservation test under pure advection/diffusion.
- Peak spectral velocity divergence: **5.218e-15**.
- Taylor–Green velocity L2 discrepancy: **2.276e-09**.
- Minimum surrogate pool: **1.023e-24**; maximum local sum-of-pools drift **8.882e-16**.
- Converted-pool MAE of same-grid QSSA against full: **4.3732e-05**.
- Coarse-grid converted-pool MAE vs sampled reference without illustrative eddy closure: **0.00066261**, with closure: **0.00067223**. These values DO NOT validate LES or DNS.
- Temporal halving on 16²: product MAE **3.270e-07**. Single pair is sensitivity only; it is not a measured convergence order.
- **Stress test** pulse strength 100, n=20, dt=0.02: min pool **5.913e-21**, conservation deviation **2.665e-15**, max chemical subcycles per half-step **7**.

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
