# PlasmaROM-H2 v0.3: independent regression checks

The values below come from executed version 0.3 code, not from a plasma experiment or real combustion validation.

## Conservative finite-volume scalar transport

Two separate positive scalar fields were advanced for 100 steps with `dt=0.005`, zero reactions and zero scalar diffusion. The velocity was divergence-free Taylor–Green on a periodic domain. Each species should conserve its **global** integral.

| Grid | Max relative species integral drift | Min scalar |
|---|---:|---:|
| 16x16 | 0.000e+00 | 3.696e-14 |
| 32x32 | 0.000e+00 | 1.030e-17 |
| 64x64 | 2.180e-16 | 8.162e-22 |

## QSSA consistency and stress errors

Reduced state is B=F1+F2; all B outflows use the algebraic fast-state fractions. Reference full states are mapped to F0/B/R/P for comparisons. Results on pulse_width=0.12:

| Pulse strength | Max absolute product error | Max absolute B error |
|---:|---:|---:|
| 1 | 0.000649452 | 0.00142337 |
| 6 | 0.0196863 | 0.0199189 |
| 15 | 0.0390022 | 0.0298535 |
| 50 | 0.0369055 | 0.0302408 |
| 100 | 0.0313673 | 0.0292703 |

**Numerical accuracy is not physical validation.** The fastest transient still violates a globally reliable QSSA assumption. This test only demonstrates model consistency and illustrates where full and reduced trajectories differ.

## Verification

- Command: `python -m pytest -q` (26 passed, one optional Cantera test skipped in the local environment).
- Command: `python -m plasmarom run --out results_highres` (generated high-resolution CSV/JSON figures and report).
- Cantera is not available here, so neither its ignition nor its conservation results were independently verified. API usage cross-checked against Cantera 3.2 docs. The GitHub optional CI job will run it when the package can be installed.
- Scope excludes state-to-state molecular plasma kinetics, a real electron-Boltzmann equation, full P-DRGEP, a validated skeletal mechanism, and reacting DNS/LES.
