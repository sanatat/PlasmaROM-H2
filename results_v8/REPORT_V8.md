# Numerical results — PlasmaROM-H2 v8

**Generated with synthetic inputs, 2026-11-15. Not experimental plasma or combustion validation.**

| Numerical check | Measured result |
|---|---:|
| Maximum relative detailed-balance ratio error (Te 1,2,5 eV) | 2.220e-16 |
| Electron energy balance relative residual (input 0 W/m³) | 7.403e-09 |
| Electron energy balance relative residual (input 3000 W/m³) | 7.089e-09 |
| Numerical peak dT/dt time (toy trajectory) | 2 s |
| Interpolated 20 K threshold crossing (same trajectory) | 3.857142857 s |

**Tests:** 112 passed, 3 skipped (Cantera not installed), 0 failures across split executions. Full scientific acceptance requires loading and replaying real mechanisms in Cantera. It has not run here.

Machine-readable source: [`VERIFIED_NUMERICAL_RESULTS.json`](VERIFIED_NUMERICAL_RESULTS.json).
