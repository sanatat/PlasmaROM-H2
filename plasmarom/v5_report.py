"""Generate a transparent v5 summary from on-disk machine-readable evidence."""
import json
from pathlib import Path


def generate_v5_report(out):
    p=Path(out)
    required=['electron_energy_summary.json','channel_reduction_summary.json',
              'external_eedf_rate.json','numerical_verification.json','parallel_comparison.json']
    data={name:json.loads((p/name).read_text()) for name in required}
    energy=data['electron_energy_summary.json']
    reduction=data['channel_reduction_summary.json']
    ext=data['external_eedf_rate.json']
    verification=data['numerical_verification.json']
    parallel=data['parallel_comparison.json']
    message=f'''# PlasmaROM-H2 v0.5 — computed verification report

**Status:** educational, independently executable numerical portfolio prototype. NOT validated plasma-assisted combustion.

## Coupled electron-energy closure (synthetic cross sections)

- Initial / final electron Te: {energy['initial_Te_eV']:.6g} / {energy['final_Te_eV']:.6g} eV
- Deposited power energy: {energy['total_power_deposited_J_m3']:.6g} J/m³
- Max relative total-energy residual: {energy['relative_energy_error']:.4g}
- Max population conservation residual: {energy['max_population_sum_error']:.4g}
- Electron Te is solved as an ODE with constant externally prescribed input power and Maxwellian approximation; fixed ne and gas bath. **No field, E/N, Boltzmann solver, or charge balance.**

## Energy-aware reaction channel pruning (synthetic)

- Result: **{reduction['status']}**
- Retained: {', '.join(reduction['retained_transitions'])}
- Removed: {', '.join(reduction['removed_transitions']) or 'none'}
- Gross excitation energy retained: {reduction['retained_energy_fraction']:.6%}
- Maximum held-out normalised trajectory error: {reduction['worst_holdout_error']:.5g} (limit {reduction['allowed_max_normalized_error']:.5g})
- **This is not published plasma P-DRGEP.** There are only two transitions and no real collision dataset.

## External EEDF CSV interface

- Test data: explicitly **SYNTHETIC** Maxwellian EEDF and invented excitation cross section
- Computed cross-section rate: {ext['rate_coefficient_m3_s']:.7g} m³/s
- EEDF normalisation integral: {ext['eedf_integral']:.8g}
- CSV SHA256 hashes and provenance are saved in `external_eedf_rate.json`. External energy distribution must already be computed; the software **does not solve EEDF**.

## Periodic numerical verification

- Upwind advection observed order: {verification['advection']['observed_order']:.4f} (≈ first order)
- Centered diffusion observed order: {verification['diffusion']['observed_order']:.4f} (≈ second order)
- These are manufactured-PDE correctness tests, not flame validation.

## Local serial versus two-process check

- Serial walltime: {parallel['serial_wall_s']:.4g} s
- Parallel walltime: {parallel['parallel_wall_s']:.4g} s
- Observed speedup: {parallel['measured_speedup']:.4g}× (may be below one owing to overhead)
- Max difference between results: {parallel['max_absolute_serial_parallel_difference']:.3g}
- One uncalibrated run of independent tasks: **NOT an HPC strong/weak scaling experiment**.

## Cannot be marked verified in this environment

- `Cantera >=3.2` was unavailable; its reactor/skeletal mechanism and 1D H₂/air flame integrations were **not executed**.
- Without a licensed real electron-impact dataset and externally computed EEDF, physical plasma coefficients remain **unverified**.
- No plasma-specific P-DRGEP, full electron Boltzmann EEDF solver, 3D DNS/LES, calibrated reacting-flow CFD, or experimental PAC validation.
- Do not invent a real reduced mechanism, ignition delay, flame speed, measured error or speedup.
'''
    (p/'REPORT_V5.md').write_text(message)
    return p/'REPORT_V5.md'
