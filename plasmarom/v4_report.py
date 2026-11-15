"""Generate an evidence report from actual computations (no canned metrics)."""
import json
from pathlib import Path
import numpy as np
from .verification import convergence_study
from .plasma_states import run_plasma
from .thermal_flow import run_thermal
from .hybrid_closure import fit_and_holdout
from .benchmark import run_benchmark


def generate_v4_report(out, with_benchmark=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    convergence=convergence_study(out=out)
    plasma=run_plasma(out)
    thermal=run_thermal(out,n=32,dt=.014,duration=1.4)
    hybrid=fit_and_holdout()
    if with_benchmark:
        run_benchmark(out)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases=hybrid['heldout'];x=np.arange(len(cases))
    fig,ax=plt.subplots(figsize=(7,4))
    ax.bar(x-.18,[d['qssa_max_abs_product_error'] for d in cases],width=.34,label='QSSA')
    ax.bar(x+.18,[d['lag_max_abs_product_error'] for d in cases],width=.34,label='Lag fit')
    ax.set(xticks=x,xticklabels=[f"S={d['pulse_strength']}, W={d['pulse_width']}" for d in cases],
           ylabel='Max absolute converted-pool error',title='Held-out toy kinetics; not plasma prediction')
    ax.legend();fig.tight_layout();fig.savefig(out/'heldout_model_comparison.png',dpi=155);plt.close(fig)
    report={'model_status':'SURROGATE NUMERICAL DEMONSTRATOR; NOT A VALIDATED PAC SOLVER',
            'verification':{k:{'observed_order':v['observed_order'],
                                'grid_errors':{str(row['n']):row['l2_error'] for row in v['rows']},
                                'max_mass_drift':max(row['mass_relative_error'] for row in v['rows'])}
                            for k,v in convergence.items()},
            'plasma':plasma,'thermal_surrogate':thermal,'lag_fit':hybrid,
            'cantera_executed_in_this_report':False,
            'not_established':['Experimental accuracy','EEDF / electric field','P-DRGEP','Real reacting CFD','DNS/LES','HPC scaling']}
    (out/'SUMMARY_V4.json').write_text(json.dumps(report,indent=2))
    lines=['# PlasmaROM-H2 v0.4 — computed verification report','',
           '> Automatically built from code execution. A numerical verification is not physical validation.','',
           '## Spatial verification',
           f"- Periodic upwind advection observed order: **{convergence['advection']['observed_order']:.3f}** (target ≈1).",
           f"- Centered diffusion observed order: **{convergence['diffusion']['observed_order']:.3f}** (target ≈2).",
           '','## Three-level H2 vibrational balance',
           f"- Data status: **{plasma['data_status']}**.",
           f"- Energy levels E(v=0,1,2) [eV]: {plasma['levels_eV']}.",
           f"- Population closure: {plasma['max_population_sum_error']:.3e}.",
           f"- Relative electron/vibrational/gas-energy residual: {plasma['relative_energy_closure']:.3e}.",
           '- Assumed Maxwellian EEDF and invented excitation cross sections are **not** calibrated plasma physics.',
           '','## Dimensionless heat source + transport',
           f"- Mean chemical+sensible energy proxy drift: {thermal['max_global_energy_proxy_drift']:.3e}.",
           f"- Maximum local species-sum drift: {thermal['max_total_species_drift']:.3e}.",
           '- Not physical enthalpy or reacting flow with variable density.',
           '','## Fitted QSSA lag correction',
           f"- Fitted dimensionless lag: {hybrid['fitted_tau_dimensionless']:.6g}.",
           '- Held-out cases were not used for parameter tuning:',
           '| Pulse strength | Pulse width | Max product error QSSA | Max product error lag |',
           '|---:|---:|---:|---:|']
    for row in hybrid['heldout']:
        lines.append('| '+ ' | '.join(f"{row[k]:.6g}" for k in (
            'pulse_strength','pulse_width','qssa_max_abs_product_error','lag_max_abs_product_error'))+' |')
    lines.extend(['','## Optional chemistry status',
                  '- The Cantera skeletal mechanism builder is implemented but was **not executed** in this environment.',
                  '- Do not state an actual reduced species count, validated ignition delay, or Cantera speedup.',
                  '','## Limitations',
                  '- No measured electron-impact cross sections supplied.',
                  '- No electric field, E/N, Boltzmann EEDF, ionization or charge conservation.',
                  '- No validated P-DRGEP or state-to-state combustion mechanism.',
                  '- No real flame CFD, 3D DNS/LES, or experimental comparison.',
                  '- Small serial CPU timing is not an HPC scaling study.',
                  '','![Spatial convergence](spatial_convergence.png)',
                  '![Synthetic H2 vibrational balance](plasma_states.png)',
                  '![Dimensionless heat proxy](thermal_surrogate.png)',
                  '![Heldout closure errors](heldout_model_comparison.png)',''])
    (out/'REPORT_V4.md').write_text('\n'.join(lines))
    return report
