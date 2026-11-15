"""Recompute the version 6 numerical demo from SYNTHETIC test inputs.

Run as: python examples/reproduce_v6.py
The __main__ guard is required by multiprocessing start method 'spawn'.
"""
from pathlib import Path
import sys
import csv
import json
import numpy as np
from scipy.integrate import trapezoid

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT))
from plasmarom.eedf import maxwellian_pdf
from plasmarom.eedf_field import run_field_rate
from plasmarom.hybrid_closure import fit_and_holdout
from plasmarom.plasma_states import levels_eV
from plasmarom.parallel import fixed_workload_scaling


def main():
    data=PROJECT/'examples'/'data'
    out=PROJECT/'results_v6'
    out.mkdir(exist_ok=True)
    E=np.linspace(0,80,1601)
    for field,te in [(50,1.),(200,3.)]:
        pdf=maxwellian_pdf(E,te)
        pdf=pdf/trapezoid(pdf,E)
        with (data/f'SYNTHETIC_EEDF_field_{field}Td.csv').open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['energy_eV','probability_per_eV'])
            w.writerows(zip(E,pdf))
    manifest=data/'SYNTHETIC_EoverN_manifest.csv'
    manifest.write_text('EoverN_Td,eedf_csv,source\n'
       '50,SYNTHETIC_EEDF_field_50Td.csv,SYNTHETIC Maxwellian Te=1eV arbitrarily mapped to 50Td\n'
       '200,SYNTHETIC_EEDF_field_200Td.csv,SYNTHETIC Maxwellian Te=3eV arbitrarily mapped to 200Td\n')
    result=run_field_rate(out,manifest,data/'SYNTHETIC_H2_v0_v1_cross_section.csv',
          'SYNTHETIC INVENTED CROSS SECTION - NOT MEASURED',float(levels_eV()[1]),100.)
    ml=fit_and_holdout()
    (out/'ml_validation_fallback.json').write_text(json.dumps(ml,indent=2))
    scaled=fixed_workload_scaling(out,powers=(100.,500.,1000.,3000.),workers=(1,2),repetitions=2)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fields=[50.,75.,100.,125.,150.,175.,200.]
    rate_points=[run_field_rate(out,manifest,data/'SYNTHETIC_H2_v0_v1_cross_section.csv',
          'SYNTHETIC INVENTED CROSS SECTION - NOT MEASURED',float(levels_eV()[1]),f)
          ['rate_coefficient_m3_s'] for f in fields]
    # Recreate the nominated 100 Td evidence JSON after exploratory grid.
    result=run_field_rate(out,manifest,data/'SYNTHETIC_H2_v0_v1_cross_section.csv',
          'SYNTHETIC INVENTED CROSS SECTION - NOT MEASURED',float(levels_eV()[1]),100.)
    fig,ax=plt.subplots(figsize=(6.3,3.8))
    ax.plot(fields,rate_points,marker='o')
    ax.set(xlabel='Invented E/N mapping [Td] (NOT physical)',
           ylabel='Collision-rate coefficient [m3/s]',
           title='Synthetic field-dependent EEDF input exercise')
    fig.tight_layout();fig.savefig(out/'synthetic_field_rates.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,3.8))
    cases=ml['heldout']
    xpos=np.arange(len(cases))
    ax.bar(xpos-.2,[x['qssa_max_abs_product_error'] for x in cases],width=.4,label='QSSA')
    ax.bar(xpos+.2,[x['lag_max_abs_product_error'] for x in cases],width=.4,label='Fitted lag')
    ax.set(xticks=xpos,xticklabels=[str(i+1) for i in range(len(cases))],
           xlabel='Held-out synthetic condition',ylabel='Maximum absolute product error',
           title='Validation selects model; held-out data only evaluates')
    ax.legend();fig.tight_layout();fig.savefig(out/'ml_heldout.png',dpi=170);plt.close(fig)
    print('Synthetic E/N rate [m3/s]:',result['rate_coefficient_m3_s'])
    print('Validation-selected toy correction:',ml['selected_model_by_validation'])
    print('Local task-parallel speedups (not HPC):',
          {k:round(v['speedup_vs_serial'],3) for k,v in scaled['by_workers'].items()})


if __name__=='__main__':
    main()
