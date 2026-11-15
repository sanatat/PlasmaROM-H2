"""Real 2-term Boltzmann EEDF integration in Cantera >=3.2.

Cross-sections are provided by the external Cantera-distributed Phelps
example_data/air-plasma-Phelps.yaml. Not bundled here. Demonstrates N2/O2,
NOT H2/NH3 or a self-consistent discharge with dynamic E/N and charges.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid

SOURCE = ('https://cantera.org/stable/examples/python/thermo/plasma-eedf.html')
CROSS_SECTIONS = ('https://github.com/Cantera/cantera-example-data/blob/main/air-plasma-Phelps.yaml')

def run_phelps_eedf(out, *, fields_Td=(50.,100.,200.), temperature_K=300.):
    if not fields_Td or not 100<=temperature_K<=1500 or any(not np.isfinite(v) or v<=0 for v in fields_Td):
        raise ValueError('Invalid E/N or temperature')
    try:
        import cantera as ct
    except ImportError as exc:
        raise RuntimeError('Cantera 3.2+ missing (needed for two-term EEDF)') from exc
    if tuple(int(x) for x in ct.__version__.split('.')[:2])<(3,2):
        raise RuntimeError('Cantera >=3.2 required for Boltzmann-two-term')
    gas=ct.Solution('example_data/air-plasma-Phelps.yaml')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=[];results=[]
    for en in fields_Td:
        # Published demonstration composition; tiny charge seed concentrations.
        gas.TPX=float(temperature_K),101325.,'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'
        gas.reduced_electric_field=float(en)*1e-21
        gas.update_electron_energy_distribution() # Cantera's REAL two-term solver
        energy=np.asarray(gas.electron_energy_levels,float)
        f0=np.asarray(gas.electron_energy_distribution,float)
        if energy.ndim!=1 or energy.shape!=f0.shape or (not np.isfinite(f0).all()) or np.any(f0<0):
            raise RuntimeError('EEDF shape, finiteness or positivity invalid')
        # Cantera f0 is normalized with sqrt(epsilon) weight, not a PDF directly.
        norm=float(trapezoid(np.sqrt(energy)*f0,energy))
        mean=float(trapezoid(energy**1.5*f0,energy)/norm)
        if not 0.95<=norm<=1.05 or mean<=0:
            raise RuntimeError(f'EEDF normalization/mean-energy check failed: {norm=}, {mean=}')
        results.append({'field_Td':float(en),'normalisation_integral':norm,'mean_energy_eV':mean})
        rows.extend({'field_Td':float(en),'energy_eV':float(e),'boltzmann_f0':float(y)} for e,y in zip(energy,f0))
    with (out/'boltzmann_eedf.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['field_Td','energy_eV','boltzmann_f0']);w.writeheader();w.writerows(rows)
    report={'status':'REAL_CANTERA_TWO_TERM_BOLTZMANN_EEDF_N2_O2_NUMERICAL_ONLY',
            'cantera_version':ct.__version__,'collision_set':CROSS_SECTIONS,
            'citation':SOURCE,'comparison_TO_experiment':'NOT_PERFORMED',
            'no_NH3_H2_claim':True,'no_discharge_model':True,'results':results}
    (out/'boltzmann_summary.json').write_text(json.dumps(report,indent=2))
    return report
