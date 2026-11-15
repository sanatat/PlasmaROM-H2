"""Optional *real* NH3/H2 gas-phase reactor demonstration, no plasma.

Uses the Alzueta 2023 modified demonstration mechanism distributed in Cantera
example_data. Any result is a numerical prediction, NOT experimental validation.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

MECHANISM = 'example_data/ammonia-CO-H2-Alzueta-2023.yaml'
MECHANISM_SOURCE = 'https://github.com/Cantera/cantera-example-data/blob/main/ammonia-CO-H2-Alzueta-2023.yaml'

def check_initial_mixture(mix):
    if not isinstance(mix,dict) or not {'NH3','H2','O2','N2'} <= set(mix):
        raise ValueError('Mixture requires NH3, H2, O2, N2')
    if any(isinstance(v, bool) or not isinstance(v,(float,int)) or not np.isfinite(v) or v < 0 for v in mix.values()):
        raise ValueError('All gas amounts must be finite nonnegative')
    # This simplified equivalence ratio is valid only for these fuels and
    # explicitly listed non-fuel species; fail closed on unknown fuel species.
    supported = {'NH3', 'H2', 'CO', 'O2', 'N2', 'AR', 'Ar', 'HE', 'He', 'H2O', 'CO2'}
    unknown = set(mix) - supported
    if unknown:
        raise ValueError(f'Cannot compute equivalence ratio for unsupported species: {sorted(unknown)}')
    if mix['NH3']<=0 or mix['H2']<=0 or mix['O2']<=0:
        raise ValueError('Mixture requires positive fuel and oxidizer')
    # NH3 + 3/4 O2 -> 1/2 N2 + 3/2 H2O
    # H2  + 1/2 O2 -> H2O; CO + 1/2 O2 -> CO2
    required = 0.75*mix['NH3'] + 0.5*mix['H2'] + 0.5*mix.get('CO', 0.)
    phi = required/mix['O2']
    return float(phi)

def _ignition(times, temperatures, threshold=20.0):
    from .cantera_demo import _ignition_assessment
    return _ignition_assessment(times, temperatures, rise_K=threshold)

def run_ammonia_h2(out, *, temperatures=(1350.,1600.), pressure_atm=1.,
                   end_time_s=.02, n_points=160,
                   mixture=None, mech=MECHANISM, phase='baseline'):
    if not 0 < end_time_s <= 2 or n_points<20 or not 0<pressure_atm<=100:
        raise ValueError('Invalid numerical settings')
    if not temperatures or any(not np.isfinite(t) or t<300 or t>3500 for t in temperatures):
        raise ValueError('Invalid initial temperatures')
    mixture = mixture or {'NH3':1.,'H2':1.,'O2':1.25,'N2':4.7}
    phi = check_initial_mixture(mixture)
    try:
        import cantera as ct
    except ImportError as exc:
        raise RuntimeError('Cantera missing: pip install cantera==3.2.0') from exc
    gas=ct.Solution(mech,phase)
    for k in mixture:
        if k not in gas.species_names:
            raise ValueError(f'{k} absent from selected reaction mechanism')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=[];cases=[]
    for T0 in temperatures:
        gas.TPX=float(T0),pressure_atm*ct.one_atm,mixture
        initial_mass=gas.Y.copy()
        reactor=ct.IdealGasConstPressureReactor(gas,energy='on',clone=True)
        net=ct.ReactorNet([reactor]);net.rtol=1e-8;net.atol=1e-17
        times=np.linspace(0.,end_time_s,n_points)
        temps=[]
        for t in times:
            net.advance(float(t))
            state=reactor.thermo
            temps.append(float(state.T))
            rows.append({'initial_T_K':T0,'time_s':float(t), 'T_K':float(state.T),
                         **{f'X_{name}':float(state[name].X[0]) for name in ('NH3','H2','O2','H2O','OH','NO','N2O') if name in state.species_names}})
        ignition=_ignition(times,temps)
        cases.append({'initial_T_K':float(T0),'pressure_atm':pressure_atm,
                      'max_temperature_rise_K':float(max(temps)-T0),**ignition})
    with (out/'ammonia_h2_reactor.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    info={'status':'NUMERICAL_CANTERA_RUN_NO_PLASMA_NO_EXPERIMENTAL_VALIDATION',
          'cantera_version':ct.__version__,'mechanism':mech,'mechanism_citation':MECHANISM_SOURCE,
          'phase':phase,'n_species':gas.n_species,'n_reactions':gas.n_reactions,
          'mixture_molar_amounts':mixture,'equivalence_ratio_from_stoichiometry':phi,
          'cases':cases,'note':'Ignition can be right-censored; no observed ignition is fabricated.'}
    (out/'ammonia_h2_summary.json').write_text(json.dumps(info,indent=2))
    return info
