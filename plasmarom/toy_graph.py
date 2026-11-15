"""Traceable DRGEP demonstration on the synthetic five-pool network.

Artificial excitation energies (1, 3 model units) are used to explain WHY
energy-target-aware mechanism reduction differs from a chemical-only target.
This is NOT the published plasma-specific P-DRGEP algorithm, which requires
calibrated electron-energy transfer and branching into electron processes.
"""
import numpy as np
from scipy.integrate import trapezoid
from .kinetics import FULL_NAMES,pulse,Parameters
from .reduction import multi_state_importance,choose_species

CHANNELS=('F0 -> F1 (source)','F1 -> F0','F1 -> F2 (source)','F2 -> F1',
          'F2 -> R','F0 -> R','R -> P','R -> F0')
NU=np.array([[-1,1,0,0,0,-1,0,1],
             [1,-1,-1,1,0,0,0,0],
             [0,0,1,-1,-1,0,0,0],
             [0,0,0,0,1,1,-1,-1],
             [0,0,0,0,0,0,1,0]],dtype=float)
PARTICIPATION=(np.abs(NU)>0)


def channel_rates(times, populations, p=Parameters()):
    f0,f1,f2,r,_=populations
    drive=pulse(times,p)
    return np.array([p.excitation*drive*f0,
                     p.relaxation10*f1,
                     p.second_excitation*drive*f1,
                     p.relaxation21*f2,
                     p.activated_oxidation*f2,
                     p.background_oxidation*f0,
                     p.radical_oxidation*r,
                     p.radical_quench*r]).T


def graph_diagnostics(sol,p=Parameters(),threshold=0.95, energy_tolerance=0.03):
    """Chemical DRGEP on full trajectories plus illustrative energy channel guard.

    The energy guard retains the two reactant/product species of any plasma
    source channel contributing >= tolerance to an artificial integrated
    power-budget diagnostic. Does not guarantee preserved plasma physics.
    """
    rates=channel_rates(sol.t,sol.y,p)
    chem=multi_state_importance(NU,rates,PARTICIPATION,targets=[0,4])
    selected=choose_species(chem,FULL_NAMES,threshold)
    energy=np.array([trapezoid(rates[:,0],sol.t)*1.,
                     trapezoid(rates[:,2],sol.t)*3.])
    fractions=energy/energy.sum() if energy.sum()>0 else np.zeros_like(energy)
    energy_guard=set()
    for j,index in enumerate((0,2)):
        if fractions[j]>=energy_tolerance:
            energy_guard.update(FULL_NAMES[k] for k in np.flatnonzero(PARTICIPATION[:,index]))
    combined=[n for n in FULL_NAMES if n in selected or n in energy_guard]
    return {'names':FULL_NAMES,'chemical_importance':dict(zip(FULL_NAMES,chem.tolist())),
            'chemical_threshold':threshold,'chemical_selected':selected,
            'energy_tolerance':energy_tolerance,'toy_energy_fractions':
                {'F0->F1':float(fractions[0]),'F1->F2':float(fractions[1])},
            'energy_guard_species':sorted(energy_guard),
            'energy_aware_candidate':combined,
            'warning':'Toy energy weights 1 and 3 are invented model units, not electron cross-section measurements. No validated skeletal mechanism or P-DRGEP is produced.'}
