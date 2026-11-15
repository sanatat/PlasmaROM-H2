"""Energy-channel-aware pruning on a 3-state toy vibrational network.

This is NOT plasma-specific P-DRGEP. It experimentally prunes a transition
channel and measures trajectory errors on held-out power levels. Energy shares
are computed on training FULL trajectories. All cross sections are synthetic.
"""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.constants import e as J_PER_EV
from scipy.integrate import trapezoid
from .electron_energy import (ElectronEnergyConfig,run_coupled_electron_energy,
                              synthetic_excitation_tables)
from .plasma_states import maxwell_rate, levels_eV
from .energy_reduction import screen_power_channels


def _powers_overlap(a,b,*,allow_same_index=False):
    """Recognize duplicates and near-identical cases, including rounded floats."""
    return any(np.isclose(x,y,rtol=1e-10,atol=1e-12)
               for i,x in enumerate(a) for j,y in enumerate(b)
               if not (allow_same_index and a is b and i==j))


def _errors(ref, case):
    return {'relative_max_Te_error':float(np.max(abs(case['te_eV']-ref['te_eV']))/max(1.,float(ref['te_eV'].max()))),
            'max_population_fraction_error':float(np.max(abs(case['fractions']-ref['fractions']))),
            'normalized_gas_heat_error':float(np.max(abs(case['gas_heat_J_m3']-ref['gas_heat_J_m3']))/
                     max(.1,float(np.max(abs(ref['gas_heat_J_m3'])))))}


def evaluate_channel_pruning(out, training_powers=(100.,3000.,6000.),
                              holdout_powers=(500.,4500.),
                              energy_fraction=0.995, allowed_error=0.05):
    training_powers=tuple(training_powers)
    holdout_powers=tuple(holdout_powers)
    if (not training_powers or not holdout_powers or
        any(isinstance(p,bool) or not np.isscalar(p) or
            not np.isfinite(p) or p < 0 for p in (*training_powers,*holdout_powers)) or
        not 0<energy_fraction<=1 or not 0<allowed_error<1):
        raise ValueError('Invalid channel-reduction configuration')
    if _powers_overlap(training_powers,holdout_powers):
        raise ValueError('Data leakage: training and holdout power levels overlap')
    if _powers_overlap(training_powers,training_powers,allow_same_index=True) or _powers_overlap(holdout_powers,holdout_powers,allow_same_index=True):
        raise ValueError('Duplicate cases within a split are not independent')
    tab01,tab12=synthetic_excitation_tables()
    tables=(tab01,tab12)
    delta=np.diff(levels_eV())
    integrated=np.zeros(2)
    training={}
    for p in training_powers:
        cfg=ElectronEnergyConfig(input_power_W_m3=float(p))
        sol=run_coupled_electron_energy(*tables,config=cfg)
        training[p]=sol
        for channel,table in enumerate(tables):
            rates=np.array([maxwell_rate(table,te,delta[channel]) for te in sol['te_eV']])
            power_channel=(cfg.ne_m3*cfg.h2_m3*J_PER_EV*delta[channel]*
                           rates*sol['fractions'][channel])
            integrated[channel]+=float(trapezoid(power_channel,sol['time_s']))
    # Energy selection is based on INTEGRATED GROSS excitation power, not net
    # exchange. Reaction reverse/quenching are retained only for retained links.
    order=np.argsort(-integrated,kind='stable')
    running=0.;keep=[]
    for i in order:
        keep.append(int(i));running+=integrated[i]
        if running >= energy_fraction*integrated.sum():break
    # Keep the v0<->v1 pathway for basic connectivity, regardless of score.
    if 0 not in keep:keep.append(0)
    keep=sorted(set(keep))
    altered=[a.copy() for a in tables]
    for i in range(2):
        if i not in keep:
            altered[i][:,1]=0

    def evaluate_cases(reference_cases):
        metrics=[]
        for power,ref in reference_cases.items():
            approx=run_coupled_electron_energy(*altered, config=ElectronEnergyConfig(input_power_W_m3=float(power)))
            metrics.append({'power_W_m3':float(power),**_errors(ref,approx)})
        return metrics

    train_errors=evaluate_cases(training)
    other={p:run_coupled_electron_energy(*tables,config=ElectronEnergyConfig(input_power_W_m3=float(p)))
           for p in holdout_powers}
    heldout=evaluate_cases(other)
    train_worst=max(max(r[k] for k in r if k!='power_W_m3') for r in train_errors)
    holdout_worst=max(max(r[k] for k in r if k!='power_W_m3') for r in heldout)
    good=len(keep)<2 and train_worst<=allowed_error and holdout_worst<=allowed_error
    result={'status':'HELDOUT_PASS' if good else 'NOT_VALIDATED_REDUCTION',
            'scientific_status':'synthetic 3-state energy-channel pruning; not P-DRGEP',
            'retained_transitions':[f'v{i}<->v{i+1}' for i in keep],
            'removed_transitions':[f'v{i}<->v{i+1}' for i in range(2) if i not in keep],
            'channel_gross_excitation_energy_J_m3':dict(zip(('v0_to_v1','v1_to_v2'),integrated.tolist())),
            'retained_energy_fraction':float(sum(integrated[i] for i in keep)/max(sum(integrated),1e-30)),
            'energy_selection_fraction_target':energy_fraction,
            'allowed_max_normalized_error':allowed_error,
            'training_errors':train_errors,'holdout_errors':heldout,
            'worst_training_error':train_worst,'worst_holdout_error':holdout_worst,
            'training_powers_W_m3':list(training_powers),'holdout_powers_W_m3':list(holdout_powers),
            'warning':'No actual PAC or published P-DRGEP; only synthetic collision tables.'}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'channel_reduction_summary.json').write_text(json.dumps(result,indent=2))
    return result
