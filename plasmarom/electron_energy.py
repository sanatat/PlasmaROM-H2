"""Zero-dimensional, power-driven Maxwellian electron-energy closure.

A deliberately constrained *plasma surrogate*: fixed ne, fixed nH2 and gas
bath temperature; no E/N, Boltzmann equation, ionization, self-consistent
charge or plasma discharge. Input cross sections are user supplied or explicitly
synthetic. Electron temperature is derived from Ue=3/2 ne e Te, NOT prescribed.
"""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.constants import e as J_PER_EV, k as K_B
from scipy.integrate import solve_ivp
from .plasma_states import levels_eV, maxwell_reversible_rates


@dataclass(frozen=True)
class ElectronEnergyConfig:
    ne_m3: float = 1e17
    h2_m3: float = 1e21
    te_initial_eV: float = 2.0
    gas_temperature_K: float = 300.0
    input_power_W_m3: float = 3000.0
    electron_gas_relaxation_s: float = 1e-4
    quench10_s: float = 6e3
    quench21_s: float = 6e3
    duration_s: float = 2e-4
    n_outputs: int = 121



def synthetic_excitation_tables():
    """Synthetic, non-measured cross sections for code-verification only."""
    gap=np.diff(levels_eV())
    energy=np.linspace(0,100,5001)
    out=[]
    for d in gap:
        sigma=np.zeros_like(energy)
        active=energy>d
        x=energy[active]-d
        sigma[active]=1e-20*(1-np.exp(-x/1.5))*np.exp(-x/35)
        out.append(np.column_stack((energy,sigma)))
    return tuple(out)

def _validate_config(c):
    vals = tuple(vars(c).values())
    if not all(np.isfinite(x) for x in vals):
        raise ValueError('All plasma parameters must be finite')
    if min(c.ne_m3, c.h2_m3, c.te_initial_eV, c.gas_temperature_K,
           c.electron_gas_relaxation_s, c.duration_s) <= 0:
        raise ValueError('Densities, temperatures, relaxation and duration must be positive')
    if min(c.quench10_s, c.quench21_s, c.input_power_W_m3) < 0 or c.n_outputs < 3:
        raise ValueError('Negative rates/power or too few output points')


def run_coupled_electron_energy(table01=None, table12=None, config=ElectronEnergyConfig()):
    """Solve populations + ln(Te) + cumulative bath heat and energy input.

    The Maxwellian assumption is imposed; inelastic collisions use tabulated
    excitation cross sections; detailed balance sets reverse rates. Electron
    energy transfers from bath and vibrations are accounted for in J/m3.
    All input energy is explicit. Does NOT claim to model a discharge.
    """
    _validate_config(config)
    levels = levels_eV()
    gaps = np.diff(levels)
    demo = table01 is None and table12 is None
    if (table01 is None) != (table12 is None):
        raise ValueError('Supply both cross-section tables or neither')
    if demo:
        table01,table12 = synthetic_excitation_tables()
    tables = (np.asarray(table01, float), np.asarray(table12, float))
    ne = config.ne_m3
    nh2 = config.h2_m3
    electron_heat_capacity_J_m3_eV = 1.5*ne*J_PER_EV
    gas_Te_eV = K_B*config.gas_temperature_K/J_PER_EV
    u_initial = electron_heat_capacity_J_m3_eV*config.te_initial_eV
    energy_levels = levels*nh2*J_PER_EV

    def rhs(t, y):
        n0, n1, n2, log_te, heat, deposited = y
        te = float(np.exp(np.clip(log_te, -18, 18)))
        rates = []
        for table, gap in zip(tables, gaps):
            k_up,k_down = maxwell_reversible_rates(table,te,gap)
            rates.append((ne*k_up,ne*k_down))
        (a,b),(c,d) = rates
        f01=a*n0-b*n1
        f12=c*n1-d*n2
        q10=config.quench10_s*n1
        q21=config.quench21_s*n2
        population_deriv = [-f01+q10, f01-f12-q10+q21, f12-q21]
        vibrational_power = nh2*J_PER_EV*(gaps[0]*f01+gaps[1]*f12)
        quench_heat = nh2*J_PER_EV*(gaps[0]*q10+gaps[1]*q21)
        elastic_power = electron_heat_capacity_J_m3_eV*(te-gas_Te_eV)/config.electron_gas_relaxation_s
        input_power=config.input_power_W_m3
        dU = input_power-vibrational_power-elastic_power
        dlogte=dU/(electron_heat_capacity_J_m3_eV*te)
        return population_deriv+[dlogte,quench_heat+elastic_power,input_power]

    t=np.linspace(0,config.duration_s,config.n_outputs)
    y0=[1.,0.,0.,np.log(config.te_initial_eV),0.,0.]
    sol=solve_ivp(rhs,(0,config.duration_s),y0,t_eval=t,method='BDF',
                  atol=1e-11,rtol=2e-8)
    if not sol.success:
        raise RuntimeError(f'Coupled electron-energy ODE failed: {sol.message}')
    p=sol.y[:3]
    te=np.exp(sol.y[3])
    ue=electron_heat_capacity_J_m3_eV*te
    uv=energy_levels@p
    ug=sol.y[4]
    injected=sol.y[5]
    total=ue+uv+ug
    residual=total-(u_initial+injected)
    scale=max(1.,float(np.max(np.abs(injected))),float(np.max(np.abs(total))))
    return {'time_s':sol.t,'fractions':p,'te_eV':te,'electron_energy_J_m3':ue,
            'vibrational_energy_J_m3':uv,'gas_heat_J_m3':ug,
            'injected_energy_J_m3':injected,'max_abs_energy_error_J_m3':float(np.max(abs(residual))),
            'relative_energy_error':float(np.max(abs(residual))/scale),
            'max_population_sum_error':float(np.max(abs(p.sum(axis=0)-1.))),
            'min_population':float(p.min()),
            'status':'Maxwellian electron-energy ODE; fixed ne and gas bath; no discharge',
            'cross_sections':'synthetic' if demo else 'caller-supplied; not independently verified'}


def run_energy_demo(out, config=ElectronEnergyConfig()):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=run_coupled_electron_energy(config=config)
    summary={k:v for k,v in d.items() if isinstance(v,(float,int,str))}
    summary.update(initial_Te_eV=float(d['te_eV'][0]),final_Te_eV=float(d['te_eV'][-1]),
                   max_Te_eV=float(d['te_eV'].max()),
                   total_power_deposited_J_m3=float(d['injected_energy_J_m3'][-1]))
    (out/'electron_energy_summary.json').write_text(json.dumps(summary,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(10,3.5))
    for i in range(3): ax[0].plot(d['time_s']*1e3,d['fractions'][i],label=f'v{i}')
    ax[0].legend(); ax[0].set(xlabel='Time [ms]',ylabel='Fraction')
    ax[1].plot(d['time_s']*1e3,d['te_eV'])
    ax[1].set(xlabel='Time [ms]',ylabel='Maxwellian electron Te [eV]')
    fig.suptitle('Synthetic excitation cross sections — not discharge validation')
    fig.tight_layout();fig.savefig(out/'electron_energy.png',dpi=155);plt.close(fig)
    return summary
